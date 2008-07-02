//-------------------------------------------------------------------------------------------
//
// Reaction.cpp
//
// Author: Struan Robertson
// Date:   23/Feb/2003
//
// This file contains the implementation of the Reaction class.
//
//-------------------------------------------------------------------------------------------
#include <limits>
#include "Reaction.h"

using namespace Constants ;
using namespace std;
using namespace mesmer;

namespace mesmer
{
  Reaction::Reaction(MoleculeManager *pMoleculeManager, const MesmerEnv& Env, const char *id)
    :m_ppPersist(),
    m_TransitionState(NULL),
    m_pMoleculeManager(pMoleculeManager),
    m_pMicroRateCalculator(NULL),
    m_pTunnelingCalculator(NULL),
    m_CellTSFlux(),
    m_GrainTSFlux(),
    m_GrainKfmc(),
    m_GrainKbmc(),
    m_Env(Env),
    m_Name(id),
    reCalcDOS(true),
    m_PreExp(0.0),
    m_NInf(0.0),
    m_kfwd(0.0),
    m_HeatOfReaction(0.0),
    m_HeatOfReactionInt(0)
  {}

  Reaction::~Reaction(){}

  /*
  Reaction::Reaction(const Reaction& reaction) {
  // Copy constructor - define later SHR 23/Feb/2003
  }

  Reaction& Reaction::operator=(const Reaction& reaction) {
  // Assignment operator - define later SHR 23/Feb/2003

  return *this ;
  }
  */


  //
  // Locate molecule in molecular map.
  //
  Molecule* Reaction::GetMolRef(PersistPtr pp)
  {
    Molecule* pMol = NULL;

    if(!pp) return NULL;
    PersistPtr ppmol = pp->XmlMoveTo("molecule");
    if(!ppmol) return NULL;

    string sRef = ppmol->XmlReadValue("ref");
    if(sRef.size()){ // if got the name of the molecule
      string sType = ppmol->XmlReadValue("me:type");
      if(sType.size()){ // initialize molecule here with the specified type (need to know m_ppIOPtr)
        PersistPtr ppMolList = m_pMoleculeManager->get_PersistPtr();
        if(!ppMolList)
        {
          cerr << "No molecules have been specified." << endl;
          return NULL;
        }
        pMol = m_pMoleculeManager->addmol(sRef, sType, ppMolList, getEnv());
      }
    }

    if(!pMol) {
      cinfo << "Failed to get a molecular reference." << endl;
      return NULL;
    }

    return pMol;
  }

  //
  // Access microcanonical rate coefficients.
  //
  void Reaction::get_MicroRateCoeffs(std::vector<double> &kmc) {
    calcGrnAvrgMicroRateCoeffs();
    kmc = m_GrainKfmc ;
  }

  //
  // Calculate grain averaged microcanonical rate coefficients.
  //
  bool Reaction::calcGrnAvrgMicroRateCoeffs() {
    if (reCalcDOS){
      if (m_CellTSFlux.size()) m_CellTSFlux.clear();

      // Calculate microcanonical rate coefficients.
      if(!m_pMicroRateCalculator->calculateMicroRateCoeffs(this))
        return false;

      // report TransitionState Flux in cells to test output
      const int MaximumCell = getEnv().MaxCell;
      if (getEnv().cellTSFluxEnabled){
        ctest << "\nTSFlux(e) cells for " << getName() << ":\n{\n";
        for (int i = 0; i < MaximumCell; ++i){
          ctest << m_CellTSFlux[i] << endl;
        }
        ctest << "}\n";
      }

      // Calculate Grain-averaged microcanonical rate coefficients.
      if (!grnAvrgMicroRateCoeffs())
        return false;

      // test grained microcanonical rate coefficients
      if (getEnv().microRateEnabled && !m_pMicroRateCalculator->testMicroRateCoeffs(this, m_ppPersist) )
        return false;
    }
    reCalcDOS = false; // reset the flag
    return true;
  }

  //
  // Access microcanonical rate coefficients - cell values are averaged
  // to give grain values. This code is similar to that in Molecule.cpp
  // and this averaging should be done there. SHR 19/Sep/2004.
  //
  bool Reaction::grnAvrgMicroRateCoeffs() {
    // This grain averaging of the microcanonical rate coefficients is
    // based on the view from the species that is
    // moving in the current reaction toward the opposite species.

    std::vector<double> shiftedTScellFlux;
    shiftTScellFlux(shiftedTScellFlux);

    // convert flux from cells to grains
    TSFluxCellToGrain(shiftedTScellFlux);

    // Calculate forward and backward grained microcanonical rate coefficients
    calcGrainRateCoeffs();

    return true;
  }

  // set the bottom energy of m_CellTSFlux
  void Reaction::setCellFluxBottom(const double fluxBottomZPE){
    m_FluxGrainZPE = (fluxBottomZPE - getEnv().EMin) / getEnv().GrainSize ; //convert to grain
    m_FluxCellOffset = int(fmod(fluxBottomZPE, getEnv().GrainSize));
  }

  // shift transitions state cell flux
  void Reaction::shiftTScellFlux(std::vector<double>& shiftedTScellFlux){
    int cellOffset = getTSFluxCellOffset();
    const int MaximumCell  = getEnv().MaxCell;
    for(int i = 0; i < cellOffset; ++i){
      shiftedTScellFlux.push_back(0.0);
    }
    for(int i = cellOffset, j = 0; i < MaximumCell; ++i, ++j){
      shiftedTScellFlux.push_back(m_CellTSFlux[j]);
    }
  }

  // calculate TSFlux in grains
  void Reaction::TSFluxCellToGrain(const std::vector<double>& shiftedTScellFlux)
  {
    const int maxGrn = getEnv().MaxGrn;
    const int grnSiz = getEnv().GrainSize;

    // resize m_GrainTSFlux to maxGrn and initialize all members to zero
    m_GrainTSFlux.clear();
    m_GrainTSFlux.resize(maxGrn, 0.0);

    int cIdx = 0; // cell iterator

    for (int i = 0; i < maxGrn ; ++i) {
      for (int j = 0; j < grnSiz; ++j, ++cIdx) {
        m_GrainTSFlux[i] += shiftedTScellFlux[cIdx];
      }
    }

    if (getEnv().grainTSFluxEnabled){
      ctest << "\nTSFlux(e) grains for " << getName() << ":\n{\n";
      for (int i = 0; i < maxGrn; ++i){
        ctest << m_GrainTSFlux[i] << endl;
      }
      ctest << "}\n";
    }
  }

  //this function retrieves the threshold energy for a reaction
  double Reaction::get_ThresholdEnergy(void) const {
    if (!m_TransitionState) {
      cinfo << "No TransitionState for " << getName() << ", threshold energy = 0." << endl;
      return 0.0;
    }

    double ThresholdEnergy = get_relative_TSZPE() - get_relative_rctZPE();
    if(IsNan(ThresholdEnergy)){
      cerr << "Reaction " << getName() << " has no threshold energy.";
      exit(1);
    }
    return ThresholdEnergy;
  } ;

  void Reaction::setHeatOfReaction(const double pdtZPE, const double rctZPE){
    m_HeatOfReaction = pdtZPE - rctZPE;
    m_HeatOfReactionInt = int(pdtZPE) - int(rctZPE);
  }

  void Reaction::setHeatOfReaction(const double value){
    m_HeatOfReaction = value;
    m_HeatOfReactionInt = int(value);
  }

}//namespace
