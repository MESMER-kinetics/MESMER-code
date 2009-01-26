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
  Reaction::Reaction(MoleculeManager *pMoleculeManager, const MesmerEnv& Env, MesmerFlags& Flags, const char *id)
    :m_ppPersist(),
    m_TransitionState(NULL),
    m_pMoleculeManager(pMoleculeManager),
    m_pMicroRateCalculator(NULL),
    m_pTunnelingCalculator(NULL),
    m_FluxCellZPE(0.0),
    m_FluxGrainZPE(0.0),
    m_FluxCellOffset(0),
    m_CellFlux(),
    m_GrainFlux(),
    m_GrainKfmc(),
    m_GrainKbmc(),
    m_fwdGrnCanonicalRate(0.0),
    m_rvsGrnCanonicalRate(0.0),
    m_fwdCellCanonicalRate(0.0),
    m_rvsCellCanonicalRate(0.0),
    m_Env(Env),
    m_Flags(Flags),
    m_Name(id),
    m_GrnFluxFirstNonZeroIdx(0),
    m_EffGrainedFwdThreshold(0),
    m_EffGrainedRvsThreshold(0),
    reCalcDOS(true),
    m_PreExp(0.0),
    m_NInf(0.0),
    m_TInf(298.0),
    m_EInf(0.0),
    m_kfwd(0.0),
    m_ERConc(0.0)
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
  Molecule* Reaction::GetMolRef(PersistPtr pp, const char* defaultType)
  {
    Molecule* pMol = NULL;

    if(!pp) return NULL;
    PersistPtr ppmol = pp->XmlMoveTo("molecule");
    if(!ppmol) return NULL;

    const char* reftxt = ppmol->XmlReadValue("ref");//using const char* in case NULL returned
    if(reftxt) // if got the name of the molecule
    {
      const char* typetxt = ppmol->XmlReadValue("me:type");
      if(!typetxt && defaultType)
        typetxt=defaultType;
      if(typetxt){ // initialize molecule here with the specified type (need to know m_ppIOPtr)
        PersistPtr ppMolList = m_pMoleculeManager->get_PersistPtr();
        if(!ppMolList)
        {
          cerr << "No molecules have been specified." << endl;
          return NULL;
        }
        pMol = m_pMoleculeManager->addmol(string(reftxt), string(typetxt), ppMolList, getEnv(), getFlags());
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
      if (m_CellFlux.size()) m_CellFlux.clear();

      // Calculate microcanonical rate coefficients.
      if(!m_pMicroRateCalculator->calculateMicroRateCoeffs(this))
        return false;

      // report Transition State Flux in cells to test output
      const int MaximumCell = getEnv().MaxCell;
      if (getFlags().cellFluxEnabled){
        ctest << "\nFlux(e) cells for " << getName() << ":\n{\n";
        for (int i = 0; i < MaximumCell; ++i){
          ctest << m_CellFlux[i] << endl;
        }
        ctest << "}\n";
      }

      // Calculate Grain-averaged microcanonical rate coefficients.
      if (!grnAvrgMicroRateCoeffs())
        return false;

      // test grained microcanonical rate coefficients
      if (getFlags().microRateEnabled && !m_pMicroRateCalculator->testMicroRateCoeffs(this, m_ppPersist) )
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

    std::vector<double> shiftedCellFlux;
    shiftCellFlux(shiftedCellFlux);

    // convert flux from cells to grains
    fluxCellToGrain(shiftedCellFlux);

    // Calculate forward and backward grained microcanonical rate coefficients
    calcGrainRateCoeffs();

    return true;
  }

  // set the bottom energy of m_CellFlux
  void Reaction::setCellFluxBottom(const double fluxBottomZPE){
    m_FluxCellZPE = fluxBottomZPE;
    m_FluxGrainZPE = fluxBottomZPE / getEnv().GrainSize ; //convert to grain
    m_FluxCellOffset = int(fmod(fluxBottomZPE, getEnv().GrainSize));
  }

  // shift transition state cell flux
  void Reaction::shiftCellFlux(std::vector<double>& shiftedCellFlux){
    int cellOffset = getFluxCellOffset();
    const int MaximumCell  = getEnv().MaxCell;
    for(int i = 0; i < cellOffset; ++i){
      shiftedCellFlux.push_back(0.0);
    }
    for(int i = cellOffset, j = 0; i < MaximumCell; ++i, ++j){
      shiftedCellFlux.push_back(m_CellFlux[j]);
    }
  }

  // calculate flux in grains
  void Reaction::fluxCellToGrain(const std::vector<double>& shiftedCellFlux)
  {
    const int maxGrn = getEnv().MaxGrn;
    const int grnSiz = getEnv().GrainSize;

    // resize m_GrainFlux to maxGrn and initialize all members to zero
    m_GrainFlux.clear();
    m_GrainFlux.resize(maxGrn, 0.0);

    int cIdx = 0; // cell iterator

    for (int i = 0; i < maxGrn ; ++i) {
      for (int j = 0; j < grnSiz; ++j, ++cIdx) {
        m_GrainFlux[i] += shiftedCellFlux[cIdx];
      }
    }

    if (getFlags().grainFluxEnabled){
      ctest << "\nFlux(e) grains for " << getName() << ":\n{\n";
      for (int i = 0; i < maxGrn; ++i){
        ctest << m_GrainFlux[i] << endl;
      }
      ctest << "}\n";
    }
  }

  //this function retrieves the activation/threshold energy for an association reaction
  double Reaction::get_ThresholdEnergy(void) {
    // ILT
    if (m_pMicroRateCalculator->getName() == "Mesmer ILT"){
      if (IsNan(m_EInf.get_value())){
        cerr << "No E_infinity provided for Reaction " << getName();
        exit(1);
      }
      if (m_EInf.get_value() < 0.0){
        cerr << "Providing negative E_infinity in Reaction " << getName() << " is invalid.";
      }
      if (m_isRvsILTpara){
        const double tempv = m_EInf.get_value() > 0.0 ? m_EInf.get_value() + getHeatOfReaction() : getHeatOfReaction();
        return tempv;
      }
      return m_EInf.get_value();
    }

    // Not ILT
    if (!m_TransitionState) {
      cerr << "No Transition State for " << getName();
      exit(1);
    }

    return (get_relative_TSZPE() - get_relative_rctZPE());
  }

  void Reaction::calcFluxFirstNonZeroIdx(void) {
    double thresh = get_ThresholdEnergy();
    double RxnHeat = getHeatOfReaction();
    if(thresh<0.0)
      m_GrnFluxFirstNonZeroIdx = int(-thresh/m_Env.GrainSize);
    else if(thresh>0.0 && thresh<RxnHeat)
      m_GrnFluxFirstNonZeroIdx = int(RxnHeat - thresh)/m_Env.GrainSize;
    else
      m_GrnFluxFirstNonZeroIdx = 0;
  };

  // Read excess reactant concentration
  bool Reaction::ReadExcessReactantConcentration(PersistPtr ppReac){
    const char* pERConctxt = ppReac->XmlReadValue("me:excessReactantConc",false);
    if (!pERConctxt){
      cerr << "Concentration of excess reactant has not been specified.";
      return false;
    } else {
      stringstream s3(pERConctxt) ;
      s3 >> m_ERConc ;
    }
    return true;
  }

  // Read ILT parameters
  bool Reaction::ReadILTParameters(PersistPtr ppReac) {
    const char* pActEnetxt = ppReac->XmlReadValue("me:activationEnergy", false);
    if (pActEnetxt)
    {
      PersistPtr ppActEne = ppReac->XmlMoveTo("me:activationEnergy") ;
      m_isRvsILTpara = ppActEne->XmlReadValue("reverse", false); //specify the direction of the following ILT parameters
      double tmpvalue = 0.0;
      stringstream s2(pActEnetxt); s2 >> tmpvalue ;
      const char* unitsTxt = ppActEne->XmlReadValue("units", false);
      string unitsInput;
      if (unitsTxt){
        unitsInput = unitsTxt;
      }
      else{
        unitsInput = "kJ/mol";
      }
      const char* pLowertxt = ppActEne->XmlReadValue("lower", false);
      const char* pUppertxt = ppActEne->XmlReadValue("upper", false);
      const char* pStepStxt = ppActEne->XmlReadValue("stepsize", false);
      double value(getConvertedEnergy(unitsInput, tmpvalue));
      if (pLowertxt && pUppertxt){
        double tmpvalueL(0.0), tmpvalueU(0.0), tmpstepsize(0.0);
        stringstream s3(pLowertxt), s4(pUppertxt), s5(pStepStxt); s3 >> tmpvalueL; s4 >> tmpvalueU; s5 >> tmpstepsize;
        double valueL(getConvertedEnergy(unitsInput, tmpvalueL));
        double valueU(getConvertedEnergy(unitsInput, tmpvalueL));
        double stepsize(getConvertedEnergy(unitsInput, tmpstepsize));
        set_EInf(valueL, valueU, stepsize);
      }
      else{
        set_EInf(value);
      }
    }
    else{
      cerr << "Specifying ILT without activation energy provided in reaction " << this->getName() << ". Please correct input file.";
      return false;
    }

    const char* pPreExptxt = ppReac->XmlReadValue("me:preExponential", false);
    if (pPreExptxt)
    {
      PersistPtr ppPreExponential = ppReac->XmlMoveTo("me:preExponential") ;
      double value = 0.0;
      stringstream s2(pPreExptxt); s2 >> value ;
      const char* pLowertxt = ppPreExponential->XmlReadValue("lower", false);
      const char* pUppertxt = ppPreExponential->XmlReadValue("upper", false);
      const char* pStepStxt = ppPreExponential->XmlReadValue("stepsize", false);
      if (pLowertxt && pUppertxt){
        double valueL(0.0), valueU(0.0), stepsize(0.0);
        stringstream s3(pLowertxt), s4(pUppertxt), s5(pStepStxt); s3 >> valueL; s4 >> valueU; s5 >> stepsize;
        set_PreExp(valueL, valueU, stepsize);
      }
      else{
        set_PreExp(value);
      }
    }
    else{
      cerr << "Specifying ILT without pre-exponential term provided in reaction " << this->getName() << ". Please correct input file.";
      return false;
    }

    const char* pNInftxt = ppReac->XmlReadValue("me:nInfinity", false);
    if (pNInftxt)
    {
      PersistPtr ppNInf = ppReac->XmlMoveTo("me:nInfinity") ;
      double value = 0.0;
      stringstream s2(pNInftxt); s2 >> value ;
      const char* pLowertxt = ppNInf->XmlReadValue("lower", false);
      const char* pUppertxt = ppNInf->XmlReadValue("upper", false);
      const char* pStepStxt = ppNInf->XmlReadValue("stepsize", false);
      if (pLowertxt && pUppertxt){
        double valueL(0.0), valueU(0.0), stepsize(0.0);
        stringstream s3(pLowertxt), s4(pUppertxt), s5(pStepStxt); s3 >> valueL; s4 >> valueU; s5 >> stepsize;
        set_NInf(valueL, valueU, stepsize);
      }
      else{
        set_NInf(value);
      }
    }

    const char* pTInftxt = ppReac->XmlReadValue("me:TInfinity");  // read XML designation
    if (!pTInftxt){           // if there is no Tinfinity
      cinfo << "Tinfinity has not been specified; set to the default value of 298 K" << endl;
    }
    else{                     // if there is a Tinfinity
      double TInf(298.0);     // set Tinf to 298
      stringstream s3(pTInftxt);    // initialize stringstream
      s3 >> TInf;             // use the stringstream to get Tinfinty from the input file
      if(TInf <= 0){          // if Tinf is <= 0, set to 298
        cinfo << "Tinfinity is less than or equal to 0; set to the default value of 298 K" << endl;
      }
      else
        set_TInf(TInf);         // else set Tinf to the value in the input
    }

    return true;
  }

  // Read parameters requires to determine reaction heats and rates.
  bool Reaction::ReadRateCoeffParameters(PersistPtr ppReac) {

    // Determine the method of MC rate coefficient calculation.
    const char* pMCRCMethodtxt = ppReac->XmlReadValue("me:MCRCMethod") ;
    if(pMCRCMethodtxt)
    {
      m_pMicroRateCalculator = MicroRateCalculator::Find(pMCRCMethodtxt);
      if(!m_pMicroRateCalculator)
      {
        cerr << "Unknown method " << pMCRCMethodtxt
          << " for the determination of Microcanonical rate coefficients in reaction "
          << getName();
        return false;
      }
    }

    // Read the transition state (if present)
    PersistPtr ppTransitionState = ppReac->XmlMoveTo("me:transitionState") ;
    if (ppTransitionState)
    {
      Molecule* pTrans = GetMolRef(ppTransitionState,"transitionState");
      if(pTrans) m_TransitionState = pTrans;
    }

    //---------------------------------------------------------
    // Microcanonical rate constants methods dependent section.
    if (m_pMicroRateCalculator->getName() == "Mesmer ILT"){
      if (ppTransitionState){
        cerr << "Reaction " << getName() << " uses ILT method, which should not have transition state.";
        return false;
      }
      cinfo << "ILT method chosen, look for ILT expressions" << endl;
      if (!ReadILTParameters(ppReac)) return false;
      const char* pTunnelingtxt = ppReac->XmlReadValue("me:tunneling") ;
      if(pTunnelingtxt)
      {
        cerr << "Tunneling parameter in Reaction " << getName() << " is invalid in ILT.";
        return false;
      }
    }
    if (m_pMicroRateCalculator->getName() == "Simple RRKM"){
      if (!ppTransitionState){
        cerr << "Reaction " << getName() << " uses RRKM method, which should have transition state.";
        return false;
      }
      // Determine the method of estimating tunneling effect.
      const char* pTunnelingtxt = ppReac->XmlReadValue("me:tunneling") ;
      if(pTunnelingtxt)
      {
        m_pTunnelingCalculator = TunnelingCalculator::Find(pTunnelingtxt);
        if(!m_pTunnelingCalculator)
        {
          cerr << "Unknown method " << pTunnelingtxt
            << " for the determination of tunneling coefficients in reaction "
            << getName();
          return false;
        }
      }
      else{
        cinfo << "No tunneling method was found for " << getName() << endl;
      }
    }
    //
    //---------------------------------------------------------

    if (!isUnimolecular()){
      cinfo << "Not a unimolecular reaction: look for excess reactant concentration." << endl;
      if (!ReadExcessReactantConcentration(ppReac)) return false;
    }

    return true ;
  }


}//namespace
