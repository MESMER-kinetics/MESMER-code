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
    Reaction::Reaction(MoleculeManager *pMoleculeManager, const MesmerEnv& Env, const char *id):
m_Env(Env),
m_Name(id),
restartCalc(true),
m_PreExp(0.0),
m_NInf(0.0),
m_kfwd(0.0),
m_HeatOfReaction(0.0),
m_pMoleculeManager(pMoleculeManager),
m_pMicroRateCalculator(NULL),
m_pTunnelingCalculator(NULL),
m_ppPersist(),
m_rct1(NULL),
m_rct2(NULL),
m_pdt1(NULL),
m_TransitionState(NULL),
m_CellKfmc(),
m_CellKbmc(),
m_GrainKfmc(),
m_GrainKbmc()
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

    string pRef = ppmol->XmlReadValue("ref");
    if(pRef.size()){ // if got the name of the molecule
        string pType = ppmol->XmlReadValue("me:type");
        if(pType.size()){ // initialize molecule here with the specified type (need to know m_ppIOPtr)
            PersistPtr ppMolList = m_pMoleculeManager->get_PersistPtr();
            if(!ppMolList)
            {
                meErrorLog.ThrowError(__FUNCTION__, string("No molecules have been specified"), obWarning);
                return false;
            }
            pMol = m_pMoleculeManager->addmol(pRef, pType, ppMolList, getEnv());
        }
    }

    if(!pMol) {
        meErrorLog.ThrowError(__FUNCTION__, string("Cannot find molecule: "), obInfo);
        return NULL;
    }

    return pMol;
}

//
// Access microcanoincal rate coefficients.
//
void Reaction::get_MicroRateCoeffs(std::vector<double> &kmc) {
    calcGrnAvrgMicroRateCoeffs();
    kmc = m_GrainKfmc ;
}

double Reaction::get_ActivationEnergy(void) {
    if (!m_TransitionState) {
        cinfo << "No TransitionState for " << getName() << ", activation energy = 0.";
        return 0.0;
    }
    double zpeReactants = m_rct2 ? m_rct1->get_zpe() + m_rct2->get_zpe() : m_rct1->get_zpe();
    double AE = m_TransitionState->get_zpe() - zpeReactants;
    if(IsNan(AE)){
        cerr << "To use ILT for reaction " << getName() << " the ZPE of the transition state needs to be set.";
        exit(1);
    }
    return AE;
} ;


//
// Calculate grain averaged microcanonical rate coefficients.
//
bool Reaction::calcGrnAvrgMicroRateCoeffs() {
    if (restartCalc){
        if (m_CellKfmc.size()) m_CellKfmc.clear();
        if (m_CellKbmc.size()) m_CellKbmc.clear();
        restartCalc = false; // reset the flag

        // Calculate microcanonical rate coefficients.
        if(!m_pMicroRateCalculator->calculateMicroRateCoeffs(this, m_CellKfmc))
            return false;

        // Calculate Grain-averaged microcanonical rate coefficients.
        if (!grnAvrgMicroRateCoeffs())
            return false;

        // test grained microcanonical rate coefficients
        if (getEnv().microRateEnabled && !m_pMicroRateCalculator->testMicroRateCoeffs(this, m_ppPersist, m_GrainKfmc) )
            return false;
    }
    return true;
}

// calculate rate constant grained average for a forward or backward reaction
void Reaction::rateConstantGrnAvg(const int _MG,
                                  const int _gsz,
                                  const std::vector<double> &CellDOS,
                                  const std::vector<double> &CellRC,
                                  std::vector<double> &GrainRC)
{
    int idx1(0), idx2(0);
    for (int i = 0; i < _MG ; ++i ) {
        int idx3(idx1);

        // Calculate the number of states in a grain.
        double gNOS = .0 ;
        for (int j = 0 ; j < _gsz ; ++j, ++idx1){
            gNOS += CellDOS[idx1] ;
        }

        // Calculate average of the grain if it contains sum states.
        if ( gNOS > 0.0 ) {
            double gSum = .0;
            for (int j= 0 ; j < _gsz ; ++j, ++idx3){
                gSum += CellRC[idx3] * CellDOS[idx3] ;
            }
            GrainRC[idx2] = gSum/gNOS ;
            idx2++;
        }
    }

    // Issue warning if number of grains produced is less that requested.
    if ( idx2 != _MG ) {
        cerr << "Number of grains produced is not equal to that is requested" << endl
            << "Number of grains requested: " << _MG << endl
            << "Number of grains produced : " << idx2 << " in " << getName();
    }
    else{
        //      cinfo << "Number of grains requested: " << MaximumGrain << endl
        //            << "Number of grains produced : " << idx2 << " in " << getName() << endl;
    }
}



//
// Add microcanonical terms to collision operator
//
void Reaction::AddMicroRates(dMatrix *CollOptr,
                             isomerMap &isomermap,
                             const double rMeanOmega)
{
    // Calculate Microcanonical rate coefficients.
    calcGrnAvrgMicroRateCoeffs() ;

    // Add microcanonical rates to the collision operator.
    AddReactionTerms(CollOptr, isomermap, rMeanOmega) ;
}

}//namespace
