//-------------------------------------------------------------------------------------------
//
// System.cpp 
//
// Author: Struan Robertson 
// Date:   11/Feb/2003
//
// This file contains the implementation of the System class.
//
//-------------------------------------------------------------------------------------------
#include "Persistence.h"
#include "Constants.h"
#include "System.h"

using namespace std ;
using namespace Constants ;

namespace mesmer
{
    System::System(): m_pMoleculeManager(0), m_pReactionManager(0) {

        m_pMoleculeManager = new MoleculeManager() ;

        m_pReactionManager = new ReactionManager(m_pMoleculeManager) ;

    }

    System::~System() { }

    //
    // Parse an input data file.
    //
    bool System::parse(PersistPtr ppIOPtr)
    {
        PersistPtr ppMolList = ppIOPtr->MoveTo("moleculeList");
        if(!ppMolList)
        {
            cerr << "No molecules have been specified" << endl;
            return false;
        }
        if(!m_pMoleculeManager->addmols(ppMolList))
            return false;;

        PersistPtr ppReacList = ppIOPtr->MoveTo("reactionList");
        if(!ppReacList)
        {
            cerr << "No reactions have been specified" << endl;
            return false;
        }if(!m_pReactionManager->addreactions(ppReacList))
            return false;

        PersistPtr ppConditions = ppIOPtr->MoveTo("me:conditions");
        if(!ppConditions)
        {
            cerr << "No conditions have been specified" << endl;
            return false;
        }
        const string Bgtxt = ppConditions->ReadValue("me:bathGas");
        if(!Bgtxt.size() || !(m_pMoleculeManager->find(Bgtxt)) )
        {
            cerr << "No bath gas specified" << endl;
            return false;
		} else {
			m_pMoleculeManager->set_BathGasMolecule(Bgtxt) ;
		}

        const char* Ttxt = ppConditions->ReadValue("me:temperature");
        const char* Ctxt = ppConditions->ReadValue("me:conc");
        if(!Ttxt || !Ctxt)
            return false;

        istringstream Tss(Ttxt);
        Tss >> temp;
        istringstream Css(Ctxt);
        Css >> conc;
        return true;
    }

    //
    // Begin calculation.
    //
    void System::calculate() 
    { 
        double beta = 1.0/(boltzmann*temp) ;

        // Build collison matrix for system.

        m_pReactionManager->BuildSystemCollisionOperator(beta, conc) ;

        m_pReactionManager->diagCollisionOperator() ;

        for (size_t i=0; i < m_pReactionManager->size() ; i++) {

            Reaction *reaction = (*m_pReactionManager)[i] ;

            reaction->CalcMicroRateCoeffs() ;

            // Work space to hold mircocanonical rates.

            vector<double> kmc(MAXGRN,0.0) ;

            reaction->get_MicroRateCoeffs(kmc, MAXGRN) ;

            CollidingMolecule *pmolecule = reaction->m_Reactant ;

            // for (int i(0) ; i < MAXCELL ; i++) 
            //     kmc[i] /= omega ;

            pmolecule->diagCollisionOperator() ;

            // Calculate matrix elements

            double kinf = pmolecule->matrixElement(MAXGRN-1,MAXGRN-1,kmc,MAXGRN) ;

            cout << endl ;
            formatFloat(cout, kinf, 6, 15) ;

        }

    }
}//namespace
