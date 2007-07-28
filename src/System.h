#ifndef GUARD_System_h
#define GUARD_System_h

//-------------------------------------------------------------------------------------------
//
// System.h 
//
// Author: Struan Robertson 
// Date:   11/Feb/2003
//
// This header file contains the declaration of the System class. This class is the route
// of all molecular and reaction data will contain information about any number of System.
// Reaction System inforamtion will be sorted in a vector of reaction maps. Molecular data
// is stored in the molecule manager.
//
//-------------------------------------------------------------------------------------------

#include <vector>
#include "MoleculeManager.h"
#include "ReactionManager.h"

namespace mesmer
{

    class System
    {
    public:

        System() ;
        ~System() ;

        //
        // Read and parse a data input file.
        //
        bool parse(PersistPtr ppIOPtr) ;

        //
        // Begin calculation.
        //
        void calculate() ;

        double igsz()                   { return GrainSize; }
        int    MAXGrn()                 { return MaxGrn; }
        int    MAXCell()                { return MaxCell; }
        double getEMin()                { return EMin; }
        bool   TestMicroRatesEnabled() { return bTestMicroRates; }

        std::vector<double> Temperatures; //Explicit emumeration of the temperatures requested in the input
        std::vector<double> Concentrations; //Bath gas concentrations: values for calculation. Use if Pressures is empty.
        std::vector<double> Pressures; //Bath gas pressures: values for calculation

    private:

        bool SetGrainParams();
        bool ReadRange(const std::string& name, std::vector<double>& vals,
                       PersistPtr ppbase, bool MustBeThere=true);
        void WriteMetadata();
        // 
        // Location of the molecule manager.
        //
        MoleculeManager *m_pMoleculeManager ;


        //
        // Location of the reaction maps.
        //
        ReactionManager *m_pReactionManager ;

        PersistPtr m_ppIOPtr; //level in XML file under <mesemer>

        double temp;
        double conc;

        double GrainSize;  //Grain size in cm-1
        int    MaxGrn;     //The number of grains
        int    MaxCell;    //The number of cells
        double MaxT;       //Maximum temperature for the purposes of setting the energy range

        double EMin, EMax; // The absolute lowest and highest energies in the system, cm-1

        bool   bTestMicroRates;  // Whether to output microcanonical rate coefficients

    } ;

    //Namespace global variable for accessing the one instance of System class
    extern System* pSys;

}//namespace
#endif // GUARD_System_h
