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
#include "tinyxml.h"

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
   bool parse(TiXmlElement* root) ;

   //
   // Begin calculation.
   //
   void calculate() ;

private:

   // 
   // Location of the molecule manager.
   //
   MoleculeManager *m_pMoleculeManager ;


   //
   // Location of the reaction maps.
   //
   ReactionManager *m_pReactionManager ;

private:
  BathGasMolecule* m_pBathGasMolecule;

  double temp;
  double conc;


} ;
}//namespace
#endif // GUARD_System_h
