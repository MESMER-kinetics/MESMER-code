//-------------------------------------------------------------------------------------------
//
// DissociationReaction.cpp
//
// Author: Struan Robertson
// Date:   30/Dec/2007
//
// This file contains the implementation of the DissociationReaction class.
//
//-------------------------------------------------------------------------------------------
#include <limits>
#include "DissociationReaction.h"

using namespace Constants ;
using namespace std;
using namespace mesmer;

namespace mesmer
{
  //
  // Add dissociation reaction terms to collision matrix.
  //
  void DissociationReaction::AddReactionTerms(dMatrix *CollOptr, isomerMap &isomermap, const double rMeanOmega) {

    // Locate reactant in system matrix.
    const int rctLocation = isomermap[dynamic_cast<CollidingMolecule*>(m_rct1)] ;

    for ( int i = 0 ; i < dynamic_cast<CollidingMolecule*>(m_rct1)->get_colloptrsize() ; ++i ) {
      int ii(rctLocation + i) ;
      (*CollOptr)[ii][ii] -= rMeanOmega * m_GrainKfmc[i] ;                            // Forward loss reaction.
    }
  }


}//namespace
