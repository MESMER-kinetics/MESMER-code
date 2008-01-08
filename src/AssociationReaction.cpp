//-------------------------------------------------------------------------------------------
//
// AssociationReaction.cpp
//
// Author: Struan Robertson
// Date:   30/Dec/2007
//
// This file contains the implementation of the AssociationReaction class.
//
//-------------------------------------------------------------------------------------------
#include <limits>
#include "AssociationReaction.h"

using namespace Constants ;
using namespace std;
using namespace mesmer;

namespace mesmer
{
  //
  // Add (REVERSIBLE) association reaction terms to collision matrix.
  //
  void AssociationReaction::AddReactionTerms(dMatrix      *CollOptr,
                                             isomerMap    &isomermap,
                                             const double rMeanOmega)
  {
    // Locate isomers in system matrix.
    const int pdtLoc =      isomermap[dynamic_cast<CollidingMolecule*>(m_pdt1)] ;
    const int sL     = (*m_sourceMap)[dynamic_cast<SuperMolecule    *>(m_srct)] ;

    // Get equilibrium constant.
    const double Keq = calcEquilibriumConstant() ;

    // Get Boltzmann distribution for detailed balance.
    const int MaximumGrain = m_Env.MaxGrn ;
    vector<double> adductBoltz(MaximumGrain, 0.0) ;
    m_pdt1->grnBoltzDist(adductBoltz) ;

    const int colloptrsize = dynamic_cast<CollidingMolecule*>(m_pdt1)->get_colloptrsize() ;

    double DissRateCoeff(0.0) ;

    const int idx = m_pdt1->get_grnZpe() - m_rct1->get_grnZpe() ;
    for ( int i = max(0,-idx) ; i < min(colloptrsize,(colloptrsize-idx)) ; ++i ) {
      int ll = i + idx ;
      int pL(pdtLoc + ll) ;

      (*CollOptr)[pL][pL] -= rMeanOmega * m_GrainKfmc[ll] ;                           // Forward loss reaction.
      (*CollOptr)[pL][sL]  = rMeanOmega * m_GrainKfmc[ll]*sqrt(adductBoltz[ll]/Keq) ; // Reactive gain.
      (*CollOptr)[sL][pL]  = (*CollOptr)[pL][sL] ;                                    // Reactive gain.
      DissRateCoeff       += m_GrainKfmc[ll]*adductBoltz[ll] ;
    }
    (*CollOptr)[sL][sL] -= DissRateCoeff/Keq ;       // Backward loss reaction from detailed balance.
  }

}//namespace
