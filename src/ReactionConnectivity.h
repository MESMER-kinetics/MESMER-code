#ifndef GUARD_ReactionConnectivity_h
#define GUARD_ReactionConnectivity_h

//-------------------------------------------------------------------------------------------
//
// ReactionConnectivity.h
//
// Author: Struan Robertson
// Date:   30/Dec/2006
//
// This header file contains the declaration of the ReactionConnectivity class.
//
//-------------------------------------------------------------------------------------------

//#include "Reaction.h"

namespace mesmer
{

class CReactionConnectivity {

public:

  CReactionConnectivity(): m_reactantmap(), m_connectivitytable(10) {};

  virtual ~CReactionConnectivity(){} ;

  void addreaction(Reaction *reaction) ;

  void printconnectivitytable() const ;

  void get_unimolecularspecies(std::vector<std::string>& species) ;

private:

  std::map<std::string, int> m_reactantmap ;

  typedef std::map<std::string, int>::const_iterator ReactantMap_itr ;

  Matrix<int> m_connectivitytable ;

} ;

}

#endif // GUARD_ReactionConnectivity_h
