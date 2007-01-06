//-------------------------------------------------------------------------------------------
//
// ReactionConnectivity.cpp
//
// Author: Struan Robertson 
// Date:   30/Dec/2006
//
// This header file contains the definition of the ReactionConnectivity class.
//
//-------------------------------------------------------------------------------------------

#include "ReactionConnectivity.h"

using namespace std ;

namespace mesmer
{

    void CReactionConnectivity::addreaction(Reaction *reaction){

        // Does table need re-sizing?

        if (m_reactantmap.size() == m_connectivitytable.size())
            m_connectivitytable.resize(2*m_connectivitytable.size()) ;

        // Get the name of the principal reactant.

        string reactant = reaction->get_ReactantName() ;

        // Search for name of reactant in reactant list.
        // If reactant not present add it.

        ReactantMap_itr itr ;
        itr = m_reactantmap.find(reactant) ;
        if ( itr == m_reactantmap.end() ){
            int idx = int(m_reactantmap.size()) ;
            m_reactantmap[reactant] = idx ;

            // Is reactant unimolecular? Mark diagonal matrix as one.

            if ( reaction->get_NumberOfReactants() == 1 ) 
                m_connectivitytable[idx][idx] = 1 ;
        }

        // Repeat for principal product.

        string product = reaction->get_ProductName() ;

        itr = m_reactantmap.find(product) ;
        if (itr == m_reactantmap.end() ){
            int idx = int(m_reactantmap.size()) ;
            m_reactantmap[product] = idx ;

            // Is product unimolecular? Mark diagonal matrix as one.

            if (reaction->get_NumberOfProducts() == 1)
                m_connectivitytable[idx][idx] = 1 ;
        }

        // Mark off-diagonal elements of table as unity to indicate connectivity.

        int indr = m_reactantmap[reactant];
        int indp = m_reactantmap[product];

        m_connectivitytable[indr][indp] = 1 ;
        m_connectivitytable[indp][indr] = 1 ;

    }

    void CReactionConnectivity::printconnectivitytable() const {

        cout << endl << "System: Connectivity table" << endl ;

        int msize = int(m_reactantmap.size()) ;

        for (int i(0) ; i < msize ; i++) {
            for (int j(0) ; j < msize ; j++) 
                cout  << " " << m_connectivitytable[i][j] ;

            cout << endl ;
        }

    }

    void CReactionConnectivity::get_unimolecularspecies(vector<string>& species) {

        for ( ReactantMap_itr itr = m_reactantmap.begin() ; itr != m_reactantmap.end() ; itr++ ){
           int idx = itr->second ;

           // Test diagonal element of connectivity table, 
           // if unity then species is unimolecular.

           if (m_connectivitytable[idx][idx] == 1) {
              species.push_back(itr->first) ;
           }
        }

    }

}
