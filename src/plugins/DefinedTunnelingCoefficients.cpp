//-------------------------------------------------------------------------------------------
//
// DefinedTunnelingCoefficients.h
//
// Author: Robin_Shannon
// Date:   _2011_02_22__10_02_55_
//
// Reads external tunneling probabilities
//
//-------------------------------------------------------------------------------------------
#include <vector>
#include <cmath>
#include <string>
#include <stdexcept>
#include "../System.h"
#include "../Spline.h"
#include "../AssociationReaction.h"
#include "../IrreversibleExchangeReaction.h"
using namespace Constants;

namespace mesmer
{
  class DefinedTunnelingCoefficients : public TunnelingCalculator
  {
  public:

    ///Constructor which registers with the list of TunnellingCalculators in the base class
    DefinedTunnelingCoefficients(const char* id) : m_id(id){ Register(); }

    virtual ~DefinedTunnelingCoefficients() {}
    virtual const char* getID()  { return m_id; }

    virtual bool calculateCellTunnelingCoeffs(Reaction* pReact, std::vector<double>& TunnelingProbability);

  private:
   // Read potential barrier details. Called from Reaction.cpp ~line 274 via ParseForPlugin
   // and store in member variables.
   virtual bool ParseData(PersistPtr pp);

   const char* m_id;
   vector<double> m_Energy;
   vector<double> m_PE;
   double m_VariationalThreshold;
  };

  //************************************************************
  //Global instance, defining its id
   DefinedTunnelingCoefficients theDefinedTunnelingCoefficients("Defined");
  //************************************************************

  bool DefinedTunnelingCoefficients::calculateCellTunnelingCoeffs(Reaction* pReact, vector<double>& TunnelingProbability){

    //TZ is the zpe of the TS
    const double TZ = pReact->get_relative_TSZPE();
    //barrier0 & barrier1 are the zpe corrected barrier heights in the forward/reverse directions
    int barrier0 = int(TZ) - int(pReact->get_relative_rctZPE());
    int barrier1 = int(TZ) - int(pReact->get_relative_pdtZPE());

  
    //Correct barrier1 and barrier0 according to maximum in vibrationally adiabatic curve
    int VaryCorrection = static_cast<int>(m_VariationalThreshold-barrier0);
    
    if(VaryCorrection < 0){
    }
    else{
      barrier0 = barrier0 + VaryCorrection;
          barrier1 = barrier1 + VaryCorrection;
    }
    
    // Spline P(E)'s.

    Spline spline ;
    spline.Initialize(m_Energy, m_PE) ;

    // Set transmission coefficients to 0 where no tunneling is possible;
    // where tunneling may occur, the transmission coefficients are calculated using a wkb formalism
    // as described by  B. C. Garrett and D. G. Truhlar, J. Phys. Chem., 1979, 83, 292
    const int MaximumCell = pReact->get_reactant()->getEnv().MaxCell;
    TunnelingProbability.clear();
    TunnelingProbability.resize(MaximumCell);
    for(int i = 0; i < MaximumCell; ++i){
      int E = i - barrier0;
      if ((E + barrier1) < 0) {
        TunnelingProbability[i] = 0.0;
      } else if (E <= 0) {
        TunnelingProbability[i] = spline.Calculate(double(i)) ;;
      } else if (E > 0 && E <= barrier0 ) { 
        //non classical reflection above the barrier
        TunnelingProbability[i] = 1.0 - spline.Calculate(double(barrier0-(i-barrier0))) ;
      } else if (E > barrier0) {
        TunnelingProbability[i] = 1.0 ;
      } else {
        // This branch should never be executed.
      }
    }

    // Calculate macroscopic transmission coefficient for testing purposes

    if (pReact->getFlags().TunnellingCoeffEnabled){
      ctest << "\nTunneling coefficients for: " << pReact->getName();

      for(size_t i(0) ; i < TunnelingProbability.size() ; ++i){
        ctest << TunnelingProbability[i] << endl;
      }
      ctest << "}\n";
	}

    return true;
  }

  bool DefinedTunnelingCoefficients::ParseData(PersistPtr pp1)
  {
    m_VariationalThreshold=0;
    PersistPtr pp = pp1->XmlMoveTo("me:DefinedTunnelingCoefficients");
    //If data is not under <me:tunneling> in <reaction>, look in TS
    if (!pp)
    {
      pp1 = getParent()->get_TransitionState()->get_PersistentPointer();
      pp = pp1->XmlMoveTo("me:DefinedTunnelingCoefficients");
    }
    if (!pp)
      return false;

    bool dataFound=false;
    while(pp = pp->XmlMoveTo("me:DefinedPE"))
    {
      double Ene = pp->XmlReadDouble("energy", optional);
      m_Energy.push_back(Ene) ;
      m_VariationalThreshold = std::max(Ene,m_VariationalThreshold);

      double prob = pp->XmlReadDouble("pE", optional);
      m_PE.push_back(prob) ;
      if(IsNan(Ene) || IsNan(prob))
      {
        cerr << "Missing energy or PE" << endl;
        return false;
      }
      else
        dataFound=true;
    }
    if(dataFound)
      cinfo << "Data for DefinedTunnelingCoefficients found" << endl;

    return true ;
  }
  }
//namespace