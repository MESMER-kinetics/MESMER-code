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
	
    // Read potential barrier details
   virtual bool ReadPE(Reaction* pReact, vector<double> &PE, vector<double> &E, double &Vary) ;
   const char* m_id;
  };

  //************************************************************
  //Global instance, defining its id (usually the only instance)
   DefinedTunnelingCoefficients theDefinedTunnelingCoefficients("Defined");
  //************************************************************

  bool DefinedTunnelingCoefficients::calculateCellTunnelingCoeffs(Reaction* pReact, vector<double>& TunnelingProbability){

    //TZ is the zpe of the TS
    const double TZ = pReact->get_relative_TSZPE();
    //barrier0 & barrier1 are the zpe corrected barrier heights in the forward/reverse directions
     int barrier0 = int(TZ) - int(pReact->get_relative_rctZPE());
     int barrier1 = int(TZ) - int(pReact->get_relative_pdtZPE());

    // Read in potential barrier details.
    	vector<double> PE, E ;
	E.clear()  ;
	PE.clear() ;
	
	//for the case where read p(E)'s are calculated variationally, therfore barrier0 may need correcting
	double VariationalThreshold;
    
    if (!ReadPE(pReact, PE, E, VariationalThreshold))
      return false ;
	//Correct barrier1 and barrier0 according to maximum in vibrationally adiabatic curve
	double VaryCorrection = VariationalThreshold-barrier0;
    
	if(VaryCorrection < 0){
	}
	else{
		barrier0 = barrier0 + VaryCorrection;
        barrier1 = barrier1 + VaryCorrection;
	}
    
	// Spline P(E)'s.

	Spline spline ;
	spline.Initialize(E, PE) ;

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

    // Calculatate macroscopic transmission coefficient for testing purposes

    if (pReact->getFlags().TunnellingCoeffEnabled){
      ctest << "\nTunneling coefficients for: " << pReact->getName();

      for(size_t i(0) ; i < TunnelingProbability.size() ; ++i){
        ctest << TunnelingProbability[i] << endl;
      }
      ctest << "}\n";
	}

    return true;
  }

  // Read potential barrier details
  bool DefinedTunnelingCoefficients::ReadPE(Reaction* pReact, vector<double> &PE, vector<double> &E, double &Vary) {

    Vary=0;
    
	// Read input data for P(E)'s

    PersistPtr pptran = pReact->get_TransitionState()->get_PersistentPointer();
    
    PersistPtr pp = pptran->XmlMoveTo("me:DefinedTunnelingCoefficients") ;
    if (!pp) {  
      // Force program to close if not p(E) dataavailable and print error message
      throw (std::runtime_error("Error: DefinedTunnelingCoefficients cannot be computed witout p(E) data")); 
    }
	
    while(pp = pp->XmlMoveTo("me:DefinedPE"))
	{
      double Ene = pp->XmlReadDouble("Energy");
      E.push_back(Ene) ;
	  Vary = max(Ene,Vary);

      double prob = pp->XmlReadDouble("pE");
      PE.push_back(prob) ;
    }

    return true ;
  }
  }
//namespace