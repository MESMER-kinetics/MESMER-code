//-------------------------------------------------------------------------------------------
//
// WKBTunnelingCoefficients.h
//
// Author: Robin_Shannon
// Date:   _2011_02_22__10_02_55_
//
// Produces WKB tunneling coefficients
//
//-------------------------------------------------------------------------------------------
#include <vector>
#include <cmath>
#include <string>
#include <stdexcept>
#include "../System.h"
#include "../AssociationReaction.h"
#include "../IrreversibleExchangeReaction.h"
using namespace Constants;
using namespace std;

namespace mesmer
{
  class WKBTunnellingCoefficients : public TunnelingCalculator
  {
  public:

    ///Constructor which registers with the list of TunnellingCalculators in the base class
    WKBTunnellingCoefficients(const std::string& id) : TunnelingCalculator(id){}

    virtual ~WKBTunnellingCoefficients() {}

    virtual bool calculateCellTunnelingCoeffs(Reaction* pReact, std::vector<double>& TunnelingProbability);

    //Integration quadrature to obtain theta
    double WKBTunnellingCoefficients::Theta(vector<double> pot , vector<double> MEP , int Ene , double mu);

    //Integation quadrature for Boltzmann averaging to obtain macroscopic transmission coefficients 
    double WKBTunnellingCoefficients::Trans(vector<double> TunProb , double Vmax, double beta );


  };

  //************************************************************
  //Global instance, defining its id (usually the only instance)
  WKBTunnellingCoefficients theWKBTunnellingCoefficients("WKB");
  //************************************************************

  bool WKBTunnellingCoefficients::calculateCellTunnelingCoeffs(Reaction* pReact, vector<double>& TunnelingProbability){





    //TZ is the zpe of the TS
    const double TZ = pReact->get_relative_TSZPE();
    //barrier0 & barrier1 are the zpe corrected barrier heights in the forward/reverse directions
    const int barrier0 = int(TZ) - int(pReact->get_relative_rctZPE());
    const int barrier1 = int(TZ) - int(pReact->get_relative_pdtZPE());


    //get properties of vectors in which to include transmission coefficients
    const int MaximumCell = pReact->get_reactant()->getEnv().MaxCell;
    TunnelingProbability.clear();
    TunnelingProbability.resize(MaximumCell);

    Molecule * p_TransitionState = pReact->get_TransitionState();

    // read input data for barrier IRC

    PersistPtr pptran = p_TransitionState->get_PersistentPointer();

    PersistPtr pp = pptran->XmlMoveTo("me:IRCPotential") ;
    if(!pp)
    {  
      // Force program to close if not potential information available and print error message
      throw (std::runtime_error("Error: WKB calculation cannot proceed without a PES")); 
    }
    const char* p = pp->XmlReadValue("units", optional);
    string units = p ? p : "kJ/mol";
    double mu = pp->XmlReadDouble("ReducedMass", optional);
    if(IsNan(mu)) mu = 1.;
    vector<double> potential ;
    vector<double> distance ;

    while(pp = pp->XmlMoveTo("me:PotentialPoint"))
    {
      double distancePoint = pp->XmlReadDouble("ReacCoord", optional);
      if(IsNan(distancePoint))
        double  distancePoint = 0.0;
      distance.push_back(distancePoint) ;

      double potentialPoint = pp->XmlReadDouble("potential", optional);
      if(IsNan(potentialPoint))
        double  potentialPoint = 0.0;
      double convertedpotentialPoint = getConvertedEnergy(units, potentialPoint) * SpeedOfLight_in_cm*PlancksConstant_in_JouleSecond ; //Convert potential point into joules
      potential.push_back(convertedpotentialPoint) ;
    }

    //set transmission coefficients to 0 where no tunneling is possible;
    //where tunneling may occur, the transmission coefficients are calculated using a wkb formalism
    // as described by  B. C. Garrett and D. G. Truhlar, J. Phys. Chem., 1979, 83, 292

    for(int i = 0; i < MaximumCell; ++i){
      int E = i - barrier0;
      if ((E + barrier1) < 0){
        TunnelingProbability[i] = 0.0;
      }
      else if (E <= 0){
        TunnelingProbability[i] = 1 / (1 + ( exp ( 2 * Theta( potential, distance , E + barrier0 , mu)))) ;

      }

      //non classical reflection above the barrier
      else if (  E > 0 && E <= barrier0 ){
        TunnelingProbability[i] =1 -( 1 / (1 + ( exp ( 2 * Theta( potential, distance , 2*barrier0 - (E + barrier0), mu))))) ;

      }
      else if ( E > barrier0){
        TunnelingProbability[i] = 1;
      }

      {

        if(IsNan(TunnelingProbability[i])) TunnelingProbability[i] = 0.0;
      }
    }
    const double beta = pReact->getEnv().beta * 1./(SpeedOfLight_in_cm * PlancksConstant_in_JouleSecond) ; //Beta in J^-1
    // Calculatate macroscopic transmission coefficient for testing purposes
    double TransmissionCoefficient = ( 1 + ( 2 * beta )*Trans( TunnelingProbability , barrier0, beta ) );

    if (pReact->getFlags().TunnellingCoeffEnabled){
      ctest << "\nTransmission Coefficient for : " << pReact->getName() << "=" << TransmissionCoefficient;
      ctest << "\nTunneling coefficients for: " << pReact->getName();

      for(int i = 0; i < MaximumCell; ++i){
        ctest << TunnelingProbability[i] << endl;
      }
      ctest << "}\n";
    }

    return true;
  }

  double WKBTunnellingCoefficients::Theta(vector<double> pot , vector<double> MEP , int Ene , double mu){
    double integral=0.0;
    //Quadrature for integration, just using a simple trapezium like rule for now, need to include quadrature described in the reference above.
    for(int i = 0; i < pot.size()-1 ; ++i){
      //Only integrate over energies > current cell energy
      if (pot[i]>= (Ene*SpeedOfLight_in_cm*PlancksConstant_in_JouleSecond)){; // convert energy into joules
      double A = MEP[i];
      double fA = sqrt ( 2  * mu * amu * fabs( pot[i]- (Ene*SpeedOfLight_in_cm*PlancksConstant_in_JouleSecond)));
      double B = MEP[i+1];
      double fB = sqrt ( 2  * mu * amu * fabs(  pot[i+1] - (Ene*SpeedOfLight_in_cm*PlancksConstant_in_JouleSecond)));

      integral += (B - A) * ((fA + fB) / 2);
      }

    }

    double Theta = (integral * 2 * M_PI) / PlancksConstant_in_JouleSecond;
    return Theta;
  }
  
  //Quadrature for Boltzmann averaging
  double WKBTunnellingCoefficients::Trans(vector<double> TunProb , double Vmax, double beta ){
    double Trans=0.0;
    for(int i = 0; i < Vmax; ++i){
      double fA = sinh( (Vmax - i)*SpeedOfLight_in_cm*PlancksConstant_in_JouleSecond* beta ) * TunProb[i];
      double fB = sinh( (Vmax - (i+1))*SpeedOfLight_in_cm*PlancksConstant_in_JouleSecond * beta) * TunProb[i + 1];

      Trans +=  SpeedOfLight_in_cm*PlancksConstant_in_JouleSecond*(fA + fB) / 2 ;
    }
    return Trans;
  }

}//namespace

