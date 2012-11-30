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
#include "../Reaction.h"

using namespace Constants;
using namespace std;

namespace mesmer
{
  class WKBTunnellingCoefficients : public TunnelingCalculator
  {
  public:

    ///Constructor which registers with the list of TunnellingCalculators in the base class
    WKBTunnellingCoefficients(const char* id) : m_id(id){ Register(); }

    virtual ~WKBTunnellingCoefficients() {}
    virtual const char* getID()  { return m_id; }

    virtual bool calculateCellTunnelingCoeffs(Reaction* pReact, std::vector<double>& TunnelingProbability);

  private:

    // Read potential barrier details
    bool ReadPotential(Reaction* pReact, vector<double> &potential, vector<double> &distance, double &mu) ;

    //Integration quadrature to obtain theta
    double Theta(const vector<double> &pot , const vector<double> &MEP , int Ene , double mu);

    //Integation quadrature for Boltzmann averaging to obtain macroscopic transmission coefficients 
    double Trans(const vector<double> &TunProb , double Vmax, double beta );
  private:
    const char* m_id;
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

    // Read in potential barrier details.
    vector<double> potential, distance ;
    double mu(0.0) ;
    if (!ReadPotential(pReact, potential, distance, mu))
      return false ;

    // Set transmission coefficients to 0 where no tunneling is possible;
    // where tunneling may occur, the transmission coefficients are calculated using a wkb formalism
    // as described by  B. C. Garrett and D. G. Truhlar, J. Phys. Chem., 1979, 83, 292
    const int MaximumCell = pReact->get_reactant()->getEnv().MaxCell;
    TunnelingProbability.clear();
    TunnelingProbability.resize(MaximumCell);
    for(int i = 0; i < MaximumCell; ++i){
      int E = i - barrier0;
      TunnelingProbability[i] = 0.0;
      if ((E + barrier1) < 0) {
        TunnelingProbability[i] = 0.0;
      } else if (E <= 0) {
        TunnelingProbability[i] = 1.0 /(1.0 + ( exp (2.0 * Theta( potential, distance , E + barrier0 , mu)))) ;
      } else if (E > 0 && E <= barrier0 ) { 
        //non classical reflection above the barrier
        TunnelingProbability[i] = 1.0 - (1.0 /(1.0 + ( exp (2.0 * Theta( potential, distance , barrier0 - E, mu))))) ;
      } else if (E > barrier0) {
        TunnelingProbability[i] = 1.0 ;
      } else {
        // This branch should never be executed.
      }
    }

    // Calculatate macroscopic transmission coefficient for testing purposes
    const double beta = pReact->getEnv().beta ; 
    double TransmissionCoefficient = ( 1.0 + ( 2.0 * beta )*Trans( TunnelingProbability , barrier0, beta ) );

    if (pReact->getFlags().TunnellingCoeffEnabled){
      ctest << "\nTransmission Coefficient for : " << pReact->getName() << "=" << TransmissionCoefficient;
      ctest << "\nTunneling coefficients for: " << pReact->getName();

      for(size_t i(0) ; i < TunnelingProbability.size() ; ++i){
        ctest << TunnelingProbability[i] << endl;
      }
      ctest << "}\n";
    }

    return true;
  }

  // Read potential barrier details
  bool WKBTunnellingCoefficients::ReadPotential(Reaction* pReact, vector<double> &potential, vector<double> &distance, double &mu) {

    // Read input data for barrier IRC.

    PersistPtr pptran = pReact->get_TransitionState()->get_PersistentPointer();
    
    PersistPtr pp = pptran->XmlMoveTo("me:IRCPotential") ;
    if (!pp) {  
      // Force program to close if not potential information available and print error message
      throw (std::runtime_error("Error: WKB calculation cannot proceed without a PES")); 
    }
    const char* p = pp->XmlReadValue("units", optional);
    string units = p ? p : "kJ/mol";
    
    mu = pp->XmlReadDouble("ReducedMass", optional);
    if(IsNan(mu)) mu = 1.0 ;

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

    return true ;
  }

  double WKBTunnellingCoefficients::Theta(const vector<double> &pot, const vector<double> &MEP, int Ene, double mu){

    // Quadrature for integration, just using a simple trapezium like rule for now,
    // need to include quadrature described in the reference above.
    double integral(0.0);
    double energy = double(Ene)*SpeedOfLight_in_cm*PlancksConstant_in_JouleSecond ;
    for(size_t i(0) ; i < pot.size()-1 ; ++i){
      //Only integrate over energies > current cell energy
      if (pot[i] >= energy) { 
        double fA = sqrt ( 2.0  * mu * amu * fabs(pot[i]   - energy));
        double fB = sqrt ( 2.0  * mu * amu * fabs(pot[i+1] - energy));

        integral += (MEP[i+1] - MEP[i]) * ((fA + fB) / 2);
      }
    }

    return (integral * 2.0 * M_PI) / PlancksConstant_in_JouleSecond;
  }

  double WKBTunnellingCoefficients::Trans(const vector<double> &TunProb , double Vmax, double beta ){

    // Trapezium quadrature for Boltzmann averaging. Units of energy are cm-1.
    double Trans(0.0); 
    double dEne(1.0) ; // cm -1. 
    for(int i = 0; i < int(Vmax) ; ++i){
      double fA = sinh( (Vmax - double(i))   * beta) * TunProb[i];
      double fB = sinh( (Vmax - double(i+1)) * beta) * TunProb[i + 1];

      Trans += (fA + fB) * 0.5 ;
    }
    return Trans*dEne ;
  }

}//namespace

