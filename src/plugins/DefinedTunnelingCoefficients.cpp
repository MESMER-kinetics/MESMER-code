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
#include <string>
#include <stdexcept>
#include "../System.h"
#include "../Spline.h"
#include "../Reaction.h"
#include "../Tunneling.h"

using namespace Constants;

namespace mesmer
{
  class DefinedTunnelingCoefficients : public TunnelingCalculator
  {
  public:

    ///Constructor which registers with the list of TunnellingCalculators in the base class
    DefinedTunnelingCoefficients(const char* id) : m_id(id),m_Energy(), m_PE(),m_VariationalThreshold(0.0){ Register(); }

    virtual ~DefinedTunnelingCoefficients() {}
    virtual const char* getID()  { return m_id; }
    virtual DefinedTunnelingCoefficients* Clone() { return new DefinedTunnelingCoefficients(*this); }

    // Read tunnelling coefficients. Called from Reaction.cpp ~line 274 via ParseForPlugin
    // and store in member variables. 
    virtual bool ParseData(PersistPtr pp);

    virtual bool calculateCellTunnelingCoeffs(Reaction* pReact, std::vector<double>& TunnelingProbability);

  private:

    const char* m_id;
    vector<double> m_Energy;
    vector<double> m_PE;
    double m_VariationalThreshold;
  };

  //************************************************************
  //Global instance, defining its id
  DefinedTunnelingCoefficients theDefinedTunnelingCoefficients("DefinedTunnelingCoefficient");
  //************************************************************

  bool DefinedTunnelingCoefficients::calculateCellTunnelingCoeffs(Reaction* pReact, vector<double>& TunnelingProbability){

    //TZ is the zpe of the TS.
    const double TZ = pReact->get_relative_TSZPE();
    //barrier0 & barrier1 are the zpe corrected barrier heights in the forward/reverse directions.
    int barrier0 = int(TZ - pReact->get_relative_rctZPE());
    int barrier1 = int(TZ - pReact->get_relative_pdtZPE());


    //Correct barrier1 and barrier0 according to maximum in vibrationally adiabatic curve
    int VaryCorrection = static_cast<int>(m_VariationalThreshold-barrier0);

    if(VaryCorrection > 0){
      barrier0 += VaryCorrection;
      barrier1 += VaryCorrection;
    }

    // Spline P(E)'s.

    Spline spline ;
    spline.Initialize(m_Energy, m_PE) ;

    // Set transmission coefficients to 0 where no tunneling is possible;
    // where tunneling may occur, the transmission coefficients are obtained
	  // from a spline fit to the user defined tunneling coefficients.
    const size_t MaximumCell = pReact->get_reactant()->getEnv().MaxCell;
    TunnelingProbability.clear();
    TunnelingProbability.resize(MaximumCell,0.0);
    for(size_t i = 0; i < TunnelingProbability.size() ; ++i){
      int E = int(i) - barrier0;
	  if (E <= 0) {
        TunnelingProbability[i] = spline.Calculate(double(i));
      } else if (E > 0 && E <= barrier0 ) { 
        //non classical reflection above the barrier
        TunnelingProbability[i] = 1.0 - spline.Calculate(double(2*barrier0-i)) ;
      } else if (E > barrier0) {
        TunnelingProbability[i] = 1.0 ;
      } 
    }

    // Writing out tunneling coefficients.

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
      m_VariationalThreshold = Ene;

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

