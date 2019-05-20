//-------------------------------------------------------------------------------------------
//
// HirstWardlawPotential.cpp
//
// Author: Struan Robertson
// Date:   29/Apr/2019
//
// Calculates the Hirst-Wardlaw potential for the CH3 + H --> CH4 reaction.
//
//-------------------------------------------------------------------------------------------
#include "../FTSTPotential.h"
#include "../MesmerMath.h"
#include "../Constants.h"
#include "../unitsConversion.h"

using namespace std;
using namespace Constants;

namespace mesmer
{

  class HirstWardlawPotential : public FTSTPotential {

  public:

    HirstWardlawPotential(const char* id) : m_id(id), 
      m_V0(0.0),
      m_ant(0.0),
      m_bnt(0.0),
      m_threshold(0.0) { Register(); }

    virtual ~HirstWardlawPotential() {}

    virtual const char* getID() { return m_id; }

    virtual HirstWardlawPotential* Clone() { return new HirstWardlawPotential(*this); }

    void Initialize() {
      m_threshold = getConvertedEnergy("kcal/mol", 103.432);
      m_ant = getConvertedEnergy("kcal/mol", 164.0);
      m_bnt = 0.43;
    }

    virtual void RxnCrdInitialize(double rxnCrd) {
      m_V0 = m_ant * exp(-m_bnt * rxnCrd*rxnCrd);
    }


    virtual double MEPPotential(double rxnCrd) {
      double potential(0.0);

      double A3 = -1.017234e-03; // Hirst Surface parameters as fit by Hase et al.
      double A2 = 7.7738886e-02; // A3, A2, A1, AO parameters of the variable
      double A1 = 7.703640e-02;  // beta of the Modified Morse function.
      double A0 = 1.686690e0;    // 

      double delR  = rxnCrd - 1.0015;
      double delR2 = delR * delR;
      double delR3 = delR * delR2;

      double bb  = A3 * delR3 + A2 * delR2 + A1 * delR + A0;
      double tmp = 1.0 - exp(-bb * delR);
      potential = m_threshold * tmp*tmp;

      return potential;

    }

    virtual double HinderingPotential(double rxnCrd, const vector<double>& angles) {
      double tmp = sin(angles[0]);
      return  m_V0 * tmp * tmp;
    }

  private:

    const char* m_id;

    double m_V0;
    double m_ant;
    double m_bnt;
    double m_threshold;

  };

  //-------------------------------------------------------------
  //Global instance, defining its id
  HirstWardlawPotential theHirstWardlawPotential("HirstWardlawPotential");
  //-------------------------------------------------------------


}
