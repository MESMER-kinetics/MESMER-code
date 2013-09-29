//-------------------------------------------------------------------------------------------
//
// EckartCoefficients.h
//
// Author: Dave Glowacki, based on fortran code written by Jeremy Harvey
// Date:   28-8-2009
//
// Produces WKB spin forbidden crossing coefficients
//
//-------------------------------------------------------------------------------------------
#include <vector>
#include <cmath>
#include <string>
#include "../System.h"

using namespace Constants;

namespace mesmer
{
  class WKBCrossingCoeff : public CrossingCalculator
  {
  public:
    virtual bool ParseData(PersistPtr pp);

    ///Constructor which registers with the list of CrossingCalculators in the base class
    WKBCrossingCoeff(const char* id) : m_id(id){ Register(); }

    virtual ~WKBCrossingCoeff() {}
    virtual const char* getID()  { return m_id; }

    virtual bool calculateCellCrossingCoeffs(Reaction* pReact, std::vector<double>& CrossingProbability);

    virtual bool ThereIsTunnellingWithCrossing(void) {return true;};

  private:
    bool ReadDoubleAndUnits(double& element, PersistPtr pp, const std::string identifier, const std::string units);
    bool ReadDoubleAndUnits2(double& element, PersistPtr pp, const std::string identifier, const std::string units);

  private:
    const char* m_id;
    double dataIsInTS;
    double m_SOCelement, m_GradDiffMagnitude, m_ReducedMass, m_AverageSlope;
  };

  //************************************************************
  //Global instance, defining its id 
  WKBCrossingCoeff theWKBCrossingCoeff("WKB");
  //************************************************************

  bool WKBCrossingCoeff::ParseData(PersistPtr pp)
  {
    //Look for data in the reaction as child of <me:crossing>
    PersistPtr ppData = pp->XmlMoveTo("me:RMS_SOC_element");
    if(ppData)
    {
      return
       ReadDoubleAndUnits2(m_SOCelement, pp, "me:RMS_SOC_element", "cm-1") &&
       ReadDoubleAndUnits2(m_GradDiffMagnitude, pp, "me:GradientDifferenceMagnitude", "a.u./Bohr") &&
       ReadDoubleAndUnits2(m_ReducedMass, pp, "me:GradientReducedMass", "a.m.u.") &&
       ReadDoubleAndUnits2(m_AverageSlope, pp, "me:AverageSlope", "a.u./Bohr");
    }
    else
    {
      // data is in TS
      pp = getParent()->get_TransitionState()->get_PersistentPointer();
      if(!pp)
        return false;
      PersistPtr ppPropList = pp->XmlMoveTo("propertyList");
      if(!ppPropList)
        ppPropList=pp; // we can get by without a propertyList element

      if(!ReadDoubleAndUnits(m_SOCelement, ppPropList, "me:RMS_SOC_element", "cm-1")){
        return false;
      }

      if(!ReadDoubleAndUnits(m_GradDiffMagnitude, ppPropList, "me:GradientDifferenceMagnitude", "a.u./Bohr")){
        return false;
      }

      if(!ReadDoubleAndUnits(m_ReducedMass, ppPropList, "me:GradientReducedMass", "a.m.u.")){
        return false;
      }
      if(!ReadDoubleAndUnits(m_AverageSlope, ppPropList, "me:AverageSlope", "a.u./Bohr")){
        return false;
      }
      return true;
    }
  }

  bool WKBCrossingCoeff::calculateCellCrossingCoeffs(Reaction* pReact, vector<double>& CrossingProbability)
  {
    double ZPE_corr_barrier_height;

    // get threshold energy:
    ZPE_corr_barrier_height  = pReact->get_ThresholdEnergy();

    const double SOCelementAU = m_SOCelement * (1/Hartree_in_RC);
    const double ReducedMassAU = m_ReducedMass * 1.822888e+3;

    //get properties of vectors in which to include Crossing coefficients
    const int MaximumCell = pReact->getEnv().MaxCell;
    CrossingProbability.clear();
    CrossingProbability.resize(MaximumCell);

    //set transmission coefficients to 0 below the ZPE corrected barrier height;
    //above the barrier, the Landau Zener transmission coefficients are calculated 
    //as described by Harvey & Aschi, Faraday Discuss, 2003 (124) 129-143

    double Ai(0);

    const double DT_1 = 4.0e+0 * pow(M_PI,2.0e+0) * pow(SOCelementAU,2.0e+0) * pow((2.0e+0 * ReducedMassAU / (m_GradDiffMagnitude * m_AverageSlope)),(2.0e+0/3.0e+0));
    const double DT_2 = pow((2.0e+0 * ReducedMassAU * pow(m_GradDiffMagnitude,2.0e+0) / pow(m_AverageSlope,4.0e+0)),(1.0e+0/3.0e+0));

    for(int i = 0; i < MaximumCell; ++i){
      double E = double(i) - ZPE_corr_barrier_height;
      double E_AU = E / Hartree_in_RC;
      double xvalue = -1.0e+0* E_AU * DT_2;
      double Aip(0), Bi(0), Bip(0); //airy returns Ai, Bi, and their derivatives, Aip & Bip
      airy2(xvalue, Ai);            //airy2 accurate over entire tunnelling regime
      //airy(xvalue, Ai, Aip, Bi, Bip); //airy accurate in shallow tunnelling regime; less so for deep tunnelling
      CrossingProbability[i] = DT_1 * pow(Ai,2.0e+0);
      // following if statement to avoid nan
      if(IsNan(CrossingProbability[i])) CrossingProbability[i] = 0.0;
    }

    if (pReact->getFlags().CrossingCoeffEnabled){
      ctest << "\nCrossing coefficients for: " << pReact->getName() << endl;
      for(int i = 0; i < MaximumCell; ++i){
        ctest << CrossingProbability[i] << endl;
      }
      ctest << "}\n";
    }

    return true;
  }

  bool WKBCrossingCoeff::ReadDoubleAndUnits(double& element, PersistPtr pp, const std::string identifier, const std::string units){

    element=pp->XmlReadPropertyDouble(identifier,true);
    string unitsTxt;

    if(element){
      unitsTxt = pp->XmlReadPropertyAttribute(identifier, "units");
      if (unitsTxt!=units){
        cinfo << "MESMER could not read units for " << identifier << "; assuming " << units << "." << endl;
      }
      cinfo << identifier << " = " << element << " " << units << endl;
    }  
    else{
      cerr << "Spin forbidden crossing: failed to read " << identifier << " (" << units << ")." << endl; 
      return false;
    }
    return true;
  }

  bool WKBCrossingCoeff::ReadDoubleAndUnits2
      (double& element, PersistPtr pp, const std::string identifier, const std::string units)
  {
    PersistPtr ppData = pp->XmlMoveTo(identifier);
    if(!ppData)
      return false;
		element = pp->XmlReadDouble(identifier); //or default
    if(!IsNan(element))
    {
      string unitsTxt = ppData->XmlReadValue("units", false);
		  if (unitsTxt!=units){
			  cinfo << "MESMER could not read units for " << identifier << "; assuming " << units << "." << endl;
  	  cinfo << identifier << " = " << element << " " << units << endl;
      }
    }
		else
    {
			cerr << "Spin forbidden crossing: failed to read " << identifier << " (" << units << ")." << endl; 
			return false;
		}
		return true;
	}

}//namespace

