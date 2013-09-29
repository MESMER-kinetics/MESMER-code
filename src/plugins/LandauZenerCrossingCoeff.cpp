//-------------------------------------------------------------------------------------------
//
// EckartCoefficients.h
//
// Author: Dave Glowacki, based on fortran code written by Jeremy Harvey
// Date:   28-8-2009
//
// Produces Landau-Zener spin forbidden crossing coefficients
//
//-------------------------------------------------------------------------------------------
#include <vector>
#include <cmath>
#include <string>
#include "../System.h"

using namespace Constants;

namespace mesmer
{
  class LandauZenerCrossingCoeff : public CrossingCalculator
  {
  public:
    virtual bool ParseData(PersistPtr pp);
    ///Constructor which registers with the list of CrossingCalculators in the base class
    LandauZenerCrossingCoeff(const char* id) : m_id(id),dataIsInTS(false){ Register(); }
  
    virtual ~LandauZenerCrossingCoeff() {}
    virtual const char* getID()  { return m_id; }
  
    virtual bool calculateCellCrossingCoeffs(Reaction* pReact, std::vector<double>& CrossingProbability);

    virtual bool ThereIsTunnellingWithCrossing(void) {return false;};

  private:
    bool ReadDoubleAndUnits(double& element, PersistPtr pp, const std::string identifier, const std::string units);
    bool ReadDoubleAndUnits2(double& element, PersistPtr pp, const std::string identifier, const std::string units);

  private:
    const char* m_id;
    double dataIsInTS;
    double m_SOCelement, m_GradDiffMagnitude, m_ReducedMass;
  };

  //************************************************************
  //Global instance, defining its id
  LandauZenerCrossingCoeff theLandauZenerCrossingCoeff("LZ");
  //************************************************************

  bool LandauZenerCrossingCoeff::ParseData(PersistPtr pp)
  {
    //Look for data in the reaction as child of <me:crossing>
    PersistPtr ppData = pp->XmlMoveTo("me:RMS_SOC_element");
    if(ppData)
    {
      return
       ReadDoubleAndUnits2(m_SOCelement, pp, "me:RMS_SOC_element", "cm-1") &&
       ReadDoubleAndUnits2(m_GradDiffMagnitude, pp, "me:GradientDifferenceMagnitude", "a.u./Bohr") &&
       ReadDoubleAndUnits2(m_ReducedMass, pp, "me:GradientReducedMass", "a.m.u.");
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

      cinfo << endl << "Spin Forbidden Crossing Data looked for in transition state molecule." << endl;

      if(!ReadDoubleAndUnits(m_SOCelement, ppPropList, "me:RMS_SOC_element", "cm-1")){
        return false;
      }

      if(!ReadDoubleAndUnits(m_GradDiffMagnitude, ppPropList, "me:GradientDifferenceMagnitude", "a.u./Bohr")){
        return false;
      }

      if(!ReadDoubleAndUnits(m_ReducedMass, ppPropList, "me:GradientReducedMass", "a.m.u.")){
        return false;
      }
      return true;
    }
  }

  bool LandauZenerCrossingCoeff::calculateCellCrossingCoeffs(Reaction* pReact, vector<double>& CrossingProbability){

    //Molecule * p_TransitionState = pReact->get_TransitionState();

    //// read input data for Landau-Zener crossing 

    //PersistPtr pp = p_TransitionState->get_PersistentPointer();
    //PersistPtr ppPropList = pp->XmlMoveTo("propertyList");
    //if(!ppPropList)
    //	ppPropList=pp; //Be forgiving; we can get by without a propertyList element

    //cinfo << endl << "Spin Forbidden Crossing Data for reaction " << pReact->getName() << endl;

    //double SOCelement(0.0), GradDiffMagnitude(0), ReducedMass(0);

    //if(!ReadDoubleAndUnits(SOCelement, ppPropList, "me:RMS_SOC_element", "cm-1")){
    //	return false;
    //}

    //if(!ReadDoubleAndUnits(GradDiffMagnitude, ppPropList, "me:GradientDifferenceMagnitude", "a.u./Bohr")){
    //	return false;
    //}

    //if(!ReadDoubleAndUnits(ReducedMass, ppPropList, "me:GradientReducedMass", "a.m.u.")){
    //	return false;
    //}

    // end read input data section

    // get threshold energy:
    double ZPE_corr_barrier_height;
    ZPE_corr_barrier_height  = pReact->get_ThresholdEnergy();

    const double SOCelementAU = m_SOCelement * (1/Hartree_in_RC);
    const double ReducedMassAU = m_ReducedMass * 1.822888e+3;

    //get properties of vectors in which to include Crossing coefficients
    const int MaximumCell = pReact->getEnv().MaxCell;
//		const int MaximumCell = pReactant->getEnv().MaxCell;
    CrossingProbability.clear();
    CrossingProbability.resize(MaximumCell);

    //set transmission coefficients to 0 below the ZPE corrected barrier height;
    //above the barrier, the Landau Zener transmission coefficients are calculated 
    //as described by Harvey & Aschi, Faraday Discuss, 2003 (124) 129-143

    for(int i = 0; i < MaximumCell; ++i){

      double E;
      if (ZPE_corr_barrier_height < 0.0){ 	//if the barrier height is negative
        E = double(i) + 0.5;                //set E to the avg energy of the cell in cm-1
      }
      else{                                 //otherwise if the barrier is gt.or.eq zero
        E = double(i) - ZPE_corr_barrier_height + 0.5;
      }

      if (E < 0){
        CrossingProbability[i] = 0.0;
      }
      else
      {
        double E_AU = E / Hartree_in_RC;
        double trans_probability = exp(-2.0e+0 * M_PI * pow(SOCelementAU,2.0e+0) / m_GradDiffMagnitude * pow(ReducedMassAU / (2.e+0 * E_AU),0.5e+0));
        CrossingProbability[i] = (1.0e+0 + trans_probability)*(1.0e+0 - trans_probability);
        // following if statement to avoid nan
        if(IsNan(CrossingProbability[i])) CrossingProbability[i] = 0.0;
      }
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

  bool LandauZenerCrossingCoeff::ReadDoubleAndUnits
    (double& element, PersistPtr pp, const std::string identifier, const std::string units){

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

  bool LandauZenerCrossingCoeff::ReadDoubleAndUnits2
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

