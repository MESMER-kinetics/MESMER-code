//-------------------------------------------------------------------------------------------
//
// GaussianModel.cpp
//
// This file contains the implementation of Gaussian energy transfer model.
//
// Author: David Glowacki (additional detail Struan Robertson).
// Date:   24/Nov/2012
//
//-------------------------------------------------------------------------------------------

#include "../EnergyTransferModel.h"
#include "../Rdouble.h"
#include "../Molecule.h"
#include "../gWellProperties.h"
#include <cmath>

using namespace std ;

namespace mesmer
{
  class Gaussian : public EnergyTransferModel
  {
  public:

    /********************************************************************************
    Constructor which registers this class with the map of energy transfer models
    kept by the base class.
    ********************************************************************************/
    Gaussian(const char* id) : m_id(id),
      m_GaussianCenter(0.0), m_GaussianWidth(0.0) { Register(); }

    virtual const char* getID()  { return m_id; }

    /******************************************************************************
    Because the class can be used for more than one molecule, a new instance is made
    for each one. This is done by EnergyTransferModel::Find() calling Clone(); a
    function of the following form is required for each derived class.
    ******************************************************************************/
    Gaussian* Clone() { return new Gaussian(*this); }

    /*************************************************************
    Read the parameters needed by the class from the XML datafile
    *************************************************************/
    virtual bool ParseData(PersistPtr pp) ;

    /*************************************************************
    This is the function which does the real work of the plugin
    *************************************************************/
    virtual double calculateTransitionProbability(double Ei, double Ej);

  private:
    const char* m_id;
    Rdouble m_GaussianCenter ;
    Rdouble m_GaussianWidth ;

  };
//*************************************
  Gaussian globalGaussian("gaussian");
//*************************************

  bool Gaussian::ParseData(PersistPtr pp) {

  bool rangeSet ;
  PersistPtr ppTop = m_parent->get_PersistentPointer();

  PersistPtr ppPropList = ppTop->XmlMoveTo("propertyList");
  ppTop = ppPropList ? ppPropList : ppTop; //At <propertyList> if exists, else at <molecule>
  PersistPtr ppProp = ppTop->XmlMoveToProperty("me:gaussianCenter"); 

  if (!ppProp)
  {
    //Data is a child of <me:energyTransferModel>
    PersistPtr ppCent = pp->XmlMoveTo("me:gaussianCenter");
    PersistPtr ppWidth = pp->XmlMoveTo("me:gaussianWidth");
    if(!ppCent || !ppWidth)
    {
      cerr << "Both <me:gaussianCenter> and <me:gaussianWidth> must be specified" << endl;
      return false;
    }
    const char* units = ppCent->XmlReadValue("units", optional);
    m_GaussianCenter = getConvertedEnergy(units ? units : "cm-1",
      pp->XmlReadDouble("me:gaussianCenter")); //read from defaults.xml if not present
    ReadRdoubleRange(string(m_parent->getName()+":gaussianCenter"), ppCent, m_GaussianCenter,
      rangeSet, getConvertedEnergy((units ? units : "cm-1"),1.0)) ;


    units = ppWidth->XmlReadValue("units");
    m_GaussianWidth = getConvertedEnergy(units ? units : "cm-1",
      pp->XmlReadDouble("me:gaussianWidth")); //read from defaults.xml if not present
    ReadRdoubleRange(string(m_parent->getName()+":gaussianWidth"), ppWidth, m_GaussianWidth,
      rangeSet, getConvertedEnergy((units ? units : "cm-1"),1.0)) ;
  }
  else
  {
    const char* units = ppProp->XmlReadPropertyAttribute("me:gaussianCenter", "units", optional);
    m_GaussianCenter = getConvertedEnergy(units ? units : "cm-1",
      ppTop->XmlReadPropertyDouble("me:gaussianCenter")); //read from defaults.xml if not present
    ReadRdoubleRange(string(m_parent->getName()+":gaussianCenter"), ppProp, m_GaussianCenter,
      rangeSet, getConvertedEnergy((units ? units : "cm-1"),1.0)) ;

    units = ppProp->XmlReadPropertyAttribute("me:gaussianWidth", "units", optional);
    m_GaussianWidth = getConvertedEnergy(units ? units : "cm-1",
      ppTop->XmlReadPropertyDouble("me:gaussianWidth")); //read from defaults.xml if not present
    ReadRdoubleRange(string(m_parent->getName()+":gaussianWidth"), ppProp, m_GaussianCenter,
      rangeSet, getConvertedEnergy((units ? units : "cm-1"),1.0)) ;
  }

  // ExponentialDown can use several bath gases with different energy transfer parameters.
  // This would also be possible for GaussianModel, but currently only the bathgas specified
  // in <me:conditions> is used.
  m_parent->getColl().addBathGas(NULL, this);

  // Issue a warning if the specified width is zero.
  if(m_GaussianWidth==0.0) {
    throw(std::runtime_error("me:gaussianWidth is equal to zero... cannot define the gaussian energy transfer function.")) ;
  }
  if (m_GaussianCenter == 0.0){
    // If m_GaussianCenter is at zero, issue a warning if the grain size is larger than the average 
  // downward energy transfer.
    double avgDownwardEnergyTransfer = 2.0*m_GaussianWidth/(sqrt(2.0*acos(-1.0)));  // this is the analytic average of a normalized gaussian centered at 0
    if( double(getParent()->getEnv().GrainSize) > avgDownwardEnergyTransfer){
      cerr << "For a gaussian centered at zero, the average downward energy transfer is" << avgDownwardEnergyTransfer << " which is larger than the grain size..." << endl;
      cerr << "Check your input data." << endl;
    }
  } else {    // issue a warning if the grain size is larger than the average downward energy transfer 
    double avgDownwardEnergyTransfer = m_GaussianCenter;  // this is the average downward energy transferred for a normalized gaussian centered at 0
    if( double(getParent()->getEnv().GrainSize) > avgDownwardEnergyTransfer){
      cerr << "for a gaussian centered at " << m_GaussianCenter << ", the average downward energy transfer is" << avgDownwardEnergyTransfer << " which is larger than the grain size..." << endl;
      cerr << "Check your input data." << endl;
    }
  }

    return true; //TODO Parse for ref attribute, when this might be the name of a bath gas
  }

  /******************************************************************************
  This is the function which does the real work of the plug-in.
  ******************************************************************************/
  double Gaussian::calculateTransitionProbability(double Ej, double Ei) {

	// The variable m_GaussianCenter makes larger energy transfers more likely
	// than small ones. However, this has the effect that transition probabilities
	// from states at energies close to zero are very high, and this can be 
	// incompatabile with the Boltzmann distribution. The manifestation of this
	// effect is the negative normalization coefficitents. To attenuate this effect,
	// a factor has been introduced, that is a function of the initial energy, which
	// gradually increases the effect of m_GaussianCenter.

	double ratio(0.0) ;
	double spread(4.0) ;
	if (m_GaussianCenter > 0.0) {
	  ratio = Ej/(spread*m_GaussianCenter) ;
	}

    double dE = (Ej - Ei - m_GaussianCenter * (ratio/(1.0 + ratio)))/m_GaussianWidth ;

    return exp(-0.5*dE*dE);
  }

}//namespace

