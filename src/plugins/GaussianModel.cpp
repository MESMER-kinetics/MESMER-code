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
#include "../MolecularComponents.h"
#include "../Molecule.h"
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
    virtual bool ReadParameters(Molecule* parent) ; 

    /*************************************************************
    This is the function which does the real work of the plugin
    *************************************************************/
    virtual double calculateTransitionProbability(double Ei, double Ej);

  private:
    const char* m_id;
    Rdouble m_GaussianCenter ;
    Rdouble m_GaussianWidth ;

  };

  /******************************************************************************
  Declaration of a global instance of the class. This makes the class known to
  Mesmer and also defines its id "Gaussian".
  ******************************************************************************/
  Gaussian globalGaussian("gaussian");

  /******************************************************************************
  The energy transfer model for each modelled molecule is specified in the XML
  datafile by a child element of <molecule>:
  <me:energyTransferModel>gaussian</me:energyTransferModel>
  If this is omitted, the default method specified in defaults.xml is used (which
  is currently "ExponentialDown", but could be changed).
  ******************************************************************************/

  /******************************************************************************
  Plugin classes usually read the data they require from the XML datafile and may
  store it in Molecule or Reaction if this is sensible or in the plugin class for
  more specialized data..
  ******************************************************************************/
  bool Gaussian::ReadParameters(Molecule* parent) { 

    setParent(parent);
    PersistPtr pp = parent->get_PersistentPointer();
    PersistPtr ppPropList = pp->XmlMoveTo("propertyList");
    if(!ppPropList)
      ppPropList=pp;

    /******************************************************************************
    The following statement reads the value of the CML property "me:GaussianCenter".
    If it is not present the default value from defaults.xml is added to the internal
    XML tree and the value returned. This mechanism, which applies for most XmlRead
    operations unless there is a 'optional' parameter, is the recommended way to
    handle default values. It allows the default to be changed by the user, logs the
    use of the default, and provides error messages, including optional exhortations
    for the user to check the default (see the manual).
    ******************************************************************************/
    const char* txt = ppPropList->XmlReadProperty("me:gaussianCenter"); //required
    if(!txt)
      return false;
    istringstream idata(txt);
    idata >> m_GaussianCenter ;

    // Needed to read the attributes.
    bool rangeSet ;
    PersistPtr ppProp = ppPropList->XmlMoveToProperty("me:gaussianCenter"); 
    ReadRdoubleRange(string(parent->getName()+":gaussianCenter"), ppProp, m_GaussianCenter, rangeSet) ;

    txt = ppPropList->XmlReadProperty("me:gaussianWidth"); //required in datafile or defaults.xml
    if(!txt)
      return false;
    istringstream iidata(txt);
    iidata >> m_GaussianWidth ;

    PersistPtr ppPropExp = ppPropList->XmlMoveToProperty("me:gaussianWidth"); 
    ReadRdoubleRange(string(parent->getName()+":gaussianWidth"), ppPropExp, m_GaussianWidth, rangeSet) ;

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
    } else {  
	  // issue a warning if the grain size is larger than the average downward energy transfer 
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