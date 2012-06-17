/********************************************************************************
This file is working code for a class which implements the exponential down model
for energy transfer during a collision, but it is heavily commented so as to be a
tutorial on the writing of classes which implement different energy transfer
models, or, more widely, models for other parts of the calculation. The parent
class in EnergyTransferModel.h is also commented.

These classes are plugin classes - they can be added or removed without having to
change any of the existing code. They derived from a base class (here it is
EnergyTransferModel) which is usually abstract. It has virtual functions which
are called by the main code and are probably redefined in the derived classes
like this one.

A plugin class can be declared in a .h file and defined in a .cpp file, as normal.
But the header file is unlikely ever to be used independently, so it may be
convenient for the declaration to be also in the .cpp file, so that the plugin is
all in one file, as here.
********************************************************************************/

#include "../EnergyTransferModel.h"
#include "../Rdouble.h"
#include "../MolecularComponents.h"
#include "../Molecule.h"
#include <cmath>

using namespace std ;

namespace mesmer
{
  class gaussian : public EnergyTransferModel
  {
  public:

    /********************************************************************************
    Constructor which registers this class with the map of energy transfer models
    kept by the base class.
    ********************************************************************************/
    gaussian(const std::string& id) : EnergyTransferModel(id),
      m_gaussianCenter(0.0), m_gaussianWidth(0.0) {}

    /******************************************************************************
    Because the class can be used for more than one molecule, a new instance is made
    for each one. This is done by EnergyTransferModel::Find() calling Clone(); a
    function of the following form is required for each derived class.
    ******************************************************************************/
    gaussian* Clone() { return new gaussian(*this); }

    /*************************************************************
    Read the parameters needed by the class from the XML datafile
    *************************************************************/
    virtual bool ReadParameters(const Molecule* parent) ; 

    /*************************************************************
    This is the function which does the real work of the plugin
    *************************************************************/
    virtual double calculateTransitionProbability(double Ei, double Ej);

  private:
    Rdouble m_gaussianCenter ;
    Rdouble m_gaussianWidth ;

  };

  /******************************************************************************
  Declaration of a global instance of the class. This makes the class known to
  Mesmer and also defines its id "ExponentialDown".
  ******************************************************************************/
  gaussian gaussian("gaussian");

  /******************************************************************************
  The energy transfer model for each modelled molecule is specified in the XML
  datafile by a child element of <molecule>:
  <me:energyTransferModel>gaussian</me:energyTransferModel>
  If this is omitted, the default method specified in defaults.xml is used (which
  is currently "ExponentialDown", but could be changed).
  ******************************************************************************/

  /******************************************************************************
  Plugin classes usually read the data they require from the XML datafile and may
  store it themselves, although less specialised data may be stored in Molecule
  or Reaction.
  ******************************************************************************/
  bool gaussian::ReadParameters(const Molecule* parent) { 

    setParent(parent);
    PersistPtr pp = parent->get_PersistentPointer();
    // There may or may not be a <propertyList> element. If not, the <property>
    //  elements are children of <molecule>
    PersistPtr ppPropList = pp->XmlMoveTo("propertyList");
    if(!ppPropList)
      ppPropList=pp;

    /******************************************************************************
    The following statement reads the value of the CML property "me:gaussianCenter". If
    it is not present the default value from defaults.xml is added to the internal
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
    double value(0.0);
    idata >> value;
    m_gaussianCenter = value;

    /******************************************************************************
    m_gaussianCenter behaves most of the time like a normal variable of type double.
    But it can be a "range variable", taking a range of values when used in grid
    search and fitting routines. The me:deltaEDown property having both "lower" and
    "upper" attributes, together with and the following code, sets this up.
    ******************************************************************************/
    // Needed to read the attributes.
    bool rangeSet ;
    PersistPtr ppProp = ppPropList->XmlMoveToProperty("me:gaussianCenter"); 
    ReadRdoubleRange(string(parent->getName()+":gaussianCenter"), ppProp, m_gaussianCenter, rangeSet) ;

    // By default, me:gaussianCenter = 0.0 
    //             me:gaussianWidth = 0.0

    txt = ppPropList->XmlReadProperty("me:gaussianWidth"); //required in datafile or defaults.xml
    if(!txt)
      return false;
    istringstream iidata(txt);
    value = 0.0;
    iidata >> value;
    m_gaussianWidth = value;

    PersistPtr ppPropExp = ppPropList->XmlMoveToProperty("me:gaussianWidth"); 
    ReadRdoubleRange(string(parent->getName()+":gaussianWidth"), ppPropExp, m_gaussianWidth, rangeSet) ;

    // issue some warning messages
		double gaussianCenter = m_gaussianCenter;
    double gaussianWidth = m_gaussianWidth;

		// issue a warning if the specified width is zero
    if(gaussianWidth==0.0) {
			cerr << "me:gaussianWidth is equal to zero... cannot define the gaussian energy transfer function" << endl;
			exit(1);
    }
    // if m_gaussianCenter is at zero, issue a warning if the grain size is larger than the average downward energy transfer
		if (gaussianCenter == 0.0){
			double avgDownwardEnergyTransfer = 2.0*gaussianWidth/(sqrt(2.0*acos(-1.0)));  // this is the analytic average of a normalized gaussian centered at 0
			if( double(getParent()->getEnv().GrainSize) > avgDownwardEnergyTransfer){
				cerr << "for a gaussian centered at zero, the average downward energy transfer is" << avgDownwardEnergyTransfer << " which is larger than the grain size..." << endl;
				cerr << "check your input and try again... exiting." << endl;
				exit(1);
			}
		}
		else{  // issue a warning if the grain size is larger than the average downward energy transfer 
			double avgDownwardEnergyTransfer = gaussianCenter;  // this is the average downward energy transferred for a normalized gaussian centered at 0
			if( double(getParent()->getEnv().GrainSize) > avgDownwardEnergyTransfer){
				cerr << "for a gaussian centered at " << gaussianCenter << ", the average downward energy transfer is" << avgDownwardEnergyTransfer << " which is larger than the grain size..." << endl;
				cerr << "check your input and try again... exiting." << endl;
				exit(1);
			}
		}

    return true ; 
  }
  /******************************************************************************
  This is the function which does the real work of the plugin
  ******************************************************************************/
  double gaussian::calculateTransitionProbability(double Ei, double Ej) {
    // return exp(-(Ei - Ej - gaussianCenter)^2/(2 * gaussianWidth)) ;

		double gaussianCenter = m_gaussianCenter;
    double gaussianWidth = m_gaussianWidth;
		double dE = Ei - Ej;

//		cout << dE << " " << exp(-0.5*(dE - gaussianCenter)*(dE - gaussianCenter) / (gaussianWidth*gaussianWidth) ) << endl;

		return exp(-0.5*(dE - gaussianCenter)*(dE - gaussianCenter) / (gaussianWidth*gaussianWidth) );
  }

  /******************************************************************************
  In summary, to make a plugin class with a new energy transfer model:
  - Copy this file, changing all the "ExponentialDown" to your class name.
  - Change the code in calculateTransitionProbability() to your model.
  - Alter ReadParameters() to input any data your model needs from the XML file.
  Any data that is essential should preferably have an entry in defaults.xml,
  The "default" attribute of this should be in uppercase if it is necessary
  for the user to review the value of the default.
  If the data does not need to be provided, add an optional parameter to the
  call to an XmlRead function.
  - In the XML file, add or change
  <me:energyTransferModel>yourID</me:energyTransferModel>
  Add your data, which should usually have an me: prefix.

  ******************************************************************************/

}//namespace