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
#include "../gWellProperties.h"
#include <cmath>
#include <cassert>

using namespace std;

namespace mesmer
{

  class ExponentialDown : public EnergyTransferModel
  {
  public:

    /********************************************************************************
    Constructor which registers this class with the list of plugins, initializes the
    ID and also does some initialization specific to this class.
    ********************************************************************************/
    ExponentialDown(const char* id) : m_id(id),
      m_deltaEDown(0.0, "cm-1"), m_refTemp(298), m_dEdExp(0.0), m_EdepParam(0.0)
    {
      Register();
    }

    virtual const char* getID() { return m_id; } //always required

    /********************************************************************************
    Optional description which will appear in a verbose listing of plugins.
    Put two spaces after \n  in Description() for good formatting in thelisting.
    ********************************************************************************/
    virtual const char* Description() {
      return
        "This calculates energy transfer probabilities on the basis of\n  "
        "the exponential down model: the probability of transition from a\n  "
        "grain of energy E to one of lower energy E1 is given by\n  "
        "  P(E|E1) = A(E1)exp(-(E - E1)/deltaEdown)\n  "
        "where A(E1) is a normalization factor. The probabilities of\n  "
        "activating collisions are found by detailed balance.\n"
        ;
    }

    /******************************************************************************
    Because the class can be used for more than one molecule, a new instance is made
    for each one. This is done by EnergyTransferModel::Find() calling Clone(); a
    function of the following form is required for each derived class.
    ******************************************************************************/
    ExponentialDown* Clone() { return new ExponentialDown(*this); }

    /*************************************************************
    Read the parameters needed by the class from the XML datafile
    *************************************************************/
    virtual bool ParseData(PersistPtr pp);
    //  virtual bool ReadParameters(Molecule* parent){return false;}//*****

      /*************************************************************
      This is the function which does the real work of the plugin
      *************************************************************/
    virtual double calculateTransitionProbability(double Ei, double Ej);

  private:
    const char* m_id; //all concrete plugin classes have this 
    Rdouble m_deltaEDown;
    double  m_refTemp;
    Rdouble m_dEdExp;
    Rdouble m_dEdAct;
    Rdouble m_EdepParam;
  };

  /******************************************************************************
  Declaration of a global instance of the class. This makes the class known to
  Mesmer and also defines its id "ExponentialDown".
  Names should be XML compatible ) so must not contain any of these characters ><"'&
  There can be additional global instances with different names (not usually necessary).
  ******************************************************************************/
  ExponentialDown exponentialDown("ExponentialDown");

  /******************************************************************************
  The energy transfer model for each modelled molecule is specified in the XML
  datafile by a child element of <molecule>:
  <me:energyTransferModel>ExponentialDown</me:energyTransferModel>
  If this is omitted, the default method specified in defaults.xml is used (which
  is currently "ExponentialDown", but could be changed).
  ******************************************************************************/

  /******************************************************************************
  Plugin classes usually read the data they require from the XML datafile and may
  store it themselves, although less specialised data may be stored in Molecule
  or Reaction.
  Returns the name of a special bath gas molecule or NULL.
  ******************************************************************************/
  bool ExponentialDown::ParseData(PersistPtr ppModel)
  {
    PersistPtr ppTop = m_parent->get_PersistentPointer();

    bool dataInProperty(true);
    PersistPtr ppPropList = ppTop->XmlMoveTo("propertyList");
    ppTop = ppPropList ? ppPropList : ppTop; //At <propertyList> if exists, else at <molecule>
    PersistPtr ppProp = ppTop->XmlMoveToProperty("me:deltaEDown");

    if (!ppProp)
    {
      //Data is a child of <me:energyTransferModel>
      dataInProperty = false;
      //PersistPtr ppModel = pp->XmlMoveTo("me:energyTransferModel");
      m_deltaEDown = ppModel->XmlReadDouble("me:deltaEDown"); //or use default
      ppProp = ppModel->XmlMoveTo("me:deltaEDown");
    }
    else //Try for data in a property
    {
      /******************************************************************************
      The following reads the content of every CML property "me:deltaEDown". If there
      is not one, the default value from defaults.xml is added to the internal XML tree
      and the value returned. This mechanism, which applies for most XmlRead
      operations unless there is a 'optional' parameter, is the recommended way to
      handle default values. It allows the default to be changed by the user, logs the
      use of the default, and provides error messages, including optional exhortations
      for the user to check the default (see the manual).
      ******************************************************************************/
      m_deltaEDown = ppTop->XmlReadPropertyDouble("me:deltaEDown");
    }
    if (IsNan(m_deltaEDown)) //unlikely failure
      return false;

    do //Loop over all <me:deltaEDown> or equivalent property.
    {
      /******************************************************************************
      The bath gas can optionally be specified (using ref or bathGas or omitted, when the
      general one specified directly under <me:conditions> is used. Each bath gas has its
      own instance of ExponentialDown (made, if necessary, in WellProperties::addBathGas()).
      ******************************************************************************/
      const char* bathGasName = NULL;
      bathGasName = ppProp->XmlReadValue("ref", optional);//ppProp is at <scalar> or <me:deltaEDown>
      if (!bathGasName)
        bathGasName = ppProp->XmlReadValue("bathGas", optional);
      ExponentialDown* pModel
        = static_cast<ExponentialDown*>(m_parent->getColl().addBathGas(bathGasName, this));
      assert(ppProp);
      istringstream ss(ppProp->XmlRead());
      ss >> pModel->m_deltaEDown;
      /******************************************************************************
      m_deltaEdown behaves most of the time like a normal variable of type double.
      But it can be a "range variable", taking a range of values when used in grid
      search and fitting routines. The me:deltaEDown property having both "lower" and
      "upper" attributes, together with and the following code, sets this up.
      ******************************************************************************/
      bool rangeSet;
      string varid = m_parent->getName() + ":deltaEDown";
      if (bathGasName)
        varid += string(":") += bathGasName;
      ReadRdoubleRange(varid, ppProp, pModel->m_deltaEDown, rangeSet);

    } while (ppProp = dataInProperty ? ppProp->XmlMoveToProperty("me:deltaEDown", true)
      : ppProp->XmlMoveTo("me:deltaEDown"));

    /******************************************************************************
    Read the temperature coefficients for all bath gases.
    The temperature dependence of <delta_E>d is accounted for as:
       <delta_E>d(T) = <delta_E>d_ref * (T / refTemp)^dEdExp * exp(-dEdAct/T)
    The default is hardwired at dEdExp = 0, dEdAct=0, so delta_E_down does not depend
    on temperature. The reference temperature of <delta_E>d, refTemp, also hardwired,
    has a default of 298K. Both of the temperature parameters can be range variables.
    ******************************************************************************/
    PersistPtr ppPropExp = dataInProperty ? ppTop : ppModel;
    while (ppPropExp = dataInProperty ?
      ppPropExp->XmlMoveToProperty("me:deltaEDownTExponent", true) :
      ppPropExp->XmlMoveTo("me:deltaEDownTExponent"))
    {
      m_refTemp = ppPropExp->XmlReadDouble("referenceTemperature", optional);
      if (IsNan(m_refTemp))
        m_refTemp = 298.;
      const char* bathGasName = ppPropExp->XmlReadValue("ref", optional);
      if (!bathGasName)
        bathGasName = ppPropExp->XmlReadValue("bathGas", optional);
      ExponentialDown* pModel
        = static_cast<ExponentialDown*>(m_parent->getColl().addBathGas(bathGasName, this));

      istringstream ss(ppPropExp->XmlRead());
      ss >> pModel->m_dEdExp;

      bool rangeSet;
      string varid = m_parent->getName() + ":deltaEDownTExponent";
      if (bathGasName)
        varid += string(":") += bathGasName;
      ReadRdoubleRange(varid, ppPropExp, pModel->m_dEdExp, rangeSet);
    }

    PersistPtr ppPropAct = dataInProperty ? ppTop : ppModel;
    while (ppPropAct = dataInProperty ?
      ppPropAct->XmlMoveToProperty("me:deltaEDownTActivation", true) :
      ppPropAct->XmlMoveTo("me:deltaEDownTActivation"))
    {
      const char* bathGasName = ppPropAct->XmlReadValue("ref", optional);
      if (!bathGasName)
        bathGasName = ppPropAct->XmlReadValue("bathGas", optional);
      ExponentialDown* pModel
        = static_cast<ExponentialDown*>(m_parent->getColl().addBathGas(bathGasName, this));

      istringstream ss(ppPropAct->XmlRead());
      ss >> pModel->m_dEdAct;

      bool rangeSet;
      string varid = m_parent->getName() + ":deltaEDownTActivation";
      if (bathGasName)
        varid += string(":") += bathGasName;
      ReadRdoubleRange(varid, ppPropAct, pModel->m_dEdAct, rangeSet);
    }

    /******************************************************************************
    Read in the energy dependence of <Delta E>d. A linear model is assumed for the 
    energy dependence, that is <Delta E>d(Ei) = <Delta E>d_ref + Ei * constant. The
    constant is hardwired to be 0.0 by default and can be fitted.
    ******************************************************************************/

    PersistPtr ppPropLinEne = dataInProperty ? ppTop : ppModel;
    while (ppPropLinEne = dataInProperty ?
      ppPropLinEne->XmlMoveToProperty("me:deltaEDownLinEne", true) :
      ppPropLinEne->XmlMoveTo("me:deltaEDownLinEne"))
    {
      const char* bathGasName = ppPropLinEne->XmlReadValue("ref", optional);
      if (!bathGasName)
        bathGasName = ppPropLinEne->XmlReadValue("bathGas", optional);
      ExponentialDown* pModel
        = static_cast<ExponentialDown*>(m_parent->getColl().addBathGas(bathGasName, this));

      istringstream ss(ppPropLinEne->XmlRead());
      ss >> pModel->m_EdepParam;

      bool rangeSet;
      string varid = m_parent->getName() + ":deltaEDownLinEne";
      if (bathGasName)
        varid += string(":") += bathGasName;
      ReadRdoubleRange(varid, ppPropLinEne, pModel->m_EdepParam, rangeSet);
    }

    return true;
  }

  /******************************************************************************
  This is the function which does the real work of the plugin
  ******************************************************************************/
  double ExponentialDown::calculateTransitionProbability(double Ei, double Ej) {
    double deltaEDown = m_deltaEDown + Ei * m_EdepParam;
    const double temperature = 1.0 / (boltzmann_RCpK * getParent()->getEnv().beta);
    if (m_dEdExp != 0.0)
      deltaEDown *= pow((temperature / m_refTemp), m_dEdExp);
    if (m_dEdAct != 0.0)
      deltaEDown *= exp(m_dEdAct * temperature);

    /******************************************************************************
    Issue a warning message if delta_E_down is smaller than grain size. When
    using cerr or cwarn, the message is displayed on the console and also added to
    the log file (usually mesmer.log). Using cinfo outputs to the log file only.
    Message may be prefixed by a "context" like "In R1:". Making an ErrorContext
    object sets the current context, and the previous context is restored when the
    object goes out of scope.
    The "once" manipulator causes the error message to be output only once if it is
    repeated (as this one would be). The message includes the context, so that if more
    than one molecule has an over-small delta_E_down there is a message for each molecule.
    ******************************************************************************/
    if (deltaEDown < double(getParent()->getEnv().GrainSize) && !getParent()->getFlags().allowSmallerDEDown) {
      ErrorContext e(getParent()->getName());
      cerr << "Delta E down is smaller than grain size: the solution may not converge." << once << endl;
    }
    return exp(-(Ei - Ej) / deltaEDown);
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
  /*
  Format for unified parsing
  <molecule>
    ...
    <me:energyTransferModel xsi:type="me:ExponentialDown">
      <me:deltaEDown units="cm-1" lower="160" upper="300" stepsize="10">190</me:deltaEDown>
      <me:deltaEDownTExponent referenceTemperature="300">1.0</me:deltaEDownTExponent>
      <me:deltaEDownTActivation units="K-1">1.2</me:deltaEDownTActivation> 
      <me:deltaEDown units="cm-1" bathGas="Ar">144.0</me:deltaEDown>
      <me:deltaEDownTExponent ref="Ar" lower="0.8" upper="1.2" stepsize="0.1">1.0</me:deltaEDownTExponent>
    </me:energyTransferModel>

  */

}//namespace

