//-------------------------------------------------------------------------------------------
//
// BiExponentialDown.cpp
//
// This file contains the implementation of Biexponential Down energy transfer model. This
// implementation differs from the simple exponential down class in that it will be assumed
// that there will be NO definition of parameters in the proprety list, i.e. parameters will
// only be defined as part of the model. 
//
// Author: Struan Robertson.
// Date:   27/Oct/2020
//
//-------------------------------------------------------------------------------------------

#include "../EnergyTransferModel.h"
#include "../Rdouble.h"
#include "../Molecule.h"
#include "../gWellProperties.h"
#include <cmath>
#include <cassert>

using namespace std;

namespace mesmer
{

  class BiExponentialDown : public EnergyTransferModel
  {
  public:

    BiExponentialDown(const char* id) : m_id(id),
      m_deltaEDown(200.0, "cm-1"), m_deltaEDown2(1000.0, "cm-1"), m_ratio(0.0), m_refTemp(298), m_dEdExp(0.0), m_dEdAct(0.0)
    {
      Register();
    }

    virtual const char* getID() { return m_id; } //always required

    virtual const char* Description() {
      return
        "This calculates energy transfer probabilities on the basis of\n  "
        "the biexponential down model: the probability of transition from a\n  "
        "grain of energy E to one of lower energy E1 is given by\n  "
        "P(E|E1) = A(E1)((1-y)exp(-(E - E1)/<DE>d1) + y exp(-(E - E1)/<DE>d2))\n  "
        "where A(E1) is a normalization factor. The probabilities of\n  "
        "activating collisions are found by detailed balance.\n"
        ;
    }

    /******************************************************************************
    Because the class can be used for more than one molecule, a new instance is made
    for each one. This is done by EnergyTransferModel::Find() calling Clone(); a
    function of the following form is required for each derived class.
    ******************************************************************************/
    BiExponentialDown* Clone() { return new BiExponentialDown(*this); }

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
    Rdouble m_deltaEDown2;
    Rdouble m_ratio;
    double  m_refTemp;
    Rdouble m_dEdExp;
    Rdouble m_dEdAct;
  };

  /******************************************************************************
  Declaration of a global instance of the class. This makes the class known to
  Mesmer and also defines its id "BiExponentialDown".
  Names should be XML compatible ) so must not contain any of these characters ><"'&
  There can be additional global instances with different names (not usually necessary).
  ******************************************************************************/
  BiExponentialDown biExponentialDown("BiExponentialDown");

  bool BiExponentialDown::ParseData(PersistPtr ppModel)
  {
    //Data is a child of <me:energyTransferModel>
    m_deltaEDown = ppModel->XmlReadDouble("me:deltaEDown"); //or use default
    m_deltaEDown2 = ppModel->XmlReadDouble("me:deltaEDown2"); //or use default
    m_ratio = ppModel->XmlReadDouble("me:ratio"); //or use default

    if (IsNan(m_deltaEDown) || IsNan(m_deltaEDown2) || IsNan(m_ratio)) //unlikely failure
      return false;

    PersistPtr ppProp = ppModel->XmlMoveTo("me:deltaEDown");
    do //Loop over all <me:deltaEDown> or equivalent property.
    {
      /******************************************************************************
      The bath gas can optionally be specified (using ref or bathGas or omitted, when
      the general one specified directly under <me:conditions> is used). Each bath gas
      has its own instance of BiExponentialDown (made, if necessary, in
      WellProperties::addBathGas()).
      ******************************************************************************/
      const char* bathGasName = NULL;
      bathGasName = ppProp->XmlReadValue("ref", optional);//ppProp is at <scalar> or <me:deltaEDown>
      if (!bathGasName)
        bathGasName = ppProp->XmlReadValue("bathGas", optional);
      BiExponentialDown* pModel
        = static_cast<BiExponentialDown*>(m_parent->getColl().addBathGas(bathGasName, this));
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

    } while (ppProp = ppProp->XmlMoveTo("me:deltaEDown"));

    // Loop over all <me:deltaEDown2>.

    PersistPtr pp = ppModel;
    while (pp = pp->XmlMoveTo("me:deltaEDown2"))
    {
      const char* bathGasName = pp->XmlReadValue("ref", optional);
      if (!bathGasName)
        bathGasName = pp->XmlReadValue("bathGas", optional);
      BiExponentialDown* pModel
        = static_cast<BiExponentialDown*>(m_parent->getColl().addBathGas(bathGasName, this));

      istringstream ss(pp->XmlRead());
      ss >> pModel->m_deltaEDown2;

      bool rangeSet;
      string varid = m_parent->getName() + ":deltaEDown2";
      if (bathGasName)
        varid += string(":") += bathGasName;
      ReadRdoubleRange(varid, pp, pModel->m_deltaEDown2, rangeSet);
    }

    // Loop over all <me:ratio>.

    pp = ppModel;
    while (pp = pp->XmlMoveTo("me:ratio"))
    {
      const char* bathGasName = pp->XmlReadValue("ref", optional);
      if (!bathGasName)
        bathGasName = pp->XmlReadValue("bathGas", optional);
      BiExponentialDown* pModel
        = static_cast<BiExponentialDown*>(m_parent->getColl().addBathGas(bathGasName, this));

      istringstream ss(pp->XmlRead());
      ss >> pModel->m_ratio;

      bool rangeSet;
      string varid = m_parent->getName() + ":ratio";
      if (bathGasName)
        varid += string(":") += bathGasName;
      ReadRdoubleRange(varid, pp, pModel->m_ratio, rangeSet);
    }

    /******************************************************************************
    Read the temperature coefficients for all bath gases.
    The temperature dependence of <delta_E_down> is accounted for as:
       <delta_E_down>(T) = <delta_E_down>_ref * (T / refTemp)^dEdExp * exp(-dEdAct/T)
    The default is hardwired at dEdExp = 0, dEdAct=0, so delta_E_down does not depend
    on temperature. The reference temperature of <DeltaEDown>, refTemp, also hardwired,
    has a default of 298K. Both of the temperature parameters can be range variables.
    ******************************************************************************/
    PersistPtr ppPropExp = ppModel;
    while (ppPropExp = ppPropExp->XmlMoveTo("me:deltaEDownTExponent"))
    {
      m_refTemp = ppPropExp->XmlReadDouble("referenceTemperature", optional);
      if (IsNan(m_refTemp))
        m_refTemp = 298.;
      const char* bathGasName = ppPropExp->XmlReadValue("ref", optional);
      if (!bathGasName)
        bathGasName = ppPropExp->XmlReadValue("bathGas", optional);
      BiExponentialDown* pModel
        = static_cast<BiExponentialDown*>(m_parent->getColl().addBathGas(bathGasName, this));

      istringstream ss(ppPropExp->XmlRead());
      ss >> pModel->m_dEdExp;

      bool rangeSet;
      string varid = m_parent->getName() + ":deltaEDownTExponent";
      if (bathGasName)
        varid += string(":") += bathGasName;
      ReadRdoubleRange(varid, ppPropExp, pModel->m_dEdExp, rangeSet);
    }

    PersistPtr ppPropAct = ppModel;
    while (ppPropAct = ppPropAct->XmlMoveTo("me:deltaEDownTActivation"))
    {
      const char* bathGasName = ppPropAct->XmlReadValue("ref", optional);
      if (!bathGasName)
        bathGasName = ppPropAct->XmlReadValue("bathGas", optional);
      BiExponentialDown* pModel
        = static_cast<BiExponentialDown*>(m_parent->getColl().addBathGas(bathGasName, this));

      istringstream ss(ppPropAct->XmlRead());
      ss >> pModel->m_dEdAct;

      bool rangeSet;
      string varid = m_parent->getName() + ":deltaEDownTActivation";
      if (bathGasName)
        varid += string(":") += bathGasName;
      ReadRdoubleRange(varid, ppPropAct, pModel->m_dEdAct, rangeSet);
    }

    return true;
  }

  /******************************************************************************
  This is the function which does the real work of the plugin
  ******************************************************************************/
  double BiExponentialDown::calculateTransitionProbability(double Ei, double Ej) {
    double deltaEDown = m_deltaEDown;
    double deltaEDown2 = m_deltaEDown2;
    const double temperature = 1.0 / (boltzmann_RCpK * getParent()->getEnv().beta);
    if (m_dEdExp != 0.0)
      deltaEDown = deltaEDown * pow((temperature / m_refTemp), m_dEdExp);
    if (m_dEdAct != 0)
      deltaEDown = deltaEDown * exp(m_dEdAct * temperature);

    /******************************************************************************
    Issue a warning message if delta_E_down is smaller than grain size.
    When using cerr or cwarn, the message is displayed on the console and also added to
    the log file (usually mesmer.log). Using cinfo outputs to the log file only.
    Message may be prefixed by a "context" like "In R1:". Making an ErrorContext object
    sets the current context, and the previous context is restored when the object goes
    out of scope.
    The "once" manipulator causes the error message to be output only once if it is
    repeated (as this one would be). The message includes the context, so that if more
    than one molecule has an over-small delta_E_down there is a message for each molecule.
    ******************************************************************************/
    if (deltaEDown < double(getParent()->getEnv().GrainSize) && !getParent()->getFlags().allowSmallerDEDown) {
      ErrorContext e(getParent()->getName());
      cerr << "Delta E down is smaller than grain size: the solution may not converge." << once << endl;
    }

    double ratio = Ei / 1000.0;
    deltaEDown *= ratio / (1 + ratio);

    return (1.0 - m_ratio) * exp(-(Ei - Ej) / deltaEDown) + m_ratio * exp(-(Ei - Ej) / deltaEDown2);
  }

}//namespace

