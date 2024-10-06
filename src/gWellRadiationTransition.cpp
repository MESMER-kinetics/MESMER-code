//-------------------------------------------------------------------------------------------
//
// gWellRadiationTransition.cpp
//
// Author: Struan Robertson
// Date:   27/Dec/2021
//
// This file contains the implementation of the gWellRadiationTransition class. This 
// class governs radiation transitions between states.
//
//-------------------------------------------------------------------------------------------
#include "Molecule.h"
#include "gWellRadiationTransition.h"
#include "System.h"

using namespace std;
using namespace Constants;

namespace mesmer
{
  //
  // Constructor, destructor and initialization
  //
  gWellRadiationTransition::gWellRadiationTransition(Molecule* pMol) : MolecularComponent(), 
    m_TransitionFrequency(),
    m_EinsteinAij(),
    m_EinsteinBij(),
    m_bActivation(false)
  {
    ErrorContext c(pMol->getName());
    m_host = pMol;
  }

  gWellRadiationTransition::~gWellRadiationTransition()
  {
  }

  bool gWellRadiationTransition::initialization() {

    PersistPtr pp = m_host->get_PersistentPointer();

    // Read in Einstein Aij and/or Bij coefficents.

    PersistPtr ppPropList = pp->XmlMoveTo("propertyList");
    if (!ppPropList)
      ppPropList = pp; // A propertyList element is not essential.

    const char* txt;
    if ((txt = ppPropList->XmlReadProperty("me:EinsteinAij", optional))) {
      istringstream idata(txt);
      double x;
      while (idata >> x)
        m_EinsteinAij.push_back(x);
    }

    if ((txt = ppPropList->XmlReadProperty("me:EinsteinBij", optional))) {
      istringstream idata(txt);
      double x;
      while (idata >> x)
        m_EinsteinBij.push_back(x);
    }

    if (m_EinsteinAij.size() == m_EinsteinBij.size()) {
      m_bActivation = false;
    } else if (m_EinsteinAij.size() == 0 && m_EinsteinBij.size() > 0) {
      m_bActivation = true;
    } else if (m_EinsteinAij.size() > 0 && m_EinsteinBij.size() == 0) {
      m_bActivation = false;
    } else {
      throw(runtime_error("Error reading Einstein Coefficients."));
    }

    return true;
  }

  // Set transition frequencies.
  void gWellRadiationTransition::setTransitionFrequencies(vector<double>& frequencies) {
    m_TransitionFrequency = frequencies;

    // Derive Bij Coefficients from Aij coefficients if required.
    if (!m_bActivation && m_EinsteinBij.size() != m_EinsteinAij.size()) {
      m_EinsteinBij.clear();
      const double h = PlancksConstant_in_JouleSecond;
      const double A2B = 1.0e-06 / (8.0 * M_PI * h);
      for (size_t i(0); i < m_EinsteinAij.size(); i++) {
        const double nu = m_TransitionFrequency[i];
        const double Bij = A2B * m_EinsteinAij[i] / (pow(nu, 3.0));
        m_EinsteinBij.push_back(Bij);
      }
    }

    return;
  }

  // Calculate the attenuated temperature. This method calculates an effective radiation
  // temperature at a distant point from the surface of a black body (star). The basis of
  // the method is the attenuation of the Planck distribution for each frequency at a
  // radial distance from the surface followed by finding the temperature of the Planck 
  // distribution that has the closest match to the attenuated distribution. 
  double gWellRadiationTransition::AttenuatedTemperature(double radTemp, double attenuation) const {

    double attnFctr = 1.0 / (attenuation * attenuation);
    double radBeta  = 1.0 / (boltzmann_RCpK * radTemp);

    vector<double> attnInten(m_TransitionFrequency.size(), 0.0);

    // Calculate and store the attenated Planck distribution for each frequency.
    for (size_t idx(0); idx < m_TransitionFrequency.size(); idx++) {
      attnInten[idx] = attnFctr * planckDistribtuion(m_TransitionFrequency[idx], radBeta);
    }

    // Calculate the effective radiation temperature by finding the temperature
    // of the Planck distribution that matches the attenuated distribution.
    double upperT(radTemp);
    double beta = 1.0 / (boltzmann_RCpK * upperT);
    double sm1 = distributionDiff(attnInten, beta);

    double lowerT(1.0);
    beta = 1.0 / (boltzmann_RCpK * lowerT);
    double sm2 = distributionDiff(attnInten, beta);

    double localT(1.0);
    size_t count(0);
    do {
      localT = (upperT + lowerT) * 0.5;
      beta = 1.0 / (boltzmann_RCpK * localT);

      double sm3 = distributionDiff(attnInten, beta) ;

      if (sm1 > sm2) {
        upperT = localT;
        sm1 = sm3;
      }
      else {
        lowerT = localT;
        sm2 = sm3;
      }
      count++;
    } while (fabs(upperT - lowerT) > 0.1 && count < 100);

    if (count == 100)
      cinfo << "Warning: Maximum iterations exceeded in calculation of local radiation temperature." << endl;

    return localT;
  }

  // Method to calculate the difference between an attenuated distribution and 
  // the Planck distribution of a certain temperature.
  double gWellRadiationTransition::distributionDiff(vector<double> &attnInten, double beta) const {

    double sum(0.0);
    for (size_t idx(0); idx < 3; idx++) {
      sum += fabs(attnInten[idx] - planckDistribtuion(m_TransitionFrequency[idx], beta));
    }

    return sum;
  }

}//namespace
