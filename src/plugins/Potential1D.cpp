//-------------------------------------------------------------------------------------------
//
// Potential1D.cpp
//
// Author: Struan Robertson
// Date:   25/Mar/2017
//
// Implementation of general 1D potential classes. 
//
//-------------------------------------------------------------------------------------------

#include "Potential1D.h"
#include "unitsConversion.h"
#include "Constants.h"
#include "MesmerConfig.h"
#include <algorithm>

using namespace std;
using namespace Constants;

namespace mesmer
{

  void AnalyticalPotential::InitializePotential(PersistPtr pp) {

    const char* p = pp->XmlReadValue("units", optional);
    m_units = p ? p : "kJ/mol";

		double val = pp->XmlReadDouble("minx", optional);
		if (!IsNan(val))
			m_minx = val;

		val = pp->XmlReadDouble("maxx", optional);
		if (!IsNan(val))
			m_maxx = val;

		vector<int> indicies;
    vector<double> coefficients;
    int maxIndex(0);
    while (pp = pp->XmlMoveTo("me:PotentialPoint"))
    {
      int index = pp->XmlReadInteger("index", optional);
      indicies.push_back(index);
      maxIndex = max(maxIndex, index);

      double coefficient = pp->XmlReadDouble("coefficient", optional);
      if (IsNan(coefficient))
        coefficient = 0.0;
      coefficients.push_back(coefficient);
    }

    // As coefficients can be supplied in any order, they are sorted here.
    m_Coeff.resize(++maxIndex, 0.0);
    for (size_t i(0); i < coefficients.size(); i++) {
      m_Coeff[indicies[i]] = coefficients[i];
    }

  }

  void AnalyticalPotential::calculatePotential(vector<double> &abscissa, vector<double> &potential) const {

    abscissa.resize(m_npoints, 0.0);
    potential.resize(m_npoints, 0.0);

    // Calculate abscissa points
    double deltax = (m_maxx - m_minx) / double(m_npoints-1);
    for (size_t i(0); i < abscissa.size(); i++) {
      abscissa[i] = m_minx + double(i)*deltax;
      double pot(m_Coeff[0]), tmp(1.0);
      for (size_t j(1); j < m_Coeff.size(); j++) {
        tmp *= abscissa[i];
        pot += m_Coeff[j] * tmp;
      }
      potential[i] = pot;
    }

    // Convert to wavenumbers.
    for (size_t i(0); i < potential.size(); i++) {
      potential[i] = getConvertedEnergy(m_units, potential[i]);
    }
  }

  void NumericalPotential::InitializePotential(PersistPtr pp) {

    const char* p = pp->XmlReadValue("units", optional);
    m_units = p ? p : "kJ/mol";

    vector<double> potential;
    vector<double> angle;

    while (pp = pp->XmlMoveTo("me:PotentialPoint"))
    {
      double anglePoint = pp->XmlReadDouble("angle", optional);
      if (IsNan(anglePoint))
        anglePoint = 0.0;
      angle.push_back(anglePoint);

      double potentialPoint = pp->XmlReadDouble("potential", optional);
      if (IsNan(potentialPoint))
        potentialPoint = 0.0;
      potential.push_back(potentialPoint);
    }

		m_spline.Initialize(angle, potential);

  }

  void NumericalPotential::calculatePotential(vector<double> &abscissa, vector<double> &potential) const {

    // Finally convert to wavenumbers.
    for (size_t i(0); i < potential.size(); i++) {
      potential[i] = getConvertedEnergy(m_units, potential[i]);
    }
  }

}//namespace

