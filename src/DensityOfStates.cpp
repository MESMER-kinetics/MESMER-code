//-------------------------------------------------------------------------------------------
//
// DensityOfStates.cpp
//
// Author: Struan Robertson
// Date:   17/Nov/2022
//
// This file contains the implementation of various base class methods.
//
//-------------------------------------------------------------------------------------------
#include "DensityOfStates.h"
#include "Constants.h"
#include "formatfloat.h"

using namespace Constants;

namespace mesmer
{
  // Method to create thermodymaic table header.
  void DensityOfStatesCalculator::ThermoContribHeader(const std::string species, const std::string id, const std::string units) {
    m_thermoContribTable = true;
    std::stringstream ss;
    ss << std::endl;
    ss << species << " : " << id << std::endl;
    ss << std::endl;
    ss << "      T/K   H(T) - H(0)/         S(T)/         G(T)/         C(T)/" << std::endl;
    if (units == "kJ/mol") {
      ss << "                  kJ/mol       J/mol/K        kJ/mol       J/mol/K" << std::endl;
    }
    else {
      ss << "                kcal/mol     cal/mol/K      kcal/mol     cal/mol/K" << std::endl;
    }
    m_thermoDynamicTable = ss.str();
  }

  // Method to add an entry to thermodymaic table header.
  void DensityOfStatesCalculator::ThermoDynamicEntry(const double beta, const double CanPrtnFn, const double internalEnergy, const double varEnergy) {
    if (m_thermoContribTable) {
      double temp = 1.0 / (beta * boltzmann_RCpK);
      double FreeEnergy = -m_unitFctr * log(CanPrtnFn) / beta;
      double enthalpy = m_unitFctr * internalEnergy;
      double entropy = (enthalpy - FreeEnergy) * beta * boltzmann_RCpK * 1000.0;
      double heatCapacity = m_unitFctr * boltzmann_RCpK * (beta * beta * varEnergy) * 1000.0;
      std::stringstream ss;
      ss << setw(10) << fixed << setprecision(2) << temp << formatFloat(enthalpy, 5, 14) << formatFloat(entropy, 5, 14)
        << formatFloat(FreeEnergy, 5, 14) << formatFloat(heatCapacity, 5, 14) << std::endl;
      m_thermoDynamicTable += ss.str();
    }
  }

  // Write out table of contribution to thermodynamic functions.
  string DensityOfStatesCalculator::ThermoDynamicWrite() const {
    return m_thermoDynamicTable;
  }


}//namespace
