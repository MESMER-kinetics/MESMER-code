//-------------------------------------------------------------------------------------------
//
// MultiHinderedRotorPotential.cpp
//
// Author: Struan Robertson
// Date:   18/Dec/2023
//
// Implementation of the  MultiHinderedRotorPotential class. 
//
//-------------------------------------------------------------------------------------------

#include "MultiHinderedRotorPotential.h"
#include "../unitsConversion.h"
#include "../Constants.h"
#include "../MesmerConfig.h"
#include <algorithm>
#include <stdexcept>
#include "../formatfloat.h"
#include "../dMatrix.h"

using namespace std;
using namespace Constants;

namespace mesmer
{

  bool MultiHinderedRotorPotential::ReadPotentialPoints(PersistPtr pp) {

    const char* p = pp->XmlReadValue("units", optional);
    m_units = p ? p : "kJ/mol";

    m_expansion = pp->XmlReadInteger("expansionSize", optional);

    // Check if sine terms are to be used.
    m_useSinTerms = pp->XmlReadBoolean("useSineTerms") || pp->XmlReadBoolean("UseSineTerms");

    while (pp = pp->XmlMoveTo("me:PotentialPoint"))
    {
      vector<double> tangles;
      const char* txt = pp->XmlReadValue("angles");
      istringstream anglesTxt(txt);
      double angle(0.0);
      while (anglesTxt >> angle) {
        tangles.push_back(angle);
      }
      m_angles.push_back(tangles);

      // Check for a consistent number of angles.
      if (m_angles[0].size() != tangles.size()) {
        cerr << "The number of angles is not the same for one row of a multi-hindered rotor potential definition." << endl;
        return false;
      }

      double potentialPoint = pp->XmlReadDouble("potential", optional);
      if (IsNan(potentialPoint))
        potentialPoint = 0.0;
      potentialPoint = getConvertedEnergy(m_units, potentialPoint);
      m_potential.push_back(potentialPoint);
    }

    m_nVar = m_angles[0].size();

    return true;

  }

  void MultiHinderedRotorPotential::initializePotential() {

    size_t nbasis(0);

    nbasis = 2 * m_nVar * m_expansion + 1;

    vector<double> wrk(nbasis, 0.0);
    dMatrix design(nbasis, 0.0);

    double degToRad(M_PI / 180.0);

    for (size_t itr(0); itr < m_potential.size(); itr++) {

      double ptnl = m_potential[itr];
      vector<double> angles = m_angles[itr];

      vector<double> basis(nbasis, 0.0);
      basis[0] = 1.0;
      wrk[0] += ptnl;

      for (size_t i(0), k(1), l(m_expansion + 1); i < angles.size(); i++, k += m_expansion, l += m_expansion) {

        // alpha coefficients.

        double angle_i = degToRad * angles[i];
        for (size_t j(1); j <= m_expansion; j++, k++, l++) {
          double nTheta_i = double(j) * angle_i;
          basis[k] = cos(nTheta_i);
          basis[l] = sin(nTheta_i);
          wrk[k] += ptnl * basis[k];
          wrk[l] += ptnl * basis[l];
        }
      }

      // Update design matrix.

      for (size_t m(0); m < basis.size(); m++) {
        design[m][m] += basis[m] * basis[m];
        for (size_t n(0); n < m; n++) {
          design[m][n] += basis[m] * basis[n];
          design[n][m]  = design[m][n];
        }
      }

    }

    design.solveLinearEquationSet(&wrk[0]);
    m_Calpha = wrk;

    // Test the result.

    double residule(0.0);
    for (size_t itr(0); itr < m_potential.size(); itr++) {

      vector<double> angles = m_angles[itr];

      double ptnl = calculatePotential(angles);
      double diff = m_potential[itr] - ptnl;
      residule += diff * diff;

      cout << formatFloat(m_potential[itr], 6, 15) << formatFloat(ptnl, 6, 15) << formatFloat(diff, 6, 15) << endl;

    }

    cout << endl;
    cout << "Residule = " << formatFloat(residule, 6, 15) << endl;

  }

  double MultiHinderedRotorPotential::calculatePotential(std::vector<double>& angles) const {

    double degToRad(M_PI / 180.0);
    double ptnl(m_Calpha[0]);
    for (size_t i(0), k(1), l(m_expansion + 1); i < angles.size(); i++, k += m_expansion, l += m_expansion) {
      double angle_i = degToRad * angles[i];
      for (size_t j(1); j <= m_expansion; j++, k++, l++) {
        double nTheta_i = double(j) * angle_i;
        ptnl += m_Calpha[k]*cos(nTheta_i);
        ptnl += m_Calpha[l]*sin(nTheta_i);
      }

    }

    return ptnl;

  }

  //void MultiHinderedRotorPotential::initializePotential_FT() {

  //  m_f0 = 0.0;
  //  m_Calpha.resize(m_nVar * m_expansion, 0.0);
  //  m_Cbeta.resize(m_nVar * (m_nVar - 1) * m_expansion * m_expansion / 2, 0.0);
  //  m_Salpha.resize(m_nVar * m_expansion, 0.0);
  //  m_Sbeta.resize(m_nVar * (m_nVar - 1) * m_expansion * m_expansion / 2, 0.0);

  //  double degToRad(M_PI / 180.0);

  //  for (size_t itr(0); itr < m_potential.size(); itr++) {

  //    double ptnl = m_potential[itr];
  //    vector<double> angles = m_angles[itr];

  //    m_f0 += ptnl;

  //    for (size_t i(0), ida(0), idb(0); i < angles.size(); i++) {

  //      // alpha coefficients.

  //      double angle_i = degToRad * angles[i];
  //      for (size_t k(1); k <= m_expansion; k++, ida++) {
  //        double nTheta_i = double(k) * angle_i;
  //        m_Calpha[ida] += ptnl * cos(nTheta_i);
  //        m_Salpha[ida] += ptnl * sin(nTheta_i);
  //      }

  //      // beta coefficients.

  //      for (size_t j(0); j < i; j++) {
  //        double angle_j = degToRad * angles[j];
  //        for (size_t k(1); k <= m_expansion; k++) {
  //          double nTheta_i = double(k) * angle_i;
  //          for (size_t l(1); l <= m_expansion; l++, idb++) {
  //            double nTheta_j = double(l) * angle_j;
  //            m_Cbeta[idb] += ptnl * cos(nTheta_i) * cos(nTheta_j);
  //            // m_Sbeta[idb] += ptnl * sin(nTheta_i) * sin(nTheta_j);
  //          }
  //        }
  //      }

  //    }
  //  }

  //  // Normalize expansion coefficients.

  //  double nrmlFctr(1.0 / double(M_PI * m_potential.size()));

  //  m_f0 *= 0.5 * nrmlFctr;

  //  for (size_t j(0); j < m_Calpha.size(); j++) {
  //    m_Calpha[j] *= nrmlFctr;
  //    m_Salpha[j] *= nrmlFctr;
  //  }

  //  for (size_t j(0); j < m_Cbeta.size(); j++) {
  //    m_Cbeta[j] *= nrmlFctr;
  //    m_Sbeta[j] *= nrmlFctr;
  //  }

  //  // Test the result.

  //  double residule(0.0);
  //  for (size_t itr(0); itr < m_potential.size(); itr++) {

  //    vector<double> angles = m_angles[itr];

  //    double ptnl = calculatePotential(angles);
  //    double diff = m_potential[itr] - ptnl;
  //    residule += diff * diff;

  //    cout << formatFloat(m_potential[itr], 6, 15) << formatFloat(ptnl, 6, 15) << formatFloat(diff, 6, 15) << endl;

  //  }

  //  cout << endl;
  //  cout << "Residule = " << formatFloat(residule, 6, 15) << endl;

  //}

  //double MultiHinderedRotorPotential::calculatePotential_FT(std::vector<double>& angles) const {

  //  double degToRad(M_PI / 180.0);

  //  // Initialize sum with the mean of the function.
  //  double sum(m_f0);

  //  for (size_t i(0), ida(0), idb(0); i < m_nVar; i++) {

  //    // alpha coefficients.

  //    double angle_i = degToRad * angles[i];
  //    for (size_t k(1); k <= m_expansion; k++, ida++) {
  //      double nTheta_i = double(k) * angle_i;
  //      sum += m_Calpha[ida] * cos(nTheta_i);
  //      sum += m_Salpha[ida] * sin(nTheta_i);
  //    }

  //    // beta coefficients.

  //    for (size_t j(0); j < i; j++) {
  //      double angle_j = degToRad * angles[j];
  //      for (size_t k(1); k <= m_expansion; k++) {
  //        double nTheta_i = double(k) * angle_i;
  //        for (size_t l(1); l <= m_expansion; l++, idb++) {
  //          double nTheta_j = double(l) * angle_j;
  //          sum += m_Cbeta[idb] * cos(nTheta_i) * cos(nTheta_j);
  //          sum += m_Sbeta[idb] * sin(nTheta_i) * sin(nTheta_j);
  //        }
  //      }
  //    }

  //  }

  //  return sum;

  //}

}//namespace

