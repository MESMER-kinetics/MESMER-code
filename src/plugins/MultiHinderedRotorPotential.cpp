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
#include <cmath>
#include "../formatfloat.h"
#include "../dMatrix.h"

using namespace std;
using namespace Constants;

namespace mesmer
{

  bool MultiHinderedRotorPotential::ReadPotentialPoints(PersistPtr pp) {

    PersistPtr pMHRP = pp;

    const char* p = pp->XmlReadValue("units", optional);
    m_units = p ? p : "kJ/mol";

    p = pp->XmlReadValue("bondRef");
    if (!p) {
      cerr << "There are no bonds defined." << endl;
      return false;
    }
    else {
      stringstream txt(p);
      string bondID;
      while (txt >> bondID) {
        m_bondIDs.push_back(bondID);
      }
    }

    m_expansion = pp->XmlReadInteger("expansionSize", optional);

    // Check if sine terms are to be used.
    // m_useSinTerms = pp->XmlReadBoolean("useSineTerms") || pp->XmlReadBoolean("UseSineTerms");

    double degToRad(M_PI / 180.0);
    while (pp = pp->XmlMoveTo("me:PotentialPoint"))
    {
      vector<double> tangles;
      const char* txt = pp->XmlReadValue("angles");
      istringstream anglesTxt(txt);
      double angle(0.0);
      while (anglesTxt >> angle) {
        tangles.push_back(angle * degToRad);
      }
      m_angles.push_back(tangles);

      // Check for a consistent number of angles.
      if (m_angles[0].size() != tangles.size()) {
        string errorMsg = "The number of angles is not the same for one row of a multi-hindered rotor potential definition.";
        throw(std::runtime_error(errorMsg));
      }

      double potentialPoint = pp->XmlReadDouble("potential", optional);
      if (IsNan(potentialPoint))
        potentialPoint = 0.0;
      potentialPoint = getConvertedEnergy(m_units, potentialPoint);
      potentialPoint = min(potentialPoint, 10000.0);
      m_potential.push_back(potentialPoint);
    }

    m_nVar = m_angles[0].size();
    if (m_nVar != m_bondIDs.size()) {
      string errorMsg = "The number of angles is different from the number of bond IDs in the definition of the multi hindered rotor potential.";
      throw(std::runtime_error(errorMsg));
    }

    // Shift data to minimum.
    size_t ii(0);
    double minEnergy(m_potential[0]);
    vector<double> tangles = m_angles[0];
    for (size_t i(0); i < m_potential.size(); i++) {
      if (m_potential[i] < minEnergy) {
        minEnergy = m_potential[i];
        tangles = m_angles[i];
        ii = i;
      }
    }

    for (size_t i(0); i < m_potential.size(); i++) {
      m_potential[i] -= minEnergy;
      for (size_t j(0); j < m_nVar; j++) {
        m_angles[i][j] -= tangles[j];
        if (m_angles[i][j] < 0.0)
          m_angles[i][j] += 2.0 * M_PI;
      }
    }

    pp = pMHRP->XmlMoveTo("me:TestLSqFit");
    if (pp) {
      m_testLSqFit = true;
    }

    return true;

  }

  void MultiHinderedRotorPotential::initializePotential() {

    size_t nbasis(0);

    nbasis = 1 + 2 * m_nVar * m_expansion + 2 * (m_nVar * (m_nVar - 1) * m_expansion * m_expansion / 2);
    // nbasis = 1 + 2 * m_nVar * m_expansion + 4 * (m_nVar * (m_nVar - 1) * m_expansion * m_expansion / 2);

    // Issue a warning if the number of basis functions exceeds the number of potential points.
    if (nbasis > m_potential.size()) {
      cwarn << "WARNING: The number of basis functions exceeds the number of potential " << endl
            << "points in the least square fit the coupled hinder rotor potential." << endl;
    }

    vector<double> wrk(nbasis, 0.0);
    dMatrix design(nbasis, 0.0);

    for (size_t itr(0); itr < m_potential.size(); itr++) {

      double ptnl = m_potential[itr];
      vector<double> angles = m_angles[itr];

      vector<double> basis(nbasis, 0.0);
      basis[0] = 1.0;
      wrk[0] += ptnl;

      for (size_t i(0), k(1); i < angles.size(); i++) {

        // alpha coefficients.

        double angle_i = angles[i];
        for (size_t n(1); n <= m_expansion; n++) {
          double nTheta_i = double(n) * angle_i;
          double cnT_i = cos(nTheta_i);
          double snT_i = sin(nTheta_i);
          basis[k] = cnT_i;
          wrk[k++] += ptnl * cnT_i;
          basis[k] = snT_i;
          wrk[k++] += ptnl * snT_i;
        }

        // beta coefficients.

        for (size_t n(1); n <= m_expansion; n++) {
          double nTheta_i = double(n) * angle_i;
          double cnT_i = cos(nTheta_i);
          double snT_i = sin(nTheta_i);
          for (size_t j(0); j < i; j++) {
            double angle_j = angles[j];
            for (size_t m(1); m <= m_expansion; m++) {
              double nTheta_j = double(m) * angle_j;
              double cnT_j = cos(nTheta_j);
              double snT_j = sin(nTheta_j);
              basis[k] = cnT_i * cnT_j;
              wrk[k++] += ptnl * cnT_i * cnT_j;
              basis[k] = snT_i * snT_j;
              wrk[k++] += ptnl * snT_i * snT_j;
              //basis[k] = cnT_i * snT_j;
              //wrk[k] += ptnl * basis[k];
            }
          }
        }

      }

      // Update design matrix.

      for (size_t i(0); i < basis.size(); i++) {
        design[i][i] += basis[i] * basis[i];
        for (size_t j(0); j < i; j++) {
          design[i][j] += basis[i] * basis[j];
          design[j][i] = design[i][j];
        }
      }

    }

    design.solveLinearEquationSet(&wrk[0]);
    m_Calpha = wrk;

    // Test the result.
    if (m_testLSqFit) {
      double residual(0.0);
      ctest << "\nCoupled hindered rotor potential (cm-1) for bonds :";
      for (size_t i(0); i < m_bondIDs.size(); i++)
        ctest << " " << m_bondIDs[i];
      ctest << endl << endl;
      for (size_t itr(0); itr < m_potential.size(); itr++) {

        vector<double> angles = m_angles[itr];

        double ptnl = calculatePotential(angles);
        double diff = m_potential[itr] - ptnl;
        residual += diff * diff;

        ctest << formatFloat(m_potential[itr], 6, 15) << formatFloat(ptnl, 6, 15) << formatFloat(diff, 6, 15) << endl;
      }

      residual = sqrt(residual/double(m_potential.size()));

      ctest << endl;
      ctest << "RMS Residual = " << formatFloat(residual, 6, 15) << endl << endl;
    }

  }

  double MultiHinderedRotorPotential::calculatePotential(std::vector<double>& angles) const {

    double ptnl(m_Calpha[0]);
    for (size_t i(0), k(1); i < angles.size(); i++) {

      // alpha coefficients.

      double angle_i = angles[i];
      for (size_t n(1); n <= m_expansion; n++) {
        double nTheta_i = double(n) * angle_i;
        double cnT_i = cos(nTheta_i);
        double snT_i = sin(nTheta_i);
        ptnl += m_Calpha[k++] * cnT_i;
        ptnl += m_Calpha[k++] * snT_i;
      }

      // beta coefficients.

      for (size_t n(1); n <= m_expansion; n++) {
        double nTheta_i = double(n) * angle_i;
        double cnT_i = cos(nTheta_i);
        double snT_i = sin(nTheta_i);
        for (size_t j(0); j < i; j++) {
          double angle_j = angles[j];
          for (size_t m(1); m <= m_expansion; m++) {
            double nTheta_j = double(m) * angle_j;
            double cnT_j = cos(nTheta_j);
            double snT_j = sin(nTheta_j);
            ptnl += m_Calpha[k++] * cnT_i * cnT_j;
            ptnl += m_Calpha[k++] * snT_i * snT_j;
            //basis[k] = cnT_i * snT_j;
          }
        }
      }

    }

    ptnl = max(ptnl, 0.0);

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

