#ifndef GUARD_MultiHinderedRotorPotential_h
#define GUARD_MultiHinderedRotorPotential_h

//-------------------------------------------------------------------------------------------
//
// MultiHinderedRotorPotential.h
//
// Author: Struan Robertson
// Date:   18/Dec/2023
//
// Definition of a the MultiHinderedRotorPotential class. The principal use of this class 
// is to support multiple hinder rotor calculations by providing a means of using multiple
// conformer potential points to produce an inierpolated potential.
//
//-------------------------------------------------------------------------------------------

#include "../XMLPersist.h"
#include <vector>
// #include "../Spline.h"

namespace mesmer
{

  // Abstract base class.

  class MultiHinderedRotorPotential
  {
  public:

		MultiHinderedRotorPotential() : m_units("kJ/mol"), m_expansion(1), m_potential(), 
      m_angles(), m_nVar(0), m_Calpha(), m_bondIDs(), m_testLSqFit(false)
    { };
		virtual ~MultiHinderedRotorPotential() {};

    bool ReadPotentialPoints(PersistPtr pp);

    void initializePotential();

    double calculatePotential(std::vector<double> &abscissa) const;

    void getBondIDs(std::vector<std::string> &bondIDSs) const {
      bondIDSs = m_bondIDs;
    };

  protected:

    std::string m_units;              // Units that potential energy is expressed in.
    size_t m_expansion;               // Number of coefficients in the cosine expansion.
    std::vector<double> m_potential;  // Potential values.
    std::vector<std::vector<double> > m_angles; // Conformer angles.
    size_t m_nVar;
    std::vector<double> m_Calpha;
    std::vector<std::string> m_bondIDs;
    bool m_testLSqFit;
  };

}  //namespace

#endif // GUARD_MultiHinderedRotorPotential_h
