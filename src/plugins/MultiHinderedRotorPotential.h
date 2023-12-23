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
#include "../Spline.h"

namespace mesmer
{

  // Abstract base class.

  class MultiHinderedRotorPotential
  {
  public:

		MultiHinderedRotorPotential() : m_units("kJ/mol"), m_expansion(1),
      m_useSinTerms(false), m_potential(), m_angles(), m_nVar(0), m_f0(0.0), m_Calpha(), m_Cbeta(), m_Salpha(), m_Sbeta()
    { };
		virtual ~MultiHinderedRotorPotential() {};

    bool ReadPotentialPoints(PersistPtr pp);

    void MultiHinderedRotorPotential::initializePotential();

    double calculatePotential(std::vector<double> &abscissa) const;

  protected:

    std::string m_units;              // Units that potential energy is expressed in.
    size_t m_expansion;               // Number of coefficients in the cosine expansion.
    bool m_useSinTerms;               // If true sine terms are used in the representation of the potential.
    std::vector<double> m_potential;  // Potential values.
    std::vector<std::vector<double> > m_angles; // Conformer angles.
    size_t m_nVar;
    double m_f0 ;
    std::vector<double> m_Calpha;
    std::vector<double> m_Cbeta;
    std::vector<double> m_Salpha;
    std::vector<double> m_Sbeta;
  };


}  //namespace

#endif // GUARD_MultiHinderedRotorPotential_h
