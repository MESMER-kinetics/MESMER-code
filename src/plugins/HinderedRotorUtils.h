#ifndef GUARD_HinderedRotorUtils_h
#define GUARD_HinderedRotorUtils_h

//-------------------------------------------------------------------------------------------
//
// HinderedRotorUtils.h
//
// Author: Struan Robertson
// Date:   1/Mar/2014
//
// Definition of a utility class that is inherited by the hindered rotor methods. 
//
//-------------------------------------------------------------------------------------------

#include "../MolecularComponents.h"
#include <vector>

namespace mesmer
{
  class gStructure ;

  class HinderedRotorUtils: public DensityOfStatesCalculator
  {
  public:

	HinderedRotorUtils(const char* id) : DensityOfStatesCalculator(),
      m_id(id),
      m_bondID(),
	  m_potentialCosCoeff(),
	  m_potentialSinCoeff(),
	  m_expansion(4),
	  m_useSinTerms(false),
	  m_calIntrlIrt(false) { Register(); } ; 
	~HinderedRotorUtils() {} ; 

    virtual const char* getID()  { return m_id; }
	void set_BondID(const char *BondID) {m_bondID = string(BondID) ; } ;
	string get_BondID () const { return m_bondID ; } ;

  protected:

	std::vector<double>& get_PotentialCosCoeff() {return m_potentialCosCoeff ; } ;
	std::vector<double>& get_PotentialSinCoeff() {return m_potentialSinCoeff ; } ;
	void set_Expansion(size_t expansion) {m_expansion = expansion ; } ;
	void set_UseSinTerms(bool useSinTerms) {m_useSinTerms = useSinTerms ; } ;
	void set_CalIntrlIrt(bool calIntrlIrt) {m_calIntrlIrt = calIntrlIrt ; } ;
	size_t get_Expansion() const {return m_expansion ; } ;
	bool get_UseSinTerms() const {return m_useSinTerms ; } ;
	bool get_CalIntrlIrt() const {return m_calIntrlIrt ; } ;

	// Calculate the Fourier coefficients from potential data points.
	void FourierCoeffs(vector<double> &angle, vector<double> &potential) ;

	// Calculate potential.
	double CalculatePotential(double angle) const ;

	// Calculate gradient.
	double CalculateGradient(double angle) const ;

  private:

    const char* m_id;

    std::string m_bondID;

	std::vector<double> m_potentialCosCoeff ; // The cosine coefficients of the hindered rotor potential.
	std::vector<double> m_potentialSinCoeff ; // The sine coefficients of the hindered rotor potential.

	size_t m_expansion ;                 // Number of coefficients in the cosine expansion.

	bool m_useSinTerms ;                 // If true sine terms are used in the representation of the potential.

	bool m_calIntrlIrt ;                 // If true the internal rotor moment of inertia will be calculated as
	                                     // function the angle of rotation about the associated bond.
  } ;

}  //namespace

#endif // GUARD_HinderedRotorUtils_h
