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
#include "../DensityOfStates.h"
#include <vector>

namespace mesmer
{
	class gStructure;

	class HinderedRotorUtils : public DensityOfStatesCalculator
	{
	public:

		HinderedRotorUtils(const char* id) : DensityOfStatesCalculator(),
			m_id(id),
			m_bondID(),
			m_expansion(1),
			m_useSinTerms(false),
			m_calIntrlIrt(false),
			m_phase(0.0)
		{
			Register();
		};
		~HinderedRotorUtils() {};

		virtual const char* getID() { return m_id; }
		void set_BondID(const char *BondID) { m_bondID = string(BondID); };
		string get_BondID() const { return m_bondID; };
		virtual string getName() { return string(m_id) + ": Bond " + m_bondID; };

	protected:

		void set_Expansion(size_t expansion) { m_expansion = expansion; };
		void set_UseSinTerms(bool useSinTerms) { m_useSinTerms = useSinTerms; };
		void set_CalIntrlIrt(bool calIntrlIrt) { m_calIntrlIrt = calIntrlIrt; };
		void set_Phase(double phase) { m_phase = phase; };
		size_t get_Expansion() const { return m_expansion; };
		bool get_UseSinTerms() const { return m_useSinTerms; };
		bool get_CalIntrlIrt() const { return m_calIntrlIrt; };
		double get_Phase() const { return m_phase; };

		// Read potential parameters
		void ReadPotentialParameters(PersistPtr ppHR, const string& SpeciesID, const string& bondID, vector<double>& ptnlCosCoeff, vector<double>& ptnlSinCoeff);

		// Calculate the Fourier coefficients from potential data points.
		void PotentialFourierCoeffs(vector<double> &angle, vector<double> &potential, vector<double> &ptnlCosCoeff, vector<double> &ptnlSinCoeff);

		// Calculate potential.
		double CalculatePotential(double angle, const vector<double> &ptnlCosCoeff, const vector<double> &ptnlSinCoeff) const;

		// Calculate gradient.
		double CalculateGradient(double angle, const vector<double> &ptnlCosCoeff, const vector<double> &ptnlSinCoeff) const;

	private:

		// Shift potential to origin.
		void ShiftPotential(vector<double> &ptnlCosCoeff, vector<double> &ptnlSinCoeff);

		const char* m_id;

		std::string m_bondID;

		size_t m_expansion;                      // Number of coefficients in the cosine expansion.

		bool m_useSinTerms;                      // If true sine terms are used in the representation of the potential.

		bool m_calIntrlIrt;                      // If true the internal rotor moment of inertia will be calculated as
												 // function the angle of rotation about the associated bond.

		double m_phase;                          // Phase difference, in degrees, between potential and configuration.
	};

}  //namespace

#endif // GUARD_HinderedRotorUtils_h
