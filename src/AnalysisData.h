#ifndef GUARD_AnalysisData_h
#define GUARD_AnalysisData_h

#include <string>
#include <sstream>
#include <vector>

namespace mesmer
{

	struct AnalysisData
	{
	public:
		AnalysisData() : m_number(0), m_selection(""), m_eigenvalues(), m_lossRef(), m_lossRateCoeff() {}

		~AnalysisData() {};

		void clear() {
			m_eigenvalues.clear();
			m_lossRef.clear();
			m_lossRateCoeff.clear();
			m_firstOrderReactionType.clear();
			m_firstOrderFromRef.clear();
			m_firstOrderToRef.clear();
			m_firstOrderRateCoeff.clear();
			m_timePoints.clear();
			m_aveEnergyRef.clear();
			m_aveEnergy.clear();
			m_PopRef.clear();
			m_Pop.clear();
		};

		// Eigenvalues.
		int m_number;
		std::string m_selection;
		std::vector<double> m_eigenvalues;

		// Loss rate coefficients.
		std::vector<std::string> m_lossRef;
		std::vector<double> m_lossRateCoeff;

		// First order rate coefficients.
		std::vector<std::string> m_firstOrderReactionType;
		std::vector<std::string> m_firstOrderFromRef;
		std::vector<std::string> m_firstOrderToRef;
		std::vector<double> m_firstOrderRateCoeff;

		// Average Energy and population.
		std::vector<double> m_timePoints;
		std::vector<std::string> m_aveEnergyRef;
		std::vector<double> m_aveEnergy;
		std::vector<std::string> m_PopRef;
		std::vector<double> m_Pop;

	};


}//namespace


#endif // GUARD_AnalysisData_h

