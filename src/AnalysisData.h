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
		AnalysisData() : m_number(0), m_selection(""), m_eigenvalues() {}

		// Eigenvalues
		int m_number;
		std::string m_selection;
		std::vector<double> m_eigenvalues;
	};


}//namespace


#endif // GUARD_AnalysisData_h

