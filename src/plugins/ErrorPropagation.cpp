//-------------------------------------------------------------------------------------------
//
// ErrorPropagation.cpp
//
// Author: Struan Robertson
// Date:   09/Jan/2016
//
// This class implements an error analysis algorithm. This is a simple Monte Carlo analysis 
// based on the assumption of guassian distributions in the model parameters. It will return,
// for a given set of temperature and pressure condiditons, estimates of the errors in the
// rate coefficient. It requires a covariance matrix of the parameters, as obtained from the
// Levenberg-Marquardt algorithm - hence this method will generally be applied in conjunction
// with a Levenberg-Marquardt fit. The error estimates generated can be used in macroscopic
// sensitivity analysis.
//
//-------------------------------------------------------------------------------------------

#include <fstream>
#include <iomanip>

#include "../System.h"
#include "../calcmethod.h"
#include "FittingUtils.h"
#include "../Sobol.h"
#include <string>
#include <cctype>

namespace {

	void ToUpper(std::string &str) {
		for (size_t i(0); i < str.size(); i++) {
			char c = std::toupper(str[i]);
			str[i] = c;
		}
	}

}

namespace mesmer
{
	class ErrorPropagation : public CalcMethod, private FittingUtils
	{
	public:

		ErrorPropagation(const char* id) :
			m_id(id),
			m_pSA(),
			m_CorrelMtx(NULL),
			m_nVar(0),
			m_nOut(0),
			m_nSample(0) {
			Register();
		}

		virtual ~ErrorPropagation() { }
		virtual const char* getID() { return m_id; }
		virtual bool ParseData(PersistPtr pp);

		// Function to do the work.
		virtual bool DoCalculation(System* pSys);

	private:

		typedef vector<vector<double> > Table;
		typedef map<string, double>::const_iterator RxnItr;

		// This method writes out the results of a sensitivity analysis to test file.
		bool WriteOutAnalysisToTest(const vector<string> &rxnId, const vector<double> &f0, const vector<double> &varf, double Temperature, double Concentration);

		// This method writes out the results of a sensitivity analysis.
		bool WriteOutAnalysis(const vector<string> &rxnId, const vector<double> &f0, const vector<double> &varf, double Temperature, double Concentration);

		// This method writes the input variable key.
		void WriteInputVariableKey(stringstream &Key) const;

		void rndLocation(const vector<double> &rndmd, const vector<double> &currentLoc, vector<double> &newLoc) {
			size_t nVar = rndmd.size();
			// Take inverse cumulative distribution of each sobol element.
			for (size_t j(0); j < nVar; j++) {
				newLoc[j] = NormalCDFInverse(rndmd[j]);
			}
			// Multiply the InvNorm with the cholesky decompostion.
			newLoc *= (*m_CorrelMtx);
			for (size_t j(0); j < nVar; j++) {
				newLoc[j] += currentLoc[j];
			}
		}

		const char* m_id;

		PersistPtr m_pSA;

		dMatrix *m_CorrelMtx;

		size_t m_nVar;             // Dimension of analysis - number of inputs.
		size_t m_nOut;             // Dimension of analysis - number of outputs.

		size_t m_nSample;

	};

	////////////////////////////////////////////////
	//Global instance
	ErrorPropagation theErrorPropagation("ErrorPropagation");
	///////////////////////////////////////////////

	bool ErrorPropagation::ParseData(PersistPtr pp)
	{
		// Read in sensitivity analysis parameters, or use values from defaults.xml.
		m_nSample = pp->XmlReadInteger("me:errorPropagationSamples");

		// Store pointer for output.
		m_pSA = pp;

		return true;
	}

	bool ErrorPropagation::DoCalculation(System* pSys) {

		m_nVar = Rdouble::withRange().size();

		if (m_nVar < 1) {
			cerr << "Error propagation requries at least one range variable to be set." << endl;
			return false;
		}

		// Read in correlation matrix.

		PersistPtr pca = pSys->getPersistPtr()->XmlMoveTo("me:analysis") ;
		PersistPtr pcm = pca->XmlMoveTo("me:covariance");
    if(!pcm)
      pcm = pca->XmlMoveTo("me:hessian");
		pcm = pcm->XmlMoveTo("matrix");
    if ((m_CorrelMtx = ReadMatrix<double>(pcm))) {
			m_CorrelMtx->cholesky();
		}
		else {
			cerr << "A correlation matrix (me:covariance or me:hessian) is required for error propagation but was not found" << endl;
			return false;
		}

		// Do not output all the intermediate results to XML
		pSys->m_Flags.overwriteXmlAnalysis = true;

		// Use the same grain numbers for for all calcuations regardless of 
		// temperature (i.e. reduce the number of times micro-rates are calculated).
		pSys->m_Flags.useTheSameCellNumber = true;

		vector<double> originalLocation(m_nVar, 0.0);
		vector<double> newLocation(m_nVar, 0.0);

		GetLocation(originalLocation);

		// Invoke SetLocation to catch any constrained parameters.
		SetLocation(originalLocation);

		vector<double> Temperature;
		vector<double> Concentration;
		pSys->getConditionsManager()->getConditions(Temperature, Concentration);

		// Instantiate a random vector generator.
		Sobol sobol;

		// Loop over condiditons. 
		size_t nConditions = Temperature.size();
		for (size_t nCnd(0); nCnd < nConditions; nCnd++) {

			//Default is to disable ctest during analysis. Restored when leaving loop.
			StopCTestOutput stop(true);

			// Test calculation to determine number of outputs.
			SetLocation(originalLocation);
			map<string, double> phenRates;
			pSys->calculate(nCnd, phenRates);

			m_nOut = phenRates.size();
			vector<double> f0(m_nOut, 0.0), varf(m_nOut, 0.0);
			vector<string> rxnId;
			RxnItr irxn = phenRates.begin();
			for (size_t nOut(0); irxn != phenRates.end(); irxn++, nOut++) {
				f0[nOut] = irxn->second ;
				rxnId.push_back(irxn->first);
			}

			// Loop over perturbed parameter values.

			long long seed(7);
			for (size_t itr(1), cnt(0) ; itr <= m_nSample ; itr++) {
				vector<double> rndmd(m_nVar, 0.0);
				sobol.sobol(rndmd.size(), &seed, rndmd);

				// Use random vector generated by sobol method to perturb parameter values.

				rndLocation(rndmd, originalLocation, newLocation);

				// Set perturbed parameters and calculate new quantities.

				SetLocation(newLocation);

				map<string, double> phenRates;

				// As some perturbations will produce unrealistic configurations the
				// diagonalizer will fail, in which case this configuration is skipped.
				try {
					pSys->calculate(nCnd, phenRates);
				}
				catch (...) {
					continue;
				}

				cnt++ ;
				irxn = phenRates.begin();
				for (size_t nOut(0); irxn != phenRates.end(); irxn++, nOut++) {
					double output = irxn->second;
					double tmp = f0[nOut] - output;
					varf[nOut] = (double(cnt-1)*varf[nOut] + tmp*tmp) / double(cnt);
				}

			}

			ctest.clear();

			// Write out results. 
			WriteOutAnalysisToTest(rxnId, f0, varf, Temperature[nCnd], Concentration[nCnd]);

			WriteOutAnalysis(rxnId, f0, varf, Temperature[nCnd], Concentration[nCnd]);

		} // End of conditions loop.

		return true;
	}

	// This method writes out the results of a sensitivity analysis.
	bool ErrorPropagation::WriteOutAnalysis(const vector<string> &rxnId, const vector<double> &f0, const vector<double> &varf, double Temperature, double Concentration) {

		// Begin table.

		PersistPtr pp = m_pSA->XmlWriteElement("me:errorPropagationTable");

		// Write out conditions

		pp->XmlWriteAttribute("temperature", Temperature, 2, true);
		pp->XmlWriteAttribute("concentration", Concentration, 2, false);

		for (size_t i(0), idx(0); i < m_nOut; i++) {

			PersistPtr ppSensInd = pp->XmlWriteElement("me:propagatedErrors");
			ppSensInd->XmlWriteAttribute("reaction", rxnId[i]);


			// Write out R^2 statistic and Standard deviation.

			stringstream ss;
			ss.str("");
			ss << sqrt(varf[i]);
			ppSensInd->XmlWriteValueElement("me:standardDeviation", ss.str());
      ss.str("");
      ss << f0[i];
			ppSensInd->XmlWriteValueElement("me:rateCoefficient", ss.str());
		}

		return true;
	}

	// Write values to .test for inspection.
	bool ErrorPropagation::WriteOutAnalysisToTest(const vector<string> &rxnId, const vector<double> &f0, const vector<double> &varf, double Temperature, double Concentration) {

		ctest << endl << "Error Propagation: Temperature = " << setw(5) << setprecision(4) << Temperature << ", Concentration = " << Concentration << endl << endl;

		stringstream Key;
		WriteInputVariableKey(Key);
		ctest << Key.str() << endl;

		for (size_t i(0), idx(0); i < m_nOut; i++) {
			ctest << " Reaction           : " << rxnId[i] << endl ;
			ctest << " Rate Coeffcient    : " << f0[i] << endl;
			ctest << " Standard Deviation : " << sqrt(varf[i]) << endl << endl;
		}

		return true;
	}

	// This method writes the input variable key.
	void ErrorPropagation::WriteInputVariableKey(stringstream &Key) const {
		Key << "Input variables:" << endl << endl;
		for (size_t iVar(0), idx(1); iVar < m_nVar; iVar++, idx++) {
			Rdouble var = *Rdouble::withRange()[iVar];
			Key << "  " << setw(30) << left << var.get_varname() << ": (" << idx << ")" << endl;
		}
		Key << endl;
	}

} //namespace
