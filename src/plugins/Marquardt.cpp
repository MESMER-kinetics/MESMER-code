//-------------------------------------------------------------------------------------------
//
// Marquardt.cpp
//
// Author: Struan Robertson
// Date:   19/Jun/2011
//
// This class implements the Levenberg-Marquardt non-linear least squares algorithm. 
//
//-------------------------------------------------------------------------------------------

#include <fstream>
#include <iomanip>

#include "../System.h"
#include "../ConditionsManager.h"
#include "../calcmethod.h"
#include "../dMatrix.h"
#include "FittingUtils.h"
#include "../TimeCounter.h"
#include "../Persistence.h"

namespace mesmer
{
	class Marquardt : public CalcMethod, private FittingUtils
	{
	public:

		Marquardt(const char* id) : FittingUtils(), m_id(id),
			m_nVar(0), m_delta(0.001), m_lambdaScale(10.0)
		{ Register(); }

		virtual ~Marquardt() {}
		virtual const char* getID()  { return m_id; }
    virtual Marquardt* Clone() { return new Marquardt(*this); }
    virtual bool ParseData(PersistPtr pp);

		//Function to do the work
		virtual bool DoCalculation(System* pSys);

	private:

		// Write out current variable values.
		void WriteVarVals(double chiSquare, double lambda) const ;

		const char* m_id;

		// Dimension of fit.
		size_t m_nVar ;

		// Numerical derivative delta.
		double m_delta ;

		// Factor for scalling lambda by during fitting.
		double m_lambdaScale ;

		unsigned m_maxIterations ;
		double m_tol ;
	};

	////////////////////////////////////////////////
	//Global instance
	Marquardt theMarquardt("marquardt");
	///////////////////////////////////////////////

	bool Marquardt::ParseData(PersistPtr pp)
	{
		// Read in Marquardt parameters, or use values from defaults.xml.
		m_delta = pp->XmlReadDouble("me:MarquardtDerivDelta");
		m_maxIterations= pp->XmlReadInteger("me:MarquardtIterations");
		m_tol = pp->XmlReadDouble("me:MarquardtTolerance");
		return true;
	}

	bool Marquardt::DoCalculation(System* pSys)
	{
		m_nVar = Rdouble::withRange().size() ;

		if (m_nVar < 1) { 
			cerr << "Marquardt requires at least one range variable to be set." << endl;
			return false ;
		}

		//Do not output all the intermediate results to XML
		pSys->m_Flags.overwriteXmlAnalysis = true;

		// Use the same grain numbers for for all calcuations regardless of 
		// temperature (i.e. reduce the number of times micro-rates are calculated).
		pSys->m_Flags.useTheSameCellNumber = true;

		// Uncomment to enable ctest output during fitting. Or use -w5 option in command.
		//ChangeErrorLevel e(obDebug); 

		//Default is to disable ctest during fitting. Restored when leaving this function.
		StopCTestOutput stop(true) ;

		//
		// Begin by finding the starting point chi-squared value.
		//

		vector<double> currentLocation(m_nVar,0.0) ; 
		vector<double> newLocation(m_nVar,0.0) ; 

		GetLocation(currentLocation) ;

		// Invoke SetLocation to catch any constrained parameters.
		SetLocation(currentLocation) ;

		double chiSquare(0.0), lambda(1.0) ;
		vector<double> residuals ;
		pSys->calculate(chiSquare, residuals) ;

		double bestChiSquare = chiSquare ;

		WriteVarVals(chiSquare, lambda) ;

		//
		// The following is slightly modified implementation of the Marquardt
		// algorithm. The modification is tha apllication of a bounds check on
		// each proposed new location of the minimum. If the the bounds check
		// fails lambda is increased with the consequence that the algorithm 
		// moves toward a short steepest decent algorithm. 
		//

		vector<double> gradient(m_nVar,0.0) ;
		qdMatrix hessian(m_nVar,0.0); 
		NumericalDerivatives(pSys, residuals, m_delta, gradient, hessian) ;

		bool converged(false) ;
		for (size_t itr(1) ; itr <= m_maxIterations && !converged ; itr++) {

			newLocation = currentLocation ;
			vector<qd_real> deltaLocation(m_nVar,0.0) ;
			qdMatrix invHessian = hessian ;

			for (size_t iVar(0) ; iVar < m_nVar ; iVar++) {
				invHessian[iVar][iVar] *= qd_real(1.0 + lambda) ;
				deltaLocation[iVar] = qd_real(gradient[iVar]) ;
			}

      try {
        invHessian.solveLinearEquationSet(&deltaLocation[0]);
      }
      catch (std::runtime_error& e)
      {
        string MtrxTitle("Marquardt Hessian Matrix: ");
        hessian.print(MtrxTitle, cerr);
        cinfo.flush();
        cwarn.flush();
        cerr.flush();
        throw(e);
      }

			for (size_t iVar(0) ; iVar < m_nVar ; iVar++) {
				newLocation[iVar] += to_double(deltaLocation[iVar]) ;
			}

			// Check bounds.    
			if (CheckBounds(newLocation)) {
				SetLocation(newLocation) ;

				pSys->calculate(chiSquare, residuals) ;

				if (chiSquare > bestChiSquare) {
					lambda *= m_lambdaScale ;
					SetLocation(currentLocation) ;
				} else {
					double relativeChange = 1.0 - chiSquare/bestChiSquare ;
					converged = (relativeChange < m_tol) ;
					lambda /= m_lambdaScale ;
					GetLocation(currentLocation) ;
					bestChiSquare = chiSquare ;
					NumericalDerivatives(pSys, residuals, m_delta, gradient, hessian) ;
				}
			} else {
				lambda *= m_lambdaScale ;
			}

      WriteVarVals(bestChiSquare, lambda) ;

			cinfo << "Iteration: " << itr << " of Marquardt. ChiSquare = " << bestChiSquare << ", Lambda = " << lambda << endl;

		}

		// Write meta data to the XML file.
		for (size_t i(0); i < m_nVar ; i++ ) {
			TimeCount events;
			std::string timeString;
			Rdouble::withRange()[i]->XmlWriteAttribute("fitted", events.setTimeStamp(timeString));
			stringstream cs;
			cs << chiSquare;
			Rdouble::withRange()[i]->XmlWriteAttribute("chiSquared", cs.str());
		}

		// Write out the results and the statisitics of the fit.
		ResultsAndStatistics(pSys, hessian) ;
		pSys->getConditionsManager()->get_generalAnalysisData()->setCovariance(hessian) ;

		return true;
	}

	//
	// Write out current variable values.
	//
	void Marquardt::WriteVarVals(double chiSquare, double lambda) const {

		cerr << endl << "Chi^2 = " << chiSquare << " Lambda = " << lambda << endl ;
		for(size_t iVar(0) ; iVar < m_nVar ; iVar++) {

			Rdouble var = *Rdouble::withRange()[iVar] ;
			cerr << var.get_varname() << "=" << setprecision(6) << double(var.originalUnits()) << "  "; 

		}
		cerr << endl ;

	}

} //namespace

