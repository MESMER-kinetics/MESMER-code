//-------------------------------------------------------------------------------------------
//
// Author: Dave Glowacki (with significant refactoring by Struan Robertson)
// Date:   2 June 2015
//
// Produces Zhu-Nakamura spin forbidden crossing coefficients
// Note that this implementation reproduces the results in Fig 5 & Tables II, III, IV in 
//    Zhu & Nakamura, JCP, 101(12), 10630, 1994
//
//-------------------------------------------------------------------------------------------
#include <vector>
#include <math.h>
#include <string>
#include "../System.h"
#include "AdiabaticCurve.h"
#include "../Rdouble.h"
#include "../gDensityOfStates.h"
#include "../Spline.h"

using namespace Constants;
using namespace std;

namespace mesmer
{

  class ZhuNakamuraCrossing : public MicroRateCalculator
  {
  public:
    virtual bool ParseData(PersistPtr pp);
    // Constructor which registers with the list of CrossingCalculators in the base class
    ZhuNakamuraCrossing(const char* id) : m_id(id), dataIsInTS(false), 
      m_H12(0.0), 
      m_ReducedMass(0.0), 
      m_harmonicX0(0.0), 
      m_harmonicDE(0.0), 
      m_harmonicFC(0.0), 
      m_exponentialA(0.0), 
      m_exponentialB(0.0), 
      m_exponentialDE(0.0) { Register(); }

    virtual ~ZhuNakamuraCrossing() {};
    virtual const char* getID() { return m_id; };

    virtual ZhuNakamuraCrossing* Clone() { return new ZhuNakamuraCrossing(*this); }

    virtual bool calculateMicroCnlFlux(Reaction* pReact);

    virtual bool calculateCellCrossingCoeffs(Reaction* pReact, std::vector<double>& CrossingProbability);

    virtual bool ThereIsTunnellingWithCrossing(void) { return true; };

  private:

    // One dimensional optimizaton method.
    bool goldenSearchOptimization(AdiabaticCurve* energyCurve, double upperLimit, double lowerLimit, double& optimalR);

    double CalculateAReasonableLimit(AdiabaticCurve*, OneDimensionalFunction*, double);

    // Function for calculating left & right turning points on the Upper & Lower Adiabatic Curves.
    double CalculateTurningPoint(AdiabaticCurve* Curve, double Rint, double TPEnergy, double sgn, double dsgn);

    // These are functions for numerical integration of an appropriately defined OneDimensionalFunction.
    double NumericalIntegration(OneDimensionalFunction*, double a, double b);
    double trapezoidRule(OneDimensionalFunction* func, double a, double b, int n);

    // This function is called by the numerical integration routines.
    void polynomialInterpolation(double xa[], double ya[], int n, double x, double &y, double &dy);

    // This function is a complex Gamma function.
    void complexGamma(double X, double Y, double &GR, double &GI);
    void argumentOfComplexNumber(double X, double Y, double &argument);

    bool ReadDoubleAndUnits2(double& element, PersistPtr pp, const std::string identifier, const std::string units);

    const char* m_id;
    double dataIsInTS;
    double m_H12, m_ReducedMass, m_harmonicX0, m_harmonicDE, m_harmonicFC, m_exponentialA, m_exponentialB, m_exponentialDE;

  };

  //************************************************************
  //Global instance, defining its id
  ZhuNakamuraCrossing theZhuNakamuraCrossing("ZhuNakamuraCrossing");
  //************************************************************

  bool ZhuNakamuraCrossing::ParseData(PersistPtr pp)
  {
    //Look for data in the reaction as child of <me:crossing>

    PersistPtr ppData = pp->XmlMoveTo("me:H12");

    //convert H12 to wavenumbers if not already in wavenumbers
    if (ppData) {

      ReadDoubleAndUnits2(m_H12, pp, "me:H12", "kJ/mol");
      ReadDoubleAndUnits2(m_ReducedMass, pp, "me:reducedMass", "a.m.u.");

      ReadDoubleAndUnits2(m_harmonicX0, pp, "me:harmonicReactantDiabat-X0", "Bohr");
      ReadDoubleAndUnits2(m_harmonicDE, pp, "me:harmonicReactantDiabat-DE", "kJ/mol");
      ReadDoubleAndUnits2(m_harmonicFC, pp, "me:harmonicReactantDiabat-FC", "kJ/mol/Bohr**2");

      ReadDoubleAndUnits2(m_exponentialA, pp, "me:exponentialProductDiabat-A", "kJ/mol");
      ReadDoubleAndUnits2(m_exponentialB, pp, "me:exponentialProductDiabat-B", "1/Bohr");
      ReadDoubleAndUnits2(m_exponentialDE, pp, "me:exponentialProductDiabat-DE", "kJ/mol");
      return true;
    }
    else
    {
      return false;
    }

  }

  bool ZhuNakamuraCrossing::calculateCellCrossingCoeffs(Reaction* pReact, vector<double>& CrossingProbability) {

    parabolicDiabat *reactantCurve;
    exponentialDiabat *productCurve;
    lowerAdiabat *lowerCurve;
    upperAdiabat *upperCurve;
    double delta(0.001), kJPerMoltoHartree(kJPerMol_in_RC / Hartree_in_RC);

    // For Zhu-Nakamura theory, we set the threshold to zero; the routines herein will calculate transmission
    // coefficients appropriate to the specified diabats, which upon convolution with the TS sum of states,
    // will give the correct Energy dependent profile in the k(E)s, and detailed balance will sort it out for
    // the reverse direction.
    //  pReact->set_EffGrnFwdThreshold(0);
    //  pReact->set_EffGrnRvsThreshold(pReact->get);

    // Set up the functions corresponding to the reactant & product diabats.
    reactantCurve = new parabolicDiabat(m_harmonicFC*kJPerMoltoHartree, m_harmonicX0, m_harmonicDE*kJPerMoltoHartree);
    productCurve = new exponentialDiabat(m_exponentialA*kJPerMoltoHartree, m_exponentialB, m_exponentialDE*kJPerMoltoHartree);

    // Set up functions corresponding to the upper & lower adiabatic states.
    lowerCurve = new lowerAdiabat(reactantCurve, productCurve, m_H12*kJPerMoltoHartree);
    upperCurve = new upperAdiabat(reactantCurve, productCurve, m_H12*kJPerMoltoHartree);

    // Find Rx & Ex, the coordinate & energy at which the two diabats cross.
    double Rx, Ex;
    diffBetweenTwoDiabats *diabatDifference;
    diabatDifference = new diffBetweenTwoDiabats(reactantCurve, productCurve);
		// Find Rt on the flipped & translated surface.
    if (!goldenSearchOptimization(diabatDifference, 0.2, -0.2, Rx)) {
      cerr << "problem finding Ex & Rx in Zhu-Nakumara Marquardt routine " << endl;
      return false;
    }
    Ex = reactantCurve->evaluateEnergy(Rx);  // calculate Ex
    if (Ex <= 0.0) {
      cerr << "Zhu-Nakumara Marquardt routine found a negative Et... check your guess value & search range " << endl;
      return false;
    }
		// Test that the energy is more or less identical on both rct & pdt curves.
    double Etest = productCurve->evaluateEnergy(Rx);
    if (abs(Ex - Etest) > (1 / Hartree_in_RC)) {
      cerr << "Zhu-Nakumara Marquardt routine could not find the crossing between the two diabats... exiting " << endl;
      return false;
    }

    // Find Rb & Eb, the coordinate and energy for the bottom of the upper Adiabat,
		// using Rx as the initial guess.
    double Rb, Eb, RbUpperLim, RbLowerLim;

    RbUpperLim = CalculateAReasonableLimit(upperCurve, reactantCurve, Rx); // get a sensible initial guess for the upper limit to Rb
    RbLowerLim = CalculateAReasonableLimit(upperCurve, productCurve, Rx);  // get a sensible initial guess for the lower limit to Rb

    if (!goldenSearchOptimization(upperCurve, RbUpperLim, RbLowerLim, Rb)) {
      cerr << "problem finding Eb & Rb in golden section search routine " << endl;
      return false;
    }
    Eb = upperCurve->evaluateEnergy(Rb);  // calcuate Eb
    if (Eb <= 0.0) {
      cerr << "Zhu-Nakumara Marquardt routine found a negative Eb... check your guess value & search range " << endl;
      return false;
    }

    // Find Et & Rt, the coordinate and energy for the top of the lower Adiabat, using Rx as the initial
		// guess. (To do this we flip the sign of the lower adiabatic curve (turning its max into a min) and
		// translate it up in energy so as to give non-negative values over the search range.)
    double Rt, Et, RtUpperLim, RtLowerLim;

    RtUpperLim = CalculateAReasonableLimit(lowerCurve, productCurve, Rx);  // get a sensible initial guess for the upper limit to Rt
    RtLowerLim = CalculateAReasonableLimit(lowerCurve, reactantCurve, Rx); // get a sensible initial guess for the lower limit to Rt

    flippedTranslatedLowerAdiabat *flippedTranslatedLowerCurve;
    flippedTranslatedLowerCurve = new flippedTranslatedLowerAdiabat(reactantCurve, productCurve, m_H12*kJPerMoltoHartree, Eb);
    if (!goldenSearchOptimization(flippedTranslatedLowerCurve, RtUpperLim, RtLowerLim, Rt)) {  // find Rt on the flipped & translated surface
      cerr << "problem finding Et & Rt in golden section search routine " << endl;
      return false;
    }
    Et = lowerCurve->evaluateEnergy(Rt);
    if (Et <= 0.0) {
      cerr << "Golden section search  routine found a negative Et... check your guess value & search range " << endl;
      return false;
    }

    // Calculate the a & b parameters (in the diabatic representation).

    double Vx = m_H12*kJPerMoltoHartree;         // Calculate coupling in A.U.
    double m = m_ReducedMass*1.822888e+3;        // Calculate mass in A.U.

    double F1 = reactantCurve->NumericalDerivatives(Rx, delta);
    double F2 = productCurve->NumericalDerivatives(Rx, delta);

    double F = sqrt(abs(F1*F2));
    double asqd = F*(F1 - F2) / (16.0*m*pow(Vx, 3.0));
    double a = sqrt(asqd);

    // In Z-N theory, b(E)**2 = (E-Ex)*(F1-F2)/(2*F*Vx) == (F1-F2)/(2*F*Vx)*E - Ex*(F1-F2)/(2*F*Vx)
    double slope = (F1 - F2) / (2 * F*Vx);
    double intercept = -Ex*(F1 - F2) / (2.0*F*Vx);
    linearDiabat *bsqd = new linearDiabat(slope, intercept);

    // Get the properties of the Crossing coefficients vector in which the Z-N probabilities will be placed. 
    const size_t MaximumGrn = pReact->getEnv().MaxGrn;

    double zeroThreshold = m_harmonicDE*kJPerMol_in_RC;
    double EinAU;
    double Lt1(Rt), Lt2(Rt);  // initialize the turning points for the top of the Lower (L) Barrier
    int lastZeroIndex, zone1StartIdx, lastIdx, zone1EndIdx, zone2StartIdx, zone2EndIdx, zone3StartIdx, zone3EndIdx;
    int gsize(pReact->getEnv().GrainSize);
    vector <double> grnCrossingProb(MaximumGrn, 0.0), grnEnergies(MaximumGrn, 0.0);
    double dErctpdt = pReact->get_relative_rctZPE() - pReact->get_relative_pdtZPE();

    for (size_t i(0); i < MaximumGrn; ++i) {
      grnEnergies[i] = double(gsize)*(double(i));
    }

    // Set Prob to zero for any cells which are below either the dE of the reactant diabat or the product zpe
    size_t ii(0);
    double Emin = min(zeroThreshold, dErctpdt);
    ii = size_t(Emin / double(gsize));

    lastZeroIndex = ii - 1;
    zone1StartIdx = ii;
		// Now we run the cell counter up to Et, to determine the range over which we will calculate zone 1 ZN
		// probabilities. This may seem an odd thing to do, but it's computationally much simpler to calculate
		// zone 1 ZN probabilities beginning from Et and then going down in energy, rather than from lower to
		// higher energies gives the last cell index which is set to zero.
		do {
			EinAU = grnEnergies[ii] / Hartree_in_RC;
      lastIdx = ii;
      ++ii;
    } while (bsqd->evaluateEnergy(EinAU) < -1.0 && EinAU < Et && ii < MaximumGrn);

    zone1EndIdx = lastIdx - 1;       // if the do/while loop conditions are satisfied & we have terminated, then we've actually gone one idx too far
    zone2StartIdx = zone1EndIdx + 1;
    do {                              // Now we run the cell counter up to Eb, to determine the range over which we will calculate zone 2 ZN probabilities
      EinAU = grnEnergies[ii] / Hartree_in_RC;
      lastIdx = ii;                  // gives the last cell index which was set to zero
      ++ii;
    } while (EinAU < Eb && ii < MaximumGrn); // we could also enforce the condition *1.0 >= bsqd->evaluateEnergy(EinAU). However this only lines up with
                                              // EinAU < Eb in the case of linear diabats; what's coded here works fine for general curved cases.
    zone2EndIdx = lastIdx - 1;       // if the do/while loop conditions are satisfied & we have terminated, then we've actually gone one idx too far
    zone3StartIdx = zone2EndIdx + 1;
    zone3EndIdx = (MaximumGrn - 1);

    // Now we evaluate the Zone 1 Probabilities, starting from Et, and working our way down in Energy
    // initialize a ftn which is the lower adiabat shifted by dE.
    lowerAdiabatShiftedByE *lowerTPCurve;
    lowerTPCurve = new lowerAdiabatShiftedByE(reactantCurve, productCurve, m_H12*kJPerMoltoHartree, 0);   // initialize the ftn with dE=0
    double leftBound, rightBound, deltaIntegral(0.0), sigma, P12, dummy;
    leftBound = rightBound = Rt;

    cout << "Calculating Zhu Nakamura Probabilities for energies less than Et..." << endl;
    lowerAdiabatActionIntegrand* actionIntegrand;
    actionIntegrand = new lowerAdiabatActionIntegrand(reactantCurve, productCurve, m_H12*kJPerMoltoHartree, m, 0.0);  // initialize the energy to zero for now

    for (int i = zone1EndIdx; i >= zone1StartIdx; --i) {  //calculate the Zhu Nakamura expression for E < Et

      EinAU = grnEnergies[i] /Hartree_in_RC;  // E = double(i) + 0.5;
      lowerTPCurve->setDE(-EinAU); // set dE to the energy where we want to find a turning point (this makes the TP into a minima on the chi-squared function)

			// Calculate Left hand turning point (LHTP). The left hand turning point becomes the right hand bound
			// in the LHTP seeach on the next iteration.
			Lt1 = CalculateTurningPoint(lowerTPCurve, rightBound, EinAU, -1.0, -1.0);
      rightBound = Lt1;

			// Calculate right hand turning point (RHTP). The right hand turning point becomes the left hand bound
			// in the RHTP seeach on the next iteration.
      Lt2 = CalculateTurningPoint(lowerTPCurve, leftBound, EinAU, 1.0, 1.0);
      leftBound = Lt2;

      actionIntegrand->setEinAU(EinAU);                                  // calculate the classical action integrand
      deltaIntegral = NumericalIntegration(actionIntegrand, Lt1, Lt2);   // evaluate the classical action integral
      double g6 = 0.32*pow(10, -2.0/asqd)*exp(-deltaIntegral);
      dummy = sqrt(1.0 - 1.0/pow(bsqd->evaluateEnergy(EinAU),2.0));
      sigma = (M_PI/(16.0*a*sqrt(abs(bsqd->evaluateEnergy(EinAU)))))*sqrt(6.0 + 10.0*dummy)/(1.0 + dummy);
      double sigc = sigma*(1.0-g6);
      double argToB = sigc/M_PI;
      dummy = MesmerGamma(argToB);
      double fB=(2.0*M_PI*pow(argToB,2.0*argToB)*exp(-2.0*argToB))/(argToB*dummy*dummy);
      dummy = fB*exp(-2.0*deltaIntegral);
      double dummy2 = 1.0 + 0.5*a*dummy/(1.0+a);
      P12 = dummy/(dummy2*dummy2 + dummy);
			grnCrossingProb[i] = (IsNan(P12)) ? 0.0 : P12;
		}

    // now we calculate the zone 2 probabilities, beginning at Et and working our way up to Eb
    // comment out Zone 2 ------------------------------------------------------------------------
    cout << "Calculating Zhu Nakamura Probabilities for energies between Et and Eb..." << endl;
    double step, g4, g5, upperLimitOft(0.0), W(0.0), previousP12(0.0);
    ZNzone2Integrand* zone2integrand;
    zone2integrand = new ZNzone2Integrand(0.0, 0.0, 0.0, 0.0, 0.0);  // ftn of the form E(t) = kk*cos(t**3/3.0 + aa*t + bb*t/(cc + 1.4842*t))
                                                                 // all parameters initialized to zero
    for (int i = zone2StartIdx; i <= zone2EndIdx; ++i) {    //calculate the Zhu Nakamura expression for Et < E < Eb

      EinAU = grnEnergies[i] / Hartree_in_RC;

      g5 = 0.38*pow(1 + bsqd->evaluateEnergy(EinAU), 1.2 - 0.4*bsqd->evaluateEnergy(EinAU)) / asqd;
      g4 = ((a - 3.0*bsqd->evaluateEnergy(EinAU)) / (a + 3.0))*sqrt(1.23 + bsqd->evaluateEnergy(EinAU));
      dummy = pow(a, 2.0 / 3.0);
      zone2integrand->setK((1.0 + g5) / dummy);                           // kk = (1.0 + g5)/dummy;                    
      zone2integrand->setA(bsqd->evaluateEnergy(EinAU) / dummy);          // aa = bsqd->evaluateEnergy(EinAU)/dummy; 
      zone2integrand->setB(g4 / dummy);                                   // bb = g4/dummy;
      zone2integrand->setC(0.61*sqrt(2.0 + bsqd->evaluateEnergy(EinAU))); // cc = 0.61*sqrt(2.0+bsqd->evaluateEnergy(EinAU));  
      zone2integrand->setD(pow(a, 1.0 / 3.0));                            // dd = pow(a, 1.0/3.0);                    

		  // The integrand W(t) has very rapid oscillations at longer times. This creates
			// a sign problem for numerical integration & extremely slow convergence. 
			// Convergence is much faster if one always begins at zero, where the integrand
			// has non-oscillatory positive values (prior to the fast oscillations)  
      int ctr(0);
      step = 1.0; upperLimitOft = step;
      P12 = 0.0; W = 0.0;
      do
      {
        previousP12 = P12;                                           
        W = NumericalIntegration(zone2integrand, 0.0, upperLimitOft);
        P12 = W*W / (1 + W*W);                                         
        upperLimitOft += step;                                       
        ++ctr;
      } while (abs((P12 - previousP12) / previousP12) > 0.01);  // converge the Probablities to within 1%

			grnCrossingProb[i] = (IsNan(P12)) ? 0.0 : P12;
		}
    // ------------------------------------------------------------------end comment out of zone 2

    // Now we evaluate the Zone 3 Probabilities, starting from Eb, and working our way up in Energy
    // initialize a ftn which is the lower adiabat shifted by dE.
    upperAdiabatShiftedByE *upperTPCurve;
    upperTPCurve = new upperAdiabatShiftedByE(reactantCurve, productCurve, m_H12*kJPerMoltoHartree, 0);   // initialize the ftn with dE=0
    leftBound = rightBound = Rb;

    cout << "Calculating Zhu Nakamura Probabilities for energies greater than Eb..." << endl;
    upperAdiabatActionIntegrand* upperActionIntegrand;
    upperActionIntegrand = new upperAdiabatActionIntegrand(reactantCurve, productCurve, m_H12*kJPerMoltoHartree, m, 0.0);  // initialize the energy to zero for now

    // Calculate the Zhu Nakamura expression for E > Eb.
    double fctr1 = (0.23*sqrt(a) / (sqrt(a) + 0.75));
    double fctr2 = -0.72 + 0.62*pow(a, 1.43);
    for (int i = zone3StartIdx; i <= zone3EndIdx; ++i) {

			// Set dE to the energy where we want to find a turning point.    
      EinAU = grnEnergies[i] / Hartree_in_RC;
      upperTPCurve->setDE(-EinAU);

      // Calculate Left hand turning point (LHTP). The left hand turning point becomes the right hand bound
			// in the LHTP search on the next iteration.
      Lt1 = CalculateTurningPoint(upperTPCurve, rightBound, EinAU, -1.0, 1.0);
      rightBound = Lt1;

			// Calculate right hand turning point (RHTP). The right hand turning point becomes the left hand bound
			// in the RHTP search on the next iteration.
      Lt2 = CalculateTurningPoint(upperTPCurve, leftBound, EinAU, 1.0, -1.0);
      leftBound = Lt2;

      upperActionIntegrand->setEinAU(EinAU);                          // calculate the classical action integrand
      sigma = NumericalIntegration(upperActionIntegrand, Lt1, Lt2);   // evaluate the classical action integral
      double bsqdval = bsqd->evaluateEnergy(EinAU);
      double g7 = fctr1*pow(40, -sigma);
      dummy = sqrt(1.0 - 1.0 / (bsqdval*bsqdval));
      double delta = (M_PI / (16.0*a*sqrt(abs(bsqdval))))*sqrt(6.0 + 10.0*dummy) / (1.0 + dummy);
			double realPart(0.0), imaginaryPart(0.0), argOfGamma(0.0);
			complexGamma(0.0, delta / M_PI, realPart, imaginaryPart);           // complex Gamma function
      argumentOfComplexNumber(realPart, imaginaryPart, argOfGamma);     // get the argument of the result
      double phi = sigma + delta / M_PI - delta / M_PI*log(delta / M_PI) + M_PI / 4.0 - g7 + argOfGamma;
      dummy = 4.0*cos(phi)*cos(phi);
      double p = exp((-M_PI / (4.0*a))*sqrt(2.0 / (bsqdval + sqrt(bsqdval*bsqdval + fctr2))));
      P12 = dummy / (dummy + p*p / (1.0 - p));
      grnCrossingProb[i] = (IsNan(P12)) ? 0.0 : P12;
    }

    cout << "Finished Calculating Zhu Nakamura Probabilities..." << endl;

		// Interpolate crossing probabilities. 

		const size_t MaximumCell = pReact->getEnv().MaxCell;
		CrossingProbability.clear();
    CrossingProbability.resize(MaximumCell, 0.0);

		Spline spline;
		spline.Initialize(grnEnergies, grnCrossingProb);

		const double cellSize = pReact->getEnv().CellSize;
		for (size_t i(0); i < MaximumCell; ++i) {
			CrossingProbability[i] = spline.Calculate(cellSize*(double(i) + 0.5));
		}

		// Output crossing probabilities if requested. 

		if (pReact->getFlags().CrossingCoeffEnabled) {
      ctest << "\nZhu Nakamura crossing coefficients for: " << pReact->getName() << endl;
      for (size_t i(0); i < MaximumGrn ; ++i) {
        ctest << grnEnergies[i] << "\t" << grnCrossingProb[i] << endl;
      }
      ctest << "}" << endl << "end of Zhu Nakamura crossing coefficients \n";
    }

    return true;
  }

  double ZhuNakamuraCrossing::CalculateAReasonableLimit(AdiabaticCurve* upperOrLowerCurve, OneDimensionalFunction* productOrReactantCurve, double Rx) {

    double Rref, Rlast, dEdR1, dEdR2, trialStep, Rlimit;

    Rref = Rx;      // establish a sensible initial guess for the upper limit to Rb
    do {
      dEdR1 = productOrReactantCurve->NumericalDerivatives(Rref, 0.001);
      dEdR2 = upperOrLowerCurve->NumericalDerivatives(Rref, 0.001);
      trialStep = (upperOrLowerCurve->evaluateEnergy(Rref) - productOrReactantCurve->evaluateEnergy(Rref)) / dEdR1;
      Rlast = Rref;
      Rref += trialStep;
    } while ((dEdR1 / dEdR2) < 0);  // for a given value of R, we assume that we can bracket the upper & lower stationary points (Rt, Eb) 
                                 // if the gradients of the adiatic curve & its corresponding diabatic curve have gradients of the same sign
    Rlimit = Rlast;
    return Rlimit;
  }

  double ZhuNakamuraCrossing::CalculateTurningPoint(AdiabaticCurve* Curve, double Rint, double TPEnergy, double sgn, double dsgn) {

    double dEdR, trialStep(0.001), TP, energy;
    // Establish a sensible initial guess for a  limit to find the turning point.
    double Rref = Rint;
    trialStep *= sgn;

    // First be sure that the gradient is of the correct sign (only really important for regions
    // of very low curvature, e.g., near the minimum or maximum of the function).
    do {
      dEdR = dsgn*Curve->NumericalDerivatives(Rref, 0.001);
      Rref += trialStep;
    } while (dEdR > 0.0);

    do {    // now we check for when the gradient changes sign
      energy = Curve->evaluateEnergy(Rref);
      Rref += trialStep;
    } while (dsgn*sgn*(energy - TPEnergy) > 0.0);

    if (sgn > 0.0) {
      double tmp = Rint;
      Rint = Rref;
      Rref = tmp;
    }

    // Find the turning point on the curve.
    if (!goldenSearchOptimization(Curve, Rint, Rref, TP)) {
      throw(std::runtime_error("CalculateTurningPoint: problem finding Ex & Rx in golden section search routine\n"));
    }

    return TP;
  }

  bool ZhuNakamuraCrossing::goldenSearchOptimization(AdiabaticCurve* energyCurve, double upperLimit, double lowerLimit, double& optimalR)
  {

    if (lowerLimit >= upperLimit) {
      cerr << "There is a problem with either the initialGuess, lowerLimit, or upperLimit you specified for the Zhu-Nakamura Marquardt routine." << endl;
      return false;
    }

    const double Gold = (3.0 - sqrt(5.0)) / 2.0;
    const double tol = 1.0e-8;

    double a = lowerLimit;
    double c = upperLimit;
    double tmp;
    tmp = energyCurve->evaluateEnergy(a);
    double eneA = tmp*tmp;
    tmp = energyCurve->evaluateEnergy(c);
    double eneC = tmp*tmp;

    static const int limit = 20;

    // Locate the initial mid-point by golden section. 
    double b = a + Gold*(c - a);
    tmp = energyCurve->evaluateEnergy(b);
    double eneB = tmp*tmp;

    double x = b + Gold*(c - b);
    tmp = energyCurve->evaluateEnergy(x);
    double eneX = tmp*tmp;

    int count = 0;
    bool converged(false);
    while (count < limit && !converged) {
      count++;

      if (eneX < eneB) {
        a = b;
        eneA = eneB;
        b = x;
        eneB = eneX;

        x = c - Gold*(c - a);
        tmp = energyCurve->evaluateEnergy(x);
        eneX = tmp*tmp;

        converged = (fabs((eneX / eneB) - 1.0) < tol);
      }
      else {
        c = x;
        eneC = eneX;
        x = b;
        eneX = eneB;

        b = a + Gold*(c - a);
        tmp = energyCurve->evaluateEnergy(b);
        eneB = tmp*tmp;

        converged = (fabs((eneB / eneX) - 1.0) < tol);
      }

    }

    // Return the location of the minimum.

    optimalR = (eneX < eneB) ? x : b;

    return true;
  }

  bool ZhuNakamuraCrossing::ReadDoubleAndUnits2(double& element, PersistPtr pp, const std::string identifier, const std::string units)
  {
    PersistPtr ppData = pp->XmlMoveTo(identifier);
    if (!ppData)
      return false;
    element = pp->XmlReadDouble(identifier); //or default
    if (!IsNan(element))
    {
      string unitsTxt = ppData->XmlReadValue("units", false);
      if (unitsTxt != units) {
        cinfo << "MESMER could not read units for " << identifier << "; assuming " << units << "." << endl;
        cinfo << identifier << " = " << element << " " << units << endl;
      }
    }
    else
    {
      cerr << "Spin forbidden crossing: failed to read " << identifier << " (" << units << ")." << endl;
      return false;
    }
    return true;
  }

  // Returns integral of the function func from a to b using Romberg's method
  double ZhuNakamuraCrossing::NumericalIntegration(OneDimensionalFunction* func, double a, double b) {

    static const int JMax(30);
    static const int JMaxP(JMax + 1);
    static const int Rombergk(5);

    double EPS(1.0e-6), ss, dss;
    int K(Rombergk);
    double s[JMaxP + 1];
    double h[JMaxP + 2];

    h[1] = 1.0;
    for (int j(1); j < JMaxP; j++) {
      s[j] = trapezoidRule(func, a, b, j);
      if (j >= K) {
        polynomialInterpolation(&h[j - K], &s[j - K], K, 0.0, ss, dss);
        if (abs(dss) <= EPS*abs(ss)) { return ss; }
      }
      h[j + 1] = 0.25*h[j];
    }
    cerr << "Too many steps in the ZhuNakamuraCrossing NumericalIntegration Routine" << endl;
    return 0.0;
  }

  double ZhuNakamuraCrossing::trapezoidRule(OneDimensionalFunction* func, double a, double b, int n) {
    static double s;

    if (n == 1) {
      return (s = 0.5*(b - a)*(func->evaluateEnergy(a) + func->evaluateEnergy(b)));
    }
    else {
      int it(1);
      for (int j(1); j < n - 1; j++)
        it <<= 1;
      double del = (b - a) / double(it);
      double x = a + 0.5*del;
      double sum = 0.0;
      for (int j(0); j < it; j++, x += del)
        sum += func->evaluateEnergy(x);
      s = 0.5*(s + del*sum);
      return s;
    }
  }

  void ZhuNakamuraCrossing::polynomialInterpolation(double xa[], double ya[], int n, double x, double &y, double &dy) {
    int i, m, ns(1);
    double den, dif, dift, ho, hp, w;
    vector<double> c(n + 1, 0.0), d(n + 1, 0.0);
    dif = abs(x - xa[1]);

    for (i = 1; i <= n; i++) {
      if ((dift = abs(x - xa[i])) < dif) {
        ns = i;
        dif = dift;
      }
      c[i] = ya[i];
      d[i] = ya[i];
    }

    y = ya[ns--];

    for (m = 1; m < n; m++) {
      for (i = 1; i <= n - m; i++) {
        ho = xa[i] - x;
        hp = xa[i + m] - x;
        w = c[i + 1] - d[i];
        if ((den = ho - hp) == 0.0) {
          cerr << "Error in ZhuNakamuraCrossing::polynomialInterpolation routine" << endl;
        }
        den = w / den;
        d[i] = hp*den;
        c[i] = ho*den;
      }
      y += (dy = (2 * ns < (n - m) ? c[ns + 1] : d[ns--]));
    }
  }

  // Function for calculation of complex Gamma - taken from 'Computation of Special functions',
  // translated by DRG (Jun 2015) from F77 to C++. Takes as input the real (X) & imaginary (Y)
  // parts of a complex number and returns the real (GR) and imaginary (GI) parts.
  void ZhuNakamuraCrossing::complexGamma(double X, double Y, double &GR, double &GI) {
    vector <double> A(11, 0.0);
    double X1, Y1, X0, Z1, TH, GR1, GI1, T, TH1, SR, SI, Z2, TH2, G0;
    int NA(0);

    A[1] = 8.333333333333333e-02;
    A[2] = -2.777777777777778e-03;
    A[3] = 7.936507936507937e-04;
    A[4] = -5.952380952380952e-04;
    A[5] = 8.417508417508418e-04;
    A[6] = -1.917526917526918e-03;
    A[7] = 6.410256410256410e-03;
    A[8] = -2.955065359477124e-02;
    A[9] = 1.796443723688307e-01;
    A[10] = -1.39243221690590;

    if (Y == 0.0 && X == int(X) && X <= 0.0) { // gamma is defined for all complex numbers
      GR = 1.0e+300;                       // except the non-positive integers
      GI = 0.0e0;
      cerr << " ZhuNakamuraCrossing::complexGamma: undefined for non-positive integers! " << endl;
      return;
    }
    else if (X < 0.0) {
      X1 = X;
      Y1 = Y;
      X = -1.0*X;
      Y = -1.0*Y;
    }

    X0 = X;
    X1 = 0.0;

    if (X <= 7.0) {
      NA = int(7 - X);
      X0 = X + NA;
    }

    Z1 = sqrt(X0*X0 + Y*Y);
    TH = atan(Y / X0);
    GR = (X0 - 0.5)*log(Z1) - TH*Y - X0 + 0.5*log(2.0*M_PI);
    GI = TH*(X0 - 0.5) + Y*log(Z1) - Y;

    for (int K = 1; K <= 10; ++K) {
      T = pow(Z1, (1.0 - 2.0*double(K)));
      GR = GR + A[K] * T*cos((2.0*K - 1.0)*TH);
      GI = GI - A[K] * T*sin((2.0*K - 1.0)*TH);
    }

    if (X <= 7.0) {
      GR1 = 0.0;
      GI1 = 0.0;
      for (int J = 0; J <= NA - 1; ++J) {
        GR1 = GR1 + 0.5*log((X + J)*(X + J) + Y*Y);
        GI1 = GI1 + atan(Y / (X + J));
      }
      GR = GR - GR1;
      GI = GI - GI1;
    }

    if (X1 < 0.0) {
      Z1 = sqrt(X*X + Y*Y);
      TH1 = atan(Y / X);
      SR = -1.0*sin(M_PI*X)*cosh(M_PI*Y);
      SI = -1.0*cos(M_PI*X)*sinh(M_PI*Y);
      Z2 = sqrt(SR*SR + SI*SI);
      TH2 = atan(SI / SR);
      if (SR < 0.0) { TH2 = M_PI + TH2; }
      GR = log(M_PI / (Z1*Z2)) - GR;
      GI = -1.0*TH1 - TH2 - GI;
      X = X1;
      Y = Y1;
    }

    G0 = exp(GR);
    GR = G0*cos(GI);
    GI = G0*sin(GI);
  }

  void ZhuNakamuraCrossing::argumentOfComplexNumber(double X, double Y, double &argument) {

    if (X == 0.0 && Y == 0.0) {
      throw(std::runtime_error("__FUNCTION__: Real & imaginary parts both equal to zero. Argument is undefined... exiting.\n"));
    }

    if (X > 0.0 || Y != 0.0) {
      argument = 2.0*atan(Y / (sqrt(X*X + Y*Y) + X));
    }
    else if (X < 0.0 && Y > 0.0) {
      argument = M_PI - 2.0*atan(Y / (sqrt(X*X + Y*Y) + fabs(X)));
    }
		else if (X < 0.0 && Y < 0.0) {
			argument = -M_PI - 2.0*atan(Y / (sqrt(X*X + Y*Y) + fabs(X)));
		}
  }

  bool ZhuNakamuraCrossing::calculateMicroCnlFlux(Reaction* pReact)
  {
    Molecule* pTS = pReact->get_TransitionState();
    if (!pTS)
    {
      cerr << "Lack of transition state in reaction " << pReact->getName() << " for Landau-Zener Crossing" << endl;
      return false;
    }

    // Allocate some work space for transistion state density of states.
    vector<double> TScellDOS;
    if (!pTS->getDOS().getCellDensityOfStates(TScellDOS))
      return false;

    // Extract densities of states from molecules.
    const size_t MaximumCell = pReact->getEnv().MaxCell;

    // Allocate space to hold transition state flux and initialize elements to zero.
    vector<double>& rxnFlux = pReact->get_CellFlux();
    rxnFlux.clear();
    rxnFlux.resize(MaximumCell, 0.0);

    vector<double> CrossingProbability;
    calculateCellCrossingCoeffs(pReact, CrossingProbability);

    vector<double> ConvolvedSumOfStates;
    FastLaplaceConvolution(TScellDOS, CrossingProbability, ConvolvedSumOfStates); // FFT convolution

    const size_t BarrierHeight = size_t(max(0.0, pReact->get_ThresholdEnergy()));

    for (size_t i = BarrierHeight; i < MaximumCell; ++i)     // Calculate k(E)s using RRKM expression.
      rxnFlux[i - BarrierHeight] = ConvolvedSumOfStates[i] * SpeedOfLight_in_cm;

    pReact->setCellFluxBottom(pReact->get_relative_rctZPE() + BarrierHeight);

    return true;
  }
}//namespace
