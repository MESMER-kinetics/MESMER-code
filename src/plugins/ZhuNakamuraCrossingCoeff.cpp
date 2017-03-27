//-------------------------------------------------------------------------------------------
//
// Author: Dave Glowacki
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
#include "../dMatrix.h"
#include "AdiabaticCurve.h"
#include "../Rdouble.h"
#include "../gDensityOfStates.h"

#define JMAX 30
#define JMAXP (JMAX+1)
#define ROMBERGK 5

using namespace Constants;
using namespace std;

namespace mesmer
{

  class ZhuNakamuraCrossingCoeff : public MicroRateCalculator
  {
  public:
    virtual bool ParseData(PersistPtr pp);
    // Constructor which registers with the list of CrossingCalculators in the base class
    ZhuNakamuraCrossingCoeff(const char* id) : m_id(id), dataIsInTS(false) { Register(); }

    virtual ~ZhuNakamuraCrossingCoeff() {};
    virtual const char* getID() { return m_id; };

    virtual ZhuNakamuraCrossingCoeff* Clone() { return new ZhuNakamuraCrossingCoeff(*this); }

    virtual bool calculateMicroCnlFlux(Reaction* pReact);

    virtual bool calculateCellCrossingCoeffs(Reaction* pReact, std::vector<double>& CrossingProbability);

    //
    virtual bool ThereIsTunnellingWithCrossing(void) { return true; };

    // this is a Levenburg-Marquardt Algorithm for finding stationary points, crossings, & turning points
    bool MarquardtOptimization(AdiabaticCurve*, double, double, double, double&);

    double CalculateAReasonableLimit(AdiabaticCurve*, OneDimensionalFunction*, double);

    // functions for calculating left & right turning points on the Upper & Lower Adiabatic Curves
    double CalculateLowerLeftTurningPoint(AdiabaticCurve*, double, double);
    double CalculateLowerRightTurningPoint(AdiabaticCurve*, double, double);

    double CalculateUpperLeftTurningPoint(AdiabaticCurve*, double, double);
    double CalculateUpperRightTurningPoint(AdiabaticCurve*, double, double);

    // these are functions for numerical integration of an appropriately defined OneDimensionalFunction
    double NumericalIntegration(OneDimensionalFunction*, double a, double b);
    double trapezoidRule(OneDimensionalFunction* func, double a, double b, int n);

    // this function is called by the numerical integration routines
    void polynomialInterpolation(double xa[], double ya[], int n, double x, double &y, double &dy);

    // this function is a complex Gamma function
    void complexGamma(double X, double Y, double &GR, double &GI);
    void argumentOfComplexNumber(double X, double Y, double &argument);

  private:
    bool ReadDoubleAndUnits2(double& element, PersistPtr pp, const std::string identifier, const std::string units);

    const char* m_id;
    double dataIsInTS;
    double m_H12, m_ReducedMass, m_harmonicX0, m_harmonicDE, m_harmonicFC, m_exponentialA, m_exponentialB, m_exponentialDE;

  };

  //************************************************************
  //Global instance, defining its id
  ZhuNakamuraCrossingCoeff theZhuNakamuraCrossingCoeff("ZhuNakamuraCrossingCoeff");
  //************************************************************

  bool ZhuNakamuraCrossingCoeff::ParseData(PersistPtr pp)
  {
    //Look for data in the reaction as child of <me:crossing>

    PersistPtr ppData = pp->XmlMoveTo("me:H12");

    //convert H12 to wavenumbers if not already in wavenumbers
    if (ppData) {
      /*
          {
            ppData->XmlMoveToProperty("me:H12");
            if (!ppData){
              cerr << "me:H12 is unspecified" << endl;
              return false;
            }
            else{
              const char* units = ppData->XmlReadValue("units", false);
              m_H12 = getConvertedEnergy(units ? units : "cm-1",pp->XmlReadDouble("me:H12"));
          }
      */

      ReadDoubleAndUnits2(m_H12, pp, "me:H12", "kJ/mol");
      ReadDoubleAndUnits2(m_ReducedMass, pp, "me:ReducedMass", "a.m.u.");

      ReadDoubleAndUnits2(m_harmonicX0, pp, "me:harmonicReactantDiabat-X0", "Bohr");
      ReadDoubleAndUnits2(m_harmonicDE, pp, "me:harmonicReactantDiabat-DE", "kJ/mol");
      ReadDoubleAndUnits2(m_harmonicFC, pp, "me:harmonicReactantDiabat-FC", "kJ/mol/Bohr**2");

      ReadDoubleAndUnits2(m_exponentialA, pp, "me:exponentialProductDiabat-A", "kJ/mol");
      ReadDoubleAndUnits2(m_exponentialB, pp, "me:exponentialProductDiabat-B", "Bohr**-1");
      ReadDoubleAndUnits2(m_exponentialDE, pp, "me:exponentialProductDiabat-DE", "kJ/mol");
      return true;
    }
    else
    {
      return false;
    }

  }

  bool ZhuNakamuraCrossingCoeff::calculateCellCrossingCoeffs(Reaction* pReact, vector<double>& CrossingProbability) {

    parabolicDiabat *reactantCurve;
    exponentialDiabat *productCurve;
    lowerAdiabat *lowerCurve;
    upperAdiabat *upperCurve;
    double delta(0.001), kJPerMoltoHartree(kJPerMol_in_RC / Hartree_in_RC);
    vector<double> gradient(1, 0.0);

    // for Zhu-Nakamura theory, we set the threshold to zero; the routines herein will calculate transmission
    // coefficients appropriate to the specified diabats, which upon convolution with the TS sum of states,
    // will give the correct Energy dependent profile in the k(E)s, and detailed balance will sort it out for
    // the reverse direction
    //  pReact->set_EffGrnFwdThreshold(0);
    //  pReact->set_EffGrnRvsThreshold(pReact->get);

/*
    //  test the Marquardt implementatation on a simple quadratic function of the form y = x^2
    parabolicDiabat *testParabola;
    testParabola = new parabolicDiabat(1,0,0);
    double Emin, Rmin;
    MarquardtOptimization(testParabola, -1.0, 5.0, -5.0, Rmin, Emin);
    MarquardtOptimization(testParabola, 1.0, 5.0, -5.0, Rmin, Emin);
*/

// set up the functions corresponding to the reactant & product diabats
    reactantCurve = new parabolicDiabat(m_harmonicFC*kJPerMoltoHartree, m_harmonicX0, m_harmonicDE*kJPerMoltoHartree);
    productCurve = new exponentialDiabat(m_exponentialA*kJPerMoltoHartree, m_exponentialB, m_exponentialDE*kJPerMoltoHartree);

    // set up functions corresponding to the upper & lower adiabatic states
    lowerCurve = new lowerAdiabat(reactantCurve, productCurve, m_H12*kJPerMoltoHartree);
    upperCurve = new upperAdiabat(reactantCurve, productCurve, m_H12*kJPerMoltoHartree);

    // Find Rx & Ex, the coordinate & energy at which the two diabats cross
    double Rx, Ex;
    diffBetweenTwoDiabats *diabatDifference;
    diabatDifference = new diffBetweenTwoDiabats(reactantCurve, productCurve);
    if (!MarquardtOptimization(diabatDifference, 0.0, 0.2, -0.2, Rx)) {  // find Rt on the flipped & translated surface
      cerr << "problem finding Ex & Rx in Zhu-Nakumara Marquardt routine " << endl;
      return false;
    }
    Ex = reactantCurve->evaluateEnergy(Rx);  // calculate Ex
    if (Ex <= 0.0) {
      cerr << "Zhu-Nakumara Marquardt routine found a negative Et... check your guess value & search range " << endl;
      return false;
    }
    double Etest = productCurve->evaluateEnergy(Rx); // test that the energy is more or less identical on both rct & pdt curves
    if (abs(Ex - Etest) > (1 / Hartree_in_RC)) {
      cerr << "Zhu-Nakumara Marquardt routine could not find the crossing between the two diabats... exiting " << endl;
      return false;
    }

    // Find Rb & Eb, the coordinate & enegy for the bottom (b) of the upper Adiabat, using Rx as the initial guess
    double Rb, Eb, RbUpperLim, RbLowerLim;

    RbUpperLim = CalculateAReasonableLimit(upperCurve, reactantCurve, Rx); // get a sensible initial guess for the upper limit to Rb
    RbLowerLim = CalculateAReasonableLimit(upperCurve, productCurve, Rx);  // get a sensible initial guess for the lower limit to Rb

    if (!MarquardtOptimization(upperCurve, Rx, RbUpperLim, RbLowerLim, Rb)) {
      cerr << "problem finding Eb & Rb in Zhu-Nakumara Marquardt routine " << endl;
      return false;
    }
    Eb = upperCurve->evaluateEnergy(Rb);  // calcuate Eb
    if (Eb <= 0.0) {
      cerr << "Zhu-Nakumara Marquardt routine found a negative Eb... check your guess value & search range " << endl;
      return false;
    }

    // Find Et & Rt, the coordinate and energy for the top (t) of the lower Adiabat, using Rx as the initial guess
    //   (to do this we flip the sign of the lower adiabatic curve (turning its max into a min) & translate it up in energy
    //    so as to give non-negative values over the search range)
    double Rt, Et, RtUpperLim, RtLowerLim;

    RtUpperLim = CalculateAReasonableLimit(lowerCurve, productCurve, Rx);  // get a sensible initial guess for the upper limit to Rt
    RtLowerLim = CalculateAReasonableLimit(lowerCurve, reactantCurve, Rx); // get a sensible initial guess for the lower limit to Rt

    flippedTranslatedLowerAdiabat *flippedTranslatedLowerCurve;
    flippedTranslatedLowerCurve = new flippedTranslatedLowerAdiabat(reactantCurve, productCurve, m_H12*kJPerMoltoHartree, Eb);
    if (!MarquardtOptimization(flippedTranslatedLowerCurve, Rx, RtUpperLim, RtLowerLim, Rt)) {  // find Rt on the flipped & translated surface
      cerr << "problem finding Et & Rt in Zhu-Nakumara Marquardt routine " << endl;
      return false;
    }
    Et = lowerCurve->evaluateEnergy(Rt);  // calculate Et
    if (Et <= 0.0) {
      cerr << "Zhu-Nakumara Marquardt routine found a negative Et... check your guess value & search range " << endl;
      return false;
    }

    // calculate the a & b parameters (in the diabatic representation)

    double Vx = m_H12*kJPerMoltoHartree;         // calculate coupling in A.U.
    double m = m_ReducedMass*1.822888e+3;        // calculate mass in A.U

    reactantCurve->setQ(Rx);
    reactantCurve->NumericalDerivatives(delta, gradient);
    double F1 = gradient[0];

    productCurve->setQ(Rx);
    productCurve->NumericalDerivatives(delta, gradient);
    double F2 = gradient[0];

    double F = sqrt(abs(F1*F2));
    double asqd = F*(F1 - F2) / (16.0*m*pow(Vx, 3.0));
    double a = sqrt(asqd);

    // in Z-N theory, b(E)**2 = (E-Ex)*(F1-F2)/(2*F*Vx) == (F1-F2)/(2*F*Vx)*E - Ex*(F1-F2)/(2*F*Vx)
    double slope = (F1 - F2) / (2 * F*Vx);
    double intercept = -1.0*Ex*(F1 - F2) / (2 * F*Vx);
    linearDiabat *bsqd;
    bsqd = new linearDiabat(slope, intercept);

    //
    // Get the properties of the Crossing coefficients vector in which the Z-N probabilities will be placed 
    const int MaximumCell = pReact->getEnv().MaxCell;

    CrossingProbability.clear();
    CrossingProbability.resize(MaximumCell);

    double zeroThreshold = m_harmonicDE*kJPerMol_in_RC;
    double EinAU, transProbability(0.0);
    double Lt1(Rt), Lt2(Rt);  // initialize the turning points for the top of the Lower (L) Barrier
    double Ut1(Rb), Ut2(Rb);  // initialize the turning points for the bottom of the Upper (U) Barrier
    int lastZeroIndex, zone1StartIdx, lastIdx, zone1EndIdx, zone2StartIdx, zone2EndIdx, zone3StartIdx, zone3EndIdx;
    int ctr(0), ii, gsize(pReact->getEnv().GrainSize);
    vector <double> grainCenteredCrossingProbabilities(MaximumCell / gsize, 0.0), grainCenteredEnergies(MaximumCell / gsize, 0.0);;
    double dErctpdt = pReact->get_relative_rctZPE() - pReact->get_relative_pdtZPE();

    double E;
    for (int i = 0; i < MaximumCell / gsize; ++i) {
      E = 0.0;
      for (int jj = 0; jj < gsize; ++jj) {
        E += double(ctr) + 0.5;
        ++ctr;
      }
      grainCenteredEnergies[i] = E / double(gsize);
    }

    ii = 0;
    do {                              // set Prob to zero for any cells which are below either the dE of the reactant diabat or the product zpe
      E = grainCenteredEnergies[ii];          // set E to the avg energy of the cell in cm-1
      EinAU = E / Hartree_in_RC;       // E in A.U.
      grainCenteredCrossingProbabilities[ii] = 0.0;
      lastIdx = ii;                  // gives the last cell index which was set to zero
      ++ii;
    } while (E <= zeroThreshold && E <= dErctpdt && ii < MaximumCell / gsize);

    lastZeroIndex = lastIdx - 1;     // if the do/while loop conditions are satisfied & we have terminated, then we've actually gone one idx too far
    zone1StartIdx = lastZeroIndex + 1;
    do {                              // Now we run the cell counter up to Et, to determine the range over which we will calculate zone 1 ZN probabilities
      E = grainCenteredEnergies[ii];          // this may seem an odd thing to do, but it's computationally much simpler to calculate zone 1 ZN probabilities
      EinAU = E / Hartree_in_RC;       //    beginning from Et and then going down in enegy, rather than from lower to higher energies.
      grainCenteredCrossingProbabilities[ii] = 0.0;
      lastIdx = ii;                  // gives the last cell index which was set to zero
      ++ii;
    } while (bsqd->evaluateEnergy(EinAU) < -1.0 && EinAU < Et && ii < MaximumCell / gsize);

    zone1EndIdx = lastIdx - 1;       // if the do/while loop conditions are satisfied & we have terminated, then we've actually gone one idx too far
    zone2StartIdx = zone1EndIdx + 1;
    do {                              // Now we run the cell counter up to Eb, to determine the range over which we will calculate zone 2 ZN probabilities
      E = grainCenteredEnergies[ii];
      EinAU = E / Hartree_in_RC;
      grainCenteredCrossingProbabilities[ii] = 0.0;
      lastIdx = ii;                  // gives the last cell index which was set to zero
      ++ii;
    } while (EinAU < Eb && ii < MaximumCell / gsize); // we could also enforce the condition *1.0 >= bsqd->evaluateEnergy(EinAU). However this only lines up with
                                              // EinAU < Eb in the case of linear diabats; what's coded here works fine for general curved cases.
    zone2EndIdx = lastIdx - 1;       // if the do/while loop conditions are satisfied & we have terminated, then we've actually gone one idx too far
    zone3StartIdx = zone2EndIdx + 1;
    zone3EndIdx = (MaximumCell / gsize - 1);

    /* // simple tests of the numerical integration routine
    OneDimensionalFunction* testCurve1;
    testCurve1 = new parabolicDiabat(1,0.0,0.0);                      // int(x**2, x=0..4)
    double testResult1 = NumericalIntegration(testCurve1,0.0,4.0);    // analytically equal to 21.33333

    OneDimensionalFunction* testCurve3;
    testCurve3 = new cosineTest2();                                   // int(0.472*cos(R**3/3.0 - 0.0575*r - 0.2436*r/(0.8895 + 1.4842*r)), x=0..10)
    double testResult3 = NumericalIntegration(testCurve3,0.0,10.0);   // MAPLE sets this equal to 0.5784069247
    */

    // Now we evaluate the Zone 1 Probabilities, starting from Et, and working our way down in Energy
    // initialize a ftn which is the lower adiabat shifted by dE; solving Chi-Square for this ftn will give turning pts
    lowerAdiabatShiftedByE *lowerTPCurve;
    lowerTPCurve = new lowerAdiabatShiftedByE(reactantCurve, productCurve, m_H12*kJPerMoltoHartree, 0);   // initialize the ftn with dE=0
    double leftBound, rightBound, deltaIntegral(0.0), sigma, P12, dummy ;
    leftBound = rightBound = Rt;

    cout << "Calculating Zhu Nakamura Probabilities for energies less than Et..." << endl;
    lowerAdiabatActionIntegrand* actionIntegrand;
    actionIntegrand = new lowerAdiabatActionIntegrand(reactantCurve, productCurve, m_H12*kJPerMoltoHartree, m, 0.0);  // initialize the energy to zero for now

    for (int i = zone1EndIdx; i >= zone1StartIdx; --i) {  //calculate the Zhu Nakamura expression for E < Et
/*
      E = grainCenteredEnergies[i];              // E = double(i) + 0.5;
      EinAU = E/Hartree_in_RC;

      lowerTPCurve->setDE(-1*EinAU); // set dE to the energy where we want to find a turning point (this makes the TP into a minima on the chi-squared function)

      Lt1 = CalculateLowerLeftTurningPoint(lowerTPCurve, rightBound, EinAU);   // calculate Left hand turning point (LHTP) by minimizing the chi-squared function
      rightBound = Lt1;           // the left hand turning point becomes the right hand bound in the LHTP seeach on the next iteration

      Lt2 = CalculateLowerRightTurningPoint(lowerTPCurve, leftBound, EinAU);   // calculate right hand turning point (RHTP) by minimizing the chi-squared function
      leftBound = Lt2;            // the right hand turning point becomes the left hand bound in the RHTP seeach on the next iteration

      actionIntegrand->setEinAU(EinAU);                                  // calculate the classical action integrand
      deltaIntegral = NumericalIntegration(actionIntegrand, Lt1, Lt2);   // evaluate the classical action integral
      g6 = 0.32*pow(10, -2.0/asqd)*exp(-1.0*deltaIntegral);
      dummy = sqrt(1.0 - 1.0/pow(bsqd->evaluateEnergy(EinAU),2.0));
      sigma = (M_PI/(16.0*a*sqrt(abs(bsqd->evaluateEnergy(EinAU)))))*sqrt(6.0 + 10.0*dummy)/(1.0 + dummy);
      sigc = sigma*(1.0-g6);
      argToB = sigc/M_PI;
      dummy = MesmerGamma(argToB);
      fB=(2.0*M_PI*pow(argToB,2.0*argToB)*exp(-2.0*argToB))/(argToB*dummy*dummy);
      dummy = fB*exp(-2.0*deltaIntegral);
      dummy2 = 1.0 + 0.5*a*dummy/(1.0+a);
      P12 = dummy/(dummy2*dummy2 + dummy);
//      cout << E << " Lt1 " << Lt1 << " Lt2 " << Lt2 << " " << EinAU/lowerCurve->evaluateEnergy(Lt1)
//                                                    << " " << EinAU/lowerCurve->evaluateEnergy(Lt2) << " " << P12 << endl;
      grainCenteredCrossingProbabilities[i] = P12;
      if(IsNan(grainCenteredCrossingProbabilities[i])) grainCenteredCrossingProbabilities[i] = 0.0;
*/
      grainCenteredCrossingProbabilities[i] = 0.0;
    }

    // now we calculate the zone 2 probabilities, beginning at Et and working our way up to Eb
    // comment out Zone 2 ------------------------------------------------------------------------
    cout << "Calculating Zhu Nakamura Probabilities for energies between Et and Eb..." << endl;
    double step, g4, g5, upperLimitOft(0.0), W(0.0), previousP12(0.0);
    ZNzone2Integrand* zone2integrand;
    zone2integrand = new ZNzone2Integrand(0.0, 0.0, 0.0, 0.0, 0.0);  // ftn of the form E(t) = kk*cos(t**3/3.0 + aa*t + bb*t/(cc + 1.4842*t))
                                                                 // all parameters initialized to zero
    for (int i = zone2StartIdx; i <= zone2EndIdx; ++i) {    //calculate the Zhu Nakamura expression for Et < E < Eb

      E = grainCenteredEnergies[i];         // E = double(i) + 0.5;         
      EinAU = E / Hartree_in_RC;

      g5 = 0.38*pow(1 + bsqd->evaluateEnergy(EinAU), 1.2 - 0.4*bsqd->evaluateEnergy(EinAU)) / asqd;
      g4 = ((a - 3.0*bsqd->evaluateEnergy(EinAU)) / (a + 3.0))*sqrt(1.23 + bsqd->evaluateEnergy(EinAU));
      dummy = pow(a, 2.0 / 3.0);
      zone2integrand->setK((1.0 + g5) / dummy);                           // kk = (1.0 + g5)/dummy;                    
      zone2integrand->setA(bsqd->evaluateEnergy(EinAU) / dummy);          // aa = bsqd->evaluateEnergy(EinAU)/dummy; 
      zone2integrand->setB(g4 / dummy);                                   // bb = g4/dummy;
      zone2integrand->setC(0.61*sqrt(2.0 + bsqd->evaluateEnergy(EinAU))); // cc = 0.61*sqrt(2.0+bsqd->evaluateEnergy(EinAU));  
      zone2integrand->setD(pow(a, 1.0 / 3.0));                            // dd = pow(a, 1.0/3.0);                    

      int ctr(0);
      step = 1.0; upperLimitOft = step;
      P12 = 0.0; W = 0.0;
      do
      {
        previousP12 = P12;                                           // The integrand W(t) has very rapid oscillations at longer times. This creates
        W = NumericalIntegration(zone2integrand, 0.0, upperLimitOft);// a sign problem for numerical integration & extremely slow convergence. 
        P12 = W*W / (1 + W*W);                                           // Convergence is much faster if one always begins at zero, where the integrand  
        upperLimitOft += step;                                       // has non-oscillatory positive values (prior to the fast oscillations)  
        ++ctr;
      } while (abs((P12 - previousP12) / previousP12) > 0.01);  // converge the Probablities to within 1%

//      cout << "E " << E << " P12 " << P12 << " ctr " << ctr << endl;
      grainCenteredCrossingProbabilities[i] = P12;
      if (IsNan(grainCenteredCrossingProbabilities[i])) grainCenteredCrossingProbabilities[i] = 0.0;
    }
    // ------------------------------------------------------------------end comment out of zone 2

    // Now we evaluate the Zone 3 Probabilities, starting from Eb, and working our way up in Energy
    // initialize a ftn which is the lower adiabat shifted by dE; solving Chi-Square for this ftn will give turning pts
    upperAdiabatShiftedByE *upperTPCurve;
    upperTPCurve = new upperAdiabatShiftedByE(reactantCurve, productCurve, m_H12*kJPerMoltoHartree, 0);   // initialize the ftn with dE=0
    double g7, argOfGamma, phi, p, bsqdval, realPart, complexPart;
    leftBound = rightBound = Rb;

    cout << "Calculating Zhu Nakamura Probabilities for energies greater than Eb..." << endl;
    upperAdiabatActionIntegrand* upperActionIntegrand;
    upperActionIntegrand = new upperAdiabatActionIntegrand(reactantCurve, productCurve, m_H12*kJPerMoltoHartree, m, 0.0);  // initialize the energy to zero for now

    for (int i = zone3StartIdx; i <= zone3EndIdx; ++i) {  //calculate the monstrosity that is the Zhu Nakamura expression for E > Eb

      E = grainCenteredEnergies[i];         //      E = double(i) + 0.5;         
      EinAU = E / Hartree_in_RC;
      realPart = 0.0; complexPart = 0.0; argOfGamma = 0.0;

      upperTPCurve->setDE(-1 * EinAU); // set dE to the energy where we want to find a turning point (this makes the TP into a minima on the chi-squared function)    

      Lt1 = CalculateUpperLeftTurningPoint(upperTPCurve, rightBound, EinAU);   // calculate Left hand turning point (LHTP) by minimizing the chi-squared function
      rightBound = Lt1;           // the left hand turning point becomes the right hand bound in the LHTP seeach on the next iteration

      Lt2 = CalculateUpperRightTurningPoint(upperTPCurve, leftBound, EinAU);   // calculate right hand turning point (RHTP) by minimizing the chi-squared function
      leftBound = Lt2;            // the right hand turning point becomes the left hand bound in the RHTP seeach on the next iteration

      sigma = 0.0;
      upperActionIntegrand->setEinAU(EinAU);                          // calculate the classical action integrand
      sigma = NumericalIntegration(upperActionIntegrand, Lt1, Lt2);   // evaluate the classical action integral
      bsqdval = bsqd->evaluateEnergy(EinAU);
      g7 = (0.23*sqrt(a) / (sqrt(a) + 0.75))*pow(40, -1.0*sigma);
      dummy = sqrt(1.0 - 1.0 / pow(bsqdval, 2.0));
      delta = (M_PI / (16.0*a*sqrt(abs(bsqdval))))*sqrt(6.0 + 10.0*dummy) / (1.0 + dummy);
      complexGamma(0.0, delta / M_PI, realPart, complexPart);           // complex Gamma function
      argumentOfComplexNumber(realPart, complexPart, argOfGamma);     // get the argument of the result
      phi = sigma + delta / M_PI - delta / M_PI*log(delta / M_PI) + M_PI / 4.0 - g7 + argOfGamma;
      dummy = 4.0*cos(phi)*cos(phi);
      p = exp((-M_PI / (4.0*a))*sqrt(2.0 / (bsqdval + sqrt(bsqdval*bsqdval - 0.72 + 0.62*pow(a, 1.43)))));
      P12 = dummy / (dummy + p*p / (1.0 - p));
      //      cout << "Zone3 E: " << E << " Lt1 " << Lt1 << " Lt2 " << Lt2 << " P12 " << P12 << endl;
      grainCenteredCrossingProbabilities[i] = P12;
      if (IsNan(grainCenteredCrossingProbabilities[i])) grainCenteredCrossingProbabilities[i] = 0.0;
    }

    cout << "Finished Calculating Zhu Nakamura Probabilities..." << endl;

    ctr = 0;
    double E0, dP12, P120;
    for (int i = 0; i < grainCenteredEnergies[0]; ++i) {
      CrossingProbability[i] = 0.0;
      ++ctr;
    }
    for (int i = 0; i < MaximumCell / gsize - 1; ++i) {
      E0 = ctr;
      P120 = grainCenteredCrossingProbabilities[i];
      dP12 = grainCenteredCrossingProbabilities[i + 1] - grainCenteredCrossingProbabilities[i];
      for (int jj = 1; jj <= gsize; ++jj) {
        CrossingProbability[ctr] = (dP12 / gsize)*(jj)+P120;
        ++ctr;
      }
    }
    for (int i = ctr; i < MaximumCell; ++i) {
      CrossingProbability[i] = 0.0;
    }


    //    for(int i=zone2StartIdx; i<=zone3EndIdx; ++i){
    //      CrossingProbability[i] = 1.0;
    //    }

    //      CrossingProbability[i] = transProbability;
    //      CrossingProbability[i] = (1.0e+0 + trans_probability)*(1.0e+0 - trans_probability);
    //      if(IsNan(CrossingProbability[i])) CrossingProbability[i] = 0.0;


    if (pReact->getFlags().CrossingCoeffEnabled) {
      ctest << "\nZhu Nakamura crossing coefficients for: " << pReact->getName() << endl;
      for (int i = 0; i < MaximumCell / gsize; ++i) {
        ctest << grainCenteredEnergies[i] << "\t" << grainCenteredCrossingProbabilities[i] << endl;
      }
      ctest << "}" << endl << "end of Zhu Nakamura crossing coefficients \n";
    }

    return true;
  }

  double ZhuNakamuraCrossingCoeff::CalculateAReasonableLimit(AdiabaticCurve* upperOrLowerCurve, OneDimensionalFunction* productOrReactantCurve, double Rx) {

    double Rref, Rlast, dEdR1, dEdR2, trialStep, Rlimit;
    vector<double> gradient(1, 0.0);

    Rref = Rx;      // establish a sensible initial guess for the upper limit to Rb
    do {
      productOrReactantCurve->setQ(Rref);
      productOrReactantCurve->NumericalDerivatives(0.001, gradient);
      dEdR1 = gradient[0];
      upperOrLowerCurve->setQ(Rref);
      upperOrLowerCurve->NumericalDerivatives(0.001, gradient);
      dEdR2 = gradient[0];
      trialStep = (upperOrLowerCurve->evaluateEnergy(Rref) - productOrReactantCurve->evaluateEnergy(Rref)) / dEdR1;
      Rlast = Rref;
      Rref += trialStep;
    } while ((dEdR1 / dEdR2) < 0);  // for a given value of R, we assume that we can bracket the upper & lower stationary points (Rt, Eb) 
                                 // if the gradients of the adiatic curve & its corresponding diabatic curve have gradients of the same sign
    Rlimit = Rlast;
    return Rlimit;
  }

  double ZhuNakamuraCrossingCoeff::CalculateLowerLeftTurningPoint(AdiabaticCurve* LowerCurve, double Rright, double TPEnergy) {

    double dEdR, leftBound, Guess, trialStep, TP, energy;
    vector<double> gradient(1, 0.0);
    double Rref = Rright;      // establish a sensible initial guess for a lower limit to find the left hand turning point

    do {                              // first be sure that the gradient is of the correct sign (only really important for regions of very low
      LowerCurve->setQ(Rref);         //   curvature, e.g., near the minimum or maximum of the function)
      LowerCurve->NumericalDerivatives(0.001, gradient);
      dEdR = gradient[0];
      trialStep = 0.001;
      Rref -= trialStep;
    } while (dEdR < 0.0);

    do {    // now we check for when the gradient changes sign
      LowerCurve->setQ(Rref);
      energy = LowerCurve->evaluateEnergy(Rref);
      trialStep = 0.001;
      Rref -= trialStep;
      //      cout << energy << " " << TPEnergy;
    } while (energy > TPEnergy);

    leftBound = Rref;
    Guess = Rright - (Rright - leftBound) / 2.0;

    if (!MarquardtOptimization(LowerCurve, Guess, Rright, leftBound, TP)) {            // find the left hand turning point on the lower surface
      cerr << "CalculateLowerLeftTurningPoint: problem finding Ex & Rx in Zhu-Nakumara Marquardt routine " << endl;
      exit(1);
    }

    return TP;
  }

  double ZhuNakamuraCrossingCoeff::CalculateLowerRightTurningPoint(AdiabaticCurve* LowerCurve, double Rleft, double TPEnergy) {

    double dEdR, rightBound, Guess, trialStep, TP, energy;
    vector<double> gradient(1, 0.0);
    double Rref = Rleft;      // establish a sensible initial guess for a lower limit to find the left hand turning point

    do {                              // first be sure that the gradient is of the correct sign (only really important for regions of very low
      LowerCurve->setQ(Rref);         //   curvature, e.g., near the minimum or maximum of the function)
      LowerCurve->NumericalDerivatives(0.001, gradient);
      dEdR = gradient[0];
      trialStep = 0.001;
      Rref += trialStep;
    } while (dEdR > 0.0);

    do {    // now we check for when the gradient changes sign
      LowerCurve->setQ(Rref);
      energy = LowerCurve->evaluateEnergy(Rref);
      trialStep = 0.001;
      Rref += trialStep;
    } while (energy > TPEnergy);

    rightBound = Rref;
    Guess = Rleft - (Rleft - rightBound) / 2.0;

    if (!MarquardtOptimization(LowerCurve, Guess, rightBound, Rleft, TP)) {            // find the left hand turning point on the lower surface
      cerr << "CalculateLowerRightTurningPoint: problem finding Ex & Rx in Zhu-Nakumara Marquardt routine " << endl;
      exit(1);
    }

    return TP;
  }

  double ZhuNakamuraCrossingCoeff::CalculateUpperLeftTurningPoint(AdiabaticCurve* UpperCurve, double Rright, double TPEnergy) {

    double dEdR, leftBound, Guess, trialStep, TP, energy;
    vector<double> gradient(1, 0.0);
    double Rref = Rright;      // establish a sensible initial guess for a lower limit to find the left hand turning point

    do {                       // first be sure that the gradient is of the correct sign (only really important for regions of very low
      UpperCurve->setQ(Rref);  //   curvature, e.g., near the minimum or maximum of the function)
      UpperCurve->NumericalDerivatives(0.001, gradient);
      dEdR = gradient[0];
      trialStep = 0.001;
      Rref -= trialStep;
    } while (dEdR > 0.0);

    do {    // now we check for when the point we're at crosses the energy threshold corresponding to the turning point we want
      UpperCurve->setQ(Rref);
      energy = UpperCurve->evaluateEnergy(Rref);
      trialStep = 0.001;
      Rref -= trialStep;
    } while (energy < TPEnergy);

    leftBound = Rref;
    Guess = Rright - (Rright - leftBound) / 2.0;

    if (!MarquardtOptimization(UpperCurve, Guess, Rright, leftBound, TP)) {            // find the left hand turning point on the lower surface
      cerr << "CalculateUpperLeftTurningPoint: problem finding Ex & Rx in Zhu-Nakumara Marquardt routine " << endl;
      exit(1);
    }

    return TP;
  }

  double ZhuNakamuraCrossingCoeff::CalculateUpperRightTurningPoint(AdiabaticCurve* UpperCurve, double Rleft, double TPEnergy) {

    double dEdR, rightBound, Guess, trialStep, TP, energy;
    vector<double> gradient(1, 0.0);
    double Rref = Rleft;      // establish a sensible initial guess for a lower limit to find the left hand turning point

    do {                      // first be sure that the gradient is of the correct sign (only really important for regions of very low
      UpperCurve->setQ(Rref); //   curvature, e.g., near the minimum or maximum of the function
      UpperCurve->NumericalDerivatives(0.001, gradient);
      dEdR = gradient[0];
      trialStep = 0.001;
      Rref += trialStep;
    } while (dEdR < 0.0);

    do {    // now we check for when the gradient changes sign
      UpperCurve->setQ(Rref);
      energy = UpperCurve->evaluateEnergy(Rref);
      trialStep = 0.001;
      Rref += trialStep;
    } while (energy < TPEnergy);

    rightBound = Rref;
    Guess = Rleft - (Rleft - rightBound) / 2.0;

    if (!MarquardtOptimization(UpperCurve, Guess, rightBound, Rleft, TP)) {            // find the left hand turning point on the lower surface
      cerr << "CalculateUpperRightTurningPoint: problem finding Ex & Rx in Zhu-Nakumara Marquardt routine " << endl;
      exit(1);
    }

    return TP;
  }

  bool ZhuNakamuraCrossingCoeff::MarquardtOptimization
  (AdiabaticCurve* energyCurve, double initialGuess, double upperLimit, double lowerLimit, double& optimalR)
  {

    if (initialGuess > upperLimit || initialGuess < lowerLimit || lowerLimit >= upperLimit) {
      cerr << "There is a problem with either the initialGuess, lowerLimit, or upperLimit you specified for the Zhu-Nakamura Marquardt routine." << endl;
      return false;
    }

    int nVar = 1;

    // Uncomment to enable ctest output during fitting. Or use -w5 option in command.
    //ChangeErrorLevel e(obDebug); 

    //Default is to disable ctest during fitting. Restored when leaving this function.
//		StopCTestOutput stop(true) ;

    vector<double> currentLocation(nVar, 0.0); // I've left these as vectors shoul dwe want to generalize this to multidimensional ftns in the future
    vector<double> newLocation(nVar, 0.0);     // for Zhu-Nakamura 1d curve optimization, it's not really necessary - DRG 4 Jun 2015

    // Set the Current Location to the initial Guess
    currentLocation[0] = initialGuess;
    energyCurve->setQ(currentLocation[0]);

    double chiSquare(0.0), lambda(0.001);

    energyCurve->setQ(currentLocation[0]);
    energyCurve->chiSquareCalculation(chiSquare);

    double bestChiSquare = chiSquare;

    //
    // The following is slightly modified implementation of the Marquardt
    // algorithm that copes with optimization bounds, copied from Marquardt.cpp 
    //

    double delta = 0.0001;
    size_t m_maxIterations = 1000;
    double tol = 1 / (10.0*Hartree_in_RC);
    double lambdaScale = 10.0;

    vector<double> gradient(nVar, 0.0);
    qdMatrix hessian(nVar, 0.0);
    energyCurve->NumericalDerivativesOfChiSquare(delta, gradient);
    energyCurve->NumericalHessianOfChiSquare(delta, hessian);

    bool converged(false);
    for (size_t itr(1); itr <= m_maxIterations && !converged; itr++) {

      newLocation = currentLocation;
      vector<qd_real> deltaLocation(nVar, 0.0);
      qdMatrix invHessian = hessian;

      for (int iVar = 0; iVar < nVar; iVar++) {
        invHessian[iVar][iVar] *= 0.5*qd_real(1.0 + lambda);
        deltaLocation[iVar] = qd_real(-0.5*gradient[iVar]);
      }

      invHessian.solveLinearEquationSet(&deltaLocation[0]);

      for (int iVar = 0; iVar < nVar; iVar++) {
        newLocation[iVar] += to_double(deltaLocation[iVar]);
      }

      // Check bounds.    
      if (newLocation[0] <= upperLimit && newLocation[0] >= lowerLimit) {

        energyCurve->setQ(newLocation[0]);
        energyCurve->chiSquareCalculation(chiSquare);

        if (chiSquare > bestChiSquare) {
          lambda *= lambdaScale;
          energyCurve->setQ(newLocation[0]);
        }
        else {
          double relativeChange = 1.0 - chiSquare / bestChiSquare;
          converged = (relativeChange < tol);
          lambda /= lambdaScale;
          currentLocation[0] = energyCurve->getQ();
          bestChiSquare = chiSquare;
          energyCurve->NumericalDerivativesOfChiSquare(delta, gradient);
          energyCurve->NumericalHessianOfChiSquare(delta, hessian);
        }
      }
      else {
        lambda *= lambdaScale;
      }
      optimalR = currentLocation[0];
      cinfo << "Iteration: " << itr << " of Marquardt. ChiSquare = " << bestChiSquare << ", Lambda = " << lambda << endl;
    }
    return true;
  }

  bool ZhuNakamuraCrossingCoeff::ReadDoubleAndUnits2(double& element, PersistPtr pp, const std::string identifier, const std::string units)
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
  double ZhuNakamuraCrossingCoeff::NumericalIntegration(OneDimensionalFunction* func, double a, double b) {
    double EPS(1.0e-6);
    int K(ROMBERGK), j;
    double ss, dss;
    double s[JMAXP + 1];
    double h[JMAXP + 2];

    h[1] = 1.0;
    for (j = 1; j <= JMAX; j++) {
      s[j] = trapezoidRule(func, a, b, j);
      if (j >= K) {
        polynomialInterpolation(&h[j - K], &s[j - K], K, 0.0, ss, dss);
        if (abs(dss) <= EPS*abs(ss)) { return ss; }
      }
      h[j + 1] = 0.25*h[j];
    }
    cerr << "Too many steps in the ZhuNakamuraCrossingCoeff NumericalIntegration Routine" << endl;
    return 0.0;
  }

  double ZhuNakamuraCrossingCoeff::trapezoidRule(OneDimensionalFunction* func, double a, double b, int n) {
    double x, tnm, sum, del;
    static double s;
    int it, j;

    if (n == 1) { return (s = 0.5*(b - a)*(func->evaluateEnergy(a) + func->evaluateEnergy(b))); }
    else {
      for (it = 1, j = 1; j < n - 1; j++) it <<= 1;
      tnm = it;
      del = (b - a) / tnm;
      x = a + 0.5*del;
      sum = 0.0;
      for (j = 1; j <= it; j++) {
        x += del;
        sum += func->evaluateEnergy(x);
      }
      s = 0.5*(s + (b - a)*sum / tnm);
      return s;
    }
  }

  void ZhuNakamuraCrossingCoeff::polynomialInterpolation(double xa[], double ya[], int n, double x, double &y, double &dy) {
    int i, m, ns(1);
    double den, dif, dift, ho, hp, w;
    double c[ROMBERGK + 1], d[ROMBERGK + 1];     //    vector <double> c(n,0.0), d(n,0.0);
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
          cerr << "Error in ZhuNakamuraCrossingCoeff::polynomialInterpolation routine" << endl;
        }
        den = w / den;
        d[i] = hp*den;
        c[i] = ho*den;
      }
      y += (dy = (2 * ns < (n - m) ? c[ns + 1] : d[ns--]));
    }
  }

  // function for calculation of complex Gamma - taken from 'Computation of Special functions',
	// translated by DRG (Jun 2015) from F77 to C++
  // takes as input the real (X) & imaginary parts (Y) of a complex number
  // and returns the real (GR) and imaginary (GI) parts
  void ZhuNakamuraCrossingCoeff::complexGamma(double X, double Y, double &GR, double &GI) {
    vector <double> A(11, 0.0);
    double X1, Y1, X0, Z1, TH, GR1, GI1, T, TH1, SR, SI, Z2, TH2, G0;
    int NA;

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
      cerr << " ZhuNakamuraCrossingCoeff::complexGamma: undefined for non-positive integers! " << endl;
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

  void ZhuNakamuraCrossingCoeff::argumentOfComplexNumber(double X, double Y, double &argument) {
    if (X > 0.0 || Y != 0.0) {
      argument = 2.0*atan(Y / (sqrt(X*X + Y*Y) + X));
    }
    else if (X < 0.0 && Y == 0.0) {
      argument = M_PI;
    }
    else if (X == 0.0 && Y == 0.0) {
      cerr << "ZhuNakamuraCrossingCoeff::argumentOfComplexNumber: real & imaginary parts both equal to zero. Argument is undefined... exiting." << endl;
      exit(1);
    }
  }
  bool ZhuNakamuraCrossingCoeff::calculateMicroCnlFlux(Reaction* pReact)
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
