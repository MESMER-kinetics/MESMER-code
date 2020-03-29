#ifndef GUARD_AdiabaticCurve_h
#define GUARD_AdiabaticCurve_h

#include <cmath>

namespace mesmer
{
  // Define an abstract base class that holds any of a range of one dimensional functions, 
  // each of which returns at energy for some input coordinate r.
  class OneDimensionalFunction {

  private:
    double Energy;

  public:
    OneDimensionalFunction(): Energy(0.0) {};
    virtual ~OneDimensionalFunction() {};

    // set Functions
    void setEnergy(double e) { Energy = e; };

    // get Functions
    double getEnergy() { return Energy; };

    // abstract virtual Energy evaluation function
    virtual double evaluateEnergy(double) = 0;

    // Function for finite central difference numerical derivative.
    double NumericalDerivatives(double r, double delta) {
      return (evaluateEnergy(r + 0.5*delta) - evaluateEnergy(r - 0.5*delta)) / delta;
    };

  };

  // exponentialDiabat is a diabatic energy curve of the form A*exp(B*r)+DE
  class exponentialDiabat : public OneDimensionalFunction {

  private:
    double A, B, DE;

  public:
    // Constructor
    exponentialDiabat(double a, double b, double dE): A(a), B(b), DE(dE) {};

    // Destructor
    virtual ~exponentialDiabat() {};

    double evaluateEnergy(double r) {
      double Energy = A*exp(B*r) + DE;
      setEnergy(Energy);
      return Energy;
    };

  };

  // parabolicDiabat is a diabatic energy curve of the form K*(r-X0)**2+DE
  class parabolicDiabat : public OneDimensionalFunction {

  private:
    double K, X0, DE;

  public:
    // Constructor
    parabolicDiabat(double k, double x0, double dE): K(k), X0(x0), DE(dE) {};

    // Destructor
    virtual ~parabolicDiabat() {};

    double evaluateEnergy(double r) {
      double Energy = K*(r - X0)*(r - X0) + DE;
      setEnergy(Energy);
      return Energy;
    };

  };

  // linearDiabat is a diabatic energy curve of the form K*r+DE
  class linearDiabat : public OneDimensionalFunction {

  private:
    double K, DE;

  public:
    // Constructor
    linearDiabat(double k, double dE): K(k), DE(dE) {};

    // Destructor
    virtual ~linearDiabat() {};

    double evaluateEnergy(double r) {
      double Energy = K*r + DE;
      setEnergy(Energy);
      return Energy;
    };
  };

  // another test function
  class cosineTest2 : public OneDimensionalFunction {

  public:
    // Constructor
    cosineTest2() {};

    // Destructor
    virtual ~cosineTest2() {};

    double evaluateEnergy(double r) {
      double Energy = 0.472*cos((r*r*r) / 3.0 - 0.0575*r - 0.2436*r / (0.8895 + 1.4842*r));
      setEnergy(Energy);
      return Energy;
    };
  };

  // ftn of the form E(t) = K*cos[t**3/3 - A*t - B*t/(C + D*t)], which are the sorts 
  // of functions required to specify the integrand in zone 2 of Zhu-Nakamura Theory
  class ZNzone2Integrand : public OneDimensionalFunction {

  private:
    double K, A, B, C, D;

  public:
    // Constructor
    ZNzone2Integrand(double k, double a, double b, double c, double d) {
      K = k; A = a; B = b; C = c; D = d;
    };

    // Destructor
    virtual ~ZNzone2Integrand() {};

    void setK(double k) { K = k; };
    void setA(double a) { A = a; };
    void setB(double b) { B = b; };
    void setC(double c) { C = c; };
    void setD(double d) { D = d; };

    double evaluateEnergy(double r) {
      double Energy = K*cos((r*r*r) / 3.0 - A*r - B*r / (C + D*r));
      setEnergy(Energy);
      return Energy;
    };
  };

  // Define an abstract base class that holds any of a range of adiabatic curves.
  class AdiabaticCurve : public OneDimensionalFunction {

  private:
    OneDimensionalFunction *R;      // reactant diabat
    OneDimensionalFunction *P;      // product diabat
    double H12; // coordinate, coupling, & Energy

  public:
    AdiabaticCurve(): R(NULL), P(NULL), H12(0.0) {};
    virtual ~AdiabaticCurve() {};

    // set Functions
    void setR(OneDimensionalFunction* rct) { R = rct; };
    void setP(OneDimensionalFunction* pdt) { P = pdt; };
    void setH12(double v12) { H12 = v12; };

    // get Functions
    OneDimensionalFunction* getR() { return R; };
    OneDimensionalFunction* getP() { return P; };
    double getH12() { return H12; };

  };

  // derived class to for the lower adiabatic curve obtained from two coupled diabats
  class lowerAdiabat : public AdiabaticCurve {

  public:
    // Constructor
    lowerAdiabat(OneDimensionalFunction* rct, OneDimensionalFunction* pdt, double v12) {
      setR(rct);  setP(pdt);  setH12(v12);
    };

    // Destructor
    virtual ~lowerAdiabat() {};

    // Energy evaluation routine
    double evaluateEnergy(double r) {
      double V1 = getR()->evaluateEnergy(r);
      double V2 = getP()->evaluateEnergy(r);
      double H12 = getH12();
      double Energy = 0.5*(V1 + V2 - sqrt((V1 - V2)*(V1 - V2) + 4.0*H12*H12));
      setEnergy(Energy);
      return Energy;
    };

  };

  // derived class for the upper adiabatic curve calculated from two coupled diabats
  class upperAdiabat : public AdiabaticCurve {

  public:
    // Constructor
    upperAdiabat(OneDimensionalFunction* rct, OneDimensionalFunction* pdt, double v12) {
      setR(rct);  setP(pdt);  setH12(v12);
    };

    // Destructor
    virtual ~upperAdiabat() {};

    // Energy evaluation routine
    double evaluateEnergy(double r) {
      double V1 = getR()->evaluateEnergy(r);
      double V2 = getP()->evaluateEnergy(r);
      double H12 = getH12();
      double Energy = 0.5*(V1 + V2 + sqrt((V1 - V2)*(V1 - V2) + 4.0*H12*H12));
      setEnergy(Energy);
      return Energy;
    };

  };

  // derived class in which the lower adiabatic is flipped (mult by -1) & translated up in Energy space
  class flippedTranslatedLowerAdiabat : public AdiabaticCurve {

  private:
    double DE;

  public:
    // Constructor
    flippedTranslatedLowerAdiabat(OneDimensionalFunction* rct, OneDimensionalFunction* pdt, double v12, double dE) {
      setR(rct);  setP(pdt);  setH12(v12); // member data on the parent class
      DE = dE; // member data on the derived class only
    };

    // Destructor
    virtual ~flippedTranslatedLowerAdiabat() {};

    // Energy evaluation routine
    double evaluateEnergy(double r) {
      double V1 = getR()->evaluateEnergy(r);
      double V2 = getP()->evaluateEnergy(r);
      double H12 = getH12();
      double Energy = -0.5*(V1 + V2 - sqrt((V1 - V2)*(V1 - V2) + 4.0*H12*H12)) + DE;
      setEnergy(Energy);
      return Energy;
    };

  };

  // derived class describing the energy difference between two diabats
  class diffBetweenTwoDiabats : public AdiabaticCurve {

  public:
    // Constructor
    diffBetweenTwoDiabats(OneDimensionalFunction* rct, OneDimensionalFunction* pdt) {
      setR(rct);  setP(pdt); // member data on the parent class
    };

    // Destructor
    virtual ~diffBetweenTwoDiabats() {};

    // Energy evaluation routine
    double evaluateEnergy(double r) {
      double V1 = getR()->evaluateEnergy(r);
      double V2 = getP()->evaluateEnergy(r);
      double Energy = V1 - V2;
      setEnergy(Energy);
      return Energy;
    };

  };

  // derived class describing the lower adiabatic curve shifted by DE
  class lowerAdiabatShiftedByE : public AdiabaticCurve {

  private:
    double DE;

  public:
    // Constructor

    lowerAdiabatShiftedByE(OneDimensionalFunction* rct, OneDimensionalFunction* pdt, double v12, double dE) {
      setR(rct);  setP(pdt);  setH12(v12); DE = dE;
    };

    // Destructor
    virtual ~lowerAdiabatShiftedByE() {};

    void setDE(double dE) { DE = dE; };

    // Energy evaluation routine
    double evaluateEnergy(double r) {
      double V1 = getR()->evaluateEnergy(r);
      double V2 = getP()->evaluateEnergy(r);
      double H12 = getH12();
      double Energy = 0.5*(V1 + V2 - sqrt((V1 - V2)*(V1 - V2) + 4.0*H12*H12)) + DE;
      setEnergy(Energy);
      return Energy;
    };

  };

  // derived class describing the lower adiabatic curve shifted by DE
  class upperAdiabatShiftedByE : public AdiabaticCurve {

  private:
    double DE;

  public:
    // Constructor

    upperAdiabatShiftedByE(OneDimensionalFunction* rct, OneDimensionalFunction* pdt, double v12, double dE) {
      setR(rct);  setP(pdt);  setH12(v12);  DE = dE;
    };

    // Destructor
    virtual ~upperAdiabatShiftedByE() {};

    void setDE(double dE) { DE = dE; };

    // Energy evaluation routine
    double evaluateEnergy(double r) {
      double V1 = getR()->evaluateEnergy(r);
      double V2 = getP()->evaluateEnergy(r);
      double H12 = getH12();
      double Energy = 0.5*(V1 + V2 + sqrt((V1 - V2)*(V1 - V2) + 4.0*H12*H12)) + DE;
      setEnergy(Energy);
      return Energy;
    };

  };

  // derived class describing the integrand of the classical action on the lower adiabatic curve
  class lowerAdiabatActionIntegrand : public AdiabaticCurve {

  private:
    double mass, EinAU;

  public:
    // Constructor

    lowerAdiabatActionIntegrand(OneDimensionalFunction* rct, OneDimensionalFunction* pdt, double v12, double m, double E) {
      setR(rct);  setP(pdt);  setH12(v12); mass = m; EinAU = E;
    };

    // Destructor
    virtual ~lowerAdiabatActionIntegrand() {};

    void setmass(double m) { mass = m; };
    void setEinAU(double E) { EinAU = E; };

    // Energy evaluation routine
    double evaluateEnergy(double r) {
      double V1, V2, H12, LowerAdiabatEnergy, integrand, Energy;

      V1 = getR()->evaluateEnergy(r);
      V2 = getP()->evaluateEnergy(r);
      H12 = getH12();
      LowerAdiabatEnergy = 0.5*(V1 + V2 - sqrt((V1 - V2)*(V1 - V2) + 4.0*H12*H12));
      integrand = sqrt(2 * mass*abs(EinAU - LowerAdiabatEnergy));

      Energy = integrand;
      setEnergy(integrand);  // in this case, the 'energy' value is the value of the action integrand at a particular r
      return Energy;
    };

  };

  // derived class describing the integrand of the classical action on the lower adiabatic curve
  class upperAdiabatActionIntegrand : public AdiabaticCurve {

  private:
    double mass, EinAU;

  public:
    // Constructor

    upperAdiabatActionIntegrand(OneDimensionalFunction* rct, OneDimensionalFunction* pdt, double v12, double m, double E) {
      setR(rct);  setP(pdt);  setH12(v12); mass = m; EinAU = E;
    };

    // Destructor
    virtual ~upperAdiabatActionIntegrand() {};

    void setmass(double m) { mass = m; };
    void setEinAU(double E) { EinAU = E; };

    // Energy evaluation routine
    double evaluateEnergy(double r) {
      double V1, V2, H12, upperAdiabatEnergy, integrand, Energy;

      V1 = getR()->evaluateEnergy(r);
      V2 = getP()->evaluateEnergy(r);
      H12 = getH12();
      upperAdiabatEnergy = 0.5*(V1 + V2 + sqrt((V1 - V2)*(V1 - V2) + 4.0*H12*H12));
      integrand = sqrt(2 * mass*abs(EinAU - upperAdiabatEnergy));

      Energy = integrand;
      setEnergy(integrand);  // in this case, the 'energy' value is the value of the action integrand at a particular r
      return Energy;
    };

  };

}//namespace

#endif // GUARD_AdiabaticCurve_h