#include "../DensityOfStates.h"
#include "../Molecule.h"
#include "../gDensityOfStates.h"

using namespace std;
namespace mesmer
{
  class ClassicalRotor : public DensityOfStatesCalculator
  {
  public:

    //Read data from XML. 
    virtual bool ReadParameters(gDensityOfStates* gdos, PersistPtr ppDOSC = NULL);

    // Function to define particular counts of the DOS of a molecule.
    virtual bool countCellDOS(gDensityOfStates* mol, const MesmerEnv& env);

    // Function to calculate contribution to canonical partition function.
    virtual void canPrtnFnCntrb(gDensityOfStates* gdos, double beta, double &PrtnFn, double &IntrlEne, double &varEne);

    // Function to return the number of degrees of freedom associated with this count.
    virtual size_t NoDegOfFreedom(gDensityOfStates* gdos);

    // Constructor which registers with the list of DensityOfStatesCalculators in TopPlugin
    // This class calculates a complete DOS: it is not an extra class. 
		ClassicalRotor(const char* id) : m_id(id), m_Sym(1.0), m_OpticalSym(1.0), m_SpinMultiplicity(1) { Register(); }

    virtual ~ClassicalRotor() {}
    virtual const char* getID() { return m_id; }
    virtual string getName() { return string(m_id); }

    // Included only in a subset of DensityOfStatesCalculators.
    // Otherwise the baseclass function returns false.
    virtual bool includesRotations() { return true; }

    virtual ClassicalRotor* Clone() { return new ClassicalRotor(*this); }

  private:

    const char* m_id;

    double m_Sym;              // Rotational Symmetry Number.
    double m_OpticalSym;       // Transition states may have 2 optical isomers, which
                               // is accounted for by an additionsl symmetry number. 
    size_t m_SpinMultiplicity; // Spin multiplicity.

  };

  //************************************************************
  //Global instance, defining its id
  ClassicalRotor theClassicalRotor("ClassicalRotors");
  //************************************************************

  //Read data from XML and store in this instance.
  bool ClassicalRotor::ReadParameters(gDensityOfStates* gdos, PersistPtr ppDOSC)
  {
    PersistPtr pp = gdos->getHost()->get_PersistentPointer();

    PersistPtr ppPropList = pp->XmlMoveTo("propertyList");
    if (!ppPropList)
      ppPropList = pp; // A propertyList element is not essential.

    // Spin multiplicity (Note spin can be defined as an attribute on a molecule).
    m_SpinMultiplicity = ppPropList->XmlReadPropertyInteger("me:spinMultiplicity", optional);
    if (m_SpinMultiplicity == 0)
      m_SpinMultiplicity = pp->XmlReadInteger("spinMultiplicity");

    // Rotational Symmetry Number
    m_Sym = ppPropList->XmlReadPropertyDouble("me:symmetryNumber");

    // Transition states may have 2 optical isomers, which
    // is accounted for by an additionsl symmetry number.
    if (gdos->getHost()->isMolType("transitionState")) {
      m_OpticalSym = ppPropList->XmlReadPropertyDouble("me:TSOpticalSymmetryNumber", optional);
      if (!IsNan(m_OpticalSym)) {
        if (m_OpticalSym > 2.0) {
          string name(gdos->getHost()->getName());
          cinfo << "Warning: The transition state " << name << " has an optical symmetry number greater than 2." << endl;
        }

        // Adjust Rotational Symmetry Number by Optical Symmetry Number if necessary.
        if (m_OpticalSym > 1.0)
          m_Sym /= m_OpticalSym;
      }
    }

    return true;
  }

  // Provide a function to define particular counts of the DOS of a molecule.
  bool ClassicalRotor::countCellDOS(gDensityOfStates* pDOS, const MesmerEnv& env)
  {
    const size_t MaximumCell = env.MaxCell;
    const double cellSize = env.CellSize;

    vector<double> cellEne;
    getCellEnergies(MaximumCell, cellSize, cellEne);
    vector<double> cellDOS(MaximumCell, 0.0);

    //
    // Initialize density of states array using calculated rotational
    // density of state from inverse Laplace transform of rotors.
    //
    vector<double> rotConst;
    RotationalTop rotorType = pDOS->get_rotConsts(rotConst);
    double qele = m_SpinMultiplicity;
    double cnt = 0.;

    switch (rotorType) {
    case NONLINEAR: //3-D symmetric/asymmetric/spherical top
      cnt = qele * sqrt(4. / (rotConst[0] * rotConst[1] * rotConst[2])) / m_Sym;
      for (size_t i(0); i < MaximumCell; ++i) {
        double Elower = cellEne[i] - 0.5*env.CellSize;
        double Eupper = cellEne[i] + 0.5*env.CellSize;
        Elower *= sqrt(Elower);
        Eupper *= sqrt(Eupper);
        cellDOS[i] = 2.0*cnt*(Eupper - Elower) / 3.0;
      }
      break;
    case LINEAR: //2-D linear
      cnt = qele / (rotConst[0] * m_Sym);
      for (size_t i(0); i < MaximumCell; ++i)
        cellDOS[i] = cnt*env.CellSize;
      break;
    default: // Assume atom.
      cellDOS[0] = qele;
      break;
    }

    // Electronic excited states.
    vector<double> eleExc;
    pDOS->getEleExcitation(eleExc);
    vector<double> tmpCellDOS(cellDOS);
    for (size_t j(0); j < eleExc.size(); ++j) {
      size_t nr = nint(eleExc[j] / cellSize);
      if (nr < MaximumCell) {
        for (size_t i(0); i < MaximumCell - nr; i++) {
          tmpCellDOS[i + nr] = tmpCellDOS[i + nr] + cellDOS[i];
        }
      }
    }

    pDOS->setCellDensityOfStates(tmpCellDOS);

    return true;
  }

  // Calculate contribution to canonical partition function.
  void  ClassicalRotor::canPrtnFnCntrb(gDensityOfStates* gdos, double beta, double &PrtnFn, double &IntrlEne, double &varEne) {

    vector<double> rotConst;
    RotationalTop rotorType = gdos->get_rotConsts(rotConst);

    double qtot(1.0), ene(0.0), var(0.0);
    qtot *= double(m_SpinMultiplicity);

    switch (rotorType) {
    case NONLINEAR://3-D symmetric/asymmetric/spherical top
      qtot *= (sqrt(M_PI / (rotConst[0] * rotConst[1] * rotConst[2]))*(pow(beta, -1.5)) / m_Sym);
      ene = 1.5 / beta;
      var = 1.5 / (beta*beta);
      break;
    case LINEAR://2-D linear
      qtot /= (rotConst[0] * m_Sym*beta);
      ene = 1.0 / beta;
      var = 1.0 / (beta*beta);
      break;
    default:
      break; // Assume atom.
    }

    // Electronic excited states.
    vector<double> eleExc;
    gdos->getEleExcitation(eleExc);
    double qelec(1.0), Eelec(0.0), varEelec(0.0);
    for (size_t j(0); j < eleExc.size(); ++j) {
      qelec += exp(-beta*eleExc[j]);
      Eelec += eleExc[j] * exp(-beta*eleExc[j]);
      varEelec += eleExc[j] * eleExc[j] * exp(-beta*eleExc[j]);
    }
    Eelec /= qelec;
    varEelec = varEelec / qelec - Eelec*Eelec;

    PrtnFn *= qtot*qelec;
    IntrlEne += ene + Eelec;
    varEne += var + varEelec;

    ThermoDynamicEntry(beta, qtot * qelec, ene + Eelec, var + varEelec);

  }

  // Function to return the number of degrees of freedom associated with this count.
  size_t ClassicalRotor::NoDegOfFreedom(gDensityOfStates* gdos) {

    vector<double> rotConst;
    RotationalTop rotorType = gdos->get_rotConsts(rotConst);

    size_t nDOF(0);
    switch (rotorType) {
    case NONLINEAR:
      nDOF = 3;
      break;
    case LINEAR:
      nDOF = 2;
      break;
    default:
      // Assume atom.
      break;
    }

    return nDOF;
  }

}//namespace
