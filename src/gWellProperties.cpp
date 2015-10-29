// MolecularComponents.cpp
//
// Author: Chi-Hsiu Liang
//
//-------------------------------------------------------------------------------------------
#include <stdexcept>
#include <numeric>
#include <cmath>
#include <iomanip>
#include "Molecule.h"
#include "System.h"
#include "ParseForPlugin.h"
#include "gWellProperties.h"
#include "gBathProperties.h"
#include "gStructure.h"

using namespace std;
using namespace Constants;
using namespace OpenBabel;

namespace mesmer
{
  //
  // Constructor, destructor and initialization
  //
  gWellProperties::gWellProperties(Molecule* pMol) : MolecularComponent(),
    m_collisionFrequency(0.0),
    m_ncolloptrsize(0),
    m_lowestBarrier(9e23),
    m_numGroupedGrains(0),
    m_pDistributionCalculator(NULL),
    m_grainDist(),
    m_egvec(NULL),
    m_egval()
  {
    ErrorContext c(pMol->getName());
    m_host = pMol;
  }

  gWellProperties::~gWellProperties()
  {
    if (m_grainDist.size()) m_grainDist.clear();
    //delete m_pEnergyTransferModel;
    for (std::map<string, EnergyTransferModel*>::iterator it = m_EnergyTransferModels.begin();
      it != m_EnergyTransferModels.end(); ++it)
      delete it->second;
  }

  bool gWellProperties::initialization(){

    // Determine the method of DOS calculation or use method from defaults.xml.
    PersistPtr pp = m_host->get_PersistentPointer();
    m_pDistributionCalculator
      = ParseForPlugin<DistributionCalculator>(m_host, "me:DistributionCalcMethod", pp);
    if (!m_pDistributionCalculator)
      return false;


    // Specify the energy transfer probability model.
    // The default value is specified in defaults.xml:
    //  <me:energyTransferModel xsi:type="ExponentialDown" default=true;/>
    EnergyTransferModel* pModel = ParseForPlugin<EnergyTransferModel>
      (m_host, "me:energyTransferModel", pp, true); //
    if (!pModel)
      return false;
    pModel->setParent(m_host);
    return true;


  }

  EnergyTransferModel* gWellProperties::addBathGas(const char* pbathGasName, EnergyTransferModel* pModel)
  {
    // Look up the energy transfer model instance for this bath gas
    bool isGeneralBG = (pbathGasName == NULL);
    if (pbathGasName == NULL) // use the general bath gas
      pbathGasName = getHost()->getMoleculeManager()->get_BathGasName().c_str();
    EnergyTransferModel* pETModel = m_EnergyTransferModels[string(pbathGasName)];
    if (pETModel == NULL)
    {
      //Unrecognized bath gas.
      //If it is not the general bath gas make a new model for it and add it to the map.
      MesmerFlags flgs = m_host->getFlags();
      getHost()->getMoleculeManager()->addmol(pbathGasName, "bathGas", m_host->getEnv(), flgs);
      pETModel = isGeneralBG ? pModel : dynamic_cast<EnergyTransferModel*>(pModel->Clone());
      m_EnergyTransferModels[string(pbathGasName)] = pETModel;
    }
    return pETModel;
  }

  const int gWellProperties::get_grnZPE(){
    double grnZpe = (m_host->getDOS().get_zpe() - m_host->getEnv().EMin) / m_host->getEnv().GrainSize; //convert to grain
    if (grnZpe < 0.0)
      cerr << "Grain zero point energy has to be a non-negative value.";

    return int(grnZpe);
  }

  //
  // Initialize the Collision Operator.
  //
  bool gWellProperties::initCollisionOperator(MesmerEnv& env, Molecule *pBathGasMolecule)
  {
    // Calculate the collision frequency.
    m_collisionFrequency = collisionFrequency(env, pBathGasMolecule);

    // Determine the energy of the reservoir grain, first find the lowest barrier 
	// associated with the current well
    // 
    PersistPtr pp = m_host->get_PersistentPointer();
    PersistPtr ppReservoirSize = pp->XmlMoveTo("me:reservoirSize");

    m_numGroupedGrains = 0; // Reset the number of grains grouped into a reservoir grain to zero.

    if (ppReservoirSize){

      // Check the size of the reservoir.
      double tmpvalue = pp->XmlReadDouble("me:reservoirSize");

      const char* unitsTxt = ppReservoirSize->XmlReadValue("units", false);
      string unitsInput("kJ/mol");
      if (unitsTxt){
        unitsInput = unitsTxt;
      }
      else {
        ctest << "No unit for reservoir size has been supplied, use kJ/mol." << endl;
      }

      const double value(getConvertedEnergy(unitsInput, tmpvalue));
      int grainLoc(int(value / double(m_host->getEnv().GrainSize)));
      int lowestBarrier = int(getLowestBarrier() / double(m_host->getEnv().GrainSize));

      if (grainLoc > 0){
        if (grainLoc > lowestBarrier){
          ctest << "The reservoir size provided is too high, corrected according to the lowest barrier height." << endl;
          grainLoc = lowestBarrier;
        }
      }
      else {
        if (abs(grainLoc) > lowestBarrier){
          ctest << "The reservoir size provided is too low, corrected to zero." << endl;
          grainLoc = 0;
        }
        else {
          grainLoc += lowestBarrier;
        }
      }

      vector<double> popDist; // Grained population distribution.
      normalizedGrnBoltzmannDistribution(popDist) ;
      double fracInRsvr(0.0) ;
      for (size_t i(0); i < size_t(grainLoc) ; ++i) {
        fracInRsvr += popDist[i] ;
      }

      m_numGroupedGrains = grainLoc;
      if (m_numGroupedGrains > 1) {
        double reservoirEnergy(m_numGroupedGrains * m_host->getEnv().GrainSize) ;
        ctest << "The reservoir for " << m_host->getName() << " is " << m_numGroupedGrains << " grains," ;
        ctest << "which is " << reservoirEnergy << " cm-1 (or " << reservoirEnergy / getConvertedEnergy("kJ/mol", 1.0) << " kJ/mol) from the well bottom." << endl;
        ctest << "At equilibrium " << 1.0 - fracInRsvr << " of the " << m_host->getName() << " population is in the active states. " << endl ;
      }

    }

    return true;
  }

  //
  // Diagonalize collision operator
  //
  void gWellProperties::diagonalizeCollisionOperator(qdMatrix *egme)
  {
    // Allocate memory.
    m_egval.clear();
    m_egval.resize(m_ncolloptrsize, 0.0);
    if (m_egvec) delete m_egvec;
    *m_egvec = *egme ;

    m_egvec->diagonalize(&m_egval[0]);

    bool reportBasisSetDetails(false);
    if (reportBasisSetDetails) {
      ctest << "\nEigenvectors of: " << m_host->getName() << endl;
      m_egvec->showFinalBits(0, m_host->getFlags().print_TabbedMatrices);

      // The eigenvalues in this series should contain a zero eigenvalue as 
      // energy transfer operators are conservative.
      ctest << "\nEigenvalues of: " << m_host->getName() << "\n{\n";
      for (size_t i(0); i < m_ncolloptrsize; ++i){
        ctest << m_egval[i] << endl;
      }
      ctest << "}\n";
    }
  }

  //
  // Calculate collision frequency.
  //
  double gWellProperties::collisionFrequency(MesmerEnv env, Molecule *pBathGasMolecule)
  {
    //
    // Lennard-Jones Collision frequency. The collision integral is calculated
    // using the formula of Neufeld et al., J.C.P. Vol. 57, Page 1100 (1972).
    // CONCentration is in molec/cm^3.
    //

    double A = 1.16145;
    double B = 0.14874;
    double C = 0.52487;
    double D = 0.77320;
    double E = 2.16178;
    double F = 2.43787;

    double temp = 1.0 / (boltzmann_RCpK*env.beta);

    // Calculate collision parameter averages.
    double bthMass = 0.0;
    bthMass = pBathGasMolecule->getStruc().getMass();

    double bthSigma = 0.0;
    bthSigma = pBathGasMolecule->getBath().getSigma();

    if (!bthSigma)
      cerr << "me:sigma is necessary for " << pBathGasMolecule->getName()
      << ". Correct input file to remove this error." << endl;

    double bthEpsilon = 0.0;
    bthEpsilon = pBathGasMolecule->getBath().getEpsilon();

    if (!bthEpsilon)
      cerr << "me:epsilon is necessary for " << pBathGasMolecule->getName()
      << ". Correct input file to remove this error.";
    double mu = amu * m_host->getStruc().getMass() * bthMass / (m_host->getStruc().getMass() + bthMass);
    double eam = sqrt(m_host->getBath().getEpsilon() * bthEpsilon);
    double sam = (m_host->getBath().getSigma() + bthSigma) * 0.5;
    double tstr = temp / eam;

    // Calculate collision integral.
    double collFrq = A * exp(-log(tstr) * B) + C * exp(-D * tstr) + E * exp(-F * tstr);

    // Calculate molecular collision frequency.
    collFrq *= (M_PI * sam * sam * 1.0e-20 * sqrt(8. * boltzmann_C * temp / (M_PI * mu)));
    // Calculate overall collision frequency.
    collFrq *= (env.conc * 1.0e6);

    return collFrq;
  }

  //
  // Calculate a reaction matrix element.
  //
  qd_real gWellProperties::matrixElement(int eigveci, int eigvecj, std::vector<double> &k) const
  {
    // Calculate matrix element starting with the higher energy
    // elements first in order to preserve precision as much as possible.
    qd_real sum = 0.0;
    for (int i = m_ncolloptrsize - 1; i >= 0; --i){
      sum += qd_real(k[i]) * ((*m_egvec)[i][eigveci] * (*m_egvec)[i][eigvecj]);
    }
    return sum;
  }

  //
  // Accessor a collision operator eigenvector.
  //
  void gWellProperties::eigenVector(int eigveci, std::vector<double> &evec) const
  {
    // evec.clear() ;
    for (size_t i(0); i < m_ncolloptrsize; ++i){
      evec[i] = to_double((*m_egvec)[i][eigveci]);
    }
  }

  //
  // Copy eigenvalues to diagonal elements of the reaction operator
  // matrix in the contracted basis representation.
  //
  void gWellProperties::copyCollisionOperatorEigenValues(qdMatrix *CollOptr,
    const size_t locate,
    const double Omega) const
  {
    // Check that the contracted basis method has been specifed.

    if (!m_host->getEnv().useBasisSetMethod)
      throw (std::runtime_error("Error: Contracted basis representation not requested."));

    // Find size of system matrix.

    size_t smsize = CollOptr->size();
    size_t nbasis = get_nbasis();

    // Check there is enough space in system matrix.

    if (locate + nbasis > smsize)
      throw (std::runtime_error("Error in the size of the reaction operator matrix in contracted basis representation."));

    // Copy collision operator eigenvalues to the diagonal elements indicated
    // by "locate" and multiply by the reduced collision frequencey.

    for (size_t i(0); i < nbasis; ++i) {
      int ii(locate + i);
      (*CollOptr)[ii][ii] = Omega * m_egval[m_ncolloptrsize - i - 1];
    }
  }

  //
  // Get normalized grain distribution.
  //
  void gWellProperties::normalizedInitialDistribution(vector<double> &grainFrac)
  {

    m_pDistributionCalculator->calculateDistribution(m_host, m_grainDist);

    double prtfn(m_grainDist[0]);
    grainFrac.push_back(prtfn);
    for (size_t i = 1; i < m_ncolloptrsize + reservoirShift(); ++i){
      prtfn += m_grainDist[i];
      if (i < m_numGroupedGrains){
        grainFrac[0] += m_grainDist[i];
      }
      else {
        grainFrac.push_back(m_grainDist[i]);
      }
    }

    for (size_t i(0); i < grainFrac.size(); ++i) {
      grainFrac[i] /= prtfn;
    }

    if (m_host->getFlags().grainBoltzmannEnabled){
      ctest << "\nGrain fraction:\n{\n";
      for (size_t i(0); i < grainFrac.size(); ++i){
        ctest << grainFrac[i] << endl;
      }
      ctest << "}\n";
    }
  }

  //
  // Get normalized grain distribution.
  //
  void gWellProperties::normalizedGrnBoltzmannDistribution(vector<double> &grainFrac)
  {
    vector<double> gDOS;
    m_host->getDOS().getGrainDensityOfStates(gDOS);

    vector<double> gEne;
    m_host->getDOS().getGrainEnergies(gEne);

    vector<double> tempGrnFrac;
    // Calculate unnormalized Boltzmann dist.
    // Note the extra 10.0 is to prevent underflow, it is removed during normalization.
    double prtfn = exp(log(gDOS[0]) - m_host->getEnv().beta * gEne[0] + 10.0);
    tempGrnFrac.push_back(prtfn);
    for (size_t i = 1; i < m_ncolloptrsize + reservoirShift(); ++i) {
      const double thisPartition = exp(log(gDOS[i]) - m_host->getEnv().beta * gEne[i] + 10.0);
      prtfn += thisPartition;
      if (i < m_numGroupedGrains){
        tempGrnFrac[0] += thisPartition;
      }
      else {
        tempGrnFrac.push_back(thisPartition);
      }
    }

    for (size_t i(0); i < tempGrnFrac.size(); ++i){
      tempGrnFrac[i] /= prtfn;
    }

    grainFrac = tempGrnFrac;

  }

  void gWellProperties::normalizedCellBoltzmannDistribution(vector<double> &CellFrac)
  {
    vector<double> gDOS;
    m_host->getDOS().getCellDensityOfStates(gDOS);

    vector<double> gEne;
    getCellEnergies(gDOS.size(), m_host->getEnv().CellSize, gEne);

    vector<double> tempCellFrac(gDOS.size(), 0.0);
    double prtfn(0.0);
    // Calculate unnormalized Boltzmann dist.
    // Note the extra 10.0 is to prevent underflow, it is removed during normalization.
    for (size_t i(0); i < tempCellFrac.size(); ++i) {
      const double thisPartition = (gDOS[i] > 0.0) ? exp(log(gDOS[i]) - m_host->getEnv().beta * gEne[i] + 10.0) : 0.0 ;
      prtfn += thisPartition;
      tempCellFrac[i] = thisPartition ;
    }

    for (size_t i(0); i < tempCellFrac.size(); ++i){
      tempCellFrac[i] /= prtfn;
    }

    CellFrac = tempCellFrac;
  }

  //
  // Accessor for number of basis functions to be used in contracted basis set method.
  //
  size_t gWellProperties::get_nbasis() const { return m_host->getEnv().nBasisSet; }

}//namespace
