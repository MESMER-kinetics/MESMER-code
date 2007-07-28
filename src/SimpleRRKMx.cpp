#include <vector>
#include "System.h"
#include "Molecule.h"
#include "MoleculeManager.h"
#include "Persistence.h"
#include "Reaction.h"
#include "Constants.h"
#include "microrate.h"

using namespace std;
using namespace Constants;
namespace mesmer
{
class SimpleRRKM : public MicroRateCalculator
{
public:

  ///Constructor which registers with the list of MicroRateCalculators in the base class
  SimpleRRKM(const string& id) : MicroRateCalculator(id){}

  virtual ~SimpleRRKM() {}

  virtual bool calculateMicroRateCoeffs(Reaction* pReact, vector<double> &kfmc);
};

//************************************************************
//Global instance, defining its id (usually the only instance)
SimpleRRKM theSimpleRRKM("Simple RRKM");
//************************************************************

bool SimpleRRKM::calculateMicroRateCoeffs(Reaction* pReact, vector<double> &kfmc)
{
  vector<CollidingMolecule *> unimolecularspecies;
  pReact->get_unimolecularspecies(unimolecularspecies);
  CollidingMolecule * pReactant = unimolecularspecies[0];

  TransitionState* pTS = pReact->get_TransitionState();
  if(!pTS)
  {
    cerr << "No transition state in Simple RRKM for reaction " << pReact->getName() << endl;
    return false;
  }

  // Allocate space to hold Micro-canonical rate coefficients.
  kfmc.resize(pSys->MAXCell());

  // Initialize microcanoincal rate coefficients.

  int i, j ;
  for (i = 0 ; i < pSys->MAXCell() ; ++i ) {
    kfmc[i] = 0.0 ;
  }

  // Allocate some work space for density of states.

  vector<double> ddosTS(pSys->MAXCell(),0.0) ; // Transistion state density of states.
  vector<double> ddos(pSys->MAXCell(),0.0) ; // Density of states of equilibrium molecule.

  // Extract densities of states from molecules.

  pReactant->cellDensityOfStates(&ddos[0]) ;
  pTS->cellDensityOfStates(&ddosTS[0]) ;

  double SumOfStates  = 0.0 ;
  int thresholdEnergy = int((pTS->get_zpe() - pReactant->get_zpe()) * KCMLTOPCM) ;
  for (i = thresholdEnergy, j = 0 ; i < pSys->MAXCell() ; ++i, ++j ) {

    // Integrate transition state density of states.

    SumOfStates += ddosTS[j] ;

    // Calculate microcanonical rate coefficients using RRKM expression.

    kfmc[i] = SumOfStates / (plancksConst*ddos[i]) ;
  }
  return true;
}

}//namespace