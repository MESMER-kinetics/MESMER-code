#include <vector>
#include <float.h>
#include "System.h"
#include "Reaction.h"
#include "Constants.h"
#include "MicroRate.h"

#if defined(WIN32)
#define IsNan _isnan
#endif

using namespace std;
using namespace Constants;
namespace mesmer
{
class SimpleILT : public MicroRateCalculator
{
public:

  ///Constructor which registers with the list of MicroRateCalculators in the base class
  SimpleILT(const string& id) : MicroRateCalculator(id){}

  virtual ~SimpleILT() {}

  virtual bool calculateMicroRateCoeffs(Reaction* pReact, vector<double> &kfmc);
};

//************************************************************
//Global instance, defining its id (usually the only instance)
SimpleILT theSimpleILT("Simple ILT");
//************************************************************

bool SimpleILT::calculateMicroRateCoeffs(Reaction* pReact, vector<double> &kfmc)
{
  vector<CollidingMolecule *> unimolecularspecies;
  pReact->get_unimolecularspecies(unimolecularspecies);
  CollidingMolecule * pReactant = unimolecularspecies[0];
  if(IsNan(pReact->get_ActivationEnergy()))
  {
    cerr << "To use SimpleILT for reaction " << pReact->getName()
      << " the Activation Energy needs to be set." << endl;
      return false;
  }

  // Allocate space to hold Micro-canonical rate coefficients.
  kfmc.resize(pSys->MAXCell());

  // Initialize microcanoincal rate coefficients.

  int i ;
  for (i = 0 ; i < pSys->MAXCell() ; ++i ) {
      kfmc[i] = 0.0 ;
  }

  // Allocate some work space for density of states.

  vector<double> ddos(pSys->MAXCell(),0.0) ; // Density of states of equilibrim molecule.

  // Extract densities of states from molecules.

  pReactant->cellDensityOfStates(&ddos[0]) ;

  // Conversion of EINF from kcal.mol^-1 to cm^-1

  int nEinf = int(pReact->get_ActivationEnergy()*KCMLTOPCM) ;

  // Calculate microcanonical rate coefficients using simple ILT expression.

  	for (i = nEinf ; i < pSys->MAXCell() ; ++i ) {
            kfmc[i] = pReact->get_PreExp()*ddos[i-nEinf] / ddos[i] ;
    }

  return true;
}

}//namespace