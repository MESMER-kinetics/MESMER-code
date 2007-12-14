#include "Reaction.h"
//
// Test the forward microcanonical rate coefficients.
//
using namespace std;
using namespace Constants ;

namespace mesmer
{

  bool MicroRateCalculator::testMicroRateCoeffs(Reaction*         pReact,
                                                vector<double> &cellKfmc,
                                                PersistPtr        ppbase,
                                                const MesmerEnv    &m_Env) const
  {
    vector<ModelledMolecule *> unimolecularspecies;
    pReact->get_unimolecularspecies(unimolecularspecies);
    ModelledMolecule * pReactant = unimolecularspecies[0];

    {
      stringstream errorMsg;
      errorMsg << "Test of microcanonical rate coefficients";
      meErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obInfo);
    }

    string comment("Microcanonical rate coefficients");
    PersistPtr ppList = ppbase->XmlWriteMainElement("me:microRateList", comment );

    // Allocate some work space for density of states.

    vector<double> CellEne(m_Env.MaxCell,0.0) ;
    vector<double> cellDOS(m_Env.MaxCell,0.0) ;

    pReactant->getCellEnergies(CellEne) ;
    pReactant->getCellDensityOfStates(cellDOS) ;

    for(int i = 0 ; i < 29 ; ++i)
    {
      double Temperature = double(i+2)*100.0 ;
      double beta = 1.0/(boltzmann_RCpK*Temperature) ;

      double sm1 = 0.0, sm2 = 0.0, tmp = 0.0;

      for ( int i = 0 ; i < m_Env.MaxCell ; ++i ) {
        tmp  = cellDOS[i] * exp(-beta * CellEne[i]) ;
        sm1 += cellKfmc[i] * tmp ;
        sm2 +=           tmp ;
      }
      sm1 /= sm2 ;
      formatFloat(ctest, Temperature, 6,  7) ;
      formatFloat(ctest, sm1,         6, 15) ;
      ctest << endl ;

      //Add to XML document
      PersistPtr ppItem = ppList->XmlWriteElement("me:microRate");
      ppItem->XmlWriteValueElement("me:T",   Temperature, 6) ;
      ppItem->XmlWriteValueElement("me:val", sm1,         6) ;
    }
    return true;
  }
}//namespace
