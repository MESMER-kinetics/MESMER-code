//-------------------------------------------------------------------------------------------
//
// SimpleBimolecularSink.cpp
//
// Authors: David Glowacki
// Date:    Dec 2011
//
// calculates a sink reaction as an irreversible association process A+B, where one of A+B is
// a modelled molecule; this routine presently assumes that the association rate 
// coefficient is independent of the internal energy of A+B (hence, "Simple")
//
//-------------------------------------------------------------------------------------------
#include <vector>
#include <string>
#include "../System.h"

using namespace std;
using namespace Constants;
namespace mesmer
{
	class SimpleBimolecularSink : public MicroRateCalculator
	{
	public:

		///Constructor which registers with the list of MicroRateCalculators in the base class
		SimpleBimolecularSink(const char* id) : m_id(id), m_AssociationRateCoeff(0.0)
    { Register(); }

    virtual const char* getID()  { return m_id; }

		virtual ~SimpleBimolecularSink() {}

		virtual SimpleBimolecularSink* Clone() { return new SimpleBimolecularSink(*this); }

		//@virtual bool ReadParameters(Reaction* pReac);
    virtual bool ParseData(PersistPtr pp);

		virtual double get_ThresholdEnergy(Reaction* pReac) ;

		virtual bool calculateMicroCnlFlux(Reaction* pReact);

	private:
		Rdouble m_AssociationRateCoeff ; // association rate coefficient

  private:
    const char* m_id;

	};

	//************************************************************
	//Global instance, defining its id (usually the only instance)
	SimpleBimolecularSink theSimpleBimolecularSink("SimpleBimolecularSink");
	//************************************************************

	bool SimpleBimolecularSink::calculateMicroCnlFlux(Reaction* pReact)
	{
		cinfo << "calculating cell flux for reaction " << pReact->getName() << " which is a SimpleBimolecularSink method"<<endl;

		Molecule* pRct = pReact->get_reactant() ;
		if(!pRct){
			cerr << "can't find reactant for reaction " << pReact->getName() << " for SimpleBimolecularSink" << endl;
			return false;
		}
		// Allocate some work space for density of states.
		vector<double> RctCellDOS; // Transistion state density of states.
		if(!pRct->getDOS().getCellDensityOfStates(RctCellDOS))
			return false; // Extract densities of states from molecules.

		// get MaxCell from MesmerEnv structure via Reaction class
		const int MaximumCell = pReact->getEnv().MaxCell;

		// Allocate space to hold transition state flux and initialize elements to zero.
		vector<double>& rxnFlux = pReact->get_CellFlux();
		rxnFlux.clear();
		rxnFlux.resize(MaximumCell, 0.0);

		// As there's presently no TS for the bimolecular sink reaction, the microcanonical flux is calculated simply as k(E)*p(E)
		// where k(E) is energy independent bimolecular rate coefficient multiplied by the excess species concentration
		// and p(E) is the DOS for the modelled molecule
		double excessConc = pReact->get_concExcessReactant();
		for (int i = 0 ; i < MaximumCell ; ++i) {
			rxnFlux[i] = m_AssociationRateCoeff * excessConc * RctCellDOS[i] ;
		}

		// the flux bottom energy is equal to the ZPE of the transition state
		pReact->setCellFluxBottom(pReact->get_relative_rctZPE() + pReact->get_ThresholdEnergy());

		return true;
	}

	double SimpleBimolecularSink::get_ThresholdEnergy(Reaction* pReact)
	{
		//	the simple bimolecular sink presently assumes a threshold energy of zero -
		//		although this will no doubt change with more sophisticated treatments
		return 0.0;
	}

	// for the simple bimolecular sink, all we do is read the BimolecularLossRateCoefficient
	// the excess species concentration lives on the Reaction class

	//@bool SimpleBimolecularSink::ReadParameters(Reaction* pReact)
  bool SimpleBimolecularSink::ParseData(PersistPtr pp)
	{
    Reaction* pReact = m_parent; //use old var name
		cinfo << "Reading parameters for Reaction " << pReact->getName() << " which uses the SimpleBimolecularSink method"<<endl;

		PersistPtr ppReac = pReact->get_PersistentPointer();

		PersistPtr ppRateCoefficient;
		const char* pRateCoefficientTxt=NULL;
		bool rangeSet(false) ;

		ppRateCoefficient = pp->XmlMoveTo("me:BimolecularLossRateCoefficient") ;
		pRateCoefficientTxt = pp->XmlReadValue("me:BimolecularLossRateCoefficient");

		if (pRateCoefficientTxt){
			double value = 0.0;
			stringstream s2(pRateCoefficientTxt); 
			s2 >> value ;

			if(value<0.0){
				cerr << "association rate coefficient should not be negative when used with Simple Bimolecular Sink" << endl;
				return false;
			}
			ReadRdoubleRange(string(pReact->getName()+":BimolecularLossRateCoefficient"), ppRateCoefficient, m_AssociationRateCoeff, rangeSet) ;
			m_AssociationRateCoeff = value ;
			if (rangeSet) {
				double valueL, valueU, stepsize ;
				m_AssociationRateCoeff.get_range(valueL,valueU,stepsize) ;
				if(valueL<0.0){
					cerr << "Lower bound of association rate coefficient should not be negative when used with Simple Bimolecular Sink.";
					return false;
				}
			}
		}
		else{
			cerr << "Requesting Simple Bimolecular Sink without providing the bimolecular loss rate coefficient for reaction "
				<< this->getID() << ". Please correct input file.";
			return false;
		}

		return true;
	}

}//namespace
