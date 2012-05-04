//-------------------------------------------------------------------------------------------
//
// DefinedSumOfStates.cpp
//
// Author: Struan Robertson
// Date:   24/Mar/2012
//
// This class implements the defined sum of states model. This class exists to allow users to
// read in their own transition states sum of states, e.g. from FTST calculations or similar.
//
//-------------------------------------------------------------------------------------------
#include <vector>
#include <string>
#include "../System.h"

using namespace std;
using namespace Constants;

namespace mesmer
{
	class DefinedSumOfStates : public MicroRateCalculator
	{
	public:

		// Constructor which registers with the list of MicroRateCalculators in the base class.
		DefinedSumOfStates(const std::string& id) : MicroRateCalculator(id) {}

		  virtual ~DefinedSumOfStates() {}
		  virtual DefinedSumOfStates* Clone() { return new DefinedSumOfStates(*this); }

		  virtual bool calculateMicroRateCoeffs(Reaction* pReac);

		  virtual bool ReadParameters(Reaction* pReac) ;

	private:

		typedef vector<pair<double,double> > WEJ ;

	};

	//************************************************************
	//Global instance, defining its id (usually the only instance) but here with an alternative name
	DefinedSumOfStates theDefinedSumOfStates("DefinedSumOfStates");
	//************************************************************

	bool DefinedSumOfStates::ReadParameters(Reaction* pReact) {

		return true ; 
	}

	//
	// This method calculates the reaction flux. 
	//

	bool DefinedSumOfStates::calculateMicroRateCoeffs(Reaction* pReact)
	{

		PersistPtr ppReac = pReact->get_PersistentPointer();

		PersistPtr pp = ppReac->XmlMoveTo("me:SumOfStates") ;

		if (!pp) {
			throw(std::runtime_error("No sum of states define for reaction"+pReact->getName()))  ;
		}

		const char* txt = pp->XmlReadValue("units", optional);
		string units = txt ? txt : "kJ/mol";

		map<int, WEJ* > SumOfStatesData ;
		while (pp = pp->XmlMoveTo("me:SumOfStatesPoint")) {

			// Read data from XML representation.

			double SumOfStates(0.0) ; 
			txt = pp->XmlRead();
			stringstream s1(txt); s1 >> SumOfStates;
			double energy = pp->XmlReadDouble("energy", false);
			energy = getConvertedEnergy(units, energy);
			int angMomMag = pp->XmlReadInteger("angMomMag", false);

			// Test to see if a new value for the angular momentum magnitude has been found. 

			if (IsNan(angMomMag)) angMomMag = 0 ;
			WEJ* pWEJ = SumOfStatesData[angMomMag] ;

			if (!pWEJ) {
				pWEJ = new WEJ ;
			  SumOfStatesData[angMomMag] = pWEJ ;
			}

			pWEJ->push_back(pair<double,double>(energy, SumOfStates)) ;

		}


		// Delete temporary data representation.
		map<int, WEJ* >::iterator itr =  SumOfStatesData.begin() ;
		for (; itr != SumOfStatesData.end() ; itr++) {
			delete itr->second ;
		}

		return true;
	}

}//namespace
