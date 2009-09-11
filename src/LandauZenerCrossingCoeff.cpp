#include "LandauZenerCrossingCoeff.h"

using namespace Constants;

namespace mesmer
{
	//************************************************************
	//Global instance, defining its id (usually the only instance)
	LandauZenerCrossingCoeff theLandauZenerCrossingCoeff("LZ");
	//************************************************************


	bool LandauZenerCrossingCoeff::calculateCellCrossingCoeffs(Reaction* pReact, vector<double>& CrossingProbability){

		Molecule * p_TransitionState = pReact->get_TransitionState();

		// read input data for Landau-Zener crossing 

		PersistPtr pp = p_TransitionState->get_PersistentPointer();
		PersistPtr ppPropList = pp->XmlMoveTo("propertyList");
		if(!ppPropList)
			ppPropList=pp; //Be forgiving; we can get by without a propertyList element

		cinfo << endl << "Spin Forbidden Crossing Data for reaction " << pReact->getName() << endl;

		double SOCelement(0.0), GradDiffMagnitude(0), ReducedMass(0);

		if(!ReadDoubleAndUnits(SOCelement, ppPropList, "me:RMS_SOC_element", "cm-1")){
			return false;
		}

		if(!ReadDoubleAndUnits(GradDiffMagnitude, ppPropList, "me:GradientDifferenceMagnitude", "a.u./Bohr")){
			return false;
		}

		if(!ReadDoubleAndUnits(ReducedMass, ppPropList, "me:GradientReducedMass", "a.m.u.")){
			return false;
		}

		// end read input data section

		std::vector<Molecule *> unimolecularspecies;
		pReact->get_unimolecularspecies(unimolecularspecies);
		Molecule * pReactant = unimolecularspecies[0];
		double ZPE_corr_barrier_height;

		// get threshold energy:
		ZPE_corr_barrier_height  = pReact->get_ThresholdEnergy();

		const double SOCelementAU = SOCelement * (1/Hartree_in_RC);
		const double ReducedMassAU = ReducedMass * 1.822888e+3;

		//get properties of vectors in which to include Crossing coefficients
		const int MaximumCell = pReactant->getEnv().MaxCell;
		CrossingProbability.clear();
		CrossingProbability.resize(MaximumCell);

		//set transmission coefficients to 0 below the ZPE corrected barrier height;
		//above the barrier, the Landau Zener transmission coefficients are calculated 
		//as described by Harvey & Aschi, Faraday Discuss, 2003 (124) 129-143

		for(int i = 0; i < MaximumCell; ++i){
			double E = double(i) - ZPE_corr_barrier_height;
			if (E <= 0){
				CrossingProbability[i] = 0.0;
			}
			else
			{
				double E_AU = E / Hartree_in_RC;
				double trans_probability = exp(-2.0e+0 * M_PI * pow(SOCelementAU,2.0e+0) / GradDiffMagnitude * pow(ReducedMassAU / (2.e+0 * E_AU),0.5e+0));
				CrossingProbability[i] = (1.0e+0 + trans_probability)*(1.0e+0 - trans_probability);
				// following if statement to avoid nan
				if(IsNan(CrossingProbability[i])) CrossingProbability[i] = 0.0;
			}
		}

		if (pReact->getFlags().CrossingCoeffEnabled){
			ctest << "\nCrossing coefficients for: " << pReact->getName() << endl;
			for(int i = 0; i < MaximumCell; ++i){
				ctest << CrossingProbability[i] << endl;
			}
			ctest << "}\n";
		}

		return true;
	}

	bool LandauZenerCrossingCoeff::ReadDoubleAndUnits(double& element, PersistPtr pp, const std::string identifier, const std::string units){

		element=pp->XmlReadPropertyDouble(identifier,true);
		string unitsTxt;

		if(element){
			unitsTxt = pp->XmlReadPropertyAttribute(identifier, "units");
			if (unitsTxt!=units){
				cinfo << "MESMER could not read units for " << identifier << "; assuming " << units << "." << endl;
			}
			cinfo << identifier << " = " << element << " " << units << endl;
		}  
		else{
			cerr << "Spin forbidden crossing: failed to read " << identifier << " (" << units << ")." << endl; 
			return false;
		}
		return true;
	}

}//namespace

