#include "ClassicalRotor.h"

using namespace std;
namespace mesmer
{
	//************************************************************
	//Global instance, defining its id (usually the only instance)
	ClassicalRotor theClassicalRotor("Classical rotors");
	//************************************************************

	// provide a function to define particular counts of the convoluted DOS of two molecules
	bool ClassicalRotor::countDimerCellDOS(ModelledMolecule* p_mol1, ModelledMolecule*  p_mol2, vector<double>& rctsCellEne, vector<double>& rctsCellDOS){
		//----------
		// Differetiating the rotors
		// <three types of rotors: (0) non-rotor (1) 2-D linear, (2) 3-D symmetric top, (3) 3-D asymmetric top>
		// hence there are 9 combinations for convolution:
		// 0-1, 0-2, 0-3, 1-1, 1-2, 1-3, 2-2, 2-3, 3-3
		// where the first three convolutions are trivial
		//

		if (!p_mol1 || !p_mol2){
			cerr << "Input molecules not defined in countrctsCellDOS.";
			return false;
		}

		vector<double> rotCs1; int rotor1Type = p_mol1->get_rotConsts(rotCs1);
		vector<double> rotCs2; int rotor2Type = p_mol2->get_rotConsts(rotCs2);

		int MaximumCell = p_mol2->getEnv().MaxCell;

		p_mol1->getCellEnergies(rctsCellEne) ; // make sure the cell energies are calculated for both molecules.
		p_mol2->getCellEnergies(rctsCellEne) ;

		double I1 = 0., I2 = 0.,                                            //for 2-D linear rotors
			I1X = 0., I2X = 0., I1Y = 0., I2Y = 0., I1Z = 0., I2Z = 0.,  //for 3-D symmetric/asymmetric/spherical top rotors
			s1 = p_mol1->get_Sym(), q1 = static_cast<double>(p_mol1->getSpinMultiplicity()),
			s2 = p_mol2->get_Sym(), q2 = static_cast<double>(p_mol2->getSpinMultiplicity());

		/*All rotational constants are sorted --- the largest one at the head, smallest one at the tail*/


		//---------------------- Assign rotational constants ----------------------
		if      (rotor1Type  <  0){}                                                   // non-rotor
		else if (rotor1Type ==  0){ I1 = rotCs1[0]; }                                  // 2-D linear
		else                      { I1X = rotCs1[0]; I1Y = rotCs1[1]; I1Z = rotCs1[2];}// 3-D symmetric/asymmetric/spherical top
		if      (rotor2Type  <  0){}                                                   // non-rotor
		else if (rotor2Type ==  0){ I2 = rotCs2[0]; }                                  // 2-D linear
		else                      { I2X = rotCs2[0]; I2Y = rotCs2[1]; I2Z = rotCs2[2];}// 3-D symmetric/asymmetric/spherical top

		//------------------------------------------------------------------------------
		//Density of states from ILT of the product of partition functions of two rotors
		int rotorType = rotor1Type + rotor2Type;
		if      (rotorType == -8) return false;              // both not rotors
		double cnt = q1 * q2 / (s1 * s2);

		if      (rotorType <  0 && rotor1Type < rotor2Type){ // only p_mol2 a rotor
			p_mol2->getCellDensityOfStates(rctsCellDOS) ;
			p_mol2->getCellEnergies(rctsCellEne) ;
		}
		else if (rotorType <  0 && rotor1Type > rotor2Type){ // only p_mol1 a rotor
			p_mol1->getCellDensityOfStates(rctsCellDOS) ;
			p_mol1->getCellEnergies(rctsCellEne) ;
		}
		else { 
			if (rotorType == 0){                            // both 2-D linear
				cnt /= (I1 * I2);
				for (int i = 0; i < MaximumCell; ++i) rctsCellDOS[i] = cnt * rctsCellEne[i];
			}
			else if (rotorType == 2 && rotor1Type < rotor2Type){ // 2-D linear vs 3-D symmetric/asymmetric/spherical top
				cnt *= (4./(3.* I1 * sqrt(I2X * I2Y * I2Z)));
				for (int i = 0; i < MaximumCell; ++i) rctsCellDOS[i] = cnt * pow(rctsCellEne[i], 1.5);
			}
			else if (rotorType == 2 && rotor1Type > rotor2Type){ // 3-D symmetric/asymmetric/spherical top vs 2-D linear
				cnt *= (4./(3.* I2 * sqrt(I1X * I1Y * I1Z)));
				for (int i = 0; i < MaximumCell; ++i) rctsCellDOS[i] = cnt * pow(rctsCellEne[i], 1.5);
			}
			else if (rotorType == 4){                            // both 3-D symmetric/asymmetric/spherical top
				cnt *= (M_PI /(2. * sqrt(I1X * I1Y * I1Z * I2X * I2Y * I2Z)));
				for (int i = 0; i < MaximumCell; ++i) rctsCellDOS[i] = cnt * (rctsCellEne[i] * rctsCellEne[i]);
			}

			//-----------------------------------------------------------------------------------------------------
			// convolution of vibrational DOS onto rotational DOS -- loop through all frequencies of both molecules
			vector<double> vfMols; p_mol1->get_VibFreq(vfMols);

			// times the scale factor (the original vibrational frequencies vector still contains the unscaled values)
			for (vector<double>::size_type i = 0; i < vfMols.size(); ++i){
				vfMols[i] *= p_mol1->get_scaleFactor();
			}
			vector<double> vfMol2; p_mol2->get_VibFreq(vfMol2);
			for (vector<double>::size_type i = 0; i < vfMol2.size(); ++i){
				vfMols.push_back(vfMol2[i] * p_mol2->get_scaleFactor());
			}

			Beyer_Swinehart(vfMols, rctsCellDOS);

			//electronic degeneracy
			vector<double> eleExc1, eleExc2;
			p_mol1->getEleExcitation(eleExc1);
			p_mol2->getEleExcitation(eleExc2);
			if (!eleExc1.empty()){
				for (int j = 0; j < static_cast<int>(eleExc1.size()); ++j){
					int iele = static_cast<int>(eleExc1[j]);
					for (int i = (MaximumCell - 1); i >= (iele - 1); --i){
						rctsCellDOS[i] += rctsCellDOS[i - iele +1];
					}
				}
			}
			if (!eleExc2.empty()){
				for (int j = 0; j < static_cast<int>(eleExc2.size()); ++j){
					int iele = static_cast<int>(eleExc2[j]);
					for (int i = (MaximumCell - 1); i >= (iele - 1); --i){
						rctsCellDOS[i] += rctsCellDOS[i - iele + 1];
					}
				}
			}
		}

		return true;
	}

	// provide a function to define particular counts of the convoluted DOS of two molecules
	bool ClassicalRotor::countCellDOS(ModelledMolecule* pMol)
	{
		vector<double> VibFreq ; 
		pMol->get_VibFreq(VibFreq) ;

		// times the scale factor
		for (vector<double>::size_type i = 0; i < VibFreq.size(); ++i){
			VibFreq[i] *= pMol->get_scaleFactor();
		}

		const int MaximumCell = pMol->getEnv().MaxCell;
        vector<double> cellEne(MaximumCell, 0.0) ;
        vector<double> cellDOS(MaximumCell, 0.0) ;

		//
		// Initialize density of states array using calculated rotational
		// density of state.
		//

		//From inverse Laplace transform of rotors
		vector<double> rotConst; int rotorType = pMol->get_rotConsts(rotConst);
		double sym = pMol->get_Sym();
		double qele = pMol->getSpinMultiplicity();
		double cnt = 0.;

		for (int i = 0 ; i < MaximumCell ; ++i ) {
			cellEne[i] = double(i) + 0.5 ;
		}

		switch (rotorType){
	  case 2: //3-D symmetric/asymmetric/spherical top
		  cnt = qele * sqrt(4./(rotConst[0] * rotConst[1] * rotConst[2]))/sym ;
		  for (int i = 0 ; i < MaximumCell ; ++i ) 
			  cellDOS[i] = cnt*sqrt(cellEne[i]) ;
		  break;
	  case 0: //2-D linear
		  cnt = qele / (rotConst[0] * sym);
		  for (int i = 0 ; i < MaximumCell ; ++i ) 
			  cellDOS[i] = cnt ;
		  break;
	  default:
		  cnt = 0.;
		  for (int i = 0 ; i < MaximumCell ; ++i ) 
			  cellDOS[i] = cnt ;
		}

		// Implementation of the Beyer-Swinehart algorithm.
		Beyer_Swinehart(VibFreq, cellDOS);

		//electronic degeneracy
		vector<double> eleExc;
		pMol->getEleExcitation(eleExc);
		if (!eleExc.empty()){
			for (int j = 0; j < eleExc.size() ; ++j){
				int iele = static_cast<int>(eleExc[j]);
				for (int i = (MaximumCell - 1); i >= (iele - 1); --i){
					cellDOS[i] += cellDOS[i - iele + 1];
				}
			}
		}

		pMol->setCellEnergies(cellEne) ;
		pMol->setCellDensityOfStates(cellDOS) ;

		return true;
	}

}//namespace
