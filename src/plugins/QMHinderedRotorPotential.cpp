#include "../Molecule.h"
#include "../MesmerMath.h"
#include "../OneDimensionalFGH.h"
#include "QMHinderedRotorPotential.h"

using namespace std;
namespace mesmer
{
  //-------------------------------------------------------------
  //Global instance, defining the id
  QMHinderedRotorPotential theQMHinderedRotorPotential("QMHinderedRotorPotential");
  //-------------------------------------------------------------

  using OpenBabel::vector3;
  //Read data from XML and store in this instance.
  bool QMHinderedRotorPotential::ReadParameters(Molecule* pMol, PersistPtr ppDOSC)
  {
    gStructure& gs = pMol->getStruc();
    if(!gs.ReadStructure())
    {
      cerr << "A complete set of atom coordinates are required for hindered rotor calculations" <<endl;
      return false;
    }

    stringstream freqStrn(ppDOSC->XmlReadValue("me:vibFreq", optional));
    freqStrn >> m_vibFreq;
    // Alvyn's comment: One mode of hindered rotation is represented by a single conjugated vibrational mode in QM packages.
    // The vibration in the vicinity obeys normal mode analysis, meaning it includes the coupled motion of all atoms.
    // The normal mode motion conserves moments of inertia (Im). Where equally, a rigid rotation of the bond with other coordinates
    // fixed also conserves Im. However, a fully relaxed calculation will introduce noise to the conservation of Im.

    const char* bondID = ppDOSC->XmlReadValue("bondRef", optional);
    if(bondID)
    {
      pair<string,string> bondats = gs.GetAtomsOfBond(bondID);
      if(bondats.first.empty())
      {
        cerr << "Unknown bond reference " << bondID << endl;
        return false;
      }
      m_bondID = bondID;
      cinfo << "Hindered rotor " << m_bondID << endl;

      vector3 coords1 = gs.GetAtomCoords(bondats.first);
      vector3 coords2 = gs.GetAtomCoords(bondats.second);

      //calc Moment of inertia about bond axis of atoms on one side of bond...
      vector<string> atomset;
      gs.GetAttachedAtoms(atomset, bondats.first, bondats.second);
      double mm1 = gs.CalcMomentAboutAxis(atomset, coords1, coords2);

      //...and the other side of the bond
      atomset.clear();
      gs.GetAttachedAtoms(atomset, bondats.second, bondats.first);
      double mm2 = gs.CalcMomentAboutAxis(atomset, coords1, coords2);

      /*
       Is the reduced moment of inertia need about the bond axis or, separately for the set of
       atoms on each side of the bond, about a parallel axis through their centre of mass?
       See:
       http://www.ccl.net/chemistry/resources/messages/2001/03/21.005-dir/index.html
       http://www.ccl.net/chemistry/resources/messages/2001/03/31.002-dir/index.html
       The bond axis is used here.
       */
      m_reducedMomentInertia = mm1 * mm2 / ( mm1 + mm2 );//units a.u.*Angstrom*Angstrom
    }

    // Section to get the cyclic pes data.
    // The cyclic data shall be input as a series of double precision number as energies, where the energies will be converted into cm-1.
    // The cyclic data are reference points for Fourier expansion.
    /*
     \documentclass[10pt,a4paper]{article}
     \usepackage[latin1]{inputenc}
     \usepackage{amsmath}
     \usepackage{amsfonts}
     \usepackage{amssymb}
     \begin{document}
     Let the function be expanded as

     \begin{equation}
     y(\theta)=a_0+a_1 \sin(\theta)+a_2 \sin(2\theta)+\ldots+a_m \sin(m\theta)+b_1 \cos(\theta)+b_2 \cos(2\theta)+\ldots+b_m \cos(m\theta)
     \end{equation}

     where m is a positive integer (the larger the better).
     The Fourier coefficients are given by

     \begin{eqnarray}
     a_0 = \frac{\sum_{i=1}^{n}y_i}{n} \\
     a_k = \frac{2}{n}\sum_{i=1}^{n} y_i \sin(k\theta_i)\\
     b_k = \frac{2}{n}\sum_{i=1}^{n} y_i \cos(k\theta_i)
     \end{eqnarray}

     where $(y_1, \theta_1), (y_2, \theta_2),\ldots, (y_n, \theta_n)$ are ordered pairs in $0 \leq \theta \leq \pi$
     \end{document}
     */

    // initialize the values for Fourier expansion
    ak.clear();
    bk.clear();
    a0 = 0.0;

    // Lines to read in the classical PES from the hindered rotor. The energies shall cover a full circle of the rotation, evenly
    // separated. More points would be required if the PES is complicated.
    vector<double> pesEnes;
    PersistPtr ppPes = ppDOSC->XmlMoveTo("me:hinderedpes");
    const char* p = ppPes->XmlReadValue("units", optional);
    string units = p ? p : "kJ/mol";
    stringstream ss(ppPes->XmlReadValue("array", false));
    double val;
    while(ss >> val)
      pesEnes.push_back(getConvertedEnergy(units, val));

    // might need some code to ground the energies to a minimum values? (Does it matter?)
    // Routine to call function in order to convert the PES into ak, bk and a0.

    // Ground the energies (make the lowest energy 0.0)
    double minimalE = pesEnes[0];
    for (size_t i(1); i < pesEnes.size(); ++i) if (minimalE > pesEnes[i]) minimalE = pesEnes[i];
    for (size_t i(0); i < pesEnes.size(); ++i) pesEnes[i] -= minimalE;
    //

    size_t expansion = 8;
    convertToFourierCoefficients(expansion, ak, bk, a0, pesEnes);

    return true;
  }

  //Adjusts DOS for the specified vibrations being treated as hindered rotors
  bool QMHinderedRotorPotential::countCellDOS(gDensityOfStates* pDOS, int MaximumCell)
  {
    cinfo << "Hindered rotor " << m_bondID << endl;

    vector<double> cellDOS;
    if(!pDOS->getCellDensityOfStates(cellDOS, 0, false)) // retrieve the DOS vector without recalculating
      return false;

    // Remove contribution from vibration and add contribution from hindered rotor

    // 1. Calculate the required energy levels (eigenvalues) from the hindered rotation.
    // How many cm-1 are required ? Estimate the number of states thru the multiple of the first hindered PES vibrational energy
    // Find maximum quantum No. for rotor.

    double bint = conMntInt2RotCnt/m_reducedMomentInertia ;
    double root = sqrt(double(MaximumCell)/bint) ;
    int kmax    = int(root + 1.0) ;

    m_numberGridPoint = 2*kmax +1 ;
    vector<double> eigenvalues;
    dMatrix wavefunctions(m_numberGridPoint, m_numberGridPoint);
    oneDimensionalFourierGridHamiltonian(m_reducedMomentInertia, &getEnergyFromFourierCoefficients,
                                         eigenvalues, wavefunctions, m_numberGridPoint, ak, bk, a0);

    vector<double> hinderedLevels(MaximumCell, 0.0);
    for(size_t i(0); i < eigenvalues.size(); ++i){
      int eneLevel = int(abs(eigenvalues[i]));
      if (eneLevel < int(MaximumCell)){
        hinderedLevels[eneLevel] += 1.0;
      }
    }

    ////----------  The following adjustments of ZPE need to be reconsidered. The energy has to be stored as classical energy in order to be clear.
    // 2. Remove the contribution of the harmonic vibration to zero point energy.
    // energy to be removed from the ZPE = m_vibFreq / 2.0
    // 3. Add the contribution of ground state hindered vibration to zero point energy.
    // need to think about how...
    ////----------

    // 4. De-convolute the harmonic vibrational energy levels from the DOS.
    //reverseBeyer_Swinehart(m_vibFreq, cellDOS);

    // 5. Convolute the hindered rotor energy levels to the DOS. The first vibrational state is ignored and treated as harmonic. It's energy
    //    is divided by two and counted as the contribution to the Zero point energy.
    vector<double> tempCellDOS;
    FastLaplaceConvolution(cellDOS, hinderedLevels, tempCellDOS);

    // 6. Set density of states
    pDOS->setCellDensityOfStates(tempCellDOS);

    return true;
  }

}//namespace
