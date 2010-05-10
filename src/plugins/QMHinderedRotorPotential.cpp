#include "Molecule.h"
#include "../OneDimensionalFGH.h"
#include "QMHinderedRotorPotential.h"

using namespace std;
namespace mesmer
{
  //-------------------------------------------------------------
  //Global instance, defining its id (usually the only instance)
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

    const char* bondID = ppDOSC->XmlReadValue("bondRef",optional);
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
    
    // Lines to read in the PES from the hindered rotor. The energies shall cover a full circle of the rotation, evenly
    // separated. More points would be required if the PES is complicated.
    vector<double> pesEnes;
    //ppDOSC->

    // Routine to call function in order to convert the PES into ak, bk and a0. 
    //convertToFourierCoefficients(ak, bk, a0, pesEnes);
    
    

    return true;
  }

  //Adjusts DOS for the specified vibrations being treated as hindered rotors
  bool QMHinderedRotorPotential::countCellDOS(gDensityOfStates* pDOS, int MaximumCell)
  {
    //bond is needed, and the PES is recorded in the bond as an array.
    //const char* bondID = ppDOSC->XmlReadValue("bondRef",optional);
    const char* bondID = "";    // temporary
    if(bondID)
    {
      PersistPtr ppMol = pDOS->getHost()->get_PersistentPointer();
      PersistPtr ppBond = ppMol->XmlMoveTo("bondArray");
      while(ppBond=ppBond->XmlMoveTo("bond"))
      {
        const char* id = ppBond->XmlReadValue("id");
        if(id && !strcmp(id, bondID))
          break;
      }
      // The PES information about the dihedral torsion of the bond should be recorded in here.
      if(ppBond){
        //ppBond->
      }
      else{
        cerr << "Cannot find the specified bond.";
      }
    }
    cinfo << "Hindered rotor " << bondID << endl;

    //double barrier  = ppDOSC->XmlReadDouble("me:barrierZPE");
    double barrier  = 0.0; // temporary
    //PersistPtr pp = ppDOSC->XmlMoveTo("me:barrierZPE");
    PersistPtr pp = NULL; // temporary
    //const char* p = pp->XmlReadValue("units", optional);
    const char* p = ""; //temporary
    string units = p ? p : "kJ/mol";
    barrier = getConvertedEnergy(units, barrier);

    //...
    vector<double> cellDOS;
    if(!pDOS->getCellDensityOfStates(cellDOS, 0, false)) // retrieve the DOS vector without recalculating
      return false;

    //***TODO remove contribution from vibration and add contribution from hindered rotor

    pDOS->setCellDensityOfStates(cellDOS) ;

    return true;
  }

}//namespace
