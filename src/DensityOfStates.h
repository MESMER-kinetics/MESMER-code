#ifndef GUARD_DensityOfStates_h
#define GUARD_DensityOfStates_h

#include <map>
#include <string>
#include <vector>
#include "XMLPersist.h"
#include "plugin.h"

namespace mesmer
{
  class Molecule;
  class gDensityOfStates;

  // Abstract base class for cell Density Of States (DOS) calculators. The
  // derived concrete classes are plugin classes, so that new classes can be 
  // added without changing any of the existing code. The constructor of a global 
  // instance registers the class with the base class. Subsequently, supplying  
  // the id (a string) to the Find function returns a pointer to a new instance.

  class DensityOfStatesCalculator : public TopPlugin
  {
  public:
    DensityOfStatesCalculator(){}
    virtual ~DensityOfStatesCalculator(){}
    virtual const char* getTypeID()  {return typeID();}
    virtual const char* typeDescription() { return
      "Normally it necessary to specify a density of states calculator that\n"
      "includes the rotations, such as ClassicalRotors, QMRotors or DefinedStates.\n"
      "In the XML datafile this is done like\n" 
      "<me:DOSCMethod>QMRotors/> or <me:DOSCMethod name=\"QMRotors\"/>\n\n"
      "If no method is specified with <me:DOSCMethod>, the value from defaults.xml\n"
      "(currently ClassicalRotors) is used. The BeyerSwinehart method for vibrations"
      " is normally automatically included.\n\n"
      "Multiple additional methods can be specified like\n"
      "<me:ExtraDOSCMethod xsi:type=\"me:HinderedRotorQM1D\">\n"
      "  data for method...\n"
      "</me:ExtraDOSCMethod>\nn";
    }
    virtual bool includesRotations(){return false;}

    //Get a pointer to a derived class by providing its id.
    static DensityOfStatesCalculator* Find(const std::string& id)
    {
      return dynamic_cast<DensityOfStatesCalculator*>(TopFind(id, typeID()));
    }

    // Read any data from XML and store in this instance. Default is do nothing.
    virtual bool ReadParameters(gDensityOfStates* gdos, PersistPtr ppDOSC=NULL){ return true; };

    // Provide a function to define particular counts of the DOS of a molecule
    virtual bool countCellDOS(gDensityOfStates* mol, size_t MaximumCell)=0;

    // Provide a function to calculate contribution to canonical partition function.
    virtual double canPrtnFnCntrb(gDensityOfStates* gdos, double beta)=0;

    // Provide a function to return the number of degrees of freedom associated with this count.
    virtual unsigned int NoDegOfFreedom(gDensityOfStates* gdos)=0;

    // Provide a function to calculate the zero point energy of a molecule.
	virtual double ZeroPointEnergy(gDensityOfStates* gdos) { return 0.0 ; } ;

    const Molecule* getParent() const {return m_parent;} ;
    void setParent(const Molecule* parent) { m_parent = parent;} ;

  private:
    const Molecule* m_parent;

    static const char* typeID() { return "Cell Density of States Calculators"; }
  };

}//namespace

#endif
