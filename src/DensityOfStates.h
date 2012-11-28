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
    virtual const char* getTypeID() override {return typeID(false);}
    virtual const char* getTypeID(bool extra) override {return typeID(extra);}


    //Get a pointer to a derived class by providing its id.
    static DensityOfStatesCalculator* Find(const std::string& id, bool extraType=false)
    {
      return dynamic_cast<DensityOfStatesCalculator*>(TopFind(id, typeID(extraType)));
    }

    // Read any data from XML and store in this instance. Default is do nothing.
    virtual bool ReadParameters(gDensityOfStates* gdos, PersistPtr ppDOSC=NULL){ return true; };

    // Provide a function to define particular counts of the DOS of a molecule
    virtual bool countCellDOS(gDensityOfStates* mol, size_t MaximumCell)=0;

    // Provide a function to calculate contribution to canonical partition function.
    virtual double canPrtnFnCntrb(gDensityOfStates* gdos, double beta)=0;

    // Provide a function to return the number of degrees of freedom associated with this count.
    virtual unsigned int NoDegOfFreedom(gDensityOfStates* gdos)=0;

    const Molecule* getParent() const {return m_parent;} ;
    void setParent(const Molecule* parent) { m_parent = parent;} ;

  private:
    const Molecule* m_parent;

  private:
    static const char* typeID(bool extraType)
    {  return extraType ? "Cell Density of States Calculators (Extra)"
        : "Cell Density of States Calculators";
    }
    /* There are separate maps in TopPlugin for normal and extra Density of
       States Calculators, which are Listed separately, but under the same
       title.
    */

  };

}//namespace

#endif
