//-------------------------------------------------------------------------------------------
//
// Fragmentation.h
//
// Author: Struan Robertson
// Date:   8th October 2019
//
// This file contains the definition of the fragmentation abstract base class.
//
// This class definition should probably be the base class to a set of plug-in classes.
//
//-------------------------------------------------------------------------------------------
#ifndef GUARD_Fragmentation_h
#define GUARD_Fragmentation_h

#include <vector>
#include "plugin.h"

namespace mesmer
{
  // Forward class definitions
  class Reaction;

//
// Abstract base class for the calculation of fragment distribution on dissocistion.
//
  class FragDist : public TopPlugin
  {
  public:

    // Constructors.
    FragDist() {};

    // Destructor.
    virtual ~FragDist() {};

    static const char* typeID() { return "Fragmentation Energy Partition Calculators"; }
    virtual const char* getTypeID() { return typeID(); }
    virtual const char* typeDescription() {
      return "This set of plugins describes how energy is distributed among the fragments \n"
             "produced when an adduct (whch might be a transition state) disintegrates \n";
    }

    //Get a pointer to a derived class by providing its id.
    static FragDist* Find(const std::string& id)
    {
      return dynamic_cast<FragDist*>(TopFind(id, typeID()));
    }
    Reaction* getParent() { return m_parent; }
    void setParent(Reaction* parent) { m_parent = parent; }

    // Read any data from XML and store in this instance. Default is do nothing.
    virtual bool ReadParameters(PersistPtr ppFragDist, std::string name) { return true; };

    // Initialize the fragment distribution.
    virtual void initialize(Reaction* pReaction) = 0;

    // Calculate distribution.
    virtual void calculate(double excessEnergy, std::vector<double>& dist) = 0;

    // Return resources
    virtual void clear() = 0;

  private:
    Reaction* m_parent;
  };

}//namespace

#endif

