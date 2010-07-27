#ifndef GUARD_MicroRate_h
#define GUARD_MicroRate_h

#include <map>
#include "XMLPersist.h"

namespace mesmer
{

  class Reaction;

  /** Abstract base class for Microcanonical rate calculators
  The derived concrete classes are plugin classes:
  -- New classes can be added without changing any of the existing code.
  A new instance is made for each use, the constructor of which registers
  the class with the base class. Subsequently, a pointer to the class is
  obtained by supplying the id (a string) to the Find function.
  **/
  class MicroRateCalculator
  {
  public:
    typedef std::map<std::string, MicroRateCalculator*> MicroRateMap;

    ///Base class constructor adds the derived class instance to the map
    MicroRateCalculator(const std::string& id) {
      get_Map()[id] = this;
      name = id;
    }

    virtual ~MicroRateCalculator() {}
    virtual MicroRateCalculator* Clone() = 0;

    //Get a pointer to a derived class by providing its id.
    static MicroRateCalculator* Find(const std::string& id)
    {
      MicroRateMap::iterator pos = get_Map().find(id);
      if(pos==get_Map().end())
        return NULL; 
      return (pos->second)->Clone();
    }

    virtual std::string getName();

    virtual bool calculateMicroRateCoeffs(Reaction* pReact) = 0 ;

    virtual bool testMicroRateCoeffs(Reaction* pReact, PersistPtr ppbase) const;
  
    virtual bool ReadParameters(Reaction* pReac) = 0 ;
    
    virtual double get_ThresholdEnergy(Reaction* pReac) ;
    
  private:
    /// Returns a reference to the map of MicroRateCalculator classes
    /// Is a function rather than a static member variable to avoid initialization problems.
    static MicroRateMap& get_Map()
    {
      static MicroRateMap m;
      return m;
    }

    std::string name;

  };

}//namespace

#endif
