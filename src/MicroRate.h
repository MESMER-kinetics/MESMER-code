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
  They have a single global instance, the constructor of which registers
  the class with the base class. Subsequently, a pointer to the class is
  obtained by supplying the id (a string) to the Find function.
  **/
  class MicroRateCalculator
  {
  public:
    typedef std::map<std::string, MicroRateCalculator*> MicroRateMap;

    ///Base class constructor adds the derived class instance to the map
    MicroRateCalculator(const std::string& id) { get_Map()[id] = this; }

    //Get a pointer to a derived class by providing its id.
    static MicroRateCalculator* Find(const std::string& id)
    {
      MicroRateMap::iterator pos = get_Map().find(id);
      return (pos==get_Map().end()) ? NULL : pos->second;
    }

    virtual bool calculateMicroRateCoeffs(Reaction* pReact, std::vector<double> &cellKfmc, const MesmerEnv &mEnv) = 0 ;

    virtual bool testMicroRateCoeffs(Reaction* pReact, std::vector<double> &cellKfmc, PersistPtr ppbase, const MesmerEnv &mEnv) const;

  private:
    /// Returns a reference to the map of MicroRateCalculator classes
    /// Is a function rather than a static member variable to avoid initialization problems.
    static MicroRateMap& get_Map()
    {
      static MicroRateMap m;
      return m;
    }
  };

}//namespace

#endif
