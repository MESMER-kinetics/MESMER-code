#ifndef GUARD_Tunneling_h
#define GUARD_Tunneling_h

#include <map>
#include "Molecule.h"

namespace mesmer
{

  /** Abstract base class for Tunneling calculators
  The derived concrete classes are plugin classes:
  -- New classes can be added without changing any of the existing code.
  They have a single global instance, the constructor of which registers
  the class with the base class. Subsequently, a pointer to the class is
  obtained by supplying the id (a string) to the Find function.
  **/
  class Reaction;
  
  class TunnelingCalculator
  {
  public:
    typedef std::map<std::string, TunnelingCalculator*> TunnelingMap;

    ///Base class constructor adds the derived class instance to the map
    TunnelingCalculator(const std::string& id) { get_Map()[id] = this; }

    //Get a pointer to a derived class by providing its id.
    static TunnelingCalculator* Find(const std::string& id)
    {
      TunnelingMap::iterator pos = get_Map().find(id);
      return (pos==get_Map().end()) ? NULL : pos->second;
    }

    virtual bool calculateCellTunnelingCoeffs(Reaction* pReact, std::vector<double>& TunnelingProbability) = 0 ;

  private:
    /// Returns a reference to the map of TunnelingCalculator classes
    /// Is a function rather than a static member variable to avoid initialization problems.

    static TunnelingMap& get_Map()
    {
      static TunnelingMap m;
      return m;
    }
  };

}//namespace

#endif
