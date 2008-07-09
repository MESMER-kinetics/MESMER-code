//-------------------------------------------------------------------------------------------
//
// DistributionCalculator.h
//
// Author: Chi-Hsiu Liang
// Date:   _2008_05_15_
//
// This file contains the definition of the distribution calculator plug-in class.
//
//-------------------------------------------------------------------------------------------
#ifndef GUARD_Distribution_h
#define GUARD_Distribution_h

#include <map>

namespace mesmer
{

  /** Abstract base class for distribution calculators
  The derived concrete classes are plugin classes:
  -- New classes can be added without changing any of the existing code.
  They have a single global instance, the constructor of which registers
  the class with the base class. Subsequently, a pointer to the class is
  obtained by supplying the id (a string) to the Find function.
  **/
  class DistributionCalculator
  {
  public:
    typedef std::map<std::string, DistributionCalculator*> DistributionMap;

    ///Base class constructor adds the derived class instance to the map
    DistributionCalculator(const std::string& id) { get_Map()[id] = this; }

    //Get a pointer to a derived class by providing its id.
    static DistributionCalculator* Find(const std::string& id)
    {
      DistributionMap::iterator pos = get_Map().find(id);
      return (pos==get_Map().end()) ? NULL : pos->second;
    }

    virtual bool calculateDistribution(std::vector<double> DOS, std::vector<double> ene, const double& beta, std::vector<double>& distribution) = 0 ;

  private:
    /// Returns a reference to the map of DistributionCalculator classes
    /// Is a function rather than a static member variable to avoid initialization problems.
    static DistributionMap& get_Map()
    {
      static DistributionMap m;
      return m;
    }
  };

}//namespace

#endif
