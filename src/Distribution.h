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
#include "plugin.h"

namespace mesmer
{

  /** Abstract base class for distribution calculators
  The derived concrete classes are plugin classes:
  -- New classes can be added without changing any of the existing code.
  They have a single global instance, the constructor of which registers
  the class with the base class. Subsequently, a pointer to the class is
  obtained by supplying the id (a string) to the Find function.
  **/
  class DistributionCalculator : public TopPlugin
  {
  public:
    virtual const char* getTypeID(){return typeID();}

    //Get a pointer to a derived class by providing its id.
    static DistributionCalculator* Find(const std::string& id)
    {
      return dynamic_cast<DistributionCalculator*>(TopFind(id, typeID()));
    }

    virtual bool calculateDistribution(std::vector<double> DOS, std::vector<double> ene, const double& beta, std::vector<double>& distribution) = 0 ;

private:
    static const char* typeID(){ return "Distribution Calculators"; }
  };

}//namespace

#endif
