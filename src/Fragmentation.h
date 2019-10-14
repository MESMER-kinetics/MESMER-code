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
#include "Reaction.h"

namespace mesmer
{

//
// Abstract base class for the calculation of fragment distribution on dissocistion.
//
  class FragDist
  {
  public:

    // Constructors.
    FragDist() {};

    // Destructor.
    virtual ~FragDist() {};

    // Initialize the fragment distribution.
    virtual void initialize(Reaction* pReaction) = 0;

    // Calculate distribution.
    virtual void calculate(double excessEnergy, std::vector<double>& dist, size_t size) = 0;

    // Return resources
    virtual void clear() = 0;

  };

}//namespace

#endif

