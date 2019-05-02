#ifndef GUARD_FTSTPotential_h
#define GUARD_FTSTPotential_h

//-------------------------------------------------------------------------------------------
//
// FTSTPotential.h
//
// Author: Struan Robertson
// Date:   29/Apr/2019
//
// This header file contains the declaration of the FTSTPotential abstract base class.
// Classes derived from this one calculate the potential, for a given relative orientation
// of fragments along the reaction coordinate, used in FTST calculations
//
//-------------------------------------------------------------------------------------------

#include "plugin.h"
#include <vector>

namespace mesmer
{
  class uFTST;

  class FTSTPotential : public TopPlugin {

  public:

    FTSTPotential() {};

    virtual ~FTSTPotential() {};

    static const char* typeID() { return "FTST Potential Energy Surfaces"; }
    virtual const char* getTypeID() { return typeID(); }

    //Get a pointer to a derived class by providing its id.
    static FTSTPotential* Find(const std::string& id)
    {
      return dynamic_cast<FTSTPotential*>(TopFind(id, typeID()));
    }

    uFTST* getParent() { return m_parent; }
    void setParent(uFTST* parent) { m_parent = parent; }

    // Initialize any general PES parameters.
    virtual void Initialize() = 0;

    // Set up parameters associated wth a value of the reaction coordinate.
    // This is usually called at the beginning of a phase integral calculation.
    virtual void RxnCrdInitialize(double rxnCrd) = 0;

    // The potential energy along the reaction coordinate.
    virtual double MEPPotential(double rxnCrd) = 0;

    // The hindering potential of the transitional modes perpendicular to 
    // the reaction coordinate.
    virtual double HinderingPotential(double rxnCrd, const std::vector<double>& angles) = 0;

  private:
    uFTST* m_parent;
  };

}

#endif