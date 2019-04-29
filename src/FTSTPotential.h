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

    virtual void initialize(double rxnCrd) = 0;

    virtual double MEPPotential(double rxnCrd) = 0;

    virtual double HinderingPotential(double rxnCrd, const std::vector<double>& angles) = 0;

  private:
    uFTST* m_parent;
  };

}

#endif