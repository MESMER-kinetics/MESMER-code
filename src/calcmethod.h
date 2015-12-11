#ifndef GUARD_CalcMethod_h
#define GUARD_CalcMethod_h

#include <map>
#include <string>
#include <vector>
#include "System.h"

namespace mesmer
{
  class System;
  /** Abstract base class as interface for plugin classes which do various
  types of execution, e.g. single calculations, gridsearch and fitting
  **/
  class CalcMethod : public TopPlugin
  {
  public:
    enum parseQuery { ALL, NOCONDITIONSOK, MODELPARAMS };
    static const char* typeID(){ return "Calculation methods"; }
    virtual ~CalcMethod(){}
    virtual const char* getTypeID(){return typeID();}

    // A CalcMethod can have its parsing behaviour customized in a number
    // of respects, which are discovered by calling this (virtual) function
    // with various parameters. Most methods return false to all questions.
    // Currently used as follows, but the way it is used is extendable.
    // - With parameter CalcMethod::ALL, a true return causes System::parse()
    //    not to parse <me:modelParameters> and other other such sections and
    //    leaves it to the CalcMethod. The CalcMethod UnitTests works like this.
    // - With parameter CalcMethod::MODELPARAMS, a true return shows that the
    //    calcMethod, e.g. ThermodynamicTable, will set its own defaults.
    // - With parameter CalcMethod::NOCONDITIONSOK, a true return prevents
    //    System::Parse() aborting when there is no <me:Conditions> section.
    virtual bool DoesOwnParsing(parseQuery q=ALL) { return false; }

    //Get a pointer to a derived class by providing its id.
    static CalcMethod* Find(const std::string& id)
    {
      return dynamic_cast<CalcMethod*>(TopFind(id, typeID()));
    }
    System* getParent() { return m_parent; }
    void setParent(System* parent) { m_parent = parent; }

    //Function to do the work
    virtual bool DoCalculation(System* pSys)=0;

  private:
    System* m_parent;
  };

}//namespace

#endif
