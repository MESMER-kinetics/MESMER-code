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
    static const char* typeID(){ return "Calculation methods"; }
    virtual ~CalcMethod(){}
    virtual const char* getTypeID(){return typeID();}

    virtual bool DoesOwnParsing() { return false; }

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
