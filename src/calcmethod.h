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

    //Parses the <me:control> section of the XML input file to find the specified method
    //For instance: <calcMethod>simpleCalc</calcMethod>
    //If there is no <calcMethod> element, the simpleCalc, set in defaults.xml, is used.
    //If there is an unrecognized method, the acceptable method are listed in an error message
    //and returns false.
    static CalcMethod* GetCalcMethod(PersistPtr ppControl, std::string& name)
    {
      const char* type = ppControl->XmlReadValue("me:calcMethod"); //or uses default
      name = std::string(type) ;
      CalcMethod* method = Find(type);
      return method;
    }
  private:
    System* m_parent;
    static const char* typeID(){ return "Calculation methods"; }

  };

}//namespace

#endif
