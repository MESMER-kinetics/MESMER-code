#ifndef GUARD_CalcMethod_h
#define GUARD_CalcMethod_h

#include <map>
#include <string>
#include <vector>

namespace mesmer
{
  class System;
  /** Abstract base class as interface for plugin classes which do various
  types of execution, e.g. single calculations, gridsearch and fitting
  **/
  class CalcMethod
  {
  public:
    typedef std::map<std::string, CalcMethod*> CalcMethodMap;

    ///Base class constructor adds the derived class instance to the map
    CalcMethod(const std::string& id){ get_Map()[id] = this; }

    virtual ~CalcMethod(){}

    //Get a pointer to a derived class by providing its id.
    static CalcMethod* Find(const std::string& id)
    {
      CalcMethodMap::iterator pos = get_Map().find(id);
      return (pos==get_Map().end()) ? NULL : pos->second;
    }

    //Function to do the work
    virtual bool DoCalculation(System* pSys)=0;

    //Parses the <me:control> section of the XML input file to find the specified method
    //For instance: <calcMethod>simpleCalc</calcMethod>
    //If there is no <calcMethod> element, the simpleCalc, set in defaults.xml, is used.
    static CalcMethod* GetCalcMethod(PersistPtr ppControl)
    {
      const char* type = ppControl->XmlReadValue("me:calcMethod"); //or uses default
      return Find(type);
    }

  private:
    /// Returns a reference to the map of CalcMethod classes
    /// Is a function rather than a static member variable to avoid initialization problems.
    static CalcMethodMap& get_Map()
    {
      static CalcMethodMap m;
      return m;
    }
  };

}//namespace

#endif
