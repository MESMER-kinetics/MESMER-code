#ifndef GUARD_MicroRate_h
#define GUARD_MicroRate_h

#include <map>
#include "XMLPersist.h"
#include "MesmerTools.h"
#include "Rdouble.h"
#include "plugin.h"

namespace mesmer
{

  class Reaction;

  /** Abstract base class for Microcanonical rate calculators
  The derived concrete classes are plugin classes:
  -- New classes can be added without changing any of the existing code.
  A new instance is made for each use, the constructor of which registers
  the class with the base class. Subsequently, a pointer to the class is
  obtained by supplying the id (a string) to the Find function.
  **/
  class MicroRateCalculator :public TopPlugin
  {
  public:

    virtual ~MicroRateCalculator() {}
    virtual const char* getTypeID(){return typeID();}
    static const char* typeID(){ return "Microcanonical Rate Calculators"; }

    //Get a pointer to a derived class by providing its id.
    static MicroRateCalculator* Find(const std::string& name)
    {
      return dynamic_cast<MicroRateCalculator*>(TopFind(name, typeID()));
    }
    Reaction* getParent() {return m_parent;} ;
    void setParent(Reaction* parent) { m_parent = parent;} ;

    //Get a list of the IDs of plugin classes derived from this one.
    static std::string List()
    {
      return TopPlugin::List(typeID(), comma);
    }

    virtual bool calculateMicroRateCoeffs(Reaction* pReact) = 0 ;

    virtual bool testMicroRateCoeffs(Reaction* pReact, PersistPtr ppbase) const;

    //@virtual bool ReadParameters(Reaction* pReac) = 0 ;

    virtual double get_ThresholdEnergy(Reaction* pReac) ;

    // Utility function to check for inconsistencies. 
    static bool ILTCheck(Reaction* pReac, PersistPtr ppReac) ;

  protected:
    Reaction* m_parent;
  };

}//namespace

#endif
