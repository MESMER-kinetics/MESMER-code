#ifndef GUARD_MicroRate_h
#define GUARD_MicroRate_h

#include <map>
#include "XMLPersist.h"
#include "Rdouble.h"

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
  class MicroRateCalculator
  {
  public:
    typedef std::map<std::string, MicroRateCalculator*> MicroRateMap;

    ///Base class constructor adds the derived class instance to the map
    MicroRateCalculator(const std::string& id) :
      m_PreExp(0.0),
      m_NInf(0.0),
      m_TInf(298.0),
      m_EInf(0.0),
      m_usesILT(false), 
      m_isRvsILTpara(false)
    {
      get_Map()[id] = this;
      name = id;
    }

    virtual ~MicroRateCalculator() {}
    virtual MicroRateCalculator* Clone() = 0;

    //Get a pointer to a derived class by providing its id.
    static MicroRateCalculator* Find(const std::string& id)
    {
      MicroRateMap::iterator pos = get_Map().find(id);
      if(pos==get_Map().end())
        return NULL; 
      return (pos->second)->Clone();
    }

    virtual std::string getName();

    virtual bool calculateMicroRateCoeffs(Reaction* pReact) = 0 ;

    virtual bool testMicroRateCoeffs(Reaction* pReact, PersistPtr ppbase) const;
  
    virtual bool ReadParameters(Reaction* pReac); //Used by MesmerILT and SimpleILT
    
    virtual double get_ThresholdEnergy(Reaction* pReac) ;
    
    bool isReverseReactionILT_Ea() {return m_isRvsILTpara;}

  private:
    /// Returns a reference to the map of MicroRateCalculator classes
    /// Is a function rather than a static member variable to avoid initialization problems.
    static MicroRateMap& get_Map()
    {
      static MicroRateMap m;
      return m;
    }

    std::string name;
    
 protected:   

    // All the parameters that follow are for an arrhenius expression of the type:
    // k(T) = Ainf*(T/Tinf)^ninf * exp(-Einf/(RT))
    Rdouble m_PreExp ;           // Preexponetial factor
    Rdouble m_NInf ;             // Modified Arrhenius parameter
    double  m_TInf ;             // T infinity
    Rdouble m_EInf ;             // E infinity
    bool    m_usesILT;
    bool    m_isRvsILTpara;      // The ILT parameters provided are for reverse direction.
  };

}//namespace

#endif
