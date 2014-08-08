#ifndef GUARD_Tunneling_h
#define GUARD_Tunneling_h

#include <map>

namespace mesmer
{

  /** Abstract base class for Tunneling calculators
  The derived concrete classes are plugin classes:
  **/
  class Reaction;
  
  class TunnelingCalculator : public TopPlugin
  {
  public:
    static const char* typeID(){ return "Tunneling Calculators"; }
    virtual const char* getTypeID(){return typeID();}

    //Get a pointer to a derived class by providing its id.
    static TunnelingCalculator* Find(const std::string& id)
    {
      return dynamic_cast<TunnelingCalculator*>(TopFind(id, typeID()));
    }
    Reaction* getParent() { return m_parent; }
    void setParent(Reaction* parent) { m_parent = parent; }

    virtual bool calculateCellTunnelingCoeffs(Reaction* pReact, std::vector<double>& TunnelingProbability) = 0 ;

private:
  Reaction* m_parent;
  };

}//namespace

#endif
