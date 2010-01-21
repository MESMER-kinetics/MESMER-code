#ifndef GUARD_EnergyTransferModels_h
#define GUARD_EnergyTransferModels_h

//-------------------------------------------------------------------------------------------
//
// EnergyTransferModel.h
//
// Author: Struan Robertson
// Date:   27/Jun/2009
//
// This header file contains the declaration of the EnergyTransferModel abstract
// base class. This class will be inherited by all energy transfer models and
// maintains a list of all models.
//
//  The derived concrete classes are plugin classes:
//  -- New classes can be added without changing any of the existing code.
//  They have a single global instance, the constructor of which registers
//  the class with the base class. Subsequently, a pointer to the class is
//  obtained by supplying the id (a string) to the Find function.
//
//-------------------------------------------------------------------------------------------

#include <map>
#include <string>

namespace mesmer
{
  // Forward declarations:
  
 class Molecule;
  
  class EnergyTransferModel
  {
  public:
    typedef std::map<std::string, EnergyTransferModel*> EnergyTransferModelMap;

    // Base class constructor adds the derived class instance to the map
    EnergyTransferModel(const std::string& id) : m_name(id),m_parent(NULL){
      get_Map()[id] = this;
    }

    virtual ~EnergyTransferModel(){}

    virtual EnergyTransferModel* Clone() = 0;

    // Get a pointer to a new instance of a derived class by providing its id.
    // Returns NULL if not recognized.
    static EnergyTransferModel* Find(const std::string& id)
    {
      EnergyTransferModelMap::iterator pos = get_Map().find(id);
      if(pos==get_Map().end())
        return NULL; 
      return (pos->second)->Clone();
    } ;

    virtual double calculateTransitionProbability(double Ei, double Ej) = 0 ;

    virtual bool ReadParameters(const Molecule* parent) = 0 ; 

    std::string getName() const {return m_name;} ;
    const Molecule* getParent() const {return m_parent;} ;
    void setParent(const Molecule* parent) { m_parent = parent;} ;


  private:
    // Returns a reference to the map of EnergyTransferModel classes
    // Is a function rather than a static member variable to avoid
    // initialization problems.
    static EnergyTransferModelMap& get_Map()
    {
      static EnergyTransferModelMap m;
      return m;
    }

  private:
    const std::string m_name;
    const Molecule* m_parent;
  };

}//namespace

#endif // GUARD_EnergyTransferModels_h
