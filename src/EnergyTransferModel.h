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
//  **/
//
//-------------------------------------------------------------------------------------------

#include <map>
#include <string>
//#include "XMLPersist.h"

namespace mesmer
{

  class EnergyTransferModel
  {
  public:
    typedef std::map<std::string, EnergyTransferModel*> EnergyTransferModelMap;

    // Base class constructor adds the derived class instance to the map
    EnergyTransferModel(const std::string& id) {
      get_Map()[id] = this;
      m_name = id;
    }

    virtual ~EnergyTransferModel(){}

    // Get a pointer to a derived class by providing its id.
    static EnergyTransferModel* Find(const std::string& id)
    {
      EnergyTransferModelMap::iterator pos = get_Map().find(id);
      return (pos==get_Map().end()) ? NULL : pos->second;
    } ;

    virtual std::string getName() const {return m_name;} ;

    virtual double calculateTransitionProbability(double Ei, double Ej) = 0 ;

    virtual bool ReadParameters() = 0 ; 

  private:
    // Returns a reference to the map of EnergyTransferModel classes
    // Is a function rather than a static member variable to avoid
    // initialization problems.
    static EnergyTransferModelMap& get_Map()
    {
      static EnergyTransferModelMap m;
      return m;
    }

    std::string m_name;
  };

}//namespace

#endif // GUARD_EnergyTransferModels_h
