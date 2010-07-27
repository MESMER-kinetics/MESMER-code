#ifndef GUARD_EnergyTransferModels_h
#define GUARD_EnergyTransferModels_h

//-------------------------------------------------------------------------------------------
//
// EnergyTransferModel.h
//
// Author: Struan Robertson
// Date:   27/Jun/2009
//
/*****************************************************************************
This header file contains the declaration of the EnergyTransferModel abstract
base class. This class will be inherited by all energy transfer models and
maintains a list of all models.

The derived concrete classes are plugin classes- new classes can be added
without changing any of the existing code. ExponentialDown.cpp is a file
containing a derived class and has lots of comments.

The constructor of the derived classes call the base class constructor which adds
the derived class to its list. Subsequently, supplying the id (a string)
to the Find function provides a pointer to a new instance. For historical reasons
not all plugin classes use a new instance on every call to Find(), but it is worth
doing because it allows them to contain data perttaining to a particular molecule
or reaction.
*****************************************************************************/

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

  /****************************************************************************
   EnergyTransferModel::Find() returns a pointer to a derived class from its id,
   or NULL if not recognized.
   In this class a new instance is returned, (a copy of a global instance)
   so that molecule-specific data can be stored in it. Some types of plugin use
   only the single global instance.
   It is the responsibility of the calling function to delete the instance.

   In this particular case the calling function is
   gWellProperties::initialization(). It parses the XML data file for a
   <me:energyTransferModel> element for each 'modelled' molecule. If not found
   the default id from defaults.xml is used, and this logged in the log file.
   EnergyTransferModel::Find(id) is called and then the derived class's
   ReadParameters(). The derived EnergyTransferModel class instance is deleted
   in the gWellProperties destructor.
  ****************************************************************************/
  static EnergyTransferModel* Find(const std::string& id)
  {
    EnergyTransferModelMap::iterator pos = get_Map().find(id);
    if(pos==get_Map().end())
      return NULL; 
    return (pos->second)->Clone();
  } ;

  /*************************************************************
  This is the function which does the real work of the plugin
  *************************************************************/
  virtual double calculateTransitionProbability(double Ei, double Ej) = 0 ;

  virtual bool ReadParameters(const Molecule* parent) = 0 ; 

  std::string getName() const {return m_name;} ;
  const Molecule* getParent() const {return m_parent;} ;
  void setParent(const Molecule* parent) { m_parent = parent;} ;


private:
  /**************************************************************************
  Returns a reference to the map of EnergyTransferModel classes.
  Is a function rather than a static member variable to avoid initialization
  problems.
  **************************************************************************/
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
