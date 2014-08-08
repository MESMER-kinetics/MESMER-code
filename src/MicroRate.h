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

/*****************************************************************************
This header file contains the declaration of the MicroRateCalculator abstract
class. This class will be inherited by all microcanonical rate models.
Plugin types like this are usually in .h files because they are mainly interface.

The derived concrete classes are plugin classes - new classes can be added
without changing any of the existing code. MesmerILT.cpp is a file
containing a derived class and has lots of comments.
Plugin classes like MesmerILT.cpp usually do not have a .h file because they
are called through their parent's interface, and never directly.

At startup the constructors of the derived classes register their presence and
their IDs, and this is passed to their grandparent class TopPlugin which maintains
a list of all plugins. This list can be displayed from the command line using:
  mesmer -t   or   mesmer -T
With the -T option the output includes descriptions of the plugins classes and
the plugin types, for those classes which have been provided with a Description()
function.

The normal way of using a plugin is to call the global function ParseForPlugin,
as in Reaction.cpp near line 250:

m_pMicroRateCalculator = ParseForPlugin<MicroRateCalculator>(this, "me:MCRCMethod", ppReac);

The template parameter is the plugin-type class and the other parameters are a
pointer to the parent, the XML type name and a PersistPtr, which is a position in
the XML file. ParseForPlugin will parse the XML input file for the specified
plugin type, making an instance of the requested plugin from its name and parsing
for its data using the plugin's virtual function. (See MesmerILT.cpp.) It can also
insert a default plugin name from defaults.xml (if appropriate) and does most of
the reporting of errors, like missing or unrecognized names, etc.

Alternatively, supplying an ID (a string like "MesmerILT") to the Find function of
the plugin type class (MicroRateCalculator here) provides a pointer to a new
instance of the plugin class.
*****************************************************************************/

  class MicroRateCalculator :public TopPlugin
  {
  public:

    virtual ~MicroRateCalculator() {}
    virtual const char* getTypeID(){return typeID();}

    /*****************************************************************************
    The name of the plugin type that appears in the list on the command line.
    *****************************************************************************/
    static const char* typeID(){ return "Microcanonical Rate Calculators"; }

    // Get a pointer to a derived class by providing its id.
    static MicroRateCalculator* Find(const std::string& name)
    {
      return dynamic_cast<MicroRateCalculator*>(TopFind(name, typeID()));
    }
    Reaction* getParent() {return m_parent;} ;
    void setParent(Reaction* parent) { m_parent = parent;} ;

    // Get a list of the IDs of plugin classes derived from this one.
    static std::string List()
    {
      return TopPlugin::List(typeID(), comma);
    }

/*****************************************************************************
Any new plugin-type class should contain the above functions. It will also need
at least one virtual function like calculateMicroCnlFlux() which does the main work.
The m_parent member variable should be Molecule*, Reaction* or System*, as
appropriate, and the functions getParent and setParent and should be compatible.
*****************************************************************************/

   // Method to calculate the microcanoical flux (W(E)/h) at the transition state.
    virtual bool calculateMicroCnlFlux(Reaction* pReact) = 0 ;

    virtual bool testMicroRateCoeffs(Reaction* pReact, PersistPtr ppbase) const;

    virtual double get_ThresholdEnergy(Reaction* pReac) ;

    // Utility function to check for inconsistencies. 
    static bool ILTCheck(Reaction* pReac, PersistPtr ppReac) ;

  protected:
    Reaction* m_parent;
  };

}//namespace

#endif
