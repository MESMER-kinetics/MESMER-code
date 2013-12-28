#ifndef GUARD_plugin_h
#define GUARD_plugin_h
//-------------------------------------------------------------------------------------------
//
// System.h
//
// Author: Chris Morley
// Date:   24 Nov 2012
//
// This header file contains the declaration of TopPlugin, an abstract class
// which is the parent of all the plugin types (e.g. DensityOfStatesCalculator)
// and the grandparent of the concrete plugin classes (e.g. Classical Rotor).
// It maintains a list of plugins and provides an instance of a plugin from an ID.
// This search is not dependent on case or the presence of leading, trailing or
// embedded whitespace in the ID.
// TopPlugin::List provides a table of plugins in short or verbose form and is
// called from the command line with  mesmer -t  and mesmer -T respectively.
//-------------------------------------------------------------------------------------------

#include <string>
#include <map>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <iostream>
#include "Persistence.h"
namespace mesmer
{

  //Parent of all plugin-type classes (e.g. DensityOfStatesCalculator)
  class TopPlugin
  {
  public:
    typedef std::map<std::string, std::map<std::string, TopPlugin*> > PluginMapType;
    typedef std::map<std::string, TopPlugin*> TypeMapType;

    //IDs on separate lines, IDs+descriptions on separate lines, comma separated
    enum format {brief, verbose, comma};

    virtual ~TopPlugin() {}

    //getID() is public because sometimes used externally
    virtual const char* getID()=0;

    //Get the plugin to parse its data. Plugins with no data use the default here.
    virtual bool ParseData(PersistPtr pp){ return true; }

    // Print out full case versions of ID and TypeID
    // and additionally Description() and typeDescription() if verbose is true.
    // For a single plugintype:
    static std::string List(const char* pluginType, format f=brief)
    {
      std::stringstream ss;
      PluginMapType::iterator titer = allPlugins().find(pluginType);
      if(titer==allPlugins().end())
        return std::string(); //empty
      else
        return ListofType(titer->second, ss, f);
    }

    //For allPlugin types:
    static std::string List(format f=brief)
    {
      std::stringstream ss;
      PluginMapType::iterator titer=allPlugins().begin();
      for(; titer!=allPlugins().end(); ++titer)
        ListofType(titer->second, ss, f);
      return ss.str();
    }

  private:
    static std::string ListofType(TypeMapType typeMap, std::stringstream& ss, format f)
    {
      for(TypeMapType::iterator piter=typeMap.begin(); piter!=typeMap.end(); ++piter)
      {
        if(piter==typeMap.begin() && f!=comma)
        {
          ss << piter->second->getTypeID() << '\n';
          if(f==verbose)
            ss << piter->second->typeDescription() << '\n';
        }
        if(f==comma)
        {
          if(piter!=typeMap.begin())
            ss << ", ";
          ss  << piter->second->getID();
        }
        else
          ss << " * " << piter->second->getID() << "\n";

        if(f==verbose)
          ss << "  " << piter->second->Description() << '\n';
      }
      return ss.str();
    }
  public:
    virtual TopPlugin* Clone(){ return NULL; }
  protected:

    static TopPlugin* TopFind(const std::string& name,const std::string& pluginType,
      bool useErrorMessage=true)
    {
      // name is converted to lower case so that the search is case insensitive.
      TypeMapType typeMap = allPlugins()[pluginType];
      TypeMapType::iterator pos = typeMap.find(toLowercase(name));
      if(pos==typeMap.end())
      {
        if(useErrorMessage)
          std::cerr << "\"" << name << "\" not recognized as one of the \n "
          << pluginType << " plugins. Possibilities are:\n "
          << List(pluginType.c_str(), comma) << std::endl;
        return NULL;
      }

      //Return a new instance, or the original instance if there is no Clone function.
      TopPlugin* newp = (pos->second)->Clone();
      return newp ? newp : pos->second;
    }

    virtual const char* getTypeID()=0;
    virtual const char* Description() { return ""; };
    virtual const char* typeDescription() { return ""; };

    void Register()
    {
      // The id strings in allPlugins are lower case/no white space versions
      std::string cid = getTypeID();
      std::string id  = getID();
      allPlugins()[cid][toLowercase(id)] = this;
    }

  private:
    static PluginMapType& allPlugins()
    {
      static PluginMapType m;
      return m;
    }

    static std::string toLowercase(const std::string& txt)
    {
      //Also removes whitespace
      std::string lctxt;
      std::remove_copy_if(txt.begin(), txt.end(), back_inserter(lctxt), ::isspace);
      std::transform(lctxt.begin(), lctxt.end(), lctxt.begin(), ::tolower);
      return lctxt;
    }

  };

}//namespace

#endif

