#ifndef GUARD_DensityOfStates_h
#define GUARD_DensityOfStates_h

#include <map>
#include <string>
#include <vector>

namespace mesmer
{
  class ModelledMolecule;
  /** Abstract base class for cell Density Of States (DOS) calculators for energy grained master equation (EGME).
  The derived concrete classes are plugin classes:
  -- New classes can be added without changing any of the existing code.
  They have a single global instance, the constructor of which registers
  the class with the base class. Subsequently, a pointer to the class is
  obtained by supplying the id (a string) to the Find function.
  **/
  class DensityOfStatesCalculator
  {
  public:
    typedef std::map<std::string, DensityOfStatesCalculator*> DensityOfStatesMap;

    ///Base class constructor adds the derived class instance to the map
    DensityOfStatesCalculator(const std::string& id) { get_Map()[id] = this; }

    //Get a pointer to a derived class by providing its id.
    static DensityOfStatesCalculator* Find(const std::string& id)
    {
      DensityOfStatesMap::iterator pos = get_Map().find(id);
      return (pos==get_Map().end()) ? NULL : pos->second;
    }

    // Provide a function to define particular counts of the convolved DOS of two molecules.
    virtual bool countDimerCellDOS(ModelledMolecule* p_mol1, ModelledMolecule*  p_mol2, std::vector<double>& rctsCellEne, std::vector<double>& rctsCellDOS) = 0; 

    // provide a function to define particular counts of the DOS of a molecule
    virtual bool countCellDOS(ModelledMolecule* mol) = 0;

  private:
    /// Returns a reference to the map of DensityOfStatesCalculator classes
    /// Is a function rather than a static member variable to avoid initialization problems.
    static DensityOfStatesMap& get_Map()
    {
      static DensityOfStatesMap m;
      return m;
    }
  };

}//namespace

#endif
