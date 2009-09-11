#ifndef GUARD_Crossing_h
#define GUARD_Crossing_h

#include <map>

namespace mesmer
{

	/** Abstract base class for Crossing calculators
	The derived concrete classes are plugin classes:
	-- New classes can be added without changing any of the existing code.
	They have a single global instance, the constructor of which registers
	the class with the base class. Subsequently, a pointer to the class is
	obtained by supplying the id (a string) to the Find function.
	**/
	class Reaction;  //to avoid compile time errors, this declaration tells the compiler that 
	//Reaction is a class; reaction objects are featured in some functions below

	class CrossingCalculator
	{
	public:
		typedef std::map<std::string, CrossingCalculator*> CrossingMap;

		///Base class constructor adds the derived class instance to the map
		CrossingCalculator(const std::string& id) { get_Map()[id] = this; }

		//Get a pointer to a derived class by providing its id.
		static CrossingCalculator* Find(const std::string& id)
		{
			CrossingMap::iterator pos = get_Map().find(id);
			return (pos==get_Map().end()) ? NULL : pos->second;
		}

		virtual bool calculateCellCrossingCoeffs(Reaction* pReact, std::vector<double>& CrossingProbability) = 0 ;

	private:
		/// Returns a reference to the map of CrossingCalculator classes
		/// Is a function rather than a static member variable to avoid initialization problems.

		static CrossingMap& get_Map()
		{
			static CrossingMap m;
			return m;
		}
	};

}//namespace

#endif
