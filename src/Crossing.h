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

	class CrossingCalculator : public TopPlugin
	{
	public:
    ~CrossingCalculator(){}
  virtual const char* getTypeID(){return typeID();}


		//Get a pointer to a derived class by providing its id.
		static CrossingCalculator* Find(const std::string& id)
		{
      return dynamic_cast<CrossingCalculator*>(TopFind(id, typeID()));
    }
    Reaction* getParent() { return m_parent; }
    void setParent(Reaction* parent) { m_parent = parent; }

		virtual bool calculateCellCrossingCoeffs(Reaction* pReact, std::vector<double>& CrossingProbability) = 0 ;

		virtual bool ThereIsTunnellingWithCrossing(void) = 0;

private:
  Reaction* m_parent;
  static const char* typeID(){ return "Crossing Calculators"; }

	};

}//namespace

#endif
