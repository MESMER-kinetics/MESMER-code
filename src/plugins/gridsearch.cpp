// gridsearch.cpp
// Calculates for all combinations of values of range variables.
// Each calculation has a range of pressure/temperature conditions.

#include "../System.h"
#include "../calcmethod.h"

namespace mesmer
{
  class GridSearch : public CalcMethod
  {
  public:
    GridSearch(const char* id) : m_id(id) { Register();}
    virtual ~GridSearch() {}
    virtual const char* getID()  { return m_id; }
    virtual bool DoCalculation(System* pSys);

  private:
    bool DoRangeCalcs(unsigned startvar, System* pSys);
    bool CalcAndReport(System* pSys);
  private:
    const char* m_id;
  };

  ////////////////////////////////////////////////
  //Global instance
  GridSearch theGridSearch("gridSearch");
  ///////////////////////////////////////////////

  inline bool GridSearch::DoCalculation(System* pSys)
  {
    //Set the values of all range variables to the bottom of their range
    for(size_t i=0;i!=Rdouble::withRange().size();++i)
      Rdouble::withRange()[i]->set_to_lower();

    return DoRangeCalcs(0, pSys);
  }

  //Calculate for all values of range variables later than startvar (recursive function)
  bool GridSearch::DoRangeCalcs(unsigned startvar, System* pSys)
  {
    //Do calculation when the last range variable is incrementing
    if(startvar >= Rdouble::withRange().size())
      return CalcAndReport(pSys);
    //Do later variables exhaustively at each iteration
    do 
    { 
      DoRangeCalcs(startvar+1, pSys);
    } while(!IsNan(++(*Rdouble::withRange()[startvar])));
    return true;
  }

  bool GridSearch::CalcAndReport(System* pSys)
  { 
    double chiSquare(1000.0);
    ctest << "Parameter Grid\n";

    for(size_t i=0;i!=Rdouble::withRange().size();++i)
    {
      cerr  << ' ' << Rdouble::withRange()[i]->get_varname() << '=' << Rdouble::withRange()[i]->originalUnits();
      ctest << ' ' << Rdouble::withRange()[i]->get_varname() << '=' << Rdouble::withRange()[i]->originalUnits();
    }
    cerr << endl;

    ctest << "\n{" << endl;
    pSys->calculate(chiSquare);
    ctest << "chiSquare = " << chiSquare << " )\n}\n";

    return true;
  }

}//namespace

