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
  GridSearch(const std::string& id) : CalcMethod(id) {}
  virtual ~GridSearch() {}
  virtual bool DoCalculation(System* pSys);

private:
  bool DoRangeCalcs(unsigned startvar, System* pSys);
  bool CalcAndReport(System* pSys);
};

////////////////////////////////////////////////
//Global instance
GridSearch theGridSearch("gridSearch");
///////////////////////////////////////////////

inline bool GridSearch::DoCalculation(System* pSys)
{
  return DoRangeCalcs(0, pSys);
}

//Calculate for all values of range variables later than startvar (recursive function)
bool GridSearch::DoRangeCalcs(unsigned startvar, System* pSys)
{
  //Do calculation when the last range variable is incrementing
  if(startvar >= Rdouble::withRange().size())
    return CalcAndReport(pSys);
  //Do later variables exhaustively at each iteration
  do { DoRangeCalcs(startvar+1, pSys); } while(!IsNan(++(*Rdouble::withRange()[startvar])));
  return true;
}

bool GridSearch::CalcAndReport(System* pSys)
{ 
  double chiSquare(1000.0);
  ctest << "Parameter Grid\n";

  for(size_t i=0;i!=Rdouble::withRange().size();++i)
  {
    cerr  << ' ' << Rdouble::withRange()[i]->get_varname() << '=' << *Rdouble::withRange()[i];
    ctest << ' ' << Rdouble::withRange()[i]->get_varname() << '=' << *Rdouble::withRange()[i];
  }
  cerr << endl;

  ctest << "\n{" << endl;
  pSys->calculate(chiSquare);
  ctest << "chiSquare = " << chiSquare << " )\n}\n";

  return true;
}

}//namespace

/*
    // produce a grid for search
    db2D gridArray; // this array grows up freely without needing dimensions

    int totalSteps = 1, dataPointSize = int(Rdouble::withRange().size());
    for (int varID(0); varID < dataPointSize; ++varID){
      // for every dimension create a series and duplicate the serie while the indices going forward
      double lower(0.0), upper(0.0), stepsize(0.0);
      Rdouble::withRange()[varID]->get_range(lower, upper, stepsize);
      int numSteps = Rdouble::withRange()[varID]->get_numsteps();
      totalSteps *= numSteps;
    }

    int spanSteps = 1;
    for (int varID(0); varID < dataPointSize; ++varID){
      double lower(0.0), upper(0.0), stepsize(0.0);
      Rdouble::withRange()[varID]->get_range(lower, upper, stepsize);
      int numSteps = Rdouble::withRange()[varID]->get_numsteps();
      int stack(0), block(0);
      while (stack < totalSteps){
        int step(0);
        while(step < spanSteps){
          double stepValue = lower + block * stepsize;
          gridArray[stack][varID] = stepValue;
          ++step;
          ++stack;
        }
        ++block;
        if (block == numSteps) block = 0;
      }
      spanSteps *= numSteps;
    }

    // TimeCount events; unsigned int timeElapsed;
    int calPoint(0);

    ofstream punchStream;
    if (m_Flags.searchMethod == GRIDSEARCHWITHPUNCH) 
      punchStream.open(m_Flags.punchFileName.c_str());

    for (int i(0); i < totalSteps; ++i){
      double chiSquare(1000.0);

      // assign values
      for (int varID(0); varID < dataPointSize; ++varID) *Rdouble::withRange()[varID] = gridArray[i][varID];

      // calculate
      cerr << "Parameter Grid " << calPoint <<endl;;
      ctest << "Parameter Grid " << calPoint << "\n{\n";
      calculate(chiSquare);

      if (dataPointSize){
        ctest << "Parameters: ( ";
        if (m_Flags.punchSymbols.size()){
          for (int varID(0); varID < dataPointSize; ++varID){
            if (m_Flags.searchMethod == GRIDSEARCHWITHPUNCH) punchStream << "Para" << varID << "\t";
          }
          if (m_Flags.searchMethod == GRIDSEARCHWITHPUNCH) punchStream << "Temperature (K)\tNumber density\t";
          if (m_Flags.searchMethod == GRIDSEARCHWITHPUNCH) punchStream << m_Flags.punchSymbols;
          m_Flags.punchSymbols.clear();
        }
        for (int varID(0); varID < dataPointSize; ++varID){
          ctest << gridArray[i][varID] << " ";
          if (m_Flags.searchMethod == GRIDSEARCHWITHPUNCH) punchStream << gridArray[i][varID] << "\t";
        }
        if (m_Flags.searchMethod == GRIDSEARCHWITHPUNCH) punchStream << m_Env.beta << "\t" << m_Env.conc << "\t";
        if (m_Flags.searchMethod == GRIDSEARCHWITHPUNCH) punchStream << m_Flags.punchNumbers;
        m_Flags.punchNumbers.clear();

        if (m_Flags.searchMethod == GRIDSEARCHWITHPUNCH) punchStream.flush();

        ctest << "chiSquare = " << chiSquare << " )\n}\n";
      }
      ++calPoint;
    }
  */
/*
    if (m_Flags.searchMethod == GRIDSEARCHWITHPUNCH)
    {
        unsigned dataPointSize = Rdouble::withRange().size();
        //ctest << "Parameters: ( ";
        if (m_Flags.punchSymbols.size()){
          for (unsigned varID(0); varID < dataPointSize; ++varID){
            punchStream << "Para" << varID << "\t";
          }
          punchStream << "Temperature (K)\tNumber density\t";
          punchStream << m_Flags.punchSymbols;
          m_Flags.punchSymbols.clear();
        }
        for (unsigned varID(0); varID < dataPointSize; ++varID){
          //ctest << gridArray[i][varID] << " ";
          //punchStream << gridArray[i][varID] << "\t";
          punchStream << *Rdouble::withRange()[varID] << "\t";
        }
        punchStream << m_Env.beta << "\t" << m_Env.conc << "\t";
        punchStream << m_Flags.punchNumbers;
        m_Flags.punchNumbers.clear();

        punchStream.flush();

    }
*/
