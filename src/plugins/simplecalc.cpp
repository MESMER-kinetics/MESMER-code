// simplecalc.cpp
// Calculates for a range of pressure/temperature conditions but
// with other model parameters held constant.

#include "../System.h"
#include "../calcmethod.h"

namespace mesmer
{
class SimpleCalc : public CalcMethod
{
public:
  SimpleCalc(const char* id) : m_id(id) { Register(); }
  virtual ~SimpleCalc() {}
  virtual const char* getID()  { return m_id; }

    //Function to do the work
  virtual bool DoCalculation(System* pSys)  ;
private:
  const char* m_id;

};

////////////////////////////////////////////////
//Global instance
SimpleCalc theSimpleCalc("simpleCalc");
///////////////////////////////////////////////

bool SimpleCalc::DoCalculation(System* pSys)
{
  double chiSquare(1000.0);
  return pSys->calculate(chiSquare) ;
}

}//namespace

