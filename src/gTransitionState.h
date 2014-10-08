#ifndef GUARD_gTransitionState_h
#define GUARD_gTransitionState_h

#include "MolecularComponents.h"

namespace mesmer
{

  class gTransitionState:public MolecularComponent
  {

    //-------------------------------------------------------------------------------------------------
    // Transition state related properties
    //-------------------------------------------------------------------------------------------------

  private:
    Rdouble m_ImFreq;            // Imaginary frequency of this barrier (For tunneling in QM calculations)
    //================================================
    // CHECK FOR INPUTFILE PARAMETERS
    int m_ImFreq_chk;
    //================================================

  public:
    //
    // Constructor, destructor and initialization
    //
    gTransitionState(Molecule* pMol);
    virtual ~gTransitionState();

    double get_ImFreq();
    void set_imFreq(const double value){ m_ImFreq = value; m_ImFreq_chk = 0;}

  };

}//namespace
#endif //GUARD_gTransitionState_h
