//-------------------------------------------------------------------------------------------
//
// PriorDistFragmentation.h
//
// Author: Struan Robertson
// Date:   8th October 2019
//
// This file contains the definition of the PriorDistFragmentation and related implementation
// classes. 
//
//-------------------------------------------------------------------------------------------
#ifndef GUARD_PriorDistFragmentation_h
#define GUARD_PriorDistFragmentation_h

#include <vector>
#include "../Reaction.h"
#include "../Fragmentation.h"

namespace mesmer
{
  //
  // Implementation class for the calculation of fragment distribution on dissocistion
  // for the prior model.
  // SHR 16/Jun/2013: This class definition should probably be a plug-in class.
  //
  class priorDist : public FragDist
  {
  public:

    // Constructors.
    priorDist() : m_pReaction(NULL),
     m_rctDOS(),
     m_upperConv(),
     m_lowerConv()
    {};

    // Destructor.
    virtual ~priorDist() {};

    // Initialize the fragment distribution.
    virtual void initialize(Reaction* pReaction);

    // Calculate distribution
    virtual void calculate(double excessEnergy, std::vector<double>& dist);

    // Return resources
    virtual void clear() {
      m_rctDOS.clear();
      m_upperConv.clear();
      m_lowerConv.clear();
    };

  protected:

    Reaction* m_pReaction;

    vector<double> m_rctDOS;

    vector<double> m_upperConv;

    vector<double> m_lowerConv;

  };

  class modPriorDist : public priorDist
  {
  public:

    // Constructors.
    modPriorDist(PersistPtr ppFragDist, std::string name) : priorDist() {
      m_order = ppFragDist->XmlReadDouble("me:modPriorOrder");
      m_nexp = ppFragDist->XmlReadDouble("me:modPriorNexp");
      m_Tref = ppFragDist->XmlReadDouble("me:modPriorTref");
      bool rangeSet(false);
      PersistPtr ppOrder = ppFragDist->XmlMoveTo("me:modPriorOrder");
      ReadRdoubleRange(name + std::string(":modPriorOrder"), ppOrder, m_order, rangeSet);
      PersistPtr ppNexp = ppFragDist->XmlMoveTo("me:modPriorNexp");
      ReadRdoubleRange(name + std::string(":modPriorNexp"), ppNexp, m_nexp, rangeSet);
    };

    // Destructor.
    virtual ~modPriorDist() {};

    // Initialize the fragment distribution.
    virtual void initialize(Reaction* pReaction);

  private:

    Rdouble m_order;
    Rdouble m_nexp;
    double m_Tref;

  };


}//namespace

#endif

