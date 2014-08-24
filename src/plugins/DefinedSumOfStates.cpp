//-------------------------------------------------------------------------------------------
//
// DefinedSumOfStates.cpp
//
// Author: Struan Robertson
// Date:   24/Mar/2012
//
// This class implements the defined sum of states model. This class exists to allow users to
// read in their own transition states sum of states, e.g. from FTST calculations or similar.
//
//-------------------------------------------------------------------------------------------
#include <vector>
#include <string>
#include "../System.h"
#include "../Spline.h"

using namespace std;
using namespace Constants;

namespace mesmer
{
  class DefinedSumOfStates : public MicroRateCalculator
  {
  public:

    // Constructor which registers with the list of MicroRateCalculators in the base class.
    DefinedSumOfStates(const char* id) :m_id(id), m_logSpline(true), m_WE(), m_E() { Register(); }

    virtual ~DefinedSumOfStates() {}
    virtual const char* getID()  { return m_id; }
    virtual DefinedSumOfStates* Clone() { return new DefinedSumOfStates(*this); }

    virtual bool calculateMicroCnlFlux(Reaction* pReac) ;

    virtual bool ReadParameters(Reaction* pReac)  { return true ; }
    virtual bool ParseData(PersistPtr pp);

  private:

    typedef vector<pair<double,double> > JData ;
    typedef map<double, JData* > SumOfStatesData ;
    typedef SumOfStatesData::const_iterator Citr ;
    typedef SumOfStatesData::iterator Itr ;

    // This method reads in an energy (only) dependent sum of states.
    bool readEneDepSOS(PersistPtr pp) ;

    // This method reads in an energy and J dependent sum of states, and interpolates over J.
    bool readEneAndJDepSOS(PersistPtr pp) ;

  private:
    const char* m_id;
    bool m_logSpline;
    vector<double> m_WE, m_E ;
  };

  //************************************************************
  DefinedSumOfStates theDefinedSumOfStates("DefinedSumOfStates");
  //************************************************************

  bool DefinedSumOfStates::ParseData(PersistPtr pp)
  {
    Molecule* pTransitionState = m_parent->get_TransitionState();

    PersistPtr ppSos = pp->XmlMoveTo("me:SumOfStates") ;
    if (!ppSos)
      ppSos = (pTransitionState->get_PersistentPointer())->XmlMoveTo("me:SumOfStates");
    if (!ppSos)
      throw(std::runtime_error("No sum of states defined in either "
      + m_parent->getName() + " or " +pTransitionState->getName()))  ;

    bool m_logSpline = !(ppSos->XmlReadBoolean("nologSpline")) ;

    // Initialize arrays, of possibly J interpolated, sums of states.
    m_E.clear()  ; m_E.push_back(0.0) ;
    m_WE.clear() ; m_WE.push_back(1.0) ;

    // Read sum of states data from XML representation, interpolating J values if required.

    if (ppSos->XmlReadBoolean("angularMomentum")) {
      readEneAndJDepSOS(ppSos) ;
    } else {
      readEneDepSOS(ppSos) ;
    }

    // Apply symmetry number.

    double symmetryNumber = pTransitionState->getDOS().get_Sym() ;
    for (size_t i(0) ; i < m_WE.size() ; i++) {
      m_WE[i] /= symmetryNumber ;
      if (m_logSpline) 
        m_WE[i] = log(m_WE[i]) ;
    }
    return true;
  }
  //
  // This method calculates the reaction flux. 
  //

  bool DefinedSumOfStates::calculateMicroCnlFlux(Reaction* pReact)
  {
    // Allocate space to hold transition state flux and initialize elements to zero.
    vector<double>& rxnFlux = pReact->get_CellFlux();
    rxnFlux.clear();
    const size_t MaximumCell = pReact->getEnv().MaxCell;
    rxnFlux.resize(MaximumCell, 0.0);

    // Spline the sum of states.

    Spline spline ;
    spline.Initialize(m_E, m_WE) ;

    double cellSize  = pReact->getEnv().CellSize ;
    double maxEnergy = m_E[m_E.size()-1] ;

    // Calculate the cell reaction flux from the spline.

    for (size_t i(0) ; i < MaximumCell ; ++i ) {
      double SumOfStates = spline.Calculate(min(double(i)*cellSize,maxEnergy)) ;
      if (m_logSpline) SumOfStates = exp(SumOfStates) ;
      rxnFlux[i] = SumOfStates * SpeedOfLight_in_cm;
    }

    // Set the location of the transition state relative to the reactants.
    pReact->setCellFluxBottom(pReact->get_relative_rctZPE() + pReact->get_ThresholdEnergy());

    return true;
  }

  // This method reads in an energy (only) dependent sum of states.
  bool DefinedSumOfStates::readEneDepSOS(PersistPtr pp)
  {
    const char* txt = pp->XmlReadValue("units", optional);
    string units = txt ? txt : "kJ/mol";

    while (pp = pp->XmlMoveTo("me:SumOfStatesPoint")) {
      double SumOfStates(0.0) ; 
      const char* txt = pp->XmlRead();
      stringstream s1(txt); s1 >> SumOfStates;
      m_WE.push_back(SumOfStates) ;
      double energy = pp->XmlReadDouble("energy", false);
      energy = getConvertedEnergy(units, energy);
      m_E.push_back(energy) ;
    }

    return true;
  } 

  // This method reads in an energy and J dependent sum of states, and interpolates over J.
  bool DefinedSumOfStates::readEneAndJDepSOS(PersistPtr pp)
  {

    const char* txt = pp->XmlReadValue("units", optional);
    string units = txt ? txt : "kJ/mol";

    int JMin(0), JMax(0) ;
    SumOfStatesData WEJ;
    while (pp = pp->XmlMoveTo("me:SumOfStatesPoint")) {

      double SumOfStates(0.0) ; 
      const char* txt = pp->XmlRead();
      stringstream s1(txt); s1 >> SumOfStates;
      double energy = pp->XmlReadDouble("energy", false);
      energy = getConvertedEnergy(units, energy);
      double angMomMag = pp->XmlReadDouble("angMomMag", false);

      JMax = max(JMax, int(angMomMag)) ;

      // Test to see if a new value of Energy has been found. 

      if (IsNan(angMomMag)) angMomMag = 0 ;
      JData* pJData = WEJ[energy] ;

      if (!pJData) {
        pJData = new JData ;
        WEJ[energy] = pJData ;
      }

      SumOfStates = (m_logSpline) ? log(max(SumOfStates,1.0)) : SumOfStates ;
      pJData->push_back(pair<double,double>(angMomMag, SumOfStates)) ;
    }

    // Use a cubic spline to interpolate and sum over the J dimension.

    Citr citr = WEJ.begin() ;
    for (size_t i(1) ; citr != WEJ.end() ; citr++, i++ ) {
      Spline spline ;
      spline.Initialize(citr->second) ;

      double sum(0.0) ;
      for (int j(JMin) ; j < JMax ; j++) {
        double tmp = spline.Calculate(double(j)) ;
        sum += (m_logSpline) ? exp(tmp) : tmp ;
      }
      m_WE.push_back(sum) ;
      m_E.push_back(citr->first) ;
    }

    // Delete temporary data representation.
    Itr itr = WEJ.begin() ;
    for (; itr != WEJ.end() ; itr++) {
      itr->second->clear() ;
      delete itr->second ;
    }

    return true;
  }

} //namespace
