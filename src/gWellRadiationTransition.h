#ifndef GUARD_gWellRadiationTransition_h
#define GUARD_gWellRadiationTransition_h

//-------------------------------------------------------------------------------------------
//
// gWellRadiationTransition.h
//
// Author: Struan Robertson
// Date:   27/Dec/2021
//
// This header file contains the declaration of the gWellRadiationTransition class. This 
// class governs radition transitions betwen states.
//
//-------------------------------------------------------------------------------------------

#include "MolecularComponents.h"
#include <vector> 
#include "gDensityOfStates.h"

namespace mesmer
{

  class gWellRadiationTransition :public MolecularComponent
  {
    //-------------------------------------------------------------------------------------------------
    // Radiation transfer related properties
    //-------------------------------------------------------------------------------------------------

  public:

    //
    // Constructor, destructor and initialization
    //
    gWellRadiationTransition(Molecule* pMol);
    virtual ~gWellRadiationTransition();
    bool initialization();

    // Set transition frequencies.
    void setTransitionFrequencies(vector<double>& frequencies) {
      m_TransitionFrequency = frequencies;
    }

    // Calculate radiation operator.
    template<class T>
    bool RadiationOperator(MesmerEnv& env, TMatrix<T>** egme) const;

    // Add radiation operator to overall transtion matrix.
    template<class T>
    void addRadiationOperator(qdMatrix* CollOptr, TMatrix<T>* egme, const size_t locate, const double RducdOmega) const;

  private:

    // Returns the Planck distribution for a frequency expressed in wavenumbers. Returned units are Js/m3.
    double planckDistribtuion(double ene, double beta) const {
      return 8.0 * M_PI * pow(ene, 3.0) * 1.0e+06 * PlancksConstant_in_JouleSecond / (exp(beta * ene) - 1.0);
    }

    const size_t reservoirShift() const { return m_numGroupedGrains == 0 ? 0 : m_numGroupedGrains - 1; }

    const double getLowestBarrier() { return m_lowestBarrier; }

    void setLowestBarrier(double value) { m_lowestBarrier = value; }

    std::vector<double> m_TransitionFrequency; // Transition freqeuncies (in most cases these are same as normal mode freqencies).
    std::vector<double> m_EinsteinBij;         // The associated Einstein Bij freqnecies. 

    size_t m_numGroupedGrains;    // Number of grains grouped into a reservoir grain.

    double m_lowestBarrier;       // lowest barrier associatied with this species
  };

  //
  // Calculate collision operator
  //
  template<class T>
  bool gWellRadiationTransition::RadiationOperator(MesmerEnv& env, TMatrix<T>** CollOp) const {

    // Are there any Einstein coeffcients?
    if (m_EinsteinBij.size() == 0)
      return false;

    vector<double> gEne;
    vector<double> gDOS;
    m_host->getDOS().getGrainEnergies(gEne);
    m_host->getDOS().getGrainDensityOfStates(gDOS);

    const size_t nradoptrsize = (*CollOp)->size() ;

    // Allocate memory.
    TMatrix<T>* transitionMatrix = new TMatrix<T>(nradoptrsize, T(0.0));

    const T beta = T(env.beta);

    // Determine the upward radiative transitions first, by filling sub diagonal elements first.
    for (size_t idx(0); idx < m_EinsteinBij.size(); idx++) {
      
      const T excitationRate = T(m_EinsteinBij[idx] * planckDistribtuion(m_TransitionFrequency[idx], env.beta));

      // Search for first transition.  
      size_t ll(0);
      for (; m_TransitionFrequency[idx] > gEne[ll]; ll++) ;

      for (size_t i = ll, j = 0; i < nradoptrsize; ++i, ++j) {
        T ei = T(gEne[i]);
        T ni = T(gDOS[i]);
        T ej = T(gEne[j]);
        T nj = T(gDOS[j]);

        // Transfer to higher energy (adsorption).
        (*transitionMatrix)[i][j] += excitationRate;

        // Transfer to lower energy (stimulated and spontaneous emission via detailed balance).
        (*transitionMatrix)[j][i] += excitationRate * (nj / ni) * exp(-beta * (ej - ei));
      }
    }

    // Mass conservation.

    for (size_t i = 0; i < nradoptrsize; ++i) {
      (*transitionMatrix)[i][i] = T(0.0);
      T sum = T(0.0);
      for (size_t j = 0; j < nradoptrsize; ++j) {
        sum += (*transitionMatrix)[j][i];
      }
      (*transitionMatrix)[i][i] = -sum;
    }

    string title("Radiation Transition Matrix:");
    transitionMatrix->print(title, stest, 20, 20, -1, -1);

    // Account for a reservoir state.

    vector<T> popDist; // Grained population distribution.
    popDist.push_back(0.0);
    T prtnFn(0.0);
    for (size_t idx(0); idx < gDOS.size(); ++idx) {
      const T tmp(T(gDOS[idx]) * exp(-beta * T(gEne[idx])));
      prtnFn += tmp;
      if (idx < std::max(m_numGroupedGrains, size_t(1))) {
        popDist[0] += tmp;
      }
      else {
        popDist.push_back(tmp);
      }
    }

    const size_t reducedCollOptrSize = nradoptrsize - reservoirShift();

    // Symmetrization of the radiation matrix.
    for (size_t i(1); i < reducedCollOptrSize; ++i) {
      for (size_t j(0); j < i; ++j) {
        (*transitionMatrix)[j][i] *= sqrt(popDist[i] / popDist[j]);
        (*transitionMatrix)[i][j]  = (*transitionMatrix)[j][i];
      }
    }

    if (*CollOp)
      delete* CollOp;  // Delete any existing matrix.
    (*CollOp) = new TMatrix<T>(reducedCollOptrSize);
    (**CollOp) = (*transitionMatrix);

    delete transitionMatrix;

    return true;
  }

  // Add radiation operator to overall transtion matrix.
  template<class T>
  void gWellRadiationTransition::addRadiationOperator(qdMatrix* CollOptr, TMatrix<T>* egme, const size_t locate, const double RducdOmega) const {

    // Find size of system matrix.

    const size_t smsize = CollOptr->size();
    const size_t nradoptrsize = egme->size();

    // Check there is enough space in system matrix.

    if (locate + nradoptrsize > smsize)
      throw (std::runtime_error("Error in the size of the system matrix."));

    // Add radiation operator to the diagonal block indicated by "locate"
    // and multiply by the reduced collision frequency.

    for (size_t i(0), ii(locate); i < nradoptrsize; ++i, ++ii) {
      for (size_t j(0), jj(locate); j < nradoptrsize; ++j, ++jj) {
        (*CollOptr)[ii][jj] += RducdOmega * (*egme)[i][j];
      }
    }

  }

}

#endif // GUARD_gWellRadiationTransition_h
