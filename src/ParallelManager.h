#ifndef GUARD_ParallelManager_h
#define GUARD_ParallelManager_h

//-------------------------------------------------------------------------------------------
//
// ParallelManager.h
//
// Author: Struan Robertson
// Date:   27/Aug/2016
//
// This header file contains the Mesmer parallel implementation based on the MPI substrate.
//
//-------------------------------------------------------------------------------------------


#ifdef PARALLEL
#include "mpi.h"
#endif

namespace mesmer
{

	class ParallelManager {

	  public: 

#ifdef PARALLEL

			// Constructor

			ParallelManager(int argc, char ** argv) : m_rank(0), m_size(1) {

				MPI_Init(&argc, &argv); 

				MPI_Comm_rank(MPI_COMM_WORLD, &m_rank);
				MPI_Comm_size(MPI_COMM_WORLD, &m_size);

				//  printf("I am %d of %d\n", m_rank + 1, m_size);
			} ;

			// Destructor

			~ParallelManager() { MPI_Finalize(); } ;

#else

			// Constructor

			ParallelManager(int argc, char ** argv) : m_rank(0), m_size(1) {} ;

			// Destructor

			~ParallelManager() 

#endif 

			int rank() { return m_rank ; } ;
			int size() { return m_size ; };

	private:

		int m_rank;
		int m_size;

	};

} //namespace mesmer


#endif // GUARD_ParallelManager_h
