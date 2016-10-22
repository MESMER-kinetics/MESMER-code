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

			// Sum vectors across all across processes and redistribute.

			void sumDouble(double *sum, int size) {
				
				vector<double> tmp(size, 0.0);

				MPI_Allreduce(sum, &tmp[0], size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

				for (size_t i(0) ; i < size_t(size) ; i++) {
					sum[i] = tmp[i];
				}

			}

			// Broadcast doubles across processes from a given route.

			void broadcastDouble(double *tmp, int size, int root) {

				MPI_Bcast(tmp, size, MPI_DOUBLE, root , MPI_COMM_WORLD);

			}

			// Brings processes to the same point.

			void barrier() {

				MPI_Barrier(MPI_COMM_WORLD);

			}

#else

			// Constructor

			ParallelManager(int argc, char ** argv) : m_rank(0), m_size(1) {} ;

			// Destructor

			~ParallelManager() 

			// Sum vectors across all across processes and redistribute.

			void sumDouble(double *sum, int size) {
				// No Op.
			}

			// Broadcast doubles across processes from a given route.

			void broadcastDouble(double *sum, int size, int root) {
				// No Op.
      }

			// Brings processes to the same point.

			void barrier() {
				// No Op.
			}

#endif 

			// Accessors 

			int rank() { return m_rank ; } ;
			int size() { return m_size ; };

	private:

		int m_rank;
		int m_size;

	};

} //namespace mesmer


#endif // GUARD_ParallelManager_h
