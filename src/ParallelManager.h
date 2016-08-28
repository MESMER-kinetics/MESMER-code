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

			// Calculate loop indicies to distribute work across processes.

			void loopDistribution(size_t loopSize, size_t &lidx, size_t &uidx) {
				size_t loopSegment   = loopSize / m_size ;
				size_t loopRemainder = loopSize % m_size ;

				if (m_rank < int(loopRemainder)) {
					lidx =  m_rank*(loopSegment + 1) ;
					uidx = (m_rank+1)*(loopSegment + 1) ;
				}
				else {
					lidx = loopRemainder + m_rank*loopSegment ;
					uidx = loopRemainder + (m_rank + 1)*loopSegment;
				}

				printf("I am %d of %d : Loopsize = %d ll = %d ul = %d \n", m_rank, m_size, loopSize, lidx, uidx);

			} ;

			// Sum vectors across all across processes and redistribute.

			void sumDouble(double *sum, double *rSum, int size) {
				MPI_Allreduce(sum, rSum, size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
			}

#else

			// Constructor

			ParallelManager(int argc, char ** argv) : m_rank(0), m_size(1) {} ;

			// Destructor

			~ParallelManager() 

			// Calculate loop indicies to distribute work across processes.

			void loopDistribution(size_t loopSize, size_t &lidx, size_t &uidx) {
				lidx = 0 ;
				uidx = loopSize ;
			};

			// Sum vectors across all across processes and redistribute.

			void sumDouble(double *sum, double * rSum, int size) {
				for (int i(0); i < size, i++) {
					rSum[i] = sum[i];
				}
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
