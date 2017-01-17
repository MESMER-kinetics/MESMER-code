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
#include "error.h"

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

		// Broadcast strings across processes from a given route.

		void broadcastString(string &tmp, int root) {
			int strSize = tmp.size() + 1;
			MPI_Bcast(&strSize, 1, MPI_LONG, root, MPI_COMM_WORLD);
			char *ctmp = new char[strSize];
			if (m_rank == root)
				strncpy(ctmp, tmp.c_str(), strSize-1);
			ctmp[strSize - 1] = '\0';
			MPI_Bcast(ctmp, strSize, MPI_CHAR, root, MPI_COMM_WORLD);
			tmp = string(ctmp);
			delete[] ctmp;
		}

		// Broadcast vector of strings across processes from a given route.

		void broadcastVecString(vector<string> &tmp, int root) {
			int size = tmp.size();
			MPI_Bcast(&size, 1, MPI_LONG, root, MPI_COMM_WORLD);
			vector<string> svtmp(size);
			if (m_rank == root)
				svtmp = tmp;
			for (int i(0); i < size; i++) {
				string tmp = svtmp[i];
				broadcastString(tmp, root);
				svtmp[i] = tmp;
			}
			tmp = svtmp;
		}

		// Broadcast integers across processes from a given route.

    void broadcastInteger(int *tmp, int size, int root) {
      MPI_Bcast(tmp, size, MPI_LONG, root, MPI_COMM_WORLD);
    }

    // Broadcast doubles across processes from a given route.

    void broadcastDouble(double *tmp, int size, int root) {
      MPI_Bcast(tmp, size, MPI_DOUBLE, root , MPI_COMM_WORLD);
    }

		// Broadcast vector of doubles across processes from a given route.

		void broadcastVecDouble(vector<double> &tmp, int root) {
			int size = tmp.size();
			MPI_Bcast(&size, 1, MPI_LONG, root, MPI_COMM_WORLD);
			vector<double> dtmp(size);
			if (m_rank == root)
				dtmp = tmp ;
			MPI_Bcast(&dtmp[0], size, MPI_DOUBLE, root, MPI_COMM_WORLD);
			tmp = dtmp;
		}

		// Brings processes to the same point.

    void barrier() {

      MPI_Barrier(MPI_COMM_WORLD);

    }

#else

    // Constructor

    ParallelManager(int argc, char ** argv) : m_rank(0), m_size(1) {
			// No Op.
		};

    // Destructor

		~ParallelManager() {} ;

    // Sum vectors across all across processes and redistribute.

    void sumDouble(double *sum, int size) {
      // No Op.
    }

		// Broadcast strings across processes from a given route.

		void broadcastString(char *tmp, int size, int root) {
			// No Op.
		}

		// Broadcast vector of strings across processes from a given route.

		void broadcastVecString(vector<string> &tmp, int root) {
			// No Op.
		}

		// Broadcast integers across processes from a given route.

    void broadcastInteger(int *tmp, int size, int root) {
      // No Op.
    }

    // Broadcast doubles across processes from a given route.

    void broadcastDouble(double *sum, int size, int root) {
      // No Op.
    }

		// Broadcast vector of doubles across processes from a given route.

		void broadcastVecDouble(vector<double> &tmp, int root) {
			// No Op.
		}

		// Brings processes to the same point.

    void barrier() {
      // No Op.
    }

#endif 

    // Accessors 

    int rank() { return m_rank; };
    int size() { return m_size; };
		bool Master() { return (m_rank == 0) ; };

  private:

    int m_rank;
    int m_size;

  };

} //namespace mesmer


#endif // GUARD_ParallelManager_h
