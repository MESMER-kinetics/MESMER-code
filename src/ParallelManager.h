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

    ParallelManager(int argc, char** argv) : m_rank(0), m_size(1), m_mutex(false) {

      MPI_Init(&argc, &argv);

      MPI_Comm_rank(MPI_COMM_WORLD, &m_rank);
      MPI_Comm_size(MPI_COMM_WORLD, &m_size);

    };

    // Destructor

    ~ParallelManager() {};

    // Finalize

    void finalize() { MPI_Finalize(); };

    // Sum vectors across all across processes and redistribute.

    void sumDouble(double* sum, int size) {

      vector<double> tmp(size, 0.0);

      MPI_Allreduce(sum, &tmp[0], size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

      for (size_t i(0); i < size_t(size); i++) {
        sum[i] = tmp[i];
      }

    }

    // Broadcast strings across processes from a given root.

    void broadcastString(string& tmp, int root) {
      int strSize(0), intSize(1);
      if (m_rank == root)
        strSize = tmp.size() + 1;
      MPI_Bcast(&strSize, intSize, MPI_INT, root, MPI_COMM_WORLD);
      char* ctmp = new char[strSize];
      if (m_rank == root)
        std::strcpy(ctmp, tmp.c_str());
      MPI_Bcast(ctmp, strSize, MPI_CHAR, root, MPI_COMM_WORLD);
      tmp = string(ctmp);
      delete[] ctmp;
    }

    // Broadcast vector of strings across processes from a given root.

    void broadcastVecString(vector<string>& tmp, int root) {
      int size = tmp.size();
      MPI_Bcast(&size, 1, MPI_INT, root, MPI_COMM_WORLD);
      if (size == 0) // Nothing to copy so return;
        return;
      if (m_rank != root)
        tmp.clear();
      for (int i(0); i < size; i++) {
        string stmp(" ");
        if (m_rank == root)
          stmp = tmp[i];
        broadcastString(stmp, root);
        if (m_rank != root)
          tmp.push_back(stmp);
      }
    }

    // Broadcast integers across processes from a given root.

    void broadcastInteger(int* tmp, int size, int root) {
      MPI_Bcast(tmp, size, MPI_INT, root, MPI_COMM_WORLD);
    }

    // Broadcast doubles across processes from a given root.

    void broadcastDouble(double* tmp, int size, int root) {
      MPI_Bcast(tmp, size, MPI_DOUBLE, root, MPI_COMM_WORLD);
    }

    // Broadcast vector of doubles across processes from a given route.

    void broadcastVecDouble(vector<double>& tmp, int root) {
      int size = tmp.size();
      MPI_Bcast(&size, 1, MPI_INT, root, MPI_COMM_WORLD);
      vector<double> dtmp(size);
      if (m_rank == root)
        dtmp = tmp;
      MPI_Bcast(&dtmp[0], size, MPI_DOUBLE, root, MPI_COMM_WORLD);
      tmp = dtmp;
    }

    // Broadcast throw across processes from an unknown rank.

    void CheckForThrow(int rank, string error_msg) {

      // First find if a  rank has thrown.
      int root;
      MPI_Allreduce(&rank, &root, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

      // If one or more ranks have thrown, combine all error messages and re-throw.
      if (root >= 0) {
        stringstream ss;
        for (int i(0); i < m_size; i++) {
          string str = error_msg;
          broadcastString(str, i);
          ss << str;
        }
        throw(runtime_error(ss.str()));
      }
    }

    // Brings processes to the same point.

    void barrier() {
      MPI_Barrier(MPI_COMM_WORLD);
    }

    // Gathers QA output from across ranks and writes to ctest.

    void writeCtest() {
      ctest << stest.str();
      for (int i(1); i < m_size; i++) {
        string tmp = stest.str();
        broadcastString(tmp, i);
        ctest << tmp;
      }
      stest.str("");
    }

#else

    // Constructor

    ParallelManager(int argc, char** argv) : m_rank(0), m_size(1), m_mutex(false) {
      // No Op.
    };

    // Destructor

    ~ParallelManager() {};

    // Finalize

    void finalize() {};

    // Sum vectors across all across processes and redistribute.

    void sumDouble(double* sum, int size) {
      // No Op.
    }

    // Broadcast strings across processes from a given route.

    void broadcastString(string& tmp, int root) {
      // No Op.
    }

    // Broadcast vector of strings across processes from a given route.

    void broadcastVecString(vector<string>& tmp, int root) {
      // No Op.
    }

    // Broadcast integers across processes from a given route.

    void broadcastInteger(int* tmp, int size, int root) {
      // No Op.
    }

    // Broadcast doubles across processes from a given route.

    void broadcastDouble(double* sum, int size, int root) {
      // No Op.
    }

    // Broadcast vector of doubles across processes from a given route.

    void broadcastVecDouble(vector<double>& tmp, int root) {
      // No Op.
    }

    // Broadcast throw across processes from an unknown rank.

    void CheckForThrow(int rank, string error_msg) {
      if (rank >= 0) {
        throw(runtime_error(error_msg));
      }
    }

    // Brings processes to the same point.

    void barrier() {
      // No Op.
    }

    // Gathers QA output from across ranks ans writes to ctest.

    void writeCtest() {
      ctest << stest.str();
      stest.str("");
    }

#endif 

    // Accessors 

    int rank() { return m_rank; };
    int size() { return m_size; };
    bool Master() { return (m_rank == 0); };
    void setMutex(const char* function) {
      if (!m_mutex) {
        m_mutex = true;
      }
      else {
        string msg = string(function) + string(": Parellel execution is recursive.\n");
        throw(runtime_error(msg));
      }
    };
    void clearMutex() { m_mutex = false; };

  private:

    int m_rank;
    int m_size;
    bool m_mutex;

  };

} //namespace mesmer


#endif // GUARD_ParallelManager_h
