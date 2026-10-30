// main.cpp
//
// Tests for newfemus::MyVector
// - If N < 20  : run verbose correctness tests (with output)
// - If N >= 20 : run timing-only test for localize()
//
// Compile, e.g.:
//   mpic++ -std=c++11 main.cpp -o test_myvector
// Run, e.g.:
//   mpirun -np 4 ./test_myvector 10
//   mpirun -np 4 ./test_myvector 5000000

#include <mpi.h>
#include <iostream>
#include <vector>
#include <cstdlib>  // std::atoi

#include "./NewMyVector.hpp"

using newfemus::MyVector;

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);

  int iproc = 0, nprocs = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &iproc);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  // --------------------------------------
  // Read N (global size) from command line or stdin
  // --------------------------------------
  unsigned N = 0;

  if (argc > 1) {
    N = static_cast<unsigned>(std::atoi(argv[1]));
  }
  else {
    if (iproc == 0) {
      std::cout << "Enter global vector size N: ";
      std::cin >> N;
    }
    // Broadcast N to all ranks
    MPI_Bcast(&N, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
  }

  if (N == 0) {
    if (iproc == 0) {
      std::cout << "N = 0, nothing to do. Exiting.\n";
    }
    MPI_Finalize();
    return 0;
  }

  // ===========================
  // Small N: correctness tests
  // ===========================
  if (N < 20) {
    if (iproc == 0) {
      std::cout << "Running VERBOSE correctness tests with N = "
                << N << " on " << nprocs << " MPI processes\n";
    }

    // 1) SERIAL constructor + fill
    if (iproc == 0) std::cout << "\n[TEST 1] Serial constructor + fill\n";

    MyVector<double> v(N, 0.0);  // serial: all ranks have full vector

    // Fill with pattern v[i] = i
    for (unsigned i = v.begin(); i < v.end(); ++i) {
      // v[i] = static_cast<double>(iproc * N + i);
      v[i] = static_cast<double>(i);
    }

    // Each rank prints what it sees for first and last entries
    for (int kp = 0; kp < nprocs; kp++) {
      if (iproc == kp) {
        std::cout << "Rank " << iproc << " (serial) v[0] = "
                  << v[0] << ", v[N-1] = " << v[N - 1] << std::endl << std::flush;
      }
      MPI_Barrier(MPI_COMM_WORLD);
    }

    // 2) scatter() to PARALLEL layout
    if (iproc == 0) std::cout << "\n[TEST 2] scatter() to parallel layout\n";

    v.scatter();  // now distributed
    //v.stack();  // now distributed

    unsigned b = v.begin();
    unsigned e = v.end();
    unsigned locsize = v.size();

    std::cout << "Rank " << iproc << " local range after scatter: ["
              << b << ", " << e << ") size = " << locsize << std::endl;

    // Check local values: v[i] should still equal i
    for (unsigned gi = b; gi < e; ++gi) {
      double val = v[gi];
      if (val != static_cast<double>(gi)) {
        std::cout << "Rank " << iproc
                  << " ERROR: v[" << gi << "] = " << val
                  << " != " << gi << std::endl;
      }
    }

    MPI_Barrier(MPI_COMM_WORLD);

    // 3) localize() global reconstruction
    if (iproc == 0) std::cout << "\n[TEST 3] localize() to std::vector\n";

    std::vector<double> global_vec;
    v.localize(global_vec);

    if (iproc == 0) {
      std::cout << "Global vector after localize():\n";
      for (unsigned i = 0; i < global_vec.size(); ++i) {
        std::cout << "  g[" << i << "] = " << global_vec[i] << "\n";
      }
    }

    MPI_Barrier(MPI_COMM_WORLD);

    // 4) broadcast / clearBroadcast
    if (iproc == 0) std::cout << "\n[TEST 4] broadcast() + clearBroadcast()\n";

    // Refill distributed vector with rank-dependent pattern: v[i] = 100*rank + global_index
    for (unsigned gi = v.begin(); gi < v.end(); ++gi) {
      v[gi] = 100.0 * iproc + gi;
    }

    // Let rank 0 override its local chunk with a special pattern: 1000 + i
    if (iproc == 0) {
      for (unsigned gi = v.begin(); gi < v.end(); ++gi) {
        v[gi] = 1000.0 + gi;
      }
    }

    // Broadcast from rank 0
    v.broadcast(0);

    // Now every rank should see rank-0’s pattern on its local range
    std::cout << "Rank " << iproc << " after broadcast, local segment:\n";
    for (unsigned gi = v.begin(); gi < v.end(); ++gi) {
      std::cout << "  v[" << gi << "] = " << v[gi] << "\n";
    }

    MPI_Barrier(MPI_COMM_WORLD);

    // Clear broadcast: non-root ranks revert, root keeps its pattern
    v.clearBroadcast();

    std::cout << "Rank " << iproc << " after clearBroadcast, local segment:\n";
    for (unsigned gi = v.begin(); gi < v.end(); ++gi) {
      std::cout << "  v[" << gi << "] = " << v[gi] << "\n";
    }

    MPI_Barrier(MPI_COMM_WORLD);

    if (iproc == 0) {
      std::cout << "\nVERBOSE small-N tests completed.\n";
    }

    MPI_Finalize();
    return 0;
  }

  // ==================================
  // Large N: timing-only for localize
  // ==================================
  const int N_REPS = 10;  // number of repetitions for timing

  if (iproc == 0) {
    std::cout << "Running TIMING test for localize() with N = "
              << N << " on " << nprocs << " MPI processes\n";
  }

  // 1) SERIAL constructor + fill
  MyVector<double> v(N, 0.0);

  for (unsigned i = v.begin(); i < v.end(); ++i) {
    v[i] = static_cast<double>(i);
  }

  // 2) scatter() to PARALLEL layout
  v.scatter();

  unsigned local_size = v.size();
  if (iproc == 0) {
    std::cout << "After scatter: local_size on rank 0 = "
              << local_size << " (other ranks similar)\n";
  }

  MPI_Barrier(MPI_COMM_WORLD);

  // 3) timing localize()
  std::vector<double> global_vec;

  // warm-up call (not timed)
  v.localize(global_vec);
  MPI_Barrier(MPI_COMM_WORLD);

  double t_start = MPI_Wtime();
  for (int rep = 0; rep < N_REPS; ++rep) {
    v.localize(global_vec);
  }
  double t_end = MPI_Wtime();

  double local_elapsed = t_end - t_start;
  double max_elapsed = 0.0;

  MPI_Reduce(&local_elapsed, &max_elapsed, 1, MPI_DOUBLE,
             MPI_MAX, 0, MPI_COMM_WORLD);

  if (iproc == 0) {
    std::cout << "\nlocalize() called " << N_REPS << " times\n"
              << "Total time (max over ranks): " << max_elapsed << " s\n"
              << "Average per call: " << (max_elapsed / N_REPS) << " s\n";
  }

  MPI_Finalize();
  return 0;
}

