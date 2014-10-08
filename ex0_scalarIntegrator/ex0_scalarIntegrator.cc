// -*- C++ -*-
// ex0_scalarIntegrator.cc
// an exercise for the sandia 2014 clinic team.
// here we integrate sin over some interval using various methods

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <array>
#include <string>
#include <chrono>

// header file for openmp
#include <omp.h>

// header files for tbb
#include <tbb/blocked_range.h>
#include <tbb/parallel_reduce.h>
#include <tbb/parallel_for.h>
#include <tbb/task_scheduler_init.h>

// header files for cuda implementation
#include "ex0_scalarIntegrator_cuda.cuh"

// header files for kokkos
#include <Kokkos_Core.hpp>

using std::string;
using std::vector;
using std::array;
using std::chrono::high_resolution_clock;
using std::chrono::duration;
using std::chrono::duration_cast;

class TbbFunctor {
public:
  
  const double _startBound;
  const double _endBound;
  const double _dx;
  double _integral;

  TbbFunctor(const double startBound, const double endBound, 
  	const double dx) : 
	_startBound(startBound), _endBound(endBound),
	_dx(dx), _integral(0) {
  }

  TbbFunctor(const TbbFunctor & other,
             tbb::split) : 
	     _startBound(other._startBound), _endBound(other._endBound),
	     _dx(other._dx), _integral(0) {
  }

  void operator()(const tbb::blocked_range<size_t> & range) {
    double sum = _integral;
    double begin = range.begin()*_dx;
    double end = range.end()*_dx;
    for (double i = begin+.5*_dx; i < end; i += _dx) {
      sum += std::sin(i);
    }
    _integral = sum;
    
  }

  void join(const TbbFunctor & other) {
    _integral += other._integral;
  }

private:
  TbbFunctor();

};

struct KokkosFunctor {

  const double _dx;

  KokkosFunctor(const double dx) : _dx(dx) {
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const unsigned int intervalIndex) const {
  }

private:
  KokkosFunctor();

};

int main(int argc, char* argv[]) {

  // a couple of inputs.  change the numberOfIntervals to control the amount
  //  of work done
  const unsigned long numberOfIntervals = 1e8;
  // the integration bounds
  const array<double, 2> bounds = {{0, 1.314}};

  // these are c++ timers...for timing
  high_resolution_clock::time_point tic;
  high_resolution_clock::time_point toc;

  const double dx = (bounds[1] - bounds[0]) / numberOfIntervals;

  // calculate analytic solution
  const double libraryAnswer =
    std::cos(bounds[0]) - std::cos(bounds[1]);

  // ===============================================================
  // ********************** < do slow serial> **********************
  // vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

  double slowSerialIntegral = 0;
  tic = high_resolution_clock::now();
  for (unsigned int intervalIndex = 0;
       intervalIndex < numberOfIntervals; ++intervalIndex) {
    if (intervalIndex % (numberOfIntervals / 10) == 0) {
      printf("serial calculation on interval %8.2e / %8.2e (%%%5.1f)\n",
             double(intervalIndex),
             double(numberOfIntervals),
             100. * intervalIndex / double(numberOfIntervals));
    }
    const double evaluationPoint =
      bounds[0] + (double(intervalIndex) + 0.5) * dx;
    slowSerialIntegral += std::sin(evaluationPoint);
  }
  slowSerialIntegral *= dx;
  toc = high_resolution_clock::now();
  const double slowSerialElapsedTime =
    duration_cast<duration<double> >(toc - tic).count();

  // check our serial answer
  const double serialRelativeError =
    std::abs(libraryAnswer - slowSerialIntegral) / std::abs(libraryAnswer);
  if (serialRelativeError > 1e-3) {
    fprintf(stderr, "our answer is too far off: %15.8e instead of %15.8e\n",
            slowSerialIntegral, libraryAnswer);
    exit(1);
  }

  // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  // ********************** </do slow serial> **********************
  // ===============================================================

  // ===============================================================
  // ********************** < do fast serial> **********************
  // vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

  const double fastSerialIntegral = slowSerialIntegral;
  const double fastSerialElapsedTime = slowSerialElapsedTime;

  // check our serial answer
  const double fastSerialRelativeError =
    std::abs(libraryAnswer - fastSerialIntegral) / std::abs(libraryAnswer);
  if (fastSerialRelativeError > 1e-3) {
    fprintf(stderr, "our answer is too far off: %15.8e instead of %15.8e\n",
            fastSerialIntegral, libraryAnswer);
    exit(1);
  }

  // output speedup
  printf("fast: time %8.2e speedup %8.2e\n",
         fastSerialElapsedTime,
         slowSerialElapsedTime / fastSerialElapsedTime);

  // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  // ********************** </do fast serial> **********************
  // ===============================================================

  // we will repeat the computation for each of the numbers of threads
  vector<unsigned int> numberOfThreadsArray;
  numberOfThreadsArray.push_back(1);
  numberOfThreadsArray.push_back(2);
  numberOfThreadsArray.push_back(4);
  numberOfThreadsArray.push_back(8);
  numberOfThreadsArray.push_back(16);
  numberOfThreadsArray.push_back(24);

  // ===============================================================
  // ********************** < do tbb> ******************************
  // vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

  printf("performing calculations with tbb\n");
  // for each number of threads
  for (const unsigned int numberOfThreads :
         numberOfThreadsArray) {

    // initialize tbb's threading system for this number of threads
    tbb::task_scheduler_init init(numberOfThreads);

    // prepare the tbb functor.
    TbbFunctor tbbFunctor(bounds[0], bounds[1], dx);

    // start timing
    tic = high_resolution_clock::now();
    // dispatch threads
    parallel_reduce(tbb::blocked_range<size_t>(0, numberOfIntervals),
                    tbbFunctor);
    // stop timing
    toc = high_resolution_clock::now();
    const double threadedElapsedTime =
      duration_cast<duration<double> >(toc - tic).count();

    // somehow get the threaded integral answer
    const double threadedIntegral = tbbFunctor._integral*dx;
    // check the answer
    const double threadedRelativeError =
      std::abs(libraryAnswer - threadedIntegral) / std::abs(libraryAnswer);
    if (threadedRelativeError > 1e-3) {
      fprintf(stderr, "our answer is too far off: %15.8e instead of %15.8e\n",
              threadedIntegral, libraryAnswer);
      exit(1);
    }

    // output speedup
    printf("%3u : time %8.2e speedup %8.2e (%%%5.1f of ideal)\n",
           numberOfThreads,
           threadedElapsedTime,
           fastSerialElapsedTime / threadedElapsedTime,
           100. * fastSerialElapsedTime / threadedElapsedTime / numberOfThreads);
  }

  // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  // ********************** </do tbb> ******************************
  // ===============================================================


  // ===============================================================
  // ********************** < do openmp> ***************************
  // vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

  printf("performing calculations with openmp\n");
  // for each number of threads
  for (const unsigned int numberOfThreads :
         numberOfThreadsArray) {

    // TODO: set the number of threads for openmp
    omp_set_num_threads(numberOfThreads);
    // start timing
    tic = high_resolution_clock::now();
    
    double threadedIntegral = 0;
    // TODO: do scalar integration with threads on openmp
    #pragma omp parallel for reduction(+: threadedIntegral)
    for (unsigned int n = 0; n < numberOfIntervals; n++) {
      threadedIntegral += std::sin(bounds[0]+.5*dx+n*dx);
    }
    threadedIntegral *= dx;
  
    // stop timing
    toc = high_resolution_clock::now();
    const double threadedElapsedTime =
      duration_cast<duration<double> >(toc - tic).count();

    // check the answer
    const double threadedRelativeError =
      std::abs(libraryAnswer - threadedIntegral) / std::abs(libraryAnswer);
    if (threadedRelativeError > 1e-3) {
      fprintf(stderr, "our answer is too far off: %15.8e instead of %15.8e\n",
              threadedIntegral, libraryAnswer);
      exit(1);
    }

    // output speedup
    printf("%3u : time %8.2e speedup %8.2e (%%%5.1f of ideal)\n",
           numberOfThreads,
           threadedElapsedTime,
           fastSerialElapsedTime / threadedElapsedTime,
           100. * fastSerialElapsedTime / threadedElapsedTime / numberOfThreads);
  }

  // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  // ********************** </do openmp> ***************************
  // ===============================================================

  // ===============================================================
  // ********************** < do cuda> *****************************
  // vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

  // we will repeat the computation for each of the numbers of threads
  vector<unsigned int> threadsPerBlockArray;
  threadsPerBlockArray.push_back(32);
  threadsPerBlockArray.push_back(64);
  threadsPerBlockArray.push_back(128);
  threadsPerBlockArray.push_back(256);
  threadsPerBlockArray.push_back(512);

  printf("performing calculations with cuda\n");
  // for each number of threads per block
  for (const unsigned int numberOfThreadsPerBlock :
         threadsPerBlockArray) {

    // start timing
    tic = high_resolution_clock::now();

    double cudaIntegral = 0;
    // TODO: do scalar integration with cuda for this number of threads per block
    cudaDoScalarIntegration(numberOfThreadsPerBlock, bounds[0], bounds[1],
			    dx, &cudaIntegral);

    // stop timing
    toc = high_resolution_clock::now();
    const double cudaElapsedTime =
      duration_cast<duration<double> >(toc - tic).count();

    // check the answer
    const double cudaRelativeError =
      std::abs(libraryAnswer - cudaIntegral) / std::abs(libraryAnswer);
    if (cudaRelativeError > 1e-3) {
      fprintf(stderr, "our answer is too far off: %15.8e instead of %15.8e\n",
              cudaIntegral, libraryAnswer);
    //  exit(1);
    }

    // output speedup
    printf("%3u : time %8.2e speedup %8.2e\n",
           numberOfThreadsPerBlock,
           cudaElapsedTime,
           fastSerialElapsedTime / cudaElapsedTime);
  }

  // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  // ********************** </do cuda> *****************************
  // ===============================================================


  // ===============================================================
  // ********************** < do kokkos> *****************************
  // vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

  printf("performing calculations with kokkos on %s\n",
         typeid(Kokkos::DefaultExecutionSpace).name());

  Kokkos::initialize();

  const KokkosFunctor kokkosFunctor(dx);

  // start timing
  tic = high_resolution_clock::now();

  const double kokkosIntegral = 0;
  // TODO: calculate integral using kokkos

  // stop timing
  toc = high_resolution_clock::now();
  const double kokkosElapsedTime =
    duration_cast<duration<double> >(toc - tic).count();

  // check the answer
  const double kokkosRelativeError =
    std::abs(libraryAnswer - kokkosIntegral) / std::abs(libraryAnswer);
  if (kokkosRelativeError > 1e-3) {
    fprintf(stderr, "our answer is too far off: %15.8e instead of %15.8e\n",
            kokkosIntegral, libraryAnswer);
    exit(1);
  }

  // output speedup
  printf("kokkos : time %8.2e speedup %8.2e\n",
         kokkosElapsedTime,
         fastSerialElapsedTime / kokkosElapsedTime);

  Kokkos::finalize();

  // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  // ********************** </do kokkos> ***************************
  // ===============================================================

  return 0;
}
