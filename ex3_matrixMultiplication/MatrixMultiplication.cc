// -*- C++ -*-
// matrixMultiplication.cc
// a huge comparison of doing naive and tiled matrix multiplication using many
//  different methods and technologies

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <array>
#include <string>
#include <chrono>
#include <algorithm>

// yucky, but for asking the system how many cores we have
#include <unistd.h>

// header file for openmp
#include <omp.h>

// header files for tbb
#include <tbb/blocked_range.h>
#include <tbb/parallel_reduce.h>
#include <tbb/parallel_for.h>
#include <tbb/task_scheduler_init.h>

// header files for cuda implementation
#include "MatrixMultiplication_cuda.cuh"

// header files for eigen
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-local-typedefs"
#include <Eigen/Core>
#pragma GCC diagnostic pop

// header files for kokkos
#include <Kokkos_Core.hpp>

using std::string;
using std::vector;
using std::array;
using std::chrono::high_resolution_clock;
using std::chrono::duration;
using std::chrono::duration_cast;

struct ColMajorMatrix {
  const unsigned int _matrixSize;
  vector<double> _data;

  ColMajorMatrix(const unsigned int matrixSize) :
    _matrixSize(matrixSize), _data(_matrixSize*_matrixSize) {
  }

  inline
  double &
  operator()(const unsigned int row, const unsigned int col) {
    return _data[row + col * _matrixSize];
  }

  inline
  double
  operator()(const unsigned int row, const unsigned int col) const {
    return _data[row + col * _matrixSize];
  }

  void
  fill(const double value) {
    std::fill(_data.begin(), _data.end(), 0);
  }
};

struct RowMajorMatrix {
  const unsigned int _matrixSize;
  vector<double> _data;

  RowMajorMatrix(const unsigned int matrixSize) :
    _matrixSize(matrixSize), _data(_matrixSize*_matrixSize) {
  }

  inline
  double &
  operator()(const unsigned int row, const unsigned int col) {
    return _data[row * _matrixSize + col];
  }

  inline
  double
  operator()(const unsigned int row, const unsigned int col) const {
    return _data[row * _matrixSize + col];
  }

  void
  fill(const double value) {
    std::fill(_data.begin(), _data.end(), 0);
  }
};

class TbbFunctorNaive {
public:

  const unsigned int _matrixSize;

  TbbFunctorNaive(const unsigned int matrixSize) :
    _matrixSize(matrixSize) {
  }

  void operator()(const tbb::blocked_range<size_t> & range) const {
    // TODO: something!
  }

private:
  TbbFunctorNaive();

};

class TbbFunctorTiled {
public:

  const unsigned int _matrixSize;
  const unsigned int _tileSize;
  const vector<double> * const _tiledLeftMatrix;
  const vector<double> * const _tiledRightMatrix;
  vector<double> * const _tiledResultMatrix;

  TbbFunctorTiled(const unsigned int matrixSize,
                  const unsigned int tileSize,
                  const vector<double> * const tiledLeftMatrix,
                  const vector<double> * const tiledRightMatrix,
                  vector<double> * const tiledResultMatrix) :
    _matrixSize(matrixSize), _tileSize(tileSize),
    _tiledLeftMatrix(tiledLeftMatrix),
    _tiledRightMatrix(tiledRightMatrix),
    _tiledResultMatrix(tiledResultMatrix) {
  }

  void operator()(const tbb::blocked_range<size_t> & range) const {
    // TODO: something!
  }

private:
  TbbFunctorTiled();

};

struct KokkosFunctor {

  const unsigned int _matrixSize;

  KokkosFunctor(const unsigned int matrixSize) :
    _matrixSize(matrixSize) {
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const unsigned int elementIndex) const {
    // TODO: something!
  }

private:
  KokkosFunctor();

};

int main(int argc, char* argv[]) {

  // a couple of inputs.  change the numberOfIntervals to control the amount
  //  of work done
  const unsigned int matrixSize = 512 * 3;
  const unsigned int numberOfRepeats = 1;

  // we will repeat the computation for each of the numbers of threads
  vector<unsigned int> numberOfThreadsArray;
  //numberOfThreadsArray.push_back(1);
  //numberOfThreadsArray.push_back(2);
  //numberOfThreadsArray.push_back(4);
  //numberOfThreadsArray.push_back(8);
  //numberOfThreadsArray.push_back(16);
  //numberOfThreadsArray.push_back(24);
  numberOfThreadsArray.push_back(sysconf(_SC_NPROCESSORS_ONLN));

  printf("using a matrix size of %u\n", matrixSize);
  char methodName[500];

  // these are c++ timers...for timing
  high_resolution_clock::time_point tic;
  high_resolution_clock::time_point toc;

  // create a c++11 random number generator
  std::mt19937 randomNumberEngine;
  std::uniform_real_distribution<double> randomNumberGenerator(0, 1);

  RowMajorMatrix leftMatrix(matrixSize);
  RowMajorMatrix rightMatrixRow(matrixSize);
  ColMajorMatrix rightMatrixCol(matrixSize);
  RowMajorMatrix resultMatrix(matrixSize);
  for (unsigned int row = 0; row < matrixSize; ++row) {
    for (unsigned int col = 0; col < matrixSize; ++col) {
      leftMatrix(row, col) = randomNumberGenerator(randomNumberEngine);
      rightMatrixRow(row, col) = randomNumberGenerator(randomNumberEngine);
      rightMatrixCol(row, col) = rightMatrixRow(row, col);
    }
  }

  // ===============================================================
  // ********************** < do cache unfriendly> *****************
  // vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

  tic = high_resolution_clock::now();
  for (unsigned int repeatIndex = 0;
       repeatIndex < numberOfRepeats; ++repeatIndex) {
    for (unsigned int row = 0; row < matrixSize; ++row) {
      for (unsigned int col = 0; col < matrixSize; ++col) {
        resultMatrix(row, col) = 0;
        for (unsigned int dummy = 0; dummy < matrixSize; ++dummy) {
          resultMatrix(row, col) +=
            leftMatrix(row, dummy) * rightMatrixRow(dummy, col);
        }
      }
    }
  }
  toc = high_resolution_clock::now();
  const double cacheUnfriendlyElapsedTime =
    duration_cast<duration<double> >(toc - tic).count();

  double cacheUnfriendlyCheckSum = 0;
  for (unsigned int row = 0; row < matrixSize; ++row) {
    for (unsigned int col = 0; col < matrixSize; ++col) {
      cacheUnfriendlyCheckSum += resultMatrix(row, col);
    }
  }
  printf("%-38s : time %6.2f seconds\n",
         "cache unfriendly", cacheUnfriendlyElapsedTime);

  // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  // ********************** </do cache unfriendly> *****************
  // ===============================================================

  resultMatrix.fill(0);

  // ===============================================================
  // ********************** < do cache friendly> *******************
  // vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv


  tic = high_resolution_clock::now();

  for (unsigned int repeatIndex = 0;
       repeatIndex < numberOfRepeats; ++repeatIndex) {
    // TODO: do cache-friendly multiplication
  }

  toc = high_resolution_clock::now();
  const double cacheFriendlyElapsedTime =
    duration_cast<duration<double> >(toc - tic).count();

  double cacheFriendlyCheckSum = 0;
  for (unsigned int row = 0; row < matrixSize; ++row) {
    for (unsigned int col = 0; col < matrixSize; ++col) {
      cacheFriendlyCheckSum += resultMatrix(row, col);
    }
  }
  sprintf(methodName, "cache friendly");
  if (std::abs(cacheUnfriendlyCheckSum - cacheFriendlyCheckSum) / cacheUnfriendlyCheckSum < 1e-3) {
    printf("%-38s : time %6.2f speedup w.r.t. unfriendly %6.2f\n",
           methodName,
           cacheFriendlyElapsedTime,
           cacheUnfriendlyElapsedTime / cacheFriendlyElapsedTime);
  } else {
    printf("%-38s : incorrect checksum %lf instead of %lf\n",
           methodName, cacheFriendlyCheckSum, cacheUnfriendlyCheckSum);
  }

  // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  // ********************** </do cache friendly> *******************
  // ===============================================================

  resultMatrix.fill(0);

  // ===============================================================
  // ********************** < do naive tbb> ************************
  // vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

  // for each number of threads
  for (const unsigned int numberOfThreads :
         numberOfThreadsArray) {

    // initialize tbb's threading system for this number of threads
    tbb::task_scheduler_init init(numberOfThreads);

    // prepare the tbb functor.
    const TbbFunctorNaive tbbFunctor(matrixSize);

    // start timing
    tic = high_resolution_clock::now();
    for (unsigned int repeatIndex = 0;
         repeatIndex < numberOfRepeats; ++repeatIndex) {
      // TODO: dispatch threads
    }
    // stop timing
    toc = high_resolution_clock::now();
    const double tbbElapsedTime =
      duration_cast<duration<double> >(toc - tic).count();

    // check the answer
    double tbbCheckSum = 0;
    for (unsigned int row = 0; row < matrixSize; ++row) {
      for (unsigned int col = 0; col < matrixSize; ++col) {
        tbbCheckSum += resultMatrix(row, col);
      }
    }
    sprintf(methodName, "naive tbb, %3u threads", numberOfThreads);
    if (std::abs(cacheUnfriendlyCheckSum - tbbCheckSum) / cacheUnfriendlyCheckSum < 1e-3) {
      printf("%-38s : time %6.2f speedup w.r.t. unfriendly %6.2f, w.r.t. friendly %6.2f (%%%5.1f of ideal)\n",
             methodName,
             tbbElapsedTime,
             cacheUnfriendlyElapsedTime / tbbElapsedTime,
             cacheFriendlyElapsedTime / tbbElapsedTime,
             100. * cacheFriendlyElapsedTime / tbbElapsedTime / numberOfThreads);
    } else {
      printf("%-38s : incorrect checksum %lf instead of %lf\n",
             methodName, tbbCheckSum, cacheUnfriendlyCheckSum);
    }
  }

  // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  // ********************** </do naive tbb> ************************
  // ===============================================================

  resultMatrix.fill(0);

  // ===============================================================
  // ********************** < do naive openmp> *********************
  // vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

  // for each number of threads
  for (const unsigned int numberOfThreads :
         numberOfThreadsArray) {

    // set the number of threads for openmp
    omp_set_num_threads(numberOfThreads);

    // start timing
    tic = high_resolution_clock::now();

    for (unsigned int repeatIndex = 0;
         repeatIndex < numberOfRepeats; ++repeatIndex) {
      // TODO: do openmp
    }

    // stop timing
    toc = high_resolution_clock::now();
    const double ompElapsedTime =
      duration_cast<duration<double> >(toc - tic).count();

    // check the answer
    double ompCheckSum = 0;
    for (unsigned int row = 0; row < matrixSize; ++row) {
      for (unsigned int col = 0; col < matrixSize; ++col) {
        ompCheckSum += resultMatrix(row, col);
      }
    }
    sprintf(methodName, "naive omp, %3u threads", numberOfThreads);
    if (std::abs(cacheUnfriendlyCheckSum - ompCheckSum) / cacheUnfriendlyCheckSum < 1e-3) {
      printf("%-38s : time %6.2f speedup w.r.t. unfriendly %6.2f, w.r.t. friendly %6.2f (%%%5.1f of ideal)\n",
             methodName,
             ompElapsedTime,
             cacheUnfriendlyElapsedTime / ompElapsedTime,
             cacheFriendlyElapsedTime / ompElapsedTime,
             100. * cacheFriendlyElapsedTime / ompElapsedTime / numberOfThreads);
    } else {
      printf("%-38s : incorrect checksum %lf instead of %lf\n",
             methodName, ompCheckSum, cacheUnfriendlyCheckSum);
    }
  }

  // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  // ********************** </do naive openmp> *********************
  // ===============================================================

  resultMatrix.fill(0);

  // ===============================================================
  // ********************** < do cuda> *****************************
  // vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

  // we will repeat the computation for each of the numbers of threads
  const vector<unsigned int> threadsPerBlockArray = {256};
  const vector<unsigned int> maxNumberOfBlocksArray = {10000};

  // warm up cuda
  {
    const unsigned int warmUpMaxNumberOfBlocks = 1e4;
    const unsigned int warmUpThreadsPerBlock   = 256;
    cudaDoMatrixMultiplication(warmUpMaxNumberOfBlocks,
                               warmUpThreadsPerBlock,
                               matrixSize);
  }

  // for each max number of blocks
  for (const unsigned int maxNumberOfBlocks :
         maxNumberOfBlocksArray) {
    // for each number of threads per block
    for (const unsigned int numberOfThreadsPerBlock :
           threadsPerBlockArray) {

      // start timing
      tic = high_resolution_clock::now();

      // do calculation with cuda for this number of threads per block
      for (unsigned int repeatIndex = 0;
           repeatIndex < numberOfRepeats; ++repeatIndex) {
        cudaDoMatrixMultiplication(maxNumberOfBlocks,
                                   numberOfThreadsPerBlock,
                                   matrixSize);
      }

      // stop timing
      toc = high_resolution_clock::now();
      const double cudaElapsedTime =
        duration_cast<duration<double> >(toc - tic).count();

      // check the answer
      double cudaCheckSum = 0;
      for (unsigned int row = 0; row < matrixSize; ++row) {
        for (unsigned int col = 0; col < matrixSize; ++col) {
          cudaCheckSum += resultMatrix(row, col);
        }
      }
      sprintf(methodName, "naive cuda %8.2e blocks %3u threads", double(maxNumberOfBlocks), numberOfThreadsPerBlock);
      if (std::abs(cacheUnfriendlyCheckSum - cudaCheckSum) / cacheUnfriendlyCheckSum < 1e-3) {
        printf("%-38s : time %6.2f speedup w.r.t. unfriendly %6.2f, w.r.t. friendly %6.2f\n",
               methodName,
               cudaElapsedTime,
               cacheUnfriendlyElapsedTime / cudaElapsedTime,
               cacheFriendlyElapsedTime / cudaElapsedTime);
      } else {
        printf("%-38s : incorrect checksum %lf instead of %lf\n",
               methodName, cudaCheckSum, cacheUnfriendlyCheckSum);
      }
    }
  }

  // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  // ********************** </do cuda> *****************************
  // ===============================================================

  resultMatrix.fill(0);

  // ===============================================================
  // ********************** < do kokkos> ***************************
  // vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

  Kokkos::initialize();

  //printf("kokkos is running on %s\n", typeid(Kokkos::DefaultExecutionSpace).name());

  // start timing
  tic = high_resolution_clock::now();

  for (unsigned int repeatIndex = 0;
       repeatIndex < numberOfRepeats; ++repeatIndex) {
    // TODO: do kokkos calculation
  }

  // stop timing
  toc = high_resolution_clock::now();
  const double kokkosElapsedTime =
    duration_cast<duration<double> >(toc - tic).count();

  // check the answer
  double kokkosCheckSum = 0;
  for (unsigned int row = 0; row < matrixSize; ++row) {
    for (unsigned int col = 0; col < matrixSize; ++col) {
      kokkosCheckSum += resultMatrix(row, col);
    }
  }
  sprintf(methodName, "naive kokkos");
  if (std::abs(cacheUnfriendlyCheckSum - kokkosCheckSum) / cacheUnfriendlyCheckSum < 1e-3) {
    printf("%-38s : time %6.2f speedup w.r.t. unfriendly %6.2f, w.r.t. friendly %6.2f\n",
           methodName,
           kokkosElapsedTime,
           cacheUnfriendlyElapsedTime / kokkosElapsedTime,
           cacheFriendlyElapsedTime / kokkosElapsedTime);
  } else {
    printf("%-38s : incorrect checksum %lf instead of %lf\n",
           methodName, kokkosCheckSum, cacheUnfriendlyCheckSum);
  }

  Kokkos::finalize();
  // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  // ********************** </do kokkos> ***************************
  // ===============================================================

  resultMatrix.fill(0);

  // ===============================================================
  // ********************** < do eigen> ****************************
  // vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

  Eigen::MatrixXd eigenLeftMatrix(matrixSize, matrixSize);
  Eigen::MatrixXd eigenRightMatrix(matrixSize, matrixSize);
  Eigen::MatrixXd eigenResultMatrix(matrixSize, matrixSize);
  for (unsigned int row = 0; row < matrixSize; ++row) {
    for (unsigned int col = 0; col < matrixSize; ++col) {
      eigenLeftMatrix(row, col) = leftMatrix(row, col);
      eigenRightMatrix(row, col) = rightMatrixRow(row, col);
    }
  }

  // warm up eigen
  eigenResultMatrix = eigenLeftMatrix * eigenRightMatrix;

  // start timing
  tic = high_resolution_clock::now();

  for (unsigned int repeatIndex = 0;
       repeatIndex < numberOfRepeats; ++repeatIndex) {
    eigenResultMatrix = eigenLeftMatrix * eigenRightMatrix;
  }

  // stop timing
  toc = high_resolution_clock::now();
  const double eigenElapsedTime =
    duration_cast<duration<double> >(toc - tic).count();

  // check the answer
  double eigenCheckSum = 0;
  for (unsigned int row = 0; row < matrixSize; ++row) {
    for (unsigned int col = 0; col < matrixSize; ++col) {
      eigenCheckSum += eigenResultMatrix(row, col);
    }
  }
  sprintf(methodName, "eigen");
  if (std::abs(cacheUnfriendlyCheckSum - eigenCheckSum) / cacheUnfriendlyCheckSum < 1e-3) {
    printf("%-38s : time %6.2f speedup w.r.t. unfriendly %6.2f, w.r.t. friendly %6.2f\n",
           methodName,
           eigenElapsedTime,
           cacheUnfriendlyElapsedTime / eigenElapsedTime,
           cacheFriendlyElapsedTime / eigenElapsedTime);
  } else {
    printf("%-38s : incorrect checksum %lf instead of %lf\n",
           methodName, eigenCheckSum, cacheUnfriendlyCheckSum);
  }

  // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  // ********************** </do eigen> ****************************
  // ===============================================================

  resultMatrix.fill(0);

  // ===============================================================
  // ********************** < do tiled> ****************************
  // vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

  //const vector<unsigned int> tileSizes = {16, 32, 64};
  const vector<unsigned int> tileSizes = {};

  for (const unsigned int tileSize : tileSizes) {

    vector<double> tiledLeftMatrix(matrixSize * matrixSize,
                                   std::numeric_limits<double>::quiet_NaN());
    vector<double> tiledRightMatrix(matrixSize * matrixSize,
                                    std::numeric_limits<double>::quiet_NaN());
    vector<double> tiledResultMatrix(matrixSize * matrixSize, 0);
    // TODO: form left matrix
    // TODO: form right matrix

    // ===============================================================
    // ********************** < do vanilla tiled> ********************
    // vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

    double tiledElapsedTime = 0;
    for (unsigned int repeatIndex = 0;
         repeatIndex < numberOfRepeats; ++repeatIndex) {
      // TODO: do tiled matrix multiplication
    }
    // check the answer
    double tiledCheckSum = 0;
    for (const double entry : tiledResultMatrix) {
      tiledCheckSum += entry;
    }
    sprintf(methodName, "tileSize %3u", tileSize);
    if (std::abs(cacheUnfriendlyCheckSum - tiledCheckSum) / cacheUnfriendlyCheckSum < 1e-3) {
      printf("%-38s : time %6.2f speedup w.r.t. unfriendly %6.2f, w.r.t. friendly %6.2f\n",
             methodName,
             tiledElapsedTime,
             cacheUnfriendlyElapsedTime / tiledElapsedTime,
             cacheFriendlyElapsedTime / tiledElapsedTime);
    } else {
      printf("%-38s : incorrect checksum %lf instead of %lf\n",
             methodName, tiledCheckSum, cacheUnfriendlyCheckSum);
    }

    // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    // ********************** </do vanilla tiled> ********************
    // ===============================================================

    std::fill(tiledResultMatrix.begin(), tiledResultMatrix.end(), 0);

    // ===============================================================
    // ********************** < do tiled tbb> ************************
    // vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

    // for each number of threads
    for (const unsigned int numberOfThreads :
           numberOfThreadsArray) {

      // initialize tbb's threading system for this number of threads
      tbb::task_scheduler_init init(numberOfThreads);

      // prepare the tbb functor.
      const TbbFunctorTiled tbbFunctor(matrixSize,
                                       tileSize,
                                       &tiledLeftMatrix,
                                       &tiledRightMatrix,
                                       &tiledResultMatrix);

      double tbbElapsedTime = 0;
      for (unsigned int repeatIndex = 0;
           repeatIndex < numberOfRepeats; ++repeatIndex) {
        // TODO: something!
      }

      // check the answer
      double tbbCheckSum = 0;
      for (const double entry : tiledResultMatrix) {
        tbbCheckSum += entry;
      }
      sprintf(methodName, "tileSize %3u, %2u tbb threads", tileSize, numberOfThreads);
      if (std::abs(cacheUnfriendlyCheckSum - tbbCheckSum) / cacheUnfriendlyCheckSum < 1e-3) {
        printf("%-38s : time %6.2f speedup w.r.t. unfriendly %6.2f, w.r.t. friendly %6.2f (%%%5.1f of ideal)\n",
               methodName,
               tbbElapsedTime,
               cacheUnfriendlyElapsedTime / tbbElapsedTime,
               cacheFriendlyElapsedTime / tbbElapsedTime,
               100. * cacheFriendlyElapsedTime / tbbElapsedTime / numberOfThreads);
      } else {
        printf("%-38s : incorrect checksum %lf instead of %lf\n",
               methodName, tbbCheckSum, cacheUnfriendlyCheckSum);
      }
    }

    // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    // ********************** </do tiled tbb> ************************
    // ===============================================================

    std::fill(tiledResultMatrix.begin(), tiledResultMatrix.end(), 0);

    // ===============================================================
    // ********************** < do tiled openmp> *********************
    // vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

    // for each number of threads
    for (const unsigned int numberOfThreads :
           numberOfThreadsArray) {

      omp_set_num_threads(numberOfThreads);

      double ompElapsedTime = 0;
      for (unsigned int repeatIndex = 0;
           repeatIndex < numberOfRepeats; ++repeatIndex) {
        // TODO: something!
      }


      double ompCheckSum = 0;
      for (const double entry : tiledResultMatrix) {
        ompCheckSum += entry;
      }
      sprintf(methodName, "tileSize %3u, %2u omp threads", tileSize, numberOfThreads);
      if (std::abs(cacheUnfriendlyCheckSum - ompCheckSum) / cacheUnfriendlyCheckSum < 1e-3) {
        printf("%-38s : time %6.2f speedup w.r.t. unfriendly %6.2f, w.r.t. friendly %6.2f (%%%5.1f of ideal)\n",
               methodName,
               ompElapsedTime,
               cacheUnfriendlyElapsedTime / ompElapsedTime,
               cacheFriendlyElapsedTime / ompElapsedTime,
               100. * cacheFriendlyElapsedTime / ompElapsedTime / numberOfThreads);
      } else {
        printf("%-38s : incorrect checksum %lf instead of %lf\n",
               methodName, ompCheckSum, cacheUnfriendlyCheckSum);
      }
    }

    // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    // ********************** </do tiled openmp> *********************
    // ===============================================================

  }

  // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  // ********************** </do tiled> ****************************
  // ===============================================================

  // ===============================================================
  // ********************** < do cublas> ***************************
  // vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

  {
#if 0
    const int cudaDeviceId = 0;
    cudaDeviceProp deviceProp;
    checkCudaErrors(cudaGetDeviceProperties(&deviceProp, cudaDeviceId));
    printf("GPU Device %d: \"%s\" with compute capability %d.%d\n\n",
           cudaDeviceId, deviceProp.name, deviceProp.major, deviceProp.minor);
#endif

    // warm up cublas
    multiplyMatricesUsingCublas(matrixSize,
                                &leftMatrix(0, 0),
                                &rightMatrixRow(0, 0),
                                &resultMatrix(0, 0));
    // start timing
    tic = high_resolution_clock::now();

    for (unsigned int repeatIndex = 0;
         repeatIndex < numberOfRepeats; ++repeatIndex) {
      multiplyMatricesUsingCublas(matrixSize,
                                  &leftMatrix(0, 0),
                                  &rightMatrixRow(0, 0),
                                  &resultMatrix(0, 0));
    }

    // stop timing
    toc = high_resolution_clock::now();
    const double cublasElapsedTime =
      duration_cast<duration<double> >(toc - tic).count();

    // check the answer
    double cublasCheckSum = 0;
    for (unsigned int row = 0; row < matrixSize; ++row) {
      for (unsigned int col = 0; col < matrixSize; ++col) {
        cublasCheckSum += resultMatrix(row, col);
      }
    }
    sprintf(methodName, "cublas");
    if (std::abs(cacheUnfriendlyCheckSum - cublasCheckSum) / cacheUnfriendlyCheckSum < 1e-3) {
      printf("%-38s : time %6.2f speedup w.r.t. unfriendly %6.2f, w.r.t. friendly %6.2f\n",
             methodName,
             cublasElapsedTime,
             cacheUnfriendlyElapsedTime / cublasElapsedTime,
             cacheFriendlyElapsedTime / cublasElapsedTime);
    } else {
      printf("%-38s : incorrect checksum %lf instead of %lf\n",
             methodName, cublasCheckSum, cacheUnfriendlyCheckSum);
    }
  }

  // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  // ********************** </do cublas> ***************************
  // ===============================================================

  return 0;
}
