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

  inline
  double &
  operator()(const unsigned int index) {
    return _data[index];
  }

  inline 
  double 
  operator()(const unsigned int index) const {
    return _data[index];
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
  double &
  operator()(const unsigned int index) {
    return _data[index];
  }

  inline 
  double 
  operator()(const unsigned int index) const {
    return _data[index];
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
  ColMajorMatrix * _rightMatrix;
  RowMajorMatrix * _leftMatrix;
  RowMajorMatrix * _resultMatrix;

  TbbFunctorNaive(const unsigned int matrixSize, ColMajorMatrix * rightMatrix,
		  RowMajorMatrix * leftMatrix, RowMajorMatrix * resultMatrix) :
    _matrixSize(matrixSize), _rightMatrix(rightMatrix), _leftMatrix(leftMatrix),
    _resultMatrix(resultMatrix) {
  }
    
  void operator()(const tbb::blocked_range<size_t> & range) const {
    // TODO: something!
    for (unsigned int i = range.begin(); i < range.end(); i++) {
	for (unsigned int j = 0; j < _matrixSize; j++) {
	    for (unsigned int k = 0; k < _matrixSize; k++) {
		(*_resultMatrix)(i, j) += (*_leftMatrix)(i, k) * (*_rightMatrix)(k, j);
	    }
	}
    }
  }
    
/*
    void operator()(const tbb::blocked_range<size_t> & range) const {
	for (unsigned int i = range.begin(); i < range.end(); i++) {
	    unsigned int bound = (i/_matrixSize)*_matrixSize;
	    unsigned int col = i%_matrixSize;
	    for (unsigned int j = 0; j < _matrixSize; j++) {
		(*_resultMatrix)(i) += (*_leftMatrix)(bound+j) * (*_rightMatrix)(col+j*_matrixSize);
	    }
	}
    }
*/
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
	unsigned int loopBound = _matrixSize/_tileSize;
	for (unsigned int i = range.begin(); i < range.end(); i++) {
	    for (unsigned int j = 0; j < loopBound; j++) {
		for (unsigned int k = 0; k < loopBound; k++) {
		    for (unsigned int i1 = 0; i1 < _tileSize; i1++) {
			for (unsigned int j1 = 0; j1 < _tileSize; j1++) {
			    for (unsigned int k1 = 0; k1 < _tileSize; k1++) {
				(*_tiledResultMatrix)[i*_tileSize*_matrixSize + j*_tileSize + i1*_matrixSize + j1] +=
				(*_tiledLeftMatrix)[i*_tileSize*_matrixSize + j*_tileSize + i1*_matrixSize + k1] *
				(*_tiledRightMatrix)[j*_tileSize*_matrixSize + i*_tileSize + k1*_matrixSize + j1];       
			    }
			}
		    }
		}
	    }
	}
  }

private:
  TbbFunctorTiled();

};

typedef Kokkos::View<double*> view_type;
typedef view_type::HostMirror host_view_type;

struct KokkosFunctor {

  const unsigned int _matrixSize;
  view_type _rightMatrix;
  view_type _leftMatrix;
  view_type _resultMatrix;


  KokkosFunctor(const unsigned int matrixSize, view_type rightMatrix,
		view_type leftMatrix, view_type resultMatrix) :
    _matrixSize(matrixSize), _rightMatrix(rightMatrix), _leftMatrix(leftMatrix),
    _resultMatrix(resultMatrix) {
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const unsigned int elementIndex) const {
    // TODO: something!
    double sum = 0.0;
    for (unsigned int k = 0; k < _matrixSize; k++) {
	sum += _leftMatrix((elementIndex/_matrixSize)*_matrixSize+k) *
	       _rightMatrix(k*_matrixSize + elementIndex%_matrixSize);
    }
    _resultMatrix(elementIndex) = sum;
  }

private:
  KokkosFunctor();

};

int main(int argc, char* argv[]) {

  // a couple of inputs.  change the numberOfIntervals to control the amount
  //  of work done
  const unsigned int matrixSize = 512*3;
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
    for (unsigned int row = 0; row < matrixSize; ++row) {
      for (unsigned int col = 0; col < matrixSize; ++col) {
        resultMatrix(row, col) = 0;
        for (unsigned int dummy = 0; dummy < matrixSize; ++dummy) {
          resultMatrix(row, col) +=
            leftMatrix(row, dummy) * rightMatrixCol(dummy, col);
        }
      }
    }
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
    const TbbFunctorNaive tbbFunctor(matrixSize, &rightMatrixCol, &leftMatrix,
				    &resultMatrix);

    // start timing
    tic = high_resolution_clock::now();
    for (unsigned int repeatIndex = 0;
         repeatIndex < numberOfRepeats; ++repeatIndex) {
      // TODO: dispatch threads
	parallel_for(tbb::blocked_range<size_t>(0, matrixSize), tbbFunctor);
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
	#pragma omp parallel for
	    //Figure out how to do it!
	    for (unsigned int i = 0; i < matrixSize; i++) {
		for (unsigned int j = 0; j < matrixSize; j++) {
		    for (unsigned int k = 0; k < matrixSize; k++) {
			resultMatrix(i, j) += leftMatrix(i,k) * rightMatrixCol(k,j);
		    }
		}
	    }
	
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
  double checkSumKokkos = 0;
  for (unsigned int repeatIndex = 0;
       repeatIndex < numberOfRepeats; ++repeatIndex) {
    // TODO: do kokkos calculation
    
    // Make views for Kokkos Matrices
    view_type rightMat("right", matrixSize*matrixSize);
    view_type leftMat("left", matrixSize*matrixSize);
    view_type resultMat("result", matrixSize*matrixSize);
    
    //Make Host view types
    host_view_type h_result = Kokkos::create_mirror_view(resultMat);
    host_view_type h_right = Kokkos::create_mirror_view(rightMat);
    host_view_type h_left = Kokkos::create_mirror_view(leftMat);
   
    for (unsigned int i = 0; i < matrixSize; i++) {
	for (unsigned int j = 0; j < matrixSize; j++) {
	    h_left(i*matrixSize + j) = leftMatrix(i, j);
	    h_right(i*matrixSize + j) = rightMatrixRow(i, j);
	}
    } 
    
    //Try to deep copy the right and left matrix over?????
    Kokkos::deep_copy(rightMat, h_right);
    Kokkos::deep_copy(leftMat, h_left);

    //Make the functor
    KokkosFunctor kokkosFunctor(matrixSize, rightMat, leftMat, resultMat);
    
    Kokkos::parallel_for(matrixSize*matrixSize, kokkosFunctor);
    Kokkos::fence();
    
    // Deep copy it back
    Kokkos::deep_copy(h_result, resultMat);
     
  for (unsigned int i = 0; i < matrixSize; i++) {
    for (unsigned int j = 0; j < matrixSize; j++) {
	checkSumKokkos += h_result(i*matrixSize + j);
    }
  }
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
  if (std::abs(cacheUnfriendlyCheckSum - checkSumKokkos) / cacheUnfriendlyCheckSum < 1e-3) {
    printf("%-38s : time %6.2f speedup w.r.t. unfriendly %6.2f, w.r.t. friendly %6.2f\n",
           methodName,
           kokkosElapsedTime,
           cacheUnfriendlyElapsedTime / kokkosElapsedTime,
           cacheFriendlyElapsedTime / kokkosElapsedTime);
  } else {
    printf("%-38s : incorrect checksum %lf instead of %lf\n",
           methodName, checkSumKokkos, cacheUnfriendlyCheckSum);
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
  const vector<unsigned int> tileSizes = {16};

  for (const unsigned int tileSize : tileSizes) {

    vector<double> tiledLeftMatrix(matrixSize * matrixSize,
                                   std::numeric_limits<double>::quiet_NaN());
    vector<double> tiledRightMatrix(matrixSize * matrixSize,
                                    std::numeric_limits<double>::quiet_NaN());
    vector<double> tiledResultMatrix(matrixSize * matrixSize, 0);
    // TODO: form left matrix
    for (unsigned int i = 0; i < matrixSize; i++) {
	for (unsigned int j = 0; j < matrixSize; j++) {
	    tiledLeftMatrix[i*matrixSize + j] = leftMatrix(i, j);
	    tiledRightMatrix[i*matrixSize + j] = rightMatrixRow(i, j);
	}
    }
    // TODO: form right matrix

    // ===============================================================
    // ********************** < do vanilla tiled> ********************
    // vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

    double tiledElapsedTime = 0;
    tic = high_resolution_clock::now();
    for (unsigned int repeatIndex = 0;
         repeatIndex < numberOfRepeats; ++repeatIndex) {
      // TODO: do tiled matrix multiplication
	unsigned int loopBound = matrixSize/tileSize;
	for (unsigned int i = 0; i < loopBound; i++) {
	    for (unsigned int j = 0; j < loopBound; j++) {
		for (unsigned int k = 0; k < loopBound; k++) {
		    for (unsigned int i1 = 0; i1 < tileSize; i1++) {
			for (unsigned int j1 = 0; j1 < tileSize; j1++) {
			    for (unsigned int k1 = 0; k1 < tileSize; k1++) {
				tiledResultMatrix[i*tileSize*matrixSize + j*tileSize + i1*matrixSize + j1] +=
				tiledLeftMatrix[i*tileSize*matrixSize + j*tileSize + i1*matrixSize + k1] *
				tiledRightMatrix[j*tileSize*matrixSize + i*tileSize + k1*matrixSize + j1];       
			    }
			}
		    }
		}
	    }
	}
    }
    toc = high_resolution_clock::now();
    tiledElapsedTime = duration_cast<duration<double>>(toc - tic).count();
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
      tic = high_resolution_clock::now();
      
      for (unsigned int repeatIndex = 0;
           repeatIndex < numberOfRepeats; ++repeatIndex) {
        // TODO: something!
	parallel_for(tbb::blocked_range<size_t>(0, matrixSize/tileSize), tbbFunctor);
      }
      toc = high_resolution_clock::now();
      tbbElapsedTime =
	duration_cast<duration<double> >(toc - tic).count();


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
      tic = high_resolution_clock::now();
      for (unsigned int repeatIndex = 0;
           repeatIndex < numberOfRepeats; ++repeatIndex) {
        // TODO: something!
	unsigned int loopBound = matrixSize/tileSize;
	#pragma omp for
	for (unsigned int i = 0; i < loopBound; i++) {
	    for (unsigned int j = 0; j < loopBound; j++) {
		for (unsigned int k = 0; k < loopBound; k++) {
		    for (unsigned int i1 = 0; i1 < tileSize; i1++) {
			for (unsigned int j1 = 0; j1 < tileSize; j1++) {
			    for (unsigned int k1 = 0; k1 < tileSize; k1++) {
				tiledResultMatrix[i*tileSize*matrixSize + j*tileSize + i1*matrixSize + j1] +=
				tiledLeftMatrix[i*tileSize*matrixSize + j*tileSize + i1*matrixSize + k1] *
				tiledRightMatrix[j*tileSize*matrixSize + i*tileSize + k1*matrixSize + j1];       
			    }
			}
		    }
		}
	    }
	}
      }

      toc = high_resolution_clock::now();
      ompElapsedTime =
	duration_cast<duration<double> >(toc - tic).count();

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
