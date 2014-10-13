// -*- C++ -*-
#include <cstdio>

#include <cuda_runtime.h>
// These come from the cublas matrix multiplication example
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wswitch"
#include <cublas_v2.h>
#include <helper_functions.h>
#include <helper_cuda.h>
#pragma GCC diagnostic pop


#include "MatrixMultiplication_cuda.cuh"

__global__
void
cudaDoNaiveMatrixMultiplication_kernel(const unsigned int matrixSize,
                                                   const double * leftMatrix,
                                                   const double * rightMatrix,
                                                   double * resultMatrix) {
  // TODO: something!
}

void
cudaDoMatrixMultiplication(const unsigned int maxNumberOfBlocks,
                           const unsigned int numberOfThreadsPerBlock,
                           const unsigned int matrixSize) {

  // TODO: something!
}


void
multiplyMatricesUsingCublas(const unsigned int matrixSize,
                            const double * leftMatrix,
                            const double * rightMatrix,
                            double * resultMatrix) {

  const unsigned int numberOfEntries = matrixSize * matrixSize;

  // allocate device memory
  double * dev_leftMatrix;
  double * dev_rightMatrix;
  double * dev_resultMatrix;
  checkCudaErrors(cudaMalloc((void **) &dev_leftMatrix, numberOfEntries * sizeof(double)));
  checkCudaErrors(cudaMalloc((void **) &dev_rightMatrix, numberOfEntries * sizeof(double)));
  checkCudaErrors(cudaMalloc((void **) &dev_resultMatrix, numberOfEntries * sizeof(double)));
  // copy matrices to the device
  checkCudaErrors(cudaMemcpy(dev_leftMatrix, leftMatrix, numberOfEntries * sizeof(double), cudaMemcpyHostToDevice));
  checkCudaErrors(cudaMemcpy(dev_rightMatrix, rightMatrix, numberOfEntries * sizeof(double), cudaMemcpyHostToDevice));

  const double alpha = 1.0f;
  const double beta  = 0.0f;
  cublasHandle_t handle;

  checkCudaErrors(cublasCreate(&handle));

  // perform the multiply
  checkCudaErrors(cublasDgemm(handle,
                              CUBLAS_OP_N, CUBLAS_OP_N, // don't transpose
                              matrixSize, matrixSize, matrixSize, // sizes
                              &alpha, // no scalar premultiply
                              dev_rightMatrix, matrixSize, // left matrix
                              dev_leftMatrix, matrixSize, // right matrix
                              &beta, // don't premultiply result by anything
                              dev_resultMatrix, matrixSize));

  // copy result from device to host
  checkCudaErrors(cudaMemcpy(resultMatrix, dev_resultMatrix,
                             numberOfEntries * sizeof(double),
                             cudaMemcpyDeviceToHost));

  // Destroy the handle
  checkCudaErrors(cublasDestroy(handle));

  // clean up memory
  checkCudaErrors(cudaFree(dev_leftMatrix));
  checkCudaErrors(cudaFree(dev_rightMatrix));
  checkCudaErrors(cudaFree(dev_resultMatrix));
}
