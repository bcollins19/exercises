// -*- C++ -*-
#ifndef MATRIXMULTIPLICATION_CUDA_CUH
#define MATRIXMULTIPLICATION_CUDA_CUH

void
cudaDoMatrixMultiplication(const unsigned int maxNumberOfBlocks,
                           const unsigned int numberOfThreadsPerBlock,
                           const unsigned int matrixSize);
void
multiplyMatricesUsingCublas(const unsigned int matrixSize,
                            const double * leftMatrix,
                            const double * rightMatrix,
                            double * resultMatrix);

#endif // MATRIXMULTIPLICATION_CUDA_CUH
