// -*- C++ -*-
#include <cstdio>

#include <cuda_runtime.h>

#include "ex0_scalarIntegrator_cuda.cuh"

#define MAX_NUMBER_OF_PUPPIES 3

__constant__ double constantPuppies[MAX_NUMBER_OF_PUPPIES];

__global__
void
cudaDoScalarIntegration_kernel(double* output) {
  // block-wide reduction storage, size is determined by third kernel
  // launch argument (thing between <<< and >>>)
  extern __shared__ double contributions[];

  // TODO: do scalar integration somehow

  // reading from global memory
  *output = 5;

}

void
cudaDoScalarIntegration(const unsigned int numberOfThreadsPerBlock,
			const double startBound, double endBound,
			const double dx, double * const output) {

  // TODO: you have to do stuff in here, the junk below is just to show syntax

  // this is how to use constant memory:
  // make some stuff that we'll copy into constant memory
  double * iLikePuppies = new double[MAX_NUMBER_OF_PUPPIES];
  iLikePuppies[0] = 0.;
  iLikePuppies[1] = 1.;
  iLikePuppies[2] = 2.;
  // copy some junk into constant memory
  cudaMemcpyToSymbol(constantPuppies, iLikePuppies,
                     sizeof(double) * MAX_NUMBER_OF_PUPPIES);
  delete[] iLikePuppies;


  // using global memory
  unsigned int amountOfJunk = 10;
  // make some array from which we'll copy to the device
  double * junk = new double[amountOfJunk];
  // fill it with junk
  for (int i = 0; i < amountOfJunk; ++i) {
    junk[i] = -i;
  }
  // this is going to be a pointer to memory *on the device*
  double * dev_junk;
  // allocate room on the device
  cudaMalloc((void **) &dev_junk, amountOfJunk*sizeof(double));
  // copy junk from host to device
  cudaMemcpy(dev_junk, junk, amountOfJunk*sizeof(double),
             cudaMemcpyHostToDevice);
  delete[] junk;

  // TODO: calculate the number of blocks
  const unsigned int numberOfBlocks = 1;

  // allocate somewhere to put our result
  double *dev_output;
  cudaMalloc((void **) &dev_output, 1*sizeof(double));

  // run the kernel
  cudaDoScalarIntegration_kernel<<<numberOfBlocks,
    numberOfThreadsPerBlock,
    numberOfThreadsPerBlock*sizeof(double)>>>(dev_output);

  // copy over the output
  cudaMemcpy(output, dev_output, 1*sizeof(double), cudaMemcpyDeviceToHost);
  // make sure that everything in flight has been completed
  cudaDeviceSynchronize();

  // clean up
  cudaFree(dev_junk);
  cudaFree(dev_output);
}
