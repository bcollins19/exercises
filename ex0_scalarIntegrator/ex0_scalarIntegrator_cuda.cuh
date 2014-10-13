// -*- C++ -*-
#ifndef EX0_SCALARINTEGRATOR_CUDA_CUH
#define EX0_SCALARINTEGRATOR_CUDA_CUH

// lots of consts for defensive programming
void
cudaDoScalarIntegration(const unsigned int numberOfThreadsPerBlock,
			const double startBound, const double endBound,
			const double dx, double * const output);

#endif // EX0_SCALARINTEGRATOR_CUDA_CUH
