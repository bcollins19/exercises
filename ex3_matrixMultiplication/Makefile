# CUDA code generation flags
GENCODE_FLAGS   := -gencode arch=compute_30,code=sm_30 -gencode arch=compute_35,code=sm_35

# OS-specific build flags
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
  # compiler
  CC              := /opt/local/bin/g++-mp-4.8
  # cuda
  CUDA_PATH       ?= /usr/local/cuda
  CUDA_LIB_PATH   ?= $(CUDA_PATH)/lib/
  LD_FLAGS        += -Xlinker -rpath $(CUDA_LIB_PATH) -L$(CUDA_LIB_PATH) -lcudart -lcufft -lcublas
  # tbb
  TBB_DIR         += /opt/tbb-4.2
  # kokkos
  KOKKOS_DIR      := /opt/kokkos
  # eigen
  EIGEN_DIR       := /opt/eigen
  # avx vectorization support
  #CC_FLAGS        += -mavx -Wa,-q
else
  # compiler
  CC              := g++
  # cuda
  CUDA_PATH       ?= /usr/local/cuda-6.5
  CUDA_LIB_PATH   ?= $(CUDA_PATH)/lib64/
  ifeq ($(OS_SIZE),32)
    CC_FLAGS      += -m32
  else
    CC_FLAGS      += -m64
  endif
  LD_FLAGS      += -L$(CUDA_LIB_PATH) -lcudart -lcufft -lcublas
  # tbb
  TBB_DIR         += /clinic/2014/sandia14/tbb-4.3
  # kokkos
  KOKKOS_DIR      := /clinic/2014/sandia14/kokkos
  # eigen
  EIGEN_DIR       := /clinic/2014/sandia14/eigen
  # tbb
  CC_INCLUDE      += -I$(TBB_DIR)/include
  LD_FLAGS        += -L$(TBB_DIR)/lib -ltbb
  # avx vectorization support
  CC_FLAGS        += -mavx
endif

# cuda helper functions
NVCC_INCLUDE      += -I$(CUDA_PATH)/samples/common/inc/

# kokkos
CC_INCLUDE        += -I$(KOKKOS_DIR)/core/src -I$(KOKKOS_DIR)/containers/src -I$(KOKKOS_DIR)/algorithms/src -I$(KOKKOS_DIR)/linalg/src -DKOKKOS_DONT_INCLUDE_CORE_CONFIG_H  -DKOKKOS_HAVE_PTHREAD -DKOKKOS_HAVE_OPENMP
LD_FLAGS          += -L$(KOKKOS_DIR) -lkokkoscore

# eigen
CC_INCLUDE        += -I$(EIGEN_DIR)

# OS-architecture specific flags
ifeq ($(OS_SIZE),32)
      NVCCFLAGS := -m32
else
      NVCCFLAGS := -m64
endif

CUDA_INC_PATH   ?= $(CUDA_PATH)/include
CUDA_BIN_PATH   ?= $(CUDA_PATH)/bin
NVCC            ?= $(CUDA_BIN_PATH)/nvcc

CC_FLAGS   += -std=c++11 -fopenmp

TARGETS = MatrixMultiplication

all: $(TARGETS)

MatrixMultiplication: MatrixMultiplication.cc MatrixMultiplication_cuda.o
	$(CC) $< -o $@ MatrixMultiplication_cuda.o -O3 $(CC_INCLUDE) $(LD_FLAGS) $(CC_FLAGS) -Wall

MatrixMultiplication_cuda.o: MatrixMultiplication_cuda.cu
	$(NVCC) $(NVCCFLAGS) -O3 $(GENCODE_FLAGS) -I$(CUDA_INC_PATH) $(NVCC_INCLUDE) -o $@ -c $<

clean:
	rm -f *.o $(TARGETS)

again: clean $(TARGETS)
