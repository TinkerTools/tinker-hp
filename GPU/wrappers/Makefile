# Compilers
CUDACC   = nvcc
FC       = pgf90

# Libraries
#CUDA_DIR  = /usr/local/cuda-10.0
#CUDA_DIR  = $(CUDAROOT)
#CUDA_DIR  = $(PGI)/linux86-64-llvm/2019/cuda/10.0

# Flags
FCFLAGS   = -O3 -Mcuda=cc60,cc70#-rc=rc4.0
CUFLAGS   = -std=c++14 -O3 -gencode arch=compute_70,code=sm_70 \
                -gencode arch=compute_80,code=sm_80 \

# Directories
#INCLUDE       = -I$(CUDA_DIR)/include
INCLUDE       = -I.
TINKER_SOURCE =../source

# Objects
OBJ = thrust_cache.o csort.o thrust_module.o

#
# Targets
#
all: libwrapper.a

libwrapper.a: $(OBJ)
	ar qc $@ $(OBJ)

#
# Generic Rules
#
thrust_cache.o: thrust_cache.cu thrust_cache.h
	$(CUDACC) $(CUFLAGS) -c $< $(INCLUDE)
%.o: %.cu
	$(CUDACC) $(CUFLAGS) -c $< $(INCLUDE)

%.o: %.f90
	$(FC) $(FCFLAGS) -c $<

# Cleaning operation
clean:
	rm *.o *.mod *.a
