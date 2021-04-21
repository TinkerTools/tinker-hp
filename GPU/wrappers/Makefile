# Compilers
CUDACC   = nvcc
FC       = pgf90

# Libraries
#CUDA_DIR  = /usr/local/cuda-10.0
#CUDA_DIR  = $(CUDAROOT)
#CUDA_DIR  = $(PGI)/linux86-64-llvm/2019/cuda/10.0

# Flags
FCFLAGS   = -O3 -Mcuda=cc60,cc70#-rc=rc4.0
CUFLAGS   = -std=c++11 -O3 -gencode arch=compute_60,code=sm_60 \
                -gencode arch=compute_70,code=sm_70 \

# Directories
#INCLUDE       = -I$(CUDA_DIR)/include
INCLUDE       = -I.
TINKER_SOURCE =../source

# Objects
OBJ = thrust_cache.o csort.o thrust_module.o

# Targets
# -------

all: libwrapper.a

libwrapper.a: $(OBJ)
	@echo '****** Making libwrapper.a ******'
	$(AT)ar -crs $@ $(OBJ)

# dependencies
# ------------

thrust_cache.o: thrust_cache.h

# Generic Rules
# -------------

%.o: %.cu
	$(CUDACC) $(CUFLAGS) -c $< $(INCLUDE)

%.o: %.f90
	$(FC) $(FCFLAGS) -c $<

# Cleaning operation
clean:
	rm *.o *.mod *.a
