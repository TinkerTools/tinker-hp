# Configurations on an ordinary Linux PCs for development
#   using GNU compiler, OpenMPI and various FFT libraries

# 2DECOMP&FFT options
#====================

# Choose pre-processing options
#   -DDOUBLE_PREC   - use double-precision (default single)
#   -DEVEN          - for evenly distributed data, use ALLTOALL
#   -DSHM	        - enable system V shared-memory implementation
#   -DOVERWRITE     - allow FFT input to be overwritten (save memory)
#   -DOCC           - use API for overlap communication and computation
#   -DGLOBAL_ARRAYS - turn on support for Global Arrays toolkit
#   -DT3PIO         - turn on LUSTRE IO optimisation code using T3PIO
OPTIONS=-DDOUBLE_PREC

# Choose one FFT engine, available options are:
#   acml       - AMD Core Math Library
#   ffte       - FFTE
#   fftpack5   - FFTPACK5
#   fftw3      - FFTW version 3.x
#   fftw3_f03  - FFTW 3.3-beta1 or later (with Fortran 2003 interface)
#   generic    - A general FFT algorithm (no 3rd-party library needed)
#   mkl        - Intel Math Kernel Library
FFT=generic

# Platform-dependent information - compiler, external library etc
#================================================================

# Inlcude path if necessary (change to your actual paths)
MKL_PATH=/opt/intel/mkl
FFTW_PATH=/opt/fftw-3.3
ifeq ($(FFT),mkl)
  INC=-I$(MKL_PATH)/include
  # Fortran compiler used to compile MKL mod files
  IFORT=gfortran
else ifeq ($(FFT),fftw3)
  INC=-I$(FFTW_PATH)/include
else ifeq ($(FFT),fftw3_f03)
  INC=-I$(FFTW_PATH)/include
else
  INC=
endif

ifneq (,$(findstring DGLOBAL_ARRAYS,$(OPTIONS)))
  GA_PATH=/opt/ga-5-0-2
  INC+= -I$(GA_PATH)/include
endif

ifneq (,$(findstring DT3PIO,$(OPTIONS)))
  T3PIO_PATH=/opt/t3pio
  INC+= -I$(T3PIO_PATH)/include
endif


#----------------------- Fortran Compiler ----------------------------
F90=mpif90

# enable preprocessing
CPPFLAGS=-cpp
# enable Cray pointer support if needed
CRAYPTR=-fcray-pointer
# optimisation or debugging flags
#OPTIM=-g -fcheck=all
OPTIM=-O3

F90FLAGS=$(OPTIM) $(CRAYPTR) $(CPPFLAGS) $(INC)
LDFLAGS=$(OPTIM)

#--------------------------C Compiler---------------------------------
CC=mpicc
CFLAGS=-O3

#-----------------------External Library------------------------------

# For FFTW
LIB_FFTW3=-L$(FFTW_PATH)/lib -lfftw3f -lfftw3

# For ACML
# This assumes that 32-bit ACML installed at /opt/. Adjust properly.
# It helps to create a symbolic link such as 'acml -> acml4.4.0'
# and update the symbolic link when future ACML version is installed. 
LIB_ACML=/opt/acml/gfortran32/lib/libacml.a -lrt

# For FFTPACK5 
FFTPACK5_PATH=$(HOME)/software/fftpack5
LIB_FFTPACK5=$(FFTPACK5_PATH)/libfftpack5.a

# For Intel MKL
MKL_LIB_PATH= $(MKL_PATH)/lib/ia32
LIB_MKL=-Wl,--start-group $(MKL_LIB_PATH)/libmkl_gf.a $(MKL_LIB_PATH)/libmkl_sequential.a $(MKL_LIB_PATH)/libmkl_core.a -Wl,--end-group -lpthread

# For FFTE
FFTE_PATH=$(HOME)/software/ffte-4.1
LIB_FFTE=$(FFTE_PATH)/libffte.a

ifeq ($(FFT),generic)
     LIBFFT=
else ifeq ($(FFT),acml)
     LIBFFT=$(LIB_ACML)
else ifeq ($(FFT),ffte)
     LIBFFT=$(LIB_FFTE)
else ifeq ($(FFT),fftpack5)
     LIBFFT=$(LIB_FFTPACK5)
else ifeq ($(FFT),fftw3)
     LIBFFT=$(LIB_FFTW3)
else ifeq ($(FFT),fftw3_f03)
     LIBFFT=$(LIB_FFTW3)
else ifeq ($(FFT),mkl)
     LIBFFT=$(LIB_MKL)
endif	
