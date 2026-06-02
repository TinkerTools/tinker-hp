include ../../src/Makefile.inc

# customise the PATHs for FFTW and libNBC as needed

FFTW3_PATH=$(HOME)/software/build/fftw-3.3
INCLUDE = -I../../include -I$(FFTW3_PATH)/include
LIBNBC_PATH = $(HOME)/software/build/libNBC-1.1.1-thread

LIBS = -L../../lib -l2decomp_fft $(LIBNBC_PATH)/lib/libnbc.a -lstdc++ -lmpi_cxx -L$(FFTW3_PATH)/lib -lfftw3

all: blocking non_blocking

blocking: blocking.o
	$(F90) -o $@ $@.o $(LIBS)

non_blocking: non_blocking.o
	$(F90) -o $@ $@.o $(LIBS)

clean:
	rm -f *.o blocking non_blocking

%.o : %.f90
	$(F90) $(INCLUDE) $(OPTIONS) $(F90FLAGS) -c $<
