include ../../src/Makefile.inc

# customise the PATHs for FFTW and Gloal Arrays toolkit as needed

FFTW3_PATH=$(HOME)/software/build/fftw-3.3
GA_PATH=$(HOME)/software/build/ga-5-0-2
INCLUDE = -I../../include -I$(FFTW3_PATH)/include -I$(GA_PATH)/include

LIBS = -L../../lib -l2decomp_fft -L$(FFTW3_PATH)/lib -lfftw3 -lfftw3f -L$(GA_PATH)/lib -lga

all: ga_test

ga_test: ga_test.o
	$(F90) -o $@ $@.o $(LIBS)

clean:
	rm -f *.o ga_test

%.o : %.f90
	$(F90) $(INCLUDE) $(OPTIONS) $(F90FLAGS) -c $<
