########################################################################
# Makefile for building the 2DECOMP&FFT library
########################################################################

# Normally, do not change anything here. Modify the platform-dependent
# Makefile.inc file instead

include Makefile.inc

SRCS = decomp_2d.f90 io.f90 glassman.f90 fft_$(FFT).f90 fft_cufft.f90

# special treatment for shared-memory code
ifneq (,$(findstring DSHM,$(OPTIONS)))
SRCS := FreeIPC.f90 $(SRCS)   
OBJS =	$(SRCS:.f90=.o) alloc_shm.o FreeIPC_c.o
else
OBJS =	$(SRCS:.f90=.o)
endif	

# special treatment for Intel MKL (need to build MKL mod files)
MKL_MOD=
ifeq ($(FFT),mkl)
   MKL_MOD=mkl_mod
endif

all: lib examples

lib: lib2decomp_fft.a

includes: $(MKL_MOD) $(OBJS)
	@mkdir -p ../include
	@[ -e decomp_2d.mod ] && mv *.mod ../include || :

lib2decomp_fft.a: $(MKL_MOD) $(OBJS)
	@mkdir -p ../include
	@ar -csr $@ $(OBJS)
	@echo '******  Installing 2decomp modules and library  ******'
	@[ -e decomp_2d.mod ] && cp *.mod ../include || :
	@cp lib2decomp_fft.a ../lib

alloc_shm.o: alloc_shm.c
	$(CC) $(CFLAGS) -c $<

FreeIPC_c.o: FreeIPC_c.c
	$(CC) $(CFLAGS) -c $<

include depend.mk

.PHONY: mkl_mod
mkl_mod:
	echo Building Intel MKL mod files...
	$(IFORT) -c $(MKL_PATH)/include/mkl_dfti.f90

.PHONY: examples
examples: lib

.PHONY: test
test: examples

.PHONY: clean
clean:
	rm -f *.o *.mod lib*.a

.PHONY: realclean
realclean: clean
	rm -f *~ \#*\#

%.o : %.f90
	$(F90) $(OPTIONS) $(F90FLAGS) -c $<
