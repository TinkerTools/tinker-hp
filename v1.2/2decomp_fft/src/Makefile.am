#
#     Sorbonne University
#     Washington University in Saint Louis
#     University of Texas at Austin
#

####################################################
# Compilation of 2decomp_fft library for Tinker-HP #
#                                                  #
# Declare the library to create                    #
####################################################

noinst_LIBRARIES = lib2decomp_fft.a

#
# Declare her sources
#

lib2decomp_fft_a_SOURCES =
lib2decomp_fft_a_SOURCES += decomp_2d.f90
lib2decomp_fft_a_SOURCES += io.f90
lib2decomp_fft_a_SOURCES += glassman.f90

#
# States that depend.mk must exist to allows parallel make to work
# It gives the dependencies. So, it gives  the order of compilation
#

BUILT_SOURCES = depend.mk

#
# States that we compile different fft_INTERFACE_.f90 source
# depending on the engine we choose with --enable-fft-INTERFACE of ./configure
#

INTERFACE=@FFTW_TYPE@

lib2decomp_fft_a_SOURCES += fft_$(INTERFACE).f90

#
# Gives the extra sources to include
# $(MKL_MOD_BIN) could be empty or mkl_dfti.o if we use MKL
#

lib2decomp_fft_a_LIBADD = @MKL_MOD_BIN@ fft_$(INTERFACE).o

#
# Gives the dependencies 
# $(MKL_MOD_BIN) could be empty mkl_dfti.o if we use MKL
#

lib2decomp_fft_a_DEPENDENCIES = @MKL_MOD_BIN@ fft_$(INTERFACE).o

#
# Gives the extra source (from mkl) we might need
# 
#

@MKL_MOD_SOURCE@ :
	- cp $(FFTW_MKL_INCLUDES_PATH)/@MKL_MOD_SOURCE@ .


# Choose pre-processing options
#   -DDOUBLE_PREC   - use double-precision (default single)
#   -DEVEN          - for evenly distributed data, use ALLTOALL
#   -DSHM           - enable system V shared-memory implementation
#   -DOVERWRITE     - allow FFT input to be overwritten (save memory)
#   -DOCC           - use API for overlap communication and computation
#   -DGLOBAL_ARRAYS - turn on support for Global Arrays toolkit
#   -DT3PIO         - turn on LUSTRE IO optimisation code using T3PIO
#

#
# Gives the extra compilation flag we use (Preprocessing, and DOUBLE_PRECISION
#

AM_FCFLAGS = -cpp -DDOUBLE_PREC -I $(FFTW_INCLUDES_PATH)

#
# Local cleaning of *.mod files
#

clean-local:
#-rm -rf *.mod mkl_dfti*
	-rm -rf *.mod @MKL_MOD_BIN@ @MKL_MOD_SOURCE@

#
# include dependencies (Key substituted by configure)
#
 
include depend.mk

#
#
# End of Makefile.am
#
