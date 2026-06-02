This folder contains two sample applications to compute multiple independent 
FFTs. The first application - 'blocking.f90' - uses the standard blocking
version of MPI communication code to transpose the data among different stages
of the computation. The second - 'non_blocking.f90' - performs the same 
computation using the non-blocking communication code supplied by 2DECOMP&FFT. 

Please note that the non-blocking communication was realised using libNBC,
a library implementing non-blocking MPI collectives (such as IALLTOALL) using 
existing MPI 1 functions. Such functionality is expected to be part of the
MPI 3 standard. 

To build the non-blocking test application, user needs to obtain a copy of 
libNBC at http://www.unixer.de/research/nbcoll/libnbc/, compile it, and link it
to 2DECOMP&FFT and the object code of the test application.

To demonstrate the idea of overlap communication and computation, the 3D FFT 
is implemented using multiple 1D FFTs (without using the advanced interface of 
FFTW) so that MPI_TEST calls can be easily inserted in the computational part 
of the code. This is required because the communication has to be explicitly 
progressed when running on the same thread as the computation.
