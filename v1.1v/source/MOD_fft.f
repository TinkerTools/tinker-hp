c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #####################################################################
c     ##                                                                 ##
c     ##  module fft  --  values and options for Fast Fourier transform  ##
c     ##                                                                 ##
c     #####################################################################
c
c
c     maxtable   maximum size of the FFT table intermediate array
c     maxprime   maximum number of prime factors of FFT dimension
c     n1mpimax   maximum number of points along the divided first dimension of the fft grid 
c     n2mpimax   maximum number of points along the divided second dimension of the fft grid 
c     n3mpimax   maximum number of points along the divided third dimension of the fft grid 
c     ngrid1     first dimension of the 2d proc grid for the 2decomp library (2d pencils decomposition)
c     ngrid2     second dimension of the 2d proc grid for the 2decomp library (2d pencils decomposition)
c     istart1    list of beginning indexes of the first axis of the grid among the procs 
c     jstart1    list of beginning indexes of the second axis of the grid among the procs
c     kstart1    list of beginning indexes of the third axis of the grid among the procs
c     iend1    list of ending indexes of the first axis of the grid among the procs
c     jend1    list of ending indexes of the second axis of the grid among the procs
c     kend1    list of ending indexes of the third axis of the grid among the procs 
c     isize1   list of size of the first axis of the grid among the procs
c     jsize1   list of size of the second axis of the grid among the procs
c     ksize1   list of size of the third axis of the grid among the procs
c     istart2    list of beginning indexes of the first axis of the grid among the procs for the transposed grid
c     jstart2    list of beginning indexes of the second axis of the grid among the procs for the transposed grid
c     kstart2    list of beginning indexes of the third axis of the grid among the procs for the transposed grid
c     iend2    list of ending indexes of the first axis of the grid among the procs for the transposed grid
c     jend2    list of ending indexes of the second axis of the grid among the procs for the transposed grid
c     kend2    list of ending indexes of the third axis of the grid among the procs for the transposed grid
c     isize2   list of size of the first axis of the transposed grid among the procs
c     jsize2   list of size of the second axis of the transposed grid among the procs
c     ksize2   list of size of the third axis of the transposed grid among the procs
c
      module fft
      use sizes
      implicit none
      integer maxtable
      integer maxprime
      integer n1mpimax,n2mpimax,n3mpimax
      parameter (maxtable=4*maxfft)
      parameter (maxprime=15)
      integer ngrid1,ngrid2
      !DIR$ ATTRIBUTES ALIGN:64:: istart1, iend1
      integer, allocatable :: istart1(:), iend1(:)
      !DIR$ ATTRIBUTES ALIGN:64:: jstart1, jend1
      integer, allocatable :: jstart1(:), jend1(:)
      !DIR$ ATTRIBUTES ALIGN:64:: kstart1, kend1
      integer, allocatable :: kstart1(:), kend1(:)
      !DIR$ ATTRIBUTES ALIGN:64:: isize1 ,jsize1,ksize1
      integer, allocatable :: isize1(:),jsize1(:),ksize1(:)
      !DIR$ ATTRIBUTES ALIGN:64:: istart2, iend2
      integer, allocatable :: istart2(:), iend2(:)
      !DIR$ ATTRIBUTES ALIGN:64:: jstart2, jend2
      integer, allocatable :: jstart2(:), jend2(:)
      !DIR$ ATTRIBUTES ALIGN:64:: kstart2, kend2
      integer, allocatable :: kstart2(:), kend2(:)
      !DIR$ ATTRIBUTES ALIGN:64:: isize2 ,jsize2,ksize2
      integer, allocatable :: isize2(:),jsize2(:),ksize2(:)
      save
      end
