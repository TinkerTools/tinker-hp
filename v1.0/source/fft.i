c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ################################################################
c     ##                                                            ##
c     ##  fft.i  --  values and options for Fast Fourier transform  ##
c     ##                                                            ##
c     ################################################################
c
c
c     maxtable   maximum size of the FFT table intermediate array
c     maxprime   maximum number of prime factors of FFT dimension
c
      integer maxtable
      integer maxprime
      integer n1mpimax,n2mpimax,n3mpimax
      parameter (maxtable=4*maxfft)
      parameter (maxprime=15)
      integer iprime
      integer ngrid1,ngrid2
      integer, pointer :: istart1(:), iend1(:)
      integer, pointer :: jstart1(:), jend1(:)
      integer, pointer :: kstart1(:), kend1(:)
      integer, pointer :: isize1(:),jsize1(:),ksize1(:)
      integer, pointer :: istart2(:), iend2(:)
      integer, pointer :: jstart2(:), jend2(:)
      integer, pointer :: kstart2(:), kend2(:)
      integer, pointer :: isize2(:),jsize2(:),ksize2(:)
      real*8 ffttable
      character*7 ffttyp
      common /fft/ ffttyp,
     $             istart1,jstart1,kstart1,
     $             iend1,jend1,kend1,
     $             isize1,jsize1,ksize1,
     $             istart2,jstart2,kstart2,
     $             iend2,jend2,kend2,
     $             isize2,jsize2,ksize2,n1mpimax,n2mpimax,n3mpimax,
     $             ngrid1,ngrid2
