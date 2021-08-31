c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine kewald  --  Ewald sum parameter assignment  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "kewald" assigns particle mesh Ewald parameters and options
c     for a periodic system
c
c
      subroutine kewald
      use atoms
      use bound
      use boxes
      use chunks
      use cutoff
      use domdec
      use ewald
      use fft
      use inform
      use iounit
      use keys
      use pme
      use potent
      use mpi
      implicit none
      integer maxpower
      parameter (maxpower=63)
      integer i,k,next
      integer ifft1,ifft2,ifft3
      integer multi(maxpower)
      integer iproc,ierr
      real*8 delta,rmax,dens
      character*20 keyword
      character*240 record
      character*240 string
c
c     allocate global arrays
c
      if (allocated(igrid)) deallocate (igrid)
      allocate (igrid(3,n))
c
c     PME grid size must be even with factors of only 2, 3 and 5
c
      data multi  /   2,   4,   6,   8,  10,  12,  16,  18,  20,
     &               24,  30,  32,  36,  40,  48,  50,  54,  60,
     &               64,  72,  80,  90,  96, 100, 108, 120, 128,
     &              144, 150, 160, 162, 180, 192, 200, 216, 240,
     &              250, 256, 270, 288, 300, 320, 324, 360, 384,
     &              400, 432, 450, 480, 486, 500, 512, 540, 576,
     &              600, 640, 648, 720, 750, 768, 800, 810, 864 /
c
c
c     default boundary treatment, B-spline order and grid density
c
      boundary = 'TINFOIL'
      bsorder = 5
      dens = 1.2d0
c
c     estimate an optimal value for the Ewald coefficient
c
      call ewaldcof (aewald,ewaldcut)
c
c     set the system extent for nonperiodic Ewald summation
c
      if (.not. use_bounds) then
         call extent (rmax)
         xbox = 2.0d0 * (rmax+ewaldcut)
         ybox = xbox
         zbox = xbox
         alpha = 90.0d0
         beta = 90.0d0
         gamma = 90.0d0
         orthogonal = .true.
         call lattice
         boundary = 'NONE'
         dens = 0.75d0
      end if
c
c     get default grid size from periodic system dimensions
c
      delta = 1.0d-8
      ifft1 = int(xbox*dens-delta) + 1
      ifft2 = int(ybox*dens-delta) + 1
      ifft3 = int(zbox*dens-delta) + 1
      ngrid1 = 0
      ngrid2 = 0
c
c     search keywords for Ewald summation commands
c
      do i = 1, nkey
         record = keyline(i)
         next = 1
         call upcase (record)
         call gettext (record,keyword,next)
         string = record(next:240)
         if (keyword(1:12) .eq. 'EWALD-ALPHA ') then
            read (string,*,err=20)  aewald
         else if (keyword(1:15) .eq. 'EWALD-BOUNDARY ') then
            boundary = 'VACUUM'
         else if (keyword(1:9) .eq. 'PME-GRID ') then
            ifft1 = 0
            ifft2 = 0
            ifft3 = 0
            read (string,*,err=10,end=10)  ifft1,ifft2,ifft3
   10       continue
            if (ifft2 .eq. 0)  ifft2 = ifft1
            if (ifft3 .eq. 0)  ifft3 = ifft1
         else if (keyword(1:10) .eq. 'PME-ORDER ') then
            read (string,*,err=20)  bsorder
         else if (keyword(1:15) .eq. 'DECOMPFFT-GRID ') then
            read (string,*,err=20) ngrid1,ngrid2
         end if
   20    continue
      end do

      nlpts = (bsorder-1) /2
      nrpts = bsorder - nlpts - 1
      grdoff = (bsorder+1)/2 +1
c
c     grid size must be even, with prime factors of 2, 3 and 5
c
      nfft1 = maxfft
      nfft2 = maxfft
      nfft3 = maxfft
      do i = maxpower, 1, -1
         k = multi(i)
         if (k .le. maxfft) then
            if (k .ge. ifft1)  nfft1 = k
            if (k .ge. ifft2)  nfft2 = k
            if (k .ge. ifft3)  nfft3 = k
         end if
      end do
c
c     check the B-spline order and charge grid dimension
c
      if (bsorder .gt. maxorder) then
         if (rank.eq.0) write (iout,30)
   30    format (/,' KEWALD  --  B-Spline Order Too Large;',
     &              ' Increase MAXORDER')
         call fatal
      end if
      if (max(nfft1,nfft2,nfft3) .gt. maxfft) then
         write (iout,40)
   40    format (/,' KEWALD  --  FFT Charge Grid Too Large;',
     &              ' Increase MAXFFT')
         call fatal
      else if (nfft1.lt.ifft1 .or. nfft2.lt.ifft2
     &             .or. nfft3.lt.ifft3) then
         if (rank.eq.0) write (iout,50)
   50    format (/,' KEWALD  --  Warning, Small Charge Grid',
     &              ' may give Poor Accuracy')
      end if
c
c     perform dynamic allocation of some pointer arrays
c
      if (allocated(thetai1)) call fft_final
      if (allocated(thetai1))  deallocate (thetai1)
      if (allocated(thetai2))  deallocate (thetai2)
      if (allocated(thetai3))  deallocate (thetai3)
c
      if (allocated(qgridin_2d))  deallocate (qgridin_2d)
      if (allocated(qgrid2in_2d))  deallocate (qgrid2in_2d)
      if (allocated(qgridout_2d))  deallocate (qgridout_2d)
      if (allocated(qgrid2out_2d))  deallocate (qgrid2out_2d)
      if (allocated(qfac_2d))  deallocate (qfac_2d)
      if (allocated(istart1)) deallocate (istart1)
      if (allocated(iend1)) deallocate (iend1)
      if (allocated(isize1)) deallocate (isize1)
      if (allocated(jstart1)) deallocate (jstart1)
      if (allocated(jend1)) deallocate (jend1)
      if (allocated(jsize1)) deallocate (jsize1)
      if (allocated(kstart1)) deallocate (kstart1)
      if (allocated(kend1)) deallocate (kend1)
      if (allocated(ksize1)) deallocate (ksize1)
      if (allocated(istart2)) deallocate (istart2)
      if (allocated(iend2)) deallocate (iend2)
      if (allocated(isize2)) deallocate (isize2)
      if (allocated(jstart2)) deallocate (jstart2)
      if (allocated(jend2)) deallocate (jend2)
      if (allocated(jsize2)) deallocate (jsize2)
      if (allocated(kstart2)) deallocate (kstart2)
      if (allocated(kend2)) deallocate (kend2)
c
c     also deallocate 2decomp stuff
c
      if (allocated(ksize2)) deallocate (ksize2)
c
      if (use_pmecore) then
        allocate (istart1(nrec))
        allocate (iend1(nrec))
        allocate (isize1(nrec))
        allocate (jstart1(nrec))
        allocate (jend1(nrec))
        allocate (jsize1(nrec))
        allocate (kstart1(nrec))
        allocate (kend1(nrec))
        allocate (ksize1(nrec))
c
        allocate (istart2(nrec))
        allocate (iend2(nrec))
        allocate (isize2(nrec))
        allocate (jstart2(nrec))
        allocate (jend2(nrec))
        allocate (jsize2(nrec))
        allocate (kstart2(nrec))
        allocate (kend2(nrec))
        allocate (ksize2(nrec))
      else
        allocate (istart1(nproc))
        allocate (iend1(nproc))
        allocate (isize1(nproc))
        allocate (jstart1(nproc))
        allocate (jend1(nproc))
        allocate (jsize1(nproc))
        allocate (kstart1(nproc))
        allocate (kend1(nproc))
        allocate (ksize1(nproc))
c
        allocate (istart2(nproc))
        allocate (iend2(nproc))
        allocate (isize2(nproc))
        allocate (jstart2(nproc))
        allocate (jend2(nproc))
        allocate (jsize2(nproc))
        allocate (kstart2(nproc))
        allocate (kend2(nproc))
        allocate (ksize2(nproc))
c
      end if
c
c     initialize the PME arrays that can be precomputed
c
      call moduli
      n1mpimax = 0
      n2mpimax = 0
      n3mpimax = 0
      if (use_pmecore) then
        if (rank.gt.ndir-1) then
          call fft_setup(nrec,rank_bis,istart1,iend1,isize1,
     $      jstart1,jend1,jsize1,kstart1,kend1,ksize1,
     $      istart2,iend2,isize2,jstart2,jend2,jsize2,
     $      kstart2,kend2,ksize2,nfft1,nfft2,nfft3,comm_rec,ngrid1,
     $      ngrid2)
        end if
        do iproc = 0, nrec-1
          call MPI_BCAST(istart1(iproc+1),1,MPI_INT,iproc+ndir,
     $     COMM_TINKER,ierr)
          call MPI_BCAST(iend1(iproc+1),1,MPI_INT,iproc+ndir,
     $     COMM_TINKER,ierr)
          call MPI_BCAST(isize1(iproc+1),1,MPI_INT,iproc+ndir,
     $     COMM_TINKER,ierr)
          call MPI_BCAST(jstart1(iproc+1),1,MPI_INT,iproc+ndir,
     $     COMM_TINKER,ierr)
          call MPI_BCAST(jend1(iproc+1),1,MPI_INT,iproc+ndir,
     $     COMM_TINKER,ierr)
          call MPI_BCAST(jsize1(iproc+1),1,MPI_INT,iproc+ndir,
     $     COMM_TINKER,ierr)
          call MPI_BCAST(kstart1(iproc+1),1,MPI_INT,iproc+ndir,
     $     COMM_TINKER,ierr)
          call MPI_BCAST(kend1(iproc+1),1,MPI_INT,iproc+ndir,
     $     COMM_TINKER,ierr)
          call MPI_BCAST(ksize1(iproc+1),1,MPI_INT,iproc+ndir,
     $     COMM_TINKER,ierr)
c
          call MPI_BCAST(istart2(iproc+1),1,MPI_INT,iproc+ndir,
     $     COMM_TINKER,ierr)
          call MPI_BCAST(iend2(iproc+1),1,MPI_INT,iproc+ndir,
     $     COMM_TINKER,ierr)
          call MPI_BCAST(isize2(iproc+1),1,MPI_INT,iproc+ndir,
     $     COMM_TINKER,ierr)
          call MPI_BCAST(jstart2(iproc+1),1,MPI_INT,iproc+ndir,
     $     COMM_TINKER,ierr)
          call MPI_BCAST(jend2(iproc+1),1,MPI_INT,iproc+ndir,
     $     COMM_TINKER,ierr)
          call MPI_BCAST(jsize2(iproc+1),1,MPI_INT,iproc+ndir,
     $     COMM_TINKER,ierr)
          call MPI_BCAST(kstart2(iproc+1),1,MPI_INT,iproc+ndir,
     $     COMM_TINKER,ierr)
          call MPI_BCAST(kend2(iproc+1),1,MPI_INT,iproc+ndir,
     $     COMM_TINKER,ierr)
          call MPI_BCAST(ksize2(iproc+1),1,MPI_INT,iproc+ndir,
     $     COMM_TINKER,ierr)
        end do
c
c    get maximum size of the grid of the neighbors
c
        do iproc = 0, nrec-1
          if (isize1(iproc+1).ge.n1mpimax) then
            n1mpimax = isize1(iproc+1)
          end if
          if (jsize1(iproc+1).ge.n2mpimax) then
            n2mpimax = jsize1(iproc+1)
          end if
          if (ksize1(iproc+1).ge.n3mpimax) then
            n3mpimax = ksize1(iproc+1)
          end if
        end do

        if (allocated (qfac_2d)) deallocate (qfac_2d)
        if (rank.gt.ndir-1) then
          allocate (qfac_2d(isize2(rank_bis+1),jsize2(rank_bis+1),
     $        ksize2(rank_bis+1)))
        end if
      else
        call fft_setup(nproc,rank,istart1,iend1,isize1,
     $   jstart1,jend1,jsize1,kstart1,kend1,ksize1,
     $   istart2,iend2,isize2,jstart2,jend2,jsize2,
     $   kstart2,kend2,ksize2,nfft1,nfft2,nfft3,COMM_TINKER,
     $   ngrid1,ngrid2)
        do iproc = 0, nproc-1
          call MPI_BCAST(istart1(iproc+1),1,MPI_INT,iproc,
     $     COMM_TINKER,ierr)
          call MPI_BCAST(iend1(iproc+1),1,MPI_INT,iproc,
     $     COMM_TINKER,ierr)
          call MPI_BCAST(isize1(iproc+1),1,MPI_INT,iproc,
     $     COMM_TINKER,ierr)
          call MPI_BCAST(jstart1(iproc+1),1,MPI_INT,iproc,
     $     COMM_TINKER,ierr)
          call MPI_BCAST(jend1(iproc+1),1,MPI_INT,iproc,
     $     COMM_TINKER,ierr)
          call MPI_BCAST(jsize1(iproc+1),1,MPI_INT,iproc,
     $     COMM_TINKER,ierr)
          call MPI_BCAST(kstart1(iproc+1),1,MPI_INT,iproc,
     $     COMM_TINKER,ierr)
          call MPI_BCAST(kend1(iproc+1),1,MPI_INT,iproc,
     $     COMM_TINKER,ierr)
          call MPI_BCAST(ksize1(iproc+1),1,MPI_INT,iproc,
     $     COMM_TINKER,ierr)
c
          call MPI_BCAST(istart2(iproc+1),1,MPI_INT,iproc,
     $     COMM_TINKER,ierr)
          call MPI_BCAST(iend2(iproc+1),1,MPI_INT,iproc,
     $     COMM_TINKER,ierr)
          call MPI_BCAST(isize2(iproc+1),1,MPI_INT,iproc,
     $     COMM_TINKER,ierr)
          call MPI_BCAST(jstart2(iproc+1),1,MPI_INT,iproc,
     $     COMM_TINKER,ierr)
          call MPI_BCAST(jend2(iproc+1),1,MPI_INT,iproc,
     $     COMM_TINKER,ierr)
          call MPI_BCAST(jsize2(iproc+1),1,MPI_INT,iproc,
     $     COMM_TINKER,ierr)
          call MPI_BCAST(kstart2(iproc+1),1,MPI_INT,iproc,
     $     COMM_TINKER,ierr)
          call MPI_BCAST(kend2(iproc+1),1,MPI_INT,iproc,
     $     COMM_TINKER,ierr)
          call MPI_BCAST(ksize2(iproc+1),1,MPI_INT,iproc,
     $     COMM_TINKER,ierr)
        end do
c
c    get maximum size of the grid of the neighbors
c
        do iproc = 0, nproc-1
          if (isize1(iproc+1).ge.n1mpimax) then
            n1mpimax = isize1(iproc+1)
          end if
          if (jsize1(iproc+1).ge.n2mpimax) then
            n2mpimax = jsize1(iproc+1)
          end if
          if (ksize1(iproc+1).ge.n3mpimax) then
            n3mpimax = ksize1(iproc+1)
          end if
        end do
c
        allocate (qfac_2d(isize2(rank+1),jsize2(rank+1),ksize2(rank+1)))
      end if
c
       call reassignpme(.true.)
c
c     print a message listing some of the Ewald parameters
c
      if (verbose.and.(rank.eq.0)) then
         write (iout,60)  aewald,nfft1,nfft2,nfft3,bsorder
   60    format (/,' Smooth Particle Mesh Ewald Parameters :',
     &           //,4x,'Ewald Coefficient',6x,'Charge Grid',
     &              ' Dimensions',6x,'B-Spline Order',
     &           //,8x,f8.4,11x,3i6,12x,i6)
      end if
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine ewaldcof  --  estimation of Ewald coefficient  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "ewaldcof" finds an Ewald coefficient such that all terms
c     beyond the specified cutoff distance will have a value less
c     than a specified tolerance
c
c
      subroutine ewaldcof (alpha,cutoff)
      implicit none
      integer i,k
      real*8 alpha,cutoff,eps
      real*8 x,xlo,xhi,y
      real*8 ratio,erfc
      external erfc
c
c
c     set the tolerance value; use of 1.0d-8 results in large
c     Ewald coefficients that ensure continuity in the gradient
c
      eps = 1.0d-8
c
c     get approximate value from cutoff and tolerance
c
      ratio = eps + 1.0d0
      x = 0.5d0
      i = 0
      do while (ratio .ge. eps)
         i = i + 1
         x = 2.0d0 * x
         y = x * cutoff
         ratio = erfc(y) / cutoff
      end do
c
c     use a binary search to refine the coefficient
c
      k = i + 60
      xlo = 0.0d0
      xhi = x
      do i = 1, k
         x = (xlo+xhi) / 2.0d0
         y = x * cutoff
         ratio = erfc(y) / cutoff
         if (ratio .ge. eps) then
            xlo = x
         else
            xhi = x
         end if
      end do
      alpha = x
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine extent  --  find maximum interatomic distance  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "extent" finds the largest interatomic distance in a system
c
c
      subroutine extent (rmax)
      use sizes
      use atoms
      implicit none
      integer i,k
      real*8 xi,yi,zi
      real*8 xk,yk,zk
      real*8 r2,rmax
c
c
c     search all atom pairs to find the largest distance
c
      rmax = 0.0d0
      do i = 1, n-1
         xi = x(i)
         yi = y(i)
         zi = z(i)
         do k = i+1, n
            xk = x(k)
            yk = y(k)
            zk = z(k)
            r2 = (xi-xk)**2 + (yi-yk)**2 + (zi-zk)**2
            rmax = max(r2,rmax)
         end do
      end do
      rmax = sqrt(rmax)
      return
      end
c
c
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine moduli  --  store the inverse DFT moduli  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     "moduli" sets the moduli of the inverse discrete Fourier
c     transform of the B-splines
c
c
      subroutine moduli
      use pme
      implicit none
      integer i
      real*8 x
      real*8 array(maxorder)
      real*8 bsarray(maxfft)
c
c
c     compute and load the moduli values
c
      x = 0.0d0
      call bspline (x,bsorder,array)
      do i = 1, maxfft
         bsarray(i) = 0.0d0
      end do
      do i = 2, bsorder+1
         bsarray(i) = array(i-1)
      end do
      call dftmod (bsmod1,bsarray,nfft1,bsorder)
      call dftmod (bsmod2,bsarray,nfft2,bsorder)
      call dftmod (bsmod3,bsarray,nfft3,bsorder)
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine bspline  --  determine B-spline coefficients  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "bspline" calculates the coefficients for an n-th order
c     B-spline approximation
c
c
      subroutine bspline (x,n,c)
      implicit none
      integer i,k,n
      real*8 x,denom
      real*8 c(*)
c
c
c     initialize the B-spline as the linear case
c
      c(1) = 1.0d0 - x
      c(2) = x
c
c     compute standard B-spline recursion to n-th order
c
      do k = 3, n
         denom = 1.0d0 / dble(k-1)
         c(k) = x * c(k-1) * denom
         do i = 1, k-2
            c(k-i) = ((x+dble(i))*c(k-i-1)
     &                  + (dble(k-i)-x)*c(k-i)) * denom
         end do
         c(1) = (1.0d0-x) * c(1) * denom
      end do
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine dftmod  --  discrete Fourier transform modulus  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "dftmod" computes the modulus of the discrete Fourier transform
c     of "bsarray" and stores it in "bsmod"
c
c
      subroutine dftmod (bsmod,bsarray,nfft,order)
      use math
      implicit none
      integer i,j,k
      integer nfft,jcut
      integer order,order2
      real*8 eps,zeta
      real*8 arg,factor
      real*8 sum1,sum2
      real*8 bsmod(*)
      real*8 bsarray(*)
c
c
c     get the modulus of the discrete Fourier transform
c
      factor = 2.0d0 * pi / dble(nfft)
      do i = 1, nfft
         sum1 = 0.0d0
         sum2 = 0.0d0
         do j = 1, nfft
            arg = factor * dble((i-1)*(j-1))
            sum1 = sum1 + bsarray(j)*cos(arg)
            sum2 = sum2 + bsarray(j)*sin(arg)
         end do
         bsmod(i) = sum1**2 + sum2**2
      end do
c
c     fix for exponential Euler spline interpolation failure
c
      eps = 1.0d-7
      if (bsmod(1) .lt. eps)  bsmod(1) = 0.5d0 * bsmod(2)
      do i = 2, nfft-1
         if (bsmod(i) .lt. eps)
     &      bsmod(i) = 0.5d0 * (bsmod(i-1)+bsmod(i+1))
      end do
      if (bsmod(nfft) .lt. eps)  bsmod(nfft) = 0.5d0 * bsmod(nfft-1)
c
c     compute and apply the optimal zeta coefficient
c
      jcut = 50
      order2 = 2 * order
      do i = 1, nfft
         k = i - 1
         if (i .gt. nfft/2)  k = k - nfft
         if (k .eq. 0) then
            zeta = 1.0d0
         else
            sum1 = 1.0d0
            sum2 = 1.0d0
            factor = pi * dble(k) / dble(nfft)
            do j = 1, jcut
               arg = factor / (factor+pi*dble(j))
               sum1 = sum1 + arg**order
               sum2 = sum2 + arg**order2
            end do
            do j = 1, jcut
               arg = factor / (factor-pi*dble(j))
               sum1 = sum1 + arg**order
               sum2 = sum2 + arg**order2
            end do
            zeta = sum2 / sum1
         end if
         bsmod(i) = bsmod(i) * zeta**2
      end do
      return
      end
