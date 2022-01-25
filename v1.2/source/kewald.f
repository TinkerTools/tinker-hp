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
      integer nbig,minfft
      integer iefft1,idfft1
      integer iefft2,idfft2
      integer iefft3,idfft3
      integer multi(maxpower)
      integer iproc,ierr,nprocloc
      real*8 delta,rmax
      real*8 edens,ddens
      real*8 size,slope
      real*8 fft1,fft2,fft3
      character*20 keyword
      character*240 record
      character*240 string
 120  format('2d parallel FFT grid : ','Ngrid1 = ',I5,2x,'Ngrid2 = ',I5)
c
      if (use_pmecore) then
        nprocloc = nrec
      else
        nprocloc = nproc
      end if
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
      bseorder = 5
      bsporder = 5
      bsdorder = 4
      edens = 1.2d0
      ddens = 0.8d0
      aeewald = 0.4d0
      apewald = 0.4d0
      adewald = 0.4d0
      minfft = 16
      nlpts = (bsorder-1) /2
      nrpts = bsorder - nlpts - 1
      grdoff = (bsorder+1)/2 +1
c
c     estimate an optimal value for the Ewald coefficient
c
      call ewaldcof (aeewald,ewaldcut)
      if (use_dewald)  call ewaldcof (adewald,dewaldcut)
c
c     modify Ewald coefficient for small unitcell dimensions
c
      if (use_polar) then
         apewald = aeewald
         size = min(xbox,ybox,zbox)
         if (size .lt. 6.0d0) then
            slope = (1.0d0-apewald) / 2.0d0
            apewald = apewald + slope*(6.0d0-size)
            minfft = 64
            if (verbose) then
               write (iout,10)
   10          format (/,' KEWALD  --  Warning, PME Grid Expanded',
     &                    ' due to Small Cell Size')
            end if
         end if
      end if
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
         edens = 0.7d0
         ddens = 0.7d0
         call drivermpi
      end if
c
c     set defaults for electrostatic and dispersion grid sizes
c
      nefft1 = 0
      nefft2 = 0
      nefft3 = 0
      ndfft1 = 0
      ndfft2 = 0
      ndfft3 = 0
c
c     get default grid size from periodic system dimensions
c
      delta = 1.0d-8
      iefft1 = int(xbox*edens-delta) + 1
      iefft2 = int(ybox*edens-delta) + 1
      iefft3 = int(zbox*edens-delta) + 1
      idfft1 = int(xbox*ddens-delta) + 1
      idfft2 = int(ybox*ddens-delta) + 1
      idfft3 = int(zbox*ddens-delta) + 1
      ngrid1 = 0
      ngrid2 = 0
cc
cc     by default, try adapting the 2d grid to the 3d real space domain decomposition
cc
c      if (mod(nprocloc,nydd).eq.0) then
c        ngrid1 = nydd
c        ngrid2 = int(nprocloc/nydd)
c      end if

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
            read (string,*,err=40)  aeewald
         else if (keyword(1:13) .eq. 'PEWALD-ALPHA ') then
            read (string,*,err=40,end=40)  apewald
         else if (keyword(1:13) .eq. 'DEWALD-ALPHA ') then
            read (string,*,err=40,end=40)  adewald
         else if (keyword(1:15) .eq. 'EWALD-BOUNDARY ') then
            boundary = 'VACUUM'
         else if (keyword(1:9) .eq. 'PME-GRID ') then
            fft1 = 0
            fft2 = 0
            fft3 = 0
            read (string,*,err=20,end=20)  fft1,fft2,fft3
   20       continue
            iefft1 = nint(fft1)
            iefft2 = nint(fft2)
            iefft3 = nint(fft3)
            if (iefft2 .eq. 0)  iefft2 = iefft1
            if (iefft3 .eq. 0)  iefft3 = iefft1
         else if (keyword(1:10) .eq. 'DPME-GRID ') then
            fft1 = 0.0d0
            fft2 = 0.0d0
            fft3 = 0.0d0
            read (string,*,err=30,end=30)  fft1,fft2,fft3
   30       continue
            idfft1 = nint(fft1)
            idfft2 = nint(fft2)
            idfft3 = nint(fft3)
            if (idfft2 .eq. 0)  idfft2 = idfft1
            if (idfft3 .eq. 0)  idfft3 = idfft1
         else if (keyword(1:10) .eq. 'PME-ORDER ') then
            read (string,*,err=40,end=40)  bseorder
         else if (keyword(1:11) .eq. 'DPME-ORDER ') then
            read (string,*,err=40,end=40)  bsdorder
         else if (keyword(1:11) .eq. 'PPME-ORDER ') then
            read (string,*,err=40,end=40)  bsporder
         else if (keyword(1:15) .eq. 'DECOMPFFT-GRID ') then
            read (string,*,err=40) ngrid1,ngrid2
         else if (keyword(1:15) .eq. 'DECOMPFFT-OPTIM ') then
c
c        keyword to let the 2decomp library optimize the 2d processor grid for parallel fft
c
           ngrid1 = 0
           ngrid2 = 0
         end if
   40    continue
      end do

c
c     grid size must be even, with prime factors of 2, 3 and 5
c
      nefft1 = maxfft
      nefft2 = maxfft
      nefft3 = maxfft
      do i = maxpower, 1, -1
         k = multi(i)
         if (k .le. maxfft) then
            if (k .ge. iefft1)  nefft1 = k
            if (k .ge. iefft2)  nefft2 = k
            if (k .ge. iefft3)  nefft3 = k
         end if
         if (nefft1 .lt. minfft)  nefft1 = minfft
         if (nefft2 .lt. minfft)  nefft2 = minfft
         if (nefft3 .lt. minfft)  nefft3 = minfft
      end do
c
c     determine dispersion grid size from allowed values
c
      if (use_dewald) then
         ndfft1 = maxfft
         ndfft2 = maxfft
         ndfft3 = maxfft
         do i = maxpower, 1, -1
            k = multi(i)
            if (k .le. maxfft) then
               if (k .ge. idfft1)  ndfft1 = k
               if (k .ge. idfft2)  ndfft2 = k
               if (k .ge. idfft3)  ndfft3 = k
            end if
         end do
         if (ndfft1 .lt. minfft)  ndfft1 = minfft
         if (ndfft2 .lt. minfft)  ndfft2 = minfft
         if (ndfft3 .lt. minfft)  ndfft3 = minfft
      end if
c
c     check the particle mesh Ewald grid dimensions
c
      nbig = max(nefft1,nefft2,nefft3,ndfft1,ndfft2,ndfft3)
      if (nbig .gt. maxfft) then
         write (iout,50)
   50    format (/,' KEWALD  --  PME Grid Size Too Large;',
     &              ' Increase MAXFFT')
         call fatal
      end if
      if (nefft1.lt.iefft1.or.
     &       nefft2.lt.iefft2.or.nefft3.lt.iefft3) then
         write (iout,60)
   60    format (/,' KEWALD  --  Warning, Small Electrostatic',
     &              'PME Grid Size')
      end if
      if (use_dewald .and. (ndfft1.lt.idfft1.or.
     &       ndfft2.lt.idfft2.or.ndfft3.lt.idfft3)) then
         write (iout,70)
   70    format (/,' KEWALD  --  Warning, Small Dispersion',
     &              'PME Grid Size')
      end if
c
c     set maximum sizes for PME grid and B-spline order
c
c
c     enforce use of max pme grid in all PME-based computations
c
      nfft1 = max(nefft1,ndfft1)
      nfft2 = max(nefft2,ndfft2)
      nfft3 = max(nefft3,ndfft3)
c
c     enforce use of max spline order in all PME-based computations
c
      bsorder = max(bseorder,bsporder,bsdorder)

      nlpts = (bsorder-1) /2
      nrpts = bsorder - nlpts - 1
      grdoff = (bsorder+1)/2 +1
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
      if (ngrid1.ne.0 .and. ngrid2.ne.0) then
        if ((use_pmecore).and.(rank.eq.ndir)) then
          write(iout,120) ngrid1,ngrid2
        else if (.not.(use_pmecore).and.(rank.eq.0)) then
          write(iout,120) ngrid1,ngrid2
        end if
      end if
c
       call reassignpme(.true.)
c
c     print a message listing some of the Ewald parameters
c
      if (verbose.and.(rank.eq.0)) then
         write (iout,80)
   80    format (/,' Particle Mesh Ewald Parameters :',
     &           //,5x,'Type',16x,'Ewald Alpha',4x,'Grid',
     &              ' Dimensions',4x,'Spline Order',/)
         write (iout,90)  aeewald,nfft1,nfft2,nfft3,bsorder
   90    format (3x,'Electrostatics',9x,f8.4,5x,3i5,7x,i5)
         if (use_polar) then
            write (iout,100)  apewald,nfft1,nfft2,nfft3,bsorder
  100       format (3x,'Polarization',11x,f8.4,5x,3i5,7x,i5)
         end if
         if (use_dewald) then
            write (iout,110)  adewald,nfft1,nfft2,nfft3,bsorder
  110       format (3x,'Dispersion',13x,f8.4,5x,3i5,7x,i5)
         end if
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
