c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ##################################################################
c     ##                                                              ##
c     ##  program minimize  --  low storage BFGS Cartesian optimizer  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "minimize" performs energy minimization in Cartesian coordinate
c     space using a low storage BFGS nonlinear optimization
c
c
#include "tinker_precision.h"
      program minimize
      use atoms
      use domdec
      use deriv ,only:info_forces,cDef
      use energi,only:info_energy
      use files
      use keys
      use inform
      use iounit
      use interfaces,only: LBFGS
      use scales
      use tinMemory ,only: prmem_requestm
      use usage
      use utils     ,only:set_to_zero1m
      use utilgpu   ,only:rec_queue
      use mpi
      implicit none
      integer i,j,imin,ierr
      integer next,freeunit,nthreadsupport
      integer iglob
      real(r_p) minimum,minimiz1
      real(r_p) grdmin
      real(r_p) grms
      real(r_p) gnorm
      real(r_p) energy
      real(r_p) eps
      real(r_p), allocatable :: xx(:)
      real(r_p), allocatable :: derivs(:,:)
      logical exist,analytic
      character*20 keyword
      character*120 minfile
      character*120 record
      character*120 string
      external energy
      external minimiz1
      external optsave
      call MPI_INIT_THREAD(MPI_THREAD_MULTIPLE,nthreadsupport,ierr)
c
c
c     set up the structure and mechanics calculation
c
      call initial
      call initmpi
      call getxyz
      call cutoffs
      call unitcell
      call lattice
c
c     setup for MPI
c
      call drivermpi
      call reinitnl(0)
c
      call mechanic
      call nblist(0)
      call allocstep
c
      if (associated(scale)) then
!$acc exit data delete(scale)
         deallocate (scale)
      end if
      allocate (scale(3*n))
!$acc enter data create(scale)
c 
c     use either analytical or numerical gradients
c
      analytic = .true.
      eps = 0.00001_re_p
c
c     search the keywords for output frequency parameters
c
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         string = record(next:120)
         if (keyword(1:9) .eq. 'PRINTOUT ') then
            read (string,*,err=10,end=10)  iprint
         else if (keyword(1:9) .eq. 'WRITEOUT ') then
            read (string,*,err=10,end=10)  iwrite
         end if
   10    continue
      end do
c
c     get termination criterion as RMS gradient per atom
c
      grdmin = -1.0_re_p
      call nextarg (string,exist)
      if (exist)  read (string,*,err=20,end=20)  grdmin
   20 continue
      if (grdmin .le. 0.0_re_p) then
         write (iout,30)
   30    format (/,' Enter RMS Gradient per Atom Criterion',
     &              ' [0.01] :  ',$)
         read (input,40)  grdmin
   40    format (f20.0)
      end if
      if (grdmin .le. 0.0_re_p)  grdmin = 0.01_re_p
c
c     write out a copy of coordinates for later update
c
      imin = freeunit ()
      minfile = filename(1:leng)//'.xyz'
      call version (minfile,'new')
      if (rank.eq.0) then
        open (unit=imin,file=minfile,status='new')
        call prtxyz (imin)
        close (unit=imin)
        outfile = minfile
      end if
c
c     set scaling parameter for function and derivative values;
c     use square root of median eigenvalue of typical Hessian
c
      set_scale = .true.
!$acc parallel loop async present(scale)
      do i = 1, 3*n
         scale(i) = 12.0_re_p
      end do
c
c     perform dynamic allocation of some local arrays
c
      allocate (xx(3*n))
!$acc enter data create(xx) async
c
c     scale the coordinates of each active atom
c
!$acc parallel loop async default(present)
      do i = 1, nloc
         iglob = glob(i)
         if (use(iglob)) then
            xx(3*(iglob-1)+1) = x(iglob) * scale(3*(iglob-1)+1)
            xx(3*(iglob-1)+2) = y(iglob) * scale(3*(iglob-1)+2)
            xx(3*(iglob-1)+3) = z(iglob) * scale(3*(iglob-1)+3)
         end if
      end do
c
c     make the call to the optimization routine
c
      call lbfgs (n,xx,minimum,grdmin,minimiz1,optsave)
c
c     unscale the final coordinates for active atoms
c
!$acc parallel loop async default(present)
      do i = 1, nloc
         iglob = glob(i)
         if (use(iglob)) then
            x(iglob) = xx(3*(iglob-1)+1) / scale(3*(iglob-1)+1)
            y(iglob) = xx(3*(iglob-1)+2) / scale(3*(iglob-1)+2)
            z(iglob) = xx(3*(iglob-1)+3) / scale(3*(iglob-1)+3)
         end if
      end do
c
c     compute the final function and RMS gradient values
c
c      if (analytic) then
c         call commstep
         call sendallpos
         call AllDirAssign
         call reassignpme(.false.)
         allocate (derivs(3,nbloc))
!$acc enter data create(derivs,minimum) async
         call set_to_zero1m(derivs,size(derivs),rec_queue)
         call reinitnl(0)
         call mechanicstep(0)
         call allocstep
         call nblist(0)
         call gradient (minimum,derivs)
!$acc wait
!$acc update host(minimum)
         call MPI_ALLREDUCE(MPI_IN_PLACE,minimum,1,MPI_RPREC,
     $        MPI_SUM,MPI_COMM_WORLD,ierr)
         call commforces(derivs)
c      else
c         minimum = energy ()
cc         call numgrad (energy,derivs,eps)
c      end if

      !Debug
      if(deb_Force) call info_forces(cDef)
      if(deb_Energy)call info_energy(rank)

      gnorm = 0.0_re_p
!$acc parallel loop async present(derivs)
      do i = 1, nloc
         do j = 1, 3
            iglob = glob(i)
            if (use(iglob)) then
               gnorm = gnorm + derivs(j,i)**2
            end if
         end do
      end do
      call MPI_ALLREDUCE(MPI_IN_PLACE,gnorm,1,MPI_RPREC,
     $     MPI_SUM,MPI_COMM_WORLD,ierr)
      gnorm = sqrt(gnorm)
      grms = gnorm / sqrt(real(3*n,t_p))
c
c     perform deallocation of some local arrays
c
!$acc exit data delete(xx,derivs,minimum) async
      deallocate (xx)
      deallocate (derivs)
c
c     write out the final function and gradient values
c
      if (rank.eq.0) then
      if (digits .ge. 8) then
         if (grms .gt. 1.0d-8) then
            write (iout,50)  minimum,grms,gnorm
   50       format (/,' Final Function Value :',2x,f20.8,
     &              /,' Final RMS Gradient :',4x,f20.8,
     &              /,' Final Gradient Norm :',3x,f20.8)
         else
            write (iout,60)  minimum,grms,gnorm
   60       format (/,' Final Function Value :',2x,f20.8,
     &              /,' Final RMS Gradient :',4x,d20.8,
     &              /,' Final Gradient Norm :',3x,d20.8)
         end if
      else if (digits .ge. 6) then
         if (grms .gt. 1.0d-6) then
            write (iout,70)  minimum,grms,gnorm
   70       format (/,' Final Function Value :',2x,f18.6,
     &              /,' Final RMS Gradient :',4x,f18.6,
     &              /,' Final Gradient Norm :',3x,f18.6)
         else
            write (iout,80)  minimum,grms,gnorm
   80       format (/,' Final Function Value :',2x,f18.6,
     &              /,' Final RMS Gradient :',4x,d18.6,
     &              /,' Final Gradient Norm :',3x,d18.6)
         end if
      else
         if (grms .gt. 1.0d-4) then
            write (iout,90)  minimum,grms,gnorm
   90       format (/,' Final Function Value :',2x,f16.4,
     &              /,' Final RMS Gradient :',4x,f16.4,
     &              /,' Final Gradient Norm :',3x,f16.4)
         else
            write (iout,100)  minimum,grms,gnorm
  100       format (/,' Final Function Value :',2x,f16.4,
     &              /,' Final RMS Gradient :',4x,d16.4,
     &              /,' Final Gradient Norm :',3x,d16.4)
         end if
      end if
      end if
c
c     write the final coordinates into a file
c
      if (rank.eq.0) then
        imin = freeunit ()
        open (unit=imin,file=minfile,status='old')
        rewind (unit=imin)
        call prtxyz (imin)
        close (unit=imin)
      end if
c
c     perform any final tasks before program exit
c
      call final
      call MPI_FINALIZE(ierr)
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  function minimiz1  --  energy and gradient for minimize  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "minimiz1" is a service routine that computes the energy and
c     gradient for a low storage BFGS optimization in Cartesian
c     coordinate space
c
c
      function minimiz1 (xx,g)
      use sizes
      use atoms
      use domdec
      use deriv ,only:info_forces,cDef
      use energi,only:info_energy
      use inform
      use scales
      use usage
      use utils,only:set_to_zero1m
      use utilgpu,only:rec_queue
      use mpi
      implicit none
      integer i,iglob,ierr
      real(r_p) minimiz1
      real(r_p) energy,e
      real(r_p) eps
      real(r_p) xx(*)
      real(r_p) g(*)
      real(r_p), allocatable :: derivs(:,:)
      logical analytic
      external energy
c
c     use either analytical or numerical gradients
c
      analytic = .true.
      eps = 0.00001_re_p
      if (rank.eq.0.and.tinkerdebug) write(*,'(x,a)') 'minimiz1'
c
c     translate optimization parameters to atomic coordinates
c
!$acc parallel loop async default(present)
      do i = 1, nloc
         iglob = glob(i)
         if (use(iglob)) then
            x(iglob) = xx(3*(iglob-1)+1) / scale(3*(iglob-1)+1)
            y(iglob) = xx(3*(iglob-1)+2) / scale(3*(iglob-1)+2)
            z(iglob) = xx(3*(iglob-1)+3) / scale(3*(iglob-1)+3)
         end if
      end do
c
c      call commstep
      call sendallpos
      call AllDirAssign
      call reassignpme(.false.)
c
c     perform dynamic allocation of some local arrays
c
      allocate (derivs(3,nbloc))
!$acc enter data create(derivs,e) async
      call set_to_zero1m(derivs,size(derivs),rec_queue)
      call reinitnl(0)
      call mechanicstep(0)
      call allocstep
      call nblist(0)
c
c     compute and store the energy and gradient
c
c      if (analytic) then
         call gradient (e,derivs)
         call commforces(derivs)
c      else
c         e = energy ()
c         call numgrad (energy,derivs,eps)
c      end if
      call allreduceen(e)

      !Debug
      if(deb_Force)  call info_forces(cDef)
      if(deb_Energy) call info_energy(rank)
!$acc update host(e) async
c
c     store Cartesian gradient as optimization gradient
c
!$acc parallel loop async default(present)
      do i = 1, nloc
         iglob = glob(i)
         if (use(iglob)) then
            g(3*(iglob-1)+1) = derivs(1,i) / scale(3*(iglob-1)+1)
            g(3*(iglob-1)+2) = derivs(2,i) / scale(3*(iglob-1)+2)
            g(3*(iglob-1)+3) = derivs(3,i) / scale(3*(iglob-1)+3)
         end if
      end do

!$acc wait
      minimiz1 = e
c
c     perform deallocation of some local arrays
c
!$acc exit data delete(derivs,e) async
      deallocate (derivs)
      end
