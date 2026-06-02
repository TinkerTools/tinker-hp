!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ##################################################################
!     ##                                                              ##
!     ##  program minimize  --  low storage BFGS Cartesian optimizer  ##
!     ##                                                              ##
!     ##################################################################
!
!
!     "minimize" performs energy minimization in Cartesian coordinate
!     space using a low storage BFGS nonlinear optimization
!
!
#include "tinker_macro.h"
program minimize
   use mpi
   implicit none
   integer ierr
   call MPI_INIT(ierr)
   call minimize_bis
   call MPI_BARRIER(MPI_COMM_WORLD,ierr)
   call MPI_FINALIZE(ierr)
end

subroutine minimize_bis
   use atomsMirror
   use domdec
   use deriv  ,only: ftot_l,cDef,dr_stride3,de_tot&
      &,info_forces,comm_forces
   use energi ,only: info_energy
   use files
   use keys
   use inform
   use iounit
   use interfaces,only: LBFGS
   use scales
   use tinMemory ,only: prmem_requestm
   use usage
   use utilgpu   ,only: rec_queue,ti_p,re_p
   use mpi
   implicit none
   integer i,j,imin,ierr
   integer next,freeunit
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
   character*240 minfile
   character*240 record
   character*240 string
   external energy
   external minimiz1
   external optsave
!
   ! Sign running program
   app_id = minimize_a
!
!     set up the structure and mechanics calculation
!
   call initial
   call initmpi
   call getxyz
   call cutoffs
   call unitcell
   call lattice
!
!     setup for MPI
!
   call drivermpi
   call reinitnl(0)
!
   call mechanic
   call nblist(0)
   call allocstep
!
   if (associated(scale)) then
!$acc exit data delete(scale)
      deallocate (scale)
   end if
   allocate (scale(3*n))
!$acc enter data create(scale)
!
!     use either analytical or numerical gradients
!
   analytic = .true.
   eps = 0.00001_re_p
!
!     search the keywords for output frequency parameters
!
   do i = 1, nkey
      next = 1
      record = keyline(i)
      call gettext (record,keyword,next)
      call upcase (keyword)
      string = record(next:240)
      if (keyword(1:9) .eq. 'PRINTOUT ') then
         read (string,*,err=10,end=10)  iprint
      else if (keyword(1:9) .eq. 'WRITEOUT ') then
         read (string,*,err=10,end=10)  iwrite
      end if
10    continue
   end do
!
!     get termination criterion as RMS gradient per atom
!
   grdmin = -1.0_re_p
   call nextarg (string,exist)
   if (.not.exist) then
     if (ranktot.eq.0) write (iout,*) 'You Need To Enter RMS Gradient per Atom Criterion'
     call MPI_BARRIER(COMM_TINKER,ierr)
     __TINKER_FATAL__
   else
     read (string,*,err=20,end=20)  grdmin
20 continue
   end if
   if (grdmin .le. 0.0_re_p)  grdmin = 0.01_re_p
!
!     write out a copy of coordinates for later update
!
   imin = freeunit ()
   minfile = filename(1:leng)//'.xyz'
   call version (minfile,'new')
   if (rank.eq.0) then
      open (unit=imin,file=minfile,status='new')
      call prtxyz (imin)
      close (unit=imin)
      outfile = minfile
   end if
!
!     set scaling parameter for function and derivative values;
!     use square root of median eigenvalue of typical Hessian
!
   set_scale = .true.
!$acc parallel loop async present(scale)
   do i = 1, 3*n
      scale(i) = 12.0_re_p
   end do
!
!     perform dynamic allocation of some local arrays
!
   allocate (xx(3*n))
!$acc enter data create(xx,minimum) async
!
!     scale the coordinates of each active atom
!
!$acc parallel loop async default(present)
   do i = 1, nloc
      iglob = glob(i)
      if (use(iglob)) then
         xx(3*(iglob-1)+1) = x(iglob) * scale(3*(iglob-1)+1)
         xx(3*(iglob-1)+2) = y(iglob) * scale(3*(iglob-1)+2)
         xx(3*(iglob-1)+3) = z(iglob) * scale(3*(iglob-1)+3)
      end if
   end do
!
!     make the call to the optimization routine
!
   call lbfgs (n,xx,minimum,grdmin,minimiz1,optsave)
!
!     unscale the final coordinates for active atoms
!
!$acc parallel loop async default(present)
   do i = 1, nloc
      iglob = glob(i)
      if (use(iglob)) then
         x(iglob) = xx(3*(iglob-1)+1) / scale(3*(iglob-1)+1)
         y(iglob) = xx(3*(iglob-1)+2) / scale(3*(iglob-1)+2)
         z(iglob) = xx(3*(iglob-1)+3) / scale(3*(iglob-1)+3)
      end if
   end do
!
!     compute the final function and RMS gradient values
!
!      if (analytic) then
!         call commstep
   call sendallpos
   call reCast_position
   call AllDirAssign
   call reassignpme(.false.)
   call reinitnl(0)
   call mechanicstep(0)
   call allocstep
   call nblist(0)
   call gradient (minimum,derivs)
!$acc wait
!$acc update host(minimum)
   call MPI_ALLREDUCE(MPI_IN_PLACE,minimum,1,MPI_RPREC,&
      &MPI_SUM,MPI_COMM_WORLD,ierr)
   call comm_forces( derivs )
!      else
!         minimum = energy ()
!c         call numgrad (energy,derivs,eps)
!      end if

   !Debug
   if(deb_Force) call info_forces(cDef)
   if(deb_Energy)call info_energy(rank)

   gnorm = 0.0_re_p
   call normp(de_tot(1:,1),dr_stride3,gnorm)
!!$acc parallel loop collapse(2) async present(de_tot)
!      do i = 1, nloc; do j = 1, 3
!         iglob = glob(i)
!         if (use(iglob)) then
!            gnorm = gnorm + mdr2md(de_tot(j,i))**2
!         end if
!      end do; end do

   call MPI_ALLREDUCE(MPI_IN_PLACE,gnorm,1,MPI_RPREC,&
      &MPI_SUM,COMM_TINKER,ierr)
   gnorm = sqrt(gnorm)
   grms  = gnorm / sqrt(real(3*n,r_p))
!
!     perform deallocation of some local arrays
!
!$acc exit data delete(xx,minimum) async
   deallocate (xx)
!
!     write out the final function and gradient values
!
   if (rank.eq.0) then
      if (digits .ge. 8) then
         if (grms .gt. 1.0d-8) then
            write (iout,50)  minimum,grms,gnorm
50          format (/,' Final Function Value :',2x,f20.8,&
               &/,' Final RMS Gradient :',4x,f20.8,&
               &/,' Final Gradient Norm :',3x,f20.8)
         else
            write (iout,60)  minimum,grms,gnorm
60          format (/,' Final Function Value :',2x,f20.8,&
               &/,' Final RMS Gradient :',4x,d20.8,&
               &/,' Final Gradient Norm :',3x,d20.8)
         end if
      else if (digits .ge. 6) then
         if (grms .gt. 1.0d-6) then
            write (iout,70)  minimum,grms,gnorm
70          format (/,' Final Function Value :',2x,f18.6,&
               &/,' Final RMS Gradient :',4x,f18.6,&
               &/,' Final Gradient Norm :',3x,f18.6)
         else
            write (iout,80)  minimum,grms,gnorm
80          format (/,' Final Function Value :',2x,f18.6,&
               &/,' Final RMS Gradient :',4x,d18.6,&
               &/,' Final Gradient Norm :',3x,d18.6)
         end if
      else
         if (grms .gt. 1.0d-4) then
            write (iout,90)  minimum,grms,gnorm
90          format (/,' Final Function Value :',2x,f16.4,&
               &/,' Final RMS Gradient :',4x,f16.4,&
               &/,' Final Gradient Norm :',3x,f16.4)
         else
            write (iout,100)  minimum,grms,gnorm
100         format (/,' Final Function Value :',2x,f16.4,&
               &/,' Final RMS Gradient :',4x,d16.4,&
               &/,' Final Gradient Norm :',3x,d16.4)
         end if
      end if
   end if
!
!     write the final coordinates into a file
!
   if (rank.eq.0) then
      imin = freeunit ()
      open (unit=imin,file=minfile,status='old')
      rewind (unit=imin)
      call prtxyz (imin)
      close (unit=imin)
   end if
!
!     perform any final tasks before program exit
!
   call final
end
!
!
!     ###############################################################
!     ##                                                           ##
!     ##  function minimiz1  --  energy and gradient for minimize  ##
!     ##                                                           ##
!     ###############################################################
!
!
!     "minimiz1" is a service routine that computes the energy and
!     gradient for a low storage BFGS optimization in Cartesian
!     coordinate space
!
!
function minimiz1 (xx,g)
   use sizes
   use atomsMirror
   use domdec
   use deriv  ,only:info_forces,cDef,ftot_l,comm_forces,get_ftot&
      &,de_tot,dr_stride3
   use energi ,only:info_energy
   use inform
   use scales
   use tinMemory,only:mipk
   use usage
   use utils  ,only:set_to_zero1m
   use utilgpu,only:rec_queue,rec_stream,mem_set,ti_p,re_p
   use mpi
   implicit none
   integer i,iglob,ierr
   real(r_p) minimiz1
   real(r_p) energy,e
   real(r_p) eps
   real(r_p) xx(*)
   real(r_p) g(*)
   real(r_p), allocatable :: derivs(:,:)
   mdyn_rtyp  zero_md
   integer(mipk) siz_
   logical analytic
   external energy
   parameter(zero_md=0)
!
!     use either analytical or numerical gradients
!
   analytic = .true.
   eps = 0.00001_re_p
   if (deb_Path) write(*,'(x,a)') 'minimiz1'
!
!     translate optimization parameters to atomic coordinates
!
!$acc parallel loop async default(present)
   do i = 1, nloc
      iglob = glob(i)
      if (use(iglob)) then
         x(iglob) = xx(3*(iglob-1)+1) / scale(3*(iglob-1)+1)
         y(iglob) = xx(3*(iglob-1)+2) / scale(3*(iglob-1)+2)
         z(iglob) = xx(3*(iglob-1)+3) / scale(3*(iglob-1)+3)
      end if
   end do
!
!      call commstep
   call sendallpos
   call reCast_position
   call AllDirAssign
   call reassignpme(.false.)
!
!     perform dynamic allocation of some local arrays
!
   allocate (derivs(3,nbloc))
!$acc enter data create(derivs,e)
   call set_to_zero1m(derivs,size(derivs),rec_queue)
   call reinitnl(0)
   call mechanicstep(0)
   call allocstep
   call nblist(0)
!
!     compute and store the energy and gradient
!
!      if (analytic) then
   call gradient (e,derivs)
   call comm_forces(derivs)
!      else
!         e = energy ()
!         call numgrad (energy,derivs,eps)
!      end if
   if (ftot_l) then
      siz_ = dr_stride3
      call get_ftot(derivs,nbloc)
      call mem_set (de_tot,zero_md,siz_,rec_stream)
   end if
   call allreduceen(e)

   !Debug
   if(deb_Force)  call info_forces(cDef)
   if(deb_Energy) call info_energy(rank)
!$acc update host(e) async
!
!     store Cartesian gradient as optimization gradient
!
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
!
!     perform deallocation of some local arrays
!
!$acc exit data delete(derivs,e)
   deallocate (derivs)
end
