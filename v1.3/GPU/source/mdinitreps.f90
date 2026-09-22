!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ###############################################################
!     ##                                                                        ##
!     ##  subroutine mdinitreps  --  initialize a multiple replicas dyn  ##
!     ##                                                                        ##
!     ###############################################################
!
!
!     "mdinitreps" tries to restart a multiple replicas dynamics from
!      restart files
!
#include "tinker_precision.h"
!
subroutine mdinitreps
   use atomsMirror, only: n, x, y, z, reCast_position, download_position
   use boxes, only: xbox, ybox, zbox, alpha, beta, gamma
   use domdec
   use files
   use inform
   use iounit
   use moldyn, only: v, a, aalt, aalt2
   use potent, only: use_lambdadyn
   use replicas
   use restart, only: read_restart
   use status, only: failed, status_t
   implicit none
   integer lext
   integer next
   integer num
   logical exist,aalt2_present
   real(r_p), allocatable :: derivs(:,:)
   type(status_t) :: restart_status
   character*7 ext
   character*240 dynfile
   character*3 numberreps
   character*20 keyword
   character*240 record
   character*240 string
!
   write(numberreps, '(i3.3)') rank_reploc
!
!     try to restart using prior velocities and accelerations
!
   dynfile = filename(1:leng)//'_reps'//numberreps//'.dyn'
   call version (dynfile,'old')
   inquire (file=dynfile,exist=exist)
   if (exist) then
      call read_restart(trim(dynfile),n,xbox,ybox,zbox,alpha,beta,gamma,&
         &x,y,z,v,a,aalt,aalt2,aalt2_present,restart_status)
      if (failed(restart_status)) then
         if (rank == 0) write (iout,'(/,a)') trim(restart_status%message)
         call fatal
      end if
      call lattice
!$acc update device(x(:),y(:),z(:),v,a,aalt,aalt2)
      if (t_p.ne.r_p) then
         call reCast_position
!$acc wait
         call download_position
      end if
!
!     Do the domain decomposition
!
      call ddpme3d
      call reinitnl(0)
      call reassignpme(.true.)
      call mechanicstep(0)
      call allocstep
      call nblist(0)
#ifdef COLVARS
!
!       for lambda-dynamics, do a "blank" colvars computation to get restart value of lambda
!
      if (use_lambdadyn) then
         allocate(derivs(3,nbloc))
         derivs = 0d0
         call allocstep
!$acc enter data create(derivs)
         call colvars_run(derivs)
!$acc    exit data delete(derivs) async
      end if
#endif
   endif
end subroutine
!
!
!
subroutine lambda_init_reps
   use keys
   use mutant
   use potent, only: use_lambdadyn
   use replicas
   implicit none
   integer i,j,num,next
   character*20 keyword
   character*240 record
   character*240 string
   real(t_p), allocatable :: list(:)
!
!     check if we want to start from different lambdas with lambdadyn
!
!     Colvars Feature Initialization
!
   if (use_lambdadyn) then
      allocate (lambdastart(nreps))
      allocate (list(nreps))
      do j = 1, nkey
         next = 1
         record = keyline(j)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:12) .eq. 'LAMBDASTART ') then
            num = 0
            list = 0d0
            string = record(next:240)
            read (string,*,err=20,end=20)  (list(i),i=1,nreps)
20          continue
            lambdastart = list
            lambda = lambdastart(rank_reploc+1)
            write(*,*) 'lambda mdinitreps = ',lambda
         endif
      end do
   end if
end
!
