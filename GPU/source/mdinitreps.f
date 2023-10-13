c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ###############################################################
c     ##                                                                        ##
c     ##  subroutine mdinitreps  --  initialize a multiple replicas dyn  ##
c     ##                                                                        ##
c     ###############################################################
c
c
c     "mdinitreps" tries to restart a multiple replicas dynamics from
c      restart files
c
#include "tinker_precision.h"
c
      subroutine mdinitreps
      use domdec
      use files
      use inform
      use iounit
      use potent, only: use_lambdadyn
      use replicas
      implicit none
      integer lext,freeunit,idyn
      integer next
      integer num
      logical exist
      real(r_p), allocatable :: derivs(:,:)
      character*7 ext
      character*240 dynfile
      character*3 numberreps
      character*20 keyword
      character*240 record
      character*240 string
c
      write(numberreps, '(i3.3)') rank_reploc
c
c     try to restart using prior velocities and accelerations
c
      dynfile = filename(1:leng)//'_reps'//numberreps//'.dyn'
      call version (dynfile,'old')
      inquire (file=dynfile,exist=exist)
      if (exist) then
         idyn = freeunit ()
         open (unit=idyn,file=dynfile,status='old')
         rewind (unit=idyn)
         call readdyn (idyn)
         close (unit=idyn)
c
c     Do the domain decomposition
c
         call ddpme3d
         call reinitnl(0)
         call reassignpme(.true.)
         call mechanicstep(0)
         call allocstep
         call nblist(0)
#ifdef COLVARS
c
c       for lambda-dynamics, do a "blank" colvars computation to get restart value of lambda
c
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
c
c     
c
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
c
c     check if we want to start from different lambdas with lambdadyn
c
c     Colvars Feature Initialization
c
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
   20         continue
              lambdastart = list
              lambda = lambdastart(rank_reploc+1)
              write(*,*) 'lambda mdinitreps = ',lambda
           endif
         end do
       end if
      end
c
