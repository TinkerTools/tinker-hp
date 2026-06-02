!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     #####################################################################
!     ##                                                                 ##
!     ##  subroutine mdinitreps  --  initialize a multiple replicas dyn  ##
!     ##                                                                 ##
!     #####################################################################
!
!
!     "mdinitreps" tries to restart a multiple replicas dynamics from
!      restart files
!
!
!> @brief 
!> tries to restart a multiple replicas dynamics from
!>  restart files
!> @param no params
subroutine mdinitreps
   use files
   use inform
   use iounit
   use keys
   use mutant
   use replicas
   use domdec, only: rank
   implicit none
   integer freeunit,idyn
   integer next
   integer i,j,num
   logical exist
   character*240 dynfile
   character*3 numberreps
   character*20 keyword
   character*240 record
   character*240 string
   real*8, allocatable :: list(:)
!
   if (deb_Path) write(iout,*), 'mdinitreps '
!
!
   write(numberreps, '(i3.3)') rank_reploc
!
!     try to restart using prior velocities and accelerations
!
   dynfile = filename(1:leng)//'_reps'//numberreps//'.dyn'
   call version (dynfile,'old')
   inquire (file=dynfile,exist=exist)
   if (exist) then
      if(rank==0) then
         write(*,*) " --- Tinker-HP loading Restart file "&
         &//trim(dynfile)//"  ---"
      endif
      idyn = freeunit ()
      open (unit=idyn,file=dynfile,status='old')
      rewind (unit=idyn)
      call readdyn (idyn)
      close (unit=idyn)
!
!     Do the domain decomposition
!
      call ddpme3d
      call reinitnl(0)
      call reassignpme(.true.)
      call mechanic_up_para(0)
      call allocstep
      call nblist(0)
   endif
!
!     check if we want to start from different lambdas with lambdadyn
!
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
20       continue
         lambdastart = list
         lambda = lambdastart(rank_reploc+1)
      endif
   end do
end
!
