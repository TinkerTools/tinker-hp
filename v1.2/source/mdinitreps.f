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
c
      subroutine mdinitreps
      use files
      use inform
      use iounit
      use keys
      use mutant
      use replicas
      implicit none
      integer lext,freeunit,idyn
      integer next
      integer i,j,k,num
      logical exist
      character*7 ext
      character*240 dynfile
      character*3 numberreps
      character*20 keyword
      character*240 record
      character*240 string
      real*8, allocatable :: list(:)
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
      endif
c
c     check if we want to start from different lambdas with lambdadyn
c
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
c
