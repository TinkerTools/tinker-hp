!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     #############################################################
!     ##                                                         ##
!     ##  subroutine getprm  --  get force field parameter file  ##
!     ##                                                         ##
!     #############################################################
!
!
!     "getprm" finds the potential energy parameter file
!     and then opens and reads the parameters
!
!
!> @brief 
!> finds the potential energy parameter file
!> and then opens and reads the parameters
!> @param no params
subroutine getprm
   use domdec
   use files
   use keys
   use inform
   use iounit
   use params
   use mpi
   implicit none
   integer i,iprm,next,length,ierr
   integer freeunit
   logical exist,useprm
   character*4 none
   character*20 keyword
   character*240 prmfile
   character*240 record
   character*240 string
   character*240 paramdir
!
   if (deb_Path) write(iout,*), 'getprm '
!
!
!     set the default name for the parameter file
!
   useprm = .true.
   prmfile = filename(1:leng)//'.prm'
!
!     search the keyword list for the parameter filename
!
   do i = 1, nkey
      next = 1
      record = keyline(i)
      call gettext (record,keyword,next)
      call upcase (keyword)
      if (keyword(1:11) .eq. 'PARAMETERS ') then
         string = record(next:240)
         next = 1
         call getstring (string,prmfile,next)
         if (next .eq. 1)  call gettext (string,prmfile,next)
      end if
   end do
!     Try to read $TINKER_PARAMS from environment
!     get_environment_variable is fortran intrinsic
   call get_environment_variable ('TINKER_PARAMS', paramdir,length)
!     if lenght is 0, the TINKER_PARAMS variable is not set
   if (length.ne.0) prmfile = trim(paramdir)//'/'//prmfile
!
!     check existence of default or specified parameter file
!
   call suffix (prmfile,'prm','old')
   inquire (file=prmfile,exist=exist)
!!
!!     test for user specified absence of a parameter file
!!
!   if (.not. exist) then
!      none = prmfile(1:4)
!      call upcase (none)
!      if (none .eq. 'NONE') then
!         exist = .true.
!         useprm = .false.
!      end if
!   end if
!
!     try to get a parameter filename from the command line
!
   if (.not. exist) then
      call nextarg (prmfile,exist)
      if (exist) then
         call suffix (prmfile,'prm','old')
         inquire (file=prmfile,exist=exist)
      end if
   end if
!
!     raise an error if no prm file given
!
   if (.not. exist) then
      if (rank.eq.0) write (iout,*) 'You need to Specify the parameter file in the keyfile or in the command line !'
      call MPI_BARRIER(COMM_TINKER,ierr)
      call fatal
   end if
   
!!     if necessary, ask for the parameter filename
!!
!   do while (.not. exist)
!      write (iout,10)
!10    format (/,' Enter Potential Parameter File Name :  ',$)
!      read (input,20)  prmfile
!20    format (a240)
!      next = 1
!      call getword (prmfile,none,next)
!      call upcase (none)
!      if (none.eq.'NONE' .and. next.eq.5) then
!         exist = .true.
!         useprm = .false.
!      else
!         call suffix (prmfile,'prm','old')
!         inquire (file=prmfile,exist=exist)
!      end if
!   end do
!
!     initialize force field control and parameter values
!
   call initprm
!
!     read the parameter file and store it for latter use
!
   nprm = 0
   if (useprm) then
      iprm = freeunit ()
      open (unit=iprm,file=prmfile,status='old')
      rewind (unit=iprm)
      do while (.true.)
         read (iprm,30,err=50,end=50)  record
30       format (a240)
         nprm = nprm + 1
         prmline(nprm) = record
         if (nprm .ge. maxprm) then
            write (iout,40)
40          format (/,' GETPRM  --  Parameter File Too Large;',&
            &' Increase MAXPRM')
            call fatal
         end if
      end do
50    continue
      close (unit=iprm)
   end if
!
!     get control and parameter values from the parameter file
!
   if (useprm)  call readprm
   return
end
