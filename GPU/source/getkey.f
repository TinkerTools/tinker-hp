c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine getkey  --  find and store contents of keyfile  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "getkey" finds a valid keyfile and stores its contents as
c     line images for subsequent keyword parameter searching
c
c
#include "tinker_precision.h"
      subroutine getkey
      use argue
      use dcdmod
      use domdec
      use files
      use keys
      use iounit
      use potent
      use replicas
      implicit none
      integer i,ikey
      integer next,length
      integer freeunit
      integer trimtext
      integer omp_get_num_threads
      logical exist,header
      character*20 keyword
      character*240 keyfile
      character*240 comment
      character*240 record
      character*240 string
c
c
c     check for a keyfile specified on the command line
c
      exist = .false.
      do i = 1, narg-1
         string = arg(i)
         call upcase (string)
         if (string(1:2) .eq. '-K') then
            keyfile = arg(i+1)
            call suffix (keyfile,'key','old')
            inquire (file=keyfile,exist=exist)
            if (.not. exist) then
               write (iout,10)
   10          format (/,' GETKEY  --  Keyfile Specified',
     &                    ' on Command Line was not Found')
               call fatal
            end if
         end if
      end do
c
c     try to get keyfile from base name of current system
c
      if (.not. exist) then
         keyfile = filename(1:leng)//'.key'
         call version (keyfile,'old')
         inquire (file=keyfile,exist=exist)
      end if
c
c     check for the existence of a generic keyfile
c
      if (.not. exist) then
         if (ldir .eq. 0) then
            keyfile = 'tinker.key'
         else
            keyfile = filename(1:ldir)//'tinker.key'
         end if
         call version (keyfile,'old')
         inquire (file=keyfile,exist=exist)
      end if
c
c     read the keyfile and store it for latter use
c
      nkey = 0
      if (exist) then
         ikey = freeunit ()
         open (unit=ikey,file=keyfile,status='old')
         rewind (unit=ikey)
         do while (.true.)
            read (ikey,20,err=40,end=40)  record
   20       format (a240)
            nkey = nkey + 1
            keyline(nkey) = record
            if (nkey .ge. maxkey) then
               write (iout,30)
   30          format (/,' GETKEY  --  Keyfile Too Large;',
     &                    ' Increase MAXKEY')
               call fatal
            end if
         end do
   40    continue
         close (unit=ikey)
      end if
c
c     check for comment lines to be echoed to the output
c
      header = .true.
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:5) .eq. 'ECHO ') then
            comment = record(next:240)
            length = trimtext (comment)
            if (header) then
               header = .false.
               write (iout,50)
   50          format ()
            end if
            if (length .eq. 0) then
               write (iout,60)
   60          format ()
            else
               write (iout,70)  comment(1:length)
   70          format (a)
            end if
         end if
      end do
c
c     set number of threads for OpenMP parallelization
c
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call upcase (record)
         call gettext (record,keyword,next)
         string = record(next:240)
         if (keyword(1:15) .eq. 'OPENMP-THREADS ') then
            read (string,*,err=80,end=80)  nthread
            !call omp_set_num_threads (nthread)
         end if
   80    continue
      end do
c
c     separate group of MPI process for PME reciprocal space
c     set number of MPI process for PME reciprocal space
c
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call upcase (record)
         call gettext (record,keyword,next)
         string = record(next:240)
         if (keyword(1:15) .eq. 'PME-PROCS ') then
            nrec = 0
            use_pmecore = .true.
c!$acc update device(use_pmecore)
            read (string,*,err=90,end=90)  nrec
   90       continue
            if (nrec.eq.0) then
              if (rank.eq.0) write (iout,*) 'no procs for PME, quitting'
              call fatal
            end if
            if (rank.eq.0) then
 1000         format(/, I10,' Separate MPI processes for PME')
              write(iout,1000) nrec
            end if
         end if
      end do
c
c     write/read trajectory in the dcd format
c
      dcdio = .false.
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call upcase (record)
         call gettext (record,keyword,next)
         string = record(next:240)
         if (keyword(1:6) .eq. 'DCDIO ') then
            dcdio = .true.
         end if
      end do
c
c     read number of replicas in multiple replicas run
c
      use_reps = .false.
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call upcase (record)
         call gettext (record,keyword,next)
         string = record(next:240)
         if (keyword(1:9) .eq. 'REPLICAS ') then
            read (string,*,err=110,end=110)  nreps
            use_reps = .true.
  110       continue
         end if
      end do
c
      return
      end

      subroutine init_keys()
      use iounit
      use files
      implicit none
      integer ixyz
      integer freeunit
      logical exist
      character*240 xyzfile
c
c
c     try to get a filename from the command line arguments
c
      call nextarg (xyzfile,exist)
      if (exist) then
         call basefile (xyzfile)
         call suffix (xyzfile,'xyz','old')
         inquire (file=xyzfile,exist=exist)
      end if
c
c     ask for the user specified input structure filename
c
      do while (.not. exist)
   10    format(/,' Enter Cartesian Coordinate File Name :  ',$)
   20    format (a240)
         write (iout,10)
         read (input,20) xyzfile
         call basefile (xyzfile)
         call suffix (xyzfile,'xyz','old')
         inquire (file=xyzfile,exist=exist)
      end do

      keys_already_read = .true.

      end subroutine init_keys
