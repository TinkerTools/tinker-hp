c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ####################################################################
c     ##                                                                ##
c     ##  module files  --  name and number of current structure files  ##
c     ##                                                                ##
c     ####################################################################
c
c
c     nprior     number of previously existing cycle files
c     ldir       length in characters of the directory name
c     leng       length in characters of the base filename
c     filename   base filename used by default for all files
c     outfile    output filename used for intermediate results
c
c
      module files
      implicit none
      integer nprior,ldir,leng
      character*240, target:: filename,outfile
      logical :: keys_already_read=.FALSE.

      contains

      subroutine init_keys()
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
         write (6,10)
         read (5,20) xyzfile
         call basefile (xyzfile)
         call suffix (xyzfile,'xyz','old')
         inquire (file=xyzfile,exist=exist)
      end do

      keys_already_read = .true.

      end subroutine init_keys
      end
