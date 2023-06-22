c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine getxyz  --  get Cartesian coordinate structure  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "getxyz" asks for a Cartesian coordinate file name,
c     then reads in the coordinates file
c
c
#include "tinker_macro.h"
      subroutine get_xyz_filename(xyzfile)
      use files
      use iounit
      implicit none
      character*240, intent(inout) :: xyzfile

      if(.not. keys_already_read) then
        call init_keys
      endif 

      xyzfile = filename
      call suffix (xyzfile,'xyz','old')
      
      end subroutine

      subroutine getxyz
      use inform
      use iounit
      use output
      use files
      implicit none
      integer ixyz
      integer freeunit
      character*240 xyzfile

      call get_xyz_filename(xyzfile)
c
c     first open and then read the Cartesian coordinates file
c
      coordtype = 'CARTESIAN'
      ixyz = freeunit ()
      open (unit=ixyz,file=xyzfile,status='old')
      rewind (unit=ixyz)
      call readxyz (ixyz)
      close (unit=ixyz)
c
c     quit if the Cartesian coordinates file contains no atoms
c
      if (abort) then
         write (iout,30)
   30    format (/,' GETXYZ  --  Cartesian Coordinates File',
     &              ' does not Contain Any Atoms')
         call fatal
      end if
      end
