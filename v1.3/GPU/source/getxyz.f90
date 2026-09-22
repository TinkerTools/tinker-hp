!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     #################################################################
!     ##                                                             ##
!     ##  subroutine getxyz  --  get Cartesian coordinate structure  ##
!     ##                                                             ##
!     #################################################################
!
!
!     "getxyz" asks for a Cartesian coordinate file name,
!     then reads in the coordinates file
!
!
#include "tinker_macro.h"
subroutine get_xyz_filename(xyzfile, key_is_loaded)
   use files
   use iounit
   implicit none
   character*240, intent(inout) :: xyzfile
   logical, intent(in), optional :: key_is_loaded

   if (present(key_is_loaded)) then
      if (.not.key_is_loaded .and. .not.keys_already_read) call init_keys
   else
      if (.not.keys_already_read) call init_keys
   end if

   xyzfile = filename
   call suffix (xyzfile,'xyz','old')

end subroutine

subroutine getxyz(key_is_loaded)
   use inform
   use interfaces, only: get_xyz_filename
   use iounit
   use output
   use files
   implicit none
   integer ixyz
   integer freeunit
   character*240 xyzfile
   logical, intent(in), optional :: key_is_loaded

   if (present(key_is_loaded)) then
      call get_xyz_filename(xyzfile,key_is_loaded)
   else
      call get_xyz_filename(xyzfile)
   end if
!
!     first open and then read the Cartesian coordinates file
!
   coordtype = 'CARTESIAN'
   ixyz = freeunit ()
   open (unit=ixyz,file=xyzfile,status='old')
   rewind (unit=ixyz)
   call readxyz (ixyz)
   close (unit=ixyz)
!
!     quit if the Cartesian coordinates file contains no atoms
!
   if (abort) then
      write (iout,30)
30    format (/,' GETXYZ  --  Cartesian Coordinates File',&
         &' does not Contain Any Atoms')
      call fatal
   end if
end
