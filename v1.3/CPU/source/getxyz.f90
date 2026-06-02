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
!> @brief 
!> asks for a Cartesian coordinate file name,
!> then reads in the coordinates file
!> @param no params
subroutine getxyz
   use inform
   use iounit
   use output
   implicit none
   integer ixyz
   integer freeunit
   character*240 xyzfile
!
   if (deb_Path) write(iout,*), 'getxyz '
!

   call get_xyz_filename(xyzfile)
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

subroutine get_xyz_filename(xyzfile)
   use iounit
   use files
   implicit none
   character*240, intent(inout) :: xyzfile

   if(.not. keys_already_read) then
      call init_keys
   endif

   xyzfile = filename
   call suffix (xyzfile,'xyz','old')

end subroutine
