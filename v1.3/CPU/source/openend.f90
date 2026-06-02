!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     #################################################################
!     ##                                                             ##
!     ##  subroutine openend  --  open a file positioned for append  ##
!     ##                                                             ##
!     #################################################################
!
!
!     "openend" opens a file on a Fortran unit such that the position
!     is set to the bottom for appending to the end of the file
!
!     note this routine is system dependent since the Fortran 90
!     standard is not supported by many Fortran 77 compilers; only
!     one of the various implementations below should be activated
!     by removing comment characters
!
!
!> @brief 
!> opens a file on a Fortran unit such that the position
!> is set to the bottom for appending to the end of the file
!> @param no params
subroutine openend (iunit,name)
   implicit none
   integer iunit
   character*240 name
!
!
!     standard Fortran 90, unavailable in some Fortran 77 compilers
!
   open (unit=iunit,file=name,status='old',position='append')
!
!     common extension supported by many Fortran 77 compilers
!
!     open (unit=iunit,file=name,status='old',access='append')
!
!     some Fortran 77 compilers open files for append by default
!
!     open (unit=iunit,file=name,status='old')
!
!     manually read to the end of file, slow but always correct
!
!     open (unit=iunit,file=name,status='old')
!     do while (.true.)
!        read (iunit,10,err=20,end=20)
!  10    format ()
!     end do
!  20 continue
   return
end
