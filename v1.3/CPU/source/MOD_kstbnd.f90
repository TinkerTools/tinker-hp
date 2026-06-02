!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     #################################################################
!     ##                                                             ##
!     ##  module kstbnd  --  forcefield parameters for stretch-bend  ##
!     ##                                                             ##
!     #################################################################
!
!
!
module kstbnd
   implicit none
   integer :: maxnsb !<maximum number of stretch-bend parameter entries
   parameter (maxnsb=2000)
   real*8 :: stbn(2,maxnsb) !<force constant parameters for stretch-bend terms
   character*12 :: ksb(maxnsb) !<string of atom classes for stretch-bend terms
   save
end
