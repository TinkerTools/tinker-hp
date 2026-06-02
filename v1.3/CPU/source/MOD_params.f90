!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     #################################################################
!     ##                                                             ##
!     ##  module params  --  contents of force field parameter file  ##
!     ##                                                             ##
!     #################################################################
!
!
module params
   use sizes
   implicit none
   integer :: nprm !<number of nonblank lines in the parameter file
   character*240 :: prmline(maxprm) !<contents of each individual parameter file line
   save
end
