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
!     nprm      number of nonblank lines in the parameter file
!     prmline   contents of each individual parameter file line
!
!
module params
   use sizes
   implicit none
   integer nprm
   character*240 prmline(maxprm)
   save
end
