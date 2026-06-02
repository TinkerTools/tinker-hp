!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ##################################################################
!     ##                                                              ##
!     ##  module kpolr  --  forcefield parameters for polarizability  ##
!     ##                                                              ##
!     ##################################################################
!
!
module kpolr
   use sizes
   implicit none
   integer :: pgrp(maxvalue,maxtyp) !<connected types in polarization group of each atom type
   real*8 :: polr(maxtyp) !<dipole polarizability parameters for each atom type
   real*8 :: athl(maxtyp) !<Thole polarizability damping value for each atom type
   real*8 :: dthl(maxtyp) !<alternate Thole direct polarization damping values
   save
end
