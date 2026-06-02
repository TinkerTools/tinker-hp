!
!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ################################################################
!     ##                                                            ##
!     ##  module kdsp  --  damped dispersion forcefield parameters  ##
!     ##                                                            ##
!     ################################################################
!
!
module kdsp
   use sizes
   implicit none
   real*8 dspsix(maxtyp) !<C6 dispersion coefficient for each atom class
   real*8 dspdmp(maxtyp) !<alpha dispersion parameter for each atom class
   save
end
