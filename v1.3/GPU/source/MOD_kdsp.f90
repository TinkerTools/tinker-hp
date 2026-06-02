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
!     dspsix   C6 dispersion coefficient for each atom class
!     dspdmp   alpha dispersion parameter for each atom class
!
!
module kdsp
   use sizes
   implicit none
   real*8 dspsix(maxtyp)
   real*8 dspdmp(maxtyp)
   save
end
