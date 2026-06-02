!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ############################################################
!     ##  COPYRIGHT (C) 2018 by Joshua Rackers & Jay W. Ponder  ##
!     ##                   All Rights Reserved                  ##
!     ############################################################
!
!     ##################################################################
!     ##                                                              ##
!     ##  module kcpen  --  charge penetration forcefield parameters  ##
!     ##                                                              ##
!     ##################################################################
!
!
!
module kcpen
   use sizes
   implicit none
   real*8 :: cpele(maxclass) !<valence electron magnitude for each atom class
   real*8 :: cpalp(maxclass) !<alpha charge penetration parameter for each atom class
   save
end
