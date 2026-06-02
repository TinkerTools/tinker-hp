!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ####################################################################
!     ##                                                                ##
!     ##  module urypot  --  specifics of Urey-Bradley functional form  ##
!     ##                                                                ##
!     ####################################################################
!
module urypot
   implicit none
   real*8 :: cury !<cubic coefficient in Urey-Bradley potential
   real*8 :: qury !<quartic coefficient in Urey-Bradley potential
   real*8 :: ureyunit !<convert Urey-Bradley energy to kcal/mole
   save
end
