!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     #####################################################################
!     ##                                                                 ##
!     ##  module chgpot  --  specifics of charge-charge functional form  ##
!     ##                                                                 ##
!     #####################################################################
!
!
module chgpot
   implicit none
   real*8 :: electric !<energy factor in kcal/mole for current force field
   real*8 :: dielec !<dielectric constant for electrostatic interactions
   real*8 :: ebuffer !<electrostatic buffering constant added to distance
   real*8 :: c2scale !<factor by which 1-2 charge interactions are scaled
   real*8 :: c3scale !<factor by which 1-3 charge interactions are scaled
   real*8 :: c4scale !<factor by which 1-4 charge interactions are scaled
   real*8 :: c5scale !<factor by which 1-5 charge interactions are scaled
   logical :: neutnbr !<logical flag governing use of neutral group neighbors
   logical :: neutcut !<logical flag governing use of neutral group cutoffs
   save
end
