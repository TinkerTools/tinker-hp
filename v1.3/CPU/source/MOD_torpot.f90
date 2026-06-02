!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ##################################################################
!     ##                                                              ##
!     ##  module torpot  --  specifics of torsional functional forms  ##
!     ##                                                              ##
!     ##################################################################
!
!
!
module torpot
   implicit none
   real*8 :: idihunit !<convert improper dihedral energy to kcal/mole
   real*8 :: itorunit !<convert improper torsion amplitudes to kcal/mole
   real*8 :: torsunit !<convert torsional parameter amplitudes to kcal/mole
   real*8 :: ptorunit !<convert pi-orbital torsion energy to kcal/mole
   real*8 :: storunit !<convert stretch-torsion energy to kcal/mole
   real*8 :: ttorunit !<convert angle-torsion energy to kcal/mole
   real*8 :: atorunit !<convert stretch-torsion energy to kcal/mole
   save
end
