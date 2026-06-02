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
!     idihunit  convert improper dihedral energy to kcal/mole
!     itorunit  convert improper torsion amplitudes to kcal/mole
!     torsunit  convert torsional parameter amplitudes to kcal/mole
!     ptorunit  convert pi-orbital torsion energy to kcal/mole
!     storunit  convert stretch-torsion energy to kcal/mole
!     atorunit  convert angle-torsion energy to kcal/mole
!     ttorunit  convert stretch-torsion energy to kcal/mole
!
!
#include "tinker_macro.h"
module torpot
   implicit none
   real(t_p) idihunit,itorunit,torsunit
   real(t_p) ptorunit,storunit,ttorunit
   real(t_p) atorunit
end
