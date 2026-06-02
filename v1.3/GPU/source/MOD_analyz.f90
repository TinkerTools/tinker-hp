!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ###################################################################
!     ##                                                               ##
!     ##  module analyz  --  energy components partitioned over atoms  ##
!     ##                                                               ##
!     ###################################################################
!
!
!     aesum   total potential energy partitioned over atoms
!     aeb     bond stretch energy partitioned over atoms
!     aea     angle bend energy partitioned over atoms
!     aeba    stretch-bend energy partitioned over atoms
!     aeub    Urey-Bradley energy partitioned over atoms
!     aeaa    angle-angle energy partitioned over atoms
!     aeopb   out-of-plane bend energy partitioned over atoms
!     aeopd   out-of-plane distance energy partitioned over atoms
!     aeid    improper dihedral energy partitioned over atoms
!     aeit    improper torsion energy partitioned over atoms
!     aet     torsional energy partitioned over atoms
!     aept    pi-orbital torsion energy partitioned over atoms
!     aeat    angle-torsion energy partitioned over atoms
!     aebt    stretch-torsion energy partitioned over atoms
!     aett    torsion-torsion energy partitioned over atoms
!     aev     van der Waals energy partitioned over atoms
!     aer     Pauli repulsion energy partitioned over atoms
!     aedsp   damped dispersion energy partitioned over atoms
!     aec     charge-charge energy partitioned over atoms
!     aem     multipole energy partitioned over atoms
!     aep     polarization energy partitioned over atoms
!     aect    charge transfer energy partitioned over atoms
!     aeg     geometric restraint energy partitioned over atoms
!     aex     extra energy term partitioned over atoms
!
!
#include "tinker_macro.h"
module analyz
   implicit none
   real(t_p), allocatable :: aesum(:),aem (:),aep (:),aev (:)&
      &, aea(:), aeba(:), aub(:), aeaa(:), aeopd(:), aeid(:)&
      &, aeit(:), aet(:), aebt(:), aett(:), aeg(:), aex(:)&
      &, aeb(:), aeopb(:), aept(:), aec(:), aer(:), aedsp(:)&
      &, aect(:), aeat(:)
end
