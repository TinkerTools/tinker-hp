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
module analyz
   implicit none
   real*8, allocatable :: aesum(:) !<total potential energy partitioned over atoms
   real*8, allocatable :: aem(:)   !<multipole energy partitioned over atoms
   real*8, allocatable :: aep(:)   !<polarization energy partitioned over atoms
   real*8, allocatable :: aev(:)   !<van der Waals energy partitioned over atoms
   real*8, allocatable :: aeb(:)   !<bond stretch energy partitioned over atoms
   real*8, allocatable :: aea(:)   !<angle bend energy partitioned over atoms
   real*8, allocatable :: aeba(:)  !<stretch-bend energy partitioned over atoms
   real*8, allocatable :: aub(:)   !<Urey-Bradley energy partitioned over atoms
   real*8, allocatable :: aeaa(:)  !<angle-angle energy partitioned over atoms
   real*8, allocatable :: aeopb(:) !<out-of-plane bend energy partitioned over atoms
   real*8, allocatable :: aeopd(:) !<out-of-plane distance energy partitioned over atoms
   real*8, allocatable :: aeid(:)  !<improper dihedral energy partitioned over atoms
   real*8, allocatable :: aeit(:)  !<improper torsion energy partitioned over atoms
   real*8, allocatable :: aet(:)   !<torsional energy partitioned over atoms
   real*8, allocatable :: aept(:)  !<pi-orbital torsion energy partitioned over atoms
   real*8, allocatable :: aebt(:)  !<stretch-torsion energy partitioned over atoms
   real*8, allocatable :: aett(:)  !<torsion-torsion energy partitioned over atoms
   real*8, allocatable :: aeg(:)   !<geometric restraint energy partitioned over atoms
   real*8, allocatable :: aex(:)   !<extra energy term partitioned over atoms
   real*8, allocatable :: aec(:)   !<charge-charge energy partitioned over atoms
   real*8, allocatable :: aer(:)   !<Pauli repulsion energy partitioned over atoms
   real*8, allocatable :: aedsp(:) !<damped dispersion energy partitioned over atoms
   real*8, allocatable :: aect(:)  !<charge transfer energy partitioned over atoms
   real*8, allocatable :: aeat(:)  !<angle-torsion energy partitioned over atoms
   save
end
