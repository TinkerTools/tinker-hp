c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ###################################################################
c     ##                                                               ##
c     ##  module analyz  --  energy components partitioned over atoms  ##
c     ##                                                               ##
c     ###################################################################
c
c
c     aesum   total potential energy partitioned over atoms
c     aeb     bond stretch energy partitioned over atoms
c     aea     angle bend energy partitioned over atoms
c     aeba    stretch-bend energy partitioned over atoms
c     aeub    Urey-Bradley energy partitioned over atoms
c     aeaa    angle-angle energy partitioned over atoms
c     aeopb   out-of-plane bend energy partitioned over atoms
c     aeopd   out-of-plane distance energy partitioned over atoms
c     aeid    improper dihedral energy partitioned over atoms
c     aeit    improper torsion energy partitioned over atoms
c     aet     torsional energy partitioned over atoms
c     aept    pi-orbital torsion energy partitioned over atoms
c     aeat    angle-torsion energy partitioned over atoms
c     aebt    stretch-torsion energy partitioned over atoms
c     aett    torsion-torsion energy partitioned over atoms
c     aev     van der Waals energy partitioned over atoms
c     aer     Pauli repulsion energy partitioned over atoms
c     aedsp   damped dispersion energy partitioned over atoms
c     aec     charge-charge energy partitioned over atoms
c     aem     multipole energy partitioned over atoms
c     aep     polarization energy partitioned over atoms
c     aect    charge transfer energy partitioned over atoms
c     aeg     geometric restraint energy partitioned over atoms
c     aex     extra energy term partitioned over atoms
c
c
#include "tinker_precision.h"
      module analyz
      implicit none
      real(t_p), allocatable :: aesum(:),aem (:),aep (:),aev (:)
     &         , aea(:), aeba(:), aub(:), aeaa(:), aeopd(:), aeid(:)
     &         , aeit(:), aet(:), aebt(:), aett(:), aeg(:), aex(:)
     &         , aeb(:), aeopb(:), aept(:), aec(:), aer(:), aedsp(:)
     &         , aect(:), aeat(:)
      end
