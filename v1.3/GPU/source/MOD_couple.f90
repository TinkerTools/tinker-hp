!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     #################################################################
!     ##                                                             ##
!     ##  module couple  --  near-neighbor atom connectivity lists   ##
!     ##                                                             ##
!     #################################################################
!
!
!     maxn12_  maximum number of atoms 1-2 connected to an atom at runTime
!     maxd12_  maximum index difference between atoms 1-2 connected to an atom at runTime
!     maxn13   maximum number of atoms 1-3 connected to an atom
!     maxn14   maximum number of atoms 1-4 connected to an atom
!     maxn15   maximum number of atoms 1-5 connected to an atom
!
!     n12      number of atoms directly bonded to each atom
!     winn12    window object corresponding to n12
!     i12      atom numbers of atoms 1-2 connected to each atom
!     wini12    window object corresponding to i12
!     n13      number of atoms in a 1-3 relation to each atom
!     winn13    window object corresponding to n13
!     i13      atom numbers of atoms 1-3 connected to each atom
!     wini13    window object corresponding to i13
!     n14      number of atoms in a 1-4 relation to each atom
!     winn14    window object corresponding to n14
!     i14      atom numbers of atoms 1-4 connected to each atom
!     wini14    window object corresponding to i14
!     n15      number of atoms in a 1-5 relation to each atom
!     winn15    window object corresponding to n15
!     i15      atom numbers of atoms 1-5 connected to each atom
!     wini15    window object corresponding to i15
!
!     ninteract_saling_n   all interactions connected atoms (1-2,3,4,5)
!     allscal_n     list of all n scaling factor of the system
!     winallscal_n  window object corresponding to allscal_n
!     numscal_n     n scaling factor's number per atom
!     winnumscal_n  window object corresponding to allnumscal_n
!     scalbeg_n     index start point of n scaling factor in allscal_n
!     winscalbeg_n  window object corresponding to allscalbeg_n
!     typscal_n     pair to allscal_n to store n scaling factor type
!     wintypscal_n  window object corresponding to alltypscal_n
!     pair_factorn  temporary pair list of n scaling factor
!     scal_facotrn  temporary scaling factor value for each atoms
!     n_factorn     number of useful scaling interaction for each atoms
!     scan_factorn  scan of n_factorn
!     sum_factorn   total number of useful n scaling factor
!     max_facorn    maximum of n scaling factor for every atoms
!
!
#include "tinker_macro.h"
module couple
   use sizes
   implicit none
   integer maxn12_,maxn13,maxn14,maxn15
   integer maxd12
   integer ninteract_scaling_n
   parameter (maxn13=4*maxvalue)
   parameter (maxn14=4*maxvalue)
   parameter (maxn15=4*maxvalue)
   integer, allocatable :: n12(:),i12(:,:)
   integer, pointer :: n13(:),i13(:,:)
   integer, pointer :: n14(:),i14(:,:),n15(:),i15(:,:)
   integer :: winn12,wini12,winn13,wini13
   integer :: winn14,wini14,winn15,wini15

   integer   ,pointer :: allscal_n(:),numscal_n(:),scalbeg_n(:)
   integer(1),pointer :: typscal_n(:)
   integer winallscal_n,winnumscal_n,winscalbeg_n,wintypscal_n
   integer   ,allocatable::pair_factorn(:,:)
   real(t_p) ,allocatable::scal_factorn(:,:)
   integer   ,allocatable::   n_factorn(:)
   integer   ,allocatable::scan_factorn(:)
   integer sum_factorn
   integer max_factorn

   ! TODO Fix Analyze: Remove this
!$acc declare create(n12,i12)
!$acc declare create(n13,i13)
!$acc declare create(n14,i14)
!$acc declare create(n15,i15)

end
