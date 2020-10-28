c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #################################################################
c     ##                                                             ##
c     ##  module couple  --  near-neighbor atom connectivity lists   ##
c     ##                                                             ##
c     #################################################################
c
c
c     maxn13   maximum number of atoms 1-3 connected to an atom
c     maxn14   maximum number of atoms 1-4 connected to an atom
c     maxn15   maximum number of atoms 1-5 connected to an atom
c
c     n12      number of atoms directly bonded to each atom
c     winn12    window object corresponding to n12
c     i12      atom numbers of atoms 1-2 connected to each atom
c     wini12    window object corresponding to i12
c     n13      number of atoms in a 1-3 relation to each atom
c     winn13    window object corresponding to n13
c     i13      atom numbers of atoms 1-3 connected to each atom
c     wini13    window object corresponding to i13
c     n14      number of atoms in a 1-4 relation to each atom
c     winn14    window object corresponding to n14
c     i14      atom numbers of atoms 1-4 connected to each atom
c     wini14    window object corresponding to i14
c     n15      number of atoms in a 1-5 relation to each atom
c     winn15    window object corresponding to n15
c     i15      atom numbers of atoms 1-5 connected to each atom
c     wini15    window object corresponding to i15
c
c     ninteract_saling_n   all interactions connected atoms (1-2,3,4,5)
c     allscal_n     list of all n scaling factor of the system
c     winallscal_n  window object corresponding to allscal_n
c     numscal_n     n scaling factor's number per atom
c     winnumscal_n  window object corresponding to allnumscal_n
c     scalbeg_n     index start point of n scaling factor in allscal_n
c     winscalbeg_n  window object corresponding to allscalbeg_n
c     typscal_n     pair to allscal_n to store n scaling factor type
c     wintypscal_n  window object corresponding to alltypscal_n
c     pair_factorn  temporary pair list of n scaling factor
c     scal_facotrn  temporary scaling factor value for each atoms
c     n_factorn     number of useful scaling interaction for each atoms
c     scan_factorn  scan of n_factorn
c     sum_factorn   total number of useful n scaling factor
c     max_facorn    maximum of n scaling factor for every atoms
c
c
#include "tinker_precision.h"
      module couple
      use sizes
      implicit none
      integer maxn13,maxn14,maxn15
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
