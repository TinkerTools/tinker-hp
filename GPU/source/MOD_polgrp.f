c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #####################################################################
c     ##                                                                 ##
c     ##  module polgrp  --  polarizable site group connectivity lists   ##
c     ##                                                                 ##
c     #####################################################################
c
c
c     maxp11   maximum number of atoms in a polarization group
c     maxp12   maximum number of atoms in groups 1-2 to an atom
c     maxp13   maximum number of atoms in groups 1-3 to an atom
c     maxp14   maximum number of atoms in groups 1-4 to an atom
c
c     np11     number of atoms in polarization group of each atom
c     winnp11    window object corresponding to np11
c     ip11     atom numbers of atoms in same group as each atom
c     winip11    window object corresponding to ip11
c     np12     number of atoms in groups 1-2 to each atom
c     winnp12    window object corresponding to np12
c     ip12     atom numbers of atoms in groups 1-2 to each atom
c     winip12    window object corresponding to ip12
c     np13     number of atoms in groups 1-3 to each atom
c     winnp13    window object corresponding to np13
c     ip13     atom numbers of atoms in groups 1-3 to each atom
c     winip13    window object corresponding to ip13
c     np14     number of atoms in groups 1-4 to each atom
c     winnp14    window object corresponding to np14
c     ip14     atom numbers of atoms in groups 1-4 to each atom
c     winip14    window object corresponding to ip14
c
c     allscal_p     list of all p scaling factor of the system
c     winallscal_p  window object corresponding to allscal_p
c     numscal_p     p scaling factor's number per atom
c     winnumscal_p  window object corresponding to allnumscal_p
c     scalbeg_p     index start point of p scaling factor in allscal_n
c     winscalbeg_p  window object corresponding to allscalbeg_p
c     typscal_p     pair to allscal_p to store p scaling factor type
c     wintypscal_p  window object corresponding to alltypscal_p
c     pair_factorp  pair list of p scaling factor
c     scal_facotrp  scaling factor value for each atoms
c     n_factorp     number of scaling useful scaling interaction for each atoms
c     n_factordp    number of scaling useful pair (d,p) scaling interaction for each atoms
c     scan_factorp  partial sum of n_factorp
c     scan_factordp partial sum of n_factordp
c     sum_factorp = sum(n_factorp)
c     max_facorp    maximum of p scaling factor for every atoms
c
c
#include "tinker_precision.h"
      module polgrp
      use sizes
      implicit none
      integer maxp11,maxp12
      integer maxp13,maxp14
      integer ninteract_scaling_p
      integer ,pointer :: np11(:),ip11(:,:),np12(:),ip12(:,:)
      integer ,pointer :: np13(:),ip13(:,:),np14(:),ip14(:,:)
      integer winnp11,winip11,winnp12,winip12
      integer winnp13,winip13,winnp14,winip14

      integer   ,pointer :: allscal_p(:),numscal_p(:),scalbeg_p(:)
      integer(1),pointer :: typscal_p(:)
      integer winallscal_p,winnumscal_p,winscalbeg_p,wintypscal_p
      integer   ,allocatable::pair_factorp(:,:)
      real(t_p) ,allocatable::scal_factorp(:,:)
      integer   ,allocatable::   n_factorp(:),   n_factordp(:)
      integer   ,allocatable::scan_factorp(:),scan_factordp(:)
      integer sum_factorp
      integer max_factorp

      !TODO Increase parameter to 120 following 1.2 version
      parameter (maxp11=120,
     &           maxp12=120,
     &           maxp13=120,
     &           maxp14=120)

      ! TODO Fix Analyze: Remove this
!$acc declare create(np11,ip11,np12,ip12)
!$acc declare create(np13,ip13,np14,ip14)
      end
