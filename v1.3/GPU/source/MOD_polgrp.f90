!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     #####################################################################
!     ##                                                                 ##
!     ##  module polgrp  --  polarizable site group connectivity lists   ##
!     ##                                                                 ##
!     #####################################################################
!
!
!     maxp11   maximum number of atoms in a polarization group
!     maxp12   maximum number of atoms in groups 1-2 to an atom
!     maxp13   maximum number of atoms in groups 1-3 to an atom
!     maxp14   maximum number of atoms in groups 1-4 to an atom
!
!     np11     number of atoms in polarization group of each atom
!     winnp11    window object corresponding to np11
!     ip11     atom numbers of atoms in same group as each atom
!     winip11    window object corresponding to ip11
!     np12     number of atoms in groups 1-2 to each atom
!     winnp12    window object corresponding to np12
!     ip12     atom numbers of atoms in groups 1-2 to each atom
!     winip12    window object corresponding to ip12
!     np13     number of atoms in groups 1-3 to each atom
!     winnp13    window object corresponding to np13
!     ip13     atom numbers of atoms in groups 1-3 to each atom
!     winip13    window object corresponding to ip13
!     np14     number of atoms in groups 1-4 to each atom
!     winnp14    window object corresponding to np14
!     ip14     atom numbers of atoms in groups 1-4 to each atom
!     winip14    window object corresponding to ip14
!
!     allscal_p     list of all p scaling factor of the system
!     winallscal_p  window object corresponding to allscal_p
!     numscal_p     p scaling factor's number per atom
!     winnumscal_p  window object corresponding to allnumscal_p
!     scalbeg_p     index start point of p scaling factor in allscal_n
!     winscalbeg_p  window object corresponding to allscalbeg_p
!     typscal_p     pair to allscal_p to store p scaling factor type
!     wintypscal_p  window object corresponding to alltypscal_p
!     pair_factorp  pair list of p scaling factor
!     scal_facotrp  scaling factor value for each atoms
!     n_factorp     number of scaling useful scaling interaction for each atoms
!     n_factordp    number of scaling useful pair (d,p) scaling interaction for each atoms
!     scan_factorp  partial sum of n_factorp
!     scan_factordp partial sum of n_factordp
!     sum_factorp = sum(n_factorp)
!     max_facorp    maximum of p scaling factor for every atoms
!
!
#include "tinker_macro.h"
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

   parameter (maxp11=200,&
      &maxp12=200,&
      &maxp13=200,&
      &maxp14=200)

   ! TODO Fix Analyze: Remove this
!$acc declare create(np11,ip11,np12,ip12)
!$acc declare create(np13,ip13,np14,ip14)
end
