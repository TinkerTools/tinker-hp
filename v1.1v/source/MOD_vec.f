c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ###################################################################
c     ##                                                               ##
c     ##  module  vec  --  General Stuff for vectorization purpose     ##
c     ##                                                               ##
c     ###################################################################
c
c
      module vec
      use sizes
      use couple
      use polgrp
c     use domdec
      implicit none
!DIR$ ATTRIBUTES ALIGN:64::itmp14
      integer itmp14(maxn14)
!DIR$ ATTRIBUTES ALIGN:64:: iptmp11
      integer iptmp11(maxp11)
!DIR$ ATTRIBUTES ALIGN:64:: kpolevec,kpolevec1
      integer kpolevec(maxvlst),kpolevec1(maxvlst)
!DIR$ ATTRIBUTES ALIGN:64::kglobvec,kglobvec1,kglobvec2
      integer kglobvec(maxvlst),kglobvec1(maxvlst),kglobvec2(maxvlst)
!DIR$ ATTRIBUTES ALIGN:64::kbisvec,kbisvec1,kbisvec2
      integer kbisvec(maxvlst),kbisvec1(maxvlst),kbisvec2(maxvlst)
!DIR$ ATTRIBUTES ALIGN:64::xposvec,yposvec,zposvec
      real*8 xposvec(maxvlst),yposvec(maxvlst),zposvec(maxvlst)
!DIR$ ATTRIBUTES ALIGN:64::xposvec1,yposvec1,zposvec1
      real*8 xposvec1(maxvlst),yposvec1(maxvlst),zposvec1(maxvlst)
!DIR$ ATTRIBUTES ALIGN:64:: term1vec ,term2vec ,term3vec
      real*8 term1vec(maxvlst) ,term2vec(maxvlst) ,term3vec(maxvlst)
!DIR$ ATTRIBUTES ALIGN:64:: term4vec, term5vec,term6vec
      real*8 term4vec(maxvlst),term5vec(maxvlst),term6vec(maxvlst)
!DIR$ ATTRIBUTES ALIGN:64:: mask, mask1,maskloop
      logical mask(maxvlst), mask1(maxvlst),maskloop(maxvlst)
      end
