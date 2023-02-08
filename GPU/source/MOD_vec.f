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
#include "tinker_macro.h"
      module vec
      use sizes  ,only: maxvlst,maxelst
      use couple ,only: maxn14
      use polgrp ,only: maxp11
c     use domdec
      implicit none
!DIR$ ATTRIBUTES ALIGN:64::itmp14
      integer itmp14(maxn14)
!DIR$ ATTRIBUTES ALIGN:64:: iptmp11
      integer iptmp11(maxp11)
!DIR$ ATTRIBUTES ALIGN:64:: kpolevec,kpolevec1
      integer,dimension(maxvlst):: kpolevec,kpolevec1
!DIR$ ATTRIBUTES ALIGN:64:: kpolevec2
      integer kpolevec2(maxelst)
!DIR$ ATTRIBUTES ALIGN:64::kglobvec,kglobvec1,kglobvec2
      integer,dimension(maxvlst):: kglobvec,kglobvec1,kglobvec2
!DIR$ ATTRIBUTES ALIGN:64::kbisvec,kbisvec1,kbisvec2
      integer,dimension(maxvlst):: kbisvec,kbisvec1,kbisvec2
!DIR$ ATTRIBUTES ALIGN:64::xposvec,yposvec,zposvec
      real(t_p),dimension(maxvlst):: xposvec,yposvec,zposvec
!DIR$ ATTRIBUTES ALIGN:64::xposvec1,yposvec1,zposvec1
      real(t_p),dimension(maxvlst):: xposvec1,yposvec1,zposvec1
!DIR$ ATTRIBUTES ALIGN:64:: term1vec ,term2vec ,term3vec
      real(t_p),dimension(maxvlst):: term1vec,term2vec,term3vec
!DIR$ ATTRIBUTES ALIGN:64:: term4vec, term5vec,term6vec
      real(t_p),dimension(maxvlst):: term4vec,term5vec,term6vec
!DIR$ ATTRIBUTES ALIGN:64:: mask, mask1,maskloop
      logical  ,dimension(maxvlst):: mask, mask1,maskloop
      end
