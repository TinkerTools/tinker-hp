c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ###################################################################
c     ##                                                               ##
c     ##   module  vec_charge  --    Stuff for vectorization purpose   ##
c     ##                                                               ##
c     ###################################################################

      module vec_charge
      use sizes
      use couple
      use vec
      use vec_elec
      implicit none
#include "tinker_precision.h"

!PGI$ ATTRIBUTES ALIGN:64::kvec,kvec1,kvec2
      integer kvec(maxelst),kvec1(maxelst),kvec2(maxelst)
!PGI$ ATTRIBUTES ALIGN:64::kkchgvec,kkchgvec1
      integer kkchgvec(maxelst),kkchgvec1(maxelst)
!PGI$ ATTRIBUTES ALIGN:64:: kkchgvec2,kkchgvec3
      integer kkchgvec2(maxelst), kkchgvec3(maxelst)
!PGI$ ATTRIBUTES ALIGN:64::iscalevec
      integer iscalevec(maxelst)
!PGI$ ATTRIBUTES ALIGN:64::rvec,r2vec,invrvec
      real(t_p)  rvec(maxelst),r2vec(maxelst),invrvec(maxelst)
!PGI$ ATTRIBUTES ALIGN:64:: rewvec
      real(t_p)  rewvec(maxelst)
!PGI$ ATTRIBUTES ALIGN:64::fikvec,scalevec
      real(t_p)  fikvec(maxelst),scalevec(maxelst)
!PGI$ ATTRIBUTES ALIGN:64::invrbvec,invrb2vec
      real(t_p)  invrbvec(maxelst),invrb2vec(maxelst)
!PGI$ ATTRIBUTES ALIGN:64:: devec,evec,efullvec
      real(t_p)  devec(maxelst),evec(maxelst),efullvec(maxelst)
!PGI$ ATTRIBUTES ALIGN:64:: erfvec
      real(t_p)  erfvec(maxelst)
!PGI$ ATTRIBUTES ALIGN:64:: decxvec
      real(t_p)  decxvec(maxelst)
!PGI$ ATTRIBUTES ALIGN:64:: decyvec
      real(t_p)  decyvec(maxelst)
!PGI$ ATTRIBUTES ALIGN:64:: deczvec
      real(t_p)  deczvec(maxelst)
      end module
