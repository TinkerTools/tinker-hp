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

!DIR$ ATTRIBUTES ALIGN:64::kvec,kvec1,kvec2
      integer kvec(maxelst),kvec1(maxelst),kvec2(maxelst)
!DIR$ ATTRIBUTES ALIGN:64::kkchgvec,kkchgvec1
      integer kkchgvec(maxelst),kkchgvec1(maxelst)
!DIR$ ATTRIBUTES ALIGN:64:: kkchgvec2,kkchgvec3
      integer kkchgvec2(maxelst), kkchgvec3(maxelst)
!DIR$ ATTRIBUTES ALIGN:64::rvec,r2vec,invrvec
      real*8  rvec(maxelst),r2vec(maxelst),invrvec(maxelst)
!DIR$ ATTRIBUTES ALIGN:64:: rewvec
      real*8  rewvec(maxelst)
!DIR$ ATTRIBUTES ALIGN:64::fikvec,scalevec
      real*8  fikvec(maxelst),scalevec(maxelst)
!DIR$ ATTRIBUTES ALIGN:64::invrbvec,invrb2vec
      real*8  invrbvec(maxelst),invrb2vec(maxelst)
!DIR$ ATTRIBUTES ALIGN:64:: devec,evec,efullvec
      real*8  devec(maxelst),evec(maxelst),efullvec(maxelst)
!DIR$ ATTRIBUTES ALIGN:64:: erfvec
      real*8  erfvec(maxelst)
!DIR$ ATTRIBUTES ALIGN:64:: decxvec
      real*8  decxvec(maxelst)
!DIR$ ATTRIBUTES ALIGN:64:: decyvec
      real*8  decyvec(maxelst)
!DIR$ ATTRIBUTES ALIGN:64:: deczvec
      real*8  deczvec(maxelst)

      end module
