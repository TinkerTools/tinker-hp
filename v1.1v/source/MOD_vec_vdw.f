c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ###################################################################
c     ##                                                               ##
c     ##  module  vec_elec  --    Stuff for vectorization purpose      ##
c     ##                                                               ##
c     ###################################################################
c
c
      module vec_vdw
      use sizes
      use couple
c     use domdec
      implicit none
!DIR$ ATTRIBUTES ALIGN:64::kvvec,kvvec1,kvvec2
      integer kvvec(maxvlst),kvvec1(maxvlst),kvvec2(maxvlst)
!DIR$ ATTRIBUTES ALIGN:64:: kvlocvec,kvlocvec1,kvlocvec2
      integer kvlocvec(maxvlst),kvlocvec1(maxvlst),kvlocvec2(maxvlst)
!DIR$ ATTRIBUTES ALIGN:64:: ktvec,ktvec2
      integer ktvec(maxvlst),ktvec2(maxvlst)
!DIR$ ATTRIBUTES ALIGN:64::vscalevec
      real*8 vscalevec(maxvlst)
!DIR$ ATTRIBUTES ALIGN:64::radminvec
      real*8 radminvec(maxvlst)
!DIR$ ATTRIBUTES ALIGN:64::radmin4vec
      real*8 radmin4vec(maxvlst)
!DIR$ ATTRIBUTES ALIGN:64::epsilonvec
      real*8 epsilonvec(maxvlst)
!DIR$ ATTRIBUTES ALIGN:64::epsilon4vec
      real*8 epsilon4vec(maxvlst)
!DIR$ ATTRIBUTES ALIGN:64::rvvec,epsvec,rvvec2
      real*8 rvvec(maxvlst),rvvec2(maxvlst),epsvec(maxvlst)
!DIR$ ATTRIBUTES ALIGN:64::rv7vec,invrhovec,epsvec2
      real*8 rv7vec(maxvlst),invrhovec(maxvlst),epsvec2(maxvlst)
!DIR$ ATTRIBUTES ALIGN:64::rik2vec,rik6vec,invrikvec
      real*8 rik2vec(maxvlst),rik6vec(maxvlst),invrikvec(maxvlst)
!DIR$ ATTRIBUTES ALIGN:64::rik2vec1
      real*8 rik2vec1(maxvlst)
!DIR$ ATTRIBUTES ALIGN:64::rikvec,rik3vec,rik4vec
      real*8 rikvec(maxvlst),rik3vec(maxvlst),rik4vec(maxvlst)
!DIR$ ATTRIBUTES ALIGN:64::rik7vec,invtmpvec
      real*8 rik7vec(maxvlst),invtmpvec(maxvlst)
!DIR$ ATTRIBUTES ALIGN:64::tapervec,evec,rik5vec
      real*8 tapervec(maxvlst),evec(maxvlst),rik5vec(maxvlst)
!DIR$ ATTRIBUTES ALIGN:64:: devec,dtapervec,rv7orhovec
      real*8 devec(maxvlst),dtapervec(maxvlst),rv7orhovec(maxvlst)
!DIR$ ATTRIBUTES ALIGN:64:: dedxvec,dedyvec,dedzvec
      real*8 dedxvec(maxvlst), dedyvec(maxvlst), dedzvec(maxvlst)
!DIR$ ATTRIBUTES ALIGN:64:: redkvec,dtauvec
      real*8 redkvec(maxvlst),dtauvec(maxvlst)
!DIR$ ATTRIBUTES ALIGN:64:: tauvec,gtauvec,tau7vec
      real*8 tauvec(maxvlst),gtauvec(maxvlst),tau7vec(maxvlst)
!DIR$ ATTRIBUTES ALIGN:64:: redivec,p6vec
      real*8 redivec(maxvlst),p6vec(maxvlst)
!DIR$ ATTRIBUTES ALIGN:64:: devxvec,devyvec,devzvec
      real*8 devxvec(maxvlst),devyvec(maxvlst),devzvec(maxvlst)
!DIR$ ATTRIBUTES ALIGN:64:: devvxvec,devvyvec,devvzvec
      real*8 devvxvec(maxvlst),devvyvec(maxvlst),devvzvec(maxvlst)

      end
