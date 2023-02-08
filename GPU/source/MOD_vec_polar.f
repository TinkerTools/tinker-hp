c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ###################################################################
c     ##                                                               ##
c     ##  module  vec_polar  --    Stuff for vectorization purpose     ##
c     ##                                                               ##
c     ###################################################################

      module vec_polar
      use sizes
      use couple
      use vec
      use vec_elec
      implicit none
#include "tinker_macro.h"
cold  !DIR$ ATTRIBUTES ALIGN:64:: kpolevec
cold  integer kpolevec(maxelst)
cold  !DIR$ ATTRIBUTES ALIGN:64::kpolevec1
cold  integer kpolevec1(maxelst)
cold !DIR$ ATTRIBUTES ALIGN:64::kpolevec2
cold !     integer kpolevec2(maxelst)
!DIR$ ATTRIBUTES ALIGN:64::iuscalevec
!DIR$ ATTRIBUTES ALIGN:64::uivec,uipvec
      real(t_p) uivec(3),uipvec(3)
      integer iuscalevec(maxelst)
!DIR$ ATTRIBUTES ALIGN:64::rdampvec
      real(t_p),dimension(maxelst):: rdampvec
!DIR$ ATTRIBUTES ALIGN:64::qritmp,qrktmp
      real(t_p),dimension(maxelst):: qritmp,qrktmp
!DIR$ ATTRIBUTES ALIGN:64::r2vec,invr2vec,rvec
      real(t_p),dimension(maxelst):: r2vec,invr2vec,rvec
!DIR$ ATTRIBUTES ALIGN:64:: invrvec
      real(t_p),dimension(maxelst):: invrvec
!DIR$ ATTRIBUTES ALIGN:64::uripvec
      real(t_p),dimension(maxelst):: uripvec
!DIR$ ATTRIBUTES ALIGN:64::ukvecx,ukvecy,ukvecz
      real(t_p),dimension(maxelst):: ukvecx,ukvecy,ukvecz
!DIR$ ATTRIBUTES ALIGN:64::ukpvecx,ukpvecy,ukpvecz
      real(t_p),dimension(maxelst):: ukpvecx,ukpvecy,ukpvecz
!DIR$ ATTRIBUTES ALIGN:64::urkpvec
      real(t_p),dimension(maxelst):: urkpvec
!DIR$ ATTRIBUTES ALIGN:64:: drivec,drkvec
      real(t_p),dimension(maxelst):: drivec,drkvec
!DIR$ ATTRIBUTES ALIGN:64:: qrivecx, qrivecy, qrivecz
      real(t_p),dimension(maxelst):: qrivecx,qrivecy,qrivecz
!DIR$ ATTRIBUTES ALIGN:64:: qrkvecx, qrkvecy, qrkvecz
      real(t_p),dimension(maxelst):: qrkvecx,qrkvecy,qrkvecz
!DIR$ ATTRIBUTES ALIGN:64:: qrimodvecx, qrimodvecy, qrimodvecz
      real(t_p),dimension(maxelst):: qrimodvecx,qrimodvecy,
     &                                          qrimodvecz
!DIR$ ATTRIBUTES ALIGN:64:: qrkmodvecx, qrkmodvecy, qrkmodvecz
      real(t_p),dimension(maxelst):: qrkmodvecx,qrkmodvecy,
     &                                          qrkmodvecz
!DIR$ ATTRIBUTES ALIGN:64:: qrrivec,qrrkvec
      real(t_p),dimension(maxelst):: qrrivec,qrrkvec
!DIR$ ATTRIBUTES ALIGN:64:: urivec,urkvec
      real(t_p),dimension(maxelst):: urivec,urkvec
!DIR$ ATTRIBUTES ALIGN:64:: duikvec,quikvec
      real(t_p),dimension(maxelst):: duikvec,quikvec
!DIR$ ATTRIBUTES ALIGN:64:: rr7vec,rr9vec,rr11vec
      real(t_p),dimension(maxelst):: rr7vec,rr9vec,rr11vec
!DIR$ ATTRIBUTES ALIGN:64:: rc3vecx ,rc3vecy ,rc3vecz
      real(t_p),dimension(maxelst):: rc3vecx,rc3vecy,rc3vecz
!DIR$ ATTRIBUTES ALIGN:64:: rc5vecx ,rc5vecy,rc5vecz
      real(t_p),dimension(maxelst):: rc5vecx,rc5vecy,rc5vecz
!DIR$ ATTRIBUTES ALIGN:64:: rc7vecx ,rc7vecy,rc7vecz
      real(t_p),dimension(maxelst):: rc7vecx,rc7vecy,rc7vecz
!DIR$ ATTRIBUTES ALIGN:64:: dsr3vec ,dsr5vec,dsr7vec
      real(t_p),dimension(maxelst):: dsr3vec,dsr5vec,dsr7vec
!DIR$ ATTRIBUTES ALIGN:64:: usr3vec,usr5vec
      real(t_p),dimension(maxelst):: usr3vec,usr5vec
!DIR$ ATTRIBUTES ALIGN:64:: usc3vec,usc5vec
      real(t_p),dimension(maxelst):: usc3vec,usc5vec
!DIR$ ATTRIBUTES ALIGN:64:: uscalevec
      real(t_p),dimension(maxelst):: uscalevec
!DIR$ ATTRIBUTES ALIGN:64:: prc3vecx ,prc3vecy,prc3vecz
      real(t_p),dimension(maxelst):: prc3vecx,prc3vecy,prc3vecz
!DIR$ ATTRIBUTES ALIGN:64:: prc5vecx ,prc5vecy,prc5vecz
      real(t_p),dimension(maxelst):: prc5vecx,prc5vecy,prc5vecz
!DIR$ ATTRIBUTES ALIGN:64:: prc7vecx ,prc7vecy,prc7vecz
      real(t_p),dimension(maxelst):: prc7vecx,prc7vecy,prc7vecz
!DIR$ ATTRIBUTES ALIGN:64:: urc3vecx ,urc3vecy ,urc3vecz
      real(t_p),dimension(maxelst):: urc3vecx,urc3vecy,urc3vecz
!DIR$ ATTRIBUTES ALIGN:64:: urc5vecx ,urc5vecy ,urc5vecz
      real(t_p),dimension(maxelst):: urc5vecx,urc5vecy,urc5vecz
!DIR$ ATTRIBUTES ALIGN:64:: drc3vecx ,drc3vecy ,drc3vecz
      real(t_p),dimension(maxelst):: drc3vecx,drc3vecy,drc3vecz
!DIR$ ATTRIBUTES ALIGN:64:: drc5vecx ,drc5vecy ,drc5vecz
      real(t_p),dimension(maxelst):: drc5vecx,drc5vecy,drc5vecz
!DIR$ ATTRIBUTES ALIGN:64:: drc7vecx ,drc7vecy ,drc7vecz
      real(t_p),dimension(maxelst):: drc7vecx,drc7vecy,drc7vecz
!DIR$ ATTRIBUTES ALIGN:64:: psr3vec ,psr5vec,psr7vec
      real(t_p),dimension(maxelst):: psr3vec,psr5vec,psr7vec
!DIR$ ATTRIBUTES ALIGN:64:: dterm1 ,dterm2 
      real(t_p),dimension(maxelst):: dterm1,dterm2
!DIR$ ATTRIBUTES ALIGN:64:: dterm3x,dterm3y,dterm3z
      real(t_p),dimension(maxelst):: dterm3x,dterm3y,dterm3z
!DIR$ ATTRIBUTES ALIGN:64:: dterm4x,dterm4y,dterm4z
      real(t_p),dimension(maxelst):: dterm4x,dterm4y,dterm4z
!DIR$ ATTRIBUTES ALIGN:64:: dterm5x,dterm5y,dterm5z
      real(t_p),dimension(maxelst):: dterm5x,dterm5y,dterm5z
!DIR$ ATTRIBUTES ALIGN:64:: dterm6x,dterm6y,dterm6z
      real(t_p),dimension(maxelst):: dterm6x,dterm6y,dterm6z
!DIR$ ATTRIBUTES ALIGN:64:: dterm7x,dterm7y,dterm7z
      real(t_p),dimension(maxelst):: dterm7x,dterm7y,dterm7z
!DIR$ ATTRIBUTES ALIGN:64:: tmp1vecx,tmp1vecy,tmp1vecz
      real(t_p),dimension(maxelst):: tmp1vecx,tmp1vecy,tmp1vecz
!DIR$ ATTRIBUTES ALIGN:64:: tmp2vecx,tmp2vecy,tmp2vecz
      real(t_p),dimension(maxelst):: tmp2vecx,tmp2vecy,tmp2vecz
!DIR$ ATTRIBUTES ALIGN:64:: tmp3vecx,tmp3vecy,tmp3vecz
      real(t_p),dimension(maxelst):: tmp3vecx,tmp3vecy,tmp3vecz
!DIR$ ATTRIBUTES ALIGN:64:: tmp4vecx,tmp4vecy,tmp4vecz
      real(t_p),dimension(maxelst):: tmp4vecx,tmp4vecy,tmp4vecz
!DIR$ ATTRIBUTES ALIGN:64:: tmp5vecx,tmp5vecy,tmp5vecz
      real(t_p),dimension(maxelst):: tmp5vecx,tmp5vecy,tmp5vecz
!DIR$ ATTRIBUTES ALIGN:64:: tmp6vecx,tmp6vecy,tmp6vecz
      real(t_p),dimension(maxelst):: tmp6vecx,tmp6vecy,tmp6vecz
!DIR$ ATTRIBUTES ALIGN:64:: tmp7vecx,tmp7vecy,tmp7vecz
      real(t_p),dimension(maxelst):: tmp7vecx,tmp7vecy,tmp7vecz
!DIR$ ATTRIBUTES ALIGN:64:: davec
      real(t_p),dimension(maxelst):: davec
!DIR$ ATTRIBUTES ALIGN:64:: tisvecx ,tisvecy,tisvecz
      real(t_p),dimension(maxelst):: tisvecx,tisvecy,tisvecz
!DIR$ ATTRIBUTES ALIGN:64:: tksvecx ,tksvecy,tksvecz
      real(t_p),dimension(maxelst):: tksvecx,tksvecy,tksvecz
!DIR$ ATTRIBUTES ALIGN:64:: ticvecx ,ticvecy,ticvecz
      real(t_p),dimension(maxelst):: ticvecx,ticvecy,ticvecz
!DIR$ ATTRIBUTES ALIGN:64:: tkcvecx ,tkcvecy,tkcvecz
      real(t_p),dimension(maxelst):: tkcvecx,tkcvecy,tkcvecz
!DIR$ ATTRIBUTES ALIGN:64:: tivec1 ,tivec2,tivec3
      real(t_p),dimension(maxelst):: tivec1,tivec2,tivec3
!DIR$ ATTRIBUTES ALIGN:64:: tivec4 ,tivec5,tivec6
      real(t_p),dimension(maxelst):: tivec4,tivec5,tivec6
!DIR$ ATTRIBUTES ALIGN:64:: tivec7 ,tivec8,tivec9
      real(t_p),dimension(maxelst):: tivec7,tivec8,tivec9
!DIR$ ATTRIBUTES ALIGN:64:: tkvec1 ,tkvec2,tkvec3
      real(t_p),dimension(maxelst):: tkvec1,tkvec2,tkvec3
!DIR$ ATTRIBUTES ALIGN:64:: tkvec4 ,tkvec5,tkvec6
      real(t_p),dimension(maxelst):: tkvec4,tkvec5,tkvec6
!DIR$ ATTRIBUTES ALIGN:64:: tkvec7 ,tkvec8,tkvec9
      real(t_p),dimension(maxelst):: tkvec7,tkvec8,tkvec9
!DIR$ ATTRIBUTES ALIGN:64:: depxvec ,depyvec,depzvec
      real(t_p),dimension(maxelst):: depxvec,depyvec,depzvec
!DIR$ ATTRIBUTES ALIGN:64:: frcxvec,frcyvec,frczvec
      real(t_p),dimension(maxelst):: frcxvec,frcyvec,frczvec
!DIR$ ATTRIBUTES ALIGN:64::ti3vecx,ti3vecy,ti3vecz
      real(t_p),dimension(maxelst):: ti3vecx,ti3vecy,ti3vecz
!DIR$ ATTRIBUTES ALIGN:64::tk3vecx,tk3vecy,tk3vecz
      real(t_p),dimension(maxelst):: tk3vecx,tk3vecy,tk3vecz
!DIR$ ATTRIBUTES ALIGN:64::ti5vecx,ti5vecy,ti5vecz
      real(t_p),dimension(maxelst):: ti5vecx,ti5vecy,ti5vecz
!DIR$ ATTRIBUTES ALIGN:64::tk5vecx,tk5vecy,tk5vecz
      real(t_p),dimension(maxelst):: tk5vecx,tk5vecy,tk5vecz
!DIR$ ATTRIBUTES ALIGN:64::turi5vec,turk5vec
      real(t_p),dimension(maxelst):: turi5vec,turk5vec
!DIR$ ATTRIBUTES ALIGN:64::turi7vec,turk7vec
      real(t_p),dimension(maxelst):: turi7vec,turk7vec
!DIR$ ATTRIBUTES ALIGN:64::efullvec
      real(t_p),dimension(maxelst):: efullvec
      end module
