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
cold  !DIR$ ATTRIBUTES ALIGN:64:: kpolevec
cold  integer kpolevec(maxelst)
cold  !DIR$ ATTRIBUTES ALIGN:64::kpolevec1
cold  integer kpolevec1(maxelst)
cold  !DIR$ ATTRIBUTES ALIGN:64::kpolevec2
cold  integer kpolevec2(maxelst)
      !DIR$ ATTRIBUTES ALIGN:64::rdampvec
      real*8 rdampvec(maxelst)
      !DIR$ ATTRIBUTES ALIGN:64::qritmp,qrktmp
      real*8 qritmp(maxelst),qrktmp(maxelst)
      !DIR$ ATTRIBUTES ALIGN:64::r2vec,invr2vec,rvec
      real*8 r2vec(maxelst),invr2vec(maxelst),rvec(maxelst)
      !DIR$ ATTRIBUTES ALIGN:64:: invrvec
      real*8 invrvec(maxelst)
      !DIR$ ATTRIBUTES ALIGN:64::uivec,uipvec,uripvec
      real*8 uivec(3),uipvec(3),uripvec(maxelst)
      !DIR$ ATTRIBUTES ALIGN:64::ukvecx,ukvecy,ukvecz
      real*8 ukvecx(maxelst),ukvecy(maxelst),ukvecz(maxelst)
      !DIR$ ATTRIBUTES ALIGN:64::ukpvecx,ukpvecy,ukpvecz
      real*8 ukpvecx(maxelst),ukpvecy(maxelst),ukpvecz(maxelst)
      !DIR$ ATTRIBUTES ALIGN:64::urkpvec
      real*8 urkpvec(maxelst)
      !DIR$ ATTRIBUTES ALIGN:64:: drivec,drkvec
      real*8 drivec(maxelst),drkvec(maxelst)
      !DIR$ ATTRIBUTES ALIGN:64:: qrivecx, qrivecy, qrivecz
      real*8 qrivecx(maxelst),qrivecy(maxelst),qrivecz(maxelst)
      !DIR$ ATTRIBUTES ALIGN:64:: qrkvecx, qrkvecy, qrkvecz
      real*8 qrkvecx(maxelst),qrkvecy(maxelst),qrkvecz(maxelst)
      !DIR$ ATTRIBUTES ALIGN:64:: qrimodvecx, qrimodvecy, qrimodvecz
      real*8 qrimodvecx(maxelst),qrimodvecy(maxelst),qrimodvecz(maxelst)
      !DIR$ ATTRIBUTES ALIGN:64:: qrkmodvecx, qrkmodvecy, qrkmodvecz
      real*8 qrkmodvecx(maxelst),qrkmodvecy(maxelst),qrkmodvecz(maxelst)
      !DIR$ ATTRIBUTES ALIGN:64:: qrrivec,qrrkvec
      real*8 qrrivec(maxelst),qrrkvec(maxelst)
      !DIR$ ATTRIBUTES ALIGN:64:: urivec,urkvec
      real*8 urivec(maxelst),urkvec(maxelst)
      !DIR$ ATTRIBUTES ALIGN:64:: duikvec,quikvec
      real*8 duikvec(maxelst),quikvec(maxelst)
      !DIR$ ATTRIBUTES ALIGN:64:: rr7vec,rr9vec,rr11vec
      real*8 rr7vec(maxelst),rr9vec(maxelst),rr11vec(maxelst)
      !DIR$ ATTRIBUTES ALIGN:64:: rc3vecx ,rc3vecy ,rc3vecz
      real*8 rc3vecx(maxelst) ,rc3vecy(maxelst) ,rc3vecz(maxelst)
      !DIR$ ATTRIBUTES ALIGN:64:: rc5vecx ,rc5vecy ,rc5vecz
      real*8 rc5vecx(maxelst) ,rc5vecy(maxelst) ,rc5vecz(maxelst)
      !DIR$ ATTRIBUTES ALIGN:64:: rc7vecx ,rc7vecy ,rc7vecz
      real*8 rc7vecx(maxelst) ,rc7vecy(maxelst) ,rc7vecz(maxelst)
      !DIR$ ATTRIBUTES ALIGN:64:: dsr3vec ,dsr5vec ,dsr7vec
      real*8 dsr3vec(maxelst),dsr5vec(maxelst),dsr7vec(maxelst)
      !DIR$ ATTRIBUTES ALIGN:64:: usr3vec,usr5vec
      real*8 usr3vec(maxelst),usr5vec(maxelst)
      !DIR$ ATTRIBUTES ALIGN:64:: usc3vec,usc5vec
      real*8 usc3vec(maxelst),usc5vec(maxelst)
      !DIR$ ATTRIBUTES ALIGN:64:: uscalevec
      real*8 uscalevec(maxelst)
      !DIR$ ATTRIBUTES ALIGN:64:: prc3vecx ,prc3vecy,prc3vecz
      real*8 prc3vecx(maxelst) ,prc3vecy(maxelst) ,prc3vecz(maxelst)
      !DIR$ ATTRIBUTES ALIGN:64:: prc5vecx ,prc5vecy,prc5vecz
      real*8 prc5vecx(maxelst) ,prc5vecy(maxelst) ,prc5vecz(maxelst)
      !DIR$ ATTRIBUTES ALIGN:64:: prc7vecx ,prc7vecy,prc7vecz
      real*8 prc7vecx(maxelst) ,prc7vecy(maxelst) ,prc7vecz(maxelst)
      !DIR$ ATTRIBUTES ALIGN:64:: urc3vecx ,urc3vecy ,urc3vecz
      real*8 urc3vecx(maxelst) ,urc3vecy(maxelst),urc3vecz(maxelst)
      !DIR$ ATTRIBUTES ALIGN:64:: urc5vecx ,urc5vecy ,urc5vecz
      real*8 urc5vecx(maxelst) ,urc5vecy(maxelst),urc5vecz(maxelst)
      !DIR$ ATTRIBUTES ALIGN:64:: drc3vecx ,drc3vecy ,drc3vecz
      real*8 drc3vecx(maxelst) ,drc3vecy(maxelst) ,drc3vecz(maxelst)
      !DIR$ ATTRIBUTES ALIGN:64:: drc5vecx ,drc5vecy ,drc5vecz
      real*8 drc5vecx(maxelst) ,drc5vecy(maxelst) ,drc5vecz(maxelst)
      !DIR$ ATTRIBUTES ALIGN:64:: drc7vecx ,drc7vecy ,drc7vecz
      real*8 drc7vecx(maxelst) ,drc7vecy(maxelst) ,drc7vecz(maxelst)
      !DIR$ ATTRIBUTES ALIGN:64:: psr3vec ,psr5vec ,psr7vec
      real*8 psr3vec(maxelst) ,psr5vec(maxelst) ,psr7vec(maxelst)
      !DIR$ ATTRIBUTES ALIGN:64:: dterm1 ,dterm2 
      real*8 dterm1(maxelst),dterm2(maxelst)
      !DIR$ ATTRIBUTES ALIGN:64:: dterm3x,dterm3y,dterm3z
      real*8 dterm3x(maxelst),dterm3y(maxelst),dterm3z(maxelst)
      !DIR$ ATTRIBUTES ALIGN:64:: dterm4x,dterm4y,dterm4z
      real*8 dterm4x(maxelst),dterm4y(maxelst),dterm4z(maxelst)
      !DIR$ ATTRIBUTES ALIGN:64:: dterm5x,dterm5y,dterm5z
      real*8 dterm5x(maxelst),dterm5y(maxelst),dterm5z(maxelst)
      !DIR$ ATTRIBUTES ALIGN:64:: dterm6x,dterm6y,dterm6z
      real*8 dterm6x(maxelst),dterm6y(maxelst),dterm6z(maxelst)
      !DIR$ ATTRIBUTES ALIGN:64:: dterm7x,dterm7y,dterm7z
      real*8 dterm7x(maxelst),dterm7y(maxelst),dterm7z(maxelst)
      !DIR$ ATTRIBUTES ALIGN:64:: tmp1vecx,tmp1vecy,tmp1vecz
      real*8 tmp1vecx(maxelst),tmp1vecy(maxelst),tmp1vecz(maxelst)
      !DIR$ ATTRIBUTES ALIGN:64:: tmp2vecx,tmp2vecy,tmp2vecz
      real*8 tmp2vecx(maxelst),tmp2vecy(maxelst),tmp2vecz(maxelst)
      !DIR$ ATTRIBUTES ALIGN:64:: tmp3vecx,tmp3vecy,tmp3vecz
      real*8 tmp3vecx(maxelst),tmp3vecy(maxelst),tmp3vecz(maxelst)
      !DIR$ ATTRIBUTES ALIGN:64:: tmp4vecx,tmp4vecy,tmp4vecz
      real*8 tmp4vecx(maxelst),tmp4vecy(maxelst),tmp4vecz(maxelst)
      !DIR$ ATTRIBUTES ALIGN:64:: tmp5vecx,tmp5vecy,tmp5vecz
      real*8 tmp5vecx(maxelst),tmp5vecy(maxelst),tmp5vecz(maxelst)
      !DIR$ ATTRIBUTES ALIGN:64:: tmp6vecx,tmp6vecy,tmp6vecz
      real*8 tmp6vecx(maxelst),tmp6vecy(maxelst),tmp6vecz(maxelst)
      !DIR$ ATTRIBUTES ALIGN:64:: tmp7vecx,tmp7vecy,tmp7vecz
      real*8 tmp7vecx(maxelst),tmp7vecy(maxelst),tmp7vecz(maxelst)
      !DIR$ ATTRIBUTES ALIGN:64:: davec
      real*8 davec(maxelst)
      !DIR$ ATTRIBUTES ALIGN:64:: tisvecx ,tisvecy,tisvecz
      real*8 tisvecx(maxelst),tisvecy(maxelst),tisvecz(maxelst)
      !DIR$ ATTRIBUTES ALIGN:64:: tksvecx ,tksvecy,tksvecz
      real*8 tksvecx(maxelst),tksvecy(maxelst),tksvecz(maxelst)
      !DIR$ ATTRIBUTES ALIGN:64:: ticvecx ,ticvecy,ticvecz
      real*8 ticvecx(maxelst),ticvecy(maxelst),ticvecz(maxelst)
      !DIR$ ATTRIBUTES ALIGN:64:: tkcvecx ,tkcvecy,tkcvecz
      real*8 tkcvecx(maxelst),tkcvecy(maxelst),tkcvecz(maxelst)
      !DIR$ ATTRIBUTES ALIGN:64:: tivec1 ,tivec2,tivec3
      real*8 tivec1(maxelst),tivec2(maxelst),tivec3(maxelst)
      !DIR$ ATTRIBUTES ALIGN:64:: tivec4 ,tivec5,tivec6
      real*8 tivec4(maxelst),tivec5(maxelst),tivec6(maxelst)
      !DIR$ ATTRIBUTES ALIGN:64:: tivec7 ,tivec8,tivec9
      real*8 tivec7(maxelst),tivec8(maxelst),tivec9(maxelst)
      !DIR$ ATTRIBUTES ALIGN:64:: tkvec1 ,tkvec2,tkvec3
      real*8 tkvec1(maxelst),tkvec2(maxelst),tkvec3(maxelst)
      !DIR$ ATTRIBUTES ALIGN:64:: tkvec4 ,tkvec5,tkvec6
      real*8 tkvec4(maxelst),tkvec5(maxelst),tkvec6(maxelst)
      !DIR$ ATTRIBUTES ALIGN:64:: tkvec7 ,tkvec8,tkvec9
      real*8 tkvec7(maxelst),tkvec8(maxelst),tkvec9(maxelst)
      !DIR$ ATTRIBUTES ALIGN:64:: depxvec ,depyvec,depzvec
      real*8 depxvec(maxelst),depyvec(maxelst),depzvec(maxelst)
      !DIR$ ATTRIBUTES ALIGN:64:: frcxvec,frcyvec,frczvec
      real*8 frcxvec(maxelst),frcyvec(maxelst),frczvec(maxelst)
      !DIR$ ATTRIBUTES ALIGN:64::ti3vecx,ti3vecy,ti3vecz
      real*8 ti3vecx(maxelst),ti3vecy(maxelst),ti3vecz(maxelst)
      !DIR$ ATTRIBUTES ALIGN:64::tk3vecx,tk3vecy,tk3vecz
      real*8 tk3vecx(maxelst),tk3vecy(maxelst),tk3vecz(maxelst)
      !DIR$ ATTRIBUTES ALIGN:64::ti5vecx,ti5vecy,ti5vecz
      real*8 ti5vecx(maxelst),ti5vecy(maxelst),ti5vecz(maxelst)
      !DIR$ ATTRIBUTES ALIGN:64::tk5vecx,tk5vecy,tk5vecz
      real*8 tk5vecx(maxelst),tk5vecy(maxelst),tk5vecz(maxelst)
      !DIR$ ATTRIBUTES ALIGN:64::turi5vec,turk5vec
      real*8 turi5vec(maxelst),turk5vec(maxelst)
      !DIR$ ATTRIBUTES ALIGN:64::turi7vec,turk7vec
      real*8 turi7vec(maxelst),turk7vec(maxelst)
      !DIR$ ATTRIBUTES ALIGN:64::efullvec
      real*8 efullvec(maxelst)
      end module
