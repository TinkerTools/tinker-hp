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
      module vec_elec
      use sizes
      use couple
      use polgrp
      use vec
      implicit none
!DIR$ ATTRIBUTES ALIGN:64:: kpolelocvec
      integer kpolelocvec(maxelst)
!DIR$ ATTRIBUTES ALIGN:64:: kpolelocvec1
      integer kpolelocvec1(maxelst)
!DIR$ ATTRIBUTES ALIGN:64:: kpolelocvec2
      integer kpolelocvec2(maxelst)
!!!DIR$ ATTRIBUTES ALIGN:64:: d2vec,d2vec1,d1vec
!!!     real*8 d2vec(maxelst),d2vec1(maxelst),d1vec(maxelst)
!DIR$ ATTRIBUTES ALIGN:64:: d1vec
      real*8 d1vec(maxelst)
!DIR$ ATTRIBUTES ALIGN:64:: invd1vec,invd2vec
      real*8 invd1vec(maxelst),invd2vec(maxelst)
!!DIR$ ATTRIBUTES ALIGN:64:: invd3vec,invd5vec,invd7vec
!!     real*8 invd3vec(maxelst),invd5vec(maxelst),invd7vec(maxelst)
!DIR$ ATTRIBUTES ALIGN:64:: duivec,puivec
      real*8 duivec(3),puivec(3)
!DIR$ ATTRIBUTES ALIGN:64:: duivecx
      real*8 duivecx(maxelst)
!DIR$ ATTRIBUTES ALIGN:64:: puivecx
      real*8 puivecx(maxelst)
!DIR$ ATTRIBUTES ALIGN:64:: duivecy
      real*8 duivecy(maxelst)
!DIR$ ATTRIBUTES ALIGN:64:: puivecy
      real*8 puivecy(maxelst)
!DIR$ ATTRIBUTES ALIGN:64:: duivecz
      real*8 duivecz(maxelst)
!DIR$ ATTRIBUTES ALIGN:64:: puivecz
      real*8 puivecz(maxelst)
!DIR$ ATTRIBUTES ALIGN:64:: dukvecx,dukvecy,dukvecz
      real*8 dukvecx(maxelst),dukvecy(maxelst),dukvecz(maxelst)
!DIR$ ATTRIBUTES ALIGN:64:: pukvecx,pukvecy,pukvecz
      real*8 pukvecx(maxelst),pukvecy(maxelst),pukvecz(maxelst)
!DIR$ ATTRIBUTES ALIGN:64:: duirvec,dukrvec
      real*8 duirvec(maxelst),dukrvec(maxelst)
!DIR$ ATTRIBUTES ALIGN:64:: puirvec,pukrvec
      real*8 puirvec(maxelst),pukrvec(maxelst)
!DIR$ ATTRIBUTES ALIGN:64:: rr1vec, rr3vec,rr5vec
      real*8 rr1vec(maxelst),rr3vec(maxelst),rr5vec(maxelst)
!DIR$ ATTRIBUTES ALIGN:64:: sc3vec,sc5vec,sc7vec
      real*8 sc3vec(maxelst),sc5vec(maxelst),sc7vec(maxelst)
!DIR$ ATTRIBUTES ALIGN:64:: ralphavec,exp2avec,bnvec
      real*8 ralphavec(maxelst),exp2avec(maxelst),bnvec(maxelst,0:5)
!DIR$ ATTRIBUTES ALIGN:64:: bn0vec
      real*8 bn0vec(maxelst)
!DIR$ ATTRIBUTES ALIGN:64:: bn1vec
      real*8 bn1vec(maxelst)
!DIR$ ATTRIBUTES ALIGN:64:: bn2vec
      real*8 bn2vec(maxelst)
!DIR$ ATTRIBUTES ALIGN:64:: bn3vec
      real*8 bn3vec(maxelst)
!DIR$ ATTRIBUTES ALIGN:64:: bn4vec
      real*8 bn4vec(maxelst)
!DIR$ ATTRIBUTES ALIGN:64:: bn5vec
      real*8 bn5vec(maxelst)
!DIR$ ATTRIBUTES ALIGN:64:: pgammavec
      real*8 pgammavec(maxelst)
!DIR$ ATTRIBUTES ALIGN:64:: dampvec,invdampvec
      real*8 dampvec(maxelst),invdampvec(maxelst)
!DIR$ ATTRIBUTES ALIGN:64:: dampvec1,expdampvec1
      real*8 dampvec1(maxelst),expdampvec1(maxelst)
!DIR$ ATTRIBUTES ALIGN:64:: tholevec
      real*8 tholevec(maxelst)
!DIR$ ATTRIBUTES ALIGN:64:: dscalevec,pscalevec
      real*8 dscalevec(maxelst),pscalevec(maxelst)
!DIR$ ATTRIBUTES ALIGN:64:: fimpvecx,fimpvecy,fimpvecz
      real*8 fimpvecx(maxelst),fimpvecy(maxelst),fimpvecz(maxelst)
!DIR$ ATTRIBUTES ALIGN:64:: fimdvecx,fimdvecy,fimdvecz
      real*8 fimdvecx(maxelst),fimdvecy(maxelst),fimdvecz(maxelst)
!DIR$ ATTRIBUTES ALIGN:64:: fidvecx,fidvecy,fidvecz
      real*8 fidvecx(maxelst),fidvecy(maxelst),fidvecz(maxelst)
!DIR$ ATTRIBUTES ALIGN:64:: fipvecx,fipvecy,fipvecz
      real*8 fipvecx(maxelst),fipvecy(maxelst),fipvecz(maxelst)
!DIR$ ATTRIBUTES ALIGN:64:: fimvecx,fimvecy,fimvecz
      real*8 fimvecx(maxelst),fimvecy(maxelst),fimvecz(maxelst)
!DIR$ ATTRIBUTES ALIGN:64:: fkmdvecx,fkmdvecy,fkmdvecz
      real*8 fkmdvecx(maxelst),fkmdvecy(maxelst),fkmdvecz(maxelst)
!DIR$ ATTRIBUTES ALIGN:64:: fkmpvecx,fkmpvecy,fkmpvecz
      real*8 fkmpvecx(maxelst),fkmpvecy(maxelst),fkmpvecz(maxelst)
!DIR$ ATTRIBUTES ALIGN:64:: fkdvecx,fkdvecy,fkdvecz
      real*8 fkdvecx(maxelst),fkdvecy(maxelst),fkdvecz(maxelst)
!DIR$ ATTRIBUTES ALIGN:64:: fkpvecx,fkpvecy,fkpvecz
      real*8 fkpvecx(maxelst),fkpvecy(maxelst),fkpvecz(maxelst)
!DIR$ ATTRIBUTES ALIGN:64:: fkmvecx,fkmvecy,fkmvecz
      real*8 fkmvecx(maxelst),fkmvecy(maxelst),fkmvecz(maxelst)
!DIR$ ATTRIBUTES ALIGN:64:: qirvecx, qirvecy,qirvecz
      real*8  qirvecx(maxelst),qirvecy(maxelst),qirvecz(maxelst)
!DIR$ ATTRIBUTES ALIGN:64:: qkrvecx, qkrvecy,qkrvecz
      real*8  qkrvecx(maxelst),qkrvecy(maxelst),qkrvecz(maxelst)
!DIR$ ATTRIBUTES ALIGN:64:: ckvec
      real*8  ckvec(maxelst)
!DIR$ ATTRIBUTES ALIGN:64:: dkvecx,dkvecy,dkvecz
      real*8  dkvecx(maxelst),dkvecy(maxelst),dkvecz(maxelst)
!DIR$ ATTRIBUTES ALIGN:64:: dirvec,dkrvec
      real*8  dirvec(maxelst),dkrvec(maxelst)
!DIR$ ATTRIBUTES ALIGN:64:: divec,qivec
      real*8  divec(3),qivec(9)
!DIR$ ATTRIBUTES ALIGN:64:: qkvec1,qkvec2,qkvec3 
      real*8  qkvec1(maxelst),  qkvec2(maxelst),  qkvec3(maxelst) 
!DIR$ ATTRIBUTES ALIGN:64:: qkvec4,qkvec5,qkvec6 
      real*8  qkvec4(maxelst),  qkvec5(maxelst),  qkvec6(maxelst) 
!DIR$ ATTRIBUTES ALIGN:64:: qkvec7,qkvec8,qkvec9 
      real*8  qkvec7(maxelst),  qkvec8(maxelst),  qkvec9(maxelst) 
!DIR$ ATTRIBUTES ALIGN:64:: qirrvec,qkrrvec
      real*8  qirrvec(maxelst),qkrrvec(maxelst)
!DIR$ ATTRIBUTES ALIGN:64:: dsc3vec,dsc5vec,dsc7vec
      real*8  dsc3vec(maxelst),dsc5vec(maxelst),dsc7vec(maxelst)
!DIR$ ATTRIBUTES ALIGN:64:: psc3vec,psc5vec,psc7vec
      real*8  psc3vec(maxelst),psc5vec(maxelst),psc7vec(maxelst)
!DIR$ ATTRIBUTES ALIGN:64:: drr3vec,drr5vec,drr7vec
      real*8  drr3vec(maxelst),drr5vec(maxelst),drr7vec(maxelst)
!DIR$ ATTRIBUTES ALIGN:64:: prr3vec,prr5vec,prr7vec
      real*8  prr3vec(maxelst),prr5vec(maxelst),prr7vec(maxelst)
!DIR$ ATTRIBUTES ALIGN:64::xpsvec,ypsvec,zpsvec
      real*8 xpsvec(maxelst),ypsvec(maxelst),zpsvec(maxelst)
      end
