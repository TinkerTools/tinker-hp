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
#include "tinker_precision.h"
!PGI$ ATTRIBUTES ALIGN:64:: kpolelocvec
      integer kpolelocvec(maxelst)
!PGI$ ATTRIBUTES ALIGN:64:: kpolelocvec1
      integer kpolelocvec1(maxelst)
!PGI$ ATTRIBUTES ALIGN:64:: kpolelocvec2
      integer kpolelocvec2(maxelst)
      integer idscalevec(maxelst),ipscalevec(maxelst)
!PGI$ ATTRIBUTES ALIGN:64:: d1vec
      real(t_p) d1vec(maxelst)
!PGI$ ATTRIBUTES ALIGN:64:: invd1vec,invd2vec
      real(t_p),dimension(maxelst)::invd1vec,invd2vec
!PGI$ ATTRIBUTES ALIGN:64:: duivec,puivec
      real(t_p) duivec(3),puivec(3)
!PGI$ ATTRIBUTES ALIGN:64:: duivecx
      real(t_p) duivecx(maxelst)
!PGI$ ATTRIBUTES ALIGN:64:: puivecx
      real(t_p) puivecx(maxelst)
!PGI$ ATTRIBUTES ALIGN:64:: duivecy
      real(t_p) duivecy(maxelst)
!PGI$ ATTRIBUTES ALIGN:64:: puivecy
      real(t_p) puivecy(maxelst)
!PGI$ ATTRIBUTES ALIGN:64:: duivecz
      real(t_p) duivecz(maxelst)
!PGI$ ATTRIBUTES ALIGN:64:: puivecz
      real(t_p) puivecz(maxelst)
!PGI$ ATTRIBUTES ALIGN:64:: dukvecx,dukvecy,dukvecz
      real(t_p),dimension(maxelst)::dukvecx,dukvecy,dukvecz
!PGI$ ATTRIBUTES ALIGN:64:: pukvecx,pukvecy,pukvecz
      real(t_p),dimension(maxelst)::pukvecx,pukvecy,pukvecz
!PGI$ ATTRIBUTES ALIGN:64:: duirvec,dukrvec
      real(t_p) duirvec(maxelst),dukrvec(maxelst)
!PGI$ ATTRIBUTES ALIGN:64:: puirvec,pukrvec
      real(t_p) puirvec(maxelst),pukrvec(maxelst)
!PGI$ ATTRIBUTES ALIGN:64:: rr1vec, rr3vec,rr5vec
      real(t_p),dimension(maxelst)::rr1vec,rr3vec,rr5vec
!PGI$ ATTRIBUTES ALIGN:64:: sc3vec,sc5vec,sc7vec
      real(t_p),dimension(maxelst)::sc3vec,sc5vec,sc7vec
!PGI$ ATTRIBUTES ALIGN:64:: ralphavec,exp2avec,bnvec
      real(t_p),dimension(maxelst)::ralphavec,exp2avec
!PGI$ ATTRIBUTES ALIGN:64:: bnvec
      real(t_p) bnvec(maxelst,0:5)
!PGI$ ATTRIBUTES ALIGN:64:: bn0vec
      real(t_p) bn0vec(maxelst)
!PGI$ ATTRIBUTES ALIGN:64:: bn1vec
      real(t_p) bn1vec(maxelst)
!PGI$ ATTRIBUTES ALIGN:64:: bn2vec
      real(t_p) bn2vec(maxelst)
!PGI$ ATTRIBUTES ALIGN:64:: bn3vec
      real(t_p) bn3vec(maxelst)
!PGI$ ATTRIBUTES ALIGN:64:: bn4vec
      real(t_p) bn4vec(maxelst)
!PGI$ ATTRIBUTES ALIGN:64:: bn5vec
      real(t_p) bn5vec(maxelst)
!PGI$ ATTRIBUTES ALIGN:64:: pgammavec
      real(t_p) pgammavec(maxelst)
!PGI$ ATTRIBUTES ALIGN:64:: dampvec,invdampvec
      real(t_p) dampvec(maxelst),invdampvec(maxelst)
!PGI$ ATTRIBUTES ALIGN:64:: dampvec1,expdampvec1
      real(t_p) dampvec1(maxelst),expdampvec1(maxelst)
!PGI$ ATTRIBUTES ALIGN:64:: tholevec
      real(t_p) tholevec(maxelst)
!PGI$ ATTRIBUTES ALIGN:64:: dscalevec,pscalevec
      real(t_p) dscalevec(maxelst),pscalevec(maxelst)
!PGI$ ATTRIBUTES ALIGN:64:: fimpvecx,fimpvecy,fimpvecz
      real(t_p),dimension(maxelst)::fimpvecx,fimpvecy,fimpvecz
!PGI$ ATTRIBUTES ALIGN:64:: fimdvecx,fimdvecy,fimdvecz
      real(t_p),dimension(maxelst)::fimdvecx,fimdvecy,fimdvecz
!PGI$ ATTRIBUTES ALIGN:64:: fidvecx,fidvecy,fidvecz
      real(t_p),dimension(maxelst)::fidvecx,fidvecy,fidvecz
!PGI$ ATTRIBUTES ALIGN:64:: fipvecx,fipvecy,fipvecz
      real(t_p),dimension(maxelst)::fipvecx,fipvecy,fipvecz
!PGI$ ATTRIBUTES ALIGN:64:: fimvecx,fimvecy,fimvecz
      real(t_p),dimension(maxelst)::fimvecx,fimvecy,fimvecz
!PGI$ ATTRIBUTES ALIGN:64:: fkmdvecx,fkmdvecy,fkmdvecz
      real(t_p),dimension(maxelst)::fkmdvecx,fkmdvecy,fkmdvecz
!PGI$ ATTRIBUTES ALIGN:64:: fkmpvecx,fkmpvecy,fkmpvecz
      real(t_p),dimension(maxelst)::fkmpvecx,fkmpvecy,fkmpvecz
!PGI$ ATTRIBUTES ALIGN:64:: fkdvecx,fkdvecy,fkdvecz
      real(t_p),dimension(maxelst)::fkdvecx,fkdvecy,fkdvecz
!PGI$ ATTRIBUTES ALIGN:64:: fkpvecx,fkpvecy,fkpvecz
      real(t_p),dimension(maxelst)::fkpvecx,fkpvecy,fkpvecz
!PGI$ ATTRIBUTES ALIGN:64:: fkmvecx,fkmvecy,fkmvecz
      real(t_p),dimension(maxelst)::fkmvecx,fkmvecy,fkmvecz
!PGI$ ATTRIBUTES ALIGN:64:: qirvecx, qirvecy,qirvecz
      real(t_p),dimension(maxelst):: qirvecx,qirvecy,qirvecz
!PGI$ ATTRIBUTES ALIGN:64:: qkrvecx, qkrvecy,qkrvecz
      real(t_p),dimension(maxelst):: qkrvecx,qkrvecy,qkrvecz
!PGI$ ATTRIBUTES ALIGN:64:: ckvec
      real(t_p)  ckvec(maxelst)
!PGI$ ATTRIBUTES ALIGN:64:: dkvecx,dkvecy,dkvecz
      real(t_p),dimension(maxelst):: dkvecx,dkvecy,dkvecz
!PGI$ ATTRIBUTES ALIGN:64:: dirvec,dkrvec
      real(t_p),dimension(maxelst):: dirvec,dkrvec
!PGI$ ATTRIBUTES ALIGN:64:: divec,qivec
      real(t_p)  divec(3),qivec(9)
!PGI$ ATTRIBUTES ALIGN:64:: qkvec1,qkvec2,qkvec3 
      real(t_p),dimension(maxelst):: qkvec1,  qkvec2,  qkvec3 
!PGI$ ATTRIBUTES ALIGN:64:: qkvec4,qkvec5,qkvec6 
      real(t_p),dimension(maxelst):: qkvec4,  qkvec5,  qkvec6 
!PGI$ ATTRIBUTES ALIGN:64:: qkvec7,qkvec8,qkvec9 
      real(t_p),dimension(maxelst):: qkvec7,  qkvec8,  qkvec9 
!PGI$ ATTRIBUTES ALIGN:64:: qirrvec,qkrrvec
      real(t_p),dimension(maxelst):: qirrvec,qkrrvec
!PGI$ ATTRIBUTES ALIGN:64:: dsc3vec,dsc5vec,dsc7vec
      real(t_p),dimension(maxelst):: dsc3vec,dsc5vec,dsc7vec
!PGI$ ATTRIBUTES ALIGN:64:: psc3vec,psc5vec,psc7vec
      real(t_p),dimension(maxelst):: psc3vec,psc5vec,psc7vec
!PGI$ ATTRIBUTES ALIGN:64:: drr3vec,drr5vec,drr7vec
      real(t_p),dimension(maxelst):: drr3vec,drr5vec,drr7vec
!PGI$ ATTRIBUTES ALIGN:64:: prr3vec,prr5vec,prr7vec
      real(t_p),dimension(maxelst):: prr3vec,prr5vec,prr7vec
!PGI$ ATTRIBUTES ALIGN:64::xpsvec,ypsvec,zpsvec
      real(t_p),dimension(maxelst)::xpsvec,ypsvec,zpsvec
      end
