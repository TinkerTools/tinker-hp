#ifndef PAIR_POLAR_INC
#define PAIR_POLAR_INC

#include "tinker_types.h"
#include "tinker_cudart.h"
#include "switch_respa.f.inc"

      M_subroutine
     &               epolar1_couple(dpui,ip,dpuk,kp,r2,pos,
     &                  aewald,alsq2,alsq2n,pgamma,damp,f,
     &                  dscale,pscale,uscale,
     &                  e,frc,frc_r,trqi,trqk,do_correct)
!$acc routine
      use tinheader ,only: ti_p
      use tinTypes  ,only: rpole_elt,real3,real6,mdyn3_r
#ifdef TINKER_CUF
      use utilcu    ,only: f_erfc
#  if defined(SINGLE)||defined(MIXED)
      use utilcu    ,only: f_sqrt,f_exp
#  endif
#endif
      implicit none

      real(t_p)  ,intent(in):: aewald,alsq2,alsq2n,pgamma,damp,f
      real(t_p)  ,intent(in):: dscale,pscale,uscale
      real(t_p)  ,intent(in):: r2
      type(rpole_elt),intent(in)::ip,kp
      type(real3),intent(in):: pos
      type(real6),intent(in):: dpuk,dpui
      logical    ,intent(in):: do_correct
      real(t_p)  ,intent(inout):: e
      type(real3),intent(inout):: frc
      type(real3),intent(inout):: trqi,trqk
      type(mdyn3_r),intent(inout)::frc_r

      real(t_p) exp2a,ralpha
      real(t_p) one,two,half
      real(t_p) ck,dkx,dky,dkz
      real(t_p) r,invr2,invr
      real(t_p) bn(0:4)
      real(t_p) rr3,rr5,rr7,rr9
      real(t_p) damp1,invdamp
      real(t_p) expdamp,expdamp1
      real(t_p) da
      real(t_p) ukx,uky,ukz,ukpx,ukpy,ukpz
      real(t_p) sc3,psc3,dsc3,usc3,psr3,dsr3,usr3
      real(t_p) sc5,psc5,dsc5,usc5,psr5,dsr5,usr5
      real(t_p) sc7,psc7,dsc7,usc7,psr7,dsr7,usr7
      real(t_p) rc3,rc3x,rc3y,rc3z
      real(t_p) prc3x,prc3y,prc3z
      real(t_p) drc3x,drc3y,drc3z
      real(t_p) urc3x,urc3y,urc3z
      real(t_p) rc5,rc5x,rc5y,rc5z
      real(t_p) prc5x,prc5y,prc5z
      real(t_p) drc5x,drc5y,drc5z
      real(t_p) urc5x,urc5y,urc5z
      real(t_p) rc7,rc7x,rc7y,rc7z
      real(t_p) prc7x,prc7y,prc7z
      real(t_p) drc7x,drc7y,drc7z
      real(t_p) urc7x,urc7y,urc7z
      real(t_p) dri,drk
      real(t_p) qrix,qriy,qriz
      real(t_p) qrkx,qrky,qrkz
      real(t_p) qrri,qrrk
      real(t_p) uri,urip,urk,urkp,duik,quik
      real(t_p) qrimodx,qrimody,qrimodz
      real(t_p) qrkmodx,qrkmody,qrkmodz
      real(t_p) term1,term2,term3
      real(t_p) dterm1,dterm2
      real(t_p) dterm3x,dterm3y,dterm3z
      real(t_p) dterm4x,dterm4y,dterm4z
      real(t_p) dterm5x,dterm5y,dterm5z
      real(t_p) dterm6x,dterm6y,dterm6z
      real(t_p) dterm7x,dterm7y,dterm7z
      real(t_p) tmp1x,tmp1y,tmp1z
      real(t_p) tisx,tisy,tisz,ticx,ticy,ticz
      real(t_p) tkcx,tkcy,tkcz,tksx,tksy,tksz
      real(t_p) ti5x,ti5y,ti5z,tk5x,tk5y,tk5z
      real(t_p) turi5,turi7,turk5,turk7
      real(t_p) depx,depy,depz
      type(real6) :: dufli,duflk
      type(real3) :: ufli,uflk
      parameter(half=0.5)
      parameter(one=1.0, two=2.0)
c
c     t reciprocal distance terms for this interaction
c
      invr2    = r2**(-1)
      r        = f_sqrt(r2)
      invr     = f_sqrt(invr2)
c
c     Calculate the real space Ewald error function terms
c
      if (do_correct) then
      bn(0)    = 0.0_ti_p
      bn(1)    = 0.0_ti_p
      bn(2)    = 0.0_ti_p
      bn(3)    = 0.0_ti_p
      bn(4)    = 0.0_ti_p
      else
      ralpha   = aewald * r
      !call erfcore_inl(ralpha, bn(0),1)
      exp2a    = f_exp( - ralpha**2)
      bn(0)    = f_erfc(ralpha)

      bn(0)    = bn(0) * invr
      bn(1)    = ( 1.0_ti_p*bn(0) + alsq2    *alsq2n*exp2a ) * invr2
      bn(2)    = ( 3.0_ti_p*bn(1) + alsq2**2 *alsq2n*exp2a ) * invr2
      bn(3)    = ( 5.0_ti_p*bn(2) + alsq2**2*alsq2    *alsq2n*exp2a )
     &         * invr2
      bn(4)    = ( 7.0_ti_p*bn(3) + alsq2**2*alsq2**2 *alsq2n*exp2a )
     &         * invr2

      bn(0)    = f * bn(0)
      bn(1)    = f * bn(1)
      bn(2)    = f * bn(2)
      bn(3)    = f * bn(3)
      bn(4)    = f * bn(4)
      end if

      rr3      = f        * invr * invr2
      rr5      = 3.0_ti_p * f*invr*invr2  * invr2
      rr7      = 5.0_ti_p * rr5  * invr2
      rr9      = 7.0_ti_p * rr7  * invr2
c
c    Apply Thole polarization damping to scale factors
c
      if (damp.ne.0.0_ti_p) then
         invdamp  = damp**(-one)
         damp1    = - pgamma * (r*invdamp)**2 * (r*invdamp)
         expdamp1 = f_exp(damp1)
         da       = damp1 * expdamp1
c
c     termediates involving Thole damping and scale factors
c
         sc3      = 1.0_ti_p - expdamp1
         sc5      = 1.0_ti_p - (1.0_ti_p - damp1)*expdamp1
         sc7      = 1.0_ti_p - (1.0_ti_p - damp1 + 0.6_ti_p*damp1**2)
     &                           *expdamp1
         rc3x     = - 3.0_ti_p *da *pos%x *invr2
         rc3y     = - 3.0_ti_p *da *pos%y *invr2
         rc3z     = - 3.0_ti_p *da *pos%z *invr2
         rc5x     = - damp1 * rc3x
         rc5y     = - damp1 * rc3y
         rc5z     = - damp1 * rc3z
         rc7x     = - (0.2_ti_p + 0.6_ti_p*damp1)*rc5x
         rc7y     = - (0.2_ti_p + 0.6_ti_p*damp1)*rc5y
         rc7z     = - (0.2_ti_p + 0.6_ti_p*damp1)*rc5z
      else
         sc3 = 1; sc5=1; sc7=1;
         rc3x = 0; rc3y=0; rc3z=0;
         rc5x = 0; rc5y=0; rc5z=0;
         rc7x = 0; rc7y=0; rc7z=0;
      end if

      if (do_correct) then
         psc3   = sc3 * pscale
         dsc3   = sc3 * dscale
         usc3   = sc3 * uscale
         psc5   = sc5 * pscale
         dsc5   = sc5 * dscale
         usc5   = sc5 * uscale
         psc7   = sc7 * pscale
         dsc7   = sc7 * dscale

         prc3x  = - rc3x * pscale
         drc3x  = - rc3x * dscale
         urc3x  = - rc3x * uscale
         prc3y  = - rc3y * pscale
         drc3y  = - rc3y * dscale
         urc3y  = - rc3y * uscale
         prc3z  = - rc3z * pscale
         drc3z  = - rc3z * dscale
         urc3z  = - rc3z * uscale

         prc5x  = - rc5x * pscale
         drc5x  = - rc5x * dscale
         urc5x  = - rc5x * uscale
         prc5y  = - rc5y * pscale
         drc5y  = - rc5y * dscale
         urc5y  = - rc5y * uscale
         prc5z  = - rc5z * pscale
         drc5z  = - rc5z * dscale
         urc5z  = - rc5z * uscale

         prc7x  = - rc7x * pscale
         drc7x  = - rc7x * dscale
         prc7y  = - rc7y * pscale
         drc7y  = - rc7y * dscale
         prc7z  = - rc7z * pscale
         drc7z  = - rc7z * dscale
      else
         psc3   = 1.0_ti_p - sc3 * pscale
         dsc3   = 1.0_ti_p - sc3 * dscale
         usc3   = 1.0_ti_p - sc3 * uscale
         psc5   = 1.0_ti_p - sc5 * pscale
         dsc5   = 1.0_ti_p - sc5 * dscale
         usc5   = 1.0_ti_p - sc5 * uscale
         psc7   = 1.0_ti_p - sc7 * pscale
         dsc7   = 1.0_ti_p - sc7 * dscale

         prc3x  = rc3x * pscale
         drc3x  = rc3x * dscale
         urc3x  = rc3x * uscale
         prc3y  = rc3y * pscale
         drc3y  = rc3y * dscale
         urc3y  = rc3y * uscale
         prc3z  = rc3z * pscale
         drc3z  = rc3z * dscale
         urc3z  = rc3z * uscale

         prc5x  = rc5x * pscale
         drc5x  = rc5x * dscale
         urc5x  = rc5x * uscale
         prc5y  = rc5y * pscale
         drc5y  = rc5y * dscale
         urc5y  = rc5y * uscale
         prc5z  = rc5z * pscale
         drc5z  = rc5z * dscale
         urc5z  = rc5z * uscale

         prc7x  = rc7x * pscale
         drc7x  = rc7x * dscale
         prc7y  = rc7y * pscale
         drc7y  = rc7y * dscale
         prc7z  = rc7z * pscale
         drc7z  = rc7z * dscale
      end if

      psr3     = bn(1) - psc3 * rr3
      dsr3     = bn(1) - dsc3 * rr3
      usr3     = bn(1) - usc3 * rr3
      psr5     = bn(2) - psc5 * rr5
      dsr5     = bn(2) - dsc5 * rr5
      usr5     = bn(2) - usc5 * rr5
      psr7     = bn(3) - psc7 * rr7
      dsr7     = bn(3) - dsc7 * rr7
c
c     termediates involving moments and distance separation
c
      dri      =    ip%dx*pos%x +   ip%dy*pos%y +  ip%dz*pos%z
      drk      =    kp%dx*pos%x +   kp%dy*pos%y +  kp%dz*pos%z

      qrix     =   ip%qxx*pos%x +  ip%qxy*pos%y +  ip%qxz*pos%z
      qriy     =   ip%qxy*pos%x +  ip%qyy*pos%y +  ip%qyz*pos%z
      qriz     =   ip%qxz*pos%x +  ip%qyz*pos%y +  ip%qzz*pos%z
      qrkx     =   kp%qxx*pos%x +  kp%qxy*pos%y +  kp%qxz*pos%z
      qrky     =   kp%qxy*pos%x +  kp%qyy*pos%y +  kp%qyz*pos%z
      qrkz     =   kp%qxz*pos%x +  kp%qyz*pos%y +  kp%qzz*pos%z
      qrri     =     qrix*pos%x +    qriy*pos%y +    qriz*pos%z
      qrrk     =     qrkx*pos%x +    qrky*pos%y +    qrkz*pos%z

      uri      =   dpui%x*pos%x +  dpui%y*pos%y +  dpui%z*pos%z
      urk      =   dpuk%x*pos%x +  dpuk%y*pos%y +  dpuk%z*pos%z
      urip     =  dpui%xx*pos%x + dpui%yy*pos%y + dpui%zz*pos%z
      urkp     =  dpuk%xx*pos%x + dpuk%yy*pos%y + dpuk%zz*pos%z

      duik     =   ip%dx*dpuk%x  +  ip%dy*dpuk%y + ip%dz*dpuk%z
     &          +  kp%dx*dpui%x  +  kp%dy*dpui%y + kp%dz*dpui%z
      quik     =    qrix*dpuk%x  +   qriy*dpuk%y +  qriz*dpuk%z
     &           -  qrkx*dpui%x  -   qrky*dpui%y -  qrkz*dpui%z
c
c     Calculate intermediate terms for polarization interaction
c
      term1    =  kp%c*uri  - ip%c*urk + duik
      term2    =   two*quik -  uri*drk - dri*urk
      term3    =   uri*qrrk -  urk*qrri
c
c     compute the energy contribution for this interaction
c
      e        = WRITE_C(e +) term1*psr3 +   term2*psr5 +  term3*psr7

      qrimodx  =                  qriy*pos%y  +  qriz*pos%z
      qrimody  =  qrix*pos%x   +                 qriz*pos%z
      qrimodz  =  qrix*pos%x   +  qriy*pos%y
      qrkmodx  =                  qrky*pos%y  +  qrkz*pos%z
      qrkmody  =  qrkx*pos%x   +                 qrkz*pos%z
      qrkmodz  =  qrkx*pos%x   +  qrky*pos%y
c
c     t the dEd/dR terms used for direct polarization force
c
      dterm1   =  bn(2)   -   dsc3*rr5
      dterm2   =  bn(3)   -   dsc5*rr7
c
c     Straight terms ( xx, yy ,zz )
c
      dterm3x  = - dsr3  + dterm1*pos%x**2 - rr3*pos%x*drc3x
      dterm3y  = - dsr3  + dterm1*pos%y**2 - rr3*pos%y*drc3y
      dterm3z  = - dsr3  + dterm1*pos%z**2 - rr3*pos%z*drc3z

      dterm4x  =  rr3*drc3x  - (dterm1 + dsr5)*pos%x
      dterm4y  =  rr3*drc3y  - (dterm1 + dsr5)*pos%y
      dterm4z  =  rr3*drc3z  - (dterm1 + dsr5)*pos%z

      dterm5x  = - dsr5      + dterm2*pos%x**2 - rr5*pos%x*drc5x
      dterm5y  = - dsr5      + dterm2*pos%y**2 - rr5*pos%y*drc5y
      dterm5z  = - dsr5      + dterm2*pos%z**2 - rr5*pos%z*drc5z

      dterm6x  =  (bn(4) - dsc7*rr9)*pos%x**2  - bn(3) - rr7*pos%x*drc7x
      dterm6y  =  (bn(4) - dsc7*rr9)*pos%y**2  - bn(3) - rr7*pos%y*drc7y
      dterm6z  =  (bn(4) - dsc7*rr9)*pos%z**2  - bn(3) - rr7*pos%z*drc7z

      dterm7x  = rr5*drc5x  - two*bn(3)*pos%x
     &           + ( dsc5   + 1.5_ti_p*dsc7 )*rr7*pos%x
      dterm7y  = rr5*drc5y  - two*bn(3)*pos%y
     &           + ( dsc5   + 1.5_ti_p*dsc7)*rr7*pos%y
      dterm7z  = rr5*drc5z  - two*bn(3)*pos%z
     &           + ( dsc5   + 1.5_ti_p*dsc7 )*rr7*pos%z
c
c     Straight terms ( xx, yy ,zz )
c
      tisx     = ip%c*dterm3x      + ip%dx*dterm4x + dri*dterm5x
     &          + two*dsr5*ip%qxx  + qrimodx*dsc7*rr7
     &          + two*qrix*dterm7x + qrri*dterm6x
      tisy     = ip%c*dterm3y      + ip%dy*dterm4y + dri*dterm5y
     &          + two*dsr5*ip%qyy  + qrimody*dsc7*rr7
     &          + two*qriy*dterm7y + qrri*dterm6y
      tisz     = ip%c*dterm3z      + ip%dz*dterm4z + dri*dterm5z
     &          + two*dsr5*ip%qzz  + qrimodz*dsc7*rr7
     &          + two*qriz*dterm7z +  qrri*dterm6z

      tksx     = kp%c*dterm3x      - kp%dx*dterm4x - drk*dterm5x
     &          + two*dsr5*kp%qxx  + qrkmodx*dsc7*rr7
     &          + two*qrkx*dterm7x + qrrk*dterm6x
      tksy     = kp%c*dterm3y      - kp%dy*dterm4y - drk*dterm5y
     &          + two*dsr5*kp%qyy  + qrkmody*dsc7*rr7
     &          + two*qrky*dterm7y + qrrk*dterm6y
      tksz     = kp%c*dterm3z      - kp%dz*dterm4z  - drk*dterm5z
     &          + two*dsr5*kp%qzz  + qrkmodz*dsc7*rr7
     &          + two*qrkz*dterm7z + qrrk*dterm6z
c
c     Crossed terms ( xy = yx , xz = zx ,yz = zy )
c
      tmp1x    = pos%x * pos%y
      tmp1y    = pos%x * pos%z
      tmp1z    = pos%y * pos%z

      dterm3x  =   dterm1*tmp1x  -  rr3*pos%y*drc3x
      dterm3y  =   dterm1*tmp1y  -  rr3*pos%z*drc3x
      dterm3z  =   dterm1*tmp1z  -  rr3*pos%z*drc3y
      dterm4x  = - dterm1*pos%x  +  rr3*drc3x
      dterm4y  = - dterm1*pos%x  +  rr3*drc3x
      dterm4z  = - dterm1*pos%y  +  rr3*drc3y

      dterm5x  =   dterm2*tmp1x  -  rr5*pos%y*drc5x
      dterm5y  =   dterm2*tmp1y  -  rr5*pos%z*drc5x
      dterm5z  =   dterm2*tmp1z  -  rr5*pos%z*drc5y
      dterm6x  =  (bn(4) - dsc7*rr9)*tmp1x - rr7*pos%y*drc7x
      dterm6y  =  (bn(4) - dsc7*rr9)*tmp1y - rr7*pos%z*drc7x
      dterm6z  =  (bn(4) - dsc7*rr9)*tmp1z - rr7*pos%z*drc7y
      dterm7x  = - dterm2*pos%x  +  rr5*drc5x
      dterm7y  = - dterm2*pos%x  +  rr5*drc5x
      dterm7z  = - dterm2*pos%y  +  rr5*drc5y
c
c     Crossed terms ( xy = yx , xz = zx ,yz = zy )
c

      ticx     = ip%c*dterm3x + ip%dy*dterm4x + dri*dterm5x
     &         + qrri*dterm6x - dsr5*ip%dx*pos%y
     &         +  two*(dsr5*ip%qxy - dsr7*pos%y*qrix + qriy*dterm7x)

      ticy     = ip%c*dterm3y + ip%dz*dterm4y + dri*dterm5y
     &         + qrri*dterm6y - dsr5*ip%dx*pos%z
     &         +  two*(dsr5*ip%qxz - dsr7*pos%z*qrix + qriz*dterm7y)

      ticz     = ip%c*dterm3z + ip%dz*dterm4z + dri*dterm5z
     &         + qrri*dterm6z - dsr5*ip%dy*pos%z
     &         +  two*(dsr5*ip%qyz - dsr7*pos%z*qriy + qriz*dterm7z)

      tkcx     = kp%c*dterm3x - kp%dy*dterm4x - drk*dterm5x
     &         + qrrk*dterm6x + dsr5*kp%dx*pos%y
     &         +  two*(dsr5*kp%qxy - dsr7*pos%y*qrkx + qrky*dterm7x)

      tkcy     = kp%c*dterm3y - kp%dz*dterm4y - drk*dterm5y
     &         + qrrk*dterm6y + dsr5*kp%dx*pos%z
     &         +  two*(dsr5*kp%qxz - dsr7*pos%z*qrkx + qrkz*dterm7y)

      tkcz     = kp%c*dterm3z- kp%dz*dterm4z - drk*dterm5z
     &         + qrrk*dterm6z + dsr5 *kp%dy*pos%z
     &         +  two*(dsr5*kp%qyz - dsr7*pos%z*qrky + qrkz*dterm7z)
c
c      Construct matrixes for dot_product
c      do Dot product
c
      frc%x     =  tisx*dpuk%xx + ticx*dpuk%yy + ticy*dpuk%zz
     &          -  tksx*dpui%xx - tkcx*dpui%yy - tkcy*dpui%zz
      frc%y     =  ticx*dpuk%xx + tisy*dpuk%yy + ticz*dpuk%zz
     &          -  tkcx*dpui%xx - tksy*dpui%yy - tkcz*dpui%zz
      frc%z     =  ticy*dpuk%xx + ticz*dpuk%yy + tisz*dpuk%zz
     &          -  tkcy*dpui%xx - tkcz*dpui%yy - tksz*dpui%zz
c
c     t the dEp/dR terms used for direct polarization force
c

      dterm1   =  bn(2) - psc3 * rr5
      dterm2   =  bn(3) - psc5 * rr7
c
c     Straight terms ( xx, yy ,zz )
c
      dterm3x  = - psr3 + dterm1*pos%x**2 - rr3*pos%x*prc3x
      dterm3y  = - psr3 + dterm1*pos%y**2 - rr3*pos%y*prc3y
      dterm3z  = - psr3 + dterm1*pos%z**2 - rr3*pos%z*prc3z

      dterm4x  = - psr5*pos%x - dterm1*pos%x + rr3*prc3x
      dterm4y  = - psr5*pos%y - dterm1*pos%y + rr3*prc3y
      dterm4z  = - psr5*pos%z - dterm1*pos%z + rr3*prc3z

      dterm5x  = - psr5 + dterm2*pos%x**2 - rr5*pos%x*prc5x
      dterm5y  = - psr5 + dterm2*pos%y**2 - rr5*pos%y*prc5y
      dterm5z  = - psr5 + dterm2*pos%z**2 - rr5*pos%z*prc5z

      dterm6x  =  (bn(4) - psc7*rr9)*pos%x**2 - bn(3) - rr7*pos%x*prc7x
      dterm6y  =  (bn(4) - psc7*rr9)*pos%y**2 - bn(3) - rr7*pos%y*prc7y
      dterm6z  =  (bn(4) - psc7*rr9)*pos%z**2 - bn(3) - rr7*pos%z*prc7z

      dterm7x  =  rr5*prc5x - two*bn(3)*pos%x
     &         +  (psc5 + 1.5_ti_p*psc7)*rr7*pos%x
      dterm7y  =  rr5*prc5y - two*bn(3)*pos%y
     &         +  (psc5 + 1.5_ti_p*psc7)*rr7*pos%y
      dterm7z  =  rr5*prc5z - two*bn(3)*pos%z
     &         +  (psc5 + 1.5_ti_p*psc7)*rr7*pos%z
c
c     Straight terms ( xx, yy ,zz )
c
      tisx = ip%c*dterm3x + ip%dx*dterm4x + dri*dterm5x
     &     + qrri*dterm6x + qrimodx*psc7*rr7
     &     + two*(psr5*ip%qxx + qrix*dterm7x)

      tisy = ip%c*dterm3y + ip%dy*dterm4y + dri*dterm5y
     &     + qrri*dterm6y + qrimody*psc7*rr7
     &     + two*(psr5*ip%qyy + qriy*dterm7y)

      tisz = ip%c*dterm3z + ip%dz*dterm4z + dri*dterm5z
     &     + qrri*dterm6z + qrimodz*psc7*rr7
     &     + two*(psr5*ip%qzz + qriz*dterm7z)

      tksx = kp%c*dterm3x - kp%dx*dterm4x - drk*dterm5x
     &     + qrrk*dterm6x + qrkmodx*psc7*rr7
     &     + two*(psr5*kp%qxx + qrkx*dterm7x)

      tksy = kp%c*dterm3y - kp%dy*dterm4y - drk*dterm5y
     &     + qrrk*dterm6y + qrkmody*psc7*rr7
     &     + two*(psr5*kp%qyy + qrky*dterm7y)

      tksz = kp%c*dterm3z - kp%dz*dterm4z - drk*dterm5z
     &     + qrrk*dterm6z + qrkmodz*psc7*rr7
     &     + two*(psr5*kp%qzz + qrkz*dterm7z)
c
c     Crossed terms ( xy = yx , xz = zx ,yz = zy )
c
      tmp1x = pos%x * pos%y
      tmp1y = pos%x * pos%z
      tmp1z = pos%y * pos%z

      dterm3x =  dterm1*tmp1x - rr3*pos%y*prc3x
      dterm3y =  dterm1*tmp1y - rr3*pos%z*prc3x
      dterm3z =  dterm1*tmp1z - rr3*pos%z*prc3y

      dterm4x =  rr3*prc3x - dterm1*pos%x
      dterm4y =  rr3*prc3x - dterm1*pos%x
      dterm4z =  rr3*prc3y - dterm1*pos%y
      dterm5x =  dterm2*tmp1x - rr5*pos%y*prc5x
      dterm5y =  dterm2*tmp1y - rr5*pos%z*prc5x
      dterm5z =  dterm2*tmp1z - rr5*pos%z*prc5y

      dterm6x =  (bn(4) - psc7*rr9)*tmp1x - rr7*pos%y*prc7x
      dterm6y =  (bn(4) - psc7*rr9)*tmp1y - rr7*pos%z*prc7x
      dterm6z =  (bn(4) - psc7*rr9)*tmp1z - rr7*pos%z*prc7y
      dterm7x =  rr5*prc5x - dterm2*pos%x
      dterm7y =  rr5*prc5x - dterm2*pos%x
      dterm7z =  rr5*prc5y - dterm2*pos%y

c
c     Crossed terms ( xy = yx , xz = zx ,yz = zy )
c
      ticx = ip%c*dterm3x + ip%dy*dterm4x + dri*dterm5x
     &     + qrri*dterm6x - psr5*ip%dx*pos%y
     &     + two*(psr5*ip%qxy - psr7*pos%y*qrix + qriy*dterm7x)

      ticy = ip%c*dterm3y + ip%dz*dterm4y + dri*dterm5y
     &     + qrri*dterm6y - psr5*ip%dx*pos%z
     &     + two*(psr5*ip%qxz - psr7*pos%z*qrix + qriz*dterm7y)

      ticz = ip%c*dterm3z + ip%dz*dterm4z + dri*dterm5z
     &     + qrri*dterm6z - psr5*ip%dy*pos%z
     &     +  two*(psr5*ip%qyz - psr7*pos%z*qriy + qriz*dterm7z)

      tkcx = kp%c*dterm3x - kp%dy*dterm4x - drk*dterm5x
     &     + qrrk*dterm6x + psr5*kp%dx*pos%y
     &     +  two*(psr5*kp%qxy - psr7*pos%y*qrkx + qrky*dterm7x)

      tkcy = kp%c*dterm3y - kp%dz*dterm4y - drk*dterm5y
     &     + qrrk*dterm6y + psr5*kp%dx*pos%z
     &     + two*(psr5*kp%qxz - psr7*pos%z*qrkx + qrkz*dterm7y)

      tkcz = kp%c*dterm3z - kp%dz*dterm4z - drk*dterm5z
     &     + qrrk*dterm6z + psr5*kp%dy*pos%z
     &     + two*(psr5*kp%qyz - psr7*pos%z*qrky + qrkz*dterm7z)
c
c     Construct matrixes for dot_product
c     Do dot product
c
      depx = tisx*dpuk%x + ticx*dpuk%y + ticy*dpuk%z
     &     - tksx*dpui%x - tkcx*dpui%y - tkcy*dpui%z
      depy = ticx*dpuk%x + tisy*dpuk%y + ticz*dpuk%z
     &     - tkcx*dpui%x - tksy*dpui%y - tkcz*dpui%z
      depz = ticy*dpuk%x + ticz*dpuk%y + tisz*dpuk%z
     &     - tkcy*dpui%x - tkcz*dpui%y - tksz*dpui%z

      frc%x = frc%x + depx
      frc%y = frc%y + depy
      frc%z = frc%z + depz

c
c     t the dtau/dr terms used for mutual polarization force
c
      dterm1 = bn(2) - usc3 * rr5
      dterm2 = bn(3) - usc5 * rr7
c
c     Straight terms ( xx, yy ,zz )
c
      dterm3x = usr5 + dterm1
      dterm3y = usr5 + dterm1
      dterm3z = usr5 + dterm1
      dterm4x = rr3  * uscale
      dterm4y = rr3  * uscale
      dterm4z = rr3  * uscale

      dterm5x = - pos%x*dterm3x + rc3x*dterm4x
      dterm5y = - pos%y*dterm3y + rc3y*dterm4y
      dterm5z = - pos%z*dterm3z + rc3z*dterm4z

      dterm6x = - usr5 + pos%x**2*dterm2 - rr5*pos%x*urc5x
      dterm6y = - usr5 + pos%y**2*dterm2 - rr5*pos%y*urc5y
      dterm6z = - usr5 + pos%z**2*dterm2 - rr5*pos%z*urc5z

      tisx =  dpui%x*dterm5x + uri*dterm6x
      tisy =  dpui%y*dterm5y + uri*dterm6y
      tisz =  dpui%z*dterm5z + uri*dterm6z
      tksx =  dpuk%x*dterm5x + urk*dterm6x
      tksy =  dpuk%y*dterm5y + urk*dterm6y
      tksz =  dpuk%z*dterm5z + urk*dterm6z
c
c     Crossed terms ( xy = yx , xz = zx ,yz = zy )
c
      dterm4x = - usr5*pos%y
      dterm4y = - usr5*pos%z
      dterm4z = - usr5*pos%z
      dterm5x = - pos%x*dterm1 + rr3*urc3x
      dterm5y = - pos%x*dterm1 + rr3*urc3x
      dterm5z = - pos%y*dterm1 + rr3*urc3y
      dterm6x =   pos%x*dterm2*pos%y - rr5*pos%y*urc5x
      dterm6y =   pos%x*dterm2*pos%z - rr5*pos%z*urc5x
      dterm6z =   pos%y*dterm2*pos%z - rr5*pos%z*urc5y

      ticx =  dpui%x*dterm4x + dpui%y*dterm5x + uri*dterm6x
      ticy =  dpui%x*dterm4y + dpui%z*dterm5y + uri*dterm6y
      ticz =  dpui%y*dterm4z + dpui%z*dterm5z + uri*dterm6z

      tkcx =  dpuk%x*dterm4x + dpuk%y*dterm5x + urk*dterm6x
      tkcy =  dpuk%x*dterm4y + dpuk%z*dterm5y + urk*dterm6y
      tkcz =  dpuk%y*dterm4z + dpuk%z*dterm5z + urk*dterm6z
c
c     Construct matrixes for dot_product
c     Dot product
c
      depx =  tisx*dpuk%xx + ticx*dpuk%yy + ticy*dpuk%zz
     &      + tksx*dpui%xx + tkcx*dpui%yy + tkcy*dpui%zz
      depy =  ticx*dpuk%xx + tisy*dpuk%yy + ticz*dpuk%zz
     &      + tkcx*dpui%xx + tksy*dpui%yy + tkcz*dpui%zz
      depz =  ticy*dpuk%xx + ticz*dpuk%yy + tisz*dpuk%zz
     &      + tkcy*dpui%xx + tkcz*dpui%yy + tksz*dpui%zz

      frc%x   = frc%x + depx
      frc%y   = frc%y + depy
      frc%z   = frc%z + depz
      frc_r%x = WRITE_C(frc_r%x +) tp2mdr(frc%x)
      frc_r%y = WRITE_C(frc_r%y +) tp2mdr(frc%y)
      frc_r%z = WRITE_C(frc_r%z +) tp2mdr(frc%z)
c
c     t the induced dipole field used for dipole torques
c
      turi5    = -psr5*urk - dsr5*urkp
      turk5    = -psr5*uri - dsr5*urip
      turi7    = - psr7*urk - dsr7*urkp
      turk7    = - psr7*uri - dsr7*urip
c
c     t induced dipole field gradient used for quadrupole torques
c
      ti5x     =  two*(psr5*dpuk%x + dsr5*dpuk%xx)
      ti5y     =  two*(psr5*dpuk%y + dsr5*dpuk%yy)
      ti5z     =  two*(psr5*dpuk%z + dsr5*dpuk%zz)
      tk5x     =  two*(psr5*dpui%x + dsr5*dpui%xx)
      tk5y     =  two*(psr5*dpui%y + dsr5*dpui%yy)
      tk5z     =  two*(psr5*dpui%z + dsr5*dpui%zz)
c
c     Torque is induced field and gradient cross permanent moments
c
      ufli%x   = psr3*dpuk%x + dsr3*dpuk%xx + pos%x*turi5
      ufli%y   = psr3*dpuk%y + dsr3*dpuk%yy + pos%y*turi5
      ufli%z   = psr3*dpuk%z + dsr3*dpuk%zz + pos%z*turi5
      uflk%x   = psr3*dpui%x + dsr3*dpui%xx + pos%x*turk5
      uflk%y   = psr3*dpui%y + dsr3*dpui%yy + pos%y*turk5
      uflk%z   = psr3*dpui%z + dsr3*dpui%zz + pos%z*turk5

      dufli%x  = pos%x*ti5x + pos%x**2*turi7
      dufli%z  = pos%y*ti5y + pos%y**2*turi7
      dufli%zz = pos%z*ti5z + pos%z**2*turi7
      dufli%y  = pos%x*ti5y + pos%y*ti5x + two*pos%x*pos%y*turi7
      dufli%xx = pos%x*ti5z + pos%z*ti5x + two*pos%x*pos%z*turi7
      dufli%yy = pos%y*ti5z + pos%z*ti5y + two*pos%y*pos%z*turi7

      duflk%x  = -(pos%x*tk5x + pos%x**2*turk7)
      duflk%z  = -(pos%y*tk5y + pos%y**2*turk7)
      duflk%zz = -(pos%z*tk5z + pos%z**2*turk7)
      duflk%y  = -(pos%x*tk5y + pos%y*tk5x + two*pos%x*pos%y*turk7)
      duflk%xx = -(pos%x*tk5z + pos%z*tk5x + two*pos%x*pos%z*turk7)
      duflk%yy = -(pos%y*tk5z + pos%z*tk5y + two*pos%y*pos%z*turk7)

      trqi%x = WRITE_C(trqi%x +) ip%dz*ufli%y - ip%dy*ufli%z
     &       + ip%qxz*dufli%y - ip%qxy*dufli%xx
     &       + (ip%qzz - ip%qyy)*dufli%yy
     &       + two*ip%qyz*( dufli%z - dufli%zz )
      trqi%y = WRITE_C(trqi%y +) ip%dx*ufli%z - ip%dz*ufli%x
     &       - ip%qyz*dufli%y
     &       + ip%qxy*dufli%yy + (ip%qxx - ip%qzz)*dufli%xx
     &       + two*ip%qxz*( dufli%zz - dufli%x )
      trqi%z = WRITE_C(trqi%z +) ip%dy*ufli%x - ip%dx*ufli%y
     &       + ip%qyz*dufli%xx
     &       - ip%qxz*dufli%yy + (ip%qyy - ip%qxx)*dufli%y
     &       + two*ip%qxy*( dufli%x - dufli%z )

      trqk%x = WRITE_C(trqk%x +) kp%dz*uflk%y - kp%dy*uflk%z
     &       + kp%qxz*duflk%y
     &       - kp%qxy*duflk%xx + (kp%qzz - kp%qyy)*duflk%yy
     &       + two*kp%qyz*( duflk%z - duflk%zz )
      trqk%y = WRITE_C(trqk%y +) kp%dx*uflk%z - kp%dz*uflk%x
     &       - kp%qyz*duflk%y
     &       + kp%qxy*duflk%yy + (kp%qxx - kp%qzz)*duflk%xx
     &       + two*kp%qxz*( duflk%zz - duflk%x )
      trqk%z = WRITE_C(trqk%z +) kp%dy*uflk%x - kp%dx*uflk%y
     &       + kp%qyz*duflk%xx
     &       - kp%qxz*duflk%yy + (kp%qyy - kp%qxx)*duflk%y
     &       + two*kp%qxy*( duflk%x - duflk%z )
      end subroutine

      M_subroutine
     &               epolar3_couple(dpui,ip,dpuk,kp,r2,pos,
     &                  aewald,alsq2,alsq2n,pgamma,damp,f,
     &                  off,shortheal,pscale,e,use_short,do_correct)
!$acc routine
      use tinheader ,only: ti_p
      use utilgpu   ,only: rpole_elt,real3,real6,real3_red
#ifdef TINKER_CUF
      use utilcu ,only: f_erfc
#  if defined(SINGLE)||defined(MIXED)
      use utilcu ,only: f_sqrt,f_exp
#  endif
#endif
      implicit none

      real(t_p)  ,intent(in):: aewald,alsq2,alsq2n,pgamma,damp,f
      real(t_p)  ,intent(in):: pscale,off,shortheal
      real(t_p)  ,intent(in):: r2
      type(rpole_elt),intent(in)::ip,kp
      type(real3),intent(in):: pos
      type(real6),intent(in):: dpuk,dpui
      logical    ,intent(in):: do_correct
      logical    ,intent(in):: use_short
#ifdef TINKER_CUF
      ener_rtyp  ,intent(inout):: e
#else
      real(t_p)  ,intent(inout):: e
#endif

      real(t_p) exp2a,ralpha
      real(t_p) one,two,half
      real(t_p) t1,t2,t3,t4,t5,t6
      real(t_p) ck,dkx,dky,dkz
      real(t_p) r,invr2,invr
      real(t_p) bn(0:4)
      real(t_p) rr3,rr5,rr7,rr9
      real(t_p) damp1,invdamp
      real(t_p) expdamp,expdamp1
      real(t_p) sc3,psc3,psr3
      real(t_p) sc5,psc5,psr5
      real(t_p) sc7,psc7,psr7
      real(t_p) rc3
      real(t_p) uri,urk,dri,drk
      real(t_p) duik,quik
      real(t_p) qrix,qriy,qriz
      real(t_p) qrkx,qrky,qrkz
      real(t_p) qrri,qrrk
      real(t_p) term1,term2,term3
      real(t_p) s,ds,e_
      parameter(half=0.5)
      parameter(one=1.0, two=2.0)

c
c     t reciprocal distance terms for this interaction
c
      invr2    = 1.0_ti_p/r2
      r        = f_sqrt(r2)
      invr     = f_sqrt(invr2)
c
c     Calculate the real space Ewald error function terms
c
c
      if (do_correct) then
      bn(0)    = 0.0_ti_p
      bn(1)    = 0.0_ti_p
      bn(2)    = 0.0_ti_p
      bn(3)    = 0.0_ti_p
      bn(4)    = 0.0_ti_p
      else
      ralpha   = aewald * r
      !call erfcore_inl(ralpha, bn(0),1)
      exp2a    = f_exp( - ralpha**2)
      bn(0)    = f_erfc(ralpha)

      bn(0)    = bn(0) * invr
      bn(1)    = ( 1.0_ti_p*bn(0) + alsq2    *alsq2n*exp2a ) * invr2
      bn(2)    = ( 3.0_ti_p*bn(1) + alsq2**2 *alsq2n*exp2a ) * invr2
      bn(3)    = ( 5.0_ti_p*bn(2) + alsq2**2*alsq2    *alsq2n*exp2a )
     &           * invr2
      bn(4)    = ( 7.0_ti_p*bn(3) + alsq2**2*alsq2**2 *alsq2n*exp2a )
     &           * invr2

      bn(0)    = f * bn(0)
      bn(1)    = f * bn(1)
      bn(2)    = f * bn(2)
      bn(3)    = f * bn(3)
      bn(4)    = f * bn(4)
      end if

      rr3      = f        * invr * invr2
      rr5      = 3.0_ti_p * f*invr*invr2  * invr2
      rr7      = 5.0_ti_p * rr5  * invr2
      rr9      = 7.0_ti_p * rr7  * invr2
c
c    Apply Thole polarization damping to scale factors
c
      if (damp.ne.0.0_ti_p) then
         invdamp  = damp**(-one)
         damp1    = - pgamma * (r*invdamp)**2 * (r*invdamp)
         expdamp1 = exp(damp1)
c
c     termediates involving Thole damping and scale factors
c
         sc3      = 1.0_ti_p - expdamp1
         sc5      = 1.0_ti_p - (1.0_ti_p - damp1)*expdamp1
         sc7      = 1.0_ti_p - (1.0_ti_p - damp1 + 0.6_ti_p*damp1**2)
     &                           *expdamp1
      else
         sc3 = 1; sc5=1; sc7=1;
      end if

      if (do_correct) then
         psr3  = bn(1) - sc3 * pscale * rr3
         psr5  = bn(2) - sc5 * pscale * rr5
         psr7  = bn(3) - sc7 * pscale * rr7
      else
         psr3  = bn(1) - (1.0_ti_p - sc3 * pscale) * rr3
         psr5  = bn(2) - (1.0_ti_p - sc5 * pscale) * rr5
         psr7  = bn(3) - (1.0_ti_p - sc7 * pscale) * rr7
      end if
c
c     termediates involving moments and distance separation
c
      dri      =    ip%dx*pos%x +   ip%dy*pos%y +  ip%dz*pos%z
      drk      =    kp%dx*pos%x +   kp%dy*pos%y +  kp%dz*pos%z

      qrix     =   ip%qxx*pos%x +  ip%qxy*pos%y +  ip%qxz*pos%z
      qriy     =   ip%qxy*pos%x +  ip%qyy*pos%y +  ip%qyz*pos%z
      qriz     =   ip%qxz*pos%x +  ip%qyz*pos%y +  ip%qzz*pos%z
      qrkx     =   kp%qxx*pos%x +  kp%qxy*pos%y +  kp%qxz*pos%z
      qrky     =   kp%qxy*pos%x +  kp%qyy*pos%y +  kp%qyz*pos%z
      qrkz     =   kp%qxz*pos%x +  kp%qyz*pos%y +  kp%qzz*pos%z
      qrri     =     qrix*pos%x +    qriy*pos%y +    qriz*pos%z
      qrrk     =     qrkx*pos%x +    qrky*pos%y +    qrkz*pos%z

      uri      =   dpui%x*pos%x +  dpui%y*pos%y +  dpui%z*pos%z
      urk      =   dpuk%x*pos%x +  dpuk%y*pos%y +  dpuk%z*pos%z

      duik     =    ip%dx*dpuk%x  +  ip%dy*dpuk%y  +  ip%dz*dpuk%z
     &          +  dpui%x*kp%dx   + dpui%y*kp%dy   + dpui%z*kp%dz
      quik     =     qrix*dpuk%x  +   qriy*dpuk%y  +   qriz*dpuk%z
     &          -  dpui%x*qrkx    - dpui%y*qrky    - dpui%z*qrkz
c
c     Calculate intermediate terms for polarization interaction
c
      term1    =        kp%c*uri  -      ip%c*urk  +      duik
      term2    =        two*quik  -       uri*drk  -   dri*urk
      term3    =        uri*qrrk  -       urk*qrri
c
c     compute the energy contribution for this interaction
c
      e_       = term1*psr3 +   term2*psr5 +  term3*psr7

      if (use_short) then
      call switch_respa_inl(r,off,shortheal,s,ds)
      e_       = s*e_
      end if

      e        = WRITE_C(e + tp2enr) (e_)

      end subroutine

      M_subroutine
     &               mpolar1_couple(dpui,ip,dpuk,kp,r2,pos
     &                  ,aewald,alsq2,alsq2n,pgamma,damp,f
     &                  ,r_cut,shortheal
     &                  ,mscale,dscale,pscale,uscale
     &                  ,e,frc,frc_r,trqi,trqk,do_correct,mode
     &     ,iglob,kglob)
!$acc routine
      use tinheader ,only: ti_p
      use tinTypes  ,only: rpole_elt,real3,real6,mdyn3_r
      use interfaces,only: m_normal,m_short,m_long
#ifdef TINKER_CUF
      use utilcu    ,only: f_erfc
#  if defined(SINGLE)||defined(MIXED)
      use utilcu    ,only: f_sqrt,f_exp
#  endif
#endif
      implicit none

      integer iglob,kglob
      real(t_p)  ,intent(in),value:: aewald,alsq2,alsq2n,pgamma,damp,f
     &           ,r_cut,shortheal
      real(t_p)  ,intent(in):: mscale,dscale,pscale,uscale
      real(t_p)  ,intent(in):: r2
      type(rpole_elt),intent(in)::ip,kp
      type(real3),intent(in):: pos
      type(real6),intent(in):: dpuk,dpui
      integer    ,intent(in):: mode
      logical    ,intent(in):: do_correct
      real(t_p)  ,intent(inout):: e
      type(real3),intent(inout):: frc
      type(real3),intent(inout):: trqi,trqk
      type(mdyn3_r),intent(inout)::frc_r

      real(t_p) exp2a,ralpha
      real(t_p) one,two,half
      real(t_p) ck,dkx,dky,dkz
      real(t_p) r,invr2,invr
      real(t_p) bn(0:5)
      real(t_p) rr1,rr3,rr5,rr7,rr9,rr11
      real(t_p) damp1,invdamp
      real(t_p) expdamp,expdamp1
      real(t_p) da
      real(t_p) ukx,uky,ukz,ukpx,ukpy,ukpz
      real(t_p) sc3,psc3,dsc3,usc3,psr3,dsr3,usr3
      real(t_p) sc5,psc5,dsc5,usc5,psr5,dsr5,usr5
      real(t_p) sc7,psc7,dsc7,usc7,psr7,dsr7,usr7
      real(t_p) rc3,rc3x,rc3y,rc3z
      real(t_p) prc3x,prc3y,prc3z
      real(t_p) drc3x,drc3y,drc3z
      real(t_p) urc3x,urc3y,urc3z
      real(t_p) rc5,rc5x,rc5y,rc5z
      real(t_p) prc5x,prc5y,prc5z
      real(t_p) drc5x,drc5y,drc5z
      real(t_p) urc5x,urc5y,urc5z
      real(t_p) rc7,rc7x,rc7y,rc7z
      real(t_p) prc7x,prc7y,prc7z
      real(t_p) drc7x,drc7y,drc7z
      real(t_p) urc7x,urc7y,urc7z
      real(t_p) dri,drk
      real(t_p) qrix,qriy,qriz
      real(t_p) qrkx,qrky,qrkz
      real(t_p) qrri,qrrk
      real(t_p) uri,urip,urk,urkp,duik,quik
      real(t_p) qrimodx,qrimody,qrimodz
      real(t_p) qrkmodx,qrkmody,qrkmodz
      real(t_p) term1,term2,term3
      real(t_p) dterm1,dterm2
      real(t_p) dterm3x,dterm3y,dterm3z
      real(t_p) dterm4x,dterm4y,dterm4z
      real(t_p) dterm5x,dterm5y,dterm5z
      real(t_p) dterm6x,dterm6y,dterm6z
      real(t_p) dterm7x,dterm7y,dterm7z
      type(real3) tmp1
      real(t_p) tisx,tisy,tisz,ticx,ticy,ticz
      real(t_p) tkcx,tkcy,tkcz,tksx,tksy,tksz
      real(t_p) ti5x,ti5y,ti5z,tk5x,tk5y,tk5z
      real(t_p) turi5,turi7,turk5,turk7
      real(t_p) depx,depy,depz
      type(real6) :: dufli,duflk
      type(real3) :: ufli,uflk
      parameter(half=0.5)
      parameter(one=1.0, two=2.0)
c
c     t reciprocal distance terms for this interaction
c
      invr2    = r2**(-1)
      r        = f_sqrt(r2)
      invr     = f_sqrt(invr2)
c
c     Calculate the real space Ewald error function terms
c
      if (do_correct) then
      bn(0)    = 0.0_ti_p
      bn(1)    = 0.0_ti_p
      bn(2)    = 0.0_ti_p
      bn(3)    = 0.0_ti_p
      bn(4)    = 0.0_ti_p
      bn(5)    = 0.0_ti_p
      else
      ralpha   = aewald * r
      !call erfcore_inl(ralpha, bn(0),1)
      exp2a    = f_exp( - ralpha**2)
      bn(0)    = f_erfc(ralpha)

      bn(0)    = bn(0) * invr
      bn(1)    = ( 1.0_ti_p*bn(0) + alsq2    *alsq2n*exp2a ) * invr2
      bn(2)    = ( 3.0_ti_p*bn(1) + alsq2**2 *alsq2n*exp2a ) * invr2
      bn(3)    = ( 5.0_ti_p*bn(2) + alsq2**2*alsq2 *alsq2n*exp2a )
     &         * invr2
      tmp1%x   = alsq2**2*alsq2**2
      bn(4)    = ( 7.0_ti_p*bn(3) + tmp1%x*alsq2n*exp2a ) *invr2
      bn(5)    = (9.0*bn(4) + tmp1%x*alsq2*alsq2n*exp2a ) *invr2

      bn(0)    = f * bn(0)
      bn(1)    = f * bn(1)
      bn(2)    = f * bn(2)
      bn(3)    = f * bn(3)
      bn(4)    = f * bn(4)
      bn(5)    = f * bn(5)
      end if

      rr1      = f        *invr
      rr3      = f        *invr*invr2
      rr5      = 3.0_ti_p *rr3 *invr2
      rr7      = 5.0_ti_p *rr5 *invr2
      rr9      = 7.0_ti_p *rr7 *invr2
      rr11     = 9.0_ti_p *rr9 *invr2
c
c    Apply Thole polarization damping to scale factors
c
      if (damp.ne.0.0_ti_p) then
         invdamp  = damp**(-one)
         damp1    = - pgamma * (r*invdamp)**2 * (r*invdamp)
         expdamp1 = f_exp(damp1)
         da       = damp1 * expdamp1
c
c     termediates involving Thole damping and scale factors
c
         sc3      = 1.0_ti_p - expdamp1
         sc5      = 1.0_ti_p - (1.0_ti_p - damp1)*expdamp1
         sc7      = 1.0_ti_p - (1.0_ti_p - damp1 + 0.6_ti_p*damp1**2)
     &                           *expdamp1
         rc3x     = - 3.0_ti_p *da *pos%x *invr2
         rc3y     = - 3.0_ti_p *da *pos%y *invr2
         rc3z     = - 3.0_ti_p *da *pos%z *invr2
         rc5x     = - damp1 * rc3x
         rc5y     = - damp1 * rc3y
         rc5z     = - damp1 * rc3z
         rc7x     = - (0.2_ti_p + 0.6_ti_p*damp1)*rc5x
         rc7y     = - (0.2_ti_p + 0.6_ti_p*damp1)*rc5y
         rc7z     = - (0.2_ti_p + 0.6_ti_p*damp1)*rc5z
      else
         sc3 = 1; sc5=1; sc7=1;
         rc3x = 0; rc3y=0; rc3z=0;
         rc5x = 0; rc5y=0; rc5z=0;
         rc7x = 0; rc7y=0; rc7z=0;
      end if

c
c     Intermediates involving moments and distance separation
c
      dri      =    ip%dx*pos%x +   ip%dy*pos%y +  ip%dz*pos%z
      drk      =    kp%dx*pos%x +   kp%dy*pos%y +  kp%dz*pos%z

      qrix     =   ip%qxx*pos%x +  ip%qxy*pos%y +  ip%qxz*pos%z
      qriy     =   ip%qxy*pos%x +  ip%qyy*pos%y +  ip%qyz*pos%z
      qriz     =   ip%qxz*pos%x +  ip%qyz*pos%y +  ip%qzz*pos%z
      qrkx     =   kp%qxx*pos%x +  kp%qxy*pos%y +  kp%qxz*pos%z
      qrky     =   kp%qxy*pos%x +  kp%qyy*pos%y +  kp%qyz*pos%z
      qrkz     =   kp%qxz*pos%x +  kp%qyz*pos%y +  kp%qzz*pos%z
      qrri     =     qrix*pos%x +    qriy*pos%y +    qriz*pos%z
      qrrk     =     qrkx*pos%x +    qrky*pos%y +    qrkz*pos%z

      !ReUse registers
      associate(
     &  dik=>dterm7x,qrrik=>dterm7y,qik=>dterm7z
     & ,qrrx=>dterm3x,qrry=>dterm3y,qrrz=>dterm3z
     & ,qikrx=>dterm6x,qikry=>dterm6y,qikrz=>dterm6z
     & ,qkirx=>qrimodx,qkiry=>qrimody,qkirz=>qrimodz
     & ,diqkx=>qrkmodx,diqky=>qrkmody,diqkz=>qrkmodz
     & ,dkqix=>dterm1,dkqiy=>dterm2,dkqiz=>duik
     & ,diqrk=>quik,dkqri=>uri
     & ,rr1_=>turk5,rr3_=>turk7, rr5_=>turi5
     & ,rr7_=>tksx,rr9_=>tksy,rr11_=>tksz
     & ,term4=>tkcx,term5=>tkcy,term6=>tkcz
     & ,s=>ticx,ds=>ticy,e_=>ticz,de=>turi7
     & ,diqkxr=>tk5x,diqkyr=>tk5y,diqkzr=>tk5z
     & ,dkqixr=>tk5x,dkqiyr=>tk5y,dkqizr=>tk5z
     & ,dqiqkx=>tisx,dqiqky=>tisy,dqiqkz=>tisz
     & ,qikrxr=>ti5x,qikryr=>ti5y,qikrzr=>ti5z
     & ,qkirxr=>ti5x,qkiryr=>ti5y,qkirzr=>ti5z
     & ,dikx=>dterm4x,diky=>dterm4y,dikz=>dterm4z
     & ,dirx=>dterm5x,diry=>dterm5y,dirz=>dterm5z
     & ,dkrx=>dterm5x,dkry=>dterm5y,dkrz=>dterm5z
     & ,qrixr=>dterm7x,qriyr=>dterm7y,qrizr=>dterm7z
     & ,qrkxr=>dterm7x,qrkyr=>dterm7y,qrkzr=>dterm7z )
c
c        intermediates involving moments and distance separation
c
         dik    = ip%dx*kp%dx + ip%dy*kp%dy + ip%dz*kp%dz

         qrrik  = qrix*qrkx + qriy*qrky + qriz*qrkz
         qik    = two*(ip%qxy*kp%qxy + ip%qxz*kp%qxz + ip%qyz*kp%qyz)
     &               + ip%qxx*kp%qxx + ip%qyy*kp%qyy + ip%qzz*kp%qzz

         qrrx   = qrky*qriz - qrkz*qriy
         qrry   = qrkz*qrix - qrkx*qriz
         qrrz   = qrkx*qriy - qrky*qrix

         qikrx  = ip%qxx*qrkx + ip%qxy*qrky + ip%qxz*qrkz
         qikry  = ip%qxy*qrkx + ip%qyy*qrky + ip%qyz*qrkz
         qikrz  = ip%qxz*qrkx + ip%qyz*qrky + ip%qzz*qrkz
         qkirx  = kp%qxx*qrix + kp%qxy*qriy + kp%qxz*qriz
         qkiry  = kp%qxy*qrix + kp%qyy*qriy + kp%qyz*qriz
         qkirz  = kp%qxz*qrix + kp%qyz*qriy + kp%qzz*qriz

         diqkx  = ip%dx*kp%qxx  + ip%dy*kp%qxy + ip%dz*kp%qxz
         diqky  = ip%dx*kp%qxy  + ip%dy*kp%qyy + ip%dz*kp%qyz
         diqkz  = ip%dx*kp%qxz  + ip%dy*kp%qyz + ip%dz*kp%qzz
         dkqix  = kp%dx*ip%qxx  + kp%dy*ip%qxy + kp%dz*ip%qxz
         dkqiy  = kp%dx*ip%qxy  + kp%dy*ip%qyy + kp%dz*ip%qyz
         dkqiz  = kp%dx*ip%qxz  + kp%dy*ip%qyz + kp%dz*ip%qzz
         diqrk  = ip%dx*qrkx  + ip%dy*qrky + ip%dz*qrkz
         dkqri  = kp%dx*qrix  + kp%dy*qriy + kp%dz*qriz

c
c        modify distances to account for Ewald and exclusions
c
         rr1_   = two*(bn(0) - mscale*rr1)
         rr3_   = two*(bn(1) - mscale*rr3)
         rr5_   = two*(bn(2) - mscale*rr5)
         rr7_   = two*(bn(3) - mscale*rr7)
         rr9_   = two*(bn(4) - mscale*rr9)
         rr11_  = two*(bn(5) - mscale*rr11)
c
c        calculate intermediate terms for multipole energy
c
         term1  = ip%c*kp%c
         term2  = kp%c*dri   - ip%c*drk   + dik
         term3  = ip%c*qrrk  + kp%c*qrri  - dri*drk
     &          + two*(dkqri - diqrk + qik)
         term4  = dri*qrrk - drk*qrri - 4.0*qrrik
         term5  = qrri*qrrk
c
c        compute the energy contributions for this interaction
c
         s      = 1.0
         if (mode.eq.m_short.or.mode.eq.m_long) then
            call switch_respa_inl(r,r_cut,shortheal,s,ds)
            if (mode.eq.m_long) s = 1.0-s
         end if

         e_     = term1*rr1_+ term2*rr3_+ term3*rr5_
     &                      + term4*rr7_+ term5*rr9_

         e      = WRITE_C(e + ) ( s*e_ )
         de     = ( term1*rr3_+ term2*rr5_+ term3*rr7_
     &          +   term4*rr9_+ term5*rr11_)
c
c        calculate intermediate terms for force and torque
c
         term1  = -kp%c*rr3_+ drk*rr5_- qrrk*rr7_
         term2  =  ip%c*rr3_+ dri*rr5_+ qrri*rr7_
         term3  = two * rr5_
         term4  = two * (-kp%c*rr5_+drk*rr7_-qrrk*rr9_)
         term5  = two * (-ip%c*rr5_-dri*rr7_-qrri*rr9_)
         term6  = 4.0 * rr7_
c
c        compute the force components for this interaction
c
         frc%x  = de*pos%x + term1*ip%dx  + term2*kp%dx
         frc%y  = de*pos%y + term1*ip%dy  + term2*kp%dy
         frc%z  = de*pos%z + term1*ip%dz  + term2*kp%dz

         frc%x  = frc%x + (term3*(diqkx-dkqix) + term4*qrix)
         frc%y  = frc%y + (term3*(diqky-dkqiy) + term4*qriy)
         frc%z  = frc%z + (term3*(diqkz-dkqiz) + term4*qriz)

         frc%x  = frc%x + (term5*qrkx          + term6*(qikrx+qkirx))
         frc%y  = frc%y + (term5*qrky          + term6*(qikry+qkiry))
         frc%z  = frc%z + (term5*qrkz          + term6*(qikrz+qkirz))

         frc%x  = s*( frc%x )
         frc%y  = s*( frc%y )
         frc%z  = s*( frc%z )

         if (mode.eq.m_short) then
            frc%x  = frc%x - ds*pos%x*e_*invr
            frc%y  = frc%y - ds*pos%y*e_*invr
            frc%z  = frc%z - ds*pos%z*e_*invr
         else if (mode.eq.m_long) then
            frc%x  = frc%x + ds*pos%x*e_*invr
            frc%y  = frc%y + ds*pos%y*e_*invr
            frc%z  = frc%z + ds*pos%z*e_*invr
         end if
c
c        compute the torque components for this interaction
c
         dqiqkx = ip%qxz*kp%qxy + ip%qyz*kp%qyy + ip%qzz*kp%qyz
         dqiqky = ip%qxx*kp%qxz + ip%qxy*kp%qyz + ip%qxz*kp%qzz
         dqiqkz = ip%qxy*kp%qxx + ip%qyy*kp%qxy + ip%qyz*kp%qxz

         dqiqkx =dqiqkx -(ip%qxy*kp%qxz + ip%qyy*kp%qyz + ip%qyz*kp%qzz)
         dqiqky =dqiqky -(ip%qxz*kp%qxx + ip%qyz*kp%qxy + ip%qzz*kp%qxz)
         dqiqkz =dqiqkz -(ip%qxx*kp%qxy + ip%qxy*kp%qyy + ip%qxz*kp%qyz)
         dqiqkx = two*dqiqkx
         dqiqky = two*dqiqky
         dqiqkz = two*dqiqkz
         dqiqkx =dqiqkx + ip%dy*qrkz + kp%dy*qriz
         dqiqky =dqiqky + ip%dz*qrkx + kp%dz*qrix
         dqiqkz =dqiqkz + ip%dx*qrky + kp%dx*qriy

         dkqixr = dkqiz*pos%y  - dkqiy*pos%z
         dkqiyr = dkqix*pos%z  - dkqiz*pos%x
         dkqizr = dkqiy*pos%x  - dkqix*pos%y

         dqiqkx =dqiqkx -(ip%dz*qrky + kp%dz*qriy)
         dqiqky =dqiqky -(ip%dx*qrkz + kp%dx*qriz)
         dqiqkz =dqiqkz -(ip%dy*qrkx + kp%dy*qrix)

         dirx   = ip%dy*pos%z - ip%dz*pos%y
         diry   = ip%dz*pos%x - ip%dx*pos%z
         dirz   = ip%dx*pos%y - ip%dy*pos%x

         depx   = + term1*dirx + term3*(dqiqkx+dkqixr)
         depy   = + term1*diry + term3*(dqiqky+dkqiyr)
         depz   = + term1*dirz + term3*(dqiqkz+dkqizr)
         qrixr  = qriz*pos%y - qriy*pos%z
         qriyr  = qrix*pos%z - qriz*pos%x
         qrizr  = qriy*pos%x - qrix*pos%y
         qikrxr = qikrz*pos%y  - qikry*pos%z
         qikryr = qikrx*pos%z  - qikrz*pos%x
         qikrzr = qikry*pos%x  - qikrx*pos%y
         dikx   = ip%dy*kp%dz - ip%dz*kp%dy
         diky   = ip%dz*kp%dx - ip%dx*kp%dz
         dikz   = ip%dx*kp%dy - ip%dy*kp%dx
         depx   =depx - ( term4*qrixr + term6*(qikrxr+qrrx)) - rr3_*dikx
         depy   =depy - ( term4*qriyr + term6*(qikryr+qrry)) - rr3_*diky
         depz   =depz - ( term4*qrizr + term6*(qikrzr+qrrz)) - rr3_*dikz
         !depx   = 0.01
         !depy   = 0.01
         !depz   = 0.01

         trqi%x = WRITE_C(trqi%x +) s*( depx )
         trqi%y = WRITE_C(trqi%y +) s*( depy )
         trqi%z = WRITE_C(trqi%z +) s*( depz )

         diqkxr = diqkz*pos%y - diqky*pos%z
         diqkyr = diqkx*pos%z - diqkz*pos%x
         diqkzr = diqky*pos%x - diqkx*pos%y
         dkrx   = kp%dy*pos%z - kp%dz*pos%y
         dkry   = kp%dz*pos%x - kp%dx*pos%z
         dkrz   = kp%dx*pos%y - kp%dy*pos%x
         depx   = term2*dkrx - term3*(dqiqkx+diqkxr)
         depy   = term2*dkry - term3*(dqiqky+diqkyr)
         depz   = term2*dkrz - term3*(dqiqkz+diqkzr)
         qkirxr = qkirz*pos%y - qkiry*pos%z
         qkiryr = qkirx*pos%z - qkirz*pos%x
         qkirzr = qkiry*pos%x - qkirx*pos%y
         qrkxr  = qrkz*pos%y - qrky*pos%z
         qrkyr  = qrkx*pos%z - qrkz*pos%x
         qrkzr  = qrky*pos%x - qrkx*pos%y
         depx   =depx - ( term5*qrkxr + term6*(qkirxr-qrrx)) + rr3_*dikx
         depy   =depy - ( term5*qrkyr + term6*(qkiryr-qrry)) + rr3_*diky
         depz   =depz - ( term5*qrkzr + term6*(qkirzr-qrrz)) + rr3_*dikz
         !depx   = 0.01
         !depy   = 0.01
         !depz   = 0.01

         trqk%x = WRITE_C(trqk%x +) s*( depx )
         trqk%y = WRITE_C(trqk%y +) s*( depy )
         trqk%z = WRITE_C(trqk%z +) s*( depz )
      end associate

      if (do_correct) then
         psc3   = sc3 * pscale
         dsc3   = sc3 * dscale
         usc3   = sc3 * uscale
         psc5   = sc5 * pscale
         dsc5   = sc5 * dscale
         usc5   = sc5 * uscale
         psc7   = sc7 * pscale
         dsc7   = sc7 * dscale

         prc3x  = - rc3x * pscale
         drc3x  = - rc3x * dscale
         urc3x  = - rc3x * uscale
         prc3y  = - rc3y * pscale
         drc3y  = - rc3y * dscale
         urc3y  = - rc3y * uscale
         prc3z  = - rc3z * pscale
         drc3z  = - rc3z * dscale
         urc3z  = - rc3z * uscale

         prc5x  = - rc5x * pscale
         drc5x  = - rc5x * dscale
         urc5x  = - rc5x * uscale
         prc5y  = - rc5y * pscale
         drc5y  = - rc5y * dscale
         urc5y  = - rc5y * uscale
         prc5z  = - rc5z * pscale
         drc5z  = - rc5z * dscale
         urc5z  = - rc5z * uscale

         prc7x  = - rc7x * pscale
         drc7x  = - rc7x * dscale
         prc7y  = - rc7y * pscale
         drc7y  = - rc7y * dscale
         prc7z  = - rc7z * pscale
         drc7z  = - rc7z * dscale
      else
         psc3   = 1.0_ti_p - sc3 * pscale
         dsc3   = 1.0_ti_p - sc3 * dscale
         usc3   = 1.0_ti_p - sc3 * uscale
         psc5   = 1.0_ti_p - sc5 * pscale
         dsc5   = 1.0_ti_p - sc5 * dscale
         usc5   = 1.0_ti_p - sc5 * uscale
         psc7   = 1.0_ti_p - sc7 * pscale
         dsc7   = 1.0_ti_p - sc7 * dscale

         prc3x  = rc3x * pscale
         drc3x  = rc3x * dscale
         urc3x  = rc3x * uscale
         prc3y  = rc3y * pscale
         drc3y  = rc3y * dscale
         urc3y  = rc3y * uscale
         prc3z  = rc3z * pscale
         drc3z  = rc3z * dscale
         urc3z  = rc3z * uscale

         prc5x  = rc5x * pscale
         drc5x  = rc5x * dscale
         urc5x  = rc5x * uscale
         prc5y  = rc5y * pscale
         drc5y  = rc5y * dscale
         urc5y  = rc5y * uscale
         prc5z  = rc5z * pscale
         drc5z  = rc5z * dscale
         urc5z  = rc5z * uscale

         prc7x  = rc7x * pscale
         drc7x  = rc7x * dscale
         prc7y  = rc7y * pscale
         drc7y  = rc7y * dscale
         prc7z  = rc7z * pscale
         drc7z  = rc7z * dscale
      end if

      psr3     = bn(1) - psc3*rr3
      dsr3     = bn(1) - dsc3*rr3
      usr3     = bn(1) - usc3*rr3
      psr5     = bn(2) - psc5*rr5
      dsr5     = bn(2) - dsc5*rr5
      usr5     = bn(2) - usc5*rr5
      psr7     = bn(3) - psc7*rr7
      dsr7     = bn(3) - dsc7*rr7

      ! Intermediates terms involving moments and distance separation
      uri      =   dpui%x*pos%x +  dpui%y*pos%y +  dpui%z*pos%z
      urk      =   dpuk%x*pos%x +  dpuk%y*pos%y +  dpuk%z*pos%z
      urip     =  dpui%xx*pos%x + dpui%yy*pos%y + dpui%zz*pos%z
      urkp     =  dpuk%xx*pos%x + dpuk%yy*pos%y + dpuk%zz*pos%z

      duik     =   ip%dx*dpuk%x  +  ip%dy*dpuk%y + ip%dz*dpuk%z
     &          +  kp%dx*dpui%x  +  kp%dy*dpui%y + kp%dz*dpui%z
      quik     =    qrix*dpuk%x  +   qriy*dpuk%y +  qriz*dpuk%z
     &           -  qrkx*dpui%x  -   qrky*dpui%y -  qrkz*dpui%z

c
c     Calculate intermediate terms for polarization interaction
c
      term1    =  kp%c*uri  - ip%c*urk + duik
      term2    =   two*quik -  uri*drk - dri*urk
      term3    =   uri*qrrk -  urk*qrri
c
c     compute the energy contribution for this interaction
c
      e        = e + term1*psr3 + term2*psr5 + term3*psr7

      qrimodx  =                  qriy*pos%y  +  qriz*pos%z
      qrimody  =  qrix*pos%x   +                 qriz*pos%z
      qrimodz  =  qrix*pos%x   +  qriy*pos%y
      qrkmodx  =                  qrky*pos%y  +  qrkz*pos%z
      qrkmody  =  qrkx*pos%x   +                 qrkz*pos%z
      qrkmodz  =  qrkx*pos%x   +  qrky*pos%y
c
c     t the dEd/dR terms used for direct polarization force
c
      dterm1   =  bn(2)   -   dsc3*rr5
      dterm2   =  bn(3)   -   dsc5*rr7
c
c     Straight terms ( xx, yy ,zz )
c
      dterm3x  = - dsr3  + dterm1*pos%x**2 - rr3*pos%x*drc3x
      dterm3y  = - dsr3  + dterm1*pos%y**2 - rr3*pos%y*drc3y
      dterm3z  = - dsr3  + dterm1*pos%z**2 - rr3*pos%z*drc3z

      dterm4x  =  rr3*drc3x  - (dterm1 + dsr5)*pos%x
      dterm4y  =  rr3*drc3y  - (dterm1 + dsr5)*pos%y
      dterm4z  =  rr3*drc3z  - (dterm1 + dsr5)*pos%z

      dterm5x  = - dsr5      + dterm2*pos%x**2 - rr5*pos%x*drc5x
      dterm5y  = - dsr5      + dterm2*pos%y**2 - rr5*pos%y*drc5y
      dterm5z  = - dsr5      + dterm2*pos%z**2 - rr5*pos%z*drc5z

      dterm6x  =  (bn(4) - dsc7*rr9)*pos%x**2  - bn(3) - rr7*pos%x*drc7x
      dterm6y  =  (bn(4) - dsc7*rr9)*pos%y**2  - bn(3) - rr7*pos%y*drc7y
      dterm6z  =  (bn(4) - dsc7*rr9)*pos%z**2  - bn(3) - rr7*pos%z*drc7z

      dterm7x  = rr5*drc5x  - two*bn(3)*pos%x
     &           + ( dsc5   + 1.5_ti_p*dsc7 )*rr7*pos%x
      dterm7y  = rr5*drc5y  - two*bn(3)*pos%y
     &           + ( dsc5   + 1.5_ti_p*dsc7)*rr7*pos%y
      dterm7z  = rr5*drc5z  - two*bn(3)*pos%z
     &           + ( dsc5   + 1.5_ti_p*dsc7 )*rr7*pos%z
c
c     Straight terms ( xx, yy ,zz )
c
      tisx     = ip%c*dterm3x      + ip%dx*dterm4x + dri*dterm5x
     &          + two*dsr5*ip%qxx  + qrimodx*dsc7*rr7
     &          + two*qrix*dterm7x + qrri*dterm6x
      tisy     = ip%c*dterm3y      + ip%dy*dterm4y + dri*dterm5y
     &          + two*dsr5*ip%qyy  + qrimody*dsc7*rr7
     &          + two*qriy*dterm7y + qrri*dterm6y
      tisz     = ip%c*dterm3z      + ip%dz*dterm4z + dri*dterm5z
     &          + two*dsr5*ip%qzz  + qrimodz*dsc7*rr7
     &          + two*qriz*dterm7z +  qrri*dterm6z

      tksx     = kp%c*dterm3x      - kp%dx*dterm4x - drk*dterm5x
     &          + two*dsr5*kp%qxx  + qrkmodx*dsc7*rr7
     &          + two*qrkx*dterm7x + qrrk*dterm6x
      tksy     = kp%c*dterm3y      - kp%dy*dterm4y - drk*dterm5y
     &          + two*dsr5*kp%qyy  + qrkmody*dsc7*rr7
     &          + two*qrky*dterm7y + qrrk*dterm6y
      tksz     = kp%c*dterm3z      - kp%dz*dterm4z  - drk*dterm5z
     &          + two*dsr5*kp%qzz  + qrkmodz*dsc7*rr7
     &          + two*qrkz*dterm7z + qrrk*dterm6z
c
c     Crossed terms ( xy = yx , xz = zx ,yz = zy )
c
      tmp1%x   = pos%x * pos%y
      tmp1%y   = pos%x * pos%z
      tmp1%z   = pos%y * pos%z

      dterm3x  =   dterm1*tmp1%x -  rr3*pos%y*drc3x
      dterm3y  =   dterm1*tmp1%y -  rr3*pos%z*drc3x
      dterm3z  =   dterm1*tmp1%z -  rr3*pos%z*drc3y
      dterm4x  = - dterm1*pos%x  +  rr3*drc3x
      dterm4y  = - dterm1*pos%x  +  rr3*drc3x
      dterm4z  = - dterm1*pos%y  +  rr3*drc3y

      dterm5x  =   dterm2*tmp1%x -  rr5*pos%y*drc5x
      dterm5y  =   dterm2*tmp1%y -  rr5*pos%z*drc5x
      dterm5z  =   dterm2*tmp1%z -  rr5*pos%z*drc5y
      dterm6x  =  (bn(4) - dsc7*rr9)*tmp1%x - rr7*pos%y*drc7x
      dterm6y  =  (bn(4) - dsc7*rr9)*tmp1%y - rr7*pos%z*drc7x
      dterm6z  =  (bn(4) - dsc7*rr9)*tmp1%z - rr7*pos%z*drc7y
      dterm7x  = - dterm2*pos%x  +  rr5*drc5x
      dterm7y  = - dterm2*pos%x  +  rr5*drc5x
      dterm7z  = - dterm2*pos%y  +  rr5*drc5y
c
c     Crossed terms ( xy = yx , xz = zx ,yz = zy )
c

      ticx     = ip%c*dterm3x + ip%dy*dterm4x + dri*dterm5x
     &         + qrri*dterm6x - dsr5*ip%dx*pos%y
     &         +  two*(dsr5*ip%qxy - dsr7*pos%y*qrix + qriy*dterm7x)

      ticy     = ip%c*dterm3y + ip%dz*dterm4y + dri*dterm5y
     &         + qrri*dterm6y - dsr5*ip%dx*pos%z
     &         +  two*(dsr5*ip%qxz - dsr7*pos%z*qrix + qriz*dterm7y)

      ticz     = ip%c*dterm3z + ip%dz*dterm4z + dri*dterm5z
     &         + qrri*dterm6z - dsr5*ip%dy*pos%z
     &         +  two*(dsr5*ip%qyz - dsr7*pos%z*qriy + qriz*dterm7z)

      tkcx     = kp%c*dterm3x - kp%dy*dterm4x - drk*dterm5x
     &         + qrrk*dterm6x + dsr5*kp%dx*pos%y
     &         +  two*(dsr5*kp%qxy - dsr7*pos%y*qrkx + qrky*dterm7x)

      tkcy     = kp%c*dterm3y - kp%dz*dterm4y - drk*dterm5y
     &         + qrrk*dterm6y + dsr5*kp%dx*pos%z
     &         +  two*(dsr5*kp%qxz - dsr7*pos%z*qrkx + qrkz*dterm7y)

      tkcz     = kp%c*dterm3z- kp%dz*dterm4z - drk*dterm5z
     &         + qrrk*dterm6z + dsr5 *kp%dy*pos%z
     &         +  two*(dsr5*kp%qyz - dsr7*pos%z*qrky + qrkz*dterm7z)
c
c      Construct matrixes for dot_product
c      do Dot product
c
      depx     =  tisx*dpuk%xx + ticx*dpuk%yy + ticy*dpuk%zz
     &         -  tksx*dpui%xx - tkcx*dpui%yy - tkcy*dpui%zz
      depy     =  ticx*dpuk%xx + tisy*dpuk%yy + ticz*dpuk%zz
     &         -  tkcx*dpui%xx - tksy*dpui%yy - tkcz*dpui%zz
      depz     =  ticy*dpuk%xx + ticz*dpuk%yy + tisz*dpuk%zz
     &         -  tkcy*dpui%xx - tkcz*dpui%yy - tksz*dpui%zz

      frc%x    = frc%x - depx
      frc%y    = frc%y - depy
      frc%z    = frc%z - depz
c
c     t the dEp/dR terms used for direct polarization force
c

      dterm1   =  bn(2) - psc3 * rr5
      dterm2   =  bn(3) - psc5 * rr7
c
c     Straight terms ( xx, yy ,zz )
c
      dterm3x  = - psr3 + dterm1*pos%x**2 - rr3*pos%x*prc3x
      dterm3y  = - psr3 + dterm1*pos%y**2 - rr3*pos%y*prc3y
      dterm3z  = - psr3 + dterm1*pos%z**2 - rr3*pos%z*prc3z

      dterm4x  = - psr5*pos%x - dterm1*pos%x + rr3*prc3x
      dterm4y  = - psr5*pos%y - dterm1*pos%y + rr3*prc3y
      dterm4z  = - psr5*pos%z - dterm1*pos%z + rr3*prc3z

      dterm5x  = - psr5 + dterm2*pos%x**2 - rr5*pos%x*prc5x
      dterm5y  = - psr5 + dterm2*pos%y**2 - rr5*pos%y*prc5y
      dterm5z  = - psr5 + dterm2*pos%z**2 - rr5*pos%z*prc5z

      dterm6x  =  (bn(4) - psc7*rr9)*pos%x**2 - bn(3) - rr7*pos%x*prc7x
      dterm6y  =  (bn(4) - psc7*rr9)*pos%y**2 - bn(3) - rr7*pos%y*prc7y
      dterm6z  =  (bn(4) - psc7*rr9)*pos%z**2 - bn(3) - rr7*pos%z*prc7z

      dterm7x  =  rr5*prc5x - two*bn(3)*pos%x
     &         +  (psc5 + 1.5_ti_p*psc7)*rr7*pos%x
      dterm7y  =  rr5*prc5y - two*bn(3)*pos%y
     &         +  (psc5 + 1.5_ti_p*psc7)*rr7*pos%y
      dterm7z  =  rr5*prc5z - two*bn(3)*pos%z
     &         +  (psc5 + 1.5_ti_p*psc7)*rr7*pos%z
c
c     Straight terms ( xx, yy ,zz )
c
      tisx = ip%c*dterm3x + ip%dx*dterm4x + dri*dterm5x
     &     + qrri*dterm6x + qrimodx*psc7*rr7
     &     + two*(psr5*ip%qxx + qrix*dterm7x)

      tisy = ip%c*dterm3y + ip%dy*dterm4y + dri*dterm5y
     &     + qrri*dterm6y + qrimody*psc7*rr7
     &     + two*(psr5*ip%qyy + qriy*dterm7y)

      tisz = ip%c*dterm3z + ip%dz*dterm4z + dri*dterm5z
     &     + qrri*dterm6z + qrimodz*psc7*rr7
     &     + two*(psr5*ip%qzz + qriz*dterm7z)

      tksx = kp%c*dterm3x - kp%dx*dterm4x - drk*dterm5x
     &     + qrrk*dterm6x + qrkmodx*psc7*rr7
     &     + two*(psr5*kp%qxx + qrkx*dterm7x)

      tksy = kp%c*dterm3y - kp%dy*dterm4y - drk*dterm5y
     &     + qrrk*dterm6y + qrkmody*psc7*rr7
     &     + two*(psr5*kp%qyy + qrky*dterm7y)

      tksz = kp%c*dterm3z - kp%dz*dterm4z - drk*dterm5z
     &     + qrrk*dterm6z + qrkmodz*psc7*rr7
     &     + two*(psr5*kp%qzz + qrkz*dterm7z)
c
c     Crossed terms ( xy = yx , xz = zx ,yz = zy )
c
      tmp1%x  = pos%x * pos%y
      tmp1%y  = pos%x * pos%z
      tmp1%z  = pos%y * pos%z

      dterm3x =  dterm1*tmp1%x - rr3*pos%y*prc3x
      dterm3y =  dterm1*tmp1%y - rr3*pos%z*prc3x
      dterm3z =  dterm1*tmp1%z - rr3*pos%z*prc3y

      dterm4x =  rr3*prc3x - dterm1*pos%x
      dterm4y =  rr3*prc3x - dterm1*pos%x
      dterm4z =  rr3*prc3y - dterm1*pos%y
      dterm5x =  dterm2*tmp1%x - rr5*pos%y*prc5x
      dterm5y =  dterm2*tmp1%y - rr5*pos%z*prc5x
      dterm5z =  dterm2*tmp1%z - rr5*pos%z*prc5y

      dterm6x =  (bn(4) - psc7*rr9)*tmp1%x - rr7*pos%y*prc7x
      dterm6y =  (bn(4) - psc7*rr9)*tmp1%y - rr7*pos%z*prc7x
      dterm6z =  (bn(4) - psc7*rr9)*tmp1%z - rr7*pos%z*prc7y
      dterm7x =  rr5*prc5x - dterm2*pos%x
      dterm7y =  rr5*prc5x - dterm2*pos%x
      dterm7z =  rr5*prc5y - dterm2*pos%y

c
c     Crossed terms ( xy = yx , xz = zx ,yz = zy )
c
      ticx = ip%c*dterm3x + ip%dy*dterm4x + dri*dterm5x
     &     + qrri*dterm6x - psr5*ip%dx*pos%y
     &     + two*(psr5*ip%qxy - psr7*pos%y*qrix + qriy*dterm7x)

      ticy = ip%c*dterm3y + ip%dz*dterm4y + dri*dterm5y
     &     + qrri*dterm6y - psr5*ip%dx*pos%z
     &     + two*(psr5*ip%qxz - psr7*pos%z*qrix + qriz*dterm7y)

      ticz = ip%c*dterm3z + ip%dz*dterm4z + dri*dterm5z
     &     + qrri*dterm6z - psr5*ip%dy*pos%z
     &     +  two*(psr5*ip%qyz - psr7*pos%z*qriy + qriz*dterm7z)

      tkcx = kp%c*dterm3x - kp%dy*dterm4x - drk*dterm5x
     &     + qrrk*dterm6x + psr5*kp%dx*pos%y
     &     +  two*(psr5*kp%qxy - psr7*pos%y*qrkx + qrky*dterm7x)

      tkcy = kp%c*dterm3y - kp%dz*dterm4y - drk*dterm5y
     &     + qrrk*dterm6y + psr5*kp%dx*pos%z
     &     + two*(psr5*kp%qxz - psr7*pos%z*qrkx + qrkz*dterm7y)

      tkcz = kp%c*dterm3z - kp%dz*dterm4z - drk*dterm5z
     &     + qrrk*dterm6z + psr5*kp%dy*pos%z
     &     + two*(psr5*kp%qyz - psr7*pos%z*qrky + qrkz*dterm7z)
c
c     Construct matrixes for dot_product
c     Do dot product
c
      depx = tisx*dpuk%x + ticx*dpuk%y + ticy*dpuk%z
     &     - tksx*dpui%x - tkcx*dpui%y - tkcy*dpui%z
      depy = ticx*dpuk%x + tisy*dpuk%y + ticz*dpuk%z
     &     - tkcx*dpui%x - tksy*dpui%y - tkcz*dpui%z
      depz = ticy*dpuk%x + ticz*dpuk%y + tisz*dpuk%z
     &     - tkcy*dpui%x - tkcz*dpui%y - tksz*dpui%z

      frc%x = frc%x - depx
      frc%y = frc%y - depy
      frc%z = frc%z - depz

c
c     t the dtau/dr terms used for mutual polarization force
c
      dterm1 = bn(2) - usc3 * rr5
      dterm2 = bn(3) - usc5 * rr7
c
c     Straight terms ( xx, yy ,zz )
c
      dterm3x = usr5 + dterm1
      dterm3y = usr5 + dterm1
      dterm3z = usr5 + dterm1
      dterm4x = rr3  * uscale
      dterm4y = rr3  * uscale
      dterm4z = rr3  * uscale

      dterm5x = - pos%x*dterm3x + rc3x*dterm4x
      dterm5y = - pos%y*dterm3y + rc3y*dterm4y
      dterm5z = - pos%z*dterm3z + rc3z*dterm4z

      dterm6x = - usr5 + pos%x**2*dterm2 - rr5*pos%x*urc5x
      dterm6y = - usr5 + pos%y**2*dterm2 - rr5*pos%y*urc5y
      dterm6z = - usr5 + pos%z**2*dterm2 - rr5*pos%z*urc5z

      tisx =  dpui%x*dterm5x + uri*dterm6x
      tisy =  dpui%y*dterm5y + uri*dterm6y
      tisz =  dpui%z*dterm5z + uri*dterm6z
      tksx =  dpuk%x*dterm5x + urk*dterm6x
      tksy =  dpuk%y*dterm5y + urk*dterm6y
      tksz =  dpuk%z*dterm5z + urk*dterm6z
c
c     Crossed terms ( xy = yx , xz = zx ,yz = zy )
c
      dterm4x = - usr5*pos%y
      dterm4y = - usr5*pos%z
      dterm4z = - usr5*pos%z
      dterm5x = - pos%x*dterm1 + rr3*urc3x
      dterm5y = - pos%x*dterm1 + rr3*urc3x
      dterm5z = - pos%y*dterm1 + rr3*urc3y
      dterm6x =   pos%x*dterm2*pos%y - rr5*pos%y*urc5x
      dterm6y =   pos%x*dterm2*pos%z - rr5*pos%z*urc5x
      dterm6z =   pos%y*dterm2*pos%z - rr5*pos%z*urc5y

      ticx =  dpui%x*dterm4x + dpui%y*dterm5x + uri*dterm6x
      ticy =  dpui%x*dterm4y + dpui%z*dterm5y + uri*dterm6y
      ticz =  dpui%y*dterm4z + dpui%z*dterm5z + uri*dterm6z

      tkcx =  dpuk%x*dterm4x + dpuk%y*dterm5x + urk*dterm6x
      tkcy =  dpuk%x*dterm4y + dpuk%z*dterm5y + urk*dterm6y
      tkcz =  dpuk%y*dterm4z + dpuk%z*dterm5z + urk*dterm6z
c
c     Construct matrixes for dot_product
c     Dot product
c
      depx =  tisx*dpuk%xx + ticx*dpuk%yy + ticy*dpuk%zz
     &      + tksx*dpui%xx + tkcx*dpui%yy + tkcy*dpui%zz
      depy =  ticx*dpuk%xx + tisy*dpuk%yy + ticz*dpuk%zz
     &      + tkcx*dpui%xx + tksy*dpui%yy + tkcz*dpui%zz
      depz =  ticy*dpuk%xx + ticz*dpuk%yy + tisz*dpuk%zz
     &      + tkcy*dpui%xx + tkcz*dpui%yy + tksz*dpui%zz

      frc%x   = frc%x - depx
      frc%y   = frc%y - depy
      frc%z   = frc%z - depz
      frc_r%x = WRITE_C(frc_r%x -) tp2mdr(frc%x)
      frc_r%y = WRITE_C(frc_r%y -) tp2mdr(frc%y)
      frc_r%z = WRITE_C(frc_r%z -) tp2mdr(frc%z)
c
c     t the induced dipole field used for dipole torques
c
      turi5    = -psr5*urk - dsr5*urkp
      turk5    = -psr5*uri - dsr5*urip
      turi7    = - psr7*urk - dsr7*urkp
      turk7    = - psr7*uri - dsr7*urip
c
c     t induced dipole field gradient used for quadrupole torques
c
      ti5x     =  two*(psr5*dpuk%x + dsr5*dpuk%xx)
      ti5y     =  two*(psr5*dpuk%y + dsr5*dpuk%yy)
      ti5z     =  two*(psr5*dpuk%z + dsr5*dpuk%zz)
      tk5x     =  two*(psr5*dpui%x + dsr5*dpui%xx)
      tk5y     =  two*(psr5*dpui%y + dsr5*dpui%yy)
      tk5z     =  two*(psr5*dpui%z + dsr5*dpui%zz)
c
c     Torque is induced field and gradient cross permanent moments
c
      ufli%x   = psr3*dpuk%x + dsr3*dpuk%xx + pos%x*turi5
      ufli%y   = psr3*dpuk%y + dsr3*dpuk%yy + pos%y*turi5
      ufli%z   = psr3*dpuk%z + dsr3*dpuk%zz + pos%z*turi5
      uflk%x   = psr3*dpui%x + dsr3*dpui%xx + pos%x*turk5
      uflk%y   = psr3*dpui%y + dsr3*dpui%yy + pos%y*turk5
      uflk%z   = psr3*dpui%z + dsr3*dpui%zz + pos%z*turk5

      dufli%x  = pos%x*ti5x + pos%x**2*turi7
      dufli%z  = pos%y*ti5y + pos%y**2*turi7
      dufli%zz = pos%z*ti5z + pos%z**2*turi7
      dufli%y  = pos%x*ti5y + pos%y*ti5x + two*pos%x*pos%y*turi7
      dufli%xx = pos%x*ti5z + pos%z*ti5x + two*pos%x*pos%z*turi7
      dufli%yy = pos%y*ti5z + pos%z*ti5y + two*pos%y*pos%z*turi7

      duflk%x  = -(pos%x*tk5x + pos%x**2*turk7)
      duflk%z  = -(pos%y*tk5y + pos%y**2*turk7)
      duflk%zz = -(pos%z*tk5z + pos%z**2*turk7)
      duflk%y  = -(pos%x*tk5y + pos%y*tk5x + two*pos%x*pos%y*turk7)
      duflk%xx = -(pos%x*tk5z + pos%z*tk5x + two*pos%x*pos%z*turk7)
      duflk%yy = -(pos%y*tk5z + pos%z*tk5y + two*pos%y*pos%z*turk7)

      trqi%x = trqi%x + ip%dz*ufli%y - ip%dy*ufli%z
     &       + ip%qxz*dufli%y - ip%qxy*dufli%xx
     &       + (ip%qzz - ip%qyy)*dufli%yy
     &       + two*ip%qyz*( dufli%z - dufli%zz )
      trqi%y = trqi%y + ip%dx*ufli%z - ip%dz*ufli%x
     &       - ip%qyz*dufli%y
     &       + ip%qxy*dufli%yy + (ip%qxx - ip%qzz)*dufli%xx
     &       + two*ip%qxz*( dufli%zz - dufli%x )
      trqi%z = trqi%z + ip%dy*ufli%x - ip%dx*ufli%y
     &       + ip%qyz*dufli%xx
     &       - ip%qxz*dufli%yy + (ip%qyy - ip%qxx)*dufli%y
     &       + two*ip%qxy*( dufli%x - dufli%z )

      trqk%x = trqk%x + kp%dz*uflk%y - kp%dy*uflk%z
     &       + kp%qxz*duflk%y
     &       - kp%qxy*duflk%xx + (kp%qzz - kp%qyy)*duflk%yy
     &       + two*kp%qyz*( duflk%z - duflk%zz )
      trqk%y = trqk%y + kp%dx*uflk%z - kp%dz*uflk%x
     &       - kp%qyz*duflk%y
     &       + kp%qxy*duflk%yy + (kp%qxx - kp%qzz)*duflk%xx
     &       + two*kp%qxz*( duflk%zz - duflk%x )
      trqk%z = trqk%z + kp%dy*uflk%x - kp%dx*uflk%y
     &       + kp%qyz*duflk%xx
     &       - kp%qxz*duflk%yy + (kp%qyy - kp%qxx)*duflk%y
     &       + two*kp%qxy*( duflk%x - duflk%z )
      end subroutine

      ! Same routine as mpolar1_couple
      ! excpept no scaling is taken under consideration
      M_subroutine
     &               mpolar1_couple_comp(dpui,ip,dpuk,kp,r2,pos
     &                  ,aewald,alsq2,alsq2n,pgamma,damp,f
     &                  ,r_cut,shortheal
     &                  ,mscale,dscale,pscale,uscale
     &                  ,e,frc,frc_r,trqi,trqk,do_correct,mode
     &     ,iglob,kglob)
!$acc routine
      use tinheader ,only: ti_p
      use tinTypes  ,only: rpole_elt,real3,real6,mdyn3_r
      use interfaces,only: m_normal,m_short,m_long
#ifdef TINKER_CUF
      use utilcu    ,only: f_erfc
#  if defined(SINGLE)||defined(MIXED)
      use utilcu    ,only: f_sqrt,f_exp
#  endif
#endif
      implicit none

      integer iglob,kglob
      real(t_p)  ,intent(in),value:: aewald,alsq2,alsq2n,pgamma,damp,f
     &           ,r_cut,shortheal
      real(t_p)  ,intent(in):: mscale,dscale,pscale,uscale
      real(t_p)  ,intent(in):: r2
      type(rpole_elt),intent(in)::ip,kp
      type(real3),intent(in):: pos
      type(real6),intent(in):: dpuk,dpui
      integer    ,intent(in):: mode
      logical    ,intent(in):: do_correct
      real(t_p)  ,intent(inout):: e
      type(real3),intent(inout):: frc
      type(real3),intent(inout):: trqi,trqk
      type(mdyn3_r),intent(inout)::frc_r

      real(t_p) exp2a,ralpha
      real(t_p) one,two,half
      real(t_p) ck,dkx,dky,dkz
      real(t_p) r,invr2,invr
      real(t_p) bn(0:5)
      real(t_p) rr1,rr3,rr5,rr7,rr9,rr11
      real(t_p) damp1,invdamp
      real(t_p) expdamp,expdamp1
      real(t_p) da
      real(t_p) ukx,uky,ukz,ukpx,ukpy,ukpz
      real(t_p) sc3,sc5,sc7
      real(t_p) sr3,sr5,sr7
      real(t_p) rc3,rc3x,rc3y,rc3z
      real(t_p) rc5,rc5x,rc5y,rc5z
      real(t_p) rc7,rc7x,rc7y,rc7z
      real(t_p) dri,drk
      real(t_p) qrix,qriy,qriz
      real(t_p) qrkx,qrky,qrkz
      real(t_p) qrri,qrrk
      real(t_p) uri,urip,urk,urkp,duik,quik
      real(t_p) qrimodx,qrimody,qrimodz
      real(t_p) qrkmodx,qrkmody,qrkmodz
      real(t_p) term1,term2,term3
      real(t_p) dterm1,dterm2
      real(t_p) dterm3x,dterm3y,dterm3z
      real(t_p) dterm4x,dterm4y,dterm4z
      real(t_p) dterm5x,dterm5y,dterm5z
      real(t_p) dterm6x,dterm6y,dterm6z
      real(t_p) dterm7x,dterm7y,dterm7z
      type(real3) tmp1
      real(t_p) tisx,tisy,tisz,ticx,ticy,ticz
      real(t_p) tkcx,tkcy,tkcz,tksx,tksy,tksz
      real(t_p) ti5x,ti5y,ti5z,tk5x,tk5y,tk5z
      real(t_p) turi5,turi7,turk5,turk7
      real(t_p) depx,depy,depz
      type(real6) :: dufli,duflk
      type(real3) :: ufli,uflk
      parameter(half=0.5)
      parameter(one=1.0, two=2.0)
c
c     t reciprocal distance terms for this interaction
c
      invr2    = r2**(-1)
      r        = f_sqrt(r2)
      invr     = f_sqrt(invr2)
c
c     Calculate the real space Ewald error function terms
c
      ralpha   = aewald * r
      !call erfcore_inl(ralpha, bn(0),1)
      exp2a    = f_exp( - ralpha**2)
      bn(0)    = f_erfc(ralpha)

      bn(0)    = bn(0) * invr
      bn(1)    = ( 1.0_ti_p*bn(0) + alsq2    *alsq2n*exp2a ) * invr2
      bn(2)    = ( 3.0_ti_p*bn(1) + alsq2**2 *alsq2n*exp2a ) * invr2
      bn(3)    = ( 5.0_ti_p*bn(2) + alsq2**2*alsq2 *alsq2n*exp2a )
     &         * invr2
      tmp1%x   = alsq2**2*alsq2**2
      bn(4)    = ( 7.0_ti_p*bn(3) + tmp1%x*alsq2n*exp2a ) *invr2
      bn(5)    = (9.0*bn(4) + tmp1%x*alsq2*alsq2n*exp2a ) *invr2

      bn(0)    = f * bn(0)
      bn(1)    = f * bn(1)
      bn(2)    = f * bn(2)
      bn(3)    = f * bn(3)
      bn(4)    = f * bn(4)
      bn(5)    = f * bn(5)

      rr1      = f        *invr
      rr3      = f        *invr*invr2
      rr5      = 3.0_ti_p *rr3 *invr2
      rr7      = 5.0_ti_p *rr5 *invr2
      rr9      = 7.0_ti_p *rr7 *invr2
      rr11     = 9.0_ti_p *rr9 *invr2
c
c    Apply Thole polarization damping to scale factors
c
      if (damp.ne.0.0_ti_p) then
         invdamp  = damp**(-one)
         damp1    = - pgamma * (r*invdamp)**2 * (r*invdamp)
         expdamp1 = f_exp(damp1)
         da       = damp1 * expdamp1
c
c     termediates involving Thole damping and scale factors
c
         sc3      = 1.0_ti_p - expdamp1
         sc5      = 1.0_ti_p - (1.0_ti_p - damp1)*expdamp1
         sc7      = 1.0_ti_p - (1.0_ti_p - damp1 + 0.6_ti_p*damp1**2)
     &                           *expdamp1
         rc3x     = - 3.0_ti_p *da *pos%x *invr2
         rc3y     = - 3.0_ti_p *da *pos%y *invr2
         rc3z     = - 3.0_ti_p *da *pos%z *invr2
         rc5x     = - damp1 * rc3x
         rc5y     = - damp1 * rc3y
         rc5z     = - damp1 * rc3z
         rc7x     = - (0.2_ti_p + 0.6_ti_p*damp1)*rc5x
         rc7y     = - (0.2_ti_p + 0.6_ti_p*damp1)*rc5y
         rc7z     = - (0.2_ti_p + 0.6_ti_p*damp1)*rc5z
      else
         sc3 = 1; sc5=1; sc7=1;
         rc3x = 0; rc3y=0; rc3z=0;
         rc5x = 0; rc5y=0; rc5z=0;
         rc7x = 0; rc7y=0; rc7z=0;
      end if

c
c     Intermediates involving moments and distance separation
c
      dri      =    ip%dx*pos%x +   ip%dy*pos%y +  ip%dz*pos%z
      drk      =    kp%dx*pos%x +   kp%dy*pos%y +  kp%dz*pos%z

      qrix     =   ip%qxx*pos%x +  ip%qxy*pos%y +  ip%qxz*pos%z
      qriy     =   ip%qxy*pos%x +  ip%qyy*pos%y +  ip%qyz*pos%z
      qriz     =   ip%qxz*pos%x +  ip%qyz*pos%y +  ip%qzz*pos%z
      qrkx     =   kp%qxx*pos%x +  kp%qxy*pos%y +  kp%qxz*pos%z
      qrky     =   kp%qxy*pos%x +  kp%qyy*pos%y +  kp%qyz*pos%z
      qrkz     =   kp%qxz*pos%x +  kp%qyz*pos%y +  kp%qzz*pos%z
      qrri     =     qrix*pos%x +    qriy*pos%y +    qriz*pos%z
      qrrk     =     qrkx*pos%x +    qrky*pos%y +    qrkz*pos%z

      !ReUse registers
      associate(
     &  dik=>dterm7x,qrrik=>dterm7y,qik=>dterm7z
     & ,qrrx=>dterm3x,qrry=>dterm3y,qrrz=>dterm3z
     & ,qikrx=>dterm6x,qikry=>dterm6y,qikrz=>dterm6z
     & ,qkirx=>qrimodx,qkiry=>qrimody,qkirz=>qrimodz
     & ,diqkx=>qrkmodx,diqky=>qrkmody,diqkz=>qrkmodz
     & ,dkqix=>dterm1,dkqiy=>dterm2,dkqiz=>duik
     & ,diqrk=>quik,dkqri=>uri
     & ,rr1_=>turk5,rr3_=>turk7, rr5_=>turi5
     & ,rr7_=>tksx,rr9_=>tksy,rr11_=>tksz
     & ,term4=>tkcx,term5=>tkcy,term6=>tkcz
     & ,s=>ticx,ds=>ticy,e_=>ticz,de=>turi7
     & ,diqkxr=>tk5x,diqkyr=>tk5y,diqkzr=>tk5z
     & ,dkqixr=>tk5x,dkqiyr=>tk5y,dkqizr=>tk5z
     & ,dqiqkx=>tisx,dqiqky=>tisy,dqiqkz=>tisz
     & ,qikrxr=>ti5x,qikryr=>ti5y,qikrzr=>ti5z
     & ,qkirxr=>ti5x,qkiryr=>ti5y,qkirzr=>ti5z
     & ,dikx=>dterm4x,diky=>dterm4y,dikz=>dterm4z
     & ,dirx=>dterm5x,diry=>dterm5y,dirz=>dterm5z
     & ,dkrx=>dterm5x,dkry=>dterm5y,dkrz=>dterm5z
     & ,qrixr=>dterm7x,qriyr=>dterm7y,qrizr=>dterm7z
     & ,qrkxr=>dterm7x,qrkyr=>dterm7y,qrkzr=>dterm7z )
c
c        intermediates involving moments and distance separation
c
         dik    = ip%dx*kp%dx + ip%dy*kp%dy + ip%dz*kp%dz

         qrrik  = qrix*qrkx + qriy*qrky + qriz*qrkz
         qik    = two*(ip%qxy*kp%qxy + ip%qxz*kp%qxz + ip%qyz*kp%qyz)
     &               + ip%qxx*kp%qxx + ip%qyy*kp%qyy + ip%qzz*kp%qzz

         qrrx   = qrky*qriz - qrkz*qriy
         qrry   = qrkz*qrix - qrkx*qriz
         qrrz   = qrkx*qriy - qrky*qrix

         qikrx  = ip%qxx*qrkx + ip%qxy*qrky + ip%qxz*qrkz
         qikry  = ip%qxy*qrkx + ip%qyy*qrky + ip%qyz*qrkz
         qikrz  = ip%qxz*qrkx + ip%qyz*qrky + ip%qzz*qrkz
         qkirx  = kp%qxx*qrix + kp%qxy*qriy + kp%qxz*qriz
         qkiry  = kp%qxy*qrix + kp%qyy*qriy + kp%qyz*qriz
         qkirz  = kp%qxz*qrix + kp%qyz*qriy + kp%qzz*qriz

         diqkx  = ip%dx*kp%qxx  + ip%dy*kp%qxy + ip%dz*kp%qxz
         diqky  = ip%dx*kp%qxy  + ip%dy*kp%qyy + ip%dz*kp%qyz
         diqkz  = ip%dx*kp%qxz  + ip%dy*kp%qyz + ip%dz*kp%qzz
         dkqix  = kp%dx*ip%qxx  + kp%dy*ip%qxy + kp%dz*ip%qxz
         dkqiy  = kp%dx*ip%qxy  + kp%dy*ip%qyy + kp%dz*ip%qyz
         dkqiz  = kp%dx*ip%qxz  + kp%dy*ip%qyz + kp%dz*ip%qzz
         diqrk  = ip%dx*qrkx  + ip%dy*qrky + ip%dz*qrkz
         dkqri  = kp%dx*qrix  + kp%dy*qriy + kp%dz*qriz

c
c        modify distances to account for Ewald and exclusions
c
         rr1_   = two*bn(0)
         rr3_   = two*bn(1)
         rr5_   = two*bn(2)
         rr7_   = two*bn(3)
         rr9_   = two*bn(4)
         rr11_  = two*bn(5)
c
c        calculate intermediate terms for multipole energy
c
         term1  = ip%c*kp%c
         term2  = kp%c*dri   - ip%c*drk   + dik
         term3  = ip%c*qrrk  + kp%c*qrri  - dri*drk
     &          + two*(dkqri - diqrk + qik)
         term4  = dri*qrrk - drk*qrri - 4.0*qrrik
         term5  = qrri*qrrk
c
c        compute the energy contributions for this interaction
c
         s      = 1.0
         if (mode.eq.m_short.or.mode.eq.m_long) then
            call switch_respa_inl(r,r_cut,shortheal,s,ds)
            if (mode.eq.m_long) s = 1.0-s
         end if

         e_     = term1*rr1_+ term2*rr3_+ term3*rr5_
     &                      + term4*rr7_+ term5*rr9_

         e      = WRITE_C(e + ) ( s*e_ )
         de     = ( term1*rr3_+ term2*rr5_+ term3*rr7_
     &          +   term4*rr9_+ term5*rr11_)
c
c        calculate intermediate terms for force and torque
c
         term1  = -kp%c*rr3_+ drk*rr5_- qrrk*rr7_
         term2  =  ip%c*rr3_+ dri*rr5_+ qrri*rr7_
         term3  = two * rr5_
         term4  = two * (-kp%c*rr5_+drk*rr7_-qrrk*rr9_)
         term5  = two * (-ip%c*rr5_-dri*rr7_-qrri*rr9_)
         term6  = 4.0 * rr7_
c
c        compute the force components for this interaction
c
         frc%x  = de*pos%x + term1*ip%dx  + term2*kp%dx
         frc%y  = de*pos%y + term1*ip%dy  + term2*kp%dy
         frc%z  = de*pos%z + term1*ip%dz  + term2*kp%dz

         frc%x  = frc%x + (term3*(diqkx-dkqix) + term4*qrix)
         frc%y  = frc%y + (term3*(diqky-dkqiy) + term4*qriy)
         frc%z  = frc%z + (term3*(diqkz-dkqiz) + term4*qriz)

         frc%x  = frc%x + (term5*qrkx          + term6*(qikrx+qkirx))
         frc%y  = frc%y + (term5*qrky          + term6*(qikry+qkiry))
         frc%z  = frc%z + (term5*qrkz          + term6*(qikrz+qkirz))

         frc%x  = s*( frc%x )
         frc%y  = s*( frc%y )
         frc%z  = s*( frc%z )

         if (mode.eq.m_short) then
            frc%x  = frc%x - ds*pos%x*e_*invr
            frc%y  = frc%y - ds*pos%y*e_*invr
            frc%z  = frc%z - ds*pos%z*e_*invr
         else if (mode.eq.m_long) then
            frc%x  = frc%x + ds*pos%x*e_*invr
            frc%y  = frc%y + ds*pos%y*e_*invr
            frc%z  = frc%z + ds*pos%z*e_*invr
         end if
c
c        compute the torque components for this interaction
c
         dqiqkx = ip%qxz*kp%qxy + ip%qyz*kp%qyy + ip%qzz*kp%qyz
         dqiqky = ip%qxx*kp%qxz + ip%qxy*kp%qyz + ip%qxz*kp%qzz
         dqiqkz = ip%qxy*kp%qxx + ip%qyy*kp%qxy + ip%qyz*kp%qxz

         dqiqkx =dqiqkx -(ip%qxy*kp%qxz + ip%qyy*kp%qyz + ip%qyz*kp%qzz)
         dqiqky =dqiqky -(ip%qxz*kp%qxx + ip%qyz*kp%qxy + ip%qzz*kp%qxz)
         dqiqkz =dqiqkz -(ip%qxx*kp%qxy + ip%qxy*kp%qyy + ip%qxz*kp%qyz)
         dqiqkx = two*dqiqkx
         dqiqky = two*dqiqky
         dqiqkz = two*dqiqkz
         dqiqkx =dqiqkx + ip%dy*qrkz + kp%dy*qriz
         dqiqky =dqiqky + ip%dz*qrkx + kp%dz*qrix
         dqiqkz =dqiqkz + ip%dx*qrky + kp%dx*qriy

         dkqixr = dkqiz*pos%y  - dkqiy*pos%z
         dkqiyr = dkqix*pos%z  - dkqiz*pos%x
         dkqizr = dkqiy*pos%x  - dkqix*pos%y

         dqiqkx =dqiqkx -(ip%dz*qrky + kp%dz*qriy)
         dqiqky =dqiqky -(ip%dx*qrkz + kp%dx*qriz)
         dqiqkz =dqiqkz -(ip%dy*qrkx + kp%dy*qrix)

         dirx   = ip%dy*pos%z - ip%dz*pos%y
         diry   = ip%dz*pos%x - ip%dx*pos%z
         dirz   = ip%dx*pos%y - ip%dy*pos%x

         depx   = + term1*dirx + term3*(dqiqkx+dkqixr)
         depy   = + term1*diry + term3*(dqiqky+dkqiyr)
         depz   = + term1*dirz + term3*(dqiqkz+dkqizr)
         qrixr  = qriz*pos%y - qriy*pos%z
         qriyr  = qrix*pos%z - qriz*pos%x
         qrizr  = qriy*pos%x - qrix*pos%y
         qikrxr = qikrz*pos%y  - qikry*pos%z
         qikryr = qikrx*pos%z  - qikrz*pos%x
         qikrzr = qikry*pos%x  - qikrx*pos%y
         dikx   = ip%dy*kp%dz - ip%dz*kp%dy
         diky   = ip%dz*kp%dx - ip%dx*kp%dz
         dikz   = ip%dx*kp%dy - ip%dy*kp%dx
         depx   =depx - ( term4*qrixr + term6*(qikrxr+qrrx)) - rr3_*dikx
         depy   =depy - ( term4*qriyr + term6*(qikryr+qrry)) - rr3_*diky
         depz   =depz - ( term4*qrizr + term6*(qikrzr+qrrz)) - rr3_*dikz
         !depx   = 0.01
         !depy   = 0.01
         !depz   = 0.01

         trqi%x = WRITE_C(trqi%x +) s*( depx )
         trqi%y = WRITE_C(trqi%y +) s*( depy )
         trqi%z = WRITE_C(trqi%z +) s*( depz )

         diqkxr = diqkz*pos%y - diqky*pos%z
         diqkyr = diqkx*pos%z - diqkz*pos%x
         diqkzr = diqky*pos%x - diqkx*pos%y
         dkrx   = kp%dy*pos%z - kp%dz*pos%y
         dkry   = kp%dz*pos%x - kp%dx*pos%z
         dkrz   = kp%dx*pos%y - kp%dy*pos%x
         depx   = term2*dkrx - term3*(dqiqkx+diqkxr)
         depy   = term2*dkry - term3*(dqiqky+diqkyr)
         depz   = term2*dkrz - term3*(dqiqkz+diqkzr)
         qkirxr = qkirz*pos%y - qkiry*pos%z
         qkiryr = qkirx*pos%z - qkirz*pos%x
         qkirzr = qkiry*pos%x - qkirx*pos%y
         qrkxr  = qrkz*pos%y - qrky*pos%z
         qrkyr  = qrkx*pos%z - qrkz*pos%x
         qrkzr  = qrky*pos%x - qrkx*pos%y
         depx   =depx - ( term5*qrkxr + term6*(qkirxr-qrrx)) + rr3_*dikx
         depy   =depy - ( term5*qrkyr + term6*(qkiryr-qrry)) + rr3_*diky
         depz   =depz - ( term5*qrkzr + term6*(qkirzr-qrrz)) + rr3_*dikz
         !depx   = 0.01
         !depy   = 0.01
         !depz   = 0.01

         trqk%x = WRITE_C(trqk%x +) s*( depx )
         trqk%y = WRITE_C(trqk%y +) s*( depy )
         trqk%z = WRITE_C(trqk%z +) s*( depz )
      end associate

      sc3   = 1.0_ti_p - sc3
      sc5   = 1.0_ti_p - sc5
      sc7   = 1.0_ti_p - sc7

      sr3     = bn(1) - sc3*rr3
      sr5     = bn(2) - sc5*rr5
      sr7     = bn(3) - sc7*rr7

      ! Intermediates terms involving moments and distance separation
      uri      =   dpui%x*pos%x +  dpui%y*pos%y +  dpui%z*pos%z
      urk      =   dpuk%x*pos%x +  dpuk%y*pos%y +  dpuk%z*pos%z
      urip     =  dpui%xx*pos%x + dpui%yy*pos%y + dpui%zz*pos%z
      urkp     =  dpuk%xx*pos%x + dpuk%yy*pos%y + dpuk%zz*pos%z

      duik     =   ip%dx*dpuk%x  +  ip%dy*dpuk%y + ip%dz*dpuk%z
     &          +  kp%dx*dpui%x  +  kp%dy*dpui%y + kp%dz*dpui%z
      quik     =    qrix*dpuk%x  +   qriy*dpuk%y +  qriz*dpuk%z
     &           -  qrkx*dpui%x  -   qrky*dpui%y -  qrkz*dpui%z
c
c     Calculate intermediate terms for polarization interaction
c
      term1    =  kp%c*uri  - ip%c*urk + duik
      term2    =   two*quik -  uri*drk - dri*urk
      term3    =   uri*qrrk -  urk*qrri
c
c     compute the energy contribution for this interaction
c
      e        = e + term1*sr3 + term2*sr5 + term3*sr7

      qrimodx  =                  qriy*pos%y  +  qriz*pos%z
      qrimody  =  qrix*pos%x   +                 qriz*pos%z
      qrimodz  =  qrix*pos%x   +  qriy*pos%y
      qrkmodx  =                  qrky*pos%y  +  qrkz*pos%z
      qrkmody  =  qrkx*pos%x   +                 qrkz*pos%z
      qrkmodz  =  qrkx*pos%x   +  qrky*pos%y
c
c     t the dEd/dR terms used for direct polarization force
c
      dterm1   =  bn(2)   -   sc3*rr5
      dterm2   =  bn(3)   -   sc5*rr7
c
c     Straight terms ( xx, yy ,zz )
c
      dterm3x  = - sr3  + dterm1*pos%x**2 - rr3*pos%x*rc3x
      dterm3y  = - sr3  + dterm1*pos%y**2 - rr3*pos%y*rc3y
      dterm3z  = - sr3  + dterm1*pos%z**2 - rr3*pos%z*rc3z

      dterm4x  =  rr3*rc3x  - (dterm1 + sr5)*pos%x
      dterm4y  =  rr3*rc3y  - (dterm1 + sr5)*pos%y
      dterm4z  =  rr3*rc3z  - (dterm1 + sr5)*pos%z

      dterm5x  = - sr5      + dterm2*pos%x**2 - rr5*pos%x*rc5x
      dterm5y  = - sr5      + dterm2*pos%y**2 - rr5*pos%y*rc5y
      dterm5z  = - sr5      + dterm2*pos%z**2 - rr5*pos%z*rc5z

      dterm6x  =  (bn(4) - sc7*rr9)*pos%x**2  - bn(3) - rr7*pos%x*rc7x
      dterm6y  =  (bn(4) - sc7*rr9)*pos%y**2  - bn(3) - rr7*pos%y*rc7y
      dterm6z  =  (bn(4) - sc7*rr9)*pos%z**2  - bn(3) - rr7*pos%z*rc7z

      dterm7x  = rr5*rc5x  - two*bn(3)*pos%x
     &           + ( sc5   + 1.5_ti_p*sc7 )*rr7*pos%x
      dterm7y  = rr5*rc5y  - two*bn(3)*pos%y
     &           + ( sc5   + 1.5_ti_p*sc7)*rr7*pos%y
      dterm7z  = rr5*rc5z  - two*bn(3)*pos%z
     &           + ( sc5   + 1.5_ti_p*sc7 )*rr7*pos%z
c
c     Straight terms ( xx, yy ,zz )
c
      tisx     = ip%c*dterm3x      + ip%dx*dterm4x + dri*dterm5x
     &          + two*sr5*ip%qxx  + qrimodx*sc7*rr7
     &          + two*qrix*dterm7x + qrri*dterm6x
      tisy     = ip%c*dterm3y      + ip%dy*dterm4y + dri*dterm5y
     &          + two*sr5*ip%qyy  + qrimody*sc7*rr7
     &          + two*qriy*dterm7y + qrri*dterm6y
      tisz     = ip%c*dterm3z      + ip%dz*dterm4z + dri*dterm5z
     &          + two*sr5*ip%qzz  + qrimodz*sc7*rr7
     &          + two*qriz*dterm7z +  qrri*dterm6z

      tksx     = kp%c*dterm3x      - kp%dx*dterm4x - drk*dterm5x
     &          + two*sr5*kp%qxx  + qrkmodx*sc7*rr7
     &          + two*qrkx*dterm7x + qrrk*dterm6x
      tksy     = kp%c*dterm3y      - kp%dy*dterm4y - drk*dterm5y
     &          + two*sr5*kp%qyy  + qrkmody*sc7*rr7
     &          + two*qrky*dterm7y + qrrk*dterm6y
      tksz     = kp%c*dterm3z      - kp%dz*dterm4z  - drk*dterm5z
     &          + two*sr5*kp%qzz  + qrkmodz*sc7*rr7
     &          + two*qrkz*dterm7z + qrrk*dterm6z
c
c     Crossed terms ( xy = yx , xz = zx ,yz = zy )
c
      tmp1%x   = pos%x * pos%y
      tmp1%y   = pos%x * pos%z
      tmp1%z   = pos%y * pos%z

      dterm3x  =   dterm1*tmp1%x -  rr3*pos%y*rc3x
      dterm3y  =   dterm1*tmp1%y -  rr3*pos%z*rc3x
      dterm3z  =   dterm1*tmp1%z -  rr3*pos%z*rc3y
      dterm4x  = - dterm1*pos%x  +  rr3*rc3x
      dterm4y  = - dterm1*pos%x  +  rr3*rc3x
      dterm4z  = - dterm1*pos%y  +  rr3*rc3y

      dterm5x  =   dterm2*tmp1%x -  rr5*pos%y*rc5x
      dterm5y  =   dterm2*tmp1%y -  rr5*pos%z*rc5x
      dterm5z  =   dterm2*tmp1%z -  rr5*pos%z*rc5y
      dterm6x  =  (bn(4) - sc7*rr9)*tmp1%x - rr7*pos%y*rc7x
      dterm6y  =  (bn(4) - sc7*rr9)*tmp1%y - rr7*pos%z*rc7x
      dterm6z  =  (bn(4) - sc7*rr9)*tmp1%z - rr7*pos%z*rc7y
      dterm7x  = - dterm2*pos%x  +  rr5*rc5x
      dterm7y  = - dterm2*pos%x  +  rr5*rc5x
      dterm7z  = - dterm2*pos%y  +  rr5*rc5y
c
c     Crossed terms ( xy = yx , xz = zx ,yz = zy )
c

      ticx     = ip%c*dterm3x + ip%dy*dterm4x + dri*dterm5x
     &         + qrri*dterm6x - sr5*ip%dx*pos%y
     &         +  two*(sr5*ip%qxy - sr7*pos%y*qrix + qriy*dterm7x)

      ticy     = ip%c*dterm3y + ip%dz*dterm4y + dri*dterm5y
     &         + qrri*dterm6y - sr5*ip%dx*pos%z
     &         +  two*(sr5*ip%qxz - sr7*pos%z*qrix + qriz*dterm7y)

      ticz     = ip%c*dterm3z + ip%dz*dterm4z + dri*dterm5z
     &         + qrri*dterm6z - sr5*ip%dy*pos%z
     &         +  two*(sr5*ip%qyz - sr7*pos%z*qriy + qriz*dterm7z)

      tkcx     = kp%c*dterm3x - kp%dy*dterm4x - drk*dterm5x
     &         + qrrk*dterm6x + sr5*kp%dx*pos%y
     &         +  two*(sr5*kp%qxy - sr7*pos%y*qrkx + qrky*dterm7x)

      tkcy     = kp%c*dterm3y - kp%dz*dterm4y - drk*dterm5y
     &         + qrrk*dterm6y + sr5*kp%dx*pos%z
     &         +  two*(sr5*kp%qxz - sr7*pos%z*qrkx + qrkz*dterm7y)

      tkcz     = kp%c*dterm3z- kp%dz*dterm4z - drk*dterm5z
     &         + qrrk*dterm6z + sr5 *kp%dy*pos%z
     &         +  two*(sr5*kp%qyz - sr7*pos%z*qrky + qrkz*dterm7z)
c
c      Construct matrixes for dot_product
c      do Dot product
c
      depx     =  tisx*dpuk%xx + ticx*dpuk%yy + ticy*dpuk%zz
     &         -  tksx*dpui%xx - tkcx*dpui%yy - tkcy*dpui%zz
      depy     =  ticx*dpuk%xx + tisy*dpuk%yy + ticz*dpuk%zz
     &         -  tkcx*dpui%xx - tksy*dpui%yy - tkcz*dpui%zz
      depz     =  ticy*dpuk%xx + ticz*dpuk%yy + tisz*dpuk%zz
     &         -  tkcy*dpui%xx - tkcz*dpui%yy - tksz*dpui%zz

      frc%x    = frc%x - depx
      frc%y    = frc%y - depy
      frc%z    = frc%z - depz
c
c     t the dEp/dR terms used for direct polarization force
c

      dterm1   =  bn(2) - sc3 * rr5
      dterm2   =  bn(3) - sc5 * rr7
c
c     Straight terms ( xx, yy ,zz )
c
      dterm3x  = - sr3 + dterm1*pos%x**2 - rr3*pos%x*rc3x
      dterm3y  = - sr3 + dterm1*pos%y**2 - rr3*pos%y*rc3y
      dterm3z  = - sr3 + dterm1*pos%z**2 - rr3*pos%z*rc3z

      dterm4x  = - sr5*pos%x - dterm1*pos%x + rr3*rc3x
      dterm4y  = - sr5*pos%y - dterm1*pos%y + rr3*rc3y
      dterm4z  = - sr5*pos%z - dterm1*pos%z + rr3*rc3z

      dterm5x  = - sr5 + dterm2*pos%x**2 - rr5*pos%x*rc5x
      dterm5y  = - sr5 + dterm2*pos%y**2 - rr5*pos%y*rc5y
      dterm5z  = - sr5 + dterm2*pos%z**2 - rr5*pos%z*rc5z

      dterm6x  =  (bn(4) - sc7*rr9)*pos%x**2 - bn(3) - rr7*pos%x*rc7x
      dterm6y  =  (bn(4) - sc7*rr9)*pos%y**2 - bn(3) - rr7*pos%y*rc7y
      dterm6z  =  (bn(4) - sc7*rr9)*pos%z**2 - bn(3) - rr7*pos%z*rc7z

      dterm7x  =  rr5*rc5x - two*bn(3)*pos%x
     &         +  (sc5 + 1.5_ti_p*sc7)*rr7*pos%x
      dterm7y  =  rr5*rc5y - two*bn(3)*pos%y
     &         +  (sc5 + 1.5_ti_p*sc7)*rr7*pos%y
      dterm7z  =  rr5*rc5z - two*bn(3)*pos%z
     &         +  (sc5 + 1.5_ti_p*sc7)*rr7*pos%z
c
c     Straight terms ( xx, yy ,zz )
c
      tisx = ip%c*dterm3x + ip%dx*dterm4x + dri*dterm5x
     &     + qrri*dterm6x + qrimodx*sc7*rr7
     &     + two*(sr5*ip%qxx + qrix*dterm7x)

      tisy = ip%c*dterm3y + ip%dy*dterm4y + dri*dterm5y
     &     + qrri*dterm6y + qrimody*sc7*rr7
     &     + two*(sr5*ip%qyy + qriy*dterm7y)

      tisz = ip%c*dterm3z + ip%dz*dterm4z + dri*dterm5z
     &     + qrri*dterm6z + qrimodz*sc7*rr7
     &     + two*(sr5*ip%qzz + qriz*dterm7z)

      tksx = kp%c*dterm3x - kp%dx*dterm4x - drk*dterm5x
     &     + qrrk*dterm6x + qrkmodx*sc7*rr7
     &     + two*(sr5*kp%qxx + qrkx*dterm7x)

      tksy = kp%c*dterm3y - kp%dy*dterm4y - drk*dterm5y
     &     + qrrk*dterm6y + qrkmody*sc7*rr7
     &     + two*(sr5*kp%qyy + qrky*dterm7y)

      tksz = kp%c*dterm3z - kp%dz*dterm4z - drk*dterm5z
     &     + qrrk*dterm6z + qrkmodz*sc7*rr7
     &     + two*(sr5*kp%qzz + qrkz*dterm7z)
c
c     Crossed terms ( xy = yx , xz = zx ,yz = zy )
c
      tmp1%x  = pos%x * pos%y
      tmp1%y  = pos%x * pos%z
      tmp1%z  = pos%y * pos%z

      dterm3x =  dterm1*tmp1%x - rr3*pos%y*rc3x
      dterm3y =  dterm1*tmp1%y - rr3*pos%z*rc3x
      dterm3z =  dterm1*tmp1%z - rr3*pos%z*rc3y

      dterm4x =  rr3*rc3x - dterm1*pos%x
      dterm4y =  rr3*rc3x - dterm1*pos%x
      dterm4z =  rr3*rc3y - dterm1*pos%y
      dterm5x =  dterm2*tmp1%x - rr5*pos%y*rc5x
      dterm5y =  dterm2*tmp1%y - rr5*pos%z*rc5x
      dterm5z =  dterm2*tmp1%z - rr5*pos%z*rc5y

      dterm6x =  (bn(4) - sc7*rr9)*tmp1%x - rr7*pos%y*rc7x
      dterm6y =  (bn(4) - sc7*rr9)*tmp1%y - rr7*pos%z*rc7x
      dterm6z =  (bn(4) - sc7*rr9)*tmp1%z - rr7*pos%z*rc7y
      dterm7x =  rr5*rc5x - dterm2*pos%x
      dterm7y =  rr5*rc5x - dterm2*pos%x
      dterm7z =  rr5*rc5y - dterm2*pos%y

c
c     Crossed terms ( xy = yx , xz = zx ,yz = zy )
c
      ticx = ip%c*dterm3x + ip%dy*dterm4x + dri*dterm5x
     &     + qrri*dterm6x - sr5*ip%dx*pos%y
     &     + two*(sr5*ip%qxy - sr7*pos%y*qrix + qriy*dterm7x)

      ticy = ip%c*dterm3y + ip%dz*dterm4y + dri*dterm5y
     &     + qrri*dterm6y - sr5*ip%dx*pos%z
     &     + two*(sr5*ip%qxz - sr7*pos%z*qrix + qriz*dterm7y)

      ticz = ip%c*dterm3z + ip%dz*dterm4z + dri*dterm5z
     &     + qrri*dterm6z - sr5*ip%dy*pos%z
     &     +  two*(sr5*ip%qyz - sr7*pos%z*qriy + qriz*dterm7z)

      tkcx = kp%c*dterm3x - kp%dy*dterm4x - drk*dterm5x
     &     + qrrk*dterm6x + sr5*kp%dx*pos%y
     &     +  two*(sr5*kp%qxy - sr7*pos%y*qrkx + qrky*dterm7x)

      tkcy = kp%c*dterm3y - kp%dz*dterm4y - drk*dterm5y
     &     + qrrk*dterm6y + sr5*kp%dx*pos%z
     &     + two*(sr5*kp%qxz - sr7*pos%z*qrkx + qrkz*dterm7y)

      tkcz = kp%c*dterm3z - kp%dz*dterm4z - drk*dterm5z
     &     + qrrk*dterm6z + sr5*kp%dy*pos%z
     &     + two*(sr5*kp%qyz - sr7*pos%z*qrky + qrkz*dterm7z)
c
c     Construct matrixes for dot_product
c     Do dot product
c
      depx = tisx*dpuk%x + ticx*dpuk%y + ticy*dpuk%z
     &     - tksx*dpui%x - tkcx*dpui%y - tkcy*dpui%z
      depy = ticx*dpuk%x + tisy*dpuk%y + ticz*dpuk%z
     &     - tkcx*dpui%x - tksy*dpui%y - tkcz*dpui%z
      depz = ticy*dpuk%x + ticz*dpuk%y + tisz*dpuk%z
     &     - tkcy*dpui%x - tkcz*dpui%y - tksz*dpui%z

      frc%x = frc%x - depx
      frc%y = frc%y - depy
      frc%z = frc%z - depz

c
c     t the dtau/dr terms used for mutual polarization force
c
      dterm1 = bn(2) - sc3 * rr5
      dterm2 = bn(3) - sc5 * rr7
c
c     Straight terms ( xx, yy ,zz )
c
      dterm3x = sr5 + dterm1
      dterm3y = sr5 + dterm1
      dterm3z = sr5 + dterm1
      dterm4x = rr3  * uscale
      dterm4y = rr3  * uscale
      dterm4z = rr3  * uscale

      dterm5x = - pos%x*dterm3x + rc3x*dterm4x
      dterm5y = - pos%y*dterm3y + rc3y*dterm4y
      dterm5z = - pos%z*dterm3z + rc3z*dterm4z

      dterm6x = - sr5 + pos%x**2*dterm2 - rr5*pos%x*rc5x
      dterm6y = - sr5 + pos%y**2*dterm2 - rr5*pos%y*rc5y
      dterm6z = - sr5 + pos%z**2*dterm2 - rr5*pos%z*rc5z

      tisx =  dpui%x*dterm5x + uri*dterm6x
      tisy =  dpui%y*dterm5y + uri*dterm6y
      tisz =  dpui%z*dterm5z + uri*dterm6z
      tksx =  dpuk%x*dterm5x + urk*dterm6x
      tksy =  dpuk%y*dterm5y + urk*dterm6y
      tksz =  dpuk%z*dterm5z + urk*dterm6z
c
c     Crossed terms ( xy = yx , xz = zx ,yz = zy )
c
      dterm4x = - sr5*pos%y
      dterm4y = - sr5*pos%z
      dterm4z = - sr5*pos%z
      dterm5x = - pos%x*dterm1 + rr3*rc3x
      dterm5y = - pos%x*dterm1 + rr3*rc3x
      dterm5z = - pos%y*dterm1 + rr3*rc3y
      dterm6x =   pos%x*dterm2*pos%y - rr5*pos%y*rc5x
      dterm6y =   pos%x*dterm2*pos%z - rr5*pos%z*rc5x
      dterm6z =   pos%y*dterm2*pos%z - rr5*pos%z*rc5y

      ticx =  dpui%x*dterm4x + dpui%y*dterm5x + uri*dterm6x
      ticy =  dpui%x*dterm4y + dpui%z*dterm5y + uri*dterm6y
      ticz =  dpui%y*dterm4z + dpui%z*dterm5z + uri*dterm6z

      tkcx =  dpuk%x*dterm4x + dpuk%y*dterm5x + urk*dterm6x
      tkcy =  dpuk%x*dterm4y + dpuk%z*dterm5y + urk*dterm6y
      tkcz =  dpuk%y*dterm4z + dpuk%z*dterm5z + urk*dterm6z
c
c     Construct matrixes for dot_product
c     Dot product
c
      depx =  tisx*dpuk%xx + ticx*dpuk%yy + ticy*dpuk%zz
     &      + tksx*dpui%xx + tkcx*dpui%yy + tkcy*dpui%zz
      depy =  ticx*dpuk%xx + tisy*dpuk%yy + ticz*dpuk%zz
     &      + tkcx*dpui%xx + tksy*dpui%yy + tkcz*dpui%zz
      depz =  ticy*dpuk%xx + ticz*dpuk%yy + tisz*dpuk%zz
     &      + tkcy*dpui%xx + tkcz*dpui%yy + tksz*dpui%zz

      frc%x   = frc%x - depx
      frc%y   = frc%y - depy
      frc%z   = frc%z - depz
      frc_r%x = WRITE_C(frc_r%x -) tp2mdr(frc%x)
      frc_r%y = WRITE_C(frc_r%y -) tp2mdr(frc%y)
      frc_r%z = WRITE_C(frc_r%z -) tp2mdr(frc%z)
c
c     t the induced dipole field used for dipole torques
c
      turi5    = -sr5*urk - sr5*urkp
      turk5    = -sr5*uri - sr5*urip
      turi7    = - sr7*urk - sr7*urkp
      turk7    = - sr7*uri - sr7*urip
c
c     t induced dipole field gradient used for quadrupole torques
c
      ti5x     =  two*(sr5*dpuk%x + sr5*dpuk%xx)
      ti5y     =  two*(sr5*dpuk%y + sr5*dpuk%yy)
      ti5z     =  two*(sr5*dpuk%z + sr5*dpuk%zz)
      tk5x     =  two*(sr5*dpui%x + sr5*dpui%xx)
      tk5y     =  two*(sr5*dpui%y + sr5*dpui%yy)
      tk5z     =  two*(sr5*dpui%z + sr5*dpui%zz)
c
c     Torque is induced field and gradient cross permanent moments
c
      ufli%x   = sr3*dpuk%x + sr3*dpuk%xx + pos%x*turi5
      ufli%y   = sr3*dpuk%y + sr3*dpuk%yy + pos%y*turi5
      ufli%z   = sr3*dpuk%z + sr3*dpuk%zz + pos%z*turi5
      uflk%x   = sr3*dpui%x + sr3*dpui%xx + pos%x*turk5
      uflk%y   = sr3*dpui%y + sr3*dpui%yy + pos%y*turk5
      uflk%z   = sr3*dpui%z + sr3*dpui%zz + pos%z*turk5

      dufli%x  = pos%x*ti5x + pos%x**2*turi7
      dufli%z  = pos%y*ti5y + pos%y**2*turi7
      dufli%zz = pos%z*ti5z + pos%z**2*turi7
      dufli%y  = pos%x*ti5y + pos%y*ti5x + two*pos%x*pos%y*turi7
      dufli%xx = pos%x*ti5z + pos%z*ti5x + two*pos%x*pos%z*turi7
      dufli%yy = pos%y*ti5z + pos%z*ti5y + two*pos%y*pos%z*turi7

      duflk%x  = -(pos%x*tk5x + pos%x**2*turk7)
      duflk%z  = -(pos%y*tk5y + pos%y**2*turk7)
      duflk%zz = -(pos%z*tk5z + pos%z**2*turk7)
      duflk%y  = -(pos%x*tk5y + pos%y*tk5x + two*pos%x*pos%y*turk7)
      duflk%xx = -(pos%x*tk5z + pos%z*tk5x + two*pos%x*pos%z*turk7)
      duflk%yy = -(pos%y*tk5z + pos%z*tk5y + two*pos%y*pos%z*turk7)

      trqi%x = trqi%x + ip%dz*ufli%y - ip%dy*ufli%z
     &       + ip%qxz*dufli%y - ip%qxy*dufli%xx
     &       + (ip%qzz - ip%qyy)*dufli%yy
     &       + two*ip%qyz*( dufli%z - dufli%zz )
      trqi%y = trqi%y + ip%dx*ufli%z - ip%dz*ufli%x
     &       - ip%qyz*dufli%y
     &       + ip%qxy*dufli%yy + (ip%qxx - ip%qzz)*dufli%xx
     &       + two*ip%qxz*( dufli%zz - dufli%x )
      trqi%z = trqi%z + ip%dy*ufli%x - ip%dx*ufli%y
     &       + ip%qyz*dufli%xx
     &       - ip%qxz*dufli%yy + (ip%qyy - ip%qxx)*dufli%y
     &       + two*ip%qxy*( dufli%x - dufli%z )

      trqk%x = trqk%x + kp%dz*uflk%y - kp%dy*uflk%z
     &       + kp%qxz*duflk%y
     &       - kp%qxy*duflk%xx + (kp%qzz - kp%qyy)*duflk%yy
     &       + two*kp%qyz*( duflk%z - duflk%zz )
      trqk%y = trqk%y + kp%dx*uflk%z - kp%dz*uflk%x
     &       - kp%qyz*duflk%y
     &       + kp%qxy*duflk%yy + (kp%qxx - kp%qzz)*duflk%xx
     &       + two*kp%qxz*( duflk%zz - duflk%x )
      trqk%z = trqk%z + kp%dy*uflk%x - kp%dx*uflk%y
     &       + kp%qyz*duflk%xx
     &       - kp%qxz*duflk%yy + (kp%qyy - kp%qxx)*duflk%y
     &       + two*kp%qxy*( duflk%x - duflk%z )
      end subroutine
#endif
