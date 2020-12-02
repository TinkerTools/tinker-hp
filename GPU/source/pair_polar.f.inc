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
      use utilcu ,only: f_erfc
#  if defined(SINGLE)||defined(MIXED)
      use utilcu ,only: f_sqrt,f_exp
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
      real(t_p) t1,t2,t3,t4,t5,t6
      real(t_p) ck,dkx,dky,dkz
      real(t_p) r,invr2,invr
      real(t_p) bn0,bn1,bn2,bn3,bn4
      real(t_p) rr3,rr5,rr7,rr9
      real(t_p) damp1,invdamp
      real(t_p) expdamp,expdamp1
      real(t_p) da
      real(t_p) qk1,qk2,qk3,qk4,qk5,qk6,qk7,qk8,qk9
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
      real(t_p) tmp2x,tmp2y,tmp2z
      real(t_p) tmp3x,tmp3y,tmp3z
      real(t_p) tmp4x,tmp4y,tmp4z
      real(t_p) tmp5x,tmp5y,tmp5z
      real(t_p) tmp6x,tmp6y,tmp6z
      real(t_p) tmp7x,tmp7y,tmp7z
      real(t_p) tisx,tisy,tisz
      real(t_p) tikx,tiky,tikz,ticx,ticy,ticz
      real(t_p) tkcx,tkcy,tkcz,tksx,tksy,tksz
      real(t_p) ti1,ti2,ti3,ti4,ti5,ti6,ti7,ti8,ti9
      real(t_p) tk1,tk2,tk3,tk4,tk5,tk6,tk7,tk8,tk9
      real(t_p) ti3x,ti3y,ti3z,tk3x,tk3y,tk3z
      real(t_p) ti5x,ti5y,ti5z,tk5x,tk5y,tk5z
      real(t_p) turi5,turi7,turk5,turk7
      real(t_p) depx,depy,depz
      real(t_p) dipolx,dipoly,dipolz
      real(t_p) qipol1,qipol2,qipol3
      real(t_p) qipol4,qipol5,qipol6
      real(t_p) qipol7,qipol8,qipol9
      type(real6) :: dufli,duflk
      type(real3) :: ufli,uflk
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
      bn0      = 0.0_ti_p
      bn1      = 0.0_ti_p
      bn2      = 0.0_ti_p
      bn3      = 0.0_ti_p
      bn4      = 0.0_ti_p
      else
      ralpha   = aewald * r
      !call erfcore_inl(ralpha, bn0,1)
      exp2a    = f_exp( - ralpha**2)
      bn0      = f_erfc(ralpha)

      bn0      = bn0 * invr
      bn1      = ( 1.0_ti_p*bn0 + alsq2    *alsq2n*exp2a ) * invr2
      bn2      = ( 3.0_ti_p*bn1 + alsq2**2 *alsq2n*exp2a ) * invr2
      bn3      = ( 5.0_ti_p*bn2 + alsq2**2*alsq2    *alsq2n*exp2a )
     &           * invr2
      bn4      = ( 7.0_ti_p*bn3 + alsq2**2*alsq2**2 *alsq2n*exp2a )
     &           * invr2

      bn0      = f * bn0
      bn1      = f * bn1
      bn2      = f * bn2
      bn3      = f * bn3
      bn4      = f * bn4
      end if

      rr3      = f        * invr * invr2
      rr5      = 3.0_ti_p * f*invr*invr2  * invr2
      rr7      = 5.0_ti_p * rr5  * invr2
      rr9      = 7.0_ti_p * rr7  * invr2
c
c    Apply Thole polarization damping to scale factors
c
      if (damp.ne.0.0_ti_p) then
         invdamp = damp**(-one)
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
         psc3     = sc3 * pscale
         dsc3     = sc3 * dscale
         usc3     = sc3 * uscale
         psc5     = sc5 * pscale
         dsc5     = sc5 * dscale
         usc5     = sc5 * uscale
         psc7     = sc7 * pscale
         dsc7     = sc7 * dscale

         psr7     = bn3 - psc7 * rr7
         dsr7     = bn3 - dsc7 * rr7
         psr3     = bn1 - psc3 * rr3
         dsr3     = bn1 - dsc3 * rr3
         usr3     = bn1 - usc3 * rr3
         psr5     = bn2 - psc5 * rr5
         dsr5     = bn2 - dsc5 * rr5
         usr5     = bn2 - usc5 * rr5

         prc3x    = - rc3x * pscale
         drc3x    = - rc3x * dscale
         urc3x    = - rc3x * uscale
         prc3y    = - rc3y * pscale
         drc3y    = - rc3y * dscale
         urc3y    = - rc3y * uscale
         prc3z    = - rc3z * pscale
         drc3z    = - rc3z * dscale
         urc3z    = - rc3z * uscale

         prc5x    = - rc5x * pscale
         drc5x    = - rc5x * dscale
         urc5x    = - rc5x * uscale
         prc5y    = - rc5y * pscale
         drc5y    = - rc5y * dscale
         urc5y    = - rc5y * uscale
         prc5z    = - rc5z * pscale
         drc5z    = - rc5z * dscale
         urc5z    = - rc5z * uscale

         prc7x    = - rc7x * pscale
         drc7x    = - rc7x * dscale
         prc7y    = - rc7y * pscale
         drc7y    = - rc7y * dscale
         prc7z    = - rc7z * pscale
         drc7z    = - rc7z * dscale
      else
         psc3     = 1.0_ti_p - sc3 * pscale
         dsc3     = 1.0_ti_p - sc3 * dscale
         usc3     = 1.0_ti_p - sc3 * uscale
         psc5     = 1.0_ti_p - sc5 * pscale
         dsc5     = 1.0_ti_p - sc5 * dscale
         usc5     = 1.0_ti_p - sc5 * uscale
         psc7     = 1.0_ti_p - sc7 * pscale
         dsc7     = 1.0_ti_p - sc7 * dscale

         psr3     = bn1 - psc3 * rr3
         dsr3     = bn1 - dsc3 * rr3
         usr3     = bn1 - usc3 * rr3
         psr5     = bn2 - psc5 * rr5
         dsr5     = bn2 - dsc5 * rr5
         usr5     = bn2 - usc5 * rr5
         psr7     = bn3 - psc7 * rr7
         dsr7     = bn3 - dsc7 * rr7

         prc3x    = rc3x * pscale
         drc3x    = rc3x * dscale
         urc3x    = rc3x * uscale
         prc3y    = rc3y * pscale
         drc3y    = rc3y * dscale
         urc3y    = rc3y * uscale
         prc3z    = rc3z * pscale
         drc3z    = rc3z * dscale
         urc3z    = rc3z * uscale

         prc5x    = rc5x * pscale
         drc5x    = rc5x * dscale
         urc5x    = rc5x * uscale
         prc5y    = rc5y * pscale
         drc5y    = rc5y * dscale
         urc5y    = rc5y * uscale
         prc5z    = rc5z * pscale
         drc5z    = rc5z * dscale
         urc5z    = rc5z * uscale

         prc7x    = rc7x * pscale
         drc7x    = rc7x * dscale
         prc7y    = rc7y * pscale
         drc7y    = rc7y * dscale
         prc7z    = rc7z * pscale
         drc7z    = rc7z * dscale
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
      urip     =  dpui%xx*pos%x + dpui%yy*pos%y + dpui%zz*pos%z
      urkp     =  dpuk%xx*pos%x + dpuk%yy*pos%y + dpuk%zz*pos%z

      duik     =    ip%dx*dpuk%x  +  ip%dy*dpuk%y  +  ip%dz*dpuk%z
     &          +  dpui%x*kp%dx  + dpui%y*kp%dy  + dpui%z*kp%dz
      quik     =     qrix*dpuk%x  +   qriy*dpuk%y  +   qriz*dpuk%z
     &          -  dpui%x*qrkx - dpui%y*qrky - dpui%z*qrkz
c
c     Calculate intermediate terms for polarization interaction
c
      term1    =         kp%c*uri  -      ip%c*urk  +      duik
      term2    =   two*quik -       uri*drk  -   dri*urk
      term3    =        uri*qrrk -       urk*qrri
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
      dterm1   =  bn2    -   dsc3*rr5
      dterm2   =  bn3    -   dsc5*rr7
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

      dterm6x  =  (bn4 - dsc7*rr9)*pos%x**2  - bn3 - rr7*pos%x*drc7x
      dterm6y  =  (bn4 - dsc7*rr9)*pos%y**2  - bn3 - rr7*pos%y*drc7y
      dterm6z  =  (bn4 - dsc7*rr9)*pos%z**2  - bn3 - rr7*pos%z*drc7z

      dterm7x  = rr5*drc5x  - two*bn3*pos%x
     &           + ( dsc5   + 1.5_ti_p*dsc7 )*rr7*pos%x
      dterm7y  = rr5*drc5y  - two*bn3*pos%y
     &           + ( dsc5   + 1.5_ti_p*dsc7)*rr7*pos%y
      dterm7z  = rr5*drc5z  - two*bn3*pos%z
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

      tmp2x    =        drc3x
      tmp2y    =        drc3x
      tmp2z    =        drc3y

      tmp3x    = pos%y * drc3x
      tmp3y    = pos%z * drc3x
      tmp3z    = pos%z * drc3y
      tmp4x    = pos%x
      tmp4y    = pos%x
      tmp4z    = pos%y

      tmp5x    = pos%y * drc5x
      tmp5y    = pos%z * drc5x
      tmp5z    = pos%z * drc5y

      tmp6x    = pos%y * drc7x
      tmp6y    = pos%z * drc7x
      tmp6z    = pos%z * drc7y

      tmp7x    =        drc5x
      tmp7y    =        drc5x
      tmp7z    =        drc5y

      dterm3x  =   dterm1*tmp1x  -  rr3*tmp3x
      dterm3y  =   dterm1*tmp1y  -  rr3*tmp3y
      dterm3z  =   dterm1*tmp1z  -  rr3*tmp3z
      dterm4x  = - dterm1*tmp4x  +  rr3*tmp2x
      dterm4y  = - dterm1*tmp4y  +  rr3*tmp2y
      dterm4z  = - dterm1*tmp4z  +  rr3*tmp2z

      dterm5x  =   dterm2*tmp1x  -  rr5*tmp5x
      dterm5y  =   dterm2*tmp1y  -  rr5*tmp5y
      dterm5z  =   dterm2*tmp1z  -  rr5*tmp5z
      dterm6x  =  (bn4           - dsc7*rr9)*tmp1x - rr7*tmp6x
      dterm6y  =  (bn4           - dsc7*rr9)*tmp1y - rr7*tmp6y
      dterm6z  =  (bn4           - dsc7*rr9)*tmp1z - rr7*tmp6z
      dterm7x  = - dterm2*tmp4x  +  rr5*tmp7x
      dterm7y  = - dterm2*tmp4y  +  rr5*tmp7y
      dterm7z  = - dterm2*tmp4z  +  rr5*tmp7z
c
c     Crossed terms ( xy = yx , xz = zx ,yz = zy )
c

      ticx     =  ip%c*dterm3x  -  dsr5*ip%dx*pos%y
     &          + ip%dy*dterm4x +  dri*dterm5x
     &          + two*(dsr5*ip%qxy - dsr7*pos%y*qrix + qriy*dterm7x)
     &          + qrri*dterm6x

      ticy     =  ip%c*dterm3y  -  dsr5*ip%dx*pos%z
     &          + ip%dz*dterm4y +         dri*dterm5y
     &          + two*(dsr5*ip%qxz - dsr7*pos%z*qrix + qriz*dterm7y)
     &          + qrri*dterm6y

      ticz     =  ip%c*dterm3z  -  dsr5*ip%dy*pos%z
     &          + ip%dz*dterm4z +  dri*dterm5z
     &          + two* (dsr5*ip%qyz - dsr7*pos%z*qriy + qriz*dterm7z)
     &          + qrri*dterm6z

      tkcx     =  kp%c*dterm3x  +  dsr5 *kp%dx*pos%y
     &          - kp%dy*dterm4x -  drk*dterm5x
     &          + two* (dsr5*kp%qxy - dsr7*pos%y*qrkx + qrky*dterm7x)
     &          + qrrk*dterm6x

      tkcy     =  kp%c*dterm3y +  dsr5 *kp%dx*pos%z
     &          - kp%dz*dterm4y - drk*dterm5y
     &          + two* (dsr5*kp%qxz - dsr7*pos%z*qrkx + qrkz*dterm7y)
     &          + qrrk*dterm6y

      tkcz     =  kp%c*dterm3z +  dsr5 *kp%dy*pos%z
     &          - kp%dz*dterm4z - drk*dterm5z
     &          + two* (dsr5*kp%qyz - dsr7*pos%z*qrky + qrkz*dterm7z)
     &          + qrrk*dterm6z
c
c      Construct matrixes for dot_product
c
      ti1      = tisx !xx
      ti5      = tisy !yy
      ti9      = tisz !zz
      ti2      = ticx !xy
      ti4      = ticx !yx
      ti3      = ticy !xz
      ti7      = ticy !zx
      ti6      = ticz !yz
      ti8      = ticz !zy
      tk1      = tksx !xx
      tk5      = tksy !yy
      tk9      = tksz !zz
      tk2      = tkcx !xy
      tk4      = tkcx !yx
      tk3      = tkcy !xz
      tk7      = tkcy !zx
      tk6      = tkcz !yz
      tk8      = tkcz !zy
c
c     do dot product
c
      frc%x     =  ti1 *dpuk%xx + ti2 *dpuk%yy + ti3 *dpuk%zz
     &          -  tk1 *dpui%xx - tk2 *dpui%yy - tk3 *dpui%zz
      frc%y     =  ti4 *dpuk%xx + ti5 *dpuk%yy + ti6 *dpuk%zz
     &          -  tk4 *dpui%xx - tk5 *dpui%yy - tk6 *dpui%zz
      frc%z     =  ti7 *dpuk%xx + ti8 *dpuk%yy + ti9 *dpuk%zz
     &          -  tk7 *dpui%xx - tk8 *dpui%yy - tk9 *dpui%zz
c
c     t the dEp/dR terms used for direct polarization force
c

      dterm1   =  bn2 - psc3 * rr5
      dterm2   =  bn3 - psc5 * rr7
c
c     Straight terms ( xx, yy ,zz )
c
      dterm3x  = - psr3
     &              + dterm1*pos%x ** 2
     &              - rr3*pos%x*prc3x
      dterm3y  = - psr3
     &              + dterm1*pos%y ** 2
     &              - rr3*pos%y*prc3y
      dterm3z  = - psr3
     &              + dterm1*pos%z ** 2
     &              - rr3*pos%z*prc3z

      dterm4x  =   rr3*prc3x
     &              - dterm1*pos%x
     &              - psr5*pos%x
      dterm4y  =   rr3*prc3y
     &              - dterm1*pos%y
     &              - psr5*pos%y
      dterm4z  =   rr3*prc3z
     &              - dterm1*pos%z
     &              - psr5*pos%z

      dterm5x  = - psr5 + dterm2*pos%x ** 2 - rr5*pos%x*prc5x
      dterm5y  = - psr5 + dterm2*pos%y ** 2 - rr5*pos%y*prc5y
      dterm5z  = - psr5 + dterm2*pos%z ** 2 - rr5*pos%z*prc5z

      dterm6x  =  (bn4 - psc7*rr9)*pos%x ** 2 - bn3 - rr7*pos%x*prc7x
      dterm6y  =  (bn4 - psc7*rr9)*pos%y ** 2 - bn3 - rr7*pos%y*prc7y
      dterm6z  =  (bn4 - psc7*rr9)*pos%z ** 2 - bn3 - rr7*pos%z*prc7z

      dterm7x  =  rr5* prc5x - two*bn3*pos%x
     &             +  (psc5 + 1.5_ti_p*psc7)*rr7*pos%x
      dterm7y  =  rr5* prc5y - two*bn3*pos%y
     &             +  (psc5 + 1.5_ti_p*psc7)*rr7*pos%y
      dterm7z  =  rr5* prc5z - two*bn3*pos%z
     &             +  (psc5 + 1.5_ti_p*psc7)*rr7*pos%z
c
c     Straight terms ( xx, yy ,zz )
c
      tisx =          ip%c            * dterm3x
     &            +         ip%dx * dterm4x
     &            +         dri     * dterm5x
     &            + two * psr5    * ip%qxx
     &            +         qrimodx * psc7
     &                                    * rr7
     &            + two * qrix    * dterm7x
     &            +         qrri    * dterm6x

      tisy =          ip%c            * dterm3y
     &            +         ip%dy * dterm4y
     &            +         dri     * dterm5y
     &            + two * psr5    * ip%qyy
     &            +         qrimody * psc7
     &                                    * rr7
     &            + two * qriy    * dterm7y
     &            +         qrri    * dterm6y

      tisz =          ip%c            * dterm3z
     &            +         ip%dz * dterm4z
     &            +         dri     * dterm5z
     &            + two * psr5    * ip%qzz
     &            +         qrimodz * psc7
     &                                    * rr7
     &            + two * qriz    * dterm7z
     &            +         qrri    * dterm6z
      tksx =          kp%c            * dterm3x
     &            -         kp%dx     * dterm4x
     &            -         drk     * dterm5x
     &            + two * psr5    * kp%qxx
     &            +         qrkmodx * psc7
     &                                    * rr7
     &            + two * qrkx    * dterm7x
     &            +         qrrk    * dterm6x

      tksy =          kp%c            * dterm3y
     &            -         kp%dy     * dterm4y
     &            -         drk     * dterm5y
     &            + two * psr5    * kp%qyy
     &            +         qrkmody * psc7
     &                                    * rr7
     &            + two * qrky    * dterm7y
     &            +         qrrk    * dterm6y

      tksz =          kp%c            * dterm3z
     &            -         kp%dz     * dterm4z
     &            -         drk     * dterm5z
     &            + two * psr5    * kp%qzz
     &            +         qrkmodz * psc7
     &                                    * rr7
     &            + two * qrkz    * dterm7z
     &            +         qrrk    * dterm6z
c
c     Crossed terms ( xy = yx , xz = zx ,yz = zy )
c
      tmp1x = pos%x * pos%y
      tmp1y = pos%x * pos%z
      tmp1z = pos%y * pos%z

      tmp2x = prc3x
      tmp2y = prc3x
      tmp2z = prc3y

      tmp3x = pos%y * prc3x
      tmp3y = pos%z * prc3x
      tmp3z = pos%z * prc3y
      tmp4x = pos%x
      tmp4y = pos%x
      tmp4z = pos%y

      tmp5x = pos%y * prc5x
      tmp5y = pos%z * prc5x
      tmp5z = pos%z * prc5y
      tmp6x = pos%y * prc7x
      tmp6y = pos%z * prc7x
      tmp6z = pos%z * prc7y

      tmp7x = prc5x
      tmp7y = prc5x
      tmp7z = prc5y

      dterm3x =  dterm1  * tmp1x
     &            - rr3  * tmp3x
      dterm3y =  dterm1  * tmp1y
     &            - rr3  * tmp3y
      dterm3z =  dterm1  * tmp1z
     &            - rr3  * tmp3z

      dterm4x =  rr3  * tmp2x
     &            - dterm1  * tmp4x
      dterm4y =  rr3  * tmp2y
     &            - dterm1  * tmp4y
      dterm4z =  rr3  * tmp2z
     &            - dterm1  * tmp4z
      dterm5x =  dterm2  * tmp1x
     &            - rr5  * tmp5x
      dterm5y =  dterm2  * tmp1y
     &            - rr5  * tmp5y
      dterm5z =  dterm2  * tmp1z
     &            - rr5  * tmp5z

      dterm6x =  (bn4 - psc7  * rr9)
     &                         * tmp1x
     &            -  rr7 * tmp6x
      dterm6y =  (bn4 - psc7  * rr9)
     &                         * tmp1y
     &            -  rr7 * tmp6y
      dterm6z =  (bn4 - psc7  * rr9)
     &                         * tmp1z
     &            -  rr7 * tmp6z
      dterm7x =  rr5  * tmp7x
     &            - dterm2  * tmp4x
      dterm7y =  rr5  * tmp7y
     &            - dterm2  * tmp4y
      dterm7z =  rr5  * tmp7z
     &            - dterm2  * tmp4z

c
c     Crossed terms ( xy = yx , xz = zx ,yz = zy )
c

      ticx =          ip%c         * dterm3x
     &            -         psr5 * ip%dx
     &                                 * pos%y
     &            +         ip%dy * dterm4x
     &            +         dri  * dterm5x
     &            + two * psr5 * ip%qxy
     &            - two * psr7 * pos%y
     &                                 * qrix
     &            + two * qriy * dterm7x
     &            +         qrri * dterm6x

      ticy =          ip%c         * dterm3y
     &            -         psr5 * ip%dx
     &                                 * pos%z
     &            +         ip%dz * dterm4y
     &            +         dri  * dterm5y
     &            + two * psr5 * ip%qxz
     &            - two * psr7 * pos%z
     &                                 * qrix
     &            + two * qriz * dterm7y
     &            +         qrri * dterm6y

      ticz =          ip%c         * dterm3z
     &            -         psr5 * ip%dy
     &                                 * pos%z
     &            +         ip%dz * dterm4z
     &            +         dri  * dterm5z
     &            + two * psr5 * ip%qyz
     &            - two * psr7 * pos%z
     &                                 * qriy
     &            + two * qriz * dterm7z
     &            +         qrri * dterm6z
      tkcx =          kp%c         * dterm3x
     &            +         psr5 * kp%dx
     &                                 * pos%y
     &            -         kp%dy  * dterm4x
     &            -         drk  * dterm5x
     &            + two * psr5 * kp%qxy
     &            - two * psr7 * pos%y
     &                                 * qrkx
     &            + two * qrky * dterm7x
     &            +         qrrk * dterm6x

      tkcy =          kp%c         * dterm3y
     &            +         psr5 * kp%dx
     &                                 * pos%z
     &            -         kp%dz  * dterm4y
     &            -         drk  * dterm5y
     &            + two * psr5 * kp%qxz
     &            - two * psr7 * pos%z
     &                                 * qrkx
     &            + two * qrkz * dterm7y
     &            +         qrrk * dterm6y

      tkcz =          kp%c         * dterm3z
     &            +         psr5 * kp%dy
     &                                 * pos%z
     &            -         kp%dz  * dterm4z
     &            -         drk  * dterm5z
     &            + two * psr5 * kp%qyz
     &            - two * psr7 * pos%z
     &                                 * qrky
     &            + two * qrkz * dterm7z
     &            +         qrrk * dterm6z
c
c      Construct matrixes for dot_product
c
      ti1 = tisx !xx
      ti5 = tisy !yy
      ti9 = tisz !zz
      ti2 = ticx !xy
      ti4 = ticx !yx
      ti3 = ticy !xz
      ti7 = ticy !zx
      ti6 = ticz !yz
      ti8 = ticz !zy
      tk1 = tksx !xx
      tk5 = tksy !yy
      tk9 = tksz !zz
      tk2 = tkcx !xy
      tk4 = tkcx !yx
      tk3 = tkcy !xz
      tk7 = tkcy !zx
      tk6 = tkcz !yz
      tk8 = tkcz !zy
c
c     Do dot product
c

      depx =  ti1 * dpuk%x
     &            + ti2 * dpuk%y
     &            + ti3 * dpuk%z
     &            - tk1 * dpui%x
     &            - tk2 * dpui%y
     &            - tk3 * dpui%z
      depy =  ti4 * dpuk%x
     &            + ti5 * dpuk%y
     &            + ti6 * dpuk%z
     &            - tk4 * dpui%x
     &            - tk5 * dpui%y
     &            - tk6 * dpui%z
      depz =  ti7 * dpuk%x
     &            + ti8 * dpuk%y
     &            + ti9 * dpuk%z
     &            - tk7 * dpui%x
     &            - tk8 * dpui%y
     &            - tk9 * dpui%z

      frc%x = frc%x + depx
      frc%y = frc%y + depy
      frc%z = frc%z + depz

c
c     t the dtau/dr terms used for mutual polarization force
c
      dterm1 = bn2 - usc3 * rr5
      dterm2 = bn3 - usc5 * rr7

c
c     Straight terms ( xx, yy ,zz )
c

      dterm3x = usr5 + dterm1
      dterm3y = usr5 + dterm1
      dterm3z = usr5 + dterm1
      dterm4x = rr3  * uscale
      dterm4y = rr3  * uscale
      dterm4z = rr3  * uscale

      dterm5x = - pos%x   * dterm3x
     &             + rc3x  * dterm4x
      dterm5y = - pos%y   * dterm3y
     &             + rc3y  * dterm4y
      dterm5z = - pos%z   * dterm3z
     &             + rc3z  * dterm4z

      dterm6x = - usr5 + pos%x ** 2 * dterm2
     &             - rr5  * pos%x      * urc5x
      dterm6y = - usr5 + pos%y ** 2 * dterm2
     &             - rr5  * pos%y      * urc5y
      dterm6z = - usr5 + pos%z ** 2 * dterm2
     &             - rr5  * pos%z      * urc5z

      tisx =  dpui%x * dterm5x
     &            + uri * dterm6x
      tisy =  dpui%y * dterm5y
     &            + uri * dterm6y
      tisz =  dpui%z * dterm5z
     &            + uri * dterm6z
      tksx =  dpuk%x * dterm5x
     &            + urk * dterm6x
      tksy =  dpuk%y * dterm5y
     &            + urk * dterm6y
      tksz =  dpuk%z * dterm5z
     &            + urk * dterm6z
c
c     Crossed terms ( xy = yx , xz = zx ,yz = zy )
c
      tmp1x  = pos%x
      tmp1y  = pos%x
      tmp1z  = pos%y

      tmp2x  = urc3x
      tmp2y  = urc3x
      tmp2z  = urc3y

      tmp4x  = pos%y
      tmp4y  = pos%z
      tmp4z  = pos%z

      tmp5x  = urc5x
      tmp5y  = urc5x
      tmp5z  = urc5y

      tmp6x  = pos%x   * dterm2
      tmp6y  = pos%x   * dterm2
      tmp6z  = pos%y   * dterm2
      dterm4x = - usr5  * tmp4x
      dterm4y = - usr5  * tmp4y
      dterm4z = - usr5  * tmp4z
      dterm5x = - tmp1x * dterm1
     &             +  rr3  * tmp2x
      dterm5y = - tmp1y * dterm1
     &             +  rr3  * tmp2y
      dterm5z = - tmp1z * dterm1
     &             +  rr3  * tmp2z
      dterm6x =   tmp6x * tmp4x
     &             - rr5   * tmp4x
     &                           * tmp5x
      dterm6y =   tmp6y * tmp4y
     &             - rr5   * tmp4y
     &                           * tmp5y
      dterm6z =   tmp6z * tmp4z
     &             - rr5   * tmp4z
     &                           * tmp5z
      ticx =  dpui%x * dterm4x
     &            + dpui%y * dterm5x
     &            + uri * dterm6x
      ticy =  dpui%x * dterm4y
     &            + dpui%z * dterm5y
     &            + uri * dterm6y
      ticz =  dpui%y * dterm4z
     &            + dpui%z * dterm5z
     &            + uri * dterm6z

      tkcx =  dpuk%x * dterm4x
     &            + dpuk%y * dterm5x
     &            + urk * dterm6x

      tkcy =  dpuk%x * dterm4y
     &            + dpuk%z * dterm5y
     &            + urk * dterm6y

      tkcz =  dpuk%y * dterm4z
     &            + dpuk%z * dterm5z
     &            + urk * dterm6z
c
c      Construct matrixes for dot_product
c
      ti1 = tisx !xx
      ti5 = tisy !yy
      ti9 = tisz !zz
      ti2 = ticx !xy
      ti4 = ticx !yx
      ti3 = ticy !xz
      ti7 = ticy !zx
      ti6 = ticz !yz
      ti8 = ticz !zy

      tk1 = tksx !xx
      tk5 = tksy !yy
      tk9 = tksz !zz
      tk2 = tkcx !xy
      tk4 = tkcx !yx
      tk3 = tkcy !xz
      tk7 = tkcy !zx
      tk6 = tkcz !yz
      tk8 = tkcz !zy
c
c     dot product
c
      depx =  ti1 * dpuk%xx
     &            + ti2 * dpuk%yy
     &            + ti3 * dpuk%zz
     &            + tk1 * dpui%xx
     &            + tk2 * dpui%yy
     &            + tk3 * dpui%zz
      depy =  ti4 * dpuk%xx
     &            + ti5 * dpuk%yy
     &            + ti6 * dpuk%zz
     &            + tk4 * dpui%xx
     &            + tk5 * dpui%yy
     &            + tk6 * dpui%zz
      depz =  ti7 * dpuk%xx
     &            + ti8 * dpuk%yy
     &            + ti9 * dpuk%zz
     &            + tk7 * dpui%xx
     &            + tk8 * dpui%yy
     &            + tk9 * dpui%zz

      frc%x   = frc%x + depx
      frc%y   = frc%y + depy
      frc%z   = frc%z + depz
      frc_r%x = WRITE_C(frc_r%x +) tp2mdr(frc%x)
      frc_r%y = WRITE_C(frc_r%y +) tp2mdr(frc%y)
      frc_r%z = WRITE_C(frc_r%z +) tp2mdr(frc%z)
c
c     t the induced dipole field used for dipole torques
c
      ti3x     =  psr3*dpuk%x  +  dsr3*dpuk%xx
      ti3y     =  psr3*dpuk%y  +  dsr3*dpuk%yy
      ti3z     =  psr3*dpuk%z  +  dsr3*dpuk%zz
      tk3x     =  psr3*dpui%x  +  dsr3*dpui%xx
      tk3y     =  psr3*dpui%y  +  dsr3*dpui%yy
      tk3z     =  psr3*dpui%z  +  dsr3*dpui%zz
      turi5    = -psr5*urk     -  dsr5*urkp
      turk5    = -psr5*uri     -  dsr5*urip
c
c     t induced dipole field gradient used for quadrupole torques
c
      ti5x     =  two * (psr5 * dpuk%x  + dsr5 * dpuk%xx)
      ti5y     =  two * (psr5 * dpuk%y  + dsr5 * dpuk%yy)
      ti5z     =  two * (psr5 * dpuk%z  + dsr5 * dpuk%zz)
      tk5x     =  two * (psr5 * dpui%x  + dsr5 * dpui%xx)
      tk5y     =  two * (psr5 * dpui%y  + dsr5 * dpui%yy)
      tk5z     =  two * (psr5 * dpui%z  + dsr5 * dpui%zz)
      turi7    = - psr7 * urk    - dsr7 * urkp
      turk7    = - psr7 * uri    - dsr7 * urip
c
c     Torque is induced field and gradient cross permanent moments
c
      ufli%x   = ti3x  + pos%x * turi5
      ufli%y   = ti3y  + pos%y * turi5
      ufli%z   = ti3z  + pos%z * turi5
      uflk%x   = tk3x  + pos%x * turk5
      uflk%y   = tk3y  + pos%y * turk5
      uflk%z   = tk3z  + pos%z * turk5

      dufli%x  = pos%x * ti5x     + pos%x**2 * turi7
      dufli%z  = pos%y * ti5y     + pos%y**2 * turi7
      dufli%zz = pos%z * ti5z     + pos%z**2 * turi7
      dufli%y  = pos%x * ti5y     + pos%y    * ti5x
     &          + two*pos%x*pos%y*turi7
      dufli%xx = pos%x * ti5z     + pos%z    * ti5x
     &          + two*pos%x*pos%z*turi7
      dufli%yy = pos%y * ti5z     + pos%z    * ti5y
     &          + two*pos%y*pos%z*turi7

      duflk%x  = pos%x * tk5x     + pos%x**2 * turk7
      duflk%z  = pos%y * tk5y     + pos%y**2 * turk7
      duflk%zz = pos%z * tk5z     + pos%z**2 * turk7
      duflk%y  = pos%x * tk5y     + pos%y    * tk5x
     &          + two*pos%x*pos%y*turk7
      duflk%xx = pos%x * tk5z     + pos%z    * tk5x
     &          + two*pos%x*pos%z*turk7
      duflk%yy = pos%y * tk5z     + pos%z    * tk5y
     &          + two*pos%y*pos%z*turk7

      duflk%x  = -duflk%x
      duflk%y  = -duflk%y
      duflk%z  = -duflk%z
      duflk%xx = -duflk%xx
      duflk%yy = -duflk%yy
      duflk%zz = -duflk%zz

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
      real(t_p) bn0,bn1,bn2,bn3,bn4
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
      bn0      = 0.0_ti_p
      bn1      = 0.0_ti_p
      bn2      = 0.0_ti_p
      bn3      = 0.0_ti_p
      bn4      = 0.0_ti_p
      else
      ralpha   = aewald * r
      !call erfcore_inl(ralpha, bn0,1)
      exp2a    = f_exp( - ralpha**2)
      bn0      = f_erfc(ralpha)

      bn0      = bn0 * invr
      bn1      = ( 1.0_ti_p*bn0 + alsq2    *alsq2n*exp2a ) * invr2
      bn2      = ( 3.0_ti_p*bn1 + alsq2**2 *alsq2n*exp2a ) * invr2
      bn3      = ( 5.0_ti_p*bn2 + alsq2**2*alsq2    *alsq2n*exp2a )
     &           * invr2
      bn4      = ( 7.0_ti_p*bn3 + alsq2**2*alsq2**2 *alsq2n*exp2a )
     &           * invr2

      bn0      = f * bn0
      bn1      = f * bn1
      bn2      = f * bn2
      bn3      = f * bn3
      bn4      = f * bn4
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
         psr3  = bn1 - sc3 * pscale * rr3
         psr5  = bn2 - sc5 * pscale * rr5
         psr7  = bn3 - sc7 * pscale * rr7
      else
         psr3  = bn1 - (1.0_ti_p - sc3 * pscale) * rr3
         psr5  = bn2 - (1.0_ti_p - sc5 * pscale) * rr5
         psr7  = bn3 - (1.0_ti_p - sc7 * pscale) * rr7
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
#endif
