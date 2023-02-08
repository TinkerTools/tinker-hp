#ifndef PAIR_POLAR_CPEN_INC
#define PAIR_POLAR_CPEN_INC

#include "tinker_macro.h"
#include "tinker_cudart.h"
#include "damping.inc.f"

      M_subroutine
     &        duo_polar_cpen(ui,uip,mi,uk,ukp,mk,p,r2,f
     &           ,corei,corek,vali,valk,alphai,alphak,dscal,pscal,wscal
     &           ,ewald,pentyp,ucflx
     &           ,poti,potk,e,frc,trqi,trqk,ver,fea)
      
      use tinTypes,only: rpole_elt,real3,real6,mdyn3_r
#ifdef TINKER_CUF
      use utilcu  ,only: f_erfc
#  if defined(SINGLE)||defined(MIXED)
      use utilcu  ,only: f_sqrt,f_exp
#  endif
#endif
      implicit none
      integer    ,intent(in ):: pentyp,ver,fea
      logical    ,intent(in ):: ucflx
      real(t_p)  ,intent(in ):: r2,f,corei,corek,vali,valk,alphai,alphak
     &           ,pscal,dscal,wscal,ewald
      type(rpole_elt),intent(in):: mi,mk
      type(real3),intent(in ):: ui,uip,uk,ukp,p
      real(t_p)  ,intent(out):: poti,potk,e
      type(real3),intent(inout):: frc,trqi,trqk

      integer     grd,ene,sca
      real(t_p)   dir,qix,qiy,qiz,qir,dkr,qkx,qky,qkz,qkr,uir,ukr
     &           ,r,rr1,rr3,rr3core,rr5,rr5core,rr7,rr9
     &           ,rr3i,rr5i,rr7i,rr9i,rr3k,rr5k,rr7k,rr9k,rr5ik,rr7ik
     &           ,tix3,tiy3,tiz3,tix5,tiy5,tiz5,tuir,tukr
     &           ,tkx3,tky3,tkz3,tkx5,tky5,tkz5
     &           ,term1,term2,term3,term1core
     &          ,term1i,term2i,term3i,term4i,term5i,term6i,term7i,term8i
     &          ,term1k,term2k,term3k,term4k,term5k,term6k,term7k,term8k
     &           ,tixx,tixy,tixz,tiyy,tiyz,tizz
     &           ,tkxx,tkxy,tkxz,tkyy,tkyz,tkzz
     &           ,depx,depy,depz
     &           ,dmpe(5),dmpi(5),dmpk(5),dmpik(5),dufi(6),dufk(6)
      type(real3) ufi,ufk
      parameter(grd=__use_grd__,ene=__use_ene__,sca=__use_sca__)

      r    = f_sqrt(r2)
c
c     intermediates involving moments and separation distance
c
      dir  = mi%dx *p%x + mi%dy *p%y + mi%dz *p%z
      qix  = mi%qxx*p%x + mi%qxy*p%y + mi%qxz*p%z
      qiy  = mi%qxy*p%x + mi%qyy*p%y + mi%qyz*p%z
      qiz  = mi%qxz*p%x + mi%qyz*p%y + mi%qzz*p%z
      qir  =    qix*p%x +    qiy*p%y +    qiz*p%z
      dkr  = mk%dx *p%x + mk%dy *p%y + mk%dz *p%z
      qkx  = mk%qxx*p%x + mk%qxy*p%y + mk%qxz*p%z
      qky  = mk%qxy*p%x + mk%qyy*p%y + mk%qyz*p%z
      qkz  = mk%qxz*p%x + mk%qyz*p%y + mk%qzz*p%z
      qkr  =    qkx*p%x +    qky*p%y +    qkz*p%z
      uir  =  ui%x *p%x +  ui%y *p%y + ui%z *p%z
      ukr  =  uk%x *p%x +  uk%y *p%y + uk%z *p%z
c
c     get reciprocal distance terms for this interaction
c
      rr1  = f / r
      rr3  = rr1 / r2
      rr5  = 3.0 * rr3 / r2
      rr7  = 5.0 * rr5 / r2
      rr9  = 7.0 * rr7 / r2
c
c     calculate real space Ewald error function damping
c
      call dampewald_inl (9,r,r2,ewald,f,dmpe)
c
c     apply charge penetration damping to scale factors
c
      call damppole_inl (r,9,pentyp,alphai,alphak,dmpi,dmpk,dmpik)
      IF (IAND(ver,sca).NE.0) THEN
         rr3core =                dscal         *rr3
         rr5core =                dscal         *rr5
         rr3i    =                dscal*dmpi (2)*rr3
         rr5i    =                dscal*dmpi (3)*rr5
         rr7i    =                dscal*dmpi (4)*rr7
         rr9i    =                dscal*dmpi (5)*rr9
         rr3k    =                dscal*dmpk (2)*rr3
         rr5k    =                dscal*dmpk (3)*rr5
         rr7k    =                dscal*dmpk (4)*rr7
         rr9k    =                dscal*dmpk (5)*rr9
         rr5ik   =                wscal*dmpik(3)*rr5
         rr7ik   =                wscal*dmpik(4)*rr7
      ELSE
         rr3core = dmpe(2) - (1.0-1.0           )*rr3
         rr5core = dmpe(3) - (1.0-1.0           )*rr5
         rr3i    = dmpe(2) - (1.0-      dmpi (2))*rr3
         rr5i    = dmpe(3) - (1.0-      dmpi (3))*rr5
         rr7i    = dmpe(4) - (1.0-      dmpi (4))*rr7
         rr9i    = dmpe(5) - (1.0-      dmpi (5))*rr9
         rr3k    = dmpe(2) - (1.0-      dmpk (2))*rr3
         rr5k    = dmpe(3) - (1.0-      dmpk (3))*rr5
         rr7k    = dmpe(4) - (1.0-      dmpk (4))*rr7
         rr9k    = dmpe(5) - (1.0-      dmpk (5))*rr9
         rr5ik   = dmpe(3) - (1.0-      dmpik(3))*rr5
         rr7ik   = dmpe(4) - (1.0-      dmpik(4))*rr7
      ENDIF
c
c     store the potential at each site for use in charge flux
c
      if (ucflx) then
         poti =-2.0 *ukr *rr3i
         potk = 2.0 *uir *rr3k
      end if
      tix3  =  2.0*rr3i*uk%x
      tiy3  =  2.0*rr3i*uk%y
      tiz3  =  2.0*rr3i*uk%z
      tkx3  =  2.0*rr3k*ui%x
      tky3  =  2.0*rr3k*ui%y
      tkz3  =  2.0*rr3k*ui%z
      tuir  = -2.0*rr5i*ukr
      tukr  = -2.0*rr5k*uir
c
c     find the field components for charge penetration damping
c
      call dampdir_inl (pentyp,r,alphai,alphak,dmpi,dmpk)
      IF (IAND(ver,sca).NE.0) THEN
      rr3i  =                pscal*dmpi(1) *rr3
      rr5i  =                pscal*dmpi(2) *rr5
      rr7i  =                pscal*dmpi(3) *rr7
      rr3k  =                pscal*dmpk(1) *rr3
      rr5k  =                pscal*dmpk(2) *rr5
      rr7k  =                pscal*dmpk(3) *rr7
      rr3   =                pscal         *rr3
      ELSE
      rr3i  = dmpe(2) - (1.0-      dmpi(1))*rr3
      rr5i  = dmpe(3) - (1.0-      dmpi(2))*rr5
      rr7i  = dmpe(4) - (1.0-      dmpi(3))*rr7
      rr3k  = dmpe(2) - (1.0-      dmpk(1))*rr3
      rr5k  = dmpe(3) - (1.0-      dmpk(2))*rr5
      rr7k  = dmpe(4) - (1.0-      dmpk(3))*rr7
      rr3   = dmpe(2)
      ENDIF

      ufi%x = -p%x*(rr3*corek + rr3k*valk - rr5k*dkr + rr7k*qkr)
     &            - rr3k*mk%dx + 2.0*rr5k*qkx
      ufi%y = -p%y*(rr3*corek + rr3k*valk - rr5k*dkr + rr7k*qkr)
     &            - rr3k*mk%dy + 2.0*rr5k*qky
      ufi%z = -p%z*(rr3*corek + rr3k*valk - rr5k*dkr + rr7k*qkr)
     &            - rr3k*mk%dz + 2.0*rr5k*qkz
      ufk%x =  p%x*(rr3*corei + rr3i*vali + rr5i*dir + rr7i*qir)
     &            - rr3i*mi%dx - 2.0*rr5i*qix
      ufk%y =  p%y*(rr3*corei + rr3i*vali + rr5i*dir + rr7i*qir)
     &            - rr3i*mi%dy - 2.0*rr5i*qiy
      ufk%z =  p%z*(rr3*corei + rr3i*vali + rr5i*dir + rr7i*qir)
     &            - rr3i*mi%dz - 2.0*rr5i*qiz

c
c     compute the energy contribution for this interaction
c
      e = -( ui%x*ufi%x + uk%x*ufk%x + ui%y*ufi%y + uk%y*ufk%y
     &      +ui%z*ufi%z + uk%z*ufk%z )
c
c     get the induced dipole field used for dipole torques
c
      ufi%x = tix3 + p%x*tuir
      ufi%y = tiy3 + p%y*tuir
      ufi%z = tiz3 + p%z*tuir
      ufk%x = tkx3 + p%x*tukr
      ufk%y = tky3 + p%y*tukr
      ufk%z = tkz3 + p%z*tukr
c
c     get induced dipole field gradient used for quadrupole torques
c
      tix5 =  4.0 *(rr5i*uk%x)
      tiy5 =  4.0 *(rr5i*uk%y)
      tiz5 =  4.0 *(rr5i*uk%z)
      tkx5 =  4.0 *(rr5k*ui%x)
      tky5 =  4.0 *(rr5k*ui%y)
      tkz5 =  4.0 *(rr5k*ui%z)
      tuir = -2.0 *rr7i*ukr 
      tukr = -2.0 *rr7k*uir 

      dufi(1)= + p%x*tix5 + p%x*p%x*tuir
      dufi(2)= + p%x*tiy5 + p%y*tix5 + 2.0*p%x*p%y*tuir
      dufi(3)= + p%y*tiy5 + p%y*p%y*tuir
      dufi(4)= + p%x*tiz5 + p%z*tix5 + 2.0*p%x*p%z*tuir
      dufi(5)= + p%y*tiz5 + p%z*tiy5 + 2.0*p%y*p%z*tuir
      dufi(6)= + p%z*tiz5 + p%z*p%z*tuir
      dufk(1)= - p%x*tkx5 - p%x*p%x*tukr
      dufk(2)= - p%x*tky5 - p%y*tkx5 - 2.0*p%x*p%y*tukr
      dufk(3)= - p%y*tky5 - p%y*p%y*tukr
      dufk(4)= - p%x*tkz5 - p%z*tkx5 - 2.0*p%x*p%z*tukr
      dufk(5)= - p%y*tkz5 - p%z*tky5 - 2.0*p%y*p%z*tukr
      dufk(6)= - p%z*tkz5 - p%z*p%z*tukr
c
c     get the field gradient for direct polarization force
c
      term1i = rr3i - rr5i*p%x*p%x
      term1core = rr3core - rr5core*p%x*p%x
      term2i = 2.0*rr5i*p%x 
      term3i = rr7i*p%x*p%x - rr5i
      term4i = 2.0*rr5i
      term5i = 5.0*rr7i*p%x
      term6i = rr9i*p%x*p%x
      term1k = rr3k - rr5k*p%x*p%x
      term2k = 2.0*rr5k*p%x
      term3k = rr7k*p%x*p%x - rr5k
      term4k = 2.0*rr5k
      term5k = 5.0*rr7k*p%x
      term6k = rr9k*p%x*p%x
      tixx   = vali*term1i + corei*term1core  
     &              + mi%dx*term2i - dir*term3i
     &              - mi%qxx*term4i + qix*term5i - qir*term6i
     &              + (qiy*p%y+qiz*p%z)*rr7i
      tkxx   = valk*term1k + corek*term1core
     &              - mk%dx*term2k + dkr*term3k
     &              - mk%qxx*term4k + qkx*term5k - qkr*term6k
     &              + (qky*p%y+qkz*p%z)*rr7k
      term1i = rr3i - rr5i*p%y*p%y
      term1core = rr3core - rr5core*p%y*p%y
      term2i = 2.0*rr5i*p%y
      term3i = rr7i*p%y*p%y - rr5i
      term4i = 2.0*rr5i
      term5i = 5.0*rr7i*p%y
      term6i = rr9i*p%y*p%y
      term1k = rr3k - rr5k*p%y*p%y
      term2k = 2.0*rr5k*p%y
      term3k = rr7k*p%y*p%y - rr5k
      term4k = 2.0*rr5k
      term5k = 5.0*rr7k*p%y
      term6k = rr9k*p%y*p%y
      tiyy   = vali*term1i + corei*term1core
     &              + mi%dy*term2i - dir*term3i
     &              - mi%qyy*term4i + qiy*term5i - qir*term6i
     &              + (qix*p%x+qiz*p%z)*rr7i
      tkyy   = valk*term1k + corek*term1core
     &              - mk%dy*term2k + dkr*term3k
     &              - mk%qyy*term4k + qky*term5k - qkr*term6k
     &              + (qkx*p%x+qkz*p%z)*rr7k
      term1i = rr3i - rr5i*p%z*p%z
      term1core = rr3core - rr5core*p%z*p%z
      term2i = 2.0*rr5i*p%z
      term3i = rr7i*p%z*p%z - rr5i
      term4i = 2.0*rr5i
      term5i = 5.0*rr7i*p%z
      term6i = rr9i*p%z*p%z
      term1k = rr3k - rr5k*p%z*p%z
      term2k = 2.0*rr5k*p%z
      term3k = rr7k*p%z*p%z - rr5k
      term4k = 2.0*rr5k
      term5k = 5.0*rr7k*p%z
      term6k = rr9k*p%z*p%z
      tizz   = vali*term1i + corei*term1core
     &              + mi%dz*term2i - dir*term3i
     &              - mi%qzz*term4i + qiz*term5i - qir*term6i
     &              + (qix*p%x+qiy*p%y)*rr7i
      tkzz   = valk*term1k + corek*term1core
     &              - mk%dz*term2k + dkr*term3k
     &              - mk%qzz*term4k + qkz*term5k - qkr*term6k
     &              + (qkx*p%x+qky*p%y)*rr7k
      term2i = rr5i*p%x 
      term1i = p%y * term2i
      term1core = rr5core*p%x*p%y
      term3i = rr5i*p%y
      term4i = p%y * (rr7i*p%x)
      term5i = 2.0*rr5i
      term6i = 2.0*rr7i*p%x
      term7i = 2.0*rr7i*p%y
      term8i = p%y*rr9i*p%x
      term2k = rr5k*p%x
      term1k = p%y * term2k
      term3k = rr5k*p%y
      term4k = p%y * (rr7k*p%x)
      term5k = 2.0*rr5k
      term6k = 2.0*rr7k*p%x
      term7k = 2.0*rr7k*p%y
      term8k = p%y*rr9k*p%x
      tixy   = -vali*term1i - corei*term1core 
     &              + mi%dy*term2i + mi%dx*term3i
     &              - dir*term4i - mi%qxy*term5i + qiy*term6i
     &              + qix*term7i - qir*term8i
      tkxy   = -valk*term1k - corek*term1core 
     &              - mk%dy*term2k - mk%dx*term3k
     &              + dkr*term4k - mk%qxy*term5k + qky*term6k
     &              + qkx*term7k - qkr*term8k
      term2i = rr5i*p%x
      term1i = p%z * term2i
      term1core = rr5core*p%x*p%z
      term3i = rr5i*p%z
      term4i = p%z * (rr7i*p%x)
      term5i = 2.0*rr5i
      term6i = 2.0*rr7i*p%x
      term7i = 2.0*rr7i*p%z
      term8i = p%z*rr9i*p%x
      term2k = rr5k*p%x
      term1k = p%z * term2k
      term3k = rr5k*p%z
      term4k = p%z * (rr7k*p%x)
      term5k = 2.0*rr5k
      term6k = 2.0*rr7k*p%x
      term7k = 2.0*rr7k*p%z
      term8k = p%z*rr9k*p%x
      tixz   = -vali*term1i - corei*term1core
     &              + mi%dz*term2i + mi%dx*term3i
     &              - dir*term4i - mi%qxz*term5i + qiz*term6i
     &              + qix*term7i - qir*term8i
      tkxz   = -valk*term1k - corek*term1core
     &              - mk%dz*term2k - mk%dx*term3k
     &              + dkr*term4k - mk%qxz*term5k + qkz*term6k
     &              + qkx*term7k - qkr*term8k
      term2i = rr5i*p%y
      term1i = p%z * term2i
      term1core = rr5core*p%y*p%z
      term3i = rr5i*p%z
      term4i = p%z * (rr7i*p%y)
      term5i = 2.0*rr5i
      term6i = 2.0*rr7i*p%y
      term7i = 2.0*rr7i*p%z
      term8i = p%z*rr9i*p%y
      term2k = rr5k*p%y
      term1k = p%z * term2k
      term3k = rr5k*p%z
      term4k = p%z * (rr7k*p%y)
      term5k = 2.0*rr5k
      term6k = 2.0*rr7k*p%y
      term7k = 2.0*rr7k*p%z
      term8k = p%z*rr9k*p%y
      tiyz   = -vali*term1i - corei*term1core
     &              + mi%dz*term2i + mi%dy*term3i
     &              - dir*term4i - mi%qyz*term5i + qiz*term6i
     &              + qiy*term7i - qir*term8i
      tkyz   = -valk*term1k - corek*term1core
     &              - mk%dz*term2k - mk%dy*term3k
     &              + dkr*term4k - mk%qyz*term5k + qkz*term6k
     &              + qky*term7k - qkr*term8k
      depx   = tixx*uk%x + tixy*uk%y + tixz*uk%z
     &              - tkxx*ui%x - tkxy*ui%y - tkxz*ui%z
      depy   = tixy*uk%x + tiyy*uk%y + tiyz*uk%z
     &              - tkxy*ui%x - tkyy*ui%y - tkyz*ui%z
      depz   = tixz*uk%x + tiyz*uk%y + tizz*uk%z
     &              - tkxz*ui%x - tkyz*ui%y - tkzz*ui%z
      frc%x  = -2.0 * depx
      frc%y  = -2.0 * depy
      frc%z  = -2.0 * depz
c
c     get the dtau/dr terms used for mutual polarization force
c
      term1  = 2.0 * rr5ik
      term2  = term1*p%x
      term3  = rr5ik - rr7ik*p%x*p%x
      tixx   = ui%x*term2 + uir*term3
      tkxx   = uk%x*term2 + ukr*term3
      term2  = term1*p%y
      term3  = rr5ik - rr7ik*p%y*p%y
      tiyy   = ui%y*term2 + uir*term3
      tkyy   = uk%y*term2 + ukr*term3
      term2  = term1*p%z
      term3  = rr5ik - rr7ik*p%z*p%z
      tizz   = ui%z*term2 + uir*term3
      tkzz   = uk%z*term2 + ukr*term3
      term1  = rr5ik*p%y
      term2  = rr5ik*p%x
      term3  = p%y*(rr7ik*p%x)
      tixy   = ui%x*term1 + ui%y*term2 - uir*term3
      tkxy   = uk%x*term1 + uk%y*term2 - ukr*term3
      term1  = rr5ik*p%z
      term3  = p%z*(rr7ik*p%x)
      tixz   = ui%x*term1 + ui%z*term2 - uir*term3
      tkxz   = uk%x*term1 + uk%z*term2 - ukr*term3
      term2  = rr5ik*p%y
      term3  = p%z*(rr7ik*p%y)
      tiyz   = ui%y*term1 + ui%z*term2 - uir*term3
      tkyz   = uk%y*term1 + uk%z*term2 - ukr*term3
      depx   = tixx*ukp%x + tixy*ukp%y + tixz*ukp%z
     &       + tkxx*uip%x + tkxy*uip%y + tkxz*uip%z
      depy   = tixy*ukp%x + tiyy*ukp%y + tiyz*ukp%z
     &       + tkxy*uip%x + tkyy*uip%y + tkyz*uip%z
      depz   = tixz*ukp%x + tiyz*ukp%y + tizz*ukp%z
     &       + tkxz*uip%x + tkyz*uip%y + tkzz*uip%z
      frc%x  = frc%x - depx
      frc%y  = frc%y - depy
      frc%z  = frc%z - depz
c
c     Torque is induced field and gradient cross permanent moments
c
      trqi%x = WRITE_C(trqi%x +) mi%dz*ufi%y - mi%dy*ufi%z
     &       + mi%qxz*dufi(2) - mi%qxy*dufi(4)
     &       + 2.0*mi%qyz*( dufi(3) - dufi(6) )
     &       + (mi%qzz - mi%qyy)*dufi(5)
      trqi%y = WRITE_C(trqi%y +) mi%dx*ufi%z - mi%dz*ufi%x
     &       - mi%qyz*dufi(2) + mi%qxy*dufi(5)
     &       + 2.0*mi%qxz*( dufi(6) - dufi(1) )
     &       + (mi%qxx - mi%qzz)*dufi(4)
      trqi%z = WRITE_C(trqi%z +) mi%dy*ufi%x - mi%dx*ufi%y
     &       + mi%qyz*dufi(4) - mi%qxz*dufi(5)
     &       + 2.0*mi%qxy*( dufi(1) - dufi(3) )
     &       + (mi%qyy - mi%qxx)*dufi(2)

      trqk%x = WRITE_C(trqk%x +) mk%dz*ufk%y - mk%dy*ufk%z
     &       + mk%qxz*dufk(2) - mk%qxy*dufk(4)
     &       + 2.0*mk%qyz*( dufk(3) - dufk(6) )
     &       + (mk%qzz - mk%qyy)*dufk(5)
      trqk%y = WRITE_C(trqk%y +) mk%dx*ufk%z - mk%dz*ufk%x
     &       - mk%qyz*dufk(2)  + mk%qxy*dufk(5)
     &       + 2.0*mk%qxz*( dufk(6) - dufk(1) )
     &       + (mk%qxx - mk%qzz)*dufk(4)
      trqk%z = WRITE_C(trqk%z +) mk%dy*ufk%x - mk%dx*ufk%y
     &       + mk%qyz*dufk(4) - mk%qxz*dufk(5)
     &       + 2.0*mk%qxy*( dufk(1) - dufk(3) )
     &       + (mk%qyy - mk%qxx)*dufk(2)

      end subroutine
#endif
