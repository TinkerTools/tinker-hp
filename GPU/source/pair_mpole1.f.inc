#include "tinker_cudart.h"
#include "tinker_types.h"
#include "convert.f.inc"
#include "switch_respa.f.inc"

      M_subroutine
     &             mpole1_couple(r2,xr,yr,zr,ip,kp,c_mscale
     &                ,aewald,f,alsq2n,alsq2
     &                ,e,frc,frc_r,ttmi,ttmk,do_correct)
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
         real(t_p)      ,intent(in ):: r2,xr,yr,zr,f,alsq2n,alsq2
         real(t_p)      ,intent(in ):: c_mscale,aewald
         type(rpole_elt),intent(in ):: ip,kp
         logical        ,intent(in ):: do_correct
#ifdef TINKER_CUF
         ener_rtyp      ,intent(out):: e
#else
         real(t_p)      ,intent(out):: e
#endif
         type(real3)    ,intent(out):: frc,ttmi,ttmk
         type(mdyn3_r)  ,intent(out):: frc_r
         real(t_p) two
         real(t_p) de,e_
         real(t_p) r,invr,invr2
         real(t_p) rr1,rr2,rr3,rr5,rr7,rr9,rr11
         real(t_p) bn0,bn1,bn2,bn3,bn4,bn5
         real(t_p) alsqt,ralpha,exp2a
         real(t_p) dikx,diky,dikz
         real(t_p) dirx,diry,dirz
         real(t_p) dkrx,dkry,dkrz
         real(t_p) qrix,qriy,qriz
         real(t_p) qrkx,qrky,qrkz
         real(t_p) qrixr,qriyr,qrizr
         real(t_p) qrkxr,qrkyr,qrkzr
         real(t_p) qrrx,qrry,qrrz
         real(t_p) qikrx,qikry,qikrz
         real(t_p) qkirx,qkiry,qkirz
         real(t_p) qikrxr,qikryr,qikrzr
         real(t_p) qkirxr,qkiryr,qkirzr
         real(t_p) diqkx,diqky,diqkz
         real(t_p) dkqix,dkqiy,dkqiz
         real(t_p) diqkxr,diqkyr,diqkzr
         real(t_p) dkqixr,dkqiyr,dkqizr
         real(t_p) dqiqkx,dqiqky,dqiqkz
         real(t_p) dri,drk,qrri,qrrk
         real(t_p) diqrk,dkqri
         real(t_p) dik,qik,qrrik
         real(t_p) term1,term2,term3
         real(t_p) term4,term5,term6
         parameter(two=2.0)
c
c        get reciprocal distance terms for this interaction
c
         invr2  = r2**(-1)
         r      = f_sqrt(r2)
         invr   = f_sqrt(invr2)
         rr1    = c_mscale*f * invr
         rr3    =       rr1 * invr2
         rr5    = 3.0 * rr3 * invr2
         rr7    = 5.0 * rr5 * invr2
         rr9    = 7.0 * rr7 * invr2
         rr11   = 9.0 * rr9 * invr2
c
c        calculate the real space Ewald error function terms
c
         if (do_correct) then
            bn0 = 0.0
            bn1 = 0.0
            bn2 = 0.0
            bn3 = 0.0
            bn4 = 0.0
            bn5 = 0.0
         else
            ralpha = aewald * r
            exp2a  = f_exp(-ralpha**2)
            !call erfcore_inl(ralpha,bn0,1)
            bn0    = f_erfc(ralpha)
            bn0    = bn0*invr

            alsqt  = alsq2 * alsq2n
            bn1    = (    bn0+alsqt*exp2a) * invr2
            alsqt  = alsq2 * alsqt
            bn2    = (3.0*bn1+alsqt*exp2a) * invr2
            alsqt  = alsq2 * alsqt
            bn3    = (5.0*bn2+alsqt*exp2a) * invr2
            alsqt  = alsq2 * alsqt
            bn4    = (7.0*bn3+alsqt*exp2a) * invr2
            alsqt  = alsq2 * alsqt
            bn5    = (9.0*bn4+alsqt*exp2a) * invr2

            bn0    = f * bn0
            bn1    = f * bn1
            bn2    = f * bn2
            bn3    = f * bn3
            bn4    = f * bn4
            bn5    = f * bn5
         end if
c
c        modify distances to account for Ewald and exclusions
c
         rr1    = bn0 - rr1
         rr3    = bn1 - rr3
         rr5    = bn2 - rr5
         rr7    = bn3 - rr7
         rr9    = bn4 - rr9
         rr11   = bn5 - rr11
c
c        intermediates involving moments and distance separation
c
         dikx   = ip%dy*kp%dz - ip%dz*kp%dy
         diky   = ip%dz*kp%dx - ip%dx*kp%dz
         dikz   = ip%dx*kp%dy - ip%dy*kp%dx
         dirx   = ip%dy*zr  - ip%dz*yr
         diry   = ip%dz*xr  - ip%dx*zr
         dirz   = ip%dx*yr  - ip%dy*xr
         dkrx   = kp%dy*zr  - kp%dz*yr
         dkry   = kp%dz*xr  - kp%dx*zr
         dkrz   = kp%dx*yr  - kp%dy*xr

         dri    = ip%dx*xr  + ip%dy*yr  + ip%dz*zr
         drk    = kp%dx*xr  + kp%dy*yr  + kp%dz*zr
         dik    = ip%dx*kp%dx + ip%dy*kp%dy + ip%dz*kp%dz
         qrix   = ip%qxx*xr + ip%qxy*yr + ip%qxz*zr
         qriy   = ip%qxy*xr + ip%qyy*yr + ip%qyz*zr
         qriz   = ip%qxz*xr + ip%qyz*yr + ip%qzz*zr
         qrkx   = kp%qxx*xr + kp%qxy*yr + kp%qxz*zr
         qrky   = kp%qxy*xr + kp%qyy*yr + kp%qyz*zr
         qrkz   = kp%qxz*xr + kp%qyz*yr + kp%qzz*zr
         qrri   = qrix*xr + qriy*yr + qriz*zr
         qrrk   = qrkx*xr + qrky*yr + qrkz*zr
         qrrik  = qrix*qrkx + qriy*qrky + qriz*qrkz
         qik    = 2.0*(ip%qxy*kp%qxy + ip%qxz*kp%qxz + ip%qyz*kp%qyz)
     &               + ip%qxx*kp%qxx + ip%qyy*kp%qyy + ip%qzz*kp%qzz
         qrixr  = qriz*yr - qriy*zr
         qriyr  = qrix*zr - qriz*xr
         qrizr  = qriy*xr - qrix*yr
         qrkxr  = qrkz*yr - qrky*zr
         qrkyr  = qrkx*zr - qrkz*xr
         qrkzr  = qrky*xr - qrkx*yr

         qrrx   = qrky*qriz - qrkz*qriy
         qrry   = qrkz*qrix - qrkx*qriz
         qrrz   = qrkx*qriy - qrky*qrix

         qikrx  = ip%qxx*qrkx + ip%qxy*qrky + ip%qxz*qrkz
         qikry  = ip%qxy*qrkx + ip%qyy*qrky + ip%qyz*qrkz
         qikrz  = ip%qxz*qrkx + ip%qyz*qrky + ip%qzz*qrkz
         qkirx  = kp%qxx*qrix + kp%qxy*qriy + kp%qxz*qriz
         qkiry  = kp%qxy*qrix + kp%qyy*qriy + kp%qyz*qriz
         qkirz  = kp%qxz*qrix + kp%qyz*qriy + kp%qzz*qriz

         qikrxr = qikrz*yr  - qikry*zr
         qikryr = qikrx*zr  - qikrz*xr
         qikrzr = qikry*xr  - qikrx*yr
         qkirxr = qkirz*yr  - qkiry*zr
         qkiryr = qkirx*zr  - qkirz*xr
         qkirzr = qkiry*xr  - qkirx*yr

         diqkx  = ip%dx*kp%qxx  + ip%dy*kp%qxy + ip%dz*kp%qxz
         diqky  = ip%dx*kp%qxy  + ip%dy*kp%qyy + ip%dz*kp%qyz
         diqkz  = ip%dx*kp%qxz  + ip%dy*kp%qyz + ip%dz*kp%qzz
         dkqix  = kp%dx*ip%qxx  + kp%dy*ip%qxy + kp%dz*ip%qxz
         dkqiy  = kp%dx*ip%qxy  + kp%dy*ip%qyy + kp%dz*ip%qyz
         dkqiz  = kp%dx*ip%qxz  + kp%dy*ip%qyz + kp%dz*ip%qzz
         diqrk  = ip%dx*qrkx  + ip%dy*qrky + ip%dz*qrkz
         dkqri  = kp%dx*qrix  + kp%dy*qriy + kp%dz*qriz

         diqkxr = diqkz*yr  - diqky*zr
         diqkyr = diqkx*zr  - diqkz*xr
         diqkzr = diqky*xr  - diqkx*yr
         dkqixr = dkqiz*yr  - dkqiy*zr
         dkqiyr = dkqix*zr  - dkqiz*xr
         dkqizr = dkqiy*xr  - dkqix*yr

         dqiqkx = ip%dy*qrkz - ip%dz*qrky +  kp%dy*qriz - kp%dz*qriy
     &     - two* (ip%qxy*kp%qxz + ip%qyy*kp%qyz + ip%qyz*kp%qzz
     &           - ip%qxz*kp%qxy - ip%qyz*kp%qyy - ip%qzz*kp%qyz)
         dqiqky = ip%dz*qrkx - ip%dx*qrkz + kp%dz*qrix - kp%dx*qriz
     &     - two*(ip%qxz*kp%qxx + ip%qyz*kp%qxy + ip%qzz*kp%qxz
     &          - ip%qxx*kp%qxz - ip%qxy*kp%qyz - ip%qxz*kp%qzz)
         dqiqkz = ip%dx*qrky - ip%dy*qrkx + kp%dx*qriy - kp%dy*qrix
     &     - two*(ip%qxx*kp%qxy + ip%qxy*kp%qyy + ip%qxz*kp%qyz
     &          - ip%qxy*kp%qxx - ip%qyy*kp%qxy - ip%qyz*kp%qxz)
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
         e_     = term1*rr1 + term2*rr3 + term3*rr5
     &                      + term4*rr7 + term5*rr9
         e      = WRITE_C(e + tp2enr) (e_)
         de     = term1*rr3 + term2*rr5 + term3*rr7
     &          + term4*rr9 + term5*rr11
c
c        calculate intermediate terms for force and torque
c
         term1  = -kp%c*rr3 + drk*rr5 - qrrk*rr7
         term2  =  ip%c*rr3 + dri*rr5 + qrri*rr7
         term3  = two * rr5
         term4  = two * (-kp%c*rr5+drk*rr7-qrrk*rr9)
         term5  = two * (-ip%c*rr5-dri*rr7-qrri*rr9)
         term6  = 4.0 * rr7
c
c        compute the force components for this interaction
c
         frc%x  = de*xr 
     &          + term1*ip%dx         + term2*kp%dx
     &          + term3*(diqkx-dkqix) + term4*qrix
     &          + term5*qrkx          + term6*(qikrx+qkirx)
         frc%y  = de*yr 
     &          + term1*ip%dy         + term2*kp%dy
     &          + term3*(diqky-dkqiy) + term4*qriy
     &          + term5*qrky          + term6*(qikry+qkiry)
         frc%z  = de*zr 
     &          + term1*ip%dz         + term2*kp%dz
     &          + term3*(diqkz-dkqiz) + term4*qriz
     &          + term5*qrkz          + term6*(qikrz+qkirz)
c
c        compute the torque components for this interaction
c
         ttmi%x = WRITE_C(ttmi%x) - rr3*dikx 
     &                            + term1*dirx  + term3*(dqiqkx+dkqixr)
     &                            - term4*qrixr - term6*(qikrxr+qrrx)
         ttmi%y = WRITE_C(ttmi%y) - rr3*diky 
     &                            + term1*diry  + term3*(dqiqky+dkqiyr)
     &                            - term4*qriyr - term6*(qikryr+qrry)
         ttmi%z = WRITE_C(ttmi%z) - rr3*dikz 
     &                            + term1*dirz  + term3*(dqiqkz+dkqizr)
     &                            - term4*qrizr - term6*(qikrzr+qrrz)
         ttmk%x = WRITE_C(ttmk%x) + rr3*dikx
     &                            + term2*dkrx  - term3*(dqiqkx+diqkxr)
     &                            - term5*qrkxr - term6*(qkirxr-qrrx)
         ttmk%y = WRITE_C(ttmk%y) + rr3*diky
     &                            + term2*dkry  - term3*(dqiqky+diqkyr)
     &                            - term5*qrkyr - term6*(qkiryr-qrry)
         ttmk%z = WRITE_C(ttmk%z) + rr3*dikz
     &                            + term2*dkrz  - term3*(dqiqkz+diqkzr)
     &                            - term5*qrkzr - term6*(qkirzr-qrrz)
         !store in large container for mixed precision
         frc_r%x  = WRITE_C( frc_r%x - ) tp2mdr( frc%x )
         frc_r%y  = WRITE_C( frc_r%y - ) tp2mdr( frc%y )
         frc_r%z  = WRITE_C( frc_r%z - ) tp2mdr( frc%z )
      end subroutine

      M_subroutine
     &             mpole3_couple(r2,xr,yr,zr,ip,kp,c_mscale
     &                ,aewald,f,alsq2n,alsq2,e,do_correct)
!$acc routine
      use tinheader ,only: ti_p
      use utilgpu   ,only: rpole_elt,real3,real6
#ifdef TINKER_CUF
      use utilcu ,only: f_erfc
#  if defined(SINGLE)||defined(MIXED)
      use utilcu ,only: f_sqrt,f_exp
#  endif
#endif
         implicit none
         real(t_p)      ,intent(in ):: r2,xr,yr,zr,f,alsq2n,alsq2
         real(t_p)      ,intent(in ):: c_mscale,aewald
         type(rpole_elt),intent(in ):: ip,kp
         logical        ,intent(in ):: do_correct
#ifdef TINKER_CUF
         ener_rtyp      ,intent(out):: e
#else
         real(t_p)      ,intent(out):: e
#endif
         real(t_p) two
         real(t_p) de,e_
         real(t_p) r,invr,invr2
         real(t_p) rr1,rr2,rr3,rr5,rr7,rr9,rr11
         real(t_p) bn0,bn1,bn2,bn3,bn4,bn5
         real(t_p) alsqt,ralpha,exp2a
         real(t_p) qrix,qriy,qriz
         real(t_p) qrkx,qrky,qrkz
         real(t_p) dri,drk,qrri,qrrk
         real(t_p) diqrk,dkqri
         real(t_p) dik,qik,qrrik
         real(t_p) term1,term2,term3
         real(t_p) term4,term5,term6
         parameter(two=2.0)
c
c        get reciprocal distance terms for this interaction
c
         invr2  = r2**(-1)
         r      = f_sqrt(r2)
         invr   = f_sqrt(invr2)
         rr1    = c_mscale*f * invr
         rr3    =       rr1 * invr2
         rr5    = 3.0 * rr3 * invr2
         rr7    = 5.0 * rr5 * invr2
         rr9    = 7.0 * rr7 * invr2
         rr11   = 9.0 * rr9 * invr2
c
c        calculate the real space Ewald error function terms
c
         if (do_correct) then
            bn0 = 0.0
            bn1 = 0.0
            bn2 = 0.0
            bn3 = 0.0
            bn4 = 0.0
            bn5 = 0.0
         else
            ralpha = aewald * r
            exp2a  = f_exp(-ralpha**2)
            !call erfcore_inl(ralpha,bn0,1)
            bn0    = f_erfc(ralpha)
            bn0    = bn0*invr

            alsqt  = alsq2 * alsq2n
            bn1    = (    bn0+alsqt*exp2a) * invr2
            alsqt  = alsq2 * alsqt
            bn2    = (3.0*bn1+alsqt*exp2a) * invr2
            alsqt  = alsq2 * alsqt
            bn3    = (5.0*bn2+alsqt*exp2a) * invr2
            alsqt  = alsq2 * alsqt
            bn4    = (7.0*bn3+alsqt*exp2a) * invr2
            alsqt  = alsq2 * alsqt
            bn5    = (9.0*bn4+alsqt*exp2a) * invr2

            bn0    = f * bn0
            bn1    = f * bn1
            bn2    = f * bn2
            bn3    = f * bn3
            bn4    = f * bn4
            bn5    = f * bn5
         end if
c
c        modify distances to account for Ewald and exclusions
c
         rr1    = bn0 - rr1
         rr3    = bn1 - rr3
         rr5    = bn2 - rr5
         rr7    = bn3 - rr7
         rr9    = bn4 - rr9
         rr11   = bn5 - rr11
c
c        intermediates involving moments and distance separation
c
         dri    = ip%dx*xr  + ip%dy*yr  + ip%dz*zr
         drk    = kp%dx*xr  + kp%dy*yr  + kp%dz*zr
         dik    = ip%dx*kp%dx + ip%dy*kp%dy + ip%dz*kp%dz
         qrix   = ip%qxx*xr + ip%qxy*yr + ip%qxz*zr
         qriy   = ip%qxy*xr + ip%qyy*yr + ip%qyz*zr
         qriz   = ip%qxz*xr + ip%qyz*yr + ip%qzz*zr
         qrkx   = kp%qxx*xr + kp%qxy*yr + kp%qxz*zr
         qrky   = kp%qxy*xr + kp%qyy*yr + kp%qyz*zr
         qrkz   = kp%qxz*xr + kp%qyz*yr + kp%qzz*zr
         qrri   = qrix*xr + qriy*yr + qriz*zr
         qrrk   = qrkx*xr + qrky*yr + qrkz*zr
         qrrik  = qrix*qrkx + qriy*qrky + qriz*qrkz
         qik    = 2.0*(ip%qxy*kp%qxy + ip%qxz*kp%qxz + ip%qyz*kp%qyz)
     &               + ip%qxx*kp%qxx + ip%qyy*kp%qyy + ip%qzz*kp%qzz

         diqrk  = ip%dx*qrkx  + ip%dy*qrky + ip%dz*qrkz
         dkqri  = kp%dx*qrix  + kp%dy*qriy + kp%dz*qriz
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
         e_     = term1*rr1 + term2*rr3 + term3*rr5
     &          + term4*rr7 + term5*rr9
         e      = WRITE_C(e + tp2enr) (e_)
      end subroutine

      M_subroutine
     &             mpole1_couple_shortlong(r2,xr,yr,zr,ip,kp,c_mscale
     &                         ,r_cut,shortheal,aewald,f,alsq2n,alsq2
     &                         ,e,frc,frc_r,ttmi,ttmk,do_correct,mode)
!$acc routine
      use interfaces,only: m_long
      use tinheader ,only: ti_p
      use utilgpu   ,only: rpole_elt,real3,real6,mdyn3_r
#ifdef TINKER_CUF
      use utilcu ,only: f_erfc
#  if defined(SINGLE)||defined(MIXED)
      use utilcu ,only: f_sqrt,f_exp
#  endif
#endif
         implicit none
         real(t_p)      ,intent(in ):: r2,xr,yr,zr,f,alsq2n,alsq2
         real(t_p)      ,intent(in ):: c_mscale,aewald,shortheal,r_cut
         type(rpole_elt),intent(in ):: ip,kp
         logical        ,intent(in ):: do_correct
         integer        ,intent(in ):: mode
#ifdef TINKER_CUF
         ener_rtyp      ,intent(out):: e
#else
         real(t_p)      ,intent(out):: e
#endif
         type(real3)    ,intent(out):: frc,ttmi,ttmk
         type(mdyn3_r)  ,intent(out):: frc_r
         real(t_p) two
         real(t_p) e_,de,s,ds
         real(t_p) r,invr,invr2
         real(t_p) rr1,rr2,rr3,rr5,rr7,rr9,rr11
         real(t_p) bn0,bn1,bn2,bn3,bn4,bn5
         real(t_p) alsqt,ralpha,exp2a
         real(t_p) dikx,diky,dikz
         real(t_p) dirx,diry,dirz
         real(t_p) dkrx,dkry,dkrz
         real(t_p) qrix,qriy,qriz
         real(t_p) qrkx,qrky,qrkz
         real(t_p) qrixr,qriyr,qrizr
         real(t_p) qrkxr,qrkyr,qrkzr
         real(t_p) qrrx,qrry,qrrz
         real(t_p) qikrx,qikry,qikrz
         real(t_p) qkirx,qkiry,qkirz
         real(t_p) qikrxr,qikryr,qikrzr
         real(t_p) qkirxr,qkiryr,qkirzr
         real(t_p) diqkx,diqky,diqkz
         real(t_p) dkqix,dkqiy,dkqiz
         real(t_p) diqkxr,diqkyr,diqkzr
         real(t_p) dkqixr,dkqiyr,dkqizr
         real(t_p) dqiqkx,dqiqky,dqiqkz
         real(t_p) dri,drk,qrri,qrrk
         real(t_p) diqrk,dkqri
         real(t_p) dik,qik,qrrik
         real(t_p) term1,term2,term3
         real(t_p) term4,term5,term6
         parameter(two=2.0)
c
c        get reciprocal distance terms for this interaction
c
         invr2  = 1.0/r2
         r      = f_sqrt(r2)
         invr   = f_sqrt(invr2)
         rr1    = c_mscale*f * invr
         rr3    =       rr1 * invr2
         rr5    = 3.0 * rr3 * invr2
         rr7    = 5.0 * rr5 * invr2
         rr9    = 7.0 * rr7 * invr2
         rr11   = 9.0 * rr9 * invr2
c
c        calculate the real space Ewald error function terms
c
         if (do_correct) then
            bn0 = 0.0
            bn1 = 0.0
            bn2 = 0.0
            bn3 = 0.0
            bn4 = 0.0
            bn5 = 0.0
         else
            ralpha = aewald * r
            exp2a  = f_exp(-ralpha**2)
            !call erfcore_inl(ralpha,bn0,1)
            bn0    = f_erfc(ralpha)
            bn0    = bn0*invr

            alsqt  = alsq2 * alsq2n
            bn1    = (    bn0+alsqt*exp2a) * invr2
            alsqt  = alsq2 * alsqt
            bn2    = (3.0*bn1+alsqt*exp2a) * invr2
            alsqt  = alsq2 * alsqt
            bn3    = (5.0*bn2+alsqt*exp2a) * invr2
            alsqt  = alsq2 * alsqt
            bn4    = (7.0*bn3+alsqt*exp2a) * invr2
            alsqt  = alsq2 * alsqt
            bn5    = (9.0*bn4+alsqt*exp2a) * invr2

            bn0    = f * bn0
            bn1    = f * bn1
            bn2    = f * bn2
            bn3    = f * bn3
            bn4    = f * bn4
            bn5    = f * bn5
         end if
c
c        modify distances to account for Ewald and exclusions
c
         rr1    = bn0 - rr1
         rr3    = bn1 - rr3
         rr5    = bn2 - rr5
         rr7    = bn3 - rr7
         rr9    = bn4 - rr9
         rr11   = bn5 - rr11
c
c        intermediates involving moments and distance separation
c
         dikx   = ip%dy*kp%dz - ip%dz*kp%dy
         diky   = ip%dz*kp%dx - ip%dx*kp%dz
         dikz   = ip%dx*kp%dy - ip%dy*kp%dx
         dirx   = ip%dy*zr  - ip%dz*yr
         diry   = ip%dz*xr  - ip%dx*zr
         dirz   = ip%dx*yr  - ip%dy*xr
         dkrx   = kp%dy*zr  - kp%dz*yr
         dkry   = kp%dz*xr  - kp%dx*zr
         dkrz   = kp%dx*yr  - kp%dy*xr

         dri    = ip%dx*xr  + ip%dy*yr  + ip%dz*zr
         drk    = kp%dx*xr  + kp%dy*yr  + kp%dz*zr
         dik    = ip%dx*kp%dx + ip%dy*kp%dy + ip%dz*kp%dz
         qrix   = ip%qxx*xr + ip%qxy*yr + ip%qxz*zr
         qriy   = ip%qxy*xr + ip%qyy*yr + ip%qyz*zr
         qriz   = ip%qxz*xr + ip%qyz*yr + ip%qzz*zr
         qrkx   = kp%qxx*xr + kp%qxy*yr + kp%qxz*zr
         qrky   = kp%qxy*xr + kp%qyy*yr + kp%qyz*zr
         qrkz   = kp%qxz*xr + kp%qyz*yr + kp%qzz*zr
         qrri   = qrix*xr + qriy*yr + qriz*zr
         qrrk   = qrkx*xr + qrky*yr + qrkz*zr
         qrrik  = qrix*qrkx + qriy*qrky + qriz*qrkz
         qik    = 2.0*(ip%qxy*kp%qxy + ip%qxz*kp%qxz + ip%qyz*kp%qyz)
     &               + ip%qxx*kp%qxx + ip%qyy*kp%qyy + ip%qzz*kp%qzz
         qrixr  = qriz*yr - qriy*zr
         qriyr  = qrix*zr - qriz*xr
         qrizr  = qriy*xr - qrix*yr
         qrkxr  = qrkz*yr - qrky*zr
         qrkyr  = qrkx*zr - qrkz*xr
         qrkzr  = qrky*xr - qrkx*yr

         qrrx   = qrky*qriz - qrkz*qriy
         qrry   = qrkz*qrix - qrkx*qriz
         qrrz   = qrkx*qriy - qrky*qrix

         qikrx  = ip%qxx*qrkx + ip%qxy*qrky + ip%qxz*qrkz
         qikry  = ip%qxy*qrkx + ip%qyy*qrky + ip%qyz*qrkz
         qikrz  = ip%qxz*qrkx + ip%qyz*qrky + ip%qzz*qrkz
         qkirx  = kp%qxx*qrix + kp%qxy*qriy + kp%qxz*qriz
         qkiry  = kp%qxy*qrix + kp%qyy*qriy + kp%qyz*qriz
         qkirz  = kp%qxz*qrix + kp%qyz*qriy + kp%qzz*qriz

         qikrxr = qikrz*yr  - qikry*zr
         qikryr = qikrx*zr  - qikrz*xr
         qikrzr = qikry*xr  - qikrx*yr
         qkirxr = qkirz*yr  - qkiry*zr
         qkiryr = qkirx*zr  - qkirz*xr
         qkirzr = qkiry*xr  - qkirx*yr

         diqkx  = ip%dx*kp%qxx  + ip%dy*kp%qxy + ip%dz*kp%qxz
         diqky  = ip%dx*kp%qxy  + ip%dy*kp%qyy + ip%dz*kp%qyz
         diqkz  = ip%dx*kp%qxz  + ip%dy*kp%qyz + ip%dz*kp%qzz
         dkqix  = kp%dx*ip%qxx  + kp%dy*ip%qxy + kp%dz*ip%qxz
         dkqiy  = kp%dx*ip%qxy  + kp%dy*ip%qyy + kp%dz*ip%qyz
         dkqiz  = kp%dx*ip%qxz  + kp%dy*ip%qyz + kp%dz*ip%qzz
         diqrk  = ip%dx*qrkx  + ip%dy*qrky + ip%dz*qrkz
         dkqri  = kp%dx*qrix  + kp%dy*qriy + kp%dz*qriz

         diqkxr = diqkz*yr  - diqky*zr
         diqkyr = diqkx*zr  - diqkz*xr
         diqkzr = diqky*xr  - diqkx*yr
         dkqixr = dkqiz*yr  - dkqiy*zr
         dkqiyr = dkqix*zr  - dkqiz*xr
         dkqizr = dkqiy*xr  - dkqix*yr

         dqiqkx = ip%dy*qrkz - ip%dz*qrky +  kp%dy*qriz - kp%dz*qriy
     &     - two* (ip%qxy*kp%qxz + ip%qyy*kp%qyz + ip%qyz*kp%qzz
     &           - ip%qxz*kp%qxy - ip%qyz*kp%qyy - ip%qzz*kp%qyz)
         dqiqky = ip%dz*qrkx - ip%dx*qrkz + kp%dz*qrix - kp%dx*qriz
     &     - two*(ip%qxz*kp%qxx + ip%qyz*kp%qxy + ip%qzz*kp%qxz
     &          - ip%qxx*kp%qxz - ip%qxy*kp%qyz - ip%qxz*kp%qzz)
         dqiqkz = ip%dx*qrky - ip%dy*qrkx + kp%dx*qriy - kp%dy*qrix
     &     - two*(ip%qxx*kp%qxy + ip%qxy*kp%qyy + ip%qxz*kp%qyz
     &          - ip%qxy*kp%qxx - ip%qyy*kp%qxy - ip%qyz*kp%qxz)
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
         call switch_respa_inl(r,r_cut,shortheal,s,ds)
         if (mode.eq.m_long) s = 1.0-s

         e_     = term1*rr1 + term2*rr3 + term3*rr5
     &                      + term4*rr7 + term5*rr9

         e      = WRITE_C(e + tp2enr) ( + s*e_ )
         de     = term1*rr3 + term2*rr5 + term3*rr7
     &          + term4*rr9 + term5*rr11
c
c        calculate intermediate terms for force and torque
c
         term1  = -kp%c*rr3 + drk*rr5 - qrrk*rr7
         term2  =  ip%c*rr3 + dri*rr5 + qrri*rr7
         term3  = two * rr5
         term4  = two * (-kp%c*rr5+drk*rr7-qrrk*rr9)
         term5  = two * (-ip%c*rr5-dri*rr7-qrri*rr9)
         term6  = 4.0 * rr7
c
c        compute the force components for this interaction
c
         frc%x  = s*(de*xr 
     &          + term1*ip%dx         + term2*kp%dx
     &          + term3*(diqkx-dkqix) + term4*qrix
     &          + term5*qrkx          + term6*(qikrx+qkirx))
     &          
         frc%y  = s*(de*yr 
     &          + term1*ip%dy         + term2*kp%dy
     &          + term3*(diqky-dkqiy) + term4*qriy
     &          + term5*qrky          + term6*(qikry+qkiry))
     &         
         frc%z  = s*(de*zr 
     &          + term1*ip%dz         + term2*kp%dz
     &          + term3*(diqkz-dkqiz) + term4*qriz
     &          + term5*qrkz          + term6*(qikrz+qkirz))
     &          
         if (mode.eq.m_long) then
            frc%x  = frc%x + ds*xr*e_*invr
            frc%y  = frc%y + ds*yr*e_*invr
            frc%z  = frc%z + ds*zr*e_*invr
         else
            frc%x  = frc%x - ds*xr*e_*invr
            frc%y  = frc%y - ds*yr*e_*invr
            frc%z  = frc%z - ds*zr*e_*invr
         end if
c
c        compute the torque components for this interaction
c
         ttmi%x = WRITE_C(ttmi%x +) s*( - rr3*dikx 
     &                            + term1*dirx  + term3*(dqiqkx+dkqixr)
     &                            - term4*qrixr - term6*(qikrxr+qrrx))
         ttmi%y = WRITE_C(ttmi%y +) s*( - rr3*diky 
     &                            + term1*diry  + term3*(dqiqky+dkqiyr)
     &                            - term4*qriyr - term6*(qikryr+qrry))
         ttmi%z = WRITE_C(ttmi%z +) s*( - rr3*dikz 
     &                            + term1*dirz  + term3*(dqiqkz+dkqizr)
     &                            - term4*qrizr - term6*(qikrzr+qrrz))
         ttmk%x = WRITE_C(ttmk%x +) s*( + rr3*dikx
     &                            + term2*dkrx  - term3*(dqiqkx+diqkxr)
     &                            - term5*qrkxr - term6*(qkirxr-qrrx))
         ttmk%y = WRITE_C(ttmk%y +) s*( + rr3*diky
     &                            + term2*dkry  - term3*(dqiqky+diqkyr)
     &                            - term5*qrkyr - term6*(qkiryr-qrry))
         ttmk%z = WRITE_C(ttmk%z +) s*( + rr3*dikz
     &                            + term2*dkrz  - term3*(dqiqkz+diqkzr)
     &                            - term5*qrkzr - term6*(qkirzr-qrrz))
         !store in large container for mixed precision
         frc_r%x  = WRITE_C( frc_r%x - ) tp2mdr( frc%x )
         frc_r%y  = WRITE_C( frc_r%y - ) tp2mdr( frc%y )
         frc_r%z  = WRITE_C( frc_r%z - ) tp2mdr( frc%z )
      end subroutine

      M_subroutine
     &             mpole3_couple_shortlong(r2,xr,yr,zr,ip,kp,c_mscale
     &                         ,r_cut,shortheal,aewald,f,alsq2n,alsq2
     &                         ,e,do_correct,mode)
!$acc routine
      use interfaces,only: m_long
      use tinheader ,only: ti_p
      use utilgpu   ,only: rpole_elt,real3,real6
#ifdef TINKER_CUF
      use utilcu ,only: f_erfc
#  if defined(SINGLE)||defined(MIXED)
      use utilcu ,only: f_sqrt,f_exp
#  endif
#endif
         implicit none
         real(t_p)      ,intent(in ):: r2,xr,yr,zr,f,alsq2n,alsq2
         real(t_p)      ,intent(in ):: c_mscale,aewald,shortheal,r_cut
         type(rpole_elt),intent(in ):: ip,kp
         logical        ,intent(in ):: do_correct
         integer        ,intent(in ):: mode
#ifdef TINKER_CUF
         ener_rtyp      ,intent(out):: e
#else
         real(t_p)      ,intent(out):: e
#endif
         real(t_p) two
         real(t_p) e_,s,ds
         real(t_p) r,invr,invr2
         real(t_p) rr1,rr2,rr3,rr5,rr7,rr9,rr11
         real(t_p) bn0,bn1,bn2,bn3,bn4,bn5
         real(t_p) alsqt,ralpha,exp2a
         real(t_p) qrix,qriy,qriz
         real(t_p) qrkx,qrky,qrkz
         real(t_p) dri,drk,qrri,qrrk
         real(t_p) diqrk,dkqri
         real(t_p) dik,qik,qrrik
         real(t_p) term1,term2,term3
         real(t_p) term4,term5,term6
         parameter(two=2.0)
c
c        get reciprocal distance terms for this interaction
c
         invr2  = 1.0/r2
         r      = f_sqrt(r2)
         invr   = f_sqrt(invr2)
         rr1    = c_mscale*f * invr
         rr3    =       rr1 * invr2
         rr5    = 3.0 * rr3 * invr2
         rr7    = 5.0 * rr5 * invr2
         rr9    = 7.0 * rr7 * invr2
         rr11   = 9.0 * rr9 * invr2
c
c        calculate the real space Ewald error function terms
c
         if (do_correct) then
            bn0 = 0.0
            bn1 = 0.0
            bn2 = 0.0
            bn3 = 0.0
            bn4 = 0.0
            bn5 = 0.0
         else
            ralpha = aewald * r
            exp2a  = exp(-ralpha**2)
            !call erfcore_inl(ralpha,bn0,1)
            bn0    = f_erfc(ralpha)
            bn0    = bn0*invr

            alsqt  = alsq2 * alsq2n
            bn1    = (    bn0+alsqt*exp2a) * invr2
            alsqt  = alsq2 * alsqt
            bn2    = (3.0*bn1+alsqt*exp2a) * invr2
            alsqt  = alsq2 * alsqt
            bn3    = (5.0*bn2+alsqt*exp2a) * invr2
            alsqt  = alsq2 * alsqt
            bn4    = (7.0*bn3+alsqt*exp2a) * invr2
            alsqt  = alsq2 * alsqt
            bn5    = (9.0*bn4+alsqt*exp2a) * invr2

            bn0    = f * bn0
            bn1    = f * bn1
            bn2    = f * bn2
            bn3    = f * bn3
            bn4    = f * bn4
            bn5    = f * bn5
         end if
c
c        modify distances to account for Ewald and exclusions
c
         rr1    = bn0 - rr1
         rr3    = bn1 - rr3
         rr5    = bn2 - rr5
         rr7    = bn3 - rr7
         rr9    = bn4 - rr9
         rr11   = bn5 - rr11
c
c        intermediates involving moments and distance separation
c
         dri    = ip%dx*xr  + ip%dy*yr  + ip%dz*zr
         drk    = kp%dx*xr  + kp%dy*yr  + kp%dz*zr
         dik    = ip%dx*kp%dx + ip%dy*kp%dy + ip%dz*kp%dz
         qrix   = ip%qxx*xr + ip%qxy*yr + ip%qxz*zr
         qriy   = ip%qxy*xr + ip%qyy*yr + ip%qyz*zr
         qriz   = ip%qxz*xr + ip%qyz*yr + ip%qzz*zr
         qrkx   = kp%qxx*xr + kp%qxy*yr + kp%qxz*zr
         qrky   = kp%qxy*xr + kp%qyy*yr + kp%qyz*zr
         qrkz   = kp%qxz*xr + kp%qyz*yr + kp%qzz*zr
         qrri   = qrix*xr + qriy*yr + qriz*zr
         qrrk   = qrkx*xr + qrky*yr + qrkz*zr
         qrrik  = qrix*qrkx + qriy*qrky + qriz*qrkz
         qik    = 2.0*(ip%qxy*kp%qxy + ip%qxz*kp%qxz + ip%qyz*kp%qyz)
     &               + ip%qxx*kp%qxx + ip%qyy*kp%qyy + ip%qzz*kp%qzz

         diqrk  = ip%dx*qrkx  + ip%dy*qrky + ip%dz*qrkz
         dkqri  = kp%dx*qrix  + kp%dy*qriy + kp%dz*qriz
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
         call switch_respa_inl(r,r_cut,shortheal,s,ds)
         if (mode.eq.m_long) s = 1.0-s

         e_     = term1*rr1 + term2*rr3 + term3*rr5
     &                      + term4*rr7 + term5*rr9

         e      = WRITE_C(e + tp2enr) ( + s*e_ )
      end subroutine

