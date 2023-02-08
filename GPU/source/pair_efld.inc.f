#ifndef PAIR_EFLD_INC_F
#define PAIR_EFLD_INC_F
#include "tinker_cudart.h"
#include "damping.inc.f"

        M_subroutine
     &           efld0_couple(d2,pos,ip,kp,alsq2,alsq2n,
     &             aewald,damp,pgamma,dscale,pscale,
     &             fid,fip,fkd,fkp,d1,bn1,bn2,sc3,sc5,do_correct)
!$acc routine
        use tinheader ,only: ti_p
        use tintypes  ,only: real3,rpole_elt
#if  defined(TINKER_CUF)
        use utilcu  ,only: f_erfc
#  if (defined(SINGLE)||defined(MIXED))
        use utilcu  ,only: f_sqrt,f_exp
#  endif
#endif
        real(t_p)      ,intent(in) :: d2
        type(real3)    ,intent(in) :: pos
        type(rpole_elt),intent(in) :: ip,kp
        real(t_p)      ,intent(in) :: aewald,damp,pgamma
     &                 ,alsq2,alsq2n,pscale,dscale
        logical        ,intent(in) :: do_correct
        type(real3)    ,intent(out)::fid,fip,fkd,fkp
        real(t_p)      ,intent(out):: d1,bn1,bn2,sc3,sc5

        real(t_p) exp2a
        real(t_p) invd1,invd2,invd3,invd5,invd7
        real(t_p) sc7,dsc3,dsc5,dsc7,psc3,psc5,psc7
        real(t_p) drr3,drr5,drr7,prr3,prr5,prr7
        real(t_p) dir,qirr,dkr,qkrr
        real(t_p) qirx,qiry,qirz,qkrx,qkry,qkrz
        real(t_p) fkmx,fkmy,fkmz,fimx,fimy,fimz
        real(t_p) invdamp,expdamp1,damp1
        real(t_p) ralpha,bn0,bn3
        real(t_p) one,two
        parameter( one=1.0_ti_p, two = 2.0_ti_p )

        damp1   = -100.0_ti_p
        invdamp = damp ** (-one)
        invd2   = d2 ** (-one)
        d1      = d2 ** 0.5_ti_p
        invd1   = d1 ** (-one)

        sc3     = one
        sc5     = one
        sc7     = one

        invd3   = invd1 * invd2
        invd5   = invd3 * invd2
        invd7   = invd5 * invd2

        if (damp.ne.0.0_ti_p) damp1 = - pgamma*(d1*invdamp)**3

        if (damp1 > -50.0_ti_p) then
           expdamp1  = f_exp(damp1)
           sc3  = one - expdamp1
           sc5  = one - expdamp1*(one - damp1)
           sc7  = one - expdamp1*(one - damp1 + 0.6_ti_p*damp1**2)
        end if

        if (do_correct) then
           ! [dp]scale equal to 1-[dp]scale in this case
           drr3    =      sc3*dscale * invd3
           drr5    =  3 * sc5*dscale * invd5
           drr7    = 15 * sc7*dscale * invd7

           prr3    =      sc3*pscale * invd3
           prr5    =  3 * sc5*pscale * invd5
           prr7    = 15 * sc7*pscale * invd7
        else
c
c     calculate the error function damping terms
c
           ralpha  = aewald * d1
           exp2a   = f_exp( -ralpha**2 )
           bn0     = f_erfc(ralpha)

           bn0     =    bn0                            * invd1
           bn1     = (  bn0  + alsq2    *alsq2n*exp2a) * invd2
           bn2     = (3*bn1  + alsq2**2 *alsq2n*exp2a) * invd2
           bn3     = (5*bn2  + alsq2**3 *alsq2n*exp2a) * invd2

           drr3    =      (one - sc3*dscale) * invd3
           drr5    =  3 * (one - sc5*dscale) * invd5
           drr7    = 15 * (one - sc7*dscale) * invd7

           prr3    =      (one - sc3*pscale) * invd3
           prr5    =  3 * (one - sc5*pscale) * invd5
           prr7    = 15 * (one - sc7*pscale) * invd7
        end if
c
c     compute some intermediate quantities
c
        dir     =  ip%dx*pos%x +  ip%dy*pos%y +  ip%dz*pos%z
        qirx    = ip%qxx*pos%x + ip%qxy*pos%y + ip%qxz*pos%z
        qiry    = ip%qxy*pos%x + ip%qyy*pos%y + ip%qyz*pos%z
        qirz    = ip%qxz*pos%x + ip%qyz*pos%y + ip%qzz*pos%z
        qirr    =   qirx*pos%x +   qiry*pos%y +   qirz*pos%z

        dkr     =  kp%dx*pos%x +  kp%dy*pos%y +   kp%dz*pos%z
        qkrx    = kp%qxx*pos%x + kp%qxy*pos%y +  kp%qxz*pos%z
        qkry    = kp%qxy*pos%x + kp%qyy*pos%y +  kp%qyz*pos%z
        qkrz    = kp%qxz*pos%x + kp%qyz*pos%y +  kp%qzz*pos%z
        qkrr    =   qkrx*pos%x +   qkry*pos%y +    qkrz*pos%z

        if (do_correct) then
           fimx = 0.0_ti_p; fimy = 0.0_ti_p; fimz = 0.0_ti_p;
           fkmx = 0.0_ti_p; fkmy = 0.0_ti_p; fkmz = 0.0_ti_p;
        else
           fimx = -( bn1*kp%c  - bn2*dkr + bn3*qkrr )*pos%x
     &            -  bn1*kp%dx + two*bn2*qkrx
           fimy = -( bn1*kp%c  - bn2*dkr + bn3*qkrr )*pos%y
     &            -  bn1*kp%dy + two*bn2*qkry
           fimz = -( bn1*kp%c  - bn2*dkr + bn3*qkrr )*pos%z
     &            -  bn1*kp%dz + two*bn2*qkrz
           fkmx =  ( bn1*ip%c  + bn2*dir + bn3*qirr )*pos%x
     &            -  bn1*ip%dx - two*bn2*qirx
           fkmy =  ( bn1*ip%c  + bn2*dir + bn3*qirr )*pos%y
     &            -  bn1*ip%dy - two*bn2*qiry
           fkmz =  ( bn1*ip%c  + bn2*dir + bn3*qirr )*pos%z
     &            -  bn1*ip%dz - two*bn2*qirz
        end if

        fid%x   =  ( drr3*kp%c  - drr5*dkr + drr7*qkrr )*pos%x
     &            +  drr3*kp%dx - two*drr5*qkrx
        fid%y   =  ( drr3*kp%c  - drr5*dkr + drr7*qkrr )*pos%y
     &            +  drr3*kp%dy - two*drr5*qkry
        fid%z   =  ( drr3*kp%c  - drr5*dkr + drr7*qkrr )*pos%z
     &            +  drr3*kp%dz - two*drr5*qkrz
        fip%x   =  ( prr3*kp%c  - prr5*dkr + prr7*qkrr )*pos%x
     &            +  prr3*kp%dx - two*prr5*qkrx
        fip%y   =  ( prr3*kp%c  - prr5*dkr + prr7*qkrr )*pos%y
     &            +  prr3*kp%dy - two*prr5*qkry
        fip%z   =  ( prr3*kp%c  - prr5*dkr + prr7*qkrr )*pos%z
     &            +  prr3*kp%dz - two*prr5*qkrz

        fkd%x   = -( drr3*ip%c  + drr5*dir + drr7*qirr )*pos%x
     &            +  drr3*ip%dx + two*drr5*qirx
        fkd%y   = -( drr3*ip%c  + drr5*dir + drr7*qirr )*pos%y
     &            +  drr3*ip%dy + two*drr5*qiry
        fkd%z   = -( drr3*ip%c  + drr5*dir + drr7*qirr )*pos%z
     &            +  drr3*ip%dz + two*drr5*qirz
        fkp%x   = -( prr3*ip%c  + prr5*dir + prr7*qirr )*pos%x
     &            +  prr3*ip%dx + two*prr5*qirx
        fkp%y   = -( prr3*ip%c  + prr5*dir + prr7*qirr )*pos%y
     &            +  prr3*ip%dy + two*prr5*qiry
        fkp%z   = -( prr3*ip%c  + prr5*dir + prr7*qirr )*pos%z
     &            +  prr3*ip%dz + two*prr5*qirz

        fid%x   =  fimx + fid%x
        fid%y   =  fimy + fid%y
        fid%z   =  fimz + fid%z
        fip%x   =  fimx + fip%x
        fip%y   =  fimy + fip%y
        fip%z   =  fimz + fip%z
        fkd%x   =  fkmx + fkd%x
        fkd%y   =  fkmy + fkd%y
        fkd%z   =  fkmz + fkd%z
        fkp%x   =  fkmx + fkp%x
        fkp%y   =  fkmy + fkp%y
        fkp%z   =  fkmz + fkp%z
        end subroutine

      M_subroutine
     &        duo_efld0(r2,p,mi,mk,pscale,dscale,aewald
     &           ,damp,pgamma,use_dirdamp
     &           ,fid,fkd,fip,fkp,tver,tfea
     &           )
      use tintypes ,only: rpole_elt
#if defined(TINKER_CUF) && (TINKER_SINGLE_PREC+TINKER_MIXED_PREC)
      use utilcu   ,only: f_sqrt
#endif
      integer    ,intent(in):: tver,tfea
      logical    ,intent(in):: use_dirdamp
      real(t_p)  ,intent(in):: r2,pscale,dscale,aewald,pgamma
      real(t_p)  ,intent(inout):: damp
      type(real3),intent(in):: p
      type(rpole_elt),intent(in):: mi,mk
      type(real3),intent(out):: fip,fkp,fid,fkd
      integer     sca
      real(t_p)   r,rr1,rr2,rr3,rr5,rr7
      real(t_p)   dir,qix,qiy,qiz,dkr,qkx,qky,qkz,qir,qkr,dmp3,dmp5,dmp7
      real(t_p)   dmpe(4),dmpik(3)
      parameter (sca=__use_sca__)
c
      r   = f_sqrt(r2)
      rr1 = r **(-1)
      rr2 = rr1 * rr1
      rr3 = rr2 * rr1
      rr5 = 3.0 * rr2 * rr3
      rr7 = 5.0 * rr2 * rr5
c
c     intermediates involving moments and separation distance
c
      dir = mi%dx *p%x + mi%dy *p%y + mi%dz *p%z
      qix = mi%qxx*p%x + mi%qxy*p%y + mi%qxz*p%z
      qiy = mi%qxy*p%x + mi%qyy*p%y + mi%qyz*p%z
      qiz = mi%qxz*p%x + mi%qyz*p%y + mi%qzz*p%z
      dkr = mk%dx *p%x + mk%dy *p%y + mk%dz *p%z
      qkx = mk%qxx*p%x + mk%qxy*p%y + mk%qxz*p%z
      qky = mk%qxy*p%x + mk%qyy*p%y + mk%qyz*p%z
      qkz = mk%qxz*p%x + mk%qyz*p%y + mk%qzz*p%z
      qir = qix*p%x + qiy*p%y + qiz*p%z
      qkr = qkx*p%x + qky*p%y + qkz*p%z
c
c     calculate real space Ewald error function damping
c
      IF (IAND(tver,sca).EQ.0)
     &   call dampewald_inl (7,r,r2,aewald,1.0,dmpe)
c
c     find the field components for Thole polarization damping
c
      call dampthole_inl(7,r,damp,pgamma,use_dirdamp,dmpik)
      IF (IAND(tver,sca).NE.0) THEN
         dmp3 = dscale*dmpik(1)*rr3
         dmp5 = dscale*dmpik(2)*rr5
         dmp7 = dscale*dmpik(3)*rr7
      ELSE
         dmp3 = dmpe(2) - (1.0-dmpik(1))*rr3
         dmp5 = dmpe(3) - (1.0-dmpik(2))*rr5
         dmp7 = dmpe(4) - (1.0-dmpik(3))*rr7
      END IF
      fid%x =-p%x*(dmp3*mk%c-dmp5*dkr+dmp7*qkr)
     &           - dmp3*mk%dx + 2.0*dmp5*qkx
      fid%y =-p%y*(dmp3*mk%c-dmp5*dkr+dmp7*qkr)
     &           - dmp3*mk%dy + 2.0*dmp5*qky
      fid%z =-p%z*(dmp3*mk%c-dmp5*dkr+dmp7*qkr)
     &            - dmp3*mk%dz + 2.0*dmp5*qkz
      fkd%x = p%x*(dmp3*mi%c+dmp5*dir+dmp7*qir)
     &           - dmp3*mi%dx - 2.0*dmp5*qix
      fkd%y = p%y*(dmp3*mi%c+dmp5*dir+dmp7*qir)
     &           - dmp3*mi%dy - 2.0*dmp5*qiy
      fkd%z = p%z*(dmp3*mi%c+dmp5*dir+dmp7*qir)
     &           - dmp3*mi%dz - 2.0*dmp5*qiz
      IF (IAND(tver,sca).NE.0) THEN
         dmp3 = pscale*dmpik(1)*rr3
         dmp5 = pscale*dmpik(2)*rr5
         dmp7 = pscale*dmpik(3)*rr7
      END IF
      fip%x =-p%x*(dmp3*mk%c-dmp5*dkr+dmp7*qkr)
     &           - dmp3*mk%dx + 2.0*dmp5*qkx
      fip%y =-p%y*(dmp3*mk%c-dmp5*dkr+dmp7*qkr)
     &           - dmp3*mk%dy + 2.0*dmp5*qky
      fip%z =-p%z*(dmp3*mk%c-dmp5*dkr+dmp7*qkr)
     &           - dmp3*mk%dz + 2.0*dmp5*qkz
      fkp%x = p%x*(dmp3*mi%c+dmp5*dir+dmp7*qir)
     &           - dmp3*mi%dx - 2.0*dmp5*qix
      fkp%y = p%y*(dmp3*mi%c+dmp5*dir+dmp7*qir)
     &           - dmp3*mi%dy - 2.0*dmp5*qiy
      fkp%z = p%z*(dmp3*mi%c+dmp5*dir+dmp7*qir)
     &           - dmp3*mi%dz - 2.0*dmp5*qiz
      end subroutine

      M_subroutine
     &        duo_efld0NEwd(r2,p,mi,mk,pscale,dscale
     &           ,damp,pgamma,use_dirdamp
     &           ,fid,fkd,fip,fkp,tver,tfea
     &           )
      use tintypes ,only: rpole_elt
#if defined(TINKER_CUF) && (TINKER_SINGLE_PREC+TINKER_MIXED_PREC)
      use utilcu   ,only: f_sqrt
#endif
      integer    ,intent(in):: tver,tfea
      logical    ,intent(in):: use_dirdamp
      real(t_p)  ,intent(in):: r2,pscale,dscale,pgamma
      real(t_p)  ,intent(inout):: damp
      type(real3),intent(in):: p
      type(rpole_elt),intent(in):: mi,mk
      type(real3),intent(out):: fip,fkp,fid,fkd
      integer     sca
      real(t_p)   r,rr1,rr2,rr3,rr5,rr7
      real(t_p)   dir,qix,qiy,qiz,dkr,qkx,qky,qkz,qir,qkr,dmp3,dmp5,dmp7
      real(t_p)   dmpik(3)
      parameter (sca=__use_sca__)
c
      r   = f_sqrt(r2)
      rr1 = r **(-1)
      rr2 = rr1 * rr1
      rr3 = rr2 * rr1
      rr5 = 3.0 * rr2 * rr3
      rr7 = 5.0 * rr2 * rr5
c
c     intermediates involving moments and separation distance
c
      dir = mi%dx *p%x + mi%dy *p%y + mi%dz *p%z
      qix = mi%qxx*p%x + mi%qxy*p%y + mi%qxz*p%z
      qiy = mi%qxy*p%x + mi%qyy*p%y + mi%qyz*p%z
      qiz = mi%qxz*p%x + mi%qyz*p%y + mi%qzz*p%z
      dkr = mk%dx *p%x + mk%dy *p%y + mk%dz *p%z
      qkx = mk%qxx*p%x + mk%qxy*p%y + mk%qxz*p%z
      qky = mk%qxy*p%x + mk%qyy*p%y + mk%qyz*p%z
      qkz = mk%qxz*p%x + mk%qyz*p%y + mk%qzz*p%z
      qir = qix*p%x + qiy*p%y + qiz*p%z
      qkr = qkx*p%x + qky*p%y + qkz*p%z
c
c     find the field components for Thole polarization damping
c
      call dampthole_inl(7,r,damp,pgamma,use_dirdamp,dmpik)
      IF (IAND(tver,sca).NE.0) THEN
         dmp3 = dscale*dmpik(1)*rr3
         dmp5 = dscale*dmpik(2)*rr5
         dmp7 = dscale*dmpik(3)*rr7
      ELSE
         dmp3 = (dmpik(1)-1.0)*rr3
         dmp5 = (dmpik(2)-1.0)*rr5
         dmp7 = (dmpik(3)-1.0)*rr7
      END IF

      fid%x =-p%x*(dmp3*mk%c-dmp5*dkr+dmp7*qkr)
     &           - dmp3*mk%dx + 2.0*dmp5*qkx
      fid%y =-p%y*(dmp3*mk%c-dmp5*dkr+dmp7*qkr)
     &           - dmp3*mk%dy + 2.0*dmp5*qky
      fid%z =-p%z*(dmp3*mk%c-dmp5*dkr+dmp7*qkr)
     &            - dmp3*mk%dz + 2.0*dmp5*qkz
      fkd%x = p%x*(dmp3*mi%c+dmp5*dir+dmp7*qir)
     &           - dmp3*mi%dx - 2.0*dmp5*qix
      fkd%y = p%y*(dmp3*mi%c+dmp5*dir+dmp7*qir)
     &           - dmp3*mi%dy - 2.0*dmp5*qiy
      fkd%z = p%z*(dmp3*mi%c+dmp5*dir+dmp7*qir)
     &           - dmp3*mi%dz - 2.0*dmp5*qiz
      IF (IAND(tver,sca).NE.0) THEN
         dmp3 = pscale*dmpik(1)*rr3
         dmp5 = pscale*dmpik(2)*rr5
         dmp7 = pscale*dmpik(3)*rr7
      END IF
      fip%x =-p%x*(dmp3*mk%c-dmp5*dkr+dmp7*qkr)
     &           - dmp3*mk%dx + 2.0*dmp5*qkx
      fip%y =-p%y*(dmp3*mk%c-dmp5*dkr+dmp7*qkr)
     &           - dmp3*mk%dy + 2.0*dmp5*qky
      fip%z =-p%z*(dmp3*mk%c-dmp5*dkr+dmp7*qkr)
     &           - dmp3*mk%dz + 2.0*dmp5*qkz
      fkp%x = p%x*(dmp3*mi%c+dmp5*dir+dmp7*qir)
     &           - dmp3*mi%dx - 2.0*dmp5*qix
      fkp%y = p%y*(dmp3*mi%c+dmp5*dir+dmp7*qir)
     &           - dmp3*mi%dy - 2.0*dmp5*qiy
      fkp%z = p%z*(dmp3*mi%c+dmp5*dir+dmp7*qir)
     &           - dmp3*mi%dz - 2.0*dmp5*qiz
      end subroutine
#endif
