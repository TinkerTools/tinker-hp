#include "tinker_cudart.h"

        M_subroutine
     &              tmatxb_couple(d2,dist,dpui,dpuk,
     &                    sdamp,pgamma,aewald,alsq2,alsq2n,uscale,
     &                          fid,fip,fkd,fkp,do_correct)
!$acc routine
#if  defined(TINKER_CUF)
           use utilcu  ,only: f_erfc
#  if (defined(SINGLE)||defined(MIXED))
           use utilcu  ,only: f_sqrt,f_exp
#  endif
#endif
           use utilgpu ,only: real3,real6
           implicit none
           type(real6),intent(in)   ::dpui,dpuk
           type(real3),intent(in)   ::dist
           real(t_p)  ,intent(in)   ::sdamp,pgamma,aewald,alsq2,alsq2n
           real(t_p)  ,intent(in)   ::uscale
           real(t_p)  ,intent(inout)::d2
           type(real3),intent(inout)::fid,fip,fkd,fkp
           logical    ,intent(in)   :: do_correct
           real(t_p) d1,ralpha,exp2a
           real(t_p) bn0,bn1,bn2
           real(t_p) sdamp1,expdamp1
           real(t_p) rr3,rr5,rr3_bn1,rr5_bn2,sc3,sc5
           real(t_p) duir,dukr,puir,pukr  ! Scalar products duir = (du(i).r)

           ! compute the distances and the scaling factors according to Thole's model.
           d1      = f_sqrt(d2)
           d2      = d2**(-1)

           ralpha  = aewald * d1
           d1      = d1**(-1)

           if (do_correct) then
              bn1     = 0.0
              bn2     = 0.0

              if( sdamp == 0.0 ) then
                sdamp1  = -100.0
                sc3      =   f_exp(sdamp1) * uscale
                sc5      =   f_exp(sdamp1) * uscale * (1 - sdamp1)
              else
                sdamp1 = - pgamma / ((d1*sdamp)*(d1*sdamp)*(d1*sdamp))
                if (sdamp1 > -50.0) then
                  expdamp1 = f_exp(sdamp1)
                  sc3      =   expdamp1 * uscale
                  sc5      =   expdamp1 * uscale * (1 - sdamp1)
                else
                  sc3     = 1
                  sc5     = 1
                end if
              end if
           else
              exp2a   = f_exp(-ralpha*ralpha)
              !call erfcore_inl(ralpha,bn0,1)
              bn0     = f_erfc(ralpha)   ! Attention to Macro f_erfc in single prec with hastings method (check tinker_cudart.h)
              bn0     = bn0 * d1
              bn1     = (     bn0 +         alsq2 * alsq2n * exp2a) * d2
              bn2     = ( 3 * bn1 + alsq2 * alsq2 * alsq2n * exp2a) * d2

              if( sdamp == 0.0 ) then
                sdamp1   = -100.0
                expdamp1 = f_exp(sdamp1)
                sc3      =   1 - expdamp1 * uscale
                sc5      =   1 - expdamp1 * uscale * (1 - sdamp1)
              else
                sdamp1 = - pgamma / ((d1*sdamp)*(d1*sdamp)*(d1*sdamp))
                if (sdamp1 > -50.0) then
                  expdamp1 = f_exp(sdamp1)
                  sc3      =   1 - expdamp1 * uscale
                  sc5      =   1 - expdamp1 * uscale * (1 - sdamp1)
                else
                  sc3     = 1
                  sc5     = 1
                end if
              end if
           end if

           ! compute the field.
           rr3     =     (1 - sc3) * (d1 * d2)
           rr5     = 3 * (1 - sc5) * (d1 * d2 * d2)
           rr3_bn1 = rr3 - bn1
           rr5_bn2 = rr5 - bn2

           duir    = dpui%x * dist%x + dpui%y * dist%y + dpui%z *dist%z
           dukr    = dpuk%x * dist%x + dpuk%y * dist%y + dpuk%z *dist%z

           puir    = dpui%xx * dist%x + dpui%yy * dist%y
     &             + dpui%zz * dist%z
           pukr    = dpuk%xx * dist%x + dpuk%yy * dist%y
     &             + dpuk%zz * dist%z

        fid%x  = WRITE_C(fid%x) - rr3_bn1*dpuk%x + rr5_bn2*dukr*dist%x
        fid%y  = WRITE_C(fid%y) - rr3_bn1*dpuk%y + rr5_bn2*dukr*dist%y
        fid%z  = WRITE_C(fid%z) - rr3_bn1*dpuk%z + rr5_bn2*dukr*dist%z

        fip%x  = WRITE_C(fip%x) - rr3_bn1*dpuk%xx + rr5_bn2*pukr*dist%x
        fip%y  = WRITE_C(fip%y) - rr3_bn1*dpuk%yy + rr5_bn2*pukr*dist%y
        fip%z  = WRITE_C(fip%z) - rr3_bn1*dpuk%zz + rr5_bn2*pukr*dist%z

        fkd%x  = WRITE_C(fkd%x) - rr3_bn1*dpui%x + rr5_bn2*duir*dist%x
        fkd%y  = WRITE_C(fkd%y) - rr3_bn1*dpui%y + rr5_bn2*duir*dist%y
        fkd%z  = WRITE_C(fkd%z) - rr3_bn1*dpui%z + rr5_bn2*duir*dist%z

        fkp%x  = WRITE_C(fkp%x) - rr3_bn1*dpui%xx + rr5_bn2*puir*dist%x
        fkp%y  = WRITE_C(fkp%y) - rr3_bn1*dpui%yy + rr5_bn2*puir*dist%y
        fkp%z  = WRITE_C(fkp%z) - rr3_bn1*dpui%zz + rr5_bn2*puir*dist%z
        end subroutine

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
        parameter( one=1.0_ti_p , two = 2.0_ti_p)

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

