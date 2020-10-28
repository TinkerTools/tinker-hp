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

