#ifndef PAIR_TMATXB_INC_F
#define PAIR_TMATXB_INC_F
#include "tinker_cudart.h"
#include "damping.inc.f"

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
           real(t_p) sdamp1,sdamp2,expdamp1,expdamp2
           real(t_p) rr3,rr5,rr3_bn1,rr5_bn2,sc3,sc5
           real(t_p) duir,dukr,puir,pukr  ! Scalar products duir = (du(i).r)
           parameter(sdamp2=-100.0,expdamp2= 0.0)

           ! compute the distances and the scaling factors according to Tholes model.
           d1      = f_sqrt(d2)
           d2      = d2**(-1)

           ralpha  = aewald * d1
           d1      = d1**(-1)

           if (do_correct) then
              bn1     = 0.0
              bn2     = 0.0

              if( sdamp == 0.0 ) then
                sc3   = (1.0 - f_exp(-100.0)) * uscale
                sc5   = (1.0 - f_exp(-100.0)*(1-sdamp2)) * uscale
              else
                sdamp1= - pgamma / ((d1*sdamp)*(d1*sdamp)*(d1*sdamp))
                if (sdamp1 > -50.0) then
                  expdamp1 = f_exp(sdamp1)
                  sc3 = (1.0- expdamp1) * uscale
                  sc5 = (1.0- expdamp1*(1.0-sdamp1)) * uscale
                else
                  sc3 = uscale
                  sc5 = uscale
                end if
              end if
           else
              exp2a   = f_exp(-ralpha*ralpha)
              !call erfcore_inl(ralpha,bn0,1)
              bn0     = f_erfc(ralpha)   ! Attention to Macro f_erfc in single prec with hastings method (check tinker_cudart.h)
              bn0     = bn0*d1
              bn1     = (   bn0 +       alsq2*alsq2n*exp2a)*d2
              bn2     = ( 3*bn1 + alsq2*alsq2*alsq2n*exp2a)*d2

              if( sdamp == 0.0 ) then
                sc3   = (1.0 - expdamp2)
                sc5   = (1.0 - expdamp2 * (1 - sdamp2))
              else
                sdamp1= - pgamma / ((d1*sdamp)*(d1*sdamp)*(d1*sdamp))
                if (sdamp1 > -50.0) then
                  expdamp1 = f_exp(sdamp1)
                  sc3 = (1.0 - expdamp1)
                  sc5 = (1.0 - expdamp1 * (1 - sdamp1))
                else
                  sc3 = 1.0
                  sc5 = 1.0
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
     &              duo_tmatxb(d2,dist,dpui,dpuk,
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
        type(real6),intent(in)   ::dpui
        type(real7),intent(in)   ::dpuk
        type(real3),intent(in)   ::dist
        real(t_p)  ,intent(in)   ::sdamp,pgamma,aewald,alsq2,alsq2n
        real(t_p)  ,intent(in)   ::uscale
        real(t_p)  ,intent(inout)::d2
        type(real3),intent(inout)::fid,fip,fkd,fkp
        logical    ,intent(in)   :: do_correct
        real(t_p) d1,ralpha,exp2a
        real(t_p) bn0,bn1,bn2
        real(t_p) sdamp1,sdamp2,expdamp1,expdamp2
        real(t_p) rr3,rr5,rr3_bn1,rr5_bn2,sc3,sc5
        real(t_p) duir,dukr,puir,pukr  ! Scalar products duir = (du(i).r)
        parameter(sdamp2=-100.0,expdamp2= 0.0)

        ! compute the distances and the scaling factors according to Thole's model.
        d1      = f_sqrt(d2)
        d2      = d2**(-1)

        ralpha  = aewald * d1
        d1      = d1**(-1)

        if (do_correct) then
           bn1     = 0.0
           bn2     = 0.0

           if( sdamp == 0.0 ) then
             sc3    = (1.0-f_exp(sdamp2)) * uscale
             sc5    = (1.0-f_exp(sdamp2) * (1-sdamp2)) * uscale
           else
             sdamp1 = - pgamma / ((d1*sdamp)*(d1*sdamp)*(d1*sdamp))
             if (sdamp1 > -50.0) then
               expdamp1 = f_exp(sdamp1)
               sc3      =  (1.0-expdamp1) * uscale
               sc5      =  (1.0-expdamp1*(1-sdamp1))* uscale
             else
               sc3  = uscale
               sc5  = uscale
             end if
           end if
        else
           exp2a    = f_exp(-ralpha*ralpha)
           bn0      = f_erfc(ralpha)   ! Attention to Macro f_erfc in single prec with hastings method (check tinker_cudart.h)
           bn0      = bn0 * d1
           bn1      = (   bn0 +       alsq2*alsq2n*exp2a) * d2
           bn2      = ( 3*bn1 + alsq2*alsq2*alsq2n*exp2a) * d2

           if( sdamp == 0.0 ) then
             sc3    = 1.0 - expdamp2
             sc5    = 1.0 - expdamp2*(1.0-sdamp2)
           else
             sdamp1 = - pgamma / ((d1*sdamp)*(d1*sdamp)*(d1*sdamp))
             if (sdamp1 > -50.0) then
               expdamp1 = f_exp(sdamp1)
               sc3  = 1.0 - expdamp1
               sc5  = 1.0 - expdamp1 * (1.0 - sdamp1)
             else
               sc3  = 1.0
               sc5  = 1.0
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
     &          + dpui%zz * dist%z
        pukr    = dpuk%xx * dist%x + dpuk%yy * dist%y
     &          + dpuk%zz * dist%z

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
     &        duo_tmatxb1a(d2,p,dpui,dpuk,damp,pgamma,aewald,scal
     &           ,fid,fip,fkd,fkp,ver,fea)
!$acc routine
#if  defined(TINKER_CUF)
      use utilcu  ,only: f_erfc
#  if (defined(SINGLE)||defined(MIXED))
      use utilcu  ,only: f_sqrt,f_exp
#  endif
#endif
      use utilgpu ,only: real3,real6
      implicit none
      integer    ,intent(in) ::ver,fea
      type(real6),intent(in) ::dpui
      type(real7),intent(in) ::dpuk
      type(real3),intent(in) ::p
      real(t_p)  ,intent(in) ::d2,scal,aewald,damp,pgamma
      type(real3),intent(out)::fid,fip,fkd,fkp
      integer     sca
      real(t_p)   ralpha,exp2a
      real(t_p)   r,rr1,rr2,rr3,rr5
      real(t_p)   dlc(6),dmpe(3),dmpik(2)
      real(t_p)   dmp3,dmp5  ! Scalar products duir = (du(i).r)
      parameter(sca=__use_sca__)
c
c     compute the distances and the scaling factors according to
c     Thole's model.
c
      r   = f_sqrt(d2)
      rr1 = r**(-1)
      rr2 = rr1 * rr1
      rr3 = rr2 * rr1
      rr5 = 3.0 * rr2 * rr3
c
c     calculate real space Ewald error function damping
c
      IF (IAND(ver,sca).EQ.0) THEN
         call dampewald_inl (5,r,d2,aewald,1.0,dmpe)
      END IF
c
c     find the field components for Thole polarization damping
c
      call dampthole2_inl (damp,pgamma,5,r,dmpik)
      IF (IAND(fea,sca).NE.0) THEN
         dmp3 = scal*dmpik(1)*rr3
         dmp5 = scal*dmpik(2)*rr5
      ELSE
         dmp3 = dmpe(2) - (1.0-dmpik(1))*rr3
         dmp5 = dmpe(3) - (1.0-dmpik(2))*rr5
      END IF
      dlc(1) = -dmp3 + dmp5*p%x*p%x
      dlc(2) =  dmp5*p%x*p%y
      dlc(3) =  dmp5*p%x*p%z
      dlc(4) = -dmp3 + dmp5*p%y*p%y
      dlc(5) =  dmp5*p%y*p%z
      dlc(6) = -dmp3 + dmp5*p%z*p%z
c
c     find the field components for charge penetration damping
c

      fid%x = -dlc(1)*dpuk%x+dlc(2)*dpuk%y+dlc(3)*dpuk%z
      fid%y = -dlc(2)*dpuk%x+dlc(4)*dpuk%y+dlc(5)*dpuk%z
      fid%z = -dlc(3)*dpuk%x+dlc(5)*dpuk%y+dlc(6)*dpuk%z
      fkd%x = -dlc(1)*dpui%x+dlc(2)*dpui%y+dlc(3)*dpui%z
      fkd%y = -dlc(2)*dpui%x+dlc(4)*dpui%y+dlc(5)*dpui%z
      fkd%z = -dlc(3)*dpui%x+dlc(5)*dpui%y+dlc(6)*dpui%z

      fip%x = -dlc(1)*dpuk%xx+dlc(2)*dpuk%yy+dlc(3)*dpuk%zz
      fip%y = -dlc(2)*dpuk%xx+dlc(4)*dpuk%yy+dlc(5)*dpuk%zz
      fip%z = -dlc(3)*dpuk%xx+dlc(5)*dpuk%yy+dlc(6)*dpuk%zz
      fkp%x = -dlc(1)*dpui%xx+dlc(2)*dpui%yy+dlc(3)*dpui%zz
      fkp%y = -dlc(2)*dpui%xx+dlc(4)*dpui%yy+dlc(5)*dpui%zz
      fkp%z = -dlc(3)*dpui%xx+dlc(5)*dpui%yy+dlc(6)*dpui%zz

      end subroutine

      M_subroutine
     &        duo_tmatxb_cpen(d2,p,dpui,dpuk,pentyp
     &           ,alphai,alphak,aewald,scal
     &           ,fid,fip,fkd,fkp,ver,fea)
!$acc routine
#if  defined(TINKER_CUF)
      use utilcu  ,only: f_erfc
#  if (defined(SINGLE)||defined(MIXED))
      use utilcu  ,only: f_sqrt,f_exp
#  endif
#endif
      use utilgpu ,only: real3,real6
      implicit none
      integer    ,intent(in) ::pentyp,ver,fea
      type(real6),intent(in) ::dpui
      type(real7),intent(in) ::dpuk
      type(real3),intent(in) ::p
      real(t_p)  ,intent(in) ::alphai,alphak,aewald,scal,d2
      type(real3),intent(out)::fid,fip,fkd,fkp
      integer   sca
      real(t_p) ralpha,exp2a
      real(t_p) bn0,bn1,bn2
      real(t_p) sdamp1,expdamp1
      real(t_p) r,rr1,rr2,rr3,rr5
      real(t_p) dlc(6),dmpe(3),dmpik(2)
      real(t_p) duir,dukr,puir,pukr  ! Scalar products duir = (du(i).r)
      parameter(sca=__use_sca__)
c
c     compute the distances and the scaling factors according to
c     Thole's model.
c
      r   = f_sqrt(d2)
      rr1 = r**(-1)
      rr2 = rr1 * rr1
      rr3 = rr2 * rr1
      rr5 = 3.0 * rr2 * rr3
c
c     calculate real space Ewald error function damping
c
      IF (IAND(ver,sca).EQ.0) THEN
         call dampewald_inl (5,r,d2,aewald,1.0,dmpe)
      END IF
c
c     find the field components for Thole polarization damping
c
c     if (use_thole) then
c        call dampthole2_inl (damp,pgamma,5,r,dmpik)
c        IF (IAND(fea,sca).NE.0) THEN
c           dmp3 = scal*dmpik(1)*rr3
c           dmp5 = scal*dmpik(2)*rr5
c        ELSE
c           dmp3 = dmpe(2) - (1.0-dmpik(1))*rr3
c           dmp5 = dmpe(3) - (1.0-dmpik(2))*rr5
c        END IF
c        dlc(1) = -dmp3 + dmp5*p%x*p%x
c        dlc(2) = dmp5*p%x*p%y
c        dlc(3) = dmp5*p%x*p%z
c        dlc(4) = -dmp3 + dmp5*p%y*p%y
c        dlc(5) = dmp5*p%y*p%z
c        dlc(6) = -dmp3 + dmp5*p%z*p%z
c
c     find the field components for charge penetration damping
c
c     else if (use_chgpen) then
         call dampmut_inl (r,alphai,alphak,pentyp,dmpik)
         if (IAND(fea,sca).NE.0) THEN
            rr3 = scal*dmpik(1)*rr3
            rr5 = scal*dmpik(2)*rr5
         ELSE
            rr3 = dmpe(2) - (1.0-dmpik(1))*rr3
            rr5 = dmpe(3) - (1.0-dmpik(2))*rr5
         END IF
         dlc(1) = -rr3 + rr5*p%x*p%x
         dlc(2) =        rr5*p%x*p%y
         dlc(3) =        rr5*p%x*p%z
         dlc(4) = -rr3 + rr5*p%y*p%y
         dlc(5) =        rr5*p%y*p%z
         dlc(6) = -rr3 + rr5*p%z*p%z
c     end if

      fid%x = -dlc(1)*dpuk%x+dlc(2)*dpuk%y+dlc(3)*dpuk%z
      fid%y = -dlc(2)*dpuk%x+dlc(4)*dpuk%y+dlc(5)*dpuk%z
      fid%z = -dlc(3)*dpuk%x+dlc(5)*dpuk%y+dlc(6)*dpuk%z
      fkd%x = -dlc(1)*dpui%x+dlc(2)*dpui%y+dlc(3)*dpui%z
      fkd%y = -dlc(2)*dpui%x+dlc(4)*dpui%y+dlc(5)*dpui%z
      fkd%z = -dlc(3)*dpui%x+dlc(5)*dpui%y+dlc(6)*dpui%z

      fip%x = -dlc(1)*dpuk%xx+dlc(2)*dpuk%yy+dlc(3)*dpuk%zz
      fip%y = -dlc(2)*dpuk%xx+dlc(4)*dpuk%yy+dlc(5)*dpuk%zz
      fip%z = -dlc(3)*dpuk%xx+dlc(5)*dpuk%yy+dlc(6)*dpuk%zz
      fkp%x = -dlc(1)*dpui%xx+dlc(2)*dpui%yy+dlc(3)*dpui%zz
      fkp%y = -dlc(2)*dpui%xx+dlc(4)*dpui%yy+dlc(5)*dpui%zz
      fkp%z = -dlc(3)*dpui%xx+dlc(5)*dpui%yy+dlc(6)*dpui%zz

      end subroutine

      ! Increment output
      M_subroutine
     &        duo_tmatxb_cpen2(d2,p,dpui,dpuk,pentyp
     &           ,alphai,alphak,aewald,scal
     &           ,fid,fip,fkd,fkp,ver,fea)
!$acc routine
#if  defined(TINKER_CUF)
      use utilcu  ,only: f_erfc
#  if (defined(SINGLE)||defined(MIXED))
      use utilcu  ,only: f_sqrt,f_exp
#  endif
#endif
      use utilgpu ,only: real3,real6
      implicit none
      integer    ,intent(in)   ::pentyp,ver,fea
      type(real6),intent(in)   ::dpui
      type(real7),intent(in)   ::dpuk
      type(real3),intent(in)   ::p
      real(t_p)  ,intent(in)   ::alphai,alphak,aewald,scal,d2
      type(real3),intent(inout)::fid,fip,fkd,fkp
      integer   sca
      real(t_p) ralpha,exp2a
      real(t_p) bn0,bn1,bn2
      real(t_p) sdamp1,expdamp1
      real(t_p) r,rr1,rr2,rr3,rr5
      real(t_p) dlc(6),dmpe(3),dmpik(2)
      real(t_p) duir,dukr,puir,pukr  ! Scalar products duir = (du(i).r)
      parameter(sca=__use_sca__)
c
c     compute the distances and the scaling factors according to
c     Thole's model.
c
      r   = f_sqrt(d2)
      rr1 = r**(-1)
      rr2 = rr1 * rr1
      rr3 = rr2 * rr1
      rr5 = 3.0 * rr2 * rr3
c
c     calculate real space Ewald error function damping
c
      IF (IAND(ver,sca).EQ.0) THEN
         call dampewald_inl (5,r,d2,aewald,1.0,dmpe)
      END IF
c
c     find the field components for Thole polarization damping
c
c     if (use_thole) then
c        call dampthole2_inl (damp,pgamma,5,r,dmpik)
c        IF (IAND(fea,sca).NE.0) THEN
c           dmp3 = scal*dmpik(1)*rr3
c           dmp5 = scal*dmpik(2)*rr5
c        ELSE
c           dmp3 = dmpe(2) - (1.0-dmpik(1))*rr3
c           dmp5 = dmpe(3) - (1.0-dmpik(2))*rr5
c        END IF
c        dlc(1) = -dmp3 + dmp5*p%x*p%x
c        dlc(2) = dmp5*p%x*p%y
c        dlc(3) = dmp5*p%x*p%z
c        dlc(4) = -dmp3 + dmp5*p%y*p%y
c        dlc(5) = dmp5*p%y*p%z
c        dlc(6) = -dmp3 + dmp5*p%z*p%z
c
c     find the field components for charge penetration damping
c
c     else if (use_chgpen) then
         call dampmut_inl (r,alphai,alphak,pentyp,dmpik)
         IF (IAND(fea,sca).NE.0) THEN
            rr3 = scal*dmpik(1)*rr3
            rr5 = scal*dmpik(2)*rr5
         ELSE
            rr3 = dmpe(2) - (1.0-dmpik(1))*rr3
            rr5 = dmpe(3) - (1.0-dmpik(2))*rr5
         END IF
         dlc(1) = -rr3 + rr5*p%x*p%x
         dlc(2) =        rr5*p%x*p%y
         dlc(3) =        rr5*p%x*p%z
         dlc(4) = -rr3 + rr5*p%y*p%y
         dlc(5) =        rr5*p%y*p%z
         dlc(6) = -rr3 + rr5*p%z*p%z
c     end if

      fid%x = fid%x - dlc(1)*dpuk%x+dlc(2)*dpuk%y+dlc(3)*dpuk%z
      fid%y = fid%y - dlc(2)*dpuk%x+dlc(4)*dpuk%y+dlc(5)*dpuk%z
      fid%z = fid%z - dlc(3)*dpuk%x+dlc(5)*dpuk%y+dlc(6)*dpuk%z
      fkd%x = fkd%x - dlc(1)*dpui%x+dlc(2)*dpui%y+dlc(3)*dpui%z
      fkd%y = fkd%y - dlc(2)*dpui%x+dlc(4)*dpui%y+dlc(5)*dpui%z
      fkd%z = fkd%z - dlc(3)*dpui%x+dlc(5)*dpui%y+dlc(6)*dpui%z

      fip%x = fip%x - dlc(1)*dpuk%xx+dlc(2)*dpuk%yy+dlc(3)*dpuk%zz
      fip%y = fip%y - dlc(2)*dpuk%xx+dlc(4)*dpuk%yy+dlc(5)*dpuk%zz
      fip%z = fip%z - dlc(3)*dpuk%xx+dlc(5)*dpuk%yy+dlc(6)*dpuk%zz
      fkp%x = fkp%x - dlc(1)*dpui%xx+dlc(2)*dpui%yy+dlc(3)*dpui%zz
      fkp%y = fkp%y - dlc(2)*dpui%xx+dlc(4)*dpui%yy+dlc(5)*dpui%zz
      fkp%z = fkp%z - dlc(3)*dpui%xx+dlc(5)*dpui%yy+dlc(6)*dpui%zz

      end subroutine

        M_subroutine
     &              tmatxb4_couple(d2,dist,dpui,dpuk,equi,equk,
     &                         sdamp,pgamma,aewald,alsq2,alsq2n,uscale,
     &                             fid,fip,fkd,fkp,
     &                             fie,fiq,fke,fkq,do_correct)
!$acc routine
#if  defined(TINKER_CUF)
           use utilcu  ,only: f_erfc
#  if (defined(SINGLE)||defined(MIXED))
           use utilcu  ,only: f_sqrt,f_exp
#  endif
#endif
           use utilgpu ,only: real3,real6
           implicit none
           type(real6),intent(in)   ::dpui,dpuk,equi,equk
           type(real3),intent(in)   ::dist
           real(t_p)  ,intent(in)   ::sdamp,pgamma,aewald,alsq2,alsq2n
           real(t_p)  ,intent(in)   ::uscale
           real(t_p)  ,intent(inout)::d2
           type(real3),intent(inout)::fid,fip,fkd,fkp,fie,fiq,fke,fkq
           logical    ,intent(in)   :: do_correct
           real(t_p) d1,ralpha,exp2a
           real(t_p) bn0,bn1,bn2
           real(t_p) sdamp1,expdamp1
           real(t_p) rr3,rr5,rr3_bn1,rr5_bn2,sc3,sc5
           real(t_p) duir,dukr,puir,pukr  ! Scalar products duir = (du(i).r)

           ! compute the distances and the scaling factors according to Thole s model.
           d1      = f_sqrt(d2)
           d2      = d2**(-1)

           ralpha  = aewald * d1
           d1      = d1**(-1)

           if (do_correct) then
              bn1     = 0.0
              bn2     = 0.0

              if( sdamp == 0.0 ) then
                sdamp1  = -100.0
                sc3    = (1-f_exp(sdamp1)) * uscale
                sc5    = (1-f_exp(sdamp1) * (1 - sdamp1)) * uscale
              else
                sdamp1 = - pgamma / ((d1*sdamp)*(d1*sdamp)*(d1*sdamp))
                if (sdamp1 > -50.0) then
                  expdamp1 = f_exp(sdamp1)
                  sc3      =  (1- expdamp1) * uscale
                  sc5      =   (1-expdamp1* (1 - sdamp1)) * uscale 
                else
                  sc3     = uscale
                  sc5     = uscale 
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
                sc3      =   (1 - expdamp1) * uscale
                sc5      =   (1 - expdamp1*(1 - sdamp1)) * uscale 
              else
                sdamp1 = - pgamma / ((d1*sdamp)*(d1*sdamp)*(d1*sdamp))
                if (sdamp1 > -50.0) then
                  expdamp1 = f_exp(sdamp1)
                  sc3      =   (1 - expdamp1) * uscale
                  sc5      =   (1 - expdamp1*(1 - sdamp1)) * uscale 
                else
                  sc3     = uscale
                  sc5     = uscale
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

           duir    = equi%x * dist%x + equi%y * dist%y + equi%z *dist%z
           dukr    = equk%x * dist%x + equk%y * dist%y + equk%z *dist%z

           puir    = equi%xx * dist%x + equi%yy * dist%y
     &             + equi%zz * dist%z
           pukr    = equk%xx * dist%x + equk%yy * dist%y
     &             + equk%zz * dist%z

        fie%x  = WRITE_C(fie%x) - rr3_bn1*equk%x + rr5_bn2*dukr*dist%x
        fie%y  = WRITE_C(fie%y) - rr3_bn1*equk%y + rr5_bn2*dukr*dist%y
        fie%z  = WRITE_C(fie%z) - rr3_bn1*equk%z + rr5_bn2*dukr*dist%z

        fiq%x  = WRITE_C(fiq%x) - rr3_bn1*equk%xx + rr5_bn2*pukr*dist%x
        fiq%y  = WRITE_C(fiq%y) - rr3_bn1*equk%yy + rr5_bn2*pukr*dist%y
        fiq%z  = WRITE_C(fiq%z) - rr3_bn1*equk%zz + rr5_bn2*pukr*dist%z

        fke%x  = WRITE_C(fke%x) - rr3_bn1*equi%x + rr5_bn2*duir*dist%x
        fke%y  = WRITE_C(fke%y) - rr3_bn1*equi%y + rr5_bn2*duir*dist%y
        fke%z  = WRITE_C(fke%z) - rr3_bn1*equi%z + rr5_bn2*duir*dist%z

        fkq%x  = WRITE_C(fkq%x) - rr3_bn1*equi%xx + rr5_bn2*puir*dist%x
        fkq%y  = WRITE_C(fkq%y) - rr3_bn1*equi%yy + rr5_bn2*puir*dist%y
        fkq%z  = WRITE_C(fkq%z) - rr3_bn1*equi%zz + rr5_bn2*puir*dist%z
        end subroutine
#endif
