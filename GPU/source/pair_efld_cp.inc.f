#ifndef PAIR_EFLD_CP
#define PAIR_EFLD_CP
#include "tinker_cudart.h"
#include "damping.inc.f"

      M_subroutine
     &        duo_efld0Cpen(mi,mk,r2,p,pscale,dscale,aewald,pentyp
     &           ,corei,corek,vali,valk,alphai,alphak
     &           ,fid,fkd,fip,fkp,tver,tfea
     &           )
      use tintypes ,only: rpole_elt
#if defined(TINKER_CUF) && (TINKER_SINGLE_PREC+TINKER_MIXED_PREC)
      use utilcu   ,only: f_sqrt
#endif
      integer  ,intent(in):: pentyp,tver,tfea
      real(t_p),intent(in):: r2,pscale,dscale,aewald,corei,corek
     &         ,vali,valk,alphai,alphak
      type(real3),intent(in):: p
      type(rpole_elt),intent(in):: mi,mk
      real(t_p),intent(out):: fip(3),fkp(3),fid(3),fkd(3)
      integer   sca
      real(t_p) r,rr1,rr2,rr3,rr5,rr7
      real(t_p) dir,qix,qiy,qiz,dkr,qkx,qky,qkz,qir,qkr,dmp3,dmp5,dmp7
      real(t_p) rr3i,rr5i,rr7i,rr3k,rr5k,rr7k
      real(t_p) dmpe(4),dmpi(3),dmpk(3),dmpik(3)
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
#if 0
      if (use_thole) then
         call dampthole_inl (7,r,damp,pgamma,dmpik)
         IF (IAND(tver,sca).NE.0) THEN
            dmp3 = dscale*dmpik(1)*rr3
            dmp5 = dscale*dmpik(2)*rr5
            dmp7 = dscale*dmpik(3)*rr7
         ELSE
            dmp3 = dmpe(2) - (1.0-dmpik(1))*rr3
            dmp5 = dmpe(3) - (1.0-dmpik(2))*rr5
            dmp7 = dmpe(4) - (1.0-dmpik(3))*rr7
         END IF
         fid(1) = -p%x*(dmp3*mk%c-dmp5*dkr+dmp7*qkr)
     &               - dmp3*mk%dx + 2.0*dmp5*qkx
         fid(2) = -p%y*(dmp3*mk%c-dmp5*dkr+dmp7*qkr)
     &               - dmp3*mk%dy + 2.0*dmp5*qky
         fid(3) = -p%z*(dmp3*mk%c-dmp5*dkr+dmp7*qkr)
     &               - dmp3*mk%dz + 2.0*dmp5*qkz
         fkd(1) = p%x*(dmp3*mi%c+dmp5*dir+dmp7*qir)
     &               - dmp3*mi%dx - 2.0*dmp5*qix
         fkd(2) = p%y*(dmp3*mi%c+dmp5*dir+dmp7*qir)
     &               - dmp3*mi%dy - 2.0*dmp5*qiy
         fkd(3) = p%z*(dmp3*mi%c+dmp5*dir+dmp7*qir)
     &               - dmp3*mi%dz - 2.0*dmp5*qiz
         IF (IAND(tver,sca).NE.0) THEN
            dmp3 = pscale*dmpik(3)*rr3
            dmp5 = pscale*dmpik(5)*rr5
            dmp7 = pscale*dmpik(7)*rr7
         ELSE
            dmp3 = dmpe(2) - (1.0-dmpik(1))*rr3
            dmp5 = dmpe(3) - (1.0-dmpik(2))*rr5
            dmp7 = dmpe(4) - (1.0-dmpik(3))*rr7
         END IF
         fip(1) = -p%x*(dmp3*mk%c-dmp5*dkr+dmp7*qkr)
     &               - dmp3*mk%dx + 2.0*dmp5*qkx
         fip(2) = -p%y*(dmp3*mk%c-dmp5*dkr+dmp7*qkr)
     &               - dmp3*mk%dy + 2.0*dmp5*qky
         fip(3) = -p%z*(dmp3*mk%c-dmp5*dkr+dmp7*qkr)
     &               - dmp3*mk%dz + 2.0*dmp5*qkz
         fkp(1) = p%x*(dmp3*mi%c+dmp5*dir+dmp7*qir)
     &               - dmp3*mi%dx - 2.0*dmp5*qix
         fkp(2) = p%y*(dmp3*mi%c+dmp5*dir+dmp7*qir)
     &               - dmp3*mi%dy - 2.0*dmp5*qiy
         fkp(3) = p%z*(dmp3*mi%c+dmp5*dir+dmp7*qir)
     &               - dmp3*mi%dz - 2.0*dmp5*qiz
c
c     find the field components for charge penetration damping
c
      else if (use_chgpen) then
#endif
         call dampdir_inl (pentyp,r,alphai,alphak,dmpi,dmpk)
         IF (IAND(tver,sca).NE.0) THEN
            rr3i = dscale*dmpi(1)*rr3
            rr5i = dscale*dmpi(2)*rr5
            rr7i = dscale*dmpi(3)*rr7
            rr3k = dscale*dmpk(1)*rr3
            rr5k = dscale*dmpk(2)*rr5
            rr7k = dscale*dmpk(3)*rr7
            rr3  = dscale*rr3
         ELSE
            rr3i = dmpe(2) - (1.0-dmpi(1))*rr3
            rr5i = dmpe(3) - (1.0-dmpi(2))*rr5
            rr7i = dmpe(4) - (1.0-dmpi(3))*rr7
            rr3k = dmpe(2) - (1.0-dmpk(1))*rr3
            rr5k = dmpe(3) - (1.0-dmpk(2))*rr5
            rr7k = dmpe(4) - (1.0-dmpk(3))*rr7
            rr3  = dmpe(2)
         END IF
         fid(1) =-p%x*(rr3*corek + rr3k*valk
     &               - rr5k*dkr + rr7k*qkr)
     &               - rr3k*mk%dx + 2.0*rr5k*qkx
         fid(2) =-p%y*(rr3*corek + rr3k*valk
     &               - rr5k*dkr + rr7k*qkr)
     &               - rr3k*mk%dy + 2.0*rr5k*qky
         fid(3) =-p%z*(rr3*corek + rr3k*valk
     &               - rr5k*dkr + rr7k*qkr)
     &               - rr3k*mk%dz + 2.0*rr5k*qkz
         fkd(1) = p%x*(rr3*corei + rr3i*vali
     &               + rr5i*dir + rr7i*qir)
     &               - rr3i*mi%dx - 2.0*rr5i*qix
         fkd(2) = p%y*(rr3*corei + rr3i*vali
     &               + rr5i*dir + rr7i*qir)
     &               - rr3i*mi%dy - 2.0*rr5i*qiy
         fkd(3) = p%z*(rr3*corei + rr3i*vali
     &               + rr5i*dir + rr7i*qir)
     &               - rr3i*mi%dz - 2.0*rr5i*qiz
         IF (IAND(tver,sca).NE.0) THEN
            rr3i = pscale*dmpi(1)*rr3
            rr5i = pscale*dmpi(2)*rr5
            rr7i = pscale*dmpi(3)*rr7
            rr3k = pscale*dmpk(1)*rr3
            rr5k = pscale*dmpk(2)*rr5
            rr7k = pscale*dmpk(3)*rr7
            rr3  = pscale*rr3
         ELSE
            rr3i = dmpe(2) - (1.0-dmpi(1))*rr3
            rr5i = dmpe(3) - (1.0-dmpi(2))*rr5
            rr7i = dmpe(4) - (1.0-dmpi(3))*rr7
            rr3k = dmpe(2) - (1.0-dmpk(1))*rr3
            rr5k = dmpe(3) - (1.0-dmpk(2))*rr5
            rr7k = dmpe(4) - (1.0-dmpk(3))*rr7
            rr3  = dmpe(2)
         END IF
         fip(1) =-p%x*(rr3*corek + rr3k*valk
     &               - rr5k*dkr + rr7k*qkr)
     &               - rr3k*mk%dx + 2.0*rr5k*qkx
         fip(2) =-p%y*(rr3*corek + rr3k*valk
     &               - rr5k*dkr + rr7k*qkr)
     &               - rr3k*mk%dy + 2.0*rr5k*qky
         fip(3) =-p%z*(rr3*corek + rr3k*valk
     &               - rr5k*dkr + rr7k*qkr)
     &               - rr3k*mk%dz + 2.0*rr5k*qkz
         fkp(1) = p%x*(rr3*corei + rr3i*vali
     &               + rr5i*dir + rr7i*qir)
     &               - rr3i*mi%dx - 2.0*rr5i*qix
         fkp(2) = p%y*(rr3*corei + rr3i*vali
     &               + rr5i*dir + rr7i*qir)
     &               - rr3i*mi%dy - 2.0*rr5i*qiy
         fkp(3) = p%z*(rr3*corei + rr3i*vali
     &               + rr5i*dir + rr7i*qir)
     &               - rr3i*mi%dz - 2.0*rr5i*qiz
#if 0
      end if
#endif
      end subroutine

      M_subroutine
     &        duo_efld0CpenNEwd(mi,mk,r2,p,pscale,dscale,pentyp
     &           ,corei,corek,vali,valk,alphai,alphak
     &           ,fid,fkd,fip,fkp,tver,tfea
     &           )
      use tintypes ,only: rpole_elt
#if defined(TINKER_CUF) && (TINKER_SINGLE_PREC+TINKER_MIXED_PREC)
      use utilcu   ,only: f_sqrt
#endif
      integer  ,intent(in):: pentyp,tver,tfea
      real(t_p),intent(in):: r2,pscale,dscale,corei,corek
     &         ,vali,valk,alphai,alphak
      type(real3),intent(in):: p
      type(rpole_elt),intent(in):: mi,mk
      real(t_p),intent(out):: fip(3),fkp(3),fid(3),fkd(3)
      integer   sca
      real(t_p) r,rr1,rr2,rr3,rr5,rr7
      real(t_p) dir,qix,qiy,qiz,dkr,qkx,qky,qkz,qir,qkr,dmp3,dmp5,dmp7
      real(t_p) rr3i,rr5i,rr7i,rr3k,rr5k,rr7k
      real(t_p) dmpe(4),dmpi(3),dmpk(3),dmpik(3)
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
      call dampdir_inl (pentyp,r,alphai,alphak,dmpi,dmpk)
      IF (IAND(tver,sca).NE.0) THEN
         rr3i = dscale*dmpi(1)*rr3
         rr5i = dscale*dmpi(2)*rr5
         rr7i = dscale*dmpi(3)*rr7
         rr3k = dscale*dmpk(1)*rr3
         rr5k = dscale*dmpk(2)*rr5
         rr7k = dscale*dmpk(3)*rr7
         rr3  = dscale*rr3
      ELSE
         rr3i = (dmpi(1)-1.0)*rr3
         rr5i = (dmpi(2)-1.0)*rr5
         rr7i = (dmpi(3)-1.0)*rr7
         rr3k = (dmpk(1)-1.0)*rr3
         rr5k = (dmpk(2)-1.0)*rr5
         rr7k = (dmpk(3)-1.0)*rr7
         rr3  = 0.0
      END IF
      fid(1) =-p%x*(rr3*corek + rr3k*valk
     &            - rr5k*dkr + rr7k*qkr)
     &            - rr3k*mk%dx + 2.0*rr5k*qkx
      fid(2) =-p%y*(rr3*corek + rr3k*valk
     &            - rr5k*dkr + rr7k*qkr)
     &            - rr3k*mk%dy + 2.0*rr5k*qky
      fid(3) =-p%z*(rr3*corek + rr3k*valk
     &            - rr5k*dkr + rr7k*qkr)
     &            - rr3k*mk%dz + 2.0*rr5k*qkz
      fkd(1) = p%x*(rr3*corei + rr3i*vali
     &            + rr5i*dir + rr7i*qir)
     &            - rr3i*mi%dx - 2.0*rr5i*qix
      fkd(2) = p%y*(rr3*corei + rr3i*vali
     &            + rr5i*dir + rr7i*qir)
     &            - rr3i*mi%dy - 2.0*rr5i*qiy
      fkd(3) = p%z*(rr3*corei + rr3i*vali
     &            + rr5i*dir + rr7i*qir)
     &            - rr3i*mi%dz - 2.0*rr5i*qiz
      IF (IAND(tver,sca).NE.0) THEN
         rr3i = pscale*dmpi(1)*rr3
         rr5i = pscale*dmpi(2)*rr5
         rr7i = pscale*dmpi(3)*rr7
         rr3k = pscale*dmpk(1)*rr3
         rr5k = pscale*dmpk(2)*rr5
         rr7k = pscale*dmpk(3)*rr7
         rr3  = pscale*rr3
      END IF
      fip(1) =-p%x*(rr3*corek + rr3k*valk
     &            - rr5k*dkr + rr7k*qkr)
     &            - rr3k*mk%dx + 2.0*rr5k*qkx
      fip(2) =-p%y*(rr3*corek + rr3k*valk
     &            - rr5k*dkr + rr7k*qkr)
     &            - rr3k*mk%dy + 2.0*rr5k*qky
      fip(3) =-p%z*(rr3*corek + rr3k*valk
     &            - rr5k*dkr + rr7k*qkr)
     &            - rr3k*mk%dz + 2.0*rr5k*qkz
      fkp(1) = p%x*(rr3*corei + rr3i*vali
     &            + rr5i*dir + rr7i*qir)
     &            - rr3i*mi%dx - 2.0*rr5i*qix
      fkp(2) = p%y*(rr3*corei + rr3i*vali
     &            + rr5i*dir + rr7i*qir)
     &            - rr3i*mi%dy - 2.0*rr5i*qiy
      fkp(3) = p%z*(rr3*corei + rr3i*vali
     &            + rr5i*dir + rr7i*qir)
     &            - rr3i*mi%dz - 2.0*rr5i*qiz
      end subroutine

      ! Increment output
      M_subroutine
     &        duo_efld0Cpen2(mi,mk,r2,p,pscale,dscale,aewald,pentyp
     &           ,corei,corek,vali,valk,alphai,alphak
     &           ,fid,fkd,fip,fkp,tver,tfea
     &           )
      use tintypes ,only: rpole_elt
#if defined(TINKER_CUF) && (TINKER_SINGLE_PREC+TINKER_MIXED_PREC)
      use utilcu   ,only: f_sqrt
#endif
      integer  ,intent(in):: pentyp,tver,tfea
      real(t_p),intent(in):: r2,pscale,dscale,aewald,corei,corek
     &         ,vali,valk,alphai,alphak
      type(real3),intent(in):: p
      type(rpole_elt),intent(in):: mi,mk
      real(t_p),intent(inout):: fip(3),fkp(3),fid(3),fkd(3)
      integer   sca
      real(t_p) r,rr1,rr2,rr3,rr5,rr7
      real(t_p) dir,qix,qiy,qiz,dkr,qkx,qky,qkz,qir,qkr,dmp3,dmp5,dmp7
      real(t_p) rr3i,rr5i,rr7i,rr3k,rr5k,rr7k
      real(t_p) dmpe(4),dmpi(3),dmpk(3),dmpik(3)
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
#if 0
      if (use_thole) then
         call dampthole_inl (7,r,damp,pgamma,dmpik)
         IF (IAND(tver,sca).NE.0) THEN
            dmp3 = dscale*dmpik(3)*rr3
            dmp5 = dscale*dmpik(5)*rr5
            dmp7 = dscale*dmpik(7)*rr7
         ELSE
            dmp3 = dmpe(2) - (1.0-dmpik(3))*rr3
            dmp5 = dmpe(3) - (1.0-dmpik(5))*rr5
            dmp7 = dmpe(4) - (1.0-dmpik(7))*rr7
         END IF
         fid(1) = -p%x*(dmp3*mk%c-dmp5*dkr+dmp7*qkr)
     &               - dmp3*mk%dx + 2.0*dmp5*qkx
         fid(2) = -p%y*(dmp3*mk%c-dmp5*dkr+dmp7*qkr)
     &               - dmp3*mk%dy + 2.0*dmp5*qky
         fid(3) = -p%z*(dmp3*mk%c-dmp5*dkr+dmp7*qkr)
     &               - dmp3*mk%dz + 2.0*dmp5*qkz
         fkd(1) = p%x*(dmp3*mi%c+dmp5*dir+dmp7*qir)
     &               - dmp3*mi%dx - 2.0*dmp5*qix
         fkd(2) = p%y*(dmp3*mi%c+dmp5*dir+dmp7*qir)
     &               - dmp3*mi%dy - 2.0*dmp5*qiy
         fkd(3) = p%z*(dmp3*mi%c+dmp5*dir+dmp7*qir)
     &               - dmp3*mi%dz - 2.0*dmp5*qiz
         IF (IAND(tver,sca).NE.0) THEN
            dmp3 = pscale*dmpik(3)*rr3
            dmp5 = pscale*dmpik(5)*rr5
            dmp7 = pscale*dmpik(7)*rr7
         ELSE
            dmp3 = dmpe(2) - (1.0-dmpik(3))*rr3
            dmp5 = dmpe(3) - (1.0-dmpik(5))*rr5
            dmp7 = dmpe(4) - (1.0-dmpik(7))*rr7
         END IF
         fip(1) = -p%x*(dmp3*mk%c-dmp5*dkr+dmp7*qkr)
     &               - dmp3*mk%dx + 2.0*dmp5*qkx
         fip(2) = -p%y*(dmp3*mk%c-dmp5*dkr+dmp7*qkr)
     &               - dmp3*mk%dy + 2.0*dmp5*qky
         fip(3) = -p%z*(dmp3*mk%c-dmp5*dkr+dmp7*qkr)
     &               - dmp3*mk%dz + 2.0*dmp5*qkz
         fkp(1) = p%x*(dmp3*mi%c+dmp5*dir+dmp7*qir)
     &               - dmp3*mi%dx - 2.0*dmp5*qix
         fkp(2) = p%y*(dmp3*mi%c+dmp5*dir+dmp7*qir)
     &               - dmp3*mi%dy - 2.0*dmp5*qiy
         fkp(3) = p%z*(dmp3*mi%c+dmp5*dir+dmp7*qir)
     &               - dmp3*mi%dz - 2.0*dmp5*qiz
c
c     find the field components for charge penetration damping
c
      else if (use_chgpen) then
#endif
         call dampdir_inl (pentyp,r,alphai,alphak,dmpi,dmpk)
         IF (IAND(tver,sca).NE.0) THEN
            rr3i = dscale*dmpi(1)*rr3
            rr5i = dscale*dmpi(2)*rr5
            rr7i = dscale*dmpi(3)*rr7
            rr3k = dscale*dmpk(1)*rr3
            rr5k = dscale*dmpk(2)*rr5
            rr7k = dscale*dmpk(3)*rr7
            rr3  = dscale*rr3
         ELSE
            rr3i = dmpe(2) - (1.0-dmpi(1))*rr3
            rr5i = dmpe(3) - (1.0-dmpi(2))*rr5
            rr7i = dmpe(4) - (1.0-dmpi(3))*rr7
            rr3k = dmpe(2) - (1.0-dmpk(1))*rr3
            rr5k = dmpe(3) - (1.0-dmpk(2))*rr5
            rr7k = dmpe(4) - (1.0-dmpk(3))*rr7
            rr3  = dmpe(2)
         END IF
         fid(1) = fid(1) -p%x*(rr3*corek + rr3k*valk
     &               - rr5k*dkr + rr7k*qkr)
     &               - rr3k*mk%dx + 2.0*rr5k*qkx
         fid(2) = fid(2) -p%y*(rr3*corek + rr3k*valk
     &               - rr5k*dkr + rr7k*qkr)
     &               - rr3k*mk%dy + 2.0*rr5k*qky
         fid(3) = fid(3) -p%z*(rr3*corek + rr3k*valk
     &               - rr5k*dkr + rr7k*qkr)
     &               - rr3k*mk%dz + 2.0*rr5k*qkz
         fkd(1) = fkd(3) +p%x*(rr3*corei + rr3i*vali
     &               + rr5i*dir + rr7i*qir)
     &               - rr3i*mi%dx - 2.0*rr5i*qix
         fkd(2) = fkd(2) +p%y*(rr3*corei + rr3i*vali
     &               + rr5i*dir + rr7i*qir)
     &               - rr3i*mi%dy - 2.0*rr5i*qiy
         fkd(3) = fkd(3) +p%z*(rr3*corei + rr3i*vali
     &               + rr5i*dir + rr7i*qir)
     &               - rr3i*mi%dz - 2.0*rr5i*qiz
         IF (IAND(tver,sca).NE.0) THEN
            rr3i = pscale*dmpi(1)*rr3
            rr5i = pscale*dmpi(2)*rr5
            rr7i = pscale*dmpi(3)*rr7
            rr3k = pscale*dmpk(1)*rr3
            rr5k = pscale*dmpk(2)*rr5
            rr7k = pscale*dmpk(3)*rr7
            rr3  = pscale*rr3
         END IF
         fip(1) = fip(1) -p%x*(rr3*corek + rr3k*valk
     &               - rr5k*dkr + rr7k*qkr)
     &               - rr3k*mk%dx + 2.0*rr5k*qkx
         fip(2) = fip(2) -p%y*(rr3*corek + rr3k*valk
     &               - rr5k*dkr + rr7k*qkr)
     &               - rr3k*mk%dy + 2.0*rr5k*qky
         fip(3) = fip(3) -p%z*(rr3*corek + rr3k*valk
     &               - rr5k*dkr + rr7k*qkr)
     &               - rr3k*mk%dz + 2.0*rr5k*qkz
         fkp(1) = fkp(1) +p%x*(rr3*corei + rr3i*vali
     &               + rr5i*dir + rr7i*qir)
     &               - rr3i*mi%dx - 2.0*rr5i*qix
         fkp(2) = fkp(2) +p%y*(rr3*corei + rr3i*vali
     &               + rr5i*dir + rr7i*qir)
     &               - rr3i*mi%dy - 2.0*rr5i*qiy
         fkp(3) = fkp(3) +p%z*(rr3*corei + rr3i*vali
     &               + rr5i*dir + rr7i*qir)
     &               - rr3i*mi%dz - 2.0*rr5i*qiz
#if 0
      end if
#endif
      end subroutine

      M_subroutine
     &        duo_efld0CpenNEwd2(mi,mk,r2,p,pscale,dscale,pentyp
     &           ,corei,corek,vali,valk,alphai,alphak
     &           ,fid,fkd,fip,fkp,tver,tfea
     &           )
      use tintypes ,only: rpole_elt
#if defined(TINKER_CUF) && (TINKER_SINGLE_PREC+TINKER_MIXED_PREC)
      use utilcu   ,only: f_sqrt
#endif
      integer  ,intent(in):: pentyp,tver,tfea
      real(t_p),intent(in):: r2,pscale,dscale,corei,corek
     &         ,vali,valk,alphai,alphak
      type(real3),intent(in):: p
      type(rpole_elt),intent(in):: mi,mk
      real(t_p),intent(inout):: fip(3),fkp(3),fid(3),fkd(3)
      integer   sca
      real(t_p) r,rr1,rr2,rr3,rr5,rr7
      real(t_p) dir,qix,qiy,qiz,dkr,qkx,qky,qkz,qir,qkr,dmp3,dmp5,dmp7
      real(t_p) rr3i,rr5i,rr7i,rr3k,rr5k,rr7k
      real(t_p) dmpe(4),dmpi(3),dmpk(3),dmpik(3)
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
c     find the field components for Thole polarization damping
c
      call dampdir_inl (pentyp,r,alphai,alphak,dmpi,dmpk)
      IF (IAND(tver,sca).NE.0) THEN
         rr3i = dscale*dmpi(1)*rr3
         rr5i = dscale*dmpi(2)*rr5
         rr7i = dscale*dmpi(3)*rr7
         rr3k = dscale*dmpk(1)*rr3
         rr5k = dscale*dmpk(2)*rr5
         rr7k = dscale*dmpk(3)*rr7
         rr3  = dscale*rr3
      ELSE
         rr3i = (dmpi(1)- 1.0)*rr3
         rr5i = (dmpi(2)- 1.0)*rr5
         rr7i = (dmpi(3)- 1.0)*rr7
         rr3k = (dmpk(1)- 1.0)*rr3
         rr5k = (dmpk(2)- 1.0)*rr5
         rr7k = (dmpk(3)- 1.0)*rr7
         rr3  = 0.0
      END IF
      fid(1) = fid(1) -p%x*(rr3*corek + rr3k*valk
     &            - rr5k*dkr + rr7k*qkr)
     &            - rr3k*mk%dx + 2.0*rr5k*qkx
      fid(2) = fid(2) -p%y*(rr3*corek + rr3k*valk
     &            - rr5k*dkr + rr7k*qkr)
     &            - rr3k*mk%dy + 2.0*rr5k*qky
      fid(3) = fid(3) -p%z*(rr3*corek + rr3k*valk
     &            - rr5k*dkr + rr7k*qkr)
     &            - rr3k*mk%dz + 2.0*rr5k*qkz
      fkd(1) = fkd(3) +p%x*(rr3*corei + rr3i*vali
     &            + rr5i*dir + rr7i*qir)
     &            - rr3i*mi%dx - 2.0*rr5i*qix
      fkd(2) = fkd(2) +p%y*(rr3*corei + rr3i*vali
     &            + rr5i*dir + rr7i*qir)
     &            - rr3i*mi%dy - 2.0*rr5i*qiy
      fkd(3) = fkd(3) +p%z*(rr3*corei + rr3i*vali
     &            + rr5i*dir + rr7i*qir)
     &            - rr3i*mi%dz - 2.0*rr5i*qiz
      IF (IAND(tver,sca).NE.0) THEN
         rr3i = pscale*dmpi(1)*rr3
         rr5i = pscale*dmpi(2)*rr5
         rr7i = pscale*dmpi(3)*rr7
         rr3k = pscale*dmpk(1)*rr3
         rr5k = pscale*dmpk(2)*rr5
         rr7k = pscale*dmpk(3)*rr7
         rr3  = pscale*rr3
      END IF
      fip(1) = fip(1) -p%x*(rr3*corek + rr3k*valk
     &            - rr5k*dkr + rr7k*qkr)
     &            - rr3k*mk%dx + 2.0*rr5k*qkx
      fip(2) = fip(2) -p%y*(rr3*corek + rr3k*valk
     &            - rr5k*dkr + rr7k*qkr)
     &            - rr3k*mk%dy + 2.0*rr5k*qky
      fip(3) = fip(3) -p%z*(rr3*corek + rr3k*valk
     &            - rr5k*dkr + rr7k*qkr)
     &            - rr3k*mk%dz + 2.0*rr5k*qkz
      fkp(1) = fkp(1) +p%x*(rr3*corei + rr3i*vali
     &            + rr5i*dir + rr7i*qir)
     &            - rr3i*mi%dx - 2.0*rr5i*qix
      fkp(2) = fkp(2) +p%y*(rr3*corei + rr3i*vali
     &            + rr5i*dir + rr7i*qir)
     &            - rr3i*mi%dy - 2.0*rr5i*qiy
      fkp(3) = fkp(3) +p%z*(rr3*corei + rr3i*vali
     &            + rr5i*dir + rr7i*qir)
     &            - rr3i*mi%dz - 2.0*rr5i*qiz
      end subroutine
#endif
