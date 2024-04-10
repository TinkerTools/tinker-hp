#ifndef PAIR_CHGTRN_INC
#define PAIR_CHGTRN_INC
#include "tinker_cudart.h"
#include "convert.f.inc"
#include "switch_respa.f.inc"

      M_subroutine duo_chgtrn
     &            (r2,xr,yr,zr,chgi,chgk,alphai,alphak,f,scal,fgrp
     &            ,rinv,off,ctrnscut,cut2,sheal,ctrntyp,u_grp,e,frc
     &            ,ver,fea)
      use ctrpot   ,only: CHGT_SEPARATE
      use tinTypes ,only: real3
      implicit none
      integer    ,intent(in):: ctrntyp,ver,fea
      logical    ,intent(in):: u_grp
      real(t_p)  ,intent(in):: r2,xr,yr,zr,chgi,chgk,alphai,alphak,fgrp
     &           ,f,scal,off,ctrnscut,cut2,sheal,rinv
      real(t_p)  ,intent(inout):: e
      type(real3),intent(out):: frc

      integer   grd,ene,sca,grp,shr,lgr
      real(t_p) rik,rr1,expi,expk,chgik,alphaik,de
      parameter(grd=__use_grd__,ene=__use_ene__,sca=__use_sca__
     &         ,grp=__use_groups__,shr=__use_shortRange__
     &         ,lgr=__use_longRange__
     &         )

      rik = f_sqrt(r2)
      !TODO when Loading check alpha[ik] value
      !TODO if (alpha[ik] .eq. 0.0) alpha[ik] = 100.0
      if (ctrntyp .eq. CHGT_SEPARATE) then
         expi = f_exp(-alphai*rik)
         expk = f_exp(-alphak*rik)
         IF (IAND(ver,ene).NE.0)  e=-chgi*expk - chgk*expi
         IF (IAND(ver,grd).NE.0) de= chgi*expk*alphak + chgk*expi*alphai
      else
         !associate(chgik=>expi,alphaik=>expk)
         chgik   = f_sqrt(f_abs(chgi*chgk))
         alphaik = 0.5* (alphai+alphak)
         IF (IAND(ver,ene).NE.0)  e= -chgik* f_exp(-alphaik*rik)
         IF (IAND(ver,grd).NE.0) de= -e* alphaik
         !end associate
      end if

      IF (IAND(ver,sca).NE.0) THEN
         IF (IAND(ver,ene).NE.0) e  = f*  e* scal
         IF (IAND(ver,grd).NE.0) de = f* de* scal
      ELSE
         IF (IAND(ver,ene).NE.0) e  = f* e
         IF (IAND(ver,grd).NE.0) de = f* de
      END IF
c
c     use energy switching if near the cutoff distance
c
      IF (IAND(fea,shr).EQ.0) THEN
         if (r2 .gt. cut2) then
            block
            real(t_p) r,r2,r3,taper,dtaper
            r      = (rik - off)*rinv
            r2     = r  * r
            r3     = r2 * r
            taper  = r3 * (6*r2 - 15*r + 10)
            dtaper = 30* (r*(1.0-r))*(r*(1.0-r)) *rinv;
            IF (IAND(ver,grd).NE.0) de = e*dtaper + de*taper
            IF (IAND(ver,ene).NE.0)  e = e * taper
            end block
        end if
      END IF
c
c     scale the interaction based on its group membership
c
      if (IAND(fea,grp).NE.0.and.u_grp) then
         IF (IAND(ver,ene).NE.0)  e =  e* fgrp
         IF (IAND(ver,grd).NE.0) de = de* fgrp
      end if
      IF (IAND(fea,shr+lgr).NE.0) THEN
         block
         real(t_p) s,ds
         call switch_respa_inl(rik,ctrnscut,sheal,s,ds)
         IF (IAND(fea,shr).NE.0) THEN
            IF (IAND(ver,grd).NE.0) de = de*s + e*ds 
            IF (IAND(ver,ene).NE.0)  e =  e*s
         ELSE
            IF (IAND(ver,grd).NE.0) de = de*(1.0-s) - e*ds 
            IF (IAND(ver,ene).NE.0)  e =  e*(1.0-s)
         END IF
         end block
      END IF
c
c     compute the force components for this interaction
c
      IF (IAND(ver,grd).NE.0) THEN
         rr1 = rik**(-1)
         frc%x = de * xr * rr1
         frc%y = de * yr * rr1
         frc%z = de * zr * rr1
      END IF
      end subroutine
#endif
