#include "image.f.inc"
#include "convert.f.inc"
#include "atomicOp.inc.f"
#include "groups.inc.f"

      M_subroutine 
     &         ker_bond(i,ia,ib,loc,ideal,force,alp,fgrp,xab,yab,zab
     &            ,bndtyp_i
     &            ,use_polymer,use_group
     &            ,cbnd,qbnd,bndunit,eb,e
#ifdef TINKER_CUF
     &            ,debx,deby,debz
#else
     &            ,ded
#endif
     &            ,g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz
     &            ,ver,fea)
!$acc routine
      use bndpot    ,only: BND_HARMONIC,BND_MORSE,BND_MORSE4
      use tinheader ,only: ti_p,zeror
      use tinTypes  ,only: real3
      implicit none
      integer  ,intent(in):: i,ver,fea,bndtyp_i
      integer  ,intent(in):: loc(*)
      real(t_p),intent(in):: cbnd,qbnd,bndunit,ideal,force,alp,fgrp
      real(t_p),intent(inout):: xab,yab,zab
      logical  ,intent(in):: use_group,use_polymer
      integer  ,intent(inout):: ia,ib
#ifdef TINKER_CUF
      ener_rtyp,intent(inout):: eb
      mdyn_rtyp,intent(inout):: debx(*),deby(*),debz(*)
      real(t_p),intent(inout):: g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz
      type(real3) ded
#else
      real(r_p),intent(inout):: eb
      real(r_p),intent(inout):: g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz
      type(real3),intent(out):: ded
#endif
      real(t_p),intent(out):: e

      integer   grd,ene,act,vir,gamd,plm,grp
      real(t_p) de,expterm,bde,dt,dt2
     &         ,deddt
     &         ,rab, ba2
      parameter(
     &         grd=__use_grd__,ene=__use_ene__,
     &         act=__use_act__,vir=__use_vir__,
     &        gamd=__use_gamd__,plm=__use_polymer__,
     &         grp=__use_groups__
     &         )

      if (use_polymer.and.IAND(fea,plm).NE.0)
     &   call image_inl (xab,yab,zab)
      rab = sqrt(xab*xab + yab*yab + zab*zab)
      dt  = rab - ideal

      if (dt.eq.0.0) then
         IF (IAND(ver,ene).NE.0)  e    =0
         IF (IAND(ver,grd+vir).NE.0) THEN
            ded%x=0; ded%y=0; ded%z=0;
         END IF
         return
      end if
c
c     harmonic potential uses Taylor expansion of Morse potential
c     through the fourth power of the bond length deviation
c
      if (bndtyp_i .eq. BND_HARMONIC) then
         dt2 = dt * dt
         IF (IAND(ver,ene).NE.0)
     &      e = bndunit *force *dt2 *(1.0 + cbnd*dt + qbnd*dt2)
         IF (IAND(ver,grd).NE.0)
     &   deddt = 2.0 *bndunit *force *dt
     &         *(1.0 + 1.5*cbnd*dt + 2.0*qbnd*dt2)
c
c     Morse potential uses energy = BDE * (1 - e**(-alpha*dt))**2)
c     with alpha = sqrt(ForceConst/BDE)
c     and BDE = Bond Dissociation Energy = ForceConst/alpha**2
c
      else if (bndtyp_i .eq. BND_MORSE) then
         expterm = exp(-alp*dt)
         bde   = bndunit * force / (alp*alp)
         IF (IAND(ver,ene).NE.0)
     &      e  = bde * (1.0_ti_p-expterm)**2
         IF (IAND(ver,grd+vir).NE.0)
     &   deddt = 2.0_ti_p * bde * alp  
     &       * (1.0_ti_p-expterm) * expterm
      else if (bndtyp_i .eq. BND_MORSE4) then
         dt2 = dt * dt
         ba2 = 7._ti_p/12._ti_p* alp*alp
         IF (IAND(ver,ene).NE.0)
     &      e  = bndunit * force * dt2 
     &        * (1.0_ti_p - alp*dt + ba2*dt2)
         IF (IAND(ver,grd+vir).NE.0)
     &   deddt = bndunit * force * dt 
     &           * (2.0_ti_p-3._ti_p*alp*dt
     &                + 4._ti_p*ba2*dt2)
      end if
c
c     compute chain rule terms needed for derivatives
c
      IF (IAND(ver,grd+vir).NE.0) THEN
         ia    = loc(ia)
         ib    = loc(ib)
         de    = merge( deddt/rab, zeror, rab.ne.zeror)
         ded%x = de* xab
         ded%y = de* yab
         ded%z = de* zab
      END IF

      if (use_group.AND.IAND(fea,grp).NE.0) then
         e     =  e    * fgrp
         deddt = deddt * fgrp
      end if
c
c     increment the total bond energy and first derivatives
c
#ifdef TINKER_CUF
      IF (IAND(ver,ene).NE.0)  eb =  eb + tp2enr(e)
      IF (IAND(ver,grd).NE.0) THEN
         call atomic_add( debx(ia),ded%x )
         call atomic_add( deby(ia),ded%y )
         call atomic_add( debz(ia),ded%z )
         call atomic_sub( debx(ib),ded%x )
         call atomic_sub( deby(ib),ded%y )
         call atomic_sub( debz(ib),ded%z )
      END IF
#else
      IF (IAND(ver,ene).NE.0)  eb =  eb + e
#endif
c
c     increment the internal virial tensor components
c
      IF (IAND(ver,vir).NE.0) THEN
         g_vxx = g_vxx + xab*ded%x
         g_vxy = g_vxy + yab*ded%x
         g_vxz = g_vxz + zab*ded%x
         g_vyy = g_vyy + yab*ded%y
         g_vyz = g_vyz + zab*ded%y
         g_vzz = g_vzz + zab*ded%z
      END IF
      end subroutine
