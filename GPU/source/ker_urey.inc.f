#include "image.f.inc"
#include "groups.inc.f"
#include "atomicOp.inc.f"

      M_subroutine
     &        ker_urey(i,ia,ic,nbloc,loc
     &            ,ideal,force,ureyunit,cury,qury,fgrp
     &            ,use_group,use_polymer
     &            ,x,y,z
     &            ,eub
#ifdef TINKER_CUF
     &            ,deubx,deuby,deubz
#else
     &            ,e,ded
     &            ,g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz
#endif
     &            ,tver,tfea)
      use tinTypes,only: real3
      implicit none
      integer  ,intent(in) :: i,nbloc,tver,tfea,loc(*)
      logical  ,intent(in) :: use_group,use_polymer
      real(t_p),intent(in) :: ideal,force,ureyunit,cury,qury,fgrp
     &         , x(*),y(*),z(*)
      integer  ,intent(inout) :: ia,ic
#ifdef TINKER_CUF
      ener_rtyp,intent(inout):: eub
      mdyn_rtyp,intent(inout):: deux(*),deuy(*),deuz(*)
      real(t_p),intent(inout):: g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz
      type(real3) ded
      real(t_p)   e
#else
      type(real3),intent(out):: ded
      real(r_p),intent(inout):: eub
     &         , g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz
      real(t_p),  intent(out):: e
#endif
      integer :: grd,ene,act,vrl,gamd,grp,plm
      real(t_p) de
      real(t_p) dt,dt2,deddt
      real(t_p) xac,yac,zac,rac
      parameter(
     &         grd=__use_grd__,
     &         ene=__use_ene__,
     &         vrl=__use_vir__,
     &         act=__use_act__,
     &        gamd=__use_gamd__,
     &         grp=__use_groups__,
     &         plm=__use_polymer__
     &         )
!$acc routine
      xac = x(ia) - x(ic)
      yac = y(ia) - y(ic)
      zac = z(ia) - z(ic)
      if (IAND(tfea,plm).NE.0.and.use_polymer)
     &   call image_inl (xac,yac,zac)

      rac = sqrt(xac*xac + yac*yac + zac*zac)
      dt  = rac - ideal
      dt2 = dt * dt
      e   = ureyunit *force *dt2 *(1.0+cury*dt+qury*dt2)
      deddt = 2.0 *ureyunit *force *dt
     &            *(1.0 + 1.5*cury*dt + 2.0*qury*dt2)
c
c     scale the interaction based on its group membership
c
      if (IAND(tfea,grp).NE.0.and.use_group) then
          e    =  e    * fgrp
         deddt = deddt * fgrp
      end if
c
c     increment the total Urey-Bradley energy
c
#ifdef TINKER_CUF
      IF (IAND(tver,ene).NE.0)  eub =  eub + tp2enr(e)
#else
      IF (IAND(tver,ene).NE.0)  eub =  eub + e
#endif

      IF (IAND(tver,grd).NE.0) THEN
      ia = loc(ia)
      ic = loc(ic)
c
c     compute chain rule terms needed for derivatives
c
      de    = deddt / rac
      ded%x = de * xac
      ded%y = de * yac
      ded%z = de * zac

#ifdef TINKER_CUF
      call atomic_add_f1( deubx(ia),ded%x )
      call atomic_add_f1( deuby(ia),ded%y )
      call atomic_add_f1( deubz(ia),ded%z )
      call atomic_sub_f1( deubx(ic),ded%x )
      call atomic_sub_f1( deuby(ic),ded%y )
      call atomic_sub_f1( deubz(ic),ded%z )
#endif
      END IF
c
c     increment the internal virial tensor components
c
      IF (IAND(tver,vrl).NE.0) THEN
         g_vxx = g_vxx + xac*ded%x
         g_vxy = g_vxy + yac*ded%x
         g_vxz = g_vxz + zac*ded%x
         g_vyy = g_vyy + yac*ded%y
         g_vyz = g_vyz + zac*ded%y
         g_vzz = g_vzz + zac*ded%z
      END IF
      end subroutine
