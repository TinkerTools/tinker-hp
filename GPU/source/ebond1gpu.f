c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine ebond1  --  bond stretch energy & derivatives  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "ebond1" calculates the bond stretching energy and
c     first derivatives with respect to Cartesian coordinates
c
c
#include "tinker_macro.h"
      module ebond1gpu_inl
#include "atomicOp.h.f"
        contains
#include "ker_bond.inc.f"
      end module

      subroutine ebond1gpu
      use atmlst
      use atoms       ,only: type,n
      use atomsMirror ,only: x,y,z
      use bndpot
      use bond
      use bound
      use deriv
      use domdec
      use energi
      use ebond1gpu_inl
      use group
      use inform   ,only: deb_Path,minmaxone
      use nvshmem
      use usage
      use virial
      use timestat ,only: timer_enter,timer_exit,timer_ebond
     &             ,quiet_timers
      use tinheader
      use tinTypes ,only: real3
      use mamd
      use potent   ,only: use_amd_wat1
      implicit none
      integer i,ia,ib,ibond,fea,ver
      real(t_p) ideal,force,e,fgrp,xab,yab,zab
      type(real3) ded
      parameter(
     &       ver=__use_grd__+__use_ene__+__use_vir__,
     &       fea=__use_gamd__+__use_polymer__+__use_groups__
     &         )

      call timer_enter( timer_ebond )
      if (deb_Path) write(*,*) 'ebond1gpu'
c
c     calculate the bond stretch energy and first derivatives
c
!$acc parallel loop present(x,y,z,use,loc,bndglob,grplist,wgrp,type
!$acc&    ,aMDwattype,deW1amd,deb,eb,eW1aMD
!$acc&    ,g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz)
!$acc&     async
#ifdef USE_NVSHMEM_CUDA
!$acc&         deviceptr(d_ibnd,d_bl,d_bk)
#else
!$acc&         present(bl,bk,ibnd)
#endif
!$acc&       reduction(+:eb,eW1aMD,g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz)
!$acc&       private(ded,fgrp)
      do ibond = 1, nbondloc
         i     = bndglob(ibond)
#ifdef USE_NVSHMEM_CUDA
         ipe   = (i-1)/nbond_pe
         ind   = mod((i-1),nbond_pe) +1
         ia    = d_ibnd(ipe)%pel(1,ind)
         ib    = d_ibnd(ipe)%pel(2,ind)
         ideal = d_bl  (ipe)%pel(ind)
         force = d_bk  (ipe)%pel(ind)
#else
         ia    = ibnd(1,i)
         ib    = ibnd(2,i)
         ideal = bl(i)
         force = bk(i)
#endif
         if (use_group.AND.IAND(fea,__use_groups__).NE.0)
     &      call groups2_inl(fgrp,ia,ib,ngrp,grplist,wgrp)
c
c     compute the value of the bond length deviation
c     decide whether to compute the current interaction
c
         if (useAll.or.use(ia).or.use(ib)) then

         xab = x(ia) - x(ib)
         yab = y(ia) - y(ib)
         zab = z(ia) - z(ib)
         call ker_bond(i,ia,ib,loc,ideal,force,fgrp,xab,yab,zab
     &           ,bndtyp_i,use_polymer,use_group
     &           ,cbnd,qbnd,bndunit,eb,e,ded
     &           ,g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz
     &           ,ver,fea)

        !Increment the total bond first derivatives
         call atomic_add( deb(1,ia),ded%x )
         call atomic_add( deb(2,ia),ded%y )
         call atomic_add( deb(3,ia),ded%z )
         call atomic_sub( deb(1,ib),ded%x )
         call atomic_sub( deb(2,ib),ded%y )
         call atomic_sub( deb(3,ib),ded%z )
c
c     aMD storage if waters are considered
c
         IF (use_amd_wat1.and.IAND(fea,__use_gamd__).NE.0) THEN
         if (type(ia)==aMDwattype(1) .or.
     $       type(ib)==aMDwattype(1)) then
            eW1aMD = eW1aMD + e
            call atomic_add_m( deW1amd(1,ia),ded%x )
            call atomic_add_m( deW1amd(2,ia),ded%y )
            call atomic_add_m( deW1amd(3,ia),ded%z )

            call atomic_sub_m( deW1amd(1,ib),ded%x )
            call atomic_sub_m( deW1amd(2,ib),ded%y )
            call atomic_sub_m( deW1amd(3,ib),ded%z )
         end if
         END IF

         end if
      end do
      call timer_exit( timer_ebond )
      end
