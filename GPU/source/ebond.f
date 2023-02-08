c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine ebond  --  bond stretch potential energy  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     "ebond" calculates the bond stretching energy
c
c
#include "tinker_macro.h"
      module ebond_inl
#include "atomicOp.h.f"
        contains
#include "ker_bond.inc.f"
      end module

      subroutine ebond
      use atmlst
      use atoms    ,only: type
      use atomsMirror
      use bndpot
      use bond
      use bound
      use deriv
      use domdec
      use ebond_inl
      use energi
      use group
      use nvshmem
      use timestat ,only: timer_enter,timer_exit,timer_ebond
     &             ,quiet_timers
      use tinheader
      use tinTypes ,only: real3
      use usage
      use virial
      implicit none
      integer     i,ia,ib,ibond,fea,ver
      real(t_p)   ideal,force,e,fgrp,xab,yab,zab
      type(real3) ded
      parameter(
     &       ver=__use_ene__,
     &       fea=__use_gamd__+__use_polymer__+__use_groups__
     &         )
c
c     zero out the bond stretching energy
c
      call timer_enter( timer_ebond )
      eb = 0.0_re_p
c
c     calculate the bond stretching energy term
c
!$acc parallel loop present(x,y,z,use,loc,bndglob,grplist,wgrp
!$acc&    ,eb,g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz)
!$acc&     async
#ifdef USE_NVSHMEM_CUDA
!$acc&         deviceptr(d_ibnd,d_bl,d_bk)
#else
!$acc&         present(bl,bk,ibnd)
#endif
!$acc&         reduction(+:eb)
!$acc&         private(fgrp)
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

         end if
      end do
      call timer_exit( timer_ebond,quiet_timers )
      end
