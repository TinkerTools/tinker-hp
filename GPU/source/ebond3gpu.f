c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine ebond3  --  bond stretch energy & analysis  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "ebond3" calculates the bond stretching energy; also
c     partitions the energy among the atoms
c
c
#include "tinker_macro.h"
      module ebond3gpu_inl
#include "atomicOp.h.f"
        contains
#include "ker_bond.inc.f"
      end module

      subroutine ebond3gpu
      use action
      use analyz
      use atmlst
      use atoms       ,only: type
      use atomsMirror ,only: x,y,z
      use bndpot
      use bond
      use bound
      use deriv
      use domdec
      use energi
      use ebond3gpu_inl
      use group
      use inform   ,only: debug,verbose,deb_Path
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
      integer     i,ia,ib,ibond,fea,ver
      logical     header,huge
      real(t_p)   ideal,force,e,fgrp,xab,yab,zab
      type(real3) ded
      parameter(
     &       ver=__use_ene__+__use_act__,
     &       fea=__use_gamd__+__use_polymer__+__use_groups__
     &         )

      call timer_enter( timer_ebond )
      if (deb_Path) write(*,*) 'ebond3gpu'
      neb    = 0
       eb    = 0
      header = rank.eq.0
c
c     calculate the bond stretch energy and first derivatives
c
!$acc parallel loop present(x,y,z,use,loc,bndglob,grplist,wgrp,type
!$acc&    ,eb,g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz)
!$acc&     async
#ifdef USE_NVSHMEM_CUDA
!$acc&         deviceptr(d_ibnd,d_bl,d_bk)
#else
!$acc&         present(bl,bk,ibnd)
#endif
!$acc&       reduction(+:eb,neb)
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
         neb = neb + 1
#if 0
         ia = loc(ia); ib = loc(ib)
         aeb(ia) = aeb(ia) + 0.5_ti_p*e
         aeb(ib) = aeb(ib) + 0.5_ti_p*e
c
c     print a message if the energy of this interaction is large
c
         huge = (e .gt. 5.0_ti_p)
         if (debug .or. (verbose.and.huge)) then
            if (header) then
               header = .false.
               write (iout,10)
   10          format (/,' Individual Bond Stretching',
     &                    ' Interactions :',
     &                 //,' Type',14x,'Atom Names',22x,'Ideal',
     &                    4x,'Actual',6x,'Energy',/)
            end if
            write (iout,20)  ia,name(ia),ib,name(ib),ideal,rab,e
   20       format (' Bond',6x,2(i7,'-',a3),13x,2f10.4,f12.4)
         end if
#endif
         end if
      end do
      call timer_exit( timer_ebond )
      end
