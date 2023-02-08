c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine estrtor1  --  stretch-torsion energy & derivs  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "estrtor1" calculates the stretch-torsion energy and first
c     derivatives with respect to Cartesian coordinates
c
c
#include "tinker_macro.h"
      module estrtor1gpu_inl
#include "atomicOp.h.f"
        contains

#include "ker_strtor.inc.f"

      subroutine estrtor1gpu_(itors,x,y,z,tors1,tors2,tors3,ist,kst,bl
     &                       ,debt)
      use atmlst   ,only: strtorglob
      use atomsMirror ,only: n
      use bound    ,only: use_polymer
      use deriv    ,only: deamdD
      use domdec   ,only: nloc,nbloc,loc
      use energi   ,only: ebt,eDaMD
      use group    ,only: use_group,ngrp,grplist,wgrp
      use nvshmem
      use potent   ,only: use_amd_dih
      use strtor   ,only: nstrtorloc
      use tinheader
      use torpot   ,only: storunit
      use usage
      use virial   ,only: g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz
      implicit none
      integer   ist(:,:),itors(:,:)
      real(t_p) tors1(:,:),tors2(:,:),tors3(:,:),kst(:,:),bl(:)
      real(r_p) x(:),y(:),z(:),debt(:,:)

      integer istrtor,iistrtor,i,ia,ib,ic,id,k,ver,fea,grp
      logical proceed
      real(t_p) fgrp
      parameter(
     &          grp=__use_groups__
     &         ,ver=__use_grd__+__use_ene__+__use_vir__
     &         ,fea=__use_groups__+__use_polymer__
     &         )
      real(t_p) e
c
c     calculate the stretch-torsion interaction energy term
c
!$acc parallel loop
#ifdef USE_NVSHMEM_CUDA
!$acc&     default(present) deviceptr(d_bl)
#else
!$acc&     present(debt,bl,ist,kst,strtorglob,itors,loc,x,y,z,use
!$acc&            ,tors1,tors2,tors3,grplist,wgrp)
#endif
!$acc&     present(ebt,g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz) async
!$acc&     reduction(+:ebt,g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz)
      do istrtor = 1, nstrtorloc
         iistrtor = strtorglob(istrtor)
         i        = ist(1,iistrtor)
#ifdef USE_NVSHMEM_CUDA
         ipe   = (i-1)/ntors_pe
         ind   = mod((i-1),ntors_pe) +1
         ia    = d_itors(ipe)%pel(1,ind)
         ib    = d_itors(ipe)%pel(2,ind)
         ic    = d_itors(ipe)%pel(3,ind)
         id    = d_itors(ipe)%pel(4,ind)
#else    
         ia    = itors(1,i)
         ib    = itors(2,i)
         ic    = itors(3,i)
         id    = itors(4,i)
#endif
c
c     decide whether to compute the current interaction
c
         if (IAND(fea,grp).NE.0.and.use_group)
     &      call groups4_inl (fgrp,ia,ib,ic,id,ngrp,grplist,wgrp)
         proceed = useAll .or.(use(ia).or. use(ib)
     &                    .or. use(ic).or. use(id))
c
c     compute the value of the torsional angle
c
         if (proceed) then
            call ker_strtor(iistrtor,i,ia,ib,ic,id,ist,loc,x,y,z
     &              ,tors1,tors2,tors3,kst,bl,storunit
     &              ,nbloc,n,use_group,use_polymer
     &              ,e,ebt,debt,g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz
     &              ,ver,fea)
         end if
      end do
c
c     increment the aMD/GaMD derivs table
c
      if (use_amd_dih) then
!$acc parallel loop collapse(2) async
!$acc&         present(deamdD,debt)
         do i = 1,nloc; do k = 1,3
            deamdD(k,i) = deamdD(k,i) + debt(k,i)
         end do; end do
!$acc serial async present(eDaMD,ebt)
         eDaMD = eDaMD + ebt
!$acc end serial
      end if
      end
      end module

      subroutine estrtor1gpu
      use atomsMirror
      use bond
      use deriv
      use estrtor1gpu_inl
      use strtor
      use tors
      call estrtor1gpu_(itors,x,y,z,tors1,tors2,tors3,ist,kst,bl,debt)
      end subroutine
