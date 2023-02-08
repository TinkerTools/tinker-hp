c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine eopbend  --   out-of-plane bending & analysis  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "eopbend3" computes the out-of-plane bend potential energy at
c     trigonal centers via a Wilson-Decius-Cross or Allinger angle
c
c
#include "tinker_macro.h"
      module eopbend_inl
#include "atomicOp.h.f"
        contains
#include "ker_opbend.inc.f"
      end module

      subroutine eopbend_(deopb)
      use action
      use angle
      use angpot
      use atmlst
      use atomsMirror
      use bound
      use domdec
      use energi
      use eopbend_inl
      use group
      use inform    ,only: deb_Path
      use math
      use opbend
      use tinheader ,only: ti_p,re_p
      use virial
      use usage
      implicit none
      real(r_p) deopb(1)
      integer i,iopbend,iopbendloc
      integer ia,ib,ic,id
      integer tver,tfea
#ifdef USE_NVSHMEM_CUDA
      integer ipe,ind
#endif
      real(t_p) e,force,fgrp
      logical   proceed
      parameter(
     &     tver=__use_ene__,
     &     tfea=__use_groups__+__use_polymer__
     &         )
c
c     zero out-of-plane energy
c
      eopb = 0

!$acc parallel loop async
#ifdef USE_NVSHMEM_CUDA
!$acc&      present(loc,use,opbendglob,opbk,iopb
#else
!$acc&      present(loc,use,opbendglob,iang,opbk,iopb
#endif
!$acc&     ,deopb,eopb,g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz)
!$acc&      reduction(+:eopb)
      do iopbendloc = 1, nopbendloc
         iopbend = opbendglob(iopbendloc)
         i       = iopb(iopbend)
#ifdef USE_NVSHMEM_CUDA
         ipe     =     (i-1)/nangle_pe
         ind     = mod((i-1),nangle_pe) +1
         ia      = d_iang(ipe)%pel(1,ind)
         ib      = d_iang(ipe)%pel(2,ind)
         ic      = d_iang(ipe)%pel(3,ind)
         id      = d_iang(ipe)%pel(4,ind)
#else
         ia      = iang(1,i)
         ib      = iang(2,i)
         ic      = iang(3,i)
         id      = iang(4,i)
#endif
         force   = opbk(iopbend)
c
c     decide whether to compute the current interaction
c
         if (use_group.and.IAND(tfea,__use_groups__).NE.0)
     &      call groups4_inl (fgrp,ia,ib,ic,id,ngrp,grplist,wgrp)
         proceed = (use(ia).or. use(ib).or. use(ic).or. use(id))
c
c     get the coordinates of the atoms at trigonal center
c
         if (proceed) then
            call ker_opbend(ia,ib,ic,id,opbtypInt,loc
     &              ,use_group,use_polymer
     &              ,opbunit,fgrp,force,copb,qopb,popb,sopb,x,y,z
     &              ,eopb,e,deopb,g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz
     &              ,tver,tfea)
         end if
      end do
      end

      subroutine eopbend
      use deriv
      implicit none
      call eopbend_(deopb)
      end subroutine
