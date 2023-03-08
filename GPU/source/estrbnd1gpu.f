c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine estrbnd1   --  stretch-bend energy and derivs  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "estrbnd1" calculates the stretch-bend potential energy and
c     first derivatives with respect to Cartesian coordinates
c
c
#include "tinker_macro.h"
      module estrbnd1gpu_inl
#include "atomicOp.h.f"
        contains
#include "ker_strbnd.inc.f"
      end module

      subroutine estrbnd1gpu_(isb,typeA,anat,bl,deba,deaMD,s)
      use angle     ,only: iang
      use angpot    ,only: stbnunit
      use atmlst    ,only: strbndglob
      use atoms     ,only: x,y,z
      !use bond
      use bound
      !use deriv
      use domdec
      use energi
      use estrbnd1gpu_inl
      use group
      use inform
      use mamd
      !use math
      use nvshmem
      use potent    ,only: use_amd_wat1
      use strbnd    ,only: nstrbndloc,sbk
      use tinheader ,only: ti_p,re_p
      use usage
      use virial
      implicit none
      integer  ,intent(in):: isb(:,:),typeA(:),s
      real(t_p),intent(in):: anat(:),bl(:)
      real(r_p),intent(out):: deba(:,:),deaMD(3,s)

      integer   grp,tver,tfea
      integer   i,ia,ib,ic,istrbnd,istrbndloc,ipe,ind
      real(t_p) e,fgrp,force1,force2
      logical proceed
      parameter(
     &      grp=__use_groups__,
     &     tver=__use_grd__+__use_ene__+__use_vir__,
     &     tfea=__use_groups__+__use_polymer__+__use_gamd__
     &         )

      if(deb_Path) write(*,*) 'estrbnd1gpu'
c
c     calculate the stretch-bend energy and first derivatives
c
!$acc parallel loop present(eba,g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz,
!$acc&     eW1aMD,deba,deaMD,isb,typeA,anat,bl,iang,sbk,strbndglob)
!$acc&     default(present) async
!$acc&     reduction(+:eba,g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz,eW1aMD)
      do istrbndloc = 1, nstrbndloc
         istrbnd = strbndglob(istrbndloc)
         i      = isb(1,istrbnd)
#ifdef USE_NVSHMEM_CUDA
         ipe     =     (i-1)/nangle_pe
         ind     = mod((i-1),nangle_pe) +1
         ia      = d_iang(ipe)%pel(1,ind)
         ib      = d_iang(ipe)%pel(2,ind)
         ic      = d_iang(ipe)%pel(3,ind)
#else
         ia     = iang(1,i)
         ib     = iang(2,i)
         ic     = iang(3,i)
#endif
         force1 = sbk(1,istrbnd)
         force2 = sbk(2,istrbnd)
         if (use_group.and.IAND(tfea,grp).NE.0)
     &      call groups3_inl(fgrp,ia,ib,ic,ngrp,grplist,wgrp)
c
c     decide whether to compute the current interaction
c
         proceed = (use(ia).or. use(ib).or. use(ic))
c
c     get the coordinates of the atoms in the angle
c
         if (proceed) then
            call ker_strbnd(i,ia,ib,ic,istrbnd,loc,isb,typeA,aMDwattype
     &              ,use_polymer,use_group,use_amd_wat1
     &              ,stbnunit,fgrp,force1,force2,anat,bl,x,y,z
     &              ,eba,eW1aMD,e,deba,deaMD
     &              ,g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz
     &              ,tver,tfea)
         end if
      end do
      end

      subroutine estrbnd1gpu
      use atoms
      use angle
      use bond
      use deriv
      use strbnd
      implicit none
      integer s
      interface
      subroutine estrbnd1gpu_(isb,typeA,anat,bl,deba,deW1aMD,s)
      integer  ,intent(in):: isb(:,:),typeA(:),s
      real(t_p),intent(in):: anat(:),bl(:)
      real(r_p),intent(out):: deba(:,:),deW1aMD(3,s)
      end subroutine
      end interface
      s = size(deW1aMD,2)
      call estrbnd1gpu_(isb,type,anat,bl,deba,deW1aMD,max(s,1))
      end subroutine
