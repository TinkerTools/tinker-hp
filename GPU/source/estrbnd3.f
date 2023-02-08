c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine estrbnd3  --  stretch-bend energy & analysis  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "estrbnd3" calculates the stretch-bend potential energy;
c     also partitions the energy among the atoms
c
c
#include "tinker_macro.h"
      module estrbnd3_inl
#include "atomicOp.h.f"
        contains
#include "ker_strbnd.inc.f"
      end module

      subroutine estrbnd3_(isb,typeA,anat,bl,deba,deaMD)
      use action
      use analyz
      use angle     ,only: iang
      use angpot    ,only: stbnunit
      use atmlst    ,only: strbndglob
      use atoms     ,only: x,y,z
      !use bond
      use bound
      use domdec
      use energi
      use estrbnd3_inl
      use group
      use inform
      use iounit
      use mamd
      use nvshmem
      use potent    ,only: use_amd_wat1
      use strbnd    ,only: nstrbndloc,sbk
      use tinheader ,only: ti_p,re_p
      use usage
      use sizes     ,only: tinkerdebug
      use virial
      implicit none
      integer  ,intent(in):: isb(:,:),typeA(:)
      real(t_p),intent(in):: anat(:),bl(:)
      real(r_p),intent(out):: deba(1),deaMD(1)

      integer   grp,tver,tfea
      integer   i,ia,ib,ic,istrbnd,istrbndloc,ipe,ind
      logical   proceed,header,huge
      real(t_p) e,fgrp,force1,force2
      parameter(
     &      grp=__use_groups__,
     &     tver=__use_ene__+__use_act__,
     &     tfea=__use_groups__+__use_polymer__+__use_gamd__
     &         )
c
c     zero out the energy component and partitioning terms
c
      if(deb_Path) write(*,*) 'estrbnd3'
      neba   = 0
      eba    = 0
      aeba   = 0
      header = (rank.eq.0)
c
c     calculate the stretch-bend energy term
c
!$acc parallel loop present(eba,g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz,
!$acc&     eW1aMD,deba,deaMD,isb,typeA,anat,bl,iang,sbk,strbndglob)
!$acc&     default(present) async
!$acc&     reduction(+:eba,neba)
      do istrbndloc = 1, nstrbndloc
         istrbnd = strbndglob(istrbndloc)
         i       = isb(1,istrbnd)
#ifdef USE_NVSHMEM
         ipe     =     (i-1)/nangle_pe
         ind     = mod((i-1),nangle_pe) +1
         ia      = d_iang(ipe)%pel(1,ind)
         ib      = d_iang(ipe)%pel(2,ind)
         ic      = d_iang(ipe)%pel(3,ind)
#else
         ia      = iang(1,i)
         ib      = iang(2,i)
         ic      = iang(3,i)
#endif
         force1  = sbk(1,istrbnd)
         force2  = sbk(2,istrbnd)
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
            if(e.ne.0.0) neba = neba + 1
#if 0
            ib = loc(ib)
            aeba(ib) = aeba(ib) + e
c
c     print a message if the energy of this interaction is large
c
            huge = (abs(e) .gt. 2.0_ti_p)
            if (debug .or. (verbose.and.huge)) then
               if (header) then
                  header = .false.
                  write (iout,10)
   10             format (/,' Individual Stretch-Bend',
     &                       ' Interactions :',
     &                    //,' Type',18x,'Atom Names',18x,'dSB 1',
     &                       5x,'dSB 2',6x,'Energy',/)
               end if
               write (iout,20)  ia,name(ia),ib,name(ib),
     &                          ic,name(ic),dr1*dt,dr2*dt,e
   20          format (' StrBend',3x,3(i7,'-',a3),2x,2f10.4,f12.4)
            end if
#endif
         end if
      end do
      end

      subroutine estrbnd3
      use atoms
      use angle
      use bond
      use deriv
      use strbnd
      implicit none
      interface
      subroutine estrbnd3_(isb,typeA,anat,bl,deba,deW1aMD)
      integer  ,intent(in):: isb(:,:),typeA(:)
      real(t_p),intent(in):: anat(:),bl(:)
      real(r_p),intent(out):: deba(*),deW1aMD(*)
      end subroutine
      end interface
      call estrbnd3_(isb,type,anat,bl,deba,deW1aMD)
      end subroutine
