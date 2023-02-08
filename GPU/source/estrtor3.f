c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine estrtor3  --  stretch-torsion energy & analysis  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "estrtor3" calculates the stretch-torsion potential energy;
c     also partitions the energy terms among the atoms
c
c
#include "tinker_macro.h"
      module estrtor3_inl
#include "atomicOp.h.f"
        contains
#include "ker_strtor.inc.f"

      subroutine estrtor3_(itors,x,y,z,tors1,tors2,tors3,ist,kst,bl
     &                    ,debt)
      use action   ,only: nebt
      use atmlst   ,only: strtorglob
      use atomsMirror ,only: n
      use bound    ,only: use_polymer
      use deriv    ,only: deamdD
      use domdec   ,only: rank,nloc,nbloc,loc
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
      real(r_p) x(:),y(:),z(:),debt(3,*)

      integer istrtor,iistrtor,i,ia,ib,ic,id,k,ver,fea,grp
      logical proceed,header
      real(t_p) fgrp
      parameter(
     &          grp=__use_groups__
     &         ,ver=__use_ene__+__use_act__
     &         ,fea=__use_groups__+__use_polymer__
     &         )
      real(t_p) e

      nebt   = 0
      header = rank.eq.0
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
!$acc&     reduction(+:ebt,nebt)
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
            nebt = nebt + 1
#if 0
            !call atomic_add( aebt(ialoc),0.5_ti_p*e1 )
            !call atomic_add( aebt(ibloc),0.5_ti_p*(e1+e2) )
            !call atomic_add( aebt(icloc),0.5_ti_p*(e3+e3) )
            !call atomic_add( aebt(idloc),0.5_ti_p*e1 )
c
c     print a message if the energy of this interaction is large
c
            huge = (abs(e) .gt. 1.0_ti_p)
            if (debug .or. (verbose.and.huge)) then
               if (header) then
                  header = .false.
                  write (iout,10)
   10             format (/,' Individual Stretch-Torsion',
     &                       ' Interactions :',
     &                    //,' Type',25x,'Atom Names',21x,'Angle',
     &                       6x,'Energy',/)
               end if
               write (iout,20)  ia,name(ia),ib,name(ib),ic,
     &                          name(ic),id,name(id),angle,e
   20          format (' StrTors',3x,4(i7,'-',a3),f11.4,f12.4)
            end if
#endif
         end if
      end do
      end
      end module

      subroutine estrtor3
      use atomsMirror
      use bond
      use deriv
      use estrtor3_inl
      use strtor
      use tors
      call estrtor3_(itors,x,y,z,tors1,tors2,tors3,ist,kst,bl,debt)
      end subroutine
