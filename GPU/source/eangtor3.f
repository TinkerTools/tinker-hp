c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine eangtor3  --  angle-torsion energy & analysis  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "eangtor3" calculates the angle-torsion potential energy;
c     also partitions the energy terms among the atoms
c
c
#include "tinker_macro.h"
      module eangtor3_inl
#include "atomicOp.h.f"
      contains
#include "ker_angtor.inc.f"
      end module

      subroutine eangtor3_(iat,anat,kant,tors1,tors2,tors3,deat)
      use action
      use analyz
      !use angle
      use angtor ,only: nangtor,nangtorloc
      use atmtyp
      use atmlst
      use atoms
      use bound
      use domdec   ,only: rank,loc
      use eangtor3_inl
      use energi
      use group
      use inform
      use iounit
      use math
      use tinheader,only: ti_p
      use torpot
      use tors     ,only: itors
      use usage
      use virial
      implicit none
      integer  ,intent(in):: iat(:,:)
      real(t_p),intent(in)::anat(:),kant(:,:)
     &         ,tors1(:,:),tors2(:,:),tors3(:,:)
      real(r_p):: deat(1)

      integer grp
      integer i,iiangtor,iangtor,tver,tfea
      integer ia,ib,ic,id
      logical proceed,header,huge
      real(t_p) fgrp,e,angle1
      parameter(
     &       grp=__use_groups__,
     &      tver=__use_ene__+__use_act__,
     &      tfea=__use_polymer__+__use_groups__
     &         )
c
      if (deb_Path) print*, "eangtor3"
#ifdef USE_NVSHMEM_CUDA
      ! TODO Remove this check
      ! Implement NVSHMEM Access to anat & itors, tors[123]
      print*, '  FATAL ERROR  '
      print*, 'NVSHMEM feature not implemented inside eangtor3'
      __TINKER_FATAL__
#endif
c
c     zero out the energy due to extra potential terms
c
      neat = 0
      eat  = 0.0_re_p
      aet  = 0.0_ti_p
c
c     print header information if debug output was requested
c
      header = (rank.eq.0)
      if (debug .and. nangtor.ne.0) then
         header = .false.
         write (iout,10)
   10    format (/,' Individual Angle-Torsion Interactions :',
     &           //,' Type',25x,'Atom Names',21x,'Angle',
     &              6x,'Energy',/)
      end if
c
c     calculate the angle-torsion interaction energy term
c
!$acc parallel loop async reduction(+:eat,neat)
!$acc&         present(eat,aeat,angtorglob,iat,anat,wgrp,grplist
!$acc&     ,itors,tors1,tors2,tors3,x,y,z,use,kant,loc)
      do iangtor = 1, nangtorloc
         iiangtor = angtorglob(iangtor)
         i        = iat(1,iiangtor)
         ia       = itors(1,i)
         ib       = itors(2,i)
         ic       = itors(3,i)
         id       = itors(4,i)
c
c     decide whether to compute the current interaction
c
         if (use_group.and.IAND(tver,grp).NE.0)
     &      call groups4_inl (fgrp,ia,ib,ic,id,ngrp,grplist,wgrp)
         proceed = (use(ia) .or. use(ib) .or.
     &                              use(ic) .or. use(id))
c
c     compute the value of the torsional angle
c
         if (proceed) then
            call ker_angtor(iiangtor,i,ia,ib,ic,id,loc,iat,radian
     &              ,atorunit,fgrp,x,y,z,anat,kant,tors1,tors2,tors3
     &              ,use_polymer,use_group,use_virial
     &              ,eat,e,deat,g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz
     &              ,tver,tfea)

            if (e.ne.0) neat = neat + 1
#ifndef _OPENACC
c           ia = loc(ia)
c           ib = loc(ib)
c           ic = loc(ic)
c           id = loc(id)
c           aeat(ia) = aeat(ia) + e1/3.0_ti_p
c           aeat(ib) = aeat(ib) + e/3.0_ti_p
c           aeat(ic) = aeat(ic) + e/3.0_ti_p
c           aeat(id) = aeat(id) + e2/3.0_ti_p
c
c     print a message if the energy of this interaction is large
c
            huge = (e .gt. 3.0d0)
            if (debug .or. (verbose.and.huge)) then
               if (header) then
                  header = .false.
                  write (iout,20)
   20             format (/,' Individual Angle-Torsion',
     &                       ' Interactions :',
     &                    //,' Type',25x,'Atom Names',21x,'Angle',
     &                       6x,'Energy',/)
               end if
               write (iout,30)  ia,name(ia),ib,name(ib),ic,
     &                          name(ic),id,name(id),angle1,e
   30          format (' AngTors',3x,4(i7,'-',a3),f11.4,f12.4)
            end if
#endif
         end if
      end do
      end

      subroutine eangtor3
      use angle
      use angtor
      use deriv
      use tors
      implicit none
      interface
      subroutine eangtor3_(iat,anat,kant,tors1,tors2,tors3,deat)
      integer  ,intent(in):: iat(:,:)
      real(t_p),intent(in)::anat(:),kant(:,:)
     &         ,tors1(:,:),tors2(:,:),tors3(:,:)
      real(r_p):: deat(*)
      end subroutine; end interface
      call eangtor3_(iat,anat,kant,tors1,tors2,tors3,deat)
      end subroutine
