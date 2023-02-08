c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine eopbend3  --  out-of-plane bending & analysis  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "eopbend3" computes the out-of-plane bend potential energy at
c     trigonal centers via a Wilson-Decius-Cross or Allinger angle;
c     also partitions the energy among the atoms
c
c
#include "tinker_macro.h"
      module eopbend3gpu_inl
#include "atomicOp.h.f"
        contains
#include "ker_opbend.inc.f"
      end module

      subroutine eopbend3gpu_(deopb)
      use action
      use analyz
      use angle
      use angpot
      use atmlst
      use atmtyp
      use atomsMirror
      use bound
      use domdec
      use energi
      use eopbend3gpu_inl
      use group
      use inform
      use iounit
      use math
      use opbend
      use tinheader ,only:ti_p,re_p
      use timestat
      use usage
      use virial
      implicit none
      real(r_p) deopb(1)
      integer i,iopbend,iopbendloc
      integer ia,ib,ic,id
      integer tver,tfea
#ifdef USE_NVSHMEM_CUDA
      integer ipe,ind
#endif
      real(t_p) e,force,fgrp
      logical   proceed,header,huge
      parameter(
     &     tver=__use_ene__+__use_act__,
     &     tfea=__use_groups__+__use_polymer__
     &         )

      if(deb_Path) write(*,*) 'eopbend3gpu'
      call timer_enter( timer_eopbend1 )
c
c     zero out the out-of-plane bend energy and partitioning
c
      neopb  = 0
       eopb  = 0
c     aeopb  = 0
      header = (rank.eq.0)

!$acc parallel loop 
#ifdef USE_NVSHMEM_CUDA
!$acc&         present(x,y,z,loc,use,opbendglob,opbk
#else
!$acc&         present(x,y,z,loc,use,opbendglob,iang,opbk
#endif
!$acc&    ,iopb,eopb) reduction(+:eopb,neopb) async
      do iopbendloc = 1, nopbendloc
         iopbend = opbendglob(iopbendloc)
         i = iopb(iopbend)
#ifdef USE_NVSHMEM_CUDA
         ipe=     (i-1)/nangle_pe
         ind= mod((i-1),nangle_pe) +1
         ia = d_iang(ipe)%pel(1,ind)
         ib = d_iang(ipe)%pel(2,ind)
         ic = d_iang(ipe)%pel(3,ind)
         id = d_iang(ipe)%pel(4,ind)
#else
         ia = iang(1,i)
         ib = iang(2,i)
         ic = iang(3,i)
         id = iang(4,i)
#endif
         force = opbk(iopbend)

         if (use_group.and.IAND(tfea,__use_groups__).NE.0)
     &      call groups4_inl(fgrp,ia,ib,ic,id,ngrp,grplist,wgrp)
c
c     decide whether to compute the current interaction
c
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

            if (e.ne.0.0) neopb = neopb + 1
#if 0
            ib = loc(ib)
            aeopb(ib) = aeopb(ib) + e
c
c     print a message if the energy of this interaction is large
c
            huge = (e .gt. 2.0_ti_p)
            if (debug .or. (verbose.and.huge)) then
               if (header) then
                  header = .false.
                  write (iout,10)
   10             format (/,' Individual Out-of-Plane Bend',
     &                       ' Interactions :',
     &                    //,' Type',25x,'Atom Names',21x,'Angle',
     &                       6x,'Energy',/)
               end if
               write (iout,20)  id,name(id),ib,name(ib),ia,
     &                          name(ia),ic,name(ic),angle1,e
   20          format (' O-P-Bend',2x,4(i7,'-',a3),f11.4,f12.4)
            end if
#endif
         end if
      end do
      call timer_exit( timer_eopbend1 )
      end

      subroutine eopbend3gpu
      use deriv
      implicit none
      call eopbend3gpu_(deopb)
      end subroutine
