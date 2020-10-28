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
#include "tinker_precision.h"
      module eopbend3gpu_inl
        contains
#include "image.f.inc"
      end module

      subroutine eopbend3gpu
      use action
      use analyz
      use angle
      use angpot
      use atmlst
      use atmtyp
      use atoms
      use bound
      use domdec
      use energi
      use eopbend3gpu_inl
      use group
      use inform
      use iounit
      use math
      use opbend
      use usage
      use tinheader ,only:ti_p,re_p
      use timestat
      implicit none
      integer i,iopbend,iopbendloc
      integer ia,ib,ic,id,ibloc
#ifdef USE_NVSHMEM_CUDA
      integer ipe,ind
#endif
      real(t_p) e,angle1,force
      real(t_p) cosine
      real(t_p) dt,dt2,dt3,dt4
      real(t_p) xia,yia,zia
      real(t_p) xib,yib,zib
      real(t_p) xic,yic,zic
      real(t_p) xid,yid,zid
      real(t_p) xab,yab,zab
      real(t_p) xcb,ycb,zcb
      real(t_p) xdb,ydb,zdb
      real(t_p) xad,yad,zad
      real(t_p) xcd,ycd,zcd
      real(t_p) rdb2,rad2,rcd2
      real(t_p) rab2,rcb2
      real(t_p) cc,ee,bkk2
      logical proceed
      logical header,huge
!$acc routine(image_acc) seq

c
c
c     zero out the out-of-plane bend energy and partitioning
c
      if(deb_Path) write(*,*) 'eopbend3gpu'
      call timer_enter( timer_eopbend3 )
      
      neopb = 0
      eopb  = 0.0_re_p
c     aeopb = 0.0_ti_p
      if (rank.eq.0) then
         header = .true.
      else
         header=.false.
      end if

!$acc parallel loop 
#ifdef USE_NVSHMEM_CUDA
!$acc&         present(x,y,z,loc,use,opbendglob,opbk
#else
!$acc&         present(x,y,z,loc,use,opbendglob,iang,opbk
#endif
!$acc&    ,iopb,aeopb) present(eopb) async
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
         ibloc = loc(ib)
         force = opbk(iopbend)
c
c     decide whether to compute the current interaction
c
         proceed = .true.

         if (proceed)  proceed = (use(ia) .or. use(ib) .or.
     &                              use(ic) .or. use(id))
c
c     get the coordinates of the atoms at trigonal center
c

         if (proceed) then
            xia = x(ia)
            yia = y(ia)
            zia = z(ia)
            xib = x(ib)
            yib = y(ib)
            zib = z(ib)
            xic = x(ic)
            yic = y(ic)
            zic = z(ic)
            xid = x(id)
            yid = y(id)
            zid = z(id)
c
c     compute the out-of-plane bending angle
c
            xab = xia - xib
            yab = yia - yib
            zab = zia - zib
            xcb = xic - xib
            ycb = yic - yib
            zcb = zic - zib
            xdb = xid - xib
            ydb = yid - yib
            zdb = zid - zib
            xad = xia - xid
            yad = yia - yid
            zad = zia - zid
            xcd = xic - xid
            ycd = yic - yid
            zcd = zic - zid
            if (use_polymer) then
               call image_inl(xab,yab,zab)
               call image_inl(xcb,ycb,zcb)
               call image_inl(xdb,ydb,zdb)
               call image_inl(xad,yad,zad)
               call image_inl(xcd,ycd,zcd)
            end if
c
c     W-D-C angle between A-B-C plane and B-D vector for D-B<AC
c
            if (opbtyp .eq. 'W-D-C') then
               rab2 = xab*xab + yab*yab + zab*zab
               rcb2 = xcb*xcb + ycb*ycb + zcb*zcb
               cc = rab2*rcb2 - (xab*xcb+yab*ycb+zab*zcb)**2
c
c     Allinger angle between A-C-D plane and D-B vector for D-B<AC
c
            else if (opbtyp .eq. 'ALLINGER') then
               rad2 = xad*xad + yad*yad + zad*zad
               rcd2 = xcd*xcd + ycd*ycd + zcd*zcd
               cc = rad2*rcd2 - (xad*xcd+yad*ycd+zad*zcd)**2
            end if
c
c     find the out-of-plane angle bending energy
c
            ee = xdb*(yab*zcb-zab*ycb) + ydb*(zab*xcb-xab*zcb)
     &              + zdb*(xab*ycb-yab*xcb)
            rdb2 = xdb*xdb + ydb*ydb + zdb*zdb
            if (rdb2.ne.0.0_ti_p .and. cc.ne.0.0_ti_p) then
               bkk2 = rdb2 - ee*ee/cc
               cosine = sqrt(bkk2/rdb2)
               cosine = min(1.0_ti_p,max(-1.0_ti_p,cosine))
               angle1 = radian * acos(cosine)
               dt = angle1
               dt2 = dt * dt
               dt3 = dt2 * dt
               dt4 = dt2 * dt2
               e = opbunit * force * dt2
     &                * (1.0_ti_p+copb*dt+qopb*dt2+popb*dt3+sopb*dt4)
c
c     increment the total out-of-plane bending energy
c
               neopb = neopb + 1
               eopb = eopb + e
!$acc atomic
               aeopb(ibloc) = aeopb(ibloc) + e
#ifndef _OPENACC
c
c     print a message if the energy of this interaction is large
c
               huge = (e .gt. 2.0_ti_p)
               if (debug .or. (verbose.and.huge)) then
                  if (header) then
                     header = .false.
                     write (iout,10)
   10                format (/,' Individual Out-of-Plane Bend',
     &                          ' Interactions :',
     &                       //,' Type',25x,'Atom Names',21x,'Angle',
     &                          6x,'Energy',/)
                  end if
                  write (iout,20)  id,name(id),ib,name(ib),ia,
     &                             name(ia),ic,name(ic),angle1,e
   20             format (' O-P-Bend',2x,4(i7,'-',a3),f11.4,f12.4)
               end if
#endif
            end if
         end if
      end do
!$acc update self(aeopb) async
      call timer_exit( timer_eopbend3 )
      end
