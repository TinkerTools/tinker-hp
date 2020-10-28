c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine eopbend1  --  out-of-plane energy and derivs  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "eopbend1" computes the out-of-plane bend potential energy and
c     first derivatives at trigonal centers via a Wilson-Decius-Cross
c     or Allinger angle
c
c
#include "tinker_precision.h"
      module eopbend1gpu_inl
        contains
#include "image.f.inc"
      end module

      subroutine eopbend1gpu
      use angle
      use angpot
      use atmlst
      use atoms
      use bound
      use deriv
      use domdec
      use energi
      use eopbend1gpu_inl
      use group
      use inform    ,only: deb_Path
      use math
      use opbend
      use tinheader ,only:ti_p,re_p
      use usage
      use virial
      use timestat
      implicit none
      integer i,iopbend,iopbendloc
      integer ia,ib,ic,id
      integer ialoc,ibloc,icloc,idloc
#ifdef USE_NVSHMEM_CUDA
      integer ipe,ind
#endif
      real(t_p) e,angle1,force
      real(t_p) dot,cosine
      real(t_p) cc,ee,bkk2,term
      real(t_p) deddt,dedcos
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
      real(t_p) dccdxia,dccdyia,dccdzia
      real(t_p) dccdxic,dccdyic,dccdzic
      real(t_p) dccdxid,dccdyid,dccdzid
      real(t_p) deedxia,deedyia,deedzia
      real(t_p) deedxic,deedyic,deedzic
      real(t_p) deedxid,deedyid,deedzid
      real(t_p) dedxia,dedyia,dedzia
      real(t_p) dedxib,dedyib,dedzib
      real(t_p) dedxic,dedyic,dedzic
      real(t_p) dedxid,dedyid,dedzid
      real(r_p) dedi(3)
      logical proceed
!$acc routine(image_acc) seq   
c
c     calculate the out-of-plane bending energy and derivatives
c
      if(deb_Path) write(*,*) 'eopbend1gpu'
      call timer_enter( timer_eopbend1 )

!$acc parallel loop private(dedi) async
#ifdef USE_NVSHMEM_CUDA
!$acc&      present(loc,use,opbendglob,opbk,iopb
#else
!$acc&      present(loc,use,opbendglob,iang,opbk,iopb
#endif
!$acc&     ,deopb,eopb,g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz)
      do iopbendloc = 1, nopbendloc
         iopbend = opbendglob(iopbendloc)
         i  = iopb(iopbend)
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
         ialoc = loc(ia)
         ibloc = loc(ib)
         icloc = loc(ic)
         idloc = loc(id)
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
               call image_inl (xab,yab,zab)
               call image_inl (xcb,ycb,zcb)
               call image_inl (xdb,ydb,zdb)
               call image_inl (xad,yad,zad)
               call image_inl (xcd,ycd,zcd)
            end if
c
c     W-D-C angle between A-B-C plane and B-D vector for D-B<AC
c
            if (opbtyp .eq. 'W-D-C') then
               rab2 = xab*xab + yab*yab + zab*zab
               rcb2 = xcb*xcb + ycb*ycb + zcb*zcb
               dot = xab*xcb+yab*ycb+zab*zcb
               cc = rab2*rcb2 - dot*dot
c
c     Allinger angle between A-C-D plane and D-B vector for D-B<AC
c
            else if (opbtyp .eq. 'ALLINGER') then
               rad2 = xad*xad + yad*yad + zad*zad
               rcd2 = xcd*xcd + ycd*ycd + zcd*zcd
               dot = xad*xcd + yad*ycd + zad*zcd
               cc = rad2*rcd2 - dot*dot
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
               deddt = opbunit * force * dt * radian
     &                * (2.0_ti_p + 3.0_ti_p*copb*dt + 4.0_ti_p*qopb*dt2
     &                 + 5.0_ti_p*popb*dt3 + 6.0_ti_p*sopb*dt4)
               dedcos = -deddt * sign(1.0_ti_p,ee) / sqrt(cc*bkk2)
c
c     chain rule terms for first derivative components
c
               if (opbtyp .eq. 'W-D-C') then
                  term = ee / cc
                  dccdxia = (xab*rcb2-xcb*dot) * term
                  dccdyia = (yab*rcb2-ycb*dot) * term
                  dccdzia = (zab*rcb2-zcb*dot) * term
                  dccdxic = (xcb*rab2-xab*dot) * term
                  dccdyic = (ycb*rab2-yab*dot) * term
                  dccdzic = (zcb*rab2-zab*dot) * term
                  dccdxid = 0.0_ti_p
                  dccdyid = 0.0_ti_p
                  dccdzid = 0.0_ti_p
               else if (opbtyp .eq. 'ALLINGER') then
                  term = ee / cc
                  dccdxia = (xad*rcd2-xcd*dot) * term
                  dccdyia = (yad*rcd2-ycd*dot) * term
                  dccdzia = (zad*rcd2-zcd*dot) * term
                  dccdxic = (xcd*rad2-xad*dot) * term
                  dccdyic = (ycd*rad2-yad*dot) * term
                  dccdzic = (zcd*rad2-zad*dot) * term
                  dccdxid = -dccdxia - dccdxic
                  dccdyid = -dccdyia - dccdyic
                  dccdzid = -dccdzia - dccdzic
               end if
               term = ee / rdb2
               deedxia = ydb*zcb - zdb*ycb
               deedyia = zdb*xcb - xdb*zcb
               deedzia = xdb*ycb - ydb*xcb
               deedxic = yab*zdb - zab*ydb
               deedyic = zab*xdb - xab*zdb
               deedzic = xab*ydb - yab*xdb
               deedxid = ycb*zab - zcb*yab + xdb*term
               deedyid = zcb*xab - xcb*zab + ydb*term
               deedzid = xcb*yab - ycb*xab + zdb*term
c
c     compute first derivative components for this angle
c
               dedxia = dedcos * (dccdxia+deedxia)
               dedyia = dedcos * (dccdyia+deedyia)
               dedzia = dedcos * (dccdzia+deedzia)
               dedxic = dedcos * (dccdxic+deedxic)
               dedyic = dedcos * (dccdyic+deedyic)
               dedzic = dedcos * (dccdzic+deedzic)
               dedxid = dedcos * (dccdxid+deedxid)
               dedyid = dedcos * (dccdyid+deedyid)
               dedzid = dedcos * (dccdzid+deedzid)
               dedxib = -dedxia - dedxic - dedxid
               dedyib = -dedyia - dedyic - dedyid
               dedzib = -dedzia - dedzic - dedzid
c
c     increment the out-of-plane bending energy and gradient
c
               eopb = eopb + e

               dedi(1) = dedxib
               dedi(2) = dedyib
               dedi(3) = dedzib
!$acc atomic update     
               deopb(1,ibloc) = deopb(1,ibloc) + dedi(1)
!$acc atomic update     
               deopb(2,ibloc) = deopb(2,ibloc) + dedi(2)
!$acc atomic update     
               deopb(3,ibloc) = deopb(3,ibloc) + dedi(3)
c
               dedi(1) = dedxia
               dedi(2) = dedyia
               dedi(3) = dedzia
!$acc atomic update
               deopb(1,ialoc) = deopb(1,ialoc) + dedi(1)
!$acc atomic update
               deopb(2,ialoc) = deopb(2,ialoc) + dedi(2)
!$acc atomic update
               deopb(3,ialoc) = deopb(3,ialoc) + dedi(3)
c
               dedi(1) = dedxic
               dedi(2) = dedyic
               dedi(3) = dedzic
!$acc atomic update
               deopb(1,icloc) = deopb(1,icloc) + dedi(1)
!$acc atomic update
               deopb(2,icloc) = deopb(2,icloc) + dedi(2)
!$acc atomic update
               deopb(3,icloc) = deopb(3,icloc) + dedi(3)
c
               dedi(1) = dedxid
               dedi(2) = dedyid
               dedi(3) = dedzid
!$acc atomic update
               deopb(1,idloc) = deopb(1,idloc) + dedi(1)
!$acc atomic update
               deopb(2,idloc) = deopb(2,idloc) + dedi(2)
!$acc atomic update
               deopb(3,idloc) = deopb(3,idloc) + dedi(3)
c
c     increment the internal virial tensor components
c
               g_vxx = g_vxx + xab*dedxia + xcb*dedxic + xdb*dedxid
               g_vxy = g_vxy + yab*dedxia + ycb*dedxic + ydb*dedxid
               g_vxz = g_vxz + zab*dedxia + zcb*dedxic + zdb*dedxid
               g_vyy = g_vyy + yab*dedyia + ycb*dedyic + ydb*dedyid
               g_vyz = g_vyz + zab*dedyia + zcb*dedyic + zdb*dedyid
               g_vzz = g_vzz + zab*dedzia + zcb*dedzic + zdb*dedzid
            end if
         end if
      end do

      call timer_exit( timer_eopbend1 )
      end
