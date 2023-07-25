
c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     "empole1" : driver for calculation of the multipole and dipole polarization
c     energy and derivatives with respect to Cartesian coordinates
c
c
#include "tinker_macro.h"
      module empole1gpu_inl
        use tinTypes , only: real3,real3_red,rpole_elt
        logical:: em1c_fi=.true.
        ener_rtyp emself
        include "erfcore_data.f.inc"
#include "atomicOp.h.f"
        contains
#include "convert.f.inc"
#include "image.f.inc"
#if defined(SINGLE) | defined(MIXED)
        include "erfcscore.f.inc"
#else
        include "erfcdcore.f.inc"
#endif
#include "groups.inc.f"
#include "atomicOp.inc.f"
#include "switch_respa.f.inc"
#include "pair_mpole1.f.inc"
      end module

      subroutine empole1gpu
      use domdec
      use energi
      use potent
      use tinheader ,only:ti_p,re_p
      use group
      implicit none

      if (use_group) then
         if (wgrp(1,2).eq.0._ti_p) return
      end if
c
c     choose the method for summing over multipole interactions
c
      if (use_lambdadyn) then
         call elambdampole1cgpu
      else
         call empole1cgpu
      end if
c
c     zero out energy and derivative terms which are not in use
c
      if (.not. use_mpole) then
         em = 0.0_ti_p
      end if
      if (.not. use_polar) then
         ep = 0.0_ti_p
      end if
      end
c
c     "elec_calc_reset" reset global data for electrostatics calculations
c
      subroutine elec_calc_reset
      use atmlst
      use cutoff
      use domdec  ,only: ndir,nproc
      use inform
      use interfaces,only: reorder_nblist,bspline_fill_sitegpu
      use mpole   ,only: ipole,npolelocnl
      use neigh
      use potent
      use shunt   ,only: off2
      use utilgpu
      use timestat
      implicit none
      integer i
c
c     Communicate poleglob
c
#ifdef _OPENACC
      call orderPole
#endif
      if (mlst2_enable) call set_ElecData_cellOrder(.false.)
c
c     check the sign of multipole components at chiral sites
c
      call timer_enter( timer_other )
      call chkpolegpu(.false.)
c
c     rotate the multipole components into the global frame
c
      call rotpolegpu
#ifdef _OPENACC
c
c     Reorder Atom-Atom neigbhor list
c
      if ((.not.(use_pmecore)).or.(use_pmecore).and.(rank.le.ndir-1))
     &   then
         if (use_mreal) then
            if (mlst_enable.and.
     &         (use_mpoleshortreal.or.use_polarshortreal)) then
               call switch('SHORTEWALD')
               !TODO Uncomment me
               call reorder_nblist(shortelst,nshortelst,nshortelstc
     &                     ,npolelocnl,off2,ipole,poleglobnl)
c!$acc parallel loop async(dir_queue) default(present)
c              do i = 1,npolelocnl
c                 nshortelstc(i) = nshortelst(i)
c              end do
            end if
            if (mlst_enable) then
               call switch('EWALD')
               !TODO Uncomment me
               call reorder_nblist(elst,nelst,nelstc
     &                     ,npolelocnl,off2,ipole,poleglobnl)
c!$acc parallel loop async(dir_queue) default(present)
c              do i = 1,npolelocnl
c                 nelstc(i) = nelst(i)
c              end do
            end if
         end if
      end if
#else
      do i = 1,npolelocnl
         if (use_mpoleshortreal.or.use_polarshortreal)
     &      nshortelstc(i) = nshortelst(i)
         nelstc(i) = nelst(i)
      end do
#endif
c
c     compute B-spline coefficients
c
      if (use_mrec.or.use_prec) then
         call bspline_fill_sitegpu
      end if

#ifdef _OPENACC
      if (dir_queue.ne.rec_queue) then
         if (app_id.eq.dynamic_a.or.app_id.eq.pimd_a)
     &      call end_dir_stream_cover
         call start_dir_stream_cover
      end if
#endif

      call timer_exit ( timer_other,quiet_timers )
      end subroutine
c
c
c     "empole1c" calculates the multipole energy and derivatives
c     with respect to Cartesian coordinates using particle mesh Ewald
c     summation and a neighbor list
c
c
      subroutine empole1cgpu
      use sizes
      use atmlst
      use atoms
      use boxes
      use chgpot
      use cflux
      use deriv
      use domdec
      use elec_wspace,only: trq=>r2Work4
      use empole1gpu_inl
      use energi
      use ewald
      use inform     ,only: deb_Path, minmaxone
      use interfaces ,only: torquegpu,emreal1c_p
      use math
      use mpole
      use mpi
      use neigh
      use potent
      use shunt
      use tinheader  ,only: zeror
      use timestat
      use tinMemory  ,only: prmem_request
      use utilcomm
      use utilgpu    ,pot=>ug_workS_r
      use virial
      implicit none
      integer i,j,ii
      integer iipole,iglob,ierr
      real(t_p) zero,one,two,three,half
      real(t_p) e,f
      real(t_p) term,fterm
      real(t_p) cii,dii,qii
      real(t_p) xd,yd,zd
      real(t_p) xq,yq,zq
      real(t_p) xv,yv,zv,vterm
      real(t_p) ci,dix,diy,diz
      real(t_p) qixx,qixy,qixz
      real(t_p) qiyy,qiyz,qizz
      real(t_p) xdfield,ydfield
      real(t_p) zdfield
      real(t_p) time0,time1
      parameter(zero=0.0_ti_p,  one=1.0_ti_p,
     &           two=2.0_ti_p,three=3.0_ti_p,
     &          half=0.5_ti_p)
c
      if (npole .eq. 0)  return
      if (em1c_fi) then
         em1c_fi = .false.
!$acc enter data create(emself)
      end if
c
      if (deb_Path) write(*,*) 'empole1cgpu'
c
c     set Ewald coefficient
c
      aewald = aeewald
c
c     zero out the atomic multipole energy and derivatives
c
!$acc serial async present(emself,emrec,em)
      em     = zero
      emself = 0
      emrec  = zero
!$acc end serial
c
c     set the energy unit conversion factor
c
      f = electric / dielec
c
c     Reset global data for electrostatic
c
      call elec_calc_reset
c
c     compute the reciprocal space part of the Ewald summation
c
      if ((.not.(use_pmecore)).or.(use_pmecore).and.(rank.le.ndir-1))
     &   then
         call timer_enter( timer_real )
         if (use_mreal) then
            call emreal1c_p
         end if

         if (use_mself) then
            if(deb_Path) print*, 'emself'
c
c     compute the Ewald self-energy term over all the atoms
c
         term  = two * aewald * aewald
         fterm = -f * aewald / sqrtpi
!$acc parallel loop async(def_queue) default(present) present(emself)
         do i = 1, npoleloc
            iipole = poleglob(i)
            ci     = rpole( 1,iipole)
            dix    = rpole( 2,iipole)
            diy    = rpole( 3,iipole)
            diz    = rpole( 4,iipole)
            qixx   = rpole( 5,iipole)
            qixy   = rpole( 6,iipole)
            qixz   = rpole( 7,iipole)
            qiyy   = rpole( 9,iipole)
            qiyz   = rpole(10,iipole)
            qizz   = rpole(13,iipole)
            cii    = ci*ci
            dii    =       dix*dix  +  diy*diy  + diz*diz
            qii    = two*(qixy*qixy + qixz*qixz + qiyz*qiyz)
     &                  + qixx*qixx + qiyy*qiyy + qizz*qizz
            e      = fterm*(cii + term*(dii/three +
     &                                  two*term*qii/5.0_ti_p))
            emself = emself + tp2enr(e)
            if (use_chgflx) pot(loc(ipole(iipole))) = 2.0 * fterm * ci
         end do
c
c     modify gradient and virial for charge flux self-energy
c
         if (use_chgflx) then
            call commDDd_ext(pot,1,ucComm+ucNeig)
            call dcflux2(pot,de_ws1)
            call adflux2(de_ws1,dem)
            call mem_set(pot,zeror,int(nbloc,mipk),def_stream)
         end if
c
c     compute the cell dipole boundary correction term
c
         if (boundary .eq. 'VACUUM') then
            write(0,*) 'FATAL ERROR VACUUM boundary unavailable !!!'
            __TINKER_FATAL__

            call prmem_request(trq,3,npoleloc,queue=def_queue)
c
!$acc data create(trq)
!$acc wait
            xd = zero
            yd = zero
            zd = zero
!$acc parallel loop async(def_queue) default(present)
            do i = 1, npoleloc
               iipole = poleglob(i)
               iglob  = ipole(iipole)
               xd = xd + rpole(2,iipole) + rpole(1,iipole)*x(iglob)
               yd = yd + rpole(3,iipole) + rpole(1,iipole)*y(iglob)
               zd = zd + rpole(4,iipole) + rpole(1,iipole)*z(iglob)
            end do
!$acc wait(def_queue)
            call MPI_ALLREDUCE(MPI_IN_PLACE,xd,1,MPI_TPREC,MPI_SUM,
     $         COMM_TINKER,ierr)
            call MPI_ALLREDUCE(MPI_IN_PLACE,yd,1,MPI_TPREC,MPI_SUM,
     $         COMM_TINKER,ierr)
            call MPI_ALLREDUCE(MPI_IN_PLACE,zd,1,MPI_TPREC,MPI_SUM,
     $         COMM_TINKER,ierr)
            if (rank.eq.0) then
!$acc serial async(def_queue) present(em)
              em = em + term*(xd*xd+yd*yd+zd*zd)
!$acc end serial
            end if
            term    = (two/3.0_ti_p) * f * (pi/volbox)
!$acc serial async(def_queue) present(emself)
            emself   = emself + tp2enr(term*(xd*xd+yd*yd+zd*zd))
!$acc end serial
            xdfield = -two * term * xd
            ydfield = -two * term * yd
            zdfield = -two * term * zd
!$acc parallel loop async(def_queue) default(present)
            do ii = 1, npoleloc
               iipole   = poleglob(ii)
               iglob    = ipole(iipole)
               i        = loc(iglob)
               dem(1,i) = dem(1,i) + two*term*rpole(1,iipole)*xd
               dem(2,i) = dem(2,i) + two*term*rpole(1,iipole)*yd
               dem(3,i) = dem(3,i) + two*term*rpole(1,iipole)*zd

               trq(1,ii)=rpole(3,iipole)*zdfield-rpole(4,iipole)*ydfield
               trq(2,ii)=rpole(4,iipole)*xdfield-rpole(2,iipole)*zdfield
               trq(3,ii)=rpole(2,iipole)*ydfield-rpole(3,iipole)*xdfield
            end do
            ! TODO Remove this comment 
            ! ---( Trouble in fixed precision )---
            !call torquegpu(npoleloc,nbloc,poleglob,loc,trq,dem,def_queue)
c
c     boundary correction to virial due to overall cell dipole
c
            xd = zero
            yd = zero
            zd = zero
            xq = zero
            yq = zero
            zq = zero
!$acc parallel loop async(def_queue) default(present)
            do i = 1, npoleloc
               iipole = poleglob(i)
               iglob  = ipole(iipole)
               xd     = xd + rpole(2,iipole)
               yd     = yd + rpole(3,iipole)
               zd     = zd + rpole(4,iipole)
               xq     = xq + rpole(1,iipole)*x(iglob)
               yq     = yq + rpole(1,iipole)*y(iglob)
               zq     = zq + rpole(1,iipole)*z(iglob)
            end do
!$acc wait(def_queue)
            call MPI_ALLREDUCE(MPI_IN_PLACE,xd,1,MPI_TPREC,MPI_SUM,
     $         COMM_TINKER,ierr)
            call MPI_ALLREDUCE(MPI_IN_PLACE,yd,1,MPI_TPREC,MPI_SUM,
     $         COMM_TINKER,ierr)
            call MPI_ALLREDUCE(MPI_IN_PLACE,zd,1,MPI_TPREC,MPI_SUM,
     $         COMM_TINKER,ierr)
            call MPI_ALLREDUCE(MPI_IN_PLACE,xq,1,MPI_TPREC,MPI_SUM,
     $         COMM_TINKER,ierr)
            call MPI_ALLREDUCE(MPI_IN_PLACE,yq,1,MPI_TPREC,MPI_SUM,
     $         COMM_TINKER,ierr)
            call MPI_ALLREDUCE(MPI_IN_PLACE,zq,1,MPI_TPREC,MPI_SUM,
     $         COMM_TINKER,ierr)
            if (rank.eq.0) then
!$acc serial default(present)
            xv = xd * xq
            yv = yd * yq
            zv = zd * zq
            vterm = term * (xd*xd + yd*yd + zd*zd + two*(xv+yv+zv)
     &                    + xq*xq + yq*yq + zq*zq)
            vir(1,1) = vir(1,1) + two*term*(xq*xq+xv) + vterm
            vir(2,1) = vir(2,1) + two*term*(xq*yq+xv)
            vir(3,1) = vir(3,1) + two*term*(xq*zq+xv)
            vir(1,2) = vir(1,2) + two*term*(yq*xq+yv)
            vir(2,2) = vir(2,2) + two*term*(yq*yq+yv) + vterm
            vir(3,2) = vir(3,2) + two*term*(yq*zq+yv)
            vir(1,3) = vir(1,3) + two*term*(zq*xq+zv)
            vir(2,3) = vir(2,3) + two*term*(zq*yq+zv)
            vir(3,3) = vir(3,3) + two*term*(zq*zq+zv) + vterm
!$acc end serial
            end if
!$acc wait
!$acc end data
c
           end if !( boundary .eq. "VACUUM" )
         end if !( use_mself )
         call timer_exit( timer_real,quiet_timers )

      end if !((.not.(use_pmecore)).or.(use_pmecore).and.(rank.le.ndir-1))

c
c     compute the reciprocal space part of the Ewald summation
c
      if ((.not.(use_pmecore)).or.(use_pmecore).and.(rank.gt.ndir-1))
     &   then
         if (use_mrec) then
         call timer_enter( timer_rec )
         call emrecip1gpu
         call timer_exit( timer_rec,quiet_timers )
         end if
      end if
c
#ifdef _OPENACC
      !Finalize async overlapping
      if (dir_queue.ne.rec_queue) call end_dir_stream_cover
#endif
c
c     Add both contribution to the energy
c
!$acc serial async(rec_queue) present(em,emself,em_r,emrec)
      em = em + enr2en(emself + em_r) + emrec
!$acc end serial
c
      end
c
c
c     "elambdampole1c" calculates the multipole energy and derivatives
c     with respect to Cartesian coordinates using particle mesh Ewald
c     summation and a neighbor list during lambda dynamics
c
c
      subroutine elambdampole1cgpu
      use sizes
      use atmlst
      use atoms
      use boxes
      use chgpot
      use cflux
      use deriv
      use domdec
      use elec_wspace,only: trq=>r2Work4
      use empole1gpu_inl
      use energi
      use ewald
      use inform     ,only: deb_Path
      use interfaces ,only: torquegpu,emreal1c_p
      use mutant
      use math
      use mpole
      use mpi
      use neigh
      use potent
      use shunt
      use timestat
      use tinheader  ,only: ti_p,zerom,zeror
      use tinMemory  ,only: prmem_request,mipk
      use utilcomm
      use utilgpu    ,pot=>ug_workS_r
      use virial
      implicit none
      integer i,j,ii
      integer iipole,iglob,ierr,altopt
      integer(mipk) siz8
      real(t_p) zero,one,two,three,half
      real(t_p) e,f,elambdatemp
      real(t_p) term,fterm
      real(t_p) cii,dii,qii
      real(t_p) xd,yd,zd,xdt,ydt,zdt
      real(t_p) xq,yq,zq
      real(t_p) xv,yv,zv,vterm
      real(t_p) ci,dix,diy,diz
      real(t_p) qixx,qixy,qixz
      real(t_p) qiyy,qiyz,qizz
      real(t_p) xdfield,ydfield
      real(t_p) zdfield
      real(t_p) time0,time1
      real(r_p),allocatable:: delambdarec0(:,:),delambdarec1(:,:)
      real(r_p) elambdarec0,elambdarec1
      real(r_p) :: g_vxx_temp,g_vxy_temp,g_vxz_temp
      real(r_p) :: g_vyy_temp,g_vyz_temp,g_vzz_temp
      real(r_p) :: g_vxx_1,g_vxy_1,g_vxz_1
      real(r_p) :: g_vyy_1,g_vyz_1,g_vzz_1
      real(r_p) :: g_vxx_0,g_vxy_0,g_vxz_0
      real(r_p) :: g_vyy_0,g_vyz_0,g_vzz_0
      parameter(zero=0.0_ti_p,  one=1.0_ti_p
     &         , two=2.0_ti_p,three=3.0_ti_p
     &         ,half=0.5_ti_p
#ifdef _OPENACC
     &         ,altopt = 0
#else
     &         ,altopt = 1
#endif
     &         )
c
      if (npole .eq. 0)  return
      if(deb_Path) write(*,*) 'elambdampole1cgpu'

      if (em1c_fi) then
         em1c_fi = .false.
!$acc enter data create(emself)
      end if
c
c     set Ewald coefficient
c
      aewald = aeewald

      allocate (delambdarec0(3,nlocrec2))
      allocate (delambdarec1(3,nlocrec2))
!$acc enter data create(elambdarec0,elambdarec1
!$acc&     ,delambdarec0,delambdarec1
!$acc&     ,g_vxx_temp,g_vxy_temp,g_vxz_temp
!$acc&     ,g_vyy_temp,g_vyz_temp,g_vzz_temp) async(rec_queue)
      elambdatemp = elambda  
c
c     zero out the atomic multipole energy and derivatives
c
!$acc serial async present(emself,emrec,em,delambdae)
      em        = zero
      emself     = 0
      emrec     = zero
      delambdae = 0
!$acc end serial
c
c     set the energy unit conversion factor
c
      f = electric / dielec
c
c     Reset global data for electrostatic
c
      call elec_calc_reset
c
c     compute the reciprocal space part of the Ewald summation
c
      if ((.not.(use_pmecore)).or.(use_pmecore).and.(rank.le.ndir-1))
     &   then
         call timer_enter( timer_real )
         if (use_mreal) then
#ifdef _OPENACC
            call emreal1c_p
#else
            if (use_chgpen) then
               call emreal1c
            else
               call emreal1c_p
            end if
#endif
         end if

         if (use_mself) then
            if(deb_Path) print*, 'emself'
c
c     compute the Ewald self-energy term over all the atoms
c
         term  = two * aewald * aewald
         fterm = -f * aewald / sqrtpi
!$acc parallel loop async(def_queue) default(present)
!$acc&         present(emself,delambdae)
!$acc&         reduction(+:emself,delambdae)
         do i = 1, npoleloc
            iipole = poleglob(i)
             iglob = ipole(iipole)
            ci     = rpole( 1,iipole)
            dix    = rpole( 2,iipole)
            diy    = rpole( 3,iipole)
            diz    = rpole( 4,iipole)
            qixx   = rpole( 5,iipole)
            qixy   = rpole( 6,iipole)
            qixz   = rpole( 7,iipole)
            qiyy   = rpole( 9,iipole)
            qiyz   = rpole(10,iipole)
            qizz   = rpole(13,iipole)
            cii    = ci*ci
            dii    =       dix*dix  +  diy*diy  + diz*diz
            qii    = two*(qixy*qixy + qixz*qixz + qiyz*qiyz)
     &                  + qixx*qixx + qiyy*qiyy + qizz*qizz
            e      = fterm*(cii + term*(dii/three +
     &                                  two*term*qii/5.0_ti_p))
            emself  = emself + tp2enr(e)
            if (mut(iglob).and.elambda.gt.0.0) then
               delambdae = delambdae + 2.0*e/elambda
            end if
            if (use_chgflx) pot(loc(iglob)) = 2.0 * fterm * ci
         end do
c
c     modify gradient and virial for charge flux self-energy
c
         if (use_chgflx) then
            call commDDd_add(pot,1,ucComm+ucBig)
            call commDDd_ext(pot,1,ucComm+ucNeig)
            call dcflux2(pot,de_ws1)
            call adflux2(de_ws1,dem)
            call mem_set(pot,zeror,int(nbloc,mipk),def_stream)
         end if
c
c     compute the cell dipole boundary correction term
c
         if (boundary .eq. 'VACUUM') then
            write(0,*) 'FATAL ERROR VACUUM boundary unavailable !!!'
            __TINKER_FATAL__

            call prmem_request(trq,3,npoleloc,queue=def_queue)
c
!$acc data create(trq)
!$acc wait
            xd  = zero
            yd  = zero
            zd  = zero
            xdt = zero
            ydt = zero
            zdt = zero
            !FIXME complete delambdae integration
!$acc parallel loop async(def_queue) default(present)
            do i = 1, npoleloc
               iipole = poleglob(i)
               iglob  = ipole(iipole)
               xd = xd + rpole(2,iipole) + rpole(1,iipole)*x(iglob)
               yd = yd + rpole(3,iipole) + rpole(1,iipole)*y(iglob)
               zd = zd + rpole(4,iipole) + rpole(1,iipole)*z(iglob)
               if (mut(iglob).and.elambda.gt.0) then
                  xdt= xdt+ (rpole(2,iipole) + rpole(1,iipole)*x(iglob))
     &                      /elambda
                  ydt= ydt+ (rpole(3,iipole) + rpole(1,iipole)*x(iglob))
     &                      /elambda
                  zdt= zdt+ (rpole(4,iipole) + rpole(1,iipole)*x(iglob))
     &                      /elambda
               end if
            end do
!$acc wait(def_queue)
            if (ndir.gt.1) then
            call MPI_ALLREDUCE(MPI_IN_PLACE,xd,1,MPI_TPREC,MPI_SUM,
     $         COMM_TINKER,ierr)
            call MPI_ALLREDUCE(MPI_IN_PLACE,yd,1,MPI_TPREC,MPI_SUM,
     $         COMM_TINKER,ierr)
            call MPI_ALLREDUCE(MPI_IN_PLACE,zd,1,MPI_TPREC,MPI_SUM,
     $         COMM_TINKER,ierr)
            call MPI_ALLREDUCE(MPI_IN_PLACE,xdt,1,MPI_TPREC,MPI_SUM,
     $         COMM_TINKER,ierr)
            call MPI_ALLREDUCE(MPI_IN_PLACE,ydt,1,MPI_TPREC,MPI_SUM,
     $         COMM_TINKER,ierr)
            call MPI_ALLREDUCE(MPI_IN_PLACE,zdt,1,MPI_TPREC,MPI_SUM,
     $         COMM_TINKER,ierr)
            end if
            if (rank.eq.0) then
            term    = (two/3.0_ti_p) * f * (pi/volbox)
!$acc serial async(def_queue) present(emself)
            emself   = emself + tp2enr(term*(xd*xd+yd*yd+zd*zd))
!$acc end serial
            end if
            xdfield = -two * term * xd
            ydfield = -two * term * yd
            zdfield = -two * term * zd
!$acc parallel loop async(def_queue) default(present)
            do ii = 1, npoleloc
               iipole   = poleglob(ii)
               iglob    = ipole(iipole)
               i        = loc(iglob)
               dem(1,i) = dem(1,i) + two*term*rpole(1,iipole)*xd
               dem(2,i) = dem(2,i) + two*term*rpole(1,iipole)*yd
               dem(3,i) = dem(3,i) + two*term*rpole(1,iipole)*zd

               trq(1,ii)=rpole(3,iipole)*zdfield-rpole(4,iipole)*ydfield
               trq(2,ii)=rpole(4,iipole)*xdfield-rpole(2,iipole)*zdfield
               trq(3,ii)=rpole(2,iipole)*ydfield-rpole(3,iipole)*xdfield
            end do
            ! TODO Remove this comment 
            ! ---( Trouble in fixed precision )---
            !call torquegpu(npoleloc,nbloc,poleglob,loc,trq,dem,def_queue)
c
c     boundary correction to virial due to overall cell dipole
c
            xd = zero
            yd = zero
            zd = zero
            xq = zero
            yq = zero
            zq = zero
!$acc parallel loop async(def_queue) default(present)
            do i = 1, npoleloc
               iipole = poleglob(i)
               iglob  = ipole(iipole)
               xd     = xd + rpole(2,iipole)
               yd     = yd + rpole(3,iipole)
               zd     = zd + rpole(4,iipole)
               xq     = xq + rpole(1,iipole)*x(iglob)
               yq     = yq + rpole(1,iipole)*y(iglob)
               zq     = zq + rpole(1,iipole)*z(iglob)
            end do
!$acc wait(def_queue)
            if (ndir.gt.1) then
            call MPI_ALLREDUCE(MPI_IN_PLACE,xd,1,MPI_TPREC,MPI_SUM,
     $         COMM_TINKER,ierr)
            call MPI_ALLREDUCE(MPI_IN_PLACE,yd,1,MPI_TPREC,MPI_SUM,
     $         COMM_TINKER,ierr)
            call MPI_ALLREDUCE(MPI_IN_PLACE,zd,1,MPI_TPREC,MPI_SUM,
     $         COMM_TINKER,ierr)
            call MPI_ALLREDUCE(MPI_IN_PLACE,xq,1,MPI_TPREC,MPI_SUM,
     $         COMM_TINKER,ierr)
            call MPI_ALLREDUCE(MPI_IN_PLACE,yq,1,MPI_TPREC,MPI_SUM,
     $         COMM_TINKER,ierr)
            call MPI_ALLREDUCE(MPI_IN_PLACE,zq,1,MPI_TPREC,MPI_SUM,
     $         COMM_TINKER,ierr)
            end if
            if (rank.eq.0) then
!$acc kernels default(present)
            xv = xd * xq
            yv = yd * yq
            zv = zd * zq
            vterm = term * (xd*xd + yd*yd + zd*zd + two*(xv+yv+zv)
     &                    + xq*xq + yq*yq + zq*zq)
            vir(1,1) = vir(1,1) + two*term*(xq*xq+xv) + vterm
            vir(2,1) = vir(2,1) + two*term*(xq*yq+xv)
            vir(3,1) = vir(3,1) + two*term*(xq*zq+xv)
            vir(1,2) = vir(1,2) + two*term*(yq*xq+yv)
            vir(2,2) = vir(2,2) + two*term*(yq*yq+yv) + vterm
            vir(3,2) = vir(3,2) + two*term*(yq*zq+yv)
            vir(1,3) = vir(1,3) + two*term*(zq*xq+zv)
            vir(2,3) = vir(2,3) + two*term*(zq*yq+zv)
            vir(3,3) = vir(3,3) + two*term*(zq*zq+zv) + vterm
!$acc end kernels
            end if
!$acc wait
!$acc end data
c
           end if !( boundary .eq. "VACUUM" )
         end if !( use_mself )
         call timer_exit( timer_real,quiet_timers )

      end if !((.not.(use_pmecore)).or.(use_pmecore).and.(rank.le.ndir-1))
c
c     compute the reciprocal space part of the Ewald summation
c
      if ((.not.(use_pmecore)).or.(use_pmecore).and.(rank.gt.ndir-1))
     &   then
         if (use_mrec) then
         call timer_enter( timer_rec )
c
c        the reciprocal part is interpolated between 0 and 1
c
         elambda = 0.0
         siz8    = 3*nlocrec2
         call altelec(altopt)
         call rotpolegpu
!$acc serial async(rec_queue)
!$acc& present(emrec,g_vxx_temp,g_vxy_temp,g_vxz_temp,
!$acc&  g_vyy_temp,g_vyz_temp,g_vzz_temp,
!$acc& g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz)
         emrec   = 0.0
         g_vxx_temp = g_vxx
         g_vxy_temp = g_vxy
         g_vxz_temp = g_vxz
         g_vyy_temp = g_vyy
         g_vyz_temp = g_vyz
         g_vzz_temp = g_vzz
         g_vxx = 0.0
         g_vxy = 0.0
         g_vxz = 0.0
         g_vyy = 0.0
         g_vyz = 0.0
         g_vzz = 0.0
!$acc end serial
         call mem_set(demrec,zerom,siz8,rec_stream)
         call emrecip1gpu
!$acc serial async(rec_queue) present(elambdarec0,emrec,
!$acc& g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz)
         elambdarec0 = emrec
         g_vxx_0 = g_vxx
         g_vxy_0 = g_vxy
         g_vxz_0 = g_vxz
         g_vyy_0 = g_vyy
         g_vyz_0 = g_vyz
         g_vzz_0 = g_vzz
         g_vxx = 0.0
         g_vxy = 0.0
         g_vxz = 0.0
         g_vyy = 0.0
         g_vyz = 0.0
         g_vzz = 0.0
         emrec       = 0.0
!$acc end serial
         call mem_move(delambdarec0,demrec,siz8,rec_stream)
         call mem_set(demrec,zerom,siz8,rec_stream)

         elambda = 1.0
         call altelec(altopt)
         call rotpolegpu
         call emrecip1gpu
!$acc serial async(rec_queue) present(elambdarec1,emrec,
!$acc& g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz)
         elambdarec1 = emrec
         g_vxx_1 = g_vxx
         g_vxy_1 = g_vxy
         g_vxz_1 = g_vxz
         g_vyy_1 = g_vyy
         g_vyz_1 = g_vyz
         g_vzz_1 = g_vzz
!$acc end serial
         call mem_move(delambdarec1,demrec,siz8,rec_stream)

         elambda = elambdatemp
!$acc wait
!$acc serial async(rec_queue)
!$acc&       present(emrec,delambdae,elambdarec0,elambdarec1,
!$acc& g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz,g_vxx_temp,g_vxy_temp,
!$acc& g_vxz_temp,g_vyy_temp,g_vyz_temp,g_vzz_temp)
         emrec     = (1-elambda)*elambdarec0 + elambda*elambdarec1
         g_vxx = g_vxx_temp + (1.0-elambda)*g_vxx_0+elambda*g_vxx_1
         g_vxy = g_vxy_temp + (1.0-elambda)*g_vxy_0+elambda*g_vxy_1
         g_vxz = g_vxz_temp + (1.0-elambda)*g_vxz_0+elambda*g_vxz_1
         g_vyy = g_vyy_temp + (1.0-elambda)*g_vyy_0+elambda*g_vyy_1
         g_vyz = g_vyz_temp + (1.0-elambda)*g_vyz_0+elambda*g_vyz_1
         g_vzz = g_vzz_temp + (1.0-elambda)*g_vzz_0+elambda*g_vzz_1
         delambdae = delambdae + elambdarec1-elambdarec0
!$acc end serial
!$acc parallel loop async(rec_queue) collapse(2) default(present)
         do i = 1,nlocrec2; do j = 1,3
            demrec(j,i) = (1-elambda)*delambdarec0(j,i)
     &                  +    elambda *delambdarec1(j,i)
         end do; end do

         !reset lambda to initial value
         call altelec(altopt)
         call rotpolegpu
!$acc update host(delambdae) async(rec_queue)

         call timer_exit( timer_rec,quiet_timers )
         end if
      end if
c
#ifdef _OPENACC
      !Finalize async overlapping
      if (dir_queue.ne.rec_queue) call end_dir_stream_cover
#endif
c
c     Add both contribution to the energy
c
!$acc serial async(rec_queue) present(em,emself,em_r,emrec,delambdae)
      em = em + enr2en(emself + em_r) + emrec
!$acc end serial
c
!$acc exit data delete(elambdarec0,elambdarec1
!$acc&    ,delambdarec0,delambdarec1
!$acc&     ,g_vxx_temp,g_vxy_temp,g_vxz_temp
!$acc&     ,g_vyy_temp,g_vyz_temp,g_vzz_temp) async(rec_queue)
      end
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine emreal1c  --  Ewald real space derivs via list  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "emreal1c" evaluates the real space portion of the Ewald
c     summation energy and gradient due to multipole interactions
c     via a neighbor list
c
c
      subroutine emreal1c_
      use atmlst     ,only: poleglobnl
      use atoms      ,only: x,y,z
      use cflux      ,only: dcflux2,adflux2
      use deriv      ,only: dem,de_ws1
      use domdec     ,only: nbloc
      use empole1gpu_inl
      use elec_wspace,only: fix=>r2Work1,fiy=>r2Work2,fiz=>r2Work3
     &               ,tem=>r2Work4
      use interfaces ,only: emreal1ca_p,torquegpu
      use inform
      use mpole      ,only: xaxis,yaxis,zaxis,npolelocnl,ipole
      use potent     ,only: use_mpoleshortreal,use_mpolelong,use_chgflx
      use tinheader  ,only: ti_p,zeror
      use tinmemory  ,only: prmem_request,mipk
      use utils      ,only: set_to_zero1
      use utilcomm
      use utilgpu    ,only: dir_queue,def_queue,rec_queue,mem_set
     &               ,rec_stream,dir_stream,def_stream
     &               ,pot=>ug_workS_r
#ifdef _OPENACC
     &               ,rec_event,stream_wait_async
#endif
      use virial     ,only: vxx=>g_vxx,vxy=>g_vxy,vxz=>g_vxz
     &               ,vyy=>g_vyy,vyz=>g_vyz,vzz=>g_vzz,use_virial
      !use mpi
      use timestat   ,only: timer_enter,timer_exit,timer_emreal
     &               ,quiet_timers
      use tinMemory  ,only: prmem_request
      implicit none
      logical*1,parameter::extract=.true.
      integer   ii,iipole,iglob,iax,iay,iaz
      real(t_p) xix,yix,zix,xiy,yiy,ziy,xiz,yiz,ziz
      character*11 mode

      call timer_enter( timer_emreal )
      def_queue = dir_queue

      call prmem_request(fix,3,npolelocnl)
      call prmem_request(fiy,3,npolelocnl)
      call prmem_request(fiz,3,npolelocnl)
      call prmem_request(tem,3,nbloc)
      if (use_chgflx) call prmem_request(pot,nbloc)

      call set_to_zero1(tem,3*nbloc,def_queue)

      if (use_mpoleshortreal) then
         call switch('SHORTEWALD')
      else
         call switch('EWALD')
      end if

#ifdef _OPENACC
      if (rec_queue.ne.dir_queue) then  !start async overlapping
         def_stream = dir_stream
         call stream_wait_async(rec_stream,dir_stream,rec_event)
      end if
#endif
c
c     calculate multipole potentials
c
      call emreal1ca_p(tem,vxx,vxy,vxz,vyy,vyz,vzz)
c
c     resolve site torques then increment forces and virial
c
      call torquegpu(tem,fix,fiy,fiz,dem,extract)
c
      if (use_virial) call torque_on_vir(fix,fiy,fiz)
c
      if (use_chgflx) then
         call commDDd_add(pot,1,ucComm+ucBig)
         call commDDd_ext(pot,1,ucComm+ucNeig)
         call dcflux2(pot,de_ws1)
         call adflux2(de_ws1,dem)
         call mem_set(pot,zeror,int(nbloc,mipk),def_stream)
      end if
c
      call timer_exit( timer_emreal )

      end

      subroutine torque_on_vir(fix,fiy,fiz)
      use atoms
      use atmlst
      use mpole
      use utilgpu ,only: def_queue
      use virial
      implicit  none
      real(t_p),intent(in),dimension(3,npolelocnl)::fix,fiy,fiz
      integer   ii,iglob,iipole,iax,iay,iaz
      real(t_p) xix,yix,zix,xiy,yiy,ziy,xiz,yiz,ziz

!$acc parallel loop vector_length(32) async(def_queue)
!$acc&         present(poleglobnl,x,y,z,ipole
!$acc&     ,xaxis,yaxis,zaxis,fix,fiy,fiz
!$acc&     ,g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz)
!$acc&         reduction(+:g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz)
      do ii = 1, npolelocnl
         iipole = poleglobnl(ii)
         iglob  = ipole(iipole)
         iaz    = zaxis(iipole)
         iax    = xaxis(iipole)
         iay    = yaxis(iipole)
         if (iaz.le.0) iaz = iglob
         if (iax.le.0) iax = iglob
         if (iay.le.0) iay = iglob
         xix    = x(iax) - x(iglob)
         yix    = y(iax) - y(iglob)
         zix    = z(iax) - z(iglob)
         xiy    = x(iay) - x(iglob)
         yiy    = y(iay) - y(iglob)
         ziy    = z(iay) - z(iglob)
         xiz    = x(iaz) - x(iglob)
         yiz    = y(iaz) - y(iglob)
         ziz    = z(iaz) - z(iglob)
         g_vxx  = g_vxx + xix*fix(1,ii) + xiy*fiy(1,ii) + xiz*fiz(1,ii)
         g_vxy  = g_vxy + yix*fix(1,ii) + yiy*fiy(1,ii) + yiz*fiz(1,ii)
         g_vxz  = g_vxz + zix*fix(1,ii) + ziy*fiy(1,ii) + ziz*fiz(1,ii)
         g_vyy  = g_vyy + yix*fix(2,ii) + yiy*fiy(2,ii) + yiz*fiz(2,ii)
         g_vyz  = g_vyz + zix*fix(2,ii) + ziy*fiy(2,ii) + ziz*fiz(2,ii)
         g_vzz  = g_vzz + zix*fix(3,ii) + ziy*fiy(3,ii) + ziz*fiz(3,ii)
      end do
      end subroutine

      subroutine torque_to_vir(fix,fiy,fiz,vxx,vxy,vxz,vyy,vyz,vzz)
      use atoms
      use atmlst
      use mpole
      use utilgpu ,only: def_queue
      implicit  none
      real(t_p),intent(in),dimension(3,npolelocnl)::fix,fiy,fiz
      real(r_p),intent(inout) :: vxx,vxy,vxz,vyy,vyz,vzz
      integer   ii,iglob,iipole,iax,iay,iaz
      real(t_p) xix,yix,zix,xiy,yiy,ziy,xiz,yiz,ziz

!$acc parallel loop vector_length(32) async(def_queue)
!$acc&         present(poleglobnl,x,y,z,ipole
!$acc&     ,xaxis,yaxis,zaxis,fix,fiy,fiz
!$acc&     ,vxx,vxy,vxz,vyy,vyz,vzz)
!$acc&         reduction(+:vxx,vxy,vxz,vyy,vyz,vzz)
      do ii = 1, npolelocnl
         iipole = poleglobnl(ii)
         iglob  = ipole(iipole)
         iaz    = zaxis(iipole)
         iax    = xaxis(iipole)
         iay    = yaxis(iipole)
         if (iaz.le.0) iaz = iglob
         if (iax.le.0) iax = iglob
         if (iay.le.0) iay = iglob
         xix    = x(iax) - x(iglob)
         yix    = y(iax) - y(iglob)
         zix    = z(iax) - z(iglob)
         xiy    = x(iay) - x(iglob)
         yiy    = y(iay) - y(iglob)
         ziy    = z(iay) - z(iglob)
         xiz    = x(iaz) - x(iglob)
         yiz    = y(iaz) - y(iglob)
         ziz    = z(iaz) - z(iglob)
         vxx    = vxx + xix*fix(1,ii) + xiy*fiy(1,ii) + xiz*fiz(1,ii)
         vxy    = vxy + yix*fix(1,ii) + yiy*fiy(1,ii) + yiz*fiz(1,ii)
         vxz    = vxz + zix*fix(1,ii) + ziy*fiy(1,ii) + ziz*fiz(1,ii)
         vyy    = vyy + yix*fix(2,ii) + yiy*fiy(2,ii) + yiz*fiz(2,ii)
         vyz    = vyz + zix*fix(2,ii) + ziy*fiy(2,ii) + ziz*fiz(2,ii)
         vzz    = vzz + zix*fix(3,ii) + ziy*fiy(3,ii) + ziz*fiz(3,ii)
      end do
      end subroutine

      ! Compute emreal interactions
      subroutine emreal1ca_ac(tem,vxx,vxy,vxz,vyy,vyz,vzz)
      use atmlst ,only:poleglobnl
      use atoms  ,only:x,y,z,n
      use cutoff ,only:ewaldshortcut,shortheal
      use chgpot ,only:electric,dielec
      use deriv  ,only:dem,delambdae
      use domdec ,only:rank,nbloc,loc
      use empole1gpu_inl ,only: image_inl,duo_mpole,groups2_inl
     &                   ,atomic_add,atomic_sub
#ifdef USE_DETERMINISTIC_REDUCTION
     &                   ,tp2enr
#endif
      use energi ,only:em=>em_r
      use ewald  ,only:aewald
      use group  ,only:use_group,ngrp,grplist,wgrp
      use inform ,only:deb_Path
      use math   ,only:sqrtpi
      !use mpi
      use mplpot ,only:n_mscale,mcorrect_ik,mcorrect_scale
      use mpole  ,only:rpole,ipole,polelocnl,npolelocnl
      use mutant ,only:elambda,mutInt
      use neigh  ,only:nelst,elst,shortelst,nshortelst
      use potent ,only:use_mpoleshortreal,use_mpolelong,use_lambdadyn
     &           ,use_chgflx
      use shunt  ,only:off,off2
      use tinheader   ,only:ti_p,zeror
      use tinMemory   ,only:prmem_request
      use utilgpu,only:dir_queue,rec_queue,def_queue,warning
     &                ,maxscaling,pot=>ug_workS_r
      use tinTypes
      use virial
      implicit none

      real(r_p),intent(inout):: vxx,vyy,vzz
      real(r_p),intent(inout):: vxy,vxz,vyz
      real(t_p),intent(inout):: tem(3,nbloc)

      integer i,j,k,iglob,kglob,kbis,ver,ver1,fea
      integer ii,kk,iipole,kkpole
      integer nnelst
      integer,pointer,save::lst(:,:),nlst(:)
#ifdef TINKER_DEBUG
      integer inter(n)
#endif
      integer(1) muti,mutk,mutik
      real(t_p) scut
      real(t_p) r2,f,e,delambdae_,fgrp,poti,potk
      real(t_p) alsq2,alsq2n
      real(t_p) xi,yi,zi
      real(t_p) xr,yr,zr
      type(rpole_elt) ip,kp
      type(mdyn3_r) frc_r
      type(real3) ttmi,ttmk,frc
      real(t_p) mscale
      real(t_p) loff2
      parameter(
     &          ver=__use_grd__+__use_ene__+__use_vir__
     &         ,ver1=ver+__use_sca__)

      if (deb_Path) write(*,'(2x,a,2L3)')
     &   'emreal1c_gpu',use_mpoleshortreal,use_mpolelong

      fea    = __use_mpi__
      if (use_group) fea = fea + __use_groups__
      if (use_lambdadyn) fea = fea + __use_lambdadyn__
      if (use_mpoleshortreal) then
         fea = fea + __use_shortRange__
      else if(use_mpolelong) then
         fea = fea + __use_longRange__
      end if
c
c     set conversion factor, cutoff and switching coefficients
c
      f      = electric / dielec
      alsq2  = 2 * aewald**2
      alsq2n = zeror
      if (aewald .gt. zeror)
     &   alsq2n = 1.0 / (sqrtpi*aewald)

      ! Configure data to be use in next loop
      scut = ewaldshortcut
      if (use_mpoleshortreal) then
         loff2 = 0.0_ti_p
          lst =>  shortelst
         nlst => nshortelst
      else if (use_mpolelong) then
         loff2 = (scut-shortheal)**2
          lst =>  elst
         nlst => nelst
      else
         loff2 = 0.0_ti_p
          lst =>  elst
         nlst => nelst
      end if
c
c     compute the real space portion of the Ewald summation
c     (Scaling interactions)
c
!$acc parallel loop
!$acc&         present(dem,tem,em,delambdae,vxx,vxy,vxz,vyy,vyz,vzz)
!$acc&         present(ipole,loc,x,y,z,rpole,
!$acc&  mcorrect_ik,mcorrect_scale,grplist,wgrp,mutInt)
!$acc&         private(ip,kp,frc,ttmi,ttmk)
!$acc&         reduction(+:em,delambdae,vxx,vxy,vxz,vyy,vyz,vzz)
!$acc&         async(def_queue)
      do ii = 1, n_mscale
         iipole = mcorrect_ik(1,ii)
         kkpole = mcorrect_ik(2,ii)
         mscale = mcorrect_scale(ii)
         iglob  = ipole(iipole)
         kglob  = ipole(kkpole)
         i      = loc(iglob)
         kbis   = loc(kglob)

         xr     = x(kglob) - x(iglob)
         yr     = y(kglob) - y(iglob)
         zr     = z(kglob) - z(iglob)
         call image_inl(xr,yr,zr)
         r2     = xr*xr + yr*yr + zr*zr
         if (r2.lt.loff2.or.r2.gt.off2) cycle  !apply cutoff

         ip%c   = rpole(01,iipole)
         ip%dx  = rpole(02,iipole)
         ip%dy  = rpole(03,iipole)
         ip%dz  = rpole(04,iipole)
         ip%qxx = rpole(05,iipole)
         ip%qxy = rpole(06,iipole)
         ip%qxz = rpole(07,iipole)
         ip%qyy = rpole(09,iipole)
         ip%qyz = rpole(10,iipole)
         ip%qzz = rpole(13,iipole)
         kp%c   = rpole(01,kkpole)
         kp%dx  = rpole(02,kkpole)
         kp%dy  = rpole(03,kkpole)
         kp%dz  = rpole(04,kkpole)
         kp%qxx = rpole(05,kkpole)
         kp%qxy = rpole(06,kkpole)
         kp%qxz = rpole(07,kkpole)
         kp%qyy = rpole(09,kkpole)
         kp%qyz = rpole(10,kkpole)
         kp%qzz = rpole(13,kkpole)
         poti = 0
         potk = 0

         if (use_lambdadyn) mutik=mutInt(iglob)+mutInt(kglob)


         if (use_group)
     &      call groups2_inl(fgrp,iglob,kglob,ngrp,grplist,wgrp)

         ! compute mpole one interaction
         call duo_mpole(r2,xr,yr,zr,ip,kp,mscale,scut
     &           ,shortheal,aewald,f,alsq2n,alsq2,use_group,fgrp
     &           ,use_lambdadyn,mutik,elambda,use_chgflx
     &           ,poti,potk,delambdae_,e,frc,ttmi,ttmk,ver1,fea)

         ! update energy
         em     = em  + tp2enr(e)

         if (use_lambdadyn) delambdae=delambdae+delambdae_

         if (use_chgflx) then
            call atomic_add( pot(i)   ,poti )
            call atomic_add( pot(kbis),potk )
         end if
c
c     increment force-based gradient and torque on first site
c
         call atomic_add(dem(1,i   ),frc%x)
         call atomic_add(dem(2,i   ),frc%y)
         call atomic_add(dem(3,i   ),frc%z)
         call atomic_sub(dem(1,kbis),frc%x)
         call atomic_sub(dem(2,kbis),frc%y)
         call atomic_sub(dem(3,kbis),frc%z)
c
c     increment force-based gradient and torque on second site
c
         call atomic_add(tem(1,i   ),ttmi%x)
         call atomic_add(tem(2,i   ),ttmi%y)
         call atomic_add(tem(3,i   ),ttmi%z)
         call atomic_add(tem(1,kbis),ttmk%x)
         call atomic_add(tem(2,kbis),ttmk%y)
         call atomic_add(tem(3,kbis),ttmk%z)
c
c     increment the virial due to pairwise Cartesian forces c
c
         vxx    = vxx  - xr * frc%x
         vxy    = vxy  - yr * frc%x
         vxz    = vxz  - zr * frc%x
         vyy    = vyy  - yr * frc%y
         vyz    = vyz  - zr * frc%y
         vzz    = vzz  - zr * frc%z
      end do
c
c     compute the real space portion of the Ewald summation
c
!$acc parallel loop gang vector_length(32)
!$acc&         present(dem,tem,em,delambdae,vxx,vxy,vxz,vyy,vyz,vzz)
!$acc&         present(poleglobnl,ipole,loc,x,y,z,rpole,
!$acc&  lst,nlst,grplist,wgrp,mutInt)
!$acc&         private(ip)
!$acc&         reduction(+:em,delambdae,vxx,vxy,vxz,vyy,vyz,vzz)
!$acc&         async(def_queue)
      do ii = 1, npolelocnl
         iipole = poleglobnl(ii)
         iglob  = ipole(iipole)
         i      = loc(iglob)
         nnelst = nlst(ii)
         if (nnelst.eq.0) cycle
         xi     = x(iglob)
         yi     = y(iglob)
         zi     = z(iglob)
         ip%c   = rpole(01,iipole)
         ip%dx  = rpole(02,iipole)
         ip%dy  = rpole(03,iipole)
         ip%dz  = rpole(04,iipole)
         ip%qxx = rpole(05,iipole)
         ip%qxy = rpole(06,iipole)
         ip%qxz = rpole(07,iipole)
         ip%qyy = rpole(09,iipole)
         ip%qyz = rpole(10,iipole)
         ip%qzz = rpole(13,iipole)
         if(use_lambdadyn) muti   = mutInt(iglob)
c
c     evaluate all sites within the cutoff distance
c
!$acc loop vector private(kp,frc,ttmi,ttmk)
!$acc&     reduction(+:em,delambdae,vxx,vxy,vxz,vyy,vyz,vzz)
         do kk = 1, nnelst
            kkpole = lst(kk,ii)
            kglob  = ipole(kkpole)
            kbis   = loc(kglob)
            xr     = x(kglob) - xi
            yr     = y(kglob) - yi
            zr     = z(kglob) - zi
            call image_inl(xr,yr,zr)
            r2     = xr*xr + yr*yr + zr*zr
            if (r2.lt.loff2.or.r2.gt.off2) cycle       !apply cutoff

            kp%c   = rpole( 1,kkpole)
            kp%dx  = rpole( 2,kkpole)
            kp%dy  = rpole( 3,kkpole)
            kp%dz  = rpole( 4,kkpole)
            kp%qxx = rpole( 5,kkpole)
            kp%qxy = rpole( 6,kkpole)
            kp%qxz = rpole( 7,kkpole)
            kp%qyy = rpole( 9,kkpole)
            kp%qyz = rpole(10,kkpole)
            kp%qzz = rpole(13,kkpole)
            poti   = 0
            potk   = 0

            if (use_lambdadyn) mutik=muti+mutInt(kglob)

            if (use_group)
     &         call groups2_inl(fgrp,iglob,kglob,ngrp,grplist,wgrp)

            ! compute mpole one interaction
            call duo_mpole(r2,xr,yr,zr,ip,kp,zeror
     &              ,scut,shortheal,aewald,f,alsq2n,alsq2,use_group,fgrp
     &              ,use_lambdadyn,mutik,elambda,use_chgflx
     &              ,poti,potk,delambdae_,e,frc,ttmi,ttmk,ver,fea)

            ! update energy
            em     = em  + tp2enr(e)

            ! update hamiltonian derivative
            if (use_lambdadyn) delambdae=delambdae+delambdae_
            if (use_chgflx) then
               call atomic_add( pot(i)   ,poti )
               call atomic_add( pot(kbis),potk )
            end if
c
c     increment force-based gradient and torque on first site
c
            call atomic_add(dem(1,i   ),frc%x)
            call atomic_add(dem(2,i   ),frc%y)
            call atomic_add(dem(3,i   ),frc%z)
            call atomic_sub(dem(1,kbis),frc%x)
            call atomic_sub(dem(2,kbis),frc%y)
            call atomic_sub(dem(3,kbis),frc%z)
c
c     increment force-based gradient and torque on second site
c
            call atomic_add(tem(1,i   ),ttmi%x)
            call atomic_add(tem(2,i   ),ttmi%y)
            call atomic_add(tem(3,i   ),ttmi%z)
            call atomic_add(tem(1,kbis),ttmk%x)
            call atomic_add(tem(2,kbis),ttmk%y)
            call atomic_add(tem(3,kbis),ttmk%z)
c
c     increment the virial due to pairwise Cartesian forces
c
            vxx    = vxx  - xr * frc%x
            vxy    = vxy  - yr * frc%x
            vxz    = vxz  - zr * frc%x
            vyy    = vyy  - yr * frc%y
            vyz    = vyz  - zr * frc%y
            vzz    = vzz  - zr * frc%z
         end do
      end do

      end

      ! CUDA Fortran routine
#ifdef _CUDA
      subroutine emreal1ca_cu(tem,vxx,vxy,vxz,vyy,vyz,vzz)
      use atmlst ,only:poleglobnl
      use atoms  ,only:x,y,z,n
      use cutoff ,only:mpoleshortcut,shortheal,ewaldshortcut
      use chgpen ,only:pcore,pval,palpha
      use chgpot ,only:electric,dielec
      use cell
      use deriv  ,only:dem,delambdae
      use domdec ,only: xbegproc,ybegproc,zbegproc
     &           ,nproc,rank,xendproc,yendproc,zendproc
     &           ,nbloc,loc
      use empole1cu
      use empole_cpencu
      use energi ,only:em=>em_r
      use ewald  ,only:aewald
      use group  ,only:use_group,grplist,wgrp
      use inform ,only: deb_Path,minmaxone
      use math   ,only:sqrtpi
      !use mpi
      use mplpot ,only:n_mscale,mcorrect_ik,mcorrect_scale,pentyp_i
      use mpole  ,only:rpole,ipole,polelocnl,npolelocnl
     &           ,npolelocnlb,npolelocnlb_pair,npolebloc
     &           ,npolelocnlb2_pair,nshortpolelocnlb2_pair,ipole
      use mutant ,only:mutInt,elambda
      use neigh  ,only:ipole_s=>celle_glob,pglob_s=>celle_pole
     &           , loc_s=>celle_loc, ieblst_s=>ieblst, eblst_s=>eblst
     &           , x_s=>celle_x, y_s=>celle_y, z_s=>celle_z
     &           , seblst_s=>shorteblst,iseblst_s=>ishorteblst
      use polpot ,only:use_thole
      use potent ,only:use_mpoleshortreal,use_mpolelong,use_lambdadyn
     &           ,use_chgpen,use_chgflx
      use shunt  ,only:off,off2
      use tinheader,only: ti_p
      use tinMemory,only: prmem_request,mipk
      use utilcu ,only:BLOCK_DIM,check_launch_kernel
      use utilgpu,only:dir_queue,rec_queue,def_queue,warning
     &           ,RED_BUFF_SIZE,BLOCK_SIZE
     &           ,pot=>ug_workS_r
     &           ,ered_buff=>ered_buf1,vred_buff,lam_buff,nred_buff
     &           ,reduce_energy_virial,reduce_buffer,get_GridDim
     &           ,zero_evir_red_buffer
     &           ,dir_stream,def_stream
      !use tinTypes
      use virial
      implicit none

      real(r_p),intent(inout):: vxx,vyy,vzz
      real(r_p),intent(inout):: vxy,vxz,vyz
      real(t_p),intent(inout):: tem(3,nbloc)

      integer i,gS1
      integer start,start1,sized
      real(t_p) r2,f,e
      real(t_p) alsq2,alsq2n
      real(t_p) p_xbeg,p_xend,p_ybeg,p_yend,p_zbeg,p_zend
      real(t_p) loff2,lcut
      integer,save::ndec
      logical,save:: first_in=.true.
      integer,save:: gS
      real(8) e0,e1,e2
      character*64 rtami

      if (deb_Path) then
         rtami = ' emreal1d_cu'//merge(' CHGPEN','',use_chgpen)
     &          //merge(' CHGFLX','',use_chgflx)
         if (use_mpoleshortreal) then
            rtami = trim(rtami)//' SHORT'
         else if (use_mpolelong) then
            rtami = trim(rtami)//' LONG'
         end if
         write(*,'(a)') trim(rtami)
      end if
c
c     set conversion factor, cutoff and switching coefficients
c
      f      = electric / dielec
      alsq2  = 2 * aewald**2
      alsq2n = 0
      if (aewald .gt. 0.0)
     &   alsq2n = 1.0 / (sqrtpi*aewald)

      p_xbeg = xbegproc(rank+1)
      p_xend = xendproc(rank+1)
      p_ybeg = ybegproc(rank+1)
      p_yend = yendproc(rank+1)
      p_zbeg = zbegproc(rank+1)
      p_zend = zendproc(rank+1)

      if (use_mpoleshortreal) then
         loff2 = 0.0_ti_p
      else if (use_mpolelong) then
         loff2 = (ewaldshortcut-shortheal)**2
      else
         loff2 = 0.0_ti_p
      end if
      lcut  = ewaldshortcut

      if (first_in) then
         call cudaMaxGridSize("emreal1_kcu",gS)
         ndec=1
         !if (dir_queue.ne.rec_queue) ndec=1
         first_in = .false.
      end if

      def_stream = dir_stream
      def_queue  = dir_queue

      if      (use_chgpen) then
!$acc host_data use_device(ipole_s,pglob_s,loc_s
!$acc&    ,ieblst_s,eblst_s,iseblst_s,seblst_s
!$acc&    ,x_s,y_s,z_s,pcore,pval,palpha,rpole,dem,tem,pot
!$acc&    ,grplist,wgrp,mutInt
!$acc&    ,ered_buff,vred_buff,lam_buff,nred_buff
!$acc&    ,mcorrect_ik,mcorrect_scale,ipole,loc,x,y,z
!$acc&    )

      start1 = 2*npolelocnlb_pair+1
         if      (use_mpoleshortreal) then
      call emreal_cpen1s_kcu<<<gS,BLOCK_DIM,0,def_stream>>>
     &    (ipole_s,pglob_s,loc_s
     &    ,iseblst_s,seblst_s(start1)
     &    ,npolelocnlb,npolelocnlb2_pair,npolebloc,n,pentyp_i
     &    ,x_s,y_s,z_s,pcore,pval,palpha,rpole
     &    ,shortheal,lcut,loff2,off2,f,aewald
     &    ,use_group,grplist,wgrp,use_lambdadyn,elambda,mutInt
     &    ,use_chgflx
     &    ,dem,tem,pot,ered_buff,vred_buff,lam_buff,nred_buff
     &    ,p_xbeg,p_xend,p_ybeg,p_yend,p_zbeg,p_zend
     &    ,mcorrect_ik,mcorrect_scale,ipole,loc,x,y,z,n_mscale
     &    )
      call check_launch_kernel(" emreal_cpen1s_kcu")
         else if (use_mpolelong) then
      call emreal_cpen1l_kcu<<<gS,BLOCK_DIM,0,def_stream>>>
     &    (ipole_s,pglob_s,loc_s
     &    ,ieblst_s,eblst_s(start1)
     &    ,npolelocnlb,npolelocnlb2_pair,npolebloc,n,pentyp_i
     &    ,x_s,y_s,z_s,pcore,pval,palpha,rpole
     &    ,shortheal,lcut,loff2,off2,f,aewald
     &    ,use_group,grplist,wgrp,use_lambdadyn,elambda,mutInt
     &    ,use_chgflx
     &    ,dem,tem,pot,ered_buff,vred_buff,lam_buff,nred_buff
     &    ,p_xbeg,p_xend,p_ybeg,p_yend,p_zbeg,p_zend
     &    ,mcorrect_ik,mcorrect_scale,ipole,loc,x,y,z,n_mscale
     &    )
      call check_launch_kernel(" emreal_cpen1l_kcu")
         else
      call emreal_cpen1_kcu<<<gS,BLOCK_DIM,0,def_stream>>>
     &    (ipole_s,pglob_s,loc_s
     &    ,ieblst_s,eblst_s(start1)
     &    ,npolelocnlb,npolelocnlb2_pair,npolebloc,n,pentyp_i
     &    ,x_s,y_s,z_s,pcore,pval,palpha,rpole
     &    ,shortheal,lcut,loff2,off2,f,aewald
     &    ,use_group,grplist,wgrp,use_lambdadyn,elambda,mutInt
     &    ,use_chgflx
     &    ,dem,tem,pot,ered_buff,vred_buff,lam_buff,nred_buff
     &    ,p_xbeg,p_xend,p_ybeg,p_yend,p_zbeg,p_zend
     &    ,mcorrect_ik,mcorrect_scale,ipole,loc,x,y,z,n_mscale
     &    )
      call check_launch_kernel(" emreal_cpen1_kcu")
         end if
!$acc end host_data
      else !(use_thole)

!$acc host_data use_device(ipole_s,pglob_s,loc_s
!$acc&    ,ieblst_s,eblst_s,iseblst_s,seblst_s
!$acc&    ,x_s,y_s,z_s,rpole,dem,tem,pot
!$acc&    ,grplist,wgrp,mutInt
!$acc&    ,ered_buff,vred_buff,lam_buff,nred_buff
!$acc&    ,mcorrect_ik,mcorrect_scale,ipole,loc,x,y,z
!$acc&    )

      !call CUDA kernel to compute the real space portion of the Ewald summation
         if (use_mpoleshortreal) then
         start1 = 2*npolelocnlb_pair+1
      call emreal1s_kcu<<<gS,BLOCK_DIM,0,def_stream>>>
     &    (ipole_s,pglob_s,loc_s,iseblst_s,seblst_s(start1)
     &    ,npolelocnlb,nshortpolelocnlb2_pair,npolebloc,n
     &    ,x_s,y_s,z_s,rpole
     &    ,shortheal,lcut,loff2,off2,f,alsq2,alsq2n,aewald
     &    ,use_group,grplist,wgrp,use_lambdadyn,elambda,mutInt
     &    ,use_chgflx,pot
     &    ,dem,tem,ered_buff,vred_buff,lam_buff,nred_buff
     &    ,p_xbeg,p_xend,p_ybeg,p_yend,p_zbeg,p_zend
     &    ,mcorrect_ik,mcorrect_scale,ipole,loc,x,y,z,n_mscale
     &    )
      call check_launch_kernel(" empole1s_kcu")

         else if (use_mpolelong) then
      ! Split long range electrostatic kernel to ease recovering process in MPI
            sized = npolelocnlb2_pair/ndec
            do i = 1,ndec
         start  = (i-1)*sized + 1
         start1 = 2*npolelocnlb_pair+1+(start-1)*BLOCK_SIZE
         if (i.eq.ndec) sized = npolelocnlb2_pair-start+1
         call emreal1l_kcu<<<gS,BLOCK_DIM,0,def_stream>>>
     &       (ipole_s,pglob_s,loc_s
     &       ,ieblst_s(start),eblst_s(start1)
     &       ,npolelocnlb,sized,npolebloc,n
     &       ,x_s,y_s,z_s,rpole
     &       ,shortheal,lcut,loff2,off2,f,alsq2,alsq2n,aewald
     &       ,use_group,grplist,wgrp,use_lambdadyn,elambda,mutInt
     &       ,use_chgflx,pot
     &       ,dem,tem,ered_buff,vred_buff,lam_buff,nred_buff
     &       ,p_xbeg,p_xend,p_ybeg,p_yend,p_zbeg,p_zend
     &       ,mcorrect_ik,mcorrect_scale,ipole,loc,x,y,z,n_mscale
     &       )
            end do
      call check_launch_kernel(" empole1l_kcu")

         else
            sized = npolelocnlb2_pair/ndec
            do i = 1,ndec
         start  = (i-1)*sized + 1
         start1 = 2*npolelocnlb_pair+1+(start-1)*BLOCK_SIZE
         if (i.eq.ndec) sized = npolelocnlb2_pair-start+1
         call emreal1_kcu<<<gS,BLOCK_DIM,0,def_stream>>>
     &       (ipole_s,pglob_s,loc_s
     &       ,ieblst_s(start),eblst_s(start1)
     &       ,npolelocnlb,sized,npolebloc,n
     &       ,x_s,y_s,z_s,rpole
     &       ,shortheal,lcut,loff2,off2,f,alsq2,alsq2n,aewald
     &       ,use_group,grplist,wgrp,use_lambdadyn,elambda,mutInt
     &       ,use_chgflx,pot
     &       ,dem,tem,ered_buff,vred_buff,lam_buff,nred_buff
     &       ,p_xbeg,p_xend,p_ybeg,p_yend,p_zbeg,p_zend
     &       ,mcorrect_ik,mcorrect_scale,ipole,loc,x,y,z,n_mscale
     &       )
            end do
      call check_launch_kernel(" empole1_kcu")

         end if  !(short/long range)

!$acc end host_data

      end if !(use_[chgpen/thole])

      call reduce_energy_virial(em,vxx,vxy,vxz,vyy,vyz,vzz
     &                         ,ered_buff,def_queue)

      if (use_lambdadyn)
     &   call reduce_buffer(lam_buff,RED_BUFF_SIZE,delambdae,def_queue)

      end

      ! CUDA C wrapper on emreal1c (check file cu_mpole1.cu)
      subroutine emreal1c_core4(tem,vxx,vxy,vxz,vyy,vyz,vzz)
      use atmlst ,only:poleglobnl
      use atoms  ,only:x,y,z,n
      use chgpot ,only:electric,dielec
      use cell
      use deriv  ,only:dem
      use domdec ,only:xbegproc,ybegproc,zbegproc
     &           ,nproc,rank,xendproc,yendproc,zendproc
     &           ,nbloc,loc
      use energi ,only:em=>em_r
      !use erf_mod
      use inform
      use ewald  ,only:aewald
      use group 
      use inform ,only:deb_Path
      use interfaces  ,only:cu_emreal1c
      use math   ,only:sqrtpi
      !use mpi
      use mplpot ,only:n_mscale,mcorrect_ik,mcorrect_scale
      use mpole  ,only:rpole,ipole,polelocnl,npolelocnl
     &           ,npolelocnlb,npolelocnlb_pair,npolebloc
     &           ,npolelocnlb2_pair
      use neigh  ,only:ipole_s=>celle_glob,pglob_s=>celle_pole
     &           ,loc_s=>celle_loc, ieblst_s=>ieblst, eblst_s=>eblst
     &           ,x_s=>celle_x, y_s=>celle_y, z_s=>celle_z
      use shunt  ,only:off2
      use utilgpu,only:dir_queue,rec_queue,def_queue,warning
     &           ,rpole_elt,RED_BUFF_SIZE
     &           ,ered_buff=>ered_buf1,vred_buff,reduce_energy_virial
     &           ,zero_evir_red_buffer
     &           ,dir_stream,def_stream
      use virial
      implicit none

      real(r_p),intent(inout):: vxx,vyy,vzz
      real(r_p),intent(inout):: vxy,vxz,vyz
      real(t_p),intent(inout):: tem(3,nbloc)

      integer i
      real(t_p) r2,f,e
      real(t_p) alsq2,alsq2n
      real(t_p) p_xbeg,p_xend,p_ybeg,p_yend,p_zbeg,p_zend

      if(deb_Path) write(*,'(2x,a)') 'emreal1c_core4'

#ifdef USE_DETERMINISTIC_REDUCTION
 13   format(3x,"ERROR : empole CUDA-C Routine is not to be with"
     &      ," fixed precision !!! ")
      write (0,13)
      __TINKER_FATAL__
#endif
 14   format(' Pairwise scaling computation is missing !!')
      write(0,14) 
      __TINKER_FATAL__
c
c     set conversion factor, cutoff and switching coefficients
c
      f      = electric / dielec
      alsq2  = 2 * aewald**2
      alsq2n = 0
      if (aewald .gt. 0.0)
     &   alsq2n = 1.0 / (sqrtpi*aewald)

      p_xbeg = xbegproc(rank+1)
      p_xend = xendproc(rank+1)
      p_ybeg = ybegproc(rank+1)
      p_yend = yendproc(rank+1)
      p_zbeg = zbegproc(rank+1)
      p_zend = zendproc(rank+1)

c
c     call C wrapper to compute the real space portion of the Ewald summation
c
      def_stream = dir_stream
      call zero_evir_red_buffer(def_queue)

!$acc host_data use_device(ipole_s,pglob_s,loc_s,ieblst_s,eblst_s,
!$acc&    x_s,y_s,z_s,rpole,wgrp,grplist,dem,tem,ered_buff,vred_buff)

      call cu_emreal1c( ipole_s,pglob_s,loc_s,ieblst_s
     &                , eblst_s(2*npolelocnlb_pair+1)
     &                , x_s,y_s,z_s,rpole
     &                , dem,tem,ered_buff,vred_buff
     &                , npolelocnlb,npolelocnlb2_pair,npolebloc,n
     &                , off2,f,alsq2,alsq2n,aewald
     &                , xcell,ycell,zcell,xcell2,ycell2,zcell2
     &                , p_xbeg,p_xend,p_ybeg,p_yend,p_zbeg,p_zend
     &                , def_stream 
     &                , ngrp,use_group,wgrp,grplist)

!$acc end host_data

      call reduce_energy_virial(em,vxx,vxy,vxz,vyy,vyz,vzz
     &                         ,ered_buff,def_queue)
      end
#endif
c
c
c     #######################################################################
c     ##                                                                   ##
c     ##  subroutine emrecip1gpu  --  PME recip multipole energy & derivs  ##
c     ##                                                                   ##
c     #######################################################################
c
c
c     "emrecip1" evaluates the reciprocal space portion of the particle
c     mesh Ewald summation energy and gradient due to multipoles
c
      module emrecip_mod
      integer :: emr_ncall=0
      real(r_p) e
      real(r_p) vxx,vyy,vzz
      real(r_p) vxy,vxz,vyz
      !indices into the electrostatic field array
      integer deriv1(10)
      integer deriv2(10)
      integer deriv3(10)
      data  deriv1  / 2, 5,  8,  9, 11, 16, 18, 14, 15, 20 /
      data  deriv2  / 3, 8,  6, 10, 14, 12, 19, 16, 20, 17 /
      data  deriv3  / 4, 9, 10,  7, 15, 17, 13, 20, 18, 19 /
c     parameter(
c    &  deriv1=(/ 2, 5,  8,  9, 11, 16, 18, 14, 15, 20 /),
c    &  deriv2=(/ 3, 8,  6, 10, 14, 12, 19, 16, 20, 17 /),
c    &  deriv3=(/ 4, 9, 10,  7, 15, 17, 13, 20, 18, 19 /))
c
      end module

      subroutine emrecip1gpu
      use atmlst
      use atoms
      use bound
      use boxes
      use chgpot
      use cflux
      use deriv
      use domdec
      use emrecip_mod
      use energi
      use ewald
      use fft
      use inform    ,only: deb_Path,minmaxone
      use interfaces,only: torquegpu,fphi_mpole_site_p
     &              ,grid_mpole_site_p
     &              ,emreal1c_cp
      use mutant
      use polar_temp,only: cmp=>fphid,fmp=>fphip !Use register pool
     &              , trqrec=>fuind
      use pme
      use pme1
      use math
      use mpole
      use potent
      use timestat
      use tinheader ,only: zeror
      use virial
      use utilcomm
      use utilgpu   ,pot=>ug_workS_r,potrec=>ug_workS_r1
      use utils     ,only:set_to_zero1,rpole_scale
     &              ,rpole_ind_extract,comput_norm
      use mpi
      implicit none
      integer status(MPI_STATUS_SIZE),tag,ierr
      integer i,j,k,ii,iipole,iglob,iloc
      integer k1,k2,k3
      integer m1,m2,m3
      integer ntot,nff
      integer nf1,nf2,nf3
      integer gridin_size,gridout_size
      real(t_p) zero,one,two,half
      real(t_p) eterm,f
      real(t_p) r1,r2,r3
      real(t_p) h1,h2,h3
      real(t_p) f1,f2,f3
      real(t_p) volterm,denom
      real(t_p) hsq,expterm
      real(t_p) term,pterm
      real(t_p) vterm,struc2
      real(t_p) trq(3),fix(3)
      real(t_p) fiy(3),fiz(3)
      integer reqsend(nproc),reqrec(nproc)
      integer req2send(nproc),req2rec(nproc)
      integer nprocloc,commloc,rankloc,proc
      real(8) time0,time1
      parameter(zero=0.0_ti_p,  one=1.0_ti_p,
     &           two=2.0_ti_p, half=0.5_ti_p)
c
c     return if the Ewald coefficient is zero
c
      if (aewald .lt. 1.0d-6)  return
c
      if (aewald .lt. 1.0d-6)  return
      f = electric / dielec

      if (deb_Path) write(*,'(2x,a)') 'emrecip1gpu'
      call timer_enter( timer_emrecip )

      emr_ncall = emr_ncall+1
      if (emr_ncall.eq.1) then
!$acc enter data create(e,vxx,vxy,vxz,vyy,vyz,vzz) 
!$acc&           copyin(deriv1,deriv2,deriv3)
      end if

      if (use_pmecore) then
        nprocloc = nrec
        commloc  = comm_rec
        rankloc  = rank_bis
      else
        nprocloc = nproc
        commloc  = COMM_TINKER
        rankloc  = rank
      end if
      f = electric / dielec
c
c     perform dynamic allocation of some global arrays
c
      j = max(npolerecloc,1)
      call mallocMpiGrid
      call prmem_request(cphirec,10,j,async=.false.)
      call prmem_request(fphirec,20,j,async=.false.)
      call prmem_request(cmp    ,10,j,async=.false.)
      call prmem_request(fmp    ,10,j,async=.false.)
      call prmem_request(trqrec , 3,j,async=.false.)
c
c     zero out the temporary virial accumulation variables
c
!$acc serial async present(e,vxx,vxy,vxz,vyy,vyz,vzz)
      vxx = zero
      vxy = zero
      vxz = zero
      vyy = zero
      vyz = zero
      vzz = zero
      e   = zero
!$acc end serial
c
      call timer_enter( timer_other )
!$acc parallel loop collapse(2) async(rec_queue) default(present)
      do i=1,npolerecloc; do j=1,20
         if (j.lt.10) cphirec(j,i) = zero
         fphirec(j,i) = zero
      end do; end do
      call timer_exit ( timer_other,quiet_timers )
c
c     zero out pme grid
c
      call timer_enter( timer_grid1 )
      call set_to_zero1(qgridin_2d(1,1,1,1,1),
     &     2*n1mpimax*n2mpimax*n3mpimax*(nrec_send+1),rec_queue)
c
c     copy multipole moments and coordinates to local storage
c
!$acc parallel loop collapse(2) async(rec_queue) default(present)
      do i = 1, npolerecloc
         do j = 1, 10
            iipole   = polerecglob(i)
            cmp(j,i) = rpole_scale(j)*rpole(rpole_ind_extract(j),iipole)
         end do
      end do

      call cmp_to_fmp_sitegpu(cmp,fmp)
c
c     assign permanent multipoles to PME grid and perform
c     the 3-D FFT forward transformation
c
      call grid_mpole_site_p(fmp)
      call timer_exit( timer_grid1,quiet_timers )

      !Exchange Grid Between reciprocal process
      call commGridFront( qgridin_2d,r_comm )

#ifdef _OPENACC
      ! Recover MPI communication with real space computations
      if (dir_queue.ne.rec_queue) then
         call start_dir_stream_cover
         call emreal1c_cp
      end if
#endif

      ! Wait for Grid communication
      call commGridFront( qgridin_2d,r_wait )
c
c     Perform 3-D FFT forward transform
c
#ifdef _OPENACC
      call cufft2d_frontmpi(qgridin_2d,qgridout_2d,n1mpimax,n2mpimax,
     $                      n3mpimax)
#else
      call   fft2d_frontmpi(qgridin_2d,qgridout_2d,n1mpimax,n2mpimax,
     $                      n3mpimax)
#endif
c
c     initialize variables required for the scalar summation
c
      pterm   = (pi/aewald)**2
      volterm = pi * volbox
      nff     = nfft1 * nfft2
      nf1     = (nfft1+1) / 2
      nf2     = (nfft2+1) / 2
      nf3     = (nfft3+1) / 2
c
      call timer_enter( timer_scalar )
!$acc serial async(rec_queue) default(present)
      if ((istart2(rankloc+1).eq.1).and.
     &    (jstart2(rankloc+1).eq.1).and.
     &    (kstart2(rankloc+1).eq.1)) then
         qfac_2d(1,1,1) = zero
      end if
!$acc end serial
c
c     make the scalar summation over reciprocal lattice
c
!$acc parallel loop collapse(3) default(present)
!!$acc&        reduction(+:vxx,vxy,vxz,vyy,vyz,vzz)
!$acc&         present(vxx,vxy,vxz,vyy,vyz,vzz,use_bounds)
!$acc&         async(rec_queue)
      do k3 = kstart2(rankloc+1),kend2(rankloc+1)
        do k2 = jstart2(rankloc+1),jend2(rankloc+1)
          do k1 = istart2(rankloc+1),iend2(rankloc+1)
            m1 = k1 - 1
            m2 = k2 - 1
            m3 = k3 - 1
            if (k1 .gt. nf1)  m1 = m1 - nfft1
            if (k2 .gt. nf2)  m2 = m2 - nfft2
            if (k3 .gt. nf3)  m3 = m3 - nfft3
            if ((m1.eq.0).and.(m2.eq.0).and.(m3.eq.0)) goto 10
            r1 = real(m1,t_p)
            r2 = real(m2,t_p)
            r3 = real(m3,t_p)
            h1 = recip(1,1)*r1 + recip(1,2)*r2 + recip(1,3)*r3
            h2 = recip(2,1)*r1 + recip(2,2)*r2 + recip(2,3)*r3
            h3 = recip(3,1)*r1 + recip(3,2)*r2 + recip(3,3)*r3
            hsq = h1*h1 + h2*h2 + h3*h3
            term = -pterm * hsq
            expterm = zero
            if (term .gt. -50.0_ti_p) then
               denom = volterm*hsq*bsmod1(k1)*bsmod2(k2)*bsmod3(k3)
               expterm = exp(term) / denom
               if (.not. use_bounds) then
                  expterm = expterm * (1.0_ti_p-cos(pi*xbox*sqrt(hsq)))
               else if (octahedron) then
                  if (mod(m1+m2+m3,2) .ne. 0)  expterm = zero
               end if
               struc2 = qgridout_2d(1,k1-istart2(rankloc+1)+1,
     &          k2-jstart2(rankloc+1)+1,k3-kstart2(rankloc+1)+1)**2
     &                + qgridout_2d(2,k1-istart2(rankloc+1)+1,
     &          k2-jstart2(rankloc+1)+1,k3-kstart2(rankloc+1)+1)**2
               eterm = half * electric * expterm * struc2
               vterm = (two/hsq) * (1.0_ti_p-term) * eterm
               vxx = vxx + h1*h1*vterm - eterm
               vxy = vxy + h1*h2*vterm
               vxz = vxz + h1*h3*vterm
               vyy = vyy + h2*h2*vterm - eterm
               vyz = vyz + h2*h3*vterm
               vzz = vzz + h3*h3*vterm - eterm
             qfac_2d(k1-istart2(rankloc+1)+1,k2-jstart2(rankloc+1)+1,k3-
     $       kstart2(rankloc+1)+1) = expterm
            end if
 10         continue
          end do
        end do
      end do
c
c     save the virial for use in polarization computation
c
!$acc serial async(rec_queue)
!$acc&       present(e,emrec,vxx,vxy,vxz,vyy,vyz,vzz)
      vmxx = vxx
      vmxy = vxy
      vmxz = vxz
      vmyy = vyy
      vmyz = vyz
      vmzz = vzz
c
c     account for zeroth grid point for nonperiodic system
c
      if ((istart2(rankloc+1).eq.1).and.
     &    (jstart2(rankloc+1).eq.1).and.
     &    (kstart2(rankloc+1).eq.1)) then
        if (.not. use_bounds) then
           expterm = half * pi / xbox
           struc2  = qgrid2in_2d(1,1,1,1,1)**2 +
     &               qgrid2in_2d(2,1,1,1,1)**2
           e       = f * expterm * struc2
           emrec   = emrec + e
        end if
      end if
!$acc end serial
!!$acc update host(vmxx,vmxy,vmxz,vmyy,vmyz,vmzz) async(rec_queue)
c
c     complete the transformation of the PME grid
c
!$acc parallel loop collapse(3) async(rec_queue) default(present)
      do k = 1, ksize2(rankloc+1)
         do j = 1, jsize2(rankloc+1)
            do i = 1, isize2(rankloc+1)
              term = qfac_2d(i,j,k)
              qgridout_2d(1,i,j,k) = term*qgridout_2d(1,i,j,k)
              qgridout_2d(2,i,j,k) = term*qgridout_2d(2,i,j,k)
            end do
         end do
      end do
      call timer_exit( timer_scalar,quiet_timers )

c     nff   = 2*isize2(rankloc+1)*jsize2(rankloc+1)*ksize2(rankloc+1)
c     time0 = comput_norm(qgridout_2d,nff,0)
c     time1 = comput_norm(qgridout_2d,nff,1)
c     call MPI_ALLREDUCE(MPI_in_place,time0,1,MPI_real8,MPI_MAX,
c    &     COMM_TINKER,i)
c     call MPI_ALLREDUCE(MPI_in_place,time1,1,MPI_real8,MPI_SUM,
c    &     COMM_TINKER,i)
c     if (rank.eq.0) print*,' gridout Li ',time0
c     if (rank.eq.0) print*,' gridout L1 ',time1
c
c     perform 3-D FFT backward transform and get potential
c
#ifdef _OPENACC
      call cufft2d_backmpi(qgridin_2d,qgridout_2d,n1mpimax,n2mpimax,
     $                     n3mpimax)
#else
      call   fft2d_backmpi(qgridin_2d,qgridout_2d,n1mpimax,n2mpimax,
     $                     n3mpimax)
#endif

      !Exchange Grid Between reciprocal process
      call commGridBack( qgridin_2d,r_comm )

      ! Wait for Grid communication
      call commGridBack( qgridin_2d,r_wait )

#ifdef _OPENACC
      ! sync streams
      if (dir_queue.ne.rec_queue) call end_dir_stream_cover
#endif

      call timer_enter ( timer_grid2 )
      call fphi_mpole_site_p
!$acc parallel loop collapse(2) async(rec_queue) default(present)
      do i = 1,npolerecloc; do j = 1,20;
         fphirec(j,i) = electric * fphirec(j,i)
      end do; end do;
      call fphi_to_cphi_sitegpu(fphirec,cphirec)
      call timer_exit( timer_grid2,quiet_timers )
c
c     increment the permanent multipole energy and gradient
c
      call timer_enter( timer_fmanage )
!$acc parallel loop async(rec_queue) reduction(+:e) default(present)
!$acc&         present(e)
      do i = 1, npolerecloc
         f1 = zero
         f2 = zero
         f3 = zero
!$acc loop seq
         do k = 1, 10
            e  = e  + fmp(k,i)*fphirec(       k ,i)
            f1 = f1 + fmp(k,i)*fphirec(deriv1(k),i)
            f2 = f2 + fmp(k,i)*fphirec(deriv2(k),i)
            f3 = f3 + fmp(k,i)*fphirec(deriv3(k),i)
         end do
         f1 = real(nfft1,t_p) * f1
         f2 = real(nfft2,t_p) * f2
         f3 = real(nfft3,t_p) * f3
         h1 = recip(1,1)*f1 + recip(1,2)*f2 + recip(1,3)*f3
         h2 = recip(2,1)*f1 + recip(2,2)*f2 + recip(2,3)*f3
         h3 = recip(3,1)*f1 + recip(3,2)*f2 + recip(3,3)*f3
         iipole = polerecglob(i)
         iglob  = ipole(iipole)
         ii     = locrec(iglob)
         demrec(1,ii) = demrec(1,ii) + h1
         demrec(2,ii) = demrec(2,ii) + h2
         demrec(3,ii) = demrec(3,ii) + h3
      end do
c
c     distribute torques into the permanent multipole gradient
c
!$acc parallel loop async(rec_queue) default(present)
      do i = 1, npolerecloc
         !iipole = polerecglob(i)
         trqrec(1,i) = cmp(4,i)*cphirec(3,i) - cmp(3,i)*cphirec(4,i)
     &               + two*(cmp(7,i)-cmp(6,i))*cphirec(10,i)
     &               + cmp(9,i)*cphirec(8,i) + cmp(10,i)*cphirec(6,i)
     &               - cmp(8,i)*cphirec(9,i) - cmp(10,i)*cphirec(7,i)
         trqrec(2,i) = cmp(2,i)*cphirec(4,i) - cmp(4,i)*cphirec(2,i)
     &               + two*(cmp(5,i)-cmp(7,i))*cphirec(9,i)
     &               + cmp(8,i)*cphirec(10,i) + cmp(9,i)*cphirec(7,i)
     &               - cmp(9,i)*cphirec(5,i) - cmp(10,i)*cphirec(8,i)
         trqrec(3,i) = cmp(3,i)*cphirec(2,i) - cmp(2,i)*cphirec(3,i)
     &               + two*(cmp(6,i)-cmp(5,i))*cphirec(8,i)
     &               + cmp(8,i)*cphirec(5,i) + cmp(10,i)*cphirec(9,i)
     &               - cmp(8,i)*cphirec(6,i) - cmp(9,i)*cphirec(10,i)
      end do

      if (use_virial) then   ! COMPUTE virial ?
c
c     permanent multipole contribution to the internal virial
c
!$acc parallel loop async(rec_queue) default(present)
!$acc&         present(vxx,vxy,vxz,vyy,vyz,vzz)
      do i = 1, npolerecloc
         vxx =vxx - cmp(2,i)*cphirec(2,i) - two*cmp(5,i)*cphirec(5,i)
     &            - cmp(8,i)*cphirec(8,i) - cmp(9,i)*cphirec(9,i)
         vxy =vxy - half*(cmp(3,i)*cphirec(2,i)+cmp(2,i)*cphirec(3,i))
     &            -(cmp(5,i)+cmp(6,i))*cphirec(8,i)
     &            -half*cmp(8,i)*(cphirec(5,i)+cphirec(6,i))
     &            -half*(cmp(9,i)*cphirec(10,i)+cmp(10,i)*cphirec(9,i))
         vxz =vxz - half*(cmp(4,i)*cphirec(2,i)+cmp(2,i)*cphirec(4,i))
     &            -(cmp(5,i)+cmp(7,i))*cphirec(9,i)
     &            -half*cmp(9,i)*(cphirec(5,i)+cphirec(7,i))
     &            -half*(cmp(8,i)*cphirec(10,i)+cmp(10,i)*cphirec(8,i))
         vyy =vyy - cmp(3,i)*cphirec(3,i) - two*cmp(6,i)*cphirec(6,i)
     &            - cmp(8,i)*cphirec(8,i) - cmp(10,i)*cphirec(10,i)
         vyz =vyz - half*(cmp(4,i)*cphirec(3,i)+cmp(3,i)*cphirec(4,i))
     &            - (cmp(6,i)+cmp(7,i))*cphirec(10,i)
     &            - half*cmp(10,i)*(cphirec(6,i)+cphirec(7,i))
     &            - half*(cmp(8,i)*cphirec(9,i)+cmp(9,i)*cphirec(8,i))
         vzz =vzz - cmp(4,i)*cphirec(4,i) - two*cmp(7,i)*cphirec(7,i)
     &            - cmp(9,i)*cphirec(9,i) - cmp(10,i)*cphirec(10,i)
      end do

      !TODO Amoebap  (Add chgflx feature)
c
c     increment the internal virial tensor components
c
!$acc serial async(rec_queue) default(present)
!$acc&       present(e,emrec,g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz)
!$acc&       present(vxx,vxy,vxz,vyy,vyz,vzz)
      emrec    = emrec + half*e
      g_vxx =  g_vxx + vxx
      g_vxy =  g_vxy + vxy
      g_vxz =  g_vxz + vxz
      g_vyy =  g_vyy + vyy
      g_vyz =  g_vyz + vyz
      g_vzz =  g_vzz + vzz
!$acc end serial
      
      else
!$acc serial async(rec_queue) present(emrec,e)
         emrec    = emrec + half*e
!$acc end serial

      end if      ! COMPUTE virial ?

      call torquegpu(npolerecloc,nlocrec2,polerecglob,locrec
     &              ,trqrec,demrec,rec_queue)

      call timer_exit(timer_fmanage,quiet_timers )
c
      if (use_chgflx) then
         def_queue = rec_queue
         if (nproc.eq.1) then
            if (npolerecloc.eq.n) then
!$acc parallel loop async(rec_queue) default(present)
               do i = 1,n; pot(i) = cphirec(1,i); end do
            else
!$acc parallel loop async(rec_queue) default(present)
               do i = 1,npolerecloc
                  iloc = loc(ipole(polerecglob(i)))
                  pot(iloc) = cphirec(1,i)
               end do
            end if
         else
            call prmem_request(potrec,nlocrec)
!$acc parallel loop async(rec_queue) default(present)
            do i = 1,npolerecloc
               iloc = locrec(ipole(polerecglob(i)))
               potrec(iloc) = cphirec(1,i)
            end do
            call commDDrd   (pot,potrec,1,ucComm)
            call commDDd_ext(pot,       1,ucNeig+ucComm)
         end if
         call dcflux2( pot,de_ws1 )
         call adflux2( de_ws1,dem )
         call mem_set(pot,zeror,int(nbloc,mipk),rec_stream)
         if (nproc.ne.1) call mem_set
     &      (potrec,zeror,int(nlocrec,mipk),rec_stream)
      end if
c
      call timer_exit(timer_emrecip,quiet_timers)
      end
