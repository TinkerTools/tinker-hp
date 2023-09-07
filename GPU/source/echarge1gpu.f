c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine echarge1  --  charge-charge energy & derivs  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "echarge1" calculates the charge-charge interaction energy
c     and first derivatives with respect to Cartesian coordinates
c
c
#include "tinker_macro.h"
      module echarge1gpu_inl
        use utilgpu , only: real3,real3_red,rpole_elt
        implicit none
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
#include "atomicOp.inc.f"
#include "switch_respa.f.inc"
#include "groups.inc.f"
#include "pair_charge.f.inc"
      end module

      subroutine echarge1gpu
      use potent ,only: use_lambdadyn
      implicit none
c
c     choose the method for summing over pairwise interactions
c
      if (use_lambdadyn) then
        call elambdacharge1cgpu
      else
        call echarge1cgpu
      end if
c
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine echarge1c  --  Ewald charge derivs via list  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "echarge1c" calculates the charge-charge interaction energy
c     and first derivatives with respect to Cartesian coordinates
c     using a particle mesh Ewald summation and a neighbor list
c
c
      subroutine echarge1cgpu
      use atmlst
      use atoms
      use bound
      use boxes
      use charge
      use chgpot
      use deriv
      use echarge1gpu_inl
      use energi
      use ewald
      use domdec
      use iounit
      use inform
      use interfaces,only: ecreal1d_p
      use inter
      use math
      use neigh     ,only: clst2_enable
      use potent
      use timestat
      use tinheader
      use usage
      use utilgpu
      use virial
      use mpi
      implicit none
      integer ii,i,iglob,iichg
      real(r_p) e
      real(t_p) de,term
      real(t_p) f,fs
      real(t_p) xd,yd,zd
      real(t_p) dedx,dedy,dedz

      if (nion.eq.0) return

      if (deb_Path)  write(*,'(1x,a)') 'echarge1cgpu'
      call timer_enter(timer_echarge)
c
c     zero out the Ewald summation energy and derivatives
c
      if (calc_e) then
!$acc enter data create(e) async(rec_queue)
!$acc serial async(rec_queue) present(e)
      e     = 0.0
!$acc end serial
      end if
c
c     compute the Ewald self-energy term over all the atoms
c
      if ((.not.(use_pmecore)).or.(use_pmecore).and.(rank.le.ndir-1))
     $   then
      if (use_cself) then
        call timer_enter(timer_other)
        f  = electric / dielec
        fs = -f * aewald / sqrtpi
        if (calc_e) then
!$acc parallel loop async(rec_queue)
!$acc&         present(chgglob,pchg,e)
           do ii = 1, nionloc
              iichg = chgglob(ii)
              e = e + fs * pchg(iichg)**2
           end do
        end if
c
c     compute the cell dipole boundary correction term
c
        if (boundary .eq. 'VACUUM') then
!$acc wait
           xd = 0.0
           yd = 0.0
           zd = 0.0
!$acc parallel loop default(present)
           do ii = 1, nionloc
             iichg = chgglob(ii)
             iglob = iion(iichg)
             i  = loc(iglob)
             xd = xd + pchg(iichg)*x(iglob)
             yd = yd + pchg(iichg)*y(iglob)
             zd = zd + pchg(iichg)*z(iglob)
           end do
           term = (2.0_ti_p/3.0_ti_p) * f * (pi/volbox)
           if (calc_e) then
!$acc serial present(e)
           e = e + term * (xd*xd+yd*yd+zd*zd)
!$acc end serial
           end if
!$acc parallel loop default(present)
           do ii = 1, nionloc
              iichg = chgglob(ii)
              iglob = iion(iichg)
              i     = loc(iglob)
              de    = 2.0 * term * pchg(iichg)
              dedx  = de * xd
              dedy  = de * yd
              dedz  = de * zd
              dec(1,i) = dec(1,i) + dedx
              dec(2,i) = dec(2,i) + dedy
              dec(3,i) = dec(3,i) + dedz
           end do
        end if
        call timer_exit( timer_other,quiet_timers )
      end if
      end if

      if (clst2_enable) call set_ChgData_CellOrder(.false.)
c
c     compute the real space part of the Ewald summation
c
      if ((.not.(use_pmecore)).or.(use_pmecore).and.(rank.le.ndir-1))
     $   then
         if (use_creal) then
            call timer_enter(timer_real)
            call ecreal1d_p
            call timer_exit( timer_real )
         end if
      end if
c
c     compute the reciprocal space part of the Ewald summation
c
      if ((.not.(use_pmecore)).or.(use_pmecore).and.(rank.gt.ndir-1))
     $   then
         if (use_crec) then
            call timer_enter(timer_rec)
            call ecrecip1gpu
            call timer_exit(timer_rec )
         end if
      end if
c
      if (calc_e) then
!$acc serial present(ecrec,e,ec,ec_r) async(rec_queue)
         ec = ec + e + ecrec + enr2en( ec_r )
!$acc end serial
!$acc exit data delete(e) async(rec_queue)
      end if

      call timer_exit(timer_echarge)
      end
c
c     subroutine elambdacharge1c : charge electrostatic interactions during lambda dynamics
c
      subroutine elambdacharge1cgpu
      use atmlst
      use atoms
      use bound
      use boxes
      use charge
      use chgpot
      use deriv
      use echarge1gpu_inl
      use energi
      use ewald
      use domdec
      use group
      use iounit
      use interfaces
      use inter
      use inform
      use math
      use mutant
      use neigh     ,only: clst2_enable
      use potent
      use timestat
      use tinMemory
      use usage
      use utilgpu
      use virial
      use mpi
      use potent
      use sizes
      implicit none
      integer ii,i,j,iglob,iichg,ierr,altopt
      integer(mipk) siz8
      real(r_p) e,de
      real(t_p) f,fs,term
      real(r_p) xd,yd,zd
      real(r_p) xdtemp,ydtemp,zdtemp
      real(r_p) dedx,dedy,dedz,zero_m
      real(t_p) elambdatemp
      real(r_p), allocatable :: delambdarec0(:,:),delambdarec1(:,:)
      real(r_p) :: elambdarec0,elambdarec1,qtemp
      real(r_p) :: g_vxx_temp,g_vxy_temp,g_vxz_temp
      real(r_p) :: g_vyy_temp,g_vyz_temp,g_vzz_temp
      real(r_p) :: g_vxx_1,g_vxy_1,g_vxz_1
      real(r_p) :: g_vyy_1,g_vyz_1,g_vzz_1
      real(r_p) :: g_vxx_0,g_vxy_0,g_vxz_0
      real(r_p) :: g_vyy_0,g_vyz_0,g_vzz_0
      parameter( zero_m=0.0 
#ifdef _OPENACC
     &         , altopt=0 
#else
     &         , altopt=1
#endif
     &         )
c
      if (nion .eq. 0)  return
c
      if (deb_Path)  write(*,'(1x,a)') 'elambdacharge1cgpu'
      call timer_enter(timer_echarge)
c
      allocate (delambdarec0(3,nlocrec2))
      allocate (delambdarec1(3,nlocrec2))
!$acc enter data create(delambdarec0,delambdarec1
!$acc&     ,elambdarec0,elambdarec1
!$acc&     ,g_vxx_temp,g_vxy_temp,g_vxz_temp
!$acc&     ,g_vyy_temp,g_vyz_temp,g_vzz_temp) async
      elambdatemp = elambda  
c
c     zero out the Ewald summation energy and derivatives
c
c!$acc serial async present(delambdae)
c      delambdae = 0.0
c!$acc end serial
c
c     set Ewald coefficient
c
      aewald = aeewald
c
c     compute the Ewald self-energy term over all the atoms
c
      if ((.not.(use_pmecore)).or.(use_pmecore).and.(rank.le.ndir-1))
     $   then
      if (use_cself) then
         f  = electric / dielec
         fs = -f * aewald / sqrtpi
!$acc parallel loop default(present) present(delambdae,ec) async
!$acc&         reduction(+:delambdae,ec)
         do ii = 1, nionloc
           iichg = chgglob(ii)
           iglob = iion(iichg)
           ec    = ec + fs * pchg(iichg)**2
           if (mut(iglob)) then
             qtemp     =  pchg_orig(iichg)
             delambdae = delambdae + fs*2.0*elambda*qtemp**2
           end if
         end do
c
c     compute the cell dipole boundary correction term
c
        if (boundary .eq. 'VACUUM') then
#ifdef _OPENACC
           __TINKER_FATAL__
#endif
           xd = 0.0
           yd = 0.0
           zd = 0.0
           xdtemp = 0.0
           ydtemp = 0.0
           zdtemp = 0.0
           do ii = 1, nionloc
             iichg = chgglob(ii)
             iglob = iion(iichg)
             i = loc(iglob)
             xd = xd + pchg(iichg)*x(iglob)
             yd = yd + pchg(iichg)*y(iglob)
             zd = zd + pchg(iichg)*z(iglob)
             if (mut(iglob)) then
               qtemp = pchg_orig(iichg)
               xdtemp = xdtemp + qtemp*x(iglob)
               ydtemp = ydtemp + qtemp*y(iglob)
               zdtemp = zdtemp + qtemp*z(iglob)
             end if
           end do
           call MPI_ALLREDUCE(MPI_IN_PLACE,xd,1,MPI_RPREC,MPI_SUM,
     $                        comm_dir,ierr)
           call MPI_ALLREDUCE(MPI_IN_PLACE,yd,1,MPI_RPREC,MPI_SUM,
     $                        comm_dir,ierr)
           call MPI_ALLREDUCE(MPI_IN_PLACE,zd,1,MPI_RPREC,MPI_SUM,
     $                        comm_dir,ierr)
           term = (2.0/3.0) * f * (pi/volbox)
           if (rank.eq.0) then
              ec = ec + term * (xd*xd+yd*yd+zd*zd)
           end if
           delambdae = delambdae + term*(xdtemp**2+ydtemp**2+zdtemp**2)
           do ii = 1, nionloc
              iichg = chgglob(ii)
              iglob = iion(iichg)
              i     = loc(iglob)
              de    = 2.0 * term * pchg(iichg)
              dedx  = de * xd
              dedy  = de * yd
              dedz  = de * zd
              dec(1,i) = dec(1,i) + dedx
              dec(2,i) = dec(2,i) + dedy
              dec(3,i) = dec(3,i) + dedz
           end do
        end if
      end if
      end if
      if (clst2_enable) call set_ChgData_CellOrder(.false.)
c
c     compute the real space part of the Ewald summation
c
      if ((.not.(use_pmecore)).or.(use_pmecore).and.(rank.le.ndir-1))
     $   then
        if (use_creal) then
           call ecreal1d_p
        end if
      end if
c
c     compute the reciprocal space part of the Ewald summation
c
      if ((.not.(use_pmecore)).or.(use_pmecore).and.(rank.gt.ndir-1))
     $  then
        call timer_enter(timer_rec)
        if (use_crec) then
c
c         the reciprocal part is interpolated between 0 and 1
c
          siz8 = 3*nlocrec2
!$acc serial async present(ecrec,g_vxx_temp,g_vxy_temp,g_vxz_temp,
!$acc&  g_vyy_temp,g_vyz_temp,g_vzz_temp,
!$acc& g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz)
          ecrec = 0.0
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
          call mem_set(decrec,zero_m,siz8,rec_stream)

          elambda = 0.0
          call altelec(altopt)
          if (elambda.lt.1.0) then
            call ecrecip1gpu
          end if
!$acc serial async present(elambdarec0,ecrec,
!$acc& g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz)
          elambdarec0  = ecrec
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

          ecrec = 0.0
!$acc end serial
          call mem_move(delambdarec0,decrec,siz8,rec_stream)
          call mem_set(decrec,zero_m,siz8,rec_stream)

          elambda = 1.0
          call altelec(altopt)
          if (elambda.gt.0.0) then
            call ecrecip1gpu
          end if
!$acc serial async present(elambdarec1,ecrec,
!$acc& g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz)
          elambdarec1  = ecrec
          g_vxx_1 = g_vxx
          g_vxy_1 = g_vxy
          g_vxz_1 = g_vxz
          g_vyy_1 = g_vyy
          g_vyz_1 = g_vyz
          g_vzz_1 = g_vzz
!$acc end serial
          call mem_move(delambdarec1,decrec,siz8,rec_stream)

          elambda   = elambdatemp
!$acc wait
!$acc serial async present(elambdarec0,elambdarec1,ecrec,delambdae,
!$acc& g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz,g_vxx_temp,g_vxy_temp,
!$acc& g_vxz_temp,g_vyy_temp,g_vyz_temp,g_vzz_temp)
          ecrec     = (1.0-elambda)*elambdarec0 + elambda*elambdarec1
          g_vxx = g_vxx_temp + (1.0-elambda)*g_vxx_0+elambda*g_vxx_1
          g_vxy = g_vxy_temp + (1.0-elambda)*g_vxy_0+elambda*g_vxy_1
          g_vxz = g_vxz_temp + (1.0-elambda)*g_vxz_0+elambda*g_vxz_1
          g_vyy = g_vyy_temp + (1.0-elambda)*g_vyy_0+elambda*g_vyy_1
          g_vyz = g_vyz_temp + (1.0-elambda)*g_vyz_0+elambda*g_vyz_1
          g_vzz = g_vzz_temp + (1.0-elambda)*g_vzz_0+elambda*g_vzz_1
          delambdae = delambdae + elambdarec1-elambdarec0
!$acc end serial
!$acc parallel loop async collapse(2) default(present)
          do i = 1,nlocrec2; do j = 1,3
             decrec(j,i) = (1-elambda)*delambdarec0(j,i)
     &                    +   elambda *delambdarec1(j,i)
          end do; end do
c
c         reset lambda to initial value
c
          call altelec(altopt)
        end if
        call timer_exit(timer_rec )
      end if
c
      if (calc_e) then
!$acc serial present(ecrec,ec,ec_r,delambdae) async(rec_queue)
c!$acc& present(elambdarec0,elambdarec1)
c         print*, ec,ec_r,ecrec,elambdarec1,elambdarec0,elambda
         ec = ec + ecrec + enr2en( ec_r )
!$acc end serial
      end if
c
!$acc update host(delambdae) async
c
!$acc exit data delete(delambdarec0,delambdarec1
!$acc&    ,elambdarec0,elambdarec1
!$acc&     ,g_vxx_temp,g_vxy_temp,g_vxz_temp
!$acc&     ,g_vyy_temp,g_vyz_temp,g_vzz_temp) async
      deallocate(delambdarec0,delambdarec1)

      call timer_exit(timer_echarge)
      end

c
c     "ecreal1dgpu" evaluates the real space portion of the Ewald sum
c     energy and forces due to atomic charge interactions, using a neighbor list
c
      subroutine ecreal1dgpu
      use atmlst
      use atoms
      use bound
      use boxes
      use charge
      use chgpot
      use couple
      use cutoff
      use deriv
      use domdec
      use echarge1gpu_inl
      use energi ,only: ec=>ec_r
      use ewald
      use group
      use iounit
      use inform
      use inter
      use math
      use molcul
      use mutant
      use neigh
      use potent
      use shunt
      use timestat
      use usage
      use tinTypes
      use utilgpu
      use virial
      use mpi
      implicit none
      integer    i,j,k,iichg,iglob,ii,kkk,kglob,kkchg,ver,ver1,fea
      integer   ,pointer:: lst(:,:),nlst(:)
      integer(1) mutik,muti,zero1
      real(t_p)  e,delambdae_,scale
      real(t_p)  f,fi,fi_,fik_,fik,pchgk,r,r2,rew,rb,rb2
      real(t_p)  xi,yi,zi,xr,yr,zr
      real(t_p)  loff2,scut,fgrp
      type(real3)   ded
      character*11  mode

      parameter( ver  = __use_grd__+__use_ene__+__use_vir__
     &         , ver1 = ver        +__use_sca__
     &         , zero1 = 0
     &         )

      if (deb_Path) write(*,'(2x,a)') 'ecreal1dgpu'

      fea = __use_mpi__
      if (use_lambdadyn) fea = fea + __use_lambdadyn__
      if (use_group)     fea = fea + __use_groups__

      !set conversion factor, cutoff and switching coefficients
      f    = electric / dielec
      if (use_cshortreal) then
         mode  = 'SHORTEWALD'
         call switch (mode)
         loff2 = 0
         scut  = off
         lst   => shortelst
         nlst  => nshortelst
         fea   = fea + __use_shortRange__
      else if (use_clong) then
         mode  = 'EWALD'
         call switch (mode)
         loff2 = (chgshortcut-shortheal)**2
         scut  = chgshortcut
         lst   => elst
         nlst  => nelst
         fea   = fea + __use_longRange__
      else
         mode  = 'EWALD'
         call switch (mode)
         loff2 = 0
         scut  = chgshortcut
         lst   => elst
         nlst  => nelst
      end if
c
c     Scaling interaction correction subroutines
c
!$acc parallel loop gang vector async(dir_queue)
!$acc&         present(ccorrect_ik,ccorrect_scale,loc,x,y,z,wgrp
!$acc&   ,grplist,mutInt,dec,ec,delambdae,chglist,pchg,pchg_orig
!$acc&   ,g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz)
!$acc&         private(ded)
!$acc&   reduction(+:ec,delambdae,g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz)
      do ii = 1, n_cscale
         iglob = ccorrect_ik(ii,1)
         kglob = ccorrect_ik(ii,2)
         scale = ccorrect_scale(2*ii+1)
         if (use_lambdadyn) then
         fi    = pchg(chglist(iglob))
         fik   = f*fi*pchg(chglist(kglob))
         fi    = pchg_orig(chglist(iglob))
         fik_  = f*fi*pchg_orig(chglist(kglob))
         mutik = mutInt(iglob)+mutInt(kglob)
         else
         fik   = f*ccorrect_scale(2*ii+2)
         end if
         xi    = x(iglob)
         yi    = y(iglob)
         zi    = z(iglob)
         if (use_group) then
            call groups2_inl(fgrp,iglob,kglob,ngrp,grplist,wgrp)
            scale = 1.0 - fgrp*(1.0-scale)
         end if
c
c     compute the energy contribution for this interaction
c
         xr    = xi - x(kglob)
         yr    = yi - y(kglob)
         zr    = zi - z(kglob)
c
c     find energy for interactions within real space cutoff
c
         call image_inl (xr,yr,zr)
         r2 = xr*xr + yr*yr + zr*zr
         if (r2.gt.loff2 .and. r2.le.off2) then
            i  = loc(iglob)
            k  = loc(kglob)
            call charge_couple(r2,xr,yr,zr,ebuffer,fik_,fik,aewald
     &                        ,scale,mutik,use_lambdadyn,shortheal
     &                        ,scut,elambda,delambdae_,e,ded,ver1,fea)
 
            if (use_lambdadyn) delambdae = delambdae+delambdae_
           !increment the overall energy and derivative expressions
            ec = ec + tp2enr(e)
            call atomic_add( dec(1,i),ded%x )
            call atomic_add( dec(2,i),ded%y )
            call atomic_add( dec(3,i),ded%z )
            call atomic_sub( dec(1,k),ded%x )
            call atomic_sub( dec(2,k),ded%y )
            call atomic_sub( dec(3,k),ded%z )
           !increment the internal virial tensor components
            g_vxx = g_vxx + xr * ded%x
            g_vxy = g_vxy + yr * ded%x
            g_vxz = g_vxz + zr * ded%x
            g_vyy = g_vyy + yr * ded%y
            g_vyz = g_vyz + zr * ded%y
            g_vzz = g_vzz + zr * ded%z
         end if
      end do
c
c     compute the real space Ewald energy and first derivatives
c
!$acc parallel loop gang vector_length(32) async(dir_queue)
!$acc&         present(chgglobnl,iion,loc,x,y,z,pchg,pchg_orig
!$acc&        ,nlst,lst,grplist,wgrp,mutInt
!$acc&        ,dec,ec,delambdae,g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz)
!$acc&   reduction(+:ec,delambdae,g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz)
      do ii = 1, nionlocnl
         iichg = chgglobnl(ii)
         iglob = iion(iichg)
         i     = loc(iglob)
         xi    = x(iglob)
         yi    = y(iglob)
         zi    = z(iglob)
         fi    = f*pchg(iichg)
         if (use_lambdadyn) then
            muti = mutInt(iglob)
            fi_  = f*pchg_orig(iichg)
         end if
!$acc loop vector private(ded)
!$acc&   reduction(+:ec,delambdae,g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz)
         do kkk = 1, nlst(ii)
            kkchg = lst(kkk,ii)
            if (kkchg.eq.0) cycle
            kglob = iion(kkchg)
            if (use_lambdadyn) then
               mutik = muti+mutInt(kglob)
            end if
c
c     compute the energy contribution for this interaction
c
            xr = xi - x(kglob)
            yr = yi - y(kglob)
            zr = zi - z(kglob)
            call image_inl (xr,yr,zr)
            r2 = xr*xr + yr*yr + zr*zr
c
c     find energy for interactions within real space cutoff
c
            if (r2.gt.loff2 .and. r2.le.off2) then
               fik  = fi*pchg(kkchg)
               k    = loc(kglob)
               if (use_lambdadyn.and.mutik.ne.zero1)
     &            fik_  = fi_*pchg_orig(kkchg)

               ! Apply pair group scaling
               if (use_group) then
                  call groups2_inl(fgrp,iglob,kglob,ngrp,grplist,wgrp)
                  fgrp = 1.0-fgrp
               else
                  fgrp = 0.0
               end if

               ! Compute charge pairwise interaction
               call charge_couple(r2,xr,yr,zr,ebuffer,fik_,fik,aewald
     &                           ,fgrp,mutik,use_lambdadyn,shortheal
     &                           ,scut,elambda,delambdae_,e,ded,ver,fea)
 
               if (use_lambdadyn) delambdae = delambdae+delambdae_
              !increment the overall energy and derivative expressions
               ec   = ec + tp2enr(e)
               call atomic_add( dec(1,i),ded%x )
               call atomic_add( dec(2,i),ded%y )
               call atomic_add( dec(3,i),ded%z )
               call atomic_sub( dec(1,k),ded%x )
               call atomic_sub( dec(2,k),ded%y )
               call atomic_sub( dec(3,k),ded%z )
              !increment the internal virial tensor components
               g_vxx = g_vxx + xr * ded%x
               g_vxy = g_vxy + yr * ded%x
               g_vxz = g_vxz + zr * ded%x
               g_vyy = g_vyy + yr * ded%y
               g_vyz = g_vyz + zr * ded%y
               g_vzz = g_vzz + zr * ded%z
            end if
         end do
      end do
      end

#ifdef _CUDA
      subroutine ecreal1d_cu
      use atmlst
      use atoms
      use bound
      use boxes
      use charge
      use chgpot
      use couple
      use cutoff
      use deriv
      use domdec
      use echargecu
      use energi , only: ec=>ec_r,calc_e
      use ewald
      use group
      use iounit
      use inform
      use inter
      use math
      use molcul
      use mpi
      use mutant
      use neigh  , iion_s=>celle_glob,chg_s=>celle_chg
     &           , loc_s=>celle_loc, ieblst_s=>ieblst, eblst_s=>eblst
     &           , x_s=>celle_x, y_s=>celle_y, z_s=>celle_z
      use potent
      use shunt
      use timestat
      use usage
      use utilcu  ,only: BLOCK_DIM,check_launch_kernel 
      use utilgpu ,only: def_queue,dir_queue,dir_stream,def_stream
     &            , vred_buff,lam_buff,ered_buff=>ered_buf1,nred_buff
     &            , reduce_energy_virial, reduce_energy_action
     &            , reduce_buffer, RED_BUFF_SIZE
      use virial
      implicit none
      integer i,szcik
      real(t_p) f
      real(t_p) loff2,scut
      real(t_p) p_xbeg,p_xend,p_ybeg,p_yend,p_zbeg,p_zend
      logical,save:: first_in=.true.
      integer,save:: gS
      character*11 mode
      logical,parameter::dyn_gS=.true.

      if (deb_Path) write(*,'(2X,A)') 'ecreal1d_cu'

      if (first_in) then
         call cudaMaxGridSize("ecreal1_kcu",gS)
         first_in=.false.
      end if
      if (dyn_gS) 
     &   gs = merge(nionlocnlb2_pair,min(nionlocnlb2_pair/8,2**14)
     &             ,nionlocnlb2_pair.lt.257)

      p_xbeg = xbegproc(rank+1)
      p_xend = xendproc(rank+1)
      p_ybeg = ybegproc(rank+1)
      p_yend = yendproc(rank+1)
      p_zbeg = zbegproc(rank+1)
      p_zend = zendproc(rank+1)

      !set conversion factor, cutoff and switching coefficients
      f    = electric / dielec
      if (use_cshortreal) then
         mode  = 'SHORTEWALD'
         call switch (mode)
         loff2 = 0
         scut  = off
      else if (use_clong) then
         mode  = 'EWALD'
         call switch (mode)
         loff2 = (chgshortcut-shortheal)**2
         scut  = chgshortcut
      else
         mode  = 'EWALD'
         call switch (mode)
         loff2 = 0
         scut  = chgshortcut
      end if

      def_stream = dir_stream
      def_queue  = dir_queue
      szcik      = size(ccorrect_ik,1)

!$acc host_data use_device(iion_s,chg_s,loc_s,ieblst_s,eblst_s
!$acc&    ,x_s,y_s,z_s,pchg,pchg_orig,mutInt,grplist,wgrp
!$acc&    ,dec,ered_buff,vred_buff,lam_buff,nred_buff
!$acc&    ,ccorrect_ik,ccorrect_scale,chglist,loc,x,y,z
!$acc&    )

      call ecreal1_kcu<<<gS,BLOCK_DIM,0,def_stream>>>
     &     ( iion_s,chg_s,loc_s,ieblst_s
     &     , eblst_s(2*nionlocnlb_pair+1)
     &     , nionlocnlb,nionlocnlb2_pair,nionbloc,n
     &     , x_s,y_s,z_s,pchg,pchg_orig,mutInt
     &     , off2,loff2,scut,shortheal,f,aewald,ebuffer
     &     , elambda,use_lambdadyn
     &     , dec,ered_buff,vred_buff,lam_buff,nred_buff
     &     , use_group, grplist, wgrp
     &     , ccorrect_ik,ccorrect_scale,chglist,x,y,z,loc,n_cscale,szcik
     &     , p_xbeg,p_xend,p_ybeg,p_yend,p_zbeg,p_zend
     &     )
      call check_launch_kernel(" ecreal1d_core_cu")

!$acc end host_data

      if (calc_e.or.use_virial) then
      call reduce_energy_virial(ec,g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz
     &                         ,ered_buff,def_queue)
      end if

      if (use_lambdadyn)
     &   call reduce_buffer(lam_buff,RED_BUFF_SIZE,delambdae,def_queue)

      end subroutine

      subroutine ecreal_lj1_cu
      use atmlst
      use atoms
      use bound
      use boxes
      use charge
      use chgpot
      use couple
      use cudafor, only: dim3
      use cutoff , only: ewaldcut,shortheal
      use deriv
      use domdec
      use eChgLjcu
      use energi , only: ec=>ec_r,calc_e
      use ewald
      use group
      use iounit
      use inform
      use inter
      use kvdws     ,only: radv,epsv
      use kvdwpr    ,only: vdwpr_l
      use math
      use molcul
      use mpi
      use neigh  , iion_s=>celle_glob,chg_s=>celle_chg
     &           , loc_s=>celle_loc, ieblst_s=>ieblst, eblst_s=>eblst
     &           , x_s=>celle_x, y_s=>celle_y, z_s=>celle_z
      use potent
      use shunt   ,only: c0,c1,c2,c3,c4,c5,off2,off,cut2,cut
      use timestat
      use usage
      use utilcu  ,only: BLOCK_DIM,check_launch_kernel 
     &            , TRP_BLOCK_DIM,transpose_z3fl
      use utilgpu ,only: def_queue,dir_queue,vred_buff
     &            , reduce_energy_virial
     &            , ered_buff=>ered_buf1,dir_stream,def_stream
     &            , transpose_az3,get_GridDim
      use vdw     ,only: ired,kred,jvdw,ivdw,radmin_c
     &            , epsilon_c,radmin,radmin4,epsilon,epsilon4
     &            , nvdwbloc,nvdwlocnl
     &            , nvdwlocnlb,nvdwclass
     &            , nvdwlocnlb_pair,nvdwlocnlb2_pair
      use vdwpot  ,only: vcorrect_ik,vcorrect_scale,n_vscale,dhal,ghal
     &            , radrule_i,epsrule_i,radepsOpt_l
      use vdw_locArray
      use virial
      implicit none
      integer i,gs1,szcik
      real(t_p) f,coff2,rinv,aewaldop
      real(t_p) p_xbeg,p_xend,p_ybeg,p_yend,p_zbeg,p_zend
      logical,save:: first_in=.true.
      integer,save:: gS
      character*11 mode
      logical,parameter::dyn_gS=.true.
      logical v0,v1

      if (deb_Path) write(*,'(2X,A)') 'ecreal_lj1_cu'
c
c     set conversion factor, cutoff and switching coefficients
c
      f    = electric / dielec
      mode = 'VDW'
      call switch (mode)
      rinv = 1.0/(cut-off)

      if (first_in) then
         call cudaMaxGridSize("echg_lj1_kcu_v0",gS)
         first_in=.false.
      end if
      if(dyn_gS) gs = get_GridDim(nionlocnlb2_pair,BLOCK_DIM)

      p_xbeg = xbegproc(rank+1)
      p_xend = xendproc(rank+1)
      p_ybeg = ybegproc(rank+1)
      p_yend = yendproc(rank+1)
      p_zbeg = zbegproc(rank+1)
      p_zend = zendproc(rank+1)
      szcik  = size(ccorrect_ik,1)
      coff2  = ewaldcut**2
      aewaldop = 2d0*aewald/sqrtpi

      def_stream = dir_stream
      def_queue  = dir_queue

      v0 = .not.calc_e.and..not.use_virial.and.ndir.eq.1
     &     .and.radepsOpt_l

!$acc host_data use_device(iion_s,chg_s,loc_s,ieblst_s,eblst_s,
!$acc&   x_s,y_s,z_s,pchg,dec,d_x,d_y,d_z,ered_buff,vred_buff,
!$acc&   cellv_jvdw,jvdw,i12,epsilon_c,radmin_c,radv,epsv,
!$acc&   radmin,epsilon,radmin4,epsilon4,vcorrect_scale,
!$acc&   ccorrect_ik,ccorrect_scale,loc,x,y,z
!$acc&   )

      if (v0) then
      call echg_lj1_kcu_v0<<<gS,BLOCK_DIM,0,def_stream>>>
     &     ( x_s,y_s,z_s,iion_s,loc_s,chg_s
     &     , ieblst_s, eblst_s(2*nionlocnlb_pair+1)
     &     , cellv_jvdw,i12,epsilon_c,radmin_c,radv,epsv
     &     , d_x,d_y,d_z,ered_buff,vred_buff
     &     , nionlocnlb2_pair,n,nionbloc,nionlocnl,nionlocnlb
     &     , nvdwclass,radrule_i,epsrule_i
     &     , cut2,cut,off2,off,rinv,shortheal
     &     , coff2,f,aewald,aewaldop,ebuffer,pchg
     &     , x,y,z,ccorrect_ik,vcorrect_scale,ccorrect_scale,szcik
     &     , n_cscale,loc,jvdw
     &     , radmin,epsilon,radmin4,epsilon4,dec
     &     , p_xbeg,p_xend,p_ybeg,p_yend,p_zbeg,p_zend
     &     )
      call check_launch_kernel(" echg_lj1_kcu_v0")
      else
      call echg_lj1_kcu_v1<<<gS,BLOCK_DIM,0,def_stream>>>
     &     ( x_s,y_s,z_s,iion_s,loc_s,chg_s
     &     , ieblst_s, eblst_s(2*nionlocnlb_pair+1)
     &     , cellv_jvdw,i12,epsilon_c,radmin_c,radv,epsv
     &     , d_x,d_y,d_z,ered_buff,vred_buff
     &     , nionlocnlb2_pair,n,nionbloc,nionlocnl,nionlocnlb
     &     , nvdwclass,radrule_i,epsrule_i
     &     , cut2,cut,off2,off,rinv,shortheal
     &     , coff2,f,aewald,aewaldop,ebuffer,pchg
     &     , x,y,z,ccorrect_ik,vcorrect_scale,ccorrect_scale,szcik
     &     , n_cscale,loc,jvdw
     &     , radmin,epsilon,radmin4,epsilon4,dec
     &     , p_xbeg,p_xend,p_ybeg,p_yend,p_zbeg,p_zend
     &     )
      call check_launch_kernel(" echg_lj1_kcu_v1")
      end if

!$acc end host_data

      if (calc_e.or.use_virial) then
      call reduce_energy_virial(ec,g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz
     &                         ,ered_buff,def_queue)
      end if

      call transpose_az3(d_x,decx,loc_s,nionlocnl,dr_stride,def_queue)
c      block
c      integer    sdx
c      type(dim3) gridS,blockS
c      sdx    = size(d_x)
c      gridS  = dim3((sdx-1)/TRP_BLOCK_DIM+1,1,1)
c      blockS = dim3(TRP_BLOCK_DIM,3,1)
c!$acc host_data use_device(d_x,dec,loc_s)
c      call transpose_z3fl<<<gridS,blockS,0,def_stream>>>
c     &                ( d_x,dec,sdx,nionlocnl,loc_s )
c      call check_launch_kernel(" transpose_z3fl ")
c!$acc end host_data
c      end block

      !call minmaxone(dec,3*nbloc,'dec')

      end subroutine
#endif
c
c
c     ####################################################################
c     ##                                                                ##
c     ##  subroutine ecrecip1gpu  --  PME recip charge energy & derivs  ##
c     ##                                                                ##
c     ####################################################################
c
c
c     "ecrecip1" evaluates the reciprocal space portion of the particle
c     mesh Ewald summation energy and gradient due to partial charges
c
c     literature reference:
c
c     U. Essmann, L. Perera, M. L Berkowitz, T. Darden, H. Lee and
c     L. G. Pedersen, "A Smooth Particle Mesh Ewald Method", Journal
c     of Chemical Physics, 103, 8577-8593 (1995)
c
c     modifications for nonperiodic systems suggested by Tom Darden
c     during May 2007
c
c
      subroutine ecrecip1gpu
      use atmlst
      use atoms
      use bound
      use boxes
      use charge
      use chgpot
      use deriv
      use domdec
      use energi
      use ewald
      use fft
      use inform
      use interfaces ,only: grid_pchg_site_p,grid_pchg_force_p
     &               , ecreal1d_cp, pme_conv_p
#ifdef _OPENACC
     &               , grid_pchg_sitecu
#endif
      use math
      use pme
      use pme1
      use potent
      use timestat
      use utils
      use utilgpu
      use virial
      use mpi
      implicit none
      integer i,j,k
      integer iichg,iglob
      integer isite,iatm
      integer i0,j0,k0
      integer it1,it2,it3
      integer k1,k2,k3
      integer m1,m2,m3
      integer nf1,nf2,nf3
      integer nff
      integer status(MPI_STATUS_SIZE),tag,ierr
      integer nprocloc,rankloc,commloc,proc
      integer ist2,jst2,kst2,ien2,jen2,ken2
      integer istart,iend,jstart,jend,kstart,kend
      integer iloc,iproc
      integer igrd0,jgrd0,kgrd0
      real(t_p) e,term,expterm
      real(t_p) vterm,pterm
      real(t_p) volterm
      real(t_p) f,fi,denom
      real(t_p) hsq,struc2
      real(t_p) de1,de2,de3
      real(t_p) dn1,dn2,dn3
      real(t_p) t1,t2,t3
      real(t_p) dt1,dt2,dt3
      real(t_p) h1,h2,h3
      real(t_p) r1,r2,r3
      integer  ,allocatable:: req(:),reqbcast(:)
c
c     return if the Ewald coefficient is zero
c
      if (aewald .lt. 1.0d-6) return

      if (deb_Path) write(*,'(2x,A)') 'ecrecip1gpu'
c
      if (use_pmecore) then
        nprocloc = nrec
        commloc  = comm_rec
        rankloc  = rank_bis
      else
        nprocloc = nproc
        commloc  = COMM_TINKER
        rankloc  = rank
      end if
c
c     dynamic allocation of local arrays
c
      call mallocMpiGrid
c
      call timer_enter( timer_grid1 )
#ifdef _OPENACC
      if (.not.associated(grid_pchg_site_p,grid_pchg_sitecu))
     &   call bspline_fill_sitegpu(1)
#else
      call bspline_fill_sitegpu(1)
#endif
      call set_to_zero1(qgridin_2d(1,1,1,1,1),
     &     2*n1mpimax*n2mpimax*n3mpimax*(nrec_send+1),rec_queue)

      call grid_pchg_site_p
      call timer_exit ( timer_grid1,quiet_timers )

      ! FFtGridin communication
      call commGridFront( qgridin_2d,r_comm )
      call commGridFront( qgridin_2d,r_wait )

#ifdef _OPENACC
      if (rec_queue.ne.dir_queue) then
         call start_dir_stream_cover
         call ecreal1d_cp
      end if
#endif
c
c     perform the 3-D FFT forward transformation
c
#ifdef _OPENACC
      call cufft2d_frontmpi(qgridin_2d,qgridout_2d,n1mpimax,n2mpimax,
     $                      n3mpimax)
#else
      call   fft2d_frontmpi(qgridin_2d,qgridout_2d,n1mpimax,n2mpimax,
     $                      n3mpimax)
#endif
#if 0
      if (rankloc.eq.0) then
!$acc update host(qgridin_2d,qgridout_2d)
      print*,'gridi norm2',comput_norm(qgridin_2d,size(qgridin_2d),2)
      print*,'grido norm2',comput_norm(qgridout_2d,size(qgridout_2d),2)
      end if
#endif
c
c     use scalar sum to get reciprocal space energy and virial
c
      call pme_conv_p(ecrec,g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz)
c
c     perform the 3-D FFT backward transformation
c
#ifdef _OPENACC
      call cufft2d_backmpi(qgridin_2d,qgridout_2d,n1mpimax,n2mpimax,
     $                     n3mpimax)
#else
      call   fft2d_backmpi(qgridin_2d,qgridout_2d,n1mpimax,n2mpimax,
     $                     n3mpimax)
#endif

      call commGridBack( qgridin_2d,r_comm )
      call commGridBack( qgridin_2d,r_wait )

#if 0
      if (rankloc.eq.0) then
!$acc update host(qgridin_2d,qgridout_2d)
      print*,'grido norm2',comput_norm(qgridout_2d,size(qgridout_2d),2)
      print*,'gridi norm2',comput_norm(qgridin_2d,size(qgridin_2d),2)
      end if
#endif
c
c     get first derivatives of the reciprocal space energy
c
      call timer_enter( timer_grid2 )
      call grid_pchg_force_p
      call timer_exit ( timer_grid2,quiet_timers )
#ifdef _OPENACC
      if (rec_queue.ne.dir_queue) then
         call end_dir_stream_cover
      end if
#endif
      end
