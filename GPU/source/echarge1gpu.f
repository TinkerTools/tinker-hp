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
c      if (use_lambdadyn) then
c        call elambdacharge1cgpu
c      else
        call echarge1cgpu
c      end if
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
      use mutant,  only: mut,elambda
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
      real(r_p) xdtemp,ydtemp,zdtemp
      real(t_p) dedx,dedy,dedz
      real(t_p) qtemp

      if (nion.eq.0) return

      if (deb_Path)  write(*,'(1x,a)') 'echarge1cgpu'
      call timer_enter(timer_echarge)
c
c     zero out the Ewald summation energy and derivatives
c
      if (calc_e) then
!$acc enter data create(e) async(rec_queue)
!$acc serial async(rec_queue) present(e,delambdae)
      e     = 0.0
      delambdae = 0.0
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
           if (use_lambdadyn) then
!$acc parallel loop async(rec_queue)
!$acc&         present(chgglob,iion,pchg,mut,delambdae)
             do ii = 1, nionloc
                iichg = chgglob(ii)
                iglob = iion(iichg)
                 if (mut(iglob)) then
                  qtemp     =  pchg_orig(iichg)
                  delambdae = delambdae + fs*2.0*elambda*qtemp**2
                 end if
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
               if (use_lambdadyn.and.mut(iglob)) then
                 qtemp = pchg_orig(iichg)
                 xdtemp = xdtemp + qtemp*x(iglob)
                 ydtemp = ydtemp + qtemp*y(iglob)
                 zdtemp = zdtemp + qtemp*z(iglob)
               end if
             end do
             term = (2.0_ti_p/3.0_ti_p) * f * (pi/volbox)
             if (calc_e) then
!$acc serial present(e)
             e = e + term * (xd*xd+yd*yd+zd*zd)
!$acc end serial
             end if
             if (use_lambdadyn) then
             delambdae =delambdae + term*(xdtemp**2+ydtemp**2+zdtemp**2)
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
            if (use_lambdadyn.and.use_cshortreal) then
c
c     save delambdae for short range computation
c
!$acc serial async(rec_queue) present(delambdae,delambdaesave)
              delambdaesave = delambdae
!$acc end serial
            end if
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

      if (use_lambdadyn) then
!$acc update host(delambdae,delambdaesave) async(rec_queue)
      end if

      call timer_exit(timer_echarge)
      end
c
c     subroutine elambdacharge1c : charge electrostatic interactions during lambda dynamics
c
c      subroutine elambdacharge1cgpu
c      use atmlst
c      use atoms
c      use bound
c      use boxes
c      use charge
c      use chgpot
c      use deriv
c      use echarge1gpu_inl
c      use energi
c      use ewald
c      use domdec
c      use group
c      use iounit
c      use interfaces
c      use inter
c      use inform
c      use math
c      use mutant
c      use neigh     ,only: clst2_enable
c      use potent
c      use timestat
c      use tinMemory
c      use usage
c      use utilgpu
c      use virial
c      use mpi
c      use potent
c      use sizes
c      implicit none
c      integer ii,i,j,iglob,iichg,ierr,altopt
c      integer(mipk) siz8
c      real(r_p) e,de
c      real(t_p) f,fs,term
c      real(r_p) xd,yd,zd
c      real(r_p) xdtemp,ydtemp,zdtemp
c      real(r_p) dedx,dedy,dedz,zero_m
c      real(t_p) elambdatemp
c      real(r_p), allocatable :: delambdarec0(:,:),delambdarec1(:,:)
c      real(r_p) :: elambdarec0,elambdarec1,qtemp
c      real(r_p) :: g_vxx_temp,g_vxy_temp,g_vxz_temp
c      real(r_p) :: g_vyy_temp,g_vyz_temp,g_vzz_temp
c      real(r_p) :: g_vxx_1,g_vxy_1,g_vxz_1
c      real(r_p) :: g_vyy_1,g_vyz_1,g_vzz_1
c      real(r_p) :: g_vxx_0,g_vxy_0,g_vxz_0
c      real(r_p) :: g_vyy_0,g_vyz_0,g_vzz_0
c      parameter( zero_m=0.0 
c#ifdef _OPENACC
c     &         , altopt=0 
c#else
c     &         , altopt=1
c#endif
c     &         )
cc
c      if (nion .eq. 0)  return
cc
c      if (deb_Path)  write(*,'(1x,a)') 'elambdacharge1cgpu'
c      call timer_enter(timer_echarge)
cc
c      allocate (delambdarec0(3,nlocrec2))
c      allocate (delambdarec1(3,nlocrec2))
c!$acc enter data create(delambdarec0,delambdarec1
c!$acc&     ,elambdarec0,elambdarec1
c!$acc&     ,g_vxx_temp,g_vxy_temp,g_vxz_temp
c!$acc&     ,g_vyy_temp,g_vyz_temp,g_vzz_temp) async
c      elambdatemp = elambda  
cc
cc     zero out the Ewald summation energy and derivatives
cc
cc!$acc serial async present(delambdae)
cc      delambdae = 0.0
cc!$acc end serial
cc
cc     set Ewald coefficient
cc
c      aewald = aeewald
cc
cc     compute the Ewald self-energy term over all the atoms
cc
c      if ((.not.(use_pmecore)).or.(use_pmecore).and.(rank.le.ndir-1))
c     $   then
c      if (use_cself) then
c         f  = electric / dielec
c         fs = -f * aewald / sqrtpi
c!$acc parallel loop default(present) present(delambdae,ec) async
c!$acc&         reduction(+:delambdae,ec)
c         do ii = 1, nionloc
c           iichg = chgglob(ii)
c           iglob = iion(iichg)
c           ec    = ec + fs * pchg(iichg)**2
c           if (mut(iglob)) then
c             qtemp     =  pchg_orig(iichg)
c             delambdae = delambdae + fs*2.0*elambda*qtemp**2
c           end if
c         end do
cc
cc     compute the cell dipole boundary correction term
cc
c        if (boundary .eq. 'VACUUM') then
c#ifdef _OPENACC
c           __TINKER_FATAL__
c#endif
c           xd = 0.0
c           yd = 0.0
c           zd = 0.0
c           xdtemp = 0.0
c           ydtemp = 0.0
c           zdtemp = 0.0
c           do ii = 1, nionloc
c             iichg = chgglob(ii)
c             iglob = iion(iichg)
c             i = loc(iglob)
c             xd = xd + pchg(iichg)*x(iglob)
c             yd = yd + pchg(iichg)*y(iglob)
c             zd = zd + pchg(iichg)*z(iglob)
c             if (mut(iglob)) then
c               qtemp = pchg_orig(iichg)
c               xdtemp = xdtemp + qtemp*x(iglob)
c               ydtemp = ydtemp + qtemp*y(iglob)
c               zdtemp = zdtemp + qtemp*z(iglob)
c             end if
c           end do
c           call MPI_ALLREDUCE(MPI_IN_PLACE,xd,1,MPI_RPREC,MPI_SUM,
c     $                        comm_dir,ierr)
c           call MPI_ALLREDUCE(MPI_IN_PLACE,yd,1,MPI_RPREC,MPI_SUM,
c     $                        comm_dir,ierr)
c           call MPI_ALLREDUCE(MPI_IN_PLACE,zd,1,MPI_RPREC,MPI_SUM,
c     $                        comm_dir,ierr)
c           term = (2.0/3.0) * f * (pi/volbox)
c           if (rank.eq.0) then
c              ec = ec + term * (xd*xd+yd*yd+zd*zd)
c           end if
c           delambdae = delambdae + term*(xdtemp**2+ydtemp**2+zdtemp**2)
c           do ii = 1, nionloc
c              iichg = chgglob(ii)
c              iglob = iion(iichg)
c              i     = loc(iglob)
c              de    = 2.0 * term * pchg(iichg)
c              dedx  = de * xd
c              dedy  = de * yd
c              dedz  = de * zd
c              dec(1,i) = dec(1,i) + dedx
c              dec(2,i) = dec(2,i) + dedy
c              dec(3,i) = dec(3,i) + dedz
c           end do
c        end if
c      end if
c      end if
c      if (clst2_enable) call set_ChgData_CellOrder(.false.)
cc
cc     compute the real space part of the Ewald summation
cc
c      if ((.not.(use_pmecore)).or.(use_pmecore).and.(rank.le.ndir-1))
c     $   then
c        if (use_creal) then
c           call ecreal1d_p
c        end if
c      end if
cc
cc     compute the reciprocal space part of the Ewald summation
cc
c      if ((.not.(use_pmecore)).or.(use_pmecore).and.(rank.gt.ndir-1))
c     $  then
c        call timer_enter(timer_rec)
c        if (use_crec) then
cc
cc         the reciprocal part is interpolated between 0 and 1
cc
c          siz8 = 3*nlocrec2
c!$acc serial async present(ecrec,g_vxx_temp,g_vxy_temp,g_vxz_temp,
c!$acc&  g_vyy_temp,g_vyz_temp,g_vzz_temp,
c!$acc& g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz)
c          ecrec = 0.0
c          g_vxx_temp = g_vxx
c          g_vxy_temp = g_vxy
c          g_vxz_temp = g_vxz
c          g_vyy_temp = g_vyy
c          g_vyz_temp = g_vyz
c          g_vzz_temp = g_vzz
c          g_vxx = 0.0
c          g_vxy = 0.0
c          g_vxz = 0.0
c          g_vyy = 0.0
c          g_vyz = 0.0
c          g_vzz = 0.0
c!$acc end serial
c          call mem_set(decrec,zero_m,siz8,rec_stream)
c
c          elambda = 0.0
c          call altelec(altopt)
c          if (elambda.lt.1.0) then
c            call ecrecip1gpu
c          end if
c!$acc serial async present(elambdarec0,ecrec,
c!$acc& g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz)
c          elambdarec0  = ecrec
c          g_vxx_0 = g_vxx
c          g_vxy_0 = g_vxy
c          g_vxz_0 = g_vxz
c          g_vyy_0 = g_vyy
c          g_vyz_0 = g_vyz
c          g_vzz_0 = g_vzz
c          g_vxx = 0.0
c          g_vxy = 0.0
c          g_vxz = 0.0
c          g_vyy = 0.0
c          g_vyz = 0.0
c          g_vzz = 0.0
c
c          ecrec = 0.0
c!$acc end serial
c          call mem_move(delambdarec0,decrec,siz8,rec_stream)
c          call mem_set(decrec,zero_m,siz8,rec_stream)
c
c          elambda = 1.0
c          call altelec(altopt)
c          if (elambda.gt.0.0) then
c            call ecrecip1gpu
c          end if
c!$acc serial async present(elambdarec1,ecrec,
c!$acc& g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz)
c          elambdarec1  = ecrec
c          g_vxx_1 = g_vxx
c          g_vxy_1 = g_vxy
c          g_vxz_1 = g_vxz
c          g_vyy_1 = g_vyy
c          g_vyz_1 = g_vyz
c          g_vzz_1 = g_vzz
c!$acc end serial
c          call mem_move(delambdarec1,decrec,siz8,rec_stream)
c
c          elambda   = elambdatemp
c!$acc wait
c!$acc serial async present(elambdarec0,elambdarec1,ecrec,delambdae,
c!$acc& g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz,g_vxx_temp,g_vxy_temp,
c!$acc& g_vxz_temp,g_vyy_temp,g_vyz_temp,g_vzz_temp)
c          ecrec     = (1.0-elambda)*elambdarec0 + elambda*elambdarec1
c          g_vxx = g_vxx_temp + (1.0-elambda)*g_vxx_0+elambda*g_vxx_1
c          g_vxy = g_vxy_temp + (1.0-elambda)*g_vxy_0+elambda*g_vxy_1
c          g_vxz = g_vxz_temp + (1.0-elambda)*g_vxz_0+elambda*g_vxz_1
c          g_vyy = g_vyy_temp + (1.0-elambda)*g_vyy_0+elambda*g_vyy_1
c          g_vyz = g_vyz_temp + (1.0-elambda)*g_vyz_0+elambda*g_vyz_1
c          g_vzz = g_vzz_temp + (1.0-elambda)*g_vzz_0+elambda*g_vzz_1
c          delambdae = delambdae + elambdarec1-elambdarec0
c!$acc end serial
c!$acc parallel loop async collapse(2) default(present)
c          do i = 1,nlocrec2; do j = 1,3
c             decrec(j,i) = (1-elambda)*delambdarec0(j,i)
c     &                    +   elambda *delambdarec1(j,i)
c          end do; end do
cc
cc         reset lambda to initial value
cc
c          call altelec(altopt)
c        end if
c        call timer_exit(timer_rec )
c      end if
cc
c      if (calc_e) then
c!$acc serial present(ecrec,ec,ec_r,delambdae) async(rec_queue)
cc!$acc& present(elambdarec0,elambdarec1)
cc         print*, ec,ec_r,ecrec,elambdarec1,elambdarec0,elambda
c         ec = ec + ecrec + enr2en( ec_r )
c!$acc end serial
c      end if
cc
c!$acc update host(delambdae) async
cc
c!$acc exit data delete(delambdarec0,delambdarec1
c!$acc&    ,elambdarec0,elambdarec1
c!$acc&     ,g_vxx_temp,g_vxy_temp,g_vxz_temp
c!$acc&     ,g_vyy_temp,g_vyz_temp,g_vzz_temp) async
c      deallocate(delambdarec0,delambdarec1)
c
c      call timer_exit(timer_echarge)
c      end

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
      use interfaces ,only: grid_pchg_site_p, grid_pchg_force_p
     &               , ecreal1d_cp, pme_conv_p, fphi_chg_site_p
#ifdef _OPENACC
     &               , grid_pchg_sitecu
#endif
      use math
      use mutant, only: mut,elambda
      use pme
      use pme1
      use potent
      use timestat
      use utils
      use utilgpu
      use virial
      use mpi
      implicit none
      integer i,j,k,ii
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
      real(t_p) f1,f2,f3
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
      if (use_lambdadyn) then
        j = max(nionrecloc,1)
c        call prmem_request(cphirec,1,j,async=.false.)
        call prmem_request(fphirec,4,j,async=.false.)
      end if
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
      if (use_lambdadyn) then
         f = electric / dielec
         call fphi_chg_site_p
!$acc parallel loop collapse(2) default(present) async(rec_queue)
         do i = 1, nionrecloc; do j = 1, 4
            fphirec(j,i) = electric * fphirec(j,i)
         end do; end do
c         call fphi_to_cphi_chg_sitegpu (fphirec,cphirec)
!$acc parallel loop default(present) async(rec_queue)
         do i = 1, nionrecloc
            iichg = chgrecglob(i)
            f1 = pchg(iichg)*fphirec(2,i)
            f2 = pchg(iichg)*fphirec(3,i)
            f3 = pchg(iichg)*fphirec(4,i)
            f1 = dble(nfft1) * f1
            f2 = dble(nfft2) * f2
            f3 = dble(nfft3) * f3
            h1 = recip(1,1)*f1 + recip(1,2)*f2 + recip(1,3)*f3
            h2 = recip(2,1)*f1 + recip(2,2)*f2 + recip(2,3)*f3
            h3 = recip(3,1)*f1 + recip(3,2)*f2 + recip(3,3)*f3
            iglob = iion(iichg)
            ii = locrec(iglob)
            decrec(1,ii) = decrec(1,ii) + h1
            decrec(2,ii) = decrec(2,ii) + h2
            decrec(3,ii) = decrec(3,ii) + h3
         end do
!$acc parallel loop async(rec_queue) 
!$acc& present(chgrecglob,iion,pchg,mut,delambdae)
         do i = 1, nionrecloc
            iichg = chgrecglob(i)
            iglob = iion(iichg)
            if (mut(iglob).and.elambda.gt.0) then
               delambdae = delambdae + (pchg(iichg)*fphirec(1,i))/
     &                      elambda
            end if
         end do
      else
         call grid_pchg_force_p
      end if
      call timer_exit ( timer_grid2,quiet_timers )
#ifdef _OPENACC
      if (rec_queue.ne.dir_queue) then
         call end_dir_stream_cover
      end if
#endif
      end
