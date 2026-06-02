!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ###################################################################
!     ##                                                               ##
!     ##  subroutine epolar3_group  --  "group" multipoles             ##
!     ##                                                               ##
!     ###################################################################
!
!
!     "empole3_group" calculates the permanent multipoles energy
!     using a subset of atoms in gaz phase
!
!
#include "tinker_precision.h"
#include "tinker_types.h"
module empole3_group_inl
   use tintypes , only: real3,rpole_elt
   include "erfcore_data.inc.f90"
contains
#include "image.inc.f90"
#if defined(SINGLE) | defined(MIXED)
   include "erfcscore.inc.f90"
#else
   include "erfcdcore.inc.f90"
#endif
#include "switch_respa.inc.f90"
#include "pair_mpole1.inc.f90"
end module

subroutine empole3_group
   use action
   use atmlst
   use atoms
   use boxes
   use chgpot
   use domdec
   use energi
   use ewald
   use group
   use iounit
   use inform
   use math
   use mpole
   use polar
   use polpot
   use potent
   use virial
   use mpi
   use utilgpu
   implicit none
   integer i,iglob,iipole,iglobgroup
   real*8 a(3,3)

   if (deb_Path) print '(2X,A)','epolar3_group'
!
!
!     zero out the polarization energy and derivatives
!
!$acc serial async present(emgroup)
   emgroup = 0.0d0
!$acc end serial

!
   if (npolegroup .eq. 0)  return

   call emreal3_group

   return
end

!     #################################################################
!     ##                                                             ##
!     ##  subroutine emreal3_group  --  group multipoles energy      ##
!     ##                                                             ##
!     #################################################################
!
!
!     "emreal3_group" evaluates the
!     multipoles energy of a subset of atoms
!
!
subroutine emreal3_group
   use action  ,only: nep_
   use analyz
   use atoms   ,only: x,y,z
   use chgpot  ,only: dielec,electric
   use atmlst  ,only: poleglobnl
   use atmtyp
   use bound
   use couple
   use cutoff  ,only: shortheal,ewaldshortcut
   use domdec
   use energi  ,only: ep=>ep_r
   use ewald
   use erf_mod
   use inform
   use inter
   use iounit
   use math
   use molcul
   use mpole
   use neigh
   use polar
   use polgrp
   use polpot
   use potent
   use shunt
   use mpi
   use utilgpu
   use timestat
   use group
   use empole3_group_inl
   implicit none
   integer j,k,iglob,kglob
   integer ii,kk,kkk,iipole,kpole,iipolegroup
   integer iglobgroup,kglobgroup,ver,fea
   integer(1) mutik
   logical use_grp,use_chgf,use_lamdyn
   real(t_p) r2,pgamma,damp,alsq2,alsq2n
   real(t_p) f,e
   real(t_p) pdi,pti,poti,potk,delambdae_
   real(t_p) one
   real(t_p) pscale,dscale,uscale
   type(rpole_elt):: ip,kp
   type(real3) :: pos,posi,frc,trqi,trqk
   character*10:: mode
   parameter (one=1.0_ti_p,alsq2=0.0,alsq2n=0.0&
      &, use_grp=.false.,use_chgf=.false.,use_lamdyn=.false.&
      &, mutik=0&
      &, ver=__use_ene__&
      &, fea=__use_mpi__&
      &)

   if (deb_Path) print '(2X,A)','>>> emreal3_group'
!     set conversion factor, cutoff and switching coefficients
!
   f = electric / dielec
!
!     compute the dipole polarization gradient components
!
!$acc parallel loop gang vector_length(32) &
!$acc         present(poleglobgroup,ipolegroup,globglobgroup &
!$acc   ,rpole,thole,pdamp, &
!$acc  loc,x,y,z,pollist,emgroup) &
!$acc         private(ip,posi) &
!$acc         async(def_queue)
   do ii = 1, npolelocgroup
      iipolegroup = poleglobgroup(ii)
      iglobgroup = ipolegroup(iipolegroup)
      iglob = globglobgroup(iglobgroup)
      iipole = pollist(iglob)

      posi%x  = x(iglob)
      posi%y  = y(iglob)
      posi%z  = z(iglob)

      pdi = pdamp(iipole)
      pti = thole(iipole)

      ip%c    = rpole( 1, iipole)
      ip%dx   = rpole( 2, iipole)
      ip%dy   = rpole( 3, iipole)
      ip%dz   = rpole( 4, iipole)
      ip%qxx  = rpole( 5, iipole)
      ip%qxz  = rpole( 7, iipole)
      ip%qxy  = rpole( 6, iipole)
      ip%qyy  = rpole( 9, iipole)
      ip%qyz  = rpole(10, iipole)
      ip%qzz  = rpole(13, iipole)
!
!$acc loop vector private(kp,pos,frc)
      do kkk = iipolegroup+1, npolegroup
         kglobgroup = ipolegroup(kkk)
         kglob = globglobgroup(kglobgroup)
         kpole = pollist(kglob)

         pos%x    = x(kglob) - posi%x
         pos%y    = y(kglob) - posi%y
         pos%z    = z(kglob) - posi%z

         !call image_inl(pos%x,pos%y,pos%z)
         r2       = pos%x**2 + pos%y**2 + pos%z**2

         kp%c     = rpole( 1, kpole)
         kp%dx    = rpole( 2, kpole)
         kp%dy    = rpole( 3, kpole)
         kp%dz    = rpole( 4, kpole)
         kp%qxx   = rpole( 5, kpole)
         kp%qxy   = rpole( 6, kpole)
         kp%qxz   = rpole( 7, kpole)
         kp%qyy   = rpole( 9, kpole)
         kp%qyz   = rpole(10, kpole)
         kp%qzz   = rpole(13, kpole)
!
!     Compute polar interaction
!
         ! compute mpole one interaction
         call duo_mpole(r2,pos%x,pos%y,pos%z,ip,kp,0.0,ewaldshortcut&
            &,shortheal,0.0,f,0.0,0.0,use_grp,one&
            &,use_lamdyn,mutik,0.0,use_chgf&
            &,poti,potk,delambdae_,e,frc,trqi,trqk,ver,fea)

         emgroup = emgroup + e

      end do
   end do

   call emreal3_correct_scale_group

   if (deb_Path) print '(2X,A)','<<< epreal3_group'
end

subroutine emreal3_correct_scale_group
   use action  ,only: nep
   use atmlst  ,only: poleglobnl
   use atoms   ,only: x,y,z
   use chgpot  ,only: dielec,electric
   use cutoff  ,only: shortheal,ewaldshortcut
   use domdec  ,only: rank,loc
   use energi  ,only: ep=>ep_r
   use empole3_group_inl
   use ewald   ,only: aewald
   use inform  ,only: deb_Path
   use math    ,only: sqrtpi
   use mpole   ,only: npolelocnl,ipole,rpole,polelocnl
   use polar   ,only: uind,uinp,thole,pdamp
   use potent  ,only: use_polarshortreal
   use shunt   ,only: off2,off
   use tinheader,only: ti_p
   use utilgpu ,only: def_queue,real3,real6,real3_red,rpole_elt
   use atoms   ,only: x,y,z
   use group

   implicit none

   integer i,k,iglob,kglob,iploc,kploc
   integer iipolegroup,kpolegroup,iglobgroup,kglobgroup
   integer nnelst,ver,fea
   integer ii,iipole,kpole,j,kbis,iga,igb
   integer(1) mutik
   logical use_grp,use_chgf,use_lamdyn
   real(t_p) alsq2,alsq2n,r2,pgamma,damp
   real(t_p) f,e,fgrp
   real(t_p) pdi,pti,poti,potk,delambdae_
   real(t_p) one,mscale
   type(rpole_elt):: ip,kp
   type(real3) :: pos,posi,frc,trqi,trqk

   parameter (one=1.0_ti_p,alsq2=0.0,alsq2n=0.0&
      &, use_grp=.false.,use_chgf=.false.,use_lamdyn=.false.&
      &, fgrp=1.0&
      &, mutik=0&
      &, ver=__use_ene__+__use_sca__&
      &, fea=__use_mpi__&
      &)

   if(deb_Path)&
      &write(*,'(2x,a)') 'epreal3c_correct_scale_group'
!
!     set conversion factor, cutoff and switching coefficients
!
   f      = electric / dielec

!$acc parallel loop gang vector_length(32) &
!$acc         present(ipole,rpole,thole,pdamp,loc, &
!$acc     x,y,z,mcorrect_ik_group, &
!$acc     mcorrect_scale_group,emgroup) &
!$acc     private(pos,ip,kp) &
!$acc         async(def_queue)
   do ii = 1, n_mscale_group
      iipole   = mcorrect_ik_group(ii,1)
      kpole    = mcorrect_ik_group(ii,2)

      mscale   = mcorrect_scale_group(ii)
      if (mscale.eq.0.0) cycle

      iglob    = ipole(iipole)
      kglob    = ipole(kpole)

      iglobgroup = loclocgroup(iglob)
      iipolegroup = pollistgroup(iglobgroup)

      kglobgroup = loclocgroup(kglob)
      kpolegroup = pollistgroup(kglobgroup)

      pos%x    = x(kglob) - x(iglob)
      pos%y    = y(kglob) - y(iglob)
      pos%z    = z(kglob) - z(iglob)

      !call image_inl(pos%x,pos%y,pos%z)
      ! cutoff
      r2       = pos%x**2 + pos%y**2 + pos%z**2


      ip%c     = rpole( 1,iipole)
      ip%dx    = rpole( 2,iipole)
      ip%dy    = rpole( 3,iipole)
      ip%dz    = rpole( 4,iipole)
      ip%qxx   = rpole( 5,iipole)
      ip%qxy   = rpole( 6,iipole)
      ip%qxz   = rpole( 7,iipole)
      ip%qyy   = rpole( 9,iipole)
      ip%qyz   = rpole(10,iipole)
      ip%qzz   = rpole(13,iipole)

      kp%c     = rpole( 1, kpole)
      kp%dx    = rpole( 2, kpole)
      kp%dy    = rpole( 3, kpole)
      kp%dz    = rpole( 4, kpole)
      kp%qxx   = rpole( 5, kpole)
      kp%qxy   = rpole( 6, kpole)
      kp%qxz   = rpole( 7, kpole)
      kp%qyy   = rpole( 9, kpole)
      kp%qyz   = rpole(10, kpole)
      kp%qzz   = rpole(13, kpole)
!
!        compute mpole one interaction
      call duo_mpole(r2,pos%x,pos%y,pos%z,ip,kp,mscale,ewaldshortcut&
         &,shortheal,0.0,f,0.0,0.0,use_grp,fgrp&
         &,use_lamdyn,mutik,0.0,use_chgf&
         &,poti,potk,delambdae_,e,frc,trqi,trqk,ver,fea)
!
!        increment energy
      emgroup       = emgroup + tp2enr(e)
   end do
end

