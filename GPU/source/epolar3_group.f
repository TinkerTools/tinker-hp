c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ###################################################################
c     ##                                                               ##
c     ##  subroutine epolar3_group  --  "group" polarization           ##
c     ##                                                               ##
c     ###################################################################
c
c
c     "epolar3_group" calculates the dipole polarization energy
c     using a subset of atoms in gaz phase
c
c
#include "tinker_precision.h"
#include "tinker_types.h"
      module epolar3group_inl
        contains
#include "convert.f.inc"
#include "image.f.inc"
#include "pair_polar.inc.f"
      end module

      subroutine epolar3_group
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
c
c
c     zero out the polarization energy and derivatives
c
!$acc serial async present(epgroup)
      epgroup = 0.0d0
!$acc end serial

c
      if (npolegroup .eq. 0)  return

      call prmem_request(uindgroup,3,npolegroup,async=.true.)
      call prmem_request(uinpgroup,3,npolegroup,async=.true.)
c
c     rotates all the multipoles of the group
c
      !do i = 1, npolegroup
      !   iglobgroup = ipolegroup(i)
      !   iglob = globglobgroup(iglobgroup)
      !   iipole = pollist(iglob)
      !   call rotmat (iipole,iglob,a)
      !   call rotsite (iipole,a)
      !end do

      call newinduce_group

      call epreal3_group

      return
      end

c     #################################################################
c     ##                                                             ##
c     ##  subroutine epreal3_group  --  group polarization energy    ##
c     ##                                                             ##
c     #################################################################
c
c
c     "epreal3_group" evaluates the 
c     polarization energy due to dipole polarization
c     of a subset of atoms
c
c
      subroutine epreal3_group
      use action  ,only: nep_
      use analyz
      use atoms   ,only: x,y,z
      use chgpot  ,only: dielec,electric
      use atmlst  ,only: poleglobnl
      use atmtyp
      use bound
      use couple
      use cutoff  ,only: shortheal
      use domdec
      use energi  ,only: ep=>ep_r
      use ewald
      use epolar3group_inl
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
      implicit none
      integer j,k,iglob,kglob
      integer ii,kk,kkk,iipole,kpole,iipolegroup
      integer iglobgroup,kglobgroup
      real(t_p) r2,pgamma,pgama1,damp
      real(t_p) f,e
      real(t_p) pdi,pti
      real(t_p) one
      real(t_p) pscale,dscale,uscale
      type(rpole_elt):: ip,kp
      type(real6) :: dpui,dpuk
      type(real3) :: pos,posi
      character*10:: mode
      parameter(one=1.0,pgama1=0.0)

      if (deb_Path) print '(2X,A)','>>> epreal3_group'
c     set conversion factor, cutoff and switching coefficients
c
      f = 0.5d0 * electric / dielec      
c
c
c     compute the dipole polarization gradient components
c
!$acc parallel loop gang vector_length(32)
!$acc&         present(poleglobgroup,ipolegroup,globglobgroup
!$acc&   ,rpole,thole,pdamp,
!$acc&  loc,x,y,z,uindgroup,uinpgroup,pollist,epgroup)
!$acc&         private(ip,dpui,posi)
!$acc&         async(def_queue)
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

         dpui%x  = uindgroup ( 1, iipolegroup)
         dpui%y  = uindgroup ( 2, iipolegroup)
         dpui%z  = uindgroup ( 3, iipolegroup)
         dpui%xx = uinpgroup ( 1, iipolegroup)
         dpui%yy = uinpgroup ( 2, iipolegroup)
         dpui%zz = uinpgroup ( 3, iipolegroup)

c
!$acc loop vector private(dpuk,kp,pos)
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

            dpuk%x   = uindgroup ( 1, kkk)
            dpuk%y   = uindgroup ( 2, kkk)
            dpuk%z   = uindgroup ( 3, kkk)
            dpuk%xx  = uinpgroup ( 1, kkk)
            dpuk%yy  = uinpgroup ( 2, kkk)
            dpuk%zz  = uinpgroup ( 3, kkk)

            pgamma   = min( pti,thole(kpole) )
            damp     = pdi * pdamp (kpole)
c
c     Compute polar interaction
c
            call epolar3_couple(dpui,ip,dpuk,kp,r2,pos,
     &               0._ti_p,0._ti_p,0._ti_p,pgamma,damp,.false.,f,
     &               off,shortheal,1.0_ti_p,
     &               e,.false.,.false.)

            epgroup = epgroup + e

         end do
      end do

      call epreal3c_correct_scale_group
      
      if (deb_Path) print '(2X,A)','<<< epreal3_group'

      return
      end

       subroutine epreal3c_correct_scale_group
      use action  ,only: nep
      use atmlst  ,only: poleglobnl
      use atoms   ,only: x,y,z
      use chgpot  ,only: dielec,electric
      use cutoff  ,only: shortheal
      use domdec  ,only: rank,loc
      use energi  ,only: ep=>ep_r
      use epolar3group_inl
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
      integer :: iipolegroup,kpolegroup
      integer :: iglobgroup,kglobgroup
      integer nnelst
      integer ii,iipole,kpole
      integer j,kbis
      integer iga,igb
      real(t_p) alsq2,alsq2n
      real(t_p) r2,pgamma,pgama1,damp
      real(t_p) f,e
      real(t_p) pdi,pti
      real(t_p) one
      real(t_p) pscale,dscale,uscale
      real(t_p) fgrp

      type(rpole_elt):: ip,kp
      type(real6) :: dpui,dpuk
      type(real3) :: pos,posi

      parameter(one=1.0_ti_p)
      parameter(pgama1=0.0)

      if(deb_Path)
     &   write(*,'(2x,a)') 'epreal3c_correct_scale_group'
c
c     set conversion factor, cutoff and switching coefficients
c
      f      = 0.5_ti_p * electric / dielec

!$acc parallel loop gang vector_length(32)
!$acc&         present(ipole,rpole,thole,pdamp,loc,
!$acc&     x,y,z,uindgroup,uinpgroup,dpucorrect_ik_group,
!$acc&     dpucorrect_scale_group,epgroup)
!$acc&     private(pos,ip,kp,dpui,dpuk)
!$acc&         async(def_queue)
      do ii = 1, n_dpuscale_group
         iipole   = dpucorrect_ik_group(2*(ii-1)+1)
         kpole    = dpucorrect_ik_group(2*(ii-1)+2)

         !dscale   = dpucorrect_scale(3*(ii-1)+1)
         pscale   = dpucorrect_scale_group(3*(ii-1)+2)
         !uscale   = dpucorrect_scale(3*(ii-1)+3)
         if (pscale.eq.0.0) cycle

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

         pdi      = pdamp(iipole)
         pti      = thole(iipole)

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

         dpui%x   = uindgroup ( 1,iipolegroup)
         dpui%y   = uindgroup ( 2,iipolegroup)
         dpui%z   = uindgroup ( 3,iipolegroup)
         dpui%xx  = uinpgroup ( 1,iipolegroup)
         dpui%yy  = uinpgroup ( 2,iipolegroup)
         dpui%zz  = uinpgroup ( 3,iipolegroup)

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

         dpuk%x   = uindgroup ( 1, kpolegroup)
         dpuk%y   = uindgroup ( 2, kpolegroup)
         dpuk%z   = uindgroup ( 3, kpolegroup)
         dpuk%xx  = uinpgroup ( 1, kpolegroup)
         dpuk%yy  = uinpgroup ( 2, kpolegroup)
         dpuk%zz  = uinpgroup ( 3, kpolegroup)

         pgamma   = min( pti,thole(kpole) )
         damp     = pdi * pdamp (kpole)
        
c
c     Compute polar interaction
c
         call epolar3_couple(dpui,ip,dpuk,kp,r2,pos,
     &            0._ti_p,0._ti_p,0._ti_p,pgamma,damp,.false.,f,
     &            off,shortheal,pscale,e,
     &            .false.,.true.)
c
c     increment energy
c
         epgroup       = epgroup + tp2enr(e)
      end do
c
      end
