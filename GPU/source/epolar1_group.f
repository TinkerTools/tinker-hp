c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ###################################################################
c     ##                                                               ##
c     ##  subroutine epolar1_group  --  "group" polarization           ##
c     ##                                                               ##
c     ###################################################################
c
c
c     "epolar1_group" calculates the dipole polarization energy and
c     derivatives with respect to Cartesian coordinates using
c     a subset of atoms in gaz phase
c
c
#include "tinker_macro.h"
      module epolar1group_inl
#include "atomicOp.h.f"
        contains
#include "convert.f.inc"
#include "image.f.inc"
#include "atomicOp.inc.f"
#include "pair_polar.inc.f"
      end module

      subroutine epolar1_group
      use action
      use atmlst
      use atoms
      use boxes
      use chgpot
      use deriv
      use domdec
      use energi
      use ewald
      use group
      use iounit
      use math
      use mpole
      use polar
      use polpot
      use potent
      use virial
      use mpi
      use utilgpu
      use utils
      use inform
      implicit none
      integer i,j,iglob,iipole,iglobgroup
      real(r_p) a(3,3)
c
c
c     zero out the polarization energy and derivatives
c
!$acc serial present(epgroup)
      epgroup = 0.0d0
!$acc end serial
      !if (allocated(depgroup)) deallocate(depgroup)
      !allocate (depgroup(3,natgroup))
      !depgroup = 0d0
      !if (allocated(uindgroup)) deallocate(uindgroup)
      !allocate (uindgroup(3,npolegroup))
      !if (allocated(uinpgroup)) deallocate(uinpgroup)
      !allocate (uinpgroup(3,npolegroup))

      call prmem_request(uindgroup,3,npolegroup, async=.true.)
      call prmem_request(uinpgroup,3,npolegroup, async=.true.)

      call prmem_requestm(depgroup,3,natgroup, async=.true.)
!$acc parallel loop collapse(2) default(present) async
      do i = 1, natgroup
         do j = 1, 3
            depgroup(j,i) = 0
         end do
      end do


c
      if (npolegroup .eq. 0)  return
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
c
      call epreal1_group
c
c     communicate the forces
c
      call commforcesgroup

      return
      end
      
c     #################################################################
c     ##                                                             ##
c     ##  subroutine epreal1_group  --  group polarization derivs    ##
c     ##                                                             ##
c     #################################################################
c
c
c     "epreal1_group" evaluates the 
c     polarization energy and gradient due to dipole polarization
c     of a subset of atoms
c
c
      subroutine epreal1_group
      use atmlst
      use atoms
      use bound
      use chgpot
      use couple
      use deriv
      use domdec
      use energi
      use elec_wspace,only: fix=>r2Work1,fiy=>r2Work2,fiz=>r2Work3
     &               ,trqvec=>r2Work4
      !use erf_mod
      use epolar1group_inl
      use ewald
      use inform     ,only: deb_Path
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
      use virial
      use mpi
      use tinMemory
      use utils
      use utilgpu
      use group
      implicit none
      integer i,j,k,iglob,kglob,kbis
      integer ii,kk,kkk,iipole,kpole
      integer iglobgroup,kglobgroup
      integer ilocgroup,klocgroup,iipolegroup
      integer iploc,kploc
      integer ix,iy,iz
      real(t_p) alsq2,alsq2n
      real(t_p) r2,pgamma,pgama1,damp
      real(t_p) f,e
      real(t_p) pdi,pti,pa,pb
      real(t_p) one
      real(t_p) pscale,dscale,uscale,fgrp
      integer iga,igb

      real(r_p),save:: vxx,vxy,vxz,vyy,vyz,vzz
      real(t_p) xix,xiy,xiz
      real(t_p) yix,yiy,yiz
      real(t_p) zix,ziy,ziz
      integer :: iax,iay,iaz
      logical, save :: f_in=.TRUE.

      type(rpole_elt):: ip,kp
      type(real6) :: dpui,dpuk
      type(real3) :: pos,posi,ufli,uflk,trqi,trqk,frc
      type(mdyn3_r):: frc_r
      parameter(one=1.0_ti_p)
      parameter(pgama1=0.0)

      interface
        subroutine epreal1c_correct_scale_group
     &   (trqvec,vxx,vxy,vxz,vyy,vyz,vzz)
          real(t_p),intent(inout):: trqvec(:,:)
          real(r_p),intent(inout):: vxx,vxy,vxz,vyy,vyz,vzz
        end subroutine epreal1c_correct_scale_group
      end interface

      if(deb_Path) then
!$acc wait
       write(*,'(2x,a)') '>>>epreal1_group'
      endif

      if (f_in) then
         f_in=.false.
         vxx = 0.0_re_p
          vxy = 0.0_re_p
          vxz = 0.0_re_p
          vyy = 0.0_re_p
          vyz = 0.0_re_p
          vzz = 0.0_re_p
!$acc enter data copyin(vxx,vxy,vxz,vyy,vyz,vzz)
      end if

      call prmem_request(fix,3,npolegroup, async=.true.)
      call prmem_request(fiy,3,npolegroup, async=.true.)
      call prmem_request(fiz,3,npolegroup, async=.true.)
      call prmem_request(trqvec,3,npolegroup, async=.true.)
      call set_to_zero1(trqvec,3*npolegroup,def_queue)

c
c     set conversion factor, cutoff and switching coefficients
c
      f = 0.5d0 * electric / dielec
c
c
c     compute the dipole polarization gradient components
c
!$acc parallel loop gang vector_length(32)
!$acc&         present(trqvec,vxx,vxy,vxz,vyy,vyz,vzz,epgroup)
!$acc&         present(poleglobgroup,ipolegroup,globglobgroup,
!$acc&  pollist,locgroup,polelocgroup,
!$acc&  rpole,thole,pdamp,
!$acc&  loc,x,y,z,uindgroup,uinpgroup,depgroup)
!$acc&         private(ip,dpui,posi)
!$acc&         async(def_queue)
      do ii = 1, npolelocgroup
         iipolegroup = poleglobgroup(ii)
         iglobgroup = ipolegroup(iipolegroup)
         iglob = globglobgroup(iglobgroup)
         iipole = pollist(iglob)
         ilocgroup = locgroup(iglobgroup)
         iploc = polelocgroup(iipolegroup)

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
!$acc loop vector private(frc,frc_r,trqi,trqk,dpuk,kp,pos)
         do kkk = iipolegroup+1, npolegroup
            kglobgroup = ipolegroup(kkk)
            kglob    = globglobgroup(kglobgroup)
            kpole    = pollist(kglob)
            klocgroup = locgroup(kglobgroup)
            kploc    = polelocgroup(kkk)

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
              frc%x=0.0;  frc%y=0.0;  frc%z=0.0;
            frc_r%x=0.0;frc_r%y=0.0;frc_r%z=0.0;
             trqi%x=0.0; trqi%y=0.0; trqi%z=0.0;
             trqk%x=0.0; trqk%y=0.0; trqk%z=0.0;
            call epolar1_couple(dpui,ip,dpuk,kp,r2,pos,
     &                  0.0_ti_p,0.0_ti_p,0.0_ti_p,pgamma,pgama1,damp,f,
     &                  1.0_ti_p,1.0_ti_p,1.0_ti_p,.false.,.false.
     &                  ,pa,pb,e,frc,frc_r,trqi,trqk,.false.)
c
c     increment energy
c
            epgroup  = epgroup + tp2enr(e)
c
c     increment gradient and virial due to Cartesian forces
c
            call atomic_sub( depgroup(1,ilocgroup),frc%x )
            call atomic_sub( depgroup(2,ilocgroup),frc%y )
            call atomic_sub( depgroup(3,ilocgroup),frc%z )
            call atomic_add( depgroup(1,klocgroup),frc%x )
            call atomic_add( depgroup(2,klocgroup),frc%y )
            call atomic_add( depgroup(3,klocgroup),frc%z )

            vxx      = vxx + pos%x*frc%x
            vxy      = vxy + pos%y*frc%x
            vxz      = vxz + pos%z*frc%x
            vyy      = vyy + pos%y*frc%y
            vyz      = vyz + pos%z*frc%y
            vzz      = vzz + pos%z*frc%z
c
c     increment torque
c
            call atomic_add( trqvec(1,iploc),trqi%x )
            call atomic_add( trqvec(2,iploc),trqi%y )
            call atomic_add( trqvec(3,iploc),trqi%z )
            call atomic_add( trqvec(1,kploc),trqk%x )
            call atomic_add( trqvec(2,kploc),trqk%y )
            call atomic_add( trqvec(3,kploc),trqk%z )
        enddo
      enddo

      call epreal1c_correct_scale_group(trqvec,vxx,vxy,vxz,vyy,vyz,vzz)

      call torquegpu_group(trqvec,fix,fiy,fiz,depgroup)

      if (use_virial) then

!$acc parallel loop async(def_queue)
!$acc&   present(globglobgroup,ipolegroup,pollist,polelocgroup
!$acc&  ,x,y,z,xaxis,yaxis,zaxis,fix,fiy,fiz)
        do i = 1, npolegroup
          iglob = globglobgroup(ipolegroup(i))
          iipole = pollist(iglob)
          ii = polelocgroup(i)
          iax    = xaxis(iipole)
          iay    = yaxis(iipole)
          iaz    = zaxis(iipole)
          if (iax.gt.0)  then
              xix = x(iax) - x(iglob) ! xix
              yix = y(iax) - y(iglob) ! yix
              zix = z(iax) - z(iglob) ! zix
          else
              xix = 0.0_ti_p ! xix
              yix = 0.0_ti_p ! yix
              zix = 0.0_ti_p ! zix
          endif
          if (iay.gt.0)  then
              xiy = x(iay) - x(iglob) ! xiy
              yiy = y(iay) - y(iglob) ! yiy
              ziy = z(iay) - z(iglob) ! ziy
          else
              xiy = 0.0_ti_p ! xiy
              yiy = 0.0_ti_p ! yiy
              ziy = 0.0_ti_p ! ziy
          endif
          if (iaz.gt.0) then
              xiz = x(iaz) - x(iglob) ! xiz
              yiz = y(iaz) - y(iglob) ! yiz
              ziz = z(iaz) - z(iglob) ! ziz
          else
              xiz = 0.0_ti_p ! xiz
              yiz = 0.0_ti_p ! yiz
              ziz = 0.0_ti_p ! ziz
          endif

          vxx = vxx + xix*fix(1,i) + xiy*fiy(1,i) + xiz*fiz(1,i)
          vxy = vxy + yix*fix(1,i) + yiy*fiy(1,i) + yiz*fiz(1,i)
          vxz = vxz + zix*fix(1,i) + ziy*fiy(1,i) + ziz*fiz(1,i)
          vyy = vyy + yix*fix(2,i) + yiy*fiy(2,i) + yiz*fiz(2,i)
          vyz = vyz + zix*fix(2,i) + ziy*fiy(2,i) + ziz*fiz(2,i)
          vzz = vzz + zix*fix(3,i) + ziy*fiy(3,i) + ziz*fiz(3,i)
        end do

!$acc serial present(vir_group,vxx,vxy,vxz,vyy,vyz,vzz) async
        vir_group(1,1) = vir_group(1,1) + vxx
        vir_group(2,1) = vir_group(2,1) + vxy
        vir_group(3,1) = vir_group(3,1) + vxz
        vir_group(1,2) = vir_group(1,2) + vxy
        vir_group(2,2) = vir_group(2,2) + vyy
        vir_group(3,2) = vir_group(3,2) + vyz
        vir_group(1,3) = vir_group(1,3) + vxz
        vir_group(2,3) = vir_group(2,3) + vyz
        vir_group(3,3) = vir_group(3,3) + vzz
        vxx = 0.0_re_p
        vxy = 0.0_re_p
        vxz = 0.0_re_p
        vyy = 0.0_re_p
        vyz = 0.0_re_p
        vzz = 0.0_re_p
!$acc end serial

      endif

      if(deb_Path) then
!$acc wait
       write(*,'(2x,a)') '<<<epreal1_group'
      endif

      return
      end


            ! Loop on scale interaction for correction
      subroutine epreal1c_correct_scale_group
     &           (trqvec,vxx,vxy,vxz,vyy,vyz,vzz)
      use atmlst  ,only: poleglobnl
      use atoms   ,only: x,y,z
      use chgpot  ,only: dielec,electric
      use deriv   ,only: dep
      use domdec  ,only: rank,loc
      use energi  ,only: ep=>ep_r
      use epolar1group_inl
      use ewald   ,only: aewald
      use inform  ,only: deb_Path
      use math    ,only: sqrtpi
      use mpole   ,only: npolelocnl,ipole,rpole,polelocnl
      use neigh   ,only: elst,nelst
      use polar   ,only: uind,uinp,thole,pdamp
      use shunt   ,only: off2
      use tinheader ,only: ti_p
      use tinTypes,only: real3,real6,mdyn3_r,rpole_elt
      use utilgpu ,only: def_queue
      use group
      use virial
      implicit none

      real(t_p),intent(inout):: trqvec(:,:)
      real(r_p),intent(inout):: vxx,vxy,vxz,vyy,vyz,vzz

      integer i,k,iglob,kglob,iploc,kploc
      integer ii,iipole,kpole,j,kbis
      integer iglobgroup,kglobgroup
      integer iipolegroup,kpolegroup
      integer nnelst
      real(t_p) alsq2,alsq2n
      real(t_p) r2,pgamma,pgama1,damp
      real(t_p) f,e,fgrp
      real(t_p) pdi,pti,pa,pb
      real(t_p) one
      real(t_p) pscale,dscale,uscale
      integer iga,igb
      type(rpole_elt):: ip,kp
      type(real6)    :: dpui,dpuk
      type(real3)    :: pos,posi,ufli,uflk,trqi,trqk,frc
      type(mdyn3_r)  :: frc_r

      parameter(one=1.0_ti_p)
      parameter(pgama1=0.0)

      if(deb_Path)
     &   write(*,'(2x,a)') 'epreal1c_correct_scale'
c
c     set conversion factor, cutoff and switching coefficients
c
      f      = 0.5_ti_p * electric / dielec

!$acc parallel loop present(trqvec,vxx,vxy,vxz,vyy,vyz,vzz)
!$acc&     present(loclocgroup,locgroup,polelocgroup,pollistgroup,
!$acc&  rpole,thole,pdamp,loc,x,y,z,
!$acc&  uindgroup,uinpgroup,depgroup,epgroup,vir,polelocnl
!$acc&  ,dpucorrect_ik_group,dpucorrect_scale_group)
!$acc&         private(pos,ip,kp,dpui,dpuk,trqk,trqi,frc_r,frc)
!$acc&         async(def_queue)
      do ii = 1, n_dpuscale_group
         iipole   = dpucorrect_ik_group(2*(ii-1)+1)
         kpole    = dpucorrect_ik_group(2*(ii-1)+2)

         iglob    = ipole(iipole)
         kglob    = ipole(kpole)

         iglobgroup = loclocgroup(iglob)
         iipolegroup= pollistgroup(iglobgroup)
         i        = locgroup(iglobgroup)
         iploc    = polelocgroup(iipolegroup)

         kglobgroup = loclocgroup(kglob)
         kpolegroup = pollistgroup(kglobgroup)
         k        = locgroup(kglobgroup)
         kploc    = polelocgroup(kpolegroup)

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

         dscale   = dpucorrect_scale_group(3*(ii-1)+1)
         pscale   = dpucorrect_scale_group(3*(ii-1)+2)
         uscale   = dpucorrect_scale_group(3*(ii-1)+3)

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
           frc%x=0.0;  frc%y=0.0;  frc%z=0.0;
         frc_r%x=0.0;frc_r%y=0.0;frc_r%z=0.0;
          trqi%x=0.0; trqi%y=0.0; trqi%z=0.0;
          trqk%x=0.0; trqk%y=0.0; trqk%z=0.0;
         call epolar1_couple(dpui,ip,dpuk,kp,r2,pos,
     &               0._ti_p,0._ti_p,0._ti_p,pgamma,pgama1,damp,f,
     &               dscale,pscale,uscale,.false.,.false.,pa,pb,
     &               e,frc,frc_r,trqi,trqk,.true.)
c
c     increment energy
c
         epgroup  = epgroup + tp2enr(e)
c
c     increment gradient and virial due to Cartesian forces
c
         call atomic_sub( depgroup(1,i),frc_r%x )
         call atomic_sub( depgroup(2,i),frc_r%y )
         call atomic_sub( depgroup(3,i),frc_r%z )
         call atomic_add( depgroup(1,k),frc_r%x )
         call atomic_add( depgroup(2,k),frc_r%y )
         call atomic_add( depgroup(3,k),frc_r%z )

         vxx      = vxx + pos%x*frc%x
         vxy      = vxy + pos%y*frc%x
         vxz      = vxz + pos%z*frc%x
         vyy      = vyy + pos%y*frc%y
         vyz      = vyz + pos%z*frc%y
         vzz      = vzz + pos%z*frc%z
c
c     increment torque
c
         call atomic_add( trqvec(1,iploc),trqi%x )
         call atomic_add( trqvec(2,iploc),trqi%y )
         call atomic_add( trqvec(3,iploc),trqi%z )
         call atomic_add( trqvec(1,kploc),trqk%x )
         call atomic_add( trqvec(2,kploc),trqk%y )
         call atomic_add( trqvec(3,kploc),trqk%z )
      end do
c
      end
