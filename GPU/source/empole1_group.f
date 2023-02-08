c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ###################################################################
c     ##                                                               ##
c     ##  subroutine empole1_group  --  "group" multipoles             ##
c     ##                                                               ##
c     ###################################################################
c
c
c     "empole1_group" calculates the permanent multipoles energy and
c     derivatives with respect to Cartesian coordinates using
c     a subset of atoms in gaz phase
c
c
#include "tinker_precision.h"
#include "tinker_types.h"
       module empole1_group_inl
        use tinTypes , only: real3,real3_red,rpole_elt
#include "atomicOp.h.f"
        contains
#include "convert.f.inc"
#include "image.f.inc"
#include "atomicOp.inc.f"
#include "switch_respa.f.inc"
#include "pair_mpole1.f.inc"
      end module


      subroutine empole1_group
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
!$acc serial async present(emgroup)
      emgroup = 0.0d0
!$acc end serial

      call prmem_requestm(demgroup,3,natgroup, async=.true.)
!$acc parallel loop collapse(2) default(present) async
      do i = 1, natgroup
         do j = 1, 3
            demgroup(j,i) = 0
         end do
      end do
c
      if (npolegroup .eq. 0)  return

c
      call emreal1_group
c
c     communicate the forces
c
      call commforcesgroup

      return
      end
      
c     #################################################################
c     ##                                                             ##
c     ##  subroutine emreal1_group  --  group multipoles derivs      ##
c     ##                                                             ##
c     #################################################################
c
c
c     "emreal1_group" evaluates the permanent
c     multipoles energy and gradient
c     of a subset of atoms
c
c
      subroutine emreal1_group
      use atmlst
      use atoms
      use bound
      use chgpot
      use couple
      use cutoff
      use deriv
      use domdec
      use energi
      use elec_wspace,only: fix=>r2Work1,fiy=>r2Work2,fiz=>r2Work3
     &               ,trqvec=>r2Work4
      !use erf_mod
      use empole1_group_inl
      use ewald
      use inform     ,only: deb_Path
      use inter
      use iounit
      use interfaces ,only: torquegpu_group
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
      integer iax,iay,iaz,ver,fea
      integer(1) mutik
      logical use_grp,use_lamdyn,use_chgf
      real(t_p) alsq2,alsq2n,fgrp
      real(t_p) r2,pgamma,damp
      real(t_p) f,e
      real(t_p) pdi,pti,poti,potk,delambdae_
      real(t_p) one
      real(t_p) xix,xiy,xiz
      real(t_p) yix,yiy,yiz
      real(t_p) zix,ziy,ziz
      type(rpole_elt) ip,kp
      type(real3)     pos,posi,ufli,uflk,trqi,trqk,frc
      type(mdyn3_r)   frc_r
      real(r_p) ,save:: vxx,vxy,vxz,vyy,vyz,vzz
      logical   ,save:: f_in=.TRUE.
      parameter (one=1.0_ti_p,alsq2=0.0,alsq2n=0.0
     &   , use_grp=.false.,use_chgf=.false.
     &   , fgrp=1.0
     &   , use_lamdyn=.false.
     &   , mutik=0
     &   , ver=__use_grd__+__use_ene__+__use_vir__
     &   , fea=__use_mpi__
     &)

      interface
        subroutine emreal1_correct_scale_group
     &   (trqvec,vxx,vxy,vxz,vyy,vyz,vzz)
          real(t_p),intent(inout):: trqvec(:,:)
          real(r_p),intent(inout):: vxx,vxy,vxz,vyy,vyz,vzz
        end subroutine emreal1_correct_scale_group
      end interface

      if(deb_Path) then
!$acc wait
       write(*,'(2x,a)') '>>>emreal1_group'
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
      f = electric / dielec
c
c
c     compute the dipole polarization gradient components
c
!$acc parallel loop gang vector_length(32)
!$acc&         present(trqvec,vxx,vxy,vxz,vyy,vyz,vzz,emgroup)
!$acc&         present(poleglobgroup,ipolegroup,globglobgroup,
!$acc&  pollist,locgroup,polelocgroup,rpole,
!$acc&  loc,x,y,z,demgroup)
!$acc&         private(ip,posi)
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
c
!$acc loop vector private(frc,trqi,trqk,kp,pos)
         do kkk = iipolegroup+1, npolegroup
            kglobgroup = ipolegroup(kkk)
            kglob = globglobgroup(kglobgroup)
            kpole = pollist(kglob)
            klocgroup = locgroup(kglobgroup)
            kploc = polelocgroup(kkk)

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
c
            ! compute mpole one interaction
            call duo_mpole(r2,pos%x,pos%y,pos%z,ip,kp,0.0,ewaldshortcut
     &              ,shortheal,0.0,f,0.0,0.0,use_grp,fgrp
     &              ,use_lamdyn,mutik,0.0,use_chgf
     &              ,poti,potk,delambdae_,e,frc,trqi,trqk,ver,fea)
c
c     increment energy
c
            emgroup          = emgroup + tp2enr(e)
c
c     increment gradient and virial due to Cartesian forces
c
            call atomic_add( demgroup(1,ilocgroup),frc%x )
            call atomic_add( demgroup(2,ilocgroup),frc%y )
            call atomic_add( demgroup(3,ilocgroup),frc%z )
            call atomic_sub( demgroup(1,klocgroup),frc%x )
            call atomic_sub( demgroup(2,klocgroup),frc%y )
            call atomic_sub( demgroup(3,klocgroup),frc%z )
c
            vxx         = vxx - pos%x * frc%x
            vxy         = vxy - pos%y * frc%x
            vxz         = vxz - pos%z * frc%x
            vyy         = vyy - pos%y * frc%y
            vyz         = vyz - pos%z * frc%y
            vzz         = vzz - pos%z * frc%z
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

      call emreal1_correct_scale_group(trqvec,vxx,vxy,vxz,vyy,vyz,vzz)

      call torquegpu_group(trqvec,fix,fiy,fiz,demgroup)

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
       write(*,'(2x,a)') '<<<emreal1_group'
      endif

      return
      end


            ! Loop on scale interaction for correction
      subroutine emreal1_correct_scale_group
     &           (trqvec,vxx,vxy,vxz,vyy,vyz,vzz)
      use atmlst  ,only: poleglobnl
      use atoms   ,only: x,y,z
      use chgpot  ,only: dielec,electric
      use cutoff  ,only: shortheal,ewaldshortcut
      use domdec  ,only: rank,loc
      use empole1_group_inl
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
      integer nnelst,ver,fea
      integer(1) mutik
      logical use_grp,use_chgf,use_lamdyn
      real(t_p) alsq2,alsq2n
      real(t_p) r2,pgamma,damp
      real(t_p) f,e,fgrp
      real(t_p) pdi,pti,poti,potk,delambdae_
      real(t_p) one
      real(t_p) mscale
      integer iga,igb
      type(rpole_elt):: ip,kp
      type(real6)    :: dpui,dpuk
      type(real3)    :: pos,posi,ufli,uflk,trqi,trqk,frc
      type(mdyn3_r)  :: frc_r

      parameter (one=1.0_ti_p,alsq2=0.0,alsq2n=0.0
     &   , use_grp=.false.,use_chgf=.false.,use_lamdyn=.false.
     &   , fgrp=1.0
     &   , mutik=0
     &   , ver=__use_grd__+__use_ene__+__use_vir__+__use_sca__
     &   , fea=__use_mpi__
     &)

      if(deb_Path)
     &   write(*,'(2x,a)') 'emreal1_correct_scale'
c
c     set conversion factor, cutoff and switching coefficients
c
      f      = electric / dielec

!$acc parallel loop present(trqvec,vxx,vxy,vxz,vyy,vyz,vzz)
!$acc&     present(loclocgroup,locgroup,polelocgroup,pollistgroup,
!$acc&  rpole,thole,pdamp,loc,x,y,z,
!$acc&  demgroup,emgroup,vir,polelocnl,
!$acc&  mcorrect_ik_group,mcorrect_scale_group)
!$acc&         private(pos,ip,kp,trqk,trqi,frc_r,frc)
!$acc&         async(def_queue)
      do ii = 1, n_mscale_group
         iipole = mcorrect_ik_group(ii,1)
         kpole = mcorrect_ik_group(ii,2)
         mscale = mcorrect_scale_group(ii)


         iglob    = ipole(iipole)
         kglob    = ipole(kpole)

         iglobgroup = loclocgroup(iglob)
         i = locgroup(iglobgroup)
         iipolegroup = pollistgroup(iglobgroup)
         iploc = polelocgroup(iipolegroup)

         kglobgroup = loclocgroup(kglob)
         k = locgroup(kglobgroup)
         kpolegroup = pollistgroup(kglobgroup)
         kploc = polelocgroup(kpolegroup)

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
c
         ! compute mpole one interaction
         call duo_mpole(r2,pos%x,pos%y,pos%z,ip,kp,mscale,ewaldshortcut
     &           ,shortheal,0.0,f,0.0,0.0,use_grp,fgrp
     &           ,use_lamdyn,mutik,0.0,use_chgf
     &           ,poti,potk,delambdae_,e,frc,trqi,trqk,ver,fea)
c
c     increment energy
c
         emgroup       = emgroup + tp2enr(e)
c
c     increment gradient and virial due to Cartesian forces
c
         call atomic_add( demgroup(1,i),frc%x )
         call atomic_add( demgroup(2,i),frc%y )
         call atomic_add( demgroup(3,i),frc%z )
         call atomic_sub( demgroup(1,k),frc%x )
         call atomic_sub( demgroup(2,k),frc%y )
         call atomic_sub( demgroup(3,k),frc%z )
c
         vxx      = vxx - pos%x*frc%x
         vxy      = vxy - pos%y*frc%x
         vxz      = vxz - pos%z*frc%x
         vyy      = vyy - pos%y*frc%y
         vyz      = vyz - pos%z*frc%y
         vzz      = vzz - pos%z*frc%z
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
      end
