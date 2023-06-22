c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine bounds  --  check periodic boundary conditions  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "bounds" finds the center of mass of each molecule and
c     translates any stray molecules back into the periodic box
c
c
#include "tinker_precision.h"
      module bounds_inl
        contains
#include "image.f.inc"
      end module

      subroutine chk_unwrap_status
      use atoms
      use inform
      implicit none
      integer i,j
      integer cpt
      integer(1),parameter:: char1m=126

      cpt = 0

!$acc parallel loop async default(present) reduction(+:cpt)
      do i = 1,n
         j = 4*(i-1)
         if (abs(pbcWrapIdx(j+1)).gt.char1m) cpt=cpt+1
         if (abs(pbcWrapIdx(j+2)).gt.char1m) cpt=cpt+1
         if (abs(pbcWrapIdx(j+3)).gt.char1m) cpt=cpt+1
      end do
!$acc wait

      if (cpt.gt.0) then
 66   format(" WARNING !! Wrapping State Is reaching full capacity !!"
     &    ,/,"  -> Dynamic will interrupt to preserve the sampling"
     &      ," Quality"
     &    ,/,"  -> Please Restart your sampling with your latest"
     &      ," restart")
         write(0,66) 
         abort = .true.
      end if
      end subroutine
c
c     #####################################################################
c     ##                                                                 ##
c     ##  subroutine boundslist  --  check periodic boundary conditions  ##
c     ##                                                                 ##
c     #####################################################################
c
c
c     "boundslist" finds the center of mass of a group of atoms and
c     translates its center of mass back into the periodic box
c
      subroutine boundslist(nlist,list)
      use sizes
      use tinheader
      use atmtyp
      use atoms       ,only: pbcWrapIdx,pbcunwrap
      use atomsMirror
      use bounds_inl
      use boxes ,box34_p=>box34
      integer i,j,k
      integer nlist,list(nlist)
      integer(1) wx,wy,wz
      real(r_p) weigh,weight
      real(r_p) xmid,ymid,zmid
      real(r_p) xfrac,yfrac,zfrac
      real(r_p) xcom,ycom,zcom
      real(r_p) box34
c
      box34 = 0.75_re_p*xbox
c
c     locate the center of mass of each molecule
c
      xmid = 0.0_re_p
      ymid = 0.0_re_p
      zmid = 0.0_re_p
      weigh = 0.0_re_p
      weight = 0.0_re_p
      do j = 1, nlist
         k = list(j)
         weigh = mass(k)
         weight = weight + weigh
         xmid = xmid + x(k)*weigh
         ymid = ymid + y(k)*weigh
         zmid = zmid + z(k)*weigh
      end do
      xmid = xmid / weight
      ymid = ymid / weight
      zmid = zmid / weight
c
c     get fractional coordinates of center of mass
c
      if (orthogonal .or. octahedron) then
         zfrac = zmid
         yfrac = ymid
         xfrac = xmid
      else if (monoclinic) then
         zfrac = zmid / beta_sin
         yfrac = ymid
         xfrac = xmid - zfrac*beta_cos
      else if (triclinic) then
         zfrac = zmid / gamma_term
         yfrac = (ymid - zfrac*beta_term) / gamma_sin
         xfrac = xmid - yfrac*gamma_cos - zfrac*beta_cos
      end if
c
c     translate center of mass into the periodic box
c
      call imagem_inl(xfrac,wx,yfrac,wy,zfrac,wz,xbox2,ybox2,zbox2
     &               ,box34)
c
c     convert translated fraction center of mass to Cartesian
c
      if (orthogonal .or. octahedron) then
         xcom = xfrac
         ycom = yfrac
         zcom = zfrac
      else if (monoclinic) then
         xcom = xfrac + zfrac*beta_cos
         ycom = yfrac
         zcom = zfrac * beta_sin
      else if (triclinic) then
         xcom = xfrac + yfrac*gamma_cos + zfrac*beta_cos
         ycom = yfrac*gamma_sin + zfrac*beta_term
         zcom = zfrac * gamma_term
      end if
c
c     translate coordinates via offset from center of mass
c     &
c     store the index of the shift in pbc cell
c
      do j = 1, nlist
         k    = list(j)
         kk   = 4*(k-1)
         x(k) = x(k) - xmid + xcom
         y(k) = y(k) - ymid + ycom
         z(k) = z(k) - zmid + zcom
         pbcWrapIdx(kk+1) = pbcWrapIdx(kk+1) + wx
         pbcWrapIdx(kk+2) = pbcWrapIdx(kk+2) + wy
         pbcWrapIdx(kk+3) = pbcWrapIdx(kk+3) + wz
      end do
      end

      ! Wrap each molecules individually (OpenACC)
      subroutine bounds_by_mol
      use sizes
      use tinheader
      use atmtyp
      use atoms       ,only: pbcWrapIdx,pbcunwrap
      use atomsMirror
      use bounds_inl
      use boxes ,box34_p=>box34
      use colvars     ,only: use_colvars
      use plumed      ,only: lplumed
      use molcul
      implicit none
      integer i,j,k,kk
      integer init,stop
      integer(1) wx,wy,wz
      real(r_p) weigh
      real(r_p) xmid,ymid,zmid
      real(r_p) xfrac,yfrac,zfrac
      real(r_p) xcom,ycom,zcom
      real(r_p) box34

      box34 = 0.75_re_p*xbox
c
c     locate the center of mass of each molecule
c
!$acc parallel loop gang worker async default(present)
      do i = 1, nmol
         init = imol(1,i)
         stop = imol(2,i)
         xmid = 0.0_re_p
         ymid = 0.0_re_p
         zmid = 0.0_re_p
         do j = init, stop
            k = kmol(j)
            weigh = mass(k)
            xmid = xmid + x(k)*weigh
            ymid = ymid + y(k)*weigh
            zmid = zmid + z(k)*weigh
         end do
         weigh = molmass(i)
         xmid = xmid / weigh
         ymid = ymid / weigh
         zmid = zmid / weigh
c
c     get fractional coordinates of center of mass
c
         if (orthogonal .or. octahedron) then
            zfrac = zmid
            yfrac = ymid
            xfrac = xmid
         else if (monoclinic) then
            zfrac = zmid / beta_sin
            yfrac = ymid
            xfrac = xmid - zfrac*beta_cos
         else if (triclinic) then
            zfrac = zmid / gamma_term
            yfrac = (ymid - zfrac*beta_term) / gamma_sin
            xfrac = xmid - yfrac*gamma_cos - zfrac*beta_cos
         end if
c
c     translate center of mass into the periodic box
c
         call imagem_inl(xfrac,wx,yfrac,wy,zfrac,wz,xbox2,ybox2,zbox2
     &                  ,box34)
c
c     convert translated fraction center of mass to Cartesian
c
         if (orthogonal .or. octahedron) then
            xcom = xfrac
            ycom = yfrac
            zcom = zfrac
         else if (monoclinic) then
            xcom = xfrac + zfrac*beta_cos
            ycom = yfrac
            zcom = zfrac * beta_sin
         else if (triclinic) then
            xcom = xfrac + yfrac*gamma_cos + zfrac*beta_cos
            ycom = yfrac*gamma_sin + zfrac*beta_term
            zcom = zfrac * gamma_term
         end if
c
c     translate coordinates via offset from center of mass
c     &
c     store the index of the shift in pbc cell
c
         do j = init, stop
            k    = kmol(j)
            kk   = 4*(k-1)
            x(k) = x(k) - xmid + xcom
            y(k) = y(k) - ymid + ycom
            z(k) = z(k) - zmid + zcom
            pbcWrapIdx(kk+1) = pbcWrapIdx(kk+1) + wx
            pbcWrapIdx(kk+2) = pbcWrapIdx(kk+2) + wy
            pbcWrapIdx(kk+3) = pbcWrapIdx(kk+3) + wz
c        if (wx.ne.0) print*,'x', k,int(wx),int(pbcWrapIdx(kk+1))
c        if (wy.ne.0) print*,'y', k,int(wy),int(pbcWrapIdx(kk+2))
c        if (wz.ne.0) print*,'z', k,int(wz),int(pbcWrapIdx(kk+3))
         end do
      end do
      end

      subroutine bounds
      use sizes
      use tinheader
      use atmtyp
      use atoms       ,only: pbcWrapIdx,pbcunwrap
      use atomsMirror
      use bounds_inl
      use boxes ,box34_p=>box34
      use colvars     ,only: use_colvars,ncvatomsmol,cvatomsmol
      use plumed      ,only: lplumed
      use molcul
      implicit none
      integer i,j,k,l,lglob
      integer init,stop
      integer nlist
      integer, allocatable :: list(:)

      if (use_colvars.or.lplumed) then
!$acc wait
!$acc update host(x,y,z,pbcWrapIdx)
         if (allocated(list)) deallocate(list)
         allocate (list(n))
      end if

      if (use_colvars) then

        !TODO Offload onto device
        if (ncvatomsmol.gt.0) then
          nlist = ncvatomsmol
          do j = 1, ncvatomsmol
            list(j) = cvatomsmol(j)
          end do
          call boundslist(nlist,list)
        end if
c
c     then the other ones, molecule by molecule
c
        do i = 1, nmol
          nlist = 0
          list = 0
          init = imol(1,i)
          stop = imol(2,i)
          do j = init, stop
            k = kmol(j)
            do l = 1, ncvatomsmol
              lglob = cvatomsmol(l)
              if (lglob.eq.k) then
                goto 20
              end if
            end do
            nlist = nlist + 1
            list(nlist) = k
 20         continue
          end do
          call boundslist(nlist,list)
        end do

      else if (lplumed) then
c
c      with PLUMED we don't know which atoms are in the cv so we wrap
c      the whole system at once
c
        !TODO Offload onto device
        nlist = n
        do i = 1, n
          list(i) = i
        end do
        call boundslist(nlist,list)

      else
c
c     wrap by molecules
c
        call bounds_by_mol

      end if

      if (use_colvars.or.lplumed) then
!$acc update device(x,y,z,pbcWrapIdx) async
         deallocate (list)
      end if
c
      call reCast_position

      call nblst_alter_bstat

      if (use_colvars.or.pbcunwrap.or.lplumed) then
         call chk_unwrap_status
      end if
      end
