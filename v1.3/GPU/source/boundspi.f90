!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     #################################################################
!     ##                                                             ##
!     ##  subroutine bounds  --  check periodic boundary conditions  ##
!     ##                                                             ##
!     #################################################################
!
!
!     "bounds" finds the center of mass of each molecule and
!     translates any stray molecules back into the periodic box
!     atoms involved in a CV (through colvars or plumed) are dealt with
!     separately
!
!
#include "tinker_precision.h"
module boundspi_inl
contains
#include "image.inc.f90"
end module

subroutine boundspi(polymer,ibead_beg,ibead_end)
   use sizes
   use atmtyp
   use atoms
   use boxes
   use molcul
   use beads
   use colvars
   use plumed
   implicit none
   type(POLYMER_COMM_TYPE), intent(inout) :: polymer
   integer, intent(in) :: ibead_beg,ibead_end
   integer i,j,k,l,lglob
   integer init,stop_
   integer nlist
   integer, allocatable :: list(:)

   if (use_colvars.or.lplumed) then
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
         call boundslistpi(polymer,ibead_beg,ibead_end,nlist,list)
      end if
!
!     then the other ones, molecule by molecule
!
      do i = 1, nmol
         nlist = 0
         list = 0
         init = imol(1,i)
         stop_ = imol(2,i)
         do j = init, stop_
            k = kmol(j)
            do l = 1, ncvatomsmol
               lglob = cvatomsmol(l)
               if (lglob.eq.k) then
                  goto 20
               end if
            end do
            nlist = nlist + 1
            list(nlist) = k
20          continue
         end do
         call boundslistpi(polymer,ibead_beg,ibead_end,nlist,list)
      end do

   else if (lplumed) then
!
!      with PLUMED we don't know which atoms are in the cv so we wrap
!      the whole system at once
!
      nlist = n
      do i = 1, n
         list(i) = i
      end do
      call boundslistpi(polymer,ibead_beg,ibead_end,nlist,list)

   else
!
!     wrap by molecules
!
      call bounds_by_mol_pi(polymer,ibead_beg,ibead_end)

   end if

   if (use_colvars.or.lplumed) then
      deallocate (list)
   end if

!$ acc update device(pbcwrap,pbcWrapIdx,polymer%eigpos(:,:,1)) async

   call nblst_alter_bstat

   if (use_colvars.or.pbcunwrap.or.lplumed) then
      call chk_unwrap_status
   end if

   return
end
!
!     "boundslist" finds the center of mass of a list of atoms
!     translates any stray atoms back into the periodic box
!
!
subroutine boundslistpi(polymer,ibead_beg,ibead_end,nlist,list)
   use sizes
   use atmtyp
   use atoms
   use boxes
   use beads
   use boundspi_inl
   use tinheader, only: re_p
   use boxes ,box34_p=>box34
   implicit none
   type(POLYMER_COMM_TYPE), intent(inout) :: polymer
   integer, intent(in) :: ibead_beg,ibead_end
   integer, intent(in) :: nlist,list(nlist)
   integer i,j,k,ibead,kk
   real*8 weigh,weight
   real*8 xmid,ymid,zmid
   real*8 dx,dy,dz
   integer(1) wx,wy,wz
   real(r_p) xcom,ycom,zcom
   real(r_p) xfrac,yfrac,zfrac
   real(r_p) box34

   box34 = 0.75_re_p*xbox

!
!
   xmid = 0.0d0
   ymid = 0.0d0
   zmid = 0.0d0
   weight = 0.0d0
   do j = 1, nlist
      k = list(j)
      weigh = mass(k)
      weight = weight + weigh
      xmid = xmid + polymer%eigpos(1,k,1)*weigh
      ymid = ymid + polymer%eigpos(2,k,1)*weigh
      zmid = zmid + polymer%eigpos(3,k,1)*weigh
   end do
   xmid = xmid / weight
   ymid = ymid / weight
   zmid = zmid / weight

!
!     get fractional coordinates of center of mass
!
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
!
!     translate center of mass into the periodic box
!
   call imagem_inl(xfrac,wx,yfrac,wy,zfrac,wz,xbox2,ybox2,zbox2&
   &,box34)
!
!     convert translated fraction center of mass to Cartesian
!
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
!
!     translate coordinates via offset from center of mass
!
   do j = 1, nlist
      k = list(j)
      kk   = 4*(k-1)
      polymer%eigpos(1,k,1) = polymer%eigpos(1,k,1) - xmid + xcom
      polymer%eigpos(2,k,1) = polymer%eigpos(2,k,1) - ymid + ycom
      polymer%eigpos(3,k,1) = polymer%eigpos(3,k,1) - zmid + zcom
      pbcWrapIdx(kk+1) = pbcWrapIdx(kk+1) + wx
      pbcWrapIdx(kk+2) = pbcWrapIdx(kk+2) + wy
      pbcWrapIdx(kk+3) = pbcWrapIdx(kk+3) + wz
   end do
   do ibead = ibead_beg,ibead_end
      do j=1,nlist
         k=list(j)
         polymer%pos(1,k,ibead) = polymer%pos(1,k,ibead) - xmid + xcom
         polymer%pos(2,k,ibead) = polymer%pos(2,k,ibead) - ymid + ycom
         polymer%pos(3,k,ibead) = polymer%pos(3,k,ibead) - zmid + zcom
      enddo
   enddo
   return
end

 ! Wrap each molecules individually (OpenACC)
subroutine bounds_by_mol_pi(polymer,ibead_beg,ibead_end)
   use sizes
   use beads
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
   type(POLYMER_COMM_TYPE), intent(inout) :: polymer
   integer, intent(in) :: ibead_beg,ibead_end
   integer i,j,k,kk,ibead
   integer init,stop_
   integer(1) wx,wy,wz
   real(r_p) weigh
   real(r_p) xmid,ymid,zmid
   real(r_p) xfrac,yfrac,zfrac
   real(r_p) xcom,ycom,zcom
   real(r_p) box34

   box34 = 0.75_re_p*xbox
!
!     locate the center of mass of each molecule
!

!!$acc update device(polymer%eigpos(:,:,1) &
!!$acc   ,polymer%pos(:,:,ibead_beg:ibead_end)) async
!!$acc parallel loop gang worker async default(present)

   do i = 1, nmol
      init = imol(1,i)
      stop_ = imol(2,i)
      xmid = 0.0_re_p
      ymid = 0.0_re_p
      zmid = 0.0_re_p
      do j = init, stop_
         k = kmol(j)
         weigh = mass(k)
         xmid = xmid + polymer%eigpos(1,k,1)*weigh
         ymid = ymid + polymer%eigpos(2,k,1)*weigh
         zmid = zmid + polymer%eigpos(3,k,1)*weigh
      end do
      weigh = molmass(i)
      xmid = xmid / weigh
      ymid = ymid / weigh
      zmid = zmid / weigh
!
!     get fractional coordinates of center of mass
!
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
!
!     translate center of mass into the periodic box
!
      call imagem_inl(xfrac,wx,yfrac,wy,zfrac,wz,xbox2,ybox2,zbox2&
      &,box34)
!
!     convert translated fraction center of mass to Cartesian
!
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
!
!     translate coordinates via offset from center of mass
!     &
!     store the index of the shift in pbc cell
!
      do j = init, stop_
         k    = kmol(j)
         kk   = 4*(k-1)
         polymer%eigpos(1,k,1) = polymer%eigpos(1,k,1)-xmid+xcom
         polymer%eigpos(2,k,1) = polymer%eigpos(2,k,1)-ymid+ycom
         polymer%eigpos(3,k,1) = polymer%eigpos(3,k,1)-zmid+zcom
         pbcWrapIdx(kk+1) = pbcWrapIdx(kk+1) + wx
         pbcWrapIdx(kk+2) = pbcWrapIdx(kk+2) + wy
         pbcWrapIdx(kk+3) = pbcWrapIdx(kk+3) + wz
!        if (wx.ne.0) print*,'x', k,int(wx),int(pbcWrapIdx(kk+1))
!        if (wy.ne.0) print*,'y', k,int(wy),int(pbcWrapIdx(kk+2))
!        if (wz.ne.0) print*,'z', k,int(wz),int(pbcWrapIdx(kk+3))

!!$acc loop seq
         do ibead = ibead_beg,ibead_end
            polymer%pos(1,k,ibead) = polymer%pos(1,k,ibead)-xmid+xcom
            polymer%pos(2,k,ibead) = polymer%pos(2,k,ibead)-ymid+ycom
            polymer%pos(3,k,ibead) = polymer%pos(3,k,ibead)-zmid+zcom
         enddo
      end do

   end do
!!$acc update host(polymer%eigpos(:,:,1),pbcWrapIdx,pbcwrap
!!$acc&   ,polymer%pos(:,:,ibead_beg:ibead_end)) async
end
