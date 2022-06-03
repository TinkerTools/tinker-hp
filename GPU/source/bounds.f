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

      subroutine bounds
      use sizes
      use tinheader
      use atmtyp
      use atomsMirror
      use bounds_inl
      use boxes ,box34_p=>box34
      use molcul
      implicit none
      integer i,j,k
      integer init,stop
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
         call imagem_inl(xfrac,yfrac,zfrac,xbox2,ybox2,zbox2,box34)
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
c
         do j = init, stop
            k = kmol(j)
            x(k) = x(k) - xmid + xcom
            y(k) = y(k) - ymid + ycom
            z(k) = z(k) - zmid + zcom
         end do
      end do
      call reCast_position

      call nblst_alter_bstat

      end
