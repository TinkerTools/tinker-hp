c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine lattice  --  setup periodic boundary conditions  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "lattice" stores the periodic box dimensions and sets angle
c     values to be used in computing fractional coordinates
c
c
#include "tinker_precision.h"
      subroutine lattice
      use boxes
      use cell
      use math
#ifdef _OPENACC
      use interfaces ,only: C_get_cell
      use utilcu,only: copy_cell_cu
      use pmestuffcu,only: copy_recip_cu
#endif
      use tinheader
      implicit none
      integer   octa
      real(r_p) alpha_cos
      real(r_p) ar1,ar2,ar3
      real(r_p) br1,br2,br3
      real(r_p) cr1,cr2,cr3
c
c     compute and store the half box length values
c
      xbox2 = 0.5_re_p * xbox
      ybox2 = 0.5_re_p * ybox
      zbox2 = 0.5_re_p * zbox
      if (octahedron)  box34 = 0.75_re_p * xbox
      octa = merge(1,0,octahedron)

#ifdef ORTHOGONAL_BOX_SHAPE_ONLY
      ! Check box shape
      if (octahedron.or.triclinic.or.monoclinic) then
         print 44
 44      format( "FATAL ERROR !!! Unsuitable Box shape for this",
     &           " configuration")
         call fatal
      end if
#endif
c
c     set replicated cell dimensions equal to the unitcell
c
      xcell   = xbox
      ycell   = ybox
      zcell   = zbox
      xcell2  = xbox2
      ycell2  = ybox2
      zcell2  = zbox2
      i_xcell = real(1.0/real(xcell,8),t_p)
      i_ycell = real(1.0/real(ycell,8),t_p)
      i_zcell = real(1.0/real(zcell,8),t_p)
#if defined(SINGLE)||defined(MIXED)
      eps_cell = 10*max(xcell,ycell,zcell)*prec1_eps
#else
      eps_cell = 1d3*max(xcell,ycell,zcell)*prec1_eps
#endif
!$acc wait
!$acc update device(xbox2,ybox2,zbox2,xcell,ycell,zcell,xcell2,
!$acc& ycell2,zcell2,i_xcell,i_ycell,i_zcell,eps_cell,box34)
#ifdef _OPENACC
      call copy_cell_cu(xcell,ycell,zcell,xcell2,ycell2,zcell2,
     &     eps_cell,octahedron,box34)
      call C_get_cell(xcell,ycell,zcell,eps_cell,octa,box34)
#endif
c
c     get values needed for fractional coordinate computations
c
      if (orthogonal .or. octahedron) then
         alpha_cos = 0.0_ti_p
         beta_sin = 1.0_ti_p
         beta_cos = 0.0_ti_p
         gamma_sin = 1.0_ti_p
         gamma_cos = 0.0_ti_p
         beta_term = 0.0_ti_p
         gamma_term = 1.0_ti_p
      else if (monoclinic) then
         alpha_cos = 0.0_ti_p
         beta_sin = sin(beta/radian)
         beta_cos = cos(beta/radian)
         gamma_sin = 1.0_ti_p
         gamma_cos = 0.0_ti_p
         beta_term = 0.0_ti_p
         gamma_term = beta_sin
      else if (triclinic) then
         alpha_cos = cos(alpha/radian)
         beta_sin = sin(beta/radian)
         beta_cos = cos(beta/radian)
         gamma_sin = sin(gamma/radian)
         gamma_cos = cos(gamma/radian)
         beta_term = (alpha_cos - beta_cos*gamma_cos) / gamma_sin
         gamma_term = sqrt(beta_sin**2 - beta_term**2)
      end if
c
c     determine the volume of the parent periodic box
c
      volbox = 0.0_re_p
      if (orthogonal .or. octahedron) then
         volbox = xbox * ybox * zbox
      else if (monoclinic) then
         volbox = beta_sin * xbox * ybox * zbox
      else if (triclinic) then
         volbox = (gamma_sin*gamma_term) * xbox * ybox * zbox
      end if
c
c     compute and store real space lattice vectors as rows
c
      ar1 = xbox
      ar2 = 0.0_ti_p
      ar3 = 0.0_ti_p
      br1 = ybox * gamma_cos
      br2 = ybox * gamma_sin
      br3 = 0.0_ti_p
      cr1 = zbox * beta_cos
      cr2 = zbox * beta_term
      cr3 = zbox * gamma_term
      lvec(1,1) = ar1
      lvec(1,2) = ar2
      lvec(1,3) = ar3
      lvec(2,1) = br1
      lvec(2,2) = br2
      lvec(2,3) = br3
      lvec(3,1) = cr1
      lvec(3,2) = cr2
      lvec(3,3) = cr3
!$acc update device(lvec(:,:))
c
c     compute and store reciprocal lattice vectors as columns
c
      if (volbox .ne. 0.0_ti_p) then
         recip(1,1) = (br2*cr3 - cr2*br3) / volbox
         recip(2,1) = (br3*cr1 - cr3*br1) / volbox
         recip(3,1) = (br1*cr2 - cr1*br2) / volbox
         recip(1,2) = (cr2*ar3 - ar2*cr3) / volbox
         recip(2,2) = (cr3*ar1 - ar3*cr1) / volbox
         recip(3,2) = (cr1*ar2 - ar1*cr2) / volbox
         recip(1,3) = (ar2*br3 - br2*ar3) / volbox
         recip(2,3) = (ar3*br1 - br3*ar1) / volbox
         recip(3,3) = (ar1*br2 - br1*ar2) / volbox
#ifdef _CUDA
         call copy_recip_cu
#endif
!$acc update device(recip(:,:))

        !get the fractional to Cartesian transformation matrix
         call frac_to_cartgpu

        !find the matrix to convert Cartesian to fractional
         call cart_to_fracgpu
      end if
c
c     volume of truncated octahedron is half of cubic parent
c
      if (octahedron)  volbox = 0.5_re_p * volbox
!$acc update device(volbox)
      end
