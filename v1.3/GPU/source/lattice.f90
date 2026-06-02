!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ##################################################################
!     ##                                                              ##
!     ##  subroutine lattice  --  setup periodic boundary conditions  ##
!     ##                                                              ##
!     ##################################################################
!
!
!     "lattice" stores the periodic box dimensions and sets angle
!     values to be used in computing fractional coordinates
!
!
#include "tinker_precision.h"
subroutine lattice
   use boxes
   use cell
   use domdec     ,only: ranktot
   use math
#ifdef _OPENACC
   use interfaces ,only: C_get_cell
   use utilcu     ,only: copy_cell_cu
   use pmestuffcu ,only: copy_recip_cu
#endif
   use tinheader
   use sizes      ,only: tinkerdebug
   implicit none
   integer   mono,tric,octa
   real(r_p) alpha_cos
   real(r_p) ar1,ar2,ar3
   real(r_p) br1,br2,br3
   real(r_p) cr1,cr2,cr3
!
!     compute and store the half box length values
!
   xbox2 = 0.5_re_p * xbox
   ybox2 = 0.5_re_p * ybox
   zbox2 = 0.5_re_p * zbox
   if (octahedron)  box34 = 0.75_re_p * xbox
   mono = merge(1,0,monoclinic)
   tric = merge(1,0,triclinic)
   octa = merge(1,0,octahedron)

#ifdef ORTHOGONAL_BOX_SHAPE_ONLY
   ! Check box shape
   if (octahedron.or.triclinic.or.monoclinic) then
      print 44
44    format( "FATAL ERROR !!! Unsuitable Box shape for this",&
         &" configuration")
      call fatal
   end if
#endif
!
!     set replicated cell dimensions equal to the unitcell
!
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
!
!     get values needed for fractional coordinate computations
!
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
!
!     determine the volume of the parent periodic box
!
   volbox = 0.0_re_p
   if (orthogonal .or. octahedron) then
      volbox = xbox * ybox * zbox
   else if (monoclinic) then
      volbox = beta_sin * xbox * ybox * zbox
   else if (triclinic) then
      volbox = (gamma_sin*gamma_term) * xbox * ybox * zbox
   end if
!
!     compute and store real space lattice vectors as rows
!
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
!
!   Transfer To Device and C environment
!
!$acc wait
!$acc update device(xbox2,ybox2,zbox2,xcell,ycell,zcell,xcell2,ycell2,zcell2&
!$acc    ,beta_sin,beta_cos,beta_term,gamma_sin,gamma_cos,gamma_term&
!$acc    ,i_xcell,i_ycell,i_zcell,eps_cell,box34 )
#ifdef _OPENACC
   call copy_cell_cu(xcell,ycell,zcell,xcell2,ycell2,zcell2,eps_cell&
            ,beta_sin,beta_cos,beta_term,gamma_sin,gamma_cos,gamma_term&
            ,monoclinic,triclinic,octahedron,box34)
   call C_get_cell(xcell,ycell,zcell,eps_cell&
         ,beta_sin,beta_cos,beta_term,gamma_sin,gamma_cos,gamma_term&
         ,mono,tric,octa,box34)
#endif
!
!     compute and store reciprocal lattice vectors as columns
!
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
#ifdef _OPENACC
      call copy_recip_cu
#endif
!$acc update device(recip(:,:))

      !get the fractional to Cartesian transformation matrix
      call frac_to_cartgpu

      !find the matrix to convert Cartesian to fractional
      call cart_to_fracgpu
   end if
   if (ranktot.eq.0.and.tinkerdebug.gt.0) then
      print '(A,3F8.2,A,3L3,/,A,2F11.5,A,2F11.5)',&
            ' Cell Info -xyx', xbox,ybox,zbox,' -shape mto',monoclinic,triclinic,octahedron,&
            '           -gamma',gamma,gamma_term,' -beta', beta,beta_term
   end if
!
!     store ctf and tfc matrices
!
   ftcmat(1,1) = ar1
   ftcmat(1,2) = ar2
   ftcmat(1,3) = ar3
   ftcmat(2,1) = br1
   ftcmat(2,2) = br2
   ftcmat(2,3) = br3
   ftcmat(3,1) = cr1
   ftcmat(3,2) = cr2
   ftcmat(3,3) = cr3

   ctfmat(1,1) = 1/xbox
   ctfmat(1,2) = -gamma_cos/(xbox*gamma_sin)
   ctfmat(1,3) = ybox*zbox*(alpha_cos*gamma_cos-beta_cos)/(volbox*gamma_sin)
   ctfmat(2,1) = 0d0
   ctfmat(2,2) = 1/(ybox*gamma_sin)
   ctfmat(2,3) =  xbox*zbox*(beta_cos*gamma_cos-alpha_cos)/(volbox*gamma_sin)
   ctfmat(3,1) = 0d0
   ctfmat(3,2) = 0d0
   ctfmat(3,3) = xbox*ybox*gamma_sin/volbox
   !$acc update device(ctfmat,ftcmat)
!
!     volume of truncated octahedron is half of cubic parent
!
   if (octahedron)  volbox = 0.5_re_p * volbox
!$acc update device(volbox)
end
