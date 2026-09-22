!
!     CUDA mirror for authoritative box and replicated-cell geometry
!
#ifdef _CUDA
#define TINKER_CUF
#include "tinker_macro.h"
module cellcu
   use cudafor
   use boxes, only: box34_h=>box34,beta_cos_h=>beta_cos,beta_sin_h=>beta_sin,&
      &beta_term_h=>beta_term,gamma_cos_h=>gamma_cos,gamma_sin_h=>gamma_sin,&
      &gamma_term_h=>gamma_term,monoclinic_h=>monoclinic,&
      &octahedron_h=>octahedron,triclinic_h=>triclinic
   use cell, only: xcell_h=>xcell,ycell_h=>ycell,zcell_h=>zcell,&
      &xcell2_h=>xcell2,ycell2_h=>ycell2,zcell2_h=>zcell2,&
      &i_xcell_h=>i_xcell,i_ycell_h=>i_ycell,i_zcell_h=>i_zcell,&
      &eps_cell_h=>eps_cell
   implicit none

   real(t_p),constant :: xcell,ycell,zcell,xcell2,ycell2,zcell2
   real(t_p),constant :: i_xcell,i_ycell,i_zcell,mincell2,eps_cell
   real(t_p),constant :: beta_sin,beta_cos,beta_term
   real(t_p),constant :: gamma_sin,gamma_cos,gamma_term,box34
   logical,constant :: octahedron,monoclinic,triclinic

contains

   subroutine copy_cell_cu
      implicit none

      xcell=xcell_h
      ycell=ycell_h
      zcell=zcell_h
      xcell2=xcell2_h
      ycell2=ycell2_h
      zcell2=zcell2_h
      mincell2=minval([xcell2_h,ycell2_h,zcell2_h])
      eps_cell=eps_cell_h
      i_xcell=i_xcell_h
      i_ycell=i_ycell_h
      i_zcell=i_zcell_h
      beta_sin=beta_sin_h
      beta_cos=beta_cos_h
      beta_term=beta_term_h
      gamma_sin=gamma_sin_h
      gamma_cos=gamma_cos_h
      gamma_term=gamma_term_h
      monoclinic=monoclinic_h
      triclinic=triclinic_h
      octahedron=octahedron_h
      box34=box34_h
   end subroutine copy_cell_cu

end module cellcu
#endif
