!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ################################################################
!     ##                                                            ##
!     ##  module pme  --  values for particle mesh Ewald summation  ##
!     ##                                                            ##
!     ################################################################
!
!
module pme
   use sizes
   implicit none
   integer :: maxorder !<maximum order of the B-spline approximation
   parameter (maxorder=20)
   integer :: nfft1 !<number of grid points along the a-axis direction
   integer :: nfft2 !<number of grid points along the b-axis direction
   integer :: nfft3 !<number of grid points along the c-axis direction
   integer :: nefft1 !<number of grid points (ele) along the a-axis direction
   integer :: nefft2 !<number of grid points (ele) along the b-axis direction
   integer :: nefft3 !<number of grid points (ele) along the c-axis direction
   integer :: ndfft1 !<number of grid points (disp) along the a-axis direction
   integer :: ndfft2 !<number of grid points (disp) along the b-axis direction
   integer :: ndfft3 !<number of grid points (disp) along the c-axis direction
   integer :: bsorder !<order of the PME B-spline approximation
   integer :: bseorder !<order of the PME B-spline approximation (ele)
   integer :: bsporder !<order of the PME B-spline approximation (polar)
   integer :: bsdorder !<order of the PME B-spline approximation (disp)
   integer, allocatable :: igrid(:,:) !<initial Ewald charge grid values for B-spline
   real*8 :: bsmod1(maxfft) !<B-spline moduli along the a-axis direction
   real*8 :: bsmod2(maxfft) !<B-spline moduli along the b-axis direction
   real*8 :: bsmod3(maxfft) !<B-spline moduli along the c-axis direction
   real*8, allocatable :: thetai1(:,:,:) !<B-spline coefficients along the a-axis
   real*8, allocatable :: thetai2(:,:,:) !<B-spline coefficients along the b-axis
   real*8, allocatable :: thetai3(:,:,:) !<B-spline coefficients along the c-axis
   real*8, allocatable :: qgridin_2d(:,:,:,:,:) !<values on the particle mesh Ewald charge grid, permanent multipoles
   real*8, allocatable :: qgridout_2d(:,:,:,:) !<values on the transposed particle mesh Ewald charge grid, permanent multipoles
   real*8, allocatable :: qgrid2in_2d(:,:,:,:,:) !<values on the particle mesh Ewald charge grid, induced dipoles
   real*8, allocatable :: qgrid2out_2d(:,:,:,:) !<values on the transposed particle mesh Ewald charge grid, induced dipoles
   real*8, allocatable :: qfac_2d(:,:,:) !<prefactors for particle mesh Ewald charge grid
   real*8, allocatable :: cphirec(:,:) !<permanent electric fields, cartesian coordinates
   real*8, allocatable :: fphirec(:,:) !<permanent electric fields, fractional coordinates
   real*8, allocatable :: cphidprec(:,:) !<dipolar electric fields, cartesian coordinates
   real*8, allocatable :: fphidprec(:,:) !<dipolar electric fields, fractional coordinates
   save
end
