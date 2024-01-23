c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ################################################################
c     ##                                                            ##
c     ##  module pme  --  values for particle mesh Ewald summation  ##
c     ##                                                            ##
c     ################################################################
c
c
c     maxorder   maximum order of the B-spline approximation
c
c     bsmod1     B-spline moduli along the a-axis direction
c     bsmod2     B-spline moduli along the b-axis direction
c     bsmod3     B-spline moduli along the c-axis direction
c     thetai1    B-spline coefficients along the a-axis
c     thetai2    B-spline coefficients along the b-axis
c     thetai3    B-spline coefficients along the c-axis
c     qgridin_2d  values on the particle mesh Ewald charge grid, permanent multipoles
c     qgridout_2d  values on the transposed particle mesh Ewald charge grid, permanent multipoles
c     qgrid2in_2d  values on the particle mesh Ewald charge grid, induced dipoles
c     qgrid2out_2d  values on the transposed particle mesh Ewald charge grid, induced dipoles
c     qfac_2d       prefactors for particle mesh Ewald charge grid
c     nfft1      number of grid points along the a-axis direction
c     nfft2      number of grid points along the b-axis direction
c     nfft3      number of grid points along the c-axis direction
c     bsorder    order of the PME B-spline approximation 
c     bseorder   order of the electrostatic PME B-spline values
c     bsporder   order of the polarization PME B-spline values
c     bsdorder   order of the dispersion PME B-spline values
c     cphirec    permanent electric fields, cartesian coordinates
c     fphirec    permanent electric fields, fractional coordinates
c     cphirec0    permanent electric fields, cartesian coordinates, elambda=0
c     fphirec0    permanent electric fields, fractional coordinates, elambda=0
c     cphi0    permanent electric fields, cartesian coordinates, elambda=0
c     cphirec1    permanent electric fields, cartesian coordinates, elambda=1
c     fphirec1    permanent electric fields, fractional coordinates, elambda=1
c     cphi1    permanent electric fields, cartesian coordinates, elambda=1
c     cphidprec    dipolar electric fields, cartesian coordinates
c     fphidprec    dipolar electric fields, fractional coordinates
c     igrid      initial Ewald charge grid values for B-spline
c
c
      module pme
      use sizes
      implicit none
      integer maxorder
      parameter (maxorder=20)
      integer nfft1,nfft2,nfft3
      integer nefft1,nefft2,nefft3
      integer ndfft1,ndfft2,ndfft3
      integer bsorder,bseorder
      integer bsporder,bsdorder
      real*8 bsmod1(maxfft),bsmod2(maxfft),bsmod3(maxfft)
      real*8, allocatable :: thetai1(:,:,:)
      real*8, allocatable :: thetai2(:,:,:)
      real*8, allocatable :: thetai3(:,:,:)
      real*8, allocatable :: qgridin_2d(:,:,:,:,:)
      real*8, allocatable :: qgridout_2d(:,:,:,:)
      real*8, allocatable :: qgrid2in_2d(:,:,:,:,:)
      real*8, allocatable :: qgrid2out_2d(:,:,:,:)
      real*8, allocatable :: qfac_2d(:,:,:)
      real*8, allocatable :: cphirec(:,:),fphirec(:,:)
      real*8, allocatable :: cphi0(:,:),cphi1(:,:)
      real*8, allocatable :: cphirec0(:,:),cphirec1(:,:)
      real*8, allocatable :: fphirec0(:,:),fphirec1(:,:)
      real*8, allocatable :: cphidprec(:,:),fphidprec(:,:)
      integer, allocatable :: igrid(:,:)
      save
      end
