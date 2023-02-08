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
c     cphidprec    dipolar electric fields, cartesian coordinates
c     fphidprec    dipolar electric fields, fractional coordinates
c     igrid      initial Ewald charge grid values for B-spline
c     pme_eps offset used to shift sites off exact lattice bounds
c
c
#include "tinker_macro.h"
      module pme
      use sizes
      implicit none
      integer maxorder
#ifdef _OPENACC
      logical,parameter::GridDecomp1d=.true.
#else
      logical,parameter::GridDecomp1d=.false.
#endif
      parameter (maxorder=20)
      integer nfft1,nfft2,nfft3
      integer nefft1,nefft2,nefft3
      integer ndfft1,ndfft2,ndfft3
      integer bsorder,bseorder,bsporder,bsdorder
#if (defined(SINGLE)||defined(MIXED))
      real(t_p) pme_eps
#else
      real(t_p),parameter:: pme_eps=1d-8
#endif
!$acc declare create(nfft1,nfft2,nfft3,bsorder)
      real(t_p) bsmod1(maxfft),bsmod2(maxfft),bsmod3(maxfft)
!DIR$ ATTRIBUTES ALIGN:64:: thetai1
      real(t_p), allocatable,target :: thetai1(:,:,:)
!DIR$ ATTRIBUTES ALIGN:64:: thetai3
      real(t_p), allocatable,target :: thetai3(:,:,:)
!DIR$ ATTRIBUTES ALIGN:64:: thetai2
      real(t_p), allocatable,target :: thetai2(:,:,:)
      real(t_p), pointer:: thetai1_p(:,:,:),thetai2_p(:,:,:)
     &         , thetai3_p(:,:,:)
!DIR$ ATTRIBUTES ALIGN:64:: qgridin_2d  
      real(t_p), allocatable,target :: qgridin_2d(:,:,:,:,:)
      real(t_p), pointer :: qgridin_p(:)
!DIR$ ATTRIBUTES ALIGN:64:: qgridout_2d 
      real(t_p), allocatable,target :: qgridout_2d(:,:,:,:)
      real(t_p), pointer :: qgridout_p(:)
!DIR$ ATTRIBUTES ALIGN:64:: qgrid2in_2d 
      real(t_p), allocatable,target :: qgrid2in_2d(:,:,:,:,:)
      real(t_p), pointer :: qgrid2in_p(:)
!DIR$ ATTRIBUTES ALIGN:64:: qgrid2out_2d
      real(t_p), allocatable,target :: qgrid2out_2d(:,:,:,:)
      real(t_p), pointer :: qgrid2out_p(:)
!DIR$ ATTRIBUTES ALIGN:64:: qfac_2d
      real(t_p), allocatable :: qfac_2d(:,:,:)
!DIR$ ATTRIBUTES ALIGN:64:: cphirec,fphirec
      real(t_p), allocatable :: cphirec(:,:),fphirec(:,:)
     &         , cphidprec(:,:),fphidprec(:,:)
!DIR$ ATTRIBUTES ALIGN:64:: igrid  
      integer, allocatable,target :: igrid(:,:)

!$acc declare create(bsmod1,bsmod2,bsmod3)

      end

c     Is_MPI_GRID1D  keeps track of MPI_GRID1D allocation
c     Is_MPI_GRID2D  keeps track of MPI_GRID2D allocation
c     MPI_GRID1D     MPI vector Type for 1d fftGrid
c     MPI_GRID2D     MPI vector Type for 2d fftGrid
      module pme1
      implicit none
      real(t_p),allocatable,target:: qgridmpi(:,:,:,:,:)
      logical :: Is_MPI_GRID1D=.false.,Is_MPI_GRID2D=.false.
      integer MPI_GRID1D,MPI_GRID2D
      integer nfloor,stride
      enum,bind(C)
         enumerator r_comm,r_wait
      end enum
      end module
