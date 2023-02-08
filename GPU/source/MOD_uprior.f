c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ####################################################################
c     ##                                                                ##
c     ##  module uprior  --  previous values of induced dipole moments  ##
c     ##                                                                ##
c     ####################################################################
c
c
c     maxualt   maximum number of sets of induced dipoles to save
c
c     gear      coefficients for Gear predictor binomial method
c     aspc      coefficients for always stable predictor-corrector
c     bpred     coefficients for induced dipole predictor polynomial
c     bpredp    coefficients for predictor polynomial in energy field
c     bpreds    coefficients for predictor for PB/GK solvation
c     bpredps   coefficients for predictor in PB/GK energy field
c     udalt     prior values for induced dipoles at each site
c     upalt     prior values for induced dipoles in energy field
c     udshortalt  prior values for short range induced dipoles at each site
c     upshortalt  prior values for short range induced dipoles in energy field
c     nualt     number of prior sets of induced dipoles in storage
c     shortnualt  number of prior sets of short range induced dipoles in storage
c     upalt_p0   Pointer to upalt or upshortalt
c     udalt_p0   Pointer to udalt or udshortalt
c     use_pred  flag to control use of induced dipole prediction
c     polpred   type of predictor polynomial (Gear, ASPC or LSQR)
c     lalt      save location index of electrical field
c     lshalt    save location index of short range electrical field
c
c
#include "tinker_macro.h"
      module uprior
#ifdef USE_NVSHMEM_CUDA
      use tinTypes ,only: r3dDPC=>Real3dDevPointerContainer
#endif
      implicit none
      integer maxualt
      parameter (maxualt=7)
      integer nualt
      integer lalt,lshalt
      !TODO Should be used with short range solver
      ! integer shortnualt
      real(t_p) gear(maxualt),aspc(maxualt)
      real(t_p) bpred(maxualt),bpredp(maxualt)
      real(t_p) bpreds(maxualt),bpredps(maxualt)
      !DIR$ ATTRIBUTES ALIGN:64:: udalt,udshortalt
      real(t_p),allocatable,target :: udalt(:,:,:),udshortalt(:,:,:)
      !DIR$ ATTRIBUTES ALIGN:64:: upalt,upshortalt
      real(t_p),allocatable,target :: upalt(:,:,:),upshortalt(:,:,:)
      real(t_p),pointer :: udalt_p0(:,:,:),upalt_p0(:,:,:)
      logical use_pred
      character*4 polpred

#ifdef USE_NVSHMEM_CUDA
      ![c][d]_*  nvshmem data structure for *
      type(r3dDPC),device,pointer::d_udalt(:),d_udshortalt(:)
      type(r3dDPC),   allocatable::c_udalt(:),c_udshortalt(:)
      type(r3dDPC),device,pointer::d_upalt(:),d_upshortalt(:)
      type(r3dDPC),   allocatable::c_upalt(:),c_upshortalt(:)
      ! Pointers to d_udalt,d_udshortalt,d_upalt,d_upshortalt
      type(r3dDPC),device,pointer::upalt_p1(:),udalt_p1(:)
      ! Temp buffer to reduce (Remote Memory Adressing)
      type(r3dDPC),device,pointer::d_altbuf(:)
      type(r3dDPC),   allocatable::c_altbuf(:)
#endif

!$acc declare create(bpred,gear,aspc)
!$acc declare create(bpreds,bpredp,bpredps)
      end
