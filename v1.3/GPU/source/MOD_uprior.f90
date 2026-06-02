!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ####################################################################
!     ##                                                                ##
!     ##  module uprior  --  previous values of induced dipole moments  ##
!     ##                                                                ##
!     ####################################################################
!
!
!     maxualt   maximum number of sets of induced dipoles to save
!     dGuSp     dimension of the Guess subspace
!
!     gear      coefficients for Gear predictor binomial method
!     aspc      coefficients for always stable predictor-corrector
!     bpred     coefficients for induced dipole predictor polynomial
!     bpredp    coefficients for predictor polynomial in energy field
!     bpreds    coefficients for predictor for PB/GK solvation
!     bpredps   coefficients for predictor in PB/GK energy field
!     udalt     prior values for induced dipoles at each site
!     upalt     prior values for induced dipoles in energy field
!     udshortalt  prior values for short range induced dipoles at each site
!     upshortalt  prior values for short range induced dipoles in energy field
!     nualt     number of prior sets of induced dipoles in storage
!     nualt1    number of prior sets of short range induced dipoles in storage
!     idxGuSp   last update index within the induce dipole guess/Solution subspace  
!     shortnualt  number of prior sets of short range induced dipoles in storage
!     upalt_p0   Pointer to upalt or upshortalt
!     udalt_p0   Pointer to udalt or udshortalt
!     use_pred  flag to control use of induced dipole prediction
!     polpred   type of predictor polynomial (Gear, ASPC or LSQR)
!     lalt      save location index of electrical field
!     lshalt    save location index of short range electrical field
!
!
#include "tinker_macro.h"
module uprior
#ifdef USE_NVSHMEM_CUDA
   use tinTypes ,only: r3dDPC=>Real3dDevPointerContainer
#endif
   implicit none
   integer :: maxualt=-1,maxualt_prm,dGuSp
   parameter (maxualt_prm=32)
   integer nualt, nualt1, idxGuSp
   integer lalt,lshalt
   !TODO Should be used with short range solver
   ! integer shortnualt
   real(t_p) gear(maxualt_prm),aspc(maxualt_prm)
   real(t_p),target:: bpred (maxualt_prm),bpredp (maxualt_prm)
   real(t_p),target:: bpreds(maxualt_prm),bpredps(maxualt_prm)
   real(t_p),allocatable,target :: udalt(:,:,:),udshortalt(:,:,:)
   real(t_p),allocatable,target :: upalt(:,:,:),upshortalt(:,:,:)
   real(t_p),allocatable,target :: udpalt(:), udpalt1(:), Audp(:)
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
