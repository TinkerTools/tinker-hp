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
module uprior
   implicit none
   integer :: maxualt !<maximum number of sets of induced dipoles to save
   parameter (maxualt=7)
   integer :: nualt !<number of prior sets of induced dipoles in storage
   real*8 :: gear(maxualt) !<coefficients for Gear predictor binomial method
   real*8 :: aspc(maxualt) !<coefficients for always stable predictor-corrector
   real*8 :: bpred(maxualt) !<coefficients for induced dipole predictor polynomial
   real*8 :: bpredp(maxualt) !<coefficients for predictor polynomial in energy field
   real*8 :: bpreds(maxualt) !<coefficients for predictor for PB/GK solvation
   real*8 :: bpredps(maxualt) !<coefficients for predictor in PB/GK energy field
   real*8, allocatable :: udalt(:,:,:) !<prior values for induced dipoles at each site
   real*8, allocatable :: udshortalt(:,:,:) !<prior values for short range induced dipoles at each site
   real*8, allocatable :: upalt(:,:,:) !<prior values for induced dipoles in energy field
   real*8, allocatable :: upshortalt(:,:,:) !<prior values for short range induced dipoles in energy field
   logical :: use_pred !<flag to control use of induced dipole prediction
   character*4 :: polpred !<type of predictor polynomial (Gear, ASPC or LSQR)
   save
end
