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
c     use_pred  flag to control use of induced dipole prediction
c     polpred   type of predictor polynomial (Gear, ASPC or LSQR)
c
c
      module uprior
      implicit none
      integer maxualt
      parameter (maxualt=7)
      integer nualt
      real*8 gear(maxualt),aspc(maxualt)
      real*8 bpred(maxualt),bpredp(maxualt)
      real*8 bpreds(maxualt),bpredps(maxualt)
      real*8, allocatable :: udalt(:,:,:),udshortalt(:,:,:)
      real*8, allocatable :: upalt(:,:,:),upshortalt(:,:,:)
      logical use_pred
      character*4 polpred
      save
      end
