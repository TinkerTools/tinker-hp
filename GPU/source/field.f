c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine field  --  get the potential energy functions  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "field" sets the force field potential energy functions from
c     a parameter file and modifications specified in a keyfile
c
c
#include "tinker_macro.h"
      subroutine field
      use ani    ,only: MLpot,use_bondorder
      use keys
      use inform ,only: app_id,dynamic_a,pimd_a
      use potent
      use uprior ,only: use_pred
      implicit none
      integer i
      character*240 record
c
c
c     set the default values for the active potentials
c
      use_bond       = .true.
      use_angle      = .true.
      use_strbnd     = .true.
      use_urey       = .true.
      use_angang     = .true.
      use_opbend     = .true.
      use_opdist     = .true.
      use_improp     = .true.
      use_imptor     = .true.
      use_tors       = .true.
      use_pitors     = .true.
      use_angtor     = .true.
      use_strtor     = .true.
      use_tortor     = .true.
      use_vdw        = .true.
      use_vdwshort   = .false.
      use_vdwlong    = .false.
      use_charge     = .true.
      use_cshortreal = .false.
      use_clong      = .false.
      use_creal      = .true.
      use_crec       = .true.
      use_cself      = .true.
      use_mpole      = .true.
      use_mpoleshortreal = .false.
      use_mpolelong      = .false.
      use_mreal      = .true.
      use_mrec       = .true.
      use_mself      = .true.
      use_polar      = .true.
      use_preal      = .true.
      use_prec       = .true.
      use_pself      = .true.
      use_polarshortreal = .false.
      use_solv       = .false.
      use_geom       = .true.
      use_extra      = .false.
      use_emtp       = .false.
      use_repulsshort    = .false.
      use_repulslong     = .false.
      use_repuls         = .true.
      use_disp           = .true.
      use_dispshortreal  = .false.
      use_displong   = .false.
      use_dispreal   = .true.
      use_disprec    = .true.
      use_dispself   = .true.
      use_chgtrn     = .true.
      use_chgflx     = .true.
      use_lambdadyn  = .false.
      use_OSRW       = .false.
      fuse_chglj     = .false.
      fuse_bonded    = .false.
      use_mlpot      = .false.
c
c     Flags for Machine learning potentials
c
      MLpot           = "ANI2X"
      use_ani_only    = .FALSE.
      use_bondorder   = .FALSE.
      use_ml_embedding= .FALSE.
      use_embd_potoff = .FALSE.
      use_embd_bond   = .TRUE.
      use_embd_angle  = .TRUE.
      use_embd_strbnd = .TRUE.
      use_embd_urey   = .TRUE.
      use_embd_angang = .TRUE.
      use_embd_opbend = .TRUE.
      use_embd_opdist = .TRUE.
      use_embd_improp = .TRUE.
      use_embd_imptor = .TRUE.
      use_embd_tors   = .TRUE.
      use_embd_pitors = .TRUE.
      use_embd_strtor = .TRUE.
      use_embd_tortor = .TRUE.
c
c     Flag for predictor-Corrector
c
      use_pred       = merge(.true.,.false.,app_id.eq.dynamic_a)
c
c     Set default values of Force field potential
c
      bonded_l             = .true.
      shortnonbonded_l     = .false.
      nonbonded_l          = .true.
      PotentialAll         = .true.
      PotentialAmoeba      = .false.
      PotentialAmoeba18    = .false.
      PotentialWaterAmoeba = .false.
      PotentialWaterCharmm = .false.
      PotentialCharmm      = .false.
c
c     read the potential energy force field parameter file
c
      call getprm
c
c     check keywords for potential function control parameters
c
      do i = 1, nkey
         record = keyline(i)
         call prmkey (record)
      end do
!$acc update device(use_mpole,use_polar)
      end
