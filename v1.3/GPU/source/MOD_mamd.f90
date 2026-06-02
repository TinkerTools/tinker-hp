!
!        ##          ##    ##    ###
!       #  #         # #  # #    #  ##
!      #    #        #  ##  #    #    #
!     ########       #      #    #    #
!    #        #      #      #    #    #
!   #          #     #      #    #  ##
!  #            #    #      #    ###
!
!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     #################################################################
!     ##                                                             ##
!     ##  mamd.f  --  parameters                                     ##
!     ##                  for AMD procedure                          ##
!     ##                                                             ##
!     #################################################################
!
!
!     amd_dih_ene       threshold energy E for dihedrals
!     amd_dih_alpha     acceleration factor alpha for dihedrals
!     amd_ep_ene        threshold energy E for potential energy
!     amd_ep_alpha      acceleration factor alpha for potential energy
!     amdoutputfreq     frequency of printing in the AMD output file
!     iamdout           index for aMD output file
!     aMDdt             timestep in the amd
!     tamd              current timestep during the aMD dih
!     tamd2             current timestep during the aMD ep
!     amdtpass          account for the passage in the etors1.f
!     amdtpass2         account for the passage in the
!     amdstate          depict the use or not of the aMD module
!     gamd_ie           energetic threshold for adding boost potential (1 or 2)
!     gamd_cmdprepsteps no of timesteps to prepare cMD (STEP 1)
!     gamd_cmdsteps     no of timesteps to produce cMD (STEP 2)
!     gamd_eqprepsteps  no of timesteps to prepare eqGaMD (STEP 3)
!     gamd_eqsteps      no of timesteps to produce eqGaMD (STEP 4)
!     gamd_sigma0P      max SD allowed for potential energy boost in GaMD
!     gamd_sigma0D      max SD allowed for dihedral boost in GaMD
!     gamd_sigmaE       max SD allowed for dual boost in GaMD
!     use_gamd_restart  specify if is it or not a restart procedure
!     gamd_restartfile  name of the restart GaMD file
!     egamd             counter for storaging of GaMD energy
!     gamdstep          localized which step of gamd we are
!     VmaxD             max GaMD potential
!     VminD             min GaMD potential
!     VmaxP
!     VminP
!     Vavg              averrage GaMD potential
!     Vdiff             storaging stuff for sigmaV
!     M2                storaging stuff for sigmaV
!     sigmaV            SD of the GaMD potential
!     sigma0            SD used in the GaMD calculation
!     cptgamd           counter for GaMD
!     gamdE
!     gamdk
!     gamdk0
!     gamdk0bis
!     gamdk0ter
!     gamd_deltaV
!     aMDwattype
!
!
#include "tinker_macro.h"
module mamd
   implicit none
   logical amddebug
   integer ncamd
   integer, pointer :: tcamd(:)
   real(r_p) amd_dih_ene
   real(t_p) amd_dih_alpha
   real(r_p) amd_ep_ene
   real(r_p) amd_ep_alpha
   real(r_p) amd_factor_dih
   real(r_p) amd_factor_tot
   logical save_amd_dih, save_amd_ene
   integer amdoutputfreq, amdoutputfreq_dih
   integer iamdout, igamdrestart
   integer amdtpass, amdtpass_dih
   real(r_p) aMDdt_dih
   real(r_p) aMDdt
   real(8) etamd, tamd, tamd_dih
   character (len=3) amdstate
   real(t_p) amdboost, amdboostavg, amdboostavg_dih
   real(t_p) amdboostavg_W1
   integer gamd_ie
   integer gamd_cmdprepsteps, gamd_cmdsteps
   integer gamd_eqprepsteps, gamd_eqsteps
   real(t_p) gamd_sigma0P, gamd_sigma0D, gamd_sigma0E
   real(t_p) gamd_sigma0W1
   logical use_gamd_restart
   character (len=100) gamd_restartfile
   integer gamdstep
   real(8) egamd
   real(r_p) VminD, VmaxD, VavgD, VdiffD, M2D
   real(r_p) VminP, VmaxP, VavgP, VdiffP, M2P
   real(r_p) VminW1, VmaxW1, VavgW1, VdiffW1, M2W1
   real(r_p) sigmaVD, sigma0D
   real(r_p) sigmaVP, sigma0P
   real(r_p) sigmaVW1, sigma0W1
   real(8)   cptgamdW1,cptgamdD,cptgamdP
   real(r_p) gamdED, gamdkD, gamdk0D, gamdk0bisD, gamdk0terD
   real(r_p) gamdEP, gamdkP, gamdk0P, gamdk0bisP, gamdk0terP
   real(r_p) gamdEW1, gamdkW1, gamdk0W1, gamdk0bisW1, gamdk0terW1
   real(r_p) gamd_deltaV
   real(t_p) amdboostD, amdboostP, amdboostavgD, amdboostavgP
   real(t_p) amdboostW1, amdboostavgW1
   real(r_p) gamd_factor_dih, gamd_factor_tot, gamd_factor_wat1
   integer aMDwattype(2)
end


