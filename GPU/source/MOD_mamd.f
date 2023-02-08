c
c        ##          ##    ##    ###
c       #  #         # #  # #    #  ##
c      #    #        #  ##  #    #    #
c     ########       #      #    #    #
c    #        #      #      #    #    #
c   #          #     #      #    #  ##
c  #            #    #      #    ###
c
c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #################################################################
c     ##                                                             ##
c     ##  mamd.f  --  parameters                                     ## 
c     ##                  for AMD procedure                          ##
c     ##                                                             ##
c     #################################################################
c
c
c     amd_dih_ene       threshold energy E for dihedrals
c     amd_dih_alpha     acceleration factor alpha for dihedrals
c     amd_ep_ene        threshold energy E for potential energy
c     amd_ep_alpha      acceleration factor alpha for potential energy
c     amdoutputfreq     frequency of printing in the AMD output file
c     iamdout           index for aMD output file
c     aMDdt             timestep in the amd
c     tamd              current timestep during the aMD dih
c     tamd2             current timestep during the aMD ep
c     amdtpass          account for the passage in the etors1.f
c     amdtpass2         account for the passage in the
c     amdstate          depict the use or not of the aMD module
c     gamd_ie           energetic threshold for adding boost potential (1 or 2)
c     gamd_cmdprepsteps no of timesteps to prepare cMD (STEP 1)
c     gamd_cmdsteps     no of timesteps to produce cMD (STEP 2)
c     gamd_eqprepsteps  no of timesteps to prepare eqGaMD (STEP 3)
c     gamd_eqsteps      no of timesteps to produce eqGaMD (STEP 4)
c     gamd_sigma0P      max SD allowed for potential energy boost in GaMD
c     gamd_sigma0D      max SD allowed for dihedral boost in GaMD
c     gamd_sigmaE       max SD allowed for dual boost in GaMD
c     use_gamd_restart  specify if is it or not a restart procedure
c     gamd_restartfile  name of the restart GaMD file
c     egamd             counter for storaging of GaMD energy
c     gamdstep          localized which step of gamd we are
c     VmaxD             max GaMD potential
c     VminD             min GaMD potential
c     VmaxP
c     VminP
c     Vavg              averrage GaMD potential
c     Vdiff             storaging stuff for sigmaV
c     M2                storaging stuff for sigmaV
c     sigmaV            SD of the GaMD potential
c     sigma0            SD used in the GaMD calculation
c     cptgamd           counter for GaMD
c     gamdE
c     gamdk
c     gamdk0
c     gamdk0bis
c     gamdk0ter
c     gamd_deltaV
c     aMDwattype
c
c
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


