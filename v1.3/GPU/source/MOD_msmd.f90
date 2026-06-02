!
!          ###    ##    ##    ###
!         #       # #  # #    #  ##
!        #        #  ##  #    #    #
!         ###     #      #    #    #
!           #     #      #    #    #
!          #      #      #    #  ##
!       ###       #      #    ###
!
!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     #################################################################
!     ##                                                             ##
!     ##  msmd.f  --  parameters, energy components                  ##
!     ##                  and derivatives                            ##
!     ##                                                             ##
!     #################################################################
!
!
!     use_smd_velconst Activation key of the constant velocity SMD in
!                      the programm
!     use_smd_forconst Activation key of the constant force SMD in
!                      the programm
!     use_smdk2 Use or not of a second k constant
!     use_atfol Print or not specific informations about SMD atoms
!     SMDk    Elastic constant (in Kcal/mol/A^2)
!     SMDk2   Transverse elastic constant (in Kcal/mol/A^2)
!     SMDVel  Correspond to the velocity of the SMD center of mass in
!             case of constant velocity SMD study
!     SMDFor  Correspond to the force of the SMD center of mass in
!             case of constant force SMD study
!     ncsmd   Number of SMD atoms assigned to the center of mass of SMD
!     tcsmd   Index of each atoms assigned to the SMD COM
!     SMDoutputFreq  Correspond to the frequency in timestseps with
!                    which the current SMD data values of printed out
!     ismdout unit of the smd output file (default = 3)
!     ismdsave unit of the smd save file (default = 40)
!     xcom    x coordinate of the initial SMD COM
!     ycom    y coordinate of the initial SMD COM
!     zcom    z coordinate of the initial SMD COM
!     cur_xcom x coordinate of the SMD COM at timestep t
!     cur_ycom y coordinate of the SMD COM at timestep t
!     cur_zcom z coordinate of the SMD COM at timestep t
!     xdir    x direction of the center of mass
!     ydir    y direction of the center of mass
!     zdir    z direction of the center of mass
!     mtotcom initial total mass of the center of mass
!     curmtotcom current total mass of the center of mass during
!                calculations
!     atfol   atom to use with use_atfol during the simulation
!     tabatfol Table containing index of atoms specified with atfol
!     tsmd    time used in the forces calculations (incrmented by 1
!             at each time step (in ps)
!     SMDdt   time step incrementation during the force calculations
!     tpass   number of passage in the routine esmd1.f (usefull for the
!             use of SMDoutputFreq)
!     naltsmd Useful to adapt the timestep for multitimestep procedure
!     dedx    1D table of x force contribution
!     dedy    1D table of y force contribution
!     dedz    1D table of z force contribution
!     com_dedx x force contribution on the SMD COM at timestep t
!     com_dedy y force contribution on the SMD COM at timestep t
!     com_dedz z force contribution on the SMD COM at timestep t
!     stock_dedx storage of the x force contribution on the SMD COM at
!                timestep t
!     stock_dedy storage of the y force contribution on the SMD COM at
!                timestep t
!     stock_dedz storage of the z force contribution on the SMD COM at
!                timestep t
!     stock_ensmd storage of the SMD energy contribution on the SMD COM
!                 at timestep t
! ### Lines added for the MPI SMD ####
!     nsmdloc number of allocated SMD atoms per processors
!     nsmdglob 1D table corresponded to the number of each SMD atoms
!     smdprocprint Decide which proc has to send the information to the master
! ####################################
!
#include "tinker_macro.h"
module msmd
   implicit none
   logical use_smdk2, use_atfol
   real(r_p) SMDk, SMDk2
   real(r_p) SMDvel, SMDfor
   integer SMDoutputFreq
   integer ismdout, ismdsave
   real(r_p) xdir, ydir, zdir
   integer ncsmd
   integer, pointer :: tcsmd(:)
   real(r_p) :: xcom, ycom, zcom
   real(r_p) :: cur_xcom, cur_ycom, cur_zcom
   real(r_p) mtotcom, curmtotcom
   integer atfol
   integer, pointer :: tabatfol(:)
   real(r_p) tsmd
   real(r_p) SMDdt
   integer tpass
   real(r_p), pointer :: dedx(:),dedy(:),dedz(:)
   real(r_p) :: com_dedx, com_dedy, com_dedz
   real(r_p) :: stock_dedx, stock_dedy, stock_dedz
   real(r_p) :: stock_ensmd
   integer  nsmdloc
   integer, pointer :: nsmdglob(:)
   integer smdprocprint
   integer naltsmd
   save
end
