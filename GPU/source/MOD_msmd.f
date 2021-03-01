c
c          ###    ##    ##    ###
c         #       # #  # #    #  ##
c        #        #  ##  #    #    #
c         ###     #      #    #    #
c           #     #      #    #    #
c          #      #      #    #  ##
c       ###       #      #    ###
c
c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #################################################################
c     ##                                                             ##
c     ##  msmd.f  --  parameters, energy components                  ## 
c     ##                  and derivatives                            ##
c     ##                                                             ##
c     #################################################################
c
c
c     use_smd_velconst Activation key of the constant velocity SMD in 
c                      the programm
c     use_smd_forconst Activation key of the constant force SMD in
c                      the programm
c     use_smdk2 Use or not of a second k constant
c     use_atfol Print or not specific informations about SMD atoms
c     SMDk    Elastic constant (in Kcal/mol/A^2) 
c     SMDk2   Transverse elastic constant (in Kcal/mol/A^2)
c     SMDVel  Correspond to the velocity of the SMD center of mass in 
c             case of constant velocity SMD study
c     SMDFor  Correspond to the force of the SMD center of mass in 
c             case of constant force SMD study 
c     ncsmd   Number of SMD atoms assigned to the center of mass of SMD
c     tcsmd   Index of each atoms assigned to the SMD COM
c     SMDoutputFreq  Correspond to the frequency in timestseps with 
c                    which the current SMD data values of printed out
c     ismdout unit of the smd output file (default = 3)
c     ismdsave unit of the smd save file (default = 40)
c     xcom    x coordinate of the initial SMD COM
c     ycom    y coordinate of the initial SMD COM
c     zcom    z coordinate of the initial SMD COM
c     cur_xcom x coordinate of the SMD COM at timestep t
c     cur_ycom y coordinate of the SMD COM at timestep t
c     cur_zcom z coordinate of the SMD COM at timestep t
c     xdir    x direction of the center of mass 
c     ydir    y direction of the center of mass
c     zdir    z direction of the center of mass
c     mtotcom initial total mass of the center of mass
c     curmtotcom current total mass of the center of mass during
c                calculations
c     atfol   atom to use with use_atfol during the simulation
c     tabatfol Table containing index of atoms specified with atfol
c     tsmd    time used in the forces calculations (incrmented by 1
c             at each time step (in ps)
c     SMDdt   time step incrementation during the force calculations
c     tpass   number of passage in the routine esmd1.f (usefull for the
c             use of SMDoutputFreq)
c     naltsmd Useful to adapt the timestep for multitimestep procedure
c     dedx    1D table of x force contribution
c     dedy    1D table of y force contribution
c     dedz    1D table of z force contribution
c     com_dedx x force contribution on the SMD COM at timestep t
c     com_dedy y force contribution on the SMD COM at timestep t
c     com_dedz z force contribution on the SMD COM at timestep t
c     stock_dedx storage of the x force contribution on the SMD COM at
c                timestep t     
c     stock_dedy storage of the y force contribution on the SMD COM at
c                timestep t
c     stock_dedz storage of the z force contribution on the SMD COM at
c                timestep t
c     stock_ensmd storage of the SMD energy contribution on the SMD COM
c                 at timestep t
c ### Lines added for the MPI SMD ####
c     nsmdloc number of allocated SMD atoms per processors
c     nsmdglob 1D table corresponded to the number of each SMD atoms
c     smdprocprint Decide which proc has to send the information to the master
c ####################################
c
#include "tinker_precision.h"
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
