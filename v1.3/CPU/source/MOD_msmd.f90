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
module msmd
   implicit none
   logical :: use_smdk2 !<Use or not of a second k constant
   logical :: use_atfol !<Print or not specific informations about SMD atoms
   integer :: SMDoutputFreq !<Correspond to the frequency in timestseps with which the current SMD data values of printed out
   integer :: ismdout !<unit of the smd output file (default = 3)
   integer :: ismdsave !<unit of the smd save file (default = 40)
   integer :: ncsmd !<Number of SMD atoms assigned to the center of mass of SMD
   integer :: atfol !<atom to use with use_atfol during the simulation
   integer :: tpass !<number of passage in the routine esmd1.f (usefull for the use of SMDoutputFreq)
   integer :: ismd !<temporary SMD atom number during the processors distrubution
   integer :: nsmdloc !<number of allocated SMD atoms per processors
   integer :: smdprocprint !<Decide which proc has to send the information to the master
   integer :: naltsmd !<Useful to adapt the timestep for multitimestep procedure
   integer, pointer :: nsmdglob(:) !<1D table corresponded to the number of each SMD atoms
   integer, pointer :: tabatfol(:) !<Table containing index of atoms specified with atfol
   integer, pointer :: tcsmd(:) !<Index of each atoms assigned to the SMD COM
   real*8 :: SMDk  !<Elastic constant (in Kcal/mol/A^2)
   real*8 :: SMDk2 !<Transverse elastic constant (in Kcal/mol/A^2)
   real*8 :: SMDvel !<Correspond to the velocity of the SMD center of mass in case of constant velocity SMD study
   real*8 :: SMDfor !<Correspond to the force of the SMD center of mass in case of constant force SMD study
   real*8 :: xdir !<x direction of the center of mass
   real*8 :: ydir !<y direction of the center of mass
   real*8 :: zdir !<z direction of the center of mass
   real*8 :: xcom !<x coordinate of the initial SMD COM
   real*8 :: ycom !<y coordinate of the initial SMD COM
   real*8 :: zcom !<z coordinate of the initial SMD COM
   real*8 :: cur_xcom !<x coordinate of the SMD COM at timestep t
   real*8 :: cur_ycom !<y coordinate of the SMD COM at timestep t
   real*8 :: cur_zcom !<z coordinate of the SMD COM at timestep t
   real*8 :: mtotcom !<initial total mass of the center of mass
   real*8 :: curmtotcom !<current total mass of the center of mass during calculations
   real*8 :: tsmd !<time used in the forces calculations (incrmented by 1 at each time step (in ps)
   real*8 :: SMDdt !<time step incrementation during the force calculations
   real*8 :: com_dedx !<x force contribution on the SMD COM at timestep t
   real*8 :: com_dedy !<y force contribution on the SMD COM at timestep t
   real*8 :: com_dedz !<z force contribution on the SMD COM at timestep t
   real*8 :: stock_dedx !<storage of the x force contribution on the SMD COM at timestep t
   real*8 :: stock_dedy !<storage of the y force contribution on the SMD COM at timestep t
   real*8 :: stock_dedz !<storage of the z force contribution on the SMD COM at timestep t
   real*8 :: stock_ensmd !<storage of the SMD energy contribution on the SMD COM at timestep t
   real*8, pointer :: dedx(:) !<1D table of x force contribution
   real*8, pointer :: dedy(:) !<1D table of y force contribution
   real*8, pointer :: dedz(:) !<1D table of z force contribution
   save
end
