!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ####################################################################
!     ##                                                                ##
!     ##  module bath  --  temperature and pressure control parameters  ##
!     ##                                                                ##
!     ####################################################################
!
!
!
module bath
   implicit none
   integer ::  voltrial !<mean number of steps between Monte Carlo moves
   real*8 :: kelvin !<target value for the system temperature (K)
   real*8 :: atmsph !<target value for the system pressure (atm)
   real*8 :: tautemp !<time constant for Berendsen thermostat (psec)
   real*8 :: taupres !<time constant for Berendsen barostat (psec)
   real*8 :: compress !<isothermal compressibility of medium (atm-1)
   real*8 :: collide !<collision frequency for Andersen thermostat
   real*8 :: eta !<eta parameter for BUSSI barostat
   real*8 :: volmove !<maximum volume move for Monte Carlo barostat (Ang**3)
   real*8 :: gammapiston !<friction for the Langevin Piston barostat
   real*8 :: masspiston !<mass of the piston for the Langevin Piston barostat
   real*8 :: extvolold !<last value of the Volume as an extended variable for the Langevin Piston barostat
   real*8 :: temppiston !<temperature of the piston extended variable (Langevin barostat)
   real*8 :: extvol !<Volume as an extended variable for the Langevin Piston barostat
   real*8 :: vextvol !<Velocity of the Volume extended variable for the Langevin Piston barostat
   real*8 :: aextvol !<acceleration of the Volume extended variable for the Langevin Piston barostat
   real*8 :: extbox(3) !<sizes of the box as extended variables for the anisotropic Langevin barostat
   real*8 :: vextbox(3) !<velocities of the sizes of the boxes for the anisotropic Langevin barostat
   real*8 :: aextbox(3) !<accelerations of the sizes of the box for the anisotropic Langevin barostat
   logical :: isothermal !<logical flag governing use of temperature control
   logical :: isobaric !<logical flag governing use of pressure control
   logical :: anisotrop !<logical flag governing use of anisotropic pressure
   logical :: use_piston !<logical flag governing use of Langevin Piston barostat
   logical :: freeze_axis(3)=.FALSE. !<info about axis to freeze during volume fluctuations (semi-isotropic barostat)
   logical :: isobaric_save=.FALSE. !<logical flag governing use of Langevin Piston barostat with QTB 
   logical :: use_piston_save=.FALSE. !<logical flag governing monitoring of stat about Langevin piston
   character*9 :: volscale !<choice of scaling method for Monte Carlo barostat
   character*11 :: barostat !<choice of pressure control method to be used
   character*11 :: thermostat !<choice of temperature control method to be used
   save
end
