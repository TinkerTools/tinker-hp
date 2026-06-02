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
!     kelvin      target value for the system temperature (K)
!     atmsph      target value for the system pressure (atm)
!     tautemp     time constant for Berendsen thermostat (psec)
!     taupres     time constant for Berendsen barostat (psec)
!     compress    isothermal compressibility of medium (atm-1)
!     collide     collision frequency for Andersen thermostat
!     volmove     maximum volume move for Monte Carlo barostat (Ang**3)
!     voltrial    mean number of steps between Monte Carlo moves
!     isothermal  logical flag governing use of temperature control
!     isobaric    logical flag governing use of pressure control
!     anisotrop   logical flag governing use of anisotropic pressure
!     thermostat  choice of temperature control method to be used
!     barostat    choice of pressure control method to be used
!     volscale    choice of scaling method for Monte Carlo barostat
!     gammapiston friction for the Langevin Piston barostat
!     masspiston  mass of the piston for the Langevin Piston barostat
!     extvol      Volume as an extended variable for the Langevin Piston barostat
!     extvolold   last value of the Volume as an extended variable for the Langevin Piston barostat
!     vextvol     Speed of the Volume extended variable for the Langevin Piston barostat
!     aextvol     Acceleration of the Volume extended variable for the Langevin Piston barostat
!     freeze_axis info about axis to freeze during volume fluctuations (semi-isotropic barostat)
!
#include "tinker_macro.h"
module bath
   implicit none
   integer voltrial
   integer semi_isotrop_dir
   real(r_p) kelvin,atmsph
   real(r_p) tautemp,taupres
   real(r_p) compress,collide
   real(r_p) vbar,qbar,gbar
   real(r_p) eta,volmove
   real(r_p) gammapiston,masspiston,extvolold
   real(r_p) extvol,vextvol,aextvol
   real(r_p) temppiston
   real(r_p) extbox(3)
   real(r_p) vextbox(3)
   real(r_p) aextbox(3)
   logical :: freeze_axis(3)=.FALSE.
   logical use_piston
   logical isothermal
   logical isobaric
   logical anisotrop
   logical :: isobaric_save=.FALSE.
   logical :: use_piston_save=.FALSE.
   character*9 volscale
   character*11 barostat
   character*11 thermostat
!$acc declare create(eta)
end
