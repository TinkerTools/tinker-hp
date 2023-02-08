c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ####################################################################
c     ##                                                                ##
c     ##  module bath  --  temperature and pressure control parameters  ##
c     ##                                                                ##
c     ####################################################################
c
c
c
c     kelvin      target value for the system temperature (K)
c     atmsph      target value for the system pressure (atm)
c     tautemp     time constant for Berendsen thermostat (psec)
c     taupres     time constant for Berendsen barostat (psec)
c     compress    isothermal compressibility of medium (atm-1)
c     collide     collision frequency for Andersen thermostat
c     volmove     maximum volume move for Monte Carlo barostat (Ang**3)
c     voltrial    mean number of steps between Monte Carlo moves
c     isothermal  logical flag governing use of temperature control
c     isobaric    logical flag governing use of pressure control
c     anisotrop   logical flag governing use of anisotropic pressure
c     thermostat  choice of temperature control method to be used
c     barostat    choice of pressure control method to be used
c     volscale    choice of scaling method for Monte Carlo barostat
c     gammapiston friction for the Langevin Piston barostat
c     masspiston  mass of the piston for the Langevin Piston barostat
c     extvol      Volume as an extended variable for the Langevin Piston barostat
c     extvolold   last value of the Volume as an extended variable for the Langevin Piston barostat
c     vextvol     Speed of the Volume extended variable for the Langevin Piston barostat
c     aextvol     Acceleration of the Volume extended variable for the Langevin Piston barostat
c
#include "tinker_macro.h"
      module bath
      implicit none
      integer voltrial
      real(r_p) kelvin,atmsph
      real(r_p) tautemp,taupres
      real(r_p) compress,collide
      real(r_p) vbar,qbar,gbar
      real(r_p) eta,volmove
      real(r_p) gammapiston,masspiston,extvolold
      real(r_p) extvol,vextvol,aextvol
      real(r_p) temppiston
      logical use_piston
      logical isothermal
      logical isobaric
      logical anisotrop
      character*9 volscale
      character*11 barostat
      character*11 thermostat
!$acc declare create(eta)
      end
