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
c
      module bath
      implicit none
      integer voltrial
      real*8 kelvin,atmsph
      real*8 tautemp,taupres
      real*8 compress,collide
      real*8 vbar,qbar,gbar
      real*8 eta,volmove
      logical isothermal
      logical isobaric
      logical anisotrop
      character*9 volscale
      character*11 barostat
      character*11 thermostat
      save 
      end
