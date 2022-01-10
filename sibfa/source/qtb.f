
c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ########################################################################
c     ##                                                                    ##
c     ##  module QTB  --  parameters and arrays for Quantum Thermal Bath    ##
c     ##                                                                    ##
c     ########################################################################
c     omegacut : cutting frequency for the QTB random force
c     domega : frequency discretization in THz (ps-1)
c     omegasmear : length of the smearing window near omegacut
c     skipseg : number to seg to be skipped (because the system has not
c     vad : past velocities used for correlation computations
c
c     Literature reference: 
c     J-L Barrat & D. Rodney JStatPhys (2011)
c     and H Dammak PRL 103, 190601 (2009)
c
      module qtb 
      implicit none
      real*8 omegacut
      real*8 domega
      real*8 omegasmear
      real*8 omegamax
      real*8, allocatable :: Htilde(:,:)
      real*8, allocatable :: rt(:,:,:)
      real*8, pointer :: noise(:,:,:)
      integer :: nseg
      integer, allocatable :: repartnoise(:)
      logical :: noQTB
      logical adaptive
      logical :: corr_pot
      real*8, allocatable :: corr_pot_ratio(:)
      logical :: register_spectra
      logical :: ir
      integer compteur, skipseg, startsavespec
      real*8, allocatable :: Cmumu_average(:,:)
      real*8, allocatable :: vad(:,:,:)
      real*8 Tseg
      integer :: nad
      save 
      end
     
