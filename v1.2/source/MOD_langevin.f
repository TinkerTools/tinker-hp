c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ########################################################################
c     ##                                                                    ##
c     ##  module langevin  --  parameters and arrays for langevin dynamics  ##
c     ##                                                                    ##
c     ########################################################################
c
c     gamma : friction parameter in ps-1
c     Rn : white noise 
c
c
      module langevin
      implicit none
      real*8 gamma
      real*8, allocatable :: gamma_friction(:)
      real*8, allocatable :: Rn(:,:)
      save
      end
