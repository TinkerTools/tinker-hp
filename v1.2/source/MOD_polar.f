c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #####################################################################
c     ##                                                                 ##
c     ##  module polar  --  polarizabilities and induced dipole moments  ##
c     ##                                                                 ##
c     #####################################################################
c
c
c     polarity  dipole polarizability for each multipole site (Ang**3)
c     winpolariy    window object corresponding to polariy
c     thole     Thole polarizability damping value for each site
c     winthole    window object corresponding to thole
c     pdamp     value of polarizability scale factor for each site
c     winpdamp    window object corresponding to pdamp
c     uind      induced dipole components at each multipole site
c     uinp      induced dipoles in field used for energy interactions
c     npolar    total number of polarizable sites in the system
c     dirdamp   direct polarization damping value for each site
c     windirdamp window object corresponding to dirdamp
c     tinypol   negligeable polarisability to allow convergence
c
c
      module polar
      implicit none
      integer npolar
      real*8, parameter :: tinypol = 1e-5
      real*8, pointer :: polarity(:),thole(:),pdamp(:)
      integer :: winpolarity,winthole,winpdamp
      real*8, allocatable :: uind(:,:),uinp(:,:)
      real*8, pointer :: dirdamp(:)
      integer :: windirdamp
      save
      end
