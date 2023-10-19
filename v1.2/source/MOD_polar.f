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
c     polarity_orig  original dipole polarizability for each multipole site (Ang**3) (lambdadyn)
c     winpolariy_orig    window object corresponding to polariy_orig
c     thole     Thole polarizability damping value for each site
c     winthole    window object corresponding to thole
c     tholed    Thole direct polarization damping value for each atom
c     wintholed    window object corresponding to tholed
c     pdamp     value of polarizability scale factor for each site
c     winpdamp    window object corresponding to pdamp
c     thlval    Thole damping parameter value for each atom type pair
c     thdval    alternate Thole direct damping value for atom type pair
c     uind      induced dipole components at each multipole site
c     uinp      induced dipoles in field used for energy interactions
c     npolar    total number of polarizable sites in the system
c     tinypol   negligeable polarisability to allow convergence
c     jpolar    index into polarization parameter matrix for each atom
c     winjpolar window object corresponding to jpolar
c
c
      module polar
      implicit none
      integer npolar
      real*8, parameter :: tinypol = 1e-5
      real*8, pointer :: polarity(:),polarity_orig(:),thole(:),pdamp(:)
      integer :: winpolarity,winpolarity_orig,winthole,winpdamp
      real*8, allocatable :: uind(:,:),uinp(:,:)
      real*8, allocatable :: thlval(:,:),thdval(:,:)
      real*8, pointer :: tholed(:)
      integer :: wintholed
      integer, pointer :: jpolar(:)
      integer :: winjpolar
      save
      end
