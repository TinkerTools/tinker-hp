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
c     tinypol   negligeable polarisability to allow convergence
c     npolar_ne_npole     (npolar.ne.npole)
c
c
#include "tinker_precision.h"
      module polar
      implicit none
      integer npolar
      logical npolar_ne_npole
      real(t_p), parameter :: tinypol = 1e-5
      integer :: winpolarity,winthole,winpdamp
      real(t_p), allocatable,target :: uind(:,:),uinp(:,:)
      real(t_p), pointer :: polarity(:),thole(:),pdamp(:)
      real(t_p), pointer :: uind_p(:,:),uinp_p(:,:)
      save
      end
 
      ! Temporary Data allocated by polarisation
      module polar_temp
      implicit none
      integer now
      real(t_p), allocatable,target :: fuind(:,:),fuinp(:,:)
      real(t_p), allocatable,target :: fphid(:,:),fphip(:,:)
      real(t_p), allocatable,target :: p_save(:)
      real(t_p), allocatable,target :: ef(:,:,:), mu(:,:,:)
     &         , murec(:,:,:)
      real(t_p), allocatable :: cphi(:,:)
      real(t_p), allocatable :: res(:,:,:), h(:,:,:)
     &         , pp(:,:,:), zr(:,:,:), diag(:)
      real(t_p), allocatable,target :: dipfield(:,:,:)
     &         , dipfieldbis(:,:,:)
      real(t_p), allocatable,target :: cmp(:,:),fmp(:,:),fphidp(:,:)
     &         , trqrec(:,:)
      end module
