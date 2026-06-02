!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     #####################################################################
!     ##                                                                 ##
!     ##  module polar  --  polarizabilities and induced dipole moments  ##
!     ##                                                                 ##
!     #####################################################################
!
!
!     polarity  dipole polarizability for each multipole site (Ang**3)
!     winpolariy    window object corresponding to polariy
!     polarity_orig  original dipole polarizability for each multipole site (Ang**3) (lambdadyn)
!     winpolariy_orig    window object corresponding to polariy_orig
!     thole     Thole polarizability damping value for each site
!     winthole    window object corresponding to thole
!     tholed    Thole direct polarization damping value for each atom
!     wintholed    window object corresponding to tholed
!     pdamp     value of polarizability scale factor for each site
!     winpdamp    window object corresponding to pdamp
!     thlval    Thole damping parameter value for each atom type pair
!     thdval    alternate Thole direct damping value for atom type pair
!     uind      induced dipole components at each multipole site
!     uinp      induced dipoles in field used for energy interactions
!     npolar    total number of polarizable sites in the system
!     tinypol   negligeable polarisability to allow convergence
!     jpolar    index into polarization parameter matrix for each atom
!     winjpolar window object corresponding to jpolar
!     npolar_ne_npole     (npolar.ne.npole)
!     use_mpolar_ker   controls use of mpolar kernel (empole+polar)
!
!
#include "tinker_macro.h"
module polar
   implicit none
   integer   npolar
   logical   npolar_ne_npole,use_mpolar_ker
   real(t_p),parameter :: tinypol = 1d-5
   real(t_p),allocatable,target :: uind(:,:),uinp(:,:)
   real(t_p),pointer :: polarity(:),polarity_orig(:),thole(:),tholed(:)&
            ,pdamp(:),uind_p(:,:),uinp_p(:,:),thlval(:,:),thdval(:,:)
   integer  ,pointer :: jpolar(:)
   integer ::winpolarity,winpolarity_orig,winjpolar&
            ,winthole,winpdamp,wintholed
end

 ! Temporary Data allocated by polarisation
module polar_temp
   implicit none
   integer now
   real(t_p),allocatable,target :: fuind(:,:),fuinp(:,:),fphid(:,:),fphip(:,:)
   real(t_p),allocatable,target :: p_save(:)
   real(t_p),allocatable,target :: ef(:,:,:), mu(:,:,:), murec(:,:,:)
   real(t_p),allocatable :: cphi(:,:)
   real(t_p),allocatable :: res(:,:,:), h(:,:,:), pp(:,:,:), zr(:,:,:), diag(:)
   real(t_p),allocatable,target :: dipfield(:,:,:), dipfieldbis(:,:,:)
   real(t_p),allocatable,target :: cmp(:,:),fmp(:,:),fphidp(:,:), trqrec(:,:)
end module
