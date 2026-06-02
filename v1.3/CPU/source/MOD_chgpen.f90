!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ##################################################################
!     ##                                                              ##
!     ##  module chgpen  --  charge penetration in current structure  ##
!     ##                                                              ##
!     ##################################################################
!
!
module chgpen
   implicit none
   integer :: ncp  !<total number of charge penetration sites in system
   integer :: winpcore !<window corresponding to pcore array
   integer :: winpval !<window corresponding to pval array
   integer :: winpval0 !<window corresponding to pval0
   integer :: winpalpha !<window corresponding to palpha array
   real*8, pointer :: pcore(:) !<number of core electrons at each multipole site
   real*8, pointer :: pval(:) !<number of valence electrons at each multipole site
   real*8, pointer :: pval0(:) !<original number of valence electrons at each multipole site for chglx
   real*8, pointer :: palpha(:) !<charge penetration damping at each multipole site
   save
end
