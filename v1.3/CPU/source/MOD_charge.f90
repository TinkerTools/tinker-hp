!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ####################################################################
!     ##                                                                ##
!     ##  module charge  --  partial charges for the current structure  ##
!     ##                                                                ##
!     ####################################################################
!
!
!     iion      !<number of the atom site for each partial charge
!     jion      !<neighbor generation site for each partial charge
!     kion      !<cutoff switching site for each partial charge
!     chglist   !<partial charge site for each atom (0=no charge)
!     nbchg     !<number of charges before each index
!     chgloc    !<global-local charge correspondance
!     chglocnl  !<global-localnl charge correspondance
!     chgrecloc  !<global-local reciprocal charge correspondance
!     pchg      !<magnitude of the partial charges (e-)
!     pchg_orig   !<original magnitude of the partial charges (e-) (lambda dyn)
!     pchg0     !<original partial charge values for charge flux
!
!
module charge
   implicit none
   integer :: nion !<total number of partial charges in system
   integer :: nionloc !<local number of partial charges in system
   integer :: nionbloc !<local+neighbors number of partial charges in system
   integer :: nionlocnl !<localnl number of partial charges in system
   integer :: nionrecloc !<local reciprocal number of partial charges in system
   integer :: winiion !<window object corresponding to iion
   integer :: winjion !<window object corresponding to jion
   integer :: winkion !<cutoff switching site for each partial charge
   integer :: winchglist !<window object corresponding to chglist
   integer :: winnbchg !<window object corresponding to nbchg
   integer :: winchgloc !<window object corresponding to chgloc
   integer :: winchglocnl !<window object corresponding to chglocnl
   integer :: winchgrecloc !<window object corresponding to chgrecloc
   integer :: winpchg !<window object corresponding to pchg
   integer :: winpchg0 !<window object corresponding to pchg0
   integer :: winpchg_orig !<window object corresponding to pchg_orig
   integer, allocatable :: chgloc(:) !<global-local charge correspondance
   integer, allocatable :: chglocnl(:) !<global-localnl charge correspondance
   integer, allocatable :: chgrecloc(:) !<global-local reciprocal charge correspondance
   integer, pointer :: chglist(:) !<partial charge site for each atom (0=no charge)
   integer, pointer :: iion(:) !<number of the atom site for each partial charge
   integer, pointer :: jion(:) !<neighbor generation site for each partial charge
   integer, pointer :: kion(:) !<cutoff switching site for each partial charge
   integer, pointer :: nbchg(:) !<number of charges before each index
   real*8, pointer ::  pchg(:) !<magnitude of the partial charges (e-)
   real*8, pointer :: pchg0(:) !<original partial charge values for charge flux
   real*8, pointer :: pchg_orig(:) !<original magnitude of the partial charges (e-) (lambda dyn)
   save
end
