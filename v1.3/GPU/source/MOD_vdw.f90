!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ######################################################################
!     ##                                                                  ##
!     ##  module vdw  --  van der Waals parameters for current structure  ##
!     ##                                                                  ##
!     ######################################################################
!
!
!     radmin     minimum energy distance for each atom class pair
!     winradmin    window object corresponding to radmin
!     epsilon    well depth parameter for each atom class pair
!     winepsilon    window object corresponding to epsilon
!     radmin4    minimum energy distance for 1-4 interaction pairs
!     winradmin4    window object corresponding to radmin4
!     epsilon4   well depth parameter for 1-4 interaction pairs
!     winepsilon4    window object corresponding to epsilon4
!     radmin_c   same as radmin but with compress class
!     winradmin_c    window object corresponding to radmin_c
!     epsilon_c  same as epsilon but with compress class
!     radmin4_c  same as radmin4 but with compress class
!     epsilon4_c same as epsilon4 but with compress class
!     radhbnd    minimum energy distance for hydrogen bonding pairs
!     winradhbnd    window object corresponding to radhbnd
!     epshbnd    well depth parameter for hydrogen bonding pairs
!     winepshbnd    window object corresponding to epshbnd
!     kred       value of reduction factor parameter for each atom
!     winkred    window object corresponding to kred
!     ired       attached atom from which reduction factor is applied
!     winired    window object corresponding to ired
!     nvdw       total number van der Waals active sites in the system
!     nvdwloc    local number van der Waals active sites in the system
!     nvdwbloc   local+neighbors number van der Waals active sites in the system
!     nvdwblocloop   First multiple of 16 after nvdwbloc (help vectorization
!     nvdwclass  total number of van der Waals different atoms class inside the system
!     ivdw       number of the atom for each van der Waals active site
!     winivdw    window object corresponding to ivdw
!     jvdw       type or class index into vdw parameters for each atom
!     winjvdw    window object corresponding to jvdw
!     jvdw_c     renumbered type or class index into vdw parameters for each atom
!     nvt        number of distinct vdw types/classes in the system
!     ivt        type/class index for each distinct vdw type or class
!     jvt        frequency of each vdw type or class in the system
!     nbvdw      number of 'vdw' atoms before each atom
!     winnbvdw    window object corresponding to nbvdw
!     vdwlocnl   glob-locnl vdw correspondance
!     nvdwlocnlb First multiple of BLOCK_SIZE after nvdwlocnl
!     nvdwlocnlb_pair  total number of vdw pair blocks interaction
!     nvdwlocnlb2_pair  total number of vdw pair blocks interaction in C2 nblist
!     nshortvdwlocnlb2_pair  total number of vdw pair blocks interaction in short range interaction list
!     skipvdw12  switch to skip vdw 1-2 Interactions computation
!     vdw_lcut2  defines a lower bound under which any vdwinteraction is neglected
!     vdweAbsurd  defines an energy limit above which any vdwinteraction is skiped
!
!
#include "tinker_macro.h"
module vdw
   implicit none
   integer nvdw,nvt,nvdwloc,nvdwbloc,nvdwlocnl
   integer nvdwblocloop
   integer nvdwlocnlb
   integer nvdwlocnlb_pair,nvdwlocnlb_pair1,nvdwlocnlb2_pair
   integer nshortvdwlocnlb2_pair
   integer nvdwclass
   integer, allocatable :: vdwlocnl(:)
   integer, pointer :: jvdw(:),ivdw(:),ired(:)
   integer, pointer :: jvdw_c(:)
   integer, pointer :: ivt(:),jvt(:),nbvdw(:)
   integer :: winjvdw,winjvdw_c,winivdw,winired,winivt
   integer :: winjvt,winnbvdw
   real(t_p), pointer:: radmin(:,:),epsilon(:,:)
   real(t_p), pointer:: radmin4(:,:),epsilon4(:,:)
   real(t_p), pointer:: radhbnd(:,:),epshbnd(:,:)
   real(t_p), pointer:: epsilon_c(:),epsilon4_c(:)
   real(t_p), pointer:: radmin_c(:),radmin4_c(:)
   real(t_p), pointer:: kred(:)
   integer winradmin,winradmin4
   integer winepsilon,winepsilon4
   integer winradmin_c,winradmin4_c
   integer winepsilon_c,winepsilon4_c
   integer winkred,winradhbnd,winepshbnd
   logical skipvdw12
   real(t_p),parameter:: vdw_lcut2=0.5**2
   real(t_p),parameter:: vdweAbsurd=2d2
end

module vdw_locArray
   implicit none
   real(t_p),allocatable,target :: xred(:),xredc(:)
   real(t_p),allocatable,target :: yred(:),yredc(:)
   real(t_p),allocatable,target :: zred(:),zredc(:)
   integer  ,allocatable,target :: loc_ired(:)
   real(t_p),allocatable,target :: loc_kred(:)
end module
