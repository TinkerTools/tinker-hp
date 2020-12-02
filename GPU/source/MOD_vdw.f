c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ######################################################################
c     ##                                                                  ##
c     ##  module vdw  --  van der Waals parameters for current structure  ##
c     ##                                                                  ##
c     ######################################################################
c
c
c     radmin     minimum energy distance for each atom class pair
c     winradmin    window object corresponding to radmin
c     epsilon    well depth parameter for each atom class pair
c     winepsilon    window object corresponding to epsilon
c     radmin4    minimum energy distance for 1-4 interaction pairs
c     winradmin4    window object corresponding to radmin4
c     epsilon4   well depth parameter for 1-4 interaction pairs
c     winepsilon4    window object corresponding to epsilon4
c     radmin_c   same as radmin but with compress class
c     winradmin_c    window object corresponding to radmin_c
c     epsilon_c  same as epsilon but with compress class
c     radmin4_c  same as radmin4 but with compress class
c     epsilon4_c same as epsilon4 but with compress class
c     radhbnd    minimum energy distance for hydrogen bonding pairs
c     winradhbnd    window object corresponding to radhbnd
c     epshbnd    well depth parameter for hydrogen bonding pairs
c     winepshbnd    window object corresponding to epshbnd
c     kred       value of reduction factor parameter for each atom
c     winkred    window object corresponding to kred
c     ired       attached atom from which reduction factor is applied
c     winired    window object corresponding to ired
c     nvdw       total number van der Waals active sites in the system
c     nvdwloc    local number van der Waals active sites in the system
c     nvdwbloc   local+neighbors number van der Waals active sites in the system
c     nvdwblocloop   First multiple of 16 after nvdwbloc (help vectorization
c     nvdwclass  total number of van der Waals different atoms class inside the system
c     ivdw       number of the atom for each van der Waals active site
c     winivdw    window object corresponding to ivdw
c     jvdw       type or class index into vdw parameters for each atom
c     winjvdw    window object corresponding to jvdw
c     jvdw_c     renumbered type or class index into vdw parameters for each atom
c     nvt        number of distinct vdw types/classes in the system
c     ivt        type/class index for each distinct vdw type or class
c     jvt        frequency of each vdw type or class in the system
c     nbvdw      number of 'vdw' atoms before each atom
c     winnbvdw    window object corresponding to nbvdw
c     vdwlocnl   glob-locnl vdw correspondance
c     nvdwlocnlb First multiple of BLOCK_SIZE after nvdwlocnl
c     nvdwlocnlb_pair  total number of vdw pair blocks interaction
c     nvdwlocnlb2_pair  total number of vdw pair blocks interaction in C2 nblist
c     nshortvdwlocnlb2_pair  total number of vdw pair blocks interaction in short range interaction list
c     skipvdw12  switch to skip vdw 1-2 Interactions computation
c
c
#include "tinker_precision.h"
      module vdw
      implicit none
      integer nvdw,nvt,nvdwloc,nvdwbloc,nvdwlocnl
      integer nvdwblocloop
      integer nvdwlocnlb
      integer nvdwlocnlb_pair,nvdwlocnlb2_pair
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
      end

      module vdw_locArray
      implicit none
      real(t_p),allocatable,target :: xred(:),xredc(:)
      real(t_p),allocatable,target :: yred(:),yredc(:)
      real(t_p),allocatable,target :: zred(:),zredc(:)
      integer  ,allocatable,target :: loc_ired(:)
      real(t_p),allocatable,target :: loc_kred(:)
      end module
