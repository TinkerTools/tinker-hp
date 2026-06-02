!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ####################################################################
!     ##                                                                ##
!     ##  module mpole  --  multipole components for current structure  ##
!     ##                                                                ##
!     ####################################################################
!
!
!     maxpole   max components (monopole=1,dipole=4,quadrupole=13)
!
!     pole      multipole values for each site in the local frame
!     winpole    window object corresponding to pole
!     pole_orig  original multipole values for each site in the local frame (lambdadyn)
!     winpole_orig    window object corresponding to pole_orig
!     rpole     multipoles rotated to the global coordinate system
!     mono0     original atomic monopole values for charge flux
!     winmono0    window object corresponding to mono0
!
!     npole     total number of multipole sites in the system
!     npoleloc  local number of multipole sites in the system
!     npolelocloop  first multiple of 16 after npoleloc (help vecto)
!     npolebloc  local+neighbor number of multipole sites in the system
!     npoleblocloop  first multiple of 16 after npolebloc (help vecto)
!     npolelocnl  localnl number of multipole sites in the system
!     npolelocnlloop  first multiple of 16 after npolelocnl (help vecto)
!     npolerecloc  local number of reciprocal multipole sites in the system
!     npolerecloc_save  save local number of reciprocal multipole sites in the system
!     npolereclocloop  first multiple of 16 after npolerecloc (help vecto)
!
!     poleloc   global-local correspondance for multipoles
!     polelocnl global-localnl correspondace for multipoles
!     polerecloc   global-local correspondance for reciprocal multipoles
!     nbpole    number of multipoles before each atom
!
!     ipole     number of the atom for each multipole site
!     winipole    window object corresponding to ipole
!     polsiz    number of multipole components at each atom
!     winpolsiz    window object corresponding to polsiz
!     pollist   multipole site for each atom (0=no multipole)
!     winpollist    window object corresponding to pollist
!     zaxis     number of the z-axis defining atom for each site
!     xaxis     number of the x-axis defining atom for each site
!     yaxis     number of the y-axis defining atom for each site
!     polaxe    local axis type for each multipole site
!     ipolaxe   local axis type for each multipole site (integer)
!     winpolaxe    window object corresponding to polaxe
!     winipolaxe    window object corresponding to ipolaxe
!
!     alphapen  alpha parameters for sibfa like charge penetration
!     winalphapen    window object corresponding to alphapen
!     betapen   beta parameters for sibfa like charge penetration
!     winbetapen    window object corresponding to betapen
!     gammapen  gamma parameters for sibfa like charge penetration
!     wingammapen    window object corresponding to gammapen
!
!     npolelocnlb First multiple of BLOCK_SIZE after npolelocnl
!     npolelocnlb_pair  total number of electrostatics pair blocks interaction
!     npolelocnlb2_pair  total number of electrostatics pair blocks interaction from C2 nblist
!     nshortpolelocnlb2_pair  total number of electrostatics pair blocks interaction in short range interaction list
!
#include "tinker_macro.h"
module mpole
   !use iso_c_binding,only:c_ptr
   implicit none
   integer maxpole
   parameter (maxpole=13)
   integer npole,npoleloc,npolebloc,npolerecloc,npolelocnl
   integer npolelocnlb,npolelocnlb_pair,npolelocnlb2_pair
   integer nshortpolelocnlb2_pair
   integer npolerecloc_old
   integer npolelocloop, npolelocnlloop,npoleblocloop,npolereclocloop
   integer :: nZ_Onlyloc=0,nZ_Onlyglob
   real(r_p) vmxx,vmyy,vmzz
   real(r_p) vmxy,vmxz,vmyz
   integer  ,pointer :: ipole(:), polsiz(:), pollist(:), iglobpole(:)
   integer  ,allocatable,target::poleloc(:),polelocnl(:)&
      &,polerecloc(:)
   integer  ,pointer :: nbpole(:)
   !TODO [xyz]axis has been turned into an allocatable Investiguate
   integer  ,pointer :: zaxis(:),xaxis(:),yaxis(:)
   real(t_p),allocatable :: rpole(:,:)
   real(t_p),pointer :: pole(:,:), pole_orig(:,:), mono0(:)
   real(t_p),pointer :: alphapen(:),betapen(:),gammapen(:)
   integer  ,pointer :: ipolaxe(:)
   character*8,pointer :: polaxe(:)
   integer   winipole,winpolsiz,winpollist,winpole,winpole_orig&
      &,winzaxis,winxaxis,winyaxis,winpolaxe,winipolaxe,winmono0&
      &,winnbpole,winalphapen,winbetapen,wingammapen

   ! Axetyp enumeration
   enum,bind(C)
      enumerator ::Ax_None=0
      enumerator ::Ax_3_Fold=1
      enumerator ::Ax_Bisector=2
      enumerator ::Ax_Z_Bisect=4
      enumerator ::Ax_Z_Only=8
      enumerator ::Ax_Z_Then_X=16
   end enum
!$acc declare create(vmxx,vmxy,vmxz,vmyy,vmyz,vmzz)
end

module elec_wspace
   real(t_p),allocatable,target:: rWork1(:),rWork2(:),rWork3(:)
   real(t_p),allocatable,target:: r2Work1(:,:),r2Work2(:,:)&
      &,r2Work3(:,:),r2Work4(:,:)
end module
