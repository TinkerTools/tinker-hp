c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ####################################################################
c     ##                                                                ##
c     ##  module mpole  --  multipole components for current structure  ##
c     ##                                                                ##
c     ####################################################################
c
c
c     maxpole   max components (monopole=1,dipole=4,quadrupole=13)
c
c     pole      multipole values for each site in the local frame
c     winpole    window object corresponding to pole
c     pole_orig  original multipole values for each site in the local frame (lambdadyn)
c     winpole_orig    window object corresponding to pole_orig
c     rpole     multipoles rotated to the global coordinate system
c     mono0     original atomic monopole values for charge flux
c     winmono0    window object corresponding to mono0
c
c     npole     total number of multipole sites in the system
c     npoleloc  local number of multipole sites in the system
c     npolelocloop  first multiple of 16 after npoleloc (help vecto)
c     npolebloc  local+neighbor number of multipole sites in the system
c     npoleblocloop  first multiple of 16 after npolebloc (help vecto)
c     npolelocnl  localnl number of multipole sites in the system
c     npolelocnlloop  first multiple of 16 after npolelocnl (help vecto)
c     npolerecloc  local number of reciprocal multipole sites in the system
c     npolerecloc_save  save local number of reciprocal multipole sites in the system
c     npolereclocloop  first multiple of 16 after npolerecloc (help vecto)
c
c     poleloc   global-local correspondance for multipoles
c     polelocnl global-localnl correspondace for multipoles
c     polerecloc   global-local correspondance for reciprocal multipoles
c     nbpole    number of multipoles before each atom
c
c     ipole     number of the atom for each multipole site
c     winipole    window object corresponding to ipole
c     polsiz    number of multipole components at each atom
c     winpolsiz    window object corresponding to polsiz
c     pollist   multipole site for each atom (0=no multipole)
c     winpollist    window object corresponding to pollist
c     zaxis     number of the z-axis defining atom for each site
c     xaxis     number of the x-axis defining atom for each site
c     yaxis     number of the y-axis defining atom for each site
c     polaxe    local axis type for each multipole site
c     ipolaxe   local axis type for each multipole site (integer)
c     winpolaxe    window object corresponding to polaxe
c     winipolaxe    window object corresponding to ipolaxe
c
c     alphapen  alpha parameters for sibfa like charge penetration
c     winalphapen    window object corresponding to alphapen
c     betapen   beta parameters for sibfa like charge penetration
c     winbetapen    window object corresponding to betapen
c     gammapen  gamma parameters for sibfa like charge penetration
c     wingammapen    window object corresponding to gammapen
c
c     npolelocnlb First multiple of BLOCK_SIZE after npolelocnl
c     npolelocnlb_pair  total number of electrostatics pair blocks interaction
c     npolelocnlb2_pair  total number of electrostatics pair blocks interaction from C2 nblist
c     nshortpolelocnlb2_pair  total number of electrostatics pair blocks interaction in short range interaction list
c
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
      integer  ,allocatable,target::poleloc(:),polelocnl(:)
     &         ,polerecloc(:)
      integer  ,pointer :: nbpole(:)
      !TODO [xyz]axis has been turned into an allocatable Investiguate
      integer  ,pointer :: zaxis(:),xaxis(:),yaxis(:)
      real(t_p),allocatable :: rpole(:,:)
      real(t_p),pointer :: pole(:,:), pole_orig(:,:), mono0(:)
      real(t_p),pointer :: alphapen(:),betapen(:),gammapen(:)
      integer  ,pointer :: ipolaxe(:)
      character*8,pointer :: polaxe(:)
      integer   winipole,winpolsiz,winpollist,winpole,winpole_orig
     &         ,winzaxis,winxaxis,winyaxis,winpolaxe,winipolaxe,winmono0
     &         ,winnbpole,winalphapen,winbetapen,wingammapen

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
      real(t_p),allocatable,target:: r2Work1(:,:),r2Work2(:,:)
     &         ,r2Work3(:,:),r2Work4(:,:)
      end module
