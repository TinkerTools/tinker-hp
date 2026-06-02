!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     #####################################################################
!     ##                                                                 ##
!     ##  module atmlst    --  local geometry terms involving each atom  ##
!     ##                                                                 ##
!     #####################################################################
!
!
!     bndlist   list of the bond numbers involving each atom
!     winbndlist window object corresponding to bndlist
!     anglist   list of the angle numbers centered on each atom
!     winanglist window object corresponding to anglist
!     balist   numbers of the bonds comprising each angle
!     winbalist window object corresponding to balist
!
!     bndglob   local - global bond  correspondance
!     angleglob local - global angle correspondance
!     torsglob  local - global torsion correspondance
!     bitorsglob local - global bitorsion correspondance
!     strbndglob local - global strech bending correspondance
!     ureyglob  local - global urey bradley correspondance
!     angangglob local - global angle correspondance
!     angtorglob local - global angle-torsion correspondance
!     opbendglob local - global out of plane bending correspondance
!     opdistglob local - global out of plane distance correspondance
!     impropglob local - global improper dihedral correspondance
!     imptorglob local - global improper torsion correspondance
!     pitorsglob local - global pi torsion correspondance
!     strtorglob local - global strech torsion correspondance
!     tortorglob local - global torsion torsion correspondance
!     vdwglob   local - global vdw correspondance
!     dispglob   local - global dispersion correspondance
!     poleglob  local - global direct multipole correspondance
!     polerecglob local - global reciprocal multipole correspondance
!     chgglob  local - global direct charge correspondance
!     chgrecglob local - global reciprocal charge correspondance
!     disprecglob local - global reciprocal dispersion correspondance
!
!     molculeglob local - global molecule correspondance
!     npfixglob local - global position restrains correspondance
!     ndfixglob local - global distance restrains correspondance
!     nafixglob local - global angle restrains correspondance
!     ntfixglob local - global torsion restrains correspondance
!     ngfixglob local - global group restrains correspondance
!     nchirfixglob local - global chiral restrains correspondance
!     ratglob local - global constrains correspondance
!
!     chgglobnl  localnl - global direct charge correspondance
!     vdwglobnl  localnl - global vdw correspondance
!     dispglobnl  localnl - global dispersion correspondance
!     poleglobbnl  localnl - global direct multipole correspondance
!
#include "tinker_macro.h"
module atmlst
#ifdef USE_NVSHMEM_CUDA
   use tinTypes,only: i2dDPC=>Int2dDevPointerContainer
#endif
   implicit none
   integer, pointer :: bndlist(:,:),anglist(:,:),balist(:,:)
   integer winbndlist,winanglist,winbalist
!DIR$ ATTRIBUTES ALIGN:64:: bndglob,angleglob,torsglob
   integer, allocatable :: bndglob(:),angleglob(:),torsglob(:)
!DIR$ ATTRIBUTES ALIGN:64:: bitorsglob,strbndglob
   integer, allocatable :: bitorsglob(:),strbndglob(:)
!DIR$ ATTRIBUTES ALIGN:64:: ureyglob,angangglob
   integer, allocatable :: ureyglob(:),angangglob(:)
!DIR$ ATTRIBUTES ALIGN:64:: opbendglob,opdistglob
   integer, allocatable :: opbendglob(:),opdistglob(:)
!DIR$ ATTRIBUTES ALIGN:64:: impropglob,imptorglob
   integer, allocatable :: impropglob(:),imptorglob(:)
!DIR$ ATTRIBUTES ALIGN:64:: pitorsglob,strtorglob
   integer, allocatable :: pitorsglob(:),strtorglob(:)
!DIR$ ATTRIBUTES ALIGN:64:: angtorglob,tortorglob
   integer, allocatable :: angtorglob(:),tortorglob(:)
!DIR$ ATTRIBUTES ALIGN:64:: vdwglob,poleglob,polerecglob
   integer, allocatable,target :: vdwglob(:),poleglob(:)&
      &, polerecglob(:)&
      &, dispglob(:),disprecglob(:)
!DIR$ ATTRIBUTES ALIGN:64:: chgglob,chgrecglob
   integer, allocatable,target :: chgglob(:),chgrecglob(:)
!DIR$ ATTRIBUTES ALIGN:64:: molculeglob
   integer, allocatable :: molculeglob(:)
!DIR$ ATTRIBUTES ALIGN:64:: npfixglob,ndfixglob,nafixglob
   integer, allocatable :: npfixglob(:),ndfixglob(:),nafixglob(:)
!DIR$ ATTRIBUTES ALIGN:64:: ntfixglob,ngfixglob,nchirglob
   integer, allocatable :: ntfixglob(:),ngfixglob(:),nchirglob(:)
!DIR$ ATTRIBUTES ALIGN:64:: ratglob
   integer, allocatable :: ratglob(:)
!DIR$ ATTRIBUTES ALIGN:64:: chgglobnl,vdwglobnl,poleglobnl
   integer, allocatable :: chgglobnl(:),vdwglobnl(:),poleglobnl(:)
   integer, allocatable :: dispglobnl(:)
   integer, pointer :: AtomKind(:)

#ifdef USE_NVSHMEM_CUDA
   !d_*     device data type container for nvshmem feature
   !c_*       host data type container for nvshmem feature
   type(i2dDPC),device,pointer::d_bndlist(:)
   type(i2dDPC),   allocatable::c_bndlist(:)
   type(i2dDPC),device,pointer::d_anglist(:)
   type(i2dDPC),   allocatable::c_anglist(:)
#endif
end
