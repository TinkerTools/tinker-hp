c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #####################################################################
c     ##                                                                 ##
c     ##  module atmlst    --  local geometry terms involving each atom  ##
c     ##                                                                 ##
c     #####################################################################
c
c
c     bndlist   list of the bond numbers involving each atom
c     winbndlist window object corresponding to bndlist 
c     anglist   list of the angle numbers centered on each atom
c     winanglist window object corresponding to anglist 
c     
c     bndglob   local - global bond  correspondance 
c     angleglob local - global angle correspondance 
c     torsglob  local - global torsion correspondance 
c     bitorsglob local - global bitorsion correspondance 
c     strbndglob local - global strech bending correspondance 
c     ureyglob  local - global urey bradley correspondance 
c     angangglob local - global angle correspondance 
c     opbendglob local - global out of plane bending correspondance 
c     opdistglob local - global out of plane distance correspondance 
c     impropglob local - global improper dihedral correspondance 
c     imptorglob local - global improper torsion correspondance 
c     pitorsglob local - global pi torsion correspondance 
c     strtorglob local - global strech torsion correspondance 
c     tortorglob local - global torsion torsion correspondance 
c     vdwglob   local - global vdw correspondance 
c     poleglob  local - global direct multipole correspondance 
c     polerecglob local - global reciprocal multipole correspondance 
c     chgglob  local - global direct charge correspondance 
c     chgrecglob local - global reciprocal charge correspondance 
c
c     molculeglob local - global molecule correspondance 
c     npfixglob local - global position restrains correspondance 
c     ndfixglob local - global distance restrains correspondance 
c     nafixglob local - global angle restrains correspondance 
c     ntfixglob local - global torsion restrains correspondance 
c     ngfixglob local - global group restrains correspondance 
c     nchirfixglob local - global chiral restrains correspondance 
c     ratglob local - global constrains correspondance 
c
c     chgglobnl  localnl - global direct charge correspondance 
c     vdwglobnl  localnl - global vdw correspondance 
c     poleglobbnl  localnl - global direct multipole correspondance 
c
#include "tinker_precision.h"
      module atmlst
#ifdef USE_NVSHMEM_CUDA
      use tinTypes,only: i2dDPC=>Int2dDevPointerContainer
#endif
      implicit none
      integer, pointer :: bndlist(:,:),anglist(:,:)
      integer winbndlist,winanglist
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
!DIR$ ATTRIBUTES ALIGN:64:: tortorglob
      integer, allocatable :: tortorglob(:)
!DIR$ ATTRIBUTES ALIGN:64:: vdwglob,poleglob,polerecglob
      integer, allocatable,target :: vdwglob(:),poleglob(:)
     &                    , polerecglob(:)
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
