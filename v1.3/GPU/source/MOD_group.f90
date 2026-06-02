!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     #################################################################
!     ##                                                             ##
!     ##  module group  --  partitioning of system into atom groups  ##
!     ##                                                             ##
!     #################################################################
!
!
!     grpmass     total mass of all the atoms in each group
!     wgrp        weight for each set of group-group interactions
!     ngrp        total number of atom groups in the system
!     kgrp        contiguous list of the atoms in each group
!     igrp        first and last atom of each group in the list
!     grplist     number of the group to which each atom belongs
!     use_group   flag to use partitioning of system into groups
!     use_intra   flag to include only intragroup interactions
!     use_inter   flag to include only intergroup interactions
!     existScaledInter flag to determine whether interactions between groups are scaled
!
!     natgroup       number of atoms involved in a group
!     globglobgroup  associated indexes
!     loclocgroup    global-group correspondance
!
!
!     nlocatgroup    number of local atoms involved in a group
!     globgroup      associated group indexes
!     locgroup       global group - local group correspondance
!
!
!     npolegroup     number of multipoles involved in a group
!     ipolegroup     multipoles-atoms acorrespondance in the group frame
!     pollistgroup   atoms-multipoles correspondance in the group frame
!     npolelocgroup  number of local multipoles involved in a group
!     globpolegroup  global/group - global/global correspondance for the multipoles
!     poleglobgroup  local-global correspondance for the multipoles in the group frame
!     polelocgroup   global-local correspondance in the group frame
!
!
!     domlengroup     number of atoms in the group per domain
!     bufbeggroup    "bufbeg" equivalent in the group frame
!     domlenpolegroup number of multipoles in the group per domain
!     bufbegpolegroup "bufbegpole" equivalent in the group frame
!     uindgroup,uinpgroup  induced dipoles associated to the group
!     epgroup         polarization energy associated to the group
!     depgroup        polarization energy derivatives associated to the group
!     vir_group       virial associated to the polarization energy of the group
!
!
#include "tinker_macro.h"

module group
   use sizes
   implicit none
   integer    ngrp
   integer  ,allocatable :: kgrp(:),grplist(:)
   integer  ,allocatable :: igrp(:,:)
   real(t_p),allocatable :: grpmass(:),wgrp(:,:)
   real(t_p),allocatable,protected:: wgrp0(:,:)

   logical   use_group,use_intra,use_inter,existScaledInter&
      &,use_group_polar,use_group_mpole
   integer   natgroup,nlocatgroup,npolegroup,npolelocgroup
   integer  ,allocatable :: globglobgroup(:),loclocgroup(:)&
      &,globgroup(:),locgroup(:),globpolegroup(:)&
      &,poleglobgroup(:),polelocgroup(:),domlengroup(:)&
      &,ipolegroup(:),domlenpolegroup(:),bufbegpolegroup(:)&
      &,bufbeggroup(:),pollistgroup(:)
   real(t_p),allocatable :: uindgroup(:,:),uinpgroup(:,:)
   mdyn_rtyp,allocatable :: depgroup(:,:),demgroup(:,:)
   ener_rtyp epgroup, emgroup
   real(r_p) vir_group(3,3)

   integer :: n_uscale_group, n_dpscale_group, n_dpuscale_group
   integer :: n_mscale_group
   integer  ,allocatable:: ucorrect_ik_group(:)&
      &, dpcorrect_ik_group(:), dpucorrect_ik_group(:)&
      &, mcorrect_ik_group(:,:)
   real(t_p),allocatable:: ucorrect_scale_group(:)&
      &, dpcorrect_scale_group(:)&
      &, dpucorrect_scale_group(:), mcorrect_scale_group(:)


contains

   subroutine save_wgrp
      implicit none
      if(.not. allocated(wgrp0)) then
         allocate(wgrp0(ngrp+1,ngrp+1))
      end if
      wgrp0 = wgrp
   end subroutine

   subroutine load_wgrp
      implicit none
      wgrp = wgrp0
!$acc update device(wgrp) async
   end subroutine

end
