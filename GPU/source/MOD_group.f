c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #################################################################
c     ##                                                             ##
c     ##  module group  --  partitioning of system into atom groups  ##
c     ##                                                             ##
c     #################################################################
c
c
c     grpmass     total mass of all the atoms in each group
c     wgrp        weight for each set of group-group interactions
c     ngrp        total number of atom groups in the system
c     kgrp        contiguous list of the atoms in each group
c     igrp        first and last atom of each group in the list
c     grplist     number of the group to which each atom belongs
c     use_group   flag to use partitioning of system into groups
c     use_intra   flag to include only intragroup interactions
c     use_inter   flag to include only intergroup interactions
c     existScaledInter flag to determine whether interactions between groups are scaled
c
c     natgroup       number of atoms involved in a group
c     globglobgroup  associated indexes
c     loclocgroup    global-group correspondance
c
c
c     nlocatgroup    number of local atoms involved in a group
c     globgroup      associated group indexes
c     locgroup       global group - local group correspondance
c
c
c     npolegroup     number of multipoles involved in a group
c     ipolegroup     multipoles-atoms acorrespondance in the group frame
c     pollistgroup   atoms-multipoles correspondance in the group frame
c     npolelocgroup  number of local multipoles involved in a group
c     globpolegroup  global/group - global/global correspondance for the multipoles
c     poleglobgroup  local-global correspondance for the multipoles in the group frame
c     polelocgroup   global-local correspondance in the group frame
c
c
c     domlengroup     number of atoms in the group per domain
c     bufbeggroup    "bufbeg" equivalent in the group frame
c     domlenpolegroup number of multipoles in the group per domain
c     bufbegpolegroup "bufbegpole" equivalent in the group frame
c     uindgroup,uinpgroup  induced dipoles associated to the group
c     epgroup         polarization energy associated to the group
c     depgroup        polarization energy derivatives associated to the group
c     vir_group       virial associated to the polarization energy of the group
c
c
#include "tinker_macro.h"

      module group
      use sizes
      implicit none
      integer    ngrp
      integer  ,allocatable :: kgrp(:),grplist(:)
      integer  ,allocatable :: igrp(:,:)
      real(t_p),allocatable :: grpmass(:),wgrp(:,:)
      real(t_p),allocatable,protected:: wgrp0(:,:)

      logical   use_group,use_intra,use_inter,existScaledInter
     &         ,use_group_polar,use_group_mpole
      integer   natgroup,nlocatgroup,npolegroup,npolelocgroup
      integer  ,allocatable :: globglobgroup(:),loclocgroup(:)
     &         ,globgroup(:),locgroup(:),globpolegroup(:)
     &         ,poleglobgroup(:),polelocgroup(:),domlengroup(:)
     &         ,ipolegroup(:),domlenpolegroup(:),bufbegpolegroup(:)
     &         ,bufbeggroup(:),pollistgroup(:)
      real(t_p),allocatable :: uindgroup(:,:),uinpgroup(:,:)
      mdyn_rtyp,allocatable :: depgroup(:,:),demgroup(:,:)
      ener_rtyp epgroup, emgroup
      real(r_p) vir_group(3,3)

      integer :: n_uscale_group, n_dpscale_group, n_dpuscale_group
      integer :: n_mscale_group
      integer  ,allocatable:: ucorrect_ik_group(:)
     &         , dpcorrect_ik_group(:), dpucorrect_ik_group(:)
     &         , mcorrect_ik_group(:,:)
      real(t_p),allocatable:: ucorrect_scale_group(:)
     &         , dpcorrect_scale_group(:)
     &         , dpucorrect_scale_group(:), mcorrect_scale_group(:)

      interface
        subroutine efld0_group(nrhs,ef)
          integer  , intent(in)    :: nrhs
          real(t_p), intent(inout) :: ef(:,:,:)
        end subroutine efld0_group
        subroutine efld0_group_correct_scaling(nrhs,ef)
          integer  , intent(in)    :: nrhs
          real(t_p), intent(inout) :: ef(:,:,:)
        end subroutine efld0_group_correct_scaling
        subroutine commfieldfull(nrhs,ef)
          integer  , intent(in)    :: nrhs
          real(t_p), intent(inout) :: ef(:,:,:)
        end subroutine commfieldfull
        subroutine commdirdirfull(nrhs,rule,mu,reqrec,reqsend)
          integer, intent(in) :: nrhs,rule
          integer, intent(inout) :: reqrec(nproc),reqsend(nproc)
          real(t_p), intent(inout) ::  mu(:,:,:)
        end subroutine commdirdirfull
        subroutine inducepcg_group(nrhs,precnd,ef,mu)
          integer  ,intent(in)   :: nrhs
          logical  ,intent(in)   :: precnd
          real(t_p),intent(in)   :: ef (:,:,:)
          real(t_p),intent(inout):: mu (:,:,:)
        end subroutine inducepcg_group
        subroutine tmatxb_group(nrhs,dodiag,mu,efi)
          integer, intent(in) ::  nrhs
          logical, intent(in) ::  dodiag
          real(t_p), intent(in) ::  mu(:,:,:)
          real(t_p), intent(inout) ::  efi(:,:,:)
        end subroutine tmatxb_group
        subroutine tmatxb_correct_interactions_group(nrhs,mu,efi)
          integer, intent(in) ::  nrhs
          real(t_p), intent(in) ::  mu(:,:,:)
          real(t_p), intent(inout) ::  efi(:,:,:)
        end subroutine tmatxb_correct_interactions_group
      end interface

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
