!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ##################################################################
!     ##                                                              ##
!     ##  module mlpot  --                                              ##
!     ##                                                              ##
!     ##################################################################
!
!
!     MLpot   machine learning potential name
!     s_save
!     MLpotcut
!     xyz_c
!
#include "tinker_macro.h"
module mlpot
   use, intrinsic :: iso_c_binding

   integer, private :: s_save=0
   integer :: naml, namloc,namloc_strict, ml_port=0, gc_stride=1000
   integer :: ml_embedding_mode=1
   integer :: ml_pos_buffer = 0
   integer, allocatable :: list(:), iaml(:)
   integer, allocatable, target :: amloc(:)
   integer(c_int32_t), allocatable :: edge_src(:), edge_dst(:)
   integer(c_int), allocatable :: pair_rank(:)
   integer(c_int32_t), allocatable :: alch_group(:)
   integer(c_int32_t), allocatable :: atomic_mlpot(:), trueat_mlpot(:)
   integer, allocatable, target :: blist1(:), blist2(:)
   logical :: ml_resources_initialized=.FALSE.
   logical :: ml_dev_active=.FALSE.
   logical :: use_dmts=.FALSE.
   logical, allocatable:: laml(:), grpmllist(:)
   real(t_p) :: mlpotcut, mlpot_rfield
   real(t_p) :: model_cutoff, model_cutoff2
   real(t_p) :: mlpotscale=1.0
   real(c_float) cell_a(3,3), vir_ml(9)
   real(c_float),allocatable,target:: xyz_c(:)
   real(c_float),allocatable:: d2(:), dp(:), aemlp(:)
   real(c_float),allocatable:: d_mlpot(:)
   real(c_float) :: dedle_ml(1), dedlv_ml(1)
   real(c_float) :: mlqtot, mlqligand
   character*240  model_file
   character*240  model_file2
   character*240  model_file_cur
   integer, parameter :: ModelSlots = 2
   character(len=240) :: model_identity(ModelSlots) = ''
   integer :: model_padding(ModelSlots) = 0

   interface
      function init_ml_ressources(rank, devID, model_file, debug, port, use_lambda, qtot, qligand, gc_stride, nat, species, nligand, index_ligand) bind(C)
         import c_int32_t, c_float, c_char
         integer(c_int32_t) :: init_ml_ressources
         integer(c_int32_t), value :: rank, devID, debug, port, use_lambda, gc_stride, nat, nligand
         real(c_float), value :: qtot, qligand
         integer(c_int32_t), dimension(*) :: species, index_ligand
         character(kind=c_char), dimension(*) :: model_file
      end function

      function init_dmts(model_file2) bind(C)
        import c_int32_t, c_char
         integer(c_int32_t) :: init_dmts
         character(kind=c_char), dimension(*), intent(in) :: model_file2
      end function

      function get_model_cutoff(model_file, cutoff) bind(C)
         import c_int32_t, c_float, c_char
         integer(c_int32_t) :: get_model_cutoff
         character(kind=c_char), dimension(*), intent(in) :: model_file
         real(c_float), intent(out) :: cutoff
      end function


      function ml_models(coord, energies, gradients, vir, cell, species, l1, l2, &
         &dist, dv, trueat, natm, model_padding, natm_full, dograd, use_lambda, &
         &elambda, vlambda, alch_group, dedle, dedlv, model_file_cur) bind(C)
         import c_float, c_int32_t,c_char
         integer(c_int32_t) :: ml_models
         real(c_float), dimension(*) :: coord, energies, gradients, vir, cell, dist, dv
         integer(c_int32_t), dimension(*) :: species, l1, l2, alch_group, trueat
         integer(c_int32_t), value :: natm, model_padding, natm_full, dograd, use_lambda
         real(c_float), value :: elambda, vlambda
         real(c_float), dimension(*) :: dedle, dedlv
         character(kind=c_char), dimension(*) :: model_file_cur
      end function

   end interface

   interface
      subroutine set_embedding_weights
      end subroutine
      subroutine init_build_ml_bond_list
      end subroutine
      subroutine build_ml_bond_list
      end subroutine
   end interface

contains

integer function find_model_slot(model_name)
   character(*), intent(in) :: model_name
   integer i

   find_model_slot = 0
   do i = 1, ModelSlots
      if (trim(model_identity(i)) == trim(model_name)) then
         find_model_slot = i
         return
      end if
   end do
end function find_model_slot

subroutine realloc_position_buffer(rsize,npairs)
   integer,intent(in):: rsize,npairs
   integer :: pair_capacity

   if (rsize.gt.s_save) then
      if (allocated(xyz_c)) then
!$acc exit data delete(xyz_c,atomic_mlpot,aemlp,d_mlpot,trueat_mlpot,alch_group)
         deallocate(xyz_c,atomic_mlpot,aemlp,d_mlpot,trueat_mlpot,alch_group)
      end if
      if (s_save == 0) then 
         s_save = rsize
      else
         s_save = rsize + ml_pos_buffer
      end if
      !write(*,*) 'realloc_position_buffer:    s_save = ',s_save,' rsize = ',rsize

      allocate(xyz_c(s_save*3),atomic_mlpot(s_save),aemlp(s_save)&
         &,d_mlpot(3*s_save),trueat_mlpot(s_save),alch_group(s_save))
!$acc enter data create(xyz_c,atomic_mlpot,aemlp,d_mlpot,trueat_mlpot,alch_group)
      ! s_save = rsize
   end if
   if (.not.allocated(pair_rank) .or. npairs.gt.size(pair_rank)) then
      if (allocated(d2)) then
!$acc exit data delete(d2,dp,edge_src,edge_dst,pair_rank)
         deallocate(d2,dp,edge_src,edge_dst,pair_rank)
      end if
      pair_capacity = max(1,ceiling(1.01*real(npairs,8)))
      allocate(pair_rank(pair_capacity))
      allocate(d2(2*pair_capacity),dp(6*pair_capacity))
      allocate(edge_src(2*pair_capacity),edge_dst(2*pair_capacity))
!$acc enter data create(d2,dp,edge_src,edge_dst,pair_rank)
   end if

end subroutine realloc_position_buffer

subroutine release_mlpot_buffers
   if (allocated(xyz_c)) then
!$acc exit data delete(xyz_c,atomic_mlpot,aemlp,d_mlpot,trueat_mlpot,alch_group)
      deallocate(xyz_c,atomic_mlpot,aemlp,d_mlpot,trueat_mlpot,alch_group)
   end if
   if (allocated(d2)) then
!$acc exit data delete(d2,dp,edge_src,edge_dst,pair_rank)
      deallocate(d2,dp,edge_src,edge_dst,pair_rank)
   end if
   s_save = 0
   model_identity = ''
   model_padding = 0
   if (ml_dev_active) then
!$acc exit data delete(dedle_ml,dedlv_ml,vir_ml)
      ml_dev_active = .false.
   end if
end subroutine release_mlpot_buffers

end module
