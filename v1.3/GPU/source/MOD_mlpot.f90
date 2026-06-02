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
   integer(c_int32_t), allocatable :: alch_group(:)
   integer(c_int32_t), allocatable :: atomic_mlpot(:), trueat_mlpot(:)
   integer, allocatable, target :: blist1(:), blist2(:)
   logical :: ml_resources_initialized=.FALSE.
   logical, allocatable:: laml(:), grpmllist(:)
   real(t_p) :: mlpotcut, mlpot_rfield
   real(t_p) :: mlpotscale=1.0
   real(c_float) cell_a(3,3), vir_ml(9)
   real(c_float),allocatable,target:: xyz_c(:)
   real(c_float),allocatable:: d2(:), dp(:), aemlp(:)
   real(c_float),allocatable:: d_mlpot(:)
   real(c_float) :: dedle_ml(1), dedlv_ml(1)
   real(c_float) :: mlqtot, mlqligand
   character*240  model_file

   interface
      function init_ml_ressources(rank, devID, model_file, debug, port, use_lambda, qtot, qligand, gc_stride, nat, species, nligand, index_ligand) bind(C)
         import c_int32_t, c_float, c_char
         integer(c_int32_t) :: init_ml_ressources
         integer(c_int32_t), value :: rank, devID, debug, port, use_lambda, gc_stride, nat, nligand
         real(c_float), value :: qtot, qligand
         integer(c_int32_t), dimension(*) :: species, index_ligand
         character(kind=c_char), dimension(*) :: model_file
      end function

      function ml_models(coord, energies, gradients, vir, cell, species, l1, l2, dist, dv, trueat, natm, npairs, natm_full, dograd, use_lambda, elambda, vlambda, alch_group, dedle, dedlv) bind(C)
         import c_float, c_int32_t
         integer(c_int32_t) :: ml_models
         real(c_float), dimension(*) :: coord, energies, gradients, vir, cell, dist, dv
         integer(c_int32_t), dimension(*) :: species, l1, l2, alch_group, trueat
         integer(c_int32_t), value :: natm, npairs, natm_full, dograd, use_lambda
         real(c_float), value :: elambda, vlambda
         real(c_float), dimension(*) :: dedle, dedlv
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

subroutine realloc_position_buffer(rsize,npairs)
   integer,intent(in):: rsize,npairs
   integer:: max_pairs

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
      write(*,*) 'realloc_position_buffer: s_save = ',s_save,' rsize = ',rsize

      allocate(xyz_c(s_save*3),atomic_mlpot(s_save),aemlp(s_save)&
         &,d_mlpot(3*s_save),trueat_mlpot(s_save),alch_group(s_save))
!$acc enter data create(xyz_c,atomic_mlpot,aemlp,d_mlpot,trueat_mlpot,alch_group)
      ! s_save = rsize
   end if
   if (2*npairs.gt.size(d2)) then
      if (allocated(d2)) then
!$acc exit data delete(d2,dp,edge_src,edge_dst)
         deallocate(d2,dp,edge_src,edge_dst)
      end if
      max_pairs = 2*ceiling(1.01*real(npairs,8))
      write(*,*) 'realloc_position_buffer: max_pairs = ',max_pairs,' npairs = ',2*npairs
      ! max_pairs = 2*npairs
      allocate(d2(max_pairs),dp(3*max_pairs))
      allocate(edge_src(max_pairs),edge_dst(max_pairs))
!$acc enter data create(d2,dp,edge_src,edge_dst)
   end if

end subroutine realloc_position_buffer

end module
