!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
#include "tinker_macro.h"
module structure
   use config, only: config_t
   use status, only: status_t
   use tinheader, only: ti_p, re_p
   implicit none

   ! structure_t is value-copy: its allocatable components own independent host storage.
   type, public :: structure_t
      character(len=240) :: source_path = ''
      character(len=240) :: title = ''
      integer :: atom_count = 0
      integer :: max_neighbor_count = 0
      logical :: has_bounds = .false.
      real(re_p), allocatable :: x(:), y(:), z(:)
      integer, allocatable :: tag(:), atom_type(:)
      character(len=3), allocatable :: atom_name(:)
      integer, allocatable :: neighbor_count(:), connectivity(:,:)
      real(re_p) :: box_length(3) = 0.0_re_p
      real(ti_p) :: box_angle(3) = 0.0_ti_p
   end type structure_t

   public :: load_initial
   public :: read_frame
   interface
      module subroutine load_initial(configuration, structure_value, status_value)
         type(config_t), intent(in) :: configuration
         type(structure_t), intent(out) :: structure_value
         type(status_t), intent(out) :: status_value
      end subroutine load_initial

      module subroutine read_frame(unit, structure_value, status_value)
         integer, intent(in) :: unit
         type(structure_t), intent(inout) :: structure_value
         type(status_t), intent(out) :: status_value
      end subroutine read_frame
   end interface
end module structure
