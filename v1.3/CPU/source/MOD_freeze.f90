!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     #####################################################################
!     ##                                                                 ##
!     ##  module freeze  --  definition of holonomic RATTLE constraints  ##
!     ##                                                                 ##
!     #####################################################################
!
!
module freeze
   implicit none
   integer :: nrat !<number of holonomic distance constraints to apply
   integer :: nratloc !<number of local constraints
   integer :: nratbloc !<number of local constraints + neighboring constraints
   integer :: winirat !<window object corresponding to irat
   integer :: winiratx !<window object corresponding to iratx
   integer :: winratx !<window object corresponding to ratx
   integer :: winkratx !<window object corresponding to kratx
   integer :: winkrat !<window object corresponding to krat
   integer :: winratimage !<window object corresponding to ratimage
   integer, pointer :: irat(:,:) !<atom numbers of atoms in a holonomic constraint
   integer, pointer :: iratx(:) !<group number of group in a holonomic constraint
   integer, pointer :: kratx(:) !<spatial constraint type (1=plane, 2=line, 3=point)
   integer, allocatable :: buflenrat1(:) !<length of buffer of constraints to send at each rattle iteration
   integer, allocatable :: buflenrat2(:) !<length of buffer of constraints to receive at each rattle iteration
   integer, allocatable :: bufbegrat1(:) !<first indexes of buffer of constraints to send at each rattle iteration
   integer, allocatable :: bufbegrat2(:) !<first indexes of buffer of constraints to receive at each rattle iteration
   integer, allocatable :: bufrat1(:) !<list of indexes of buffer of constraints to send at each rattle iteration
   integer, allocatable :: bufrat2(:) !<list of indexes of buffer of constraints to receive at each rattle iteration
   logical, pointer :: ratimage(:) !<flag to use minimum image for holonomic constraint
   logical :: use_rattle !<logical flag to set use of holonomic contraints
   real*8, pointer :: krat(:) !<ideal distance value for holonomic constraint
   real*8 :: rateps !<convergence tolerance for holonomic constraints
   save
end
