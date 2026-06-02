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
!     rateps       convergence tolerance for holonomic constraints
!     krat         ideal distance value for holonomic constraint
!     winkrat    window object corresponding to krat
!     nrat         number of holonomic distance constraints to apply
!     nratx        number of atom group holonomic constraints to apply
!     irat         atom numbers of atoms in a holonomic constraint
!     winirat    window object corresponding to irat
!     iratx        group number of group in a holonomic constraint
!     winiratx    window object corresponding to iratx
!     kratx        spatial constraint type (1=plane, 2=line, 3=point)
!     winkratx    window object corresponding to kratx
!     ratimage     flag to use minimum image for holonomic constraint
!     winratimage    window object corresponding to ratimage
!     use_rattle   logical flag to set use of holonomic contraints
!
!     nratloc      number of local constraints
!     nratbloc      number of local constraints + neighboring constraints
!
!     buflenrat1   length of buffer of constraints to send at each rattle iteration
!     buflenrat2   length of buffer of constraints to receive at each rattle iteration
!     buflenrat1   first indexes of buffer of constraints to send at each rattle iteration
!     buflenrat2   first indexes of buffer of constraints to receive at each rattle iteration
!     bufrat1   list of indexes of buffer of constraints to send at each rattle iteration
!     bufrat2   list of indexes of buffer of constraints to receive at each rattle iteration
!
!
#include "tinker_macro.h"
module freeze
   implicit none
   integer nrat,nratx
   !integer, pointer :: irat(:,:),iratx(:)
   !integer, pointer ::  kratx(:)
   integer,allocatable :: irat(:,:), iratx(:)
   integer,allocatable ::  kratx(:)
   integer :: winirat,winiratx,winratx,winkratx
   !DIR$ ATTRIBUTES ALIGN:64:: buflenrat1,buflenrat2
   integer, allocatable :: buflenrat1(:),buflenrat2(:)
   !DIR$ ATTRIBUTES ALIGN:64:: bufbegrat1,bufbegrat2
   integer, allocatable :: bufbegrat1(:),bufbegrat2(:)
   !DIR$ ATTRIBUTES ALIGN:64:: bufrat1,bufrat2
   integer, allocatable :: bufrat1(:),bufrat2(:)
   integer nratloc,nratbloc
   !   real(t_p), pointer :: krat(:)
   real(t_p),allocatable :: krat(:)
   real(t_p) rateps
   !   logical, pointer :: ratimage(:)
   logical, allocatable :: ratimage(:)
   integer :: winkrat,winratimage
   logical use_rattle
   integer,allocatable :: nratmol(:),iratmol(:,:)
   save
end
