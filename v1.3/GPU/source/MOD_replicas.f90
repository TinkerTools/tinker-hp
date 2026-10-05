!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     #####################################################################
!     ##                                                                 ##
!     ##  module replicass   --  replicas variables                      ##
!     ##                                                                 ##
!     #####################################################################
!
!     nreps: number of replicas
!     rank_reploc: index of the local replica
!     use_repls: logical flag to determine if multiple replicas are used
!     lambdastart: starting values of lambda for multiple replicas lambda dyn
!
!
#include "tinker_precision.h"
module replicas
   implicit none
   integer :: nreps
   integer :: rank_reploc
   logical :: use_reps
   real(t_p), allocatable :: lambdastart(:)
   save
end module replicas
