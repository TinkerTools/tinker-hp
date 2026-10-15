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
   integer, protected, pointer :: rank_reploc => null()
   logical :: use_reps
   real(t_p), allocatable :: lambdastart(:)
   save
contains

   subroutine bind_replica_view(index)
      integer, target, intent(inout) :: index
      rank_reploc => index
   end subroutine bind_replica_view

   subroutine clear_replica_view()
      nullify(rank_reploc)
   end subroutine clear_replica_view

end module replicas
