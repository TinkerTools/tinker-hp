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
!
module replicas
   implicit none
   integer :: nreps !<number of replicas
   integer :: rank_reploc !<index of the local replica
   integer :: COMM_ROOT2ROOT !<root to roor communicator
   logical :: use_reps=.FALSE. !<logical flag to determine if multiple replicas are used
   real*8, allocatable :: lambdastart(:) !<array of initial lambda values (lambda dynamics)
   save
end module replicas
