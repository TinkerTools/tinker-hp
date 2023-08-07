c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #####################################################################
c     ##                                                                 ##
c     ##  module replicass   --  replicas variables                      ##
c     ##                                                                 ##
c     #####################################################################
c
c     nreps: number of replicas
c     rank_reploc: index of the local replica
c     COMM_ROOT2ROOT: root to roor communicator
c     use_repls: logical flag to determine if multiple replicas are used
c
c
      module replicas
      implicit none
      integer :: nreps
      integer :: rank_reploc
      integer :: COMM_ROOT2ROOT
      logical :: use_reps=.FALSE.
      real*8, allocatable :: lambdastart(:)
      save
      end module replicas
