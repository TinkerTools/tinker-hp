c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ############################################################
c     ##                                                        ##
c     ##  bound.i  --  control of periodic boundary conditions  ##
c     ##                                                        ##
c     ############################################################
c
c
c     polycut       cutoff distance for infinite polymer nonbonds
c     polycut2      square of infinite polymer nonbond cutoff
c     use_bounds    flag to use periodic boundary conditions
c     use_replica   flag to use replicates for periodic system
c     use_polymer   flag to mark presence of infinite polymer
c
c
      real*8 polycut,polycut2
      logical use_bounds
      logical use_replica
      logical use_polymer
      common /bound/ use_bounds,use_replica,use_polymer
