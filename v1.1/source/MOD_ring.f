c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #####################################################################
c     ##                                                                 ##
c     ##  module ring  --  number and location of small ring structures  ##
c     ##                                                                 ##
c     #####################################################################
c
c
c     nring3   total number of 3-membered rings in the system
c     iring3   numbers of the atoms involved in each 3-ring
c     nring4   total number of 4-membered rings in the system
c     iring4   numbers of the atoms involved in each 4-ring
c     nring5   total number of 5-membered rings in the system
c     iring5   numbers of the atoms involved in each 5-ring
c     nring6   total number of 6-membered rings in the system
c     iring6   numbers of the atoms involved in each 6-ring
c
c
      module ring
      use sizes
      implicit none
      integer nring3,iring3(3,maxring)
      integer nring4,iring4(4,maxring)
      integer nring5,iring5(5,maxring)
      integer nring6,iring6(6,maxring)
      save
      end
