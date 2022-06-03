c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ###################################################################
c     ##                                                               ##
c     ##  module output  --  control of coordinate output file format  ##
c     ##                                                               ##
c     ###################################################################
c
c
c     archive    logical flag to save structures in an archive
c     noversion  logical flag governing use of filename versions
c     overwrite  logical flag to overwrite intermediate files inplace
c     cyclesave  logical flag to mark use of numbered cycle files
c     coordtype  selects Cartesian, internal, rigid body or none
c     new_restart logical flag to write restart  in a new file
c     f_mdsave    logical flag to force md traj writing
c
c
      module output
      implicit none
      logical archive,noversion
      logical overwrite,cyclesave
      logical new_restart,f_mdsave
      character*9 coordtype
      end
