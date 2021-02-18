c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ####################################################################
c     ##                                                                ##
c     ##  module deriv  --  Cartesian coordinate derivative components  ##
c     ##                                                                ##
c     ####################################################################
c
c
c     desum   total energy Cartesian coordinate derivatives
c     deb     bond stretch Cartesian coordinate derivatives
c     dea     angle bend Cartesian coordinate derivatives
c     deba    stretch-bend Cartesian coordinate derivatives
c     deub    Urey-Bradley Cartesian coordinate derivatives
c     deaa    angle-angle Cartesian coordinate derivatives
c     deopb   out-of-plane bend Cartesian coordinate derivatives
c     deopd   out-of-plane distance Cartesian coordinate derivatives
c     deid    improper dihedral Cartesian coordinate derivatives
c     deit    improper torsion Cartesian coordinate derivatives
c     det     torsional Cartesian coordinate derivatives
c     dept    pi-orbital torsion Cartesian coordinate derivatives
c     debt    stretch-torsion Cartesian coordinate derivatives
c     deat    angle-torsion Cartesian coordinate derivatives
c     dett    torsion-torsion Cartesian coordinate derivatives
c     dev     van der Waals Cartesian coordinate derivatives
c     dec     charge-charge Cartesian coordinate derivatives
c     decrec  reciprocal charge-charge Cartesian coordinate derivatives
c     dem     multipole Cartesian coordinate derivatives
c     demrec  reciprocal multipole Cartesian coordinate derivatives
c     dep     polarization Cartesian coordinate derivatives
c     dep     reciprocal polarization Cartesian coordinate derivatives
c     deg     geometric restraint Cartesian coordinate derivatives
c     dex     extra energy term Cartesian coordinate derivatives
c     desave  stored Cartesian coordinate derivatives
c     desmd   extra smd energy term Cartesian coordinate derivatives
c
c     dotstgrad : flag when the main program is testgrad (communication
c      of the forces one by one)
c
c
      module deriv
      implicit none
      real*8, allocatable :: desum(:,:),deb(:,:),dea(:,:),deba(:,:)
      real*8, allocatable :: deub(:,:),deaa(:,:),deopb(:,:),deopd(:,:)
      real*8, allocatable :: deid(:,:),det(:,:),dept(:,:),deit(:,:)
      real*8, allocatable :: deat(:,:)
      real*8, allocatable :: debt(:,:),dett(:,:),dev(:,:),dec(:,:)
      real*8, allocatable :: dem(:,:),dep(:,:)
      real*8, allocatable :: deg(:,:),dex(:,:)
      real*8, allocatable :: decrec(:,:),demrec(:,:),deprec(:,:)
      real*8, allocatable :: debond(:,:),desave(:,:)
      real*8, allocatable :: desmd(:,:)
      logical dotstgrad
      save
      end
