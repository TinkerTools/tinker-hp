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
c     der     repulsion Cartesian coordinate derivatives
c     dedsp   dispersion Cartesian coordinate derivatives
c     dedsprec   reciprocal dispersion Cartesian coordinate derivatives
c     dect    charge transfer Cartesian coordinate derivatives
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
c     Lambda-dynamics derivatives
c
c     delambda           hamiltonian derivative with respect to lambda (to be sent to colvar)
c     delambdae          hamiltonian derivative with respect to elambda
c     delambdav          hamiltonian derivative with respect to vlambda
c     delambdaesave      stored hamiltonian derivative with respect to elambda
c     delambdavsave      stored hamiltonian derivative with respect to vlambda
c     dlambdaelambda     derivative of elambda with respect to lambda
c     dlambdavlambda     derivative of vlambda with respect to lambda     
c
c     Orthogonal Space Random Walk - note x stands for Cartesian coordinates
c     dxdelambda         hamiltonian double derivative with respect to x and lambda (to be sent to colvar)
c     dxdelambdae        hamiltonian double derivative with respect to x and elambda (electrostatic interactions)
c     dxdelambdav        hamiltonian double derivative with respect to x and vlambda (vdw interactions)
c     d2edlambda2         hamiltonian double derivative with respect to lambda (to be sent to colvar)
c     d2edlambdae2        hamiltonian double derivative with respect to elambda (electrostatic interactions)
c     d2edlambdav2        hamiltonian double derivative with respect to vlambda (vdw interactions)
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
      real*8, allocatable :: der(:,:),dedsp(:,:),dect(:,:)
      real*8, allocatable :: dedsprec(:,:)
      real*8, allocatable :: dem(:,:),dep(:,:)
      real*8, allocatable :: deg(:,:),dex(:,:)
      real*8, allocatable :: decrec(:,:),demrec(:,:),deprec(:,:)
      real*8, allocatable :: debond(:,:),desave(:,:)
      real*8, allocatable :: desmd(:,:)
      real*8 :: delambda,delambdae,delambdav
      real*8 :: delambdaesave,delambdavsave
      real*8 :: d2edlambda2,d2edlambdae2,d2edlambdav2
      real*8 :: dlambdaelambda, dlambdavlambda
      real*8, allocatable  :: dxdelambda(:,:)
      real*8, allocatable :: dxdelambdae(:,:), dxdelambdav(:,:)
      logical dotstgrad
      save

      contains

      subroutine resetForcesRec
      implicit none
        if(allocated(demrec)) demrec = 0.0d0
        if(allocated(decrec)) decrec = 0.0d0
        if(allocated(deprec)) deprec = 0.0d0
        if(allocated(dedsprec)) dedsprec = 0.0d0
      end subroutine resetForcesRec

      end
