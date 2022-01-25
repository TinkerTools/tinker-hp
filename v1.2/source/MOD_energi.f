c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #################################################################
c     ##                                                             ##
c     ##  module energi  --  individual potential energy components  ##
c     ##                                                             ##
c     #################################################################
c
c
c     esum   total potential energy of the system
c     eb     bond stretch potential energy of the system
c     ea     angle bend potential energy of the system
c     eba    stretch-bend potential energy of the system
c     eub    Urey-Bradley potential energy of the system
c     eaa    angle-angle potential energy of the system
c     eopb   out-of-plane bend potential energy of the system
c     eopd   out-of-plane distance potential energy of the system
c     eid    improper dihedral potential energy of the system
c     eit    improper torsion potential energy of the system
c     et     torsional potential energy of the system
c     ept    pi-orbital torsion potential energy of the system
c     eat    angle-torsion potential energy of the system
c     ebt    stretch-torsion potential energy of the system
c     ett    torsion-torsion potential energy of the system
c     ev     van der Waals potential energy of the system
c     er     Pauli repulsion potential energy of the system
c     edsp   dispersion potential energy of the system
c     ec     charge-charge potential energy of the system
c     em     atomic multipole potential energy of the system
c     ep     polarization potential energy of the system
c     ect    charge transfer potential energy of the system
c     eg     geometric restraint potential energy of the system
c     ex     extra term potential energy of the system
c     esave  stored potential energy of the system
c     ensmd  extra term smd potential energy of the system
c
c
      module energi
      implicit none
      real*8 esum,eb,ea,eba
      real*8 eub,eaa,eopb,eopd
      real*8 eid,eit,et,ept
      real*8 eat,ebt,ett,ev,ec
      real*8 em,ep
      real*8 er,edsp,ect
      real*8 eg,ex
      real*8 esave
      real*8 ensmd
      save
      end
