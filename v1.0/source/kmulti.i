c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #################################################################
c     ##                                                             ##
c     ##  kmulti.i  --  forcefield parameters for atomic multipoles  ##
c     ##                                                             ##
c     #################################################################
c
c
c     maxnmp   maximum number of atomic multipole parameter entries
c
c     multip   atomic monopole, dipole and quadrupole values
c     mpaxis   type of local axis definition for atomic multipoles
c     kmp      string of atom types for atomic multipoles
c
c
      integer maxnmp
      parameter (maxnmp=2000)
      real*8 multip,sibfacp
      character*8 mpaxis
      character*16 kmp
      common /kmulti/ multip(13,maxnmp),mpaxis(maxnmp),kmp(maxnmp),
     $                sibfacp(3,maxtyp)
