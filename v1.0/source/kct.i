c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #################################################################
c     ##                                                             ##
c     ##  kct.i  --  forcefield parameters for charge transfer       ##
c     ##                                                             ##
c     #################################################################
c
c     vdwct1 : donnor vdw radii for charge transfer
c     vdwct2 : acceptor vdw radii for charge transfer
c     hybrid : hybridation coefficients for charge transfer
c
c
      integer maxnct
      parameter (maxnct=2000)
      real*8 sibfact1,sibfact2
      real*8 hybrid,tas,tap,ma
      common /kct/ sibfact1(maxtyp),sibfact2(maxtyp),hybrid(2,maxtyp),
     $             tas(5,200),tap(5,200),ma(5,200)
