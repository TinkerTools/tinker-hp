c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c
c     #################################################################
c     ##                                                             ##
c     ##  module  kct  --  forcefield parameters for charge transfer ##
c     ##                                                             ##
c     #################################################################
c
c     sibfact1  atom type indexed vdw radius for electron donnnor 
c     sibfact2  atom type indexed vdw radius for electron acceptor 
c     hybrid    atom type indexed lone pair hybridization coefficients
c     tas       atom type indexed coefficients to compute integrals 
c     tap       atom type indexed coefficients to compute integrals 
c     ma        atom type indexed coefficients to compute integrals
c     ae        atom type indexed ionization potential
c     ah        atom type indexed electronic affinity
c
c
      module kct
      use sizes
      implicit none
      real*8 sibfact1(maxtyp),sibfact2(maxtyp)
      real*8 hybrid(2,maxtyp)
      real*8 tas(100,100),tap(100,100),ma(100,100)
      real*8 ah(maxtyp),ae(maxtyp)
      save
      end
