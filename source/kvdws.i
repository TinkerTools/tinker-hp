c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ##################################################################
c     ##                                                              ##
c     ##  kvdws.i  --  forcefield parameters for van der Waals terms  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     rad      van der Waals radius parameter for each atom type
c     eps      van der Waals well depth parameter for each atom type
c     rad4     van der Waals radius parameter in 1-4 interactions
c     eps4     van der Waals well depth parameter in 1-4 interactions
c     reduct   van der Waals reduction factor for each atom type
c
c
      real*8 rad,eps
      real*8 rad4,eps4
      real*8 reduct
      common /kvdws/ rad(maxtyp),eps(maxtyp),rad4(maxtyp),eps4(maxtyp),
     &               reduct(maxtyp)
