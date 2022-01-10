c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #################################################################################
c     ##                                                                             ##
c     ##  module dispersion  --  dispersion caracteristics for the current structure ##
c     ##                                                                             ##
c     #################################################################################
c
c     admp6 damping parameter for exponential coeff of the 1/(r6) term
c     admp8 damping parameter for exponential coeff of the 1/(r8) term
c     admp8 damping parameter for exponential coeff of the 1/(r10) term
c     bdmp  coeff used to compute reduced distance
c     c6disp  1/r6 coefficient
c     c8disp  1/r8 coefficient
c     c10disp  1/r10 coefficient
c     cxd  atom-atom exchange dispersion coefficient 
c     axd  atom-atom exchange dispersion damping parameter
c     cxdla  atom-lp exchange dispersion coefficient 
c     axdla  atom-lp exchange dispersion damping parameter
c     cxdlp  lp-lp exchange dispersion coefficient 
c     axdlp  lp-lp exchange dispersion damping parameter
c     scdp   atom-atom dispersion coefficient
c     facdispij global dispersion coefficient 
c     discoff  dispersion coefficient 
c     colpa  atom-lp dispersion coefficient 
c     colp  lp-lp dispersion coefficient 
c
c     sibfadisp atom-type indexed dispersion vdw radii
c     vdwdisp   dispersion vdw radii for the current structure
c
c      
      module dispersion
      use sizes
      real*8 admp6,admp8,admp10,bdmp
      real*8 c6disp,c8disp,c10disp
      real*8 cxd,axd,cxdlp,axdlp,cxdla,axdla
      real*8 scdp,colpa,colp
      real*8 facdispij,discof
      real*8 sibfadisp(3,maxtyp)
      real*8, allocatable :: vdwdisp1(:),vdwdisp2(:),vdwdisp3(:)
      save
      end 
