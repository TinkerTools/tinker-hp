c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ################################################################################
c     ##                                                                            ##
c     ##  module repulsion  --   repulsion caracteristics for the current structure ##
c     ##                                                                            ##
c     ################################################################################
c
c     alpha  repulsion square overlap/r damping parameter
c     alpha2  repulsion square overlap/r2 damping parameter
c     cvrep11 conversion factor for bond-bond square overlap/r term
c     cvrep12 conversion factor for bond-bond square overlap/r2 term
c     cvrep21 conversion factor for bond-lp square overlap/r term
c     cvrep22 conversion factor for bond-lp square overlap/r2 term
c     cvrep31 conversion factor for lp-lp square overlap/r term
c     cvrep32 conversion factor for lp-lp square overlap/r2 term
c     
c     sibfarep atom type indexed repulsion-vdw radii
c     forb  atom type indexed prefactor to compute overlap
c     gorb  atom type indexed prefactor to compute overlap
c     forbrep  prefactor to compute overlap for the current structure
c     gorbrep  prefactor to compute overlap for the current structure
c
c     vdwrep repulsion-vdw radii for the current structure
c     chyb     hybridization coefficients for the lone pairs of the current structure
c     
c
c      
      module repulsion
      use sizes
      implicit none
      real*8 alpha,alpha2
      real*8 cvrep11,cvrep12
      real*8 cvrep21,cvrep22
      real*8 cvrep31,cvrep32
      real*8 sibfarep(maxtyp),forb(100,100),gorb(maxtyp)
      real*8, allocatable :: vdwrep(:),forbrep(:,:),gorbrep(:)
      real*8, allocatable :: dincr_lprep(:)
      real*8, allocatable :: chyb(:,:)
      save
      end
