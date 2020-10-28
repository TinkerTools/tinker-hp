c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #####################################################################
c     ##                                                                 ##
c     ##  module vdwpot  --  specifics of van der Waals functional form  ##
c     ##                                                                 ##
c     #####################################################################
c
c
c     abuck       value of "A" constant in Buckingham vdw potential
c     bbuck       value of "B" constant in Buckingham vdw potential
c     cbuck       value of "C" constant in Buckingham vdw potential
c     ghal        value of "gamma" in buffered 14-7 vdw potential
c     dhal        value of "delta" in buffered 14-7 vdw potential
c     v2scale     factor by which 1-2 vdw interactions are scaled
c     v3scale     factor by which 1-3 vdw interactions are scaled
c     v4scale     factor by which 1-4 vdw interactions are scaled
c     v5scale     factor by which 1-5 vdw interactions are scaled
c     igauss      coefficients of Gaussian fit to vdw potential
c     ngauss      number of Gaussians used in fit to vdw potential
c     use_vcorr   flag to use long range vdw der Waals correction
c     vdwindex    indexing mode (atom type or class) for vdw parameters
c     vdwtyp      type of van der Waals potential energy function
c     radtyp      type of parameter (sigma or R-min) for atomic size
c     radsiz      atomic size provided as radius or diameter
c     radrule     combining rule for atomic size parameters
c     epsrule     combining rule for vdw well depth parameters
c     gausstyp    type of Gaussian fit to van der Waals potential
c     vcorrect_ik      pair mscale interactions container
c     vcorrect_scale   vscale value of vcorrect_ik interaction
c     n_vscale         number of vscale interactions
c
c
#include "tinker_precision.h"
      module vdwpot
      implicit none
      integer maxgauss
      parameter (maxgauss=10)
      integer ngauss
      real(t_p) abuck,bbuck,cbuck
      real(t_p) ghal,dhal
      real(t_p) v2scale,v3scale
      real(t_p) v4scale,v5scale
      real(t_p) igauss(2,maxgauss)
      integer  n_vscale
      integer  ,allocatable:: vcorrect_ik(:,:)
      integer  ,allocatable:: vcorrect_type(:)
      real(t_p),allocatable:: vcorrect_scale(:)
      logical  use_vcorr
      character*5 vdwindex
      character*5 radtyp
      character*8 radsiz,gausstyp
      character*10 radrule,epsrule
      character*13 vdwtyp
!$acc declare create(dhal,ghal)
!$acc declare create (v2scale,v3scale,v4scale,v5scale)
      end
