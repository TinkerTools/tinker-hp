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
c     radrule_i   combining rule for atomic size parameters
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
      enum, bind(C)
      enumerator ARITHMETIC_RL, GEOMETRIC_RL, CUBIC_MEAN_RL, MMFF94_RL
      enumerator   HARMONIC_RL,       HHG_RL,        W_H_RL,DEFAULT_RL
      end enum
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
      logical  use_vcorr, radepsOpt_l
      integer radrule_i, epsrule_i
      character*5 vdwindex
      character*5 radtyp
      character*8 radsiz,gausstyp
      character*10 radrule,epsrule
      character*13 vdwtyp
!$acc declare create (v2scale,v3scale,v4scale,v5scale)

      contains

      subroutine set_vdwradepsrule( arule,irule )
      implicit none
      character*10,intent(in) :: arule
      integer     ,intent(out):: irule

      if (arule(1:6).eq.'MMFF94') then
         irule = MMFF94_RL
      else if (arule(1:10).eq.'ARITHMETIC') then
         irule = ARITHMETIC_RL
      else if (arule(1:9).eq.'GEOMETRIC') then
         irule = GEOMETRIC_RL
      else if (arule(1:9).eq.'HARMONIC') then
         irule = HARMONIC_RL
      else if (arule(1:10).eq.'CUBIC-MEAN') then
         irule = CUBIC_MEAN_RL
      else if (arule(1:3).eq.'HHG') then
         irule = HHG_RL
      else if (arule(1:3).eq.'W-H') then
         irule = W_H_RL
      else
         irule = DEFAULT_RL
      end if
      end subroutine

      end
