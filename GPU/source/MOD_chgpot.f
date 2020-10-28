c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #####################################################################
c     ##                                                                 ##
c     ##  module chgpot  --  specifics of charge-charge functional form  ##
c     ##                                                                 ##
c     #####################################################################
c
c
c     electric   energy factor in kcal/mole for current force field
c     dielec     dielectric constant for electrostatic interactions
c     ebuffer    electrostatic buffering constant added to distance
c     c2scale    factor by which 1-2 charge interactions are scaled
c     c3scale    factor by which 1-3 charge interactions are scaled
c     c4scale    factor by which 1-4 charge interactions are scaled
c     c5scale    factor by which 1-5 charge interactions are scaled
c     neutnbr    logical flag governing use of neutral group neighbors
c     neutcut    logical flag governing use of neutral group cutoffs
c     ccorrect_ik      pair cscale interactions container
c     ccorrect_scale   vscale value of ccorrect_ik interaction
c     n_cscale         number of cscale interactions
c
c
#include "tinker_precision.h"
      module chgpot
      implicit none
      real(t_p) electric
      real(t_p) dielec,ebuffer
      real(t_p) c2scale,c3scale
      real(t_p) c4scale,c5scale
      integer n_cscale
      integer  ,allocatable:: ccorrect_ik(:,:)
      real(t_p),allocatable:: ccorrect_scale(:)
      logical neutnbr,neutcut
      end
