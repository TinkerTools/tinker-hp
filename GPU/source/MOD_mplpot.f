c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ##################################################################
c     ##                                                              ##
c     ##  module mplpot  --  specifics of atomic multipole functions  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     m2scale   scale factor for 1-2 multipole energy interactions
c     m3scale   scale factor for 1-3 multipole energy interactions
c     m4scale   scale factor for 1-4 multipole energy interactions
c     m5scale   scale factor for 1-5 multipole energy interactions
c     mcorrect_ik      pair mscale interactions container
c     mcorrect_scale   mscale value of mcorrect_ik interaction
c     n_mscale         number of mscale interactions
c
c
#include "tinker_precision.h"
      module mplpot
      implicit none
      real(t_p) m2scale,m3scale
      real(t_p) m4scale,m5scale
      integer n_mscale
      integer  ,allocatable:: mcorrect_ik(:,:)
      real(t_p),allocatable:: mcorrect_scale(:)
!$acc declare create(m2scale,m3scale,m4scale,m5scale)
      end
