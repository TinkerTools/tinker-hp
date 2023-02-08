c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  module dsppot  --  dispersion interaction scale factors  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     dsp2scale   scale factor for 1-2 dispersion energy interactions
c     dsp3scale   scale factor for 1-3 dispersion energy interactions
c     dsp4scale   scale factor for 1-4 dispersion energy interactions
c     dsp5scale   scale factor for 1-5 dispersion energy interactions
c     use_dcorr   flag to use long range dispersion correction
c     dspscal_ik  pair rscale interactions container
c     dspscal_val rscale value of rscak_ik interaction
c     n_dspscal   number of rscaled interactions
c
c
#include "tinker_macro.h"
      module dsppot
      implicit none
      real(t_p) dsp2scale
      real(t_p) dsp3scale
      real(t_p) dsp4scale
      real(t_p) dsp5scale
      logical use_dcorr
      integer   n_dspscal
      integer  ,allocatable:: dspscal_ik(:,:)
      real(t_p),allocatable:: dspscal_val(:)
      end
