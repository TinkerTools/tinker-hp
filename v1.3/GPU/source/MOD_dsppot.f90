!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!
!     ###############################################################
!     ##                                                           ##
!     ##  module dsppot  --  dispersion interaction scale factors  ##
!     ##                                                           ##
!     ###############################################################
!
!
!     dsp2scale   scale factor for 1-2 dispersion energy interactions
!     dsp3scale   scale factor for 1-3 dispersion energy interactions
!     dsp4scale   scale factor for 1-4 dispersion energy interactions
!     dsp5scale   scale factor for 1-5 dispersion energy interactions
!     use_dcorr   flag to use long range dispersion correction
!     dspscal_ik  pair rscale interactions container
!     dspscal_val rscale value of rscak_ik interaction
!     n_dspscal   number of rscaled interactions
!
!
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
