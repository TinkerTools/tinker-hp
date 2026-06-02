!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ##############################################################
!     ##                                                          ##
!     ##  module reppot  --  repulsion interaction scale factors  ##
!     ##                                                          ##
!     ##############################################################
!
!
!     r2scale   scale factor for 1-2 repulsion energy interactions
!     r3scale   scale factor for 1-3 repulsion energy interactions
!     r4scale   scale factor for 1-4 repulsion energy interactions
!     r5scale   scale factor for 1-5 repulsion energy interactions
!     rscal_ik  pair rscale interactions container
!     rscal_val rscale value of rscak_ik interaction
!     n_rscal   number of rscaled interactions
!
!
#include "tinker_macro.h"
module reppot
   implicit none
   real(t_p) r2scale
   real(t_p) r3scale
   real(t_p) r4scale
   real(t_p) r5scale
   integer   n_rscal
   integer  ,allocatable:: rscal_ik(:,:)
   real(t_p),allocatable:: rscal_val(:)
end
