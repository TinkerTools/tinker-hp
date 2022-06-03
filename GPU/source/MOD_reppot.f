c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ##############################################################
c     ##                                                          ##
c     ##  module reppot  --  repulsion interaction scale factors  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     r2scale   scale factor for 1-2 repulsion energy interactions
c     r3scale   scale factor for 1-3 repulsion energy interactions
c     r4scale   scale factor for 1-4 repulsion energy interactions
c     r5scale   scale factor for 1-5 repulsion energy interactions
c
c
#include "tinker_precision.h"
      module reppot
      implicit none
      real(t_p) r2scale
      real(t_p) r3scale
      real(t_p) r4scale
      real(t_p) r5scale
      end
