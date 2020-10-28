#ifndef SWITCH_RESPA
#define SWITCH_RESPA
#include "tinker_cudart.h"
      M_subroutine 
     &             switch_respa_inl(r,rc,lambda,s,ds)
!$acc routine
      use tinheader ,only: ti_p
      implicit none
      real(t_p),intent(in)   :: r,rc,lambda
      real(t_p),intent(inout):: s,ds
      real(t_p) u,du
c
      if (r.le.(rc-lambda)) then
        s = 1.0_ti_p
        ds = 0.0_ti_p
      else if (((rc-lambda).le.r).and.(r.le.rc)) then
        u = (1.0_ti_p/lambda)*(r-rc+lambda)
        du = 1.0_ti_p/lambda
        s = 1 + u**3*(15*u - 6*u**2 - 10)
        ds = 3*du*u**2*(15*u - 6*u**2 - 10) + 
     $   u**3*(15*du - 12*du*u)
      else if (r.ge.rc) then
        s = 0.0_ti_p
        ds = 0.0_ti_p
      end if
      end
#endif
