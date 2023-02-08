#ifndef DAMPING_INC_F
#define DAMPING_INC_F
#include "tinker_macro.h"

c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine dampewald  --  find Ewald damping coefficients  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "dampewald" finds coefficients for error function damping used
c     for Ewald real space interactions
c
c
      M_subroutine dampewald_inl 
     &            (rorder,r,r2,ewald,scale,dmpe)
      use math ,only: sqrtpi
#ifdef TINKER_CUF
      use utilcu ,only: f_erfc
#  if  TINKER_SINGLE_PREC + TINKER_MIXED_PREC
     &           , f_exp
#  endif
#endif
!$acc routine
      implicit none
      integer  ,intent(in):: rorder
      real(t_p),intent(in):: r,r2,scale,ewald
      real(t_p),intent(out):: dmpe(*)
      integer i,niter
      real(t_p) bfac,aesq2,afac,exp2a,ra,ir2,bni,bni_1
c
c     Match with reference
c       1  -->  1
c       2  -->  3
c       3  -->  5
c       4  -->  7
c       5  -->  9
c     compute the successive Ewald damping factors
c
      ra      = ewald* r
      exp2a   = f_exp(-ra*ra)
      bni_1   = f_erfc(ra) / r
      dmpe(1) = scale* bni_1
      aesq2   = 2.0* ewald* ewald
      afac    = merge((sqrtpi*ewald)**(-1), 0.0, ewald.gt.0.0)
      ir2     = r2**(-1)
      niter   = ishft(rorder-1,-1)
      do i = 1, niter
         bfac      = real(2*i-1,t_p)
         afac      = aesq2 * afac
         bni       = (bfac*bni_1+afac*exp2a) * ir2
         dmpe(i+1) = scale * bni
         bni_1     = bni
      end do
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine dampthole  --  find Thole damping coefficients  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "dampthole" generates coefficients for the Thole damping
c     function for powers of the interatomic distance
c
c     literature reference:
c
c     B. T. Thole, "Molecular Polarizabilities Calculated with a
c     Modified Dipole Interaction", Chemical Physics, 59, 341-350 (1981)
c
c
      M_subroutine dampthole_inl
     &            (rorder,r,damp,pgamma,use_dirdamp,dmpik)
#ifdef TINKER_CUF
      use utilcu ,only: f_erfc
#  if  TINKER_SINGLE_PREC + TINKER_MIXED_PREC
     &           , f_exp
#  endif
#endif
!$acc routine
      implicit none
      integer  ,intent(in):: rorder
      logical  ,intent(in):: use_dirdamp
      real(t_p),intent(in):: r,pgamma
      real(t_p),intent(inout):: damp
      real(t_p),intent(out):: dmpik(*)
      real(t_p) damp2,damp3,expdamp
c
c     use alternate Thole model for AMOEBA+ direct polarization
c
c     Match with Refi
c          1 --> 3   i&k
c          2 --> 5   i&k
c          3 --> 7   i&k
c          4 --> 9   i&k
      dmpik(1) = 1.0
      dmpik(2) = 1.0
      dmpik(3) = 1.0
      if (rorder.gt.9) dmpik(4) = 1.0
      if (use_dirdamp) then
         if (damp.ne.0.0 .and. pgamma.ne.0.0) then
            damp = pgamma * (r/damp)**(1.5)
            if (damp .lt. 50.0) then
               expdamp = exp(-damp)
               dmpik(1) = 1.0 - expdamp
               dmpik(2) = 1.0 - expdamp*(1.0+0.5*damp)
               if (rorder .ge. 7) then
                  damp2 = damp * damp
                  dmpik(3) = 1.0 - expdamp*(1.0+0.65*damp
     &                                  +0.15*damp2)
               end if
            end if
         end if
c
c     use original AMOEBA Thole polarization damping factors
c
      else
         if (damp.ne.0.0 .and. pgamma.ne.0.0) then
            damp = pgamma * (r/damp)**3
            if (damp .lt. 50.0) then
               expdamp = f_exp(-damp)
               dmpik(1) = 1.0 - expdamp
               dmpik(2) = 1.0 - expdamp*(1.0+damp)
               if (rorder .ge. 7) then
                  damp2 = damp * damp
                  dmpik(3) = 1.0 - expdamp*(1.0+damp+0.6*damp2)
                  if (rorder .ge. 9) then
                     damp3 = damp * damp2
                     dmpik(4) = 1.0 - expdamp*(1.0+damp
     &                                     +(18.0/35.0)*damp2
     &                                     +(9.0/35.0)*damp3)
                  end if
               end if
            end if
         end if
      end if
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine dampthole2  --  original Thole damping values  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "dampthole2" finds coefficients for the original Thole damping
c     function used by AMOEBA and for mutual polarization by AMOEBA+
c
c     literature reference:
c
c     B. T. Thole, "Molecular Polarizabilities Calculated with a
c     Modified Dipole Interaction", Chemical Physics, 59, 341-350 (1981)
c
c
      M_subroutine dampthole2_inl 
     &            (damp,pgamma,rorder,r,dmpik)
#ifdef TINKER_CUF
      use utilcu  ,only: f_abs
#  if TINKER_SINGLE_PREC + TINKER_MIXED_PREC
     &            , f_sqrt,f_exp
#  endif
#endif
!$acc routine
      implicit none
      integer  ,intent(in) :: rorder
      real(t_p),intent(in) :: r,damp,pgamma
      real(t_p),intent(out):: dmpik(*)
      real(t_p) sdamp, damp2, damp3, expdamp
c
c     initialize the Thole damping factors to a value of one
c
c     Match with Refi
c          1 --> 3   i&k
c          2 --> 5   i&k
c          3 --> 7   i&k
c          4 --> 9   i&k
      dmpik(1) = 1.0
      dmpik(2) = 1.0
      if (rorder.ge.7) dmpik(3) = 1.0
      if (rorder.ge.9) dmpik(4) = 1.0
c
c     assign original Thole polarization model damping factors
c
      if (damp.ne.0.0 .and. pgamma.ne.0.0) then
         sdamp = pgamma * (r/damp)**3
         if (sdamp .lt. 50.0) then
            expdamp = f_exp(-sdamp)
            dmpik(1) = 1.0 - expdamp
            dmpik(2) = 1.0 - expdamp*(1.0+sdamp)
            if (rorder.ge.7) then
               damp2 = sdamp * sdamp
               dmpik(3) = 1.0 - expdamp*(1.0+sdamp+0.6*damp2)
               if (rorder.ge.9) then
                  damp3 = sdamp * damp2
                  dmpik(4) = 1.0 - expdamp*(
     &                     1.0+sdamp+(18.0/35.0)*damp2+(9.0/35.0)*damp3)
               end if
            end if
         end if
      end if
      end
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine damppole  --  penetration damping coefficents  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "damppole" generates coefficients for the charge penetration
c     damping function for powers of the interatomic distance
c
c     literature references:
c
c     L. V. Slipchenko and M. S. Gordon, "Electrostatic Energy in the
c     Effective Fragment Potential Method: Theory and Application to
c     the Benzene Dimer", Journal of Computational Chemistry, 28,
c     276-291 (2007)  [Gordon f1 and f2 models]
c
c     J. A. Rackers, Q. Wang, C. Liu, J.-P. Piquemal, P. Ren and
c     J. W. Ponder, "An Optimized Charge Penetration Model for Use with
c     the AMOEBA Force Field", Physical Chemistry Chemical Physics, 19,
c     276-291 (2017)
c
c
      M_subroutine damppole_inl 
     &            (r,rorder,pentyp,alphai,alphak,dmpi,dmpk,dmpik)
      use mplpot  ,only: PT_GORDON1,PT_GORDON2
#ifdef TINKER_CUF
      use utilcu ,only: f_erfc, f_abs
#  if  TINKER_SINGLE_PREC + TINKER_MIXED_PREC
     &           , f_exp
#  endif
#endif
      implicit none
      integer  ,intent(in ):: rorder,pentyp
      real(t_p),intent(in ):: r,alphai,alphak
      real(t_p),intent(out):: dmpi(*),dmpk(*),dmpik(*)
      real(t_p) termi,termk
      real(t_p) termi2,termk2
      real(t_p) alphai2,alphak2
      real(t_p) eps,diff
      real(t_p) expi,expk
      real(t_p) dampi,dampk
      real(t_p) dampi2,dampi3,dampi4,dampi5,dampi6,dampi7,dampi8
      real(t_p) dampk2,dampk3,dampk4,dampk5,dampk6
      parameter( eps=0.001 )
c
c     Match with Refi
c          1 --> 1   i&k
c          2 --> 3   i&k
c          3 --> 5   i&k
c          4 --> 7   i&k
c          5 --> 9   i&k
c          6 --> 11  ik
c
c     compute tolerance and exponential damping factors
c
      diff  = f_abs(alphai-alphak)
      dampi = alphai * r
      dampk = alphak * r
      expi  = f_exp(-dampi)
      expk  = f_exp(-dampk)
c
c     core-valence charge penetration damping for Gordon f1
c
      if (pentyp .eq. PT_GORDON1) then
         dampi2 = dampi * dampi
         dampi3 = dampi * dampi2
         dampi4 = dampi2 * dampi2
         dampi5 = dampi2 * dampi3
         dmpi(1) = 1.0 - (1.0 + 0.5*dampi)*expi
         dmpi(2) = 1.0 - (1.0 + dampi + 0.5*dampi2)*expi
         dmpi(3) = 1.0 - (1.0 + dampi + 0.5*dampi2
     &                + dampi3/6.0)*expi
         dmpi(4) = 1.0 - (1.0 + dampi + 0.5*dampi2
     &                + dampi3/6.0 + dampi4/30.0)*expi
         dmpi(5) = 1.0 - (1.0 + dampi + 0.5*dampi2
     &                + dampi3/6.0 + 4.0*dampi4/105.0
     &                + dampi5/210.0)*expi
         if (diff .lt. eps) then
            dmpk(1) = dmpi(1)
            dmpk(2) = dmpi(2)
            dmpk(3) = dmpi(3)
            dmpk(4) = dmpi(4)
            dmpk(5) = dmpi(5)
         else
            dampk2 = dampk * dampk
            dampk3 = dampk * dampk2
            dampk4 = dampk2 * dampk2
            dampk5 = dampk2 * dampk3
            dmpk(1) = 1.0 - (1.0 + 0.5*dampk)*expk
            dmpk(2) = 1.0 - (1.0 + dampk + 0.5*dampk2)*expk
            dmpk(3) = 1.0 - (1.0 + dampk + 0.5*dampk2
     &                   + dampk3/6.0)*expk
            dmpk(4) = 1.0 - (1.0 + dampk + 0.5*dampk2
     &                   + dampk3/6.0 + dampk4/30.0)*expk
            dmpk(5) = 1.0 - (1.0 + dampk + 0.5*dampk2
     &                   + dampk3/6.0 + 4.0*dampk4/105.0
     &                   + dampk5/210.0)*expk
         end if
c
c     valence-valence charge penetration damping for Gordon f1
c
         if (diff .lt. eps) then
            dampi6   = dampi3 * dampi3
            dampi7   = dampi3 * dampi4
            dmpik(1) = 1.0 - (1.0 + 11.0*dampi/16.0
     &                    + 3.0*dampi2/16.0
     &                    + dampi3/48.0)*expi
            dmpik(2) = 1.0 - (1.0 + dampi + 0.5*dampi2
     &                    + 7.0*dampi3/48.0
     &                    + dampi4/48.0)*expi
            dmpik(3) = 1.0 - (1.0 + dampi + 0.5*dampi2
     &                    + dampi3/6.0 + dampi4/24.0
     &                    + dampi5/144.0)*expi
            dmpik(4) = 1.0 - (1.0 + dampi + 0.5*dampi2
     &                    + dampi3/6.0 + dampi4/24.0
     &                    + dampi5/120.0 + dampi6/720.0)*expi
            dmpik(5) = 1.0 - (1.0 + dampi + 0.5*dampi2
     &                    + dampi3/6.0 + dampi4/24.0
     &                    + dampi5/120.0 + dampi6/720.0
     &                    + dampi7/5040.0)*expi
            if (rorder .ge. 11) then
               dampi8 = dampi4 * dampi4
               dmpik(6) = 1.0 - (1.0 + dampi + 0.5*dampi2
     &                        + dampi3/6.0 + dampi4/24.0
     &                        + dampi5/120.0 + dampi6/720.0
     &                        + dampi7/5040.0 + dampi8/45360.0)*expi
            end if
         else
            alphai2 = alphai * alphai
            alphak2 = alphak * alphak
            termi = alphak2 / (alphak2-alphai2)
            termk = alphai2 / (alphai2-alphak2)
            termi2 = termi * termi
            termk2 = termk * termk
            dmpik(1) = 1.0 - termi2*(1.0 + 2.0*termk
     &                    + 0.5*dampi)*expi
     &                 - termk2*(1.0 + 2.0*termi
     &                      + 0.5*dampk)*expk
            dmpik(2) = 1.0 - termi2*(1.0+dampi+0.5*dampi2)*expi
     &                    - termk2*(1.0+dampk+0.5*dampk2)*expk
     &                    - 2.0*termi2*termk*(1.0+dampi)*expi
     &                    - 2.0*termk2*termi*(1.0+dampk)*expk
            dmpik(3) = 1.0 - termi2*(1.0 + dampi + 0.5*dampi2
     &                    + dampi3/6.0)*expi
     &                 - termk2*(1.0 + dampk + 0.5*dampk2
     &                      + dampk3/6.0)*expk
     &                 - 2.0*termi2*termk
     &                      *(1.0 + dampi + dampi2/3.0)*expi
     &                 - 2.0*termk2*termi
     &                      *(1.0 + dampk + dampk2/3.0)*expk
            dmpik(4) = 1.0 - termi2*(1.0 + dampi + 0.5*dampi2
     &                    + dampi3/6.0 + dampi4/30.0)*expi
     &                 - termk2*(1.0 + dampk + 0.5*dampk2
     &                      + dampk3/6.0 + dampk4/30.0)*expk
     &                 - 2.0*termi2*termk*(1.0 + dampi
     &                      + 2.0*dampi2/5.0 + dampi3/15.0)*expi
     &                 - 2.0*termk2*termi*(1.0 + dampk
     &                      + 2.0*dampk2/5.0 + dampk3/15.0)*expk
            dmpik(5) = 1.0 - termi2*(1.0 + dampi + 0.5*dampi2
     &                    + dampi3/6.0 + 4.0*dampi4/105.0
     &                    + dampi5/210.0)*expi
     &                 - termk2*(1.0 + dampk + 0.5*dampk2
     &                      + dampk3/6.0 + 4.0*dampk4/105.0
     &                      + dampk5/210.0)*expk
     &                 - 2.0*termi2*termk*(1.0 + dampi
     &                      + 3.0*dampi2/7.0
     &                      + 2.0*dampi3/21.0
     &                      + dampi4/105.0)*expi 
     &                 - 2.0*termk2*termi*(1.0 + dampk
     &                      + 3.0*dampk2/7.0
     &                      + 2.0*dampk3/21.0
     &                      + dampk4/105.0)*expk
            if (rorder .ge. 11) then
               dampi6 = dampi3 * dampi3
               dampk6 = dampk3 * dampk3
               dmpik(6) = 1.0 - termi2*(1.0 + dampi
     &                        + 0.5*dampi2 + dampi3/6.0
     &                        + 5.0*dampi4/126.0
     &                        + 2.0*dampi5/315.0
     &                        + dampi6/1890.0)*expi
     &                     - termk2*(1.0 + dampk
     &                          + 0.5*dampk2 + dampk3/6.0
     &                          + 5.0*dampk4/126.0
     &                          + 2.0*dampk5/315.0
     &                          + dampk6/1890.0)*expk
     &                     - 2.0*termi2*termk*(1.0 + dampi
     &                          + 4.0*dampi2/9.0 + dampi3/9.0
     &                          + dampi4/63.0 + dampi5/945.0)*expi
     &                     - 2.0*termk2*termi*(1.0 + dampk
     &                          + 4.0*dampk2/9.0 + dampk3/9.0
     &                          + dampk4/63.0 + dampk5/945.0)*expk 
            end if
         end if
c
c     core-valence charge penetration damping for Gordon f2
c
      else if (pentyp .eq. PT_GORDON2) then
         dampi2 = dampi * dampi
         dampi3 = dampi * dampi2
         dmpi(1) = 1.0 - expi
         dmpi(2) = 1.0 - (1.0 + dampi)*expi
         dmpi(3) = 1.0 - (1.0 + dampi + dampi2/3.0)*expi
         dmpi(4) = 1.0 - (1.0 + dampi + 0.4*dampi2
     &                + dampi3/15.0)*expi
         if (diff .lt. eps) then
            dmpk(1) = dmpi(1)
            dmpk(2) = dmpi(2)
            dmpk(3) = dmpi(3)
            dmpk(4) = dmpi(4)
         else
            dampk2 = dampk * dampk
            dampk3 = dampk * dampk2
            dmpk(1) = 1.0 - expk
            dmpk(2) = 1.0 - (1.0 + dampk)*expk
            dmpk(3) = 1.0 - (1.0 + dampk + dampk2/3.0)*expk
            dmpk(4) = 1.0 - (1.0 + dampk + 0.4*dampk2
     &                   + dampk3/15.0)*expk
         end if
c
c     valence-valence charge penetration damping for Gordon f2
c
         dampi4 = dampi2 * dampi2
         dampi5 = dampi2 * dampi3
         if (diff .lt. eps) then
            dampi6 = dampi3 * dampi3
            dmpik(1) = 1.0 - (1.0 + 0.5*dampi)*expi
            dmpik(2) = 1.0 - (1.0 + dampi + 0.5*dampi2)*expi
            dmpik(3) = 1.0 - (1.0 + dampi + 0.5*dampi2
     &                    + dampi3/6.0)*expi
            dmpik(4) = 1.0 - (1.0 + dampi + 0.5*dampi2
     &                    + dampi3/6.0 + dampi4/30.0)*expi
            dmpik(5) = 1.0 - (1.0 + dampi + 0.5*dampi2
     &                    + dampi3/6.0 + 4.0*dampi4/105.0
     &                    + dampi5/210.0)*expi
            if (rorder .ge. 11) then
               dmpik(6) = 1.0 - (1.0 + dampi + 0.5*dampi2
     &                        + dampi3/6.0 + 5.0*dampi4/126.0
     &                        + 2.0*dampi5/315.0
     &                        + dampi6/1890.0)*expi
            end if
         else
            dampk4 = dampk2 * dampk2
            dampk5 = dampk2 * dampk3
            alphai2 = alphai * alphai
            alphak2 = alphak * alphak
            termi = alphak2 / (alphak2-alphai2)
            termk = alphai2 / (alphai2-alphak2)
            dmpik(1) = 1.0 - termi*expi - termk*expk
            dmpik(2) = 1.0 - termi*(1.0 + dampi)*expi
     &                    - termk*(1.0 + dampk)*expk
            dmpik(3) = 1.0 - termi*(1.0 + dampi + dampi2/3.0)*expi
     &                    - termk*(1.0 + dampk + dampk2/3.0)*expk
            dmpik(4) = 1.0 - termi*(1.0 + dampi + 0.4*dampi2
     &                    + dampi3/15.0)*expi
     &                 - termk*(1.0 + dampk + 0.4*dampk2
     &                      + dampk3/15.0)*expk
            dmpik(5) = 1.0 - termi*(1.0 + dampi + 3.0*dampi2/7.0
     &                    + 2.0*dampi3/21.0 + dampi4/105.0)*expi
     &                 - termk*(1.0 + dampk + 3.0*dampk2/7.0
     &                      + 2.0*dampk3/21.0 + dampk4/105.0)*expk
            if (rorder .ge. 11) then
               dmpik(6) = 1.0 - termi*(1.0 + dampi
     &                        + 4.0*dampi2/9.0 + dampi3/9.0
     &                        + dampi4/63.0 + dampi5/945.0)*expi
     &                     - termk*(1.0 + dampk + 4.0*dampk2/9.0
     &                        + dampk3/9.0
     &                        + dampk4/63.0 + dampk5/945.0)*expk
            end if
         end if
      end if
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine dampdir  --  direct field damping coefficents  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "dampdir" generates coefficients for the direct field damping
c     function for powers of the interatomic distance
c
c
      M_subroutine dampdir_inl
     &            (pentyp,r,alphai,alphak,dmpi,dmpk)
      use mplpot ,only: PT_GORDON1,PT_GORDON2
#ifdef TINKER_CUF
      use utilcu  ,only: f_abs
#  if TINKER_SINGLE_PREC + TINKER_MIXED_PREC
     &            , f_sqrt,f_exp
#  endif
#endif
!$acc routine
      implicit none
      integer  ,intent(in):: pentyp
      real(t_p),intent(in):: r,alphai,alphak
      real(t_p),intent(out):: dmpi(*),dmpk(*)
      real(t_p) eps,diff
      real(t_p) expi,expk,dampi,dampk,dampi2,dampk2
     &         ,dampi3,dampk3,dampi4,dampk4
      parameter( eps=0.001 )
c
c     Match with reference
c       1 --> 3
c       2 --> 5
c       3 --> 7
c     compute tolerance and exponential damping factors
c
      diff  = f_abs(alphai-alphak)
      dampi = alphai * r
      dampk = alphak * r
      expi  = f_exp(-dampi)
      expk  = f_exp(-dampk)
c
c     core-valence charge penetration damping for Gordon f1
c
      if (pentyp .eq. PT_GORDON1) then
         dampi2  = dampi * dampi
         dampi3  = dampi * dampi2
         dampi4  = dampi2 * dampi2
         dmpi(1) = 1.0 - (1.0 + dampi + 0.5*dampi2)*expi
         dmpi(2) = 1.0 - (1.0 + dampi + 0.5*dampi2 
     &                + dampi3/6.0)*expi
         dmpi(3) = 1.0 - (1.0 + dampi + 0.5*dampi2
     &                + dampi3/6.0 + dampi4/30.0)*expi
         if (diff .lt. eps) then
            dmpk(1) = dmpi(1)
            dmpk(2) = dmpi(2)
            dmpk(3) = dmpi(3)
         else
            dampk2  = dampk * dampk
            dampk3  = dampk * dampk2
            dampk4  = dampk2 * dampk2
            dmpk(1) = 1.0 - (1.0 + dampk + 0.5*dampk2)*expk
            dmpk(2) = 1.0 - (1.0 + dampk + 0.5*dampk2
     &                   + dampk3/6.0)*expk
            dmpk(3) = 1.0 - (1.0 + dampk + 0.5*dampk2
     &                   + dampk3/6.0 + dampk4/30.0)*expk
         end if
c
c     core-valence charge penetration damping for Gordon f2
c
      else if (pentyp .eq. PT_GORDON2) then
         dampi2  = dampi * dampi
         dampi3  = dampi * dampi2
         dmpi(1) = 1.0 - (1.0 + dampi)*expi
         dmpi(2) = 1.0 - (1.0 + dampi + dampi2/3.0)*expi
         dmpi(3) = 1.0 - (1.0 + dampi + 0.4*dampi2
     &                + dampi3/15.0)*expi
         if (diff .lt. eps) then
            dmpk(1) = dmpi(1)
            dmpk(2) = dmpi(2)
            dmpk(3) = dmpi(3)
         else
            dampk2 = dampk * dampk
            dampk3 = dampk * dampk2
            dmpk(1) = 1.0 - (1.0 + dampk)*expk
            dmpk(2) = 1.0 - (1.0 + dampk + dampk2/3.0)*expk
            dmpk(3) = 1.0 - (1.0 + dampk + 0.4*dampk2
     &                   + dampk3/15.0)*expk
         end if
      else
         dmpi(1) = 1.0
      end if
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine dampmut  --  mutual field damping coefficents  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "dampmut" generates coefficients for the mutual field damping
c     function for powers of the interatomic distance
c
c
      M_subroutine dampmut_inl
     &            (r,alphai,alphak,pentyp,dmpik)
      use mplpot  ,only: PT_GORDON1,PT_GORDON2
#ifdef TINKER_CUF
      use utilcu  ,only: f_abs
#  if TINKER_SINGLE_PREC + TINKER_MIXED_PREC
     &            , f_sqrt,f_exp
#  endif
#endif
!$acc routine
      implicit none
      integer  ,intent(in):: pentyp
      real(t_p),intent(in):: r,alphai,alphak
      real(t_p),intent(out):: dmpik(*)
      real(t_p) termi,termk,termi2,termk2
      real(t_p) alphai2,alphak2,diff
      real(t_p) eps,expi,expk
      real(t_p) dampi,dampk,dampi2,dampi3
      real(t_p) dampi4,dampi5,dampk2,dampk3
      parameter(eps=0.001)
c
c     compute tolerance and exponential damping factors
c
      diff = f_abs(alphai-alphak)
      dampi = alphai * r
      dampk = alphak * r
      expi = f_exp(-dampi)
      expk = f_exp(-dampk)
c
c     valence-valence charge penetration damping for Gordon f1
c
      if (pentyp .eq. PT_GORDON1) then
         dampi2 = dampi * dampi
         dampi3 = dampi * dampi2
         if (diff .lt. eps) then
            dampi4 = dampi2 * dampi2
            dampi5 = dampi2 * dampi3
            dmpik(1) = 1.0 - (1.0 + dampi + 0.5*dampi2
     &                    + 7.0*dampi3/48.0
     &                    + dampi4/48.0)*expi
            dmpik(2) = 1.0 - (1.0 + dampi + 0.5*dampi2
     &                    + dampi3/6.0 + dampi4/24.0
     &                    + dampi5/144.0)*expi
         else
            dampk2 = dampk * dampk
            dampk3 = dampk * dampk2
            alphai2 = alphai * alphai
            alphak2 = alphak * alphak
            termi = alphak2 / (alphak2-alphai2)
            termk = alphai2 / (alphai2-alphak2)
            termi2 = termi * termi
            termk2 = termk * termk
            dmpik(1) = 1.0 - termi2*(1.0+dampi+0.5*dampi2)*expi
     &                     - termk2*(1.0+dampk+0.5*dampk2)*expk
     &                     - 2.0*termi2*termk*(1.0+dampi)*expi
     &                     - 2.0*termk2*termi*(1.0+dampk)*expk
            dmpik(2) = 1.0 - termi2*(1.0+dampi+0.5*dampi2
     &                            +dampi3/6.0)*expi
     &                    - termk2*(1.0+dampk+0.5*dampk2
     &                         +dampk3/6.00)*expk
     &                    - 2.0*termi2*termk
     &                         *(1.0+dampi+dampi2/3.0)*expi
     &                    - 2.0*termk2*termi*(1.0+dampk+dampk2/3.0)*expk
         end if
c
c     valence-valence charge penetration damping for Gordon f2
c
      else if (pentyp .eq. PT_GORDON2) then
         dampi2 = dampi * dampi
         if (diff .lt. eps) then
            dampi3   = dampi * dampi2
            dmpik(1) = 1.0 - (1.0 + dampi + 0.5*dampi2)*expi
            dmpik(2) = 1.0 - (1.0 + dampi + 0.5*dampi2
     &                    + dampi3/6.0)*expi
         else
            dampk2   = dampk * dampk
            alphai2  = alphai * alphai
            alphak2  = alphak * alphak
            termi    = alphak2 / (alphak2-alphai2)
            termk    = alphai2 / (alphai2-alphak2)
            dmpik(1) = 1.0 - termi*(1.0 + dampi)*expi
     &                     - termk*(1.0 + dampk)*expk
            dmpik(2) = 1.0 - termi*(1.0 + dampi + dampi2/3.0)*expi
     &                     - termk*(1.0 + dampk + dampk2/3.0)*expk
         end if
      end if
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine damprep  --  find repulsion damping coefficents  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "damprep" generates coefficients for the Pauli repulsion
c     damping function for powers of the interatomic distance
c
c     literature reference:
c
c     J. A. Rackers and J. W. Ponder, "Classical Pauli Repulsion: An
c     Anisotropic, Atomic Multipole Model", Journal of Chemical Physics,
c     150, 084104 (2019)
c
c
      M_subroutine damprep_inl
     &            (r,r2,rr1,rr3,rr5,rr7,rr9,rr11
     &            ,rorder,dmpi,dmpk,dmpik)
#ifdef TINKER_CUF
      use utilcu  ,only: f_abs
#  if TINKER_SINGLE_PREC + TINKER_MIXED_PREC
     &            , f_sqrt,f_exp
#  endif
#endif
!$acc routine
      implicit none
      integer  ,intent(in):: rorder
      real(t_p),intent(in):: r,r2,rr1,rr3,rr5,rr7,rr9,rr11,dmpi,dmpk
      real(t_p),intent(out):: dmpik(*)

      real(t_p) r3,r4,r5,r6,r7,r8,d3s,d4s,d5s,s,ds,d2s
      real(t_p) dmpi2 ,dmpk2 ,dmpi22,dmpi23,dmpi24,dmpi25
     &         ,dmpi26,dmpi27,dmpk22,dmpk23,dmpk24,dmpk25,dmpk26
      real(t_p) eps,diff
      real(t_p) expi,expk,dampi,dampk,pre,term,tmp
      parameter( eps=0.001 )
c
c     compute tolerance value for damping exponents
c
      diff= f_abs(dmpi-dmpk)
c
c     treat the case where alpha damping exponents are equal
c
      if (diff .lt. eps) then
         r3 = r2 * r
         r4 = r3 * r
         r5 = r4 * r
         r6 = r5 * r
         dmpi2 = 0.5 * dmpi
         dampi = dmpi2 * r
         expi = f_exp(-dampi)
         dmpi22 = dmpi2 * dmpi2
         dmpi23 = dmpi22 * dmpi2
         dmpi24 = dmpi23 * dmpi2
         dmpi25 = dmpi24 * dmpi2
         pre = 2.0
         s = (r + dmpi2*r2 + dmpi22*r3/3.0) * expi
         ds = (dmpi22*r3 + dmpi23*r4) * expi / 3.0
         d2s = dmpi24 * expi * r5 / 9.0
         d3s = dmpi25 * expi * r6 / 45.0
         if (rorder .ge. 9) then
            r7 = r6 * r
            dmpi26 = dmpi25 * dmpi2
            d4s = (dmpi25*r6 + dmpi26*r7) * expi / 315.0
            if (rorder .ge. 11) then
               r8 = r7 * r
               dmpi27 = dmpi2 * dmpi26
               d5s = (dmpi25*r6 + dmpi26*r7 + dmpi27*r8/3.0)
     &                   * expi / 945.0
            end if
         end if
c
c     treat the case where alpha damping exponents are unequal
c
      else
         r3 = r2 * r
         r4 = r3 * r
         dmpi2 = 0.5 * dmpi
         dmpk2 = 0.5 * dmpk
         dampi = dmpi2 * r
         dampk = dmpk2 * r
         expi = f_exp(-dampi)
         expk = f_exp(-dampk)
         dmpi22 = dmpi2 * dmpi2
         dmpi23 = dmpi22 * dmpi2
         dmpi24 = dmpi23 * dmpi2
         dmpk22 = dmpk2 * dmpk2
         dmpk23 = dmpk22 * dmpk2
         dmpk24 = dmpk23 * dmpk2
         term = dmpi22 - dmpk22
         pre = 128.0 * dmpi23 * dmpk23 / term**4
         tmp = 4.0 * dmpi2 * dmpk2 / term
         s = (dampi-tmp)*expk + (dampk+tmp)*expi
         ds = (dmpi2*dmpk2*r2 - 4.0*dmpi2*dmpk22*r/term
     &            - 4.0*dmpi2*dmpk2/term) * expk
     &      + (dmpi2*dmpk2*r2 + 4.0*dmpi22*dmpk2*r/term
     &            + 4.0*dmpi2*dmpk2/term) * expi
         d2s = (dmpi2*dmpk2*r2/3.0
     &             + dmpi2*dmpk22*r3/3.0
     &             - (4.0/3.0)*dmpi2*dmpk23*r2/term
     &             - 4.0*dmpi2*dmpk22*r/term
     &             - 4.0*dmpi2*dmpk2/term) * expk
     &       + (dmpi2*dmpk2*r2/3.0
     &             + dmpi22*dmpk2*r3/3.0
     &             + (4.0/3.0)*dmpi23*dmpk2*r2/term
     &             + 4.0*dmpi22*dmpk2*r/term
     &             + 4.0*dmpi2*dmpk2/term) * expi
         d3s = (dmpi2*dmpk23*r4/15.0
     &             + dmpi2*dmpk22*r3/5.0
     &             + dmpi2*dmpk2*r2/5.0
     &             - (4.0/15.0)*dmpi2*dmpk24*r3/term
     &             - (8.0/5.0)*dmpi2*dmpk23*r2/term
     &             - 4.0*dmpi2*dmpk22*r/term
     &             - 4.0/term*dmpi2*dmpk2) * expk
     &       + (dmpi23*dmpk2*r4/15.0 
     &             + dmpi22*dmpk2*r3/5.0
     &             + dmpi2*dmpk2*r2/5.0 
     &             + (4.0/15.0)*dmpi24*dmpk2*r3/term
     &             + (8.0/5.0)*dmpi23*dmpk2*r2/term
     &             + 4.0*dmpi22*dmpk2*r/term
     &             + 4.0/term*dmpi2*dmpk2) * expi
         if (rorder .ge. 9) then
            r5 = r4 * r
            dmpi25 = dmpi24 * dmpi2
            dmpk25 = dmpk24 * dmpk2
            d4s = (dmpi2*dmpk24*r5/105.0
     &                + (2.0/35.0)*dmpi2*dmpk23*r4
     &                + dmpi2*dmpk22*r3/7.0
     &                + dmpi2*dmpk2*r2/7.0
     &                - (4.0/105.0)*dmpi2*dmpk25*r4/term
     &                - (8.0/21.0)*dmpi2*dmpk24*r3/term
     &                - (12.0/7.0)*dmpi2*dmpk23*r2/term
     &                - 4.0*dmpi2*dmpk22*r/term
     &                - 4.0*dmpi2*dmpk2/term) * expk
     &          + (dmpi24*dmpk2*r5/105.0
     &                + (2.0/35.0)*dmpi23*dmpk2*r4
     &                + dmpi22*dmpk2*r3/7.0
     &                + dmpi2*dmpk2*r2/7.0
     &                + (4.0/105.0)*dmpi25*dmpk2*r4/term
     &                + (8.0/21.0)*dmpi24*dmpk2*r3/term
     &                + (12.0/7.0)*dmpi23*dmpk2*r2/term
     &                + 4.0*dmpi22*dmpk2*r/term
     &                + 4.0*dmpi2*dmpk2/term) * expi
            if (rorder .ge. 11) then
               r6 = r5 * r
               dmpi26 = dmpi25 * dmpi2
               dmpk26 = dmpk25 * dmpk2
               d5s = (dmpi2*dmpk25*r6/945.0
     &                   + (2.0/189.0)*dmpi2*dmpk24*r5
     &                   + dmpi2*dmpk23*r4/21.0
     &                   + dmpi2*dmpk22*r3/9.0
     &                   + dmpi2*dmpk2*r2/9.0
     &                   - (4.0/945.0)*dmpi2*dmpk26*r5/term
     &                   - (4.0/63.0)*dmpi2*dmpk25*r4/term
     &                   - (4.0/9.0)*dmpi2*dmpk24*r3/term
     &                   - (16.0/9.0)*dmpi2*dmpk23*r2/term
     &                   - 4.0*dmpi2*dmpk22*r/term
     &                   - 4.0*dmpi2*dmpk2/term) * expk
     &             + (dmpi25*dmpk2*r6/945.0
     &                   + (2.0/189.0)*dmpi24*dmpk2*r5
     &                   + dmpi23*dmpk2*r4/21.0
     &                   + dmpi22*dmpk2*r3/9.0
     &                   + dmpi2*dmpk2*r2/9.0
     &                   + (4.0/945.0)*dmpi26*dmpk2*r5/term
     &                   + (4.0/63.0)*dmpi25*dmpk2*r4/term
     &                   + (4.0/9.0)*dmpi24*dmpk2*r3/term
     &                   + (16.0/9.0)*dmpi23*dmpk2*r2/term
     &                   + 4.0*dmpi22*dmpk2*r/term
     &                   + 4.0*dmpi2*dmpk2/term) * expi
            end if
         end if
      end if
c
c     convert partial derivatives into full derivatives
c
        s =   s* rr1
       ds =  ds* rr3
      d2s = d2s* rr5
      d3s = d3s* rr7
      dmpik(1) = 0.5* pre*  s*s
      dmpik(2) =      pre*  s*ds
      dmpik(3) =      pre* (s*d2s + ds*ds)
      dmpik(4) =      pre* (s*d3s + 3.0*ds*d2s)

      if (rorder .ge. 9) then
         d4s      = d4s * rr9
         dmpik(5) = pre * (s*d4s + 4.0*ds*d3s + 3.0*d2s*d2s)
         if (rorder .ge. 11) then
            d5s = d5s * rr11
            dmpik(6) = pre * (s*d5s + 5.0*ds*d4s + 10.0*d2s*d3s)
         end if
      end if
      end
#endif
