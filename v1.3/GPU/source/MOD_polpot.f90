!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ####################################################################
!     ##                                                                ##
!     ##  module polpot  --  specifics of polarization functional form  ##
!     ##                                                                ##
!     ####################################################################
!
!
!     poleps    induced dipole convergence criterion (rms Debyes/atom)
!     p2scale   scale factor for 1-2 polarization energy interactions
!     p3scale   scale factor for 1-3 polarization energy interactions
!     p4scale   scale factor for 1-4 polarization energy interactions
!     p5scale   scale factor for 1-5 polarization energy interactions
!     p41scale  additional factor for 1-4 intragroup polarization
!     p2iscale  scale factor for 1-2 intragroup polarization energy
!     p3iscale  scale factor for 1-3 intragroup polarization energy
!     p4iscale  scale factor for 1-4 intragroup polarization energy
!     p5iscale  scale factor for 1-5 intragroup polarization energy
!     d1scale   scale factor for intra-group direct induction
!     d2scale   scale factor for 1-2 group direct induction
!     d3scale   scale factor for 1-3 group direct induction
!     d4scale   scale factor for 1-4 group direct induction
!     u1scale   scale factor for intra-group mutual induction
!     u2scale   scale factor for 1-2 group mutual induction
!     u3scale   scale factor for 1-3 group mutual induction
!     u4scale   scale factor for 1-4 group mutual induction
!     w2scale   scale factor for 1-2 induced dipole interactions
!     w3scale   scale factor for 1-3 induced dipole interactions
!     w4scale   scale factor for 1-4 induced dipole interactions
!     w5scale   scale factor for 1-5 induced dipole interactions
!     politer   maximum number of induced dipole SCF iterations
!     poltyp    type of polarization potential (direct or mutual)
!     polalg    algorithm to be used to solve the induced dipoles
!               1) Preconditioned Conjugate Gradient
!               2) Jacobi/DIIS
!               3) TCG
!               5) DC-Jacobi/DIIS
!     polalgshort algorithm to be used to solve the short range induced dipoles in respa1 integrators
!               1) Preconditioned Conjugate Gradient
!               2) Jacobi/DIIS
!               3) TCG
!               5) DC-Jacobi/DIIS
!     polprt    printing flag for induce.
!               1) Print convergence information after the final iteration
!               2) Print convergence information at each iteration
!               3) Print the converged induced dipoles
!               4) Also print the auxiliary dipoles (uinp, uint)
!     polgsf    whether to use dipoles from a previous iteration as a
!               guess (1) or to use direct field dipoles (0)
!
!     tcg related keywords :
!               - tcgorder : self explanatory
!               - tcgprec : use a (diagonal) preconditioner
!               - tcgguess : use 'alpha*E' as a guess
!               - tcgpeek : use a peek-step
!               - tcgomega : omega value for the peek-step
!               -*short : idem for short range tcg polarization in respa1 integrators
!               - tcgomegafit : true if omega has to be fitted
!               - omegafitstep : if TRUE, tcgomega has to be refitted at the current step
!               - residue : 3,N vector; contains residue of the previous
!                    tcg iteration
!               - munp : TCG induced dipoles without the peek step part.
!                    Useful for omega refitting.
!               - efres : electric field used to compute E_TCG, used for
!                    omega refitting.
!               - omegafitfreq : each 'omegafitfreq', refit the
!                    peek-step's omega to the energy from a CG
!               - epCG : polarization energy issued from CG
!
!
#include "tinker_macro.h"
module polpot
   implicit none

   ! Polarisation Solvers present in Tinker-HP
   enum,bind(C)
      enumerator :: pcg_SId=1
      enumerator jacobi_SId
      enumerator tcg_SId
      enumerator :: dc_diis_SId=5
      enumerator :: step_pcg_SId=60
      enumerator step_pcg_short_SId
   end enum

   integer politer,polalg,polprt,polgsf,polff
   integer tcgorder,polalgshort
   logical tcgprec,tcgguess,tcgpeek
   integer tcgordershort
   logical tcgprecshort,tcgguessshort,tcgpeekshort
   integer n_uscale,n_dpscale,n_dpuscale
   real(r_p) poleps
   real(t_p) p2scale
   real(t_p) p3scale,p4scale
   real(t_p) p5scale,p41scale
   real(t_p) p2iscale,p3iscale
   real(t_p) p4iscale,p5iscale
   real(t_p) d1scale,d2scale
   real(t_p) d3scale,d4scale
   real(t_p) u1scale,u2scale
   real(t_p) u3scale,u4scale
   real(t_p) w2scale,w3scale
   real(t_p) w4scale,w5scale
   real(t_p) tcgomega,tcgomegashort
   character(6) poltyp
   logical omegafitstep,tcgomegafit
   integer omegafitfreq
   real(t_p), allocatable :: residue(:,:), munp(:,:), efres(:,:)
   real(t_p) epCG
   logical use_thole,use_tholed,dpequal

   integer  ,allocatable:: ucorrect_ik(:),dpcorrect_ik(:),&
      &dpucorrect_ik(:)
   real(t_p),allocatable:: ucorrect_scale(:),dpcorrect_scale(:),&
      &dpucorrect_scale(:)

!$acc declare create(p2scale,p3scale,p4scale,p41scale,p5scale, &
!$acc u1scale,u2scale,u3scale,u4scale,d1scale,d2scale,d3scale,d4scale)
end
