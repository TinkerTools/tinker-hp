c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ####################################################################
c     ##                                                                ##
c     ##  module polpot  --  specifics of polarization functional form  ##
c     ##                                                                ##
c     ####################################################################
c
c
c     poleps    induced dipole convergence criterion (rms Debyes/atom)
c     p2scale   scale factor for 1-2 polarization energy interactions
c     p3scale   scale factor for 1-3 polarization energy interactions
c     p4scale   scale factor for 1-4 polarization energy interactions
c     p5scale   scale factor for 1-5 polarization energy interactions
c     p41scale  additional factor for 1-4 intragroup polarization
c     d1scale   scale factor for intra-group direct induction
c     d2scale   scale factor for 1-2 group direct induction
c     d3scale   scale factor for 1-3 group direct induction
c     d4scale   scale factor for 1-4 group direct induction
c     u1scale   scale factor for intra-group mutual induction
c     u2scale   scale factor for 1-2 group mutual induction
c     u3scale   scale factor for 1-3 group mutual induction
c     u4scale   scale factor for 1-4 group mutual induction
c     politer   maximum number of induced dipole SCF iterations
c     poltyp    type of polarization potential (direct or mutual)
c     polalg    algorithm to be used to solve the induced dipoles
c               1) Preconditioned Conjugate Gradient
c               2) Jacobi/DIIS
c               3) TCG
c               5) DC-Jacobi/DIIS
c     polalgshort algorithm to be used to solve the short range induced dipoles in respa1 integrators
c               1) Preconditioned Conjugate Gradient
c               2) Jacobi/DIIS
c               3) TCG
c               5) DC-Jacobi/DIIS
c     polprt    printing flag for induce.
c               1) Print convergence information after the final iteration
c               2) Print convergence information at each iteration
c               3) Print the converged induced dipoles
c               4) Also print the auxiliary dipoles (uinp, uint)
c     polgsf    whether to use dipoles from a previous iteration as a
c               guess (1) or to use direct field dipoles (0)
c
c     tcg related keywords :
c               - tcgorder : self explanatory
c               - tcgprec : use a (diagonal) preconditioner
c               - tcgguess : use 'alpha*E' as a guess 
c               - tcgpeek : use a peek-step
c               - tcgomega : omega value for the peek-step
c               -*short : idem for short range tcg polarization in respa1 integrators
c               - tcgomegafit : true if omega has to be fitted
c               - omegafitstep : if TRUE, tcgomega has to be refitted at the current step
c               - residue : 3,N vector; contains residue of the previous
c                    tcg iteration
c               - munp : TCG induced dipoles without the peek step part.
c                    Useful for omega refitting.
c               - efres : electric field used to compute E_TCG, used for
c                    omega refitting.
c               - omegafitfreq : each 'omegafitfreq', refit the
C                    peek-step's omega to the energy from a CG
c               - epCG : polarization energy issued from CG
c
c
#include "tinker_precision.h"
      module polpot
      implicit none

      ! Polarisation Solvers present in Tinker-HP
      enum,bind(C)
      enumerator pcg_SId=1
      enumerator jacobi_SId
      enumerator tcg_SId
      enumerator dc_diis_SId=5
      enumerator step_pcg_SId=60
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
      real(t_p) d1scale,d2scale
      real(t_p) d3scale,d4scale
      real(t_p) u1scale,u2scale
      real(t_p) u3scale,u4scale
      real(t_p) tcgomega,tcgomegashort
      character*6 poltyp
      logical :: omegafitstep,tcgomegafit
      integer :: omegafitfreq
      real(t_p), allocatable :: residue(:,:), munp(:,:), efres(:,:)
      real(t_p):: epCG

      integer  ,allocatable:: ucorrect_ik(:),dpcorrect_ik(:),
     &                        dpucorrect_ik(:)
      real(t_p),allocatable:: ucorrect_scale(:),dpcorrect_scale(:),
     &                        dpucorrect_scale(:)

!$acc declare create(p2scale,p3scale,p4scale,p41scale,p5scale,
!$acc& u1scale,u2scale,u3scale,u4scale,d1scale,d2scale,d3scale,d4scale)
      end
