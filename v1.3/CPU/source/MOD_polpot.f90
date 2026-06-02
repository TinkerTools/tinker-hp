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
module polpot
   implicit none
   integer :: politer !<maximum number of induced dipole SCF iterations
   integer :: polalg !<algorithm to be used to solve the induced dipoles
                !<1) Preconditioned Conjugate Gradient
                !<2) Jacobi/DIIS
                !<3) TCG
                !<5) DC-Jacobi/DIIS
   integer :: polprt !<printing flag for induce.
               !<1) Print convergence information after the final iteration
               !<2) Print convergence information at each iteration
               !<3) Print the converged induced dipoles
               !<4) Also print the auxiliary dipoles (uinp, uint)
   integer :: polgsf !<whether to use dipoles from a previous iteration as a guess (1) or to use direct field dipoles (0)
   integer :: tcgorder !<tcgorder : order (1 or 2) of TCG solver
   integer :: polalgshort !<algorithm to be used to solve the short range induced dipoles in respa1 integrators
                !<1) Preconditioned Conjugate Gradient
                !<2) Jacobi/DIIS
                !<3) TCG
                !<5) DC-Jacobi/DIIS
   logical :: tcgprec !<tcgprec : use a (diagonal) preconditioner
   logical :: tcgguess !<tcgguess : use 'alpha*E' as a guess
   logical :: tcgpeek !<tcgpeek : use a peek-step
   integer :: tcgordershort !<tcgorder : order (1 or 2) of short range TCG solver
   logical :: tcgprecshort !<tcgprec : use a (diagonal) preconditioner for short range TCG solver
   logical :: tcgguessshort !<tcgguess : use 'alpha*E' as a guess for short range TCG solver
   logical :: tcgpeekshort !<tcgpeek : use a peek-step for short range TCG solver
   real*8 :: poleps !<induced dipole convergence criterion (rms Debyes/atom)
   real*8 :: p2scale !<scale factor for 1-2 polarization energy interactions
   real*8 :: p3scale !<scale factor for 1-3 polarization energy interactions
   real*8 :: p4scale !<scale factor for 1-4 polarization energy interactions
   real*8 :: p5scale !<scale factor for 1-5 polarization energy interactions
   real*8 :: p41scale !<additional factor for 1-4 intragroup polarization    e
   real*8 :: p2iscale !<scale factor for 1-2 intragroup polarization energy
   real*8 :: p3iscale !<scale factor for 1-3 intragroup polarization energy
   real*8 :: p4iscale !<scale factor for 1-4 intragroup polarization energy
   real*8 :: p5iscale !<scale factor for 1-5 intragroup polarization energy
   real*8 :: d1scale !<scale factor for intra-group direct induction
   real*8 :: d2scale !<scale factor for 1-2 group direct induction
   real*8 :: d3scale !<scale factor for 1-3 group direct induction
   real*8 :: d4scale !<scale factor for 1-4 group direct induction
   real*8 :: u1scale !<scale factor for intra-group mutual induction
   real*8 :: u2scale !<scale factor for 1-2 group mutual induction
   real*8 :: u3scale !<scale factor for 1-3 group mutual induction
   real*8 :: u4scale !<scale factor for 1-4 group mutual induction
   real*8 :: w2scale !<scale factor for 1-2 induced dipole interactions
   real*8 :: w3scale !<scale factor for 1-3 induced dipole interactions
   real*8 :: w4scale !<scale factor for 1-4 induced dipole interactions
   real*8 :: w5scale !<scale factor for 1-5 induced dipole interactions
   real*8 :: tcgomega !<tcgomega : omega value for the peek-step
   real*8 :: tcgomegashort !<tcgomega : omega value for the peek-step (short range)
   character*6 :: poltyp !<type of polarization (only MUTUAL implemented)
   logical :: omegafitstep !<if TRUE, tcgomega has to be refitted at the current step
   logical :: tcgomegafit !<true if omega has to be fitted
   integer :: omegafitfreq !<each 'omegafitfreq', refit the peek-step's omega to the energy from a CG
   real*8, allocatable :: residue(:,:) !<3,N vector; contains residue of the previous tcg iteration
   real*8, allocatable :: munp(:,:) !<TCG induced dipoles without the peek step part. Useful for omega refitting.
   real*8, allocatable :: efres(:,:) !<electric field used to compute E_TCG, used for omega refitting.
   real*8  :: epCG !<polarization energy issued from CG
   logical :: use_thole !<flag governing use of thole damping
   logical :: use_tholed !<flag governing use of direct thole damping
   logical :: dpequal !<flag governing single set of induced dipoles (d and p are equal)
   save
end
