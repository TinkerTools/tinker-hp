!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     #####################################################################
!     ##                                                                 ##
!     ##  module ktrtor  --  forcefield parameters for torsion-torsions  ##
!     ##                                                                 ##
!     #####################################################################
!
!
!
module ktrtor
   implicit none
   integer :: maxntt !<maximum number of torsion-torsion parameter entries
   integer :: maxtgrd !<maximum dimension of torsion-torsion spline grid
   integer :: maxtgrd2 !<maximum dimension of torsion-torsion spline grid points
   parameter (maxntt=100)
   parameter (maxtgrd=30)
   parameter (maxtgrd2=maxtgrd*maxtgrd)
   integer :: tnx(maxntt) !<number of columns in torsion-torsion spline grid
   integer :: tny(maxntt) !<number of rows in torsion-torsion spline grid
   real*8 :: ttx(maxtgrd,maxntt) !<angle values for first torsion of spline grid
   real*8 :: tty(maxtgrd,maxntt) !<angle values for second torsion of spline grid
   real*8 :: tbf(maxtgrd2,maxntt) !<function values at points on spline grid
   real*8 :: tbx(maxtgrd2,maxntt) !<gradient over first torsion of spline grid
   real*8 :: tby(maxtgrd2,maxntt) !<gradient over second torsion of spline grid
   real*8 :: tbxy(maxtgrd2,maxntt) !<Hessian cross components over spline grid
   character*20 :: ktt(maxntt) !<string of torsion-torsion atom classes
   save
end
