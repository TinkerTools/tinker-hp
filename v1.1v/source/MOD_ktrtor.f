c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #####################################################################
c     ##                                                                 ##
c     ##  module ktrtor  --  forcefield parameters for torsion-torsions  ##
c     ##                                                                 ##
c     #####################################################################
c
c
c     maxntt    maximum number of torsion-torsion parameter entries
c     maxtgrd   maximum dimension of torsion-torsion spline grid
c     maxtgrd2  maximum number of torsion-torsion spline grid points
c
c     ttx       angle values for first torsion of spline grid
c     tty       angle values for second torsion of spline grid
c     tbf       function values at points on spline grid
c     tbx       gradient over first torsion of spline grid
c     tby       gradient over second torsion of spline grid
c     tbxy      Hessian cross components over spline grid
c     tnx       number of columns in torsion-torsion spline grid
c     tny       number of rows in torsion-torsion spline grid
c     ktt       string of torsion-torsion atom classes
c
c
      module ktrtor
      implicit none
      integer maxntt,maxtgrd,maxtgrd2
      parameter (maxntt=100)
      parameter (maxtgrd=30)
      parameter (maxtgrd2=maxtgrd*maxtgrd)
      integer tnx(maxntt),tny(maxntt)
      real*8 ttx(maxtgrd,maxntt),tty(maxtgrd,maxntt)
      real*8 tbf(maxtgrd2,maxntt)
      real*8 tbx(maxtgrd2,maxntt),tby(maxtgrd2,maxntt)
      real*8 tbxy(maxtgrd2,maxntt)
      character*20 ktt(maxntt)
      save 
      end
