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
!     maxntt    maximum number of torsion-torsion parameter entries
!     maxtgrd   maximum dimension of torsion-torsion spline grid
!     maxtgrd2  maximum number of torsion-torsion spline grid points
!
!     ttx       angle values for first torsion of spline grid
!     tty       angle values for second torsion of spline grid
!     tbf       function values at points on spline grid
!     tbx       gradient over first torsion of spline grid
!     tby       gradient over second torsion of spline grid
!     tbxy      Hessian cross components over spline grid
!     tnx       number of columns in torsion-torsion spline grid
!     tny       number of rows in torsion-torsion spline grid
!     ktt       signature of torsion-torsion atom classes
!     ktt_sys   system signature of torsion-torsion atom classes
!
!
#include "tinker_macro.h"
module ktrtor
   implicit none
   integer maxntt,maxtgrd,maxtgrd2
   parameter (maxntt=100)
   parameter (maxtgrd=30)
   parameter (maxtgrd2=maxtgrd*maxtgrd)
   integer tnx(maxntt),tny(maxntt)
   real(t_p) ttx(maxtgrd,maxntt),tty(maxtgrd,maxntt)
   real(t_p) tbf(maxtgrd2,maxntt)
   real(t_p) tbx(maxtgrd2,maxntt),tby(maxtgrd2,maxntt)
   real(t_p) tbxy(maxtgrd2,maxntt)
   integer(8) ktt(maxntt)
   integer(8) ktt_sys(0:maxntt)
   logical,private:: data_on_device=.false.
contains

   subroutine malloc_device_data
      if (.not.data_on_device) then
!$acc enter data create(tbxy,tbf,tbx,tby,ttx,tty,tnx,tny)
         data_on_device=.true.
      end if
   end subroutine
   subroutine free_device_data
      if (data_on_device) then
!$acc exit data delete(tbxy,tbf,tbx,tby,ttx,tty,tnx,tny)
         data_on_device=.false.
      end if
   end subroutine

end
