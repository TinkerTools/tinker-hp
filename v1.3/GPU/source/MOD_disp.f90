!
!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ################################################################
!     ##                                                            ##
!     ##  module disp  --  damped dispersion for current structure  ##
!     ##                                                            ##
!     ################################################################
!
!
!     ndisp     total number of dispersion sites in the system
!     ndisploc  local number of dispersion sites in the system
!     ndispbloc local+neighbor number of dispersion sites in the system
!     ndisploc localnl number of dispersion sites in the system
!     ndisprecloc local reciprocal number of dispersion sites in the system
!     idisp     number of the atom for each dispersion site
!     displist     number of the dispersion site for each atom
!     nbdisp  number of dispersion sites before each atom
!     csixpr    pairwise sum of C6 dispersion coefficients
!     csix      C6 dispersion coefficient value at each site
!     adisp     alpha dispersion damping value at each site
!     winidisp  window object corresponding to idisp
!     wincsix   window object corresponding to csix
!     winadisp  window object corresponding to adisp
!     windisplist  window object corresponding to displist
!     winnbdisp  window object corresponding to nbdisp
!     displocnl glob-locnl dispersion correspondance
!     disprecloc  global-local reciprocal dispersion correspondance
!
!
#include "tinker_macro.h"
module disp
   implicit none
   integer ndisp,ndisploc,ndispbloc,ndisplocnl,ndisprecloc
   integer,pointer :: idisp(:),displist(:),nbdisp(:)
   integer,allocatable,target :: displocnl(:),disprecloc(:)
   real(t_p) csixpr
   real(t_p), pointer :: csix(:)
   real(t_p), pointer :: adisp(:)
   integer winidisp,wincsix,winadisp,windisplist,winnbdisp
end
