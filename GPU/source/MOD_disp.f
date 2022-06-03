c
c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ################################################################
c     ##                                                            ##
c     ##  module disp  --  damped dispersion for current structure  ##        
c     ##                                                            ##
c     ################################################################
c
c
c     ndisp     total number of dispersion sites in the system
c     ndisploc  local number of dispersion sites in the system
c     ndispbloc local+neighbor number of dispersion sites in the system
c     ndisploc localnl number of dispersion sites in the system
c     ndisprecloc local reciprocal number of dispersion sites in the system
c     idisp     number of the atom for each dispersion site
c     displist     number of the dispersion site for each atom
c     nbdisp  number of dispersion sites before each atom 
c     csixpr    pairwise sum of C6 dispersion coefficients
c     csix      C6 dispersion coefficient value at each site
c     adisp     alpha dispersion damping value at each site
c     winidisp  window object corresponding to idisp
c     wincsix   window object corresponding to csix
c     winadisp  window object corresponding to adisp
c     windisplist  window object corresponding to displist
c     winnbdisp  window object corresponding to nbdisp
c     displocnl glob-locnl dispersion correspondance
c     disprecloc  global-local reciprocal dispersion correspondance
c
c
#include "tinker_precision.h"
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
