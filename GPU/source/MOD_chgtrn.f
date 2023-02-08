c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ################################################################
c     ##                                                            ##
c     ##  module chgtrn  --  charge transfer for current structure  ##        
c     ##                                                            ##
c     ################################################################
c
c
c     nct       total number of dispersion sites in the system
c     chgct     charge for charge transfer at each multipole site
c     dmpct     charge transfer damping factor at each multipole site
c     winchgct  window associated to chgct array
c     windmpct  window associated to dmpct array
c
c
#include "tinker_macro.h"
      module chgtrn
      implicit none
      integer nct
      real(t_p), pointer :: chgct(:)
      real(t_p), pointer :: dmpct(:)
      integer :: winchgct,windmpct
      save
      end
