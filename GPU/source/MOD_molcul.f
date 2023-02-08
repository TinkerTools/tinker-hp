c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #####################################################################
c     ##                                                                 ##
c     ##  module molcul  --  individual molecules within current system  ##
c     ##                                                                 ##
c     #####################################################################
c
c
c     molmass   molecular weight for each molecule in the system
c     winmolmass    window object corresponding to molmass
c     totmass   total weight of all the molecules in the system
c     nmol      total number of separate molecules in the system
c     nmoleloc  local number of separate molecules in the system
c     kmol      contiguous list of the atoms in each molecule
c     winkmol    window object corresponding to kmol
c     imol      first and last atom of each molecule in the list
c     winimol    window object corresponding to imol
c     molcule   number of the molecule to which each atom belongs
c     winmolcule    window object corresponding to molcule
c
c
#include "tinker_macro.h"
      module molcul
      implicit none
      integer nmol,nmoleloc
      integer :: winmolcule,winkmol,winimol
      integer :: winmolmass
      integer, pointer :: molcule(:),kmol(:),imol(:,:)
      real(r_p) totmass
      real(t_p), pointer :: molmass(:)
      end
