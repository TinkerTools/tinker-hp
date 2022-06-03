c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #################################################################
c     ##                                                             ##
c     ##  module group  --  partitioning of system into atom groups  ##
c     ##                                                             ##
c     #################################################################
c
c
c     grpmass     total mass of all the atoms in each group
c     wgrp        weight for each set of group-group interactions
c     ngrp        total number of atom groups in the system
c     kgrp        contiguous list of the atoms in each group
c     igrp        first and last atom of each group in the list
c     grplist     number of the group to which each atom belongs
c     use_group   flag to use partitioning of system into groups
c     use_intra   flag to include only intragroup interactions
c     use_inter   flag to include only intergroup interactions
c
c     natgroup    number of atoms involved in a group
c     globglobgroup   associated indexes
c     loclocgroup   global-group correspondance
c
c     nlocatgroup    number of local atoms involved in a group
c     globgroup   associated group indexes
c     locgroup   global group - local group correspondance
c
c
c     npolegroup  number of multipoles involved in a group
c     ipolegroup  associated multipoles indexes
c     pollistgroup  aotms-multipoles correspondance in the group frame
c     npolelocgroup  number of local multipoles involved in a group
c     poleglobgroup   local-global correspondance for the multipoles in the group frame
c     polelocgroup  global-local correspondance in the group frame
c      
c
c
c     domlengroup number of atoms in the group per domain
c     bufbeggroup "bufbeg" equivalent in the group frame
c     domlenpolegroup number of multipoles in the group per domain
c     bufbegpolegroup "bufbegpole" equivalent in the group frame
c     uindgroup,uinpgroup induced dipoles associated to the group
c     epgroup  polarization energy associated to the group
c     depgroup  polarization energy derivatives associated to the group
c     vir_group virial associated to the polarization energy of the group
c
c
#include "tinker_precision.h"
#include "tinker_types.h"

      module group
      use sizes
      implicit none
      integer ngrp
      integer, allocatable :: kgrp(:),grplist(:)
      integer, allocatable :: igrp(:,:)
      real(t_p), allocatable :: grpmass(:),wgrp(:,:)

      logical use_group
      logical use_intra
      logical use_inter
      integer :: natgroup,nlocatgroup
      integer :: npolegroup, npolelocgroup
      integer, allocatable :: globglobgroup(:),loclocgroup(:)
      integer, allocatable :: globgroup(:),locgroup(:)
      integer, allocatable :: poleglobgroup(:),polelocgroup(:)
      integer, allocatable :: domlengroup(:),ipolegroup(:)
      integer, allocatable :: domlenpolegroup(:),bufbegpolegroup(:)
      integer, allocatable :: bufbeggroup(:),pollistgroup(:)
      real(t_p), allocatable :: uindgroup(:,:),uinpgroup(:,:)
      mdyn_rtyp, pointer :: depgroup(:,:)
      real(r_p) vir_group(3,3)
      ener_rtyp epgroup
      end
