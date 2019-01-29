c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ######################################################################
c     ##                                                                  ##
c     ##  module vdw  --  van der Waals parameters for current structure  ##
c     ##                                                                  ##
c     ######################################################################
c
c
c     radmin     minimum energy distance for each atom class pair
c     epsilon    well depth parameter for each atom class pair
c     radmin4    minimum energy distance for 1-4 interaction pairs
c     epsilon4   well depth parameter for 1-4 interaction pairs
c     radhbnd    minimum energy distance for hydrogen bonding pairs
c     epshbnd    well depth parameter for hydrogen bonding pairs
c     kred       value of reduction factor parameter for each atom
c     ired       attached atom from which reduction factor is applied
c     nvdw       total number van der Waals active sites in the system
c     nvdwloc    local number van der Waals active sites in the system
c     nvdwbloc   local+neighbors number van der Waals active sites in the system
c     ivdw       number of the atom for each van der Waals active site
c     jvdw       type or class index into vdw parameters for each atom
c     nvt        number of distinct vdw types/classes in the system
c     ivt        type/class index for each distinct vdw type or class
c     jvt        frequency of each vdw type or class in the system
c     nbvdw      number of 'vdw' atoms before each atom
c     vdwlocnl   glob-locnl vdw correspondance
c
c
      module vdw
      implicit none
      integer nvdw,nvt,nvdwloc,nvdwbloc,nvdwlocnl
      integer, pointer :: jvdw(:),ivdw(:),ired(:)
      integer, pointer :: ivt(:),jvt(:),nbvdw(:)
      integer, pointer :: vdwlocnl(:)
      real*8, pointer :: radmin(:,:),epsilon(:,:)
      real*8, pointer :: radmin4(:,:),epsilon4(:,:)
      real*8, pointer :: radhbnd(:,:),epshbnd(:,:)
      real*8, pointer :: kred(:)
      save
      end
