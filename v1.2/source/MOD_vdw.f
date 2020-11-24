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
c     winradmin    window object corresponding to radmin
c     epsilon    well depth parameter for each atom class pair
c     winepsilon    window object corresponding to epsilon
c     radmin4    minimum energy distance for 1-4 interaction pairs
c     winradmin4    window object corresponding to radmin4
c     epsilon4   well depth parameter for 1-4 interaction pairs
c     winepsilon4    window object corresponding to epsilon4
c     radhbnd    minimum energy distance for hydrogen bonding pairs
c     winradhbnd    window object corresponding to radhbnd
c     epshbnd    well depth parameter for hydrogen bonding pairs
c     winepshbnd    window object corresponding to epshbnd
c     kred       value of reduction factor parameter for each atom
c     winkred    window object corresponding to kred
c     ired       attached atom from which reduction factor is applied
c     winired    window object corresponding to ired
c     nvdw       total number van der Waals active sites in the system
c     nvdwloc    local number van der Waals active sites in the system
c     nvdwbloc   local+neighbors number van der Waals active sites in the system
c     ivdw       number of the atom for each van der Waals active site
c     winivdw    window object corresponding to ivdw
c     jvdw       type or class index into vdw parameters for each atom
c     winjvdw    window object corresponding to jvdw
c     nvt        number of distinct vdw types/classes in the system
c     ivt        type/class index for each distinct vdw type or class
c     jvt        frequency of each vdw type or class in the system
c     nbvdw      number of 'vdw' atoms before each atom
c     winnbvdw    window object corresponding to nbvdw
c     vdwlocnl   glob-locnl vdw correspondance
c
c
      module vdw
      implicit none
      integer nvdw,nvt,nvdwloc,nvdwbloc,nvdwlocnl
      integer, pointer :: jvdw(:),ivdw(:),ired(:)
      integer, pointer :: ivt(:),jvt(:),nbvdw(:)
      integer, allocatable :: vdwlocnl(:)
      integer :: winjvdw,winivdw,winired,winivt
      integer :: winjvt,winnbvdw
      real*8, pointer :: radmin(:,:),epsilon(:,:)
      real*8, pointer :: radmin4(:,:),epsilon4(:,:)
      real*8, pointer :: radhbnd(:,:),epshbnd(:,:)
      real*8, pointer :: kred(:)
      integer :: winradmin,winepsilon,winradmin4
      integer :: winepsilon4,winradhbnd,winepshbnd
      integer :: winkred
      save
      end
