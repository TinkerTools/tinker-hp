c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #####################################################################
c     ##                                                                 ##
c     ##  module freeze  --  definition of holonomic RATTLE constraints  ##
c     ##                                                                 ##
c     #####################################################################
c
c
c     rateps       convergence tolerance for holonomic constraints
c     krat         ideal distance value for holonomic constraint
c     nrat         number of holonomic distance constraints to apply
c     nratx        number of atom group holonomic constraints to apply
c     irat         atom numbers of atoms in a holonomic constraint
c     iratx        group number of group in a holonomic constraint
c     kratx        spatial constraint type (1=plane, 2=line, 3=point)
c     ratimage     flag to use minimum image for holonomic constraint
c     use_rattle   logical flag to set use of holonomic contraints
c
c     nratloc      number of local constraints
c     nratbloc      number of local constraints + neighboring constraints
c
c     buflenrat1   length of buffer of constraints to send at each rattle iteration
c     buflenrat2   length of buffer of constraints to receive at each rattle iteration
c     buflenrat1   first indexes of buffer of constraints to send at each rattle iteration
c     buflenrat2   first indexes of buffer of constraints to receive at each rattle iteration
c     bufrat1   list of indexes of buffer of constraints to send at each rattle iteration
c     bufrat2   list of indexes of buffer of constraints to receive at each rattle iteration
c
c
      module freeze
      implicit none
      integer nrat,nratx
      integer, pointer :: irat(:,:),iratx(:)
      integer, pointer ::  kratx(:)
      !DIR$ ATTRIBUTES ALIGN:64:: buflenrat1,buflenrat2
      integer, allocatable :: buflenrat1(:),buflenrat2(:)
      !DIR$ ATTRIBUTES ALIGN:64:: bufbegrat1,bufbegrat2
      integer, allocatable :: bufbegrat1(:),bufbegrat2(:)
      !DIR$ ATTRIBUTES ALIGN:64:: bufrat1   ,bufrat2   
      integer, allocatable :: bufrat1   (:),bufrat2   (:)
      integer nratloc,nratbloc
      real*8, pointer :: krat(:)
      real*8 rateps
      logical, pointer :: ratimage(:)
      logical use_rattle
      save
      end
