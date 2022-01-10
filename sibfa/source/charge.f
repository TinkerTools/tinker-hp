c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ####################################################################
c     ##                                                                ##
c     ##  module charge  --  partial charges for the current structure  ##
c     ##                                                                ##
c     ####################################################################
c
c
c     nion      total number of partial charges in system
c     nionloc   local number of partial charges in system
c     nionbloc   local+neighbors number of partial charges in system
c     nionlocnl  localnl number of partial charges in system
c     nionrecloc   local reciprocal number of partial charges in system
c     iion      number of the atom site for each partial charge
c     jion      neighbor generation site for each partial charge
c     kion      cutoff switching site for each partial charge
c     chglist   partial charge site for each atom (0=no charge)
c     nbchg     number of charges before each index
c     chgloc    global-local charge correspondance
c     chglocnl  global-localnl charge correspondance
c     chgrecloc  global-local reciprocal charge correspondance
c     pchg      magnitude of the partial charges (e-)
c
c
      module charge
      implicit none
      integer nion,nionloc,nionbloc,nionlocnl,nionrecloc
      integer, pointer :: iion(:)
      integer, pointer :: jion(:),kion(:)
      integer, pointer :: chglist(:)
      integer, pointer :: nbchg(:)
      integer, pointer :: chgloc(:),chglocnl(:)
      integer, pointer :: chgrecloc(:)
      real*8, pointer ::  pchg(:)
      save 
      end
