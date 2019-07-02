c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #####################################################################
c     ##                                                                 ##
c     ##  module atmlst    --  local geometry terms involving each atom  ##
c     ##                                                                 ##
c     #####################################################################
c
c
c     bndlist   list of the bond numbers involving each atom
c     anglist   list of the angle numbers centered on each atom
c     
c     bndglob   local - global bond  correspondance 
c     angleglob local - global angle correspondance 
c     torsglob  local - global torsion correspondance 
c     bitorsglob local - global bitorsion correspondance 
c     strbndglob local - global strech bending correspondance 
c     ureyglob  local - global urey bradley correspondance 
c     angangglob local - global angle correspondance 
c     opbendglob local - global out of plane bending correspondance 
c     opdistglob local - global out of plane distance correspondance 
c     impropglob local - global improper dihedral correspondance 
c     imptorglob local - global improper torsion correspondance 
c     pitorsglob local - global pi torsion correspondance 
c     strtorglob local - global strech torsion correspondance 
c     tortorglob local - global torsion torsion correspondance 
c     vdwglob   local - global vdw correspondance 
c     poleglob  local - global direct multipole correspondance 
c     polerecglob local - global reciprocal multipole correspondance 
c     chgglob  local - global direct charge correspondance 
c     chgrecglob local - global reciprocal charge correspondance 
c
c     molculeglob local - global molecule correspondance 
c     npfixglob local - global position restrains correspondance 
c     ndfixglob local - global distance restrains correspondance 
c     nafixglob local - global angle restrains correspondance 
c     ntfixglob local - global torsion restrains correspondance 
c     ngfixglob local - global group restrains correspondance 
c     nchirfixglob local - global chiral restrains correspondance 
c     ratglob local - global constrains correspondance 
c
c     chgglobnl  localnl - global direct charge correspondance 
c     vdwglobnl  localnl - global vdw correspondance 
c     poleglobbnl  localnl - global direct multipole correspondance 
c
      module atmlst
      implicit none
      integer, pointer :: bndlist(:,:),anglist(:,:)
      integer, allocatable :: bndglob(:),angleglob(:),torsglob(:)
      integer, allocatable :: bitorsglob(:),strbndglob(:)
      integer, allocatable :: ureyglob(:),angangglob(:)
      integer, allocatable :: opbendglob(:),opdistglob(:)
      integer, allocatable :: impropglob(:),imptorglob(:)
      integer, allocatable :: pitorsglob(:),strtorglob(:)
      integer, allocatable :: tortorglob(:)
      integer, allocatable :: vdwglob(:),poleglob(:),polerecglob(:)
      integer, allocatable :: chgglob(:),chgrecglob(:)
      integer, allocatable :: molculeglob(:)
      integer, allocatable :: npfixglob(:),ndfixglob(:),nafixglob(:)
      integer, allocatable :: ntfixglob(:),ngfixglob(:),nchirglob(:)
      integer, allocatable :: ratglob(:)
      integer, allocatable :: chgglobnl(:),vdwglobnl(:),poleglobnl(:)
      save
      end
