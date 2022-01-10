c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ####################################################################
c     ##                                                                ##
c     ##  module neigh  --  pairwise neighbor list indices and storage  ## c     ##                                                                ##
c     ####################################################################
c
c
c     ineigup     number of time steps between each neighbor list update
c     ncell_tot   number of cells corresponding to the unit cell: cell-list method
c
c     nvlst       number of sites in list for each vdw site
c     vlst        site numbers in neighbor list of each vdw site
c     nelst       number of sites in list for each electrostatic site
c     elst        site numbers in list of each electrostatic site
c     nbondlst    number of bonds in list for each bond (repulsion bond-bond interaction)
c     bondlst     bond numbers in list of each bond (repulsion bond-bond interaction)
c     nbondlplst  number of sites in list for each bond (repulsion bond-lp interaction)
c     bondlplst   site numbers in list of each bond (repulsion bond-lp interaction)
c     nlplplst    number of lp in list for each lp (repulsion lp-lp interaction)
c     lplplst     lp numbers in list of each lp (repulsion lp-lp interaction)
c     natatlst    number of sites in list for each site (dispersion atom-atom interaction)
c     atatlst     site numbers in list of each site (dispersion atom-atom interaction)
c     nlpacclst   number of acceptors in list for each lp (charge transfer interaction)
c     lpacclst    acceptor numbers in list of each lp (charge transfer interaction)
c     naccpotlst  number of sites in list for each acceptor (elec potential for ct interaction)
c     accpotlst   sites numbers in list of each acceptor (elec potential for ct interaction)
c     lppotlst    number of sites in list for each lp (elec potential for ct interaction)
c     lppotlst    sites numbers in list of each lp (elec potential for ct interaction)
c
c
c     ineignl     localnl-global correspondance
c     locnl     global-localnlcorrespondance
c     neigecell   neighboring cells, cell-list method
c     numneigcell  number of neighboring cells, cell-list method
c     repartcell  index of the cell corresponding to each atom, cell-list method
c     cell_len    number of atoms in each cell, cell-list method
c     indcell, bufbegcell localcell-global correspondance
c     xbegcell, ybegcell, zbegcell   x,y,z coordinates of the beginning of  the cells
c     xendcell, yendcell, zendcell   x,y,z coordinates of the ending of  the cells
c
c     lbuffer     width of the neighbor list buffer region
c     lbuf2       square of half the neighbor list buffer width
c     vbuf2       square of vdw cutoff plus neighbor list buffer
c     mbuf2       square of multipole cutoff plus neighbor list buffer
c     torquebuf2  square of torque cutoff plus neighbor list buffer
c     repbuf2     square of repulsion cutoff plus neighbor list buffer
c     dispbuf2    square of dispersion cutoff plus neighbor list buffer
c     ctransferbuf2  square of charge transfer cutoff plus neighbor list buffer
c     mpolectbuf2  square of electrostatic potential for charge transfer cutoff plus neighbor list buffer
c
c
c
      module neigh
      implicit none
      integer ineigup
      integer ncell_tot
      integer, allocatable :: nvlst(:)
      integer, allocatable :: vlst(:,:)
      integer, allocatable :: nelst(:)
      integer, allocatable :: elst(:,:)
      integer, allocatable :: nbondlst(:),nbondlplst(:),nlplplst(:)
      integer, allocatable :: bondlst(:,:),bondlplst(:,:),lplplst(:,:)
      integer, allocatable :: natatlst(:),natlplst(:),nlpatlst(:)
      integer, allocatable :: atatlst(:,:),atlplst(:,:),lpatlst(:,:)
      integer, allocatable :: nlpacclst(:)
      integer, allocatable :: lpacclst(:,:)
      integer, allocatable :: naccpotlst(:,:),accpotlst(:,:,:)
      integer, allocatable :: nlppotlst(:),lppotlst(:,:)
      integer, allocatable :: ineignl(:),locnl(:)
      integer, allocatable :: neigcell(:,:),numneigcell(:),repartcell(:)
      integer, allocatable :: cell_len(:),indcell(:),bufbegcell(:)
      real*8 lbuffer,lbuf2
      real*8 vbuf2,cbuf2,mbuf2,torquebuf2
      real*8 repbuf2,dispbuf2,ctransferbuf2,mpolectbuf2
      real*8, allocatable :: xbegcell(:),ybegcell(:),zbegcell(:)
      real*8, allocatable :: xendcell(:),yendcell(:),zendcell(:)
      save
      end
