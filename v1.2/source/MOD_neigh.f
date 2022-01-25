c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ####################################################################
c     ##                                                                ##
c     ##  module neigh  --  pairwise neighbor list indices and storage  ##
c     ##                                                                ##
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
c     ineignl     localnl-global correspondance
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
c     cbuf2       square of short range charge cutoff plus neighbor list buffer
c     vbuf2       square of vdw cutoff plus neighbor list buffer
c     mbuf2       square of multipole cutoff plus neighbor list buffer
c     dbuf2       square of dispersion cutoff plus neighbor list buffer
c     cshortbuf2       square of charge cutoff plus neighbor list buffer
c     vshortbuf2       square of short range vdw cutoff plus neighbor list buffer
c     mshortbuf2       square of short range multipole cutoff plus neighbor list buffer
c     dshortbuf2       square of short range dispersion cutoff plus neighbor list buffer
c     torquebuf2  square of torque cutoff plus neighbor list buffer
c     torqueshortbuf2  square of short range torque cutoff plus neighbor list buffer
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
      integer, allocatable :: nshortvlst(:)
      integer, allocatable :: shortvlst(:,:)
      integer, allocatable :: nshortelst(:)
      integer, allocatable :: shortelst(:,:)
      integer, allocatable :: ineignl(:)
      integer, allocatable :: neigcell(:,:),numneigcell(:),repartcell(:)
      integer, allocatable :: cell_len(:),indcell(:),bufbegcell(:)
      real*8 lbuffer,lbuf2
      real*8 vbuf2,cbuf2,mbuf2,torquebuf2,dbuf2
      real*8 vshortbuf2,cshortbuf2,mshortbuf2,torqueshortbuf2
      real*8 dshortbuf2
      real*8, allocatable :: xbegcell(:),ybegcell(:),zbegcell(:)
      real*8, allocatable :: xendcell(:),yendcell(:),zendcell(:)
      save
      end
