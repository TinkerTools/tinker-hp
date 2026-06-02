!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ####################################################################
!     ##                                                                ##
!     ##  module neigh  --  pairwise neighbor list indices and storage  ##
!     ##                                                                ##
!     ####################################################################
!
!
!     ineigup     !<number of time steps between each neighbor list update
!     ncell_tot   !<number of cells corresponding to the unit cell: cell-list method
!
!     nvlst       !<number of sites in list for each vdw site
!     vlst        !<site numbers in neighbor list of each vdw site
!     nelst       !<number of sites in list for each electrostatic site
!     elst        !<site numbers in list of each electrostatic site
!     ineignl     !<localnl-global correspondance
!     neigecell   !<neighboring cells, cell-list method
!     numneigcell  !<number of neighboring cells, cell-list method
!     repartcell  !<index of the cell corresponding to each atom, cell-list method
!     cell_len    !<number of atoms in each cell, cell-list method
!     lcell_nl    !<list of cells to be included for neighbor search, cell-list method
!     ncell_nl    !<number of cells to be included for neighbor search, cell-list method
!     indcell, bufbegcell localcell-global correspondance
!     xbegcell, ybegcell, zbegcell   x,y,z coordinates of the beginning of  the cells
!     xendcell, yendcell, zendcell   x,y,z coordinates of the ending of  the cells
!
!     lbuffer     !<width of the neighbor list buffer region
!     lbuf2       !<square of half the neighbor list buffer width
!     cbuf2       !<square of short range charge cutoff plus neighbor list buffer
!     vbuf2       !<square of vdw cutoff plus neighbor list buffer
!     mbuf2       !<square of multipole cutoff plus neighbor list buffer
!     dbuf2       !<square of dispersion cutoff plus neighbor list buffer
!     cshortbuf2       !<square of charge cutoff plus neighbor list buffer
!     vshortbuf2       !<square of short range vdw cutoff plus neighbor list buffer
!     mshortbuf2       !<square of short range multipole cutoff plus neighbor list buffer
!     dshortbuf2       !<square of short range dispersion cutoff plus neighbor list buffer
!     torquebuf2  !<square of torque cutoff plus neighbor list buffer
!     torqueshortbuf2  !<square of short range torque cutoff plus neighbor list buffer
!
!
module neigh
   implicit none
   integer :: ineigup !<number of time steps between each neighbor list update
   integer :: ncell_tot !<number of cells corresponding to the unit cell: cell-list method
   integer :: ncell_nl !<number of cells to be included for neighbor search, cell-list method
   integer, allocatable :: nvlst(:) !<number of sites in list for each vdw site
   integer, allocatable :: vlst(:,:) !<site numbers in neighbor list of each vdw site
   integer, allocatable :: nelst(:) !<number of sites in list for each electrostatic site
   integer, allocatable :: elst(:,:) !<site numbers in list of each electrostatic site
   integer, allocatable :: nshortvlst(:)!<number of sites in list for each short vdw site
   integer, allocatable :: shortvlst(:,:) !<site numbers in neighbor list of each short vdw site)
   integer, allocatable :: nshortelst(:)!<number of sites in list for each short electrostatic site
   integer, allocatable :: shortelst(:,:) !<site numbers in list of each short electrostatic site   )
   integer, allocatable :: ineignl(:)!<localnl-global correspondance
   integer, allocatable :: neigcell(:,:) !<neighboring cells, cell-list method
   integer, allocatable :: numneigcell(:) !<number of neighboring cells, cell-list method
   integer, allocatable :: repartcell(:) !<index of the cell corresponding to each atom, cell-list method
   integer, allocatable :: cell_len(:) !<number of atoms in each cell, cell-list method
   integer, allocatable :: indcell(:) !<localcell-global correspondance
   integer, allocatable :: bufbegcell(:) !<localcell-global correspondance
   integer, allocatable :: lcell_nl(:) !<list of cells to be included for neighbor search, cell-list method
   real*8 :: lbuffer !<width of the neighbor list buffer region
   real*8 :: lbuf2 !<square of half the neighbor list buffer width
   real*8 :: vbuf2 !<square of vdw cutoff plus neighbor list buffer
   real*8 :: cbuf2 !<square of short range charge cutoff plus neighbor list buffer
   real*8 :: mbuf2 !<square of multipole cutoff plus neighbor list buffer
   real*8 :: torquebuf2 !<square of torque cutoff plus neighbor list buffer
   real*8 :: dbuf2 !<square of dispersion cutoff plus neighbor list buffer
   real*8 :: vshortbuf2 !<square of short range vdw cutoff plus neighbor list buffer
   real*8 :: cshortbuf2 !<square of charge cutoff plus neighbor list buffer
   real*8 :: mshortbuf2 !<square of short range multipole cutoff plus neighbor list buffer
   real*8 :: torqueshortbuf2 !<square of short range torque cutoff plus neighbor list buffer
   real*8 :: dshortbuf2 !<square of short range dispersion cutoff plus neighbor list buffer
   save
end
