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
c     vblst       site numbers block in neighbor list of each vdw site block
c     ivblst      first block pair of neighbor list C2
c     vb_tag      tag on each vdw block atoms  (to be use in mpi midpointimage)
c     nelst       number of sites in list for each electrostatic site
c     nelstc      number of sites under the cutoff in list for each electrostatic site
c     nshortelstc      number of sites under the cutoff in short range list for each electrostatic site
c     elst        site numbers in list of each electrostatic site
c     ineignl     localnl-global correspondance
c     neigecell   neighboring cells, cell-list method
c     numneigcell number of neighboring cells, cell-list method
c     repartcell  index of the cell corresponding to each atom, cell-list method
c     cell_len    number of atoms in each cell, cell-list method
c     cell_len1   number of atoms in each cell, cell-list method when cells are small enough
c     cell_lenr   number of atoms in each cell, use in couple to cell_len for cell spliting
c     cell_len2   number of atoms in each cell, a partner to cell_len
c     cell_len2r   number of atoms in each cell, a partner to cell_len
c     cell_scan   partial sum of cell_len
c     indcell,    bufbegcell localcell-global correspondance
c     xbegcell, ybegcell, zbegcell   x,y,z coordinates of the beginning of  the cells
c     xendcell, yendcell, zendcell   x,y,z coordinates of the ending of  the cells
c
c     cell*_order cell numbering array
c     cell*_glob    local atoms index following cell numbering
c     cell*_x      x position atoms renumber from cell_reorder
c     cell*_y      y position atoms renumber from cell_reorder
c     cell*_z      z position atoms renumber from cell_reorder
c
c
c     nx_cell     number of boxes in x dimension
c     ny_cell     number of boxes in y dimension
c     nz_cell     number of boxes in z dimension
c     bvx_cell     number of boxes in x dimension to reach vdw_cutoff+buffer edge
c     bvy_cell     number of boxes in y dimension to reach vdw_cutoff+buffer edge
c     bvz_cell     number of boxes in z dimension to reach vdw_cutoff+buffer edge
c     bex_cell     number of boxes in x dimension to reach elec_cutoff+buffer edge
c     bey_cell     number of boxes in y dimension to reach elec_cutoff+buffer edge
c     bez_cell     number of boxes in z dimension to reach elec_cutoff+buffer edge
c     bcx_cell     number of boxes in x dimension to reach chg_cutoff+buffer edge
c     bcy_cell     number of boxes in y dimension to reach chg_cutoff+buffer edge
c     bcz_cell     number of boxes in z dimension to reach chg_cutoff+buffer edge
c     lbuffer     width of the neighbor list buffer region
c     lbuf2       square of half the neighbor list buffer width
c     vbuf2       square of vdw cutoff plus neighbor list buffer
c     mbuf2       square of multipole cutoff plus neighbor list buffer
c     cshortbuf2       square of charge cutoff plus neighbor list buffer
c     vshortbuf2       square of short range vdw cutoff plus neighbor list buffer
c     mshortbuf2       square of short range multipole cutoff plus neighbor list buffer
c     cbuf2       square of short range charge cutoff plus neighbor list buffer
c     torquebuf2  square of torque cutoff plus neighbor list buffer
c     torqueshortbuf2  square of short range torque cutoff plus neighbor list buffer
c     auto_lst    logical flag to switch from automatic to manual neighbor list construction
c     rebuild_lst logical flag to trigger neighbor list reconstruction
c     dovlst      logical flag to rebuild vdw neighbor list
c     domlst      logical flag to rebuild multipole neighbor list
c     doclst      logical flag to rebuild charge neighbor list
c
c     matb_lst    adjacency matrix of pairwise blocks interactions
c     buffMatb    adjacency matrix maximum storage size (Go)
c     szoMatb     sizeof(matb_lst)
c     szMatb      number of element in matb_lst's line
c     nbMatb      number of blocks stored by matb_lst
c     niterMatb   number of iteration required to build Adjacency list
c     offSet[l][r]Mb  BlockId offset left and right stored by Matb
c
c     vlst_enable  logical flag to control vlst  construct
c     vlst2_enable logical flag to control vblst construct 
c     mlst_enable  logical flag to control elst  construct
c     mlst2_enable logical flag to control eblst construct 
c
c
c
#include "tinker_precision.h"
      module neigh
      implicit none
      integer ineigup
      integer ncell_tot,max_cell_len,max_numneigcell
      integer::bvx_cell=1,bvy_cell=1,bvz_cell=1
      integer::bex_cell=1,bey_cell=1,bez_cell=1
      integer::bcx_cell=1,bcy_cell=1,bcz_cell=1
      integer nx_cell,ny_cell,nz_cell
      real(t_p) defaultlbuffer,defaultlbuffer1,lbuffer,lbuf2
      real(t_p) vbuf2,cbuf2,mbuf2,torquebuf2
      real(t_p) vshortbuf2,cshortbuf2,mshortbuf2,torqueshortbuf2
      logical rebuild_lst,auto_lst,dovlst,doclst,domlst
      logical vlst_enable,mlst_enable,clst_enable
      logical vlst2_enable,mlst2_enable,clst2_enable
      integer bit_si,bit_sh

!DIR$ ATTRIBUTES ALIGN:64 :: nvlst
      integer, allocatable,target:: nvlst(:),nshortvlst(:)
!DIR$ ATTRIBUTES ALIGN:64 :: vlst
      integer, allocatable,target:: vlst(:,:),shortvlst(:,:)
      integer, allocatable,target:: vblst(:),ivblst(:),vb_tag(:)
      integer, allocatable,target:: shortvblst(:),ishortvblst(:)
     &                    ,shortvb_tag(:)

      ! Adjacency Matrix attributes
      integer, allocatable, target:: matb_lst(:)
      integer(int_ptr_kind()):: szoMatb=0
      real(8) :: buffMatb=2.0d0
      integer szMatb
      integer nbMatb, niterMatb, offsetlMb, offsetrMb

!DIR$ ATTRIBUTES ALIGN:64 :: nelst
      integer, allocatable,target:: nelst(:),nelstc(:),nshortelst(:)
     &                    ,nshortelstc(:)
!DIR$ ATTRIBUTES ALIGN:64 :: elst 
      integer, allocatable,target:: elst(:,:),shortelst(:,:)
      integer, allocatable,target:: eblst(:),ieblst(:),eb_tag(:)
      integer, allocatable,target:: shorteblst(:),ishorteblst(:)
     &                    ,shorteb_tag(:)
!DIR$ ATTRIBUTES ALIGN:64 :: ineignl

      integer, allocatable :: ineignl(:)
!DIR$ ATTRIBUTES ALIGN:64 :: neigcell,numneigcell,repartcell
      integer, allocatable :: neigcell(:,:),numneigcell(:),repartcell(:)
!DIR$ ATTRIBUTES ALIGN:64 :: cell_len,indcell,bufbegcell
      integer, allocatable :: cell_len(:),cell_lenr(:),cell_len2(:)
     &                       ,cell_len2r(:)
      integer, allocatable :: indcell(:),bufbegcell(:)
      integer, allocatable :: cell_scan(:),cell_scan1(:)
      integer(1), allocatable :: cell_len1(:)

!DIR$ ATTRIBUTES ALIGN:64:: xbegcell,ybegcell,zbegcell
      real(t_p), allocatable :: xbegcell(:),ybegcell(:),zbegcell(:)
!DIR$ ATTRIBUTES ALIGN:64:: xendcell,yendcell,zendcell
      real(t_p), allocatable :: xendcell(:),yendcell(:),zendcell(:)

      integer,target,allocatable :: cellv_key(:),cellv_glob(:)
     &              ,cellv_loc(:),cellv_jvdw(:)
      integer,target,allocatable :: celle_key(:),celle_glob(:)
     &              ,celle_pole(:),celle_loc(:),celle_ploc(:)
     &              ,celle_plocnl(:)
      integer,target,allocatable :: celle_chg(:)
      real(t_p),target,allocatable:: celle_x(:),celle_y(:),celle_z(:)

      parameter( bit_si=8*sizeof(bvx_cell) )
      parameter( bit_sh=5 )
      parameter(
#if (defined(SINGLE) || defined(MIXED))
     &  defaultlbuffer  = 0.7
     & ,defaultlbuffer1 = 1.0
#else
     &  defaultlbuffer  = 2.0
     & ,defaultlbuffer1 = 2.0
#endif
     & )

!$acc declare create(nvlst,nelst,nelstc)
      end
