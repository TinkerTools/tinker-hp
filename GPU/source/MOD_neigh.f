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
c     boxPart contains box partitionning information
c     pair_atlst  a derived type that gather all attributes of a pairwise
c                neighbor list stored in a list structure 
c
c
#include "tinker_macro.h"
      module neigh
      implicit none
      type boxPart
         integer ipart,bx_c,by_c,bz_c,nx_c,ny_c,nz_c
         logical nx_l,ny_l,nz_l
      end type
      type pair_atlst
         integer   natmnl,natmnlb,nbpairs,npairs
         ! Spatial reordering key and global atoms
         integer  ,allocatable:: c_key(:),c_glob(:)
         ! pairwise block-list and atom-list
         integer  ,allocatable:: blist(:),list(:)
         ! (cutoff+lbuffer)**2
         real(t_p) cut_buff2
         ! Spatial reordering of positions
         real(t_p),allocatable:: cell_x(:),cell_y(:),cell_z(:)
         type(boxPart) bPar
      end type

      ! Type Element nomenclature
      ! a   : atom
      ! b   : block
      ! ct  : cut
      ! l   : list
      ! n   : number
      ! o   : order
      ! p   : pair
      ! gl  : global
      ! lc  : local
      ! s   : spatial
      ! _1  : short
      ! *n  : ***** n times  (b2 : block-block)
      ! ---------------------------------------
      ! Some Details
      ! ctoff   :  cutoff
      ! ctoff2  :  cutoff**2
      ! ctbuf   :  cutoff + buffer
      ! ctbuf2  :  (cutoff + buffer)**2
      ! s_key   :  spatial key issue from radix sort
      ! b_rmid  :  block satellite data (radius & center)
      type nblst_t
         integer idx,na,nab,nb,n_tot
         integer nb2p_0,nb2p,nbap,nb2p_1,nbap_1
         integer BLOCK_SIZE
         logical use2lists
         real(t_p)  ctoff, ctoff2, ctbuf, ctbuf2
         real(t_p)  ctoff_beg, ctoff2_beg, ctbuf_beg, ctbuf2_beg
         real(t_p)  ctoff_1, ctoff2_1, ctbuf_1, ctbuf2_1
         integer,pointer    :: a_kind(:), ag_id(:)
         integer,allocatable:: na2l(:), na2l_1(:), a2l(:,:),a2l_1(:,:)
         integer,allocatable:: b2pl(:),  bapl  (:), b_stat(:)
         integer,pointer    :: b2pl_1(:),bapl_1(:)
         integer,allocatable:: s_key(:), s_kind(:), sgl_id(:),slc_id(:)
     &          , gl_key(:)
         real(t_p),allocatable:: b_rmid(:), so_x(:),so_y(:),so_z(:)
         type(boxPart) bPar
      end type

      ! b_stat values
      enum,bind(C)
         enumerator normal_block,inside_block,disc_block
      end enum
      ! nblist ID
      enum,bind(C)
         enumerator vdwl_id,displ_id
         enumerator elecl_id,chrl_id
      end enum

      integer ineigup
      integer ncell_tot,max_cell_len,max_numneigcell
      integer nx_cell,ny_cell,nz_cell
      real(t_p) defaultlbuffer,defaultlbuffer1,lbuffer,lbuf2
      real(t_p) vbuf2,cbuf2,mbuf2,torquebuf2,dbuf2
      real(t_p) vshortbuf2,cshortbuf2,mshortbuf2,torqueshortbuf2
     &        , dshortbuf2
      logical rebuild_lst,auto_lst,dovlst,doclst,domlst
      logical vlst_enable,mlst_enable,clst_enable
      logical vlst2_enable,mlst2_enable,clst2_enable
      logical only_bsat
      integer bit_si,bit_sh
      type(boxPart) e_bP,v_bP,c_bP

      ! Structure list for ani feature
      type(pair_atlst),target:: list_ani
!DIR$ ATTRIBUTES ALIGN:64 :: nvlst
      integer, allocatable,target:: nvlst(:),nshortvlst(:)
!DIR$ ATTRIBUTES ALIGN:64 :: vlst
      integer, allocatable,target:: vlst(:,:),shortvlst(:,:)
      integer, allocatable,target:: vblst(:),ivblst(:),vb_tag(:)
      integer, allocatable,target:: shortvblst(:),ishortvblst(:)
     &                    ,shortvb_tag(:)

      ! Adjacency Matrix & List attributes
      integer, allocatable,target:: matb_lst(:), bb_lst(:)
      integer(8):: szoMatb=0
      real(8)   :: buffMatb=2.0d0
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

      integer,target,allocatable :: cellv_jvdw(:)
      integer,target,allocatable :: celle_key(:),celle_glob(:)
     &              ,celle_pole(:),celle_loc(:),celle_ploc(:)
     &              ,celle_plocnl(:)
      integer,target,allocatable :: celle_chg(:)
      real(t_p),target,allocatable:: celle_x(:),celle_y(:),celle_z(:)

      integer,allocatable,target:: ngh_works_i(:)

      ! ---------------------------------------------------------
      !    nblst_t struct attibutes  ( to be use with load_nbl )
      ! ---------------------------------------------------------
      integer   na,nab,nb
      integer   nb2p,nb2p_0,nbap,nb2p_1,nbap_1
      integer  ,pointer:: na2l(:),a2l(:,:), na2l_1(:),a2l_1(:,:)
      integer  ,pointer:: b2pl  (:), bapl  (:), abpl  (:), b_stat(:)
      integer  ,pointer:: b2pl_1(:), bapl_1(:), abpl_1(:)
      integer  ,pointer:: s_key(:),s_kind(:),sgl_id(:),slc_id(:)
      real(t_p),pointer:: b_rmid(:),so_x(:),so_y(:),so_z(:)
      type(boxPart)    :: bPar
      ! =========================================================

      type(nblst_t),target:: disp_nbl, vdw_nbl, elec_nbl

      parameter( bit_si=8*sizeof(ineigup) )
      parameter( bit_sh=5 )

      interface
        subroutine build_cell_list2(ineignl,cell_order,
     &             nlocnl,buf,bP)
        import boxPart
        integer  ,intent(in)   :: nlocnl
        integer  ,intent(in)   :: ineignl(:)
        integer  ,intent(inout):: cell_order(:)
        real(t_p),intent(in)   :: buf
        type(boxPart),intent(inout)::bP
        end subroutine
      end interface

      interface
        subroutine build_pairwise_list(atList,natmnl,cut2,a_plist)
          import pair_atlst
          integer  ,intent(in)::natmnl
          integer  ,intent(in)::atList(:)
          real(t_p),intent(in):: cut2
          type(pair_atlst),target:: a_plist
        end subroutine
      end interface

      interface
        subroutine nblist_build1(l)
          import nblst_t
          type(nblst_t),intent(inout),target :: l
        end subroutine
      end interface

!$acc declare create(nvlst,nelst,nelstc)
      contains

      subroutine init_aboxPart(a_bP)
      implicit none
      type(boxPart),intent(inout):: a_bP
      a_bP%ipart=0
      a_bP%bx_c=1
      a_bP%by_c=1
      a_bP%bz_c=1
      a_bP%nx_l=.true.
      a_bP%ny_l=.true.
      a_bP%nz_l=.true.
      end subroutine

      subroutine init_boxPart
      implicit none
      e_bP%ipart=0
      e_bP%bx_c=1
      e_bP%by_c=1
      e_bP%bz_c=1
      e_bP%nx_l=.true.
      e_bP%ny_l=.true.
      e_bP%nz_l=.true.
      v_bP%ipart=0
      v_bP%bx_c=1
      v_bP%by_c=1
      v_bP%bz_c=1
      v_bP%nx_l=.true.
      v_bP%ny_l=.true.
      v_bP%nz_l=.true.
      c_bP%ipart=0
      c_bP%bx_c=1
      c_bP%by_c=1
      c_bP%bz_c=1
      c_bP%nx_l=.true.
      c_bP%ny_l=.true.
      c_bP%nz_l=.true.
      call init_aboxPart(list_ani%bPar)
      call init_aboxPart(vdw_nbl%bPar)
      call init_aboxPart(elec_nbl%bPar)
      call init_aboxPart(disp_nbl%bPar)
      end subroutine

      subroutine init_nbl(l,id,typ,glob,na_,ntot,b_siz
     &                   ,lcut,cut,cut_1,lbuf,use2lists_)
      implicit none
      type(nblst_t)              :: l
      integer  ,target,intent(in):: typ(:),glob(:)
      integer  ,intent(in) :: na_,ntot,b_siz,id
      logical  ,intent(in) :: use2lists_
      real(t_p),intent(in) :: cut,cut_1,lcut,lbuf
      integer siz_t,siz_g

      siz_t = size(typ)
      siz_g = size(glob)

      l%a_kind (1:siz_t) => typ (1:)
      l%ag_id  (1:siz_g) => glob(1:)
      l%idx        = id
      l%BLOCK_SIZE = b_siz
      l%n_tot      = ntot
      l%na         = na_
      l%nb         = (na_-1)/b_siz+1
      l%nab        = b_siz*l%nb
      l%nab        = merge(na_+b_siz,l%nab,na_.eq.l%nab)
      na           = l%na
      nb           = l%nb
      nab          = l%nab

      l%ctoff      = cut
      l%ctoff2     = l%ctoff*l%ctoff
      l%ctbuf      = cut+lbuf
      l%ctbuf2     = l%ctbuf*l%ctbuf
      l%ctoff_beg  = lcut
      l%ctoff2_beg = l%ctoff_beg*l%ctoff_beg
      l%ctbuf_beg  = lcut-lbuf
      l%ctbuf2_beg = l%ctbuf_beg*l%ctbuf_beg
      l%ctoff_1    = cut_1
      l%ctoff2_1   = l%ctoff_1*l%ctoff_1
      l%ctbuf_1    = cut_1+lbuf
      l%ctbuf2_1   = l%ctbuf_1*l%ctbuf_1
      l%use2lists  = use2lists_
      end subroutine

      subroutine load_nbl2mod(l,option)
      implicit none
      type(nblst_t),target:: l
      integer,intent(in),optional:: option
      integer,parameter:: deflt=0,short=1
      integer(8) siz1,siz2
      integer n_nbl

      na     = l%na
      nb     = l%nb
      nab    = l%nab
      nb2p_0 = l%nb2p_0
      nb2p   = l%nb2p
      nbap   = l%nbap
      n_nbl  = merge(2,1,l%use2lists)

      if (allocated(l%na2l)) then
         siz1 = size(l%na2l)
         na2l(1:siz1) => l%na2l(1:)
         siz1 = size(l%a2l,1)
         siz2 = size(l%a2l,2)
         a2l(1:siz1,1:siz2) => l%a2l(1:siz1*siz2,1)
      end if
      if (allocated(l%b2pl)) then
         siz1 = size(l%b2pl)
         b2pl(1:siz1)     => l%b2pl(1:siz1)
         bapl(1:l%nb2p_0) => l%bapl(1:)
         abpl(1:l%nb2p_0*l%BLOCK_SIZE) => l%bapl(n_nbl*l%nb2p_0+1:)
      end if
      if (allocated(l%b_stat)) then
         b_stat(1:  nb) => l%b_stat(1:)
         b_rmid(1:4*nb) => l%b_rmid(1:)
      end if
      if (allocated(l%s_key)) then
         siz1 = size(l%s_key)
         s_key(1:siz1) => l%s_key(1:siz1)
         sgl_id(1:nab) => l%sgl_id(1:)
         slc_id(1:nab) => l%slc_id(1:)
         so_x(1:nab) => l%so_x(1:)
         so_y(1:nab) => l%so_y(1:)
         so_z(1:nab) => l%so_z(1:)
      end if
      if (allocated(l%s_kind)) then
         siz1 = size(l%s_kind)
         s_kind(1:siz1) => l%s_kind(1:siz1)
      end if

      if (l%use2lists) then
         nbap_1 = l%nbap_1
         nb2p_1 = l%nb2p_1
         if (allocated(l%a2l_1)) then
            siz1 = size(l%na2l_1)
            na2l_1(1:siz1) => l%na2l_1(1:siz1)
            siz1 = size(l%a2l_1,1)
            siz2 = size(l%a2l_1,2)
            a2l_1(1:siz1,1:siz2) => l%a2l_1(1:siz1*siz2,1)
         end if
         if (associated(l%b2pl_1)) then
            siz1 = size(l%b2pl_1)
            b2pl_1(1:nb2p_1) => l%b2pl(nb2p_0+1:)
            siz1 = size(l%bapl_1)
            bapl_1(1:nbap_1) => l%bapl(nb2p_0+1:)
            abpl_1(1:nb2p_0*l%BLOCK_SIZE) => l%bapl
     &                     ((2+l%BLOCK_SIZE)*l%nb2p_0+1:)
         end if
      end if
      end subroutine

      end
