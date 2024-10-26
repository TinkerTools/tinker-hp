c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #######################################################################
c     ##                                                                   ##
c     ##  subroutine nblist_vdw  -- set of vdw neighbor lists subroutines  ##
c     ##                                                                   ##
c     #######################################################################
c
c
c     "nblist" constructs and maintains nonbonded pair neighbor lists
c     for vdw
c
c
#include "tinker_macro.h"
      module nblist_b_inl
        contains
#include "image.f.inc"
#include "midpointimage.f.inc"
      end module

c
c     spatial reordering of vdw data
c
      subroutine set_VdwData_CellOrder
      use atoms  ,only: n
      use atmlst ,only: vdwglobnl
      use domdec
      use kvdws  ,only: rad,eps,radv,epsv
      use neigh
      use utilgpu
      use vdw
      implicit none
      integer i,it,iglob,nab_,na_

      nab_ = vdw_nbl%nab
      na_  = vdw_nbl%na
      call prmem_request(slc_id    ,nab_,async=.true.)
      call prmem_request(cellv_jvdw,nab_,async=.true.)
      call prmem_request(radv      ,nab_,async=.true.)
      call prmem_request(epsv      ,nab_,async=.true.)

!$acc parallel loop async(def_queue)
!$acc&         present(sgl_id,slc_id,cellv_jvdw,
!$acc&   s_key,loc,jvdw_c,jvdw,vdwglobnl,rad,eps,radv,epsv)
      do i = 1, nab_
         if (i.le.na_) then
            iglob         = ivdw(vdwglobnl(s_key(i)))
            sgl_id(i)     = iglob
            cellv_jvdw(i) = jvdw_c(iglob)
            radv(i)       = rad(jvdw(iglob))
            epsv(i)       = eps(jvdw(iglob))
         else
            sgl_id(i)     = n
            cellv_jvdw(i) = 1
            radv(i)       = inf
            epsv(i)       = inf
         end if
      end do
      end subroutine

      subroutine spatialOrder_pos
      use atoms   ,only: n,x,y,z
      use domdec
      use neigh
      use nblist_b_inl
      use utilgpu
      use vdw
      implicit none
      integer i,iglob
      real(t_p) xr,yr,zr

!$acc parallel loop async(def_queue) default(present)
      do i = 1, nab
         if (i.le.na) then
            iglob     = sgl_id(i)
            slc_id(i) = loc(iglob)
            xr        = x(iglob)
            yr        = y(iglob)
            zr        = z(iglob)
            call image_inl(xr,yr,zr)
            so_x(i)   = xr
            so_y(i)   = yr
            so_z(i)   = zr
         else
            slc_id(i) = nbloc
            so_x(i)   = inf
            so_y(i)   = inf
            so_z(i)   = inf
         end if
      end do
      end subroutine
c
c     build vdw block-block &  block-atoms neighbor list
c
      subroutine vdw_nblist_build1
      use atmlst  ,only: vdwglobnl
      use cutoff  ,only: use_shortvlist,vdwcut,vdwshortcut,shortheal
      use inform  ,only: deb_Path
      use neigh   ,only: vdw_nbl,vdwl_id,lbuffer,init_nbl,nblist_build1
      use utilgpu ,only: BLOCK_SIZE
      use vdw     ,only: nvdwlocnl,ired,kred
     &            ,nvdwlocnlb_pair,nvdwlocnlb_pair1,nvdwlocnlb2_pair
     &            ,nshortvdwlocnlb2_pair,nvdw,ivdw

      if (deb_Path) write(*,*) 'vdw_nblist_build1'

      call init_nbl(vdw_nbl,vdwl_id,vdwglobnl,ivdw,nvdwlocnl,nvdw
     &             ,BLOCK_SIZE,vdwshortcut-shortheal,vdwcut
     &             ,vdwshortcut,lbuffer,use_shortvlist)
      call nblist_build1(vdw_nbl)
      call set_VdwData_CellOrder
      end subroutine

      subroutine nblist_build1(l)
      use atoms   ,only: x,y,z,n
      use boxes   ,only: octahedron
      use cell    ,only: xcell2,ycell2,zcell2,xcell,ycell,zcell,eps_cell
      use domdec  ,only: rank,nproc
      use neigh   ,only: nblst_t,nb2p_0,na,nab,nb
     &            ,nbMatb,niterMatb,buffMatb,offsetlMb,offsetrMb,szMatb
     &            ,matb_lst,bb_lst
     &            ,cell_scan,cell_len,cell_lenr,cell_len2,cell_len2r 
     &            ,s_key,sgl_id,slc_id
     &            ,bit_si,bit_sh,disc_block,inside_block,vshortbuf2
     &            ,lbuffer,build_cell_list2
     &            ,v_bP,b2pl,b2pl_1,bapl,bapl_1,abpl,abpl_1
     &            ,b_stat,b_rmid,trackb=>ngh_workS_i,cellv_jvdw
     &            ,so_x,so_y,so_z
      use utils   ,only: set_to_zero1_int1,set_to_zero1_int
     &            ,associate_ptr
      use inform  ,only: deb_Path, app_id,pimd_a,analyze_beads_a,pibar_a
      use nblist_b_inl
      use interfaces ,only: pre_process_adjacency_matrix

#ifdef _OPENACC
      use thrust  ,only: thrust_exclusive_scan
      use nblistcu,only: sztb,filter_lst_sparse,filter_lsts_sparse
     &            ,blocksat_kcu,InMID_BS
      use utilcu  ,only: check_launch_kernel,mem_download_i
     &            ,prdmem_request,dmem_set,chk_cuAPI_O
      use cudafor
#endif

      use tinheader,only: ti_p
      use tinMemory,only: prmem_mvrequest,prmem_request,szoi,mipk
      use utilgpu ,only: rec_queue,BLOCK_SIZE,inf,mem_set
     &            ,rec_stream,mem_move_i4,get_atom_glob_id
      implicit none
      type(nblst_t),intent(inout),target :: l
      integer i,iblock,inblock,icell
      integer k,kblock,knblock,kcell
      integer iglob,iv,ierrSync
      integer it1,it2,block1,block2,temp,bit2,pos2
      integer nx,ny,nz
      integer nxy_cell,icell_scan,kcell_scan
      integer posi,nGS,n_nbl
      integer iter ,siztr,sizba,sizb2,nb2p_old
      integer nkblock,complete_b,nb_pair,nneig
      integer(1) icell_len,kcell_len
      integer(mipk) i8_
      integer,pointer:: atomk(:),atomgid(:)
      real(t_p) l_xcell,l_ycell,l_zcell,rdn,rdn1
      real(t_p) ctbuf,ctbuf2,ctbuf2_1,ctbuf2_b

c     real(t_p) xi,yi,zi
c     integer j,ibl,iblv(32)

      if (deb_Path) write(*,*) ' nblist_build1'

      ! Get parameters and set variables from structure
      ctbuf    = l%ctbuf
      ctbuf2   = l%ctbuf2
      ctbuf2_1 = l%ctbuf2_1
      ctbuf2_b = l%ctbuf2_beg
      l%nb2p_0 = 0
      l%nb2p   = 0
      l%nb2p_1 = 0
      l%nbap   = 0
      l%nbap_1 = 0
      nb2p_0   = 0
      ! Get data from struct via associate_ptr
      call associate_ptr(atomk,l%a_kind,int(l%na,mipk))
      call associate_ptr(atomgid,l%ag_id,int(n,mipk))

      ! Add margin to serve as out-of-bound exclusion interactions for C2 nblist
      szMatb    = ((nb-1)/bit_si) + 1
      nbMatb    = max(1,min(int(buffMatb*1d9/(4*szMatb)),nb))
      niterMatb = (nb-1)/nbMatb + 1

      n_nbl     = merge(2,1,l%use2lists)
#ifdef _OPENACC
      siztr     = n_nbl*(sztb+1)*nb+32
#endif

      ! Allocate buffers
      call prmem_request(l%so_x,   nab,async=.false.)
      call prmem_request(l%so_y,   nab,async=.false.)
      call prmem_request(l%so_z,   nab,async=.false.)
      call prmem_request(l%sgl_id, nab,async=.false.)
      call prmem_request(l%slc_id, nab,async=.false.)
      call prmem_request(l%s_key,    n,async=.false.)
      call prmem_request(l%b_stat,  nb,async=.false.)
      call prmem_request(l%b_rmid,4*nb,async=.false.)
      call prmem_request(trackb, siztr,async=.false.)

      ! Attach module pointers to structure
      b_stat(1:nb)   => l%b_stat(1:)
      b_rmid(1:nb*4) => l%b_rmid(1:)
      sgl_id(1:nab)  => l%sgl_id(1:)
      slc_id(1:nab)  => l%slc_id(1:)
      so_x(1:nab)    => l%so_x  (1:)
      so_y(1:nab)    => l%so_y  (1:)
      so_z(1:nab)    => l%so_z  (1:)
      s_key(1:n)     => l%s_key (1:)

      ! Box partitionning and spatial reodering
      if (l%n_tot.eq.n) then
         call build_cell_list2(l%a_kind,s_key,na,l%ctbuf,l%bPar)
      else
         call get_atom_glob_id(l%a_kind,l%ag_id,sgl_id,na,n)
         call build_cell_list2(sgl_id,s_key,na,l%ctbuf,l%bPar)
      end if
c
c     Build Adjacency list
c
      do iter = 0,niterMatb-1
        offsetlMb = iter*nbMatb
        offsetrMb = min((iter+1)*nbMatb,nb)
        if (iter.eq.niterMatb-1) nbMatb=nb-offsetlMb

        if (octahedron) then
           call build_adjacency_matrix_octahedron(na,l%bPar)
        else
           call build_adjacency_matrix(na,l%bPar)
        end if

        call pre_process_adjacency_matrix_ (nb,nb_pair)
        nb2p_old = 2*nb2p_0 + 1
        nb2p_0 = nb2p_0 + nb_pair

        if (niterMatb.eq.1) then
           call prmem_request (bb_lst,2*nb2p_0,async=.true.)
        else
           call prmem_mvrequest (bb_lst,2*nb2p_0,async=.true.)
        end if

        call fill_adjacency_list_(nb,nb_pair,bb_lst(nb2p_old:))
      end do
      l%nb2p_0 = nb2p_0  ! Set nb2p_0 to struct
 
      ! Allocate memory for lists variables & Attach associated pointers
      sizba = n_nbl*nb2p_0*(BLOCK_SIZE+1)
      sizb2 = n_nbl*nb2p_0
      call prmem_request (l%bapl,sizba,async=.false.)
      call prmem_request (l%b2pl,sizb2,async=.false.)
      bapl(1:sizba)             => l%bapl(1:) 
      abpl(1:nb2p_0*BLOCK_SIZE) => l%bapl(n_nbl*nb2p_0+1:)
      b2pl(1:sizb2)             => l%b2pl(1:) 
      if (l%use2lists) then
         l%bapl_1(1:nb2p_0) => l%bapl(nb2p_0+1:)
           bapl_1(1:nb2p_0) => l%bapl(nb2p_0+1:) 
           abpl_1(1:nb2p_0*BLOCK_SIZE) => l%bapl
     &                         (n_nbl*nb2p_0+nb2p_0*BLOCK_SIZE+1:)
         l%b2pl_1(1:nb2p_0) => l%b2pl(nb2p_0+1:)
           b2pl_1(1:nb2p_0) => l%b2pl(nb2p_0+1:) 
      end if

      ! Zero list & workspace
      i8_ = n_nbl*nb2p_0*int(BLOCK_SIZE+1,mipk)
      if (i8_.gt.huge(sizba)) then
 17      format('ERROR',/
     &         ,'Required size for nblist is above 32 bits numerical '
     &         ,'limits. ','Crash Idetified or system is to largea'
     &         ,2I16)
         write(0,17) i8_,sizba
         call fatal
      end if
      call mem_set(bapl,0,i8_,rec_stream)
      call mem_set(trackb,0,int(siztr,mipk),rec_stream)

      ! Debug info
      if (deb_Path) then
 13   format( I3," iteration(s) required to build Adj lst" )
 14   format( I10," nblocks process at once over",I10 )
 15   format(1x,'Size Adj Mat (Mo)',F12.3,';',4x,'buffer (Go)',F10.6
     &    ,/,' Cutbuf2',F9.4,' shortcut',F9.4,' nbl buffer', F9.4 )
 16   format(1x,'natoms',I8,';',3x,'vlist mem space (Mo)',F12.3)
         nbMatb = max(1,min(int(buffMatb*1d9/(4*szMatb)),nb))
         i8_ = nb2p_0*(BLOCK_SIZE+4)*szoi
         if (l%use2lists) i8_ = 2*i8_
         print 13, niterMatb
         print 14, nbMatb,nb
         print 15, int(szMatb*4,8)*nbMatb/1d6, buffMatb,l%ctbuf2
     &           ,l%ctoff_beg,lbuffer
         print 16, na,i8_/1d6
      end if
c
      l_xcell = xcell/l%bPar%nx_c
      l_ycell = ycell/l%bPar%ny_c
      l_zcell = zcell/l%bPar%nz_c
      nGS     = min(2**12,max(nb,4)/4)
#ifdef _CUDA
!$acc host_data use_device(sgl_id,cell_scan,cell_len
!$acc&    ,x,y,z,matb_lst,bb_lst,atomk,atomgid,s_key
!$acc&    ,sgl_id,slc_id,cellv_jvdw,trackb,b2pl,abpl,bapl
!$acc&    ,so_x,so_y,so_z,b_stat,b_rmid)

      ! Build block satellite data ( middle point and locality )
      call blocksat_kcu<<<nGS,128,0,rec_stream>>>
     &     (atomk,atomgid,s_key,x,y,z
     &     ,nb,na,nab,n,ctbuf
     &     ,l_xcell,l_ycell,l_zcell,eps_cell
     &     ,slc_id,sgl_id,so_x,so_y,so_z,b_stat,b_rmid)
      call check_launch_kernel(" blocksat_kcu")

      ! Set irregular flag to all blocks if box is to small
      !TODO Lambda : Set to true
      if(app_id==pimd_a .or. app_id==analyze_beads_a
     &     .or. app_id==pibar_a.or.octahedron.or.l%use2lists) then
!$acc parallel loop async(rec_queue) deviceptr(b_stat)
         do i = 1,nb
            b_stat(i)=disc_block
         end do
      else if (3*ctbuf.gt.xcell) then
      !  call dmem_set(b_stat,disc_block,int(nb,8),rec_stream)
!$acc parallel loop async(rec_queue) deviceptr(b_stat)
         do i = 1,nb
            if (b_stat(i).ne.inside_block) b_stat(i)=disc_block
         end do
      end if
c
c     From Block-Block to compressed Block-Atoms list
c
      if (l%use2lists) then
!$acc host_data use_device(abpl_1)
        call filter_lsts_sparse <<<*,128,0,rec_stream>>>
     &       (sgl_id,cell_scan,so_x,so_y,so_z,b_stat,b_rmid,matb_lst
     &       ,nab,nb,szMatb,nb2p_0
     &       ,ctbuf2_1,ctbuf2_b,ctbuf2
     &       ,bb_lst,bapl,abpl,trackb(1)
     &       ,b2pl,trackb(3),abpl_1
     &       ,trackb(17),trackb(n_nbl*nb+17))
        call check_launch_kernel(" filter_lsts_sparse")
!$acc end host_data
      else

        call filter_lst_sparse <<<*,128,0,rec_stream>>>
     &       (sgl_id,cell_scan,so_x,so_y,so_z,b_stat,b_rmid,matb_lst
     &       ,nab,nb,szMatb,nb2p_0,ctbuf2
     &       ,bb_lst,bapl,abpl,trackb(1),b2pl,trackb(3)
     &       ,trackb(17),trackb(nb+17))
        call check_launch_kernel(" filter_lst_sparse")

      end if

!$acc end host_data
!$acc update host(trackb(1:16)) async(rec_queue)
#else
      !TODO Add an OpenACC version of this kernel
#endif
      if (deb_Path) then
         it1=0; it2=0
!$acc parallel loop async
         do i = 1,nb;
            if (b_stat(i).eq.2) it1=it1+1;
            if (b_stat(i).eq.1) it2=it2+1;
         end do
!$acc wait(rec_queue)
         print '(2A,3I6)', 'number of disc/irregular & '
     &        ,'internal & normal blocks'
     &        ,it1,it2,nb-(it1+it2)
!$acc wait
         call get_vdw_nblst_sizes
         print '(A,6I12)','C1 vs C2 vdw nblist pairs blocks number '
     &         , nb2p_0  ,trackb(1:4)
      end if
c!$acc update host(b2pl,bapl,abpl)
c      do i = 0,trackb(3)-1
c         print*, i+1,b2pl(2*i+1),b2pl(2*i+2)
c      end do
c      do i = 0,trackb(1)-1
c         ibl = bapl(i+1)
c         if (ibl.eq.13.or.ibl.eq.14) then
c         iblv(:) = [(abpl(32*i+j),j=1,32)] 
c 43      format( I6,32I5 )
c         print 43, ibl, iblv
c         end if
c      end do
      end subroutine
c
      subroutine get_vdw_nblst_sizes
      use neigh
      use vdw
      implicit none

      nvdwlocnlb       = vdw_nbl%nab
      nvdwlocnlb2_pair = ngh_workS_i(1)
      nvdwlocnlb_pair1 = ngh_workS_i(3)
      vdw_nbl%nbap     = ngh_workS_i(1)
      vdw_nbl%nb2p     = ngh_workS_i(3)

      if (vdw_nbl%use2lists) then
         vdw_nbl%nbap_1 = ngh_workS_i(2)
         vdw_nbl%nb2p_1 = ngh_workS_i(4)
      else
         vdw_nbl%nb2p_1 = 0
         vdw_nbl%nbap_1 = 0
      end if
      end subroutine
c
      subroutine nblst_alter_bstat
      use neigh   ,only: vdw_nbl,disp_nbl,elec_nbl,d_bl=>disc_block
      use utilgpu ,only: mem_set,rec_stream
      implicit none
      integer i
      integer(ipk_) nb_

      if (allocated(vdw_nbl%b_stat)) then
         nb_ = vdw_nbl%nb
         call mem_set(vdw_nbl%b_stat,d_bl,nb_,rec_stream)
      end if
      if (allocated(disp_nbl%b_stat)) then
         nb_ = disp_nbl%nb
         call mem_set(disp_nbl%b_stat,d_bl,nb_,rec_stream)
      end if
      if (allocated(elec_nbl%b_stat)) then
         nb_ = elec_nbl%nb
         call mem_set(elec_nbl%b_stat,d_bl,nb_,rec_stream)
      end if
      end subroutine

      subroutine mlist_block
      use atoms   ,only: x,y,z,n
      use atmlst  ,only: poleglobnl
      use boxes   ,only: octahedron
      use cutoff  ,only: use_shortmlist
      use cell    ,only: xcell2,ycell2,zcell2,xcell,ycell,zcell
      use domdec  ,only: rank,nproc
      use mpole   ,only: npolelocnl,npolelocnlb
     &            ,npolelocnlb_pair,npolelocnlb2_pair
     &            ,nshortpolelocnlb2_pair,ipole,npole
      use neigh   ,only: cell_scan,cell_len,matb_lst,celle_glob,eblst
     &            ,mbuf2,mshortbuf2
     &            ,nbMatb,niterMatb,buffMatb,offsetlMb,offsetrMb,szMatb
     &            ,celle_key,celle_x,celle_y,celle_z,cell_lenr
     &            ,ieblst,shorteblst,ishorteblst,bit_si,bit_sh
     &            ,cell_len2,cell_len2r,build_cell_list2,e_bP
      use utils   ,only: set_to_zero1_int1,set_to_zero1_int
      use inform  ,only: deb_Path
      use interfaces,only:pre_process_adjacency_matrix
     &               ,set_ElecData_CellOrder
#ifdef _OPENACC
     &               ,cu_filter_lst_sparse
      use thrust  ,only: thrust_exclusive_scan
      use nblistcu,only: filter_lst_sparse
      use utilcu  ,only: check_launch_kernel
      use utilgpu ,only: rec_stream
      use cudafor
#endif
      use utilgpu ,only: rec_queue,BLOCK_SIZE,inf
      use tinheader,only: ti_p
      use tinMemory,only: prmem_request,prmem_mvrequest,mipk,szoi
      implicit none
      integer i,iblock,inblock,icell
      integer k,kblock,knblock,kcell
      integer iglob,iv,ierrSync
      integer it1,it2,block1,block2,temp,nblock,bit2,pos2
      integer nx,ny,nz
      integer nxy_cell,icell_scan,kcell_scan
      integer posi
      integer iter
      integer nkblock,complete_b,nb_pair,nneig
      integer(1) icell_len,kcell_len
      integer(mipk) szolst
      real(t_p) lenx_cell,leny_cell,lenz_cell,rdn,rdn1

      if (deb_Path) print*, 'mlist_neiglist_block'
      npolelocnlb_pair = 0
      npolelocnlb = BLOCK_SIZE*((npolelocnl-1)/BLOCK_SIZE + 1)
      ! Add margin to serve as out-of-bound exclusion interactions for C2 nblist
      if (npolelocnlb == npolelocnl)
     &   npolelocnlb = npolelocnl + BLOCK_SIZE

      nblock     = ((npolelocnl-1)/BLOCK_SIZE) + 1
      szMatb     = ((nblock-1)/bit_si) + 1

      nbMatb     = max(1,min(int(buffMatb*1d9/(4*szMatb)),nblock))
      niterMatb  = (nblock-1)/nbMatb + 1
      npolelocnlb_pair = 0

      call prmem_request(celle_key,n,async=.true.)

      if (npole.eq.n) then
         call build_cell_list2(poleglobnl,celle_key,npolelocnl,
     &               sqrt(mbuf2),e_bP)
      else
         call prmem_request(celle_glob,npolelocnlb,async=.true.)
!$acc parallel loop async default(present)
         do i = 1,npolelocnl
            celle_glob(i) = ipole(poleglobnl(i))
         end do
         call build_cell_list2(celle_glob,celle_key,npolelocnl,
     &               sqrt(mbuf2),e_bP)
      end if
c
c     Build Adjacency list piece by piece
c
      do iter = 0,niterMatb-1

        offsetlMb = iter*nbMatb
        offsetrMb = min((iter+1)*nbMatb,nblock)
        if (iter.eq.niterMatb-1) nbMatb=nblock-offsetlMb

        if (octahedron) then
           call build_adjacency_matrix_octahedron(npolelocnl,e_bP)
        else
           call build_adjacency_matrix(npolelocnl,e_bP)
        end if

        call pre_process_adjacency_matrix (nblock,nb_pair)
        npolelocnlb_pair = npolelocnlb_pair + nb_pair

        call prmem_mvrequest (eblst,npolelocnlb_pair*(BLOCK_SIZE+2),
     &          async=.true.)
        call fill_adjacency_list(nblock,npolelocnlb_pair,eblst)

      end do

      call prmem_request (ieblst,npolelocnlb_pair,async=.true.)

      ! zero list and ilist
      call set_to_zero1_int(eblst(npolelocnlb_pair*2+1),
     &          npolelocnlb_pair*(BLOCK_SIZE),rec_queue)
      call set_to_zero1_int(ieblst,npolelocnlb_pair,rec_queue)

      if (use_shortmlist) then
         call prmem_request (shorteblst,npolelocnlb_pair*
     &                         (BLOCK_SIZE+2),async=.true.)
         call prmem_request (ishorteblst,npolelocnlb_pair,
     &                         async=.true.)
         call set_to_zero1_int(shorteblst,npolelocnlb_pair*
     &                         (2+BLOCK_SIZE),rec_queue)
         call set_to_zero1_int(ishorteblst,npolelocnlb_pair,
     &                         rec_queue)
      end if

      call Init_blocklen(nblock)
      call set_ElecData_CellOrder(.true.)

c
c     Debug print
c
      if (deb_Path) then
c13   format( I3," iteration(s) required to build Adj lst" )
c14   format( I10," nblocks process at once over",I10 )
c15   format(1x,'Size Adj Mat (Mo)',F12.3,';',4x,'buffer (Go)',F10.6)
 16   format(1x,'npolelocnl',I8,';',3x,'mlist mem space (Mo)',F12.3)
         !nbMatb = max(1,min(int(buffMatb*1d9/(4*szMatb)),nblock))
         szolst = npolelocnlb_pair*(BLOCK_SIZE+3)*szoi
         if (use_shortmlist) szolst = 2*szolst
         !print 13, niterMatb
         !print 14, nbMatb,nblock
         !print 15, int(szMatb*4,8)*nbMatb/1d6, buffMatb
         print 16, npolelocnl, szolst/1d6
      end if

c
c     From Block-Block to Block-Atoms list
c
#ifdef _CUDA
      if (use_shortmlist) then
!$acc host_data use_device(celle_glob,cell_scan,cell_len,cell_lenr,
!$acc&    celle_x,celle_y,celle_z,matb_lst,eblst,shorteblst,
!$acc&    cell_len2,cell_len2r)
      call filter_shortlst_sparse1 <<<*,128,0,rec_stream>>>
     &     (celle_glob,cell_scan,celle_x,celle_y,celle_z,matb_lst,
     &      npolelocnlb,nblock,szMatb,npolelocnlb_pair,
     &      mshortbuf2,0.0_ti_p,mbuf2,
     &      eblst,cell_len,cell_lenr,shorteblst,cell_len2,cell_len2r)
      call check_launch_kernel(" mlist:filter_shortlst_sparse1")
!$acc end host_data
      else
!$acc host_data use_device(celle_glob,cell_scan,cell_len,cell_lenr,
!$acc&    celle_x,celle_y,celle_z,matb_lst,eblst)

c     call cu_filter_lst_sparse
c    &     (celle_glob,cell_scan,celle_x,celle_y,celle_z,matb_lst,
c    &      npolelocnlb,nblock,szMatb,npolelocnlb_pair,mbuf2,
c    &      eblst,cell_len,xcell,ycell,zcell,xcell2,ycell2,zcell2,
c    &      rec_stream)
c     call filter_lst_sparse<<<*,128,0,rec_stream>>>
c    &     (celle_glob,cell_scan,celle_x,celle_y,celle_z,matb_lst,
c    &      npolelocnlb,nblock,szMatb,npolelocnlb_pair,mbuf2,
c    &      eblst,cell_len)
      call filter_lst_sparse1<<<*,128,0,rec_stream>>>
     &     (celle_glob,cell_scan,celle_x,celle_y,celle_z,matb_lst,
     &      npolelocnlb,nblock,szMatb,npolelocnlb_pair,mbuf2,
     &      eblst,cell_len,cell_lenr)
      call check_launch_kernel(" filter_lst_sparse1")

!$acc end host_data
      end if
#else
      !TODO Add an OpenACC version of this kernel
#endif

      call finalize_list_C2
     &     (nblock,npolelocnlb_pair,npolelocnlb
     &     ,npolelocnlb2_pair,ieblst,eblst,cell_len,cell_scan
     &     ,cell_lenr)

      if (use_shortmlist) call finalize_list_C2
     &        (nblock,npolelocnlb_pair,npolelocnlb
     &        ,nshortpolelocnlb2_pair,ishorteblst,shorteblst
     &        ,cell_len2,cell_scan,cell_len2r)

      end subroutine

c
c    "clist_block" performs a complete rebuild of the block
c     electrostatic neighbor lists for charges using linked cells method
c
      subroutine clist_block
      use atoms   ,only: x,y,z,n
      use atmlst  ,only: chgglobnl
      use boxes   ,only: octahedron
      use charge  ,only: nionlocnl,nionlocnlb
     &            ,nionlocnlb_pair,nionlocnlb2_pair
     &            ,nshortionlocnlb2_pair,iion,nion
      use cell    ,only: xcell2,ycell2,zcell2,xcell,ycell,zcell
      use cutoff  ,only: use_shortclist
      use domdec  ,only: rank,nproc
      use neigh   ,only: cell_scan,cell_len,matb_lst,celle_glob,eblst
     &            ,cbuf2,cshortbuf2,vbuf2,vshortbuf2
     &            ,nbMatb,niterMatb,buffMatb,offsetlMb,offsetrMb,szMatb
     &            ,celle_key,celle_x,celle_y,celle_z,cell_lenr
     &            ,ieblst,shorteblst,ishorteblst,bit_si,bit_sh
     &            ,cell_len2,cell_len2r,build_cell_list2,c_bP
      use utils   ,only: set_to_zero1_int1,set_to_zero1_int
      use inform  ,only: deb_Path
      use interfaces,only: pre_process_adjacency_matrix
     &              ,set_ElecData_CellOrder
#ifdef _OPENACC
     &              ,cu_filter_lst_sparse
      use thrust  ,only: thrust_exclusive_scan
      use nblistcu,only: filter_lst_sparse
      use utilcu  ,only: check_launch_kernel
      use utilgpu ,only: rec_stream
      use cudafor
#endif
      use potent  ,only: fuse_chglj
      use utilgpu ,only: rec_queue,BLOCK_SIZE,inf
      use tinheader,only: ti_p
      use tinMemory,only: prmem_request,prmem_mvrequest,mipk,szoi
      implicit none
      integer i,iblock,inblock,icell
      integer k,kblock,knblock,kcell
      integer iglob,iv,ierrSync
      integer it1,it2,block1,block2,temp,nblock,bit2,pos2
      integer nx,ny,nz
      integer nxy_cell,icell_scan,kcell_scan
      integer posi
      integer iter
      integer nkblock,complete_b,nb_pair,nneig
      integer(1) icell_len,kcell_len
      integer(mipk) szolst
      real(t_p) lenx_cell,leny_cell,lenz_cell,rdn,rdn1
      real(t_p) cutb2,cutb2_s
c
      if (deb_Path) print*, 'clist_neiglist_block'
      nionlocnlb_pair = 0
      nionlocnlb = BLOCK_SIZE*((nionlocnl-1)/BLOCK_SIZE + 1)
      ! Add margin to serve as out-of-bound exclusion interactions for C2 nblist
      if (nionlocnlb == nionlocnl)
     &   nionlocnlb = nionlocnl + BLOCK_SIZE

      nblock     = ((nionlocnl-1)/BLOCK_SIZE) + 1
      szMatb     = ((nblock-1)/bit_si) + 1

      nbMatb     = max(1,min(int(buffMatb*1d9/(4*szMatb)),nblock))
      niterMatb  = (nblock-1)/nbMatb + 1
      nionlocnlb_pair = 0
      cutb2      = merge(max(cbuf2,vbuf2),cbuf2,fuse_chglj)
      !cutb2_s    = merge(max(cshortbuf2,vshortbuf2),cshortbuf2,fuse_chglj)

      call prmem_request(celle_key,n,async=.true.)

      if (nion.eq.n) then
         call build_cell_list2(chgglobnl,celle_key,nionlocnl,
     &               sqrt(cutb2),c_bP)
      else
         call prmem_request(celle_glob,nionlocnlb,async=.true.)
!$acc    parallel loop async(rec_queue) default(present)
         do i = 1,nionlocnl
            celle_glob(i) = iion(chgglobnl(i))
         end do
         call build_cell_list2(celle_glob,celle_key,nionlocnl,
     &               sqrt(cutb2),c_bP)
      end if

c
c     Build Adjacency list piece by piece
c
      do iter = 0,niterMatb-1

        offsetlMb = iter*nbMatb
        offsetrMb = min((iter+1)*nbMatb,nblock)
        if (iter.eq.niterMatb-1) nbMatb=nblock-offsetlMb

        if (octahedron) then
           call build_adjacency_matrix_octahedron(nionlocnl,c_bP)
        else
           call build_adjacency_matrix(nionlocnl,c_bP)
        end if

        call pre_process_adjacency_matrix (nblock,nb_pair)
        nionlocnlb_pair = nionlocnlb_pair + nb_pair

        call prmem_mvrequest (eblst,nionlocnlb_pair*(BLOCK_SIZE+2),
     &          async=.true.)
        call fill_adjacency_list(nblock,nionlocnlb_pair,eblst)

      end do

      call prmem_request (ieblst,nionlocnlb_pair,async=.true.)

      ! zero list and ilist
      call set_to_zero1_int(eblst(nionlocnlb_pair*2+1),
     &          nionlocnlb_pair*(BLOCK_SIZE),rec_queue)
      call set_to_zero1_int(ieblst,nionlocnlb_pair,rec_queue)

      if (use_shortclist) then
         call prmem_request (shorteblst,nionlocnlb_pair*
     &                         (BLOCK_SIZE+2),async=.true.)
         call prmem_request (ishorteblst,nionlocnlb_pair,
     &                         async=.true.)
         call set_to_zero1_int(shorteblst,nionlocnlb_pair*
     &                         (2+BLOCK_SIZE),rec_queue)
         call set_to_zero1_int(ishorteblst,nionlocnlb_pair,
     &                         rec_queue)
      end if

      call Init_blocklen(nblock)
      call set_ChgData_CellOrder(.true.)
c
c     Debug print
c
      if (deb_Path) then
c13   format( I3," iteration(s) required to build Adj lst" )
c14   format( I10," nblocks process at once over",I10 )
c15   format(1x,'Size Adj Mat (Mo)',F12.3,';',4x,'buffer (Go)',F10.6)
 16   format(1x,'nionlocnl',I8,';',3x,'clist mem space (Mo)',F12.3)
         !nbMatb = max(1,min(int(buffMatb*1d9/(4*szMatb)),nblock))
         szolst = nionlocnlb_pair*(BLOCK_SIZE+3)*szoi
         if (use_shortclist) szolst = 2*szolst
         !print 13, niterMatb
         !print 14, nbMatb,nblock
         !print 15, int(szMatb*4,8)*nbMatb/1d6, buffMatb
         print 16, nionlocnl, szolst/1d6
      end if

c
c     From Block-Block to Block-Atoms list
c
#ifdef _CUDA
      if (use_shortclist) then
!$acc host_data use_device(celle_glob,cell_scan,cell_len,cell_lenr,
!$acc&    celle_x,celle_y,celle_z,matb_lst,eblst,shorteblst,
!$acc&    cell_len2,cell_len2r)
      call filter_shortlst_sparse1 <<<*,128,0,rec_stream>>>
     &     (celle_glob,cell_scan,celle_x,celle_y,celle_z,matb_lst,
     &      nionlocnlb,nblock,szMatb,nionlocnlb_pair,
     &      cshortbuf2,0.0_ti_p,cutb2,
     &      eblst,cell_len,cell_lenr,shorteblst,cell_len2,cell_len2r)
      call check_launch_kernel(" clist:filter_shortlst_sparse1")
!$acc end host_data
      else
!$acc host_data use_device(celle_glob,cell_scan,cell_len,cell_lenr,
!$acc&    celle_x,celle_y,celle_z,matb_lst,eblst)

c     call cu_filter_lst_sparse
c    &     (celle_glob,cell_scan,celle_x,celle_y,celle_z,matb_lst,
c    &      nionlocnlb,nblock,szMatb,nionlocnlb_pair,cutb2,
c    &      eblst,cell_len,xcell,ycell,zcell,xcell2,ycell2,zcell2,
c    &      rec_stream)
c     call filter_lst_sparse<<<*,128,0,rec_stream>>>
c    &     (celle_glob,cell_scan,celle_x,celle_y,celle_z,matb_lst,
c    &      nionlocnlb,nblock,szMatb,nionlocnlb_pair,mbuf2,
c    &      eblst,cell_len)
      call filter_lst_sparse1<<<*,128,0,rec_stream>>>
     &     (celle_glob,cell_scan,celle_x,celle_y,celle_z,matb_lst,
     &      nionlocnlb,nblock,szMatb,nionlocnlb_pair,cutb2,
     &      eblst,cell_len,cell_lenr)
      call check_launch_kernel(" clist:filter_lst_sparse1")

!$acc end host_data
      end if
#else
      !TODO Add an OpenACC version of this kernel
#endif

      call finalize_list_C2
     &     (nblock,nionlocnlb_pair,nionlocnlb
     &     ,nionlocnlb2_pair,ieblst,eblst,cell_len,cell_scan
     &     ,cell_lenr)

      if (use_shortclist) call finalize_list_C2
     &        (nblock,nionlocnlb_pair,nionlocnlb
     &        ,nshortionlocnlb2_pair,ishorteblst,shorteblst
     &        ,cell_len2,cell_scan,cell_len2r)

      end subroutine

c
c    "vlistcell" performs a complete rebuild of the
c     vdw neighbor lists for charges using linked cells method
c
      subroutine vlistcellgpu
      use atmlst
      use atoms
      use domdec
      use inform ,only: deb_Path
      use iounit
      use kvdws
      use nblist_b_inl
      use neigh
      use vdw
      use utilgpu
      use mpi
      implicit none
      integer iglob,iloc
      integer i,ii,icell,j,k,nneigloc
      integer ineig,iivdw,iv
      integer inumneig,ibufbeg,cap
      integer kcell,kloc,kglob,kbis
      integer ncell_loc,kcell_len
      real(t_p) xr,yr,zr,xi,yi,zi,xk,yk,zk,r2,rdn
      real(t_p) zbeg,zend,ybeg,yend,xbeg,xend
      real(t_p) pos1,pos2,pos3
      real(t_p) xred(nbloc)
      real(t_p) yred(nbloc)
      real(t_p) zred(nbloc)
      logical docompute
!$acc routine(fatal_acc)
c
      if(deb_Path) write(*,*) 'vlistcellgpu'
c
!$acc data create(xred,yred,zred)
!$acc&     present(vdwglob,vdwglobnl,ivdw,jvdw,x,y,z,
!$acc&  vlst,nvlst,ired,kred,loc,rad,indcell,neigcell,
!$acc&  cell_len,bufbegcell,numneigcell,repartcell)
!$acc&     async(rec_queue)

      zbeg = zbegproc(rank+1)
      zend = zendproc(rank+1)
      ybeg = ybegproc(rank+1)
      yend = yendproc(rank+1)
      xbeg = xbegproc(rank+1)
      xend = xendproc(rank+1)
c
c     apply reduction factors to find coordinates for each site
c

!$acc parallel loop async(rec_queue)
      do ii = 1, nvdwbloc
         iivdw   = vdwglob(ii)
         iglob   = ivdw(iivdw)
         i       = loc (iglob)
         iv      = ired(iglob)
         rdn     = kred(iglob)
         xred(i) = rdn*(x(iglob)-x(iv)) + x(iv)
         yred(i) = rdn*(y(iglob)-y(iv)) + y(iv)
         zred(i) = rdn*(z(iglob)-z(iv)) + z(iv)
      end do
c
c     perform a complete list build
c
!$acc parallel loop vector_length(32) private(j) async(rec_queue)
      do i = 1, nvdwlocnl
        iivdw    = vdwglobnl(i)
        iglob    = ivdw(iivdw)
        icell    = repartcell(iglob)
        iloc     = loc(iglob)
        inumneig = numneigcell(icell)
        j        = 0
        xi       = xred(iloc)
        yi       = yred(iloc)
        zi       = zred(iloc)
c
c     align data of the local cell and the neighboring ones
c
!$acc loop seq
        do ineig = 0, inumneig
           if (ineig.ne.0) then
              kcell = neigcell(ineig,icell)
           else
              kcell = icell
           end if
           kcell_len = cell_len(kcell)
c          if (ii.gt.cell_len(kcell)) cycle
           ibufbeg   = bufbegcell(kcell)
c
c     do the neighbor search
c
!$acc loop vector
           do ii = 1,kcell_len
              kglob = indcell(ibufbeg+ii-1)
              if (kglob.le.iglob) cycle
              if (rad(jvdw(kglob)).eq.0) cycle
              kbis  = loc(kglob)
              xk    = xred(kbis)
              yk    = yred(kbis)
              zk    = zred(kbis)
              pos1  = xi - xk
              pos2  = yi - yk
              pos3  = zi - zk
              call midpointimage_inl(xk,yk,zk,pos1,pos2,pos3)
c             docompute = .true.
              if    ((zk.lt.zbeg).or.(zk.ge.zend)
     $           .or.(yk.lt.ybeg).or.(yk.ge.yend)
     $           .or.(xk.lt.xbeg).or.(xk.ge.xend)) cycle
c
c     compute the distances and build the list accordingly
c
              r2  = pos1**2 + pos2**2 + pos3**2
              if (r2.gt.vbuf2) cycle
!$acc atomic capture
              j   = j + 1
              cap = j
!$acc end atomic
              vlst(cap,i) = kglob
           end do
        end do
        nvlst(i) = j
c
c     check to see if the neighbor list is too long
c
        if (j .ge. maxvlst) then
           if (rank.eq.0) then
             print*,' VBUILD  --  Too many Neighbors;',
     &                  ' Increase MAXVLST'
             call fatal_acc
           end if
        end if
      end do

!$acc end data
      end
c
c    "vlistcell2gpu" performs a complete rebuild of the
c     short range and regular vdw neighbor lists for charges 
c      using linked cells method
c
      subroutine vlistcell2gpu
      use atmlst
      use atoms
      use cutoff
      use domdec
      use inform ,only: deb_Path
      use iounit
      use kvdws
      use nblist_b_inl
      use neigh
      use vdw
      use utilgpu
      use mpi
      implicit none
      integer iglob,iloc
      integer i,ii,icell,j,k,nneigloc
      integer ineig,iivdw,iv
      integer inumneig,ibufbeg,cap
      integer kcell,kloc,kglob,kbis
      integer ncell_loc,kcell_len
      real(t_p) xr,yr,zr,xi,yi,zi,xk,yk,zk,r2,rdn
      real(t_p) zbeg,zend,ybeg,yend,xbeg,xend
      real(t_p) pos1,pos2,pos3
      real(t_p) xred(nbloc)
      real(t_p) yred(nbloc)
      real(t_p) zred(nbloc)
      real(t_p) vbufbeg2
      logical docompute
!$acc routine(fatal_acc)
c
      if(deb_Path) write(*,*) 'vlistcell2gpu'
c
!$acc data create(xred,yred,zred)
!$acc&     present(vdwglob,vdwglobnl,ivdw,jvdw,x,y,z,
!$acc&  vlst,nvlst,ired,kred,loc,rad,indcell,neigcell,
!$acc&  cell_len,bufbegcell,numneigcell,repartcell,
!$acc&  shortvlst,nshortvlst)
!$acc&     async(rec_queue)

      zbeg = zbegproc(rank+1)
      zend = zendproc(rank+1)
      ybeg = ybegproc(rank+1)
      yend = yendproc(rank+1)
      xbeg = xbegproc(rank+1)
      xend = xendproc(rank+1)
c
c     starting distances for long range real space interactions
c
      vbufbeg2 = (vdwshortcut-lbuffer-shortheal)**2
c
c     apply reduction factors to find coordinates for each site
c
!$acc parallel loop async(rec_queue)
      do ii = 1, nvdwbloc
         iivdw   = vdwglob(ii)
         iglob   = ivdw(iivdw)
         i       = loc (iglob)
         iv      = ired(iglob)
         rdn     = kred(iglob)
         xred(i) = rdn*(x(iglob)-x(iv)) + x(iv)
         yred(i) = rdn*(y(iglob)-y(iv)) + y(iv)
         zred(i) = rdn*(z(iglob)-z(iv)) + z(iv)
      end do
c
c     perform a complete list build
c
!$acc parallel loop vector_length(32) private(j,k) async(rec_queue)
      do i = 1, nvdwlocnl
        iivdw    = vdwglobnl(i)
        iglob    = ivdw(iivdw)
        icell    = repartcell(iglob)
        iloc     = loc(iglob)
        inumneig = numneigcell(icell)
        j        = 0
        k        = 0
        xi       = xred(iloc)
        yi       = yred(iloc)
        zi       = zred(iloc)
c
c     align data of the local cell and the neighboring ones
c
!$acc loop seq
        do ineig = 0, inumneig
           if (ineig.ne.0) then
              kcell = neigcell(ineig,icell)
           else
              kcell = icell
           end if
           kcell_len = cell_len(kcell)
c          if (ii.gt.cell_len(kcell)) cycle
           ibufbeg   = bufbegcell(kcell)
c
c     do the neighbor search
c
!$acc loop vector
           do ii = 1,kcell_len
              kglob = indcell(ibufbeg+ii-1)
              if (kglob.le.iglob) cycle
              if (rad(jvdw(kglob)).eq.0) cycle
              kbis  = loc(kglob)
              xk    = xred(kbis)
              yk    = yred(kbis)
              zk    = zred(kbis)
              pos1  = xi - xk
              pos2  = yi - yk
              pos3  = zi - zk
              call midpointimage_inl(xk,yk,zk,pos1,pos2,pos3)
c             docompute = .true.
              if    ((zk.lt.zbeg).or.(zk.ge.zend)
     $           .or.(yk.lt.ybeg).or.(yk.ge.yend)
     $           .or.(xk.lt.xbeg).or.(xk.ge.xend)) cycle
c
c     compute the distances and build the list accordingly
c
              r2  = pos1**2 + pos2**2 + pos3**2
              if (r2.le.vshortbuf2) then
!$acc atomic capture
                 k   = k + 1
                 cap = k
!$acc end atomic
                 shortvlst(cap,i) = kglob
              end if
              if (r2.le.vbuf2) then
!$acc atomic capture
                 j   = j + 1
                 cap = j
!$acc end atomic
                 vlst(cap,i) = kglob
              end if
           end do
        end do
        nshortvlst(i) = k
        nvlst(i)      = j
c
c     check to see if the neighbor list is too long
c
        if (j .ge. maxvlst) then
           if (rank.eq.0) then
             print*,' VBUILD  --  Too many Neighbors;',
     &                  ' Increase MAXVLST'
             call fatal_acc
           end if
        end if
      end do

!$acc end data
      end
c
c    "mlistcellgpu" rebuilds the regular electrostatic neighbor lists
c     for  multipoles using linked cells method
c
      subroutine mlistcellgpu
      use sizes
      use atmlst
      use atoms
      use cell
      use domdec
      use inform ,only: deb_Path
      use iounit
      use mpole
      use nblist_b_inl
      use neigh
      use mpi
      use utilgpu
      implicit none
      integer iglob
      integer i,ii,j,jj,k,icell,nneigloc
      integer ineig,iipole,kkpole
      integer kcell,kloc,kglob
      integer ncell_loc
      integer ibufbeg,kbufbeg,kcell_len
      integer inumneig
      real(t_p) pos1,pos2,pos3
      real(t_p) xi,yi,zi,xk,yk,zk,r2
      real(t_p) zbeg,zend,ybeg,yend,xbeg,xend
      logical docompute
!$acc routine(fatal_acc)

      if(deb_Path) write(*,*) 'mlistcellgpu'
!$acc data present(poleglobnl,ipole,repartcell,cell_len,bufbegcell,
!$acc&  numneigcell,neigcell,indcell,elst,nelst,pollist,x,y,z)

      xbeg = xbegproc(rank+1)
      xend = xendproc(rank+1)
      ybeg = ybegproc(rank+1)
      yend = yendproc(rank+1)
      zbeg = zbegproc(rank+1)
      zend = zendproc(rank+1)
c
c     perform a complete list build
c
c     print*,'mean cell_len ',sum(cell_len)/size(cell_len)
!$acc parallel loop vector_length(32) private(j) async(rec_queue)
      do i = 1, npolelocnl
         iipole    = poleglobnl(i)
         iglob     = ipole(iipole)
         icell     = repartcell(iglob)
         inumneig  = numneigcell(icell)
         j         = 0
         xi        = x(iglob)
         yi        = y(iglob)
         zi        = z(iglob)
c
c        align data of the local cell and the neighboring ones
c
         do ineig = 0,inumneig
            if (ineig.ne.0) then
               kcell  = neigcell(ineig,icell)
            else
               kcell  = icell
            end if
            kcell_len = cell_len(kcell)
c           if (ii.gt.cell_len(kcell)) cycle
            ibufbeg   = bufbegcell(kcell)
c
c     do the neighbor search
c 
!$acc loop vector
            do ii=1,kcell_len
               kglob  = indcell(ibufbeg+ii-1)
               kkpole = pollist(kglob)
c
c     skip atom if it is not in the multipole list
c
               if (kkpole.eq.0.or.kglob.le.iglob) cycle
               xk     =  x(kglob)
               yk     =  y(kglob)
               zk     =  z(kglob)
               pos1   =  xi - xk
               pos2   =  yi - yk
               pos3   =  zi - zk
               call midpointimage_inl(xk,yk,zk,pos1,pos2,pos3)
               if ((zk.lt.zbeg).or.(zk.ge.zend)
     &         .or.(yk.lt.ybeg).or.(yk.ge.yend)
     &         .or.(xk.lt.xbeg).or.(xk.ge.xend)) cycle
c
c     compute the distances and build the list accordingly
c
               r2 = pos1**2 + pos2**2 + pos3**2
               if (r2.le.mbuf2) then
!$acc atomic capture
                  j          = j + 1
                  jj         = j
!$acc end atomic
c                 kkpole     = pollist(kglob)
                  elst(jj,i) = kkpole
               end if
            end do
         end do
         nelst(i) = j
c
c     check to see if the neighbor list is too long
c
         if (j .ge. maxelst) then
            if (rank.eq.0) then
              print*,' MBUILD  --  Too many Neighbors;',
     &                   ' Increase MAXELST'
              call fatal_acc
            end if
         end if
      end do

!$acc end data
      return
      end
c
c    "mlistcell2gpu" performs a complete rebuild of the
c     short range and regular electrostatic neighbor lists for 
c     multipoles using linked cells method
c
      subroutine mlistcell2gpu
      use sizes
      use atmlst
      use atoms
      use cell
      use domdec
      use inform ,only: deb_Path
      use iounit
      use mpole
      use nblist_b_inl
      use neigh
      use mpi
      use utilgpu
      implicit none
      integer iglob
      integer i,ii,j,cap,k,kk,icell,nneigloc
      integer ineig,iipole,kkpole
      integer kcell,kloc,kglob
      integer ncell_loc
      integer ibufbeg,kbufbeg,kcell_len
      integer inumneig
      real(t_p) pos1,pos2,pos3
      real(t_p) xi,yi,zi,xk,yk,zk,r2
      real(t_p) zbeg,zend,ybeg,yend,xbeg,xend
      real(t_p) mbufbeg2
      logical docompute
!$acc routine(fatal_acc)

      if(deb_Path) write(*,*) 'mlistcell2gpu'
!$acc data present(poleglobnl,ipole,repartcell,cell_len,bufbegcell,
!$acc&  numneigcell,neigcell,indcell,elst,nelst,pollist,x,y,z,
!$acc&  nshortelst,shortelst)

      xbeg = xbegproc(rank+1)
      xend = xendproc(rank+1)
      ybeg = ybegproc(rank+1)
      yend = yendproc(rank+1)
      zbeg = zbegproc(rank+1)
      zend = zendproc(rank+1)
c
c     perform a complete list build
c
c     print*,'mean cell_len ',sum(cell_len)/size(cell_len)
!$acc parallel loop vector_length(32) private(j,k) async(rec_queue)
      do i = 1, npolelocnl
         iipole    = poleglobnl(i)
         iglob     = ipole(iipole)
         icell     = repartcell(iglob)
         inumneig  = numneigcell(icell)
         j         = 0
         k         = 0
         xi        = x(iglob)
         yi        = y(iglob)
         zi        = z(iglob)
c
c        align data of the local cell and the neighboring ones
c
         do ineig = 0,inumneig
            if (ineig.ne.0) then
               kcell  = neigcell(ineig,icell)
            else
               kcell  = icell
            end if
            kcell_len = cell_len(kcell)
c           if (ii.gt.cell_len(kcell)) cycle
            ibufbeg   = bufbegcell(kcell)
c
c     do the neighbor search
c 
!$acc loop vector
            do ii=1,kcell_len
               kglob  = indcell(ibufbeg+ii-1)
               kkpole = pollist(kglob)
c
c     skip atom if it is not in the multipole list
c
               if (kkpole.eq.0.or.kglob.le.iglob) cycle
               xk     =  x(kglob)
               yk     =  y(kglob)
               zk     =  z(kglob)
               pos1   =  xi - xk
               pos2   =  yi - yk
               pos3   =  zi - zk
               call midpointimage_inl(xk,yk,zk,pos1,pos2,pos3)
               if ((zk.lt.zbeg).or.(zk.ge.zend)
     &         .or.(yk.lt.ybeg).or.(yk.ge.yend)
     &         .or.(xk.lt.xbeg).or.(xk.ge.xend)) cycle
c
c     compute the distances and build the list accordingly
c
               r2 = pos1**2 + pos2**2 + pos3**2
               if (r2.le.mshortbuf2) then
!$acc atomic capture
                  k          = k + 1
                  cap        = k
!$acc end atomic
                  shortelst(cap,i) = kkpole
               end if
               if (r2.le.mbuf2) then
!$acc atomic capture
                  j          = j + 1
                  cap        = j
!$acc end atomic
c                 kkpole     = pollist(kglob)
                  elst(cap,i)= kkpole
               end if
            end do
         end do
         nelst(i)      = j
         nshortelst(i) = k
c
c     check to see if the neighbor list is too long
c
         if (j .ge. maxelst) then
            if (rank.eq.0) then
              print*,' MBUILD  --  Too many Neighbors;',
     &                   ' Increase MAXELST'
              call fatal_acc
            end if
         end if
      end do

!$acc end data
      end
c
c    "clistcell" performs a complete rebuild of the
c     electrostatic neighbor lists for charges using linked cells method
c
      subroutine clistcellgpu
      use sizes
      use atmlst
      use atoms
      use charge
      use domdec
      use inform ,only: deb_Path
      use iounit
      use nblist_b_inl
      use neigh
      use utilgpu
      use mpi
      implicit none
      integer iglob
      integer i,ii,icell,j,k,nneigloc,cap
      integer ineig,inumneig,ibufbeg,iichg,kkchg
      integer kcell,kcell_len,kloc,kglob
      integer ncell_loc
      real(t_p) xr,yr,zr,xi,yi,zi,xk,yk,zk,r2
      real(t_p) zbeg,zend,ybeg,yend,xbeg,xend
      real(t_p) pos1,pos2,pos3
      logical docompute
!$acc routine(fatal_acc)
c
      if(deb_Path) write(*,*) 'clistcellgpu'

      zbeg = zbegproc(rank+1)
      zend = zendproc(rank+1)
      ybeg = ybegproc(rank+1)
      yend = yendproc(rank+1)
      xbeg = xbegproc(rank+1)
      xend = xendproc(rank+1)
c
c     perform a complete list build
c
!$acc parallel loop vector_length(32) private(j) async(rec_queue)
!$acc&         present(chgglobnl,chglist,repartcell,cell_len
!$acc&  ,bufbegcell,numneigcell,neigcell,indcell,elst,nelst,x,y,z)
      do i = 1, nionlocnl
         iichg  = chgglobnl(i)
         iglob  = iion(iichg)
         icell  = repartcell(iglob)
         xi     = x(iglob)
         yi     = y(iglob)
         zi     = z(iglob)
         j      = 0
c
c        align data of the local cell and the neighboring ones
c
         do ineig = 0, numneigcell(icell)
            if (ineig.ne.0) then
               kcell = neigcell(ineig,icell)
            else
               kcell = icell
            end if
            kcell_len = cell_len(kcell)
            ibufbeg   = bufbegcell(kcell)
c
c     do the neighbor search
c
!$acc loop vector
            do ii = 1,kcell_len
               kglob = indcell(ibufbeg+ii-1)
               kkchg = chglist(kglob)
               if (kkchg.eq.0.or.kglob.le.iglob) cycle
               xk    = x(kglob)
               yk    = y(kglob)
               zk    = z(kglob)
               pos1  = xi - xk
               pos2  = yi - yk
               pos3  = zi - zk
               call midpointimage_inl(xk,yk,zk,pos1,pos2,pos3)

               if  ((zk.lt.zbeg).or.(zk.ge.zend)
     $          .or.(yk.lt.ybeg).or.(yk.ge.yend)
     $          .or.(xk.lt.xbeg).or.(xk.ge.xend)) cycle
c                docompute = .false.
c
c        compute the distances and build the list accordingly
c
               r2  = pos1**2 + pos2**2 + pos3**2
               if (r2 .gt. cbuf2) cycle
!$acc atomic capture
               j           = j + 1
               cap         = j
!$acc end atomic
               elst(cap,i) = kkchg
            end do
         end do
         nelst(i) = j
c
c     check to see if the neighbor list is too long
c
         if (j .ge. maxelst) then
            if (rank.eq.0) then
               print*,' CBUILD  --  Too many Neighbors;',
     &                ' Increase MAXELST'
               call fatal_acc
            end if
         end if
      end do
      end
c
c    "clistcell2gpu" performs a complete rebuild of the
c     electrostatic short range and regular neighbor lists for charges 
c     using linked cells method
c
      subroutine clistcell2gpu
      use sizes
      use atmlst
      use atoms
      use charge
      use cutoff
      use domdec
      use inform ,only: deb_Path
      use iounit
      use nblist_b_inl
      use neigh
      use utilgpu
      use mpi
      implicit none
      integer iglob
      integer i,ii,icell,j,k,nneigloc,cap
      integer ineig,inumneig,ibufbeg,iichg,kkchg
      integer kcell,kcell_len,kloc,kglob
      integer ncell_loc
      real(t_p) xr,yr,zr,xi,yi,zi,xk,yk,zk,r2
      real(t_p) zbeg,zend,ybeg,yend,xbeg,xend
      real(t_p) pos1,pos2,pos3
      real(t_p) cbufbeg2
      logical docompute
!$acc routine(fatal_acc)
c
      if(deb_Path) write(*,*) 'clistcell2gpu'

      zbeg = zbegproc(rank+1)
      zend = zendproc(rank+1)
      ybeg = ybegproc(rank+1)
      yend = yendproc(rank+1)
      xbeg = xbegproc(rank+1)
      xend = xendproc(rank+1)
c
c     starting distances for long range real space interactions
c
      cbufbeg2 = (chgshortcut-lbuffer-shortheal)**2
c
c     perform a complete list build
c
!$acc parallel loop vector_length(32) private(j,k) async(rec_queue)
!$acc&         present(chgglobnl,chglist,repartcell,cell_len,bufbegcell,
!$acc&  numneigcell,neigcell,indcell,elst,nelst,nshortelst,
!$acc&  shortelst,x,y,z)
      do i = 1, nionlocnl
         iichg  = chgglobnl(i)
         iglob  = iion(iichg)
         icell  = repartcell(iglob)
         xi     = x(iglob)
         yi     = y(iglob)
         zi     = z(iglob)
         k      = 0
         j      = 0
c
c        align data of the local cell and the neighboring ones
c
         do ineig = 0, numneigcell(icell)
            if (ineig.ne.0) then
               kcell = neigcell(ineig,icell)
            else
               kcell = icell
            end if
            kcell_len = cell_len(kcell)
            ibufbeg   = bufbegcell(kcell)
c
c     do the neighbor search
c
!$acc loop vector
            do ii = 1,kcell_len
               kglob = indcell(ibufbeg+ii-1)
               kkchg = chglist(kglob)
               if (kkchg.eq.0.or.kglob.le.iglob) cycle
               xk    = x(kglob)
               yk    = y(kglob)
               zk    = z(kglob)
               pos1  = xi - xk
               pos2  = yi - yk
               pos3  = zi - zk
               call midpointimage_inl(xk,yk,zk,pos1,pos2,pos3)

               if  ((zk.lt.zbeg).or.(zk.ge.zend)
     $          .or.(yk.lt.ybeg).or.(yk.ge.yend)
     $          .or.(xk.lt.xbeg).or.(xk.ge.xend)) cycle
c                docompute = .false.
c
c        compute the distances and build the list accordingly
c
               r2  = pos1**2 + pos2**2 + pos3**2
               if (r2 .le. cshortbuf2) then
!$acc atomic capture
                  k           = k + 1
                  cap         = k
!$acc end atomic
                  shortelst(cap,i) = kkchg
               end if
               if (r2.le.cbuf2.and.r2.ge.cbufbeg2) then
!$acc atomic capture
                  j           = j + 1
                  cap         = j
!$acc end atomic
                  elst(cap,i) = kkchg
               end if
            end do
         end do
         nshortelst(i) = k
         nelst     (i) = j
c
c     check to see if the neighbor list is too long
c
         if (j .ge. maxelst) then
            if (rank.eq.0) then
               print*,' CBUILD  --  Too many Neighbors;',
     &                    ' Increase MAXELST'
               call fatal_acc
            end if
         end if
      end do
      end
