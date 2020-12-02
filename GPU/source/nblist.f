c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine nblist  --  maintain pairwise neighbor lists  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "nblist" constructs and maintains nonbonded pair neighbor lists
c     for vdw and electrostatic interactions
c
c
#include "tinker_precision.h"
      module nblist_inl
        use tinheader,only:ti_p
        contains
#include "image.f.inc"
      end module

      subroutine nblist(istep)
      use atoms ,only:n
      use boxes,only: xbox2,ybox2,zbox2
      use cutoff,only:use_clist,use_mlist,use_vlist
     &          ,use_shortvlist,use_shortclist,use_shortmlist
      use cutoff,only:vdwcut,chgcut,mpolecut
      use domdec
      use inform,only: deb_Path,tindPath
      use mpole ,only: npolelocnl
      use neigh ,only: ineigup
     &          ,mlst_en=>mlst_enable,mlst2_en=>mlst2_enable
     &          ,clst_en=>clst_enable,clst2_en=>clst2_enable
     &          ,vlst_en=>vlst_enable,vlst2_en=>vlst2_enable
      use neigh , only: mbuf2,cbuf2,vbuf2,lbuffer
      use potent
      use sizes ,only:tinkerdebug
      use timestat
      use tinMemory,only:print_memory_usage
      use mpi
      use utilgpu,only:rec_queue
      implicit none
      integer istep,modnl
      integer i,j,k
      real(t_p)  time0,time1
      real(t_p) boxedge2
c
c     check the consistency of the cutoffs+buffer compared to the size
c     of the box (minimum image convention)
c 
      boxedge2 = min(xbox2,ybox2,zbox2)
      if (use_clist) then
        if (cbuf2.gt.boxedge2*boxedge2) then
          if (rank.eq.0) then
            print*,'Error in neigbor list: max cutoff + ',
     $   'buffer should be less than half one edge of the box'
            print*,'Charge cutoff = ', chgcut
            print*,'List buffer = ',lbuffer
          end if
          call fatal
        end if
      end if
      if (use_mlist) then
        if (mbuf2.gt.boxedge2*boxedge2) then
          if (rank.eq.0) then
            print*,'Error in neigbor list: max cutoff + ',
     $   'buffer should be less than half one edge of the box'
            print*,'Multipole cutoff = ', mpolecut
            print*,'List buffer = ',lbuffer
          end if
          call fatal
        end if
      end if
      if (use_vlist) then
        if (vbuf2.gt.boxedge2*boxedge2) then
          if (rank.eq.0) then
            print*,'Error in neigbor list: max cutoff + ',
     $   'buffer should be less than half one edge of the box'
            print*,'VDW cutoff = ', vdwcut
            print*,'List buffer = ',lbuffer
          end if
          call fatal
        end if
      end if
c
c     check number of steps between nl updates
c
      modnl = mod(istep,ineigup)
      if (modnl.ne.0) return
      if (btest(tinkerdebug,tindPath)) call print_memory_usage
      !call mHalte

      if (use_mpole.or.use_charge) call SpatialReorder
      if ((use_pmecore).and.(rank.gt.ndir-1)) return
c
c     build the cells at the beginning and assign the particules to them
c
      call timer_enter( timer_nl )
      ! Init classical list build
      if (mlst_en.or.vlst_en.or.clst_en) then
         call init_cell_list_verlet
         call build_cell_listgpu(istep)
      end if
c
      ! Classical list construct
      if (use_shortclist) then
         call timer_enter( timer_clistcell )
         if (clst_en) call clistcell2gpu
         call timer_exit ( timer_clistcell )
      else if (use_clist) then
         call timer_enter( timer_clistcell )
         if (clst_en) call clistcellgpu
         call timer_exit ( timer_clistcell )
      endif
      if (use_shortvlist) then
         call timer_enter( timer_vlistcell )
         if (vlst_en) call vlistcell2gpu
         call timer_exit ( timer_vlistcell )
      else if (use_vlist) then
         call timer_enter( timer_vlistcell )
         if (vlst_en) call vlistcellgpu
         call timer_exit ( timer_vlistcell )
      endif
      if (use_shortmlist) then
         call timer_enter( timer_mlistcell )
         if (mlst_en) call mlistcell2gpu
         call timer_exit ( timer_mlistcell )
      else if (use_mlist) then
         call timer_enter( timer_mlistcell )
         if (mlst_en) call mlistcellgpu
         call timer_exit ( timer_mlistcell )
      endif        

      ! C2 list (Atom-Blocks)
      if (use_vlist.or.use_shortvlist) then
         call timer_enter( timer_vlistcell )
         if (vlst2_en) call vlist_block
         call timer_exit ( timer_vlistcell )
      end if
      if (use_mlist.or.use_shortmlist) then
         call timer_enter( timer_mlistcell )
         if (mlst2_en) call mlist_block
         call timer_exit ( timer_mlistcell )
      end if
      if (use_clist.or.use_shortclist) then
         call timer_enter( timer_clistcell )
         if (clst2_en) call clist_block
         call timer_exit ( timer_clistcell )
      end if

      !call save_atoms_nl

!$acc wait
      call timer_exit( timer_nl )
      end

      subroutine init_cell_list_verlet
      use sizes
      use domdec
      use cutoff
      use mpole     ,only: npolelocnl
      use neigh
      use potent
      use mpi
      use tinMemory ,only:prmem_request,s_prmem,sd_prmem
      use utils     ,only: set_to_zero1_int8
     &              ,set_to_zero1_int
      use utilgpu   ,only:rec_queue
      implicit none
      integer i,j
      integer(8),save:: size_elst,size_vlst
      integer(8) size_selst,size_svlst

      ! Reallocate space fo list
      if (mlst_enable.or.clst_enable) then
         if (.not.allocated(nelst)) then
            allocate (nelst (nlocnl))
            allocate (nelstc(nlocnl))
            allocate ( elst(maxelst,nlocnl))
            size_elst = int(nlocnl,8)*int(maxelst,8)
            s_prmem = s_prmem + 4*nlocnl*int(maxelst+2,8)
#ifdef _OPENACC
           sd_prmem =sd_prmem + 4*nlocnl*int(maxelst+2,8)
#endif
!$acc    enter data create(elst) async
         else if (nlocnl > size(nelst)) then
            s_prmem = s_prmem - sizeof(nelst)*int(maxelst+2,8)
#ifdef _OPENACC
           sd_prmem =sd_prmem - sizeof(nelst)*int(maxelst+2,8)
#endif
!$acc    exit data delete(elst) async
            deallocate (nelst)
            deallocate ( elst)
            deallocate (nelstc)
            allocate (nelst (nlocnl))
            allocate (nelstc(nlocnl))
            allocate ( elst (maxelst,nlocnl))
            size_elst = int(nlocnl,8)*int(maxelst,8)
            s_prmem = s_prmem + 4*nlocnl*int(maxelst+2,8)
#ifdef _OPENACC
           sd_prmem =sd_prmem + 4*nlocnl*int(maxelst+2,8)
#endif
!$acc    enter data create(elst) async
         end if
      end if

      if (vlst_enable) then
         if (.not.allocated(nvlst)) then
            allocate (nvlst (nlocnl))
            allocate ( vlst(maxvlst,nlocnl))
            size_vlst = int(nlocnl,8)*int(maxvlst,8)
            s_prmem = s_prmem + 4*nlocnl*int(maxvlst+2,8)
#ifdef _OPENACC
           sd_prmem =sd_prmem + 4*nlocnl*int(maxvlst+2,8)
#endif
!$acc    enter data create(vlst) async
         else if (nlocnl > size(nvlst)) then
            s_prmem = s_prmem - sizeof(nelst)*int(maxvlst+2,8)
#ifdef _OPENACC
           sd_prmem =sd_prmem - sizeof(nelst)*int(maxvlst+2,8)
#endif
!$acc    exit data delete(vlst) async
            deallocate (nvlst)
            deallocate ( vlst)
            allocate ( nvlst(nlocnl))
            allocate ( vlst(maxvlst,nlocnl))
            size_vlst = int(nlocnl,8)*int(maxvlst,8)
!$acc    enter data create(vlst) async
            s_prmem = s_prmem + 4*nlocnl*int(maxvlst+2,8)
#ifdef _OPENACC
           sd_prmem =sd_prmem + 4*nlocnl*int(maxvlst+2,8)
#endif
         end if
      end if

      if ((mlst_enable.or.clst_enable).and.
     &    (use_shortclist.or.use_shortmlist))
     &   then
        size_elst = int(nlocnl,8)*int(maxelst,8)
        call prmem_request(nshortelst,        nlocnl,async=.true.)
        call prmem_request(nshortelstc,       nlocnl,async=.true.)
        call prmem_request( shortelst,maxelst,nlocnl,async=.true.)
        call set_to_zero1_int (nshortelst,   nlocnl,rec_queue)
        call set_to_zero1_int8( shortelst,size_elst,rec_queue)
      end if

      if (vlst_enable.and.use_shortvlist) then
        size_vlst = int(nlocnl,8)*int(maxvlst,8)
        call prmem_request(nshortvlst,        nlocnl,async=.true.)
        call prmem_request( shortvlst,maxvlst,nlocnl,async=.true.)
        call set_to_zero1_int (nshortvlst,   nlocnl,rec_queue)
        call set_to_zero1_int8( shortvlst,size_vlst,rec_queue)
      end if

      ! set data to zero
      if ((mlst_enable.or.clst_enable).and.vlst_enable) then
!$acc parallel loop async
         do i = 1,nlocnl
            nelst(i) = 0
            nvlst(i) = 0
         end do
      else if (mlst_enable.or.clst_enable) then
         call set_to_zero1_int(nelst,nlocnl,rec_queue)
      else if (vlst_enable) then
         call set_to_zero1_int(nvlst,nlocnl,rec_queue)
      end if

      if (mlst_enable.or.clst_enable)
     &   call set_to_zero1_int8(elst,size_elst,rec_queue)
      if (vlst_enable) call set_to_zero1_int8(vlst,size_vlst,rec_queue)
      end subroutine
c
      subroutine mlistcell
      use sizes
      use atmlst
      use atoms
      use domdec
      use inform ,only: deb_Path
      use iounit
      use mpole
      use neigh
      use mpi
      implicit none
      integer iglob
      integer i,icell,j,k,nneigloc
      integer ineig,iipole,kkpole
      integer kcell,kglob
      integer ncell_loc
      !DIR$ ATTRIBUTES ALIGN:64:: index,indcell_loc
      integer, allocatable :: index(:),indcell_loc(:)
      real(t_p) xr,yr,zr,xi,yi,zi,xk,yk,zk,r2
      !DIR$ ATTRIBUTES ALIGN:64:: pos,r2vec
      real(t_p), allocatable :: pos(:,:),r2vec(:)
      logical docompute
c
!$acc data present(elst,nelst)
!$acc update host(elst(:,:),nelst(:))
!$acc end data
      allocate (index(nbloc))
      allocate (indcell_loc(nbloc))
      allocate(pos(3,nbloc))
      allocate(r2vec(nbloc))

      if(deb_Path) write(*,*) 'mlistcell'
c
c     perform a complete list build
c
      do i = 1, npolelocnl
        iipole = poleglobnl(i)
        iglob  = ipole(iipole)
        icell = repartcell(iglob)
c
c       align data of the local cell and the neighboring ones
c
        ncell_loc = cell_len(icell)
        !DIR$ ASSUME_ALIGNED indcell:64
        indcell_loc(1:ncell_loc) = 
     $  indcell(bufbegcell(icell):(bufbegcell(icell)+cell_len(icell)-1))
        do ineig = 1, numneigcell(icell)
          kcell = neigcell(ineig,icell)
          indcell_loc(ncell_loc+1:(ncell_loc+cell_len(kcell))) = 
     $  indcell(bufbegcell(kcell):(bufbegcell(kcell)+cell_len(kcell)-1))
          ncell_loc = ncell_loc + cell_len(kcell)
        end do
c
c       do the neighbor search
c
        nneigloc = 0 
        xi = x(iglob)
        yi = y(iglob)
        zi = z(iglob)
        do k = 1, ncell_loc
          kglob = indcell_loc(k)
          kkpole = pollist(kglob)
c
c   skip atom if it is not in the multipole list
c
          if (kkpole.eq.0) cycle
          if (kglob.le.iglob) cycle
          xk = x(kglob)
          yk = y(kglob)
          zk = z(kglob)
          pos(1,nneigloc+1) = xi - xk
          pos(2,nneigloc+1) = yi - yk
          pos(3,nneigloc+1) = zi - zk
          call midpointimage(xi,yi,zi,xk,yk,zk,pos(1,nneigloc+1),
     $       pos(2,nneigloc+1),pos(3,nneigloc+1),docompute)
          if (docompute) then
            nneigloc = nneigloc + 1
            index(nneigloc) = kglob
          end if
        end do
c
c       compute the distances and build the list accordingly
c
        r2vec(1:nneigloc) = pos(1,1:nneigloc)*pos(1,1:nneigloc) + 
     $      pos(2,1:nneigloc)*pos(2,1:nneigloc) + 
     $      pos(3,1:nneigloc)*pos(3,1:nneigloc)
        
        j = 0
        do k = 1, nneigloc
          r2 = r2vec(k)
          kglob = index(k)
          if (r2 .le. mbuf2) then
             j = j + 1
             kkpole = pollist(kglob)
             elst(j,i) = kkpole
          end if
        end do
        nelst(i) = j
c
c     check to see if the neighbor list is too long
c
        if (nelst(i) .ge. maxelst) then
           if (rank.eq.0) then
             write (iout,10)
   10        format (/,' MBUILD  --  Too many Neighbors;',
     &                  ' Increase MAXELST')
             call fatal
           end if
        end if
      end do
!$acc data present(elst,nelst)
!$acc update device(elst(:,:),nelst(:)) async
!$acc end data
c
      deallocate (pos)
      deallocate (index)
      deallocate (indcell_loc)
      deallocate (r2vec)
!$acc wait
      return
      end
c
c    "mlistcell2" performs a complete rebuild of the
c     short range and regular electrostatic neighbor lists for 
c     multipoles using linked cells method
c
      subroutine mlistcell2
      use sizes
      use atmlst
      use atoms
      use cutoff
      use domdec
      use inform ,only: deb_Path
      use iounit
      use mpole
      use neigh
      use mpi
      implicit none
      integer iglob
      integer i,icell,j,j1,k,nneigloc
      integer ineig,iipole,kkpole
      integer kcell,kglob
      integer ncell_loc
      integer, allocatable :: index(:),indcell_loc(:)
      real(t_p) xi,yi,zi,xk,yk,zk,r2
      real(t_p) mbufbeg2
      real(t_p), allocatable :: pos(:,:),r2vec(:)
      logical docompute
c
      allocate (index(nbloc))
      allocate (indcell_loc(nbloc))
      allocate(pos(3,nbloc))
      allocate(r2vec(nbloc))

      if(deb_Path) write(*,*) 'mlistcell2'
c
c     starting distances for long range real space interactions
c
      mbufbeg2 = (mpoleshortcut-lbuffer-shortheal)**2
c
c     perform a complete list build
c
      do i = 1, npolelocnl
        iipole = poleglobnl(i)
        iglob  = ipole(iipole)
        icell = repartcell(iglob)
c
c       align data of the local cell and the neighboring ones
c
        ncell_loc = cell_len(icell)
        indcell_loc(1:ncell_loc) = 
     $  indcell(bufbegcell(icell):(bufbegcell(icell)+cell_len(icell)-1))
        do ineig = 1, numneigcell(icell)
          kcell = neigcell(ineig,icell)
          indcell_loc(ncell_loc+1:(ncell_loc+cell_len(kcell))) = 
     $  indcell(bufbegcell(kcell):(bufbegcell(kcell)+cell_len(kcell)-1))
          ncell_loc = ncell_loc + cell_len(kcell)
        end do
c
c       do the neighbor search
c
        nneigloc = 0 
        xi = x(iglob)
        yi = y(iglob)
        zi = z(iglob)
        do k = 1, ncell_loc
          kglob = indcell_loc(k)
          kkpole = pollist(kglob)
c
c   skip atom if it is not in the multipole list
c
          if (kkpole.eq.0) cycle
          if (kglob.le.iglob) cycle
          xk = x(kglob)
          yk = y(kglob)
          zk = z(kglob)
          pos(1,nneigloc+1) = xi - xk
          pos(2,nneigloc+1) = yi - yk
          pos(3,nneigloc+1) = zi - zk
          call midpointimage(xi,yi,zi,xk,yk,zk,pos(1,nneigloc+1),
     $       pos(2,nneigloc+1),pos(3,nneigloc+1),docompute)
          if (docompute) then
            nneigloc = nneigloc + 1
            index(nneigloc) = kglob
          end if
        end do
c
c       compute the distances and build the list accordingly
c
        r2vec(1:nneigloc) = pos(1,1:nneigloc)*pos(1,1:nneigloc) + 
     $      pos(2,1:nneigloc)*pos(2,1:nneigloc) + 
     $      pos(3,1:nneigloc)*pos(3,1:nneigloc)
        
        j = 0
        j1 = 0
        do k = 1, nneigloc
          r2 = r2vec(k)
          kglob = index(k)
          if (r2 .le. mshortbuf2) then
             j1 = j1 + 1
             kkpole = pollist(kglob)
             shortelst(j1,i) = kkpole
          end if
          if (r2.le.mbuf2) then
             j = j + 1
             kkpole = pollist(kglob)
             elst(j,i) = kkpole
          end if
        end do
        nelst(i) = j
        nshortelst(i) = j1
c
c     check to see if the neighbor list is too long
c
        if (nelst(i) .ge. maxelst) then
           if (rank.eq.0) then
             write (iout,10)
   10        format (/,' MBUILD  --  Too many Neighbors;',
     &                  ' Increase MAXELST')
             call fatal
           end if
        end if
      end do
c     
      deallocate (pos)
      deallocate (index)
      deallocate (indcell_loc)
      deallocate (r2vec)
      return
      end
c
c    subroutine initmpipme : build the arrays to communicate direct and reciprocal fields
c    during the calculation of the induced dipoles
c
c
      subroutine initmpipme
      use atmlst
      use domdec
      use mpole
      use pme
      use mpi
      use tinMemory
      implicit none
      integer ierr,iipole
      integer i,iproc,tag,iglob
      integer count1
      integer status(MPI_STATUS_SIZE)
      integer, allocatable :: req(:),req2(:),count(:)
      allocate (req(nproc*nproc))
      allocate (req2(nproc*nproc))
      allocate (count(nproc))
c
      count = 0 
c
c     deal with Direct-Recip communications
c
      call prmem_request(buf1,nblocrecdir,async=.true.)
      call prmem_request(buf2,nblocrecdir,async=.true.)

      call prmem_request(buflen1,nproc,config=mhostonly)
      call prmem_request(buflen2,nproc,config=mhostonly)
      call prmem_request(bufbeg1,nproc,config=mhostonly)
      call prmem_request(bufbeg2,nproc,config=mhostonly)
      buflen1 = 0
      buflen2 = 0
      bufbeg1 = 0
      bufbeg2 = 0

!$acc data present(buf1,buf2)
c
      do i = 1, npolerecloc
        iipole = polerecglob(i)
        iglob = ipole(iipole)
        if (repart(iglob).ne.rank) then
          buflen2(repart(iglob)+1) = buflen2(repart(iglob)+1)+1
        end if
      end do
      count1 = 0
      do iproc = 1, nrecdir_recep1
        if (precdir_recep1(iproc).ne.rank) then
          if (buflen2(precdir_recep1(iproc)+1).ne.0) then
            bufbeg2(precdir_recep1(iproc)+1) = count1 + 1
          else
            bufbeg2(precdir_recep1(iproc)+1) = 1
          end if
          count1 = count1 + buflen2(precdir_recep1(iproc)+1)
        end if
      end do
c
      do i = 1, npolerecloc
        iipole = polerecglob(i)
        iglob = ipole(iipole)
        if (repart(iglob).ne.rank) then
          buf2(bufbeg2(repart(iglob)+1)+count(repart(iglob)+1))=
     $      iipole
          count(repart(iglob)+1) = count(repart(iglob)+1) + 1
        end if
      end do
!$acc update device(buf2,bufbeg2) async
c
c     send and receive sizes of the buffers
c
       do i = 1, nrecdir_send1
         if (precdir_send1(i).ne.rank) then
          tag = nproc*rank + precdir_send1(i) + 1
          call MPI_IRECV(buflen1(precdir_send1(i)+1),1,MPI_INT,
     $   precdir_send1(i),tag,COMM_TINKER,req(tag),ierr)
        end if
      end do
      do i = 1, nrecdir_recep1
        if (precdir_recep1(i).ne.rank) then
          tag = nproc*precdir_recep1(i) + rank + 1
          call MPI_ISEND(buflen2(precdir_recep1(i)+1),1,MPI_INT,
     $     precdir_recep1(i),tag,COMM_TINKER,req(tag),ierr)
        end if
      end do
c
      do i = 1, nrecdir_send1
        if (precdir_send1(i).ne.rank) then
          tag = nproc*rank + precdir_send1(i) + 1
          call MPI_WAIT(req(tag),status,ierr)
        end if
      end do
      do i = 1, nrecdir_recep1
        if (precdir_recep1(i).ne.rank) then
          tag = nproc*precdir_recep1(i) + rank + 1
          call MPI_WAIT(req(tag),status,ierr)
        end if
      end do
      count1 = 0
      do iproc = 1, nrecdir_send1
        if (precdir_send1(iproc).ne.rank) then
          if (buflen1(precdir_send1(iproc)+1).ne.0) then
            bufbeg1(precdir_send1(iproc)+1) = count1 + 1
          else
            bufbeg1(precdir_send1(iproc)+1) = 1
          end if
          count1 = count1 + buflen1(precdir_send1(iproc)+1)
        end if
      end do
!$acc update device(bufbeg1) async
c
c     send and receive list of corresponding indexes
c
      do i = 1, nrecdir_send1
        if (precdir_send1(i).ne.rank) then
          tag = nproc*rank + precdir_send1(i) + 1
          call MPI_IRECV(buf1(bufbeg1(precdir_send1(i)+1)),
     $     buflen1(precdir_send1(i)+1),
     $     MPI_INT,precdir_send1(i),tag,COMM_TINKER,req2(tag),ierr)
        end if
      end do
      do i = 1, nrecdir_recep1
        if (precdir_recep1(i).ne.rank) then
          tag = nproc*precdir_recep1(i) + rank + 1
          call MPI_ISEND(buf2(bufbeg2(precdir_recep1(i)+1)),
     $     buflen2(precdir_recep1(i)+1),MPI_INT,precdir_recep1(i),tag,
     $     COMM_TINKER,req2(tag),ierr)
        end if
      end do
c
      do i = 1, nrecdir_send1
        if (precdir_send1(i).ne.rank) then
          tag = nproc*rank + precdir_send1(i) + 1
          call MPI_WAIT(req2(tag),status,ierr)
        end if
      end do
      do i = 1, nrecdir_recep1
        if (precdir_recep1(i).ne.rank) then
          tag = nproc*precdir_recep1(i) + rank + 1
          call MPI_WAIT(req2(tag),status,ierr)
        end if
      end do
!$acc update device(buf1) async
c
!$acc end data 

      !TODO Remove declare directive on thetai
c     call prmem_request(thetai1,4,bsorder,nlocrec,async=.true.)
c     call prmem_request(thetai2,4,bsorder,nlocrec,async=.true.)
c     call prmem_request(thetai3,4,bsorder,nlocrec,async=.true.)

      if (allocated(thetai1)) deallocate (thetai1)
      if (allocated(thetai2)) deallocate (thetai2)
      if (allocated(thetai3)) deallocate (thetai3)
      allocate (thetai1(4,bsorder,nlocrec))
      allocate (thetai2(4,bsorder,nlocrec))
      allocate (thetai3(4,bsorder,nlocrec))

      deallocate (req)
      deallocate (req2)
      deallocate (count)

      return
      end

      subroutine check_nl_rebuild()
      use atoms
      use domdec
      use mpi
      use neigh
      use nblist_inl
      use usage
      integer i,iglob
      real(t_p) :: move=0,red_move=0
      real(t_p) xpos,ypos,zpos
      move = 0.0

!$acc parallel loop async copy(j)
!$acc&         present(glob,use,x,y,z,xold_nl,yold_nl,zold_nl)
      do i = 1, nloc
         iglob = glob(i)
         if (use(iglob)) then
            xpos = xold_nl(iglob) - x(iglob)
            ypos = yold_nl(iglob) - y(iglob)
            zpos = zold_nl(iglob) - z(iglob)
            move = max(move, xpos**2+ypos**2+zpos**2)
         end if
      end do
!$acc wait
      if (nproc.ne.1) then
         call MPI_ALLREDUCE(move,red_move,1,MPI_TPREC,MPI_MAX,
     &                      COMM_TINKER,i)
         move=red_move
      end if

      if (move .gt. lbuffer**2/4) then
         rebuild_lst = .true.
         if (rank.eq.0) then
            print*,'REBUILD_NBLIST !!!'
            print*,'maximum atom move ',sqrt(move)
         end if
      else
         rebuild_lst = .false.
      end if

      end subroutine

      subroutine save_atoms_nl()
      use atoms
      use domdec
      use neigh
      use usage
      use tinheader,only: ti_p
      use tinMemory
      implicit none
      integer i,iglob
      logical ::first_in=.true.

      if (first_in) then
         call prmem_request(xold_nl,n,async=.true.)
         call prmem_request(yold_nl,n,async=.true.)
         call prmem_request(zold_nl,n,async=.true.)
!$acc parallel loop present(xold_nl,yold_nl,zold_nl) async
         do i = 1,n
            xold_nl(i) = 0.0_ti_p
            yold_nl(i) = 0.0_ti_p
            zold_nl(i) = 0.0_ti_p
         end do
         first_in=.false.
      end if

!$acc parallel loop async
!$acc&         present(glob,x,y,z,xold_nl,yold_nl,zold_nl)
      do i = 1,nbloc
         iglob = glob(i)
         xold_nl(iglob) = x(iglob)
         yold_nl(iglob) = y(iglob)
         zold_nl(iglob) = z(iglob)
      end do

      end subroutine
c
c     subroutine reinitnl : get the number of particules whose nl has to be computed
c     and the associated indexes
c
      subroutine reinitnl(istep)
      use atoms
      use domdec
      use cell    ,only:xcell,ycell,zcell
      use inform  ,only:deb_Path,tindPath,abort
      use neigh
      use utilgpu ,only:rec_queue
      use tinMemory
      use timestat,only: timer_enter,timer_exit,timer_nl
     &            ,quiet_timers
      use sizes   ,only: tinkerdebug
      implicit none
      real(t_p) d,mbuf,vbuf,torquebuf,bigbuf,mcell
      integer iproc,i,iglob,modnl
      integer iloc,istep,idomlen,nloc_cap
      integer ibufbeg
      integer:: s_nlocnl=0
!$acc routine(distprocpart1)
c
      !if (istep.ne.0) call check_nl_rebuild
      if (mod(istep,ineigup).ne.0) return
      if (abort) call fatal
c
      if (deb_Path) write(*,*) ' reinitnl - istep(',istep,')'
      call timer_enter(timer_nl)
      mbuf      = sqrt(mbuf2)
      vbuf      = sqrt(vbuf2)
      torquebuf = mbuf + lbuffer + 2.0
      bigbuf    = max( vbuf,torquebuf )
      mcell     = max(xcell,max(ycell,zcell))
c
      nlocnl    = nloc
c
      if (.not.allocated(ineignl)) then
         call prmem_request(ineignl,n)
!$acc wait
      end if
c
!$acc data copy(nlocnl,abort)
!$acc&     present(ineignl,glob)
!$acc&     async(rec_queue)
c
!$acc parallel loop async(rec_queue)
      do i = 1,n
         if (i.le.nloc) then
            ineignl(i) = glob(i)
         else
            ineignl(i) = 0
         end if
      end do
c
      do iproc = 1, nbig_recep
        idomlen = domlen(pbig_recep(iproc)+1)
        ibufbeg = bufbeg(pbig_recep(iproc)+1)
!$acc parallel loop async(rec_queue)
        do i = 1, idomlen
          iloc  = ibufbeg+i-1
          iglob = glob(iloc)
          call distprocpart1(iglob,rank,d,.true.,x,y,z)
          if (d.le.(bigbuf/2)) then
!$acc atomic capture
            nlocnl = nlocnl + 1
            nloc_cap = nlocnl
!$acc end atomic
            ineignl(nloc_cap) = iglob
c           locnl(iglob) = nloc_cap
          else if (d.ge.mcell) then
             !print*,'distprocpart error',iglob,d
!$acc atomic write
             abort = .true.
          end if
        end do
      end do
!$acc end data
!$acc wait(rec_queue)

      if (btest(tinkerdebug,tindPath).and.nlocnl.gt.s_nlocnl) then
         s_nlocnl = nlocnl
 15      format('iRank',I4,'; get nlocnl increased to',I10,'; nloc',I10)
         write(*,15) rank,nlocnl,nloc
      end if

      if (abort) then
         print*, "errors found in distprocpart1 routine"
      end if

      call timer_exit(timer_nl)
      end
c
c
c     subroutine build_cell_list : build the cells in order to build the non bonded neighbor
c     lists with the cell-list method
c
      subroutine build_cell_list
      use atoms
      use bound
      use domdec
      use neigh
      use mpi
      implicit none
      integer i,proc,icell,j,k,p,q,r,iglob
      integer count,iloc
      integer temp_x,temp_y,temp_z

      integer numneig,tempcell
      real(t_p) xmin,xmax,ymin,ymax,zmin,zmax
      real(t_p) lenx,leny,lenz
      real(t_p) mbuf,vbuf,bigbuf
      real(t_p) lenx_cell,leny_cell,lenz_cell
      real(t_p) xr,yr,zr
      real(t_p) eps1,eps2
      real(t_p), allocatable :: xbegcelltemp(:),ybegcelltemp(:)
      real(t_p), allocatable :: zbegcelltemp(:)
      real(t_p), allocatable :: xendcelltemp(:),yendcelltemp(:)
      real(t_p), allocatable :: zendcelltemp(:)
      integer, allocatable :: filledcell(:),indcelltemp(:)

c
      eps1 = 1.0d-10
      eps2 = 1.0d-8
      mbuf = sqrt(mbuf2)
      vbuf = sqrt(vbuf2)+2.0
      if ((mbuf+lbuffer).gt.vbuf) then
        bigbuf = (mbuf+lbuffer)/3
      else
        bigbuf = vbuf/3
      end if
c
c
c     divide the searching domain in cells of size the multipole cutoff
c
      xmin = xbegproc(rank+1)
      xmax = xendproc(rank+1)
      ymin = ybegproc(rank+1)
      ymax = yendproc(rank+1)
      zmin = zbegproc(rank+1)
      zmax = zendproc(rank+1)
      do i = 1, nbig_recep
        proc = pbig_recep(i)
        if (xbegproc(proc+1).le.xmin) xmin = xbegproc(proc+1)
        if (xendproc(proc+1).ge.xmax) xmax = xendproc(proc+1)
        if (ybegproc(proc+1).le.ymin) ymin = ybegproc(proc+1)
        if (yendproc(proc+1).ge.ymax) ymax = yendproc(proc+1)
        if (zbegproc(proc+1).le.zmin) zmin = zbegproc(proc+1)
        if (zendproc(proc+1).ge.zmax) zmax = zendproc(proc+1)
      end do
c
      lenx = abs(xmax-xmin)
      nx_cell = max(1,int(lenx/(bigbuf)))
      lenx_cell = lenx/nx_cell
      leny = abs(ymax-ymin)
      ny_cell = max(1,int(leny/(bigbuf)))
      leny_cell = leny/ny_cell
      lenz = abs(zmax-zmin)
      nz_cell = max(1,int(lenz/(bigbuf)))
      lenz_cell = lenz/nz_cell
      ncell_tot = nx_cell*ny_cell*nz_cell
c
      allocate (xbegcelltemp(nx_cell))
      allocate (xendcelltemp(nx_cell))
      allocate (ybegcelltemp(ny_cell))
      allocate (yendcelltemp(ny_cell))
      allocate (zbegcelltemp(nz_cell))
      allocate (zendcelltemp(nz_cell))
      if (allocated(xbegcell)) deallocate (xbegcell)
      allocate (xbegcell(ncell_tot))
      if (allocated(ybegcell)) deallocate (ybegcell)
      allocate (ybegcell(ncell_tot))
      if (allocated(zbegcell)) deallocate (zbegcell)
      allocate (zbegcell(ncell_tot))
      if (allocated(xendcell)) deallocate (xendcell)
      allocate (xendcell(ncell_tot))
      if (allocated(yendcell)) deallocate (yendcell)
      allocate (yendcell(ncell_tot))
      if (allocated(zendcell)) deallocate (zendcell)
      allocate (zendcell(ncell_tot))
      if (allocated(neigcell)) deallocate (neigcell)
      allocate (neigcell(400,ncell_tot))
      if (allocated(numneigcell)) deallocate (numneigcell)
      allocate (numneigcell(ncell_tot))
      allocate (filledcell(ncell_tot))
c
      do i = 0, nx_cell-1
        xbegcelltemp(i+1) = xmin + i*lenx_cell
        xendcelltemp(i+1) = xmin + (i+1)*lenx_cell
      end do
      do i = 0, ny_cell-1
        ybegcelltemp(i+1) = ymin + i*leny_cell
        yendcelltemp(i+1) = ymin + (i+1)*leny_cell
      end do
      do i = 0, nz_cell-1
        zbegcelltemp(i+1) = zmin + i*lenz_cell
        zendcelltemp(i+1) = zmin + (i+1)*lenz_cell
      end do
c
c     assign cell
c
      do k = 1, nz_cell
        do j = 1, ny_cell
          do i = 1, nx_cell
              icell = (k-1)*ny_cell*nx_cell+(j-1)*nx_cell+i
              xbegcell(icell) = xbegcelltemp(i)
              xendcell(icell) = xendcelltemp(i)
              ybegcell(icell) = ybegcelltemp(j)
              yendcell(icell) = yendcelltemp(j)
              zbegcell(icell) = zbegcelltemp(k)
              zendcell(icell) = zendcelltemp(k)
              numneig = 0
              filledcell = 0
              filledcell(icell) = 1
c
c              do p = -1,1
c                do q = -1,1
c                  do r = -1,1
              do p = -3,3
                do q = -3,3
                  do r = -3,3
                    if ((p.eq.0).and.(q.eq.0).and.(r.eq.0)) goto 10
c
                    temp_x = p+i
                    temp_y = q+j-1
                    temp_z = r+k-1
c                    if ((i.eq.1).and.(p.eq.-1)) temp_x = nx_cell
c                    if ((i.eq.nx_cell).and.(p.eq.1)) temp_x = 1
c                    if ((j.eq.1).and.(q.eq.-1)) temp_y = ny_cell-1
c                    if ((j.eq.ny_cell).and.(q.eq.1)) temp_y = 0
c                    if ((k.eq.1).and.(r.eq.-1)) temp_z = nz_cell-1
c                    if ((k.eq.nz_cell).and.(r.eq.1)) temp_z = 0

                    if ((i.eq.1).and.(p.eq.-3)) temp_x = nx_cell-2
                    if ((i.eq.1).and.(p.eq.-2)) temp_x = nx_cell-1
                    if ((i.eq.1).and.(p.eq.-1)) temp_x = nx_cell
                    if ((i.eq.2).and.(p.eq.-3)) temp_x = nx_cell-1
                    if ((i.eq.2).and.(p.eq.-2)) temp_x = nx_cell
                    if ((i.eq.3).and.(p.eq.-3)) temp_x = nx_cell

                    if ((i.eq.nx_cell).and.(p.eq.1)) temp_x = 1
                    if ((i.eq.nx_cell).and.(p.eq.2)) temp_x = 2
                    if ((i.eq.nx_cell).and.(p.eq.3)) temp_x = 3
                    if ((i.eq.nx_cell-1).and.(p.eq.2)) temp_x = 1
                    if ((i.eq.nx_cell-1).and.(p.eq.3)) temp_x = 2
                    if ((i.eq.nx_cell-2).and.(p.eq.3)) temp_x = 1

                    if ((j.eq.1).and.(q.eq.-3)) temp_y = ny_cell-3
                    if ((j.eq.1).and.(q.eq.-2)) temp_y = ny_cell-2
                    if ((j.eq.1).and.(q.eq.-1)) temp_y = ny_cell-1
                    if ((j.eq.2).and.(q.eq.-3)) temp_y = ny_cell-2
                    if ((j.eq.2).and.(q.eq.-2)) temp_y = ny_cell-1
                    if ((j.eq.3).and.(q.eq.-3)) temp_y = ny_cell-1

                    if ((j.eq.ny_cell).and.(q.eq.1)) temp_y = 0
                    if ((j.eq.ny_cell).and.(q.eq.2)) temp_y = 1
                    if ((j.eq.ny_cell).and.(q.eq.3)) temp_y = 2
                    if ((j.eq.ny_cell-1).and.(q.eq.2)) temp_y = 0
                    if ((j.eq.ny_cell-1).and.(q.eq.3)) temp_y = 1
                    if ((j.eq.ny_cell-2).and.(q.eq.3)) temp_y = 0

                    if ((k.eq.1).and.(r.eq.-3)) temp_z = nz_cell-3
                    if ((k.eq.1).and.(r.eq.-2)) temp_z = nz_cell-2
                    if ((k.eq.1).and.(r.eq.-1)) temp_z = nz_cell-1
                    if ((k.eq.2).and.(r.eq.-3)) temp_z = nz_cell-2
                    if ((k.eq.2).and.(r.eq.-2)) temp_z = nz_cell-1
                    if ((k.eq.3).and.(r.eq.-3)) temp_z = nz_cell-1
                    if ((k.eq.nz_cell).and.(r.eq.1)) temp_z = 0
                    if ((k.eq.nz_cell).and.(r.eq.2)) temp_z = 1
                    if ((k.eq.nz_cell).and.(r.eq.3)) temp_z = 2
                    if ((k.eq.nz_cell-1).and.(r.eq.2)) temp_z = 0
                    if ((k.eq.nz_cell-1).and.(r.eq.3)) temp_z = 1
                    if ((k.eq.nz_cell-2).and.(r.eq.3)) temp_z = 0

                    tempcell = temp_z*ny_cell*nx_cell+temp_y*nx_cell+
     $                temp_x
                    if (filledcell(tempcell).eq.1) goto 10
                    filledcell(tempcell) = 1
                    numneig = numneig+1
                    neigcell(numneig,icell) = tempcell
 10               continue
                  end do
                end do
              end do
              numneigcell(icell) = numneig
          end do
        end do
      end do
!$acc update device(neigcell,numneigcell) async
      deallocate (filledcell)
      deallocate (xbegcelltemp)
      deallocate (xendcelltemp)
      deallocate (ybegcelltemp)
      deallocate (yendcelltemp)
      deallocate (zbegcelltemp)
      deallocate (zendcelltemp)
c
c     assign the atoms to the cells
c
      if (allocated(cell_len)) deallocate (cell_len)
      allocate (cell_len(ncell_tot))
      if (allocated(indcell)) deallocate (indcell)
      allocate (indcell(n))
      if (allocated(bufbegcell)) deallocate (bufbegcell)
      allocate (bufbegcell(ncell_tot))
      if (allocated(repartcell)) deallocate (repartcell)
      allocate (repartcell(n))
      allocate (indcelltemp(n))
      cell_len = 0
c
      do i = 1, nlocnl
        iglob = ineignl(i)
        xr = x(iglob)
        yr = y(iglob)
        zr = z(iglob)
        if (use_bounds) call image(xr,yr,zr)
        if (abs(xr-xmax).lt.eps1) xr = xr-eps2
        if (abs(yr-ymax).lt.eps1) yr = yr-eps2
        if (abs(zr-zmax).lt.eps1) zr = zr-eps2
        do icell = 1, ncell_tot
          if ((zr.ge.zbegcell(icell)).and.
     $     (zr.lt.zendcell(icell)).and.(yr.ge.ybegcell(icell))
     $    .and.(yr.lt.yendcell(icell)).and.(xr.ge.xbegcell(icell))
     $    .and.(xr.lt.xendcell(icell))) then
            repartcell(iglob) = icell
            cell_len(icell) = cell_len(icell) + 1
            indcelltemp(iglob) = cell_len(icell)
          end if
        end do
      end do

c
      bufbegcell(1) = 1
      count = cell_len(1)
      do icell = 2, ncell_tot
        if (cell_len(icell).ne.0) then
          bufbegcell(icell) = count + 1
        else
          bufbegcell(icell) = 1
        end if
        count = count + cell_len(icell)
      end do
c
      do i = 1, nlocnl
        iglob = ineignl(i)
        icell = repartcell(iglob)
        iloc  = bufbegcell(icell) + indcelltemp(iglob) - 1
        indcell(iloc) = iglob
      end do
!$acc update device(cell_len,bufbegcell,indcell,repartcell) async
      deallocate (indcelltemp)
      return
      end
c
c    "clistcell" performs a complete rebuild of the
c     electrostatic neighbor lists for charges using linked cells method
c
      subroutine clistcell
      use sizes
      use atmlst
      use atoms
      use charge
      use domdec
      use inform ,only: deb_Path
      use iounit
      use neigh
      use mpi
      implicit none
      integer iglob
      integer i,icell,j,k,nneigloc
      integer ineig,iichg,kkchg
      integer kcell,kglob
      integer ncell_loc
      integer, allocatable :: index(:),indcell_loc(:)
      real(t_p) xr,yr,zr,xi,yi,zi,xk,yk,zk,r2
      real(t_p), allocatable :: pos(:,:),r2vec(:)
      logical docompute
!$acc data present(elst,nelst)
!$acc update host(elst(:,:),nelst(:))
!$acc end data
c
      allocate (index(nbloc))
      allocate (indcell_loc(nbloc))
      allocate(pos(3,nbloc))
      allocate(r2vec(nbloc))

      if(deb_Path) write(*,*) 'clistcell'

c
c     perform a complete list build
c
      do i = 1, nionlocnl
        iichg = chgglobnl(i)
        iglob  = iion(iichg)
        icell = repartcell(iglob)
c
c       align data of the local cell and the neighboring ones
c
        ncell_loc = cell_len(icell)
        indcell_loc(1:ncell_loc) = 
     $  indcell(bufbegcell(icell):(bufbegcell(icell)+cell_len(icell)-1))
        do ineig = 1, numneigcell(icell)
          kcell = neigcell(ineig,icell)
          indcell_loc(ncell_loc+1:(ncell_loc+cell_len(kcell))) = 
     $  indcell(bufbegcell(kcell):(bufbegcell(kcell)+cell_len(kcell)-1))
          ncell_loc = ncell_loc + cell_len(kcell)
        end do
c
c       do the neighbor search
c
        nneigloc = 0 
        xi = x(iglob)
        yi = y(iglob)
        zi = z(iglob)
        do k = 1, ncell_loc
          kglob = indcell_loc(k)
          kkchg = chglist(kglob)
          if (kkchg.eq.0) cycle
          if (kglob.le.iglob) cycle
          xk = x(kglob)
          yk = y(kglob)
          zk = z(kglob)
          pos(1,nneigloc+1) = xi - xk
          pos(2,nneigloc+1) = yi - yk
          pos(3,nneigloc+1) = zi - zk
          call midpointimage(xi,yi,zi,xk,yk,zk,pos(1,nneigloc+1),
     $       pos(2,nneigloc+1),pos(3,nneigloc+1),docompute)
          if (docompute) then
            nneigloc = nneigloc + 1
            index(nneigloc) = kglob
          end if
        end do
c
c       compute the distances and build the list accordingly
c
        r2vec(1:nneigloc) = pos(1,1:nneigloc)*pos(1,1:nneigloc) + 
     $      pos(2,1:nneigloc)*pos(2,1:nneigloc) + 
     $      pos(3,1:nneigloc)*pos(3,1:nneigloc)
        
        j = 0
        do k = 1, nneigloc
          r2 = r2vec(k)
          kglob = index(k)
          if (r2 .le. cbuf2) then
             j = j + 1
             kkchg = chglist(kglob)
             elst(j,i) = kkchg
          end if
        end do
        nelst(i) = j
c
c     check to see if the neighbor list is too long
c
        if (nelst(i) .ge. maxelst) then
           if (rank.eq.0) then
             write (iout,10)
   10        format (/,' MBUILD  --  Too many Neighbors;',
     &                  ' Increase MAXELST')
             call fatal
           end if
        end if
      end do
!$acc data present(elst,nelst)
!$acc update device(elst(:,:),nelst(:)) async
!$acc end data
c     
      deallocate (pos)
      deallocate (index)
      deallocate (indcell_loc)
      deallocate (r2vec)
!$acc wait
      return
      end
c
c    "clistcell2" performs a complete rebuild of the
c     electrostatic short range and regular neighbor lists for charges 
c     using linked cells method
c
      subroutine clistcell2
      use sizes
      use atmlst
      use atoms
      use charge
      use cutoff
      use domdec
      use iounit
      use neigh
      use mpi
      implicit none
      integer iglob
      integer i,icell,j,j1,k,nneigloc
      integer ineig,iichg,kkchg
      integer kcell,kglob
      integer ncell_loc
      integer, allocatable :: index(:),indcell_loc(:)
      real(t_p) xi,yi,zi,xk,yk,zk,r2
      real(t_p) cbufbeg2
      real(t_p), allocatable :: pos(:,:),r2vec(:)
      logical docompute
c
      allocate (index(nbloc))
      allocate (indcell_loc(nbloc))
      allocate(pos(3,nbloc))
      allocate(r2vec(nbloc))
c
c     starting distances for long range real space interactions
c
      cbufbeg2 = (chgshortcut-lbuffer-shortheal)**2
c
c     perform a complete list build
c
      do i = 1, nionlocnl
        iichg = chgglobnl(i)
        iglob  = iion(iichg)
        icell = repartcell(iglob)
c
c       align data of the local cell and the neighboring ones
c
        ncell_loc = cell_len(icell)
        indcell_loc(1:ncell_loc) = 
     $  indcell(bufbegcell(icell):(bufbegcell(icell)+cell_len(icell)-1))
        do ineig = 1, numneigcell(icell)
          kcell = neigcell(ineig,icell)
          indcell_loc(ncell_loc+1:(ncell_loc+cell_len(kcell))) = 
     $  indcell(bufbegcell(kcell):(bufbegcell(kcell)+cell_len(kcell)-1))
          ncell_loc = ncell_loc + cell_len(kcell)
        end do
c
c       do the neighbor search
c
        nneigloc = 0 
        xi = x(iglob)
        yi = y(iglob)
        zi = z(iglob)
        do k = 1, ncell_loc
          kglob = indcell_loc(k)
          kkchg = chglist(kglob)
          if (kkchg.eq.0) cycle
          if (kglob.le.iglob) cycle
          xk = x(kglob)
          yk = y(kglob)
          zk = z(kglob)
          pos(1,nneigloc+1) = xi - xk
          pos(2,nneigloc+1) = yi - yk
          pos(3,nneigloc+1) = zi - zk
          call midpointimage(xi,yi,zi,xk,yk,zk,pos(1,nneigloc+1),
     $       pos(2,nneigloc+1),pos(3,nneigloc+1),docompute)
          if (docompute) then
            nneigloc = nneigloc + 1
            index(nneigloc) = kglob
          end if
        end do
c
c       compute the distances and build the list accordingly
c
        r2vec(1:nneigloc) = pos(1,1:nneigloc)*pos(1,1:nneigloc) + 
     $      pos(2,1:nneigloc)*pos(2,1:nneigloc) + 
     $      pos(3,1:nneigloc)*pos(3,1:nneigloc)
        
        j = 0
        j1 = 0
        do k = 1, nneigloc
          r2 = r2vec(k)
          kglob = index(k)
          if (r2 .le. cshortbuf2) then
             j1 = j1 + 1
             kkchg = chglist(kglob)
             shortelst(j1,i) = kkchg
          end if
          if ((r2 .le. cbuf2).and.(r2.ge.cbufbeg2)) then
             j = j + 1
             kkchg = chglist(kglob)
             elst(j,i) = kkchg
          end if
        end do
        nelst(i) = j
        nshortelst(i) = j1
c
c     check to see if the neighbor list is too long
c
        if ((nelst(i) .ge. maxelst).or.(nshortelst(i).ge.maxelst)) then
           if (rank.eq.0) then
             write (iout,10)
   10        format (/,' MBUILD  --  Too many Neighbors;',
     &                  ' Increase MAXELST')
             call fatal
           end if
        end if
      end do
c     
      deallocate (pos)
      deallocate (index)
      deallocate (indcell_loc)
      deallocate (r2vec)
      return
      end
c
c    "vlistcell" performs a complete rebuild of the
c     vdw neighbor lists for charges using linked cells method
c
      subroutine vlistcell
      use atmlst
      use atoms
      use domdec
      use inform ,only: deb_Path
      use iounit
      use kvdws
      use neigh
      use vdw
      use mpi
      implicit none
      integer iglob,iloc
      integer i,ii,icell,j,k,nneigloc
      integer ineig,iivdw,iv
      integer kcell,kglob,kbis
      integer ncell_loc
      integer, allocatable :: index(:),indcell_loc(:)
      real(t_p) xr,yr,zr,xi,yi,zi,xk,yk,zk,r2,rdn
      real(t_p), allocatable :: pos(:,:),r2vec(:)
      real(t_p), allocatable :: xred(:)
      real(t_p), allocatable :: yred(:)
      real(t_p), allocatable :: zred(:)
      logical docompute
c
!$acc data present(vlst,nvlst)
!$acc update host(vlst(:,:),nvlst(:))
!$acc end data
      allocate (xred(nbloc))
      allocate (yred(nbloc))
      allocate (zred(nbloc))
c
      allocate (index(nbloc))
      allocate (indcell_loc(nbloc))
      allocate(pos(3,nbloc))
      allocate(r2vec(nbloc))
      if(deb_Path) write(*,*) 'vlistcell'
c
c     apply reduction factors to find coordinates for each site
c
      do ii = 1, nvdwbloc
         iivdw = vdwglob(ii)
         iglob = ivdw(iivdw)
         i = loc(iglob)
         iv = ired(iglob)
         rdn = kred(iglob)
         xred(i) = rdn*(x(iglob)-x(iv)) + x(iv)
         yred(i) = rdn*(y(iglob)-y(iv)) + y(iv)
         zred(i) = rdn*(z(iglob)-z(iv)) + z(iv)
      end do
c
c     perform a complete list build
c
      do i = 1, nvdwlocnl
        iivdw = vdwglobnl(i)
        iglob  = ivdw(iivdw)
        icell = repartcell(iglob)
        iloc = loc(iglob)
c
c       align data of the local cell and the neighboring ones
c
        ncell_loc = cell_len(icell)
        indcell_loc(1:ncell_loc) = 
     $  indcell(bufbegcell(icell):(bufbegcell(icell)+cell_len(icell)-1))
        do ineig = 1, numneigcell(icell)
          kcell = neigcell(ineig,icell)
          indcell_loc(ncell_loc+1:(ncell_loc+cell_len(kcell))) = 
     $  indcell(bufbegcell(kcell):(bufbegcell(kcell)+cell_len(kcell)-1))
          ncell_loc = ncell_loc + cell_len(kcell)
        end do
c
c       do the neighbor search
c
        nneigloc = 0 
        xi = xred(iloc)
        yi = yred(iloc)
        zi = zred(iloc)
        do k = 1, ncell_loc
          kglob = indcell_loc(k)
          if (kglob.le.iglob) cycle
          if (rad(jvdw(kglob)).eq.0) cycle
          kbis = loc(kglob)
          xk = xred(kbis)
          yk = yred(kbis)
          zk = zred(kbis)
          pos(1,nneigloc+1) = xi - xk
          pos(2,nneigloc+1) = yi - yk
          pos(3,nneigloc+1) = zi - zk
          call midpointimage(xi,yi,zi,xk,yk,zk,pos(1,nneigloc+1),
     $       pos(2,nneigloc+1),pos(3,nneigloc+1),docompute)
          if (docompute) then
            nneigloc = nneigloc + 1
            index(nneigloc) = kglob
          end if
        end do
c
c       compute the distances and build the list accordingly
c
        r2vec(1:nneigloc) = pos(1,1:nneigloc)*pos(1,1:nneigloc) + 
     $      pos(2,1:nneigloc)*pos(2,1:nneigloc) + 
     $      pos(3,1:nneigloc)*pos(3,1:nneigloc)
        
        j = 0
        do k = 1, nneigloc
          r2 = r2vec(k)
          kglob = index(k)
          if (r2 .le. vbuf2) then
             j = j + 1
             vlst(j,i) = kglob
          end if
        end do
        nvlst(i) = j
c
c     check to see if the neighbor list is too long
c
        if (nvlst(i) .ge. maxvlst) then
           if (rank.eq.0) then
             write (iout,10)
   10        format (/,' VBUILD  --  Too many Neighbors;',
     &                  ' Increase MAXVLST')
             call fatal
           end if
        end if
      end do
!$acc data present(vlst,nvlst)
!$acc update device(vlst(:,:),nvlst(:)) async
!$acc end data
c     
      deallocate (xred)
      deallocate (yred)
      deallocate (zred)
c
      deallocate (pos)
      deallocate (index)
      deallocate (indcell_loc)
      deallocate (r2vec)
!$acc wait
      return
      end
c
c    "vlistcell2" performs a complete rebuild of the
c     short range and regular vdw neighbor lists for charges 
c      using linked cells method
c
      subroutine vlistcell2
      use atmlst
      use atoms
      use cutoff
      use domdec
      use iounit
      use kvdws
      use neigh
      use vdw
      use mpi
      implicit none
      integer iglob,iloc
      integer i,ii,icell,j,j1,k,nneigloc
      integer ineig,iivdw,iv
      integer kcell,kglob,kbis
      integer ncell_loc
      integer, allocatable :: index(:),indcell_loc(:)
      real(t_p) xi,yi,zi,xk,yk,zk,r2,rdn
      real(t_p) vbufbeg2
      real(t_p), allocatable :: pos(:,:),r2vec(:)
      real(t_p), allocatable :: xred(:)
      real(t_p), allocatable :: yred(:)
      real(t_p), allocatable :: zred(:)
      logical docompute
c
      allocate (xred(nbloc))
      allocate (yred(nbloc))
      allocate (zred(nbloc))
c
      allocate (index(nbloc))
      allocate (indcell_loc(nbloc))
      allocate(pos(3,nbloc))
      allocate(r2vec(nbloc))
c
c     starting distances for long range real space interactions
c
c
c     apply reduction factors to find coordinates for each site
c
      do ii = 1, nvdwbloc
         iivdw = vdwglob(ii)
         iglob = ivdw(iivdw)
         i = loc(iglob)
         iv = ired(iglob)
         rdn = kred(iglob)
         xred(i) = rdn*(x(iglob)-x(iv)) + x(iv)
         yred(i) = rdn*(y(iglob)-y(iv)) + y(iv)
         zred(i) = rdn*(z(iglob)-z(iv)) + z(iv)
      end do
c
c     perform a complete list build
c
      do i = 1, nvdwlocnl
        iivdw = vdwglobnl(i)
        iglob  = ivdw(iivdw)
        icell = repartcell(iglob)
        iloc = loc(iglob)
c
c       align data of the local cell and the neighboring ones
c
        ncell_loc = cell_len(icell)
        indcell_loc(1:ncell_loc) = 
     $  indcell(bufbegcell(icell):(bufbegcell(icell)+cell_len(icell)-1))
        do ineig = 1, numneigcell(icell)
          kcell = neigcell(ineig,icell)
          indcell_loc(ncell_loc+1:(ncell_loc+cell_len(kcell))) = 
     $  indcell(bufbegcell(kcell):(bufbegcell(kcell)+cell_len(kcell)-1))
          ncell_loc = ncell_loc + cell_len(kcell)
        end do
c
c       do the neighbor search
c
        nneigloc = 0 
        xi = xred(iloc)
        yi = yred(iloc)
        zi = zred(iloc)
        do k = 1, ncell_loc
          kglob = indcell_loc(k)
          if (kglob.le.iglob) cycle
          if (rad(jvdw(kglob)).eq.0) cycle
          kbis = loc(kglob)
          xk = xred(kbis)
          yk = yred(kbis)
          zk = zred(kbis)
          pos(1,nneigloc+1) = xi - xk
          pos(2,nneigloc+1) = yi - yk
          pos(3,nneigloc+1) = zi - zk
          call midpointimage(xi,yi,zi,xk,yk,zk,pos(1,nneigloc+1),
     $       pos(2,nneigloc+1),pos(3,nneigloc+1),docompute)
          if (docompute) then
            nneigloc = nneigloc + 1
            index(nneigloc) = kglob
          end if
        end do
c
c       compute the distances and build the list accordingly
c
        r2vec(1:nneigloc) = pos(1,1:nneigloc)*pos(1,1:nneigloc) + 
     $      pos(2,1:nneigloc)*pos(2,1:nneigloc) + 
     $      pos(3,1:nneigloc)*pos(3,1:nneigloc)
        
        j = 0
        j1 = 0
        do k = 1, nneigloc
          r2 = r2vec(k)
          kglob = index(k)
          if (r2 .le. vshortbuf2) then
             j1 = j1 + 1
             shortvlst(j1,i) = kglob
          end if
          if (r2.le.vbuf2) then
             j = j + 1
             vlst(j,i) = kglob
          end if
        end do
        nvlst(i) = j
        nshortvlst(i) = j1
c
c     check to see if the neighbor list is too long
c
        if ((nvlst(i).ge.maxvlst).or.(nshortvlst(i).ge.maxvlst)) then
           if (rank.eq.0) then
             write (iout,10)
   10        format (/,' VBUILD  --  Too many Neighbors;',
     &                  ' Increase MAXVLST')
             call fatal
           end if
        end if
      end do
c     
      deallocate (xred)
      deallocate (yred)
      deallocate (zred)
c
      deallocate (pos)
      deallocate (index)
      deallocate (indcell_loc)
      deallocate (r2vec)
      return
      end
