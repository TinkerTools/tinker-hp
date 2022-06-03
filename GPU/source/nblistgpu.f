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
      module nblistgpu_inl
        contains
#include "image.f.inc"
#include "midpointimage.f.inc"
        subroutine ExchangeScalar(var1,var2)
        real(t_p) var1,var2
        real(t_p) temp
!$acc routine
        temp=var1;var1=var2;var2=temp;
        end subroutine

        function spbox2cmap(a,b,c,nx,ny,nxy) result(id)
        integer,intent(in )::a,b,c,nx,ny,nxy
        integer,intent(out):: id
        id = a + b*nx + c*nxy
        end function

        function spbox2smap(a,b,c,nx,ny,nxy) result(id)
        integer,intent(in )::a,b,c,nx,ny,nxy
        integer,intent(out):: id
        if (btest(c,0)) then
           id = merge(0,a + b*nx + c*nxy,btest(b,0))
        else
           id = merge(0,a + b*nx + c*nxy,btest(b,0))
        end if
        end function

        function ixyz2box(xr,yr,zr,xw,yw,zw
     &                   ,lenx_cell,leny_cell,lenz_cell
     &                   ,nx_cell,ny_cell,nz_cell,xmin,ymin,zmin)
        use cell ,only: xcell2,ycell2,zcell2,eps_cell
        implicit none
        real(t_p),intent(inout) :: xr,yr,zr
        real(t_p),intent(in)::lenx_cell,leny_cell,lenz_cell
     &           ,xmin,ymin,zmin
        integer  ,intent(in):: nx_cell,ny_cell,nz_cell
        integer  ,intent(out)::xw,yw,zw
        integer  ixyz2box
!$acc routine

        ! When to close to superior edge of cubix box
        if ((xcell2-xr).lt.eps_cell) xr = xr-0.05*lenx_cell
        if ((ycell2-yr).lt.eps_cell) yr = yr-0.05*leny_cell
        if ((zcell2-zr).lt.eps_cell) zr = zr-0.05*lenz_cell

        ! Find cell coordinate
        xw  = int((xr-xmin)/lenx_cell)
        yw  = int((yr-ymin)/leny_cell)
        zw  = int((zr-zmin)/lenz_cell)

        ! From cell coodinate to cell actual index
        ixyz2box = (xw + nx_cell*yw + nx_cell*ny_cell*zw) + 1
        end function

        function imageOutCellOcta (bx,by,bz,bxcell2,bycell2,bzcell2
     &                 ,nx_cell,ny_cell,nz_cell) result(boxid)
        implicit none
        integer,intent(in)::bxcell2,bycell2,bzcell2
        integer,intent(in)::nx_cell,ny_cell,nz_cell
        integer boxid
        integer bx,by,bz
!$acc routine

        if (bx >= bxcell2) then
           bx = bx - bxcell2
        else
           bx = bx + bxcell2
        end if
        if (by >= bycell2) then
           by = by - bycell2
        else
           by = by + bycell2
        end if
        if (bz >= bzcell2) then
           bz = bz - bzcell2
        else
           bz = bz + bzcell2
        end if
        boxid = (bx + nx_cell*by + nx_cell*ny_cell*bz) + 1
        end function
      end module
c
c     Subroutine : Reorder the neigbhor list to obtain atoms under cutoff
c     At the begining of the column
c
      subroutine reorder_nblist(nblist,nneig,nneig_cut,nloc,cutoff2,
     &                          ktype,glob)
      use atoms
      use inform ,only: deb_Path
      use nblistgpu_inl
      use utilgpu,only:rank,dir_queue

      implicit none
      integer,intent(in)   :: nloc
      integer,intent(inout):: nblist(:,:)
      integer,intent(out)  :: nneig_cut(:)
      integer,intent(in)   :: nneig(:)
      integer,intent(in)   :: glob(:),ktype(:)
      real(t_p),intent(in) :: cutoff2

      integer,parameter :: SMEM_SIZE=400
      integer i,j,k,kk,iglob,kglob,kbis
      integer write_index,read_index,read_index_s,read_cap
      integer read_pos(SMEM_SIZE)
      integer read_nei(SMEM_SIZE)
      integer read_posl(SMEM_SIZE)
      integer nneig_interact,nneig_cut_tot,nnneig
      integer neighbor
      integer count
      real(t_p) xi,yi,zi,xk,yk,zk
      logical restart

      if (deb_Path) write(*,'(2x,a)') "reorder_nblist"
      nneig_interact = 0
      write_index    = 0
!$acc parallel loop gang private(read_index,read_nei,read_pos,
!$acc&  read_posl)
!$acc&         present(nblist,nneig,nneig_cut,glob,ktype)
!$acc&         vector_length(32) async(dir_queue)
      do i = 1,nloc
         nnneig     = nneig(i)
         iglob      = ktype(glob(i))
         xi         = x(iglob)
         yi         = y(iglob)
         zi         = z(iglob)
         restart    = .true.

         ! as long as storing space is full
         ! Device Shared memory purpose
         do while (restart)
            read_index = 0
c
c           Look for interaction over the cutoff
c           And store them
c
!$acc loop vector
            do k = 1,nnneig
               kk     = nblist(k,i)
               kglob  = ktype(kk)
               xk     = x(kglob) - xi
               yk     = y(kglob) - yi
               zk     = z(kglob) - zi
               call image_inl (xk,yk,zk)
               if ((xk**2+yk**2+zk**2).gt.cutoff2) then
!$acc atomic capture
                  read_index = read_index + 1
                  read_cap   = read_index
!$acc end atomic
                  if (read_cap.gt.SMEM_SIZE) cycle  ! read_nei is full
                  !save idex and neighbor to transfer right
                  read_pos(read_cap) = k
                  read_nei(read_cap) = kk
                  nblist(k,i)        = -1
               end if
            end do
            !if (i.eq.1) then
            !   count = 0
            !   print*,nnneig
            !   do k=1,nnneig
            !      if (nblist(k,i).eq.-1) count=count+1
            !      print*, k,nblist(k,i)
            !   end do
            !   print*,'count = ',count
            !end if

            if (read_index.le.SMEM_SIZE) then 
               restart = .false.
               nneig_interact = nnneig - read_index
               read_index_s   = read_index
               nneig_cut(i)   = nneig_interact
            else
               nneig_interact = nnneig - SMEM_SIZE
               nnneig         = nnneig - SMEM_SIZE
               read_index_s   = SMEM_SIZE
            end if
            read_index   = 0
            !if (rank.eq.1) print*,i,read_index_s,nnneig
            !if (i.eq.1) then
            !   do k=1,nnneig
            !      print*, read_index, k
            !   end do
            !end if
c
c           Replace all interactions over index cutoff (read_index_s)
c
!$acc loop vector
            do k = nneig_interact+1,nneig_interact+read_index_s
               kk = nblist(k,i)
               nblist(k,i) = read_nei(k-nneig_interact)
               write_index = read_pos(k-nneig_interact)
               !save or not interaction to transfer left
               if (kk.ne.-1) then
                  read_nei(k-nneig_interact) = kk
               else
                  read_nei(k-nneig_interact) = -1
               end if
               !save naighbor index left of cutoff
               if (write_index.le.nneig_interact) then
!$acc atomic capture
                  read_index = read_index + 1
                  read_cap   = read_index
!$acc end atomic
                  read_posl(read_cap) = write_index
               end if
            end do
            if (i.eq.1.and.rank.eq.1) then
               count=0
               do k=1,nnneig
                  !print*, k,nblist(k,i)
                  if (nblist(k,i).eq.-1) count=count+1
               end do
               !print*,count
            end if

c
c           Retrive all missing interactions under the cutoff index
c
            read_index = 0
!$acc loop vector
            do k = 1,read_index_s
               kk = read_nei(k)
               if (kk.ne.-1) then
!$acc atomic capture
                  read_index = read_index + 1
                  read_cap   = read_index
!$acc end atomic
                  write_index= read_posl(read_cap)
                  nblist(write_index,i) = kk
               end if
            end do
            !if (rank.eq.1) print*, read_index

            !do k = 1, nneig(i)
            !   if (nblist(k,i).eq.-1) then
            !      print*,k,i
            !   end if
            !end do
         end do
      end do

#ifdef TINKER_DEBUG
      read_index =0
!$acc parallel loop vector_length(32) async(dir_queue)
!$acc&         present(nblist,nneig,nneig_cut) private(read_index)
      do i = 1,nloc
!$acc loop vector
         do k = 1,nneig(i)
            if (nblist(k,i).eq.-1) then
!$acc atomic update
               read_index = read_index + 1
            end if
         end do
         if (read_index.ne.0) then
            print*,'-----error-nblist-reordering-----'
            print*,i,read_index,nneig(i),nneig_cut(i)
            stop
         end if
      end do
!$acc serial
      !print*, 'Check nblist suceed'
!$acc end serial
#endif
      end subroutine
c
c     subroutine initmpipme : build the arrays to communicate direct and reciprocal fields
c     during the calculation of the induced dipoles
c
      subroutine initmpipmegpu
      use atmlst
      use domdec
      use mpole
      use pme
      use mpi
      use tinMemory
      use utilgpu ,only: rec_queue
      implicit none
      integer ierr,iipole
      integer i,iproc,tag,iglob
      integer ipre_rec1,ipre_rec2
      integer irepart,sbuf1
      integer count1,count_cap
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

c     buf1 = 0
c     buf2 = 0

      call prmem_request(buflen1,nproc,config=mhostonly)
      call prmem_request(buflen2,nproc,async=.true.)
      if (allocated(bufbeg1)) deallocate (bufbeg1)
      allocate (bufbeg1(nproc))
      if (allocated(bufbeg2)) deallocate (bufbeg2)
      allocate (bufbeg2(nproc))

!$acc data present(polerecglob,ipole,repart,buf1,buf2,
!$acc&  bufbeg1,bufbeg2,buflen2)
!$acc&     copyin(precdir_recep1,count)

!$acc parallel loop async(rec_queue)
      do i = 1,nproc
         buflen2(i) = 0
         bufbeg2(i) = 0
      end do
      buflen1 = 0
      bufbeg1 = 0

!$acc parallel loop async(rec_queue)
      do i = 1, npolerecloc
         iipole  = polerecglob(i)
         iglob   = ipole(iipole)
         irepart = repart(iglob)
         if (irepart.ne.rank) then
!$acc atomic update
            buflen2(irepart+1) = buflen2(irepart+1)+1
         end if
      end do

      count1 = 0
!$acc serial loop async(rec_queue)
      do iproc = 1, nrecdir_recep1
         ipre_rec1 = precdir_recep1(iproc)
         if (ipre_rec1.ne.rank) then
            if (buflen2(ipre_rec1+1).ne.0) then
                bufbeg2(ipre_rec1+1) = count1 + 1
            else
               bufbeg2(ipre_rec1+1) = 1
            end if
            count1 = count1 + buflen2(ipre_rec1+1)
         end if
      end do
!$acc update host(bufbeg2,buflen2) async(rec_queue)
c
!$acc parallel loop async(rec_queue)
      do i = 1, npolerecloc
         iipole  = polerecglob(i)
         iglob   = ipole(iipole)
         irepart = repart(iglob)
         if (irepart.ne.rank) then
!$acc atomic capture
            count_cap = count(irepart+1)
            count(irepart+1) = count(irepart+1) + 1
!$acc end atomic
            buf2(bufbeg2(irepart+1)+count_cap)= iipole
         end if
      end do
!$acc wait
c
c     send and receive sizes of the buffers
c
      do i = 1, nrecdir_send1
         if (precdir_send1(i).ne.rank) then
            tag = nproc*rank + precdir_send1(i) + 1
            call MPI_IRECV(buflen1(precdir_send1(i)+1),1,MPI_INT,
     $           precdir_send1(i),tag,COMM_TINKER,req(tag),ierr)
         end if
      end do
      do i = 1, nrecdir_recep1
         if (precdir_recep1(i).ne.rank) then
            tag = nproc*precdir_recep1(i) + rank + 1
            call MPI_ISEND(buflen2(precdir_recep1(i)+1),1,MPI_INT,
     $          precdir_recep1(i),tag,COMM_TINKER,req(tag),ierr)
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
!$acc update device(bufbeg1) async(rec_queue)
c
c     send and receive list of corresponding indexes
c
!$acc host_data use_device(buf1,buf2)
      do i = 1, nrecdir_send1
        if (precdir_send1(i).ne.rank) then
           tag = nproc*rank + precdir_send1(i) + 1
           call MPI_IRECV(buf1(bufbeg1(precdir_send1(i)+1)),
     $          buflen1(precdir_send1(i)+1),MPI_INT,precdir_send1(i),
     $          tag,COMM_TINKER,req2(tag),ierr)
        end if
      end do
      do i = 1, nrecdir_recep1
         if (precdir_recep1(i).ne.rank) then
            tag = nproc*precdir_recep1(i) + rank + 1
            call MPI_ISEND(buf2(bufbeg2(precdir_recep1(i)+1)),
     $           buflen2(precdir_recep1(i)+1),MPI_INT,precdir_recep1(i),
     $           tag,COMM_TINKER,req2(tag),ierr)
         end if
      end do
!$acc end host_data
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
c
!$acc end data

      !TODO Remove declare directive on thetai
      call prmem_request(thetai1,4,bsorder,nlocrec,async=.false.)
      call prmem_request(thetai2,4,bsorder,nlocrec,async=.false.)
      call prmem_request(thetai3,4,bsorder,nlocrec,async=.false.)

      if (.not.associated(thetai1_p)) then
c        sbuf1 = merge(nlocrec+int(real(nlocrec,t_p)*mem_inc),nlocrec
c    &                ,extra_alloc)
c        allocate (thetai1(4,bsorder,sbuf1))
c        allocate (thetai2(4,bsorder,sbuf1))
c        allocate (thetai3(4,bsorder,sbuf1))
c        s_prmem =  s_prmem + 12*bsorder*sbuf1*szoTp
c       sd_prmem = sd_prmem + 12*bsorder*sbuf1*szoTp
         call AssociateThetai_p
#ifdef _OPENACC
         call attach_pmecu_pointer(1)
#endif
!$acc parallel loop async default(present)
         do i =1,size(thetai1)
            thetai1(i,1,1) =0
            thetai2(i,1,1) =0
            thetai3(i,1,1) =0
         end do
      else if (nlocrec>size(thetai1,dim=3)) then
c        sbuf1 = merge(nlocrec+int(real(nlocrec,t_p)*mem_inc),nlocrec
c    &                ,extra_alloc)
c        tag = size(thetai1,dim=3)
c        s_prmem =  s_prmem - 12*bsorder*tag*szoTp
c       sd_prmem = sd_prmem - 12*bsorder*tag*szoTp
c        deallocate (thetai1)
c        deallocate (thetai2)
c        deallocate (thetai3)
c        allocate (thetai1(4,bsorder,sbuf1))
c        allocate (thetai2(4,bsorder,sbuf1))
c        allocate (thetai3(4,bsorder,sbuf1))
c        s_prmem =  s_prmem + 12*bsorder*sbuf1*szoTp
c       sd_prmem = sd_prmem + 12*bsorder*sbuf1*szoTp
         call AssociateThetai_p
#ifdef _OPENACC
         call attach_pmecu_pointer(1)
#endif
      end if

      deallocate (req)
      deallocate (req2)
      deallocate (count)
      end
c
c     check Atom Inboxing algorithm
c
      subroutine check_atoms_inboxing(cell_len,ncell_tot,nlocnl)
      use domdec ,only: rank
      use sizes  ,only: tinkerdebug
      use utilgpu,only: rec_queue
      implicit none
      integer ncell_tot,nlocnl
      integer,intent(in)::cell_len(ncell_tot)
      integer i, k

      k = 0
!$acc parallel loop default(present) async(rec_queue)
      do i = 1,ncell_tot
         k = k + cell_len(i)
      end do
!$acc wait
      if (k.ne.nlocnl) then
32       format( ' build_cell_listgpu : Atoms are missing from',
     &   ' inboxing procedure  (rank ',I5,')',
     &   /,' total count',I10,5x,'nlocnl',I10)
         print 32, rank,k,nlocnl
      end if
      end subroutine
c
c
c     subroutine build_cell_list : build the cells in order to build the non bonded neighbor
c     lists with the cell-list method
c
      subroutine build_cell_listgpu(istep)
      use atoms
      use boxes
      use bound
      use cell
      use domdec
      use inform    ,only: deb_Path
      use nblistgpu_inl
      use neigh
      use mpi
      use tinheader
      use tinMemory ,only:prmem_request
      use utilgpu   ,only: nSMP, cores_SMP
     &              , openacc_abort,rec_queue
      use utils     ,only: set_to_zero1_int8
      use sizes     ,only: tinkerdebug
      implicit none
      integer i,proc,icell,j,k,p,q,r,istep,iglob,ii,kk
      integer count,iloc
      integer,parameter:: ncell2buff=2
      integer ncell2buffb
      integer num_neigcell
      integer Nblock,Nrest,kstart,kend,buffersize,bufloc
      integer temp_x,temp_y,temp_z,locat
      integer temi_x,temi_y,temi_z
      integer nx_cell2,ny_cell2,nz_cell2
      integer xw,yw,zw
      integer box1,box2,box2i,setb,setbi
      integer icell_len,imageCell,distImage2Cell
      integer temp,numneig,numneig_cap,tempcell
      real(t_p) xmin,xmax,ymin,ymax,zmin,zmax
      real(t_p) lenx,leny,lenz
      real(t_p) mbuf,vbuf,bigbuf
      real(t_p) lenx_cell,leny_cell,lenz_cell
      real(t_p) eps_x,eps_y,eps_z
      real(t_p) xr,yr,zr
      real(t_p) cellXmin,cellXmax,cellYmin,cellYmax,cellZmin,cellZmax
     &         ,cellXin,cellYin,cellZin
      logical   close2Zbox,far2Zbox

      real(t_p), allocatable :: xbegcelltemp(:),ybegcelltemp(:)
      real(t_p), allocatable :: zbegcelltemp(:)
      real(t_p), allocatable :: xendcelltemp(:),yendcelltemp(:)
      real(t_p), allocatable :: zendcelltemp(:)
      integer, allocatable :: filledcell(:,:),indcelltemp(:)
      logical docompute
!$acc routine(fatal_acc)
c
      if (deb_Path) write(*,'(2x,a)') 'build_cell_listgpu'

      mbuf = sqrt(mbuf2)
      vbuf = sqrt(vbuf2)
      !FIXME over buffer on vbuf when done

      ncell2buffb = ncell2buff
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

 20   continue
      if ((mbuf).gt.vbuf) then
        bigbuf = (mbuf)/ncell2buffb
      else
        bigbuf = vbuf/ncell2buffb
      end if
c
      lenx       = abs(xmax-xmin)
      nx_cell    = max(2*ncell2buffb+1,int(lenx/(bigbuf)))
      lenx_cell  = lenx/nx_cell
      leny       = abs(ymax-ymin)
      ny_cell    = max(2*ncell2buffb+1,int(leny/(bigbuf)))
      leny_cell  = leny/ny_cell
      lenz       = abs(zmax-zmin)
      nz_cell    = max(2*ncell2buffb+1,int(lenz/(bigbuf)))
      lenz_cell  = lenz/nz_cell

      ! Force even cell number along 3 dimensions
      if (octahedron) then
         if (nx_cell.ne.ny_cell) then
24          format ('build_cell_listgpu -- Improper cell number',/,
     &              ' nx',I5,' ny',I5,' nz',I5)
            print 24, nx_cell,ny_cell,nz_cell
            call fatal
         end if
         if (btest(nx_cell,0)) then
            nx_cell     = nx_cell-1
            if (nx_cell <= 2*ncell2buffb+1 ) then
               ncell2buffb = ncell2buffb + 1
               goto 20
            end if
            ny_cell     = nx_cell
            nz_cell     = nx_cell
            lenx_cell   = lenx/nx_cell
            leny_cell   = leny/ny_cell
            lenz_cell   = lenz/nz_cell
         end if
         nx_cell2  = nx_cell/2
         ny_cell2  = ny_cell/2
         nz_cell2  = nz_cell/2
      end if

      ncell_tot  = nx_cell*ny_cell*nz_cell
      eps_x      = lenx*prec_eps
      eps_y      = leny*prec_eps
      eps_z      = lenz*prec_eps

      ! compute max neighbor cell
      if (octahedron) then
         num_neigcell = 2*((2*ncell2buffb+1)**3)
      else
         num_neigcell = (2*ncell2buffb+1)**3-1
      end if
      call prmem_request (numneigcell, ncell_tot,async=.true.)
      call prmem_request(neigcell,num_neigcell,ncell_tot
     &                         ,async=.true.)

c
c     Assign atoms to cells
c
      if (.not.allocated(repartcell)) then
         allocate(repartcell(n))
!$acc enter data create(repartcell) async(rec_queue)
      end if
      call prmem_request(bufbegcell,ncell_tot,async=.true.)
      call prmem_request(cell_len,ncell_tot,async=.true.)

      allocate (indcelltemp(n))
!$acc data create(indcelltemp)
!$acc&     present(cell_len,indcell,bufbegcell,repartcell,
!$acc&  use_bounds) async(rec_queue)

!$acc parallel loop async(rec_queue)
      do i = 1,ncell_tot
         cell_len(i) = 0
      end do
!$acc parallel loop async(rec_queue)
      do i = 1,n
         indcelltemp(i) = 0
      end do

!$acc parallel loop async(rec_queue)
      do i = 1,nlocnl
         iglob = ineignl(i)
         xr    = x(iglob)
         yr    = y(iglob)
         zr    = z(iglob)
         if (use_bounds) call image_inl(xr,yr,zr)
         if ((xcell2-xr).lt.eps_cell) xr = xr-0.05*lenx_cell
         if ((ycell2-yr).lt.eps_cell) yr = yr-0.05*leny_cell
         if ((zcell2-zr).lt.eps_cell) zr = zr-0.05*lenz_cell
         xw  = int((xr-xmin)/lenx_cell)
         yw  = int((yr-ymin)/leny_cell)
         zw  = int((zr-zmin)/lenz_cell)
         icell = (xw + nx_cell*yw + nx_cell*ny_cell*zw) + 1
         repartcell (iglob) = icell
!$acc atomic capture
         cell_len (icell) = cell_len(icell)+1
         icell_len        = cell_len(icell)
!$acc end atomic
         indcelltemp (iglob) = icell_len
#ifdef TINKER_DEBUG
            if (xw<0.or.xw>=nx_cell) print*,'xb out',xw,xr,
     &         real(0.05*lenx_cell,4),eps_cell
            if (yw<0.or.yw>=ny_cell) print*,'yb out',xw,yr,
     &         real(0.05*leny_cell,4),eps_cell
            if (zw<0.or.zw>=nz_cell) print*,'zb out',xw,zr,
     &         real(0.05*lenz_cell,4),eps_cell
#endif
      end do
!$acc update host(cell_len) async(rec_queue)

      if (tinkerdebug.gt.0)
     &   call check_atoms_inboxing(cell_len,ncell_tot,nlocnl)
c
c     Find neighbor cells for octahedron
c
      if (octahedron) then

c
c     Build Octahedron Adjacency matrix on cells
c
      szMatb  = ((ncell_tot-1)/bit_si) + 1
      nbMatb  = max(1,ncell_tot)
      call prmem_request (matb_lst,szMatb*nbMatb,async=.true.)
      call set_to_zero1_int8(matb_lst(1),int(szMatb,8)*nbMatb,
     &     rec_queue)

      ! loop on every cell
!$acc parallel loop async(rec_queue) private(locat)
!$acc&         vector_length(32)
!$acc&         present(matb_lst,neigcell,numneigcell)
      do ii = 1,ncell_tot
         k  = (ii-1)/(nx_cell*ny_cell)
         j  = (ii-1 - ((ii-1)/(nx_cell*ny_cell))*(nx_cell*ny_cell))
     &        /nx_cell
         i  = mod(ii-1,nx_cell)
         locat = 0
         numneigcell(ii) = 0
         if (cell_len(ii).eq.0) then
            cycle
         end if

!$acc loop vector collapse(3)
         do p = -ncell2buffb,ncell2buffb
            do q = -ncell2buffb,ncell2buffb
               do r = -ncell2buffb,ncell2buffb
                  temp_x = i + r
                  temp_y = j + q
                  temp_z = k + p


                  ! Find cell image in cubic box
                  if      (temp_x.lt.0) then
                     temp_x = temp_x + nx_cell
                  else if (temp_x.ge.nx_cell) then
                     temp_x = temp_x - nx_cell
                  end if
                  if      (temp_y.lt.0) then
                     temp_y = temp_y + ny_cell
                  else if (temp_y.ge.ny_cell) then
                     temp_y = temp_y - ny_cell
                  end if
                  if      (temp_z.lt.0) then
                     temp_z = temp_z + nz_cell
                  else if (temp_z.ge.nz_cell) then
                     temp_z = temp_z - nz_cell
                  end if
                  tempcell = (temp_z*ny_cell*nx_cell+temp_y*nx_cell+
     &                        temp_x) + 1

                  ! Where is neighbor cell located compare to octahedron
                  cellXmin = xmin +  temp_x   *lenx_cell
                  cellXmax = xmin + (temp_x+1)*lenx_cell
                  if (abs(cellXmin)>abs(cellXmax))
     &               call ExchangeScalar(cellXmin,cellXmax)
                  cellYmin = ymin +  temp_y   *leny_cell
                  cellYmax = ymin + (temp_y+1)*leny_cell
                  if (abs(cellYmin)>abs(cellYmax))
     &               call ExchangeScalar(cellYmin,cellYmax)
                  cellZmin = zmin +  temp_z   *lenz_cell
                  cellZmax = zmin + (temp_z+1)*lenz_cell
                  if (abs(cellZmin)>abs(cellZmax))
     &               call ExchangeScalar(cellZmin,cellZmax)

                  close2Zbox =(abs(cellXmin)+abs(cellYmin)+abs(cellZmin)
     &                         .le.box34)
                  far2Zbox   =(abs(cellXmax)+abs(cellYmax)+abs(cellZmax)
     &                         .le.box34)

                  ! Find (tempcell,ii) location in matb_lst
                  box1 = (ii-1)*szMatb
                  box2 = box1 + ishft(tempcell-1,-bit_sh) + 1
                  setb = ishft( 1, iand(tempcell-1,bit_si-1) )

                  ! Find neighbor cell attributes
                  if (close2Zbox.and.far2Zbox) then !Inside box
                     ! Avoid ii cell
                     if (p==0.and.q==0.and.r==0) cycle
                     !if (cell_len(tempcell).eq.0) cycle !might be useful with MPI
!$acc atomic
                     matb_lst(box2) = ior( matb_lst(box2),setb )
                  else if (close2Zbox.and..not.far2Zbox) then !Partially inside box
                     imageCell = imageOutCellOcta(temp_x,temp_y,temp_z
     &                              ,nx_cell2,ny_cell2,nz_cell2,nx_cell
     &                              ,ny_cell,nz_cell)
                     ! Max Distance between imageCell and cell i
                     distImage2Cell = max(abs(temp_x-i),
     &                                max(abs(temp_y-j),abs(temp_z-k)))

                     ! Find (imageCell,ii) location in matb_lst
                     box2i = box1 + ishft(imageCell-1,-bit_sh) + 1
                     setbi = ishft( 1, iand(imageCell-1,bit_si-1) )

                     ! Avoid ii cell, only add his image
                     if (p==0.and.q==0.and.r==0.and.
     &                   distImage2Cell>ncell2buffb) then
!$acc atomic
                        matb_lst(box2i) = ior( matb_lst(box2i),setbi )
                     else if (distImage2Cell>ncell2buffb) then
c!$acc atomic capture
c                        locat = locat+2
c                        count = locat
c!$acc end atomic
c                        !if (cell_len(tempcell).eq.0) cycle !!! Pay attention to 'locat'
c                        neigcell(count-1,ii) = tempcell
c                        neigcell(count  ,ii) = imageCell
!$acc atomic
                        matb_lst(box2 ) = ior( matb_lst(box2 ),setb  )
!$acc atomic
                        matb_lst(box2i) = ior( matb_lst(box2i),setbi )
                     else
                        !if (cell_len(tempcell).eq.0) cycle
!$acc atomic
                        matb_lst(box2 ) = ior( matb_lst(box2 ),setb  )
                     end if
                  else if (.not.close2Zbox.and..not.far2Zbox) then !Outside box
                     !if (cell_len(imageCell).eq.0) cycle !Might be useful with MPI
                     imageCell = imageOutCellOcta(temp_x,temp_y,temp_z
     &                              ,nx_cell2,ny_cell2,nz_cell2,nx_cell
     &                              ,ny_cell,nz_cell)
                     ! Max Distance between imageCell and cell i
                     distImage2Cell = max(abs(temp_x-i),
     &                                max(abs(temp_y-j),abs(temp_z-k)))

                     ! Find (imageCell,ii) location in matb_lst
                     box2i = box1 + ishft(imageCell-1,-bit_sh) + 1
                     setbi = ishft( 1, iand(imageCell-1,bit_si-1) )

                     if (cell_len(imageCell).eq.0.or.
     &                   distImage2Cell<ncell2buffb+1) cycle
!$acc atomic
                     matb_lst(box2i) = ior( matb_lst(box2i),setbi )
                  else           ! Absurd
                     print*,':build_cell_listgpu'
                     print*,tempcell,ii,close2Zbox,far2Zbox
                     print*,cellXmin,cellYmin,cellZmin
                     print*,cellXmax,cellYmax,cellZmax
                     print*,'impossible case in neighbor cell'
                     call fatal_acc
                  end if
               end do
            end do
         end do
         if (locat>num_neigcell) then
            print*,'build_cell_listgpu: max neighbor cell reached'
            call fatal_acc
         end if
         !numneigcell(ii) = locat
      end do
c
c     From Adjacency matrix to neighbor cells list
c
!$acc parallel loop gang vector collapse(2) async(rec_queue)
!$acc&         present(matb_lst,neigcell,numneigcell) private(locat)
      do i = 1, ncell_tot
         do j = 0, ncell_tot-1
            box2 = ishft( j,-bit_sh ) + 1
            setb = iand ( j,bit_si-1 )
            if (btest( matb_lst((i-1)*szMatb+box2),setb )) then
!$acc atomic capture
               numneigcell(i) = numneigcell(i) + 1
               count = numneigcell(i)
!$acc end atomic
               neigcell(count,i) = j+1
            end if
         end do
      end do
c     print*, 'num neighbor box',nx_cell,ny_cell,nz_cell,
c    &        num_neigcell
c     do i = 1,ncell_tot
c        if (cell_len(i).eq.0) cycle
c        print*, i, cell_len(i), numneigcell(i)
c     end do
c     print*, '---------------------------------'

      else
      !TODO Integrate ii cell to neigcell and update [vcm]list2?gpu
c
c     Find neighbor cells for cubic box
c
!$acc parallel loop async(rec_queue) private(locat)
!$acc&         vector_length(32)
!$acc&         present(neigcell,numneigcell)
      do ii = 1, ncell_tot
         k  = (ii-1)/(nx_cell*ny_cell)
         j  = (ii-1 - ((ii-1)/(nx_cell*ny_cell))*(nx_cell*ny_cell))
     &        /nx_cell
         i  = mod(ii-1,nx_cell)
         locat = 0
!$acc loop vector collapse(3)
         do p = -ncell2buff,ncell2buff
            do q = -ncell2buff,ncell2buff
               do r = -ncell2buff,ncell2buff
                  temp_x = i + r
                  temp_y = j + q
                  temp_z = k + p

                  ! Avoid ii cell
                  if (p==0.and.q==0.and.r==0) cycle

                  ! Find cell image in cubic box
                  if      (temp_x.lt.0) then
                     temp_x = temp_x + nx_cell
                  else if (temp_x.ge.nx_cell) then
                     temp_x = temp_x - nx_cell
                  end if
                  if      (temp_y.lt.0) then
                     temp_y = temp_y + ny_cell
                  else if (temp_y.ge.ny_cell) then
                     temp_y = temp_y - ny_cell
                  end if
                  if      (temp_z.lt.0) then
                     temp_z = temp_z + nz_cell
                  else if (temp_z.ge.nz_cell) then
                     temp_z = temp_z - nz_cell
                  end if
                  tempcell = (temp_z*ny_cell*nx_cell+temp_y*nx_cell+
     &                        temp_x) + 1
c                 locat =((r+ncell2buff) 
c    &                   +(q+ncell2buff)*(2*ncell2buff+1)
c    &                   +(p+ncell2buff)*(2*ncell2buff+1)**2) + 1
!$acc atomic capture
                  locat = locat +1
                  count = locat
!$acc end atomic
                  neigcell (count,ii) = tempcell
               end do
            end do
         end do
         numneigcell(ii) = num_neigcell
      end do

      end if
c
c     Scan cell_len and build index for repartcell
c
      if (.not.allocated(indcell)) then
         allocate(indcell(n))
!$acc enter data create(indcell) async(rec_queue)
      end if
c
!$acc wait
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
!$acc update device(bufbegcell) async(rec_queue)
c
!$acc parallel loop async(rec_queue)
      do i = 1, nlocnl
        iglob = ineignl(i)
        icell = repartcell(iglob)
        iloc  = bufbegcell(icell) + indcelltemp(iglob) - 1
        indcell(iloc) = iglob
      end do

!$acc end data
      deallocate (indcelltemp)
      end
c
c     subroutine build_cell_list2 : build the cells small enough in order to build
c     the non bonded neighbor lists with the cell-list method
c
      subroutine build_cell_list2(ineignl,cell_order,
     &           nlocnl,buf,bP)
      use atoms ,only: x,y,z,n
      use bound ,only: use_bounds
      use boxes ,only: orthogonal,octahedron,xbox,ybox,zbox
      use cell  ,only: xcell2,ycell2,zcell2,eps_cell
      use domdec,only: nproc,rank,COMM_TINKER,comm_dir,comm_rec,ndir
      use inform,only: deb_Path,tindPath
      use nblistgpu_inl
      use neigh ,only: repartcell,cell_scan,cell_len,cell_len1,boxPart
     &          ,nx_cell,ny_cell,nz_cell,ncell_tot,max_cell_len
      use mpi
      use potent,only: use_pmecore
#ifdef _OPENACC
      use thrust
#endif
      use tinheader,only: ti_eps,prec_eps
      use tinMemory,only: prmem_request
      use sizes    ,only: tinkerdebug
      use utils    ,only: set_to_zero1_int
      use utilgpu  ,only: nSMP,cores_SMP
     &             ,rec_queue,BLOCK_SIZE
#ifdef _OPENACC
     &             ,rec_stream
#endif
      implicit none

      integer  ,intent(in)   :: nlocnl
      integer  ,intent(in)   :: ineignl(:)
      integer  ,intent(inout):: cell_order(:)
      real(t_p),intent(in)   :: buf
      type(boxPart),intent(inout)::bP

      integer i,j,k,proc,icell,istep,iglob
      integer bx_cell,by_cell,bz_cell
      integer ierr
      integer xw,yw,zw,miter
      integer icell_len
      integer commloc
      integer levelUp,incr,iter
      integer,save:: save_ncell
      integer min_cell_len
      real(t_p) xmin,xmax,ymin,ymax,zmin,zmax
      real(t_p) lenx,leny,lenz,xl
      real(t_p) mbuf,vbuf
      real(t_p) lenx_cell,leny_cell,lenz_cell
      real(t_p) eps_x,eps_y,eps_z
      real(t_p) xr,yr,zr
      logical xb,yb,zb
      logical first_loop
      integer,save:: first_in=0

      if (use_pmecore) then  !Deal with pme_core
         if (rank.lt.ndir) then
            commloc = comm_dir
         else
            commloc = comm_rec
         end if
      else
         commloc = COMM_TINKER
      end if

      if (deb_Path) write(*,'(2x,a)') 'build_cell_listgpu2'
      first_in=first_in+1

      ! Fetch Partition information
      iter    = bp%ipart
      bx_cell = bp%bx_c
      by_cell = bp%by_c
      bz_cell = bp%bz_c
      nx_cell = bp%nx_c
      ny_cell = bp%ny_c
      nz_cell = bp%nz_c
      xb      = bp%nx_l
      yb      = bp%ny_l
      zb      = bp%nz_l

      first_loop   = .true.
      incr  = merge(2,1,octahedron)
      max_cell_len = BLOCK_SIZE+1
      min_cell_len = huge(min_cell_len)
      xmin  = -xcell2
      xmax  =  xcell2
      ymin  = -ycell2
      ymax  =  ycell2
      zmin  = -zcell2
      zmax  =  zcell2
      lenx  = abs(xmax-xmin)
      leny  = abs(ymax-ymin)
      lenz  = abs(zmax-zmin)

      levelUp = 1

      if (first_in.eq.1) then
         xl = real((((xbox*ybox*zbox)*BLOCK_SIZE)/n)**(1.0d0/3),t_p)
         if(deb_Path)print'(A,2F10.6)','box axe 3diag ',xl,xl*sqrt(3.0)
c        nx_cell = max(1,int(lenx/xl))
c        ny_cell = max(1,int(leny/xl)+1)
c        nz_cell = max(1,int(lenz/xl)+1)
         save_ncell = 0
c        nx_cell = max(1,int(floor(lenx/7)))
c        ny_cell = max(1,int(floor(leny/7)))
c        nz_cell = max(1,int(floor(lenz/7)))
      end if

      call prmem_request(repartcell,n)

c
c     Put Atoms into small boxes
c
      do while (max_cell_len .gt. BLOCK_SIZE)
         if (.not.first_loop) then
            miter = mod(iter,4)
            if (abs(max_cell_len-BLOCK_SIZE).lt.4) then
               nx_cell = nx_cell + incr
               xb = .false.
            else if (iter.le.2) then
               bx_cell = bx_cell + 1
               xb = .true.
            else if (iter.le.5) then
               nz_cell = nz_cell + incr
               ny_cell = ny_cell + incr
               yb = .false.
               zb = .false.
            else if (miter.eq.0.and.iter.eq.6.and.bx_cell.eq.1) then
               bz_cell = bz_cell + 1
               by_cell = by_cell + 1
               yb = .true.
               zb = .true.
            else if (miter.gt.0) then
               i = merge(2,incr,levelUp.eq.1); levelUp=i;
               nx_cell = nx_cell + levelUp
               xb = .false.
            else if (miter.eq.0) then
               nz_cell = nz_cell + incr
               ny_cell = ny_cell + incr
               yb = .false.
               zb = .false.
            end if
c           nx_cell = nx_cell + 1
         end if

         iter       = iter + 1
         first_loop = .false.

         if (xb) then
            nx_cell = max(1,int(bx_cell*lenx/buf))
         else
            bx_cell = ceiling(buf*nx_cell/lenx)
         endif
         if (yb) then
            ny_cell = max(1,int(by_cell*leny/buf))
         else
            by_cell = ceiling(buf*ny_cell/leny)
         endif
         if (zb) then
            nz_cell = max(1,int(bz_cell*lenz/buf))
         else
            bz_cell = ceiling(buf*nz_cell/lenz)
         endif

         if (octahedron) then
            ! Force even cell number along 3 dimensions
            if (btest(nx_cell,0)) nx_cell = nx_cell-1
            if (btest(ny_cell,0)) ny_cell = ny_cell-1
            if (btest(nz_cell,0)) nz_cell = nz_cell-1
            !nx_cell2  = nx_cell/2
            !ny_cell2  = ny_cell/2
            !nz_cell2  = nz_cell/2
         end if

         lenx_cell  = lenx/nx_cell
         leny_cell  = leny/ny_cell
         lenz_cell  = lenz/nz_cell
         ncell_tot  = nx_cell*ny_cell*nz_cell
         save_ncell = max(ncell_tot,save_ncell)
         max_cell_len = 0
         min_cell_len = huge(min_cell_len)

         call prmem_request(cell_len ,save_ncell,async=.false.)
         call prmem_request(cell_len1,save_ncell,async=.false.)
         call prmem_request(cell_scan,save_ncell,async=.false.)

         call set_to_zero1_int(cell_len,size(cell_len),rec_queue)
!$acc parallel loop async(rec_queue)
!$acc&         default(present)
         do i = 1,nlocnl
            iglob = ineignl(i)
            xr    = x(iglob)
            yr    = y(iglob)
            zr    = z(iglob)
            if (use_bounds) call image_inl(xr,yr,zr)
            if ((xcell2-xr).lt.eps_cell) xr = xr-0.05*lenx_cell
            if ((ycell2-yr).lt.eps_cell) yr = yr-0.05*leny_cell
            if ((zcell2-zr).lt.eps_cell) zr = zr-0.05*lenz_cell
            xw    = int((xr-xmin)/lenx_cell)
            yw    = int((yr-ymin)/leny_cell)
            zw    = int((zr-zmin)/lenz_cell)
#ifdef TINKER_DEBUG
            if (xw<0.or.xw>=nx_cell) print*,'xb out',xw,xr,
     &         real(0.05*lenx_cell,4),eps_cell
            if (yw<0.or.yw>=ny_cell) print*,'yb out',xw,yr,
     &         real(0.05*leny_cell,4),eps_cell
            if (zw<0.or.zw>=nz_cell) print*,'zb out',xw,zr,
     &         real(0.05*lenz_cell,4),eps_cell
#endif
            icell = (xw + nx_cell*yw + nx_cell*ny_cell*zw) + 1
            repartcell (i)   = icell
!$acc atomic capture
            cell_len (icell) = cell_len(icell)+1
            icell_len        = cell_len(icell)
!$acc end atomic
            max_cell_len     = max(max_cell_len,icell_len)
         end do

         if (btest(tinkerdebug,tindPath)) then
!$acc parallel loop async(rec_queue) present(cell_len)
            do i = 1,ncell_tot
               min_cell_len     = min(min_cell_len,cell_len(i))
            end do
         end if
!$acc wait
         if(tinkerdebug.gt.0)
     &      call check_atoms_inboxing(cell_len,ncell_tot,nlocnl)

         if (nproc.ne.1) then
            call MPI_AllReduce(MPI_IN_PLACE,max_cell_len,1,MPI_INT,
     &           MPI_MAX,commloc,ierr)
            if (btest(tinkerdebug,tindPath))
     &      call MPI_AllReduce(MPI_IN_PLACE,min_cell_len,1,MPI_INT,
     &           MPI_MIN,commloc,ierr)
         end if
 20      format(a,7I5,2X,I8,2x,3I3)
         if (deb_Path)
     &   print 20,'box iter ',iter,nx_cell,ny_cell,nz_cell,max_cell_len
     &                       ,nlocnl/ncell_tot,min_cell_len,ncell_tot
     &                       ,bx_cell,by_cell,bz_cell
      end do

      ! Save box partition info
      bp%ipart = iter-1
      bp%bx_c  = bx_cell
      bp%by_c  = by_cell
      bp%bz_c  = bz_cell
      bp%nx_c  = nx_cell
      bp%ny_c  = ny_cell
      bp%nz_c  = nz_cell
      bp%nx_l  = xb
      bp%ny_l  = yb
      bp%nz_l  = zb


#ifdef _OPENACC
!$acc parallel loop async present(cell_order)
      do i = 1,nlocnl
         cell_order(i) = i
      end do
!$acc wait
!$acc host_data use_device(repartcell,cell_order,cell_len,cell_scan)
      call thrust_stable_sort_by_key(repartcell,nlocnl,cell_order
     &                              ,rec_stream)
      call thrust_exclusive_scan(cell_len,ncell_tot,cell_scan
     &                              ,rec_stream)
!$acc end host_data
#else
      call sort3 (nlocnl,repartcell,cell_order)
      cell_scan(1) = 0
      do i = 2,ncell_tot
         cell_scan(i) = cell_len(i-1)+cell_scan(i-1)
      end do
#endif

!$acc parallel loop present(cell_len1,cell_len) async(rec_queue)
      do i = 1,ncell_tot
         !cell_len1(i) = int(cell_len(i),1)  !unsupported by pgi
         cell_len1(i) = cell_len(i)
      end do

      end subroutine

      ! Change Interface module if modify calling arguments
      subroutine set_pairlist_Cellorder(atList,a_plist,buildnl)
      use atoms  ,only: n,x,y,z
      use domdec
      use neigh
      use utilgpu
      use mpole
      implicit none
      integer atList(*)
      type(pair_atlst),target:: a_plist
      logical,intent(in):: buildnl

      integer i,iglob,iipole,iploc
      integer natmnlb,natmnl
      integer  ,pointer::c_glob(:),c_key(:) 
      real(t_p),pointer::cell_x(:),cell_y(:),cell_z(:)

      ! Fetch a_plist attributes
      natmnlb =  a_plist%natmnlb
      natmnl  =  a_plist%natmnl 

      if (buildnl) then

         call prmem_request(a_plist%c_glob,natmnlb,async=.true.)
         call prmem_request(a_plist%cell_x,natmnlb,async=.true.)
         call prmem_request(a_plist%cell_y,natmnlb,async=.true.)
         call prmem_request(a_plist%cell_z,natmnlb,async=.true.)
         c_glob  => a_plist%c_glob
         c_key   => a_plist%c_key
         cell_x  => a_plist%cell_x
         cell_y  => a_plist%cell_y
         cell_z  => a_plist%cell_z

!$acc parallel loop async(rec_queue)
!$acc&         present(atList,x,y,z
!$acc&   ,c_glob,c_key,cell_x,cell_y,cell_z)
         do i = 1, natmnlb
            if (i.le.natmnl) then
               iglob     = atList(c_key(i))
               c_glob(i) = iglob
               cell_x(i) = x(iglob)
               cell_y(i) = y(iglob)
               cell_z(i) = z(iglob)
            else
               c_glob(i) = n+1
               cell_x(i) = inf
               cell_y(i) = inf
               cell_z(i) = inf
            end if
         end do

      else

         c_glob  => a_plist%c_glob
         cell_x  => a_plist%cell_x
         cell_y  => a_plist%cell_y
         cell_z  => a_plist%cell_z
!$acc parallel loop async(rec_queue)
!$acc&         present(c_glob
!$acc&   ,celle_x,celle_y,celle_z,x,y,z)
         do i = 1, natmnlb
            if (i.le.natmnl) then
               iglob     = c_glob(i)
               cell_x(i) = x(iglob)
               cell_y(i) = y(iglob)
               cell_z(i) = z(iglob)
            else
               cell_x(i) = inf
               cell_y(i) = inf
               cell_z(i) = inf
            end if
         end do

      end if
      end subroutine

      ! Change Interface module if modify calling arguments
      subroutine set_ElecData_CellOrder(rebuild_nl)
      use atoms  ,only: n,x,y,z
      use atmlst ,only: poleglobnl
      use domdec
      use neigh
      use utilgpu
      use mpole
      implicit none
      logical,intent(in):: rebuild_nl
      integer i,iglob,iipole,iploc

      if (rebuild_nl) then

         call prmem_request(celle_glob,npolelocnlb,async=.true.)
         call prmem_request(celle_pole,npolelocnlb,async=.true.)
         call prmem_request(celle_loc ,npolelocnlb,async=.true.)
         call prmem_request(celle_ploc,npolelocnlb,async=.true.)
         call prmem_request(celle_plocnl,npolelocnlb,async=.true.)
         call prmem_request(celle_x   ,npolelocnlb,async=.true.)
         call prmem_request(celle_y   ,npolelocnlb,async=.true.)
         call prmem_request(celle_z   ,npolelocnlb,async=.true.)

!$acc parallel loop async(rec_queue)
!$acc&         present(celle_glob,celle_pole,celle_key,
!$acc&   poleglobnl,celle_x,celle_y,celle_z,celle_plocnl,polelocnl)
         do i = 1, npolelocnlb
            if (i.le.npolelocnl) then
               iipole        = poleglobnl(celle_key(i))
               iglob         = ipole(iipole)
               celle_pole(i) = iipole
               celle_plocnl(i) = polelocnl(iipole)
               celle_glob(i) = iglob
               celle_x   (i) = x(iglob)
               celle_y   (i) = y(iglob)
               celle_z   (i) = z(iglob)
            else
               celle_pole(i) = n
               celle_glob(i) = n
               celle_plocnl(i)= npolelocnl
               celle_x   (i) = inf
               celle_y   (i) = inf
               celle_z   (i) = inf
            end if
         end do

      else

!$acc parallel loop async(rec_queue)
!$acc&         present(celle_glob,celle_pole,celle_key,
!$acc&   celle_loc,celle_ploc,celle_x,celle_y,celle_z,x,y,z,
!$acc&   loc,poleloc)
         do i = 1, npolelocnlb
            if (i.le.npolelocnl) then
               iglob        = celle_glob(i)
               celle_loc(i) = loc(iglob)
               celle_ploc(i)= poleloc(celle_pole(i))
               celle_x  (i) = x(iglob)
               celle_y  (i) = y(iglob)
               celle_z  (i) = z(iglob)
            else
               celle_loc(i) = nbloc
               celle_ploc(i)= npolebloc
               celle_x  (i) = inf
               celle_y  (i) = inf
               celle_z  (i) = inf
            end if
         end do

      end if
      end subroutine

      subroutine set_ChgData_CellOrder(rebuild_nl)
      use atoms  ,only: n,x,y,z
      use atmlst ,only: chgglobnl
      use charge
      use domdec
      use kvdws  ,only: rad,eps,radv,epsv
      use neigh
      use potent ,only: fuse_chglj
      use utilgpu
      use vdw
      implicit none
      logical,intent(in):: rebuild_nl
      integer i,iglob,iichg

      if (rebuild_nl) then

         call prmem_request(celle_glob,nionlocnlb,async=.true.)
         call prmem_request(celle_chg ,nionlocnlb,async=.true.)
         call prmem_request(celle_loc ,nionlocnlb,async=.true.)
         call prmem_request(celle_ploc,nionlocnlb,async=.true.)
         call prmem_request(celle_plocnl,nionlocnlb,async=.true.)
         call prmem_request(celle_x   ,nionlocnlb,async=.true.)
         call prmem_request(celle_y   ,nionlocnlb,async=.true.)
         call prmem_request(celle_z   ,nionlocnlb,async=.true.)
         if (fuse_chglj) then
            call prmem_request(cellv_jvdw,nionlocnlb,async=.true.)
            call prmem_request(radv   ,nionlocnlb,async=.true.)
            call prmem_request(epsv   ,nionlocnlb,async=.true.)
            nvdwlocnlb = nionlocnlb
         end if

         if (fuse_chglj) then
!$acc parallel loop async(rec_queue)
!$acc&         present(celle_glob,celle_chg,celle_key,
!$acc&   chgglobnl,celle_x,celle_y,celle_z,celle_plocnl,chglocnl)
         do i = 1, nionlocnlb
            if (i.le.nionlocnl) then
               iichg         = chgglobnl(celle_key(i))
               iglob         = iion(iichg)
               celle_chg(i)  = iichg
               !celle_plocnl(i) = chglocnl(iichg)
               celle_glob(i) = iglob
               celle_x   (i) = x(iglob)
               celle_y   (i) = y(iglob)
               celle_z   (i) = z(iglob)
               cellv_jvdw(i) = jvdw_c(iglob)
               radv(i)       = rad(jvdw(iglob))
               epsv(i)       = eps(jvdw(iglob))
            else
               celle_chg(i)  = n
               celle_glob(i) = n
               celle_plocnl(i)= nionlocnl
               celle_x   (i) = inf
               celle_y   (i) = inf
               celle_z   (i) = inf
               cellv_jvdw(i) = 1
               radv(i)       = inf
               epsv(i)       = inf
            end if
         end do

         else
!$acc parallel loop async(rec_queue)
!$acc&         present(celle_glob,celle_chg,celle_key,
!$acc&   chgglobnl,celle_x,celle_y,celle_z,celle_plocnl,chglocnl)
         do i = 1, nionlocnlb
            if (i.le.nionlocnl) then
               iichg         = chgglobnl(celle_key(i))
               iglob         = iion(iichg)
               celle_chg(i)  = iichg
               !celle_plocnl(i) = chglocnl(iichg)
               celle_glob(i) = iglob
               celle_x   (i) = x(iglob)
               celle_y   (i) = y(iglob)
               celle_z   (i) = z(iglob)
            else
               celle_chg(i)  = n
               celle_glob(i) = n
               celle_plocnl(i)= nionlocnl
               celle_x   (i) = inf
               celle_y   (i) = inf
               celle_z   (i) = inf
            end if
         end do

         end if ! fuse_chglj

      else

!$acc parallel loop async(rec_queue)
!$acc&         present(celle_glob,celle_chg,celle_key,
!$acc&   celle_loc,celle_ploc,celle_x,celle_y,celle_z,x,y,z,
!$acc&   loc)
         do i = 1, nionlocnlb
            if (i.le.nionlocnl) then
               iglob        = celle_glob(i)
               celle_loc(i) = loc(iglob)
c              celle_ploc(i)= chgloc(celle_chg(i))
               celle_x  (i) = x(iglob)
               celle_y  (i) = y(iglob)
               celle_z  (i) = z(iglob)
            else
               celle_loc(i) = nbloc
c              celle_ploc(i)= nionbloc
               celle_x  (i) = inf
               celle_y  (i) = inf
               celle_z  (i) = inf
            end if
         end do

      end if
      end subroutine

      subroutine build_adjacency_matrix
     &           (nlocnl,bP)
      use domdec  ,only: xbegproc,ybegproc,zbegproc,
     &                   xendproc,yendproc,zendproc,
     &                   rank,nproc
      use cutoff  ,only: use_shortclist,use_shortmlist,
     &                   use_shortvlist
      use cell    ,only: xcell2,ycell2,zcell2
      use inform  ,only: deb_Path
      use neigh   ,only: matb_lst,cell_scan,cell_len1,cell_len2,boxPart
     &            , offl=>offsetlMb,offr=>offsetrMb,szMatb,nbMatb
     &            , ncell_tot,nx_cell,ny_cell,nz_cell
     &            , bit_si,bit_sh,cell_lenr,cell_len2r
      use nblistgpu_inl ,only: midpointimage_inl
      use utils   ,only: set_to_zero1_int1,set_to_zero1_int8
      use utilgpu ,only: prmem_request,rec_queue,BLOCK_SIZE
      implicit none
      integer      ,intent(in):: nlocnl
      type(boxPart),intent(in):: bP

      integer i,iblock,inblock,icell
      integer k,kblock,knblock,kcell
      integer bx_cell,by_cell,bz_cell
      integer it1,it2,block1,block2,nblock,bit2
      integer nx,ny,nz
      integer icell_x,icell_y,icell_z
      integer kcell_x,kcell_y,kcell_z
      integer nxy_cell,icell_scan,kcell_scan
      integer(1) icell_len,kcell_len
      integer setb,capt
      real(t_p) lenx_cell,leny_cell,lenz_cell
      real(t_p) xmin,xmax,ymin,ymax,zmin,zmax
      real(t_p) xbeg,xend,ybeg,yend,zbeg,zend
      real(t_p) xmid_ik,ymid_ik,zmid_ik
      real(t_p) xbi,ybi,zbi,xbk,ybk,zbk,posx,posy,posz

      xmin = -xcell2
      xmax =  xcell2
      ymin = -ycell2
      ymax =  ycell2
      zmin = -zcell2
      zmax =  zcell2
      if (deb_Path) write(*,'(2x,a)') "build_adjacency_matrix"

      ! Get box partition info
      bx_cell = bP%bx_c
      by_cell = bP%by_c
      bz_cell = bP%bz_c

c     zbeg = zbegproc(rank+1)
c     zend = zendproc(rank+1)
c     ybeg = ybegproc(rank+1)
c     yend = yendproc(rank+1)
c     xbeg = xbegproc(rank+1)
c     xend = xendproc(rank+1)

      lenx_cell = abs(xmax-xmin)/nx_cell
      leny_cell = abs(ymax-ymin)/ny_cell
      lenz_cell = abs(zmax-zmin)/nz_cell
      nxy_cell  = nx_cell*ny_cell

      nblock    = ((nlocnl-1)/BLOCK_SIZE) + 1

      call prmem_request (matb_lst,szMatb*nbMatb,async=.true.)
      call prmem_request (cell_lenr,nblock,async=.true.)
      if (use_shortvlist.or.use_shortmlist.or.use_shortclist) then
         call prmem_request (cell_len2 ,nblock,async=.true.)
         call prmem_request (cell_len2r,nblock,async=.true.)
      end if
      call set_to_zero1_int8(matb_lst(1),int(szMatb,8)*nbMatb,
     &     rec_queue)

      ! Find bocks to compute interactions with
!$acc parallel loop present(cell_len1,cell_scan,matb_lst)
!$acc&              async(rec_queue)
      do icell = 1,ncell_tot
         icell_len = cell_len1(icell)
         if (icell_len.eq.0) cycle
         icell_scan = cell_scan(icell)
         icell_x    =  mod(icell-1,nx_cell)
         icell_y    = (icell-1 -((icell-1)/nxy_cell)*nxy_cell)/nx_cell
         icell_z    = (icell-1)/nxy_cell

         ! midpoint box i coordinates
         ! xbi        = xmin + (icell_x+0.5)*lenx_cell
         ! ybi        = ymin + (icell_y+0.5)*leny_cell
         ! zbi        = zmin + (icell_z+0.5)*lenz_cell

         ! retrieve first block number of box i
         iblock     = icell_scan / BLOCK_SIZE
         inblock    = 1
         ! Find number of blocks contained in box i
         if (((icell_scan+icell_len) / BLOCK_SIZE).gt.iblock)
     &      inblock = 2

         ! Skip kblock search if iblock[+1] is outside offsets
         if (inblock.eq.2) then
            if ((iblock  .lt.offl.or.iblock  .ge.offr)
     &     .and.(iblock+1.lt.offl.or.iblock+1.ge.offr)) cycle
         else
            if (iblock.lt.offl.or.iblock.ge.offr) cycle
         end if

         ! loop on neighbors box
!$acc loop vector collapse(3)
         do nz = -bz_cell,bz_cell
            do ny = -by_cell,by_cell
               do nx = -bx_cell,bx_cell

                  kcell_x = icell_x + nx
                  kcell_y = icell_y + ny
                  kcell_z = icell_z + nz

                  if      (kcell_x.lt.0) then
                     kcell_x = kcell_x + nx_cell
                  else if (kcell_x.gt.nx_cell-1) then
                     kcell_x = kcell_x - nx_cell
                  end if
                  if      (kcell_y.lt.0) then
                     kcell_y = kcell_y + ny_cell
                  else if (kcell_y.gt.ny_cell-1) then
                     kcell_y = kcell_y - ny_cell
                  end if
                  if      (kcell_z.lt.0) then
                     kcell_z = kcell_z + nz_cell
                  else if (kcell_z.gt.nz_cell-1) then
                     kcell_z = kcell_z - nz_cell
                  end if

                  kcell = ( kcell_x + nx_cell*kcell_y + 
     &                         nxy_cell*kcell_z ) + 1
                  kcell_len  = cell_len1(kcell)
                  if (kcell_len.eq.0.or.kcell.lt.icell) cycle

                  kcell_scan = cell_scan(kcell)

                  kblock     = kcell_scan/BLOCK_SIZE
                  knblock    = 1
                  ! Find number of blocks contained in box k
                  if (((kcell_scan+kcell_len) / BLOCK_SIZE).gt.kblock)
     &               knblock = 2

                  ! midpoint box k coordinates
c                 xbk    = xmin + (kcell_x+0.5)*lenx_cell
c                 ybk    = ymin + (kcell_y+0.5)*leny_cell
c                 zbk    = zmin + (kcell_z+0.5)*lenz_cell
c                 posx   =  xbi - xbk
c                 posy   =  ybi - ybk
c                 posz   =  zbi - zbk
                  ! Apply mdipointimage rule
                  ! TODO add test when single mpi process
                  ! TODO add test on vuf2 on mid point distance
c                 call midpointimage_inl(xbk,ybk,zbk,posx,posy,posz)
c                 if   ((zbk.ge.zbeg).and.(zbk.lt.zend)
c    &             .and.(ybk.ge.ybeg).and.(ybk.lt.yend)
c    &             .and.(xbk.ge.xbeg).and.(xbk.lt.xend)) then

                  ! Authorize blocks interactions
!$acc loop seq
                  do it1 = 0,inblock-1
                     if (iblock+it1.lt.offl.or.iblock+it1.ge.offr) cycle
                     block1 = (iblock+it1-offl)*szMatb
                     do it2 = 0,knblock-1
                        block2 = block1 +ishft(kblock+it2,-bit_sh) +1
                        bit2   =  iand( kblock+it2,bit_si-1 )
                        setb   = ishft( 1, bit2 )
!$acc atomic
                        matb_lst(block2) = ior (matb_lst(block2),(setb))
                     end do
                  end do

               end do
            end do
         end do
      end do
      end

      subroutine build_adjacency_matrix_octahedron
     &           (nlocnl,bP)
      use boxes   ,only: box34
      use domdec  ,only: xbegproc,ybegproc,zbegproc,
     &                   xendproc,yendproc,zendproc,
     &                   rank,nproc
      use cutoff  ,only: use_shortclist,use_shortmlist,
     &                   use_shortvlist
      use cell    ,only: xcell2,ycell2,zcell2
      use inform  ,only: deb_Path
      use neigh   ,only: matb_lst,cell_scan,cell_len1,cell_len2,boxPart
     &            , offl=>offsetlMb,offr=>offsetrMb,szMatb,nbMatb
     &            , ncell_tot,nx_cell,ny_cell,nz_cell
     &            , bit_si,bit_sh,cell_lenr,cell_len2r
      use nblistgpu_inl
      use utils   ,only: set_to_zero1_int1,set_to_zero1_int8
      use utilgpu ,only: prmem_request,rec_queue,BLOCK_SIZE
      implicit none
      integer      ,intent(in):: nlocnl
      type(boxPart),intent(in):: bP

      integer bx_cell,by_cell,bz_cell
      integer i,iblock,inblock,icell,ilcell,nlcell
      integer k,kblock,knblock,kcell
      integer it1,it2,block1,block2,nblock,bit2
      integer nx,ny,nz
      integer icell_x,icell_y,icell_z
      integer kcell_x,kcell_y,kcell_z
      integer nx_cell2,ny_cell2,nz_cell2
      integer imageCell
      integer nxy_cell,icell_scan,kcell_scan
      integer(1) icell_len,kcell_len
      integer setb,capt
      real(t_p) lenx_cell,leny_cell,lenz_cell
      real(t_p) xmin,xmax,ymin,ymax,zmin,zmax
      real(t_p) xbeg,xend,ybeg,yend,zbeg,zend
      real(t_p) xmid_ik,ymid_ik,zmid_ik
      real(t_p) xbi,ybi,zbi,xbk,ybk,zbk,posx,posy,posz
      real(t_p) cellXmin,cellXmax,cellYmin,cellYmax,cellZmin,cellZmax
     &         ,cellXin,cellYin,cellZin
      logical   close2Zbox,far2Zbox,outcell
!$acc routine(fatal_acc)

      xmin = -xcell2
      xmax =  xcell2
      ymin = -ycell2
      ymax =  ycell2
      zmin = -zcell2
      zmax =  zcell2

      if(deb_Path) write(*,'(2x,a)')"build_adjacency_matrix_octahedron"

      ! Get box partition info
      bx_cell = bP%bx_c
      by_cell = bP%by_c
      bz_cell = bP%bz_c

c     zbeg = zbegproc(rank+1)
c     zend = zendproc(rank+1)
c     ybeg = ybegproc(rank+1)
c     yend = yendproc(rank+1)
c     xbeg = xbegproc(rank+1)
c     xend = xendproc(rank+1)

      lenx_cell = abs(xmax-xmin)/nx_cell
      leny_cell = abs(ymax-ymin)/ny_cell
      lenz_cell = abs(zmax-zmin)/nz_cell
      nxy_cell  = nx_cell*ny_cell
      nx_cell2  = nx_cell/2
      ny_cell2  = ny_cell/2
      nz_cell2  = nz_cell/2

      nblock = ((nlocnl-1)/BLOCK_SIZE) + 1

      call prmem_request (matb_lst,szMatb*nbMatb,async=.true.)
      call prmem_request (cell_lenr,nblock,async=.true.)
      if (use_shortvlist.or.use_shortmlist.or.use_shortclist) then
         call prmem_request (cell_len2 ,nblock,async=.true.)
         call prmem_request (cell_len2r,nblock,async=.true.)
      end if
      call set_to_zero1_int8(matb_lst(1),int(szMatb,8)*nbMatb,
     &     rec_queue)

      ! Find bocks to compute interactions with
!$acc parallel loop present(cell_len1,cell_scan,matb_lst)
!$acc&              async(rec_queue)
      do icell = 1,ncell_tot
         icell_len = cell_len1(icell)
         if (icell_len.eq.0) cycle
         icell_scan = cell_scan(icell)
         icell_x    =  mod(icell-1,nx_cell)
         icell_y    = (icell-1 -((icell-1)/nxy_cell)*nxy_cell)/nx_cell
         icell_z    = (icell-1)/nxy_cell

         ! midpoint box i coordinates
         ! xbi        = xmin + (icell_x+0.5)*lenx_cell
         ! ybi        = ymin + (icell_y+0.5)*leny_cell
         ! zbi        = zmin + (icell_z+0.5)*lenz_cell

         ! retrieve first block number of box i
         iblock     = icell_scan / BLOCK_SIZE
         inblock    = 1
         ! Find number of blocks contained in box i
         if (((icell_scan+icell_len) / BLOCK_SIZE).gt.iblock)
     &      inblock = 2

         ! Skip kblock search if iblock[+1] is outside offsets
         if (inblock.eq.2) then
            if ((iblock  .lt.offl.or.iblock  .ge.offr)
     &     .and.(iblock+1.lt.offl.or.iblock+1.ge.offr)) cycle
         else
            if (iblock.lt.offl.or.iblock.ge.offr) cycle
         end if

         ! loop on neighbors box
!$acc loop vector collapse(3)
         do nz = -bz_cell,bz_cell
            do ny = -by_cell,by_cell
               do nx = -bx_cell,bx_cell

                  kcell_x = icell_x + nx
                  kcell_y = icell_y + ny
                  kcell_z = icell_z + nz

                  ! Find cell image in cubic box
                  if      (kcell_x.lt.0) then
                     kcell_x = kcell_x + nx_cell
                  else if (kcell_x.gt.nx_cell-1) then
                     kcell_x = kcell_x - nx_cell
                  end if
                  if      (kcell_y.lt.0) then
                     kcell_y = kcell_y + ny_cell
                  else if (kcell_y.gt.ny_cell-1) then
                     kcell_y = kcell_y - ny_cell
                  end if
                  if      (kcell_z.lt.0) then
                     kcell_z = kcell_z + nz_cell
                  else if (kcell_z.gt.nz_cell-1) then
                     kcell_z = kcell_z - nz_cell
                  end if

                  ! Where is neighbor cell located compare to octahedron
                  cellXmin = xmin +  kcell_x   *lenx_cell
                  cellXmax = xmin + (kcell_x+1)*lenx_cell
                  if (abs(cellXmin)>abs(cellXmax))
     &               call ExchangeScalar(cellXmin,cellXmax)
                  cellYmin = ymin +  kcell_y   *leny_cell
                  cellYmax = ymin + (kcell_y+1)*leny_cell
                  if (abs(cellYmin)>abs(cellYmax))
     &               call ExchangeScalar(cellYmin,cellYmax)
                  cellZmin = zmin +  kcell_z   *lenz_cell
                  cellZmax = zmin + (kcell_z+1)*lenz_cell
                  if (abs(cellZmin)>abs(cellZmax))
     &               call ExchangeScalar(cellZmin,cellZmax)

                  close2Zbox =(abs(cellXmin)+abs(cellYmin)+abs(cellZmin)
     &                         .le.box34)
                  far2Zbox   =(abs(cellXmax)+abs(cellYmax)+abs(cellZmax)
     &                         .le.box34)


                  ! Find neighbor cell attributes
                  if (close2Zbox.and.far2Zbox) then !Inside box
                     nlcell=1
                     outcell=.false.
                  else if (close2Zbox.and..not.far2Zbox) then !Partially inside box
                     nlcell=2
                     outcell=.false.
                  else if (.not.close2Zbox.and..not.far2Zbox) then !Outside box
                     nlcell=1
                     outcell=.true.
                  else           ! Absurd
                     print*,':build_cell_listgpu'
                     print*,icell,close2Zbox,far2Zbox
                     print*,cellXmin,cellYmin,cellZmin
                     print*,cellXmax,cellYmax,cellZmax
                     print*,'impossible case in neighbor cell'
                     call fatal_acc
                  end if

!$acc loop seq
                  do ilcell = 1,nlcell
                     if (ilcell==1) then
                        if (outcell) then
                           kcell = imageOutCellOcta
     &                               (kcell_x,kcell_y,kcell_z
     &                               ,nx_cell2,ny_cell2,nz_cell2,nx_cell
     &                               ,ny_cell,nz_cell)
                        else
                           kcell = ( kcell_x + nx_cell*kcell_y + 
     &                                  nxy_cell*kcell_z ) + 1
                        end if
                     else
                        kcell = imageOutCellOcta(kcell_x,kcell_y,kcell_z
     &                               ,nx_cell2,ny_cell2,nz_cell2,nx_cell
     &                               ,ny_cell,nz_cell)
                     end if

                     kcell_len  = cell_len1(kcell)
                     if (kcell_len.eq.0.or.kcell.lt.icell) cycle

                     kcell_scan = cell_scan(kcell)

                     kblock     = kcell_scan/BLOCK_SIZE
                     knblock    = 1
                     ! Find number of blocks contained in box k
                     if (((kcell_scan+kcell_len)/BLOCK_SIZE).gt.kblock)
     &                  knblock = 2

                     ! Authorize blocks interactions
!$acc loop seq
                     do it1 = 0,inblock-1
                     if (iblock+it1.lt.offl.or.iblock+it1.ge.offr) cycle
                        block1 = (iblock+it1-offl)*szMatb
                        do it2 = 0,knblock-1
                           block2 = block1 +ishft(kblock+it2,-bit_sh) +1
                           bit2   =  iand( kblock+it2,bit_si-1 )
                           setb   = ishft( 1, bit2 )
!$acc atomic
                        matb_lst(block2) = ior (matb_lst(block2),(setb))
                        end do
                     end do
                  end do

               end do
            end do
         end do
      end do
      end
c
c     Preprocess adjacency matrix to find optimise memory allocation
c     Count every neighors blocks interaction
c
      subroutine pre_process_adjacency_matrix
     &           (nblock,nlocnlb_pair)
      use cutoff  ,only: use_shortclist,use_shortmlist,
     &                   use_shortvlist
      use neigh   ,only: cell_len,cell_scan,matb_lst
     &            ,offl=>offsetlMb,szMatb,nbMatb,cell_scan1
     &            ,bit_si,bit_sh
      use utils   ,only: set_to_zero1_int
      use utilgpu ,only: rec_queue,prmem_request
#ifdef _OPENACC
      use thrust  ,only: thrust_exclusive_scan_int_async
      use utilgpu ,only: rec_stream
#endif
      implicit none
      integer,intent(in):: nblock
      integer,intent(out):: nlocnlb_pair
      integer it1,it2,pos2,bit2

      !Re use cell_len & cell_scan1 for neighbor blocks
      call prmem_request(cell_scan1,nblock,async=.true.)
      call set_to_zero1_int(cell_len(offl+1),nbMatb,rec_queue)

      nlocnlb_pair = 0

!$acc parallel loop async(rec_queue)
!$acc&    present(cell_len,matb_lst)
      do it1 = 0,nbMatb-1
!$acc loop vector
         do it2 = offl+it1,nblock-1
            pos2 = ishft( it2,-bit_sh ) + 1
            bit2 = iand ( it2,bit_si-1)
            if (btest(matb_lst(it1*szMatb+pos2),bit2)) then
!$acc atomic
               cell_len(offl+it1+1) = cell_len(offl+it1+1) + 1
            end if
         end do
      end do
      ! Non blocking reduction with pgi-20

!$acc parallel loop async(rec_queue) present(cell_len)
      do it1 = 1,nbMatb
         nlocnlb_pair=nlocnlb_pair+cell_len(offl+it1)
      end do
!$acc wait(rec_queue)
      ! Apply partial sum to cell_scan1
#ifdef _OPENACC
!$acc host_data use_device(cell_len,cell_scan1)
      call thrust_exclusive_scan_int_async(cell_len(offl+1),nbMatb,
     &                      cell_scan1(offl+1),rec_stream)
!$acc end host_data
#else
      cell_scan1(offl+1) = 0
      do it1 = 2,nbMatb
         cell_scan1(offl+it1) = cell_scan1(offl+it1-1) +
     &                          cell_len(offl+it1-1)
      end do
#endif
c  Wait for nlocnlb_pair to be updated on host memory
c  Implicit Wait seems to be removed with PGI-20

      ! Complete scan when useful
      if (nbMatb.ne.nblock.and.offl.ne.0) then
!$acc parallel loop async(rec_queue) present(cell_scan1,cell_len)
         do it1 = 1, nbMatb
            cell_scan1(offl+it1) = cell_scan1(offl+it1) +
     &                         cell_scan1(offl) + cell_len(offl)
         end do
      end if

      end subroutine
      subroutine pre_process_adjacency_matrix_
     &           (nblock,nlocnlb_pair)
      use cutoff  ,only: use_shortclist,use_shortmlist,
     &                   use_shortvlist
      use neigh   ,only: cell_len,cell_scan,matb_lst
     &            ,offl=>offsetlMb,szMatb,nbMatb,cell_scan1
     &            ,bit_si,bit_sh
      use utils   ,only: set_to_zero1_int
      use utilgpu ,only: rec_queue,prmem_request
#ifdef _OPENACC
      use thrust  ,only: thrust_exclusive_scan_int_async
      use utilgpu ,only: rec_stream
#endif
      implicit none
      integer,intent(in):: nblock
      integer,intent(out):: nlocnlb_pair
      integer it1,it2,pos2,bit2

      !Re use cell_len & cell_scan1 for neighbor blocks
      !call prmem_request(cell_scan1,nblock,async=.true.)
      !call set_to_zero1_int(cell_len(offl+1),nbMatb,rec_queue)

      nlocnlb_pair = 0

!$acc parallel loop async(rec_queue) reduction(+:nlocnlb_pair)
!$acc&         present(matb_lst)
      do it1 = 0,nbMatb-1
!$acc loop vector
         do it2 = offl+it1,nblock-1
            pos2 = ishft( it2,-bit_sh ) + 1
            bit2 = iand ( it2,bit_si-1)
            if (btest(matb_lst(it1*szMatb+pos2),bit2)) then
               nlocnlb_pair=nlocnlb_pair+1
            end if
         end do
      end do
      ! Non blocking reduction with pgi-20
!$acc wait(rec_queue)
      end subroutine
c
c     Build nblist with optimal space according to matb_lst
c     From full storage to adjacency list
c
      subroutine fill_adjacency_list
     &           (nblock,nlocnlb_pair,lst)
      use neigh  ,only: cell_len,cell_lenr,cell_scan1,matb_lst
     &           ,offl=>offsetlMb,offr=>offsetrMb,szMatb,nbMatb
     &           ,bit_si,bit_sh,cell_len2,cell_len2r
      use utilgpu,only: BLOCK_SIZE,rec_queue
      implicit none
      integer,intent(in)   :: nblock,nlocnlb_pair
      integer,intent(inout):: lst(2*nlocnlb_pair)
      integer k,i,iblock,it1,it2,pos2,bit2

!$acc parallel loop async(rec_queue) private(k)
!$acc&         present(cell_scan1,matb_lst,lst)
      do it1 = 0,nbMatb-1
         iblock = cell_scan1(offl+it1+1)
         k      = 0
!$acc loop vector
         do it2 = offl+it1,nblock-1
            pos2 = ishft( it2,-bit_sh ) + 1
            bit2 = iand ( it2,bit_si-1)
            if (btest(matb_lst(it1*szMatb+pos2),bit2)) then
!$acc atomic capture
               i = k
               k = k + 1
!$acc end atomic
               lst(2*(iblock+i)+1) = offl+it1+1
               lst(2*(iblock+i)+2) =      it2+1
            end if
         end do
      end do
      end subroutine
      subroutine fill_adjacency_list_
     &           (nblock,nlocnlb_pair,lst)
      use neigh  ,only: cell_len,cell_lenr,cell_scan1,matb_lst
     &           ,offl=>offsetlMb,offr=>offsetrMb,szMatb,nbMatb
     &           ,bit_si,bit_sh,cell_len2,cell_len2r
      use utilgpu,only: BLOCK_SIZE,rec_queue
      implicit none
      integer,intent(in)   :: nblock,nlocnlb_pair
      integer,intent(inout):: lst(2*nlocnlb_pair)
      integer k,i,iblock,it1,it2,pos2,bit2

      k = 0
!$acc parallel loop async(rec_queue) copyin(k)
!$acc&         present(cell_scan1,matb_lst,lst)
      do it1 = 0,nbMatb-1
!$acc loop vector
         do it2 = offl+it1,nblock-1
            pos2 = ishft( it2,-bit_sh ) + 1
            bit2 = iand ( it2,bit_si-1)
            if (btest(matb_lst(it1*szMatb+pos2),bit2)) then
!$acc atomic capture
               i = k
               k = k + 1
!$acc end atomic
               lst(2*i+1) = offl+it1+1
               lst(2*i+2) =      it2+1
            end if
         end do
      end do
      end subroutine
c
c     Re use cell_len for blocks
c     initiate to block_size
c
      subroutine Init_blocklen(nblock)
      use cutoff  ,only: use_shortclist,use_shortmlist,
     &                   use_shortvlist
      use neigh  ,only: cell_len,cell_lenr,cell_scan,cell_scan1
     &           ,offl=>offsetlMb,offr=>offsetrMb,szMatb,nbMatb
     &           ,bit_si,bit_sh,cell_len2,cell_len2r,matb_lst
      use utilgpu,only: BLOCK_SIZE,rec_queue
      implicit none
      integer,intent(in)::nblock
      integer i

      if (use_shortvlist.or.use_shortmlist.or.use_shortclist) then
!$acc parallel loop async(rec_queue) present(cell_len,cell_lenr,
!$acc&         cell_len2,cell_len2r,cell_scan,cell_scan1)
         do i = 1,nblock
            cell_len  (i) = BLOCK_SIZE
            cell_lenr (i) = 0
            cell_len2 (i) = BLOCK_SIZE
            cell_len2r(i) = 0
            cell_scan (i) = cell_scan1(i)
         end do
      else
!$acc parallel loop async(rec_queue) present(cell_len,cell_lenr,
!$acc&         cell_scan,cell_scan1)
         do i = 1,nblock
            cell_len (i) = BLOCK_SIZE
            cell_lenr(i) = 0
            cell_scan(i) = cell_scan1(i)
         end do
      end if

      end subroutine
c
c     finalize the pair block-atoms list
c     complete lst to fit a multiple of BLOCK_SIZE
c     building his pair block vector
c
      subroutine finalize_list_C2 (nblock,nlocnlb_pair,nlocnlb
     &           ,nlocnlb2_pair,ilst,lst,cell_len,cell_scan,cell_lenr)
      use domdec ,only: rank
      use inform ,only: deb_Path
      use neigh  ,only: cell_scan1
      use utilgpu,only: BLOCK_SIZE,rec_queue,prmem_request
#ifdef _OPENACC
     &           ,rec_stream
      use thrust
#endif
      implicit none
      integer,intent(in):: nblock,nlocnlb,nlocnlb_pair
      integer,intent(out):: nlocnlb2_pair
      integer,intent(inout):: ilst(nlocnlb_pair)
     &       ,lst(nlocnlb_pair*(BLOCK_SIZE+2))
     &       ,cell_len(nblock),cell_scan(nblock)
     &       ,cell_lenr(nblock)
      integer i,k,nneig,nneigr,iv,ivr,nkblock,complete_b
      integer kcell_scanr,kcell_scan
      integer buff_i

      call prmem_request(cell_scan1,nblock,async=.true.)

      nlocnlb2_pair=0
!$acc parallel loop async(rec_queue)
!$acc&         present(ilst,cell_scan,lst,cell_len,cell_lenr)
!$acc&         vector_length(BLOCK_SIZE)
!!$acc&        private(buff)
      do k = 1, nblock
         nneig       = cell_len(k)
         iv          = cell_scan(k)
         nneigr      = cell_lenr(k)
         ! Last block has only itself as a neigbhor
         if (k.le.nblock) ivr = cell_scan(k+1)
         nkblock     = ((nneig+nneigr-1)/BLOCK_SIZE) + 1
         complete_b  = nkblock*BLOCK_SIZE - (nneig+nneigr)
         kcell_scan  = 2*nlocnlb_pair+BLOCK_SIZE*iv
         kcell_scanr = 2*nlocnlb_pair+BLOCK_SIZE*ivr
         cell_len(k) = nkblock
!$acc loop vector
         do i = 1, complete_b
          !buff_i = lst(kcell_scanr-complete_b+i)
          !lst(kcell_scanr-complete_b+i) = nlocnlb
            lst(kcell_scan+nneig+i) = nlocnlb
         end do
c!$acc loop vector
c         do i = 1,nkblock
c            ilst(iv+i) = k
c         end do
           nlocnlb2_pair = nlocnlb2_pair+nkblock
      end do
c Wait for nlocnl_pair
c Implicit wait removed with PGI-20.4
!$acc wait(rec_queue)

#ifdef _OPENACC
!$acc host_data use_device(cell_len,cell_scan1,lst)
      call thrust_exclusive_scan(cell_len,nblock,cell_scan1,rec_stream)
      call thrust_remove_zero_async_int(lst(2*nlocnlb_pair+1)
     &            ,nlocnlb_pair*(BLOCK_SIZE),rec_stream)
!$acc end host_data
#endif

!$acc parallel loop async(rec_queue)
!$acc&         present(cell_len,cell_scan1,ilst)
      do i = 1, nblock
         iv = cell_len(i)
         kcell_scan = cell_scan1(i)
!$acc loop vector
         do k = 1,iv
            ilst(kcell_scan+k)=i
         end do
      end do

      if (deb_Path)
     &   print*,'C1 vs C2 nblist pairs blocks number ',nlocnlb_pair
     &      , nlocnlb2_pair

      end subroutine
      subroutine finalize_list2_C2 (nblock,nlocnlb_pair,nlocnlb
     &           ,nlocnlb2_pair,ilst,lst,cell_len,cell_scan,cell_lenr)
      use domdec ,only: rank
      use inform ,only: deb_Path
      use neigh  ,only: cell_scan1
      use utilgpu,only: BLOCK_SIZE,rec_queue,prmem_request
#ifdef _OPENACC
     &           ,rec_stream
      use thrust
#endif
      implicit none
      integer,intent(in):: nblock,nlocnlb,nlocnlb_pair
      integer,intent(out):: nlocnlb2_pair
      integer,intent(inout):: ilst(nlocnlb_pair)
     &       ,lst(nlocnlb_pair*BLOCK_SIZE)
     &       ,cell_len(nblock),cell_scan(nblock)
     &       ,cell_lenr(nblock)
      integer i,k,nneig,nneigr,iv,ivr,nkblock,complete_b
      integer kcell_scanr,kcell_scan
      integer buff_i

      call prmem_request(cell_scan1,nblock,async=.true.)

      nlocnlb2_pair=0
!$acc parallel loop async(rec_queue)
!$acc&         present(ilst,cell_scan,lst,cell_len,cell_lenr)
!$acc&         vector_length(BLOCK_SIZE)
!!$acc&        private(buff)
      do k = 1, nblock
         nneig       = cell_len(k)
         iv          = cell_scan(k)
         nneigr      = cell_lenr(k)
         ! Last block has only itself as a neigbhor
         if (k.le.nblock) ivr = cell_scan(k+1)
         nkblock     = ((nneig+nneigr-1)/BLOCK_SIZE) + 1
         complete_b  = nkblock*BLOCK_SIZE - (nneig+nneigr)
         kcell_scan  = BLOCK_SIZE*iv
         kcell_scanr = BLOCK_SIZE*ivr
         cell_len(k) = nkblock
!$acc loop vector
         do i = 1, complete_b
          !buff_i = lst(kcell_scanr-complete_b+i)
          !lst(kcell_scanr-complete_b+i) = nlocnlb
            lst(kcell_scan+nneig+i) = nlocnlb
         end do
c!$acc loop vector
c         do i = 1,nkblock
c            ilst(iv+i) = k
c         end do
           nlocnlb2_pair = nlocnlb2_pair+nkblock
      end do
c Wait for nlocnl_pair
c Implicit wait removed with PGI-20.4
!$acc wait(rec_queue)

#ifdef _OPENACC
!$acc host_data use_device(cell_len,cell_scan1,lst)
      call thrust_exclusive_scan(cell_len,nblock,cell_scan1,rec_stream)
      call thrust_remove_zero_async_int(lst,nlocnlb_pair*(BLOCK_SIZE)
     &                                 ,rec_stream)
!$acc end host_data
#endif

!$acc parallel loop async(rec_queue)
!$acc&         present(cell_len,cell_scan1,ilst)
      do i = 1, nblock
         iv = cell_len(i)
         kcell_scan = cell_scan1(i)
!$acc loop vector
         do k = 1,iv
            ilst(kcell_scan+k)=i
         end do
      end do

      if (deb_Path)
     &   print*,'C1 vs C2 nblist pairs blocks number ',nlocnlb_pair
     &      , nlocnlb2_pair

      end subroutine

      subroutine SpatialReorder
      use atoms ,only: n
      use atmlst,only: poleglobnl,poleglob,polerecglob
     &          ,chgglob,chgrecglob
      use charge,only: nionrecloc,nionloc,iion,nion
      use domdec,only: nproc,rank,nrec,ndir
      use mpole ,only: npolelocnl,ipole,npoleloc,npolerecloc,npole
      use inform,only: deb_Path
#ifdef _OPENACC
      use interfaces,only:grid_uind_site_p,grid_uind_sitecu
     &              ,grid_pchg_site_p,grid_pchg_sitecu
#endif
      use neigh
      use potent,only: use_pmecore,use_mpole,use_charge
      use utilgpu
      implicit none
      integer i,iipole,iichg,iglob
      integer num
      integer,pointer:: poleArray(:),Array1(:)

#ifndef _OPENACC
      return
#endif

      if (use_pmecore) then
         if (rank.lt.ndir) return
         if (nrec.ne.1) return

         if (use_mpole) then
#ifdef _OPENACC
         if (.not.associated(grid_uind_site_p,grid_uind_sitecu)) return
#endif
            num = npolerecloc
            poleArray => polerecglob
         else if (use_charge) then
#ifdef _OPENACC
         if (.not.associated(grid_pchg_site_p,grid_pchg_sitecu)) return
#endif
            num = nionrecloc
            poleArray => chgrecglob
         end if
      else
         if (nproc.ne.1) return
         if (use_mpole) then
#ifdef _OPENACC
            if (.not.(.not.mlst2_enable
     &          .and.associated(grid_uind_site_p,grid_uind_sitecu)))
     &         return
#endif
            num = npoleloc
            poleArray => poleglob
         else if (use_charge) then
#ifdef _OPENACC
            if (.not.(.not.clst2_enable
     &          .and.associated(grid_pchg_site_p,grid_pchg_sitecu)))
     &         return
#endif
            num = nionloc
            poleArray => chgglob
         end if

      end if
      Array1 => poleArray

      if (deb_Path) write(*,*) 'SpatialReorder'

      call prmem_request(celle_key,num,async=.true.)
      call prmem_request(celle_glob,num,async=.true.)

      if (use_mpole) then
         if (n.ne.num) then
         call prmem_request(celle_glob,num,async=.true.)
!$acc    parallel loop async default(present)
         do i = 1,num
            celle_glob(i) = ipole(poleArray(i))
         end do
         poleArray => celle_glob
         end if
         call build_cell_list2(poleArray,celle_key,num,
     &              sqrt(mbuf2),e_bP)
         call prmem_request(celle_pole,num,async=.true.)
!$acc    parallel loop default(present) async
         do i = 1, npolerecloc
            iipole = Array1(celle_key(i))
            iglob  = ipole(iipole)
            celle_pole(i) = iipole
            celle_glob(i) = iglob
         end do
      else if (use_charge) then
         if (n.ne.num) then
         call prmem_request(celle_glob,num,async=.true.)
!$acc    parallel loop async default(present)
         do i = 1,num
            celle_glob(i) = iion(poleArray(i))
         end do
         poleArray => celle_glob
         end if
         call build_cell_list2(poleArray,celle_key,num,
     &              sqrt(cbuf2),c_bP)
         call prmem_request(celle_chg,nionrecloc,async=.true.)
!$acc    parallel loop default(present) async
         do i = 1, nionrecloc
            iichg  = Array1(celle_key(i))
            iglob  = iion(iichg)
            celle_chg (i) = iichg
            celle_glob(i) = iglob
         end do
      end if
      end subroutine

c
c     "build_pairwise_list" constructs a pairwise atom nblist
c          inside 'a_plist' structure
c
      subroutine build_pairwise_list(atList,natmnl,cut2,a_plist)
      use atoms   ,only: x,y,z,n
      use atmlst  ,only: poleglobnl
      use boxes   ,only: octahedron
      use cutoff  ,only: use_shortmlist
      use cell    ,only: xcell2,ycell2,zcell2,xcell,ycell,zcell
      use domdec  ,only: rank,nproc,xbegproc,xendproc
     &            ,ybegproc,yendproc,zbegproc,zendproc
      use neigh   ,only: pair_atlst,cell_scan,cell_len,matb_lst
     &            ,nbMatb,niterMatb,buffMatb,offsetlMb,offsetrMb,szMatb
     &            ,cell_lenr
     &            ,ieblst,bit_si,bit_sh
     &            ,cell_len2,cell_len2r,build_cell_list2
      use utils   ,only: set_to_zero1_int1,set_to_zero1_int
      use inform  ,only: deb_Path
      use interfaces,only:pre_process_adjacency_matrix
     &               ,set_ElecData_CellOrder
#ifdef _OPENACC
     &               ,cu_filter_lst_sparse
      use thrust  ,only: thrust_exclusive_scan,thrust_remove
      use nblistcu,only: filter_pairwise_atom_lst
      use utilcu  ,only: check_launch_kernel
      use utilgpu ,only: rec_stream
      use cudafor
#endif
      use utilgpu ,only: rec_queue,BLOCK_SIZE,inf
      use tinheader,only: ti_p
      use tinMemory,only: prmem_request,prmem_mvrequest,mipk,szoi
      implicit none
      integer  ,intent(in)::natmnl
      integer  ,intent(in)::atList(:)
      real(t_p),intent(in):: cut2
      type(pair_atlst),target:: a_plist

      integer natmnlb,nbpairs,npairs
      integer nb_pair,szlist,nblock
      integer iglob,iv,ierrSync
      integer i,j,iter
      integer(1) icell_len,kcell_len
      integer(mipk) szolst
      real(t_p) lenx_cell,leny_cell,lenz_cell,rdn,rdn1
      integer  ,pointer:: list(:),blist(:),c_glob(:)
      real(t_p),pointer:: cell_x(:),cell_y(:),cell_z(:)
      real(t_p) xbeg,xend,ybeg,yend,zbeg,zend

      if (deb_Path) print*, 'build_pairwise_list'

      nbpairs = 0
      natmnlb = BLOCK_SIZE*((natmnl-1)/BLOCK_SIZE + 1)
      ! Add margin to serve as out-of-bound exclusion interactions for C2 nblist
      if (natmnlb == natmnl) natmnlb = natmnl + BLOCK_SIZE

      nblock     = ((natmnl-1)/BLOCK_SIZE) + 1
      szMatb     = ((nblock-1)/bit_si) + 1

      nbMatb     = max(1,min(int(buffMatb*1d9/(4*szMatb)),nblock))
      niterMatb  = (nblock-1)/nbMatb + 1
      nbpairs = 0

      call prmem_request(a_plist%c_key,n,async=.true.)

      call build_cell_list2(atList,a_plist%c_key,natmnl,
     &                      sqrt(cut2),a_plist%bPar)
c
c     Build Adjacency list piece by piece
c
      do iter = 0,niterMatb-1

        offsetlMb = iter*nbMatb
        offsetrMb = min((iter+1)*nbMatb,nblock)
        if (iter.eq.niterMatb-1) nbMatb=nblock-offsetlMb

        if (octahedron) then
           call build_adjacency_matrix_octahedron(natmnl,a_plist%bPar)
        else
           call build_adjacency_matrix(natmnl,a_plist%bPar)
        end if

        call pre_process_adjacency_matrix (nblock,nb_pair)
        nbpairs = nbpairs + nb_pair

        call prmem_mvrequest (a_plist%blist,nbpairs*2
     &                       ,async=.true.)
        call fill_adjacency_list(nblock,nbpairs,a_plist%blist)

      end do

      szlist = 2*nbpairs*(BLOCK_SIZE**2)
      call prmem_request(a_plist%list,szlist,async=.true.)
      list  => a_plist%list
      blist => a_plist%blist

      ! zero list and ilist
!$acc parallel loop async present(list)
      do i = 1,nbpairs*2*(BLOCK_SIZE)**2
         list(i) = -1
      end do

      ! Set pairwise list attributes
      a_plist%nbpairs   = nbpairs
      a_plist%natmnl    = natmnl
      a_plist%natmnlb   = natmnlb
      a_plist%cut_buff2 = cut2

      call Init_blocklen(nblock)
      call set_pairlist_Cellorder(atList,a_plist,.true.)
c
c     Debug print
c
      if (deb_Path) then
 13   format( I3," iteration(s) required to build Adj lst" )
c14   format( I10," nblocks process at once over",I10 )
c15   format(1x,'Size Adj Mat (Mo)',F12.3,';',4x,'buffer (Go)',F10.6)
 16   format(1x,'natmnl',I8,';',3x,'mlist mem space (Mo)',F12.3)
         !nbMatb = max(1,min(int(buffMatb*1d9/(4*szMatb)),nblock))
         szolst = nbpairs*(BLOCK_SIZE**2+2)*szoi
         !if (use_shortmlist) szolst = 2*szolst
         print 13, niterMatb
         !print 14, nbMatb,nblock
         !print 15, int(szMatb*4,8)*nbMatb/1d6, buffMatb
         print 16, natmnl, szolst/1d6
         print*, 'sz list',szlist,'; cutbuf2',cut2
      end if

c
c     From Block-Block to Block-Atoms list
c
#ifdef _CUDA
      c_glob=> a_plist%c_glob
      cell_x=> a_plist%cell_x
      cell_y=> a_plist%cell_y
      cell_z=> a_plist%cell_z
      xbeg  =  xbegproc(rank+1)
      xend  =  xendproc(rank+1)
      ybeg  =  ybegproc(rank+1)
      yend  =  yendproc(rank+1)
      zbeg  =  zbegproc(rank+1)
      zend  =  zendproc(rank+1)

!$acc host_data use_device(c_glob,cell_scan,
!$acc&    cell_x,cell_y,cell_z,
!$acc&    matb_lst,list,blist)

      call filter_pairwise_atom_lst<<<*,128,0,rec_stream>>>
     &     (c_glob,cell_scan,cell_x,cell_y,cell_z,matb_lst,
     &      n,natmnlb,nblock,szMatb,nbpairs,cut2,
     &      list,blist,xbeg,xend,ybeg,yend,zbeg,zend)
      call check_launch_kernel(" filter_pairwise_atom_lst")

!$acc end host_data
#else
      !TODO Add an OpenACC version of this kernel
#endif

#ifdef _OPENACC
!$acc wait(rec_queue)
!$acc host_data use_device(list)
      call thrust_remove(list,szlist,-1,nb_pair,rec_stream)
!$acc end host_data
#endif

      a_plist%npairs = nb_pair/2

      if (deb_Path) then
      print*, 'C1 v pair-atom numbers',a_plist%nbpairs,a_plist%npairs
      end if

c     call finalize_list_C2
c    &     (nblock,nbpairs,natmnlb
c    &     ,npolelocnlb2_pair,ieblst,eblst,cell_len,cell_scan
c    &     ,cell_lenr)

c     if (use_shortmlist) call finalize_list_C2
c    &        (nblock,nbpairs,natmnlb
c    &        ,nshortpolelocnlb2_pair,ishorteblst,shorteblst
c    &        ,cell_len2,cell_scan,cell_len2r)


      end subroutine
