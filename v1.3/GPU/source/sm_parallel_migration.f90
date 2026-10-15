!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
#include "tinker_macro.h"
submodule (parallel) parallel_migration
   implicit none

   ! ReassignRespa Exchange space
   type, private :: t_elt
      real(r_p) :: x,y,z,vx,vy,vz,ax,ay,az
   end type
   integer, save :: max_atoms_send=0
   integer, save :: max_atoms_recv=0
   type(t_elt), allocatable, target, save :: data_send(:,:), data_recv(:,:) ! Data (x,y,...) to exchange

contains
#include "image.inc.f90"
#include "convert.inc.f90"

module subroutine release_migration_workspace()
   if (allocated(data_send)) then
!$acc exit data delete(data_send)
      deallocate(data_send)
   end if
   if (allocated(data_recv)) then
!$acc exit data delete(data_recv)
      deallocate(data_recv)
   end if
   max_atoms_send = 0
   max_atoms_recv = 0
end subroutine release_migration_workspace

module subroutine migrate_atoms_legacy()

   if (associated(active_domain)) call check_domain_layout()
   call reassign
   if (associated(active_domain)) call refresh_domain_layout(active_domain)
end subroutine migrate_atoms_legacy

module subroutine reassign_respa(ialt, nalt)
   integer, intent(in) :: ialt, nalt

   call reassignrespa(ialt, nalt)
   if (associated(active_domain)) then
      call refresh_domain_layout(active_domain)
   end if
end subroutine reassign_respa
!
!     subroutine reassign: deal with particles changing of domain between two time steps
!
#ifdef _OPENACC
subroutine reassign
   implicit none
   call reassignrespa(0,0)
end subroutine
#else
subroutine reassign
   use atomsMirror
   use bound
   use cell
   use domdec
   use parallel, only: get_domain_comm, get_direct_comm
   use freeze
   use moldyn
   use neigh
   use potent
   use timestat ,only: timer_enter,timer_exit,timer_eneig&
      &,quiet_timers
   use tinheader,only: ti_eps
   use mpi
   use qtb, only: qtb_thermostat, reassignqtb
   implicit none
   integer i
   integer iproc,ierr,iglob

   real(t_p) xr,yr,zr,eps1,eps2
   integer tag,status(MPI_STATUS_SIZE)
   integer, allocatable :: reqsend(:),reqrec(:)
   integer, allocatable :: count2(:)
   integer, allocatable :: buf3(:), buflen3(:), bufbeg3(:)
   integer, allocatable :: buflen4(:),bufbeg4(:)
   real(r_p), allocatable :: buffers(:,:),buffer(:,:)
   real(r_p), allocatable :: buffers1(:,:),buffer1(:,:)
   integer, allocatable :: buf3bis(:,:)
   integer, allocatable :: glob1(:)
   integer rankloc,nprocloc,proc,j,jglob
   integer nloc1
   integer nloc_save,max_atoms_recv,max_atoms_send
   integer, allocatable :: iglob_send(:,:), iglob_recv(:,:)
   integer, allocatable :: n_data_send(:), n_data_recv(:)
#if (defined(SINGLE)||defined(MIXED))
   parameter(eps1 = 1d1*ti_eps, eps2 = 1d3*ti_eps)
#else
   parameter(eps1 = 1d2*ti_eps, eps2 = 1d4*ti_eps)
#endif
!
   if (nproc.eq.1) return

   call timer_enter( timer_eneig )
   if (use_pmecore) then
      nprocloc = ndir
      rankloc  = rank_bis
   else
      nprocloc = nproc
      rankloc  = rank
   end if
   allocate (glob1(n))
!
   allocate (reqsend(nproc))
   allocate (reqrec(nproc))
!
   allocate (count2(nproc))
   allocate (buflen3(nproc))
   allocate (bufbeg3(nproc))
   allocate (buf3bis(nloc,nproc))

   count2  = 0
   buflen3 = 0
   bufbeg3 = 0
!
!     get particules that changed of domain
!
   do i = 1, nloc
      iglob = glob(i)
      xr    = x(iglob)
      yr    = y(iglob)
      zr    = z(iglob)
      call image(xr,yr,zr)
      if ((xcell2-abs(xr)).lt.eps_cell) xr = xr-sign(4*eps_cell,xr)
      if ((ycell2-abs(yr)).lt.eps_cell) yr = yr-sign(4*eps_cell,yr)
      if ((zcell2-abs(zr)).lt.eps_cell) zr = zr-sign(4*eps_cell,zr)
      do iproc = 0, nprocloc-1
         if (iproc.eq.rank) cycle
         if ((zr.ge.zbegproc(iproc+1)).and.(zr.lt.zendproc(iproc+1))&
            &.and.(yr.ge.ybegproc(iproc+1)).and.(yr.lt.yendproc(iproc+1))&
            &.and.(xr.ge.xbegproc(iproc+1)).and.(xr.lt.xendproc(iproc+1)))&
            &then
            buflen3(iproc+1) = buflen3(iproc+1)+1
            buf3bis(buflen3(iproc+1),iproc+1) = iglob
         end if
      end do
   end do
   bufbeg3(1) = 1
   do iproc = 1, nproc-1
      bufbeg3(iproc+1) = bufbeg3(iproc)+buflen3(iproc)
   end do
!
   allocate (buf3(bufbeg3(nproc)+buflen3(nproc)))
   buf3 = 0
!
   do iproc = 1, nproc
      buf3(bufbeg3(iproc):(bufbeg3(iproc)+buflen3(iproc)-1)) =&
         &buf3bis(1:buflen3(iproc),iproc)
   end do
!
   if (nprocloc.eq.1) then
      call timer_exit( timer_eneig,quiet_timers )
      return
   end if
!
   if ((use_pmecore).and.(rank.ge.ndir)) goto 20
!
!     get size of buffers
!
   allocate (buflen4(nprocloc))
   buflen4 = 0
   allocate (bufbeg4(nprocloc))
   bufbeg4 = 0
   allocate(n_data_send(0:nneig_send))
   n_data_send(:)=0
   allocate(n_data_recv(nneig_recep))
   n_data_recv(:)=0
!
   do iproc = 1, nneig_recep
      proc = pneig_recep(iproc)
      tag  = nprocloc*rankloc + proc + 1
      call MPI_IRECV(buflen4(proc+1),1,MPI_INT,&
         &proc,tag,merge(get_direct_comm(),get_domain_comm(),use_pmecore),reqrec(iproc),ierr)
   end do
   do iproc = 1, nneig_send
      proc = pneig_send(iproc)
      tag  = nprocloc*proc + rankloc + 1
      call MPI_ISEND(buflen3(proc+1),1,MPI_INT,&
         &proc,tag,merge(get_direct_comm(),get_domain_comm(),use_pmecore),reqsend(iproc),ierr)
   end do
!
   do iproc = 1, nneig_recep
      proc = pneig_recep(iproc)
      call MPI_WAIT(reqrec(iproc),status,ierr)
      n_data_recv(iproc)=buflen4(proc+1)
   end do
   do iproc = 1, nneig_send
      proc = pneig_recep(iproc)
      call MPI_WAIT(reqsend(iproc),status,ierr)
      n_data_send(iproc)=buflen3(proc+1)
   end do
!
   bufbeg4(pneig_recep(1)+1) = 1
   do iproc = 2, nneig_recep
      bufbeg4(pneig_recep(iproc)+1) = bufbeg4(pneig_recep(iproc-1)+1)&
         &+buflen4(pneig_recep(iproc-1)+1)
   end do
!
!     communicate the corresponding indexes, positions, speed and accelerations
!
   proc = pneig_send(nneig_send)
   allocate (buffers(10,bufbeg3(proc+1)+buflen3(proc+1)))
   proc = pneig_recep(nneig_recep)
   allocate (buffer(10,bufbeg4(proc+1)+buflen4(proc+1)))
!
!     Begin reception
!
   do i = 1, nneig_recep
      proc = pneig_recep(i)
      tag = nprocloc*rankloc + proc + 1
      call MPI_IRECV(buffer(1,bufbeg4(proc+1)),10*buflen4(proc+1),&
         &MPI_RPREC,proc,tag,get_domain_comm(),reqrec(i),ierr)
   end do

   max_atoms_send=maxval(n_data_send)
   allocate(iglob_send(nloc,0:nneig_send))
   iglob_send(:,:)=0
   do i = 1, nneig_send
      proc = pneig_send(i)
      do j = 0, buflen3(proc+1)-1
         jglob = buf3(bufbeg3(proc+1)+j)
         iglob_send(j+1,i)=jglob
         buffers(1,bufbeg3(proc+1)+j) = x(jglob)
         buffers(2,bufbeg3(proc+1)+j) = y(jglob)
         buffers(3,bufbeg3(proc+1)+j) = z(jglob)
         buffers(4,bufbeg3(proc+1)+j) = v(1,jglob)
         buffers(5,bufbeg3(proc+1)+j) = v(2,jglob)
         buffers(6,bufbeg3(proc+1)+j) = v(3,jglob)
         buffers(7,bufbeg3(proc+1)+j) = a(1,jglob)
         buffers(8,bufbeg3(proc+1)+j) = a(2,jglob)
         buffers(9,bufbeg3(proc+1)+j) = a(3,jglob)
         buffers(10,bufbeg3(proc+1)+j) = jglob
      end do
   end do
!
!     send the positions, the velocities and the accelerations
!
   do i = 1, nneig_send
      proc = pneig_send(i)
      tag = nprocloc*proc + rankloc + 1
      call MPI_ISEND(buffers(1,bufbeg3(proc+1)),10*buflen3(proc+1),&
         &MPI_RPREC,proc,tag,get_domain_comm(),reqsend(i),ierr)
   end do
   do i = 1, nneig_send
      proc = pneig_send(i)
      call MPI_WAIT(reqsend(i),status,ierr)
   end do

   max_atoms_recv=maxval(n_data_recv)
   allocate(iglob_recv(max_atoms_recv,nneig_recep))
   iglob_recv(:,:)=0
   do i = 1, nneig_recep
      proc = pneig_recep(i)
      call MPI_WAIT(reqrec(i),status,ierr)
      do j = 0, buflen4(proc+1)-1
         iglob = int(buffer(10,bufbeg4(proc+1)+j))
         iglob_recv(j+1,i)=iglob
         x(iglob) = buffer(1,bufbeg4(proc+1)+j)
         y(iglob) = buffer(2,bufbeg4(proc+1)+j)
         z(iglob) = buffer(3,bufbeg4(proc+1)+j)
         v(1,iglob) = buffer(4,bufbeg4(proc+1)+j)
         v(2,iglob) = buffer(5,bufbeg4(proc+1)+j)
         v(3,iglob) = buffer(6,bufbeg4(proc+1)+j)
         a(1,iglob) = buffer(7,bufbeg4(proc+1)+j)
         a(2,iglob) = buffer(8,bufbeg4(proc+1)+j)
         a(3,iglob) = buffer(9,bufbeg4(proc+1)+j)
      end do
   end do
!$acc data present(x,y,z,v,a)
!$acc update device(x,y,z,v,a) async
!$acc end data
!
!     if rattle is used, also send the old coordinates
!
   if (use_rattle) then

      proc = pneig_send(nneig_send)
      allocate (buffers1(4,bufbeg3(proc+1)+buflen3(proc+1)))
      proc = pneig_recep(nneig_recep)
      allocate (buffer1(4,bufbeg4(proc+1)+buflen4(proc+1)))
!
!     Begin reception
!
      do i = 1, nneig_recep
         proc = pneig_recep(i)
         tag = nprocloc*rankloc + proc + 1
         call MPI_IRECV(buffer1(1,bufbeg4(proc+1)),4*buflen4(proc+1),&
            &MPI_RPREC,proc,tag,get_domain_comm(),reqrec(i),ierr)
      end do
      do i = 1, nneig_send
         proc = pneig_send(i)
         do j = 0, buflen3(proc+1)-1
            jglob = buf3(bufbeg3(proc+1)+j)
            buffers1(1,bufbeg3(proc+1)+j) = jglob
            buffers1(2,bufbeg3(proc+1)+j) = xold(jglob)
            buffers1(3,bufbeg3(proc+1)+j) = yold(jglob)
            buffers1(4,bufbeg3(proc+1)+j) = zold(jglob)
         end do
      end do
!
!       send the old positions
!
      do i = 1, nneig_send
         proc = pneig_send(i)
         tag = nprocloc*proc + rankloc + 1
         call MPI_ISEND(buffers1(1,bufbeg3(proc+1)),4*buflen3(proc+1),&
            &MPI_RPREC,proc,tag,get_domain_comm(),reqsend(i),ierr)
      end do
      do i = 1, nneig_send
         proc = pneig_send(i)
         call MPI_WAIT(reqsend(i),status,ierr)
      end do
      do i = 1, nneig_recep
         proc = pneig_recep(i)
         call MPI_WAIT(reqrec(i),status,ierr)
         do j = 0, buflen4(proc+1)-1
            iglob = int(buffer1(1,bufbeg4(proc+1)+j))
            xold(iglob) = buffer1(2,bufbeg4(proc+1)+j)
            yold(iglob) = buffer1(3,bufbeg4(proc+1)+j)
            zold(iglob) = buffer1(4,bufbeg4(proc+1)+j)
         end do
      end do
      deallocate (buffer1)
      deallocate (buffers1)
   end if

!
!     reorder indexes accordingly and build local repart array
!
   repart = -1
   repartrec = -1
!
!     remove atoms that left the domain
!
   nloc1 = 0
   glob1 = 0
   do i = 1, nloc
      iglob = glob(i)
      do j = 1, bufbeg3(nproc)+buflen3(nproc)-1
         jglob = buf3(j)
         if (iglob.eq.jglob) goto 10
      end do
      nloc1 = nloc1 + 1
      glob1(nloc1)  = iglob
      loc(iglob)    = nloc1
      repart(iglob) = rank
10    continue
   end do
!
!     add atoms that entered the domain
!
   proc = pneig_recep(nneig_recep)+1
   do j = 1, bufbeg4(proc)+buflen4(proc)-1
      jglob         = int(buffer(10,j))
      nloc1         = nloc1 + 1
      glob1(nloc1)  = jglob
      loc(jglob)    = nloc1
      repart(jglob) = rank
   end do
   deallocate (buffer)
   deallocate (buffers)
!
   if(qtb_thermostat) call reassignqtb(nloc,max_atoms_recv&
      &,nneig_recep, n_data_recv, iglob_recv, pneig_recep&
      &,nneig_send , n_data_send, iglob_send, pneig_send&
      &,max_atoms_recv, max_atoms_send)

   nloc = nloc1
   domlen(rank+1) = nloc
   glob = glob1
!
20 call orderbuffer(.false.)
!
   call timer_exit( timer_eneig,quiet_timers )
   deallocate (buf3bis)
   deallocate (glob1)
   deallocate (count2)
   deallocate (buf3)
   deallocate (buflen3)
   deallocate (bufbeg3)
   deallocate (reqsend)
   deallocate (reqrec)
end
#endif

subroutine reassignrespa(ialt,nalt)
   use atoms       ,only: pbcWrap
   use atomsMirror ,only: x,y,z,xold,yold,zold,n
   use bound       ,only: use_bounds
   use boxes
   use cell
   use domdec,only : nloc,nproc,glob,rank,rank_bis,ndir,nxdd,nydd,nzdd&
             ,xbegproc,ybegproc,zbegproc,xendproc,yendproc,zendproc&
             ,repart,repartrec,loc,locrec,domlen&
             ,nneig_send,nneig_recep,pneig_recep,pneig_send
   use parallel, only: get_domain_comm, get_direct_comm
   use freeze   ,only : use_rattle
   use inform   ,only : deb_Path,tindPath
   use iso_c_binding
   use moldyn   ,only : a,v
   use potent   ,only: use_pmecore
   use timestat ,only: timer_enter,timer_exit,timer_eneig&
      &,quiet_timers
   use tinTypes ,only: real3,real6
   use tinheader,only: ti_eps
   use tinMemory,only: prmem_request,debMem,mem_inc,extra_alloc
   use mpi
   use sizes    ,only: tinkerdebug
   use qtb      ,only: qtb_thermostat, reassignqtb
   use mdstuf   ,only: mts
   use moldyn   ,only: aalt,aalt2
   use colvars  ,only: use_colvars
   implicit none
   integer, intent(in) :: ialt, nalt

   integer,pointer :: iglob_send(:,:) ! Global index to send to each proc
   integer,pointer :: iglob_recv(:,:) ! Global index to recieve from procs
   real(t_p) xr, yr, zr ! adjusted position
   type(t_elt),pointer :: d
   type(real3),pointer :: b
   type(real3),pointer :: old_send(:,:),old_recv(:,:)
   type(real6),pointer :: aalt_send(:,:),aalt_recv(:,:),paa
   type(c_ptr) void_p
   integer, pointer :: pool_i1(:), pool_i2(:)

   integer  n_data_send(0:nneig_send),   n_data_recv(nneig_recep)
   integer req_iglob_send(nneig_send),req_iglob_recv(nneig_recep)
   integer  req_data_send(nneig_send), req_data_recv(nneig_recep)
   integer  mpi_status1(MPI_STATUS_SIZE)

   integer nprocloc, rankloc ! local (pme) mpi info
   integer i, ierr, iglob, iloc, iproc, ineighbor
   integer nloc_save, n_data_send_capture,max_data_recv
   integer n_ar, n_br, n_cr
   integer nloc_capture
   integer s_bufi,pbcw
   integer :: max_data_recv_save=0
   logical ifind
   character(len=40) ::fmt
   real(t_p) lena_box, lenb_box, lenc_box
   real(t_p) ixbeg,ixend,iybeg,iyend,izbeg,izend,xtemp
#if  (defined(SINGLE)||defined(MIXED))
   real(t_p), parameter:: eps1= 10*ti_eps, eps2= 1d2*ti_eps
#else
   real(t_p), parameter:: eps1= 1d2*ti_eps, eps2= 1d4*ti_eps
#endif

   if (ialt/=nalt) return
   if (use_pmecore) then
      nprocloc = ndir
      rankloc  = rank_bis
      if (rank >= ndir) then
         call timer_enter( timer_eneig )
         goto 20
      end if
   else
      nprocloc = nproc
      rankloc  = rank
   end if
   if (nproc.eq.1) return
   if (use_pmecore.and.ndir.eq.1) then
      call timer_enter( timer_eneig )
      goto 20
   end if

   if (deb_Path) write(*,'(A,2I4)') '  -- reassignrespa ',ialt,nalt
   call timer_enter( timer_eneig )

   s_bufi = 2*nloc*(nneig_send+1)
   call domain_workspace_i(1,s_bufi,pool_i1)
   !Associate iglob_send pointer to buffer
   iglob_send(1:2*nloc,0:nneig_send) => pool_i1(1:s_bufi)

!$acc wait
!$acc data copyin(pneig_recep, pneig_send) async &
!$acc     create(n_data_send) &
!$acc   present(xbegproc,xendproc,ybegproc,yendproc,zbegproc,zendproc) &
!$acc     present(repart,loc,glob,x,y,z,a,v,iglob_send)

   ! Gather process domain of rank
   ixbeg = xbegproc(rank+1)
   ixend = xendproc(rank+1)
   iybeg = ybegproc(rank+1)
   iyend = yendproc(rank+1)
   izbeg = zbegproc(rank+1)
   izend = zendproc(rank+1)
   lena_box= xcell/nxdd
   lenb_box= ycell/nydd
   lenc_box= zcell/nzdd

!$acc parallel loop async
   do i = 0,nneig_send
      n_data_send(i) = 0
   end do

   ! Get process domain of each atoms
!$acc parallel loop async
   do i = 1, nloc
      iglob = glob(i)
      xr    = x(iglob)
      yr    = y(iglob)
      zr    = z(iglob)
      pbcw  = pbcWrap(iglob)
      if (monoclinic.or.triclinic) call ctf_inl( xr,yr,zr)
      call image_orth_inl(xr,yr,zr)
      ! Check box limits (SP Issue)
      xtemp = xr
      ifind = .false.
      if ((xbox2-abs(xr)).lt.eps_cell) xr = xr-sign(4*eps_cell,xr)
      if ((ybox2-abs(yr)).lt.eps_cell) yr = yr-sign(4*eps_cell,yr)
      if ((zbox2-abs(zr)).lt.eps_cell) zr = zr-sign(4*eps_cell,zr)
      n_ar  = int((xr+xcell2)/lena_box)
      n_br  = int((yr+ycell2)/lenb_box)
      n_cr  = int((zr+zcell2)/lenc_box)
      iproc = n_ar + n_br*nxdd + n_cr*nydd*nxdd
      if (iproc.eq.rankloc) then
!$acc atomic capture
         n_data_send(0)      = n_data_send(0) + 1
         n_data_send_capture = n_data_send(0)
!$acc end atomic
         iglob_send(2*(n_data_send_capture-1)+1,0)=iglob
         iglob_send(2*(n_data_send_capture-1)+2,0)= pbcw
         ifind=.true.
      else
         do ineighbor = 1, nneig_send
            if (iproc.eq.pneig_send(ineighbor)) then
!$acc atomic capture
               n_data_send(ineighbor) = n_data_send(ineighbor) + 1
               n_data_send_capture    = n_data_send(ineighbor)
!$acc end atomic
               iglob_send(2*(n_data_send_capture-1)+1,ineighbor)= iglob
               iglob_send(2*(n_data_send_capture-1)+2,ineighbor)= pbcw
               ifind=.true.
               exit
            end if
         end do
      end if

!      if  (     (xr>=ixbeg).and.(xr<ixend)&
!         &.and. (yr>=iybeg).and.(yr<iyend)&
!         &.and. (zr>=izbeg).and.(zr<izend) ) then  ! (Rank domain)
!!$acc atomic capture
!         n_data_send(0) = n_data_send(0) + 1
!         n_data_send_capture = n_data_send(0)
!!$acc end atomic
!         iglob_send(2*(n_data_send_capture-1)+1,0)=iglob
!         iglob_send(2*(n_data_send_capture-1)+2,0)= pbcw
!         ifind = .true.
!      else  ! (not in rank domain)
!!$acc loop seq
!         do ineighbor = 1, nneig_send   ! Look among neighbour
!            iproc = pneig_send(ineighbor)+1
!            if( (iproc /= rank+1)&
!               &.and. (xr>=xbegproc(iproc)) .and. (xr<xendproc(iproc))&
!               &.and. (yr>=ybegproc(iproc)) .and. (yr<yendproc(iproc))&
!               &.and. (zr>=zbegproc(iproc)) .and. (zr<zendproc(iproc)) )then
!!$acc atomic capture
!               n_data_send(ineighbor) = n_data_send(ineighbor) + 1
!               n_data_send_capture = n_data_send(ineighbor)
!!$acc end atomic
!               iglob_send(2*(n_data_send_capture-1)+1,ineighbor)= iglob
!               iglob_send(2*(n_data_send_capture-1)+2,ineighbor)= pbcw
!               ifind = .true.
!               exit ! if found
!            end if
!         end do
!      end if ! (Find atom domain)

   end do
!$acc update host(n_data_send) async
!Wait for n_data_send
!$acc wait

   ! Allocate and fill data to send (Optimize reallocation)
   s_bufi = maxval(n_data_send(1:nneig_send))
   if (s_bufi.gt.max_atoms_send) then
      max_atoms_send =&
         &merge( int(s_bufi*(1+mem_inc)),s_bufi,extra_alloc )
      if (debMem.ne.0) print*,'reassign::realloc s',s_bufi,rank
      if (allocated(data_send)) then
!$acc exit data delete(data_send)
         deallocate(data_send)
      end if
      allocate( data_send( max_atoms_send,nneig_send ) )
!$acc enter data create(data_send)
   end if

!$acc parallel loop present(iglob_send) async &
!$acc         vector_length(512)
   do ineighbor = 1, nneig_send  ! Gather Data for each neighbor domain
!$acc loop vector
      do i = 1, n_data_send(ineighbor)
         iglob = iglob_send(2*(i-1)+1,ineighbor)
         d    => data_send( i, ineighbor )
         d%x  = x(iglob)
         d%y  = y(iglob)
         d%z  = z(iglob)
         d%vx = v(1,iglob)
         d%vy = v(2,iglob)
         d%vz = v(3,iglob)
         d%ax = a(1,iglob)
         d%ay = a(2,iglob)
         d%az = a(3,iglob)
      end do
   end do

   ! Initiate send data
   req_iglob_send(:) = MPI_REQUEST_NULL
   req_data_send (:) = MPI_REQUEST_NULL
!Wait for iglob_send and data_send
!$acc wait
!$acc host_data use_device(iglob_send)
   do ineighbor = 1, nneig_send
      call MPI_Isend(iglob_send(1,ineighbor)&
         &,2*n_data_send(ineighbor)&
         &,MPI_INT,   pneig_send(ineighbor), 0, get_domain_comm(),&
         &req_iglob_send(ineighbor), ierr)
   end do
!$acc end host_data

   ! Probe for sizes
   n_data_recv(:) = 0
   do ineighbor = 1, nneig_recep
      call MPI_Probe(pneig_recep(ineighbor), 0, get_domain_comm(),&
         &mpi_status1, ierr)
      call MPI_Get_count( mpi_status1, MPI_INT,&
         &n_data_recv(ineighbor), ierr )
   end do

   ! Allocate and recv data depending on message size
   n_data_recv(:)= n_data_recv(:)/2
   max_data_recv = maxval(n_data_recv)
   s_bufi        = 2*max_data_recv*nneig_recep
   call domain_workspace_i(2,s_bufi,pool_i2)
   iglob_recv(1:2*max_data_recv,1:nneig_recep)=>pool_i2(1:s_bufi)
   if ( max_data_recv.gt.max_atoms_recv ) then
      max_atoms_recv = merge(int(max_data_recv*(1+mem_inc))&
         &,max_data_recv,extra_alloc)
      if (debMem.ne.0) print*,'reassign::realloc r',max_data_recv&
         &,rank
      if (allocated(data_recv)) then
!$acc exit data delete(data_recv)
         deallocate(data_recv)
      end if
      allocate( data_recv( max_atoms_recv, nneig_recep) )
!$acc enter data create(data_recv)
   end if

   req_iglob_recv(:) = MPI_REQUEST_NULL
   req_data_recv(:)  = MPI_REQUEST_NULL
!Wait for iglob_recv and data_recv
!$acc wait
!$acc host_data use_device(iglob_recv,data_recv,data_send)
   do ineighbor = 1, nneig_recep
      call MPI_Irecv(iglob_recv(1,ineighbor),2*n_data_recv(ineighbor)&
         &,MPI_INT, pneig_recep(ineighbor), 0, get_domain_comm(),&
         &req_iglob_recv(ineighbor), ierr)
      call MPI_Irecv(data_recv(1,ineighbor),9*n_data_recv(ineighbor)&
         &,MPI_RPREC, pneig_recep(ineighbor), 1, get_domain_comm(),&
         &req_data_recv(ineighbor), ierr)
   end do

   ! Wait for iglob exchange before sending datas
   ! Overlaping those comms like before disrupts mpi with
   ! multi-node running
   call MPI_Waitall(nneig_recep,req_iglob_recv,MPI_STATUSES_IGNORE,&
      &ierr)
   call MPI_Waitall(nneig_send ,req_iglob_send,MPI_STATUSES_IGNORE,&
      &ierr)

   do ineighbor = 1, nneig_send
      call MPI_Isend(data_send(1,ineighbor),9*n_data_send(ineighbor)&
         &,MPI_RPREC, pneig_send(ineighbor), 1, get_domain_comm(),&
         &req_data_send(ineighbor), ierr)
   end do
!$acc end host_data

   nloc_save = nloc
   ! Rebuild glob(:) with old local atoms
   nloc = n_data_send(0)

   if (allocated(repartrec)) then
!$acc parallel loop present(iglob_send) async
      do i = 1, n
         if (i.lt.nloc+1)&
            &glob(i) = iglob_send(2*(i-1)+1, 0)
         loc(i)     = 0
         repart(i)  =-1
         locrec(i)  = 0
         repartrec(i) = -1
      end do
   else
!$acc parallel loop present(iglob_send) async
      do i = 1, n
         if (i.lt.nloc+1)&
            &glob(i) = iglob_send(2*(i-1)+1, 0)
         loc(i)     = 0
         repart(i)  =-1
         locrec(i)  = 0
      end do
   end if

   ! Wait all communications
   call MPI_Waitall(nneig_recep,req_data_recv ,MPI_STATUSES_IGNORE,&
      &ierr)
   call MPI_Waitall(nneig_send ,req_data_send ,MPI_STATUSES_IGNORE,&
      &ierr)

   ! Add new atoms to glob(:) and add data
!$acc parallel loop vector_length(512) &
!$acc         present(iglob_recv,nloc) copyin(n_data_recv) copy(nloc)
   do ineighbor = 1, nneig_recep
!$acc loop vector
      do i = 1, n_data_recv(ineighbor)
         iglob = iglob_recv(2*(i-1)+1,ineighbor)
         pbcw  = iglob_recv(2*(i-1)+2,ineighbor)
!$acc atomic capture
         nloc = nloc+1
         nloc_capture = nloc
!$acc end atomic
         glob(nloc_capture) = iglob
         d => data_recv(i,ineighbor)
         pbcWrap(iglob)     = pbcw
         x(iglob)   = d%x
         y(iglob)   = d%y
         z(iglob)   = d%z
         v(1,iglob) = d%vx
         v(2,iglob) = d%vy
         v(3,iglob) = d%vz
         a(1,iglob) = d%ax
         a(2,iglob) = d%ay
         a(3,iglob) = d%az
      end do
   end do

   if (use_rattle) then
      ! Recast workspace
      void_p = c_loc(data_send)
      call c_f_pointer(void_p,old_send,[max_atoms_send,nneig_send])
      void_p = c_loc(data_recv)
      call c_f_pointer(void_p,old_recv,[max_atoms_recv,nneig_recep])

!$acc parallel loop present(iglob_send,old_send,xold,yold,zold) async &
!$acc         vector_length(512)
      do ineighbor = 1, nneig_send  ! Gather Data for each neighbor domain
!$acc loop vector
         do i = 1, n_data_send(ineighbor)
            iglob = iglob_send(2*(i-1)+1,ineighbor)
            b    =>   old_send(i,ineighbor)
            b%x  = xold(iglob)
            b%y  = yold(iglob)
            b%z  = zold(iglob)
         end do
      end do

!$acc wait
!$acc host_data use_device(old_send, old_recv)
      do ineighbor = 1, nneig_send
         call MPI_Isend(old_send(1,ineighbor),3*n_data_send(ineighbor)&
            &,MPI_TPREC, pneig_send(ineighbor), 0, get_domain_comm(),&
            &req_data_send(ineighbor), ierr)
      end do
      do ineighbor = 1, nneig_recep
         call MPI_Irecv(old_recv(1,ineighbor),3*n_data_recv(ineighbor)&
            &,MPI_TPREC, pneig_recep(ineighbor), 0, get_domain_comm(),&
            &req_data_recv(ineighbor), ierr)
      end do
!$acc end host_data

      ! Wait all communications
      call MPI_Waitall(nneig_recep,req_data_recv ,MPI_STATUSES_IGNORE,&
         &ierr)
      call MPI_Waitall(nneig_send ,req_data_send ,MPI_STATUSES_IGNORE,&
         &ierr)

!$acc parallel loop vector_length(512) &
!$acc         present(iglob_recv) copyin(n_data_recv) copy(nloc)
      do ineighbor = 1, nneig_recep
!$acc loop vector
         do i = 1, n_data_recv(ineighbor)
            iglob   = iglob_recv(2*(i-1)+1,ineighbor)
            b       =>  old_recv(i,ineighbor)
            xold(iglob) = b%x
            yold(iglob) = b%y
            zold(iglob) = b%z
         end do
      end do

   end if  ! (use_rattle)

   if (mts .and. use_colvars) then
   !! send also aalt and aalt2 to have complete forces (including slow, int and fast) for ABF
     
      ! Recast workspace
      void_p = c_loc(data_send)
      call c_f_pointer(void_p,aalt_send,[max_atoms_send,nneig_send])
      void_p = c_loc(data_recv)
      call c_f_pointer(void_p,aalt_recv,[max_atoms_recv,nneig_recep])

!$acc parallel loop present(iglob_send,aalt_send,aalt,aalt2) async &
!$acc         vector_length(512)
      do ineighbor = 1, nneig_send  ! Gather Data for each neighbor domain
!$acc loop vector
         do i = 1, n_data_send(ineighbor)
            iglob = iglob_send(2*(i-1)+1,ineighbor)
            paa    =>   aalt_send(i,ineighbor)
            paa%x  = real(aalt(1,iglob),t_p)
            paa%y  = real(aalt(2,iglob),t_p)
            paa%z  = real(aalt(3,iglob),t_p)
            paa%xx  = real(aalt2(1,iglob),t_p)
            paa%yy  = real(aalt2(2,iglob),t_p)
            paa%zz  = real(aalt2(3,iglob),t_p)
         end do
      end do

!$acc wait
!$acc host_data use_device(aalt_send, aalt_recv)
      do ineighbor = 1, nneig_send 
        call MPI_Isend(aalt_send(1,ineighbor),6*n_data_send (ineighbor) &
            &,MPI_TPREC, pneig_send(ineighbor), 0, get_domain_comm(), &
            & req_data_send(ineighbor), ierr)
      end do
      do ineighbor = 1, nneig_recep
        call MPI_Irecv(aalt_recv(1,ineighbor),6*n_data_recv (ineighbor)&
              & ,MPI_TPREC, pneig_recep(ineighbor), 0, get_domain_comm(),&
              & req_data_recv(ineighbor), ierr)
      end do
!$acc end host_data

       ! Wait all communications
      call MPI_Waitall(nneig_recep,req_data_recv ,MPI_STATUSES_IGNORE,ierr)
      call MPI_Waitall(nneig_send ,req_data_send ,MPI_STATUSES_IGNORE,ierr)

!$acc parallel loop vector_length(512)&
!$acc         present(iglob_recv) copyin(n_data_recv) copy(nloc)
      do ineighbor = 1, nneig_recep
!$acc loop vector
         do i = 1, n_data_recv(ineighbor)
            iglob   = iglob_recv(2*(i-1)+1,ineighbor)
            paa       =>  aalt_recv(i,ineighbor)
            aalt(1,iglob) = real(paa%x,t_p)
            aalt(2,iglob) = real(paa%y,t_p)
            aalt(3,iglob) = real(paa%z,t_p)
            aalt2(1,iglob) = real(paa%xx,t_p)
            aalt2(2,iglob) = real(paa%yy,t_p)
            aalt2(3,iglob) = real(paa%zz,t_p)
         end do
      end do

   end if ! (mts and use_colvars)

   if (btest(tinkerdebug,tindPath).and.&
      &max_data_recv.gt.max_data_recv_save) then
      max_data_recv_save = max_data_recv
16    format(' iRank',I4,'; ',I6,' atoms(max) reassigned; nloc=',I10)
      print 16,rank,max_data_recv_save,nloc
   end if
17 format(' reassign nloc ',I9,' max_data_recv',I9)
18 format('(a,I3,a,',I0,'I6,a,',I0,'I6,a,I7,a)')
19 format(I3,6F12.6)
   if (deb_Path) write(*,17) nloc,max_data_recv
   if (tinkerdebug.gt.0.and.nproc.gt.1) then
      call AtomDebRepart(ierr)
      if (ierr.ne.0) then
         write(fmt,18) nneig_send+1,nneig_recep
         write(*,19) rank,ixbeg,ixend,iybeg,iyend&
            &,izbeg,izend
         write(* ,fmt) 'r',rank, ' send(',n_data_send,&
            &') recv(',n_data_recv, ') nloc_s(',nloc_save,')'
         call emergency_save
         __TINKER_FATAL__
      end if
   end if

!$acc parallel loop
   do i = 1,nloc
      iglob         = glob(i)
      loc(iglob)    = i
      repart(iglob) = rank
   end do

   domlen(rank+1) = nloc

   if(qtb_thermostat) call reassignqtb(nloc_save,max_data_recv&
      &,nneig_recep, n_data_recv, pool_i2, pneig_recep&
      &,nneig_send , n_data_send, pool_i1, pneig_send&
      &,max_atoms_recv, max_atoms_send)

!$acc end data
!copyin(pneig_recep, pneig_send)

20 call orderbuffer_gpu(.false.)
   call timer_exit( timer_eneig,quiet_timers )
end subroutine

end submodule parallel_migration
