c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
#include "tinker_macro.h"
      module commstuffpi_inl
        ! ReassignRespa Exchange space
        type t_elt
          real(r_p) :: x,y,z,vx,vy,vz,ax,ay,az
        end type
        integer:: max_atoms_send=0
        integer:: max_atoms_recv=0
        type(t_elt),allocatable,target :: data_send(:,:), data_recv(:,:) ! Data (x,y,...) to exchange

      contains
#include "image.f.inc"
#include "convert.f.inc"
      end module

      module commstuffpi
      implicit none

      contains

      subroutine full_gradient_prepare_pi(polymer,polymer_ctr
     &             ,istep,respa)
      use beads
      use domdec
      use atomsMirror
      use cutoff, only: use_list
      implicit none
      type(POLYMER_COMM_TYPE), intent(inout) :: polymer
      type(POLYMER_COMM_TYPE), intent(inout) :: polymer_ctr
      integer, intent(in) :: istep
      logical, intent(in) :: respa
      integer :: iloc,i

      !! UPDATE DIRECT SPACE BEFORE GRADIENT CALCULATION !!
      call update_direct_space_pi(polymer)

      !! COMMUNICATE FOR GRADIENT CALCULATIONS !!
         !! PROPAGATION REPARTITION -> GRADIENT REPARTITION !!
      call comm_for_gradient(polymer,.TRUE.,send_vel=nproc>1)
      if(contract .AND. polymer_ctr%nbeads>1) then
        !! PREPARE AND COMMUNICATE CONTRACTED POLYMER !!
        call contract_polymer(polymer,polymer_ctr)
        call comm_for_gradient(polymer_ctr,.FALSE.)
        !! REASSIGN DOFs ACCORDING TO CENTROID POSITION !!
          !! -> COMMUNICATE ALSO CONTRACTED POLYMER !!
        call reassignpi(polymer,polymer_ctr)
        !! GET NECESSARY CONTRACTED DOFs FROM NEIGHBOURING PROCS !!
        call commpospi(polymer_ctr,.FALSE.,.FALSE.)
      else
        !! REASSIGN DOFs ACCORDING TO CENTROID POSITION !!
        call reassignpi(polymer)
      endif
      !! GET NECESSARY DOFs FROM NEIGHBOURING PROCS !!
      call commpospi(polymer,.TRUE.,.FALSE.)

      !! LOAD CENTROID !! 
      call load_bead(polymer,0)

      !! REASSIGN AND COMMUNICATE POSITIONS FOR RECIPROCAL SPACE !!
      call reassignpme(.false.)
      call commposrec
      if(nproc>1) then
        !SYNCHRONIZE CENTROID
!$acc parallel loop async default(present)
        DO iloc=1,nbloc
          i=glob(iloc)
          polymer%eigpos(1,i,1)=x(i)
          polymer%eigpos(2,i,1)=y(i)
          polymer%eigpos(3,i,1)=z(i)
        ENDDO
        ! REORGANIZE POSITIONS AND VELOCITIES AFTER REASSIGN
        if(respa) then
          call comm_for_normal_modes(polymer
     &       ,polymer%pos,polymer%vel,polymer%forces)
        else
          call comm_for_normal_modes(polymer
     &       ,polymer%pos,polymer%vel)
        endif
        call update_normal_modes_pi(polymer)
      endif
      
      if(.not. centroid_recip) then
        if(contract .AND. nbeads_ctr>1) then
          call commposrecpi(polymer_ctr,.FALSE.)
        else
          call commposrecpi(polymer,.FALSE.)
        endif
      endif

      !! CENTROID IS LOADED !!
      call reCast_position
      call reinitnl(istep)
      if(respa) then
        call mechanicsteprespa(istep,.FALSE.)
        call allocsteprespa(.false.)
      else
        call mechanicstep(istep)
        call allocstep
      endif
      if (use_list) call nblist(istep)

      end subroutine full_gradient_prepare_pi

      subroutine fast_gradient_prepare_pi(polymer,istep)
      use beads
      use domdec
      use atomsMirror
      use cutoff, only: use_list
      implicit none
      type(POLYMER_COMM_TYPE), intent(inout) :: polymer
      integer, intent(in) :: istep
      integer :: iloc,i

      !! UPDATE DIRECT SPACE BEFORE GRADIENT CALCULATION !!
      call update_direct_space_pi(polymer)
      
      !! COMMUNICATE FOR GRADIENT CALCULATIONS !!
        !! PROPAGATION REPARTITION -> GRADIENT REPARTITION !!
        !! (don't send velocities as there is no reassign) !!
      call comm_for_gradient(polymer,.TRUE.,.false.)
      !! GET NECESSARY DOFs FROM NEIGHBOURING PROCS (RESPA FAST) !!
      call commpospi(polymer,.TRUE.,.TRUE.)

      !! LOAD CENTROID !! 
      call load_bead(polymer,0)
      call mechanicsteprespa(istep,.true.)
      call allocsteprespa(.true.)

      end subroutine fast_gradient_prepare_pi

      subroutine reassignpi(polymer,polymer_ctr)
        use atoms       ,only : pbcWrap
        use atomsMirror ,only : x,y,z,xold,yold,zold,n
        use beads
        use bound ,only : use_bounds
        use cell
        use domdec
        use freeze,only : use_rattle
        use inform,only : deb_Path,tindPath
        use iso_c_binding
        use moldyn,only : a,v
        use commstuffpi_inl
        use potent,only : use_pmecore
        use timestat ,only: timer_enter,timer_exit,timer_eneig
     &               ,quiet_timers
        use tinTypes ,only: real3,real3_red
        use tinheader,only: ti_eps
        use tinMemory,only: prmem_request,debMem,mem_inc,extra_alloc
        use utilcomm ,only: buffMpi_i1,buffMpi_i2
        use mpi
        use sizes    ,only: tinkerdebug
        use qtb
        use utilgpu
        implicit none
        type(POLYMER_COMM_TYPE), intent(inout) :: polymer
        type(POLYMER_COMM_TYPE), intent(inout),optional :: polymer_ctr
        integer,pointer :: iglob_send(:,:) ! Global index to send to each proc
        integer,pointer :: iglob_recv(:,:) ! Global index to recieve from procs
        real(t_p) xr, yr, zr ! adjusted position
        type(t_elt),pointer :: d
        type(real3),pointer :: b
        type(real3_red),pointer :: d_ctr
        type(real3),pointer :: old_send(:,:),old_recv(:,:)
        type(c_ptr) void_p

        integer  n_data_send(0:nneig_send),   n_data_recv(nneig_recep)
        integer req_iglob_send(nneig_send),req_iglob_recv(nneig_recep)
        integer  req_data_send(nneig_send), req_data_recv(nneig_recep)
        integer  mpi_status1(MPI_STATUS_SIZE)
        type(real3_red), allocatable,target :: data_send_ctr(:,:)
     &                                       , data_recv_ctr(:,:)

        integer ibeadbeg,ibeadend, nbeadsproc, nbeadssend 
        integer ibeadbeg_ctr,ibeadend_ctr, nbeadsproc_ctr
        integer nprocloc, commloc, rankloc ! local (pme) mpi info
        integer i,k,j, ierr, iglob, iloc, iproc, ineighbor
        integer nloc_save, n_data_send_capture,max_data_recv
        integer nloc_capture
        integer s_bufi,pbcw
        integer :: max_data_recv_save=0
        character(len=40) ::fmt
        real(r_p) :: sqrtnu
        real(t_p) ixbeg,ixend,iybeg,iyend,izbeg,izend,xtemp
        logical :: send_ctr,ifind
#if  (defined(SINGLE)||defined(MIXED))
        real(t_p), parameter:: eps1= 10*ti_eps, eps2= 1d2*ti_eps
#else
        real(t_p), parameter:: eps1= 1d2*ti_eps, eps2= 1d4*ti_eps
#endif
!$acc routine(image_acc) seq

        send_ctr = present(polymer_ctr)

        if (use_pmecore) then
           nprocloc = ndir
           commloc  = comm_dir
           rankloc  = rank_bis
           if(rank >= ndir) then
              call timer_enter( timer_eneig )
           end if
        else
           nprocloc = nproc
           commloc  = COMM_TINKER
           rankloc  = rank
        end if
        if (nproc.eq.1) return
        if (use_pmecore.and.ndir.eq.1) then
           call timer_enter( timer_eneig )
           goto 20
        end if

        if (deb_Path) write(*,*) '  -- reassignpi ',send_ctr
        call timer_enter( timer_eneig )

        s_bufi = 2*nloc*(nneig_send+1)
        call prmem_request( buffMpi_i1,max(s_bufi,size(buffMpi_i1)) )
        !Associate iglob_send pointer to buffer
        iglob_send(1:2*nloc,0:nneig_send) => buffMpi_i1(1:s_bufi)

!$acc wait
!$acc data copyin(pneig_recep, pneig_send) async
!$acc&     create(n_data_send)
!$acc&   present(xbegproc,xendproc,ybegproc,yendproc,zbegproc,zendproc)
!$acc&     present(repart,loc,glob,iglob_send)

        ! Gather process domain of rank
        ixbeg = xbegproc(rank+1)
        ixend = xendproc(rank+1)
        iybeg = ybegproc(rank+1)
        iyend = yendproc(rank+1)
        izbeg = zbegproc(rank+1)
        izend = zendproc(rank+1)

        ibeadbeg=polymer%ibead_beg(rank_polymer+1)
        ibeadend=polymer%ibead_end(rank_polymer+1)
        nbeadsproc = ibeadend-ibeadbeg+1

        if(send_ctr) then
          ibeadbeg_ctr=polymer_ctr%ibead_beg(rank_polymer+1)
          ibeadend_ctr=polymer_ctr%ibead_end(rank_polymer+1)
          nbeadsproc_ctr = ibeadend_ctr-ibeadbeg_ctr+1
        else
          nbeadsproc_ctr = 0
        endif

!$acc parallel loop async
        do i = 0,nneig_send 
           n_data_send(i) = 0
        end do

        ! Get process domain of each atoms
!$acc parallel loop async
        do i = 1, nloc
          iglob = glob(i)
          xr = polymer%eigpos(1,iglob,1)
          yr = polymer%eigpos(2,iglob,1)
          zr = polymer%eigpos(3,iglob,1)
          pbcw = pbcWrap(iglob)
          call image_inl(xr,yr,zr)
          ! Check box limits (SP Issue)
          xtemp = xr
          ifind=.FALSE.
          if ((xcell2-abs(xr)).lt.eps_cell) xr = xr-sign(4*eps_cell,xr)
          if ((ycell2-abs(yr)).lt.eps_cell) yr = yr-sign(4*eps_cell,yr)
          if ((zcell2-abs(zr)).lt.eps_cell) zr = zr-sign(4*eps_cell,zr)
          if  ( (xr>=ixbeg).and.(xr<ixend)
     &    .and. (yr>=iybeg).and.(yr<iyend)
     &    .and. (zr>=izbeg).and.(zr<izend) ) then  ! (Rank domain)
!$acc atomic capture
              n_data_send(0) = n_data_send(0) + 1
              n_data_send_capture = n_data_send(0)
!$acc end atomic
              iglob_send(2*(n_data_send_capture-1)+1,0) = iglob
              iglob_send(2*(n_data_send_capture-1)+2,0) = pbcw
              ifind=.TRUE.
          else  ! (not in rank domain)
!$acc loop seq
            do ineighbor = 1, nneig_send   ! Look among neighbour
              iproc = pneig_send(ineighbor)+1
              if( (iproc /= rank+1)
     &  .and. (xr>=xbegproc(iproc)) .and. (xr<xendproc(iproc))
     &  .and. (yr>=ybegproc(iproc)) .and. (yr<yendproc(iproc))
     &  .and. (zr>=zbegproc(iproc)) .and. (zr<zendproc(iproc)) ) then
!$acc atomic capture
                n_data_send(ineighbor) = n_data_send(ineighbor) + 1
                n_data_send_capture = n_data_send(ineighbor)
!$acc end atomic
                iglob_send(2*(n_data_send_capture-1)+1,ineighbor)=iglob
                iglob_send(2*(n_data_send_capture-1)+2,ineighbor)=pbcw
                ifind=.TRUE.
                exit ! if found
              endif
            enddo
          endif ! (Find atom domain)
          if (.not.ifind) then
             print*, 'Issue reassign',rank,iglob,xr,yr,zr
     &        ,polymer%eigpos(1,iglob,1)
     &        ,polymer%eigpos(2,iglob,1)
     &        ,polymer%eigpos(3,iglob,1)
     &        ,eps_cell,xcell2>xr,ycell2>yr,zcell2>zr
     &        ,xcell2,ycell2,zcell2,xcell,ycell,zcell
          end if
        enddo
!$acc update host(n_data_send) async
!Wait for n_data_send
!$acc wait

        ! Allocate and fill data to send (Optimzed reallocation)
        s_bufi = maxval(n_data_send(1:nneig_send))
        if (s_bufi.gt.max_atoms_send) then
          max_atoms_send =
     &            merge( int(s_bufi*(1+mem_inc)),s_bufi,extra_alloc )
          if (debMem) print*,'reassign::realloc s',s_bufi,rank
          if (allocated(data_send)) then
!$acc exit data delete(data_send) async
             deallocate(data_send)
          end if
          allocate(data_send((nbeadsproc+1)*max_atoms_send,nneig_send))
!$acc enter data create(data_send) async
        end if

        
!$acc parallel loop present(iglob_send) async
!$acc&         vector_length(512)
        do ineighbor = 1, nneig_send  ! Gather Data for each neighbor domain
!$acc loop vector
          do i = 1, n_data_send(ineighbor)
            iglob = iglob_send(2*(i-1)+1,ineighbor)
            d => data_send((i-1)*(nbeadsproc+1)+1, ineighbor )
            d%x  = polymer%eigpos(1,iglob,1)
            d%y  = polymer%eigpos(2,iglob,1)
            d%z  = polymer%eigpos(3,iglob,1)
            d%vx = polymer%eigvel(1,iglob,1)
            d%vy = polymer%eigvel(2,iglob,1)
            d%vz = polymer%eigvel(3,iglob,1)
            d%ax = 0.d0
            d%ay = 0.d0
            d%az = 0.d0
!$acc loop seq
            do k=1,nbeadsproc
              d => data_send((i-1)*(nbeadsproc+1)+k+1, ineighbor)
              d%x  = polymer%pos(1,iglob,k+ibeadbeg-1)
              d%y  = polymer%pos(2,iglob,k+ibeadbeg-1)
              d%z  = polymer%pos(3,iglob,k+ibeadbeg-1)
              d%vx = polymer%vel(1,iglob,k+ibeadbeg-1)
              d%vy = polymer%vel(2,iglob,k+ibeadbeg-1)
              d%vz = polymer%vel(3,iglob,k+ibeadbeg-1)
              d%ax = polymer%forces(1,iglob,k+ibeadbeg-1)
              d%ay = polymer%forces(2,iglob,k+ibeadbeg-1)
              d%az = polymer%forces(3,iglob,k+ibeadbeg-1)
            enddo
          end do
        end do

        ! Initiate send data
        req_iglob_send(:) = MPI_REQUEST_NULL
        req_data_send (:) = MPI_REQUEST_NULL
!Wait for iglob_send and data_send
!$acc wait
!$acc host_data use_device(iglob_send)
        do ineighbor = 1, nneig_send 
          call MPI_Isend(iglob_send(1,ineighbor)
     &            ,2*n_data_send(ineighbor)
     &            ,MPI_INT,   pneig_send(ineighbor), 0, COMM_TINKER,
     &            req_iglob_send(ineighbor), ierr)
        end do
!$acc end host_data

        ! Probe for sizes
        n_data_recv(:) = 0
        do ineighbor = 1, nneig_recep
          call MPI_Probe(pneig_recep(ineighbor), 0, COMM_TINKER,
     &                   mpi_status1, ierr)
          call MPI_Get_count( mpi_status1, MPI_INT,
     &                        n_data_recv(ineighbor), ierr )
        end do

        ! Allocate and recv data depending on message size
        n_data_recv(:)= n_data_recv(:)/2
        max_data_recv = maxval(n_data_recv)
        s_bufi        = 2*max_data_recv*nneig_recep
        call prmem_request( buffMpi_i2,max(s_bufi,size(buffMpi_i2)) )
        iglob_recv(1:2*max_data_recv,1:nneig_recep)
     &     =>buffMpi_i2(1:s_bufi)
        if ( max_data_recv.gt.max_atoms_recv ) then
          max_atoms_recv = merge(int(max_data_recv*(1+mem_inc))
     &                          ,max_data_recv,extra_alloc)
          if (debMem) print*,'reassign::realloc r',max_data_recv,rank
          if (allocated(data_recv)) then
!$acc exit data delete(data_recv)
             deallocate(data_recv)
          end if
          allocate(data_recv((nbeadsproc+1)*max_atoms_recv,nneig_recep))
!$acc enter data create(data_recv) async
        end if

        req_iglob_recv(:) = MPI_REQUEST_NULL
        req_data_recv(:)  = MPI_REQUEST_NULL
!Wait for iglob_recv and data_recv
!$acc wait
!$acc host_data use_device(iglob_recv,data_recv,data_send)
        do ineighbor = 1, nneig_recep
          call MPI_Irecv(iglob_recv(1,ineighbor)
     &            ,2*n_data_recv(ineighbor)
     &            ,MPI_INT, pneig_recep(ineighbor), 0, COMM_TINKER,
     &            req_iglob_recv(ineighbor), ierr)
          call MPI_Irecv(data_recv(1,ineighbor)
     &            ,9*n_data_recv(ineighbor)*(nbeadsproc+1)
     &            ,MPI_RPREC, pneig_recep(ineighbor), 1, COMM_TINKER,
     &            req_data_recv(ineighbor), ierr)
        end do

        ! Wait for iglob exchange before sending datas
        ! Overlaping those comms like before disrupts mpi with
        ! multi-node running
        call MPI_Waitall(nneig_recep,req_iglob_recv,MPI_STATUSES_IGNORE,
     &                   ierr)
        call MPI_Waitall(nneig_send ,req_iglob_send,MPI_STATUSES_IGNORE,
     &                   ierr)

        do ineighbor = 1, nneig_send 
          call MPI_Isend(data_send(1,ineighbor)
     &            ,9*n_data_send(ineighbor)*(nbeadsproc+1)
     &            ,MPI_RPREC, pneig_send(ineighbor), 1, COMM_TINKER,
     &            req_data_send(ineighbor), ierr)
        end do
!$acc end host_data

        nloc_save = nloc
        ! Rebuild glob(:) with old local atoms
        nloc = n_data_send(0)

!$acc parallel loop present(iglob_send) async
        do i = 1, n
           if (i.lt.nloc+1)
     &        glob(i) = iglob_send(2*(i-1)+1, 0)
           loc(i)    = 0
           repart(i) = -1
           locrec(i)    = 0
           repartrec(i) = -1
        end do

        ! Wait all communications
        call MPI_Waitall(nneig_recep,req_data_recv ,MPI_STATUSES_IGNORE,
     &                   ierr)
        call MPI_Waitall(nneig_send ,req_data_send ,MPI_STATUSES_IGNORE,
     &                   ierr)

        sqrtnu = sqrt(real(nbeads,r_p))

        ! Add new atoms to glob(:) and add data
!$acc parallel loop vector_length(512)
!$acc&         present(iglob_recv,nloc) copyin(n_data_recv) copy(nloc)
        do ineighbor = 1, nneig_recep
!$acc loop vector
          do i = 1, n_data_recv(ineighbor)
            iglob = iglob_recv(2*(i-1)+1,ineighbor)
            pbcwrap(iglob) = iglob_recv(2*(i-1)+2,ineighbor)
!$acc atomic capture
            nloc = nloc+1
            nloc_capture = nloc
!$acc end atomic
            glob(nloc_capture) = iglob
            d => data_recv((i-1)*(nbeadsproc+1)+1,ineighbor)
            polymer%eigpos(1,iglob,1) = d%x
            polymer%eigpos(2,iglob,1) = d%y
            polymer%eigpos(3,iglob,1) = d%z
            polymer%eigvel(1,iglob,1) = d%vx
            polymer%eigvel(2,iglob,1) = d%vy
            polymer%eigvel(3,iglob,1) = d%vz
!$acc loop seq
            do k=1,nbeadsproc
              d => data_recv((i-1)*(nbeadsproc+1)+k+1,ineighbor)
              polymer%pos(1,iglob,k+ibeadbeg-1)   = d%x
              polymer%pos(2,iglob,k+ibeadbeg-1)   = d%y
              polymer%pos(3,iglob,k+ibeadbeg-1)   = d%z
              polymer%vel(1,iglob,k+ibeadbeg-1)  = d%vx
              polymer%vel(2,iglob,k+ibeadbeg-1)  = d%vy
              polymer%vel(3,iglob,k+ibeadbeg-1)  = d%vz
              polymer%forces(1,iglob,k+ibeadbeg-1)   = d%ax
              polymer%forces(2,iglob,k+ibeadbeg-1)   = d%ay
              polymer%forces(3,iglob,k+ibeadbeg-1)   = d%az
            enddo
          end do
        end do
      ! sort glob to ensure the same ordering across COMM_POLYMER
#ifdef _OPENACC
      call dev_sort(nloc,glob,rec_stream)
#endif
        if (btest(tinkerdebug,tindPath).and.
     &     max_data_recv.gt.max_data_recv_save) then
           max_data_recv_save = max_data_recv
 16     format('iRank',I4,'; ',I6,' atoms(max) reassigned; nloc=',I10)
           print 16,rank,max_data_recv_save,nloc
        end if
 17     format('reassign nloc ',I9,' max_data_recv',I9)
 18     format('(a,I3,a,',I0,'I6,a,',I0,'I6,a,I7,")")')
 19     format(I3,12F11.6)
        if (deb_Path) write(*,17) nloc,max_data_recv
        if (tinkerdebug.gt.0.and.nproc.gt.1) then
           call AtomDebRepart(ierr)
           if (ierr.ne.0) then
              write(fmt,18) nneig_send+1,nneig_recep
              write(*,19) rank,xbegproc,xendproc,ybegproc
     &                   ,yendproc,zbegproc,zendproc
              write(* ,fmt), 'r',rank, ' send(',n_data_send,
     &           ') recv(',n_data_recv, ') nloc_s(',nloc_save
           end if
        end if

!$acc parallel loop async
        do i = 1,nloc
          iglob         = glob(i)
          loc(iglob)    = i
          repart(iglob) = rank
        end do

        domlen(rank+1) = nloc


        if(send_ctr) then
        ! SEND AND RECEIVE DATA FOR CONTRACTED POLYMER
          allocate(data_send_ctr(nbeadsproc_ctr*max_atoms_send 
     &        ,nneig_send))
          allocate(data_recv_ctr(nbeadsproc_ctr*max_atoms_recv
     &        ,nneig_recep))
!$acc data create(data_send_ctr,data_recv_ctr) async

          req_data_recv(:)  = MPI_REQUEST_NULL
          req_data_send(:)  = MPI_REQUEST_NULL
!$acc wait
!$acc host_data use_device(data_recv_ctr)
          do ineighbor = 1, nneig_recep
            call MPI_Irecv(data_recv_ctr(1,ineighbor)
     &              ,3*n_data_recv(ineighbor)*nbeadsproc_ctr
     &              ,MPI_RPREC, pneig_recep(ineighbor), 1, COMM_TINKER,
     &              req_data_recv(ineighbor), ierr)
          end do
!$acc end host_data


!$acc parallel loop present(iglob_send,polymer_ctr%pos) async
!$acc&         vector_length(512)
          do ineighbor = 1, nneig_send  ! Gather Data for each neighbor domain
!$acc loop vector
            do i = 1, n_data_send(ineighbor)
              iglob = iglob_send(2*(i-1)+1,ineighbor)
!$acc loop seq
              do k=1,nbeadsproc_ctr
                d_ctr => data_send_ctr((i-1)*nbeadsproc_ctr+k
     &                     , ineighbor)
                d_ctr%x  = polymer_ctr%pos(1,iglob,k+ibeadbeg_ctr-1)
                d_ctr%y  = polymer_ctr%pos(2,iglob,k+ibeadbeg_ctr-1)
                d_ctr%z  = polymer_ctr%pos(3,iglob,k+ibeadbeg_ctr-1)
              enddo
            end do
          end do


!$acc wait
!$acc host_data use_device(data_send_ctr)
          do ineighbor = 1, nneig_send 
            call MPI_Isend(data_send_ctr(1,ineighbor)
     &              ,3*n_data_send(ineighbor)*nbeadsproc_ctr
     &              ,MPI_RPREC, pneig_send(ineighbor), 1, COMM_TINKER,
     &              req_data_send(ineighbor), ierr)
          end do
!$acc end host_data

          call MPI_Waitall(nneig_recep,req_data_recv 
     &                     ,MPI_STATUSES_IGNORE,ierr)
          call MPI_Waitall(nneig_send ,req_data_send 
     &                  ,MPI_STATUSES_IGNORE, ierr)

!$acc parallel loop vector_length(512)
!$acc&         present(iglob_recv) copyin(n_data_recv)
          do ineighbor = 1, nneig_recep
!$acc loop vector
            do i = 1, n_data_recv(ineighbor)
              iglob = iglob_recv(2*(i-1)+1,ineighbor)
!$acc loop seq
              do k=1,nbeadsproc_ctr
                d_ctr => data_recv_ctr((i-1)*nbeadsproc_ctr+k
     &             ,ineighbor)
                polymer_ctr%pos(1,iglob,k+ibeadbeg_ctr-1)   = d_ctr%x
                polymer_ctr%pos(2,iglob,k+ibeadbeg_ctr-1)   = d_ctr%y
                polymer_ctr%pos(3,iglob,k+ibeadbeg_ctr-1)   = d_ctr%z
              enddo
            end do
          end do

!$acc end data
        endif

        if(piqtb) call reassignqtb(nloc_save,max_data_recv
     &      ,nneig_recep, n_data_recv, buffMpi_i2, pneig_recep 
     &      ,nneig_send , n_data_send, buffMpi_i1, pneig_send 
     &      ,max_atoms_recv, max_atoms_send)

!$acc end data

20      call orderbuffer_gpu(.false.)
        call update_nlocpi(nloc)

        call timer_exit( timer_eneig,quiet_timers )
        
      end subroutine reassignpi


      subroutine commpospi(polymer,send_centroid,fast)
      use atomsMirror
      use beads
      use domdec
      use inform   ,only: deb_Path
      use freeze
      use mpi
      use potent
      use timestat ,only: timer_enter,timer_exit,timer_commstep
     &             ,quiet_timers
      use tinMemory,only: prmem_requestm
      use utilcomm ,only: buffMpi_p1,buffMpi_p2
      use beads
      implicit none
      type(POLYMER_COMM_TYPE), intent(inout) :: polymer
      logical, intent(in) :: send_centroid
      logical, intent(in) :: fast
      integer i,k,j,tag,ierr,iproc
      integer iloc,iglob,idomlen,ibufbeg
      real(r_p),pointer:: buffer(:,:,:),buffers(:,:,:)
      integer, allocatable :: reqrec(:),reqsend(:)
      integer,dimension(nproc):: ptemp_recep,ptemp_send
      integer status(MPI_STATUS_SIZE)
      integer ntemp_recep,ntemp_send 
      integer nbeadsproc,ibeadbeg,ibeadend,nbeadssend 
      integer cshift   

      if (nproc.eq.1.or.(use_pmecore).and.(rank.gt.ndir-1)) return

      if (fast) then
        ptemp_recep = pneig_recep
        ptemp_send  = pneig_send 
        ntemp_recep = nneig_recep
        ntemp_send  = nneig_send 
      else
        ptemp_recep = pbig_recep
        ptemp_send  = pbig_send 
        ntemp_recep = nbig_recep
        ntemp_send  = nbig_send 
      end if   

      if(send_centroid) then
        cshift=1
      else
        cshift=0
      endif

      ibeadbeg=polymer%ibead_beg(rank_polymer+1)
      ibeadend=polymer%ibead_end(rank_polymer+1)
      nbeadsproc = ibeadend-ibeadbeg+1
      nbeadssend = nbeadsproc+cshift
c
c
c     communicate positions
c
      call timer_enter( timer_commstep )
      if (deb_Path) write(*,*) '   >> commpos'
 
      !Fetch memory from pool
      allocate (reqsend(nproc))
      allocate (reqrec(nproc))
      call prmem_requestm(buffMpi_p1,3*nbloc*nbeadssend)
      call prmem_requestm(buffMpi_p2,3*nloc*nbeadssend )
      buffer (1:nbeadssend,1:3,1:nbloc)
     &   =>buffMpi_p1(1:3*nbloc*nbeadssend)
      buffers(1:nbeadssend,1:3,1:nloc)
     &   =>buffMpi_p2(1:3*nloc*nbeadssend)
c
c     MPI : begin reception in buffer
c
!$acc host_data use_device(buffer)
      do i = 1, ntemp_recep
        tag = nproc*rank + ptemp_recep(i) + 1
        call MPI_IRECV(buffer(1,1,bufbeg(ptemp_recep(i)+1)),
     $   3*domlen(ptemp_recep(i)+1)*nbeadssend,
     $   MPI_RPREC,ptemp_recep(i),tag,COMM_TINKER,reqrec(i),ierr)
      end do
!$acc end host_data
c
c     MPI : move in buffer
c
!$acc parallel loop async default(present) collapse(3)
      do i = 1, nloc; do j=1,3; do k=1,nbeadsproc
        iglob = glob(i)
        if(k==cshift) then
          buffers(1,j,i) = polymer%eigpos(j,iglob,1)
        endif
        buffers(k+cshift,j,i) = polymer%pos(j,iglob,k+ibeadbeg-1)
      end do ; enddo ; enddo
c
c     MPI : begin sending
c
!$acc wait
!$acc host_data use_device(buffers)
      do i = 1, ntemp_send 
        tag = nproc*ptemp_send(i) + rank + 1
        call MPI_ISEND(buffers,3*nloc*nbeadssend 
     $       ,MPI_RPREC,
     $       ptemp_send(i),tag,COMM_TINKER,
     $       reqsend(i),ierr)
      end do
!$acc end host_data
c
      do i = 1, ntemp_recep
        call MPI_WAIT(reqrec(i),status,ierr)
      end do
      do i = 1, ntemp_send 
        call MPI_WAIT(reqsend(i),status,ierr)
      end do
c
c     MPI : move in global arrays
c
      do iproc = 1, ntemp_recep
         idomlen = domlen(ptemp_recep(iproc)+1)
         ibufbeg = bufbeg(ptemp_recep(iproc)+1)
!$acc parallel loop async default(present)  collapse(2)
         do i = 1, idomlen; do j=1,3
           iloc   = ibufbeg+i-1
           iglob  = glob(iloc)
           if(cshift==1) then
            polymer%eigpos(j,iglob,1)=buffer(1,j,iloc)
           endif
!$acc loop seq
            do k = 1, nbeadsproc
              polymer%pos(j,iglob,k+ibeadbeg-1) 
     &           = buffer(k+cshift,j,iloc)
            enddo
         end do; enddo
      end do
c
      nullify(buffer)
      nullify(buffers)
      deallocate (reqsend)
      deallocate (reqrec)

      if (deb_Path) write(*,*) '   << commpos'
      call timer_exit( timer_commstep,quiet_timers )
      end subroutine commpospi

      subroutine commposrecpi(polymer,send_centroid)
      use atomsMirror
      use domdec
      use inform   ,only: deb_Path
      use mpi
      use potent   ,only: use_pmecore
      use pme      ,only: GridDecomp1d
      use timestat ,only: timer_enter,timer_exit,timer_commstep
     &             ,quiet_timers
      use tinMemory,only: prmem_requestm
      use utilcomm ,only: buffMpi_p1,buffMpi_p2
     &             ,reqsend=>reqs_recdir,reqrec=>reqr_recdir
      use beads
      implicit none
      type(POLYMER_COMM_TYPE), intent(inout) :: polymer
      logical, intent(in) :: send_centroid
      integer i,iproc,iprec
      integer iglob,iloc,j,k
      integer tag,ierr,ibufbeg,idomlen
      integer rankloc,commloc,nprocloc
      integer status(MPI_STATUS_SIZE)
      real(r_p),pointer :: buffer(:,:,:),buffers(:,:,:)
      integer nbeadsproc,ibeadbeg,ibeadend,nbeadssend 
      integer cshift      

      if (nproc.eq.1) return
      if(send_centroid) then
        cshift=1
      else
        cshift=0
      endif

c
      if (GridDecomp1d.and.Bdecomp1d.and..not.use_pmecore) return

      call timer_enter( timer_commstep )
      if (deb_Path) write(*,*) '   >> commposrecpi'

      ibeadbeg=polymer%ibead_beg(rank_polymer+1)
      ibeadend=polymer%ibead_end(rank_polymer+1)
      nbeadsproc = ibeadend-ibeadbeg+1
      nbeadssend = nbeadsproc+cshift
c
      if (use_pmecore) then
        nprocloc = nrec
        commloc  = comm_rec
        rankloc  = rank_bis
      else
        nprocloc = nproc
        commloc  = COMM_TINKER
        rankloc  = rank
      end if
c
c     allocate some arrays
c
      call prmem_requestm(buffMpi_p1
     &   ,3*max(nblocrecdir,nlocrec2)*nbeadssend )
      call prmem_requestm(buffMpi_p2
     &   ,3*max(nloc,nlocrec)*nbeadssend)
      buffer (1:nbeadssend,1:3,1:nblocrecdir)
     &   =>buffMpi_p1(1:3*nblocrecdir*nbeadssend)
      buffers(1:nbeadssend,1:3,1:nloc)   
     &   =>buffMpi_p2(1:3*nloc*nbeadssend)
c
c     first do dir rec communications
c
c
c     MPI : begin reception in buffer
c
!$acc data present(x,y,z,glob,globrec)
!$acc data present(buffer,buffers)

!$acc host_data use_device(buffer,buffers)
      do i = 1, nrecdir_recep2
        if (precdir_recep2(i).ne.rank) then
          tag = nproc*rank + precdir_recep2(i) + 1
          call MPI_IRECV(buffer(1,1,bufbeg(precdir_recep2(i)+1)),
     $         3*domlen(precdir_recep2(i)+1)*nbeadssend,MPI_RPREC,
     $         precdir_recep2(i),tag,COMM_TINKER,reqrec(i),ierr)
        end if
      end do
c
c     MPI : move in buffer
c
!$acc parallel loop async default(present) collapse(3)
      do i = 1, nloc; do j=1,3; do k=1,nbeadsproc
        iglob = glob(i)
        if(k==cshift) then
          buffers(1,j,i) = polymer%eigpos(j,iglob,1)
        endif
        buffers(k+cshift,j,i) = polymer%pos(j,iglob,k+ibeadbeg-1)
      end do ; enddo ; enddo
c
c     MPI : begin sending
c
!$acc wait
      do i = 1, nrecdir_send2
        if (precdir_send2(i).ne.rank) then
          tag = nproc*precdir_send2(i) + rank + 1
          call MPI_ISEND(buffers,3*nloc*nbeadssend 
     &       ,MPI_RPREC,precdir_send2(i),
     &         tag,COMM_TINKER,reqsend(i),ierr)
        end if
      end do
!$acc end host_data
c
c     MPI : move in global arrays
c
      do iproc = 1, nrecdir_recep2
        if (precdir_recep2(iproc).ne.rank) then
          call MPI_WAIT(reqrec(iproc),status,ierr)
          ibufbeg = bufbeg(precdir_recep2(iproc)+1)
          idomlen = domlen(precdir_recep2(iproc)+1)
!$acc parallel loop async collapse(3) default(present)
          do i = 1, idomlen; do j=1,3; do k=1,nbeadsproc
            iloc     = ibufbeg+i-1
            iglob    = glob(iloc)
            if(k==cshift) then
              polymer%eigpos(j,iglob,1) = buffer(1,j,iloc)
            endif
            polymer%pos(j,iglob,k+ibeadbeg-1) = buffer(k+cshift,j,iloc) 
          enddo; enddo; enddo
        end if
      end do
c
      do i = 1, nrecdir_send2
         if (precdir_send2(i).ne.rank) then
            call MPI_WAIT(reqsend(i),status,ierr)
         end if
      end do
c
!$acc end data
!$acc wait
      nullify(buffer)
      nullify(buffers)
c
c     then do rec rec communications
c
      buffer (1:nbeadssend,1:3,1:nlocrec2)
     &      =>buffMpi_p1(1:3*nlocrec2*nbeadssend)
      buffers(1:nbeadssend,1:3,1:nlocrec)  
     &      =>buffMpi_p2(1:3*nlocrec*nbeadssend)

!$acc data present(buffer,buffers)
c
c     MPI : begin reception in buffer
c
!$acc host_data use_device(buffer,buffers)
      do i = 1, nrec_recep
        tag = nprocloc*rankloc + prec_recep(i) + 1
        call MPI_IRECV(buffer(1,1,bufbegrec(prec_recep(i)+1))
     $      ,3*domlenrec(prec_recep(i)+1)*nbeadssend,MPI_RPREC
     $      ,prec_recep(i),tag,commloc,reqrec(i),ierr)
      end do
c
c     MPI : move in buffer
c
!$acc parallel loop async default(present) collapse(3)
      do i = 1, nlocrec; do j=1,3; do k=1,nbeadsproc
        iglob = globrec(i)
        if(k==cshift) then
          buffers(1,j,i) = polymer%eigpos(j,iglob,1)
        endif
        buffers(k+cshift,j,i) = polymer%pos(j,iglob,k+ibeadbeg-1)
      end do ; enddo ; enddo
c
c     MPI : begin sending
c
!$acc wait
      do i = 1, nrec_send 
         tag = nprocloc*prec_send(i) + rankloc+1
         call MPI_ISEND(buffers,3*nlocrec*nbeadssend 
     $      ,MPI_RPREC,prec_send(i)
     $       ,tag,commloc,reqsend(i),ierr)
      end do
!$acc end host_data
c
c     MPI : move in global arrays
c
      do iproc = 1, nrec_recep
         call MPI_WAIT(reqrec(iproc),status,ierr)
         ibufbeg = bufbegrec(prec_recep(iproc)+1)
         idomlen = domlenrec(prec_recep(iproc)+1)
!$acc parallel loop async collapse(3) default(present)
         do i = 1, idomlen; do j=1,3; do k=1,nbeadsproc
            iloc     = ibufbeg+i-1
            iglob    = globrec(iloc)
            if(k==cshift) then
              polymer%eigpos(j,iglob,1) = buffer(1,j,iloc)
            endif
            polymer%pos(j,iglob,k+ibeadbeg-1) = buffer(k+cshift,j,iloc) 
          enddo; enddo; enddo
      end do
      ! Wait for sending requests
      do i = 1, nrec_send
        call MPI_WAIT(reqsend(i),status,ierr)
      end do

!$acc end data
!$acc end data
      nullify(buffer)
      nullify(buffers)

      if (deb_Path) write(*,*) '   << commposrecpi'
      call timer_exit( timer_commstep,quiet_timers )
      end subroutine commposrecpi


      subroutine commforcesblocpi(derivs,nbeadsproc,fast)
      use beads
      use deriv
      use domdec
      use inform,only: deb_Path
      use mpi
      use potent
      use sizes
      use timestat
      use tinMemory ,only: prmem_requestm
      use utilcomm  ,only: buffMpi_p1,buffMpi_p2
      implicit none
      integer, intent(in) :: nbeadsproc
      real(r_p), intent(inout) :: derivs(3,nbloc,*)
      logical, intent(in) :: fast
      integer i,j,k,tag,ierr,sz,sz2,ii
      integer reqsend(nproc),reqrec(nproc)
      real(r_p), pointer :: buffer(:,:,:,:),buffers(:,:,:)
!      real(r_p) derivs(3,nbloc,nbeads)
      integer,dimension(nproc):: ptemp_recep,ptemp_send
      integer status(MPI_STATUS_SIZE)
      integer ntemp_recep,ntemp_send 

c
      if (nproc.eq.1.or.((use_pmecore).and.(rank.gt.ndir-1))) return
      if (deb_Path) write(*,*) '   >> commforcesblocpi'

      if (fast) then
        ptemp_recep = pneig_recep
        ptemp_send  = pneig_send 
        ntemp_recep = nneig_recep
        ntemp_send  = nneig_send 
      else
        ptemp_recep = pbig_recep
        ptemp_send  = pbig_send 
        ntemp_recep = nbig_recep
        ntemp_send  = nbig_send 
      end if   
c
c     allocate some arrays
c
!$acc wait
      sz = 3*max(1,nloc)*ntemp_send*nbeadsproc
      sz2 = 3*nbloc*nbeadsproc
      call prmem_requestm(buffMpi_p1,sz,async=.false.)
      call prmem_requestm(buffMpi_p2,sz2,async=.false.)
      buffer(1:nbeadsproc,1:3,1:max(1,nloc),1:ntemp_send)
     &    => buffMpi_p1(1:sz)
      buffers(1:nbeadsproc,1:3,1:nbloc) 
     &     =>buffMpi_p2(1:sz2)
!$acc enter data attach(buffer,buffers)

!$acc parallel loop collapse(3) default(present)
      do k=1,nbeadsproc; do i=1,nbloc; do j=1,3
        buffers(k,j,i) = derivs(j,i,k)
      end do ; enddo ; enddo
c
c     communicate forces
c
c     MPI : begin reception in buffer
c

!$acc host_data use_device(buffer)
      do i = 1, ntemp_send 
         tag = nproc*rank + ptemp_send(i) + 1
         call MPI_IRECV(buffer(1,1,1,i),3*nloc*nbeadsproc,
     $        MPI_RPREC,ptemp_send(i),tag,
     $        COMM_TINKER,reqrec(i),ierr)
      end do
!$acc end host_data
c
c     MPI : begin sending
c

!$acc wait
!$acc host_data use_device(buffers)
      do i = 1, ntemp_recep
         tag = nproc*ptemp_recep(i) + rank + 1
         call MPI_ISEND(buffers(1,1,bufbeg(ptemp_recep(i)+1)),
     $        3*domlen(ptemp_recep(i)+1)*nbeadsproc
     $        ,MPI_RPREC,ptemp_recep(i),tag,
     $        COMM_TINKER,reqsend(i),ierr)
      end do
!$acc end host_data
c
      do i = 1, ntemp_send 
        call MPI_WAIT(reqrec(i),status,ierr)
      end do
      do i = 1, ntemp_recep
        call MPI_WAIT(reqsend(i),status,ierr)
      end do
c
c     MPI : move in global arrays
c
!$acc parallel loop gang vector collapse(3) async
!$acc&         present(derivs,buffer)
      do k=1,nbeadsproc; do i = 1, nloc; do j = 1, 3
!$acc loop seq
        do ii = 1, ntemp_send 
          derivs(j,i,k) = derivs(j,i,k) + buffer(k,j,i,ii)
        end do
      enddo ; enddo; enddo
c
!$acc exit data detach(buffer,buffers) async
      nullify(buffer)
      nullify(buffers)

      if (deb_Path) write(*,*) '   << commforcesblocpi'
      end subroutine commforcesblocpi


      subroutine comm_for_gradient(polymer,send_centroid,send_vel)
      use beads
      use mpi
      use beads
      use domdec
      use inform, only: deb_path
      implicit none
      type(POLYMER_COMM_TYPE),intent(inout) :: polymer
      logical, intent(in) :: send_centroid
      logical, intent(in), optional :: send_vel
      integer :: ibead,k,i,j,iglob,iproc,maxsend,nsend,maxloc
      integer :: maxrecv,nproc_recv
      integer :: isend,ii,ierr,ns,nr,iloc,kk
      integer, parameter :: mptag_size=0,mptag_data=1
      integer :: reqs(2*nproc_polymer)
      LOGICAL :: flag
      real(r_p) :: sqrtnu
      integer :: nbeadslocmax, nlocproc
      integer :: ibeadbeg,ibeadend,nu
      integer :: ilocend,ilocbeg
      integer :: cshift,nelt
      integer :: send_vel_
      real(r_p), allocatable :: bufferpi_s(:,:,:)
      real(r_p), allocatable :: bufferpi_r(:,:,:)

      if(nproc_polymer==1) return
      reqs(:)=MPI_REQUEST_NULL

      if(deb_path) then
          write(*,*) '>> comm_for_gradient'
      endif
      if(send_centroid) then
        cshift=1
      else
        cshift=0
      endif

      send_vel_=.true.
      if(present(send_vel)) send_vel_=send_vel
      if(allocated(polymer%vel) .and. send_vel_) then
        nelt=2
      else
        nelt=1
      endif

      nu=polymer%nbeads
      sqrtnu=sqrt(real(nu))

      if(pi_comm_report) then
        write(pi_comm_unit,*) "  entering comm_for_gradient"
        FLUSH(pi_comm_unit)
      endif
     
      maxloc=nlocpis(nproc_polymer)
      nbeadslocmax=polymer%nbeadsloc(nproc_polymer)
      nlocproc=int(nloc/nproc_polymer)

      allocate(bufferpi_s(nelt*3,maxloc*(nbeadslocmax+cshift)
     &     ,nproc_polymer))
      allocate(bufferpi_r(nelt*3,maxloc*(nbeadslocmax+cshift)
     &     ,nproc_polymer))
      
!$acc enter data create(bufferpi_r,bufferpi_s) async

        !PUT DATA INTO BUFFER ARRAYS

      ilocbeg=ilocpi_beg(rank_polymer+1)
      ilocend=ilocpi_end(rank_polymer+1)

!$acc parallel loop collapse(2) async
!$acc& present(bufferpi_s)
        DO k=1,nu; DO iloc=0,ilocend-ilocbeg !; DO j=1,3
          iproc=polymer%iproc(k)
          if(rank_polymer==iproc-1) CYCLE
          ibeadbeg=polymer%ibead_beg(iproc)
          i=glob(ilocbeg+iloc)
          if(k-ibeadbeg==0 .and. send_centroid) then
            kk=iloc +1
!$acc loop seq
            DO j=1,3
            bufferpi_s(j,kk,iproc) = polymer%eigpos(j,i,1)
            ENDDO
            if(nelt==2) then
!$acc loop seq
              DO j=1,3
                bufferpi_s(3+j,kk,iproc) = polymer%eigvel(j,i,1)
              ENDDO
            endif
          endif
          kk=(k-ibeadbeg + cshift)*nlocpi+ iloc +1
!$acc loop seq
          DO j=1,3
          bufferpi_s(j,kk,iproc)=polymer%pos(j,i,k)
          ENDDO
          if(nelt==2) then
!$acc loop seq
            DO j=1,3
              bufferpi_s(3+j,kk,iproc)=polymer%vel(j,i,k)
            ENDDO
          endif
        ENDDO; ENDDO !; ENDDO

!$acc wait
!$acc host_data use_device(bufferpi_s,bufferpi_r)
      ii=0
      DO iproc=1,nproc_polymer
        if(rank_polymer==iproc-1) CYCLE
        ii=ii+1
        ibeadbeg = polymer%ibead_beg(iproc)
        ibeadend = polymer%ibead_end(iproc)
        if(pi_comm_report) then
          write(pi_comm_unit,*) "  send"
     &       ,nelt*3*(ilocend-ilocbeg+1)*(ibeadend-ibeadbeg+1+cshift)
     &         ,"ats to",iproc-1
        endif
        call MPI_ISEND(bufferpi_s(1,1,iproc)
     &     ,nelt*3*(ilocend-ilocbeg+1)*(ibeadend-ibeadbeg+1+cshift)
     &      ,MPI_RPREC
     &     ,iproc-1, mptag_data, COMM_POLYMER, reqs(ii),ierr)
      ENDDO

      ibeadbeg = polymer%ibead_beg(rank_polymer+1)
      ibeadend = polymer%ibead_end(rank_polymer+1)
      DO iproc=1,nproc_polymer
        if(rank_polymer==iproc-1) CYCLE
        ii=ii+1
        ilocbeg=ilocpi_beg(iproc)
        ilocend=ilocpi_end(iproc)
        if(pi_comm_report) then
          write(pi_comm_unit,*) "  receive"
     &       ,nelt*3*(ilocend-ilocbeg+1)*(ibeadend-ibeadbeg+1+cshift)
     &         ,"ats from",iproc-1
        endif
        call MPI_IRECV(bufferpi_r(1,1,iproc)
     &   ,nelt*3*(ilocend-ilocbeg+1)*(ibeadend-ibeadbeg+1+cshift)
     &   ,MPI_RPREC,iproc-1,mptag_data, COMM_POLYMER
     &   ,reqs(ii),ierr)
      ENDDO
!$acc end host_data

      !WAIT FOR COMMUNICATIONS
      if(pi_comm_report) then
        write(pi_comm_unit,*) "  Wait for data"
        FLUSH(pi_comm_unit)
      endif
      call MPI_WAITALL(ii,reqs,MPI_STATUSES_IGNORE,ierr)
      

      !PUT DATA INTO POLYMER ARRAYS

!$acc parallel loop collapse(2) async
!$acc& present(bufferpi_r,ilocpi_beg,ilocpi_end)
        DO k=1-cshift,ibeadend-ibeadbeg+1; DO iloc=1,nloc !; DO j=1,3
          iproc=min(nproc_polymer,int((iloc-1)/nlocproc)+1)
          if(rank_polymer==iproc-1) CYCLE
          ilocbeg=ilocpi_beg(iproc)
          ilocend=ilocpi_end(iproc)
          i=glob(iloc)
          if(k==0) then
            kk=iloc-ilocbeg + 1
!$acc loop seq
            DO j=1,3
            polymer%eigpos(j,i,1)=bufferpi_r(j,kk,iproc)
            ENDDO
            if(nelt==2) then
!$acc loop seq
              DO j=1,3
                polymer%eigvel(j,i,1)=bufferpi_r(3+j,kk,iproc)
              ENDDO
            endif
          else
            ibead=k+ibeadbeg-1
            kk=(k+cshift-1)*(ilocend-ilocbeg+1)+ iloc-ilocbeg +1
!$acc loop seq
            DO j=1,3
            polymer%pos(j,i,ibead)=bufferpi_r(j,kk,iproc)
            ENDDO
            if(nelt==2) then
!$acc loop seq
              DO j=1,3
                polymer%vel(j,i,ibead)=bufferpi_r(3+j,kk,iproc)
              ENDDO
            endif
          endif
        ENDDO; ENDDO !; ENDDO

c!$acc wait
c      call MPI_BARRIER(COMM_POLYMER,ierr)

!$acc exit data delete(bufferpi_s,bufferpi_r) async
      deallocate(bufferpi_s,bufferpi_r)

      if(pi_comm_report) then
        write(pi_comm_unit,*) "  exiting comm_for_gradient"
        FLUSH(pi_comm_unit)
      endif

      if(deb_path) then
          write(*,*) '<< comm_for_gradient'
      endif

      end subroutine comm_for_gradient

      subroutine comm_for_normal_modes(polymer
     &   ,array1,array2,array3,array4)
      use beads
      use mpi
      use beads
      use domdec
      use inform, only: deb_path
      implicit none
      type(POLYMER_COMM_TYPE),intent(inout) :: polymer
      real(r_p), intent(inout) :: array1(:,:,:)
      real(r_p), intent(inout), optional :: array2(:,:,:)
      real(r_p), intent(inout), optional :: array3(:,:,:)
      real(r_p), intent(inout), optional :: array4(:,:,:)
      logical :: send_pos_, send_vel_,slow_
      logical :: correctorder
      integer :: ibead,k,i,j,iglob,iproc,maxsend,nsend,maxloc
      integer :: maxrecv,nproc_recv
      integer :: isend,ii,ierr,ns,nr,iloc,kk
      integer, parameter :: mptag_size=0,mptag_data=1
      integer :: reqs(2*nproc_polymer)
      LOGICAL :: flag
      real(r_p) :: sqrtnu
      integer :: nbeadslocmax,nlocproc
      integer :: ibeadbeg,ibeadend,nu
      integer :: ilocend,ilocbeg,nelt,ielt
      real(r_p), allocatable :: bufferpi_s(:,:,:)
      real(r_p), allocatable :: bufferpi_r(:,:,:)

      if(nproc_polymer==1) return

      if(deb_path) then
!$acc wait
          write(*,*) '>> comm_for_normal_modes'
      endif

      reqs(:)=MPI_REQUEST_NULL

      nu=polymer%nbeads
      sqrtnu=sqrt(real(nu,r_p))

      if(pi_comm_report) then
        write(pi_comm_unit,*) "  entering comm_for_normal_modes"
        FLUSH(pi_comm_unit)
      endif
     
      maxloc=nlocpis(nproc_polymer)
      nbeadslocmax=polymer%nbeadsloc(nproc_polymer)
      nlocproc=int(nloc/nproc_polymer)

      nelt=1
      if(present(array2)) then
        nelt=nelt+1
      endif
      if(present(array3)) then
        if(.not. present(array2)) then
          write(*,*) 'Error in comm_for_normal_modes: array3'
     &     ,' is present but array2 is not'
          call fatal
        endif
        nelt=nelt+1
      endif
      if(present(array4)) then
        if(.not. present(array3)) then
          write(*,*) 'Error in comm_for_normal_modes: array4'
     &     ,' is present but array3 is not'
          call fatal
        endif
        nelt=nelt+1
      endif
      allocate(bufferpi_s(nelt*3,maxloc*nbeadslocmax,nproc_polymer))
      allocate(bufferpi_r(nelt*3,maxloc*nbeadslocmax,nproc_polymer))

c      allocate(bufferpi_s(3*3,maxloc*nbeadslocmax,nproc_polymer))
c      allocate(bufferpi_r(3*3,maxloc*nbeadslocmax,nproc_polymer))
      
!$acc enter data create(bufferpi_r,bufferpi_s) async

        !PUT DATA INTO BUFFER ARRAYS
      ibeadbeg = polymer%ibead_beg(rank_polymer+1)
      ibeadend = polymer%ibead_end(rank_polymer+1)
!$acc parallel loop collapse(2) async
!$acc& present(bufferpi_s,ilocpi_beg,ilocpi_end)
        DO k=0,ibeadend-ibeadbeg; DO iloc=1,nloc !; DO j=1,3
          iproc=min(nproc_polymer,int((iloc-1)/nlocproc)+1)
          if(rank_polymer==iproc-1) cycle
          ilocbeg=ilocpi_beg(iproc)
          ilocend=ilocpi_end(iproc)
          i=glob(iloc)
          ibead=k+ibeadbeg
          kk=k*(ilocend-ilocbeg+1)+ iloc-ilocbeg + 1
!$acc loop seq
          do j=1,3
          bufferpi_s(  j,kk,iproc)=array1(j,i,ibead)
          enddo
          if(nelt==1) cycle
!$acc loop seq
          do j=1,3
          bufferpi_s(3+j,kk,iproc)=array2(j,i,ibead)
          enddo
          if(nelt==2) cycle
!$acc loop seq
          do j=1,3
          bufferpi_s(6+j,kk,iproc)=array3(j,i,ibead)
          enddo
          if(nelt==3) cycle
!$acc loop seq
          do j=1,3
          bufferpi_s(9+j,kk,iproc)=array4(j,i,ibead)
          enddo
        ENDDO; ENDDO !; ENDDO

!$acc wait
!$acc host_data use_device(bufferpi_s,bufferpi_r)
      ii=0
      DO iproc=1,nproc_polymer
        if(rank_polymer==iproc-1) CYCLE
        ii=ii+1
        ilocbeg=ilocpi_beg(iproc)
        ilocend=ilocpi_end(iproc)
        if(pi_comm_report) then
          write(pi_comm_unit,*) "  send"
     &       ,nelt*3*(ilocend-ilocbeg+1)*(ibeadend-ibeadbeg+1)
     &         ,"ats to",iproc-1
        endif
        call MPI_ISEND(bufferpi_s(1,1,iproc)
     &     ,nelt*3*(ilocend-ilocbeg+1)*(ibeadend-ibeadbeg+1)
     &      ,MPI_RPREC
     &     ,iproc-1, mptag_data, COMM_POLYMER, reqs(ii),ierr)
      ENDDO

      !RECEIVE DATA
      ilocbeg=ilocpi_beg(rank_polymer+1)
      ilocend=ilocpi_end(rank_polymer+1)
      DO iproc=1,nproc_polymer
        if(rank_polymer==iproc-1) CYCLE
        ii=ii+1
        ibeadbeg = polymer%ibead_beg(iproc)
        ibeadend = polymer%ibead_end(iproc)
        if(pi_comm_report) then
          write(pi_comm_unit,*) "  receive"
     &       ,nelt*3*(ilocend-ilocbeg+1)*(ibeadend-ibeadbeg+1)
     &         ,"ats to",iproc-1
        endif
        call MPI_IRECV(bufferpi_r(1,1,iproc)
     &   ,nelt*3*(ilocend-ilocbeg+1)*(ibeadend-ibeadbeg+1)
     &   ,MPI_RPREC,iproc-1,mptag_data, COMM_POLYMER
     &   ,reqs(ii),ierr)
      ENDDO
!$acc end host_data

      !WAIT FOR COMMUNICATIONS
      call MPI_WAITALL(ii,reqs,MPI_STATUSES_IGNORE,ierr)

      !PUT DATA INTO POLYMER ARRAYS

!$acc parallel loop collapse(2) async
!$acc& present(polymer,bufferpi_r,polymer%ibead_beg,polymer%iproc)
        DO k=1,nu; DO iloc=0,ilocend-ilocbeg !; DO j=1,3
          iproc=polymer%iproc(k)
          if(rank_polymer==iproc-1) cycle
          ibeadbeg=polymer%ibead_beg(iproc)
          i=glob(ilocbeg+iloc)
          kk=(k-ibeadbeg)*nlocpi + iloc + 1
!$acc loop seq
          do j=1,3
          array1(j,i,k)=bufferpi_r(j,kk,iproc)
          enddo
          if(nelt==1) cycle
!$acc loop seq
          do j=1,3
          array2(j,i,k)=bufferpi_r(3+j,kk,iproc)
          enddo
          if(nelt==2) cycle
!$acc loop seq
          do j=1,3
          array3(j,i,k)=bufferpi_r(6+j,kk,iproc)
          enddo
          if(nelt==3) cycle
!$acc loop seq
          do j=1,3
          array4(j,i,k)=bufferpi_r(9+j,kk,iproc)
          enddo
        ENDDO; ENDDO !; ENDDO

c!$acc wait
c      call MPI_BARRIER(COMM_POLYMER,ierr)

!$acc exit data delete(bufferpi_s,bufferpi_r) async
      deallocate(bufferpi_s,bufferpi_r)

      if(pi_comm_report) then
        write(pi_comm_unit,*) "  exiting comm_for_normal_modes"
        FLUSH(pi_comm_unit)
      endif

      if(deb_path) then
          write(*,*) '<< comm_for_normal_modes'
      endif

      end subroutine comm_for_normal_modes

      subroutine save_deriv_rec(temprec)
      use potent ,pa=>PotentialAll
      use domdec
      use deriv
      implicit none
      real(r_p), intent(inout) :: temprec(3,nlocrec2)
      integer :: i,j

      if (pa) then
!$acc parallel loop collapse(2) async
!$acc&   present(demrec,deprec,decrec,dedsprec)
         do i = 1, nlocrec2
            do j = 1, 3
              temprec(j,i)= demrec(j,i) 
     &                      + deprec(j,i)
     &                      + decrec(j,i)
     &                      + dedsprec(j,i)
              demrec(j,i)=0.d0
              deprec(j,i)=0.d0
              decrec(j,i)=0.d0
              dedsprec(j,i)=0.d0
            end do
         end do
      else
         if (use_polar) then
!$acc parallel loop collapse(2) async
!$acc&      present(demrec,deprec)
            do i = 1, nlocrec2
               do j = 1, 3
                  temprec(j,i)= demrec(j,i) + deprec(j,i)
                  demrec(j,i)=0.d0
                  deprec(j,i)=0.d0
               end do
            end do
         else if (use_mpole) then
!$acc parallel loop collapse(2) async
!$acc&      present(demrec)
            do i = 1, nlocrec2
               do j = 1, 3
                  temprec(j,i)= demrec(j,i)
                  demrec(j,i)=0.d0
               end do
            end do
         end if
         if (use_charge) then
!$acc parallel loop collapse(2) async
!$acc&      present(decrec)
            do i = 1, nlocrec2
               do j = 1, 3
                  temprec(j,i)=  decrec(j,i)
                  decrec(j,i)=0.d0
               end do
            end do
         end if
         if (use_disprec) then
!$acc parallel loop collapse(2) async
!$acc&      present(decrec)
            do i = 1, nlocrec2
               do j = 1, 3
                  temprec(j,i)=temprec(j,i) + dedsprec(j,i)
                  dedsprec(j,i)=0.d0
               end do
            end do
         end if
      end if


      end subroutine save_deriv_rec

      subroutine commforcesrecpi(derivs,temprec,nbeadsproc)
      use deriv
      use domdec
      use mpi
      use potent ,pa=>PotentialAll
      use pme    ,only: GridDecomp1d
      use inform ,only: deb_Path
      use sizes
      use timestat ,only: timer_enter,timer_exit,quiet_timers
     &             ,timer_dirreccomm,timer_fmanage
      use tinMemory,only: prmem_requestm
      use utilcomm ,only: buffMpi_p1,buffMpi_p2,buffMpi_p3
      implicit none
      integer, intent(in) :: nbeadsproc
      real(r_p), intent(inout) :: derivs(3,nbloc,*)
      real(r_p), intent(inout) :: temprec(3,nlocrec2,*)
      integer i,j,k,tag,ierr,iproc,iglob,iloc,ilocrec,ibufbeg
      integer jloc,jglob,jlocrec,ii
      integer sz1,sz2,sz3
      integer rankloc,commloc,nprocloc
      integer reqsend(nproc),reqrec(nproc)
      integer status(MPI_STATUS_SIZE)
      real(r_p),pointer:: buffer(:,:,:,:)
      real(r_p),pointer:: buffers(:,:,:),temprecbuf(:,:,:)

      !case sequential
      if (nproc.eq.1) then
!$acc parallel loop collapse(3)
!$acc&         default(present) async
         do k= 1,nbeadsproc; do i = 1, nloc
            do j = 1, 3
               derivs(j,i,k) = derivs(j,i,k) + temprec(j,i,k)
            end do
         end do; end do

         return
      end if

      if (GridDecomp1d.and.Bdecomp1d.and..not.use_pmecore) then
!$acc parallel loop collapse(3)
!$acc&         default(present) async
         do k= 1,nbeadsproc; do i = 1, nbloc
            do j = 1, 3
              ilocrec = locrec(glob(i))
              if (ilocrec.eq.0) cycle
              derivs(j,i,k) = derivs(j,i,k) + temprec(j,ilocrec,k)
            end do
         end do; end do

         return
      end if

      if (deb_Path) write(*,*) '   >> commforcesrecpi'
      call timer_enter( timer_dirreccomm )
c
      if (use_pmecore) then
        nprocloc = nrec
        commloc  = comm_rec
        rankloc  = rank_bis
      else
        nprocloc = nproc
        commloc  = COMM_TINKER
        rankloc  = rank
      end if
c
c     allocate some arrays
c
      sz1=3*max(max(1,nloc)*nrecdir_send1,nlocrec*nrec_send)*nbeadsproc
      sz2=3*max(1,nblocrecdir)*nbeadsproc
      sz3=3*max(1,nlocrec2)*nbeadsproc
      call prmem_requestm(buffMpi_p1,max(size(buffMpi_p1),sz1))
      call prmem_requestm(buffMpi_p2,max(size(buffMpi_p2),sz2))
      call prmem_requestm(buffMpi_p3,max(size(buffMpi_p3),sz3))

      sz1 = 3*max(1,nlocrec)*nrec_send
      buffer(1:nbeadsproc,1:3,1:max(1,nlocrec),1:nrec_send) 
     &   => buffMpi_p1(1:sz1)
      buffers(1:nbeadsproc,1:3,1:max(1,nblocrecdir))        
     &   => buffMpi_p2(1:sz2)
      temprecbuf(1:nbeadsproc,1:3,1:max(1,nlocrec2))        
     &   => buffMpi_p3(1:sz3)
!$acc enter data attach(buffer,buffers)
c
c     first do rec rec communications
c
c
c     MPI : begin reception in buffer
c
!$acc host_data use_device(buffer)
      do i = 1, nrec_send
        tag = nprocloc*rankloc + prec_send(i) + 1
        call MPI_IRECV(buffer(1,1,1,i),3*nlocrec*nbeadsproc,MPI_RPREC
     $      ,prec_send(i),tag,commloc,reqrec(i),ierr)
      end do
!$acc end host_data
c
c     MPI : move in buffer
c
!$acc data present(temprec,temprecbuf,buffers,derivs)

!$acc parallel loop collapse(3)
!$acc&         default(present) async
         do k= 1,nbeadsproc; do i = 1, nlocrec2
            do j = 1, 3
              temprecbuf(k,j,i) = temprec(j,i,k)
            end do
         end do; end do
c
!$acc wait
!$acc host_data use_device(temprec)
      do i = 1, nrec_recep
        tag = nprocloc*prec_recep(i) + rankloc + 1
        call MPI_ISEND(temprecbuf(1,1,bufbegrec(prec_recep(i)+1))
     $      ,3*domlenrec(prec_recep(i)+1)*nbeadsproc
     $      ,MPI_RPREC
     $      ,prec_recep(i),tag,commloc,reqsend(i),ierr)
      end do
!$acc end host_data
c
      do i = 1, nrec_send
        call MPI_WAIT(reqrec(i),status,ierr)
      end do
      do i = 1, nrec_recep
        call MPI_WAIT(reqsend(i),status,ierr)
      end do
c
c     mpi : move in global arrays
c
!$acc parallel loop collapse(3) present(buffer)
      do k=1,nbeadsproc; do i = 1,nlocrec
        do j = 1,3
!$acc loop seq
          do ii = 1,nrec_send
             temprec(j,i,k) = temprec(j,i,k) + buffer(k,j,i,ii)
          end do
        end do
      end do; enddo

!$acc exit data detach(buffer)
      sz1 = 3*max(1,nloc)*nrecdir_send1*nbeadsproc
      nullify(buffer)
      buffer(1:nbeadsproc,1:3,1:max(1,nloc),1:nrecdir_send1) 
     &   => buffMpi_p1(1:sz1)
!$acc enter data attach(buffer)
c
c     then do rec dir communications
c
c     MPI : begin reception in buffer
c
!$acc host_data use_device(buffer)
      do i = 1, nrecdir_send1
         if (precdir_send1(i).ne.rank) then
           tag = nproc*rank + precdir_send1(i) + 1
           call MPI_IRECV(buffer(1,1,1,i),3*nloc*nbeadsproc,MPI_RPREC
     $         ,precdir_send1(i),tag,COMM_TINKER,reqrec(i),ierr)
         end if
      end do
!$acc end host_data
c
c     Move in array
c
      do ii = 1, nrecdir_recep1
        if (precdir_recep1(ii).ne.rank) then
          ibufbeg =bufbeg(precdir_recep1(ii)+1)-1 
!$acc parallel loop collapse(3)
!$acc&    present(repartrec,locrec,glob)
          do k=1,nbeadsproc; do i = 1, domlen(precdir_recep1(ii)+1)
            do j = 1,3
              iloc = ibufbeg+i
              iglob = glob(iloc)
              if (repartrec(iglob).eq.rankloc) then
                ilocrec = locrec(iglob)
                buffers(k,j,iloc) = temprec(j,ilocrec,k)
              else
                buffers(k,j,iloc) = 0
              end if
            end do
          end do; enddo
        end if
      end do
c
!$acc host_data use_device(buffers)
      do i = 1, nrecdir_recep1
        if (precdir_recep1(i).ne.rank) then
          tag = nproc*precdir_recep1(i) + rank + 1
          call MPI_ISEND(buffers(1,1,bufbeg(precdir_recep1(i)+1))
     $        ,3*domlen(precdir_recep1(i)+1)*nbeadsproc,MPI_RPREC
     $        ,precdir_recep1(i),tag,COMM_TINKER,reqsend(i),ierr)
        end if
      end do
!$acc end host_data
c
      do i = 1, nrecdir_send1
        if (precdir_send1(i).ne.rank) then
        call MPI_WAIT(reqrec(i),status,ierr)
        end if
      end do
c
c     MPI : move in global arrays
c
      do ii = 1, nrecdir_send1
        if (precdir_send1(ii).ne.rank) then
!$acc parallel loop collapse(3) async present(buffer)
          do k=1,nbeadsproc; do i = 1,nloc
            do j = 1,3
              derivs(j,i,k) = derivs(j,i,k) + buffer(k,j,i,ii)
            end do
          end do; enddo
        else
!$acc parallel loop collapse(3) async
          do k=1,nbeadsproc; do i = 1, nlocrec
            do j = 1,3
              iglob = globrec(i)
              iloc  = loc(iglob)
              if (repart(iglob).eq.rank) then
                 derivs(j,iloc,k) = derivs(j,iloc,k) + temprec(j,i,k)
              end if
            end do
          end do; enddo
        end if
      end do

      do i = 1, nrecdir_recep1
        if (precdir_recep1(i).ne.rank) then
        call MPI_WAIT(reqsend(i),status,ierr)
        end if
      end do

!$acc wait
!$acc end data
!$acc exit data detach(buffer,buffers,temprecbuf)
      nullify(buffer)
      nullify(buffers)
      nullify(temprecbuf)

      call timer_exit( timer_dirreccomm,quiet_timers )
      if (deb_Path) write(*,*) '   << commforcesrecpi'


      end subroutine commforcesrecpi


      subroutine commforces_ctr(polymer_ctr,derivs_ctr)
      use beads
      use mpi
      use beads
      use domdec
      implicit none
      type(POLYMER_COMM_TYPE),intent(inout) :: polymer_ctr
      real(r_p), intent(inout) :: derivs_ctr(3,nbloc,*)
      integer :: ibead,k,i,j,iglob,iproc,maxsend,nsend,maxloc
      integer :: maxrecv,nproc_recv
      integer :: isend,ii,ierr,ns,nr,iloc,kk
      integer, parameter :: mptag_size=0,mptag_data=1
      integer :: reqs(2*nproc_polymer)
      LOGICAL :: flag
      real(r_p) :: sqrtnu
      integer :: nbeadslocmax,nlocproc,nbeadsproc
      integer :: ibeadbeg,ibeadend,nu
      integer :: ilocend,ilocbeg
      real(r_p), allocatable :: bufferpi_s(:,:,:)
      real(r_p), allocatable :: bufferpi_r(:,:,:)

      nu=polymer_ctr%nbeads

      ilocbeg=ilocpi_beg(rank_polymer+1)
      ilocend=ilocpi_end(rank_polymer+1)
      ibeadbeg = polymer_ctr%ibead_beg(rank_polymer+1)
      ibeadend = polymer_ctr%ibead_end(rank_polymer+1)
      nbeadsproc = ibeadend-ibeadbeg+1

!$acc parallel loop collapse(3) async
!$acc& present(polymer,glob,derivs_ctr)
      DO k=1,nbeadsproc; DO iloc=ilocbeg,ilocend ; DO j=1,3
        i=glob(iloc)
        polymer_ctr%forces_slow(j,i,k+ibeadbeg-1)
     &     =-derivs_ctr(j,iloc,k)
      ENDDO; ENDDO ; ENDDO
        
      if(nproc_polymer==1)   return

      reqs(:)=MPI_REQUEST_NULL

      if(pi_comm_report) then
        write(pi_comm_unit,*) "  entering comm_for_normal_modes"
        FLUSH(pi_comm_unit)
      endif
     
      maxloc=nlocpis(nproc_polymer)
      nbeadslocmax=polymer_ctr%nbeadsloc(nproc_polymer)
      nlocproc=int(nloc/nproc_polymer)

      allocate(bufferpi_s(3,maxloc*nbeadslocmax,nproc_polymer))
      allocate(bufferpi_r(3,maxloc*nbeadslocmax,nproc_polymer))
      
!$acc enter data create(bufferpi_r,bufferpi_s) async

        !PUT DATA INTO BUFFER ARRAYS
      
!$acc parallel loop collapse(2) async
!$acc& present(polymer_ctr,bufferpi_s,ilocpi_beg,ilocpi_end)
        DO k=1,ibeadend-ibeadbeg+1; DO iloc=1,nloc !; DO j=1,3
          iproc=min(nproc_polymer,int((iloc-1)/nlocproc)+1)
          if(rank_polymer/=iproc-1) then
            ilocbeg=ilocpi_beg(iproc)
            ilocend=ilocpi_end(iproc)
            kk=(k-1)*(ilocend-ilocbeg+1)+ iloc-ilocbeg + 1
!$acc loop seq
            do j=1,3
            bufferpi_s(j,kk,iproc)=derivs_ctr(j,iloc,k)
            enddo
          endif
        ENDDO; ENDDO !; ENDDO

!$acc wait
!$acc host_data use_device(bufferpi_s,bufferpi_r)
      ii=0
      DO iproc=1,nproc_polymer
        if(rank_polymer==iproc-1) CYCLE
        ii=ii+1
        ilocbeg=ilocpi_beg(iproc)
        ilocend=ilocpi_end(iproc)
        if(pi_comm_report) then
          write(pi_comm_unit,*) "  send"
     &       ,3*(ilocend-ilocbeg+1)*(ibeadend-ibeadbeg+1)
     &         ,"ats to",iproc-1
        endif
        call MPI_ISEND(bufferpi_s(1,1,iproc)
     &     ,3*(ilocend-ilocbeg+1)*(ibeadend-ibeadbeg+1)
     &      ,MPI_RPREC
     &     ,iproc-1, mptag_data, COMM_POLYMER, reqs(ii),ierr)
      ENDDO

      !RECEIVE DATA
      ilocbeg=ilocpi_beg(rank_polymer+1)
      ilocend=ilocpi_end(rank_polymer+1)
      DO iproc=1,nproc_polymer
        if(rank_polymer==iproc-1) CYCLE
        ii=ii+1
        ibeadbeg = polymer_ctr%ibead_beg(iproc)
        ibeadend = polymer_ctr%ibead_end(iproc)
        if(pi_comm_report) then
          write(pi_comm_unit,*) "  receive"
     &       ,3*(ilocend-ilocbeg+1)*(ibeadend-ibeadbeg+1)
     &         ,"ats to",iproc-1
        endif
        call MPI_IRECV(bufferpi_r(1,1,iproc)
     &   ,3*(ilocend-ilocbeg+1)*(ibeadend-ibeadbeg+1)
     &   ,MPI_RPREC,iproc-1,mptag_data, COMM_POLYMER
     &   ,reqs(ii),ierr)
      ENDDO
!$acc end host_data

      !WAIT FOR COMMUNICATIONS
      call MPI_WAITALL(ii,reqs,MPI_STATUSES_IGNORE,ierr)

      !PUT DATA INTO POLYMER ARRAYS

!$acc parallel loop collapse(2) async
!$acc& present(polymer,bufferpi_r,polymer%ibead_beg,polymer%iproc)
        DO k=1,nu; DO iloc=0,ilocend-ilocbeg !; DO j=1,3
          iproc=polymer_ctr%iproc(k)
          if(rank_polymer/=iproc-1) then
            ibeadbeg=polymer_ctr%ibead_beg(iproc)
            i=glob(ilocbeg+iloc)
            kk=(k-ibeadbeg)*nlocpi + iloc + 1
!$acc loop seq
            do j=1,3
            polymer_ctr%forces_slow(j,i,k)=-bufferpi_r(j,kk,iproc)
            enddo
          endif
        ENDDO; ENDDO !; ENDDO

c!$acc wait
c      call MPI_BARRIER(COMM_POLYMER,ierr)

!$acc exit data delete(bufferpi_s,bufferpi_r) async
      deallocate(bufferpi_s,bufferpi_r)

      if(pi_comm_report) then
        write(pi_comm_unit,*) "  exiting comm_for_normal_modes"
        FLUSH(pi_comm_unit)
      endif

      end subroutine commforces_ctr

      end module commstuffpi
