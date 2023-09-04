c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
      module commstuffpi
      implicit none

      type real3
        real*8 x,y,z
      end type real3
      type t_elt
        real*8 :: x,y,z,vx,vy,vz,ax,ay,az
      end type
      integer:: max_atoms_send=0
      integer:: max_atoms_recv=0
      integer, allocatable :: loc_save(:)
      save

      contains

      subroutine full_gradient_prepare_pi(polymer,polymer_ctr
     &             ,istep,respa)
      use beads
      use domdec
      use atoms
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
      use atoms
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
        use atoms ,only : x,y,z,xold,yold,zold,n
        use beads
        use bound ,only : use_bounds
        use cell  ,only : xcell2, ycell2, zcell2
        use domdec,only : nloc,nproc,glob,comm_dir,rank,rank_bis,ndir,
     &                    xbegproc,ybegproc,zbegproc,
     &                    xendproc,yendproc,zendproc,
     &                    repart, loc, domlen, COMM_TINKER,
     &                    nneig_send, nneig_recep,
     &                    pneig_recep, pneig_send, 
     &                    nproc_polymer, COMM_POLYMER, rank_polymer
        use freeze,only : use_rattle
        use moldyn,only : a,v
        use potent,only : use_pmecore
        use mpi
        use qtb
        implicit none

        type(POLYMER_COMM_TYPE), intent(inout) :: polymer
        type(POLYMER_COMM_TYPE), intent(inout),optional :: polymer_ctr
        integer,allocatable :: iglob_send(:,:) ! Global index to send to each proc
        integer,allocatable :: iglob_recv(:,:) ! Global index to recieve from procs
        
        real*8 xr, yr, zr ! adjusted position
        type(t_elt),pointer :: d
        type(real3),pointer :: b
        type(real3),pointer :: d_ctr
        type(real3),pointer :: old_send(:,:),old_recv(:,:)

        integer  n_data_send(0:nneig_send),   n_data_recv(nneig_recep)
        integer req_iglob_send(nneig_send),req_iglob_recv(nneig_recep)
        integer  req_data_send(nneig_send), req_data_recv(nneig_recep)
        integer  mpi_status1(MPI_STATUS_SIZE)
        type(t_elt), allocatable,target :: data_send(:,:)
     &                                       , data_recv(:,:)
        type(real3), allocatable,target :: data_send_ctr(:,:)
     &                                       , data_recv_ctr(:,:)

        integer ibeadbeg,ibeadend, nbeadsproc, nbeadssend 
        integer ibeadbeg_ctr,ibeadend_ctr, nbeadsproc_ctr
        integer nprocloc, commloc, rankloc ! local (pme) mpi info
        integer i,k,j, ierr, iglob, iloc, iproc, ineighbor
        integer nloc_save, n_data_send_capture,max_data_recv
        integer nloc_capture
        integer s_bufi
        integer :: max_data_recv_save=0
        character(len=40) ::fmt
        real*8 :: sqrtnu
        real*8 ixbeg,ixend,iybeg,iyend,izbeg,izend,xtemp
        logical :: send_ctr
        real*8, parameter:: eps1= 1.0d-10, eps2= 1.0d-8

        send_ctr = present(polymer_ctr)

        if (use_pmecore) then
           nprocloc = ndir
           commloc  = comm_dir
           rankloc  = rank_bis
        else
           nprocloc = nproc
           commloc  = COMM_TINKER
           rankloc  = rank
        end if
        if (nproc.eq.1) return
        if (use_pmecore.and.ndir.eq.1) then
           goto 20
        end if

        s_bufi = nloc*(nneig_send+1)
        !Associate iglob_send pointer to buffer
        allocate(iglob_send(1:nloc,0:nneig_send))

        if(.not. allocated(loc_save)) then
          allocate(loc_save(n))
        endif

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

        n_data_send(:) = 0
        ! Get process domain of each atoms
        do i = 1, nloc
          iglob = glob(i)
          xr = polymer%eigpos(1,iglob,1)
          yr = polymer%eigpos(2,iglob,1)
          zr = polymer%eigpos(3,iglob,1)
          if (use_bounds) call image(xr,yr,zr)
          ! Check box limits (SP Issue)
          if (abs(xr-xcell2).lt.eps1) xr = xr-eps2
          if (abs(yr-ycell2).lt.eps1) yr = yr-eps2
          if (abs(zr-zcell2).lt.eps1) zr = zr-eps2
          if  ( (xr>=ixbeg).and.(xr<ixend)
     &    .and. (yr>=iybeg).and.(yr<iyend)
     &    .and. (zr>=izbeg).and.(zr<izend) ) then  ! (Rank domain)
              n_data_send(0) = n_data_send(0) + 1
              iglob_send(n_data_send(0),0) = glob(i)
          else  ! (not in rank domain)
            do ineighbor = 1, nneig_send   ! Look among neighbour
              iproc = pneig_send(ineighbor)+1
              if( (iproc /= rank+1)
     &  .and. (xr>=xbegproc(iproc)) .and. (xr<xendproc(iproc))
     &  .and. (yr>=ybegproc(iproc)) .and. (yr<yendproc(iproc))
     &  .and. (zr>=zbegproc(iproc)) .and. (zr<zendproc(iproc)) ) then
                n_data_send(ineighbor) = n_data_send(ineighbor) + 1
                iglob_send(n_data_send(ineighbor),ineighbor) = glob(i)
                exit ! if found
              endif
            enddo
          endif ! (Find atom domain)
        enddo

        ! Allocate and fill data to send (Optimzed reallocation)
        max_atoms_send = maxval(n_data_send(1:nneig_send))
        allocate(data_send((nbeadsproc+1)*max_atoms_send,nneig_send))

        
        do ineighbor = 1, nneig_send  ! Gather Data for each neighbor domain
          do i = 1, n_data_send(ineighbor)
            iglob = iglob_send(i,ineighbor)
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
        do ineighbor = 1, nneig_send 
          call MPI_Isend(iglob_send(1,ineighbor)
     &            ,n_data_send(ineighbor)
     &            ,MPI_INT,   pneig_send(ineighbor), 0, COMM_TINKER,
     &            req_iglob_send(ineighbor), ierr)
        end do

        ! Probe for sizes
        n_data_recv(:) = 0
        do ineighbor = 1, nneig_recep
          call MPI_Probe(pneig_recep(ineighbor), 0, COMM_TINKER,
     &                   mpi_status1, ierr)
          call MPI_Get_count( mpi_status1, MPI_INT,
     &                        n_data_recv(ineighbor), ierr )
        end do

        ! Allocate and recv data depending on message size
        max_atoms_recv = maxval(n_data_recv)
        allocate(iglob_recv(1:max_atoms_recv,1:nneig_recep))
        allocate(data_recv((nbeadsproc+1)*max_atoms_recv,nneig_recep))

        req_iglob_recv(:) = MPI_REQUEST_NULL
        req_data_recv(:)  = MPI_REQUEST_NULL
!Wait for iglob_recv and data_recv
        do ineighbor = 1, nneig_recep
          call MPI_Irecv(iglob_recv(1,ineighbor),n_data_recv(ineighbor)
     &            ,MPI_INT, pneig_recep(ineighbor), 0, COMM_TINKER,
     &            req_iglob_recv(ineighbor), ierr)
          call MPI_Irecv(data_recv(1,ineighbor)
     &            ,9*n_data_recv(ineighbor)*(nbeadsproc+1)
     &            ,MPI_REAL8, pneig_recep(ineighbor), 1, COMM_TINKER,
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
     &            ,MPI_REAL8, pneig_send(ineighbor), 1, COMM_TINKER,
     &            req_data_send(ineighbor), ierr)
        end do

        nloc_save = nloc
        ! Rebuild glob(:) with old local atoms
        nloc = n_data_send(0)
        do i=1,nloc
          glob(i) = iglob_send(i, 0)
        enddo

        loc_save(:)=loc(:)
        loc(:)=0
        repart(:)=-1

        ! Wait all communications
        call MPI_Waitall(nneig_recep,req_data_recv ,MPI_STATUSES_IGNORE,
     &                   ierr)
        call MPI_Waitall(nneig_send ,req_data_send ,MPI_STATUSES_IGNORE,
     &                   ierr)

        sqrtnu = sqrt(real(nbeads,8))

        ! Add new atoms to glob(:) and add data
        do ineighbor = 1, nneig_recep
          do i = 1, n_data_recv(ineighbor)
            iglob = iglob_recv(i,ineighbor)
            nloc = nloc+1
            glob(nloc) = iglob
            d => data_recv((i-1)*(nbeadsproc+1)+1,ineighbor)
            polymer%eigpos(1,iglob,1) = d%x
            polymer%eigpos(2,iglob,1) = d%y
            polymer%eigpos(3,iglob,1) = d%z
            polymer%eigvel(1,iglob,1) = d%vx
            polymer%eigvel(2,iglob,1) = d%vy
            polymer%eigvel(3,iglob,1) = d%vz
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

          req_data_recv(:)  = MPI_REQUEST_NULL
          req_data_send(:)  = MPI_REQUEST_NULL
          do ineighbor = 1, nneig_recep
            call MPI_Irecv(data_recv_ctr(1,ineighbor)
     &              ,3*n_data_recv(ineighbor)*nbeadsproc_ctr
     &              ,MPI_REAL8, pneig_recep(ineighbor), 1, COMM_TINKER,
     &              req_data_recv(ineighbor), ierr)
          end do


          do ineighbor = 1, nneig_send  ! Gather Data for each neighbor domain
            do i = 1, n_data_send(ineighbor)
              iglob = iglob_send(i,ineighbor)
              do k=1,nbeadsproc_ctr
                d_ctr => data_send_ctr((i-1)*nbeadsproc_ctr+k
     &                     , ineighbor)
                d_ctr%x  = polymer_ctr%pos(1,iglob,k+ibeadbeg_ctr-1)
                d_ctr%y  = polymer_ctr%pos(2,iglob,k+ibeadbeg_ctr-1)
                d_ctr%z  = polymer_ctr%pos(3,iglob,k+ibeadbeg_ctr-1)
              enddo
            end do
          end do

          do ineighbor = 1, nneig_send 
            call MPI_Isend(data_send_ctr(1,ineighbor)
     &              ,3*n_data_send(ineighbor)*nbeadsproc_ctr
     &              ,MPI_REAL8, pneig_send(ineighbor), 1, COMM_TINKER,
     &              req_data_send(ineighbor), ierr)
          end do

          call MPI_Waitall(nneig_recep,req_data_recv 
     &                     ,MPI_STATUSES_IGNORE,ierr)
          call MPI_Waitall(nneig_send ,req_data_send 
     &                  ,MPI_STATUSES_IGNORE, ierr)

          do ineighbor = 1, nneig_recep
            do i = 1, n_data_recv(ineighbor)
              iglob = iglob_recv(i,ineighbor)
              do k=1,nbeadsproc_ctr
                d_ctr => data_recv_ctr((i-1)*nbeadsproc_ctr+k
     &             ,ineighbor)
                polymer_ctr%pos(1,iglob,k+ibeadbeg_ctr-1)   = d_ctr%x
                polymer_ctr%pos(2,iglob,k+ibeadbeg_ctr-1)   = d_ctr%y
                polymer_ctr%pos(3,iglob,k+ibeadbeg_ctr-1)   = d_ctr%z
              enddo
            end do
          end do

        endif

        if(piqtb) call reassignqtb(nloc_save,max_atoms_recv
     &      ,nneig_recep, n_data_recv, iglob_recv, pneig_recep 
     &      ,nneig_send , n_data_send, iglob_send, pneig_send 
     &      ,max_atoms_recv, max_atoms_send)


20      call orderbuffer(.false.)
        call update_nlocpi(nloc)

      end subroutine reassignpi


      subroutine commpospi(polymer,send_centroid,fast)
      use atoms
      use beads
      use domdec
      use freeze
      use mpi
      use potent
      use beads
      implicit none
      type(POLYMER_COMM_TYPE), intent(inout) :: polymer
      logical, intent(in) :: send_centroid
      logical, intent(in) :: fast
      integer i,k,j,tag,ierr,iproc
      integer iloc,iglob,idomlen,ibufbeg
      real*8, allocatable:: buffer(:,:,:),buffers(:,:,:)
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
c     communicate positions
c
      !Fetch memory from pool
      allocate(reqsend(nproc))
      allocate(reqrec(nproc))
      allocate(buffer (1:3,1:nbeadssend,1:nbloc))
      allocate(buffers(1:3,1:nbeadssend,1:nloc))
c
c     MPI : begin reception in buffer
c
      do i = 1, ntemp_recep
        tag = nproc*rank + ptemp_recep(i) + 1
        call MPI_IRECV(buffer(1,1,bufbeg(ptemp_recep(i)+1)),
     $   3*domlen(ptemp_recep(i)+1)*nbeadssend,
     $   MPI_REAL8,ptemp_recep(i),tag,COMM_TINKER,reqrec(i),ierr)
      end do
c
c     MPI : move in buffer
c
      if(send_centroid) then
        do i = 1, nloc
          iglob = glob(i)
          buffers(:,1,i) = polymer%eigpos(:,iglob,1)
        end do
      endif
      do i = 1, nloc
        iglob = glob(i)
        do k=1,nbeadsproc
          buffers(:,k+cshift,i) = polymer%pos(:,iglob,k+ibeadbeg-1)
        enddo
      end do 
c
c     MPI : begin sending
c
      do i = 1, ntemp_send 
        tag = nproc*ptemp_send(i) + rank + 1
        call MPI_ISEND(buffers,3*nloc*nbeadssend 
     $       ,MPI_REAL8,
     $       ptemp_send(i),tag,COMM_TINKER,
     $       reqsend(i),ierr)
      end do
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
         if(send_centroid) then
          do i = 1, idomlen
            iloc   = ibufbeg+i-1
            iglob  = glob(iloc)
            polymer%eigpos(:,iglob,1)=buffer(:,1,iloc)
          end do
         endif
         do i = 1, idomlen
          iloc   = ibufbeg+i-1
          iglob  = glob(iloc)
          do k = 1, nbeadsproc
            polymer%pos(:,iglob,k+ibeadbeg-1) 
     &           = buffer(:,k+cshift,iloc)
          enddo
         end do
      end do
c
      end subroutine commpospi

      subroutine commposrecpi(polymer,send_centroid)
      use atoms
      use domdec
      use mpi
      use potent   ,only: use_pmecore
      use beads
      implicit none
      type(POLYMER_COMM_TYPE), intent(inout) :: polymer
      logical, intent(in) :: send_centroid
      integer i,iproc,iprec
      integer iglob,iloc,j,k
      integer tag,ierr,ibufbeg,idomlen
      integer rankloc,commloc,nprocloc
      integer status(MPI_STATUS_SIZE)
      real*8,allocatable :: buffer(:,:,:),buffers(:,:,:)
      integer nbeadsproc,ibeadbeg,ibeadend,nbeadssend 
      integer cshift    
      integer, allocatable :: reqrec(:),reqsend(:) 

      if (nproc.eq.1) return

      if(send_centroid) then
        cshift=1
      else
        cshift=0
      endif
c
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

      allocate (reqsend(nproc))
      allocate (reqrec(nproc))
c
c     allocate some arrays
c
      allocate(buffer (1:3,1:nbeadssend,1:nblocrecdir))
      allocate(buffers(1:3,1:nbeadssend,1:nloc))   
c
c     first do dir rec communications
c
c
c     MPI : begin reception in buffer
c
      do i = 1, nrecdir_recep2
        if (precdir_recep2(i).ne.rank) then
          tag = nproc*rank + precdir_recep2(i) + 1
          call MPI_IRECV(buffer(1,1,bufbeg(precdir_recep2(i)+1)),
     $         3*domlen(precdir_recep2(i)+1)*nbeadssend,MPI_REAL8,
     $         precdir_recep2(i),tag,COMM_TINKER,reqrec(i),ierr)
        end if
      end do
c
c     MPI : move in buffer
c
      if(send_centroid) then
        do i = 1, nloc
          iglob = glob(i)
          buffers(:,1,i) = polymer%eigpos(:,iglob,1)
        enddo
      endif
      do i = 1, nloc
        iglob = glob(i)
        do k=1,nbeadsproc
        buffers(:,k+cshift,i) = polymer%pos(:,iglob,k+ibeadbeg-1)
        enddo
      end do
c
c     MPI : begin sending
c
      do i = 1, nrecdir_send2
        if (precdir_send2(i).ne.rank) then
          tag = nproc*precdir_send2(i) + rank + 1
          call MPI_ISEND(buffers,3*nloc*nbeadssend 
     &       ,MPI_REAL8,precdir_send2(i),
     &         tag,COMM_TINKER,reqsend(i),ierr)
        end if
      end do
c
c     MPI : move in global arrays
c
      do iproc = 1, nrecdir_recep2
        if (precdir_recep2(iproc)==rank) cycle
        call MPI_WAIT(reqrec(iproc),status,ierr)
        ibufbeg = bufbeg(precdir_recep2(iproc)+1)
        idomlen = domlen(precdir_recep2(iproc)+1)
        if(send_centroid) then
          do i = 1, idomlen
            iloc     = ibufbeg+i-1
            iglob    = glob(iloc)
            polymer%eigpos(:,iglob,1) = buffer(:,1,iloc)
          enddo
        endif
        do i = 1, idomlen
          iloc     = ibufbeg+i-1
          iglob    = glob(iloc)
          do k=1,nbeadsproc
            polymer%pos(:,iglob,k+ibeadbeg-1) 
     &          = buffer(:,k+cshift,iloc) 
          enddo
        enddo
      end do
c
      do i = 1, nrecdir_send2
         if (precdir_send2(i).ne.rank) then
            call MPI_WAIT(reqsend(i),status,ierr)
         end if
      end do
c
c     then do rec rec communications
c
      deallocate(buffer,buffers)
      allocate(buffer (1:3,1:nbeadssend,1:nlocrec2))
      allocate(buffers(1:3,1:nbeadssend,1:nlocrec))

c
c     MPI : begin reception in buffer
c
      do i = 1, nrec_recep
        tag = nprocloc*rankloc + prec_recep(i) + 1
        call MPI_IRECV(buffer(1,1,bufbegrec(prec_recep(i)+1))
     $      ,3*domlenrec(prec_recep(i)+1)*nbeadssend,MPI_REAL8
     $      ,prec_recep(i),tag,commloc,reqrec(i),ierr)
      end do
c
c     MPI : move in buffer
c
      if(send_centroid) then
        do i = 1, nlocrec
          iglob = globrec(i)
          buffers(:,1,i) = polymer%eigpos(:,iglob,1)
        end do
      endif
      do i = 1, nlocrec
        iglob = globrec(i)
        do k=1,nbeadsproc
          buffers(:,k+cshift,i) = polymer%pos(:,iglob,k+ibeadbeg-1)
        enddo
      end do 
c
c     MPI : begin sending
c
      do i = 1, nrec_send 
         tag = nprocloc*prec_send(i) + rankloc+1
         call MPI_ISEND(buffers,3*nlocrec*nbeadssend 
     $      ,MPI_REAL8,prec_send(i)
     $       ,tag,commloc,reqsend(i),ierr)
      end do
c
c     MPI : move in global arrays
c
      do iproc = 1, nrec_recep
        call MPI_WAIT(reqrec(iproc),status,ierr)
        ibufbeg = bufbegrec(prec_recep(iproc)+1)
        idomlen = domlenrec(prec_recep(iproc)+1)
        if(send_centroid) then
          do i = 1, idomlen
            iloc     = ibufbeg+i-1
            iglob    = globrec(iloc)
            polymer%eigpos(:,iglob,1) = buffer(:,1,iloc)
          enddo
        endif
        do i = 1, idomlen
          iloc     = ibufbeg+i-1
          iglob    = globrec(iloc)
          do k=1,nbeadsproc
            polymer%pos(:,iglob,k+ibeadbeg-1) 
     &           = buffer(:,k+cshift,iloc) 
          enddo
        enddo
      end do
      ! Wait for sending requests
      do i = 1, nrec_send 
        call MPI_WAIT(reqsend(i),status,ierr)
      end do

      end subroutine commposrecpi


      subroutine commforcesblocpi(derivs,nbeadsproc,fast)
      use beads
      use deriv
      use domdec
      use mpi
      use potent
      use sizes
      implicit none
      integer, intent(in) :: nbeadsproc
      real*8, intent(inout) :: derivs(3,nbloc,*)
      logical, intent(in) :: fast
      integer i,j,k,tag,ierr,sz,sz2,ii
      integer reqsend(nproc),reqrec(nproc)
      real*8, allocatable :: buffer(:,:,:,:),buffers(:,:,:)
!      real*8 derivs(3,nbloc,nbeads)
      integer,dimension(nproc):: ptemp_recep,ptemp_send
      integer status(MPI_STATUS_SIZE)
      integer ntemp_recep,ntemp_send 

c
      if (nproc.eq.1.or.((use_pmecore).and.(rank.gt.ndir-1))) return

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
      sz = 3*max(1,nloc)*ntemp_send*nbeadsproc
      sz2 = 3*nbloc*nbeadsproc
      allocate(buffer(1:3,1:nbeadsproc,1:max(1,nloc),1:ntemp_send))
      allocate(buffers(1:3,1:nbeadsproc,1:nbloc))

      do k=1,nbeadsproc; do i=1,nbloc; do j=1,3
        buffers(j,k,i) = derivs(j,i,k)
      end do ; enddo ; enddo
c
c     communicate forces
c
c     MPI : begin reception in buffer
c
      do i = 1, ntemp_send 
         tag = nproc*rank + ptemp_send(i) + 1
         call MPI_IRECV(buffer(1,1,1,i),3*nloc*nbeadsproc,
     $        MPI_REAL8,ptemp_send(i),tag,
     $        COMM_TINKER,reqrec(i),ierr)
      end do
c
c     MPI : begin sending
c

      do i = 1, ntemp_recep
         tag = nproc*ptemp_recep(i) + rank + 1
         call MPI_ISEND(buffers(1,1,bufbeg(ptemp_recep(i)+1)),
     $        3*domlen(ptemp_recep(i)+1)*nbeadsproc
     $        ,MPI_REAL8,ptemp_recep(i),tag,
     $        COMM_TINKER,reqsend(i),ierr)
      end do
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
      do ii = 1, ntemp_send 
        do k=1,nbeadsproc; do i = 1, nloc; do j = 1, 3
          derivs(j,i,k) = derivs(j,i,k) + buffer(j,k,i,ii)
        enddo ; enddo; enddo
      end do
c
      end subroutine commforcesblocpi


      subroutine comm_for_gradient(polymer,send_centroid,send_vel)
      use beads
      use mpi
      use beads
      use domdec
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
      real*8 :: sqrtnu
      integer :: nbeadslocmax, nlocproc
      integer :: ibeadbeg,ibeadend,nu
      integer :: ilocend,ilocbeg
      integer :: cshift,nelt
      logical :: send_vel_
      real*8, allocatable :: bufferpi_s(:,:,:)
      real*8, allocatable :: bufferpi_r(:,:,:)

      if(nproc_polymer==1) return
      reqs(:)=MPI_REQUEST_NULL

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
      
        !PUT DATA INTO BUFFER ARRAYS

      ilocbeg=ilocpi_beg(rank_polymer+1)
      ilocend=ilocpi_end(rank_polymer+1)

      DO k=1,nu; DO iloc=0,ilocend-ilocbeg !; DO j=1,3
        iproc=polymer%iproc(k)
        if(rank_polymer==iproc-1)  CYCLE
        ibeadbeg=polymer%ibead_beg(iproc)
        i=glob(ilocbeg+iloc)
        if(k-ibeadbeg==0 .and. send_centroid) then
          kk=iloc +1
          DO j=1,3
          bufferpi_s(j,kk,iproc) = polymer%eigpos(j,i,1)
          ENDDO
          if(nelt==2) then
            DO j=1,3
              bufferpi_s(3+j,kk,iproc) = polymer%eigvel(j,i,1)
            ENDDO
          endif
        endif
        kk=(k-ibeadbeg + cshift)*nlocpi+ iloc +1
        DO j=1,3
        bufferpi_s(j,kk,iproc)=polymer%pos(j,i,k)
        ENDDO
        if(nelt==2) then
          DO j=1,3
            bufferpi_s(3+j,kk,iproc)=polymer%vel(j,i,k)
          ENDDO
        endif
      ENDDO; ENDDO !; ENDDO

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
     &      ,MPI_REAL8
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
     &   ,MPI_REAL8,iproc-1,mptag_data, COMM_POLYMER
     &   ,reqs(ii),ierr)
      ENDDO

      !WAIT FOR COMMUNICATIONS
      if(pi_comm_report) then
        write(pi_comm_unit,*) "  Wait for data"
        FLUSH(pi_comm_unit)
      endif
      call MPI_WAITALL(ii,reqs,MPI_STATUSES_IGNORE,ierr)
      

      !PUT DATA INTO POLYMER ARRAYS

        if(send_centroid) then
          DO iloc=1,nloc
            iproc=min(nproc_polymer,int((iloc-1)/nlocproc)+1)
            if(rank_polymer==iproc-1) CYCLE
            ilocbeg=ilocpi_beg(iproc)
            ilocend=ilocpi_end(iproc)
            i=glob(iloc)
            kk=iloc-ilocbeg + 1
            DO j=1,3
            polymer%eigpos(j,i,1)=bufferpi_r(j,kk,iproc)
            ENDDO
            if(nelt==2) then
              DO j=1,3
                polymer%eigvel(j,i,1)=bufferpi_r(3+j,kk,iproc)
              ENDDO
            endif
          enddo
        endif
        DO k=1,ibeadend-ibeadbeg+1; DO iloc=1,nloc !; DO j=1,3
          iproc=min(nproc_polymer,int((iloc-1)/nlocproc)+1)
          if(rank_polymer==iproc-1) CYCLE
          ilocbeg=ilocpi_beg(iproc)
          ilocend=ilocpi_end(iproc)
          i=glob(iloc)
          
          ibead=k+ibeadbeg-1
          kk=(k+cshift-1)*(ilocend-ilocbeg+1)+ iloc-ilocbeg +1
          DO j=1,3
          polymer%pos(j,i,ibead)=bufferpi_r(j,kk,iproc)
          ENDDO
          if(nelt==2) then
            DO j=1,3
              polymer%vel(j,i,ibead)=bufferpi_r(3+j,kk,iproc)
            ENDDO
          endif
        ENDDO; ENDDO !; ENDDO

      deallocate(bufferpi_s,bufferpi_r)

      if(pi_comm_report) then
        write(pi_comm_unit,*) "  exiting comm_for_gradient"
        FLUSH(pi_comm_unit)
      endif

      end subroutine comm_for_gradient

      subroutine comm_for_normal_modes(polymer
     &   ,array1,array2,array3,array4)
      use beads
      use mpi
      use beads
      use domdec
      implicit none
      type(POLYMER_COMM_TYPE),intent(inout) :: polymer
      real*8, intent(inout) :: array1(:,:,:)
      real*8, intent(inout), optional :: array2(:,:,:)
      real*8, intent(inout), optional :: array3(:,:,:)
      real*8, intent(inout), optional :: array4(:,:,:)
      logical :: send_pos_, send_vel_,slow_
      logical :: correctorder
      integer :: ibead,k,i,j,iglob,iproc,maxsend,nsend,maxloc
      integer :: maxrecv,nproc_recv
      integer :: isend,ii,ierr,ns,nr,iloc,kk
      integer, parameter :: mptag_size=0,mptag_data=1
      integer :: reqs(2*nproc_polymer)
      LOGICAL :: flag
      real*8 :: sqrtnu
      integer :: nbeadslocmax,nlocproc
      integer :: ibeadbeg,ibeadend,nu
      integer :: ilocend,ilocbeg,nelt,ielt
      real*8, allocatable :: bufferpi_s(:,:,:)
      real*8, allocatable :: bufferpi_r(:,:,:)

      if(nproc_polymer==1) return

      reqs(:)=MPI_REQUEST_NULL

      nu=polymer%nbeads
      sqrtnu=sqrt(real(nu,8))

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

        !PUT DATA INTO BUFFER ARRAYS
      ibeadbeg = polymer%ibead_beg(rank_polymer+1)
      ibeadend = polymer%ibead_end(rank_polymer+1)
      DO k=0,ibeadend-ibeadbeg; DO iloc=1,nloc !; DO j=1,3
        iproc=min(nproc_polymer,int((iloc-1)/nlocproc)+1)
        if(rank_polymer==iproc-1) cycle
        ilocbeg=ilocpi_beg(iproc)
        ilocend=ilocpi_end(iproc)
        i=glob(iloc)
        ibead=k+ibeadbeg
        kk=k*(ilocend-ilocbeg+1)+ iloc-ilocbeg + 1
        do j=1,3
          bufferpi_s(  j,kk,iproc)=array1(j,i,ibead)
        enddo
        if(nelt==1) cycle
        do j=1,3
          bufferpi_s(3+j,kk,iproc)=array2(j,i,ibead)
        enddo
        if(nelt==2) cycle
        do j=1,3
          bufferpi_s(6+j,kk,iproc)=array3(j,i,ibead)
        enddo
        if(nelt==3) cycle
        do j=1,3
          bufferpi_s(9+j,kk,iproc)=array4(j,i,ibead)
        enddo
      ENDDO; ENDDO !; ENDDO

      ii=0
      DO iproc=1,nproc_polymer
        if(rank_polymer==iproc-1) CYCLE
        ii=ii+1
        ilocbeg=ilocpi_beg(iproc)
        ilocend=ilocpi_end(iproc)
        if(pi_comm_report) then
          write(pi_comm_unit,*) rank_polymer,"  send"
     &       ,nelt*3*(ilocend-ilocbeg+1)*(ibeadend-ibeadbeg+1)
     &         ,"ats to",iproc-1
        endif
        call MPI_ISEND(bufferpi_s(1,1,iproc)
     &     ,nelt*3*(ilocend-ilocbeg+1)*(ibeadend-ibeadbeg+1)
     &      ,MPI_REAL8
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
          write(pi_comm_unit,*) rank_polymer,"  receive"
     &       ,nelt*3*(ilocend-ilocbeg+1)*(ibeadend-ibeadbeg+1)
     &         ,"ats from",iproc-1
        endif
        call MPI_IRECV(bufferpi_r(1,1,iproc)
     &   ,nelt*3*(ilocend-ilocbeg+1)*(ibeadend-ibeadbeg+1)
     &   ,MPI_REAL8,iproc-1,mptag_data, COMM_POLYMER
     &   ,reqs(ii),ierr)
      ENDDO

      !WAIT FOR COMMUNICATIONS
      call MPI_WAITALL(ii,reqs,MPI_STATUSES_IGNORE,ierr)

      !PUT DATA INTO POLYMER ARRAYS

        DO k=1,nu; 
          iproc=polymer%iproc(k)
          if(rank_polymer==iproc-1) cycle

          DO iloc=0,ilocend-ilocbeg !; DO j=1,3
            ibeadbeg=polymer%ibead_beg(iproc)
            i=glob(ilocbeg+iloc)
            kk=(k-ibeadbeg)*nlocpi + iloc + 1
            do j=1,3
              array1(j,i,k)=bufferpi_r(j,kk,iproc)
            enddo
            if(nelt==1) cycle
            do j=1,3
              array2(j,i,k)=bufferpi_r(3+j,kk,iproc)
            enddo
            if(nelt==2) cycle
            do j=1,3
              array3(j,i,k)=bufferpi_r(6+j,kk,iproc)
            enddo
            if(nelt==3) cycle
            do j=1,3
              array4(j,i,k)=bufferpi_r(9+j,kk,iproc)
            enddo
          ENDDO
        ENDDO  !; ENDDO

      deallocate(bufferpi_s,bufferpi_r)

      if(pi_comm_report) then
        write(pi_comm_unit,*) "  exiting comm_for_normal_modes"
        FLUSH(pi_comm_unit)
      endif


      end subroutine comm_for_normal_modes

      subroutine save_deriv_rec(temprec)
      use domdec
      use deriv
      implicit none
      real*8, intent(inout) :: temprec(3,nlocrec2)
      integer :: i,j

      do i = 1, nlocrec2
        do j = 1, 3
            temprec(j,i)= demrec(j,i) 
     &                  + deprec(j,i) 
     &                  + decrec(j,i)
     &                  + dedsprec(j,i)
            demrec(j,i)=0.d0
            deprec(j,i)=0.d0
            decrec(j,i)=0.d0
            dedsprec(j,i)=0.d0
        end do
      end do


      end subroutine save_deriv_rec

      subroutine commforcesrecpi(derivs,temprec,nbeadsproc)
      use deriv
      use domdec
      use mpi
      use potent
      use sizes
      implicit none
      integer, intent(in) :: nbeadsproc
      real*8, intent(inout) :: derivs(3,nbloc,*)
      real*8, intent(inout) :: temprec(3,nlocrec2,*)
      integer i,j,k,tag,ierr,iproc,iglob,iloc,ilocrec,ibufbeg
      integer jloc,jglob,jlocrec,ii
      integer sz1,sz2,sz3
      integer rankloc,commloc,nprocloc
      integer reqsend(nproc),reqrec(nproc)
      integer status(MPI_STATUS_SIZE)
      real*8,allocatable:: buffer(:,:,:,:)
      real*8,allocatable:: buffers(:,:,:),temprecbuf(:,:,:)

      !case sequential
      if (nproc.eq.1) then

         do k= 1,nbeadsproc; do i = 1, nloc
            derivs(:,i,k) = derivs(:,i,k) + temprec(:,i,k)
         end do; end do

         return
      end if
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

      allocate(buffer(1:3,1:nbeadsproc,1:max(1,nlocrec),1:nrec_send))
      allocate(buffers(1:3,1:nbeadsproc,1:max(1,nblocrecdir)))        
      allocate(temprecbuf(1:3,1:nbeadsproc,1:max(1,nlocrec2)))       
c
c     first do rec rec communications
c
c
c     MPI : begin reception in buffer
c
      do i = 1, nrec_send 
        tag = nprocloc*rankloc + prec_send(i) + 1
        call MPI_IRECV(buffer(1,1,1,i),3*nlocrec*nbeadsproc,MPI_REAL8
     $      ,prec_send(i),tag,commloc,reqrec(i),ierr)
      end do
c
c     MPI : move in buffer
c

      do k= 1,nbeadsproc; do i = 1, nlocrec2
        temprecbuf(:,k,i) = temprec(:,i,k)
      end do; end do
c
      do i = 1, nrec_recep
        tag = nprocloc*prec_recep(i) + rankloc + 1
        call MPI_ISEND(temprecbuf(1,1,bufbegrec(prec_recep(i)+1))
     $      ,3*domlenrec(prec_recep(i)+1)*nbeadsproc
     $      ,MPI_REAL8
     $      ,prec_recep(i),tag,commloc,reqsend(i),ierr)
      end do
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
      do ii = 1,nrec_send 
        do k=1,nbeadsproc; do i = 1,nlocrec
          temprec(:,i,k) = temprec(:,i,k) + buffer(:,k,i,ii)
        end do; enddo
      end do

      sz1 = 3*max(1,nloc)*nrecdir_send1*nbeadsproc
      deallocate(buffer)
      allocate(buffer(1:3,1:nbeadsproc,1:max(1,nloc),1:nrecdir_send1))
c
c     then do rec dir communications
c
c     MPI : begin reception in buffer
c
      do i = 1, nrecdir_send1
         if (precdir_send1(i).ne.rank) then
           tag = nproc*rank + precdir_send1(i) + 1
           call MPI_IRECV(buffer(1,1,1,i),3*nloc*nbeadsproc,MPI_REAL8
     $         ,precdir_send1(i),tag,COMM_TINKER,reqrec(i),ierr)
         end if
      end do
c
c     Move in array
c
      do ii = 1, nrecdir_recep1
        if (precdir_recep1(ii)==rank) cycle
        ibufbeg =bufbeg(precdir_recep1(ii)+1)-1 
        do k=1,nbeadsproc; do i = 1, domlen(precdir_recep1(ii)+1)
          iloc = ibufbeg+i
          iglob = glob(iloc)
          if (repartrec(iglob).eq.rankloc) then
            ilocrec = locrec(iglob)
            buffers(:,k,iloc) = temprec(:,ilocrec,k)
          else
            buffers(:,k,iloc) = 0
          end if
        end do; enddo
      end do
c
      do i = 1, nrecdir_recep1
        if (precdir_recep1(i).ne.rank) then
          tag = nproc*precdir_recep1(i) + rank + 1
          call MPI_ISEND(buffers(1,1,bufbeg(precdir_recep1(i)+1))
     $        ,3*domlen(precdir_recep1(i)+1)*nbeadsproc,MPI_REAL8
     $        ,precdir_recep1(i),tag,COMM_TINKER,reqsend(i),ierr)
        end if
      end do
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
          do k=1,nbeadsproc; do i = 1,nloc
            derivs(:,i,k) = derivs(:,i,k) + buffer(:,k,i,ii)
          end do; enddo
        else
          do k=1,nbeadsproc; do i = 1, nlocrec
            iglob = globrec(i)
            iloc  = loc(iglob)
            if (repart(iglob).eq.rank) then
              derivs(:,iloc,k) = derivs(:,iloc,k) + temprec(:,i,k)
            end if
          end do; enddo
        end if
      end do

      do i = 1, nrecdir_recep1
        if (precdir_recep1(i).ne.rank) then
        call MPI_WAIT(reqsend(i),status,ierr)
        end if
      end do

      end subroutine commforcesrecpi


      subroutine commforces_ctr(polymer_ctr,derivs_ctr)
      use beads
      use mpi
      use beads
      use domdec
      implicit none
      type(POLYMER_COMM_TYPE),intent(inout) :: polymer_ctr
      real*8, intent(inout) :: derivs_ctr(3,nbloc,*)
      integer :: ibead,k,i,j,iglob,iproc,maxsend,nsend,maxloc
      integer :: maxrecv,nproc_recv
      integer :: isend,ii,ierr,ns,nr,iloc,kk
      integer, parameter :: mptag_size=0,mptag_data=1
      integer :: reqs(2*nproc_polymer)
      LOGICAL :: flag
      real*8 :: sqrtnu
      integer :: nbeadslocmax,nlocproc,nbeadsproc
      integer :: ibeadbeg,ibeadend,nu
      integer :: ilocend,ilocbeg
      real*8, allocatable :: bufferpi_s(:,:,:)
      real*8, allocatable :: bufferpi_r(:,:,:)

      nu=polymer_ctr%nbeads

      ilocbeg=ilocpi_beg(rank_polymer+1)
      ilocend=ilocpi_end(rank_polymer+1)
      ibeadbeg = polymer_ctr%ibead_beg(rank_polymer+1)
      ibeadend = polymer_ctr%ibead_end(rank_polymer+1)
      nbeadsproc = ibeadend-ibeadbeg+1

      DO k=1,nbeadsproc; DO iloc=1,nloc
        i=glob(iloc)
        polymer_ctr%forces_slow(:,i,k+ibeadbeg-1)
     &     =-derivs_ctr(:,iloc,k)
      ENDDO; ENDDO
        
      if(nproc_polymer==1)   return

      reqs(:)=MPI_REQUEST_NULL

      if(pi_comm_report) then
        write(pi_comm_unit,*) "  entering commforces_ctr"
        FLUSH(pi_comm_unit)
      endif
     
      maxloc=nlocpis(nproc_polymer)
      nbeadslocmax=polymer_ctr%nbeadsloc(nproc_polymer)
      nlocproc=int(nloc/nproc_polymer)

      allocate(bufferpi_s(3,maxloc*nbeadslocmax,nproc_polymer))
      allocate(bufferpi_r(3,maxloc*nbeadslocmax,nproc_polymer))
      
        !PUT DATA INTO BUFFER ARRAYS
      
      DO k=1,ibeadend-ibeadbeg+1; DO iloc=1,nloc !; DO j=1,3
        iproc=min(nproc_polymer,int((iloc-1)/nlocproc)+1)
        if(rank_polymer==iproc-1) cycle
        ilocbeg=ilocpi_beg(iproc)
        ilocend=ilocpi_end(iproc)
        kk=(k-1)*(ilocend-ilocbeg+1)+ iloc-ilocbeg + 1
        bufferpi_s(:,kk,iproc)=derivs_ctr(:,iloc,k)
      ENDDO; ENDDO !; ENDDO

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
     &      ,MPI_REAL8
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
     &   ,MPI_REAL8,iproc-1,mptag_data, COMM_POLYMER
     &   ,reqs(ii),ierr)
      ENDDO

      !WAIT FOR COMMUNICATIONS
      call MPI_WAITALL(ii,reqs,MPI_STATUSES_IGNORE,ierr)

      !PUT DATA INTO POLYMER ARRAYS
      DO k=1,nu
        iproc=polymer_ctr%iproc(k)
        if(rank_polymer==iproc-1) cycle
        DO iloc=0,ilocend-ilocbeg
          ibeadbeg=polymer_ctr%ibead_beg(iproc)
          i=glob(ilocbeg+iloc)
          kk=(k-ibeadbeg)*nlocpi + iloc + 1
          polymer_ctr%forces_slow(:,i,k)=-bufferpi_r(:,kk,iproc)
        enddo
      ENDDO

      deallocate(bufferpi_s,bufferpi_r)

      if(pi_comm_report) then
        write(pi_comm_unit,*) "  exiting commforces_ctr"
        FLUSH(pi_comm_unit)
      endif

      end subroutine commforces_ctr

      end module commstuffpi
