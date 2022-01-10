c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #######################################################################################
c     ##                                                                                   ##
c     ##  subroutine baoabopi  --  BAOAB Path Integral Langevin molecular dynamics step    ##
c     ##                                                                                   ##
c     #######################################################################################
c
c
c     "baoabpi" performs a single molecular dynamics time step
c     via the BAOAB recursion formula
c
c     literature reference:
c
c     Efficient molecular dynamics using geodesic integration 
c     and solvent-solute splitting B. Leimkuhler and C. Matthews,
c     Proceedings of the Royal Society A, 472: 20160138, 2016
c
c
      subroutine baoabpi4site1 (istep,dt)
      use atmtyp
      use atoms
      use bath
      use cutoff
      use domdec
      use energi
      use freeze
      use langevin
      use mdstuf
      use moldyn
      use timestat
      use units
      use usage
      use mpi
      implicit none
      integer i,j,istep,iglob,ibead
      real*8 dt,dt_x,factor
      real*8 etot,eksum,epot
      real*8 temp,pres
      real*8 part1,part2
      real*8 a1,a2,normal
      real*8 ekin(3,3)
      real*8 stress(3,3)
      real*8 time0,time1
c      real*8, allocatable :: derivs(:,:)
c
c
c     find quarter step velocities and half step positions via BAOAB recursion
c
      do i = 1, nloc
         iglob = glob(i)
         if (atomic(iglob).eq.0) then
           v(:,iglob) = 0d0
         else if (use(iglob)) then
            do j = 1, 3
               v(j,iglob) = v(j,iglob) + 0.5*dt*a(j,iglob)
            end do
         end if
      end do
c
c      if (use_rattle) call rattle2(dt)
c
c      do i = 1, nloc
c        iglob = glob(i)
c        if (use(iglob)) then
c          xold(iglob) = x(iglob)
c          yold(iglob) = y(iglob)
c          zold(iglob) = z(iglob)
c          x(iglob) = x(iglob) + v(1,iglob)*0.5*dt
c          y(iglob) = y(iglob) + v(2,iglob)*0.5*dt
c          z(iglob) = z(iglob) + v(3,iglob)*0.5*dt
c        end if
c      end do
c
c      if (use_rattle) call rattle(0.5*dt,xold,yold,zold)
c      if (use_rattle) call rattle2(0.5*dt)
c
c      if (use_rattle) call rattle2(dt)
c
c     compute random part
c
c      deallocate (Rn)
c      allocate (Rn(3,nloc))
c      do i = 1, nloc
c        do j = 1, 3
c          Rn(j,i) = normal()
c        end do
c      end do
c      do i = 1, nloc
c         iglob = glob(i)
c         if (use(iglob)) then
c            do j = 1, 3
c               v(j,iglob) = a1*v(j,iglob) + 
c     $            a2*Rn(j,i)/sqrt(mass(iglob))
c            end do
c         end if
c      end do
c
c      if (use_rattle) call rattle2(dt)
c
c     find full step positions via BAOAB recursion
c     
c      do i = 1, nloc
c         iglob = glob(i)
c         if (use(iglob)) then
c            xold(iglob) = x(iglob)
c            yold(iglob) = y(iglob)
c            zold(iglob) = z(iglob)
c            x(iglob) = x(iglob) + v(1,iglob)*0.5*dt
c            y(iglob) = y(iglob) + v(2,iglob)*0.5*dt
c            z(iglob) = z(iglob) + v(3,iglob)*0.5*dt
c         end if
c      end do
c
c      if (use_rattle) call rattle(0.5*dt,xold,yold,zold)
c      if (use_rattle) call rattle2(0.5*dt)
c
c
c     make half-step temperature and pressure corrections
c
c      call temper (dt,eksum,ekin,temp)
c      call pressure2 (epot,temp)
c
      end

      subroutine baoabpi4site2 (istep,dt)
      use atmtyp
      use atoms
      use bath
      use beads
      use cutoff
      use deriv
      use domdec
      use energi
      use freeze
      use langevin
      use mdstuf
      use moldyn
      use neigh
      use timestat
      use units
      use usage
      use mpi
      implicit none
      integer i,j,istep,iglob,ibead,ierr
      integer iloc,iloc1,iloc2,iloc3,ilocrec1,ilocrec2,ilocrec3
      real*8 oterm,hterm
      real*8 dt,dt_x,factor
      real*8 pres
      real*8 part1,part2
      real*8 stress(3,3)
      real*8 time0,time1
      real*8, allocatable :: derivs(:,:)

      oterm = 0.73612d0
      hterm = 0.13194d0
c
c     move M center to linear combination of H-O-H coordinates
c
c      do i = 1, nlocnl
c         iglob = ineignl(i)
      do i = 1, nbloc
         iglob = glob(i)
         if (atomic(iglob) .eq. 0) then
            x(iglob) = oterm*x(iglob-3) + hterm*x(iglob-2) + 
     $           hterm*x(iglob-1)
            y(iglob) = oterm*y(iglob-3) + hterm*y(iglob-2) +
     $           hterm*y(iglob-1)
            z(iglob) = oterm*z(iglob-3) + hterm*z(iglob-2) + 
     $           hterm*z(iglob-1)
         end if
      end do
c
c     rebuild the neighbor lists
c
      if (use_list) call nblist(istep)

c      call pressure2 (epot,temp)
c
      allocate (derivs(3,nbloc))
      derivs = 0d0
c
c     get the potential energy and atomic forces
c
      call gradient (epotpi_loc,derivs)
c
c     transfer the forces on the M center to the H-O-H sites
c
      do i = 1, nlocnl
         iglob = ineignl(i)
         iloc = loc(iglob)
         if (atomic(iglob) .eq. 0) then
            iloc1 = loc(iglob-1)
            iloc2 = loc(iglob-2)
            iloc3 = loc(iglob-3)
            do j = 1, 3
               derivs(j,iloc3) = derivs(j,iloc3) + oterm*derivs(j,iloc)
               derivs(j,iloc2) = derivs(j,iloc2) + hterm*derivs(j,iloc)
               derivs(j,iloc1) = derivs(j,iloc1) + hterm*derivs(j,iloc)
            end do
         end if
      end do
      do i = 1, nlocrec
         iglob = globrec(i)
         iloc = locrec1(iglob)
         if (atomic(iglob) .eq. 0) then
            iloc1 = locrec1(iglob-1)
            iloc2 = locrec1(iglob-2)
            iloc3 = locrec1(iglob-3)
            do j = 1, 3
               decrec(j,iloc3) = decrec(j,iloc3) + oterm*decrec(j,iloc)
               decrec(j,iloc2) = decrec(j,iloc2) + hterm*decrec(j,iloc)
               decrec(j,iloc1) = decrec(j,iloc1) + hterm*decrec(j,iloc)
            end do
         end if
      end do
c
c     MPI : get total energy (and allreduce the virial)
c
      call reduceen(epotpi_loc)
      call MPI_BCAST(epotpi_loc,1,MPI_REAL8,0,COMM_BEAD,ierr)
c
c     communicate forces
c
      call commforces(derivs)
c
c     use Newton's second law to get the next accelerations;
c     find the full-step velocities using the BAOAB recursion
c
c      write(*,*) 'x 1 = ',x(1),y(1),v(1,1),a(1,1)
      do i = 1, nloc
         iglob = glob(i)
         if (atomic(iglob) .eq. 0) then
            do j = 1, 3
               v(j,iglob) = 0.0d0
               a(j,iglob) = 0.0d0
            end do
         else if (use(iglob)) then
            do j = 1, 3
               a(j,iglob) = -convert * derivs(j,i)/mass(iglob)
               v(j,iglob) = v(j,iglob) + 0.5*dt*a(j,iglob)
            end do
         end if
      end do
c      write(*,*) 'x 2 = ',x(1),y(1),v(1,1),a(1,1)
c
c     perform deallocation of some local arrays
c
      deallocate(derivs)
c
c     find the constraint-corrected full-step velocities
c
c      if (use_rattle)  call rattle2 (dt)
c
c      call temper (dt,eksum,ekin,temp)
c     call pressure (dt,epot,ekin,temp,pres,stress,istep)
      call kinetic (eksumpi_loc,ekinpi_loc,temppi)
c
c     total energy is sum of kinetic and potential energies
c
      etotpi_loc = eksumpi_loc + epotpi_loc

c
c     compute statistics and save trajectory for this step
c
c      call mdstat (istep,dt,etotpi,epotpi,eksumpi,temppi,pres)
      call mdsavebeads (istep,dt)
c      call mdrest (istep)
      return
      end

c     !subroutine commposbead: communicate all positions to global master to be able
c     to make normal mode transformation and propagation,
c     then communicate back the positions after NM propagation
c
      subroutine aoapi4site(istep,dt)
      use atoms
      use boxes, only: volbox,xbox,ybox,zbox
      use domdec
      use mpi
      use units
      use beads
      use langevin 
      use bath
      use cutoff
      use energi
      use freeze
      use mdstuf
      use moldyn
      use qtb
      use timestat
      use units
      use usage
      use atmtyp
      use math
      use iounit
      use inform
      use molcul
      implicit none
      real*8, intent(in) :: dt
      integer, intent(in) :: istep
      integer, allocatable :: reqsend(:),reqrec(:),nbeadscomm(:)
      integer, allocatable :: nloccomm(:)
      real*8, allocatable :: indexposcomm(:,:,:,:),buffer(:,:,:)
      integer i,j,k,tagmpi,status(MPI_STATUS_SIZE),ierr,iproc, ibead,l
      integer :: modstep
      !real*8, allocatable :: eigmat(:,:),eigmattr(:,:),omkpi(:),WORK(:)
      real*8 :: eigx0,eigv0
      real*8 :: dt2, a1,a2, sqrtmass
      real*8, allocatable :: pos_full(:,:,:),vel_full(:,:,:)
      real*8, allocatable :: forces_full(:,:,:)
      real*8 :: normal
      !real*8 :: hbar
      real*8 :: buffer_energy(3)
      real*8 :: ekprim,ekvir,dedv_mean,presvir
      real*8 :: nrealdof
      real*8 :: mpitime1, mpitime2
      real*8 :: time_tot, time_com
      real*8 :: acentroid,a1p,a2p,scale,third,gammak
      real*8 :: dens

      allocate (reqsend(nproctot))
      allocate (reqrec(nproctot))

      !hbar=(planck*1.d11*avogadro)/(2.*pi)

      time_tot=0.d0
      time_com=0.d0
c
c     reduce potential and kinetic energies
c
      mpitime1=mpi_wtime()
      buffer_energy=0
      DO ibead=1,nbeadsloc
        buffer_energy(1)=buffer_energy(1) + epotpi(ibead)
        buffer_energy(2)=buffer_energy(2) + eksumpi(ibead)
        buffer_energy(3)=buffer_energy(3) + dedv_pi(ibead)
      ENDDO
      if (ranktot.eq.0) then
        
        call MPI_REDUCE(MPI_IN_PLACE,buffer_energy,3,MPI_REAL8
     $     ,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        epotpi_loc=buffer_energy(1)/(nbeads*nproc)
        eksumpi_loc=buffer_energy(2)/(nbeads*nproc)
        dedv_mean=buffer_energy(3)/(nbeads*nproc)

      else
        call MPI_REDUCE(buffer_energy,buffer_energy,3,MPI_REAL8
     $     ,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      end if

c
c     first get the number of beads per process
c
      if (ranktot.eq.0) allocate(nbeadscomm(nproctot))
   
      if (ranktot.eq.0) then
        do i = 1, nproctot-1
          tagmpi = nproctot*ranktot + i + 1
          call MPI_IRECV(nbeadscomm(i+1),1,
     $     MPI_INT,i,tagmpi,MPI_COMM_WORLD,reqrec(i),ierr)
          call MPI_WAIT(reqrec(i),status,ierr)
c          write(*,*) 'nbeads of ',i,' = ',nbeadscomm(i+1)
        end do
      else
        tagmpi = ranktot + 1
        call MPI_ISEND(nbeadsloc,1,MPI_INT,0,tagmpi,MPI_COMM_WORLD,
     $   reqsend(1),ierr)
        call MPI_WAIT(reqsend(1),status,ierr)
      end if
c
c     get the number of atoms per process
c
      if (ranktot.eq.0) allocate(nloccomm(nproctot))
   
      if (ranktot.eq.0) then
        do i = 1, nproctot-1
          tagmpi = nproctot*ranktot + i + 1
          call MPI_IRECV(nloccomm(i+1),1,
     $     MPI_INT,i,tagmpi,MPI_COMM_WORLD,reqrec(i),ierr)
          call MPI_WAIT(reqrec(i),status,ierr)
c          write(*,*) 'nloc of ',i,' = ',nloccomm(i+1)
        end do
      else
        tagmpi = ranktot + 1
c
c     the proc having more than 1 beads have necessarily the whole system for all these beads
c
        call MPI_ISEND(nlocpi(1),1,MPI_INT,0,tagmpi,MPI_COMM_WORLD,
     $   reqsend(1),ierr)
        call MPI_WAIT(reqsend(1),status,ierr)
      end if
c
c     get their indexes and positions, velocities and forces
c
      if (ranktot.eq.0) then
        allocate(indexposcomm(10,n,nbeads,nproctot))
      else
        allocate(buffer(10,nlocpi(1),nbeadsloc))
      end if
   
      if (ranktot.eq.0) then
        do i = 1, nproctot-1
          tagmpi = nproctot*ranktot + i + 1
          call MPI_IRECV(indexposcomm(1,1,1,i+1),
     $     10*nloccomm(i+1)*nbeadscomm(i+1),
     $     MPI_REAL8,i,tagmpi,MPI_COMM_WORLD,reqrec(i),ierr)
          call MPI_WAIT(reqrec(i),status,ierr)
c       do k = 1, nbeadscomm(i+1)
c         do j = 1, nloccomm(i+1)
c           write(*,*) 'beads ',k,'indexes of ',i,' = ',
c  $      indexposcomm(1,j,k,i+1),
c  $      indexposcomm(2,j,k,i+1),indexposcomm(3,j,k,i+1),
c  $      indexposcomm(4,j,k,i+1)
c         end do
c       end do
        end do
      else
        tagmpi = ranktot + 1
c
c     the proc having more than 1 beads have necessarily the whole system for all these beads
c
        do k = 1, nbeadsloc
          do j = 1, nlocpi(1)
            buffer(1,j,k) = real(globpi(j,k))

            buffer(2,j,k) = pospi(1,glob(j),k)
            buffer(3,j,k) = pospi(2,glob(j),k)
            buffer(4,j,k) = pospi(3,glob(j),k)

            buffer(5,j,k) = velpi(1,glob(j),k)
            buffer(6,j,k) = velpi(2,glob(j),k)
            buffer(7,j,k) = velpi(3,glob(j),k)

            !send forces
            buffer(8,j,k) = api(1,glob(j),k)*mass(glob(j))/convert
            buffer(9,j,k) = api(2,glob(j),k)*mass(glob(j))/convert
            buffer(10,j,k) = api(3,glob(j),k)*mass(glob(j))/convert
          end do
        end do
        call MPI_ISEND(buffer,10*nlocpi(1)*nbeadsloc,MPI_REAL8,0,tagmpi,
     $   MPI_COMM_WORLD,reqsend(1),ierr)
        call MPI_WAIT(reqsend(1),status,ierr)
      end if
      
      mpitime2=mpi_wtime()

      time_com=mpitime2-mpitime1

      !if (ranktot/=0) then
       ! deallocate(buffer)
      !end if
c
c
c     DO THE NM PROPAGATION
c
c      indexposcomm = 0d0
c
c     propagate AOA
      if(ranktot .eq. 0) then

c       ORDER POSITIONS AND VELOCITES IN pos_full AND vel_full
          allocate(pos_full(3,n,nbeads),vel_full(3,n,nbeads))
          allocate(forces_full(3,n,nbeads))
          DO l=1,nbeadsloc
            DO i=1,nloc ; DO j=1,3
              k=glob(i)
              pos_full(j,k,l) =pospi(j,k,l)
              vel_full(j,k,l) =velpi(j,k,l)
              forces_full(j,k,l) =api(j,k,l)*mass(k)/convert
            ENDDO  ; ENDDO
          ENDDO
        IF(nproctot>nbeads) THEN
c           EACH PROC TREATS ONLY 1 BEAD (OR SPATIAL PART OF IT)
          ibead=1
          do iproc = 1, nproctot-1
            if(mod(iproc,nproc)==0) ibead=ibead+1
            DO i=1,nloccomm(iproc+1) ; DO j=1,3
              k=nint(indexposcomm(1,i,1,iproc+1))
              pos_full(j,k,ibead)=indexposcomm(j+1,i,1,iproc+1)  
              vel_full(j,k,ibead)=indexposcomm(j+4,i,1,iproc+1)  
              forces_full(j,k,ibead)=indexposcomm(j+7,i,1,iproc+1)  
            ENDDO  ; ENDDO
          enddo
        ELSE
c           EACH PROC CAN TREAT MULTIPLE BEADS BUT NO SPATIAL PARALELLIZATION
          ibead=nbeadsloc
          do iproc = 1, nproctot-1
            DO l=1,nbeadscomm(iproc+1)
              ibead=ibead+1
              DO i=1,nloccomm(iproc+1) ; DO j=1,3
                k=nint(indexposcomm(1,i,l,iproc+1))
                pos_full(j,k,ibead)= indexposcomm(j+1,i,l,iproc+1)  
                vel_full(j,k,ibead)= indexposcomm(j+4,i,l,iproc+1)  
                forces_full(j,k,ibead)= indexposcomm(j+7,i,l,iproc+1)  
              ENDDO  ; ENDDO
            ENDDO         
          enddo
        ENDIF

        !deallocate(indexposcomm)

      call compute_observables_pi_4site(pos_full,forces_full
     &          ,dedv_mean,Ekcentroid,ekprim,ekvir,presvir)
    

        temppi = 2.0d0 * ekvir / (nfree * gasconst)
        temppi_cl = 2.0d0 * eksumpi_loc / (nfree * gasconst)
c        dens = (1.0d24/extvol) * (totmass/avogadro)
c
c        modstep = mod(istep,iprint)
c        if (verbose) then
c          if(isobaric) then
c            if (modstep .eq. 1) then
c                write (iout,10)
c10            format (/,4x,'MD Step',6x,'E Total',3x,'E Potential'
c     &          ,7x,'Ek vir',7x,'Ek prim',7x,'Temp',5x,'Temp_cl'
c     &          ,5x,'Pres',7x,'Density',7x,'Volume'/)
c            end if
c            write (iout,40)  istep,ekvir+epotpi_loc,epotpi_loc,ekvir
c     &             ,ekprim,temppi,temppi_cl,presvir,dens,extvol
c40         format (i10,4f14.4,3f11.2,f10.4,f14.4)
c          else
c            if (modstep .eq. 1) then
c                write (iout,60)
c60            format (/,4x,'MD Step',6x,'E Total',3x,'E Potential'
c     &          ,7x,'Ek vir',7x,'Ek prim',7x,'Temp',5x,'Temp_cl'
c     &          ,5x,'Pres'/)
c            end if
c            write (iout,50)  istep,ekvir+epotpi_loc,epotpi_loc,ekvir
c     &             ,ekprim,temppi,temppi_cl,presvir
c50         format (i10,4f14.4,3f11.2)
c          endif
c        end if
      call mdstatpi(istep,dt,ekvir,ekprim,presvir)

c       COMPUTE POLYMER TRANSFER MATRIX (TO BE REPLACED BY FFT !!!)
!        allocate(eigmat(nbeads,nbeads),eigmattr(nbeads,nbeads))
!        allocate(omkpi(nbeads),WORK(3*nbeads))
!        eigmat=0
!        DO i=1,nbeads-1
!          eigmat(i,i)=2
!          eigmat(i+1,i)=-1
!          eigmat(i,i+1)=-1
!        ENDDO
!        eigmat(1,nbeads)=-1
!        eigmat(nbeads,1)=-1
!        eigmat(nbeads,nbeads)=2
!        call DSYEV('V','U',nbeads,eigMat,nbeads, 
!     $      omkpi,WORK,3*nbeads,ierr)		
!        omkpi(1)=0
!        omkpi(:)=sqrt(omkpi)*(nbeads*boltzmann*kelvin/hbar)
!        eigmattr=transpose(eigmat)
!        deallocate(WORK)

        dt2=0.5*dt
        a1 = exp(-gamma*dt)
        a2 = sqrt((1.-a1**2)*nbeads*boltzmann*kelvin)
c       transform to normal modes
        DO i=1,n ; if(atomic(i)==0) CYCLE ; DO j=1,3
          vel_full(j,i,:)=matmul(eigmattr,vel_full(j,i,:))
          pos_full(j,i,:)=matmul(eigmattr,pos_full(j,i,:))
        ENDDO ; ENDDO

        
        if (ir) then
          k = mod(istep-1,nseg)+1
          do i=1,n
            do j=1,3
              vad(j,i,k)=vel_full(j,i,1)/sqrt(real(nbeads))
            enddo
          enddo 
c              write(*,*) 'compteur=', compteur
c              write(*,*) 'nseg=', nseg
c              write(*,*) 'Tseg=', Tseg
          if ((mod(istep,nseg).eq.0)) then
c             write(*,*) 'compteur=', compteur
c             write(*,*) 'skipseg=', skipseg
              call irspectra_pimd
              compteur=compteur+1
          endif
        endif

        if(isobaric) then
c         propagate volume velocity from previous step
          aextvol = 3.d0*nbeads*convert*(
     &        extvol*(presvir-atmsph)/prescon 
     &        +gasconst*kelvin 
     &      )/masspiston
          vextvol = vextvol + dt2*aextvol
          
c         update Ekcentroid and propagate volume velocity           
          presvir=presvir-prescon*(Ekcentroid/(3*nbeads*volbox))
          Ekcentroid=0
          DO i=1,n ;if(atomic(i)==0) CYCLE ; DO j=1,3          
            Ekcentroid=Ekcentroid+mass(i)*vel_full(j,i,1)**2
          ENDDO ; ENDDO
          Ekcentroid=Ekcentroid/convert
          
          !write(*,*) Ekcentroid,nfree*nbeads*kelvin*gasconst
          
          presvir=presvir+prescon*(Ekcentroid/(3*nbeads*volbox))
          aextvol = 3.d0*nbeads*convert*(
     &        volbox*(presvir-atmsph)/prescon 
     &        +gasconst*kelvin 
     &      )/masspiston
          vextvol = vextvol + dt2*aextvol

c         propagate centroid isobaric (half step)
          acentroid = sinh(dt2*vextvol)/vextvol
          scale=exp(dt2*vextvol)  
          extvolold = extvol
          extvol = extvol*exp(3.d0*dt2*vextvol)
          DO i=1,n ; if(atomic(i)==0) CYCLE ; DO j=1,3
            pos_full(j,i,1)=pos_full(j,i,1)*scale
     &        +acentroid*vel_full(j,i,1)
            vel_full(j,i,1)=vel_full(j,i,1)/scale
          ENDDO ; ENDDO 
        else
c         propagate centroid (half step)
          DO i=1,n ; if(atomic(i)==0) CYCLE ; DO j=1,3
            pos_full(j,i,1)=pos_full(j,i,1) + dt2*vel_full(j,i,1)
          ENDDO ; ENDDO
        endif

c         propagate springs (half step)
        DO k=2,nbeads
          DO i=1,n ; if(atomic(i)==0) CYCLE ; DO j=1,3           
            eigx0=pos_full(j,i,k)
            eigv0=vel_full(j,i,k)
            pos_full(j,i,k)=eigx0*cos(omkpi(k)*dt2)
     $            +eigv0*sin(omkpi(k)*dt2)/omkpi(k)
            vel_full(j,i,k)=eigv0*cos(omkpi(k)*dt2)
     $           -eigx0*sin(omkpi(k)*dt2)*omkpi(k)
          ENDDO ; ENDDO
        ENDDO

c         propagate langevin (full step) (all beads same gamma for now !)
        DO k=1,nbeads
          gammak=max(gamma,omkpi(k))
          a1 = exp(-gammak*dt)
          a2 = sqrt((1.-a1**2)*nbeads*boltzmann*kelvin)
          DO i=1,n ; if(atomic(i)==0) CYCLE ; DO j=1,3 
            vel_full(j,i,k)=vel_full(j,i,k)*a1 
     &              + a2*normal()/sqrt(mass(i))
          ENDDO; ENDDO
        ENDDO

        if(isobaric) then
c         langevin piston (full step)
          a1p = exp(-gammapiston*dt)
          a2p = sqrt((1-a1p**2)*nbeads*boltzmann*kelvin/masspiston)
          vextvol = vextvol*a1p + a2p*normal()

c         propagate centroid isobaric (half step)
          acentroid = sinh(dt2*vextvol)/vextvol
          scale=exp(dt2*vextvol)  
          extvol = extvol*exp(3.d0*dt2*vextvol)
          !volbox=extvol
          
          DO i=1,n ; if(atomic(i)==0) CYCLE ; DO j=1,3 
            pos_full(j,i,1)=pos_full(j,i,1)*scale
     &        +acentroid*vel_full(j,i,1)
            vel_full(j,i,1)=vel_full(j,i,1)/scale
          ENDDO ; ENDDO 
        else
c         propagate centroid (half step)
          DO i=1,n ; if(atomic(i)==0) CYCLE ; DO j=1,3 
            pos_full(j,i,1)=pos_full(j,i,1) + dt2*vel_full(j,i,1)
          ENDDO ; ENDDO
        endif

c         propagate springs (half step)
        DO k=2,nbeads
          DO i=1,n ; if(atomic(i)==0) CYCLE ; DO j=1,3 
            eigx0=pos_full(j,i,k)
            eigv0=vel_full(j,i,k)
            pos_full(j,i,k)=eigx0*cos(omkpi(k)*dt2)
     $            +eigv0*sin(omkpi(k)*dt2)/omkpi(k)
            vel_full(j,i,k)=eigv0*cos(omkpi(k)*dt2)
     $           -eigx0*sin(omkpi(k)*dt2)*omkpi(k)
          ENDDO ; ENDDO
        ENDDO
        
        Ekcentroid=0
        DO i=1,n ; if(atomic(i)==0) CYCLE ;DO j=1,3          
          Ekcentroid=Ekcentroid+mass(i)*vel_full(j,i,1)**2
        ENDDO ; ENDDO
        Ekcentroid=Ekcentroid/convert

c         transform back to coordinates
        DO i=1,n ; if(atomic(i)==0) CYCLE ; DO j=1,3 
          vel_full(j,i,:)=matmul(eigmat,vel_full(j,i,:))
          pos_full(j,i,:)=matmul(eigmat,pos_full(j,i,:))          
        ENDDO ; ENDDO

c       PUT BACK POSITIONS AND VELOCITES AT CORRECT INDICES
       ! allocate(indexposcomm(7,n,nbeads,nproctot))
        DO l=1,nbeadsloc
          DO i=1,nloc ; DO j=1,3
                k=glob(i)
                if(atomic(k)==0) cycle
                pospi(j,k,l)=pos_full(j,k,l)
                velpi(j,k,l)=vel_full(j,k,l)
          ENDDO  ; ENDDO
        ENDDO
        IF(nproctot>nbeads) THEN
c           EACH PROC TREATS ONLY 1 BEAD (OR SPATIAL PART OF IT)
          ibead=1
          do iproc = 1, nproctot-1
            if(mod(iproc,nproc)==0) ibead=ibead+1
            DO i=1,nloccomm(iproc+1) ; DO j=1,3
              k=nint(indexposcomm(1,i,1,iproc+1))
              indexposcomm(j+1,i,1,iproc+1)=pos_full(j,k,ibead) 
              indexposcomm(j+4,i,1,iproc+1) =vel_full(j,k,ibead)
            ENDDO  ; ENDDO
          enddo
        ELSE
c           EACH PROC CAN TREAT MULTIPLE BEADS BUT NO SPATIAL PARALELLIZATION
          ibead=nbeadsloc
          do iproc = 1, nproctot-1
            DO l=1,nbeadscomm(iproc+1)
              ibead=ibead+1
              DO i=1,nloccomm(iproc+1) ; DO j=1,3
                k=nint(indexposcomm(1,i,l,iproc+1))
                indexposcomm(j+1,i,l,iproc+1)=pos_full(j,k,ibead)
                indexposcomm(j+4,i,l,iproc+1)=vel_full(j,k,ibead)
              ENDDO  ; ENDDO
            ENDDO         
          enddo
        ENDIF
        deallocate(pos_full,vel_full,forces_full)
      endif

      mpitime1=mpi_wtime()

      time_tot=mpitime1-mpitime2


!      if (ranktot/=0) then
!         allocate(buffer(7,nlocpi(1),nbeadsloc))
!       end if
c
c     communicate back positions
c
      if (ranktot.eq.0) then
        do i = 1, nproctot-1
          tagmpi = nproctot*ranktot + i + 1
          call MPI_ISEND(indexposcomm(1,1,1,i+1),10*nloccomm(i+1)*
     $     nbeadscomm(i+1),MPI_REAL8,i,tagmpi,MPI_COMM_WORLD,reqsend(i),
     $     ierr)
          call MPI_WAIT(reqsend(i),status,ierr)
        end do
      else
        tagmpi = ranktot + 1
c
c     the proc having more than 1 beads have necessarily the whole system for all these beads
c
        call MPI_IRECV(buffer,10*nlocpi(1)*nbeadsloc,MPI_REAL8,0,tagmpi,
     $   MPI_COMM_WORLD,reqrec(1),ierr)
        call MPI_WAIT(reqrec(1),status,ierr)
        do k = 1, nbeadsloc
          do j = 1, nlocpi(1)
            pospi(1,glob(j),k) = buffer(2,j,k)
            pospi(2,glob(j),k) = buffer(3,j,k)
            pospi(3,glob(j),k) = buffer(4,j,k)
            velpi(1,glob(j),k) = buffer(5,j,k)
            velpi(2,glob(j),k) = buffer(6,j,k)
            velpi(3,glob(j),k) = buffer(7,j,k)
           ! write(*,*) 'j = ',glob(j),'k = ',k,'x = ',pospi(1,glob(j),k)
          end do
        end do
      end if
c
      if (ranktot.eq.0) then
        deallocate(nbeadscomm)
        deallocate(nloccomm)
        deallocate(indexposcomm)
      else
        deallocate(buffer)
      end if
      
      if(isobaric) then 
        if (ranktot.ne.0) extvolold=extvol
c       broadcast new volume  
        call MPI_BCAST(extvol,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
        !call MPI_BARRIER(MPI_COMM_WORLD,ierr)
c       rescale box
        third=1.d0/3.d0
        scale=(extvol/extvolold)**third
        xbox = (extvol)**third
        ybox = (extvol)**third
        zbox = (extvol)**third  
        
        !xbox = xbox*scale
        !ybox = ybox*scale
        !zbox = zbox*scale
        !call MPI_BARRIER(MPI_COMM_WORLD,ierr)
c
c     propagate the new box dimensions to other lattice values
c 
        call lattice  
        !call MPI_BARRIER(MPI_COMM_WORLD,ierr)   
        call ddpme3dnpt(scale,istep)
      endif
      
      mpitime2=mpi_wtime()

c      if (ranktot.eq.0) then
c            write(*,*) 'Time spent in com', mpitime2-mpitime1
c      endif
      end subroutine
