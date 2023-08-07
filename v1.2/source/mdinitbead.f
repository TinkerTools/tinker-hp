c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine mdinitbead  --  initialize a pimd trajectory  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "mdinitbead" initializes the velocities and positions
c     for a molecular pimd trajectory, including restarts
c
c
      subroutine mdinitbead(dt,polymer)
      use atmtyp
      use atoms
      use bath
      use beads
      use bound
      use commstuffpi
      use couple
      use domdec
      use files
      use keys
      use langevin, only: gamma_friction,gamma
      use freeze
      use inform
      use iounit
      use math
      use mdstuf
      use molcul
      use moldyn
      use mpole
      use neigh
      use qtb
      use units 
      use uprior
      use usage
      implicit none
      real*8, intent(in) :: dt
      type(POLYMER_COMM_TYPE), intent(inout) :: polymer
      integer i,j,k,idyn,iglob,nh,ibead
      integer next
      integer lext,freeunit
      integer ierr
      real*8 e,wt
      real*8 maxwell,speed
      real*8 hmax,hmass
      real*8 sum,dmass
      real*8 vec(3),Rn,Rnvar
      real*8, allocatable :: derivs(:,:)
      logical exist,heavy
      character*7 ext
      character*20 keyword
      character*240 dynfile
      character*240 record
      character*240 string
      character*3 numberbeads
      logical restart
      interface
         function normal ()
         real*8 normal
         end function
      end interface
     
c
c     try to restart using prior velocities and accelerations
c
      restart=.FALSE.
      do ibead=1,nbeads
        write(numberbeads, '(i3.3)') ibead
        dynfile = filename(1:leng)//'_beads'//numberbeads//'.dyn'
        call version (dynfile,'old')
        inquire (file=dynfile,exist=exist)

        if((.NOT.exist) .and. restart) then
          if(ranktot==0) then
            write(0,*) "Error: not all beads can been restarted!"
            call fatal
          endif
        endif
        restart=exist

        if (exist) then
           if(ranktot==0) then
            write(*,*) " --- Tinker-HP loading Restart file "
     &           //trim(dynfile)//"  ---"
           endif
           idyn = freeunit ()
           open (unit=idyn,file=dynfile,status='old')
           rewind (unit=idyn)
           call readdyn (idyn)
           close (unit=idyn)
        else
c 
c       set velocities and accelerations for cartesian dynamics
c 
c           allocate (derivs(3,nbloc))
c           derivs = 0d0
c          ! call allocstep
c           call gradient (epotpi_loc,derivs)
c           call commforces(derivs)
c 
           a(:,:) = 0.d0
           if(allocated(aalt)) aalt(:,:) = 0.d0
           do i = 1, nloc
              iglob = glob(i)
              if (use(iglob)) then
c                 speed = maxwell (mass(iglob),nbeads*kelvin)
c                 call ranvec (vec)
                 do j = 1, 3
                   v(j,iglob)=normal()*sqrt(nbeads
     &                  *kelvin*boltzmann/mass(iglob))
c                    v(j,iglob) = speed * vec(j)
c                   if(i==1 .AND. j==1) write(0,*) ranktot,v(j,iglob)
                 end do
              else
                 do j = 1, 3
                    v(j,iglob) = 0.0d0
                 end do
              end if
           end do
c           deallocate (derivs)
c           if (nuse .eq. n)  call mdrest (0)
        end if

        do iglob=1,n
          polymer%pos(1,iglob,ibead)=x(iglob)
          polymer%pos(2,iglob,ibead)=y(iglob)
          polymer%pos(3,iglob,ibead)=z(iglob)
          
          wt=mass(iglob)/convert
          polymer%vel(:,iglob,ibead)=v(:,iglob)
          polymer%forces(:,iglob,ibead)=a(:,iglob)*wt
        enddo

        if(allocated(aalt)) then
           do iglob=1,n
              wt=mass(iglob)/convert
              polymer%forces_slow(:,iglob,ibead)=aalt(:,iglob)*wt
           enddo
        else
           polymer%forces_slow(:,:,:)=0.d0
        endif
c 
c       check for any prior dynamics coordinate sets
c 
        i = 0
        exist = .true.
        do while (exist)
           i = i + 1
           lext = 3
           call numeral (i,ext,lext)
           dynfile = filename(1:leng)//'_beads'//numberbeads//
     $                          '.'//ext(1:lext)
           inquire (file=dynfile,exist=exist)
           if (.not.exist .and. i.lt.100) then
              lext = 2
              call numeral (i,ext,lext)
              dynfile = filename(1:leng)//'_beads'//numberbeads//
     $                          '.'//ext(1:lext)
              inquire (file=dynfile,exist=exist)
           end if
           if (.not.exist .and. i.lt.10) then
              lext = 1
              call numeral (i,ext,lext)
              dynfile = filename(1:leng)//'_beads'//numberbeads//
     $                          '.'//ext(1:lext)
              inquire (file=dynfile,exist=exist)
           end if
        end do
        nprior = i - 1
      enddo


      if(restart) then
        do iglob=1,n
          !iglob=glob(i)
          x(iglob)=0.d0
          y(iglob)=0.d0
          z(iglob)=0.d0
          do ibead=1,nbeads
            x(iglob)=x(iglob)+polymer%pos(1,iglob,ibead)
            y(iglob)=y(iglob)+polymer%pos(2,iglob,ibead)
            z(iglob)=z(iglob)+polymer%pos(3,iglob,ibead)
          enddo
          x(iglob)=x(iglob)/real(nbeads,8)
          y(iglob)=y(iglob)/real(nbeads,8)
          z(iglob)=z(iglob)/real(nbeads,8)
        enddo
        call ddpme3d
        call reinitnl(0)
        call reassignpme(.true.)
        call mechanicstep(0)
        call allocstep
        call nblist(0)
        call update_nlocpi(nloc)
      endif

      call comm_for_normal_modes(polymer,polymer%pos
     &   ,polymer%vel,polymer%forces,polymer%forces_slow )
      call update_normal_modes_pi(polymer)

      if(restart) then
        if(integrate.eq.'BAOABRESPA') then
          call set_eigforces_pi(polymer,polymer%forces_slow)
        else
          call set_eigforces_pi(polymer,polymer%forces)
        endif
      endif

      if(.not. piqtb) then
        if(allocated(gamma_friction)) then
          deallocate(gamma_friction)
        endif
        allocate(gamma_friction(nbeads))
        if(default_lambda_trpmd) lambda_trpmd = 1.0d0
        do ibead=1,nbeads
          gamma_friction(ibead) = max(gamma,lambda_trpmd*omkpi(ibead))
        enddo
      endif



      end
c
