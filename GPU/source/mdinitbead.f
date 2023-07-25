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
#include "tinker_macro.h"
      subroutine mdinitbead(dt,polymer)
      use atmtyp
      use atomsMirror
      use bath
      use beads
      use bound
      use couple
      use domdec
      use files
      use keys
      use freeze
      use inform
      use iounit
      use langevin
      use qtb
      use math
      use mdstuf
      use molcul
      use moldyn
      use mpole
      use neigh
      use units 
      use uprior
      use usage
      use random_mod
      use commstuffpi, only: comm_for_normal_modes

      implicit none
      integer :: ibead
      real(r_p), intent(in) :: dt
      type(POLYMER_COMM_TYPE), intent(inout) :: polymer
      integer i,j,k,idyn,iglob,nh
      integer next
      integer lext,freeunit
      integer ierr
      real(r_p) :: e
      real(r_p) :: maxwell,speed
      real(r_p) :: hmax,hmass
      real(r_p) :: sum,dmass
      real(t_p) :: vec(3)
      real(r_p), allocatable :: speedvec(:)
      real(r_p) :: wt
      logical exist,heavy,restart
      character*7 ext
      character*20 keyword
      character*120 dynfile
      character*120 record
      character*120 string
      character*3 numberbeads
     
c
c     try to restart using prior velocities and accelerations
c
      restart = .false.
      do ibead = 1,nbeads
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
        restart = exist

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
           a(:,:) = 0.0_re_p
           aalt(:,:) = 0.0_re_p
!$acc update device(a,aalt) async
!$acc data present(glob,mass,v,samplevec) async
c
#ifdef _OPENACC
          if (.not.host_rand_platform) then
             allocate(speedvec(nloc))
!$acc data create(speedvec) async
              call rand_unitgpu(samplevec(1),nloc)
              call maxwellgpu(mass,nbeads*kelvin,nloc,speedvec)
!$acc parallel loop collapse(2) async
              do i = 1, nloc
                 do j = 1, 3
                    iglob = glob(i)
                    if (use(iglob)) then
                       v(j,iglob) = speedvec(i) 
     &                            * real(samplevec(3*(i-1)+j),r_p)
                    else
                       v(j,iglob)    = 0.0_re_p
                    end if
                 end do
              end do
!$acc end data
              deallocate(speedvec)
           end if
#endif
          if (host_rand_platform) then
!$acc wait
!$acc update host(glob)
              do i = 1, nloc
                 iglob = glob(i)
                 if (use(iglob)) then
                    speed = maxwell (mass(iglob),nbeads*kelvin)
                    call ranvec (vec)
                    do j = 1, 3
                       v(j,iglob) = speed * real(vec(j),r_p)
                    end do
                 else
                    do j = 1, 3
                       v(j,iglob)    = 0.0_re_p
                    end do
                 end if
              end do
!$acc update device(v)
           end if
!$acc end data
        end if

!$acc wait
!$acc parallel loop async default(present) 
        do i=1,n
          iglob=glob(i)
          wt=mass(iglob)/convert
          polymer%pos(1,iglob,ibead)=x(iglob)
          polymer%pos(2,iglob,ibead)=y(iglob)
          polymer%pos(3,iglob,ibead)=z(iglob)
!$acc loop seq
          do j=1,3
            polymer%vel(j,iglob,ibead)=v(j,iglob)
            polymer%forces(j,iglob,ibead)=a(j,iglob)*wt
            polymer%forces_slow(j,iglob,ibead)=aalt(j,iglob)*wt
          enddo
        enddo
c
c     check for any prior dynamics coordinate sets
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
!$acc parallel loop default(present) async
        do iglob=1,n
          x(iglob)=0.d0
          y(iglob)=0.d0
          z(iglob)=0.d0
!$acc loop seq
          do ibead=1,nbeads
            x(iglob)=x(iglob)+polymer%pos(1,iglob,ibead)
            y(iglob)=y(iglob)+polymer%pos(2,iglob,ibead)
            z(iglob)=z(iglob)+polymer%pos(3,iglob,ibead)
          enddo
          x(iglob)=x(iglob)/real(nbeads,r_p)
          y(iglob)=y(iglob)/real(nbeads,r_p)
          z(iglob)=z(iglob)/real(nbeads,r_p)
        enddo
        call reCast_position
        call ddpme3d
        call reinitnl(0)
        call reassignpme(.false.)
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
!$acc exit data delete(gamma_friction)
          deallocate(gamma_friction)
        endif
        allocate(gamma_friction(nbeads))
        if(default_lambda_trpmd) lambda_trpmd = 1.0
        do ibead=1,nbeads
          gamma_friction(ibead) = max(gamma,lambda_trpmd*omkpi(ibead))
        enddo
!$acc enter data copyin(gamma_friction(:)) async
      endif


      end
c
