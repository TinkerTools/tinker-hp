c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ###################################################################
c     ##                                                               ##
c     ##  program pimd  --  run pimd molecular or stochastic dynamics  ##
c     ##                                                               ##
c     ###################################################################
c
c
c     "pimd" computes a molecular dynamics trajectory
c     in one of the standard statistical mechanical ensembles and using
c     any of several possible integration methods
c
c
      program pimd
      use mpi
      implicit none
      integer ierr,nthreadsupport
c      call MPI_INIT(ierr)
      call MPI_INIT_THREAD(MPI_THREAD_MULTIPLE,nthreadsupport,ierr)
      call pimd_bis
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      call MPI_FINALIZE(ierr)
      end
c
      subroutine pimd_bis
      use atoms
      use bath
      use beads
      use bound
      use boxes
      use domdec
      use keys
      use inform
      use iounit
      use mdstuf
      use commstuffpi
      use moldyn
      use neigh
      use timestat
      use mpi
      use cutoff
      use ascii, only: int_to_str
      use qtb, only: adaptive_qtb,piqtb
      use units
      implicit none
      integer istep,nstep,ierr
      integer mode
      real*8 dt,dtdump,time0,time1,bufpi
      logical exist,query
      character*240 string
      logical pi_isobaric
      logical only_long
      TYPE(POLYMER_COMM_TYPE), save :: polymer,polymer_ctr
      INTERFACE
        subroutine baoabpi(polymer,polymer_ctr,istep,dt)
          import :: POLYMER_COMM_TYPE 
          TYPE(POLYMER_COMM_TYPE), intent(inout) :: polymer,polymer_ctr
          integer, intent(in) :: istep
          real*8, intent(in) :: dt
        end subroutine baoabpi
        subroutine baoabrespapi(polymer,polymer_ctr,istep,dt)
          import :: POLYMER_COMM_TYPE 
          TYPE(POLYMER_COMM_TYPE), intent(inout) :: polymer,polymer_ctr
          integer, intent(in) :: istep
          real*8, intent(in) :: dt
        end subroutine baoabrespapi
      END INTERFACE

      path_integral_md=.TRUE.
c
c
 1000 Format(' Time for 100 Steps: ',f15.4,/,
     $  ' Ave. Time per step: ',f15.4)
 1010 Format(' ns per day: ',f15.4)
c
c     set up the structure and molecular mechanics calculation
c
      call initial
      call init_keys
      call read_bead_config
      call getxyz
      call initmpi
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      if(nproctot==1) pi_comm_report=.FALSE.
      if(pi_comm_report) then
        open(newunit=pi_comm_unit,
     &    file="comm_report_rank"//int_to_str(ranktot)//".out")
      endif

      call cutoffs
      call unitcell
      call lattice
c
c     get parameters
c
      call mechanic
c
c     setup for MPI
c
c
c     get nprocforce to do correct allocation of MPI arrays
c
      call drivermpi
      call reinitnl(0)
      call mechanic_init_para
c
c     allocate some arrays
c
      if (allocated(v)) deallocate (v)
      allocate (v(3,n))
      if (allocated(a)) deallocate (a)
      allocate (a(3,n))
      if (allocated(aalt)) deallocate (aalt)
      allocate (aalt(3,n))
c
      a = 0d0
      v = 0d0
      aalt = 0d0
c
      call nblist(0)
c
c     initialize the temperature, pressure and coupling baths
c
      kelvin = 0.0d0
      atmsph = 0.0d0
      isothermal = .false.
      isobaric = .false.
c
c     initialize the simulation length as number of time steps
c
      query = .true.
      call nextarg (string,exist)
      if (exist) then
         read (string,*,err=10,end=10)  nstep
         query = .false.
      end if
   10 continue
      if (query) then
         write (iout,20)
   20    format (/,' Enter the Number of Dynamics Steps to be',
     &              ' Taken :  ',$)
         read (input,30)  nstep
   30    format (i10)
      end if

c
c     get the length of the dynamics time step in picoseconds
c
      dt = -1.0d0
      call nextarg (string,exist)
      if (exist)  read (string,*,err=40,end=40)  dt
   40 continue
      do while (dt .lt. 0.0d0)
         write (iout,50)
   50    format (/,' Enter the Time Step Length in Femtoseconds',
     &              ' [1.0] :  ',$)
         read (input,60,err=70)  dt
   60    format (f20.0)
         if (dt .le. 0.0d0)  dt = 1.0d0
   70    continue
      end do
      dt = 0.001d0 * dt
c
c     enforce bounds on thermostat and barostat coupling times
c
      tautemp = max(tautemp,dt)
      taupres = max(taupres,dt)
c
c     set the time between trajectory snapshot coordinate dumps
c
      dtdump = -1.0d0
      call nextarg (string,exist)
      if (exist)  read (string,*,err=80,end=80)  dtdump
   80 continue
      do while (dtdump .lt. 0.0d0)
         write (iout,90)
   90    format (/,' Enter Time between Dumps in Picoseconds',
     &              ' [0.1] :  ',$)
         read (input,100,err=110)  dtdump
  100    format (f20.0)
         if (dtdump .le. 0.0d0)  dtdump = 0.1d0
  110    continue
      end do
      iwrite = nint(dtdump/dt)
c
c     get choice of statistical ensemble for periodic system
c
      if (use_bounds) then
         mode = -1
         call nextarg (string,exist)
         if (exist)  read (string,*,err=120,end=120)  mode
  120    continue
         do while (mode.lt.1 .or. mode.gt.4)
            write (iout,130)
  130       format (/,' Available Statistical Mechanical Ensembles :',
     &              //,4x,'(1) Microcanonical (NVE)',
     &              /,4x,'(2) Canonical (NVT)',
     &              /,4x,'(3) Isoenthalpic-Isobaric (NPH)',
     &              /,4x,'(4) Isothermal-Isobaric (NPT)',
     &              //,' Enter the Number of the Desired Choice',
     &                 ' [1] :  ',$)
            read (input,140,err=150)  mode
  140       format (i10)
            if (mode .le. 0)  mode = 1
  150       continue
         end do
         if (mode.eq.2 .or. mode.eq.4) then
            isothermal = .true.
            kelvin = -1.0d0
            call nextarg (string,exist)
            if (exist)  read (string,*,err=170,end=170)  kelvin
  170       continue
            do while (kelvin .lt. 0.0d0)
               write (iout,180)
  180          format (/,' Enter the Desired Temperature in Degrees',
     &                    ' K [298] :  ',$)
               read (input,190,err=200)  kelvin
  190          format (f20.0)
               if (kelvin .le. 0.0d0)  kelvin = 298.0d0
  200          continue
            end do
         end if
        if (mode.eq.1 .or. mode.eq.3) then
            if(ranktot.eq.0) then
               write(iout,*) 'NVE, NPH and NPT not implemented yet!'
               call fatal
            endif
        endif
        if(mode.eq.4) then
            isobaric = .true.
            atmsph = -1.0d0
            call nextarg (string,exist)
            if (exist)  read (string,*,err=210,end=210)  atmsph
  210       continue
            do while (atmsph .lt. 0.0d0)
               write (iout,220)
  220          format (/,' Enter the Desired Pressure in Atm',
     &                    ' [1.0] :  ',$)
               read (input,230,err=240)  atmsph
  230          format (f20.0)
               if (atmsph .le. 0.0d0)  atmsph = 1.0d0
  240          continue
            end do
         end if
      end if
c
c     use constant energy or temperature for nonperiodic system
c
      if (.not. use_bounds) then
         mode = -1
         call nextarg (string,exist)
         if (exist)  read (string,*,err=250,end=250)  mode
  250    continue
         do while (mode.lt.1 .or. mode.gt.2)
            write (iout,260)
  260       format (/,' Available Simulation Control Modes :',
     &              //,4x,'(1) Constant Total Energy Value (E)',
     &              /,4x,'(2) Constant Temperature via Thermostat (T)',
     &              //,' Enter the Number of the Desired Choice',
     &                 ' [1] :  ',$)
            read (input,270,err=280)  mode
  270       format (i10)
            if (mode .le. 0)  mode = 1
  280       continue
         end do
         if (mode .eq. 2) then
            isothermal = .true.
            kelvin = -1.0d0
            call nextarg (string,exist)
            if (exist)  read (string,*,err=290,end=290)  kelvin
  290       continue
            do while (kelvin .lt. 0.0d0)
               write (iout,300)
  300          format (/,' Enter the Desired Temperature in Degrees',
     &                    ' K [298] :  ',$)
               read (input,310,err=320)  kelvin
  310          format (f20.0)
               if (kelvin .le. 0.0d0)  kelvin = 298.0d0
  320          continue
            end do
         end if
      end if
c
c     setup dynamics
c
      call reset_observables(.true.)
      call allocpi(polymer,polymer_ctr)     
      if(centroid_longrange) then
        use_shortmlist = use_mlist
        use_shortclist = use_clist
        use_shortvlist = use_vlist
      endif

      if(add_buffer_pi) then
        bufpi = 2.d0*hbar_planck/sqrt(12.d0*boltzmann*kelvin)
        call update_lbuffer(lbuffer+bufpi)
        write(*,*) 'Added buffer of ',bufpi,' Angstroms'
        write(*,*) 'The buffer is now ',lbuffer,' Angstroms'
        write(*,*) 'with nlupdate = ',ineigup
      endif

      call mdinit(dt)    
      call mdinitbead(polymer)

      !if(qtb_thermostat) then
      !  if(ranktot==0) then
      !    write(0,*) 'PIQTB not available yet!'
      !    call fatal
      !  endif
      !endif  
      

      pi_isobaric=isobaric
      if(isobaric) then
        if(pi_start_isobaric>dt) isobaric=.FALSE.
        !masspiston=masspiston*nbeads
        !vextvol=vextvol*sqrt(real(nbeads,8))
        if(ranktot.eq.0 .and. barostat/='LANGEVIN') then
            write(0,*) "Error: PIMD only compatible"
     &         //" with langevin barostat yet"
            call fatal
        endif
      endif   

      only_long = centroid_longrange .and.
     &      (.not. (contract.and.nbeads_ctr==1) )
      if(only_long .and. polar_allbeads) then
        if(ranktot==0) then
            write(0,*) "Warning: centroid polarization contribution "//
     &    "to the virial might not be accurate !"//
     &    " Pressure could be WRONG !"
        endif
      endif       


c
c     print out a header line for the dynamics computation
c

      if(ranktot==0) then
        write(iout,*)
        write(iout,'(A)') 
     &     ' Path Integral MD trajectory via the '
     &      //trim(integrate)//' algorithm'
        write(iout,'(A,i3,A)', advance="no") 
     &     '  - with ',nbeads,' beads'
        if(contract) then
          write(iout,'(A,i3,A)', advance="no") 
     &     ' and ',nbeads_ctr,' contracted beads'
        endif
        write(iout,*)
        if (piqtb) then
          if (adaptive_qtb) then
            write(iout,'(A)')
     &     '  - with adaptive PI-QTB thermostat'
          else
            write(iout,'(A)')
     &     '  - with PI-QTB thermostat'
          endif
        endif
        if(cay_correction) then
          write(iout,'(A)') 
     &     '  - using CAY correction for fluctuation modes'
        endif
        if(centroid_longrange) then
           write(iout,'(A)') 
     &      '  - long-range interactions on centroid only'
        endif
        write(iout,*)
      endif

c
c     integrate equations of motion to take a time step
c
      time0=mpi_wtime()
      do istep = 1, nstep

        if(istep*dt>pi_start_isobaric) isobaric=pi_isobaric

        ! perform a step of PIMD integration
        if(integrate.eq.'BAOAB') then
          call baoabpi(polymer,polymer_ctr,istep,dt)
        elseif(integrate.eq.'BAOABRESPA') then
          call baoabrespapi(polymer,polymer_ctr,istep,dt)
        elseif(ranktot==0) then
          write(0,*) "Unknown integrator "//integrate//" for PIMD"
          call fatal
        endif

        if (mod(istep,iprint).eq.0) then
           time1 = mpi_wtime()
           timestep = time1-time0
           if (ranktot.eq.0) then
             call MPI_REDUCE(MPI_IN_PLACE,timestep,1,MPI_REAL8,MPI_SUM,
     $          0,MPI_COMM_WORLD,ierr)
           else
             call MPI_REDUCE(timestep,timestep,1,MPI_REAL8,MPI_SUM,0,
     $          MPI_COMM_WORLD,ierr)
           end if
           if (ranktot.eq.0) then
            write(6,1000) (timestep)/nproctot
     &         ,(timestep)/dble(nproctot*iprint)
            write(6,1010) 86400*dt*dble(iprint*nproctot)/(1000*timestep)
           end if
           time0 = time1
         end if
      end do

c     perform any final tasks before program exit
c
      call final

      end subroutine
        


      
