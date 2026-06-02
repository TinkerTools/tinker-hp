!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ###################################################################
!     ##                                                               ##
!     ##  program pimd  --  run pimd molecular or stochastic dynamics  ##
!     ##                                                               ##
!     ###################################################################
!
!
!     "pimd" computes a molecular dynamics trajectory
!     in one of the standard statistical mechanical ensembles and using
!     any of several possible integration methods
!
!
#include "tinker_macro.h"
program pimd
   use mpi
#ifdef _OPENACC
   use utilgpu,only: bind_gpu
#endif
   implicit none
   integer ierr,nthreadsupport
#ifdef _OPENACC
   call bind_gpu
#endif
   call MPI_INIT(ierr)
!      call MPI_INIT_THREAD(MPI_THREAD_MULTIPLE,nthreadsupport,ierr)
   call pimd_bis
   call MPI_BARRIER(MPI_COMM_WORLD,ierr)
   call MPI_FINALIZE(ierr)
end
!
subroutine pimd_bis
   use atomsMirror
   use bath
   use beads
   use bound
   use boxes, only: volbox
   use cutoff
   use domdec
   use keys
   use inform
   use iounit
   use mdstuf
   use moldyn
   use neigh
   use timestat
   use mpi
   use neigh, only: lbuffer
   use utilgpu ,only: rec_queue
   use tinMemory
   use tinheader
   use ascii, only: int_to_str
   use utils, only: set_to_zero2m
   use units
   use qtb, only: qtb_thermostat,adaptive_qtb,piqtb
   use commstuffpi
   use virial, only: use_virial
   implicit none
   integer i,istep,nstep,ierr,k,ibead,ibeadglob!,nstep_therm
   integer iglob
   integer mode,next,nseg
   real(r_p) dt,dtdump,time0,time1,timestep
   logical exist,query,restart,restart_bead
   character*20 keyword
   character*240 record
   character*240 string
   real(r_p) maxwell
   integer nbeads_para,j,iloc
   logical :: pi_isobaric,only_long
   real(t_p) :: bufpi

   app_id = pimd_a
   path_integral_md=.TRUE.
!
!
1000 Format(' Time for ',I4,' Steps: ',f15.4,/,&
      &' Ave. Time per step: ',f15.4)
1010 Format(' ns per day: ',f15.4)
!
!     set up the structure and molecular mechanics calculation
!
   call initial
   call init_keys
   call read_bead_config
   call initmpi
   call MPI_BARRIER(MPI_COMM_WORLD,ierr)
   call getxyz

   write(*,*) ranktot,'polymer',nproc_polymer,rank_polymer
   write(*,*) ranktot,'spatial',nproc,rank

   if(nproctot==1) pi_comm_report=.FALSE.
   if(pi_comm_report) then
      open(newunit=pi_comm_unit,&
         &file="comm_report_rank"//int_to_str(ranktot)//".out")
   endif

   call unitcell
   call cutoffs
   call lattice
!
!     setup for MPI
!
!
!     get nprocforce to do correct allocation of MPI arrays
!
   call drivermpi
   call reinitnl(0)
!
!     allocate some arrays
!
   call prmem_requestm(    v,3,n)
   call prmem_requestm(    a,3,n)
   call prmem_requestm( aalt,3,n)
   call prmem_requestm(aalt2,3,n)
!
   call set_to_zero2m(   a,    v,3*n,rec_queue)
   call set_to_zero2m(aalt,aalt2,3*n,rec_queue)
!
   call mechanic
   call nblist(0)
!
!     initialize the temperature, pressure and coupling baths
!
   kelvin = 0.0
   atmsph = 0.0
   isothermal = .false.
   isobaric = .false.
!
!     initialize the simulation length as number of time steps
!
   query = .true.
   call nextarg (string,exist)
   if (.not.exist) then
      if (ranktot.eq.0) write (iout,20)
20    format (/,' You Need To Enter the Number of Dynamics Steps')
      __TINKER_FATAL__
   else
      read (string,*,err=10,end=10)  nstep
      query = .false.
   end if
10 continue
!
!     get the length of the dynamics time step in picoseconds
!
   dt = -1.0
   call nextarg (string,exist)
   if (.not.exist) then
      if (ranktot.eq.0) write (iout,50)
50    format (/,' You Need To Enter the Time Step Length in Femtoseconds')
      __TINKER_FATAL__
   else
      read (string,*,err=40,end=40)  dt
   end if
40 continue
   dt = 0.001 * dt
!
!     enforce bounds on thermostat and barostat coupling times
!
   tautemp = max(tautemp,dt)
   taupres = max(taupres,dt)
!
!     set the time between trajectory snapshot coordinate dumps
!
   dtdump = -1.0
   call nextarg (string,exist)
   if (.not.exist) then
      if (rank.eq.0) write (iout,90)
90    format (/,' You Need To Enter Time between Dumps in Picoseconds')
   else
      read (string,*,err=80,end=80)  dtdump
   end if
80 continue
   iwrite = nint(dtdump/dt)
!
!     get choice of statistical ensemble for periodic system
!
!  if (use_bounds) then
   mode = -1
   call nextarg (string,exist)
   if (.not.exist) then
      if (ranktot.eQ.0) write (iout,130)
 130     format (/,' Available Statistical Mechanical Ensembles :'&
         ,//,4x,'(1) Microcanonical (NVE)'&
          ,/,4x,'(2) Canonical (NVT)'&
          ,/,4x,'(3) Isoenthalpic-Isobaric (NPH)'&
          ,/,4x,'(4) Isothermal-Isobaric (NPT)')
      __TINKER_FATAL__
   else
      read (string,*,err=120,end=120)  mode
   end if
 120  continue
   if (mode.eq.2 .or. mode.eq.4) then
      isothermal = .true.
      kelvin = -1.0
      call nextarg (string,exist)
      if (.not.exist) then
         write (iout,180)
 180        format (/,' You Need To Enter the Desired Temperature in Kelvin')
         __TINKER_FATAL__
      else
         read (string,*,err=170,end=170)  kelvin
      end if
 170  continue
      if (kelvin .le. 0.0)  kelvin = 298.0
   end if
   if (mode.eq.1 .or. mode.eq.3) then
      if (ranktot.eq.0) &
         write(iout,*) ' ERROR! NVE, NPH and NPT not implemented yet!'
      __TINKER_FATAL__
   endif
   if (mode.eq.4) then
      isobaric = .true.
      atmsph   = -1.0
      call nextarg (string,exist)
      if (.not.exist) then
         if (ranktot.eq.0) write (iout,220)
 220        format (/,' You Need To Enter the Desired Pressure in Atm')
      else
         read (string,*,err=210,end=210)  atmsph
      end if
 210     continue
      if (atmsph .le. 0.0)  atmsph = 1.0
   end if
!  end if
!
!     use constant energy or temperature for nonperiodic system
!
!   if (.not. use_bounds) then
!      mode = -1
!      call nextarg (string,exist)
!      if (exist)  read (string,*,err=250,end=250)  mode
!250   continue
!      do while (mode.lt.1 .or. mode.gt.2)
!         write (iout,260)
!260      format (/,' Available Simulation Control Modes :',&
!            &//,4x,'(1) Constant Total Energy Value (E)',&
!            &/,4x,'(2) Constant Temperature via Thermostat (T)',&
!            &//,' Enter the Number of the Desired Choice',&
!            &' [1] :  ',$)
!         read (input,270,err=280)  mode
!270      format (i10)
!         if (mode .le. 0)  mode = 1
!280      continue
!      end do
!      if (mode .eq. 2) then
!         isothermal = .true.
!         kelvin = -1.0
!         call nextarg (string,exist)
!         if (exist)  read (string,*,err=290,end=290)  kelvin
!290      continue
!         do while (kelvin .lt. 0.0)
!            write (iout,300)
!300         format (/,' Enter the Desired Temperature in Degrees',&
!               &' K [298] :  ',$)
!            read (input,310,err=320)  kelvin
!310         format (f20.0)
!            if (kelvin .le. 0.0)  kelvin = 298.0
!320         continue
!         end do
!      end if
!   end if
!
!     setup dynamics
!
   call reset_observables(.true.)
   call allocpi()
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
   call mdinitbead(dt,polymer_main)

   pi_isobaric=isobaric
   if(isobaric) then
      if(pi_start_isobaric>dt) isobaric=.FALSE.
      !masspiston=masspiston*nbeads
      !vextvol=vextvol*sqrt(real(nbeads,r_p))
      if(ranktot.eq.0 .and. barostat/='LANGEVIN') then
         write(0,*) "Error: PIMD only compatible"&
            &//" with langevin barostat yet"
         call fatal
      endif
   endif

   only_long = centroid_longrange .and.&
      &(.not. (contract.and.nbeads_ctr==1) )
   if(only_long .and. polar_allbeads) then
      if(ranktot==0 .and. use_virial) then
         write(0,*) "Warning: centroid polarization contribution "//&
            &"to the virial might not be accurate !"//&
            &" Pressure could be WRONG !"
      endif
   endif

!
!     print out a header line for the dynamics computation
!

   if(ranktot==0) then
      write(iout,*)
      write(iout,'(A)')&
         &' Path Integral MD trajectory via the '&
         &//trim(integrate)//' algorithm'
      write(iout,'(A,i3,A)', advance="no")&
         &'  - with ',nbeads,' beads'
      if(contract) then
         write(iout,'(A,i3,A)', advance="no")&
            &' and ',nbeads_ctr,' contracted beads'
      endif
      write(iout,*)
      if (piqtb) then
         if (adaptive_qtb) then
            write(iout,'(A)')&
               &'  - with adaptive PI-QTB thermostat'
         else
            write(iout,'(A)')&
               &'  - with PI-QTB thermostat'
         endif
      endif
      if(cay_correction) then
         write(iout,'(A)')&
            &'  - using CAY correction for fluctuation modes'
      endif
      if(centroid_longrange) then
         write(iout,'(A)')&
            &'  - long-range interactions on centroid only'
         if(.not. polar_allbeads) then
            write(iout,'(A)') '  - polarization on centroid only'
         endif
      endif
      write(iout,*)
   endif

!
!     integrate equations of motion to take a time step
!

   time0 = 0 !mpi_wtime()
   do istep = 1, nstep
      step_c = istep
      call timer_enter( timer_timestep )

      if(istep*dt>pi_start_isobaric) isobaric=pi_isobaric

      ! perform a step of PIMD integration
      if(integrate.eq.'BAOAB') then
         call baoabpi(istep,dt)
      elseif(integrate.eq.'BAOABRESPA') then
         call baoabrespapi(istep,dt)
      elseif(ranktot==0) then
         write(0,*) "Unknown integrator "//integrate//" for PIMD"
         call fatal
      endif

      call timer_exit( timer_timestep )

      if (mod(istep,iprint).eq.0) then
         time1 = timer_get_total( timer_timestep )
         timestep = time1-time0
         time0 = time1
         if (ranktot.eq.0) then
            call MPI_REDUCE(MPI_IN_PLACE,timestep,1,MPI_REAL8,MPI_SUM,&
               &0,MPI_COMM_WORLD,ierr)
         else
            call MPI_REDUCE(timestep,timestep,1,MPI_REAL8,MPI_SUM,0,&
               &MPI_COMM_WORLD,ierr)
         end if
         if (verbose.and.ranktot.eq.0) then
            write(6,1000) iprint, (timestep)/nproctot,(timestep)/&
               &dble(nproctot*iprint)
            write(6,1010) 86400*dt*real(iprint*nproctot,t_p)/&
               &(1000*timestep)
         end if
      end if

      ! Abort if any problem detected
      if (abort) call fatal
   end do

!     perform any final tasks before program exit
!
   call final

end subroutine




