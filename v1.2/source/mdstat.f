c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine mdstat  --  compute averages over a trajectory  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "mdstat" is called at each molecular dynamics time step to
c     form statistics on various average values and fluctuations,
c     and to periodically save the state of the trajectory
c
c
      subroutine mdstat (istep,dt,etot,epot,ekin,temp,pres)
      use sizes
      use atoms
      use bound
      use boxes
      use cutoff
      use domdec
      use inform
      use inter
      use iounit
      use mdstuf
      use molcul
      use potent
      use timestat
      use units
      use mpi
      implicit none
      integer istep,modstep,ierr,i
      real*8 dt,temp,pres
      real*8 etot,epot,ekin
      real*8 pico,dens
      real*8 fluctuate,fluctuate2
      real*8 intfluct,intfluct2
      real*8 potfluct,potfluct2
      real*8 kinfluct,kinfluct2
      real*8 tfluct,pfluct,dfluct
      real*8 tfluct2,pfluct2,dfluct2
      real*8 etot_sum,etot2_sum
      real*8 eint_sum,eint2_sum
      real*8 etot_ave,etot2_ave
      real*8 eint_ave,eint2_ave
      real*8 epot_sum,epot2_sum
      real*8 ekin_sum,ekin2_sum
      real*8 epot_ave,epot2_ave
      real*8 ekin_ave,ekin2_ave
      real*8 temp_sum,temp2_sum
      real*8 temp_ave,temp2_ave
      real*8 pres_sum,pres2_sum
      real*8 pres_ave,pres2_ave
      real*8 dens_sum,dens2_sum
      real*8 dens_ave,dens2_ave
      save etot_sum,etot2_sum
      save eint_sum,eint2_sum
      save epot_sum,epot2_sum
      save ekin_sum,ekin2_sum
      save temp_sum,temp2_sum
      save pres_sum,pres2_sum
      save dens_sum,dens2_sum
c
c
c     set number of steps for block averages of properties
c
      modstep = mod(istep,iprint)
c
c     zero out summation variables for new averaging period
c
      if (modstep.eq.1 .or. iprint.eq.1) then
         etot_sum = 0.0d0
         etot2_sum = 0.0d0
         epot_sum = 0.0d0
         epot2_sum = 0.0d0
         ekin_sum = 0.0d0
         ekin2_sum = 0.0d0
         eint_sum = 0.0d0
         eint2_sum = 0.0d0
         temp_sum = 0.0d0
         temp2_sum = 0.0d0
         pres_sum = 0.0d0
         pres2_sum = 0.0d0
         dens_sum = 0.0d0
         dens2_sum = 0.0d0
      end if
      if (rank.eq.0) then
        if (verbose) then
           if (modstep .eq. 1) then
               write (iout,10)
   10          format (/,4x,'MD Step',6x,'E Total',3x,'E Potential',
     &                    5x,'E Kinetic',7x,'Temp',7x,'Pres',/)
           end if
           write (iout,30)  istep,etot,epot,ekin,temp,pres
   30      format (i10,3f14.4,2f11.2)
        end if
c
c       print header for the averages over a group of recent steps
c
        if (modstep .eq. 0) then
           pico = dble(istep) * dt
           write (iout,50)  iprint,istep
   50      format (/,' Average Values for the last',i6,' out of',
     &                i9,' Dynamics Steps')
           write (iout,60)  pico
   60      format (/,' Simulation Time',5x,f15.4,' Picosecond')
        end if
c
c       compute total energy and fluctuation for recent steps
c
        etot_sum = etot_sum + etot
        etot2_sum = etot2_sum + etot**2
        if (modstep .eq. 0) then
           etot_ave = etot_sum / dble(iprint)
           etot2_ave = etot2_sum / dble(iprint)
           fluctuate2 = etot2_ave - etot_ave**2
           if (fluctuate2 .gt. 0.0d0) then
              fluctuate = sqrt(fluctuate2)
           else
              fluctuate = 0.0d0
           end if
           write (iout,70)  etot_ave,fluctuate
   70      format (' Total Energy',8x,f15.4,' Kcal/mole',3x,
     &                '(+/-',f9.4,')')
        end if
c
c       compute average potential energy and its fluctuation
c
        epot_sum = epot_sum + epot
        epot2_sum = epot2_sum + epot**2
        if (modstep .eq. 0) then
           epot_ave = epot_sum / dble(iprint)
           epot2_ave = epot2_sum / dble(iprint)
           potfluct2 = epot2_ave - epot_ave**2
           if (potfluct2 .gt. 0.0d0) then
              potfluct = sqrt(potfluct2)
           else
              potfluct = 0.0d0
           end if
           write (iout,80)  epot_ave,potfluct
   80      format (' Potential Energy',4x,f15.4,' Kcal/mole',3x,
     &                '(+/-',f9.4,')')
        end if
c
c       compute average kinetic energy and its fluctuation
c
        ekin_sum = ekin_sum + ekin
        ekin2_sum = ekin2_sum + ekin**2
        if (modstep .eq. 0) then
           ekin_ave = ekin_sum / dble(iprint)
           ekin2_ave = ekin2_sum / dble(iprint)
           kinfluct2 = ekin2_ave - ekin_ave**2
           if (kinfluct2 .gt. 0.0d0) then
              kinfluct = sqrt(kinfluct2)
           else
              kinfluct = 0.0d0
           end if
           write (iout,90)  ekin_ave,kinfluct
   90      format (' Kinetic Energy',6x,f15.4,' Kcal/mole',3x,
     &                '(+/-',f9.4,')')
        end if
c
c       compute average intermolecular energy and its fluctuation
c
        if (nmol.ne.1 .and. nmol.ne.n .and. .not.use_ewald) then
           eint_sum = eint_sum + einter
           eint2_sum = eint2_sum + einter**2
           if (modstep .eq. 0) then
              eint_ave = eint_sum / dble(iprint)
              eint2_ave = eint2_sum / dble(iprint)
              intfluct2 = eint2_ave - eint_ave**2
              if (intfluct2 .gt. 0.0d0) then
                 intfluct = sqrt(intfluct2)
              else
                 intfluct = 0.0d0
              end if
              write (iout,100)  eint_ave,intfluct
  100         format (' Intermolecular',6x,f15.4,' Kcal/mole',3x,
     &                   '(+/-',f9.4,')')
           end if
        end if
c
c       compute the average temperature and its fluctuation
c
        temp_sum = temp_sum + temp
        temp2_sum = temp2_sum + temp**2
        if (modstep .eq. 0) then
           temp_ave = temp_sum / dble(iprint)
           temp2_ave = temp2_sum / dble(iprint)
           tfluct2 = temp2_ave - temp_ave**2
           if (tfluct2 .gt. 0.0d0) then
              tfluct = sqrt(tfluct2)
           else
              tfluct = 0.0d0
           end if
           write (iout,110)  temp_ave,tfluct
  110      format (' Temperature',9x,f15.2,' Kelvin',6x,
     &                '(+/-',f9.2,')')
        end if
c
c       compute the average pressure and its fluctuation
c
        if (use_bounds) then
           pres_sum = pres_sum + pres
           pres2_sum = pres2_sum + pres**2
           if (modstep .eq. 0) then
              pres_ave = pres_sum / dble(iprint)
              pres2_ave = pres2_sum / dble(iprint)
              pfluct2 = pres2_ave - pres_ave**2
              if (pfluct2 .gt. 0.0d0) then
                 pfluct = sqrt(pfluct2)
              else
                 pfluct = 0.0d0
              end if
              write (iout,120)  pres_ave,pfluct
  120         format (' Pressure',12x,f15.2,' Atmosphere',2x,
     &                   '(+/-',f9.2,')')
           end if
c
c       compute the average density and its fluctuation
c
           dens = (1.0d24/volbox) * (totmass/avogadro)
           dens_sum = dens_sum + dens
           dens2_sum = dens2_sum + dens**2
           if (modstep .eq. 0) then
              dens_ave = dens_sum / dble(iprint)
              dens2_ave = dens2_sum / dble(iprint)
              dfluct2 = dens2_ave - dens_ave**2
              if (dfluct2 .gt. 0.0d0) then
                 dfluct = sqrt(dfluct2)
              else
                 dfluct = 0.0d0
              end if
              write (iout,130)  dens_ave,dfluct
  130         format (' Density',13x,f15.4,' Grams/cc',4x,
     &                   '(+/-',f9.4,')')
           end if
        end if
      end if
      return
      end
c
c     subroutine mdstattime: compute averages of various timings and
c      output them if necessary
c
      subroutine mdstattime(istep,dt)
      use domdec
      use inform
      use iounit
      use timestat
      use mpi
      use potent
      implicit none
      integer istep,modstep,ierr,i
      real*8 dt
      real*8 buffer(24)
c
c     set number of steps for block averages of properties
c
      modstep = mod(istep,iprint)
c
c     sum the timing values to get average over the cores
c
      if (modstep.eq.0) then
        buffer(1) = timeinte
        buffer(2) = timereneig
        buffer(3) = timecommpos
        buffer(4) = timeparam
        buffer(5) = timenl
        buffer(6) = timegrad
        buffer(7) = timered
        buffer(8) = timetp
        buffer(9) = timecommforces
        buffer(10) = timestep
        buffer(11) = timebonded
        buffer(12) = timevdw
        buffer(13) = timeelec
        buffer(14) = timepolar
        buffer(15) = timecleargrad
        buffer(16) = timereal
        buffer(17) = timerec
        buffer(18) = timecommforcesrec
        buffer(19) = timecommforcesreal
        buffer(20) = timegrid
        buffer(21) = timefft
        buffer(22) = timescalar
        buffer(23) = timegrid2
        buffer(24) = timerecreccomm
c
        if (rank.eq.0) then
          call MPI_REDUCE(MPI_IN_PLACE,buffer,24,MPI_REAL8,MPI_SUM,0,
     $       COMM_TINKER,ierr)
        else
          call MPI_REDUCE(buffer,buffer,24,MPI_REAL8,MPI_SUM,0,
     $       COMM_TINKER,ierr)
        end if
c
        if (rank.eq.0) then
          do i = 1, 24
            buffer(i) = buffer(i)/nproc
          end do
          timeinte = buffer(1)
          timereneig = buffer(2)
          timecommpos = buffer(3)
          timeparam = buffer(4)
          timenl = buffer(5)
          timegrad = buffer(6)
          timered = buffer(7)
          timetp = buffer(8)
          timecommforces = buffer(9)
          timestep = buffer(10)
          timebonded = buffer(11)
          timevdw = buffer(12)
          timeelec = buffer(13)
          timepolar = buffer(14)
          timecleargrad = buffer(15)
          timereal = buffer(16)
          timerec = buffer(17)
          timecommforcesrec = buffer(18)
          timecommforcesreal = buffer(19)
          timegrid = buffer(20)
          timefft = buffer(21)
          timescalar = buffer(22)
          timegrid2 = buffer(23)
          timerecreccomm = buffer(24)
        end if
      end if
c
c     display the different times of the computation
c
      if ((rank.eq.0).and.(verbose)) then
        if (modstep.eq.0) then
          write(iout,*) '=================================='
          write(iout,160) timeinte/dble(iprint),timeinte*100/timestep
 160      format('Ave time for pos./vel. update     =  ',F14.6,
     $      ' seconds',2x,F14.3,'%')
          write(iout,170) timereneig/dble(iprint),
     $      timereneig*100/timestep
 170      format('Ave time for reneigh              =  ',F14.6,
     $      ' seconds',2x,F14.3,'%')
          write(iout,180) timecommpos/dble(iprint),
     $     timecommpos*100/timestep
 180      format('Ave time for pos. comms.          =  ',F14.6,
     $      ' seconds',2x,F14.3,'%')
          write(iout,190) timeparam/dble(iprint),
     $     timeparam*100/timestep
 190      format('Ave time for param regeneration   =  ',F14.6,
     $      ' seconds',2x,F14.3,'%')
          write(iout,200) timegrad/dble(iprint),
     $     timegrad*100/timestep
 200      format('Ave time in gradient routines     =  ',F14.6,
     $      ' seconds',2x,F14.3,'%')
          write(iout,210) timenl/dble(iprint),
     $     timenl*100/timestep
 210      format('Ave time for neighbor lists       =  ',F14.6,
     $      ' seconds',2x,F14.3,'%')
          write(iout,220) timered/dble(iprint),
     $     timered*100/timestep
 220      format('Ave time for energy reduction     =  ',F14.6,
     $      ' seconds',2x,F14.3,'%')
          write(iout,230) timetp/dble(iprint),
     $     timetp*100/timestep
 230      format('Ave time for Temp/Press control   =  ',F14.6,
     $      ' seconds',2x,F14.3,'%')
          write(iout,240) timecommforces/dble(iprint),
     $     timecommforces*100/timestep
 240      format('Ave time for forces comms.        =  ',F14.6,
     $      ' seconds',2x,F14.3,'%')
          write(iout,250) timestep/dble(iprint),
     $    timestep*100/timestep
 250      format('Ave time per step                 =  ',F14.6,
     $      ' seconds',2x,F14.3,'%')
          write(iout,260) 86400*dt*dble(iprint)/(1000*timestep)
 260      format(' ns per day: ',f15.4)
          write(iout,*) '=================================='
        end if
      end if

      if ((rank.eq.0).and.(verboseforcestime)) then
        if (modstep.eq.0) then
          write(iout,*) '=================================='
          write(iout,270) timebonded/dble(iprint),
     $     timebonded*100/timegrad
 270      format('Ave time for bonded forces        =  ',F14.6,
     $      ' seconds',2x,F14.3,'%')
          write(iout,280) timevdw/dble(iprint),
     $     timevdw*100/timegrad
 280      format('Ave time for vdw forces           =  ',F14.6,
     $      ' seconds',2x,F14.3,'%')
          write(iout,290) timeelec/dble(iprint),
     $     timeelec*100/timegrad
 290      format('Ave time for electrostatic forces =  ',F14.6,
     $      ' seconds',2x,F14.3,'%')
          write(iout,300) timepolar/dble(iprint),
     $     timepolar*100/timegrad
 300      format('Ave time for polarization forces  =  ',F14.6,
     $      ' seconds',2x,F14.3,'%')
          write(iout,310) timecleargrad/dble(iprint),
     $     timecleargrad*100/timegrad
 310      format('Ave time for zero./sum. forces    =  ',F14.6,
     $      ' seconds',2x,F14.3,'%')
          write(iout,*) '=================================='
        end if
      end if

      if ((rank.eq.0).and.(verbosepmetime)) then
        if (modstep.eq.0) then
          write(iout,*) '=================================='
          write(iout,320) timereal/dble(iprint)
 320      format('Ave time for real space forces (permanent)        =  '
     $      ,F14.6,' seconds')
          write(iout,330) timerec/dble(iprint)
 330      format('Ave time for rec. space forces (permanent)        =  '
     $      ,F14.6,' seconds')
          write(iout,340) timecommforcesreal/dble(iprint)
 340      format('Ave time for real space forces comms. (permanent) =  '
     $      ,F14.6,' seconds')
          write(iout,350) timecommforcesrec/dble(iprint)
 350      format('Ave time for rec. space forces comms. (permanent) =  '
     $      ,F14.6,' seconds')
          write(iout,360) timegrid/dble(iprint)
 360      format('Ave time fill pme grid (permanent)                =  '
     $      ,F14.6,' seconds')
          write(iout,370) timefft/dble(iprint)
 370      format('Ave time to do ffts (permanent)                   =  '
     $      ,F14.6,' seconds')
          write(iout,380) timescalar/dble(iprint)
 380      format('Ave time to do pme scalar prod. (permanent)       =  '
     $      ,F14.6,' seconds')
          write(iout,390) timegrid2/dble(iprint)
 390      format('Ave time to get rec. potential (permanent)        =  '
     $      ,F14.6,' seconds')
          write(iout,400) timerecreccomm/dble(iprint)
 400      format('Ave time for rec. rec. comms. (permanent)         =  '
     $      ,F14.6,' seconds')
          write(iout,*) '=================================='
        end if
      end if
c
      if (modstep.eq.0) then
c
c     put timing variables to zero
c
        timeinte  = 0d0
        timereneig = 0.0d0
        timecommpos = 0.0d0
        timeparam = 0.0d0
        timenl = 0d0
        timegrad = 0.0d0
        timered = 0.0d0
        timetp = 0.0d0
        timecommforces = 0.0d0
        timestep = 0d0
        timebonded = 0d0
        timevdw = 0d0
        timeelec = 0d0
        timepolar = 0d0
        timecleargrad = 0d0
        timereal = 0d0
        timerec = 0d0
        timecommforcesrec = 0d0
        timecommforcesreal = 0d0
        timegrid = 0d0
        timefft = 0d0
        timescalar = 0d0
        timegrid2 = 0d0
        timerecreccomm = 0d0
      end if
      return
      end

