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
      use adqtb
      use sizes
      use atoms
      use bath
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
      use stat
      use timestat
      use units
      use mpi
      implicit none
      integer istep,modstep,ierr,i
      real*8 dt,temp,pres
      real*8 etot,epot,ekin,vol
      real*8 pico,dens,pistontemp
      real*8 fluctuate,fluctuate2
      real*8 intfluct,intfluct2
      real*8 potfluct,potfluct2
      real*8 kinfluct,kinfluct2
      real*8 tfluct,pfluct,dfluct,volfluct
      real*8 tfluct2,pfluct2,dfluct2,volfluct2
c      real*8 etot_sum,etot2_sum
c      real*8 eint_sum,eint2_sum
      real*8 etot_ave,etot2_ave
      real*8 eint_ave,eint2_ave
c      real*8 epot_sum,epot2_sum
c      real*8 ekin_sum,ekin2_sum
      real*8 epot_ave,epot2_ave
      real*8 ekin_ave,ekin2_ave
c      real*8 temp_sum,temp2_sum
      real*8 temp_ave,temp2_ave
c      real*8 pres_sum,pres2_sum
      real*8 pres_ave,pres2_ave
c      real*8 dens_sum,dens2_sum
      real*8 dens_ave,dens2_ave
      real*8 vol_ave,vol2_ave
      real*8 pistontemp_ave,pistontemp2_ave
      real*8 buffer(18)
c      save etot_sum,etot2_sum
c      save eint_sum,eint2_sum
c      save epot_sum,epot2_sum
c      save ekin_sum,ekin2_sum
c      save temp_sum,temp2_sum
c      save pres_sum,pres2_sum
c      save dens_sum,dens2_sum
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
         vol_sum = 0.0d0
         vol2_sum = 0.0d0
         pistontemp_sum = 0.0d0
         pistontemp2_sum = 0.0d0
      end if
c
c     print energy, temperature and pressure for current step
c
      if (modstep.eq.0) then
c
c     sum the timing values to get average over the cores
c
        buffer(1) = timereneig
        buffer(2) = timecommstep
        buffer(3) = timeparam
        buffer(4) = timeforcescomm
        buffer(5) = timedirreccomm
        buffer(6) = timebondedvdw
        buffer(7) = timebonded
        buffer(8) = timenonbonded
        buffer(9) = timereal
        buffer(10) = timerealdip
        buffer(11) = timegrid1
        buffer(12) = timeffts
        buffer(13) = timescalar
        buffer(14) = timegrid2
        buffer(15) = timerecreccomm
        buffer(16) = timerec
        buffer(17) = timerecdip
c
        if (rank.eq.0) then
          call MPI_REDUCE(MPI_IN_PLACE,buffer,17,MPI_REAL8,MPI_SUM,0,
     $       COMM_BEAD,ierr)
        else
          call MPI_REDUCE(buffer,buffer,17,MPI_REAL8,MPI_SUM,0,
     $       COMM_BEAD,ierr)
        end if
c
        if (use_pmecore) then
          if (rank.eq.0) then
            do i = 1, 5
              buffer(i) = buffer(i)/ndir
            end do
            do i = 6, 8
              buffer(i) = buffer(i)/nproc
            end do
            do i = 9, 10
              buffer(i) = buffer(i)/ndir
            end do
            do i = 11, 17
              buffer(i) = buffer(i)/nrec
            end do
            timereneig = buffer(1)
            timecommstep = buffer(2)
            timeparam = buffer(3)
            timeforcescomm = buffer(4)
            timedirreccomm = buffer(5)
            timebondedvdw = buffer(6)
            timebonded = buffer(7)
            timenonbonded = buffer(8)
            timereal = buffer(9)
            timerealdip = buffer(10)
            timegrid1 = buffer(11)
            timeffts = buffer(12)
            timescalar = buffer(13) 
            timegrid2 = buffer(14) 
            timerecreccomm = buffer(15)
            timerec = buffer(16)
            timerecdip = buffer(17)
          end if
        else
          if (rank.eq.0) then
            do i = 1, 17
              buffer(i) = buffer(i)/nproc
            end do
            timereneig = buffer(1)
            timecommstep = buffer(2)
            timeparam = buffer(3)
            timeforcescomm = buffer(4)
            timedirreccomm = buffer(5)
            timebondedvdw = buffer(6)
            timebonded = buffer(7)
            timenonbonded = buffer(8)
            timereal = buffer(9)
            timerealdip = buffer(10)
            timegrid1 = buffer(11)
            timeffts = buffer(12)
            timescalar = buffer(13) 
            timegrid2 = buffer(14) 
            timerecreccomm = buffer(15)
            timerec = buffer(16)
            timerecdip = buffer(17)
          end if
        end if
      end if
      if (rank.eq.0) then
        if (verbose) then
           if (modstep .eq. 1) then
              if (use_bounds) then
                 write (iout,10)
   10            format (/,4x,'MD Step',6x,'E Total',3x,'E Potential',
     &                      5x,'E Kinetic',7x,'Temp',7x,'Pres',
     &                      5x,'Volume',/)
              else
                 write (iout,20)
   20            format (/,4x,'MD Step',6x,'E Total',3x,'E Potential',
     &                      5x,'E Kinetic',7x,'Temp',5x,'Volume'/)
              end if
           end if
           if (use_bounds .and. integrate.ne.'STOCHASTIC') then
              write (iout,30)  istep,etot,epot,ekin/corr_fact_qtb,temp
     &                          ,pres,volbox
   30         format (i10,3f14.4,3f11.2)
           else
              write (iout,40)  istep,etot,epot,ekin/corr_fact_qtb,temp
     &                          ,volbox
   40         format (i10,3f14.4,2f11.2)
           end if
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
        open(44,file='Total_energie.out')
        write(44,*) etot
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

c       compute the average volume and its fluctuation
c
          open(45,file='Volume.out')
          write(45,*) volbox
           vol_sum = vol_sum + volbox
           vol2_sum = vol2_sum + volbox**2
           if (modstep .eq. 0) then
              vol_ave = vol_sum / dble(iprint)
              vol2_ave = vol2_sum / dble(iprint)
              volfluct2 = vol2_ave - vol_ave**2
              if (volfluct2 .gt. 0.0d0) then
                 volfluct = sqrt(volfluct2)
              else
                 volfluct = 0.0d0
              end if
              write (iout,150)  vol_ave,volfluct
  150         format (' Volume',13x,f15.4,' Angstrom^3',3x,
     &                   '(+/-',f9.4,')')
           end if

c       compute the average piston temperature and its fluctuation
c
           if ((barostat.eq.'LANGEVIN').and.(rank.eq.0)) then
             pistontemp = (masspiston*vextvol*vextvol)/boltzmann
             pistontemp_sum = pistontemp_sum + pistontemp
             pistontemp2_sum = pistontemp2_sum + pistontemp**2
             if (modstep .eq. 0) then
                pistontemp_ave = pistontemp_sum / dble(iprint)
                pistontemp2_ave = pistontemp2_sum / dble(iprint)
                dfluct2 = pistontemp2_ave - pistontemp_ave**2
                if (dfluct2 .gt. 0.0d0) then
                   dfluct = sqrt(dfluct2)
                else
                   dfluct = 0.0d0
                end if
                write (iout,131)  pistontemp_ave,dfluct
  131           format (' Piston temperature',5x,f15.4,' K',4x,
     &                     '(+/-',f9.4,')')
             end if
           end if
        end if
      end if
c
c
c     display the different times of the computation
c
      if ((rank.eq.0).and.(verbose)) then
        if (modstep.eq.0) then
c          write(iout,150) timeclear/dble(iprint)
c 150      format('Ave time for aclear comms            =  ',F24.15)
c          timeclear = 0.0d0
          write(iout,160) timereneig/dble(iprint)
 160      format('Ave time for reneig                  =  ',F24.15)
          timereneig = 0.0d0
          write(iout,170) timecommstep/dble(iprint)
 170      format('Ave time for positions comm          =  ',F24.15)
          timecommstep = 0.0d0
          write(iout,180) timeparam/dble(iprint)
 180      format('Ave time for param                   =  ',F24.15)
          timeparam = 0.0d0
          write(iout,190) timeforcescomm/dble(iprint)
 190      format('Ave time for forces comm             =  ',F24.15)
          timeforcescomm = 0.0d0
          write(iout,210) timedirreccomm/dble(iprint)
 210      format('Ave time for reciprocal forces comm  =  ',F24.15)
          timedirreccomm = 0.0d0
          write(iout,220) timebondedvdw/dble(iprint)
 220      format('Ave time for real space forces comm  =  ',F24.15)
          timebondedvdw = 0.0d0
          write(iout,230) timebonded/dble(iprint)
 230      format('Ave time for bonded forces           =  ',F24.15)
          timebonded = 0.0d0
          write(iout,240) timenonbonded/dble(iprint)
 240      format('Ave time for non bonded forces       =  ',F24.15)
          timenonbonded = 0.0d0
          write(iout,250) timenl/dble(iprint)
 250      format('Ave time for neighbor lists          =  ',F24.15)
          timenl = 0.0d0
          write(iout,260) timereal/dble(iprint)
 260      format('Ave time for real space (permanent)  =  ',F24.15)
          timereal = 0.0d0
          write(iout,270) timerealdip/dble(iprint)
 270      format('Ave time for real space (polar)      =  ',F24.15)
          timerealdip = 0.0d0
          write(iout,280) timegrid1/dble(iprint)
 280      format('Ave time for fill grid (permanent)   =  ',F24.15)
          timegrid1 = 0.0d0
          write(iout,290) timeffts/dble(iprint)
 290      format('Ave time for ffts (permanent)        =  ',F24.15)
          timeffts = 0.0d0
          write(iout,300) timescalar/dble(iprint)
 300      format('Ave time for scalar prod (permanent) =  ',F24.15)
          timescalar = 0.0d0
          write(iout,310) timegrid2/dble(iprint)
 310      format('Ave time for extract grid (permanent)=  ',F24.15)
          timegrid2 = 0.0d0
          write(iout,320) timerecreccomm/dble(iprint)
 320      format('Ave time for rec-rec comm (permanent)=  ',F24.15)
          timerecreccomm = 0.0d0
          write(iout,330) timerec/dble(iprint)
 330      format('Ave time for recip space (permanent) =  ',F24.15)
          timerec = 0.0d0
          write(iout,340) timerecdip/dble(iprint)
 340      format('Ave time for recip space (polar)     =  ',F24.15)
          timerecdip = 0.0d0
        end if
      end if
c
      if (modstep.eq.0) then
c
c     put timing variables to zero
c
        timeclear = 0.0d0
        timereneig = 0.0d0
        timecommstep = 0.0d0
        timeparam = 0.0d0
        timeforcescomm = 0.0d0
        timedirreccomm = 0.0d0
        timebondedvdw = 0.0d0
        timebonded = 0.0d0
        timenonbonded = 0.0d0
        timereal = 0.0d0
        timerealdip = 0.0d0
        timegrid1 = 0.0d0
        timeffts = 0.0d0
        timescalar = 0.0d0
        timegrid2 = 0.0d0
        timerecreccomm = 0.0d0
        timerec = 0.0d0
        timerecdip = 0.0d0
      end if
      return
      end
