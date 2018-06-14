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
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'bound.i'
      include 'boxes.i'
      include 'cutoff.i'
      include 'inform.i'
      include 'inter.i'
      include 'iounit.i'
      include 'mdstuf.i'
      include 'molcul.i'
      include 'units.i'
      include 'openmp.i'
      include 'potent.i'
      include 'timestat.i'
      include 'mpif.h'
      integer istep,modstep
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
c
c     print energy, temperature and pressure for current step
c
      if (rank.eq.0) then
        if (verbose) then
           if (modstep .eq. 1) then
              if (use_bounds .and. integrate.ne.'STOCHASTIC') then
                 write (iout,10)
   10            format (/,4x,'MD Step',6x,'E Total',3x,'E Potential',
     &                      5x,'E Kinetic',7x,'Temp',7x,'Pres',/)
              else
                 write (iout,20)
   20            format (/,4x,'MD Step',6x,'E Total',3x,'E Potential',
     &                      5x,'E Kinetic',7x,'Temp',/)
              end if
           end if
           if (use_bounds .and. integrate.ne.'STOCHASTIC') then
              write (iout,30)  istep,etot,epot,ekin,temp,pres
   30         format (i10,3f14.4,2f11.2)
           else
              write (iout,40)  istep,etot,epot,ekin,temp
   40         format (i10,3f14.4,f11.2)
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
cc
cc       note deformation value for potential energy smoothing
cc
c        if (use_smooth) then
c           if (modstep .eq. 0) then
c              write (iout,140)  deform
c  140         format (' Deformation',9x,f15.3,' Sqr Angs')
c           end if
c        end if
      end if
c
c     display the different times of the computation
c
      if ((rank.eq.0).and.(verbose)) then
        if (modstep.eq.0) then
          write(iout,150) timeclear/dble(iprint)
 150      format('Ave time for aclear comms            =  ',F24.15)
          timeclear = 0.0d0
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
        end if
      end if
      if (((use_pmecore).and.(rank.eq.ndir)).or.
     $ (.not.(use_pmecore).and.(rank.eq.0))) then
        if ((modstep.eq.0).and.(verbose)) then
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
      return
      end
