c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine mdstatpi  --  compute averages over a pimd      ##
c     ##                                           trajectory        ##
c     #################################################################
c
c
c     "mdstatpi" is called at each molecular dynamics time step to
c     form statistics on various average values and fluctuations,
c     and to periodically save the state of the trajectory
c
c
      subroutine mdstatpi (istep,dt,ekvir,ekprim,presvir)
      use sizes
      use atoms
      use bath
      use beads
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
      integer mode
      real*8 dt,temp,pres
      real*8 ekvir,ekprim,presvir
      real*8 etotppi,epot,ekin,vol
      real*8 pico,dens,pistontemppi
      real*8 fluctuatepi,fluctuate2pi
      real*8 intfluctpi,intfluct2pi
      real*8 potfluctpi,potfluct2pi
      real*8 kinfluctpi,kinfluct2pi
      real*8 tfluctpi,pfluctpi,dfluctpi,volfluctpi
      real*8 tfluct2pi,pfluct2pi,dfluct2pi,volfluct2pi
      real*8 etotpi_ave,etot2pi_ave
c      real*8 eint_sum,eint2_sum
      real*8 eintpi_ave,eint2pi_ave
c      real*8 epot_sum,epot2_sum
c      real*8 ekin_sum,ekin2_sum
      real*8 epotppi_ave,epot2pi_ave
      real*8 ekinpi_ave,ekin2pi_ave
c      real*8 temp_sum,temp2_sum
      real*8 tempppi_ave,temp2pi_ave
c      real*8 pres_sum,pres2_sum
      real*8 prespi_ave,pres2pi_ave
c      real*8 dens_sum,dens2_sum
      real*8 denspi_ave,dens2pi_ave
      real*8 volpi_ave,vol2pi_ave
      real*8 pistontemppi_ave,pistontemp2pi_ave
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
         etotpi_sum = 0.0d0
         etot2pi_sum = 0.0d0
         epotpi_sum = 0.0d0
         epot2pi_sum = 0.0d0
         ekinpi_sum = 0.0d0
         ekin2pi_sum = 0.0d0
         eintpi_sum = 0.0d0
         eint2pi_sum = 0.0d0
         temppi_sum = 0.0d0
         temp2pi_sum = 0.0d0
         prespi_sum = 0.0d0
         pres2pi_sum = 0.0d0
         denspi_sum = 0.0d0
         dens2pi_sum = 0.0d0
         volpi_sum = 0.0d0
         vol2pi_sum = 0.0d0
         pistontemppi_sum = 0.0d0
         pistontemp2pi_sum = 0.0d0
      end if

        dens = (1.0d24/extvol) * (totmass/avogadro)
        if (verbose) then
          if(isobaric) then
            if (modstep .eq. 1) then
                write (iout,110)
110            format (/,4x,'MD Step',6x,'E Total',3x,'E Potential'
     &          ,7x,'Ek vir',7x,'Ek prim',7x,'Temp',5x,'Temp_cl'
     &          ,5x,'Pres',7x,'Density',7x,'Volume'/)
            end if
            write (iout,140)  istep,ekvir+epotpi_loc,epotpi_loc,ekvir
     &             ,ekprim,temppi,temppi_cl,presvir,dens,extvol
140         format (i10,4f14.4,3f11.2,f10.4,f14.4)
          else
            if (modstep .eq. 1) then
                write (iout,160)
160            format (/,4x,'MD Step',6x,'E Total',3x,'E Potential'
     &          ,7x,'Ek vir',7x,'Ek prim',7x,'Temp',5x,'Temp_cl'
     &          ,5x,'Pres'/)
            end if
            write (iout,150)  istep,ekvir+epotpi_loc,epotpi_loc,ekvir
     &             ,ekprim,temppi,temppi_cl,presvir
150         format (i10,4f14.4,3f11.2)
          endif

        end if

c       COMPUTE POLYMER TRANSFER MATRIX (TO BE REPLACED BY FFT !!!)
!        allocate(eigmat(nbeads,nbeads),eigmattr(nbeads,nbeads))
!        allocate(omkpi(nbeads),WORK(3*nbeads))
!        eigmat=0
!        DO i=1,nbeads-1
!          eigmat(i,i)=2
!          eigmat(i+1,i)=-1
c
c       print header for the averages over a group of recent steps
c
        if (modstep .eq. 0) then
           pico = dble(istep) * dt
           write (iout,250)  iprint,istep
  250      format (/,' Average Values for the last',i6,' out of',
     &                i9,' Dynamics Steps')
           write (iout,260)  pico
  260      format (/,' Simulation Time',5x,f15.4,' Picosecond')
        end if
c
c       compute total energy and fluctuation for recent steps
c
        etotppi=ekvir+epotpi_loc
        etotpi_sum = etotpi_sum + etotppi 
        etot2pi_sum = etot2pi_sum +  etotppi**2
        if (modstep .eq. 0) then
           etotpi_ave = etotpi_sum / dble(iprint)
           etot2pi_ave = etot2pi_sum / dble(iprint)
           fluctuate2pi = etot2pi_ave - etotpi_ave**2
           if (fluctuate2pi .gt. 0.0d0) then
              fluctuatepi = sqrt(fluctuate2pi)
           else
              fluctuatepi = 0.0d0
           end if
           write (iout,70)  etotpi_ave,fluctuatepi
   70      format (' Total Energy',8x,f15.4,' Kcal/mole',3x,
     &                '(+/-',f9.4,')')
        end if

c       compute average potential energy and its fluctuation

        epotpi_sum = epotpi_sum + epotpi_loc
        epot2pi_sum = epot2pi_sum + epotpi_loc**2
        if (modstep .eq. 0) then
           epotpi_ave = epotpi_sum / dble(iprint)
           epot2pi_ave = epot2pi_sum / dble(iprint)
           potfluct2pi = epot2pi_ave - epotpi_ave**2
           if (potfluct2pi .gt. 0.0d0) then
              potfluctpi = sqrt(potfluct2pi)
           else
              potfluctpi = 0.0d0
           end if
           write (iout,80)  epotpi_ave,potfluctpi
   80      format (' Potential Energy',4x,f15.4,' Kcal/mole',3x,
     &                '(+/-',f9.4,')')
        end if
c
c       compute average kinetic energy and its fluctuation
c
        ekinpi_sum = ekinpi_sum + ekvir
        ekin2pi_sum = ekin2pi_sum + ekvir**2
        if (modstep .eq. 0) then
           ekinpi_ave = ekinpi_sum / dble(iprint)
           ekin2pi_ave = ekin2pi_sum / dble(iprint)
           kinfluct2pi = ekin2pi_ave - ekinpi_ave**2
           if (kinfluct2pi .gt. 0.0d0) then
              kinfluctpi = sqrt(kinfluct2pi)
           else
              kinfluctpi = 0.0d0
           end if
           write (iout,90)  ekinpi_ave,kinfluctpi
   90      format (' Kinetic Energy',6x,f15.4,' Kcal/mole',3x,
     &                '(+/-',f9.4,')')
        end if
c
c       compute average intermolecular energy and its fluctuation
c
c        if (nmol.ne.1 .and. nmol.ne.n .and. .not.use_ewald) then
c           eint_sum = eint_sum + einter
c           eint2_sum = eint2_sum + einter**2
c           if (modstep .eq. 0) then
c              eint_ave = eintpi_sum / dble(iprint)
c              eint2_ave = eint2_sum / dble(iprint)
c              intfluct2 = eint2_ave - eint_ave**2
c              if (intfluct2 .gt. 0.0d0) then
c                 intfluct = sqrt(intfluct2)
c              else
c                 intfluct = 0.0d0
c              end if
c              write (iout,100)  eint_ave,intfluct
c  100         format (' Intermolecular',6x,f15.4,' Kcal/mole',3x,
c     &                   '(+/-',f9.4,')')
c           end if
c        end if
c
c       compute the average temperature and its fluctuation
c
        temppi_sum = temppi_sum + temppi
        temp2pi_sum = temp2pi_sum + temppi**2
        if (modstep .eq. 0) then
           temppi_ave = temppi_sum / dble(iprint)
           temp2pi_ave = temp2pi_sum / dble(iprint)
           tfluct2pi = temp2pi_ave - temppi_ave**2
           if (tfluct2pi .gt. 0.0d0) then
              tfluctpi = sqrt(tfluct2pi)
           else
              tfluctpi = 0.0d0
           end if
           write (iout,310)  temppi_ave,tfluctpi
  310      format (' Temperature',9x,f15.2,' Kelvin',6x,
     &                '(+/-',f9.2,')')
        end if

c       compute the average pressure and its fluctuation
c
        if (use_bounds) then
           prespi_sum = prespi_sum + presvir
           pres2pi_sum = pres2pi_sum + presvir**2
           if (modstep .eq. 0) then
              prespi_ave = prespi_sum / dble(iprint)
              pres2pi_ave = pres2pi_sum / dble(iprint)
              pfluct2pi = pres2pi_ave - prespi_ave**2
              if (pfluct2pi .gt. 0.0d0) then
                 pfluctpi = sqrt(pfluct2pi)
              else
                 pfluctpi = 0.0d0
              end if
              write (iout,120)  prespi_ave,pfluctpi
  120         format (' Pressure',12x,f15.2,' Atmosphere',2x,
     &                   '(+/-',f9.2,')')
           end if

c       compute the average density and its fluctuation
c
           denspi_sum = denspi_sum + dens
           dens2pi_sum = dens2pi_sum + dens**2
           if (modstep .eq. 0) then
              denspi_ave = denspi_sum / dble(iprint)
              dens2pi_ave = dens2pi_sum / dble(iprint)
              dfluct2pi = dens2pi_ave - denspi_ave**2
              if (dfluct2pi .gt. 0.0d0) then
                 dfluctpi = sqrt(dfluct2pi)
              else
                 dfluctpi = 0.0d0
              end if
              write (iout,130)  denspi_ave,dfluctpi
  130         format (' Density',13x,f15.4,' Grams/cc',4x,
     &                   '(+/-',f9.4,')')
           end if

c       compute the average volume and its fluctuation
c
           volpi_sum = volpi_sum + extvol
           vol2pi_sum = vol2pi_sum + extvol**2
           if (modstep .eq. 0) then
              volpi_ave = volpi_sum / dble(iprint)
              vol2pi_ave = vol2pi_sum / dble(iprint)
              volfluct2pi = vol2pi_ave - volpi_ave**2
              if (volfluct2pi .gt. 0.0d0) then
                 volfluctpi = sqrt(volfluct2pi)
              else
                 volfluctpi = 0.0d0
              end if
              write (iout,350)  volpi_ave,volfluctpi
  350         format (' Volume',13x,f15.4,' Angstrom^3',3x,
     &                   '(+/-',f9.4,')')
           end if

c       compute the average piston temperature and its fluctuation
c
           if ((barostat.eq.'LANGEVIN').and.(rank.eq.0)) then
             pistontemppi = (masspiston*vextvol*vextvol)/boltzmann
             pistontemppi_sum = pistontemppi_sum + pistontemppi
             pistontemp2pi_sum = pistontemp2pi_sum + pistontemppi**2
             if (modstep .eq. 0) then
                pistontemppi_ave = pistontemppi_sum / dble(iprint)
                pistontemp2pi_ave = pistontemp2pi_sum / dble(iprint)
                dfluct2pi = pistontemp2pi_ave - pistontemppi_ave**2
                if (dfluct2pi .gt. 0.0d0) then
                   dfluctpi = sqrt(dfluct2pi)
                else
                   dfluctpi = 0.0d0
                end if
                write (iout,131)  pistontemppi_ave,dfluctpi
  131           format (' Piston temperature',1x,f15.4,' K',12x,
     &                     '(+/-',f9.4,')')
             end if
           end if
        end if
c
c
c     display the different times of the computation
cc
c      if ((rank.eq.0).and.(verbose)) then
c        if (modstep.eq.0) then
cc          write(iout,150) timeclear/dble(iprint)
cc 150      format('Ave time for aclear comms            =  ',F24.15)
cc          timeclear = 0.0d0
c          write(iout,160) timereneig/dble(iprint)
c 160      format('Ave time for reneig                  =  ',F24.15)
c          timereneig = 0.0d0
c          write(iout,170) timecommstep/dble(iprint)
c 170      format('Ave time for positions comm          =  ',F24.15)
c          timecommstep = 0.0d0
c          write(iout,180) timeparam/dble(iprint)
c 180      format('Ave time for param                   =  ',F24.15)
c          timeparam = 0.0d0
c          write(iout,190) timeforcescomm/dble(iprint)
c 190      format('Ave time for forces comm             =  ',F24.15)
c          timeforcescomm = 0.0d0
c          write(iout,210) timedirreccomm/dble(iprint)
c 210      format('Ave time for reciprocal forces comm  =  ',F24.15)
c          timedirreccomm = 0.0d0
c          write(iout,220) timebondedvdw/dble(iprint)
c 220      format('Ave time for real space forces comm  =  ',F24.15)
c          timebondedvdw = 0.0d0
c          write(iout,230) timebonded/dble(iprint)
c 230      format('Ave time for bonded forces           =  ',F24.15)
c          timebonded = 0.0d0
c          write(iout,240) timenonbonded/dble(iprint)
c 240      format('Ave time for non bonded forces       =  ',F24.15)
c          timenonbonded = 0.0d0
c          write(iout,250) timenl/dble(iprint)
c 250      format('Ave time for neighbor lists          =  ',F24.15)
c          timenl = 0.0d0
c          write(iout,260) timereal/dble(iprint)
c 260      format('Ave time for real space (permanent)  =  ',F24.15)
c          timereal = 0.0d0
c          write(iout,270) timerealdip/dble(iprint)
c 270      format('Ave time for real space (polar)      =  ',F24.15)
c          timerealdip = 0.0d0
c          write(iout,280) timegrid1/dble(iprint)
c 280      format('Ave time for fill grid (permanent)   =  ',F24.15)
c          timegrid1 = 0.0d0
c          write(iout,290) timeffts/dble(iprint)
c 290      format('Ave time for ffts (permanent)        =  ',F24.15)
c          timeffts = 0.0d0
c          write(iout,300) timescalar/dble(iprint)
c 300      format('Ave time for scalar prod (permanent) =  ',F24.15)
c          timescalar = 0.0d0
c          write(iout,310) timegrid2/dble(iprint)
c 310      format('Ave time for extract grid (permanent)=  ',F24.15)
c          timegrid2 = 0.0d0
c          write(iout,320) timerecreccomm/dble(iprint)
c 320      format('Ave time for rec-rec comm (permanent)=  ',F24.15)
c          timerecreccomm = 0.0d0
c          write(iout,330) timerec/dble(iprint)
c 330      format('Ave time for recip space (permanent) =  ',F24.15)
c          timerec = 0.0d0
c          write(iout,340) timerecdip/dble(iprint)
c 340      format('Ave time for recip space (polar)     =  ',F24.15)
c          timerecdip = 0.0d0
c        end if
c      end if
cc
c      if (modstep.eq.0) then
cc
cc     put timing variables to zero
cc
c        timeclear = 0.0d0
c        timereneig = 0.0d0
c        timecommstep = 0.0d0
c        timeparam = 0.0d0
c        timeforcescomm = 0.0d0
c        timedirreccomm = 0.0d0
c        timebondedvdw = 0.0d0
c        timebonded = 0.0d0
c        timenonbonded = 0.0d0
c        timereal = 0.0d0
c        timerealdip = 0.0d0
c        timegrid1 = 0.0d0
c        timeffts = 0.0d0
c        timescalar = 0.0d0
c        timegrid2 = 0.0d0
c        timerecreccomm = 0.0d0
c        timerec = 0.0d0
c        timerecdip = 0.0d0
     
      return
      end
