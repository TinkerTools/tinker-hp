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
      subroutine mdstatpi (istep,dt)
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
      use timestat
      use units
      use mpi
      use virial
      implicit none
      integer, intent(in) ::  istep
      real*8, intent(in) :: dt
      integer :: mode,modstep,ierr,i
      real*8 :: temp,pres
      real*8 :: etotppi,epot,ekin,vol
      real*8 :: pico,dens,pistontemppi
      real*8 :: fluctuatepi,fluctuate2pi
      real*8 :: intfluctpi,intfluct2pi
      real*8 :: potfluctpi,potfluct2pi
      real*8 :: kinfluctpi,kinfluct2pi
      real*8 :: tfluctpi,pfluctpi,dfluctpi,volfluctpi
      real*8 :: tfluct2pi,pfluct2pi,dfluct2pi,volfluct2pi
      real*8 :: gyrfluctpi,gyrfluct2pi
      real*8 :: etotpi_ave,etot2pi_ave
      real*8 :: Tcentroidfluct,Tcentroidfluct2,Tcentroid
      real*8 :: eintpi_ave,eint2pi_ave
      real*8 :: epotpi_ave,epot2pi_ave
      real*8 :: ekinpi_ave,ekin2pi_ave
      real*8 :: temppi_ave,temp2pi_ave
      real*8 :: prespi_ave,pres2pi_ave
      real*8 :: denspi_ave,dens2pi_ave
      real*8 :: volpi_ave,vol2pi_ave
      real*8 :: gyrpi_ave,gyr2pi_ave
      real*8 :: Tcentroid_ave,Tcentroid2_ave
      real*8 :: pistontemppi_ave,pistontemp2pi_ave
      real*8,save :: box_ave(3),box2_ave(3)
      real*8,save :: stress_ave(3,3),stress2_ave(3,3)
      real*8,save :: eintpi_sum,eint2pi_sum
      real*8,save :: epotpi_sum,epot2pi_sum
      real*8,save :: etotpi_sum,etot2pi_sum
      real*8,save :: ekinpi_sum,ekin2pi_sum
      real*8,save :: temppi_sum,temp2pi_sum
      real*8,save :: prespi_sum,pres2pi_sum
      real*8,save :: denspi_sum,dens2pi_sum
      real*8,save :: volpi_sum,vol2pi_sum
      real*8,save :: gyrpi_sum,gyr2pi_sum
      real*8,save :: Tcentroid_sum,Tcentroid2_sum
      real*8,save :: pistontemppi_sum,pistontemp2pi_sum
      real*8,save :: box_sum(3),box2_sum(3)
      real*8,save :: stress_sum(3,3),stress2_sum(3,3)
      real*8 :: lth 

      if(ranktot/=0) return
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
         gyrpi_sum = 0.0d0
         gyr2pi_sum = 0.0d0
         Tcentroid_sum = 0.0d0
         Tcentroid2_sum = 0.0d0
         pistontemppi_sum = 0.0d0
         pistontemp2pi_sum = 0.0d0
         box_sum = 0.0d0
         box2_sum = 0.0d0
         stress_sum = 0.0d0
         stress2_sum = 0.0d0
      end if

        dens = (1.0d24/real(volbox,8))
     &          * (real(totmass,8)/real(avogadro,8))
        if (verbose) then
          if(modstep==1) then
            write(iout,'(A)',advance="no") "    MD Step"
     &                  //"      E total"
     &                  //"   E Potential"            
     &                  //"     E Kinetic"
c     &                  //"     Ek prim"
     &                  //"       Temp"
     &                  //"    Temp_cl"
            write(iout,'(A)',advance="no") "       Pres"
c            if(isobaric) then
c              write(iout,'(A)',advance="no") 
c     &                    "     Density"
c     &                  //"      Volume"
c            endif
            write(iout,*)
          endif

          write(iout,'(i10,3f14.4,2f11.2)',advance="no") istep
     &                ,ekvir+epotpi,epotpi,ekvir !,ekprim
     &                ,temppi,temppi_cl/nbeads
          write(iout,'(f11.2)',advance="no") presvir
c           if(isobaric) then
c             write(iout,'(f12.4,f12.2)',advance="no") dens,volbox
c           endif
          write(iout,*)
        end if

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
        etotppi=ekvir+epotpi
        etotpi_sum = etotpi_sum + etotppi 
        etot2pi_sum = etot2pi_sum +  etotppi**2
        if (modstep .eq. 0) then
           etotpi_ave = etotpi_sum / dble(iprint)
           etot2pi_ave = etot2pi_sum / dble(iprint)
           fluctuate2pi = etot2pi_ave - etotpi_ave**2
           if (fluctuate2pi .gt. 0.0_8) then
              fluctuatepi = sqrt(fluctuate2pi)
           else
              fluctuatepi = 0.0_8
           end if
           write (iout,70)  etotpi_ave,fluctuatepi
   70      format (' Total Energy',8x,f15.4,' Kcal/mole',3x,
     &                '(+/-',f9.4,')')
        end if

c       compute average potential energy and its fluctuation

        epotpi_sum = epotpi_sum + epotpi
        epot2pi_sum = epot2pi_sum + epotpi**2
        if (modstep .eq. 0) then
           epotpi_ave = epotpi_sum / dble(iprint)
           epot2pi_ave = epot2pi_sum / dble(iprint)
           potfluct2pi = epot2pi_ave - epotpi_ave**2
           if (potfluct2pi .gt. 0.0_8) then
              potfluctpi = sqrt(potfluct2pi)
           else
              potfluctpi = 0.0_8
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
           if (kinfluct2pi .gt. 0.0_8) then
              kinfluctpi = sqrt(kinfluct2pi)
           else
              kinfluctpi = 0.0_8
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
c              if (intfluct2 .gt. 0.0_8) then
c                 intfluct = sqrt(intfluct2)
c              else
c                 intfluct = 0.0_8
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
           if (tfluct2pi .gt. 0.0_8) then
              tfluctpi = sqrt(tfluct2pi)
           else
              tfluctpi = 0.0_8
           end if
           write (iout,310)  temppi_ave,tfluctpi
  310      format (' Temperature',9x,f15.2,' Kelvin',6x,
     &                '(+/-',f9.2,')')
        end if

        gyrpi_sum = gyrpi_sum + gyrpi
        gyr2pi_sum = gyr2pi_sum + gyrpi**2
        if (modstep .eq. 0) then
           gyrpi_ave = gyrpi_sum  / dble(iprint)
           gyr2pi_ave = gyr2pi_sum / dble(iprint)
           gyrfluct2pi = gyr2pi_ave - gyrpi_ave**2
           if (gyrfluct2pi .gt. 0.0d0) then
              gyrfluctpi = sqrt(gyrfluct2pi)
           else
              gyrfluctpi = 0.0d0
           end if
           lth=hbar_planck/sqrt(boltzmann*kelvin)
           write (iout,360)  gyrpi_ave,gyrfluctpi
  360      format (' Gyration radius',5x,f15.4,' Angstrom',4x,
     &                '(+/-',f9.4,')')
        end if

        Tcentroid = ekcentroid*2.d0/(3*n*gasconst)
        Tcentroid_sum = Tcentroid_sum + Tcentroid
        Tcentroid2_sum = Tcentroid2_sum + Tcentroid**2
        if (modstep .eq. 0) then
           Tcentroid_ave = Tcentroid_sum / dble(iprint)
           Tcentroid2_ave = Tcentroid2_sum / dble(iprint)
           Tcentroidfluct2 = Tcentroid2_ave - Tcentroid_ave**2
           if (gyrfluct2pi .gt. 0.0d0) then
              Tcentroidfluct = sqrt(Tcentroidfluct2)
           else
              Tcentroidfluct = 0.0d0
           end if
           write (iout,370)  Tcentroid_ave,Tcentroidfluct
  370      format (' Centroid Temperature',5x,f15.4,' K',4x,
     &                '(+/-',f9.4,')')
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
              if (pfluct2pi .gt. 0.0_8) then
                 pfluctpi = sqrt(pfluct2pi)
              else
                 pfluctpi = 0.0_8
              end if
              write (iout,120)  prespi_ave,pfluctpi
  120         format (' Pressure',12x,f15.2,' Atmosphere',2x,
     &                   '(+/-',f9.2,')')
          endif

c       compute the average density and its fluctuation
c
           denspi_sum = denspi_sum + dens
           dens2pi_sum = dens2pi_sum + dens**2
           if (modstep .eq. 0) then
              denspi_ave = denspi_sum / dble(iprint)
              dens2pi_ave = dens2pi_sum / dble(iprint)
              dfluct2pi = dens2pi_ave - denspi_ave**2
              if (dfluct2pi .gt. 0.0_8) then
                 dfluctpi = sqrt(dfluct2pi)
              else
                 dfluctpi = 0.0_8
              end if
              write (iout,130)  denspi_ave,dfluctpi
  130         format (' Density',13x,f15.4,' Grams/cc',4x,
     &                   '(+/-',f9.4,')')
           end if

c       compute the average volume and its fluctuation
c
           volpi_sum = volpi_sum + volbox 
           vol2pi_sum = vol2pi_sum + volbox**2
           if (modstep .eq. 0) then
              volpi_ave = volpi_sum / dble(iprint)
              vol2pi_ave = vol2pi_sum / dble(iprint)
              volfluct2pi = vol2pi_ave - volpi_ave**2
              if (volfluct2pi .gt. 0.0_8) then
                 volfluctpi = sqrt(volfluct2pi)
              else
                 volfluctpi = 0.0_8
              end if
              write (iout,350)  volpi_ave,volfluctpi
  350         format (' Volume',13x,f15.4,' Angstrom^3',3x,
     &                   '(+/-',f9.4,')')
           end if

c       compute the average piston temperature and its fluctuation
c
           if (isobaric.and.(barostat.eq.'LANGEVIN')) then
             ! pistontemppi = masspiston*vextvol**2/boltzmann
             pistontemppi_sum = pistontemppi_sum + pistontemppi
             pistontemp2pi_sum = pistontemp2pi_sum + pistontemppi**2
             if (modstep .eq. 0) then
                pistontemppi_ave = pistontemppi_sum / dble(iprint)
                pistontemp2pi_ave = pistontemp2pi_sum / dble(iprint)
                dfluct2pi = pistontemp2pi_ave - pistontemppi_ave**2
                if (dfluct2pi .gt. 0.0_8) then
                   dfluctpi = sqrt(dfluct2pi)
                else
                   dfluctpi = 0.0_8
                end if
                write (iout,131)  pistontemppi_ave,dfluctpi
  131           format (' Piston Temperature',1x,f15.4,' K',12x,
     &                     '(+/-',f9.4,')')
             end if
           end if

           if (isobaric .and. anisotrop) then
             stress_sum = stress_sum + stresspi
             stress2_sum = stress2_sum + stresspi**2
             box_sum(1) = box_sum(1) + xbox
             box_sum(2) = box_sum(2) + ybox
             box_sum(3) = box_sum(3) + zbox
             box2_sum(1) = box2_sum(1) + xbox**2
             box2_sum(2) = box2_sum(2) + ybox**2
             box2_sum(3) = box2_sum(3) + zbox**2
             if (modstep .eq. 0) then
              box_ave =  box_sum / dble(iprint)
              box2_ave = box2_sum / dble(iprint)
              write (iout,'(A,6x,3f15.4,A)') ' Box parameters'
     &                ,box_ave,'  Angstrom'  
              stress_ave =  stress_sum / dble(iprint)
              stress2_ave = stress2_sum / dble(iprint)
              write (iout,'(A,A,2x,6f15.4,A)') ' Stress Tensor '
     &                ,'(xx yy zz xy yz xz)'
     &                ,stress_ave(1,1),stress_ave(2,2),stress_ave(3,3)
     &                ,stress_ave(1,2),stress_ave(2,3),stress_ave(1,3)
     &                ,' Atmosphere'
             end if
           end if
        end if

      return
      end
