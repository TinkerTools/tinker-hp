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
#include "tinker_precision.h"
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
      real(r_p), intent(in) :: dt
      integer,parameter:: d_prec= kind(1.0d0)
      integer mode,modstep,ierr,i
      integer :: period,freq
      logical :: display
      real(r_p) temp,pres
      real(r_p) etotppi,epot,ekin,vol
      real(d_prec) pico,dens,pistontemppi
      real(d_prec) fluctuatepi,fluctuate2pi
      real(d_prec) intfluctpi,intfluct2pi
      real(d_prec) potfluctpi,potfluct2pi
      real(d_prec) kinfluctpi,kinfluct2pi
      real(d_prec) tfluctpi,pfluctpi,dfluctpi,volfluctpi
      real(d_prec) tfluct2pi,pfluct2pi,dfluct2pi,volfluct2pi
      real(d_prec) gyrfluctpi,gyrfluct2pi
      real(d_prec) Tcentroidfluct,Tcentroidfluct2,Tcentroid
      real(d_prec) etotpi_ave,etot2pi_ave
      real(d_prec) eintpi_ave,eint2pi_ave
      real(d_prec) epotpi_ave,epot2pi_ave
      real(d_prec) ekinpi_ave,ekin2pi_ave
      real(d_prec) temppi_ave,temp2pi_ave
      real(d_prec) prespi_ave,pres2pi_ave
      real(d_prec) denspi_ave,dens2pi_ave
      real(d_prec) volpi_ave,vol2pi_ave
      real(d_prec) gyrpi_ave,gyr2pi_ave
      real(d_prec) Tcentroid_ave,Tcentroid2_ave
      real(d_prec) pistontemppi_ave,pistontemp2pi_ave
      real(d_prec),save :: box_ave(3),box2_ave(3)
      real(d_prec),save :: stress_ave(3,3),stress2_ave(3,3)
      real(d_prec),save :: etotpi_sum,etot2pi_sum
      real(d_prec),save :: eintpi_sum,eint2pi_sum
      real(d_prec),save :: epotpi_sum,epot2pi_sum
      real(d_prec),save :: ekinpi_sum,ekin2pi_sum
      real(d_prec),save :: temppi_sum,temp2pi_sum
      real(d_prec),save :: prespi_sum,pres2pi_sum
      real(d_prec),save :: denspi_sum,dens2pi_sum
      real(d_prec),save :: volpi_sum,vol2pi_sum
      real(d_prec),save :: gyrpi_sum,gyr2pi_sum
      real(d_prec),save :: Tcentroid_sum,Tcentroid2_sum
      real(d_prec),save :: pistontemppi_sum,pistontemp2pi_sum
      real(d_prec),save :: box_sum(3),box2_sum(3)
      real(d_prec),save :: stress_sum(3,3),stress2_sum(3,3)
      integer, save  :: count=0
      real(d_prec) :: lth 
#ifdef SINGLE
      character(1) :: edcim = '2'
#else
      character(1) :: edcim = '4'
#endif
      character(2) :: esize
      character(5) :: f_ener

      esize='14'
      if (n>5d6) esize='17'
      f_ener = 'f'//esize//'.'//edcim

c
c
c     set number of steps for block averages of properties
c
#if (defined(SINGLE) | defined(MIXED))
      freq    = 100
      period  = max(1,iprint/freq)
#else
      freq    = iprint
      period  = 1
#endif
      modstep = mod(istep,iprint)
      display = (mod(istep,period).eq.0.and.verbose)

      call inquire_calc_e(istep+1,period)

      if(ranktot/=0) return

c
c     zero out summation variables for new averaging period
c
      if (modstep.eq.1 .or. iprint.eq.1) then
         count=0
         etotpi_sum = 0.0_d_prec
         etot2pi_sum = 0.0_d_prec
         epotpi_sum = 0.0_d_prec
         epot2pi_sum = 0.0_d_prec
         ekinpi_sum = 0.0_d_prec
         ekin2pi_sum = 0.0_d_prec
         eintpi_sum = 0.0_d_prec
         eint2pi_sum = 0.0_d_prec
         temppi_sum = 0.0_d_prec
         temp2pi_sum = 0.0_d_prec
         prespi_sum = 0.0_d_prec
         pres2pi_sum = 0.0_d_prec
         denspi_sum = 0.0_d_prec
         dens2pi_sum = 0.0_d_prec
         volpi_sum = 0.0_d_prec
         vol2pi_sum = 0.0_d_prec
         gyrpi_sum = 0.0_d_prec
         gyr2pi_sum = 0.0_d_prec
         Tcentroid_sum = 0.0_d_prec
         Tcentroid2_sum = 0.0_d_prec
         pistontemppi_sum = 0.0_d_prec
         pistontemp2pi_sum = 0.0_d_prec
         box_sum = 0.0_d_prec
         box2_sum = 0.0_d_prec
         stress_sum = 0.0_d_prec
         stress2_sum = 0.0_d_prec
      end if
      if(display) then
        count=count+1
        dens = (1.0d24/real(volbox,d_prec))
     &          * (real(totmass,d_prec)/real(avogadro,d_prec))
!$acc update host(epotpi,eintrapi,einterpi,presvir
!$acc&   ,ekvir,ekprim,temppi,temppi_cl,gyrpi,stresspi,ekcentroid) async
!$acc wait
      endif

        !dens = (1.0d24/extvol) * (totmass/avogadro)
c        dens = (1.0d24/real(volbox,d_prec))
c     &          * (real(totmass,d_prec)/real(avogadro,d_prec))
        if (verbose) then
          if(modstep==1) then
            write(iout,'(A)',advance="no") "    MD Step"
     &                  //"      E total"
     &                  //"   E Potential"            
     &                  //"     E Kinetic"
c     &                  //"     Ek prim"
     &                  //"       Temp"
     &                  //"    Temp_cl"
            if(use_virial) then
              write(iout,'(A)',advance="no") "       Pres"
            endif
c            if(isobaric) then
c              write(iout,'(A)',advance="no") 
c     &                    "     Density"
c     &                  //"      Volume"
c            endif
            write(iout,*)
          endif
          if (display) then
            write(iout,'(i10,3'//f_ener//',2f11.2)',advance="no") istep
     &                 ,ekvir+epotpi,epotpi,ekvir !,ekprim
     &                 ,temppi,temppi_cl/nbeads
            if(use_virial) then
              write(iout,'(f11.2)',advance="no") presvir
            endif
c            if(isobaric) then
c              write(iout,'(f12.4,f12.2)',advance="no") dens,volbox
c            endif
            write(iout,*)
          endif
        end if

c
c       print header for the averages over a group of recent steps
c
        if (verbose .and. modstep .eq. 0) then
!$acc wait
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
        if (display) then
          etotppi=ekvir+epotpi
          etotpi_sum = etotpi_sum + etotppi 
          etot2pi_sum = etot2pi_sum +  etotppi**2
        endif
        if (verbose.and.modstep .eq. 0) then
           etotpi_ave = etotpi_sum / real(count,d_prec)
           etot2pi_ave = etot2pi_sum / real(count,d_prec)
           fluctuate2pi = etot2pi_ave - etotpi_ave**2
           if (fluctuate2pi .gt. 0.0_d_prec) then
              fluctuatepi = sqrt(fluctuate2pi)
           else
              fluctuatepi = 0.0_d_prec
           end if
           write (iout,70)  etotpi_ave,fluctuatepi
   70      format (' Total Energy',8x,f15.4,' Kcal/mole',3x,
     &                '(+/-',f9.4,')')
        end if

c       compute average potential energy and its fluctuation

        if (display) then
          epotpi_sum = epotpi_sum + epotpi
          epot2pi_sum = epot2pi_sum + epotpi**2
        endif
        if (verbose.and.modstep .eq. 0) then
           epotpi_ave = epotpi_sum / real(count,d_prec)
           epot2pi_ave = epot2pi_sum / real(count,d_prec)
           potfluct2pi = epot2pi_ave - epotpi_ave**2
           if (potfluct2pi .gt. 0.0_d_prec) then
              potfluctpi = sqrt(potfluct2pi)
           else
              potfluctpi = 0.0_d_prec
           end if
           write (iout,80)  epotpi_ave,potfluctpi
   80      format (' Potential Energy',4x,f15.4,' Kcal/mole',3x,
     &                '(+/-',f9.4,')')
        end if
c
c       compute average kinetic energy and its fluctuation
c
        if (display) then
        ekinpi_sum = ekinpi_sum + ekvir
        ekin2pi_sum = ekin2pi_sum + ekvir**2
        endif
        if (verbose.and.modstep .eq. 0) then
           ekinpi_ave = ekinpi_sum / real(count,d_prec)
           ekin2pi_ave = ekin2pi_sum / real(count,d_prec)
           kinfluct2pi = ekin2pi_ave - ekinpi_ave**2
           if (kinfluct2pi .gt. 0.0_d_prec) then
              kinfluctpi = sqrt(kinfluct2pi)
           else
              kinfluctpi = 0.0_d_prec
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
c              if (intfluct2 .gt. 0.0_d_prec) then
c                 intfluct = sqrt(intfluct2)
c              else
c                 intfluct = 0.0_d_prec
c              end if
c              write (iout,100)  eint_ave,intfluct
c  100         format (' Intermolecular',6x,f15.4,' Kcal/mole',3x,
c     &                   '(+/-',f9.4,')')
c           end if
c        end if
c
c       compute the average temperature and its fluctuation
c
        if (display) then
        temppi_sum = temppi_sum + temppi
        temp2pi_sum = temp2pi_sum + temppi**2
        endif
        if (verbose.and.modstep .eq. 0) then
           temppi_ave = temppi_sum / real(count,d_prec)
           temp2pi_ave = temp2pi_sum / real(count,d_prec)
           tfluct2pi = temp2pi_ave - temppi_ave**2
           if (tfluct2pi .gt. 0.0_d_prec) then
              tfluctpi = sqrt(tfluct2pi)
           else
              tfluctpi = 0.0_d_prec
           end if
           write (iout,310)  temppi_ave,tfluctpi
  310      format (' Temperature',9x,f15.2,' Kelvin',6x,
     &                '(+/-',f9.2,')')
        end if

        if (display) then
        gyrpi_sum = gyrpi_sum + gyrpi
        gyr2pi_sum = gyr2pi_sum + gyrpi**2
        endif
        if (verbose.and.modstep .eq. 0) then
           gyrpi_ave = gyrpi_sum / real(count,d_prec)
           gyr2pi_ave = gyr2pi_sum / real(count,d_prec)
           gyrfluct2pi = gyr2pi_ave - gyrpi_ave**2
           if (gyrfluct2pi .gt. 0.0_d_prec) then
              gyrfluctpi = sqrt(gyrfluct2pi)
           else
              gyrfluctpi = 0.0_d_prec
           end if
           lth=hbar_planck/sqrt(boltzmann*kelvin)
           write (iout,360)  gyrpi_ave,gyrfluctpi
  360      format (' Gyration radius',5x,f15.4,' Angstrom',4x,
     &                '(+/-',f9.4,')')
        end if

        if (display) then
        Tcentroid = ekcentroid*2.d0/(3*n*gasconst)
        Tcentroid_sum = Tcentroid_sum + Tcentroid
        Tcentroid2_sum = Tcentroid2_sum + Tcentroid**2
        endif
        if (verbose.and.modstep .eq. 0) then
           Tcentroid_ave = Tcentroid_sum / real(count,d_prec)
           Tcentroid2_ave = Tcentroid2_sum / real(count,d_prec)
           Tcentroidfluct2 = Tcentroid2_ave - Tcentroid_ave**2
           if (gyrfluct2pi .gt. 0.0_d_prec) then
              Tcentroidfluct = sqrt(Tcentroidfluct2)
           else
              Tcentroidfluct = 0.0_d_prec
           end if
           write (iout,370)  Tcentroid_ave,Tcentroidfluct
  370      format (' Centroid Temperature',5x,f15.4,' K',4x,
     &                '(+/-',f9.4,')')
        end if

c       compute the average pressure and its fluctuation
c
        if (use_bounds) then
          if(use_virial) then
            if (display) then
            prespi_sum = prespi_sum + presvir
            pres2pi_sum = pres2pi_sum + presvir**2
            endif
            if (verbose.and.modstep .eq. 0) then
              prespi_ave = prespi_sum / real(count,d_prec)
              pres2pi_ave = pres2pi_sum / real(count,d_prec)
              pfluct2pi = pres2pi_ave - prespi_ave**2
              if (pfluct2pi .gt. 0.0_d_prec) then
                 pfluctpi = sqrt(pfluct2pi)
              else
                 pfluctpi = 0.0_d_prec
              end if
              write (iout,120)  prespi_ave,pfluctpi
  120         format (' Pressure',12x,f15.2,' Atmosphere',2x,
     &                   '(+/-',f9.2,')')
            end if
          endif

c       compute the average density and its fluctuation
c
            if (display) then
            denspi_sum = denspi_sum + dens
            dens2pi_sum = dens2pi_sum + dens**2
            endif
            if (verbose.and.modstep .eq. 0) then
              denspi_ave = denspi_sum / real(count,d_prec)
              dens2pi_ave = dens2pi_sum / real(count,d_prec)
              dfluct2pi = dens2pi_ave - denspi_ave**2
              if (dfluct2pi .gt. 0.0_d_prec) then
                 dfluctpi = sqrt(dfluct2pi)
              else
                 dfluctpi = 0.0_d_prec
              end if
              write (iout,130)  denspi_ave,dfluctpi
  130         format (' Density',13x,f15.4,' Grams/cc',4x,
     &                   '(+/-',f9.4,')')
           end if

c       compute the average volume and its fluctuation
c
            if (display) then
           volpi_sum = volpi_sum + volbox 
           vol2pi_sum = vol2pi_sum + volbox**2
           endif
           if (verbose.and.modstep .eq. 0) then
              volpi_ave = volpi_sum / real(count,d_prec)
              vol2pi_ave = vol2pi_sum / real(count,d_prec)
              volfluct2pi = vol2pi_ave - volpi_ave**2
              if (volfluct2pi .gt. 0.0_d_prec) then
                 volfluctpi = sqrt(volfluct2pi)
              else
                 volfluctpi = 0.0_d_prec
              end if
              write (iout,350)  volpi_ave,volfluctpi
  350         format (' Volume',13x,f15.4,' Angstrom^3',3x,
     &                   '(+/-',f9.4,')')
           end if

c       compute the average piston temperature and its fluctuation
c
           if (isobaric .and.(barostat.eq.'LANGEVIN')) then
              if (display) then
             !pistontemppi = masspiston*vextvol**2/boltzmann
             pistontemppi_sum = pistontemppi_sum + temppiston
             pistontemp2pi_sum = pistontemp2pi_sum + temppiston**2
             endif
             if (verbose.and.modstep .eq. 0) then
              pistontemppi_ave=pistontemppi_sum/real(count,d_prec)
              pistontemp2pi_ave=pistontemp2pi_sum/real(count,d_prec)
              dfluct2pi = pistontemp2pi_ave - pistontemppi_ave**2
              if (dfluct2pi .gt. 0.0_d_prec) then
                 dfluctpi = sqrt(dfluct2pi)
              else
                 dfluctpi = 0.0_d_prec
              end if
              write (iout,131)  pistontemppi_ave,dfluctpi
  131          format (' Piston Temperature',1x,f15.4,' K',12x,
     &                    '(+/-',f9.4,')')
             end if
           end if
          
          if (isobaric .and. anisotrop) then
              if (display) then
             stress_sum = stress_sum + stresspi
             stress2_sum = stress2_sum + stresspi**2
             box_sum(1) = box_sum(1) + xbox
             box_sum(2) = box_sum(2) + ybox
             box_sum(3) = box_sum(3) + zbox
             box2_sum(1) = box2_sum(1) + xbox**2
             box2_sum(2) = box2_sum(2) + ybox**2
             box2_sum(3) = box2_sum(3) + zbox**2
             endif
             if (verbose.and.modstep .eq. 0) then
              box_ave =  box_sum / real(count,d_prec)
              box2_ave = box2_sum / real(count,d_prec)
              write (iout,'(A,6x,3f15.4,A)') ' Box parameters'
     &                ,box_ave,'  Angstrom'  
              stress_ave =  stress_sum / real(count,d_prec)
              stress2_ave = stress2_sum / real(count,d_prec)
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
