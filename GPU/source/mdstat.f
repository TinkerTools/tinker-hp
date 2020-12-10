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
#include "tinker_precision.h"
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
      use virial ,only:use_virial
      implicit none
      integer istep,modstep,ierr,i
      integer,parameter:: d_prec= kind(1.0d0)
      integer period,freq
      real(r_p) dt,temp,pres
      real(r_p) etot,epot,ekin
      real(d_prec) pico,dens
      real(d_prec) fluctuate,fluctuate2
      real(d_prec) intfluct,intfluct2
      real(d_prec) potfluct,potfluct2
      real(d_prec) kinfluct,kinfluct2
      real(d_prec) tfluct,pfluct,dfluct
      real(d_prec) tfluct2,pfluct2,dfluct2
      real(d_prec) etot_sum,etot2_sum
      real(d_prec) eint_sum,eint2_sum
      real(d_prec) etot_ave,etot2_ave
      real(d_prec) eint_ave,eint2_ave
      real(d_prec) epot_sum,epot2_sum
      real(d_prec) ekin_sum,ekin2_sum
      real(d_prec) epot_ave,epot2_ave
      real(d_prec) ekin_ave,ekin2_ave
      real(d_prec) temp_sum,temp2_sum
      real(d_prec) temp_ave,temp2_ave
      real(d_prec) pres_sum,pres2_sum
      real(d_prec) pres_ave,pres2_ave
      real(d_prec) dens_sum,dens2_sum
      real(d_prec) dens_ave,dens2_ave
      real*8 buffer(18)
      logical display
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
#if (defined(SINGLE) | defined(MIXED))
      freq    = 100
      period  = max(1,iprint/freq)
#else
      freq    = iprint
      period  = 1
#endif

      modstep = mod(istep,iprint)
      display = (mod(istep,period).eq.0.and.verbose)
c
c     zero out summation variables for new averaging period
c
      if (modstep.eq.1 .or. iprint.eq.1) then
         etot_sum  = 0.0_d_prec
         etot2_sum = 0.0_d_prec
         epot_sum  = 0.0_d_prec
         epot2_sum = 0.0_d_prec
         ekin_sum  = 0.0_d_prec
         ekin2_sum = 0.0_d_prec
         eint_sum  = 0.0_d_prec
         eint2_sum = 0.0_d_prec
         temp_sum  = 0.0_d_prec
         temp2_sum = 0.0_d_prec
         pres_sum  = 0.0_d_prec
         pres2_sum = 0.0_d_prec
         dens_sum  = 0.0_d_prec
         dens2_sum = 0.0_d_prec
      end if
c
c     print energy, temperature and pressure for current step
c
      if (rank.eq.0) then
        if (modstep.eq.0.or.display) then
!$acc data present(etot,epot,ekin,temp,pres) async
!$acc update host(etot,epot,ekin,temp,pres) async
!$acc end data
        end if
        if (verbose) then
           if (modstep .eq. 1) then
              if (use_bounds .and. integrate.ne.'STOCHASTIC') then
                 if (n>5d6) then
                  if  (use_virial) then
                  write (iout,11)
   11             format (/,4x,'MD Step',9x,'E Total',6x,'E Potential',
     &                      8x,'E Kinetic',7x,'Temp',7x,'Pres',/)
                  else
                  write (iout,13)
   13             format (/,4x,'MD Step',9x,'E Total',6x,'E Potential',
     &                      8x,'E Kinetic',7x,'Temp'/)
                  end  if
                 else
                  if (use_virial) then
                  write (iout,10)
   10             format (/,4x,'MD Step',6x,'E Total',3x,'E Potential',
     &                      5x,'E Kinetic',7x,'Temp',7x,'Pres',/)
                  else
                  write (iout,20)
   20             format (/,4x,'MD Step',6x,'E Total',3x,'E Potential',
     &                      5x,'E Kinetic',7x,'Temp',/)
                  end  if
                 end if
              end if
           end if
           if (display) then
           if  (use_bounds .and. integrate.ne.'STOCHASTIC') then
!$acc wait
              if (n>5d6) then
                 if (use_virial) then
                 write (iout,33)  istep,etot,epot,ekin,temp,pres
   33            format (i10,3f17.4,2f11.2)
                 else
                 write (iout,34)  istep,etot,epot,ekin,temp
   34            format (i10,3f17.4,f11.2)
                 end if
              else
                 if (use_virial) then
                 write (iout,30)  istep,etot,epot,ekin,temp,pres
                 !if(n.lt.5d4) write (iout,*) epot
#ifdef SINGLE
   30            format (i10,3f14.2,2f11.2)
#else
   30            format (i10,3f14.4,2f11.2)
#endif
                 else
                 write (iout,31)  istep,etot,epot,ekin,temp
#ifdef SINGLE
   31            format (i10,3f14.2,f11.2)
#else
   31            format (i10,3f14.4,f11.2)
#endif
                 end if
              end if
           else
!$acc wait
              write (iout,40)  istep,etot,epot,ekin,temp
   40         format (i10,3f14.4,f11.2)
           end  if
           end if
        end if
c
c       print  header for the averages over a group of recent steps
c
        if (verbose.and.modstep.eq.0) then
!$acc wait
           pico = real(istep,d_prec) * dt
           write (iout,50)  iprint,istep
   50      format (/,' Average Values for the last',i6,' out of',
     &                i9,' Dynamics Steps')
           write (iout,60)  pico
   60      format (/,' Simulation Time',5x,f15.4,' Picosecond')
        end if
c
c       compute total energy and fluctuation for recent steps
c
        if (display) then
           etot_sum = etot_sum + real(etot,d_prec)
           etot2_sum = etot2_sum + real(etot,d_prec)**2
        end if
        if (verbose.and.modstep.eq.0) then
           etot_ave = etot_sum / real(freq,d_prec)
           etot2_ave = etot2_sum / real(freq,d_prec)
           fluctuate2 = etot2_ave - etot_ave**2
           if (fluctuate2 .gt. 0.0_ti_p) then
              fluctuate = sqrt(fluctuate2)
           else
              fluctuate = 0.0_ti_p
           end if
           write (iout,70)  etot_ave,fluctuate
   70      format (' Total Energy',8x,f15.4,' Kcal/mole',3x,
     &                '(+/-',f10.4,')')
        end if
c
c       compute average potential energy and its fluctuation
c
        if (display) then
        epot_sum = epot_sum + real(epot,d_prec)
        epot2_sum = epot2_sum + real(epot,d_prec)**2
        end if
        if (verbose.and.modstep.eq.0) then
           epot_ave = epot_sum / real(freq,d_prec)
           epot2_ave = epot2_sum / real(freq,d_prec)
           potfluct2 = epot2_ave - epot_ave**2
           if (potfluct2 .gt. 0.0_ti_p) then
              potfluct = sqrt(potfluct2)
           else
              potfluct = 0.0_ti_p
           end if
           write (iout,80)  epot_ave,potfluct
   80      format (' Potential Energy',4x,f15.4,' Kcal/mole',3x,
     &                '(+/-',f10.4,')')
        end if
c
c       compute average kinetic energy and its fluctuation
c
        if (display) then
        ekin_sum = ekin_sum + real(ekin,d_prec)
        ekin2_sum = ekin2_sum + real(ekin,d_prec)**2
        end if
        if (verbose.and.modstep.eq.0) then
           ekin_ave = ekin_sum / real(freq,d_prec)
           ekin2_ave = ekin2_sum / real(freq,d_prec)
           kinfluct2 = ekin2_ave - ekin_ave**2
           if (kinfluct2 .gt. 0.0_ti_p) then
              kinfluct = sqrt(kinfluct2)
           else
              kinfluct = 0.0_ti_p
           end if
           write (iout,90)  ekin_ave,kinfluct
   90      format (' Kinetic Energy',6x,f15.4,' Kcal/mole',3x,
     &                '(+/-',f10.4,')')
        end if
c
c       compute average intermolecular energy and its fluctuation
c
        if (nmol.ne.1 .and. nmol.ne.n .and. .not.use_ewald) then
           if (display) then
           eint_sum = eint_sum + real(einter,d_prec)
           eint2_sum = eint2_sum + real(einter,d_prec)**2
           end if
           if (verbose.and.modstep.eq.0) then
              eint_ave = eint_sum / real(freq,d_prec)
              eint2_ave = eint2_sum / real(freq,d_prec)
              intfluct2 = eint2_ave - eint_ave**2
              if (intfluct2 .gt. 0.0_ti_p) then
                 intfluct = sqrt(intfluct2)
              else
                 intfluct = 0.0_ti_p
              end if
              write (iout,100)  eint_ave,intfluct
  100         format (' Intermolecular',6x,f15.4,' Kcal/mole',3x,
     &                   '(+/-',f10.4,')')
           end if
        end if
c
c       compute the average temperature and its fluctuation
c
        if (display) then
        temp_sum = temp_sum + real(temp,d_prec)
        temp2_sum = temp2_sum + real(temp,d_prec)**2
        end if
        if (verbose.and.modstep.eq.0) then
           temp_ave = temp_sum / real(freq,d_prec)
           temp2_ave = temp2_sum / real(freq,d_prec)
           tfluct2 = temp2_ave - temp_ave**2
           if (tfluct2 .gt. 0.0_ti_p) then
              tfluct = sqrt(tfluct2)
           else
              tfluct = 0.0_ti_p
           end if
           write (iout,110)  temp_ave,tfluct
  110      format (' Temperature',9x,f15.2,' Kelvin',6x,
     &                '(+/-',f10.2,')')
        end if
c
c       compute the average pressure and its fluctuation
c
        if (use_bounds) then
        if (display) then
           pres_sum = pres_sum + real(pres,d_prec)
           pres2_sum = pres2_sum + real(pres,d_prec)**2
        end if
           if (verbose.and.modstep.eq.0) then
              pres_ave = pres_sum / real(freq,d_prec)
              pres2_ave = pres2_sum / real(freq,d_prec)
              pfluct2 = pres2_ave - pres_ave**2
              if (pfluct2 .gt. 0.0_ti_p) then
                 pfluct = sqrt(pfluct2)
              else
                 pfluct = 0.0_ti_p
              end if
              if (use_virial) then
              write (iout,120)  pres_ave,pfluct
  120         format (' Pressure',12x,f15.2,' Atmosphere',2x,
     &                   '(+/-',f10.2,')')
              end if
           end if
c
c       compute the average density and its fluctuation
c
        if (display) then
           dens = (1.0d24/real(volbox,d_prec))
     &          * (real(totmass,d_prec)/real(avogadro,d_prec))
           dens_sum = dens_sum + dens
           dens2_sum = dens2_sum + dens**2
        end if
           if (verbose.and.modstep.eq.0) then
              dens_ave = dens_sum / real(freq,d_prec)
              dens2_ave = dens2_sum / real(freq,d_prec)
              dfluct2 = dens2_ave - dens_ave**2
              if (dfluct2 .gt. 0.0_ti_p) then
                 dfluct = sqrt(dfluct2)
              else
                 dfluct = 0.0_ti_p
              end if
              write (iout,130)  dens_ave,dfluct
  130         format (' Density',13x,f15.4,' Grams/cc',4x,
     &                   '(+/-',f10.4,')')
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
      if (modstep.eq.0.and.btest(tinkertime,sumy_time)) then
         call timer_exit( timer_timestep,quiet_timers )
         call timer_enter( timer_timestep )
         call display_timers( stat_timers,config=ave_disp,slot=1,
     &        iter=iprint )
         call timer_save( stat_timers,slot=1 )
      end if
      end
