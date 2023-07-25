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
#include "tinker_macro.h"

      ! Determine wether or not to compute energy
      subroutine inquire_calc_e(istep,period)
      use bath   ,only: isobaric,barostat
      use energi ,only: calc_e
      use inform ,only: iprint,verbose,deb_Path,deb_Energy
      use plumed ,only: lplumed
      implicit none
      integer,intent(in):: istep,period

      ! Compute energy for next step ?
      calc_e =   (mod(istep,iprint).eq.0)
     &       .or.(mod(istep,period).eq.0.and.verbose)
     &       .or.(isobaric.and.barostat.eq.'MONTECARLO')
     &       .or.lplumed

      if (deb_Path.or.(deb_Energy.and.calc_e))
     &    print*, 'inquire compute energy',calc_e
      end subroutine

      subroutine emergency_save
      use argue    ,only:arg
      use dcdmod   ,only:dcdio
      use inform   ,only:deb_Path,deb_Force,abort,n_fwriten
      use moldyn   ,only:stepfast,stepint,nalt,step_c
      use mdstuf1  ,only:epot
      use output   ,only:archive,new_restart,f_mdsave
      implicit none
      logical save_dcdio,save_arc
      real(r_p) dt
      read(arg(3),*) dt
      dt = 1d-3*dt
      new_restart = .true.
      f_mdsave    = .true.
      save_dcdio  = dcdio
      save_arc    = archive
      dcdio       = .false.
      archive     = .false.
      n_fwriten   = - 1
      call mdsave(step_c,dt,epot)
      new_restart = .false.
      f_mdsave    = .false.
      dcdio       = save_dcdio
      archive     = save_arc
      end subroutine

      subroutine mdstat (istep,dt,etot,epot,ekin,temp,pres)
      use sizes
      use atoms
      use bound
      use boxes
      use bath
      use cutoff
      use domdec
      use energi ,only: calc_e,epot_mean,epot_std10
     &           , etot_ave,etot_std
      use inform
      use inter
      use iounit
      use mdstuf
      use mdstuf1
      use molcul
      use potent
      use mdstate
      use timestat
      use units
      use mpi
      use virial ,only:use_virial
      implicit none
      integer,parameter:: d_prec= kind(1.0d0)
      integer istep,modstep
      integer period,freq
      integer, save :: count=0
      real(r_p) dt,temp,pres
      real(r_p) etot,epot,ekin
      real(d_prec) pico,dens
      real(d_prec) fluctuate2
      real(d_prec) intfluct,intfluct2
      real(d_prec) potfluct,potfluct2
      real(d_prec) kinfluct,kinfluct2
      real(d_prec) tfluct,pfluct,dfluct,vfluct
      real(d_prec) tfluct2,pfluct2,dfluct2,vfluct2
      real(d_prec) tpistonfluct
      real(d_prec) tpistonfluct2
      real(d_prec) etot_sum,etot2_sum
      real(d_prec) eint_sum,eint2_sum
      real(d_prec) etot2_ave
      real(d_prec) eint_ave,eint2_ave
      real(d_prec) epot_sum,epot2_sum
      real(d_prec) ekin_sum,ekin2_sum
      real(d_prec) epot_ave,epot2_ave
      real(d_prec) ekin_ave,ekin2_ave
      real(d_prec) temp_sum,temp2_sum
      real(d_prec) temp_ave,temp2_ave
      real(d_prec) tpiston_sum,tpiston2_sum
      real(d_prec) tpiston_ave,tpiston2_ave
      real(d_prec) pres_sum,pres2_sum
      real(d_prec) pres_ave,pres2_ave
      real(d_prec) dens_sum,dens2_sum
      real(d_prec) dens_ave,dens2_ave
      real(d_prec)  vol_sum, vol2_sum
      real(d_prec)  vol_ave, vol2_ave
      real(d_prec)  box_sum(3), box2_sum(3)
      real(d_prec)  box_ave(3), box2_ave(3)
      real(d_prec)  stress_sum(3,3), stress2_sum(3,3)
      real(d_prec)  stress_ave(3,3), stress2_ave(3,3)
      real(d_prec) buffer(18)
      logical display
      save etot_sum,etot2_sum
      save eint_sum,eint2_sum
      save epot_sum,epot2_sum
      save ekin_sum,ekin2_sum
      save temp_sum,temp2_sum
      save tpiston_sum,tpiston2_sum
      save pres_sum,pres2_sum
      save dens_sum,dens2_sum
      save  vol_sum, vol2_sum
      save  box_sum, box2_sum
      save  stress_sum, stress2_sum
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
      if (iprint.eq.1) freq=1
#else
      freq    = iprint
      period  = 1
#endif

      if (track_mds.and.mod(istep,ms_back_p).eq.0) then
         call mds_save
      end if

      modstep =  mod(istep,iprint)
      display = (mod(istep,period).eq.0.and.verbose)

c
c     zero out summation variables for new averaging period
c
      if (modstep.eq.1 .or. iprint.eq.1) then
         count=0
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
         tpiston_sum  = 0.0_d_prec
         tpiston2_sum = 0.0_d_prec
         pres_sum  = 0.0_d_prec
         pres2_sum = 0.0_d_prec
         dens_sum  = 0.0_d_prec
         dens2_sum = 0.0_d_prec
         vol_sum   = 0.0_d_prec
         vol2_sum  = 0.0_d_prec
         box_sum   = 0.0_d_prec
         box2_sum  = 0.0_d_prec
         stress_sum   = 0.0_d_prec
         stress2_sum  = 0.0_d_prec
      end if
      if(display) then
        count=count+1
        dens = (1.0d24/real(volbox,d_prec))
     &          * (real(totmass,d_prec)/real(avogadro,d_prec))
      endif
c
c     print energy, temperature and pressure for current step
c
      if (rank.eq.0) then
        if (modstep.eq.0.or.display) then
!$acc update host(etot,epot,ekin,temp,pres,stress) async
           if (.not.calc_e) then
   09         format(
     &        " -- ERROR -- mdstat ",/
     &       ,"    Energy computation is disabled !!"
     &       ," istep-",I0,1X,"period-",I0,2X,"verbose-",L1)
              write(0,09) istep,period,verbose
              call fatal
           end if
        end if
        if (verbose) then
           if (modstep .eq. 1) then
              if (use_bounds .and. integrate.ne.'STOCHASTIC') then
                 write(iout,'(A)',advance="no") 
     &                         "     MD Step"
     &                       //"     E total"
     &                       //"   E Potential"            
     &                       //"     E Kinetic"
     &                       //"       Temp"
                 if(use_virial) then
                   write(iout,'(A)',advance="no") 
     &                         "       Pres"
                 endif
c                 if(isobaric) then
c                   write(iout,'(A)',advance="no") 
c     &                         "     Density"
c     &                       //"      Volume"
c                 endif
                 write(iout,*)
              end if
           end if
           if (display) then
           if  (use_bounds .and. integrate.ne.'STOCHASTIC') then
!$acc wait
               write(iout,'(i10,3'//f_ener//',f11.2)',advance="no")
     &             istep,etot,epot,ekin,temp
               if(use_virial) then
                 write(iout,'(f11.2)',advance="no") pres
               endif
               !if(isobaric) then
               !  write(iout,'(f12.4,f12.2)',advance="no") dens,volbox
               !endif
               write(iout,*)
           else
!$acc wait
              write (iout,'(i10,3f14.4,f11.2)')  
     &             istep,etot,epot,ekin,temp
           end if
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
        if (display.and.etot_ave.eq.0.0) then
           etot_ave = etot
           etot_std = abs(etot)
        end if
c
c       compute total energy and fluctuation for recent steps
c
        if (display) then
           etot_sum  = etot_sum + real(etot,d_prec)
           etot2_sum = etot2_sum + real(etot,d_prec)**2
           if (.not.isothermal.and..not.isobaric) then
              if (abs(etot-etot_ave).gt.9*etot_std) then
  600            format(/,"Detected brutal shift in Total Energy"
     &                   ," (Etot|mean|std)",/,3F14.4)
                 write(0,600) etot,etot_ave,etot_std
                 __TINKER_FATAL__
              end if
           end if
        end if
        if (verbose.and.modstep.eq.0) then
           etot_ave = etot_sum / real(count,d_prec)
           etot2_ave = etot2_sum / real(count,d_prec)
           fluctuate2 = etot2_ave - etot_ave**2
           if (fluctuate2 .gt. 0.0_ti_p) then
              etot_std = sqrt(fluctuate2)
           else
              etot_std = 0.0_ti_p
           end if
           write (iout,70)  etot_ave,etot_std
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
           epot_ave = epot_sum / real(count,d_prec)
           epot2_ave = epot2_sum / real(count,d_prec)
           potfluct2 = epot2_ave - epot_ave**2
           if (potfluct2 .gt. 0.0_ti_p) then
              potfluct = sqrt(potfluct2)
           else
              potfluct = 0.0_ti_p
           end if
           epot_mean  =    epot_ave
           epot_std10 = 10*potfluct
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
           ekin_ave = ekin_sum / real(count,d_prec)
           ekin2_ave = ekin2_sum / real(count,d_prec)
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
              eint_ave = eint_sum / real(count,d_prec)
              eint2_ave = eint2_sum / real(count,d_prec)
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
           temp_ave = temp_sum / real(count,d_prec)
           temp2_ave = temp2_sum / real(count,d_prec)
           tfluct2 = temp2_ave - temp_ave**2
           if (tfluct2 .gt. 0.0_ti_p) then
              tfluct = sqrt(tfluct2)
           else
              tfluct = 0.0_ti_p
           end if
           write (iout,110)  temp_ave,tfluct
  110      format (' Temperature',7x,f15.2,'   Kelvin',6x,
     &                '(+/-',f8.2,'  )')
        end if

        if(use_piston .or. use_piston_save) then
          if (display) then
          tpiston_sum = tpiston_sum + (temppiston-kelvin)
          tpiston2_sum = tpiston2_sum + (temppiston -kelvin)**2
          end if
          if (verbose.and.modstep.eq.0) then
             tpistonfluct2 = (tpiston2_sum  
     &          - tpiston_sum**2/real(count,d_prec))
     &              /(real(count-1,d_prec))
             tpiston_ave = kelvin + tpiston_sum / real(count,d_prec)
             if (tpistonfluct2 .gt. 0.0_ti_p) then
                tpistonfluct = sqrt(tpistonfluct2)
             else
                tpistonfluct = 0.0_ti_p
             end if
             write (iout,'(A,12x,f15.2,A,6x,A,f10.2,A)')  ' T_piston'
     &            ,tpiston_ave,' Kelvin' ,'(+/-',tpistonfluct,')'
          end if
        endif
c
c       compute the average pressure and its fluctuation
c
        if (use_bounds) then
        if (display) then
           pres_sum = pres_sum + real(pres,d_prec)
           pres2_sum = pres2_sum + real(pres,d_prec)**2
        end if
           if (verbose.and.modstep.eq.0) then
              pres_ave = pres_sum / real(count,d_prec)
              pres2_ave = pres2_sum / real(count,d_prec)
              pfluct2 = pres2_ave - pres_ave**2
              if (pfluct2 .gt. 0.0_ti_p) then
                 pfluct = sqrt(pfluct2)
              else
                 pfluct = 0.0_ti_p
              end if
              if (use_virial) then
              write (iout,120)  pres_ave,pfluct
  120         format (' Pressure',10x,f15.2,'   Atmosphere',2x,
     &                   '(+/-',f8.2,'  )')
              end if
           end if
c
c       compute the average density and its fluctuation
c
        if (display) then
           dens_sum = dens_sum + dens
           dens2_sum = dens2_sum + dens**2
        end if
           if (verbose.and.modstep.eq.0) then
              dens_ave = dens_sum / real(count,d_prec)
              dens2_ave = dens2_sum / real(count,d_prec)
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
c
c       compute the average volume and its fluctuation
c
        if (display.and.(isobaric.or.isobaric_save)) then
            vol_sum =  vol_sum + volbox*1d-3
           vol2_sum = vol2_sum +(volbox*1d-3)**2
           if (anisotrop) then
             stress_sum = stress_sum + stress
             stress2_sum = stress2_sum + stress**2
             box_sum(1) = box_sum(1) + xbox
             box_sum(2) = box_sum(2) + ybox
             box_sum(3) = box_sum(3) + zbox
             box2_sum(1) = box2_sum(1) + xbox**2
             box2_sum(2) = box2_sum(2) + ybox**2
             box2_sum(3) = box2_sum(3) + zbox**2
           endif
        end if
           if (verbose.and.(isobaric.or.isobaric_save)
     &         .and.modstep.eq.0) then
              vol_ave =  vol_sum / real(freq,d_prec)
              vol2_ave = vol2_sum / real(freq,d_prec)
              vfluct2  = vol2_ave - vol_ave**2
              if (vfluct2 .gt. 0.0_ti_p) then
                 vfluct = sqrt(vfluct2)
              else
                 vfluct = 0.0_ti_p
              end if
              write (iout,140)  vol_ave,vfluct
  140         format (' Volume',14x,f15.4,' nm^3',8x,
     &                   '(+/-',f10.4,')')
              if (anisotrop) then
                box_ave =  box_sum / real(freq,d_prec)
                box2_ave = box2_sum / real(freq,d_prec)
                write (iout,'(A,6x,3f15.4,A)') ' Box parameters'
     &                  ,box_ave,'  Angstrom'  
                stress_ave =  stress_sum / real(freq,d_prec)
                stress2_ave = stress2_sum / real(freq,d_prec)
                write (iout,'(A,A,2x,6f15.4,A)') ' Stress Tensor '
     &                  ,'(xx yy zz xy yz xz)'
     &                  ,stress_ave(1,1),stress_ave(2,2),stress_ave(3,3)
     &                  ,stress_ave(1,2),stress_ave(2,3),stress_ave(1,3)
     &                  ,' Atmosphere'
              endif
           end if
           if (verbose.and.isobaric.and.mtc_nacc.ne.0.and.modstep.eq.0)
     &        then
              write (iout,150) mtc_nacc
  150         format(' Montecarlo Barostat applied',3x,I5,' times')
              mtc_nacc = 0
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

      call inquire_calc_e(istep+1,period)

      end
