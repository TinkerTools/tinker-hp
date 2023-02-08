c
c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine temper  --  thermostat applied at half step  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "temper" applies a velocity correction at the half time step
c     as needed for the Nose-Hoover extended system thermostat
c
c     literature references:
c
c     D. Frenkel and B. Smit, "Understanding Molecular Simulation,
c     2nd Edition", Academic Press, San Diego, CA, 2002; see Appendix
c     E.2 for implementation details
c
c     G. J. Martyna, M. E. Tuckerman, D. J. Tobias and M. L. Klein,
c     "Explicit Reversible Integrators for Extended Systems Dynamics",
c     Molecular Physics, 87, 1117-1157 (1996)
c
c
#include "tinker_macro.h"
      subroutine temper (dt,eksum,ekin,temp)
      use atoms
      use atmlst
      use atmtyp
      use bath
      use domdec
      use energi  ,only: calc_e
      use group
      use mdstuf
      use molcul
      use moldyn
      use mpi
      use random_mod
      use units
      use usage
      use virial  ,only: use_virial
      use timestat
      implicit none
      integer i,j,nc,ns,iglob
      integer k,m,ierr
      real(r_p) dt,dtc,dts
      real(r_p) dt2,dt4,dt8
      real(r_p) eksum
      real(r_p) temp,ekt
      real(r_p),save:: scale
      real(r_p) expterm,speed
      real(r_p) w(3)
      real(r_p) ekin(3,3)
      real(r_p),save:: c,d
      real(t_p) s,si
      real(r_p) kt,rate,trial
      logical,save::f_in=.true.
c
      call timer_enter( timer_other )
      if (.not. isothermal)  goto 10
      if (f_in) then
         f_in=.false.
         pick1=0
!$acc enter data copyin(pick1,c,d,pick,scale)
      end if
c
      if (thermostat .eq. 'BERENDSEN' .or.
     &    thermostat .eq. 'BUSSI'    ) then
         call kineticgpu (eksum,ekin,temp)
      end if
c
c     couple to external temperature bath via Berendsen scaling
c
      if (thermostat .eq. 'BERENDSEN') then
!$acc serial async present(scale,temp)
         scale = 1.0_re_p
         if (temp .ne. 0.0_re_p)
     &      scale = sqrt(1.0_re_p + (dt/tautemp)*(kelvin/temp-1.0_re_p))
!$acc end serial
!$acc parallel loop collapse(2) async present(scale,glob,use,v)
            do i = 1, nloc; do j = 1, 3
               iglob = glob(i)
               if (use(iglob)) then
                  v(j,iglob) = scale * v(j,iglob)
               end if
            end do; end do
c
c     couple to external temperature bath via Bussi scaling
c
      else if (thermostat .eq. 'BUSSI') then
c
#ifdef _OPENACC
         if (.not.host_rand_platform) then
            call normalgpu(samplevec(1),nfree)
!$acc parallel loop async present(pick1,samplevec)
            do i = 2, nfree
               pick1 = pick1 + real(samplevec(i)**2,r_p)
            end do
         end if
#endif
         if (host_rand_platform) then
            pick1 = 0.0_re_p
            pick = normal ()
            do i = 1, nfree-1
               si = normal ()
               pick1  = pick1 + real(si*si,r_p)
            end do
!$acc update device(pick1,pick) async
         end if
c
!$acc serial async present(pick,pick1,c,d,samplevec,temp,scale,eta)
#ifdef _OPENACC
         if (.not.host_rand_platform) pick = samplevec(1)
#endif
         if (temp .eq. 0.0_re_p)  temp = 0.1_re_p
         c = exp(-dt/tautemp)
         d = (1.0_re_p-c) * (kelvin/temp) / real(nfree,r_p)
         scale = c + (pick1+pick*pick)*d + 2.0_re_p*pick*sqrt(c*d)
         scale = sqrt(scale)
         if (pick+sqrt(c/d) .lt. 0.0_re_p)  scale = -scale
         pick1 = 0.0_re_p
!$acc end serial
         if (nproc.ne.1) then
!$acc wait
!$acc host_data use_device(scale)
            call MPI_BCAST(scale,1,MPI_RPREC,0,COMM_TINKER,ierr)
!$acc end host_data
         end if
c        if(rank.eq.0) print*,scale,pick
         if (barostat.eq.'BUSSI') then
!$acc serial async present(eta,scale)
         eta = eta * scale
!$acc end serial
         end if
c
!$acc parallel loop collapse(2) async
!$acc&         present(glob,use,v,scale)
         do i = 1, nloc; do j = 1, 3
            iglob = glob(i)
            if (use(iglob)) then
               v(j,iglob) = scale * v(j,iglob)
            end if
         end do; end do
c
c     select random velocities via Andersen stochastic collisions
c
      else if (thermostat .eq. 'ANDERSEN') then
#ifdef _OPENACC
 34      format('ANDERSEN Thermostat is unavailable on device platform'
     &       ,/,3x,'Use host application to benefit from it')
         write(0,34)
         __TINKER_FATAL__
#endif
         kt = boltzmann * kelvin
         rate = 1000.0_re_p * dt * collide
         if (barostat.eq.'MONTECARLO' .and.
     &            volscale.eq.'MOLECULAR') then
            call molecule(.false.)
            rate = rate / real(nmol,t_p)**(2.0_re_p/3.0_re_p)
            do i = 1, nmoleloc
               iglob = molculeglob(iglob)
               trial = random ()
               if (trial .lt. rate) then
                  do j = imol(1,iglob), imol(2,iglob)
                     k = kmol(j)
                     speed = sqrt(kt/mass(k))
                     do m = 1, 3
                        v(m,k) = speed * normal ()
                     end do
                  end do
               end if
            end do
         else
            rate = rate / real(nuse,r_p)**(2.0_re_p/3.0_re_p)
            do i = 1, nloc
               iglob = glob(i)
               if (use(iglob)) then
                  trial = random ()
                  if (trial .lt. rate) then
                     speed = sqrt(kt/mass(iglob))
                     do j = 1, 3
                        v(j,iglob) = speed * normal ()
                     end do
                  end if
               end if
            end do
         end if
      end if
c
c     recompute kinetic energy and instantaneous temperature
c
  10  continue
      if (calc_e.or.use_virial) call kineticgpu ( eksum,ekin,temp )
      call timer_exit ( timer_other,quiet_timers )

      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine temper2  --  thermostat applied at full step  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "temper2" computes the instantaneous temperature and applies a
c     thermostat via Berendsen or Bussi-Parrinello velocity scaling,
c     Andersen stochastic collisions or Nose-Hoover extended system
c
c     literature references:
c
c     H. J. C. Berendsen, J. P. M. Postma, W. F. van Gunsteren,
c     A. DiNola and J. R. Hauk, "Molecular Dynamics with Coupling
c     to an External Bath", Journal of Chemical Physics, 81,
c     3684-3690 (1984)
c
c     G. Bussi and M. Parrinello, "Stochastic Thermostats: Comparison
c     of Local and Global Schemes", Computer Physics Communications,
c     179, 26-29 (2008)
c
c     H. C. Andersen, "Molecular Dynamics Simulations at Constant
c     Pressure and/or Temperature", Journal of Chemical Physics,
c     72, 2384-2393 (1980)
c
c
      subroutine temper2 (temp)
      use atmtyp
      use bath
      use domdec
      use group
      use mdstuf
      use moldyn
      use units
      use usage
      use timestat
      implicit none

      real(r_p),save:: eksum
      real(r_p) temp
      real(r_p),save:: ekin(3,3)
      logical  ,save::f_in=.true.

      if (f_in) then
         f_in=.false.
!$acc enter data create(ekin,eksum)
      end if
c
c     get instantaneous temperature from the kinetic energy
c
      call timer_enter( timer_other )
      call kineticgpu (eksum,ekin,temp)

      call timer_exit( timer_other,quiet_timers )
      if (.not. isothermal)  return
      return
cc
cc     couple to external temperature bath via Berendsen scaling
cc
c      if (thermostat .eq. 'BERENDSEN') then
c         call kinetic (eksum,ekin)
c         temp = 2.0_re_p * eksum / (real(nfree,t_p) * gasconst)
c         if (temp .eq. 0.0_re_p)  temp = 0.1_re_p
c         scale = sqrt(1.0_re_p + (dt/tautemp)*(kelvin/temp-1.0_re_p))
c         do i = 1, nloc
c            iglob = glob(i)
c            if (use(iglob)) then
c               do j = 1, 3
c                  v(j,iglob) = scale * v(j,iglob)
c               end do
c            end if
c         end do
cc
cc     couple to external temperature bath via Bussi scaling
cc
c      else if (thermostat .eq. 'BUSSI') then
c         call kinetic (eksum,ekin)
c         temp = 2.0_re_p * eksum / (real(nfree,t_p) * gasconst)
c         if (temp .eq. 0.0_re_p)  temp = 0.1_re_p
c         c = exp(-dt/tautemp)
c         d = (1.0_re_p-c) * (kelvin/temp) / real(nfree,t_p)
c         r = normal ()
c         s = 0.0_re_p
c         do i = 1, nfree-1
c            si = normal ()
c            s = s + si*si
c         end do
c         scale = c + (s+r*r)*d + 2.0_re_p*r*sqrt(c*d)
c         scale = sqrt(scale)
c         if (r+sqrt(c/d) .lt. 0.0_re_p)  scale = -scale
c         eta = eta * scale
c         do i = 1, nloc
c            iglob = glob(i)
c            if (use(iglob)) then
c               do j = 1, 3
c                  v(j,iglob) = scale * v(j,iglob)
c               end do
c            end if
c         end do
cc
cc     select random velocities via Andersen stochastic collisions
cc
c      else if (thermostat .eq. 'ANDERSEN') then
c         kt = boltzmann * kelvin
c         rate = 1000.0_re_p * dt * collide
c         rate = rate / real(nuse,t_p)**(2.0_re_p/3.0_re_p)
c         do i = 1, nloc
c            iglob = glob(i)
c            if (use(iglob)) then
c               trial = random ()
c               if (trial .lt. rate) then
c                  speed = sqrt(kt/mass(iglob))
c                  do j = 1, 3
c                     v(j,iglob) = speed * normal ()
c                  end do
c               end if
c            end if
c         end do
cc
cc     make full-step velocity correction for Nose-Hoover system
cc
c      else if (thermostat .eq. 'NOSE-HOOVER') then
c         ekt = gasconst * kelvin
c         nc = 5
c         ns = 3
c         dtc = dt / real(nc,t_p)
c         w(1) = 1.0_re_p / (2.0_re_p-2.0_re_p**(1.0_re_p/3.0_re_p))
c         w(2) = 1.0_re_p - 2.0_re_p*w(1)
c         w(3) = w(1)
c         scale = 1.0_re_p
c         do i = 1, nc
c            do j = 1, ns
c               dts = w(j) * dtc
c               dt2 = 0.5_re_p * dts
c               dt4 = 0.25_re_p * dts
c               dt8 = 0.125_re_p * dts
c               gnh(4) = (qnh(3)*vnh(3)*vnh(3)-ekt) / qnh(4)
c               vnh(4) = vnh(4) + gnh(4)*dt4
c               gnh(3) = (qnh(2)*vnh(2)*vnh(2)-ekt) / qnh(3)
c               expterm = exp(-vnh(4)*dt8)
c               vnh(3) = expterm * (vnh(3)*expterm+gnh(3)*dt4)
c               gnh(2) = (qnh(1)*vnh(1)*vnh(1)-ekt) / qnh(2)
c               expterm = exp(-vnh(3)*dt8)
c               vnh(2) = expterm * (vnh(2)*expterm+gnh(2)*dt4)
c               gnh(1) = (2.0_re_p*eksum-real(nfree,t_p)*ekt) / qnh(1)
c               expterm = exp(-vnh(2)*dt8)
c               vnh(1) = expterm * (vnh(1)*expterm+gnh(1)*dt4)
c               scale = scale * exp(-vnh(1)*dt2)
c               eksum = eksum * scale * scale
c               gnh(1) = (2.0_re_p*eksum-real(nfree,t_p)*ekt) / qnh(1)
c               expterm = exp(-vnh(2)*dt8)
c               vnh(1) = expterm * (vnh(1)*expterm+gnh(1)*dt4)
c               gnh(2) = (qnh(1)*vnh(1)*vnh(1)-ekt) / qnh(2)
c               expterm = exp(-vnh(3)*dt8)
c               vnh(2) = expterm * (vnh(2)*expterm+gnh(2)*dt4)
c               gnh(3) = (qnh(2)*vnh(2)*vnh(2)-ekt) / qnh(3)
c               expterm = exp(-vnh(4)*dt8)
c               vnh(3) = expterm * (vnh(3)*expterm+gnh(3)*dt4)
c               gnh(4) = (qnh(3)*vnh(3)*vnh(3)-ekt) / qnh(4)
c               vnh(4) = vnh(4) + gnh(4)*dt4
c            end do
c         end do
c         do i = 1, nloc
c            iglob = glob(i)
c            if (use(iglob)) then
c               do j = 1, 3
c                  v(j,iglob) = scale * v(j,iglob)
c               end do
c            end if
c         end do
c      end if
cc
cc     recompute kinetic energy and instantaneous temperature
cc
c      call kinetic (eksum,ekin)
c      temp = 2.0_re_p * eksum / (real(nfree,t_p) * gasconst)
      end
