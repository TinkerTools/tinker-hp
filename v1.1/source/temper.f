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
      subroutine temper (dt,eksum,ekin,temp)
      use atmlst
      use atmtyp
      use bath
      use domdec
      use group
      use mdstuf
      use molcul
      use moldyn
      use units
      use usage
      use mpi
      implicit none
      integer i,j,nc,ns,iglob
      integer k,m
      real*8 dt,dtc,dts
      real*8 dt2,dt4,dt8
      real*8 eksum,ekt
      real*8 temp
      real*8 scale,expterm,speed
      real*8 w(3)
      real*8 ekin(3,3)
      real*8 c,d,r,s,si
      real*8 random,normal
      real*8 kt,rate,trial
c
      call kinetic (eksum,ekin,temp)
      if (.not. isothermal)  return
c
c     couple to external temperature bath via Berendsen scaling
c
      if (thermostat .eq. 'BERENDSEN') then
         scale = 1.0d0
         if (temp .ne. 0.0d0)
     &      scale = sqrt(1.0d0 + (dt/tautemp)*(kelvin/temp-1.0d0))
            do i = 1, nloc
               iglob = glob(i)
               if (use(iglob)) then
                  do j = 1, 3
                     v(j,iglob) = scale * v(j,iglob)
                  end do
               end if
            end do
c
c     couple to external temperature bath via Bussi scaling
c
      else if (thermostat .eq. 'BUSSI') then
         if (temp .eq. 0.0d0)  temp = 0.1d0
         c = exp(-dt/tautemp)
         d = (1.0d0-c) * (kelvin/temp) / dble(nfree)
         r = normal ()
         s = 0.0d0
         do i = 1, nfree-1
            si = normal ()
            s = s + si*si
         end do
         scale = c + (s+r*r)*d + 2.0d0*r*sqrt(c*d)
         scale = sqrt(scale)
         if (r+sqrt(c/d) .lt. 0.0d0)  scale = -scale
         eta = eta * scale
         do i = 1, nloc
            iglob = glob(i)
            if (use(iglob)) then
               do j = 1, 3
                  v(j,iglob) = scale * v(j,iglob)
               end do
            end if
         end do
c
c     select random velocities via Andersen stochastic collisions
c
      else if (thermostat .eq. 'ANDERSEN') then
         kt = boltzmann * kelvin
         rate = 1000.0d0 * dt * collide
         if (barostat.eq.'MONTECARLO' .and.
     &            volscale.eq.'MOLECULAR') then
            call molecule(.false.)
            rate = rate / dble(nmol)**(2.0d0/3.0d0)
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
            rate = rate / dble(nuse)**(2.0d0/3.0d0)
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
      call kinetic (eksum,ekin,temp)
      temp = 2.0d0 * eksum / (dble(nfree) * gasconst)
c
      return
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
      subroutine temper2 (dt,temp)
      use atmtyp
      use bath
      use domdec
      use group
      use mdstuf
      use moldyn
      use units
      use usage
      implicit none
      integer i,j,k,m,iglob
      integer nc,ns
      real*8 dt,dtc,dts
      real*8 dt2,dt4,dt8
      real*8 eksum,ekt
      real*8 scale,speed
      real*8 c,d,r,s,si
      real*8 random,normal
      real*8 kt,rate,trial
      real*8 temp,expterm
      real*8 w(3)
      real*8 ekin(3,3)
c
c
c     get instantaneous temperature from the kinetic energy
c
      call kinetic (eksum,ekin,temp)
      if (.not. isothermal)  return
      return
cc
cc     couple to external temperature bath via Berendsen scaling
cc
c      if (thermostat .eq. 'BERENDSEN') then
c         call kinetic (eksum,ekin)
c         temp = 2.0d0 * eksum / (dble(nfree) * gasconst)
c         if (temp .eq. 0.0d0)  temp = 0.1d0
c         scale = sqrt(1.0d0 + (dt/tautemp)*(kelvin/temp-1.0d0))
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
c         temp = 2.0d0 * eksum / (dble(nfree) * gasconst)
c         if (temp .eq. 0.0d0)  temp = 0.1d0
c         c = exp(-dt/tautemp)
c         d = (1.0d0-c) * (kelvin/temp) / dble(nfree)
c         r = normal ()
c         s = 0.0d0
c         do i = 1, nfree-1
c            si = normal ()
c            s = s + si*si
c         end do
c         scale = c + (s+r*r)*d + 2.0d0*r*sqrt(c*d)
c         scale = sqrt(scale)
c         if (r+sqrt(c/d) .lt. 0.0d0)  scale = -scale
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
c         rate = 1000.0d0 * dt * collide
c         rate = rate / dble(nuse)**(2.0d0/3.0d0)
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
c         dtc = dt / dble(nc)
c         w(1) = 1.0d0 / (2.0d0-2.0d0**(1.0d0/3.0d0))
c         w(2) = 1.0d0 - 2.0d0*w(1)
c         w(3) = w(1)
c         scale = 1.0d0
c         do i = 1, nc
c            do j = 1, ns
c               dts = w(j) * dtc
c               dt2 = 0.5d0 * dts
c               dt4 = 0.25d0 * dts
c               dt8 = 0.125d0 * dts
c               gnh(4) = (qnh(3)*vnh(3)*vnh(3)-ekt) / qnh(4)
c               vnh(4) = vnh(4) + gnh(4)*dt4
c               gnh(3) = (qnh(2)*vnh(2)*vnh(2)-ekt) / qnh(3)
c               expterm = exp(-vnh(4)*dt8)
c               vnh(3) = expterm * (vnh(3)*expterm+gnh(3)*dt4)
c               gnh(2) = (qnh(1)*vnh(1)*vnh(1)-ekt) / qnh(2)
c               expterm = exp(-vnh(3)*dt8)
c               vnh(2) = expterm * (vnh(2)*expterm+gnh(2)*dt4)
c               gnh(1) = (2.0d0*eksum-dble(nfree)*ekt) / qnh(1)
c               expterm = exp(-vnh(2)*dt8)
c               vnh(1) = expterm * (vnh(1)*expterm+gnh(1)*dt4)
c               scale = scale * exp(-vnh(1)*dt2)
c               eksum = eksum * scale * scale
c               gnh(1) = (2.0d0*eksum-dble(nfree)*ekt) / qnh(1)
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
c      temp = 2.0d0 * eksum / (dble(nfree) * gasconst)
      return
      end
