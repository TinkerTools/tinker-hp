!
!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ##############################################################
!     ##                                                          ##
!     ##  subroutine temper  --  thermostat applied at half step  ##
!     ##                                                          ##
!     ##############################################################
!
!
!     "temper" applies a velocity correction at the half time step
!     as needed for the Nose-Hoover extended system thermostat
!
!     literature references:
!
!     D. Frenkel and B. Smit, "Understanding Molecular Simulation,
!     2nd Edition", Academic Press, San Diego, CA, 2002; see Appendix
!     E.2 for implementation details
!
!     G. J. Martyna, M. E. Tuckerman, D. J. Tobias and M. L. Klein,
!     "Explicit Reversible Integrators for Extended Systems Dynamics",
!     Molecular Physics, 87, 1117-1157 (1996)
!
!
!> @brief 
!> computes and displays the total potential energy
!> @param[in] dt: duration of a timestep
!> @param[in] eksum: kinetic energy
!> @param[in] ekin: kinetic energy tensor
!> @param[in] temp: instantaneous temperature
subroutine temper (dt,eksum,ekin,temp)
   use atmlst
   use atmtyp
   use bath
   use domdec
   use group
   use inform
   use iounit
   use mdstuf
   use molcul
   use moldyn
   use units
   use usage
   use mpi
   implicit none
   integer i,j,iglob,ierr
   integer k,m
   real*8 dt
   real*8 eksum
   real*8 temp
   real*8 scale,speed
   real*8 ekin(3,3)
   real*8 c,d,r,s,si
   real*8 random,normal
   real*8 kt,rate,trial
!
   if (deb_Path) write(iout,*), 'temper '
!
!
   call kinetic (eksum,ekin,temp)
   if (.not. isothermal)  return
!
!     couple to external temperature bath via Berendsen scaling
!
   if (thermostat .eq. 'BERENDSEN') then
      scale = 1.0d0
      if (temp .ne. 0.0d0)&
      &scale = sqrt(1.0d0 + (dt/tautemp)*(kelvin/temp-1.0d0))
      do i = 1, nloc
         iglob = glob(i)
         if (use(iglob)) then
            do j = 1, 3
               v(j,iglob) = scale * v(j,iglob)
            end do
         end if
      end do
!
!     couple to external temperature bath via Bussi scaling
!
   else if (thermostat .eq. 'BUSSI') then
      if (rank.eq.0) then
         if (temp .eq. 0.0d0)  temp = 0.1d0
         c = exp(-dt/tautemp)
         d = (1.0d0-c) * (kelvin/temp) / dble(nfree)
         if (rank.eq.0) then
            r = normal ()
         end if
         s = 0.0d0
         do i = 1, nfree-1
            si = normal ()
            s = s + si*si
         end do
         scale = c + (s+r*r)*d + 2.0d0*r*sqrt(c*d)
         scale = sqrt(scale)
         if (r+sqrt(c/d) .lt. 0.0d0)  scale = -scale
      end if
      call MPI_BCAST(scale,1,MPI_REAL8,0,COMM_TINKER,ierr)
      eta = eta * scale
      do i = 1, nloc
         iglob = glob(i)
         if (use(iglob)) then
            do j = 1, 3
               v(j,iglob) = scale * v(j,iglob)
            end do
         end if
      end do
!
!     select random velocities via Andersen stochastic collisions
!
   else if (thermostat .eq. 'ANDERSEN') then
      kt = boltzmann * kelvin
      rate = 1000.0d0 * dt * collide
      if (barostat.eq.'MONTECARLO' .and.&
      &volscale.eq.'MOLECULAR') then
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
!
!     recompute kinetic energy and instantaneous temperature
!
   call kinetic (eksum,ekin,temp)
   temp = 2.0d0 * eksum / (dble(nfree) * gasconst)
!
   return
end
!
!
!     ###############################################################
!     ##                                                           ##
!     ##  subroutine temper2  --  thermostat applied at full step  ##
!     ##                                                           ##
!     ###############################################################
!
!
!     "temper2" computes the instantaneous temperature and applies a
!     thermostat via Berendsen or Bussi-Parrinello velocity scaling,
!     Andersen stochastic collisions or Nose-Hoover extended system
!
!     literature references:
!
!     H. J. C. Berendsen, J. P. M. Postma, W. F. van Gunsteren,
!     A. DiNola and J. R. Hauk, "Molecular Dynamics with Coupling
!     to an External Bath", Journal of Chemical Physics, 81,
!     3684-3690 (1984)
!
!     G. Bussi and M. Parrinello, "Stochastic Thermostats: Comparison
!     of Local and Global Schemes", Computer Physics Communications,
!     179, 26-29 (2008)
!
!     H. C. Andersen, "Molecular Dynamics Simulations at Constant
!     Pressure and/or Temperature", Journal of Chemical Physics,
!     72, 2384-2393 (1980)
!
!
!> @brief 
!> computes the instantaneous temperature and applies a
!> thermostat via Berendsen or Bussi-Parrinello velocity scaling,
!> Andersen stochastic collisions or Nose-Hoover extended system
!> @param[in] temp: instaneous temperature
subroutine temper2 (temp)
   use atmtyp
   use bath
   use domdec
   use group
   use mdstuf
   use moldyn
   use units
   use usage
   implicit none
   real*8 eksum
   real*8 temp
   real*8 ekin(3,3)
!
!
!     get instantaneous temperature from the kinetic energy
!
   call kinetic (eksum,ekin,temp)
   if (.not. isothermal)  return
   return
!c
!c     couple to external temperature bath via Berendsen scaling
!c
!      if (thermostat .eq. 'BERENDSEN') then
!         call kinetic (eksum,ekin)
!         temp = 2.0d0 * eksum / (dble(nfree) * gasconst)
!         if (temp .eq. 0.0d0)  temp = 0.1d0
!         scale = sqrt(1.0d0 + (dt/tautemp)*(kelvin/temp-1.0d0))
!         do i = 1, nloc
!            iglob = glob(i)
!            if (use(iglob)) then
!               do j = 1, 3
!                  v(j,iglob) = scale * v(j,iglob)
!               end do
!            end if
!         end do
!c
!c     couple to external temperature bath via Bussi scaling
!c
!      else if (thermostat .eq. 'BUSSI') then
!         call kinetic (eksum,ekin)
!         temp = 2.0d0 * eksum / (dble(nfree) * gasconst)
!         if (temp .eq. 0.0d0)  temp = 0.1d0
!         c = exp(-dt/tautemp)
!         d = (1.0d0-c) * (kelvin/temp) / dble(nfree)
!         r = normal ()
!         s = 0.0d0
!         do i = 1, nfree-1
!            si = normal ()
!            s = s + si*si
!         end do
!         scale = c + (s+r*r)*d + 2.0d0*r*sqrt(c*d)
!         scale = sqrt(scale)
!         if (r+sqrt(c/d) .lt. 0.0d0)  scale = -scale
!         eta = eta * scale
!         do i = 1, nloc
!            iglob = glob(i)
!            if (use(iglob)) then
!               do j = 1, 3
!                  v(j,iglob) = scale * v(j,iglob)
!               end do
!            end if
!         end do
!c
!c     select random velocities via Andersen stochastic collisions
!c
!      else if (thermostat .eq. 'ANDERSEN') then
!         kt = boltzmann * kelvin
!         rate = 1000.0d0 * dt * collide
!         rate = rate / dble(nuse)**(2.0d0/3.0d0)
!         do i = 1, nloc
!            iglob = glob(i)
!            if (use(iglob)) then
!               trial = random ()
!               if (trial .lt. rate) then
!                  speed = sqrt(kt/mass(iglob))
!                  do j = 1, 3
!                     v(j,iglob) = speed * normal ()
!                  end do
!               end if
!            end if
!         end do
!c
!c     make full-step velocity correction for Nose-Hoover system
!c
!      else if (thermostat .eq. 'NOSE-HOOVER') then
!         ekt = gasconst * kelvin
!         nc = 5
!         ns = 3
!         dtc = dt / dble(nc)
!         w(1) = 1.0d0 / (2.0d0-2.0d0**(1.0d0/3.0d0))
!         w(2) = 1.0d0 - 2.0d0*w(1)
!         w(3) = w(1)
!         scale = 1.0d0
!         do i = 1, nc
!            do j = 1, ns
!               dts = w(j) * dtc
!               dt2 = 0.5d0 * dts
!               dt4 = 0.25d0 * dts
!               dt8 = 0.125d0 * dts
!               gnh(4) = (qnh(3)*vnh(3)*vnh(3)-ekt) / qnh(4)
!               vnh(4) = vnh(4) + gnh(4)*dt4
!               gnh(3) = (qnh(2)*vnh(2)*vnh(2)-ekt) / qnh(3)
!               expterm = exp(-vnh(4)*dt8)
!               vnh(3) = expterm * (vnh(3)*expterm+gnh(3)*dt4)
!               gnh(2) = (qnh(1)*vnh(1)*vnh(1)-ekt) / qnh(2)
!               expterm = exp(-vnh(3)*dt8)
!               vnh(2) = expterm * (vnh(2)*expterm+gnh(2)*dt4)
!               gnh(1) = (2.0d0*eksum-dble(nfree)*ekt) / qnh(1)
!               expterm = exp(-vnh(2)*dt8)
!               vnh(1) = expterm * (vnh(1)*expterm+gnh(1)*dt4)
!               scale = scale * exp(-vnh(1)*dt2)
!               eksum = eksum * scale * scale
!               gnh(1) = (2.0d0*eksum-dble(nfree)*ekt) / qnh(1)
!               expterm = exp(-vnh(2)*dt8)
!               vnh(1) = expterm * (vnh(1)*expterm+gnh(1)*dt4)
!               gnh(2) = (qnh(1)*vnh(1)*vnh(1)-ekt) / qnh(2)
!               expterm = exp(-vnh(3)*dt8)
!               vnh(2) = expterm * (vnh(2)*expterm+gnh(2)*dt4)
!               gnh(3) = (qnh(2)*vnh(2)*vnh(2)-ekt) / qnh(3)
!               expterm = exp(-vnh(4)*dt8)
!               vnh(3) = expterm * (vnh(3)*expterm+gnh(3)*dt4)
!               gnh(4) = (qnh(3)*vnh(3)*vnh(3)-ekt) / qnh(4)
!               vnh(4) = vnh(4) + gnh(4)*dt4
!            end do
!         end do
!         do i = 1, nloc
!            iglob = glob(i)
!            if (use(iglob)) then
!               do j = 1, 3
!                  v(j,iglob) = scale * v(j,iglob)
!               end do
!            end if
!         end do
!      end if
!c
!c     recompute kinetic energy and instantaneous temperature
!c
!      call kinetic (eksum,ekin)
!      temp = 2.0d0 * eksum / (dble(nfree) * gasconst)
   return
end
