!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ###############################################################
!     ##                                                           ##
!     ##  subroutine pressure  --  constant pressure via barostat  ##
!     ##                                                           ##
!     ###############################################################
!
!
!     "pressure" uses the internal virial to find the pressure
!     in a periodic box and maintains a constant desired pressure
!     via a barostat method
!
!
!> @brief 
!> uses the internal virial to find the pressure
!> in a periodic box and maintains a constant desired pressure
!> via a barostat method
!> @param[in] ekin: kinetic energy
!> @param[in] vir: virial
!> @param[in] pres: instataneous pressure
!> @param[in] stress: stress tensor
subroutine stress_press(ekin,vir,pres,stress)
   use boxes
   use inform
   use iounit
   use units
   implicit none
   real*8, intent(in) :: ekin(3,3),vir(3,3)
   real*8, intent(out) :: pres,stress(3,3)
   real*8 factor
   integer i,j
!
   if (deb_Path) write(iout,*), 'stress_press '
!

!
!     calculate the stress tensor for anisotropic systems
!
   factor = prescon / volbox
   do i = 1, 3
      do j = 1, 3
         stress(j,i) = factor * (2.0d0*ekin(j,i)-vir(j,i))
      end do
   end do
!
!     set isotropic pressure to the average of tensor diagonal
!
   pres = (stress(1,1)+stress(2,2)+stress(3,3)) / 3.0d0
!
!     compute pressure via numerical virial if necessary
!
   ! call dedvcalc()
   ! if(kin_instant) then
   !    ekin_trace=(ekin(1,1)+ekin(2,2)+ekin(3,3))/corr_fact_qtb
   !    pres=prescon*(2.0d0*ekin_trace/(3.d0*volbox)-dedv)
   ! else
   !    pres=prescon*(nfree*kelvin*gasconst/(3.d0*volbox)-dedv)
   ! endif

end subroutine stress_press

!> @brief 
!> uses the internal virial to find the pressure
!> in a periodic box and maintains a constant desired pressure
!> via a barostat method
!> @param[in] dt: timestep
!> @param[in] ekin: kinetic energy
!> @param[in] vir: virial
!> @param[in] pres: instataneous pressure
!> @param[in] istep: index of timestep
subroutine pressure (dt,ekin,pres,stress,istep)
   use bath
   use bound
   use boxes
   use domdec
   use inform
   use iounit
   use mdstuf
   use units
   use virial
   implicit none
   integer istep
   real*8 dt
   real*8 pres
   real*8 ekin(3,3)
   real*8 stress(3,3)
!
   if (deb_Path) write(iout,*), 'pressure '
!
!
!
!     only necessary if periodic boundaries are in use
!
   if (.not. use_bounds)  return

   call stress_press(ekin,vir,pres,stress)

!     use either the Berendsen or Monte Carlo barostat method
!
   if (isobaric) then
      if (barostat .eq. 'BERENDSEN') then
         call pscale (dt,pres,stress,istep)
!         else if (barostat .eq. 'MONTECARLO') then
!            call pmonte(epot,temp)
      end if
   end if
   return
end
!
!
!     "pressure2" applies a box size and velocity correction at
!     the half time step as needed for the Monte Carlo barostat
!
!> @brief 
!> applies a box size and velocity correction at
!> the half time step as needed for the Monte Carlo barostat
!> @param[in] epot: potential energy
!> @param[in] temp: temperature in K
subroutine pressure2 (epot,temp)
   use bath
   use bound
   use boxes
   use domdec
   use inform
   use iounit
   use units
   use virial
   implicit none
   real*8 epot
   real*8 temp
!
   if (deb_Path) write(iout,*), 'pressure2 '
!
!
!     only necessary if periodic boundaries are in use
!
   if (.not. use_bounds)  return
!
!
!     use either the Berendsen or Monte Carlo barostat method
!
   if (isobaric) then
!         if (barostat .eq. 'BERENDSEN') then
!            call pscale (dt,pres,stress,istep)
      if (barostat .eq. 'MONTECARLO') then
         call pmonte(epot,temp)
      end if
   end if
   return
end
!
!
!     #############################################################
!     ##                                                         ##
!     ##  subroutine pscale  --  Berendsen barostat via scaling  ##
!     ##                                                         ##
!     #############################################################
!
!
!     "pscale" implements a Berendsen barostat by scaling the
!     coordinates and box dimensions via coupling to an external
!     constant pressure bath
!
!     literature references:
!
!     H. J. C. Berendsen, J. P. M. Postma, W. F. van Gunsteren,
!     A. DiNola and J. R. Hauk, "Molecular Dynamics with Coupling
!     to an External Bath", Journal of Chemical Physics, 81,
!     3684-3690 (1984)
!
!     S. E. Feller, Y. Zhang, R. W. Pastor, B. R. Brooks, "Constant
!     Pressure Molecular Dynamics Simulation: The Langevin Piston
!     Method", Journal of Chemical Physics, 103, 4613-4621 (1995)
!
!     code for anisotropic pressure coupling was provided by Guido
!     Raos, Dipartimento di Chimica, Politecnico di Milano, Italy
!
!
!> @brief 
!> implements a Berendsen barostat by scaling the
!> coordinates and box dimensions via coupling to an external
!> constant pressure bath
!> @param[in] dt: value of timestep
!> @param[in] pres: instantaneous pressure
!> @param[in] stress: stress tensor
!> @param[in] istep: index of timestep 
subroutine pscale (dt,pres,stress,istep)
   use atoms
   use bath
   use boxes
   use domdec
   use inform
   use iounit
   use math
   use usage
   implicit none
   integer i,iglob,istep

   real*8 :: dt,pres,stress(3,3)

   real*8 third
   real*8 :: scale(3)
!
   if (deb_Path) write(iout,*), 'pscale '
   third = 1.0d0 / 3.0d0

   if (anisotrop) then
!     find the anisotropic scale factor for constant pressure
      scale(:) = 1.0d0
      do i=1,3
         if (.not. freeze_axis(i)) then
            scale(i) = 1.0d0 + (third*dt*compress/taupres)*(stress(i,i)-atmsph)
         end if
      enddo

   else
!     find the isotropic scale factor for constant pressure
      scale(:) = (1.0d0 + (dt*compress/taupres)*(pres-atmsph))**third
   end if


   call rescale_box(istep,scale)
!
!     couple to pressure bath via atom scaling in Cartesian space
!
   do i = 1, nbloc
      iglob = glob(i)
      if (use(iglob)) then
         x(iglob) = x(iglob) * scale(1)
         y(iglob) = y(iglob) * scale(2)
         z(iglob) = z(iglob) * scale(3)
      end if
   end do

   return
end
!
!     propagate Volume with Langevin equation with a BAOAB propagator
!
!> @brief 
!> initialization of Langevin piston
!> @param no params
subroutine initialize_langevin_piston()
   use bath
   use boxes
   use mpi
   use domdec
   use inform
   use iounit
   use units, only: boltzmann
   implicit none
   integer ierr
   interface
      function maxwell (mass,temper)
         real*8 maxwell
         real*8 mass
         real*8 temper
      end function
   end interface
!
   if (deb_Path) write(iout,*), 'initialize_langevin_piston '
!

   extvol = volbox
   extvolold = volbox
   if (rank.eq.0) then
      if (.not. anisotrop) then
         vextvol  = maxwell(masspiston,kelvin)
      else
         vextbox(1) = maxwell(masspiston,kelvin)
         vextbox(2) = maxwell(masspiston,kelvin)
         vextbox(3) = maxwell(masspiston,kelvin)
      endif
   end if
   if (.not. anisotrop) then
      call MPI_BCAST(vextvol,1,MPI_REAL8,0,COMM_TINKER,ierr)
   else
      call MPI_BCAST(vextbox,3,MPI_REAL8,0,COMM_TINKER,ierr)
   endif
   aextvol = 0d0
   aextbox = 0d0
end subroutine initialize_langevin_piston
!
!
!     #################################################################
!     ##                                                             ##
!     ##  subroutine dedvcalc  --  find virial tensor via finite
!                                                        difference  ##
!     ##                                                             ##
!     #################################################################
!
!
!> @brief 
!> find virial tensor via finite differences
!> @param no params
subroutine dedvcalc()
   use atoms
   use bath
   use bound
   use boxes
   use domdec
   use inform
   use iounit
   use units
   use virial
   implicit none
   integer i,iglob
   real*8 energy,third
   real*8 delta,step,scale
   real*8 vold,xboxold
   real*8 yboxold,zboxold
   real*8 epos,eneg
   real*8, allocatable :: xoldloc(:)
   real*8, allocatable :: yoldloc(:)
   real*8, allocatable :: zoldloc(:)
!
   if (deb_Path) write(iout,*), 'dedvcalv '
!
!
!     set relative volume change for finite-differences
!
   if (.not. use_bounds)  return

   if(virnum) then
      delta = 0.000001d0
      step = volbox * delta
!
!     perform dynamic allocation of some local arrays
!
      allocate (xoldloc(n))
      allocate (yoldloc(n))
      allocate (zoldloc(n))
!
!     store original box dimensions and coordinate values
!
      xboxold = xbox
      yboxold = ybox
      zboxold = zbox
      vold = volbox
      do i = 1, nbloc
         iglob = glob(i)
         xoldloc(iglob) = x(iglob)
         yoldloc(iglob) = y(iglob)
         zoldloc(iglob) = z(iglob)
      end do
!
!     get scale factor to reflect a negative volume change
!
      volbox = vold - step
      third = 1.0d0 / 3.0d0
      scale = (volbox/vold)**third
!
!     set new box dimensions and coordinate values
!
      xbox = xboxold * scale
      ybox = yboxold * scale
      zbox = zboxold * scale
      call lattice
      do i = 1, nbloc
         iglob = glob(i)
         x(iglob) = xoldloc(iglob) * scale
         y(iglob) = yoldloc(iglob) * scale
         z(iglob) = zoldloc(iglob) * scale
      end do
!
!     compute potential energy for negative volume change
!
      eneg = energy ()
      call allreduceen(eneg)
!
!     get scale factor to reflect a positive volume change
!
      volbox = vold + step
      third = 1.0d0 / 3.0d0
      scale = (volbox/vold)**third
!
!     set new box dimensions and coordinate values
!
      xbox = xboxold * scale
      ybox = yboxold * scale
      zbox = zboxold * scale
      call lattice
      do i = 1, nbloc
         iglob = glob(i)
         x(iglob) = xoldloc(iglob) * scale
         y(iglob) = yoldloc(iglob) * scale
         z(iglob) = zoldloc(iglob) * scale
      end do
!
!     compute potential energy for positive volume change
!
      epos = energy ()
      call allreduceen(epos)
!
!     restore original box dimensions and coordinate values
!
      xbox = xboxold
      ybox = yboxold
      zbox = zboxold
      call lattice
      do i = 1, nbloc
         iglob = glob(i)
         x(iglob) = xoldloc(iglob)
         y(iglob) = yoldloc(iglob)
         z(iglob) = zoldloc(iglob)
      end do
!
!     perform deallocation of some local arrays
!
      deallocate (xoldloc)
      deallocate (yoldloc)
      deallocate (zoldloc)
!
!     get virial and finite difference values of dE/dV
!
      dedv = (epos-eneg) / (2.0d0*step)
   else
      dedv = (vir(1,1)+vir(2,2)+vir(3,3)) / (3.0d0*volbox)
   endif

end subroutine dedvcalc
!

!
!
!
!
!
!     #################################################################
!     ##                                                             ##
!     ##  subroutine ptest  --  find pressure via finite-difference  ##
!     ##                                                             ##
!     #################################################################
!
!
!     "ptest" compares the virial-based value of dE/dV to an estimate
!     from finite-difference volume changes; also finds the isotropic
!     pressure via finite-differences
!
!     original version written by John D. Chodera, University of
!     California, Berkeley, December 2010
!
!
!> @brief 
!> "ptest" compares the virial-based value of dE/dV to an estimate
!> from finite-difference volume changes; also finds the isotropic
!> pressure via finite-differences
!> @param no params
subroutine ptest
   use atoms
   use bath
   use bound
   use boxes
   use domdec
   use inform
   use iounit
   use units
   use virial
   implicit none
   integer i,iglob
   real*8 energy,third
   real*8 delta,step,scale
   real*8 vold,xboxold
   real*8 yboxold,zboxold
   real*8 epos,eneg
   real*8 dedv_vir,dedv_fd
   real*8 pres_vir,pres_fd
   real*8, allocatable :: xoldloc(:)
   real*8, allocatable :: yoldloc(:)
   real*8, allocatable :: zoldloc(:)
!
   if (deb_Path) write(iout,*), 'ptest '
!
!
!
!     set relative volume change for finite-differences
!
   if (.not. use_bounds)  return
   delta = 0.000001d0
   step = volbox * delta
!
!     perform dynamic allocation of some local arrays
!
   allocate (xoldloc(n))
   allocate (yoldloc(n))
   allocate (zoldloc(n))
!
!     store original box dimensions and coordinate values
!
   xboxold = xbox
   yboxold = ybox
   zboxold = zbox
   vold = volbox
   do i = 1, nbloc
      iglob = glob(i)
      xoldloc(iglob) = x(iglob)
      yoldloc(iglob) = y(iglob)
      zoldloc(iglob) = z(iglob)
   end do
!
!     get scale factor to reflect a negative volume change
!
   volbox = vold - step
   third = 1.0d0 / 3.0d0
   scale = (volbox/vold)**third
!
!     set new box dimensions and coordinate values
!
   xbox = xboxold * scale
   ybox = yboxold * scale
   zbox = zboxold * scale
   call lattice
   do i = 1, nbloc
      iglob = glob(i)
      x(iglob) = xoldloc(iglob) * scale
      y(iglob) = yoldloc(iglob) * scale
      z(iglob) = zoldloc(iglob) * scale
   end do
!
!     compute potential energy for negative volume change
!
   eneg = energy ()
   call allreduceen(eneg)
!
!     get scale factor to reflect a positive volume change
!
   volbox = vold + step
   third = 1.0d0 / 3.0d0
   scale = (volbox/vold)**third
!
!     set new box dimensions and coordinate values
!
   xbox = xboxold * scale
   ybox = yboxold * scale
   zbox = zboxold * scale
   call lattice
   do i = 1, nbloc
      iglob = glob(i)
      x(iglob) = xoldloc(iglob) * scale
      y(iglob) = yoldloc(iglob) * scale
      z(iglob) = zoldloc(iglob) * scale
   end do
!
!     compute potential energy for positive volume change
!
   epos = energy ()
   call allreduceen(epos)
!
!     restore original box dimensions and coordinate values
!
   xbox = xboxold
   ybox = yboxold
   zbox = zboxold
   call lattice
   do i = 1, nbloc
      iglob = glob(i)
      x(iglob) = xoldloc(iglob)
      y(iglob) = yoldloc(iglob)
      z(iglob) = zoldloc(iglob)
   end do
!
!     perform deallocation of some local arrays
!
   deallocate (xoldloc)
   deallocate (yoldloc)
   deallocate (zoldloc)
!
!     get virial and finite difference values of dE/dV
!
   dedv_vir = (vir(1,1)+vir(2,2)+vir(3,3)) / (3.0d0*volbox)
   dedv_fd = (epos-eneg) / (2.0d0*delta*volbox)
   if (rank.eq.0) then
      write (iout,10)  dedv_vir
10    format (/,' dE/dV (Virial-based) :',11x,f15.6,' Kcal/mole/A**3')
      write (iout,20)  dedv_fd
20    format (' dE/dV (Finite Diff) :',12x,f15.6,' Kcal/mole/A**3')
   end if
!
!     compute analytical and finite-difference isotropic pressure
!
   pres_vir = prescon * (dble(n)*gasconst*kelvin/volbox-dedv_vir)
   pres_fd = prescon * (dble(n)*gasconst*kelvin/volbox-dedv_fd)
   if (rank.eq.0) then
      if (kelvin .eq. 0.0d0) then
         write (iout,30)  pres_vir
         write (iout,40)  pres_fd
30       format (/,' Pressure (Analytical, 0 K) :',5x,f15.3,&
         &' Atmospheres')
40       format (' Pressure (Numerical, 0 K) :',6x,f15.3,&
         &' Atmospheres')
      else
         write (iout,50)  nint(kelvin),pres_vir
         write (iout,60)  nint(kelvin),pres_fd
50       format (/,' Pressure (Analytical,',i4,' K) :',3x,f15.3,&
         &' Atmospheres')
60       format (' Pressure (Numerical,',i4,' K) :',4x,f15.3,&
         &' Atmospheres')
      end if
   end if
   return
end
!
!> @brief 
!> applies monte-carlo barostat
!> @param[in] epot: potential energy
!> @param[in] temp: temperature
subroutine pmonte (epot,temp)
   use atmlst
   use atmtyp
   use atoms
   use bath
   use boxes
   use domdec
   use energi
   use group
   use inform
   use iounit
   use math
   use mdstuf
   use molcul
   use moldyn
   use units
   use usage
   use mpi
   implicit none
   integer i,j,k,iglob,ierr
   integer start,stop
   real*8 epot,temp,term
   real*8 energy
   real*8 kt,expterm
   real*8 third,weigh
   real*8 step,scale
   real*8 eold
   real*8 xcm,ycm,zcm
   real*8 vxcm,vycm,vzcm
   real*8 volold
   real*8 dpot,dpv,dkin
   real*8 xmove,ymove,zmove
   real*8 vxmove,vymove,vzmove
   real*8 xboxold,yboxold,zboxold
   real*8 alphaold,betaold,gammaold



   real*8, allocatable :: xold1(:)
   real*8, allocatable :: yold1(:)
   real*8, allocatable :: zold1(:)
   real*8, allocatable :: vold(:,:)
   real*8 valrand
   logical dotrial
   logical isotropic
   real*8 random
   external random
!
   if (deb_Path) write(iout,*), 'pmonte '
!
!
!
!     decide whether to attempt a box size change at this step
!
   dotrial = .false.
   if (rank.eq.0) then
      valrand = random()
   end if
   call MPI_BCAST(valrand,1,MPI_REAL8,0,COMM_TINKER,ierr)
   call MPI_BCAST(epot,1,MPI_REAL8,0,COMM_TINKER,ierr)
!
   if (valrand .lt. 1.0d0/dble(voltrial))  dotrial = .true.
!
!     set constants and decide on type of trial box size change
!
   if (dotrial) then
      third = 1.0d0 / 3.0d0
      kt = gasconst * temp
      if (isothermal)  kt = gasconst * kelvin
      isotropic = .true.
!
!     perform dynamic allocation of some local arrays
!
      allocate (xold1(n))
      allocate (yold1(n))
      allocate (zold1(n))
      allocate (vold(3,n))
!
!     save the system state prior to trial box size change
!
      xboxold = xbox
      yboxold = ybox
      zboxold = zbox
      alphaold = alpha
      betaold  = beta
      gammaold = gamma
      volold = volbox
      eold = epot
      do i = 1, nbloc
         iglob = glob(i)
         if (use(iglob)) then
            xold1(iglob) = x(iglob)
            yold1(iglob) = y(iglob)
            zold1(iglob) = z(iglob)
            vold(1,iglob) = v(1,iglob)
            vold(2,iglob) = v(2,iglob)
            vold(3,iglob) = v(3,iglob)
         end if
      end do
!
!     for the isotropic case, change the lattice lengths uniformly
!
      if (isotropic) then
         if (rank.eq.0) then
            valrand = random()
         end if
         call MPI_BCAST(valrand,1,MPI_REAL8,0,COMM_TINKER,ierr)
         step = volmove * (2.0d0*valrand-1.0d0)
         volbox = volbox + step
         scale = (volbox/volold)**third
         xbox = xbox * scale
         ybox = ybox * scale
         zbox = zbox * scale
         call lattice
         if (volscale .eq. 'MOLECULAR') then
            call molecule(.false.)
            scale = scale - 1.0d0
            do i = 1, nmoleloc
               iglob = molculeglob(i)
               xcm = 0.0d0
               ycm = 0.0d0
               zcm = 0.0d0
               vxcm = 0.0d0
               vycm = 0.0d0
               vzcm = 0.0d0
               start = imol(1,iglob)
               stop = imol(2,iglob)
               do j = start, stop
                  k = kmol(j)
                  weigh = mass(k)
                  xcm = xcm + x(k)*weigh
                  ycm = ycm + y(k)*weigh
                  zcm = zcm + z(k)*weigh
                  vxcm = vxcm + v(1,k)*weigh
                  vycm = vycm + v(2,k)*weigh
                  vzcm = vzcm + v(3,k)*weigh
               end do
               xmove = scale * xcm/molmass(iglob)
               ymove = scale * ycm/molmass(iglob)
               zmove = scale * zcm/molmass(iglob)
               vxmove = scale * vxcm/molmass(iglob)
               vymove = scale * vycm/molmass(iglob)
               vzmove = scale * vzcm/molmass(iglob)
               do j = start, stop
                  k = kmol(j)
                  if (use(k)) then
                     x(k) = x(k) + xmove
                     y(k) = y(k) + ymove
                     z(k) = z(k) + zmove
                     v(1,k) = v(1,k) - vxmove
                     v(2,k) = v(2,k) - vymove
                     v(3,k) = v(3,k) - vzmove
                  end if
               end do
            end do
         else
            do i = 1, nbloc
               iglob = glob(i)
               if (use(iglob)) then
                  x(iglob) = x(iglob) * scale
                  y(iglob) = y(iglob) * scale
                  z(iglob) = z(iglob) * scale
                  v(1,iglob) = v(1,iglob) / scale
                  v(2,iglob) = v(2,iglob) / scale
                  v(3,iglob) = v(3,iglob) / scale
               end if
            end do
!            call ddpme3dnpt(scale,0)
         end if
      end if
!
!     get the potential energy and PV work changes for trial move
!
      epot = energy ()
      call allreduceen(epot)
      dpot = epot - eold
      dpv = atmsph * (volbox-volold) / prescon
!
!     estimate the kinetic energy change as an ideal gas term
!
      if (volscale .eq. 'MOLECULAR') then
         dkin = dble(nmol) * kt * log(volold/volbox)
      else
         dkin = dble(nuse) * kt * log(volold/volbox)
      end if
!
!     compute the kinetic energy change from the actual velocities;
!     scale the kinetic energy change to match virial pressure
!
!        dkin = 0.0d0
!        do i = 1, n
!           if (use(i)) then
!              term = 0.5d0 * mass(i) / convert
!              do j = 1, 3
!                 dkin = dkin + term*(v(j,i)**2-vold(j,i)**2)
!              end do
!           end if
!        end do
!        dkin = 0.907d0 * dkin
!
!     acceptance ratio from Epot change, Ekin change and PV work
!
      term = -(dpot+dpv+dkin) / kt
      expterm = exp(term)
!
!     reject the step, and restore values prior to trial change
!
      if (rank.eq.0) then
         valrand = random()
      end if
      call MPI_BCAST(valrand,1,MPI_REAL8,0,COMM_TINKER,ierr)
      if (valrand .gt. expterm) then
         epot = eold
         xbox = xboxold
         ybox = yboxold
         zbox = zboxold
         call lattice
         do i = 1, nbloc
            iglob = glob(i)
            if (use(iglob)) then
               x(iglob) = xold1(iglob)
               y(iglob) = yold1(iglob)
               z(iglob) = zold1(iglob)
               v(1,iglob) = vold(1,iglob)
               v(2,iglob) = vold(2,iglob)
               v(3,iglob) = vold(3,iglob)
            end if
         end do
         return
      end if
!
!        rescale domain decomposition related stuff
!
      if (volscale.eq.'MOLECULAR') then
         call ddpme3dnpt(scale+1,0)
      else
         call ddpme3dnpt(scale,0)
      end if

!
!     perform deallocation of some local arrays
!
      deallocate (xold1)
      deallocate (yold1)
      deallocate (zold1)
      deallocate (vold)
   end if
   return
end
!
!     subroutine rescale_box: rescale simulation box (isotropic)
!
!> @brief 
!> rescale simulation box (isotropic)
!> @param[in] istep: index of timestep
!> @param[in] scale: scaling factor
subroutine rescale_box(istep,scale)
   use atoms
   use bath
   use boxes
   use domdec
   use inform
   use iounit
   use moldyn
   use usage
   implicit none
   integer, intent(in) :: istep
   real*8, intent(in) :: scale(3)
!
   if (deb_Path) write(iout,*), 'rescale_box '
!
!
!     modify the current periodic box dimension values
!
   xbox = xbox * scale(1)
   ybox = ybox * scale(2)
   zbox = zbox * scale(3)
!
!     propagate the new box dimensions to other lattice values
!
   call lattice
!
!   also rescale xbegproc, xendproc...
!
   call ddpme3dnptaniso(scale,istep)

end subroutine rescale_box
