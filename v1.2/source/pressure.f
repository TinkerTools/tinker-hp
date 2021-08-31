c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine pressure  --  constant pressure via barostat  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "pressure" uses the internal virial to find the pressure
c     in a periodic box and maintains a constant desired pressure
c     via a barostat method
c
c
      subroutine pressure (dt,ekin,pres,stress,istep)
      use bath
      use bound
      use boxes
      use domdec
      use units
      use virial
      implicit none
      integer i,j,istep
      real*8 dt
      real*8 pres
      real*8 factor
      real*8 ekin(3,3)
      real*8 stress(3,3)
c
c
c     only necessary if periodic boundaries are in use
c
      if (.not. use_bounds)  return
c
c     calculate the stress tensor for anisotropic systems
c
      factor = prescon / volbox
      do i = 1, 3
         do j = 1, 3
            stress(j,i) = factor * (2.0d0*ekin(j,i)-vir(j,i))
         end do
      end do
c
c     set isotropic pressure to the average of tensor diagonal
c
      pres = (stress(1,1)+stress(2,2)+stress(3,3)) / 3.0d0
c
c     use either the Berendsen or Monte Carlo barostat method
c
      if (isobaric) then
         if (barostat .eq. 'BERENDSEN') then
c           call pscale (dt,pres,stress,istep)
            call pscale (dt,pres,istep)
c         else if (barostat .eq. 'MONTECARLO') then
c            call pmonte(epot,temp)
         end if
      end if
      return
      end
c
c
c     "pressure2" applies a box size and velocity correction at
c     the half time step as needed for the Monte Carlo barostat
c
      subroutine pressure2 (epot,temp)
      use bath
      use bound
      use boxes
      use domdec
      use units
      use virial
      implicit none
      real*8 epot
      real*8 temp
c
c     only necessary if periodic boundaries are in use
c
      if (.not. use_bounds)  return
c
c
c     use either the Berendsen or Monte Carlo barostat method
c
      if (isobaric) then
c         if (barostat .eq. 'BERENDSEN') then
c            call pscale (dt,pres,stress,istep)
         if (barostat .eq. 'MONTECARLO') then
            call pmonte(epot,temp)
         end if
      end if
      return
      end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine pscale  --  Berendsen barostat via scaling  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "pscale" implements a Berendsen barostat by scaling the
c     coordinates and box dimensions via coupling to an external
c     constant pressure bath
c
c     literature references:
c
c     H. J. C. Berendsen, J. P. M. Postma, W. F. van Gunsteren,
c     A. DiNola and J. R. Hauk, "Molecular Dynamics with Coupling
c     to an External Bath", Journal of Chemical Physics, 81,
c     3684-3690 (1984)
c
c     S. E. Feller, Y. Zhang, R. W. Pastor, B. R. Brooks, "Constant
c     Pressure Molecular Dynamics Simulation: The Langevin Piston
c     Method", Journal of Chemical Physics, 103, 4613-4621 (1995)
c
c     code for anisotropic pressure coupling was provided by Guido
c     Raos, Dipartimento di Chimica, Politecnico di Milano, Italy
c
c
c     subroutine pscale (dt,pres,stress,istep)
      subroutine pscale (dt,pres,istep)
      use atoms
      use bath
      use boxes
      use domdec
      use math
      use usage
      implicit none
      integer i,iglob,istep

      real*8 dt,pres

      real*8 scale,third

c
c
c     find the isotropic scale factor for constant pressure
c
c      if (.not. anisotrop) then
         scale = 1.0d0
         third = 1.0d0 / 3.0d0
         scale = (1.0d0 + (dt*compress/taupres)*(pres-atmsph))**third
c
c     modify the current periodic box dimension values
c
         xbox = xbox * scale
         ybox = ybox * scale
         zbox = zbox * scale
c
c     propagate the new box dimensions to other lattice values
c
         call lattice
c
c     couple to pressure bath via atom scaling in Cartesian space
c
            do i = 1, nbloc
               iglob = glob(i)
               if (use(iglob)) then
                  x(iglob) = x(iglob) * scale
                  y(iglob) = y(iglob) * scale
                  z(iglob) = z(iglob) * scale
               end if
            end do
c
c   also rescale xbegproc, xendproc...
c
            call ddpme3dnpt(scale,istep)
cc
      return
      end
c
c     propagate Volume with Langevin equation with a BAOAB propagator
c
      subroutine plangevin ()
      use atoms
      use bath
      use boxes
      use domdec
      use math
      use moldyn
      use units
      use usage
      use mpi
      implicit none
      integer ierr


      real*8 third
      real*8 maxwell
      third = 1.0d0 / 3.0d0

      extvol = volbox
      extvolold = volbox
      if (rank.eq.0) then
        vextvol = maxwell(masspiston,kelvin)
      end if
      call MPI_BCAST(vextvol,1,MPI_REAL8,0,COMM_TINKER,ierr)
      aextvol = 0d0
      return
      end
c
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine ptest  --  find pressure via finite-difference  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "ptest" compares the virial-based value of dE/dV to an estimate
c     from finite-difference volume changes; also finds the isotropic
c     pressure via finite-differences
c
c     original version written by John D. Chodera, University of
c     California, Berkeley, December 2010
c
c
      subroutine ptest
      use atoms
      use bath
      use bound
      use boxes
      use domdec
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
c
c
c     set relative volume change for finite-differences
c
      if (.not. use_bounds)  return
      delta = 0.000001d0
      step = volbox * delta
c
c     perform dynamic allocation of some local arrays
c
      allocate (xoldloc(n))
      allocate (yoldloc(n))
      allocate (zoldloc(n))
c
c     store original box dimensions and coordinate values
c
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
c
c     get scale factor to reflect a negative volume change
c
      volbox = vold - step
      third = 1.0d0 / 3.0d0
      scale = (volbox/vold)**third
c
c     set new box dimensions and coordinate values
c
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
c
c     compute potential energy for negative volume change
c
      eneg = energy ()
      call allreduceen(eneg)
c
c     get scale factor to reflect a positive volume change
c
      volbox = vold + step
      third = 1.0d0 / 3.0d0
      scale = (volbox/vold)**third
c
c     set new box dimensions and coordinate values
c
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
c
c     compute potential energy for positive volume change
c
      epos = energy ()
      call allreduceen(epos)
c
c     restore original box dimensions and coordinate values
c
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
c
c     perform deallocation of some local arrays
c
      deallocate (xoldloc)
      deallocate (yoldloc)
      deallocate (zoldloc)
c
c     get virial and finite difference values of dE/dV
c
      dedv_vir = (vir(1,1)+vir(2,2)+vir(3,3)) / (3.0d0*volbox)
      dedv_fd = (epos-eneg) / (2.0d0*delta*volbox)
      if (rank.eq.0) then
        write (iout,10)  dedv_vir
   10   format (/,' dE/dV (Virial-based) :',11x,f15.6,' Kcal/mole/A**3')
        write (iout,20)  dedv_fd
   20   format (' dE/dV (Finite Diff) :',12x,f15.6,' Kcal/mole/A**3')
      end if
c
c     compute analytical and finite-difference isotropic pressure
c
      pres_vir = prescon * (dble(n)*gasconst*kelvin/volbox-dedv_vir)
      pres_fd = prescon * (dble(n)*gasconst*kelvin/volbox-dedv_fd)
      if (rank.eq.0) then
      if (kelvin .eq. 0.0d0) then
         write (iout,30)  pres_vir
         write (iout,40)  pres_fd
   30    format (/,' Pressure (Analytical, 0 K) :',5x,f15.3,
     &              ' Atmospheres')
   40    format (' Pressure (Numerical, 0 K) :',6x,f15.3,
     &              ' Atmospheres')
      else
         write (iout,50)  nint(kelvin),pres_vir
         write (iout,60)  nint(kelvin),pres_fd
   50    format (/,' Pressure (Analytical,',i4,' K) :',3x,f15.3,
     &              ' Atmospheres')
   60    format (' Pressure (Numerical,',i4,' K) :',4x,f15.3,
     &              ' Atmospheres')
      end if
      end if
      return
      end
c
      subroutine pmonte (epot,temp)
      use atmlst
      use atmtyp
      use atoms
      use bath
      use boxes
      use domdec
      use energi
      use group
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
      real*8 energy,random
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
      external random
c
c
c     decide whether to attempt a box size change at this step
c
      dotrial = .false.
      if (rank.eq.0) then
        valrand = random()
      end if
      call MPI_BCAST(valrand,1,MPI_REAL8,0,COMM_TINKER,ierr)
      call MPI_BCAST(epot,1,MPI_REAL8,0,COMM_TINKER,ierr)
c
      if (valrand .lt. 1.0d0/dble(voltrial))  dotrial = .true.
c
c     set constants and decide on type of trial box size change
c
      if (dotrial) then
         third = 1.0d0 / 3.0d0
         kt = gasconst * temp
         if (isothermal)  kt = gasconst * kelvin
         isotropic = .true.
c
c     perform dynamic allocation of some local arrays
c
         allocate (xold1(n))
         allocate (yold1(n))
         allocate (zold1(n))
         allocate (vold(3,n))
c
c     save the system state prior to trial box size change
c
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
c
c     for the isotropic case, change the lattice lengths uniformly
c
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
c            call ddpme3dnpt(scale,0)
           end if
         end if
c
c     get the potential energy and PV work changes for trial move
c
         epot = energy ()
         call allreduceen(epot)
         dpot = epot - eold
         dpv = atmsph * (volbox-volold) / prescon
c
c     estimate the kinetic energy change as an ideal gas term
c
         if (volscale .eq. 'MOLECULAR') then
           dkin = dble(nmol) * kt * log(volold/volbox)
         else
           dkin = dble(nuse) * kt * log(volold/volbox)
         end if
c
c     compute the kinetic energy change from the actual velocities;
c     scale the kinetic energy change to match virial pressure
c
c        dkin = 0.0d0
c        do i = 1, n
c           if (use(i)) then
c              term = 0.5d0 * mass(i) / convert
c              do j = 1, 3
c                 dkin = dkin + term*(v(j,i)**2-vold(j,i)**2)
c              end do
c           end if
c        end do
c        dkin = 0.907d0 * dkin
c
c     acceptance ratio from Epot change, Ekin change and PV work
c
         term = -(dpot+dpv+dkin) / kt
         expterm = exp(term)
c
c     reject the step, and restore values prior to trial change
c
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
c
c        rescale domain decomposition related stuff
c
         if (volscale.eq.'MOLECULAR') then
           call ddpme3dnpt(scale+1,0)
         else
           call ddpme3dnpt(scale,0)
         end if

c
c     perform deallocation of some local arrays
c
         deallocate (xold1)
         deallocate (yold1)
         deallocate (zold1)
         deallocate (vold)
      end if
      return
      end
c
c     subroutine rescale: rescale positions and speeds after a change of volume
c
      subroutine rescale(istep)
      use atoms
      use bath
      use boxes
      use domdec
      use moldyn
      use usage
      implicit none
      real*8 third,scale
      integer iglob,i,istep
      third = 1.0d0 / 3.0d0

      scale  =  (extvol/extvolold)**third
c
c     modify the current periodic box dimension values
c
      xbox = xbox * scale
      ybox = ybox * scale
      zbox = zbox * scale
c
c     propagate the new box dimensions to other lattice values
c
c      if (rank.eq.0) then
c      write(*,*) 'volume rescale = ',xbox*ybox*zbox
c      end if
      call lattice

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
c      if (rank.eq.0) then
c      write(*,*) 'scale = ',scale
c      end if
c
c   also rescale xbegproc, xendproc...
c
      call ddpme3dnpt(scale,istep)
      return
      end
