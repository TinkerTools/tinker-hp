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
      subroutine pressure (dt,epot,ekin,temp,pres,stress,istep)
      implicit none
      include 'sizes.i'
      include 'bath.i'
      include 'boxes.i'
      include 'bound.i'
      include 'units.i'
      include 'virial.i'
      include 'openmp.i'
      integer i,j,istep
      real*8 dt,epot
      real*8 temp,pres
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
            call pscale (dt,pres,stress,istep)
         else if (barostat .eq. 'MONTECARLO') then
            if (rank.eq.0) write(*,*) 'Monte-Carlo barostat not 
     $       implemented'
            call fatal
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
      subroutine pscale (dt,pres,stress,istep)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'bath.i'
      include 'boxes.i'
      include 'math.i'
      include 'openmp.i'
      include 'usage.i'
      integer i,j,k,iglob,istep
      integer start,stop
      real*8 dt,pres
      real*8 weigh,cosine
      real*8 scale,third
      real*8 xcm,xmove
      real*8 ycm,ymove
      real*8 zcm,zmove
      real*8 stress(3,3)
      real*8 temp(3,3)
      real*8 hbox(3,3)
      real*8 ascale(3,3)
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
c         if (integrate .ne. 'RIGIDBODY') then
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
cc     couple to pressure bath via center of mass of rigid bodies
cc
c         else
c            scale = scale - 1.0d0
c            do i = 1, ngrp
c               start = igrp(1,i)
c               stop = igrp(2,i)
c               xcm = 0.0d0
c               ycm = 0.0d0
c               zcm = 0.0d0
c               do j = start, stop
c                  k = kgrp(j)
c                  weigh = mass(k)
c                  xcm = xcm + x(k)*weigh
c                  ycm = ycm + y(k)*weigh
c                  zcm = zcm + z(k)*weigh
c               end do
c               xmove = scale * xcm/grpmass(i)
c               ymove = scale * ycm/grpmass(i)
c               zmove = scale * zcm/grpmass(i)
c               do j = start, stop
c                  k = kgrp(j)
c                  x(k) = x(k) + xmove
c                  y(k) = y(k) + ymove
c                  z(k) = z(k) + zmove
c               end do
c            end do
c         end if
cc
cc     find the anisotropic scale factors for constant pressure
cc
c      else
c         scale = dt*compress / (3.0d0*taupres)
c         do i = 1, 3
c            do j = 1, 3
c               if (j. eq. i) then
c                  ascale(j,i) = 1.0d0 + scale*(stress(i,i)-atmsph)
c               else
c                  ascale(j,i) = scale*stress(j,i)
c               end if
c            end do
c         end do
cc
cc     modify the current periodic box dimension values
cc
c         temp(1,1) = xbox
c         temp(2,1) = 0.0d0
c         temp(3,1) = 0.0d0
c         temp(1,2) = ybox * gamma_cos
c         temp(2,2) = ybox * gamma_sin
c         temp(3,2) = 0.0d0
c         temp(1,3) = zbox * beta_cos
c         temp(2,3) = zbox * beta_term
c         temp(3,3) = zbox * gamma_term
c         do i = 1, 3
c            do j = 1, 3
c               hbox(j,i) = 0.0d0
c               do k = 1, 3
c                  hbox(j,i) = hbox(j,i) + ascale(j,k)*temp(k,i)
c               end do
c            end do
c         end do
c         xbox = sqrt(hbox(1,1)**2 + hbox(2,1)**2 + hbox(3,1)**2)
c         ybox = sqrt(hbox(1,2)**2 + hbox(2,2)**2 + hbox(3,2)**2)
c         zbox = sqrt(hbox(1,3)**2 + hbox(2,3)**2 + hbox(3,3)**2)
c         if (monoclinic) then
c            cosine = (hbox(1,1)*hbox(1,3) + hbox(2,1)*hbox(2,3)
c     &                  + hbox(3,1)*hbox(3,3)) / (xbox*zbox)
c            beta = radian * acos(cosine)
c         else if (triclinic) then
c            cosine = (hbox(1,2)*hbox(1,3) + hbox(2,2)*hbox(2,3)
c     &                  + hbox(3,2)*hbox(3,3)) / (ybox*zbox)
c            alpha = radian * acos(cosine)
c            cosine = (hbox(1,1)*hbox(1,3) + hbox(2,1)*hbox(2,3)
c     &                  + hbox(3,1)*hbox(3,3)) / (xbox*zbox)
c            beta = radian * acos(cosine)
c            cosine = (hbox(1,1)*hbox(1,2) + hbox(2,1)*hbox(2,2)
c     &                  + hbox(3,1)*hbox(3,2)) / (xbox*ybox)
c            gamma = radian * acos(cosine)
c         end if
cc
cc     propagate the new box dimensions to other lattice values
cc
c         call lattice
cc
cc     couple to pressure bath via atom scaling in Cartesian space
cc
c         if (integrate .ne. 'RIGIDBODY') then
c            do i = 1, n
c               if (use(i)) then
c                  x(i) = ascale(1,1)*x(i) + ascale(1,2)*y(i)
c     &                      + ascale(1,3)*z(i)
c                  y(i) = ascale(2,1)*x(i) + ascale(2,2)*y(i)
c     &                      + ascale(2,3)*z(i)
c                  z(i) = ascale(3,1)*x(i) + ascale(3,2)*y(i)
c     &                      + ascale(3,3)*z(i)
c               end if
c            end do
cc
cc     couple to pressure bath via center of mass of rigid bodies
cc
c         else
c            ascale(1,1) = ascale(1,1) - 1.0d0
c            ascale(2,2) = ascale(2,2) - 1.0d0
c            ascale(3,3) = ascale(3,3) - 1.0d0
c            do i = 1, ngrp
c               start = igrp(1,i)
c               stop = igrp(2,i)
c               xcm = 0.0d0
c               ycm = 0.0d0
c               zcm = 0.0d0
c               do j = start, stop
c                  k = kgrp(j)
c                  weigh = mass(k)
c                  xcm = xcm + x(k)*weigh
c                  ycm = xcm + y(k)*weigh
c                  zcm = xcm + z(k)*weigh
c               end do
c               xcm = xcm / grpmass(i)
c               ycm = ycm / grpmass(i)
c               zcm = zcm / grpmass(i)
c               xmove = ascale(1,1)*xcm + ascale(1,2)*ycm
c     &                    + ascale(1,3)*zcm
c               ymove = ascale(2,1)*xcm + ascale(2,2)*ycm
c     &                    + ascale(2,3)*zcm
c               zmove = ascale(3,1)*xcm + ascale(3,2)*ycm
c     &                    + ascale(3,3)*zcm
c               do j = start, stop
c                  k = kgrp(j)
c                  x(k) = x(k) + xmove
c                  y(k) = y(k) + ymove
c                  z(k) = z(k) + zmove
c               end do
c            end do
c         end if
c      end if
      return
      end
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
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'bath.i'
      include 'bound.i'
      include 'boxes.i'
      include 'iounit.i'
      include 'units.i'
      include 'virial.i'
      integer i
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
      do i = 1, n
         xoldloc(i) = x(i)
         yoldloc(i) = y(i)
         zoldloc(i) = z(i)
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
      do i = 1, n
         x(i) = xoldloc(i) * scale
         y(i) = yoldloc(i) * scale
         z(i) = zoldloc(i) * scale
      end do
c
c     compute potential energy for negative volume change
c
c      eneg = energy ()
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
      do i = 1, n
         x(i) = xoldloc(i) * scale
         y(i) = yoldloc(i) * scale
         z(i) = zoldloc(i) * scale
      end do
c
c     compute potential energy for positive volume change
c
c      epos = energy ()
c
c     restore original box dimensions and coordinate values
c
      xbox = xboxold
      ybox = yboxold
      zbox = zboxold
      call lattice
      do i = 1, n
         x(i) = xoldloc(i)
         y(i) = yoldloc(i)
         z(i) = zoldloc(i)
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
      write (iout,10)  dedv_vir
   10 format (/,' dE/dV (Virial-based) :',11x,f15.6,' Kcal/mole/A**3')
      write (iout,20)  dedv_fd
   20 format (' dE/dV (Finite Diff) :',12x,f15.6,' Kcal/mole/A**3')
c
c     compute analytical and finite-difference isotropic pressure
c
      pres_vir = prescon * (dble(n)*gasconst*kelvin/volbox-dedv_vir)
      pres_fd = prescon * (dble(n)*gasconst*kelvin/volbox-dedv_fd)
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
      return
      end
