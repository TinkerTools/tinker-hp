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
#include "tinker_precision.h"
subroutine pressure (dt,ekin,pres,stress,istep)
   use bath
   use bound
   use boxes
   use domdec
   use tinheader ,only:ti_p,re_p
   use units
   use virial
   implicit none
   integer i,j,istep
   real(r_p) dt
   real(r_p) ekin(3,3)
   real(r_p) pres
   real(r_p) factor
   real(r_p) stress(3,3)
   real(r_p) stres1,stres2,stres3
!
!     only necessary if periodic boundaries are in use
!     and isobaric simulation
!
   if (.not.(use_bounds.and.use_virial))  return
!
!     calculate the stress tensor for anisotropic systems
!
   factor = prescon / volbox
   if (anisotrop) then
!$acc parallel loop collapse(2) default(present) async
      do i = 1, 3
         do j = 1, 3
            stress(j,i) = factor * (2*ekin(j,i)-vir(j,i))
         end do
      end do
   endif
!
!     set isotropic pressure to the average of tensor diagonal
!
!$acc host_data use_device(ekin,stress,pres)
!$acc serial async deviceptr(ekin,stress,pres)
   !print*,'ek ',ekin(1,1),ekin(2,2),ekin(3,3)
   !print*,'vir', vir(1,1), vir(2,2), vir(3,3)
   stres1 = factor * (2*ekin(1,1)-vir(1,1))
   stres2 = factor * (2*ekin(2,2)-vir(2,2))
   stres3 = factor * (2*ekin(3,3)-vir(3,3))
   pres   = (stres1+stres2+stres3) / 3.0_re_p
!$acc end serial
!$acc end host_data
!
!     use the Berendsen barostat method
!
   if (isobaric) then
      if (barostat.eq.'BERENDSEN') call pscale(dt,pres,stress,istep)
   end if
end
!
!     "pressure2" applies a box size and velocity correction at
!     the half time step as needed for the Monte Carlo barostat
!
subroutine pressure2 (epot,temp)
   use bath
   use bound
   use boxes
   use domdec
   use units
   use virial
   implicit none
   real(r_p) epot,temp
!
!     only necessary if periodic boundaries are in use
!     and isobaric simulation
!
   if (.not.(use_bounds.and.isobaric))  return

   ! --- Monte Carlo barostat method --- !
   if (barostat.eq.'MONTECARLO') call pmonte(epot,temp)

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
subroutine pscale (dt,pres,stress,istep)
   use atomsMirror
   use bath
   use boxes
   use domdec
   use math
   use inform    ,only:deb_Path
   use tinheader ,only:ti_p,re_p
   use usage
   implicit none
   integer i,j,k,iglob,istep
   integer start,stop
   real(r_p) dt,pres
   real(t_p) weigh,cosine
   real(r_p) scale(3),third,scalex,scaley,scalez
   real(r_p) xcm,xmove
   real(r_p) ycm,ymove
   real(r_p) zcm,zmove
   real(r_p) stress(3,3)
   real(t_p) temp(3,3)
   real(t_p) hbox(3,3)
   real(t_p) ascale(3,3)
   parameter(third = 1.0_re_p / 3.0_re_p)
!
!
!     find the isotropic scale factor for constant pressure
!
!      if (.not. anisotrop) then
!$acc wait
!$acc update host(pres)

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
!
!     modify the current periodic box dimension values
   call rescale_box(istep,scale)

   if (deb_Path) then
13    format(A,10F16.6)
      print 13,'pscale',scale(1),scale(2),scale(3),dt,compress,pres,atmsph,xbox,ybox,zbox
   end if
!
!     couple to pressure bath via atom scaling in Cartesian space
!
   scalex = scale(1)
   scaley = scale(2)
   scalez = scale(3)
!$acc parallel loop async &
!$acc         present(glob,x,y,z)
   do i = 1, nbloc
      iglob = glob(i)
      if (use(iglob)) then
         x(iglob) = x(iglob) * scalex
         y(iglob) = y(iglob) * scaley
         z(iglob) = z(iglob) * scalez
      end if
   end do
   call reCast_position
!
end
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
subroutine ptest
   use atoms
   use bath
   use bound
   use boxes
   use domdec
   use iounit
   use tinheader ,only:ti_p,re_p
   use units
   use virial
   implicit none
   integer i,iglob
   real(r_p) energy
   real(r_p) third
   real(r_p) delta,step,scale
   real(r_p) vold,xboxold
   real(r_p) yboxold,zboxold
   real(r_p) epos,eneg
   real(r_p) dedv_vir,dedv_fd
   real(r_p) pres_vir,pres_fd
   real(r_p), allocatable :: xoldloc(:)
   real(r_p), allocatable :: yoldloc(:)
   real(r_p), allocatable :: zoldloc(:)
!
!
!     set relative volume change for finite-differences
!
   if (.not. use_bounds)  return
   delta = 0.000001_re_p
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
   third = 1.0_re_p / 3.0_re_p
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
   scale = (volbox/vold)**third
!$acc update device(volbox) async
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
!$acc update device(xbox,ybox,zbox) async
   call lattice
   do i = 1, nbloc
      iglob = glob(i)
      x(iglob) = xoldloc(iglob)
      y(iglob) = yoldloc(iglob)
      z(iglob) = zoldloc(iglob)
   end do
!$acc update device(x(:),y(:),z(:))
!
!     perform deallocation of some local arrays
!
   deallocate (xoldloc)
   deallocate (yoldloc)
   deallocate (zoldloc)
!
!     get virial and finite difference values of dE/dV
!
   dedv_vir = (vir(1,1)+vir(2,2)+vir(3,3)) / (3.0_re_p*volbox)
   dedv_fd = real(epos-eneg,r_p) / (2.0_re_p*delta*volbox)
   if (rank.eq.0) then
      write (iout,10)  dedv_vir
10    format (/,' dE/dV (Virial-based) :',11x,f15.6,' Kcal/mole/A**3')
      write (iout,20)  dedv_fd
20    format (' dE/dV (Finite Diff) :',12x,f15.6,' Kcal/mole/A**3')
   end if
!
!     compute analytical and finite-difference isotropic pressure
!
   pres_vir = prescon * (real(n,r_p)*gasconst*kelvin/&
      &volbox-dedv_vir)
   pres_fd = prescon * (real(n,r_p)*gasconst*kelvin/&
      &volbox-dedv_fd)
   if (rank.eq.0) then
      if (kelvin .eq. 0.0_re_p) then
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
end
!
subroutine pmonte (epot,temp)
   use atmlst
   use atmtyp
   use atomsMirror
   use bath
   use boxes
   use domdec
   use energi
   use group
   use math
   use mdstuf
   use molcul
   use moldyn
   use inform     ,only: mtc_nacc,deb_Path,deb_Energy,deb_Force&
      &,verbose
   use random_mod
   use units
   use usage
   use sizes      ,only: tinkerdebug
   use mpi
   implicit none
   integer i,j,k,iglob,ierr
   integer start,stop
   real(r_p) epot
   real(r_p) temp,term
   real(r_p) kt,expterm
   real(r_p) third,weigh
   real(r_p) step
   real(r_p) scale
   real(r_p) eold
   real(r_p) rnd6
   real(r_p) xcm,ycm,zcm
   real(r_p) vxcm,vycm,vzcm
   real(r_p) volold,cosine
   real(r_p) dpot,dpv,dkin
   real(r_p) xmove,ymove,zmove
   real(r_p) vxmove,vymove,vzmove
   real(r_p) xboxold,yboxold,zboxold
   real(t_p) alphaold,betaold,gammaold
   real(t_p) temp3(3,3)
   real(t_p) hbox(3,3)
   real(t_p) ascale(3,3)
   real(r_p), allocatable :: xold1(:)
   real(r_p), allocatable :: yold1(:)
   real(r_p), allocatable :: zold1(:)
   real(r_p), allocatable :: vold(:,:)
   real(t_p) valrand
   logical dotrial
   logical isotropic
   parameter(third = 1.0_re_p / 3.0_re_p)
   interface
      function energy ()
         real(r_p) energy
      end function
   end interface
!
!
!     decide whether to attempt a box size change at this step
!
   dotrial = .false.
   if (rank.eq.0) then
      valrand = random()
   end if
   call MPI_BCAST(valrand,1,MPI_TPREC,0,COMM_TINKER,ierr)
!
   if (valrand .lt. 1.0_ti_p/real(voltrial,t_p)) dotrial=.true.
!
!     set constants and decide on type of trial box size change
!
   if (dotrial) then
14    format(A,2x,2F12.6)
      if (deb_Path) write(*,14) ' __montecarlo barostat__'&
         &,valrand,1.0_ti_p/real(voltrial,t_p)
      isotropic = .true.
!
!     perform dynamic allocation of some local arrays
!
      allocate (xold1(n))
      allocate (yold1(n))
      allocate (zold1(n))
      allocate (vold(3,n))
!$acc data create(xold1,yold1,zold1,vold,eold,dpot) &
!$acc     present(x,y,z,v,glob,mass,molmass,kmol,use,imol &
!$acc            ,epot,temp) async

      !Broadcast potential
!$acc update host(epot,temp) async
!$acc wait
      call MPI_BCAST(epot,1,MPI_RPREC,0,COMM_TINKER,ierr)

      kt = gasconst * temp
      if (isothermal)  kt = gasconst * kelvin
!
!     save the system state prior to trial box size change
!
      xboxold  = xbox
      yboxold  = ybox
      zboxold  = zbox
      alphaold = alpha
      betaold  = beta
      gammaold = gamma
      volold   = volbox
      eold     = epot
!$acc parallel loop async
      do i = 1, nbloc
         iglob = glob(i)
         if (use(iglob)) then
            xold1(iglob)  = x(iglob)
            yold1(iglob)  = y(iglob)
            zold1(iglob)  = z(iglob)
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
         call MPI_BCAST(valrand,1,MPI_TPREC,0,COMM_TINKER,ierr)
         step   = volmove * (2.0_re_p*valrand-1.0_re_p)
         volbox = volbox + step
         scale  = (volbox/volold)**third
         xbox   = xbox * scale
         ybox   = ybox * scale
         zbox   = zbox * scale
!           print*," _scale_ ",scale
!$acc update device(volbox,xbox,ybox,zbox) async
         call lattice
         if (volscale .eq. 'MOLECULAR') then
            if (nproc.gt.1) then
               call molecule(.false.)
!$acc wait
            end if
            scale = scale - 1.0_re_p
!$acc parallel loop gang vector async
            do i = 1, nmoleloc
               iglob = molculeglob(i)
               xcm   = 0.0_re_p
               ycm   = 0.0_re_p
               zcm   = 0.0_re_p
               vxcm  = 0.0_re_p
               vycm  = 0.0_re_p
               vzcm  = 0.0_re_p
               start = imol(1,iglob)
               stop  = imol(2,iglob)
               do j = start, stop
                  k     = kmol(j)
                  weigh = mass(k)
                  xcm   =  xcm +   x(k)*weigh
                  ycm   =  ycm +   y(k)*weigh
                  zcm   =  zcm +   z(k)*weigh
                  vxcm  = vxcm + v(1,k)*weigh
                  vycm  = vycm + v(2,k)*weigh
                  vzcm  = vzcm + v(3,k)*weigh
               end do
               xmove  = scale *  xcm/molmass(iglob)
               ymove  = scale *  ycm/molmass(iglob)
               zmove  = scale *  zcm/molmass(iglob)
               vxmove = scale * vxcm/molmass(iglob)
               vymove = scale * vycm/molmass(iglob)
               vzmove = scale * vzcm/molmass(iglob)
               do j = start, stop
                  k = kmol(j)
                  if (use(k)) then
                     x(k)   =   x(k) +  xmove
                     y(k)   =   y(k) +  ymove
                     z(k)   =   z(k) +  zmove
                     v(1,k) = v(1,k) - vxmove
                     v(2,k) = v(2,k) - vymove
                     v(3,k) = v(3,k) - vzmove
                  end if
               end do
            end do
         else
!$acc parallel loop async
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
!              call ddpme3dnpt(scale,0)
         end if
         call reCast_position
      end if
!
!     get the potential energy and PV work changes for trial move
!
      epot = energy ()

!$acc update device(epot)
      if (nproc.gt.1) then
         call allreduceen(epot)
!$acc update host(epot)
      end if
      if (deb_Energy.or.deb_Force) then
         if(ranktot.eq.0) write(*,*) 'montecarlo energy'
         call info_energy(ranktot)
      end if

      dpot = epot - eold
      dpv = atmsph * (volbox-volold) / prescon
!
!     estimate the kinetic energy change as an ideal gas term
!
      if (volscale .eq. 'MOLECULAR') then
         dkin = real(nmol,r_p) * kt * log(volold/volbox)
      else
         dkin = real(nuse,r_p) * kt * log(volold/volbox)
      end if
!
!     compute the kinetic energy change from the actual velocities;
!     scale the kinetic energy change to match virial pressure
!
!        dkin = 0.0_re_p
!        do i = 1, n
!           if (use(i)) then
!              term = 0.5_re_p * mass(i) / convert
!              do j = 1, 3
!                 dkin = dkin + term*(v(j,i)**2-vold(j,i)**2)
!              end do
!           end if
!        end do
!        dkin = 0.907_re_p * dkin
!
!     acceptance ratio from Epot change, Ekin change and PV work
!
      term = -(dpot+dpv+dkin) / kt
      expterm = exp(term)
!
!     reject the step, and restore values prior to trial change
!
15    format(A,F20.6,D24.12,4F20.6)
      if (rank.eq.0) then
         valrand = random()
      end if
      call MPI_BCAST(valrand,1,MPI_TPREC,0,COMM_TINKER,ierr)
      if (valrand .gt. expterm) then
         if (rank.eq.0.and.tinkerdebug.gt.0)&
            &write(*,15) ' Reject montecarlo',valrand&
            &,expterm,epot,eold,dpot,dkin
         epot = eold
         xbox = xboxold
         ybox = yboxold
         zbox = zboxold
!$acc update device(xbox,ybox,zbox,epot)
         call lattice
!$acc parallel loop async
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
         call reCast_position
         goto 66
      else
         if (rank.eq.0.and.tinkerdebug.gt.0)&
            &write(*,15) ' Accept montecarlo',valrand&
            &,expterm,epot,eold,dpot,dkin
         if (rank.eq.0.and.verbose.and.tinkerdebug.eq.0)&
            &mtc_nacc = mtc_nacc + 1
!    &         write(*,*) 'Applied montecarlo barostat'
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
66    continue
!$acc end data
      deallocate (xold1)
      deallocate (yold1)
      deallocate (zold1)
      deallocate (vold)
   end if
end


subroutine dedvcalc()
   use atmlst
   use atmtyp
   use atomsMirror
   use bath
   use boxes
   use bound
   use domdec
   use energi
   use group
   use math
   use mdstuf
   use molcul
   use moldyn
   use inform     ,only: deb_Path,deb_Energy,deb_Force,verbose
   use random_mod
   use units
   use usage
   use mpi
   use virial
   implicit none
   integer i,iglob
   real(r_p) third
   real(r_p) delta,step,scale
   real(r_p) pres
   real(r_p) vold,xboxold
   real(r_p) yboxold,zboxold
   real(r_p) epos,eneg
   real(r_p), allocatable :: xoldloc(:)
   real(r_p), allocatable :: yoldloc(:)
   real(r_p), allocatable :: zoldloc(:)
   interface
      function energy ()
         real(r_p) energy
      end function
   end interface

   dedv=0.d0
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
!$acc data create(xoldloc,yoldloc,zoldloc,epos,eneg) &
!$acc     present(x,y,z,v,glob,mass,molmass,kmol,use,imol) async
!
!     store original box dimensions and coordinate values
!
      xboxold = xbox
      yboxold = ybox
      zboxold = zbox
      vold = volbox
!$acc parallel loop async
      do i = 1, nbloc
         iglob = glob(i)
         if (use(iglob)) then
            xoldloc(iglob) = x(iglob)
            yoldloc(iglob) = y(iglob)
            zoldloc(iglob) = z(iglob)
         end if
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
!$acc update device(volbox,xbox,ybox,zbox) async
      call lattice
!$acc parallel loop async
      do i = 1, nbloc
         iglob = glob(i)
         x(iglob) = xoldloc(iglob) * scale
         y(iglob) = yoldloc(iglob) * scale
         z(iglob) = zoldloc(iglob) * scale
      end do
      call reCast_position
!
!     compute potential energy for negative volume change
!
      eneg = energy ()
      if (nproc.gt.1) then
!$acc update device(eneg)
         call allreduceen(eneg)
!$acc update host(eneg)
      end if
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
!$acc update device(volbox,xbox,ybox,zbox) async
      call lattice
!$acc parallel loop async
      do i = 1, nbloc
         iglob    = glob(i)
         x(iglob) = xoldloc(iglob) * scale
         y(iglob) = yoldloc(iglob) * scale
         z(iglob) = zoldloc(iglob) * scale
      end do
      call reCast_position
!
!     compute potential energy for positive volume change
!
      epos = energy ()
      if (nproc.gt.1) then
!$acc update device(epos)
         call allreduceen(epos)
!$acc update host(epos)
      end if
!
!     restore original box dimensions and coordinate values
!
      xbox = xboxold
      ybox = yboxold
      zbox = zboxold
      volbox = vold
!$acc update device(volbox,xbox,ybox,zbox) async
      call lattice
!$acc parallel loop async
      do i = 1, nbloc
         iglob = glob(i)
         x(iglob) = xoldloc(iglob)
         y(iglob) = yoldloc(iglob)
         z(iglob) = zoldloc(iglob)
      end do
      call reCast_position
!$acc end data
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
!$acc serial async copyout(dedv) present(vir,volbox)
      dedv = (vir(1,1)+vir(2,2)+vir(3,3)) / (3.0d0*volbox)
!$acc end serial
   endif

   !write(*,*) dedv
!
!
!     compute analytical and finite-difference isotropic pressure
!
!      pres_vir = prescon * (dble(n)*gasconst*kelvin/volbox-dedv_vir)
!      pres= prescon * (dble(n)*gasconst*kelvin/volbox-dedv_fd)
!
!     for 4site water model
!
!      pres_vir = prescon * (dble(n-n/4)*gasconst*kelvin/volbox-dedv_vir)
!      pres= prescon * (dble(n-n/4)*gasconst*kelvin/volbox-dedv_fd)
!      write(*,*) 'pres=', pres, pres_vir,dble(n-n/4)
   return
end subroutine dedvcalc
!
!     subroutine rescale: rescale positions and speeds after a change of volume
!
subroutine rescale(istep)
   use atomsMirror
   use bath
   use boxes
   use domdec
   use moldyn
   use usage
   implicit none
   real(r_p) third,scale
   integer iglob,i,istep
   parameter(third=1.0/3.0)

   scale =  (extvol/extvolold)**third
!
!     modify the current periodic box dimension values
!
   xbox = xbox * scale
   ybox = ybox * scale
   zbox = zbox * scale
!$acc update device(xbox,ybox,zbox) async
!
!     propagate the new box dimensions to other lattice values
!
   call lattice

!$acc parallel loop present(x,y,z,v,use,glob) async
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
   call reCast_position
!
!     also rescale xbegproc, xendproc...
!
   call ddpme3dnpt(scale,istep)
end

!
!     propagate Volume with Langevin equation with a BAOAB propagator
!
subroutine initialize_langevin_piston()
   use bath
   use boxes
   use mpi
   use domdec
   implicit none
   integer ierr
   interface
      function maxwell (mass,temper)
         real(r_p) maxwell
         real(r_p) mass
         real(r_p) temper
      end function
   end interface

   extvol     = volbox
   extvolold  = volbox
   temppiston = kelvin
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
      call MPI_BCAST(vextvol,1,MPI_RPREC,0,COMM_TINKER,ierr)
   else
      call MPI_BCAST(vextbox,3,MPI_RPREC,0,COMM_TINKER,ierr)
   endif
   aextvol    = 0.0d0
   aextbox    = 0.0d0
end subroutine initialize_langevin_piston

subroutine rescale_box(istep,scale)
   use atoms
   use bath
   use boxes
   use domdec
   use moldyn
   use usage
   implicit none
   integer, intent(in) :: istep
   real(r_p), intent(in) :: scale(3)
   real*8 third
   integer iglob,i
!
!     modify the current periodic box dimension values
!
   xbox = xbox * scale(1)
   ybox = ybox * scale(2)
   zbox = zbox * scale(3)
!$acc update device(xbox,ybox,zbox) async
!
!     propagate the new box dimensions to other lattice values
!
   call lattice
!
!   also rescale xbegproc, xendproc...
!
   call ddpme3dnptaniso(scale,istep)
end subroutine rescale_box
