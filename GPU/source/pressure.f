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
c
c     only necessary if periodic boundaries are in use
c     and isobaric simulation
c
      if (.not.(use_bounds.and.use_virial))  return
c
c     calculate the stress tensor for anisotropic systems
c
      factor = prescon / volbox
c!$acc parallel loop collapse(2) default(present) async
c      do i = 1, 3
c         do j = 1, 3
c            stress(j,i) = factor * (2*ekin(j,i)-vir(j,i))
c         end do
c      end do
c
c     set isotropic pressure to the average of tensor diagonal
c
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
c
c     use the Berendsen barostat method
c
      if (isobaric) then
         if (barostat.eq.'BERENDSEN') call pscale(dt,pres,stress,istep)
      end if
      end
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
      real(r_p) epot,temp
c
c     only necessary if periodic boundaries are in use 
c     and isobaric simulation
c
      if (.not.(use_bounds.and.isobaric))  return

      ! --- Monte Carlo barostat method --- !
      if (barostat.eq.'MONTECARLO') call pmonte(epot,temp)

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
      real(r_p) scale,third
      real(r_p) xcm,xmove
      real(r_p) ycm,ymove
      real(r_p) zcm,zmove
      real(r_p) stress(3,3)
      real(t_p) temp(3,3)
      real(t_p) hbox(3,3)
      real(t_p) ascale(3,3)
      parameter(third = 1.0_re_p / 3.0_re_p)
c
c
c     find the isotropic scale factor for constant pressure
c
c      if (.not. anisotrop) then
!$acc wait
!$acc update host(pres)
         scale = 1.0_re_p
         scale = (1.0_re_p + (dt*compress/taupres)*(pres-atmsph))**third
c
c     modify the current periodic box dimension values
c
         xbox = xbox * scale
         ybox = ybox * scale
         zbox = zbox * scale
!$acc update device(xbox,ybox,zbox) async
c
c     propagate the new box dimensions to other lattice values
c
         call lattice

         if (deb_Path) then
 13      format(A,8F16.6)
         print 13,'pscale',scale,dt,compress,pres,atmsph,xbox,ybox,zbox
         end if
c
c     couple to pressure bath via atom scaling in Cartesian space
c
!$acc parallel loop async
!$acc&         present(glob,x,y,z)
         do i = 1, nbloc
            iglob = glob(i)
            if (use(iglob)) then
               x(iglob) = x(iglob) * scale
               y(iglob) = y(iglob) * scale
               z(iglob) = z(iglob) * scale
            end if
         end do
         call reCast_position
c
c   also rescale xbegproc, xendproc...
c
         call ddpme3dnpt(scale,istep)
c
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
c
c
c     set relative volume change for finite-differences
c
      if (.not. use_bounds)  return
      delta = 0.000001_re_p
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
      third = 1.0_re_p / 3.0_re_p
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
      scale = (volbox/vold)**third
!$acc update device(volbox) async
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
!$acc update device(xbox,ybox,zbox) async
      call lattice
      do i = 1, nbloc
         iglob = glob(i)
         x(iglob) = xoldloc(iglob)
         y(iglob) = yoldloc(iglob)
         z(iglob) = zoldloc(iglob)
      end do
!$acc update device(x(:),y(:),z(:))
c
c     perform deallocation of some local arrays
c
      deallocate (xoldloc)
      deallocate (yoldloc)
      deallocate (zoldloc)
c
c     get virial and finite difference values of dE/dV
c
      dedv_vir = (vir(1,1)+vir(2,2)+vir(3,3)) / (3.0_re_p*volbox)
      dedv_fd = real(epos-eneg,r_p) / (2.0_re_p*delta*volbox)
      if (rank.eq.0) then
        write (iout,10)  dedv_vir
   10   format (/,' dE/dV (Virial-based) :',11x,f15.6,' Kcal/mole/A**3')
        write (iout,20)  dedv_fd
   20   format (' dE/dV (Finite Diff) :',12x,f15.6,' Kcal/mole/A**3')
      end if
c
c     compute analytical and finite-difference isotropic pressure
c
      pres_vir = prescon * (real(n,r_p)*gasconst*kelvin/
     &                                     volbox-dedv_vir)
      pres_fd = prescon * (real(n,r_p)*gasconst*kelvin/
     &                                    volbox-dedv_fd)
      if (rank.eq.0) then
      if (kelvin .eq. 0.0_re_p) then
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
      end
c
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
      use inform     ,only: mtc_nacc,deb_Path,deb_Energy,deb_Force
     &               ,verbose
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
c
c
c     decide whether to attempt a box size change at this step
c
      dotrial = .false.
      if (rank.eq.0) then
        valrand = random()
      end if
      call MPI_BCAST(valrand,1,MPI_TPREC,0,COMM_TINKER,ierr)
c
      if (valrand .lt. 1.0_ti_p/real(voltrial,t_p)) dotrial=.true.
c
c     set constants and decide on type of trial box size change
c
      if (dotrial) then
 14      format(A,2x,2F12.6)
         if (deb_Path) write(*,14) ' __montecarlo barostat__'
     &      ,valrand,1.0_ti_p/real(voltrial,t_p)
         isotropic = .true.
c
c     perform dynamic allocation of some local arrays
c
         allocate (xold1(n))
         allocate (yold1(n))
         allocate (zold1(n))
         allocate (vold(3,n))
!$acc data create(xold1,yold1,zold1,vold,eold,dpot)
!$acc&     present(x,y,z,v,glob,mass,molmass,kmol,use,imol
!$acc&            ,epot,temp) async

         !Broadcast potential
!$acc update host(epot,temp) async
!$acc wait
         call MPI_BCAST(epot,1,MPI_RPREC,0,COMM_TINKER,ierr)

         kt = gasconst * temp
         if (isothermal)  kt = gasconst * kelvin
c
c     save the system state prior to trial box size change
c
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
c
c     for the isotropic case, change the lattice lengths uniformly
c
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
c           print*," _scale_ ",scale
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
c              call ddpme3dnpt(scale,0)
            end if
            call reCast_position
         end if
c
c     get the potential energy and PV work changes for trial move
c
         epot = energy ()

!$acc update device(epot)
         if (nproc.gt.1) then
            call allreduceen(epot)
!$acc update host(epot)
         end if
         if (deb_Energy.or.deb_Force) then
            if(rank.eq.0) write(*,*) 'montecarlo energy'
            call info_energy(rank)
         end if

         dpot = epot - eold
         dpv = atmsph * (volbox-volold) / prescon
c
c     estimate the kinetic energy change as an ideal gas term
c
         if (volscale .eq. 'MOLECULAR') then
           dkin = real(nmol,r_p) * kt * log(volold/volbox)
         else
           dkin = real(nuse,r_p) * kt * log(volold/volbox)
         end if
c
c     compute the kinetic energy change from the actual velocities;
c     scale the kinetic energy change to match virial pressure
c
c        dkin = 0.0_re_p
c        do i = 1, n
c           if (use(i)) then
c              term = 0.5_re_p * mass(i) / convert
c              do j = 1, 3
c                 dkin = dkin + term*(v(j,i)**2-vold(j,i)**2)
c              end do
c           end if
c        end do
c        dkin = 0.907_re_p * dkin
c
c     acceptance ratio from Epot change, Ekin change and PV work
c
         term = -(dpot+dpv+dkin) / kt
         expterm = exp(term)
c
c     reject the step, and restore values prior to trial change
c
 15      format(A,F20.6,D24.12,4F20.6)
         if (rank.eq.0) then
           valrand = random()
         end if
         call MPI_BCAST(valrand,1,MPI_TPREC,0,COMM_TINKER,ierr)
         if (valrand .gt. expterm) then
            if (rank.eq.0.and.tinkerdebug.gt.0)
     &      write(*,15) ' Reject montecarlo',valrand
     &                 ,expterm,epot,eold,dpot,dkin
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
            if (rank.eq.0.and.tinkerdebug.gt.0)
     &      write(*,15) ' Accept montecarlo',valrand
     &                 ,expterm,epot,eold,dpot,dkin
            if (rank.eq.0.and.verbose.and.tinkerdebug.eq.0)
     &         mtc_nacc = mtc_nacc + 1
c    &         write(*,*) 'Applied montecarlo barostat'
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
  66     continue
!$acc end data
         deallocate (xold1)
         deallocate (yold1)
         deallocate (zold1)
         deallocate (vold)
      end if
      end
c
c     subroutine rescale: rescale positions and speeds after a change of volume
c
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
c
c     modify the current periodic box dimension values
c
      xbox = xbox * scale
      ybox = ybox * scale
      zbox = zbox * scale
!$acc update device(xbox,ybox,zbox) async
c
c     propagate the new box dimensions to other lattice values
c
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
c
c     also rescale xbegproc, xendproc...
c
      call ddpme3dnpt(scale,istep)
      end

c
c     propagate Volume with Langevin equation with a BAOAB propagator
c
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
      temppiston = 0.0
      aextvol    = 0.0
      if ( rank.eq.0 ) vextvol  = maxwell(masspiston,kelvin)
      if (nproc.gt.1 ) 
     &   call MPI_BCAST(vextvol,1,MPI_RPREC,0,COMM_TINKER,ierr)
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
      real(r_p), intent(in) :: scale
      real*8 third
      integer iglob,i
c
c     modify the current periodic box dimension values
c
      xbox = xbox * scale
      ybox = ybox * scale
      zbox = zbox * scale
!$acc update device(xbox,ybox,zbox) async
c
c     propagate the new box dimensions to other lattice values
c
      call lattice
c
c   also rescale xbegproc, xendproc...
c
      call ddpme3dnpt(scale,istep)
      end subroutine rescale_box

      subroutine pressure_iso(eksum,dedv,pres)
      use boxes
      use units
      implicit none
      real(r_p), intent(inout) :: pres
      real(r_p), intent(in) :: dedv,eksum
!$acc serial async present(pres,dedv,eksum,volbox)
      pres = prescon*( -dedv + 2.0*eksum
     &                       /(3.0*volbox) )
!$acc end serial
      end subroutine pressure_iso
