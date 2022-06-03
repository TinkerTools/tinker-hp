c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine mdrest  --  stop system translation & rotation  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "mdrest" finds and removes any translational or rotational
c     kinetic energy of the overall system center of mass
c
c
#include "tinker_precision.h"
      subroutine mdrest (istep)
      use atmtyp
      use atoms
      use bound
      use domdec
      use group
      use inform
      use iounit
      use mdstuf
      use moldyn
      use units
      use mpi
      implicit none
      integer i,j,k,istep,iglob,ierr
      real(t_p) etrans,erot
      real(t_p) weigh,totmass,eps
      real(t_p) xx,yy,zz,xy,xz,yz
      real(t_p) xtot,ytot,ztot
      real(t_p) xdel,ydel,zdel
      real(t_p) mang(3),vang(3)
      real(t_p) vtot(3),tensor(3,3)
      real(t_p), allocatable :: xcm(:)
      real(t_p), allocatable :: ycm(:)
      real(t_p), allocatable :: zcm(:)
c
c
c     check steps between center of mass motion removal
c
      if (.not.dorest)  return
      if (mod(istep,irest) .ne. 0)  return
      if (deb_Path) write(*,*) ' mdrestgpu'
!$acc wait
!$acc update host(v)
c
c     zero out the total mass and overall linear velocity
c
      totmass = 0.0_ti_p
      do j = 1, 3
         vtot(j) = 0.0_ti_p
      end do
c
c     compute linear velocity of the system center of mass
c
      do i = 1, nloc
         iglob = glob(i)
         weigh = mass(iglob)
         totmass = totmass + weigh
         do j = 1, 3
            vtot(j) = vtot(j) + v(j,iglob)*weigh
         end do
      end do
      call MPI_ALLREDUCE(MPI_IN_PLACE,vtot,3,MPI_TPREC,MPI_SUM,
     $   COMM_TINKER,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,totmass,1,MPI_TPREC,MPI_SUM,
     $   COMM_TINKER,ierr)
c
c     compute translational kinetic energy of overall system
c
      etrans = 0.0_ti_p
      do j = 1, 3
         vtot(j) = vtot(j) / totmass
         etrans = etrans + vtot(j)**2
      end do
      etrans = 0.5_ti_p * etrans * totmass / convert
c
c     find the center of mass coordinates of the overall system
c
      if (.not. use_bounds) then
         xtot = 0.0_ti_p
         ytot = 0.0_ti_p
         ztot = 0.0_ti_p
            do i = 1, nloc
               iglob = glob(i)
               weigh = mass(iglob)
               xtot = xtot + x(iglob)*weigh
               ytot = ytot + y(iglob)*weigh
               ztot = ztot + z(iglob)*weigh
            end do
         call MPI_ALLREDUCE(MPI_IN_PLACE,xtot,1,MPI_TPREC,MPI_SUM,
     $      COMM_TINKER,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,ytot,1,MPI_TPREC,MPI_SUM,
     $      COMM_TINKER,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,ztot,1,MPI_TPREC,MPI_SUM,
     $      COMM_TINKER,ierr)
         xtot = xtot / totmass
         ytot = ytot / totmass
         ztot = ztot / totmass
c
c     compute the angular momentum of the overall system
c
         do j = 1, 3
            mang(j) = 0.0_ti_p
         end do
         do i = 1, nloc
               iglob = glob(i)
               weigh = mass(iglob)
               mang(1) = mang(1) + (y(iglob)*v(3,iglob)-z(iglob)*
     $          v(2,iglob))*weigh
               mang(2) = mang(2) + (z(iglob)*v(1,iglob)-x(iglob)*
     $          v(3,iglob))*weigh
               mang(3) = mang(3) + (x(iglob)*v(2,iglob)-y(iglob)*
     $          v(1,iglob))*weigh
         end do
         call MPI_ALLREDUCE(MPI_IN_PLACE,mang,3,MPI_TPREC,MPI_SUM,
     $      COMM_TINKER,ierr)
c
         mang(1) = mang(1) - (ytot*vtot(3)-ztot*vtot(2))*totmass
         mang(2) = mang(2) - (ztot*vtot(1)-xtot*vtot(3))*totmass
         mang(3) = mang(3) - (xtot*vtot(2)-ytot*vtot(1))*totmass
c
c     calculate the moment of inertia tensor
c
         xx = 0.0_ti_p
         xy = 0.0_ti_p
         xz = 0.0_ti_p
         yy = 0.0_ti_p
         yz = 0.0_ti_p
         zz = 0.0_ti_p
         do i = 1, nloc
            iglob = glob(i)
            weigh = mass(iglob)
            xdel = x(iglob) - xtot
            ydel = y(iglob) - ytot
            zdel = z(iglob) - ztot
            xx = xx + xdel*xdel*weigh
            xy = xy + xdel*ydel*weigh
            xz = xz + xdel*zdel*weigh
            yy = yy + ydel*ydel*weigh
            yz = yz + ydel*zdel*weigh
            zz = zz + zdel*zdel*weigh
         end do
         call MPI_ALLREDUCE(MPI_IN_PLACE,xx,1,MPI_TPREC,MPI_SUM,
     $      COMM_TINKER,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,xy,1,MPI_TPREC,MPI_SUM,
     $      COMM_TINKER,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,xz,1,MPI_TPREC,MPI_SUM,
     $      COMM_TINKER,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,yy,1,MPI_TPREC,MPI_SUM,
     $      COMM_TINKER,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,yz,1,MPI_TPREC,MPI_SUM,
     $      COMM_TINKER,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,zz,1,MPI_TPREC,MPI_SUM,
     $      COMM_TINKER,ierr)
         tensor(1,1) = yy + zz
         tensor(2,1) = -xy
         tensor(3,1) = -xz
         tensor(1,2) = -xy
         tensor(2,2) = xx + zz
         tensor(3,2) = -yz
         tensor(1,3) = -xz
         tensor(2,3) = -yz
         tensor(3,3) = xx + yy
c
c     fix to avoid singularity for one- or two-body systems
c
         if (n .le. 2) then
            eps = 0.000001_ti_p
            tensor(1,1) = tensor(1,1) + eps
            tensor(2,2) = tensor(2,2) + eps
            tensor(3,3) = tensor(3,3) + eps
         end if
c
c     diagonalize the moment of inertia tensor
c
         call invert (3,tensor)
c
c     compute angular velocity and rotational kinetic energy
c
         erot = 0.0_ti_p
         do i = 1, 3
            vang(i) = 0.0_ti_p
            do j = 1, 3
               vang(i) = vang(i) + tensor(i,j)*mang(j)
            end do
            erot = erot + vang(i)*mang(i)
         end do
         erot = 0.5_ti_p * erot / convert
      end if
c
c     eliminate any translation of the overall system
c
      do i = 1, nloc
         iglob = glob(i)
         do j = 1, 3
            v(j,iglob) = v(j,iglob) - vtot(j)
         end do
      end do
c
c     print the translational velocity of the overall system
c
      if (debug) then
         write (iout,10)  (vtot(i),i=1,3),etrans
   10    format (' System Linear Velocity :  ',3d12.2,
     &           /,' Translational Kinetic Energy :',10x,f12.4,
     &              ' Kcal/mole')
      end if
!$acc update device(v)
      return
      end
c
c
c
      module mdrest_inl
      real(r_p) etrans,erot,totmass
      real(r_p) vtot(3),vtot_deb(3)
      real(r_p) vtot1,vtot2,vtot3
      logical:: mdrest_l=.true.
      end module

      subroutine mdrestgpu (istep)
      use atmtyp
      use atoms
      use bound
      use domdec
      use group
      use inform
      use iounit
      use mdstuf
      use mdrest_inl
      use moldyn
      use units
      use mpi
      use timestat
      implicit none
      integer i,j,k,istep,iglob,ierr
      real(r_p) weigh
      real(t_p) eps
      real(t_p) xx,yy,zz,xy,xz,yz
      real(t_p) xtot,ytot,ztot
      real(t_p) xdel,ydel,zdel
      real(t_p) mang(3),vang(3)
      real(t_p) tensor(3,3)
      real(t_p) mang1,mang2,mang3
      real(t_p) vang1,vang2,vang3
      real(t_p), allocatable :: xcm(:)
      real(t_p), allocatable :: ycm(:)
      real(t_p), allocatable :: zcm(:)
c
c
c     check steps between center of mass motion removal
c
      if (.not.dorest)  return
      if (mod(istep,irest) .ne. 0)  return
      if (deb_Path) write(*,*) ' mdrestgpu'
      call timer_enter(timer_other)

      if (mdrest_l) then
!$acc enter data create(vtot1,vtot2,vtot3,vtot,etrans,totmass)
         mdrest_l=.false.
      end if
c
c     zero out the total mass and overall linear velocity
c
!$acc host_data use_device(vtot1,vtot2,vtot3,vtot,etrans,totmass
!$acc&         ,glob,mass,v)

!$acc serial async deviceptr(totmass,vtot1,vtot2,vtot3,etrans)
      totmass = 0.0_re_p
      vtot1   = 0.0_re_p
      vtot2   = 0.0_re_p
      vtot3   = 0.0_re_p
      etrans  = 0.0_re_p
!$acc end serial
c
c     compute linear velocity of the system center of mass
c
!$acc parallel loop gang vector collapse(2) async
!$acc&         deviceptr(vtot1,vtot2,vtot3,totmass,glob,mass,v)
!$acc&         reduction(vtot1,vtot2,vtot3,totmass)
      do i = 1, nloc
         do j = 1, 3
            iglob = glob(i)
            weigh = mass(iglob)
            if (j.eq.1) then
               totmass = totmass + weigh
               vtot1 = vtot1 + v(1,iglob)*weigh
            else if (j.eq.2) then
               vtot2 = vtot2 + v(2,iglob)*weigh
            else
               vtot3 = vtot3 + v(3,iglob)*weigh
            end if
         end do
      end do
c
      if (nproc.eq.1) then
!$acc serial async deviceptr(vtot,vtot1,vtot2,vtot3,totmass,etrans)
         !compute translational kinetic energy of overall system
         vtot(1) = vtot1 / totmass
         vtot(2) = vtot2 / totmass
         vtot(3) = vtot3 / totmass
         etrans  = etrans + vtot(1)**2 + vtot(2)**2 + vtot(3)**2
         etrans  = 0.5_re_p * etrans * totmass / convert
!$acc end serial
      else
         !compute translational kinetic energy of overall system
!$acc serial async deviceptr(vtot,vtot1,vtot2,vtot3)
         vtot(1) = vtot1
         vtot(2) = vtot2
         vtot(3) = vtot3
!$acc end serial
!$acc wait
         call MPI_ALLREDUCE(MPI_IN_PLACE,vtot,3,MPI_RPREC,MPI_SUM,
     $        COMM_TINKER,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,totmass,1,MPI_RPREC,MPI_SUM,
     $        COMM_TINKER,ierr)
 
 
!$acc serial async deviceptr(vtot,vtot1,vtot2,vtot3,totmass,etrans)
         vtot(1) = vtot(1) / totmass
         vtot(2) = vtot(2) / totmass
         vtot(3) = vtot(3) / totmass
         etrans  = etrans + vtot(1)**2 + vtot(2)**2 + vtot(3)**2
         etrans  = 0.5_re_p * etrans * totmass / convert
!$acc end serial
      end if
c
c     find the center of mass coordinates of the overall system
c
      if (.not. use_bounds) then
         xtot = 0.0_ti_p
         ytot = 0.0_ti_p
         ztot = 0.0_ti_p
            do i = 1, nloc
               iglob = glob(i)
               weigh = mass(iglob)
               xtot = xtot + x(iglob)*weigh
               ytot = ytot + y(iglob)*weigh
               ztot = ztot + z(iglob)*weigh
            end do
         call MPI_ALLREDUCE(MPI_IN_PLACE,xtot,1,MPI_TPREC,MPI_SUM,
     $      COMM_TINKER,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,ytot,1,MPI_TPREC,MPI_SUM,
     $      COMM_TINKER,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,ztot,1,MPI_TPREC,MPI_SUM,
     $      COMM_TINKER,ierr)
         xtot = xtot / totmass
         ytot = ytot / totmass
         ztot = ztot / totmass
c
c     compute the angular momentum of the overall system
c
         do j = 1, 3
            mang(j) = 0.0_ti_p
         end do
         do i = 1, nloc
               iglob = glob(i)
               weigh = mass(iglob)
               mang(1) = mang(1) + (y(iglob)*v(3,iglob)-z(iglob)*
     $          v(2,iglob))*weigh
               mang(2) = mang(2) + (z(iglob)*v(1,iglob)-x(iglob)*
     $          v(3,iglob))*weigh
               mang(3) = mang(3) + (x(iglob)*v(2,iglob)-y(iglob)*
     $          v(1,iglob))*weigh
         end do
         call MPI_ALLREDUCE(MPI_IN_PLACE,mang,3,MPI_TPREC,MPI_SUM,
     $      COMM_TINKER,ierr)
c
         mang(1) = mang(1) - (ytot*vtot(3)-ztot*vtot(2))*totmass
         mang(2) = mang(2) - (ztot*vtot(1)-xtot*vtot(3))*totmass
         mang(3) = mang(3) - (xtot*vtot(2)-ytot*vtot(1))*totmass
c
c     calculate the moment of inertia tensor
c
         xx = 0.0_ti_p
         xy = 0.0_ti_p
         xz = 0.0_ti_p
         yy = 0.0_ti_p
         yz = 0.0_ti_p
         zz = 0.0_ti_p
         do i = 1, nloc
            iglob = glob(i)
            weigh = mass(iglob)
            xdel = x(iglob) - xtot
            ydel = y(iglob) - ytot
            zdel = z(iglob) - ztot
            xx = xx + xdel*xdel*weigh
            xy = xy + xdel*ydel*weigh
            xz = xz + xdel*zdel*weigh
            yy = yy + ydel*ydel*weigh
            yz = yz + ydel*zdel*weigh
            zz = zz + zdel*zdel*weigh
         end do
         call MPI_ALLREDUCE(MPI_IN_PLACE,xx,1,MPI_TPREC,MPI_SUM,
     $      COMM_TINKER,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,xy,1,MPI_TPREC,MPI_SUM,
     $      COMM_TINKER,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,xz,1,MPI_TPREC,MPI_SUM,
     $      COMM_TINKER,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,yy,1,MPI_TPREC,MPI_SUM,
     $      COMM_TINKER,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,yz,1,MPI_TPREC,MPI_SUM,
     $      COMM_TINKER,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,zz,1,MPI_TPREC,MPI_SUM,
     $      COMM_TINKER,ierr)
         tensor(1,1) = yy + zz
         tensor(2,1) = -xy
         tensor(3,1) = -xz
         tensor(1,2) = -xy
         tensor(2,2) = xx + zz
         tensor(3,2) = -yz
         tensor(1,3) = -xz
         tensor(2,3) = -yz
         tensor(3,3) = xx + yy
c
c     fix to avoid singularity for one- or two-body systems
c
         if (n .le. 2) then
            eps = 0.000001_ti_p
            tensor(1,1) = tensor(1,1) + eps
            tensor(2,2) = tensor(2,2) + eps
            tensor(3,3) = tensor(3,3) + eps
         end if
c
c     diagonalize the moment of inertia tensor
c
         call invert (3,tensor)
c
c     compute angular velocity and rotational kinetic energy
c
         erot = 0.0_re_p
         do i = 1, 3
            vang(i) = 0.0_ti_p
            do j = 1, 3
               vang(i) = vang(i) + tensor(i,j)*mang(j)
            end do
            erot = erot + vang(i)*mang(i)
         end do
         erot = 0.5_re_p * erot / convert
      end if
c
c     eliminate any translation of the overall system
c
!$acc parallel loop collapse(2) async deviceptr(glob,v,vtot)
      do i = 1, nloc
         do j = 1, 3
            iglob = glob(i)
            v(j,iglob) = v(j,iglob) - vtot(j)
         end do
      end do
!$acc end host_data
c
c     print the translational velocity of the overall system
c
      if (debug) then
!$acc update host(vtot,etrans) async
!$acc wait
         vtot_deb(:) = vtot(:)
         write (iout,10)  (vtot_deb(i),i=1,3),etrans
   10    format (' System Linear Velocity :  ',3d12.2,
     &           /,' Translational Kinetic Energy :',10x,f12.4,
     &              ' Kcal/mole')
      end if
c
      call timer_exit(timer_other)
      end
