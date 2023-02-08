c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     "empole1" : driver for calculation of the multipole and dipole polarization
c     energy and derivatives with respect to Cartesian coordinates
c
c
#include "tinker_macro.h"
      subroutine empole1
      use domdec
      use energi
      use potent
      use mpi
      implicit none
c
c     choose the method for summing over multipole interactions
c
      if (use_lambdadyn) then
        call elambdampole1c
      else
        call empole1c
      end if
c
c     zero out energy and derivative terms which are not in use
c
      if (.not. use_mpole) then
         em = 0.0d0
      end if
      if (.not. use_polar) then
         ep = 0.0d0
      end if
      return
      end
c
c
c     "elambdampole1c" calculates the multipole energy and derivatives
c     with respect to Cartesian coordinates using particle mesh Ewald
c     summation and a neighbor list during lambda dynamics
c
c
      subroutine elambdampole1c
      use sizes
      use atmlst
      use atoms
      use boxes
      use chgpot
      use deriv
      use domdec
      use energi
      use ewald
      use math
      use mpole
      use mutant
      use potent
      use timestat
      use virial
      use mpi
      implicit none
      integer i,ii
      integer iipole,iglob,iloc,ierr
      real*8 e,f
      real*8 term,fterm
      real*8 cii,dii,qii
      real*8 xi,yi,zi
      real*8 xd,yd,zd
      real*8 xdtemp,ydtemp,zdtemp
      real*8 qtemp,muxtemp,muytemp,muztemp
      real*8 xq,yq,zq
      real*8 xv,yv,zv,vterm
      real*8 ci,dix,diy,diz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
      real*8 xdfield,ydfield
      real*8 zdfield
      real*8 fx,fy,fz
      real*8 vxx,vyy,vzz
      real*8 vxy,vxz,vyz
      real*8 trq(3),frcx(3)
      real*8 frcy(3),frcz(3)
      real*8 time0,time1
      real*8 citemp,dixtemp,diytemp,diztemp
      real*8 qixxtemp,qixytemp,qixztemp
      real*8 qiyytemp,qiyztemp,qizztemp
      real*8 elambdatemp
      real*8, allocatable :: delambdarec0(:,:),delambdarec1(:,:)
      real*8 :: elambdarec0,elambdarec1
      real*8, allocatable :: pot(:)
      real*8, allocatable :: decfx(:)
      real*8, allocatable :: decfy(:)
      real*8, allocatable :: decfz(:)
c
      allocate (delambdarec0(3,nlocrec2))
      allocate (delambdarec1(3,nlocrec2))
      elambdatemp = elambda  
c
c     zero out the atomic multipole energy and derivatives
c
      em = 0.0d0
      dem = 0d0
      if (npole .eq. 0)  return
      delambdae = 0d0

      aewald = aeewald
c
c     set the energy unit conversion factor
c
      f = electric / dielec
c
c     check the sign of multipole components at chiral sites
c
      call chkpole(.false.)
c
c     rotate the multipole components into the global frame
c
      call rotpole
c
c     compute the reciprocal space part of the Ewald summation
c
      if ((.not.(use_pmecore)).or.(use_pmecore).and.(rank.gt.ndir-1))
     $   then
        if (use_mrec) then 
          time0 = mpi_wtime()
c
c         the reciprocal part is interpolated between 0 and 1
c
          elambda = 0d0
          call MPI_BARRIER(hostcomm,ierr)
          if (hostrank.eq.0) call altelec
          call MPI_BARRIER(hostcomm,ierr)
          call rotpole
          em = 0d0
          demrec = 0d0
          if (elambda.lt.1d0) then
            call emrecip1
          end if
          elambdarec0  = em
          delambdarec0 = demrec

          elambda = 1d0
          call MPI_BARRIER(hostcomm,ierr)
          if (hostrank.eq.0) call altelec
          call MPI_BARRIER(hostcomm,ierr)
          call rotpole
          em = 0d0
          demrec = 0d0
          if (elambda.gt.0d0) then
            call emrecip1
          end if
          elambdarec1  = em
          delambdarec1 = demrec

          elambda = elambdatemp 
          em = (1-elambda)*elambdarec0 + elambda*elambdarec1
          demrec = (1-elambda)*delambdarec0+elambda*delambdarec1
          delambdae = delambdae + elambdarec1-elambdarec0
c
c         reset lambda to initial value
c
          call MPI_BARRIER(hostcomm,ierr)
          if (hostrank.eq.0) call altelec
          call MPI_BARRIER(hostcomm,ierr)
          call rotpole
          time1 = mpi_wtime()
          timerec = timerec + time1 - time0
        end if
      end if
c
c     compute the real space part of the Ewald summation
c
      if ((.not.(use_pmecore)).or.(use_pmecore).and.(rank.le.ndir-1))
     $   then
        if (use_mreal) then
          time0 = mpi_wtime()
          call emreal1c
          time1 = mpi_wtime()
          timereal = timereal + time1 - time0
        end if
c
        if (use_mself) then
c
c     compute the Ewald self-energy term over all the atoms
c
c
c     perform dynamic allocation of some local arrays
c
          allocate (pot(nbloc))
          allocate (decfx(nbloc))
          allocate (decfy(nbloc))
          allocate (decfz(nbloc))
          pot = 0d0
          term = 2.0d0 * aewald * aewald
          fterm = -f * aewald / sqrtpi
          do i = 1, npoleloc
             iipole = poleglob(i)
             iglob = ipole(iipole)
             iloc = loc(iglob)
             ci = rpole(1,iipole)
             dix = rpole(2,iipole)
             diy = rpole(3,iipole)
             diz = rpole(4,iipole)
             qixx = rpole(5,iipole)
             qixy = rpole(6,iipole)
             qixz = rpole(7,iipole)
             qiyy = rpole(9,iipole)
             qiyz = rpole(10,iipole)
             qizz = rpole(13,iipole)
             cii = ci*ci
             dii = dix*dix + diy*diy + diz*diz
             qii = 2.0d0*(qixy*qixy+qixz*qixz+qiyz*qiyz)
     &                + qixx*qixx + qiyy*qiyy + qizz*qizz
             e = fterm * (cii + term*(dii/3.0d0+2.0d0*term*qii/5.0d0))
             em = em + e
             pot(iloc) = 2.0d0 * fterm * ci
             if (mut(iglob).and.elambda.gt.0) then
               citemp = ci/elambda
               dixtemp = rpole(2,iipole)/elambda
               diytemp = rpole(3,iipole)/elambda
               diztemp = rpole(4,iipole)/elambda
               qixxtemp = rpole(5,iipole)/elambda
               qixytemp = rpole(6,iipole)/elambda
               qixztemp = rpole(7,iipole)/elambda
               qiyytemp = rpole(9,iipole)/elambda
               qiyztemp = rpole(10,iipole)/elambda
               qizztemp = rpole(13,iipole)/elambda
               delambdae = delambdae + fterm*2d0*elambda*citemp**2
               delambdae = delambdae + fterm*term*2d0*elambda*
     $           (dixtemp**2 + diytemp**2 + diztemp**2)/3d0
               delambdae = delambdae + fterm*term**2*8d0*elambda*
     $           (qixytemp**2 + qixztemp**2 + qiyztemp**2)/5d0
               delambdae = delambdae + fterm*term**2*4d0*elambda*
     $           (qixxtemp**2 + qiyytemp**2 + qizztemp**2)/5d0
             end if
          end do
c
c     modify gradient and virial for charge flux self-energy
c
          if (use_chgflx) then
             call commpot(pot,1)
             call dcflux (pot,decfx,decfy,decfz)
             do i = 1, npolebloc
                iipole = poleglob(i)
                iglob = ipole(iipole)
                iloc = loc(iglob)
                xi = x(iglob)
                yi = y(iglob)
                zi = z(iglob)
                fx = decfx(iloc)
                fy = decfy(iloc)
                fz = decfz(iloc)
                dem(1,iloc) = dem(1,iloc) + fx
                dem(2,iloc) = dem(2,iloc) + fy
                dem(3,iloc) = dem(3,iloc) + fz
                vxx = xi * fx
                vxy = yi * fx
                vxz = zi * fx
                vyy = yi * fy
                vyz = zi * fy
                vzz = zi * fz
                vir(1,1) = vir(1,1) + vxx
                vir(2,1) = vir(2,1) + vxy
                vir(3,1) = vir(3,1) + vxz
                vir(1,2) = vir(1,2) + vxy
                vir(2,2) = vir(2,2) + vyy
                vir(3,2) = vir(3,2) + vyz
                vir(1,3) = vir(1,3) + vxz
                vir(2,3) = vir(2,3) + vyz
                vir(3,3) = vir(3,3) + vzz
             end do
          end if
c
c     perform deallocation of some local arrays
c
          deallocate (pot)
          deallocate (decfx)
          deallocate (decfy)
          deallocate (decfz)
c
c       compute the cell dipole boundary correction term
c
          if (boundary .eq. 'VACUUM') then
             xd = 0.0d0
             yd = 0.0d0
             zd = 0.0d0
             xdtemp = 0.0d0
             ydtemp = 0.0d0
             zdtemp = 0.0d0
             do i = 1, npoleloc
                iipole = poleglob(i)
                iglob = ipole(iipole)
                xd = xd + rpole(2,iipole) + rpole(1,iipole)*x(iglob)
                yd = yd + rpole(3,iipole) + rpole(1,iipole)*y(iglob)
                zd = zd + rpole(4,iipole) + rpole(1,iipole)*z(iglob)
                if (mut(iglob)) then
                  qtemp = rpole(1,iipole)/elambda
                  muxtemp = rpole(2,iipole)/elambda
                  muytemp = rpole(3,iipole)/elambda
                  muztemp = rpole(4,iipole)/elambda
                  xdtemp = xdtemp + muxtemp + qtemp*x(iglob)
                  ydtemp = ydtemp + muytemp + qtemp*y(iglob)
                  zdtemp = zdtemp + muztemp + qtemp*z(iglob)
                end if
             end do
             call MPI_ALLREDUCE(MPI_IN_PLACE,xd,1,MPI_REAL8,MPI_SUM,
     $          COMM_TINKER,ierr)
             call MPI_ALLREDUCE(MPI_IN_PLACE,yd,1,MPI_REAL8,MPI_SUM,
     $          COMM_TINKER,ierr)
             call MPI_ALLREDUCE(MPI_IN_PLACE,zd,1,MPI_REAL8,MPI_SUM,
     $          COMM_TINKER,ierr)
             call MPI_ALLREDUCE(MPI_IN_PLACE,xdtemp,1,MPI_REAL8,MPI_SUM,
     $          COMM_TINKER,ierr)
             call MPI_ALLREDUCE(MPI_IN_PLACE,ydtemp,1,MPI_REAL8,MPI_SUM,
     $          COMM_TINKER,ierr)
             call MPI_ALLREDUCE(MPI_IN_PLACE,zdtemp,1,MPI_REAL8,MPI_SUM,
     $          COMM_TINKER,ierr)
             term = (2.0d0/3.0d0) * f * (pi/volbox)
             if (rank.eq.0) then
               em = em + term*(xd*xd+yd*yd+zd*zd)
             end if
             do ii = 1, npoleloc
                iipole = poleglob(ii)
                iglob = ipole(iipole)
                i = loc(iglob)
                dem(1,i) = dem(1,i) + 2.0d0*term*rpole(1,iipole)*xd
                dem(2,i) = dem(2,i) + 2.0d0*term*rpole(1,iipole)*yd
                dem(3,i) = dem(3,i) + 2.0d0*term*rpole(1,iipole)*zd
             end do
             xdfield = -2.0d0 * term * xd
             ydfield = -2.0d0 * term * yd
             zdfield = -2.0d0 * term * zd
             do i = 1, npoleloc
                iipole = poleglob(i)
              trq(1) = rpole(3,iipole)*zdfield - rpole(4,iipole)*ydfield
              trq(2) = rpole(4,iipole)*xdfield - rpole(2,iipole)*zdfield
              trq(3) = rpole(2,iipole)*ydfield - rpole(3,iipole)*xdfield
                call torque(iipole,trq,frcx,frcy,frcz,dem)
             end do
c
c       boundary correction to virial due to overall cell dipole
c
             xd = 0.0d0
             yd = 0.0d0
             zd = 0.0d0
             xq = 0.0d0
             yq = 0.0d0
             zq = 0.0d0
             do i = 1, npoleloc
                iipole = poleglob(i)
                iglob = ipole(iipole)
                xd = xd + rpole(2,iipole)
                yd = yd + rpole(3,iipole)
                zd = zd + rpole(4,iipole)
                xq = xq + rpole(1,iipole)*x(iglob)
                yq = yq + rpole(1,iipole)*y(iglob)
                zq = zq + rpole(1,iipole)*z(iglob)
             end do
             call MPI_ALLREDUCE(MPI_IN_PLACE,xd,1,MPI_REAL8,MPI_SUM,
     $          COMM_TINKER,ierr)
             call MPI_ALLREDUCE(MPI_IN_PLACE,yd,1,MPI_REAL8,MPI_SUM,
     $          COMM_TINKER,ierr)
             call MPI_ALLREDUCE(MPI_IN_PLACE,zd,1,MPI_REAL8,MPI_SUM,
     $          COMM_TINKER,ierr)
             call MPI_ALLREDUCE(MPI_IN_PLACE,xq,1,MPI_REAL8,MPI_SUM,
     $          COMM_TINKER,ierr)
             call MPI_ALLREDUCE(MPI_IN_PLACE,yq,1,MPI_REAL8,MPI_SUM,
     $          COMM_TINKER,ierr)
             call MPI_ALLREDUCE(MPI_IN_PLACE,zq,1,MPI_REAL8,MPI_SUM,
     $          COMM_TINKER,ierr)
             if (rank.eq.0) then
               xv = xd * xq
               yv = yd * yq
               zv = zd * zq
               vterm = term * (xd*xd + yd*yd + zd*zd + 2.0d0*(xv+yv+zv)
     &                            + xq*xq + yq*yq + zq*zq)
               vir(1,1) = vir(1,1) + 2.0d0*term*(xq*xq+xv) + vterm
               vir(2,1) = vir(2,1) + 2.0d0*term*(xq*yq+xv)
               vir(3,1) = vir(3,1) + 2.0d0*term*(xq*zq+xv)
               vir(1,2) = vir(1,2) + 2.0d0*term*(yq*xq+yv)
               vir(2,2) = vir(2,2) + 2.0d0*term*(yq*yq+yv) + vterm
               vir(3,2) = vir(3,2) + 2.0d0*term*(yq*zq+yv)
               vir(1,3) = vir(1,3) + 2.0d0*term*(zq*xq+zv)
               vir(2,3) = vir(2,3) + 2.0d0*term*(zq*yq+zv)
               vir(3,3) = vir(3,3) + 2.0d0*term*(zq*zq+zv) + vterm
             end if
          end if
        end if
      end if
      return
      end
c
c     "empole1c" calculates the multipole energy and derivatives
c     with respect to Cartesian coordinates using particle mesh Ewald
c     summation and a neighbor list
c
c
      subroutine empole1c
      use sizes
      use atmlst
      use atoms
      use boxes
      use chgpot
      use deriv
      use domdec
      use energi
      use ewald
      use math
      use mpole
      use potent
      use timestat
      use virial
      use mpi
      implicit none
      integer i,ii
      integer iipole,iglob,iloc,ierr
      real*8 e,f
      real*8 term,fterm
      real*8 cii,dii,qii
      real*8 xi,yi,zi
      real*8 xd,yd,zd
      real*8 xq,yq,zq
      real*8 xv,yv,zv,vterm
      real*8 ci,dix,diy,diz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
      real*8 xdfield,ydfield
      real*8 zdfield
      real*8 fx,fy,fz
      real*8 vxx,vyy,vzz
      real*8 vxy,vxz,vyz
      real*8 trq(3),frcx(3)
      real*8 frcy(3),frcz(3)
      real*8 time0,time1
      real*8, allocatable :: pot(:)
      real*8, allocatable :: decfx(:)
      real*8, allocatable :: decfy(:)
      real*8, allocatable :: decfz(:)
c
c
c     zero out the atomic multipole energy and derivatives
c
      em = 0.0d0
      dem = 0d0
      if (npole .eq. 0)  return
      aewald = aeewald
c
c     set the energy unit conversion factor
c
      f = electric / dielec
c
c     check the sign of multipole components at chiral sites
c
      call chkpole(.false.)
c
c     rotate the multipole components into the global frame
c
      call rotpole
c
c     compute the reciprocal space part of the Ewald summation
c
      if ((.not.(use_pmecore)).or.(use_pmecore).and.(rank.gt.ndir-1))
     $   then
        if (use_mrec) then 
          time0 = mpi_wtime()
          call emrecip1
          time1 = mpi_wtime()
          timerec = timerec + time1 - time0
        end if
      end if
c
c     compute the reciprocal space part of the Ewald summation
c
      if ((.not.(use_pmecore)).or.(use_pmecore).and.(rank.le.ndir-1))
     $   then
        if (use_mreal) then
          time0 = mpi_wtime()
          call emreal1c
          time1 = mpi_wtime()
          timereal = timereal + time1 - time0
        end if
c
        if (use_mself) then
c
c     compute the Ewald self-energy term over all the atoms
c
c
c     perform dynamic allocation of some local arrays
c
          allocate (pot(nbloc))
          allocate (decfx(nbloc))
          allocate (decfy(nbloc))
          allocate (decfz(nbloc))
          pot = 0d0
          term = 2.0d0 * aewald * aewald
          fterm = -f * aewald / sqrtpi
          do i = 1, npoleloc
             iipole = poleglob(i)
             iglob = ipole(iipole)
             iloc = loc(iglob)
             ci = rpole(1,iipole)
             dix = rpole(2,iipole)
             diy = rpole(3,iipole)
             diz = rpole(4,iipole)
             qixx = rpole(5,iipole)
             qixy = rpole(6,iipole)
             qixz = rpole(7,iipole)
             qiyy = rpole(9,iipole)
             qiyz = rpole(10,iipole)
             qizz = rpole(13,iipole)
             cii = ci*ci
             dii = dix*dix + diy*diy + diz*diz
             qii = 2.0d0*(qixy*qixy+qixz*qixz+qiyz*qiyz)
     &                + qixx*qixx + qiyy*qiyy + qizz*qizz
             e = fterm * (cii + term*(dii/3.0d0+2.0d0*term*qii/5.0d0))
             em = em + e
             pot(iloc) = 2.0d0 * fterm * ci
          end do
c
c     modify gradient and virial for charge flux self-energy
c
          if (use_chgflx) then
             call commpot(pot,1)
             call dcflux (pot,decfx,decfy,decfz)
             do i = 1, npolebloc
                iipole = poleglob(i)
                iglob = ipole(iipole)
                iloc = loc(iglob)
                xi = x(iglob)
                yi = y(iglob)
                zi = z(iglob)
                fx = decfx(iloc)
                fy = decfy(iloc)
                fz = decfz(iloc)
                dem(1,iloc) = dem(1,iloc) + fx
                dem(2,iloc) = dem(2,iloc) + fy
                dem(3,iloc) = dem(3,iloc) + fz
                vxx = xi * fx
                vxy = yi * fx
                vxz = zi * fx
                vyy = yi * fy
                vyz = zi * fy
                vzz = zi * fz
                vir(1,1) = vir(1,1) + vxx
                vir(2,1) = vir(2,1) + vxy
                vir(3,1) = vir(3,1) + vxz
                vir(1,2) = vir(1,2) + vxy
                vir(2,2) = vir(2,2) + vyy
                vir(3,2) = vir(3,2) + vyz
                vir(1,3) = vir(1,3) + vxz
                vir(2,3) = vir(2,3) + vyz
                vir(3,3) = vir(3,3) + vzz
             end do
          end if
c
c     perform deallocation of some local arrays
c
          deallocate (pot)
          deallocate (decfx)
          deallocate (decfy)
          deallocate (decfz)
c
c       compute the cell dipole boundary correction term
c
          if (boundary .eq. 'VACUUM') then
             xd = 0.0d0
             yd = 0.0d0
             zd = 0.0d0
             do i = 1, npoleloc
                iipole = poleglob(i)
                iglob = ipole(iipole)
                xd = xd + rpole(2,iipole) + rpole(1,iipole)*x(iglob)
                yd = yd + rpole(3,iipole) + rpole(1,iipole)*y(iglob)
                zd = zd + rpole(4,iipole) + rpole(1,iipole)*z(iglob)
             end do
             call MPI_ALLREDUCE(MPI_IN_PLACE,xd,1,MPI_REAL8,MPI_SUM,
     $          COMM_TINKER,ierr)
             call MPI_ALLREDUCE(MPI_IN_PLACE,yd,1,MPI_REAL8,MPI_SUM,
     $          COMM_TINKER,ierr)
             call MPI_ALLREDUCE(MPI_IN_PLACE,zd,1,MPI_REAL8,MPI_SUM,
     $          COMM_TINKER,ierr)
             term = (2.0d0/3.0d0) * f * (pi/volbox)
             if (rank.eq.0) then
               em = em + term*(xd*xd+yd*yd+zd*zd)
             end if
             do ii = 1, npoleloc
                iipole = poleglob(ii)
                iglob = ipole(iipole)
                i = loc(iglob)
                dem(1,i) = dem(1,i) + 2.0d0*term*rpole(1,iipole)*xd
                dem(2,i) = dem(2,i) + 2.0d0*term*rpole(1,iipole)*yd
                dem(3,i) = dem(3,i) + 2.0d0*term*rpole(1,iipole)*zd
             end do
             xdfield = -2.0d0 * term * xd
             ydfield = -2.0d0 * term * yd
             zdfield = -2.0d0 * term * zd
             do i = 1, npoleloc
                iipole = poleglob(i)
              trq(1) = rpole(3,iipole)*zdfield - rpole(4,iipole)*ydfield
              trq(2) = rpole(4,iipole)*xdfield - rpole(2,iipole)*zdfield
              trq(3) = rpole(2,iipole)*ydfield - rpole(3,iipole)*xdfield
                call torque(iipole,trq,frcx,frcy,frcz,dem)
             end do
c
c       boundary correction to virial due to overall cell dipole
c
             xd = 0.0d0
             yd = 0.0d0
             zd = 0.0d0
             xq = 0.0d0
             yq = 0.0d0
             zq = 0.0d0
             do i = 1, npoleloc
                iipole = poleglob(i)
                iglob = ipole(iipole)
                xd = xd + rpole(2,iipole)
                yd = yd + rpole(3,iipole)
                zd = zd + rpole(4,iipole)
                xq = xq + rpole(1,iipole)*x(iglob)
                yq = yq + rpole(1,iipole)*y(iglob)
                zq = zq + rpole(1,iipole)*z(iglob)
             end do
             call MPI_ALLREDUCE(MPI_IN_PLACE,xd,1,MPI_REAL8,MPI_SUM,
     $          COMM_TINKER,ierr)
             call MPI_ALLREDUCE(MPI_IN_PLACE,yd,1,MPI_REAL8,MPI_SUM,
     $          COMM_TINKER,ierr)
             call MPI_ALLREDUCE(MPI_IN_PLACE,zd,1,MPI_REAL8,MPI_SUM,
     $          COMM_TINKER,ierr)
             call MPI_ALLREDUCE(MPI_IN_PLACE,xq,1,MPI_REAL8,MPI_SUM,
     $          COMM_TINKER,ierr)
             call MPI_ALLREDUCE(MPI_IN_PLACE,yq,1,MPI_REAL8,MPI_SUM,
     $          COMM_TINKER,ierr)
             call MPI_ALLREDUCE(MPI_IN_PLACE,zq,1,MPI_REAL8,MPI_SUM,
     $          COMM_TINKER,ierr)
             if (rank.eq.0) then
               xv = xd * xq
               yv = yd * yq
               zv = zd * zq
               vterm = term * (xd*xd + yd*yd + zd*zd + 2.0d0*(xv+yv+zv)
     &                            + xq*xq + yq*yq + zq*zq)
               vir(1,1) = vir(1,1) + 2.0d0*term*(xq*xq+xv) + vterm
               vir(2,1) = vir(2,1) + 2.0d0*term*(xq*yq+xv)
               vir(3,1) = vir(3,1) + 2.0d0*term*(xq*zq+xv)
               vir(1,2) = vir(1,2) + 2.0d0*term*(yq*xq+yv)
               vir(2,2) = vir(2,2) + 2.0d0*term*(yq*yq+yv) + vterm
               vir(3,2) = vir(3,2) + 2.0d0*term*(yq*zq+yv)
               vir(1,3) = vir(1,3) + 2.0d0*term*(zq*xq+zv)
               vir(2,3) = vir(2,3) + 2.0d0*term*(zq*yq+zv)
               vir(3,3) = vir(3,3) + 2.0d0*term*(zq*zq+zv) + vterm
             end if
          end if
        end if
      end if
      return
      end
c
c
c     ####################################################################
c     ##                                                                ##
c     ##  subroutine emrecip1  --  PME recip multipole energy & derivs  ##
c     ##                                                                ##
c     ####################################################################
c
c
c     "emrecip1" evaluates the reciprocal space portion of the particle
c     mesh Ewald summation energy and gradient due to multipoles
c
      subroutine emrecip1
      use atmlst
      use atoms
      use bound
      use boxes
      use chgpot
      use deriv
      use domdec
      use energi
      use ewald
      use fft
      use pme
      use math
      use mpole
      use potent
      use timestat
      use virial
      use mpi
      implicit none
      integer status(MPI_STATUS_SIZE),tag,ierr
      integer i,j,k,ii,iipole,iglob,iloc,ilocrec
      integer k1,k2,k3
      integer m1,m2,m3
      integer nff
      integer nf1,nf2,nf3
      integer deriv1(10)
      integer deriv2(10)
      integer deriv3(10)
      real(t_p) e,eterm,f
      real(t_p) r1,r2,r3
      real(t_p) h1,h2,h3
      real(t_p) f1,f2,f3
      real(t_p) vxx,vyy,vzz
      real(t_p) vxy,vxz,vyz
      real(t_p) volterm,denom
      real(t_p) hsq,expterm
      real(t_p) term,pterm
      real(t_p) vterm,struc2
      real(t_p) xi,yi,zi,frcx,frcy,frcz
      real(t_p) trq(3),fix(3)
      real(t_p) fiy(3),fiz(3)
      real(t_p), allocatable :: cmp(:,:),fmp(:,:)
      real(t_p), allocatable :: qgridmpi(:,:,:,:,:)
      real(t_p), allocatable :: pot(:),potrec(:)
      real(t_p), allocatable :: decfx(:)
      real(t_p), allocatable :: decfy(:)
      real(t_p), allocatable :: decfz(:)
      integer, allocatable :: reqsend(:),reqrec(:)
      integer, allocatable :: req2send(:),req2rec(:)
      integer nprocloc,commloc,rankloc,proc
      real*8 time0,time1
      time0 = mpi_wtime()
c
      if (use_pmecore) then
        nprocloc = nrec
        commloc  = comm_rec
        rankloc  = rank_bis
      else
        nprocloc = nproc
        commloc  = COMM_TINKER
        rankloc  = rank
      end if
c
c     set Ewald coefficient
c
      aewald = aeewald

      allocate (qgridmpi(2,n1mpimax,n2mpimax,n3mpimax,nrec_recep))
      allocate (reqsend(nproc))
      allocate (reqrec(nproc))
      allocate (req2send(nproc))
      allocate (req2rec(nproc))
c
c     indices into the electrostatic field array
c
      data deriv1  / 2, 5,  8,  9, 11, 16, 18, 14, 15, 20 /
      data deriv2  / 3, 8,  6, 10, 14, 12, 19, 16, 20, 17 /
      data deriv3  / 4, 9, 10,  7, 15, 17, 13, 20, 18, 19 /
c
c
c     return if the Ewald coefficient is zero
c
      if (aewald .lt. 1.0d-6)  return
      f = electric / dielec
c
c     perform dynamic allocation of some global arrays
c
      if (allocated(cphirec)) deallocate (cphirec)
      allocate (cphirec(10,max(npolerecloc,1)))
      cphirec = 0d0
      if (allocated(fphirec)) deallocate(fphirec)
      allocate (fphirec(20,max(npolerecloc,1)))
      fphirec = 0d0
c
c     dynamic allocation of local arrays
c
      allocate (cmp(10,npolerecloc))
      allocate (fmp(10,npolerecloc))
c
c     zero out the temporary virial accumulation variables
c
      vxx = 0.0d0
      vxy = 0.0d0
      vxz = 0.0d0
      vyy = 0.0d0
      vyz = 0.0d0
      vzz = 0.0d0
c
c     zero out pme grid
c
      qgridin_2d = 0d0
c
c     copy multipole moments and coordinates to local storage
c
      do i = 1, npolerecloc
         iipole = polerecglob(i)
         iglob = ipole(iipole)
         cmp(1,i) = rpole(1,iipole)
         cmp(2,i) = rpole(2,iipole)
         cmp(3,i) = rpole(3,iipole)
         cmp(4,i) = rpole(4,iipole)
         cmp(5,i) = rpole(5,iipole)
         cmp(6,i) = rpole(9,iipole)
         cmp(7,i) = rpole(13,iipole)
         cmp(8,i) = 2.0d0 * rpole(6,iipole)
         cmp(9,i) = 2.0d0 * rpole(7,iipole)
         cmp(10,i) = 2.0d0 * rpole(10,iipole)
         call cmp_to_fmp_site(cmp(1,i),fmp(1,i))
         call bspline_fill_site(iglob,i)
      end do
c
c     assign permanent multipoles to PME grid and perform
c     the 3-D FFT forward transformation
c
      do i = 1, npolerecloc
        iipole = polerecglob(i)
        iglob = ipole(iipole)
        call grid_mpole_site(iglob,i,fmp(1,i))
      end do
      time1 = mpi_wtime()
      timegrid1 = timegrid1 + time1-time0
c
c     MPI : Begin reception
c
      time0 = mpi_wtime()
      do i = 1, nrec_recep
        tag = nprocloc*rankloc + prec_recep(i) + 1
        call MPI_IRECV(qgridmpi(1,1,1,1,i),2*n1mpimax*n2mpimax*
     $   n3mpimax,MPI_REAL8,prec_recep(i),tag,
     $   commloc,reqrec(i),ierr)
      end do
c
c     MPI : begin sending
c
      do i = 1, nrec_send
        proc = prec_send(i)
        tag = nprocloc*prec_send(i) + rankloc + 1
        call MPI_ISEND(qgridin_2d(1,1,1,1,i+1),
     $   2*n1mpimax*n2mpimax*n3mpimax,MPI_REAL8,proc,tag,
     $   commloc,reqsend(i),ierr)
      end do
c
      do i = 1, nrec_recep
        call MPI_WAIT(reqrec(i),status,ierr)
      end do
      do i = 1, nrec_send
        call MPI_WAIT(reqsend(i),status,ierr)
      end do
c
c     do the reduction 'by hand'

      do i = 1, nrec_recep
        qgridin_2d(:,:,:,:,1) = qgridin_2d(:,:,:,:,1) + 
     $    qgridmpi(:,:,:,:,i) 
      end do
      time1 = mpi_wtime()
      timerecreccomm = timerecreccomm + time1 - time0
c
c     Perform 3-D FFT forward transform
c
      time0 = mpi_wtime()
      call fft2d_frontmpi(qgridin_2d,qgridout_2d,n1mpimax,n2mpimax,
     $ n3mpimax)
      time1 = mpi_wtime()
      timeffts = timeffts + time1-time0
c
c     initialize variables required for the scalar summation
c
      time0 = mpi_wtime()
      pterm = (pi/aewald)**2
      volterm = pi * volbox
      nff = nfft1 * nfft2
      nf1 = (nfft1+1) / 2
      nf2 = (nfft2+1) / 2
      nf3 = (nfft3+1) / 2
c
      if ((istart2(rankloc+1).eq.1).and.(jstart2(rankloc+1).eq.1).and.
     $   (kstart2(rankloc+1).eq.1)) then
           qfac_2d(1,1,1) = 0.0d0
      end if
c
c     make the scalar summation over reciprocal lattice
c
      do k3 = kstart2(rankloc+1),kend2(rankloc+1)
        do k2 = jstart2(rankloc+1),jend2(rankloc+1)
          do k1 = istart2(rankloc+1),iend2(rankloc+1)
            m1 = k1 - 1
            m2 = k2 - 1
            m3 = k3 - 1
            if (k1 .gt. nf1)  m1 = m1 - nfft1
            if (k2 .gt. nf2)  m2 = m2 - nfft2
            if (k3 .gt. nf3)  m3 = m3 - nfft3
            if ((m1.eq.0).and.(m2.eq.0).and.(m3.eq.0)) goto 10
            r1 = dble(m1)
            r2 = dble(m2)
            r3 = dble(m3)
            h1 = recip(1,1)*r1 + recip(1,2)*r2 + recip(1,3)*r3
            h2 = recip(2,1)*r1 + recip(2,2)*r2 + recip(2,3)*r3
            h3 = recip(3,1)*r1 + recip(3,2)*r2 + recip(3,3)*r3
            hsq = h1*h1 + h2*h2 + h3*h3
            term = -pterm * hsq
            expterm = 0.0d0
            if (term .gt. -50.0d0) then
               denom = volterm*hsq*bsmod1(k1)*bsmod2(k2)*bsmod3(k3)
               expterm = exp(term) / denom
               if (.not. use_bounds) then
                  expterm = expterm * (1.0d0-cos(pi*xbox*sqrt(hsq)))
               else if (octahedron) then
                  if (mod(m1+m2+m3,2) .ne. 0)  expterm = 0.0d0
               end if
               struc2 = qgridout_2d(1,k1-istart2(rankloc+1)+1,
     $          k2-jstart2(rankloc+1)+1,k3-kstart2(rankloc+1)+1)**2
     &          + qgridout_2d(2,k1-istart2(rankloc+1)+1,k2-
     $          jstart2(rankloc+1)+1,k3-kstart2(rankloc+1)+1)**2
               eterm = 0.5d0 * electric * expterm * struc2
               vterm = (2.0d0/hsq) * (1.0d0-term) * eterm
               vxx = vxx + h1*h1*vterm - eterm
               vxy = vxy + h1*h2*vterm
               vxz = vxz + h1*h3*vterm
               vyy = vyy + h2*h2*vterm - eterm
               vyz = vyz + h2*h3*vterm
               vzz = vzz + h3*h3*vterm - eterm
             qfac_2d(k1-istart2(rankloc+1)+1,k2-jstart2(rankloc+1)+1,k3-
     $       kstart2(rankloc+1)+1) = expterm
           end if
 10        continue
          end do
        end do
      end do
c
c     save the virial for use in polarization computation
c
      vmxx = vxx
      vmxy = vxy
      vmxz = vxz
      vmyy = vyy
      vmyz = vyz
      vmzz = vzz
c
c     account for zeroth grid point for nonperiodic system
c
      if ((istart2(rankloc+1).eq.1).and.(jstart2(rankloc+1).eq.1)
     $   .and.(kstart2(rankloc+1).eq.1)) then
        if (.not. use_bounds) then
           expterm = 0.5d0 * pi / xbox
           struc2 = qgrid2in_2d(1,1,1,1,1)**2 +
     $       qgrid2in_2d(2,1,1,1,1)**2
           e = f * expterm * struc2
           em = em + e
        end if
      end if
c
c     complete the transformation of the PME grid
c
      do k = 1, ksize2(rankloc+1)
         do j = 1, jsize2(rankloc+1)
           do i = 1, isize2(rankloc+1)
              term = qfac_2d(i,j,k)
              qgridout_2d(1,i,j,k) = term*qgridout_2d(1,i,j,k)
              qgridout_2d(2,i,j,k) = term*qgridout_2d(2,i,j,k)
            end do
         end do
      end do
      time1 = mpi_wtime()
      timescalar = timescalar + time1-time0
c
c     perform 3-D FFT backward transform and get potential
c
      time0 = mpi_wtime()
      call fft2d_backmpi(qgridin_2d,qgridout_2d,n1mpimax,n2mpimax,
     $ n3mpimax)
      time1 = mpi_wtime()
      timeffts = timeffts + time1-time0
c
c     MPI : Begin reception
c
      time0 = mpi_wtime()
      do i = 1, nrec_send
        proc = prec_send(i)
        tag = nprocloc*rankloc + prec_send(i) + 1
        call MPI_IRECV(qgridin_2d(1,1,1,1,i+1),
     $   2*n1mpimax*n2mpimax*n3mpimax,MPI_REAL8,
     $   prec_send(i),tag,commloc,req2rec(i),ierr)
      end do
c
c     MPI : begin sending
c
      do i = 1, nrec_recep
        tag = nprocloc*prec_recep(i) + rankloc + 1
        call MPI_ISEND(qgridin_2d,
     $   2*n1mpimax*n2mpimax*n3mpimax,
     $   MPI_REAL8,prec_recep(i),tag,commloc,req2send(i),
     $   ierr)
      end do
c
      do i = 1, nrec_send
        call MPI_WAIT(req2rec(i),status,ierr)
      end do
      do i = 1, nrec_recep
        call MPI_WAIT(req2send(i),status,ierr)
      end do
      time1 = mpi_wtime()
      timerecreccomm = timerecreccomm + time1-time0

      time0 = mpi_wtime()
      do i = 1, npolerecloc
        iipole = polerecglob(i)
        iglob = ipole(iipole)
        call fphi_mpole_site(iglob,i)
      end do
      do i = 1, npolerecloc
         do j = 1, 20
            fphirec(j,i) = electric * fphirec(j,i)
         end do
        call fphi_to_cphi_site (fphirec(1,i),cphirec(1,i))
      end do
c
c     increment the permanent multipole energy and gradient
c
      e = 0.0d0
      do i = 1, npolerecloc
         f1 = 0.0d0
         f2 = 0.0d0
         f3 = 0.0d0
         do k = 1, 10
            e = e + fmp(k,i)*fphirec(k,i)
            f1 = f1 + fmp(k,i)*fphirec(deriv1(k),i)
            f2 = f2 + fmp(k,i)*fphirec(deriv2(k),i)
            f3 = f3 + fmp(k,i)*fphirec(deriv3(k),i)
         end do
         f1 = dble(nfft1) * f1
         f2 = dble(nfft2) * f2
         f3 = dble(nfft3) * f3
         h1 = recip(1,1)*f1 + recip(1,2)*f2 + recip(1,3)*f3
         h2 = recip(2,1)*f1 + recip(2,2)*f2 + recip(2,3)*f3
         h3 = recip(3,1)*f1 + recip(3,2)*f2 + recip(3,3)*f3
         iipole = polerecglob(i)
         iglob = ipole(iipole)
         ii = locrec(iglob)
         demrec(1,ii) = demrec(1,ii) + h1
         demrec(2,ii) = demrec(2,ii) + h2
         demrec(3,ii) = demrec(3,ii) + h3
      end do
      e = 0.5d0 * e
      em = em + e
c
c     distribute torques into the permanent multipole gradient
c
      do i = 1, npolerecloc
         iipole = polerecglob(i)
         trq(1) = cmp(4,i)*cphirec(3,i) - cmp(3,i)*cphirec(4,i)
     &               + 2.0d0*(cmp(7,i)-cmp(6,i))*cphirec(10,i)
     &               + cmp(9,i)*cphirec(8,i) + cmp(10,i)*cphirec(6,i)
     &               - cmp(8,i)*cphirec(9,i) - cmp(10,i)*cphirec(7,i)
         trq(2) = cmp(2,i)*cphirec(4,i) - cmp(4,i)*cphirec(2,i)
     &               + 2.0d0*(cmp(5,i)-cmp(7,i))*cphirec(9,i)
     &               + cmp(8,i)*cphirec(10,i) + cmp(9,i)*cphirec(7,i)
     &               - cmp(9,i)*cphirec(5,i) - cmp(10,i)*cphirec(8,i)
         trq(3) = cmp(3,i)*cphirec(2,i) - cmp(2,i)*cphirec(3,i)
     &               + 2.0d0*(cmp(6,i)-cmp(5,i))*cphirec(8,i)
     &               + cmp(8,i)*cphirec(5,i) + cmp(10,i)*cphirec(9,i)
     &               - cmp(8,i)*cphirec(6,i) - cmp(9,i)*cphirec(10,i)
         call torque_rec (iipole,trq,fix,fiy,fiz,demrec)
      end do
c
c     permanent multipole contribution to the internal virial
c
      do i = 1, npolerecloc
         vxx = vxx - cmp(2,i)*cphirec(2,i) - 2.0d0*cmp(5,i)*cphirec(5,i)
     &            - cmp(8,i)*cphirec(8,i) - cmp(9,i)*cphirec(9,i)
         vxy = vxy - 0.5d0*(cmp(3,i)*cphirec(2,i)+cmp(2,i)*cphirec(3,i))
     &            -(cmp(5,i)+cmp(6,i))*cphirec(8,i)
     &            -0.5d0*cmp(8,i)*(cphirec(5,i)+cphirec(6,i))
     &            -0.5d0*(cmp(9,i)*cphirec(10,i)+cmp(10,i)*cphirec(9,i))
         vxz = vxz - 0.5d0*(cmp(4,i)*cphirec(2,i)+cmp(2,i)*cphirec(4,i))
     &            -(cmp(5,i)+cmp(7,i))*cphirec(9,i)
     &            -0.5d0*cmp(9,i)*(cphirec(5,i)+cphirec(7,i))
     &            -0.5d0*(cmp(8,i)*cphirec(10,i)+cmp(10,i)*cphirec(8,i))
         vyy = vyy - cmp(3,i)*cphirec(3,i) - 2.0d0*cmp(6,i)*cphirec(6,i)
     &            - cmp(8,i)*cphirec(8,i) - cmp(10,i)*cphirec(10,i)
         vyz = vyz - 0.5d0*(cmp(4,i)*cphirec(3,i)+cmp(3,i)*cphirec(4,i))
     &            - (cmp(6,i)+cmp(7,i))*cphirec(10,i)
     &            - 0.5d0*cmp(10,i)*(cphirec(6,i)+cphirec(7,i))
     &            - 0.5d0*(cmp(8,i)*cphirec(9,i)+cmp(9,i)*cphirec(8,i))
         vzz = vzz - cmp(4,i)*cphirec(4,i) - 2.0d0*cmp(7,i)*cphirec(7,i)
     &            - cmp(9,i)*cphirec(9,i) - cmp(10,i)*cphirec(10,i)
      end do
c
c     perform dynamic allocation of some local arrays
c
      if (use_chgflx) then
         allocate (pot(nbloc))
         allocate (potrec(nlocrec))
         allocate (decfx(nbloc))
         allocate (decfy(nbloc))
         allocate (decfz(nbloc))
c
c     modify the gradient and virial for charge flux
c
         pot(:) = 0.0
         do i = 1, npolerecloc
           iipole  = polerecglob(i)
           iglob   = ipole(iipole)
           ilocrec = locrec(iglob)
           potrec(ilocrec) = cphirec(1,i)
         end do
c
c     communicate reciprocal potential to get local values
c
         call commpotrec(pot,potrec)
         call commpot(pot,1)
c
         call dcflux (pot,decfx,decfy,decfz)
         do i = 1, npolebloc
            iipole = poleglob(i)
            ii = ipole(iipole)
            iloc = loc(ii)
            xi = x(ii)
            yi = y(ii)
            zi = z(ii)
            frcx = decfx(iloc)
            frcy = decfy(iloc)
            frcz = decfz(iloc)
            dem(1,iloc) = dem(1,iloc) + frcx
            dem(2,iloc) = dem(2,iloc) + frcy
            dem(3,iloc) = dem(3,iloc) + frcz
            vxx = vxx + xi*frcx
            vxy = vxy + yi*frcx
            vxz = vxz + zi*frcx
            vyy = vyy + yi*frcy
            vyz = vyz + zi*frcy
            vzz = vzz + zi*frcz
         end do
c
c     perform deallocation of some local arrays
c
         deallocate (potrec)
         deallocate (pot)
         deallocate (decfx)
         deallocate (decfy)
         deallocate (decfz)
      end if
c
c     increment the internal virial tensor components
c
      vir(1,1) = vir(1,1) + vxx
      vir(2,1) = vir(2,1) + vxy
      vir(3,1) = vir(3,1) + vxz
      vir(1,2) = vir(1,2) + vxy
      vir(2,2) = vir(2,2) + vyy
      vir(3,2) = vir(3,2) + vyz
      vir(1,3) = vir(1,3) + vxz
      vir(2,3) = vir(2,3) + vyz
      vir(3,3) = vir(3,3) + vzz
c
      deallocate (reqsend)
      deallocate (reqrec)
      deallocate (req2send)
      deallocate (req2rec)
      time1 = mpi_wtime()
      timegrid2 = timegrid2 + time1-time0
      return
      end
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine emreal1c  --  Ewald real space derivs via list  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "emreal1c" evaluates the real space portion of the Ewald
c     summation energy and gradient due to multipole interactions
c     via a neighbor list
c
c     if longrange, calculates just the long range part
c     if shortrange, calculates just the short range part
      subroutine emreal1c
      use atmlst
      use atoms
      use bound
      use chgpen
      use chgpot
      use couple
      use cutoff
      use deriv
      use domdec
      use energi
      use group
      use iounit
      use math
      use mplpot
      use mpole
      use mutant
      use neigh
      use potent
      use shunt
      use virial
      implicit none
      integer i,j,iglob,kglob,kbis,nnelst
      integer ii,kkk,iipole,kkpole
      integer ix,iy,iz
      real(t_p) e,de,f
      real(t_p) scalek
      real(t_p) xi,yi,zi
      real(t_p) xr,yr,zr
      real(t_p) xix,yix,zix
      real(t_p) xiy,yiy,ziy
      real(t_p) xiz,yiz,ziz
      real(t_p) r,r2,rr1,rr3
      real(t_p) rr5,rr7,rr9,rr11
      real(t_p) rr1i,rr3i,rr5i,rr7i
      real(t_p) rr1k,rr3k,rr5k,rr7k
      real(t_p) rr1ik,rr3ik,rr5ik
      real(t_p) rr7ik,rr9ik,rr11ik
      real(t_p) ci,dix,diy,diz
      real(t_p) qixx,qixy,qixz
      real(t_p) qiyy,qiyz,qizz
      real(t_p) ck,dkx,dky,dkz
      real(t_p) qkxx,qkxy,qkxz
      real(t_p) qkyy,qkyz,qkzz
      real(t_p) dir,dkr,dik,qik
      real(t_p) qix,qiy,qiz,qir
      real(t_p) qkx,qky,qkz,qkr
      real(t_p) diqk,dkqi,qiqk
      real(t_p) dirx,diry,dirz
      real(t_p) dkrx,dkry,dkrz
      real(t_p) dikx,diky,dikz
      real(t_p) qirx,qiry,qirz
      real(t_p) qkrx,qkry,qkrz
      real(t_p) qikx,qiky,qikz
      real(t_p) qixk,qiyk,qizk
      real(t_p) qkxi,qkyi,qkzi
      real(t_p) qikrx,qikry,qikrz
      real(t_p) qkirx,qkiry,qkirz
      real(t_p) diqkx,diqky,diqkz
      real(t_p) dkqix,dkqiy,dkqiz
      real(t_p) diqkrx,diqkry,diqkrz
      real(t_p) dkqirx,dkqiry,dkqirz
      real(t_p) dqikx,dqiky,dqikz
      real(t_p) corei,corek
      real(t_p) vali,valk
      real(t_p) alphai,alphak
      real(t_p) term1,term2,term3
      real(t_p) term4,term5,term6
      real(t_p) term1i,term2i,term3i
      real(t_p) term1k,term2k,term3k
      real(t_p) term1ik,term2ik,term3ik
      real(t_p) term4ik,term5ik
      real(t_p) citemp,cktemp,term1temp
      real(t_p) dirtemp,dkrtemp,diktemp,term2temp
      real(t_p) qkrtemp,qirtemp,qiktemp,dkqitemp,diqktemp,qiqktemp
      real(t_p) term3temp,term4temp,term5temp
      real(t_p) poti,potk
      real(t_p) frcx,frcy,frcz
      real(t_p) vxx,vyy,vzz
      real(t_p) vxy,vxz,vyz
      real(t_p) ttmi(3),ttmk(3)
      real(t_p) fix(3),fiy(3),fiz(3)
      real(t_p) dmpi(9),dmpk(9)
      real(t_p) dmpik(11),dmpe(11)
      real(t_p) fgrp
      real(t_p) s,ds,mpoleshortcut2
      real(t_p) facts,factds
      logical testcut,shortrange,longrange,fullrange
      real(t_p), allocatable :: mscale(:)
      real(t_p), allocatable :: tem(:,:)
      real(t_p), allocatable :: pot(:)
      mdyn_rtyp, allocatable :: decfx(:)
      mdyn_rtyp, allocatable :: decfy(:)
      mdyn_rtyp, allocatable :: decfz(:)
      character*11 mode
      character*80 :: RoutineName
      external erfc
 1000 format(' Warning, system moved too much since last neighbor list'
     $   ' update, try lowering nlupdate')
c     compute the short, long, or full real space part of the summation
      shortrange = use_mpoleshortreal
      longrange  = use_mpolelong
      fullrange  = .not.(shortrange.or.longrange)
      if (shortrange) then 
         RoutineName = 'emrealshort1c'
         mode        = 'SHORTEWALD'
      else if (longrange) then
         RoutineName = 'emreallong1c'
         mode        = 'EWALD'
      else
         RoutineName = 'emreal1c'
         mode        = 'EWALD'
      endif
c
c
c     perform dynamic allocation of some local arrays
c
      allocate (mscale(n))
      allocate (tem(3,nbloc))
      allocate (pot(nbloc))
      allocate (decfx(nbloc))
      allocate (decfy(nbloc))
      allocate (decfz(nbloc))
c
c     initialize connected atom scaling and torque arrays
c
      mscale = 1.0d0
      tem = 0.0d0
      pot = 0d0
c
c     set conversion factor, cutoff and switching coefficients
c
      f = electric / dielec
      call switch (mode)
      mpoleshortcut2 = (mpoleshortcut-shortheal)**2
c
c     compute the real space portion of the Ewald summation
c
      do ii = 1, npolelocnl
         iipole = poleglobnl(ii)
         iglob  = ipole(iipole)
         i      = loc  (iglob)
         if (i.eq.0.or.i.gt.nbloc) then
           write(iout,1000)
           cycle
         endif
         xi = x(iglob)
         yi = y(iglob)
         zi = z(iglob)
         ci = rpole(1,iipole)
         dix = rpole(2,iipole)
         diy = rpole(3,iipole)
         diz = rpole(4,iipole)
         qixx = rpole(5,iipole)
         qixy = rpole(6,iipole)
         qixz = rpole(7,iipole)
         qiyy = rpole(9,iipole)
         qiyz = rpole(10,iipole)
         qizz = rpole(13,iipole)
         if (use_chgpen) then
            corei = pcore(iipole)
            vali = pval(iipole)
            alphai = palpha(iipole)
         end if
         do j = 1, n12(iglob)
            mscale(i12(j,iglob)) = m2scale
         end do
         do j = 1, n13(iglob)
            mscale(i13(j,iglob)) = m3scale
         end do
         do j = 1, n14(iglob)
            mscale(i14(j,iglob)) = m4scale
         end do
         do j = 1, n15(iglob)
            mscale(i15(j,iglob)) = m5scale
         end do
c
c     evaluate all sites within the cutoff distance
c
         if (shortrange) then
           nnelst = nshortelst(ii)
         else
           nnelst = nelst(ii)
         end if
         do kkk = 1, nnelst
            if (shortrange) then
              kkpole = shortelst(kkk,ii)
            else
              kkpole = elst(kkk,ii)
            end if
            kglob = ipole(kkpole)
            if (use_group)  call groups (fgrp,iglob,kglob,0,0,0,0)
            kbis = loc(kglob)
            if (kbis.eq.0) then
              write(iout,1000)
              cycle
            end if
            xr = x(kglob) - xi
            yr = y(kglob) - yi
            zr = z(kglob) - zi
            if (use_bounds)  call image (xr,yr,zr)
            r2 = xr*xr + yr*yr + zr*zr
            testcut = merge(r2 .le. off2.and.r2.ge.mpoleshortcut2,
     &                      r2 .le. off2,
     &                      longrange
     &                     )
            if (testcut) then
               r = sqrt(r2)
               ck = rpole(1,kkpole)
               dkx = rpole(2,kkpole)
               dky = rpole(3,kkpole)
               dkz = rpole(4,kkpole)
               qkxx = rpole(5,kkpole)
               qkxy = rpole(6,kkpole)
               qkxz = rpole(7,kkpole)
               qkyy = rpole(9,kkpole)
               qkyz = rpole(10,kkpole)
               qkzz = rpole(13,kkpole)
c
c     intermediates involving moments and separation distance
c
               dir = dix*xr + diy*yr + diz*zr
               qix = qixx*xr + qixy*yr + qixz*zr
               qiy = qixy*xr + qiyy*yr + qiyz*zr
               qiz = qixz*xr + qiyz*yr + qizz*zr
               qir = qix*xr + qiy*yr + qiz*zr
               dkr = dkx*xr + dky*yr + dkz*zr
               qkx = qkxx*xr + qkxy*yr + qkxz*zr
               qky = qkxy*xr + qkyy*yr + qkyz*zr
               qkz = qkxz*xr + qkyz*yr + qkzz*zr
               qkr = qkx*xr + qky*yr + qkz*zr
               dik = dix*dkx + diy*dky + diz*dkz
               qik = qix*qkx + qiy*qky + qiz*qkz
               diqk = dix*qkx + diy*qky + diz*qkz
               dkqi = dkx*qix + dky*qiy + dkz*qiz
               qiqk = 2.0d0*(qixy*qkxy+qixz*qkxz+qiyz*qkyz)
     &                   + qixx*qkxx + qiyy*qkyy + qizz*qkzz
c
c     additional intermediates involving moments and distance
c
               dirx = diy*zr - diz*yr
               diry = diz*xr - dix*zr
               dirz = dix*yr - diy*xr
               dkrx = dky*zr - dkz*yr
               dkry = dkz*xr - dkx*zr
               dkrz = dkx*yr - dky*xr
               dikx = diy*dkz - diz*dky
               diky = diz*dkx - dix*dkz
               dikz = dix*dky - diy*dkx
               qirx = qiz*yr - qiy*zr
               qiry = qix*zr - qiz*xr
               qirz = qiy*xr - qix*yr
               qkrx = qkz*yr - qky*zr
               qkry = qkx*zr - qkz*xr
               qkrz = qky*xr - qkx*yr
               qikx = qky*qiz - qkz*qiy
               qiky = qkz*qix - qkx*qiz
               qikz = qkx*qiy - qky*qix
               qixk = qixx*qkx + qixy*qky + qixz*qkz
               qiyk = qixy*qkx + qiyy*qky + qiyz*qkz
               qizk = qixz*qkx + qiyz*qky + qizz*qkz
               qkxi = qkxx*qix + qkxy*qiy + qkxz*qiz
               qkyi = qkxy*qix + qkyy*qiy + qkyz*qiz
               qkzi = qkxz*qix + qkyz*qiy + qkzz*qiz
               qikrx = qizk*yr - qiyk*zr
               qikry = qixk*zr - qizk*xr
               qikrz = qiyk*xr - qixk*yr
               qkirx = qkzi*yr - qkyi*zr
               qkiry = qkxi*zr - qkzi*xr
               qkirz = qkyi*xr - qkxi*yr
               diqkx = dix*qkxx + diy*qkxy + diz*qkxz
               diqky = dix*qkxy + diy*qkyy + diz*qkyz
               diqkz = dix*qkxz + diy*qkyz + diz*qkzz
               dkqix = dkx*qixx + dky*qixy + dkz*qixz
               dkqiy = dkx*qixy + dky*qiyy + dkz*qiyz
               dkqiz = dkx*qixz + dky*qiyz + dkz*qizz
               diqkrx = diqkz*yr - diqky*zr
               diqkry = diqkx*zr - diqkz*xr
               diqkrz = diqky*xr - diqkx*yr
               dkqirx = dkqiz*yr - dkqiy*zr
               dkqiry = dkqix*zr - dkqiz*xr
               dkqirz = dkqiy*xr - dkqix*yr
               dqikx = diy*qkz - diz*qky + dky*qiz - dkz*qiy
     &                 - 2.0d0*(qixy*qkxz+qiyy*qkyz+qiyz*qkzz
     &                         -qixz*qkxy-qiyz*qkyy-qizz*qkyz)
               dqiky = diz*qkx - dix*qkz + dkz*qix - dkx*qiz
     &                 - 2.0d0*(qixz*qkxx+qiyz*qkxy+qizz*qkxz
     &                         -qixx*qkxz-qixy*qkyz-qixz*qkzz)
               dqikz = dix*qky - diy*qkx + dkx*qiy - dky*qix
     &                 - 2.0d0*(qixx*qkxy+qixy*qkyy+qixz*qkyz
     &                         -qixy*qkxx-qiyy*qkxy-qiyz*qkxz)
c
c     get reciprocal distance terms for this interaction
c
               rr1 = f / r
               rr3 = rr1 / r2
               rr5 = 3.0d0 * rr3 / r2
               rr7 = 5.0d0 * rr5 / r2
               rr9 = 7.0d0 * rr7 / r2
               rr11 = 9.0d0 * rr9 / r2
c
c     calculate real space Ewald error function damping
c
               call dampewald (11,r,r2,f,dmpe)
c
c     find damped multipole intermediates and energy value
c
               if (use_chgpen) then
                  corek = pcore(kkpole)
                  valk = pval(kkpole)
                  alphak = palpha(kkpole)
                  term1 = corei*corek
                  term1i = corek*vali
                  term2i = corek*dir
                  term3i = corek*qir
                  term1k = corei*valk
                  term2k = -corei*dkr
                  term3k = corei*qkr
                  term1ik = vali*valk
                  term2ik = valk*dir - vali*dkr + dik
                  term3ik = vali*qkr + valk*qir - dir*dkr
     &                         + 2.0d0*(dkqi-diqk+qiqk)
                  term4ik = dir*qkr - dkr*qir - 4.0d0*qik
                  term5ik = qir*qkr
                  call damppole (r,11,alphai,alphak,
     &                            dmpi,dmpk,dmpik)
                  scalek = mscale(kglob)
                  if (use_group)  scalek = scalek * fgrp
                  rr1i = dmpe(1) - (1.0d0-scalek*dmpi(1))*rr1
                  rr3i = dmpe(3) - (1.0d0-scalek*dmpi(3))*rr3
                  rr5i = dmpe(5) - (1.0d0-scalek*dmpi(5))*rr5
                  rr7i = dmpe(7) - (1.0d0-scalek*dmpi(7))*rr7
                  rr1k = dmpe(1) - (1.0d0-scalek*dmpk(1))*rr1
                  rr3k = dmpe(3) - (1.0d0-scalek*dmpk(3))*rr3
                  rr5k = dmpe(5) - (1.0d0-scalek*dmpk(5))*rr5
                  rr7k = dmpe(7) - (1.0d0-scalek*dmpk(7))*rr7
                  rr1ik = dmpe(1) - (1.0d0-scalek*dmpik(1))*rr1
                  rr3ik = dmpe(3) - (1.0d0-scalek*dmpik(3))*rr3
                  rr5ik = dmpe(5) - (1.0d0-scalek*dmpik(5))*rr5
                  rr7ik = dmpe(7) - (1.0d0-scalek*dmpik(7))*rr7
                  rr9ik = dmpe(9) - (1.0d0-scalek*dmpik(9))*rr9
                  rr11ik = dmpe(11) - (1.0d0-scalek*dmpik(11))*rr11
                  rr1 = dmpe(1) - (1.0d0-scalek)*rr1
                  rr3 = dmpe(3) - (1.0d0-scalek)*rr3
                  e = term1*rr1 + term4ik*rr7ik + term5ik*rr9ik
     &                   + term1i*rr1i + term1k*rr1k + term1ik*rr1ik
     &                   + term2i*rr3i + term2k*rr3k + term2ik*rr3ik
     &                   + term3i*rr5i + term3k*rr5k + term3ik*rr5ik
c
c     find damped multipole intermediates for force and torque
c
                  de = term1*rr3 + term4ik*rr9ik + term5ik*rr11ik
     &                    + term1i*rr3i + term1k*rr3k + term1ik*rr3ik
     &                    + term2i*rr5i + term2k*rr5k + term2ik*rr5ik
     &                    + term3i*rr7i + term3k*rr7k + term3ik*rr7ik
                  term1 = -corek*rr3i - valk*rr3ik
     &                       + dkr*rr5ik - qkr*rr7ik
                  term2 = corei*rr3k + vali*rr3ik
     &                       + dir*rr5ik + qir*rr7ik
                  term3 = 2.0d0 * rr5ik
                  term4 = -2.0d0 * (corek*rr5i+valk*rr5ik
     &                                -dkr*rr7ik+qkr*rr9ik)
                  term5 = -2.0d0 * (corei*rr5k+vali*rr5ik
     &                                +dir*rr7ik+qir*rr9ik)
                  term6 = 4.0d0 * rr7ik
c
c     find standard multipole intermediates and energy value
c
               else
                  term1 = ci*ck
                  term2 = ck*dir - ci*dkr + dik
                  term3 = ci*qkr + ck*qir - dir*dkr
     &                       + 2.0d0*(dkqi-diqk+qiqk)
                  term4 = dir*qkr - dkr*qir - 4.0d0*qik
                  term5 = qir*qkr
                  scalek = 1.0d0 - mscale(kglob)
                  if (use_group)  then
                    scalek = 1.0d0 - mscale(kglob) * fgrp
                  end if
                  rr1 = dmpe(1) - scalek*rr1
                  rr3 = dmpe(3) - scalek*rr3
                  rr5 = dmpe(5) - scalek*rr5
                  rr7 = dmpe(7) - scalek*rr7
                  rr9 = dmpe(9) - scalek*rr9
                  rr11 = dmpe(11) - scalek*rr11
                  e = term1*rr1 + term2*rr3 + term3*rr5
     &                   + term4*rr7 + term5*rr9
c
                  if ((use_lambdadyn).and.(elambda.gt.0).and.
     $                (elambda.lt.1)) then
                     if (mut(iglob).and.mut(kglob)) then
                       citemp = ci/elambda
                       cktemp = ck/elambda
                       term1temp = 2d0*citemp*cktemp*elambda
                       dirtemp = dir/elambda
                       dkrtemp = dkr/elambda
                       diktemp = dik/(elambda**2)
                       term2temp = cktemp*dirtemp*2d0*elambda - 
     $                  citemp*dkrtemp*2d0*elambda + diktemp*2d0*elambda
                       qkrtemp = qkr/elambda
                       qirtemp = qir/elambda
                       qiktemp = qik/(elambda**2)
                       dkqitemp = dkqi/(elambda**2)
                       diqktemp = diqk/(elambda**2)
                       qiqktemp = qiqk/(elambda**2)
                       term3temp = citemp*qkrtemp*2d0*elambda +
     $                   cktemp*qirtemp*2d0*elambda - 
     $                   dirtemp*dkrtemp*2d0*elambda + 
     $                2.0d0*(dkqitemp*2d0*elambda-diqktemp*2d0*elambda
     $                 +qiqktemp*2d0*elambda)
                       term4temp = dirtemp*qkrtemp*2d0*elambda 
     $                   - dkrtemp*qirtemp*2d0*elambda 
     $                   - 4.0d0*qiktemp*2d0*elambda
                       term5temp = qirtemp*qkrtemp*2d0*elambda
                       delambdae = delambdae +  term1temp*rr1 
     $                 + term2temp*rr3 + term3temp*rr5  + term4temp*rr7 
     $                 + term5temp*rr9
                     else if ((mut(iglob).and..not.mut(kglob)).or.
     $                       (mut(kglob).and..not.mut(iglob))) then
                       delambdae = delambdae +  e/elambda
                     end if
                  else
                    delambdae = 0d0
                  end if
c
c     find standard multipole intermediates for force and torque
c
                  de = term1*rr3 + term2*rr5 + term3*rr7
     &                    + term4*rr9 + term5*rr11
                  term1 = -ck*rr3 + dkr*rr5 - qkr*rr7
                  term2 = ci*rr3 + dir*rr5 + qir*rr7
                  term3 = 2.0d0 * rr5
                  term4 = -2.0d0 * (ck*rr5-dkr*rr7+qkr*rr9)
                  term5 = -2.0d0 * (ci*rr5+dir*rr7+qir*rr9)
                  term6 = 4.0d0 * rr7
               end if
c
c     store the potential at each site for use in charge flux
c
               if (use_chgflx) then
                  if (use_chgpen) then
                     term1i = corek*rr1i + valk*rr1ik
                     term1k = corei*rr1k + vali*rr1ik
                     term2i = -dkr * rr3ik
                     term2k = dir * rr3ik
                     term3i = qkr * rr5ik
                     term3k = qir * rr5ik
                     poti = term1i + term2i + term3i
                     potk = term1k + term2k + term3k
                  else
                     poti = ck*rr1 - dkr*rr3 + qkr*rr5
                     potk = ci*rr1 + dir*rr3 + qir*rr5
                  end if
                  pot(i) = pot(i) + poti
                  pot(kbis) = pot(kbis) + potk
               end if 
c
c     compute the energy contributions for this interaction
c
               if(shortrange .or. longrange)
     &            call switch_respa(r,mpoleshortcut,shortheal,s,ds)

c
c     fix the s and ds factors, depending on the range
c
               if(shortrange) then
                  facts  =         s
                  factds =      - ds
               else if(longrange) then
                  facts  = 1.0d0 - s
                  factds =        ds
               else
                  facts  = 1.0d0
                  factds = 0.0d0
               endif

               em = em + facts * e
c
c     compute the force components for this interaction
c
               frcx = facts*(de*xr + term1*dix + term2*dkx
     &                   + term3*(diqkx-dkqix) + term4*qix
     &                   + term5*qkx + term6*(qixk+qkxi))
     &              + factds*xr*e/r
               frcy = facts*(de*yr + term1*diy + term2*dky
     &                   + term3*(diqky-dkqiy) + term4*qiy
     &                   + term5*qky + term6*(qiyk+qkyi))
     &                + factds*yr*e/r
               frcz = facts*(de*zr + term1*diz + term2*dkz
     &                   + term3*(diqkz-dkqiz) + term4*qiz
     &                   + term5*qkz + term6*(qizk+qkzi))
     &                + factds*zr*e/r

c              if (iglob.eq.1214.and.kglob.eq.1224) then
c14   format('empole',f4.1,/,4F16.6,/,13F10.4,/,13F10.4,/,4F16.6)
c                 print 14, mscale(kglob)
c    &                    , xr,yr,zr,r2
c    &                    , rpole(1:13,iipole),rpole(1:13,kkpole) 
c    &                    , e,frcx,frcy,frcz
c              end if
c
c     compute the torque components for this interaction
c
               if (use_chgpen)  rr3 = rr3ik
               ttmi(1) = facts*(-rr3*dikx + term1*dirx
     &                      + term3*(dqikx+dkqirx)
     &                      - term4*qirx - term6*(qikrx+qikx))
               ttmi(2) = facts*(-rr3*diky + term1*diry
     &                      + term3*(dqiky+dkqiry)
     &                      - term4*qiry - term6*(qikry+qiky))
               ttmi(3) = facts*(-rr3*dikz + term1*dirz
     &                      + term3*(dqikz+dkqirz)
     &                      - term4*qirz - term6*(qikrz+qikz))
               ttmk(1) = facts*(rr3*dikx + term2*dkrx
     &                      - term3*(dqikx+diqkrx)
     &                      - term5*qkrx - term6*(qkirx-qikx))
               ttmk(2) = facts*(rr3*diky + term2*dkry
     &                      - term3*(dqiky+diqkry)
     &                      - term5*qkry - term6*(qkiry-qiky))
               ttmk(3) = facts*(rr3*dikz + term2*dkrz
     &                      - term3*(dqikz+diqkrz)
     &                      - term5*qkrz - term6*(qkirz-qikz))
c
c     increment force-based gradient and torque on first site
c
               dem(1,i) = dem(1,i) + frcx
               dem(2,i) = dem(2,i) + frcy
               dem(3,i) = dem(3,i) + frcz
               tem(1,i) = tem(1,i) + ttmi(1)
               tem(2,i) = tem(2,i) + ttmi(2)
               tem(3,i) = tem(3,i) + ttmi(3)
c
c     increment force-based gradient and torque on second site
c
               dem(1,kbis) = dem(1,kbis) - frcx
               dem(2,kbis) = dem(2,kbis) - frcy
               dem(3,kbis) = dem(3,kbis) - frcz
               tem(1,kbis) = tem(1,kbis) + ttmk(1)
               tem(2,kbis) = tem(2,kbis) + ttmk(2)
               tem(3,kbis) = tem(3,kbis) + ttmk(3)
c
c     increment the virial due to pairwise Cartesian forces
c
               vxx = -xr * frcx
               vxy = -0.5d0 * (yr*frcx+xr*frcy)
               vxz = -0.5d0 * (zr*frcx+xr*frcz)
               vyy = -yr * frcy
               vyz = -0.5d0 * (zr*frcy+yr*frcz)
               vzz = -zr * frcz
               vir(1,1) = vir(1,1) + vxx
               vir(2,1) = vir(2,1) + vxy
               vir(3,1) = vir(3,1) + vxz
               vir(1,2) = vir(1,2) + vxy
               vir(2,2) = vir(2,2) + vyy
               vir(3,2) = vir(3,2) + vyz
               vir(1,3) = vir(1,3) + vxz
               vir(2,3) = vir(2,3) + vyz
               vir(3,3) = vir(3,3) + vzz
            end if
         end do
c
c     reset exclusion coefficients for connected atoms
c
         do j = 1, n12(iglob)
            mscale(i12(j,iglob)) = 1.0d0
         end do
         do j = 1, n13(iglob)
            mscale(i13(j,iglob)) = 1.0d0
         end do
         do j = 1, n14(iglob)
            mscale(i14(j,iglob)) = 1.0d0
         end do
         do j = 1, n15(iglob)
            mscale(i15(j,iglob)) = 1.0d0
         end do
      end do
c
c     resolve site torques then increment forces and virial
c
      do ii = 1, npolelocnl
         iipole = poleglobnl(ii)
         iglob = ipole(iipole)
         i = loc(iglob)
         call torque (iipole,tem(1,i),fix,fiy,fiz,dem)
         iz = zaxis(iipole)
         ix = xaxis(iipole)
         iy = abs(yaxis(iipole))
         if (iz .eq. 0)  iz = iglob
         if (ix .eq. 0)  ix = iglob
         if (iy .eq. 0)  iy = iglob
         xiz = x(iz) - x(iglob)
         yiz = y(iz) - y(iglob)
         ziz = z(iz) - z(iglob)
         xix = x(ix) - x(iglob)
         yix = y(ix) - y(iglob)
         zix = z(ix) - z(iglob)
         xiy = x(iy) - x(iglob)
         yiy = y(iy) - y(iglob)
         ziy = z(iy) - z(iglob)
         vxx = xix*fix(1) + xiy*fiy(1) + xiz*fiz(1)
         vxy = 0.5d0 * (yix*fix(1) + yiy*fiy(1) + yiz*fiz(1)
     &                    + xix*fix(2) + xiy*fiy(2) + xiz*fiz(2))
         vxz = 0.5d0 * (zix*fix(1) + ziy*fiy(1) + ziz*fiz(1)
     &                    + xix*fix(3) + xiy*fiy(3) + xiz*fiz(3)) 
         vyy = yix*fix(2) + yiy*fiy(2) + yiz*fiz(2)
         vyz = 0.5d0 * (zix*fix(2) + ziy*fiy(2) + ziz*fiz(2)
     &                    + yix*fix(3) + yiy*fiy(3) + yiz*fiz(3))
         vzz = zix*fix(3) + ziy*fiy(3) + ziz*fiz(3)
         vir(1,1) = vir(1,1) + vxx
         vir(2,1) = vir(2,1) + vxy
         vir(3,1) = vir(3,1) + vxz
         vir(1,2) = vir(1,2) + vxy
         vir(2,2) = vir(2,2) + vyy
         vir(3,2) = vir(3,2) + vyz
         vir(1,3) = vir(1,3) + vxz
         vir(2,3) = vir(2,3) + vyz
         vir(3,3) = vir(3,3) + vzz
      end do
c
c     modify the gradient and virial for charge flux
c
      if (use_chgflx) then
c
c      first communicate potential
c
         call commpotsum(pot,0)        
         call commpot(pot,1)        
c
         call dcflux (pot,decfx,decfy,decfz)
         do ii = 1, npolebloc
            iipole = poleglob(ii)
            iglob = ipole(iipole)
            i = loc(iglob)
            xi = x(iglob)
            yi = y(iglob)
            zi = z(iglob)
            frcx = decfx(i)
            frcy = decfy(i)
            frcz = decfz(i)
            dem(1,i) = dem(1,i) + frcx
            dem(2,i) = dem(2,i) + frcy
            dem(3,i) = dem(3,i) + frcz
            vxx = xi * frcx
            vxy = yi * frcx
            vxz = zi * frcx
            vyy = yi * frcy
            vyz = zi * frcy
            vzz = zi * frcz
            vir(1,1) = vir(1,1) + vxx
            vir(2,1) = vir(2,1) + vxy
            vir(3,1) = vir(3,1) + vxz
            vir(1,2) = vir(1,2) + vxy
            vir(2,2) = vir(2,2) + vyy
            vir(3,2) = vir(3,2) + vyz
            vir(1,3) = vir(1,3) + vxz
            vir(2,3) = vir(2,3) + vyz
            vir(3,3) = vir(3,3) + vzz
         end do
      end if
c
c     perform deallocation of some local arrays
c
      deallocate (mscale)
      deallocate (tem)
      deallocate (pot)
      deallocate (decfx)
      deallocate (decfy)
      deallocate (decfz)
      return
      end
