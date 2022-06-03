c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #######################################################################
c     ##                                                                   ##
c     ##  module subAtoms  --  number, position and type of current atoms  ##
c     ##                                                                   ##
c     #######################################################################
c
c
#include "tinker_precision.h"
      submodule(atomsMirror) subAtoms
      use atmtyp    ,only: mass
      use atoms     , nr=>n,xr=>x,yr=>y,zr=>z
     &              , xor=>xold,yor=>yold,zor=>zold
      use domdec    ,only: rank,nproc,nloc,nbloc,glob
      use deriv     ,only: ftot_l,tdes_l,de_tot,derivx,derivy,derivz
      use inform    ,only: deb_path
      use moldyn    ,only: v
      use tinMemory
      use units     ,only: convert
      use utilgpu   ,only: mem_set,rec_stream
      use usage
      implicit none

      contains
#include "convert.f.inc"

      module subroutine atomsmirror_init
      implicit none
      integer i,j
      logical,save:: f_in=.true.

      if (nr.eq.0) then
         write(*,*) "ERROR Atoms Mirror Initialisation Failed !!"
         call fatal
      end if

      if (.not.allocated(xr).or..not.allocated(xor)) then
         write(*,*) " Positions Data not allocated "
         write(*,*) " Cannot Initialize mirror module "
         call fatal
      end if

      n = nr

#ifdef MIXED
      call prmem_requestm(x,n)
      call prmem_requestm(y,n)
      call prmem_requestm(z,n)

!$acc parallel loop default(present)
      do i = 1,n
         x(i) = 0.0_re_p
         y(i) = 0.0_re_p
         z(i) = 0.0_re_p
      end do
#else
      x => xr
      y => yr
      z => zr
#endif
      xold => xor
      yold => yor
      zold => zor

      end subroutine

      module subroutine reCast_position
      implicit none
      integer i

      if (t_p.eq.r_p) return
      if (deb_path) print*, "reCast_position"

!$acc parallel loop default(present) async
      do i = 1,n
         xr(i) = real(x(i),t_p)
         yr(i) = real(y(i),t_p)
         zr(i) = real(z(i),t_p)
      end do

      end subroutine

      module subroutine download_position(queue)
      integer,optional::queue
      if (present(queue)) then
!$acc update host(xr,yr,zr) async(queue)
      else
!$acc update host(xr,yr,zr)
      end if
      end subroutine

      module subroutine download_mirror_position(queue)
      integer,optional::queue
      if (present(queue)) then
!$acc update host(x,y,z) async(queue)
      else
!$acc update host(x,y,z)
      end if
      end subroutine

      subroutine atomsmirror_finalize
      implicit none
      end subroutine

      module subroutine integrate_vel0(a,dt)
      implicit none
      real(r_p),intent(in):: a(:,:)
      real(r_p),intent(in):: dt
      integer i,j,iglob

      if (useAll.and.nproc.eq.1) then ! Use all atoms in the box & sequential run
!$acc parallel loop collapse(2) async default(present)
         do i = 1, nloc; do j = 1, 3
            v(j,i) = v(j,i) + a(j,i)*dt
         end do; end do
      else if (useAll) then ! Use all atoms in the box
!$acc parallel loop collapse(2) async default(present)
         do i = 1, nloc; do j = 1, 3
            iglob = glob(i)
            v(j,iglob) = v(j,iglob) + a(j,iglob)*dt
         end do; end do
      else
!$acc parallel loop collapse(2) async default(present)
         do i = 1, nloc; do j = 1, 3
            iglob = glob(i)
            if (use(iglob))
     &         v(j,iglob) = v(j,iglob) + a(j,iglob)*dt
         end do; end do
      end if
      end subroutine

      module subroutine integrate_vel1(aalt,dta,a,dt)
      implicit none
      real(r_p),intent(in):: a(:,:),aalt(:,:)
      real(r_p),intent(in):: dt,dta
      integer i,j,iglob

      if (useAll.and.nproc.eq.1) then ! Use all atoms in the box & sequential run
!$acc parallel loop collapse(2) async default(present)
         do i = 1, nloc; do j = 1, 3
            v(j,i) = v(j,i) + a(j,i)*dt + aalt(j,i)*dta
         end do; end do
      else if (useAll) then ! Use all atoms in the box
!$acc parallel loop collapse(2) async default(present)
         do i = 1, nloc; do j = 1, 3
            iglob = glob(i)
            v(j,iglob) = v(j,iglob) + a(j,iglob)*dt + aalt(j,iglob)*dta
         end do; end do
      else
!$acc parallel loop collapse(2) async default(present)
         do i = 1, nloc; do j = 1, 3
            iglob = glob(i)
            if (use(iglob))
     &         v(j,iglob) =v(j,iglob) +a(j,iglob)*dt 
     &                                +aalt(j,iglob)*dta
         end do; end do
      end if
      end subroutine

      module subroutine integrate_acc_vel0(derivs,a,dt)
      implicit none
      real(r_p),intent( in):: derivs(:,:)
      real(r_p),intent(out):: a(:,:)
      real(r_p),intent( in):: dt
      integer i,j,iglob
      integer(mipk) len,nloc1
      mdyn_rtyp,parameter:: zer=0
      real(r_p) ac,acx,acy,acz

!$acc host_data use_device(glob,a,v,derivs,de_tot,mass
!$acc&         ,derivx,derivy,derivz)
      if (ftot_l) then  ! Use de_tot buffer

      if (tdes_l) then  ! de_tot transposition

      if (useAll) then ! Use all atoms in the box
!$acc parallel loop async 
!$acc&         deviceptr(v,a,mass,glob,derivx,derivy,derivz)
         do i = 1, nbloc;
            if (i.le.nloc) then
               iglob      = glob(i)
               acx        = -convert * mdr2md(derivx(i)) / mass(iglob)
               acy        = -convert * mdr2md(derivy(i)) / mass(iglob)
               acz        = -convert * mdr2md(derivz(i)) / mass(iglob)
               v(1,iglob) = v(1,iglob) + acx*dt
               v(2,iglob) = v(2,iglob) + acy*dt
               v(3,iglob) = v(3,iglob) + acz*dt
               a(1,iglob) = acx
               a(2,iglob) = acy
               a(3,iglob) = acz
               derivx(i) = 0
               derivy(i) = 0
               derivz(i) = 0
            else
               derivx(i) = 0
               derivy(i) = 0
               derivz(i) = 0
            end if
         end do;
      else
!$acc parallel loop async
!$acc&         deviceptr(v,a,mass,glob,derivx,derivy,derivz)
         do i = 1, nbloc;
            iglob = glob(i)
            if (i.le.nloc.and.use(iglob)) then
               acx        = -convert * mdr2md(derivx(i)) / mass(iglob)
               acy        = -convert * mdr2md(derivy(i)) / mass(iglob)
               acz        = -convert * mdr2md(derivz(i)) / mass(iglob)
               v(1,iglob) = v(1,iglob) + acx*dt
               v(2,iglob) = v(2,iglob) + acy*dt
               v(3,iglob) = v(3,iglob) + acz*dt
               a(1,iglob) = acx
               a(2,iglob) = acy
               a(3,iglob) = acz
               derivx(i) = 0
               derivy(i) = 0
               derivz(i) = 0
            else
               derivx(i) = 0
               derivy(i) = 0
               derivz(i) = 0
            end if
         end do;
      end if

      else

      if (useAll) then ! Use all atoms in the box
!$acc parallel loop collapse(2) async
!$acc&         deviceptr(v,a,de_tot,mass,glob)
         do i = 1, nloc; do j = 1, 3
            iglob = glob(i)
            ac         = -convert * mdr2md(de_tot(j,i)) / mass(iglob)
            v(j,iglob) = v(j,iglob) + ac*dt
            a(j,iglob) = ac
            de_tot(j,i)= 0
         end do; end do
      else
!$acc parallel loop collapse(2) async
!$acc&         deviceptr(v,a,de_tot,mass,glob)
         do i = 1, nloc; do j = 1, 3
            iglob = glob(i)
            if (use(iglob)) then
               ac         = -convert * mdr2md(de_tot(j,i)) / mass(iglob)
               v(j,iglob) = v(j,iglob) + ac*dt
               a(j,iglob) = ac
               de_tot(j,i)= 0
            end if
         end do; end do
      end if
      len = 3*(nbloc-nloc); nloc1 = 3*nloc

      end if  ! de_tot transposition

      else

      if (useAll) then ! Use all atoms in the box
!$acc parallel loop collapse(2) async
!$acc&         deviceptr(v,a,derivs,mass,glob)
         do i = 1, nloc; do j = 1, 3
            iglob = glob(i)
            ac         = -convert * derivs(j,i) / mass(iglob)
            v(j,iglob) = v(j,iglob) + ac*dt
            a(j,iglob) = ac
         end do; end do
      else
!$acc parallel loop collapse(2) async
!$acc&         deviceptr(v,a,derivs,mass,glob)
         do i = 1, nloc; do j = 1, 3
            iglob = glob(i)
            if (use(iglob)) then
               ac         = -convert * derivs(j,i) / mass(iglob)
               v(j,iglob) = v(j,iglob) + ac*dt
               a(j,iglob) = ac
            end if
         end do; end do
      end if

      end if  ! Use de_tot buffer
!$acc end host_data

      if (ftot_l.and..not.tdes_l.and.len>0)
     &   call mem_set(de_tot,zer,len,rec_stream,nloc1)
      end subroutine

      module subroutine integrate_acc_vel1(aalt,dta,derivs,a,dt)
      implicit none
      real(r_p),intent( in):: derivs(:,:),aalt(:,:)
      real(r_p),intent(out):: a(:,:)
      real(r_p),intent( in):: dt,dta
      integer i,j,iglob
      integer(mipk) slen,nloc1
      mdyn_rtyp,parameter:: zer=0
      real(r_p) ac,acx,acy,acz

!$acc host_data use_device(glob,a,aalt,v,derivs,de_tot,mass
!$acc&         ,derivx,derivy,derivz)
      if (ftot_l) then

      if (tdes_l) then  ! de_tot transposition

      if (useAll) then ! Use all atoms in the box
!$acc parallel loop async
!$acc&         deviceptr(v,a,de_tot,mass,glob
!$acc&                  ,derivx,derivy,derivz)
         do i = 1, nbloc;
            if (i.le.nloc) then
               iglob      = glob(i)
               acx        = -convert * mdr2md(derivx(i)) / mass(iglob)
               acy        = -convert * mdr2md(derivy(i)) / mass(iglob)
               acz        = -convert * mdr2md(derivz(i)) / mass(iglob)
               v(1,iglob) = v(1,iglob) + acx*dt + aalt(1,iglob)*dta
               v(2,iglob) = v(2,iglob) + acy*dt + aalt(2,iglob)*dta
               v(3,iglob) = v(3,iglob) + acz*dt + aalt(3,iglob)*dta
               a(1,iglob) = acx
               a(2,iglob) = acy
               a(3,iglob) = acz
               derivx(i) = 0
               derivy(i) = 0
               derivz(i) = 0
            else
               derivx(i) = 0
               derivy(i) = 0
               derivz(i) = 0
            end if
         end do;
      else
!$acc parallel loop async
!$acc&         deviceptr(v,a,de_tot,mass,glob
!$acc&                  ,derivx,derivy,derivz)
         do i = 1, nbloc;
            iglob = glob(i)
            if (i.le.nloc.and.use(iglob)) then
               acx        = -convert * mdr2md(derivx(i)) / mass(iglob)
               acy        = -convert * mdr2md(derivy(i)) / mass(iglob)
               acz        = -convert * mdr2md(derivz(i)) / mass(iglob)
               v(1,iglob) = v(1,iglob) + acx*dt + aalt(1,iglob)*dta
               v(2,iglob) = v(2,iglob) + acy*dt + aalt(2,iglob)*dta
               v(3,iglob) = v(3,iglob) + acz*dt + aalt(3,iglob)*dta
               a(1,iglob) = acx
               a(2,iglob) = acy
               a(3,iglob) = acz
               derivx(i) = 0
               derivy(i) = 0
               derivz(i) = 0
            else
               derivx(i) = 0
               derivy(i) = 0
               derivz(i) = 0
            end if
         end do;
      end if

      else

      if (useAll) then ! Use all atoms in the box
!$acc parallel loop collapse(2) async
!$acc&         deviceptr(v,a,aalt,de_tot,mass,glob)
         do i = 1, nloc; do j = 1, 3
            iglob = glob(i)
            ac         = -convert * mdr2md(de_tot(j,i)) / mass(iglob)
            v(j,iglob) = v(j,iglob) + ac*dt +aalt(j,iglob)*dta
            a(j,iglob) = ac
            de_tot(j,i)= 0
         end do; end do
      else
!$acc parallel loop collapse(2) async
!$acc&         deviceptr(v,a,aalt,de_tot,mass,glob)
         do i = 1, nloc; do j = 1, 3
            iglob = glob(i)
            if (use(iglob)) then
               ac         = -convert * mdr2md(de_tot(j,i)) / mass(iglob)
               v(j,iglob) = v(j,iglob) + ac*dt +aalt(j,iglob)*dta
               a(j,iglob) = ac
               de_tot(j,i)= 0
            end if
         end do; end do
      end if
      slen = 3*(nbloc-nloc); nloc1 = 3*nloc 

      end if  ! de_tot transposition

      else

      if (useAll) then ! Use all atoms in the box
!$acc parallel loop collapse(2) async 
!$acc&         deviceptr(v,a,aalt,derivs,mass,glob)
         do i = 1, nloc; do j = 1, 3
            iglob = glob(i)
            ac         = -convert * derivs(j,i) / mass(iglob)
            v(j,iglob) = v(j,iglob) + ac*dt +aalt(j,iglob)*dta
            a(j,iglob) = ac
         end do; end do
      else
!$acc parallel loop collapse(2) async
!$acc&         deviceptr(v,a,aalt,derivs,mass,glob)
         do i = 1, nloc; do j = 1, 3
            iglob = glob(i)
            if (use(iglob)) then
               ac         = -convert * derivs(j,i) / mass(iglob)
               v(j,iglob) = v(j,iglob) + ac*dt +aalt(j,iglob)*dta
               a(j,iglob) = ac
            end if
         end do; end do
      end if

      end if
!$acc end host_data

      if (ftot_l.and..not.tdes_l.and.slen>0)
     &   call mem_set(de_tot,zer,slen,rec_stream,nloc1)
      end subroutine

      module subroutine integrate_pos(dt)
      implicit none
      real(r_p),intent(in):: dt
      integer i,j,iglob

!$acc host_data use_device(x,y,z,xr,yr,zr,v)
      if (useAll.and.nproc.eq.1) then
!$acc parallel loop async
!$acc&         deviceptr(x,y,z,xr,yr,zr,v)
         do iglob = 1, n
            x(iglob) = x(iglob) + v(1,iglob)*dt
            y(iglob) = y(iglob) + v(2,iglob)*dt
            z(iglob) = z(iglob) + v(3,iglob)*dt
            if (t_p.ne.r_p) then
               xr(iglob) = real( x(iglob),t_p )
               yr(iglob) = real( y(iglob),t_p )
               zr(iglob) = real( z(iglob),t_p )
            end if
         end do
      else
!$acc parallel loop async default(present)
!$acc&         deviceptr(x,y,z,xr,yr,zr,v)
         do i = 1, nloc
            iglob = glob(i)
            if (use(iglob)) then
               x(iglob) = x(iglob) + v(1,iglob)*dt
               y(iglob) = y(iglob) + v(2,iglob)*dt
               z(iglob) = z(iglob) + v(3,iglob)*dt
               if (t_p.ne.r_p.and.nproc.eq.1) then
                  xr(iglob) = real( x(iglob),t_p )
                  yr(iglob) = real( y(iglob),t_p )
                  zr(iglob) = real( z(iglob),t_p )
               end if
            end if
         end do
      end if
!$acc end host_data

      end subroutine

      end submodule
