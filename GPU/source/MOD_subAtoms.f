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
      use atoms, nr=>n,xr=>x,yr=>y,zr=>z
     &         ,xor=>xold,yor=>yold,zor=>zold
      use tinheader
      use tinMemory
      use inform ,only: deb_path
      implicit none

      contains

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

      end submodule
