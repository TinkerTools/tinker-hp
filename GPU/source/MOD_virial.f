c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ###############################################################
c     ##                                                           ##
c     ##  module virial  --  components of internal virial tensor  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     vir    total internal virial Cartesian tensor components
c     g_vx*  scalar component of virial
c     virsave stored internal virial Cartesian tensor components
c     g_svx*  scalar component of stored virial Cartesian
c     use_virial  switch on virial computing
c
c
#include "tinker_precision.h"
      module virial
      implicit none
      real(r_p) vir(3,3)
      real(r_p) g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz
      real(r_p) virsave(3,3)
      real(r_p) g_svxx,g_svxy,g_svxz,g_svyy,g_svyz,g_svzz
      real(r_p) viramdD(3,3)

      logical :: use_virial=.false.

!$acc declare create(vir)
      contains

      subroutine info_virial(rank)
      implicit none
      integer,intent(in):: rank
!$acc wait
!$acc update host(vir)

      if (rank.eq.0) then
 20   format( 80('~'))
         print 20
         print*,'virial ', vir(1,:)
         print*,'       ', vir(2,:)
         print*,'       ', vir(3,:)
      end if
      end subroutine

      subroutine create_vir_device()
!$acc enter data create(g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz)
!$acc enter data create(g_svxx,g_svxy,g_svxz,g_svyy,g_svyz,g_svzz)
!$acc enter data create(virsave,viramdD)
      end subroutine

      subroutine delete_vir_device()
!$acc exit data delete(g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz)
!$acc exit data delete(g_svxx,g_svxy,g_svxz,g_svyy,g_svyz,g_svzz)
!$acc exit data delete(virsave,viramdD)
      end subroutine

      end
