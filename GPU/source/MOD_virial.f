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
      logical:: use_virial=.false.

      real(r_p) vir(3,3),virsave(3,3)
     &         ,viramdD(3,3)
      real(r_p) g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz
     &         ,g_svxx,g_svxy,g_svxz,g_svyy,g_svyz,g_svzz

!$acc declare create(vir)
      contains

      subroutine zero_virial(viri)
      real(r_p) viri(*)
      integer i
!$acc parallel loop async present(viri)
      do i = 1,9; viri(i)=0; end do
      end subroutine

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

      subroutine info_virial1(rank)
      implicit none
      integer,intent(in):: rank
!$acc update host(g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz) async
!$acc wait

      if (rank.eq.0) then
 20   format( 80('~'))
 21   format( A,3F17.7 )
         print 20
         print 21,'virial ', g_vxx,g_vxy,g_vxz
         print 21,'       ', g_vyy,g_vyz
         print 21,'       ', g_vzz
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
