!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ###################################################################
!     ##                                                               ##
!     ##  module angpot  --  specifics of angle bend functional forms  ##
!     ##                                                               ##
!     ###################################################################
!
!
!     angunit    convert angle bending energy to kcal/mole
!     stbnunit   convert stretch-bend energy to kcal/mole
!     aaunit     convert angle-angle energy to kcal/mole
!     opbunit    convert out-of-plane bend energy to kcal/mole
!     opdunit    convert out-of-plane distance energy to kcal/mole
!     cang       cubic coefficient in angle bending potential
!     qang       quartic coefficient in angle bending potential
!     pang       quintic coefficient in angle bending potential
!     sang       sextic coefficient in angle bending potential
!     copb       cubic coefficient in out-of-plane bend potential
!     qopb       quartic coefficient in out-of-plane bend potential
!     popb       quintic coefficient in out-of-plane bend potential
!     sopb       sextic coefficient in out-of-plane bend potential
!     copd       cubic coefficient in out-of-plane distance potential
!     qopd       quartic coefficient in out-of-plane distance potential
!     popd       quintic coefficient in out-of-plane distance potential
!     sopd       sextic coefficient in out-of-plane distance potential
!     angtyp     type of angle bending function for each bond angle
!     angtypI    Integer type of angle bending function for each bond angle
!     winangtyp  window object corresponding to angtyp
!     winangtypI window object corresponding to angtypI
!     opbtyp     type of out-of-plane bend potential energy function
!     opbtypI    Integer type of out-of-plane bend potential energy function
!
!
#include "tinker_macro.h"
module angpot
   implicit none
   enum, bind(C)
      enumerator OPB_W_D_C
      enumerator OPB_ALLINGER
   end enum
   enum,bind(C)
      enumerator ANG_HARMONIC, ANG_IN_PLANE
      enumerator ANG_FOURIER,  ANG_LINEAR
      enumerator ANG_PS
   end enum
   real(t_p) angunit,stbnunit,aaunit
   real(t_p) opbunit,opdunit
   real(t_p) cang,qang,pang,sang
   real(t_p) copb,qopb,popb,sopb
   real(t_p) copd,qopd,popd,sopd
   integer     opbtypInt
   character*8 opbtyp
   character*8, pointer :: angtyp(:)
   integer    , pointer :: angtypI(:)
   integer :: winangtyp, winangtypI
   real(t_p) :: c5z_ps(245)
   integer :: idx_ps(245,3)
   real(t_p), allocatable :: fmat_ps(:,:,:), dfmat_ps(:,:,:)

!$acc declare create(c5z_ps,idx_ps)
end
