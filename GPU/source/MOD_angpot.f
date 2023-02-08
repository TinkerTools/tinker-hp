c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ###################################################################
c     ##                                                               ##
c     ##  module angpot  --  specifics of angle bend functional forms  ##
c     ##                                                               ##
c     ###################################################################
c
c
c     angunit    convert angle bending energy to kcal/mole
c     stbnunit   convert stretch-bend energy to kcal/mole
c     aaunit     convert angle-angle energy to kcal/mole
c     opbunit    convert out-of-plane bend energy to kcal/mole
c     opdunit    convert out-of-plane distance energy to kcal/mole
c     cang       cubic coefficient in angle bending potential
c     qang       quartic coefficient in angle bending potential
c     pang       quintic coefficient in angle bending potential
c     sang       sextic coefficient in angle bending potential
c     copb       cubic coefficient in out-of-plane bend potential
c     qopb       quartic coefficient in out-of-plane bend potential
c     popb       quintic coefficient in out-of-plane bend potential
c     sopb       sextic coefficient in out-of-plane bend potential
c     copd       cubic coefficient in out-of-plane distance potential
c     qopd       quartic coefficient in out-of-plane distance potential
c     popd       quintic coefficient in out-of-plane distance potential
c     sopd       sextic coefficient in out-of-plane distance potential
c     angtyp     type of angle bending function for each bond angle
c     angtypI    Integer type of angle bending function for each bond angle
c     winangtyp  window object corresponding to angtyp
c     winangtypI window object corresponding to angtypI
c     opbtyp     type of out-of-plane bend potential energy function
c     opbtypI    Integer type of out-of-plane bend potential energy function
c
c
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
      end
