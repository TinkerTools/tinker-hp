#include "tinker_macro.h"
      module utilbaoab
      implicit none

      interface
      module subroutine apply_b_piston(dt,pres,stress)
      real(r_p), intent(in) :: dt,pres,stress(3,3)
      end subroutine

      module subroutine apply_o_piston(dt)
      real(r_p), intent(in) :: dt
      end subroutine

      module subroutine apply_a_piston(dt,istep,A_full)
      integer  ,intent(in):: istep
      logical  ,intent(in):: A_full
      real(r_p),intent(in):: dt
      end subroutine

      module subroutine set_langevin_thermostat_coeff(dta)
      real(r_p),intent(in):: dta
      end subroutine

      module subroutine apply_langevin_thermostat(dt)
      real(r_p),intent(in):: dt
      end subroutine
      end interface

      end module utilbaoab
