#ifndef CONVERT_F_INC
#define CONVERT_F_INC
#include "tinker_types.h"

#  ifdef USE_DETERMINISTIC_REDUCTION
#  define _FRAC_BITS_E Z'040000000'
#  define _FRAC_BITS_F Z'400000000'
#    ifdef TINKER_CUF
      attributes(device)
     &mdyn_rtyp function tp2mdr( input )
      use utilcu,only: __tp2ll_rz
      implicit none
      real(t_p),intent(in):: input
      mdyn_rtyp,parameter :: frac_bits=_FRAC_BITS_F
      tp2mdr = __tp2ll_rz(input*frac_bits)
      end function

      attributes(device)
     &ener_rtyp function tp2enr( input )
      use utilcu,only: __tp2ll_rz
      implicit none
      real(t_p),intent(in):: input
      ener_rtyp,parameter :: frac_bits=_FRAC_BITS_E
      tp2enr = __tp2ll_rz(input*frac_bits)
      end function

      attributes(device)
     &real(md_p) function mdr2md( input )
      use cudadevice,only: __longlong_as_double
      implicit none
      mdyn_rtyp ,intent(in):: input
      ener_rtyp ,parameter :: frac_bits=_FRAC_BITS_F
      real(md_p),parameter ::   fl_bits=real(frac_bits,md_p)
      mdr2md = __longlong_as_double(input)/fl_bits
      end function

      attributes(device)
     &real(en_p) function enr2en( input )
      use cudadevice,only: __longlong_as_double
      implicit none
      ener_rtyp ,intent(in):: input
      ener_rtyp ,parameter :: frac_bits=_FRAC_BITS_E
      real(en_p),parameter ::   fl_bits=real(frac_bits,en_p)
      enr2en =__longlong_as_double(input)/fl_bits
      end function
#    else
      mdyn_rtyp function tp2mdr( input )
!$acc routine
      implicit none
      real(t_p),intent(in):: input
      mdyn_rtyp,parameter :: frac_bits=_FRAC_BITS_F
      tp2mdr = int(input*frac_bits,8)
      !tp2mdr = nint(input*frac_bits,8) !Not Recognize intrinsic with OpenACC
      end function

      ener_rtyp function tp2enr( input )
!$acc routine
      implicit none
      real(t_p),intent(in):: input
      ener_rtyp,parameter :: frac_bits=_FRAC_BITS_E
      tp2enr = int(input*frac_bits,8)
      end function

      real(md_p) function mdr2md( input )
!$acc routine
      implicit none
      mdyn_rtyp ,intent(in):: input
      mdyn_rtyp ,parameter :: frac_bits=_FRAC_BITS_F
      real(md_p),parameter ::   fl_bits=real(frac_bits,md_p)
      mdr2md = real(input,md_p)/fl_bits
      end function

      real(en_p) function enr2en( input )
!$acc routine
      implicit none
      ener_rtyp ,intent(in):: input
      ener_rtyp ,parameter :: frac_bits=_FRAC_BITS_E
      real(en_p),parameter ::   fl_bits=real(frac_bits,en_p)
      enr2en = real(input,en_p)/fl_bits
      end function
#    endif
#  elif TINKER_MIXED_PREC
#     define tp2mdr(x) real(x,md_p)
#     define mdr2md(x) x
#     define tp2enr(x) real(x,en_p)
#     define enr2en(x) x
#  else
#     define tp2mdr(x) x
#     define mdr2md(x) x
#     define tp2enr(x) x
#     define enr2en(x) x
#  endif

#endif
