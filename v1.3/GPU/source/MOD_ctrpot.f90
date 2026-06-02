!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ##################################################################
!     ##                                                              ##
!     ##  module ctrpot  --  charge transfer functional form details  ##
!     ##                                                              ##
!     ##################################################################
!
!
!     ctrntyp   type of charge transfer term (SEPARATE or COMBINED)
!
!
module ctrpot
   implicit none
   enum, bind(C)
      enumerator CHGT_SEPARATE,CHGT_COMBINED
   end enum
   character*8 ctrntyp
   integer ctrntyp_ID
end
