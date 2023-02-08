c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ##################################################################
c     ##                                                              ##
c     ##  module ctrpot  --  charge transfer functional form details  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     ctrntyp   type of charge transfer term (SEPARATE or COMBINED)
c
c
      module ctrpot
      implicit none
      enum, bind(C)
      enumerator CHGT_SEPARATE,CHGT_COMBINED
      end enum
      character*8 ctrntyp
      integer ctrntyp_ID
      end
