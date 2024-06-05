c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #############################################################
c     ##                                                         ##
c     ##  function energy  --  evaluates energy terms and total  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "energy" calls the subroutines to calculate the potential
c     energy terms and sums up to form the total energy
c
c
      function energy ()
      use sizes
      use energi
      use inform
      use iounit
      use potent
      use vdwpot
      implicit none
      real*8 energy
      logical isnan
c
      if (deb_Path) write(iout,*), 'energy '
c
c
c
c     zero out each of the potential energy components
c
      eb = 0.0d0
      ea = 0.0d0
      eba = 0.0d0
      eub = 0.0d0
      eaa = 0.0d0
      eopb = 0.0d0
      eopd = 0.0d0
      eid = 0.0d0
      eit = 0.0d0
      et = 0.0d0
      ept = 0.0d0
      eat = 0.0d0
      ebt = 0.0d0
      ett = 0.0d0
      ev = 0.0d0
      er = 0.0d0
      edsp = 0.0d0
      ect = 0.0d0
      ec = 0.0d0
      em = 0.0d0
      ep = 0.0d0
      eg = 0.0d0
      ex = 0.0d0
c
c     alter partial charges and multipoles for charge flux
c
      if (use_chgflx)  call alterchg
c
c     call the local geometry energy component routines
c
      if (use_bond)  call ebond
      if (use_angle)  call eangle
      if (use_strbnd)  call estrbnd
      if (use_urey)  call eurey
      if (use_angang)  call eangang
      if (use_opbend)  call eopbend
      if (use_opdist)  call eopdist
      if (use_improp)  call eimprop
      if (use_imptor)  call eimptor
      if (use_tors)  call etors
      if (use_pitors)  call epitors
      if (use_strtor)  call estrtor
      if (use_angtor)  call eangtor
      if (use_tortor)  call etortor
c
c     call the van der Waals energy component routines
c
      if (use_vdw) then
         if (vdwtyp .eq. 'LENNARD-JONES')  call elj
         if (vdwtyp .eq. 'BUFFERED-14-7')  call ehal
      end if
      if (use_repuls)  call erepel
      if (use_disp)  call edisp
c
c     call the electrostatic energy component routines
c
      if (use_charge) call echarge
      if (use_mpole)  call empole0
      if (use_polar)  call epolar

      if (use_chgtrn)  call echgtrn
c
c
c     call any miscellaneous energy component routines
c
      if (use_geom)  call egeom
      if (use_extra)  call extra
c
c     sum up to give the total potential energy
c
      esum = eb + ea + eba + eub + eaa + eopb + eopd + eid + eit
     &          + et + ept + eat + ebt + ett + ev + er + edsp + ec
     &          + em + ect + ep + eg + ex 
      energy = esum
c
c     check for an illegal value for the total energy
c
      if (isnan(esum)) then
         write (iout,10)
   10    format (/,' ENERGY  --  Illegal Value for the Total',
     &              ' Potential Energy')
         call fatal
      end if
      return
      end
