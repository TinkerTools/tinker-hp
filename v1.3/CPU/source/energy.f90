!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     #############################################################
!     ##                                                         ##
!     ##  function energy  --  evaluates energy terms and total  ##
!     ##                                                         ##
!     #############################################################
!
!
!     "energy" calls the subroutines to calculate the potential
!     energy terms and sums up to form the total energy
!
!
!> @brief 
!> calls the subroutines to calculate the potential
!> energy terms and sums up to form the total energy
!> @param no params
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
!
   if (deb_Path) write(iout,*), 'energy '
!
!
!
!     zero out each of the potential energy components
!
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
!
!     alter partial charges and multipoles for charge flux
!
   if (use_chgflx)  call alterchg
!
!     call the local geometry energy component routines
!
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
!
!     call the van der Waals energy component routines
!
   if (use_vdw) then
      if (vdwtyp .eq. 'LENNARD-JONES')  call elj
      if (vdwtyp .eq. 'BUFFERED-14-7')  call ehal
   end if
   if (use_repuls)  call erepel
   if (use_disp)  call edisp
!
!     call the electrostatic energy component routines
!
   if (use_charge) call echarge
   if (use_mpole)  call empole0
   if (use_polar)  call epolar

   if (use_chgtrn)  call echgtrn
!
!
!     call any miscellaneous energy component routines
!
   if (use_geom)  call egeom
   if (use_extra)  call extra
!
!     sum up to give the total potential energy
!
   esum = eb + ea + eba + eub + eaa + eopb + eopd + eid + eit&
   &+ et + ept + eat + ebt + ett + ev + er + edsp + ec&
   &+ em + ect + ep + eg + ex
   energy = esum
!
!     check for an illegal value for the total energy
!
   if (isnan(esum)) then
      write (iout,10)
10    format (/,' ENERGY  --  Illegal Value for the Total',&
      &' Potential Energy')
      call fatal
   end if
   return
end
