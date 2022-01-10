c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine gradient  --  find energy & gradient components  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "gradient" calls subroutines to calculate the potential energy
c     and first derivatives with respect to Cartesian coordinates
c
c
      subroutine gradient (energy,derivs)
      use sizes
      use deriv
      use domdec
      use energi
      use inter
      use iounit
      use potent
      use timestat
      use vdwpot
      use virial
      use mpi
      implicit none
      integer i
      real*8 energy
      real*8 derivs(3,nbloc)
      real*8 time0,time1,time2
      logical isnan
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
      ebt = 0.0d0
      ett = 0.0d0
      ev = 0.0d0
      ec = 0.0d0
      em = 0.0d0
      ep = 0.0d0
      eg = 0.0d0
      erep = 0.0d0
      exdisp = 0.0d0
      ect = 0.0d0
      ex = 0.0d0
c
c     zero out each of the first derivative components
c
      deb = 0d0
      dea = 0d0
      deba = 0d0
      deub = 0d0
      deaa = 0d0
      deopb = 0d0
      deopd = 0d0
      deid = 0d0
      deit = 0d0
      det = 0d0
      dept = 0d0
      debt = 0d0
      dett = 0d0
      dev = 0d0
      dec = 0.0d0
      decrec = 0.0d0
      dem = 0d0
      demrec = 0d0
      dep = 0d0
      deprec= 0d0
      deg = 0d0
      derep = 0d0
      dexdisp = 0d0
      dect = 0d0
      dex = 0d0
c
c     zero out the virial and the intermolecular energy
c
      vir = 0d0
      einter = 0.0d0
c
c     call the local geometry energy and gradient routines
c
      time0 = mpi_wtime()
      if (use_bond)  call ebond1
      if (use_strbnd)  call estrbnd1
      if (use_urey)  call eurey1
      if (use_angang)  call eangang1
      if (use_opbend)  call eopbend1
      if (use_opdist)  call eopdist1
      if (use_improp)  call eimprop1
      if (use_imptor)  call eimptor1
      if (use_tors)  call etors1
      if (use_pitors)  call epitors1
      if (use_strtor)  call estrtor1
      if (use_tortor)  call etortor1
      if (use_angle)  call eangle1
      time1 = mpi_wtime()
      timebonded = timebonded + time1-time0
c
c     call the van der Waals energy and gradient routines
c
      if (use_vdw) then
         if (vdwtyp .eq. 'LENNARD-JONES')  call elj1
         if (vdwtyp .eq. 'BUFFERED-14-7')  call ehal1
      end if
c
c     call the electrostatic energy and gradient routines
c
      if (use_charge) call echarge1
      if (use_mpole)  call empole1
      if (use_polar)  call epolar1
c
      if (use_repulsion) call erepulsion1
      if (use_dispersion) call edispersion1
      if (use_ctransfer) call ectransfer1
c
c     call any miscellaneous energy and gradient routines
c
      if (use_geom)  call egeom1
      if (use_extra)  call extra1
      time2 = mpi_wtime()
      timenonbonded = timenonbonded + time2-time1
c
c     sum up to get the total energy and first derivatives
c
      esum = eit + eopd + eopb + eaa + eub + eba + ea + eb + em + ep
     $       + ec + ev + et + ept + ebt + ett + eg + ex + eid 
     $       + erep + exdisp + ect 
      energy = esum
c
      desum = deaa + deba + dea + deb + dec + dem + deopb + deopd + deid
     $ + dep + dev + deub + deit + det + dept + debt + dett + deg + dex
     $ + derep + dexdisp + dect
c
      debond = deaa + deba + dea + deb + deopb + deopd + deid 
     $ + deub + deit + det + dept + debt + dett + deg

c
c      derivs(:,1:nloc) = derivs(:,1:nloc) + desum(:,1:nloc)
      derivs(:,1:nbloc) = derivs(:,1:nbloc) + desum(:,1:nbloc)
c
c     check for an illegal value for the total energy
c
      if (isnan(esum)) then
         write (iout,10)
   10    format (/,' GRADIENT  --  Illegal Value for the Total',
     &              ' Potential Energy')
         call fatal
      end if
      return
      end
