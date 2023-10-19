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
      use mutant
      use potent
      use timestat
      use vdwpot
      use virial
      use mpi
#ifdef PLUMED
      use atoms
      use atmtyp
      use boxes
      use plumed
#endif
#ifdef COLVARS
      use colvars
#endif
      implicit none
#ifdef PLUMED
      integer iglob,iloc
#endif
      real*8 energy
      real*8 derivs(3,nbloc)
      real*8 time0,time1
      logical isnan
      time0 = mpi_wtime()
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
      er = 0d0
      edsp = 0d0
      ect = 0d0
      ec = 0.0d0
      em = 0.0d0
      ep = 0.0d0
      eg = 0.0d0
      ensmd = 0.0d0
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
      deat = 0d0
      debt = 0d0
      dett = 0d0
      dev = 0d0
      der = 0d0
      dedsp = 0d0
      dect = 0d0
      dec = 0.0d0
      dem = 0d0
      dep = 0d0
      deg = 0d0
      desmd = 0d0
      dex = 0d0
#ifdef COLVARS
      delambdae = 0d0
      delambdav = 0d0
      delambda = 0d0
#endif
c
c     zero out the virial and the intermolecular energy
c
      vir = 0d0
      virsave = 0d0
      einter = 0.0d0
      time1 = mpi_wtime()
      timecleargrad = timecleargrad + time1-time0
c
c     alter partial charges and multipoles for charge flux
c
      if (use_chgflx)  call alterchg
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
      if (use_angtor)  call eangtor1
      if (use_tortor)  call etortor1
      if (use_angle)  call eangle1
      time1 = mpi_wtime()
      timebonded = timebonded + time1-time0
c
c     call the van der Waals energy and gradient routines
c
      time0 = mpi_wtime()
      if (use_vdw) then
         if (vdwtyp .eq. 'LENNARD-JONES')  call elj1
         if (vdwtyp .eq. 'BUFFERED-14-7')  call ehal1
      end if
      if (use_repuls)  call erepel1
      if (use_disp)  call edisp1
      time1 = mpi_wtime()
      timevdw = timevdw + time1-time0
c
c     call the electrostatic energy and gradient routines
c
      time0 = mpi_wtime()
      if (use_charge) call echarge1
      if (use_mpole)  call empole1
      time1 = mpi_wtime()
      timeelec = timeelec + time1-time0
      time0 = mpi_wtime()
      if (use_polar)  call epolar1
      time1 = mpi_wtime()
      timepolar = timepolar + time1-time0

      if (use_chgtrn)  call echgtrn1
c
c     call any miscellaneous energy and gradient routines
c
      if (use_geom)  call egeom1
      if (use_extra)  call extra1
c
      if (use_smd_velconst .or. use_smd_forconst) call esmd1
c
c     sum up to get the total energy and first derivatives
c
      time0 = mpi_wtime()
      esum = eit + eopd + eopb + eaa + eub + eba + ea + eb + em + ep
     $       + ec + ev + er + edsp + ect +et + ept + eat + ebt + ett 
     $       + eg + ex + eid + ensmd
      energy = esum
c
      desum = deaa + deba + dea + deb + dec + dem + deopb + deopd + deid
     $ + dep + dev + der + dedsp + dect +deub + deit + det + dept + deat
     $ + debt + dett + deg + dex + desmd
c
      debond = deaa + deba + dea + deb + deopb + deopd + deid 
     $ + deub + deit + det + dept + deat + debt + dett + deg

c
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
#ifdef COLVARS
      if (use_colvars) then
        call prepare_colvars 
c
c       only the master does colvars computations
c
        if (rank.eq.0) call compute_colvars_tinker()
        energy = esum
        call distrib_colvars(derivs)
      end if
#endif

#ifdef PLUMED
      if (lplumed) then
        pl_force = 0.0d0
        pl_virial = 0.d0 ! replace with actual virial
        pl_epot = energy

        ncount = ncount + 1
        do iloc = 1, nloc
          iglob = glob(iloc)
          pl_glob(iloc) = glob(iloc) - 1
c
c  get the unwrapped coordinates
c
          pl_pos(1,iloc) = x(iglob) + pbcwrapindex(1,iglob)*xbox
          pl_pos(2,iloc) = y(iglob) + pbcwrapindex(2,iglob)*ybox
          pl_pos(3,iloc) = z(iglob) + pbcwrapindex(3,iglob)*zbox
          pl_mass( iloc) = mass(iglob)
        enddo
c
c local number of atoms and their global indices
        call plumed_f_gcmd("setAtomsNlocal"//char(0),nloc)
        call plumed_f_gcmd("setAtomsGatindex"//char(0),pl_glob)
c
c local counter of the step (could be replaces)
        call plumed_f_gcmd("setStep"//char(0),ncount)
c
c local masses, positions and forces
        call plumed_f_gcmd("setMasses"//char(0),pl_mass)
        call plumed_f_gcmd("setPositions"//char(0),pl_pos)
        call plumed_f_gcmd("setForces"//char(0),pl_force)
c
c cell vectors
c check for non orthorhombic systems whether it's by rows or by columns)
        call plumed_f_gcmd("setBox"//char(0),lvec)
c
c virial should be both input and out
c unclear if the PLUMED virial is correct
        call plumed_f_gcmd("setVirial"//char(0),pl_virial)
c
c potential energy as collective variable
        call plumed_f_gcmd("setEnergy"//char(0),pl_epot)
c
c actual calculation
        call plumed_f_gcmd("calc"//char(0),0)
        vir(1,1) = vir(1,1) + pl_virial(1,1)
        vir(2,1) = vir(2,1) + pl_virial(1,2)
        vir(3,1) = vir(3,1) + pl_virial(1,3)
        vir(1,2) = vir(1,2) + pl_virial(2,1)
        vir(2,2) = vir(2,2) + pl_virial(2,2)
        vir(3,2) = vir(3,2) + pl_virial(2,3)
        vir(1,3) = vir(1,3) + pl_virial(3,1)
        vir(2,3) = vir(2,3) + pl_virial(3,2)
        vir(3,3) = vir(3,3) + pl_virial(3,3)
c
c unpdate local derivatives
c pl_force is in kcal/mol/A; no conversion should be needed
        derivs(:,1:nloc) = derivs(:,1:nloc) - pl_force(:,1:nloc)
      endif
#endif
      time1 = mpi_wtime()
      timecleargrad = timecleargrad + time1-time0
      return
      end
