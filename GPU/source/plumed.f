c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c
c     ###################################################################
c     ##                                                               ##
c     ##  Plumed utility functions                                     ##
c     ##                                                               ##
c     ###################################################################
c
#include "tinker_macro.h"
      submodule(plumed) PlumedImp
      use atomsMirror
      use atmtyp
      use boxes
      use deriv    ,only: ftot_l, de_tot
      use tinMemory
      use tinheader,only:ti_p,re_p
      use timestat
      use inform
      use domdec
      use virial
      implicit none
      contains
#include "convert.f.inc"

      module subroutine plumed_init(dt)
      implicit none
      real(r_p),intent(in):: dt
#ifdef PLUMED
      integer i
      ! Initiate plumed
      if (deb_Path) write(*,*) "plumed_init"

      ! Allocate memory
      call prmem_requestm(pl_force,3,n)
      call prmem_requestm(pl_pos,3,n)
      call prmem_requestm(pl_mass,n)
      call prmem_request(pl_glob,n)
!$acc enter data create(pl_virial)

c
c conversion factors to go from TINKER's units to PLUMED
c
c kcal/mol -> kJ/mol
      energyUnits = 4.184_re_p
c angstrom -> nm
      lengthUnits = 0.1_re_p
c fs -> ps
      timeUnits   = 0.001_re_p

      call plumed_f_gcreate()
      call plumed_f_gcmd("setRealPrecision"//char(0),8)
      call plumed_f_gcmd("setMPIFComm"//char(0),COMM_TINKER)
      call plumed_f_gcmd("setMDEnergyUnits"//char(0),energyUnits)
      call plumed_f_gcmd("setMDLengthUnits"//char(0),lengthUnits)
      call plumed_f_gcmd("setMDTimeUnits"//char(0),timeUnits)
      call plumed_f_gcmd("setPlumedDat"//char(0),trim(pl_input)//
     $    char(0))
      call plumed_f_gcmd("setLogFile"//char(0),trim(pl_output)//char(0))
      call plumed_f_gcmd("setNatoms"//char(0),n)
      call plumed_f_gcmd("setMDEngine"//char(0),"TinkerHP");
      call plumed_f_gcmd("setTimestep"//char(0),dt);
      call plumed_f_gcmd("init"//char(0),0);

      do i = 1,nloc
         pl_glob(i) = glob(i) - 1
         pl_mass(i) = mass(glob(i))
      end do
#else
 15   format(/,
     &"                ***  FATAL ERROR *** ",/,
     &" ----------------------------------------------------",/,
     &" Plumed Feature is unavailable with this build !",/,
     &" Either Remove the keyword or rebuild Tinker-HP with ",/,
     &" PLUMED_SUPPORT=1 ",/,
     &" ----------------------------------------------------")
      if (rank.eq.0) write(0,15)
      call fatal
#endif
      end subroutine

      module subroutine eplumed(energy,derivs)
      implicit none
      real(r_p) ,intent(in)   :: energy
      real(r_p) ,intent(inout):: derivs(:,:)
#ifdef PLUMED
      integer i,j,iloc,iglob

      call timer_enter(timer_plumed)
      if (deb_Path) print*,'eplumed'

!$acc data present(pl_force,pl_pos,pl_glob,pl_mass
!$acc&    ,glob,x,y,z,mass,vir,glob,energy,derivs)

      ncount = ncount + 1
      if (nproc.eq.1) then
!$acc parallel loop async
         do i = 1,n
           pl_pos(1,i) = x(i)
           pl_pos(2,i) = y(i)
           pl_pos(3,i) = z(i)
         end do
!$acc update host(pl_pos,energy) async
      else
!$acc parallel loop async
         do iloc = 1, nloc
           iglob          = glob(iloc)
           pl_glob (iloc) = iglob - 1
           pl_pos(1,iloc) = x(iglob)
           pl_pos(2,iloc) = y(iglob)
           pl_pos(3,iloc) = z(iglob)
           pl_mass( iloc) = mass(iglob)
         enddo
!$acc update host(pl_glob,pl_pos,pl_mass,energy) async
      end if

      do i = 1,3*nloc
         pl_force(i,1) = 0.0_re_p
      end do
      pl_virial      = 0.0_re_p ! replace with actual virial
!$acc wait
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
      pl_epot        = energy
      call plumed_f_gcmd("setEnergy"//char(0),pl_epot)
c
c actual calculation
      call plumed_f_gcmd("calc"//char(0),0)

c     if (iand(ncount,511).eq.0)
c    &print '(A,f16.8,d19.10,I7)'
c    &     ,'plu',sum(pl_virial(:,:)**2),sum(pl_force(:,1:nloc)**2)
c    &     ,ncount

      if (use_virial) then
!$acc update device(pl_virial,pl_force) async
!$acc serial async
         vir(1,1) = vir(1,1) + pl_virial(1,1)
         vir(2,1) = vir(2,1) + pl_virial(1,2)
         vir(3,1) = vir(3,1) + pl_virial(1,3)
         vir(1,2) = vir(1,2) + pl_virial(2,1)
         vir(2,2) = vir(2,2) + pl_virial(2,2)
         vir(3,2) = vir(3,2) + pl_virial(2,3)
         vir(1,3) = vir(1,3) + pl_virial(3,1)
         vir(2,3) = vir(2,3) + pl_virial(3,2)
         vir(3,3) = vir(3,3) + pl_virial(3,3)
!$acc end serial
      else
!$acc update device(pl_force) async
      end if
c
c     update local derivatives
c     !pl_force is in kcal/mol/A; no conversion should be needed
c
      if (ftot_l) then
!$acc parallel loop collapse(2) async
         do i = 1,nloc; do j = 1,3
            de_tot(j,i) = de_tot(j,i) - rp2mdr(pl_force(j,i))
         end do; end do
      else
!$acc parallel loop collapse(2) async
         do i = 1,nloc; do j = 1,3
            derivs(j,i) = derivs(j,i) - pl_force(j,i)
         end do; end do
      end if
!$acc end data

      call timer_exit(timer_plumed,quiet_timers)
#endif
      end subroutine

      end submodule
