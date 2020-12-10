c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ####################################################################
c     ##                                                                ##
c     ##  module deriv  --  Cartesian coordinate derivative components  ##
c     ##                                                                ##
c     ####################################################################
c
c
c     desum   total energy Cartesian coordinate derivatives
c     deb     bond stretch Cartesian coordinate derivatives
c     dea     angle bend Cartesian coordinate derivatives
c     deba    stretch-bend Cartesian coordinate derivatives
c     deub    Urey-Bradley Cartesian coordinate derivatives
c     deaa    angle-angle Cartesian coordinate derivatives
c     deopb   out-of-plane bend Cartesian coordinate derivatives
c     deopd   out-of-plane distance Cartesian coordinate derivatives
c     deid    improper dihedral Cartesian coordinate derivatives
c     deit    improper torsion Cartesian coordinate derivatives
c     det     torsional Cartesian coordinate derivatives
c     dept    pi-orbital torsion Cartesian coordinate derivatives
c     debt    stretch-torsion Cartesian coordinate derivatives
c     dett    torsion-torsion Cartesian coordinate derivatives
c     dev     van der Waals Cartesian coordinate derivatives
c     dec     charge-charge Cartesian coordinate derivatives
c     decrec  reciprocal charge-charge Cartesian coordinate derivatives
c     dem     multipole Cartesian coordinate derivatives
c     demrec  reciprocal multipole Cartesian coordinate derivatives
c     dep     polarization Cartesian coordinate derivatives
c     deprec  reciprocal polarization Cartesian coordinate derivatives
c     deg     geometric restraint Cartesian coordinate derivatives
c     dex     extra energy term Cartesian coordinate derivatives
c     desave  stored Cartesian coordinate derivatives
c     desmd   extra smd energy term Cartesian coordinate derivatives
c     deamdD  derivatives from aMD on dihedrals
c     deamdP  derivatives from aMD on potential energy
c     deW1aMD
c     deW2aMD
c
c     dotstgrad : flag when the main program is testgrad (communication
c      of the forces one by one)
c
c     cBond       Bonded force communication parameter 
c                 (check info_forces implementation )
c     cNBond  Non Bonded force communication parameter
c     cSNBond Short Non Bonded force communication parameter
c     cDef        Default force communication parameter
c
c
#include "tinker_precision.h"
#include "tinker_types.h"

      module deriv
      implicit none
      !DIR$ ATTRIBUTES ALIGN:64:: desum,deb,dea
      real(r_p),allocatable :: desum(:,:)
      real(r_p),allocatable :: deb(:,:),dea(:,:)
      !DIR$ ATTRIBUTES ALIGN:64:: deba,deit,dec
      real(r_p),allocatable :: deba(:,:),deit(:,:)
      mdyn_rtyp,allocatable :: dec(:,:)
      !DIR$ ATTRIBUTES ALIGN:64:: deub,deaa,deopb
      real(r_p),allocatable :: deub(:,:),deaa(:,:),deopb(:,:)
      !DIR$ ATTRIBUTES ALIGN:64:: deid,det,dept
      real(r_p),allocatable :: deid(:,:),det(:,:),dept(:,:)
      !DIR$ ATTRIBUTES ALIGN:64:: debt,dett,dev
      real(r_p),allocatable :: debt(:,:),dett(:,:)
      real(r_p),allocatable :: dev(:,:)
      !DIR$ ATTRIBUTES ALIGN:64:: dem,dep,deopd
      mdyn_rtyp,allocatable :: dem(:,:)
      mdyn_rtyp,allocatable :: dep(:,:)
      real(r_p),allocatable :: deopd(:,:)
      !DIR$ ATTRIBUTES ALIGN:64:: deg,dex
      real(r_p),allocatable :: deg(:,:),dex(:,:)
      !DIR$ ATTRIBUTES ALIGN:64:: decrec,demrec
      real(r_p),allocatable :: decrec(:,:),demrec(:,:)
      !DIR$ ATTRIBUTES ALIGN:64:: deprec
      real(r_p),allocatable :: deprec(:,:)
      !DIR$ ATTRIBUTES ALIGN:64:: debond
      mdyn_rtyp,allocatable :: debond(:,:)
      !DIR$ ATTRIBUTES ALIGN:64:: desave
      real(r_p),allocatable :: desave(:,:)
      !DIR$ ATTRIBUTES ALIGN:64:: desmd
      real(r_p),allocatable :: desmd(:,:)
      !DIR$ ATTRIBUTES ALIGN:64:: deamdD,deamdP
      real(r_p),allocatable :: deamdD(:,:),deamdP(:,:)
      !DIR$ ATTRIBUTES ALIGN:64:: deW1aMD,deW2aMD
      real(r_p),allocatable :: deW1aMD(:,:),deW2aMD(:,:)      
      logical dotstgrad
      integer cBond,cNBond,cSNBond,cDef
      parameter( cBond=1,cNBond=4,cSNBond=2,cDef=5 )

      procedure(resetForcesAmoeba),pointer::resetForces_p
      procedure(addForcesAmoeba),pointer::addForces_p

      interface
      module subroutine ConfigPotential
      end subroutine
      end interface

      interface
      module subroutine resetForcesAmoeba
      end subroutine
      module subroutine resetForcesAmoeba18
      end subroutine
      module subroutine resetForcesWaterAm
      end subroutine
      module subroutine resetForcesCharmm
      end subroutine
      module subroutine resetForcesSMD
      end subroutine
      module subroutine resetForcesAMD
      end subroutine
      module subroutine resetForcesDefault
      end subroutine
      end interface

      interface
      module subroutine addForcesAmoeba
      end subroutine
      module subroutine addForcesAmoeba18
      end subroutine
      module subroutine addForcesCharmm
      end subroutine
      module subroutine addForcesWaterAm
      end subroutine
      module subroutine addForcesSMD
      end subroutine
      module subroutine addForcesDefault
      end subroutine
      end interface

      interface
      module subroutine info_forces(rule)
      integer,intent(in)::rule
      end subroutine
      end interface

      interface
      module subroutine prtEForces(des,etot)
      real(r_p) des(:,:)
      real(r_p) etot
      end subroutine
      end interface

      contains

      subroutine deriv_void()
      implicit none
      end subroutine


      end
