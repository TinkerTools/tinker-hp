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
c     derivs  total energy Cartesian coordinate derivatives
c     deb     bond stretch Cartesian coordinate derivatives
c     dea     angle bend Cartesian coordinate derivatives
c     dmlpot  Machine learning potential Cartesian coordinate derivatives
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
c     deat    angle-torsion Cartesian coordinate derivatives
c     dett    torsion-torsion Cartesian coordinate derivatives
c     dev     van der Waals Cartesian coordinate derivatives
c     der     repulsion Cartesian coordinate derivatives
c     dedsp   dispersion Cartesian coordinate derivatives
c     dedsprec   reciprocal dispersion Cartesian coordinate derivatives
c     dect    charge transfer Cartesian coordinate derivatives
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
c     de_mpi  forces communication buffer
c     de_ws   work space buffer
c     deW1aMD
c     deW2aMD
c
c     Lambda-dynamics derivatives
c
c     delambda           hamiltonian derivative with respect to lambda (to be sent to colvar)
c     delambdae          hamiltonian derivative with respect to elambda
c     delambdav          hamiltonian derivative with respect to vlambda
c     dlambdaelambda     derivative of elambda with respect to lambda
c     dlambdavlambda     derivative of vlambda with respect to lambda     
c
c     Orthogonal Space Random Walk - note x stands for Cartesian coordinates
c     dxdelambda         hamiltonian double derivative with respect to x and lambda (to be sent to colvar)
c     dxdelambdae        hamiltonian double derivative with respect to x and elambda (electrostatic interactions)
c     dxdelambdav        hamiltonian double derivative with respect to x and vlambda (vdw interactions)
c     d2edlambda2         hamiltonian double derivative with respect to lambda (to be sent to colvar)
c     d2edlambdae2        hamiltonian double derivative with respect to elambda (electrostatic interactions)
c     d2edlambdav2        hamiltonian double derivative with respect to vlambda (vdw interactions)
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
c     ftot_l      use /de_tot/ buffer for partial summation
c     fdebs_l     use /de_tot/ for all bonded buffers (switch)
c     tdes_l      transpose force buffers in a (n,3) shape (switch) 
c     dr_stride   inner stride between one set of coordinates of an atom (x -- y -- z)
c     dr_stride3  stride between two consecutive buffers ( deb -- dea )
c     dr_nb0      number of buffers used in dr_buf0
c     dr_nbb      number of buffers used for bonded forces
c     dr_nbnb     number of buffers used for non bonded forces
c     dr_nbnbr    number of buffers used for non bonded reciprocal forces
c     dr_obb      offset buffer for bonded forces
c     dr_obnb     offset buffer for non bonded forces
c     dr_obnbr    offset buffer for non bonded reciprocal forces
c
c
#include "tinker_macro.h"

      module deriv
      implicit none
      logical fdebs_l, tdes_l, ftot_l
      integer dr_stride,dr_stride3,dr_strider,dr_strider3
      integer dr_nb0,dr_nb1,dr_nbb,dr_nbnb,dr_nbnbr
      integer(ipk_) dr_obb,dr_obnb,dr_obnbr,dr_obws

      real(r_p),allocatable,target :: de_buff0(:)
      real(r_p),allocatable,target :: de_buffr(:)
      mdyn_rtyp,allocatable,target :: de_buff1(:)

      mdyn_rtyp,pointer :: d_x(:),d_y(:),d_z(:)
      mdyn_rtyp,pointer :: de1x(:),de1y(:),de1z(:),de1d(:)

      mdyn_rtyp,pointer :: de_tot(:,:),de_tot1(:)
      mdyn_rtyp,pointer :: de_mpi(:,:)
      real(r_p),pointer :: de_ws0(:,:)
      mdyn_rtyp,pointer :: de_ws1(:,:)
      mdyn_rtyp,pointer :: de_ws2(:,:)

      !DIR$ ATTRIBUTES ALIGN:64:: deb,dea
      real(r_p),pointer :: deb(:,:),dea(:,:),deub(:,:)
      !DIR$ ATTRIBUTES ALIGN:64:: deba,deit,dec
      real(r_p),pointer :: deba(:,:),deopb(:,:),deopd(:,:),deaa(:,:)
      !DIR$ ATTRIBUTES ALIGN:64:: deub,deaa
      real(r_p),pointer :: deit(:,:),deid(:,:)
      !DIR$ ATTRIBUTES ALIGN:64:: debt,dett,deat,deopd,deopb
      real(r_p),pointer :: det(:,:),dept(:,:),dett(:,:)
     &         ,debt(:,:),deat(:,:)
      !DIR$ ATTRIBUTES ALIGN:64:: deg,dex
      real(r_p),pointer :: deg(:,:),dex(:,:)
      real(r_p),pointer :: dmlpot(:,:)

      mdyn_rtyp,pointer,dimension(:) :: derivx,derivy,derivz
     &         ,debx,deby,debz,deax,deay,deaz,deubx,deuby,deubz

      !DIR$ ATTRIBUTES ALIGN:64:: desmd,dev,dem,dep,deopd
      mdyn_rtyp,pointer :: desmd(:,:),dev(:,:),der(:,:),dedsp(:,:)
     &         ,dem(:,:),dep(:,:),dec(:,:),dect(:,:)
      !DIR$ ATTRIBUTES ALIGN:64:: decrec,demrec,deprec
      real(r_p),pointer :: demrec(:,:),deprec(:,:),decrec(:,:)
     &         ,dedsprec(:,:)

      mdyn_rtyp,pointer,dimension(:) :: devx,devy,devz,decx,decy,decz
     &         ,demx,demy,demz,depx,depy,depz
     &         ,decrecx,decrecy,decrecz,demrecx,demrecy,demrecz
     &         ,deprecx,deprecy,deprecz

      !DIR$ ATTRIBUTES ALIGN:64:: deamdD,deamdP,deW1aMD,deW2aMD
      real(r_p),pointer :: deamdD(:,:),deamdP(:,:), deW1aMD(:,:)
     &         ,deW2aMD(:,:)
      !DIR$ ATTRIBUTES ALIGN:64:: desave
      mdyn_rtyp,pointer ::desave(:,:)

      real(r_p) delambda,delambdae,delambdav
      real(r_p) d2edlambda2,d2edlambdae2,d2edlambdav2
      real(r_p) dlambdaelambda, dlambdavlambda
      real(r_p), allocatable :: dxdelambda(:,:)
      real(r_p), allocatable :: dxdelambdae(:,:), dxdelambdav(:,:)

      logical dotstgrad
      integer cBond,cNBond,cSNBond,cDef
      enum,bind(C)
      enumerator idBond,idSNBond,idNBond
      end enum
      parameter( cBond=1,cSNBond=2,cNBond=4,cDef=5 )

      interface
      module subroutine ConfigPotential
      end subroutine
      end interface

      interface
      module subroutine mem_alloc_deriv(opt)
      integer,optional:: opt
      end subroutine
      module subroutine mem_free_deriv
      end subroutine
      end interface

      interface
      module subroutine zero_forces
      end subroutine
      module subroutine zero_forces_host
      end subroutine
      module subroutine zero_forces_rec
      end subroutine
      module subroutine check_nzero
      end subroutine
      module subroutine get_ftot(derivs,nbloc_)
      integer   nbloc_
      real(r_p) derivs(3,nbloc_)
      end subroutine
      module subroutine add_forces(deadd)
      mdyn_rtyp deadd(*)
      end subroutine
      module subroutine add_forces1(deadd)
      real(r_p) deadd(*)
      end subroutine
      module subroutine add_forces_rec(deadd)
      mdyn_rtyp deadd(dr_strider3)
      end subroutine
      module subroutine add_forces_rec1(deadd)
      real(r_p) deadd(dr_strider3)
      end subroutine
      module subroutine sum_forces_rec1(desum)
      real(r_p) desum(dr_stride3)
      end subroutine
      module subroutine add_forces_rec_1d(deadd)
      mdyn_rtyp deadd(dr_stride3)
      end subroutine
      module subroutine add_forces_rec_1d1(deadd,derec)
      mdyn_rtyp deadd(dr_stride3)
      real(r_p) derec(dr_strider3)
      end subroutine
      module subroutine add_forces_rec1_1d(deadd)
      real(r_p) deadd(dr_stride3)
      end subroutine
      module subroutine add_forces_rec1_1d1(deadd,derec)
      real(r_p) deadd(dr_stride3), derec(dr_strider3)
      end subroutine
      module subroutine remove_desave(derivs)
      real(r_p) derivs(3,dr_stride)
      end subroutine
      end interface

      interface
      module subroutine resetForcesAMD
      end subroutine
      end interface

      ! Force Communications routines
      interface
      module subroutine comm_forces_dd(derivs,opt)
      mdyn_rtyp derivs(3,*)
      integer,optional:: opt
      end subroutine
      module subroutine comm_forces_dd1(derivs,opt)
      real(r_p) derivs(3,*)
      integer,optional:: opt
      end subroutine
      module subroutine comm_forces_rec(de_rec,opt)
      real(r_p) de_rec(3,dr_strider)
      integer,optional:: opt
      end subroutine
      module subroutine comm_forces_recdir(derivs,de_rec,opt)
      implicit none
      mdyn_rtyp,intent(inout):: derivs(dr_stride3)
      real(r_p),intent(in)::    de_rec(3,dr_strider)
      integer  ,optional:: opt
      end subroutine
      module subroutine comm_forces_recdir1(derivs,de_rec,opt)
      implicit none
      real(r_p),intent(inout):: derivs(3,dr_stride)
      real(r_p),intent(in)::    de_rec(3,dr_strider)
      integer  ,optional:: opt
      end subroutine
      module subroutine comm_forces(derivs,opt)
      real(r_p) derivs(3,*)
      integer  ,optional:: opt
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
