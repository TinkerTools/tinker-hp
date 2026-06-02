!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ####################################################################
!     ##                                                                ##
!     ##  module deriv  --  Cartesian coordinate derivative components  ##
!     ##                                                                ##
!     ####################################################################
!
!
!     derivs  total energy Cartesian coordinate derivatives
!     deb     bond stretch Cartesian coordinate derivatives
!     dea     angle bend Cartesian coordinate derivatives
!     dmlpot  Machine learning potential Cartesian coordinate derivatives
!     deba    stretch-bend Cartesian coordinate derivatives
!     deub    Urey-Bradley Cartesian coordinate derivatives
!     deaa    angle-angle Cartesian coordinate derivatives
!     deopb   out-of-plane bend Cartesian coordinate derivatives
!     deopd   out-of-plane distance Cartesian coordinate derivatives
!     deid    improper dihedral Cartesian coordinate derivatives
!     deit    improper torsion Cartesian coordinate derivatives
!     det     torsional Cartesian coordinate derivatives
!     dept    pi-orbital torsion Cartesian coordinate derivatives
!     debt    stretch-torsion Cartesian coordinate derivatives
!     deat    angle-torsion Cartesian coordinate derivatives
!     dett    torsion-torsion Cartesian coordinate derivatives
!     dev     van der Waals Cartesian coordinate derivatives
!     der     repulsion Cartesian coordinate derivatives
!     dedsp   dispersion Cartesian coordinate derivatives
!     dedsprec   reciprocal dispersion Cartesian coordinate derivatives
!     dect    charge transfer Cartesian coordinate derivatives
!     dec     charge-charge Cartesian coordinate derivatives
!     decrec  reciprocal charge-charge Cartesian coordinate derivatives
!     dem     multipole Cartesian coordinate derivatives
!     demrec  reciprocal multipole Cartesian coordinate derivatives
!     dep     polarization Cartesian coordinate derivatives
!     deprec  reciprocal polarization Cartesian coordinate derivatives
!     deg     geometric restraint Cartesian coordinate derivatives
!     dex     extra energy term Cartesian coordinate derivatives
!     desave  stored Cartesian coordinate derivatives
!     desmd   extra smd energy term Cartesian coordinate derivatives
!     deamdD  derivatives from aMD on dihedrals
!     deamdP  derivatives from aMD on potential energy
!     de_mpi  forces communication buffer
!     de_ws   work space buffer
!     deW1aMD
!     deW2aMD
!
!     Lambda-dynamics derivatives
!
!     delambda           hamiltonian derivative with respect to lambda (to be sent to colvar)
!     delambdae          hamiltonian derivative with respect to elambda
!     delambdav          hamiltonian derivative with respect to vlambda
!     delambdaesave      stored hamiltonian derivative with respect to elambda
!     delambdavsave      stored hamiltonian derivative with respect to vlambda
!     dlambdaelambda     derivative of elambda with respect to lambda
!     dlambdavlambda     derivative of vlambda with respect to lambda
!
!     Orthogonal Space Random Walk - note x stands for Cartesian coordinates
!     dxdelambda         hamiltonian double derivative with respect to x and lambda (to be sent to colvar)
!     dxdelambdae        hamiltonian double derivative with respect to x and elambda (electrostatic interactions)
!     dxdelambdav        hamiltonian double derivative with respect to x and vlambda (vdw interactions)
!     d2edlambda2         hamiltonian double derivative with respect to lambda (to be sent to colvar)
!     d2edlambdae2        hamiltonian double derivative with respect to elambda (electrostatic interactions)
!     d2edlambdav2        hamiltonian double derivative with respect to vlambda (vdw interactions)
!
!     dotstgrad : flag when the main program is testgrad (communication
!      of the forces one by one)
!
!     cBond       Bonded force communication parameter
!                 (check info_forces implementation )
!     cNBond  Non Bonded force communication parameter
!     cSNBond Short Non Bonded force communication parameter
!     cDef        Default force communication parameter
!
!     ftot_l      use /de_tot/ buffer for partial summation
!     fdebs_l     use /de_tot/ for all bonded buffers (switch)
!     tdes_l      transpose force buffers in a (n,3) shape (switch)
!     dr_stride   inner stride between one set of coordinates of an atom (x -- y -- z)
!     dr_stride3  stride between two consecutive buffers ( deb -- dea )
!     dr_nb0      number of buffers used in dr_buf0
!     dr_nbb      number of buffers used for bonded forces
!     dr_nbnb     number of buffers used for non bonded forces
!     dr_nbnbr    number of buffers used for non bonded reciprocal forces
!     dr_obb      offset buffer for bonded forces
!     dr_obnb     offset buffer for non bonded forces
!     dr_obnbr    offset buffer for non bonded reciprocal forces
!
!
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
   real(r_p),pointer :: det(:,:),dept(:,:),dett(:,:)&
      &,debt(:,:),deat(:,:)
   !DIR$ ATTRIBUTES ALIGN:64:: deg,dex
   real(r_p),pointer :: deg(:,:),dex(:,:)
   real(r_p),pointer :: dmlpot(:,:)

   mdyn_rtyp,pointer,dimension(:) :: derivx,derivy,derivz&
      &,debx,deby,debz,deax,deay,deaz,deubx,deuby,deubz

   !DIR$ ATTRIBUTES ALIGN:64:: desmd,dev,dem,dep,deopd
   mdyn_rtyp,pointer :: desmd(:,:),dev(:,:),der(:,:),dedsp(:,:)&
      &,dem(:,:),dep(:,:),dec(:,:),dect(:,:)
   !DIR$ ATTRIBUTES ALIGN:64:: decrec,demrec,deprec
   real(r_p),pointer :: demrec(:,:),deprec(:,:),decrec(:,:)&
      &,dedsprec(:,:)

   mdyn_rtyp,pointer,dimension(:) :: devx,devy,devz,decx,decy,decz&
      &,demx,demy,demz,depx,depy,depz&
      &,decrecx,decrecy,decrecz,demrecx,demrecy,demrecz&
      &,deprecx,deprecy,deprecz

   !DIR$ ATTRIBUTES ALIGN:64:: deamdD,deamdP,deW1aMD,deW2aMD
   real(r_p),pointer :: deamdD(:,:),deamdP(:,:), deW1aMD(:,:)&
      &,deW2aMD(:,:)
   !DIR$ ATTRIBUTES ALIGN:64:: desave
   mdyn_rtyp,pointer ::desave(:,:)

   real(r_p) delambda,delambdae,delambdav
   real(r_p) delambdaesave,delambdavsave
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
