!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     #####################################################################
!     ##                                                                 ##
!     ##  module imptor  --  improper torsions in the current structure  ##
!     ##                                                                 ##
!     #####################################################################
!
!
!     itors1   !<1-fold amplitude and phase for each improper torsion
!     winitors1 !<window corresponding to itors1
!     itors2   !<2-fold amplitude and phase for each improper torsion
!     winitors2 !<window corresponding to itors2
!     itors3   !<3-fold amplitude and phase for each improper torsion
!     winitors3 !<window corresponding to itors3
!     nitors   !<total number of improper torsional angles in the system
!     winbimptors !<window corresponding to nbimptors
!     iitors   !<numbers of the atoms in each improper torsional angle
!     winiitors !<window corresponding to iitors
!
!     nbimptor !<number of improper torsions before each atom
!
!
module imptor
   implicit none
   integer :: nitors
   integer :: nitorsloc
   integer :: winiitors !<window corresponding to iitors
   integer :: winnbimptor !<window corresponding to nbimptors
   integer :: winitors1 !<window corresponding to itors1
   integer :: winitors2 !<window corresponding to itors2
   integer :: winitors3 !<window corresponding to itors3
   integer, pointer :: iitors(:,:) !<numbers of the atoms in each improper torsional angle
   integer, pointer :: nbimptor(:) !<number of improper torsions before each atom
   real*8, pointer :: itors1(:,:) !<1-fold amplitude and phase for each improper torsion
   real*8, pointer :: itors2(:,:) !<2-fold amplitude and phase for each improper torsion
   real*8, pointer :: itors3(:,:) !<3-fold amplitude and phase for each improper torsion
   save
end
