!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ######################################################################
!     ##                                                                  ##
!     ##  module tors  --  torsional angles within the current structure  ##
!     ##                                                                  ##
!     ######################################################################
!
!
!     tors1   !<1-fold amplitude and phase for each torsional angle
!     wintors1    !<window object corresponding to tors1
!     tors2   !<2-fold amplitude and phase for each torsional angle
!     wintors2    !<window object corresponding to tors2
!     tors3   !<3-fold amplitude and phase for each torsional angle
!     wintors3    !<window object corresponding to tors3
!     tors4   !<4-fold amplitude and phase for each torsional angle
!     wintors4    !<window object corresponding to tors4
!     tors5   !<5-fold amplitude and phase for each torsional angle
!     wintors5    !<window object corresponding to tors5
!     tors6   !<6-fold amplitude and phase for each torsional angle
!     wintors6    !<window object corresponding to tors6
!     ntors   !<total number of torsional angles in the system
!     ntorsloc   !<local number of torsional angles in the system
!     nbtors   !<number of torsional angles before each atom
!     winnbtors    !<window object corresponding to nbtors
!     itors   !<numbers of the atoms in each torsional angle
!     winitors    !<window object corresponding to itors
!
!
module tors
   implicit none
   integer :: ntors !<total number of torsional angles in the system
   integer :: ntorsloc !<local number of torsional angles in the system
   integer, pointer :: nbtors(:) !<number of torsional angles before each atom
   integer, pointer :: itors(:,:) !<numbers of the atoms in each torsional angle
   integer :: winnbtors !<window object corresponding to nbtors
   integer :: winitors !<window object corresponding to ibtors
   integer :: wintors1 !<window object corresponding to tors1
   integer :: wintors2 !<window object corresponding to tors2
   integer :: wintors3 !<window object corresponding to tors3
   integer :: wintors4 !<window object corresponding to tors4
   integer :: wintors5 !<window object corresponding to tors5
   integer :: wintors6 !<window object corresponding to tors6
   real*8, pointer :: tors1(:,:) !<1-fold amplitude and phase for each torsional angle
   real*8, pointer :: tors2(:,:) !<2-fold amplitude and phase for each torsional angle
   real*8, pointer :: tors3(:,:) !<3-fold amplitude and phase for each torsional angle
   real*8, pointer :: tors4(:,:) !<4-fold amplitude and phase for each torsional angle
   real*8, pointer :: tors5(:,:) !<5-fold amplitude and phase for each torsional angle
   real*8, pointer :: tors6(:,:) !<6-fold amplitude and phase for each torsional angle
   save
end
