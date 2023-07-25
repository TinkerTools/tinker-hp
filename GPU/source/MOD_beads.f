c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #####################################################################
c     ##                                                                 ##
c     ##  module beads   --  pimd variables                              ##
c     ##                                                                 ##
c     #####################################################################
c
c
#include "tinker_precision.h"
      module beads
      use dcdmod
      implicit none
      logical :: path_integral_md=.FALSE.
      integer :: nbeads,nbeadsloc   

      character(len=80) :: save_beads

      logical :: contract
      logical :: centroid_longrange
      logical :: centroid_recip
      logical :: polar_allbeads
      logical :: cay_correction
      integer :: nbeads_ctr, nbeadsloc_ctr  

      logical :: pi_comm_report
      integer :: pi_comm_unit
      logical :: PITIMER
      
      integer :: max_polymer_index=0
      integer :: ipolymer_loaded=-1
      integer :: ibead_loaded=-1
      logical :: add_buffer_pi

      integer, allocatable :: ilocpi_beg(:)
      integer, allocatable :: ilocpi_end(:)
      integer, allocatable :: nlocpis(:)
      integer :: nlocpi,nlocpiproc

      real(r_p) :: pi_start_isobaric
      real(r_p) :: lambda_trpmd
      logical :: default_lambda_trpmd=.TRUE.
      real(r_p) :: gyrpi
      real(r_p) :: temppi,temppi_cl,temppi_centroid
      real(r_p) :: epotpi_loc,etotpi_loc
      real(r_p) :: eksumpi_loc,ekinpi_loc(3,3)
      real(r_p) :: ekdynpi, epotpi, eintrapi, einterpi
      real(r_p) :: ekprim,ekvir,presvir
      real(r_p) :: ekvir_fast
      real(r_p) :: ekcentroid
      real(r_p) :: eintrapi_loc
      real(r_p) :: einterpi_loc
      real(r_p) :: virpi(3,3),ekinpi(3,3),ekinpi_fast(3,3)
      real(r_p) :: stresspi(3,3)
      real(r_p) :: dedvpi,dedvintrapi,dedvinterpi
      real(r_p) :: dippi(3),dipindpi(3)
      logical :: pot_gathered=.TRUE.

      real(r_p), allocatable :: bufferpi(:)

      
      real(8), allocatable :: eigmat(:,:)
      real(8), allocatable ::  omkpi(:)
      
      real(8), allocatable :: contractor_mat(:,:)

      TYPE POLYMER_COMM_TYPE
        integer :: nbeads
        integer :: nbeadsproc
        integer,allocatable :: nbeadsloc(:)
        integer,allocatable :: ibead_beg(:)
        integer,allocatable :: ibead_end(:)
        integer,allocatable :: iproc(:)
        integer :: index

        logical :: allocated=.FALSE.

        real(r_p), allocatable :: pos(:,:,:),vel(:,:,:)
        real(r_p), allocatable :: eigpos(:,:,:),eigvel(:,:,:)
        real(r_p), allocatable :: forces(:,:,:)
        real(r_p), allocatable :: eigforces(:,:,:)
        real(r_p), allocatable :: forces_slow(:,:,:)

        real(r_p), allocatable :: epot(:)
        real(r_p), allocatable :: vir(:,:,:)

        real(r_p), allocatable :: dip(:,:)
        real(r_p), allocatable :: dipind(:,:)

        real(r_p) :: gyr

        type(dcdinfo_t),allocatable :: dcdinfo(:)
      END TYPE     

      TYPE(POLYMER_COMM_TYPE) :: polymer,polymer_ctr

      save

      interface

      module subroutine allocate_polymer(polymer,nu,allocate_vel
     &      ,allocate_forces,allocate_slow)
      implicit none
      TYPE(POLYMER_COMM_TYPE), intent(inout) :: polymer
      integer, intent(in) :: nu
      LOGICAL, intent(in) :: allocate_vel, allocate_slow
      LOGICAL, intent(in) :: allocate_forces
      end subroutine allocate_polymer

      module subroutine deallocate_polymer(polymer)
      implicit none
      TYPE(POLYMER_COMM_TYPE), intent(inout) :: polymer
      end subroutine deallocate_polymer

      module subroutine initialize_normal_modes(polymer,verbose)
      implicit none
      TYPE(POLYMER_COMM_TYPE), intent(inout) :: polymer
      logical, intent(in) :: verbose
      end subroutine initialize_normal_modes

      module subroutine update_nlocpi(nloc)
      implicit none
      integer, intent(in) :: nloc
      end subroutine update_nlocpi

      module subroutine allocpi()
      implicit none
      end subroutine allocpi

      module subroutine update_normal_modes_pi(polymer)
      implicit none
      type(POLYMER_COMM_TYPE), intent(inout) :: polymer
      end subroutine update_normal_modes_pi

      module subroutine set_eigforces_pi(polymer,forces)
      implicit none
      type(POLYMER_COMM_TYPE), intent(inout) :: polymer
      real(r_p), intent(in) :: forces(:,:,:)
      end subroutine set_eigforces_pi

      module subroutine update_direct_space_pi(polymer)
      implicit none
      type(POLYMER_COMM_TYPE), intent(inout) :: polymer
      end subroutine update_direct_space_pi

      module subroutine contract_polymer(polymer,polymer_ctr)
      implicit none
      type(POLYMER_COMM_TYPE), intent(in) :: polymer
      type(POLYMER_COMM_TYPE), intent(inout) :: polymer_ctr
      end subroutine contract_polymer

      module subroutine project_forces_ctr(polymer,polymer_ctr,slow)
      implicit none
      type(POLYMER_COMM_TYPE), intent(inout) :: polymer,polymer_ctr
      logical, intent(in), optional :: slow
      end subroutine project_forces_ctr

      module subroutine load_bead(polymer,ibead,force_load)
      implicit none
      type(POLYMER_COMM_TYPE), intent(in) :: polymer
      integer, intent(in) :: ibead
      logical, intent(in), optional :: force_load
      end subroutine load_bead

      module subroutine stopwatchpi(timer,barrier,reset)
      implicit none
      real(8), intent(inout) :: timer
      logical, intent(in) :: barrier,reset
      end subroutine stopwatchpi

      module subroutine kinetic_pi_fast(polymer)
      implicit none
      type(POLYMER_COMM_TYPE), intent(inout) :: polymer
      end subroutine kinetic_pi_fast


      module subroutine kinetic_pi(polymer)
      implicit none
      type(POLYMER_COMM_TYPE), intent(inout) :: polymer
      end subroutine kinetic_pi

      module subroutine reduce_observables_pi(polymer)
      implicit none
      type(POLYMER_COMM_TYPE), intent(in) :: polymer
      end subroutine reduce_observables_pi

      module subroutine reset_observables(init)
      implicit none
      logical, intent(in) :: init
      end subroutine reset_observables

      module subroutine read_bead_config()
      implicit none
      end subroutine read_bead_config


      end interface

      end module beads
