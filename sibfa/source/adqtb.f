

c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #################################################################################
c     ##                                                                             ##
c     ##  module ADQTB  --  parameters and arrays for Adaptive Quantum Thermal Bath  ##
c     ##                                                                             ##
c     #################################################################################
c
c     nseg : length in time steps of the segments for correlation
c     computations
c     vad : past velocities used for correlation computations
c     fad : past random forces used for correlation computations
c     skipseg : number to seg to be skipped (because the system has not
c     thermalized yet) before to start the adQTB.
c     ntype :  number of different type in the system
c     typemax : maximum of different types in the system
c     adqtb_type : array with the type of each atoms
c     coo_fact_qtb : value of the corrected factor while computing the
c     corrected kinetic and pressure values
c
c     Literature reference: 
c     Mangaud et al
c
      module adqtb
      use qtb
      implicit none
      real*8 A_gamma
      real*8 corr_fact_qtb
      real*8, allocatable :: fad(:,:,:)
      real*8, allocatable :: gamma_type(:,:)
      real*8, allocatable :: mCvv_average_type(:,:)
      real*8, allocatable :: Cvf_average_type(:,:)
      real*8, allocatable :: Cff_average_type(:,:)
      real*8, allocatable :: dFDR_average_type(:,:)
      integer, allocatable :: adqtb_type(:)
      integer, allocatable :: ntype(:)   
      integer :: typemax
      real*8, allocatable :: vad_piston(:),fad_piston(:)
      real*8, allocatable :: gamma_piston(:)
      real*8, allocatable :: mCvv_average_piston(:)
      real*8, allocatable :: Cvf_average_piston(:)
      real*8, allocatable :: Cff_average_piston(:)
      real*8, allocatable :: dFDR_average_piston(:)
      real*8 A_gamma_piston
      save
      end
