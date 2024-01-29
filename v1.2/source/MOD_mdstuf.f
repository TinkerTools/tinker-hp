c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ###################################################################
c     ##                                                               ##
c     ##  module mdstuf  --  control of molecular dynamics trajectory  ##
c     ##                                                               ##
c     ###################################################################
c
c
c     nfree       total number of degrees of freedom for a system
c     irest       steps between removal of COM inertia (0=no removal)
c     bmnmix      mixing coefficient for use with Beeman integrator
c     dorest      logical flag to remove center of mass inertia
c     velsave     logical flag to save velocity vector components
c     frcsave     logical flag to save force vector components
c     uindsave    logical flag to save induced atomic dipoles
c     integrate   type of molecular dynamics integration algorithm
c     mts         logical flag regarding multi-timestep integration
c
c
      module mdstuf
      implicit none
      integer nfree,irest
      integer bmnmix
      logical dorest
      logical velsave
      logical frcsave
      logical uindsave
      logical mts
      character*11 integrate
      save
      end
