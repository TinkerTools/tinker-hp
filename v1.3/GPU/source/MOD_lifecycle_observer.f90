!
!     Small lifecycle event ledger used by focused teardown contracts.
!
module lifecycle_observer
   implicit none
   private

   integer, parameter, public :: LifecycleCflux = 1
   integer, parameter, public :: LifecyclePolar = 2
   integer, parameter, public :: LifecycleTmatxb = 3
   integer, parameter, public :: LifecycleReduction = 4
   integer, parameter, public :: LifecyclePme = 5
   integer, parameter, public :: LifecycleDevice = 6
   integer, parameter, public :: LifecycleTopology = 7
   integer, parameter, public :: LifecycleSharedWindow = 8
   integer, parameter, public :: LifecycleTicket = 9
   integer, parameter :: LifecycleCount = LifecycleTicket
   integer, save :: event_counts(LifecycleCount) = 0
   integer, save :: event_serial(LifecycleCount) = 0
   integer, save :: next_serial = 0

   public :: reset_lifecycle_events, record_lifecycle_event, lifecycle_event_count, &
      &lifecycle_event_order

contains

   subroutine reset_lifecycle_events()
      event_counts = 0
      event_serial = 0
      next_serial = 0
   end subroutine reset_lifecycle_events

   subroutine record_lifecycle_event(resource)
      integer, intent(in) :: resource

      if (resource < 1 .or. resource > LifecycleCount) error stop 'invalid lifecycle resource'
      next_serial = next_serial + 1
      event_counts(resource) = event_counts(resource) + 1
      event_serial(resource) = next_serial
   end subroutine record_lifecycle_event

   integer function lifecycle_event_count(resource)
      integer, intent(in) :: resource

      if (resource < 1 .or. resource > LifecycleCount) error stop 'invalid lifecycle resource'
      lifecycle_event_count = event_counts(resource)
   end function lifecycle_event_count

   integer function lifecycle_event_order(resource)
      integer, intent(in) :: resource

      if (resource < 1 .or. resource > LifecycleCount) error stop 'invalid lifecycle resource'
      lifecycle_event_order = event_serial(resource)
   end function lifecycle_event_order

end module lifecycle_observer
