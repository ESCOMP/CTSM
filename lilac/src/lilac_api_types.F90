module lilac_api_types

   implicit none

   use lilac_constants, only : STRING_128
   
contains

   type :: lilac_init_data_t
      character(len=STRING_32) :: component_name
      integer :: mpicom_lilac
      integer :: mpicom_component
      integer :: output_unit_lilac
      integer :: output_unit_global_shared ! this should be the same for all instances of lilac!
      integer :: output_unit_component

   end type lilac_init_data_t

   type :: lilac_clock_data_t
      logical :: calendar_is_leap
      integer :: start_year
      integer :: start_month
      integer :: start_day
      integer :: start_second ! seconds since midnight

      integer :: stop_year
      integer :: stop_month
      integer :: stop_day
      integer :: stop_second ! seconds since midnight

      integer :: time_step_seconds
   end type lilac_clock_data_t


   type :: lilac_exchange_fields_t
      character(len=STRING_128) :: long_name
      character(len=STRING_128) :: short_name
      character(len=STRING_128) :: field_name
      character(len=STRING_128) :: units
      integer :: field_type
   end type lilac_exchange_fields_t
   
end module lilac_api_types
