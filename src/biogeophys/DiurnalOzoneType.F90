module DiurnalOzoneType

    !-----------------------------------------------------------------------
    ! !DESCRIPTION:
    ! Sets up type for converting multi-day ozone input to sub-daily values using an input
    ! ozone anomaly file
    !
    ! !USES:
  #include "shr_assert.h"
    use shr_kind_mod           , only : r8 => shr_kind_r8
    use decompMod              , only : bounds_type
    use clm_varcon             , only : spval
    use clm_varctl             , only : iulog
    use abortutils             , only : endrun
    use shr_ozone_coupling_mod , only : atm_ozone_frequency_unset, atm_ozone_frequency_subdaily, & 
                                        atm_ozone_frequency_multiday_average
    implicit none
    save
    private
  
    ! !PUBLIC TYPES:
    type, public :: diurnal_ozone_anom_type
       private
       ! Private data members
       integer               :: ozone_input_frequency  ! Which ozone input frequency are we receiving?
       integer               :: ntimes 
       real(r8), allocatable :: sec(:)                 ! seconds of day (size ntimes)
       real(r8), allocatable :: o3_anomaly_grc(:,:)    ! o3 anomaly data [grc, ntimes]
  
     contains
       ! Public routines
       procedure, public :: Init
       procedure, public :: Interp

    end type diurnal_ozone_anom_type
  
    ! !PRIVATE TYPES:
    integer, parameter :: ozone_frequency_unset = 0
    integer, parameter :: ozone_frequency_subdaily = 1
    integer, parameter :: ozone_frequency_multiday_average = 2

    character(len=*), parameter, private :: sourcefile = &
    __FILE__


  contains
  
    ! ========================================================================
    ! Infrastructure routines (initialization, restart, etc.)
    ! ========================================================================
  
    !-----------------------------------------------------------------------
    subroutine Init(this)
      !
      ! !DESCRIPTION:
      ! Initialize ozone data structure
      !
      !
      ! !ARGUMENTS:
      class(diurnal_ozone_anom_type), intent(inout) :: this
      !-----------------------------------------------------------------------
  
      if (atm_ozone_frequency_val == atm_ozone_frequency_unset) then
         this%ozone_input_frequency = ozone_frequency_unset
      else if (atm_ozone_frequency_val == atm_ozone_frequency_subdaily) then
         this%ozone_input_frequency = ozone_frequency_subdaily
      else if (atm_ozone_frequency_val == atm_ozone_frequency_multiday_average) then
        this%ozone_input_frequency = ozone_frequency_multiday_average
      else
         call endrun('unknown ozone frequency')
      end if
  
    end subroutine Init
  
  
    !-----------------------------------------------------------------------

  
  end module DiurnalOzoneType
  