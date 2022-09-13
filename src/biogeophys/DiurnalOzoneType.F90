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

    implicit none
    save
    private
  
    ! !PUBLIC TYPES:
    type, public :: diurnal_ozone_anom_type
       private
       ! Private data members
       integer               :: ntimes 
       real(r8), allocatable :: sec(:)                 ! seconds of day (size ntimes)
       real(r8), allocatable :: o3_anomaly_grc(:,:)    ! o3 anomaly data [grc, ntimes]
  
     contains
       ! Public routines
       procedure, public  :: Init
       procedure, private :: InitAllocate
       !procedure, public  :: Interp - will add this eventually

    end type diurnal_ozone_anom_type
  
    character(len=*), parameter, private :: sourcefile = &
    __FILE__


  contains
  
    ! ========================================================================
    ! Infrastructure routines (initialization, etc.)
    ! ========================================================================
  
    !-----------------------------------------------------------------------
    subroutine Init(this, bounds, nsec)
      !
      ! DESCRIPTION:
      ! Initialize ozone data structure
      !
      !
      ! ARGUMENTS:
      class(diurnal_ozone_anom_type), intent(inout) :: this
      type(bounds_type),              intent(in)    :: bounds 
      integer,                        intent(in)    :: nsec
      !-----------------------------------------------------------------------

      call this%InitAllocate(bounds, nsec)
  
    end subroutine Init
  
  
    !-----------------------------------------------------------------------
    subroutine InitAllocate(this, bounds, nsec)
      !
      ! !DESCRIPTION:
      ! Allocate module variables and data structures
      !
      ! !USES:
      use shr_infnan_mod, only: nan => shr_infnan_nan, assignment(=)
      !
      ! !ARGUMENTS:
      class(diurnal_ozone_anom_type), intent(inout) :: this
      type(bounds_type),              intent(in)    :: bounds
      integer,                        intent(in)    :: nsec 
      !
      ! LOCAL VARIABLES:
      integer  :: begg, endg
      !---------------------------------------------------------------------

      begg = bounds%begg; endg = bounds%endg

      allocate(this%o3_anomaly_grc(begg:endg,1:nsec))
      this%o3_anomaly_grc(:,:) = nan

      end subroutine InitAllocate

  !-----------------------------------------------------------------------

  end module DiurnalOzoneType
  