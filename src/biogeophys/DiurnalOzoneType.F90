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
       integer                       :: ntimes = 24
       real(r8), public, allocatable :: o3_anomaly_grc(:,:)    ! o3 anomaly data [grc, ntimes]
  
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
    subroutine Init(this, bounds)
      !
      ! DESCRIPTION:
      ! Initialize ozone data structure
      !
      !
      ! ARGUMENTS:
      class(diurnal_ozone_anom_type), intent(inout) :: this
      type(bounds_type),              intent(in)    :: bounds 
      !-----------------------------------------------------------------------

      call this%InitAllocate(bounds)

    end subroutine Init
  
  
    !-----------------------------------------------------------------------
    subroutine InitAllocate(this, bounds)
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
      !
      ! LOCAL VARIABLES:
      integer  :: begg, endg
      !---------------------------------------------------------------------

      begg = bounds%begg; endg = bounds%endg

      allocate(this%o3_anomaly_grc(begg:endg,1:this%ntimes)); this%o3_anomaly_grc(:,:) = nan

      end subroutine InitAllocate

  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  subroutine Interp(this, bounds, forc_o3, forc_o3_down)
    !
    ! !DESCRIPTION:
    ! Allocate module variables and data structures
    !
    ! !USES:
    use shr_infnan_mod,    only : nan => shr_infnan_nan, assignment(=)
    use clm_time_manager , only : get_curr_time
    !
    ! !ARGUMENTS:
    class(diurnal_ozone_anom_type), intent(in)  :: this   
    type(bounds_type),              intent(in)  :: bounds                       
    real(r8),                       intent(in)  :: forc_o3( bounds%begg: )      ! ozone partial pressure (mol/mol)
    real(r8),                       intent(out) :: forc_o3_down( bounds%begg: ) ! ozone partial pressure, downscaled (mol/mol)
    !
    ! LOCAL VARIABLES:
    integer :: j     ! time stamp to grab
    integer :: yr    ! year
    integer :: mon   ! month
    integer :: day   ! day of month
    integer :: tod   ! time of day (seconds past 0Z)
    integer :: begg, endg
    !-----------------------------------------------------------------------

    begg = bounds%begg; endg = bounds%endg

    ! Get current date/time - we really only need seconds
    call get_curr_date(yr, mon, day, tod)

    !! interpolate here!
   
    ! apply anomaly
    forc_o3_down(begg:endg) = forc_o3(begg:eng)*this%o3_anomaly_grc(begg:endg, j)

    end subroutine Interp

  !-----------------------------------------------------------------------

  end module DiurnalOzoneType
  