module FATESFireFactoryMod

  !---------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Factory to create an instance of fire_method_type. This module figures
  ! out the particular type to return.
  !
  ! !USES:
  use abortutils          , only : endrun
  use shr_log_mod         , only : errMsg => shr_log_errMsg
  use clm_varctl          , only : iulog

  implicit none
  save
  private
  !
  ! !PUBLIC ROUTINES:
  public :: create_fates_fire_data_method  ! create an object of class fates_fire_base_type

  ! These parameters set the ranges of the cases in subroutine
  ! create_fates_fire_data_method. We declare them public in order to
  ! use them as flags elsewhere in the CTSM and FATES-SPITFIRE.
  ! They correspond one-to-one to the fates_spitfire_mode options listed
  ! in bld/namelist_files/namelist_definition_clm4_5.xml
  integer, public, parameter :: no_fire = 0              ! value of no_fire mode
  integer, public, parameter :: scalar_lightning = 1     ! value of scalar_lightning mode
  integer, public, parameter :: lightning_from_data = 2  ! value of lightning_from_data mode
  integer, public, parameter :: successful_ignitions = 3 ! value of successful_ignitions mode
  integer, public, parameter :: anthro_ignitions = 4     ! value of anthro_ignitions mode
  integer, public, parameter :: anthro_suppression = 5   ! value of anthro_suppression mode

  character(len=*), parameter, private :: sourcefile = &
       __FILE__

contains

  !-----------------------------------------------------------------------
  subroutine create_fates_fire_data_method( fates_fire_data_method )
    !
    ! !DESCRIPTION:
    ! Create and return an object of fates_fire_data_method_type.
    ! The particular type is determined based on a namelist parameter.
    !
    ! !USES:
    use clm_varctl, only: fates_spitfire_mode, use_fates_sp, use_fates_ed_st3
    use shr_fire_emis_mod, only : shr_fire_emis_mechcomps_n
    use FATESFireBase,      only: fates_fire_base_type
    use FATESFireNoDataMod, only: fates_fire_no_data_type
    use FATESFireDataMod, only: fates_fire_data_type
    !
    ! !ARGUMENTS:
    class(fates_fire_base_type), allocatable, intent(inout) :: fates_fire_data_method  ! function result
    !
    ! !LOCAL VARIABLES:
    !-----------------------------------------------------------------------

    !
    ! For FATES options that bypass fire...
    !
    if ( use_fates_sp .or. use_fates_ed_st3 )then
      !
      ! Make sure fire-emissions is NOT on
      !
      if ( shr_fire_emis_mechcomps_n > 0 )then
         if ( use_fates_sp )then
            write(iulog,*) "Fire emissions can NOT be on with FATES-SP mode: ", &
                        errMsg(sourcefile, __LINE__)
            call endrun(msg="Fire emission with FATES requires FATES to NOT be in Satellite Phenology (SP) mode" )
         else if ( use_fates_ed_st3 )then
            write(iulog,*) "Fire emissions can NOT be on with FATES ST3 mode: ", &
                        errMsg(sourcefile, __LINE__)
            call endrun(msg="Fire emission with FATES requires FATES to NOT be in Static Stand Structure mode" )
         end if
         ! For unit-testing return with a FATESFireData type, so there isn't a run-time error
         ! Also do the FATESFireData type, as using FATESFireNoData type will fail with an error
         allocate(fates_fire_data_type :: fates_fire_data_method)
         return
      end if
      allocate(fates_fire_no_data_type :: fates_fire_data_method)
    else
      !
      ! For regular FATES options that include fire
      !
      select case (fates_spitfire_mode)

         ! No-fire, scalar-lightning and successful_ignitions ALL do NOT need input data from the base class
         case (no_fire:scalar_lightning)
            allocate(fates_fire_no_data_type :: fates_fire_data_method)
         case (successful_ignitions)
            allocate(fates_fire_no_data_type :: fates_fire_data_method)
         ! Lightning from data, and the anthro types (ignition and suppression) need lightning data from the base class
         case (lightning_from_data)
            allocate(fates_fire_data_type :: fates_fire_data_method)
         case (anthro_ignitions:anthro_suppression)
            allocate(fates_fire_data_type :: fates_fire_data_method)

         case default
            write(iulog,*) 'Unrecognized fates_spitfire_mode option = ', fates_spitfire_mode, ' in: ', &
                           errMsg(sourcefile, __LINE__)
            call endrun(msg="Unknown option for namelist item fates_spitfire_mode:")
            ! For unit-testing, make sure a valid fates_fire_data_method is set and return, otherwise it fails with a seg-fault
            allocate(fates_fire_no_data_type :: fates_fire_data_method)

      end select
      ! -------------------------------------------------------------------------------------------------------
      ! For now we die with a error whenever fire-emissions are turned on -- because this isn't setup in FATES
      !
      if ( fates_spitfire_mode /= no_fire ) then
         if ( shr_fire_emis_mechcomps_n > 0 )then
            write(iulog,*) "Fire emissions can NOT be on with FATES currently: ", &
                        errMsg(sourcefile, __LINE__)
            call endrun(msg="Fire emission with FATES can NOT currently be turned on (see issue #1045)" )
            return
         end if
      end if
      ! -------------------------------------------------------------------------------------------------------
    end if

  end subroutine create_fates_fire_data_method

end module FATESFireFactoryMod
