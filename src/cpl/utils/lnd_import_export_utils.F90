module lnd_import_export_utils

  use shr_kind_mod          , only : r8 => shr_kind_r8
  use shr_infnan_mod        , only : isnan => shr_infnan_isnan
  use shr_sys_mod           , only : shr_sys_abort
  use clm_varctl            , only : iulog
  use decompmod             , only : bounds_type
  use atm2lndType           , only : atm2lnd_type
  use Wateratm2lndBulkType  , only : wateratm2lndbulk_type

  implicit none
  private ! except

  public :: derive_quantities
  public :: check_for_errors
  public :: check_for_nans

!=============================================================================
contains
!=============================================================================

  !===========================================================================

  subroutine derive_quantities( bounds, atm2lnd_inst, wateratm2lndbulk_inst, &
    forc_rainc, forc_rainl, forc_snowc, forc_snowl )

    !-------------------------------------------------------------------------
    ! Convert the input data from the mediator to the land model
    !-------------------------------------------------------------------------

    use clm_varcon, only: rair, o2_molar_const
    use QSatMod, only: QSat

    ! input/output variabes
    type(bounds_type), intent(in) :: bounds  ! bounds
    type(atm2lnd_type), intent(inout) :: atm2lnd_inst ! clm internal input data type
    type(wateratm2lndbulk_type), intent(inout) :: wateratm2lndbulk_inst
    real(r8), intent(in) :: forc_rainc(bounds%begg:bounds%endg)  ! convective rain (mm/s)
    real(r8), intent(in) :: forc_rainl(bounds%begg:bounds%endg)  ! large scale rain (mm/s)
    real(r8), intent(in) :: forc_snowc(bounds%begg:bounds%endg)  ! convective snow (mm/s)
    real(r8), intent(in) :: forc_snowl(bounds%begg:bounds%endg)  ! large scale snow (mm/s)

    ! local variables
    integer :: g  ! indices
    integer :: begg, endg  ! bounds
    real(r8) :: qsat_kg_kg  ! saturation specific humidity (kg/kg)
    real(r8) :: forc_t  ! atmospheric temperature (Kelvin)
    real(r8) :: forc_q  ! atmospheric specific humidity (kg/kg)
    real(r8) :: forc_pbot  ! atmospheric pressure (Pa)
    character(len=*), parameter :: subname='(cpl:utils:derive_quantities)'

    !-------------------------------------------------------------------------

    ! Set bounds
    begg = bounds%begg; endg=bounds%endg

    !--------------------------
    ! Derived quantities
    !--------------------------

    do g = begg, endg
       forc_t    = atm2lnd_inst%forc_t_not_downscaled_grc(g)
       forc_q    = wateratm2lndbulk_inst%forc_q_not_downscaled_grc(g)
       forc_pbot = atm2lnd_inst%forc_pbot_not_downscaled_grc(g)

       atm2lnd_inst%forc_hgt_u_grc(g) = atm2lnd_inst%forc_hgt_grc(g)  !observational height of wind [m]
       atm2lnd_inst%forc_hgt_t_grc(g) = atm2lnd_inst%forc_hgt_grc(g)  !observational height of temperature [m]
       atm2lnd_inst%forc_hgt_q_grc(g) = atm2lnd_inst%forc_hgt_grc(g)  !observational height of humidity [m]

       atm2lnd_inst%forc_vp_grc(g) = forc_q * forc_pbot  / (0.622_r8 + 0.378_r8 * forc_q)

       atm2lnd_inst%forc_rho_not_downscaled_grc(g) = &
            (forc_pbot - 0.378_r8 * atm2lnd_inst%forc_vp_grc(g)) / (rair * forc_t)

       atm2lnd_inst%forc_po2_grc(g) = o2_molar_const * forc_pbot

       atm2lnd_inst%forc_wind_grc(g) = sqrt(atm2lnd_inst%forc_u_grc(g)**2 + atm2lnd_inst%forc_v_grc(g)**2)

       atm2lnd_inst%forc_solar_not_downscaled_grc(g) = &
              atm2lnd_inst%forc_solad_not_downscaled_grc(g,1) &
            + atm2lnd_inst%forc_solai_grc(g,1) &
            + atm2lnd_inst%forc_solad_not_downscaled_grc(g,2) &
            + atm2lnd_inst%forc_solai_grc(g,2)

       wateratm2lndbulk_inst%forc_rain_not_downscaled_grc(g)  = forc_rainc(g) + forc_rainl(g)
       wateratm2lndbulk_inst%forc_snow_not_downscaled_grc(g)  = forc_snowc(g) + forc_snowl(g)

       call QSat(forc_t, forc_pbot, qsat_kg_kg)

       wateratm2lndbulk_inst%forc_rh_grc(g) = 100.0_r8*(forc_q / qsat_kg_kg)
    end do

  end subroutine derive_quantities

  !===========================================================================

  subroutine check_for_errors( bounds, atm2lnd_inst, wateratm2lndbulk_inst )

    ! input/output variabes
    type(bounds_type), intent(in) :: bounds  ! bounds
    type(atm2lnd_type), intent(inout) :: atm2lnd_inst ! clm internal input data type
    type(wateratm2lndbulk_type), intent(inout) :: wateratm2lndbulk_inst

    ! local variables
    integer :: g  ! indices
    integer :: begg, endg  ! bounds
    character(len=*), parameter :: subname='(cpl:utils:check_for_errors)'

    !-------------------------------------------------------------------------

    ! Set bounds
    begg = bounds%begg; endg=bounds%endg

    !--------------------------
    ! Error checks
    !--------------------------

    ! Check that solar, specific-humidity, and LW downward aren't negative
    do g = begg, endg
       if ( atm2lnd_inst%forc_lwrad_not_downscaled_grc(g) <= 0.0_r8 ) then
          call shr_sys_abort( subname//&
               ' ERROR: Longwave down sent from the atmosphere model is negative or zero' )
       end if
       if ( (atm2lnd_inst%forc_solad_not_downscaled_grc(g,1) < 0.0_r8) .or. &
            (atm2lnd_inst%forc_solad_not_downscaled_grc(g,2) < 0.0_r8) .or. &
            (atm2lnd_inst%forc_solai_grc(g,1) < 0.0_r8) .or. &
            (atm2lnd_inst%forc_solai_grc(g,2) < 0.0_r8) ) then
          call shr_sys_abort( subname//&
               ' ERROR: One of the solar fields (indirect/diffuse, vis or near-IR)'// &
               ' from the atmosphere model is negative or zero' )
       end if
       if ( wateratm2lndbulk_inst%forc_q_not_downscaled_grc(g) < 0.0_r8 )then
          call shr_sys_abort( subname//&
               ' ERROR: Bottom layer specific humidty sent from the atmosphere model is less than zero' )
       end if
    end do

    ! Make sure relative humidity is properly bounded
    ! atm2lnd_inst%forc_rh_grc(g) = min( 100.0_r8, atm2lnd_inst%forc_rh_grc(g) )
    ! atm2lnd_inst%forc_rh_grc(g) = max(   0.0_r8, atm2lnd_inst%forc_rh_grc(g) )

  end subroutine check_for_errors

  !=============================================================================

  subroutine check_for_nans(array, fname, begg, direction)
    use GridcellType    , only : grc                

    ! input/output variables
    real(r8)         , intent(in) :: array(:)
    character(len=*) , intent(in) :: fname
    integer          , intent(in) :: begg
    character(len=*) , intent(in) :: direction

    ! local variables
    integer :: i
    !---------------------------------------------------------------------------

    ! Check if any input from mediator or output to mediator is NaN

    if (any(isnan(array))) then
       write(iulog,*) '# of NaNs = ', count(isnan(array))
       write(iulog,*) 'Which are NaNs = ', isnan(array)
       do i = 1, size(array)
          if (isnan(array(i))) then
             write(iulog,*) "NaN found in field ", trim(fname), ' at gridcell index/lon/lat: ',begg+i-1,grc%londeg(begg+i-1),grc%latdeg(begg+i-1)
          end if
       end do
       call shr_sys_abort(' ERROR: One or more of the CTSM cap '//direction//' fields are NaN ' )
    end if
  end subroutine check_for_nans

end module lnd_import_export_utils
