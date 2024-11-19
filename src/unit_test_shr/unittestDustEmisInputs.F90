module unittestDustEmisInputs

  ! unit testing dust emission inputs

  use unittestSubgridMod, only : bounds, unittest_subgrid_teardown
  use unittestSimpleSubgridSetupsMod, only : setup_single_veg_patch
  use clm_varpar, only : nlevsoi, nlevgrnd, clm_varpar_init
  use clm_varctl, only : soil_layerstruct_predefined, create_crop_landunit, use_crop, create_crop_landunit
  use shr_kind_mod , only : r8 => shr_kind_r8
  use unittestFilterBuilderMod, only : filter_from_range
  use atm2lndType, only : atm2lnd_type, atm2lnd_params_type
  use SoilStateType, only : soilstate_type
  use CanopyStateType, only : canopystate_type
  use TemperatureType, only : temperature_type
  use WaterType, only : water_type
  use FrictionVelocityMod, only : frictionvel_type
  use unittestWaterTypeFactory, only : unittest_water_type_factory_type
  use SoilStateInitTimeConstMod, only : ThresholdSoilMoistZender2003, MassFracClay
  use SoilStateInitTimeConstMod, only : MassFracClayLeung2023
  use abortutils, only : endrun

  implicit none
  private

  type, public :: unittest_dust_emis_input_type
     integer, allocatable :: filter_nolakep(:)      ! non-lake filter (patches)
     integer :: num_nolakep              ! number of patches in non-lake filter
     type(atm2lnd_type) :: atm2lnd_inst
     type(soilstate_type) :: soilstate_inst
     type(canopystate_type) :: canopystate_inst
     type(temperature_type), private :: temperature_inst
     type(unittest_water_type_factory_type), private :: water_factory
     type(water_type) :: water_inst
     type(frictionvel_type) :: frictionvel_inst
   contains
     procedure, public :: setUp
     procedure, public :: tearDown
     procedure, private:: setupSoilState
     procedure, public :: create_atm2lnd
     procedure, public :: create_fv
     procedure, public :: print_values
  end type unittest_dust_emis_input_type

contains

  !-----------------------------------------------------------------------

  subroutine setUp(this)
    use ColumnType, only : col
    class(unittest_dust_emis_input_type), intent(inout) :: this
    ! Allocate and initiatlize the test object for input objects dust-emission needs

    character(len=5) :: NLFilename = 'none'
    real(r8), allocatable :: urb_em(:)
    real(r8), allocatable :: exice_init_conc_col(:)
    integer :: begl, endl, begc, endc
    integer :: c
    type(atm2lnd_params_type) :: atm2lnd_params
    integer, parameter :: snl = 3

    !-----------------------------------------------------------------------
    ! Settings needed for clm_varpar
    soil_layerstruct_predefined = '20SL_8.5m'
    create_crop_landunit = .true.
    use_crop = .false.
    call clm_varpar_init( actual_maxsoil_patches=17, surf_numpft=15, surf_numcft=2, actual_nlevurb=5 )
    ! Water object initialization -- before subgrid
    call this%water_factory%init()
    call this%water_factory%setup_before_subgrid( &
        my_nlevsoi = nlevsoi, &
        nlevgrnd_additional = nlevgrnd - nlevsoi, &
        my_nlevsno = snl)
    ! Setup the subgrid structure for a single bare-soil patch
    call setup_single_veg_patch( pft_type=0 )
    begl = bounds%begl
    endl = bounds%endl
    begc = bounds%begc
    endc = bounds%endc

    ! Create the nolake filter
    call filter_from_range(start=bounds%begp, end=bounds%endp, numf=this%num_nolakep, filter=this%filter_nolakep)
    ! atmosphere to land parameter object
    atm2lnd_params = atm2lnd_params_type( repartition_rain_snow = .false., &
                                          glcmec_downscale_longwave = .false., &
                                          lapse_rate = 0.01_r8 &  ! arbitrary (this is unused for these tests)
    )

    call this%atm2lnd_inst%InitForTesting(bounds, atm2lnd_params)
    ! Water and soil state -- after the subgrid setup
    call this%water_factory%setup_after_subgrid(snl = snl)
    call this%setupSoilState( )   ! This needs to happen before the water_type object creation
    call this%water_factory%create_water_type(this%water_inst, watsat_col=this%soilstate_inst%watsat_col)
    ! Canopy state, friction velocity, and temperature state ojects
    call this%canopystate_inst%SetNMLForTesting()
    call this%canopystate_inst%Init(bounds)
    call this%frictionvel_inst%InitForTesting(bounds)
    allocate( urb_em(begl:endl) )
    urb_em(begl:endl) = 0.99_r8  ! Arbitrary won't matter here
    allocate( exice_init_conc_col(begc:endc) )
    exice_init_conc_col(begc:endc) = 0.0_r8 ! zero, so it doesn't affect anything.
    call this%temperature_inst%Init(bounds,           &
                               em_roof_lun=urb_em(begl:endl), &
                               em_wall_lun=urb_em(begl:endl), &
                               em_improad_lun=urb_em(begl:endl), &
                               em_perroad_lun=urb_em(begl:endl), &
                               is_simple_buildtemp=.true., is_prog_buildtemp=.false., &
                               exice_init_conc_col = exice_init_conc_col(begc:endc))
    deallocate( urb_em )
  end subroutine setUp

  !-----------------------------------------------------------------------

  subroutine tearDown(this)
    class(unittest_dust_emis_input_type), intent(inout) :: this

    call this%water_factory%teardown(this%water_inst)
    call unittest_subgrid_teardown()
    call this%atm2lnd_inst%Clean()
    deallocate( this%filter_nolakep )
  end subroutine tearDown

  !-----------------------------------------------------------------------

  subroutine setupSoilState(this)
    !
    ! !DESCRIPTION:
    ! Sets up the external environment used by Dust emissions - i.e., things accessed via
    ! 'use' statements.
    !
    ! Assumes nlevgrnd and nlevsoi have been set, and that all necessary subgrid setup has
    ! been completed.
    !
    use ColumnType, only : col
    use GridcellType, only : grc
    use shr_dust_emis_mod , only : is_dust_emis_zender, is_dust_emis_leung
    class(unittest_dust_emis_input_type), intent(in) :: this
    !
    integer :: c,j
    real(r8), parameter :: clay = 10.0_r8

    !-----------------------------------------------------------------------
    ! Setting the soil structrure is needed for water_state types
    col%dz(:,1:nlevgrnd) = 1.0_r8
    do j = 1, nlevgrnd
       do c = bounds%begc, bounds%endc
          col%z(c,j) = sum(col%dz(c,1:j-1)) + 0.5_r8*col%dz(c,j)
       end do
    end do

    call this%soilstate_inst%Init(bounds)
    do c = bounds%begc, bounds%endc
       this%soilstate_inst%watsat_col(c,:) = 0.05_r8 * (c - bounds%begc + 1)
    end do
    ! These are needed for dust emissions initialization
    do c = bounds%begc, bounds%endc
       this%soilstate_inst%gwc_thr_col(c) = ThresholdSoilMoistZender2003( clay )
       if ( is_dust_emis_zender() )then
          this%soilstate_inst%mss_frc_cly_vld_col(c) = MassFracClay( clay )
       else if ( is_dust_emis_leung() )then
          this%soilstate_inst%mss_frc_cly_vld_col(c) = MassFracClayLeung2023( clay )
       else
          call endrun("ERROR: do NOT know about this dust_emis_method")
       end if
    end do

  end subroutine setupSoilState

  !-----------------------------------------------------------------------

  subroutine create_atm2lnd(this, forc_t, forc_pbot, forc_rho )
    ! Initializes some fields needed for dust emissions in this%atm2lnd_inst, and sets
    ! forcing fields based on inputs. Excluded inputs are given a default value
    class(unittest_dust_emis_input_type), intent(inout) :: this
    real(r8), intent(in), optional :: forc_t(:)
    real(r8), intent(in), optional :: forc_pbot(:)
    real(r8), intent(in), optional :: forc_rho(:)

    real(r8), parameter :: forc_t_default = 301._r8
    real(r8), parameter :: forc_pbot_default = 100000._r8
    real(r8), parameter :: forc_rho_default = 1.1_r8
    ! ------------------------------------------------------------------------

    if (present(forc_t)) then
       this%atm2lnd_inst%forc_t_downscaled_col(bounds%begc:bounds%endc) = forc_t(:)
    else
       this%atm2lnd_inst%forc_t_downscaled_col(bounds%begc:bounds%endc) = forc_t_default
    end if

    if (present(forc_pbot)) then
       this%atm2lnd_inst%forc_pbot_downscaled_col(bounds%begc:bounds%endc) = forc_pbot(:)
    else
       this%atm2lnd_inst%forc_pbot_downscaled_col(bounds%begc:bounds%endc) = forc_pbot_default
    end if

    if (present(forc_rho)) then
       this%atm2lnd_inst%forc_rho_downscaled_col(bounds%begc:bounds%endc) = forc_rho(:)
    else
       this%atm2lnd_inst%forc_rho_downscaled_col(bounds%begc:bounds%endc) = forc_rho_default
    end if

  end subroutine create_atm2lnd

  !-----------------------------------------------------------------------

  subroutine create_fv(this, fv, u10, ram1, obu)
    ! Initializes some fields needed for dust emissions in this%frictionvel_inst, and sets
    ! fields based on inputs. Excluded inputs are given a default value
    class(unittest_dust_emis_input_type), intent(inout) :: this
    real(r8), intent(in), optional :: fv
    real(r8), intent(in), optional :: u10
    real(r8), intent(in), optional :: ram1
    real(r8), intent(in), optional :: obu

    real(r8), parameter :: fv_default = 2.0_r8
    real(r8), parameter :: u10_default = 4._r8
    real(r8), parameter :: ram1_default = 200._r8
    real(r8), parameter :: obu_default = -100._r8
    ! ------------------------------------------------------------------------

    if (present(fv)) then
       this%frictionvel_inst%fv_patch(bounds%begp:bounds%endp) = fv
    else
       this%frictionvel_inst%fv_patch(bounds%begp:bounds%endp) = fv_default
    end if

    if (present(u10)) then
       this%frictionvel_inst%u10_patch(bounds%begp:bounds%endp) = u10
    else
       this%frictionvel_inst%u10_patch(bounds%begp:bounds%endp) = u10_default
    end if

    if (present(ram1)) then
       this%frictionvel_inst%ram1_patch(bounds%begp:bounds%endp) = ram1
    else
       this%frictionvel_inst%ram1_patch(bounds%begp:bounds%endp) = ram1_default
    end if

    if (present(obu)) then
       this%frictionvel_inst%obu_patch(bounds%begp:bounds%endp) = obu
    else
       this%frictionvel_inst%obu_patch(bounds%begp:bounds%endp) = obu_default
    end if

  end subroutine create_fv

  !-----------------------------------------------------------------------

  subroutine print_values(this)
    use LandunitType, only : lun
    use PatchType, only : patch
    class(unittest_dust_emis_input_type), intent(inout) :: this
    integer :: p, c, l

    do l = bounds%begl, bounds%endl
       print *, 'landunit type= ', lun%itype(l)
    end do
    do c = bounds%begc, bounds%endc
       print *, 'watsat = ', this%soilstate_inst%watsat_col(c,1)
       print *, 'h2osoi_vol = ', this%water_inst%waterstatebulk_inst%h2osoi_vol_col(c,1)
       print *, 'frac_sno = ', this%water_inst%waterdiagnosticbulk_inst%frac_sno_col(c)
       print *, 'mss_frac_clay_vld = ', this%soilstate_inst%mss_frc_cly_vld_col(c)
    end do
    do p = bounds%begp, bounds%endp
       print *, 'patch type= ', patch%itype(p)
       print *, 'patch weight= ', patch%wtgcell(p)
       print *, 'patch active= ', patch%active(p)
       print *, 'tlai = ', this%canopystate_inst%tlai_patch(p)
       print *, 'tsai = ', this%canopystate_inst%tsai_patch(p)
    end do
  end subroutine print_values

  !-----------------------------------------------------------------------

end module unittestDustEmisInputs