module SoilWaterMovementMod

#include "shr_assert.h"

  !-----------------------------------------------------------------------
  ! DESCRIPTION
  ! module contains different subroutines to couple soil and root water interactions
  !
  ! created by Jinyun Tang, Mar 12, 2014
  use shr_log_mod    , only : errMsg => shr_log_errMsg
  use shr_kind_mod      , only : r8 => shr_kind_r8
  use shr_sys_mod         , only : shr_sys_flush
  use clm_instMod    , only : clm_fates
 
  !
  implicit none
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: SoilWater            ! Calculate soil hydrology   
  public :: init_soilwater_movement
  private :: soilwater_zengdecker2009
  private :: soilwater_moisture_form
!  private :: soilwater_mixed_form
!  private :: soilwater_head_form
  private :: compute_hydraulic_properties
  private :: compute_moisture_fluxes_and_derivs
  private :: compute_RHS_moisture_form
  private :: compute_LHS_moisture_form
  private :: compute_qcharge
  private :: IceImpedance
  private :: TridiagonalCol
  !
  ! The following is only public for the sake of unit testing; it should not be called
  ! directly by CLM code outside this module
  public :: BaseflowSink
  public :: use_aquifer_layer
  !
  ! !PUBLIC DATA MEMBERS:
!  logical, public :: use_bedrock = .false. ! true => run with spatially variable soil depth

  ! !PRIVATE DATA MEMBERS:

  ! Solution method 
  integer, parameter :: zengdecker_2009 = 0
  integer, parameter :: moisture_form = 1
  integer, parameter :: mixed_form = 2
  integer, parameter :: head_form = 3

  ! Boundary conditions
  integer, parameter :: bc_head  = 0
  integer, parameter :: bc_flux  = 1
  integer, parameter :: bc_zero_flux  = 2
  integer, parameter :: bc_waterTable = 3
  integer, parameter :: bc_aquifer    = 4

  ! Soil hydraulic properties
  integer, parameter :: soil_hp_clapphornberg_1978=0
  integer, parameter :: soil_hp_vanGenuchten_1980=1

  real(r8),parameter :: m_to_mm = 1.e3_r8 !convert meters to mm

  integer :: soilwater_movement_method    ! method for solving richards equation
  integer :: upper_boundary_condition     ! named variable for the boundary condition
  integer :: lower_boundary_condition     ! named variable for the boundary condition

  ! Adaptive time stepping algorithmic control parameters
  real(r8) :: dtmin             ! minimum time step length (seconds)
  real(r8) :: verySmall         ! a very small number: used to check for sub step completion
  real(r8) :: xTolerUpper       ! tolerance to halve length of substep
  real(r8) :: xTolerLower       ! tolerance to double length of substep
  integer  :: expensive
  integer  :: inexpensive
  integer  :: flux_calculation

  character(len=*), parameter, private :: sourcefile = &
       __FILE__

  !-----------------------------------------------------------------------

contains

!#1
  !-----------------------------------------------------------------------
  subroutine init_soilwater_movement()
    !
    !DESCRIPTION
    !specify method for doing soil&root water interactions
    !
    ! !USES:
    use abortutils      , only : endrun   
    use fileutils       , only : getavu, relavu
    use spmdMod         , only : mpicom, masterproc
    use shr_mpi_mod     , only : shr_mpi_bcast
    use clm_varctl      , only : iulog, use_bedrock
    use controlMod      , only : NLFilename
    use clm_nlUtilsMod  , only : find_nlgroup_name

    ! !ARGUMENTS:
    !------------------------------------------------------------------------------
    implicit none
    integer            :: nu_nml                     ! unit for namelist file
    integer            :: nml_error                  ! namelist i/o error flag
    character(*), parameter    :: subName = "('init_soilwater_movement')"

    !-----------------------------------------------------------------------

! MUST agree with name in namelist and read statement
    namelist /soilwater_movement_inparm/      &
         soilwater_movement_method,    &
         upper_boundary_condition,     &
         lower_boundary_condition,     &
         dtmin,                        &
         verySmall,                    &
         xTolerUpper,                  &
         xTolerLower,                  &
         expensive,                    &
         inexpensive,                  &
         flux_calculation

    ! Default values for namelist

    soilwater_movement_method = zengdecker_2009
    upper_boundary_condition = bc_flux
    lower_boundary_condition = bc_aquifer

    dtmin=60._r8          
    verySmall=1.e-8_r8    
    xTolerUpper=1.e-1_r8  
    xTolerLower=1.e-2_r8  
    expensive=42
    inexpensive=1
    flux_calculation=inexpensive  

    ! Read soilwater_movement namelist
    if (masterproc) then
       nu_nml = getavu()
       open( nu_nml, file=trim(NLFilename), status='old', iostat=nml_error )
       call find_nlgroup_name(nu_nml, 'soilwater_movement_inparm', status=nml_error)
       if (nml_error == 0) then
          read(nu_nml, nml=soilwater_movement_inparm,iostat=nml_error)
          if (nml_error /= 0) then
             call endrun(subname // ':: ERROR reading soilwater_movement namelist')
          end if
       else
          call endrun(subname // ':: ERROR reading soilwater_movement namelist')
       end if
       close(nu_nml)
       call relavu( nu_nml )

!  test for namelist consistency
       if((soilwater_movement_method == zengdecker_2009) .and. &
            (lower_boundary_condition /= bc_aquifer)) then
          call endrun(subname // ':: ERROR inconsistent soilwater_movement namelist: ZD09 must use bc_aquifer lbc')
       endif
       if((use_bedrock) .and. (lower_boundary_condition /= bc_zero_flux)) then
          call endrun(subname // ':: ERROR inconsistent soilwater_movement namelist: use_bedrock requires bc_zero_flux lbc')
       endif
    endif

    call shr_mpi_bcast(soilwater_movement_method, mpicom)
    call shr_mpi_bcast(upper_boundary_condition, mpicom)
    call shr_mpi_bcast(lower_boundary_condition, mpicom)
    call shr_mpi_bcast(dtmin, mpicom)
    call shr_mpi_bcast(verySmall, mpicom)
    call shr_mpi_bcast(xTolerUpper, mpicom)
    call shr_mpi_bcast(xTolerLower, mpicom)
    call shr_mpi_bcast(expensive, mpicom)
    call shr_mpi_bcast(inexpensive, mpicom)
    call shr_mpi_bcast(flux_calculation, mpicom)


    if (masterproc) then

       write(iulog,*) ' '
       write(iulog,*) 'soilwater_movement settings:'
       write(iulog,*) '  soilwater_movement_method  = ',soilwater_movement_method
       write(iulog,*) '  upper_boundary_condition   = ',upper_boundary_condition
       write(iulog,*) '  lower_boundary_condition   = ',lower_boundary_condition

       write(iulog,*) '  use_bedrock                = ',use_bedrock
       write(iulog,*) '  dtmin                      = ',dtmin
       write(iulog,*) '  verySmall                  = ',verySmall
       write(iulog,*) '  xTolerUpper                = ',xTolerUpper
       write(iulog,*) '  xTolerLower                = ',xTolerLower
       write(iulog,*) '  expensive                  = ',expensive
       write(iulog,*) '  inexpensive                = ',inexpensive
       write(iulog,*) '  flux_calculation           = ',flux_calculation
    endif

  end subroutine init_soilwater_movement
  

!#2
   !------------------------------------------------------------------------------   
   function use_aquifer_layer() result(lres)
     !
     !DESCRIPTION
     ! return true if an aquifer layer is used 
     ! otherwise false
     implicit none
     logical :: lres

     if(lower_boundary_condition == bc_aquifer .or. lower_boundary_condition == bc_watertable)then
        lres=.true.
     else
        lres=.false.
     endif
     return

   end function use_aquifer_layer

!#3
  !-----------------------------------------------------------------------
  subroutine SoilWater(bounds, num_hydrologyc, filter_hydrologyc, &
       num_urbanc, filter_urbanc, soilhydrology_inst, soilstate_inst, &
       waterfluxbulk_inst, waterstatebulk_inst, temperature_inst, &
       canopystate_inst, energyflux_inst, soil_water_retention_curve)
    !
    ! DESCRIPTION
    ! select one subroutine to do the soil and root water coupling
    !
    !USES
    use shr_kind_mod      , only : r8 => shr_kind_r8
    use clm_varpar        , only : nlevsoi
    use decompMod         , only : bounds_type   
    use abortutils        , only : endrun   
    use clm_varpar        , only : nlevsoi
    use SoilHydrologyType , only : soilhydrology_type
    use SoilStateType     , only : soilstate_type
    use TemperatureType   , only : temperature_type
    use WaterFluxBulkType     , only : waterfluxbulk_type
    use EnergyFluxType    , only : energyflux_type
    use WaterStateBulkType    , only : waterstatebulk_type
    use CanopyStateType   , only : canopystate_type
    use ColumnType        , only : col
    use SoilWaterRetentionCurveMod, only : soil_water_retention_curve_type
    use clm_varcon        , only : denh2o, denice
    use clm_varctl,  only : use_flexibleCN   
    !
    ! !ARGUMENTS:
    type(bounds_type)        , intent(in)    :: bounds                ! bounds
    integer                  , intent(in)    :: num_hydrologyc        ! number of column soil points in column filter
    integer                  , intent(in)    :: filter_hydrologyc(:)  ! column filter for soil points
    integer                  , intent(in)    :: num_urbanc            ! number of column urban points in column filter
    integer                  , intent(in)    :: filter_urbanc(:)      ! column filter for urban points
    type(soilhydrology_type) , intent(inout) :: soilhydrology_inst
    type(soilstate_type)     , intent(inout) :: soilstate_inst
    type(waterfluxbulk_type)     , intent(inout) :: waterfluxbulk_inst
    type(energyflux_type)    , intent(in)    :: energyflux_inst
    type(waterstatebulk_type)    , intent(inout) :: waterstatebulk_inst
    type(canopystate_type)   , intent(inout) :: canopystate_inst
    type(temperature_type)   , intent(in)    :: temperature_inst
    class(soil_water_retention_curve_type), intent(in) :: soil_water_retention_curve
    !
    ! !LOCAL VARIABLES:
    character(len=32)              :: subname = 'SoilWater' ! subroutine name
    real(r8) :: xs(bounds%begc:bounds%endc)                !excess soil water above urban ponding limit
    real(r8) :: watmin
    integer  :: fc, c, j
    
    !------------------------------------------------------------------------------

    associate(                                                         &
      wa                 =>    soilhydrology_inst%wa_col             , & ! Input:  [real(r8) (:)   ] water in the unconfined aquifer (mm)
      dz                 =>    col%dz                                , & ! Input:  [real(r8) (:,:) ]  layer thickness (m)    
      h2osoi_ice         =>    waterstatebulk_inst%h2osoi_ice_col        , & ! Output: [real(r8) (:,:) ] liquid water (kg/m2)
      h2osoi_vol         =>    waterstatebulk_inst%h2osoi_vol_col        , & ! Output: [real(r8) (:,:) ] liquid water (kg/m2)
      h2osoi_liq         =>    waterstatebulk_inst%h2osoi_liq_col          & ! Output: [real(r8) (:,:) ] liquid water (kg/m2)
    )      

    select case(soilwater_movement_method)

    case (zengdecker_2009)

       call soilwater_zengdecker2009(bounds, num_hydrologyc, filter_hydrologyc, &
            num_urbanc, filter_urbanc, soilhydrology_inst, soilstate_inst, &
            waterfluxbulk_inst, waterstatebulk_inst, temperature_inst, &
            canopystate_inst, energyflux_inst, soil_water_retention_curve)

    case (moisture_form)

       call soilwater_moisture_form(bounds, num_hydrologyc, filter_hydrologyc, &
            num_urbanc, filter_urbanc, soilhydrology_inst, soilstate_inst, &
            waterfluxbulk_inst, waterstatebulk_inst, temperature_inst, &
            canopystate_inst, energyflux_inst, soil_water_retention_curve)

    case (mixed_form)

!!$       call soilwater_mixed_form(bounds, num_hydrologyc, filter_hydrologyc, &
!!$            num_urbanc, filter_urbanc, soilhydrology_inst, soilstate_inst, &
!!$            waterfluxbulk_inst, waterstate_inst, temperature_inst)

    case (head_form)

!!$       call soilwater_head_form(bounds, num_hydrologyc, filter_hydrologyc, &
!!$            num_urbanc, filter_urbanc, soilhydrology_inst, soilstate_inst, &
!!$            waterfluxbulk_inst, waterstate_inst, temperature_inst)

    case default

       call endrun(subname // ':: a SoilWater implementation must be specified!')          

    end select

    if (use_flexibleCN) then
       !a work around of the negative liquid water. Jinyun Tang, Jan 14, 2015
       watmin = 0.001_r8
       
       do j = 1, nlevsoi-1
          do fc = 1, num_hydrologyc
             c = filter_hydrologyc(fc)
             if (h2osoi_liq(c,j) < 0._r8) then
                xs(c) = watmin - h2osoi_liq(c,j)
             else
                xs(c) = 0._r8
             end if
             h2osoi_liq(c,j  ) = h2osoi_liq(c,j  ) + xs(c)
             h2osoi_liq(c,j+1) = h2osoi_liq(c,j+1) - xs(c)
          end do
       end do
       
       j = nlevsoi
       do fc = 1, num_hydrologyc
          c = filter_hydrologyc(fc)
          if (h2osoi_liq(c,j) < watmin) then
             xs(c) = watmin-h2osoi_liq(c,j)
          else
             xs(c) = 0._r8
          end if
          h2osoi_liq(c,j) = h2osoi_liq(c,j) + xs(c)
          wa(c) = wa(c) - xs(c)
       end do
       
       !update volumetric soil moisture for bgc calculation
       do j = 1, nlevsoi
          do fc = 1, num_hydrologyc
             c = filter_hydrologyc(fc)
             h2osoi_vol(c,j) = h2osoi_liq(c,j)/(dz(c,j)*denh2o) &
                  + h2osoi_ice(c,j)/(dz(c,j)*denice)
          enddo
       enddo
    end if
   end associate 
  end subroutine SoilWater

!#5
  !-----------------------------------------------------------------------   
  subroutine BaseflowSink(bounds, num_hydrologyc, &
       filter_hydrologyc, baseflow_sink, waterfluxbulk_inst, soilstate_inst)
    !
    ! Generic routine to apply baseflow as a sink condition that
    ! is vertically distributed over the soil column. 
    !
    !USES:
    use decompMod        , only : bounds_type
    use shr_kind_mod     , only : r8 => shr_kind_r8
    use clm_varpar       , only : nlevsoi, max_patch_per_col
    use SoilStateType    , only : soilstate_type
    use WaterFluxBulkType    , only : waterfluxbulk_type
    use PatchType        , only : patch
    use ColumnType       , only : col
    !
    ! !ARGUMENTS:
    type(bounds_type)    , intent(in)    :: bounds                          ! bounds
    integer              , intent(in)    :: num_hydrologyc                  ! number of column soil points in column filter
    integer              , intent(in)    :: filter_hydrologyc(:)            ! column filter for soil points
    real(r8)             , intent(out)   :: baseflow_sink(bounds%begc:,1:) ! vertically distributed baseflow sink (mm H2O/s) (+ = to rof)
    type(waterfluxbulk_type) , intent(inout) :: waterfluxbulk_inst
    type(soilstate_type) , intent(inout) :: soilstate_inst
    !
    ! !LOCAL VARIABLES:
    integer  :: p,c,fc,j                                              ! do loop indices
    integer  :: pi                                                    ! patch index
    real(r8) :: temp(bounds%begc:bounds%endc)                         ! accumulator for rootr weighting
    !-----------------------------------------------------------------------   

    ! Enforce expected array sizes
    SHR_ASSERT_ALL((ubound(baseflow_sink)  == (/bounds%endc, nlevsoi/)), errMsg(sourcefile, __LINE__))

!this is just a placeholder for now
    baseflow_sink = 0.

  end subroutine BaseflowSink


!#6
  !-----------------------------------------------------------------------
  subroutine soilwater_zengdecker2009(bounds, num_hydrologyc, filter_hydrologyc, &
       num_urbanc, filter_urbanc, soilhydrology_inst, soilstate_inst, &
       waterfluxbulk_inst, waterstatebulk_inst, temperature_inst, &
       canopystate_inst, energyflux_inst, soil_water_retention_curve)
    !
    ! !DESCRIPTION:
    ! Soil hydrology
    ! Soil moisture is predicted from a 10-layer model (as with soil
    ! temperature), in which the vertical soil moisture transport is governed
    ! by infiltration, runoff, gradient diffusion, gravity, and root
    ! extraction through canopy transpiration.  The net water applied to the
    ! surface layer is the snowmelt plus precipitation plus the throughfall
    ! of canopy dew minus surface runoff and evaporation.
    ! CLM3.5 uses a zero-flow bottom boundary condition.
    !
    ! The vertical water flow in an unsaturated porous media is described by
    ! Darcy's law, and the hydraulic conductivity and the soil negative
    ! potential vary with soil water content and soil texture based on the work
    ! of Clapp and Hornberger (1978) and Cosby et al. (1984). The equation is
    ! integrated over the layer thickness, in which the time rate of change in
    ! water mass must equal the net flow across the bounding interface, plus the
    ! rate of internal source or sink. The terms of water flow across the layer
    ! interfaces are linearly expanded by using first-order Taylor expansion.
    ! The equations result in a tridiagonal system equation.
    !
    ! Note: length units here are all millimeter
    ! (in temperature subroutine uses same soil layer
    ! structure required but lengths are m)
    !
    ! Richards equation:
    !
    ! d wat      d     d wat d psi
    ! ----- = - -- [ k(----- ----- - 1) ] + S
    !   dt      dz       dz  d wat
    !
    ! where: wat = volume of water per volume of soil (mm**3/mm**3)
    ! psi = soil matrix potential (mm)
    ! dt  = time step (s)
    ! z   = depth (mm)
    ! dz  = thickness (mm)
    ! qin = inflow at top (mm h2o /s)
    ! qout= outflow at bottom (mm h2o /s)
    ! s   = source/sink flux (mm h2o /s)
    ! k   = hydraulic conductivity (mm h2o /s)
    !
    !                       d qin                  d qin
    ! qin[n+1] = qin[n] +  --------  d wat(j-1) + --------- d wat(j)
    !                       d wat(j-1)             d wat(j)
    !                ==================|=================
    !                                  < qin
    !
    !                 d wat(j)/dt * dz = qin[n+1] - qout[n+1] + S(j)
    !
    !                                  > qout
    !                ==================|=================
    !                        d qout               d qout
    ! qout[n+1] = qout[n] + --------- d wat(j) + --------- d wat(j+1)
    !                        d wat(j)             d wat(j+1)
    !
    !
    ! Solution: linearize k and psi about d wat and use tridiagonal
    ! system of equations to solve for d wat,
    ! where for layer j
    !
    !
    ! r_j = a_j [d wat_j-1] + b_j [d wat_j] + c_j [d wat_j+1]
    !
    ! !USES:
    use shr_kind_mod               , only : r8 => shr_kind_r8     
    use shr_const_mod              , only : SHR_CONST_TKFRZ, SHR_CONST_LATICE, SHR_CONST_G
    use decompMod                  , only : bounds_type        
    use clm_varcon                 , only : wimp,grav,hfus,tfrz
    use clm_varcon                 , only : e_ice,denh2o, denice
    use clm_varpar                 , only : nlevsoi, max_patch_per_col, nlevgrnd
    use clm_time_manager           , only : get_step_size, get_nstep
    use column_varcon              , only : icol_roof, icol_road_imperv
    use clm_varctl                 , only : use_flexibleCN, use_hydrstress
    use TridiagonalMod             , only : Tridiagonal
    use abortutils                 , only : endrun     
    use SoilStateType              , only : soilstate_type
    use SoilHydrologyType          , only : soilhydrology_type
    use TemperatureType            , only : temperature_type
    use WaterFluxBulkType              , only : waterfluxbulk_type
    use EnergyFluxType             , only : energyflux_type
    use WaterStateBulkType             , only : waterstatebulk_type
    use CanopyStateType            , only : canopystate_type
    use SoilWaterRetentionCurveMod , only : soil_water_retention_curve_type
    use PatchType                  , only : patch
    use ColumnType                 , only : col
    use clm_varctl                 , only : iulog
    use SoilWaterPlantSinkMod      , only : COmpute_EffecRootFrac_And_VertTranSink
    !
    ! !ARGUMENTS:
    type(bounds_type)       , intent(in)    :: bounds               ! bounds
    integer                 , intent(in)    :: num_hydrologyc       ! number of column soil points in column filter
    integer                 , intent(in)    :: filter_hydrologyc(:) ! column filter for soil points
    integer                 , intent(in)    :: num_urbanc           ! number of column urban points in column filter
    integer                 , intent(in)    :: filter_urbanc(:)     ! column filter for urban points
    type(soilhydrology_type), intent(inout) :: soilhydrology_inst
    type(soilstate_type)    , intent(inout) :: soilstate_inst
    type(waterfluxbulk_type)    , intent(inout) :: waterfluxbulk_inst
    type(waterstatebulk_type)   , intent(inout) :: waterstatebulk_inst
    type(canopystate_type)  , intent(inout) :: canopystate_inst
    type(temperature_type)  , intent(in)    :: temperature_inst
    type(energyflux_type)   , intent(in)    :: energyflux_inst
    class(soil_water_retention_curve_type), intent(in) :: soil_water_retention_curve
    !
    ! !LOCAL VARIABLES:
    character(len=32) :: subname = 'soilwater_zengdecker2009' ! subroutine name 
    integer  :: p,c,fc,j                                     ! do loop indices
    integer  :: jtop(bounds%begc:bounds%endc)                ! top level at each column
    real(r8) :: dtime                                        ! land model time step (sec)
    real(r8) :: hk(bounds%begc:bounds%endc,1:nlevsoi)        ! hydraulic conductivity [mm h2o/s]
    real(r8) :: dhkdw(bounds%begc:bounds%endc,1:nlevsoi)     ! d(hk)/d(vol_liq)
    real(r8) :: amx(bounds%begc:bounds%endc,1:nlevsoi+1)     ! "a" left off diagonal of tridiagonal matrix
    real(r8) :: bmx(bounds%begc:bounds%endc,1:nlevsoi+1)     ! "b" diagonal column for tridiagonal matrix
    real(r8) :: cmx(bounds%begc:bounds%endc,1:nlevsoi+1)     ! "c" right off diagonal tridiagonal matrix
    real(r8) :: rmx(bounds%begc:bounds%endc,1:nlevsoi+1)     ! "r" forcing term of tridiagonal matrix
    real(r8) :: zmm(bounds%begc:bounds%endc,1:nlevsoi+1)     ! layer depth [mm]
    real(r8) :: dzmm(bounds%begc:bounds%endc,1:nlevsoi+1)    ! layer thickness [mm]
    real(r8) :: den                                          ! used in calculating qin, qout
    real(r8) :: dqidw0(bounds%begc:bounds%endc,1:nlevsoi+1)  ! d(qin)/d(vol_liq(i-1))
    real(r8) :: dqidw1(bounds%begc:bounds%endc,1:nlevsoi+1)  ! d(qin)/d(vol_liq(i))
    real(r8) :: dqodw1(bounds%begc:bounds%endc,1:nlevsoi+1)  ! d(qout)/d(vol_liq(i))
    real(r8) :: dqodw2(bounds%begc:bounds%endc,1:nlevsoi+1)  ! d(qout)/d(vol_liq(i+1))
    real(r8) :: dsmpdw(bounds%begc:bounds%endc,1:nlevsoi+1)  ! d(smp)/d(vol_liq)
    real(r8) :: num                                          ! used in calculating qin, qout
    real(r8) :: qin(bounds%begc:bounds%endc,1:nlevsoi+1)     ! flux of water into soil layer [mm h2o/s]
    real(r8) :: qout(bounds%begc:bounds%endc,1:nlevsoi+1)    ! flux of water out of soil layer [mm h2o/s]
    real(r8) :: s_node                                       ! soil wetness
    real(r8) :: s1                                           ! "s" at interface of layer
    real(r8) :: s2                                           ! k*s**(2b+2)
    real(r8) :: smp(bounds%begc:bounds%endc,1:nlevsoi)       ! soil matrix potential [mm]
    real(r8) :: sdamp                                        ! extrapolates soiwat dependence of evaporation
    integer  :: pi                                           ! patch index
    real(r8) :: temp(bounds%begc:bounds%endc)                ! accumulator for rootr weighting
    integer  :: jwt(bounds%begc:bounds%endc)                 ! index of the soil layer right above the water table (-)
    real(r8) :: smp1,dsmpdw1,wh,wh_zwt,ka
    real(r8) :: dwat2(bounds%begc:bounds%endc,1:nlevsoi+1)
    real(r8) :: dzq                                          ! used in calculating qin, qout (difference in equilbirium matric potential)
    real(r8) :: zimm(bounds%begc:bounds%endc,0:nlevsoi)      ! layer interface depth [mm]
    real(r8) :: zq(bounds%begc:bounds%endc,1:nlevsoi+1)      ! equilibrium matric potential for each layer [mm]
    real(r8) :: vol_eq(bounds%begc:bounds%endc,1:nlevsoi+1)  ! equilibrium volumetric water content
    real(r8) :: tempi                                        ! temp variable for calculating vol_eq
    real(r8) :: temp0                                        ! temp variable for calculating vol_eq
    real(r8) :: voleq1                                       ! temp variable for calculating vol_eq
    real(r8) :: zwtmm(bounds%begc:bounds%endc)               ! water table depth [mm]
    real(r8) :: imped(bounds%begc:bounds%endc,1:nlevsoi)             
    real(r8) :: vol_ice(bounds%begc:bounds%endc,1:nlevsoi)
    real(r8) :: z_mid
    real(r8) :: vwc_zwt(bounds%begc:bounds%endc)
    real(r8) :: vwc_liq(bounds%begc:bounds%endc,1:nlevsoi+1) ! liquid volumetric water content
    real(r8) :: smp_grad(bounds%begc:bounds%endc,1:nlevsoi+1)
    real(r8) :: dsmpds                                       !temporary variable
    real(r8) :: dhkds                                        !temporary variable
    real(r8) :: hktmp                                        !temporary variable
    integer :: nstep
    !-----------------------------------------------------------------------

    associate(& 
         z                 =>    col%z                              , & ! Input:  [real(r8) (:,:) ]  layer depth (m)                                 
         zi                =>    col%zi                             , & ! Input:  [real(r8) (:,:) ]  interface level below a "z" level (m)           
         dz                =>    col%dz                             , & ! Input:  [real(r8) (:,:) ]  layer thickness (m)                             

         origflag          =>    soilhydrology_inst%origflag        , & ! Input:  constant
         qcharge           =>    soilhydrology_inst%qcharge_col     , & ! Input:  [real(r8) (:)   ]  aquifer recharge rate (mm/s)                      
         zwt               =>    soilhydrology_inst%zwt_col         , & ! Input:  [real(r8) (:)   ]  water table depth (m)                             
         fracice           =>    soilhydrology_inst%fracice_col     , & ! Input:  [real(r8) (:,:) ]  fractional impermeability (-)                   
         icefrac           =>    soilhydrology_inst%icefrac_col     , & ! Input:  [real(r8) (:,:) ]  fraction of ice                                 
         hkdepth           =>    soilhydrology_inst%hkdepth_col     , & ! Input:  [real(r8) (:)   ]  decay factor (m)                                  

         smpmin            =>    soilstate_inst%smpmin_col          , & ! Input:  [real(r8) (:)   ]  restriction for min of soil potential (mm)        
         watsat            =>    soilstate_inst%watsat_col          , & ! Input:  [real(r8) (:,:) ]  volumetric soil water at saturation (porosity)  
         hksat             =>    soilstate_inst%hksat_col           , & ! Input:  [real(r8) (:,:) ]  hydraulic conductivity at saturation (mm H2O /s)
         bsw               =>    soilstate_inst%bsw_col             , & ! Input:  [real(r8) (:,:) ]  Clapp and Hornberger "b"                        
         sucsat            =>    soilstate_inst%sucsat_col          , & ! Input:  [real(r8) (:,:) ]  minimum soil suction (mm)                       
         eff_porosity      =>    soilstate_inst%eff_porosity_col    , & ! Input:  [real(r8) (:,:) ]  effective porosity = porosity - vol_ice         
         smp_l             =>    soilstate_inst%smp_l_col           , & ! Input:  [real(r8) (:,:) ]  soil matrix potential [mm]                      
         hk_l              =>    soilstate_inst%hk_l_col            , & ! Input:  [real(r8) (:,:) ]  hydraulic conductivity (mm/s)                   

         h2osoi_ice        =>    waterstatebulk_inst%h2osoi_ice_col     , & ! Input:  [real(r8) (:,:) ]  ice water (kg/m2)                               
         h2osoi_liq        =>    waterstatebulk_inst%h2osoi_liq_col     , & ! Input:  [real(r8) (:,:) ]  liquid water (kg/m2)                            
         h2osoi_vol        =>    waterstatebulk_inst%h2osoi_vol_col     , & ! Input:  [real(r8) (:,:) ]  volumetric soil water (0<=h2osoi_vol<=watsat) [m3/m3]

         qflx_deficit      =>    waterfluxbulk_inst%qflx_deficit_col    , & ! Input:  [real(r8) (:)   ]  water deficit to keep non-negative liquid water content
         qflx_infl         =>    waterfluxbulk_inst%qflx_infl_col       , & ! Input:  [real(r8) (:)   ]  infiltration (mm H2O /s)                          

         qflx_rootsoi_col  =>    waterfluxbulk_inst%qflx_rootsoi_col    , & ! Output: [real(r8) (:,:) ]  vegetation/soil water exchange (mm H2O/s) (+ = to atm)
         qflx_tran_veg_col =>    waterfluxbulk_inst%qflx_tran_veg_col   , & ! Input:  [real(r8) (:)   ]  vegetation transpiration (mm H2O/s) (+ = to atm)
         rootr_col         =>    soilstate_inst%rootr_col           , & ! Input:  [real(r8) (:,:) ]  effective fraction of roots in each soil layer  
         t_soisno          =>    temperature_inst%t_soisno_col        & ! Input:  [real(r8) (:,:) ]  soil temperature (Kelvin)                       
         )

      ! Get time step
      
      nstep = get_nstep()
      dtime = get_step_size()


      ! Because the depths in this routine are in mm, use local
      ! variable arrays instead of pointers

      do j = 1, nlevsoi
         do fc = 1, num_hydrologyc
            c = filter_hydrologyc(fc)
            zmm(c,j) = z(c,j)*1.e3_r8
            dzmm(c,j) = dz(c,j)*1.e3_r8
            zimm(c,j) = zi(c,j)*1.e3_r8

            ! calculate icefrac up here
            vol_ice(c,j) = min(watsat(c,j), h2osoi_ice(c,j)/(dz(c,j)*denice))
            icefrac(c,j) = min(1._r8,vol_ice(c,j)/watsat(c,j))
            vwc_liq(c,j) = max(h2osoi_liq(c,j),1.0e-6_r8)/(dz(c,j)*denh2o)
         end do
      end do

      do fc = 1, num_hydrologyc 
         c = filter_hydrologyc(fc)
         zimm(c,0) = 0.0_r8
         zwtmm(c)  = zwt(c)*1.e3_r8
      end do

      ! compute jwt index
      ! The layer index of the first unsaturated layer, i.e., the layer right above
      ! the water table

      do fc = 1, num_hydrologyc
         c = filter_hydrologyc(fc)
         jwt(c) = nlevsoi
         ! allow jwt to equal zero when zwt is in top layer
         do j = 1,nlevsoi
            if(zwt(c) <= zi(c,j)) then
               jwt(c) = j-1
               exit
            end if
         enddo

         ! compute vwc at water table depth (mainly for case when t < tfrz)
         !     this will only be used when zwt is below the soil column
         vwc_zwt(c) = watsat(c,nlevsoi)
         if(t_soisno(c,jwt(c)+1) < tfrz) then
            vwc_zwt(c) = vwc_liq(c,nlevsoi)
            do j = nlevsoi,nlevgrnd
               if(zwt(c) <= zi(c,j)) then
                  smp1 = hfus*(tfrz-t_soisno(c,j))/(grav*t_soisno(c,j)) * 1000._r8  !(mm)
                  !smp1 = max(0._r8,smp1)
                  smp1 = max(sucsat(c,nlevsoi),smp1)
                  vwc_zwt(c) = watsat(c,nlevsoi)*(smp1/sucsat(c,nlevsoi))**(-1._r8/bsw(c,nlevsoi))
                  ! for temperatures close to tfrz, limit vwc to total water content 
                  vwc_zwt(c) = min(vwc_zwt(c), 0.5*(watsat(c,nlevsoi) + h2osoi_vol(c,nlevsoi)) )
                  exit
               endif
            enddo
         endif
      end do

      ! calculate the equilibrium water content based on the water table depth

      do j=1,nlevsoi 
         do fc=1, num_hydrologyc
            c = filter_hydrologyc(fc)
            if ((zwtmm(c) <= zimm(c,j-1))) then 
               vol_eq(c,j) = watsat(c,j)

               ! use the weighted average from the saturated part (depth > wtd) and the equilibrium solution for the
               ! rest of the layer, the equilibrium solution is based on Clapp-Hornberg parameterization
               ! and no extension to full range swrc is needed

            else if ((zwtmm(c) .lt. zimm(c,j)) .and. (zwtmm(c) .gt. zimm(c,j-1))) then
               tempi = 1.0_r8
               temp0 = (((sucsat(c,j)+zwtmm(c)-zimm(c,j-1))/sucsat(c,j)))**(1._r8-1._r8/bsw(c,j))
               voleq1 = -sucsat(c,j)*watsat(c,j)/(1._r8-1._r8/bsw(c,j))/(zwtmm(c)-zimm(c,j-1))*(tempi-temp0)
               vol_eq(c,j) = (voleq1*(zwtmm(c)-zimm(c,j-1)) + watsat(c,j)*(zimm(c,j)-zwtmm(c)))/(zimm(c,j)-zimm(c,j-1))
               vol_eq(c,j) = min(watsat(c,j),vol_eq(c,j))
               vol_eq(c,j) = max(vol_eq(c,j),0.0_r8)
            else
               tempi = (((sucsat(c,j)+zwtmm(c)-zimm(c,j))/sucsat(c,j)))**(1._r8-1._r8/bsw(c,j))
               temp0 = (((sucsat(c,j)+zwtmm(c)-zimm(c,j-1))/sucsat(c,j)))**(1._r8-1._r8/bsw(c,j))
               vol_eq(c,j) = -sucsat(c,j)*watsat(c,j)/(1._r8-1._r8/bsw(c,j))/(zimm(c,j)-zimm(c,j-1))*(tempi-temp0)
               vol_eq(c,j) = max(vol_eq(c,j),0.0_r8)
               vol_eq(c,j) = min(watsat(c,j),vol_eq(c,j))
            endif
            zq(c,j) = -sucsat(c,j)*(max(vol_eq(c,j)/watsat(c,j),0.01_r8))**(-bsw(c,j))
            zq(c,j) = max(smpmin(c), zq(c,j))
         end do
      end do

      ! If water table is below soil column calculate zq for the 11th layer
      j = nlevsoi
      do fc=1, num_hydrologyc
         c = filter_hydrologyc(fc)
         if(jwt(c) == nlevsoi) then 
            tempi = 1._r8
            temp0 = (((sucsat(c,j)+zwtmm(c)-zimm(c,j))/sucsat(c,j)))**(1._r8-1._r8/bsw(c,j))
            vol_eq(c,j+1) = -sucsat(c,j)*watsat(c,j)/(1._r8-1._r8/bsw(c,j))/(zwtmm(c)-zimm(c,j))*(tempi-temp0)
            vol_eq(c,j+1) = max(vol_eq(c,j+1),0.0_r8)
            vol_eq(c,j+1) = min(watsat(c,j),vol_eq(c,j+1))
            zq(c,j+1) = -sucsat(c,j)*(max(vol_eq(c,j+1)/watsat(c,j),0.01_r8))**(-bsw(c,j))
            zq(c,j+1) = max(smpmin(c), zq(c,j+1))
         end if
      end do

      ! Hydraulic conductivity and soil matric potential and their derivatives

      sdamp = 0._r8
      do j = 1, nlevsoi
         do fc = 1, num_hydrologyc
            c = filter_hydrologyc(fc)
            ! compute hydraulic conductivity based on liquid water content only

            if (origflag == 1) then
               s1 = 0.5_r8*(h2osoi_vol(c,j) + h2osoi_vol(c,min(nlevsoi, j+1))) / &
                    (0.5_r8*(watsat(c,j)+watsat(c,min(nlevsoi, j+1))))
            else
               s1 = 0.5_r8*(vwc_liq(c,j) + vwc_liq(c,min(nlevsoi, j+1))) / &
                    (0.5_r8*(watsat(c,j)+watsat(c,min(nlevsoi, j+1))))
            endif
            s1 = min(1._r8, s1)
            s2 = hksat(c,j)*s1**(2._r8*bsw(c,j)+2._r8)

            ! replace fracice with impedance factor, as in zhao 97,99
            if (origflag == 1) then
               imped(c,j)=(1._r8-0.5_r8*(fracice(c,j)+fracice(c,min(nlevsoi, j+1))))
            else
               imped(c,j)=10._r8**(-e_ice*(0.5_r8*(icefrac(c,j)+icefrac(c,min(nlevsoi, j+1)))))
            endif
            hk(c,j) = imped(c,j)*s1*s2
            dhkdw(c,j) = imped(c,j)*(2._r8*bsw(c,j)+3._r8)*s2* &
                 (1._r8/(watsat(c,j)+watsat(c,min(nlevsoi, j+1))))

            !compute un-restricted hydraulic conductivity
            !call soil_water_retention_curve%soil_hk(hksat(c,j), imped(c,j), s1, bsw(c,j), hktmp, dhkds)
            !if(hktmp/=hk(c,j))write(10,*)'diff',hktmp,hk(c,j)
            !    call endrun('bad in hk')
            !endif    
            !apply ice impedance
            !hk(c,j) = imped(c,j)*hk(c,j)          
            !dhkdw(c,j) = imped(c,j) * dhkds * (1._r8/(watsat(c,j)+watsat(c,min(nlevsoi, j+1))))


            ! compute matric potential and derivative based on liquid water content only
            if (origflag == 1) then
               s_node = max(h2osoi_vol(c,j)/watsat(c,j), 0.01_r8)
            else
               s_node = max(vwc_liq(c,j)/watsat(c,j), 0.01_r8)
            endif
            s_node = min(1.0_r8, s_node)

            !call soil_water_retention_curve%soil_suction(sucsat(c,j), s_node, bsw(c,j), smp(c,j), dsmpds)

            smp(c,j) = -sucsat(c,j)*s_node**(-bsw(c,j))
            smp(c,j) = max(smpmin(c), smp(c,j))
            !do not turn on the line below, which will cause bit to bit error, jyt, 2014 Mar 6
            !dsmpdw(c,j) = dsmpds/watsat(c,j)

            if (origflag == 1) then             
               dsmpdw(c,j) = -bsw(c,j)*smp(c,j)/(s_node*watsat(c,j))
            else
               dsmpdw(c,j) = -bsw(c,j)*smp(c,j)/vwc_liq(c,j)
            endif

            smp_l(c,j) = smp(c,j)
            hk_l(c,j) = hk(c,j)

         end do
      end do

      ! aquifer (11th) layer
      do fc = 1, num_hydrologyc
         c = filter_hydrologyc(fc)
         zmm(c,nlevsoi+1) = 0.5*(1.e3_r8*zwt(c) + zmm(c,nlevsoi))
         if(jwt(c) < nlevsoi) then
            dzmm(c,nlevsoi+1) = dzmm(c,nlevsoi)
         else
            dzmm(c,nlevsoi+1) = (1.e3_r8*zwt(c) - zmm(c,nlevsoi))
         end if
      end do

      ! Set up r, a, b, and c vectors for tridiagonal solution

      ! Node j=1 (top)

      j = 1
      do fc = 1, num_hydrologyc
         c = filter_hydrologyc(fc)
         qin(c,j)    = qflx_infl(c)
         den    = (zmm(c,j+1)-zmm(c,j))
         dzq    = (zq(c,j+1)-zq(c,j))
         num    = (smp(c,j+1)-smp(c,j)) - dzq
         qout(c,j)   = -hk(c,j)*num/den
         dqodw1(c,j) = -(-hk(c,j)*dsmpdw(c,j)   + num*dhkdw(c,j))/den
         dqodw2(c,j) = -( hk(c,j)*dsmpdw(c,j+1) + num*dhkdw(c,j))/den
         rmx(c,j) =  qin(c,j) - qout(c,j) - qflx_rootsoi_col(c,j)
         amx(c,j) =  0._r8
         bmx(c,j) =  dzmm(c,j)*(sdamp+1._r8/dtime) + dqodw1(c,j)
         cmx(c,j) =  dqodw2(c,j)
         
      end do

      ! Nodes j=2 to j=nlevsoi-1

      do j = 2, nlevsoi - 1
         do fc = 1, num_hydrologyc
            c = filter_hydrologyc(fc)
            den    = (zmm(c,j) - zmm(c,j-1))
            dzq    = (zq(c,j)-zq(c,j-1))
            num    = (smp(c,j)-smp(c,j-1)) - dzq
            qin(c,j)    = -hk(c,j-1)*num/den
            dqidw0(c,j) = -(-hk(c,j-1)*dsmpdw(c,j-1) + num*dhkdw(c,j-1))/den
            dqidw1(c,j) = -( hk(c,j-1)*dsmpdw(c,j)   + num*dhkdw(c,j-1))/den
            den    = (zmm(c,j+1)-zmm(c,j))
            dzq    = (zq(c,j+1)-zq(c,j))
            num    = (smp(c,j+1)-smp(c,j)) - dzq
            qout(c,j)   = -hk(c,j)*num/den
            dqodw1(c,j) = -(-hk(c,j)*dsmpdw(c,j)   + num*dhkdw(c,j))/den
            dqodw2(c,j) = -( hk(c,j)*dsmpdw(c,j+1) + num*dhkdw(c,j))/den
            rmx(c,j)    =  qin(c,j) - qout(c,j) - qflx_rootsoi_col(c,j)
            amx(c,j)    = -dqidw0(c,j)
            bmx(c,j)    =  dzmm(c,j)/dtime - dqidw1(c,j) + dqodw1(c,j)
            cmx(c,j)    =  dqodw2(c,j)
            
         end do
      end do

      ! Node j=nlevsoi (bottom)

      j = nlevsoi
      do fc = 1, num_hydrologyc
         c = filter_hydrologyc(fc)
         if(j > jwt(c)) then !water table is in soil column
            den    = (zmm(c,j) - zmm(c,j-1))
            dzq    = (zq(c,j)-zq(c,j-1))
            num    = (smp(c,j)-smp(c,j-1)) - dzq
            qin(c,j)    = -hk(c,j-1)*num/den
            dqidw0(c,j) = -(-hk(c,j-1)*dsmpdw(c,j-1) + num*dhkdw(c,j-1))/den
            dqidw1(c,j) = -( hk(c,j-1)*dsmpdw(c,j)   + num*dhkdw(c,j-1))/den
            qout(c,j)   =  0._r8
            dqodw1(c,j) =  0._r8
            rmx(c,j)    =  qin(c,j) - qout(c,j) - qflx_rootsoi_col(c,j)
            amx(c,j)    = -dqidw0(c,j)
            bmx(c,j)    =  dzmm(c,j)/dtime - dqidw1(c,j) + dqodw1(c,j)
            cmx(c,j)    =  0._r8

            ! next set up aquifer layer; hydrologically inactive
            rmx(c,j+1) = 0._r8
            amx(c,j+1) = 0._r8
            bmx(c,j+1) = dzmm(c,j+1)/dtime
            cmx(c,j+1) = 0._r8
         else ! water table is below soil column

            ! compute aquifer soil moisture as average of layer 10 and saturation
            if(origflag == 1) then
               s_node = max(0.5*(1.0_r8+h2osoi_vol(c,j)/watsat(c,j)), 0.01_r8)
            else
               s_node = max(0.5*((vwc_zwt(c)+vwc_liq(c,j))/watsat(c,j)), 0.01_r8)
            endif
            s_node = min(1.0_r8, s_node)

            ! compute smp for aquifer layer
            !call soil_water_retention_curve%soil_suction(sucsat(c,j), s_node, bsw(c,j), smp1, dsmpds)
            smp1 = -sucsat(c,j)*s_node**(-bsw(c,j))
            smp1 = max(smpmin(c), smp1)

            ! compute dsmpdw for aquifer layer
            !dsmpdw1 = dsmpds/watsat(c,j)
            dsmpdw1 = -bsw(c,j)*smp1/(s_node*watsat(c,j))

            ! first set up bottom layer of soil column
            den    = (zmm(c,j) - zmm(c,j-1))
            dzq    = (zq(c,j)-zq(c,j-1))
            num    = (smp(c,j)-smp(c,j-1)) - dzq
            qin(c,j)    = -hk(c,j-1)*num/den
            dqidw0(c,j) = -(-hk(c,j-1)*dsmpdw(c,j-1) + num*dhkdw(c,j-1))/den
            dqidw1(c,j) = -( hk(c,j-1)*dsmpdw(c,j)   + num*dhkdw(c,j-1))/den
            den    = (zmm(c,j+1)-zmm(c,j))
            dzq    = (zq(c,j+1)-zq(c,j))
            num    = (smp1-smp(c,j)) - dzq
            qout(c,j)   = -hk(c,j)*num/den
            dqodw1(c,j) = -(-hk(c,j)*dsmpdw(c,j)   + num*dhkdw(c,j))/den
            dqodw2(c,j) = -( hk(c,j)*dsmpdw1 + num*dhkdw(c,j))/den

            rmx(c,j) =  qin(c,j) - qout(c,j) - qflx_rootsoi_col(c,j)
            amx(c,j) = -dqidw0(c,j)
            bmx(c,j) =  dzmm(c,j)/dtime - dqidw1(c,j) + dqodw1(c,j)
            cmx(c,j) =  dqodw2(c,j)

            ! next set up aquifer layer; den/num unchanged, qin=qout
            qin(c,j+1)    = qout(c,j)
            dqidw0(c,j+1) = -(-hk(c,j)*dsmpdw(c,j) + num*dhkdw(c,j))/den
            dqidw1(c,j+1) = -( hk(c,j)*dsmpdw1   + num*dhkdw(c,j))/den
            qout(c,j+1)   =  0._r8  ! zero-flow bottom boundary condition
            dqodw1(c,j+1) =  0._r8  ! zero-flow bottom boundary condition
            rmx(c,j+1) =  qin(c,j+1) - qout(c,j+1)
            amx(c,j+1) = -dqidw0(c,j+1)
            bmx(c,j+1) =  dzmm(c,j+1)/dtime - dqidw1(c,j+1) + dqodw1(c,j+1)
            cmx(c,j+1) =  0._r8
         endif
      end do

      ! Solve for dwat

      jtop(bounds%begc : bounds%endc) = 1
      call Tridiagonal(bounds, 1, nlevsoi+1, &
           jtop(bounds%begc:bounds%endc), &
           num_hydrologyc, filter_hydrologyc, &
           amx(bounds%begc:bounds%endc, :), &
           bmx(bounds%begc:bounds%endc, :), &
           cmx(bounds%begc:bounds%endc, :), &
           rmx(bounds%begc:bounds%endc, :), &
           dwat2(bounds%begc:bounds%endc, :) )

      ! Renew the mass of liquid water
      ! also compute qcharge from dwat in aquifer layer
      ! update in drainage for case jwt < nlevsoi

      do fc = 1,num_hydrologyc
         c = filter_hydrologyc(fc)

         do j = 1, nlevsoi
            h2osoi_liq(c,j) = h2osoi_liq(c,j) + dwat2(c,j)*dzmm(c,j)
         end do

         ! calculate qcharge for case jwt < nlevsoi
         if(jwt(c) < nlevsoi) then
            wh_zwt = 0._r8   !since wh_zwt = -sucsat - zq_zwt, where zq_zwt = -sucsat

            ! Recharge rate qcharge to groundwater (positive to aquifer)
            s_node = max(h2osoi_vol(c,jwt(c)+1)/watsat(c,jwt(c)+1), 0.01_r8)
            s1 = min(1._r8, s_node)

            !scs: this is the expression for unsaturated hk
            ka = imped(c,jwt(c)+1)*hksat(c,jwt(c)+1) &
                 *s1**(2._r8*bsw(c,jwt(c)+1)+3._r8)

            !compute unsaturated hk, this shall be tested later, because it
            !is not bit for bit
            !call soil_water_retention_curve%soil_hk(hksat(c,jwt(c)+1), s1, bsw(c,jwt(c)+1), ka)
            !apply ice impedance
            !ka = imped(c,jwt(c)+1) * ka 
            ! Recharge rate qcharge to groundwater (positive to aquifer)
            smp1 = max(smpmin(c), smp(c,max(1,jwt(c))))
            wh      = smp1 - zq(c,max(1,jwt(c)))

            !scs: original formulation
            if(jwt(c) == 0) then
               qcharge(c) = -ka * (wh_zwt-wh)  /((zwt(c)+1.e-3)*1000._r8)
            else
               !             qcharge(c) = -ka * (wh_zwt-wh)/((zwt(c)-z(c,jwt(c)))*1000._r8)
               !scs: 1/2, assuming flux is at zwt interface, saturation deeper than zwt
               qcharge(c) = -ka * (wh_zwt-wh)/((zwt(c)-z(c,jwt(c)))*1000._r8*2.0)
            endif

            ! To limit qcharge  (for the first several timesteps)
            qcharge(c) = max(-10.0_r8/dtime,qcharge(c))
            qcharge(c) = min( 10.0_r8/dtime,qcharge(c))
         else
            ! if water table is below soil column, compute qcharge from dwat2(11)
            qcharge(c) = dwat2(c,nlevsoi+1)*dzmm(c,nlevsoi+1)/dtime
         endif
      end do

      ! compute the water deficit and reset negative liquid water content
      !  Jinyun Tang
      do fc = 1, num_hydrologyc
         c = filter_hydrologyc(fc)
         qflx_deficit(c) = 0._r8
         do j = 1, nlevsoi
            if(h2osoi_liq(c,j)<0._r8)then
               qflx_deficit(c) = qflx_deficit(c) - h2osoi_liq(c,j)
            endif
         enddo
      enddo

    end associate 
         
  end subroutine soilwater_zengdecker2009

!#7
!-----------------------------------------------------------------------
   subroutine soilwater_moisture_form(bounds, num_hydrologyc, &
        filter_hydrologyc, num_urbanc, filter_urbanc, soilhydrology_inst, &
        soilstate_inst, waterfluxbulk_inst, waterstatebulk_inst, temperature_inst, &
        canopystate_inst, energyflux_inst, soil_water_retention_curve)
    !
    ! !DESCRIPTION:
    ! Soil hydrology
    ! Soil moisture is predicted from a n-layer model (as with soil
    ! temperature), in which the vertical soil moisture transport is governed
    ! by infiltration, runoff, gradient diffusion, gravity, and root
    ! extraction through canopy transpiration.  The net water applied to the
    ! surface layer is the snowmelt plus precipitation plus the throughfall
    ! of canopy dew minus surface runoff and evaporation.
    !
    ! Options are included for head conditions at the upper and lower boundary,
    ! a flux boundary condition at the upper boundary, and a free drainage
    ! condition
    ! at the lower boundary
    !
    ! The vertical water flow in an unsaturated porous media is described by
    ! Darcy's law, and the hydraulic conductivity and the soil negative
    ! potential vary with soil water content and soil texture based on the work
    ! of van Genuchten (1980). The equation is integrated over the layer
    ! thickness,
    ! in which the time rate of change in water mass must equal the net flow
    ! across
    ! the bounding interface, plus the rate of internal source or sink. The
    ! terms of
    ! water flow across the layer interfaces are linearly expanded by using
    ! first-order Taylor expansion. The equations result in a tridiagonal system
    ! equation.
    !
    ! Note: length units here are all millimeter
    ! (in temperature subroutine uses same soil layer
    ! structure required but lengths are m)
    !
    ! Richards equation:
    !
    ! d wat      d     d wat d psi
    ! ----- = - -- [ k(----- ----- - 1) ] + S
    !   dt      dz       dz  d wat
    !
    ! where: wat = volume of water per volume of soil (mm**3/mm**3)
    ! psi = soil matrix potential (mm)
    ! dt  = time step (s)
    ! z   = depth (mm)
    ! dz  = thickness (mm)
    ! qin = inflow at top (mm h2o /s)
    ! qout= outflow at bottom (mm h2o /s)
    ! s   = source/sink flux (mm h2o /s)
    ! k   = hydraulic conductivity (mm h2o /s)
    !
    !                       d qin                  d qin
    ! qin[n+1] = qin[n] +  --------  d wat(j-1) + --------- d wat(j)
    !                       d wat(j-1)             d wat(j)
    !                ==================|=================
    !                                  < qin
    !
    !                 d wat(j)/dt * dz = qin[n+1] - qout[n+1] + S(j)
    !
    !                                  > qout
    !                ==================|=================
    !                        d qout               d qout
    ! qout[n+1] = qout[n] + --------- d wat(j) + --------- d wat(j+1)
    !                        d wat(j)             d wat(j+1)
    !
    !
    ! Solution: linearize k and psi about d wat and use tridiagonal
    ! system of equations to solve for d wat,
    ! system of equations to solve for d wat,
    ! where for layer j
    !
    !
    ! r_j = a_j [d wat_j-1] + b_j [d wat_j] + c_j [d wat_j+1]
    !
    ! !USES:
    use shr_kind_mod         , only : r8 => shr_kind_r8
    use shr_const_mod        , only : SHR_CONST_TKFRZ, SHR_CONST_LATICE,SHR_CONST_G
    use abortutils           , only : endrun
    use decompMod            , only : bounds_type
    use clm_varctl           , only : iulog, use_hydrstress
    use clm_varcon           , only : denh2o, denice
    use clm_varpar           , only : nlevsoi
    use clm_time_manager     , only : get_step_size, get_nstep
    use SoilStateType        , only : soilstate_type
    use SoilHydrologyType    , only : soilhydrology_type
    use TemperatureType      , only : temperature_type
    use WaterFluxBulkType        , only : waterfluxbulk_type
    use WaterStateBulkType       , only : waterstatebulk_type
    use EnergyFluxType       , only : energyflux_type
    use CanopyStateType      , only : canopystate_type
    use SoilWaterRetentionCurveMod , only : soil_water_retention_curve_type
    use PatchType            , only : patch
    use ColumnType           , only : col
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type)       , intent(in)    :: bounds               ! bounds
    integer                 , intent(in)    :: num_hydrologyc       ! number of column soil points in column filter
    integer                 , intent(in)    :: filter_hydrologyc(:) ! column filter for soil points
    integer                 , intent(in)    :: num_urbanc           ! number of column urban points in column filter
    integer                 , intent(in)    :: filter_urbanc(:)     ! column filter for urban points
    type(soilhydrology_type), intent(inout) :: soilhydrology_inst
    type(soilstate_type)    , intent(inout) :: soilstate_inst
    type(waterfluxbulk_type)    , intent(inout) :: waterfluxbulk_inst
    type(waterstatebulk_type)   , intent(inout) :: waterstatebulk_inst
    type(temperature_type)  , intent(in)    :: temperature_inst
    type(canopystate_type)  , intent(inout) :: canopystate_inst
    type(energyflux_type)   , intent(in)    :: energyflux_inst
    class(soil_water_retention_curve_type), intent(in) :: soil_water_retention_curve

    ! !LOCAL VARIABLES:
    integer  :: nstep
    integer  :: nlayers
    integer  :: p,c,fc,j                                     ! do loop indices
    integer  :: jtop                                         ! top level at each column
    real(r8) :: dtime                                        ! land model time step (sec)
    real(r8) :: hk(bounds%begc:bounds%endc,1:nlevsoi)        ! hydraulic conductivity [mm h2o/s]
    real(r8) :: smp(bounds%begc:bounds%endc,1:nlevsoi)       ! soil matrix potential [mm]
    real(r8) :: dhkdw(bounds%begc:bounds%endc,1:nlevsoi)     ! d(hk)/d(vol_liq)
    real(r8) :: dsmpdw(bounds%begc:bounds%endc,1:nlevsoi)    ! d(smp)/d(vol_liq)
    real(r8) :: imped(bounds%begc:bounds%endc,1:nlevsoi)
    real(r8) :: rmx_save(bounds%begc:bounds%endc,1:nlevsoi)  ! residual vector
    real(r8) :: vwc_save(bounds%begc:bounds%endc,1:nlevsoi)  ! liquid volumetric water content
    real(r8) :: rhs(1:nlevsoi)                               ! RHS vector for input to lapack

    real(r8) :: vwc_liq(bounds%begc:bounds%endc,1:nlevsoi)   ! liquid volumetric water content
    real(r8) :: dt_dz(bounds%begc:bounds%endc,1:nlevsoi)     !temporary variable
    real(r8) :: qin(bounds%begc:bounds%endc,1:nlevsoi)       ! flux of water into soil layer [mm h2o/s]
    real(r8) :: qout(bounds%begc:bounds%endc,1:nlevsoi)      ! flux of water out of soil layer [mm h2o/s]
    real(r8) :: dqidw0(bounds%begc:bounds%endc,1:nlevsoi)    ! d(qin)/d(vol_liq(i-1))
    real(r8) :: dqidw1(bounds%begc:bounds%endc,1:nlevsoi)    ! d(qin)/d(vol_liq(i))
    real(r8) :: dqodw1(bounds%begc:bounds%endc,1:nlevsoi)    ! d(qout)/d(vol_liq(i))
    real(r8) :: dqodw2(bounds%begc:bounds%endc,1:nlevsoi)    ! d(qout)/d(vol_liq(i+1))
    real(r8) :: dwat(bounds%begc:bounds%endc,1:nlevsoi)
    real(r8) :: amx(bounds%begc:bounds%endc,1:nlevsoi)       ! "a" left off diagonal of tridiagonal matrix
    real(r8) :: bmx(bounds%begc:bounds%endc,1:nlevsoi)       ! "b" diagonal column for tridiagonal matrix
    real(r8) :: cmx(bounds%begc:bounds%endc,1:nlevsoi)       ! "c" right off diagonal tridiagonal matrix
    real(r8) :: rmx(bounds%begc:bounds%endc,1:nlevsoi)       ! "r" forcing term of tridiagonal matrix
    real(r8) :: dLow(1:nlevsoi-1)                            ! lower diagonal vector
    real(r8) :: dUpp(1:nlevsoi-1)                            ! upper diagonal vector
    real(r8) :: diag(1:nlevsoi)                              ! diagonal vector
    integer  :: err                                          ! error code from the lapack routines
    character(len=32) :: subname = 'soilwater_moisture_form' ! subroutine name   
    integer,parameter :: soil_hydraulic_properties_method = 0 ! 0 for CH78
    integer,parameter :: ice_impedance_method = 1            ! 1 for power form
    integer,parameter :: lapack=1
    integer,parameter :: clmOrig=2
!    integer,parameter :: tridiagSolution=clmOrig
    integer,parameter :: tridiagSolution=lapack

    integer  :: nSubstep  ! substep index
    real(r8) :: dtsub     ! length of the substep (seconds)
    real(r8) :: dtdone    ! substep completed (seconds)
!!$    real(r8),parameter :: dtmin=60._r8          ! minimum time step length (seconds)
!!$    real(r8),parameter :: verySmall=1.e-8_r8    ! a very small number: used to check for sub step completion
!!$    real(r8),parameter :: xTolerUpper=1.e-1_r8  ! tolerance to halve length of substep
!!$    real(r8),parameter :: xTolerLower=1.e-2_r8  ! tolerance to double length of substep
!!$    integer,parameter :: expensive=42
!!$    integer,parameter :: inexpensive=1
!!$    integer,parameter :: flux_calculation=inexpensive  ! used to calculate first order taylor series approximation of the net flux
    real(r8) :: qin_test,qout_test ! fluxes into and out of each layer
    real(r8) :: fluxNet0(1:nlevsoi)  ! net flux in each layer based on the first order taylor series approximation
    real(r8) :: fluxNet1(1:nlevsoi)  ! net flux in each layer based on the flux at the start of the substep
    real(r8) :: errorVec(1:nlevsoi)  ! error vector
    real(r8) :: errorMax(1)        ! max error
    real(r8) :: qcTemp             ! drainage flux for a given substep

! temporarily use local variables for the following
    real(r8) :: vwc_liq_ub(bounds%begc:bounds%endc)   ! liquid volumetric water content at upper boundary
    real(r8) :: vwc_liq_lb(bounds%begc:bounds%endc)   ! liquid volumetric water content at lower boundary
    real(r8) :: vLiqIter(bounds%begc:bounds%endc,1:nlevsoi)   !  iteration increment for the volumetric liquid water content (v/v)
    real(r8) :: vLiqRes(bounds%begc:bounds%endc,1:nlevsoi)   ! residual for the volumetric liquid water content (v/v)

    real(r8) :: dwat_temp
    !-----------------------------------------------------------------------

    associate(&
         nbedrock          =>    col%nbedrock                       , & ! Input:  [real(r8) (:,:) ]  depth to bedrock (m)                                 
         z                 =>    col%z                              , & ! Input:  [real(r8) (:,:) ]  layer depth (m)                                 
         zi                =>    col%zi                             , & ! Input:  [real(r8) (:,:) ]  interface level below a "z" level (m)           
         dz                =>    col%dz                             , & ! Input:  [real(r8) (:,:) ]  layer thickness (m)                             

         nsubsteps         =>    soilhydrology_inst%num_substeps_col, & ! Input:  [real(r8) (:)   ]  adaptive timestep counter

         qcharge           =>    soilhydrology_inst%qcharge_col     , & ! Input:  [real(r8) (:)   ]  aquifer recharge rate (mm/s)                      
         zwt               =>    soilhydrology_inst%zwt_col         , & ! Input:  [real(r8) (:)   ]  water table depth (m)                             

         smp_l             =>    soilstate_inst%smp_l_col           , & ! Input:  [real(r8) (:,:) ]  soil matrix potential [mm]                      
         hk_l              =>    soilstate_inst%hk_l_col            , & ! Input:  [real(r8) (:,:) ]  hydraulic conductivity (mm/s)                   
         h2osoi_ice        =>    waterstatebulk_inst%h2osoi_ice_col     , & ! Input:  [real(r8) (:,:) ]  ice water (kg/m2)                               
         h2osoi_liq        =>    waterstatebulk_inst%h2osoi_liq_col     , & ! Input:  [real(r8) (:,:) ]  liquid water (kg/m2)                            
         qflx_rootsoi_col  =>    waterfluxbulk_inst%qflx_rootsoi_col      &
         )  ! end associate statement

      ! Get time step

      nstep = get_nstep()
      dtime = get_step_size()

      ! main spatial loop
      do fc = 1, num_hydrologyc
         c = filter_hydrologyc(fc)

         ! set number of layers over which to solve soilwater movement
         nlayers = nbedrock(c)

         dwat_temp = 0.

         ! initialize the number of substeps
         nsubstep=0

         ! initialize substeps
         dtsub = dtime   ! length of the substep
         dtdone = 0._r8  ! substep completed

         ! initialize qcharge
         qcharge(c) = 0._r8

         ! adaptive sub-step loop
         do ! continuous do loop with exit clause

            ! keep track of the number of trials
            nsubstep=nsubstep+1

            ! calculate commonly used variables
            do j = 1, nlayers
               vwc_liq(c,j) = max(h2osoi_liq(c,j),1.0e-6_r8)/(dz(c,j)*denh2o)
               dt_dz(c,j)       = dtsub/(m_to_mm * dz(c,j))
            end do

            ! Hydraulic conductivity and soil matric potential and their derivatives
            call compute_hydraulic_properties(c, nlayers, &
                 soilhydrology_inst, soilstate_inst, &
                 soil_water_retention_curve, vwc_liq(c,1:nlayers), &
                 hk(c,1:nlayers) ,smp(c,1:nlayers), &
                 dhkdw(c,1:nlayers), dsmpdw(c,1:nlayers), &
                 imped(c,1:nlayers))

            ! Soil moisture fluxes and their derivatives
            call compute_moisture_fluxes_and_derivs(c, nlayers, &
                 soilhydrology_inst, soilstate_inst, &
                 temperature_inst, waterfluxbulk_inst, &
                 soil_water_retention_curve, &
                 vwc_liq(c,1:nlayers), &
                 hk(c,1:nlayers), &
                 smp(c,1:nlayers), &
                 dhkdw(c,1:nlayers), &
                 dsmpdw(c,1:nlayers), &
                 imped(c,1:nlayers), &
                 qin(c,1:nlayers), &
                 qout(c,1:nlayers), &
                 dqidw0(c,1:nlayers), &
                 dqidw1(c,1:nlayers), &
                 dqodw1(c,1:nlayers), &
                 dqodw2(c,1:nlayers))

            ! RHS of system of equations
            call compute_RHS_moisture_form(c, nlayers, &           
                 qflx_rootsoi_col(c,1:nlayers), &
                 vwc_liq(c,1:nlayers), &
                 qin(c,1:nlayers), &
                 qout(c,1:nlayers), &
                 dt_dz(c,1:nlayers), &
                 rmx(c,1:nlayers))
            
            ! LHS of system of equations
            call compute_LHS_moisture_form(c, nlayers, &
                 soilhydrology_inst, &
                 dt_dz(c,1:nlayers), &
                 dqidw0(c,1:nlayers), &
                 dqidw1(c,1:nlayers), &
                 dqodw1(c,1:nlayers), &
                 dqodw2(c,1:nlayers), &
                 amx(c,1:nlayers), &
                 bmx(c,1:nlayers), &
                 cmx(c,1:nlayers))

            ! solve A.X=B (done above, but using lapack as a test)

            ! solve tridiagonal matrix using the original clm solution
            if(tridiagSolution==clmOrig)then

               ! Solve for dwat

               jtop = 1
               call TridiagonalCol(c, 1, nlayers, &
                    jtop, &
                    !num_hydrologyc, filter_hydrologyc, &
                    amx(c, 1:nlayers), &
                    bmx(c, 1:nlayers), &
                    cmx(c, 1:nlayers), &
                    rmx(c, 1:nlayers), &
                    dwat(c, 1:nlayers) )

            ! solve tridiagonal matrix using the lapack solution
            else

               ! NOTE: check for efficiency gains by passing vectors directly (no
               ! copy); however, intent(inout)

               ! get a copy of the input vectors
               dLow(1:nlayers-1) = amx(filter_hydrologyc(fc),2:nlayers)
               diag(1:nlayers)   = bmx(filter_hydrologyc(fc),1:nlayers)
               dUpp(1:nlayers-1) = cmx(filter_hydrologyc(fc),1:nlayers-1)

               ! get a copy of the residual vector
               rhs(1:nlayers) = rmx(filter_hydrologyc(fc),1:nlayers)

               ! call the lapack tri-diagonal solver
               call dgtsv(nlayers,   & ! intent(in):    [integer]       number of state variables
                          1,         & ! intent(in):    [integer]       number of columns of the matrix B
                          dLow,      & ! intent(inout): [r8(nlayers-1)] sub-diagonal elements of A
                          diag,      & ! intent(inout): [r8(nlayers  )] diagonal elements of A
                          dUpp,      & ! intent(inout): [r8(nlayers-1)] super-diagonal elements of A
                          rhs,       & ! intent(inout): [r8(nlayers  )] RHS vector; becomes the solution vector on output
                          nlayers,   & ! intent(in):    [integer]       the leading dimension of matrix rhs
                          err)
               if(err/=0) call endrun(subname // ':: problem with the lapack solver')

               ! save the iteration increment
               dwat(filter_hydrologyc(fc),1:nlayers) = rhs(1:nlayers)

            endif  ! solution method for the tridiagonal solution

            ! **********
            ! error estimation segment
            ! (could be a separate subroutine)

            fluxNet1(1:nlevsoi) = 0._r8
            fluxNet0(1:nlevsoi) = 0._r8

            ! estimate errors for each layer
            do j = 1, nLayers

               ! compute the net flux using the first order taylor series approximation

               ! different options just for clarity
               ! NOTE: we would never use the expensive option; just helps clarify what is happenning
               if(flux_calculation==expensive)then

                  ! compute influx (mm s-1)
                  if(j==1)      then; qin_test = qin(c,j) + dqidw1(c,j)*dwat(c,j)  ! upper layer
                                else; qin_test = qin(c,j) + dqidw0(c,j)*dwat(c,j-1) + dqidw1(c,j)*dwat(c,j)
                  endif

                  ! compute outflux (mm s-1)
                  if(j==nLayers)then; qout_test = qout(c,j) + dqodw1(c,j)*dwat(c,j) ! lower layer
                                else; qout_test = qout(c,j) + dqodw1(c,j)*dwat(c,j) + dqodw2(c,j)*dwat(c,j+1)
                  endif

                  ! compute the net flux
                  fluxNet0(j) = qin_test - qout_test - qflx_rootsoi_col(c,j) 

               ! flux calculation is inexpensive
               else

                  ! convert iteration increment to a net flux (mm s-1)
                  fluxNet0(j) = dwat(c,j)/dt_dz(c,j)

               endif  ! switch between the expensive and inexpensive fluxcalculations

               ! compute the net flux at the start of the sub-step
               fluxNet1(j) = qin(c,j) - qout(c,j) - qflx_rootsoi_col(c,j)

            end do  ! looping through layers

            ! compute absolute errors
            errorVec(:nlayers) = abs(fluxNet1(:nlayers) - fluxNet0(:nlayers))*dtsub*0.5_r8
            errorMax = maxval(errorVec(:nlayers))

            ! check if the error is above the upper tolerance
            if(errorMax(1) > xTolerUpper .and. dtsub > dtmin)then
               dtsub = max(dtsub/2._r8, dtmin)  ! halve the length of the sub-step
               cycle                            ! sub-step rejected; try again
            endif


            ! end of error estimation segment
            ! **********

            ! Renew the mass of liquid water
            do j = 1, nlayers
               h2osoi_liq(c,j) = h2osoi_liq(c,j) + dwat(c,j) * (m_to_mm * dz(c,j))
            end do

             ! compute drainage from the bottom of the soil column
            select case(lower_boundary_condition)

               ! flux boundary condition
               case(bc_flux)
                  qcTemp = hk(c,nlayers) + dhkdw(c,nlayers)*dwat(c,nlayers)

               ! zero flux lower boundary
               case(bc_zero_flux)
                  qcTemp = 0._r8

               ! coupling with a water table
               case(bc_waterTable)
                  qcTemp = qout(c,nlayers) + dqodw1(c,nlayers)*dwat(c,nlayers)

               case default
                  call endrun(subname // ':: the lower boundary condition must be specified!')

            end select  ! case for the lower boundary condition

            ! increment the qcharge flux
            qcharge(c) = qcharge(c) + qcTemp *(dtsub / dtime)

            ! increment substep and check for completion of the substep
            dtdone = dtdone + dtsub
            if(abs(dtime - dtdone) < verySmall) exit  ! time step completed so exit the do loop

            ! check if the error is below the lower tolerance
            if(errorMax(1) < xTolerLower)then
                dtsub = dtsub*2._r8
            endif

            ! ensure that the substep does not exceed the time remaining
            dtsub = min(dtsub, dtime - dtdone)

         end do  ! substep loop

!  save number of adaptive substeps used during time step
         nsubsteps(c) = nsubstep

! check for negative moisture values
         do j = 2, nlayers
            if(h2osoi_liq(c,j) < -1e-6_r8) then
               write(*,*) 'layer, h2osoi_liq: ', c,j,h2osoi_liq(c,j)
               !      call endrun(subname // ':: negative soil moisture values found!')
            endif
         end do

      end do  ! spatial loop


! calculate qcharge when water table is in soil column and bc_watertable
      if(use_aquifer_layer()) then
      ! compute flux of water to aquifer

         call compute_qcharge(bounds, &
              num_hydrologyc, filter_hydrologyc, soilhydrology_inst, &
              soilstate_inst, waterstatebulk_inst, &
              soil_water_retention_curve, &
              dwat(bounds%begc:bounds%endc,1:nlevsoi), &
              smp(bounds%begc:bounds%endc,1:nlevsoi), &
              imped(bounds%begc:bounds%endc,1:nlevsoi), &
              vwc_liq(bounds%begc:bounds%endc,1:nlevsoi))

      endif

    end associate

   end subroutine soilwater_moisture_form

!#8
!-----------------------------------------------------------------------
   subroutine compute_hydraulic_properties(c, nlayers, &
        soilhydrology_inst, soilstate_inst, &
        soil_water_retention_curve, vwc_liq, hk ,smp, &
        dhkdw, dsmpdw, imped)
    !
    ! !DESCRIPTION:

    ! Calculate soil hydraulic conductivity and matric potential
    !
    ! !USES:
    use shr_kind_mod         , only : r8 => shr_kind_r8
    use shr_const_mod        , only : SHR_CONST_TKFRZ, SHR_CONST_LATICE, SHR_CONST_G
    use abortutils           , only : endrun
    use decompMod            , only : bounds_type
    use clm_varcon           , only : e_ice
    use clm_varpar           , only : nlevsoi
    use SoilWaterRetentionCurveMod, only : soil_water_retention_curve_type
    use SoilStateType        , only : soilstate_type
    use SoilHydrologyType    , only : soilhydrology_type
    use ColumnType           , only : col
    !
    ! !ARGUMENTS:
    implicit none
    integer                 , intent(in)    :: c                    ! column index
    integer                 , intent(in)    :: nlayers              ! lower boundary index

    type(soilhydrology_type), intent(in) :: soilhydrology_inst
    type(soilstate_type)    , intent(in) :: soilstate_inst

    class(soil_water_retention_curve_type), intent(in) :: soil_water_retention_curve
    real(r8), intent(in)  :: vwc_liq(1:nlayers)

    real(r8), intent(out) :: hk(1:nlayers)
    real(r8), intent(out) :: smp(1:nlayers)
    real(r8), intent(out) :: dhkdw(1:nlayers)
    real(r8), intent(out) :: dsmpdw(1:nlayers)
    real(r8), intent(out) :: imped(1:nlayers)
    !
    ! !LOCAL VARIABLES:
    integer  :: j                              ! do loop indices
    real(r8) :: s1                             ! "s" at interface of layer
 real(r8) :: s2(1:nlayers)                     ! "s" at layer node
    real(r8) :: dsmpds                         !temporary variable
    real(r8) :: dhkds                          !temporary variable
    character(len=32)  :: subname = 'calculate_hydraulic_properties'     ! subroutine name   
    !-----------------------------------------------------------------------

!scs: originally, associate statements selected sections rather than 
!     entire arrays, but due to pgi bug, removed array section selections
!     using array sections allowed consistent 1d indexing throughout
    associate(&
         icefrac           =>    soilhydrology_inst%icefrac_col     , & ! Input:  [real(r8) (:,:) ]  fraction of ice                                 
         watsat            =>    soilstate_inst%watsat_col          , & ! Input:  [real(r8) (:,:) ]  volumetric soil water at saturation (porosity)  
         smp_l             =>    soilstate_inst%smp_l_col           , & ! Input:  [real(r8) (:,:) ]  soil matrix potential [mm]                      
         hk_l              =>    soilstate_inst%hk_l_col              & ! Input:  [real(r8) (:,:) ]  hydraulic conductivity (mm/s)                   
         )  ! end associate statement

         hk     = 0._r8
         smp    = 0._r8
         dhkdw  = 0._r8
         dsmpdw = 0._r8
         imped  = 0._r8
      
         ! Hydraulic conductivity and soil matric potential and their derivatives

         ! compute the relative saturation at each layer first b/c
         ! it is used later in calculation of s1
         do j = 1, nlayers
            s2(j) = vwc_liq(j)/watsat(c,j)
            ! impose constraints on relative saturation at the layer node
            s2(j) = min(s2(j), 1._r8)
            s2(j) = max(0.01_r8, s2(j))
         enddo

         do j = 1, nlayers
            ! s1 is interface value, s2 is node value
            if(j==nlayers)then
             s1 = s2(j)
             call IceImpedance(icefrac(c,j), e_ice, imped(j) )
            else
             s1 = 0.5_r8 * (s2(j) + s2(j+1))
             call IceImpedance(0.5_r8*(icefrac(c,j) + icefrac(c,j+1)), e_ice, imped(j) )
            endif

  ! impose constraints on relative saturation at the layer interface
            s1 = min(s1, 1._r8)
            s1 = max(0.01_r8, s1)

            call soil_water_retention_curve%soil_hk(c, j, s1, &
                 imped(j), soilstate_inst, hk(j), dhkds=dhkds)

            call soil_water_retention_curve%soil_suction(c, j, s2(j), &
                 soilstate_inst,smp(j), dsmpds=dsmpds)


            ! NOTE: here just save derivative in hydraulic conductivity
            !        w.r.t relative raturation at the layer interface
            !          --> the derivatives at the nodes are calculated
            !               in compute_moisture_fluxes_and_derivs
            dhkdw(j) = dhkds  ! NOTE: the variable name does not make sense 
            
            ! compute derivative w.r.t. volumetric liquid water content
            ! NOTE: derivative for the layer
            dsmpdw(j) = dsmpds /  watsat(c,j) 
            
            ! save for the output
            smp_l(c,j) = smp(j)
            hk_l(c,j) = hk(j)

         end do

    end associate

   end subroutine compute_hydraulic_properties

!#9
!-----------------------------------------------------------------------
   subroutine compute_moisture_fluxes_and_derivs(c, nlayers, &
        soilhydrology_inst, soilstate_inst, temperature_inst, waterfluxbulk_inst, &
        soil_water_retention_curve, vwc_liq, hk ,smp, dhkdw, dsmpdw, &
        imped, qin, qout, dqidw0, dqidw1, dqodw1, dqodw2)
    !
    ! !DESCRIPTION:

    ! Calculate fluxes at the boundary of each layer and derivatives w.r.t.
    ! relevant state variables
    !
    ! !USES:
    use shr_kind_mod         , only : r8 => shr_kind_r8
    use shr_const_mod        , only : SHR_CONST_TKFRZ, SHR_CONST_LATICE, SHR_CONST_G
    use abortutils           , only : endrun
    use decompMod            , only : bounds_type
    use clm_varpar           , only : nlevsoi
    use SoilWaterRetentionCurveMod, only : soil_water_retention_curve_type
    use SoilStateType        , only : soilstate_type
    use SoilHydrologyType    , only : soilhydrology_type
    use TemperatureType      , only : temperature_type
    use WaterFluxBulkType        , only : waterfluxbulk_type
    use ColumnType           , only : col
    !
    ! !ARGUMENTS:
    implicit none
    integer                 , intent(in)    :: c                    ! column index
    integer                 , intent(in)    :: nlayers              ! lower boundary index

    type(soilhydrology_type), intent(in) :: soilhydrology_inst
    type(soilstate_type)    , intent(in) :: soilstate_inst
    type(temperature_type)  , intent(in) :: temperature_inst
    type(waterfluxbulk_type)    , intent(in) :: waterfluxbulk_inst

    class(soil_water_retention_curve_type), intent(in) :: soil_water_retention_curve
    real(r8), intent(in)  :: vwc_liq(1:nlayers)
    real(r8), intent(in)  :: hk(1:nlayers)
    real(r8), intent(in)  :: smp(1:nlayers)
    real(r8), intent(in)  :: dhkdw(1:nlayers)
    real(r8), intent(in)  :: dsmpdw(1:nlayers)
    real(r8), intent(in)  :: imped(1:nlayers)

    real(r8), intent(out) :: qin(1:nlayers)
    real(r8), intent(out) :: qout(1:nlayers)
    real(r8), intent(out) :: dqidw0(1:nlayers)
    real(r8), intent(out) :: dqidw1(1:nlayers)
    real(r8), intent(out) :: dqodw1(1:nlayers)
    real(r8), intent(out) :: dqodw2(1:nlayers)
    !
    ! !LOCAL VARIABLES:
    integer  :: j                                            ! do loop indices
    real(r8) :: s1                                           ! "s" at interface of layer
    real(r8) :: s2                                           ! "s" at layer node
    real(r8) :: smp1          !temporary variable
    real(r8) :: hk1           !temporary variable
    real(r8) :: num, den      ! used in calculating qin, qout
    real(r8) :: dhkds1, dhkds2                                        !temporary variable
    real(r8),parameter :: m_to_mm = 1.e3_r8  !convert meters to mm
!scs: temporarily use local variables for the following
    real(r8) :: vwc_liq_ub   ! liquid volumetric water content at upper boundary
    real(r8) :: vwc_liq_lb   ! liquid volumetric water content at lower boundary
    character(len=32)  :: subname = 'calculate_moisture_fluxes_and_derivs'     ! subroutine name   
    integer  :: jwt              !index of layer above water table
    real(r8) :: dhkdw1, dsmpdw1, dsmpds1
    integer,parameter :: zdflag = 0
    !-----------------------------------------------------------------------

    associate(&
         z                 =>    col%z                              , & ! Input:  [real(r8) (:,:) ]  layer depth (m)                                 
         zi                =>    col%zi                             , & ! Input:  [real(r8) (:,:) ]  interface level below a "z" level (m)           
         dz                =>    col%dz                            , & ! Input:  [real(r8) (:,:) ]  layer thickness (m)                             
         zwt               =>    soilhydrology_inst%zwt_col         , & ! Input:  [real(r8) (:)   ]  water table depth (m)                             
         t_soisno          =>    temperature_inst%t_soisno_col      , & ! Input:  [real(r8) (:,:) ]  soil temperature (Kelvin)                      
         qflx_infl         =>    waterfluxbulk_inst%qflx_infl_col       , & ! Input:  [real(r8) (:)   ]  infiltration (mm H2O /s)                          
         watsat            =>    soilstate_inst%watsat_col            & ! Input:  [real(r8) (:,:) ]  volumetric soil water at saturation (porosity)  
         )  ! end associate statement

      qin    = 0._r8
      qout   = 0._r8
      dqidw0 = 0._r8
      dqidw1 = 0._r8
      dqodw1 = 0._r8
      dqodw2 = 0._r8

    ! Soil moisture fluxes and their derivatives

    j = 1
       
    ! select type of boundary condition
    select case(upper_boundary_condition)
          
       ! head boundary conditions
       case(bc_head)
          
          ! compute the relative saturation at the upper boundary
          s1 = vwc_liq_ub / watsat(c,j)
          !               s1 = (vwc_liq_ub - watres(c,j)) / (watsat(c,j) - watres(c,j))
          s1 = min(s1, 1._r8)
          s1 = max(0.01_r8, s1)
          
          ! compute the hydraulic conductivity and matric potential at the boundary
            call soil_water_retention_curve%soil_hk(c, j, s1, &
                 imped(j), soilstate_inst, hk1)

            call soil_water_retention_curve%soil_suction(c, j, s1, &
                 soilstate_inst,smp1)
          
          ! compute the infiltration to the top of the layer
          den      = z(c,j)* m_to_mm * 0.5_r8
          qin(j) = -hk1*(smp(j) - smp1)/den + hk1
          
          ! compute the flux derivative w.r.t influx
          dqidw1(j) = -hk1*dsmpdw(j)/den
          
       ! flux boundary condition
       case(bc_flux)
          
          ! specify infiltration
          qin(j) = qflx_infl(c)
          
          ! compute the flux derivative w.r.t influx
          dqidw1(j) = 0._r8
          
       ! check the case is specified correctly
       case default
          
          call endrun(subname // ':: the upper boundary condition must be specified!')
          
    end select       

    ! NOTE: the rest of the code for j=1 is identical to the code for
    !         the interior nodes (j=2,...,nlevsoi-1)
    ! Not sure why this is split out in this way....

    ! compute derivatives in hydraulic conductivity at the interface w.r.t.
    !  volumetric liquid water content in the layer above and the layer  below
    ! NOTE: dhkdw(j) is the derivative in hydrualic conductivity at the
    !        layer interface w.r.t relative saturation at the interface
    dhkds1 = 0.5_r8 * dhkdw(j) / watsat(c,j)   ! derivative w.r.t. volumetric liquid water in the upper layer
    dhkds2 = 0.5_r8 * dhkdw(j) / watsat(c,j+1) ! derivative w.r.t. volumetric liquid water in the lower layer

!scs: this is how zd is done
    if (zdflag == 1) then 
       dhkds1 = dhkdw(j)/(watsat(c,j)+watsat(c,min(nlevsoi, j+1)))
       dhkds2 = dhkds1
    endif
!scs

    ! compute flux at the bottom of the j-th layer
    ! NOTE: hk(j) is hydraulic conductivity at the bottom of the j-th
    ! layer
    num    = (smp(j+1)-smp(j))
    den    = m_to_mm * (z(c,j+1)-z(c,j))
    qout(j)   = -hk(j)*num/den + hk(j)

    ! compute flux derivatives
    dqodw1(j) = (hk(j)*dsmpdw(j) - dhkds1*num)/den + dhkds1
    dqodw2(j) = (-hk(j)*dsmpdw(j+1) - dhkds2*num)/den + dhkds2

    ! Interior nodes

    do j = 2, nlayers - 1
          
       ! get the flux from above
       qin(j) = qout(j-1)
          
       ! get the derivatives from above
       dqidw0(j) = dqodw1(j-1)
       dqidw1(j) = dqodw2(j-1)
          
       ! compute derivatives in hydraulic conductivity at the interface w.r.t.
       !  volumetric liquid water content in the layer above and the layer below
       ! NOTE: dhkdw(j) is the derivative in hydrualic conductivity at the
       !        layer interface w.r.t relative saturation at the interface
       dhkds1 = 0.5_r8 * dhkdw(j) / watsat(c,j)   ! derivative w.r.t. volumetric liquid water in the upper layer
       dhkds2 = 0.5_r8 * dhkdw(j) / watsat(c,j+1) ! derivative w.r.t. volumetric liquid water in the lower layer
!scs: this is how zd is done
             if (zdflag == 1) then 
                dhkds1 = dhkdw(j)/(watsat(c,j)+watsat(c,min(nlevsoi, j+1)))
                dhkds2 = dhkds1
             endif
!scs
          
       ! compute flux at the bottom of the j-th layer
       ! NOTE: hk(j) is hydraulic conductivity at the bottom of the j-th  layer
       num    = (smp(j+1)-smp(j))
       den    = m_to_mm * (z(c,j+1)-z(c,j))
       qout(j)   = -hk(j)*num/den + hk(j)
          
       ! compute flux derivatives
       dqodw1(j) = (hk(j)*dsmpdw(j) - dhkds1*num)/den + dhkds1 
       dqodw2(j) = (-hk(j)*dsmpdw(j+1) - dhkds2*num)/den + dhkds2

    end do

    ! Node j=nlevsoi (bottom)
    
    j = nlayers
       
    ! get the flux from above
    qin(j) = qout(j-1)
       
    ! get the derivatives from above
    dqidw0(j) = dqodw1(j-1)
    dqidw1(j) = dqodw2(j-1)
       
    ! select type of boundary condition
    select case(lower_boundary_condition)

       ! water table boundary conditions (moving lower boundary)
       case(bc_waterTable)

          jwt = nlevsoi
          !  Locate index of layer above water table
          do j = 1,nlevsoi
             if(zwt(c) <= zi(c,j)) then
                jwt = j-1
                exit
             end if
          enddo

          !  water table within soil column
          j=nlevsoi
          if(j > jwt) then

            ! no drainage
            qout(j) = 0.
            dqodw1(j) = 0.
            dqodw2(j) = 0.

          ! water table below soil column
          else

             ! following Niu et al. (JGR 2007), assume hydraulic conductivity at the bottom of the soil column
             ! is equal to hydraulic conductivty at the lowest node. This provides a free-drainage lower boundary
             ! condition when the water table is a long way below the soil column
             dhkds1 = dhkdw(j) / watsat(c,j)

!scs: this is how zd is done
             if (zdflag == 1) then 
                dhkds1 = dhkdw(j)/(watsat(c,j)+watsat(c,min(nlevsoi, j+1)))
                dhkds2 = dhkds1
             endif
!scs
             ! compute flux
             num    = -smp(j)  ! NOTE: assume saturation at water table depth (smp=0)
             den    = m_to_mm * (zwt(c) - z(c,j))
             qout(j) = -hk(j)*num/den + hk(j)  ! using hydraulic conductivity at the midpoint of the two nodes

             ! compute derivatives
             dqodw1(j) = (hk(j)*dsmpdw(j) - dhkds1*num)/den + dhkds1
             dqodw2(j) = 0._r8  ! no layer below

          endif  ! switch between water table above/below soil column

          
       ! head boundary conditions
       case(bc_head)
          
          ! compute the relative saturation at the lower boundary
          s1 = vwc_liq_lb / watsat(c,j)
!scs: mc's original expression          s1 = (vwc_liq_lb - watres(c,j)) / (watsat(c,j) - watres(c,j))
          s1 = min(s1, 1._r8)
          s1 = max(0.01_r8, s1)
          
          ! compute the hydraulic conductivity and matric potential at the boundary
            call soil_water_retention_curve%soil_hk(c, j, s1, &
                 imped(j), soilstate_inst, hk1)

            call soil_water_retention_curve%soil_suction(c, j, s1, &
                 soilstate_inst,smp1)
       
          ! compute the flux at the bottom of the layer
          den      = m_to_mm * z(c,j) * 0.5_r8
          qout(j) = -hk1*(smp1 - smp(j))/den + hk1
          
          ! compute the flux derivative
          dqodw1(j) = hk1*dsmpdw(j)/den

          ! stop, because not tested yet
          call endrun(subname // ':: the specified lower boundary condition is not fully implemented/tested yet')
          
       ! flux boundary condition
       case(bc_flux)
          
          ! free drainage
          qout(j) = hk(j)
          
          ! compute the flux derivative
          dqodw1(j) = dhkdw(j) / watsat(c,j)

       case(bc_zero_flux)
    
          ! no drainage
          qout(j) = 0.
          
          ! compute the flux derivative
          dqodw1(j) = 0.

       ! check the case is specified correctly
       case default
          
       call endrun(subname // ':: the lower boundary condition must be specified!')
       
    end select
       
    
    end associate

   end subroutine compute_moisture_fluxes_and_derivs


!#10
!-----------------------------------------------------------------------
   subroutine compute_RHS_moisture_form(c, nlayers, vert_trans_sink, &
        vwc_liq, qin, qout, dt_dz, rmx)
    !
    ! !DESCRIPTION:

    ! Calculate RHS of moisture-based form of Richards equation
    !
    ! !USES:
    use shr_kind_mod         , only : r8 => shr_kind_r8
    use shr_const_mod        , only : SHR_CONST_TKFRZ, SHR_CONST_LATICE, SHR_CONST_G
    use abortutils           , only : endrun
    use decompMod            , only : bounds_type
    use clm_varpar           , only : nlevsoi
    use ColumnType           , only : col
    !
    ! !ARGUMENTS:
    implicit none
    integer                 , intent(in)    :: c                    ! column index
    integer                 , intent(in)    :: nlayers              ! lower boundary index

    real(r8), intent(in)  :: vert_trans_sink(1:nlayers) ! vertically distributed transpiration sink (mm H2O/s) (+ = to atm)
    real(r8), intent(in)  :: vwc_liq(1:nlayers)
    real(r8), intent(in)  :: qin(1:nlayers)
    real(r8), intent(in)  :: qout(1:nlayers)
    real(r8), intent(in)  :: dt_dz(1:nlayers)
    real(r8), intent(out) :: rmx(1:nlayers)       ! "r" forcing term of tridiagonal matrix
    !
    ! !LOCAL VARIABLES:
    integer  :: j                                     ! do loop indices
    real(r8) :: fluxNet       !temporary variable
    integer,parameter :: ixIterate=1
    integer,parameter :: ixNonIterate=2
    integer           :: iterateOption=ixNonIterate
    character(len=32)  :: subname = 'compute_RHS_moisture_form'     ! subroutine name   
    !-----------------------------------------------------------------------

    associate(&
         dz                 =>    col%dz                               & ! Input:  [real(r8) (:,:) ]  layer depth (m)                                 
         )  ! end associate statement

         rmx = 0._r8

!  Compute RHS of moisture-based form of Richards equation
         do j = 1, nlayers

            fluxNet     = qin(j) - qout(j) - vert_trans_sink(j)
!            rmx(j)    = vwc_liq(j) - (vwc_liq(j) + fluxNet*dt_dz(j))

 ! different options for iterative and non-iterative solutions
            select case(iterateOption)
             case(ixIterate)
              ! for iterative solution entire soil hydrology routine should be
              ! called multiple times using updated versions of vwc_trial
              ! vwc_trial = wwc_liq for iteration=1
              !rmx(j)    = vwc_trial(j) - (vwc_liq(j) + fluxNet*dt_dz(j))
              call endrun(subname // ':: ERROR iterative option not yet implemented')
             case(ixNonIterate)
                rmx(j)    = -fluxNet*dt_dz(j)
             case default
                call endrun(subname // ':: ERROR unknown iterate option')
            end select

         end do
   
    end associate

   end subroutine compute_RHS_moisture_form

!#11
!-----------------------------------------------------------------------
   subroutine compute_LHS_moisture_form(c, nlayers, &
        soilhydrology_inst, dt_dz, dqidw0, dqidw1, dqodw1, dqodw2, &
        amx, bmx, cmx)
    !
    ! !DESCRIPTION:

    ! Calculate LHS of moisture-based form of Richards equation
    ! LHS is a tridiagonal matrix 
    !
    ! !USES:
    use shr_kind_mod         , only : r8 => shr_kind_r8
    use shr_const_mod        , only : SHR_CONST_TKFRZ, SHR_CONST_LATICE, SHR_CONST_G
    use abortutils           , only : endrun
    use decompMod            , only : bounds_type
    use clm_varpar           , only : nlevsoi
    use ColumnType           , only : col
    use SoilHydrologyType    , only : soilhydrology_type
    !
    ! !ARGUMENTS:
    implicit none
    integer                 , intent(in)    :: c                    ! column index
    integer                 , intent(in)    :: nlayers              ! lower boundary index
    type(soilhydrology_type), intent(in)    :: soilhydrology_inst

    real(r8), intent(in)  :: dt_dz(1:nlayers)
    real(r8), intent(in)  :: dqidw0(1:nlayers)
    real(r8), intent(in)  :: dqidw1(1:nlayers)
    real(r8), intent(in)  :: dqodw1(1:nlayers)
    real(r8), intent(in)  :: dqodw2(1:nlayers)
    real(r8), intent(out) :: amx(1:nlayers)       ! "a" left off diagonal of tridiagonal matrix
    real(r8), intent(out) :: bmx(1:nlayers)       ! "b" diagonal column for tridiagonal matrix
    real(r8), intent(out) :: cmx(1:nlayers)       ! "c" right off diagonal tridiagonal matrix

    !
    ! !LOCAL VARIABLES:
    integer  :: j                                     ! do loop indices
    character(len=32)  :: subname = 'compute_LHS_moisture_form'     ! subroutine name   
    integer  :: jwt              !index of layer above water table
    !-----------------------------------------------------------------------

    associate(&
         zi                =>    col%zi                             , & ! Input:  [real(r8) (:,:) ]  layer interface depth (m)                                 
         zwt               =>    soilhydrology_inst%zwt_col           & ! Input:  [real(r8) (:)   ]  water table depth (m)                             
         )  ! end associate statement


      amx = 0._r8
      bmx = 0._r8
      cmx = 0._r8

!  Compute LHS of moisture-based form of Richards equation
!  Top soil layer
      j=1
      amx(j)    = 0._r8
      bmx(j)    = -1._r8 - (-dqidw1(j) + dqodw1(j))*dt_dz(j)
      cmx(j)    = -dqodw2(j)*dt_dz(j)         

!  Interior soil layers      
      do j = 2, nlayers-1
         amx(j)    = dqidw0(j)*dt_dz(j)
         bmx(j)    = -1._r8 - (-dqidw1(j) + dqodw1(j))*dt_dz(j)
         cmx(j)    = -dqodw2(j)*dt_dz(j)            
      enddo

         j=nlayers
         amx(j)    = dqidw0(j)*dt_dz(j)
         bmx(j)    = -1._r8 - (-dqidw1(j) + dqodw1(j))*dt_dz(j)
         cmx(j)    = 0._r8         

    end associate

   end subroutine compute_LHS_moisture_form

!#12
!-----------------------------------------------------------------------
   subroutine compute_qcharge(bounds, num_hydrologyc, &
        filter_hydrologyc, soilhydrology_inst, soilstate_inst, &
        waterstatebulk_inst, soil_water_retention_curve, &
        dwat, smp, imped, vwc_liq)
    !
    ! !DESCRIPTION:

    ! Calculate additional terms due to presence of aquifer layer
    !
    ! !USES:
    use shr_kind_mod         , only : r8 => shr_kind_r8
    use shr_const_mod        , only : SHR_CONST_TKFRZ, SHR_CONST_LATICE, SHR_CONST_G
    use abortutils           , only : endrun
    use decompMod            , only : bounds_type
    use clm_time_manager     , only : get_step_size
    use clm_varpar           , only : nlevsoi
    use SoilWaterRetentionCurveMod, only : soil_water_retention_curve_type
    use SoilStateType        , only : soilstate_type
    use SoilHydrologyType    , only : soilhydrology_type
    use TemperatureType      , only : temperature_type
    use WaterFluxBulkType        , only : waterfluxbulk_type
    use WaterStateBulkType       , only : waterstatebulk_type
    use ColumnType           , only : col
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type)       , intent(in)    :: bounds               ! bounds
    integer                 , intent(in)    :: num_hydrologyc       ! number of column soil points in column filter
    integer                 , intent(in)    :: filter_hydrologyc(:) ! column filter for soil points

    type(soilhydrology_type), intent(in) :: soilhydrology_inst
    type(soilstate_type)    , intent(in) :: soilstate_inst
    type(waterstatebulk_type)    , intent(in) :: waterstatebulk_inst

!    integer,  intent(in)  :: soil_hydraulic_properties_method
    class(soil_water_retention_curve_type), intent(in) :: soil_water_retention_curve
    real(r8), intent(in)  :: dwat(bounds%begc:bounds%endc,1:nlevsoi)
    real(r8), intent(in)  :: imped(bounds%begc:bounds%endc,1:nlevsoi)
    real(r8), intent(in)  :: smp(bounds%begc:bounds%endc,1:nlevsoi)
    real(r8), intent(in)  :: vwc_liq(bounds%begc:bounds%endc,1:nlevsoi)
!    real(r8), intent(in)  :: hksat(bounds%begc:bounds%endc,1:nlevsoi)

    !
    ! !LOCAL VARIABLES:
    integer  :: c,fc,j                                     ! do loop indices
    integer  :: jwt
    real(r8) :: ka
    real(r8) :: wh_zwt
    real(r8) :: wh
    real(r8) :: dtime         ! land model time step (sec)
    real(r8) :: smp1           !temporary variable
    real(r8) :: s1            !temporary variable
    real(r8) :: dz_aquifer        !temporary variable

    character(len=32)  :: subname = 'compute_qcharge'     ! subroutine name   
    !-----------------------------------------------------------------------

    associate(&
         qcharge           =>    soilhydrology_inst%qcharge_col     , & ! Input:  [real(r8) (:)   ]  aquifer recharge rate (mm/s)                      
         wa                => soilhydrology_inst%wa_col             , & ! Input:  [real(r8) (:)   ]  water in the unconfined aquifer (mm)              
         zwt               =>    soilhydrology_inst%zwt_col         , & ! Input:  [real(r8) (:)   ]  water table depth (m)                             
         sucsat            =>    soilstate_inst%sucsat_col          , & ! Input:  [real(r8) (:,:) ]  minimum soil suction (mm)                       
         watsat            =>    soilstate_inst%watsat_col          , & ! Input:  [real(r8) (:,:) ]  volumetric soil water at saturation (porosity)  
         smpmin            =>    soilstate_inst%smpmin_col          , & ! Input:  [real(r8) (:)   ]  restriction for min of soil potential (mm)        
         h2osoi_vol        =>    waterstatebulk_inst%h2osoi_vol_col     , & ! Input:  [real(r8) (:,:) ]  volumetric soil water (0<=h2osoi_vol<=watsat) [m3/m3]
         z                 =>    col%z                              , & ! Input:  [real(r8) (:,:) ]  layer depth (m)                                 
         zi                =>    col%zi                               & ! Input:  [real(r8) (:,:) ]  layer interface depth (m)                                 
         )  ! end associate statement

      ! Get time step

      dtime = get_step_size()

      ! compute flux of water to aquifer
      do fc = 1, num_hydrologyc
         c = filter_hydrologyc(fc)
!  Locate index of layer above water table
         jwt = nlevsoi
         do j = 1,nlevsoi
            if(zwt(c) <= zi(c,j)) then
               jwt = j-1
               exit
            end if
         enddo
         ! calculate qcharge for case jwt < nlevsoi
         if(jwt < nlevsoi) then
            wh_zwt = - sucsat(c,jwt+1) - zwt(c)*m_to_mm

            ! Recharge rate qcharge to groundwater (positive to aquifer)
!            s1 = max(h2osoi_vol(c,jwt+1)/watsat(c,jwt+1), 0.01_r8)
            s1 = max(vwc_liq(c,jwt+1)/watsat(c,jwt+1), 0.01_r8)
            s1 = min(1._r8, s1)

            !this is the expression for unsaturated hk
            call soil_water_retention_curve%soil_hk(c, jwt+1, s1, &
                 imped(c,jwt+1), soilstate_inst, ka)

            ! Recharge rate qcharge to groundwater (positive to aquifer)
            smp1 = max(smpmin(c), smp(c,max(1,jwt)))
            wh      = smp1 - z(c,max(1,jwt))*m_to_mm

            ! original formulation
            if(jwt == 0) then
               qcharge(c) = -ka * (wh_zwt-wh)  /((zwt(c)+1.e-3)*m_to_mm)
            else
               !             qcharge(c) = -ka * (wh_zwt-wh)/((zwt(c)-z(c,jwt))*1000._r8)
               ! 1/2, assuming flux is at zwt interface, saturation deeper than zwt
               qcharge(c) = -ka * (wh_zwt-wh)/((zwt(c)-z(c,jwt))*m_to_mm*2.0)
            endif

            ! To limit qcharge  (for the first several timesteps)
            qcharge(c) = max(-10.0_r8/dtime,qcharge(c))
            qcharge(c) = min( 10.0_r8/dtime,qcharge(c))
         endif

      end do
   
    end associate

   end subroutine compute_qcharge

!#13
  !-----------------------------------------------------------------------
  subroutine IceImpedance(icefrac, e_ice, imped)
    !
    !DESCRIPTION
    ! compute soil suction potential
    !
    ! !USES
    use shr_kind_mod  , only : r8 => shr_kind_r8
    use abortutils    , only : endrun
    !
    ! !ARGUMENTS:
    implicit none
    real(r8), intent(in)  :: icefrac    !fraction of pore space filled with ice
    real(r8), intent(in)  :: e_ice      !shape parameter

    real(r8), intent(out) :: imped      !hydraulic conductivity reduction due to the presence of ice in pore space
    !
    ! !LOCAL VARIABLES:
    character(len=32) :: subname = 'IceImpedance'  ! subroutine name
    !------------------------------------------------------------------------------

    imped = 10._r8**(-e_ice*icefrac)

  end subroutine IceImpedance

!#14
  !-----------------------------------------------------------------------
  subroutine TridiagonalCol (ci, lbj, ubj, jtop, a, b, c, r, u)
    !
    ! !DESCRIPTION:
    ! Tridiagonal matrix solution
    !
    ! !USES:
    use shr_kind_mod   , only : r8 => shr_kind_r8
    use shr_log_mod    , only : errMsg => shr_log_errMsg
    use clm_varpar     , only : nlevurb
    use column_varcon  , only : icol_roof, icol_sunwall, icol_shadewall
    use clm_varctl     , only : iulog
    use decompMod      , only : bounds_type
    use ColumnType     , only : col                
    !
    ! !ARGUMENTS:
    implicit none
    integer , intent(in)    :: ci                 ! lbinning and ubing level indices
    integer , intent(in)    :: lbj, ubj                 ! lbinning and ubing level indices
    integer , intent(in)    :: jtop     ! top level for each column [col]
    real(r8), intent(in)    :: a(  lbj: ) ! "a" left off diagonal of tridiagonal matrix [col, j]
    real(r8), intent(in)    :: b(  lbj: ) ! "b" diagonal column for tridiagonal matrix [col, j]
    real(r8), intent(in)    :: c(  lbj: ) ! "c" right off diagonal tridiagonal matrix [col, j]
    real(r8), intent(in)    :: r(  lbj: ) ! "r" forcing term of tridiagonal matrix [col, j]
    real(r8), intent(inout) :: u(  lbj: ) ! solution [col, j]
    !
    integer  :: j,fc                   !indices
    real(r8) :: gam(lbj:ubj)      !temporary
    real(r8) :: bet               !temporary
    !-----------------------------------------------------------------------

    ! Enforce expected array sizes
    SHR_ASSERT_ALL((ubound(a)    == (/ubj/)), errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(b)    == (/ubj/)), errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(c)    == (/ubj/)), errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(r)    == (/ubj/)), errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(u)    == (/ubj/)), errMsg(sourcefile, __LINE__))

    ! Solve the matrix

    bet = b(jtop)

    do j = lbj, ubj
       if ((col%itype(ci) == icol_sunwall .or. col%itype(ci) == icol_shadewall &
            .or. col%itype(ci) == icol_roof) .and. j <= nlevurb) then
          if (j >= jtop) then
             if (j == jtop) then
                u(j) = r(j) / bet
             else
                gam(j) = c(j-1) / bet
                bet = b(j) - a(j) * gam(j)
                u(j) = (r(j) - a(j)*u(j-1)) / bet
             end if
          end if
       else if (col%itype(ci) /= icol_sunwall .and. col%itype(ci) /= icol_shadewall &
            .and. col%itype(ci) /= icol_roof) then
          if (j >= jtop) then
             if (j == jtop) then
                u(j) = r(j) / bet
             else
                gam(j) = c(j-1) / bet
                bet = b(j) - a(j) * gam(j)
                u(j) = (r(j) - a(j)*u(j-1)) / bet
             end if
          end if
       end if
    end do

    do j = ubj-1,lbj,-1
       if ((col%itype(ci) == icol_sunwall .or. col%itype(ci) == icol_shadewall &
            .or. col%itype(ci) == icol_roof) .and. j <= nlevurb-1) then
          if (j >= jtop) then
             u(j) = u(j) - gam(j+1) * u(j+1)
          end if
       else if (col%itype(ci) /= icol_sunwall .and. col%itype(ci) /= icol_shadewall &
            .and. col%itype(ci) /= icol_roof) then
          if (j >= jtop) then
             u(j) = u(j) - gam(j+1) * u(j+1)
          end if
       end if
    end do
    
  end subroutine TridiagonalCol

 end module SoilWaterMovementMod
