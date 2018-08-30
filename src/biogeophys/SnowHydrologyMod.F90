! -*- mode: f90; indent-tabs-mode: nil; f90-do-indent:3; f90-if-indent:3; f90-type-indent:3; f90-program-indent:2; f90-associate-indent:0; f90-continuation-indent:5  -*-
module SnowHydrologyMod

#include "shr_assert.h"

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Calculate snow hydrology.
  ! - Using as input aerosol deposition from atmosphere model calculate
  !   aerosol fluxes and masses in each layer - need for surface albedo calculation
  ! - Change of snow mass and the snow water onto soil
  ! - Change in snow layer thickness due to compaction
  ! - Combine snow layers less than a min thickness
  ! - Subdivide snow layers if they exceed maximum thickness
  ! - Construct snow/no-snow filters

  !
  ! !USES:
  use shr_kind_mod    , only : r8 => shr_kind_r8
  use shr_log_mod     , only : errMsg => shr_log_errMsg
  use decompMod       , only : bounds_type
  use abortutils      , only : endrun
  use clm_varpar      , only : nlevsno
  use clm_varctl      , only : iulog
  use clm_varcon      , only : namec, h2osno_max
  use atm2lndType     , only : atm2lnd_type
  use AerosolMod      , only : aerosol_type
  use TemperatureType , only : temperature_type
  use WaterFluxBulkType   , only : waterfluxbulk_type
  use WaterStateBulkType  , only : waterstatebulk_type
  use WaterDiagnosticBulkType  , only : waterdiagnosticbulk_type
  use LandunitType    , only : lun
  use TopoMod, only : topo_type
  use ColumnType      , only : col
  use landunit_varcon , only : istsoil, istdlak, istsoil, istwet, istice_mec, istcrop
  use clm_time_manager, only : get_step_size, get_nstep
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: SnowHydrology_readnl       ! Read namelist
  public :: SnowWater                  ! Change of snow mass and the snow water onto soil
  public :: SnowCompaction             ! Change in snow layer thickness due to compaction
  public :: CombineSnowLayers          ! Combine snow layers less than a min thickness
  public :: DivideSnowLayers           ! Subdivide snow layers if they exceed maximum thickness
  public :: InitSnowLayers             ! Initialize cold-start snow layer thickness
  public :: BuildSnowFilter            ! Construct snow/no-snow filters
  public :: SnowCapping                ! Remove snow mass for capped columns
  public :: NewSnowBulkDensity         ! Compute bulk density of any newly-fallen snow

  ! The following are public just for the sake of unit testing:
  public :: SnowCappingExcess          ! Determine the excess snow that needs to be capped
  public :: SnowHydrologySetControlForTesting ! Set some of the control settings
  !
  ! !PRIVATE MEMBER FUNCTIONS:
  private :: Combo                  ! Returns the combined variables: dz, t, wliq, wice.
  private :: MassWeightedSnowRadius ! Mass weighted snow grain size
  !
  ! !PUBLIC DATA MEMBERS:
  !  Aerosol species indices:
  !  1= hydrophillic black carbon
  !  2= hydrophobic black carbon
  !  3= hydrophilic organic carbon
  !  4= hydrophobic organic carbon
  !  5= dust species 1
  !  6= dust species 2
  !  7= dust species 3
  !  8= dust species 4
  !
  real(r8), public, parameter :: scvng_fct_mlt_bcphi = 0.20_r8 ! scavenging factor for hydrophillic BC inclusion in meltwater [frc]
  real(r8), public, parameter :: scvng_fct_mlt_bcpho = 0.03_r8 ! scavenging factor for hydrophobic BC inclusion in meltwater  [frc]
  real(r8), public, parameter :: scvng_fct_mlt_ocphi = 0.20_r8 ! scavenging factor for hydrophillic OC inclusion in meltwater [frc]
  real(r8), public, parameter :: scvng_fct_mlt_ocpho = 0.03_r8 ! scavenging factor for hydrophobic OC inclusion in meltwater  [frc]
  real(r8), public, parameter :: scvng_fct_mlt_dst1  = 0.02_r8 ! scavenging factor for dust species 1 inclusion in meltwater  [frc]
  real(r8), public, parameter :: scvng_fct_mlt_dst2  = 0.02_r8 ! scavenging factor for dust species 2 inclusion in meltwater  [frc]
  real(r8), public, parameter :: scvng_fct_mlt_dst3  = 0.01_r8 ! scavenging factor for dust species 3 inclusion in meltwater  [frc]
  real(r8), public, parameter :: scvng_fct_mlt_dst4  = 0.01_r8 ! scavenging factor for dust species 4 inclusion in meltwater  [frc]

  ! The following are public for the sake of unit testing
  integer, parameter, public :: LoTmpDnsSlater2017            = 2    ! For temperature below -15C use equation from Slater 2017
  integer, parameter, public :: LoTmpDnsTruncatedAnderson1976 = 1    ! Truncate low temp. snow density from the Anderson-1976 version at -15C

  ! Definition of snow pack vertical structure
  ! Hardcoded maximum of 12 snowlayers, this is checked elsewhere (controlMod.F90)
  ! The bottom layer has no limit on thickness, hence the last element of the dzmax_*
  ! arrays is 'huge'.
  real(r8), parameter :: dzmin(12) = &       ! minimum of top snow layer
               (/ 0.010_r8, 0.015_r8, 0.025_r8, 0.055_r8, 0.115_r8, 0.235_r8, &
                  0.475_r8, 0.955_r8, 1.915_r8, 3.835_r8, 7.675_r8, 15.355_r8 /)
  real(r8), parameter :: dzmax_l(12) = &     ! maximum thickness of layer when no layers beneath
               (/ 0.03_r8, 0.07_r8, 0.18_r8, 0.41_r8, 0.88_r8, 1.83_r8, &
                  3.74_r8, 7.57_r8, 15.24_r8, 30.59_r8, 61.3_r8, huge(1._r8)  /)
  real(r8), parameter :: dzmax_u(12) = &     ! maximum thickness of layer when layers beneath
               (/ 0.02_r8, 0.05_r8, 0.11_r8, 0.23_r8, 0.47_r8, 0.95_r8, &
                  1.91_r8, 3.83_r8, 7.67_r8, 15.35_r8, 30.71_r8, huge(1._r8)  /)

  !
  ! !PRIVATE DATA MEMBERS:

  integer, parameter :: OverburdenCompactionMethodAnderson1976 = 1
  integer, parameter :: OverburdenCompactionMethodVionnet2012  = 2

  ! If true, the density of new snow depends on wind speed, and there is also
  ! wind-dependent snow compaction
  logical  :: wind_dependent_snow_density                      ! If snow density depends on wind or not
  integer  :: overburden_compaction_method = -1
  integer  :: new_snow_density            = LoTmpDnsSlater2017 ! Snow density type
  real(r8) :: upplim_destruct_metamorph   = 100.0_r8           ! Upper Limit on Destructive Metamorphism Compaction [kg/m3]
  real(r8) :: overburden_compress_Tfactor = 0.08_r8            ! snow compaction overburden exponential factor (1/K)
  real(r8) :: min_wind_snowcompact        = 5._r8              ! minimum wind speed tht results in compaction (m/s)

  ! ------------------------------------------------------------------------
  ! Parameters controlling the resetting of the snow pack
  ! ------------------------------------------------------------------------

  logical  :: reset_snow      = .false.  ! If set to true, we reset the non-glc snow pack, based on the following parameters
  logical  :: reset_snow_glc  = .false.  ! If set to true, we reset the glc snow pack, based reset_snow_glc_ela

  ! Default for reset_snow_glc_ela implies that snow will be reset for all glacier columns if reset_snow_glc = .true.
  real(r8) :: reset_snow_glc_ela = 1.e9_r8  ! equilibrium line altitude (m); snow is reset for glacier columns below this elevation if reset_snow_glc = .true. (ignored if reset_snow_glc = .false.)

  ! The following are public simply to support unit testing

  ! 35 mm was chosen by Raymond Sellevold, based on finding the location in Greenland
  ! with the least amount of snowfall from Sept. 1 (roughly the end of the melt season)
  ! and Jan. 1 (when we typically start simulations). This location with the least amount
  ! of snowfall had an average of 35 mm snow fall over this 4-month period.
  real(r8), parameter, public :: reset_snow_h2osno = 35._r8  ! mm SWE to reset the snow pack to

  ! We scale the number of reset time steps with the number of snow layers, since we can
  ! remove up to one layer per time step. In the absence of snow accumulation, we might
  ! be able to get away with 1 reset time step per layer. However, we specify a larger
  ! number to be more robust.
  real(r8), parameter, public :: reset_snow_timesteps_per_layer = 4

  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine SnowHydrology_readnl( NLFilename)
    !
    ! !DESCRIPTION:
    ! Read the namelist for SnowHydrology
    !
    ! !USES:
    use fileutils      , only : getavu, relavu, opnfil
    use shr_nl_mod     , only : shr_nl_find_group_name
    use spmdMod        , only : masterproc, mpicom
    use shr_mpi_mod    , only : shr_mpi_bcast
    !
    ! !ARGUMENTS:
    character(len=*), intent(in) :: NLFilename ! Namelist filename
    !
    ! !LOCAL VARIABLES:
    integer :: ierr                 ! error code
    integer :: unitn                ! unit for namelist file
    character(len=64) :: snow_overburden_compaction_method
    character(len=25) :: lotmp_snowdensity_method

    character(len=*), parameter :: subname = 'SnowHydrology_readnl'
    !-----------------------------------------------------------------------

    namelist /clm_snowhydrology_inparm/ &
         wind_dependent_snow_density, snow_overburden_compaction_method, &
         lotmp_snowdensity_method, upplim_destruct_metamorph, &
         overburden_compress_Tfactor, min_wind_snowcompact, &
         reset_snow, reset_snow_glc, reset_snow_glc_ela

    ! Initialize options to default values, in case they are not specified in the namelist
    wind_dependent_snow_density = .false.
    snow_overburden_compaction_method = ' '

    if (masterproc) then
       unitn = getavu()
       write(iulog,*) 'Read in clm_SnowHydrology_inparm  namelist'
       call opnfil (NLFilename, unitn, 'F')
       call shr_nl_find_group_name(unitn, 'clm_SnowHydrology_inparm', status=ierr)
       if (ierr == 0) then
          read(unitn, clm_snowhydrology_inparm, iostat=ierr)
          if (ierr /= 0) then
             call endrun(msg="ERROR reading clm_snowhydrology_inparm namelist"//errmsg(sourcefile, __LINE__))
          end if
       else
          call endrun(msg="ERROR finding clm_snowhydrology_inparm namelist"//errmsg(sourcefile, __LINE__))
       end if
       call relavu( unitn )
    end if

    call shr_mpi_bcast (wind_dependent_snow_density, mpicom)
    call shr_mpi_bcast (snow_overburden_compaction_method, mpicom)
    call shr_mpi_bcast (lotmp_snowdensity_method   , mpicom)
    call shr_mpi_bcast (upplim_destruct_metamorph  , mpicom)
    call shr_mpi_bcast (overburden_compress_Tfactor, mpicom)
    call shr_mpi_bcast (min_wind_snowcompact       , mpicom)
    call shr_mpi_bcast (reset_snow                 , mpicom)
    call shr_mpi_bcast (reset_snow_glc             , mpicom)
    call shr_mpi_bcast (reset_snow_glc_ela         , mpicom)

    if (masterproc) then
       write(iulog,*) ' '
       write(iulog,*) 'SnowHydrology settings:'
       write(iulog,nml=clm_snowhydrology_inparm)
       write(iulog,*) ' '
    end if

    if (      trim(lotmp_snowdensity_method) == 'Slater2017' ) then
       new_snow_density = LoTmpDnsSlater2017
    else if ( trim(lotmp_snowdensity_method) == 'TruncatedAnderson1976' ) then
       new_snow_density = LoTmpDnsTruncatedAnderson1976
    else
       call endrun(msg="ERROR bad lotmp_snowdensity_method name"//errmsg(sourcefile, __LINE__))
    end if

    if (trim(snow_overburden_compaction_method) == 'Anderson1976') then
       overburden_compaction_method = OverburdenCompactionMethodAnderson1976
    else if (trim(snow_overburden_compaction_method) == 'Vionnet2012') then
       overburden_compaction_method = OverburdenCompactionMethodVionnet2012
    else
       call endrun(msg="ERROR bad snow_overburden_compaction_method name"// &
            errMsg(sourcefile, __LINE__))
    end if

  end subroutine SnowHydrology_readnl


  !-----------------------------------------------------------------------
  subroutine SnowWater(bounds, &
       num_snowc, filter_snowc, num_nosnowc, filter_nosnowc, &
       atm2lnd_inst, waterfluxbulk_inst, waterstatebulk_inst, waterdiagnosticbulk_inst, aerosol_inst)
    !
    ! !DESCRIPTION:
    ! Evaluate the change of snow mass and the snow water onto soil.
    ! Water flow within snow is computed by an explicit and non-physical
    ! based scheme, which permits a part of liquid water over the holding
    ! capacity (a tentative value is used, i.e. equal to 0.033*porosity) to
    ! percolate into the underlying layer.  Except for cases where the
    ! porosity of one of the two neighboring layers is less than 0.05, zero
    ! flow is assumed. The water flow out of the bottom of the snow pack will
    ! participate as the input of the soil water and runoff.  This subroutine
    ! uses a filter for columns containing snow which must be constructed prior
    ! to being called.
    !
    ! !USES:
    use clm_varcon        , only : denh2o, denice, wimp, ssi
    use AerosolMod        , only : AerosolFluxes
    !
    ! !ARGUMENTS:
    type(bounds_type)     , intent(in)    :: bounds
    integer               , intent(in)    :: num_snowc         ! number of snow points in column filter
    integer               , intent(in)    :: filter_snowc(:)   ! column filter for snow points
    integer               , intent(in)    :: num_nosnowc       ! number of non-snow points in column filter
    integer               , intent(in)    :: filter_nosnowc(:) ! column filter for non-snow points
    type(atm2lnd_type)    , intent(in)    :: atm2lnd_inst
    type(waterfluxbulk_type)  , intent(inout) :: waterfluxbulk_inst
    type(waterstatebulk_type) , intent(inout) :: waterstatebulk_inst
    type(waterdiagnosticbulk_type) , intent(inout) :: waterdiagnosticbulk_inst
    type(aerosol_type)    , intent(inout) :: aerosol_inst
    !
    ! !LOCAL VARIABLES:
    integer  :: g                                                  ! gridcell loop index
    integer  :: c, j, fc, l                                        ! do loop/array indices
    real(r8) :: dtime                                              ! land model time step (sec)
    real(r8) :: qin(bounds%begc:bounds%endc)                       ! water flow into the element (mm/s)
    real(r8) :: qout(bounds%begc:bounds%endc)                      ! water flow out of the elmement (mm/s)
    real(r8) :: qin_bc_phi  (bounds%begc:bounds%endc)              ! flux of hydrophilic BC into   layer [kg]
    real(r8) :: qout_bc_phi (bounds%begc:bounds%endc)              ! flux of hydrophilic BC out of layer [kg]
    real(r8) :: qin_bc_pho  (bounds%begc:bounds%endc)              ! flux of hydrophobic BC into   layer [kg]
    real(r8) :: qout_bc_pho (bounds%begc:bounds%endc)              ! flux of hydrophobic BC out of layer [kg]
    real(r8) :: qin_oc_phi  (bounds%begc:bounds%endc)              ! flux of hydrophilic OC into   layer [kg]
    real(r8) :: qout_oc_phi (bounds%begc:bounds%endc)              ! flux of hydrophilic OC out of layer [kg]
    real(r8) :: qin_oc_pho  (bounds%begc:bounds%endc)              ! flux of hydrophobic OC into   layer [kg]
    real(r8) :: qout_oc_pho (bounds%begc:bounds%endc)              ! flux of hydrophobic OC out of layer [kg]
    real(r8) :: qin_dst1    (bounds%begc:bounds%endc)              ! flux of dust species 1 into   layer [kg]
    real(r8) :: qout_dst1   (bounds%begc:bounds%endc)              ! flux of dust species 1 out of layer [kg]
    real(r8) :: qin_dst2    (bounds%begc:bounds%endc)              ! flux of dust species 2 into   layer [kg]
    real(r8) :: qout_dst2   (bounds%begc:bounds%endc)              ! flux of dust species 2 out of layer [kg]
    real(r8) :: qin_dst3    (bounds%begc:bounds%endc)              ! flux of dust species 3 into   layer [kg]
    real(r8) :: qout_dst3   (bounds%begc:bounds%endc)              ! flux of dust species 3 out of layer [kg]
    real(r8) :: qin_dst4    (bounds%begc:bounds%endc)              ! flux of dust species 4 into   layer [kg]
    real(r8) :: qout_dst4   (bounds%begc:bounds%endc)              ! flux of dust species 4 out of layer [kg]
    real(r8) :: wgdif                                              ! ice mass after minus sublimation
    real(r8) :: vol_liq(bounds%begc:bounds%endc,-nlevsno+1:0)      ! partial volume of liquid water in layer
    real(r8) :: vol_ice(bounds%begc:bounds%endc,-nlevsno+1:0)      ! partial volume of ice lens in layer
    real(r8) :: eff_porosity(bounds%begc:bounds%endc,-nlevsno+1:0) ! effective porosity = porosity - vol_ice
    real(r8) :: mss_liqice(bounds%begc:bounds%endc,-nlevsno+1:0)   ! mass of liquid+ice in a layer
    !-----------------------------------------------------------------------

    associate( &
         dz             => col%dz                            , & ! Input:  [real(r8) (:,:) ] layer depth (m)
         snl            => col%snl                           , & ! Input:  [integer  (:)   ] number of snow layers

         frac_sno_eff   => waterdiagnosticbulk_inst%frac_sno_eff_col  , & ! Input:  [real(r8) (:)   ] eff. fraction of ground covered by snow (0 to 1)
         frac_sno       => waterdiagnosticbulk_inst%frac_sno_col      , & ! Input:  [real(r8) (:)   ] fraction of ground covered by snow (0 to 1)
         h2osno         => waterstatebulk_inst%h2osno_col        , & ! Input:  [real(r8) (:)   ] snow water (mm H2O)
         int_snow       => waterstatebulk_inst%int_snow_col      , & ! Output: [real(r8) (:)   ] integrated snowfall [mm]
         h2osoi_ice     => waterstatebulk_inst%h2osoi_ice_col    , & ! Output: [real(r8) (:,:) ] ice lens (kg/m2)
         h2osoi_liq     => waterstatebulk_inst%h2osoi_liq_col    , & ! Output: [real(r8) (:,:) ] liquid water (kg/m2)

         qflx_snomelt   => waterfluxbulk_inst%qflx_snomelt_col   , & ! Input:  [real(r8) (:)   ] snow melt (mm H2O /s)
         qflx_rain_grnd => waterfluxbulk_inst%qflx_rain_grnd_col , & ! Input:  [real(r8) (:)   ] rain on ground after interception (mm H2O/s) [+]
         qflx_sub_snow  => waterfluxbulk_inst%qflx_sub_snow_col  , & ! Input:  [real(r8) (:)   ] sublimation rate from snow pack (mm H2O /s) [+]
         qflx_dew_snow  => waterfluxbulk_inst%qflx_dew_snow_col  , & ! Input:  [real(r8) (:)   ] surface dew added to snow pack (mm H2O /s) [+]
         qflx_evap_grnd => waterfluxbulk_inst%qflx_evap_grnd_col , & ! Input:  [real(r8) (:)   ] ground surface evaporation rate (mm H2O/s) [+]
         qflx_dew_grnd  => waterfluxbulk_inst%qflx_dew_grnd_col  , & ! Input:  [real(r8) (:)   ] ground surface dew formation (mm H2O /s) [+]
         qflx_snow_drain => waterfluxbulk_inst%qflx_snow_drain_col,& ! Output: [real(r8) (:)   ] net snow melt
         qflx_rain_plus_snomelt => waterfluxbulk_inst%qflx_rain_plus_snomelt_col , & ! Output: [real(r8) (:)   ] rain plus snow melt falling on the soil (mm/s)
         snow_depth     => waterdiagnosticbulk_inst%snow_depth_col    , & ! Output: [real(r8) (:)   ] snow height (m)

         mss_bcphi      => aerosol_inst%mss_bcphi_col        , & ! Output: [real(r8) (:,:) ] hydrophillic BC mass in snow (col,lyr) [kg]
         mss_bcpho      => aerosol_inst%mss_bcpho_col        , & ! Output: [real(r8) (:,:) ] hydrophobic  BC mass in snow (col,lyr) [kg]
         mss_ocphi      => aerosol_inst%mss_ocphi_col        , & ! Output: [real(r8) (:,:) ] hydrophillic OC mass in snow (col,lyr) [kg]
         mss_ocpho      => aerosol_inst%mss_ocpho_col        , & ! Output: [real(r8) (:,:) ] hydrophobic  OC mass in snow (col,lyr) [kg]
         mss_dst1       => aerosol_inst%mss_dst1_col         , & ! Output: [real(r8) (:,:) ] mass of dust species 1 in snow (col,lyr) [kg]
         mss_dst2       => aerosol_inst%mss_dst2_col         , & ! Output: [real(r8) (:,:) ] mass of dust species 2 in snow (col,lyr) [kg]
         mss_dst3       => aerosol_inst%mss_dst3_col         , & ! Output: [real(r8) (:,:) ] mass of dust species 3 in snow (col,lyr) [kg]
         mss_dst4       => aerosol_inst%mss_dst4_col         , & ! Output: [real(r8) (:,:) ] mass of dust species 4 in snow (col,lyr) [kg]

         begc           => bounds%begc                       , &
         endc           => bounds%endc                         &
    )

    ! Determine model time step

    dtime = get_step_size()

    ! Renew the mass of ice lens (h2osoi_ice) and liquid (h2osoi_liq) in the
    ! surface snow layer resulting from sublimation (frost) / evaporation (condense)

    do fc = 1,num_snowc
       c = filter_snowc(fc)
       l=col%landunit(c)

       wgdif = h2osoi_ice(c,snl(c)+1) &
            + frac_sno_eff(c) * (qflx_dew_snow(c) - qflx_sub_snow(c)) * dtime
       h2osoi_ice(c,snl(c)+1) = wgdif
       if (wgdif < 0._r8) then
          h2osoi_ice(c,snl(c)+1) = 0._r8
          h2osoi_liq(c,snl(c)+1) = h2osoi_liq(c,snl(c)+1) + wgdif
       end if
       h2osoi_liq(c,snl(c)+1) = h2osoi_liq(c,snl(c)+1) +  &
            frac_sno_eff(c) * (qflx_rain_grnd(c) + qflx_dew_grnd(c) &
            - qflx_evap_grnd(c)) * dtime

       ! if negative, reduce deeper layer's liquid water content sequentially
       if(h2osoi_liq(c,snl(c)+1) < 0._r8) then
          do j = snl(c)+1, 1
             wgdif=h2osoi_liq(c,j)
             if (wgdif >= 0._r8) exit
             h2osoi_liq(c,j) = 0._r8
             h2osoi_liq(c,j+1) = h2osoi_liq(c,j+1) + wgdif
          enddo
       end if
    end do

    ! Porosity and partial volume

    do j = -nlevsno+1, 0
       do fc = 1, num_snowc
          c = filter_snowc(fc)
          if (j >= snl(c)+1) then
             ! need to scale dz by frac_sno to convert to grid cell average depth
             vol_ice(c,j)      = min(1._r8, h2osoi_ice(c,j)/(dz(c,j)*frac_sno_eff(c)*denice))
             eff_porosity(c,j) = 1._r8 - vol_ice(c,j)
             vol_liq(c,j)      = min(eff_porosity(c,j),h2osoi_liq(c,j)/(dz(c,j)*frac_sno_eff(c)*denh2o))
          end if
       end do
    end do

    ! Capillary forces within snow are usually two or more orders of magnitude
    ! less than those of gravity. Only gravity terms are considered.
    ! the genernal expression for water flow is "K * ss**3", however,
    ! no effective parameterization for "K".  Thus, a very simple consideration
    ! (not physically based) is introduced:
    ! when the liquid water of layer exceeds the layer's holding
    ! capacity, the excess meltwater adds to the underlying neighbor layer.

    ! Also compute aerosol fluxes through snowpack in this loop:
    ! 1) compute aerosol mass in each layer
    ! 2) add aerosol mass flux from above layer to mass of this layer
    ! 3) qout_xxx is mass flux of aerosol species xxx out bottom of
    !    layer in water flow, proportional to (current) concentration
    !    of aerosol in layer multiplied by a scavenging ratio.
    ! 4) update mass of aerosol in top layer, accordingly
    ! 5) update mass concentration of aerosol accordingly

    do c = bounds%begc,bounds%endc
       qin(c)         = 0._r8
       qin_bc_phi (c) = 0._r8
       qin_bc_pho (c) = 0._r8
       qin_oc_phi (c) = 0._r8
       qin_oc_pho (c) = 0._r8
       qin_dst1   (c) = 0._r8
       qin_dst2   (c) = 0._r8
       qin_dst3   (c) = 0._r8
       qin_dst4   (c) = 0._r8
    end do

    do j = -nlevsno+1, 0
       do fc = 1, num_snowc
          c = filter_snowc(fc)
          if (j >= snl(c)+1) then

             h2osoi_liq(c,j) = h2osoi_liq(c,j) + qin(c)

             mss_bcphi(c,j) = mss_bcphi(c,j) + qin_bc_phi(c)
             mss_bcpho(c,j) = mss_bcpho(c,j) + qin_bc_pho(c)
             mss_ocphi(c,j) = mss_ocphi(c,j) + qin_oc_phi(c)
             mss_ocpho(c,j) = mss_ocpho(c,j) + qin_oc_pho(c)

             mss_dst1(c,j)  = mss_dst1(c,j) + qin_dst1(c)
             mss_dst2(c,j)  = mss_dst2(c,j) + qin_dst2(c)
             mss_dst3(c,j)  = mss_dst3(c,j) + qin_dst3(c)
             mss_dst4(c,j)  = mss_dst4(c,j) + qin_dst4(c)

             if (j <= -1) then
                ! No runoff over snow surface, just ponding on surface
                if (eff_porosity(c,j) < wimp .OR. eff_porosity(c,j+1) < wimp) then
                   qout(c) = 0._r8
                else
                   ! dz must be scaled by frac_sno to obtain gridcell average value
                   qout(c) = max(0._r8,(vol_liq(c,j) &
                        - ssi*eff_porosity(c,j))*dz(c,j)*frac_sno_eff(c))
                   qout(c) = min(qout(c),(1._r8-vol_ice(c,j+1) &
                        - vol_liq(c,j+1))*dz(c,j+1)*frac_sno_eff(c))
                end if
             else
                qout(c) = max(0._r8,(vol_liq(c,j) &
                     - ssi*eff_porosity(c,j))*dz(c,j)*frac_sno_eff(c))
             end if
             qout(c) = qout(c)*1000._r8
             h2osoi_liq(c,j) = h2osoi_liq(c,j) - qout(c)
             qin(c) = qout(c)

             ! mass of ice+water: in extremely rare circumstances, this can
             ! be zero, even though there is a snow layer defined. In
             ! this case, set the mass to a very small value to
             ! prevent division by zero.

             mss_liqice(c,j) = h2osoi_liq(c,j)+h2osoi_ice(c,j)
             if (mss_liqice(c,j) < 1E-30_r8) then
                mss_liqice(c,j) = 1E-30_r8
             endif

             ! BCPHI:
             ! 1. flux with meltwater:
             qout_bc_phi(c) = qout(c)*scvng_fct_mlt_bcphi*(mss_bcphi(c,j)/mss_liqice(c,j))
             if (qout_bc_phi(c) > mss_bcphi(c,j)) then
                qout_bc_phi(c) = mss_bcphi(c,j)
             endif
             mss_bcphi(c,j) = mss_bcphi(c,j) - qout_bc_phi(c)
             qin_bc_phi(c) = qout_bc_phi(c)

             ! BCPHO:
             ! 1. flux with meltwater:
             qout_bc_pho(c) = qout(c)*scvng_fct_mlt_bcpho*(mss_bcpho(c,j)/mss_liqice(c,j))
             if (qout_bc_pho(c) > mss_bcpho(c,j)) then
                qout_bc_pho(c) = mss_bcpho(c,j)
             endif
             mss_bcpho(c,j) = mss_bcpho(c,j) - qout_bc_pho(c)
             qin_bc_pho(c) = qout_bc_pho(c)

             ! OCPHI:
             ! 1. flux with meltwater:
             qout_oc_phi(c) = qout(c)*scvng_fct_mlt_ocphi*(mss_ocphi(c,j)/mss_liqice(c,j))
             if (qout_oc_phi(c) > mss_ocphi(c,j)) then
                qout_oc_phi(c) = mss_ocphi(c,j)
             endif
             mss_ocphi(c,j) = mss_ocphi(c,j) - qout_oc_phi(c)
             qin_oc_phi(c) = qout_oc_phi(c)

             ! OCPHO:
             ! 1. flux with meltwater:
             qout_oc_pho(c) = qout(c)*scvng_fct_mlt_ocpho*(mss_ocpho(c,j)/mss_liqice(c,j))
             if (qout_oc_pho(c) > mss_ocpho(c,j)) then
                qout_oc_pho(c) = mss_ocpho(c,j)
             endif
             mss_ocpho(c,j) = mss_ocpho(c,j) - qout_oc_pho(c)
             qin_oc_pho(c) = qout_oc_pho(c)

             ! DUST 1:
             ! 1. flux with meltwater:
             qout_dst1(c) = qout(c)*scvng_fct_mlt_dst1*(mss_dst1(c,j)/mss_liqice(c,j))
             if (qout_dst1(c) > mss_dst1(c,j)) then
                qout_dst1(c) = mss_dst1(c,j)
             endif
             mss_dst1(c,j) = mss_dst1(c,j) - qout_dst1(c)
             qin_dst1(c) = qout_dst1(c)

             ! DUST 2:
             ! 1. flux with meltwater:
             qout_dst2(c) = qout(c)*scvng_fct_mlt_dst2*(mss_dst2(c,j)/mss_liqice(c,j))
             if (qout_dst2(c) > mss_dst2(c,j)) then
                qout_dst2(c) = mss_dst2(c,j)
             endif
             mss_dst2(c,j) = mss_dst2(c,j) - qout_dst2(c)
             qin_dst2(c) = qout_dst2(c)

             ! DUST 3:
             ! 1. flux with meltwater:
             qout_dst3(c) = qout(c)*scvng_fct_mlt_dst3*(mss_dst3(c,j)/mss_liqice(c,j))
             if (qout_dst3(c) > mss_dst3(c,j)) then
                qout_dst3(c) = mss_dst3(c,j)
             endif
             mss_dst3(c,j) = mss_dst3(c,j) - qout_dst3(c)
             qin_dst3(c) = qout_dst3(c)

             ! DUST 4:
             ! 1. flux with meltwater:
             qout_dst4(c) = qout(c)*scvng_fct_mlt_dst4*(mss_dst4(c,j)/mss_liqice(c,j))
             if (qout_dst4(c) > mss_dst4(c,j)) then
                qout_dst4(c) = mss_dst4(c,j)
             endif
             mss_dst4(c,j) = mss_dst4(c,j) - qout_dst4(c)
             qin_dst4(c) = qout_dst4(c)

          end if
       end do
    end do

    ! Compute aerosol fluxes through snowpack and aerosol deposition fluxes into top layere

    call AerosolFluxes(bounds, num_snowc, filter_snowc, &
         atm2lnd_inst, aerosol_inst)

    ! Adjust layer thickness for any water+ice content changes in excess of previous
    ! layer thickness. Strictly speaking, only necessary for top snow layer, but doing
    ! it for all snow layers will catch problems with older initial files.
    ! Layer interfaces (zi) and node depths (z) do not need adjustment here because they
    ! are adjusted in CombineSnowLayers and are not used up to that point.

    do j = -nlevsno+1, 0
       do fc = 1, num_snowc
          c = filter_snowc(fc)
          if (j >= snl(c)+1) then
             dz(c,j) = max(dz(c,j),h2osoi_liq(c,j)/denh2o + h2osoi_ice(c,j)/denice)
          end if
       end do
    end do

    do fc = 1, num_snowc
       c = filter_snowc(fc)
       ! Qout from snow bottom
       qflx_snow_drain(c) = qflx_snow_drain(c) + (qout(c) / dtime)

       qflx_rain_plus_snomelt(c) = (qout(c) / dtime) &
            + (1.0_r8 - frac_sno_eff(c)) * qflx_rain_grnd(c)
       int_snow(c) = int_snow(c) + frac_sno_eff(c) &
                     * (qflx_dew_snow(c) + qflx_dew_grnd(c) + qflx_rain_grnd(c)) * dtime
    end do

    do fc = 1, num_nosnowc
       c = filter_nosnowc(fc)
       qflx_snow_drain(c) = qflx_snomelt(c)

       qflx_rain_plus_snomelt(c) = qflx_rain_grnd(c) + qflx_snomelt(c)
       ! reset accumulated snow when no snow present
       if (h2osno(c) <= 0._r8) then
          int_snow(c) = 0._r8
          frac_sno(c) = 0._r8
          snow_depth(c) = 0._r8
       end if
    end do

    end associate
  end subroutine SnowWater

  !-----------------------------------------------------------------------
  subroutine SnowCompaction(bounds, num_snowc, filter_snowc, &
       temperature_inst, waterstatebulk_inst, waterdiagnosticbulk_inst, atm2lnd_inst)
    !
    ! !DESCRIPTION:
    ! Determine the change in snow layer thickness due to compaction and
    ! settling.
    ! Three metamorphisms of changing snow characteristics are implemented,
    ! i.e., destructive, overburden, and melt. The treatments of the former
    ! two are from SNTHERM.89 and SNTHERM.99 (1991, 1999). The contribution
    ! due to melt metamorphism is simply taken as a ratio of snow ice
    ! fraction after the melting versus before the melting.
    !
    ! !USES:
    use clm_varcon      , only : denice, denh2o, tfrz, rpi, int_snow_max
    use clm_varctl      , only : subgridflag
    !
    ! !ARGUMENTS:
    type(bounds_type)      , intent(in) :: bounds
    integer                , intent(in) :: num_snowc       ! number of column snow points in column filter
    integer                , intent(in) :: filter_snowc(:) ! column filter for snow points
    type(temperature_type) , intent(in) :: temperature_inst
    type(waterstatebulk_type)  , intent(in) :: waterstatebulk_inst
    type(waterdiagnosticbulk_type)  , intent(in) :: waterdiagnosticbulk_inst
    type(atm2lnd_type)     , intent(in) :: atm2lnd_inst
    !
    ! !LOCAL VARIABLES:
    integer :: j, l, c, fc                      ! indices
    integer :: g                                ! gridcell index
    real(r8):: dtime                            ! land model time step (sec)
    ! parameters
    real(r8), parameter :: c3 = 2.777e-6_r8     ! [1/s]
    real(r8), parameter :: c4 = 0.04_r8         ! [1/K]
    real(r8), parameter :: c5 = 2.0_r8          !
    !
    real(r8) :: burden(bounds%begc:bounds%endc)  ! pressure of overlying snow [kg/m2]
    real(r8) :: zpseudo(bounds%begc:bounds%endc) ! wind drift compaction / pseudo depth (only valid if wind_dependent_snow_density is .true.)
    logical  :: mobile(bounds%begc:bounds%endc)  ! current snow layer is mobile, i.e. susceptible to wind drift (only valid if wind_dependent_snow_density is .true.)
    real(r8) :: ddz1   ! Rate of settling of snowpack due to destructive metamorphism.
    real(r8) :: ddz2   ! Rate of compaction of snowpack due to overburden.
    real(r8) :: ddz3   ! Rate of compaction of snowpack due to melt [1/s]
    real(r8) :: dexpf  ! expf=exp(-c4*(273.15-t_soisno)).
    real(r8) :: fi     ! Fraction of ice relative to the total water content at current time step
    real(r8) :: td     ! t_soisno - tfrz [K]
    real(r8) :: pdzdtc ! Nodal rate of change in fractional-thickness due to compaction [fraction/s]
    real(r8) :: void   ! void (1 - vol_ice - vol_liq)
    real(r8) :: wx     ! water mass (ice+liquid) [kg/m2]
    real(r8) :: bi     ! partial density of ice [kg/m3]
    real(r8) :: wsum   ! snowpack total water mass (ice+liquid) [kg/m2]
    real(r8) :: fsno_melt
    real(r8) :: ddz4   ! Rate of compaction of snowpack due to wind drift.
    real(r8) :: int_snow_limited ! integrated snowfall, limited to be no greater than int_snow_max [mm]
    !-----------------------------------------------------------------------

    associate( &
         snl          => col%snl                          , & ! Input:  [integer (:)    ] number of snow layers
         n_melt       => col%n_melt                       , & ! Input:  [real(r8) (:)   ] SCA shape parameter
         lakpoi       => lun%lakpoi                       , & ! Input:  [logical  (:)   ] true => landunit is a lake point
         urbpoi       => lun%urbpoi                       , & ! Input:  [logical  (:)   ] true => landunit is an urban point
         forc_wind    => atm2lnd_inst%forc_wind_grc       , & ! Input:  [real(r8) (:)   ]  atmospheric wind speed (m/s)

         t_soisno     => temperature_inst%t_soisno_col    , & ! Input:  [real(r8) (:,:) ] soil temperature (Kelvin)
         imelt        => temperature_inst%imelt_col       , & ! Input:  [integer (:,:)  ] flag for melting (=1), freezing (=2), Not=0

         frac_sno     => waterdiagnosticbulk_inst%frac_sno_eff_col , & ! Input:  [real(r8) (:)   ] snow covered fraction
         swe_old      => waterdiagnosticbulk_inst%swe_old_col      , & ! Input:  [real(r8) (:,:) ] initial swe values
         int_snow     => waterstatebulk_inst%int_snow_col     , & ! Input:  [real(r8) (:)   ] integrated snowfall [mm]
         frac_iceold  => waterdiagnosticbulk_inst%frac_iceold_col  , & ! Input:  [real(r8) (:,:) ] fraction of ice relative to the tot water
         h2osoi_ice   => waterstatebulk_inst%h2osoi_ice_col   , & ! Input:  [real(r8) (:,:) ] ice lens (kg/m2)
         h2osoi_liq   => waterstatebulk_inst%h2osoi_liq_col   , & ! Input:  [real(r8) (:,:) ] liquid water (kg/m2)

         dz           => col%dz                             & ! Output: [real(r8) (: ,:) ] layer depth (m)
    )

    ! Get time step

    dtime = get_step_size()

    ! Begin calculation - note that the following column loops are only invoked if snl(c) < 0

    do fc = 1, num_snowc
       c = filter_snowc(fc)
       
       burden(c)  = 0._r8
       zpseudo(c) = 0._r8
       mobile(c)  = .true.
    end do

    do j = -nlevsno+1, 0
       do fc = 1, num_snowc
          c = filter_snowc(fc)
          g = col%gridcell(c)
          if (j >= snl(c)+1) then

             wx = (h2osoi_ice(c,j) + h2osoi_liq(c,j))
             void = 1._r8 - (h2osoi_ice(c,j)/denice + h2osoi_liq(c,j)/denh2o)&
                  /(frac_sno(c) * dz(c,j))

             ! Allow compaction only for non-saturated node and higher ice lens node.
             if (void > 0.001_r8 .and. h2osoi_ice(c,j) > .1_r8) then

                bi = h2osoi_ice(c,j) / (frac_sno(c) * dz(c,j))
                fi = h2osoi_ice(c,j) / wx
                td = tfrz-t_soisno(c,j)
                dexpf = exp(-c4*td)

                ! Settling as a result of destructive metamorphism

                ddz1 = -c3*dexpf
                if (bi > upplim_destruct_metamorph) ddz1 = ddz1*exp(-46.0e-3_r8*(bi-upplim_destruct_metamorph))

                ! Liquid water term

                if (h2osoi_liq(c,j) > 0.01_r8*dz(c,j)*frac_sno(c)) ddz1=ddz1*c5

                select case (overburden_compaction_method)
                case (OverburdenCompactionMethodAnderson1976)
                   ddz2 = OverburdenCompactionAnderson1976( &
                        burden = burden(c), &
                        wx = wx, &
                        td = td, &
                        bi = bi)

                case (OverburdenCompactionMethodVionnet2012)
                   ddz2 = OverburdenCompactionVionnet2012( &
                        h2osoi_liq = h2osoi_liq(c,j), &
                        dz = dz(c,j), &
                        burden = burden(c), &
                        wx = wx, &
                        td = td, &
                        bi = bi)

                case default
                   call endrun(msg="Unknown overburden_compaction_method")
                end select

                ! Compaction occurring during melt

                if (imelt(c,j) == 1) then
                   l = col%landunit(c)
                   ! For consistency with other uses of subgridflag==1 (e.g., in
                   ! CanopyHydrologyMod), we apply this code over all landunits other
                   ! than lake and urban. (In CanopyHydrologyMod, the uses of subgridflag
                   ! are in a nolake filter, and check .not. urbpoi.)
                   if(subgridflag==1 .and. (.not. lakpoi(l) .and. .not. urbpoi(l))) then
                      ! first term is delta mass over mass
                      ddz3 = max(0._r8,min(1._r8,(swe_old(c,j) - wx)/wx))

                      ! 2nd term is delta fsno over fsno, allowing for negative values for ddz3
                      if((swe_old(c,j) - wx) > 0._r8) then
                         wsum = sum(h2osoi_liq(c,snl(c)+1:0)+h2osoi_ice(c,snl(c)+1:0))
                         int_snow_limited = min(int_snow(c), int_snow_max)
                         fsno_melt = 1. - (acos(2.*min(1._r8,wsum/int_snow_limited) - 1._r8)/rpi)**(n_melt(c))
                         
                         ddz3 = ddz3 - max(0._r8,(fsno_melt - frac_sno(c))/frac_sno(c))
                      endif
                      ddz3 = -1._r8/dtime * ddz3
                   else
                      ddz3 = - 1._r8/dtime * max(0._r8,(frac_iceold(c,j) - fi)/frac_iceold(c,j))
                   endif
                else
                   ddz3 = 0._r8
                end if

                ! Compaction occurring due to wind drift
                if (wind_dependent_snow_density) then
                   call WindDriftCompaction( &
                        bi = bi, &
                        forc_wind = forc_wind(g), &
                        dz = dz(c,j), &
                        zpseudo = zpseudo(c), &
                        mobile = mobile(c), &
                        compaction_rate = ddz4)
                else
                   ddz4 = 0.0_r8
                end if

                ! Time rate of fractional change in dz (units of s-1)
                pdzdtc = ddz1 + ddz2 + ddz3 + ddz4

                ! The change in dz due to compaction
                ! Limit compaction to be no greater than fully saturated layer thickness
                dz(c,j) = max(dz(c,j) * (1._r8+pdzdtc*dtime),(h2osoi_ice(c,j)/denice+ h2osoi_liq(c,j)/denh2o)/frac_sno(c))

             else
                ! saturated node is immobile
                !
                ! This is only needed if wind_dependent_snow_density is true, but it's
                ! simplest just to update mobile always
                mobile(c) = .false.
             end if

             ! Pressure of overlying snow

             burden(c) = burden(c) + wx

          end if
       end do
    end do

    end associate
  end subroutine SnowCompaction

  !-----------------------------------------------------------------------
  subroutine CombineSnowLayers(bounds, num_snowc, filter_snowc, &
        aerosol_inst, temperature_inst, waterfluxbulk_inst, waterstatebulk_inst, waterdiagnosticbulk_inst)
    !
    ! !DESCRIPTION:
    ! Combine snow layers that are less than a minimum thickness or mass
    ! If the snow element thickness or mass is less than a prescribed minimum,
    ! then it is combined with a neighboring element.  The subroutine
    ! clm\_combo.f90 then executes the combination of mass and energy.
    !
    ! !USES:
    use LakeCon          , only : lsadz
    !
    ! !ARGUMENTS:
    type(bounds_type)      , intent(in)    :: bounds
    integer                , intent(inout) :: num_snowc       ! number of column snow points in column filter
    integer                , intent(inout) :: filter_snowc(:) ! column filter for snow points
    type(aerosol_type)     , intent(inout) :: aerosol_inst
    type(temperature_type) , intent(inout) :: temperature_inst
    type(waterfluxbulk_type)   , intent(inout) :: waterfluxbulk_inst
    type(waterstatebulk_type)  , intent(inout) :: waterstatebulk_inst
    type(waterdiagnosticbulk_type)  , intent(inout) :: waterdiagnosticbulk_inst
    !
    ! !LOCAL VARIABLES:
    integer :: c, fc                            ! column indices
    integer :: i,k                              ! loop indices
    integer :: j,l                              ! node indices
    integer :: msn_old(bounds%begc:bounds%endc) ! number of top snow layer
    integer :: mssi(bounds%begc:bounds%endc)    ! node index
    integer :: neibor                           ! adjacent node selected for combination
    real(r8):: zwice(bounds%begc:bounds%endc)   ! total ice mass in snow
    real(r8):: zwliq (bounds%begc:bounds%endc)  ! total liquid water in snow
    real(r8):: dzminloc(size(dzmin))            ! minimum of top snow layer (local)
    real(r8):: dtime                            !land model time step (sec)

    !-----------------------------------------------------------------------

    associate( &
         ltype            => lun%itype                           , & ! Input:  [integer  (:)   ] landunit type
         urbpoi           => lun%urbpoi                          , & ! Input:  [logical  (:)   ] true => landunit is an urban point

         t_soisno         => temperature_inst%t_soisno_col       , & ! Output: [real(r8) (:,:) ] soil temperature (Kelvin)

         mss_bcphi        => aerosol_inst%mss_bcphi_col          , & ! Output: [real(r8) (:,:) ] hydrophilic BC mass in snow (col,lyr) [kg]
         mss_bcpho        => aerosol_inst%mss_bcpho_col          , & ! Output: [real(r8) (:,:) ] hydrophobic BC mass in snow (col,lyr) [kg]
         mss_ocphi        => aerosol_inst%mss_ocphi_col          , & ! Output: [real(r8) (:,:) ] hydrophilic OC mass in snow (col,lyr) [kg]
         mss_ocpho        => aerosol_inst%mss_ocpho_col          , & ! Output: [real(r8) (:,:) ] hydrophobic OC mass in snow (col,lyr) [kg]
         mss_dst1         => aerosol_inst%mss_dst1_col           , & ! Output: [real(r8) (:,:) ] dust species 1 mass in snow (col,lyr) [kg]
         mss_dst2         => aerosol_inst%mss_dst2_col           , & ! Output: [real(r8) (:,:) ] dust species 2 mass in snow (col,lyr) [kg]
         mss_dst3         => aerosol_inst%mss_dst3_col           , & ! Output: [real(r8) (:,:) ] dust species 3 mass in snow (col,lyr) [kg]
         mss_dst4         => aerosol_inst%mss_dst4_col           , & ! Output: [real(r8) (:,:) ] dust species 4 mass in snow (col,lyr) [kg]

         frac_sno         => waterdiagnosticbulk_inst%frac_sno_col        , & ! Input:  [real(r8) (:)   ] fraction of ground covered by snow (0 to 1)
         frac_sno_eff     => waterdiagnosticbulk_inst%frac_sno_eff_col    , & ! Input:  [real(r8) (:)   ] fraction of ground covered by snow (0 to 1)
         snow_depth       => waterdiagnosticbulk_inst%snow_depth_col      , & ! Output: [real(r8) (:)   ] snow height (m)
         int_snow         => waterstatebulk_inst%int_snow_col        , & ! Output:  [real(r8) (:)   ] integrated snowfall [mm]
         h2osno           => waterstatebulk_inst%h2osno_col          , & ! Output: [real(r8) (:)   ] snow water (mm H2O)
         h2osoi_ice       => waterstatebulk_inst%h2osoi_ice_col      , & ! Output: [real(r8) (:,:) ] ice lens (kg/m2)
         h2osoi_liq       => waterstatebulk_inst%h2osoi_liq_col      , & ! Output: [real(r8) (:,:) ] liquid water (kg/m2)
         snw_rds          => waterdiagnosticbulk_inst%snw_rds_col         , & ! Output: [real(r8) (:,:) ] effective snow grain radius (col,lyr) [microns, m^-6]

         qflx_sl_top_soil => waterfluxbulk_inst%qflx_sl_top_soil_col , & ! Output: [real(r8) (:)   ] liquid water + ice from layer above soil to top soil layer or sent to qflx_qrgwl (mm H2O/s)

         snl              => col%snl                             , & ! Output: [integer  (:)   ] number of snow layers
         dz               => col%dz                              , & ! Output: [real(r8) (:,:) ] layer depth (m)
         zi               => col%zi                              , & ! Output: [real(r8) (:,:) ] interface level below a "z" level (m)
         z                => col%z                                 & ! Output: [real(r8) (:,:) ] layer thickness (m)
    )

    ! Determine model time step

    dtime = get_step_size()

    ! Check the mass of ice lens of snow, when the total is less than a small value,
    ! combine it with the underlying neighbor.

    dzminloc(:) = dzmin(:) ! dzmin will stay constant between timesteps

    ! Add lsadz to dzmin for lakes
    ! Determine whether called from LakeHydrology
    ! Note: this assumes that this function is called separately with the lake-snow and non-lake-snow filters.
    if (num_snowc > 0) then
       c = filter_snowc(1)
       l = col%landunit(c)
       if (ltype(l) == istdlak) then ! Called from LakeHydrology
          dzminloc(:) = dzmin(:) + lsadz
       end if
    end if

    do fc = 1, num_snowc
       c = filter_snowc(fc)

       msn_old(c) = snl(c)
       qflx_sl_top_soil(c) = 0._r8
    end do

    ! The following loop is NOT VECTORIZED

    do fc = 1, num_snowc
       c = filter_snowc(fc)
       l = col%landunit(c)
       do j = msn_old(c)+1,0
          ! use 0.01 to avoid runaway ice buildup
          if (h2osoi_ice(c,j) <= .01_r8) then
             if (ltype(l) == istsoil .or. urbpoi(l) .or. ltype(l) == istcrop) then
                h2osoi_liq(c,j+1) = h2osoi_liq(c,j+1) + h2osoi_liq(c,j)
                h2osoi_ice(c,j+1) = h2osoi_ice(c,j+1) + h2osoi_ice(c,j)

                if (j == 0) then
                   qflx_sl_top_soil(c) = (h2osoi_liq(c,j) + h2osoi_ice(c,j))/dtime
                end if

                if (j /= 0) dz(c,j+1) = dz(c,j+1) + dz(c,j)

                ! NOTE: Temperature, and similarly snw_rds, of the
                ! underlying snow layer are NOT adjusted in this case.
                ! Because the layer being eliminated has a small mass,
                ! this should not make a large difference, but it
                ! would be more thorough to do so.
                if (j /= 0) then
                   mss_bcphi(c,j+1) = mss_bcphi(c,j+1)  + mss_bcphi(c,j)
                   mss_bcpho(c,j+1) = mss_bcpho(c,j+1)  + mss_bcpho(c,j)
                   mss_ocphi(c,j+1) = mss_ocphi(c,j+1)  + mss_ocphi(c,j)
                   mss_ocpho(c,j+1) = mss_ocpho(c,j+1)  + mss_ocpho(c,j)
                   mss_dst1(c,j+1)  = mss_dst1(c,j+1)   + mss_dst1(c,j)
                   mss_dst2(c,j+1)  = mss_dst2(c,j+1)   + mss_dst2(c,j)
                   mss_dst3(c,j+1)  = mss_dst3(c,j+1)   + mss_dst3(c,j)
                   mss_dst4(c,j+1)  = mss_dst4(c,j+1)   + mss_dst4(c,j)
                end if

             else if (ltype(l) /= istsoil .and. .not. urbpoi(l) .and. ltype(l) /= istcrop .and. j /= 0) then

                h2osoi_liq(c,j+1) = h2osoi_liq(c,j+1) + h2osoi_liq(c,j)
                h2osoi_ice(c,j+1) = h2osoi_ice(c,j+1) + h2osoi_ice(c,j)
                dz(c,j+1) = dz(c,j+1) + dz(c,j)

                mss_bcphi(c,j+1) = mss_bcphi(c,j+1)  + mss_bcphi(c,j)
                mss_bcpho(c,j+1) = mss_bcpho(c,j+1)  + mss_bcpho(c,j)
                mss_ocphi(c,j+1) = mss_ocphi(c,j+1)  + mss_ocphi(c,j)
                mss_ocpho(c,j+1) = mss_ocpho(c,j+1)  + mss_ocpho(c,j)
                mss_dst1(c,j+1)  = mss_dst1(c,j+1)   + mss_dst1(c,j)
                mss_dst2(c,j+1)  = mss_dst2(c,j+1)   + mss_dst2(c,j)
                mss_dst3(c,j+1)  = mss_dst3(c,j+1)   + mss_dst3(c,j)
                mss_dst4(c,j+1)  = mss_dst4(c,j+1)   + mss_dst4(c,j)

             end if

             ! shift all elements above this down one.
             if (j > snl(c)+1 .and. snl(c) < -1) then
                do i = j, snl(c)+2, -1
                   ! If the layer closest to the surface is less than 0.1 mm and the ltype is not
                   ! urban, soil or crop, the h2osoi_liq and h2osoi_ice associated with this layer is sent
                   ! to qflx_qrgwl later on in the code.  To keep track of this for the snow balance
                   ! error check, we add this to qflx_sl_top_soil here
                   if (ltype(l) /= istsoil .and. ltype(l) /= istcrop .and. .not. urbpoi(l) .and. i == 0) then
                      qflx_sl_top_soil(c) = (h2osoi_liq(c,i) + h2osoi_ice(c,i))/dtime
                   end if

                   t_soisno(c,i)   = t_soisno(c,i-1)
                   h2osoi_liq(c,i) = h2osoi_liq(c,i-1)
                   h2osoi_ice(c,i) = h2osoi_ice(c,i-1)

                   mss_bcphi(c,i)   = mss_bcphi(c,i-1)
                   mss_bcpho(c,i)   = mss_bcpho(c,i-1)
                   mss_ocphi(c,i)   = mss_ocphi(c,i-1)
                   mss_ocpho(c,i)   = mss_ocpho(c,i-1)
                   mss_dst1(c,i)    = mss_dst1(c,i-1)
                   mss_dst2(c,i)    = mss_dst2(c,i-1)
                   mss_dst3(c,i)    = mss_dst3(c,i-1)
                   mss_dst4(c,i)    = mss_dst4(c,i-1)
                   snw_rds(c,i)     = snw_rds(c,i-1)

                   dz(c,i)         = dz(c,i-1)
                end do
             end if
             snl(c) = snl(c) + 1
          end if
       end do
    end do

    do fc = 1, num_snowc
       c = filter_snowc(fc)
       h2osno(c) = 0._r8
       snow_depth(c) = 0._r8
       zwice(c)  = 0._r8
       zwliq(c)  = 0._r8
    end do

    do j = -nlevsno+1,0
       do fc = 1, num_snowc
          c = filter_snowc(fc)
          if (j >= snl(c)+1) then
             h2osno(c) = h2osno(c) + h2osoi_ice(c,j) + h2osoi_liq(c,j)
             snow_depth(c) = snow_depth(c) + dz(c,j)
             zwice(c)  = zwice(c) + h2osoi_ice(c,j)
             zwliq(c)  = zwliq(c) + h2osoi_liq(c,j)
          end if
       end do
    end do

    ! Check the snow depth - all snow gone
    ! The liquid water assumes ponding on soil surface.

    do fc = 1, num_snowc
       c = filter_snowc(fc)
       l = col%landunit(c)
       if (snow_depth(c) > 0._r8) then
          if ((ltype(l) == istdlak .and. snow_depth(c) < 0.01_r8 + lsadz ) .or. &
               ((ltype(l) /= istdlak) .and. ((frac_sno_eff(c)*snow_depth(c) < 0.01_r8)  &
               .or. (h2osno(c)/(frac_sno_eff(c)*snow_depth(c)) < 50._r8)))) then

             snl(c) = 0
             h2osno(c) = zwice(c)

             mss_bcphi(c,:) = 0._r8
             mss_bcpho(c,:) = 0._r8
             mss_ocphi(c,:) = 0._r8
             mss_ocpho(c,:) = 0._r8
             mss_dst1(c,:)  = 0._r8
             mss_dst2(c,:)  = 0._r8
             mss_dst3(c,:)  = 0._r8
             mss_dst4(c,:)  = 0._r8

             if (h2osno(c) <= 0._r8) snow_depth(c) = 0._r8
             ! this is where water is transfered from layer 0 (snow) to layer 1 (soil)
             if (ltype(l) == istsoil .or. urbpoi(l) .or. ltype(l) == istcrop) then
                h2osoi_liq(c,0) = 0.0_r8
                h2osoi_liq(c,1) = h2osoi_liq(c,1) + zwliq(c)
             end if
             if (ltype(l) == istwet) then
                h2osoi_liq(c,0) = 0.0_r8
             endif
             if (ltype(l)==istice_mec) then
                h2osoi_liq(c,0) = 0.0_r8
             endif
          endif
       end if
       if (h2osno(c) <= 0._r8) then
          snow_depth(c) = 0._r8
          frac_sno(c) = 0._r8
          frac_sno_eff(c) = 0._r8
          int_snow(c) = 0._r8
       endif
    end do

    ! Check the snow depth - snow layers combined
    ! The following loop IS NOT VECTORIZED

    do fc = 1, num_snowc
       c = filter_snowc(fc)

       ! Two or more layers

       if (snl(c) < -1) then

          msn_old(c) = snl(c)
          mssi(c) = 1

          do i = msn_old(c)+1,0
             if ((frac_sno_eff(c)*dz(c,i) < dzminloc(mssi(c))) .or. &
                  ((h2osoi_ice(c,i) + h2osoi_liq(c,i))/(frac_sno_eff(c)*dz(c,i)) < 50._r8)) then
                if (i == snl(c)+1) then
                   ! If top node is removed, combine with bottom neighbor.
                   neibor = i + 1
                else if (i == 0) then
                   ! If the bottom neighbor is not snow, combine with the top neighbor.
                   neibor = i - 1
                else
                   ! If none of the above special cases apply, combine with the thinnest neighbor
                   neibor = i + 1
                   if ((dz(c,i-1)+dz(c,i)) < (dz(c,i+1)+dz(c,i))) neibor = i-1
                end if

                ! Node l and j are combined and stored as node j.
                if (neibor > i) then
                   j = neibor
                   l = i
                else
                   j = i
                   l = neibor
                end if

                ! this should be included in 'Combo' for consistency,
                ! but functionally it is the same to do it here
                mss_bcphi(c,j)=mss_bcphi(c,j)+mss_bcphi(c,l)
                mss_bcpho(c,j)=mss_bcpho(c,j)+mss_bcpho(c,l)
                mss_ocphi(c,j)=mss_ocphi(c,j)+mss_ocphi(c,l)
                mss_ocpho(c,j)=mss_ocpho(c,j)+mss_ocpho(c,l)
                mss_dst1(c,j)=mss_dst1(c,j)+mss_dst1(c,l)
                mss_dst2(c,j)=mss_dst2(c,j)+mss_dst2(c,l)
                mss_dst3(c,j)=mss_dst3(c,j)+mss_dst3(c,l)
                mss_dst4(c,j)=mss_dst4(c,j)+mss_dst4(c,l)

                ! mass-weighted combination of effective grain size:
                snw_rds(c,j) = (snw_rds(c,j)*(h2osoi_liq(c,j)+h2osoi_ice(c,j)) + &
                     snw_rds(c,l)*(h2osoi_liq(c,l)+h2osoi_ice(c,l))) / &
                     (h2osoi_liq(c,j)+h2osoi_ice(c,j)+h2osoi_liq(c,l)+h2osoi_ice(c,l))

                call Combo (dz(c,j), h2osoi_liq(c,j), h2osoi_ice(c,j), &
                     t_soisno(c,j), dz(c,l), h2osoi_liq(c,l), h2osoi_ice(c,l), t_soisno(c,l) )

                ! Now shift all elements above this down one.
                if (j-1 > snl(c)+1) then

                   do k = j-1, snl(c)+2, -1
                      t_soisno(c,k) = t_soisno(c,k-1)
                      h2osoi_ice(c,k) = h2osoi_ice(c,k-1)
                      h2osoi_liq(c,k) = h2osoi_liq(c,k-1)

                      mss_bcphi(c,k) = mss_bcphi(c,k-1)
                      mss_bcpho(c,k) = mss_bcpho(c,k-1)
                      mss_ocphi(c,k) = mss_ocphi(c,k-1)
                      mss_ocpho(c,k) = mss_ocpho(c,k-1)
                      mss_dst1(c,k)  = mss_dst1(c,k-1)
                      mss_dst2(c,k)  = mss_dst2(c,k-1)
                      mss_dst3(c,k)  = mss_dst3(c,k-1)
                      mss_dst4(c,k)  = mss_dst4(c,k-1)
                      snw_rds(c,k)   = snw_rds(c,k-1)

                      dz(c,k) = dz(c,k-1)
                   end do
                end if

                ! Decrease the number of snow layers
                snl(c) = snl(c) + 1
                if (snl(c) >= -1) EXIT

             else

                ! The layer thickness is greater than the prescribed minimum value
                mssi(c) = mssi(c) + 1

             end if
          end do

       end if

    end do

    ! Reset the node depth and the depth of layer interface

    do j = 0, -nlevsno+1, -1
       do fc = 1, num_snowc
          c = filter_snowc(fc)
          if (j >= snl(c) + 1) then
             z(c,j) = zi(c,j) - 0.5_r8*dz(c,j)
             zi(c,j-1) = zi(c,j) - dz(c,j)
          end if
       end do
    end do

    end associate
  end subroutine CombineSnowLayers

  !-----------------------------------------------------------------------
  subroutine DivideSnowLayers(bounds, num_snowc, filter_snowc, &
        aerosol_inst, temperature_inst, waterstatebulk_inst, waterdiagnosticbulk_inst, is_lake)
    !
    ! !DESCRIPTION:
    ! Subdivides snow layers if they exceed their prescribed maximum thickness.
    !
    ! !USES:
    use clm_varcon,  only : tfrz
    use LakeCon   ,  only : lsadz
    !
    ! !ARGUMENTS:
    type(bounds_type)      , intent(in)    :: bounds
    integer                , intent(in)    :: num_snowc       ! number of column snow points in column filter
    integer                , intent(in)    :: filter_snowc(:) ! column filter for snow points
    type(aerosol_type)     , intent(inout) :: aerosol_inst
    type(temperature_type) , intent(inout) :: temperature_inst
    type(waterstatebulk_type)  , intent(inout) :: waterstatebulk_inst
    type(waterdiagnosticbulk_type)  , intent(inout) :: waterdiagnosticbulk_inst
    logical                , intent(in)    :: is_lake  !TODO - this should be examined and removed in the future
    !
    ! !LOCAL VARIABLES:
    integer  :: j, c, fc, k                              ! indices
    real(r8) :: drr                                      ! thickness of the combined [m]
    integer  :: msno                                     ! number of snow layer 1 (top) to msno (bottom)
    real(r8) :: dzsno(bounds%begc:bounds%endc,nlevsno)   ! Snow layer thickness [m]
    real(r8) :: swice(bounds%begc:bounds%endc,nlevsno)   ! Partial volume of ice [m3/m3]
    real(r8) :: swliq(bounds%begc:bounds%endc,nlevsno)   ! Partial volume of liquid water [m3/m3]
    real(r8) :: tsno(bounds%begc:bounds%endc ,nlevsno)   ! Nodel temperature [K]
    real(r8) :: zwice                                    ! temporary
    real(r8) :: zwliq                                    ! temporary
    real(r8) :: propor                                   ! temporary
    real(r8) :: dtdz                                     ! temporary
    ! temporary variables mimicking the structure of other layer division variables
    real(r8) :: mbc_phi(bounds%begc:bounds%endc,nlevsno) ! mass of BC in each snow layer
    real(r8) :: zmbc_phi                                 ! temporary
    real(r8) :: mbc_pho(bounds%begc:bounds%endc,nlevsno) ! mass of BC in each snow layer
    real(r8) :: zmbc_pho                                 ! temporary
    real(r8) :: moc_phi(bounds%begc:bounds%endc,nlevsno) ! mass of OC in each snow layer
    real(r8) :: zmoc_phi                                 ! temporary
    real(r8) :: moc_pho(bounds%begc:bounds%endc,nlevsno) ! mass of OC in each snow layer
    real(r8) :: zmoc_pho                                 ! temporary
    real(r8) :: mdst1(bounds%begc:bounds%endc,nlevsno)   ! mass of dust 1 in each snow layer
    real(r8) :: zmdst1                                   ! temporary
    real(r8) :: mdst2(bounds%begc:bounds%endc,nlevsno)   ! mass of dust 2 in each snow layer
    real(r8) :: zmdst2                                   ! temporary
    real(r8) :: mdst3(bounds%begc:bounds%endc,nlevsno)   ! mass of dust 3 in each snow layer
    real(r8) :: zmdst3                                   ! temporary
    real(r8) :: mdst4(bounds%begc:bounds%endc,nlevsno)   ! mass of dust 4 in each snow layer
    real(r8) :: zmdst4                                   ! temporary
    real(r8) :: rds(bounds%begc:bounds%endc,nlevsno)
    ! Variables for consistency check
    real(r8) :: dztot(bounds%begc:bounds%endc)
    real(r8) :: snwicetot(bounds%begc:bounds%endc)
    real(r8) :: snwliqtot(bounds%begc:bounds%endc)
    real(r8) :: offset ! temporary
    !-----------------------------------------------------------------------

    associate( &
         t_soisno   => temperature_inst%t_soisno_col    , & ! Output: [real(r8) (:,:) ] soil temperature (Kelvin)

         h2osoi_ice => waterstatebulk_inst%h2osoi_ice_col   , & ! Output: [real(r8) (:,:) ] ice lens (kg/m2)
         h2osoi_liq => waterstatebulk_inst%h2osoi_liq_col   , & ! Output: [real(r8) (:,:) ] liquid water (kg/m2)
         frac_sno   => waterdiagnosticbulk_inst%frac_sno_eff_col , & ! Output: [real(r8) (:)   ] fraction of ground covered by snow (0 to 1)
         snw_rds    => waterdiagnosticbulk_inst%snw_rds_col      , & ! Output: [real(r8) (:,:) ] effective snow grain radius (col,lyr) [microns, m^-6]

         mss_bcphi  => aerosol_inst%mss_bcphi_col       , & ! Output: [real(r8) (:,:) ] hydrophilic BC mass in snow (col,lyr) [kg]
         mss_bcpho  => aerosol_inst%mss_bcpho_col       , & ! Output: [real(r8) (:,:) ] hydrophobic BC mass in snow (col,lyr) [kg]
         mss_ocphi  => aerosol_inst%mss_ocphi_col       , & ! Output: [real(r8) (:,:) ] hydrophilic OC mass in snow (col,lyr) [kg]
         mss_ocpho  => aerosol_inst%mss_ocpho_col       , & ! Output: [real(r8) (:,:) ] hydrophobic OC mass in snow (col,lyr) [kg]
         mss_dst1   => aerosol_inst%mss_dst1_col        , & ! Output: [real(r8) (:,:) ] dust species 1 mass in snow (col,lyr) [kg]
         mss_dst2   => aerosol_inst%mss_dst2_col        , & ! Output: [real(r8) (:,:) ] dust species 2 mass in snow (col,lyr) [kg]
         mss_dst3   => aerosol_inst%mss_dst3_col        , & ! Output: [real(r8) (:,:) ] dust species 3 mass in snow (col,lyr) [kg]
         mss_dst4   => aerosol_inst%mss_dst4_col        , & ! Output: [real(r8) (:,:) ] dust species 4 mass in snow (col,lyr) [kg]

         snl        => col%snl                          , & ! Output: [integer  (:)   ] number of snow layers
         dz         => col%dz                           , & ! Output: [real(r8) (:,:) ] layer depth (m)
         zi         => col%zi                           , & ! Output: [real(r8) (:,:) ] interface level below a "z" level (m)
         z          => col%z                              & ! Output: [real(r8) (:,:) ] layer thickness (m)
    )

    if ( is_lake ) then
       ! Initialize for consistency check
       do j = -nlevsno+1,0
          do fc = 1, num_snowc
             c = filter_snowc(fc)

             if (j == -nlevsno+1) then
                dztot(c) = 0._r8
                snwicetot(c) = 0._r8
                snwliqtot(c) = 0._r8
             end if

             if (j >= snl(c)+1) then
                dztot(c) = dztot(c) + dz(c,j)
                snwicetot(c) = snwicetot(c) + h2osoi_ice(c,j)
                snwliqtot(c) = snwliqtot(c) + h2osoi_liq(c,j)
             end if
          end do
       end do
    end if

    ! Begin calculation - note that the following column loops are only invoked
    ! for snow-covered columns

    do j = 1,nlevsno
       do fc = 1, num_snowc
          c = filter_snowc(fc)
          if (j <= abs(snl(c))) then
             if (is_lake) then
                dzsno(c,j) = dz(c,j+snl(c))
             else
                dzsno(c,j) = frac_sno(c)*dz(c,j+snl(c))
             end if
             swice(c,j) = h2osoi_ice(c,j+snl(c))
             swliq(c,j) = h2osoi_liq(c,j+snl(c))
             tsno(c,j)  = t_soisno(c,j+snl(c))

             mbc_phi(c,j) = mss_bcphi(c,j+snl(c))
             mbc_pho(c,j) = mss_bcpho(c,j+snl(c))
             moc_phi(c,j) = mss_ocphi(c,j+snl(c))
             moc_pho(c,j) = mss_ocpho(c,j+snl(c))
             mdst1(c,j)   = mss_dst1(c,j+snl(c))
             mdst2(c,j)   = mss_dst2(c,j+snl(c))
             mdst3(c,j)   = mss_dst3(c,j+snl(c))
             mdst4(c,j)   = mss_dst4(c,j+snl(c))
             rds(c,j)     = snw_rds(c,j+snl(c))
          end if
       end do
    end do

    loop_snowcolumns: do fc = 1, num_snowc
       c = filter_snowc(fc)

       msno = abs(snl(c))

       ! Now traverse layers from top to bottom in a dynamic way, as the total
       ! number of layers (msno) may increase during the loop.
       ! Impose k < nlevsno; the special case 'k == nlevsno' is not relevant,
       ! as it is neither allowed to subdivide nor does it have layers below.
       k = 1
       loop_layers: do while( k <= msno .and. k < nlevsno )

          ! Current layer is bottom layer
          if (k == msno) then

             if (is_lake) then
                offset = 2._r8 * lsadz
             else
                offset = 0._r8
             end if

             if (dzsno(c,k) > dzmax_l(k) + offset) then
                ! Subdivide layer into two layers with equal thickness, water
                ! content, ice content and temperature
                msno = msno + 1
                dzsno(c,k)     = dzsno(c,k) / 2.0_r8
                dzsno(c,k+1)   = dzsno(c,k)
                swice(c,k)     = swice(c,k) / 2.0_r8
                swice(c,k+1)   = swice(c,k)
                swliq(c,k)     = swliq(c,k) / 2.0_r8
                swliq(c,k+1)   = swliq(c,k)

                if (k == 1) then
                   ! special case
                   tsno(c,k+1)    = tsno(c,k)
                else
                   ! use temperature gradient
                   dtdz           = (tsno(c,k-1) - tsno(c,k))/((dzsno(c,k-1)+2*dzsno(c,k))/2.0_r8)
                   tsno(c,k+1) = tsno(c,k) - dtdz*dzsno(c,k)/2.0_r8
                   if (tsno(c,k+1) >= tfrz) then
                      tsno(c,k+1)  = tsno(c,k)
                   else
                      tsno(c,k) = tsno(c,k) + dtdz*dzsno(c,k)/2.0_r8
                   endif
                end if

                mbc_phi(c,k)   = mbc_phi(c,k) / 2.0_r8
                mbc_phi(c,k+1) = mbc_phi(c,k)
                mbc_pho(c,k)   = mbc_pho(c,k) / 2.0_r8
                mbc_pho(c,k+1) = mbc_pho(c,k)
                moc_phi(c,k)   = moc_phi(c,k) / 2.0_r8
                moc_phi(c,k+1) = moc_phi(c,k)
                moc_pho(c,k)   = moc_pho(c,k) / 2.0_r8
                moc_pho(c,k+1) = moc_pho(c,k)
                mdst1(c,k)     = mdst1(c,k) / 2.0_r8
                mdst1(c,k+1)   = mdst1(c,k)
                mdst2(c,k)     = mdst2(c,k) / 2.0_r8
                mdst2(c,k+1)   = mdst2(c,k)
                mdst3(c,k)     = mdst3(c,k) / 2.0_r8
                mdst3(c,k+1)   = mdst3(c,k)
                mdst4(c,k)     = mdst4(c,k) / 2.0_r8
                mdst4(c,k+1)   = mdst4(c,k)

                rds(c,k+1)     = rds(c,k)
             end if
          end if

          ! There are layers below (note this is not exclusive with previous
          ! if-statement, since msno may have increased in the previous if-statement)
          if (k < msno) then

             if (is_lake) then
                offset = lsadz
             else
                offset = 0._r8
             end if

             if (dzsno(c,k) > dzmax_u(k) + offset ) then
                ! Only dump excess snow to underlying layer in a conservative fashion.
                ! Other quantities will depend on the height of the excess snow: a ratio is used for this.
                drr      = dzsno(c,k) - dzmax_u(k) - offset

                propor   = drr/dzsno(c,k)
                zwice    = propor*swice(c,k)
                zwliq    = propor*swliq(c,k)
                zmbc_phi = propor*mbc_phi(c,k)
                zmbc_pho = propor*mbc_pho(c,k)
                zmoc_phi = propor*moc_phi(c,k)
                zmoc_pho = propor*moc_pho(c,k)
                zmdst1   = propor*mdst1(c,k)
                zmdst2   = propor*mdst2(c,k)
                zmdst3   = propor*mdst3(c,k)
                zmdst4   = propor*mdst4(c,k)

                propor         = (dzmax_u(k)+offset)/dzsno(c,k)
                swice(c,k)     = propor*swice(c,k)
                swliq(c,k)     = propor*swliq(c,k)
                mbc_phi(c,k)   = propor*mbc_phi(c,k)
                mbc_pho(c,k)   = propor*mbc_pho(c,k)
                moc_phi(c,k)   = propor*moc_phi(c,k)
                moc_pho(c,k)   = propor*moc_pho(c,k)
                mdst1(c,k)     = propor*mdst1(c,k)
                mdst2(c,k)     = propor*mdst2(c,k)
                mdst3(c,k)     = propor*mdst3(c,k)
                mdst4(c,k)     = propor*mdst4(c,k)

                ! Set depth layer k to maximum allowed value
                dzsno(c,k)  = dzmax_u(k)  + offset

                mbc_phi(c,k+1) = mbc_phi(c,k+1)+zmbc_phi  ! (combo)
                mbc_pho(c,k+1) = mbc_pho(c,k+1)+zmbc_pho  ! (combo)
                moc_phi(c,k+1) = moc_phi(c,k+1)+zmoc_phi  ! (combo)
                moc_pho(c,k+1) = moc_pho(c,k+1)+zmoc_pho  ! (combo)
                mdst1(c,k+1)   = mdst1(c,k+1)+zmdst1  ! (combo)
                mdst2(c,k+1)   = mdst2(c,k+1)+zmdst2  ! (combo)
                mdst3(c,k+1)   = mdst3(c,k+1)+zmdst3  ! (combo)
                mdst4(c,k+1)   = mdst4(c,k+1)+zmdst4  ! (combo)

                ! Mass-weighted combination of radius
                rds(c,k+1) = MassWeightedSnowRadius( rds(c,k), rds(c,k+1), &
                     (swliq(c,k+1)+swice(c,k+1)), (zwliq+zwice) )

                call Combo (dzsno(c,k+1), swliq(c,k+1), swice(c,k+1), tsno(c,k+1), drr, &
                     zwliq, zwice, tsno(c,k))
             end if
          end if
          k = k+1
       end do loop_layers

       snl(c) = -msno

    end do loop_snowcolumns

    do j = -nlevsno+1,0
       do fc = 1, num_snowc
          c = filter_snowc(fc)
          if (j >= snl(c)+1) then
             if (is_lake) then
                dz(c,j) = dzsno(c,j-snl(c))
             else
                dz(c,j) = dzsno(c,j-snl(c))/frac_sno(c)
             end if
             h2osoi_ice(c,j) = swice(c,j-snl(c))
             h2osoi_liq(c,j) = swliq(c,j-snl(c))
             t_soisno(c,j)   = tsno(c,j-snl(c))
             mss_bcphi(c,j)   = mbc_phi(c,j-snl(c))
             mss_bcpho(c,j)   = mbc_pho(c,j-snl(c))
             mss_ocphi(c,j)   = moc_phi(c,j-snl(c))
             mss_ocpho(c,j)   = moc_pho(c,j-snl(c))
             mss_dst1(c,j)    = mdst1(c,j-snl(c))
             mss_dst2(c,j)    = mdst2(c,j-snl(c))
             mss_dst3(c,j)    = mdst3(c,j-snl(c))
             mss_dst4(c,j)    = mdst4(c,j-snl(c))
             snw_rds(c,j)     = rds(c,j-snl(c))

          end if
       end do
    end do

    ! Consistency check
    if (is_lake) then
       do j = -nlevsno + 1, 0
          do fc = 1, num_snowc
             c = filter_snowc(fc)

             if (j >= snl(c)+1) then
                dztot(c) = dztot(c) - dz(c,j)
                snwicetot(c) = snwicetot(c) - h2osoi_ice(c,j)
                snwliqtot(c) = snwliqtot(c) - h2osoi_liq(c,j)
             end if

             if (j == 0) then
                if ( abs(dztot(c)) > 1.e-10_r8 .or. abs(snwicetot(c)) > 1.e-7_r8 .or. &
                     abs(snwliqtot(c)) > 1.e-7_r8 ) then
                   write(iulog,*)'Inconsistency in SnowDivision_Lake! c, remainders', &
                        'dztot, snwicetot, snwliqtot = ',c,dztot(c),snwicetot(c),snwliqtot(c)
                   call endrun(decomp_index=c, clmlevel=namec, msg=errmsg(sourcefile, __LINE__))
                end if
             end if
          end do
       end do
    end if

    ! Reset the node depth and the depth of layer interface

    do j = 0, -nlevsno+1, -1
       do fc = 1, num_snowc
          c = filter_snowc(fc)
          if (j >= snl(c)+1) then
             z(c,j)    = zi(c,j) - 0.5_r8*dz(c,j)
             zi(c,j-1) = zi(c,j) - dz(c,j)
          end if
       end do
    end do

    end associate
  end subroutine DivideSnowLayers

  !-----------------------------------------------------------------------
  subroutine InitSnowLayers (bounds, snow_depth)
    !
    ! !DESCRIPTION:
    ! Initialize snow layer depth from specified total depth.
    !
    ! !USES:
    use clm_varcon         , only : spval
    !
    ! !ARGUMENTS:
    type(bounds_type)      , intent(in)    :: bounds
    real(r8)               , intent(in)    :: snow_depth(bounds%begc:)
    !
    !
    ! LOCAL VARAIBLES:
    integer               :: c,l,j              ! indices
    real(r8)              :: minbound, maxbound ! helper variables
    !------------------------------------------------------------------------

    SHR_ASSERT_ALL((ubound(snow_depth)  == (/bounds%endc/)), errMsg(sourcefile, __LINE__))

    associate( &
         snl => col%snl,   & ! Output: [integer (:)    ]  number of snow layers
         dz  => col%dz,    & ! Output: [real(r8) (:,:) ]  layer thickness (m)  (-nlevsno+1:nlevgrnd)
         z   => col%z,     & ! Output: [real(r8) (:,:) ]  layer depth (m) (-nlevsno+1:nlevgrnd)
         zi  => col%zi     & ! Output: [real(r8) (:,:) ]  interface level below a "z" level (m) (-nlevsno+0:nlevgrnd)
    )

    loop_columns: do c = bounds%begc,bounds%endc
       l = col%landunit(c)

       dz(c,-nlevsno+1: 0) = spval
       z (c,-nlevsno+1: 0) = spval
       zi(c,-nlevsno  :-1) = spval

       ! Special case: lake
       if (lun%lakpoi(l)) then
          snl(c)              = 0
          dz(c,-nlevsno+1:0)  = 0._r8
          z(c,-nlevsno+1:0)   = 0._r8
          zi(c,-nlevsno+0:0)  = 0._r8
          cycle
       end if

       ! LvK 9-JUN-2015: in CanopyHydrologyMod , snow_depth is scaled with frac_sno
       ! Here we do not apply scaling to snow_depth, so inconsistent? TODO

       ! Special case: too little snow for snowpack existence
       if (snow_depth(c) < dzmin(1)) then
          snl(c)              = 0
          dz(c,-nlevsno+1:0)  = 0._r8
          z(c,-nlevsno+1:0)   = 0._r8
          zi(c,-nlevsno+0:0)  = 0._r8
          cycle
       end if

       ! There has to be at least one snow layer
       snl(c)   = -1
       minbound = dzmin(1)
       maxbound = dzmax_l(1)

       if (snow_depth(c) >= minbound .and. snow_depth(c) <= maxbound) then
          ! Special case: single layer
          dz(c,0) = snow_depth(c)

       else
          ! Search for appropriate number of layers (snl) by increasing the number
          ! the number of layers and check for matching bounds.
          snl(c) = snl(c) - 1
          minbound = maxbound
          maxbound = sum(dzmax_u(1:-snl(c)))

          do while(snow_depth(c) > maxbound .and. -snl(c) < nlevsno )
             snl(c) = snl(c) - 1
             minbound = maxbound
             maxbound = sum(dzmax_u(1:-snl(c)))
          end do

          ! Set thickness of all layers except bottom two
          do j = 1, -snl(c)-2
             dz(c,j+snl(c))  = dzmax_u(j)
          enddo

          ! Determine whether the two bottom layers should be equal in size,
          ! or not. The rule here is: always create equal size when possible.
          if (snow_depth(c) <= sum(dzmax_u(1:-snl(c)-2)) + 2 * dzmax_u(-snl(c)-1)) then
             dz(c,-1) = (snow_depth(c) - sum(dzmax_u(1:-snl(c)-2))) / 2._r8
             dz(c,0)  = dz(c,-1)
          else
             dz(c,-1) = dzmax_u(-snl(c)-1)
             dz(c,0)  = snow_depth(c) - sum(dzmax_u(1:-snl(c)-1))
          endif
       endif

       ! Initialize the node depth and the depth of layer interface
       do j = 0, snl(c)+1, -1
          z(c,j)    = zi(c,j) - 0.5_r8*dz(c,j)
          zi(c,j-1) = zi(c,j) - dz(c,j)
       end do

    end do loop_columns

    end associate
  end subroutine InitSnowLayers

  !-----------------------------------------------------------------------
  subroutine SnowCapping(bounds, num_initc, filter_initc, num_snowc, filter_snowc, &
       aerosol_inst, waterfluxbulk_inst, waterstatebulk_inst, topo_inst )
    !
    ! !DESCRIPTION:
    ! Removes mass from bottom snow layer for columns that exceed the maximum snow depth.
    ! This routine is called twice: once for non-lake columns and once for lake columns. 
    ! The initialization of the snow capping fluxes should only be done ONCE for each group,
    ! therefore they are a passed as an extra argument (filter_initc). 
    ! Density and temperature of the layer are conserved (density needs some work, temperature is a state
    ! variable)
    !
    ! !ARGUMENTS:
    type(bounds_type)      , intent(in)    :: bounds          
    integer                , intent(in)    :: num_initc       ! number of column points that need to be initialized
    integer                , intent(in)    :: filter_initc(:) ! column filter for points that need to be initialized
    integer                , intent(in)    :: num_snowc       ! number of column snow points in column filter
    integer                , intent(in)    :: filter_snowc(:) ! column filter for snow points
    type(aerosol_type)     , intent(inout) :: aerosol_inst
    type(waterfluxbulk_type)   , intent(inout) :: waterfluxbulk_inst 
    type(waterstatebulk_type)  , intent(inout) :: waterstatebulk_inst
    class(topo_type)   , intent(in)    :: topo_inst
    !
    ! !LOCAL VARIABLES:
    real(r8)   :: dtime                            ! land model time step (sec)
    real(r8)   :: mss_snwcp_tot                    ! total snow capping mass [kg/m2] 
    real(r8)   :: mss_snow_bottom_lyr              ! total snow mass (ice+liquid) in bottom layer [kg/m2]
    real(r8)   :: snwcp_flux_ice                   ! snow capping flux (ice) [kg/m2]
    real(r8)   :: snwcp_flux_liq                   ! snow capping flux (liquid) [kg/m2]
    real(r8)   :: icefrac                          ! fraction of ice mass w.r.t. total mass [unitless]
    real(r8)   :: frac_adjust                      ! fraction of mass remaining after capping
    real(r8)   :: rho                              ! partial density of ice (not scaled with frac_sno) [kg/m3]
    integer    :: fc, c                            ! counters
    real(r8)   :: h2osno_excess(bounds%begc:bounds%endc) ! excess snow that needs to be capped [mm H2O]
    logical    :: apply_runoff(bounds%begc:bounds%endc)  ! for columns with capping, whether the capping flux should be sent to runoff
    ! Always keep at least this fraction of the bottom snow layer when doing snow capping
    ! This needs to be slightly greater than 0 to avoid roundoff problems
    real(r8), parameter :: min_snow_to_keep = 1.e-9  ! fraction of bottom snow layer to keep with capping

    !-----------------------------------------------------------------------
    associate( &
        qflx_snwcp_ice     => waterfluxbulk_inst%qflx_snwcp_ice_col   , & ! Output: [real(r8) (:)   ]  excess solid h2o due to snow capping (outgoing) (mm H2O /s) [+]
        qflx_snwcp_liq     => waterfluxbulk_inst%qflx_snwcp_liq_col   , & ! Output: [real(r8) (:)   ]  excess liquid h2o due to snow capping (outgoing) (mm H2O /s) [+]
        qflx_snwcp_discarded_ice => waterfluxbulk_inst%qflx_snwcp_discarded_ice_col, & ! Output: [real(r8) (:)   ]  excess solid h2o due to snow capping, which we simply discard in order to reset the snow pack (mm H2O /s) [+]
        qflx_snwcp_discarded_liq => waterfluxbulk_inst%qflx_snwcp_discarded_liq_col, & ! Output: [real(r8) (:)   ]  excess liquid h2o due to snow capping, which we simply discard in order to reset the snow pack (mm H2O /s) [+]
        h2osoi_ice         => waterstatebulk_inst%h2osoi_ice_col      , & ! In/Out: [real(r8) (:,:) ] ice lens (kg/m2)                       
        h2osoi_liq         => waterstatebulk_inst%h2osoi_liq_col      , & ! In/Out: [real(r8) (:,:) ] liquid water (kg/m2)                   
        h2osno             => waterstatebulk_inst%h2osno_col          , & ! Input:  [real(r8) (:)   ] snow water (mm H2O)
        mss_bcphi          => aerosol_inst%mss_bcphi_col          , & ! In/Out: [real(r8) (:,:) ] hydrophilic BC mass in snow (col,lyr) [kg]
        mss_bcpho          => aerosol_inst%mss_bcpho_col          , & ! In/Out: [real(r8) (:,:) ] hydrophobic BC mass in snow (col,lyr) [kg]
        mss_ocphi          => aerosol_inst%mss_ocphi_col          , & ! In/Out: [real(r8) (:,:) ] hydrophilic OC mass in snow (col,lyr) [kg]
        mss_ocpho          => aerosol_inst%mss_ocpho_col          , & ! In/Out: [real(r8) (:,:) ] hydrophobic OC mass in snow (col,lyr) [kg]
        mss_dst1           => aerosol_inst%mss_dst1_col           , & ! In/Out: [real(r8) (:,:) ] dust species 1 mass in snow (col,lyr) [kg]
        mss_dst2           => aerosol_inst%mss_dst2_col           , & ! In/Out: [real(r8) (:,:) ] dust species 2 mass in snow (col,lyr) [kg]
        mss_dst3           => aerosol_inst%mss_dst3_col           , & ! In/Out: [real(r8) (:,:) ] dust species 3 mass in snow (col,lyr) [kg]
        mss_dst4           => aerosol_inst%mss_dst4_col           , & ! In/Out: [real(r8) (:,:) ] dust species 4 mass in snow (col,lyr) [kg]
        topo               => topo_inst%topo_col                  , & ! Input : [real(r8) (:)   ] column surface height (m)
        dz                 => col%dz                                & ! In/Out: [real(r8) (:,:) ] layer depth (m)
    )

    ! Determine model time step
    dtime = get_step_size()

    ! Initialize capping fluxes for all columns in domain (lake or non-lake)
    do fc = 1, num_initc
       c = filter_initc(fc)
       qflx_snwcp_ice(c) = 0.0_r8
       qflx_snwcp_liq(c) = 0.0_r8
       qflx_snwcp_discarded_ice(c) = 0.0_r8
       qflx_snwcp_discarded_liq(c) = 0.0_r8
    end do

    call SnowCappingExcess(bounds, num_snowc, filter_snowc, &
         h2osno = h2osno(bounds%begc:bounds%endc), &
         topo = topo(bounds%begc:bounds%endc), &
         h2osno_excess = h2osno_excess(bounds%begc:bounds%endc), &
         apply_runoff = apply_runoff(bounds%begc:bounds%endc))

    loop_columns: do fc = 1, num_snowc
       c = filter_snowc(fc)

       if (h2osno_excess(c) > 0._r8) then
          mss_snow_bottom_lyr = h2osoi_ice(c,0) + h2osoi_liq(c,0) 
          mss_snwcp_tot = min(h2osno_excess(c), mss_snow_bottom_lyr * (1._r8 - min_snow_to_keep)) ! Can't remove more mass than available

          ! Ratio of snow/liquid in bottom layer determines partitioning of runoff fluxes
          icefrac = h2osoi_ice(c,0) / mss_snow_bottom_lyr
          snwcp_flux_ice = mss_snwcp_tot/dtime * icefrac
          snwcp_flux_liq = mss_snwcp_tot/dtime * (1._r8 - icefrac)
          if (apply_runoff(c)) then
             qflx_snwcp_ice(c) = snwcp_flux_ice
             qflx_snwcp_liq(c) = snwcp_flux_liq
          else
             qflx_snwcp_discarded_ice(c) = snwcp_flux_ice
             qflx_snwcp_discarded_liq(c) = snwcp_flux_liq
          end if

          rho = h2osoi_ice(c,0) / dz(c,0) ! ice only

          ! Adjust water content
          h2osoi_ice(c,0) = h2osoi_ice(c,0) - snwcp_flux_ice*dtime
          h2osoi_liq(c,0) = h2osoi_liq(c,0) - snwcp_flux_liq*dtime

          ! Scale dz such that ice density (or: pore space) is conserved
          !
          ! Avoid scaling dz for very low ice densities. This can occur, in principle, if
          ! the layer is mostly liquid water. Furthermore, this check is critical in the
          ! unlikely event that rho is 0, which can happen if the layer is entirely liquid
          ! water.
          if (rho > 1.0_r8) then
            dz(c,0) = h2osoi_ice(c,0) / rho 
          end if

          ! Check that water capacity is still positive
          if (h2osoi_ice(c,0) < 0._r8 .or. h2osoi_liq(c,0) < 0._r8 ) then
             write(iulog,*)'ERROR: capping procedure failed (negative mass remaining) c = ',c
             write(iulog,*)'h2osoi_ice = ', h2osoi_ice(c,0), ' h2osoi_liq = ', h2osoi_liq(c,0)
             call endrun(decomp_index=c, clmlevel=namec, msg=errmsg(sourcefile, __LINE__))
          end if

          ! Correct the top layer aerosol mass to account for snow capping.
          ! This approach conserves the aerosol mass concentration but not aerosol mass. 
          frac_adjust = (mss_snow_bottom_lyr - mss_snwcp_tot) / mss_snow_bottom_lyr
          mss_bcphi(c,0)   = mss_bcphi(c,0) * frac_adjust 
          mss_bcpho(c,0)   = mss_bcpho(c,0) * frac_adjust
          mss_ocphi(c,0)   = mss_ocphi(c,0) * frac_adjust
          mss_ocpho(c,0)   = mss_ocpho(c,0) * frac_adjust
          mss_dst1(c,0)    = mss_dst1(c,0) * frac_adjust
          mss_dst2(c,0)    = mss_dst2(c,0) * frac_adjust
          mss_dst3(c,0)    = mss_dst3(c,0) * frac_adjust
          mss_dst4(c,0)    = mss_dst4(c,0) * frac_adjust
       end if

    end do loop_columns

    end associate
  end subroutine SnowCapping

  !-----------------------------------------------------------------------
  subroutine SnowCappingExcess(bounds, num_snowc, filter_snowc, &
       h2osno, topo, h2osno_excess, apply_runoff)
    !
    ! !DESCRIPTION:
    ! Determine the amount of excess snow that needs to be capped
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds          
    integer  , intent(in)  :: num_snowc                     ! number of column snow points in column filter
    integer  , intent(in)  :: filter_snowc(:)               ! column filter for snow points
    real(r8) , intent(in)  :: h2osno( bounds%begc: )        ! snow water (mm H2O)
    real(r8) , intent(in)  :: topo( bounds%begc: )          ! column surface height (m)
    real(r8) , intent(out) :: h2osno_excess( bounds%begc: ) ! excess snow that needs to be capped (mm H2O)
    logical  , intent(out) :: apply_runoff( bounds%begc: )  ! whether capped snow should be sent to runoff; only valid where h2osno_excess > 0
    !
    ! !LOCAL VARIABLES:
    integer :: fc, c, l
    integer :: reset_snow_timesteps
    logical :: is_reset_snow_active  ! whether snow resetting is active in this time step for at least some points

    character(len=*), parameter :: subname = 'SnowCappingExcess'
    !-----------------------------------------------------------------------

    SHR_ASSERT_ALL((ubound(h2osno) == (/bounds%endc/)), errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(topo) == (/bounds%endc/)), errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(h2osno_excess) == (/bounds%endc/)), errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(apply_runoff) == (/bounds%endc/)), errMsg(sourcefile, __LINE__))

    do fc = 1, num_snowc
       c = filter_snowc(fc)
       h2osno_excess(c) = 0._r8
       if (h2osno(c) > h2osno_max) then
          h2osno_excess(c) = h2osno(c) - h2osno_max
          apply_runoff(c) = .true.
       end if
    end do

    ! Implement snow resetting (i.e., resetting points that have h2osno greater than some
    ! value) by applying the snow capping scheme to a value that's smaller than
    ! h2osno_max, but NOT sending the resulting capping flux to the coupler. It is easier
    ! to implement the resetting this way than to try to manually reset the snow pack,
    ! because there are so many snow pack variables that need to be kept consistent if
    ! doing this resetting manually. Note that we need to continue to apply the resetting
    ! for some number of time steps, because we can remove at most one snow layer per
    ! time step.

    ! It is important that this snow resetting comes after the standard check (h2osno(c) >
    ! h2osno_max), so that we override any standard capping.
    is_reset_snow_active = .false.
    if (reset_snow .or. reset_snow_glc) then
       reset_snow_timesteps = reset_snow_timesteps_per_layer * nlevsno
       if (get_nstep() <= reset_snow_timesteps) then
          is_reset_snow_active = .true.
       end if
    end if

    if (is_reset_snow_active) then
       do fc = 1, num_snowc
          c = filter_snowc(fc)
          l = col%landunit(c)
          if ((lun%itype(l) /= istice_mec) .and. &
               reset_snow .and. &
               (h2osno(c) > reset_snow_h2osno)) then
             h2osno_excess(c) = h2osno(c) - reset_snow_h2osno
             apply_runoff(c) = .false.
          else if ((lun%itype(l) == istice_mec) .and. &
               reset_snow_glc .and. &
               (h2osno(c) > reset_snow_h2osno) .and. &
               (topo(c) <= reset_snow_glc_ela)) then
             h2osno_excess(c) = h2osno(c) - reset_snow_h2osno
             apply_runoff(c) = .false.
          end if
       end do
    end if

  end subroutine SnowCappingExcess

  !-----------------------------------------------------------------------
  subroutine NewSnowBulkDensity(bounds, num_c, filter_c, atm2lnd_inst, bifall)
    !
    ! !DESCRIPTION:
    ! Compute the bulk density of any newly-fallen snow.
    !
    ! The return value is placed in bifall. Only columns within the given filter are set:
    ! all other columns remain at their original values.
    !
    ! !USES:
    use clm_varcon,  only : tfrz
    !
    ! !ARGUMENTS:
    type(bounds_type)  , intent(in)    :: bounds
    integer            , intent(in)    :: num_c                ! number of columns in filterc
    integer            , intent(in)    :: filter_c(:)          ! column-level filter to operate on
    type(atm2lnd_type) , intent(in)    :: atm2lnd_inst
    real(r8)           , intent(inout) :: bifall(bounds%begc:) ! bulk density of newly fallen dry snow [kg/m3]
    !
    ! !LOCAL VARIABLES:
    integer :: fc, c, g
    real(r8) :: t_for_bifall_degC  ! temperature to use in bifall equation (deg C)

    character(len=*), parameter :: subname = 'NewSnowBulkDensity'
    !-----------------------------------------------------------------------

    SHR_ASSERT_ALL((ubound(bifall) == (/bounds%endc/)), errMsg(sourcefile, __LINE__))

    associate( &
         forc_t      => atm2lnd_inst%forc_t_downscaled_col , & ! Input:  [real(r8) (:)   ]  atmospheric temperature (Kelvin)        
         forc_wind   => atm2lnd_inst%forc_wind_grc           & ! Input:  [real(r8) (:)   ]  atmospheric wind speed (m/s)
         )

    do fc = 1, num_c
       c = filter_c(fc)
       g = col%gridcell(c)

       if (forc_t(c) > tfrz + 2._r8) then
          bifall(c) = 50._r8 + 1.7_r8*(17.0_r8)**1.5_r8
       else if (forc_t(c) > tfrz - 15._r8) then
          bifall(c) = 50._r8 + 1.7_r8*(forc_t(c) - tfrz + 15._r8)**1.5_r8
       else if ( new_snow_density == LoTmpDnsTruncatedAnderson1976 ) then
          bifall(c) = 50._r8
       else if (new_snow_density == LoTmpDnsSlater2017) then 
          ! Andrew Slater: A temp of about -15C gives the nicest
          ! "blower" powder, but as you get colder the flake size decreases so
          ! density goes up. e.g. the smaller snow crystals from the Arctic and Antarctic
          ! winters
          if (forc_t(c) > tfrz - 57.55_r8) then
             t_for_bifall_degC = (forc_t(c)-tfrz)
          else
             ! Below -57.55 deg C, the following function starts to decrease with
             ! decreasing temperatures. Limit the function to avoid this turning over.
             t_for_bifall_degC = -57.55_r8
          end if
          bifall(c) = -(50._r8/15._r8 + 0.0333_r8*15_r8)*t_for_bifall_degC - 0.0333_r8*t_for_bifall_degC**2
       end if

       if (wind_dependent_snow_density .and. forc_wind(g) > 0.1_r8 ) then
       ! Density offset for wind-driven compaction, initial ideas based on Liston et. al (2007) J. Glaciology,
       ! 53(181), 241-255. Modified for a continuous wind impact and slightly more sensitive
       ! to wind - Andrew Slater, 2016
          bifall(c) = bifall(c) + (266.861_r8 * ((1._r8 + TANH(forc_wind(g)/5.0_r8))/2._r8)**8.8_r8)
       end if

    end do

    end associate

  end subroutine NewSnowBulkDensity

  !-----------------------------------------------------------------------
  pure function OverburdenCompactionAnderson1976(burden, wx, td, bi) &
       result(compaction_rate)
    !
    ! !DESCRIPTION:
    ! Compute snow overburden compaction for a single column and level using the Anderson
    ! 1976 formula
    !
    ! From Anderson 1976: A point energy and mass balance model of a snow cover, NOAA
    ! Technical Report NWS 19
    !
    ! !ARGUMENTS:
    real(r8) :: compaction_rate ! function result
    real(r8) , intent(in) :: burden ! pressure of overlying snow in this column [kg/m2]
    real(r8) , intent(in) :: wx     ! water mass (ice+liquid) [kg/m2]
    real(r8) , intent(in) :: td     ! t_soisno - tfrz [K]
    real(r8) , intent(in) :: bi     ! partial density of ice [kg/m3]
    !
    ! !LOCAL VARIABLES:
    real(r8), parameter :: c2 = 23.e-3_r8       ! [m3/kg]
    real(r8), parameter :: eta0 = 9.e+5_r8      ! The Viscosity Coefficient Eta0 [kg-s/m2]

    character(len=*), parameter :: subname = 'OverburdenCompactionAnderson1976'
    !-----------------------------------------------------------------------

    compaction_rate = -(burden+wx/2._r8)*exp(-overburden_compress_Tfactor*td - c2*bi)/eta0

  end function OverburdenCompactionAnderson1976

  !-----------------------------------------------------------------------
  function OverburdenCompactionVionnet2012(h2osoi_liq, dz, burden, wx, td, bi) &
       result(compaction_rate)
    !
    ! !DESCRIPTION:
    ! Compute snow overburden compaction for a single column and level using the Vionnet
    ! et al. 2012 formula
    !
    ! From Vionnet V et al. 2012, "The detailed snowpack scheme Crocus and its
    ! implementation in SURFEX v7.2", Geosci. Model Dev. 5, 773–791.
    !
    ! Preconditions (required to avoid divide by 0):
    ! - dz > 0
    ! - bi > 0
    !
    ! !USES:
    use clm_varcon, only : denh2o
    !
    ! !ARGUMENTS:
    real(r8) :: compaction_rate ! function result
    real(r8) , intent(in) :: h2osoi_liq ! liquid water in this column and level [kg/m2]
    real(r8) , intent(in) :: dz         ! layer depth for this column and level [m]
    real(r8) , intent(in) :: burden     ! pressure of overlying snow in this column [kg/m2]
    real(r8) , intent(in) :: wx         ! water mass (ice+liquid) [kg/m2]
    real(r8) , intent(in) :: td         ! t_soisno - tfrz [K]
    real(r8) , intent(in) :: bi         ! partial density of ice [kg/m3]
    !
    ! !LOCAL VARIABLES:
    real(r8) :: f1, f2                          ! overburden compaction modifiers to viscosity
    real(r8) :: eta                             ! Viscosity

    real(r8), parameter :: ceta = 450._r8       ! overburden compaction constant [kg/m3]
    real(r8), parameter :: aeta = 0.1_r8        ! overburden compaction constant [1/K]
    real(r8), parameter :: beta = 0.023_r8      ! overburden compaction constant [m3/kg]
    real(r8), parameter :: eta0 = 7.62237e6_r8  ! The Viscosity Coefficient Eta0 [kg-s/m2]

    character(len=*), parameter :: subname = 'OverburdenCompactionVionnet2012'
    !-----------------------------------------------------------------------

    f1 = 1._r8 / (1._r8 + 60._r8 * h2osoi_liq / (denh2o * dz))
    f2 = 4.0_r8 ! currently fixed to maximum value, holds in absence of angular grains
    eta = f1*f2*(bi/ceta)*exp(aeta*td + beta*bi)*eta0
    compaction_rate = -(burden+wx/2._r8) / eta

  end function OverburdenCompactionVionnet2012

  !-----------------------------------------------------------------------
  subroutine WindDriftCompaction(bi, forc_wind, dz, &
       zpseudo, mobile, compaction_rate)
    !
    ! !DESCRIPTION:
    !
    ! Compute wind drift compaction for a single column and level.
    !
    ! Also updates zpseudo and mobile for this column. However, zpseudo remains unchanged
    ! if mobile is already false or becomes false within this subroutine.
    !
    ! The structure of the updates done here for zpseudo and mobile requires that this
    ! subroutine be called first for the top layer of snow, then for the 2nd layer down,
    ! etc. - and finally for the bottom layer. Before beginning the loops over layers,
    ! mobile should be initialized to .true. and zpseudo should be initialized to 0.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    real(r8) , intent(in)    :: bi              ! partial density of ice [kg/m3]
    real(r8) , intent(in)    :: forc_wind       ! atmospheric wind speed [m/s]
    real(r8) , intent(in)    :: dz              ! layer depth for this column and level [m]
    real(r8) , intent(inout) :: zpseudo         ! wind drift compaction / pseudo depth for this column at this layer
    logical  , intent(inout) :: mobile          ! whether this snow column is still mobile at this layer (i.e., susceptible to wind drift)
    real(r8) , intent(out)   :: compaction_rate ! rate of compaction of snowpack due to wind drift, for the current column and layer
    !
    ! !LOCAL VARIABLES:
    real(r8) :: Frho        ! Mobility density factor [-]
    real(r8) :: MO          ! Mobility index [-]
    real(r8) :: SI          ! Driftability index [-]
    real(r8) :: gamma_drift ! Scaling factor for wind drift time scale [-]
    real(r8) :: tau_inverse ! Inverse of the effective time scale [1/s]

    real(r8), parameter :: rho_min = 50._r8      ! wind drift compaction / minimum density [kg/m3]
    real(r8), parameter :: rho_max = 350._r8     ! wind drift compaction / maximum density [kg/m3]
    real(r8), parameter :: drift_gs = 0.35e-3_r8 ! wind drift compaction / grain size (fixed value for now)
    real(r8), parameter :: drift_sph = 1.0_r8    ! wind drift compaction / sphericity
    real(r8), parameter :: tau_ref = 48._r8 * 3600._r8  ! wind drift compaction / reference time [s]

    character(len=*), parameter :: subname = 'WindDriftCompaction'
    !-----------------------------------------------------------------------

    if (mobile) then
       Frho = 1.25_r8 - 0.0042_r8*(max(rho_min, bi)-rho_min)
       ! assuming dendricity = 0, sphericity = 1, grain size = 0.35 mm Non-dendritic snow
       MO = 0.34_r8 * (-0.583_r8*drift_gs - 0.833_r8*drift_sph + 0.833_r8) + 0.66_r8*Frho
       SI = -2.868_r8 * exp(-0.085_r8*forc_wind) + 1._r8 + MO

       if (SI > 0.0_r8) then
          SI = min(SI, 3.25_r8)
          ! Increase zpseudo (wind drift / pseudo depth) to the middle of
          ! the pseudo-node for the sake of the following calculation
          zpseudo = zpseudo + 0.5_r8 * dz * (3.25_r8 - SI)
          gamma_drift = SI*exp(-zpseudo/0.1_r8)
          tau_inverse = gamma_drift / tau_ref
          compaction_rate = -max(0.0_r8, rho_max-bi) * tau_inverse
          ! Further increase zpseudo to the bottom of the pseudo-node for
          ! the sake of calculations done on the underlying layer (i.e.,
          ! the next time through the j loop).
          zpseudo = zpseudo + 0.5_r8 * dz * (3.25_r8 - SI)
       else  ! SI <= 0
          mobile = .false.
          compaction_rate = 0._r8
       end if
    else  ! .not. mobile
       compaction_rate = 0._r8
    end if

  end subroutine WindDriftCompaction


  !-----------------------------------------------------------------------
  subroutine Combo(dz,  wliq,  wice, t, dz2, wliq2, wice2, t2)
    !
    ! !DESCRIPTION:
    ! Combines two elements and returns the following combined
    ! variables: dz, t, wliq, wice.
    ! The combined temperature is based on the equation:
    ! the sum of the enthalpies of the two elements =
    ! that of the combined element.
    !
    ! !USES:
    use clm_varcon,  only : cpice, cpliq, tfrz, hfus
    !
    ! !ARGUMENTS:
    implicit none
    real(r8), intent(in)    :: dz2   ! nodal thickness of 2 elements being combined [m]
    real(r8), intent(in)    :: wliq2 ! liquid water of element 2 [kg/m2]
    real(r8), intent(in)    :: wice2 ! ice of element 2 [kg/m2]
    real(r8), intent(in)    :: t2    ! nodal temperature of element 2 [K]
    real(r8), intent(inout) :: dz    ! nodal thickness of 1 elements being combined [m]
    real(r8), intent(inout) :: wliq  ! liquid water of element 1
    real(r8), intent(inout) :: wice  ! ice of element 1 [kg/m2]
    real(r8), intent(inout) :: t     ! nodel temperature of element 1 [K]
    !
    ! !LOCAL VARIABLES:
    real(r8) :: dzc   ! Total thickness of nodes 1 and 2 (dzc=dz+dz2).
    real(r8) :: wliqc ! Combined liquid water [kg/m2]
    real(r8) :: wicec ! Combined ice [kg/m2]
    real(r8) :: tc    ! Combined node temperature [K]
    real(r8) :: h     ! enthalpy of element 1 [J/m2]
    real(r8) :: h2    ! enthalpy of element 2 [J/m2]
    real(r8) :: hc    ! temporary
    !-----------------------------------------------------------------------

    dzc = dz+dz2
    wicec = (wice+wice2)
    wliqc = (wliq+wliq2)
    h = (cpice*wice+cpliq*wliq) * (t-tfrz)+hfus*wliq
    h2= (cpice*wice2+cpliq*wliq2) * (t2-tfrz)+hfus*wliq2

    hc = h + h2
    tc = tfrz + (hc - hfus*wliqc) / (cpice*wicec + cpliq*wliqc)

    dz = dzc
    wice = wicec
    wliq = wliqc
    t = tc

  end subroutine Combo

  !-----------------------------------------------------------------------
  function MassWeightedSnowRadius( rds1, rds2, swtot, zwtot ) result(mass_weighted_snowradius)
    !
    ! !DESCRIPTION:
    ! Calculate the mass weighted snow radius when two layers are combined
    !
    ! !USES:
    use AerosolMod   , only : snw_rds_min
    use SnowSnicarMod, only : snw_rds_max
    implicit none
    ! !ARGUMENTS:
    real(r8), intent(IN) :: rds1         ! Layer 1 radius
    real(r8), intent(IN) :: rds2         ! Layer 2 radius
    real(r8), intent(IN) :: swtot        ! snow water total layer 2
    real(r8), intent(IN) :: zwtot        ! snow water total layer 1
    real(r8) :: mass_weighted_snowradius ! resulting bounded mass weighted snow radius

    SHR_ASSERT( (swtot+zwtot > 0.0_r8), errMsg(sourcefile, __LINE__))
    mass_weighted_snowradius = (rds2*swtot + rds1*zwtot)/(swtot+zwtot)

    if (      mass_weighted_snowradius > snw_rds_max ) then
       mass_weighted_snowradius = snw_rds_max
    else if ( mass_weighted_snowradius < snw_rds_min ) then
       mass_weighted_snowradius = snw_rds_min
    end if
  end function MassWeightedSnowRadius

  !-----------------------------------------------------------------------
  subroutine BuildSnowFilter(bounds, num_nolakec, filter_nolakec, &
       num_snowc, filter_snowc, num_nosnowc, filter_nosnowc)
    !
    ! !DESCRIPTION:
    ! Constructs snow filter for use in vectorized loops for snow hydrology.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    type(bounds_type) , intent(in)  :: bounds
    integer           , intent(in)  :: num_nolakec       ! number of column non-lake points in column filter
    integer           , intent(in)  :: filter_nolakec(:) ! column filter for non-lake points
    integer           , intent(out) :: num_snowc         ! number of column snow points in column filter
    integer           , intent(out) :: filter_snowc(:)   ! column filter for snow points
    integer           , intent(out) :: num_nosnowc       ! number of column non-snow points in column filter
    integer           , intent(out) :: filter_nosnowc(:) ! column filter for non-snow points
    !
    ! !LOCAL VARIABLES:
    integer  :: fc, c
    !-----------------------------------------------------------------------

    ! Build snow/no-snow filters for other subroutines

    num_snowc = 0
    num_nosnowc = 0
    do fc = 1, num_nolakec
       c = filter_nolakec(fc)
       if (col%snl(c) < 0) then
          num_snowc = num_snowc + 1
          filter_snowc(num_snowc) = c
       else
          num_nosnowc = num_nosnowc + 1
          filter_nosnowc(num_nosnowc) = c
       end if
    end do
  end subroutine BuildSnowFilter

  subroutine SnowHydrologySetControlForTesting( set_winddep_snowdensity, set_new_snow_density, &
       set_reset_snow, set_reset_snow_glc, set_reset_snow_glc_ela)
    !
    ! !DESCRIPTION:
    ! Sets some of the control settings for SnowHydrologyMod
    ! NOTE: THIS IS JUST HERE AS AN INTERFACE FOR UNIT TESTING.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    logical, intent(in), optional :: set_winddep_snowdensity  ! Set wind dependent snow density
    integer, intent(in), optional :: set_new_snow_density     ! snow density method
    logical, intent(in), optional :: set_reset_snow           ! whether to reset the snow pack, non-glc_mec points
    logical, intent(in), optional :: set_reset_snow_glc       ! whether to reset the snow pack, glc_mec points
    real(r8), intent(in), optional :: set_reset_snow_glc_ela  ! elevation below which to reset the snow pack if set_reset_snow_glc is true (m)
    !-----------------------------------------------------------------------
    if (present(set_winddep_snowdensity)) then
       wind_dependent_snow_density = set_winddep_snowdensity
    end if
    if (present(set_new_snow_density)) then
       new_snow_density            = set_new_snow_density
    end if
    if (present(set_reset_snow)) then
       reset_snow = set_reset_snow
    end if
    if (present(set_reset_snow_glc)) then
       reset_snow_glc = set_reset_snow_glc
    end if
    if (present(set_reset_snow_glc_ela)) then
       reset_snow_glc_ela = set_reset_snow_glc_ela
    end if

  end subroutine SnowHydrologySetControlForTesting

end module SnowHydrologyMod
