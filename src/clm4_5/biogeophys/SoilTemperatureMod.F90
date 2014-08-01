module SoilTemperatureMod

#include "shr_assert.h"

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Calculates snow and soil temperatures including phase change
  !
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use abortutils     ,  only: endrun
  use perf_mod       ,  only: t_startf, t_stopf
  use decompMod      ,  only: bounds_type
  use shr_log_mod    , only : errMsg => shr_log_errMsg
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: SoilTemperature                     ! Snow and soil temperatures including phase change
  !
  ! The following is only public for the sake of unit testing; it should not be called
  ! directly by CLM code outside this module
  public :: ComputeGroundHeatFluxAndDeriv       ! Computes G and dG/dT on surface of standing water, snow and soil
  public :: ComputeHeatDiffFluxAndFactor        ! Heat diffusion at layer interface and factor used in setting up of banded matrix
  public :: SetRHSVec                           ! Sets up the RHS vector for the numerical solution of temperature for snow/standing-water/soil
  public :: SetRHSVec_Snow                      ! Sets up the RHS vector corresponding to snow layers for Urban+Non-Urban columns
  public :: SetRHSVec_SnowUrban                 ! Sets up the RHS vector corresponding to snow layers for Urban columns
  public :: SetRHSVec_SnowUrbanNonRoad          ! Sets up the RHS vector corresponding to snow layers for Urban columns that are sunwall, shadewall, and roof columns
  public :: SetRHSVec_SnowUrbanRoad             ! Sets up the RHS vector corresponding to snow layers for Urban columns that are pervious, and impervious columns
  public :: SetRHSVec_SnowNonUrban              ! Sets up the RHS vector corresponding to snow layers for Non-Urban columns
  public :: SetRHSVec_StandingSurfaceWater      ! Sets up the RHS vector corresponding to standing water layers for Urban+Non-Urban columns
  public :: SetRHSVec_Soil                      ! Sets up the RHS vector corresponding to soil layers for Urban+Non-Urban columns
  public :: SetRHSVec_SoilUrban                 ! Sets up the RHS vector corresponding to soil layers for Urban columns
  public :: SetRHSVec_SoilUrbanNonRoad          ! Sets up the RHS vector corresponding to soil layers for Urban columns that are pervious, and impervious columns
  public :: SetRHSVec_SoilUrbanRoad             ! Sets up the RHS vector corresponding to soil layers for Urban columns that are pervious, and impervious columns
  public :: SetRHSVec_SoilNonUrban              ! Sets up the RHS vector corresponding to soil layers for Non-Urban columns
  public :: SetRHSVec_Soil_StandingSurfaceWater ! Adds contribution from standing water in the RHS vector corresponding to soil layers
  public :: SetMatrix                           ! Sets up the matrix for the numerical solution of temperature for snow/standing-water/soil
  public :: AssembleMatrixFromSubmatrices       ! Assemble the full matrix from submatrices.
  public :: SetMatrix_Snow                      ! Set up the matrix entries correspodning to snow layers for Urban+Non-Urban columns
  public :: SetMatrix_SnowUrban                 ! Set up the matrix entries correspodning to snow layers for Urban column
  public :: SetMatrix_SnowUrbanNonRoad          ! Set up the matrix entries correspodning to snow layers for Urban column that are sunwall, shadewall, and roof columns
  public :: SetMatrix_SnowUrbanRoad             ! Set up the matrix entries correspodning to snow layers for Urban column that are pervious, and impervious columns
  public :: SetMatrix_SnowNonUrban              ! Set up the matrix entries correspodning to snow layers for Non-Urban column
  public :: SetMatrix_Snow_Soil                 ! Set up the matrix entries correspodning to snow-soil interaction
  public :: SetMatrix_Snow_SoilUrban            ! Set up the matrix entries correspodning to snow-soil interaction for Urban column
  public :: SetMatrix_Snow_SoilUrbanNonRoad     ! Set up the matrix entries correspodning to snow-soil interaction for Urban column that are sunwall, shadewall, and roof columns
  public :: SetMatrix_Snow_SoilUrbanRoad        ! Set up the matrix entries correspodning to snow-soil interaction for Urban column that are pervious, and impervious columns
  public :: SetMatrix_Snow_SoilNonUrban         ! Set up the matrix entries correspodning to snow-soil interaction for Non-Urban column
  public :: SetMatrix_Soil                      ! Set up the matrix entries correspodning to soil layers for Urban+Non-Urban columns
  public :: SetMatrix_SoilUrban                 ! Set up the matrix entries correspodning to soil layers for Urban column
  public :: SetMatrix_SoilUrbanNonRoad          ! Set up the matrix entries correspodning to soil layers for Urban column that are sunwall, shadewall, and roof columns
  public :: SetMatrix_SoilUrbanRoad             ! Set up the matrix entries correspodning to soil layers for Urban column that are pervious, and impervious columns
  public :: SetMatrix_SoilNonUrban              ! Set up the matrix entries correspodning to soil layers for Non-Urban column
  public :: SetMatrix_Soil_Snow                 ! Set up the matrix entries correspodning to soil-snow interction for Urban+Non-Urban columns
  public :: SetMatrix_Soil_SnowUrban            ! Set up the matrix entries correspodning to soil-snow interction for Urban column
  public :: SetMatrix_Soil_SnowUrbanNonRoad     ! Set up the matrix entries correspodning to soil-snow interction for Urban column that are sunwall, shadewall, and roof columns
  public :: SetMatrix_Soil_SnowUrbanRoad        ! Set up the matrix entries correspodning to soil-snow interction for Urban column that are pervious, and impervious columns
  public :: SetMatrix_Soil_SnowNonUrban         ! Set up the matrix entries correspodning to soil-snow interction for Non-Urban column
  public :: SetMatrix_StandingSurfaceWater      ! Set up the matrix entries correspodning to standing surface water
  public :: SetMatrix_StandingSurfaceWater_Soil ! Set up the matrix entries correspodning to standing surface water-soil interaction
  public :: SetMatrix_Soil_StandingSurfaceWater ! Set up the matrix entries correspodning to soil-standing surface water interction
  !
  ! !PRIVATE MEMBER FUNCTIONS:
  private :: SoilThermProp                    ! Set therm conduct. and heat cap of snow/soil layers
  private :: PhaseChangeH2osfc                ! When surface water freezes move ice to bottom snow layer
  private :: PhaseChange_beta                 ! Calculation of the phase change within snow and soil layers
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine SoilTemperature(bounds, num_urbanl, filter_urbanl, &
       num_nolakec, filter_nolakec, xmf, fact, c_h2osfc, xmf_h2osfc)
    !
    ! !DESCRIPTION:
    ! Snow and soil temperatures including phase change
    ! o The volumetric heat capacity is calculated as a linear combination
    !   in terms of the volumetric fraction of the constituent phases.
    ! o The thermal conductivity of soil is computed from
    !   the algorithm of Johansen (as reported by Farouki 1981), and the
    !   conductivity of snow is from the formulation used in
    !   SNTHERM (Jordan 1991).
    ! o Boundary conditions:
    !   F = Rnet - Hg - LEg (top),  F= 0 (base of the soil column).
    ! o Soil / snow temperature is predicted from heat conduction
    !   in 10 soil layers and up to 5 snow layers.
    !   The thermal conductivities at the interfaces between two
    !   neighboring layers (j, j+1) are derived from an assumption that
    !   the flux across the interface is equal to that from the node j
    !   to the interface and the flux from the interface to the node j+1.
    !   The equation is solved using the Crank-Nicholson method and
    !   results in a tridiagonal system equation.
    !
    ! !USES:
    use clmtype
    use clm_time_manager, only : get_step_size
    use shr_infnan_mod  , only : nan => shr_infnan_nan, assignment(=)
    use clm_varcon      , only : cnfac, cpice, cpliq, denh2o
    use column_varcon  , only : icol_roof, icol_sunwall, icol_shadewall
    use landunit_varcon, only : istwet, istice, istice_mec, istsoil, istcrop
    use clm_varpar     , only : nlevsno, nlevgrnd, nlevurb
    use BandDiagonalMod, only : BandDiagonal
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds                               ! bounds
    integer , intent(in)  :: num_nolakec                                  ! number of column non-lake points in column filter
    integer , intent(in)  :: filter_nolakec(:)                            ! column filter for non-lake points
    integer , intent(in)  :: num_urbanl                                   ! number of urban landunits in clump
    integer , intent(in)  :: filter_urbanl(:)                             ! urban landunit filter
    real(r8), intent(out) :: xmf( bounds%begc: )                          ! total latent heat of phase change of ground water [col]
    real(r8), intent(out) :: fact( bounds%begc: , -nlevsno+1: )           ! used in computing tridiagonal matrix [col, lev]
    real(r8), intent(out) :: xmf_h2osfc( bounds%begc: )                   ! latent heat of phase change of surface water [col]
    real(r8), intent(out) :: c_h2osfc( bounds%begc: )                     ! heat capacity of surface water [col]
    !
    ! !LOCAL VARIABLES:
    integer  :: j,c,l,g,pi                                                ! indices
    integer  :: fc                                                        ! lake filtered column indices
    integer  :: fl                                                        ! urban filtered landunit indices
    integer  :: jtop(bounds%begc:bounds%endc)                             ! top level at each column
    real(r8) :: dtime                                                     ! land model time step (sec)
    real(r8) :: cv (bounds%begc:bounds%endc,-nlevsno+1:nlevgrnd)          ! heat capacity [J/(m2 K)]
    real(r8) :: tk (bounds%begc:bounds%endc,-nlevsno+1:nlevgrnd)          ! thermal conductivity [W/(m K)]
    real(r8) :: fn (bounds%begc:bounds%endc,-nlevsno+1:nlevgrnd)          ! heat diffusion through the layer interface [W/m2]
    real(r8) :: fn1(bounds%begc:bounds%endc,-nlevsno+1:nlevgrnd)          ! heat diffusion through the layer interface [W/m2]
    real(r8) :: dzm                                                       ! used in computing tridiagonal matrix
    real(r8) :: dzp                                                       ! used in computing tridiagonal matrix
    real(r8) :: sabg_lyr_col(bounds%begc:bounds%endc,-nlevsno+1:1)        ! absorbed solar radiation (col,lyr) [W/m2]
    real(r8) :: hs_top(bounds%begc:bounds%endc)                           ! net energy flux into surface layer (col) [W/m2]
    logical  :: cool_on(bounds%begl:bounds%endl)                          ! is urban air conditioning on?
    logical  :: heat_on(bounds%begl:bounds%endl)                          ! is urban heating on?
    real(r8) :: fn_h2osfc(bounds%begc:bounds%endc)                        ! heat diffusion through standing-water/soil interface [W/m2]
    real(r8) :: dz_h2osfc(bounds%begc:bounds%endc)                        ! height of standing surface water [m]
    integer, parameter :: nband=5                                         ! number of non-zero entries of the banded matrix
    real(r8) :: bmatrix(bounds%begc:bounds%endc,nband,-nlevsno:nlevgrnd)  ! banded matrix for numerical solution of temperature
    real(r8) :: tvector(bounds%begc:bounds%endc,-nlevsno:nlevgrnd)        ! initial temperature solution [Kelvin]
    real(r8) :: rvector(bounds%begc:bounds%endc,-nlevsno:nlevgrnd)        ! RHS vector for numerical solution of temperature
    real(r8) :: tk_h2osfc(bounds%begc:bounds%endc)                        ! thermal conductivity of h2osfc [W/(m K)] [col]
    real(r8) :: dhsdT(bounds%begc:bounds%endc)                            ! temperature derivative of "hs" [col]
    real(r8) :: hs_soil(bounds%begc:bounds%endc)                          ! heat flux on soil [W/m2]
    real(r8) :: hs_top_snow(bounds%begc:bounds%endc)                      ! heat flux on top snow layer [W/m2]
    real(r8) :: hs_h2osfc(bounds%begc:bounds%endc)                        ! heat flux on standing water [W/m2]
    integer  :: jbot(bounds%begc:bounds%endc)                             ! bottom level at each column
    !-----------------------------------------------------------------------

    ! Enforce expected array sizes
    SHR_ASSERT_ALL((ubound(xmf)        == (/bounds%endc/)),           errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(fact)       == (/bounds%endc, nlevgrnd/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(xmf_h2osfc) == (/bounds%endc/)),           errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(c_h2osfc)   == (/bounds%endc/)),           errMsg(__FILE__, __LINE__))

   associate(&
   ltype                     =>    lun%itype                     , & ! Input:  [integer (:)]  landunit type
   t_building                =>    lps%t_building                , & ! Input:  [real(r8) (:)]  internal building temperature (K)
   t_building_max            =>    lps%t_building_max            , & ! Input:  [real(r8) (:)]  maximum internal building temperature (K)
   t_building_min            =>    lps%t_building_min            , & ! Input:  [real(r8) (:)]  minimum internal building temperature (K)
   urbpoi                    =>    lun%urbpoi                    , & ! Input:  [logical (:)]  true => landunit is an urban point
   frac_sno_eff              =>    cps%frac_sno_eff              , & ! Input:  [real(r8) (:)]  eff. fraction of ground covered by snow (0 to 1)
   frac_h2osfc               =>    cps%frac_h2osfc               , & ! Input:  [real(r8) (:)]  fraction of ground covered by surface water (0 to 1)
   h2osfc                    =>    cws%h2osfc                    , & ! Input:  [real(r8) (:)]  surface water (mm)                      
   t_h2osfc                  =>    ces%t_h2osfc                  , & ! Input:  [real(r8) (:)]  surface water temperature
   ctype                     =>    col%itype                     , & ! Input:  [integer (:)]  column type
   clandunit                 =>    col%landunit                  , & ! Input:  [integer (:)]  column's landunit
   snl                       =>    cps%snl                       , & ! Input:  [integer (:)]  number of snow layers
   t_grnd                    =>    ces%t_grnd                    , & ! Input:  [real(r8) (:)]  ground surface temperature [K]
   hc_soi                    =>    ces%hc_soi                    , & ! Input:  [real(r8) (:)]  soil heat content (MJ/m2)
   hc_soisno                 =>    ces%hc_soisno                 , & ! Input:  [real(r8) (:)]  soil plus snow plus lake heat content (MJ/m2)
   eflx_fgr12                =>    cef%eflx_fgr12                , & ! Input:  [real(r8) (:)]  heat flux between soil layer 1 and 2 (W/m2)
   eflx_fgr                  =>    cef%eflx_fgr                  , & ! Input:  [real(r8) (:,:)]  (rural) soil downward heat flux (W/m2) (1:nlevgrnd)
   zi                        =>    cps%zi                        , & ! Input:  [real(r8) (:,:)]  interface level below a "z" level (m)
   dz                        =>    cps%dz                        , & ! Input:  [real(r8) (:,:)]  layer depth (m)
   z                         =>    cps%z                         , & ! Input:  [real(r8) (:,:)]  layer thickness (m)
   t_soisno                  =>    ces%t_soisno                  , & ! Input:  [real(r8) (:,:)]  soil temperature (Kelvin)
   eflx_building_heat        =>    cef%eflx_building_heat        , & ! Output: [real(r8) (:)]  heat flux from urban building interior to walls, roof (W/m**2)
   tssbef                    =>    ces%tssbef                    , & ! Input:  [real(r8) (:,:)]  temperature at previous time step [K]
   eflx_urban_ac             =>    cef%eflx_urban_ac             , & ! Output: [real(r8) (:)]  urban air conditioning flux (W/m**2)
   eflx_urban_heat           =>    cef%eflx_urban_heat           , & ! Output: [real(r8) (:)]  urban heating flux (W/m**2)
   begc                      =>    bounds%begc                   , & ! Input:  [integer ] beginning column index
   endc                      =>    bounds%endc                     & ! Input:  [integer ] ending column index
   )

    ! Get step size

    dtime = get_step_size()

    ! Restrict internal building temperature to between min and max
    ! and determine if heating or air conditioning is on
    do fl = 1,num_urbanl
       l = filter_urbanl(fl)
       if (urbpoi(l)) then
          cool_on(l) = .false.
          heat_on(l) = .false.
          if (t_building(l) > t_building_max(l)) then
            t_building(l) = t_building_max(l)
            cool_on(l) = .true.
            heat_on(l) = .false.
          else if (t_building(l) < t_building_min(l)) then
            t_building(l) = t_building_min(l)
            cool_on(l) = .false.
            heat_on(l) = .true.
          end if
       end if
    end do

    ! set up compact matrix for band diagonal solver, requires additional
    !     sub/super diagonals (1 each), and one additional row for t_h2osfc
    jtop = -9999
    do fc = 1,num_nolakec
       c = filter_nolakec(fc)
       jtop(c) = snl(c)
       ! compute jbot
       if ((ctype(c) == icol_sunwall .or. ctype(c) == icol_shadewall &
            .or. ctype(c) == icol_roof) ) then
          jbot(c) = nlevurb
       else
          jbot(c) = nlevgrnd
       endif
    end do

    ! Compute ground surface and soil temperatures

    ! Thermal conductivity and Heat capacity

    tk_h2osfc(begc:endc) = nan
    call SoilThermProp(bounds, num_nolakec, filter_nolakec, &
         tk( begc:endc, -nlevsno+1: ),                      &
         cv( begc:endc, -nlevsno+1:),                       &
         tk_h2osfc( begc:endc))

    ! Net ground heat flux into the surface and its temperature derivative
    ! Added a pfts loop here to get the average of hs and dhsdT over
    ! all PFTs on the column. Precalculate the terms that do not depend on PFT.

    call ComputeGroundHeatFluxAndDeriv(bounds, num_nolakec, filter_nolakec, &
         hs_h2osfc( begc:endc ),                                            &
         hs_top_snow( begc:endc ),                                          &
         hs_soil( begc:endc ),                                              &
         hs_top( begc:endc ),                                               &
         dhsdT( begc:endc ),                                                &
         sabg_lyr_col( begc:endc, -nlevsno+1: ))

    ! Determine heat diffusion through the layer interface and factor used in computing
    ! banded diagonal matrix and set up vector r and vectors a, b, c that define banded
    ! diagonal matrix and solve system

    call ComputeHeatDiffFluxAndFactor(bounds, num_nolakec, filter_nolakec, &
         dtime,                                                            &
         tk( begc:endc, -nlevsno+1: ),                                     &
         cv( begc:endc, -nlevsno+1: ),                                     &
         fn( begc:endc, -nlevsno+1: ),                                     &
         fact( begc:endc, -nlevsno+1: ))

    ! compute thermal properties of h2osfc

    do fc = 1,num_nolakec
       c = filter_nolakec(fc)
       dz_h2osfc(c) = max(1.0e-6_r8,1.0e-3*h2osfc(c))
       c_h2osfc(c)  = cpliq*denh2o*dz_h2osfc(c) !"areametric" heat capacity [J/K/m^2]
    enddo


    ! Set up vector r and vectors a, b, c that define tridiagonal

    call SetRHSVec(bounds, num_nolakec, filter_nolakec, &
         dtime,                                         &
         hs_h2osfc( begc:endc ),                        &
         hs_top_snow( begc:endc ),                      &
         hs_soil( begc:endc ),                          &
         hs_top( begc:endc ),                           &
         dhsdT( begc:endc ),                            &
         sabg_lyr_col (begc:endc, -nlevsno+1: ),        &
         tk( begc:endc, -nlevsno+1: ),                  &
         tk_h2osfc( begc:endc ),                        &
         fact( begc:endc, -nlevsno+1: ),                &
         fn( begc:endc, -nlevsno+1: ),                  &
         c_h2osfc( begc:endc ),                         &
         dz_h2osfc( begc:endc ),                        &
         rvector( begc:endc, -nlevsno: ))

    ! Set up the banded diagonal matrix

    call SetMatrix(bounds, num_nolakec, filter_nolakec, &
         dtime,                                         &
         nband,                                         &
         dhsdT( begc:endc ),                            &
         tk( begc:endc, -nlevsno+1: ),                  &
         tk_h2osfc( begc:endc ),                        &
         fact( begc:endc, -nlevsno+1: ),                &
         c_h2osfc( begc:endc ),                         &
         dz_h2osfc( begc:endc ),                        &
         bmatrix( begc:endc, 1:, -nlevsno: ))

    ! initialize initial temperature vector

    tvector(begc:endc, :) = nan
    do fc = 1,num_nolakec
       c = filter_nolakec(fc)
       do j = snl(c)+1, 0
          tvector(c,j-1) = t_soisno(c,j)
       end do

       ! surface water layer has two coefficients
       tvector(c,0) = t_h2osfc(c)

       ! soil layers; top layer will have one offset and one extra coefficient
       tvector(c,1:nlevgrnd) = t_soisno(c,1:nlevgrnd)

    enddo

    call t_startf( 'SoilTempBandDiag')

    ! Solve the system

    call BandDiagonal(bounds, -nlevsno, nlevgrnd, jtop(begc:endc), jbot(begc:endc), &
         num_nolakec, filter_nolakec, nband, bmatrix(begc:endc, :, :), &
         rvector(begc:endc, :), tvector(begc:endc, :))
     call t_stopf( 'SoilTempBandDiag')

    ! return temperatures to original array
    do fc = 1,num_nolakec
       c = filter_nolakec(fc)
       do j = snl(c)+1, 0
          t_soisno(c,j) = tvector(c,j-1) !snow layers
       end do
       t_soisno(c,1:nlevgrnd)   = tvector(c,1:nlevgrnd)  !soil layers

       if (frac_h2osfc(c) == 0._r8) then
          t_h2osfc(c)=t_soisno(c,1)
       else
          t_h2osfc(c)              = tvector(c,0)           !surface water
       endif
    enddo

    ! Melting or Freezing

    do j = -nlevsno+1,nlevgrnd
       do fc = 1,num_nolakec
          c = filter_nolakec(fc)
          l = clandunit(c)
          if ((ctype(c) == icol_sunwall .or. ctype(c) == icol_shadewall &
              .or. ctype(c) == icol_roof) .and. j <= nlevurb) then
             if (j >= snl(c)+1) then
                if (j <= nlevurb-1) then
                   fn1(c,j) = tk(c,j)*(t_soisno(c,j+1)-t_soisno(c,j))/(z(c,j+1)-z(c,j))
                else if (j == nlevurb) then
                   ! For urban sunwall, shadewall, and roof columns, there is a non-zero heat flux across
                   ! the bottom "soil" layer and the equations are derived assuming a prescribed internal
                   ! building temperature. (See Oleson urban notes of 6/18/03).
                   ! Note new formulation for fn, this will be used below in net energey flux computations
                   fn1(c,j) = tk(c,j) * (t_building(l) - t_soisno(c,j))/(zi(c,j) - z(c,j))
                   fn(c,j)  = tk(c,j) * (t_building(l) - tssbef(c,j))/(zi(c,j) - z(c,j))
                end if
             end if
          else if (ctype(c) /= icol_sunwall .and. ctype(c) /= icol_shadewall &
                   .and. ctype(c) /= icol_roof) then
             if (j >= snl(c)+1) then
                if (j <= nlevgrnd-1) then
                   fn1(c,j) = tk(c,j)*(t_soisno(c,j+1)-t_soisno(c,j))/(z(c,j+1)-z(c,j))
                else if (j == nlevgrnd) then
                   fn1(c,j) = 0._r8
                end if
             end if
          end if
       end do
    end do

    do fc = 1,num_nolakec
       c = filter_nolakec(fc)
       l = clandunit(c)
       if (urbpoi(l)) then
         if (ctype(c) == icol_sunwall .or. ctype(c) == icol_shadewall .or. ctype(c) == icol_roof) then
           eflx_building_heat(c) = cnfac*fn(c,nlevurb) + (1-cnfac)*fn1(c,nlevurb)
         else
           eflx_building_heat(c) = 0._r8
         end if
         if (cool_on(l)) then
           eflx_urban_ac(c) = abs(eflx_building_heat(c))
           eflx_urban_heat(c) = 0._r8
         else if (heat_on(l)) then
           eflx_urban_ac(c) = 0._r8
           eflx_urban_heat(c) = abs(eflx_building_heat(c))
         else
           eflx_urban_ac(c) = 0._r8
           eflx_urban_heat(c) = 0._r8
         end if
       end if
    end do

    xmf_h2osfc=0.
    ! compute phase change of h2osfc
    call PhaseChangeH2osfc (bounds, num_nolakec, filter_nolakec, &
         fact(bounds%begc:bounds%endc, :),                       &
         dhsdT(bounds%begc:bounds%endc),                         &
         c_h2osfc(bounds%begc:bounds%endc),                      &
         xmf_h2osfc(bounds%begc:bounds%endc))

    call Phasechange_beta (bounds, num_nolakec, filter_nolakec, &
         fact(bounds%begc:bounds%endc, :),                      &
         dhsdT(bounds%begc:bounds%endc),                        &
         xmf(bounds%begc:bounds%endc))

    do fc = 1,num_nolakec
       c = filter_nolakec(fc)
       ! this expression will (should) work whether there is snow or not
       if (snl(c) < 0) then
          if(frac_h2osfc(c) /= 0._r8) then
              t_grnd(c) = frac_sno_eff(c) * t_soisno(c,snl(c)+1) &
                   + (1.0_r8 - frac_sno_eff(c) - frac_h2osfc(c)) * t_soisno(c,1) &
                   + frac_h2osfc(c) * t_h2osfc(c)
          else
              t_grnd(c) = frac_sno_eff(c) * t_soisno(c,snl(c)+1) &
                   + (1.0_r8 - frac_sno_eff(c)) * t_soisno(c,1)
          end if

       else
          if(frac_h2osfc(c) /= 0._r8) then
             t_grnd(c) = (1 - frac_h2osfc(c)) * t_soisno(c,1) + frac_h2osfc(c) * t_h2osfc(c)
          else
             t_grnd(c) = t_soisno(c,1)
          end if
       endif
    end do

    ! Initialize soil heat content
    do fc = 1,num_nolakec
       c = filter_nolakec(fc)
       l = clandunit(c)
       if (.not. urbpoi(l)) then
         hc_soisno(c) = 0._r8
         hc_soi(c)    = 0._r8
       end if
       eflx_fgr12(c)= 0._r8
    end do

    ! Calculate soil heat content and soil plus snow heat content
    do j = -nlevsno+1,nlevgrnd
       do fc = 1,num_nolakec
          c = filter_nolakec(fc)
          l = clandunit(c)

          if (j == 1) then ! this only needs to be done once
             eflx_fgr12(c) = -cnfac*fn(c,1) - (1._r8-cnfac)*fn1(c,1)
          end if
          if (j > 0 .and. j < nlevgrnd .and. (ltype(l) == istsoil .or. ltype(l) == istcrop)) then
             eflx_fgr(c,j) = -cnfac*fn(c,j) - (1._r8-cnfac)*fn1(c,j)
          else if (j == nlevgrnd .and. (ltype(l) == istsoil .or. ltype(l) == istcrop)) then
             eflx_fgr(c,j) = 0._r8
          end if

          if (.not. urbpoi(l)) then
            if (j >= snl(c)+1) then
               hc_soisno(c) = hc_soisno(c) + cv(c,j)*t_soisno(c,j) / 1.e6_r8
            endif
            if (j >= 1) then
               hc_soi(c) = hc_soi(c) + cv(c,j)*t_soisno(c,j) / 1.e6_r8
            end if
          end if
       end do
    end do

    end associate

   end subroutine SoilTemperature

   !-----------------------------------------------------------------------
   subroutine SoilThermProp (bounds,  num_nolakec, filter_nolakec, tk, cv, &
        tk_h2osfc)
     !
     ! !DESCRIPTION:
     ! Calculation of thermal conductivities and heat capacities of
     ! snow/soil layers
     ! (1) The volumetric heat capacity is calculated as a linear combination
     !     in terms of the volumetric fraction of the constituent phases.
     !
     ! (2) The thermal conductivity of soil is computed from the algorithm of
     !     Johansen (as reported by Farouki 1981), and of snow is from the
     !     formulation used in SNTHERM (Jordan 1991).
     ! The thermal conductivities at the interfaces between two neighboring
     ! layers (j, j+1) are derived from an assumption that the flux across
     ! the interface is equal to that from the node j to the interface and the
     ! flux from the interface to the node j+1.
     !
     ! !USES:
     use clmtype
     use clm_varcon  , only : denh2o, denice, tfrz, tkwat, tkice, tkair, &
                              cpice,  cpliq, thk_bedrock
     use landunit_varcon,only:istice, istice_mec, istwet
     use column_varcon,only : icol_roof, icol_sunwall, icol_shadewall, icol_road_perv, icol_road_imperv
     use clm_varpar  , only : nlevsno, nlevgrnd, nlevurb, nlevsoi
     use clm_varctl  , only : iulog
     !
     ! !ARGUMENTS:
     implicit none
     type(bounds_type), intent(in) :: bounds                   ! bounds
     integer , intent(in)  :: num_nolakec                      ! number of column non-lake points in column filter
     integer , intent(in)  :: filter_nolakec(:)                ! column filter for non-lake points
     real(r8), intent(out) :: cv( bounds%begc: , -nlevsno+1: ) ! heat capacity [J/(m2 K)] [col, lev]
     real(r8), intent(out) :: tk( bounds%begc: , -nlevsno+1: ) ! thermal conductivity at the layer interface [W/(m K)] [col, lev]
     real(r8), intent(out) :: tk_h2osfc( bounds%begc: )        ! thermal conductivity of h2osfc [W/(m K)] [col]
     !
     ! !LOCAL VARIABLES:
     integer  :: l,c,j                     ! indices
     integer  :: fc                        ! lake filtered column indices
     real(r8) :: dksat                     ! thermal conductivity for saturated soil (j/(k s m))
     real(r8) :: dke                       ! kersten number
     real(r8) :: fl                        ! volume fraction of liquid or unfrozen water to total water
     real(r8) :: satw                      ! relative total water content of soil.
     real(r8) :: zh2osfc
     !-----------------------------------------------------------------------

    call t_startf( 'SoilThermProp' )

    ! Enforce expected array sizes
    SHR_ASSERT_ALL((ubound(cv)        == (/bounds%endc, nlevgrnd/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(tk)        == (/bounds%endc, nlevgrnd/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(tk_h2osfc) == (/bounds%endc/)),           errMsg(__FILE__, __LINE__))

   associate(& 
   ltype                     =>    lun%itype               , & ! Input:  [integer (:)]  landunit type                            
   h2osfc                    =>    cws%h2osfc              , & ! Input:  [real(r8) (:)]  surface (mm H2O)                        
   frac_sno                  =>    cps%frac_sno_eff        , & ! Input:  [real(r8) (:)]  fractional snow covered area            
   ctype                     =>    col%itype               , & ! Input:  [integer (:)]  column type                              
   clandunit                 =>   col%landunit             , & ! Input:  [integer (:)]  column's landunit                        
   snl                       =>    cps%snl                 , & ! Input:  [integer (:)]  number of snow layers                    
   h2osno                    =>    cws%h2osno              , & ! Input:  [real(r8) (:)]  snow water (mm H2O)                     
   watsat                    =>    cps%watsat              , & ! Input:  [real(r8) (:,:)]  volumetric soil water at saturation (porosity)
   tksatu                    =>    cps%tksatu              , & ! Input:  [real(r8) (:,:)]  thermal conductivity, saturated soil [W/m-K]
   tkmg                      =>    cps%tkmg                , & ! Input:  [real(r8) (:,:)]  thermal conductivity, soil minerals  [W/m-K]
   tkdry                     =>    cps%tkdry               , & ! Input:  [real(r8) (:,:)]  thermal conductivity, dry soil (W/m/Kelvin)
   csol                      =>    cps%csol                , & ! Input:  [real(r8) (:,:)]  heat capacity, soil solids (J/m**3/Kelvin)
   thk                       =>    cps%thk                 , & ! Output: [real(r8) (:,:)]  thermal conductivity of each layer  [W/m-K] (-nlevsno+1:nlevgrnd)
   dz                        =>    cps%dz                  , & ! Input:  [real(r8) (:,:)]  layer depth (m)                       
   zi                        =>    cps%zi                  , & ! Input:  [real(r8) (:,:)]  interface level below a "z" level (m) 
   z                         =>    cps%z                   , & ! Input:  [real(r8) (:,:)]  layer thickness (m)                   
   t_soisno                  =>    ces%t_soisno            , & ! Input:  [real(r8) (:,:)]  soil temperature (Kelvin)             
   h2osoi_liq                =>    cws%h2osoi_liq          , & ! Input:  [real(r8) (:,:)]  liquid water (kg/m2)                  
   h2osoi_ice                =>    cws%h2osoi_ice          , & ! Input:  [real(r8) (:,:)]  ice lens (kg/m2)                      
   bw                        =>    cws%bw                  , & ! Output: [real(r8) (:,:)]  partial density of water in the snow pack (ice + liquid) [kg/m3] (-nlevsno+1:0)
   tk_wall                   =>    lps%tk_wall             , & ! Input:  [real(r8) (:,:)]  thermal conductivity of urban wall    
   tk_roof                   =>    lps%tk_roof             , & ! Input:  [real(r8) (:,:)]  thermal conductivity of urban roof    
   tk_improad                =>    lps%tk_improad          , & ! Input:  [real(r8) (:,:)]  thermal conductivity of urban impervious road
   cv_wall                   =>    lps%cv_wall             , & ! Input:  [real(r8) (:,:)]  thermal conductivity of urban wall    
   cv_roof                   =>    lps%cv_roof             , & ! Input:  [real(r8) (:,:)]  thermal conductivity of urban roof    
   cv_improad                =>    lps%cv_improad          , & ! Input:  [real(r8) (:,:)]  thermal conductivity of urban impervious road
   nlev_improad              =>    lps%nlev_improad          & ! Input:  [integer (:)]  number of impervious road layers         
   )

    ! Thermal conductivity of soil from Farouki (1981)
    ! Urban values are from Masson et al. 2002, Evaluation of the Town Energy Balance (TEB)
    ! scheme with direct measurements from dry districts in two cities, J. Appl. Meteorol.,
    ! 41, 1011-1026.

    do j = -nlevsno+1,nlevgrnd
       do fc = 1, num_nolakec
          c = filter_nolakec(fc)

          ! Only examine levels from 1->nlevgrnd
          if (j >= 1) then    
             l = clandunit(c)
             if ((ctype(c) == icol_sunwall .OR. ctype(c) == icol_shadewall) .and. j <= nlevurb) then
                thk(c,j) = tk_wall(l,j)
             else if (ctype(c) == icol_roof .and. j <= nlevurb) then
                thk(c,j) = tk_roof(l,j)
             else if (ctype(c) == icol_road_imperv .and. j >= 1 .and. j <= nlev_improad(l)) then
                thk(c,j) = tk_improad(l,j)
             else if (ltype(l) /= istwet .AND. ltype(l) /= istice .AND. ltype(l) /= istice_mec &
                      .AND. ctype(c) /= icol_sunwall .AND. ctype(c) /= icol_shadewall .AND. &
                      ctype(c) /= icol_roof) then

                satw = (h2osoi_liq(c,j)/denh2o + h2osoi_ice(c,j)/denice)/(dz(c,j)*watsat(c,j))
                satw = min(1._r8, satw)
                if (satw > .1e-6_r8) then
                   if (t_soisno(c,j) >= tfrz) then       ! Unfrozen soil
                      dke = max(0._r8, log10(satw) + 1.0_r8)
                   else                               ! Frozen soil
                      dke = satw
                   end if
                   fl = (h2osoi_liq(c,j)/(denh2o*dz(c,j))) / (h2osoi_liq(c,j)/(denh2o*dz(c,j)) + &
                                                              h2osoi_ice(c,j)/(denice*dz(c,j)))
                   dksat = tkmg(c,j)*tkwat**(fl*watsat(c,j))*tkice**((1._r8-fl)*watsat(c,j))
                   thk(c,j) = dke*dksat + (1._r8-dke)*tkdry(c,j)
                else
                   thk(c,j) = tkdry(c,j)
                endif
                if (j > nlevsoi) thk(c,j) = thk_bedrock
             else if (ltype(l) == istice .OR. ltype(l) == istice_mec) then
                thk(c,j) = tkwat
                if (t_soisno(c,j) < tfrz) thk(c,j) = tkice
             else if (ltype(l) == istwet) then                         
                if (j > nlevsoi) then 
                   thk(c,j) = thk_bedrock
                else
                   thk(c,j) = tkwat
                   if (t_soisno(c,j) < tfrz) thk(c,j) = tkice
                endif
             endif
          endif

          ! Thermal conductivity of snow, which from Jordan (1991) pp. 18
          ! Only examine levels from snl(c)+1 -> 0 where snl(c) < 1
          if (snl(c)+1 < 1 .AND. (j >= snl(c)+1) .AND. (j <= 0)) then  
             bw(c,j) = (h2osoi_ice(c,j)+h2osoi_liq(c,j))/(frac_sno(c)*dz(c,j))
             thk(c,j) = tkair + (7.75e-5_r8 *bw(c,j) + 1.105e-6_r8*bw(c,j)*bw(c,j))*(tkice-tkair)
          end if

       end do
    end do

    ! Thermal conductivity at the layer interface

    do j = -nlevsno+1,nlevgrnd
       do fc = 1,num_nolakec
          c = filter_nolakec(fc)
          if ((ctype(c) == icol_sunwall .or. ctype(c) == icol_shadewall &
               .or. ctype(c) == icol_roof) .and. j <= nlevurb) then
            if (j >= snl(c)+1 .AND. j <= nlevurb-1) then
               tk(c,j) = thk(c,j)*thk(c,j+1)*(z(c,j+1)-z(c,j)) &
                         /(thk(c,j)*(z(c,j+1)-zi(c,j))+thk(c,j+1)*(zi(c,j)-z(c,j)))
            else if (j == nlevurb) then

               ! For urban sunwall, shadewall, and roof columns, there is a non-zero heat flux across
               ! the bottom "soil" layer and the equations are derived assuming a prescribed internal
               ! building temperature. (See Oleson urban notes of 6/18/03).
               tk(c,j) = thk(c,j)
            end if
          else if (ctype(c) /= icol_sunwall .and. ctype(c) /= icol_shadewall &
                   .and. ctype(c) /= icol_roof) then
            if (j >= snl(c)+1 .AND. j <= nlevgrnd-1) then
               tk(c,j) = thk(c,j)*thk(c,j+1)*(z(c,j+1)-z(c,j)) &
                         /(thk(c,j)*(z(c,j+1)-zi(c,j))+thk(c,j+1)*(zi(c,j)-z(c,j)))
            else if (j == nlevgrnd) then
               tk(c,j) = 0._r8
            end if
          end if
       end do
    end do

    ! calculate thermal conductivity of h2osfc
    do fc = 1, num_nolakec
       c = filter_nolakec(fc)
       zh2osfc=1.0e-3*(0.5*h2osfc(c)) !convert to [m] from [mm]
       tk_h2osfc(c)= tkwat*thk(c,1)*(z(c,1)+zh2osfc) &
                       /(tkwat*z(c,1)+thk(c,1)*zh2osfc)
    enddo

    ! Soil heat capacity, from de Vires (1963)
    ! Urban values are from Masson et al. 2002, Evaluation of the Town Energy Balance (TEB)
    ! scheme with direct measurements from dry districts in two cities, J. Appl. Meteorol.,
    ! 41, 1011-1026.

    do j = 1, nlevgrnd
       do fc = 1,num_nolakec
          c = filter_nolakec(fc)
          l = clandunit(c)
          if ((ctype(c) == icol_sunwall .OR. ctype(c) == icol_shadewall) .and. j <= nlevurb) then
             cv(c,j) = cv_wall(l,j) * dz(c,j)
          else if (ctype(c) == icol_roof .and. j <= nlevurb) then
             cv(c,j) = cv_roof(l,j) * dz(c,j)
          else if (ctype(c) == icol_road_imperv .and. j >= 1 .and. j <= nlev_improad(l)) then
             cv(c,j) = cv_improad(l,j) * dz(c,j)
          else if (ltype(l) /= istwet .AND. ltype(l) /= istice .AND. ltype(l) /= istice_mec &
                   .AND. ctype(c) /= icol_sunwall .AND. ctype(c) /= icol_shadewall .AND. &
                   ctype(c) /= icol_roof) then
             cv(c,j) = csol(c,j)*(1-watsat(c,j))*dz(c,j) + (h2osoi_ice(c,j)*cpice + h2osoi_liq(c,j)*cpliq)
          else if (ltype(l) == istwet) then 
             cv(c,j) = (h2osoi_ice(c,j)*cpice + h2osoi_liq(c,j)*cpliq)
             if (j > nlevsoi) cv(c,j) = csol(c,j)*dz(c,j)
          else if (ltype(l) == istice .OR. ltype(l) == istice_mec) then
             cv(c,j) = (h2osoi_ice(c,j)*cpice + h2osoi_liq(c,j)*cpliq)
          endif
          if (j == 1) then
             if (snl(c)+1 == 1 .AND. h2osno(c) > 0._r8) then
                cv(c,j) = cv(c,j) + cpice*h2osno(c)
             end if
          end if
       enddo
    end do

    ! Snow heat capacity

    do j = -nlevsno+1,0
       do fc = 1,num_nolakec
          c = filter_nolakec(fc)
          if (snl(c)+1 < 1 .and. j >= snl(c)+1) then
             cv(c,j) = cpliq*h2osoi_liq(c,j) + cpice*h2osoi_ice(c,j)
          end if
       end do
    end do
    call t_stopf( 'SoilThermProp' )

    end associate 
   end subroutine SoilThermProp

   !-----------------------------------------------------------------------
   subroutine PhaseChangeH2osfc (bounds, num_nolakec, filter_nolakec, fact, &
        dhsdT,c_h2osfc,xmf_h2osfc)
     !
     ! !DESCRIPTION:
     ! Only freezing is considered.  When water freezes, move ice to bottom snow layer.
     !
     ! !USES:
     use clmtype
     use clm_time_manager, only : get_step_size
     use clm_varcon  , only : tfrz, hfus, grav,denice,cnfac,cpice
     use clm_varpar  , only : nlevsno, nlevgrnd
     use clm_varctl  , only : iulog
     !
     ! !ARGUMENTS:
     implicit none
     type(bounds_type), intent(in) :: bounds                         ! bounds
     integer , intent(in)    :: num_nolakec                          ! number of column non-lake points in column filter
     integer , intent(in)    :: filter_nolakec(:)                    ! column filter for non-lake points
     real(r8), intent(inout) :: fact  ( bounds%begc: , -nlevsno+1: ) ! temporary [col, lev]
     real(r8), intent(in)    :: dhsdT ( bounds%begc: )               ! temperature derivative of "hs" [col]
     real(r8), intent(in)    :: c_h2osfc( bounds%begc: )             ! heat capacity of surface water [col]
     real(r8), intent(out)   :: xmf_h2osfc( bounds%begc: )           ! latent heat of phase change of surface water [col]
     !
     ! !LOCAL VARIABLES:
     integer  :: j,c,g                              !do loop index
     integer  :: fc                                 !lake filtered column indices
     real(r8) :: dtime                              !land model time step (sec)
     real(r8) :: heatr                              !energy residual or loss after melting or freezing
     real(r8) :: temp1                              !temporary variables [kg/m2]
     real(r8) :: hm(bounds%begc:bounds%endc)                        !energy residual [W/m2]
     real(r8) :: xm(bounds%begc:bounds%endc)                        !melting or freezing within a time step [kg/m2]
     real(r8) :: tinc                               !t(n+1)-t(n) (K)
     real(r8) :: smp                                !frozen water potential (mm)
     real(r8) :: rho_avg
     real(r8) :: z_avg
     real(r8) :: dcv(bounds%begc:bounds%endc) 
     real(r8) :: t_h2osfc_new
     real(r8) :: c1
     real(r8) :: c2
     real(r8) :: h_excess
     real(r8) :: c_h2osfc_ice
     !-----------------------------------------------------------------------

    call t_startf( 'PhaseChangeH2osfc' )

    ! Enforce expected array sizes
    SHR_ASSERT_ALL((ubound(fact)       == (/bounds%endc, nlevgrnd/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(dhsdT)      == (/bounds%endc/)),           errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(c_h2osfc)   == (/bounds%endc/)),           errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(xmf_h2osfc) == (/bounds%endc/)),           errMsg(__FILE__, __LINE__))

   associate(& 
   frac_sno                  =>    cps%frac_sno_eff        , & ! Input:  [real(r8) (:)]  fraction of ground covered by snow (0 to 1)
   frac_h2osfc               =>    cps%frac_h2osfc         , & ! Input:  [real(r8) (:)]  fraction of ground covered by surface water (0 to 1)
   t_h2osfc                  =>    ces%t_h2osfc            , & ! Input:  [real(r8) (:)]  surface water temperature               
   t_h2osfc_bef              =>    ces%t_h2osfc_bef        , & ! Input:  [real(r8) (:)]  saved surface water temperature         
   h2osfc                    =>    cws%h2osfc              , & ! Input:  [real(r8) (:)]  surface water (mm)                      
   int_snow                  =>    cws%int_snow            , & ! Input:  [real(r8) (:)]  integrated snowfall [mm]                
   qflx_h2osfc_to_ice        =>    cwf%qflx_h2osfc_to_ice  , & ! Input:  [real(r8) (:)]  conversion of h2osfc to ice             
   snl                       =>    cps%snl                 , & ! Input:  [integer (:)]  number of snow layers                    
   h2osno                    =>    cws%h2osno              , & ! Input:  [real(r8) (:)]  snow water (mm H2O)                     
   snow_depth                =>    cps%snow_depth          , & ! Input:  [real(r8) (:)] snow height (m)                          
   h2osoi_ice                =>    cws%h2osoi_ice          , & ! Input:  [real(r8) (:,:)] ice lens (kg/m2) (new)                 
   t_soisno                  =>    ces%t_soisno            , & ! Input:  [real(r8) (:,:)] soil temperature (Kelvin)              
   tssbef                    =>    ces%tssbef              , & ! Input:  [real(r8) (:,:)] temperature at previous time step [K]  
   dz                        =>    cps%dz                    & ! Input:  [real(r8) (:,:)] layer thickness (m)                    
   )

    ! Get step size

    dtime = get_step_size()

    ! Initialization

    do fc = 1,num_nolakec
       c = filter_nolakec(fc)
       xmf_h2osfc(c) = 0._r8
       hm(c) = 0._r8
       xm(c) = 0._r8
       qflx_h2osfc_to_ice(c) = 0._r8
    end do

    ! Freezing identification
    do fc = 1,num_nolakec
       c = filter_nolakec(fc)
       ! If liquid exists below melt point, freeze some to ice.
       if ( frac_h2osfc(c) > 0._r8 .AND. t_h2osfc(c) <= tfrz) then
          tinc = t_h2osfc(c)-tfrz
          t_h2osfc(c) = tfrz
          ! energy absorbed beyond freezing temperature
          hm(c) = dhsdT(c)*tinc - tinc*c_h2osfc(c)/dtime

          ! mass of water converted from liquid to ice
          xm(c) = hm(c)*dtime/hfus  
          temp1 = h2osfc(c) - xm(c)    

          ! compute change in cv due to additional ice
          dcv(c)=cpice*min(xm(c),h2osfc(c))

          z_avg=frac_sno(c)*snow_depth(c)
          if (z_avg > 0._r8) then 
             rho_avg=min(800._r8,h2osno(c)/z_avg)
          else
             rho_avg=200._r8
          endif
!=====================  xm < h2osfc  ====================================
          if(temp1 >= 0._r8) then ! add some frozen water to snow column
             ! add ice to snow column
             h2osno(c) = h2osno(c) + xm(c)
             int_snow(c) = int_snow(c) + xm(c)

             if(snl(c) < 0) h2osoi_ice(c,0) = h2osoi_ice(c,0) + xm(c)

             ! remove ice from h2osfc
             h2osfc(c) = h2osfc(c) - xm(c)
             
             xmf_h2osfc(c) = -frac_h2osfc(c)*hm(c)

             qflx_h2osfc_to_ice(c) = xm(c)/dtime

             ! update snow depth
             if (frac_sno(c) > 0 .and. snl(c) < 0) then 
                snow_depth(c)=h2osno(c)/(rho_avg*frac_sno(c))
             else
                snow_depth(c)=h2osno(c)/denice
             endif
!=========================  xm > h2osfc  =============================
          else !all h2osfc converted to ice, apply residual heat to top soil layer

             rho_avg=(h2osno(c)*rho_avg + h2osfc(c)*denice)/(h2osno(c) + h2osfc(c))
             h2osno(c) = h2osno(c) + h2osfc(c)
             int_snow(c) = int_snow(c) + h2osfc(c)

             qflx_h2osfc_to_ice(c) = h2osfc(c)/dtime

             ! excess energy is used to cool ice layer
             if(snl(c) < 0) h2osoi_ice(c,0) = h2osoi_ice(c,0) + h2osfc(c)

             ! compute heat capacity of frozen h2osfc layer
             c_h2osfc_ice=cpice*denice*(1.0e-3*h2osfc(c)) !h2osfc in [m]

             ! cool frozen h2osfc layer with extra heat
             t_h2osfc_new = t_h2osfc(c) - temp1*hfus/(dtime*dhsdT(c) - c_h2osfc_ice)

             ! next, determine equilibrium temperature of combined ice/snow layer
             xmf_h2osfc(c) = -frac_h2osfc(c)*hm(c)
             if (snl(c) == 0) then
                t_soisno(c,0) = t_h2osfc_new
             else if (snl(c) == -1) then
                c1=frac_sno(c)/fact(c,0) - dhsdT(c)*dtime
                if ( frac_h2osfc(c) /= 0.0_r8 )then
                   c2=frac_h2osfc(c)*(c_h2osfc_ice/dtime)
                else
                   c2=0.0_r8
                end if
                ! account for the change in t_soisno(c,0) via xmf_h2osfc(c)
                xmf_h2osfc(c) = xmf_h2osfc(c) + frac_sno(c)*t_soisno(c,0)/fact(c,0)
                t_soisno(c,0) = (c1*t_soisno(c,0)+ c2*t_h2osfc_new) &
                                   /(c1 + c2)             
                xmf_h2osfc(c) = xmf_h2osfc(c) - frac_sno(c)*t_soisno(c,0)/fact(c,0)

             else
                c1=frac_sno(c)/fact(c,0)
                if ( frac_h2osfc(c) /= 0.0_r8 )then
                   c2=frac_h2osfc(c)*(c_h2osfc_ice/dtime)
                else
                   c2=0.0_r8
                end if
                xmf_h2osfc(c) = xmf_h2osfc(c) + c1*t_soisno(c,0)
                t_soisno(c,0) = (c1*t_soisno(c,0)+ c2*t_h2osfc_new) &
                               /(c1 + c2)             
                xmf_h2osfc(c) = xmf_h2osfc(c) - c1*t_soisno(c,0)
             endif

             ! set h2osfc to zero (all liquid converted to ice)
             h2osfc(c) = 0._r8

             ! update snow depth
             if (frac_sno(c) > 0 .and. snl(c) < 0) then 
                snow_depth(c)=h2osno(c)/(rho_avg*frac_sno(c))
             else
                snow_depth(c)=h2osno(c)/denice
             endif

          endif
       endif
    enddo
    call t_stopf( 'PhaseChangeH2osfc' )
    end associate 
   end subroutine PhaseChangeH2osfc

   !-----------------------------------------------------------------------
   subroutine Phasechange_beta (bounds, num_nolakec, filter_nolakec, fact, dhsdT, xmf)
     !
     ! !DESCRIPTION:
     ! Calculation of the phase change within snow and soil layers:
     ! (1) Check the conditions for which the phase change may take place,
     !     i.e., the layer temperature is great than the freezing point
     !     and the ice mass is not equal to zero (i.e. melting),
     !     or the layer temperature is less than the freezing point
     !     and the liquid water mass is greater than the allowable supercooled 
     !     liquid water calculated from freezing point depression (i.e. freezing).
     ! (2) Assess the rate of phase change from the energy excess (or deficit)
     !     after setting the layer temperature to freezing point.
     ! (3) Re-adjust the ice and liquid mass, and the layer temperature
     !
     ! !USES:
     use clmtype
     use clm_time_manager, only : get_step_size
     use clm_varcon  , only : tfrz, hfus, grav
     use column_varcon,only : icol_roof, icol_sunwall, icol_shadewall, icol_road_perv
     use landunit_varcon,only:istsoil, istcrop, istice_mec
     use clm_varpar  , only : nlevsno, nlevgrnd,nlevurb
     use clm_varctl  , only : iulog
     !
     ! !ARGUMENTS:
     implicit none
     type(bounds_type), intent(in) :: bounds                      ! bounds
     integer , intent(in) :: num_nolakec                          ! number of column non-lake points in column filter
     integer , intent(in) :: filter_nolakec(:)                    ! column filter for non-lake points
     real(r8), intent(in) :: fact  ( bounds%begc: , -nlevsno+1: ) ! temporary [col, lev]
     real(r8), intent(in) :: dhsdT ( bounds%begc: )               ! temperature derivative of "hs" [col]
     real(r8), intent(out):: xmf   ( bounds%begc: )               ! total latent heat of phase change [col]
     !
     ! !LOCAL VARIABLES:
     integer  :: j,c,g,l                            !do loop index
     integer  :: fc                                 !lake filtered column indices
     real(r8) :: dtime                              !land model time step (sec)
     real(r8) :: heatr                              !energy residual or loss after melting or freezing
     real(r8) :: temp1                              !temporary variables [kg/m2]
     real(r8) :: hm(bounds%begc:bounds%endc,-nlevsno+1:nlevgrnd)    !energy residual [W/m2]
     real(r8) :: xm(bounds%begc:bounds%endc,-nlevsno+1:nlevgrnd)    !melting or freezing within a time step [kg/m2]
     real(r8) :: wmass0(bounds%begc:bounds%endc,-nlevsno+1:nlevgrnd)!initial mass of ice and liquid (kg/m2)
     real(r8) :: wice0 (bounds%begc:bounds%endc,-nlevsno+1:nlevgrnd)!initial mass of ice (kg/m2)
     real(r8) :: wliq0 (bounds%begc:bounds%endc,-nlevsno+1:nlevgrnd)!initial mass of liquid (kg/m2)
     real(r8) :: supercool(bounds%begc:bounds%endc,nlevgrnd)        !supercooled water in soil (kg/m2) 
     real(r8) :: propor                             !proportionality constant (-)
     real(r8) :: tinc(bounds%begc:bounds%endc,-nlevsno+1:nlevgrnd)  !t(n+1)-t(n) (K)
     real(r8) :: smp                                !frozen water potential (mm)
     !-----------------------------------------------------------------------

    call t_startf( 'PhaseChangebeta' )

    ! Enforce expected array sizes
    SHR_ASSERT_ALL((ubound(fact)  == (/bounds%endc, nlevgrnd/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(dhsdT) == (/bounds%endc/)),           errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(xmf)   == (/bounds%endc/)),           errMsg(__FILE__, __LINE__))

   associate(& 
   qflx_snow_melt            =>    cwf%qflx_snow_melt      , & ! Input:  [real(r8) (:)]  net snow melt                           
   frac_sno_eff              =>    cps%frac_sno_eff        , & ! Input:  [real(r8) (:)]  eff. fraction of ground covered by snow (0 to 1)
   frac_sno                  =>    cps%frac_sno            , & ! Input:  [real(r8) (:)]  fraction of ground covered by snow (0 to 1)
   frac_h2osfc               =>    cps%frac_h2osfc         , & ! Input:  [real(r8) (:)]  fraction of ground covered by surface water (0 to 1)
   ltype                     =>    lun%itype               , & ! Input:  [integer (:)]  landunit type                            
   urbpoi                    =>    lun%urbpoi              , & ! Input:  [logical (:)]  true => landunit is an urban point       
   ctype                     =>    col%itype               , & ! Input:  [integer (:)] column type                               
   snl                       =>    cps%snl                 , & ! Input:  [integer (:)]  number of snow layers                    
   h2osno                    =>    cws%h2osno              , & ! Input:  [real(r8) (:)]  snow water (mm H2O)                     
   snow_depth                =>    cps%snow_depth          , & ! Input:  [real(r8) (:)]  snow height (m)                         
   qflx_snomelt              =>    cwf%qflx_snomelt        , & ! Output: [real(r8) (:)]  snow melt (mm H2O /s)                   
   eflx_snomelt              =>    cef%eflx_snomelt        , & ! Output: [real(r8) (:)] snow melt heat flux (W/m**2)             
   eflx_snomelt_u            =>    cef%eflx_snomelt_u      , & ! Output: [real(r8) (:)] urban snow melt heat flux (W/m**2)       
   eflx_snomelt_r            =>    cef%eflx_snomelt_r      , & ! Output: [real(r8) (:)] rural snow melt heat flux (W/m**2)       
   h2osoi_liq                =>    cws%h2osoi_liq          , & ! Input:  [real(r8) (:,:)] liquid water (kg/m2) (new)             
   h2osoi_ice                =>    cws%h2osoi_ice          , & ! Input:  [real(r8) (:,:)] ice lens (kg/m2) (new)                 
   imelt                     =>    cps%imelt               , & ! Output: [integer (:,:)] flag for melting (=1), freezing (=2), Not=0 (new)
   t_soisno                  =>    ces%t_soisno            , & ! Input:  [real(r8) (:,:)] soil temperature (Kelvin)              
   tssbef                    =>    ces%tssbef              , & ! Input:  [real(r8) (:,:)] temperature at previous time step [K]  
   bsw                       =>    cps%bsw                 , & ! Input:  [real(r8) (:,:)] Clapp and Hornberger "b"               
   sucsat                    =>    cps%sucsat              , & ! Input:  [real(r8) (:,:)] minimum soil suction (mm)              
   watsat                    =>    cps%watsat              , & ! Input:  [real(r8) (:,:)] volumetric soil water at saturation (porosity)
   dz                        =>    cps%dz                  , & ! Input:  [real(r8) (:,:)] layer thickness (m)                    
   qflx_snofrz_lyr           =>    cwf%qflx_snofrz_lyr     , & ! Output: [real(r8) (:,:)] snow freezing rate (positive definite) (col,lyr) [kg m-2 s-1]
   qflx_snofrz_col           =>    cwf%qflx_snofrz_col     , & ! Output: [real(r8) (:)] column-integrated snow freezing rate (positive definite) [kg m-2 s-1]
   qflx_glcice               =>    cwf%qflx_glcice         , & ! Output: [real(r8) (:)] flux of new glacier ice (mm H2O/s) [+ = ice grows]
   qflx_glcice_melt          =>    cwf%qflx_glcice_melt      & ! Output: [real(r8) (:)] ice melt (positive definite) (mm H2O/s)  
   )

    ! Get step size

    dtime = get_step_size()

    ! Initialization

    do fc = 1,num_nolakec
       c = filter_nolakec(fc)
       l = col%landunit(c)

       qflx_snomelt(c) = 0._r8
       xmf(c) = 0._r8
       qflx_snofrz_lyr(c,-nlevsno+1:0) = 0._r8
       qflx_snofrz_col(c) = 0._r8
       qflx_glcice_melt(c) = 0._r8
       qflx_snow_melt(c) = 0._r8
    end do

    do j = -nlevsno+1,nlevgrnd       ! all layers
       do fc = 1,num_nolakec
          c = filter_nolakec(fc)
          if (j >= snl(c)+1) then

             ! Initialization
             imelt(c,j) = 0
             hm(c,j) = 0._r8
             xm(c,j) = 0._r8
             wice0(c,j) = h2osoi_ice(c,j)
             wliq0(c,j) = h2osoi_liq(c,j)
             wmass0(c,j) = h2osoi_ice(c,j) + h2osoi_liq(c,j)
          endif   ! end of snow layer if-block
       end do   ! end of column-loop
    enddo   ! end of level-loop

!--  snow layers  --------------------------------------------------- 
    do j = -nlevsno+1,0             
       do fc = 1,num_nolakec
          c = filter_nolakec(fc)
          if (j >= snl(c)+1) then

             ! Melting identification
             ! If ice exists above melt point, melt some to liquid.
             if (h2osoi_ice(c,j) > 0._r8 .AND. t_soisno(c,j) > tfrz) then
                imelt(c,j) = 1
!                tinc(c,j) = t_soisno(c,j) - tfrz 
                tinc(c,j) = tfrz - t_soisno(c,j) 
                t_soisno(c,j) = tfrz
             endif

             ! Freezing identification
             ! If liquid exists below melt point, freeze some to ice.
             if (h2osoi_liq(c,j) > 0._r8 .AND. t_soisno(c,j) < tfrz) then
                imelt(c,j) = 2
!                tinc(c,j) = t_soisno(c,j) - tfrz 
                tinc(c,j) = tfrz - t_soisno(c,j) 
                t_soisno(c,j) = tfrz
             endif
          endif   ! end of snow layer if-block
       end do   ! end of column-loop
    enddo   ! end of level-loop

!-- soil layers   --------------------------------------------------- 
    do j = 1,nlevgrnd             
       do fc = 1,num_nolakec
          c = filter_nolakec(fc)
          l = col%landunit(c)
          supercool(c,j) = 0.0_r8
          ! add in urban condition if-block
          if ((ctype(c) /= icol_sunwall .and. ctype(c) /= icol_shadewall &
               .and. ctype(c) /= icol_roof) .or. ( j <= nlevurb)) then



          if (h2osoi_ice(c,j) > 0. .AND. t_soisno(c,j) > tfrz) then
             imelt(c,j) = 1
!             tinc(c,j) = t_soisno(c,j) - tfrz 
             tinc(c,j) = tfrz - t_soisno(c,j) 
             t_soisno(c,j) = tfrz
          endif

          ! from Zhao (1997) and Koren (1999)
          supercool(c,j) = 0.0_r8
          if (ltype(l) == istsoil .or. ltype(l) == istcrop .or. ctype(c) == icol_road_perv) then
             if(t_soisno(c,j) < tfrz) then
                smp = hfus*(tfrz-t_soisno(c,j))/(grav*t_soisno(c,j)) * 1000._r8  !(mm)
                supercool(c,j) = watsat(c,j)*(smp/sucsat(c,j))**(-1._r8/bsw(c,j))
                supercool(c,j) = supercool(c,j)*dz(c,j)*1000._r8       ! (mm)
             endif
          endif

          if (h2osoi_liq(c,j) > supercool(c,j) .AND. t_soisno(c,j) < tfrz) then
             imelt(c,j) = 2
!             tinc(c,j) = t_soisno(c,j) - tfrz
             tinc(c,j) = tfrz - t_soisno(c,j) 
             t_soisno(c,j) = tfrz
          endif

          ! If snow exists, but its thickness is less than the critical value (0.01 m)
          if (snl(c)+1 == 1 .AND. h2osno(c) > 0._r8 .AND. j == 1) then
             if (t_soisno(c,j) > tfrz) then
                imelt(c,j) = 1
!                tincc,j) = t_soisno(c,j) - tfrz
                tinc(c,j) = tfrz - t_soisno(c,j)
                t_soisno(c,j) = tfrz
             endif
          endif

       endif

       end do
    enddo


    do j = -nlevsno+1,nlevgrnd       ! all layers
       do fc = 1,num_nolakec
          c = filter_nolakec(fc)

          if ((ctype(c) /= icol_sunwall .and. ctype(c) /= icol_shadewall &
          .and. ctype(c) /= icol_roof) .or. ( j <= nlevurb)) then

             if (j >= snl(c)+1) then

                ! Calculate the energy surplus and loss for melting and freezing
                if (imelt(c,j) > 0) then
              
                   ! added unique cases for this calculation,
                   ! to account for absorbed solar radiation in each layer

                   !==================================================================
                   if (j == snl(c)+1) then ! top layer                   
                      hm(c,j) = dhsdT(c)*tinc(c,j) - tinc(c,j)/fact(c,j)

                      if ( j==1 .and. frac_h2osfc(c) /= 0.0_r8 ) then
                         hm(c,j) = hm(c,j) - frac_h2osfc(c)*(dhsdT(c)*tinc(c,j))
                      end if
                   else if (j == 1) then
                      hm(c,j) = (1.0_r8 - frac_sno_eff(c) - frac_h2osfc(c)) &
                           *dhsdT(c)*tinc(c,j) - tinc(c,j)/fact(c,j)
                   else ! non-interfacial snow/soil layers                   
                      hm(c,j) = - tinc(c,j)/fact(c,j)
                   endif
                endif

                ! These two errors were checked carefully (Y. Dai).  They result from the
                ! computed error of "Tridiagonal-Matrix" in subroutine "thermal".
                if (imelt(c,j) == 1 .AND. hm(c,j) < 0._r8) then
                   hm(c,j) = 0._r8
                   imelt(c,j) = 0
                endif
                if (imelt(c,j) == 2 .AND. hm(c,j) > 0._r8) then
                   hm(c,j) = 0._r8
                   imelt(c,j) = 0
                endif

                ! The rate of melting and freezing
   
                if (imelt(c,j) > 0 .and. abs(hm(c,j)) > 0._r8) then
                   xm(c,j) = hm(c,j)*dtime/hfus                           ! kg/m2
   
                   ! If snow exists, but its thickness is less than the critical value
                   ! (1 cm). Note: more work is needed to determine how to tune the
                   ! snow depth for this case
                   if (j == 1) then
                      if (snl(c)+1 == 1 .AND. h2osno(c) > 0._r8 .AND. xm(c,j) > 0._r8) then
                         temp1 = h2osno(c)                           ! kg/m2
                         h2osno(c) = max(0._r8,temp1-xm(c,j))
                         propor = h2osno(c)/temp1
                         snow_depth(c) = propor * snow_depth(c)
                         heatr = hm(c,j) - hfus*(temp1-h2osno(c))/dtime   ! W/m2
                         if (heatr > 0._r8) then
                            xm(c,j) = heatr*dtime/hfus                    ! kg/m2
                            hm(c,j) = heatr                               ! W/m2
                         else
                            xm(c,j) = 0._r8
                            hm(c,j) = 0._r8
                         endif
                         qflx_snomelt(c) = max(0._r8,(temp1-h2osno(c)))/dtime   ! kg/(m2 s)
                         xmf(c) = hfus*qflx_snomelt(c)
                         qflx_snow_melt(c) = qflx_snomelt(c) 
                      endif
                   endif

                   heatr = 0._r8
                   if (xm(c,j) > 0._r8) then
                      h2osoi_ice(c,j) = max(0._r8, wice0(c,j)-xm(c,j))
                      heatr = hm(c,j) - hfus*(wice0(c,j)-h2osoi_ice(c,j))/dtime
                   else if (xm(c,j) < 0._r8) then
                      if (j <= 0) then
                         h2osoi_ice(c,j) = min(wmass0(c,j), wice0(c,j)-xm(c,j))  ! snow
                      else
                         if (wmass0(c,j) < supercool(c,j)) then
                            h2osoi_ice(c,j) = 0._r8
                         else
                            h2osoi_ice(c,j) = min(wmass0(c,j) - supercool(c,j),wice0(c,j)-xm(c,j))
                         endif
                      endif
                      heatr = hm(c,j) - hfus*(wice0(c,j)-h2osoi_ice(c,j))/dtime
                   endif
   
                   h2osoi_liq(c,j) = max(0._r8,wmass0(c,j)-h2osoi_ice(c,j))

                   if (abs(heatr) > 0._r8) then
                      if (j == snl(c)+1) then

                         if(j==1) then
                            t_soisno(c,j) = t_soisno(c,j) + fact(c,j)*heatr &
                                 /(1._r8-(1.0_r8 - frac_h2osfc(c))*fact(c,j)*dhsdT(c))
                         else
                            t_soisno(c,j) = t_soisno(c,j) + fact(c,j)*heatr &
                                 /(1._r8-fact(c,j)*dhsdT(c))
                         endif

                      else if (j == 1) then
   
                         t_soisno(c,j) = t_soisno(c,j) + fact(c,j)*heatr &
                              /(1._r8-(1.0_r8 - frac_sno_eff(c) - frac_h2osfc(c))*fact(c,j)*dhsdT(c))
                      else
                         t_soisno(c,j) = t_soisno(c,j) + fact(c,j)*heatr
                      endif

                      if (j <= 0) then    ! snow
                         if (h2osoi_liq(c,j)*h2osoi_ice(c,j)>0._r8) t_soisno(c,j) = tfrz
                      end if
                   endif  ! end of heatr > 0 if-block

                   if (j >= 1) then 
                      xmf(c) = xmf(c) + hfus*(wice0(c,j)-h2osoi_ice(c,j))/dtime
                   else
                      xmf(c) = xmf(c) + frac_sno_eff(c)*hfus*(wice0(c,j)-h2osoi_ice(c,j))/dtime
                   endif

                   if (imelt(c,j) == 1 .AND. j < 1) then
                      qflx_snomelt(c) = qflx_snomelt(c) + max(0._r8,(wice0(c,j)-h2osoi_ice(c,j)))/dtime


                   endif

                   ! layer freezing mass flux (positive):
                   if (imelt(c,j) == 2 .AND. j < 1) then
                      qflx_snofrz_lyr(c,j) = max(0._r8,(h2osoi_ice(c,j)-wice0(c,j)))/dtime
                   endif

                endif

             endif   ! end of snow layer if-block

          endif

          ! For glacier_mec columns, compute negative ice flux from melted ice.
          ! Note that qflx_glcice can also include a positive component from excess snow,
          !  as computed in Hydrology2Mod.F90.  
          
          ! Note also here: There is no need to calculate a negative ice flux for
          ! qflx_glcice here over soil columns (analagous to the + ice flux potentially
          ! calculated in Hydrology2Mod.F90), since by definition there is no ice to melt
          ! over bare land.

          l = col%landunit(c)
          if (ltype(l)==istice_mec) then

             if (j>=1 .and. h2osoi_liq(c,j) > 0._r8) then   ! ice layer with meltwater
                ! melting corresponds to a negative ice flux
                qflx_glcice_melt(c) = qflx_glcice_melt(c) + h2osoi_liq(c,j)/dtime
                qflx_glcice(c) =      qflx_glcice(c)      - h2osoi_liq(c,j)/dtime

                ! convert layer back to pure ice by "borrowing" ice from below the column
                h2osoi_ice(c,j) = h2osoi_ice(c,j) + h2osoi_liq(c,j)
                h2osoi_liq(c,j) = 0._r8

             endif  ! liquid water is present
          endif     ! istice_mec

       end do   ! end of column-loop
    enddo   ! end of level-loop

    ! Needed for history file output

    do fc = 1,num_nolakec
       c = filter_nolakec(fc)
       eflx_snomelt(c) = qflx_snomelt(c) * hfus
       l = col%landunit(c)
       if (urbpoi(l)) then
         eflx_snomelt_u(c) = eflx_snomelt(c)
       else if (ltype(l) == istsoil .or. ltype(l) == istcrop) then
         eflx_snomelt_r(c) = eflx_snomelt(c)
       end if
    end do

    call t_stopf( 'PhaseChangebeta' )
    do j = -nlevsno+1,0
       do fc = 1,num_nolakec
          c = filter_nolakec(fc)
          qflx_snofrz_col(c) = qflx_snofrz_col(c) + qflx_snofrz_lyr(c,j)
       end do
    end do

    end associate 
   end subroutine Phasechange_beta

   !-----------------------------------------------------------------------
   subroutine ComputeGroundHeatFluxAndDeriv(bounds, num_nolakec, filter_nolakec, &
        hs_h2osfc, hs_top_snow, hs_soil, hs_top, dhsdT, sabg_lyr_col)

    !
    ! !DESCRIPTION:
    ! Computes ground heat flux on:
    ! (1) The surface of standing water,
    ! (2) The surface of snow,
    ! (3) The surface of soil, and
    ! (4) Net energy flux into soil surface.
    ! Additionally, derivative of ground heat flux w.r.t to temeprature
    !
    ! !USES:
    use clmtype
    use clm_atmlnd      , only : a2l_downscaled_col
    use shr_infnan_mod  , only : nan => shr_infnan_nan, assignment(=)
    use clm_varcon      , only : sb, hvap
    use column_varcon  , only : icol_road_perv, icol_road_imperv
    use clm_varpar     , only : nlevsno, max_pft_per_col
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds                            ! bounds
    integer , intent(in)  :: num_nolakec                               ! number of column non-lake points in column filter
    integer , intent(in)  :: filter_nolakec( : )                       ! column filter for non-lake points
    real(r8), intent(out) :: hs_h2osfc( bounds%begc: )                 ! heat flux on standing water [W/m2]
    real(r8), intent(out) :: hs_top_snow( bounds%begc: )               ! heat flux on top snow layer [W/m2]
    real(r8), intent(out) :: hs_soil( bounds%begc: )                   ! heat flux on soil [W/m2]
    real(r8), intent(out) :: hs_top (bounds%begc: )                    ! net energy flux into surface layer (col) [W/m2]
    real(r8), intent(out) :: dhsdT( bounds%begc: )                     ! temperature derivative of "hs" [col]
    real(r8), intent(out) :: sabg_lyr_col( bounds%begc:, -nlevsno+1: ) ! absorbed solar radiation (col,lyr) [W/m2]
    !
    ! !LOCAL VARIABLES:
    integer  :: j,c,p,l,g,pi                                           ! indices
    integer  :: fc                                                     ! lake filtered column indices
    real(r8) :: hs(bounds%begc:bounds%endc)                            ! net energy flux into the surface (w/m2)
    real(r8) :: lwrad_emit(bounds%begc:bounds%endc)                    ! emitted longwave radiation
    real(r8) :: dlwrad_emit(bounds%begc:bounds%endc)                   ! time derivative of emitted longwave radiation
    integer  :: lyr_top                                                ! index of top layer of snowpack (-4 to 0) [idx]
    real(r8) :: eflx_gnet_top                                          ! net energy flux into surface layer, pft-level [W/m2]
    real(r8) :: lwrad_emit_snow(bounds%begc:bounds%endc)               !
    real(r8) :: lwrad_emit_soil(bounds%begc:bounds%endc)               !
    real(r8) :: lwrad_emit_h2osfc(bounds%begc:bounds%endc)             !
    real(r8) :: eflx_gnet_snow                                         !
    real(r8) :: eflx_gnet_soil                                         !
    real(r8) :: eflx_gnet_h2osfc                                       !
    !-----------------------------------------------------------------------

    ! Enforce expected array sizes
    SHR_ASSERT_ALL((ubound(hs_h2osfc)     == (/bounds%endc/)),   errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(hs_top_snow)   == (/bounds%endc/)),   errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(hs_soil)       == (/bounds%endc/)),   errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(hs_top)        == (/bounds%endc/)),   errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(dhsdT)         == (/bounds%endc/)),   errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(sabg_lyr_col)  == (/bounds%endc,1/)), errMsg(__FILE__, __LINE__))

   associate(&
   forc_lwrad                =>    a2l_downscaled_col%forc_lwrad , & ! Input:  [real(r8) (:)]  downward infrared (longwave) radiation (W/m**2)
   eflx_traffic              =>    lef%eflx_traffic              , & ! Input:  [real(r8) (:)]  traffic sensible heat flux (W/m**2)
   eflx_wasteheat            =>    lef%eflx_wasteheat            , & ! Input:  [real(r8) (:)]  sensible heat flux from urban heating/cooling sources of waste heat (W/m**2)
   eflx_heat_from_ac         =>    lef%eflx_heat_from_ac         , & ! Input:  [real(r8) (:)] sensible heat flux put back into canyon due to removal by AC (W/m**2)
   wtlunit_roof              =>    lun%wtlunit_roof              , & ! Input:  [real(r8) (:)]  weight of roof with respect to landunit
   urbpoi                    =>    lun%urbpoi                    , & ! Input:  [logical (:)]  true => landunit is an urban point
   frac_sno_eff              =>    cps%frac_sno_eff              , & ! Input:  [real(r8) (:)]  eff. fraction of ground covered by snow (0 to 1)
   t_h2osfc                  =>    ces%t_h2osfc                  , & ! Input:  [real(r8) (:)]  surface water temperature
   eflx_sh_snow              =>    pef%eflx_sh_snow              , & ! Input:  [real(r8) (:)]  sensible heat flux from snow (W/m**2) [+ to atm]
   eflx_sh_soil              =>    pef%eflx_sh_soil              , & ! Input:  [real(r8) (:)]  sensible heat flux from soil (W/m**2) [+ to atm]
   eflx_sh_h2osfc            =>    pef%eflx_sh_h2osfc            , & ! Input:  [real(r8) (:)]  sensible heat flux from surface water (W/m**2) [+ to atm]
   qflx_ev_snow              =>    pwf%qflx_ev_snow              , & ! Input:  [real(r8) (:)]  evaporation flux from snow (W/m**2) [+ to atm]
   qflx_ev_soil              =>    pwf%qflx_ev_soil              , & ! Input:  [real(r8) (:)]  evaporation flux from soil (W/m**2) [+ to atm]
   qflx_ev_h2osfc            =>    pwf%qflx_ev_h2osfc            , & ! Input:  [real(r8) (:)]  evaporation flux from h2osfc (W/m**2) [+ to atm]
   sabg_soil                 =>    pef%sabg_soil                 , & ! Input:  [real(r8) (:)]  solar radiation absorbed by soil (W/m**2)
   sabg_snow                 =>    pef%sabg_snow                 , & ! Input:  [real(r8) (:)]  solar radiation absorbed by snow (W/m**2)
   sabg_chk                  =>    pef%sabg_chk                  , & ! Input:  [real(r8) (:)]  sum of soil/snow using current fsno, for balance check
   ctype                     =>    col%itype                     , & ! Input:  [integer (:)]  column type
   npfts                     =>    col%npfts                     , & ! Input:  [integer (:)]  column's number of pfts
   pfti                      =>    col%pfti                      , & ! Input:  [integer (:)]  column's beginning pft index
   snl                       =>    cps%snl                       , & ! Input:  [integer (:)]  number of snow layers
   htvp                      =>    cps%htvp                      , & ! Input:  [real(r8) (:)]  latent heat of vapor of water (or sublimation) [j/kg]
   emg                       =>    cps%emg                       , & ! Input:  [real(r8) (:)]  ground emissivity
   t_grnd                    =>    ces%t_grnd                    , & ! Input:  [real(r8) (:)]  ground surface temperature [K]
   zi                        =>    cps%zi                        , & ! Input:  [real(r8) (:,:)]  interface level below a "z" level (m)
   dz                        =>    cps%dz                        , & ! Input:  [real(r8) (:,:)]  layer depth (m)
   z                         =>    cps%z                         , & ! Input:  [real(r8) (:,:)]  layer thickness (m)
   t_soisno                  =>    ces%t_soisno                  , & ! Input:  [real(r8) (:,:)]  soil temperature (Kelvin)
   pactive                   =>    pft%active                    , & ! Input:  [logical (:)]  true=>do computations on this pft
   pgridcell                 =>    pft%gridcell                  , & ! Input:  [integer (:)]  pft's gridcell index
   plandunit                 =>    pft%landunit                  , & ! Input:  [integer (:)]  pft's landunit index
   pwtcol                    =>    pft%wtcol                     , & ! Input:  [real(r8) (:)]  weight of pft relative to column
   frac_veg_nosno            =>    pps%frac_veg_nosno            , & ! Input:  [integer (:)]  fraction of vegetation not covered by snow (0 OR 1 now) [-] (new)
   cgrnd                     =>    pef%cgrnd                     , & ! Input:  [real(r8) (:)]  deriv. of soil energy flux wrt to soil temp [w/m2/k]
   dlrad                     =>    pef%dlrad                     , & ! Input:  [real(r8) (:)]  downward longwave radiation blow the canopy [W/m2]
   sabg                      =>    pef%sabg                      , & ! Input:  [real(r8) (:)]  solar radiation absorbed by ground (W/m**2)
   eflx_sh_grnd              =>    pef%eflx_sh_grnd              , & ! Input:  [real(r8) (:)]  sensible heat flux from ground (W/m**2) [+ to atm]
   qflx_evap_soi             =>    pwf%qflx_evap_soi             , & ! Input:  [real(r8) (:)]  soil evaporation (mm H2O/s) (+ = to atm)
   qflx_tran_veg             =>    pwf%qflx_tran_veg             , & ! Input:  [real(r8) (:)]  vegetation transpiration (mm H2O/s) (+ = to atm)
   eflx_gnet                 =>    pef%eflx_gnet                 , & ! Output: [real(r8) (:)]  net ground heat flux into the surface (W/m**2)
   dgnetdT                   =>    pef%dgnetdT                   , & ! Output: [real(r8) (:)]  temperature derivative of ground net heat flux
   eflx_lwrad_net            =>    pef%eflx_lwrad_net            , & ! Input:  [real(r8) (:)]  net infrared (longwave) rad (W/m**2) [+ = to atm]
   eflx_wasteheat_pft        =>    pef%eflx_wasteheat_pft        , & ! Input:  [real(r8) (:)]  sensible heat flux from urban heating/cooling sources of waste heat (W/m**2)
   eflx_heat_from_ac_pft     =>    pef%eflx_heat_from_ac_pft     , & ! Input:  [real(r8) (:)] sensible heat flux put back into canyon due to removal by AC (W/m**2)
   eflx_traffic_pft          =>    pef%eflx_traffic_pft          , & ! Input:  [real(r8) (:)]  traffic sensible heat flux (W/m**2)
   eflx_anthro               =>    pef%eflx_anthro               , & ! Input:  [real(r8) (:)]  total anthropogenic heat flux (W/m**2)
   sabg_lyr                  =>    pef%sabg_lyr                  , & ! Output: [real(r8) (:,:)]  absorbed solar radiation (pft,lyr) [W/m2]
   begc                      =>    bounds%begc                   , & ! Input:  [integer ] beginning column index
   endc                      =>    bounds%endc                     & ! Input:  [integer ] ending column index
   )

    ! Net ground heat flux into the surface and its temperature derivative
    ! Added a pfts loop here to get the average of hs and dhsdT over
    ! all PFTs on the column. Precalculate the terms that do not depend on PFT.

    do fc = 1,num_nolakec
       c = filter_nolakec(fc)
       lwrad_emit(c)  =    emg(c) * sb * t_grnd(c)**4
       dlwrad_emit(c) = 4._r8*emg(c) * sb * t_grnd(c)**3

       ! fractionate lwrad_emit; balanced in CanopyFluxes & Biogeophysics2
       lwrad_emit_snow(c)    =    emg(c) * sb * t_soisno(c,snl(c)+1)**4
       lwrad_emit_soil(c)    =    emg(c) * sb * t_soisno(c,1)**4
       lwrad_emit_h2osfc(c)  =    emg(c) * sb * t_h2osfc(c)**4
    end do

    hs_soil(begc:endc)   = 0._r8
    hs_h2osfc(begc:endc) = 0._r8
    hs(begc:endc)        = 0._r8
    dhsdT(begc:endc)     = 0._r8
    do pi = 1,max_pft_per_col
       do fc = 1,num_nolakec
          c = filter_nolakec(fc)
          if ( pi <= npfts(c) ) then
             p = pfti(c) + pi - 1
             l = plandunit(p)
             g = pgridcell(p)

             if (pactive(p)) then
                if (.not. urbpoi(l)) then
                   eflx_gnet(p) = sabg(p) + dlrad(p) &
                        + (1-frac_veg_nosno(p))*emg(c)*forc_lwrad(c) - lwrad_emit(c) &
                        - (eflx_sh_grnd(p)+qflx_evap_soi(p)*htvp(c))
                   ! save sabg for balancecheck, in case frac_sno is set to zero later
                   sabg_chk(p) = frac_sno_eff(c) * sabg_snow(p) + (1._r8 - frac_sno_eff(c) ) * sabg_soil(p)

                   eflx_gnet_snow = sabg_snow(p) + dlrad(p) &
                        + (1-frac_veg_nosno(p))*emg(c)*forc_lwrad(c) - lwrad_emit_snow(c) &
                        - (eflx_sh_snow(p)+qflx_ev_snow(p)*htvp(c))

                   eflx_gnet_soil = sabg_soil(p) + dlrad(p) &
                        + (1-frac_veg_nosno(p))*emg(c)*forc_lwrad(c) - lwrad_emit_soil(c) &
                        - (eflx_sh_soil(p)+qflx_ev_soil(p)*htvp(c))

                   eflx_gnet_h2osfc = sabg_soil(p) + dlrad(p) &
                        + (1-frac_veg_nosno(p))*emg(c)*forc_lwrad(c) - lwrad_emit_h2osfc(c) &
                        - (eflx_sh_h2osfc(p)+qflx_ev_h2osfc(p)*htvp(c))
                else
                   ! For urban columns we use the net longwave radiation (eflx_lwrad_net) because of
                   ! interactions between urban columns.

                   ! All wasteheat and traffic flux goes into canyon floor
                   if (ctype(c) == icol_road_perv .or. ctype(c) == icol_road_imperv) then
                      eflx_wasteheat_pft(p) = eflx_wasteheat(l)/(1._r8-wtlunit_roof(l))
                      eflx_heat_from_ac_pft(p) = eflx_heat_from_ac(l)/(1._r8-wtlunit_roof(l))
                      eflx_traffic_pft(p) = eflx_traffic(l)/(1._r8-wtlunit_roof(l))
                   else
                      eflx_wasteheat_pft(p) = 0._r8
                      eflx_heat_from_ac_pft(p) = 0._r8
                      eflx_traffic_pft(p) = 0._r8
                   end if
                   ! Include transpiration term because needed for previous road
                   ! and include wasteheat and traffic flux
                   eflx_gnet(p) = sabg(p) + dlrad(p)  &
                        - eflx_lwrad_net(p) &
                        - (eflx_sh_grnd(p) + qflx_evap_soi(p)*htvp(c) + qflx_tran_veg(p)*hvap) &
                        + eflx_wasteheat_pft(p) + eflx_heat_from_ac_pft(p) + eflx_traffic_pft(p)
                   eflx_anthro(p)   = eflx_wasteheat_pft(p) + eflx_traffic_pft(p)
                   eflx_gnet_snow   = eflx_gnet(p)
                   eflx_gnet_soil   = eflx_gnet(p)
                   eflx_gnet_h2osfc = eflx_gnet(p)
                end if
                dgnetdT(p) = - cgrnd(p) - dlwrad_emit(c)
                hs(c) = hs(c) + eflx_gnet(p) * pwtcol(p)
                dhsdT(c) = dhsdT(c) + dgnetdT(p) * pwtcol(p)
                ! separate surface fluxes for soil/snow
                hs_soil(c) = hs_soil(c) + eflx_gnet_soil * pwtcol(p)
                hs_h2osfc(c) = hs_h2osfc(c) + eflx_gnet_h2osfc * pwtcol(p)

             end if
          end if
       end do
    end do

    !       Additional calculations with SNICAR:
    !       Set up tridiagonal matrix in a new manner. There is now
    !       absorbed solar radiation in each snow layer, instead of
    !       only the surface. Following the current implementation,
    !       absorbed solar flux should be: S + ((delS/delT)*dT),
    !       where S is absorbed radiation, and T is temperature. Now,
    !       assume delS/delT is zero, then it is OK to just add S
    !       to each layer

    ! Initialize:
    sabg_lyr_col(begc:endc,-nlevsno+1:1) = 0._r8
    hs_top(begc:endc)                    = 0._r8
    hs_top_snow(begc:endc)               = 0._r8

    do pi = 1,max_pft_per_col
       do fc = 1,num_nolakec
          c = filter_nolakec(fc)
          lyr_top = snl(c) + 1
          if ( pi <= npfts(c) ) then
             p = pfti(c) + pi - 1
             if (pactive(p)) then
                g = pgridcell(p)
                l = plandunit(p)
                if (.not. urbpoi(l)) then

                   eflx_gnet_top = sabg_lyr(p,lyr_top) + dlrad(p) + (1-frac_veg_nosno(p))*emg(c)*forc_lwrad(c) &
                        - lwrad_emit(c) - (eflx_sh_grnd(p)+qflx_evap_soi(p)*htvp(c))

                   hs_top(c) = hs_top(c) + eflx_gnet_top*pwtcol(p)

                   eflx_gnet_snow = sabg_lyr(p,lyr_top) + dlrad(p) + (1-frac_veg_nosno(p))*emg(c)*forc_lwrad(c) &
                        - lwrad_emit_snow(c) - (eflx_sh_snow(p)+qflx_ev_snow(p)*htvp(c))

                   eflx_gnet_soil = sabg_lyr(p,lyr_top) + dlrad(p) + (1-frac_veg_nosno(p))*emg(c)*forc_lwrad(c) &
                        - lwrad_emit_soil(c) - (eflx_sh_soil(p)+qflx_ev_soil(p)*htvp(c))

                   hs_top_snow(c) = hs_top_snow(c) + eflx_gnet_snow*pwtcol(p)

                   do j = lyr_top,1,1
                      sabg_lyr_col(c,j) = sabg_lyr_col(c,j) + sabg_lyr(p,j) * pwtcol(p)
                   enddo
                else

                   hs_top(c)      = hs_top(c) + eflx_gnet(p)*pwtcol(p)
                   hs_top_snow(c) = hs_top_snow(c) + eflx_gnet(p)*pwtcol(p)
                   sabg_lyr_col(c,lyr_top) = sabg_lyr_col(c,lyr_top) + sabg(p) * pwtcol(p)

                endif
             endif

          endif
       enddo
    enddo

    end associate

  end subroutine ComputeGroundHeatFluxAndDeriv

  !-----------------------------------------------------------------------
  subroutine ComputeHeatDiffFluxAndFactor(bounds, num_nolakec, filter_nolakec, dtime, &
       tk, cv, fn, fact)

    !
    ! !DESCRIPTION:
    ! Computes:
    ! (1) Heat diffusion at the interface of layers.
    ! (2) Factor used in computing tridiagonal matrix
    !
    ! !USES:
    use clmtype
    use shr_infnan_mod  , only : nan => shr_infnan_nan, assignment(=)
    use clm_varcon      , only : capr, cnfac
    use column_varcon  , only : icol_roof, icol_sunwall, icol_shadewall
    use clm_varpar     , only : nlevsno, nlevgrnd, nlevurb
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds                     ! bounds
    integer , intent(in)  :: num_nolakec                        ! number of column non-lake points in column filter
    integer , intent(in)  :: filter_nolakec(:)                  ! column filter for non-lake points
    real(r8), intent(in)  :: dtime                              ! land model time step (sec)
    real(r8), intent(in)  :: tk(bounds%begc: ,-nlevsno+1: )     ! thermal conductivity [W/(m K)]
    real(r8), intent(in)  :: cv (bounds%begc: ,-nlevsno+1: )    ! heat capacity [J/(m2 K)]
    real(r8), intent(out) :: fn (bounds%begc: ,-nlevsno+1: )    ! heat diffusion through the layer interface [W/m2]
    real(r8), intent(out) :: fact( bounds%begc: , -nlevsno+1: ) ! used in computing tridiagonal matrix [col, lev]
    !
    ! !LOCAL VARIABLES:
    integer  :: j,c,l                                           ! indices
    integer  :: fc                                              ! lake filtered column indices
    real(r8) :: dzm                                             ! used in computing tridiagonal matrix
    !-----------------------------------------------------------------------

    ! Enforce expected array sizes
    SHR_ASSERT_ALL((ubound(tk)   == (/bounds%endc, nlevgrnd/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(cv)   == (/bounds%endc, nlevgrnd/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(fact) == (/bounds%endc, nlevgrnd/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(fn)   == (/bounds%endc, nlevgrnd/)), errMsg(__FILE__, __LINE__))

   associate(&
   t_building                =>    lps%t_building                , & ! Input:  [real(r8) (:)]  internal building temperature (K)
   ctype                     =>    col%itype                     , & ! Input:  [integer (:)]  column type
   clandunit                 =>    col%landunit                  , & ! Input:  [integer (:)]  column's landunit
   snl                       =>    cps%snl                       , & ! Input:  [integer (:)]  number of snow layers
   zi                        =>    cps%zi                        , & ! Input:  [real(r8) (:,:)]  interface level below a "z" level (m)
   dz                        =>    cps%dz                        , & ! Input:  [real(r8) (:,:)]  layer depth (m)
   z                         =>    cps%z                         , & ! Input:  [real(r8) (:,:)]  layer thickness (m)
   t_soisno                  =>    ces%t_soisno                  , & ! Input:  [real(r8) (:,:)]  soil temperature (Kelvin)
   eflx_bot                  =>    cef%eflx_bot                  , & ! Input:  [real(r8) (:)]  heat flux from beneath column (W/m**2) [+ = upward]
   begc                      =>    bounds%begc                   , & ! Input:  [integer ] beginning column index
   endc                      =>    bounds%endc                     & ! Input:  [integer ] ending column index
   )

    ! Determine heat diffusion through the layer interface and factor used in computing
    ! tridiagonal matrix and set up vector r and vectors a, b, c that define tridiagonal
    ! matrix and solve system

    do j = -nlevsno+1,nlevgrnd
       do fc = 1,num_nolakec
          c = filter_nolakec(fc)
          l = clandunit(c)
          if ((ctype(c) == icol_sunwall .or. ctype(c) == icol_shadewall &
              .or. ctype(c) == icol_roof) .and. j <= nlevurb) then
             if (j >= snl(c)+1) then
                if (j == snl(c)+1) then
                   fact(c,j) = dtime/cv(c,j)
                   fn(c,j) = tk(c,j)*(t_soisno(c,j+1)-t_soisno(c,j))/(z(c,j+1)-z(c,j))
                else if (j <= nlevurb-1) then
                   fact(c,j) = dtime/cv(c,j)
                   fn(c,j) = tk(c,j)*(t_soisno(c,j+1)-t_soisno(c,j))/(z(c,j+1)-z(c,j))
                   dzm     = (z(c,j)-z(c,j-1))
                else if (j == nlevurb) then
                   fact(c,j) = dtime/cv(c,j)
                   ! For urban sunwall, shadewall, and roof columns, there is a non-zero heat flux across
                   ! the bottom "soil" layer and the equations are derived assuming a prescribed internal
                   ! building temperature. (See Oleson urban notes of 6/18/03).
                   fn(c,j) = tk(c,j) * (t_building(l) - cnfac*t_soisno(c,j))/(zi(c,j) - z(c,j))
                end if
             end if
          else if (ctype(c) /= icol_sunwall .and. ctype(c) /= icol_shadewall &
                   .and. ctype(c) /= icol_roof) then
             if (j >= snl(c)+1) then
                if (j == snl(c)+1) then
                   fact(c,j) = dtime/cv(c,j) * dz(c,j) / (0.5_r8*(z(c,j)-zi(c,j-1)+capr*(z(c,j+1)-zi(c,j-1))))
                   fn(c,j) = tk(c,j)*(t_soisno(c,j+1)-t_soisno(c,j))/(z(c,j+1)-z(c,j))
                else if (j <= nlevgrnd-1) then
                   fact(c,j) = dtime/cv(c,j)
                   fn(c,j) = tk(c,j)*(t_soisno(c,j+1)-t_soisno(c,j))/(z(c,j+1)-z(c,j))
                   dzm     = (z(c,j)-z(c,j-1))
                else if (j == nlevgrnd) then
                   fact(c,j) = dtime/cv(c,j)
                   fn(c,j) = eflx_bot(c)
                end if
             end if
          end if
       end do
    end do

  end associate

  end subroutine ComputeHeatDiffFluxAndFactor

  !-----------------------------------------------------------------------
  subroutine SetRHSVec(bounds, num_nolakec, filter_nolakec, dtime, &
       hs_h2osfc, hs_top_snow, hs_soil, hs_top, dhsdT, sabg_lyr_col, tk, &
       tk_h2osfc, fact, fn, c_h2osfc, dz_h2osfc, rvector)

    !
    ! !DESCRIPTION:
    ! Setup the RHS-Vector for the numerical solution of temperature for snow,
    ! standing surface water and soil layers.
    !
    !           |===========|
    !           |   Snow    |
    !           !===========|
    ! rvector = |    SSW    |
    !           !===========|
    !           !   Soil    |
    !           !===========|
    !
    ! !USES:
    use clmtype
    use shr_infnan_mod  , only : nan => shr_infnan_nan, assignment(=)
    use clm_varcon      , only : cnfac, denh2o, cpliq
    use column_varcon  , only : icol_roof, icol_sunwall, icol_shadewall
    use clm_varpar     , only : nlevsno, nlevgrnd, nlevurb
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds                             ! bounds
    integer , intent(in)  :: num_nolakec                                ! number of column non-lake points in column filter
    integer , intent(in)  :: filter_nolakec(:)                          ! column filter for non-lake points
    real(r8), intent(in)  :: dtime                                      ! land model time step (sec)
    real(r8), intent(in)  :: hs_h2osfc( bounds%begc: )                  ! heat flux on standing water [W/m2]
    real(r8), intent(in)  :: hs_top_snow( bounds%begc: )                ! heat flux on top snow layer [W/m2]
    real(r8), intent(in)  :: hs_soil( bounds%begc: )                    ! heat flux on soil [W/m2]
    real(r8), intent(in)  :: hs_top( bounds%begc: )                     ! net energy flux into surface layer (col) [W/m2]
    real(r8), intent(in)  :: dhsdT( bounds%begc: )                      ! temperature derivative of "hs" [col]
    real(r8), intent(in)  :: sabg_lyr_col( bounds%begc: , -nlevsno+1: ) ! absorbed solar radiation (col,lyr) [W/m2]
    real(r8), intent(in)  :: tk( bounds%begc: , -nlevsno+1: )           ! thermal conductivity [W/(m K)]
    real(r8), intent(in)  :: tk_h2osfc( bounds%begc: )                  ! thermal conductivity of h2osfc [W/(m K)] [col]
    real(r8), intent(in)  :: fact( bounds%begc: , -nlevsno+1: )         ! used in computing tridiagonal matrix [col, lev]
    real(r8), intent(in)  :: fn( bounds%begc: , -nlevsno+1: )           ! heat diffusion through the layer interface [W/m2]
    real(r8), intent(in)  :: c_h2osfc( bounds%begc: )                   ! heat capacity of surface water [col]
    real(r8), intent(in)  :: dz_h2osfc( bounds%begc: )                  ! Thickness of standing water [m]
    real(r8), intent(out) :: rvector( bounds%begc: , -nlevsno: )        ! RHS vector used in numerical solution of temperature
    !
    ! !LOCAL VARIABLES:
    integer  :: j,c                                                     ! indices
    integer  :: fc                                                      ! lake filtered column indices
    real(r8) :: rt (bounds%begc:bounds%endc,-nlevsno+1:nlevgrnd)        ! "r" vector for tridiagonal solution
    real(r8) :: dzm                                                     ! used in computing tridiagonal matrix
    real(r8) :: dzp                                                     ! used in computing tridiagonal matrix
    real(r8) :: fn_h2osfc(bounds%begc:bounds%endc)                      ! heat diffusion through standing-water/soil interface [W/m2]
    real(r8) :: rt_snow(bounds%begc:bounds%endc,-nlevsno:-1)            ! RHS vector corresponding to snow layers
    real(r8) :: rt_ssw(bounds%begc:bounds%endc,1)                       ! RHS vector corresponding to standing surface water
    real(r8) :: rt_soil(bounds%begc:bounds%endc,1:nlevgrnd)             ! RHS vector corresponding to soil layer
    !-----------------------------------------------------------------------

    ! Enforce expected array sizes
    SHR_ASSERT_ALL((ubound(hs_h2osfc)    == (/bounds%endc/)),           errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(hs_top_snow)  == (/bounds%endc/)),           errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(hs_soil)      == (/bounds%endc/)),           errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(hs_top)       == (/bounds%endc/)),           errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(dhsdT)        == (/bounds%endc/)),           errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(sabg_lyr_col) == (/bounds%endc, 1/)),        errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(tk)           == (/bounds%endc, nlevgrnd/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(tk_h2osfc)    == (/bounds%endc/)),           errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(fact)         == (/bounds%endc, nlevgrnd/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(fn)           == (/bounds%endc, nlevgrnd/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(c_h2osfc)     == (/bounds%endc/)),           errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(dz_h2osfc)    == (/bounds%endc/)),           errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(rvector)      == (/bounds%endc, nlevgrnd/)), errMsg(__FILE__, __LINE__))

   associate(&
   frac_sno_eff              =>    cps%frac_sno_eff              , & ! Input:  [real(r8) (:)]  eff. fraction of ground covered by snow (0 to 1)
   frac_h2osfc               =>    cps%frac_h2osfc               , & ! Input:  [real(r8) (:)]  fraction of ground covered by surface water (0 to 1)
   h2osfc                    =>    cws%h2osfc                    , & ! Input:  [real(r8) (:)]  surface water (mm)
   t_h2osfc                  =>    ces%t_h2osfc                  , & ! Input:  [real(r8) (:)]  surface water temperature
   ctype                     =>    col%itype                     , & ! Input:  [integer (:)]  column type
   snl                       =>    cps%snl                       , & ! Input:  [integer (:)]  number of snow layers
   zi                        =>    cps%zi                        , & ! Input:  [real(r8) (:,:)]  interface level below a "z" level (m)
   z                         =>    cps%z                         , & ! Input:  [real(r8) (:,:)]  layer thickness (m)
   t_soisno                  =>    ces%t_soisno                  , & ! Input:  [real(r8) (:,:)]  soil temperature (Kelvin)
   begc                      =>    bounds%begc                   , & ! Input:  [integer ] beginning column index
   endc                      =>    bounds%endc                     & ! Input:  [integer ] ending column index
   )

    ! Initialize
    rvector(begc:endc, :) = nan

    ! Set entries in RHS vector for snow layers
    call SetRHSVec_Snow(bounds, num_nolakec, filter_nolakec, &
         hs_top_snow( begc:endc ),                           &
         hs_top( begc:endc ),                                &
         dhsdT( begc:endc ),                                 &
         sabg_lyr_col (begc:endc, -nlevsno+1: ),             &
         fact( begc:endc, -nlevsno+1: ),                     &
         fn( begc:endc, -nlevsno+1: ),                       &
         rt_snow( begc:endc, -nlevsno:))

    ! Set entries in RHS vector for surface water layer
    call SetRHSVec_StandingSurfaceWater(bounds, num_nolakec, filter_nolakec, &
         dtime,                                                              &
         hs_h2osfc( begc:endc ),                                             &
         dhsdT( begc:endc ),                                                 &
         tk_h2osfc( begc:endc ),                                             &
         c_h2osfc( begc:endc ),                                              &
         dz_h2osfc( begc:endc ),                                             &
         fn_h2osfc( begc:endc ),                                             &
         rt_ssw( begc:endc, 1:1))

    ! Set entries in RHS vector for soil layers
    call SetRHSVec_Soil(bounds, num_nolakec, filter_nolakec, &
         hs_top_snow( begc:endc ),                           &
         hs_soil( begc:endc ),                               &
         hs_top( begc:endc ),                                &
         dhsdT( begc:endc ),                                 &
         sabg_lyr_col (begc:endc, -nlevsno+1: ),             &
         fact( begc:endc, -nlevsno+1: ),                     &
         fn( begc:endc, -nlevsno+1: ),                       &
         fn_h2osfc( begc:endc ),                             &
         c_h2osfc( begc:endc ),                              &
         rt_soil( begc:endc, 1: ))

    ! Combine the RHS vector
    do fc = 1,num_nolakec
       c = filter_nolakec(fc)
       rvector(c, -nlevsno:-1) = rt_snow(c, -nlevsno:-1)
       rvector(c, 0         )  = rt_ssw(c, 1          )
       rvector(c, 1:nlevgrnd)  = rt_soil(c, 1:nlevgrnd )
    end do

    end associate

  end subroutine SetRHSVec

  !-----------------------------------------------------------------------
  subroutine SetRHSVec_Snow(bounds, num_nolakec, filter_nolakec, &
       hs_top_snow, hs_top, dhsdT, sabg_lyr_col, &
       fact, fn, rt)

    !
    ! !DESCRIPTION:
    ! Sets up RHS vector corresponding to snow layers.
    !
    ! !USES:
    use clmtype
    use shr_infnan_mod  , only : nan => shr_infnan_nan, assignment(=)
    use clm_varpar     , only : nlevsno, nlevgrnd
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds                             ! bounds
    integer , intent(in)  :: num_nolakec                                ! number of column non-lake points in column filter
    integer , intent(in)  :: filter_nolakec(:)                          ! column filter for non-lake points
    real(r8), intent(in)  :: hs_top_snow( bounds%begc: )                ! heat flux on top snow layer [W/m2]
    real(r8), intent(in)  :: hs_top( bounds%begc: )                     ! net energy flux into surface layer (col) [W/m2]
    real(r8), intent(in)  :: dhsdT( bounds%begc: )                      ! temperature derivative of "hs" [col]
    real(r8), intent(in)  :: sabg_lyr_col( bounds%begc: , -nlevsno+1: ) ! absorbed solar radiation (col,lyr) [W/m2]
    real(r8), intent(in)  :: fact( bounds%begc: , -nlevsno+1: )         ! used in computing tridiagonal matrix [col, lev]
    real(r8), intent(in)  :: fn (bounds%begc: , -nlevsno+1: )           ! heat diffusion through the layer interface [W/m2]
    real(r8), intent(out) :: rt(bounds%begc: , -nlevsno: )              ! rhs vector entries
    !-----------------------------------------------------------------------

    ! Enforce expected array sizes
    SHR_ASSERT_ALL((ubound(hs_top_snow)  == (/bounds%endc/)),           errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(hs_top)       == (/bounds%endc/)),           errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(dhsdT)        == (/bounds%endc/)),           errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(sabg_lyr_col) == (/bounds%endc, 1/)),        errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(fact)         == (/bounds%endc, nlevgrnd/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(fn)           == (/bounds%endc, nlevgrnd/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(rt)           == (/bounds%endc, -1/)),       errMsg(__FILE__, __LINE__))

   associate(&
   begc                      =>    bounds%begc                   , & ! Input:  [integer ] beginning column index
   endc                      =>    bounds%endc                     & ! Input:  [integer ] ending column index
   )

    ! Initialize
    rt(begc:endc, : ) = nan

    call SetRHSVec_SnowUrban(bounds, num_nolakec, filter_nolakec, &
         hs_top_snow( begc:endc ),                                &
         hs_top( begc:endc ),                                     &
         dhsdT( begc:endc ),                                      &
         sabg_lyr_col (begc:endc, -nlevsno+1: ),                  &
         fact( begc:endc, -nlevsno+1: ),                          &
         fn( begc:endc, -nlevsno+1: ),                            &
         rt( begc:endc, -nlevsno:))

    call SetRHSVec_SnowNonUrban(bounds, num_nolakec, filter_nolakec, &
         hs_top_snow( begc:endc ),                                   &
         hs_top( begc:endc ),                                        &
         dhsdT( begc:endc ),                                         &
         sabg_lyr_col (begc:endc, -nlevsno+1: ),                     &
         fact( begc:endc, -nlevsno+1: ),                             &
         fn( begc:endc, -nlevsno+1: ),                               &
         rt( begc:endc, -nlevsno:))

    end associate

  end subroutine SetRHSVec_Snow

  !-----------------------------------------------------------------------
  subroutine SetRHSVec_SnowUrban(bounds, num_nolakec, filter_nolakec, &
       hs_top_snow, hs_top, dhsdT, sabg_lyr_col, &
       fact, fn, rt)

    !
    ! !DESCRIPTION:
    ! Sets up RHS vector corresponding to snow layers for urban columns
    !
    ! !USES:
    use clmtype
    use shr_infnan_mod  , only : nan => shr_infnan_nan, assignment(=)
    use clm_varcon      , only : cnfac
    use column_varcon  , only : icol_roof, icol_sunwall, icol_shadewall
    use clm_varpar     , only : nlevsno, nlevgrnd
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds                             ! bounds
    integer , intent(in)  :: num_nolakec                                ! number of column non-lake points in column filter
    integer , intent(in)  :: filter_nolakec(:)                          ! column filter for non-lake points
    real(r8), intent(in)  :: hs_top_snow( bounds%begc: )                ! heat flux on top snow layer [W/m2]
    real(r8), intent(in)  :: hs_top( bounds%begc: )                     ! net energy flux into surface layer (col) [W/m2]
    real(r8), intent(in)  :: dhsdT( bounds%begc: )                      ! temperature derivative of "hs" [col]
    real(r8), intent(in)  :: sabg_lyr_col( bounds%begc: , -nlevsno+1: ) ! absorbed solar radiation (col,lyr) [W/m2]
    real(r8), intent(in)  :: fact( bounds%begc: , -nlevsno+1: )         ! used in computing tridiagonal matrix [col, lev]
    real(r8), intent(in)  :: fn (bounds%begc: , -nlevsno+1: )           ! heat diffusion through the layer interface [W/m2]
    real(r8), intent(inout) :: rt(bounds%begc: , -nlevsno: )            ! rhs vector entries
    !
    ! !LOCAL VARIABLES:
    integer  :: j,c,l                                                   ! indices
    integer  :: fc                                                      ! lake filtered column indices
    real(r8) :: dzm                                                     ! used in computing tridiagonal matrix
    real(r8) :: dzp                                                     ! used in computing tridiagonal matrix
    !-----------------------------------------------------------------------

    ! Enforce expected array sizes
    SHR_ASSERT_ALL((ubound(hs_top_snow)  == (/bounds%endc/)),           errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(hs_top)       == (/bounds%endc/)),           errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(dhsdT)        == (/bounds%endc/)),           errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(sabg_lyr_col) == (/bounds%endc, 1/)),        errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(fact)         == (/bounds%endc, nlevgrnd/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(fn)           == (/bounds%endc, nlevgrnd/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(rt)           == (/bounds%endc, -1/)),       errMsg(__FILE__, __LINE__))

   associate(&
   ctype                     =>    col%itype                     , & ! Input:  [integer (:)]  column type
   urbpoi                    =>    lun%urbpoi                    , & ! Input:  [logical (:)]  true => landunit is an urban point
   clandunit                 =>    col%landunit                  , & ! Input:  [integer (:)]  column's landunit
   snl                       =>    cps%snl                       , & ! Input:  [integer (:)]  number of snow layers
   z                         =>    cps%z                         , & ! Input:  [real(r8) (:,:)]  layer thickness (m)
   t_soisno                  =>    ces%t_soisno                  , & ! Input:  [real(r8) (:,:)]  soil temperature (Kelvin)
   begc                      =>    bounds%begc                   , & ! Input:  [integer ] beginning column index
   endc                      =>    bounds%endc                     & ! Input:  [integer ] ending column index
   )

    call SetRHSVec_SnowUrbanNonRoad(bounds, num_nolakec, filter_nolakec, &
         hs_top_snow( begc:endc ),                                       &
         hs_top( begc:endc ),                                            &
         dhsdT( begc:endc ),                                             &
         sabg_lyr_col (begc:endc, -nlevsno+1: ),                         &
         fact( begc:endc, -nlevsno+1: ),                                 &
         fn( begc:endc, -nlevsno+1: ),                                   &
         rt( begc:endc, -nlevsno:))

    call SetRHSVec_SnowUrbanRoad(bounds, num_nolakec, filter_nolakec, &
         hs_top_snow( begc:endc ),                                    &
         hs_top( begc:endc ),                                         &
         dhsdT( begc:endc ),                                          &
         sabg_lyr_col (begc:endc, -nlevsno+1: ),                      &
         fact( begc:endc, -nlevsno+1: ),                              &
         fn( begc:endc, -nlevsno+1: ),                                &
         rt( begc:endc, -nlevsno:))

    end associate

  end subroutine SetRHSVec_SnowUrban

  !-----------------------------------------------------------------------
  subroutine SetRHSVec_SnowUrbanNonRoad(bounds, num_nolakec, filter_nolakec, &
       hs_top_snow, hs_top, dhsdT, sabg_lyr_col, &
       fact, fn, rt)

    !
    ! !DESCRIPTION:
    ! Sets up RHS vector corresponding to snow layers for urban sunwall/shadewall/roof columns
    !
    ! !USES:
    use clmtype
    use shr_infnan_mod  , only : nan => shr_infnan_nan, assignment(=)
    use clm_varcon      , only : cnfac
    use column_varcon  , only : icol_roof, icol_sunwall, icol_shadewall
    use clm_varpar     , only : nlevsno, nlevgrnd
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds                             ! bounds
    integer , intent(in)  :: num_nolakec                                ! number of column non-lake points in column filter
    integer , intent(in)  :: filter_nolakec(:)                          ! column filter for non-lake points
    real(r8), intent(in)  :: hs_top_snow( bounds%begc: )                ! heat flux on top snow layer [W/m2]
    real(r8), intent(in)  :: hs_top( bounds%begc: )                     ! net energy flux into surface layer (col) [W/m2]
    real(r8), intent(in)  :: dhsdT( bounds%begc: )                      ! temperature derivative of "hs" [col]
    real(r8), intent(in)  :: sabg_lyr_col( bounds%begc: , -nlevsno+1: ) ! absorbed solar radiation (col,lyr) [W/m2]
    real(r8), intent(in)  :: fact( bounds%begc: , -nlevsno+1: )         ! used in computing tridiagonal matrix [col, lev]
    real(r8), intent(in)  :: fn (bounds%begc: , -nlevsno+1: )           ! heat diffusion through the layer interface [W/m2]
    real(r8), intent(inout) :: rt(bounds%begc: , -nlevsno: )            ! rhs vector entries
    !
    ! !LOCAL VARIABLES:
    integer  :: j,c,l                                                   ! indices
    integer  :: fc                                                      ! lake filtered column indices
    real(r8) :: dzm                                                     ! used in computing tridiagonal matrix
    real(r8) :: dzp                                                     ! used in computing tridiagonal matrix
    real(r8) :: rt_snow_urban(bounds%begc:bounds%endc,-nlevsno:-1)      ! rhs vector entries for urban columns
    real(r8) :: rt_snow_nonurban(bounds%begc:bounds%endc,-nlevsno:-1)   ! rhs vector entries for non-urban columns
    !-----------------------------------------------------------------------

    ! Enforce expected array sizes
    SHR_ASSERT_ALL((ubound(hs_top_snow)  == (/bounds%endc/)),           errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(hs_top)       == (/bounds%endc/)),           errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(dhsdT)        == (/bounds%endc/)),           errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(sabg_lyr_col) == (/bounds%endc, 1/)),        errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(fact)         == (/bounds%endc, nlevgrnd/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(fn)           == (/bounds%endc, nlevgrnd/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(rt)           == (/bounds%endc, -1/)),       errMsg(__FILE__, __LINE__))

   associate(&
   ctype                     =>    col%itype                     , & ! Input:  [integer (:)]  column type
   urbpoi                    =>    lun%urbpoi                    , & ! Input:  [logical (:)]  true => landunit is an urban point
   clandunit                 =>    col%landunit                  , & ! Input:  [integer (:)]  column's landunit
   snl                       =>    cps%snl                       , & ! Input:  [integer (:)]  number of snow layers
   z                         =>    cps%z                         , & ! Input:  [real(r8) (:,:)]  layer thickness (m)
   t_soisno                  =>    ces%t_soisno                  , & ! Input:  [real(r8) (:,:)]  soil temperature (Kelvin)
   begc                      =>    bounds%begc                   , & ! Input:  [integer ] beginning column index
   endc                      =>    bounds%endc                     & ! Input:  [integer ] ending column index
   )

    !
    ! urban columns ------------------------------------------------------------------
    !
    do j = -nlevsno+1,0
       do fc = 1,num_nolakec
          c = filter_nolakec(fc)
          l = clandunit(c)
          if (urbpoi(l)) then
            if ((ctype(c) == icol_sunwall .or. ctype(c) == icol_shadewall &
                .or. ctype(c) == icol_roof)) then
               if (j >= snl(c)+1) then
                  if (j == snl(c)+1) then
                     dzp     = z(c,j+1)-z(c,j)
                     ! changed hs to hs_top
                     rt(c,j-1) = t_soisno(c,j) +  fact(c,j)*( hs_top(c) - dhsdT(c)*t_soisno(c,j) + cnfac*fn(c,j) )
                  else
                     dzm     = (z(c,j)-z(c,j-1))
                     dzp     = (z(c,j+1)-z(c,j))
                     rt(c,j-1) = t_soisno(c,j) + cnfac*fact(c,j)*( fn(c,j) - fn(c,j-1) )
                     rt(c,j-1) = rt(c,j-1) + (fact(c,j)*sabg_lyr_col(c,j))
                  end if
               end if
            end if
          end if
       enddo
    end do

    end associate

  end subroutine SetRHSVec_SnowUrbanNonRoad

  !-----------------------------------------------------------------------
  subroutine SetRHSVec_SnowUrbanRoad(bounds, num_nolakec, filter_nolakec, &
       hs_top_snow, hs_top, dhsdT, sabg_lyr_col, &
       fact, fn, rt)

    !
    ! !DESCRIPTION:
    ! Sets up RHS vector corresponding to snow layers for urban road
    ! (impervious + pervious) columns
    !
    ! !USES:
    use clmtype
    use shr_infnan_mod  , only : nan => shr_infnan_nan, assignment(=)
    use clm_varcon      , only : cnfac
    use column_varcon  , only : icol_road_perv, icol_road_imperv
    use clm_varpar     , only : nlevsno, nlevgrnd
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds                             ! bounds
    integer , intent(in)  :: num_nolakec                                ! number of column non-lake points in column filter
    integer , intent(in)  :: filter_nolakec(:)                          ! column filter for non-lake points
    real(r8), intent(in)  :: hs_top_snow( bounds%begc: )                ! heat flux on top snow layer [W/m2]
    real(r8), intent(in)  :: hs_top( bounds%begc: )                     ! net energy flux into surface layer (col) [W/m2]
    real(r8), intent(in)  :: dhsdT( bounds%begc: )                      ! temperature derivative of "hs" [col]
    real(r8), intent(in)  :: sabg_lyr_col( bounds%begc: , -nlevsno+1: ) ! absorbed solar radiation (col,lyr) [W/m2]
    real(r8), intent(in)  :: fact( bounds%begc: , -nlevsno+1: )         ! used in computing tridiagonal matrix [col, lev]
    real(r8), intent(in)  :: fn (bounds%begc: , -nlevsno+1: )           ! heat diffusion through the layer interface [W/m2]
    real(r8), intent(inout) :: rt(bounds%begc: , -nlevsno: )            ! rhs vector entries
    !
    ! !LOCAL VARIABLES:
    integer  :: j,c,l                                                   ! indices
    integer  :: fc                                                      ! lake filtered column indices
    real(r8) :: dzm                                                     ! used in computing tridiagonal matrix
    real(r8) :: dzp                                                     ! used in computing tridiagonal matrix
    real(r8) :: rt_snow_urban(bounds%begc:bounds%endc,-nlevsno:-1)      !
    real(r8) :: rt_snow_nonurban(bounds%begc:bounds%endc,-nlevsno:-1)   !
    !-----------------------------------------------------------------------

    ! Enforce expected array sizes
    SHR_ASSERT_ALL((ubound(hs_top_snow)  == (/bounds%endc/)),           errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(hs_top)       == (/bounds%endc/)),           errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(dhsdT)        == (/bounds%endc/)),           errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(sabg_lyr_col) == (/bounds%endc, 1/)),        errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(fact)         == (/bounds%endc, nlevgrnd/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(fn)           == (/bounds%endc, nlevgrnd/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(rt)           == (/bounds%endc, -1/)),       errMsg(__FILE__, __LINE__))

   associate(&
   ctype                     =>    col%itype                     , & ! Input:  [integer (:)]  column type
   urbpoi                    =>    lun%urbpoi                    , & ! Input:  [logical (:)]  true => landunit is an urban point
   clandunit                 =>    col%landunit                  , & ! Input:  [integer (:)]  column's landunit
   snl                       =>    cps%snl                       , & ! Input:  [integer (:)]  number of snow layers
   z                         =>    cps%z                         , & ! Input:  [real(r8) (:,:)]  layer thickness (m)
   t_soisno                  =>    ces%t_soisno                  , & ! Input:  [real(r8) (:,:)]  soil temperature (Kelvin)
   begc                      =>    bounds%begc                   , & ! Input:  [integer ] beginning column index
   endc                      =>    bounds%endc                     & ! Input:  [integer ] ending column index
   )

    !
    ! urban road columns -------------------------------------------------------------
    !
    do j = -nlevsno+1,0
       do fc = 1,num_nolakec
          c = filter_nolakec(fc)
          l = clandunit(c)
          if (urbpoi(l)) then
            if (ctype(c) == icol_road_imperv .or. ctype(c) == icol_road_perv) then
               if (j >= snl(c)+1) then
                  if (j == snl(c)+1) then
                     dzp     = z(c,j+1)-z(c,j)
                     rt(c,j-1) = t_soisno(c,j) +  fact(c,j)*( hs_top_snow(c) &
                          - dhsdT(c)*t_soisno(c,j) + cnfac*fn(c,j) )
                  else
                     dzm     = (z(c,j)-z(c,j-1))
                     dzp     = (z(c,j+1)-z(c,j))

                     rt(c,j-1) = t_soisno(c,j) + cnfac*fact(c,j)*( fn(c,j) - fn(c,j-1) )
                     rt(c,j-1) = rt(c,j-1) + fact(c,j)*sabg_lyr_col(c,j)

                  end if
               end if
            end if
          end if
       enddo
    end do

    end associate

  end subroutine SetRHSVec_SnowUrbanRoad

  !-----------------------------------------------------------------------
  subroutine SetRHSVec_SnowNonUrban(bounds, num_nolakec, filter_nolakec, &
       hs_top_snow, hs_top, dhsdT, sabg_lyr_col, &
       fact, fn, rt)

    !
    ! !DESCRIPTION:
    ! Sets up RHS vector corresponding to snow layers for non-urban columns
    !
    ! !USES:
    use clmtype
    use shr_infnan_mod  , only : nan => shr_infnan_nan, assignment(=)
    use clm_varcon      , only : cnfac
    use column_varcon  , only : icol_roof, icol_sunwall, icol_shadewall
    use clm_varpar     , only : nlevsno, nlevgrnd
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds                             ! bounds
    integer , intent(in)  :: num_nolakec                                ! number of column non-lake points in column filter
    integer , intent(in)  :: filter_nolakec(:)                          ! column filter for non-lake points
    real(r8), intent(in)  :: hs_top_snow( bounds%begc: )                ! heat flux on top snow layer [W/m2]
    real(r8), intent(in)  :: hs_top( bounds%begc: )                     ! net energy flux into surface layer (col) [W/m2]
    real(r8), intent(in)  :: dhsdT( bounds%begc: )                      ! temperature derivative of "hs" [col]
    real(r8), intent(in)  :: sabg_lyr_col( bounds%begc: , -nlevsno+1: ) ! absorbed solar radiation (col,lyr) [W/m2]
    real(r8), intent(in)  :: fact( bounds%begc: , -nlevsno+1: )         ! used in computing tridiagonal matrix [col, lev]
    real(r8), intent(in)  :: fn (bounds%begc: , -nlevsno+1: )           ! heat diffusion through the layer interface [W/m2]
    real(r8), intent(inout) :: rt(bounds%begc: , -nlevsno: )            ! rhs vector entries
    !
    ! !LOCAL VARIABLES:
    integer  :: j,c,l                                                   ! indices
    integer  :: fc                                                      ! lake filtered column indices
    real(r8) :: dzm                                                     ! used in computing tridiagonal matrix
    real(r8) :: dzp                                                     ! used in computing tridiagonal matrix
    !-----------------------------------------------------------------------

    ! Enforce expected array sizes
    SHR_ASSERT_ALL((ubound(hs_top_snow)  == (/bounds%endc/)),           errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(hs_top)       == (/bounds%endc/)),           errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(dhsdT)        == (/bounds%endc/)),           errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(sabg_lyr_col) == (/bounds%endc, 1/)),        errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(fact)         == (/bounds%endc, nlevgrnd/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(fn)           == (/bounds%endc, nlevgrnd/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(rt)           == (/bounds%endc, -1/)),       errMsg(__FILE__, __LINE__))

   associate(&
   ctype                     =>    col%itype                     , & ! Input:  [integer (:)]  column type
   urbpoi                    =>    lun%urbpoi                    , & ! Input:  [logical (:)]  true => landunit is an urban point
   clandunit                 =>    col%landunit                  , & ! Input:  [integer (:)]  column's landunit
   snl                       =>    cps%snl                       , & ! Input:  [integer (:)]  number of snow layers
   z                         =>    cps%z                         , & ! Input:  [real(r8) (:,:)]  layer thickness (m)
   t_soisno                  =>    ces%t_soisno                  , & ! Input:  [real(r8) (:,:)]  soil temperature (Kelvin)
   begc                      =>    bounds%begc                   , & ! Input:  [integer ] beginning column index
   endc                      =>    bounds%endc                     & ! Input:  [integer ] ending column index
   )

    !
    ! non-urban columns --------------------------------------------------------------
    !
    do j = -nlevsno+1,0
       do fc = 1,num_nolakec
          c = filter_nolakec(fc)
          l = clandunit(c)
          if (.not. urbpoi(l)) then
             if (j >= snl(c)+1) then
                if (j == snl(c)+1) then
                   dzp     = z(c,j+1)-z(c,j)
                   rt(c,j-1) = t_soisno(c,j) +  fact(c,j)*( hs_top_snow(c) &
                        - dhsdT(c)*t_soisno(c,j) + cnfac*fn(c,j) )

                else
                   dzm     = (z(c,j)-z(c,j-1))
                   dzp     = (z(c,j+1)-z(c,j))

                   rt(c,j-1) = t_soisno(c,j) + cnfac*fact(c,j)*( fn(c,j) - fn(c,j-1) )
                   rt(c,j-1) = rt(c,j-1) + fact(c,j)*sabg_lyr_col(c,j)

                end if
             end if
          end if
       enddo
    end do

    end associate

  end subroutine SetRHSVec_SnowNonUrban

  !-----------------------------------------------------------------------
  subroutine SetRHSVec_StandingSurfaceWater(bounds, num_nolakec, filter_nolakec, dtime, &
       hs_h2osfc, dhsdT, tk_h2osfc, c_h2osfc, dz_h2osfc, fn_h2osfc, rt)

    !
    ! !DESCRIPTION:
    ! Sets up RHS vector corresponding to standing surface water
    !
    ! !USES:
    use clmtype
    use shr_infnan_mod  , only : nan => shr_infnan_nan, assignment(=)
    use clm_varcon      , only : cnfac
    use column_varcon  , only : icol_roof, icol_sunwall, icol_shadewall
    use clm_varpar     , only : nlevsno, nlevgrnd
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds                     ! bounds
    integer , intent(in)  :: num_nolakec                        ! number of column non-lake points in column filter
    integer , intent(in)  :: filter_nolakec(:)                  ! column filter for non-lake points
    real(r8), intent(in)  :: dtime                              ! land model time step (sec)
    real(r8), intent(in)  :: hs_h2osfc(bounds%begc: )           !
    real(r8), intent(in)  :: dhsdT(bounds%begc: )               ! temperature derivative of "hs" [col]
    real(r8), intent(in)  :: tk_h2osfc(bounds%begc: )           !
    real(r8), intent(in)  :: c_h2osfc( bounds%begc: )           ! heat capacity of surface water [col]
    real(r8), intent(in)  :: dz_h2osfc(bounds%begc: )           ! Thickness of standing water [m]
    real(r8), intent(out) :: fn_h2osfc (bounds%begc: )          ! heat diffusion through standing-water/soil interface [W/m2]
    real(r8), intent(out) :: rt(bounds%begc:bounds%endc, 1:1 )  ! rhs vector entries
    !
    ! !LOCAL VARIABLES:
    integer  :: j,c                                             ! indices
    integer  :: fc                                              ! lake filtered column indices
    real(r8) :: dzm                                             ! used in computing tridiagonal matrix
    !-----------------------------------------------------------------------

    ! Enforce expected array sizes
    SHR_ASSERT_ALL((ubound(hs_h2osfc)    == (/bounds%endc/)),   errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(dhsdT)        == (/bounds%endc/)),   errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(tk_h2osfc)    == (/bounds%endc/)),   errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(c_h2osfc)     == (/bounds%endc/)),   errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(dz_h2osfc)    == (/bounds%endc/)),   errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(fn_h2osfc)    == (/bounds%endc/)),   errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(rt)           == (/bounds%endc,1/)), errMsg(__FILE__, __LINE__))

   associate(&
   t_h2osfc                  =>    ces%t_h2osfc                  , & ! Input:  [real(r8) (:)]  surface water temperature
   ctype                     =>    col%itype                     , & ! Input:  [integer (:)]  column type
   z                         =>    cps%z                         , & ! Input:  [real(r8) (:,:)]  layer thickness (m)
   t_soisno                  =>    ces%t_soisno                  , & ! Input:  [real(r8) (:,:)]  soil temperature (Kelvin)
   begc                      =>    bounds%begc                   , & ! Input:  [integer ] beginning column index
   endc                      =>    bounds%endc                     & ! Input:  [integer ] ending column index
   )

    ! Initialize
    rt(begc:endc, : ) = nan

    !
    ! surface water ------------------------------------------------------------------
    !
    do fc = 1,num_nolakec
       c = filter_nolakec(fc)

       ! surface water layer has two coefficients
       dzm=(0.5*dz_h2osfc(c)+z(c,1))

       fn_h2osfc(c)=tk_h2osfc(c)*(t_soisno(c,1)-t_h2osfc(c))/dzm
       rt(c,1)= t_h2osfc(c) +  (dtime/c_h2osfc(c)) &
            *( hs_h2osfc(c) - dhsdT(c)*t_h2osfc(c) + cnfac*fn_h2osfc(c) )!rhs for h2osfc

    enddo

    end associate

  end subroutine SetRHSVec_StandingSurfaceWater

  !-----------------------------------------------------------------------
  subroutine SetRHSVec_Soil(bounds, num_nolakec, filter_nolakec, &
       hs_top_snow, hs_soil, hs_top, dhsdT, sabg_lyr_col, fact, fn, fn_h2osfc, c_h2osfc, rt)

    !
    ! !DESCRIPTION:
    ! Sets up RHS vector corresponding to soil layers
    !
    ! !USES:
    use clmtype
    use shr_infnan_mod  , only : nan => shr_infnan_nan, assignment(=)
    use clm_varcon      , only : cnfac
    use column_varcon  , only : icol_roof, icol_sunwall, icol_shadewall
    use clm_varpar     , only : nlevsno, nlevgrnd, nlevurb
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds                                     ! bounds
    integer , intent(in)  :: num_nolakec                                        ! number of column non-lake points in column filter
    integer , intent(in)  :: filter_nolakec(:)                                  ! column filter for non-lake points
    real(r8), intent(in)  :: hs_top_snow(bounds%begc: )                         ! heat flux on top snow layer [W/m2]
    real(r8), intent(in)  :: hs_soil(bounds%begc: )                             ! heat flux on soil [W/m2]
    real(r8), intent(in)  :: hs_top(bounds%begc: )                              ! net energy flux into surface layer (col) [W/m2]
    real(r8), intent(in)  :: dhsdT(bounds%begc: )                               ! temperature derivative of "hs" [col]
    real(r8), intent(in)  :: sabg_lyr_col(bounds%begc:, -nlevsno+1: )           ! absorbed solar radiation (col,lyr) [W/m2]
    real(r8), intent(in)  :: fact( bounds%begc: , -nlevsno+1: )                 ! used in computing tridiagonal matrix [col, lev]
    real(r8), intent(in)  :: fn (bounds%begc: ,-nlevsno+1: )                    ! heat diffusion through the layer interface [W/m2]
    real(r8), intent(in)  :: fn_h2osfc (bounds%begc: )                          ! heat diffusion through standing-water/soil interface [W/m2]
    real(r8), intent(in)  :: c_h2osfc( bounds%begc: )                           ! heat capacity of surface water [col]
    real(r8), intent(out) :: rt(bounds%begc: ,1: )                              ! rhs vector entries
    !-----------------------------------------------------------------------

    ! Enforce expected array sizes
    SHR_ASSERT_ALL((ubound(hs_soil)      == (/bounds%endc/)),           errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(dhsdT)        == (/bounds%endc/)),           errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(sabg_lyr_col) == (/bounds%endc, 1/)),        errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(fact)         == (/bounds%endc, nlevgrnd/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(fn)           == (/bounds%endc, nlevgrnd/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(fn_h2osfc)    == (/bounds%endc/)),           errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(c_h2osfc)     == (/bounds%endc/)),           errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(rt)           == (/bounds%endc, nlevgrnd/)), errMsg(__FILE__, __LINE__))

   associate(&
   begc                      =>    bounds%begc                   , & ! Input:  [integer ] beginning column index
   endc                      =>    bounds%endc                     & ! Input:  [integer ] ending column index
   )

    ! Initialize
    rt(begc:endc, : ) = nan

    call SetRHSVec_SoilUrban(bounds, num_nolakec, filter_nolakec, &
         hs_top_snow( begc:endc ),                                &
         hs_soil( begc:endc ),                                    &
         hs_top( begc:endc ),                                     &
         dhsdT( begc:endc ),                                      &
         sabg_lyr_col (begc:endc, -nlevsno+1: ),                  &
         fact( begc:endc, -nlevsno+1: ),                          &
         fn( begc:endc, -nlevsno+1: ),                            &
         fn_h2osfc( begc:endc ),                                  &
         c_h2osfc( begc:endc ),                                   &
         rt( begc:endc, 1: ))

    call SetRHSVec_SoilNonUrban(bounds, num_nolakec, filter_nolakec, &
         hs_top_snow( begc:endc ),                                   &
         hs_soil( begc:endc ),                                       &
         hs_top( begc:endc ),                                        &
         dhsdT( begc:endc ),                                         &
         sabg_lyr_col (begc:endc, -nlevsno+1: ),                     &
         fact( begc:endc, -nlevsno+1: ),                             &
         fn( begc:endc, -nlevsno+1: ),                               &
         fn_h2osfc( begc:endc ),                                     &
         c_h2osfc( begc:endc ),                                      &
         rt( begc:endc, 1: ))

    call SetRHSVec_Soil_StandingSurfaceWater(bounds, num_nolakec, filter_nolakec, &
         hs_top_snow( begc:endc ),                                                &
         hs_soil( begc:endc ),                                                    &
         hs_top( begc:endc ),                                                     &
         dhsdT( begc:endc ),                                                      &
         sabg_lyr_col (begc:endc, -nlevsno+1: ),                                  &
         fact( begc:endc, -nlevsno+1: ),                                          &
         fn( begc:endc, -nlevsno+1: ),                                            &
         fn_h2osfc( begc:endc ),                                                  &
         c_h2osfc( begc:endc ),                                                   &
         rt( begc:endc, 1: ))

    end associate

  end subroutine SetRHSVec_Soil

  !-----------------------------------------------------------------------
  subroutine SetRHSVec_SoilUrban(bounds, num_nolakec, filter_nolakec, &
       hs_top_snow, hs_soil, hs_top, dhsdT, sabg_lyr_col, fact, fn, fn_h2osfc, c_h2osfc, rt)

    !
    ! !DESCRIPTION:
    ! Sets up RHS vector corresponding to soil layers for urban columns
    !
    ! !USES:
    use clmtype
    use shr_infnan_mod  , only : nan => shr_infnan_nan, assignment(=)
    use clm_varcon      , only : cnfac
    use column_varcon  , only : icol_roof, icol_sunwall, icol_shadewall
    use clm_varpar     , only : nlevsno, nlevgrnd, nlevurb
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds                                     ! bounds
    integer , intent(in)  :: num_nolakec                                        ! number of column non-lake points in column filter
    integer , intent(in)  :: filter_nolakec(:)                                  ! column filter for non-lake points
    real(r8), intent(in)  :: hs_top_snow(bounds%begc: )                         ! heat flux on top snow layer [W/m2]
    real(r8), intent(in)  :: hs_soil(bounds%begc: )                             ! heat flux on soil [W/m2]
    real(r8), intent(in)  :: hs_top(bounds%begc: )                              ! net energy flux into surface layer (col) [W/m2]
    real(r8), intent(in)  :: dhsdT(bounds%begc: )                               ! temperature derivative of "hs" [col]
    real(r8), intent(in)  :: sabg_lyr_col(bounds%begc:, -nlevsno+1: )           ! absorbed solar radiation (col,lyr) [W/m2]
    real(r8), intent(in)  :: fact( bounds%begc: , -nlevsno+1: )                 ! used in computing tridiagonal matrix [col, lev]
    real(r8), intent(in)  :: fn (bounds%begc: ,-nlevsno+1: )                    ! heat diffusion through the layer interface [W/m2]
    real(r8), intent(in)  :: fn_h2osfc (bounds%begc: )                          ! heat diffusion through standing-water/soil interface [W/m2]
    real(r8), intent(in)  :: c_h2osfc( bounds%begc: )                           ! heat capacity of surface water [col]
    real(r8), intent(inout) :: rt(bounds%begc: ,1: )                            ! rhs vector entries
    !
    ! !LOCAL VARIABLES:
    integer  :: j,c,l                                                           ! indices
    integer  :: fc                                                              ! lake filtered column indices
    !-----------------------------------------------------------------------

    ! Enforce expected array sizes
    SHR_ASSERT_ALL((ubound(hs_soil)      == (/bounds%endc/)),           errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(dhsdT)        == (/bounds%endc/)),           errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(sabg_lyr_col) == (/bounds%endc, 1/)),        errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(fact)         == (/bounds%endc, nlevgrnd/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(fn)           == (/bounds%endc, nlevgrnd/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(fn_h2osfc)    == (/bounds%endc/)),           errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(c_h2osfc)     == (/bounds%endc/)),           errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(rt)           == (/bounds%endc, nlevgrnd/)), errMsg(__FILE__, __LINE__))

   associate(&
   snl                       =>    cps%snl                       , & ! Input:  [integer (:)]  number of snow layers
   urbpoi                    =>    lun%urbpoi                    , & ! Input:  [logical (:)]  true => landunit is an urban point
   clandunit                 =>    col%landunit                  , & ! Input:  [integer (:)]  column's landunit
   frac_sno_eff              =>    cps%frac_sno_eff              , & ! Input:  [real(r8) (:)]  eff. fraction of ground covered by snow (0 to 1)
   frac_h2osfc               =>    cps%frac_h2osfc               , & ! Input:  [real(r8) (:)]  fraction of ground covered by surface water (0 to 1)
   ctype                     =>    col%itype                     , & ! Input:  [integer (:)]  column type
   zi                        =>    cps%zi                        , & ! Input:  [real(r8) (:,:)]  interface level below a "z" level (m)
   z                         =>    cps%z                         , & ! Input:  [real(r8) (:,:)]  layer thickness (m)
   t_soisno                  =>    ces%t_soisno                  , & ! Input:  [real(r8) (:,:)]  soil temperature (Kelvin)
   begc                      =>    bounds%begc                   , & ! Input:  [integer ] beginning column index
   endc                      =>    bounds%endc                     & ! Input:  [integer ] ending column index
   )

    call SetRHSVec_SoilUrbanNonRoad(bounds, num_nolakec, filter_nolakec, &
         hs_top_snow( begc:endc ),                                       &
         hs_soil( begc:endc ),                                           &
         hs_top( begc:endc ),                                            &
         dhsdT( begc:endc ),                                             &
         sabg_lyr_col (begc:endc, -nlevsno+1: ),                         &
         fact( begc:endc, -nlevsno+1: ),                                 &
         fn( begc:endc, -nlevsno+1: ),                                   &
         fn_h2osfc( begc:endc ),                                         &
         c_h2osfc( begc:endc ),                                          &
         rt( begc:endc, 1: ))

    call SetRHSVec_SoilUrbanRoad(bounds, num_nolakec, filter_nolakec, &
         hs_top_snow( begc:endc ),                                    &
         hs_soil( begc:endc ),                                        &
         hs_top( begc:endc ),                                         &
         dhsdT( begc:endc ),                                          &
         sabg_lyr_col (begc:endc, -nlevsno+1: ),                      &
         fact( begc:endc, -nlevsno+1: ),                              &
         fn( begc:endc, -nlevsno+1: ),                                &
         fn_h2osfc( begc:endc ),                                      &
         c_h2osfc( begc:endc ),                                       &
         rt( begc:endc, 1: ))

    end associate

  end subroutine SetRHSVec_SoilUrban

  !-----------------------------------------------------------------------
  subroutine SetRHSVec_SoilUrbanNonRoad(bounds, num_nolakec, filter_nolakec, &
       hs_top_snow, hs_soil, hs_top, dhsdT, sabg_lyr_col, fact, fn, fn_h2osfc, c_h2osfc, rt)

    !
    ! !DESCRIPTION:
    ! Sets up RHS vector corresponding to soil layers for urban sunwall/shadewall/roof columns
    !
    ! !USES:
    use clmtype
    use shr_infnan_mod  , only : nan => shr_infnan_nan, assignment(=)
    use clm_varcon      , only : cnfac
    use column_varcon  , only : icol_roof, icol_sunwall, icol_shadewall
    use clm_varpar     , only : nlevsno, nlevgrnd, nlevurb
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds                                     ! bounds
    integer , intent(in)  :: num_nolakec                                        ! number of column non-lake points in column filter
    integer , intent(in)  :: filter_nolakec(:)                                  ! column filter for non-lake points
    real(r8), intent(in)  :: hs_top_snow(bounds%begc: )                         ! heat flux on top snow layer [W/m2]
    real(r8), intent(in)  :: hs_soil(bounds%begc: )                             ! heat flux on soil [W/m2]
    real(r8), intent(in)  :: hs_top(bounds%begc: )                              ! net energy flux into surface layer (col) [W/m2]
    real(r8), intent(in)  :: dhsdT(bounds%begc: )                               ! temperature derivative of "hs" [col]
    real(r8), intent(in)  :: sabg_lyr_col(bounds%begc:, -nlevsno+1: )           ! absorbed solar radiation (col,lyr) [W/m2]
    real(r8), intent(in)  :: fact( bounds%begc: , -nlevsno+1: )                 ! used in computing tridiagonal matrix [col, lev]
    real(r8), intent(in)  :: fn (bounds%begc: ,-nlevsno+1: )                    ! heat diffusion through the layer interface [W/m2]
    real(r8), intent(in)  :: fn_h2osfc (bounds%begc: )                          ! heat diffusion through standing-water/soil interface [W/m2]
    real(r8), intent(in)  :: c_h2osfc( bounds%begc: )                           ! heat capacity of surface water [col]
    real(r8), intent(inout) :: rt(bounds%begc: ,1: )                            ! rhs vector entries
    !
    ! !LOCAL VARIABLES:
    integer  :: j,c,l                                                           ! indices
    integer  :: fc                                                              ! lake filtered column indices
    !-----------------------------------------------------------------------

    ! Enforce expected array sizes
    SHR_ASSERT_ALL((ubound(hs_soil)      == (/bounds%endc/)),           errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(dhsdT)        == (/bounds%endc/)),           errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(sabg_lyr_col) == (/bounds%endc, 1/)),        errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(fact)         == (/bounds%endc, nlevgrnd/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(fn)           == (/bounds%endc, nlevgrnd/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(fn_h2osfc)    == (/bounds%endc/)),           errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(c_h2osfc)     == (/bounds%endc/)),           errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(rt)           == (/bounds%endc, nlevgrnd/)), errMsg(__FILE__, __LINE__))

   associate(&
   snl                       =>    cps%snl                       , & ! Input:  [integer (:)]  number of snow layers
   urbpoi                    =>    lun%urbpoi                    , & ! Input:  [logical (:)]  true => landunit is an urban point
   clandunit                 =>    col%landunit                  , & ! Input:  [integer (:)]  column's landunit
   frac_sno_eff              =>    cps%frac_sno_eff              , & ! Input:  [real(r8) (:)]  eff. fraction of ground covered by snow (0 to 1)
   frac_h2osfc               =>    cps%frac_h2osfc               , & ! Input:  [real(r8) (:)]  fraction of ground covered by surface water (0 to 1)
   ctype                     =>    col%itype                     , & ! Input:  [integer (:)]  column type
   zi                        =>    cps%zi                        , & ! Input:  [real(r8) (:,:)]  interface level below a "z" level (m)
   z                         =>    cps%z                         , & ! Input:  [real(r8) (:,:)]  layer thickness (m)
   t_soisno                  =>    ces%t_soisno                  , & ! Input:  [real(r8) (:,:)]  soil temperature (Kelvin)
   begc                      =>    bounds%begc                   , & ! Input:  [integer ] beginning column index
   endc                      =>    bounds%endc                     & ! Input:  [integer ] ending column index
   )

    !
    ! urban columns ------------------------------------------------------------------
    !
    do j = 1,nlevurb
       do fc = 1,num_nolakec
          c = filter_nolakec(fc)
          l = clandunit(c)
          if (urbpoi(l)) then
            if ((ctype(c) == icol_sunwall .or. ctype(c) == icol_shadewall &
                .or. ctype(c) == icol_roof)) then
               if (j >= snl(c)+1) then
                  if (j == snl(c)+1) then
                     ! changed hs to hs_top
                     rt(c,j) = t_soisno(c,j) +  fact(c,j)*( hs_top(c) - dhsdT(c)*t_soisno(c,j) + cnfac*fn(c,j) )
                  else if (j <= nlevurb-1) then
                     ! if this is a snow layer or the top soil layer,
                     ! add absorbed solar flux to factor 'rt'
                     if (j == 1) then
                        rt(c,j) = t_soisno(c,j) + cnfac*fact(c,j)*( fn(c,j) - fn(c,j-1) )
                        rt(c,j) = rt(c,j) + (fact(c,j)*sabg_lyr_col(c,j))
                     else
                        rt(c,j) = t_soisno(c,j) + cnfac*fact(c,j)*( fn(c,j) - fn(c,j-1) )
                     endif

                  else if (j == nlevurb) then
                     ! For urban sunwall, shadewall, and roof columns, there is a non-zero heat flux across
                     ! the bottom "soil" layer and the equations are derived assuming a prescribed internal
                     ! building temperature. (See Oleson urban notes of 6/18/03).
                     rt(c,j) = t_soisno(c,j) + fact(c,j)*( fn(c,j) - cnfac*fn(c,j-1) )
                  end if
               end if
            end if
          end if
       enddo
    end do

    end associate

  end subroutine SetRHSVec_SoilUrbanNonRoad

  !-----------------------------------------------------------------------
  subroutine SetRHSVec_SoilUrbanRoad(bounds, num_nolakec, filter_nolakec, &
       hs_top_snow, hs_soil, hs_top, dhsdT, sabg_lyr_col, fact, fn, fn_h2osfc, c_h2osfc, rt)

    !
    ! !DESCRIPTION:
    ! Sets up RHS vector corresponding to soil layers for urban road
    ! (impervious + pervious) columns
    !
    ! !USES:
    use clmtype
    use shr_infnan_mod  , only : nan => shr_infnan_nan, assignment(=)
    use clm_varcon      , only : cnfac
    use column_varcon  , only : icol_road_perv, icol_road_imperv
    use clm_varpar     , only : nlevsno, nlevgrnd, nlevurb
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds                                     ! bounds
    integer , intent(in)  :: num_nolakec                                        ! number of column non-lake points in column filter
    integer , intent(in)  :: filter_nolakec(:)                                  ! column filter for non-lake points
    real(r8), intent(in)  :: hs_top_snow(bounds%begc: )                         ! heat flux on top snow layer [W/m2]
    real(r8), intent(in)  :: hs_soil(bounds%begc: )                             ! heat flux on soil [W/m2]
    real(r8), intent(in)  :: hs_top(bounds%begc: )                              ! net energy flux into surface layer (col) [W/m2]
    real(r8), intent(in)  :: dhsdT(bounds%begc: )                               ! temperature derivative of "hs" [col]
    real(r8), intent(in)  :: sabg_lyr_col(bounds%begc:, -nlevsno+1: )           ! absorbed solar radiation (col,lyr) [W/m2]
    real(r8), intent(in)  :: fact( bounds%begc: , -nlevsno+1: )                 ! used in computing tridiagonal matrix [col, lev]
    real(r8), intent(in)  :: fn (bounds%begc: ,-nlevsno+1: )                    ! heat diffusion through the layer interface [W/m2]
    real(r8), intent(in)  :: fn_h2osfc (bounds%begc: )                          ! heat diffusion through standing-water/soil interface [W/m2]
    real(r8), intent(in)  :: c_h2osfc( bounds%begc: )                           ! heat capacity of surface water [col]
    real(r8), intent(inout) :: rt(bounds%begc: ,1: )                            ! rhs vector entries
    !
    ! !LOCAL VARIABLES:
    integer  :: j,c,l                                                           ! indices
    integer  :: fc                                                              ! lake filtered column indices
    !-----------------------------------------------------------------------

    ! Enforce expected array sizes
    SHR_ASSERT_ALL((ubound(hs_soil)      == (/bounds%endc/)),           errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(dhsdT)        == (/bounds%endc/)),           errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(sabg_lyr_col) == (/bounds%endc, 1/)),        errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(fact)         == (/bounds%endc, nlevgrnd/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(fn)           == (/bounds%endc, nlevgrnd/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(fn_h2osfc)    == (/bounds%endc/)),           errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(c_h2osfc)     == (/bounds%endc/)),           errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(rt)           == (/bounds%endc, nlevgrnd/)), errMsg(__FILE__, __LINE__))

   associate(&
   snl                       =>    cps%snl                       , & ! Input:  [integer (:)]  number of snow layers
   urbpoi                    =>    lun%urbpoi                    , & ! Input:  [logical (:)]  true => landunit is an urban point
   clandunit                 =>    col%landunit                  , & ! Input:  [integer (:)]  column's landunit
   frac_sno_eff              =>    cps%frac_sno_eff              , & ! Input:  [real(r8) (:)]  eff. fraction of ground covered by snow (0 to 1)
   frac_h2osfc               =>    cps%frac_h2osfc               , & ! Input:  [real(r8) (:)]  fraction of ground covered by surface water (0 to 1)
   ctype                     =>    col%itype                     , & ! Input:  [integer (:)]  column type
   zi                        =>    cps%zi                        , & ! Input:  [real(r8) (:,:)]  interface level below a "z" level (m)
   z                         =>    cps%z                         , & ! Input:  [real(r8) (:,:)]  layer thickness (m)
   t_soisno                  =>    ces%t_soisno                  , & ! Input:  [real(r8) (:,:)]  soil temperature (Kelvin)
   begc                      =>    bounds%begc                   , & ! Input:  [integer ] beginning column index
   endc                      =>    bounds%endc                     & ! Input:  [integer ] ending column index
   )

    !
    ! urban road columns -------------------------------------------------------------
    !
    do j = 1,nlevgrnd
       do fc = 1,num_nolakec
          c = filter_nolakec(fc)
          l = clandunit(c)
          if (urbpoi(l)) then
            if (ctype(c) == icol_road_imperv .or. ctype(c) == icol_road_perv) then
               if (j == snl(c)+1) then
                     rt(c,j) = t_soisno(c,j) +  fact(c,j)*( hs_top_snow(c) &
                          - dhsdT(c)*t_soisno(c,j) + cnfac*fn(c,j) )
               else if (j == 1) then
                  ! this is the snow/soil interface layer
                  rt(c,j) = t_soisno(c,j) + fact(c,j) &
                       *((1._r8-frac_sno_eff(c))*(hs_soil(c) - dhsdT(c)*t_soisno(c,j)) &
                       + cnfac*(fn(c,j) - frac_sno_eff(c) * fn(c,j-1)))

                  rt(c,j) = rt(c,j) +  frac_sno_eff(c)*fact(c,j)*sabg_lyr_col(c,j)

               else if (j <= nlevgrnd-1) then
                  rt(c,j) = t_soisno(c,j) + cnfac*fact(c,j)*( fn(c,j) - fn(c,j-1) )

               else if (j == nlevgrnd) then
                  rt(c,j) = t_soisno(c,j) - cnfac*fact(c,j)*fn(c,j-1) + fact(c,j)*fn(c,j)
               end if
            end if
          end if
       enddo
    end do

    end associate

  end subroutine SetRHSVec_SoilUrbanRoad

  !-----------------------------------------------------------------------
  subroutine SetRHSVec_SoilNonUrban(bounds, num_nolakec, filter_nolakec, &
       hs_top_snow, hs_soil, hs_top, dhsdT, sabg_lyr_col, fact, fn, fn_h2osfc, c_h2osfc, rt)

    !
    ! !DESCRIPTION:
    ! Sets up RHS vector corresponding to soil layers.
    !
    ! !USES:
    use clmtype
    use shr_infnan_mod  , only : nan => shr_infnan_nan, assignment(=)
    use clm_varcon      , only : cnfac
    use column_varcon  , only : icol_roof, icol_sunwall, icol_shadewall
    use clm_varpar     , only : nlevsno, nlevgrnd, nlevurb
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds                                     ! bounds
    integer , intent(in)  :: num_nolakec                                        ! number of column non-lake points in column filter
    integer , intent(in)  :: filter_nolakec(:)                                  ! column filter for non-lake points
    real(r8), intent(in)  :: hs_top_snow(bounds%begc: )                         ! heat flux on top snow layer [W/m2]
    real(r8), intent(in)  :: hs_soil(bounds%begc: )                             ! heat flux on soil [W/m2]
    real(r8), intent(in)  :: hs_top(bounds%begc: )                              ! net energy flux into surface layer (col) [W/m2]
    real(r8), intent(in)  :: dhsdT(bounds%begc: )                               ! temperature derivative of "hs" [col]
    real(r8), intent(in)  :: sabg_lyr_col(bounds%begc:, -nlevsno+1: )           ! absorbed solar radiation (col,lyr) [W/m2]
    real(r8), intent(in)  :: fact( bounds%begc: , -nlevsno+1: )                 ! used in computing tridiagonal matrix [col, lev]
    real(r8), intent(in)  :: fn (bounds%begc: ,-nlevsno+1: )                    ! heat diffusion through the layer interface [W/m2]
    real(r8), intent(in)  :: fn_h2osfc (bounds%begc: )                          ! heat diffusion through standing-water/soil interface [W/m2]
    real(r8), intent(in)  :: c_h2osfc( bounds%begc: )                           ! heat capacity of surface water [col]
    real(r8), intent(inout) :: rt(bounds%begc: ,1: )                            ! rhs vector entries
    !
    ! !LOCAL VARIABLES:
    integer  :: j,c,l                                                           ! indices
    integer  :: fc                                                              ! lake filtered column indices
    !-----------------------------------------------------------------------

    ! Enforce expected array sizes
    SHR_ASSERT_ALL((ubound(hs_soil)      == (/bounds%endc/)),           errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(dhsdT)        == (/bounds%endc/)),           errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(sabg_lyr_col) == (/bounds%endc, 1/)),        errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(fact)         == (/bounds%endc, nlevgrnd/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(fn)           == (/bounds%endc, nlevgrnd/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(fn_h2osfc)    == (/bounds%endc/)),           errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(c_h2osfc)     == (/bounds%endc/)),           errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(rt)           == (/bounds%endc, nlevgrnd/)), errMsg(__FILE__, __LINE__))

   associate(&
   snl                       =>    cps%snl                       , & ! Input:  [integer (:)]  number of snow layers
   urbpoi                    =>    lun%urbpoi                    , & ! Input:  [logical (:)]  true => landunit is an urban point
   clandunit                 =>    col%landunit                  , & ! Input:  [integer (:)]  column's landunit
   frac_sno_eff              =>    cps%frac_sno_eff              , & ! Input:  [real(r8) (:)]  eff. fraction of ground covered by snow (0 to 1)
   frac_h2osfc               =>    cps%frac_h2osfc               , & ! Input:  [real(r8) (:)]  fraction of ground covered by surface water (0 to 1)
   ctype                     =>    col%itype                     , & ! Input:  [integer (:)]  column type
   zi                        =>    cps%zi                        , & ! Input:  [real(r8) (:,:)]  interface level below a "z" level (m)
   z                         =>    cps%z                         , & ! Input:  [real(r8) (:,:)]  layer thickness (m)
   t_soisno                  =>    ces%t_soisno                  , & ! Input:  [real(r8) (:,:)]  soil temperature (Kelvin)
   begc                      =>    bounds%begc                   , & ! Input:  [integer ] beginning column index
   endc                      =>    bounds%endc                     & ! Input:  [integer ] ending column index
   )

    !
    ! non-urban columns --------------------------------------------------------------
    !
    do j = 1,nlevgrnd
       do fc = 1,num_nolakec
          c = filter_nolakec(fc)
          l = clandunit(c)
          if (.not. urbpoi(l)) then
             if (j == snl(c)+1) then
                   rt(c,j) = t_soisno(c,j) +  fact(c,j)*( hs_top_snow(c) &
                        - dhsdT(c)*t_soisno(c,j) + cnfac*fn(c,j) )
             else if (j == 1) then
                ! this is the snow/soil interface layer
                rt(c,j) = t_soisno(c,j) + fact(c,j) &
                     *((1._r8-frac_sno_eff(c))*(hs_soil(c) - dhsdT(c)*t_soisno(c,j)) &
                     + cnfac*(fn(c,j) - frac_sno_eff(c) * fn(c,j-1)))

                rt(c,j) = rt(c,j) +  frac_sno_eff(c)*fact(c,j)*sabg_lyr_col(c,j)

             else if (j <= nlevgrnd-1) then
                rt(c,j) = t_soisno(c,j) + cnfac*fact(c,j)*( fn(c,j) - fn(c,j-1) )

             else if (j == nlevgrnd) then
                rt(c,j) = t_soisno(c,j) - cnfac*fact(c,j)*fn(c,j-1) + fact(c,j)*fn(c,j)
             end if
          end if
       enddo
    end do

    end associate

  end subroutine SetRHSVec_SoilNonUrban

  !-----------------------------------------------------------------------
  subroutine SetRHSVec_Soil_StandingSurfaceWater(bounds, num_nolakec, filter_nolakec, &
       hs_top_snow, hs_soil, hs_top, dhsdT, sabg_lyr_col, fact, fn, fn_h2osfc, c_h2osfc, rt)

    !
    ! !DESCRIPTION:
    ! Sets up RHS vector corresponding to soil layers.
    !
    ! !USES:
    use clmtype
    use shr_infnan_mod  , only : nan => shr_infnan_nan, assignment(=)
    use clm_varcon      , only : cnfac
    use column_varcon  , only : icol_roof, icol_sunwall, icol_shadewall
    use clm_varpar     , only : nlevsno, nlevgrnd, nlevurb
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds                                     ! bounds
    integer , intent(in)  :: num_nolakec                                        ! number of column non-lake points in column filter
    integer , intent(in)  :: filter_nolakec(:)                                  ! column filter for non-lake points
    real(r8), intent(in)  :: hs_top_snow(bounds%begc: )                         ! heat flux on top snow layer [W/m2]
    real(r8), intent(in)  :: hs_soil(bounds%begc: )                             ! heat flux on soil [W/m2]
    real(r8), intent(in)  :: hs_top(bounds%begc: )                              ! net energy flux into surface layer (col) [W/m2]
    real(r8), intent(in)  :: dhsdT(bounds%begc: )                               ! temperature derivative of "hs" [col]
    real(r8), intent(in)  :: sabg_lyr_col(bounds%begc:, -nlevsno+1: )           ! absorbed solar radiation (col,lyr) [W/m2]
    real(r8), intent(in)  :: fact( bounds%begc: , -nlevsno+1: )                 ! used in computing tridiagonal matrix [col, lev]
    real(r8), intent(in)  :: fn (bounds%begc: ,-nlevsno+1: )                    ! heat diffusion through the layer interface [W/m2]
    real(r8), intent(in)  :: fn_h2osfc (bounds%begc: )                          ! heat diffusion through standing-water/soil interface [W/m2]
    real(r8), intent(in)  :: c_h2osfc( bounds%begc: )                           ! heat capacity of surface water [col]
    real(r8), intent(inout) :: rt(bounds%begc: ,1: )                            ! rhs vector entries
    !
    ! !LOCAL VARIABLES:
    integer  :: j,c,l                                                           ! indices
    integer  :: fc                                                              ! lake filtered column indices
    !-----------------------------------------------------------------------

    ! Enforce expected array sizes
    SHR_ASSERT_ALL((ubound(hs_soil)      == (/bounds%endc/)),           errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(dhsdT)        == (/bounds%endc/)),           errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(sabg_lyr_col) == (/bounds%endc, 1/)),        errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(fact)         == (/bounds%endc, nlevgrnd/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(fn)           == (/bounds%endc, nlevgrnd/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(fn_h2osfc)    == (/bounds%endc/)),           errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(c_h2osfc)     == (/bounds%endc/)),           errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(rt)           == (/bounds%endc, nlevgrnd/)), errMsg(__FILE__, __LINE__))

   associate(&
   snl                       =>    cps%snl                       , & ! Input:  [integer (:)]  number of snow layers
   urbpoi                    =>    lun%urbpoi                    , & ! Input:  [logical (:)]  true => landunit is an urban point
   clandunit                 =>    col%landunit                  , & ! Input:  [integer (:)]  column's landunit
   frac_sno_eff              =>    cps%frac_sno_eff              , & ! Input:  [real(r8) (:)]  eff. fraction of ground covered by snow (0 to 1)
   frac_h2osfc               =>    cps%frac_h2osfc               , & ! Input:  [real(r8) (:)]  fraction of ground covered by surface water (0 to 1)
   ctype                     =>    col%itype                     , & ! Input:  [integer (:)]  column type
   zi                        =>    cps%zi                        , & ! Input:  [real(r8) (:,:)]  interface level below a "z" level (m)
   z                         =>    cps%z                         , & ! Input:  [real(r8) (:,:)]  layer thickness (m)
   t_soisno                  =>    ces%t_soisno                  , & ! Input:  [real(r8) (:,:)]  soil temperature (Kelvin)
   begc                      =>    bounds%begc                   , & ! Input:  [integer ] beginning column index
   endc                      =>    bounds%endc                     & ! Input:  [integer ] ending column index
   )

    !
    ! surface water  -----------------------------------------------------------------
    !
    do fc = 1,num_nolakec
       c = filter_nolakec(fc)
       if ( frac_h2osfc(c) /= 0.0_r8 )then
          rt(c,1)=rt(c,1) &
               -frac_h2osfc(c)*fact(c,1)*((hs_soil(c) - dhsdT(c)*t_soisno(c,1)) &
               +cnfac*fn_h2osfc(c))
       end if
    end do

    end associate

  end subroutine SetRHSVec_Soil_StandingSurfaceWater

  !-----------------------------------------------------------------------
  subroutine SetMatrix(bounds, num_nolakec, filter_nolakec, dtime, nband, &
       dhsdT, tk, tk_h2osfc, fact, c_h2osfc, dz_h2osfc, bmatrix)

    !
    ! !DESCRIPTION:
    ! Setup the matrix for the numerical solution of temperature for snow,
    ! standing surface water and soil layers.
    !
    !
    !           |===========|===========|===========|
    !           |   Snow    |           | Snow-Soil |
    !           !===========|===========|===========|
    ! bmatrix = |           |    SSW    |  SSW-Soil |
    !           !===========|===========|===========|
    !           ! Soil-Snow | Soil-SSW  |    Soil   |
    !           !===========|===========|===========|
    !
    !
    !
    ! !USES:
    use clmtype
    use shr_infnan_mod  , only : nan => shr_infnan_nan, assignment(=)
    use clm_varcon      , only : cnfac
    use column_varcon  , only : icol_roof, icol_sunwall, icol_shadewall
    use clm_varpar     , only : nlevsno, nlevgrnd, nlevurb
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds                                    ! bounds
    integer , intent(in)  :: num_nolakec                                       ! number of column non-lake points in column filter
    integer , intent(in)  :: filter_nolakec(:)                                 ! column filter for non-lake points
    real(r8), intent(in)  :: dtime                                             ! land model time step (sec)
    integer , intent(in)  :: nband                                             ! number of bands of the tridigonal matrix
    real(r8), intent(in)  :: dhsdT(bounds%begc: )                              ! temperature derivative of "hs" [col]
    real(r8), intent(in)  :: tk(bounds%begc: ,-nlevsno+1: )                    ! thermal conductivity [W/(m K)]
    real(r8), intent(in)  :: tk_h2osfc(bounds%begc: )                          ! thermal conductivity [W/(m K)]
    real(r8), intent(in)  :: fact( bounds%begc: , -nlevsno+1: )                ! used in computing tridiagonal matrix [col, lev]
    real(r8), intent(in)  :: c_h2osfc( bounds%begc: )                          ! heat capacity of surface water [col]
    real(r8), intent(in)  :: dz_h2osfc(bounds%begc: )                          ! Thickness of standing water [m]
    real(r8), intent(out) :: bmatrix(bounds%begc: , 1:,-nlevsno: )             ! matrix for numerical solution of temperature
    !
    ! !LOCAL VARIABLES:
    integer  :: j,c                                                            ! indices
    integer  :: fc                                                             ! lake filtered column indices
    real(r8) :: at (bounds%begc:bounds%endc,-nlevsno+1:nlevgrnd)               ! "a" vector for tridiagonal matrix
    real(r8) :: bt (bounds%begc:bounds%endc,-nlevsno+1:nlevgrnd)               ! "b" vector for tridiagonal matrix
    real(r8) :: ct (bounds%begc:bounds%endc,-nlevsno+1:nlevgrnd)               ! "c" vector for tridiagonal matrix
    real(r8) :: dzm                                                            ! used in computing tridiagonal matrix
    real(r8) :: dzp                                                            ! used in computing tridiagonal matrix
    real(r8) :: bmatrix_snow(bounds%begc:bounds%endc,nband,-nlevsno:-1      )  ! block-diagonal matrix for snow layers
    real(r8) :: bmatrix_ssw(bounds%begc:bounds%endc,nband,       0:0       )   ! block-diagonal matrix for standing surface water
    real(r8) :: bmatrix_soil(bounds%begc:bounds%endc,nband,       1:nlevgrnd)  ! block-diagonal matrix for soil layers
    real(r8) :: bmatrix_snow_soil(bounds%begc:bounds%endc,nband,-1:-1)         ! off-diagonal matrix for snow-soil interaction
    real(r8) :: bmatrix_ssw_soil(bounds%begc:bounds%endc,nband, 0:0 )          ! off-diagonal matrix for standing surface water-soil interaction
    real(r8) :: bmatrix_soil_snow(bounds%begc:bounds%endc,nband, 1:1 )         ! off-diagonal matrix for soil-snow interaction
    real(r8) :: bmatrix_soil_ssw(bounds%begc:bounds%endc,nband, 1:1 )          ! off-diagonal matrix for soil-standing surface water interaction
    !-----------------------------------------------------------------------

    ! Enforce expected array sizes
    SHR_ASSERT_ALL((ubound(dhsdT)     == (/bounds%endc/)),                  errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(tk)        == (/bounds%endc, nlevgrnd/)),        errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(tk_h2osfc) == (/bounds%endc/)),                  errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(fact)      == (/bounds%endc, nlevgrnd/)),        errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(c_h2osfc)  == (/bounds%endc/)),                  errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(dz_h2osfc) == (/bounds%endc/)),                  errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(bmatrix)   == (/bounds%endc, nband, nlevgrnd/)), errMsg(__FILE__, __LINE__))

   associate(&
   frac_sno_eff              =>    cps%frac_sno_eff              , & ! Input:  [real(r8) (:)]  eff. fraction of ground covered by snow (0 to 1)
   frac_h2osfc               =>    cps%frac_h2osfc               , & ! Input:  [real(r8) (:)]  fraction of ground covered by surface water (0 to 1)
   h2osfc                    =>    cws%h2osfc                    , & ! Input:  [real(r8) (:)]  surface water (mm)
   t_h2osfc                  =>    ces%t_h2osfc                  , & ! Input:  [real(r8) (:)]  surface water temperature
   ctype                     =>    col%itype                     , & ! Input:  [integer (:)]  column type
   snl                       =>    cps%snl                       , & ! Input:  [integer (:)]  number of snow layers
   t_grnd                    =>    ces%t_grnd                    , & ! Input:  [real(r8) (:)]  ground surface temperature [K]
   zi                        =>    cps%zi                        , & ! Input:  [real(r8) (:,:)]  interface level below a "z" level (m)
   z                         =>    cps%z                         , & ! Input:  [real(r8) (:,:)]  layer thickness (m)
   t_soisno                  =>    ces%t_soisno                  , & ! Input:  [real(r8) (:,:)]  soil temperature (Kelvin)
   begc                      =>    bounds%begc                   , & ! Input:  [integer ] beginning column index
   endc                      =>    bounds%endc                     & ! Input:  [integer ] ending column index
   )

    ! Assemble smaller matrices

    call SetMatrix_Snow(bounds, num_nolakec, filter_nolakec, nband, &
         dhsdT( begc:endc ),                                        &
         tk( begc:endc, -nlevsno+1: ),                              &
         fact( begc:endc, -nlevsno+1: ),                            &
         bmatrix_snow( begc:endc, 1:, -nlevsno: ))

    call SetMatrix_Snow_Soil(bounds, num_nolakec, filter_nolakec, nband, &
         tk( begc:endc, -nlevsno+1: ),                                   &
         fact( begc:endc, -nlevsno+1: ),                                 &
         bmatrix_snow_soil( begc:endc, 1:, -1: ))

    call SetMatrix_Soil(bounds, num_nolakec, filter_nolakec, nband, &
         dhsdT( begc:endc ),                                        &
         tk( begc:endc, -nlevsno+1: ),                              &
         tk_h2osfc( begc:endc ),                                    &
         dz_h2osfc( begc:endc ),                                    &
         fact( begc:endc, -nlevsno+1: ),                            &
         bmatrix_soil( begc:endc, 1:, 1: ))

    call SetMatrix_Soil_Snow(bounds, num_nolakec, filter_nolakec, nband, &
         tk( begc:endc, -nlevsno+1: ),                                   &
         fact( begc:endc, -nlevsno+1: ),                                 &
         bmatrix_soil_snow( begc:endc, 1:, 1: ))

    call SetMatrix_StandingSurfaceWater(bounds, num_nolakec, filter_nolakec, dtime, nband, &
         dhsdT( begc:endc ),                                                               &
         tk( begc:endc, -nlevsno+1: ),                                                     &
         tk_h2osfc( begc:endc ),                                                           &
         fact( begc:endc, -nlevsno+1: ),                                                   &
         c_h2osfc( begc:endc ),                                                            &
         dz_h2osfc( begc:endc ),                                                           &
         bmatrix_ssw( begc:endc, 1:, 0: ))

    call SetMatrix_StandingSurfaceWater_Soil(bounds, num_nolakec, filter_nolakec, dtime, nband, &
         tk( begc:endc, -nlevsno+1: ),                                                          &
         tk_h2osfc( begc:endc ),                                                                &
         fact( begc:endc, -nlevsno+1: ),                                                        &
         c_h2osfc( begc:endc ),                                                                 &
         dz_h2osfc( begc:endc ),                                                                &
         bmatrix_ssw_soil( begc:endc, 1:, 0: ))

    call SetMatrix_Soil_StandingSurfaceWater(bounds, num_nolakec, filter_nolakec, nband, &
         tk_h2osfc( begc:endc ),                                                         &
         fact( begc:endc, -nlevsno+1: ),                                                 &
         dz_h2osfc( begc:endc ),                                                         &
         bmatrix_soil_ssw( begc:endc, 1:, 1: ))

    call AssembleMatrixFromSubmatrices(bounds, num_nolakec, filter_nolakec, nband, &
         bmatrix_snow( begc:endc, 1:, -nlevsno: ),                                 &
         bmatrix_ssw( begc:endc, 1:, 0: ),                                         &
         bmatrix_soil( begc:endc, 1:, 1: ),                                        &
         bmatrix_snow_soil( begc:endc, 1:, -1: ),                                  &
         bmatrix_ssw_soil( begc:endc, 1:, 0: ),                                    &
         bmatrix_soil_snow( begc:endc, 1:, 1: ),                                   &
         bmatrix_soil_ssw( begc:endc, 1:, 1: ),                                    &
         bmatrix( begc:endc, 1:, -nlevsno: ))

    end associate

  end subroutine SetMatrix

  !-----------------------------------------------------------------------
  subroutine AssembleMatrixFromSubmatrices(bounds, num_nolakec, filter_nolakec, nband, &
       bmatrix_snow, bmatrix_ssw, bmatrix_soil, bmatrix_snow_soil, &
       bmatrix_ssw_soil, bmatrix_soil_snow, bmatrix_soil_ssw, bmatrix)

    !
    ! !DESCRIPTION:
    ! Assemble the full matrix from submatrices.
    !
    ! Non-zero pattern of bmatrix:
    !
    !        SNOW-LAYERS
    !            |
    !            |  STANDING-SURFACE-WATER
    !            |         |
    !            |         |              SOIL-LAYERS
    !            |         |                  |
    !            v         v                  v
    !
    !      -5 -4 -3 -2 -1| 0| 1  2  3  4  5  6  7  8  9 10 11 12 13 14 15
    !      ==============================================================
    !  -5 | x  x         |  |                                            |
    !  -4 | x  x  x      |  |                                            |
    !  -3 |    x  x  x   |  |                                            |
    !  -2 |       x  x  x|  |                                            |
    !  -1 |          x  x|  | x                                          |
    !      ==============================================================
    !   0 |              | x| x                                          |
    !      ==============================================================
    !   1 |             x| x| x  x                                       |
    !   2 |              |  | x  x  x                                    |
    !   3 |              |  |    x  x  x                                 |
    !   4 |              |  |       x  x  x                              |
    !   5 |              |  |          x  x  x                           |
    !   6 |              |  |             x  x  x                        |
    !   7 |              |  |                x  x  x                     |
    !   8 |              |  |                   x  x  x                  |
    !   9 |              |  |                      x  x  x               |
    !  10 |              |  |                         x  x  x            |
    !  11 |              |  |                            x  x  x         |
    !  12 |              |  |                               x  x  x      |
    !  13 |              |  |                                  x  x  x   |
    !  14 |              |  |                                     x  x  x|
    !  15 |              |  |                                        x  x|
    !      ==============================================================
    !
    !
    ! !USES:
    use clmtype
    use shr_infnan_mod  , only : nan => shr_infnan_nan, assignment(=)
    use clm_varcon      , only : cnfac
    use column_varcon  , only : icol_roof, icol_sunwall, icol_shadewall
    use clm_varpar     , only : nlevsno, nlevgrnd, nlevurb
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds                                                 ! bounds
    integer , intent(in)  :: num_nolakec                                                    ! number of column non-lake points in column filter
    integer , intent(in)  :: filter_nolakec(:)                                              ! column filter for non-lake points
    integer , intent(in)  :: nband                                                          ! number of bands of the tridigonal matrix
    real(r8), intent(in)  :: bmatrix_snow(bounds%begc:bounds%endc,nband,-nlevsno:-1      )  ! block-diagonal matrix for snow layers
    real(r8), intent(in)  :: bmatrix_ssw(bounds%begc:bounds%endc,nband,       0:0       )   ! block-diagonal matrix for standing surface water
    real(r8), intent(in)  :: bmatrix_soil(bounds%begc:bounds%endc,nband,       1:nlevgrnd)  ! block-diagonal matrix for soil layers
    real(r8), intent(in)  :: bmatrix_snow_soil(bounds%begc:bounds%endc,nband,-1:-1)         ! off-diagonal matrix for snow-soil interaction
    real(r8), intent(in)  :: bmatrix_ssw_soil(bounds%begc:bounds%endc,nband, 0:0 )          ! off-diagonal matrix for standing surface water-soil interaction
    real(r8), intent(in)  :: bmatrix_soil_snow(bounds%begc:bounds%endc,nband, 1:1 )         ! off-diagonal matrix for soil-snow interaction
    real(r8), intent(in)  :: bmatrix_soil_ssw(bounds%begc:bounds%endc,nband, 1:1 )          ! off-diagonal matrix for soil-standing surface water interaction
    real(r8), intent(out) :: bmatrix(bounds%begc: , 1:,-nlevsno: )                          ! full matrix used in numerical solution of temperature
    !
    ! !LOCAL VARIABLES:
    integer  :: j,c                                                                         ! indices
    integer  :: fc                                                                          ! lake filtered column indices
    !-----------------------------------------------------------------------

    ! Enforce expected array sizes
    SHR_ASSERT_ALL((ubound(bmatrix_snow)        == (/bounds%endc, nband, -1/)),       errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(bmatrix_ssw)         == (/bounds%endc, nband, 0/)),        errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(bmatrix_soil)        == (/bounds%endc, nband, nlevgrnd/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(bmatrix_snow_soil)   == (/bounds%endc, nband, -1/)),       errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(bmatrix_ssw_soil)    == (/bounds%endc, nband, 0/)),        errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(bmatrix_soil_snow)   == (/bounds%endc, nband, 1/)),        errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(bmatrix_soil_ssw)    == (/bounds%endc, nband, 1/)),        errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(bmatrix)             == (/bounds%endc, nband, nlevgrnd/)), errMsg(__FILE__, __LINE__))

   associate(&
   frac_sno_eff              =>    cps%frac_sno_eff              , & ! Input:  [real(r8) (:)]  eff. fraction of ground covered by snow (0 to 1)
   frac_h2osfc               =>    cps%frac_h2osfc               , & ! Input:  [real(r8) (:)]  fraction of ground covered by surface water (0 to 1)
   h2osfc                    =>    cws%h2osfc                    , & ! Input:  [real(r8) (:)]  surface water (mm)
   t_h2osfc                  =>    ces%t_h2osfc                  , & ! Input:  [real(r8) (:)]  surface water temperature
   ctype                     =>    col%itype                     , & ! Input:  [integer (:)]  column type
   snl                       =>    cps%snl                       , & ! Input:  [integer (:)]  number of snow layers
   t_grnd                    =>    ces%t_grnd                    , & ! Input:  [real(r8) (:)]  ground surface temperature [K]
   zi                        =>    cps%zi                        , & ! Input:  [real(r8) (:,:)]  interface level below a "z" level (m)
   z                         =>    cps%z                         , & ! Input:  [real(r8) (:,:)]  layer thickness (m)
   t_soisno                  =>    ces%t_soisno                  , & ! Input:  [real(r8) (:,:)]  soil temperature (Kelvin)
   begc                      =>    bounds%begc                   , & ! Input:  [integer ] beginning column index
   endc                      =>    bounds%endc                     & ! Input:  [integer ] ending column index
   )

    ! Assemble the full matrix

    bmatrix(begc:endc, :, :) = 0.0_r8
    do fc = 1,num_nolakec
       c = filter_nolakec(fc)

       ! Snow
       bmatrix(c,2:3,-5   ) = bmatrix_snow(c,2:3,-5   )
       bmatrix(c,2:4,-4:-2) = bmatrix_snow(c,2:4,-4:-2)
       bmatrix(c,3:4,-1   ) = bmatrix_snow(c,3:4,-1   )

       ! Snow-Soil
       bmatrix(c,1,-1) = bmatrix_snow_soil(c,1,-1)

       ! StandingSurfaceWater
       bmatrix(c,3,0) = bmatrix_ssw(c,3,0)

       ! StandingSurfaceWater-Soil
       bmatrix(c,2,0) = bmatrix_ssw_soil(c,2,0)

       ! Soil
       bmatrix(c,2:3,1           )  = bmatrix_soil(c,2:3,1           )
       bmatrix(c,2:4,2:nlevgrnd-1)  = bmatrix_soil(c,2:4,2:nlevgrnd-1)
       bmatrix(c,3:4,nlevgrnd    )  = bmatrix_soil(c,3:4,nlevgrnd    )

       ! Soil-Snow
       bmatrix(c,5,1)  = bmatrix_soil_snow(c,5,1)

       ! Soil-StandingSurfaceWater
       bmatrix(c,4,1)  = bmatrix_soil_ssw(c,4,1)

    end do

    end associate

  end subroutine AssembleMatrixFromSubmatrices

  !-----------------------------------------------------------------------
  subroutine SetMatrix_Snow(bounds, num_nolakec, filter_nolakec, nband, &
       dhsdT, tk, fact, bmatrix_snow)

    !
    ! !DESCRIPTION:
    ! Setup the matrix entries correspodning to internal snow layers
    !
    ! !USES:
    use clmtype
    use clm_varcon     , only : cnfac
    use column_varcon  , only : icol_roof, icol_sunwall, icol_shadewall
    use clm_varpar     , only : nlevsno, nlevgrnd, nlevurb
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds                               ! bounds
    integer , intent(in)  :: num_nolakec                                  ! number of column non-lake points in column filter
    integer , intent(in)  :: filter_nolakec(:)                            ! column filter for non-lake points
    integer , intent(in)  :: nband                                        ! number of bands of the tridigonal matrix
    real(r8), intent(in)  :: dhsdT(bounds%begc: )                         ! temperature derivative of "hs" [col]
    real(r8), intent(in)  :: tk(bounds%begc: ,-nlevsno+1: )               ! thermal conductivity [W/(m K)]
    real(r8), intent(in)  :: fact( bounds%begc: , -nlevsno+1: )           ! used in computing tridiagonal matrix [col, lev]
    real(r8), intent(out) :: bmatrix_snow(bounds%begc: , 1:, -nlevsno: )  ! matrix enteries
    !-----------------------------------------------------------------------

    ! Enforce expected array sizes
    SHR_ASSERT_ALL((ubound(dhsdT)          == (/bounds%endc/)),             errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(tk)             == (/bounds%endc, nlevgrnd/)),   errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(fact)           == (/bounds%endc, nlevgrnd/)),   errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(bmatrix_snow)   == (/bounds%endc, nband, -1/)),  errMsg(__FILE__, __LINE__))

   associate(&
   begc                      =>    bounds%begc                   , & ! Input:  [integer ] beginning column index
   endc                      =>    bounds%endc                     & ! Input:  [integer ] ending column index
   )

    ! Initialize
    bmatrix_snow(begc:endc, :, :) = 0.0_r8

    call SetMatrix_SnowUrban(bounds, num_nolakec, filter_nolakec, nband, &
         dhsdT( begc:endc ),                                             &
         tk( begc:endc, -nlevsno+1: ),                                   &
         fact( begc:endc, -nlevsno+1: ),                                 &
         bmatrix_snow( begc:endc, 1:, -nlevsno: ))

    call SetMatrix_SnowNonUrban(bounds, num_nolakec, filter_nolakec, nband, &
         dhsdT( begc:endc ),                                                &
         tk( begc:endc, -nlevsno+1: ),                                      &
         fact( begc:endc, -nlevsno+1: ),                                    &
         bmatrix_snow( begc:endc, 1:, -nlevsno: ))

  end associate

  end subroutine SetMatrix_Snow

  !-----------------------------------------------------------------------
  subroutine SetMatrix_SnowUrban(bounds, num_nolakec, filter_nolakec, nband, &
       dhsdT, tk, fact, bmatrix_snow)

    !
    ! !DESCRIPTION:
    ! Setup the matrix entries correspodning to internal snow layers for
    ! urban soil columns
    !
    ! !USES:
    use clmtype
    use clm_varcon     , only : cnfac
    use column_varcon  , only : icol_roof, icol_sunwall, icol_shadewall
    use clm_varpar     , only : nlevsno, nlevgrnd, nlevurb
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds                                 ! bounds
    integer , intent(in)  :: num_nolakec                                    ! number of column non-lake points in column filter
    integer , intent(in)  :: filter_nolakec(:)                              ! column filter for non-lake points
    integer , intent(in)  :: nband                                          ! number of bands of the tridigonal matrix
    real(r8), intent(in)  :: dhsdT(bounds%begc: )                           ! temperature derivative of "hs" [col]
    real(r8), intent(in)  :: tk(bounds%begc: ,-nlevsno+1: )                 ! thermal conductivity [W/(m K)]
    real(r8), intent(in)  :: fact( bounds%begc: , -nlevsno+1: )             ! used in computing tridiagonal matrix [col, lev]
    real(r8), intent(inout) :: bmatrix_snow(bounds%begc: , 1:, -nlevsno: )  ! matrix enteries
    !-----------------------------------------------------------------------

    ! Enforce expected array sizes
    SHR_ASSERT_ALL((ubound(dhsdT)          == (/bounds%endc/)),                  errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(tk)             == (/bounds%endc, nlevgrnd/)),        errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(fact)           == (/bounds%endc, nlevgrnd/)),        errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(bmatrix_snow)   == (/bounds%endc, nband, -1/)), errMsg(__FILE__, __LINE__))

   associate(& 
   begc                      =>    bounds%begc                   , & ! Input:  [integer ] beginning column index
   endc                      =>    bounds%endc                     & ! Input:  [integer ] ending column index
   )

    call SetMatrix_SnowUrbanNonRoad(bounds, num_nolakec, filter_nolakec, nband, &
         dhsdT( begc:endc ),                                                    &
         tk( begc:endc, -nlevsno+1: ),                                          &
         fact( begc:endc, -nlevsno+1: ),                                        &
         bmatrix_snow( begc:endc, 1:, -nlevsno: ))

    call SetMatrix_SnowUrbanRoad(bounds, num_nolakec, filter_nolakec, nband, &
         dhsdT( begc:endc ),                                                 &
         tk( begc:endc, -nlevsno+1: ),                                       &
         fact( begc:endc, -nlevsno+1: ),                                     &
         bmatrix_snow( begc:endc, 1:, -nlevsno: ))

  end associate

  end subroutine SetMatrix_SnowUrban

  !-----------------------------------------------------------------------
  subroutine SetMatrix_SnowUrbanNonRoad(bounds, num_nolakec, filter_nolakec, nband, &
       dhsdT, tk, fact, bmatrix_snow)

    !
    ! !DESCRIPTION:
    ! Setup the matrix entries correspodning to internal snow layers for
    ! urban sunwall/shadewall/roof columns
    !
    ! !USES:
    use clmtype
    use clm_varcon     , only : cnfac
    use column_varcon  , only : icol_roof, icol_sunwall, icol_shadewall
    use clm_varpar     , only : nlevsno, nlevgrnd, nlevurb
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds                                 ! bounds
    integer , intent(in)  :: num_nolakec                                    ! number of column non-lake points in column filter
    integer , intent(in)  :: filter_nolakec(:)                              ! column filter for non-lake points
    integer , intent(in)  :: nband                                          ! number of bands of the tridigonal matrix
    real(r8), intent(in)  :: dhsdT(bounds%begc: )                           ! temperature derivative of "hs" [col]
    real(r8), intent(in)  :: tk(bounds%begc: ,-nlevsno+1: )                 ! thermal conductivity [W/(m K)]
    real(r8), intent(in)  :: fact( bounds%begc: , -nlevsno+1: )             ! used in computing tridiagonal matrix [col, lev]
    real(r8), intent(inout) :: bmatrix_snow(bounds%begc: , 1:, -nlevsno: )  ! matrix enteries
    !
    ! !LOCAL VARIABLES:
    integer  :: j,c,l                                                       ! indices
    integer  :: fc                                                          ! lake filtered column indices
    real(r8) :: dzm                                                         ! used in computing tridiagonal matrix
    real(r8) :: dzp                                                         ! used in computing tridiagonal matrix
    !-----------------------------------------------------------------------

    ! Enforce expected array sizes
    SHR_ASSERT_ALL((ubound(dhsdT)          == (/bounds%endc/)),            errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(tk)             == (/bounds%endc, nlevgrnd/)),  errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(fact)           == (/bounds%endc, nlevgrnd/)),  errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(bmatrix_snow)   == (/bounds%endc, nband, -1/)), errMsg(__FILE__, __LINE__))

   associate(& 
   snl                       =>    cps%snl                       , & ! Input:  [integer (:)]  number of snow layers                    
   ctype                     =>    col%itype                     , & ! Input:  [integer (:)]  column type
   urbpoi                    =>    lun%urbpoi                    , & ! Input:  [logical (:)]  true => landunit is an urban point
   clandunit                 =>    col%landunit                  , & ! Input:  [integer (:)]  column's landunit
   z                         =>    cps%z                         , & ! Input:  [real(r8) (:,:)]  layer thickness (m)
   zi                        =>    cps%zi                        , & ! Input:  [real(r8) (:,:)]  interface level below a "z" level (m)
   begc                      =>    bounds%begc                   , & ! Input:  [integer ] beginning column index
   endc                      =>    bounds%endc                     & ! Input:  [integer ] ending column index
   )

    !
    ! urban non-road columns ---------------------------------------------------------
    !
    do j = -nlevsno+1,0
       do fc = 1,num_nolakec
          c = filter_nolakec(fc)
          l = clandunit(c)
          if (urbpoi(l)) then
             if ((ctype(c) == icol_sunwall .or. ctype(c) == icol_shadewall &
                 .or. ctype(c) == icol_roof)) then
                if (j >= snl(c)+1) then
                   if (j == snl(c)+1) then
                      dzp     = z(c,j+1)-z(c,j)
                      bmatrix_snow(c,4,j-1) = 0._r8
                      bmatrix_snow(c,3,j-1) = 1+(1._r8-cnfac)*fact(c,j)*tk(c,j)/dzp-fact(c,j)*dhsdT(c)
                      if ( j /= 0) then
                         bmatrix_snow(c,2,j-1) =  -(1._r8-cnfac)*fact(c,j)*tk(c,j)/dzp
                      end if
                   else if (j <= nlevurb-1) then
                      dzm     = (z(c,j)-z(c,j-1))
                      dzp     = (z(c,j+1)-z(c,j))
                      bmatrix_snow(c,4,j-1) =   - (1._r8-cnfac)*fact(c,j)* tk(c,j-1)/dzm
                      bmatrix_snow(c,3,j-1) = 1._r8+ (1._r8-cnfac)*fact(c,j)*(tk(c,j)/dzp + tk(c,j-1)/dzm)
                      if (j /= 0) then
                         bmatrix_snow(c,2,j-1) =   - (1._r8-cnfac)*fact(c,j)* tk(c,j)/dzp
                      end if
                   end if
               end if
            end if
          end if
       enddo
    end do

  end associate

  end subroutine SetMatrix_SnowUrbanNonRoad

  !-----------------------------------------------------------------------
  subroutine SetMatrix_SnowUrbanRoad(bounds, num_nolakec, filter_nolakec, nband, &
       dhsdT, tk, fact, bmatrix_snow)

    !
    ! !DESCRIPTION:
    ! Setup the matrix entries correspodning to internal snow layers for
    ! urban road (impervious + pervious) columns
    !
    ! !USES:
    use clmtype
    use clm_varcon     , only : cnfac
    use column_varcon  , only : icol_road_perv, icol_road_imperv
    use clm_varpar     , only : nlevsno, nlevgrnd, nlevurb
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds                                 ! bounds
    integer , intent(in)  :: num_nolakec                                    ! number of column non-lake points in column filter
    integer , intent(in)  :: filter_nolakec(:)                              ! column filter for non-lake points
    integer , intent(in)  :: nband                                          ! number of bands of the tridigonal matrix
    real(r8), intent(in)  :: dhsdT(bounds%begc: )                           ! temperature derivative of "hs" [col]
    real(r8), intent(in)  :: tk(bounds%begc: ,-nlevsno+1: )                 ! thermal conductivity [W/(m K)]
    real(r8), intent(in)  :: fact( bounds%begc: , -nlevsno+1: )             ! used in computing tridiagonal matrix [col, lev]
    real(r8), intent(inout) :: bmatrix_snow(bounds%begc: , 1:, -nlevsno: )  ! matrix enteries
    !
    ! !LOCAL VARIABLES:
    integer  :: j,c,l                                                       ! indices
    integer  :: fc                                                          ! lake filtered column indices
    real(r8) :: dzm                                                         ! used in computing tridiagonal matrix
    real(r8) :: dzp                                                         ! used in computing tridiagonal matrix
    !-----------------------------------------------------------------------

    ! Enforce expected array sizes
    SHR_ASSERT_ALL((ubound(dhsdT)          == (/bounds%endc/)),            errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(tk)             == (/bounds%endc, nlevgrnd/)),  errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(fact)           == (/bounds%endc, nlevgrnd/)),  errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(bmatrix_snow)   == (/bounds%endc, nband, -1/)), errMsg(__FILE__, __LINE__))

   associate(& 
   snl                       =>    cps%snl                       , & ! Input:  [integer (:)]  number of snow layers                    
   ctype                     =>    col%itype                     , & ! Input:  [integer (:)]  column type
   urbpoi                    =>    lun%urbpoi                    , & ! Input:  [logical (:)]  true => landunit is an urban point
   clandunit                 =>    col%landunit                  , & ! Input:  [integer (:)]  column's landunit
   z                         =>    cps%z                         , & ! Input:  [real(r8) (:,:)]  layer thickness (m)
   zi                        =>    cps%zi                        , & ! Input:  [real(r8) (:,:)]  interface level below a "z" level (m)
   begc                      =>    bounds%begc                   , & ! Input:  [integer ] beginning column index
   endc                      =>    bounds%endc                     & ! Input:  [integer ] ending column index
   )

    !
    ! urban road columns -------------------------------------------------------------
    !
    do j = -nlevsno+1,0
       do fc = 1,num_nolakec
          c = filter_nolakec(fc)
          l = clandunit(c)
          if (urbpoi(l)) then
            if (ctype(c) == icol_road_imperv .or. ctype(c) == icol_road_perv) then
               if (j >= snl(c)+1) then
                  if (j == snl(c)+1) then
                     dzp     = z(c,j+1)-z(c,j)
                     bmatrix_snow(c,4,j-1) = 0._r8
                     bmatrix_snow(c,3,j-1) = 1+(1._r8-cnfac)*fact(c,j)*tk(c,j)/dzp-fact(c,j)*dhsdT(c)
                     if ( j /= 0) then
                        bmatrix_snow(c,2,j-1) =  -(1._r8-cnfac)*fact(c,j)*tk(c,j)/dzp
                     end if
                  else if (j <= nlevgrnd-1) then
                     dzm     = (z(c,j)-z(c,j-1))
                     dzp     = (z(c,j+1)-z(c,j))
                     bmatrix_snow(c,4,j-1) =   - (1._r8-cnfac)*fact(c,j)* tk(c,j-1)/dzm
                     bmatrix_snow(c,3,j-1) = 1._r8+ (1._r8-cnfac)*fact(c,j)*(tk(c,j)/dzp + tk(c,j-1)/dzm)
                     if ( j /= 0) then
                        bmatrix_snow(c,2,j-1) =   - (1._r8-cnfac)*fact(c,j)* tk(c,j)/dzp
                     end if
                  end if
               end if
            end if
          end if
       enddo
    end do

  end associate

  end subroutine SetMatrix_SnowUrbanRoad

  !-----------------------------------------------------------------------
  subroutine SetMatrix_SnowNonUrban(bounds, num_nolakec, filter_nolakec, nband, &
       dhsdT, tk, fact, bmatrix_snow)

    !
    ! !DESCRIPTION:
    ! Setup the matrix entries correspodning to internal snow layers for non-urban columns
    !
    ! !USES:
    use clmtype
    use clm_varcon     , only : cnfac
    use column_varcon  , only : icol_roof, icol_sunwall, icol_shadewall
    use clm_varpar     , only : nlevsno, nlevgrnd, nlevurb
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds                                 ! bounds
    integer , intent(in)  :: num_nolakec                                    ! number of column non-lake points in column filter
    integer , intent(in)  :: filter_nolakec(:)                              ! column filter for non-lake points
    integer , intent(in)  :: nband                                          ! number of bands of the tridigonal matrix
    real(r8), intent(in)  :: dhsdT(bounds%begc: )                           ! temperature derivative of "hs" [col]
    real(r8), intent(in)  :: tk(bounds%begc: ,-nlevsno+1: )                 ! thermal conductivity [W/(m K)]
    real(r8), intent(in)  :: fact( bounds%begc: , -nlevsno+1: )             ! used in computing tridiagonal matrix [col, lev]
    real(r8), intent(inout) :: bmatrix_snow(bounds%begc: , 1:, -nlevsno: )  ! matrix enteries
    !
    ! !LOCAL VARIABLES:
    integer  :: j,c,l                                                       ! indices
    integer  :: fc                                                          ! lake filtered column indices
    real(r8) :: dzm                                                         ! used in computing tridiagonal matrix
    real(r8) :: dzp                                                         ! used in computing tridiagonal matrix
    !-----------------------------------------------------------------------

    ! Enforce expected array sizes
    SHR_ASSERT_ALL((ubound(dhsdT)          == (/bounds%endc/)),            errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(tk)             == (/bounds%endc, nlevgrnd/)),  errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(fact)           == (/bounds%endc, nlevgrnd/)),  errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(bmatrix_snow)   == (/bounds%endc, nband, -1/)), errMsg(__FILE__, __LINE__))

   associate(& 
   snl                       =>    cps%snl                       , & ! Input:  [integer (:)]  number of snow layers                    
   ctype                     =>    col%itype                     , & ! Input:  [integer (:)]  column type
   urbpoi                    =>    lun%urbpoi                    , & ! Input:  [logical (:)]  true => landunit is an urban point
   clandunit                 =>    col%landunit                  , & ! Input:  [integer (:)]  column's landunit
   z                         =>    cps%z                         , & ! Input:  [real(r8) (:,:)]  layer thickness (m)
   zi                        =>    cps%zi                        , & ! Input:  [real(r8) (:,:)]  interface level below a "z" level (m)
   begc                      =>    bounds%begc                   , & ! Input:  [integer ] beginning column index
   endc                      =>    bounds%endc                     & ! Input:  [integer ] ending column index
   )

    !
    ! non-urban landunits ------------------------------------------------------------
    !
    do j = -nlevsno+1,0
       do fc = 1,num_nolakec
          c = filter_nolakec(fc)
          l = clandunit(c)
          if (.not. urbpoi(l)) then
             if (j >= snl(c)+1) then
                if (j == snl(c)+1) then
                   dzp     = z(c,j+1)-z(c,j)
                   bmatrix_snow(c,4,j-1) = 0._r8
                   bmatrix_snow(c,3,j-1) = 1+(1._r8-cnfac)*fact(c,j)*tk(c,j)/dzp-fact(c,j)*dhsdT(c)
                   if ( j /= 0) then
                      bmatrix_snow(c,2,j-1) =  -(1._r8-cnfac)*fact(c,j)*tk(c,j)/dzp
                   end if
                else if (j <= nlevgrnd-1) then
                   dzm     = (z(c,j)-z(c,j-1))
                   dzp     = (z(c,j+1)-z(c,j))
                   bmatrix_snow(c,4,j-1) =   - (1._r8-cnfac)*fact(c,j)* tk(c,j-1)/dzm
                   bmatrix_snow(c,3,j-1) = 1._r8+ (1._r8-cnfac)*fact(c,j)*(tk(c,j)/dzp + tk(c,j-1)/dzm)
                   if ( j /= 0) then
                      bmatrix_snow(c,2,j-1) =   - (1._r8-cnfac)*fact(c,j)* tk(c,j)/dzp
                   end if
                end if
             end if
          end if
       enddo
    end do

  end associate

  end subroutine SetMatrix_SnowNonUrban

  !-----------------------------------------------------------------------
  subroutine SetMatrix_Snow_Soil(bounds, num_nolakec, filter_nolakec, nband, &
       tk, fact, bmatrix_snow_soil)

    !
    ! !DESCRIPTION:
    ! Setup the matrix entries correspodning to snow-soil interaction
    !
    ! !USES:
    use clmtype
    use clm_varcon     , only : cnfac
    use column_varcon  , only : icol_roof, icol_sunwall, icol_shadewall
    use clm_varpar     , only : nlevsno, nlevgrnd, nlevurb

    implicit none
    type(bounds_type), intent(in) :: bounds                               ! bounds
    integer , intent(in)  :: num_nolakec                                  ! number of column non-lake points in column filter
    integer , intent(in)  :: filter_nolakec(:)                            ! column filter for non-lake points
    integer , intent(in)  :: nband                                        ! number of bands of the tridigonal matrix
    real(r8), intent(in)  :: tk(bounds%begc: ,-nlevsno+1: )               ! thermal conductivity [W/(m K)]
    real(r8), intent(in)  :: fact( bounds%begc: , -nlevsno+1: )           ! used in computing tridiagonal matrix [col, lev]
    real(r8), intent(out) :: bmatrix_snow_soil(bounds%begc: , 1:,-1: )    ! matrix enteries
    !-----------------------------------------------------------------------

    ! Enforce expected array sizes
    SHR_ASSERT_ALL((ubound(tk)                  == (/bounds%endc, nlevgrnd/)),  errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(fact)                == (/bounds%endc, nlevgrnd/)),  errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(bmatrix_snow_soil)   == (/bounds%endc, nband, -1/)), errMsg(__FILE__, __LINE__))

   associate(& 
   begc                      =>    bounds%begc                   , & ! Input:  [integer ] beginning column index
   endc                      =>    bounds%endc                     & ! Input:  [integer ] ending column index
   )

    ! Initialize
    bmatrix_snow_soil(begc:endc, :, :) = 0.0_r8

    call SetMatrix_Snow_SoilUrban(bounds, num_nolakec, filter_nolakec, nband, &
         tk( begc:endc, -nlevsno+1: ),                                        &
         fact( begc:endc, -nlevsno+1: ),                                      &
         bmatrix_snow_soil( begc:endc, 1:, -1: ))

    call SetMatrix_Snow_SoilNonUrban(bounds, num_nolakec, filter_nolakec, nband, &
         tk( begc:endc, -nlevsno+1: ),                                           &
         fact( begc:endc, -nlevsno+1: ),                                         &
         bmatrix_snow_soil( begc:endc, 1:, -1: ))

    end associate

  end subroutine SetMatrix_Snow_Soil

  !-----------------------------------------------------------------------
  subroutine SetMatrix_Snow_SoilUrban(bounds, num_nolakec, filter_nolakec, nband, &
       tk, fact, bmatrix_snow_soil)
    !
    ! !DESCRIPTION:
    ! Setup the matrix entries correspodning to snow-soil interaction for urban columns
    !
    ! !USES:
    use clmtype
    use clm_varcon     , only : cnfac
    use column_varcon  , only : icol_roof, icol_sunwall, icol_shadewall
    use clm_varpar     , only : nlevsno, nlevgrnd, nlevurb
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds                               ! bounds
    integer , intent(in)  :: num_nolakec                                  ! number of column non-lake points in column filter
    integer , intent(in)  :: filter_nolakec(:)                            ! column filter for non-lake points
    integer , intent(in)  :: nband                                        ! number of bands of the tridigonal matrix
    real(r8), intent(in)  :: tk(bounds%begc: ,-nlevsno+1: )               ! thermal conductivity [W/(m K)]
    real(r8), intent(in)  :: fact( bounds%begc: , -nlevsno+1: )           ! used in computing tridiagonal matrix [col, lev]
    real(r8), intent(inout) :: bmatrix_snow_soil(bounds%begc: , 1:,-1: )  ! matrix enteries
    !-----------------------------------------------------------------------

    ! Enforce expected array sizes
    SHR_ASSERT_ALL((ubound(tk)                  == (/bounds%endc, nlevgrnd/)),  errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(fact)                == (/bounds%endc, nlevgrnd/)),  errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(bmatrix_snow_soil)   == (/bounds%endc, nband, -1/)), errMsg(__FILE__, __LINE__))

   associate(& 
   begc                      =>    bounds%begc                   , & ! Input:  [integer ] beginning column index
   endc                      =>    bounds%endc                     & ! Input:  [integer ] ending column index
   )

    call SetMatrix_Snow_SoilUrbanNonRoad(bounds, num_nolakec, filter_nolakec, nband, &
         tk( begc:endc, -nlevsno+1: ),                                               &
         fact( begc:endc, -nlevsno+1: ),                                             &
         bmatrix_snow_soil( begc:endc, 1:, -1: ))

    call SetMatrix_Snow_SoilUrbanRoad(bounds, num_nolakec, filter_nolakec, nband, &
         tk( begc:endc, -nlevsno+1: ),                                            &
         fact( begc:endc, -nlevsno+1: ),                                          &
         bmatrix_snow_soil( begc:endc, 1:, -1: ))

    end associate

  end subroutine SetMatrix_Snow_SoilUrban

  !-----------------------------------------------------------------------
  subroutine SetMatrix_Snow_SoilUrbanNonRoad(bounds, num_nolakec, filter_nolakec, nband, &
       tk, fact, bmatrix_snow_soil)
    !
    ! !DESCRIPTION:
    ! Setup the matrix entries correspodning to snow-soil interaction for
    ! urban sunwall/shadewall/roof columns
    !
    ! !USES:
    use clmtype
    use clm_varcon     , only : cnfac
    use column_varcon  , only : icol_roof, icol_sunwall, icol_shadewall
    use clm_varpar     , only : nlevsno, nlevgrnd, nlevurb
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds                               ! bounds
    integer , intent(in)  :: num_nolakec                                  ! number of column non-lake points in column filter
    integer , intent(in)  :: filter_nolakec(:)                            ! column filter for non-lake points
    integer , intent(in)  :: nband                                        ! number of bands of the tridigonal matrix
    real(r8), intent(in)  :: tk(bounds%begc: ,-nlevsno+1: )               ! thermal conductivity [W/(m K)]
    real(r8), intent(in)  :: fact( bounds%begc: , -nlevsno+1: )           ! used in computing tridiagonal matrix [col, lev]
    real(r8), intent(inout) :: bmatrix_snow_soil(bounds%begc: , 1:,-1: )  ! matrix enteries
    !
    ! !LOCAL VARIABLES:
    integer  :: j,c,l                                                     ! indices
    integer  :: fc                                                        ! lake filtered column indices
    real(r8) :: dzm                                                       ! used in computing tridiagonal matrix
    real(r8) :: dzp                                                       ! used in computing tridiagonal matrix
    !-----------------------------------------------------------------------

    ! Enforce expected array sizes
    SHR_ASSERT_ALL((ubound(tk)                  == (/bounds%endc, nlevgrnd/)),  errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(fact)                == (/bounds%endc, nlevgrnd/)),  errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(bmatrix_snow_soil)   == (/bounds%endc, nband, -1/)), errMsg(__FILE__, __LINE__))

   associate(&
   snl                       =>    cps%snl                       , & ! Input:  [integer (:)]  number of snow layers                    
   ctype                     =>    col%itype                     , & ! Input:  [integer (:)]  column type
   urbpoi                    =>    lun%urbpoi                    , & ! Input:  [logical (:)]  true => landunit is an urban point
   clandunit                 =>    col%landunit                  , & ! Input:  [integer (:)]  column's landunit
   z                         =>    cps%z                         , & ! Input:  [real(r8) (:,:)]  layer thickness (m)
   zi                        =>    cps%zi                        , & ! Input:  [real(r8) (:,:)]  interface level below a "z" level (m)
   begc                      =>    bounds%begc                   , & ! Input:  [integer ] beginning column index
   endc                      =>    bounds%endc                     & ! Input:  [integer ] ending column index
   )

    !
    ! urban non-road columns ---------------------------------------------------------
    !
    do j = 0,0
       do fc = 1,num_nolakec
          c = filter_nolakec(fc)
          l = clandunit(c)
          if (urbpoi(l)) then
             if ((ctype(c) == icol_sunwall .or. ctype(c) == icol_shadewall &
                 .or. ctype(c) == icol_roof)) then
                if (j >= snl(c)+1) then
                   if (j == snl(c)+1) then
                      dzp     = z(c,j+1)-z(c,j)
                      bmatrix_snow_soil(c,1,j-1) =  -(1._r8-cnfac)*fact(c,j)*tk(c,j)/dzp
                   else if (j <= nlevurb-1) then
                      dzm     = (z(c,j)-z(c,j-1))
                      dzp     = (z(c,j+1)-z(c,j))
                      bmatrix_snow_soil(c,1,j-1) =   - (1._r8-cnfac)*fact(c,j)* tk(c,j)/dzp
                   end if
                end if
             end if
          end if
       enddo
    end do

    end associate

  end subroutine SetMatrix_Snow_SoilUrbanNonRoad

  !-----------------------------------------------------------------------
  subroutine SetMatrix_Snow_SoilUrbanRoad(bounds, num_nolakec, filter_nolakec, nband, &
       tk, fact, bmatrix_snow_soil)
    !
    ! !DESCRIPTION:
    ! Setup the matrix entries correspodning to snow-soil interaction for
    ! urban road (impervious + pervious) columns
    !
    ! !USES:
    use clmtype
    use clm_varcon     , only : cnfac
    use column_varcon  , only : icol_road_perv, icol_road_imperv
    use clm_varpar     , only : nlevsno, nlevgrnd, nlevurb
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds                               ! bounds
    integer , intent(in)  :: num_nolakec                                  ! number of column non-lake points in column filter
    integer , intent(in)  :: filter_nolakec(:)                            ! column filter for non-lake points
    integer , intent(in)  :: nband                                        ! number of bands of the tridigonal matrix
    real(r8), intent(in)  :: tk(bounds%begc: ,-nlevsno+1: )               ! thermal conductivity [W/(m K)]
    real(r8), intent(in)  :: fact( bounds%begc: , -nlevsno+1: )           ! used in computing tridiagonal matrix [col, lev]
    real(r8), intent(inout) :: bmatrix_snow_soil(bounds%begc: , 1:,-1: )  ! matrix enteries
    !
    ! !LOCAL VARIABLES:
    integer  :: j,c,l                                                     ! indices
    integer  :: fc                                                        ! lake filtered column indices
    real(r8) :: dzm                                                       ! used in computing tridiagonal matrix
    real(r8) :: dzp                                                       ! used in computing tridiagonal matrix
    !-----------------------------------------------------------------------

    ! Enforce expected array sizes
    SHR_ASSERT_ALL((ubound(tk)                  == (/bounds%endc, nlevgrnd/)),  errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(fact)                == (/bounds%endc, nlevgrnd/)),  errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(bmatrix_snow_soil)   == (/bounds%endc, nband, -1/)), errMsg(__FILE__, __LINE__))

   associate(& 
   snl                       =>    cps%snl                       , & ! Input:  [integer (:)]  number of snow layers                    
   ctype                     =>    col%itype                     , & ! Input:  [integer (:)]  column type
   urbpoi                    =>    lun%urbpoi                    , & ! Input:  [logical (:)]  true => landunit is an urban point
   clandunit                 =>    col%landunit                  , & ! Input:  [integer (:)]  column's landunit
   z                         =>    cps%z                         , & ! Input:  [real(r8) (:,:)]  layer thickness (m)
   zi                        =>    cps%zi                        , & ! Input:  [real(r8) (:,:)]  interface level below a "z" level (m)
   begc                      =>    bounds%begc                   , & ! Input:  [integer ] beginning column index
   endc                      =>    bounds%endc                     & ! Input:  [integer ] ending column index
   )

    !
    ! urban road columns -------------------------------------------------------------
    !
    do j = 0,0
       do fc = 1,num_nolakec
          c = filter_nolakec(fc)
          l = clandunit(c)
          if (urbpoi(l)) then
             if (ctype(c) == icol_road_imperv .or. ctype(c) == icol_road_perv) then
                if (j >= snl(c)+1) then
                   if (j == snl(c)+1) then
                      dzp     = z(c,j+1)-z(c,j)
                      bmatrix_snow_soil(c,1,j-1) =  -(1._r8-cnfac)*fact(c,j)*tk(c,j)/dzp
                   else if (j <= nlevgrnd-1) then
                     dzm     = (z(c,j)-z(c,j-1))
                     dzp     = (z(c,j+1)-z(c,j))
                    bmatrix_snow_soil(c,1,j-1) =   - (1._r8-cnfac)*fact(c,j)* tk(c,j)/dzp
                   end if
                end if
             end if
          end if
       enddo
    end do

    end associate

  end subroutine SetMatrix_Snow_SoilUrbanRoad

  !-----------------------------------------------------------------------
  subroutine SetMatrix_Snow_SoilNonUrban(bounds, num_nolakec, filter_nolakec, nband, &
       tk, fact, bmatrix_snow_soil)
    !
    ! !DESCRIPTION:
    ! Setup the matrix entries correspodning to snow-soil interaction for
    ! non-urban columns
    !
    ! !USES:
    use clmtype
    use clm_varcon     , only : cnfac
    use column_varcon  , only : icol_roof, icol_sunwall, icol_shadewall
    use clm_varpar     , only : nlevsno, nlevgrnd, nlevurb
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds                               ! bounds
    integer , intent(in)  :: num_nolakec                                  ! number of column non-lake points in column filter
    integer , intent(in)  :: filter_nolakec(:)                            ! column filter for non-lake points
    integer , intent(in)  :: nband                                        ! number of bands of the tridigonal matrix
    real(r8), intent(in)  :: tk(bounds%begc: ,-nlevsno+1: )               ! thermal conductivity [W/(m K)]
    real(r8), intent(in)  :: fact( bounds%begc: , -nlevsno+1: )           ! used in computing tridiagonal matrix [col, lev]
    real(r8), intent(inout) :: bmatrix_snow_soil(bounds%begc: , 1:,-1: )  ! matrix enteries
    !
    ! !LOCAL VARIABLES:
    integer  :: j,c,l                                                     ! indices
    integer  :: fc                                                        ! lake filtered column indices
    real(r8) :: dzm                                                       ! used in computing tridiagonal matrix
    real(r8) :: dzp                                                       ! used in computing tridiagonal matrix
    !-----------------------------------------------------------------------

    ! Enforce expected array sizes
    SHR_ASSERT_ALL((ubound(tk)                  == (/bounds%endc, nlevgrnd/)),  errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(fact)                == (/bounds%endc, nlevgrnd/)),  errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(bmatrix_snow_soil)   == (/bounds%endc, nband, -1/)), errMsg(__FILE__, __LINE__))

   associate(&
   snl                       =>    cps%snl                       , & ! Input:  [integer (:)]  number of snow layers                    
   ctype                     =>    col%itype                     , & ! Input:  [integer (:)]  column type
   urbpoi                    =>    lun%urbpoi                    , & ! Input:  [logical (:)]  true => landunit is an urban point
   clandunit                 =>    col%landunit                  , & ! Input:  [integer (:)]  column's landunit
   z                         =>    cps%z                         , & ! Input:  [real(r8) (:,:)]  layer thickness (m)
   zi                        =>    cps%zi                        , & ! Input:  [real(r8) (:,:)]  interface level below a "z" level (m)
   begc                      =>    bounds%begc                   , & ! Input:  [integer ] beginning column index
   endc                      =>    bounds%endc                     & ! Input:  [integer ] ending column index
   )

    !
    ! non-urban columns --------------------------------------------------------------
    !
    do j = 0,0
       do fc = 1,num_nolakec
          c = filter_nolakec(fc)
          l = clandunit(c)
          if (.not. urbpoi(l)) then
             if (j >= snl(c)+1) then
                if (j == snl(c)+1) then
                   dzp     = z(c,j+1)-z(c,j)
                   bmatrix_snow_soil(c,1,j-1) =  -(1._r8-cnfac)*fact(c,j)*tk(c,j)/dzp
                else if (j <= nlevgrnd-1) then
                   dzm     = (z(c,j)-z(c,j-1))
                   dzp     = (z(c,j+1)-z(c,j))
                   bmatrix_snow_soil(c,1,j-1) =   - (1._r8-cnfac)*fact(c,j)* tk(c,j)/dzp
                end if
             end if
          end if
       enddo
    end do

    end associate

  end subroutine SetMatrix_Snow_SoilNonUrban

  !-----------------------------------------------------------------------
  subroutine SetMatrix_Soil(bounds, num_nolakec, filter_nolakec, nband, &
       dhsdT, tk, tk_h2osfc, dz_h2osfc, fact, bmatrix_soil)
    !
    ! !DESCRIPTION:
    ! Setup the matrix entries correspodning to internal soil layers.
    !
    ! !USES:
    use clmtype
    use clm_varcon     , only : cnfac
    use column_varcon  , only : icol_roof, icol_sunwall, icol_shadewall
    use clm_varpar     , only : nlevsno, nlevgrnd, nlevurb
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds                        ! bounds
    integer , intent(in)  :: num_nolakec                           ! number of column non-lake points in column filter
    integer , intent(in)  :: filter_nolakec(:)                     ! column filter for non-lake points
    integer , intent(in)  :: nband                                 ! number of bands of the tridigonal matrix
    real(r8), intent(in)  :: dhsdT(bounds%begc: )                  ! temperature derivative of "hs" [col]
    real(r8), intent(in)  :: tk(bounds%begc: ,-nlevsno+1: )        ! thermal conductivity [W/(m K)]
    real(r8), intent(in)  :: tk_h2osfc(bounds%begc: )              ! thermal conductivity [W/(m K)]
    real(r8), intent(in)  :: dz_h2osfc(bounds%begc: )              ! Thickness of standing water [m]
    real(r8), intent(in)  :: fact( bounds%begc: , -nlevsno+1: )    ! used in computing tridiagonal matrix [col, lev]
    real(r8), intent(out) :: bmatrix_soil(bounds%begc: , 1:, 1: )  ! matrix enteries
    !
    ! !LOCAL VARIABLES:
    integer  :: j,c,l                                              ! indices
    integer  :: fc                                                 ! lake filtered column indices
    real(r8) :: dzm                                                ! used in computing tridiagonal matrix
    real(r8) :: dzp                                                ! used in computing tridiagonal matrix
    !-----------------------------------------------------------------------

    ! Enforce expected array sizes
    SHR_ASSERT_ALL((ubound(dhsdT)          == (/bounds%endc/)),                  errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(tk)             == (/bounds%endc, nlevgrnd/)),        errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(tk_h2osfc)      == (/bounds%endc/)),                  errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(dz_h2osfc)      == (/bounds%endc/)),                  errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(fact)           == (/bounds%endc, nlevgrnd/)),        errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(bmatrix_soil)   == (/bounds%endc, nband, nlevgrnd/)), errMsg(__FILE__, __LINE__))

   associate(&
   snl                       =>    cps%snl                       , & ! Input:  [integer (:)]  number of snow layers                    
   frac_h2osfc               =>    cps%frac_h2osfc               , & ! Input:  [real(r8) (:)]  fraction of ground covered by surface water (0 to 1)
   frac_sno_eff              =>    cps%frac_sno_eff              , & ! Input:  [real(r8) (:)]  eff. fraction of ground covered by snow (0 to 1)
   ctype                     =>    col%itype                     , & ! Input:  [integer (:)]  column type
   urbpoi                    =>    lun%urbpoi                    , & ! Input:  [logical (:)]  true => landunit is an urban point
   clandunit                 =>    col%landunit                  , & ! Input:  [integer (:)]  column's landunit
   zi                        =>    cps%zi                        , & ! Input:  [real(r8) (:,:)]  interface level below a "z" level (m)
   z                         =>    cps%z                         , & ! Input:  [real(r8) (:,:)]  layer thickness (m)
   begc                      =>    bounds%begc                   , & ! Input:  [integer ] beginning column index
   endc                      =>    bounds%endc                     & ! Input:  [integer ] ending column index
   )

    ! Initialize
    bmatrix_soil(begc:endc, :, :) = 0.0_r8

    call SetMatrix_SoilUrban(bounds, num_nolakec, filter_nolakec, nband, &
         dhsdT( begc:endc ),                                             &
         tk( begc:endc, -nlevsno+1: ),                                   &
         tk_h2osfc( begc:endc ),                                         &
         dz_h2osfc( begc:endc ),                                         &
         fact( begc:endc, -nlevsno+1: ),                                 &
         bmatrix_soil( begc:endc, 1:, 1: ))

    call SetMatrix_SoilNonUrban(bounds, num_nolakec, filter_nolakec, nband, &
         dhsdT( begc:endc ),                                                &
         tk( begc:endc, -nlevsno+1: ),                                      &
         tk_h2osfc( begc:endc ),                                            &
         dz_h2osfc( begc:endc ),                                            &
         fact( begc:endc, -nlevsno+1: ),                                    &
         bmatrix_soil( begc:endc, 1:, 1: ))

    ! the solution will be organized as (snow:h2osfc:soil) to minimize
    !     bandwidth; this requires a 5-element band instead of 3
    do fc = 1,num_nolakec
       c = filter_nolakec(fc)

       ! surface water layer has two coefficients
       dzm=(0.5*dz_h2osfc(c)+z(c,1))

       ! diagonal element correction for presence of h2osfc
       if ( frac_h2osfc(c) /= 0.0_r8 ) then
          bmatrix_soil(c,3,1)=bmatrix_soil(c,3,1)+ frac_h2osfc(c) &
               *((1._r8-cnfac)*fact(c,1)*tk_h2osfc(c)/dzm + fact(c,1)*dhsdT(c))
       end if

    enddo

    end associate

  end subroutine SetMatrix_Soil

  !-----------------------------------------------------------------------
  subroutine SetMatrix_SoilUrban(bounds, num_nolakec, filter_nolakec, nband, &
       dhsdT, tk, tk_h2osfc, dz_h2osfc, fact, bmatrix_soil)
    !
    ! !DESCRIPTION:
    ! Setup the matrix entries correspodning to internal soil layers for
    ! urban columns
    !
    ! !USES:
    use clmtype
    use clm_varcon     , only : cnfac
    use column_varcon  , only : icol_roof, icol_sunwall, icol_shadewall
    use clm_varpar     , only : nlevsno, nlevgrnd, nlevurb
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds                          ! bounds
    integer , intent(in)  :: num_nolakec                             ! number of column non-lake points in column filter
    integer , intent(in)  :: filter_nolakec(:)                       ! column filter for non-lake points
    integer , intent(in)  :: nband                                   ! number of bands of the tridigonal matrix
    real(r8), intent(in)  :: dhsdT(bounds%begc: )                    ! temperature derivative of "hs" [col]
    real(r8), intent(in)  :: tk(bounds%begc: ,-nlevsno+1: )          ! thermal conductivity [W/(m K)]
    real(r8), intent(in)  :: tk_h2osfc(bounds%begc: )                ! thermal conductivity [W/(m K)]
    real(r8), intent(in)  :: dz_h2osfc(bounds%begc: )                ! Thickness of standing water [m]
    real(r8), intent(in)  :: fact( bounds%begc: , -nlevsno+1: )      ! used in computing tridiagonal matrix [col, lev]
    real(r8), intent(inout) :: bmatrix_soil(bounds%begc: , 1:, 1: )  ! matrix enteries
    !-----------------------------------------------------------------------

    ! Enforce expected array sizes
    SHR_ASSERT_ALL((ubound(dhsdT)          == (/bounds%endc/)),                  errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(tk)             == (/bounds%endc, nlevgrnd/)),        errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(tk_h2osfc)      == (/bounds%endc/)),                  errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(dz_h2osfc)      == (/bounds%endc/)),                  errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(fact)           == (/bounds%endc, nlevgrnd/)),        errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(bmatrix_soil)   == (/bounds%endc, nband, nlevgrnd/)), errMsg(__FILE__, __LINE__))

   associate(& 
   begc                      =>    bounds%begc                   , & ! Input:  [integer ] beginning column index
   endc                      =>    bounds%endc                     & ! Input:  [integer ] ending column index
   )

    call SetMatrix_SoilUrbanNonRoad(bounds, num_nolakec, filter_nolakec, nband, &
         dhsdT( begc:endc ),                                                    &
         tk( begc:endc, -nlevsno+1: ),                                          &
         tk_h2osfc( begc:endc ),                                                &
         dz_h2osfc( begc:endc ),                                                &
         fact( begc:endc, -nlevsno+1: ),                                        &
         bmatrix_soil( begc:endc, 1:, 1: ))

    call SetMatrix_SoilUrbanRoad(bounds, num_nolakec, filter_nolakec, nband, &
         dhsdT( begc:endc ),                                                 &
         tk( begc:endc, -nlevsno+1: ),                                       &
         tk_h2osfc( begc:endc ),                                             &
         dz_h2osfc( begc:endc ),                                             &
         fact( begc:endc, -nlevsno+1: ),                                     &
         bmatrix_soil( begc:endc, 1:, 1: ))

    end associate

  end subroutine SetMatrix_SoilUrban

  !-----------------------------------------------------------------------
  subroutine SetMatrix_SoilUrbanNonRoad(bounds, num_nolakec, filter_nolakec, nband, &
       dhsdT, tk, tk_h2osfc, dz_h2osfc, fact, bmatrix_soil)
    !
    ! !DESCRIPTION:
    ! Setup the matrix entries correspodning to internal soil layers for
    ! urban sunwall/shadewall/roof columns
    !
    ! !USES:
    use clmtype
    use clm_varcon     , only : cnfac
    use column_varcon  , only : icol_roof, icol_sunwall, icol_shadewall
    use clm_varpar     , only : nlevsno, nlevgrnd, nlevurb
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds                          ! bounds
    integer , intent(in)  :: num_nolakec                             ! number of column non-lake points in column filter
    integer , intent(in)  :: filter_nolakec(:)                       ! column filter for non-lake points
    integer , intent(in)  :: nband                                   ! number of bands of the tridigonal matrix
    real(r8), intent(in)  :: dhsdT(bounds%begc: )                    ! temperature derivative of "hs" [col]
    real(r8), intent(in)  :: tk(bounds%begc: ,-nlevsno+1: )          ! thermal conductivity [W/(m K)]
    real(r8), intent(in)  :: tk_h2osfc(bounds%begc: )                ! thermal conductivity [W/(m K)]
    real(r8), intent(in)  :: dz_h2osfc(bounds%begc: )                ! Thickness of standing water [m]
    real(r8), intent(in)  :: fact( bounds%begc: , -nlevsno+1: )      ! used in computing tridiagonal matrix [col, lev]
    real(r8), intent(inout) :: bmatrix_soil(bounds%begc: , 1:, 1: )  ! matrix enteries
    !
    ! !LOCAL VARIABLES:
    integer  :: j,c,l                                                ! indices
    integer  :: fc                                                   ! lake filtered column indices
    real(r8) :: dzm                                                  ! used in computing tridiagonal matrix
    real(r8) :: dzp                                                  ! used in computing tridiagonal matrix
    !-----------------------------------------------------------------------

    ! Enforce expected array sizes
    SHR_ASSERT_ALL((ubound(dhsdT)          == (/bounds%endc/)),                  errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(tk)             == (/bounds%endc, nlevgrnd/)),        errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(tk_h2osfc)      == (/bounds%endc/)),                  errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(dz_h2osfc)      == (/bounds%endc/)),                  errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(fact)           == (/bounds%endc, nlevgrnd/)),        errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(bmatrix_soil)   == (/bounds%endc, nband, nlevgrnd/)), errMsg(__FILE__, __LINE__))

   associate(&
   snl                       =>    cps%snl                       , & ! Input:  [integer (:)]  number of snow layers                    
   frac_h2osfc               =>    cps%frac_h2osfc               , & ! Input:  [real(r8) (:)]  fraction of ground covered by surface water (0 to 1)
   frac_sno_eff              =>    cps%frac_sno_eff              , & ! Input:  [real(r8) (:)]  eff. fraction of ground covered by snow (0 to 1)
   ctype                     =>    col%itype                     , & ! Input:  [integer (:)]  column type
   urbpoi                    =>    lun%urbpoi                    , & ! Input:  [logical (:)]  true => landunit is an urban point
   clandunit                 =>    col%landunit                  , & ! Input:  [integer (:)]  column's landunit
   zi                        =>    cps%zi                        , & ! Input:  [real(r8) (:,:)]  interface level below a "z" level (m)
   z                         =>    cps%z                         , & ! Input:  [real(r8) (:,:)]  layer thickness (m)
   begc                      =>    bounds%begc                   , & ! Input:  [integer ] beginning column index
   endc                      =>    bounds%endc                     & ! Input:  [integer ] ending column index
   )

    !
    ! urban non-road columns ---------------------------------------------------------
    !
    do j = 1,nlevurb
       do fc = 1,num_nolakec
          c = filter_nolakec(fc)
          l = clandunit(c)
          if (urbpoi(l)) then
             if ((ctype(c) == icol_sunwall .or. ctype(c) == icol_shadewall &
                 .or. ctype(c) == icol_roof)) then
                if (j >= snl(c)+1) then
                   if (j == snl(c)+1) then
                      dzp     = z(c,j+1)-z(c,j)
                      if (j /= 1) then
                         bmatrix_soil(c,4,j) = 0._r8
                       end if
                      bmatrix_soil(c,3,j) = 1+(1._r8-cnfac)*fact(c,j)*tk(c,j)/dzp-fact(c,j)*dhsdT(c)
                      bmatrix_soil(c,2,j) =  -(1._r8-cnfac)*fact(c,j)*tk(c,j)/dzp
                  else if (j <= nlevurb-1) then
                      dzm     = (z(c,j)-z(c,j-1))
                      dzp     = (z(c,j+1)-z(c,j))
                      if (j /= 1) then
                         bmatrix_soil(c,4,j) =   - (1._r8-cnfac)*fact(c,j)* tk(c,j-1)/dzm
                      end if
                      bmatrix_soil(c,3,j) = 1._r8+ (1._r8-cnfac)*fact(c,j)*(tk(c,j)/dzp + tk(c,j-1)/dzm)
                      bmatrix_soil(c,2,j) =   - (1._r8-cnfac)*fact(c,j)* tk(c,j)/dzp
                   else if (j == nlevurb) then
                      ! For urban sunwall, shadewall, and roof columns, there is a non-zero heat flux across
                      ! the bottom "soil" layer and the equations are derived assuming a prescribed internal
                      ! building temperature. (See Oleson urban notes of 6/18/03).
                      dzm     = ( z(c,j)-z(c,j-1))
                      dzp     = (zi(c,j)-z(c,j))
                      bmatrix_soil(c,4,j) =   - (1._r8-cnfac)*fact(c,j)*(tk(c,j-1)/dzm)
                      bmatrix_soil(c,3,j) = 1._r8+ (1._r8-cnfac)*fact(c,j)*(tk(c,j-1)/dzm + tk(c,j)/dzp)
                      bmatrix_soil(c,2,j) = 0._r8
                   end if
                end if
             end if
          end if
       enddo
    end do

    end associate

  end subroutine SetMatrix_SoilUrbanNonRoad

  !-----------------------------------------------------------------------
  subroutine SetMatrix_SoilUrbanRoad(bounds, num_nolakec, filter_nolakec, nband, &
       dhsdT, tk, tk_h2osfc, dz_h2osfc, fact, bmatrix_soil)
    !
    ! !DESCRIPTION:
    ! Setup the matrix entries correspodning to internal soil layers for
    ! urban road (impervious + pervious) columns
    !
    ! !USES:
    use clmtype
    use clm_varcon     , only : cnfac
    use column_varcon  , only : icol_road_perv, icol_road_imperv
    use clm_varpar     , only : nlevsno, nlevgrnd, nlevurb
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds                          ! bounds
    integer , intent(in)  :: num_nolakec                             ! number of column non-lake points in column filter
    integer , intent(in)  :: filter_nolakec(:)                       ! column filter for non-lake points
    integer , intent(in)  :: nband                                   ! number of bands of the tridigonal matrix
    real(r8), intent(in)  :: dhsdT(bounds%begc: )                    ! temperature derivative of "hs" [col]
    real(r8), intent(in)  :: tk(bounds%begc: ,-nlevsno+1: )          ! thermal conductivity [W/(m K)]
    real(r8), intent(in)  :: tk_h2osfc(bounds%begc: )                ! thermal conductivity [W/(m K)]
    real(r8), intent(in)  :: dz_h2osfc(bounds%begc: )                ! Thickness of standing water [m]
    real(r8), intent(in)  :: fact( bounds%begc: , -nlevsno+1: )      ! used in computing tridiagonal matrix [col, lev]
    real(r8), intent(inout) :: bmatrix_soil(bounds%begc: , 1:, 1: )  ! matrix enteries
    !
    ! !LOCAL VARIABLES:
    integer  :: j,c,l                                                ! indices
    integer  :: fc                                                   ! lake filtered column indices
    real(r8) :: dzm                                                  ! used in computing tridiagonal matrix
    real(r8) :: dzp                                                  ! used in computing tridiagonal matrix
    !-----------------------------------------------------------------------

    ! Enforce expected array sizes
    SHR_ASSERT_ALL((ubound(dhsdT)          == (/bounds%endc/)),                  errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(tk)             == (/bounds%endc, nlevgrnd/)),        errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(tk_h2osfc)      == (/bounds%endc/)),                  errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(dz_h2osfc)      == (/bounds%endc/)),                  errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(fact)           == (/bounds%endc, nlevgrnd/)),        errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(bmatrix_soil)   == (/bounds%endc, nband, nlevgrnd/)), errMsg(__FILE__, __LINE__))

   associate(&
   snl                       =>    cps%snl                       , & ! Input:  [integer (:)]  number of snow layers                    
   frac_h2osfc               =>    cps%frac_h2osfc               , & ! Input:  [real(r8) (:)]  fraction of ground covered by surface water (0 to 1)
   frac_sno_eff              =>    cps%frac_sno_eff              , & ! Input:  [real(r8) (:)]  eff. fraction of ground covered by snow (0 to 1)
   ctype                     =>    col%itype                     , & ! Input:  [integer (:)]  column type
   urbpoi                    =>    lun%urbpoi                    , & ! Input:  [logical (:)]  true => landunit is an urban point
   clandunit                 =>    col%landunit                  , & ! Input:  [integer (:)]  column's landunit
   zi                        =>    cps%zi                        , & ! Input:  [real(r8) (:,:)]  interface level below a "z" level (m)
   z                         =>    cps%z                         , & ! Input:  [real(r8) (:,:)]  layer thickness (m)
   begc                      =>    bounds%begc                   , & ! Input:  [integer ] beginning column index
   endc                      =>    bounds%endc                     & ! Input:  [integer ] ending column index
   )

    !
    ! urban road columns -------------------------------------------------------------
    !
    do j = 1,nlevgrnd
       do fc = 1,num_nolakec
          c = filter_nolakec(fc)
          l = clandunit(c)
          if (urbpoi(l)) then
             if (ctype(c) == icol_road_imperv .or. ctype(c) == icol_road_perv) then
                if (j >= snl(c)+1) then
                  if (j == snl(c)+1) then
                     dzp     = z(c,j+1)-z(c,j)
                     if (j /= 1) then
                        bmatrix_soil(c,4,j) = 0._r8
                     end if
                     bmatrix_soil(c,3,j) = 1+(1._r8-cnfac)*fact(c,j)*tk(c,j)/dzp-fact(c,j)*dhsdT(c)
                     bmatrix_soil(c,2,j) =  -(1._r8-cnfac)*fact(c,j)*tk(c,j)/dzp
                  else if (j == 1) then
                     ! this is the snow/soil interface layer
                     dzm     = (z(c,j)-z(c,j-1))
                     dzp     = (z(c,j+1)-z(c,j))
                     if (j /= 1) then
                        bmatrix_soil(c,4,j) =   - frac_sno_eff(c) * (1._r8-cnfac) * fact(c,j) &
                            * tk(c,j-1)/dzm
                     end if
                     bmatrix_soil(c,3,j) = 1._r8 + (1._r8-cnfac)*fact(c,j)*(tk(c,j)/dzp &
                          + frac_sno_eff(c) * tk(c,j-1)/dzm) &
                          - (1._r8 - frac_sno_eff(c))*fact(c,j)*dhsdT(c)
                     bmatrix_soil(c,2,j) = - (1._r8-cnfac)*fact(c,j)*tk(c,j)/dzp
                  else if (j <= nlevgrnd-1) then
                     dzm     = (z(c,j)-z(c,j-1))
                     dzp     = (z(c,j+1)-z(c,j))
                     bmatrix_soil(c,4,j) =   - (1._r8-cnfac)*fact(c,j)* tk(c,j-1)/dzm
                     bmatrix_soil(c,3,j) = 1._r8+ (1._r8-cnfac)*fact(c,j)*(tk(c,j)/dzp + tk(c,j-1)/dzm)
                     bmatrix_soil(c,2,j) =   - (1._r8-cnfac)*fact(c,j)* tk(c,j)/dzp
                  else if (j == nlevgrnd) then
                     dzm     = (z(c,j)-z(c,j-1))
                     bmatrix_soil(c,4,j) =   - (1._r8-cnfac)*fact(c,j)*tk(c,j-1)/dzm
                     bmatrix_soil(c,3,j) = 1._r8+ (1._r8-cnfac)*fact(c,j)*tk(c,j-1)/dzm
                     bmatrix_soil(c,2,j) = 0._r8
                  end if
               end if
            end if
          end if
       enddo
    end do

    end associate

  end subroutine SetMatrix_SoilUrbanRoad

  !-----------------------------------------------------------------------
  subroutine SetMatrix_SoilNonUrban(bounds, num_nolakec, filter_nolakec, nband, &
       dhsdT, tk, tk_h2osfc, dz_h2osfc, fact, bmatrix_soil)
    !
    ! !DESCRIPTION:
    ! Setup the matrix entries correspodning to internal soil layers.
    !
    ! !USES:
    use clmtype
    use clm_varcon     , only : cnfac
    use column_varcon  , only : icol_roof, icol_sunwall, icol_shadewall
    use clm_varpar     , only : nlevsno, nlevgrnd, nlevurb
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds                          ! bounds
    integer , intent(in)  :: num_nolakec                             ! number of column non-lake points in column filter
    integer , intent(in)  :: filter_nolakec(:)                       ! column filter for non-lake points
    integer , intent(in)  :: nband                                   ! number of bands of the tridigonal matrix
    real(r8), intent(in)  :: dhsdT(bounds%begc: )                    ! temperature derivative of "hs" [col]
    real(r8), intent(in)  :: tk(bounds%begc: ,-nlevsno+1: )          ! thermal conductivity [W/(m K)]
    real(r8), intent(in)  :: tk_h2osfc(bounds%begc: )                ! thermal conductivity [W/(m K)]
    real(r8), intent(in)  :: dz_h2osfc(bounds%begc: )                ! Thickness of standing water [m]
    real(r8), intent(in)  :: fact( bounds%begc: , -nlevsno+1: )      ! used in computing tridiagonal matrix [col, lev]
    real(r8), intent(inout) :: bmatrix_soil(bounds%begc: , 1:, 1: )  ! matrix enteries
    !
    ! !LOCAL VARIABLES:
    integer  :: j,c,l                                                ! indices
    integer  :: fc                                                   ! lake filtered column indices
    real(r8) :: dzm                                                  ! used in computing tridiagonal matrix
    real(r8) :: dzp                                                  ! used in computing tridiagonal matrix
    !------------------------------------------------------------------------------

    ! Enforce expected array sizes
    SHR_ASSERT_ALL((ubound(dhsdT)          == (/bounds%endc/)),                  errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(tk)             == (/bounds%endc, nlevgrnd/)),        errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(tk_h2osfc)      == (/bounds%endc/)),                  errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(dz_h2osfc)      == (/bounds%endc/)),                  errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(fact)           == (/bounds%endc, nlevgrnd/)),        errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(bmatrix_soil)   == (/bounds%endc, nband, nlevgrnd/)), errMsg(__FILE__, __LINE__))

   associate(&
   snl                       =>    cps%snl                       , & ! Input:  [integer (:)]  number of snow layers                    
   frac_h2osfc               =>    cps%frac_h2osfc               , & ! Input:  [real(r8) (:)]  fraction of ground covered by surface water (0 to 1)
   frac_sno_eff              =>    cps%frac_sno_eff              , & ! Input:  [real(r8) (:)]  eff. fraction of ground covered by snow (0 to 1)
   ctype                     =>    col%itype                     , & ! Input:  [integer (:)]  column type
   urbpoi                    =>    lun%urbpoi                    , & ! Input:  [logical (:)]  true => landunit is an urban point
   clandunit                 =>    col%landunit                  , & ! Input:  [integer (:)]  column's landunit
   zi                        =>    cps%zi                        , & ! Input:  [real(r8) (:,:)]  interface level below a "z" level (m)
   z                         =>    cps%z                         , & ! Input:  [real(r8) (:,:)]  layer thickness (m)
   begc                      =>    bounds%begc                   , & ! Input:  [integer ] beginning column index
   endc                      =>    bounds%endc                     & ! Input:  [integer ] ending column index
   )

    !
    ! non-urban columns --------------------------------------------------------------
    !
    do j = 1,nlevgrnd
       do fc = 1,num_nolakec
          c = filter_nolakec(fc)
          l = clandunit(c)
          if (.not. urbpoi(l)) then
             if (j >= snl(c)+1) then
                if (j == snl(c)+1) then
                   dzp     = z(c,j+1)-z(c,j)
                   if (j /= 1) then
                      bmatrix_soil(c,4,j) = 0._r8
                   end if
                   bmatrix_soil(c,3,j) = 1+(1._r8-cnfac)*fact(c,j)*tk(c,j)/dzp-fact(c,j)*dhsdT(c)
                   bmatrix_soil(c,2,j) =  -(1._r8-cnfac)*fact(c,j)*tk(c,j)/dzp
                else if (j == 1) then
                   ! this is the snow/soil interface layer
                   dzm     = (z(c,j)-z(c,j-1))
                   dzp     = (z(c,j+1)-z(c,j))
                   if (j /= 1) then
                      bmatrix_soil(c,4,j) =   - frac_sno_eff(c) * (1._r8-cnfac) * fact(c,j) &
                          * tk(c,j-1)/dzm
                   end if
                   bmatrix_soil(c,3,j) = 1._r8 + (1._r8-cnfac)*fact(c,j)*(tk(c,j)/dzp &
                        + frac_sno_eff(c) * tk(c,j-1)/dzm) &
                        - (1._r8 - frac_sno_eff(c))*fact(c,j)*dhsdT(c)
                   bmatrix_soil(c,2,j) = - (1._r8-cnfac)*fact(c,j)*tk(c,j)/dzp
                else if (j <= nlevgrnd-1) then
                   dzm     = (z(c,j)-z(c,j-1))
                   dzp     = (z(c,j+1)-z(c,j))
                   bmatrix_soil(c,4,j) =   - (1._r8-cnfac)*fact(c,j)* tk(c,j-1)/dzm
                   bmatrix_soil(c,3,j) = 1._r8+ (1._r8-cnfac)*fact(c,j)*(tk(c,j)/dzp + tk(c,j-1)/dzm)
                   bmatrix_soil(c,2,j) =   - (1._r8-cnfac)*fact(c,j)* tk(c,j)/dzp
                else if (j == nlevgrnd) then
                   dzm     = (z(c,j)-z(c,j-1))
                   bmatrix_soil(c,4,j) =   - (1._r8-cnfac)*fact(c,j)*tk(c,j-1)/dzm
                   bmatrix_soil(c,3,j) = 1._r8+ (1._r8-cnfac)*fact(c,j)*tk(c,j-1)/dzm
                   bmatrix_soil(c,2,j) = 0._r8
                end if
             end if
          end if
       enddo
    end do

    end associate

  end subroutine SetMatrix_SoilNonUrban


  !-----------------------------------------------------------------------
  subroutine SetMatrix_Soil_Snow(bounds, num_nolakec, filter_nolakec, nband, &
       tk, fact, bmatrix_soil_snow)
    !
    ! !DESCRIPTION:
    ! Setup the matrix entries correspodning to soil-snow interaction
    !
    ! !USES:
    use clmtype
    use clm_varcon     , only : cnfac
    use column_varcon  , only : icol_roof, icol_sunwall, icol_shadewall
    use clm_varpar     , only : nlevsno, nlevgrnd, nlevurb
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds                             ! bounds
    integer , intent(in)  :: num_nolakec                                ! number of column non-lake points in column filter
    integer , intent(in)  :: filter_nolakec(:)                          ! column filter for non-lake points
    integer , intent(in)  :: nband                                      ! number of bands of the tridigonal matrix
    real(r8), intent(in)  :: tk(bounds%begc: ,-nlevsno+1: )             ! thermal conductivity [W/(m K)]
    real(r8), intent(in)  :: fact( bounds%begc: , -nlevsno+1: )         ! used in computing tridiagonal matrix [col, lev]
    real(r8), intent(out) :: bmatrix_soil_snow(bounds%begc: , 1: ,1: )  ! matrix enteries
    !------------------------------------------------------------------------------

    ! Enforce expected array sizes
    SHR_ASSERT_ALL((ubound(tk)                  == (/bounds%endc, nlevgrnd/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(fact)                == (/bounds%endc, nlevgrnd/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(bmatrix_soil_snow)   == (/bounds%endc, nband, 1/)), errMsg(__FILE__, __LINE__))

   associate(&
   begc                      =>    bounds%begc                   , & ! Input:  [integer ] beginning column index
   endc                      =>    bounds%endc                     & ! Input:  [integer ] ending column index
   )

    ! Initialize
    bmatrix_soil_snow(begc:endc, :, :) = 0.0_r8

    call SetMatrix_Soil_SnowUrban(bounds, num_nolakec, filter_nolakec, nband, &
         tk( begc:endc, -nlevsno+1: ),                                        &
         fact( begc:endc, -nlevsno+1: ),                                      &
         bmatrix_soil_snow( begc:endc, 1:, 1: ))

    call SetMatrix_Soil_SnowNonUrban(bounds, num_nolakec, filter_nolakec, nband, &
         tk( begc:endc, -nlevsno+1: ),                                           &
         fact( begc:endc, -nlevsno+1: ),                                         &
         bmatrix_soil_snow( begc:endc, 1:, 1: ))

    end associate

  end subroutine SetMatrix_Soil_Snow

  !-----------------------------------------------------------------------
  subroutine SetMatrix_Soil_SnowUrban(bounds, num_nolakec, filter_nolakec, nband, &
       tk, fact, bmatrix_soil_snow)
    !
    ! !DESCRIPTION:
    ! Setup the matrix entries correspodning to soil-snow interaction for
    ! urban columns
    !
    ! !USES:
    use clmtype
    use clm_varcon     , only : cnfac
    use column_varcon  , only : icol_roof, icol_sunwall, icol_shadewall
    use clm_varpar     , only : nlevsno, nlevgrnd, nlevurb
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds                               ! bounds
    integer , intent(in)  :: num_nolakec                                  ! number of column non-lake points in column filter
    integer , intent(in)  :: filter_nolakec(:)                            ! column filter for non-lake points
    integer , intent(in)  :: nband                                        ! number of bands of the tridigonal matrix
    real(r8), intent(in)  :: tk(bounds%begc: ,-nlevsno+1: )               ! thermal conductivity [W/(m K)]
    real(r8), intent(in)  :: fact( bounds%begc: , -nlevsno+1: )           ! used in computing tridiagonal matrix [col, lev]
    real(r8), intent(inout) :: bmatrix_soil_snow(bounds%begc: , 1: ,1: )  ! matrix enteries
    !-----------------------------------------------------------------------

    ! Enforce expected array sizes
    SHR_ASSERT_ALL((ubound(tk)                  == (/bounds%endc, nlevgrnd/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(fact)                == (/bounds%endc, nlevgrnd/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(bmatrix_soil_snow)   == (/bounds%endc, nband, 1/)), errMsg(__FILE__, __LINE__))

   associate(& 
   begc                      =>    bounds%begc                   , & ! Input:  [integer ] beginning column index
   endc                      =>    bounds%endc                     & ! Input:  [integer ] ending column index
   )

    call SetMatrix_Soil_SnowUrbanNonRoad(bounds, num_nolakec, filter_nolakec, nband, &
         tk( begc:endc, -nlevsno+1: ),                                               &
         fact( begc:endc, -nlevsno+1: ),                                             &
         bmatrix_soil_snow( begc:endc, 1:, 1: ))

    call SetMatrix_Soil_SnowUrbanRoad(bounds, num_nolakec, filter_nolakec, nband, &
         tk( begc:endc, -nlevsno+1: ),                                            &
         fact( begc:endc, -nlevsno+1: ),                                          &
         bmatrix_soil_snow( begc:endc, 1:, 1: ))

    end associate

  end subroutine SetMatrix_Soil_SnowUrban


  !-----------------------------------------------------------------------
  subroutine SetMatrix_Soil_SnowUrbanNonRoad(bounds, num_nolakec, filter_nolakec, nband, &
       tk, fact, bmatrix_soil_snow)
    !
    ! !DESCRIPTION:
    ! Setup the matrix entries correspodning to soil-snow interaction for
    ! urban sunwall/shadewall/roof columns
    !
    ! !USES:
    use clmtype
    use clm_varcon     , only : cnfac
    use column_varcon  , only : icol_roof, icol_sunwall, icol_shadewall
    use clm_varpar     , only : nlevsno, nlevgrnd, nlevurb
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds                               ! bounds
    integer , intent(in)  :: num_nolakec                                  ! number of column non-lake points in column filter
    integer , intent(in)  :: filter_nolakec(:)                            ! column filter for non-lake points
    integer , intent(in)  :: nband                                        ! number of bands of the tridigonal matrix
    real(r8), intent(in)  :: tk(bounds%begc: ,-nlevsno+1: )               ! thermal conductivity [W/(m K)]
    real(r8), intent(in)  :: fact( bounds%begc: , -nlevsno+1: )           ! used in computing tridiagonal matrix [col, lev]
    real(r8), intent(inout) :: bmatrix_soil_snow(bounds%begc: , 1: ,1: )  ! matrix enteries
    !
    ! !LOCAL VARIABLES:
    integer  :: j,c,l                                                     ! indices
    integer  :: fc                                                        ! lake filtered column indices
    real(r8) :: dzm                                                       ! used in computing tridiagonal matrix
    real(r8) :: dzp                                                       ! used in computing tridiagonal matrix
    !-----------------------------------------------------------------------

    ! Enforce expected array sizes
    SHR_ASSERT_ALL((ubound(tk)                  == (/bounds%endc, nlevgrnd/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(fact)                == (/bounds%endc, nlevgrnd/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(bmatrix_soil_snow)   == (/bounds%endc, nband, 1/)), errMsg(__FILE__, __LINE__))

   associate(&
   snl                       =>    cps%snl                       , & ! Input:  [integer (:)]  number of snow layers                    
   frac_sno_eff              =>    cps%frac_sno_eff              , & ! Input:  [real(r8) (:)]  eff. fraction of ground covered by snow (0 to 1)
   ctype                     =>    col%itype                     , & ! Input:  [integer (:)]  column type
   urbpoi                    =>    lun%urbpoi                    , & ! Input:  [logical (:)]  true => landunit is an urban point
   clandunit                 =>    col%landunit                  , & ! Input:  [integer (:)]  column's landunit
   z                         =>    cps%z                         , & ! Input:  [real(r8) (:,:)]  layer thickness (m)
   zi                        =>    cps%zi                        , & ! Input:  [real(r8) (:,:)]  interface level below a "z" level (m)
   begc                      =>    bounds%begc                   , & ! Input:  [integer ] beginning column index
   endc                      =>    bounds%endc                     & ! Input:  [integer ] ending column index
   )

    !
    !
    do j = 1,1
       do fc = 1,num_nolakec
          c = filter_nolakec(fc)
          l = clandunit(c)
          if (urbpoi(l)) then
             if ((ctype(c) == icol_sunwall .or. ctype(c) == icol_shadewall &
                 .or. ctype(c) == icol_roof)) then
                if (j >= snl(c)+1) then
                   if (j == snl(c)+1) then
                      dzp     = z(c,j+1)-z(c,j)
                      bmatrix_soil_snow(c,5,j) = 0._r8
                   else if (j <= nlevurb-1) then
                      dzm     = (z(c,j)-z(c,j-1))
                      dzp     = (z(c,j+1)-z(c,j))
                      bmatrix_soil_snow(c,5,j) =   - (1._r8-cnfac)*fact(c,j)* tk(c,j-1)/dzm
                   end if
                end if
             end if
          end if
       enddo
    end do

    end associate

  end subroutine SetMatrix_Soil_SnowUrbanNonRoad

  !-----------------------------------------------------------------------
  subroutine SetMatrix_Soil_SnowUrbanRoad(bounds, num_nolakec, filter_nolakec, nband, &
       tk, fact, bmatrix_soil_snow)
    !
    ! !DESCRIPTION:
    ! Setup the matrix entries correspodning to soil-snow interaction for
    ! urban road (impervious + pervious) columns
    !
    ! !USES:
    use clmtype
    use clm_varcon     , only : cnfac
    use column_varcon  , only : icol_road_imperv, icol_road_perv
    use clm_varpar     , only : nlevsno, nlevgrnd, nlevurb
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds                               ! bounds
    integer , intent(in)  :: num_nolakec                                  ! number of column non-lake points in column filter
    integer , intent(in)  :: filter_nolakec(:)                            ! column filter for non-lake points
    integer , intent(in)  :: nband                                        ! number of bands of the tridigonal matrix
    real(r8), intent(in)  :: tk(bounds%begc: ,-nlevsno+1: )               ! thermal conductivity [W/(m K)]
    real(r8), intent(in)  :: fact( bounds%begc: , -nlevsno+1: )           ! used in computing tridiagonal matrix [col, lev]
    real(r8), intent(inout) :: bmatrix_soil_snow(bounds%begc: , 1: ,1: )  ! matrix enteries
    !
    ! !LOCAL VARIABLES:
    integer  :: j,c,l                                                     ! indices
    integer  :: fc                                                        ! lake filtered column indices
    real(r8) :: dzm                                                       ! used in computing tridiagonal matrix
    real(r8) :: dzp                                                       ! used in computing tridiagonal matrix
    !-----------------------------------------------------------------------

    ! Enforce expected array sizes
    SHR_ASSERT_ALL((ubound(tk)                  == (/bounds%endc, nlevgrnd/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(fact)                == (/bounds%endc, nlevgrnd/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(bmatrix_soil_snow)   == (/bounds%endc, nband, 1/)), errMsg(__FILE__, __LINE__))

   associate(& 
   snl                       =>    cps%snl                       , & ! Input:  [integer (:)]  number of snow layers                    
   frac_sno_eff              =>    cps%frac_sno_eff              , & ! Input:  [real(r8) (:)]  eff. fraction of ground covered by snow (0 to 1)
   ctype                     =>    col%itype                     , & ! Input:  [integer (:)]  column type
   urbpoi                    =>    lun%urbpoi                    , & ! Input:  [logical (:)]  true => landunit is an urban point
   clandunit                 =>    col%landunit                  , & ! Input:  [integer (:)]  column's landunit
   z                         =>    cps%z                         , & ! Input:  [real(r8) (:,:)]  layer thickness (m)
   zi                        =>    cps%zi                        , & ! Input:  [real(r8) (:,:)]  interface level below a "z" level (m)
   begc                      =>    bounds%begc                   , & ! Input:  [integer ] beginning column index
   endc                      =>    bounds%endc                     & ! Input:  [integer ] ending column index
   )

    !
    ! urban road columns -------------------------------------------------------------
    !
    do j = 1,1
       do fc = 1,num_nolakec
          c = filter_nolakec(fc)
          l = clandunit(c)
          if (urbpoi(l)) then
            if (ctype(c) == icol_road_imperv .or. ctype(c) == icol_road_perv) then
               if (j >= snl(c)+1) then
                  if (j == snl(c)+1) then
                     dzp     = z(c,j+1)-z(c,j)
                     bmatrix_soil_snow(c,5,j) = 0._r8
                  else if (j == 1) then
                     ! this is the snow/soil interface layer
                     dzm     = (z(c,j)-z(c,j-1))
                     dzp     = (z(c,j+1)-z(c,j))

                     bmatrix_soil_snow(c,5,j) =   - frac_sno_eff(c) * (1._r8-cnfac) * fact(c,j) &
                          * tk(c,j-1)/dzm
                  end if
               end if
            end if
          end if
       end do
    end do

    end associate

  end subroutine SetMatrix_Soil_SnowUrbanRoad

  !-----------------------------------------------------------------------
  subroutine SetMatrix_Soil_SnowNonUrban(bounds, num_nolakec, filter_nolakec, nband, &
       tk, fact, bmatrix_soil_snow)
    !
    ! !DESCRIPTION:
    ! Setup the matrix entries correspodning to soil-snow interaction for
    ! non urban columns
    !
    ! !USES:
    use clmtype
    use clm_varcon     , only : cnfac
    use column_varcon  , only : icol_roof, icol_sunwall, icol_shadewall
    use clm_varpar     , only : nlevsno, nlevgrnd, nlevurb
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds                               ! bounds
    integer , intent(in)  :: num_nolakec                                  ! number of column non-lake points in column filter
    integer , intent(in)  :: filter_nolakec(:)                            ! column filter for non-lake points
    integer , intent(in)  :: nband                                        ! number of bands of the tridigonal matrix
    real(r8), intent(in)  :: tk(bounds%begc: ,-nlevsno+1: )               ! thermal conductivity [W/(m K)]
    real(r8), intent(in)  :: fact( bounds%begc: , -nlevsno+1: )           ! used in computing tridiagonal matrix [col, lev]
    real(r8), intent(inout) :: bmatrix_soil_snow(bounds%begc: , 1: ,1: )  ! matrix enteries
    !
    ! !LOCAL VARIABLES:
    integer  :: j,c,l                                                     ! indices
    integer  :: fc                                                        ! lake filtered column indices
    real(r8) :: dzm                                                       ! used in computing tridiagonal matrix
    real(r8) :: dzp                                                       ! used in computing tridiagonal matrix
    !-----------------------------------------------------------------------

    ! Enforce expected array sizes
    SHR_ASSERT_ALL((ubound(tk)                  == (/bounds%endc, nlevgrnd/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(fact)                == (/bounds%endc, nlevgrnd/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(bmatrix_soil_snow)   == (/bounds%endc, nband, 1/)), errMsg(__FILE__, __LINE__))

   associate(&
   snl                       =>    cps%snl                       , & ! Input:  [integer (:)]  number of snow layers                    
   frac_sno_eff              =>    cps%frac_sno_eff              , & ! Input:  [real(r8) (:)]  eff. fraction of ground covered by snow (0 to 1)
   ctype                     =>    col%itype                     , & ! Input:  [integer (:)]  column type
   urbpoi                    =>    lun%urbpoi                    , & ! Input:  [logical (:)]  true => landunit is an urban point
   clandunit                 =>    col%landunit                  , & ! Input:  [integer (:)]  column's landunit
   z                         =>    cps%z                         , & ! Input:  [real(r8) (:,:)]  layer thickness (m)
   zi                        =>    cps%zi                        , & ! Input:  [real(r8) (:,:)]  interface level below a "z" level (m)
   begc                      =>    bounds%begc                   , & ! Input:  [integer ] beginning column index
   endc                      =>    bounds%endc                     & ! Input:  [integer ] ending column index
   )

    !
    ! non-urban columns --------------------------------------------------------------
    !
    do j = 1,1
       do fc = 1,num_nolakec
          c = filter_nolakec(fc)
          l = clandunit(c)
          if (.not. urbpoi(l)) then
             if (j >= snl(c)+1) then
                if (j == snl(c)+1) then
                   dzp     = z(c,j+1)-z(c,j)
                   bmatrix_soil_snow(c,5,j) = 0._r8
                else if (j == 1) then
                   ! this is the snow/soil interface layer
                   dzm     = (z(c,j)-z(c,j-1))
                   dzp     = (z(c,j+1)-z(c,j))

                   bmatrix_soil_snow(c,5,j) =  -frac_sno_eff(c) * (1._r8-cnfac) * fact(c,j) &
                        * tk(c,j-1)/dzm
                end if
             end if
          end if
       end do
    end do

    end associate

  end subroutine SetMatrix_Soil_SnowNonUrban

  !-----------------------------------------------------------------------
  subroutine SetMatrix_StandingSurfaceWater(bounds, num_nolakec, filter_nolakec, dtime, nband, &
       dhsdT, tk, tk_h2osfc, fact, c_h2osfc, dz_h2osfc, bmatrix_ssw)
    !
    ! !DESCRIPTION:
    ! Setup the matrix entries correspodning to internal standing water layer
    !
    ! !USES:
    use clmtype
    use clm_varcon     , only : cnfac
    use column_varcon  , only : icol_roof, icol_sunwall, icol_shadewall
    use clm_varpar     , only : nlevsno, nlevgrnd
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds                        ! bounds
    integer , intent(in)  :: num_nolakec                           ! number of column non-lake points in column filter
    integer , intent(in)  :: filter_nolakec(:)                     ! column filter for non-lake points
    real(r8), intent(in)  :: dtime                                 ! land model time step (sec)
    integer , intent(in)  :: nband                                 ! number of bands of the tridigonal matrix
    real(r8), intent(in)  :: dhsdT(bounds%begc: )                  ! temperature derivative of "hs" [col]
    real(r8), intent(in)  :: tk(bounds%begc: ,-nlevsno+1: )        ! thermal conductivity [W/(m K)]
    real(r8), intent(in)  :: tk_h2osfc(bounds%begc: )              ! thermal conductivity [W/(m K)]
    real(r8), intent(in)  :: fact( bounds%begc: , -nlevsno+1: )    ! used in computing tridiagonal matrix [col, lev]
    real(r8), intent(in)  :: c_h2osfc( bounds%begc: )              ! heat capacity of surface water [col]
    real(r8), intent(in)  :: dz_h2osfc(bounds%begc: )              ! Thickness of standing water [m]
    real(r8), intent(out) :: bmatrix_ssw(bounds%begc: , 1:, 0: )   ! matrix enteries
    !
    ! !LOCAL VARIABLES:
    integer  :: c                                                  ! indices
    integer  :: fc                                                 ! lake filtered column indices
    real(r8) :: dzm                                                ! used in computing tridiagonal matrix
    !-----------------------------------------------------------------------

    ! Enforce expected array sizes
    SHR_ASSERT_ALL((ubound(dhsdT)          == (/bounds%endc/)),           errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(tk)             == (/bounds%endc, nlevgrnd/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(tk_h2osfc)      == (/bounds%endc/)),           errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(fact)           == (/bounds%endc, nlevgrnd/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(c_h2osfc)       == (/bounds%endc/)),           errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(dz_h2osfc)      == (/bounds%endc/)),           errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(bmatrix_ssw)   == (/bounds%endc, nband, 0/)), errMsg(__FILE__, __LINE__))

   associate(& 
   z                         =>    cps%z                         , & ! Input:  [real(r8) (:,:)]  layer thickness (m)
   begc                      =>    bounds%begc                   , & ! Input:  [integer ] beginning column index
   endc                      =>    bounds%endc                     & ! Input:  [integer ] ending column index
   )

    ! Initialize
    bmatrix_ssw(begc:endc, :, :) = 0.0_r8

    do fc = 1,num_nolakec
       c = filter_nolakec(fc)

       ! surface water layer has two coefficients
       dzm=(0.5*dz_h2osfc(c)+z(c,1))

       bmatrix_ssw(c,3,0)= 1+(1._r8-cnfac)*(dtime/c_h2osfc(c)) &
            *tk_h2osfc(c)/dzm -(dtime/c_h2osfc(c))*dhsdT(c) !interaction from atm

    enddo

    end associate

  end subroutine SetMatrix_StandingSurfaceWater

  !-----------------------------------------------------------------------
  subroutine SetMatrix_StandingSurfaceWater_Soil(bounds, num_nolakec, filter_nolakec, dtime, nband, &
       tk, tk_h2osfc, fact, c_h2osfc, dz_h2osfc, bmatrix_ssw_soil)
    !
    ! !DESCRIPTION:
    ! Setup the matrix entries correspodning to standing surface water-soil layer interaction
    !
    ! !USES:
    use clmtype
    use clm_varcon     , only : cnfac
    use column_varcon  , only : icol_roof, icol_sunwall, icol_shadewall
    use clm_varpar     , only : nlevsno, nlevgrnd
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds                             ! bounds
    integer , intent(in)  :: num_nolakec                                ! number of column non-lake points in column filter
    integer , intent(in)  :: filter_nolakec(:)                          ! column filter for non-lake points
    real(r8), intent(in)  :: dtime                                      ! land model time step (sec)
    integer , intent(in)  :: nband                                      ! number of bands of the tridigonal matrix
    real(r8), intent(in)  :: tk(bounds%begc: ,-nlevsno+1: )             ! thermal conductivity [W/(m K)]
    real(r8), intent(in)  :: tk_h2osfc(bounds%begc: )                   ! thermal conductivity [W/(m K)]
    real(r8), intent(in)  :: fact( bounds%begc: , -nlevsno+1: )         ! used in computing tridiagonal matrix [col, lev]
    real(r8), intent(in)  :: c_h2osfc( bounds%begc: )                   ! heat capacity of surface water [col]
    real(r8), intent(in)  :: dz_h2osfc(bounds%begc: )                   ! Thickness of standing water [m]
    real(r8), intent(out) :: bmatrix_ssw_soil(bounds%begc: , 1: ,0: )   ! matrix enteries
    !
    ! !LOCAL VARIABLES:
    integer  :: c                                                       ! indices
    integer  :: fc                                                      ! lake filtered column indices
    real(r8) :: dzm                                                     ! used in computing tridiagonal matrix
    !-----------------------------------------------------------------------

    ! Enforce expected array sizes
    SHR_ASSERT_ALL((ubound(tk)                  == (/bounds%endc, nlevgrnd/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(tk_h2osfc)           == (/bounds%endc/)),           errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(fact)                == (/bounds%endc, nlevgrnd/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(c_h2osfc)            == (/bounds%endc/)),           errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(dz_h2osfc)           == (/bounds%endc/)),           errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(bmatrix_ssw_soil)   == (/bounds%endc, nband, 0/)), errMsg(__FILE__, __LINE__))

   associate(& 
   z                         =>    cps%z                         , & ! Input:  [real(r8) (:,:)]  layer thickness (m)
   begc                      =>    bounds%begc                   , & ! Input:  [integer ] beginning column index
   endc                      =>    bounds%endc                     & ! Input:  [integer ] ending column index
   )

    ! Initialize
    bmatrix_ssw_soil(begc:endc, :, :) = 0.0_r8

    do fc = 1,num_nolakec
       c = filter_nolakec(fc)

       ! surface water layer has two coefficients
       dzm=(0.5*dz_h2osfc(c)+z(c,1))

       bmatrix_ssw_soil(c,2,0)= -(1._r8-cnfac)*(dtime/c_h2osfc(c))*tk_h2osfc(c)/dzm !flux to top soil layer

    enddo

  end associate

  end subroutine SetMatrix_StandingSurfaceWater_Soil

  !-----------------------------------------------------------------------
  subroutine SetMatrix_Soil_StandingSurfaceWater(bounds, num_nolakec, filter_nolakec, nband, &
       tk_h2osfc, fact, dz_h2osfc, bmatrix_soil_ssw)
    !
    ! !DESCRIPTION:
    ! Setup the matrix entries correspodning to soil layer-standing surface water interaction
    !
    ! !USES:
    use clmtype
    use clm_varcon     , only : cnfac
    use column_varcon  , only : icol_roof, icol_sunwall, icol_shadewall
    use clm_varpar     , only : nlevsno, nlevgrnd
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds                             ! bounds
    integer , intent(in)  :: num_nolakec                                ! number of column non-lake points in column filter
    integer , intent(in)  :: filter_nolakec(:)                          ! column filter for non-lake points
    integer , intent(in)  :: nband                                      ! number of bands of the tridigonal matrix
    real(r8), intent(in)  :: tk_h2osfc(bounds%begc: )                   ! thermal conductivity [W/(m K)]
    real(r8), intent(in)  :: fact( bounds%begc: , -nlevsno+1: )         ! used in computing tridiagonal matrix [col, lev]
    real(r8), intent(in)  :: dz_h2osfc(bounds%begc: )                   ! Thickness of standing water [m]
    real(r8), intent(out) :: bmatrix_soil_ssw(bounds%begc: , 1:, 1: )   ! matrix enteries
    !
    ! !LOCAL VARIABLES:
    integer  :: c                                                       ! indices
    integer  :: fc                                                      ! lake filtered column indices
    real(r8) :: dzm                                                     ! used in computing tridiagonal matrix
    !-----------------------------------------------------------------------

    ! Enforce expected array sizes
    SHR_ASSERT_ALL((ubound(tk_h2osfc)           == (/bounds%endc/)),           errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(fact)                == (/bounds%endc, nlevgrnd/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(dz_h2osfc)           == (/bounds%endc/)),           errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(bmatrix_soil_ssw)   == (/bounds%endc, nband, 1/)), errMsg(__FILE__, __LINE__))

   associate(& 
   z                         =>    cps%z                         , & ! Input:  [real(r8) (:,:)]  layer thickness (m)
   frac_h2osfc               =>    cps%frac_h2osfc               , & ! Input:  [real(r8) (:)]  fraction of ground covered by surface water (0 to 1)
   begc                      =>    bounds%begc                   , & ! Input:  [integer ] beginning column index
   endc                      =>    bounds%endc                     & ! Input:  [integer ] ending column index
   )

    ! Initialize
    bmatrix_soil_ssw(begc:endc, :, :) = 0.0_r8

    do fc = 1,num_nolakec
       c = filter_nolakec(fc)

       ! surface water layer has two coefficients
       dzm=(0.5*dz_h2osfc(c)+z(c,1))

       ! top soil layer has sub coef shifted to 2nd super diagonal
       if ( frac_h2osfc(c) /= 0.0_r8 )then
          bmatrix_soil_ssw(c,4,1)=  - frac_h2osfc(c) * (1._r8-cnfac) * fact(c,1) &
               * tk_h2osfc(c)/dzm !flux from h2osfc
       end if
    enddo

    end associate

  end subroutine SetMatrix_Soil_StandingSurfaceWater

end module SoilTemperatureMod
