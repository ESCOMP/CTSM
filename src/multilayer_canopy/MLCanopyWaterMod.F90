module MLCanopyWaterMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Update canopy water
  !
  ! !USES:
  use abortutils, only : endrun
  use clm_varctl, only : iulog
  use shr_kind_mod, only : r8 => shr_kind_r8
  !
  ! !PUBLIC TYPES:
  implicit none
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: CanopyInterception     ! Interception and throughfall
  public :: CanopyEvaporation      ! Update canopy intercepted water for evaporation
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine CanopyInterception (num_filter, filter, mlcanopy_inst)
    !
    ! !DESCRIPTION:
    ! Interception and throughfall
    !
    ! !USES:
    use MLclm_varcon, only : dewmx, maximum_leaf_wetted_fraction, interception_fraction, fwet_exponent
    use MLclm_varcon, only : clm45_interception_p1, clm45_interception_p2
    use MLclm_varctl, only : dtime_substep, fpi_type
    use MLCanopyFluxesType, only : mlcanopy_type
    !
    ! !ARGUMENTS:
    integer, intent(in) :: num_filter            ! Number of patches in filter
    integer, intent(in) :: filter(:)             ! Patch filter
    type(mlcanopy_type), intent(inout) :: mlcanopy_inst
    !
    ! !LOCAL VARIABLES:
    integer  :: fp                               ! Filter index
    integer  :: p                                ! Patch index for CLM g/l/c/p hierarchy
    integer  :: ic                               ! Aboveground layer index
    integer  :: n                                ! Number of leaf layers
    real(r8) :: dtime                            ! Model time step (s)
    real(r8) :: fracrain                         ! Fraction of precipitation that is rain
    real(r8) :: fracsnow                         ! Fraction of precipitation that is snow
    real(r8) :: fpi                              ! Fraction of precipitation intercepted
    real(r8) :: qflx_through_rain                ! Rain precipitation direct through canopy (kg H2O/m2/s)
    real(r8) :: qflx_through_snow                ! Snow precipitation direct through canopy (kg H2O/m2/s)
    real(r8) :: qflx_candrip                     ! Flux of water falling off canopy (kg H2O/m2/s)
    real(r8) :: h2ocanmx                         ! Maximum allowed water on canopy layer (kg H2O/m2)
    real(r8) :: xrun                             ! Excess water that exceeds the leaf capacity (kg H2O/m2/s)
    !---------------------------------------------------------------------

    associate ( &
                                                           ! *** Input ***
    qflx_rain    => mlcanopy_inst%qflx_rain_forcing   , &  ! Rainfall (mm H2O/s = kg H2O/m2/s)
    qflx_snow    => mlcanopy_inst%qflx_snow_forcing   , &  ! Snowfall (mm H2O/s = kg H2O/m2/s)
    lai          => mlcanopy_inst%lai_canopy          , &  ! Leaf area index of canopy (m2/m2)
    sai          => mlcanopy_inst%sai_canopy          , &  ! Stem area index of canopy (m2/m2)
    ncan         => mlcanopy_inst%ncan_canopy         , &  ! Number of aboveground layers
    dlai         => mlcanopy_inst%dlai_profile        , &  ! Canopy layer leaf area index (m2/m2)
    dpai         => mlcanopy_inst%dpai_profile        , &  ! Canopy layer plant area index (m2/m2)
                                                           ! *** Input/Output ***
    h2ocan       => mlcanopy_inst%h2ocan_profile      , &  ! Canopy layer intercepted water (kg H2O/m2)
                                                           ! *** Output ***
    qflx_intr    => mlcanopy_inst%qflx_intr_canopy    , &  ! Intercepted precipitation (kg H2O/m2/s)
    qflx_tflrain => mlcanopy_inst%qflx_tflrain_canopy , &  ! Total rain throughfall onto ground (kg H2O/m2/s)
    qflx_tflsnow => mlcanopy_inst%qflx_tflsnow_canopy , &  ! Total snow throughfall onto ground (kg H2O/m2/s)
    fwet         => mlcanopy_inst%fwet_profile        , &  ! Canopy layer fraction of plant area index that is wet
    fdry         => mlcanopy_inst%fdry_profile          &  ! Canopy layer fraction of plant area index that is green and dry
    )

    ! Time step

    dtime = dtime_substep

    do fp = 1, num_filter
       p = filter(fp)

       ! Fraction of precipitation that is rain and snow

       if ((qflx_snow(p) + qflx_rain(p)) > 0._r8) then
          fracrain = qflx_rain(p) / (qflx_snow(p) + qflx_rain(p))
          fracsnow = qflx_snow(p) / (qflx_snow(p) + qflx_rain(p))
       else
          fracrain = 0._r8
          fracsnow = 0._r8
       end if

       ! Fraction of precipitation that is intercepted

       select case (fpi_type)
       case (1) ! CLM4.5
          fpi = clm45_interception_p1 * (1._r8 - exp(clm45_interception_p2*(lai(p) + sai(p))))
       case (2) ! CLM5
          fpi = interception_fraction * tanh(lai(p) + sai(p))
       case default
          call endrun (msg=' ERROR: CanopyInterception: fpi_type not valid')
       end select

       ! Direct throughfall

       qflx_through_rain = qflx_rain(p) * (1._r8 - fpi)
       qflx_through_snow = qflx_snow(p) * (1._r8 - fpi)

       ! Intercepted precipitation

       qflx_intr(p) = (qflx_snow(p) + qflx_rain(p)) * fpi

       ! Find number of layers with lai+sai

       n = 0
       do ic = 1, ncan(p)
          if (dpai(p,ic) > 0._r8) n = n + 1
       end do

       ! Loop through layers for water balance calculation

       qflx_candrip = 0._r8
       do ic = 1, ncan(p)

          if (dpai(p,ic) > 0._r8) then

             ! Maximum external water held in layer

             h2ocanmx = dewmx * dpai(p,ic)

             ! Water storage of intercepted precipitation. Intercepted water
             ! is applied equally to all layers.

             h2ocan(p,ic) = h2ocan(p,ic) + qflx_intr(p) * dtime / float(n)

             ! Excess water that exceeds the maximum capacity. If xrun > 0
             ! then h2ocan is set to h2ocanmx and excess water is added to
             ! throughfall.

             xrun = (h2ocan(p,ic) - h2ocanmx) / dtime
             if (xrun > 0._r8) then
                qflx_candrip = qflx_candrip + xrun
                h2ocan(p,ic) = h2ocanmx
             end if

             ! Wetted fraction of canopy

             fwet(p,ic) = max((h2ocan(p,ic)/h2ocanmx),0._r8)**fwet_exponent
             fwet(p,ic) = min (fwet(p,ic), maximum_leaf_wetted_fraction)

             ! Fraction of canopy that is green and dry 

             fdry(p,ic) = (1._r8 - fwet(p,ic)) * (dlai(p,ic) / dpai(p,ic))

          else

             h2ocan(p,ic) = 0._r8
             fwet(p,ic) = 0._r8
             fdry(p,ic) = 0._r8

          end if

       end do

       ! Total throughfall onto ground

       qflx_tflrain(p) = qflx_through_rain + qflx_candrip * fracrain
       qflx_tflsnow(p) = qflx_through_snow + qflx_candrip * fracsnow

    end do

    end associate
  end subroutine CanopyInterception

  !-----------------------------------------------------------------------
  subroutine CanopyEvaporation (num_filter, filter, mlcanopy_inst)
    !
    ! !DESCRIPTION:
    ! Update canopy intercepted water for evaporation and dew
    !
    ! !USES:
    use MLclm_varcon, only : mmh2o
    use MLclm_varctl, only : dtime_substep
    use MLclm_varpar, only : isun, isha
    use MLCanopyFluxesType, only : mlcanopy_type
    !
    ! !ARGUMENTS:
    integer, intent(in) :: num_filter               ! Number of patches in filter
    integer, intent(in) :: filter(:)                ! Patch filter
    type(mlcanopy_type), intent(inout) :: mlcanopy_inst
    !
    ! !LOCAL VARIABLES:
    integer  :: fp                                  ! Filter index
    integer  :: p                                   ! Patch index for CLM g/l/c/p hierarchy
    integer  :: ic                                  ! Aboveground layer index
    real(r8) :: dtime                               ! Model time step (s)
    real(r8) :: dew                                 ! Water (kg H2O/m2)
    !---------------------------------------------------------------------

    associate ( &
                                                    ! *** Input ***
    ncan      => mlcanopy_inst%ncan_canopy     , &  ! Number of aboveground layers
    dpai      => mlcanopy_inst%dpai_profile    , &  ! Canopy layer plant area index (m2/m2)
    fracsun   => mlcanopy_inst%fracsun_profile , &  ! Canopy layer sunlit fraction (-)
    trleaf    => mlcanopy_inst%trleaf_leaf     , &  ! Leaf transpiration flux (mol H2O/m2 leaf/s)
    evleaf    => mlcanopy_inst%evleaf_leaf     , &  ! Leaf evaporation flux (mol H2O/m2 leaf/s)
                                                    ! *** Input/Output ***
    h2ocan    => mlcanopy_inst%h2ocan_profile    &  ! Canopy layer intercepted water (kg H2O/m2)
    )

    dtime = dtime_substep

    do fp = 1, num_filter
       p = filter(fp)
       do ic = 1, ncan(p)

          if (dpai(p,ic) > 0._r8) then

             ! Add dew from both evaporation and transpiration

             dew = (evleaf(p,ic,isun) + trleaf(p,ic,isun)) * fracsun(p,ic) * dpai(p,ic) * mmh2o * dtime
             if (dew < 0._r8) then
                h2ocan(p,ic) = h2ocan(p,ic) - dew
             end if

             dew = (evleaf(p,ic,isha) + trleaf(p,ic,isha)) * (1._r8 - fracsun(p,ic)) * dpai(p,ic) * mmh2o * dtime
             if (dew < 0._r8) then
                h2ocan(p,ic) = h2ocan(p,ic) - dew
             end if

             ! Evaporate intercepted water

             if (evleaf(p,ic,isun) > 0._r8) then
                h2ocan(p,ic) = h2ocan(p,ic) - evleaf(p,ic,isun) * fracsun(p,ic) * dpai(p,ic) * mmh2o * dtime
             end if

             if (evleaf(p,ic,isha) > 0._r8) then
                h2ocan(p,ic) = h2ocan(p,ic) - evleaf(p,ic,isha) * (1._r8 - fracsun(p,ic)) * dpai(p,ic) * mmh2o * dtime
             end if

             ! Do not allow h2ocan to go negative (not used)

!            h2ocan(p,ic) = max (0._r8, h2ocan(p,ic))

          end if

       end do
    end do

    end associate
  end subroutine CanopyEvaporation

end module MLCanopyWaterMod
