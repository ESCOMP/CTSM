module CNAnnualUpdateMod
  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module for updating annual summation variables
  !
  ! !USES:
  use shr_kind_mod        , only : r8 => shr_kind_r8
  use decompMod           , only : bounds_type
  use CNVegCarbonFluxType , only : cnveg_carbonflux_type
  use CNvegStateType      , only : cnveg_state_type
  use PatchType           , only : patch
  use filterColMod        , only : filter_col_type, col_filter_from_filter_and_logical_array
  use ColumnType          , only : col
  !
  implicit none
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public:: CNAnnualUpdate
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine CNAnnualUpdate(bounds, num_bgc_soilc, filter_bgc_soilc, num_bgc_vegp, filter_bgc_vegp, &
       cnveg_state_inst, cnveg_carbonflux_inst)
    !
    ! !DESCRIPTION:
    ! On the radiation time step, update annual summation variables
    !
    ! !USES:
    use clm_time_manager, only: get_step_size_real, get_curr_days_per_year
    use clm_varcon      , only: secspday
    use SubgridAveMod   , only: p2c
    !
    ! !ARGUMENTS:
    type(bounds_type)           , intent(in)    :: bounds  
    integer                     , intent(in)    :: num_bgc_soilc         ! number of bgc soil columns in filter
    integer                     , intent(in)    :: filter_bgc_soilc(:)   ! filter for bgc soil columns
    integer                     , intent(in)    :: num_bgc_vegp         ! number of bgc veg patches in filter
    integer                     , intent(in)    :: filter_bgc_vegp(:)   ! filter for bgc veg patches
    type(cnveg_state_type)      , intent(inout) :: cnveg_state_inst
    type(cnveg_carbonflux_type) , intent(inout) :: cnveg_carbonflux_inst
    !
    ! !LOCAL VARIABLES:
    integer :: c,p          ! indices
    integer :: fp,fc        ! lake filter indices
    real(r8):: secspyear
    real(r8):: dt           ! radiation time step (seconds)
    logical :: end_of_year(bounds%begc:bounds%endc) ! whether each column has reached the end of the year, according to its own annsum_counter
    type(filter_col_type) :: filter_endofyear_c
    !-----------------------------------------------------------------------

    dt = get_step_size_real()
    secspyear = get_curr_days_per_year() * secspday

    do fc = 1,num_bgc_soilc
       c = filter_bgc_soilc(fc)
       if(.not.col%is_fates(c))then
          cnveg_state_inst%annsum_counter_col(c) = cnveg_state_inst%annsum_counter_col(c) + dt
          if (cnveg_state_inst%annsum_counter_col(c) >= secspyear) then
             end_of_year(c) = .true.
             cnveg_state_inst%annsum_counter_col(c) = 0._r8
          else
             end_of_year(c) = .false.
          end if
       end if
    end do
    

    do fp = 1,num_bgc_vegp
       p = filter_bgc_vegp(fp)
       c = patch%column(p)

       if (end_of_year(c) .and. .not.col%is_fates(c)) then

          ! update annual plant ndemand accumulator
          cnveg_state_inst%annsum_potential_gpp_patch(p)  = cnveg_state_inst%tempsum_potential_gpp_patch(p)
          cnveg_state_inst%tempsum_potential_gpp_patch(p) = 0._r8

          ! update annual total N retranslocation accumulator
          cnveg_state_inst%annmax_retransn_patch(p)  = cnveg_state_inst%tempmax_retransn_patch(p)
          cnveg_state_inst%tempmax_retransn_patch(p) = 0._r8

          ! update annual average 2m air temperature accumulator
          cnveg_state_inst%annavg_t2m_patch(p)  = cnveg_state_inst%tempavg_t2m_patch(p)
          cnveg_state_inst%tempavg_t2m_patch(p) = 0._r8

          ! update annual NPP accumulator, convert to annual total
          cnveg_carbonflux_inst%annsum_npp_patch(p) = cnveg_carbonflux_inst%tempsum_npp_patch(p) * dt
          cnveg_carbonflux_inst%tempsum_npp_patch(p) = 0._r8

          ! update annual litfall accumulator, convert to annual total
          cnveg_carbonflux_inst%annsum_litfall_patch(p) = cnveg_carbonflux_inst%tempsum_litfall_patch(p) * dt
          cnveg_carbonflux_inst%tempsum_litfall_patch(p) = 0._r8

       end if
    end do

    ! Get column-level averages, just for the columns that have reached their personal end-of-year
    if(num_bgc_vegp>0)then
       filter_endofyear_c = col_filter_from_filter_and_logical_array( &
            bounds = bounds, &
            num_orig = num_bgc_soilc, &
            filter_orig = filter_bgc_soilc, &
            logical_col = end_of_year(bounds%begc:bounds%endc))
       
       call p2c(bounds, filter_endofyear_c%num, filter_endofyear_c%indices, &
            cnveg_carbonflux_inst%annsum_npp_patch(bounds%begp:bounds%endp), &
            cnveg_carbonflux_inst%annsum_npp_col(bounds%begc:bounds%endc))
       
       call p2c(bounds, filter_endofyear_c%num, filter_endofyear_c%indices, &
            cnveg_state_inst%annavg_t2m_patch(bounds%begp:bounds%endp), &
            cnveg_state_inst%annavg_t2m_col(bounds%begc:bounds%endc))
    end if
    
 end subroutine CNAnnualUpdate

end module CNAnnualUpdateMod
