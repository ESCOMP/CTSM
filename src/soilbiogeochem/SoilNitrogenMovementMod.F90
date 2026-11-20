module SoilNitrogenMovementMod
 
  !------------------------------------------------------------------------
  ! DESCRIPTION
  ! implementation of Soil-Water-Atmosphere-Plant (SWAP3.2) 
  ! and Pantakar 1980 algorithm to Community Land Model 
  ! SWAP3.2 website: https://edepot.wur.nl/39776
  ! Author: Jinmu Luo, Cornell EAS, April 1 2024
  
  use decompMod                        , only : bounds_type
  use shr_kind_mod                     , only : r8 => shr_kind_r8
  use shr_infnan_mod                   , only : isnan => shr_infnan_isnan
  use shr_infnan_mod                   , only : isinf => shr_infnan_isinf
  use clm_varctl                       , only : iulog
  use spmdMod                          , only : masterproc
  use abortutils                       , only : endrun
  use clm_time_manager                 , only : get_step_size_real
  use clm_time_manager                 , only : get_curr_date
  use SoilBiogeochemNitrogenFluxType   , only : soilbiogeochem_nitrogenflux_type
  use SoilBiogeochemNitrogenStateType  , only : soilbiogeochem_nitrogenstate_type
  use WaterStateBulkType               , only : waterstatebulk_type
  use SoilStatetype                    , only : soilstate_type
  use SoilHydrologyType                , only : soilhydrology_type
  use ColumnType                       , only : col
 
  !
  implicit none
  private
  !
  ! PUBLIC MEMBER FUNCTIONS 
  public SoilNitrogenMovement

  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  !------------------------------------------------------------------------

  contains

  !------------------------------------------------------------------------
  subroutine SoilNitrogenMovement(bounds, num_bgc_soilc, filter_bgc_soilc, waterstatebulk_inst, &
             soilstate_inst, soilhydrology_inst, soilbiogeochem_nitrogenflux_inst, soilbiogeochem_nitrogenstate_inst)
    !
    ! implementation of the advection-diffusion algorithm in Patankar1980
    ! 
    ! This part of code is only designed for fast aqueous transport and leaching of inorganic nitrate 
    ! no sources and other sinks are included in this module  
    ! leaching flux is taken out of soil pool
    ! Author: Jinmu Luo, Cornell EAS, Apr 11 2024                             
    !                                     
    !USES:
    use decompMod         , only : bounds_type
    use clm_varpar        , only : nlevdecomp, nlevgrnd
    use clm_time_manager  , only : get_step_size_real, get_curr_date
    use clm_varcon        , only : zsoi, zisoi, dzsoi_decomp, mmh2o_to_m3h2o_per_m2
    use ColumnType        , only : col
    use clm_varctl        , only : use_bedrock
    use TridiagonalMod    , only : Tridiagonal
    !ARGUMENTS: 
    type(bounds_type)                       , intent(in)    :: bounds               ! bounds
    integer                                 , intent(in)    :: num_bgc_soilc        ! number of soil columns in filter
    integer                                 , intent(in)    :: filter_bgc_soilc(:)  ! filter for soil columns
    type(waterstatebulk_type)               , intent(in)    :: waterstatebulk_inst
    type(soilstate_type)                    , intent(in)    :: soilstate_inst
    type(soilhydrology_type)                , intent(in)    :: soilhydrology_inst
    type(soilbiogeochem_nitrogenflux_type)  , intent(inout) :: soilbiogeochem_nitrogenflux_inst
    type(soilbiogeochem_nitrogenstate_type) , intent(inout) :: soilbiogeochem_nitrogenstate_inst

    !LOCAL VARIABLES:
    integer  :: c,fc,j                                            ! do loop indices
    integer  :: year, mon, day, tod 
    integer  :: jtop(bounds%begc:bounds%endc)                     ! top level at each column
    real(r8) :: dtime                                             ! land model time step (sec)
    real(r8) :: wafc, wafc2                                       ! Fraction of water that is liquid by mass
    real(r8) :: dispersion_length = 0.1_r8                        ! dispersion length (m), Jury et al., 1991
    real(r8) :: theta, thetasat                                   ! soil water and soil water at the saturation level
    real(r8) :: no3_diffusivity_in_water = 1.7e-9_r8              ! Molecular diffusivity of NO3- in water, m2/s
    real(r8) :: dissolve_frac = 1.0_r8                            ! dissolve fraction
    real(r8) :: flux_component_gridpoint_ahead                    ! A function in Patankar 1980, figure 5.6
    real(r8) :: peclet_num                                        ! Peclet number in Patankar 1980, foumula 5.18
    real(r8) :: peclet_num_in, peclet_num_out                     ! temporary Peclet numbers
    real(r8) :: qflx_in, qflx_out                                 ! water fluxes same as qin, qout but in m3 H2O/m2/s
    real(r8) :: dz_node(1:nlevdecomp+1)                           ! difference between nodes
    real(r8) :: mass_old(bounds%begc:bounds%endc)                 ! Temporal column mass, g/m2
    real(r8) :: mass_new(bounds%begc:bounds%endc)                 ! Temporal column mass, g/m2
    real(r8) :: swliq(bounds%begc:bounds%endc,1:nlevdecomp)       ! volumetric liquid soil water [m3/m3], hardwired to 1 for non-transport layers and layers below bedrock
    real(r8) :: total_diffusivity(bounds%begc:bounds%endc,1:nlevdecomp+1)  ! Total diffusivity
    real(r8) :: a_tri(bounds%begc:bounds%endc,0:nlevdecomp+1)     ! "a" vector for tridiagonal matrix
    real(r8) :: b_tri(bounds%begc:bounds%endc,0:nlevdecomp+1)     ! "b" vector for tridiagonal matrix
    real(r8) :: c_tri(bounds%begc:bounds%endc,0:nlevdecomp+1)     ! "c" vector for tridiagonal matrix
    real(r8) :: r_tri(bounds%begc:bounds%endc,0:nlevdecomp+1)     ! "r" vector for tridiagonal solution 
    real(r8) :: conc_trcr(bounds%begc:bounds%endc,0:nlevdecomp+1) ! temporary for concentration, g/m3H2O
  
    ! set up the A function, table 5.2 in Patankar 1980 has multiple A function. 
    ! Notes: According to Table 5.2, here we use the "Power Law" version of the function
    !        A is a dimensionless coefficient described in equation 5.37 of Patankar (1980)
    !        The same identical function appears in CLM's SoilBiogeochemLittVertTranspMod.F90 as aaa(pe)
    !        Patankar (1980) is posted here: https://github.com/ESCOMP/CTSM/pull/2992#discussion_r2294809728
    flux_component_gridpoint_ahead(peclet_num) = max(0._r8, ( 1._r8 - 0.1_r8 * abs(peclet_num) )**5 )

    associate(&
         h2osoi_vol          => waterstatebulk_inst%h2osoi_vol_col                       , & ! Input:  [real(r8) (:,:) ]  volumetric soil water (0<=h2osoi_vol<=watsat) [m3/m3]
         h2osoi_liq          => waterstatebulk_inst%h2osoi_liq_col                       , & ! Input:  [real(r8) (:,:) ]  col liquid water (kg/m2) 
         h2osoi_ice          => waterstatebulk_inst%h2osoi_ice_col                       , & ! Input:  [real(r8) (:,:) ]  col ice lens (kg/m2)
         watsat              => soilstate_inst%watsat_col                                , & ! Input:  [real(r8) (:,:) ]  volumetric soil water at saturation (porosity)
         qout                => soilhydrology_inst%qout_col                              , & ! Input:  [real(r8) (:,:) ]  soil water out of the bottom, mm h2o/s 
         qin                 => soilhydrology_inst%qin_col                               , & ! Input:  [real(r8) (:,:) ]  soil water into the bottom, mm h2o/s  
         smin_no3_vr         => soilbiogeochem_nitrogenstate_inst%smin_no3_vr_col        , & ! Inout:  [real(r8) (:,:) ]  soil nitrate concentration, gN/m3
         smin_no3_leached_vr => soilbiogeochem_nitrogenflux_inst%smin_no3_leached_vr_col   & ! Output: [real(r8) (:,:) ]  rate of mineral NO3 leaching (gN/m3/s)  
         )
   
    !Get the size of model time step  
    dtime = get_step_size_real()
    call get_curr_date(year, mon, day, tod)

    ! Preparing for the necessary parameters, like Delta z, Jtop, and the total diffusivity coefficients.
    do j = 1, nlevdecomp + 1
       if (j < nlevdecomp) then
          dz_node(j) = zsoi(j+1) - zsoi(j)
       end if
       do fc = 1, num_bgc_soilc
          c = filter_bgc_soilc(fc)
          jtop(c) = 0
          if (j < nlevdecomp) then
             wafc = h2osoi_liq(c,j)/(h2osoi_liq(c,j) + h2osoi_ice(c,j)) 
             wafc2 = h2osoi_liq(c,j+1)/(h2osoi_liq(c,j+1) + h2osoi_ice(c,j+1))
             swliq(c,j) = wafc * h2osoi_vol(c,j) 
             swliq(c,j+1) = wafc2 * h2osoi_vol(c,j+1) 
             if (swliq(c,j) == 0._r8 .or. swliq(c,j+1) == 0._r8) then
                theta = 0._r8
                thetasat = 1._r8
             else
                theta = swliq(c,j) + dzsoi_decomp(j)/2 * (swliq(c,j+1) - swliq(c,j))/dz_node(j)
                thetasat = watsat(c,j) + dzsoi_decomp(j)/2 * (watsat(c,j+1) - watsat(c,j))/dz_node(j)
             end if
             ! here we refer the j as the interface of j + zj/2 
             total_diffusivity(c,j) = theta * (no3_diffusivity_in_water * theta**(7/3) * thetasat**(-2))
             total_diffusivity(c,j) = total_diffusivity(c,j) + dispersion_length * mmh2o_to_m3h2o_per_m2 * abs(qout(c,j))
          else
             !no gradient for the last layer
             total_diffusivity(c,j) = total_diffusivity(c,j-1)
          end if
       end do ! Loop for columns 
    end do ! Loop for depths    
    dz_node(nlevdecomp) = dz_node(nlevdecomp-1)
    dz_node(nlevdecomp+1) = dz_node(nlevdecomp) 

    ! Calculate the tridiagonal matrix, using Crank-Nicholson soluiton 
    ! dc/dt = (1-alpha)dF(t+1)/dx + alpha*F(t)/dx, 
    ! alpha=0
    do j = 0, nlevdecomp + 1
       do fc = 1, num_bgc_soilc
          c = filter_bgc_soilc(fc)
          if ( j==0 .or. j==nlevdecomp+1) then
             !atmosphere and bottom layer, no concentration gradient here
             conc_trcr(c,j) = 0._r8
             a_tri(c,j) = 0._r8
             b_tri(c,j) = 1._r8
             c_tri(c,j) = 0._r8
             r_tri(c,j) = 0._r8
          elseif (swliq(c,j) == 0._r8 .or. j > col%nbedrock(c)) then
             ! extremely dry condition and layers beneath the bedrock, no aqueous transport of nitrate
             conc_trcr(c,j) = dissolve_frac * smin_no3_vr(c,j)
             a_tri(c,j) = 0._r8
             b_tri(c,j) = 1._r8
             c_tri(c,j) = 0._r8
             r_tri(c,j) = conc_trcr(c,j)
             swliq(c,j) = 1.0_r8      ! change swliq into 1 to be used in the update session below 
          elseif ( j == 1) then
             ! topmost soil layer, flux only interacts with the layer below it
             conc_trcr(c,j) = dissolve_frac * smin_no3_vr(c,j)/swliq(c,j)
             qflx_out = qout(c,j) * mmh2o_to_m3h2o_per_m2
             peclet_num_out = qflx_out * dz_node(j) / total_diffusivity(c,j)
             a_tri(c,j) = 0._r8
             c_tri(c,j) = -total_diffusivity(c,j) / dz_node(j) * flux_component_gridpoint_ahead(peclet_num_out) - max(-qflx_out, 0._r8)
             b_tri(c,j) = total_diffusivity(c,j) / dz_node(j) * flux_component_gridpoint_ahead(peclet_num_out) + max(qflx_out, 0._r8) + swliq(c,j) / dtime * dzsoi_decomp(j)
             r_tri(c,j) = conc_trcr(c,j)/dtime*swliq(c,j)*dzsoi_decomp(j)
          elseif ( j == col%nbedrock(c)) then
             ! Assume the bottom layer concentration is always zero
             ! This method count the loss at this layer as the leaching flux at the bottom
             conc_trcr(c,j) = 0._r8 
             a_tri(c,j) = 0._r8
             b_tri(c,j) = 1._r8
             c_tri(c,j) = 0._r8
             r_tri(c,j) = conc_trcr(c,j)/dtime*swliq(c,j)*dzsoi_decomp(j)
          else
             ! Active layers from second one to bedrock-1,  concentration should be in gN/m3Water
             conc_trcr(c,j) = dissolve_frac * smin_no3_vr(c,j)/swliq(c,j)
             qflx_in = qin(c,j) * mmh2o_to_m3h2o_per_m2  ! mm H2O/s to m3 H2O/m2/s
             qflx_out = qout(c,j) * mmh2o_to_m3h2o_per_m2
             peclet_num_in = qflx_in * dz_node(j-1) / total_diffusivity(c,j-1)
             peclet_num_out = qflx_out * dz_node(j) / total_diffusivity(c,j)
             a_tri(c,j) = -total_diffusivity(c,j-1) / dz_node(j-1) * flux_component_gridpoint_ahead(peclet_num_in) - max(qflx_in, 0._r8)
             c_tri(c,j) = -total_diffusivity(c,j) / dz_node(j) * flux_component_gridpoint_ahead(peclet_num_out) - max(-qflx_out, 0._r8)
             b_tri(c,j) = total_diffusivity(c,j-1) / dz_node(j-1) * flux_component_gridpoint_ahead(peclet_num_in) + max(-qflx_in, 0._r8) + &
                          total_diffusivity(c,j) / dz_node(j) * flux_component_gridpoint_ahead(peclet_num_out) + max(qflx_out, 0._r8) + swliq(c,j) / dtime * dzsoi_decomp(j)
             r_tri(c,j) = conc_trcr(c,j)/dtime*swliq(c,j)*dzsoi_decomp(j) 
          end if
       end do ! Loop for columns 
    end do ! Loop for depths

    ! solve the tridiagonal matrix 
    call Tridiagonal(bounds, 0, nlevdecomp+1, &
                     jtop(bounds%begc:bounds%endc), &
                     num_bgc_soilc, filter_bgc_soilc, &
                     a_tri(bounds%begc:bounds%endc, :), &
                     b_tri(bounds%begc:bounds%endc, :), &
                     c_tri(bounds%begc:bounds%endc, :), &
                     r_tri(bounds%begc:bounds%endc, :), &
                     conc_trcr(bounds%begc:bounds%endc,0:nlevdecomp+1))

    ! Calculate the leaching flux
    mass_old(bounds%begc:bounds%endc) = 0._r8
    mass_new(bounds%begc:bounds%endc) = 0._r8 
    do j = 1, nlevdecomp
       do fc = 1, num_bgc_soilc
          c = filter_bgc_soilc(fc)
          smin_no3_leached_vr(c,j) = 0._r8
          mass_old(c) = mass_old(c) + smin_no3_vr(c,j)*dzsoi_decomp(j)
          mass_new(c) = mass_new(c) + (smin_no3_vr(c,j) * (1._r8 - dissolve_frac) + conc_trcr(c,j) * swliq(c,j)) * dzsoi_decomp(j)
       end do 
    end do 

    do fc = 1, num_bgc_soilc
       c = filter_bgc_soilc(fc)
       ! g/m3/sec, leaching mass is at the layer above the bedrock
       smin_no3_leached_vr(c, col%nbedrock(c)) = max(0._r8, (mass_old(c) - mass_new(c))/dzsoi_decomp(col%nbedrock(c))/dtime)
    end do 

    ! Update the pools of interest
    do j = 1, nlevdecomp 
       do fc = 1, num_bgc_soilc
          c = filter_bgc_soilc(fc)
          smin_no3_vr(c,j) = smin_no3_vr(c,j) - smin_no3_vr(c,j) * dissolve_frac + conc_trcr(c,j) * swliq(c,j)
          ! Return this leaching flux back to smin_no3 pool, and update will be finished in CNNStateUpdate3Mod 
          if( j == col%nbedrock(c) ) then
            smin_no3_vr(c,j) = smin_no3_vr(c,j) + smin_no3_leached_vr(c,j)*dtime
          end if 
       end do ! loop for columns
    end do !loop for depths 

    end associate


  end subroutine SoilNitrogenMovement


end module SoilNitrogenMovementMod
