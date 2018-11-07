module SoilWaterPlantSinkMod

   use clm_varctl       , only : use_hydrstress
   use decompMod        , only : bounds_type
   use shr_kind_mod          , only : r8 => shr_kind_r8
   use shr_log_mod           , only : errMsg => shr_log_errMsg
   use abortutils            , only : endrun
   use clm_varctl            , only : iulog
   use landunit_varcon       , only : istsoil,istcrop
   use column_varcon         , only : icol_road_perv
   implicit none

   character(len=*), parameter, private :: sourcefile = &
         __FILE__

contains
   
   subroutine Compute_EffecRootFrac_And_VertTranSink(bounds, num_hydrologyc, &
         filter_hydrologyc, soilstate_inst, canopystate_inst, waterfluxbulk_inst, energyflux_inst)
      
      ! ---------------------------------------------------------------------------------
      ! This is a wrapper for calculating the effective root fraction and soil
      ! water sink due to plant transpiration. 
      ! Calculate Soil Water Sink to Roots over different types
      ! of columns and for different process modules
      ! The super-set of all columns that should have a root water sink
      ! is filter_hydrologyc
      ! There are three groups of columns:
      ! 1) impervious roads, 2) non-natural vegetation and 3) natural vegetation
      ! There are several methods available.
      ! 1) the default version, 2) hydstress version and 3) fates boundary conditions
      !
      ! There are only two quantities that are the result of this routine, and its
      ! children:
      !   waterfluxbulk_inst%qflx_rootsoi_col(c,j)
      !   soilstate_inst%rootr_col(c,j)
      !
      !
      ! ---------------------------------------------------------------------------------

      use SoilStateType       , only : soilstate_type
      use WaterFluxBulkType       , only : waterfluxbulk_type
      use CanopyStateType     , only : canopystate_type
      use EnergyFluxType      , only : energyflux_type
      use ColumnType          , only : col 
      use LandunitType        , only : lun

      ! Arguments
      type(bounds_type)       , intent(in)    :: bounds               ! bounds
      integer                 , intent(in)    :: num_hydrologyc       ! number of column soil points in column filter
      integer                 , intent(in)    :: filter_hydrologyc(num_hydrologyc) ! column filter for soil points
      type(soilstate_type)    , intent(inout) :: soilstate_inst
      type(waterfluxbulk_type)    , intent(inout) :: waterfluxbulk_inst
      type(canopystate_type)  , intent(in)    :: canopystate_inst
      type(energyflux_type)   , intent(in)    :: energyflux_inst

      ! Local Variables
      integer  :: filterc(bounds%endc-bounds%begc+1)           !column filter
      integer  :: num_filterc
      integer  :: num_filterc_tot
      integer  :: fc
      integer  :: c
      integer  :: l

      num_filterc_tot = 0

      ! 1) pervious roads
      num_filterc = 0
      do fc = 1, num_hydrologyc
         c = filter_hydrologyc(fc)
         if (col%itype(c) == icol_road_perv) then
            num_filterc = num_filterc + 1
            filterc(num_filterc) = c
         end if
      end do
      num_filterc_tot = num_filterc_tot+num_filterc
      if(use_hydrstress) then
         call Compute_EffecRootFrac_And_VertTranSink_HydStress_Roads(bounds, &
               num_filterc,filterc, soilstate_inst, waterfluxbulk_inst)
      else
         call Compute_EffecRootFrac_And_VertTranSink_Default(bounds, &
               num_filterc,filterc, soilstate_inst, waterfluxbulk_inst)
      end if


      ! Note: 2 and 3 really don't need to be split.  But I am leaving
      ! it split in case someone wants to calculate uptake in a special
      ! way for a specific LU or coverage type (RGK 04/2017).  Feel
      ! free to consolidate if there are no plans to do such a thing.

         
      ! 2) not ( pervious road or natural vegetation) , everything else
      num_filterc = 0
      do fc = 1, num_hydrologyc
         c = filter_hydrologyc(fc)
         l = col%landunit(c)
         if ( (col%itype(c) /= icol_road_perv) .and. (lun%itype(l) /= istsoil) ) then
            num_filterc = num_filterc + 1
            filterc(num_filterc) = c
         end if
      end do
      num_filterc_tot = num_filterc_tot+num_filterc
      if(use_hydrstress) then
         call Compute_EffecRootFrac_And_VertTranSink_HydStress(bounds, &
               num_filterc, filterc, waterfluxbulk_inst, soilstate_inst, &
               canopystate_inst, energyflux_inst)
      else
         call Compute_EffecRootFrac_And_VertTranSink_Default(bounds, &
               num_filterc,filterc, soilstate_inst, waterfluxbulk_inst)
      end if
      

      ! 3) Natural vegetation
      num_filterc = 0
      do fc = 1, num_hydrologyc
         c = filter_hydrologyc(fc)
         l = col%landunit(c)
         if ( (lun%itype(l) == istsoil) ) then
            num_filterc = num_filterc + 1
            filterc(num_filterc) = c
         end if
      end do
      num_filterc_tot = num_filterc_tot+num_filterc
      if (use_hydrstress) then
         call Compute_EffecRootFrac_And_VertTranSink_HydStress(bounds, &
              num_filterc, filterc, waterfluxbulk_inst, soilstate_inst, &
              canopystate_inst,energyflux_inst)
      else
         call Compute_EffecRootFrac_And_VertTranSink_Default(bounds, &
              num_filterc,filterc, soilstate_inst, waterfluxbulk_inst)
      end if

      if (num_hydrologyc /= num_filterc_tot) then
          write(iulog,*) 'The total number of columns flagged to root water uptake'
          write(iulog,*) 'did not match the total number calculated'
          write(iulog,*) 'This is likely a problem with the interpretation of column/lu filters.'
          call endrun(msg=errMsg(sourcefile, __LINE__))
      end if


      return
   end subroutine Compute_EffecRootFrac_And_VertTranSink

   ! ====================================================================================

   subroutine Compute_EffecRootFrac_And_VertTranSink_HydStress_Roads(bounds, &
         num_filterc,filterc, soilstate_inst, waterfluxbulk_inst)
      
      use SoilStateType    , only : soilstate_type
      use WaterFluxBulkType    , only : waterfluxbulk_type
      use clm_varpar       , only : nlevsoi
      use clm_varpar       , only : max_patch_per_col
      use PatchType        , only : patch
      use ColumnType       , only : col

      ! Arguments
      type(bounds_type)       , intent(in)    :: bounds      
      integer                 , intent(in)    :: num_filterc
      integer                 , intent(in)    :: filterc(:) 
      type(soilstate_type)    , intent(inout) :: soilstate_inst
      type(waterfluxbulk_type)    , intent(inout) :: waterfluxbulk_inst

      ! Locals
      integer :: j
      integer :: c
      integer :: fc
      integer :: pi
      integer :: p
      real(r8) :: temp(bounds%begc:bounds%endc)                         ! accumulator for rootr weighting


      associate(& 
            qflx_rootsoi_col    => waterfluxbulk_inst%qflx_rootsoi_col    , & ! Output: [real(r8) (:,:) ]  
                                                                          ! vegetation/soil water exchange (mm H2O/s) (+ = to atm)
            qflx_tran_veg_patch => waterfluxbulk_inst%qflx_tran_veg_patch , & ! Input:  [real(r8) (:)   ]  
                                                                          ! vegetation transpiration (mm H2O/s) (+ = to atm) 
            qflx_tran_veg_col   => waterfluxbulk_inst%qflx_tran_veg_col   , & ! Input:  [real(r8) (:)   ]  
                                                                          ! vegetation transpiration (mm H2O/s) (+ = to atm)
            rootr_patch         => soilstate_inst%rootr_patch         , & ! Input:  [real(r8) (:,:) ]  
                                                                          ! effective fraction of roots in each soil layer  
            rootr_col           => soilstate_inst%rootr_col             & ! Output: [real(r8) (:,:) ]  
                                                                          !effective fraction of roots in each soil layer  
            )

        ! First step is to calculate the column-level effective rooting
        ! fraction in each soil layer. This is done outside the usual
        ! PATCH-to-column averaging routines because it is not a simple
        ! weighted average of the PATCH level rootr arrays. Instead, the
        ! weighting depends on both the per-unit-area transpiration
        ! of the PATCH and the PATCHEs area relative to all PATCHES.
        
        temp(bounds%begc : bounds%endc) = 0._r8
    

        do j = 1, nlevsoi
           do fc = 1, num_filterc
              c = filterc(fc)
              rootr_col(c,j) = 0._r8
           end do
        end do
        
        do pi = 1,max_patch_per_col
           do j = 1,nlevsoi
              do fc = 1, num_filterc
                 c = filterc(fc)
                 if (pi <= col%npatches(c)) then
                    p = col%patchi(c) + pi - 1
                    if (patch%active(p)) then
                       rootr_col(c,j) = rootr_col(c,j) + rootr_patch(p,j) * &
                             qflx_tran_veg_patch(p) * patch%wtcol(p)
                    end if
                 end if
              end do
           end do
           do fc = 1, num_filterc
              c = filterc(fc)
              if (pi <= col%npatches(c)) then
                 p = col%patchi(c) + pi - 1
                 if (patch%active(p)) then
                    temp(c) = temp(c) + qflx_tran_veg_patch(p) * patch%wtcol(p)
                 end if
              end if
           end do
        end do

   
        do j = 1, nlevsoi
           do fc = 1, num_filterc
              c = filterc(fc)
              if (temp(c) /= 0._r8) then
                 rootr_col(c,j) = rootr_col(c,j)/temp(c)
              end if
              qflx_rootsoi_col(c,j) = rootr_col(c,j)*qflx_tran_veg_col(c)
           end do
        end do
      end associate
      return
   end subroutine Compute_EffecRootFrac_And_VertTranSink_HydStress_Roads
   
   ! ==================================================================================
   
   subroutine Compute_EffecRootFrac_And_VertTranSink_HydStress( bounds, &
           num_filterc, filterc, waterfluxbulk_inst, soilstate_inst, &
           canopystate_inst, energyflux_inst)


        !
        !USES:
        use decompMod        , only : bounds_type
        use clm_varpar       , only : nlevsoi
        use clm_varpar       , only : max_patch_per_col
        use SoilStateType    , only : soilstate_type
        use WaterFluxBulkType    , only : waterfluxbulk_type
        use CanopyStateType  , only : canopystate_type
        use PatchType        , only : patch
        use ColumnType       , only : col
        use clm_varctl       , only : iulog
        use PhotosynthesisMod, only : plc, params_inst
        use column_varcon    , only : icol_road_perv
        use shr_infnan_mod   , only : isnan => shr_infnan_isnan
        use EnergyFluxType   , only : energyflux_type
        !
        ! !ARGUMENTS:
        type(bounds_type)    , intent(in)    :: bounds          ! bounds
        integer              , intent(in)    :: num_filterc     ! number of column soil points in column filter
        integer              , intent(in)    :: filterc(:)      ! column filter for soil points
        type(waterfluxbulk_type) , intent(inout) :: waterfluxbulk_inst
        type(soilstate_type) , intent(inout) :: soilstate_inst
        type(canopystate_type) , intent(in)  :: canopystate_inst
        type(energyflux_type), intent(in)    :: energyflux_inst
        !
        ! !LOCAL VARIABLES:
        integer  :: p,c,fc,j                                              ! do loop indices
        integer  :: pi                                                    ! patch index
        real(r8) :: temp(bounds%begc:bounds%endc)                         ! accumulator for rootr weighting
        real(r8) :: grav2                 ! soil layer gravitational potential relative to surface (mm H2O)
        integer , parameter :: soil=1,root=4  ! index values
        !-----------------------------------------------------------------------   
        
        associate(&
              k_soil_root         => soilstate_inst%k_soil_root_patch   , & ! Input:  [real(r8) (:,:) ]  
                                                                            ! soil-root interface conductance (mm/s)
              qflx_phs_neg_col    => waterfluxbulk_inst%qflx_phs_neg_col    , & ! Input:  [real(r8) (:)   ]  n
                                                                            ! net neg hydraulic redistribution flux(mm H2O/s)
              qflx_tran_veg_col   => waterfluxbulk_inst%qflx_tran_veg_col   , & ! Input:  [real(r8) (:)   ]  
                                                                            ! vegetation transpiration (mm H2O/s) (+ = to atm)
              qflx_tran_veg_patch => waterfluxbulk_inst%qflx_tran_veg_patch , & ! Input:  [real(r8) (:)   ]  
                                                                            ! vegetation transpiration (mm H2O/s) (+ = to atm)
              qflx_rootsoi_col    => waterfluxbulk_inst%qflx_rootsoi_col    , & ! Output: [real(r8) (:)   ]
                                                                            ! col root and soil water 
                                                                            ! exchange [mm H2O/s] [+ into root]
              rootr_col           => soilstate_inst%rootr_col           , & ! Input:  [real(r8) (:,:) ]
                                                                            ! effective fraction of roots in each soil layer
              rootr_patch         => soilstate_inst%rootr_patch         , & ! Input:  [real(r8) (:,:) ]  
                                                                            ! effective fraction of roots in each soil layer
              smp                 => soilstate_inst%smp_l_col           , & ! Input:  [real(r8) (:,:) ]  soil matrix pot. [mm]
              frac_veg_nosno      => canopystate_inst%frac_veg_nosno_patch , & ! Input:  [integer  (:)  ] 
                                                                            ! fraction of vegetation not 
                                                                            ! covered by snow (0 OR 1) [-]  
              z                   => col%z                              , & ! Input: [real(r8) (:,:) ]  layer node depth (m)
              vegwp               => canopystate_inst%vegwp_patch         & ! Input: [real(r8) (:,:) ]  vegetation water 
                                                                            ! matric potential (mm)
              )
          
          do fc = 1, num_filterc
             c = filterc(fc)
             qflx_phs_neg_col(c) = 0._r8
             
             do j = 1, nlevsoi
                grav2 = z(c,j) * 1000._r8
                temp(c) = 0._r8
                do pi = 1,max_patch_per_col
                   if (pi <= col%npatches(c)) then
                      p = col%patchi(c) + pi - 1
                      if (patch%active(p).and.frac_veg_nosno(p)>0) then 
                         if (patch%wtcol(p) > 0._r8) then
                            temp(c) = temp(c) + k_soil_root(p,j) &
                                  * (smp(c,j) - vegwp(p,4) - grav2)* patch%wtcol(p)
                         endif
                      end if
                   end if
                end do
                qflx_rootsoi_col(c,j)= temp(c)
                
                if (temp(c) < 0._r8) qflx_phs_neg_col(c) = qflx_phs_neg_col(c) + temp(c)
             end do
             
             ! Back out the effective root density
             if( sum(qflx_rootsoi_col(c,:))>0.0_r8 ) then
                do j = 1, nlevsoi
                   rootr_col(c,j) = qflx_rootsoi_col(c,j)/sum( qflx_rootsoi_col(c,:))
                end do
             else
                rootr_col(c,:) = 0.0_r8
             end if
          end do
          
        end associate

        return
     end subroutine Compute_EffecRootFrac_And_VertTranSink_HydStress
     
     ! ==================================================================================

     subroutine Compute_EffecRootFrac_And_VertTranSink_Default(bounds, num_filterc, &
           filterc, soilstate_inst, waterfluxbulk_inst)

    !
    ! Generic routine to apply transpiration as a sink condition that
    ! is vertically distributed over the soil column. Should be
    ! applicable to any Richards solver that is not coupled to plant
    ! hydraulics.
    !
    !USES:
    use decompMod        , only : bounds_type
    use shr_kind_mod     , only : r8 => shr_kind_r8
    use clm_varpar       , only : nlevsoi, max_patch_per_col
    use SoilStateType    , only : soilstate_type
    use WaterFluxBulkType    , only : waterfluxbulk_type
    use PatchType        , only : patch
    use ColumnType       , only : col
    use clm_varctl       , only : use_hydrstress
    use column_varcon    , only : icol_road_perv
    !
    ! !ARGUMENTS:
    type(bounds_type)    , intent(in)    :: bounds                          ! bounds
    integer              , intent(in)    :: num_filterc                     ! number of column soil points in column filter
    integer              , intent(in)    :: filterc(num_filterc)            ! column filter for soil points
    type(waterfluxbulk_type) , intent(inout) :: waterfluxbulk_inst
    type(soilstate_type) , intent(inout) :: soilstate_inst
    !
    ! !LOCAL VARIABLES:
    integer  :: p,c,fc,j                                              ! do loop indices
    integer  :: pi                                                    ! patch index
    real(r8) :: temp(bounds%begc:bounds%endc)                         ! accumulator for rootr weighting
    associate(& 
          qflx_rootsoi_col    => waterfluxbulk_inst%qflx_rootsoi_col    , & ! Output: [real(r8) (:,:) ]  
                                                                        ! vegetation/soil water exchange (m H2O/s) (+ = to atm)
          qflx_tran_veg_patch => waterfluxbulk_inst%qflx_tran_veg_patch , & ! Input:  [real(r8) (:)   ]  
                                                                        ! vegetation transpiration (mm H2O/s) (+ = to atm) 
          qflx_tran_veg_col   => waterfluxbulk_inst%qflx_tran_veg_col   , & ! Input:  [real(r8) (:)   ]  
                                                                        ! vegetation transpiration (mm H2O/s) (+ = to atm)
          rootr_patch         => soilstate_inst%rootr_patch         , & ! Input: [real(r8) (:,:) ]
                                                                        ! effective fraction of roots in each soil layer  
          rootr_col           => soilstate_inst%rootr_col             & ! Output: [real(r8) (:,:) ]  
                                                                        ! effective fraction of roots in each soil layer  
          )
      
      ! First step is to calculate the column-level effective rooting
      ! fraction in each soil layer. This is done outside the usual
      ! PATCH-to-column averaging routines because it is not a simple
      ! weighted average of the PATCH level rootr arrays. Instead, the
      ! weighting depends on both the per-unit-area transpiration
      ! of the PATCH and the PATCHEs area relative to all PATCHES.
      
      temp(bounds%begc : bounds%endc) = 0._r8
      
      do j = 1, nlevsoi
         do fc = 1, num_filterc
            c = filterc(fc)
            rootr_col(c,j) = 0._r8
         end do
      end do
      
      do pi = 1,max_patch_per_col
         do j = 1,nlevsoi
            do fc = 1, num_filterc
               c = filterc(fc)
               if (pi <= col%npatches(c)) then
                  p = col%patchi(c) + pi - 1
                  if (patch%active(p)) then
                     rootr_col(c,j) = rootr_col(c,j) + rootr_patch(p,j) * &
                           qflx_tran_veg_patch(p) * patch%wtcol(p)
                  end if
               end if
            end do
         end do
         do fc = 1, num_filterc
            c = filterc(fc)
            if (pi <= col%npatches(c)) then
               p = col%patchi(c) + pi - 1
               if (patch%active(p)) then
                  temp(c) = temp(c) + qflx_tran_veg_patch(p) * patch%wtcol(p)
               end if
            end if
         end do
      end do
      
      do j = 1, nlevsoi
         do fc = 1, num_filterc
            c = filterc(fc)
            if (temp(c) /= 0._r8) then
               rootr_col(c,j) = rootr_col(c,j)/temp(c)
            end if
            qflx_rootsoi_col(c,j) = rootr_col(c,j)*qflx_tran_veg_col(c)

         end do
      end do
    end associate
    return
 end subroutine Compute_EffecRootFrac_And_VertTranSink_Default

end module SoilWaterPlantSinkMod

