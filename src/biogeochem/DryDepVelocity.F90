Module DryDepVelocity                                              

  !-----------------------------------------------------------------------  
  !  
  ! Purpose:  
  ! Deposition velocity (m/s) 
  !  
  ! Method:  
  ! This code simulates dry deposition velocities using the Wesely scheme.   
  ! Details of this method can be found in:  
  ! 
  ! M.L Wesely. Parameterization of surface resistances to gaseous dry deposition 
  ! in regional-scale numericl models. 1989. Atmospheric Environment vol.23 No.6  
  ! pp. 1293-1304. 
  ! 
  ! In Wesely (1998) "the magnitude of the dry deposition velocity can be found 
  ! as: 
  ! 
  !  |vd|=(ra+rb+rc)^-1 
  ! 
  ! where ra is the aerodynamic resistance (common to all gases) between a  
  ! specific height and the surface, rb is the quasilaminar sublayer resistance 
  ! (whose only dependence on the porperties of the gas of interest is its 
  ! molecular diffusivity in air), and rc is the bulk surface resistance". 
  ! 
  ! In this subroutine both ra and rb are calculated elsewhere in CLM.  
  ! 
  ! In Wesely (1989) rc is estimated for five seasonal categories and 11 landuse 
  ! types.  For each season and landuse type, Wesely compiled data into a  
  ! look-up-table for several parameters used to calculate rc. In this subroutine 
  ! the same values are used as found in wesely's look-up-tables, the only  
  ! difference is that this subroutine uses a CLM generated LAI to select values 
  ! from the look-up-table instead of seasonality.  Inaddition, Wesely(1989)  
  ! land use types are "mapped" into CLM patch types.
  ! 
  ! Subroutine written to operate at the patch level. 
  ! 
  ! Output: 
  ! 
  ! vd(n_species) !Dry deposition velocity [m s-1] for each molecule or species 
  !  
  ! Author: Beth Holland and  James Sulzman 
  ! 
  ! Modified: Francis Vitt -- 30 Mar 2007
  ! Modified: Maria Val Martin -- 15 Jan 2014
  !  Corrected major bugs in the leaf and stomatal resitances. The code is now
  !  coupled to LAI and Rs uses the Ball-Berry Scheme. Also, corrected minor 
  !  bugs in rlu and rcl calculations. Added 
  !  no vegetation removal for CO. See README for details and 
  !     Val Martin et al., 2014 GRL for major corrections
  ! Modified: Louisa Emmons -- 30 November 2017
  !  Corrected the equation calculating stomatal resistance from rssun and rssha,
  !     and removed factor that scaled Rs to match observations
  !
  !----------------------------------------------------------------------- 

  use shr_log_mod          , only : errMsg => shr_log_errMsg
  use shr_kind_mod         , only : r8 => shr_kind_r8 
  use abortutils           , only : endrun
  use clm_time_manager     , only : get_nstep, get_curr_date, get_curr_time 
  use spmdMod              , only : masterproc
  use seq_drydep_mod       , only : n_drydep, drydep_list 
  use seq_drydep_mod       , only : drydep_method, DD_XLND
  use seq_drydep_mod       , only : index_o3=>o3_ndx, index_o3a=>o3a_ndx, index_so2=>so2_ndx, index_h2=>h2_ndx
  use seq_drydep_mod       , only : index_co=>co_ndx, index_ch4=>ch4_ndx, index_pan=>pan_ndx
  use seq_drydep_mod       , only : index_xpan=>xpan_ndx
  use decompMod            , only : bounds_type
  use clm_varcon           , only : namep
  use atm2lndType          , only : atm2lnd_type
  use CanopyStateType      , only : canopystate_type
  use FrictionVelocityMod  , only : frictionvel_type
  use PhotosynthesisMod    , only : photosyns_type
  use WaterStateBulkType       , only : waterstatebulk_type
  use WaterDiagnosticBulkType       , only : waterdiagnosticbulk_type
  use GridcellType         , only : grc                
  use LandunitType         , only : lun                
  use PatchType            , only : patch                
  !
  implicit none 
  private
  !
  public :: depvel_compute
  !
  type, public :: drydepvel_type

     real(r8), pointer, public  :: velocity_patch (:,:) ! Dry Deposition Velocity
     real(r8), pointer, private :: rs_drydep_patch (:)  ! Stomatal resistance associated with dry deposition velocity for Ozone

   contains

     procedure , public  :: Init 
     procedure , private :: InitAllocate 
     procedure , private :: InitHistory

  end type drydepvel_type
  !----------------------------------------------------------------------- 

  character(len=*), parameter, private :: sourcefile = &
       __FILE__

CONTAINS 

  !------------------------------------------------------------------------
  subroutine Init(this, bounds)

    class(drydepvel_type) :: this
    type(bounds_type), intent(in) :: bounds  

    call this%InitAllocate(bounds)
    call this%InitHistory(bounds)

  end subroutine Init

  !------------------------------------------------------------------------
  subroutine InitAllocate(this, bounds)
    !
    ! !USES:
    use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
    use seq_drydep_mod , only : n_drydep, drydep_method, DD_XLND
    !
    ! !ARGUMENTS:
    class(drydepvel_type) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer :: begp, endp
    !------------------------------------------------------------------------

    begp = bounds%begp; endp= bounds%endp

    ! Dry Deposition Velocity 
    if ( n_drydep > 0 .and. drydep_method == DD_XLND )then
       allocate(this%velocity_patch(begp:endp, n_drydep));  this%velocity_patch(:,:) = nan 
       allocate(this%rs_drydep_patch(begp:endp))         ;  this%rs_drydep_patch(:)  = nan 
    end if

  end subroutine InitAllocate

  !-----------------------------------------------------------------------
  subroutine InitHistory(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize history output fields for dry deposition diagnositics
    !
    ! !USES 
    use clm_varcon     , only : spval
    use histFileMod    , only : hist_addfld1d
    use seq_drydep_mod , only : mapping
    !
    ! !ARGUMENTS:
    class(drydepvel_type) :: this
    type(bounds_type), intent(in) :: bounds
    real(r8), pointer  :: ptr_1d(:)  ! pointer to 1d patch array
    !
    ! !LOCAL VARIABLES
    integer :: ispec 
    integer :: begp, endp
    !---------------------------------------------------------------------

    begp = bounds%begp; endp = bounds%endp

    if ( n_drydep == 0 .or. drydep_method /= DD_XLND ) return

    do ispec=1,n_drydep 
       if(mapping(ispec) <= 0) cycle 

       this%velocity_patch(begp:endp,ispec)= spval
       ptr_1d => this%velocity_patch(begp:endp,ispec)
       call hist_addfld1d ( fname='DRYDEPV_'//trim(drydep_list(ispec)), units='cm/sec',  &
            avgflag='A', long_name='Dry Deposition Velocity', &
            ptr_patch=ptr_1d, default='inactive' )
    end do

    this%rs_drydep_patch(begp:endp)= spval
    call hist_addfld1d ( fname='RS_DRYDEP_O3', units='s/m',  &
         avgflag='A', long_name='Stomatal Resistance Associated with Ozone Dry Deposition Velocity', &
         ptr_patch=this%rs_drydep_patch, default='inactive' )

  end subroutine InitHistory

  !----------------------------------------------------------------------- 
  subroutine depvel_compute( bounds, &
       atm2lnd_inst, canopystate_inst, waterstatebulk_inst, waterdiagnosticbulk_inst, frictionvel_inst, &
       photosyns_inst, drydepvel_inst)
    !
    ! !DESCRIPTION:
    ! computes the dry deposition velocity of tracers
    !
    ! !USES:
    use shr_const_mod  , only : tmelt => shr_const_tkfrz
    use seq_drydep_mod , only : seq_drydep_setHCoeff, mapping, drat, foxd
    use seq_drydep_mod , only : rcls, h2_a, h2_b, h2_c, ri, rac, rclo, rlu, rgss, rgso
    use landunit_varcon, only : istsoil, istice_mec, istdlak, istwet
    use clm_varctl     , only : iulog
    use pftconMod      , only : noveg, ndllf_evr_tmp_tree, ndllf_evr_brl_tree
    use pftconMod      , only : ndllf_dcd_brl_tree, nbrdlf_evr_trp_tree
    use pftconMod      , only : nbrdlf_evr_tmp_tree, nbrdlf_dcd_trp_tree
    use pftconMod      , only : nbrdlf_dcd_tmp_tree, nbrdlf_dcd_brl_tree
    use pftconMod      , only : nbrdlf_evr_shrub, nbrdlf_dcd_tmp_shrub
    use pftconMod      , only : nbrdlf_dcd_brl_shrub,nc3_arctic_grass
    use pftconMod      , only : nc3_nonarctic_grass, nc4_grass, nc3crop
    use pftconMod      , only : nc3irrig, npcropmin, npcropmax
    use clm_varcon     , only : spval

    !
    ! !ARGUMENTS:
    type(bounds_type)      , intent(in)    :: bounds  
    type(atm2lnd_type)     , intent(in)    :: atm2lnd_inst
    type(canopystate_type) , intent(in)    :: canopystate_inst
    type(waterstatebulk_type)  , intent(in)    :: waterstatebulk_inst
    type(waterdiagnosticbulk_type)  , intent(in)    :: waterdiagnosticbulk_inst
    type(frictionvel_type) , intent(in)    :: frictionvel_inst
    type(photosyns_type)   , intent(in)    :: photosyns_inst
    type(drydepvel_type)   , intent(inout) :: drydepvel_inst
    !
    ! !LOCAL VARIABLES:
    integer  :: c
    real(r8) :: soilw, var_soilw, fact_h2, dv_soil_h2
    integer  :: pi,g, l
    integer  :: ispec 
    integer  :: length 
    integer  :: wesveg       !wesely vegegation index  
    integer  :: clmveg       !clm veg index from ivegtype 
    integer  :: i 
    integer  :: index_season !seasonal index based on LAI.  This indexs wesely data tables 
    integer  :: nstep        !current step 
    integer  :: indexp 

    real(r8) :: pg           ! surface pressure 
    real(r8) :: tc           ! temperature in celsius  
    real(r8) :: es           ! saturation vapor pressur 
    real(r8) :: ws           ! saturation mixing ratio 
    real(r8) :: rmx          ! resistance by vegetation 
    real(r8) :: qs           ! saturation specific humidity 
    real(r8) :: dewm         ! multiplier for rs when dew occurs 
    real(r8) :: crs          ! multiplier to calculate crs 
    real(r8) :: rdc          ! part of lower canopy resistance 
    real(r8) :: rain         ! rain fall 
    real(r8) :: spec_hum     ! specific humidity 
    real(r8) :: solar_flux   ! solar radiation(direct beam) W/m2 
    real(r8) :: lat          ! latitude in degrees 
    real(r8) :: lon          ! longitude in degrees 
    real(r8) :: sfc_temp     ! surface temp 
    real(r8) :: minlai       ! minimum of monthly lai 
    real(r8) :: maxlai       ! maximum of monthly lai 
    real(r8) :: rds         ! resistance for aerosols

    !mvm 11/30/2013
    real(r8) :: rlu_lai      ! constant to calculate rlu over bulk canopy
    
    logical  :: has_dew
    logical  :: has_rain
    real(r8), parameter :: rain_threshold      = 1.e-7_r8  ! of the order of 1cm/day expressed in m/s

    ! local arrays: dependent on species only 
    real(r8), dimension(n_drydep) :: rsmx  !vegetative resistance (plant mesophyll) 
    real(r8), dimension(n_drydep) :: rclx !lower canopy resistance  
    real(r8), dimension(n_drydep) :: rlux !vegetative resistance (upper canopy) 
    real(r8), dimension(n_drydep) :: rgsx !gournd resistance 
    real(r8), dimension(n_drydep) :: heff   
    real(r8) :: rs    ! stomatal resistance associated with dry deposition velocity (s/m)
    real(r8) :: rc    !combined surface resistance 
    real(r8) :: cts   !correction to flu rcl and rgs for frost 
    real(r8) :: rlux_o3 !to calculate O3 leaf resistance in dew/rain conditions

    ! constants 
    real(r8), parameter :: slope = 0._r8      ! Used to calculate  rdc in (lower canopy resistance) 
    integer, parameter :: wveg_unset = -1     ! Unset Wesley vegetation type
    character(len=32), parameter :: subname = "depvel_compute"

    ! jfl : mods for PAN
    real(r8)                  :: dv_pan
    real(r8) :: c0_pan(11) = (/ 0.000_r8, 0.006_r8, 0.002_r8, 0.009_r8, 0.015_r8, &
                                0.006_r8, 0.000_r8, 0.000_r8, 0.000_r8, 0.002_r8, 0.002_r8 /)
    real(r8) :: k_pan (11) = (/ 0.000_r8, 0.010_r8, 0.005_r8, 0.004_r8, 0.003_r8, &
                                0.005_r8, 0.000_r8, 0.000_r8, 0.000_r8, 0.075_r8, 0.002_r8 /)
    !----------------------------------------------------------------------- 

    if ( n_drydep == 0 .or. drydep_method /= DD_XLND ) return

    associate(                                                    & 
         forc_solai =>    atm2lnd_inst%forc_solai_grc           , & ! Input:  [real(r8) (:,:) ] direct beam radiation (visible only)             
         forc_solad =>    atm2lnd_inst%forc_solad_grc           , & ! Input:  [real(r8) (:,:) ] direct beam radiation (visible only)             
         forc_t     =>    atm2lnd_inst%forc_t_downscaled_col    , & ! Input:  [real(r8) (:)   ] downscaled atmospheric temperature (Kelvin)                   
         forc_q     =>    atm2lnd_inst%forc_q_downscaled_col    , & ! Input:  [real(r8) (:)   ] downscaled atmospheric specific humidity (kg/kg)              
         forc_psrf  =>    atm2lnd_inst%forc_pbot_downscaled_col , & ! Input:  [real(r8) (:)   ] downscaled surface pressure (Pa)                              
         forc_rain  =>    atm2lnd_inst%forc_rain_downscaled_col , & ! Input:  [real(r8) (:)   ] downscaled rain rate [mm/s]                                   

         h2osoi_vol =>    waterstatebulk_inst%h2osoi_vol_col        , & ! Input:  [real(r8) (:,:) ] volumetric soil water (0<=h2osoi_vol<=watsat)   
         snow_depth =>    waterdiagnosticbulk_inst%snow_depth_col        , & ! Input:  [real(r8) (:)   ] snow height (m)                                   

         ram1       =>    frictionvel_inst%ram1_patch           , & ! Input:  [real(r8) (:)   ] aerodynamical resistance                           
         rb1        =>    frictionvel_inst%rb1_patch            , & ! Input:  [real(r8) (:)   ] leaf boundary layer resistance [s/m]               
         vds        =>    frictionvel_inst%vds_patch            , & ! Input:  [real(r8) (:)   ] aerodynamical resistance                           

         rssun      =>    photosyns_inst%rssun_patch            , & ! Input:  [real(r8) (:)   ] stomatal resistance                                
         rssha      =>    photosyns_inst%rssha_patch            , & ! Input:  [real(r8) (:)   ] shaded stomatal resistance (s/m)                   

         fsun       =>    canopystate_inst%fsun_patch           , & ! Input:  [real(r8) (:)   ] sunlit fraction of canopy                          
         elai       =>    canopystate_inst%elai_patch           , & ! Input:  [real(r8) (:)   ] one-sided leaf area index with burying by snow     
         mlaidiff   =>    canopystate_inst%mlaidiff_patch       , & ! Input:  [real(r8) (:)   ] difference in lai between month one and month two  
         annlai     =>    canopystate_inst%annlai_patch         , & ! Input:  [real(r8) (:,:) ] 12 months of monthly lai from input data set     
         
         velocity   =>    drydepvel_inst%velocity_patch         , & ! Output: [real(r8) (:,:) ] cm/sec
         rs_drydep  =>    drydepvel_inst%rs_drydep_patch          & ! Output: [real(r8) (:)   ]  stomatal resistance associated with Ozone dry deposition velocity (s/m)
         )

      !_________________________________________________________________ 
      ! Begin loop through patches

      pft_loop: do pi = bounds%begp,bounds%endp
         l = patch%landunit(pi)
         
         active: if (patch%active(pi)) then

            c = patch%column(pi)
            g = patch%gridcell(pi)
            pg         = forc_psrf(c)  
            spec_hum   = forc_q(c)
            rain       = forc_rain(c) 
            sfc_temp   = forc_t(c) 
            solar_flux = forc_solad(g,1) 
            lat        = grc%latdeg(g) 
            lon        = grc%londeg(g) 
            clmveg     = patch%itype(pi) 
            soilw      = h2osoi_vol(c,1)

            !map CLM veg type into Wesely veg type  
            wesveg = wveg_unset 
            if (clmveg == noveg                               ) wesveg = 8 
            if (clmveg == ndllf_evr_tmp_tree                  ) wesveg = 5 
            if (clmveg == ndllf_evr_brl_tree                  ) wesveg = 5 
            if (clmveg == ndllf_dcd_brl_tree                  ) wesveg = 5 
            if (clmveg == nbrdlf_evr_trp_tree                 ) wesveg = 4 
            if (clmveg == nbrdlf_evr_tmp_tree                 ) wesveg = 4 
            if (clmveg == nbrdlf_dcd_trp_tree                 ) wesveg = 4 
            if (clmveg == nbrdlf_dcd_tmp_tree                 ) wesveg = 4 
            if (clmveg == nbrdlf_dcd_brl_tree                 ) wesveg = 4 
            if (clmveg == nbrdlf_evr_shrub                    ) wesveg = 11 
            if (clmveg == nbrdlf_dcd_tmp_shrub                ) wesveg = 11 
            if (clmveg == nbrdlf_dcd_brl_shrub                ) wesveg = 11 
            if (clmveg == nc3_arctic_grass                    ) wesveg = 3 
            if (clmveg == nc3_nonarctic_grass                 ) wesveg = 3 
            if (clmveg == nc4_grass                           ) wesveg = 3 
            if (clmveg == nc3crop                             ) wesveg = 2 
            if (clmveg == nc3irrig                            ) wesveg = 2 
            if (clmveg >= npcropmin .and. clmveg <= npcropmax ) wesveg = 2 
            if (wesveg == wveg_unset )then
               write(iulog,*) 'clmveg = ', clmveg, 'lun%itype = ', lun%itype(l)
               call endrun(decomp_index=pi, clmlevel=namep, &
                    msg='ERROR: Not able to determine Wesley vegetation type'//&
                    errMsg(sourcefile, __LINE__))
            end if

            ! create seasonality index used to index wesely data tables from LAI,  Bascially 
            !if elai is between max lai from input data and half that max the index_season=1 


            !mail1j and mlai2j are the two monthly lai values pulled from a CLM input data set 
            !/fs/cgd/csm/inputdata/lnd/clm2/rawdata/mksrf_lai.nc.  lai for dates in the middle  
            !of the month are interpolated using using these values and stored in the variable  
            !elai (done elsewhere).  If the difference between mlai1j and mlai2j is greater 
            !than zero it is assumed to be fall and less than zero it is assumed to be spring. 

            !wesely seasonal "index_season" 
            ! 1 - midsummer with lush vegetation 
            ! 2 - Autumn with unharvested cropland 
            ! 3 - Late autumn after frost, no snow 
            ! 4 - Winter, snow on ground and subfreezing 
            ! 5 - Transitional spring with partially green short annuals 


            !mlaidiff=jan-feb 
            minlai=minval(annlai(:,pi)) 
            maxlai=maxval(annlai(:,pi)) 

            index_season = -1

            if ( lun%itype(l) /= istsoil )then
               if ( lun%itype(l) == istice_mec ) then
                  wesveg       = 8
                  index_season = 4
               elseif ( lun%itype(l) == istdlak ) then
                  wesveg       = 7
                  index_season = 4
               elseif ( lun%itype(l) == istwet ) then
                  wesveg       = 9
                  index_season = 2
               elseif ( lun%urbpoi(l) ) then
                  wesveg       = 1
                  index_season = 2
               end if
            else if ( snow_depth(c) > 0 ) then
               index_season = 4
            else if(elai(pi) > 0.5_r8*maxlai) then  
               index_season = 1  
            endif

            if (index_season<0) then 
               if (elai(pi) < (minlai+0.05*(maxlai-minlai))) then  
                  index_season = 3 
               endif
            endif

            if (index_season<0) then 
               if (mlaidiff(pi) > 0.0_r8) then 
                  index_season = 2 
               elseif (mlaidiff(pi) < 0.0_r8) then 
                  index_season = 5 
               elseif (mlaidiff(pi).eq.0.0_r8) then 
                  index_season = 3 
               endif
            endif

            if (index_season<0) then 
               call endrun('ERROR: not able to determine season'//errmsg(sourcefile, __LINE__))
            endif

            ! saturation specific humidity 
            ! 
            es = 611_r8*exp(5414.77_r8*((1._r8/tmelt)-(1._r8/sfc_temp))) 
            ws = .622_r8*es/(pg-es) 
            qs = ws/(1._r8+ws) 

            has_dew = .false.
            if( qs <= spec_hum ) then
               has_dew = .true.
            end if
            if( sfc_temp < tmelt ) then
               has_dew = .false.
            end if

            has_rain = rain > rain_threshold

            if ( has_dew .or. has_rain ) then
               dewm = 3._r8
            else
               dewm = 1._r8
            end if

            !Define tc
            tc = sfc_temp - tmelt

            ! 
            ! rdc (lower canopy res) 
            ! 
            rdc=100._r8*(1._r8+1000._r8/(solar_flux+10._r8))/(1._r8+1000._r8*slope) 

            ! surface resistance : depends on both land type and species 
            ! land types are computed seperately, then resistance is computed as average of values 
            ! following wesely rc=(1/(rs+rm) + 1/rlu +1/(rdc+rcl) + 1/(rac+rgs))**-1 

            !******************************************************* 
            call seq_drydep_setHCoeff( sfc_temp, heff(:n_drydep) )
            !********************************************************* 

            species_loop1: do ispec=1, n_drydep
               if(mapping(ispec) <= 0) cycle 

               if(ispec.eq.index_o3.or.ispec.eq.index_o3a.or.ispec.eq.index_so2) then 
                  rmx=0._r8
               else 
                  rmx=1._r8/((heff(ispec)/3000._r8)+(100._r8*foxd(ispec))) 
               endif

               ! correction for frost 
               cts = 1000._r8*exp( -tc - 4._r8 ) 

               !ground resistance  
               rgsx(ispec) = 1._r8/((heff(ispec)/(1.e5_r8*(rgss(index_season,wesveg)+cts))) + & 
                    (foxd(ispec)/(rgso(index_season,wesveg)+cts))) 

               !-------------------------------------------------------------------------------------
               ! special case for H2 and CO;; CH4 is set ot a fraction of dv(H2)
               !-------------------------------------------------------------------------------------
               if( ispec == index_h2 .or. ispec == index_co .or. ispec == index_ch4 ) then

                  if( ispec == index_co ) then
                     fact_h2 = 1.0_r8
                  elseif ( ispec == index_h2 ) then
                     fact_h2 = 0.5_r8
                  elseif ( ispec == index_ch4 ) then
                     fact_h2 = 50.0_r8
                  end if

                  !-------------------------------------------------------------------------------------
                  ! no deposition on snow, ice, desert, and water
                  !-------------------------------------------------------------------------------------
                  if( wesveg == 1 .or. wesveg == 7 .or. wesveg == 8 .or. index_season == 4 ) then
                     rgsx(ispec) = spval
                  else
                     var_soilw = max( .1_r8,min( soilw,.3_r8 ) )
                     if( wesveg == 3 ) then
                        var_soilw = log( var_soilw )
                     end if
                     dv_soil_h2 = h2_c(wesveg) + var_soilw*(h2_b(wesveg) + var_soilw*h2_a(wesveg))
                     if( dv_soil_h2 > 0._r8 ) then
                        rgsx(ispec) = fact_h2/(dv_soil_h2*1.e-4_r8)
                     end if
                  end if
               end if

               !-------------------------------------------------------------------------------------
               ! no deposition on water or no vegetation or snow (elai<=0)
               !-------------------------------------------------------------------------------------

               no_dep: if( wesveg == 7 .or. elai(pi).le.0_r8 ) then !mvm 11/26/2013
                  rclx(ispec)  = spval
                  rsmx(ispec)  = spval
                  rlux(ispec)  = spval
                  rs           = spval
               else
                  
                  !Stomatal resistance  
                  
                  ! fvitt -- at midnight rssun and/or rssha can be zero in some places which sets rs to zero 
                  ! --- this fix prevents divide by zero error (when rsmx is zero)
                  if (rssun(pi)>0._r8 .and. rssun(pi)<1.e30 .and. rssha(pi)>0._r8 .and. rssha(pi)<1.e30 ) then
                     !LKE: corrected rs to add rssun and rssha in parallel (11/30/2017)
                     rs=1._r8/(fsun(pi)*elai(pi)/rssun(pi) + (1.-fsun(pi))*elai(pi)/rssha(pi))
                  else
                     rs=spval
                  endif

                  rsmx(ispec) = rs*drat(ispec)+rmx 

                  ! Leaf resistance
                  !MVM: adjusted rlu by LAI to get leaf resistance over bulk canopy (gao and wesely, 1995)
                  rlu_lai=cts+rlu(index_season,wesveg)/elai(pi)
                  rlux(ispec) = rlu_lai/(1.e-5_r8*heff(ispec)+foxd(ispec))      

                  !Lower canopy resistance 
                  rclx(ispec) = 1._r8/((heff(ispec)/(1.e5_r8*(rcls(index_season,wesveg)+cts))) + &
                       (foxd(ispec)/(rclo(index_season,wesveg)+cts)))

                  !-----------------------------------
                  !mvm 11/30/2013: special case for CO
                  !Dry deposition of CO and hydrocarbons is negligibly 
                  !small in vegetation [Mueller and Brasseur, 1995].
                  !------------------------------------
                  if( ispec == index_co ) then
                     rclx(ispec)   = spval
                     rsmx(ispec)   = spval
                     rlux(ispec)   = spval
                  endif

                  !--------------------------------------------
                  ! jfl : special case for PAN
                  !--------------------------------------------
                  if( ispec == index_pan ) then
                     dv_pan =  c0_pan(wesveg) * (1._r8 - exp(-k_pan(wesveg)*(rs*drat(ispec))*1.e-2_r8 ))

                     if( dv_pan > 0._r8 .and. index_season /= 4 ) then
                        rsmx(ispec) = ( 1._r8/dv_pan )
                     end if
                  end if

               endif no_dep
               if ( ispec == index_o3 )then
                  rs_drydep(pi) = rs
               end if

            end do species_loop1


            !----------------------------------------------
            !Adjustment for dew and rain in leaf resitances
            !---------------------------------------------
            ! no effect over water            
            no_water: if( wesveg.ne.7 ) then
               !MVM: effect only on vegetated areas (elai> 0)
               with_LAI: if (elai(pi).gt.0._r8) then

                  ! 
                  ! no effect if sfc_temp < O C 
                  ! 
                  non_freezing: if(sfc_temp.gt.tmelt) then
                     if( has_dew ) then 
                        rlu_lai=cts+rlu(index_season,wesveg)/elai(pi)
                        rlux_o3 = 1._r8/((1._r8/3000._r8)+(1._r8/(3._r8*rlu_lai))) 

                        if (index_o3 > 0) then
                           rlux(index_o3) = rlux_o3
                        endif
                        if (index_o3a > 0) then
                           rlux(index_o3a) = rlux_o3
                        endif
                     endif

                     if(has_rain) then 
                        rlu_lai=cts+rlu(index_season,wesveg)/elai(pi)
                        rlux_o3 = 1._r8/((1._r8/1000._r8)+(1._r8/(3._r8*rlu_lai)))

                        if (index_o3 > 0) then
                           rlux(index_o3) = rlux_o3
                        endif
                        if (index_o3a > 0) then
                           rlux(index_o3a) = rlux_o3
                        endif
                     endif

                     species_loop2: do ispec=1,n_drydep 
                        if(mapping(ispec).le.0) cycle 
                        if(ispec.ne.index_o3.and.ispec.ne.index_o3a.and.ispec.ne.index_so2) then 

                           if( has_dew .or. has_rain) then
                              rlu_lai=cts+rlu(index_season,wesveg)/elai(pi)
                              rlux(ispec)=1._r8/((1._r8/(3._r8*rlu_lai))+ & 
                                   (1.e-7_r8*heff(ispec))+(foxd(ispec)/rlux_o3))
                           endif

                        elseif(ispec.eq.index_so2) then 

                           if( has_dew ) then
                              rlux(ispec) = 100._r8
                           endif

                           if(has_rain) then 
                              rlu_lai=cts+rlu(index_season,wesveg)/elai(pi)
                              rlux(ispec) = 1._r8/((1._r8/5000._r8)+(1._r8/(3._r8*rlu_lai))) 
                           endif

                           if( has_dew .or. has_rain ) then
                              !MVM:rlux=50 for SO2 in dew or rain only for *urban land* type surfaces.
                              if (wesveg.eq.1) then
                                 rlux(ispec)=50._r8
                              endif
                           endif
                        end if
                        !mvm 11/30/2013: special case for CO
                        if( ispec.eq.index_co ) then
                           rlux(ispec) = spval
                        endif
                     end do species_loop2
                  endif non_freezing
               endif with_LAI
            endif no_water

            ! resistance for aerosols 
            rds = 1._r8/vds(pi)

            species_loop3: do ispec=1,n_drydep 
               if(mapping(ispec) <= 0) cycle 

               ! 
               ! compute rc 
               ! 
               rc = 1._r8/((1._r8/rsmx(ispec))+(1._r8/rlux(ispec)) + & 
                    (1._r8/(rdc+rclx(ispec)))+(1._r8/(rac(index_season,wesveg)+rgsx(ispec))))
               rc = max( 10._r8, rc)
               !
               ! assume no surface resistance for SO2 over water
               !
               if ( drydep_list(ispec) == 'SO2' .and. wesveg == 7 ) then
                  rc = 0._r8
               end if

               select case( drydep_list(ispec) )
               case ( 'SO4' )
                  velocity(pi,ispec) = (1._r8/(ram1(pi)+rds))*100._r8
               case ( 'NH4','NH4NO3','XNH4NO3' )
                  velocity(pi,ispec) = (1._r8/(ram1(pi)+0.5_r8*rds))*100._r8
               case ( 'Pb' )
                  velocity(pi,ispec) = 0.2_r8
               case ( 'CB1', 'CB2', 'OC1', 'OC2', 'SOAM', 'SOAI', 'SOAT', 'SOAB', 'SOAX' )
                  velocity(pi,ispec) = 0.10_r8
               case ( 'SO2' )
                  velocity(pi,ispec) = (1._r8/(ram1(pi)+rb1(pi)+rc))*200._r8
               case default
                  velocity(pi,ispec) = (1._r8/(ram1(pi)+rb1(pi)+rc))*100._r8
               end select
            end do species_loop3
         endif active
      end do pft_loop

     end associate 

  end subroutine depvel_compute

end module DryDepVelocity
