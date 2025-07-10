module DistParamType

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Spatially distributed parameter data type allocation and initialization
  ! -------------------------------------------------------- 
  !
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
  use shr_sys_mod    , only : shr_sys_abort
  use clm_varcon     , only : spval, ispval, grlnd
  use clm_varctl     , only : iulog
  use clm_varctl     , only : paramfile, distributed_paramfile
  use ColumnType     , only : col
  use spmdMod        , only : masterproc, mpicom
  use fileutils      , only : getfil, opnfil, getavu, relavu
  use shr_nl_mod     , only : shr_nl_find_group_name
  use shr_mpi_mod    , only : shr_mpi_bcast
  use decompMod      , only : bounds_type
  use ncdio_pio      , only : ncd_pio_closefile, ncd_pio_openfile
  use ncdio_pio      , only : file_desc_t, ncd_inqdid, ncd_inqdlen
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  private
  character(len=*), parameter, private :: sourcefile = __FILE__
  !
  type, public :: distparam_type
     ! CanopyHydrologyMod
     real(r8), pointer :: liq_canopy_storage_scalar    (:)   ! Canopy-storage-of-liquid-water parameter (kg/m2)
     real(r8), pointer :: snow_canopy_storage_scalar   (:)   ! Canopy-storage-of-snow parameter (kg/m2)
     real(r8), pointer :: snowcan_unload_temp_fact     (:)   ! Temperature canopy snow unload scaling (C2 in Eq. 14, Roesch et al. 2001) (K*s)
     real(r8), pointer :: snowcan_unload_wind_fact     (:)   ! Wind canopy snow unload scaling (modifies 1.56e5, where 1.56e5 is C3 in Eq. 15, Roesch et al. 2001) (-)
     real(r8), pointer :: interception_fraction        (:)   ! Fraction of intercepted precipitation (-)
     real(r8), pointer :: maximum_leaf_wetted_fraction (:)   ! Maximum fraction of leaf that may be wet (-)

     ! SoilHydrologyMod
     real(r8), pointer :: aq_sp_yield_min              (:)   ! Minimum aquifer specific yield (unitless)
     real(r8), pointer :: n_baseflow                   (:)   ! Drainage power law exponent (unitless)
     real(r8), pointer :: perched_baseflow_scalar      (:)   ! Scalar multiplier for perched base flow rate (kg/m2/s)
     real(r8), pointer :: e_ice                        (:)   ! Soil ice impedance factor (unitless)
     real(r8), pointer :: baseflow_scalar              (:)   ! Scalar multiplier for base flow rate ()

     ! SaturatedExcessRunoff
     real(r8), pointer :: fff                          (:)   ! Decay factor for fractional saturated area (1/m)
     
     ! initVerticalMod
     real(r8), pointer :: slopebeta                    (:)   ! exponent for microtopography pdf sigma (unitless)
     real(r8), pointer :: slopemax                     (:)   ! max topographic slope for microtopography pdf sigma (unitless)
     real(r8), pointer :: zbedrock_sf                  (:)   ! parameter to scale zbedrock (m)

     ! SnowCoverFractionSwensonLawrence2012Mod
     real(r8), pointer :: n_melt_coef                  (:)   ! SCA shape parameter
     real(r8), pointer :: accum_factor                 (:)   ! Accumulation constant for fractional snow covered area (unitless)

     ! SnowHydrologyMod
     real(r8), pointer :: upplim_destruct_metamorph    (:)   ! Upper limit on destructive metamorphism compaction (kg/m3)

     ! SoilHydrologyInitTimeConstMod
     real(r8), pointer :: pc                           (:)   ! Threshold probability for surface water (unitless)
     real(r8), pointer :: om_frac_sf                   (:)   ! Scale factor for organic matter fraction (unitless)

     ! SoilStateInitTimeConstMod
     real(r8), pointer :: bsw_sf                       (:)   ! Scale factor for bsw (unitless)
     real(r8), pointer :: hksat_sf                     (:)   ! Scale factor for hksat (unitless)
     real(r8), pointer :: sucsat_sf                    (:)   ! Scale factor for sucsat (unitless)
     real(r8), pointer :: watsat_sf                    (:)   ! Scale factor for watsat (unitless)

     ! SurfaceResistanceMod
     real(r8), pointer :: d_max                        (:)   ! Dry surface layer parameter (mm)
     real(r8), pointer :: frac_sat_soil_dsl_init       (:)   ! Fraction of saturated soil for moisture value at which DSL initiates (unitless)

     ! SurfaceWaterMod
     real(r8), pointer :: mu                           (:)   ! Connectivity exponent for surface water (unitless)

     ! WaterDiagnosticBulkType
     real(r8), pointer :: zlnd                         (:)   ! Momentum roughness length for soil, glacier, wetland (m)
     real(r8), pointer :: snw_rds_min                  (:)   ! minimum allowed snow effective radius (also cold "fresh snow" value) [microns]

     ! atm2lndType
     real(r8), pointer :: precip_repartition_nonglc_all_rain_t   (:)   ! Rain temperature threshold for non-glacier landunits (C)
     real(r8), pointer :: precip_repartition_nonglc_all_snow_t   (:)   ! Snow temperature threshold for non-glacier landunits (C)
     
   contains

     procedure, public :: Init
     procedure, public :: readDistributedParams
     procedure, public :: Clean

  end type distparam_type

  type(distparam_type), public, target :: distparams !
  
  !------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  subroutine Init(this, bounds)
    !
    ! !ARGUMENTS:
    class(distparam_type)  :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer :: begc, endc
    integer :: begg, endg
    !------------------------------------------------------------------------

    begc = bounds%begc; endc = bounds%endc
    begg = bounds%begg; endg = bounds%endg

    ! these should be patch; todo
    ! CanopyHydrology
    allocate(this%liq_canopy_storage_scalar    (begc:endc))                     ; this%liq_canopy_storage_scalar    (:) = nan
    allocate(this%snow_canopy_storage_scalar   (begc:endc))                     ; this%snow_canopy_storage_scalar   (:) = nan
    allocate(this%snowcan_unload_temp_fact     (begc:endc))                     ; this%snowcan_unload_temp_fact     (:) = nan
    allocate(this%snowcan_unload_wind_fact     (begc:endc))                     ; this%snowcan_unload_wind_fact     (:) = nan
    allocate(this%interception_fraction        (begc:endc))                     ; this%interception_fraction        (:) = nan
    allocate(this%maximum_leaf_wetted_fraction (begc:endc))                     ; this%maximum_leaf_wetted_fraction (:) = nan

    ! SoilHydrology
    allocate(this%aq_sp_yield_min           (begc:endc))                     ; this%aq_sp_yield_min           (:) = nan
    allocate(this%n_baseflow                (begc:endc))                     ; this%n_baseflow                (:) = nan
    allocate(this%perched_baseflow_scalar   (begc:endc))                     ; this%perched_baseflow_scalar   (:) = nan
    allocate(this%e_ice                     (begc:endc))                     ; this%e_ice                     (:) = nan
    allocate(this%baseflow_scalar           (begc:endc))                     ; this%baseflow_scalar           (:) = nan

    ! SaturatedExcessRunoff
    allocate(this%fff                       (begc:endc))                     ; this%fff                       (:) = nan

    ! initVertical
    allocate(this%slopebeta                 (begg:endc))                     ; this%slopebeta                 (:) = nan
    allocate(this%slopemax                  (begg:endc))                     ; this%slopemax                  (:) = nan
    allocate(this%zbedrock_sf               (begg:endg))                     ; this%zbedrock_sf               (:) = nan

    ! SnowCoverFractionSwensonLawrence2012Mod
    allocate(this%n_melt_coef               (begg:endc))                     ; this%n_melt_coef               (:) = nan
    allocate(this%accum_factor              (begg:endc))                     ; this%accum_factor              (:) = nan

     ! SnowHydrologyMod
    allocate(this%upplim_destruct_metamorph (begg:endc))                     ; this%upplim_destruct_metamorph (:) = nan

    ! SoilHydrologyInitTimeConstMod
    allocate(this%pc                        (begg:endc))                     ; this%pc                        (:) = nan
    allocate(this%om_frac_sf                (begg:endc))                     ; this% om_frac_sf               (:) = nan

    ! SoilStateInitTimeConstMod
    allocate(this%bsw_sf                    (begg:endc))                     ; this%bsw_sf                    (:) = nan
    allocate(this%hksat_sf                  (begg:endc))                     ; this%hksat_sf                  (:) = nan
    allocate(this%sucsat_sf                 (begg:endc))                     ; this%sucsat_sf                 (:) = nan
    allocate(this%watsat_sf                 (begg:endc))                     ; this%watsat_sf                 (:) = nan

    ! SurfaceResistanceMod
    allocate(this%d_max                     (begg:endc))                     ; this%d_max                     (:) = nan
    allocate(this%frac_sat_soil_dsl_init    (begg:endc))                     ; this%frac_sat_soil_dsl_init    (:) = nan

    ! SurfaceWaterMod
    allocate(this%mu                        (begg:endc))                     ; this%mu                        (:) = nan

    ! WaterDiagnosticBulkType
    allocate(this%zlnd                      (begg:endc))                     ; this%zlnd                      (:) = nan
    allocate(this%snw_rds_min               (begg:endc))                     ; this%snw_rds_min               (:) = nan

    ! atm2lndType
    allocate(this%precip_repartition_nonglc_all_rain_t   (begg:endc))        ; this%precip_repartition_nonglc_all_rain_t                (:) = nan
    allocate(this%precip_repartition_nonglc_all_snow_t   (begg:endc))        ; this%precip_repartition_nonglc_all_snow_t                (:) = nan

    !    allocate(this% (begg:endc))                     ; this%             (:) = nan

  end subroutine Init

  !-----------------------------------------------------------------------
  subroutine readDistributedParams(this,bounds)
    !
    ! !USES:
    use ncdio_pio       , only : check_var, ncd_io
    use paramUtilMod    , only : readNcdioScalar
    use abortutils      , only : endrun
    use shr_log_mod     , only : errMsg => shr_log_errMsg
    use domainMod       , only : ldomain
    use clm_varctl      , only : NLFilename => NLFilename_in

    !
    ! !ARGUMENTS:
    implicit none
    class(distparam_type)              :: this
    type(bounds_type), intent(in)      :: bounds
    !
    ! !LOCAL VARIABLES:
    logical            :: readvar               ! whether the variable was found
    integer            :: c, g
    integer            :: dimid ! netCDF dimension id
    integer            :: dim1, dim2            ! dimension to compare
    ! namelist variable names must match file
    real(r8)           :: baseflow_scalar       ! read in namelist
    real(r8)           :: precip_repartition_nonglc_all_rain_t  ! read in namelist
    real(r8)           :: precip_repartition_nonglc_all_snow_t  ! read in namelist
    real(r8)           :: fscalar_in            ! read in - scalar - float
    real(r8), pointer  :: fparam_in(:)          ! read in - 1D - float
    type(file_desc_t)  :: ncids, ncidd          ! pio netCDF file ids
    character(len=256) :: locfns, locfnd        ! local filenames
    integer            :: ierr                  ! error code
    integer            :: unitn                 ! unit for namelist file
    character(len=256) :: nmlname               ! namelist
    character(len=*), parameter :: subname = 'readDistributedParams'
    !--------------------------------------------------------------------

    if (masterproc) then
       write(iulog,*) trim(subname)//' :: reading CLM spatially distributed parameters'
    end if

    ! global (scalar) parameters
    call getfil (paramfile, locfns, 0)
    call ncd_pio_openfile (ncids, trim(locfns), 0)
    ! spatially distributed parameters
    call getfil (distributed_paramfile, locfnd, 0)
    call ncd_pio_openfile (ncidd, trim(locfnd), 0)

    ! should a check be added to make sure dimensions match surface data / mesh files?
    call ncd_inqdid(ncidd,'lon',dimid) 
    call ncd_inqdlen(ncidd,dimid,dim1)
    call ncd_inqdid(ncidd,'lat',dimid) 
    call ncd_inqdlen(ncidd,dimid,dim2)
    if (dim1 /= ldomain%ni .or. dim2 /= ldomain%nj) then
       call endrun(msg="ERROR dimensions of distributed parameter file incompatible with domain dimensions"//errmsg(sourcefile, __LINE__))
    endif
    
    !-----------------------------------------------------------
    ! CanopyHydrology !
    !-----------------------------------------------------------

    ! these should be patch/pft not column... leaving in for now.
    
    ! Canopy-storage-of-liquid-water parameter (kg/m2)
    call check_var(ncidd, 'liq_canopy_storage_scalar', readvar)
    if (.false.) then
       ! if distributed parameter found, overwrite global value
       allocate(fparam_in(bounds%begg:bounds%endg))
       call ncd_io(ncid=ncidd, varname='', flag='read', data=fparam_in, dim1name=grlnd, readvar=readvar)

       do c = bounds%begc,bounds%endc
          g = col%gridcell(c)
          this%liq_canopy_storage_scalar(c) = fparam_in(g)
       enddo

       deallocate(fparam_in)
    else ! use global parameter value
       call readNcdioScalar(ncids, 'liq_canopy_storage_scalar', subname, fscalar_in)
       this%liq_canopy_storage_scalar(:) = fscalar_in
    endif

    ! Canopy-storage-of-snow parameter (kg/m2)
    call check_var(ncidd, 'snow_canopy_storage_scalar', readvar)
    if (.false.) then
       ! if distributed parameter found, overwrite global value
       allocate(fparam_in(bounds%begg:bounds%endg))
       call ncd_io(ncid=ncidd, varname='', flag='read', data=fparam_in, dim1name=grlnd, readvar=readvar)

       do c = bounds%begc,bounds%endc
          g = col%gridcell(c)
          this%snow_canopy_storage_scalar(c) = fparam_in(g)
       enddo

       deallocate(fparam_in)
    else ! use global parameter value
       call readNcdioScalar(ncids, 'snow_canopy_storage_scalar', subname, fscalar_in)
       this%snow_canopy_storage_scalar(:) = fscalar_in
    endif

    ! Temperature canopy snow unload scaling (C2 in Eq. 14, Roesch et al. 2001) (K*s)
    call check_var(ncidd, 'snowcan_unload_temp_fact', readvar)
    if (.false.) then
       ! if distributed parameter found, overwrite global value
       allocate(fparam_in(bounds%begg:bounds%endg))
       call ncd_io(ncid=ncidd, varname='snowcan_unload_temp_fact', flag='read', data=fparam_in, dim1name=grlnd, readvar=readvar)

       do c = bounds%begc,bounds%endc
          g = col%gridcell(c)
          this%snowcan_unload_temp_fact(c) = fparam_in(g)
       enddo

       deallocate(fparam_in)
    else ! use global parameter value
       call readNcdioScalar(ncids, 'snowcan_unload_temp_fact', subname, fscalar_in)
       this%snowcan_unload_temp_fact(:) = fscalar_in
    endif

    ! Wind canopy snow unload scaling (modifies 1.56e5, where 1.56e5 is C3 in Eq. 15, Roesch et al. 2001) (-)
    call check_var(ncidd, 'snowcan_unload_wind_fact', readvar)
    if (.false.) then
       ! if distributed parameter found, overwrite global value
       allocate(fparam_in(bounds%begg:bounds%endg))
       call ncd_io(ncid=ncidd, varname='snowcan_unload_wind_fact', flag='read', data=fparam_in, dim1name=grlnd, readvar=readvar)

       do c = bounds%begc,bounds%endc
          g = col%gridcell(c)
          this%snowcan_unload_wind_fact(c) = fparam_in(g)
       enddo

       deallocate(fparam_in)
    else ! use global parameter value
       call readNcdioScalar(ncids, 'snowcan_unload_wind_fact', subname, fscalar_in)
       this%snowcan_unload_wind_fact(:) = fscalar_in
    endif

    ! Fraction of intercepted precipitation (-)
    call check_var(ncidd, 'interception_fraction', readvar)
    if (.false.) then
       ! if distributed parameter found, overwrite global value
       allocate(fparam_in(bounds%begg:bounds%endg))
       call ncd_io(ncid=ncidd, varname='interception_fraction', flag='read', data=fparam_in, dim1name=grlnd, readvar=readvar)

       do c = bounds%begc,bounds%endc
          g = col%gridcell(c)
          this%interception_fraction(c) = fparam_in(g)
       enddo

       deallocate(fparam_in)
    else ! use global parameter value
       call readNcdioScalar(ncids, 'interception_fraction', subname, fscalar_in)
       this%interception_fraction(:) = fscalar_in
    endif

    ! Maximum fraction of leaf that may be wet (-)
    call check_var(ncidd, 'maximum_leaf_wetted_fraction', readvar)
    if (.false.) then
       ! if distributed parameter found, overwrite global value
       allocate(fparam_in(bounds%begg:bounds%endg))
       call ncd_io(ncid=ncidd, varname='maximum_leaf_wetted_fraction', flag='read', data=fparam_in, dim1name=grlnd, readvar=readvar)

       do c = bounds%begc,bounds%endc
          g = col%gridcell(c)
          this%maximum_leaf_wetted_fraction(c) = fparam_in(g)
       enddo

       deallocate(fparam_in)
    else ! use global parameter value
       call readNcdioScalar(ncids, 'maximum_leaf_wetted_fraction', subname, fscalar_in)
       this%maximum_leaf_wetted_fraction(:) = fscalar_in
    endif

    !-----------------------------------------------------------
    ! SoilHydrology !
    !-----------------------------------------------------------

    ! Minimum aquifer specific yield (unitless)
    call check_var(ncidd, 'aq_sp_yield_min', readvar)
    if (readvar) then
       ! if distributed parameter found, overwrite global value
       allocate(fparam_in(bounds%begg:bounds%endg))
       call ncd_io(ncid=ncidd, varname='aq_sp_yield_min', flag='read', data=fparam_in, dim1name=grlnd, readvar=readvar)

       do c = bounds%begc,bounds%endc
          g = col%gridcell(c)
          this%aq_sp_yield_min(c) = fparam_in(g)
       enddo

       deallocate(fparam_in)
    else ! use global parameter value
       call readNcdioScalar(ncids, 'aq_sp_yield_min', subname, fscalar_in)
       this%aq_sp_yield_min(:) = fscalar_in
    endif
       
    ! Drainage power law exponent (unitless)
    call check_var(ncidd, 'n_baseflow', readvar)
    if (readvar) then
       ! if distributed parameter found, overwrite global value
       allocate(fparam_in(bounds%begg:bounds%endg))
       call ncd_io(ncid=ncidd, varname='n_baseflow', flag='read', data=fparam_in, dim1name=grlnd, readvar=readvar)

       do c = bounds%begc,bounds%endc
          g = col%gridcell(c)
          this%n_baseflow(c) = fparam_in(g)
       enddo

       deallocate(fparam_in)
    else ! use global parameter value
    call readNcdioScalar(ncids, 'n_baseflow', subname, fscalar_in)
       this%n_baseflow(:) = fscalar_in
    endif
    
    ! Scalar multiplier for perched base flow rate (kg/m2/s)
    call check_var(ncidd, 'perched_baseflow_scalar', readvar)
    if (readvar) then
       ! if distributed parameter found, overwrite global value
       allocate(fparam_in(bounds%begg:bounds%endg))
       call ncd_io(ncid=ncidd, varname='perched_baseflow_scalar', flag='read', data=fparam_in, dim1name=grlnd, readvar=readvar)

       do c = bounds%begc,bounds%endc
          g = col%gridcell(c)
          this%perched_baseflow_scalar(c) = fparam_in(g)
       enddo

       deallocate(fparam_in)
    else ! use global parameter value
       call readNcdioScalar(ncids, 'perched_baseflow_scalar', subname, fscalar_in)
       this%perched_baseflow_scalar(:) = fscalar_in
    endif
    
    ! Soil ice impedance factor (unitless)
    call check_var(ncidd, 'e_ice', readvar)
    if (readvar) then
       ! if distributed parameter found, overwrite global value
       allocate(fparam_in(bounds%begg:bounds%endg))
       call ncd_io(ncid=ncidd, varname='e_ice', flag='read', data=fparam_in, dim1name=grlnd, readvar=readvar)

       do c = bounds%begc,bounds%endc
          g = col%gridcell(c)
          this%e_ice(c) = fparam_in(g)
       enddo

       deallocate(fparam_in)
    else ! use global parameter value
       call readNcdioScalar(ncids, 'e_ice', subname, fscalar_in)
       this%e_ice(:) = fscalar_in
    endif

    ! Scalar multiplier for base flow rate ()
    call check_var(ncidd, 'baseflow_scalar', readvar)
    if (readvar) then
       ! if distributed parameter found, overwrite global value
       allocate(fparam_in(bounds%begg:bounds%endg))
       call ncd_io(ncid=ncidd, varname='baseflow_scalar', flag='read', data=fparam_in, dim1name=grlnd, readvar=readvar)

       do c = bounds%begc,bounds%endc
          g = col%gridcell(c)
          this%baseflow_scalar(c) = fparam_in(g)
       enddo

       deallocate(fparam_in)
    else ! global parameter value from namelist
       nmlname = 'soilhydrology_inparm'
       namelist /soilhydrology_inparm/ baseflow_scalar

       if (masterproc) then
          unitn = getavu()
          write(iulog,*) 'Read in '//nmlname//'  namelist'
          call opnfil (NLFilename, unitn, 'F')
          call shr_nl_find_group_name(unitn, nmlname, status=ierr)
          if (ierr == 0) then
             read(unitn, nml=soilhydrology_inparm, iostat=ierr)
             if (ierr /= 0) then
                call endrun(msg="ERROR reading "//nmlname//"namelist"//errmsg(sourcefile, __LINE__))
             end if
          else
             call endrun(msg="ERROR could NOT find "//nmlname//"namelist"//errmsg(sourcefile, __LINE__))
          end if
          call relavu( unitn )
       end if

       call shr_mpi_bcast (baseflow_scalar, mpicom)
       this%baseflow_scalar(:) = baseflow_scalar
    endif

    !-----------------------------------------------------------
    ! SaturatedExcessRunoff !
    !-----------------------------------------------------------
    ! Decay factor for fractional saturated area (1/m)
    call check_var(ncidd, 'fff', readvar)
    if (readvar) then
       ! if distributed parameter found, overwrite global value
       allocate(fparam_in(bounds%begg:bounds%endg))
       call ncd_io(ncid=ncidd, varname='fff', flag='read', data=fparam_in, dim1name=grlnd, readvar=readvar)

       do c = bounds%begc,bounds%endc
          g = col%gridcell(c)
          this%fff(c) = fparam_in(g)
       enddo

       deallocate(fparam_in)
    else ! use global parameter value
       call readNcdioScalar(ncids, 'fff', subname, fscalar_in)
       this%fff(:) = fscalar_in
    endif

    
    !-----------------------------------------------------------
    ! initVertical !
    !-----------------------------------------------------------

    ! exponent for microtopography pdf sigma (unitless)
    call check_var(ncidd, 'slopebeta', readvar)
    if (readvar) then
       ! if distributed parameter found, overwrite global value
       allocate(fparam_in(bounds%begg:bounds%endg))
       call ncd_io(ncid=ncidd, varname='slopebeta', flag='read', data=fparam_in, dim1name=grlnd, readvar=readvar)

       do c = bounds%begc,bounds%endc
          g = col%gridcell(c)
          this%slopebeta(c) = fparam_in(g)
       enddo

       deallocate(fparam_in)
    else ! use global parameter value
       call readNcdioScalar(ncids, 'slopebeta', subname, fscalar_in)
       this%slopebeta(:) = fscalar_in
    endif

    ! max topographic slope for microtopography pdf sigma (unitless)
    call check_var(ncidd, 'slopemax', readvar)
    if (readvar) then
       ! if distributed parameter found, overwrite global value
       allocate(fparam_in(bounds%begg:bounds%endg))
       call ncd_io(ncid=ncidd, varname='slopemax', flag='read', data=fparam_in, dim1name=grlnd, readvar=readvar)

       do c = bounds%begc,bounds%endc
          g = col%gridcell(c)
          this%slopemax(c) = fparam_in(g)
       enddo

       deallocate(fparam_in)
    else ! use global parameter value
       call readNcdioScalar(ncids, 'slopemax', subname, fscalar_in)
       this%slopemax(:) = fscalar_in
    endif

    ! parameter to scale zbedrock (m)
    call check_var(ncidd, 'zbedrock_sf', readvar)
    if (readvar) then
       ! if distributed parameter found, overwrite global value
       allocate(fparam_in(bounds%begg:bounds%endg))
       call ncd_io(ncid=ncidd, varname='zbedrock_sf', flag='read', data=fparam_in, dim1name=grlnd, readvar=readvar)

       do g = bounds%begg,bounds%endg
          this%zbedrock_sf(g) = fparam_in(g)
       enddo

       deallocate(fparam_in)
    else ! use global parameter value
       call readNcdioScalar(ncids, 'zbedrock_sf', subname, fscalar_in)
       this%zbedrock_sf(:) = fscalar_in
    endif
       
    !-----------------------------------------------------------
    ! SnowCoverFractionSwensonLawrence2012 !
    !-----------------------------------------------------------

    ! SCA shape parameter
    call check_var(ncidd, 'n_melt_coef', readvar)
    if (readvar) then
       ! if distributed parameter found, overwrite global value
       allocate(fparam_in(bounds%begg:bounds%endg))
       call ncd_io(ncid=ncidd, varname='n_melt_coef', flag='read', data=fparam_in, dim1name=grlnd, readvar=readvar)

       do c = bounds%begc,bounds%endc
          g = col%gridcell(c)
          this%n_melt_coef(c) = fparam_in(g)
       enddo

       deallocate(fparam_in)
    else ! use global parameter value
       call readNcdioScalar(ncids, 'n_melt_coef', subname, fscalar_in)
       this%n_melt_coef(:) = fscalar_in
    endif
       
    ! Accumulation constant for fractional snow covered area (unitless) 
    call check_var(ncidd, 'accum_factor', readvar)
    if (readvar) then
       ! if distributed parameter found, overwrite global value
       allocate(fparam_in(bounds%begg:bounds%endg))
       call ncd_io(ncid=ncidd, varname='accum_factor', flag='read', data=fparam_in, dim1name=grlnd, readvar=readvar)

       do c = bounds%begc,bounds%endc
          g = col%gridcell(c)
          this%accum_factor(c) = fparam_in(g)
       enddo

       deallocate(fparam_in)
    else ! use global parameter value
       call readNcdioScalar(ncids, 'accum_factor', subname, fscalar_in)
       this%accum_factor(:) = fscalar_in
    endif

    !-----------------------------------------------------------
    ! SnowHydrologyMod !
    !-----------------------------------------------------------

    ! Upper limit on destructive metamorphism compaction (kg/m3)
    call check_var(ncidd, 'upplim_destruct_metamorph', readvar)
    if (readvar) then
       ! if distributed parameter found, overwrite global value
       allocate(fparam_in(bounds%begg:bounds%endg))
       call ncd_io(ncid=ncidd, varname='upplim_destruct_metamorph', flag='read', data=fparam_in, dim1name=grlnd, readvar=readvar)

       do c = bounds%begc,bounds%endc
          g = col%gridcell(c)
          this%upplim_destruct_metamorph(c) = fparam_in(g)
       enddo

       deallocate(fparam_in)
    else ! use global parameter value
       call readNcdioScalar(ncids, 'upplim_destruct_metamorph', subname, fscalar_in)
       this%upplim_destruct_metamorph(:) = fscalar_in
    endif

    !-----------------------------------------------------------
    ! SoilHydrologyInitTimeConstMod !
    !-----------------------------------------------------------

    ! Scale factor for organic matter fraction (unitless)
    call check_var(ncidd, 'om_frac_sf', readvar)
    if (readvar) then
       ! if distributed parameter found, overwrite global value
       allocate(fparam_in(bounds%begg:bounds%endg))
       call ncd_io(ncid=ncidd, varname='om_frac_sf', flag='read', data=fparam_in, dim1name=grlnd, readvar=readvar)

       do c = bounds%begc,bounds%endc
          g = col%gridcell(c)
          this%om_frac_sf(c) = fparam_in(g)
       enddo

       deallocate(fparam_in)
    else ! use global parameter value
       call readNcdioScalar(ncids, 'om_frac_sf', subname, fscalar_in)
       this%om_frac_sf(:) = fscalar_in
    endif

    !-----------------------------------------------------------
    ! SoilStateInitTimeConstMod !
    !-----------------------------------------------------------

    call check_var(ncidd, 'bsw_sf', readvar)
    if (readvar) then
       ! if distributed parameter found, overwrite global value
       allocate(fparam_in(bounds%begg:bounds%endg))
       call ncd_io(ncid=ncidd, varname='bsw_sf', flag='read', data=fparam_in, dim1name=grlnd, readvar=readvar)

       do c = bounds%begc,bounds%endc
          g = col%gridcell(c)
          this%bsw_sf(c) = fparam_in(g)
       enddo

       deallocate(fparam_in)
    else ! use global parameter value
       call readNcdioScalar(ncids, 'bsw_sf', subname, fscalar_in)
       this%bsw_sf(:) = fscalar_in
    endif

    call check_var(ncidd, 'hksat_sf', readvar)
    if (readvar) then
       ! if distributed parameter found, overwrite global value
       allocate(fparam_in(bounds%begg:bounds%endg))
       call ncd_io(ncid=ncidd, varname='hksat_sf', flag='read', data=fparam_in, dim1name=grlnd, readvar=readvar)

       do c = bounds%begc,bounds%endc
          g = col%gridcell(c)
          this%hksat_sf(c) = fparam_in(g)
       enddo

       deallocate(fparam_in)
    else ! use global parameter value
       call readNcdioScalar(ncids, 'hksat_sf', subname, fscalar_in)
       this%hksat_sf(:) = fscalar_in
    endif

    call check_var(ncidd, 'sucsat_sf', readvar)
    if (readvar) then
       ! if distributed parameter found, overwrite global value
       allocate(fparam_in(bounds%begg:bounds%endg))
       call ncd_io(ncid=ncidd, varname='sucsat_sf', flag='read', data=fparam_in, dim1name=grlnd, readvar=readvar)

       do c = bounds%begc,bounds%endc
          g = col%gridcell(c)
          this%sucsat_sf(c) = fparam_in(g)
       enddo

       deallocate(fparam_in)
    else ! use global parameter value
       call readNcdioScalar(ncids, 'sucsat_sf', subname, fscalar_in)
       this%sucsat_sf(:) = fscalar_in
    endif

    call check_var(ncidd, 'watsat_sf', readvar)
    if (readvar) then
       ! if distributed parameter found, overwrite global value
       allocate(fparam_in(bounds%begg:bounds%endg))
       call ncd_io(ncid=ncidd, varname='watsat_sf', flag='read', data=fparam_in, dim1name=grlnd, readvar=readvar)

       do c = bounds%begc,bounds%endc
          g = col%gridcell(c)
          this%watsat_sf(c) = fparam_in(g)
       enddo

       deallocate(fparam_in)
    else ! use global parameter value
       call readNcdioScalar(ncids, 'watsat_sf', subname, fscalar_in)
       this%watsat_sf(:) = fscalar_in
    endif

    !-----------------------------------------------------------
    ! SurfaceResistanceMod !
    !-----------------------------------------------------------

    ! Dry surface layer parameter (mm)
    call check_var(ncidd, 'd_max', readvar)
    if (readvar) then
       ! if distributed parameter found, overwrite global value
       allocate(fparam_in(bounds%begg:bounds%endg))
       call ncd_io(ncid=ncidd, varname='d_max', flag='read', data=fparam_in, dim1name=grlnd, readvar=readvar)

       do c = bounds%begc,bounds%endc
          g = col%gridcell(c)
          this%d_max(c) = fparam_in(g)
       enddo

       deallocate(fparam_in)
    else ! use global parameter value
       call readNcdioScalar(ncids, 'd_max', subname, fscalar_in)
       this%d_max(:) = fscalar_in
    endif

    ! Fraction of saturated soil for moisture value at which DSL initiates (unitless)
    call check_var(ncidd, 'frac_sat_soil_dsl_init', readvar)
    if (readvar) then
       ! if distributed parameter found, overwrite global value
       allocate(fparam_in(bounds%begg:bounds%endg))
       call ncd_io(ncid=ncidd, varname='frac_sat_soil_dsl_init', flag='read', data=fparam_in, dim1name=grlnd, readvar=readvar)

       do c = bounds%begc,bounds%endc
          g = col%gridcell(c)
          this%frac_sat_soil_dsl_init(c) = fparam_in(g)
       enddo

       deallocate(fparam_in)
    else ! use global parameter value
       call readNcdioScalar(ncids, 'frac_sat_soil_dsl_init', subname, fscalar_in)
       this%frac_sat_soil_dsl_init(:) = fscalar_in
    endif

    !-----------------------------------------------------------
    ! SurfaceWaterMod !
    !-----------------------------------------------------------


    ! Threshold probability for surface water (unitless)
    call check_var(ncidd, 'pc', readvar)
    if (readvar) then
       ! if distributed parameter found, overwrite global value
       allocate(fparam_in(bounds%begg:bounds%endg))
       call ncd_io(ncid=ncidd, varname='pc', flag='read', data=fparam_in, dim1name=grlnd, readvar=readvar)

       do c = bounds%begc,bounds%endc
          g = col%gridcell(c)
          this%pc(c) = fparam_in(g)
       enddo

       deallocate(fparam_in)
    else ! use global parameter value
       call readNcdioScalar(ncids, 'pc', subname, fscalar_in)
       this%pc(:) = fscalar_in
    endif

    ! Connectivity exponent for surface water (unitless)
    call check_var(ncidd, 'mu', readvar)
    if (readvar) then
       ! if distributed parameter found, overwrite global value
       allocate(fparam_in(bounds%begg:bounds%endg))
       call ncd_io(ncid=ncidd, varname='mu', flag='read', data=fparam_in, dim1name=grlnd, readvar=readvar)

       do c = bounds%begc,bounds%endc
          g = col%gridcell(c)
          this%mu(c) = fparam_in(g)
       enddo

       deallocate(fparam_in)
    else ! use global parameter value
       call readNcdioScalar(ncids, 'mu', subname, fscalar_in)
       this%mu(:) = fscalar_in
    endif

    !-----------------------------------------------------------
    ! WaterDiagnosticBulkType !
    !-----------------------------------------------------------

    ! Momentum roughness length for soil, glacier, wetland (m) 
    call check_var(ncidd, 'zlnd', readvar)
    if (readvar) then
       ! if distributed parameter found, overwrite global value
       allocate(fparam_in(bounds%begg:bounds%endg))
       call ncd_io(ncid=ncidd, varname='zlnd', flag='read', data=fparam_in, dim1name=grlnd, readvar=readvar)

       do c = bounds%begc,bounds%endc
          g = col%gridcell(c)
          this%zlnd(c) = fparam_in(g)
       enddo

       deallocate(fparam_in)
    else ! use global parameter value
       call readNcdioScalar(ncids, 'zlnd', subname, fscalar_in)
       this%zlnd(:) = fscalar_in
    endif

    ! minimum allowed snow effective radius (also cold "fresh snow" value) [microns] 
    call check_var(ncidd, 'snw_rds_min', readvar)
    if (.false.) then
       ! if distributed parameter found, overwrite global value
       allocate(fparam_in(bounds%begg:bounds%endg))
       call ncd_io(ncid=ncidd, varname='snw_rds_min', flag='read', data=fparam_in, dim1name=grlnd, readvar=readvar)

       do c = bounds%begc,bounds%endc
          g = col%gridcell(c)
          this%snw_rds_min(c) = fparam_in(g)
       enddo

       deallocate(fparam_in)
    else ! use global parameter value
       call readNcdioScalar(ncids, 'snw_rds_min', subname, fscalar_in)
       this%snw_rds_min(:) = fscalar_in
    endif

    !-----------------------------------------------------------
    ! atm2lndType !
    !-----------------------------------------------------------
    call check_var(ncidd, 'precip_repartition_nonglc_all_snow_t', readvar)
    if (readvar) then
       ! if distributed parameter found, overwrite global value
       allocate(fparam_in(bounds%begg:bounds%endg))
       call ncd_io(ncid=ncidd, varname='precip_repartition_nonglc_all_snow_t', flag='read', data=fparam_in, dim1name=grlnd, readvar=readvar)

       do c = bounds%begc,bounds%endc
          g = col%gridcell(c)
          this%precip_repartition_nonglc_all_snow_t(c) = fparam_in(g)
       enddo

       deallocate(fparam_in)
    else ! global parameter value from namelist
       nmlname = 'atm2lnd_inparm'
       namelist /atm2lnd_inparm/ precip_repartition_nonglc_all_snow_t

       if (masterproc) then
          unitn = getavu()
          write(iulog,*) 'Read in '//nmlname//'  namelist'
          call opnfil (NLFilename, unitn, 'F')
          call shr_nl_find_group_name(unitn, nmlname, status=ierr)
          if (ierr == 0) then
             read(unitn, nml=atm2lnd_inparm, iostat=ierr)
             if (ierr /= 0) then
                call endrun(msg="ERROR reading "//nmlname//"namelist"//errmsg(sourcefile, __LINE__))
             end if
          else
             call endrun(msg="ERROR could NOT find "//nmlname//"namelist"//errmsg(sourcefile, __LINE__))
          end if
          call relavu( unitn )
       end if

       call shr_mpi_bcast (precip_repartition_nonglc_all_snow_t, mpicom)
       this%precip_repartition_nonglc_all_snow_t(:) = precip_repartition_nonglc_all_snow_t
    endif
    
    call check_var(ncidd, 'precip_repartition_nonglc_all_rain_t', readvar)
    if (readvar) then
       ! if distributed parameter found, overwrite global value
       allocate(fparam_in(bounds%begg:bounds%endg))
       call ncd_io(ncid=ncidd, varname='precip_repartition_nonglc_all_rain_t', flag='read', data=fparam_in, dim1name=grlnd, readvar=readvar)

       do c = bounds%begc,bounds%endc
          g = col%gridcell(c)
          this%precip_repartition_nonglc_all_rain_t(c) = fparam_in(g)
       enddo

       deallocate(fparam_in)
    else ! global parameter value from namelist
       nmlname = 'atm2lnd_inparm'
       namelist /atm2lnd_inparm/ precip_repartition_nonglc_all_rain_t

       if (masterproc) then
          unitn = getavu()
          write(iulog,*) 'Read in '//nmlname//'  namelist'
          call opnfil (NLFilename, unitn, 'F')
          call shr_nl_find_group_name(unitn, nmlname, status=ierr)
          if (ierr == 0) then
             read(unitn, nml=atm2lnd_inparm, iostat=ierr)
             if (ierr /= 0) then
                call endrun(msg="ERROR reading "//nmlname//"namelist"//errmsg(sourcefile, __LINE__))
             end if
          else
             call endrun(msg="ERROR could NOT find "//nmlname//"namelist"//errmsg(sourcefile, __LINE__))
          end if
          call relavu( unitn )
       end if

       call shr_mpi_bcast (precip_repartition_nonglc_all_rain_t, mpicom)
       this%precip_repartition_nonglc_all_rain_t(:) = precip_repartition_nonglc_all_rain_t
    endif

    ! close parameter files
    call ncd_pio_closefile(ncids)
    call ncd_pio_closefile(ncidd)
    
  end subroutine readDistributedParams

  !------------------------------------------------------------------------
  subroutine Clean(this)
    !
    ! !ARGUMENTS:
    class(distparam_type) :: this
    !------------------------------------------------------------------------

    deallocate(this%liq_canopy_storage_scalar)
    deallocate(this%snow_canopy_storage_scalar)
    deallocate(this%snowcan_unload_temp_fact)
    deallocate(this%snowcan_unload_wind_fact)
    deallocate(this%interception_fraction)
    deallocate(this%maximum_leaf_wetted_fraction)

    deallocate(this%aq_sp_yield_min)
    deallocate(this%n_baseflow)
    deallocate(this%perched_baseflow_scalar)
    deallocate(this%e_ice)

    deallocate(this%fff)

    deallocate(this%slopebeta)
    deallocate(this%slopemax)
    deallocate(this%zbedrock_sf)

    deallocate(this%n_melt_coef)
    deallocate(this%accum_factor)

    deallocate(this%upplim_destruct_metamorph)

    deallocate(this%pc)
    deallocate(this%mu)

    deallocate(this%bsw_sf)
    deallocate(this%hksat_sf)
    deallocate(this%sucsat_sf)
    deallocate(this%watsat_sf)
    deallocate(this%om_frac_sf)

    deallocate(this%d_max)
    deallocate(this%frac_sat_soil_dsl_init)
    
    deallocate(this%zlnd)
    deallocate(this%snw_rds_min)

    deallocate(this%precip_repartition_nonglc_all_rain_t)
    deallocate(this%precip_repartition_nonglc_all_snow_t)

    ! deallocate(this%)


  end subroutine Clean

end module Distparamtype
