Module SoilHydrologyType

  use shr_kind_mod          , only : r8 => shr_kind_r8
  use shr_log_mod           , only : errMsg => shr_log_errMsg
  use abortutils            , only : endrun
  use decompMod             , only : bounds_type
  use clm_varpar            , only : nlevgrnd, nlayer, nlayert, nlevsoi 
  use clm_varcon            , only : spval
  use clm_varctl            , only : iulog
  use LandunitType          , only : lun                
  use ColumnType            , only : col                
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  !
  type, public :: soilhydrology_type

     integer :: h2osfcflag              ! true => surface water is active (namelist)       
     integer :: origflag                ! used to control soil hydrology properties (namelist)

     real(r8), pointer :: num_substeps_col   (:)    ! col adaptive timestep counter     
     ! NON-VIC
     real(r8), pointer :: frost_table_col   (:)     ! col frost table depth                    
     real(r8), pointer :: zwt_col           (:)     ! col water table depth
     real(r8), pointer :: zwts_col          (:)     ! col water table depth, the shallower of the two water depths
     real(r8), pointer :: zwt_perched_col   (:)     ! col perched water table depth
     real(r8), pointer :: wa_col            (:)     ! col water in the unconfined aquifer (mm)
     real(r8), pointer :: qcharge_col       (:)     ! col aquifer recharge rate (mm/s) 
     real(r8), pointer :: fracice_col       (:,:)   ! col fractional impermeability (-)
     real(r8), pointer :: icefrac_col       (:,:)   ! col fraction of ice       
     real(r8), pointer :: h2osfc_thresh_col (:)     ! col level at which h2osfc "percolates"   (time constant)
     real(r8), pointer :: xs_urban_col      (:)     ! col excess soil water above urban ponding limit

     ! VIC 
     real(r8), pointer :: hkdepth_col       (:)     ! col VIC decay factor (m) (time constant)                    
     real(r8), pointer :: b_infil_col       (:)     ! col VIC b infiltration parameter (time constant)                    
     real(r8), pointer :: ds_col            (:)     ! col VIC fracton of Dsmax where non-linear baseflow begins (time constant)                    
     real(r8), pointer :: dsmax_col         (:)     ! col VIC max. velocity of baseflow (mm/day) (time constant)
     real(r8), pointer :: Wsvic_col         (:)     ! col VIC fraction of maximum soil moisutre where non-liear base flow occurs (time constant)
     real(r8), pointer :: porosity_col      (:,:)   ! col VIC porosity (1-bulk_density/soil_density)
     real(r8), pointer :: vic_clm_fract_col (:,:,:) ! col VIC fraction of VIC layers in CLM layers 
     real(r8), pointer :: depth_col         (:,:)   ! col VIC layer depth of upper layer  
     real(r8), pointer :: c_param_col       (:)     ! col VIC baseflow exponent (Qb) 
     real(r8), pointer :: expt_col          (:,:)   ! col VIC pore-size distribution related paramter(Q12) 
     real(r8), pointer :: ksat_col          (:,:)   ! col VIC Saturated hydrologic conductivity 
     real(r8), pointer :: phi_s_col         (:,:)   ! col VIC soil moisture dissusion parameter 
     real(r8), pointer :: moist_col         (:,:)   ! col VIC soil moisture (kg/m2) for VIC soil layers 
     real(r8), pointer :: moist_vol_col     (:,:)   ! col VIC volumetric soil moisture for VIC soil layers 
     real(r8), pointer :: max_moist_col     (:,:)   ! col VIC max layer moist + ice (mm) 
     real(r8), pointer :: top_moist_col     (:)     ! col VIC soil moisture in top layers
     real(r8), pointer :: top_max_moist_col (:)     ! col VIC maximum soil moisture in top layers
     real(r8), pointer :: top_ice_col       (:)     ! col VIC ice len in top layers
     real(r8), pointer :: top_moist_limited_col(:)  ! col VIC soil moisture in top layers, limited to no greater than top_max_moist_col
     real(r8), pointer :: ice_col           (:,:)   ! col VIC soil ice (kg/m2) for VIC soil layers

   contains

     ! Public routines
     procedure, public  :: Init
     procedure, public  :: Restart

     ! Private routines
     procedure, private :: InitAllocate
     procedure, private :: InitHistory
     procedure, private :: InitCold
     procedure, private :: ReadNL

  end type soilhydrology_type

  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  !-----------------------------------------------------------------------

contains
  
  !------------------------------------------------------------------------
  subroutine Init(this, bounds, NLFilename)

    class(soilhydrology_type) :: this
    type(bounds_type), intent(in)    :: bounds  
    character(len=*), intent(in) :: NLFilename

    call this%ReadNL(NLFilename)
    call this%InitAllocate(bounds) 
    call this%InitHistory(bounds)
    call this%InitCold(bounds)

  end subroutine Init

  !------------------------------------------------------------------------
  subroutine InitAllocate(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize module data structure
    !
    ! !USES:
    use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
    !
    ! !ARGUMENTS:
    class(soilhydrology_type) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer :: begp, endp
    integer :: begc, endc
    integer :: begg, endg
    !------------------------------------------------------------------------

    begp = bounds%begp; endp= bounds%endp
    begc = bounds%begc; endc= bounds%endc
    begg = bounds%begg; endg= bounds%endg

    allocate(this%num_substeps_col   (begc:endc))                ; this%num_substeps_col   (:)     = nan
    allocate(this%frost_table_col   (begc:endc))                 ; this%frost_table_col   (:)     = nan
    allocate(this%zwt_col           (begc:endc))                 ; this%zwt_col           (:)     = nan
    allocate(this%zwt_perched_col   (begc:endc))                 ; this%zwt_perched_col   (:)     = nan
    allocate(this%zwts_col          (begc:endc))                 ; this%zwts_col          (:)     = nan

    allocate(this%wa_col            (begc:endc))                 ; this%wa_col            (:)     = nan
    allocate(this%qcharge_col       (begc:endc))                 ; this%qcharge_col       (:)     = nan
    allocate(this%fracice_col       (begc:endc,nlevgrnd))        ; this%fracice_col       (:,:)   = nan
    allocate(this%icefrac_col       (begc:endc,nlevgrnd))        ; this%icefrac_col       (:,:)   = nan
    allocate(this%h2osfc_thresh_col (begc:endc))                 ; this%h2osfc_thresh_col (:)     = nan
    allocate(this%xs_urban_col      (begc:endc))                 ; this%xs_urban_col      (:)     = nan

    allocate(this%hkdepth_col       (begc:endc))                 ; this%hkdepth_col       (:)     = nan
    allocate(this%b_infil_col       (begc:endc))                 ; this%b_infil_col       (:)     = nan
    allocate(this%ds_col            (begc:endc))                 ; this%ds_col            (:)     = nan
    allocate(this%dsmax_col         (begc:endc))                 ; this%dsmax_col         (:)     = nan
    allocate(this%Wsvic_col         (begc:endc))                 ; this%Wsvic_col         (:)     = nan
    allocate(this%depth_col         (begc:endc,nlayert))         ; this%depth_col         (:,:)   = nan
    allocate(this%porosity_col      (begc:endc,nlayer))          ; this%porosity_col      (:,:)   = nan
    allocate(this%vic_clm_fract_col (begc:endc,nlayer, nlevsoi)) ; this%vic_clm_fract_col (:,:,:) = nan
    allocate(this%c_param_col       (begc:endc))                 ; this%c_param_col       (:)     = nan
    allocate(this%expt_col          (begc:endc,nlayer))          ; this%expt_col          (:,:)   = nan
    allocate(this%ksat_col          (begc:endc,nlayer))          ; this%ksat_col          (:,:)   = nan
    allocate(this%phi_s_col         (begc:endc,nlayer))          ; this%phi_s_col         (:,:)   = nan
    allocate(this%moist_col         (begc:endc,nlayert))         ; this%moist_col         (:,:)   = nan
    allocate(this%moist_vol_col     (begc:endc,nlayert))         ; this%moist_vol_col     (:,:)   = nan
    allocate(this%max_moist_col     (begc:endc,nlayer))          ; this%max_moist_col     (:,:)   = nan
    allocate(this%top_moist_col     (begc:endc))                 ; this%top_moist_col     (:)     = nan
    allocate(this%top_max_moist_col (begc:endc))                 ; this%top_max_moist_col (:)     = nan
    allocate(this%top_ice_col       (begc:endc))                 ; this%top_ice_col       (:)     = nan
    allocate(this%top_moist_limited_col(begc:endc))              ; this%top_moist_limited_col(:)  = nan
    allocate(this%ice_col           (begc:endc,nlayert))         ; this%ice_col           (:,:)   = nan

  end subroutine InitAllocate

  !------------------------------------------------------------------------
  subroutine InitHistory(this, bounds)
    !
    ! !USES:
    use histFileMod    , only : hist_addfld1d
    !
    ! !ARGUMENTS:
    class(soilhydrology_type) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer           :: begc, endc
    integer           :: begg, endg
    !------------------------------------------------------------------------

    begc = bounds%begc; endc= bounds%endc
    begg = bounds%begg; endg= bounds%endg

    this%wa_col(begc:endc) = spval
    call hist_addfld1d (fname='WA',  units='mm',  &
         avgflag='A', long_name='water in the unconfined aquifer (vegetated landunits only)', &
         ptr_col=this%wa_col, l2g_scale_type='veg')

    this%qcharge_col(begc:endc) = spval
    call hist_addfld1d (fname='QCHARGE',  units='mm/s',  &
         avgflag='A', long_name='aquifer recharge rate (vegetated landunits only)', &
         ptr_col=this%qcharge_col, l2g_scale_type='veg')

    this%num_substeps_col(begc:endc) = spval
    call hist_addfld1d (fname='NSUBSTEPS',  units='unitless',  &
         avgflag='A', long_name='number of adaptive timesteps in CLM timestep', &
         ptr_col=this%num_substeps_col, l2g_scale_type='veg', &
         default='inactive')

    this%frost_table_col(begc:endc) = spval
    call hist_addfld1d (fname='FROST_TABLE',  units='m',  &
         avgflag='A', long_name='frost table depth (vegetated landunits only)', &
         ptr_col=this%frost_table_col, l2g_scale_type='veg', default='inactive')

    this%zwt_col(begc:endc) = spval
    call hist_addfld1d (fname='ZWT',  units='m',  &
         avgflag='A', long_name='water table depth (vegetated landunits only)', &
         ptr_col=this%zwt_col, l2g_scale_type='veg')

    this%zwt_perched_col(begc:endc) = spval
    call hist_addfld1d (fname='ZWT_PERCH',  units='m',  &
         avgflag='A', long_name='perched water table depth (vegetated landunits only)', &
         ptr_col=this%zwt_perched_col, l2g_scale_type='veg')

  end subroutine InitHistory

  !-----------------------------------------------------------------------
  subroutine InitCold(this, bounds)
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(soilhydrology_type) :: this
    type(bounds_type) , intent(in)    :: bounds
    ! !LOCAL VARIABLES:
    integer :: c ! indices

    !-----------------------------------------------------------------------

    ! Nothing for now

    ! needs to be initialized to spval to avoid problems when 
    ! averaging for the accum field
    do c = bounds%begc, bounds%endc
       this%num_substeps_col(c) = spval
    end do

  end subroutine InitCold

  !------------------------------------------------------------------------
  subroutine Restart(this, bounds, ncid, flag)
    ! 
    ! !USES:
    use ncdio_pio  , only : file_desc_t, ncd_io, ncd_double
    use restUtilMod
    !
    ! !ARGUMENTS:
    class(soilhydrology_type) :: this
    type(bounds_type) , intent(in)    :: bounds 
    type(file_desc_t) , intent(inout) :: ncid   ! netcdf id
    character(len=*)  , intent(in)    :: flag   ! 'read' or 'write'
    !
    ! !LOCAL VARIABLES:
    integer :: j,c ! indices
    logical :: readvar      ! determine if variable is on initial file
    !-----------------------------------------------------------------------

    call restartvar(ncid=ncid, flag=flag, varname='FROST_TABLE', xtype=ncd_double,  & 
         dim1name='column', &
         long_name='frost table depth', units='m', &
         interpinic_flag='interp', readvar=readvar, data=this%frost_table_col)
    if (flag == 'read' .and. .not. readvar) then
       this%frost_table_col(bounds%begc:bounds%endc) = col%zi(bounds%begc:bounds%endc,nlevsoi)
    end if

    call restartvar(ncid=ncid, flag=flag, varname='WA', xtype=ncd_double,  & 
         dim1name='column', &
         long_name='water in the unconfined aquifer', units='mm', &
         interpinic_flag='interp', readvar=readvar, data=this%wa_col)

    call restartvar(ncid=ncid, flag=flag, varname='ZWT', xtype=ncd_double,  & 
         dim1name='column', &
         long_name='water table depth', units='m', &
         interpinic_flag='interp', readvar=readvar, data=this%zwt_col)

    call restartvar(ncid=ncid, flag=flag, varname='ZWT_PERCH', xtype=ncd_double,  & 
         dim1name='column', &
         long_name='perched water table depth', units='m', &
         interpinic_flag='interp', readvar=readvar, data=this%zwt_perched_col)
    if (flag == 'read' .and. .not. readvar) then
       this%zwt_perched_col(bounds%begc:bounds%endc) = col%zi(bounds%begc:bounds%endc,nlevsoi)
    end if

  end subroutine Restart

   !-----------------------------------------------------------------------
   subroutine ReadNL( this, NLFilename )
     !
     ! !DESCRIPTION:
     ! Read namelist for SoilHydrology
     !
     ! !USES:
     use shr_mpi_mod    , only : shr_mpi_bcast
     use shr_log_mod    , only : errMsg => shr_log_errMsg
     use spmdMod        , only : masterproc, mpicom
     use fileutils      , only : getavu, relavu, opnfil
     use clm_nlUtilsMod , only : find_nlgroup_name
     use clm_varctl     , only : iulog 
     use abortutils     , only : endrun
     !
     ! !ARGUMENTS:
     class(soilhydrology_type) :: this
     character(len=*), intent(IN) :: NLFilename ! Namelist filename
     !
     ! !LOCAL VARIABLES:
     integer :: ierr                 ! error code
     integer :: unitn                ! unit for namelist file
     integer :: origflag=0            !use to control soil hydraulic properties
     integer :: h2osfcflag=1          !If surface water is active or not
     character(len=32) :: subname = 'SoilHydrology_readnl'  ! subroutine name
     !-----------------------------------------------------------------------

     namelist / clm_soilhydrology_inparm / h2osfcflag, origflag

     ! preset values

     origflag = 0          
     h2osfcflag = 1        

     if ( masterproc )then

        unitn = getavu()
        write(iulog,*) 'Read in clm_soilhydrology_inparm  namelist'
        call opnfil (NLFilename, unitn, 'F')
        call find_nlgroup_name(unitn, 'clm_soilhydrology_inparm', status=ierr)
        if (ierr == 0) then
           read(unitn, clm_soilhydrology_inparm, iostat=ierr)
           if (ierr /= 0) then
              call endrun(msg="ERROR reading clm_soilhydrology_inparm namelist"//errmsg(sourcefile, __LINE__))
           end if
        else
           call endrun(msg="ERROR finding clm_soilhydrology_inparm namelist"//errmsg(sourcefile, __LINE__))
        end if
        call relavu( unitn )

     end if

     call shr_mpi_bcast(h2osfcflag, mpicom)
     call shr_mpi_bcast(origflag,   mpicom)

     this%h2osfcflag = h2osfcflag
     this%origflag   = origflag

   end subroutine ReadNL

end Module SoilHydrologyType
