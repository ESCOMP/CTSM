module ActiveLayerMod
  
  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module holding routines for calculation of active layer dynamics
  !
  ! !USES:
  use shr_kind_mod    , only : r8 => shr_kind_r8
  use decompMod       , only : bounds_type
  use shr_const_mod   , only : SHR_CONST_TKFRZ
  use clm_varctl      , only : iulog, use_cn
  use clm_varcon      , only : spval  
  use landunit_varcon , only : istsoil, istcrop
  use TemperatureType , only : temperature_type
  use GridcellType    , only : grc
  use LandunitType    , only : lun
  use ColumnType      , only : col
  !
  implicit none
  save
  private
  !
  ! !PUBLIC TYPES:
  type, public :: active_layer_type
     private
     ! Public data members:
     real(r8) , pointer, public :: altmax_col               (:)   ! col maximum annual depth of thaw
     real(r8) , pointer, public :: altmax_lastyear_col      (:)   ! col prior year maximum annual depth of thaw
     integer  , pointer, public :: altmax_indx_col          (:)   ! col maximum annual depth of thaw
     integer  , pointer, public :: altmax_lastyear_indx_col (:)   ! col prior year maximum annual depth of thaw

     ! Private data members:
     real(r8) , pointer :: alt_col                  (:)   ! col current depth of thaw
     integer  , pointer :: alt_indx_col             (:)   ! col current depth of thaw
   contains
     ! Public routines
     procedure, public :: alt_calc
     procedure, public :: Init
     procedure, public :: Restart

     ! Private routines
     procedure, private :: InitAllocate
     procedure, private :: InitHistory
     procedure, private :: InitCold
  end type active_layer_type

  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public:: alt_calc
  !-----------------------------------------------------------------------
  
contains

  ! ========================================================================
  ! Science routines
  ! ========================================================================

  !-----------------------------------------------------------------------
  subroutine alt_calc(this, num_soilc, filter_soilc, &
       temperature_inst)
    !
    ! !DESCRIPTION:
    !  define active layer thickness similarly to frost_table, except set as deepest thawed layer and define on nlevgrnd
    !  also update annual maxima, and keep track of prior year for rooting memory
    !
    ! BUG(wjs, 2014-12-15, bugz 2107) Because of this routine's placement in the driver
    ! sequence (it is called very early in each timestep, before weights are adjusted and
    ! filters are updated), it may be necessary for this routine to compute values over
    ! inactive as well as active points (since some inactive points may soon become
    ! active) - so that's what is done now. Currently, it seems to be okay to do this,
    ! because the variables computed here seem to only depend on quantities that are valid
    ! over inactive as well as active points.
    !
    ! !USES:
    use shr_const_mod    , only : SHR_CONST_TKFRZ
    use clm_varpar       , only : nlevgrnd
    use clm_time_manager , only : get_curr_date, get_step_size
    use clm_varctl       , only : iulog
    use clm_varcon       , only : zsoi
    !
    ! !ARGUMENTS:
    class(active_layer_type), intent(inout) :: this
    integer                , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                , intent(in)    :: filter_soilc(:) ! filter for soil columns
    type(temperature_type) , intent(in)    :: temperature_inst
    !
    ! !LOCAL VARIABLES:
    integer  :: c, j, fc, g                     ! counters
    integer  :: alt_ind                         ! index of base of activel layer
    integer  :: year                            ! year (0, ...) for nstep+1
    integer  :: mon                             ! month (1, ..., 12) for nstep+1
    integer  :: day                             ! day of month (1, ..., 31) for nstep+1
    integer  :: sec                             ! seconds into current date for nstep+1
    integer  :: dtime                           ! time step length in seconds
    integer  :: k_frz                           ! index of first nonfrozen soil layer
    logical  :: found_thawlayer                 ! used to break loop when first unfrozen layer reached
    real(r8) :: t1, t2, z1, z2                  ! temporary variables
    !-----------------------------------------------------------------------

    associate(                                                                & 
         t_soisno             =>    temperature_inst%t_soisno_col , & ! Input:   [real(r8) (:,:) ]  soil temperature (Kelvin)  (-nlevsno+1:nlevgrnd)

         alt                  =>    this%alt_col                  , & ! Output:  [real(r8) (:)   ]  current depth of thaw
         altmax               =>    this%altmax_col               , & ! Output:  [real(r8) (:)   ]  maximum annual depth of thaw
         altmax_lastyear      =>    this%altmax_lastyear_col      , & ! Output:  [real(r8) (:)   ]  prior year maximum annual depth of thaw
         alt_indx             =>    this%alt_indx_col             , & ! Output:  [integer  (:)   ]  current depth of thaw
         altmax_indx          =>    this%altmax_indx_col          , & ! Output:  [integer  (:)   ]  maximum annual depth of thaw
         altmax_lastyear_indx =>    this%altmax_lastyear_indx_col   & ! Output:  [integer  (:)   ]  prior year maximum annual depth of thaw
         )

      ! on a set annual timestep, update annual maxima
      ! make this 1 January for NH columns, 1 July for SH columns
      call get_curr_date(year, mon, day, sec)
      dtime =  get_step_size()
      if ( (mon .eq. 1) .and. (day .eq. 1) .and. ( sec / dtime .eq. 1) ) then
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            g = col%gridcell(c)
            if ( grc%lat(g) > 0._r8 ) then 
               altmax_lastyear(c) = altmax(c)
               altmax_lastyear_indx(c) = altmax_indx(c)
               altmax(c) = 0._r8
               altmax_indx(c) = 0
            endif
         end do
      endif
      if ( (mon .eq. 7) .and. (day .eq. 1) .and. ( sec / dtime .eq. 1) ) then
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            g = col%gridcell(c)
            if ( grc%lat(g) <= 0._r8 ) then 
               altmax_lastyear(c) = altmax(c)
               altmax_lastyear_indx(c) = altmax_indx(c)
               altmax(c) = 0._r8
               altmax_indx(c) = 0
            endif
         end do
      endif

      do fc = 1,num_soilc
         c = filter_soilc(fc)

         ! calculate alt for a given timestep
         ! start from base of soil and search upwards for first thawed layer.
         ! note that this will put talik in with active layer
         ! a different way of doing this could be to keep track of how long a given layer has ben frozen for, and define ALT as the first layer that has been frozen for less than 2 years.
         if (t_soisno(c,nlevgrnd) > SHR_CONST_TKFRZ ) then
            alt(c) = zsoi(nlevgrnd)
            alt_indx(c) = nlevgrnd
         else
            k_frz=0
            found_thawlayer = .false.
            do j=nlevgrnd-1,1,-1
               if ( ( t_soisno(c,j) > SHR_CONST_TKFRZ ) .and. .not. found_thawlayer ) then
                  k_frz=j
                  found_thawlayer = .true.
               endif
            end do

            if ( k_frz > 0 ) then
               ! define active layer as the depth at which the linearly interpolated temperature line intersects with zero
               z1 = zsoi(k_frz)
               z2 = zsoi(k_frz+1)
               t1 = t_soisno(c,k_frz)
               t2 = t_soisno(c,k_frz+1)
               alt(c) = z1 + (t1-SHR_CONST_TKFRZ)*(z2-z1)/(t1-t2)
               alt_indx(c) = k_frz
            else
               alt(c)=0._r8
               alt_indx(c) = 0
            endif
         endif


         ! if appropriate, update maximum annual active layer thickness
         if (alt(c) > altmax(c)) then
            altmax(c) = alt(c)
            altmax_indx(c) = alt_indx(c)
         endif

      end do

    end associate 

  end subroutine alt_calc

  ! ========================================================================
  ! Infrastructure routines (for initialization & restart)
  ! ========================================================================

  !-----------------------------------------------------------------------
  subroutine Init(this, bounds)
    !
    ! !ARGUMENTS:
    class(active_layer_type), intent(inout) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'Init'
    !-----------------------------------------------------------------------

    call this%InitAllocate(bounds)
    call this%InitHistory(bounds)
    call this%InitCold(bounds)

  end subroutine Init

  !-----------------------------------------------------------------------
  subroutine InitAllocate(this, bounds)
    !
    ! !ARGUMENTS:
    class(active_layer_type), intent(inout) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'InitAllocate'
    !-----------------------------------------------------------------------

    associate( &
         begc => bounds%begc, &
         endc => bounds%endc  &
         )

    allocate(this%alt_col                  (begc:endc))           ; this%alt_col                  (:)   = spval     
    allocate(this%altmax_col               (begc:endc))           ; this%altmax_col               (:)   = spval
    allocate(this%altmax_lastyear_col      (begc:endc))           ; this%altmax_lastyear_col      (:)   = spval
    allocate(this%alt_indx_col             (begc:endc))           ; this%alt_indx_col             (:)   = huge(1)
    allocate(this%altmax_indx_col          (begc:endc))           ; this%altmax_indx_col          (:)   = huge(1)
    allocate(this%altmax_lastyear_indx_col (begc:endc))           ; this%altmax_lastyear_indx_col (:)   = huge(1)

    end associate

  end subroutine InitAllocate

  !-----------------------------------------------------------------------
  subroutine InitHistory(this, bounds)
    !
    ! !USES:
    use histFileMod   , only: hist_addfld1d
    !
    ! !ARGUMENTS:
    class(active_layer_type), intent(inout) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    character(len=:), allocatable :: active_if_cn

    character(len=*), parameter :: subname = 'InitHistory'
    !-----------------------------------------------------------------------

    associate( &
         begc => bounds%begc, &
         endc => bounds%endc  &
         )

    if (use_cn) then
       active_if_cn = 'active'
    else
       active_if_cn = 'inactive'
    end if

    this%alt_col(begc:endc) = spval
    call hist_addfld1d (fname='ALT', units='m', &
         avgflag='A', long_name='current active layer thickness', &
         ptr_col=this%alt_col, default=active_if_cn)

    this%altmax_col(begc:endc) = spval
    call hist_addfld1d (fname='ALTMAX', units='m', &
         avgflag='A', long_name='maximum annual active layer thickness', &
         ptr_col=this%altmax_col, default=active_if_cn)

    this%altmax_lastyear_col(begc:endc) = spval
    call hist_addfld1d (fname='ALTMAX_LASTYEAR', units='m', &
         avgflag='A', long_name='maximum prior year active layer thickness', &
         ptr_col=this%altmax_lastyear_col, default='inactive')

    end associate

  end subroutine InitHistory

  !-----------------------------------------------------------------------
  subroutine InitCold(this, bounds)
    !
    ! !ARGUMENTS:
    class(active_layer_type), intent(inout) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: l, c

    character(len=*), parameter :: subname = 'InitCold'
    !-----------------------------------------------------------------------

    do c = bounds%begc, bounds%endc
       l = col%landunit(c)

       if (lun%itype(l) == istsoil .or. lun%itype(l) == istcrop) then
          this%alt_col(c)               = 0._r8 !iniitialized to spval for all columns
          this%altmax_col(c)            = 0._r8 !iniitialized to spval for all columns
          this%altmax_lastyear_col(c)   = 0._r8 !iniitialized to spval for all columns
          this%alt_indx_col(c)          = 0     !initiialized to huge  for all columns
          this%altmax_indx_col(c)       = 0     !initiialized to huge  for all columns
          this%altmax_lastyear_indx_col = 0     !initiialized to huge  for all columns
       end if
    end do

  end subroutine InitCold

  !-----------------------------------------------------------------------
  subroutine Restart(this, bounds, ncid, flag)
    !
    ! !USES:
    use ncdio_pio  , only : file_desc_t, ncd_double, ncd_int
    use restUtilMod
    !
    ! !ARGUMENTS:
    class(active_layer_type), intent(inout) :: this
    type(bounds_type) , intent(in)    :: bounds 
    type(file_desc_t) , intent(inout) :: ncid   ! netcdf id
    character(len=*)  , intent(in)    :: flag   ! 'read' or 'write'
    !
    ! !LOCAL VARIABLES:
    logical :: readvar      ! determine if variable is on initial file

    character(len=*), parameter :: subname = 'Restart'
    !-----------------------------------------------------------------------

    call restartvar(ncid=ncid, flag=flag, varname='altmax', xtype=ncd_double,  &
         dim1name='column', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%altmax_col) 

    call restartvar(ncid=ncid, flag=flag, varname='altmax_lastyear', xtype=ncd_double,  &
         dim1name='column', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%altmax_lastyear_col) 

    call restartvar(ncid=ncid, flag=flag, varname='altmax_indx', xtype=ncd_int,  &
         dim1name='column', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%altmax_indx_col) 

    call restartvar(ncid=ncid, flag=flag, varname='altmax_lastyear_indx', xtype=ncd_int,  &
         dim1name='column', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%altmax_lastyear_indx_col) 

  end subroutine Restart



end module ActiveLayerMod
