module TopoMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Handles topographic height of each column
  !
  ! !USES:
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
  use decompMod      , only : bounds_type
  use PatchType      , only : patch
  use ColumnType     , only : col
  use LandunitType   , only : lun
  use glc2lndMod     , only : glc2lnd_type
  use glcBehaviorMod , only : glc_behavior_type
  use landunit_varcon, only : istice_mec, istsoil
  use filterColMod   , only : filter_col_type, col_filter_from_logical_array_active_only
  use clm_varctl     , only : use_hillslope
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  private

  type, public :: topo_type
     private

     ! Public member data

     real(r8), pointer, public :: topo_col(:)  ! surface elevation (m)

     ! Private member data

     logical, pointer :: needs_downscaling_col(:)  ! whether a column needs to be downscaled
   contains
     procedure, public :: Init
     procedure, public :: Restart
     procedure, public :: Clean
     procedure, public :: UpdateTopo  ! Update topographic height each time step
     procedure, public :: DownscaleFilterc  ! Returns column-level filter: which columns need downscaling

     procedure, private :: InitAllocate
     procedure, private :: InitHistory
     procedure, private :: InitCold
  end type topo_type

contains

  !-----------------------------------------------------------------------
  subroutine Init(this, bounds)
    ! !ARGUMENTS:
    class(topo_type), intent(inout) :: this
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
    ! !ARGUMENTS:
    class(topo_type), intent(inout) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: begc, endc

    character(len=*), parameter :: subname = 'InitAllocate'
    !-----------------------------------------------------------------------

    begc = bounds%begc
    endc = bounds%endc

    allocate(this%topo_col(begc:endc))
    this%topo_col(:) = nan

    allocate(this%needs_downscaling_col(begc:endc))
    this%needs_downscaling_col(:) = .false.

  end subroutine InitAllocate

  !-----------------------------------------------------------------------
  subroutine InitHistory(this, bounds)
    ! !USES:
    use histFileMod  , only : hist_addfld1d
    !
    ! !ARGUMENTS:
    class(topo_type), intent(inout) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'InitHistory'
    !-----------------------------------------------------------------------

    call hist_addfld1d(fname='TOPO_COL', units='m', &
         avgflag='A', long_name='column-level topographic height', &
         ptr_col=this%topo_col, default='inactive')

    call hist_addfld1d(fname='TOPO_COL_ICE', units='m', &
         avgflag='A', long_name='column-level topographic height (ice landunits only)', &
         ptr_col=this%topo_col, l2g_scale_type='ice', default='inactive')

  end subroutine InitHistory

  !-----------------------------------------------------------------------
  subroutine InitCold(this, bounds)
    ! !USES:
    use column_varcon    , only: col_itype_to_icemec_class
    use clm_instur, only : topo_glc_mec
    ! !ARGUMENTS:
    class(topo_type), intent(inout) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: c, l, g
    integer :: icemec_class            ! current icemec class (1..maxpatch_glcmec)

    character(len=*), parameter :: subname = 'InitCold'
    !-----------------------------------------------------------------------

    do c = bounds%begc, bounds%endc
       l = col%landunit(c)
       g = col%gridcell(c)

       if (lun%itype(l) == istice_mec) then
          ! For ice_mec landunits, initialize topo_col based on surface dataset; this
          ! will get overwritten in the run loop by values sent from CISM
          icemec_class = col_itype_to_icemec_class(col%itype(c))
          this%topo_col(c) = topo_glc_mec(g, icemec_class)
          this%needs_downscaling_col(c) = .true.
       else
          ! For other landunits, arbitrarily initialize topo_col to 0 m; for landunits
          ! where this matters, this will get overwritten in the run loop by values sent
          ! from CISM
          if (lun%itype(l) == istsoil .and. use_hillslope) then
             this%topo_col(c) = col%hill_elev(c)
             this%needs_downscaling_col(c) = .true.
          else
             this%topo_col(c) = 0._r8
             this%needs_downscaling_col(c) = .false.
          endif

       end if
    end do

  end subroutine InitCold

  !-----------------------------------------------------------------------
  subroutine Restart(this, bounds, ncid, flag)
    ! !USES:
    use ncdio_pio, only : file_desc_t, ncd_double
    use restUtilMod
    !
    ! !ARGUMENTS:
    class(topo_type), intent(inout) :: this
    type(bounds_type), intent(in)    :: bounds 
    type(file_desc_t), intent(inout) :: ncid   ! netcdf id
    character(len=*) , intent(in)    :: flag   ! 'read', 'write' or 'define'
    !
    ! !LOCAL VARIABLES:
    integer :: p, c
    real(r8), pointer :: rparr(:)
    logical :: readvar

    character(len=*), parameter :: subname = 'Restart'
    !-----------------------------------------------------------------------

    allocate(rparr(bounds%begp:bounds%endp))

    ! TODO(wjs, 2016-04-05) Rename these restart variables to get rid of 'glc' in their
    ! names. However, this will require some changes to init_interp, too.

    ! This one is not actually an area, but has interpinic_flag='area' because we want to
    ! interpolate it under the same conditions under which we interpolate areas.
    call restartvar(ncid=ncid, flag=flag, varname='cols1d_topoglc', xtype=ncd_double,   &
         dim1name='column',                                                             &
         long_name='mean elevation on glacier elevation classes', units='m',            &
         interpinic_flag='area', readvar=readvar, data=this%topo_col)

    if (flag /= 'read') then
       do p=bounds%begp,bounds%endp
          c = patch%column(p)
          rparr(p) = this%topo_col(c)
       enddo
       ! This one has interpinic_flag = 'skip' because it isn't read back in
       call restartvar(ncid=ncid, flag=flag, varname='pfts1d_topoglc', xtype=ncd_double,   &
            dim1name='pft',                                                             &
            long_name='mean elevation on glacier elevation classes', units='m',            &
            interpinic_flag='skip', readvar=readvar, data=rparr)
    end if

    deallocate(rparr)

  end subroutine Restart


  !-----------------------------------------------------------------------
  subroutine UpdateTopo(this, bounds, num_icemecc, filter_icemecc, &
       glc2lnd_inst, glc_behavior, atm_topo)
    !
    ! !DESCRIPTION:
    ! Update topographic heights
    !
    ! Should be called each time step.
    !
    ! Should be called after glc2lndMod:update_glc2lnd_fracs, and before
    ! atm2lndMod:downscale_forcings
    !
    ! !ARGUMENTS:
    class(topo_type)        , intent(inout) :: this
    type(bounds_type)       , intent(in)    :: bounds
    integer                 , intent(in)    :: num_icemecc       ! number of points in filter_icemecc
    integer                 , intent(in)    :: filter_icemecc(:) ! col filter for ice_mec
    type(glc2lnd_type)      , intent(in)    :: glc2lnd_inst
    type(glc_behavior_type) , intent(in)    :: glc_behavior
    real(r8)                , intent(in)    :: atm_topo( bounds%begg: ) ! atmosphere topographic height [m]
    !
    ! !LOCAL VARIABLES:
    integer :: begc, endc
    integer :: c, l, g
    real(r8), allocatable :: mean_hillslope_elevation(:)
    real(r8):: mhe_norm

    character(len=*), parameter :: subname = 'UpdateTopo'
    !-----------------------------------------------------------------------

    begc = bounds%begc
    endc = bounds%endc

    ! Reset needs_downscaling_col each time step, because this is potentially
    ! time-varying for some columns. It's simplest just to reset it everywhere, rather
    ! than trying to figure out where it does and does not need to be reset.
    this%needs_downscaling_col(begc:endc) = .false.

    call glc_behavior%icemec_cols_need_downscaling(bounds, num_icemecc, filter_icemecc, &
         this%needs_downscaling_col(begc:endc))

    ! In addition to updating topo_col, this also sets some additional elements of
    ! needs_downscaling_col to .true. (but leaves the already-.true. values as is.)
    call glc2lnd_inst%update_glc2lnd_topo(bounds, &
         this%topo_col(begc:endc), &
         this%needs_downscaling_col(begc:endc))

    ! calculate area-weighted mean hillslope elevation on each landunit
    if (use_hillslope) then
       allocate(mean_hillslope_elevation(bounds%begl:bounds%endl))
       mean_hillslope_elevation(:) = 0._r8
       do l = bounds%begl, bounds%endl
          if (lun%itype(l) == istsoil) then
             mhe_norm = 0._r8
             do c = lun%coli(l), lun%colf(l)
                mean_hillslope_elevation(l) = mean_hillslope_elevation(l) &
                     + col%hill_elev(c)*col%hill_area(c)
                mhe_norm = mhe_norm + col%hill_area(c)
             enddo
             if (mhe_norm > 0) then
                mean_hillslope_elevation(l) = mean_hillslope_elevation(l)/mhe_norm
             endif
          endif
       enddo
    endif
       
    ! For any point that isn't downscaled, set its topo value to the atmosphere's
    ! topographic height. This shouldn't matter, but is useful if topo_col is written to
    ! the history file.
    !
    ! This could operate over a filter like 'allc' in order to just operate over active
    ! points, but I'm not sure that would speed things up much, and would require passing
    ! in this additional filter.
    do c = bounds%begc, bounds%endc
       if (.not. this%needs_downscaling_col(c)) then
          g = col%gridcell(c)
          l = col%landunit(c)
          if (lun%itype(l) == istsoil .and. use_hillslope) then
             this%topo_col(c) = atm_topo(g) &
                  + (col%hill_elev(c) - mean_hillslope_elevation(l))
             this%needs_downscaling_col(c) = .true.
          else
             this%topo_col(c) = atm_topo(g)
          endif
       end if
    end do

    call glc_behavior%update_glc_classes(bounds, this%topo_col(begc:endc))

  end subroutine UpdateTopo

  !-----------------------------------------------------------------------
  function DownscaleFilterc(this, bounds) result(filter)
    !
    ! !DESCRIPTION:
    ! Returns a column-level filter: which columns need downscaling.
    !
    ! This filter only contains active points.
    !
    ! The main reason it's important to have this filter (as opposed to just doing the
    ! downscaling for all columns) is because of downscaled fields that are normalized
    ! (like longwave radiation): Consider a gridcell with a glc_mec column and a
    ! vegetated column (outside of the icemask, so the vegetated column doesn't have its
    ! topographic height explicitly set). If we called the downscaling code for all
    ! columns, the longwave radiation would get adjusted over the vegetated column. This
    ! is undesirable, because it means that adding a downscaled column in a gridcell can
    ! change answers for all other columns in that gridcell.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    type(filter_col_type) :: filter  ! function result
    class(topo_type), intent(in) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'DownscaleFilterc'
    !-----------------------------------------------------------------------

    ! Currently this creates the filter on the fly, recreating it every time this function
    ! is called. In principle, we should be able to compute and save this filter when
    ! UpdateTopo is called, returning the already-computed filter when this function is
    ! called. However, the problem with that is the need to have a different filter for
    ! each clump (and potentially another filter for calls from outside a clump
    ! loop). This will become easier to handle if we rework CLM's threading so that there
    ! is a separate instance of each object for each clump: in that case, we'll have
    ! multiple instances of topo_type, each corresponding to one clump, each with its own
    ! filter.

    filter = col_filter_from_logical_array_active_only(bounds, &
         this%needs_downscaling_col(bounds%begc:bounds%endc))

  end function DownscaleFilterc


  !-----------------------------------------------------------------------
  subroutine Clean(this)
    ! !ARGUMENTS:
    class(topo_type), intent(inout) :: this
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'Clean'
    !-----------------------------------------------------------------------

    deallocate(this%topo_col)
    deallocate(this%needs_downscaling_col)

  end subroutine Clean

end module TopoMod
