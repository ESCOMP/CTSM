module lnd2glcMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Handle arrays used for exchanging data from land model to glc
  ! For now glc datais send and received on the lnd grid and decomposition.
  !
  ! The fields sent from the lnd component to the glc component via
  !  the coupler are labeled 's2x', or sno to coupler.
  ! The fields received by the lnd component from the glc component
  !  via the coupler are labeled 'x2s', or coupler to sno.
  ! 'Sno' is a misnomer in that the exchanged data are related to
  !  the ice beneath the snow, not the snow itself.  But by CESM convention,
  ! 'ice' refers to sea ice, not land ice.
  !
  ! !USES:
  use shr_kind_mod    , only : r8 => shr_kind_r8
  use shr_infnan_mod  , only : nan => shr_infnan_nan, assignment(=)
  use shr_log_mod     , only : errMsg => shr_log_errMsg
  use decompMod       , only : get_proc_bounds, bounds_type, subgrid_level_column
  use domainMod       , only : ldomain
  use clm_varpar      , only : maxpatch_glc
  use clm_varctl      , only : iulog, use_hillslope
  use clm_varcon      , only : spval, tfrz
  use column_varcon   , only : col_itype_to_ice_class
  use landunit_varcon , only : istice, istsoil
  use abortutils      , only : endrun
  use TemperatureType , only : temperature_type
  use WaterFluxBulkType   , only : waterfluxbulk_type
  use LandunitType    , only : lun                
  use ColumnType      , only : col
  use TopoMod         , only : topo_type
  !
  ! !PUBLIC TYPES:
  implicit none
  private
  save

  ! land -> glc variables structure
  type, public :: lnd2glc_type
     real(r8), pointer :: tsrf_grc(:,:) => null()
     real(r8), pointer :: topo_grc(:,:) => null()
     real(r8), pointer :: qice_grc(:,:) => null()

   contains

     procedure, public  :: Init
     procedure, public  :: update_lnd2glc
     procedure, private :: InitAllocate
     procedure, private :: InitHistory

  end type lnd2glc_type

  ! !PUBLIC MEMBER FUNCTIONS:
  
  ! The following is public simply to support unit testing, and should not generally be
  ! called from outside this module.
  !
  ! Note that it is not a type-bound procedure, because it doesn't actually involve the
  ! lnd2glc_type. This suggests that perhaps it belongs in some other module.
  public :: bareland_normalization ! compute normalization factor for fluxes from the bare land portion of the grid cell

  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  !------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  subroutine Init(this, bounds)

    class(lnd2glc_type) :: this
    type(bounds_type), intent(in) :: bounds  

    call this%InitAllocate(bounds)
    call this%InitHistory(bounds)
    
  end subroutine Init

  !------------------------------------------------------------------------
  subroutine InitAllocate(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize land variables required by glc
    !
    ! !USES:
    use clm_varcon , only : spval
    use histFileMod, only : hist_addfld1d
    !
    ! !ARGUMENTS:
    class(lnd2glc_type) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer :: begg,endg 
    !------------------------------------------------------------------------

    begg = bounds%begg; endg = bounds%endg

    allocate(this%tsrf_grc(begg:endg,0:maxpatch_glc)) ; this%tsrf_grc(:,:)=0.0_r8
    allocate(this%topo_grc(begg:endg,0:maxpatch_glc)) ; this%topo_grc(:,:)=0.0_r8
    allocate(this%qice_grc(begg:endg,0:maxpatch_glc)) ; this%qice_grc(:,:)=0.0_r8

  end subroutine InitAllocate

  !------------------------------------------------------------------------
  subroutine InitHistory(this, bounds)
    !
    ! !USES:
    use histFileMod, only : hist_addfld1d,hist_addfld2d 
    !
    ! !ARGUMENTS:
    class(lnd2glc_type) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    real(r8), pointer :: data2dptr(:,:)
    integer  :: begg, endg
    !---------------------------------------------------------------------

    begg = bounds%begg; endg = bounds%endg

    this%qice_grc(begg:endg,0:maxpatch_glc) = spval
    ! For this and the following fields, set up a pointer to the field simply for the
    ! sake of changing the indexing, so that levels start with an index of 1, as is
    ! assumed by histFileMod - so levels go 1:(nec+1) rather than 0:nec
    data2dptr => this%qice_grc(:,0:maxpatch_glc)
    call hist_addfld2d (fname='QICE_FORC', units='mm/s', type2d='elevclas', &
         avgflag='A', long_name='qice forcing sent to GLC', &
         ptr_lnd=data2dptr, default='inactive')

    this%tsrf_grc(begg:endg,0:maxpatch_glc) = spval
    data2dptr => this%tsrf_grc(:,0:maxpatch_glc)
    call hist_addfld2d (fname='TSRF_FORC', units='K', type2d='elevclas', &
         avgflag='A', long_name='surface temperature sent to GLC', &
         ptr_lnd=data2dptr, default='inactive')

    this%topo_grc(begg:endg,0:maxpatch_glc) = spval
    data2dptr => this%topo_grc(:,0:maxpatch_glc)
    call hist_addfld2d (fname='TOPO_FORC', units='m', type2d='elevclas', &
         avgflag='A', long_name='topograephic height sent to GLC', &
         ptr_lnd=data2dptr, default='inactive')

  end subroutine InitHistory


  !------------------------------------------------------------------------------
  subroutine update_lnd2glc(this, bounds, num_do_smb_c, filter_do_smb_c, &
       temperature_inst, waterfluxbulk_inst, topo_inst, init)
    !
    ! !DESCRIPTION:
    ! Assign values to lnd2glc+
    !
    ! !ARGUMENTS:
    class(lnd2glc_type)    , intent(inout) :: this
    type(bounds_type)      , intent(in)    :: bounds  
    integer                , intent(in)    :: num_do_smb_c       ! number of columns in filter_do_smb_c
    integer                , intent(in)    :: filter_do_smb_c(:) ! column filter: columns where smb calculations are performed
    type(temperature_type) , intent(in)    :: temperature_inst
    type(waterfluxbulk_type)   , intent(in)    :: waterfluxbulk_inst
    type(topo_type)        , intent(in)    :: topo_inst
    logical                , intent(in)    :: init               ! if true=>only set a subset of fields
    !
    ! !LOCAL VARIABLES:
    integer  :: c, l, g, n, fc                   ! indices
    logical, allocatable :: fields_assigned(:,:) ! tracks whether fields have already been assigned for each index [begg:endg, 0:maxpatch_glc]
    real(r8) :: flux_normalization               ! factor by which fluxes should be normalized

    character(len=*), parameter :: subname = 'update_lnd2glc'
    !------------------------------------------------------------------------------

    ! Initialize to reasonable defaults. These values will be sent for gridcells /
    ! columns outside the do_smb filter.

    ! NOTE(wjs, 2018-07-03) qice should be 0 outside the do_smb filter to ensure conservation
    this%qice_grc(bounds%begg : bounds%endg, :) = 0._r8

    ! NOTE(wjs, 2018-07-03) tsrf can be anything outside the do_smb filter; 0 deg C seems
    ! as reasonable as anything (based on input from Bill Lipscomb and Gunter Leguy)
    this%tsrf_grc(bounds%begg : bounds%endg, :) = tfrz

    ! NOTE(wjs, 2018-07-03) The topo values outside the do_smb filter could matter for
    ! gridcells where we compute SMB for some but not all elevation classes (possibly
    ! because we haven't even allocated memory for some elevation classes - i.e., if we're
    ! not using the 'virtual' behavior in that gridcell). In glc2lndMod, we ensure that
    ! this cannot occur for gridcells within the icemask (i.e., within the icemask, we
    ! ensure that there are no points that have (non-virtual and compute-SMB)), so this
    ! isn't a conservation issue, but it could still be important, e.g., for generating
    ! appropriate forcings for a later dlnd-driven T compset. I'm not sure what is "right"
    ! here. We've historically used 0 for this, and maybe that's as good as anything,
    ! because it may lead to the 0 SMB values being ignored for the sake of vertical
    ! interpolation, but I'm not sure about this. But maybe it would be better to use
    ! topo at the center of each elevation class?
    this%topo_grc(bounds%begg : bounds%endg, :) = 0._r8     
  
    ! Fill the lnd->glc data on the clm grid

    allocate(fields_assigned(bounds%begg:bounds%endg, 0:maxpatch_glc))
    fields_assigned(:,:) = .false.

    do fc = 1, num_do_smb_c
      c = filter_do_smb_c(fc)
      l = col%landunit(c)
      g = col%gridcell(c) 

      ! Set vertical index and a flux normalization, based on whether the column in question is glacier or vegetated.  
      if (lun%itype(l) == istice) then
         n = col_itype_to_ice_class(col%itype(c))
         flux_normalization = 1.0_r8
      else if (lun%itype(l) == istsoil) then
         n = 0  !0-level index (bareland information)
         flux_normalization = bareland_normalization(c)
      else
         ! Other landunit types do not pass information in the lnd2glc fields.
         ! Note: for this to be acceptable, we need virtual vegetated columns in any grid
         ! cell that is made up solely of glacier plus some other special landunit (e.g.,
         ! glacier + lake) -- otherwise CISM wouldn't have any information for the non-
         ! glaciated portion of the grid cell.
         cycle
      end if

      ! Make sure we haven't already assigned the coupling fields for this point
      ! (this could happen, for example, if there were multiple columns in the
      ! istsoil landunit, which we aren't prepared to handle)
      !
      ! BUG(wjs, 2022-07-17, ESCOMP/CTSM#204) We have a known bug in the handling of bare
      ! land fluxes when we potentially have multiple vegetated columns in a grid cell.
      ! The most common configuration where this is the case is when use_hillslope is
      ! true. In order to allow hillslope hydrology runs to work for now, we are
      ! bypassing this error check when use_hillslope is true - under the assumption
      ! that, for now, people aren't going to be interested in SMB in a run with
      ! hillslope hydrology. Once we resolve ESCOMP/CTSM#204, we should remove the '.and.
      ! .not. use_hillslope' part of this conditional.
      if (fields_assigned(g,n) .and. .not. use_hillslope) then
         write(iulog,*) subname//' ERROR: attempt to assign coupling fields twice for the same index.'
         write(iulog,*) 'One possible cause is having multiple columns in the istsoil landunit,'
         write(iulog,*) 'which this routine cannot handle.'
         write(iulog,*) 'g, n = ', g, n
         call endrun(subgrid_index=c, subgrid_level=subgrid_level_column, msg=errMsg(sourcefile, __LINE__))
      end if

      ! Send surface temperature, topography, and SMB flux (qice) to coupler.
      ! t_soisno and topo_col are valid even in initialization, so tsrf and topo
      ! are set here regardless of the value of init. But qflx_glcice is not valid
      ! until the run loop; thus, in initialization, we will use the default value
      ! for qice, as set above.
      fields_assigned(g,n) = .true.
      this%tsrf_grc(g,n) = temperature_inst%t_soisno_col(c,1)
      this%topo_grc(g,n) = topo_inst%topo_col(c)
      if (.not. init) then
         this%qice_grc(g,n) = waterfluxbulk_inst%qflx_glcice_col(c) * flux_normalization

         ! Check for bad values of qice
         if ( abs(this%qice_grc(g,n)) > 1.0_r8) then
            write(iulog,*) 'WARNING: qice out of bounds: g, n, qice =', g, n, this%qice_grc(g,n)
         end if
      end if

    end do

    deallocate(fields_assigned)
                
  end subroutine update_lnd2glc

  !-----------------------------------------------------------------------
  real(r8) function bareland_normalization(c)
    !
    ! !DESCRIPTION:
    ! Compute normalization factor for fluxes from the bare land portion of the grid
    ! cell. Fluxes should be multiplied by this factor before being sent to CISM.
    !
    ! The point of this is: CISM effectively has two land cover types: glaciated and
    ! bare. CLM, on the other hand, subdivides the bare land portion of the grid cell into
    ! multiple landunits. However, we currently don't do any sort of averaging of
    ! quantities computed in the different "bare land" landunits - instead, we simply send
    ! the values computed in the natural vegetated landunit - these fluxes (like SMB) are
    ! 0 in the other landunits. To achieve conservation, we need to normalize these
    ! natural veg. fluxes by the fraction of the "bare land" area accounted for by the
    ! natural veg. landunit.
    !
    ! For example, consider a grid cell that is:
    !   60% glacier_mec
    !   30% natural veg
    !   10% lake
    !
    ! According to CISM, this grid cell is 60% icesheet, 40% "bare land". Now suppose CLM
    ! has an SMB flux of 1m in the natural veg landunit. If we simply sent 1m of ice to
    ! CISM, conservation would be broken, since it would also apply 1m of ice to the 10%
    ! of the grid cell that CLM says is lake. So, instead, we must multiply the 1m of ice
    ! by (0.3/0.4), thus "spreading out" the SMB from the natural veg. landunit, so that
    ! 0.75m of ice is grown throughout the bare land portion of CISM.
    !
    ! Note: If the non-glaciated area of the grid cell is 0, then we arbitrarily return a
    ! normalization factor of 1.0, in order to avoid divide-by-zero errors.
    !
    ! Note: We currently aren't careful about how we would handle things if there are
    ! multiple columns within the vegetated landunit. If that possibility were introduced,
    ! this code - as well as the code in update_clm_s2x - may need to be reworked somewhat.
    !
    ! !USES:
    use subgridWeightsMod , only : get_landunit_weight
    !
    ! !ARGUMENTS:
    integer, intent(in) :: c  ! column index
    !
    ! !LOCAL VARIABLES:
    integer  :: g             ! grid cell index
    real(r8) :: area_glacier  ! fractional area of the glacier_mec landunit in this grid cell
    real(r8) :: area_this_col ! fractional area of column c in the grid cell

    real(r8), parameter :: tol = 1.e-13_r8  ! tolerance for checking subgrid weight equality
    character(len=*), parameter :: subname = 'bareland_normalization'
    !-----------------------------------------------------------------------

    g = col%gridcell(c)

    area_glacier = get_landunit_weight(g, istice)

    if (abs(area_glacier - 1.0_r8) < tol) then
       ! If the whole grid cell is glacier, then the normalization factor is arbitrary;
       ! set it to 1 so we don't do any normalization in this case
       bareland_normalization = 1.0_r8
    else
       area_this_col = col%wtgcell(c)
       bareland_normalization = area_this_col / (1.0_r8 - area_glacier)
    end if

  end function bareland_normalization

end module lnd2glcMod

