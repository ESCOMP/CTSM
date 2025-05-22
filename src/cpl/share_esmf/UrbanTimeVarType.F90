module UrbanTimeVarType

  !------------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Urban Time Varying Data
  !
  ! !USES:
  use ESMF            , only : ESMF_LogFoundError, ESMF_LOGERR_PASSTHRU, ESMF_Finalize, ESMF_END_ABORT
  use dshr_strdata_mod, only : shr_strdata_type
  use shr_kind_mod    , only : r8 => shr_kind_r8, CL => shr_kind_CL
  use shr_log_mod     , only : errMsg => shr_log_errMsg
  use abortutils      , only : endrun
  use decompMod       , only : bounds_type, subgrid_level_landunit
  use clm_varctl      , only : iulog, FL => fname_len
  use landunit_varcon , only : isturb_MIN, isturb_MAX
  use clm_varcon      , only : spval
  use LandunitType    , only : lun
  use GridcellType    , only : grc
  !
  implicit none
  private
  !
  ! !PUBLIC TYPE
  type, public :: urbantv_type
     !
     real(r8), public, pointer :: t_building_max(:)    ! lun maximum internal building air temperature (K)
     real(r8), public, pointer :: p_ac(:)              ! lun air-conditioning adoption rate (unitless, between 0 and 1)
     type(shr_strdata_type)    :: sdat_urbantv         ! urban time varying input data stream
   contains
     ! !PUBLIC MEMBER FUNCTIONS:
     procedure, public :: Init              ! Allocate and initialize urbantv
     procedure, public :: urbantv_init      ! Initialize urban time varying stream
     procedure, public :: urbantv_interp    ! Interpolate urban time varying stream
  end type urbantv_type
  
  integer      , private              :: stream_varname_MIN       ! minimum index for stream_varnames
  integer      , private              :: stream_varname_MAX       ! maximum index for stream_varnames
  character(15), private, pointer     :: stream_varnames(:)       ! urban time varying variable names

  character(len=*), parameter, private :: sourcefile = &
       __FILE__

!==============================================================================
contains
!==============================================================================

  subroutine Init(this, bounds, NLFilename)
    !
    ! Allocate module variables and data structures
    !
    ! !USES:
    use shr_infnan_mod  , only : nan => shr_infnan_nan, assignment(=)
    use histFileMod     , only : hist_addfld1d
    use UrbanParamsType , only : urban_explicit_ac
    !
    ! !ARGUMENTS:
    class(urbantv_type) :: this
    type(bounds_type) , intent(in) :: bounds
    character(len=*)  , intent(in) :: NLFilename   ! Namelist filename
    !
    ! !LOCAL VARIABLES:
    integer		:: begl, endl
    !---------------------------------------------------------------------

    begl = bounds%begl; endl = bounds%endl

    ! Determine the minimum and maximum indices for stream_varnames
    stream_varname_MIN = 1
    ! Get value for the maximum index for stream_varnames: if using explicit AC adoption scheme, 
    ! then set maximum index to 6 for reading in tbuildmax and p_ac for three urban density classes;
    ! otherwise, set to 3 to only read in tbuildmax for three urban density classes. 
    if (urban_explicit_ac) then
       stream_varname_MAX = 6
    else
       stream_varname_MAX = 3
    end if

    ! Allocate urbantv data structure

    allocate(this%t_building_max(begl:endl)); this%t_building_max(:) = nan
    allocate(this%p_ac(begl:endl)); this%p_ac(:) = nan
    allocate(stream_varnames(stream_varname_MIN:stream_varname_MAX))

    call this%urbantv_init(bounds, NLFilename)
    call this%urbantv_interp(bounds)

    ! Add history fields
    call hist_addfld1d (fname='TBUILD_MAX', units='K',      &
          avgflag='A', long_name='prescribed maximum interior building temperature',   &
          ptr_lunit=this%t_building_max, default='inactive', set_nourb=spval, &
          l2g_scale_type='unity')
    if (urban_explicit_ac) then
       call hist_addfld1d (fname='P_AC', units='a fraction between 0 and 1',      &
             avgflag='A', long_name='prescribed air-conditioning ownership rate',   &
             ptr_lunit=this%p_ac, default='inactive', set_nourb=spval, &
             l2g_scale_type='unity')
    end if

  end subroutine Init

  !==============================================================================
  subroutine urbantv_init(this, bounds, NLFilename)
    !
    ! !DESCRIPTION:
    ! Initialize data stream information for urban time varying data
    !
    ! !USES:
    use clm_nlUtilsMod   , only : find_nlgroup_name
    use spmdMod          , only : masterproc, mpicom, iam
    use shr_mpi_mod      , only : shr_mpi_bcast
    use landunit_varcon  , only : isturb_tbd, isturb_hd, isturb_md
    use dshr_strdata_mod , only : shr_strdata_init_from_inline
    use lnd_comp_shr     , only : mesh, model_clock
    use UrbanParamsType  , only : urban_explicit_ac
    !
    ! !ARGUMENTS:
    implicit none
    class(urbantv_type)           :: this
    type(bounds_type), intent(in) :: bounds
    character(len=*),  intent(in) :: NLFilename   ! Namelist filename
    !
    ! !LOCAL VARIABLES:
    integer            :: n
    integer            :: stream_year_first_urbantv         ! first year in urban tv stream to use
    integer            :: stream_year_last_urbantv          ! last year in urban tv stream to use
    integer            :: model_year_align_urbantv          ! align stream_year_first_urbantv with this model year
    integer            :: nu_nml                            ! unit for namelist file
    integer            :: nml_error                         ! namelist i/o error flag
    character(len=FL)  :: stream_fldFileName_urbantv        ! urban tv streams filename
    character(len=FL)  :: stream_meshfile_urbantv           ! urban tv streams filename
    character(len=CL)  :: urbantvmapalgo = 'nn'             ! mapping alogrithm for urban ac
    character(len=CL)  :: urbantv_tintalgo = 'linear'       ! time interpolation alogrithm
    integer            :: rc                                ! error code
    character(*), parameter :: subName = "('urbantv_init')"
    !-----------------------------------------------------------------------

    namelist /urbantv_streams/       &
         stream_year_first_urbantv,  &
         stream_year_last_urbantv,   &
         model_year_align_urbantv,   &
         urbantvmapalgo,             &
         stream_fldFileName_urbantv, &
         stream_meshfile_urbantv, &
         urbantv_tintalgo

    ! Default values for namelist
    stream_year_first_urbantv  = 1       ! first year in stream to use
    stream_year_last_urbantv   = 1       ! last  year in stream to use
    model_year_align_urbantv   = 1       ! align stream_year_first_urbantv with this model year
    stream_fldFileName_urbantv = ' '
    stream_meshfile_urbantv    = ' '
    stream_varnames(1) = "tbuildmax_TBD"
    stream_varnames(2) = "tbuildmax_HD"
    stream_varnames(3) = "tbuildmax_MD"
    if (urban_explicit_ac) then
       stream_varnames(4) = "p_ac_TBD"
       stream_varnames(5) = "p_ac_HD"
       stream_varnames(6) = "p_ac_MD"
    end if

    ! Read urbantv_streams namelist
    if (masterproc) then
       open( newunit=nu_nml, file=trim(NLFilename), status='old', iostat=nml_error )
       call find_nlgroup_name(nu_nml, 'urbantv_streams', status=nml_error)
       if (nml_error == 0) then
          read(nu_nml, nml=urbantv_streams,iostat=nml_error)
          if (nml_error /= 0) then
             call endrun(msg='ERROR reading urbantv_streams namelist'//errMsg(sourcefile, __LINE__))
          end if
       end if
       close(nu_nml)
    endif

    call shr_mpi_bcast(stream_year_first_urbantv  , mpicom)
    call shr_mpi_bcast(stream_year_last_urbantv   , mpicom)
    call shr_mpi_bcast(model_year_align_urbantv   , mpicom)
    call shr_mpi_bcast(stream_fldFileName_urbantv , mpicom)
    call shr_mpi_bcast(stream_meshfile_urbantv    , mpicom)
    call shr_mpi_bcast(urbantv_tintalgo           , mpicom)

    if (masterproc) then
       write(iulog,*) ' '
       write(iulog,'(a)') 'urbantv_streams settings:'
       write(iulog,'(a,i8)') '  stream_year_first_urbantv  = ',stream_year_first_urbantv
       write(iulog,'(a,i8)') '  stream_year_last_urbantv   = ',stream_year_last_urbantv
       write(iulog,'(a,i8)') '  model_year_align_urbantv   = ',model_year_align_urbantv
       write(iulog,'(a,a)' ) '  stream_fldFileName_urbantv = ',stream_fldFileName_urbantv
       write(iulog,'(a,a)' ) '  stream_meshfile_urbantv    = ',stream_meshfile_urbantv
       write(iulog,'(a,a)' ) '  urbantv_tintalgo           = ',urbantv_tintalgo
       do n = stream_varname_MIN,stream_varname_MAX
          write(iulog,'(a,a)' ) '  stream_varname         = ',trim(stream_varnames(n))
       end do
       write(iulog,*) ' '
    endif

    ! Initialize the cdeps data type this%sdat_urbantv
    call shr_strdata_init_from_inline(this%sdat_urbantv,             &
         my_task             = iam,                                  &
         logunit             = iulog,                                &
         compname            = 'LND',                                &
         model_clock         = model_clock,                          &
         model_mesh          = mesh,                                 &
         stream_meshfile     = trim(stream_meshfile_urbantv),        &
         stream_lev_dimname  = 'null',                               &
         stream_mapalgo      = trim(urbantvmapalgo),                 &
         stream_filenames    = (/trim(stream_fldfilename_urbantv)/), &
         stream_fldlistFile  = stream_varnames(stream_varname_MIN:stream_varname_MAX), &
         stream_fldListModel = stream_varnames(stream_varname_MIN:stream_varname_MAX), &
         stream_yearFirst    = stream_year_first_urbantv,            &
         stream_yearLast     = stream_year_last_urbantv,             &
         stream_yearAlign    = model_year_align_urbantv,             &
         stream_offset       = 0,                                    &
         stream_taxmode      = 'extend',                             &
         stream_dtlimit      = 1.0e30_r8,                            &
         stream_tintalgo     = urbantv_tintalgo,                     &
         stream_name         = 'Urban time varying data',            &
         rc                  = rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
       call ESMF_Finalize(endflag=ESMF_END_ABORT)
    end if

  end subroutine urbantv_init

  !==============================================================================
  subroutine urbantv_interp(this, bounds)
    !
    ! !DESCRIPTION:
    ! Interpolate data stream information for urban time varying data.
    !
    ! !USES:
    use clm_time_manager , only : get_curr_date
    use clm_instur       , only : urban_valid
    use dshr_methods_mod , only : dshr_fldbun_getfldptr
    use dshr_strdata_mod , only : shr_strdata_advance
    use shr_infnan_mod   , only : nan => shr_infnan_nan, assignment(=)
    use UrbanParamsType  , only : urban_explicit_ac
    !
    ! !ARGUMENTS:
    class(urbantv_type)           :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    logical :: found
    integer :: l, ig, g, ip, n
    integer :: year    ! year (0, ...) for nstep+1
    integer :: mon     ! month (1, ..., 12) for nstep+1
    integer :: day     ! day of month (1, ..., 31) for nstep+1
    integer :: sec     ! seconds into current date for nstep+1
    integer :: mcdate  ! Current model date (yyyymmdd)
    integer :: lindx   ! landunit index
    integer :: gindx   ! gridcell index
    integer :: lsize
    integer :: rc
    real(r8), pointer :: dataptr1d(:)
    real(r8), pointer :: dataptr2d(:,:)
    !-----------------------------------------------------------------------

    ! Advance sdat stream
    call get_curr_date(year, mon, day, sec)
    mcdate = year*10000 + mon*100 + day
    call shr_strdata_advance(this%sdat_urbantv, ymd=mcdate, tod=sec, logunit=iulog, istr='hdmdyn', rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
       call ESMF_Finalize(endflag=ESMF_END_ABORT)
    end if

    ! Create 2d array for all stream variable data
    lsize = bounds%endg - bounds%begg + 1
    allocate(dataptr2d(lsize, stream_varname_MIN:stream_varname_MAX))
    do n = stream_varname_MIN,stream_varname_MAX
       call dshr_fldbun_getFldPtr(this%sdat_urbantv%pstrm(1)%fldbun_model, trim(stream_varnames(n)), &
            fldptr1=dataptr1d, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
       end if
       ! Note that the size of dataptr1d includes ocean points so it will be around 3x larger than lsize
       ! So an explicit loop is required here
       do g = 1,lsize
          dataptr2d(g,n) = dataptr1d(g)
       end do
    end do

    ! Determine this%tbuilding_max (and this%p_ac, if applicable) for all landunits
    do l = bounds%begl,bounds%endl
       if (lun%urbpoi(l)) then
          ! Note that since l is within [begl, endl] bounds, we can assume
          ! lun%gricell(l) is within [begg, endg]
          ig = lun%gridcell(l) - bounds%begg + 1

          do n = stream_varname_MIN,stream_varname_MAX
             if (stream_varnames((lun%itype(l)-6)) == stream_varnames(n)) then
                this%t_building_max(l) = dataptr2d(ig,n)
             end if
             if (urban_explicit_ac) then
                if (stream_varnames((lun%itype(l)-3)) == stream_varnames(n)) then
                   this%p_ac(l) = dataptr2d(ig,n)
                end if
             end if
          end do
       else
          this%t_building_max(l) = spval
          this%p_ac(l) = spval
       end if
    end do
    deallocate(dataptr2d)

    ! Error check
    found = .false.
    do l = bounds%begl,bounds%endl
       if (lun%urbpoi(l)) then
          do g = bounds%begg,bounds%endg
             if (g == lun%gridcell(l)) exit
          end do
          ! Check for valid urban data
          if ( .not. urban_valid(g) .or. (this%t_building_max(l) <= 0._r8)) then
             found = .true.
             gindx = g
             lindx = l
             exit
          else if (urban_explicit_ac .and. (this%p_ac(l) < 0._r8 .or. this%p_ac(l) > 1._r8)) then
             found = .true.
             gindx = g
             lindx = l
             exit
          end if
       end if
    end do
    if ( found ) then
       write(iulog,*)'ERROR: no valid urban data for g= ',gindx
       write(iulog,*)'landunit type:   ',lun%itype(lindx)
       write(iulog,*)'urban_valid:     ',urban_valid(gindx)
       write(iulog,*)'t_building_max:  ',this%t_building_max(lindx)
       write(iulog,*)'p_ac:            ',this%p_ac(lindx)
       call endrun(subgrid_index=lindx, subgrid_level=subgrid_level_landunit, &
            msg=errmsg(sourcefile, __LINE__))
    end if

  end subroutine urbantv_interp

end module UrbanTimeVarType
