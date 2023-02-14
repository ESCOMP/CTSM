module UrbanTimeVarType

  !------------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Urban Time Varying Data
  !
  ! !USES:
  use shr_kind_mod    , only : r8 => shr_kind_r8, CL => shr_kind_CL
  use shr_log_mod     , only : errMsg => shr_log_errMsg
  use abortutils      , only : endrun
  use decompMod       , only : bounds_type
  use clm_varctl      , only : iulog, inst_name
  use landunit_varcon , only : isturb_MIN, isturb_MAX
  use clm_varcon      , only : spval
  use LandunitType    , only : lun
  use GridcellType    , only : grc
  use mct_mod
  use shr_strdata_mod , only : shr_strdata_type
  !
  implicit none
  save
  private
  !
  !

  ! !PUBLIC TYPE
  type, public :: urbantv_type

     real(r8), public, pointer     :: t_building_max(:)    ! lun maximum internal building air temperature (K)
     type(shr_strdata_type)        :: sdat_urbantv         ! urban time varying input data stream
   contains

     ! !PUBLIC MEMBER FUNCTIONS:
     procedure, public :: Init              ! Allocate and initialize urbantv
     procedure, public :: urbantv_init      ! Initialize urban time varying stream
     procedure, public :: urbantv_interp    ! Interpolate urban time varying stream

  end type urbantv_type

  !-----------------------------------------------------------------------
  character(15), private :: stream_var_name(isturb_MIN:isturb_MAX)

  character(len=*), parameter, private :: sourcefile = &
       __FILE__

contains

  !-----------------------------------------------------------------------
  subroutine Init(this, bounds, NLFilename)
    !
    ! Allocate module variables and data structures
    !
    ! !USES:
    use shr_infnan_mod  , only : nan => shr_infnan_nan, assignment(=)
    use histFileMod     , only : hist_addfld1d
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

    ! Allocate urbantv data structure

    allocate(this%t_building_max      (begl:endl))          ; this%t_building_max      (:)   = nan

    call this%urbantv_init(bounds, NLFilename)
    call this%urbantv_interp(bounds)

    ! Add history fields
    call hist_addfld1d (fname='TBUILD_MAX', units='K',      &
          avgflag='A', long_name='prescribed maximum interior building temperature',   &
          ptr_lunit=this%t_building_max, default='inactive', set_nourb=spval, &
          l2g_scale_type='unity')


  end subroutine Init

  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  subroutine urbantv_init(this, bounds, NLFilename)
   !
   ! !DESCRIPTION:
   ! Initialize data stream information for urban time varying data
   !
   ! !USES:
   use clm_time_manager , only : get_calendar
   use ncdio_pio        , only : pio_subsystem
   use shr_pio_mod      , only : shr_pio_getiotype
   use clm_nlUtilsMod   , only : find_nlgroup_name
   use ndepStreamMod    , only : clm_domain_mct
   use spmdMod          , only : masterproc, mpicom, comp_id
   use fileutils        , only : getavu, relavu
   use shr_mpi_mod      , only : shr_mpi_bcast
   use shr_string_mod   , only : shr_string_listAppend
   use shr_strdata_mod  , only : shr_strdata_create, shr_strdata_print
   use decompMod        , only : gsmap_lnd_gdc2glo
   use domainMod        , only : ldomain
   use shr_infnan_mod   , only : nan => shr_infnan_nan, assignment(=)
   use landunit_varcon  , only : isturb_TBD, isturb_HD, isturb_MD
   !
   ! !ARGUMENTS:
   implicit none
   class(urbantv_type)           :: this
   type(bounds_type), intent(in) :: bounds
   character(len=*),  intent(in) :: NLFilename   ! Namelist filename
   !
   ! !LOCAL VARIABLES:
   integer            :: begl, endl                           ! landunits
   integer            :: ifield                               ! field index
   integer            :: stream_year_first_urbantv            ! first year in urban tv stream to use
   integer            :: stream_year_last_urbantv             ! last year in urban tv stream to use
   integer            :: model_year_align_urbantv             ! align stream_year_first_urbantv
                                                              !  with this model year
   integer            :: nu_nml                               ! unit for namelist file
   integer            :: nml_error                            ! namelist i/o error flag
   type(mct_ggrid)    :: dom_clm                              ! domain information
   character(len=CL)  :: stream_fldFileName_urbantv           ! urban tv streams filename
   character(len=CL)  :: urbantvmapalgo = 'nn'                ! mapping alogrithm for urban ac
   character(len=CL)  :: urbantv_tintalgo = 'linear'          ! time interpolation alogrithm
   character(len=CL)  :: fldList                              ! field string
   character(*), parameter :: urbantvString = "tbuildmax_"    ! base string for field string
   character(*), parameter :: subName = "('urbantv_init')"
   character(*), parameter :: F00 = "('(urbantv_init) ',4a)"
   !-----------------------------------------------------------------------
   namelist /urbantv_streams/       &
        stream_year_first_urbantv,  &
        stream_year_last_urbantv,   &
        model_year_align_urbantv,   &
        urbantvmapalgo,             &
        stream_fldFileName_urbantv, &
        urbantv_tintalgo
   !-----------------------------------------------------------------------

   begl = bounds%begl; endl = bounds%endl

   ! Default values for namelist
   stream_year_first_urbantv  = 1       ! first year in stream to use
   stream_year_last_urbantv   = 1       ! last  year in stream to use
   model_year_align_urbantv   = 1       ! align stream_year_first_urbantv with this model year
   stream_fldFileName_urbantv = ' '

   ! Read urbantv_streams namelist
   if (masterproc) then
      nu_nml = getavu()
      open( nu_nml, file=trim(NLFilename), status='old', iostat=nml_error )
      call find_nlgroup_name(nu_nml, 'urbantv_streams', status=nml_error)
      if (nml_error == 0) then
         read(nu_nml, nml=urbantv_streams,iostat=nml_error)
         if (nml_error /= 0) then
            call endrun(msg='ERROR reading urbantv_streams namelist'//errMsg(sourcefile, __LINE__))
         end if
      end if
      close(nu_nml)
      call relavu( nu_nml )
   endif

   call shr_mpi_bcast(stream_year_first_urbantv, mpicom)
   call shr_mpi_bcast(stream_year_last_urbantv, mpicom)
   call shr_mpi_bcast(model_year_align_urbantv, mpicom)
   call shr_mpi_bcast(stream_fldFileName_urbantv, mpicom)
   call shr_mpi_bcast(urbantv_tintalgo, mpicom)

   if (masterproc) then
      write(iulog,*) ' '
      write(iulog,*) 'urbantv_streams settings:'
      write(iulog,*) '  stream_year_first_urbantv  = ',stream_year_first_urbantv
      write(iulog,*) '  stream_year_last_urbantv   = ',stream_year_last_urbantv
      write(iulog,*) '  model_year_align_urbantv   = ',model_year_align_urbantv
      write(iulog,*) '  stream_fldFileName_urbantv = ',stream_fldFileName_urbantv
      write(iulog,*) '  urbantv_tintalgo           = ',urbantv_tintalgo
      write(iulog,*) ' '
   endif

   call clm_domain_mct (bounds, dom_clm)

   ! create the field list for these urbantv fields...use in shr_strdata_create
   stream_var_name(:)          = "NOT_SET"
   stream_var_name(isturb_TBD) = urbantvString//"TBD"
   stream_var_name(isturb_HD)  = urbantvString//"HD"
   stream_var_name(isturb_MD)  = urbantvString//"MD"
   fldList = ""
   do ifield = isturb_MIN, isturb_MAX
      call shr_string_listAppend( fldList, stream_var_name(ifield) )
   end do

   call shr_strdata_create(this%sdat_urbantv,name="clmurbantv",     &
        pio_subsystem=pio_subsystem,                   &
        pio_iotype=shr_pio_getiotype(inst_name),           &
        mpicom=mpicom, compid=comp_id,                 &
        gsmap=gsmap_lnd_gdc2glo, ggrid=dom_clm,        &
        nxg=ldomain%ni, nyg=ldomain%nj,                &
        yearFirst=stream_year_first_urbantv,           &
        yearLast=stream_year_last_urbantv,             &
        yearAlign=model_year_align_urbantv,            &
        offset=0,                                      &
        domFilePath='',                                &
        domFileName=trim(stream_fldFileName_urbantv),  &
        domTvarName='time',                            &
        domXvarName='lon' ,                            &
        domYvarName='lat' ,                            &
        domAreaName='area',                            &
        domMaskName='LANDMASK',                        &
        filePath='',                                   &
        filename=(/trim(stream_fldFileName_urbantv)/) , &
        fldListFile=fldList,                           &
        fldListModel=fldList,                          &
        fillalgo='none',                               &
        mapalgo=urbantvmapalgo,                        &
        calendar=get_calendar(),                       &
        tintalgo=urbantv_tintalgo,                     &
        taxmode='extend'                                 )

   if (masterproc) then
      call shr_strdata_print(this%sdat_urbantv,'urban time varying data')
   endif


  end subroutine urbantv_init

  !-----------------------------------------------------------------------
  subroutine urbantv_interp(this, bounds)
  !
  ! !DESCRIPTION:
  ! Interpolate data stream information for urban time varying data.
  !
  ! !USES:
  use clm_time_manager, only : get_curr_date
  use spmdMod         , only : mpicom
  use shr_strdata_mod , only : shr_strdata_advance
  use clm_instur      , only : urban_valid
  !
  ! !ARGUMENTS:
  class(urbantv_type)           :: this
  type(bounds_type), intent(in) :: bounds
  !
  ! !LOCAL VARIABLES:
  logical :: found
  integer :: l, glun, ig, g, ip
  integer :: year    ! year (0, ...) for nstep+1
  integer :: mon     ! month (1, ..., 12) for nstep+1
  integer :: day     ! day of month (1, ..., 31) for nstep+1
  integer :: sec     ! seconds into current date for nstep+1
  integer :: mcdate  ! Current model date (yyyymmdd)
  integer :: lindx   ! landunit index
  integer :: gindx   ! gridcell index
  !-----------------------------------------------------------------------

   call get_curr_date(year, mon, day, sec)
   mcdate = year*10000 + mon*100 + day

   call shr_strdata_advance(this%sdat_urbantv, mcdate, sec, mpicom, 'urbantvdyn')

   do l = bounds%begl,bounds%endl
      if (lun%urbpoi(l)) then
         glun  = lun%gridcell(l)
         ip = mct_aVect_indexRA(this%sdat_urbantv%avs(1),trim(stream_var_name(lun%itype(l))))
         !
         ! Determine vector index corresponding to glun
         !
         ig = 0
         do g = bounds%begg,bounds%endg
            ig = ig+1
            if (g == glun) exit
         end do

         this%t_building_max(l) = this%sdat_urbantv%avs(1)%rAttr(ip,ig)
      else
         this%t_building_max(l) = spval
      end if
   end do

   found = .false.
   do l = bounds%begl,bounds%endl
      if (lun%urbpoi(l)) then
         glun  = lun%gridcell(l)
         !
         ! Determine vector index corresponding to glun
         !
         ig = 0
         do g = bounds%begg,bounds%endg
            ig = ig+1
            if (g == glun) exit
         end do

         if ( .not. urban_valid(g) .or. (this%t_building_max(l) <= 0._r8)) then
            found = .true.
            gindx = g
            lindx = l
            exit
         end if
      end if
   end do
   if ( found ) then
      write(iulog,*)'ERROR: no valid urban data for g= ',gindx
      write(iulog,*)'landunit type:   ',lun%itype(l)
      write(iulog,*)'urban_valid:     ',urban_valid(gindx)
      write(iulog,*)'t_building_max:  ',this%t_building_max(lindx)
      call endrun(msg=errmsg(sourcefile, __LINE__))
   end if


  end subroutine urbantv_interp

  !-----------------------------------------------------------------------

end module UrbanTimeVarType
