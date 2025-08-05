module DistParamsStreamMod


  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Contains methods for reading from a spatially distributed parameter file

  ! !USES
  use ESMF
  use dshr_strdata_mod , only : shr_strdata_type
  use shr_kind_mod     , only : r8 => shr_kind_r8, CL => shr_kind_cl
  use shr_infnan_mod   , only : nan => shr_infnan_nan, assignment(=)
  use shr_log_mod      , only : errMsg => shr_log_errMsg
  use spmdMod          , only : mpicom, masterproc
  use clm_varctl       , only : iulog, FL => fname_len
  use clm_varctl       , only : use_distributed_parameters
  use abortutils       , only : endrun
  use decompMod        , only : bounds_type
  use DistParamType    , only : distributed_parameter_type, distparam_class

  ! !PUBLIC TYPES:
  implicit none
  private

  !
  type, public :: distributed_parameter_stream_type
     !
   contains

      ! !PUBLIC MEMBER FUNCTIONS:
      procedure, public :: Init            ! Initialize and read data in
      procedure, public :: Clean           ! Clean and deallocate the object

  end type distributed_parameter_stream_type

  type(distributed_parameter_stream_type), public, target :: distributed_parameter_stream 

  ! ! PRIVATE DATA:
  type, private :: streamcontrol_type
     character(len=FL)  :: stream_fldFileName_distparams   ! data Filename
     character(len=FL)  :: stream_meshfile_distparams      ! mesh Filename
     character(len=CL)  :: mapalgo_distparams              ! map algo
  contains
     procedure, private :: ReadNML      ! Read in control namelist
  end type streamcontrol_type

  type(streamcontrol_type), private :: control             ! Stream control data
  logical                 , private :: NMLRead = .false.   ! If namelist has been read
  logical                 , private :: InitDone = .false.  ! If initialization of streams has been done

  character(len=*), parameter, private :: sourcefile = &
       __FILE__

!==============================================================================
contains
!==============================================================================

  !------------------------------------------------------------------------
  subroutine AssignDistributedParameter(this,bounds,dataptr1d,stream_varname)
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    implicit none
    class(distparam_class), intent(inout)       :: this           ! distributed parameter
    type(bounds_type)     , intent(in)          :: bounds
    real(r8)              , pointer, intent(in) :: dataptr1d(:)   ! stream data
    character(len=*)      , intent(in)          :: stream_varname ! name of parameter on stream file
   !
    ! local variables
    integer                     :: g, ig
    character(len=*), parameter :: subname = 'AssignDistributedParameter'
    !--------------------------------------------------------------------

    if (trim(stream_varname) == trim(this%name)) then
       allocate(this%val(bounds%begg:bounds%endg))
       this%val(:) = nan
       ig = 0
       do g = bounds%begg,bounds%endg
          ig = ig+1
          this%val(g) = dataptr1d(ig)
       end do
       this%is_distributed = .true.
    end if
    
  end subroutine AssignDistributedParameter
  
  !------------------------------------------------------------------------
  subroutine Init(this, bounds, NLFilename, distparams)
   !
   ! Initialize the distributed parameters stream object
   !
   ! Uses:
   use spmdMod          , only : iam
   use lnd_comp_shr     , only : mesh, model_clock
   use dshr_strdata_mod , only : shr_strdata_init_from_inline, shr_strdata_print
   use dshr_strdata_mod , only : shr_strdata_advance
   use dshr_methods_mod , only : dshr_fldbun_getfldptr
   use fileutils        , only : getfil
   use ncdio_pio        , only : ncd_pio_closefile, ncd_pio_openfile
   use ncdio_pio        , only : file_desc_t, var_desc_t
   use ncdio_pio        , only : ncd_inqvname, ncd_inqnvars
   !
   ! arguments
   implicit none
   class(distributed_parameter_stream_type) :: this
   type(bounds_type), intent(in)   :: bounds
   character(len=*),  intent(in)   :: NLFilename   ! Namelist filename
   type(distributed_parameter_type), intent(inout) :: distparams
   !
   ! local variables
   integer                        :: ig, g, n           ! Indices
   integer                        :: begg, endg
   integer                        :: year               ! year (0, ...) for nstep+1
   integer                        :: mon                ! month (1, ..., 12) for nstep+1
   integer                        :: day                ! day of month (1, ..., 31) for nstep+1
   integer                        :: sec                ! seconds into current date for nstep+1
   integer                        :: mcdate             ! Current model date (yyyymmdd)
   integer                        :: nparams            ! number of parameters on file
   integer                        :: nvariables         ! number of variables on file
   integer                        :: dimid              ! netCDF dimension id
   type(file_desc_t)              :: ncid               ! pio netCDF file id
   character(len=256)             :: locfn              ! local filename
   type(shr_strdata_type)         :: sdat_distparams    ! input data stream
   character(len=256)             :: varname            ! variable name
   type(var_desc_t)               :: var_desc           ! variable descriptor
   character(len=CL), allocatable :: stream_varnames_in(:) ! array of stream field names
   character(len=CL), allocatable :: stream_varnames(:) ! array of stream field names
   integer                        :: rc                 ! error code
   real(r8), pointer              :: dataptr1d(:)       ! temporary pointer
   character(len=*), parameter    :: stream_name = 'Distributed parameters'
   character(len=*), parameter    :: subname     = 'DistParamsStream::Init'
   !-----------------------------------------------------------------------

      call control%ReadNML( bounds, NLFileName )

      if ( use_distributed_parameters )then

         ! query parameter file
         call getfil (control%stream_fldFileName_distparams, locfn, 0)
         call ncd_pio_openfile (ncid, trim(locfn), 0)

         call ncd_inqnvars(ncid,nvariables=nvariables)

         allocate(stream_varnames_in(nvariables))
         nparams = 0
         do n=1,nvariables
            call ncd_inqvname(ncid,n,varname,var_desc)
            if (.not. any(varname == (/'time','lat','lon'/))) then 
               nparams = nparams + 1
               write(stream_varnames_in(nparams),'(a)') trim(varname)
            endif
         enddo
         
         ! close parameter file
         call ncd_pio_closefile(ncid)

         ! keep valid variable names (not entirely satisfying, but avoids having to explicitly index array later
         allocate(stream_varnames(nparams))
         stream_varnames = stream_varnames_in(1:nparams)
         deallocate(stream_varnames_in)

         if (masterproc) then
            write(iulog,*) '  stream_varnames                  = ',stream_varnames
         end if

         ! Initialize the cdeps data type sdat_distparams
         call shr_strdata_init_from_inline(sdat_distparams,                                   &
              my_task             = iam,                                                &
              logunit             = iulog,                                              &
              compname            = 'LND',                                              &
              model_clock         = model_clock,                                        &
              model_mesh          = mesh,                                               &
              stream_meshfile     = control%stream_meshfile_distparams,                 &
              stream_lev_dimname  = 'null',                                             &
              stream_mapalgo      = control%mapalgo_distparams,                         &
              stream_filenames    = (/trim(control%stream_fldFileName_distparams)/),    &
!              stream_fldlistFile  = stream_varnames(:nparams),                                    &
!              stream_fldListModel = stream_varnames(:nparams),                                    &
              stream_fldlistFile  = stream_varnames,                                    &
              stream_fldListModel = stream_varnames,                                    &
              stream_yearFirst    = 2000,                                               &
              stream_yearLast     = 2000,                                               &
              stream_yearAlign    = 1,                                                  &
              stream_offset       = 0,                                                  &
              stream_taxmode      = 'extend',                                           &
              stream_dtlimit      = 1.0e30_r8,                                          &
              stream_tintalgo     = 'linear',                                           &
              stream_name         = 'Distributed parameters',                           &
              rc                  = rc)
         if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
            call ESMF_Finalize(endflag=ESMF_END_ABORT)
         end if

         ! Explicitly set current date to a hardcoded constant value. Otherwise
         ! using the real date can cause roundoff differences that are
         ! detrected as issues with exact restart.  EBK M05/20/2017
         ! call get_curr_date(year, mon, day, sec)
         year = 2000
         mon  = 12
         day  = 31
         sec  = 0
         mcdate = year*10000 + mon*100 + day

         call shr_strdata_advance(sdat_distparams, ymd=mcdate, tod=sec, logunit=iulog, istr='distparams', rc=rc)
         if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
            call ESMF_Finalize(endflag=ESMF_END_ABORT)
         end if

         ! Get pointer for stream data that is time and spatially interpolate to model time and grid
         begg = bounds%begg; endg = bounds%endg

         do n = 1,size(stream_varnames)
            call dshr_fldbun_getFldPtr(sdat_distparams%pstrm(1)%fldbun_model, stream_varnames(n), fldptr1=dataptr1d, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
               call ESMF_Finalize(endflag=ESMF_END_ABORT)
            end if

            ! Assign values to distributed parameter arrays if stream data exist
            call AssignDistributedParameter(distparams%aq_sp_yield_min,bounds,dataptr1d,stream_varnames(n))
            call AssignDistributedParameter(distparams%n_baseflow,bounds,dataptr1d,stream_varnames(n))
            call AssignDistributedParameter(distparams%perched_baseflow_scalar,bounds,dataptr1d,stream_varnames(n))
            call AssignDistributedParameter(distparams%e_ice,bounds,dataptr1d,stream_varnames(n))
            call AssignDistributedParameter(distparams%baseflow_scalar,bounds,dataptr1d,stream_varnames(n))
            call AssignDistributedParameter(distparams%fff,bounds,dataptr1d,stream_varnames(n))
            call AssignDistributedParameter(distparams%slopebeta,bounds,dataptr1d,stream_varnames(n))
            call AssignDistributedParameter(distparams%slopemax,bounds,dataptr1d,stream_varnames(n))
            call AssignDistributedParameter(distparams%n_melt_coef,bounds,dataptr1d,stream_varnames(n))
            call AssignDistributedParameter(distparams%accum_factor,bounds,dataptr1d,stream_varnames(n))
            call AssignDistributedParameter(distparams%upplim_destruct_metamorph,bounds,dataptr1d,stream_varnames(n))
            call AssignDistributedParameter(distparams%pc,bounds,dataptr1d,stream_varnames(n))
            call AssignDistributedParameter(distparams%mu,bounds,dataptr1d,stream_varnames(n))
            call AssignDistributedParameter(distparams%bsw_sf,bounds,dataptr1d,stream_varnames(n))
            call AssignDistributedParameter(distparams%hksat_sf,bounds,dataptr1d,stream_varnames(n))
            call AssignDistributedParameter(distparams%sucsat_sf,bounds,dataptr1d,stream_varnames(n))
            call AssignDistributedParameter(distparams%watsat_sf,bounds,dataptr1d,stream_varnames(n))
            call AssignDistributedParameter(distparams%d_max,bounds,dataptr1d,stream_varnames(n))
            call AssignDistributedParameter(distparams%frac_sat_soil_dsl_init,bounds,dataptr1d,stream_varnames(n))
            call AssignDistributedParameter(distparams%zlnd,bounds,dataptr1d,stream_varnames(n))
            call AssignDistributedParameter(distparams%snw_rds_min,bounds,dataptr1d,stream_varnames(n))
            call AssignDistributedParameter(distparams%precip_repartition_nonglc_all_rain_t,bounds,dataptr1d,stream_varnames(n))
            call AssignDistributedParameter(distparams%precip_repartition_nonglc_all_snow_t,bounds,dataptr1d,stream_varnames(n))

         end do
         deallocate(stream_varnames)
         InitDone = .true.
      end if

  end subroutine Init

  !==============================================================================
  subroutine Clean(this)
   !
   ! Deallocate and clean the object
   !
   ! Uses:
   !
   ! arguments
   implicit none
   class(distributed_parameter_stream_type) :: this
   !
   ! local variables
   !-----------------------------------------------------------------------
   !deallocate(this%distparams)
   !this%distparams => NULL()
   InitDone = .false.

  end subroutine Clean

  !==============================================================================
  subroutine ReadNML(this, bounds, NLFilename)
   !
   ! Read the namelist data stream information.
   !
   ! Uses:
   use shr_nl_mod       , only : shr_nl_find_group_name
   use shr_log_mod      , only : errMsg => shr_log_errMsg
   use shr_mpi_mod      , only : shr_mpi_bcast
   !
   ! arguments
   implicit none
   class(streamcontrol_type)     :: this
   type(bounds_type), intent(in) :: bounds
   character(len=*),  intent(in) :: NLFilename   ! Namelist filename
   !
   ! local variables
   integer            :: nu_nml    ! unit for namelist file
   integer            :: nml_error ! namelist i/o error flag
   character(len=FL)  :: stream_fldFileName_distparams = ' '
   character(len=FL)  :: stream_meshfile_distparams = ' '
   character(len=CL)  :: mapalgo_distparams = 'bilinear'
   character(len=*), parameter :: namelist_name = 'distparams_streams'    ! MUST agree with group name in namelist definition to read.
   character(len=*), parameter :: subName = "('distparams_streams::ReadNML')"
   !-----------------------------------------------------------------------

   namelist /distparams_streams/ & ! MUST agree with namelist_name above
        use_distributed_parameters, stream_fldFileName_distparams, &
        stream_meshfile_distparams

   ! Default values for namelist

   ! Read distparams_streams namelist
   if (masterproc) then
      open( newunit=nu_nml, file=trim(NLFilename), status='old', iostat=nml_error )
      call shr_nl_find_group_name(nu_nml, namelist_name, status=nml_error)
      if (nml_error == 0) then
         read(nu_nml, nml=distparams_streams,iostat=nml_error)   ! MUST agree with namelist_name above
         if (nml_error /= 0) then
            call endrun(msg=' ERROR reading '//namelist_name//' namelist'//errMsg(sourcefile, __LINE__))
         end if
      else
         call endrun(msg=' ERROR finding '//namelist_name//' namelist'//errMsg(sourcefile, __LINE__))
      end if
      close(nu_nml)
   endif

   call shr_mpi_bcast(use_distributed_parameters    , mpicom)
   call shr_mpi_bcast(mapalgo_distparams            , mpicom)
   call shr_mpi_bcast(stream_fldFileName_distparams , mpicom)
   call shr_mpi_bcast(stream_meshfile_distparams    , mpicom)

   ! Error checking
   if ( .not. use_distributed_parameters )then
      if ( len_trim(stream_fldFileName_distparams) /= 0 )then
         call endrun(msg=' ERROR stream_fldFileName_distparams is set, but use_distributed_parameters is FALSE' &
                     //errMsg(sourcefile, __LINE__))
      end if
      if ( len_trim(stream_meshfile_distparams) /= 0 )then
         call endrun(msg=' ERROR stream_meshfile_distparams is set, but use_distributed_parameters is FALSE' &
                     //errMsg(sourcefile, __LINE__))
      end if
   else
      if ( len_trim(stream_fldFileName_distparams) == 0 )then
         call endrun(msg=' ERROR stream_fldFileName_distparams is NOT set, but use_distributed_parameters is TRUE' &
                     //errMsg(sourcefile, __LINE__))
      end if
      if ( len_trim(stream_meshfile_distparams) == 0 )then
         call endrun(msg=' ERROR stream_meshfile_distparams is NOT set, but use_distributed_parameters is TRUE' &
                    //errMsg(sourcefile, __LINE__))
      end if
   end if

   if (masterproc) then
      write(iulog,*) ' '
      write(iulog,*) namelist_name, ' stream settings:'
      write(iulog,*) '  use_distributed_parameters               = ',use_distributed_parameters
      if ( use_distributed_parameters )then
         write(iulog,*) '  stream_fldFileName_distparams = ',trim(stream_fldFileName_distparams)
         write(iulog,*) '  stream_meshfile_distparams    = ',trim(stream_meshfile_distparams)
         write(iulog,*) '  mapalgo_distparams             = ',trim(mapalgo_distparams)
      end if
   endif
   this%stream_fldFileName_distparams = stream_fldFileName_distparams
   this%stream_meshfile_distparams    = stream_meshfile_distparams
   this%mapalgo_distparams             = mapalgo_distparams

   ! Mark namelist read as having been done
   NMLRead = .true.

 end subroutine ReadNML

end module DistParamsStreamMod
