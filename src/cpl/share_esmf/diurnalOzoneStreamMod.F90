module diurnalOzoneStreamMod

#include "shr_assert.h"

   !-----------------------------------------------------------------------
   ! !DESCRIPTION:
   ! Read in file to convert input ozone to sub-daily values from stream
   !
   ! !USES:
   use ESMF
   use dshr_strdata_mod      , only : shr_strdata_type
   use shr_kind_mod          , only : r8 => shr_kind_r8, CL => shr_kind_CL, CS => shr_kind_CS
   use shr_log_mod           , only : errMsg => shr_log_errMsg
   use spmdMod               , only : masterproc, mpicom
   use clm_varctl            , only : iulog
   use abortutils            , only : endrun
   use decompMod             , only : bounds_type
   use DiurnalOzoneType      , only : diurnal_ozone_anom_type

   !
   ! PUBLIC TYPES:
   implicit none
   private

   ! PUBLIC MEMBER FUNCTIONS:
   public :: read_O3_stream    ! read dataset for diurnal ozone anomaly

   character(len=*), parameter :: sourcefile = &
      __FILE__

   !==============================================================================
   contains

   !==============================================================================

   subroutine read_O3_stream(diurnalOzoneAnomInst, bounds)
      !
      ! Initialize data stream information for LAI.
      !
      ! USES:
      use spmdMod          , only : iam
      use clm_nlUtilsMod   , only : find_nlgroup_name
      use shr_log_mod      , only : errMsg => shr_log_errMsg
      use shr_mpi_mod      , only : shr_mpi_bcast
      use controlMod       , only : NLFilename
      use lnd_comp_shr     , only : mesh, model_clock
      use dshr_strdata_mod , only : shr_strdata_init_from_inline, shr_strdata_print
      use dshr_strdata_mod , only : shr_strdata_advance
      use dshr_methods_mod , only : dshr_fldbun_getfldptr
      use ncdio_pio        , only : file_desc_t, ncd_pio_openfile, ncd_io

      !
      ! ARGUMENTS:
      type(diurnal_ozone_anom_type), intent(inout) :: diurnalOzoneAnomInst ! instance of diurnal ozone anomaly type
      type(bounds_type),             intent(in)    :: bounds               ! bounds
      !
      ! !LOCAL VARIABLES:
      type(shr_strdata_type)      :: sdat_dO3                              ! input data stream
      type(file_desc_t)           :: ncid                                  ! netcdf file id
      real(r8), pointer           :: dataptr2d(:,:)                        ! first dimension is level, second is data on that level
      real(r8), pointer           :: secs(:)                               ! time-of-day (units = seconds) dimension on file
      integer                     :: ig, g, j                              ! indices
      integer                     :: nu_nml                                ! unit for namelist file
      integer                     :: nml_error                             ! namelist i/o error flag
      integer                     :: year                                  ! year (0, ...) for nstep+1
      integer                     :: mon                                   ! month (1, ..., 12) for nstep+1
      integer                     :: day                                   ! day of month (1, ..., 31) for nstep+1
      integer                     :: sec                                   ! seconds into current date for nstep+1
      integer                     :: mcdate                                ! current model date (yyyymmdd)
      integer                     :: rc                                    ! error code
      integer                     :: nlevsec                               ! dimension of 'secs' variable
      integer, allocatable        :: g_to_ig(:)                            ! array matching gridcell index to data index
      logical                     :: readvar 
      character(len=CL)           :: stream_fldFileName_dO3 = ' '          ! diurnal ozone stream filename to read
      character(len=CL)           :: stream_meshfile_dO3 = ' '             ! diurnal ozone stream meshfile
      character(len=CL)           :: dO3_mapalgo = ' '                     ! mapping alogrithm
      character(len=CL)           :: stream_lev_dimname = 'sec'            ! name of vertical layer dimension
      character(*), parameter     :: stream_var_name = "O3_diurnal_factor" ! base string for field string
      character(len=*), parameter :: subName = "('dO3_init')"

      !-----------------------------------------------------------------------
      !
      ! namelist variables
      !
      namelist /dO3_streams/     &
         dO3_mapalgo,            &
         stream_fldFileName_dO3, &
         stream_meshfile_dO3

      ! Read dO3_streams namelist
      if (masterproc) then
         open( newunit=nu_nml, file=trim(NLFilename), status='old', iostat=nml_error )
         call find_nlgroup_name(nu_nml, 'do3_streams', status=nml_error)
         if (nml_error == 0) then
            read(nu_nml, nml=dO3_streams,iostat=nml_error)
            if (nml_error /= 0) then
               call endrun(subname // ':: ERROR reading do3_streams namelist')
            end if
         else
            call endrun(subname // ':: ERROR finding do3_streams namelist')
         end if
         close(nu_nml)
      endif

      call shr_mpi_bcast(stream_fldFileName_dO3 , mpicom)
      call shr_mpi_bcast(stream_meshfile_dO3    , mpicom)
      call shr_mpi_bcast(dO3_mapalgo            , mpicom)

      if (masterproc) then
         write(iulog,*)
         write(iulog,'(a)') 'do3_stream settings:'
         write(iulog,'(a,a)' ) '  stream_fldFileName_do3 = ',trim(stream_fldFileName_dO3)
         write(iulog,'(a,a)' ) '  stream_meshfile_do3    = ',trim(stream_meshfile_dO3)
         write(iulog,'(a,a)' ) '  stream varname         = ',trim(stream_var_name)
         write(iulog,*)
      endif

      ! Initialize the cdeps data type sdat_dO3
      call shr_strdata_init_from_inline(sdat_dO3,                &
         my_task             = iam,                              &
         logunit             = iulog,                            &
         compname            = 'LND',                            &
         model_clock         = model_clock,                      &
         model_mesh          = mesh,                             &
         stream_meshfile     = trim(stream_meshfile_dO3),        &
         stream_lev_dimname  = trim(stream_lev_dimname),         & 
         stream_mapalgo      = trim(dO3_mapalgo),                &
         stream_filenames    = (/trim(stream_fldfilename_dO3)/), &
         stream_fldlistFile  = (/trim(stream_var_name)/),        &
         stream_fldListModel = (/trim(stream_var_name)/),        &
         stream_yearFirst    = 2000,                             &
         stream_yearLast     = 2000,                             &
         stream_yearAlign    = 1,                                &
         stream_offset       = 0,                                &
         stream_taxmode      = 'extend',                         &
         stream_dtlimit      = 1.0e30_r8,                        &
         stream_tintalgo     = 'linear',                         &
         stream_name         = 'do3 data',                       &
         rc                  = rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
         call ESMF_Finalize(endflag=ESMF_END_ABORT)
      end if

      ! Explicitly set current date to a hardcoded constant value
      ! as in ch4FInundatedStreamType
      year = 2000
      mon  = 12
      day  = 31
      sec  = 0
      mcdate = year*10000 + mon*100 + day

      call shr_strdata_advance(sdat_dO3, ymd=mcdate, tod=sec, logunit=iulog, istr='do3', rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
         call ESMF_Finalize(endflag=ESMF_END_ABORT)
      end if

      ! Map gridcell to 1->local_size 
      if ( .not. allocated(g_to_ig) )then
         allocate (g_to_ig(bounds%begg:bounds%endg) )
         ig = 0
         do g = bounds%begg,bounds%endg
            ig = ig+1
            g_to_ig(g) = ig
         end do
      end if

      ! Get pointer for stream data that is time and spatially interpolated to model time and grid
      call dshr_fldbun_getFldPtr(sdat_dO3%pstrm(1)%fldbun_model, trim(stream_var_name), fldptr2=dataptr2d, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
         call ESMF_Finalize(endflag=ESMF_END_ABORT)
      end if

      ! set the nlevsec size
      nlevsec = size(dataptr2d, dim=1)

      ! read in the seconds array as well
      allocate(secs(nlevsec))
      call ncd_pio_openfile(ncid, trim(stream_fldFileName_dO3), 0)
      call ncd_io(ncid=ncid, varname=trim(stream_lev_dimname), flag='read', data=secs, readvar=readvar)
      if (.not. readvar) then
         call endrun(msg=' ERROR: secs NOT on diurnal ozone file'//errMsg(sourcefile, __LINE__))
      end if

      ! initialize arrays
      call diurnalOzoneAnomInst%Init(bounds, nlevsec)
      
      ! set the diurnal ozone anomaly
      do g = bounds%begg, bounds%endg
         ig = g_to_ig(g)
         do j = 1, nlevsec
            diurnalOzoneAnomInst%o3_anomaly_grc(g,j) = dataptr2d(j,ig)
         end do
      end do

      ! set the seconds array
      do j = 1, nlevsec
         diurnalOzoneAnomInst%time_arr(j) = secs(j)
      end do

   end subroutine read_O3_stream

!==============================================================================

end module diurnalOzoneStreamMod
