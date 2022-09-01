module diurnalOzoneStreamMod

#include "shr_assert.h"

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Read in file to convert input ozone to sub-daily values from stream
  !
  ! !USES:
  use ESMF
  use shr_kind_mod     , only : r8 => shr_kind_r8, CL => shr_kind_CL, CS => shr_kind_CS
  use dshr_strdata_mod , only : shr_strdata_type
  use decompMod        , only : bounds_type
  use abortutils       , only : endrun
  use clm_varctl       , only : iulog
  use perf_mod         , only : t_startf, t_stopf
  use spmdMod          , only : masterproc, mpicom, iam
  !
  ! !PUBLIC TYPES:
  implicit none
  private

  ! !PUBLIC MEMBER FUNCTIONS:
  public :: dO3_init    ! position dataset for diurnal ozone anomaly

  ! !PRIVATE MEMBER DATA:
  integer, allocatable        :: g_to_ig(:)         ! Array matching gridcell index to data index
  type(shr_strdata_type)      :: sdat_dO3           ! diurnal ozone anomaly input data stream
  character(*), parameter     :: dO3String = "dO3_" ! base string for field string
  integer     , parameter     :: numdO3Fields = 16  ! number of fields to build field string
  character(len=CS)           :: stream_varnames(numdO3Fields)

  character(len=*), parameter :: sourcefile = &
       __FILE__

!==============================================================================
contains
!==============================================================================

  subroutine dO3_init(bounds)
    !
    ! Initialize data stream information for LAI.
    !
    ! !USES:
    use shr_mpi_mod      , only : shr_mpi_bcast
    use clm_nlUtilsMod   , only : find_nlgroup_name
    use lnd_comp_shr     , only : mesh, model_clock
    use dshr_strdata_mod , only : shr_strdata_init_from_inline
    use controlMod       , only : NLFilename
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds          ! bounds
    !
    ! !LOCAL VARIABLES:
    integer                 :: i,n                        ! index
    integer                 :: nu_nml                     ! unit for namelist file
    integer                 :: nml_error                  ! namelist i/o error flag
    character(len=CL)       :: stream_fldFileName_dO3     ! diurnal ozone stream filename to read
    character(len=CL)       :: stream_meshfile_dO3        ! diurnal ozone stream meshfile
    character(len=CL)       :: dO3_mapalgo = 'bilinear'   ! Mapping alogrithm
    integer                 :: rc
    character(*), parameter :: subName = "('dO3dyn_init')"
    !-----------------------------------------------------------------------
    !
    ! deal with namelist variables here in init
    !
    namelist /dO3_streams/         &
         dO3_mapalgo,              &
         stream_fldFileName_dO3,   &
         stream_meshfile_dO3

    ! Default values for namelist
    stream_fldFileName_dO3 = ''
    stream_meshfile_dO3    = ''
    do n = 1,numdO3Fields
       write(stream_varnames(n),'(a,i0)') dO3String,n
    end do

    ! Read lai_streams namelist
    if (masterproc) then
       open( newunit=nu_nml, file=trim(NLFilename), status='old', iostat=nml_error )
       call find_nlgroup_name(nu_nml, 'dO3_streams', status=nml_error)
       if (nml_error == 0) then
          read(nu_nml, nml=dO3_streams,iostat=nml_error)
          if (nml_error /= 0) then
             call endrun(subname // ':: ERROR reading dO3_streams namelist')
          end if
       else
          call endrun(subname // ':: ERROR finding dO3_streams namelist')
       end if
       close(nu_nml)
    endif
    call shr_mpi_bcast(stream_fldFileName_dO3 , mpicom)
    call shr_mpi_bcast(stream_meshfile_dO3    , mpicom)
    call shr_mpi_bcast(dO3_tintalgo           , mpicom)

    if (masterproc) then
       write(iulog,*)
       write(iulog,'(a)') 'dO3_stream settings:'
       write(iulog,'(a,a)' ) '  stream_fldFileName_dO3 = ',trim(stream_fldFileName_dO3)
       write(iulog,'(a,a)' ) '  stream_meshfile_dO3    = ',trim(stream_meshfile_dO3)
       do n = 1,numdO3Fields
          write(iulog,'(a,a)' ) '  stream_varname         = ',trim(stream_varnames(n))
       end do
       write(iulog,*)
    endif

    ! Initialize the cdeps data type sdat_lai
    call shr_strdata_init_from_inline(sdat_dO3,                  &
         my_task             = iam,                              &
         logunit             = iulog,                            &
         compname            = 'LND',                            &
         model_clock         = model_clock,                      &
         model_mesh          = mesh,                             &
         stream_meshfile     = trim(stream_meshfile_dO3),        &
         stream_lev_dimname  = 'null',                           &
         stream_mapalgo      = trim(dO3_mapalgo),                &
         stream_filenames    = (/trim(stream_fldfilename_dO3)/), &
         stream_fldlistFile  = stream_varnames,                  &
         stream_fldListModel = stream_varnames,                  &
         stream_yearFirst    = 'null',                           &
         stream_yearLast     = 'null',                           &
         stream_yearAlign    = 'null',                           &
         stream_offset       = 'null',                           &
         stream_taxmode      = 'null',                           &
         stream_dtlimit      = 'null',                           &
         stream_tintalgo     = 'null',                           &
         stream_name         = 'dO3 data',                       &
         rc                  = rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
       call ESMF_Finalize(endflag=ESMF_END_ABORT)
    end if

  end subroutine dO3_init

end module diurnalOzoneStreamMod
