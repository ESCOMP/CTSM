module dynpftFileMod

#include "shr_assert.h"

  !---------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Handle reading of the pftdyn dataset, which specifies transient areas of natural Patches
  !
  ! !USES:
  use shr_kind_mod          , only : r8 => shr_kind_r8
  use shr_log_mod           , only : errMsg => shr_log_errMsg
  use decompMod             , only : bounds_type, BOUNDS_LEVEL_PROC
  use dynFileMod            , only : dyn_file_type
  use dynVarTimeUninterpMod , only : dyn_var_time_uninterp_type
  use clm_varctl            , only : iulog
  use abortutils            , only : endrun
  use spmdMod               , only : masterproc, mpicom
  use clm_varcon            , only : grlnd, nameg
  use LandunitType          , only : lun                
  use ColumnType            , only : col                
  use PatchType             , only : patch                
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  implicit none
  private
  save
  public :: dynpft_init     ! initialize information read from pftdyn dataset
  public :: dynpft_interp   ! interpolate pftdyn information to current time step
  !
  ! !PRIVATE MEMBER FUNCTIONS:
  private :: dynpft_check_consistency   ! check consistency with surface dataset
  private :: dynpft_read_consistency_nl ! read namelist associated with consistency checks
  !
  ! ! PRIVATE TYPES
  type(dyn_file_type), target      :: dynpft_file   ! information for the pftdyn file
  type(dyn_var_time_uninterp_type) :: wtpatch       ! weight of each patch relative to the natural veg landunit

  character(len=*), parameter      :: varname = 'PCT_NAT_PFT'  ! name of variable on file

  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  !---------------------------------------------------------------------------

contains


  !-----------------------------------------------------------------------
  subroutine dynpft_init(bounds, dynpft_filename)
    !
    ! !DESCRIPTION:
    ! Initialize dynamic pft dataset (position it to the right time samples
    ! that bound the initial model date)
    !
    ! !USES:
    use clm_varpar     , only : maxveg, natpft_size
    use ncdio_pio
    use dynTimeInfoMod , only : YEAR_POSITION_START_OF_TIMESTEP
    !
    ! !ARGUMENTS:
    type(bounds_type) , intent(in) :: bounds          ! proc-level bounds
    character(len=*)  , intent(in) :: dynpft_filename ! name of file containing transient pft information
    !
    ! !LOCAL VARIABLES:
    integer  :: wtpatch_shape(2)                  ! shape of the wtpatch data

    character(len= 32)     :: subname='dynpft_init'! subroutine name
    !-----------------------------------------------------------------------

    SHR_ASSERT_ALL(bounds%level == BOUNDS_LEVEL_PROC, subname // ': argument must be PROC-level bounds')

    if (masterproc) then
       write(iulog,*) 'Attempting to read pft dynamic landuse data .....'
    end if

    ! Get the year from the START of the timestep; this way, we'll update PFT areas
    ! starting after the year boundary. This is consistent with the timing of other area
    ! updates.
    dynpft_file = dyn_file_type(dynpft_filename, YEAR_POSITION_START_OF_TIMESTEP)

    ! Consistency checks
    call check_dim(dynpft_file, 'natpft', natpft_size)
    call dynpft_check_consistency(bounds)

    ! read data PCT_NAT_PFT corresponding to correct year
    !
    ! Note: if you want to change PCT_NAT_PFT so that it is interpolated, rather than
    ! jumping to each year's value on Jan 1 of that year, simply change wtpatch to be of type
    ! dyn_var_time_interp_type (rather than dyn_var_time_uninterp_type), and change the
    ! following constructor to construct a variable of dyn_var_time_interp_type. That's
    ! all you need to do.

    wtpatch_shape = [(bounds%endg-bounds%begg+1), natpft_size]
    wtpatch = dyn_var_time_uninterp_type( &
         dyn_file=dynpft_file, varname=varname, &
         dim1name=grlnd, conversion_factor=100._r8, &
         do_check_sums_equal_1=.true., data_shape=wtpatch_shape)

  end subroutine dynpft_init

  !-----------------------------------------------------------------------
  subroutine dynpft_check_consistency(bounds)
    !
    ! !DESCRIPTION:
    ! Check consistency between dynpft file and surface dataset.
    !
    ! This is done by assuming that PCT_NAT_PFT at time 1 in the pftdyn file agrees with
    ! PCT_NAT_PFT on the surface dataset.
    !
    ! !USES:
    use clm_instur, only : wt_nat_patch
    use clm_varpar, only : natpft_size
    use ncdio_pio
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds  ! proc-level bounds
    !
    ! !LOCAL VARIABLES:
    logical             :: check_dynpft_consistency ! whether to do the consistency check in this routine
    integer             :: g                        ! index
    real(r8), pointer   :: wtpatch_time1(:,:)       ! weight of each pft in each grid cell at first time
    logical             :: readvar                  ! whether variable was read
    real(r8), parameter :: tol = 1.e-6_r8           ! tolerance for checking equality

    character(len=*), parameter :: subname = 'dynpft_check_consistency'
    !-----------------------------------------------------------------------
    
    call dynpft_read_consistency_nl(check_dynpft_consistency)

    if (check_dynpft_consistency) then

       ! Read first time slice of PCT_NAT_PFT

       allocate(wtpatch_time1(bounds%begg:bounds%endg, natpft_size))
       call ncd_io(ncid=dynpft_file, varname=varname, flag='read', data=wtpatch_time1, &
            dim1name=grlnd, nt=1, readvar=readvar)
       if (.not. readvar) then
          call endrun(msg=' ERROR: ' // trim(varname) // ' NOT on landuse_timeseries file'//&
               errMsg(sourcefile, __LINE__))
       end if

       ! Convert from PCT to weight on grid cell
       wtpatch_time1(bounds%begg:bounds%endg,:) = wtpatch_time1(bounds%begg:bounds%endg,:) / 100._r8
    
       ! Compare with values read from surface dataset
       do g = bounds%begg, bounds%endg
          if (any(abs(wtpatch_time1(g,:) - wt_nat_patch(g,:)) > tol)) then
             write(iulog,*) subname//' mismatch between PCT_NAT_PFT at initial time and that obtained from surface dataset'
             write(iulog,*) 'On landuse_timeseries file: ', wtpatch_time1(g,:)
             write(iulog,*) 'On surface dataset: ', wt_nat_patch(g,:)
             write(iulog,*) ' '
             write(iulog,*) 'Confirm that the year of your surface dataset'
             write(iulog,*) 'corresponds to the first year of your landuse_timeseries file'
             write(iulog,*) '(e.g., for a landuse_timeseries file starting at year 1850, which is typical,'
             write(iulog,*) 'you should be using an 1850 surface dataset),'
             write(iulog,*) 'and that your landuse_timeseries file is compatible with the surface dataset.'
             write(iulog,*) ' '
             write(iulog,*) 'If you are confident that you are using the correct landuse_timeseries file'
             write(iulog,*) 'and the correct surface dataset, then you can bypass this check by setting:'
             write(iulog,*) '  check_dynpft_consistency = .false.'
             write(iulog,*) 'in user_nl_clm'
             write(iulog,*) ' '
             call endrun(decomp_index=g, clmlevel=nameg, msg=errMsg(sourcefile, __LINE__))
          end if
       end do

       deallocate(wtpatch_time1)

    end if

  end subroutine dynpft_check_consistency

  !-----------------------------------------------------------------------
  subroutine dynpft_read_consistency_nl(check_dynpft_consistency)
    !
    ! !DESCRIPTION:
    ! Read namelist settings related to pftdyn consistency checks
    !
    ! !USES:
    use fileutils      , only : getavu, relavu
    use clm_nlUtilsMod , only : find_nlgroup_name
    use controlMod     , only : NLFilename
    use shr_mpi_mod    , only : shr_mpi_bcast
    !
    ! !ARGUMENTS:
    logical, intent(out) :: check_dynpft_consistency ! whether to do the consistency check
    !
    ! !LOCAL VARIABLES:
    integer :: nu_nml    ! unit for namelist file
    integer :: nml_error ! namelist i/o error flag

    character(len=*), parameter :: subname = 'dynpft_read_consistency_nl'
    !-----------------------------------------------------------------------

    namelist /dynpft_consistency_checks/ &
         check_dynpft_consistency

    ! Set default namelist values
    check_dynpft_consistency = .true.

    ! Read namelist
    if (masterproc) then
       nu_nml = getavu()
       open( nu_nml, file=trim(NLFilename), status='old', iostat=nml_error )
       call find_nlgroup_name(nu_nml, 'dynpft_consistency_checks', status=nml_error)
       if (nml_error == 0) then
          read(nu_nml, nml=dynpft_consistency_checks,iostat=nml_error)
          if (nml_error /= 0) then
             call endrun(msg='ERROR reading dynpft_consistency_checks namelist'//errMsg(sourcefile, __LINE__))
          end if
       else
          call endrun(msg='ERROR finding dynpft_consistency_checks namelist'//errMsg(sourcefile, __LINE__))
       end if
       close(nu_nml)
       call relavu( nu_nml )
    endif

    call shr_mpi_bcast (check_dynpft_consistency, mpicom)

    if (masterproc) then
       write(iulog,*) ' '
       write(iulog,*) 'dynpft_consistency_checks settings:'
       write(iulog,nml=dynpft_consistency_checks)
       write(iulog,*) ' '
    end if

  end subroutine dynpft_read_consistency_nl



  !-----------------------------------------------------------------------
  subroutine dynpft_interp(bounds)
    !
    ! !DESCRIPTION:
    ! Get pft weights for current model time
    !
    ! Sets pft%wtcol
    !
    ! Note that PFT weights currently jump to their new value at the start of the year.
    ! However, as mentioned above, this behavior can be changed to time interpolation
    ! simply by making wtpatch a dyn_var_time_interp_type variable rather than
    ! dyn_var_time_uninterp_type.
    !
    ! !USES:
    use landunit_varcon , only : istsoil
    use clm_varpar      , only : natpft_lb, natpft_ub
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds  ! proc-level bounds
    !
    ! !LOCAL VARIABLES:
    integer               :: m,p,l,g          ! indices
    real(r8), allocatable :: wtpatch_cur(:,:)   ! current pft weights
    character(len=32) :: subname='dynpft_interp' ! subroutine name
    !-----------------------------------------------------------------------

    ! assumes that each landunit has only 1 column, 
    ! and SCAM and CNDV have not been defined
    !
    ! NOTE(wjs, 2014-12-10) I'm not sure if there is still the requirement that SCAM
    ! hasn't been defined

    SHR_ASSERT_ALL(bounds%level == BOUNDS_LEVEL_PROC, subname // ': argument must be PROC-level bounds')

    ! Get pft weights for this time step

    call dynpft_file%time_info%set_current_year()

    allocate(wtpatch_cur(bounds%begg:bounds%endg, natpft_lb:natpft_ub))
    call wtpatch%get_current_data(wtpatch_cur)

    do p = bounds%begp,bounds%endp
       g = patch%gridcell(p)
       l = patch%landunit(p)

       if (lun%itype(l) == istsoil) then
          m = patch%itype(p)

          ! Note that the following assignment assumes that all Patches share a single column
          patch%wtcol(p) = wtpatch_cur(g,m)
       end if

    end do

    deallocate(wtpatch_cur)

  end subroutine dynpft_interp

end module dynpftFileMod
