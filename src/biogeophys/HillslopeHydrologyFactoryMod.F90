module HillslopeHydrologyFactoryMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Calculate geomorphological quantities for hillslope columns.
  !
  ! !USES:
#include "shr_assert.h"
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use shr_log_mod    , only : errMsg => shr_log_errMsg
  use spmdMod        , only : masterproc
  use abortutils     , only : endrun
  use clm_varctl     , only : iulog, use_hillslope
  use decompMod      , only : bounds_type
  implicit none
  private   
  save

! namelist variable specifying geomorphology defining equation set
  character(len=256) :: hillslope_geomorphology   

  !-----------------------------------------------------------------------

  ! !PUBLIC ROUTINES
  public  :: create_and_init_hillslope_geomorphology_type

  ! !PRIVATE ROUTINES
  private :: hillslope_geomorphology_readNL

contains
     
  function create_and_init_hillslope_geomorphology_type(bounds) result(hg)
    !
    ! !DESCRIPTION:
    ! Create and initialize an object of hillslope_geomorphology_type, 
    ! and return this object. The particular type is determined based on 
    ! the use_hillslope namelist parameter.
       !
       ! !USES:
    use controlMod                       , only : NLFilename
    use clm_varctl                       , only : use_hillslope, fsurdat
    use HillslopeHydrologyBaseMod        , only : hillslope_geomorphology_type
    use HillslopeHydrologyIndependentMod , only : hillslope_geomorphology_independent_type
    use HillslopeHydrologyTroch02Mod     , only : hillslope_geomorphology_troch02_type

    !
    ! !ARGUMENTS:
    class(hillslope_geomorphology_type), allocatable :: hg ! function result
    type(bounds_type), intent(in) :: bounds

    !-----------------------------------------------------------------

    call hillslope_geomorphology_readNL(NLFilename)

    select case (trim(hillslope_geomorphology))
    case ('independent')
       allocate(hg, source = hillslope_geomorphology_independent_type())
    case ('troch02')
       allocate(hg, source = hillslope_geomorphology_troch02_type())

    end select

    call hg%Init(bounds, fsurdat)

  end function create_and_init_hillslope_geomorphology_type
  
  !-----------------------------------------------------------------------
  subroutine hillslope_geomorphology_readNL(NLFilename)
   !
   !DESCRIPTIONS
   ! Read the namelist for soil resistance method
    !
    ! !USES:
    use abortutils      , only : endrun   
    use fileutils       , only : getavu, relavu
    use spmdMod         , only : mpicom, masterproc
    use shr_mpi_mod     , only : shr_mpi_bcast
    use clm_varctl      , only : iulog
    use clm_nlUtilsMod  , only : find_nlgroup_name

    ! !ARGUMENTS:
    !------------------------------------------------------------------------------
    implicit none
    character(len=*), intent(IN) :: NLFilename ! Namelist filename
    integer                      :: nu_nml     ! unit for namelist file
    integer                      :: nml_error  ! namelist i/o error flag
    character(*), parameter      :: subName = "('hillslope_geomorphology_readNL')"

    !-----------------------------------------------------------------------

! MUST agree with name in namelist and read statement
    namelist /hillslope_hydrology_inparm/ hillslope_geomorphology

    ! Default values for namelist

   hillslope_geomorphology = 'independent'

    ! Read soil_resis namelist
    if (masterproc) then
       nu_nml = getavu()
       open( nu_nml, file=trim(NLFilename), status='old', iostat=nml_error )
       call find_nlgroup_name(nu_nml, 'hillslope_hydrology_inparm', status=nml_error)
       if (nml_error == 0) then
          read(nu_nml, nml=hillslope_hydrology_inparm,iostat=nml_error)
          if (nml_error /= 0) then
             call endrun(subname // ':: ERROR reading hillslope_hydrology namelist')
          end if
       else
          call endrun(subname // ':: ERROR reading hillslope_hydrology namelist')
       end if
       close(nu_nml)
       call relavu( nu_nml )

    endif

    call shr_mpi_bcast(hillslope_geomorphology, mpicom)

    if (masterproc) then
       write(iulog,*) ' '
       write(iulog,*) 'hillslope hydrology settings:'
       write(iulog,*) '  hillslope_geomorphology  = ',hillslope_geomorphology
    endif

  end subroutine hillslope_geomorphology_readNL

end module HillslopeHydrologyFactoryMod
