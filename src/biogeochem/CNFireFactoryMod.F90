module CNFireFactoryMod

  !---------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Factory to create an instance of fire_method_type. This module figures
  ! out the particular type to return.
  !
  ! !USES:
  use abortutils          , only : endrun
  use shr_log_mod         , only : errMsg => shr_log_errMsg
  use clm_varctl          , only : iulog

  implicit none
  save
  private
  !
  ! !PUBLIC ROUTINES:
  public :: CNFireReadNML         ! read the fire namelist
  public :: create_cnfire_method  ! create an object of class fire_method_type

  ! !PRIVATE DATA MEMBERS:
  character(len=80), private :: fire_method = "li2014qianfrc"

  character(len=*), parameter, private :: sourcefile = &
       __FILE__

contains

  !-----------------------------------------------------------------------
  subroutine CNFireReadNML( NLFilename )
    !
    ! !DESCRIPTION:
    ! Read the namelist for cnfire
    !
    ! !USES:
    use fileutils      , only : getavu, relavu, opnfil
    use shr_nl_mod     , only : shr_nl_find_group_name
    use spmdMod        , only : masterproc, mpicom
    use shr_mpi_mod    , only : shr_mpi_bcast
    !
    ! !ARGUMENTS:
    character(len=*), intent(in) :: NLFilename ! Namelist filename
    !
    ! !LOCAL VARIABLES:
    integer :: ierr                 ! error code
    integer :: unitn                ! unit for namelist file

    character(len=*), parameter :: subname = 'CNFireReadNML'
    character(len=*), parameter :: nmlname = 'cnfire_inparm'
    !-----------------------------------------------------------------------

    namelist /cnfire_inparm/ fire_method

    ! Initialize options to default values, in case they are not specified in
    ! the namelist

    if (masterproc) then
       unitn = getavu()
       write(iulog,*) 'Read in '//nmlname//'  namelist'
       call opnfil (NLFilename, unitn, 'F')
       call shr_nl_find_group_name(unitn, nmlname, status=ierr)
       if (ierr == 0) then
          read(unitn, nml=cnfire_inparm, iostat=ierr)
          if (ierr /= 0) then
             call endrun(msg="ERROR reading "//nmlname//"namelist"//errmsg(sourcefile, __LINE__))
          end if
       else
          call endrun(msg="ERROR finding "//nmlname//"namelist"//errmsg(sourcefile, __LINE__))
       end if
       call relavu( unitn )
    end if

    call shr_mpi_bcast (fire_method, mpicom)

    if (masterproc) then
       write(iulog,*) ' '
       write(iulog,*) nmlname//' settings:'
       write(iulog,nml=cnfire_inparm)
       write(iulog,*) ' '
    end if
  end subroutine CNFireReadNML
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  subroutine create_cnfire_method( NLFilename, cnfire_method )
    !
    ! !DESCRIPTION:
    ! Create and return an object of fire_method_type. The particular type
    ! is determined based on a namelist parameter.
    !
    ! !USES:
    use FireMethodType   , only : fire_method_type
    use CNFireNoFireMod  , only : cnfire_nofire_type
    use CNFireLi2014Mod  , only : cnfire_li2014_type
    use CNFireLi2016Mod  , only : cnfire_li2016_type
    use CNFireLi2021Mod  , only : cnfire_li2021_type
    use CNFireLi2024Mod  , only : cnfire_li2024_type
    use decompMod        , only : bounds_type
    !
    ! !ARGUMENTS:
    character(len=*), intent(in) :: NLFilename ! Namelist filename
    class(fire_method_type), allocatable, intent(inout) :: cnfire_method
    !
    ! !LOCAL VARIABLES:
    character(len=*), parameter :: subname = 'create_cnfire_method'
    !-----------------------------------------------------------------------
    
    select case (trim(fire_method))
       
    case ("nofire")
       allocate(cnfire_nofire_type :: cnfire_method)
    case ("li2014qianfrc")
       allocate(cnfire_li2014_type :: cnfire_method)
    case ("li2016crufrc")
       allocate(cnfire_li2016_type :: cnfire_method)
    case ("li2021gswpfrc")
       allocate(cnfire_li2021_type :: cnfire_method)
    case ("li2024gswpfrc", "li2024crujra")
       allocate(cnfire_li2024_type :: cnfire_method)

    case default
       write(iulog,*) subname//' ERROR: unknown method: ', fire_method
       call endrun(msg=errMsg(sourcefile, __LINE__))

    end select
    call cnfire_method%FireReadNML( NLFilename )

  end subroutine create_cnfire_method
  !-----------------------------------------------------------------------

end module CNFireFactoryMod
