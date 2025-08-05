module DistParamMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Spatially distributed parameter data type allocation and initialization
  ! -------------------------------------------------------- 
  !
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
  use shr_sys_mod    , only : shr_sys_abort
  use clm_varcon     , only : spval, ispval, grlnd
  use clm_varctl     , only : iulog
  use clm_varctl     , only : use_distributed_parameters
  use spmdMod        , only : masterproc, mpicom
  use shr_mpi_mod    , only : shr_mpi_bcast
  use decompMod      , only : bounds_type

  !
  ! !PUBLIC TYPES:
  implicit none
  save
  private
  character(len=*), parameter, private :: sourcefile = __FILE__
  !

  public :: InitDistributedParameters

  !------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  subroutine InitDistributedParameters(bounds)
    !
    ! !USES:
    use DistParamType   , only : distributed_parameters, InitGlobalParameters
    use DistParamsStreamMod, only : distributed_parameter_stream
    use clm_varctl      , only : NLFilename => NLFilename_in
    use abortutils      , only : endrun

    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in)      :: bounds
    !
    ! !LOCAL VARIABLES:
    logical            :: readvar               ! whether the variable was found
    character(len=*), parameter :: subname = 'InitDistributedParameters'
    !--------------------------------------------------------------------

    if (masterproc) then
       write(iulog,*) trim(subname)//' :: reading CLM parameters'
    end if

    ! Initialize distributed_parameters object
    call distributed_parameters%Init()

    ! Initialize distributed parameters based on stream data
    call distributed_parameter_stream%Init(bounds, NLFilename, distributed_parameters)

    ! Initialize global parameters
    call InitGlobalParameters(bounds)

  end subroutine InitDistributedParameters

end module DistParamMod
