module CropPoolsMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! This module contains variables giving information about the different crop-specific
  ! pools (e.g., the number of reproductive grain and structure pools), and subroutines
  ! for working with these variables.
  !
  ! !USES:
  !
  ! !PUBLIC TYPES:
  implicit none
  private

  integer, public :: ngrain ! number of crop reproductive grain pools
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: crop_pools_init  ! initialize module-level data

  character(len=*), parameter, private :: sourcefile = &
       __FILE__

contains
  !-----------------------------------------------------------------------
  subroutine crop_pools_init()
    !
    ! !DESCRIPTION:
    ! Initialize module-level data
    !
    ! !ARGUMENTS:
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'crop_pools_init'
    !-----------------------------------------------------------------------

    ! TODO(wjs, 2022-02-14) This will eventually depend on whether we're using AgSys
    ngrain = 1

  end subroutine crop_pools_init

end module CropPoolsMod
