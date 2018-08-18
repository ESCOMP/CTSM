module lilac

#include "ESMF.h"
  use ESMF

  use atmos_comp, only : atmos_setvm, atmos_register
  use land_comp,  only : land_setvm,  land_register
  use coupler_comp, only : usercpl_setvm, usercpl_register

  implicit none

  !--------------------------------------------------------------------------
  ! Public interfaces
  !--------------------------------------------------------------------------
  public :: lilac_init
  public :: lilac_run
  public :: lilac_final

  type(LilacType), save         :: lilac_obj

contains

  type, public :: LilacType
    private

    type(ESMFInfoType)             :: esmf_info

  contains
    procedure, public  :: init  => init
    procedure, public  :: run => run
    procedure, public  :: final => final
  end type LilacType

contains

  subroutine lilac_init(self)
    implicit none

    print *, "lilac_init()"

    ! Initialize ESMF structures
    call self%esmf_info%init("lilac")

  end subroutine lilac_init

  subroutine lilac_run(self)
    implicit none

    call self%esmf_info%run()

  end subroutine lilac_run

  subroutine lilac_final(self)
    implicit none

    call self%esmf_info%final()

  end subroutine lilac_final

end module lilac
