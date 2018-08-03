module lilac

  implicit none

  !--------------------------------------------------------------------------
  ! Public interfaces
  !--------------------------------------------------------------------------
  public :: lilac_init
  public :: lilac_run
  public :: lilac_final

contains

  subroutine lilac_init()
    implicit none
    print *, "lilac_init()"
    flush(6)

  end subroutine lilac_init

  subroutine lilac_run()
    implicit none
    print *, "lilac_run()"
  end subroutine lilac_run

  subroutine lilac_final()
    implicit none
    print *, "lilac_final()"
  end subroutine lilac_final

end module lilac
