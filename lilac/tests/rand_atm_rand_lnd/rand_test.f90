module rand_test

  use lilac, only : LilacType
  use ESMF

  implicit none

  character(*), parameter :: modname = "(rand_test)"
  integer, parameter      :: num_timesteps = 10
  type(LilacType), save   :: lilac_obj

  !--------------------------------------------------------------------------
  ! Public interfaces
  !--------------------------------------------------------------------------
  public :: atm_driver

  private  :: atm_init 
  private  :: lnd_init 
  private  :: atm_run 
  private  :: lnd_run 
  private  :: atm_final 
  private  :: lnd_final

contains

  subroutine atm_driver(rc)

    integer, intent(out)  :: rc


    call atm_init(rc) 
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    call atm_run(rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return


    call atm_final(rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

  end subroutine atm_driver

  subroutine atm_init(rc)

    integer, intent(out)  :: rc

    ! Initialize atmosphere
    ! TODO

    ! Initialize land via lilac
    call lnd_init(rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

  end subroutine atm_init


  subroutine lnd_init(rc)

    integer, intent(out)  :: rc
    character(len=ESMF_MAXSTR), parameter :: lilac_name="lilac_rand_test"

    call lilac_obj%init(lilac_name)


  end subroutine lnd_init


  subroutine atm_run(rc)

    integer, intent(out)  :: rc

    integer :: n

    ! Run atm for num_timesteps
    do n = 1, num_timesteps, 1
       print *, "------ Running land -------"
       ! Run land via lilac
       call lnd_run(rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    end do

  end subroutine atm_run


  subroutine lnd_run(rc)

    integer, intent(out)  :: rc

    call lilac_obj%run()


  end subroutine lnd_run


  subroutine atm_final(rc)

    integer, intent(out)  :: rc


  end subroutine atm_final


  subroutine lnd_final(rc)

    integer, intent(out)  :: rc

    call lilac_obj%final()
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

  end subroutine lnd_final

end module rand_test
