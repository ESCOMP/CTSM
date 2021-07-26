module abortutils

  !-----------------------------------------------------------------------
  ! !MODULE: abortutils
  !
  ! !DESCRIPTION:
  ! Functions to abort the model for abnormal termination.
  !
  ! Also a related function to write information about a point, since this is often done
  ! in conjunction with aborting the model, or at least issuing a warning.
  !-----------------------------------------------------------------------

  private
  save

  public :: endrun              ! Abort the model for abnormal termination
  public :: write_point_context ! Write context for the given index, including global index information and more

  interface endrun
     module procedure endrun_vanilla
     module procedure endrun_write_point_context
  end interface

  character(len=*), parameter, private :: sourcefile = &
       __FILE__

contains

  !-----------------------------------------------------------------------
  subroutine endrun_vanilla(msg, additional_msg)

    !-----------------------------------------------------------------------
    ! !DESCRIPTION:
    ! Abort the model for abnormal termination
    !
    use shr_sys_mod , only: shr_sys_abort
    use clm_varctl  , only: iulog
    !
    ! !ARGUMENTS:
    implicit none

    ! Generally you want to at least provide msg. The main reason to separate msg from
    ! additional_msg is to supported expected-exception unit testing: you can put
    ! volatile stuff in additional_msg, as in:
    !   call endrun(msg='Informative message', additional_msg=errmsg(__FILE__, __LINE__))
    ! and then just assert against msg.
    character(len=*), intent(in), optional :: msg            ! string to be passed to shr_sys_abort
    character(len=*), intent(in), optional :: additional_msg ! string to be printed, but not passed to shr_sys_abort
    !-----------------------------------------------------------------------

    if (present (additional_msg)) then
       write(iulog,*)'ENDRUN: ', trim(additional_msg)
    else
       write(iulog,*)'ENDRUN:'
    end if

    call shr_sys_abort(msg)

  end subroutine endrun_vanilla

  !-----------------------------------------------------------------------
  subroutine endrun_write_point_context(subgrid_index, subgrid_level, msg, additional_msg)

    !-----------------------------------------------------------------------
    ! Description:
    ! Abort the model for abnormal termination
    !
    ! This version also prints additional information about the point causing the error.
    !
    use shr_sys_mod , only: shr_sys_abort
    use clm_varctl  , only: iulog
    !
    ! Arguments:
    implicit none
    integer          , intent(in)           :: subgrid_index  ! index of interest (can be at any subgrid level or gridcell level)
    character(len=*) , intent(in)           :: subgrid_level  ! one of nameg, namel, namec or namep

    ! Generally you want to at least provide msg. The main reason to separate msg from
    ! additional_msg is to supported expected-exception unit testing: you can put
    ! volatile stuff in additional_msg, as in:
    !   call endrun(msg='Informative message', additional_msg=errmsg(__FILE__, __LINE__))
    ! and then just assert against msg.
    character(len=*), intent(in), optional :: msg            ! string to be passed to shr_sys_abort
    character(len=*), intent(in), optional :: additional_msg ! string to be printed, but not passed to shr_sys_abort
    !
    ! Local Variables:
    integer :: igrc, ilun, icol 
    !-----------------------------------------------------------------------

    call write_point_context(subgrid_index, subgrid_level)

    if (present (additional_msg)) then
       write(iulog,*)'ENDRUN: ', additional_msg
    else
       write(iulog,*)'ENDRUN:'
    end if

    call shr_sys_abort(msg)

  end subroutine endrun_write_point_context

  !-----------------------------------------------------------------------
  subroutine write_point_context(subgrid_index, subgrid_level)

    !-----------------------------------------------------------------------
    ! Description:
    ! Write various information giving context for the given index at the given subgrid
    ! level, including global index information and more.
    !
    use shr_sys_mod  , only : shr_sys_flush, shr_sys_abort
    use shr_log_mod  , only : errMsg => shr_log_errMsg
    use clm_varcon   , only : nameg, namel, namec, namep
    use decompMod    , only : get_global_index
    use GridcellType , only : grc
    use LandunitType , only : lun
    use ColumnType   , only : col
    use PatchType    , only : patch
    use spmdMod      , only : iam
    !
    ! Arguments:
    integer          , intent(in) :: subgrid_index  ! index of interest (can be at any subgrid level or gridcell level)
    character(len=*) , intent(in) :: subgrid_level  ! one of nameg, namel, namec or namep
    !
    ! Local Variables:
    integer :: igrc, ilun, icol, ipft
    !-----------------------------------------------------------------------

    if (trim(subgrid_level) == nameg) then

       igrc = subgrid_index
       write(iulog,'(a, i0, a, i0)') 'iam = ', iam, ': local  gridcell index = ', igrc
       write(iulog,'(a, i0, a, i0)') 'iam = ', iam, ': global gridcell index = ', &
            get_global_index(subgrid_index=igrc, subgrid_level=nameg)
       write(iulog,'(a, i0, a, f12.7)') 'iam = ', iam, ': gridcell longitude    = ', grc%londeg(igrc)
       write(iulog,'(a, i0, a, f12.7)') 'iam = ', iam, ': gridcell latitude     = ', grc%latdeg(igrc)

    else if (trim(subgrid_level) == namel) then

       ilun = subgrid_index
       igrc = lun%gridcell(ilun)
       write(iulog,'(a, i0, a, i0)') 'iam = ', iam, ': local  landunit index = ', ilun
       write(iulog,'(a, i0, a, i0)') 'iam = ', iam, ': global landunit index = ', &
            get_global_index(subgrid_index=ilun, subgrid_level=namel)
       write(iulog,'(a, i0, a, i0)') 'iam = ', iam, ': global gridcell index = ', &
            get_global_index(subgrid_index=igrc, subgrid_level=nameg)
       write(iulog,'(a, i0, a, f12.7)') 'iam = ', iam, ': gridcell longitude    = ', grc%londeg(igrc)
       write(iulog,'(a, i0, a, f12.7)') 'iam = ', iam, ': gridcell latitude     = ', grc%latdeg(igrc)
       write(iulog,'(a, i0, a, i0)') 'iam = ', iam, ': landunit type         = ', lun%itype(subgrid_index)

    else if (trim(subgrid_level) == namec) then

       icol = subgrid_index
       ilun = col%landunit(icol)
       igrc = col%gridcell(icol)
       write(iulog,'(a, i0, a, i0)') 'iam = ', iam, ': local  column   index = ', icol
       write(iulog,'(a, i0, a, i0)') 'iam = ', iam, ': global column   index = ', &
            get_global_index(subgrid_index=icol, subgrid_level=namec)
       write(iulog,'(a, i0, a, i0)') 'iam = ', iam, ': global landunit index = ', &
            get_global_index(subgrid_index=ilun, subgrid_level=namel)
       write(iulog,'(a, i0, a, i0)') 'iam = ', iam, ': global gridcell index = ', &
            get_global_index(subgrid_index=igrc, subgrid_level=nameg)
       write(iulog,'(a, i0, a, f12.7)') 'iam = ', iam, ': gridcell longitude    = ', grc%londeg(igrc)
       write(iulog,'(a, i0, a, f12.7)') 'iam = ', iam, ': gridcell latitude     = ', grc%latdeg(igrc)
       write(iulog,'(a, i0, a, i0)') 'iam = ', iam, ': column   type         = ', col%itype(icol)
       write(iulog,'(a, i0, a, i0)') 'iam = ', iam, ': landunit type         = ', lun%itype(ilun)

    else if (trim(subgrid_level) == namep) then

       ipft = subgrid_index
       icol = patch%column(ipft)
       ilun = patch%landunit(ipft)
       igrc = patch%gridcell(ipft)
       write(iulog,'(a, i0, a, i0)') 'iam = ', iam, ': local  patch    index = ', ipft
       write(iulog,'(a, i0, a, i0)') 'iam = ', iam, ': global patch    index = ', &
            get_global_index(subgrid_index=ipft, subgrid_level=namep)
       write(iulog,'(a, i0, a, i0)') 'iam = ', iam, ': global column   index = ', &
            get_global_index(subgrid_index=icol, subgrid_level=namec)
       write(iulog,'(a, i0, a, i0)') 'iam = ', iam, ': global landunit index = ', &
            get_global_index(subgrid_index=ilun, subgrid_level=namel)
       write(iulog,'(a, i0, a, i0)') 'iam = ', iam, ': global gridcell index = ', &
            get_global_index(subgrid_index=igrc, subgrid_level=nameg)
       write(iulog,'(a, i0, a, f12.7)') 'iam = ', iam, ': gridcell longitude    = ', grc%londeg(igrc)
       write(iulog,'(a, i0, a, f12.7)') 'iam = ', iam, ': gridcell latitude     = ', grc%latdeg(igrc)
       write(iulog,'(a, i0, a, i0)') 'iam = ', iam, ': pft      type         = ', patch%itype(ipft)
       write(iulog,'(a, i0, a, i0)') 'iam = ', iam, ': column   type         = ', col%itype(icol)
       write(iulog,'(a, i0, a, i0)') 'iam = ', iam, ': landunit type         = ', lun%itype(ilun)

    else
       call shr_sys_abort('subgrid_level '//trim(subgrid_level)//'not supported '//errmsg(sourcefile, __LINE__))

    end if

    call shr_sys_flush(iulog)

  end subroutine write_point_context

end module abortutils
