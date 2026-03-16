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

  use shr_log_mod  , only : errMsg => shr_log_errMsg

  implicit none
  private

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
  subroutine endrun_vanilla(msg, additional_msg, line, file)

    !-----------------------------------------------------------------------
    ! !DESCRIPTION:
    ! Abort the model for abnormal termination
    !
    use shr_abort_mod , only: shr_abort_abort
    use clm_varctl  , only: iulog
    !
    ! !ARGUMENTS:
    ! Generally you want to at least provide msg. The main reason to separate msg from
    ! additional_msg is to supported expected-exception unit testing: you can put
    ! volatile stuff in additional_msg, as in:
    !   call endrun(msg='Informative message', additional_msg=datetime )
    ! and then just assert against msg.
    character(len=*), intent(in), optional :: msg            ! string to be passed to shr_abort
    character(len=*), intent(in), optional :: additional_msg ! string to be printed, but not passed to shr_abort
    integer , intent(in), optional :: line                   ! Line number for the endrun call
    character(len=*), intent(in), optional :: file           ! file for the endrun call
    !-----------------------------------------------------------------------

    if (present (additional_msg)) then
       write(iulog,*)'ENDRUN: ', trim(additional_msg)
    else
       write(iulog,*)'ENDRUN:'
    end if

    ! Don't pass file and line to shr_abort_abort since the PFUNIT test version doesn't have those options
    if ( present(file) .and. present(line) ) then
       write(iulog,*) errMsg(file, line)
    end if
    call shr_abort_abort(string=msg)

  end subroutine endrun_vanilla

  !-----------------------------------------------------------------------
  subroutine endrun_write_point_context(subgrid_index, subgrid_level, msg, additional_msg, line, file)

    !-----------------------------------------------------------------------
    ! Description:
    ! Abort the model for abnormal termination
    !
    ! This version also prints additional information about the point causing the error.
    !
    use clm_varctl  , only: iulog
    use decompMod   , only: subgrid_level_unspecified
    !
    ! Arguments:
    integer , intent(in) :: subgrid_index ! index of interest (can be at any subgrid level or gridcell level)
    integer , intent(in) :: subgrid_level ! one of the subgrid_level_* constants defined in decompMod; subgrid_level_unspecified is allowed here, in which case the additional information will not be printed
    integer , intent(in), optional :: line                   ! Line number for the endrun call
    character(len=*), intent(in), optional :: file           !file for the endrun call
    character(len=*), intent(in), optional :: msg            ! string to be passed to shr_abort
    character(len=*), intent(in), optional :: additional_msg ! string to be printed, but not passed to shr_abort
    !
    ! Local Variables:
    integer :: igrc, ilun, icol 
    !-----------------------------------------------------------------------

    if (subgrid_level /= subgrid_level_unspecified) then
       call write_point_context(subgrid_index, subgrid_level)
    end if

    call endrun_vanilla(msg=msg, additional_msg=additional_msg, line=line, file=file)

  end subroutine endrun_write_point_context

  !-----------------------------------------------------------------------
  subroutine write_point_context(subgrid_index, subgrid_level)

    !-----------------------------------------------------------------------
    ! Description:
    ! Write various information giving context for the given index at the given subgrid
    ! level, including global index information and more.
    !
    ! NOTE: DO NOT CALL AN ABORT FROM HERE AS THAT WOULD SHORT CIRUIT THE ERROR REPORTING!!
    !
    use clm_varctl   , only : iulog
    use decompMod    , only : subgrid_level_gridcell, subgrid_level_landunit, subgrid_level_column, subgrid_level_patch
    use decompMod    , only : get_global_index
    use GridcellType , only : grc
    use LandunitType , only : lun
    use ColumnType   , only : col
    use PatchType    , only : patch
    use spmdMod      , only : iam
    !
    ! Arguments:
    integer , intent(in) :: subgrid_index ! index of interest (can be at any subgrid level or gridcell level)
    integer , intent(in) :: subgrid_level ! one of the subgrid_level_* constants defined in decompMod
    !
    ! Local Variables:
    integer, parameter :: unset = -9999 ! Unset value for an index
    integer :: igrc=unset, ilun=unset, icol=unset, ipft=unset   ! Local index for grid-cell, landunit, column, and patch
    integer :: ggrc=unset, glun=unset, gcol=unset, gpft=unset   ! Global index for grid-cell, landunit, column, and patch
    logical :: bad_point = .false. ! Flag to indicate if the point is bad (i.e., global index is -1)
    !-----------------------------------------------------------------------

    if (subgrid_level == subgrid_level_gridcell) then

       igrc = subgrid_index
       ggrc = get_global_index(subgrid_index=igrc, subgrid_level=subgrid_level_gridcell, donot_abort_on_badindex=.true.)

    else if (subgrid_level == subgrid_level_landunit) then

       ilun = subgrid_index
       glun = get_global_index(subgrid_index=ilun, subgrid_level=subgrid_level_landunit, donot_abort_on_badindex=.true.)
       if ( glun /= -1 ) then
          igrc = lun%gridcell(ilun)
          ggrc = get_global_index(subgrid_index=igrc, subgrid_level=subgrid_level_gridcell, donot_abort_on_badindex=.true.)
       else
          bad_point = .true.
       end if

    else if (subgrid_level == subgrid_level_column) then

       icol = subgrid_index
       gcol = get_global_index(subgrid_index=icol, subgrid_level=subgrid_level_column, donot_abort_on_badindex=.true.)
       if ( gcol /= -1 ) then
          ilun = col%landunit(icol)
          igrc = col%gridcell(icol)
          ggrc = get_global_index(subgrid_index=igrc, subgrid_level=subgrid_level_gridcell, donot_abort_on_badindex=.true.)
          glun = get_global_index(subgrid_index=ilun, subgrid_level=subgrid_level_landunit, donot_abort_on_badindex=.true.)
       else
          bad_point = .true.
       end if

    else if (subgrid_level == subgrid_level_patch) then

       ipft = subgrid_index
       gpft = get_global_index(subgrid_index=ipft, subgrid_level=subgrid_level_patch, donot_abort_on_badindex=.true.)
       if ( gpft /= -1 ) then
          icol = patch%column(ipft)
          ilun = patch%landunit(ipft)
          igrc = patch%gridcell(ipft)
          ggrc = get_global_index(subgrid_index=igrc, subgrid_level=subgrid_level_gridcell, donot_abort_on_badindex=.true.)
          glun = get_global_index(subgrid_index=ilun, subgrid_level=subgrid_level_landunit, donot_abort_on_badindex=.true.)
          gcol = get_global_index(subgrid_index=icol, subgrid_level=subgrid_level_column, donot_abort_on_badindex=.true.)
       else
          bad_point = .true.
       end if

    end if

    !
    ! Badpoint should already be determined, but check again in case one of the subsequent
    ! calls to get_global_index returns -1
    ! If one of the global indices is -1 then this is a bad point, so flag a bad-point
    if ( igrc /= unset) then
      if ( ggrc == -1 ) bad_point = .true.
    end if
    if ( ilun /= unset) then
      if ( glun == -1 ) bad_point = .true.
    end if
    if ( icol /= unset) then
      if ( gcol == -1 ) bad_point = .true.
    end if
    if ( ipft /= unset) then
      if ( gpft == -1 ) bad_point = .true.
    end if

    if (bad_point) then
       write(iulog,*) 'A bad input point was given: subgrid_index = ', subgrid_index, &
            ', subgrid_level = ', subgrid_level
       write(iulog,*) errMsg(sourcefile, __LINE__)
       write(iulog,*) 'Continuing the endrun without writing point context information'
       return
    end if

    if (subgrid_level == subgrid_level_gridcell) then

       write(iulog,'(a, i0, a, i0)') 'iam = ', iam, ': local  gridcell index = ', igrc
       write(iulog,'(a, i0, a, i0)') 'iam = ', iam, ': global gridcell index = ', ggrc
       write(iulog,'(a, i0, a, f12.7)') 'iam = ', iam, ': gridcell longitude    = ', grc%londeg(igrc)
       write(iulog,'(a, i0, a, f12.7)') 'iam = ', iam, ': gridcell latitude     = ', grc%latdeg(igrc)

    else if (subgrid_level == subgrid_level_landunit) then

       write(iulog,'(a, i0, a, i0)') 'iam = ', iam, ': local  landunit index = ', ilun
       write(iulog,'(a, i0, a, i0)') 'iam = ', iam, ': global landunit index = ', glun
       write(iulog,'(a, i0, a, i0)') 'iam = ', iam, ': global gridcell index = ', ggrc
       write(iulog,'(a, i0, a, f12.7)') 'iam = ', iam, ': gridcell longitude    = ', grc%londeg(igrc)
       write(iulog,'(a, i0, a, f12.7)') 'iam = ', iam, ': gridcell latitude     = ', grc%latdeg(igrc)
       write(iulog,'(a, i0, a, i0)') 'iam = ', iam, ': landunit type         = ', lun%itype(subgrid_index)

    else if (subgrid_level == subgrid_level_column) then

       write(iulog,'(a, i0, a, i0)') 'iam = ', iam, ': local  column   index = ', icol
       write(iulog,'(a, i0, a, i0)') 'iam = ', iam, ': global column   index = ', gcol
       write(iulog,'(a, i0, a, i0)') 'iam = ', iam, ': global landunit index = ', glun
       write(iulog,'(a, i0, a, i0)') 'iam = ', iam, ': global gridcell index = ', ggrc
       write(iulog,'(a, i0, a, f12.7)') 'iam = ', iam, ': gridcell longitude    = ', grc%londeg(igrc)
       write(iulog,'(a, i0, a, f12.7)') 'iam = ', iam, ': gridcell latitude     = ', grc%latdeg(igrc)
       write(iulog,'(a, i0, a, i0)') 'iam = ', iam, ': column   type         = ', col%itype(icol)
       write(iulog,'(a, i0, a, i0)') 'iam = ', iam, ': landunit type         = ', lun%itype(ilun)

    else if (subgrid_level == subgrid_level_patch) then

       write(iulog,'(a, i0, a, i0)') 'iam = ', iam, ': local  patch    index = ', ipft
       write(iulog,'(a, i0, a, i0)') 'iam = ', iam, ': global patch    index = ', gpft
       write(iulog,'(a, i0, a, i0)') 'iam = ', iam, ': global column   index = ', gcol
       write(iulog,'(a, i0, a, i0)') 'iam = ', iam, ': global landunit index = ', glun
       write(iulog,'(a, i0, a, i0)') 'iam = ', iam, ': global gridcell index = ', ggrc
       write(iulog,'(a, i0, a, f12.7)') 'iam = ', iam, ': gridcell longitude    = ', grc%londeg(igrc)
       write(iulog,'(a, i0, a, f12.7)') 'iam = ', iam, ': gridcell latitude     = ', grc%latdeg(igrc)
       write(iulog,'(a, i0, a, i0)') 'iam = ', iam, ': pft      type         = ', patch%itype(ipft)
       write(iulog,'(a, i0, a, i0)') 'iam = ', iam, ': column   type         = ', col%itype(icol)
       write(iulog,'(a, i0, a, i0)') 'iam = ', iam, ': landunit type         = ', lun%itype(ilun)

    else
       write(iulog,*) 'subgrid_level not supported: ', subgrid_level
       write(iulog,*) errMsg(sourcefile, __LINE__)
       write(iulog,*) 'Continuing the endrun without writing point context information'
       return
    end if

  end subroutine write_point_context

end module abortutils
