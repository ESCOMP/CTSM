module GetGlobalValuesMod

  !-----------------------------------------------------------------------
  ! Obtain and Write Global Index information
  !-----------------------------------------------------------------------
  implicit none
  private

  ! PUBLIC MEMBER FUNCTIONS:

  public :: get_global_index
  public :: get_global_index_array
  public :: write_point_context

  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  integer function get_global_index(subgrid_index, subgrid_level)

    !----------------------------------------------------------------
    ! Description
    ! Determine global index space value for target point at given subgrid level
    !
    ! Uses:
    use shr_log_mod, only: errMsg => shr_log_errMsg
    use decompMod  , only: bounds_type, get_subgrid_level_gindex, get_proc_bounds
    use spmdMod    , only: iam
    use clm_varcon , only: nameg, namel, namec, namep
    use clm_varctl , only: iulog
    use shr_sys_mod, only: shr_sys_abort
    !
    ! Arguments 
    integer          , intent(in) :: subgrid_index  ! index of interest (can be at any subgrid level or gridcell level)
    character(len=*) , intent(in) :: subgrid_level  ! one of nameg, namel, namec or namep
    !
    ! Local Variables:
    type(bounds_type) :: bounds_proc   ! processor bounds
    integer           :: beg_index     ! beginning proc index for subgrid_level
    integer, pointer  :: gindex(:)
    !----------------------------------------------------------------

    call get_proc_bounds(bounds_proc, allow_call_from_threaded_region=.true.)

    if (trim(subgrid_level) == nameg) then
       beg_index = bounds_proc%begg
    else if (trim(subgrid_level) == namel) then
       beg_index = bounds_proc%begl
    else if (trim(subgrid_level) == namec) then
       beg_index = bounds_proc%begc
    else if (trim(subgrid_level) == namep) then
       beg_index = bounds_proc%begp
    else
       call shr_sys_abort('subgrid_level of '//trim(subgrid_level)//' not supported' // &
            errmsg(sourcefile, __LINE__))
    end if

    call get_subgrid_level_gindex(subgrid_level=trim(subgrid_level), gindex=gindex)
    get_global_index = gindex(subgrid_index - beg_index + 1)

  end function get_global_index

  !-----------------------------------------------------------------------
  function get_global_index_array(subgrid_index, bounds1, bounds2, subgrid_level)

    !----------------------------------------------------------------
    ! Description
    ! Determine global index space value for target array at given subgrid level
    !
    ! Example from histFileMod.F90:
    ! ilarr = get_global_index_array(lun%gridcell(bounds%begl:bounds%endl), bounds%begl, bounds%endl, subgrid_level=nameg)
    ! Note that the last argument (subgrid_level) is set to nameg, which corresponds
    ! to the "gridcell" not the "lun" of the first argument.
    !
    ! Uses:
#include "shr_assert.h"
    use shr_log_mod, only: errMsg => shr_log_errMsg
    use decompMod  , only: bounds_type, get_subgrid_level_gindex, get_proc_bounds
    use spmdMod    , only: iam
    use clm_varcon , only: nameg, namel, namec, namep
    use clm_varctl , only: iulog
    use shr_sys_mod, only: shr_sys_abort
    !
    ! Arguments 
    integer          , intent(in) :: bounds1  ! lower bound of the input & returned arrays
    integer          , intent(in) :: bounds2  ! upper bound of the input & returned arrays
    integer          , intent(in) :: subgrid_index(bounds1:)  ! array of indices of interest (can be at any subgrid level or gridcell level)
    character(len=*) , intent(in) :: subgrid_level  ! one of nameg, namel, namec or namep
    integer                       :: get_global_index_array(bounds1:bounds2)
    !
    ! Local Variables:
    type(bounds_type) :: bounds_proc   ! processor bounds
    integer           :: beg_index     ! beginning proc index for subgrid_level
    integer           :: i
    integer , pointer :: gindex(:)
    !----------------------------------------------------------------

    SHR_ASSERT_ALL_FL((ubound(subgrid_index) == (/bounds2/)), sourcefile, __LINE__)
    call get_proc_bounds(bounds_proc, allow_call_from_threaded_region=.true.)

    if (trim(subgrid_level) == nameg) then
       beg_index = bounds_proc%begg
    else if (trim(subgrid_level) == namel) then
       beg_index = bounds_proc%begl
    else if (trim(subgrid_level) == namec) then
       beg_index = bounds_proc%begc
    else if (trim(subgrid_level) == namep) then
       beg_index = bounds_proc%begp
    else
       call shr_sys_abort('subgrid_level of '//trim(subgrid_level)//' not supported' // &
            errmsg(__FILE__, __LINE__))
    end if

    call get_subgrid_level_gindex(subgrid_level=trim(subgrid_level), gindex=gindex)
    do i=bounds1,bounds2
       get_global_index_array(i) = gindex(subgrid_index(i) - beg_index + 1)
    enddo

  end function get_global_index_array

  !-----------------------------------------------------------------------
  subroutine write_point_context(subgrid_index, subgrid_level)

    !-----------------------------------------------------------------------
    ! Description:
    ! Write various information giving context for the given index at the given subgrid
    ! level, including global index information and more.
    !
    use shr_sys_mod  , only : shr_sys_flush
    use shr_sys_mod  , only : shr_sys_abort
    use shr_log_mod  , only : errMsg => shr_log_errMsg
    use clm_varctl   , only : iulog
    use clm_varcon   , only : nameg, namel, namec, namep
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

end module GetGlobalValuesMod
