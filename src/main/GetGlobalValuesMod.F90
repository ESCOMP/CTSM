module GetGlobalValuesMod

  !-----------------------------------------------------------------------
  ! Obtain and Write Global Index information
  !-----------------------------------------------------------------------
  implicit none
  private

  ! PUBLIC MEMBER FUNCTIONS:

  public :: GetGlobalIndex
  public :: GetGlobalIndexArray
  public :: GetGlobalWrite

  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  integer function GetGlobalIndex(decomp_index, clmlevel)

    !----------------------------------------------------------------
    ! Description
    ! Determine global index space value for target point at given clmlevel
    !
    ! Uses:
    use shr_log_mod, only: errMsg => shr_log_errMsg
    use decompMod  , only: bounds_type, get_clmlevel_gindex, get_proc_bounds
    use spmdMod    , only: iam
    use clm_varcon , only: nameg, namel, namec, namep
    use clm_varctl , only: iulog
    use shr_sys_mod, only: shr_sys_abort
    !
    ! Arguments 
    integer          , intent(in) :: decomp_index
    character(len=*) , intent(in) :: clmlevel
    !
    ! Local Variables:
    type(bounds_type) :: bounds_proc   ! processor bounds
    integer           :: beg_index     ! beginning proc index for clmlevel
    integer, pointer  :: gindex(:)
    !----------------------------------------------------------------

    call get_proc_bounds(bounds_proc)

    if (trim(clmlevel) == nameg) then
       beg_index = bounds_proc%begg
    else if (trim(clmlevel) == namel) then
       beg_index = bounds_proc%begl
    else if (trim(clmlevel) == namec) then
       beg_index = bounds_proc%begc
    else if (trim(clmlevel) == namep) then
       beg_index = bounds_proc%begp
    else
       call shr_sys_abort('clmlevel of '//trim(clmlevel)//' not supported' // &
            errmsg(sourcefile, __LINE__))
    end if

    call get_clmlevel_gindex(clmlevel=trim(clmlevel), gindex=gindex)
    GetGlobalIndex = gindex(decomp_index - beg_index + 1)

  end function GetGlobalIndex

  !-----------------------------------------------------------------------
  function GetGlobalIndexArray(decomp_index, bounds1, bounds2, clmlevel)

    !----------------------------------------------------------------
    ! Description
    ! Determine global index space value for target array at given clmlevel
    !
    ! Example from histFileMod.F90:
    ! ilarr = GetGlobalIndexArray(lun%gridcell(bounds%begl:bounds%endl), bounds%begl, bounds%endl, clmlevel=nameg)
    ! Note that the last argument (clmlevel) is set to nameg, which corresponds
    ! to the "gridcell" not the "lun" of the first argument.
    !
    ! Uses:
#include "shr_assert.h"
    use shr_log_mod, only: errMsg => shr_log_errMsg
    use decompMod  , only: bounds_type, get_clmlevel_gindex, get_proc_bounds
    use spmdMod    , only: iam
    use clm_varcon , only: nameg, namel, namec, namep
    use clm_varctl , only: iulog
    use shr_sys_mod, only: shr_sys_abort
    !
    ! Arguments 
    integer          , intent(in) :: bounds1  ! lower bound of the input & returned arrays
    integer          , intent(in) :: bounds2  ! upper bound of the input & returned arrays
    integer          , intent(in) :: decomp_index(bounds1:)
    character(len=*) , intent(in) :: clmlevel
    integer                       :: GetGlobalIndexArray(bounds1:bounds2)
    !
    ! Local Variables:
    type(bounds_type) :: bounds_proc   ! processor bounds
    integer           :: beg_index     ! beginning proc index for clmlevel
    integer           :: i
    integer , pointer :: gindex(:)
    !----------------------------------------------------------------

    SHR_ASSERT_ALL_FL((ubound(decomp_index) == (/bounds2/)), sourcefile, __LINE__)
    call get_proc_bounds(bounds_proc)

    if (trim(clmlevel) == nameg) then
       beg_index = bounds_proc%begg
    else if (trim(clmlevel) == namel) then
       beg_index = bounds_proc%begl
    else if (trim(clmlevel) == namec) then
       beg_index = bounds_proc%begc
    else if (trim(clmlevel) == namep) then
       beg_index = bounds_proc%begp
    else
       call shr_sys_abort('clmlevel of '//trim(clmlevel)//' not supported' // &
            errmsg(__FILE__, __LINE__))
    end if

    call get_clmlevel_gindex(clmlevel=trim(clmlevel), gindex=gindex)
    do i=bounds1,bounds2
       GetGlobalIndexArray(i) = gindex(decomp_index(i) - beg_index + 1)
    enddo

  end function GetGlobalIndexArray

  !-----------------------------------------------------------------------
  subroutine GetGlobalWrite(decomp_index, clmlevel)

    !-----------------------------------------------------------------------
    ! Description:
    ! Write global index information for input local indices
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
    integer          , intent(in) :: decomp_index
    character(len=*) , intent(in) :: clmlevel
    !
    ! Local Variables:
    integer :: igrc, ilun, icol, ipft 
    !-----------------------------------------------------------------------

    if (trim(clmlevel) == nameg) then

       igrc = decomp_index
       write(iulog,'(a, i0, a, i0)') 'iam = ', iam, ': local  gridcell index = ', igrc
       write(iulog,'(a, i0, a, i0)') 'iam = ', iam, ': global gridcell index = ', &
            GetGlobalIndex(decomp_index=igrc, clmlevel=nameg)
       write(iulog,'(a, i0, a, f12.7)') 'iam = ', iam, ': gridcell longitude    = ', grc%londeg(igrc)
       write(iulog,'(a, i0, a, f12.7)') 'iam = ', iam, ': gridcell latitude     = ', grc%latdeg(igrc)

    else if (trim(clmlevel) == namel) then

       ilun = decomp_index
       igrc = lun%gridcell(ilun)
       write(iulog,'(a, i0, a, i0)') 'iam = ', iam, ': local  landunit index = ', ilun
       write(iulog,'(a, i0, a, i0)') 'iam = ', iam, ': global landunit index = ', &
            GetGlobalIndex(decomp_index=ilun, clmlevel=namel)
       write(iulog,'(a, i0, a, i0)') 'iam = ', iam, ': global gridcell index = ', &
            GetGlobalIndex(decomp_index=igrc, clmlevel=nameg)
       write(iulog,'(a, i0, a, f12.7)') 'iam = ', iam, ': gridcell longitude    = ', grc%londeg(igrc)
       write(iulog,'(a, i0, a, f12.7)') 'iam = ', iam, ': gridcell latitude     = ', grc%latdeg(igrc)
       write(iulog,'(a, i0, a, i0)') 'iam = ', iam, ': landunit type         = ', lun%itype(decomp_index)

    else if (trim(clmlevel) == namec) then

       icol = decomp_index
       ilun = col%landunit(icol)
       igrc = col%gridcell(icol)
       write(iulog,'(a, i0, a, i0)') 'iam = ', iam, ': local  column   index = ', icol
       write(iulog,'(a, i0, a, i0)') 'iam = ', iam, ': global column   index = ', &
            GetGlobalIndex(decomp_index=icol, clmlevel=namec)
       write(iulog,'(a, i0, a, i0)') 'iam = ', iam, ': global landunit index = ', &
            GetGlobalIndex(decomp_index=ilun, clmlevel=namel)
       write(iulog,'(a, i0, a, i0)') 'iam = ', iam, ': global gridcell index = ', &
            GetGlobalIndex(decomp_index=igrc, clmlevel=nameg)
       write(iulog,'(a, i0, a, f12.7)') 'iam = ', iam, ': gridcell longitude    = ', grc%londeg(igrc)
       write(iulog,'(a, i0, a, f12.7)') 'iam = ', iam, ': gridcell latitude     = ', grc%latdeg(igrc)
       write(iulog,'(a, i0, a, i0)') 'iam = ', iam, ': column   type         = ', col%itype(icol)
       write(iulog,'(a, i0, a, i0)') 'iam = ', iam, ': landunit type         = ', lun%itype(ilun)

    else if (trim(clmlevel) == namep) then

       ipft = decomp_index
       icol = patch%column(ipft)
       ilun = patch%landunit(ipft)
       igrc = patch%gridcell(ipft)
       write(iulog,'(a, i0, a, i0)') 'iam = ', iam, ': local  patch    index = ', ipft
       write(iulog,'(a, i0, a, i0)') 'iam = ', iam, ': global patch    index = ', &
            GetGlobalIndex(decomp_index=ipft, clmlevel=namep)
       write(iulog,'(a, i0, a, i0)') 'iam = ', iam, ': global column   index = ', &
            GetGlobalIndex(decomp_index=icol, clmlevel=namec)
       write(iulog,'(a, i0, a, i0)') 'iam = ', iam, ': global landunit index = ', &
            GetGlobalIndex(decomp_index=ilun, clmlevel=namel)
       write(iulog,'(a, i0, a, i0)') 'iam = ', iam, ': global gridcell index = ', &
            GetGlobalIndex(decomp_index=igrc, clmlevel=nameg)
       write(iulog,'(a, i0, a, f12.7)') 'iam = ', iam, ': gridcell longitude    = ', grc%londeg(igrc)
       write(iulog,'(a, i0, a, f12.7)') 'iam = ', iam, ': gridcell latitude     = ', grc%latdeg(igrc)
       write(iulog,'(a, i0, a, i0)') 'iam = ', iam, ': pft      type         = ', patch%itype(ipft)
       write(iulog,'(a, i0, a, i0)') 'iam = ', iam, ': column   type         = ', col%itype(icol)
       write(iulog,'(a, i0, a, i0)') 'iam = ', iam, ': landunit type         = ', lun%itype(ilun)

    else		       
       call shr_sys_abort('clmlevel '//trim(clmlevel)//'not supported '//errmsg(sourcefile, __LINE__))

    end if

    call shr_sys_flush(iulog)

  end subroutine GetGlobalWrite

end module GetGlobalValuesMod
