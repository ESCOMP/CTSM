#include <misc.h>
#include <preproc.h>

module RunoffMod

#if (defined RTM)
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: RunoffMod
!
! !DESCRIPTION:
! Module containing utilities for history file and coupler runoff data
!
! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use abortutils  , only : endrun
  use clm_mct_mod
  use clmtype     , only : lndrof, ocnrof, allrof

! !PUBLIC TYPES:
  implicit none
  private

  public :: runoff_flow
  type runoff_flow
!    - local
     real(r8), pointer :: runoff(:)      ! RTM flow (m**3 H2O/s)
     real(r8), pointer :: dvolrdt(:)     ! RTM change in storage (m**3/s)
     real(r8), pointer :: lonc(:)        ! lon of cell
     real(r8), pointer :: latc(:)        ! lat of cell
     real(r8), pointer :: area(:)        ! area of cell
     integer , pointer :: mask(:)        ! mask of cell 0=none, 1=lnd, 2=ocn
     integer , pointer :: dsi(:)         ! downstream index
!    - global
     real(r8), pointer :: rlon(:)        ! rtm longitude list, 1d
     real(r8), pointer :: rlat(:)        ! rtm latitude list, 1d
     integer , pointer :: glo2gdc(:)     ! rtm global to gdc index
     integer , pointer :: gdc2glo(:)     ! rtm gdc to global index
     integer , pointer :: gdc2gsn(:)     ! rtm gdc to gsn index
     integer , pointer :: gdc2i(:)       ! rtm gdc to i index
     integer , pointer :: gdc2j(:)       ! rtm gdc to j index
     integer , pointer :: num_rtm(:)     ! num of cells on each pe
!    - local
     integer           :: begr,endr      ! local start/stop indices
     integer           :: lnumr          ! rtm gdc local number of cells
     integer           :: begrl,endrl    ! local start/stop indices
     integer           :: lnumrl         ! rtm gdc local number of lnd cells
     integer           :: begro,endro    ! local start/stop indices
     integer           :: lnumro         ! rtm gdc local number of ocn cells
     integer           :: numr           ! rtm gdc global number of cells
     integer           :: numrl          ! rtm gdc global number of lnd cells
     integer           :: numro          ! rtm gdc global number of ocn cells
  end type runoff_flow
!
  type (runoff_flow) ,public :: runoff
  type(mct_gsMap)    ,public :: gsMap_rtm_gdc2glo
  integer,allocatable,public ::  perm_rtm_gdc2glo(:)
  type(mct_sMatP)    ,public :: sMatP_l2r
!
! !PUBLIC MEMBER FUNCTIONS:
  public get_proc_rof_bounds
  public get_proc_rof_total
  public get_proc_rof_global

  interface map_rof_dc2sn
     module procedure map_rof_dc2sn_sl_int
     module procedure map_rof_dc2sn_sl_real
  end interface
  public map_rof_dc2sn

  interface map_rof_sn2dc
     module procedure map_rof_sn2dc_sl_int
     module procedure map_rof_sn2dc_sl_real
  end interface
  public map_rof_sn2dc
!
! !REVISION HISTORY:
! Mariana Vertenstein: Created 10/2003
!
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_proc_rof_total
!
! !INTERFACE:
  subroutine get_proc_rof_total(pid, num_rtm)
!
! !DESCRIPTION:
! Determine number of land and ocean runoff points for this process
!
! !USES:
    use shr_kind_mod, only : r8 => shr_kind_r8
!
! !ARGUMENTS:
    implicit none
    integer, intent(in)  :: pid
!    integer, intent(out) :: num_lnd
!    integer, intent(out) :: num_ocn
    integer, intent(out) :: num_rtm
!
! !REVISION HISTORY:
! Mariana Vertenstein: Created 10/2003
!
!EOP
!-----------------------------------------------------------------------

     num_rtm = runoff%num_rtm(pid)

  end subroutine get_proc_rof_total

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_proc_rof_bounds
!
! !INTERFACE:
  subroutine get_proc_rof_bounds(beg_rtm, end_rtm)
!
! !DESCRIPTION:
! Determine beginning and ending indices of land and ocean runoff
! for this processor.
!
! !USES:
    use shr_kind_mod, only : r8 => shr_kind_r8
!
! !ARGUMENTS:
    implicit none
    integer, intent(out) :: beg_rtm
    integer, intent(out) :: end_rtm
!
! !REVISION HISTORY:
! Mariana Vertenstein: Created 10/2003
!
!EOP
!-----------------------------------------------------------------------

    beg_rtm = runoff%begr
    end_rtm = runoff%endr

  end subroutine get_proc_rof_bounds
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_proc_rof_global
!
! !INTERFACE:
  subroutine get_proc_rof_global(num_rtm, num_lnd, num_ocn)
!
! !DESCRIPTION:
! Determine number of land and ocean runoff points across all processors.
!
! !USES:
    use shr_kind_mod, only : r8 => shr_kind_r8
    use spmdMod     , only : npes
!
! !ARGUMENTS:
    implicit none
    integer, intent(out) :: num_rtm
    integer, intent(out) :: num_lnd
    integer, intent(out) :: num_ocn
!
! !REVISION HISTORY:
! Mariana Vertenstein: Created 10/2003
!
!EOP
!
! LOCAL VARIABLES:
    integer :: pid
!-----------------------------------------------------------------------

    num_rtm = runoff%numr
    num_lnd = runoff%numrl
    num_ocn = runoff%numro

  end subroutine get_proc_rof_global

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: map_rof_dc2sn_sl_int
!
! !INTERFACE:
  subroutine map_rof_dc2sn_sl_int(arraydc,arraysn,type1d)
!
! !DESCRIPTION:
! Determine beginning and ending indices of land and ocean runoff
! for this processor.
!
! !USES:
    use shr_kind_mod, only : r8 => shr_kind_r8
!
! !ARGUMENTS:
    implicit none
    integer, pointer :: arraydc(:)
    integer, pointer :: arraysn(:)
    character(len=*), intent(in) :: type1d
!
! !REVISION HISTORY:
! Mariana Vertenstein: Created 10/2003
!
!EOP
! !LOCAL VARIABLES:
    integer n
!-----------------------------------------------------------------------

    if (type1d == allrof .or. type1d == lndrof .or. type1d == ocnrof) then
       arraysn = 0
       do n = 1,runoff%numr
          arraysn(runoff%gdc2gsn(n)) = arraydc(n)
       enddo
    else
       write(6,*) 'map_rof_dc2sn type1d invalid ',type1d
       call endrun()
    endif

  end subroutine map_rof_dc2sn_sl_int
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: map_rof_dc2sn_sl_real
!
! !INTERFACE:
  subroutine map_rof_dc2sn_sl_real(arraydc,arraysn,type1d)
!
! !DESCRIPTION:
! Determine beginning and ending indices of land and ocean runoff
! for this processor.
!
! !USES:
    use shr_kind_mod, only : r8 => shr_kind_r8
!
! !ARGUMENTS:
    implicit none
    real(r8), pointer :: arraydc(:)
    real(r8), pointer :: arraysn(:)
    character(len=*), intent(in) :: type1d
!
! !REVISION HISTORY:
! Mariana Vertenstein: Created 10/2003
!
!EOP
! !LOCAL VARIABLES:
    integer n
!-----------------------------------------------------------------------

    if (type1d == allrof .or. type1d == lndrof .or. type1d == ocnrof) then
       arraysn = 0._r8
       do n = 1,runoff%numr
          arraysn(runoff%gdc2gsn(n)) = arraydc(n)
       enddo
    else
       write(6,*) 'map_rof_dc2sn type1d invalid ',type1d
       call endrun()
    endif

  end subroutine map_rof_dc2sn_sl_real
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: map_rof_sn2dc_sl_int
!
! !INTERFACE:
  subroutine map_rof_sn2dc_sl_int(arraysn,arraydc,type1d)
!
! !DESCRIPTION:
! Determine beginning and ending indices of land and ocean runoff
! for this processor.
!
! !USES:
    use shr_kind_mod, only : r8 => shr_kind_r8
!
! !ARGUMENTS:
    implicit none
    integer, pointer :: arraysn(:)
    integer, pointer :: arraydc(:)
    character(len=*), intent(in) :: type1d
!
! !REVISION HISTORY:
! Mariana Vertenstein: Created 10/2003
!
!EOP
! !LOCAL VARIABLES:
    integer n
!-----------------------------------------------------------------------

    if (type1d == allrof .or. type1d == lndrof .or. type1d == ocnrof) then
       arraydc = 0
       do n = 1,runoff%numr
          arraydc(n) = arraysn(runoff%gdc2gsn(n))
       enddo
    else
       write(6,*) 'map_rof_sn2dc type1d invalid ',type1d
       call endrun()
    endif

  end subroutine map_rof_sn2dc_sl_int
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: map_rof_sn2dc_sl_real
!
! !INTERFACE:
  subroutine map_rof_sn2dc_sl_real(arraysn,arraydc,type1d)
!
! !DESCRIPTION:
! Determine beginning and ending indices of land and ocean runoff
! for this processor.
!
! !USES:
    use shr_kind_mod, only : r8 => shr_kind_r8
!
! !ARGUMENTS:
    implicit none
    real(r8), pointer :: arraysn(:)
    real(r8), pointer :: arraydc(:)
    character(len=*), intent(in) :: type1d
!
! !REVISION HISTORY:
! Mariana Vertenstein: Created 10/2003
!
!EOP
! !LOCAL VARIABLES:
    integer n
!-----------------------------------------------------------------------

    if (type1d == allrof .or. type1d == lndrof .or. type1d == ocnrof) then
       arraydc = 0._r8
       do n = 1,runoff%numr
          arraydc(n) = arraysn(runoff%gdc2gsn(n))
       enddo
    else
       write(6,*) 'map_rof_sn2dc type1d invalid ',type1d
       call endrun()
    endif

   end subroutine map_rof_sn2dc_sl_real
!-----------------------------------------------------------------------
#endif

end module RunoffMod
