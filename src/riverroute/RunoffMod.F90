
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
  use clmtype     , only : allrof

! !PUBLIC TYPES:
  implicit none
  private

  integer,parameter,public :: nt_rtm = 2    ! number of tracers
  character(len=3),parameter,public :: rtm_tracers(nt_rtm) = &
     (/'LIQ','ICE'/)

  public :: runoff_flow
  type runoff_flow
!    - local
     real(r8), pointer :: runoff(:,:)      ! RTM flow (m**3 H2O/s)
     real(r8), pointer :: runofflnd(:,:)   ! runoff masked for land (m**3 H2O/s)
     real(r8), pointer :: runoffocn(:,:)   ! runoff masked for ocn  (m**3 H2O/s)
     real(r8), pointer :: dvolrdt(:,:)     ! RTM change in storage (mm/s)
     real(r8), pointer :: dvolrdtlnd(:,:)  ! dvolrdt masked for land (mm/s)
     real(r8), pointer :: dvolrdtocn(:,:)  ! dvolrdt masked for ocn  (mm/s)
     real(r8), pointer :: volr(:,:)        ! RTM storage (m**3)
     real(r8), pointer :: volrlnd(:,:)     ! RTM storage masked for land (m**3)
     real(r8), pointer :: lonc(:)          ! lon of cell
     real(r8), pointer :: latc(:)          ! lat of cell
     real(r8), pointer :: area(:)          ! area of cell
     integer , pointer :: gindex(:)        ! global index
     integer , pointer :: mask(:)          ! mask of cell 0=none, 1=lnd, 2=ocn
     integer , pointer :: dsi(:)           ! downstream index
!    - global
     real(r8), pointer :: rlon(:)          ! rtm longitude list, 1d
     real(r8), pointer :: rlat(:)          ! rtm latitude list, 1d
     integer , pointer :: num_rtm(:)       ! num of cells on each pe
!    - local
     integer           :: begr,endr        ! local start/stop indices
     integer           :: lnumr            ! rtm gdc local number of cells
     integer           :: begrl,endrl      ! local start/stop indices
     integer           :: lnumrl           ! rtm gdc local number of lnd cells
     integer           :: begro,endro      ! local start/stop indices
     integer           :: lnumro           ! rtm gdc local number of ocn cells
     integer           :: numr             ! rtm gdc global number of cells
     integer           :: numrl            ! rtm gdc global number of lnd cells
     integer           :: numro            ! rtm gdc global number of ocn cells
!    - need 1d field pointers for history files
     real(r8), pointer :: runofflnd_nt1(:)
     real(r8), pointer :: runofflnd_nt2(:)
     real(r8), pointer :: runoffocn_nt1(:)
     real(r8), pointer :: runoffocn_nt2(:)
     real(r8), pointer :: dvolrdtlnd_nt1(:)
     real(r8), pointer :: dvolrdtlnd_nt2(:)
     real(r8), pointer :: dvolrdtocn_nt1(:)
     real(r8), pointer :: dvolrdtocn_nt2(:)
     real(r8), pointer :: volr_nt1(:)
     real(r8), pointer :: volr_nt2(:)
  end type runoff_flow
!
  type (runoff_flow)         ,public :: runoff
  type(mct_gsMap)    ,target ,public :: gsMap_rtm_gdc2glo
!
! !PUBLIC MEMBER FUNCTIONS:
  public get_proc_rof_bounds
  public get_proc_rof_total
  public get_proc_rof_global

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
    integer, optional, intent(out) :: num_lnd
    integer, optional, intent(out) :: num_ocn
!
! !REVISION HISTORY:
! Mariana Vertenstein: Created 10/2003
!
!
! !LOCAL VARIABLES:
!EOP
    integer :: pid
!-----------------------------------------------------------------------

    num_rtm = runoff%numr
    if (present(num_lnd)) then
       num_lnd = runoff%numrl
    endif
    if (present(num_ocn)) then
       num_ocn = runoff%numro
    endif

  end subroutine get_proc_rof_global

!-----------------------------------------------------------------------
#endif

end module RunoffMod
