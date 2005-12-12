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
!
! !PUBLIC TYPES:
  implicit none
  save
!
  type runoff_flow
     real(r8), pointer :: lnd(:)         ! RTM river (channel) flow (m**3 H2O /s)
     real(r8), pointer :: ocn(:)         ! RTM river flow into ocean (m**3/s)
     real(r8), pointer :: lnd_dvolrdt(:) ! RTM river flow into ocean change in storage (mm/s)
     real(r8), pointer :: ocn_dvolrdt(:) ! RTM river (channel) change in storage (mm/s)
     integer , pointer :: lnd_ixy(:)     ! RTM longitude index of channel runoff point
     integer , pointer :: ocn_ixy(:)     ! RTM longitude index of ocean runoff point
     integer , pointer :: lnd_jxy(:)     ! RTM latitude index of channel runoff point
     integer , pointer :: ocn_jxy(:)     ! RTM latitude index of ocean runoff point
     real(r8), pointer :: lnd_area(:)    ! RTM gridcell area (km^2)
     real(r8), pointer :: ocn_area(:)    ! RTM gridcell area (km^2)
     integer :: beg_ocn                  ! RTM beginning ocn runoff indices on this processor
     integer :: end_ocn                  ! RTM ending ocn runoff indices on this processor
     integer :: beg_lnd                  ! RTM beginning land runoff indices on this processor
     integer :: end_lnd                  ! RTM ending land runoff indices on this processor
     integer , pointer :: num_ocn(:)     ! RTM per- proc ending ocn runoff indices
     integer , pointer :: num_lnd(:)     ! RTM per- proc ending ocn runoff indices
  end type runoff_flow
!
  type (runoff_flow) :: runoff
  public runoff_flow
!
! !PUBLIC MEMBER FUNCTIONS:
  public set_roflnd
  public set_rofocn
  public set_proc_rof_bounds
  public get_proc_rof_total
  public get_proc_rof_global
  public UpdateRunoff
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
! !IROUTINE: set_roflnd
!
! !INTERFACE:
  subroutine set_roflnd(nrlon, nrlat, mask, area, nsize)
!
! !DESCRIPTION:
! Allocate memory and initialize ocean runoff
!
! !USES:
    use shr_kind_mod, only : r8 => shr_kind_r8
!
! !ARGUMENTS:
    implicit none
    integer , intent(in) :: nrlon
    integer , intent(in) :: nrlat
    integer , intent(in) :: mask(nrlon,nrlat)
    real(r8), intent(in) :: area(nrlon,nrlat)
    integer , intent(in) :: nsize
!
! !REVISION HISTORY:
! Mariana Vertenstein: Created 10/2003
!
!EOP
!
! LOCAL VARIABLES:
    integer :: i,j,n  ! indices
    integer :: ier    ! error status
!-----------------------------------------------------------------------

    allocate(runoff%lnd(nsize), runoff%lnd_ixy(nsize), runoff%lnd_jxy(nsize), &
         runoff%lnd_area(nsize), runoff%lnd_dvolrdt(nsize), stat=ier)
    if (ier /= 0) then
       write(6,*)'set_roflnd: allocation error for runoff components lnd, lnd_ixy, lnd_jxy, and lnd_area'
       call endrun
    end if

    n = 0
    do j = 1,nrlat
       do i = 1,nrlon
          if (mask(i,j) == 1) then
             n = n + 1
             runoff%lnd_ixy(n) = i
             runoff%lnd_jxy(n) = j
             runoff%lnd_area(n) = area(i,j)
             runoff%lnd(n) = 0._r8
             runoff%lnd_dvolrdt(n) = 0._r8
          endif
       end do
    end do

  end subroutine set_roflnd

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: set_rofocn
!
! !INTERFACE:
  subroutine set_rofocn(nrlon, nrlat, mask, area, nsize)
!
! !DESCRIPTION:
! Allocate memory and initialize ocean runoff
!
! !USES:
    use shr_kind_mod, only : r8 => shr_kind_r8
!
! !ARGUMENTS:
    implicit none
    integer , intent(in) :: nrlon
    integer , intent(in) :: nrlat
    integer , intent(in) :: mask(nrlon,nrlat)
    real(r8), intent(in) :: area(nrlon,nrlat)
    integer , intent(in) :: nsize

!
! !REVISION HISTORY:
! Mariana Vertenstein: Created 10/2003
!
!EOP
!
! LOCAL VARIABLES:
    integer :: i,j,n  ! indices
    integer :: ier    ! error status
!-----------------------------------------------------------------------

    allocate(runoff%ocn(nsize), runoff%ocn_ixy(nsize), runoff%ocn_jxy(nsize), &
         runoff%ocn_area(nsize), runoff%ocn_dvolrdt(nsize), stat=ier)
    if (ier /= 0) then
       write(6,*)'set_rofocn: allocation error for runoff components ocn, ocn_ixy, ocn_jxy, and ocn_area'
       call endrun
    end if

    n = 0
    do j = 1,nrlat
       do i = 1,nrlon
          if (mask(i,j) == 1) then
             n = n + 1
             runoff%ocn_ixy(n) = i
             runoff%ocn_jxy(n) = j
             runoff%ocn_area(n) = area(i,j)
             runoff%ocn(n) = 0._r8
             runoff%ocn_dvolrdt(n) = 0._r8
          end if
       end do
    end do

  end subroutine set_rofocn

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: set_proc_rof_bounds
!
! !INTERFACE:
  subroutine set_proc_rof_bounds()
!
! !DESCRIPTION:
! Set the beginning and ending indices of the per-process
! beginning and ending indices for land and ocean runoff
!
! !USES:
    use shr_kind_mod, only : r8 => shr_kind_r8
    use spmdMod     , only : npes, iam
!
! !ARGUMENTS:
    implicit none
!
! !REVISION HISTORY:
! Mariana Vertenstein: Created 10/2003
!
!EOP
!
! LOCAL VARIABLES:
    integer :: pid               ! processor id
    integer :: nlnd,nocn         ! number of land, ocean runoff points
    integer :: ier               ! error status
    integer, allocatable :: beg_lnd(:), end_lnd(:)  ! beginning, ending land runoff indices
    integer, allocatable :: beg_ocn(:), end_ocn(:)  ! beginning, ending ocean runoff indices
!-----------------------------------------------------------------------

    allocate (beg_lnd(0:npes-1), end_lnd(0:npes-1), beg_ocn(0:npes-1), end_ocn(0:npes-1), &
         stat=ier)
    if (ier /= 0) then
       write(6,*)'set_proc_rof_bounds: alloc err for temporaries beg_lnd, end_lnd, beg_ocn, end_ocn'
       call endrun()
    end if

    allocate (runoff%num_lnd(0:npes-1), runoff%num_ocn(0:npes-1), stat=ier)
    if (ier /= 0) then
       write(6,*)'set_proc_rof_bounds: allocation error for runoff%num_lnd, runoff%num_ocn'
       call endrun()
    end if

    nlnd = size(runoff%lnd)/npes
    nocn = size(runoff%ocn)/npes

    if (npes > nlnd) then
       write (6,*) 'set_proc_rof_bounds: number of processes exceeds number of land runoff cells'
       call endrun
    else if (npes > nocn) then
       write (6,*) 'set_proc_rof_bounds: number of processes exceeds number of ocean runoff cells'
       call endrun
    end if

    do pid = 0,npes-1

       beg_lnd(pid) = pid*nlnd + 1
       beg_ocn(pid) = pid*nocn + 1
       if (pid < npes-1) then
          end_lnd(pid) = beg_lnd(pid) + nlnd -1
          end_ocn(pid) = beg_ocn(pid) + nocn -1
       else
          end_lnd(pid) = size(runoff%lnd)
          end_ocn(pid) = size(runoff%ocn)
       end if
       if (pid == iam) then
          runoff%beg_lnd = beg_lnd(pid)
          runoff%end_lnd = end_lnd(pid)
          runoff%beg_ocn = beg_ocn(pid)
          runoff%end_ocn = end_ocn(pid)
       end if
       runoff%num_lnd(pid) = end_lnd(pid) - beg_lnd(pid) + 1
       runoff%num_ocn(pid) = end_ocn(pid) - beg_ocn(pid) + 1

       if (runoff%num_lnd(pid) <= 0) then
          write(6,*)'set_proc_rof_bounds_error: num_lnd ',runoff%num_lnd(pid), &
               ' is invalid for pid ',pid
          call endrun()
       else if (runoff%num_ocn(pid) <= 0) then
          write(6,*)'set_proc_rof_bounds_error: num_ocn ',runoff%num_ocn(pid), &
               ' is invalid for pid ',pid
          call endrun()
       end if

    end do

    write(6,*)'iam= ',iam,' beg_ocn= ',beg_ocn(iam),' end_ocn= ',end_ocn(iam), &
         ' num_ocn= ',runoff%num_ocn(iam)
    write(6,*)'iam= ',iam,' beg_lnd= ',beg_lnd(iam),' end_lnd= ',end_lnd(iam), &
         ' num_lnd= ',runoff%num_lnd(iam)

    deallocate(beg_lnd, end_lnd, beg_ocn, end_ocn)

  end subroutine set_proc_rof_bounds

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: UpdateRunoff
!
! !INTERFACE:
  subroutine UpdateRunoff(nrloni, nrlonf, nrlati, nrlatf, flxocn_r, flxlnd_r, &
                          dvolrdt_ocn_r, dvolrdt_lnd_r)
!
! !DESCRIPTION:
! Update the land and ocean runoff vectors. This determine the
! ocean runoff vectors to send to the coupler as well as the
! ocean and land runoff vectors for history output.
!
! !USES:
    use shr_kind_mod, only : r8 => shr_kind_r8
!
! !ARGUMENTS:
    implicit none
    integer , intent(in) :: nrloni
    integer , intent(in) :: nrlonf
    integer , intent(in) :: nrlati
    integer , intent(in) :: nrlatf
    real(r8), intent(in) :: flxocn_r(nrloni:nrlonf, nrlati:nrlatf)
    real(r8), intent(in) :: flxlnd_r(nrloni:nrlonf, nrlati:nrlatf)
    real(r8), intent(in) :: dvolrdt_ocn_r(nrloni:nrlonf, nrlati:nrlatf)
    real(r8), intent(in) :: dvolrdt_lnd_r(nrloni:nrlonf, nrlati:nrlatf)
!
! !REVISION HISTORY:
! Mariana Vertenstein: Created 10/2003
!
!EOP
!
! LOCAL VARIABLES:
    integer :: i,j,n   ! indices
!-----------------------------------------------------------------------

    do n = 1,size(runoff%ocn)
       i = runoff%ocn_ixy(n)
       j = runoff%ocn_jxy(n)
       runoff%ocn(n) = flxocn_r(i,j)
       runoff%ocn_dvolrdt(n) = dvolrdt_ocn_r(i,j)
    end do

    do n = 1,size(runoff%lnd)
       i = runoff%lnd_ixy(n)
       j = runoff%lnd_jxy(n)
       runoff%lnd(n) = flxlnd_r(i,j)
       runoff%lnd_dvolrdt(n) = dvolrdt_lnd_r(i,j)
    end do

  end subroutine UpdateRunoff

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_proc_rof_total
!
! !INTERFACE:
  subroutine get_proc_rof_total(pid, num_lnd, num_ocn)
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
    integer, intent(out) :: num_lnd
    integer, intent(out) :: num_ocn
!
! !REVISION HISTORY:
! Mariana Vertenstein: Created 10/2003
!
!EOP
!-----------------------------------------------------------------------

    num_lnd = runoff%num_lnd(pid)
    num_ocn = runoff%num_ocn(pid)

  end subroutine get_proc_rof_total

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_proc_rof_bounds
!
! !INTERFACE:
  subroutine get_proc_rof_bounds(beg_lnd, end_lnd, beg_ocn, end_ocn)
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
    integer, intent(out) :: beg_lnd
    integer, intent(out) :: end_lnd
    integer, intent(out) :: beg_ocn
    integer, intent(out) :: end_ocn
!
! !REVISION HISTORY:
! Mariana Vertenstein: Created 10/2003
!
!EOP
!-----------------------------------------------------------------------

    beg_lnd = runoff%beg_lnd
    end_lnd = runoff%end_lnd
    beg_ocn = runoff%beg_ocn
    end_ocn = runoff%end_ocn

  end subroutine get_proc_rof_bounds

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_proc_rof_global
!
! !INTERFACE:
  subroutine get_proc_rof_global(num_lnd, num_ocn)
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

    num_lnd = 0
    num_ocn = 0
    do pid = 0,npes-1
       num_lnd = num_lnd + runoff%num_lnd(pid)
       num_ocn = num_ocn + runoff%num_ocn(pid)
    end do

  end subroutine get_proc_rof_global

#endif

end module RunoffMod
