#include <misc.h>
#include <preproc.h>

module iobinary

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: iobinary
!
! !DESCRIPTION:
! Set of wrappers to write binary I/O
!
! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
#if (defined SPMD)
  use spmdMod        , only : masterproc, mpicom, MPI_REAL8, MPI_INTEGER, &
                              MPI_LOGICAL
  use spmdGathScatMod, only : scatter_data_from_master, gather_data_to_master
#else
  use spmdMod        , only : masterproc
#endif
  use decompMod      , only : map_sn2dc, map_dc2sn
  use abortutils     , only : endrun
  use clmtype        , only : nameg, namel, namec, namep, lndrof, ocnrof
!
! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: readin, wrtout  ! read and write bindary I/O
  interface readin
     module procedure readin_1darray_int
     module procedure readin_2darray_int
     module procedure readin_1darray_real
     module procedure readin_2darray_real
  end interface
  interface wrtout
     module procedure wrtout_1darray_int
     module procedure wrtout_2darray_int
     module procedure wrtout_1darray_real
     module procedure wrtout_2darray_real
  end interface
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!
! !PRIVATE MEMBER FUNCTIONS:
  private :: getnum  !get 1d type
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: readin_1d_array_int
!
! !INTERFACE:
  subroutine readin_1darray_int (iu, iarr, clmlevel)
!
! !DESCRIPTION:
! Wrapper routine to read integer 1d array from restart binary file
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: iu                  !input unit
    integer, pointer    :: iarr(:)             !input data
    character(len=*), intent(in) :: clmlevel   !type of input data
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: ier                             !return code
    integer :: nsize                           !size of first dimension
    integer, pointer :: iglobdc(:)             !temporary
    integer, pointer :: iglobsn(:)             !temporary
!-----------------------------------------------------------------------
    if (masterproc) then
       nsize = getnum (clmlevel)
       allocate (iglobdc(nsize), iglobsn(nsize), stat=ier)
       if (ier /=0) then
          write(6,*)'readin_1d_array_int allocation error'; call endrun()
       end if
       read (iu,iostat=ier) iglobsn
       if (ier /= 0 ) then
          write (6,*)'readin_1darray_int error ',ier,' on i/o unit = ',iu
          call endrun()
       endif
       if (clmlevel == nameg .or. clmlevel == namel .or. &
           clmlevel == namec .or. clmlevel == namep) then
          call map_sn2dc(iglobsn, iglobdc, clmlevel)
       else
          iglobdc(:) = iglobsn(:)
       end if
    end if
#if (defined SPMD)
    call scatter_data_from_master(iarr, iglobdc, clmlevel=clmlevel)
#else
    iarr(:) = iglobdc(:)
#endif
    if (masterproc) deallocate(iglobdc, iglobsn)

  end subroutine readin_1darray_int

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: readin_1darray_real
!
! !INTERFACE:
  subroutine readin_1darray_real (iu, rarr, clmlevel)
!
! !DESCRIPTION:
! Wrapper routine to read real 1d array from restart binary file

! !ARGUMENTS:
    implicit none
    integer, intent(in) :: iu                  !input unit
    real(r8), pointer   :: rarr(:)             !input data
    character(len=*), intent(in) :: clmlevel   !input data type
!
! !REVISION HISTORY:
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: ier                              !return code
    integer :: nsize                            !size of first dimension
    real(r8), pointer :: rglobdc(:)             !temporary
    real(r8), pointer :: rglobsn(:)             !temporary
!-----------------------------------------------------------------------
    if (masterproc) then
       nsize = getnum (clmlevel)
       allocate (rglobdc(nsize), rglobsn(nsize), stat=ier )
       if (ier /=0) then
          write(6,*)'readin_1d_array_real allocation error'; call endrun()
       end if
       read (iu,iostat=ier) rglobsn
       if (ier /= 0 ) then
          write(6,*)'readin_1darray_real error ',ier,' on i/o unit = ',iu
          call endrun()
       endif
       if (clmlevel == nameg .or. clmlevel == namel .or. &
           clmlevel == namec .or. clmlevel == namep) then
          call map_sn2dc(rglobsn, rglobdc, clmlevel)
       else
          rglobdc(:) = rglobsn(:)
       end if
    endif
#if (defined SPMD)
    call scatter_data_from_master(rarr, rglobdc, clmlevel=clmlevel)
#else
    rarr(:) = rglobdc(:)
#endif
    if (masterproc)deallocate(rglobdc, rglobsn)

  end subroutine readin_1darray_real

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: readin_2d_arrayint
!
! !INTERFACE:
  subroutine readin_2darray_int (iu, iarr, clmlevel)
!
! !DESCRIPTION:
! Wrapper routine to read integer 2d array from restart binary file
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: iu                  !input unit
    integer, pointer    :: iarr(:,:)           !input data
    character(len=*), intent(in) :: clmlevel   !type of input data
!
! !REVISION HISTORY:
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: ier                               !error status
    integer :: lb,ub                             !bounds of first dimension
    integer :: nsize                             !size of second dimension
    integer, pointer :: iglobdc(:,:)             !temporary
    integer, pointer :: iglobsn(:,:)             !temporary
!-----------------------------------------------------------------------
    if (masterproc) then
       nsize = getnum (clmlevel)
       lb = lbound(iarr, dim=1)
       ub = ubound(iarr, dim=1)
       allocate (iglobdc(lb:ub,nsize), iglobsn(lb:ub,nsize), stat=ier)
       if (ier /=0) then
          write(6,*)'readin_2d_array_int allocation error'; call endrun()
       end if
       read (iu,iostat=ier) iglobsn
       if (ier /= 0 ) then
          write(6,*)'readin_2darray_int error ',ier,' on i/o unit = ',iu
          call endrun()
       endif
       if (clmlevel == nameg .or. clmlevel == namel .or. &
           clmlevel == namec .or. clmlevel == namep) then
          call map_sn2dc(iglobsn, iglobdc, clmlevel, lb, ub)
       else
          iglobdc(:,:) = iglobsn(:,:)
       end if
    endif
#if (defined SPMD)
    call scatter_data_from_master(iarr, iglobdc, clmlevel=clmlevel)
#else
    iarr(:,:) = iglobdc(:,:)
#endif
    if (masterproc) deallocate(iglobdc, iglobsn)

  end subroutine readin_2darray_int

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: readin_2darray_real
!
! !INTERFACE:
  subroutine readin_2darray_real (iu, rarr, clmlevel)
!
! !DESCRIPTION:
! Wrapper routine to read real 2d array from restart binary file
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: iu                 !input unit
    real(r8), pointer   :: rarr(:,:)          !input data
    character(len=*), intent(in) :: clmlevel  !type of input data
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: ier                               !return code
    integer :: lb,ub                             !bounds of first dimension
    integer :: nsize                             !size of second dimension
    real(r8), pointer :: rglobdc(:,:)             !temporary
    real(r8), pointer :: rglobsn(:,:)             !temporary
!-----------------------------------------------------------------------
    if (masterproc) then
       nsize = getnum (clmlevel)
       lb = lbound(rarr, dim=1)
       ub = ubound(rarr, dim=1)
       allocate (rglobdc(lb:ub,nsize), rglobsn(lb:ub,nsize), stat=ier)
       if (ier /=0) then
          write(6,*)'readin_2d_array_real allocation error'; call endrun()
       end if
       read (iu,iostat=ier) rglobsn
       if (ier /= 0 ) then
          write(6,*)'readin_2darray_real error ',ier,' on i/o unit = ',iu
          call endrun()
       endif
       if (clmlevel == nameg .or. clmlevel == namel .or. &
           clmlevel == namec .or. clmlevel == namep) then
          call map_sn2dc(rglobsn, rglobdc, clmlevel, lb, ub)
       else
          rglobdc(:,:) = rglobsn(:,:)
       end if
    endif
#if (defined SPMD)
    call scatter_data_from_master(rarr, rglobdc, clmlevel=clmlevel)
#else
    rarr(:,:) = rglobdc(:,:)
#endif
    if (masterproc) deallocate(rglobdc, rglobsn)

  end subroutine readin_2darray_real

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: wrtout_1d_array_int
!
! !INTERFACE:
  subroutine wrtout_1darray_int (iu, iarr, clmlevel)
!
! !DESCRIPTION:
! Wrapper routine to write integer array to restart binary file
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: iu                 !output unit
    integer, pointer    :: iarr(:)            !output data
    character(len=*), intent(in) :: clmlevel  !output 1d type
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: ier                             !errorcode
    integer :: nsize                           !size of first dimension
    integer, pointer :: iglobdc(:)             !temporary
    integer, pointer :: iglobsn(:)             !temporary
!-----------------------------------------------------------------------
    if (masterproc) then
       nsize = getnum (clmlevel)
       allocate (iglobdc(nsize), iglobsn(nsize), stat=ier)
       if (ier /=0) then
          write(6,*)'wrtout_1d_array_int allocation error'; call endrun()
       end if
    endif
#if (defined SPMD)
    call gather_data_to_master(iarr, iglobdc, clmlevel=clmlevel)
#else
    iglobdc(:) = iarr(:)
#endif
    if (masterproc) then
       if (clmlevel == nameg .or. clmlevel == namel .or. &
           clmlevel == namec .or. clmlevel == namep) then
          call map_dc2sn(iglobdc, iglobsn, clmlevel)
       else
          iglobsn(:) = iglobdc(:)
       end if
       write (iu, iostat=ier) iglobsn
       if (ier /= 0 ) then
          write(6,*)'wrtout_1darray_int error ',ier,' on i/o unit = ',iu
          call endrun
       end if
       deallocate(iglobdc, iglobsn)
    end if

  end subroutine wrtout_1darray_int

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: wrtout_1d_array_real
!
! !INTERFACE:
  subroutine wrtout_1darray_real (iu, rarr, clmlevel)
!
! !DESCRIPTION:
! Wrapper routine to write real array to restart binary file
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: iu                 !output unit
    real(r8), pointer   :: rarr(:)            !output data
    character(len=*), intent(in) :: clmlevel  !input data type
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: ier                            !return code
    integer :: nsize                          !size of first dimension
    real(r8), pointer :: rglobdc(:)           !temporary
    real(r8), pointer :: rglobsn(:)           !temporary
!-----------------------------------------------------------------------
    if (masterproc) then
       nsize = getnum (clmlevel)
       allocate (rglobdc(nsize), rglobsn(nsize), stat=ier)
       if (ier /=0) then
          write(6,*)'wrtout_1d_array_real allocation error'; call endrun()
       end if
    endif
#if (defined SPMD)
    call gather_data_to_master(rarr, rglobdc, clmlevel=clmlevel)
#else
    rglobdc(:) = rarr(:)
#endif
    if (masterproc) then
       if (clmlevel == nameg .or. clmlevel == namel .or. &
           clmlevel == namec .or. clmlevel == namep) then
          call map_dc2sn(rglobdc, rglobsn, clmlevel)
       else
          rglobsn(:) = rglobdc(:)
       end if
       write (iu, iostat=ier) rglobsn
       if (ier /= 0 ) then
          write(6,*)'wrtout_1darray_real error ',ier,' on i/o unit = ',iu
          call endrun
       end if
       deallocate(rglobdc, rglobsn)
    end if

  end subroutine wrtout_1darray_real

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: wrtout_2d_array_int
!
! !INTERFACE:
  subroutine wrtout_2darray_int (iu, iarr, clmlevel)
!
! !DESCRIPTION:
! Wrapper routine to write integer array to restart binary file
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: iu                !output unit
    integer, pointer    :: iarr(:,:)         !output data
    character(len=*), intent(in) :: clmlevel !output data 1d type
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: ier                           !return code
    integer :: lb,ub                         !bounds of first dimension
    integer :: nsize                         !size of second dimension
    integer, pointer :: iglobdc(:,:)         !temporary
    integer, pointer :: iglobsn(:,:)         !temporary
!-----------------------------------------------------------------------
    if (masterproc) then
       nsize = getnum (clmlevel)
       lb = lbound(iarr, dim=1)
       ub = ubound(iarr, dim=1)
       allocate (iglobdc(lb:ub,nsize), iglobsn(lb:ub,nsize), stat=ier)
       if (ier /=0) then
          write(6,*)'wrtout_2d_array_int allocation error'; call endrun()
       end if
    endif
#if (defined SPMD)
    call gather_data_to_master(iarr, iglobdc, clmlevel=clmlevel)
#else
    iglobdc(:,:) = iarr(:,:)
#endif
    if (masterproc) then
       if (clmlevel == nameg .or. clmlevel == namel .or. &
           clmlevel == namec .or. clmlevel == namep) then
          call map_dc2sn(iglobdc, iglobsn, clmlevel, lb, ub)
       else
          iglobsn(:,:) = iglobdc(:,:)
       end if
       write (iu, iostat=ier) iglobsn
       if (ier /= 0 ) then
          write(6,*)'wrtout_2darray_int error ',ier,' on i/o unit = ',iu
          call endrun()
       end if
       deallocate(iglobdc, iglobsn)
    endif

  end subroutine wrtout_2darray_int

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: wrtout_2darray_real
!
! !INTERFACE:
  subroutine wrtout_2darray_real (iu, rarr, clmlevel)
!
! !DESCRIPTION:
! Wrapper routine to write real array to restart binary file
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: iu                 !input unit
    real(r8), pointer   :: rarr(:,:)          !output data
    character(len=*), intent(in) :: clmlevel  !type of input data
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: ier                            !return code
    integer :: lb,ub                          !bounds of first dimension
    integer :: nsize                          !size of second dimension
    real(r8), pointer :: rglobdc(:,:)         !temporary
    real(r8), pointer :: rglobsn(:,:)         !temporary
!-----------------------------------------------------------------------
    if (masterproc) then
       nsize = getnum (clmlevel)
       lb = lbound(rarr, dim=1)
       ub = ubound(rarr, dim=1)
       allocate (rglobdc(lb:ub,nsize), rglobsn(lb:ub,nsize), stat=ier)
       if (ier /=0) then
          write(6,*)'readin_2d_array_real allocation error'; call endrun()
       end if
    endif
#if (defined SPMD)
    call gather_data_to_master(rarr, rglobdc, clmlevel=clmlevel)
#else
    rglobdc(:,:) = rarr(:,:)
#endif
    if (masterproc) then
       if (clmlevel == nameg .or. clmlevel == namel .or. &
           clmlevel == namec .or. clmlevel == namep) then
          call map_dc2sn(rglobdc, rglobsn, clmlevel, lb, ub)
       else
          rglobsn(:,:) = rglobdc(:,:)
       end if
       write (iu, iostat=ier) rglobsn
       if (ier /= 0 ) then
          write(6,*)'wrtout_2darray_real error ',ier,' on i/o unit = ',iu
          call endrun
       end if
       deallocate(rglobdc, rglobsn)
    endif

  end subroutine wrtout_2darray_real

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: getnum
!
! !INTERFACE:
  integer function getnum (type1d)
!
! !DESCRIPTION:
! Determines size (across all processors) from 1d type
!
! !USES:
  use decompMod, only : get_proc_global
#if (defined RTM)
  use RunoffMod, only : get_proc_rof_global
#endif
!
! !ARGUMENTS:
  implicit none
  character(len=*), intent(in) :: type1d    ! type of 1d array
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: numg        ! total number of gridcells across all processors
    integer :: numl        ! total number of landunits across all processors
    integer :: numc        ! total number of columns   across all processors
    integer :: nump        ! total number of pfts      across all processors
#if (defined RTM)
    integer :: num_roflnd  ! total number of land  runoff points across all procs
    integer :: num_rofocn  ! total number of ocean runoff points across all procs
#endif
!-----------------------------------------------------------------------
    call get_proc_global(numg, numl, numc, nump)
#if (defined RTM)
    call get_proc_rof_global(num_roflnd, num_rofocn)
#endif

    select case (type1d)
    case(nameg)
       getnum = numg
    case(namel)
       getnum = numl
    case(namec)
       getnum = numc
    case(namep)
       getnum = nump
#if (defined RTM)
    case(lndrof)
       getnum = num_roflnd
    case(ocnrof)
       getnum = num_rofocn
#endif
    case default
       write(6,*) 'getnum errror: no match for type ',trim(type1d)
       call endrun()
    end select

  end function getnum

end module iobinary
