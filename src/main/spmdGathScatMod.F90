#include <misc.h>
#include <preproc.h>

module spmdGathScatMod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: spmdGathScatMod
!
! !DESCRIPTION:
! Perform SPMD gather and scatter operations.
!
! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
  use spmdMod
  use clm_mct_mod
  use abortutils, only : endrun
!
! !PUBLIC TYPES:
  implicit none
!
! !PUBLIC MEMBER FUNCTIONS:
  public  scatter_data_from_master, gather_data_to_master, allgather_data

  interface scatter_data_from_master
     module procedure scatter_1darray_int
     module procedure scatter_1darray_real
     module procedure scatter_2darray_int
     module procedure scatter_2darray_real
     module procedure scatter_gs_1darray_int
     module procedure scatter_gs_1darray_real
  end interface

  interface gather_data_to_master
     module procedure gather_1darray_int
     module procedure gather_1darray_real
     module procedure gather_2darray_int
     module procedure gather_2darray_real
     module procedure gather_gs_1darray_int
     module procedure gather_gs_1darray_real
  end interface

  interface allgather_data
     module procedure allgather_1darray_int
     module procedure allgather_1darray_real
     module procedure allgather_2darray_int
     module procedure allgather_2darray_real
  end interface
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
!
  private spmd_compute_mpigs

!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: spmd_compute_mpigs
!
! !INTERFACE:
  subroutine spmd_compute_mpigs (clmlevel, nfact, numtot, numperproc, &
                                 displs, indexi)
!
! !DESCRIPTION:
! Compute arguments for gatherv, scatterv for vectors
!
! !USES:
    use clmtype  , only : nameg, namel, namec, namep, ocnrof, lndrof, allrof
    use decompMod, only : get_proc_bounds, get_proc_total
#if (defined RTM)
    use RunoffMod, only : get_proc_rof_bounds, get_proc_rof_total
#endif
!
! !ARGUMENTS:
    implicit none
    character(len=*), intent(in) :: clmlevel     ! type of input data
    integer, intent(in ) :: nfact                ! multiplicative factor for patches
    integer, intent(out) :: numtot               ! total number of elements (to send or recv)
    integer, intent(out) :: numperproc(0:npes-1) ! per-PE number of items to receive
    integer, intent(out) :: displs(0:npes-1)     ! per-PE displacements
    integer, intent(out) :: indexi               ! beginning array index (grid,land,col or pft)
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: pid          ! processor id
    integer :: begg, endg   ! beginning and ending gridcell index in processor
    integer :: begl, endl   ! beginning and ending landunit index in processor
    integer :: begc, endc   ! beginning and ending column index in processor
    integer :: begp, endp   ! beginning and ending pft index in processor
    integer :: ncells       ! total number of gridcells on the processor
    integer :: nlunits      ! total number of landunits on the processor
    integer :: ncols        ! total number of columns on the processor
    integer :: npfts        ! total number of pfts on the processor
    integer :: begr, endr   ! beginning and ending rtm index in processor
    integer :: nroff        ! total number of rtm cells on the processor
!----------------------------------------------------------------------

    select case (clmlevel)
    case(nameg)

       call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)
       call get_proc_total(iam, ncells, nlunits, ncols, npfts)
       numtot = ncells * nfact
       do pid = 0,npes-1
          call get_proc_total(pid, ncells, nlunits, ncols, npfts)
          numperproc(pid) = ncells * nfact
       end do
       indexi = begg

    case(namel)

       call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)
       call get_proc_total(iam, ncells, nlunits, ncols, npfts)
       numtot = nlunits * nfact
       do pid = 0,npes-1
          call get_proc_total(pid, ncells, nlunits, ncols, npfts)
          numperproc(pid) = nlunits * nfact
       end do
       indexi = begl

    case(namec)

       call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)
       call get_proc_total(iam, ncells, nlunits, ncols, npfts)
       numtot = ncols * nfact
       do pid = 0,npes-1
          call get_proc_total(pid, ncells, nlunits, ncols, npfts)
          numperproc(pid) = ncols * nfact
       end do
       indexi = begc

    case(namep)

       call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)
       call get_proc_total(iam, ncells, nlunits, ncols, npfts)
       numtot = npfts * nfact
       do pid = 0,npes-1
          call get_proc_total (pid, ncells, nlunits, ncols, npfts)
          numperproc(pid) = npfts * nfact
       end do
       indexi = begp

#if (defined RTM)
    case(allrof)

       call get_proc_rof_bounds(begr, endr)
       call get_proc_rof_total(iam, nroff)
       numtot = nroff * nfact
       do pid = 0,npes-1
          call get_proc_rof_total(pid, nroff)
          numperproc(pid) = nroff * nfact
       end do
       indexi = begr

    case(lndrof)

       call get_proc_rof_bounds(begr, endr)
       call get_proc_rof_total(iam, nroff)
       numtot = nroff * nfact
       do pid = 0,npes-1
          call get_proc_rof_total(pid, nroff)
          numperproc(pid) = nroff * nfact
       end do
       indexi = begr

    case(ocnrof)

       call get_proc_rof_bounds(begr, endr)
       call get_proc_rof_total(iam, nroff)
       numtot = nroff * nfact
       do pid = 0,npes-1
          call get_proc_rof_total(pid, nroff)
          numperproc(pid) = nroff * nfact
       end do
       indexi = begr
#endif

    case default

       write(6,*) 'COMPUTE_MPIGS: Invalid expansion character: ',trim(clmlevel)
       call endrun

    end select

    displs(0) = 0
    do pid = 1,npes-1
       displs(pid) = displs(pid-1) + numperproc(pid-1)
    end do

  end subroutine spmd_compute_mpigs

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: scatter_gs_1darray_int
!
! !INTERFACE:
  subroutine scatter_gs_1darray_int (ilocal, iglobal, gsmap, perm, lbeg, lend)
!
! !DESCRIPTION:
! Wrapper routine to scatter integer 1d array
!
! !USES:
!
! !ARGUMENTS:
    implicit none
    integer, pointer              :: ilocal(:)   !local read data (out)
    integer, pointer              :: iglobal(:)  !global read data (in)
    type(mct_gsMap) , intent(in ) :: gsmap       !global seg map
    integer, pointer              :: perm(:)     !gsmap permuter
    integer, optional,intent(in ) :: lbeg,lend   !beg/end indices of rlocal
!
! !REVISION HISTORY:
! Author: T Craig
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: n,lsize,lb,ub
    type(mct_aVect) :: AVi, AVo   ! attribute vectors
    integer,pointer :: ivect(:) ! local vector

!-----------------------------------------------------------------------

  if (masterproc) then
     lsize = size(iglobal)
     call mct_aVect_init(AVi,iList='array',lsize=lsize)
     call mct_aVect_importIattr(AVi,"array",iglobal,lsize)
  endif
  call mct_aVect_scatter(AVi, AVo, gsmap, 0, mpicom)
  call mct_aVect_unpermute(AVo, perm, dieWith='scatter_gs_1darray_real')

  lsize = size(ilocal)
  lb = 1
  ub = lsize
  if (present(lbeg)) then
     lb = lbeg
  endif
  if (present(lend)) then
     ub = lend
  endif

  allocate(ivect(lsize))
  call mct_aVect_exportIattr(AVo,"array",ivect,lsize)

  do n = lb,ub
     ilocal(n) = ivect(n-lb+1)
  enddo

  deallocate(ivect)
  if (masterproc) then
     call mct_aVect_clean(AVi)
  endif
  call mct_aVect_clean(AVo)

  end subroutine scatter_gs_1darray_int

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: scatter_gs_1darray_real
!
! !INTERFACE:
  subroutine scatter_gs_1darray_real (rlocal, rglobal, gsmap, perm, lbeg, lend)
!
! !DESCRIPTION:
! Wrapper routine to scatter real 1d array
!
! !USES:
!
! !ARGUMENTS:
    implicit none
    real(r8), pointer             :: rlocal(:)   !local read data (out)
    real(r8), pointer             :: rglobal(:)  !global read data (in)
    type(mct_gsMap)  ,intent(in ) :: gsmap       !global seg map
    integer,  pointer             :: perm(:)     !gsmap permuter
    integer, optional,intent(in ) :: lbeg,lend   !beg/end indices of rlocal
!
! !REVISION HISTORY:
! Author: T Craig
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: n,lsize,lb,ub
    type(mct_aVect) :: AVi, AVo   ! attribute vectors
    real(r8),pointer :: rvect(:) ! local vector

!-----------------------------------------------------------------------

  if (masterproc) then
     lsize = size(rglobal)
     call mct_aVect_init(AVi,rList='array',lsize=lsize)
     call mct_aVect_importRattr(AVi,"array",rglobal,lsize)
  endif
  call mct_aVect_scatter(AVi, AVo, gsmap, 0, mpicom)
  call mct_aVect_unpermute(AVo, perm, dieWith='scatter_gs_1darray_real')

  lsize = size(rlocal)
  lb = 1
  ub = lsize
  if (present(lbeg)) then
     lb = lbeg
  endif
  if (present(lend)) then
     ub = lend
  endif

  allocate(rvect(lsize))
  call mct_aVect_exportRattr(AVo,"array",rvect,lsize)

  do n = lb,ub
     rlocal(n) = rvect(n-lb+1)
  enddo

  deallocate(rvect)
  if (masterproc) then
     call mct_aVect_clean(AVi)
  endif
  call mct_aVect_clean(AVo)

  end subroutine scatter_gs_1darray_real

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: scatter_1darray_int
!
! !INTERFACE:
  subroutine scatter_1darray_int (ilocal, iglobal, clmlevel)
!
! !DESCRIPTION:
! Wrapper routine to scatter integer 1d array
!
! !USES:
!
! !ARGUMENTS:
    implicit none
    integer,pointer,  dimension(:) :: ilocal(:)   !local read data
    integer,pointer,  dimension(:) :: iglobal(:)  !global read data
    character(len=*), intent(in) :: clmlevel      !type of input data
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: ier                    !return code
    integer :: beg                    !temporary
    integer :: numsendv(0:npes-1)     !vector of items to be sent
    integer :: displsv(0:npes-1)      !displacement vector
    integer :: numrecv                !number of items to be received
!-----------------------------------------------------------------------
    call spmd_compute_mpigs (clmlevel, 1, numrecv, numsendv, displsv, beg)
    if (masterproc) then
       call mpi_scatterv (iglobal, numsendv, displsv, MPI_INTEGER, &
            ilocal(beg), numrecv , MPI_INTEGER , 0, mpicom, ier)
    else
       call mpi_scatterv (0, numsendv, displsv, MPI_INTEGER, &
            ilocal(beg), numrecv , MPI_INTEGER , 0, mpicom, ier)
    endif
    if (ier /= 0) then
       write(6,*)'scatter_1darray_int error: ',ier
       call endrun
    endif
  end subroutine scatter_1darray_int

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: scatter_1darray_real
!
! !INTERFACE:
  subroutine scatter_1darray_real (rlocal, rglobal, clmlevel)
!
! !DESCRIPTION:
! Wrapper routine to scatter 1d real array from master processor
!
! !ARGUMENTS:
    implicit none
    real(r8), pointer, dimension(:) :: rlocal  !local read data
    real(r8), pointer, dimension(:) :: rglobal !global read data
    character(len=*), intent(in) :: clmlevel   !input data type
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: ier                             !return code
    integer :: beg                             !temporaries
    integer :: numsendv(0:npes-1)              !vector of items to be sent
    integer :: displsv(0:npes-1)               !displacement vector
    integer :: numrecv                         !number of items to be received
!-----------------------------------------------------------------------
    call spmd_compute_mpigs (clmlevel, 1, numrecv, numsendv, displsv, beg)
    if (masterproc) then
       call mpi_scatterv (rglobal, numsendv, displsv, MPI_REAL8, &
            rlocal(beg), numrecv , MPI_REAL8 , 0, mpicom, ier)
    else
       call mpi_scatterv (0._r8, numsendv, displsv, MPI_REAL8, &
            rlocal(beg), numrecv , MPI_REAL8 , 0, mpicom, ier)
    endif
    if (ier/=0 ) then
       write(6,*)'scatter_1darray_real error: ',ier
       call endrun
    endif
  end subroutine scatter_1darray_real

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: scatter_2darray_int
!
! !INTERFACE:
  subroutine scatter_2darray_int (ilocal, iglobal, clmlevel)
!
! !DESCRIPTION:
! Wrapper routine to scatter 2d integer array from master processor
!
! !ARGUMENTS:
    implicit none
    integer, pointer, dimension(:,:) :: ilocal  !local read data
    integer, pointer, dimension(:,:) :: iglobal !global read data
    character(len=*), intent(in) :: clmlevel    !type of input data
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: ier                   !error status
    integer :: ndim1                 !size of first dimension
    integer :: l1                    !lower bound of first dimension
    integer :: beg                   !temporaries
    integer :: numsendv(0:npes-1)    !vector of items to be sent
    integer :: displsv(0:npes-1)     !displacement vector
    integer :: numrecv               !number of items to be received
!-----------------------------------------------------------------------
    if (masterproc) then
       if (lbound(ilocal,dim=1) /= lbound(iglobal,dim=1)) then
          write(6,*)' lower bounds of global and local input arrays do not match'
          write(6,*)' l1 local = ',lbound(ilocal,dim=1)
          write(6,*)' l1 global= ',lbound(iglobal,dim=1)
          call endrun
       endif
    endif
    l1 = lbound(ilocal,dim=1)
    ndim1 = size(ilocal,dim=1)
    call spmd_compute_mpigs (clmlevel, ndim1, numrecv, numsendv, displsv, beg)
    if (masterproc) then
       call mpi_scatterv (iglobal(l1,1), numsendv, displsv, MPI_INTEGER, &
            ilocal(l1,beg), numrecv ,  MPI_INTEGER, 0, mpicom, ier)
    else
       call mpi_scatterv (0, numsendv, displsv, MPI_INTEGER, &
            ilocal(l1,beg), numrecv ,  MPI_INTEGER, 0, mpicom, ier)
    endif
    if (ier /= 0) then
       write(6,*)'scatter_2darray_int error: ',ier
       call endrun
    endif
  end subroutine scatter_2darray_int

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: scatter_2darray_real
!
! !INTERFACE:
  subroutine scatter_2darray_real (rlocal, rglobal, clmlevel)
!
! !DESCRIPTION:
! Wrapper routine to scatter 2d integer array from master processor
!
! !ARGUMENTS:
    implicit none
    real(r8), pointer, dimension(:,:) :: rlocal  !local read data
    real(r8), pointer, dimension(:,:) :: rglobal !global read data
    character(len=*), intent(in) :: clmlevel     !type of input data
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: ier                            !return code
    integer :: ndim1                          !size of second dimension
    integer :: l1                             !lower,upper bounds of input arrays
    integer :: beg                            !temporaries
    integer :: numsendv(0:npes-1)             !vector of items to be sent
    integer :: displsv(0:npes-1)              !displacement vector
    integer :: numrecv                        !number of items to be received
!-----------------------------------------------------------------------
    if (masterproc) then
       if (lbound(rlocal,dim=1) /= lbound(rglobal,dim=1)) then
          write(6,*)' lower bounds of global and local input arrays do not match'
          write(6,*)' l1 local = ',lbound(rlocal,dim=1)
          write(6,*)' l1 global= ',lbound(rglobal,dim=1)
          call endrun
       endif
    endif
    l1 = lbound(rlocal,dim=1)
    ndim1 = size(rlocal,dim=1)
    call spmd_compute_mpigs (clmlevel, ndim1, numrecv, numsendv, displsv, beg)
    if (masterproc) then
       call mpi_scatterv (rglobal(l1,1), numsendv, displsv, MPI_REAL8, &
            rlocal(l1,beg), numrecv ,  MPI_REAL8, 0, mpicom, ier)
    else
       call mpi_scatterv (0._r8, numsendv, displsv, MPI_REAL8, &
            rlocal(l1,beg), numrecv ,  MPI_REAL8, 0, mpicom, ier)
    endif
    if (ier /= 0) then
       write(6,*)'scatter_2darray_real error: ',ier
       call endrun
    endif
  end subroutine scatter_2darray_real

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: gather_gs_1darray_int
!
! !INTERFACE:
  subroutine gather_gs_1darray_int (ilocal, iglobal, gsmap, perm, lbeg, lend)
!
! !DESCRIPTION:
! Wrapper routine to gather integer 1d array
!
! !USES:
!
! !ARGUMENTS:
    implicit none
    integer, pointer              :: ilocal(:)   !local read data (in)
    integer, pointer              :: iglobal(:)  !global read data (out)
    type(mct_gsMap) , intent(in ) :: gsmap       !global seg map
    integer, pointer              :: perm(:)     !gsmap permuter
    integer, optional,intent(in ) :: lbeg,lend   !beg/end indices of rlocal
!
! !REVISION HISTORY:
! Author: T Craig
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: n,lsize,lb,ub
    type(mct_aVect) :: AVi, AVo   ! attribute vectors
    integer,pointer :: ivect(:) ! local vector

!-----------------------------------------------------------------------

  lsize = size(ilocal)

  allocate(ivect(lsize))
  lb = 1
  ub = lsize
  if (present(lbeg)) then
     lb = lbeg
  endif
  if (present(lend)) then
     ub = lend
  endif

  do n = lb,ub
     ivect(n-lb+1) = ilocal(n)
  enddo

  call mct_aVect_init(AVi,iList='array',lsize=lsize)
  call mct_aVect_importIattr(AVi,"array",ivect,lsize)

  call mct_aVect_permute(AVi, perm, dieWith='gather_gs_1darray_real')
  call mct_aVect_gather(AVi, AVo, gsmap, 0, mpicom)

  lsize = size(iglobal)

  if (masterproc) then
     call mct_aVect_exportIattr(AVo,"array",iglobal,lsize)
  endif

  deallocate(ivect)
  call mct_aVect_clean(AVi)
  call mct_aVect_clean(AVo)

  end subroutine gather_gs_1darray_int

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: gather_gs_1darray_real
!
! !INTERFACE:
  subroutine gather_gs_1darray_real (rlocal, rglobal, gsmap, perm, lbeg, lend)
!
! !DESCRIPTION:
! Wrapper routine to gather real 1d array
!
! !USES:
!
! !ARGUMENTS:
    implicit none
    real(r8), pointer             :: rlocal(:)   !local read data (in)
    real(r8), pointer             :: rglobal(:)  !global read data (out)
    type(mct_gsMap)  ,intent(in ) :: gsmap       !global seg map
    integer,  pointer             :: perm(:)     !gsmap permuter
    integer, optional,intent(in ) :: lbeg,lend   !beg/end indices of rlocal
!
! !REVISION HISTORY:
! Author: T Craig
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: n,lsize,lb,ub
    type(mct_aVect) :: AVi, AVo   ! attribute vectors
    real(r8),pointer :: rvect(:) ! local vector

!-----------------------------------------------------------------------


  lsize = size(rlocal)

  allocate(rvect(lsize))
  lb = 1
  ub = lsize
  if (present(lbeg)) then
     lb = lbeg
  endif
  if (present(lend)) then
     ub = lend
  endif

  do n = lb,ub
     rvect(n-lb+1) = rlocal(n)
  enddo

  call mct_aVect_init(AVi,rList='array',lsize=lsize)
  call mct_aVect_importRattr(AVi,"array",rvect,lsize)

  call mct_aVect_permute(AVi, perm, dieWith='gather_gs_1darray_real')
  call mct_aVect_gather(AVi, AVo, gsmap, 0, mpicom)

  lsize = size(rglobal)

  if (masterproc) then
     call mct_aVect_exportRattr(AVo,"array",rglobal,lsize)
  endif

  deallocate(rvect)
  call mct_aVect_clean(AVi)
  call mct_aVect_clean(AVo)

  end subroutine gather_gs_1darray_real

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: gather_1darray_int
!
! !INTERFACE:
  subroutine gather_1darray_int (ilocal, iglobal, clmlevel)
!
! !DESCRIPTION:
! Wrapper routine to gather 1d integer array on master processor
!
! !ARGUMENTS:
    implicit none
    integer, pointer, dimension(:) :: ilocal     !output data
    integer, pointer, dimension(:) :: iglobal    !output data
    character(len=*), intent(in) :: clmlevel     !input data type
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: ier                    !errorcode
    integer :: beg                    !temporary
    integer :: numrecvv(0:npes-1)     !vector of items to be received
    integer :: displsv(0:npes-1)      !displacement vector
    integer :: numsend                !number of items to be sent
!-----------------------------------------------------------------------
    call spmd_compute_mpigs (clmlevel, 1, numsend, numrecvv, displsv, beg)
    if (masterproc) then
       call mpi_gatherv (ilocal(beg), numsend , MPI_INTEGER, &
            iglobal, numrecvv, displsv, MPI_INTEGER, 0, mpicom, ier)
    else
       call mpi_gatherv (ilocal(beg), numsend , MPI_INTEGER, &
            0, numrecvv, displsv, MPI_INTEGER, 0, mpicom, ier)
    endif
    if (ier/=0 ) then
       write(6,*)'gather_1darray_int error: ',ier
       call endrun
    endif
  end subroutine gather_1darray_int

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: gather_1darray_real
!
! !INTERFACE:
  subroutine gather_1darray_real (rlocal, rglobal, clmlevel)
!
! !DESCRIPTION:
! Wrapper routine to gather 1d real array on master processor
!
! !ARGUMENTS:
    implicit none
    real(r8), pointer, dimension(:) :: rlocal    !output data
    real(r8), pointer, dimension(:) :: rglobal   !output data
    character(len=*), intent(in) :: clmlevel     !input data type
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: ier                    !return code
    integer :: beg                    !temporary
    integer :: numrecvv(0:npes-1)     !vector of items to be received
    integer :: displsv(0:npes-1)      !displacement vector
    integer :: numsend                !number of items to be sent
!-----------------------------------------------------------------------
    call spmd_compute_mpigs (clmlevel, 1, numsend, numrecvv, displsv, beg)
    if (masterproc) then
       call mpi_gatherv (rlocal(beg), numsend , MPI_REAL8, &
            rglobal, numrecvv, displsv, MPI_REAL8, 0, mpicom, ier)
    else
       call mpi_gatherv (rlocal(beg), numsend , MPI_REAL8, &
            0._r8, numrecvv, displsv, MPI_REAL8, 0, mpicom, ier)
    endif
    if (ier /= 0) then
       write(6,*)'gather_1darray_real error: ',ier
       call endrun
    endif
  end subroutine gather_1darray_real

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: gather_2darray_int
!
! !INTERFACE:
  subroutine gather_2darray_int (ilocal, iglobal, clmlevel)
!
! !DESCRIPTION:
! Wrapper routine to gather 2d integer array on master processor
!
! !ARGUMENTS:
    implicit none
    integer, pointer, dimension(:,:) :: ilocal   !read data
    integer, pointer, dimension(:,:) :: iglobal  !global data
    character(len=*), intent(in) :: clmlevel     !type of input data
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: ier                   !return code
    integer :: ndim1                 !size of second dimension
    integer :: l1                    !lower bounds of input arrays
    integer :: beg                   !temporaries
    integer :: numrecvv(0:npes-1)    !vector of items to be received
    integer :: displsv(0:npes-1)     !displacement vector
    integer :: numsend               !number of items to be sent
!-----------------------------------------------------------------------
    if (masterproc) then
       if (lbound(ilocal,dim=1) /= lbound(iglobal,dim=1)) then
          write(6,*)' lower bounds of global and local input arrays do not match'
          write(6,*)' l1 local = ',lbound(ilocal,dim=1)
          write(6,*)' l1 global= ',lbound(iglobal,dim=1)
          call endrun
       endif
    endif
    l1 = lbound(ilocal,dim=1)
    ndim1 = size(ilocal,dim=1)
    call spmd_compute_mpigs (clmlevel, ndim1, numsend, numrecvv, displsv, beg)
    if (masterproc) then
       call mpi_gatherv (ilocal(l1,beg), numsend , MPI_INTEGER, &
            iglobal(l1,1), numrecvv, displsv, MPI_INTEGER, 0, mpicom, ier)
    else
       call mpi_gatherv (ilocal(l1,beg), numsend , MPI_INTEGER, &
            0, numrecvv, displsv, MPI_INTEGER, 0, mpicom, ier)
    endif
    if (ier /= 0) then
       write(6,*)'gather_2darray_int error: ',ier
       call endrun
    endif
  end subroutine gather_2darray_int

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: gather_2darray_real
!
! !INTERFACE:
  subroutine gather_2darray_real (rlocal, rglobal, clmlevel)
!
! !DESCRIPTION:
! Wrapper routine to gather 2d real array
!
! !ARGUMENTS:
    implicit none
    real(r8), pointer, dimension(:,:) :: rlocal   !local data
    real(r8), pointer, dimension(:,:) :: rglobal  !global data
    character(len=*), intent(in) :: clmlevel      !type of input data
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: ier                    !return code
    integer :: ndim1                  !size of second dimension
    integer :: l1                     !lower bounds of input arrays
    integer :: beg                    !temporary
    integer :: numrecvv(0:npes-1)     !vector of items to be received
    integer :: displsv(0:npes-1)      !displacement vector
    integer :: numsend                !number of items to be sent
!-----------------------------------------------------------------------
    if (masterproc) then
       if (lbound(rlocal,dim=1) /= lbound(rglobal,dim=1)) then
          write(6,*)' lower bounds of global and local input arrays do not match'
          write(6,*)' l1 local = ',lbound(rlocal,dim=1)
          write(6,*)' l1 global= ',lbound(rglobal,dim=1)
          call endrun
       endif
    endif
    l1 = lbound(rlocal,dim=1)
    ndim1 = size(rlocal,dim=1)
    call spmd_compute_mpigs (clmlevel, ndim1, numsend, numrecvv, displsv, beg)
    if (masterproc) then
       call mpi_gatherv (rlocal(l1,beg), numsend , MPI_REAL8, &
            rglobal(l1,1), numrecvv, displsv, MPI_REAL8, 0, mpicom, ier)
    else
       call mpi_gatherv (rlocal(l1,beg), numsend , MPI_REAL8, &
            0._r8, numrecvv, displsv, MPI_REAL8, 0, mpicom, ier)
    endif
    if (ier /= 0) then
       write(6,*)'gather_2darray_real error: ',ier
       call endrun
    endif
  end subroutine gather_2darray_real

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: allgather_1darray_int
!
! !INTERFACE:
  subroutine allgather_1darray_int (ilocal, iglobal, clmlevel)
!
! !DESCRIPTION:
! Wrapper routine to perform an allgatherv of 1d integer array
!
! !ARGUMENTS:
    implicit none
    integer, pointer, dimension(:) :: ilocal     !output data
    integer, pointer, dimension(:) :: iglobal    !output data
    character(len=*), intent(in) :: clmlevel     !input data type
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: ier                    !errorcode
    integer :: beg                    !temporary
    integer :: numrecvv(0:npes-1)     !vector of items to be received
    integer :: displsv(0:npes-1)      !displacement vector
    integer :: numsend                !number of items to be sent
!-----------------------------------------------------------------------
    call spmd_compute_mpigs (clmlevel, 1, numsend, numrecvv, displsv, beg)
    call mpi_allgatherv (ilocal(beg), numsend , MPI_INTEGER, &
                         iglobal, numrecvv, displsv, MPI_INTEGER, mpicom, ier)
    if (ier/=0 ) then
       write(6,*)'gather_1darray_int error: ',ier
       call endrun
    endif
  end subroutine allgather_1darray_int

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: allgather_1darray_real
!
! !INTERFACE:
  subroutine allgather_1darray_real (rlocal, rglobal, clmlevel)
!
! !DESCRIPTION:
! Wrapper routine to perform an allgatherv of 1d real array
!
! !ARGUMENTS:
    implicit none
    real(r8), pointer, dimension(:) :: rlocal    !output data
    real(r8), pointer, dimension(:) :: rglobal   !output data
    character(len=*), intent(in) :: clmlevel     !input data type
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: ier                    !return code
    integer :: beg                    !temporary
    integer :: numrecvv(0:npes-1)     !vector of items to be received
    integer :: displsv(0:npes-1)      !displacement vector
    integer :: numsend                !number of items to be sent
!-----------------------------------------------------------------------
    call spmd_compute_mpigs (clmlevel, 1, numsend, numrecvv, displsv, beg)
    call mpi_allgatherv (rlocal(beg), numsend , MPI_REAL8, &
                         rglobal, numrecvv, displsv, MPI_REAL8, mpicom, ier)
    if (ier /= 0) then
       write(6,*)'gather_1darray_real error: ',ier
       call endrun
    endif
  end subroutine allgather_1darray_real

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: allgather_2darray_int
!
! !INTERFACE:
  subroutine allgather_2darray_int (ilocal, iglobal, clmlevel)
!
! !DESCRIPTION:
! Wrapper routine to perform an allgatherv of 2d integer array
!
! !ARGUMENTS:
    implicit none
    integer, pointer, dimension(:,:) :: ilocal   !read data
    integer, pointer, dimension(:,:) :: iglobal  !global data
    character(len=*), intent(in) :: clmlevel     !type of input data
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: ier                   !return code
    integer :: ndim1                 !size of second dimension
    integer :: l1                    !lower bounds of input arrays
    integer :: beg                   !temporaries
    integer :: numrecvv(0:npes-1)    !vector of items to be received
    integer :: displsv(0:npes-1)     !displacement vector
    integer :: numsend               !number of items to be sent
!-----------------------------------------------------------------------
    if (lbound(ilocal,dim=1) /= lbound(iglobal,dim=1)) then
       write(6,*)' lower bounds of global and local input arrays do not match'
       write(6,*)' l1 local = ',lbound(ilocal,dim=1)
       write(6,*)' l1 global= ',lbound(iglobal,dim=1)
       call endrun
    endif
    l1 = lbound(ilocal,dim=1)
    ndim1 = size(ilocal,dim=1)
    call spmd_compute_mpigs (clmlevel, ndim1, numsend, numrecvv, displsv, beg)
    call mpi_allgatherv (ilocal(l1,beg), numsend , MPI_INTEGER, &
                         iglobal(l1,1), numrecvv, displsv, MPI_INTEGER, mpicom, ier)
    if (ier /= 0) then
       write(6,*)'gather_2darray_int error: ',ier
       call endrun
    endif
  end subroutine allgather_2darray_int

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: allgather_2darray_real
!
! !INTERFACE:
  subroutine allgather_2darray_real (rlocal, rglobal, clmlevel)
!
! !DESCRIPTION:
! Wrapper routine to perform an allgatherv of 2d real array
!
! !ARGUMENTS:
    implicit none
    real(r8), pointer, dimension(:,:) :: rlocal   !local data
    real(r8), pointer, dimension(:,:) :: rglobal  !global data
    character(len=*), intent(in) :: clmlevel      !type of input data
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: ier                    !return code
    integer :: ndim1                  !size of second dimension
    integer :: l1                     !lower bounds of input arrays
    integer :: beg                    !temporary
    integer :: numrecvv(0:npes-1)     !vector of items to be received
    integer :: displsv(0:npes-1)      !displacement vector
    integer :: numsend                !number of items to be sent
!-----------------------------------------------------------------------
    if (lbound(rlocal,dim=1) /= lbound(rglobal,dim=1)) then
       write(6,*)' lower bounds of global and local input arrays do not match'
       write(6,*)' l1 local = ',lbound(rlocal,dim=1)
       write(6,*)' l1 global= ',lbound(rglobal,dim=1)
       call endrun
    endif
    l1 = lbound(rlocal,dim=1)
    ndim1 = size(rlocal,dim=1)
    call spmd_compute_mpigs (clmlevel, ndim1, numsend, numrecvv, displsv, beg)
    call mpi_allgatherv (rlocal(l1,beg), numsend , MPI_REAL8, &
                         rglobal(l1,1), numrecvv, displsv, MPI_REAL8, mpicom, ier)
    if (ier /= 0) then
       write(6,*)'gather_2darray_real error: ',ier
       call endrun
    endif
  end subroutine allgather_2darray_real

end module spmdGathScatMod
