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
  use clm_varcon, only: spval, ispval
  use decompMod, only : get_clmlevel_gsmap
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
  end interface

  interface gather_data_to_master
     module procedure gather_1darray_int
     module procedure gather_1darray_real
     module procedure gather_2darray_int
     module procedure gather_2darray_real
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
  private scatter_gs_1darray_int
  private scatter_gs_1darray_real
  private gather_gs_1darray_int
  private gather_gs_1darray_real

!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: scatter_gs_1darray_int
!
! !INTERFACE:
  subroutine scatter_gs_1darray_int (ilocal, iglobal, gsmap, perm)
!
! !DESCRIPTION:
! Wrapper routine to scatter integer 1d array
!
! !USES:
!
! !ARGUMENTS:
    implicit none
    integer, pointer              :: ilocal(:)   !local readl2ddata (out)
    integer, pointer              :: iglobal(:)  !global read data (in)
    type(mct_gsMap) , intent(in ) :: gsmap       !global seg map
    integer, pointer              :: perm(:)     !gsmap permuter
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
    character(len=*),parameter :: subname = 'scatter_gs_1darray_int'

!-----------------------------------------------------------------------

  if (masterproc) then
     lsize = size(iglobal)
     call mct_aVect_init(AVi,iList='array',lsize=lsize)
     call mct_aVect_importIattr(AVi,"array",iglobal,lsize)
  endif
  call mct_aVect_scatter(AVi, AVo, gsmap, 0, mpicom)
  call mct_aVect_unpermute(AVo, perm, dieWith=subname)

  lsize = size(ilocal)
  lb = lbound(ilocal, dim=1)
  ub = ubound(ilocal, dim=1)

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
  subroutine scatter_gs_1darray_real (rlocal, rglobal, gsmap, perm)
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
    character(len=*),parameter :: subname = 'scatter_gs_1darray_real'

!-----------------------------------------------------------------------

  if (masterproc) then
     lsize = size(rglobal)
     call mct_aVect_init(AVi,rList='array',lsize=lsize)
     call mct_aVect_importRattr(AVi,"array",rglobal,lsize)
  endif
  call mct_aVect_scatter(AVi, AVo, gsmap, 0, mpicom)
  call mct_aVect_unpermute(AVo, perm, dieWith=subname)

  lsize = size(rlocal)
  lb = lbound(rlocal, dim=1)
  ub = ubound(rlocal, dim=1)

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
    integer, pointer              :: ilocal(:)   !local read data (in)
    integer, pointer              :: iglobal(:)  !global read data (out)
    character(len=*), intent(in) :: clmlevel     !type of input data
!
! !REVISION HISTORY:
! Author: T Craig
!
!EOP
!
! !LOCAL VARIABLES:
    type(mct_gsMap) , pointer :: gsmap       !global seg map
    integer, pointer,dimension(:) :: perm
    character(len=*),parameter :: subname = 'scatter_1darray_int'

!-----------------------------------------------------------------------

  call get_clmlevel_gsmap(clmlevel,gsmap,perm)
  call scatter_gs_1darray_int(ilocal,iglobal,gsmap,perm)

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
! Wrapper routine to scatter real 1d array
!
! !USES:
!
! !ARGUMENTS:
    implicit none
    real(r8), pointer             :: rlocal(:)   !local read data (in)
    real(r8), pointer             :: rglobal(:)  !global read data (out)
    character(len=*), intent(in) :: clmlevel     !type of input data
!
! !REVISION HISTORY:
! Author: T Craig
!
!EOP
!
! !LOCAL VARIABLES:
    type(mct_gsMap),pointer       :: gsmap       !global seg map
    integer, pointer,dimension(:) :: perm
    character(len=*),parameter :: subname = 'scatter_1darray_real'

!-----------------------------------------------------------------------

  call get_clmlevel_gsmap(clmlevel,gsmap,perm)
  call scatter_gs_1darray_real(rlocal,rglobal,gsmap,perm)

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
! Wrapper routine to scatter integer 2d array
!
! !USES:
!
! !ARGUMENTS:
    implicit none
    integer, pointer              :: ilocal(:,:)   !local read data (in)
    integer, pointer              :: iglobal(:,:)  !global read data (out)
    character(len=*), intent(in) :: clmlevel     !type of input data
!
! !REVISION HISTORY:
! Author: T Craig
!
!EOP
!
! !LOCAL VARIABLES:
    type(mct_gsMap) , pointer :: gsmap       !global seg map
    integer, pointer,dimension(:) :: perm
    integer, pointer :: locarr(:)
    integer, pointer :: gloarr(:)
    integer :: lbg,ubg,lbl,ubl,sizg,sizl,n
    character(len=*),parameter :: subname = 'scatter_2darray_int'

!-----------------------------------------------------------------------

  call get_clmlevel_gsmap(clmlevel,gsmap,perm)

  lbl = lbound(ilocal, dim=2)
  ubl = ubound(ilocal, dim=2)
  sizl = size(ilocal ,dim=1)
  allocate(locarr(sizl))

  if (masterproc) then
     lbg = lbound(iglobal,dim=2)
     ubg = ubound(iglobal,dim=2)
     if (ubg-lbg /= ubl-lbl) then
        write(6,*) trim(subname),' error second dim ',lbg,ubg,lbl,ubl
        call endrun()
     endif
     sizg = size(iglobal,dim=1)
     allocate(gloarr(sizg))
  endif

  do n = lbl,ubl
     if (masterproc) gloarr(:) = iglobal(:,lbg+n-lbl)
     call scatter_gs_1darray_int(locarr,gloarr,gsmap,perm)
     ilocal(:,n) = locarr(:)
  enddo

  deallocate(locarr)
  if (masterproc) deallocate(gloarr)

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
! Wrapper routine to scatter real 2d array
!
! !USES:
!
! !ARGUMENTS:
    implicit none
    real(r8), pointer             :: rlocal(:,:)   !local read data (in)
    real(r8), pointer             :: rglobal(:,:)  !global read data (out)
    character(len=*), intent(in) :: clmlevel     !type of input data
!
! !REVISION HISTORY:
! Author: T Craig
!
!EOP
!
! !LOCAL VARIABLES:
    type(mct_gsMap),pointer       :: gsmap       !global seg map
    integer, pointer,dimension(:) :: perm
    real(r8), pointer :: locarr(:)
    real(r8), pointer :: gloarr(:)
    integer :: lbg,ubg,lbl,ubl,sizg,sizl,n
    character(len=*),parameter :: subname = 'scatter_2darray_real'

!-----------------------------------------------------------------------

  call get_clmlevel_gsmap(clmlevel,gsmap,perm)
  lbl = lbound(rlocal, dim=2)
  ubl = ubound(rlocal, dim=2)
  sizl = size(rlocal ,dim=1)
  allocate(locarr(sizl))

  if (masterproc) then
     lbg = lbound(rglobal,dim=2)
     ubg = ubound(rglobal,dim=2)
     if (ubg-lbg /= ubl-lbl) then
        write(6,*) trim(subname),' error second dim ',lbg,ubg,lbl,ubl
        call endrun()
     endif
     sizg = size(rglobal,dim=1)
     allocate(gloarr(sizg))
  endif

  do n = lbl,ubl
     if (masterproc) gloarr(:) = rglobal(:,lbg+n-lbl)
     call scatter_gs_1darray_real(locarr,gloarr,gsmap,perm)
     rlocal(:,n) = locarr(:)
  enddo

  deallocate(locarr)
  if (masterproc) deallocate(gloarr)

  end subroutine scatter_2darray_real

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: gather_gs_1darray_int
!
! !INTERFACE:
  subroutine gather_gs_1darray_int (ilocal, iglobal, gsmap, perm, imissing)
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
    integer, optional,intent(in ) :: imissing    !missing value
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
    character(len=*),parameter :: subname = 'gather_gs_1darray_int'

!-----------------------------------------------------------------------

  lsize = size(ilocal)
  lb = lbound(ilocal, dim=1)
  ub = ubound(ilocal, dim=1)

  allocate(ivect(lsize))

  if (present(imissing)) then
     call mct_aVect_init(AVi,iList='array:mask',lsize=lsize)
  else
     call mct_aVect_init(AVi,iList='array',lsize=lsize)
  endif

  do n = lb,ub
     ivect(n-lb+1) = ilocal(n)
  enddo
  call mct_aVect_importIattr(AVi,"array",ivect,lsize)

  if (present(imissing)) then
     do n = lb,ub
        ivect(n-lb+1) = 1
     enddo
     call mct_aVect_importIattr(AVi,"mask",ivect,lsize)
  endif

  deallocate(ivect)

  call mct_aVect_permute(AVi, perm, dieWith=subname)
  if (present(imissing)) then
! tcx wait for update in mct, then get rid of "mask"
!     call mct_aVect_gather(AVi, AVo, gsmap, 0, mpicom, imissing = imissing)
     call mct_aVect_gather(AVi, AVo, gsmap, 0, mpicom)
  else
     call mct_aVect_gather(AVi, AVo, gsmap, 0, mpicom)
  endif

  if (masterproc) then
     lsize = size(iglobal)
     call mct_aVect_exportIattr(AVo,"array",iglobal,lsize)
     if (present(imissing)) then
        allocate(ivect(lsize))
        call mct_aVect_exportIattr(AVo,"mask",ivect,lsize)
        do n = 1,lsize
           if (ivect(n) == 0) iglobal(n) = imissing
        enddo
        deallocate(ivect)
     endif
     call mct_aVect_clean(AVo)
  endif

  call mct_aVect_clean(AVi)

  end subroutine gather_gs_1darray_int

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: gather_gs_1darray_real
!
! !INTERFACE:
  subroutine gather_gs_1darray_real (rlocal, rglobal, gsmap, perm, missing)
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
    real(r8),optional,intent(in ) :: missing     !missing value
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
    integer ,pointer :: ivect(:) ! local vector
    character(len=*),parameter :: subname = 'gather_gs_1darray_real'

!-----------------------------------------------------------------------

  lsize = size(rlocal)
  lb = lbound(rlocal, dim=1)
  ub = ubound(rlocal, dim=1)

  allocate(rvect(lsize))

  if (present(missing)) then
     call mct_aVect_init(AVi,rList='array',iList='mask',lsize=lsize)
  else
     call mct_aVect_init(AVi,rList='array',lsize=lsize)
  endif

  do n = lb,ub
     rvect(n-lb+1) = rlocal(n)
  enddo
  call mct_aVect_importRattr(AVi,"array",rvect,lsize)

  if (present(missing)) then
     allocate(ivect(lsize))
     do n = lb,ub
        ivect(n-lb+1) = 1
     enddo
     call mct_aVect_importIattr(AVi,"mask",ivect,lsize)
     deallocate(ivect)
  endif

  deallocate(rvect)

  call mct_aVect_permute(AVi, perm, dieWith=subname)
  if (present(missing)) then
! tcx wait for update in mct, then get rid of "mask"
!     call mct_aVect_gather(AVi, AVo, gsmap, 0, mpicom, missing = missing)
     call mct_aVect_gather(AVi, AVo, gsmap, 0, mpicom)
  else
     call mct_aVect_gather(AVi, AVo, gsmap, 0, mpicom)
  endif

  if (masterproc) then
     lsize = size(rglobal)
     call mct_aVect_exportRattr(AVo,"array",rglobal,lsize)
     if (present(missing)) then
        allocate(ivect(lsize))
        call mct_aVect_exportIattr(AVo,"mask",ivect,lsize)
        do n = 1,lsize
           if (ivect(n) == 0) rglobal(n) = missing
        enddo
        deallocate(ivect)
     endif
     call mct_aVect_clean(AVo)
  endif

  call mct_aVect_clean(AVi)

  end subroutine gather_gs_1darray_real

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: gather_1darray_int
!
! !INTERFACE:
  subroutine gather_1darray_int (ilocal, iglobal, clmlevel, imissing)
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
    character(len=*),  intent(in) :: clmlevel    !type of input data
    integer, optional, intent(in) :: imissing    !missing value
!
! !REVISION HISTORY:
! Author: T Craig
!
!EOP
!
! !LOCAL VARIABLES:
    type(mct_gsMap) , pointer :: gsmap       !global seg map
    integer, pointer,dimension(:) :: perm
    character(len=*),parameter :: subname = 'gather_1darray_int'

!-----------------------------------------------------------------------

  call get_clmlevel_gsmap(clmlevel,gsmap,perm)
  if (present(imissing)) then
     call gather_gs_1darray_int(ilocal,iglobal,gsmap,perm,imissing)
  else
     call gather_gs_1darray_int(ilocal,iglobal,gsmap,perm)
  endif

  end subroutine gather_1darray_int

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: gather_1darray_real
!
! !INTERFACE:
  subroutine gather_1darray_real (rlocal, rglobal, clmlevel, missing)
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
    character(len=*)  , intent(in):: clmlevel    !type of input data
    real(r8), optional, intent(in):: missing     !missing value
!
! !REVISION HISTORY:
! Author: T Craig
!
!EOP
!
! !LOCAL VARIABLES:
    type(mct_gsMap),pointer       :: gsmap       !global seg map
    integer, pointer,dimension(:) :: perm
    character(len=*),parameter :: subname = 'gather_1darray_real'

!-----------------------------------------------------------------------

  call get_clmlevel_gsmap(clmlevel,gsmap,perm)
  if (present(missing)) then
     call gather_gs_1darray_real(rlocal,rglobal,gsmap,perm,missing)
  else
     call gather_gs_1darray_real(rlocal,rglobal,gsmap,perm)
  endif

  end subroutine gather_1darray_real

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: gather_2darray_int
!
! !INTERFACE:
  subroutine gather_2darray_int (ilocal, iglobal, clmlevel, imissing)
!
! !DESCRIPTION:
! Wrapper routine to gather integer 2d array
! Assume iglobal only defined on root pe, if not, it's still ok
!
! !USES:
!
! !ARGUMENTS:
    implicit none
    integer, pointer              :: ilocal(:,:)   !local read data (in)
    integer, pointer              :: iglobal(:,:)  !global read data (out)
    character(len=*) , intent(in) :: clmlevel      !type of input data
    integer, optional, intent(in) :: imissing      !missing value
!
! !REVISION HISTORY:
! Author: T Craig
!
!EOP
!
! !LOCAL VARIABLES:
    type(mct_gsMap) , pointer :: gsmap       !global seg map
    integer, pointer,dimension(:) :: perm
    integer, pointer :: locarr(:)
    integer, pointer :: gloarr(:)
    integer :: lbg,ubg,lbl,ubl,sizg,sizl,n
    character(len=*),parameter :: subname = 'gather_2darray_int'

!-----------------------------------------------------------------------

  call get_clmlevel_gsmap(clmlevel,gsmap,perm)
  lbl = lbound(ilocal, dim=2)
  ubl = ubound(ilocal, dim=2)
  sizl = size(ilocal ,dim=1)
  allocate(locarr(sizl))

  if (masterproc) then
     lbg = lbound(iglobal,dim=2)
     ubg = ubound(iglobal,dim=2)
     if (ubg-lbg /= ubl-lbl) then
        write(6,*) trim(subname),' error second dim ',lbg,ubg,lbl,ubl
        call endrun()
     endif
     sizg = size(iglobal,dim=1)
     allocate(gloarr(sizg))
  endif

  do n = lbl,ubl
     locarr(:) = ilocal(:,n)
     if (present(imissing)) then
        call gather_gs_1darray_int(locarr,gloarr,gsmap,perm,imissing)
     else
        call gather_gs_1darray_int(locarr,gloarr,gsmap,perm)
     endif
     if (masterproc) iglobal(:,lbg+n-lbl) = gloarr(:)
  enddo

  deallocate(locarr)
  if (masterproc) deallocate(gloarr)

  end subroutine gather_2darray_int

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: gather_2darray_real
!
! !INTERFACE:
  subroutine gather_2darray_real (rlocal, rglobal, clmlevel, missing)
!
! !DESCRIPTION:
! Wrapper routine to gather real 2d array
! Assume rglobal allocated only on root pe
!
! !USES:
!
! !ARGUMENTS:
    implicit none
    real(r8), pointer             :: rlocal(:,:)   !local read data (in)
    real(r8), pointer             :: rglobal(:,:)  !global read data (out)
    character(len=*)  , intent(in):: clmlevel      !type of input data
    real(r8), optional, intent(in):: missing       !missing value
!
! !REVISION HISTORY:
! Author: T Craig
!
!EOP
!
! !LOCAL VARIABLES:
    type(mct_gsMap),pointer       :: gsmap       !global seg map
    integer, pointer,dimension(:) :: perm
    real(r8), pointer :: locarr(:)
    real(r8), pointer :: gloarr(:)
    integer :: lbg,ubg,lbl,ubl,sizg,sizl,n
    character(len=*),parameter :: subname = 'gather_2darray_real'

!-----------------------------------------------------------------------

  call get_clmlevel_gsmap(clmlevel,gsmap,perm)
  lbl = lbound(rlocal, dim=2)
  ubl = ubound(rlocal, dim=2)
  sizl = size(rlocal ,dim=1)
  allocate(locarr(sizl))

  if (masterproc) then
     lbg = lbound(rglobal,dim=2)
     ubg = ubound(rglobal,dim=2)
     if (ubg-lbg /= ubl-lbl) then
        write(6,*) trim(subname),' error second dim ',lbg,ubg,lbl,ubl
        call endrun()
     endif
     sizg = size(rglobal,dim=1)
     allocate(gloarr(sizg))
  endif

  do n = lbl,ubl
     locarr(:) = rlocal(:,n)
     if (present(missing)) then
        call gather_gs_1darray_real(locarr,gloarr,gsmap,perm,missing)
     else
        call gather_gs_1darray_real(locarr,gloarr,gsmap,perm)
     endif
     if (masterproc) rglobal(:,lbg+n-lbl) = gloarr(:)
  enddo

  deallocate(locarr)
  if (masterproc) deallocate(gloarr)

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
    character(len=*),parameter :: subname = 'allgather_1darray_int'

!-----------------------------------------------------------------------

    call gather_data_to_master(ilocal,iglobal,clmlevel)
    call mpi_bcast (iglobal, size(iglobal), MPI_INTEGER, 0, mpicom, ier)
    if (ier/=0 ) then
       write(6,*) trim(subname),ier
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
    integer :: ier                    !errorcode
    character(len=*),parameter :: subname = 'allgather_1darray_real'

!-----------------------------------------------------------------------

    call gather_data_to_master(rlocal,rglobal,clmlevel)
    call mpi_bcast (rglobal, size(rglobal), MPI_REAL8, 0, mpicom, ier)
    if (ier/=0 ) then
       write(6,*) trim(subname),ier
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
    integer :: ier                    !errorcode
    character(len=*),parameter :: subname = 'allgather_2darray_int'

!-----------------------------------------------------------------------

    call gather_data_to_master(ilocal,iglobal,clmlevel)
    call mpi_bcast (iglobal, size(iglobal), MPI_INTEGER, 0, mpicom, ier)
    if (ier/=0 ) then
       write(6,*) trim(subname),ier
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
    integer :: ier                    !errorcode
    character(len=*),parameter :: subname = 'allgather_2darray_real'

!-----------------------------------------------------------------------

    call gather_data_to_master(rlocal,rglobal,clmlevel)
    call mpi_bcast (rglobal, size(rglobal), MPI_REAL8, 0, mpicom, ier)
    if (ier/=0 ) then
       write(6,*) trim(subname),ier
       call endrun
    endif

  end subroutine allgather_2darray_real

end module spmdGathScatMod
