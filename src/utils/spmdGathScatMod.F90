module spmdGathScatMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Perform SPMD gather and scatter operations.
  !
  ! !USES:
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use spmdMod        , only : masterproc, mpicom
  use mct_mod        , only : mct_aVect, mct_gsMap
  use mct_mod        , only : mct_aVect_init, mct_aVect_importIattr, mct_aVect_scatter
  use mct_mod        , only : mct_aVect_gather, mct_aVect_exportIattr, mct_aVect_clean
  use mct_mod        , only : mct_aVect_exportRattr, mct_aVect_importRattr
  use abortutils     , only : endrun
  use clm_varctl     , only : iulog
  !
  ! !PUBLIC TYPES:
  implicit none
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: scatter_data_from_master
  public :: gather_data_to_master
  public :: gsmap_global_init
  
  type(mct_gsMap), target, public  :: gsmap_global  ! global seg map

!-----------------------------------------------------------------------
contains
!-----------------------------------------------------------------------

  subroutine gsmap_global_init(gindex_global)
    !
    ! !USES:
    use spmdMod   , only : mpicom, comp_id
    use decompMod , only : nglob_x, nglob_y
    use mct_mod   , only : mct_gsMap_init
    !
    ! !ARGUMENTS:
    integer, intent(in) :: gindex_global(:)
    !
    ! !LOCAL VARIABLES:
    integer :: lsize, gsize
    !-----------------------------------------------------------------------

    lsize = size(gindex_global)
    gsize = nglob_x * nglob_y

    call mct_gsMap_init(gsmap_global, gindex_global, mpicom, comp_id, lsize, gsize)

  end subroutine gsmap_global_init

  !-----------------------------------------------------------------------
  subroutine scatter_data_from_master (aglobal, alocal)
    !
    ! !DESCRIPTION:
    ! Wrapper routine to scatter int 1d array
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    integer , pointer :: aglobal(:) ! global data (input)
    integer , pointer :: alocal(:)  ! local  data (output)

    ! !LOCAL VARIABLES:
    integer                    :: n,lb,ub ! indices
    integer                    :: lsize      ! size of local array
    type(mct_aVect)            :: AVi, AVo   ! attribute vectors
    integer ,pointer           :: adata(:)   ! local data array
    character(len=*),parameter :: subname = 'scatter_1darray_int'
    !-----------------------------------------------------------------------

    if (masterproc) then
       call mct_aVect_init(AVi, rList="", iList="f1", lsize=lsize)
       call mct_aVect_importIattr(AVi, 'f1', aglobal, size(aglobal,dim=1))
    endif
    call mct_aVect_scatter(AVi, AVo, gsmap_global, 0, mpicom)
    lsize = size(alocal, dim=1)
    allocate(adata(lsize))
    call mct_aVect_exportIattr(AVo, 'f1', adata, lsize)
    lb = lbound(alocal, dim=1); ub = ubound(alocal, dim=1)
    do n = lb,ub
       alocal(n) = adata(n-lb+1)
    enddo
    deallocate(adata)
    if (masterproc) then
       call mct_aVect_clean(AVi)
    endif
    call mct_aVect_clean(AVo)

  end subroutine scatter_data_from_master

  !-----------------------------------------------------------------------
  subroutine gather_data_to_master (alocal, aglobal)
    !
    ! !DESCRIPTION:
    ! Wrapper routine to gather int 1d array
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    integer , pointer :: alocal(:)  ! local  data (output)
    integer , pointer :: aglobal(:) ! global data (input)
    !
    ! !LOCAL VARIABLES:
    integer            :: n,lb,ub    ! indices
    integer            :: lsize      ! size of local array
    type(mct_aVect)    :: AVi, AVo   ! attribute vectors
    integer ,pointer   :: adata(:)   ! temporary data array
    character(len=*),parameter :: subname = 'gather_1darray_int'
    !-----------------------------------------------------------------------

    lsize = size(alocal, dim=1)
    lb = lbound(alocal, dim=1); ub = ubound(alocal, dim=1)
    call mct_aVect_init(AVi, rList="", iList='f1', lsize=lsize)
    allocate(adata(lsize))
    do n = lb,ub
       adata(n-lb+1) = alocal(n)
    enddo
    call mct_aVect_importIattr(AVi, 'f1', adata, lsize)
    deallocate(adata)
    call mct_aVect_gather(AVi, AVo, gsmap_global, 0, mpicom)
    if (masterproc) then
       lsize = size(aglobal,dim=1)
       call mct_aVect_exportIattr(AVo, 'f1', aglobal, lsize)
       call mct_aVect_clean(AVo)
    endif
    call mct_aVect_clean(AVi)

  end subroutine gather_data_to_master

end module spmdGathScatMod
