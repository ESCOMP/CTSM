module mkdomainMod

  !-------------------------------
  ! Determine lon/lat of model 
  !-------------------------------

  use ESMF
  use shr_kind_mod, only : r8 => shr_kind_r8
  use mkutilsMod  , only : chkerr

  implicit none
  private

  public :: mkdomain

  character(len=*) , parameter :: u_FILE_u = &
       __FILE__

!=================================================================================
contains
!=================================================================================

  subroutine mkdomain(mesh_o, lon_o, lat_o, rc)
    
    ! input output variables
    type(ESMF_Mesh) , intent(in)  :: mesh_o
    real(r8)        , intent(out) :: lon_o(:)
    real(r8)        , intent(out) :: lat_o(:)
    integer         , intent(out) :: rc

    ! local variables:
    integer               :: no
    integer               :: ns_o
    integer               :: spatialDim
    real(r8), allocatable :: ownedElemCoords(:)
    character(len=*), parameter :: subname = 'mkdomain'
    !-----------------------------------------------------------------------

    rc = ESMF_SUCCESS

    call ESMF_MeshGet(mesh_o, spatialDim=spatialDim, numOwnedElements=ns_o, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    allocate(ownedElemCoords(spatialDim*ns_o))
    call ESMF_MeshGet(mesh_o, ownedElemCoords=ownedElemCoords, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    do no = 1,ns_o
       lon_o(no) = ownedElemCoords(2*no-1)
       lat_o(no) = ownedElemCoords(2*no)
    end do

  end subroutine mkdomain

end module mkdomainMod
