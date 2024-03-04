module mkdomainMod

  !-------------------------------
  ! Determine lon/lat of model 
  !-------------------------------

  use ESMF
  use pio
  use shr_kind_mod, only : r8 => shr_kind_r8
  use shr_sys_mod , only : shr_sys_abort 
  use mkutilsMod  , only : chkerr
  use mkvarctl    , only : root_task, ndiag
  use mkfileMod   , only : mkfile_output

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
    type(ESMF_Mesh)   , intent(in)    :: mesh_o
    real(r8)          , intent(out)   :: lon_o(:)
    real(r8)          , intent(out)   :: lat_o(:)
    integer           , intent(out)   :: rc

    ! local variables:
    integer               :: no
    integer               :: ns_o
    integer               :: spatialDim
    integer               :: k
    real(r8), allocatable :: ownedElemCoords(:)
    character(len=*), parameter :: subname = 'mkdomain'
    !-----------------------------------------------------------------------

    rc = ESMF_SUCCESS

    if (root_task) then
       write(ndiag,*)
       write(ndiag,'(1x,80a1)') ('=',k=1,80)
       write(ndiag,*)
       write(ndiag,'(a)') 'Attempting to create model lats and lons from model mesh .....'
       flush(ndiag)
    end if

    call ESMF_MeshGet(mesh_o, spatialDim=spatialDim, numOwnedElements=ns_o, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    allocate(ownedElemCoords(spatialDim*ns_o))
    call ESMF_MeshGet(mesh_o, ownedElemCoords=ownedElemCoords, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    do no = 1,ns_o
       lon_o(no) = ownedElemCoords(2*no-1)
       lat_o(no) = ownedElemCoords(2*no)
    end do

    if (root_task) then
       write (ndiag,'(a)') 'Successfully made model lats and lons'
       flush(ndiag)
    end if

  end subroutine mkdomain

end module mkdomainMod
