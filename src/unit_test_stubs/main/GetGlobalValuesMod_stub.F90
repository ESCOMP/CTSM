module GetGlobalValuesMod

  ! Stub of GetGlobalValuesMod, which satisfies routine signatures with minimal
  ! dependencies

  implicit none
  
  public :: write_point_context
  public :: get_global_index

contains

  subroutine write_point_context(subgrid_index, subgrid_level)
    integer, intent(in) :: subgrid_index
    character(len=*), intent(in) :: subgrid_level

    ! do nothing
  end subroutine write_point_context

  !-----------------------------------------------------------------------
  integer function get_global_index(subgrid_index, subgrid_level)

    !----------------------------------------------------------------
    ! Description
    ! Determine global index space value for target point at given subgrid_level
    !
    ! Uses:
    !
    ! Arguments 
    integer          , intent(in) :: subgrid_index
    character(len=*) , intent(in) :: subgrid_level

    ! De essentially nothing, just set the index to a negative value to signal it's not real
    get_global_index = -1
  end function get_global_index

end module GetGlobalValuesMod
