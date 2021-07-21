module GetGlobalValuesMod

  ! Stub of GetGlobalValuesMod, which satisfies routine signatures with minimal
  ! dependencies

  implicit none
  
  public :: write_point_context
  public :: get_global_index

contains

  subroutine write_point_context(decomp_index, clmlevel)
    integer, intent(in) :: decomp_index
    character(len=*), intent(in) :: clmlevel

    ! do nothing
  end subroutine write_point_context

  !-----------------------------------------------------------------------
  integer function get_global_index(decomp_index, clmlevel)

    !----------------------------------------------------------------
    ! Description
    ! Determine global index space value for target point at given clmlevel
    !
    ! Uses:
    !
    ! Arguments 
    integer          , intent(in) :: decomp_index
    character(len=*) , intent(in) :: clmlevel

    ! De essentially nothing, just set the index to a negative value to signal it's not real
    get_global_index = -1
  end function get_global_index

end module GetGlobalValuesMod
