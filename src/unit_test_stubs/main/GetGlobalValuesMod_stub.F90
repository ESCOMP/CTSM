module GetGlobalValuesMod

  ! Stub of GetGlobalValuesMod, which satisfies routine signatures with minimal
  ! dependencies

  implicit none
  
  public :: GetGlobalWrite
  public :: GetGlobalIndex

contains

  subroutine GetGlobalWrite(decomp_index, clmlevel)
    integer, intent(in) :: decomp_index
    character(len=*), intent(in) :: clmlevel

    ! do nothing
  end subroutine GetGlobalWrite

  !-----------------------------------------------------------------------
  integer function GetGlobalIndex(decomp_index, clmlevel)

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
    GetGlobalIndex = -1
  end function GetGlobalIndex

end module GetGlobalValuesMod
