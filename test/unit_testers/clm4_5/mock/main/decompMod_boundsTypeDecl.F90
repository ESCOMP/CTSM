module decompMod

  ! This is a mock replacement for decompMod, which only contains the type declaration
  ! for the bounds_type. The purpose of this is to remove a bunch of dependencies of the
  ! true decompMod

  !
  ! !PUBLIC TYPES:
  implicit none

  type bounds_type
     integer :: begg, endg       ! beginning and ending gridcell index
     integer :: begl, endl       ! beginning and ending landunit index
     integer :: begc, endc       ! beginning and ending column index
     integer :: begp, endp       ! beginning and ending pft index

     integer :: level            ! whether defined on the proc or clump level
     integer :: clump_index      ! if defined on the clump level, this gives the clump index
  end type bounds_type
  public bounds_type


end module decompMod
