module FilterDataBaseType
  implicit none
  private
  public :: filtered_data_base_type
type , abstract :: filtered_data_base_type
  private
  integer, private :: num_filter
  Integer, pointer, private :: filter(:)
contains
   procedure , public :: InitDataBase
   procedure(InitData_interface) , public, deferred :: InitData
   procedure , public :: ReInitDataBase
   procedure(ReInitData_interface) , public, deferred :: ReInitData
   procedure , public :: ValidateFiltersMatch
   procedure(CompressData_interface) , public, deferred :: CompressData
   procedure(DecompressData_interface) , public, deferred :: DecompressData
end type filtered_data_base_type

  abstract interface

  subroutine InitData_interface(this)
    import :: filtered_data_base_type
    class(filtered_data_base_type), intent(inout) :: this
  end subroutine InitData_interface

  subroutine ReInitData_interface(this)
    import :: filtered_data_base_type
    class(filtered_data_base_type), intent(inout) :: this
  end subroutine ReInitData_interface

  subroutine CompressData_interface(this)
    import :: filtered_data_base_type
    class(filtered_data_base_type), intent(inout) :: this
  end subroutine CompressData_interface

  subroutine DecompressData_interface(this)
    import :: filtered_data_base_type
    class(filtered_data_base_type), intent(inout) :: this
  end subroutine DecompressData_interface

  end interface

contains

  subroutine ValidateFiltersMatch(this)
    class(filtered_data_base_type), intent(in) :: this
  end subroutine ValidateFiltersMatch

  subroutine InitDataBase(this)
    class(filtered_data_base_type), intent(in) :: this
  end subroutine InitDataBase

  subroutine ReInitDataBase(this)
    class(filtered_data_base_type), intent(in) :: this
  end subroutine ReInitDataBase

end module FilterDataBaseType
