module column_varcon

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module containing column indices and associated variables and routines.
  !
  ! !USES:
#include "shr_assert.h"
  use landunit_varcon, only : isturb_MIN
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  private

  !------------------------------------------------------------------
  ! Initialize column type constants
  !------------------------------------------------------------------

  ! urban column types

  integer, parameter, public :: icol_roof        = isturb_MIN*10 + 1
  integer, parameter, public :: icol_sunwall     = isturb_MIN*10 + 2
  integer, parameter, public :: icol_shadewall   = isturb_MIN*10 + 3
  integer, parameter, public :: icol_road_imperv = isturb_MIN*10 + 4
  integer, parameter, public :: icol_road_perv   = isturb_MIN*10 + 5

  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: is_hydrologically_active   ! returns true if the given column type is hydrologically active
  public :: ice_class_to_col_itype     ! convert an ice class (1..maxpatch_glc) into col%itype
  public :: col_itype_to_ice_class     ! convert col%itype into an ice class (1..maxpatch_glc)
  public :: write_coltype_metadata     ! write column type metadata to a netcdf file

  character(len=*), parameter, private :: sourcefile = &
       __FILE__

contains

  !-----------------------------------------------------------------------
  function is_hydrologically_active(col_itype, lun_itype) &
       result(hydrologically_active)
    !
    ! !DESCRIPTION:
    ! Returns a logical value saying whether the given column type is hydrologically
    ! active
    !
    ! Note that calling this can be bad for performance, because it operates on a single
    ! point rather than a loop. So in performance-critical parts of the code (or just
    ! about anywhere, really), you should use the pre-set col%hydrologically_active(c).
    !
    ! !USES:
    use landunit_varcon, only : istsoil, istcrop
    !
    ! !ARGUMENTS:
    logical :: hydrologically_active  ! function result
    integer, intent(in) :: col_itype  ! col%itype value
    integer, intent(in) :: lun_itype  ! lun%itype value for the landunit on which this column sits
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'is_hydrologically_active'
    !-----------------------------------------------------------------------

    ! If we had an easy way to figure out which landunit a column was on based on
    ! col_itype (which would be very helpful!), then we wouldn't need lun_itype.

    if (lun_itype == istsoil .or. lun_itype == istcrop) then
       hydrologically_active = .true.
    else if (col_itype == icol_road_perv) then
       hydrologically_active = .true.
    else
       hydrologically_active = .false.
    end if

  end function is_hydrologically_active

  
  !-----------------------------------------------------------------------
  function ice_class_to_col_itype(ice_class) result(col_itype)
    !
    ! !DESCRIPTION:
    ! Convert an ice class (1..maxpatch_glc) into col%itype
    !
    ! !USES:
    use clm_varpar, only : maxpatch_glc
    use landunit_varcon, only : istice
    !
    ! !ARGUMENTS:
    integer :: col_itype             ! function result
    integer, intent(in) :: ice_class ! ice class, between 1 and maxpatch_glc
    !
    ! !LOCAL VARIABLES:
    
    character(len=*), parameter :: subname = 'ice_class_to_col_itype'
    !-----------------------------------------------------------------------
    
    SHR_ASSERT_FL((1 <= ice_class .and. ice_class <= maxpatch_glc), sourcefile, __LINE__)

    col_itype = istice*100 + ice_class

  end function ice_class_to_col_itype

  !-----------------------------------------------------------------------
  function col_itype_to_ice_class(col_itype) result(ice_class)
    !
    ! !DESCRIPTION:
    ! Convert a col%itype value (for an ice landunit) into an ice class (1..maxpatch_glc)
    !
    ! !USES:
    use clm_varpar, only : maxpatch_glc
    use landunit_varcon, only : istice
    !
    ! !ARGUMENTS:
    integer :: ice_class             ! function result
    integer, intent(in) :: col_itype ! col%itype value for an ice landunit
    !
    ! !LOCAL VARIABLES:
    
    character(len=*), parameter :: subname = 'col_itype_to_ice_class'
    !-----------------------------------------------------------------------
    
    ice_class = col_itype - istice*100

    ! The following assertion is here to ensure that col_itype is really from an
    ! istice landunit
    SHR_ASSERT_FL((1 <= ice_class .and. ice_class <= maxpatch_glc), sourcefile, __LINE__)

  end function col_itype_to_ice_class

  !-----------------------------------------------------------------------
  subroutine write_coltype_metadata(att_prefix, ncid)
    !
    ! !DESCRIPTION:
    ! Writes column type metadata to a netcdf file.
    !
    ! Note that, unlike pft and landunit metadata, this column type metadata is NOT
    ! stored in an array. This is because of the trickiness of encoding column values for
    ! crop & ice. So instead, other code must call this routine to do the work of
    ! adding the appropriate metadata directly to a netcdf file.
    !
    ! !USES:
    use ncdio_pio, only : file_desc_t, ncd_global, ncd_putatt
    !
    ! !ARGUMENTS:
    character(len=*)  , intent(in)    :: att_prefix ! prefix for attributes (e.g., 'icol_')
    type(file_desc_t) , intent(inout) :: ncid       ! local file id
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'write_coltype_metadata'
    !-----------------------------------------------------------------------

    call ncd_putatt(ncid, ncd_global, att_prefix // 'vegetated_or_bare_soil', 1)
    call ncd_putatt(ncid, ncd_global, att_prefix // 'crop'                  , 2) 
    call ncd_putatt(ncid, ncd_global, att_prefix // 'crop_noncompete'       , '2*100+m, m=cft_lb,cft_ub')
    call ncd_putatt(ncid, ncd_global, att_prefix // 'landice'               , 3) 
    call ncd_putatt(ncid, ncd_global, att_prefix // 'landice', '4*100+m, m=1,glcnec')
    call ncd_putatt(ncid, ncd_global, att_prefix // 'deep_lake'             , 5) 
    call ncd_putatt(ncid, ncd_global, att_prefix // 'wetland'               , 6) 
    call ncd_putatt(ncid, ncd_global, att_prefix // 'urban_roof'            , icol_roof)
    call ncd_putatt(ncid, ncd_global, att_prefix // 'urban_sunwall'         , icol_sunwall)
    call ncd_putatt(ncid, ncd_global, att_prefix // 'urban_shadewall'       , icol_shadewall)
    call ncd_putatt(ncid, ncd_global, att_prefix // 'urban_impervious_road' , icol_road_imperv)
    call ncd_putatt(ncid, ncd_global, att_prefix // 'urban_pervious_road'   , icol_road_perv)

  end subroutine write_coltype_metadata


end module column_varcon
