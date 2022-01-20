module mkesmfMod

  use ESMF
  use shr_kind_mod , only : r8 => shr_kind_r8, r4 => shr_kind_r4
  use shr_sys_mod  , only : shr_sys_abort
  use mkUtilsMod   , only : chkerr

  implicit none
  private

  public :: regrid_rawdata

  character(len=*) , parameter :: u_FILE_u = &
       __FILE__

!===============================================================
contains
!===============================================================

  subroutine regrid_rawdata(field_i, field_o,  routehandle, data_i, data_o, rc)

    ! input/output variables
    type(ESMF_Field)       , intent(in)    :: field_i
    type(ESMF_Field)       , intent(inout) :: field_o
    type(ESMF_RouteHandle) , intent(inout) :: routehandle
    real(r8)               , intent(in)    :: data_i(:)
    real(r8)               , intent(inout) :: data_o(:)
    integer                , intent(out)   :: rc

    ! local variables
    logical           :: checkflag = .false.
    real(r8), pointer :: dataptr(:)
    ! --------------------------------------------

    rc = ESMF_SUCCESS

    ! Interpolate data_i to data_o
    call ESMF_FieldGet(field_i, farrayptr=dataptr, rc=rc)
    dataptr(:) = data_i(:)
    call ESMF_FieldGet(field_o, farrayptr=dataptr, rc=rc)
    dataptr(:) = 0._r8
    call ESMF_FieldRegrid(field_i, field_o, routehandle=routehandle, &
         termorderflag=ESMF_TERMORDER_SRCSEQ, checkflag=checkflag, zeroregion=ESMF_REGION_TOTAL, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field_o, farrayptr=dataptr, rc=rc)
    data_o(:) = dataptr(:)
    call ESMF_VMLogMemInfo("After field regrid in regrid_data")

  end subroutine regrid_rawdata

end module mkesmfMod
