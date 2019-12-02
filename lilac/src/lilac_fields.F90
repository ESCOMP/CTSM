module lilac_fields

  ! This module is used by both CTSM and the lilac atmcap to ensure that the field bundles
  ! exchanged between components are identical

  use ESMF
  use lilac_methods, only : chkerr
  use lilac_utils  , only : atm2lnd, lnd2atm

  implicit none
  public

  character(*),parameter :: u_FILE_u = &
       __FILE__

!===============================================================================
contains
!===============================================================================

  subroutine lilac_field_bundle_to_land(mesh, fieldbundle, rc)

    type(ESMF_Mesh)        :: mesh
    type(ESMF_FieldBundle) :: fieldbundle
    integer,  intent(out)  :: rc

    integer :: n

    rc = ESMF_SUCCESS

    ! Add empty fields to field bundle
    do n = 1, size(atm2lnd)
       call fldbundle_add(trim(atm2lnd(n)%fldname), mesh, fieldbundle, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end do

  end subroutine lilac_field_bundle_to_land

  !===============================================================================

  subroutine lilac_field_bundle_from_land(mesh, fieldbundle, rc)

    type(ESMF_Mesh)        :: mesh 
    type(ESMF_FieldBundle) :: fieldbundle
    integer, intent(out)   :: rc

    integer :: n

    rc = ESMF_SUCCESS

    ! Add empty fields to field bundle
    do n = 1, size(atm2lnd)
       call fldbundle_add( trim(lnd2atm(n)%fldname), mesh, fieldbundle, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end do

  end subroutine lilac_field_bundle_from_land

  !===============================================================================

  subroutine fldbundle_add(fldname, mesh, fieldbundle, rc)

    !---------------------------
    ! Create an empty input field with name 'stdname' to add to fieldbundle
    !---------------------------

    ! input/output variables
    character(len=*)       , intent(in)    :: fldname
    type(ESMF_Mesh)        , intent(in)    :: mesh
    type(ESMF_FieldBundle) , intent(inout) :: fieldbundle
    integer                , intent(out)   :: rc

    ! local variables
    type(ESMF_Field) :: field
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    field = ESMF_FieldCreate(mesh, ESMF_TYPEKIND_R8 , meshloc=ESMF_MESHLOC_ELEMENT , name=trim(fldname), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_FieldBundleAdd(fieldbundle, (/field/), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

  end subroutine fldbundle_add


end module lilac_fields
