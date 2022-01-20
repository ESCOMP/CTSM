module mkesmfMod

  use ESMF
  use pio          , only : file_desc_t, iosystem_desc_t, io_desc_t, var_desc_t
  use pio          , only : pio_openfile, pio_closefile, pio_nowrite
  use pio          , only : pio_double, pio_real, pio_int, pio_offset_kind, pio_get_var
  use pio          , only : pio_read_darray, pio_setframe, pio_fill_double, pio_get_att
  use pio          , only : io_desc_t, var_desc_t 
  use pio          , only : PIO_BCAST_ERROR, PIO_RETURN_ERROR, PIO_NOERR, PIO_INTERNAL_ERROR
  use pio          , only : PIO_REAL, PIO_INT, PIO_DOUBLE, PIO_SHORT
  use shr_kind_mod , only : r8 => shr_kind_r8, r4 => shr_kind_r4
  use shr_sys_mod  , only : shr_sys_abort
  use mkUtilsMod   , only : chkerr

  implicit none
  private

  public :: regrid_data

  character(len=*) , parameter :: u_FILE_u = &
       __FILE__

!===============================================================
contains
!===============================================================

  subroutine regrid_data(field_i, field_o, varname, filename, data_i, data_o, rc)

    use mkvarctl , only : ndiag, root_task
    use mkpioMod , only : mkpio_iodesc_rawdata, pio_iotype, pio_ioformat, io_subsystem

    ! input/output variables
    type(ESMF_Field) , intent(inout) :: field_i   ! raw data field
    type(ESMF_Field) , intent(inout) :: field_o   ! model field  
    character(len=*) , intent(in)    :: varname   ! field name in rawdata file
    character(len=*) , intent(in)    :: filename  ! file name of rawdata file
    real(r8)         , intent(inout) :: data_i(:) ! input raw data
    real(r8)         , intent(inout) :: data_o(:) ! mapped data
    integer          , intent(out)   :: rc

    ! local variables
    type(ESMF_Mesh)        :: mesh_i
    type(ESMF_RouteHandle) :: routehandle
    integer                :: srcMaskValue = 0
    integer                :: dstMaskValue = -987987 ! spval for RH mask values
    integer                :: srcTermProcessing_Value = 0
    type(file_desc_t)      :: pioid
    type(var_desc_t)       :: pio_varid
    type(io_desc_t)        :: pio_iodesc 
    integer                :: pio_vartype
    real(r4), allocatable  :: data_real(:)
    real(r8), allocatable  :: data_double(:)
    real(r8), pointer      :: dataptr(:)
    integer                :: lsize 
    integer                :: rcode
    logical                :: checkflag = .false.
    character(len=*), parameter :: subname = 'mklakwat'
    !-------------------------------------------------

    rc = ESMF_SUCCESS

    ! create route handle to map field_model to field_data
    call ESMF_FieldRegridStore(field_i, field_o, routehandle=routehandle, &
         srcMaskValues=(/srcMaskValue/), dstMaskValues=(/dstMaskValue/), &
         regridmethod=ESMF_REGRIDMETHOD_CONSERVE, normType=ESMF_NORMTYPE_DSTAREA, &
         srcTermProcessing=srcTermProcessing_Value, &
         ignoreDegenerate=.true., unmappedaction=ESMF_UNMAPPEDACTION_IGNORE, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! Get data_i - Read in varname from filename 
    rcode = pio_openfile(io_subsystem, pioid, pio_iotype, trim(filename), pio_nowrite)
    call ESMF_FieldGet(field_i, mesh=mesh_i, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call mkpio_iodesc_rawdata(mesh_i, trim(varname), pioid, pio_varid, pio_vartype, pio_iodesc, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (pio_vartype == PIO_REAL) then
       allocate(data_real(lsize))
       call pio_read_darray(pioid, pio_varid, pio_iodesc, data_real, rcode)
       data_i(:) = real(data_real(:), kind=r8)
       deallocate(data_real)
    else if (pio_vartype == PIO_DOUBLE) then
       allocate(data_double(lsize))
       call pio_read_darray(pioid, pio_varid, pio_iodesc, data_double, rcode)
       data_i(:) = data_double(:)
       deallocate(data_double)
    else
       call shr_sys_abort(subName//"ERROR: only real and double types are supported")
    end if

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

  end subroutine regrid_data

end module mkesmfMod
