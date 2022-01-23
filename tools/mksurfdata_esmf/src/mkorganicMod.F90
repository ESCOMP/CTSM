module mkorganicMod

  ! make organic matter dataset

  use ESMF
  use pio
  use shr_kind_mod , only : r8 => shr_kind_r8, r4 => shr_kind_r4
  use shr_sys_mod  , only : shr_sys_abort
  use mkpioMod     , only : mkpio_get_rawdata, mkpio_get_dimlengths
  use mkpioMod     , only : mkpio_iodesc_rawdata, pio_iotype, pio_ioformat, pio_iosystem
  use mkesmfMod    , only : regrid_rawdata
  use mkutilsMod   , only : chkerr
  use mkvarctl     , only : root_task, ndiag
  use mkvarpar     , only : nlevsoi

  implicit none
  private

  public mkorganic      ! Set organic soil

  character(len=*) , parameter :: u_FILE_u = &
       __FILE__

!=================================================================================
contains
!=================================================================================

  subroutine mkorganic(file_mesh_i, file_data_i, mesh_o, organic_o, rc)

    ! input/output variables
    character(len=*)  , intent(in)    :: file_mesh_i      ! input mesh file name
    character(len=*)  , intent(in)    :: file_data_i      ! input data file name
    type(ESMF_Mesh)   , intent(in)    :: mesh_o
    real(r8)          , intent(out)   :: organic_o(:,:)   ! output grid: %lake
    integer           , intent(out)   :: rc

    ! local variables
    type(ESMF_RouteHandle) :: routehandle
    type(ESMF_Mesh)        :: mesh_i
    type(ESMF_Field)       :: field_i
    type(ESMF_Field)       :: field_o
    type(file_desc_t)      :: pioid
    type(var_desc_t)       :: pio_varid
    type(io_desc_t)        :: pio_iodesc
    integer                :: pio_vartype
    integer                :: ni,no
    integer                :: ns_i, ns_o
    integer                :: nlay
    integer                :: n, l  ! indices
    integer                :: rcode, ier             ! error status
    integer                :: srcMaskValue = -987987 ! spval for RH mask values
    integer                :: dstMaskValue = -987987 ! spval for RH mask values
    integer                :: srcTermProcessing_Value = 0
    integer                :: ndims
    integer , allocatable  :: dimlengths(:)
    real(r4), allocatable  :: data_real(:,:)
    real(r8), allocatable  :: data_double(:,:)
    real(r8), allocatable  :: data_i(:,:)
    real(r8), allocatable  :: data_o(:,:)
    character(len=256)     :: varname
    character(len=*), parameter :: subname = 'mkorganic'
    !-----------------------------------------------------------------------

    rc = ESMF_SUCCESS

    if (root_task) then
       write (ndiag,'(a)') ' Attempting to make organic mater dataset .....'
    end if

    ! Set variable name in raw data file
    varname = 'ORGANIC'

    ! Open raw data file - need to do this first to obtain ungridded dimension size
    call ESMF_VMLogMemInfo("Before pio_openfile for "//trim(file_data_i))
    rcode = pio_openfile(pio_iosystem, pioid, pio_iotype, trim(file_data_i), pio_nowrite)
    call ESMF_VMLogMemInfo("After pio_openfile "//trim(file_data_i))

    ! Get dimensions of raw data. 
    ! NOTE:
    !  - raw data is dimensions (lon,lat,lev) 
    !  - input read from pio has dimensions(n,lev)
    !  - esmf field dataptr has dimensions (lev,n)
    allocate(dimlengths(3))
    call mkpio_get_dimlengths(pioid, trim(varname), ndims, dimlengths)
    nlay = dimlengths(3)
    if (nlay /= nlevsoi) then
       write(6,*)'nlay, nlevsoi= ',nlay,nlevsoi,' do not match'
       call shr_sys_abort()
    end if
    deallocate(dimlengths)

    ! Read in input mesh
    call ESMF_VMLogMemInfo("Before create mesh_i in "//trim(subname))
    mesh_i = ESMF_MeshCreate(filename=trim(file_mesh_i), fileformat=ESMF_FILEFORMAT_ESMFMESH, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Create io descriptor for input raw data 
    ! This will query the raw data file for the dimensions of the variable varname and 
    ! create iodesc for either single or multi level input data
    call ESMF_VMLogMemInfo("Before mkpio_iodesc in "//trim(subname))
    call mkpio_iodesc_rawdata(mesh_i, trim(varname), pioid, pio_varid, pio_vartype, pio_iodesc, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMLogMemInfo("After mkpio_iodesc in "//trim(subname))

    ! Determine ns_i and allocate data_i
    call ESMF_MeshGet(mesh_i, numOwnedElements=ns_i, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMLogMemInfo("After create mesh_i in "//trim(subname))
    allocate(data_i(nlay,ns_i),stat=ier)
    if (ier/=0) call shr_sys_abort()

    ! Determine ns_o and allocate data_o
    ns_o = size(organic_o, dim=1)
    allocate(data_o(nlay,ns_o),stat=ier)
    if (ier/=0) call shr_sys_abort()

    ! Read the input raw data (all levels read at once)
    ! - levels are the innermost dimension for esmf fields
    ! - levels are the outermost dimension in pio reads
    ! Input data is read into (ns_i,nlay) array and then transferred to data_i(nlay,ns_i)
    if (pio_vartype == PIO_REAL) then
       write(6,*)'DEBUG: ns_i, nlay = ',ns_i, nlay
       allocate(data_real(ns_i,nlay))
       call pio_read_darray(pioid, pio_varid, pio_iodesc, data_real, rcode)
       do l = 1,nlay
          do n = 1,ns_i
             data_i(l,n) = real(data_real(n,l), kind=r8)
          end do
       end do
       deallocate(data_real)
    else if (pio_vartype == PIO_DOUBLE) then
       allocate(data_double(ns_i,nlay))
       call pio_read_darray(pioid, pio_varid, pio_iodesc, data_double, rcode)
       do l = 1,nlay
          do n = 1,ns_i
             data_i(l,n) = data_double(n,l)
          end do
       end do
       deallocate(data_double)
    else
       call shr_sys_abort(subName//"ERROR: only real and double types are supported")
    end if
    call ESMF_VMLogMemInfo("After read_darray in "//trim(subname))

    ! Create field on input mesh - using dimension information from raw data file
    field_i = ESMF_FieldCreate(mesh_i, ESMF_TYPEKIND_R8, meshloc=ESMF_MESHLOC_ELEMENT, &
               ungriddedLbound=(/1/), ungriddedUbound=(/nlay/), gridToFieldMap=(/2/), rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMLogMemInfo("After field_i creation in "//trim(subname))

    ! Create field on model model
    field_o = ESMF_FieldCreate(mesh_o, ESMF_TYPEKIND_R8, meshloc=ESMF_MESHLOC_ELEMENT, &
               ungriddedLbound=(/1/), ungriddedUbound=(/nlay/), gridToFieldMap=(/2/), rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMLogMemInfo("After field_o creation in "//trim(subname))

    ! Create route handle to map field_model to field_data
    call ESMF_FieldRegridStore(field_i, field_o, routehandle=routehandle, &
        !srcMaskValues=(/srcMaskValue/), dstMaskValues=(/dstMaskValue/), &
         regridmethod=ESMF_REGRIDMETHOD_CONSERVE, normType=ESMF_NORMTYPE_FRACAREA, &
         srcTermProcessing=srcTermProcessing_Value, &
         ignoreDegenerate=.true., unmappedaction=ESMF_UNMAPPEDACTION_IGNORE, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMLogMemInfo("After regridstore in "//trim(subname))

    ! Regrid raw data to model resolution
    call regrid_rawdata(field_i, field_o, routehandle, data_i, data_o, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMLogMemInfo("After regrid_data in "//trim(subname))

    ! Close the file and free the iodesc and deallocate memory
    call pio_closefile(pioid)
    call pio_freedecomp(pio_iosystem, pio_iodesc)
    call ESMF_VMLogMemInfo("After pio_freedecomp in "//trim(subname))

    ! Now compute organic_o - need to remember that nlay is innermost dimension in data_o
    ! but outermost dimension in organic_o
    do l = 1,nlevsoi
       do n = 1,size(organic_o, dim=1)
          organic_o(n,l) = data_o(l,n)
          if ((organic_o(n,l)) > 130.000001_r8) then
             write (6,*) trim(subname)//' error: organic = ',organic_o(n,l), &
                  ' greater than 130.000001 for n,lev ',n,l
             call shr_sys_abort()
          end if
       enddo
    end do

    ! Release memory
    ! TODO: something is hanging here
    ! call ESMF_RouteHandleDestroy(routehandle, nogarbage = .true., rc=rc)
    ! if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort()
    ! call ESMF_MeshDestroy(mesh_i, nogarbage = .true., rc=rc)
    ! if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort()
    ! call ESMF_VMLogMemInfo("After destroy operations in "//trim(subname))

    if (root_task) then
       write (ndiag,'(a,i8)') 'Successfully made organic matter '
    end if

  end subroutine mkorganic

end module mkorganicMod
