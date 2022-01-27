module mksoilcolMod

  use ESMF
  use pio
  use shr_kind_mod , only : r8 => shr_kind_r8, r4 => shr_kind_r4
  use shr_sys_mod  , only : shr_sys_abort
  use mkpioMod     , only : mkpio_get_rawdata
  use mkpioMod     , only : mkpio_iodesc_rawdata, pio_iotype, pio_ioformat, pio_iosystem
  use mkesmfMod    , only : regrid_rawdata
  use mkutilsMod   , only : chkerr, mkrank
  use mkvarctl     , only : root_task, ndiag, mpicom, MPI_INTEGER, MPI_MAX

  implicit none
  private

  public  :: mksoilcol      ! Set soil colors

  integer , parameter :: unsetcol  = -999      ! flag to indicate soil color NOT set
  integer , private   :: soil_color= unsetcol  ! soil color to override with

  character(len=*) , parameter :: u_FILE_u = &
       __FILE__

!=================================================================================
contains
!=================================================================================

  subroutine mksoilcol(file_data_i, file_mesh_i, mesh_o, soil_color_o, nsoilcol, rc)

    ! input/output variables
    character(len=*)  , intent(in)    :: file_mesh_i     ! input mesh file name
    character(len=*)  , intent(in)    :: file_data_i     ! input data file name
    type(ESMF_Mesh)   , intent(in)    :: mesh_o
    integer           , pointer       :: soil_color_o(:) ! soil color classes
    integer           , intent(out)   :: nsoilcol        ! number of soil colors 
    integer           , intent(out)   :: rc

    ! local variables:
    type(ESMF_RouteHandle) :: routehandle
    type(ESMF_Mesh)        :: mesh_i
    type(ESMF_Field)       :: field_i
    type(ESMF_Field)       :: field_o
    type(file_desc_t)      :: pioid
    integer                :: ni,no
    integer                :: ns_i, ns_o
    integer                :: n,l
    integer                :: rcode, ier             ! error status
    integer                :: srcMaskValue = -987987 ! spval for RH mask values
    integer                :: dstMaskValue = -987987 ! spval for RH mask values
    integer                :: srcTermProcessing_Value = 0
    real(r4), allocatable  :: soilcol_i(:)
    real(r4), allocatable  :: data_i(:,:)
    real(r4), allocatable  :: data_o(:,:)
    real(r4), allocatable  :: mask_i(:)
    logical                :: has_color    ! whether this grid cell has non-zero color
    integer                :: nsoilcol_local
    integer                :: maxindex(1)
    character(len=*), parameter :: subname = 'mksoilcol'
    !-----------------------------------------------------------------------

    rc = ESMF_SUCCESS

    ! Error check soil_color if it is set
    if ( soil_color /= unsetcol )then
       if ( soil_color < 0 .or. soil_color > 20 )then
          write(6,*)'soil_color is out of range = ', soil_color
          call shr_sys_abort()
       end if
       write(6,*) 'Replace soil color for all points with: ', soil_color
       do no = 1,size(soil_color_o)
          soil_color_o(no) = soil_color
       end do
       RETURN
    end if

    if (root_task) then
       write (ndiag,'(a)') 'Attempting to make soil color classes .....'
    end if

    ! Open raw data file - need to do this first to obtain ungridded dimension size
    call ESMF_VMLogMemInfo("Before pio_openfile for "//trim(file_data_i))
    rcode = pio_openfile(pio_iosystem, pioid, pio_iotype, trim(file_data_i), pio_nowrite)

    ! Determine ns_i and allocate data_i
    call ESMF_VMLogMemInfo("Before create mesh_i in "//trim(subname))
    mesh_i = ESMF_MeshCreate(filename=trim(file_mesh_i), fileformat=ESMF_FILEFORMAT_ESMFMESH, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMLogMemInfo("After create mesh_i in "//trim(subname))
    call ESMF_MeshGet(mesh_i, numOwnedElements=ns_i, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    allocate(soilcol_i(ns_i),stat=ier)
    if (ier/=0) call shr_sys_abort()

    ! Read in input data
    call mkpio_get_rawdata(pioid, 'SOIL_COLOR', mesh_i, soilcol_i, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMLogMemInfo("After mkpio_getrawdata in "//trim(subname))
 
    ! Determine maximum number of soil colors across all processors
    ! This will be used for the ungridded dimension of data_o below
    nsoilcol_local = maxval(soilcol_i)
    call mpi_allreduce(nsoilcol_local, nsoilcol, 1, MPI_INTEGER, MPI_MAX, mpicom, rcode)
    ! TODO: MPI_SUCCESS could not be accessed from mkvarctl
    !if (rcode /= MPI_SUCCESS) call shr_sys_abort('error for mpi_allredice in '//trim(subname))

    ! Determine ns_o and allocate data_o
    ns_o = size(soil_color_o)
    allocate(data_o(0:nsoilcol, ns_o),stat=ier)
    if (ier/=0) call shr_sys_abort()

    ! Determine input landmask (frac_i)
    allocate(mask_i(ns_i))
    call mkpio_get_rawdata(pioid, 'LANDMASK', mesh_i, mask_i, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! Now determine data_i as a real 2d array
    allocate(data_i(0:nsoilcol,ns_i))
    data_i(:,:) = 0._r4
    do l = 1,nsoilcol
       do n = 1,ns_i
          if (int(soilcol_i(n)) == l) then
             data_i(l,n) = 1._r4 * mask_i(n) 
          end if
       end do
    end do

    ! Create field on input mesh
    field_i = ESMF_FieldCreate(mesh_i, ESMF_TYPEKIND_R4, meshloc=ESMF_MESHLOC_ELEMENT, &
         ungriddedLbound=(/0/), ungriddedUbound=(/nsoilcol/), gridToFieldMap=(/2/), rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMLogMemInfo("After field_i creation in "//trim(subname))

    ! Create field on model model
    field_o = ESMF_FieldCreate(mesh_o, ESMF_TYPEKIND_R4, meshloc=ESMF_MESHLOC_ELEMENT, &
         ungriddedLbound=(/0/), ungriddedUbound=(/nsoilcol/), gridToFieldMap=(/2/), rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMLogMemInfo("After field_o creation in "//trim(subname))

    ! Create route handle to map field_model to field_data
    call ESMF_FieldRegridStore(field_i, field_o, routehandle=routehandle, &
         srcMaskValues=(/srcMaskValue/), dstMaskValues=(/dstMaskValue/), &
         regridmethod=ESMF_REGRIDMETHOD_CONSERVE, normType=ESMF_NORMTYPE_FRACAREA, &
         srcTermProcessing=srcTermProcessing_Value, &
         ignoreDegenerate=.true., unmappedaction=ESMF_UNMAPPEDACTION_IGNORE, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMLogMemInfo("After regridstore in "//trim(subname))

    ! Regrid data_i to data_o
    call regrid_rawdata(field_i, field_o, routehandle, data_i, data_o, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_LogWrite(subname//'after regrid rawdata')

    ! Determine dominant soil color in each output grid cell
    !call dominant_soil_color(data_o, nsoilcol, soil_color_o)

    soil_color_o(:) = 0
    do no = 1,ns_o
       ! If the output cell has any non-zero-colored inputs, then set the weight of
       ! zero-colored inputs to 0, to ensure that the zero-color is NOT dominant.
       if (any(data_o(1:nsoilcol,no) > 0.)) then
          has_color = .true.
          data_o(0,no) = 0.0
       else
          has_color = .false.
       end if

       ! Rank non-zero weights by color type. wsti(1) is the most extensive color type. 
       if (has_color) then
          !maxindex = maxloc(data_o(0:nsoilcol,no))
          !soil_color_o(no) = maxindex
          call mkrank (nsoilcol, data_o(0:nsoilcol,no), 9999, 1, maxindex)
          soil_color_o(no) = maxindex(1)
       end if

       ! If land but no color, set color to 15 (in older dataset generic soil color 4)
       if (nsoilcol == 8) then
          if (soil_color_o(no)==0) then
             soil_color_o(no) = 4
          end if
       else if (nsoilcol == 20) then
          if (soil_color_o(no)==0) then
             soil_color_o(no) = 15
          end if
       else
          write(6,*) 'MKSOILCOL error: unhandled nsoilcol: ', nsoilcol
          call shr_sys_abort()
       end if

       ! Error checks
       if (soil_color_o(no) < 0 .or. soil_color_o(no) > nsoilcol) then
          write (6,*) 'MKSOILCOL error: land model soil color = ', &
               soil_color_o(no),' is not valid for lon,lat = ',no
          call shr_sys_abort()
       end if

    end do

    if (root_task) then
       write (ndiag,'(a)') 'Successfully made soil color classes'
       write (ndiag,'(a)')
    end if

    deallocate(data_i)
    deallocate(data_o)
    deallocate(soilcol_i)
    deallocate(mask_i)
    call ESMF_FieldDestroy(field_i, nogarbage = .true., rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort()
    call ESMF_FieldDestroy(field_o, nogarbage = .true., rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort()
    call ESMF_RouteHandleDestroy(routehandle, nogarbage = .true., rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort()

    call ESMF_LogWrite(subname//' finished routine mksoilcol')

  end subroutine mksoilcol

end module mksoilcolMod
