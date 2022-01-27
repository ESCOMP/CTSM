module mksoilcolMod

  use ESMF
  use pio
  use shr_kind_mod , only : r8 => shr_kind_r8, r4 => shr_kind_r4
  use shr_sys_mod  , only : shr_sys_abort
  use mkpioMod     , only : mkpio_get_rawdata
  use mkpioMod     , only : mkpio_iodesc_rawdata, pio_iotype, pio_ioformat, pio_iosystem
  use mkesmfMod    , only : regrid_rawdata
  use mkutilsMod   , only : chkerr, mkrank
  use mkvarctl     , only : root_task, ndiag

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
    integer                :: ndims
    real(r8), allocatable  :: soilcol_i(:)
    real(r8), allocatable  :: data_i(:,:)
    real(r8), allocatable  :: data_o(:,:)
    real(r8), allocatable  :: mask_i(:)
    integer, parameter          :: num = 2      ! set soil mapunit number
    integer                     :: wsti(num)    ! index to 1st and 2nd largest wst
    real(r8)                    :: wt           ! map overlap weight
    logical                     :: has_color    ! whether this grid cell has non-zero color
    integer, parameter          :: miss = 99999 ! missing data indicator
    character(len=35), allocatable :: col(:)  ! name of each color
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
 
    ! ** TODO: determine the maximum value of soilcol_i - across all processors - doing an all reduce
    ! For now hardwire nsoilcol to 20
    nsoilcol = 20

    ! Determine ns_o and allocate data_o
    ns_o = size(soil_color_o)
    allocate(data_o(0:nsoilcol, ns_o),stat=ier)
    if (ier/=0) call shr_sys_abort()

    ! Define the model color classes: 0 to nsoilcol
    allocate(col(0:nsoilcol))
    if (nsoilcol == 20) then
       col(0)  = 'no soil                            '
       col(1)  = 'class 1: light                     '
       col(2)  = 'class 2:                           '
       col(3)  = 'class 3:                           '
       col(4)  = 'class 4:                           '
       col(5)  = 'class 5:                           '
       col(6)  = 'class 6:                           '
       col(7)  = 'class 7:                           '
       col(8)  = 'class 8:                           '
       col(9)  = 'class 9:                           '
       col(10) = 'class 10:                          '
       col(11) = 'class 11:                          '
       col(12) = 'class 12:                          '
       col(13) = 'class 13:                          '
       col(14) = 'class 14:                          '
       col(15) = 'class 15:                          '
       col(16) = 'class 16:                          '
       col(17) = 'class 17:                          '
       col(18) = 'class 18:                          '
       col(19) = 'class 19:                          '
       col(20) = 'class 20: dark                     '
    else if (nsoilcol == 8) then
       col(0) = 'no soil                            '
       col(1) = 'class 1: light                     '
       col(2) = 'class 2:                           '
       col(3) = 'class 3:                           '
       col(4) = 'class 4:                           '
       col(5) = 'class 5:                           '
       col(6) = 'class 6:                           '
       col(7) = 'class 7:                           '
       col(8) = 'class 8: dark                      '
    else
       write(6,*)'nsoilcol value of ',nsoilcol,' is not currently supported'
       call shr_sys_abort()
    end if

    ! Create field on input mesh
    field_i = ESMF_FieldCreate(mesh_i, ESMF_TYPEKIND_R8, meshloc=ESMF_MESHLOC_ELEMENT, &
         ungriddedLbound=(/0/), ungriddedUbound=(/nsoilcol/), gridToFieldMap=(/2/), rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMLogMemInfo("After field_i creation in "//trim(subname))

    ! Create field on model model
    field_o = ESMF_FieldCreate(mesh_o, ESMF_TYPEKIND_R8, meshloc=ESMF_MESHLOC_ELEMENT, &
         ungriddedLbound=(/0/), ungriddedUbound=(/nsoilcol/), gridToFieldMap=(/2/), rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMLogMemInfo("After field_o creation in "//trim(subname))

    ! Determine input landmask (frac_i)
    allocate(mask_i(ns_i))
    call mkpio_get_rawdata(pioid, 'LANDMASK', mesh_i, mask_i, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! Now determine data_i as a real 2d array
    allocate(data_i(0:nsoilcol,ns_i))
    data_i(:,:) = 0._r8
    do l = 1,nsoilcol
       do n = 1,ns_i
          if (int(soilcol_i(n)) == l) then
             data_i(l,n) = 1._r8 * mask_i(n) 
          end if
       end do
    end do

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
          call mkrank (nsoilcol, data_o(0:nsoilcol,no), miss, num, wsti)
          soil_color_o(no) = wsti(1)
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

    deallocate(col)
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
