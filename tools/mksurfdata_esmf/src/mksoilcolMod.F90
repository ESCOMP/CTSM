module mksoilcolMod

  use ESMF
  use pio              , only : file_desc_t, pio_openfile, pio_closefile, pio_nowrite
  use pio              , only : pio_syncfile, pio_inq_varid, pio_put_var, var_desc_t
  use shr_kind_mod     , only : r8 => shr_kind_r8, r4 => shr_kind_r4
  use shr_sys_mod      , only : shr_sys_abort
  use mkpioMod         , only : mkpio_get_rawdata, pio_iotype, pio_iosystem
  use mkvarctl         , only : root_task, ndiag, mpicom, unsetcol
  use mkdiagnosticsMod , only : output_diagnostics_index
  use mkutilsMod       , only : chkerr
  use mkfileMod        , only : mkfile_output

  implicit none
  private

#include <mpif.h>

  public  :: mksoilcol      ! Set soil colors
  private :: get_dominant_soilcol

  integer                :: num_soilcolors
  type(ESMF_DynamicMask) :: dynamicMask

  character(len=*) , parameter :: u_FILE_u = &
       __FILE__

!=================================================================================
contains
!=================================================================================

  subroutine mksoilcol(file_data_i, file_mesh_i, mesh_o, pioid_o, rc)

    ! input/output variables
    character(len=*)  , intent(in)    :: file_mesh_i     ! input mesh file name
    character(len=*)  , intent(in)    :: file_data_i     ! input data file name
    type(ESMF_Mesh)   , intent(in)    :: mesh_o          ! model mesho
    type(file_desc_t) , intent(inout) :: pioid_o
    integer           , intent(out)   :: rc

    ! local variables:
    type(ESMF_RouteHandle) :: routehandle
    type(ESMF_Mesh)        :: mesh_i
    type(ESMF_Field)       :: field_i
    type(ESMF_Field)       :: field_o
    type(ESMF_Field)       :: field_dstfrac
    type(file_desc_t)      :: pioid_i
    type(var_desc_t)       :: pio_varid
    integer                :: ni,no, k
    integer                :: ns_i, ns_o
    integer , allocatable  :: mask_i(:)
    real(r4), allocatable  :: rmask_i(:)
    real(r8), allocatable  :: frac_o(:)
    real(r4), allocatable  :: soil_color_i(:)
    integer , allocatable  :: soil_color_o(:) ! soil color classes
    integer                :: nsoilcol        ! number of soil colors
    real(r4), pointer      :: dataptr(:)
    real(r8), pointer      :: dataptr_r8(:)
    integer                :: nsoilcol_local
    integer                :: rcode, ier
    integer                :: srcTermProcessing_Value = 0
    character(len=*), parameter :: subname = 'mksoilcol'
    !-----------------------------------------------------------------------

    rc = ESMF_SUCCESS

    ! Note soil_color_override has been removed - instead should now use tools
    ! subset_data and modify_fsurdat

    if (root_task) then
       write(ndiag,*)
       write(ndiag,'(1x,80a1)') ('=',k=1,80)
       write(ndiag,*)
       write(ndiag,'(a)') 'Attempting to make soil color classes .....'
       write(ndiag,'(a)') ' Input file is '//trim(file_data_i)
       write(ndiag,'(a)') ' Input mesh file is '//trim(file_mesh_i)
    end if

    ! Open soil color data file
    call ESMF_VMLogMemInfo("Before pio_openfile for "//trim(file_data_i))
    rcode = pio_openfile(pio_iosystem, pioid_i, pio_iotype, trim(file_data_i), pio_nowrite)

    ! Read in input mesh
    call ESMF_VMLogMemInfo("Before create mesh_i in "//trim(subname))
    mesh_i = ESMF_MeshCreate(filename=trim(file_mesh_i), fileformat=ESMF_FILEFORMAT_ESMFMESH, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMLogMemInfo("After create mesh_i in "//trim(subname))

    ! Determine ns_i
    call ESMF_MeshGet(mesh_i, numOwnedElements=ns_i, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Determine ns_o and allocate output data
    call ESMF_MeshGet(mesh_o, numOwnedElements=ns_o, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    allocate (soil_color_o(ns_o)); soil_color_o(:) = -999

    ! Get the landmask from the file and reset the mesh mask based on that
    allocate(rmask_i(ns_i), stat=ier)
    if (ier/=0) call shr_sys_abort()
    allocate(mask_i(ns_i), stat=ier)
    if (ier/=0) call shr_sys_abort()
    call mkpio_get_rawdata(pioid_i, 'LANDMASK', mesh_i, rmask_i, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    do ni = 1,ns_i
       if (rmask_i(ni) > 0._r4) then
          mask_i(ni) = 1
       else
          mask_i(ni) = 0
       end if
    end do
    deallocate(rmask_i)
    call ESMF_MeshSet(mesh_i, elementMask=mask_i, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! Read in input soil color data
    allocate(soil_color_i(ns_i),stat=ier)
    if (ier/=0) call shr_sys_abort()
    call mkpio_get_rawdata(pioid_i, 'SOIL_COLOR', mesh_i, soil_color_i, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMLogMemInfo("After mkpio_getrawdata in "//trim(subname))

    ! Scale the input soil color by the input mask
    do ni = 1,ns_i
       if (mask_i(ni) == 0) then
          soil_color_i(ni) = 0._r4
       end if
    end do

    ! Determine maximum number of soil colors across all processors
    nsoilcol_local = maxval(soil_color_i)
    call mpi_allreduce(nsoilcol_local, nsoilcol, 1, MPI_INTEGER, MPI_MAX, mpicom, rcode)

    ! Set module variable (used in the get_dominant_soilcol routine)
    num_soilcolors = nsoilcol

     ! Create ESMF fields that will be used below
    field_i = ESMF_FieldCreate(mesh_i, ESMF_TYPEKIND_R4, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    field_o = ESMF_FieldCreate(mesh_o, ESMF_TYPEKIND_R4, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    field_dstfrac = ESMF_FieldCreate(mesh_o, ESMF_TYPEKIND_R8, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! Create a route handle
    call ESMF_FieldRegridStore(field_i, field_o, routehandle=routehandle, &
         regridmethod=ESMF_REGRIDMETHOD_CONSERVE, srcTermProcessing=srcTermProcessing_Value, &
         ignoreDegenerate=.true., unmappedaction=ESMF_UNMAPPEDACTION_IGNORE, &
         dstFracField= field_dstfrac, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMLogMemInfo("After regridstore in "//trim(subname))

    ! Determine frac_o
    call ESMF_FieldGet(field_dstfrac, farrayptr=dataptr_r8, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    allocate(frac_o(ns_o))
    frac_o(:) = dataptr_r8(:)

    ! Create a dynamic mask object
    ! The dynamic mask object further holds a pointer to the routine that will be called in order to
    ! handle dynamically masked elements - in this case its DynMaskProc (see below)
    call ESMF_DynamicMaskSetR4R8R4(dynamicMask, dynamicMaskRoutine=get_dominant_soilcol,  &
         dynamicSrcMaskValue=0._r4, handleAllElements=.true., rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! Determine dominant soil color in the field regrid call below
    call ESMF_FieldGet(field_i, farrayptr=dataptr, rc=rc)
    dataptr(:) = soil_color_i(:)
    call ESMF_FieldGet(field_o, farrayptr=dataptr, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    dataptr(:) = 0._r4

    call ESMF_FieldRegrid(field_i, field_o, routehandle=routehandle, dynamicMask=dynamicMask, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call ESMF_FieldGet(field_o, farrayptr=dataptr, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    do no = 1,ns_o
       soil_color_o(no) = int(dataptr(no))
       if (soil_color_o(no) < 0 .or. soil_color_o(no) > nsoilcol) then
          write (6,*) 'MKSOILCOL error: land model soil color = ', &
               soil_color_o(no),' is not valid for lon,lat = ',no
          call shr_sys_abort()
       end if
    end do

    ! Write output data
    if (root_task)  write(ndiag, '(a)') trim(subname)//" writing out soil color"
    call mkfile_output(pioid_o,  mesh_o, 'SOIL_COLOR', soil_color_o,  rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in calling mkfile_output')
    if (root_task)  write(ndiag, '(a)') trim(subname)//" writing out mksoil_color"
    rcode = pio_inq_varid(pioid_o, 'mxsoil_color', pio_varid)
    rcode = pio_put_var(pioid_o, pio_varid, nsoilcol)
    call pio_syncfile(pioid_o)

    ! Compare global area of each soil color on input and output grids
    call output_diagnostics_index(mesh_i, mesh_o, mask_i, frac_o, &
         0, nsoilcol, int(soil_color_i), soil_color_o, 'soil color type',  ndiag, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! Close the input file 
    call pio_closefile(pioid_i)

    ! Clean up memory
    call ESMF_RouteHandleDestroy(routehandle, nogarbage = .true., rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort()
    call ESMF_MeshDestroy(mesh_i, nogarbage = .true., rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort()
    call ESMF_FieldDestroy(field_i, nogarbage = .true., rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort()
    call ESMF_FieldDestroy(field_o, nogarbage = .true., rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort()
    call ESMF_FieldDestroy(field_dstfrac, nogarbage = .true., rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort()
    call ESMF_VMLogMemInfo("After destroy operations in "//trim(subname))

    if (root_task) then
       write (ndiag,'(a)') 'Successfully made soil color classes'
       write (ndiag,'(a)')
    end if

  end subroutine mksoilcol

  !================================================================================================
  subroutine get_dominant_soilcol(dynamicMaskList, dynamicSrcMaskValue, dynamicDstMaskValue, rc)

    ! input/output arguments
    type(ESMF_DynamicMaskElementR4R8R4) , pointer              :: dynamicMaskList(:)
    real(ESMF_KIND_R4)                  , intent(in), optional :: dynamicSrcMaskValue
    real(ESMF_KIND_R4)                  , intent(in), optional :: dynamicDstMaskValue
    integer                             , intent(out)          :: rc

    ! local variables
    integer            :: ni, no, n
    real(ESMF_KIND_R4) :: wts_o(0:num_soilcolors)
    logical            :: has_color
    integer            :: soil_color_o
    integer            :: maxindex(1)
    !---------------------------------------------------------------

    rc = ESMF_SUCCESS

    if (associated(dynamicMaskList)) then
       do no = 1, size(dynamicMaskList)
          wts_o(:) = 0.d0
          do ni = 1, size(dynamicMaskList(no)%factor)
             if (dynamicSrcMaskValue /= dynamicMaskList(no)%srcElement(ni)) then
                do n = 0,num_soilcolors
                   if (dynamicMaskList(no)%srcElement(ni) == n) then
                      wts_o(n) = wts_o(n) + dynamicMaskList(no)%factor(ni)
                   end if
                enddo
             end if
          end do

          ! If the output cell has any non-zero-colored inputs, then set the weight of
          ! zero-colored inputs to 0, to ensure that the zero-color is NOT dominant.
          soil_color_o = 0
          if (any(wts_o(1:num_soilcolors) > 0.)) then
             has_color = .true.
             wts_o(0) = 0.0
          else
             has_color = .false.
          end if

          ! Find index of maximum weight
          if (has_color) then
             call mkrank (num_soilcolors, wts_o(0:num_soilcolors), maxindex)
             soil_color_o = maxindex(1)
          end if

          ! If no color, set color to 15 (in older dataset generic soil color 4)
          if (num_soilcolors == 8) then
             if (soil_color_o == 0) then
                soil_color_o = 4
             end if
          else if (num_soilcolors == 20) then
             if (soil_color_o == 0) then
                soil_color_o = 15
             end if
          end if
          dynamicMaskList(no)%dstElement = real(soil_color_o, kind=r4)

       end do
    end if

    contains

      subroutine mkrank (n, a, iv)
        ! Return indices of largest [num] values in array [a].

        ! input/output variables
        integer , intent(in) :: n      !array length
        real(r4), intent(in) :: a(0:n) !array to be ranked
        integer , intent(out):: iv(1)  !index to [num] largest values in array [a]

        ! local variables:
        real(r4) :: a_max  !maximum value in array
        real(r4) :: delmax !tolerance for finding if larger value
        integer  :: i      !array index
        integer  :: m      !do loop index
        integer  :: k      !do loop index
        integer  :: miss   !missing data value
        !-----------------------------------------------------------------------

        ! Find index of largest non-zero number
        delmax = 1.e-06
        miss = 9999
        iv(1) = miss

        a_max = -9999.
        do i = 0, n
           if (a(i)>0. .and. (a(i)-a_max)>delmax) then
              a_max = a(i)
              iv(1)  = i
           end if
        end do
        ! iv(1) = miss indicates no values > 0. this is an error
        if (iv(1) == miss) then
           write (6,*) 'MKRANK error: iv(1) = missing'
           call shr_sys_abort()
        end if
      end subroutine mkrank

  end subroutine get_dominant_soilcol

end module mksoilcolMod
