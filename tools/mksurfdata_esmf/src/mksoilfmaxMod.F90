module mksoilfmaxMod

  !-----------------------------------------------------------------------
  ! Make soil fmax
  !-----------------------------------------------------------------------

  use ESMF
  use shr_kind_mod     , only : r8 => shr_kind_r8, r4=>shr_kind_r4
  use shr_sys_mod      , only : shr_sys_abort
  use pio              , only : file_desc_t, pio_openfile, pio_closefile, pio_nowrite, pio_syncfile
  use mkpioMod         , only : mkpio_get_rawdata, pio_iotype, pio_iosystem
  use mkesmfMod        , only : regrid_rawdata, create_routehandle_r8
  use mkdiagnosticsMod , only : output_diagnostics_area
  use mkvarctl         , only : ndiag, root_task, spval
  use mkutilsMod       , only : chkerr
  use mkfileMod        , only : mkfile_output

  implicit none
  private           ! By default make data private

  public :: mksoilfmax         ! Make percent fmax

  character(len=*) , parameter :: u_FILE_u = &
       __FILE__

!=================================================================================
contains
!=================================================================================

  subroutine mksoilfmax(file_mesh_i, file_data_i, mesh_o, pioid_o, rc)
    !
    ! make percent fmax
    !
    ! input/output variables
    character(len=*)  , intent(in)    :: file_mesh_i ! input mesh file name
    character(len=*)  , intent(in)    :: file_data_i ! input data file name
    type(ESMF_Mesh)   , intent(in)    :: mesh_o      ! output mesh
    type(file_desc_t) , intent(inout) :: pioid_o
    integer           , intent(out)   :: rc

    ! local variables:
    type(ESMF_RouteHandle) :: routehandle
    type(ESMF_Mesh)        :: mesh_i
    type(file_desc_t)      :: pioid_i
    integer                :: ni,no
    integer                :: ns_i, ns_o
    integer                :: n,l,k
    integer , allocatable  :: mask_i(:)
    real(r8), allocatable  :: frac_i(:)
    real(r8), allocatable  :: frac_o(:)
    real(r8), allocatable  :: fmax_i(:)  ! input grid: percent fmax
    real(r8), allocatable  :: fmax_o(:)  ! output grid: %fmax
    integer                :: ier, rcode ! error status
    character(len=32)      :: subname = 'mksoilfmax'
    !-----------------------------------------------------------------------

    rc = ESMF_SUCCESS

    if (root_task) then
       write(ndiag,*)
       write(ndiag,'(1x,80a1)') ('=',k=1,80)
       write(ndiag,*)
       write(ndiag,'(a)') 'Attempting to make %fmax .....'
       write(ndiag,'(a)') ' Input file is '//trim(file_data_i)
       write(ndiag,'(a)') ' Input mesh file is '//trim(file_mesh_i)
    end if
    call ESMF_VMLogMemInfo("At start of "//trim(subname))

    ! Open input data file
    rcode = pio_openfile(pio_iosystem, pioid_i, pio_iotype, trim(file_data_i), pio_nowrite)

    ! Read in input mesh
    mesh_i = ESMF_MeshCreate(filename=trim(file_mesh_i), fileformat=ESMF_FILEFORMAT_ESMFMESH, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMLogMemInfo("After create mesh_i in "//trim(subname))

    ! Determine ns_i
    call ESMF_MeshGet(mesh_i, numOwnedElements=ns_i, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Determine ns_o and allocate output data
    call ESMF_MeshGet(mesh_o, numOwnedElements=ns_o, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    allocate(fmax_o(ns_o)); fmax_o(:) = spval

    ! Get the landmask from the file and reset the mesh mask based on that
    allocate(frac_i(ns_i), stat=ier)
    if (ier/=0) call shr_sys_abort()
    allocate(mask_i(ns_i), stat=ier)
    if (ier/=0) call shr_sys_abort()
    call mkpio_get_rawdata(pioid_i, 'LANDMASK', mesh_i, frac_i, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    do ni = 1,ns_i
       if (frac_i(ni) > 0._r8) then
          mask_i(ni) = 1
       else
          mask_i(ni) = 0
       end if
    end do
    call ESMF_MeshSet(mesh_i, elementMask=mask_i, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! Create a route handle between the input and output mesh
    allocate(frac_o(ns_o))
    call create_routehandle_r8(mesh_i=mesh_i, mesh_o=mesh_o, norm_by_fracs=.true., &
         routehandle=routehandle, frac_o=frac_o, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMLogMemInfo("After create routehandle in "//trim(subname))
    do n = 1, ns_o
       if ((frac_o(n) < 0.0) .or. (frac_o(n) > 1.0001)) then
          write(6,*) "ERROR:: frac_o out of range: ", frac_o(n),n
       end if
    end do

    ! Read in input data
    allocate(fmax_i(ns_i), stat=ier)
    if (ier/=0) call shr_sys_abort()
    call mkpio_get_rawdata(pioid_i, 'FMAX', mesh_i, fmax_i, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMLogMemInfo("After mkpio_getrawdata in "//trim(subname))

    ! Regrid fmax_i to fmax_o, in points with no data, use globalAvg
    fmax_i(:) = fmax_i(:) * frac_i(:)
    fmax_o(:) = 0._r8
    call regrid_rawdata(mesh_i, mesh_o, routehandle, fmax_i, fmax_o, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_LogWrite(subname//'after regrid rawdata in '//trim(subname))
    do n = 1,ns_o
       if (frac_o(n) == 0._r8) then
          fmax_o(n) = .365783_r8
       end if
    end do

    ! Check for conservation
    do no = 1, ns_o
       if ((fmax_o(no)) > 1.000001_r8) then
          write (6,*) 'MKFMAX error: fmax = ',fmax_o(no),' greater than 1.000001 for no = ',no
          call shr_sys_abort()
       end if
    enddo

    ! Compare global areas on input and output grids
    call output_diagnostics_area(mesh_i, mesh_o, mask_i, frac_o, &
         fmax_i*0.01_r8, fmax_o*0.01_r8, "Max Fractional Sataturated Area", &
         percent=.false., ndiag=ndiag, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Write output data
    if (root_task)  write(ndiag, '(a)') trim(subname)//" writing out soil fmax (maximum fraction saturated area)"
    call mkfile_output (pioid_o,  mesh_o,  'FMAX', fmax_o, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in calling mkfile_output')
    call pio_syncfile(pioid_o)

    ! Close the input file
    call pio_closefile(pioid_i)

    ! Release memory
    call ESMF_RouteHandleDestroy(routehandle, nogarbage = .true., rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort()
    call ESMF_MeshDestroy(mesh_i, nogarbage = .true., rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort()

    if (root_task) then
       write (ndiag,'(a)') 'Successfully made %fmax'
    end if
    call ESMF_VMLogMemInfo("At end of "//trim(subname))

  end subroutine mksoilfmax

end module mksoilfmaxMod
