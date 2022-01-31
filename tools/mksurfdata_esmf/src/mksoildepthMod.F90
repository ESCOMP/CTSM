module mksoildepthMod

  !-----------------------------------------------------------------------
  ! make fraction soildepth from input soildepth data
  !-----------------------------------------------------------------------
  !
  use ESMF
  use pio
  use shr_kind_mod   , only : r8 => shr_kind_r8, r4=>shr_kind_r4, cs => shr_kind_cs
  use shr_sys_mod    , only : shr_sys_abort
  use mkpioMod       , only : mkpio_get_rawdata
  use mkpioMod       , only : pio_iotype, pio_ioformat, pio_iosystem
  use mkesmfMod      , only : regrid_rawdata, create_routehandle_r8
  use mkutilsMod     , only : chkerr
  use mkchecksMod    , only : min_bad, max_bad  
  use mkvarctl

  implicit none
  private

  public mksoildepth           ! regrid soildepth data

  character(len=*) , parameter :: u_FILE_u = &
       __FILE__

!===============================================================
contains
!===============================================================

  subroutine mksoildepth(file_mesh_i, file_data_i, mesh_o, soildepth_o, rc)
    !
    ! make soildepth
    !
    ! input/output variables
    character(len=*)  , intent(in)  :: file_mesh_i ! input mesh file name
    character(len=*)  , intent(in)  :: file_data_i ! input data file name
    type(ESMF_Mesh)   , intent(in)  :: mesh_o      ! output mesh
    real(r8)          , intent(out) :: soildepth_o(:) ! output grid: fraction soildepth
    integer           , intent(out) :: rc

    ! local variables:
    type(ESMF_RouteHandle)      :: routehandle
    type(ESMF_Mesh)             :: mesh_i
    type(file_desc_t)           :: pioid
    integer                     :: ns_i, ns_o
    integer                     :: n
    real(r8), allocatable       :: soildepth_i(:)  
    real(r8), allocatable       :: frac_i(:)
    real(r8), allocatable       :: frac_o(:)
    integer                     :: ier, rcode                ! error status
    character(len=CS)           :: varname 
    integer                     :: varnum
    real(r8), parameter         :: min_valid = 0._r8         ! minimum valid value
    real(r8), parameter         :: max_valid = 100.000001_r8 ! maximum valid value
    character(len=*), parameter :: subname = 'mksoildepth'
    !-----------------------------------------------------------------------

    rc = ESMF_SUCCESS

    if (root_task) then
       write (ndiag,'(a)') 'Attempting to make soildepth .....'
    end if

    ! Open input data file
    call ESMF_VMLogMemInfo("Before pio_openfile for "//trim(file_data_i))
    rcode = pio_openfile(pio_iosystem, pioid, pio_iotype, trim(file_data_i), pio_nowrite)

    ! Read in input mesh
    call ESMF_VMLogMemInfo("Before create mesh_i in "//trim(subname))
    mesh_i = ESMF_MeshCreate(filename=trim(file_mesh_i), fileformat=ESMF_FILEFORMAT_ESMFMESH, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMLogMemInfo("After create mesh_i in "//trim(subname))

    ! Create a route handle between the input and output mesh
    call create_routehandle_r8(mesh_i, mesh_o, routehandle, rc)
    call ESMF_VMLogMemInfo("After create routehandle in "//trim(subname))

    ! Determine ns_i and allocate soildepth_i
    call ESMF_MeshGet(mesh_i, numOwnedElements=ns_i, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    allocate(soildepth_i(ns_i), stat=ier)
    if (ier/=0) call shr_sys_abort()

    ! Determine frac_o (regrid frac_i to frac_o)
    allocate(frac_i(ns_i))
    allocate(frac_o(ns_o))
    call mkpio_get_rawdata(pioid, 'LANDMASK', mesh_i, frac_i, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call regrid_rawdata(mesh_i, mesh_o, routehandle, frac_i, frac_o, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    do n = 1, ns_o
       if ((frac_o(n) < 0.0) .or. (frac_o(n) > 1.0001)) then
          write(6,*) "ERROR:: frac_o out of range: ", frac_o(n),n
          call shr_sys_abort ()
       end if
    end do
    call ESMF_VMLogMemInfo("After regrid landmask in  "//trim(subname))

    ! Determine variable name to read in
    varnum = 1
    select case (varnum)
    case(1)
       varname = 'Avg_Depth_Median'
    case(2)
       varname = 'Avg_Depth_Mean'
    case(3)
       varname = 'Upland_Valley_Depth_Median'
    case(4)
       varname = 'Upland_Valley_Depth_Mean'
    case(5)
       varname = 'Upland_Hillslope_Depth_Median'
    case(6)
       varname = 'Upland_Hillslope_Depth_Mean'
    case(7)
       varname = 'Lowland_Depth_Mean'
    case(8)
       varname = 'Lowland_Depth_Mean'
    end select

    ! Read in input data
    call mkpio_get_rawdata(pioid, trim(varname), mesh_i, soildepth_i, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMLogMemInfo("After mkpio_getrawdata in "//trim(subname))

    ! Regrid soildepth_i to soildepth_o and scale by 1/frac_o
    call regrid_rawdata(mesh_i, mesh_o, routehandle, soildepth_i, soildepth_o, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_LogWrite(subname//'after regrid rawdata in '//trim(subname))
    do n = 1,ns_o
       if (frac_o(n) > 0._r8) then
          soildepth_o(n) = soildepth_o(n) / frac_o(n)
       else
          soildepth_o(n) = 0._r8
       end if
    end do

    ! Check validity of output data
    if ( min_bad(soildepth_o, min_valid, 'soildepth') .or. max_bad(soildepth_o, max_valid, 'soildepth')) then
       call shr_sys_abort()
    end if

    ! Close the input file
    call pio_closefile(pioid)
    call ESMF_VMLogMemInfo("After pio_closefile in "//trim(subname))

    ! Release memory
    deallocate(soildepth_i)
    deallocate(frac_i)
    deallocate(frac_o)
    call ESMF_RouteHandleDestroy(routehandle, nogarbage = .true., rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort()
    call ESMF_MeshDestroy(mesh_i, nogarbage = .true., rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort()
    call ESMF_VMLogMemInfo("After destroy operations in "//trim(subname))

    if (root_task) then
       write (ndiag,'(a)') 'Successfully made soildepth'
       write (ndiag,*)
    end if

  end subroutine mksoildepth

end module mksoildepthMod
