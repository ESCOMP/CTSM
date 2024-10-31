module mklanwatMod

  !-----------------------------------------------------------------------
  ! make %lake and %wetland from input lake / wetland data
  ! also make lake parameters
  !-----------------------------------------------------------------------

  use ESMF
  use pio              , only : file_desc_t, pio_openfile, pio_closefile, pio_nowrite, pio_syncfile
  use shr_kind_mod     , only : r8 => shr_kind_r8, r4 => shr_kind_r4
  use shr_sys_mod      , only : shr_sys_abort
  use mkpioMod         , only : mkpio_get_rawdata, pio_iotype, pio_iosystem
  use mkesmfMod        , only : regrid_rawdata, create_routehandle_r8
  use mkdiagnosticsMod , only : output_diagnostics_continuous, output_diagnostics_area
  use mkchecksMod      , only : min_bad
  use mkutilsMod       , only : chkerr
  use mkvarctl         , only : root_task, ndiag, spval, no_inlandwet
  use mkfileMod        , only : mkfile_output

  implicit none
  private

  public :: mkpctlak    ! make %lake
  public :: mklakdep    ! make lake depth
  public :: mkwetlnd    ! make %wetland
  public :: update_max_array_lake  ! Update the maximum lake percent

  real(r8), allocatable  :: frac_o_mklak_nonorm(:)
  type(ESMF_RouteHandle) :: routehandle_mklak_nonorm

  character(len=*) , parameter :: u_FILE_u = &
       __FILE__

!===============================================================
contains
!===============================================================

  subroutine mkpctlak(file_mesh_i, file_data_i, mesh_o, lake_o, pioid_o, rc)

    ! -------------------
    ! make %lake
    ! PCT_LAKE is written to fsurdat in mksurfdata after adjustments are made
    ! -------------------

    ! uses
    use mkinputMod, only: mksrf_fdynuse

    ! input/output variables
    character(len=*)  , intent(in)    :: file_mesh_i      ! input mesh file name
    character(len=*)  , intent(in)    :: file_data_i      ! input data file name
    type(ESMF_Mesh)   , intent(in)    :: mesh_o
    type(file_desc_t) , intent(inout) :: pioid_o
    real(r8)          , intent(out)   :: lake_o(:)        ! output grid: %lake
    integer           , intent(out)   :: rc

    ! local variables
    type(ESMF_Mesh)        :: mesh_i
    type(file_desc_t)      :: pioid_i
    integer , allocatable  :: mask_i(:)
    real(r8), allocatable  :: rmask_i(:)
    real(r8), allocatable  :: lake_i(:)      ! input grid: percent lake
    integer                :: ni,no,k        ! indices
    integer                :: ns_i,ns_o      ! local sizes
    integer                :: ier,rcode      ! error status
    character(len=*), parameter :: subname = ' mkpctlak '
    !-----------------------------------------------------------------------

    rc = ESMF_SUCCESS

    if (root_task) then
       write(ndiag,*)
       write(ndiag,'(1x,80a1)') ('=',k=1,80)
       write(ndiag,*)
       write(ndiag,'(a)')'Attempting to make %lake'
       write(ndiag,'(a)') ' Input file is '//trim(file_data_i)
       write(ndiag,'(a)') ' Input mesh file is '//trim(file_mesh_i)
    end if

    ! ----------------------------------------
    ! Read i input data and input mesh and create route handle
    ! ----------------------------------------

    ! Open raw data file
    rcode = pio_openfile(pio_iosystem, pioid_i, pio_iotype, trim(file_data_i), pio_nowrite)

    ! Read in input mesh
    mesh_i = ESMF_MeshCreate(filename=trim(file_mesh_i), fileformat=ESMF_FILEFORMAT_ESMFMESH, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMLogMemInfo("After create mesh_i in "//trim(subname))

    ! Determine ns_i
    call ESMF_MeshGet(mesh_i, numOwnedElements=ns_i, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Determine ns_o
    call ESMF_MeshGet(mesh_o, numOwnedElements=ns_o, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Get the landmask from the input data file and reset the mesh mask based on that
    allocate(rmask_i(ns_i), stat=ier)
    if (ier/=0) call shr_sys_abort()
    allocate(mask_i(ns_i), stat=ier)
    if (ier/=0) call shr_sys_abort()
    call mkpio_get_rawdata(pioid_i, 'LANDMASK', mesh_i, rmask_i, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    do ni = 1,ns_i
       if (rmask_i(ni) > 0._r8) then
          mask_i(ni) = 1
       else
          mask_i(ni) = 0
       end if
    end do
    call ESMF_MeshSet(mesh_i, elementMask=mask_i, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! Create a route handle between the input and output mesh
    if (.not. ESMF_RouteHandleIsCreated(routehandle_mklak_nonorm)) then
       allocate(frac_o_mklak_nonorm(ns_o))
       ! Note that norm_by_fracs is false in the following because this routehandle is
       ! used to map fields that are expressed in terms of % of the grid cell.
       call create_routehandle_r8(mesh_i=mesh_i, mesh_o=mesh_o, norm_by_fracs=.false., &
            routehandle=routehandle_mklak_nonorm, frac_o=frac_o_mklak_nonorm, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_VMLogMemInfo("After create routehandle in "//trim(subname))
    end if

    ! ----------------------------------------
    ! Create %lake
    ! ----------------------------------------

    lake_o(:) = 0._r8
    if (root_task) then
       write (ndiag,*) 'Attempting to make %lake .....'
    end if

    ! Read in lake_i
    allocate(lake_i(ns_i), stat=rcode)
    if (rcode/=0) call shr_sys_abort()
    call mkpio_get_rawdata(pioid_i, 'PCT_LAKE', mesh_i, lake_i, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Regrid lake_i to lake_o
    call regrid_rawdata(mesh_i, mesh_o, routehandle_mklak_nonorm, lake_i, lake_o, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    do no = 1,ns_o
       if (lake_o(no) < 1.) lake_o(no) = 0.
    enddo

    ! Check global areas
    call output_diagnostics_area(mesh_i, mesh_o, mask_i, frac_o_mklak_nonorm, &
         lake_i, lake_o, "pct lake", percent=.true., ndiag=ndiag, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    deallocate (lake_i)

    ! ----------------------------------------
    ! Wrap things up
    ! ----------------------------------------

    call pio_closefile(pioid_i)

    ! Release memory
    if (mksrf_fdynuse == ' ') then  ! ...else we will reuse it
       deallocate(frac_o_mklak_nonorm)
       call ESMF_RouteHandleDestroy(routehandle_mklak_nonorm, nogarbage = .true., rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort()
    end if
    call ESMF_MeshDestroy(mesh_i, nogarbage = .true., rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort()

    if (root_task) then
       write (ndiag,*) 'Successfully made %lake'
    end if
    call ESMF_VMLogMemInfo("At end of "//trim(subname))

  end subroutine mkpctlak

!===============================================================
  subroutine mklakdep(file_mesh_i, file_data_i, mesh_o, pioid_o, fsurdat, rc)

    ! -------------------
    ! make lake depth
    ! LAKE_DEPTH is written out to fsurdat here
    ! -------------------

    ! uses

    ! input/output variables
    character(len=*)  , intent(in)    :: file_mesh_i      ! input mesh file name
    character(len=*)  , intent(in)    :: file_data_i      ! input data file name
    type(ESMF_Mesh)   , intent(in)    :: mesh_o
    type(file_desc_t) , intent(inout) :: pioid_o
    character(len=*)  , intent(in)    :: fsurdat
    integer           , intent(out)   :: rc

    ! local variables
    type(ESMF_Mesh)        :: mesh_i
    type(ESMF_RouteHandle) :: routehandle
    type(file_desc_t)      :: pioid_i
    integer , allocatable  :: mask_i(:)
    real(r8), allocatable  :: rmask_i(:)
    real(r8), allocatable  :: frac_o(:)
    real(r8), allocatable  :: lakedepth_i(:) ! iput grid: lake depth (m)
    real(r8), allocatable  :: lakedepth_o(:) ! output grid: lake depth (m)
    integer                :: ni,no,k        ! indices
    integer                :: ns_i,ns_o      ! local sizes
    integer                :: ier,rcode      ! error status
    real(r8), parameter    :: min_valid_lakedepth = 0._r8
    character(len=*), parameter :: subname = ' mklakdep '
    !-----------------------------------------------------------------------

    rc = ESMF_SUCCESS

    if (root_task) then
       write(ndiag,*)
       write(ndiag,'(1x,80a1)') ('=',k=1,80)
       write(ndiag,*)
       write(ndiag,'(a)')'Attempting to make lake depth'
       write(ndiag,'(a)') ' Input file is '//trim(file_data_i)
       write(ndiag,'(a)') ' Input mesh file is '//trim(file_mesh_i)
    end if

    ! ----------------------------------------
    ! Read i input data and input mesh and create route handle
    ! ----------------------------------------

    ! Open raw data file
    rcode = pio_openfile(pio_iosystem, pioid_i, pio_iotype, trim(file_data_i), pio_nowrite)

    ! Read in input mesh
    mesh_i = ESMF_MeshCreate(filename=trim(file_mesh_i), fileformat=ESMF_FILEFORMAT_ESMFMESH, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMLogMemInfo("After create mesh_i in "//trim(subname))

    ! Determine ns_i
    call ESMF_MeshGet(mesh_i, numOwnedElements=ns_i, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Determine ns_o
    call ESMF_MeshGet(mesh_o, numOwnedElements=ns_o, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Get the landmask from the input data file and reset the mesh mask based on that
    allocate(rmask_i(ns_i), stat=ier)
    if (ier/=0) call shr_sys_abort()
    allocate(mask_i(ns_i), stat=ier)
    if (ier/=0) call shr_sys_abort()
    call mkpio_get_rawdata(pioid_i, 'LANDMASK', mesh_i, rmask_i, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    do ni = 1,ns_i
       if (rmask_i(ni) > 0._r8) then
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

    ! ----------------------------------------
    ! Create lake parameter (lakedepth)
    ! ----------------------------------------

    if (root_task) then
       write (ndiag,*) 'Attempting to make lake depth .....'
    end if

    ! lakedepth
    allocate(lakedepth_i(ns_i), stat=rcode)
    if (rcode/=0) call shr_sys_abort()
    call mkpio_get_rawdata(pioid_i, 'LAKEDEPTH', mesh_i, lakedepth_i, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! regrid lakedepth_i to lakedepth_o - this also returns lakedepth_i to be used in the global sums below
    allocate (lakedepth_o(ns_o)); lakedepth_o(:) = spval
    call regrid_rawdata(mesh_i, mesh_o, routehandle, lakedepth_i, lakedepth_o, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMLogMemInfo("After regrid_rawdata for lakedepth in "//trim(subname))
    do no = 1,ns_o
       if (frac_o(no) == 0._r8) then
          lakedepth_o(no) = 10._r8
       end if
    enddo
    if (min_bad(lakedepth_o, min_valid_lakedepth, 'lakedepth')) then
       call shr_sys_abort()
    end if

    if (fsurdat /= ' ') then
       if (root_task) write(ndiag, '(a)') trim(subname)//" writing out lakedepth"
       call mkfile_output(pioid_o,  mesh_o,  'LAKEDEPTH', lakedepth_o, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in calling mkfile_output')
       call pio_syncfile(pioid_o)
    end if

    ! Check global areas for lake depth
    call output_diagnostics_continuous(mesh_i, mesh_o, &
         lakedepth_i, lakedepth_o, "lake depth", "m", &
         ndiag=ndiag, rc=rc, mask_i=mask_i, frac_o=frac_o)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! ----------------------------------------
    ! Wrap things up
    ! ----------------------------------------

    call pio_closefile(pioid_i)

    ! Release memory
    deallocate(frac_o)
    call ESMF_RouteHandleDestroy(routehandle, nogarbage = .true., rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort()
    call ESMF_MeshDestroy(mesh_i, nogarbage = .true., rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort()

    if (root_task) then
       write (ndiag,*) 'Successfully made lake depth'
    end if
    call ESMF_VMLogMemInfo("At end of "//trim(subname))

  end subroutine mklakdep

!===============================================================
  subroutine mkwetlnd(file_mesh_i, file_data_i, mesh_o, swmp_o, rc)

    ! ----------------------------------------
    ! Create %wetland
    ! Note PCT_WETLAND is written out of mksurfdata after adjustments are made
    ! ----------------------------------------

    ! input/output variables
    character(len=*)  , intent(in)    :: file_mesh_i      ! input mesh file name
    character(len=*)  , intent(in)    :: file_data_i      ! input data file name
    type(ESMF_Mesh)   , intent(in)    :: mesh_o
    real(r8)          , intent(out)   :: swmp_o(:)        ! output grid: %lake
    integer           , intent(out)   :: rc

    ! local variables
    type(ESMF_RouteHandle) :: routehandle_nonorm
    type(ESMF_Mesh)        :: mesh_i
    type(file_desc_t)      :: pioid_i
    integer , allocatable  :: mask_i(:) 
    real(r8), allocatable  :: rmask_i(:)
    real(r8), allocatable  :: frac_o_nonorm(:)
    real(r8), allocatable  :: swmp_i(:)      ! input grid: percent wetland
    integer                :: ni,no,k        ! indices
    integer                :: ns_i,ns_o      ! local sizes
    integer                :: ier,rcode      ! error status
    character(len=*), parameter :: subname = ' mkwetlnd '
    !-----------------------------------------------------------------------

    rc = ESMF_SUCCESS

    if ( no_inlandwet) then
       if (root_task) then
          write (ndiag,*) 'Attempting to make %wetland .....'
          if (root_task)  write(ndiag, '(a)') trim(subname)//" setting PCT_WETLAND to zero"
       end if
       swmp_o(:) = 0._r8
       RETURN
    end if

    if (root_task) then
       write(ndiag,*) 'Attempting to make %wetland .....'
       write(ndiag,'(a)') ' Input file is '//trim(file_data_i)
       write(ndiag,'(a)') ' Input mesh file is '//trim(file_mesh_i)
    end if

    ! Open input data file
    rcode = pio_openfile(pio_iosystem, pioid_i, pio_iotype, trim(file_data_i), pio_nowrite)

    ! Read in input mesh
    mesh_i = ESMF_MeshCreate(filename=trim(file_mesh_i), fileformat=ESMF_FILEFORMAT_ESMFMESH, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMLogMemInfo("After create mesh_i in "//trim(subname))

    ! Determine ns_i
    call ESMF_MeshGet(mesh_i, numOwnedElements=ns_i, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Determine ns_o
    call ESMF_MeshGet(mesh_o, numOwnedElements=ns_o, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Get the landmask from the input data file and reset the mesh mask based on that
    allocate(rmask_i(ns_i), stat=ier)
    if (ier/=0) call shr_sys_abort()
    allocate(mask_i(ns_i), stat=ier)
    if (ier/=0) call shr_sys_abort()
    call mkpio_get_rawdata(pioid_i, 'LANDMASK', mesh_i, rmask_i, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    do ni = 1,ns_i
       if (rmask_i(ni) > 0._r8) then
          mask_i(ni) = 1
       else
          mask_i(ni) = 0
       end if
    end do
    call ESMF_MeshSet(mesh_i, elementMask=mask_i, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! Create a route handle between the input and output mesh
    allocate(frac_o_nonorm(ns_o))
    ! Note that norm_by_fracs is false in the following because this routehandle is
    ! used to map fields that are expressed in terms of % of the grid cell.
    call create_routehandle_r8(mesh_i=mesh_i, mesh_o=mesh_o, norm_by_fracs=.false., &
         routehandle=routehandle_nonorm, frac_o=frac_o_nonorm, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMLogMemInfo("After create routehandle in "//trim(subname))

    ! read in swmp_i
    allocate(swmp_i(ns_i), stat=rcode)
    if (rcode/=0) call shr_sys_abort()
    call mkpio_get_rawdata(pioid_i, 'PCT_WETLAND', mesh_i, swmp_i, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! regrid swmp_i to swmp_o - this also returns swmp_i to be used in the global sums below
    call regrid_rawdata(mesh_i, mesh_o, routehandle_nonorm, swmp_i, swmp_o, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMLogMemInfo("After regrid_data for wetland")
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    do no = 1,ns_o
       if (swmp_o(no) < 1.) swmp_o(no) = 0.
    enddo

    ! Check global areas
    call output_diagnostics_area(mesh_i, mesh_o, mask_i, frac_o_nonorm, &
         swmp_i, swmp_o, "pct wetland", percent=.true., ndiag=ndiag, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Close the single input data file
    call pio_closefile(pioid_i)

    ! Release memory
    call ESMF_RouteHandleDestroy(routehandle_nonorm, nogarbage = .true., rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort()
    call ESMF_MeshDestroy(mesh_i, nogarbage = .true., rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort()

    if (root_task) then
       write (ndiag,*) 'Successfully made %wetland'
    end if
    call ESMF_VMLogMemInfo("At end of "//trim(subname))

  end subroutine mkwetlnd

!===============================================================
subroutine update_max_array_lake(pct_lakmax_arr,pct_lake_arr)
  !
  ! !DESCRIPTION:
  ! Update the maximum lake percent for landuse.timeseries file
  !
  ! !ARGUMENTS:
  real(r8)         , intent(inout):: pct_lakmax_arr(:)       ! max lake percent
  real(r8)         , intent(in):: pct_lake_arr(:)            ! lake percent that is used to update the old pct_lakmax_arr
  !
  ! !LOCAL VARIABLES:
  integer :: n,ns              ! indices

  character(len=*), parameter :: subname = 'update_max_array_lake'
  !-----------------------------------------------------------------------
  ns = size(pct_lake_arr,1)
  do n = 1, ns
     if (pct_lake_arr(n) > pct_lakmax_arr(n)) then
        pct_lakmax_arr(n) = pct_lake_arr(n)
     end if
  end do

end subroutine update_max_array_lake

end module mklanwatMod
