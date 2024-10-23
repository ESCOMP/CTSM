module mkvocefMod

  !-----------------------------------------------------------------------
  ! Make VOC percentage emissions for surface dataset
  !-----------------------------------------------------------------------

  use ESMF
  use pio
  use shr_kind_mod , only : r8 => shr_kind_r8, r4 => shr_kind_r4
  use shr_sys_mod  , only : shr_sys_abort
  use mkpioMod     , only : mkpio_get_rawdata
  use mkpioMod     , only : mkpio_iodesc_rawdata, pio_iotype, pio_iosystem
  use mkesmfMod    , only : regrid_rawdata, create_routehandle_r8
  use mkutilsMod   , only : chkerr
  use mkvarctl     , only : root_task, ndiag, mpicom
  use mkfileMod    , only : mkfile_output

  implicit none
  private

  public :: mkvocef  ! Get the percentage emissions for VOC for different land cover types

  character(len=*) , parameter :: u_FILE_u = &
       __FILE__

!=================================================================================
contains
!=================================================================================

  subroutine mkvocef(file_mesh_i, file_data_i, mesh_o, pioid_o, lat_o, rc)
    !
    ! make volatile organic coumpunds (VOC) emission factors.
    !
    ! input/output variables
    character(len=*)  , intent(in)    :: file_mesh_i ! input mesh file name
    character(len=*)  , intent(in)    :: file_data_i ! input data file name
    type(ESMF_Mesh)   , intent(in)    :: mesh_o      ! model mesho
    type(file_desc_t) , intent(inout) :: pioid_o     ! output file descripter
    real(r8)          , intent(in)    :: lat_o(:)    ! output latitudes
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
    real(r8), allocatable  :: ef_btr_i(:) ! input grid: EFs for broadleaf trees
    real(r8), allocatable  :: ef_fet_i(:) ! input grid: EFs for fineleaf evergreen
    real(r8), allocatable  :: ef_fdt_i(:) ! input grid: EFs for fineleaf deciduous
    real(r8), allocatable  :: ef_shr_i(:) ! input grid: EFs for shrubs
    real(r8), allocatable  :: ef_grs_i(:) ! input grid: EFs for grasses
    real(r8), allocatable  :: ef_crp_i(:) ! input grid: EFs for crops
    real(r8), allocatable  :: ef_btr_o(:) ! output grid: EFs for broadleaf trees
    real(r8), allocatable  :: ef_fet_o(:) ! output grid: EFs for fineleaf evergreen
    real(r8), allocatable  :: ef_fdt_o(:) ! output grid: EFs for fineleaf deciduous
    real(r8), allocatable  :: ef_shr_o(:) ! output grid: EFs for shrubs
    real(r8), allocatable  :: ef_grs_o(:) ! output grid: EFs for grasses
    real(r8), allocatable  :: ef_crp_o(:) ! output grid: EFs for crops
    integer                :: ier, rcode  ! error status
    real(r8)               :: relerr = 0.00001_r8 ! max error: sum overlap wts ne 1
    character(len=*), parameter :: subname = 'mkvocef'
    !-----------------------------------------------------------------------

    rc = ESMF_SUCCESS
    call ESMF_VMLogMemInfo("At start of "//trim(subname))

    if (root_task) then
       write(ndiag,*)
       write(ndiag,'(1x,80a1)') ('=',k=1,80)
       write(ndiag,*)
       write(ndiag,'(a)') 'Attempting to make VOC emission factors .....'
       write(ndiag,'(a)') ' Input data file is '//trim(file_data_i)
       write(ndiag,'(a)') ' Input mesh file is '//trim(file_mesh_i)
    end if

    ! Open input data file
    rcode = pio_openfile(pio_iosystem, pioid_i, pio_iotype, trim(file_data_i), pio_nowrite)
    call ESMF_VMLogMemInfo("After pio_openfile "//trim(file_data_i))

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

    ! Allocate output variables
    allocate (ef_btr_o(ns_o)) ; ef_btr_o(:) = 0._r8
    allocate (ef_fet_o(ns_o)) ; ef_fet_o(:) = 0._r8
    allocate (ef_fdt_o(ns_o)) ; ef_fdt_o(:) = 0._r8
    allocate (ef_shr_o(ns_o)) ; ef_shr_o(:) = 0._r8
    allocate (ef_grs_o(ns_o)) ; ef_grs_o(:) = 0._r8
    allocate (ef_crp_o(ns_o)) ; ef_crp_o(:) = 0._r8

    ! Get the landmask from the input data file and reset the mesh mask based on that
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

    ! Read input Emission Factors
    allocate (ef_btr_i(ns_i))
    if (ier/=0) call shr_sys_abort()
    allocate (ef_fet_i(ns_i))
    if (ier/=0) call shr_sys_abort()
    allocate (ef_fdt_i(ns_i))
    if (ier/=0) call shr_sys_abort()
    allocate (ef_shr_i(ns_i))
    if (ier/=0) call shr_sys_abort()
    allocate (ef_grs_i(ns_i))
    if (ier/=0) call shr_sys_abort()
    allocate (ef_crp_i(ns_i))
    if (ier/=0) call shr_sys_abort()

    call mkpio_get_rawdata(pioid_i, 'ef_btr', mesh_i, ef_btr_i, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call mkpio_get_rawdata(pioid_i, 'ef_fet', mesh_i, ef_fet_i, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call mkpio_get_rawdata(pioid_i, 'ef_fdt', mesh_i, ef_fdt_i, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call mkpio_get_rawdata(pioid_i, 'ef_shr', mesh_i, ef_shr_i, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call mkpio_get_rawdata(pioid_i, 'ef_grs', mesh_i, ef_grs_i, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call mkpio_get_rawdata(pioid_i, 'ef_crp', mesh_i, ef_crp_i, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! Regrid input data to model resolution
    call regrid_rawdata(mesh_i, mesh_o, routehandle, ef_btr_i, ef_btr_o, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call regrid_rawdata(mesh_i, mesh_o, routehandle, ef_fet_i, ef_fet_o, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call regrid_rawdata(mesh_i, mesh_o, routehandle, ef_fdt_i, ef_fdt_o, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call regrid_rawdata(mesh_i, mesh_o, routehandle, ef_shr_i, ef_shr_o, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call regrid_rawdata(mesh_i, mesh_o, routehandle, ef_grs_i, ef_grs_o, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call regrid_rawdata(mesh_i, mesh_o, routehandle, ef_crp_i, ef_crp_o, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Check for conservation
    do no = 1, ns_o
       if ( ef_btr_o(no) < 0._r8 ) then
          write (6,*) 'MKVOCEF error: EF btr = ',ef_btr_o(no), ' is negative for no = ',no
          call shr_sys_abort()
       end if
       if ( ef_fet_o(no) < 0._r8 ) then
          write (6,*) 'MKVOCEF error: EF fet = ',ef_fet_o(no), ' is negative for no = ',no
          call shr_sys_abort()
       end if
       if ( ef_fdt_o(no) < 0._r8 ) then
          write (6,*) 'MKVOCEF error: EF fdt = ',ef_fdt_o(no), ' is negative for no = ',no
          call shr_sys_abort()
       end if
       if ( ef_shr_o(no) < 0._r8 ) then
          write (6,*) 'MKVOCEF error: EF shr = ',ef_shr_o(no), ' is negative for no = ',no
          call shr_sys_abort()
       end if
       if ( ef_grs_o(no) < 0._r8 ) then
          write (6,*) 'MKVOCEF error: EF grs = ',ef_grs_o(no), ' is negative for no = ',no
          call shr_sys_abort()
       end if
       if ( ef_crp_o(no) < 0._r8 ) then
          write (6,*) 'MKVOCEF error: EF crp = ',ef_crp_o(no), ' is negative for no = ',no
          call shr_sys_abort()
       end if
    enddo


    ! If have pole points on grid - set south pole to glacier
    ! north pole is assumed as non-land
    do no = 1,ns_o
       if (abs((lat_o(no) - 90._r8)) < 1.e-6_r8) then
          ef_btr_o(no) = 0._r8
          ef_fet_o(no) = 0._r8
          ef_fdt_o(no) = 0._r8
          ef_shr_o(no) = 0._r8
          ef_grs_o(no) = 0._r8
          ef_crp_o(no) = 0._r8
       end if
    end do

    if (root_task)  write(ndiag, '(a)') trim(subname)//" writing out voc emission factors"
    call mkfile_output(pioid_o,  mesh_o,  'EF1_BTR', ef_btr_o, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in calling mkfile_output')
    call mkfile_output(pioid_o,  mesh_o,  'EF1_FET', ef_fet_o, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in calling mkfile_output')
    call mkfile_output(pioid_o,  mesh_o,  'EF1_FDT', ef_fdt_o, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in calling mkfile_output')
    call mkfile_output(pioid_o,  mesh_o,  'EF1_SHR', ef_shr_o, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in calling mkfile_output')
    call mkfile_output(pioid_o,  mesh_o,  'EF1_GRS', ef_grs_o, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in calling mkfile_output')
    call mkfile_output(pioid_o,  mesh_o,  'EF1_CRP', ef_crp_o, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in calling mkfile_output')
    call pio_syncfile(pioid_o)

    ! -----------------------------------------------------------------
    ! Wrap up
    ! -----------------------------------------------------------------

    ! Close the input file
    call pio_closefile(pioid_i)
    call ESMF_VMLogMemInfo("After pio_closefile in "//trim(subname))

    ! Release memory
    call ESMF_RouteHandleDestroy(routehandle, nogarbage = .true., rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort()
    call ESMF_MeshDestroy(mesh_i, nogarbage = .true., rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort()
    call ESMF_VMLogMemInfo("After destroy operations in "//trim(subname))

    if (root_task) then
       write (ndiag,'(a)') 'Successfully made VOC Emission Factors'
    end if

  end subroutine mkvocef

end module mkvocefMod
