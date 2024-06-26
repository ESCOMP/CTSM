module mkVICparamsMod

  !-----------------------------------------------------------------------
  ! make parameters for VIC
  !-----------------------------------------------------------------------

  use ESMF
  use pio
  use shr_kind_mod     , only : r8 => shr_kind_r8, r4 => shr_kind_r4, cl => shr_kind_cl
  use shr_sys_mod      , only : shr_sys_abort
  use mkpioMod         , only : mkpio_get_rawdata, pio_iotype, pio_iosystem
  use mkesmfMod        , only : regrid_rawdata, create_routehandle_r8
  use mkutilsMod       , only : chkerr
  use mkvarctl         , only : root_task, ndiag, spval
  use mkchecksMod      , only : min_bad
  use mkfileMod        , only : mkfile_output
  use mkdiagnosticsMod , only : output_diagnostics_continuous

  implicit none
  private

  public :: mkVICparams            ! make VIC parameters

  character(len=*) , parameter :: u_FILE_u = &
       __FILE__

!===============================================================
contains
!===============================================================

  subroutine mkVICparams(file_mesh_i, file_data_i, mesh_o, pioid_o, rc)
    !
    ! make VIC parameters
    !
    ! input/output variables
    character(len=*) , intent(in)    :: file_mesh_i ! input mesh file name
    character(len=*) , intent(in)    :: file_data_i ! input data file name
    type(ESMF_Mesh)  , intent(in)    :: mesh_o      ! model mesho
    type(file_desc_t), intent(inout) :: pioid_o     ! output file descripter
    integer          , intent(out)   :: rc

    ! local variables:
    type(ESMF_RouteHandle) :: routehandle
    type(ESMF_Mesh)        :: mesh_i
    type(file_desc_t)      :: pioid_i
    integer                :: ni,no
    integer                :: ns_i, ns_o
    integer                :: n,l,k
    real(r8), allocatable  :: binfl_o(:)  ! output VIC b parameter for the Variable Infiltration Capacity Curve (unitless)
    real(r8), allocatable  :: ws_o(:)     ! output VIC Ws parameter for the ARNO curve (unitless)
    real(r8), allocatable  :: dsmax_o(:)  ! output VIC Dsmax parameter for the ARNO curve (mm/day)
    real(r8), allocatable  :: ds_o(:)     ! output VIC Ds parameter for the ARNO curve (unitless)
    integer , allocatable  :: mask_i(:) 
    real(r8), allocatable  :: rmask_i(:)
    real(r8), allocatable  :: frac_o(:)
    real(r8), allocatable  :: data_i(:) ! data on input grid
    real(r8), parameter    :: min_valid_binfl = 0._r8
    real(r8), parameter    :: min_valid_ws    = 0._r8
    real(r8), parameter    :: min_valid_dsmax = 0._r8
    real(r8), parameter    :: min_valid_ds    = 0._r8
    integer                :: ier,rcode ! error status
    character(len=cl)      :: errmsg
    character(len=*), parameter :: subname = 'mkVICparams'
    !-----------------------------------------------------------------------

    rc = ESMF_SUCCESS

    if (root_task) then
       write(ndiag,*)
       write(ndiag,'(1x,80a1)') ('=',k=1,80)
       write(ndiag,*)
       write(ndiag,'(a)')'Attempting to make VIC parameters.....'
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

    ! Determine ns_o
    call ESMF_MeshGet(mesh_o, numOwnedElements=ns_o, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Allocate output variables
    allocate (binfl_o(ns_o)) ; binfl_o(:) = spval
    allocate (ws_o(ns_o))    ; ws_o(:)    = spval
    allocate (dsmax_o(ns_o)) ; dsmax_o(:) = spval
    allocate (ds_o(ns_o))    ; ds_o(:)    = spval

    ! Get the landmask from the file and reset the mesh mask based on that
    allocate(rmask_i(ns_i), stat=ier)
    if (ier/=0) call shr_sys_abort()
    allocate(mask_i(ns_i), stat=ier)
    if (ier/=0) call shr_sys_abort()
    call mkpio_get_rawdata(pioid_i, 'mask', mesh_i, rmask_i, rc=rc)
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
    do n = 1, ns_o
       if ((frac_o(n) < 0.0) .or. (frac_o(n) > 1.0001)) then
          write(errmsg,'(a,f13.5,2x,i4)') "ERROR:: frac_o out of range: ", frac_o(n),n
          call shr_sys_abort(trim(errmsg),u_FILE_u,__LINE__)
       end if
    end do

    ! -----------------------------------------------------------------
    ! Determine binfl
    ! -----------------------------------------------------------------

    ! Read in input data_i for a variety of inputs
    allocate(data_i(ns_i), stat=ier)
    if (ier/=0) call shr_sys_abort('allocation error for binfl_i',u_FILE_u,__LINE__)

    ! Read in binfl_i into data_i
    call mkpio_get_rawdata(pioid_i, 'binfl', mesh_i, data_i, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMLogMemInfo("After mkpio_getrawdata in "//trim(subname))

    ! Regrid binfl_i to binfl_o and check validity of output data
    call regrid_rawdata(mesh_i, mesh_o, routehandle, data_i, binfl_o, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    do no = 1,ns_o
       if (frac_o(no) == 0._r8) then
          binfl_o(no) = 0.1_r8
       end if
    end do
    if (min_bad(binfl_o, min_valid_binfl, 'binfl')) then
       call shr_sys_abort('error for min_bad',u_FILE_u,__LINE__)
    end if

    ! Calculate global diagnostics for binfl
    call output_diagnostics_continuous(mesh_i, mesh_o, data_i, binfl_o, &
       "VIC b parameter", "unitless", ndiag=ndiag, rc=rc, mask_i=mask_i, frac_o=frac_o)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! -----------------------------------------------------------------
    ! Determine Ws
    ! -----------------------------------------------------------------

    ! Read in Ws into data_i
    call mkpio_get_rawdata(pioid_i, 'Ws', mesh_i, data_i, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMLogMemInfo("After mkpio_getrawdata in "//trim(subname))

    ! Regrid Ws_i to Ws_o and check validity of output data
    call regrid_rawdata(mesh_i, mesh_o, routehandle, data_i, ws_o, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    do no = 1,ns_o
       if (frac_o(no) == 0._r8) then
          ws_o(no) = 0.75_r8
       end if
    end do
    if (min_bad(ws_o, min_valid_ws, 'Ws')) then
       call shr_sys_abort()
    end if

    ! Calculate global diagnostics for Ws
    call output_diagnostics_continuous(mesh_i, mesh_o, data_i, ws_o, &
       "VIC Ws parameter", "unitless", ndiag=ndiag, rc=rc, mask_i=mask_i, frac_o=frac_o)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! -----------------------------------------------------------------
    ! Determine DsMax
    ! -----------------------------------------------------------------

    ! Read in Dsmax into data_i
    call mkpio_get_rawdata(pioid_i, 'Dsmax', mesh_i, data_i, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMLogMemInfo("After mkpio_getrawdata in "//trim(subname))

    ! Regrid Dsmax_i to Dsmax_o and check validity of output data
    call regrid_rawdata(mesh_i, mesh_o, routehandle, data_i, dsmax_o, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    do no = 1,ns_o
       if (frac_o(no) == 0._r8) then
          dsmax_o(no) = 10._r8
       end if
    end do
    if (min_bad(dsmax_o, min_valid_dsmax, 'Dsmax')) then
       call shr_sys_abort()
    end if

    ! Calculate global diagnostics for Dsmax
    call output_diagnostics_continuous(mesh_i, mesh_o, data_i, dsmax_o, &
       "VIC Dsmax parameter", "mm/day", ndiag=ndiag, rc=rc, mask_i=mask_i, frac_o=frac_o)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! -----------------------------------------------------------------
    ! Regrid Ds
    ! -----------------------------------------------------------------

    ! Read in Ds into data_i
    call mkpio_get_rawdata(pioid_i, 'Ds', mesh_i, data_i, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMLogMemInfo("After mkpio_getrawdata in "//trim(subname))

    ! Regrid Ds_i to Ds_o and check validity of output data
    call regrid_rawdata(mesh_i, mesh_o, routehandle, data_i, ds_o, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    do no = 1,ns_o
       if (frac_o(no) == 0._r8) then
          ds_o(no) = 0.1_r8
       end if
    end do
    if (min_bad(ds_o, min_valid_ds, 'Ds')) then
       call shr_sys_abort()
    end if

    ! Calculate global diagnostics for Ws
    call output_diagnostics_continuous(mesh_i, mesh_o, data_i, ds_o, &
       "VIC Ds parameter", "unitless", ndiag=ndiag, rc=rc, mask_i=mask_i, frac_o=frac_o)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! -----------------------------------------------------------------
    ! Write output
    ! -----------------------------------------------------------------
    
    if (root_task)  write(ndiag, '(a)') trim(subname)//" writing out VIC parameters"
    call mkfile_output(pioid_o,  mesh_o,  'binfl', binfl_o, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in calling mkfile_output for binfl')
    call mkfile_output(pioid_o,  mesh_o,  'Ws', ws_o, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in calling mkfile_output for Ws')
    call mkfile_output(pioid_o,  mesh_o,  'Dsmax', dsmax_o, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in calling mkfile_output for Dsmax')
    call mkfile_output(pioid_o,  mesh_o,  'Ds', ds_o, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in calling mkfile_output for Ds')
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
       write (ndiag,'(a)') 'Successfully made VIC parameters'
    end if

  end subroutine mkVICparams

end module mkVICparamsMod
