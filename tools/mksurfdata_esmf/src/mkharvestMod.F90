module mkharvestMod

  !-----------------------------------------------------------------------
  ! Make harvest and grazing data to add to the dynamic PFT file.
  !-----------------------------------------------------------------------

  use ESMF
  use pio
  use shr_kind_mod      , only : r8 => shr_kind_r8, r4=>shr_kind_r4, cs => shr_kind_cs, cl => shr_kind_cl
  use shr_sys_mod       , only : shr_sys_abort
  use mkpioMod          , only : pio_iotype, pio_ioformat, pio_iosystem
  use mkpioMod          , only : mkpio_get_rawdata, mkpio_get_rawdata_level, mkpio_get_dimlengths
  use mkpioMod          , only : mkpio_def_spatial_var, mkpio_iodesc_rawdata
  use mkpioMod          , only : mkpio_put_time_slice, mkpio_iodesc_output
  use mkfileMod         , only : mkfile_output
  use mkesmfMod         , only : regrid_rawdata, create_routehandle_r8, get_meshareas
  use mkutilsMod        , only : chkerr
  use mkvarctl          , only : root_task, ndiag, mpicom
  use mkdiagnosticsMod  , only : output_diagnostics_area

  implicit none
  private

  ! public member functions
  public :: mkharvest                 ! Calculate the harvest values on output grid

  ! private data members:
  integer, parameter :: numharv =  9  ! number of harvest and grazing fields
  integer, parameter :: harlen  = 25  ! length of strings for harvest fieldnames
  character(len=harlen), parameter  :: harvest_fieldnames(numharv) = (/ &
       'HARVEST_VH1            ',  &
       'HARVEST_VH2            ',  &
       'HARVEST_SH1            ',  &
       'HARVEST_SH2            ',  &
       'HARVEST_SH3            ',  &
       'GRAZING                ',  &
       'FERTNITRO_CFT          ',  &
       'UNREPRESENTED_PFT_LULCC',  &
       'UNREPRESENTED_CFT_LULCC'   &
       /)
  character(len=harlen), parameter  :: harvest_const_fieldnames(numharv) = (/ &
       'CONST_HARVEST_VH1      ',  &
       'CONST_HARVEST_VH2      ',  &
       'CONST_HARVEST_SH1      ',  &
       'CONST_HARVEST_SH2      ',  &
       'CONST_HARVEST_SH3      ',  &
       'CONST_GRAZING          ',  &
       'CONST_FERTNITRO_CFT    ',  &
       'UNREPRESENTED_PFT_LULCC',  &
       'UNREPRESENTED_CFT_LULCC'   &
       /)

  type(ESMF_Mesh)        :: mesh_i
  type(ESMF_RouteHandle) :: routehandle_r8
  real(r8), allocatable  :: frac_o(:)

  character(len=*) , parameter :: u_FILE_u = &
       __FILE__

!=================================================================================
contains
!=================================================================================

  subroutine mkharvest(file_mesh_i, file_data_i, mesh_o, pioid_o, ntime, rc)
    !
    ! Make harvest data for the dynamic PFT dataset.
    ! This dataset consists of the normalized harvest or grazing fraction (0-1) of
    ! the model.
    !
    ! input/output variables:
    character(len=*)      , intent(in)    :: file_mesh_i
    character(len=*)      , intent(in)    :: file_data_i ! input data file name
    type(ESMF_Mesh)       , intent(in)    :: mesh_o      ! model mesh
    type(file_desc_t)     , intent(inout) :: pioid_o
    integer, optional     , intent(in)    :: ntime
    integer               , intent(out)   :: rc          ! return code

    ! local variables:
    type(file_desc_t)      :: pioid_i
    type(var_desc_t)       :: pio_varid_i
    type(var_desc_t)       :: pio_varid_o
    type(io_desc_t)        :: pio_iodesc_i
    type(io_desc_t)        :: pio_iodesc_o
    integer                :: ns_i, ns_o         ! input/output sizes
    integer                :: ni,no              ! indices
    integer                :: k,l,m              ! indices
    integer                :: ifld               ! indices
    integer                :: dims2nd            ! variable dimension lengths of 3rd dimension
    character(len=cs)      :: name2nd            ! name of 2nd dimension
    logical                :: varexists          ! true if variable exists on file
    character(len=cs)      :: varname_i          ! input variable name
    character(len=cs)      :: varname_o          ! output variable name
    real(r8) , allocatable :: rmask_i(:)         ! real value of input mask (read in)
    integer  , allocatable :: mask_i(:)          ! integer value of rmask_i
    real(r8) , allocatable :: data1d_i(:)        ! input 1d data
    real(r8) , allocatable :: data1d_o(:)        ! otuput 1d data
    real(r8) , allocatable :: data2d_o(:,:)      ! output 2d data
    real(r8) , allocatable :: read_data2d_i(:,:)
    real(r8) , allocatable :: read_data2d_o(:,:)
    real(r8) , allocatable :: area_i(:)          ! areas on input mesh
    real(r8) , allocatable :: area_o(:)          ! areas on output mesh
    integer                :: rcode, ier         ! error status
    character(len=*), parameter :: subname = 'mkharvest'
    !-----------------------------------------------------------------------

    rc = ESMF_SUCCESS

    if (root_task) then
       write(ndiag,*)
       write(ndiag,'(1x,80a1)') ('=',k=1,80)
       write(ndiag,*)
       write(ndiag,'(a)') 'Attempting to make harvest fields .....'
       write(ndiag,'(a)') ' Input file is '//trim(file_data_i)
       write(ndiag,'(a)') ' Input mesh file is '//trim(file_mesh_i)
    end if

    ! Determine ns_o
    call ESMF_MeshGet(mesh_o, numOwnedElements=ns_o, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Normally read in the harvesting file, and then regrid to output grid
    rcode = pio_openfile(pio_iosystem, pioid_i, pio_iotype, trim(file_data_i), pio_nowrite)

    ! Read in input mesh
    if (.not. ESMF_MeshIsCreated(mesh_i)) then
       mesh_i = ESMF_MeshCreate(filename=trim(file_mesh_i), fileformat=ESMF_FILEFORMAT_ESMFMESH, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    ! Determine ns_i
    call ESMF_MeshGet(mesh_i, numOwnedElements=ns_i, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

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

    ! Create a route handle between the input and output mesh
    ! NOTE: this must be done after mask_i is set in mesh_i
    if (.not. ESMF_RouteHandleIsCreated(routehandle_r8)) then
       allocate(frac_o(ns_o))
       call create_routehandle_r8(mesh_i=mesh_i, mesh_o=mesh_o, norm_by_fracs=.true., &
            routehandle=routehandle_r8, frac_o=frac_o, rc=rc)
       call ESMF_VMLogMemInfo("After create routehandle in "//trim(subname))
    end if

    ! Read in input 1d fields if they exists and map to output grid
    do ifld = 1,numharv
       varname_i = trim(harvest_fieldnames(ifld))
       if (.not. present(ntime)) then  ! not transient, i.e. constant
          varname_o = trim(harvest_const_fieldnames(ifld))
       else
          varname_o = varname_i
       end if
       call mkharvest_check_input_var(pioid_i, trim(varname_i), varexists, dims2nd, name2nd)
       if (varexists) then
          if (dims2nd == 0) then

             ! 1d output
             allocate(data1d_i(ns_i))
             allocate(data1d_o(ns_o))

             ! read in input 1d variable
             call mkpio_get_rawdata(pioid_i, varname_i, mesh_i, data1d_i, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return

             ! regrid input variable
             call regrid_rawdata(mesh_i, mesh_o, routehandle_r8, data1d_i, data1d_o, rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return

             ! write out mapped variable
             if (present(ntime)) then
                if (root_task)  write(ndiag, '(a,i8)') subname//" writing out 1d "//trim(varname_o)//' at time ',ntime
                rcode = pio_inq_varid(pioid_o, trim(varname_o), pio_varid_o)
                call mkpio_iodesc_output(pioid_o, mesh_o, trim(varname_o), pio_iodesc_o, rc)
                if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in making an iodesc for '//trim(varname_o))
                call pio_setframe(pioid_o, pio_varid_o, int(ntime, kind=Pio_Offset_Kind))
                call pio_write_darray(pioid_o, pio_varid_o, pio_iodesc_o, data1d_o, rcode)
                call pio_freedecomp(pioid_o, pio_iodesc_o)
             else
                if (root_task)  write(ndiag, '(a)') subname//" writing out 1d "//trim(varname_o)
                call mkfile_output(pioid_o, mesh_o, trim(varname_o), data1d_o, rc=rc)
                if (ChkErr(rc,__LINE__,u_FILE_u)) return
             end if
             call pio_syncfile(pioid_o)

             ! Compare global areas on input and output grids for 1d variables
             call output_diagnostics_area(mesh_i, mesh_o, mask_i, frac_o, &
                  data1d_i*0.01_r8, data1d_o*0.01_r8, trim(varname_o), percent=.false., ndiag=ndiag, rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
                  
             deallocate(data1d_i)
             deallocate(data1d_o)

          else ! 2d output

             ! Read in input data 
             allocate(read_data2d_i(dims2nd, ns_i))
             allocate(read_data2d_o(dims2nd, ns_o))
             call mkpio_get_rawdata(pioid_i, trim(varname_i), mesh_i, read_data2d_i, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return

             ! Regrid input input data
             call regrid_rawdata(mesh_i, mesh_o, routehandle_r8, read_data2d_i, read_data2d_o, 1, dims2nd, rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
             deallocate(read_data2d_i)

             ! Fill in output 2d array
             allocate(data2d_o(ns_o,dims2nd))
             do l = 1,dims2nd 
                do no = 1,ns_o
                   data2d_o(no,l) = read_data2d_o(l,no)
                end do
             end do
             deallocate(read_data2d_o)

             ! write out variable
             if (present(ntime)) then
                if (root_task)  write(ndiag, '(a,i8)') subname//" writing out 2d "//trim(varname_o)//' at time ',ntime
                rcode = pio_inq_varid(pioid_o, trim(varname_o), pio_varid_o)
                call mkpio_iodesc_output(pioid_o, mesh_o, trim(varname_o), pio_iodesc_o, rc)
                if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in making an iodesc for '//trim(varname_o))
                call pio_setframe(pioid_o, pio_varid_o, int(ntime, kind=Pio_Offset_Kind))
                call pio_write_darray(pioid_o, pio_varid_o, pio_iodesc_o, data2d_o, rcode)
                call pio_freedecomp(pioid_o, pio_iodesc_o)
             else
                if (root_task)  write(ndiag, '(a)') subname//" writing out 2d "//trim(varname_o)
                call mkfile_output(pioid_o, mesh_o, trim(varname_o), data2d_o, lev1name=trim(name2nd), rc=rc)
                if (ChkErr(rc,__LINE__,u_FILE_u)) return
             end if
             call pio_syncfile(pioid_o)
             deallocate(data2d_o)

          end if
       end if
    end do

    if (.not. present(ntime)) then  ! ...else we will reuse it
       deallocate(frac_o)
       call ESMF_RouteHandleDestroy(routehandle_r8, nogarbage = .true., rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort()
       call ESMF_MeshDestroy(mesh_i, nogarbage = .true., rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort()
    end if

    if (root_task) then
       write (ndiag,'(a)') 'Successfully made harvest and grazing'
    end if

  end subroutine mkharvest

  !=================================================================================
  subroutine mkharvest_check_input_var(pioid, varname, varexists, dims2nd, name2nd)

    ! input/output variables
    type(file_desc_t) , intent(inout) :: pioid
    character(len=*)  , intent(in)    :: varname
    logical           , intent(out)   :: varexists
    integer           , intent(out)   :: dims2nd
    character(len=*)  , intent(out)   :: name2nd

    ! local variable
    type(var_desc_t) :: pio_varid
    integer          :: ndims          ! number of variable dimension
    integer          :: dimids(3)      ! variable dimension dim ids
    integer          :: dimlens(3)     ! variable dimensions sizes
    integer          :: ifld           ! index
    integer          :: rcode          ! error status
    character(len=*)  , parameter :: subname = 'mkharvest_check_input_var'
    !-----------------------------------------------------------------------

    dims2nd = -999
    name2nd = 'unset'

    call pio_seterrorhandling(pioid, PIO_BCAST_ERROR)
    rcode = pio_inq_varid(pioid, varname, pio_varid)
    call pio_seterrorhandling(pioid, PIO_INTERNAL_ERROR)
    if (rcode == PIO_NOERR) then
       varexists = .true.
    else
       varexists = .false.
    end if
    if (varexists) then
       rcode = pio_inq_varndims(pioid, pio_varid, ndims)
       if ( ndims == 2 )then
          dims2nd = 0
       else if (ndims == 3 )then
          rcode = pio_inq_vardimid(pioid, pio_varid, dimids)
          rcode = pio_inq_dimlen(pioid, dimids(3), dimlens(3))
          dims2nd = dimlens(3)
          rcode = pio_inq_vardimid(pioid, pio_varid, dimids)
          rcode = pio_inq_dimname(pioid, dimids(3), name2nd)
       else
          write(*,*) 'ERROR:: bad dimensionality for variable = ',trim(varname)
          call shr_sys_abort()
       end if
       if (root_task) then
          write(ndiag,'(a)') " reading: "//trim(varname)
       end if
    else
       if (root_task) then
          write(ndiag,'(a)') " skipping: "//trim(varname)
       end if
    end if

  end subroutine mkharvest_check_input_var

end module mkharvestMod
