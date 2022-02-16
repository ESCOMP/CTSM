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

  implicit none
  private

#include <mpif.h>

  ! public member functions
  public :: mkharvest                 ! Calculate the harvest values on output grid
  public :: mkharvest_parse_oride     ! Parse the over-ride string

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

  character(len=CL), parameter :: string_undef = 'UNSET'
  real(r8),          parameter :: real_undef   = -999.99
  real(r8),  pointer           :: oride_harv(:)        ! array that can override harvesting

  type(ESMF_Mesh)        :: mesh_i
  type(ESMF_RouteHandle) :: routehandle_r8
  real(r8), allocatable  :: frac_o(:)
  logical                :: initialized = .false.

  character(len=*) , parameter :: u_FILE_u = &
       __FILE__

!=================================================================================
contains
!=================================================================================

  subroutine mkharvest(file_mesh_i, file_data_i, mesh_o, pioid_o, all_veg, constant, ntime, rc)
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
    logical               , intent(in)    :: all_veg
    logical               , intent(in)    :: constant
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

    ! If all veg then write out output data with values set to harvest_initval and return
    if (all_veg) then
       varname_i = trim(harvest_fieldnames(ifld))
       if (constant) then
          varname_o = trim(harvest_const_fieldnames(ifld))
       else
          varname_o = varname_i
       end if
       call mkharvest_check_input_var(pioid_i, trim(varname_i), varexists, dims2nd, name2nd)
       if (varexists) then
          if (dims2nd == 0) then
             allocate(data1d_o(ns_o))
             data1d_o(:) = 0._r8
             call mkfile_output(pioid_o, mesh_o, trim(varname_o), data1d_o, rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
             deallocate(data1d_o)
          else
             allocate(data2d_o(ns_o,dims2nd))
             data2d_o(:,:) = 0._r8
             call mkfile_output(pioid_o, mesh_o, trim(varname_o), data2d_o, lev1name=name2nd, rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
             deallocate(data2d_o)
          end if
       end if
       call pio_closefile(pioid_i)
       RETURN
    end if

    ! Read in input mesh
    if (.not. ESMF_MeshIsCreated(mesh_i)) then
       mesh_i = ESMF_MeshCreate(filename=trim(file_mesh_i), fileformat=ESMF_FILEFORMAT_ESMFMESH, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    ! Determine ns_i
    call ESMF_MeshGet(mesh_i, numOwnedElements=ns_i, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! TODO: imlement oride_harv funcitionality - this is hard-wired for now
    allocate( oride_harv(numharv) )
    oride_harv(:) = real_undef

    if ( all(oride_harv == real_undef ) )then

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
          call create_routehandle_r8(mesh_i, mesh_o, routehandle_r8, frac_o=frac_o, rc=rc)
          call ESMF_VMLogMemInfo("After create routehandle in "//trim(subname))
       end if

       ! Read in input 1d fields if they exists and map to output grid
       do ifld = 1,numharv
          varname_i = trim(harvest_fieldnames(ifld))
          if (constant) then
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
                   rcode = pio_inq_varid(pioid_o, trim(varname_o), pio_varid_o)
                   call mkpio_iodesc_output(pioid_o, mesh_o, trim(varname_o), pio_iodesc_o, rc)
                   if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in making an iodesc for '//trim(varname_o))
                   call pio_setframe(pioid_o, pio_varid_o, int(ntime, kind=Pio_Offset_Kind))
                   call pio_write_darray(pioid_o, pio_varid_o, pio_iodesc_o, data1d_o, rcode)
                   call pio_freedecomp(pioid_o, pio_iodesc_o)
                else
                   if (root_task)  write(ndiag, '(a)') " writing out 1d "//trim(varname_o)
                   call mkfile_output(pioid_o, mesh_o, trim(varname_o), data1d_o, rc=rc)
                   if (ChkErr(rc,__LINE__,u_FILE_u)) return
                   call pio_syncfile(pioid_o)
                end if

                ! TODO: uncomment the following and validate
                ! Compare global areas on input and output grids for 1d variables
                ! call mkharvest check_global_sums('harvest type '//trim(varname_o), ns_i, ns_o, &
                !      mesh_i, mesh_o, mask_i, frac_o, data1d_i, data1d_o, rc)
                ! if (ChkErr(rc,__LINE__,u_FILE_u)) return

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
                   rcode = pio_inq_varid(pioid_o, trim(varname_o), pio_varid_o)
                   call mkpio_iodesc_output(pioid_o, mesh_o, trim(varname_o), pio_iodesc_o, rc)
                   if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in making an iodesc for '//trim(varname_o))
                   call pio_setframe(pioid_o, pio_varid_o, int(ntime, kind=Pio_Offset_Kind))
                   call pio_write_darray(pioid_o, pio_varid_o, pio_iodesc_o, data2d_o, rcode)
                   call pio_freedecomp(pioid_o, pio_iodesc_o)
                else
                   if (root_task)  write(ndiag, '(a)') " writing out 2d "//trim(varname_o)
                   call mkfile_output(pioid_o, mesh_o, trim(varname_o), data2d_o, lev1name=trim(name2nd), rc=rc)
                   if (ChkErr(rc,__LINE__,u_FILE_u)) return
                end if
                deallocate(data2d_o)

             end if
          end if

       end do

    else

       ! Otherwise override the harvesting with the input harvest values
       ! TODO: implement this
       ! if ( any(oride_harv == real_undef ) )then
       !    write(6,*) subname, ' error some override harvesting fields set and others are not = ', oride_harv
       !    call shr_sys_abort()
       ! end if
       ! do k = 1, harvdata%num1Dfields()
       !    m = ind1D(k)
       !    if ( oride_harv(m) < 0.0_r8 .or. oride_harv(m) > 100.0_r8 )then
       !       write(6,*) subname, ' error override harvesting field out of range', &
       !            oride_harv(m), ' field = ', mkharvest_fieldname(m)
       !       call shr_sys_abort()
       !    end if
       ! end do
       ! do no = 1,ns_o
       !    do k = 1, harvdata%num1Dfields()
       !       m = ind1D(k)
       !       data1D_o(no) = oride_harv(m)
       !    end do
       ! end do

    end if

    ! If constant model, clean up the mapping
    if (constant) then
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
  subroutine mkharvest_parse_oride( string )
    !
    ! Parse the string with harvest and grazing information on it, to override
    ! the file with this information rather than reading from a file.
    !
    use shr_string_mod, only: shr_string_betweenTags

    ! input/output variables:
    character(len=CS), intent(in) :: string  ! String to parse with harvest and grazing data

    ! local variables:
    character(len=CS) :: substring       ! substring between tags
    integer           :: rc              ! error return code
    character(len=*), parameter :: harv_start = "<harv>"
    character(len=*), parameter :: harv_end   = "</harv>"
    character(len=*), parameter :: graz_start = "<graz>"
    character(len=*), parameter :: graz_end   = "</graz>"
    character(len=*), parameter :: subname = 'mkharvest_parse_oride'
    !-----------------------------------------------------------------------

    call shr_string_betweenTags( string, harv_start, harv_end, substring, rc )
    if ( rc /= 0 )then
       write(6,*) subname//'Trouble finding harvest start end tags'
       call shr_sys_abort()
    end if
    read(substring,*) oride_harv(1:numharv-1)
    call shr_string_betweenTags( string, graz_start, graz_end, substring, rc )
    if ( rc /= 0 )then
       write(6,*) subname//'Trouble finding grazing start end tags'
       call shr_sys_abort()
    end if
    read(substring,*) oride_harv(numharv)
    if ( harvest_fieldnames(numharv) /= 'GRAZING' )then
       write(6,*) subname, ' grazing is NOT last field as was expected'
       call shr_sys_abort()
    end if

  end subroutine mkharvest_parse_oride

  !=================================================================================
  subroutine mkharvest_check_global_sums(name, ns_i, ns_o, mesh_i, mesh_o, mask_i, frac_o, &
       data_i, data_o, rc)

    ! input/otuput variables
    character(len=*) , intent(in)  :: name
    integer          , intent(in)  :: ns_i
    integer          , intent(in)  :: ns_o
    type(ESMF_Mesh)  , intent(in)  :: mesh_i
    type(ESMF_Mesh)  , intent(in)  :: mesh_o
    integer          , intent(in)  :: mask_i(:)
    real(r8)         , intent(in)  :: frac_o(:)
    real(r8)         , intent(in)  :: data_i(:) ! 1D input data
    real(r8)         , intent(in)  :: data_o(:) ! 1D output data
    integer          , intent(out) :: rc

    ! local variables
    integer               :: ni, no, k, m
    integer               :: ier
    real(r8), allocatable :: area_i(:)
    real(r8), allocatable :: area_o(:)
    real(r8)              :: local_i(numharv)  ! input  grid: global area harvesting
    real(r8)              :: local_o(numharv)  ! output grid: global area harvesting
    real(r8)              :: global_i(numharv) ! input  grid: global area harvesting
    real(r8)              :: global_o(numharv) ! output grid: global area harvesting
    real(r8), parameter   :: fac = 1.e-06_r8   ! Output factor
    real(r8), parameter   :: rat = fac/100._r8 ! Output factor divided by 100%
    character(len=*) , parameter :: unit = '10**6 km**2' ! Output units
    !-----------------------------------------------------------------------

    rc = ESMF_SUCCESS

    allocate(area_i(ns_i), stat=ier)
    if (ier/=0) call shr_sys_abort()
    allocate(area_o(ns_o), stat=ier)
    if (ier/=0) call shr_sys_abort()
    call get_meshareas(mesh_i, area_i, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call get_meshareas(mesh_o, area_o, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! Input grid global area
    local_i(:) = 0.
    do ni = 1, ns_i
       local_i(m) = local_i(m) + data_i(ni)  *area_i(ni) * mask_i(ni)
    end do
    call mpi_reduce(local_i, global_i, numharv, MPI_REAL8, MPI_SUM, 0, mpicom, ier)

    ! Output grid global area
    local_o(:) = 0.
    do no = 1,ns_o
       local_o(m) = local_o(m) + data_o(no) * area_o(no) * frac_o(no)
    end do
    call mpi_reduce(local_o, global_o, numharv, MPI_REAL8, MPI_SUM, 0, mpicom, ier)

    ! Comparison
    write (ndiag,*)
    write (ndiag,*)
    write (ndiag,'(1x,70a1)') ('.',k=1,70)
    write (ndiag,101) unit, unit
101 format (1x,'harvest type' ,20x,' input grid area',' output grid area',/ &
            1x,33x,'     ',A,'      ',A)
    write (ndiag,'(1x,70a1)') ('.',k=1,70)
    write (ndiag,*)
    write (ndiag,102) trim(name), global_i*rat, global_o*rat
102 format (1x,a35,f16.3,f17.3)

  end subroutine mkharvest_check_global_sums

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
