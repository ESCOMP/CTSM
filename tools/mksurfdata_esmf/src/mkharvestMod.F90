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
  use mkpioMod          , only : mkpio_def_spatial_var, mkpio_iodesc_output, mkpio_iodesc_rawdata
  use mkpioMod          , only : mkpio_put_time_slice
  use mkesmfMod         , only : regrid_rawdata, create_routehandle_r8, get_meshareas
  use mkutilsMod        , only : chkerr
  use mkvarctl          , only : root_task, ndiag, mpicom

  implicit none
  private

#include <mpif.h>

  integer, private, parameter :: numharv =  9  ! number of harvest and grazing fields

  ! public member functions
  public :: mkharvest                 ! Calculate the harvest values on output grid
  public :: mkharvest_parse_oride     ! Parse the over-ride string

  ! private data members:
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

  logical :: is_field_1d(num_harv)

  character(len=CL), parameter :: string_undef = 'UNSET'
  real(r8),          parameter :: real_undef   = -999.99
  character(len=CL)            :: harvest_longnames(numharv) = string_undef
  character(len=CL)            :: harvest_units(numharv)     = string_undef
  real(r8),  pointer           :: oride_harv(:)        ! array that can override harvesting
  logical                      :: initialized = .false.

  type(ESMF_Mesh)       , private :: mesh_i
  type(ESMF_RouteHandle), private :: routehandle
  real(r8), allocatable , private :: frac_o(:)

  character(len=*) , parameter :: u_FILE_u = &
       __FILE__

!=================================================================================
contains
!=================================================================================
    !-----------------------------------------------------------------------

    rc = ESMF_SUCCESS

    ! TODO: check that number of global elements in mesh is identical
    ! to number of global elements in input data

    if (root_task) then
       write (ndiag,'(a)') 'Attempting to initialize harvest module .... '
    end if

    lconstant = .false.
    if ( present(constant) ) lconstant = constant

    initialized = .true.

    ! Open input harvest file
    rcode = pio_openfile(pio_iosystem, pioid, pio_iotype, trim(file_data_i), pio_nowrite)

    do ifld = 1, numharv
       call pio_seterrorhandling(pioid, PIO_BCAST_ERROR)
       ! Determine if variable is on input dataset (fharvest)
       varname = mkharvest_fieldname(ifld, constant=lconstant)
       rCode = pio_inq_varid(pioid, varname, pio_varid)
       call pio_seterrorhandling(pioid, PIO_INTERNAL_ERROR)
       if (rcode == PIO_NOERR) then
          varexists = .true.
       else
          varexists = .false.
       end if
       if (varexists) then
          ! get dims2nd for variable
          call mkpio_get_dimlengths(pioid, varname, ndims, dim_lengths)
          if ( ndims == 2 )then
             dims2nd(ifld) = 0
          else if ( ndims == 3 )then
             dims2nd(ifld) = dim_lengths(3)
          else
             write(*,*) 'ERROR:: bad dimensionality for variable = ', mkharvest_fieldname(ifld, constant=lconstant)
             call shr_sys_abort()
          end if
          if (root_task) then
             write(ndiag,'(a)') "Will Read: "//mkharvest_fieldname(ifld, constant=lconstant)
          end if
       else
          if (root_task) then
             write(ndiag,'(a)') "Will Skip: "//mkharvest_fieldname(ifld, constant=lconstant)
          end if
       end if
    end do
    call pio_closefile(pioid)

    if (root_task) then
       write (ndiag,'(a)') 'finished mkharvest_init'
    end if

  end subroutine mkharvest_init

  !=================================================================================
  subroutine mkharvest(file_data_i, mesh_o, pioid_o, constant, rc)
    !
    ! Make harvest data for the dynamic PFT dataset.
    ! This dataset consists of the normalized harvest or grazing fraction (0-1) of
    ! the model.
    !
    ! input/output variables:
    character(len=*)      , intent(in)    :: file_data_i ! input data file name
    type(ESMF_Mesh)       , intent(in)    :: mesh_o      ! model mesh
    type(file_desc_t)     , intent(in)    :: pioid_o 
    logical, optional     , intent(in)    :: constant
    integer               , intent(out)   :: rc          ! return code

    ! local variables:
    type(file_desc_t)      :: pioid_i
    type(var_desc_t)       :: pio_varid
    integer                :: ni,no,ns_i,ns_o    ! indices
    integer                :: k,l,n,m            ! indices
    integer                :: ndims
    integer                :: nlev               ! inner dimension of input 2d data read in
    integer                :: ifld               ! indices
    character(len=cs)      :: varname
    logical                :: varexists          ! If variable exists or not
    integer  , allocatable :: mask_i(:)
    real(r8) , allocatable :: frac_i(:)
    real(r8) , allocatable :: frac_o(:)
    real(r8) , allocatable :: data_i(:,:)
    real(r8) , allocatable :: data_o(:,:)
    real(r8) , allocatable :: read_data2d_i(:,:) ! input 2d data read in
    real(r8) , allocatable :: read_data2d_o(:,:) ! regridded input 2d data
    real(r8) , allocatable :: area_i(:)
    real(r8) , allocatable :: area_o(:)
    integer                :: dims2nd(numharv)   ! Dimension lengths of 3rd dimension for each variable on file
    logical                :: lconstant          ! local version of constant flag
    integer                :: dim_lengths(3)     ! Dimension lengths on file
    integer                :: ndims              ! Number of dimensions on file
    logical                :: varexists          ! If variable exists on file
    integer                :: ifld               ! indices
    integer                :: rcode, ier         ! error status
    logical                :: lconstant
    character(len=*), parameter :: subname = 'mkharvest'
    !-----------------------------------------------------------------------

    rc = ESMF_SUCCESS

    ! Normally read in the harvesting file, and then regrid to output grid

    if (root_task) then
       write (ndiag,'(a)') 'Attempting to make harvest fields .....'
    end if

    rcode = pio_openfile(pio_iosystem, pioid, pio_iotype, trim(file_data_i), pio_nowrite)

    ! Read in input mesh
    if (.not. ESMF_MeshCreated(mesh_i)) then
       mesh_i = ESMF_MeshCreate(filename=trim(file_mesh_i), fileformat=ESMF_FILEFORMAT_ESMFMESH, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if
    
    ! Create a route handle between the input and output mesh
    if (.not. RouteHandleisCreated(routehandle_r8)) then
       allocate(frac_o(ns_o))
       call create_routehandle_r8(mesh_i, mesh_o, routehandle, frac_o=frac_o, rc=rc)
       call ESMF_VMLogMemInfo("After create routehandle in "//trim(subname))
    end if

    ! Determine ns_i
    call ESMF_MeshGet(mesh_i, numOwnedElements=ns_i, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    
    ! Determine ns_o
    call ESMF_MeshGet(mesh_o, numOwnedElements=ns_o, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    
    allocate( oride_harv(numharv) )
    oride_harv(:) = real_undef

    if ( all(oride_harv == real_undef ) )then

       ! Get the landmask from the file and reset the mesh mask based on that
       allocate(frac_i(ns_i), stat=ier)
       if (ier/=0) call shr_sys_abort()
       allocate(mask_i(ns_i), stat=ier)
       if (ier/=0) call shr_sys_abort()
       call mkpio_get_rawdata(pioid, 'LANDMASK', mesh_i, frac_i, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       do ni = 1,ns_i
          if (frac_i(ni) > 0._r4) then
             mask_i(ni) = 1
          else
             mask_i(ni) = 0
          end if
       end do
       call ESMF_MeshSet(mesh_i, elementMask=mask_i, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       ! Read in input 1d fields if they exists and map to output grid
       do ifld = 1,numharv
          varname_i = trim(harvest_fieldnames(ifld))
          if (lconstant) then
             varname_o = trim(harvest_const_fieldnames(ifld)
          else
             varname_o = varname_i
          end if
          ! Check if the variable is on the input file
          call pio_seterrorhandling(pioid_i, PIO_BCAST_ERROR)
          rCode = pio_inq_varid(pioid, varname, pio_varid)
          call pio_seterrorhandling(pioid, PIO_INTERNAL_ERROR)
          if (rcode == PIO_NOERR) then
             varexists = .true.
          else
             varexists = .false.
          end if
          if (varexists) then
             call mkpio_get_dimlengths(pioid, varname, ndims, dim_lengths)
             if ( ndims == 2 )then
                dims2nd(ifld) = 0
             else if (ndims == 3 )then
                dims2nd(ifld) = dim_lengths(3)
             else
                write(*,*) 'ERROR:: bad dimensionality for variable = ',trim(varname_i)
                call shr_sys_abort()
             end if
             if (root_task) then
                write(ndiag,'(a)') "Will Read: "//trim(varname_i)
             end if
          else
             if (root_task) then
                write(ndiag,'(a)') "Will Skip: "//trim(varname_i)
             end if
          end if
          if (varexists) then
             if (dim2nd(ifld) == 0) then
                allocate(data1d_i(ns_i))
                allocate(data1d_o(ns_o))
                ! read in input 1d variable
                call mkpio_get_rawdata(pioid, varname, mesh_i, data1d_i, rc=rc)
                if (chkerr(rc,__LINE__,u_FILE_u)) return

                ! regrid input variable
                call regrid_rawdata(mesh_i, mesh_o, routehandle, data1d_i, data1d_o, rc=rc)
                if (ChkErr(rc,__LINE__,u_FILE_u)) return
                deallocate(data1d_i)
                if (root_task)  write(ndiag, '(a)') trim(subname)//" writing out "//trim(varname_o)

                ! write out mapped variable
                call mkfile_output(pioid_o, mesh_o, trim(varname_o), data1d_o, rc=rc)
                if (ChkErr(rc,__LINE__,u_FILE_u)) return
                deallocate(data1d_o))
             else
                nlev = dims2nd(ifld)
                allocate(read_data2d_i(nlev, ns_i))
                allocate(read_data2d_o(nlev, ns_o))
                ! read in input variable
                call mkpio_get_rawdata(pioid, varname, mesh_i, read_data2d_i, rc=rc)
                if (chkerr(rc,__LINE__,u_FILE_u)) return

                ! regrid input variable
                call regrid_rawdata(mesh_i, mesh_o, routehandle, read_data2d_i, read_data2d_o, 1, nlev, rc)
                if (ChkErr(rc,__LINE__,u_FILE_u)) return
                deallocate(read_data2d_i)

                ! write out variable
                allocate(data2d_o(ns_o, nlev))
                do l = 1,nlev
                   do n = 1,ns_o
                      data2d_o(n,l) = read_data2d_o(l,n)
                   end do
                end do
                call mkfile_output(pioid_o, mesh_o, trim(varname_o), data2d_o, rc=rc)
                if (ChkErr(rc,__LINE__,u_FILE_u)) return
                deallocate(data2d_o)
             end if
          end if
       end do

       ! Compare global areas on input and output grids
       ! call check_global_sums('harvest type', ns_i, ns_o, mesh_i, mesh_o, mask_i, frac_o, IND1d, rc)
       ! if (ChkErr(rc,__LINE__,u_FILE_u)) return

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
    lconstant = .false.
    if ( present(constant) ) lconstant = constant
    if (lconstant) then
       deallocate(frac_o)
       call ESMF_RouteHandleDestroy(routehandle, nogarbage = .true., rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort()
       call ESMF_MeshDestroy(mesh_i, nogarbage = .true., rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort()
    end if

    if (root_task) then
       write (ndiag,'(a)') 'Successfully made harvest and grazing'
       write (ndiag,*)
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
  subroutine check_global_sums(this, name, ns_i, ns_o, mesh_i, mesh_o, mask_i, frac_o, ind1D, rc)

    ! input/otuput variables
    class(harvestDataType), intent(inout) :: this       ! harvestData object
    character(len=*)      , intent(in)    :: name
    integer               , intent(in)    :: ns_i
    integer               , intent(in)    :: ns_o
    type(ESMF_Mesh)       , intent(in)    :: mesh_i
    type(ESMF_Mesh)       , intent(in)    :: mesh_o
    integer               , intent(in)    :: mask_i(:)
    real(r8)              , intent(in)    :: frac_o(:)
    integer               , intent(in)    :: ind1D(:)
    integer               , intent(out)   :: rc

    ! local variables
    integer               :: ni, no, k, m
    integer               :: ier
    real(r8), allocatable :: area_i(:)
    real(r8), allocatable :: area_o(:)
    real(r8)              :: local_i(numharv)  ! input  grid: global area harvesting
    real(r8)              :: local_o(numharv)  ! output grid: global area harvesting
    real(r8)              :: global_i(numharv) ! input  grid: global area harvesting
    real(r8)              :: global_o(numharv) ! output grid: global area harvesting
    real(r8), pointer     :: data1D_i(:)            ! 1D input data
    real(r8), pointer     :: data1D_o(:)            ! 1D output data
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
       do k = 1, this%num1Dfields()
          m = ind1D(k)
          data1D_i => this%get1DFieldPtr( m )
          local_i(m) = local_i(m) + data1D_i(ni)  *area_i(ni) * mask_i(ni)
       end do
    end do
    call mpi_reduce(local_i, global_i, numharv, MPI_REAL8, MPI_SUM, 0, mpicom, ier)

    ! Output grid global area
    local_o(:) = 0.
    do no = 1,ns_o
       do k = 1, this%num1Dfields()
          m = ind1D(k)
          data1D_o => this%get1DFieldPtr( m, output=.true. )
          local_o(m) = local_o(m) + data1D_o(no) * area_o(no) * frac_o(no)
       end do
    end do
    call mpi_reduce(local_o, global_o, numharv, MPI_REAL8, MPI_SUM, 0, mpicom, ier)

    ! Comparison
    write (ndiag,*)
    write (ndiag,*)
    write (ndiag,'(1x,70a1)') ('.',k=1,70)
    write (ndiag,101) unit, unit
101 format (1x,'harvest type   ',20x,' input grid area',' output grid area',/ &
            1x,33x,'     ',A,'      ',A)
    write (ndiag,'(1x,70a1)') ('.',k=1,70)
    write (ndiag,*)
    do k = 1, this%num1Dfields()
       m = ind1D(k)
       write (ndiag,102) mkharvest_fieldname(m), global_i(m)*rat, global_o(m)*rat
    end do
102 format (1x,a35,f16.3,f17.3)

  end subroutine check_global_sums

end module mkharvestMod
