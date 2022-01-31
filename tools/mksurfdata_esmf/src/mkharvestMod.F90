module mkharvestMod

  !-----------------------------------------------------------------------
  ! Make harvest and grazing data to add to the dynamic PFT file.
  !-----------------------------------------------------------------------

  use ESMF
  use pio
  use shr_kind_mod , only : r8 => shr_kind_r8, CL => shr_kind_CL
  use mkpioMod
  use mkvarpar
  use mkvarctl

  implicit none
  private

  public :: harvestDataType
  integer, private, parameter :: numharv =  9  ! number of harvest and grazing fields

  ! public types
  type :: harvestDataType
     private
     type(ESMF_Routehandle) :: routehandle
     type(ESMF_Mesh)        :: mesh_i
     integer                :: ns_i
     integer                :: ns_o
     integer                :: dims2nd(numharv)    ! 2nd dimension size
     integer                :: CFTdimsize          ! Size of CFT dimension
     integer                :: PFTdimsize          ! Size of PFT dimension
     integer                :: indices1D(numharv)  ! Field indices for CFT variables
     integer                :: indicesCFT(numharv) ! Field indices for CFT variables
     integer                :: indicesPFT(numharv) ! Field indices for PFT variables
     real(r8), pointer      :: data1D(:,:)         ! Input 1D data
     real(r8), pointer      :: data2DCFT(:,:,:)    ! Input 2D data with CFT's
     real(r8), pointer      :: data2DPFT(:,:,:)    ! Input 2D data with PFT's
     real(r8), pointer      :: OutData1D(:,:)      ! Output 1D data
     real(r8), pointer      :: OutData2DCFT(:,:,:) ! Output 2D data with CFT's
     real(r8), pointer      :: OutData2DPFT(:,:,:) ! Output 2D data with natural PFT's
   contains
     procedure :: init           ! Initialization
     procedure :: get1DFieldPtr  ! Get a pointer to a 1D field
     procedure :: get2DFieldPtr  ! Get a pointer to a 2D field
     procedure :: getFieldsIdx   ! Get field indexes to 1D and 2D fields
     procedure :: getFieldsDim   ! Get dimension names for this field
     procedure :: isField1D      ! Return true if field is a 1D field
     procedure :: isField2D      ! Return true if field is a 2D field
     procedure :: num1DFields    ! Return the number of 1D fields
     procedure :: num2DFields    ! Return the number of 2D fields
     procedure :: clean          ! Clean and deallocate everything
  end type harvestDataType

  ! public member functions
  public :: mkharvest_init            ! Initialization
  public :: mkharvest                 ! Calculate the harvest values on output grid
  public :: mkharvest_fieldname       ! Field name for harvest fields on landuse.timeseries
  public :: mkharvest_longname        ! Long name
  public :: mkharvest_units           ! units
  public :: mkharvest_numtypes        ! Number of harvest types
  public :: mkharvest_parse_oride     ! Parse the over-ride string

  ! private member functions: (but public because unit test uses them)
  public mkharvest_fieldInBounds    ! Check that field index is within bounds

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

  character(len=CL), parameter :: string_undef = 'UNSET'
  real(r8),          parameter :: real_undef   = -999.99
  character(len=CL)            :: harvest_longnames(numharv) = string_undef
  character(len=CL)            :: harvest_units(numharv)     = string_undef
  real(r8),  pointer           :: oride_harv(:)        ! array that can override harvesting
  logical                      :: initialized = .false.

!=================================================================================
contains
!=================================================================================

  subroutine init( this, dims2nd, ns_i, ns_o, init_value )
    !
    ! Initialization of the harvestData object
    !
    ! input/output variables
    class(harvestDataType), intent(inout) :: this       ! harvestData object
    integer,                intent(in)    :: dims2nd(:) ! 2nd Dimension sizes
    integer,                intent(in)    :: ns_i       ! Input dimension size
    integer,                intent(in)    :: ns_o       ! Output dimension size
    real(r8),               intent(in)    :: init_value ! Initial value
    !
    ! !local variables:
    character(len=*),  parameter :: subname = 'harvestData::init'
    integer :: num2nd           ! number of non 1D variables
    integer :: numCFT, numPFT   ! number of CFT and PFT variables respectively
    integer :: num1D            ! number of 1D variables
    integer :: n                ! index
    !-----------------------------------------------------------------------

    if ( size(dims2nd) /= numharv )then
       write(*,*) subname//':ERROR:: dims2nd given to init is not the right size'
       call shr_sys_abort()
    end if

    this%CFTdimsize = 64
    this%PFTdimsize = 15
    this%dims2nd = dims2nd
    num2nd = 0
    numCFT = 0
    numPFT = 0
    num1D  = 0
    this%indices1D  = -1
    this%indicesPFT = -1
    this%indicesCFT = -1

    do n = 1, numharv
       if ( dims2nd(n) == 0 )then
          num1D = num1D + 1
          this%indices1D(n) = num1D
       else
          num2nd = num2nd + 1
          if (dims2nd(n) == this%CFTdimsize) then
             numCFT = numCFT + 1
             this%indicesCFT(n) = numCFT
          else if ( dims2nd(n) == this%PFTdimsize )then
             numPFT = numPFT + 1
             this%indicesPFT(n) = numPFT
          else
             write(*,*) 'ERROR:: dims2nd is not the right size (should be 0, 15, or 64) = ', dims2nd(n)
             call shr_sys_abort()
          end if
       end if
    end do

    allocate( this%data1D(ns_i,num1D) )
    allocate( this%OutData1D(ns_o,num1D) )

    this%OutData1D(:,:) = init_value

    if ( num2nd > 0 ) then
       allocate( this%data2DCFT   (ns_i,this%CFTdimsize,numCFT) )
       allocate( this%OutData2DCFT(ns_o,this%CFTdimsize,numCFT) )

       this%OutData2DCFT(:,:,:) = init_value

       allocate( this%data2DPFT   (ns_i,this%PFTdimsize,numPFT) )
       allocate( this%OutData2DPFT(ns_o,this%PFTdimsize,numPFT) )

       this%OutData2DPFT(:,:,:) = init_value
    end if

  end subroutine init

  !=================================================================================
  subroutine mkharvest_init(file_mesh_i, file_data_i, mesh_o, init_val, harvdata, pioid, constant)
    !
    ! Initialization of mkharvest module.
    !
    ! input/output variables:
    real(r8)              , intent(in)    :: init_val ! initial value to set to
    character(len=*)      , intent(in)    :: fharvest ! input harvest dataset file name
    type(file_desc_t)     , intent(out)   :: pioid
    type(harvestDataType) , intent(inout) :: harvdata ! Harvest data
    logical, optional     , intent(in)    :: constant ! Flag if variables are CONST_ version for surface dataset

    ! local variables:
    character(len=CL) :: lunits           ! local units read in
    integer           :: ifld             ! indices
    integer           :: ret              ! return code
    logical           :: lconstant        ! local version of constant flag
    logical           :: varexists        ! If variable exists on file
    integer           :: dim_lengths(3)   ! Dimension lengths on file
    integer           :: dims2nd(numharv) ! Dimension lengths of 3rd dimension for each variable on file
    integer           :: ndims            ! Number of dimensions on file
    integer           :: rcode
    character(len=*),  parameter :: subname = 'mkharvest_init'
    !-----------------------------------------------------------------------

    ! TODO: check that number of global elements in mesh is identical to number of global elements in input data

    lconstant = .false.
    if ( present(constant) ) lconstant = constant
    initialized = .true.

    ! Open fharvest
    call ESMF_VMLogMemInfo("Before pio_openfile for "//trim(file_data_i))
    rcode = pio_openfile(pio_iosystem, pioid, pio_iotype, trim(fharvest), pio_nowrite)
    call ESMF_VMLogMemInfo("After pio_openfile "//trim(file_data_i))

    ! Read in input mesh
    mesh_i = ESMF_MeshCreate(filename=trim(file_mesh_i), fileformat=ESMF_FILEFORMAT_ESMFMESH, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Determine ns_i 
    call ESMF_MeshGet(mesh_i, numOwnedElements=ns_i, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Determine ns_o and allocate data_o
    call ESMF_MeshGet(mesh_o, numOwnedElements=ns_o, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Create a route handle between the input and output mesh
    call create_routehandle_r8(mesh_i, mesh_o, routehandle, rc)
    call ESMF_VMLogMemInfo("After create routehandle in "//trim(subname))

    do ifld = 1, numharv
       call pio_seterrorhandling(pioid, PIO_BCAST_ERROR)
       ! Determine if variable is on input dataset (fharvest)
       varname = mkharvest_fieldname(ifld, constant=lconstant)
       rCode = pio_inq_varid(varname, pio_varid)
       call pio_seterrorhandling(pioid, PIO_INTERNAL_ERROR)
       if (rcode == PIO_NOERR) then
          varexists = .true.
       else
          varexists = .false.
       end if
       if (varexists) then
          ! Get variable attributes if they exist
          call pio_seterrorhandling(pioid, PIO_BCAST_ERROR)
          rcode = pio_inq_attname(pioid, pio_varid, attnum, attname)
          call pio_seterrorhandling(pioid, PIO_INTERNAL_ERROR)
          if (rcode == 'PIO_NOERR') then
             rcode = pio_get_att_text( pioid, pio_varid, 'long_name', harvest_longnames(ifld))
             rcode = pio_get_att_text( pioid, pio_varid, 'units',     harvest_units(ifld))
          end if

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
       else
          if (root_task) then
             write(ndiag,'(a)') "SKIPPING: "//mkharvest_fieldname(ifld, constant=lconstant)
          end if
          harvest_longnames(ifld) = trim(mkharvest_fieldname(ifld, constant=lconstant)) // " (zeroed out)"
          harvest_units(ifld) = "not_read_in"
       end if
    end do
    call pio_closefile(pioid)

    ! Initialize harvest datatype
    call harvdata%init( dims2nd, ns_i, ns_o, init_val, routehandle, mesh_i )

    allocate( oride_harv(numharv) )
    oride_harv(:) = real_undef

  end subroutine mkharvest_init

  !=================================================================================
  subroutine mkharvest(file_data_i, mesh_o, harvdata, rc)
    !
    ! Make harvest data for the dynamic PFT dataset.
    ! This dataset consists of the normalized harvest or grazing fraction (0-1) of
    ! the model.
    !
    ! input/output variables:
    character(len=*)      , intent(in)    :: file_data_i ! input data file name
    type(ESMF_Mesh)       , intent(in)    :: mesh_o      ! model mesh
    type(harvestDataType) , intent(inout) :: harvdata    ! Harvest data
    integer               , intent(out)   :: rc          ! return code

    ! local variables:
    type(file_desc_t)           :: pioid
    type(var_desc_t)            :: pio_varid
    integer                     :: ni,no,ns_i,ns_o       ! indices
    integer                     :: k,l,n,m               ! indices
    integer                     :: rcode, ier             ! error status
    integer                     :: ndims
    integer , allocatable       :: dimlengths(:)
    real(r8), allocatable       :: data_i(:,:)
    real(r8), allocatable       :: data_o(:,:)
    real(r8), allocatable       :: frac_i(:)
    real(r8), allocatable       :: frac_o(:)
    real(r8)                    :: gharv_o(numharv)      ! output grid: global area harvesting
    real(r8)                    :: garea_o               ! output grid: global area
    real(r8)                    :: gharv_i(numharv)      ! input grid: global area harvesting
    real(r8)                    :: garea_i               ! input grid: global area
    integer                     :: ifld                  ! indices
    logical                     :: varexists             ! If variable exists or not
    integer, allocatable        :: ind1D(:)              ! Index of 1D harvest fields
    integer, allocatable        :: ind2D(:)              ! Index of 2D harvest fields
    real(r8), pointer           :: data1D_i(:)           ! 1D input data
    real(r8), pointer           :: data2D_i(:,:)         ! 2D output data
    real(r8), pointer           :: data1D_o(:)           ! 1D output data
    real(r8), pointer           :: data2D_o(:,:)         ! 2D output data
    real(r8), allocatable       :: read_data2d_i(:,:)    ! input 2d data read in
    real(r8), allocatable       :: read_data2d_o(:,:)    ! regridded input 2d data
    character(len=CS)           :: varname
    integer                     :: nlev                  ! inner dimension of input 2d data read in
    character(len=*), parameter :: unit = '10**6 km**2'  ! Output units
    real(r8), parameter         :: fac = 1.e-06_r8       ! Output factor
    real(r8), parameter         :: rat = fac/100._r8     ! Output factor divided by 100%
    character(len=*), parameter :: subname = 'mkharvest'
    !-----------------------------------------------------------------------

    if (root_task) then
       write (ndiag,'(a)') 'Attempting to make harvest fields .....'
    end if

    ! -----------------------------------------------------------------
    ! Normally read in the harvesting file, and then regrid to output grid
    ! -----------------------------------------------------------------

    call harvdata%getFieldsIdx( ind1D, ind2D )

    if ( all(oride_harv == real_undef ) )then

       ! -----------------------------------------------------------------
       ! Read input harvesting file
       ! -----------------------------------------------------------------
       
       ! Determine frac_o (regrid frac_i to frac_o)
       allocate(frac_i(ns_i))
       allocate(frac_o(ns_o))
       call mkpio_get_rawdata(pioid, 'LANDMASK', mesh_i, frac_i, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call regrid_rawdata(mesh_i, mesh_o, routehandle, frac_i, frac_o, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_VMLogMemInfo("After regrid landmask in  "//trim(subname))

       ! Read in input 1d fields if they exists and map to output grid
       do k = 1, harvdata%num1Dfields()
          ifld = ind1D(k)
          data1d_i => harvdata%get1DFieldPtr( ifld )
          data1d_o => harvdata%get1DFieldPtr( ifld, output=.true. )
          data1d_i(:) = 0.0_r8
          call pio_seterrorhandling(pioid, PIO_BCAST_ERROR)
          varname = mkharvest_fieldname(ifld)
          rCode = pio_inq_varid(varname, pio_varid)
          call pio_seterrorhandling(pioid, PIO_INTERNAL_ERROR)
          if (rcode == PIO_NOERR) then
             call mkpio_get_rawdata(pioid, varname, mesh_i, data1d_i, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
             call regrid_rawdata(mesh_i, mesh_o, routehandle, data1d_i, data1d_o, rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
             do n = 1,ns_o
                if (frac_o(n) > 0._r8) then
                   data1d_o(n) = datda1d_o(n) / frac_o(n)
                else
                   data1d_o(n) = 0._r8
                end if
             end do
          else
             write(*,*) "SKIP: "//mkharvest_fieldname(ifld)
          end if
       end do

       ! Read in input 1d fields if they exists and map to output grid
       do k = 1, harvdata%num2Dfields()
          ifld = ind2D(k)
          data2d_i => harvdata%get2DFieldPtr( ifld )
          data2d_o => harvdata%get2DFieldPtr( ifld, output=.true. )
          data2d_i(:,:) = 0.0_r8
          call pio_seterrorhandling(pioid, PIO_BCAST_ERROR)
          varname = mkharvest_fieldname(ifld)
          rCode = pio_inq_varid(varname, pio_varid)
          call pio_seterrorhandling(pioid, PIO_INTERNAL_ERROR)
          if (rcode == PIO_NOERR) then
             nlev = size(data2d_i, dim=2)
             allocate(read_data2d_i(nlev, ns_i))
             allocate(read_data2d_o(nlev, ns_o))
             call mkpio_get_rawdata(pioid, varname, mesh_i, read_data2d_i, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
             call regrid_rawdata(mesh_i, mesh_o, routehandle, read_data2d_i, read_data2d_o, rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
             do n = 1,ns_o
                if (frac_o(n) > 0._r8) then
                   read_data2d_o(:,n) = datda2d_o(:,n) / frac_o(n)
                else
                   read_data2d_o(:,n) = 0._r8
                end if
             end do
             do l = 1,nlev
                do n = 1,ns_o
                   data2d_o(n,l) = read_data2d_o(l,n)
                end do
             end do
             deallocate(read_data2d_i)
             deallocate(read_data2d_o)
          else
             write(*,*) "SKIP: "//mkharvest_fieldname(ifld)
          end if
       end do

       call pio_closefile(pioid)
       
       ! -----------------------------------------------------------------
       ! Error check
       ! Compare global areas on input and output grids
       ! -----------------------------------------------------------------

       ! gharv_i(:) = 0.
       ! garea_i = 0.
       ! do ni = 1, ns_i
       !    garea_i = garea_i + area_src(ni)*re**2
       !    do k = 1, harvdata%num1Dfields()
       !       m = ind1D(k)
       !       data1D_i => harvdata%get1DFieldPtr( m )
       !       gharv_i(m) = gharv_i(m) + data1D_i(ni)*area_src(ni) * mask(ni)*re**2
       !    end do
       ! end do
       ! gharv_o(:) = 0.
       ! garea_o = 0.
       ! do no = 1,ns_o
       !    garea_o = garea_o + area_dst(no)*re**2
       !    do k = 1, harvdata%num1Dfields()
       !       m = ind1D(k)
       !       data1D_o => harvdata%get1DFieldPtr( m, output=.true. )
       !       gharv_o(m) = gharv_o(m) + data1D_o(no)*area_dst(no)* frac_dst(no)*re**2
       !    end do
       ! end do

       ! Write out to diagnostic output file
!        write (ndiag,*)
!        write (ndiag,'(1x,70a1)') ('=',k=1,70)
!        write (ndiag,*) 'Harvesting Output'
!        write (ndiag,'(1x,70a1)') ('=',k=1,70)

!        write (ndiag,*)
!        write (ndiag,'(1x,70a1)') ('.',k=1,70)
!        write (ndiag,1001) unit, unit
! 1001   format (1x,'harvest type   ',20x,' input grid area',' output grid area',/ &
!             1x,33x,'     ',A,'      ',A)
!        write (ndiag,'(1x,70a1)') ('.',k=1,70)
!        write (ndiag,*)
!        do k = 1, harvdata%num1Dfields()
!           m = ind1D(k)
!           write (ndiag,1002) mkharvest_fieldname(m), gharv_i(m)*rat,gharv_o(m)*rat
!        end do
! 1002   format (1x,a35,f16.3,f17.3)

    else

       ! -----------------------------------------------------------------
       ! Otherwise override the harvesting with the input harvest values
       ! -----------------------------------------------------------------

       if ( any(oride_harv == real_undef ) )then
          write(6,*) subname, ' error some override harvesting fields set and others are not = ', oride_harv
          call shr_sys_abort()
       end if
       do k = 1, harvdata%num1Dfields()
          m = ind1D(k)
          if ( oride_harv(m) < 0.0_r8 .or. oride_harv(m) > 100.0_r8 )then
             write(6,*) subname, ' error override harvesting field out of range', &
                  oride_harv(m), ' field = ', mkharvest_fieldname(m)
             call shr_sys_abort()
          end if
       end do
       do no = 1,ns_o
          do k = 1, harvdata%num1Dfields()
             m = ind1D(k)
             data1D_o => harvdata%get1DFieldPtr( m, output=.true. )
             data1D_o(no) = oride_harv(m)
          end do
       end do

    end if

    deallocate( ind1D, ind2D )
    write (6,*) 'Successfully made harvest and grazing'
    write (6,*)

  end subroutine mkharvest

  !=================================================================================
  subroutine mkharvest_parse_oride( string )
    !
    ! Parse the string with harvest and grazing information on it, to override
    ! the file with this information rather than reading from a file.
    !
    use shr_string_mod, only: shr_string_betweenTags

    ! input/output variables:
    character(len=256), intent(in) :: string  ! String to parse with harvest and grazing data

    ! local variables:
    integer :: rc                         ! error return code
    character(len=256) :: substring       ! substring between tags
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
  function get1DFieldPtr( this, nfield, output ) result(ptr1D)
    !
    ! Returns 2D pointer to field data for this index
    !
    ! input/output variables
    class(harvestDataType) , intent(in) :: this     ! harvestData object
    integer                , intent(in) :: nfield   ! field index
    real(r8)               , pointer    :: ptr1D(:) ! Pointer to 1D data
    logical, optional      , intent(in) :: output   ! Flag if this is the output pointer or not (input)
    !
    ! !local variables:
    logical :: loutput    ! Local output flag
    character(len=*),  parameter :: subname = 'harvestData::get1DFieldPtr'
    !-----------------------------------------------------------------------

    loutput = .false.
    if (present(output)) loutput = output

    if ( mkharvest_fieldInBounds( nfield ) .and. this%isField1D(nfield) )then
       if ( .not. loutput ) then
          ptr1D => this%data1D(:,this%indices1D(nfield))
       else
          ptr1D => this%OutData1D(:,this%indices1D(nfield))
       end if
    else
       call shr_sys_abort()
    end if
  end function get1DFieldPtr

  !=================================================================================
  function get2DFieldPtr( this, nfield, output ) result(ptr2D)
    !
    ! Returns 2D pointer to field data for this index
    !
    ! input/output variables:
    class(harvestDataType) , intent(in) :: this       ! harvestData object
    integer                , intent(in) :: nfield     ! field index
    real(r8)               , pointer    :: ptr2D(:,:) ! Pointer to 2D data
    logical, optional      , intent(in) :: output     ! Flag if this is the output pointer or not (input)
    !
    ! local variables:
    logical :: loutput    ! Local output flag
    character(len=*),  parameter :: subname = 'harvestData::get2DFieldPtr'
    !-----------------------------------------------------------------------

    loutput = .false.
    if ( present(output) ) loutput = output
    if ( mkharvest_fieldInBounds( nfield ) .and. this%isField2D(nfield) )then
       if ( .not. loutput ) then
          if ( this%dims2nd(nfield) == this%CFTdimsize )then
             ptr2D => this%data2DCFT(:,:,this%indicesCFT(nfield))
          else
             ptr2D => this%data2DPFT(:,:,this%indicesPFT(nfield))
          end if
       else
          if ( this%dims2nd(nfield) == this%CFTdimsize )then
             ptr2D => this%OutData2DCFT(:,:,this%indicesCFT(nfield))
          else
             ptr2D => this%OutData2DPFT(:,:,this%indicesPFT(nfield))
          end if
       end if
    else
       call shr_sys_abort()
    end if
  end function get2DFieldPtr

  !=================================================================================
  subroutine getFieldsIdx( this, indices1D, indices2D )
    !
    ! Returns list of 1D and 2D fields indices
    !
    ! input/output variables:
    class(harvestDataType), intent(in) :: this  ! harvestData object
    integer, allocatable :: indices1D(:)        ! List of 1D indices
    integer, allocatable :: indices2D(:)        ! List of 2D indices
    !
    ! local variables:
    character(len=*),  parameter :: subname = 'harvestData::getFieldsIdx'
    integer :: ifld, n1, n2  ! field index and field index
    !-----------------------------------------------------------------------

    allocate( indices1D(max(1,this%num1DFields()) ) )
    allocate( indices2D(max(1,this%num2DFields()) ) )
    indices1D = -1
    indices2D = -1
    n1 = 0
    n2 = 0
    do ifld = 1, mkharvest_numtypes()
       if ( this%isField1D(ifld) )then
          n1 = n1 + 1
          indices1D(n1) = ifld
       else if ( this%isField2D(ifld) )then
          n2 = n2 + 1
          indices2D(n2) = ifld
       end if
    end do
  end subroutine getFieldsIdx

  !=================================================================================
  function getFieldsDim( this, nfield ) result(dimname)
    !
    ! Returns list of 1D and 2D fields indices
    !
    ! input/output variables:
    class(harvestDataType), intent(in) :: this  ! harvestData object
    integer, intent(in) :: nfield               ! field index
    character(len=10)   :: dimname              ! Dimension names
    !
    ! local variables:
    character(len=*),  parameter :: subname = 'harvestData::getFieldsDim'
    !-----------------------------------------------------------------------

    if ( this%dims2nd(nfield) == this%CFTdimsize )then
       dimname = "cft"
    else if ( this%dims2nd(nfield) == this%PFTdimsize )then
       dimname = "natpft"
    else
       dimname = "none"
    end if
  end function  getFieldsDim

  !=================================================================================
  logical function isField1D( this, nfield )
    !
    ! Returns true if this field index is a 1D field
    !
    ! input/output variables:
    class(harvestDataType), intent(in) :: this  ! harvestData object
    integer, intent(in) :: nfield               ! field index
    !
    ! local variables:
    character(len=*),  parameter :: subname = 'harvestData::isField1D'
    !-----------------------------------------------------------------------

    isField1D = .false.
    if ( mkharvest_fieldInBounds( nfield ) )then
       if ( this%dims2nd(nfield) == 0 ) isField1D = .true.
    else
       call shr_sys_abort()
    end if
  end function isField1D

  !=================================================================================
  logical function isField2D( this, nfield )
    !
    ! Returns true if this field index is a 2D field
    !
    ! input/output variables:
    class(harvestDataType), intent(in) :: this  ! harvestData object
    integer, intent(in) :: nfield              ! field index
    !
    ! local variables:
    character(len=*),  parameter :: subname = 'harvestData::isField2D'
    !-----------------------------------------------------------------------

    isField2D = .false.
    if ( mkharvest_fieldInBounds( nfield ) )then
       if ( this%dims2nd(nfield) /= 0 ) isField2D = .true.
    else
       call shr_sys_abort()
    end if
  end function isField2D

  !=================================================================================
  integer function num1DFields( this )
    !
    ! Returns the number of 1D fields
    !
    ! input/output variables:
    class(harvestDataType), intent(in) :: this  ! harvestData object
    !
    ! local variables:
    character(len=*),  parameter :: subname = 'harvestData::num1DFields'
    !-----------------------------------------------------------------------
    num1DFields = count( this%dims2nd == 0)
  end function num1DFields

  !=================================================================================
  integer function num2DFields( this )
    !
    ! Returns the number of 2D fields
    !
    ! input/output variables:
    class(harvestDataType), intent(in) :: this  ! harvestData object
    !
    ! local variables:
    character(len=*),  parameter :: subname = 'harvestData::num2DFields'
    !-----------------------------------------------------------------------

    num2DFields = count( this%dims2nd /= 0)
  end function num2DFields


  !=================================================================================
  logical function  mkharvest_fieldInBounds( nfield )
    !
    ! Return true if field index is in bounds and initialization done
    !
    ! input/output variables:
    integer, intent(in) :: nfield    ! field index
    !
    ! local variables:
    character(len=*), parameter :: subname = 'mkharvest_fieldInBounds'
    !-----------------------------------------------------------------------

    if (      nfield < 1       )then
       write(6,*) subname, ' ERROR nfield < 1'
       mkharvest_fieldInBounds = .false.
    else if ( nfield > numharv )then
       write(6,*) subname, ' ERROR nfield > max fields'
       mkharvest_fieldInBounds = .false.
    else if ( .not. initialized ) then
       write(6,*) subname, ' ERROR mkharvest NOT initialized yet!'
       mkharvest_fieldInBounds = .false.
    else
       mkharvest_fieldInBounds = .true.
    end if

  end function mkharvest_fieldInBounds

  !=================================================================================
  character(len=harlen) function mkharvest_fieldname( nfield, constant )
    !
    ! Return harvest fieldname of input field number.
    !
    ! input/output variables:
    integer, intent(in) :: nfield
    logical, intent(in), optional :: constant
    !
    ! local variables:
    character(len=*), parameter :: subname = 'mkharvest_fieldname'
    logical :: lconstant     ! local version of constant flag
    !-----------------------------------------------------------------------
    lconstant = .false.
    if ( present(constant) ) lconstant = constant

    if ( mkharvest_fieldInBounds( nfield ) )then
       if ( .not. lconstant )then
          mkharvest_fieldname = harvest_fieldnames(nfield)
       else
          mkharvest_fieldname = harvest_const_fieldnames(nfield)
       end if
    else
       call shr_sys_abort()
    end if

  end function mkharvest_fieldname

  !=================================================================================
  character(len=CL) function mkharvest_units( nfield )
    !
    !              Return units description of harvest fields
    !
    ! input/output variables:
    integer, intent(in) :: nfield

    ! local variables:
    character(len=*), parameter :: subname = 'mkharvest_units'
    !-----------------------------------------------------------------------

    if ( mkharvest_fieldInBounds( nfield ) )then
       mkharvest_units    = harvest_units(nfield)
    else
       call shr_sys_abort()
    end if

  end function mkharvest_units

  !=================================================================================
  character(len=CL) function mkharvest_longname( nfield )
    !
    !              Return longname description of given input field number.

    ! input/output variables:
    integer, intent(in) :: nfield
    !
    ! local variables:
    character(len=*), parameter :: subname = 'mkharvest_longname'
    !-----------------------------------------------------------------------

    if ( mkharvest_fieldInBounds( nfield ) )then
       mkharvest_longname = harvest_longnames(nfield)
    else
       call shr_sys_abort()
    end if

  end function mkharvest_longname

  !=================================================================================
  integer function mkharvest_numtypes( )
    !
    !              Return number of different harvest field types.
    !
    ! input/output variables:
    character(len=*), parameter :: subname = 'mkharvest_numtypes'

    ! local variables:
    !-----------------------------------------------------------------------
    mkharvest_numtypes = numharv

  end function mkharvest_numtypes

  !=================================================================================
  subroutine clean( this )
    !
    !   Clean and deallocate the harvestData object
    !
    ! input/output variables:
    class(harvestDataType), intent(inout) :: this       ! harvestData object
    !
    ! local variables:
    character(len=*),  parameter :: subname = 'harvestData::clean'
    !-----------------------------------------------------------------------

    this%CFTdimsize = -1
    this%PFTdimsize = -1

    if ( associated(this%data1D)    ) deallocate( this%data1D )
    if ( associated(this%Outdata1D) ) deallocate( this%OutData1D )

    if ( associated(this%data2DCFT)   ) deallocate( this%data2DCFT    )
    if ( associated(this%OutData2DCFT)) deallocate( this%OutData2DCFT )
    if ( associated(this%data2DPFT )  ) deallocate( this%data2DPFT    )
    if ( associated(this%OutData2DPFT)) deallocate( this%OutData2DPFT )
    this%data2DCFT    => null()
    this%OutData2DCFT => null()
    this%data2DPFT    => null()
    this%OutData2DPFT => null()
  end subroutine clean

end module mkharvestMod
