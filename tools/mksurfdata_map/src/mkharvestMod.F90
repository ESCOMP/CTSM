module mkharvestMod
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: mkharvest
!
! !DESCRIPTION:
! Make harvest and grazing data to add to the dynamic PFT file.
!
! !REVISION HISTORY:
! Author: Erik Kluzek
!
!-----------------------------------------------------------------------
! !USES:
  use shr_kind_mod , only : r8 => shr_kind_r8, CL => shr_kind_CL
  use shr_sys_mod  , only : shr_sys_flush
  use mkdomainMod  , only : domain_checksame

  implicit none

  private

! !PUBLIC DATA MEMBERS:

  public :: harvestDataType
  integer, private, parameter :: numharv =  9  ! number of harvest and grazing fields

  type :: harvestDataType
     private
     real(r8), pointer :: data1D(:,:)            ! Input 1D data
     real(r8), pointer :: data2DCFT(:,:,:)       ! Input 2D data with CFT's
     real(r8), pointer :: data2DPFT(:,:,:)       ! Input 2D data with PFT's
     real(r8), pointer :: OutData1D(:,:)         ! Output 1D data
     real(r8), pointer :: OutData2DCFT(:,:,:)    ! Output 2D data with CFT's
     real(r8), pointer :: OutData2DPFT(:,:,:)    ! Output 2D data with natural PFT's
     integer              :: dims2nd(numharv)    ! 2nd dimension size
     integer              :: CFTdimsize          ! Size of CFT dimension
     integer              :: PFTdimsize          ! Size of PFT dimension
     integer              :: indices1D(numharv)  ! Field indices for CFT variables
     integer              :: indicesCFT(numharv) ! Field indices for CFT variables
     integer              :: indicesPFT(numharv) ! Field indices for PFT variables
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

! !PUBLIC MEMBER FUNCTIONS:
  public mkharvest_init            ! Initialization
  public mkharvest                 ! Calculate the harvest values on output grid
  public mkharvest_fieldname       ! Field name for harvest fields on landuse.timeseries
  public mkharvest_longname        ! Long name
  public mkharvest_units           ! units
  public mkharvest_numtypes        ! Number of harvest types
  public mkharvest_parse_oride     ! Parse the over-ride string

! !PRIVATE MEMBER FUNCTIONS: (but public because unit test uses them)
  public mkharvest_fieldInBounds    ! Check that field index is within bounds

! !PRIVATE DATA MEMBERS:

  integer, parameter :: harlen  = 25  ! length of strings for harvest fieldnames
  character(len=harlen), parameter  :: harvest_fieldnames(numharv) = (/ &
                                                        'HARVEST_VH1  ',  &
                                                        'HARVEST_VH2  ',  &
                                                        'HARVEST_SH1  ',  &
                                                        'HARVEST_SH2  ',  &
                                                        'HARVEST_SH3  ',  &
                                                        'GRAZING      ',  &
                                                        'FERTNITRO_CFT',  &
                                                        'PFT_LULCC    ',  &
                                                        'CFT_LULCC    '   &
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
  character(len=CL), parameter :: string_undef = 'STRING_UNDEFINED'
  real(r8),          parameter :: real_undef   = -999.99
  character(len=CL), save :: harvest_longnames(numharv) = string_undef
  character(len=CL), save :: harvest_units(numharv)     = string_undef
  real(r8),  pointer     :: oride_harv(:)        ! array that can override harvesting
  logical          , save :: initialized = .false.

!EOP
!-----------------------------------------------------------------------
contains
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: init
!
! !INTERFACE:
  subroutine init( this, dims2nd, ns_i, ns_o, init_value )
!
! !DESCRIPTION:
!   Initialization of the harvestData object
!
! !USES:
    implicit none
!
! !ARGUMENTS:
    class(harvestDataType), intent(INOUT) :: this       ! harvestData object
    integer,                intent(IN)    :: dims2nd(:) ! 2nd Dimension sizes
    integer,                intent(IN)    :: ns_i       ! Input dimension size
    integer,                intent(IN)    :: ns_o       ! Output dimension size
    real(r8),               intent(IN)    :: init_value ! Initial value
!
! !REVISION HISTORY:
! Author: Erik Kluzek
!
! !LOCAL VARIABLES:
    character(len=*),  parameter :: subname = 'harvestData::init'
    integer :: num2nd           ! number of non 1D variables
    integer :: numCFT, numPFT   ! number of CFT and PFT variables respectively
    integer :: num1D            ! number of 1D variables
    integer :: n                ! index
!EOP
!-----------------------------------------------------------------------
    if ( size(dims2nd) /= numharv )then
       write(*,*) subname//':ERROR:: dims2nd given to init is not the right size'
       call abort()
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
          if (      dims2nd(n) == this%CFTdimsize )then
             numCFT = numCFT + 1
             this%indicesCFT(n) = numCFT
          else if ( dims2nd(n) == this%PFTdimsize )then
             numPFT = numPFT + 1
             this%indicesPFT(n) = numPFT
          else
             write(*,*) 'ERROR:: dims2nd is not the right size (should be 0, 15, or 64) = ', dims2nd(n)
             call abort()
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

!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: get1DFieldPtr
!
! !INTERFACE:
  function get1DFieldPtr( this, nfield, output ) result(ptr1D)
!
! !DESCRIPTION:
!   Returns 2D pointer to field data for this index
!
! !USES:
    implicit none
!
! !ARGUMENTS:
    class(harvestDataType), intent(IN) :: this  ! harvestData object
    integer, intent(in) :: nfield               ! field index
    real(r8), pointer :: ptr1D(:)               ! Pointer to 1D data
    logical, optional, intent(in) :: output     ! Flag if this is the output pointer or not (input)
!
! !REVISION HISTORY:
! Author: Erik Kluzek
!
! !LOCAL VARIABLES:
    character(len=*),  parameter :: subname = 'harvestData::get1DFieldPtr'
    logical :: loutput    ! Local output flag
!EOP
!-----------------------------------------------------------------------
    loutput = .false.
    if ( present(output) ) loutput = output
    if ( mkharvest_fieldInBounds( nfield ) .and. this%isField1D(nfield) )then
       if ( .not. loutput ) then
          ptr1D => this%data1D(:,this%indices1D(nfield))
       else
          ptr1D => this%OutData1D(:,this%indices1D(nfield))
       end if
    else
       call abort()
    end if
  end function get1DFieldPtr

!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: get2DFieldPtr
!
! !INTERFACE:
  function get2DFieldPtr( this, nfield, output ) result(ptr2D)
!
! !DESCRIPTION:
!   Returns 2D pointer to field data for this index
!
! !USES:
    implicit none
!
! !ARGUMENTS:
    class(harvestDataType), intent(IN) :: this  ! harvestData object
    integer, intent(in) :: nfield               ! field index
    real(r8), pointer :: ptr2D(:,:)             ! Pointer to 2D data
    logical, optional, intent(in) :: output     ! Flag if this is the output pointer or not (input)
!
! !REVISION HISTORY:
! Author: Erik Kluzek
!
! !LOCAL VARIABLES:
    character(len=*),  parameter :: subname = 'harvestData::get2DFieldPtr'
    logical :: loutput    ! Local output flag
!EOP
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
       call abort()
    end if
  end function get2DFieldPtr

!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: getFieldsIdx
!
! !INTERFACE:
  subroutine getFieldsIdx( this, indices1D, indices2D )
!
! !DESCRIPTION:
!   Returns list of 1D and 2D fields indices
!
! !USES:
    implicit none
!
! !ARGUMENTS:
    class(harvestDataType), intent(IN) :: this  ! harvestData object
    integer, allocatable :: indices1D(:)        ! List of 1D indices
    integer, allocatable :: indices2D(:)        ! List of 2D indices
!
! !REVISION HISTORY:
! Author: Erik Kluzek
!
! !LOCAL VARIABLES:
    character(len=*),  parameter :: subname = 'harvestData::getFieldsIdx'
    integer :: ifld, n1, n2  ! field index and field index
!EOP
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

!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: getFieldsDim
!
! !INTERFACE:
  function getFieldsDim( this, nfield ) result(dimname)
!
! !DESCRIPTION:
!   Returns list of 1D and 2D fields indices
!
! !USES:
    implicit none
!
! !ARGUMENTS:
    class(harvestDataType), intent(IN) :: this  ! harvestData object
    integer, intent(in) :: nfield               ! field index
    character(len=10)   :: dimname              ! Dimension names
!
! !REVISION HISTORY:
! Author: Erik Kluzek
!
! !LOCAL VARIABLES:
    character(len=*),  parameter :: subname = 'harvestData::getFieldsDim'
!EOP
!-----------------------------------------------------------------------
    if ( this%dims2nd(nfield) == this%CFTdimsize )then
       dimname = "cft"
    else if ( this%dims2nd(nfield) == this%PFTdimsize )then
       dimname = "natpft"
    else
       dimname = "none"
    end if
  end function  getFieldsDim

!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: isField1D
!
! !INTERFACE:
  logical function isField1D( this, nfield )
!
! !DESCRIPTION:
!   Returns true if this field index is a 1D field
!
! !USES:
    implicit none
!
! !ARGUMENTS:
    class(harvestDataType), intent(IN) :: this  ! harvestData object
    integer, intent(in) :: nfield               ! field index
!
! !REVISION HISTORY:
! Author: Erik Kluzek
!
! !LOCAL VARIABLES:
    character(len=*),  parameter :: subname = 'harvestData::isField1D'
!EOP
!-----------------------------------------------------------------------
    isField1D = .false.
    if ( mkharvest_fieldInBounds( nfield ) )then
       if ( this%dims2nd(nfield) == 0 ) isField1D = .true.
    else
       call abort()
    end if
  end function isField1D

!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: isField2D
!
! !INTERFACE:
  logical function isField2D( this, nfield )
!
! !DESCRIPTION:
!   Returns true if this field index is a 2D field
!
! !USES:
    implicit none
!
! !ARGUMENTS:
    class(harvestDataType), intent(IN) :: this  ! harvestData object
    integer, intent(in) :: nfield              ! field index
!
! !REVISION HISTORY:
! Author: Erik Kluzek
!
! !LOCAL VARIABLES:
    character(len=*),  parameter :: subname = 'harvestData::isField2D'
!EOP
!-----------------------------------------------------------------------
    isField2D = .false.
    if ( mkharvest_fieldInBounds( nfield ) )then
       if ( this%dims2nd(nfield) /= 0 ) isField2D = .true.
    else
       call abort()
    end if
  end function isField2D

!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: num1DFields
!
! !INTERFACE:
  integer function num1DFields( this )
!
! !DESCRIPTION:
!   Returns the number of 1D fields
!
! !USES:
    implicit none
!
! !ARGUMENTS:
    class(harvestDataType), intent(IN) :: this  ! harvestData object
!
! !REVISION HISTORY:
! Author: Erik Kluzek
!
! !LOCAL VARIABLES:
    character(len=*),  parameter :: subname = 'harvestData::num1DFields'
!EOP
!-----------------------------------------------------------------------
    num1DFields = count( this%dims2nd == 0)
  end function num1DFields

!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: num2DFields
!
! !INTERFACE:
  integer function num2DFields( this )
!
! !DESCRIPTION:
!   Returns the number of 2D fields
!
! !USES:
    implicit none
!
! !ARGUMENTS:
    class(harvestDataType), intent(IN) :: this  ! harvestData object
!
! !REVISION HISTORY:
! Author: Erik Kluzek
!
! !LOCAL VARIABLES:
    character(len=*),  parameter :: subname = 'harvestData::num2DFields'
!EOP
!-----------------------------------------------------------------------
    num2DFields = count( this%dims2nd /= 0)
  end function num2DFields

!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: mkharvest_init
!
! !INTERFACE:
  subroutine mkharvest_init( ns_o, init_val, harvdata, fharvest, constant )
!
! !DESCRIPTION:
!              Initialization of mkharvest module.
!
! !USES:
    use mkncdio 
    implicit none
!
! !ARGUMENTS:
    integer              , intent(in) :: ns_o          ! clm output grid resolution
    real(r8)             , intent(in) :: init_val      ! initial value to set to
    type(harvestDataType), intent(INOUT) :: harvdata   ! Harvest data
    character(len=*)     , intent(in) :: fharvest      ! input harvest dataset file name
    logical, intent(in), optional :: constant          ! Flag if variables are CONST_ version for surface dataset 
                                                       ! rather than landuse.timeseries
!
! !REVISION HISTORY:
! Author: Erik Kluzek
!
!
! !LOCAL VARIABLES:
    character(len=*),  parameter :: subname = 'mkharvest_init'
    character(len=CL) :: lunits                 ! local units read in
    integer  :: ncid,varid                      ! input netCDF id's
    integer  :: ifld                            ! indices
    integer  :: ret                             ! return code
    logical  :: lconstant                       ! local version of constant flag
    logical  :: varexists                       ! If variable exists on file
    integer  :: dim_lengths(3)                  ! Dimension lengths on file
    integer  :: dims2nd(numharv)                ! Dimension lengths of 3rd dimension for each variable on file
    integer  :: ndims                           ! Number of dimensions on file
    integer  :: ns_i                            ! clm input grid resolution (nlat*nlon)
!EOP
!-----------------------------------------------------------------------
    lconstant = .false.
    if ( present(constant) ) lconstant = constant

    initialized = .true.
    call check_ret(nf_open(fharvest, 0, ncid), subname)
    dims2nd(:) = 0
    ns_i       = 0
    do ifld = 1, numharv
       call check_ret(nf_inq_varid (   ncid, mkharvest_fieldname(ifld, constant=lconstant), varid), subname, varexists=varexists)
       if ( .not. varexists )then
          write(*,*) "SKIP: "//mkharvest_fieldname(ifld, constant=lconstant)
          harvest_longnames(ifld) = trim(mkharvest_fieldname(ifld, constant=lconstant)) // " (zeroed out)"
          harvest_units(ifld)     = "not_read_in"
       else
          call check_ret(nf_get_att_text( ncid, varid, 'long_name', harvest_longnames(ifld)), subname )
          ret = nf_get_att_text( ncid, varid, 'units',     harvest_units(ifld))
          if ( ret == nf_enotatt )then
             harvest_units(ifld) = "unitless"
          else if ( ret == nf_noerr )then
          else
             write(*,*) 'ERROR:: bad return code from NetCDF get attribute= '// nf_strerror(ret)
             call abort()
          end if
          call get_dim_lengths(ncid, mkharvest_fieldname(ifld, constant=lconstant), ndims, dim_lengths)
          if ( ns_i == 0 )then
             ns_i = dim_lengths(1)*dim_lengths(2)
          else if ( ns_i /= dim_lengths(1)*dim_lengths(2) )then
             write(*,*) 'ERROR:: bad dimension sizes for variable = ', mkharvest_fieldname(ifld, constant=lconstant)
             call abort()
          end if
          if (      ndims == 2 )then
             dims2nd(ifld) = 0
          else if ( ndims == 3 )then
             dims2nd(ifld) = dim_lengths(3)
          else
             write(*,*) 'ERROR:: bad dimensionality for variable = ', mkharvest_fieldname(ifld, constant=lconstant)
             call abort()
          end if
             
       end if
    end do
    call harvdata%init( dims2nd, ns_i, ns_o, init_val )

    call check_ret(nf_close(ncid), subname)

    allocate( oride_harv(numharv) )
    oride_harv(:) = real_undef

  end subroutine mkharvest_init

!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: mkharvest_fieldInBounds
!
! !INTERFACE:
   logical function  mkharvest_fieldInBounds( nfield )
!
! !DESCRIPTION:
!   Return true if field index is in bounds and initialization done
!
! !USES:
    implicit none
!
! !ARGUMENTS:
    integer, intent(in) :: nfield    ! field index
!
! !REVISION HISTORY:
! Author: Erik Kluzek
!
! !LOCAL VARIABLES:
    character(len=*), parameter :: subname = 'mkharvest_fieldInBounds'
!EOP
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

!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: mkharvest_fieldname
!
! !INTERFACE:
  character(len=harlen) function mkharvest_fieldname( nfield, constant )
!
! !DESCRIPTION:
!              Return harvest fieldname of input field number.
!
! !USES:
    implicit none
!
! !ARGUMENTS:
    integer, intent(in) :: nfield
    logical, intent(in), optional :: constant
!
! !REVISION HISTORY:
! Author: Erik Kluzek
!
! !LOCAL VARIABLES:
    character(len=*), parameter :: subname = 'mkharvest_fieldname'
    logical :: lconstant     ! local version of constant flag
!EOP
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
         call abort()
      end if

  end function mkharvest_fieldname

!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: mkharvest_units
!
! !INTERFACE:
  character(len=CL) function mkharvest_units( nfield )
!
! !DESCRIPTION:
!              Return units description of harvest fields
!
! !USES:
    implicit none
!
! !ARGUMENTS:
    integer, intent(in) :: nfield
!
! !REVISION HISTORY:
! Author: Erik Kluzek
!
!
! !LOCAL VARIABLES:
    character(len=*), parameter :: subname = 'mkharvest_units'
!EOP
!-----------------------------------------------------------------------

      if ( mkharvest_fieldInBounds( nfield ) )then
         mkharvest_units    = harvest_units(nfield)
      else
         call abort()
      end if

  end function mkharvest_units

!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: mkharvest_longname
!
! !INTERFACE:
  character(len=CL) function mkharvest_longname( nfield )
!
! !DESCRIPTION:
!              Return longname description of given input field number.
!
! !USES:
    implicit none
!
! !ARGUMENTS:
    integer, intent(in) :: nfield
!
! !REVISION HISTORY:
! Author: Erik Kluzek
!
!
! !LOCAL VARIABLES:
    character(len=*), parameter :: subname = 'mkharvest_longname'
!EOP
!-----------------------------------------------------------------------

      if ( mkharvest_fieldInBounds( nfield ) )then
         mkharvest_longname = harvest_longnames(nfield)
      else
         call abort()
      end if

  end function mkharvest_longname

!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: mkharvest_numtypes
!
! !INTERFACE:
  integer function mkharvest_numtypes( )
!
! !DESCRIPTION:
!              Return number of different harvest field types.
!
! !USES:
    implicit none
!
! !ARGUMENTS:
    character(len=*), parameter :: subname = 'mkharvest_numtypes'
!
! !REVISION HISTORY:
! Author: Erik Kluzek
!
!
! !LOCAL VARIABLES:
!EOP
!-----------------------------------------------------------------------
      mkharvest_numtypes = numharv

  end function mkharvest_numtypes

!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: clean
!
! !INTERFACE:
  subroutine clean( this )
!
! !DESCRIPTION:
!   Clean and deallocate the harvestData object
!
! !USES:
    implicit none
!
! !ARGUMENTS:
    class(harvestDataType), intent(INOUT) :: this       ! harvestData object
!
! !REVISION HISTORY:
! Author: Erik Kluzek
!
! !LOCAL VARIABLES:
    character(len=*),  parameter :: subname = 'harvestData::clean'
!EOP
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

!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: mkharvest
!
! !INTERFACE:
subroutine mkharvest(ldomain, mapfname, datfname, ndiag, harvdata)
!
! !DESCRIPTION:
! Make harvest data for the dynamic PFT dataset.
! This dataset consists of the normalized harvest or grazing fraction (0-1) of
! the model.
!
! !USES:
  use mkdomainMod, only : domain_type, domain_clean, domain_read
  use mkgridmapMod
  use mkvarpar	
  use mkvarctl    
  use mkncdio
!
! !ARGUMENTS:
  implicit none
  type(domain_type)    , intent(in) :: ldomain     ! 
  character(len=*)     , intent(in) :: mapfname    ! input mapping file name
  character(len=*)     , intent(in) :: datfname    ! input data file name
  integer              , intent(in) :: ndiag       ! unit number for diag out
  type(harvestDataType), intent(INOUT) :: harvdata ! Harvest data
!
! !CALLED FROM:
! subroutine mksrfdat in module mksrfdatMod
!
! !REVISION HISTORY:
! Author: Erik Kluzek
!
!
! !LOCAL VARIABLES:
!EOP
  type(gridmap_type)    :: tgridmap
  type(domain_type)     :: tdomain            ! local domain
  real(r8) :: gharv_o(numharv)                ! output grid: global area harvesting
  real(r8) :: garea_o                         ! output grid: global area
  real(r8) :: gharv_i(numharv)                ! input grid: global area harvesting
  real(r8) :: garea_i                         ! input grid: global area
  integer  :: ifld                            ! indices
  integer  :: k,n,m,ni,no,ns_i,ns_o           ! indices
  integer  :: ncid,varid                      ! input netCDF id's
  logical  :: varexists                       ! If variable exists or not
  integer  :: ier                             ! error status
  integer, allocatable :: ind1D(:)            ! Index of 1D harvest fields
  integer, allocatable :: ind2D(:)            ! Index of 2D harvest fields
  real(r8), pointer :: data1D_i(:)            ! 1D input data
  real(r8), pointer :: data2D_i(:,:)          ! 2D output data
  real(r8), pointer :: data1D_o(:)            ! 1D output data
  real(r8), pointer :: data2D_o(:,:)          ! 2D output data

  character(len=*), parameter :: unit = '10**6 km**2' ! Output units
  real(r8), parameter :: fac = 1.e-06_r8              ! Output factor
  real(r8), parameter :: rat = fac/100._r8            ! Output factor divided by 100%
  character(len=*), parameter :: subname = 'mkharvest'
!-----------------------------------------------------------------------

  write (6,*) 'Attempting to make harvest fields .....'
  call shr_sys_flush(6)

  ! -----------------------------------------------------------------
  ! Normally read in the harvesting file, and then regrid to output grid
  ! -----------------------------------------------------------------
  call harvdata%getFieldsIdx( ind1D, ind2D )
  
  if ( all(oride_harv == real_undef ) )then

     ! -----------------------------------------------------------------
     ! Read input harvesting file
     ! -----------------------------------------------------------------

     ! Obtain input grid info, read HARVEST_VH1, HARVEST_VH2, ... GRAZING etc.

     call domain_read(tdomain,datfname)
     ns_i = tdomain%ns
     ns_o = ldomain%ns

     write (6,*) 'Open harvest file: ', trim(datfname)
     call check_ret(nf_open(datfname, 0, ncid), subname)
     do k = 1, harvdata%num1Dfields()
        ifld = ind1D(k)
        call check_ret( nf_inq_varid(ncid, mkharvest_fieldname(ifld), varid), subname, varexists=varexists )
        data1D_i => harvdata%get1DFieldPtr( ifld )
        if ( .not. varexists )then
           write(*,*) "SKIP: "//mkharvest_fieldname(ifld)
           data1D_i(:) = 0.0_r8
        else
           call check_ret(nf_get_var_double (ncid, varid, data1D_i), subname)
        end if 
     end do
     do k = 1, harvdata%num2Dfields()
        ifld = ind2D(k)
        call check_ret( nf_inq_varid(ncid, mkharvest_fieldname(ifld), varid), subname, varexists=varexists )
        data2D_i => harvdata%get2DFieldPtr( ifld )
        if ( .not. varexists )then
           write(*,*) "SKIP: "//mkharvest_fieldname(ifld)
           data2D_i(:,:) = 0.0_r8
        else
           call check_ret(nf_get_var_double (ncid, varid, data2D_i), subname)
        end if 
     end do
     call check_ret(nf_close(ncid), subname)

     ! Area-average normalized harvest on input grid [data*_i] to output grid [data*_o]

     call gridmap_mapread(tgridmap, mapfname )

     ! Error checks for domain and map consistencies
     
     call domain_checksame( tdomain, ldomain, tgridmap )

     ! Determine data* on output grid

     do k = 1, harvdata%num1Dfields()
        ifld = ind1D(k)
        data1D_i => harvdata%get1DFieldPtr( ifld )
        data1D_o => harvdata%get1DFieldPtr( ifld, output=.true. )
        call gridmap_areaave(tgridmap, data1D_i, data1D_o, nodata=0._r8)
     end do
     do k = 1, harvdata%num2Dfields()
        ifld = ind2D(k)
        data2D_i => harvdata%get2DFieldPtr( ifld )
        data2D_o => harvdata%get2DFieldPtr( ifld, output=.true. )
        do m = lbound(data2D_i(:,:),dim=2), ubound(data2D_i(:,:),dim=2)
           call gridmap_areaave(tgridmap, data2D_i(:,m), data2D_o(:,m), nodata=0._r8)
        end do
     end do

     ! -----------------------------------------------------------------
     ! Error check
     ! Compare global areas on input and output grids
     ! -----------------------------------------------------------------

     gharv_i(:) = 0.
     garea_i = 0.
     do ni = 1, ns_i
        garea_i = garea_i + tgridmap%area_src(ni)*re**2
        do k = 1, harvdata%num1Dfields()
           m = ind1D(k)
           data1D_i => harvdata%get1DFieldPtr( m )
           gharv_i(m) = gharv_i(m) + data1D_i(ni)*tgridmap%area_src(ni)* &
                                                tgridmap%frac_src(ni)*re**2
        end do
     end do

     gharv_o(:) = 0.
     garea_o = 0.
     do no = 1,ns_o
        garea_o = garea_o + tgridmap%area_dst(no)*re**2
        do k = 1, harvdata%num1Dfields()
           m = ind1D(k)
           data1D_o => harvdata%get1DFieldPtr( m, output=.true. )
           gharv_o(m) = gharv_o(m) + data1D_o(no)*tgridmap%area_dst(no)* &
                                                tgridmap%frac_dst(no)*re**2
        end do
     end do

     ! Write out to diagnostic output file
     !

     write (ndiag,*)
     write (ndiag,'(1x,70a1)') ('=',k=1,70)
     write (ndiag,*) 'Harvesting Output'
     write (ndiag,'(1x,70a1)') ('=',k=1,70)

     write (ndiag,*)
     write (ndiag,'(1x,70a1)') ('.',k=1,70)
     write (ndiag,1001) unit, unit
1001 format (1x,'harvest type   ',20x,' input grid area',' output grid area',/ &
             1x,33x,'     ',A,'      ',A)
     write (ndiag,'(1x,70a1)') ('.',k=1,70)
     write (ndiag,*)
     do k = 1, harvdata%num1Dfields()
        m = ind1D(k)
        write (ndiag,1002) mkharvest_fieldname(m), gharv_i(m)*rat,gharv_o(m)*rat
     end do
1002 format (1x,a35,f16.3,f17.3)

     ! Deallocate dynamic memory

     call domain_clean(tdomain) 
     call gridmap_clean(tgridmap)

  else

     ! -----------------------------------------------------------------
     ! Otherwise override the harvesting with the input harvest values
     ! -----------------------------------------------------------------
  
     if ( any(oride_harv == real_undef ) )then
         write(6,*) subname, ' error some override harvesting fields set ', &
                    'and others are not = ', oride_harv
         call abort()
     end if
     do k = 1, harvdata%num1Dfields()
        m = ind1D(k)
        if ( oride_harv(m) < 0.0_r8 .or. oride_harv(m) > 100.0_r8 )then
            write(6,*) subname, ' error override harvesting field out of range', &
                       oride_harv(m), ' field = ', mkharvest_fieldname(m)
            call abort()
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

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: mkharvest_parse_oride
!
! !INTERFACE:
subroutine mkharvest_parse_oride( string )
!
! !DESCRIPTION:
! Parse the string with harvest and grazing information on it, to override
! the file with this information rather than reading from a file.
!
! !USES:
   use shr_string_mod, only: shr_string_betweenTags
! !ARGUMENTS:
   character(len=256), intent(IN) :: string  ! String to parse with harvest and grazing data
!
! !CALLED FROM:
! subroutine mksrfdat in module mksrfdatMod
!
! !REVISION HISTORY:
! Author: Erik Kluzek
!
!
! !LOCAL VARIABLES:
!EOP
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
     call abort()
  end if
  read(substring,*) oride_harv(1:numharv-1)
  call shr_string_betweenTags( string, graz_start, graz_end, substring, rc )
  if ( rc /= 0 )then
     write(6,*) subname//'Trouble finding grazing start end tags'
     call abort()
  end if
  read(substring,*) oride_harv(numharv)
  if ( harvest_fieldnames(numharv) /= 'GRAZING' )then
     write(6,*) subname, ' grazing is NOT last field as was expected'
     call abort()
  end if

!-----------------------------------------------------------------------

end subroutine mkharvest_parse_oride

!-----------------------------------------------------------------------

end module mkharvestMod
