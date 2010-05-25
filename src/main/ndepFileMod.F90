module ndepFileMod
#ifdef CN

!----------------------------------------------------------------------- 
!BOP
!
! !MODULE: ndepFileMod
! 
! !DESCRIPTION: 
! Contains methods for reading in nitrogen deposition data file
! Also includes functions for dynamic ndep file handling and 
! interpolation.
!
! !USES
  use abortutils, only : endrun
  use ncdio
  use clmtype
  use spmdMod     
  use clm_varpar,   only: lsmlat, lsmlon
  use clm_varctl,   only: iulog
  use shr_kind_mod, only: r8 => shr_kind_r8
!
! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: ndeprd  ! Read nitrogen deposition dataset
  public :: ndepdyn_init  ! position datasets for dynamic ndep
  public :: ndepdyn_interp ! interpolates between two years of ndep file data
!
! !REVISION HISTORY:
! Created by Peter Thornton, 1 June 2004
! 2/5/05, PET: Added ndepdyn_init and ndepdyn_interp
!
!EOP
!
! ! PRIVATE TYPES
  real(r8), parameter :: days_per_year = 365._r8
  integer , pointer   :: yearsndep(:)
  real(r8), pointer   :: ndepdyn1(:)   
  real(r8), pointer   :: ndepdyn2(:)   
  real(r8), pointer   :: ndepdyn(:)
  integer :: nt1
  integer :: nt2
  integer :: ncid
!----------------------------------------------------------------------- 

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ndeprd
!
! !INTERFACE:
  subroutine ndeprd(fndepdat)
!
! !DESCRIPTION: 
! Read the nitrogen deposition dataset.
!
! !USES:
    use shr_kind_mod, only: r8 => shr_kind_r8
    use clm_varctl  , only : single_column
    use fileutils   , only : getfil
    use decompMod   , only : get_proc_bounds
    use clm_atmlnd  , only : clm_a2l
!
! !ARGUMENTS:
    implicit none
    include 'netcdf.inc'
    character(len=*), intent(in), optional :: fndepdat	
!
! !CALLED FROM:
! subroutine initialize in module initializeMod
!
! !REVISION HISTORY:
! Created by Peter Thornton, 1 June 2004
!
!
! !LOCAL VARIABLES:
!EOP
    character(len=256) :: locfn              ! local file name
    character(len=256) :: units = ' '        ! units of ndep_year variable
    integer  :: ncid,dimid,varid             ! netCDF id's
    integer  :: g,begg,endg                  ! start/stop gridcells
    integer  :: ier, ret                     ! error status 
    real(r8), pointer :: ndep(:)             ! annual nitrogen deposition rate (gN/m2/yr)
    character(len=32) :: subname = 'ndeprd'  ! subroutine name
!-----------------------------------------------------------------------

   ! Initialize data to zero - no nitrogen deposition

    call get_proc_bounds(begg,endg)
    allocate(ndep(begg:endg))
    ndep(:) = 0._r8

    ! Obtain netcdf file and read surface data if appropriate
    if (present(fndepdat)) then	
       if (masterproc) then

          write(iulog,*) 'Attempting to read nitrogen deposition data .....'
          call getfil (fndepdat, locfn, 0)
          call check_ret(nf_open(locfn, 0, ncid), subname)

          if ( .not. single_column )then
             call check_dim(ncid, 'lon' , lsmlon)
             call check_dim(ncid, 'lat' , lsmlat)
          else
             lsmlon = 1
             lsmlat = 1
          end if
          
          ! Check that units are correct on NDEP_year variable
          call check_ret(nf_inq_varid(ncid, 'NDEP_year', varid), subname)
          call check_ret(nf_get_att_text(ncid, varid, 'units', units), subname)
          if ( trim(units) .ne. "g(N)/m2/yr" .and. trim(units) .ne. "gN/m2/yr" )then
             call endrun( 'ERROR: units are NOT what is expected on ndepdat file = '//trim(units) )
          end if
       endif

       call ncd_iolocal(ncid,'NDEP_year','read', ndep, grlnd, status=ret)
       if (ret /= 0) call endrun( trim(subname)//' ERROR: NDEP_year NOT on ndepdat file' )
       write(iulog,*) 'Successfully read nitrogen deposition data'
       write(iulog,*)
    end if

    do g = begg,endg
       clm_a2l%forc_ndep(g) = ndep(g)/(86400._r8 * 365._r8)
    end do

  end subroutine ndeprd

!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: ndepdyn_init
!
! !INTERFACE:
  subroutine ndepdyn_init()
!
! !DESCRIPTION:
! Initialize dynamic nitrogen deposition dataset
!
! !USES:
    use shr_kind_mod, only : r8 => shr_kind_r8
    use decompMod   , only : get_proc_bounds
    use clm_time_manager, only : get_curr_date
    use clm_varctl  , only : fndepdyn
    use fileutils   , only : getfil
!
! !ARGUMENTS:
    implicit none
    include 'netcdf.inc'
!
!
! !LOCAL VARIABLES:
!EOP
    integer  :: i,j,m,n,g                        ! indices
    integer  :: ntimes                           ! number of input time samples
    real(r8) :: sumpct                           ! sum for error check
    integer  :: varid                            ! netcdf ids
    integer  :: year                             ! year (0, ...) for nstep+1
    integer  :: mon                              ! month (1, ..., 12) for nstep+1
    integer  :: day                              ! day of month (1, ..., 31) for nstep+1
    integer  :: sec                              ! seconds into current date for nstep+1
    integer  :: ier                              ! error status
    logical  :: found                            ! true => input dataset bounding dates found
    integer  :: begg,endg                        ! local beg/end indices
    type(gridcell_type), pointer :: gptr         ! pointer to gridcell derived subtype
    character(len=256) :: locfn                  ! local file name
    character(len=256) :: units = ' '            ! units of ndep_year variable
    character(len= 32) :: subname='ndepdyn_init' ! subroutine name
 !-----------------------------------------------------------------------

    ! Set pointers into derived type

    gptr => clm3%g

    ! Get relevant sizes

    call get_proc_bounds(begg,endg)

    allocate(ndepdyn1(begg:endg), ndepdyn2(begg:endg), ndepdyn(begg:endg), stat=ier)
    if (ier /= 0) then
       write(iulog,*)'ndepdyn_init allocation error for ndepdyn1, ndepdyn2, ndepdyn'
       call endrun()
    end if

    if (masterproc) then
       
       ! Obtain file

       write(iulog,*) 'Attempting to read dynamic ndep data .....'
       call getfil (fndepdyn, locfn, 0)
       call check_ret(nf_open(locfn, 0, ncid), subname)

       ! Obtain ndep years from dynamic ndep file

       call check_ret(nf_inq_dimid(ncid, 'time', varid), subname)
       call check_ret(nf_inq_dimlen(ncid, varid, ntimes), subname)

       allocate (yearsndep(ntimes), stat=ier)
       if (ier /= 0) then
          write(iulog,*)'ndepdyn_init allocation error for yearsndep'; call endrun()
       end if

       call check_ret(nf_inq_varid(ncid, 'YEAR', varid), subname)
       call check_ret(nf_get_var_int(ncid, varid, yearsndep), subname)

       ! Determine if current date spans the years
       ! If current year is less than first dynamic PFT timeseries year,
       ! then use the first year from dynamic pft file for both nt1 and nt2,
       ! forcing constant weights until the model year enters the dynamic
       ! pft dataset timeseries range.
       ! If current year is equal to or greater than the last dynamic pft
       ! timeseries year, then use the last year for both nt1 and nt2, 
       ! forcing constant weights for the remainder of the simulation.
       ! This mechanism permits the introduction of a dynamic pft period in the middle
       ! of a simulation, with constant weights before and after the dynamic period.

       call get_curr_date(year, mon, day, sec)

       if (year < yearsndep(1)) then
          nt1 = 1
          nt2 = 1
       else if (year >= yearsndep(ntimes)) then
          nt1 = ntimes
          nt2 = ntimes
       else
          found = .false.
          do n = 1,ntimes-1 
             if (year == yearsndep(n)) then
                nt1 = n
                nt2 = nt1 + 1
                found = .true.
             end if   
          end do
          if (.not. found) then
             write(iulog,*)'ndepdyn_init error: model year not found in ndepdyn timeseries'
             write(iulog,*)'model year = ',year
             call endrun()
          end if
       end if
       !
       ! Check that units are correct on NDEP_year variable
       !
       call check_ret(nf_inq_varid(ncid, 'NDEP_year', varid), subname)
       call check_ret(nf_get_att_text(ncid, varid, 'units', units), subname)
       if ( trim(units) .ne. "g(N)/m2/yr" )then
          call endrun( 'ERROR: units are NOT what is expected on fndepdyn file' )
       end if
            
    end if   ! end of if-masterproc block

    call mpi_bcast (nt1   , 1, MPI_INTEGER, 0, mpicom, ier)
    call mpi_bcast (nt2   , 1, MPI_INTEGER, 0, mpicom, ier)
    call mpi_bcast (ntimes, 1, MPI_INTEGER, 0, mpicom, ier)
    if (.not. masterproc) then
       allocate(yearsndep(ntimes))
    end if
    call mpi_bcast (yearsndep, size(yearsndep), MPI_INTEGER, 0, mpicom, ier)

    ! Get ndep time samples bracketing the current time
    call ndepdyn_getdata(nt1, ndepdyn1)
    call ndepdyn_getdata(nt2, ndepdyn2)

  end subroutine ndepdyn_init

!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: ndepdyn_interp
!
! !INTERFACE:
  subroutine ndepdyn_interp()
!
! !DESCRIPTION:
! Time interpolate dynamic ndep data to get ndep for model time
!
! !USES:
    use shr_kind_mod, only : r8 => shr_kind_r8
    use clm_time_manager, only : get_curr_date, get_curr_calday
    use decompMod   , only : get_proc_bounds
    use clm_atmlnd  , only : clm_a2l
!
! !ARGUMENTS:
    implicit none
!
!
! !LOCAL VARIABLES:
!EOP
    integer  :: i,j,m,p,l,g      ! indices
    integer  :: year             ! year (0, ...) for nstep+1
    integer  :: mon              ! month (1, ..., 12) for nstep+1
    integer  :: day              ! day of month (1, ..., 31) for nstep+1
    integer  :: sec              ! seconds into current date for nstep+1
    real(r8) :: cday             ! current calendar day (1.0 = 0Z on Jan 1)
    integer  :: begg, endg       ! per-proc gridcell ending gridcell indices
    integer  :: ier              ! error status
    real(r8) :: wt1              ! time interpolation weights
    real(r8), allocatable :: ndep(:,:)       ! input ndep
    type(gridcell_type), pointer :: gptr         ! pointer to gridcell derived subtype
    character(len=32) :: subname='ndepdyn_interp' ! subroutine name
!-----------------------------------------------------------------------

    ! Set pointers into derived type

    gptr => clm3%g

    ! Get relevant sizes

    call get_proc_bounds(begg, endg)

    ! Interpolat ndep to current time step - output in ndep

    ! If necessary, obtain new time sample

    ! Get current date

    call get_curr_date(year, mon, day, sec)

    ! Obtain new time sample if necessary.
    ! The first condition is the regular crossing of a year boundary
    ! when within the dynndep timeseries range. The second condition is
    ! the case of the first entry into the dynndep timeseries range from
    ! an earlier period of constant weights.
    
    if (year > yearsndep(nt1) .or. (nt1 == 1 .and. nt2 == 1 .and. year == yearsndep(1))) then

       if (year >= yearsndep(size(yearsndep))) then
          nt1 = size(yearsndep)
          nt2 = size(yearsndep)
       else
          nt1 = nt2
          nt2 = nt1 + 1
       end if
       
       if (nt2 > size(yearsndep)) then
          write(iulog,*)subname,' error - current year is past input data boundary'
       end if
       
       do g = begg,endg
          ndepdyn1(g) = ndepdyn2(g)
       end do

       call ndepdyn_getdata(nt2, ndepdyn2)

    end if  ! end of need new data if-block 

    ! Interpolate ndep to current time

    cday = get_curr_calday() 

    wt1 = ((days_per_year + 1._r8) - cday)/days_per_year
    
    ! assign interpolated flux field to forc_ndep
    ! convert units from gN/yr -> gN/s
    
    do g = begg,endg
!      --- recoded for roundoff error, tcraig 3/07 from k.lindsay
!      clm_a2l%forc_ndep(g) = (ndepdyn1(g)*wt1 + ndepdyn2(g)* wt2)/(86400._r8 * 365._r8)
       clm_a2l%forc_ndep(g) = (ndepdyn2(g) + wt1*(ndepdyn1(g)-ndepdyn2(g)))/(86400._r8 * 365._r8)
    end do

  end subroutine ndepdyn_interp

!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: ndepdyn_getdata
!
! !INTERFACE:
  subroutine ndepdyn_getdata(ntime, ndep)
!
! !DESCRIPTION:
! Obtain dynamic ndep 
!
! !USES:
    use clmtype     , only : grlnd
    use shr_kind_mod, only : r8 => shr_kind_r8
!
! !ARGUMENTS:
    implicit none
    include 'netcdf.inc'
    integer , intent(in)  :: ntime
    real(r8), pointer     :: ndep(:)
!
!
! !LOCAL VARIABLES:
!EOP
    integer  :: beg3d(3), len3d(3)      ! input sizes
    integer  :: ret                     ! error status
    character(len=32) :: subname='ndepdyn_getdata' ! subroutine name
!-----------------------------------------------------------------------
    
    beg3d(1) = 1     ;  len3d(1) = lsmlon
    beg3d(2) = 1     ;  len3d(2) = lsmlat
    beg3d(3) = ntime ;  len3d(3) = 1
    
    call ncd_iolocal(ncid,'NDEP_year','read',ndep,grlnd,beg3d,len3d,status=ret)
    if (ret /= 0) call endrun( trim(subname)//' ERROR: NDEP_year NOT on ndepdyn file' )

  end subroutine ndepdyn_getdata
#endif

end module ndepFileMod
