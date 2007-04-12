#include <misc.h>
#include <preproc.h>

module ndepFileMod

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
  subroutine ndeprd(ndep)
!
! !DESCRIPTION: 
! Read the nitrogen deposition dataset.
!
! !USES:
    use shr_kind_mod, only: r8 => shr_kind_r8
    use clm_varctl  , only : fndepdat, single_column
    use fileutils   , only : getfil
    use decompMod   , only : get_proc_bounds,gsMap_lnd_gdc2glo,perm_lnd_gdc2glo
!
! !ARGUMENTS:
    implicit none
    include 'netcdf.inc'
    real(r8), pointer :: ndep(:)         ! annual nitrogen deposition rate (gN/m2/yr)
!
! !CALLED FROM:
! subroutine initialize in module initializeMod
!
! !REVISION HISTORY:
! Created by Peter Thornton, 1 June 2004
!
!EOP
!
! !LOCAL VARIABLES:
    character(len=256) :: locfn                          ! local file name
    integer  :: ncid,dimid,varid                         ! netCDF id's
    integer  :: begg,endg                                ! start/stop gridcells
    integer  :: ier                                      ! error status 
    character(len=32) :: subname = 'ndeprd'              ! subroutine name
!-----------------------------------------------------------------------

    call get_proc_bounds(begg,endg)

    ! Initialize data to zero - no nitrogen deposition

    ndep(:)   = 0._r8
       
    ! read data if file was specified in namelist
       
    if (fndepdat /= ' ') then

       ! Obtain netcdf file and read surface data

       if (masterproc) then

          write (6,*) 'Attempting to read nitrogen deposition data .....'

          call getfil (fndepdat, locfn, 0)
          call check_ret(nf_open(locfn, 0, ncid), subname)

          if ( .not. single_column )then
             call check_dim(ncid, 'lon' , lsmlon)
             call check_dim(ncid, 'lat' , lsmlat)
          else
             lsmlon = 1
             lsmlat = 1
          end if
       endif 

       if ( .not. single_column )then
          call ncd_iolocal(ncid,'NDEP_year','read',ndep,begg,endg,gsMap_lnd_gdc2glo,perm_lnd_gdc2glo)
       else
          call endrun('ndeprd not implemented for SCAM' )
       end if

    endif

    if ( masterproc )then
       write (6,*) 'Successfully read nitrogen deposition data'
       write (6,*)
    end if

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
    use decompMod   , only : get_proc_global, ldecomp
    use clm_time_manager, only : get_curr_date
    use clm_varctl  , only : fndepdyn
    use fileutils   , only : getfil
!
! !ARGUMENTS:
    implicit none
    include 'netcdf.inc'
!
!EOP
!
! !LOCAL VARIABLES:
    integer  :: i,j,m,n,g                       ! indices
    integer  :: ntimes                          ! number of input time samples
    real(r8) :: sumpct                          ! sum for error check
    integer  :: varid                           ! netcdf ids
    integer  :: year                            ! year (0, ...) for nstep+1
    integer  :: mon                             ! month (1, ..., 12) for nstep+1
    integer  :: day                             ! day of month (1, ..., 31) for nstep+1
    integer  :: sec                             ! seconds into current date for nstep+1
    integer  :: ier                             ! error status
    logical  :: found                           ! true => input dataset bounding dates found
    integer  :: numg                            ! total number of gridcells across all procs
    integer  :: numl                            ! total number of landunits across all procs
    integer  :: numc                            ! total number of columns across all procs
    integer  :: nump                            ! total number of pfts across all procs
    integer , allocatable :: landmask_ndepdyn(:,:) ! input landmask
    real(r8), allocatable :: ndep(:,:)          ! input ndep
    type(gridcell_type), pointer :: gptr        ! pointer to gridcell derived subtype
    character(len=256) :: locfn                 ! local file name
    character(len= 32) :: subname='ndepdyn_init' ! subroutine name
 !-----------------------------------------------------------------------

    allocate(landmask_ndepdyn(lsmlon,lsmlat))

    ! Set pointers into derived type

    gptr => clm3%g

    ! Get relevant sizes

    call get_proc_global(numg, numl, numc, nump)


    allocate(ndepdyn1(numg), ndepdyn2(numg), ndepdyn(numg), stat=ier)
    if (ier /= 0) then
       write(6,*)'ndepdyn_init allocation error for ndepdyn1, ndepdyn2, ndepdyn'
       call endrun()
    end if

    if (masterproc) then
       
       ! Obtain file

       write (6,*) 'Attempting to read dynamic ndep data .....'
       call getfil (fndepdyn, locfn, 0)
       call check_ret(nf_open(locfn, 0, ncid), subname)

       ! Obtain ndep years from dynamic ndep file

       call check_ret(nf_inq_dimid(ncid, 'time', varid), subname)
       call check_ret(nf_inq_dimlen(ncid, varid, ntimes), subname)

       allocate (yearsndep(ntimes), stat=ier)
       if (ier /= 0) then
          write(6,*)'ndepdyn_init allocation error for yearsndep'; call endrun()
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
             write(6,*)'ndepdyn_init error: model year not found in ndepdyn timeseries'
             write(6,*)'model year = ',year
             call endrun()
          end if
       end if
            
       ! Get ndep time samples bracketing the current time

       allocate (ndep(lsmlon,lsmlat), stat=ier)
       if (ier /= 0) then
          write(6,*)subname,' allocation error for ndep'; call endrun()
       end if

       call ndepdyn_getdata(nt1, ndep)

!dir$ concurrent
!cdir nodep
       do g = 1,numg
          i = ldecomp%gdc2i(g)
          j = ldecomp%gdc2j(g)
          ndepdyn1(g) = ndep(i,j)
       end do
       
       call ndepdyn_getdata(nt2, ndep)

!dir$ concurrent
!cdir nodep
       do g = 1,numg
          i = ldecomp%gdc2i(g)
          j = ldecomp%gdc2j(g)
          ndepdyn2(g) = ndep(i,j)
       end do

       deallocate(ndep)

    end if   ! end of if-masterproc block

    call mpi_bcast (nt1   , 1, MPI_INTEGER, 0, mpicom, ier)
    call mpi_bcast (nt2   , 1, MPI_INTEGER, 0, mpicom, ier)
    call mpi_bcast (ntimes, 1, MPI_INTEGER, 0, mpicom, ier)
    if (.not. masterproc) then
       allocate(yearsndep(ntimes))
    end if
    call mpi_bcast (yearsndep, size(yearsndep), MPI_INTEGER, 0, mpicom, ier)
    call mpi_bcast (ndepdyn1  , size(ndepdyn1)  , MPI_REAL8, 0, mpicom, ier)
    call mpi_bcast (ndepdyn2  , size(ndepdyn2)  , MPI_REAL8, 0, mpicom, ier)

    deallocate(landmask_ndepdyn)

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
    use decompMod   , only : get_proc_global, get_proc_bounds, ldecomp
    use clm_atmlnd  , only : clm_a2l
    use clm_varcon  , only : istsoil
!
! !ARGUMENTS:
    implicit none
!
!EOP
!
! !LOCAL VARIABLES:
    integer  :: i,j,m,p,l,g      ! indices
    integer  :: year             ! year (0, ...) for nstep+1
    integer  :: mon              ! month (1, ..., 12) for nstep+1
    integer  :: day              ! day of month (1, ..., 31) for nstep+1
    integer  :: sec              ! seconds into current date for nstep+1
    real(r8) :: cday             ! current calendar day (1.0 = 0Z on Jan 1)
    integer  :: numg             ! total number of gridcells across all processors
    integer  :: numl             ! total number of landunits across all processors
    integer  :: numc             ! total number of columns across all processors
    integer  :: nump             ! total number of pfts across all processors
    integer  :: begp, endp       ! per-proc beginning and ending pft indices
    integer  :: begc, endc       ! per-proc beginning and ending column indices
    integer  :: begl, endl       ! per-proc beginning and ending landunit indices
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

    call get_proc_global(numg, numl, numc, nump)
    call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)

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
          write(6,*)subname,' error - current year is past input data boundary'
       end if
       
       if (masterproc) then
!dir$ concurrent
!cdir nodep
          do g = 1,numg
             ndepdyn1(g) = ndepdyn2(g)
          end do

          allocate(ndep(lsmlon,lsmlat), stat=ier)
          if (ier /= 0) then
             write(6,*)subname,' allocation error for ndep'; call endrun()
          end if
          call ndepdyn_getdata(nt2, ndep)

!dir$ concurrent
!cdir nodep
          do g = 1,numg
             i = ldecomp%gdc2i(g)
             j = ldecomp%gdc2j(g)
             ndepdyn2(g) = ndep(i,j)
          end do
          deallocate(ndep)
       end if  ! end of if-masterproc if-block

       call mpi_bcast (ndepdyn1, size(ndepdyn1), MPI_REAL8, 0, mpicom, ier)
       call mpi_bcast (ndepdyn2, size(ndepdyn2), MPI_REAL8, 0, mpicom, ier)
    
    end if  ! end of need new data if-block 

    ! Interpolate ndep to current time

    cday = get_curr_calday() 

    wt1 = ((days_per_year + 1._r8) - cday)/days_per_year
    
    ! assign interpolated flux field to forc_ndep
    ! convert units from gN/yr -> gN/s
    
!dir$ concurrent
!cdir nodep
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
    use shr_kind_mod, only : r8 => shr_kind_r8
!
! !ARGUMENTS:
    implicit none
    include 'netcdf.inc'
    integer , intent(in)  :: ntime
    real(r8), intent(out) :: ndep(:,:)
!
!EOP
!
! !LOCAL VARIABLES:
    integer  :: i,j,m
    integer  :: err, ierr, jerr, sumerr
    integer  :: varid                             ! netcdf variable id
    integer  :: beg3d(3), end3d(3), len3d(3)      ! input sizes
    character(len=32) :: subname='ndepdyn_getdata' ! subroutine name
!-----------------------------------------------------------------------
    
    beg3d(1) = 1     ;  len3d(1) = lsmlon
    beg3d(2) = 1     ;  len3d(2) = lsmlat
    beg3d(3) = ntime ;  len3d(3) = 1
    
    call check_ret(nf_inq_varid(ncid, 'NDEP_year', varid), subname)
    call check_ret(nf_get_vara_double(ncid, varid, beg3d, len3d, ndep), subname)

  end subroutine ndepdyn_getdata

end module ndepFileMod
