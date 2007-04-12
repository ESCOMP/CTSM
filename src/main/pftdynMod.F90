#include <misc.h>
#include <preproc.h>

module pftdynMod

!---------------------------------------------------------------------------
!BOP
!
! !MODULE: pftdynMod
!
! !USES:
  use spmdMod
  use clmtype
  use decompMod   , only : gsmap_lnd_gdc2glo,perm_lnd_gdc2glo
  use decompMod   , only : get_proc_bounds
  use ncdio
  use clm_varsur  , only : pctspec
  use clm_varpar  , only : max_pft_per_col
  use shr_kind_mod, only : r8 => shr_kind_r8
  use abortutils  , only : endrun
!
! !DESCRIPTION:
! Determine pft weights at current time using dynamic landuse datasets.
! ASSUMES that only have one dynamic landuse dataset.
!
! !PUBLIC TYPES:
  implicit none
  private
  save
  public :: pftdyn_init
  public :: pftdyn_interp
  public :: pftdyn_wbal_init
  public :: pftdyn_wbal
!
! !REVISION HISTORY:
! Created by Gordon Bonan, Sam Levis and Mariana Vertenstein
!
!EOP
!
! ! PRIVATE TYPES
  real(r8), parameter :: days_per_year = 365._r8
  integer , pointer   :: yearspft(:)
  real(r8), pointer   :: wtpft1(:,:)   
  real(r8), pointer   :: wtpft2(:,:)   
  real(r8), pointer   :: wtcol_old(:)
  integer :: nt1
  integer :: nt2
  integer :: ncid
!---------------------------------------------------------------------------

contains
  
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: pftdyn_init
!
! !INTERFACE:
  subroutine pftdyn_init()
!
! !DESCRIPTION:
! Initialize dynamic landuse dataset (position it to the right time samples to 
! that bound the initial model date
!
! !USES:
    use clm_time_manager, only : get_curr_date
    use clm_varctl  , only : fpftdyn
    use clm_varpar  , only : lsmlon, lsmlat, numpft
    use fileutils   , only : getfil
    use spmdGathScatMod, only : gather_data_to_master
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
    integer  :: begg,endg                       ! beg/end indices for land gridcells
    integer  :: begl,endl                       ! beg/end indices for land landunits
    integer  :: begc,endc                       ! beg/end indices for land columns
    integer  :: begp,endp                       ! beg/end indices for land pfts
    real(r8), pointer :: pctgla(:)          ! percent of gcell is glacier
    real(r8), pointer :: pctlak(:)          ! percent of gcell is lake
    real(r8), pointer :: pctwet(:)          ! percent of gcell is wetland
    real(r8), pointer :: pcturb(:)          ! percent of gcell is urbanized
    type(gridcell_type), pointer :: gptr        ! pointer to gridcell derived subtype
    character(len=256) :: locfn                 ! local file name
    character(len= 32) :: subname='pftdyn_init' ! subroutine name
 !-----------------------------------------------------------------------

    call get_proc_bounds(begg,endg,begl,endl,begc,endc,begp,endp)

    allocate(pctgla(begg:endg),pctlak(begg:endg))
    allocate(pctwet(begg:endg),pcturb(begg:endg))

    ! Set pointers into derived type

    gptr => clm3%g

    ! pctspec must be saved between time samples
    ! position to first time sample - assume that first time sample must match starting date
    ! check consistency -  special landunits, grid, frac and mask
    ! only do this once

    ! read data PCT_PFT corresponding to correct year

    allocate(wtpft1(begg:endg,0:numpft), wtpft2(begg:endg,0:numpft), stat=ier)
    if (ier /= 0) then
       write(6,*)'pctpft_dyn_init allocation error for wtpft1, wtpft2'
       call endrun()
    end if

    allocate(wtcol_old(begp:endp),stat=ier)
    if (ier /= 0) then
       write(6,*)'pctpft_dyn_init allocation error for wtcol_old'
       call endrun()
    end if

    if (masterproc) then

       ! Obtain file

       write (6,*) 'Attempting to read pft dynamic landuse data .....'
       call getfil (fpftdyn, locfn, 0)
       call check_ret(nf_open(locfn, 0, ncid), subname)

       ! Obtain pft years from dynamic landuse file

       call check_ret(nf_inq_dimid(ncid, 'time', varid), subname)
       call check_ret(nf_inq_dimlen(ncid, varid, ntimes), subname)

       ! Consistency check

       call check_dim(ncid, 'lsmpft', numpft+1)

    endif

    call mpi_bcast(ntimes,1,MPI_INTEGER,0,mpicom,ier)

    allocate (yearspft(ntimes), stat=ier)
    if (ier /= 0) then
       write(6,*)'pctpft_dyn_init allocation error for yearspft'; call endrun()
    end if

    if (masterproc) then
       call check_ret(nf_inq_varid(ncid, 'YEAR', varid), subname)
       call check_ret(nf_get_var_int(ncid, varid, yearspft), subname)
    endif

    call mpi_bcast(yearspft,ntimes,MPI_INTEGER,0,mpicom,ier)

    call ncd_iolocal(ncid, 'PCT_WETLAND', 'read', pctwet, begg, endg, gsMap_lnd_gdc2glo, perm_lnd_gdc2glo)
    call ncd_iolocal(ncid, 'PCT_LAKE'   , 'read', pctlak, begg, endg, gsMap_lnd_gdc2glo, perm_lnd_gdc2glo)
    call ncd_iolocal(ncid, 'PCT_GLACIER', 'read', pctgla, begg, endg, gsMap_lnd_gdc2glo, perm_lnd_gdc2glo)
    call ncd_iolocal(ncid, 'PCT_URBAN'  , 'read', pcturb, begg, endg, gsMap_lnd_gdc2glo, perm_lnd_gdc2glo)

    ! Consistency check
    do g = begg,endg
       if (pctlak(g)+pctwet(g)+pcturb(g)+pctgla(g) /= pctspec(g)) then 
          write(6,*)'mismatch between input pctspec = ',&
                     pctlak(g)+pctwet(g)+pcturb(g)+pctgla(g),&
                    ' and that obtained from surface dataset ', pctspec(g),' at g= ',g
           call endrun()
       end if
    end do

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

    if (year < yearspft(1)) then
       nt1 = 1
       nt2 = 1
    else if (year >= yearspft(ntimes)) then
       nt1 = ntimes
       nt2 = ntimes
    else
       found = .false.
       do n = 1,ntimes-1 
          if (year == yearspft(n)) then
             nt1 = n
             nt2 = nt1 + 1
             found = .true.
          end if   
       end do
       if (.not. found) then
          write(6,*)'pftdyn_init error: model year not found in pftdyn timeseries'
          write(6,*)'model year = ',year
          call endrun()
       end if
    end if

    ! Get pctpft time samples bracketing the current time

    call pftdyn_getdata(nt1, wtpft1, begg,endg,0,numpft)
    call pftdyn_getdata(nt2, wtpft2, begg,endg,0,numpft)
    do m = 0,numpft
!dir$ concurrent
!cdir nodep
       do g = begg,endg
          wtpft1(g,m) = wtpft1(g,m)/100._r8
          wtpft2(g,m) = wtpft2(g,m)/100._r8
       end do
    end do
       
    deallocate(pctgla,pctlak,pctwet,pcturb)

  end subroutine pftdyn_init

!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: pftdyn_interp
!
! !INTERFACE:
  subroutine pftdyn_interp()
!
! !DESCRIPTION:
! Time interpolate dynamic landuse data to get pft weights for model time
!
! !USES:
    use clm_time_manager, only : get_curr_date, get_curr_calday
    use clm_varcon  , only : istsoil
    use clm_varpar  , only : numpft, lsmlon, lsmlat
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
    integer  :: ier              ! error status
    real(r8) :: wt1              ! time interpolation weights
    integer  :: begg,endg                       ! beg/end indices for land gridcells
    integer  :: begl,endl                       ! beg/end indices for land landunits
    integer  :: begc,endc                       ! beg/end indices for land columns
    integer  :: begp,endp                       ! beg/end indices for land pfts
    type(gridcell_type), pointer :: gptr         ! pointer to gridcell derived subtype
    type(landunit_type), pointer :: lptr         ! pointer to landunit derived subtype
    type(pft_type)     , pointer :: pptr         ! pointer to pft derived subtype
    character(len=32) :: subname='pftdyn_interp' ! subroutine name
!-----------------------------------------------------------------------

    ! Set pointers into derived type

    gptr => clm3%g
    lptr => clm3%g%l
    pptr => clm3%g%l%c%p

    call get_proc_bounds(begg,endg,begl,endl,begc,endc,begp,endp)

    ! Interpolat pctpft to current time step - output in pctpft_intp
    ! Map interpolated pctpft to subgrid weights
    ! assumes that maxpatch_pft = numpft + 1, that each landunit has only 1 column, 
    ! SCAM and DGVM have not been defined, and that create_croplandunit = .false.

    ! If necessary, obtain new time sample

    ! Get current date

    call get_curr_date(year, mon, day, sec)

    ! Obtain new time sample if necessary.
    ! The first condition is the regular crossing of a year boundary
    ! when within the dynpft timeseries range. The second condition is
    ! the case of the first entry into the dynpft timeseries range from
    ! an earlier period of constant weights.

    if (year > yearspft(nt1) .or. (nt1 == 1 .and. nt2 == 1 .and. year == yearspft(1))) then

       if (year >= yearspft(size(yearspft))) then
          nt1 = size(yearspft)
          nt2 = size(yearspft)
       else
          nt1 = nt2
          nt2 = nt1 + 1
       end if
       
       if (nt2 > size(yearspft)) then
          write(6,*)subname,' error - current year is past input data boundary'
       end if
       
       do m = 0,numpft
!dir$ concurrent
!cdir nodep
          do g = begg,endg
             wtpft1(g,m) = wtpft2(g,m)
          end do
       end do

       call pftdyn_getdata(nt2, wtpft2, begg,endg,0,numpft)

       do m = 0,numpft
!dir$ concurrent
!cdir nodep
          do g = begg,endg
             wtpft2(g,m) = wtpft2(g,m)/100._r8
          end do
       end do
    
    end if  ! end of need new data if-block 

    ! Interpolate pft weight to current time

    cday = get_curr_calday() 

    wt1 = ((days_per_year + 1._r8) - cday)/days_per_year

!dir$ concurrent
!cdir nodep
    do p = begp,endp
       g = pptr%gridcell(p)
       l = pptr%landunit(p)
       if (lptr%itype(l) == istsoil) then
          m = pptr%itype(p)
          wtcol_old(p)      = pptr%wtcol(p)
!         --- recoded for roundoff performance, tcraig 3/07 from k.lindsay
!         pptr%wtgcell(p)   = wtpft1(g,m)*wt1 + wtpft2(g,m)*wt2
          pptr%wtgcell(p)   = wtpft2(g,m) + wt1*(wtpft1(g,m)-wtpft2(g,m))
          pptr%wtlunit(p)   = pptr%wtgcell(p) / lptr%wtgcell(l)
          pptr%wtcol(p)     = pptr%wtgcell(p) / lptr%wtgcell(l)
       end if
    end do
    
  end subroutine pftdyn_interp

!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: pftdyn_get_data
!
! !INTERFACE:
  subroutine pftdyn_getdata(ntime, pctpft, lb1,ub1,lb2,ub2)
!
! !DESCRIPTION:
! Obtain dynamic landuse data (pctpft) and make sure that
! percentage of PFTs sum to 100% cover for vegetated landunit
!
! !USES:
    use clm_varpar  , only : numpft, lsmlon, lsmlat
!
! !ARGUMENTS:
    implicit none
    include 'netcdf.inc'
    integer , intent(in)  :: ntime
    integer , intent(in)  :: lb1,ub1,lb2,ub2
    real(r8), intent(out) :: pctpft(lb1:ub1,lb2:ub2)
!
!EOP
!
! !LOCAL VARIABLES:
    integer  :: i,j,m,n
    integer  :: begg,endg         
    integer  :: err, ierr
    real(r8) :: sumpct,sumerr                     ! temporary
    integer  :: start(4), count(4)                ! input sizes
    real(r8),pointer :: arrayl(:)                 ! temporary array
    character(len=32) :: subname='pftdyn_getdata' ! subroutine name
!-----------------------------------------------------------------------
    
    call get_proc_bounds(begg,endg)

    allocate(arrayl(begg:endg))
    do n = 0,numpft
       start(1) = 1
       count(1) = lsmlon
       start(2) = 1
       count(2) = lsmlat
       start(3) = n+1      ! dataset is 1:numpft+1, not 0:numpft
       count(3) = 1
       start(4) = ntime 
       count(4) = 1
       call ncd_iolocal(ncid, 'PCT_PFT', 'read', arrayl, begg, endg, gsMap_lnd_gdc2glo, perm_lnd_gdc2glo, start, count)
       pctpft(begg:endg,n) = arrayl(begg:endg)
    enddo
    deallocate(arrayl)

    err = 0
    do n = begg,endg
       if (pctspec(n) < 100._r8) then
          sumpct = 0._r8
          do m = 0, numpft
             sumpct = sumpct + pctpft(n,m) * 100._r8/(100._r8-pctspec(n))
          end do
          if (abs(sumpct - 100._r8) > 0.1_r8) then
             err = 1; ierr = n; sumerr = sumpct
          end if
          if (sumpct < -0.000001_r8) then
             err = 2; ierr = n; sumerr = sumpct
          end if
       end if
    end do
    if (err == 1) then
       write(6,*) subname,' error: sum(pct) over numpft+1 is not = 100.',sumerr,ierr,pctspec(ierr),pctpft(ierr,:)
       call endrun()
    else if (err == 2) then
       write(6,*)subname,' error: sum(pct) over numpft+1 is < 0.',sumerr,ierr,pctspec(ierr),pctpft(ierr,:)
       call endrun()
    end if
    
  end subroutine pftdyn_getdata

!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: pftdyn_wbal_init
!
! !INTERFACE:
  subroutine pftdyn_wbal_init()
!
! !DESCRIPTION:
! initialize the column-level mass-balance correction term.
! Called in every timestep.
!
! !USES:
!
! !ARGUMENTS:
    implicit none
!
!EOP
!
! !LOCAL VARIABLES:
    integer  :: begp, endp    ! proc beginning and ending pft indices
    integer  :: begc, endc    ! proc beginning and ending column indices
    integer  :: begl, endl    ! proc beginning and ending landunit indices
    integer  :: begg, endg    ! proc beginning and ending gridcell indices
    integer  :: c             ! indices
    type(column_type),   pointer :: cptr         ! pointer to column derived subtype
!-----------------------------------------------------------------------

    ! Set pointers into derived type

    cptr => clm3%g%l%c

    ! Get relevant sizes

    call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)

    ! set column-level canopy water mass balance correction flux
    ! term to 0 at the beginning of every timestep
    
!dir$ concurrent
!cdir nodep
    do c = begc,endc
       cptr%cwf%h2ocan_loss(c) = 0._r8
    end do
    
  end subroutine pftdyn_wbal_init

!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: pftdyn_wbal
!
! !INTERFACE:
  subroutine pftdyn_wbal()
!
! !DESCRIPTION:
! modify pft-level state and flux variables to maintain water balance with
! dynamic pft-weights.
!
! !USES:
    use clm_varcon  , only : istsoil
    use clm_time_manager, only : get_step_size
!
! !ARGUMENTS:
    implicit none
!
!EOP
!
! !LOCAL VARIABLES:
    integer  :: begp, endp    ! proc beginning and ending pft indices
    integer  :: begc, endc    ! proc beginning and ending column indices
    integer  :: begl, endl    ! proc beginning and ending landunit indices
    integer  :: begg, endg    ! proc beginning and ending gridcell indices
    integer  :: pi,p,c,l,g    ! indices
    integer  :: ier           ! error code
    real(r8) :: dwt           ! change in pft weight (relative to column)
    real(r8) :: dtime         ! land model time step (sec)
    real(r8) :: init_h2ocan   ! initial canopy water mass
    real(r8) :: new_h2ocan    ! canopy water mass after weight shift
    real(r8), allocatable :: loss_h2ocan(:) ! canopy water mass loss due to weight shift
    type(landunit_type), pointer :: lptr         ! pointer to landunit derived subtype
    type(column_type),   pointer :: cptr         ! pointer to column derived subtype
    type(pft_type)   ,   pointer :: pptr         ! pointer to pft derived subtype
    character(len=32) :: subname='pftdyn_wbal' ! subroutine name
!-----------------------------------------------------------------------

    ! Set pointers into derived type

    lptr => clm3%g%l
    cptr => clm3%g%l%c
    pptr => clm3%g%l%c%p

    ! Get relevant sizes

    call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)

    ! Allocate loss_h2ocan
    allocate(loss_h2ocan(begp:endp), stat=ier)
    if (ier /= 0) then
          write(6,*)subname,' allocation error for loss_h2ocan'; call endrun()
    end if

    ! Get time step

    dtime = get_step_size()

    ! set column-level canopy water mass balance correction flux
    ! term to 0 at the beginning of every weight-shifting timestep

!dir$ concurrent
!cdir nodep
    do c = begc,endc
       cptr%cwf%h2ocan_loss(c) = 0._r8
    end do

!dir$ concurrent
!cdir nodep
    do p = begp,endp
       l = pptr%landunit(p)
       loss_h2ocan(p) = 0._r8

       if (lptr%itype(l) == istsoil) then

          ! calculate the change in weight for the timestep
          dwt = pptr%wtcol(p)-wtcol_old(p)
  
          if (dwt > 0._r8) then
          
             ! if the pft gained weight, then the 
             ! initial canopy water state is redistributed over the
             ! new (larger) area, conserving mass.

             pptr%pws%h2ocan(p) = pptr%pws%h2ocan(p) * (wtcol_old(p)/pptr%wtcol(p))
          
          else if (dwt < 0._r8) then
          
             ! if the pft lost weight on the timestep, then the canopy water
             ! mass associated with the lost weight is directed to a 
             ! column-level flux term that gets added to the precip flux
             ! for every pft calculation in Hydrology1()
             
             init_h2ocan = pptr%pws%h2ocan(p) * wtcol_old(p)
             loss_h2ocan(p) = pptr%pws%h2ocan(p) * (-dwt)
             new_h2ocan = init_h2ocan - loss_h2ocan(p)
             if (abs(new_h2ocan) < 1e-8_r8) then
                new_h2ocan = 0._r8
                loss_h2ocan(p) = init_h2ocan
             end if
             if (pptr%wtcol(p) /= 0._r8) then  
                pptr%pws%h2ocan(p) = new_h2ocan/pptr%wtcol(p)
             else
                pptr%pws%h2ocan(p) = 0._r8
                loss_h2ocan(p) = init_h2ocan
             end if 
          

          end if
          
       end if
       
    end do

!dir$ nointerchange
    do pi = 1,max_pft_per_col
!dir$ concurrent
!cdir nodep
       do c = begc,endc
          if (pi <= cptr%npfts(c)) then
             p = cptr%pfti(c) + pi - 1
             cptr%cwf%h2ocan_loss(c) = cptr%cwf%h2ocan_loss(c) + loss_h2ocan(p)/dtime
          end if
       end do
    end do

    ! Deallocate loss_h2ocan
    deallocate(loss_h2ocan)
    
  end subroutine pftdyn_wbal

end module pftdynMod
