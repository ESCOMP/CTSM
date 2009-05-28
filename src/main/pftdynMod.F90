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
  use decompMod   , only : gsmap_lnd_gdc2glo
  use decompMod   , only : get_proc_bounds
  use ncdio
  use clm_varsur  , only : pctspec
  use clm_varpar  , only : max_pft_per_col
  use clm_varctl  , only : iulog
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
#ifdef CN
  public :: pftdyn_cnbal
#endif
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
    use clm_varpar  , only : lsmlon, lsmlat, numpft, maxpatch_pft
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
    integer  :: ier, ret                        ! error status
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

    ! Error check

    if ( maxpatch_pft /= numpft+1 )then
       call endrun( subname//' maxpatch_pft does NOT equal numpft+1 -- this is invalid for dynamic PFT case' )
    end if

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
       call endrun( subname//' allocation error for wtpft1, wtpft2' )
    end if

    allocate(wtcol_old(begp:endp),stat=ier)
    if (ier /= 0) then
       call endrun( subname//' allocation error for wtcol_old' )
    end if

    if (masterproc) then

       ! Obtain file

       write(iulog,*) 'Attempting to read pft dynamic landuse data .....'
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
       write(iulog,*) subname//' allocation error for yearspft'; call endrun()
    end if

    if (masterproc) then
       call check_ret(nf_inq_varid(ncid, 'YEAR', varid), subname)
       call check_ret(nf_get_var_int(ncid, varid, yearspft), subname)
    endif

    call mpi_bcast(yearspft,ntimes,MPI_INTEGER,0,mpicom,ier)

    call ncd_iolocal(ncid, 'PCT_WETLAND', 'read', pctwet, grlnd, status=ret)
    if (ret /= 0) call endrun( trim(subname)//' ERROR: PCT_WETLAND NOT on pftdyn file' )
    call ncd_iolocal(ncid, 'PCT_LAKE'   , 'read', pctlak, grlnd, status=ret)
    if (ret /= 0) call endrun( trim(subname)//' ERROR: PCT_LAKE NOT on pftdyn file' )
    call ncd_iolocal(ncid, 'PCT_GLACIER', 'read', pctgla, grlnd, status=ret)
    if (ret /= 0) call endrun( trim(subname)//' ERROR: PCT_GLACIER NOT on pftdyn file' )
    call ncd_iolocal(ncid, 'PCT_URBAN'  , 'read', pcturb, grlnd, status=ret)
    if (ret /= 0) call endrun( trim(subname)//' ERROR: PCT_URBAN NOT on pftdyn file' )

    ! Consistency check
    do g = begg,endg
       if (pctlak(g)+pctwet(g)+pcturb(g)+pctgla(g) /= pctspec(g)) then 
          write(iulog,*) subname//'mismatch between input pctspec = ',&
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
          write(iulog,*) subname//' error: model year not found in pftdyn timeseries'
          write(iulog,*)'model year = ',year
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
    integer  :: i,j,m,p,l,g,c    ! indices
    integer  :: year             ! year (0, ...) for nstep+1
    integer  :: mon              ! month (1, ..., 12) for nstep+1
    integer  :: day              ! day of month (1, ..., 31) for nstep+1
    integer  :: sec              ! seconds into current date for nstep+1
    real(r8) :: cday             ! current calendar day (1.0 = 0Z on Jan 1)
    integer  :: ier              ! error status
    integer  :: lbc,ubc
    real(r8) :: wt1              ! time interpolation weights
    real(r8), pointer :: wtpfttot1(:)           ! summation of pft weights for renormalization
    real(r8), pointer :: wtpfttot2(:)           ! summation of pft weights for renormalization
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

    allocate(wtpfttot1(begc:endc),wtpfttot2(begc:endc))
    wtpfttot1(:) = 0._r8
    wtpfttot2(:) = 0._r8

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
          write(iulog,*)subname,' error - current year is past input data boundary'
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
       c = pptr%column(p)
       g = pptr%gridcell(p)
       l = pptr%landunit(p)
       if (lptr%itype(l) == istsoil) then
          m = pptr%itype(p)
          wtcol_old(p)      = pptr%wtcol(p)
!         --- recoded for roundoff performance, tcraig 3/07 from k.lindsay
!         pptr%wtgcell(p)   = wtpft1(g,m)*wt1 + wtpft2(g,m)*wt2
          wtpfttot1(c) = wtpfttot1(c)+pptr%wtgcell(p)    
          pptr%wtgcell(p)   = wtpft2(g,m) + wt1*(wtpft1(g,m)-wtpft2(g,m))
          pptr%wtlunit(p)   = pptr%wtgcell(p) / lptr%wtgcell(l)
          pptr%wtcol(p)     = pptr%wtgcell(p) / lptr%wtgcell(l)
          wtpfttot2(c) = wtpfttot2(c)+pptr%wtgcell(p)
       end if

    end do

!   Renormalize pft weights so that sum of pft weights relative to grid cell 
!   remain constant even as land cover changes.  Doing this eliminates 
!   soil balance error warnings.  (DML, 4/8/2009)
!dir$ concurrent
!cdir nodep
    do p = begp,endp
       c = pptr%column(p)
       g = pptr%gridcell(p)
       l = pptr%landunit(p)
       if (lptr%itype(l) == istsoil) then
          if (wtpfttot2(c) .ne. 0) then
             pptr%wtgcell(p)   = (wtpfttot1(c)/wtpfttot2(c))*pptr%wtgcell(p)
             pptr%wtlunit(p)   = pptr%wtgcell(p) / lptr%wtgcell(l)
             pptr%wtcol(p)     = pptr%wtgcell(p) / lptr%wtgcell(l)
          end if
       end if

    end do
   
    deallocate(wtpfttot1,wtpfttot2) 
    
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
    integer  :: err, ierr, ret
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
       call ncd_iolocal(ncid, 'PCT_PFT', 'read', arrayl, grlnd, start, count, status=ret)
       if (ret /= 0) call endrun( trim(subname)//' ERROR: PCT_PFT NOT on pftdyn file' )
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
       write(iulog,*) subname,' error: sum(pct) over numpft+1 is not = 100.',sumerr,ierr,pctspec(ierr),pctpft(ierr,:)
       call endrun()
    else if (err == 2) then
       write(iulog,*)subname,' error: sum(pct) over numpft+1 is < 0.',sumerr,ierr,pctspec(ierr),pctpft(ierr,:)
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
    real(r8) :: dtime         ! land model time step (sec)
    real(r8) :: dwt           ! change in pft weight (relative to column)
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
          write(iulog,*)subname,' allocation error for loss_h2ocan'; call endrun()
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
  
#ifdef CN
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: pftdyn_cnbal
!
! !INTERFACE:
  subroutine pftdyn_cnbal()
!
! !DESCRIPTION:
! modify pft-level state and flux variables to maintain carbon and nitrogen balance with
! dynamic pft-weights.
!
! !USES:
    use shr_kind_mod, only : r8 => shr_kind_r8
    use shr_const_mod,only : SHR_CONST_PDB
    use decompMod   , only : get_proc_bounds
    use clm_varcon  , only : istsoil, c13ratio
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
    real(r8) :: dt            ! land model time step (sec)
    real(r8) :: init_h2ocan   ! initial canopy water mass
    real(r8) :: new_h2ocan    ! canopy water mass after weight shift
    real(r8), allocatable :: dwt_leafc_seed(:)       ! pft-level mass gain due to seeding of new area
    real(r8), allocatable :: dwt_leafn_seed(:)       ! pft-level mass gain due to seeding of new area
    real(r8), allocatable :: dwt_leafc13_seed(:)     ! pft-level mass gain due to seeding of new area
    real(r8), allocatable :: dwt_deadstemc_seed(:)       ! pft-level mass gain due to seeding of new area
    real(r8), allocatable :: dwt_deadstemn_seed(:)       ! pft-level mass gain due to seeding of new area
    real(r8), allocatable :: dwt_deadstemc13_seed(:)     ! pft-level mass gain due to seeding of new area
    real(r8), allocatable :: dwt_frootc_to_litter(:)       ! pft-level mass loss due to weight shift
    real(r8), allocatable :: dwt_livecrootc_to_litter(:)   ! pft-level mass loss due to weight shift
    real(r8), allocatable :: dwt_deadcrootc_to_litter(:)   ! pft-level mass loss due to weight shift
    real(r8), allocatable, target :: dwt_frootc13_to_litter(:)     ! pft-level mass loss due to weight shift
    real(r8), allocatable, target :: dwt_livecrootc13_to_litter(:) ! pft-level mass loss due to weight shift
    real(r8), allocatable, target :: dwt_deadcrootc13_to_litter(:) ! pft-level mass loss due to weight shift
    real(r8), allocatable, target :: dwt_frootn_to_litter(:)       ! pft-level mass loss due to weight shift
    real(r8), allocatable, target :: dwt_livecrootn_to_litter(:)   ! pft-level mass loss due to weight shift
    real(r8), allocatable, target :: dwt_deadcrootn_to_litter(:)   ! pft-level mass loss due to weight shift
    real(r8), allocatable :: conv_cflux(:)         ! pft-level mass loss due to weight shift
    real(r8), allocatable :: prod10_cflux(:)       ! pft-level mass loss due to weight shift
    real(r8), allocatable :: prod100_cflux(:)      ! pft-level mass loss due to weight shift
    real(r8), allocatable, target :: conv_c13flux(:)       ! pft-level mass loss due to weight shift
    real(r8), allocatable, target :: prod10_c13flux(:)     ! pft-level mass loss due to weight shift
    real(r8), allocatable, target :: prod100_c13flux(:)    ! pft-level mass loss due to weight shift
    real(r8), allocatable, target :: conv_nflux(:)         ! pft-level mass loss due to weight shift
    real(r8), allocatable, target :: prod10_nflux(:)       ! pft-level mass loss due to weight shift
    real(r8), allocatable, target :: prod100_nflux(:)      ! pft-level mass loss due to weight shift
    real(r8) :: c3_del13c     ! typical del13C for C3 photosynthesis (permil, relative to PDB)
    real(r8) :: c4_del13c     ! typical del13C for C4 photosynthesis (permil, relative to PDB)
    real(r8) :: c3_r1         ! isotope ratio (13c/12c) for C3 photosynthesis
    real(r8) :: c4_r1         ! isotope ratio (13c/12c) for C4 photosynthesis
    real(r8) :: c3_r2         ! isotope ratio (13c/[12c+13c]) for C3 photosynthesis
    real(r8) :: c4_r2         ! isotope ratio (13c/[12c+13c]) for C4 photosynthesis
    real(r8) :: kprod10       ! decay constant for 10-year product pool
    real(r8) :: kprod100      ! decay constant for 100-year product pool
    real(r8) :: t1,t2,wt_new,wt_old
    real(r8) :: init_state, change_state, new_state
	real(r8) :: tot_leaf, pleaf, pstor, pxfer
	real(r8) :: leafc_seed, leafn_seed, leafc13_seed
	real(r8) :: deadstemc_seed, deadstemn_seed, deadstemc13_seed
    real(r8), pointer :: dwt_ptr0, dwt_ptr1, dwt_ptr2, dwt_ptr3, ptr
    real(r8) :: pconv(0:16)     ! proportion of deadstem to conversion flux
    real(r8) :: pprod10(0:16)   ! proportion of deadstem to 10-yr product pool
    real(r8) :: pprod100(0:16)  ! proportion of deadstem to 100-yr product pool
    type(landunit_type), pointer :: lptr         ! pointer to landunit derived subtype
    type(column_type),   pointer :: cptr         ! pointer to column derived subtype
    type(pft_type)   ,   pointer :: pptr         ! pointer to pft derived subtype
    character(len=32) :: subname='pftdyn_cbal' ! subroutine name
!-----------------------------------------------------------------------
    
    ! set deadstem proportions
    ! veg type:      0       1       2       3       4       5       6       7       8       9      10      11      12     &
    !                13      14      15      16
    pconv   =    (/0.0_r8, 0.6_r8, 0.6_r8, 0.6_r8, 0.6_r8, 0.6_r8, 0.6_r8, 0.6_r8, 0.6_r8, 0.8_r8, 0.8_r8, 0.8_r8, 1.0_r8, &
                   1.0_r8, 1.0_r8, 1.0_r8, 1.0_r8/)
    pprod10 =    (/0.0_r8, 0.3_r8, 0.3_r8, 0.3_r8, 0.4_r8, 0.3_r8, 0.4_r8, 0.3_r8, 0.3_r8, 0.2_r8, 0.2_r8, 0.2_r8, 0.0_r8, &
                   0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8/)
    pprod100 =   (/0.0_r8, 0.1_r8, 0.1_r8, 0.1_r8, 0.0_r8, 0.1_r8, 0.0_r8, 0.1_r8, 0.1_r8, 0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, &
                   0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8/)
    
    ! Set pointers into derived type

    lptr => clm3%g%l
    cptr => clm3%g%l%c
    pptr => clm3%g%l%c%p

    ! Get relevant sizes

    call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)

    ! Allocate pft-level mass loss arrays
    allocate(dwt_leafc_seed(begp:endp), stat=ier)
    if (ier /= 0) then
          write(iulog,*)subname,' allocation error for dwt_leafc_seed'; call endrun()
    end if
    allocate(dwt_leafn_seed(begp:endp), stat=ier)
    if (ier /= 0) then
          write(iulog,*)subname,' allocation error for dwt_leafn_seed'; call endrun()
    end if
    allocate(dwt_leafc13_seed(begp:endp), stat=ier)
    if (ier /= 0) then
          write(iulog,*)subname,' allocation error for dwt_leafc13_seed'; call endrun()
    end if
    allocate(dwt_deadstemc_seed(begp:endp), stat=ier)
    if (ier /= 0) then
          write(iulog,*)subname,' allocation error for dwt_deadstemc_seed'; call endrun()
    end if
    allocate(dwt_deadstemn_seed(begp:endp), stat=ier)
    if (ier /= 0) then
          write(iulog,*)subname,' allocation error for dwt_deadstemn_seed'; call endrun()
    end if
    allocate(dwt_deadstemc13_seed(begp:endp), stat=ier)
    if (ier /= 0) then
          write(iulog,*)subname,' allocation error for dwt_deadstemc13_seed'; call endrun()
    end if
    allocate(dwt_frootc_to_litter(begp:endp), stat=ier)
    if (ier /= 0) then
          write(iulog,*)subname,' allocation error for dwt_frootc_to_litter'; call endrun()
    end if
    allocate(dwt_livecrootc_to_litter(begp:endp), stat=ier)
    if (ier /= 0) then
          write(iulog,*)subname,' allocation error for dwt_livecrootc_to_litter'; call endrun()
    end if
    allocate(dwt_deadcrootc_to_litter(begp:endp), stat=ier)
    if (ier /= 0) then
          write(iulog,*)subname,' allocation error for dwt_deadcrootc_to_litter'; call endrun()
    end if
    allocate(dwt_frootc13_to_litter(begp:endp), stat=ier)
    if (ier /= 0) then
          write(iulog,*)subname,' allocation error for dwt_frootc13_to_litter'; call endrun()
    end if
    allocate(dwt_livecrootc13_to_litter(begp:endp), stat=ier)
    if (ier /= 0) then
          write(iulog,*)subname,' allocation error for dwt_livecrootc13_to_litter'; call endrun()
    end if
    allocate(dwt_deadcrootc13_to_litter(begp:endp), stat=ier)
    if (ier /= 0) then
          write(iulog,*)subname,' allocation error for dwt_deadcrootc13_to_litter'; call endrun()
    end if
    allocate(dwt_frootn_to_litter(begp:endp), stat=ier)
    if (ier /= 0) then
          write(iulog,*)subname,' allocation error for dwt_frootn_to_litter'; call endrun()
    end if
    allocate(dwt_livecrootn_to_litter(begp:endp), stat=ier)
    if (ier /= 0) then
          write(iulog,*)subname,' allocation error for dwt_livecrootn_to_litter'; call endrun()
    end if
    allocate(dwt_deadcrootn_to_litter(begp:endp), stat=ier)
    if (ier /= 0) then
          write(iulog,*)subname,' allocation error for dwt_deadcrootn_to_litter'; call endrun()
    end if
    allocate(conv_cflux(begp:endp), stat=ier)
    if (ier /= 0) then
          write(iulog,*)subname,' allocation error for conv_cflux'; call endrun()
    end if
    allocate(prod10_cflux(begp:endp), stat=ier)
    if (ier /= 0) then
          write(iulog,*)subname,' allocation error for prod10_cflux'; call endrun()
    end if
    allocate(prod100_cflux(begp:endp), stat=ier)
    if (ier /= 0) then
          write(iulog,*)subname,' allocation error for prod100_cflux'; call endrun()
    end if
    allocate(conv_c13flux(begp:endp), stat=ier)
    if (ier /= 0) then
          write(iulog,*)subname,' allocation error for conv_c13flux'; call endrun()
    end if
    allocate(prod10_c13flux(begp:endp), stat=ier)
    if (ier /= 0) then
          write(iulog,*)subname,' allocation error for prod10_c13flux'; call endrun()
    end if
    allocate(prod100_c13flux(begp:endp), stat=ier)
    if (ier /= 0) then
          write(iulog,*)subname,' allocation error for prod100_c13flux'; call endrun()
    end if
    allocate(conv_nflux(begp:endp), stat=ier)
    if (ier /= 0) then
          write(iulog,*)subname,' allocation error for conv_nflux'; call endrun()
    end if
    allocate(prod10_nflux(begp:endp), stat=ier)
    if (ier /= 0) then
          write(iulog,*)subname,' allocation error for prod10_nflux'; call endrun()
    end if
    allocate(prod100_nflux(begp:endp), stat=ier)
    if (ier /= 0) then
          write(iulog,*)subname,' allocation error for prod100_nflux'; call endrun()
    end if

    ! Get time step
    dt = real( get_step_size(), r8 )

!dir$ concurrent
!cdir nodep
	do p = begp,endp
		! initialize all the pft-level local flux arrays
		dwt_leafc_seed(p) = 0._r8
		dwt_leafn_seed(p) = 0._r8
		dwt_leafc13_seed(p) = 0._r8
		dwt_deadstemc_seed(p) = 0._r8
		dwt_deadstemn_seed(p) = 0._r8
		dwt_deadstemc13_seed(p) = 0._r8
		dwt_frootc_to_litter(p) = 0._r8
		dwt_livecrootc_to_litter(p) = 0._r8
		dwt_deadcrootc_to_litter(p) = 0._r8
		dwt_frootc13_to_litter(p) = 0._r8
		dwt_livecrootc13_to_litter(p) = 0._r8
		dwt_deadcrootc13_to_litter(p) = 0._r8
		dwt_frootn_to_litter(p) = 0._r8
		dwt_livecrootn_to_litter(p) = 0._r8
		dwt_deadcrootn_to_litter(p) = 0._r8
		conv_cflux(p) = 0._r8
		prod10_cflux(p) = 0._r8
		prod100_cflux(p) = 0._r8
		conv_c13flux(p) = 0._r8
		prod10_c13flux(p) = 0._r8
		prod100_c13flux(p) = 0._r8
		conv_nflux(p) = 0._r8
		prod10_nflux(p) = 0._r8
		prod100_nflux(p) = 0._r8
       
		l = pptr%landunit(p)
		c = pptr%column(p)
		if (lptr%itype(l) == istsoil) then

			! calculate the change in weight for the timestep
			dwt = pptr%wtcol(p)-wtcol_old(p)

			! PFTs for which weight increases on this timestep
			if (dwt > 0._r8) then

				! first identify PFTs that are initiating on this timestep
				! and set all the necessary state and flux variables
				if (wtcol_old(p) == 0._r8) then

					! set initial conditions for PFT that is being initiated
					! in this time step.  Based on the settings in cnIniTimeVar.

					! pft-level carbon state variables
					pptr%pcs%leafc(p)              = 0._r8
					pptr%pcs%leafc_storage(p)      = 0._r8
					pptr%pcs%leafc_xfer(p)         = 0._r8
					pptr%pcs%frootc(p)             = 0._r8
					pptr%pcs%frootc_storage(p)     = 0._r8
					pptr%pcs%frootc_xfer(p)        = 0._r8
					pptr%pcs%livestemc(p)          = 0._r8
					pptr%pcs%livestemc_storage(p)  = 0._r8
					pptr%pcs%livestemc_xfer(p)     = 0._r8
					pptr%pcs%deadstemc(p)          = 0._r8
					pptr%pcs%deadstemc_storage(p)  = 0._r8
					pptr%pcs%deadstemc_xfer(p)     = 0._r8
					pptr%pcs%livecrootc(p)         = 0._r8
					pptr%pcs%livecrootc_storage(p) = 0._r8
					pptr%pcs%livecrootc_xfer(p)    = 0._r8
					pptr%pcs%deadcrootc(p)         = 0._r8
					pptr%pcs%deadcrootc_storage(p) = 0._r8
					pptr%pcs%deadcrootc_xfer(p)    = 0._r8
					pptr%pcs%gresp_storage(p)      = 0._r8
					pptr%pcs%gresp_xfer(p)         = 0._r8
					pptr%pcs%cpool(p)              = 0._r8
					pptr%pcs%xsmrpool(p)           = 0._r8
					pptr%pcs%pft_ctrunc(p)         = 0._r8
					pptr%pcs%dispvegc(p)           = 0._r8
					pptr%pcs%storvegc(p)           = 0._r8
					pptr%pcs%totvegc(p)            = 0._r8
					pptr%pcs%totpftc(p)            = 0._r8

					! pft-level carbon-13 state variables
					pptr%pc13s%leafc(p)              = 0._r8
					pptr%pc13s%leafc_storage(p)      = 0._r8
					pptr%pc13s%leafc_xfer(p)         = 0._r8
					pptr%pc13s%frootc(p)             = 0._r8
					pptr%pc13s%frootc_storage(p)     = 0._r8
					pptr%pc13s%frootc_xfer(p)        = 0._r8
					pptr%pc13s%livestemc(p)          = 0._r8
					pptr%pc13s%livestemc_storage(p)  = 0._r8
					pptr%pc13s%livestemc_xfer(p)     = 0._r8
					pptr%pc13s%deadstemc(p)          = 0._r8
					pptr%pc13s%deadstemc_storage(p)  = 0._r8
					pptr%pc13s%deadstemc_xfer(p)     = 0._r8
					pptr%pc13s%livecrootc(p)         = 0._r8
					pptr%pc13s%livecrootc_storage(p) = 0._r8
					pptr%pc13s%livecrootc_xfer(p)    = 0._r8
					pptr%pc13s%deadcrootc(p)         = 0._r8
					pptr%pc13s%deadcrootc_storage(p) = 0._r8
					pptr%pc13s%deadcrootc_xfer(p)    = 0._r8
					pptr%pc13s%gresp_storage(p)      = 0._r8
					pptr%pc13s%gresp_xfer(p)         = 0._r8
					pptr%pc13s%cpool(p)              = 0._r8
					pptr%pc13s%xsmrpool(p)           = 0._r8
					pptr%pc13s%pft_ctrunc(p)         = 0._r8
					pptr%pc13s%dispvegc(p)           = 0._r8
					pptr%pc13s%storvegc(p)           = 0._r8
					pptr%pc13s%totvegc(p)            = 0._r8
					pptr%pc13s%totpftc(p)            = 0._r8

					! pft-level nitrogen state variables
					pptr%pns%leafn(p)	           = 0._r8
					pptr%pns%leafn_storage(p)      = 0._r8
					pptr%pns%leafn_xfer(p)         = 0._r8
					pptr%pns%frootn(p)	           = 0._r8
					pptr%pns%frootn_storage(p)     = 0._r8
					pptr%pns%frootn_xfer(p)        = 0._r8
					pptr%pns%livestemn(p)	       = 0._r8
					pptr%pns%livestemn_storage(p)  = 0._r8
					pptr%pns%livestemn_xfer(p)     = 0._r8
					pptr%pns%deadstemn(p)	       = 0._r8
					pptr%pns%deadstemn_storage(p)  = 0._r8
					pptr%pns%deadstemn_xfer(p)     = 0._r8
					pptr%pns%livecrootn(p)         = 0._r8
					pptr%pns%livecrootn_storage(p) = 0._r8
					pptr%pns%livecrootn_xfer(p)    = 0._r8
					pptr%pns%deadcrootn(p)         = 0._r8
					pptr%pns%deadcrootn_storage(p) = 0._r8
					pptr%pns%deadcrootn_xfer(p)    = 0._r8
					pptr%pns%retransn(p)	       = 0._r8
					pptr%pns%npool(p)	           = 0._r8
					pptr%pns%pft_ntrunc(p)         = 0._r8
					pptr%pns%dispvegn(p)           = 0._r8
					pptr%pns%storvegn(p)           = 0._r8
					pptr%pns%totvegn(p)            = 0._r8
					pptr%pns%totpftn (p)           = 0._r8

					! initialize same flux and epv variables that are set
					! in CNiniTimeVar
					pptr%pcf%psnsun(p) = 0._r8
					pptr%pcf%psnsha(p) = 0._r8
					pptr%pc13f%psnsun(p) = 0._r8
					pptr%pc13f%psnsha(p) = 0._r8
					pptr%pps%laisun(p) = 0._r8
					pptr%pps%laisha(p) = 0._r8
					pptr%pps%lncsun(p) = 0._r8
					pptr%pps%lncsha(p) = 0._r8
					pptr%pps%vcmxsun(p) = 0._r8
					pptr%pps%vcmxsha(p) = 0._r8
					pptr%pps%alphapsnsun(p) = 0._r8
					pptr%pps%alphapsnsha(p) = 0._r8

					pptr%pepv%dormant_flag(p) = 1._r8
					pptr%pepv%days_active(p) = 0._r8
					pptr%pepv%onset_flag(p) = 0._r8
					pptr%pepv%onset_counter(p) = 0._r8
					pptr%pepv%onset_gddflag(p) = 0._r8
					pptr%pepv%onset_fdd(p) = 0._r8
					pptr%pepv%onset_gdd(p) = 0._r8
					pptr%pepv%onset_swi(p) = 0.0_r8
					pptr%pepv%offset_flag(p) = 0._r8
					pptr%pepv%offset_counter(p) = 0._r8
					pptr%pepv%offset_fdd(p) = 0._r8
					pptr%pepv%offset_swi(p) = 0._r8
					pptr%pepv%lgsf(p) = 0._r8
					pptr%pepv%bglfr(p) = 0._r8
					pptr%pepv%bgtr(p) = 0._r8
					! difference from CNiniTimeVar: using column-level
					! information to initialize annavg_t2m.
					pptr%pepv%annavg_t2m(p) = cptr%cps%cannavg_t2m(c)
					pptr%pepv%tempavg_t2m(p) = 0._r8
					pptr%pepv%gpp(p) = 0._r8
					pptr%pepv%availc(p) = 0._r8
					pptr%pepv%xsmrpool_recover(p) = 0._r8
					pptr%pepv%xsmrpool_c13ratio(p) = c13ratio
					pptr%pepv%alloc_pnow(p) = 1._r8
					pptr%pepv%c_allometry(p) = 0._r8
					pptr%pepv%n_allometry(p) = 0._r8
					pptr%pepv%plant_ndemand(p) = 0._r8
					pptr%pepv%tempsum_potential_gpp(p) = 0._r8
					pptr%pepv%annsum_potential_gpp(p) = 0._r8
					pptr%pepv%tempmax_retransn(p) = 0._r8
					pptr%pepv%annmax_retransn(p) = 0._r8
					pptr%pepv%avail_retransn(p) = 0._r8
					pptr%pepv%plant_nalloc(p) = 0._r8
					pptr%pepv%plant_calloc(p) = 0._r8
					pptr%pepv%excess_cflux(p) = 0._r8
					pptr%pepv%downreg(p) = 0._r8
					pptr%pepv%prev_leafc_to_litter(p) = 0._r8
					pptr%pepv%prev_frootc_to_litter(p) = 0._r8
					pptr%pepv%tempsum_npp(p) = 0._r8
					pptr%pepv%annsum_npp(p) = 0._r8
					pptr%pepv%rc13_canair(p) = 0._r8
					pptr%pepv%rc13_psnsun(p) = 0._r8
					pptr%pepv%rc13_psnsha(p) = 0._r8

				end if  ! end initialization of new pft

				! (still in dwt > 0 block)

				! set the seed sources for leaf and deadstem
				! leaf source is split later between leaf, leaf_storage, leaf_xfer
				leafc_seed   = 0._r8
				leafn_seed   = 0._r8
				leafc13_seed = 0._r8
				deadstemc_seed   = 0._r8
				deadstemn_seed   = 0._r8
				deadstemc13_seed = 0._r8
				if (pptr%itype(p) /= 0) then
					leafc_seed = 1._r8
					leafn_seed  = leafc_seed / pftcon%leafcn(pptr%itype(p))
					if (pftcon%woody(pptr%itype(p)) == 1._r8) then
						deadstemc_seed = 0.1_r8
						deadstemn_seed = deadstemc_seed / pftcon%deadwdcn(pptr%itype(p))
					end if

					! 13c state is initialized assuming del13c = -28 permil for C3, and -13 permil for C4.
					! That translates to ratios of (13c/(12c+13c)) of 0.01080455 for C3, and 0.01096945 for C4
					! based on the following formulae: 
					! r1 (13/12) = PDB + (del13c * PDB)/1000.0
					! r2 (13/(13+12)) = r1/(1+r1)
					! PDB = 0.0112372_R8  (ratio of 13C/12C in Pee Dee Belemnite, C isotope standard)
					c3_del13c = -28._r8
					c4_del13c = -13._r8
					c3_r1 = SHR_CONST_PDB + ((c3_del13c*SHR_CONST_PDB)/1000._r8)
					c3_r2 = c3_r1/(1._r8 + c3_r1)
					c4_r1 = SHR_CONST_PDB + ((c4_del13c*SHR_CONST_PDB)/1000._r8)
					c4_r2 = c4_r1/(1._r8 + c4_r1)

					if (pftcon%c3psn(pptr%itype(p)) == 1._r8) then
						leafc13_seed     = leafc_seed     * c3_r2
						deadstemc13_seed = deadstemc_seed * c3_r2
					else
						leafc13_seed     = leafc_seed     * c4_r2
						deadstemc13_seed = deadstemc_seed * c4_r2
					end if 
				end if

				! When PFT area expands (dwt > 0), the pft-level mass density 
				! is modified to conserve the original pft mass distributed
				! over the new (larger) area, plus a term to account for the 
				! introduction of new seed source for leaf and deadstem
				t1 = wtcol_old(p)/pptr%wtcol(p)
				t2 = dwt/pptr%wtcol(p)

				tot_leaf = pptr%pcs%leafc(p) + pptr%pcs%leafc_storage(p) + pptr%pcs%leafc_xfer(p)
				pleaf = 0._r8
				pstor = 0._r8
				pxfer = 0._r8
				if (tot_leaf /= 0._r8) then
					! when adding seed source to non-zero leaf state, use current proportions
					pleaf = pptr%pcs%leafc(p)/tot_leaf
					pstor = pptr%pcs%leafc_storage(p)/tot_leaf
					pxfer = pptr%pcs%leafc_xfer(p)/tot_leaf
				else
					! when initiating from zero leaf state, use evergreen flag to set proportions
					if (pftcon%evergreen(pptr%itype(p)) == 1._r8) then
						pleaf = 1._r8
					else
						pstor = 1._r8
					end if
				end if 
				pptr%pcs%leafc(p)         = pptr%pcs%leafc(p)*t1         + leafc_seed*pleaf*t2
				pptr%pcs%leafc_storage(p) = pptr%pcs%leafc_storage(p)*t1 + leafc_seed*pstor*t2
				pptr%pcs%leafc_xfer(p)    = pptr%pcs%leafc_xfer(p)*t1    + leafc_seed*pxfer*t2
				pptr%pcs%frootc(p)  		   = pptr%pcs%frootc(p) 			* t1
				pptr%pcs%frootc_storage(p)     = pptr%pcs%frootc_storage(p) 	* t1
				pptr%pcs%frootc_xfer(p) 	   = pptr%pcs%frootc_xfer(p)		* t1
				pptr%pcs%livestemc(p)		   = pptr%pcs%livestemc(p)  		* t1
				pptr%pcs%livestemc_storage(p)  = pptr%pcs%livestemc_storage(p)  * t1
				pptr%pcs%livestemc_xfer(p)     = pptr%pcs%livestemc_xfer(p) 	* t1
				pptr%pcs%deadstemc(p)     = pptr%pcs%deadstemc(p)*t1     + deadstemc_seed*t2
				pptr%pcs%deadstemc_storage(p)  = pptr%pcs%deadstemc_storage(p)  * t1
				pptr%pcs%deadstemc_xfer(p)     = pptr%pcs%deadstemc_xfer(p) 	* t1
				pptr%pcs%livecrootc(p)  	   = pptr%pcs%livecrootc(p) 		* t1
				pptr%pcs%livecrootc_storage(p) = pptr%pcs%livecrootc_storage(p) * t1
				pptr%pcs%livecrootc_xfer(p)    = pptr%pcs%livecrootc_xfer(p)	* t1
				pptr%pcs%deadcrootc(p)  	   = pptr%pcs%deadcrootc(p) 		* t1
				pptr%pcs%deadcrootc_storage(p) = pptr%pcs%deadcrootc_storage(p) * t1
				pptr%pcs%deadcrootc_xfer(p)    = pptr%pcs%deadcrootc_xfer(p)	* t1
				pptr%pcs%gresp_storage(p)	   = pptr%pcs%gresp_storage(p)  	* t1
				pptr%pcs%gresp_xfer(p)  	   = pptr%pcs%gresp_xfer(p) 		* t1
				pptr%pcs%cpool(p)			   = pptr%pcs%cpool(p)  			* t1
				pptr%pcs%xsmrpool(p)		   = pptr%pcs%xsmrpool(p)			* t1
				pptr%pcs%pft_ctrunc(p)  	   = pptr%pcs%pft_ctrunc(p) 		* t1
				pptr%pcs%dispvegc(p)		   = pptr%pcs%dispvegc(p)			* t1
				pptr%pcs%storvegc(p)		   = pptr%pcs%storvegc(p)			* t1
				pptr%pcs%totvegc(p) 		   = pptr%pcs%totvegc(p)			* t1
				pptr%pcs%totpftc(p) 		   = pptr%pcs%totpftc(p)			* t1

				! pft-level carbon-13 state variables 
				tot_leaf = pptr%pc13s%leafc(p) + pptr%pc13s%leafc_storage(p) + pptr%pc13s%leafc_xfer(p)
				pleaf = 0._r8
				pstor = 0._r8
				pxfer = 0._r8
				if (tot_leaf /= 0._r8) then
					pleaf = pptr%pc13s%leafc(p)/tot_leaf
					pstor = pptr%pc13s%leafc_storage(p)/tot_leaf
					pxfer = pptr%pc13s%leafc_xfer(p)/tot_leaf
				else
					! when initiating from zero leaf state, use evergreen flag to set proportions
					if (pftcon%evergreen(pptr%itype(p)) == 1._r8) then
						pleaf = 1._r8
					else
						pstor = 1._r8
					end if
				end if 
				pptr%pc13s%leafc(p)         = pptr%pc13s%leafc(p)*t1         + leafc13_seed*pleaf*t2
				pptr%pc13s%leafc_storage(p) = pptr%pc13s%leafc_storage(p)*t1 + leafc13_seed*pstor*t2
				pptr%pc13s%leafc_xfer(p)    = pptr%pc13s%leafc_xfer(p)*t1    + leafc13_seed*pxfer*t2
				pptr%pc13s%frootc(p)			 = pptr%pc13s%frootc(p) 			* t1
				pptr%pc13s%frootc_storage(p)	 = pptr%pc13s%frootc_storage(p) 	* t1
				pptr%pc13s%frootc_xfer(p)		 = pptr%pc13s%frootc_xfer(p)		* t1
				pptr%pc13s%livestemc(p) 		 = pptr%pc13s%livestemc(p)  		* t1
				pptr%pc13s%livestemc_storage(p)  = pptr%pc13s%livestemc_storage(p)  * t1
				pptr%pc13s%livestemc_xfer(p)	 = pptr%pc13s%livestemc_xfer(p) 	* t1
				pptr%pc13s%deadstemc(p)     = pptr%pc13s%deadstemc(p)*t1     + deadstemc13_seed*t2
				pptr%pc13s%deadstemc_storage(p)  = pptr%pc13s%deadstemc_storage(p)  * t1
				pptr%pc13s%deadstemc_xfer(p)	 = pptr%pc13s%deadstemc_xfer(p) 	* t1
				pptr%pc13s%livecrootc(p)		 = pptr%pc13s%livecrootc(p) 		* t1
				pptr%pc13s%livecrootc_storage(p) = pptr%pc13s%livecrootc_storage(p) * t1
				pptr%pc13s%livecrootc_xfer(p)	 = pptr%pc13s%livecrootc_xfer(p)	* t1
				pptr%pc13s%deadcrootc(p)		 = pptr%pc13s%deadcrootc(p) 		* t1
				pptr%pc13s%deadcrootc_storage(p) = pptr%pc13s%deadcrootc_storage(p) * t1
				pptr%pc13s%deadcrootc_xfer(p)	 = pptr%pc13s%deadcrootc_xfer(p)	* t1
				pptr%pc13s%gresp_storage(p) 	 = pptr%pc13s%gresp_storage(p)  	* t1
				pptr%pc13s%gresp_xfer(p)		 = pptr%pc13s%gresp_xfer(p) 		* t1
				pptr%pc13s%cpool(p) 			 = pptr%pc13s%cpool(p)  			* t1
				pptr%pc13s%xsmrpool(p)  		 = pptr%pc13s%xsmrpool(p)			* t1
				pptr%pc13s%pft_ctrunc(p)		 = pptr%pc13s%pft_ctrunc(p) 		* t1
				pptr%pc13s%dispvegc(p)  		 = pptr%pc13s%dispvegc(p)			* t1
				pptr%pc13s%storvegc(p)  		 = pptr%pc13s%storvegc(p)			* t1
				pptr%pc13s%totvegc(p)			 = pptr%pc13s%totvegc(p)			* t1
				pptr%pc13s%totpftc(p)			 = pptr%pc13s%totpftc(p)			* t1

				tot_leaf = pptr%pns%leafn(p) + pptr%pns%leafn_storage(p) + pptr%pns%leafn_xfer(p)
				pleaf = 0._r8
				pstor = 0._r8
				pxfer = 0._r8
				if (tot_leaf /= 0._r8) then
					pleaf = pptr%pns%leafn(p)/tot_leaf
					pstor = pptr%pns%leafn_storage(p)/tot_leaf
					pxfer = pptr%pns%leafn_xfer(p)/tot_leaf
				else
					! when initiating from zero leaf state, use evergreen flag to set proportions
					if (pftcon%evergreen(pptr%itype(p)) == 1._r8) then
						pleaf = 1._r8
					else
						pstor = 1._r8
					end if
				end if 
				! pft-level nitrogen state variables
				pptr%pns%leafn(p)         = pptr%pns%leafn(p)*t1         + leafn_seed*pleaf*t2
				pptr%pns%leafn_storage(p) = pptr%pns%leafn_storage(p)*t1 + leafn_seed*pstor*t2
				pptr%pns%leafn_xfer(p)    = pptr%pns%leafn_xfer(p)*t1    + leafn_seed*pxfer*t2
				pptr%pns%frootn(p)  		   = pptr%pns%frootn(p) 			* t1
				pptr%pns%frootn_storage(p)     = pptr%pns%frootn_storage(p) 	* t1
				pptr%pns%frootn_xfer(p) 	   = pptr%pns%frootn_xfer(p)		* t1
				pptr%pns%livestemn(p)		   = pptr%pns%livestemn(p)  		* t1
				pptr%pns%livestemn_storage(p)  = pptr%pns%livestemn_storage(p)  * t1
				pptr%pns%livestemn_xfer(p)     = pptr%pns%livestemn_xfer(p) 	* t1
				pptr%pns%deadstemn(p)     = pptr%pns%deadstemn(p)*t1     + deadstemn_seed*t2
				pptr%pns%deadstemn_storage(p)  = pptr%pns%deadstemn_storage(p)  * t1
				pptr%pns%deadstemn_xfer(p)     = pptr%pns%deadstemn_xfer(p) 	* t1
				pptr%pns%livecrootn(p)  	   = pptr%pns%livecrootn(p) 		* t1
				pptr%pns%livecrootn_storage(p) = pptr%pns%livecrootn_storage(p) * t1
				pptr%pns%livecrootn_xfer(p)    = pptr%pns%livecrootn_xfer(p)	* t1
				pptr%pns%deadcrootn(p)  	   = pptr%pns%deadcrootn(p) 		* t1
				pptr%pns%deadcrootn_storage(p) = pptr%pns%deadcrootn_storage(p) * t1
				pptr%pns%deadcrootn_xfer(p)    = pptr%pns%deadcrootn_xfer(p)	* t1
				pptr%pns%retransn(p)		   = pptr%pns%retransn(p)			* t1
				pptr%pns%npool(p)			   = pptr%pns%npool(p)  			* t1
				pptr%pns%pft_ntrunc(p)  	   = pptr%pns%pft_ntrunc(p) 		* t1
				pptr%pns%dispvegn(p)		   = pptr%pns%dispvegn(p)			* t1
				pptr%pns%storvegn(p)		   = pptr%pns%storvegn(p)			* t1
				pptr%pns%totvegn(p) 		   = pptr%pns%totvegn(p)			* t1
				pptr%pns%totpftn(p) 		   = pptr%pns%totpftn(p)			* t1

				! update temporary seed source arrays
				! These are calculated in terms of the required contributions from
				! column-level seed source
				dwt_leafc_seed(p)   = leafc_seed   * dwt
				dwt_leafc13_seed(p) = leafc13_seed * dwt
				dwt_leafn_seed(p)   = leafn_seed   * dwt
				dwt_deadstemc_seed(p)   = deadstemc_seed   * dwt
				dwt_deadstemc13_seed(p) = deadstemc13_seed * dwt
				dwt_deadstemn_seed(p)   = deadstemn_seed   * dwt

			else if (dwt < 0._r8) then

				! if the pft lost weight on the timestep, then the carbon and nitrogen state
				! variables are directed to litter, CWD, and wood product pools.

				! N.B. : the conv_cflux, prod10_cflux, and prod100_cflux fluxes are accumulated
				! as negative values, but the fluxes for pft-to-litter are accumulated as 
				! positive values

				! set local weight variables for this pft
				wt_new = pptr%wtcol(p)
				wt_old = wtcol_old(p)

				!---------------
				! C state update
				!---------------

				! leafc 
				ptr => pptr%pcs%leafc(p)
				init_state = ptr*wt_old
				change_state = ptr*dwt
				new_state = init_state+change_state
				if (wt_new /= 0._r8) then
					ptr = new_state/wt_new
					conv_cflux(p) = conv_cflux(p) + change_state
				else
					ptr = 0._r8
					conv_cflux(p) = conv_cflux(p) - init_state
				end if

				! leafc_storage 
				ptr => pptr%pcs%leafc_storage(p)
				init_state = ptr*wt_old
				change_state = ptr*dwt
				new_state = init_state+change_state
				if (wt_new /= 0._r8) then
					ptr = new_state/wt_new
					conv_cflux(p) = conv_cflux(p) + change_state
				else
					ptr = 0._r8
					conv_cflux(p) = conv_cflux(p) - init_state
				end if

				! leafc_xfer 
				ptr => pptr%pcs%leafc_xfer(p)
				init_state = ptr*wt_old
				change_state = ptr*dwt
				new_state = init_state+change_state
				if (wt_new /= 0._r8) then
					ptr = new_state/wt_new
					conv_cflux(p) = conv_cflux(p) + change_state
				else
					ptr = 0._r8
					conv_cflux(p) = conv_cflux(p) - init_state
				end if

				! frootc 
				ptr => pptr%pcs%frootc(p)
				init_state = ptr*wt_old
				change_state = ptr*dwt
				new_state = init_state+change_state
				if (wt_new /= 0._r8) then
					ptr = new_state/wt_new
					dwt_frootc_to_litter(p) = dwt_frootc_to_litter(p) - change_state
				else
					ptr = 0._r8
					dwt_frootc_to_litter(p) = dwt_frootc_to_litter(p) + init_state
				end if

				! frootc_storage 
				ptr => pptr%pcs%frootc_storage(p)
				init_state = ptr*wt_old
				change_state = ptr*dwt
				new_state = init_state+change_state
				if (wt_new /= 0._r8) then
					ptr = new_state/wt_new
					conv_cflux(p) = conv_cflux(p) + change_state
				else
					ptr = 0._r8
					conv_cflux(p) = conv_cflux(p) - init_state
				end if

				! frootc_xfer 
				ptr => pptr%pcs%frootc_xfer(p)
				init_state = ptr*wt_old
				change_state = ptr*dwt
				new_state = init_state+change_state
				if (wt_new /= 0._r8) then
					ptr = new_state/wt_new
					conv_cflux(p) = conv_cflux(p) + change_state
				else
					ptr = 0._r8
					conv_cflux(p) = conv_cflux(p) - init_state
				end if

				! livestemc 
				ptr => pptr%pcs%livestemc(p)
				init_state = ptr*wt_old
				change_state = ptr*dwt
				new_state = init_state+change_state
				if (wt_new /= 0._r8) then
					ptr = new_state/wt_new
					conv_cflux(p) = conv_cflux(p) + change_state
				else
					ptr = 0._r8
					conv_cflux(p) = conv_cflux(p) - init_state
				end if

				! livestemc_storage 
				ptr => pptr%pcs%livestemc_storage(p)
				init_state = ptr*wt_old
				change_state = ptr*dwt
				new_state = init_state+change_state
				if (wt_new /= 0._r8) then
					ptr = new_state/wt_new
					conv_cflux(p) = conv_cflux(p) + change_state
				else
					ptr = 0._r8
					conv_cflux(p) = conv_cflux(p) - init_state
				end if

				! livestemc_xfer 
				ptr => pptr%pcs%livestemc_xfer(p)
				init_state = ptr*wt_old
				change_state = ptr*dwt
				new_state = init_state+change_state
				if (wt_new /= 0._r8) then
					ptr = new_state/wt_new
					conv_cflux(p) = conv_cflux(p) + change_state
				else
					ptr = 0._r8
					conv_cflux(p) = conv_cflux(p) - init_state
				end if

				! deadstemc 
				ptr => pptr%pcs%deadstemc(p)
				init_state = ptr*wt_old
				change_state = ptr*dwt
				new_state = init_state+change_state
				if (wt_new /= 0._r8) then
					ptr = new_state/wt_new
					conv_cflux(p) = conv_cflux(p) + change_state*pconv(pptr%itype(p))
					prod10_cflux(p) = prod10_cflux(p) + change_state*pprod10(pptr%itype(p))
					prod100_cflux(p) = prod100_cflux(p) + change_state*pprod100(pptr%itype(p))
				else
					ptr = 0._r8
					conv_cflux(p) = conv_cflux(p) - init_state*pconv(pptr%itype(p))
					prod10_cflux(p) = prod10_cflux(p) - init_state*pprod10(pptr%itype(p))
					prod100_cflux(p) = prod100_cflux(p) - init_state*pprod100(pptr%itype(p))
				end if

				! deadstemc_storage 
				ptr => pptr%pcs%deadstemc_storage(p)
				init_state = ptr*wt_old
				change_state = ptr*dwt
				new_state = init_state+change_state
				if (wt_new /= 0._r8) then
					ptr = new_state/wt_new
					conv_cflux(p) = conv_cflux(p) + change_state
				else
					ptr = 0._r8
					conv_cflux(p) = conv_cflux(p) - init_state
				end if

				! deadstemc_xfer 
				ptr => pptr%pcs%deadstemc_xfer(p)
				init_state = ptr*wt_old
				change_state = ptr*dwt
				new_state = init_state+change_state
				if (wt_new /= 0._r8) then
					ptr = new_state/wt_new
					conv_cflux(p) = conv_cflux(p) + change_state
				else
					ptr = 0._r8
					conv_cflux(p) = conv_cflux(p) - init_state
				end if

				! livecrootc 
				ptr => pptr%pcs%livecrootc(p)
				init_state = ptr*wt_old
				change_state = ptr*dwt
				new_state = init_state+change_state
				if (wt_new /= 0._r8) then
					ptr = new_state/wt_new
					dwt_livecrootc_to_litter(p) = dwt_livecrootc_to_litter(p) - change_state
				else
					ptr = 0._r8
					dwt_livecrootc_to_litter(p) = dwt_livecrootc_to_litter(p) + init_state
				end if

				! livecrootc_storage 
				ptr => pptr%pcs%livecrootc_storage(p)
				init_state = ptr*wt_old
				change_state = ptr*dwt
				new_state = init_state+change_state
				if (wt_new /= 0._r8) then
					ptr = new_state/wt_new
					conv_cflux(p) = conv_cflux(p) + change_state
				else
					ptr = 0._r8
					conv_cflux(p) = conv_cflux(p) - init_state
				end if

				! livecrootc_xfer 
				ptr => pptr%pcs%livecrootc_xfer(p)
				init_state = ptr*wt_old
				change_state = ptr*dwt
				new_state = init_state+change_state
				if (wt_new /= 0._r8) then
					ptr = new_state/wt_new
					conv_cflux(p) = conv_cflux(p) + change_state
				else
					ptr = 0._r8
					conv_cflux(p) = conv_cflux(p) - init_state
				end if

				! deadcrootc 
				ptr => pptr%pcs%deadcrootc(p)
				init_state = ptr*wt_old
				change_state = ptr*dwt
				new_state = init_state+change_state
				if (wt_new /= 0._r8) then
					ptr = new_state/wt_new
					dwt_deadcrootc_to_litter(p) = dwt_deadcrootc_to_litter(p) - change_state
				else
					ptr = 0._r8
					dwt_deadcrootc_to_litter(p) = dwt_deadcrootc_to_litter(p) + init_state
				end if

				! deadcrootc_storage 
				ptr => pptr%pcs%deadcrootc_storage(p)
				init_state = ptr*wt_old
				change_state = ptr*dwt
				new_state = init_state+change_state
				if (wt_new /= 0._r8) then
					ptr = new_state/wt_new
					conv_cflux(p) = conv_cflux(p) + change_state
				else
					ptr = 0._r8
					conv_cflux(p) = conv_cflux(p) - init_state
				end if

				! deadcrootc_xfer 
				ptr => pptr%pcs%deadcrootc_xfer(p)
				init_state = ptr*wt_old
				change_state = ptr*dwt
				new_state = init_state+change_state
				if (wt_new /= 0._r8) then
					ptr = new_state/wt_new
					conv_cflux(p) = conv_cflux(p) + change_state
				else
					ptr = 0._r8
					conv_cflux(p) = conv_cflux(p) - init_state
				end if

				! gresp_storage 
				ptr => pptr%pcs%gresp_storage(p)
				init_state = ptr*wt_old
				change_state = ptr*dwt
				new_state = init_state+change_state
				if (wt_new /= 0._r8) then
					ptr = new_state/wt_new
					conv_cflux(p) = conv_cflux(p) + change_state
				else
					ptr = 0._r8
					conv_cflux(p) = conv_cflux(p) - init_state
				end if

				! gresp_xfer 
				ptr => pptr%pcs%gresp_xfer(p)
				init_state = ptr*wt_old
				change_state = ptr*dwt
				new_state = init_state+change_state
				if (wt_new /= 0._r8) then
					ptr = new_state/wt_new
					conv_cflux(p) = conv_cflux(p) + change_state
				else
					ptr = 0._r8
					conv_cflux(p) = conv_cflux(p) - init_state
				end if

				! cpool 
				ptr => pptr%pcs%cpool(p)
				init_state = ptr*wt_old
				change_state = ptr*dwt
				new_state = init_state+change_state
				if (wt_new /= 0._r8) then
					ptr = new_state/wt_new
					conv_cflux(p) = conv_cflux(p) + change_state
				else
					ptr = 0._r8
					conv_cflux(p) = conv_cflux(p) - init_state
				end if

				! xsmrpool 
				ptr => pptr%pcs%xsmrpool(p)
				init_state = ptr*wt_old
				change_state = ptr*dwt
				new_state = init_state+change_state
				if (wt_new /= 0._r8) then
					ptr = new_state/wt_new
					conv_cflux(p) = conv_cflux(p) + change_state
				else
					ptr = 0._r8
					conv_cflux(p) = conv_cflux(p) - init_state
				end if

				! pft_ctrunc 
				ptr => pptr%pcs%pft_ctrunc(p)
				init_state = ptr*wt_old
				change_state = ptr*dwt
				new_state = init_state+change_state
				if (wt_new /= 0._r8) then
					ptr = new_state/wt_new
					conv_cflux(p) = conv_cflux(p) + change_state
				else
					ptr = 0._r8
					conv_cflux(p) = conv_cflux(p) - init_state
				end if

				!-----------------
				! C13 state update
				!-----------------

				! set pointers to the conversion and product pool fluxes for this pft
				! dwt_ptr0 is reserved for local assignment to dwt_xxx_to_litter fluxes
				dwt_ptr1 => conv_c13flux(p)
				dwt_ptr2 => prod10_c13flux(p)
				dwt_ptr3 => prod100_c13flux(p)

				! leafc 
				ptr => pptr%pc13s%leafc(p)
				init_state = ptr*wt_old
				change_state = ptr*dwt
				new_state = init_state+change_state
				if (wt_new /= 0._r8) then
					ptr = new_state/wt_new
					dwt_ptr1 = dwt_ptr1 + change_state
				else
					ptr = 0._r8
					dwt_ptr1 = dwt_ptr1 - init_state
				end if

				! leafc_storage 
				ptr => pptr%pc13s%leafc_storage(p)
				init_state = ptr*wt_old
				change_state = ptr*dwt
				new_state = init_state+change_state
				if (wt_new /= 0._r8) then
					ptr = new_state/wt_new
					dwt_ptr1 = dwt_ptr1 + change_state
				else
					ptr = 0._r8
					dwt_ptr1 = dwt_ptr1 - init_state
				end if

				! leafc_xfer 
				ptr => pptr%pc13s%leafc_xfer(p)
				init_state = ptr*wt_old
				change_state = ptr*dwt
				new_state = init_state+change_state
				if (wt_new /= 0._r8) then
					ptr = new_state/wt_new
					dwt_ptr1 = dwt_ptr1 + change_state
				else
					ptr = 0._r8
					dwt_ptr1 = dwt_ptr1 - init_state
				end if

				! frootc 
				ptr => pptr%pc13s%frootc(p)
				dwt_ptr0 => dwt_frootc13_to_litter(p)
				init_state = ptr*wt_old
				change_state = ptr*dwt
				new_state = init_state+change_state
				if (wt_new /= 0._r8) then
					ptr = new_state/wt_new
					dwt_ptr0 = dwt_ptr0 - change_state
				else
					ptr = 0._r8
					dwt_ptr0 = dwt_ptr0 + init_state
				end if

				! frootc_storage 
				ptr => pptr%pc13s%frootc_storage(p)
				init_state = ptr*wt_old
				change_state = ptr*dwt
				new_state = init_state+change_state
				if (wt_new /= 0._r8) then
					ptr = new_state/wt_new
					dwt_ptr1 = dwt_ptr1 + change_state
				else
					ptr = 0._r8
					dwt_ptr1 = dwt_ptr1 - init_state
				end if

				! frootc_xfer 
				ptr => pptr%pc13s%frootc_xfer(p)
				init_state = ptr*wt_old
				change_state = ptr*dwt
				new_state = init_state+change_state
				if (wt_new /= 0._r8) then
					ptr = new_state/wt_new
					dwt_ptr1 = dwt_ptr1 + change_state
				else
					ptr = 0._r8
					dwt_ptr1 = dwt_ptr1 - init_state
				end if

				! livestemc 
				ptr => pptr%pc13s%livestemc(p)
				init_state = ptr*wt_old
				change_state = ptr*dwt
				new_state = init_state+change_state
				if (wt_new /= 0._r8) then
					ptr = new_state/wt_new
					dwt_ptr1 = dwt_ptr1 + change_state
				else
					ptr = 0._r8
					dwt_ptr1 = dwt_ptr1 - init_state
				end if

				! livestemc_storage 
				ptr => pptr%pc13s%livestemc_storage(p)
				init_state = ptr*wt_old
				change_state = ptr*dwt
				new_state = init_state+change_state
				if (wt_new /= 0._r8) then
					ptr = new_state/wt_new
					dwt_ptr1 = dwt_ptr1 + change_state
				else
					ptr = 0._r8
					dwt_ptr1 = dwt_ptr1 - init_state
				end if

				! livestemc_xfer 
				ptr => pptr%pc13s%livestemc_xfer(p)
				init_state = ptr*wt_old
				change_state = ptr*dwt
				new_state = init_state+change_state
				if (wt_new /= 0._r8) then
					ptr = new_state/wt_new
					dwt_ptr1 = dwt_ptr1 + change_state
				else
					ptr = 0._r8
					dwt_ptr1 = dwt_ptr1 - init_state
				end if

				! deadstemc 
				ptr => pptr%pc13s%deadstemc(p)
				init_state = ptr*wt_old
				change_state = ptr*dwt
				new_state = init_state+change_state
				if (wt_new /= 0._r8) then
					ptr = new_state/wt_new
					dwt_ptr1 = dwt_ptr1 + change_state*pconv(pptr%itype(p))
					dwt_ptr2 = dwt_ptr2 + change_state*pprod10(pptr%itype(p))
					dwt_ptr3 = dwt_ptr3 + change_state*pprod100(pptr%itype(p))
				else
					ptr = 0._r8
					dwt_ptr1 = dwt_ptr1 - init_state*pconv(pptr%itype(p))
					dwt_ptr2 = dwt_ptr2 - init_state*pprod10(pptr%itype(p))
					dwt_ptr3 = dwt_ptr3 - init_state*pprod100(pptr%itype(p))
				end if

				! deadstemc_storage 
				ptr => pptr%pc13s%deadstemc_storage(p)
				init_state = ptr*wt_old
				change_state = ptr*dwt
				new_state = init_state+change_state
				if (wt_new /= 0._r8) then
					ptr = new_state/wt_new
					dwt_ptr1 = dwt_ptr1 + change_state
				else
					ptr = 0._r8
					dwt_ptr1 = dwt_ptr1 - init_state
				end if

				! deadstemc_xfer 
				ptr => pptr%pc13s%deadstemc_xfer(p)
				init_state = ptr*wt_old
				change_state = ptr*dwt
				new_state = init_state+change_state
				if (wt_new /= 0._r8) then
					ptr = new_state/wt_new
					dwt_ptr1 = dwt_ptr1 + change_state
				else
					ptr = 0._r8
					dwt_ptr1 = dwt_ptr1 - init_state
				end if

				! livecrootc 
				ptr => pptr%pc13s%livecrootc(p)
				dwt_ptr0 => dwt_livecrootc13_to_litter(p)
				init_state = ptr*wt_old
				change_state = ptr*dwt
				new_state = init_state+change_state
				if (wt_new /= 0._r8) then
					ptr = new_state/wt_new
					dwt_ptr0 = dwt_ptr0 - change_state
				else
					ptr = 0._r8
					dwt_ptr0 = dwt_ptr0 + init_state
				end if

				! livecrootc_storage 
				ptr => pptr%pc13s%livecrootc_storage(p)
				init_state = ptr*wt_old
				change_state = ptr*dwt
				new_state = init_state+change_state
				if (wt_new /= 0._r8) then
					ptr = new_state/wt_new
					dwt_ptr1 = dwt_ptr1 + change_state
				else
					ptr = 0._r8
					dwt_ptr1 = dwt_ptr1 - init_state
				end if

				! livecrootc_xfer 
				ptr => pptr%pc13s%livecrootc_xfer(p)
				init_state = ptr*wt_old
				change_state = ptr*dwt
				new_state = init_state+change_state
				if (wt_new /= 0._r8) then
					ptr = new_state/wt_new
					dwt_ptr1 = dwt_ptr1 + change_state
				else
					ptr = 0._r8
					dwt_ptr1 = dwt_ptr1 - init_state
				end if

				! deadcrootc 
				ptr => pptr%pc13s%deadcrootc(p)
				dwt_ptr0 => dwt_deadcrootc13_to_litter(p)
				init_state = ptr*wt_old
				change_state = ptr*dwt
				new_state = init_state+change_state
				if (wt_new /= 0._r8) then
					ptr = new_state/wt_new
					dwt_ptr0 = dwt_ptr0 - change_state
				else
					ptr = 0._r8
					dwt_ptr0 = dwt_ptr0 + init_state
				end if

				! deadcrootc_storage 
				ptr => pptr%pc13s%deadcrootc_storage(p)
				init_state = ptr*wt_old
				change_state = ptr*dwt
				new_state = init_state+change_state
				if (wt_new /= 0._r8) then
					ptr = new_state/wt_new
					dwt_ptr1 = dwt_ptr1 + change_state
				else
					ptr = 0._r8
					dwt_ptr1 = dwt_ptr1 - init_state
				end if

				! deadcrootc_xfer 
				ptr => pptr%pc13s%deadcrootc_xfer(p)
				init_state = ptr*wt_old
				change_state = ptr*dwt
				new_state = init_state+change_state
				if (wt_new /= 0._r8) then
					ptr = new_state/wt_new
					dwt_ptr1 = dwt_ptr1 + change_state
				else
					ptr = 0._r8
					dwt_ptr1 = dwt_ptr1 - init_state
				end if

				! gresp_storage 
				ptr => pptr%pc13s%gresp_storage(p)
				init_state = ptr*wt_old
				change_state = ptr*dwt
				new_state = init_state+change_state
				if (wt_new /= 0._r8) then
					ptr = new_state/wt_new
					dwt_ptr1 = dwt_ptr1 + change_state
				else
					ptr = 0._r8
					dwt_ptr1 = dwt_ptr1 - init_state
				end if

				! gresp_xfer 
				ptr => pptr%pc13s%gresp_xfer(p)
				init_state = ptr*wt_old
				change_state = ptr*dwt
				new_state = init_state+change_state
				if (wt_new /= 0._r8) then
					ptr = new_state/wt_new
					dwt_ptr1 = dwt_ptr1 + change_state
				else
					ptr = 0._r8
					dwt_ptr1 = dwt_ptr1 - init_state
				end if

				! cpool 
				ptr => pptr%pc13s%cpool(p)
				init_state = ptr*wt_old
				change_state = ptr*dwt
				new_state = init_state+change_state
				if (wt_new /= 0._r8) then
					ptr = new_state/wt_new
					dwt_ptr1 = dwt_ptr1 + change_state
				else
					ptr = 0._r8
					dwt_ptr1 = dwt_ptr1 - init_state
				end if

				! pft_ctrunc 
				ptr => pptr%pc13s%pft_ctrunc(p)
				init_state = ptr*wt_old
				change_state = ptr*dwt
				new_state = init_state+change_state
				if (wt_new /= 0._r8) then
					ptr = new_state/wt_new
					dwt_ptr1 = dwt_ptr1 + change_state
				else
					ptr = 0._r8
					dwt_ptr1 = dwt_ptr1 - init_state
				end if

				!---------------
				! N state update
				!---------------

				! set pointers to the conversion and product pool fluxes for this pft
				! dwt_ptr0 is reserved for local assignment to dwt_xxx_to_litter fluxes
				dwt_ptr1 => conv_nflux(p)
				dwt_ptr2 => prod10_nflux(p)
				dwt_ptr3 => prod100_nflux(p)

				! leafn 
				ptr => pptr%pns%leafn(p)
				init_state = ptr*wt_old
				change_state = ptr*dwt
				new_state = init_state+change_state
				if (wt_new /= 0._r8) then
					ptr = new_state/wt_new
					dwt_ptr1 = dwt_ptr1 + change_state
				else
					ptr = 0._r8
					dwt_ptr1 = dwt_ptr1 - init_state
				end if

				! leafn_storage  
				ptr => pptr%pns%leafn_storage(p)
				init_state = ptr*wt_old
				change_state = ptr*dwt
				new_state = init_state+change_state
				if (wt_new /= 0._r8) then
					ptr = new_state/wt_new
					dwt_ptr1 = dwt_ptr1 + change_state
				else
					ptr = 0._r8
					dwt_ptr1 = dwt_ptr1 - init_state
				end if

				! leafn_xfer  
				ptr => pptr%pns%leafn_xfer(p)
				init_state = ptr*wt_old
				change_state = ptr*dwt
				new_state = init_state+change_state
				if (wt_new /= 0._r8) then
					ptr = new_state/wt_new
					dwt_ptr1 = dwt_ptr1 + change_state
				else
					ptr = 0._r8
					dwt_ptr1 = dwt_ptr1 - init_state
				end if

				! frootn 
				ptr => pptr%pns%frootn(p)
				dwt_ptr0 => dwt_frootn_to_litter(p)
				init_state = ptr*wt_old
				change_state = ptr*dwt
				new_state = init_state+change_state
				if (wt_new /= 0._r8) then
					ptr = new_state/wt_new
					dwt_ptr0 = dwt_ptr0 - change_state
				else
					ptr = 0._r8
					dwt_ptr0 = dwt_ptr0 + init_state
				end if

				! frootn_storage 
				ptr => pptr%pns%frootn_storage(p)
				init_state = ptr*wt_old
				change_state = ptr*dwt
				new_state = init_state+change_state
				if (wt_new /= 0._r8) then
					ptr = new_state/wt_new
					dwt_ptr1 = dwt_ptr1 + change_state
				else
					ptr = 0._r8
					dwt_ptr1 = dwt_ptr1 - init_state
				end if

				! frootn_xfer  
				ptr => pptr%pns%frootn_xfer(p)
				init_state = ptr*wt_old
				change_state = ptr*dwt
				new_state = init_state+change_state
				if (wt_new /= 0._r8) then
					ptr = new_state/wt_new
					dwt_ptr1 = dwt_ptr1 + change_state
				else
					ptr = 0._r8
					dwt_ptr1 = dwt_ptr1 - init_state
				end if

				! livestemn  
				ptr => pptr%pns%livestemn(p)
				init_state = ptr*wt_old
				change_state = ptr*dwt
				new_state = init_state+change_state
				if (wt_new /= 0._r8) then
					ptr = new_state/wt_new
					dwt_ptr1 = dwt_ptr1 + change_state
				else
					ptr = 0._r8
					dwt_ptr1 = dwt_ptr1 - init_state
				end if

				! livestemn_storage 
				ptr => pptr%pns%livestemn_storage(p)
				init_state = ptr*wt_old
				change_state = ptr*dwt
				new_state = init_state+change_state
				if (wt_new /= 0._r8) then
					ptr = new_state/wt_new
					dwt_ptr1 = dwt_ptr1 + change_state
				else
					ptr = 0._r8
					dwt_ptr1 = dwt_ptr1 - init_state
				end if

				! livestemn_xfer 
				ptr => pptr%pns%livestemn_xfer(p)
				init_state = ptr*wt_old
				change_state = ptr*dwt
				new_state = init_state+change_state
				if (wt_new /= 0._r8) then
					ptr = new_state/wt_new
					dwt_ptr1 = dwt_ptr1 + change_state
				else
					ptr = 0._r8
					dwt_ptr1 = dwt_ptr1 - init_state
				end if

				! deadstemn 
				ptr => pptr%pns%deadstemn(p)
				init_state = ptr*wt_old
				change_state = ptr*dwt
				new_state = init_state+change_state
				if (wt_new /= 0._r8) then
					ptr = new_state/wt_new
					dwt_ptr1 = dwt_ptr1 + change_state*pconv(pptr%itype(p))
					dwt_ptr2 = dwt_ptr2 + change_state*pprod10(pptr%itype(p))
					dwt_ptr3 = dwt_ptr3 + change_state*pprod100(pptr%itype(p))
				else
					ptr = 0._r8
					dwt_ptr1 = dwt_ptr1 - init_state*pconv(pptr%itype(p))
					dwt_ptr2 = dwt_ptr2 - init_state*pprod10(pptr%itype(p))
					dwt_ptr3 = dwt_ptr3 - init_state*pprod100(pptr%itype(p))
				end if

				! deadstemn_storage 
				ptr => pptr%pns%deadstemn_storage(p)
				init_state = ptr*wt_old
				change_state = ptr*dwt
				new_state = init_state+change_state
				if (wt_new /= 0._r8) then
					ptr = new_state/wt_new
					dwt_ptr1 = dwt_ptr1 + change_state
				else
					ptr = 0._r8
					dwt_ptr1 = dwt_ptr1 - init_state
				end if

				! deadstemn_xfer 
				ptr => pptr%pns%deadstemn_xfer(p)
				init_state = ptr*wt_old
				change_state = ptr*dwt
				new_state = init_state+change_state
				if (wt_new /= 0._r8) then
					ptr = new_state/wt_new
					dwt_ptr1 = dwt_ptr1 + change_state
				else
					ptr = 0._r8
					dwt_ptr1 = dwt_ptr1 - init_state
				end if

				! livecrootn 
				ptr => pptr%pns%livecrootn(p)
				dwt_ptr0 => dwt_livecrootn_to_litter(p)
				init_state = ptr*wt_old
				change_state = ptr*dwt
				new_state = init_state+change_state
				if (wt_new /= 0._r8) then
					ptr = new_state/wt_new
					dwt_ptr0 = dwt_ptr0 - change_state
				else
					ptr = 0._r8
					dwt_ptr0 = dwt_ptr0 + init_state
				end if

				! livecrootn_storage  
				ptr => pptr%pns%livecrootn_storage(p)
				init_state = ptr*wt_old
				change_state = ptr*dwt
				new_state = init_state+change_state
				if (wt_new /= 0._r8) then
					ptr = new_state/wt_new
					dwt_ptr1 = dwt_ptr1 + change_state
				else
					ptr = 0._r8
					dwt_ptr1 = dwt_ptr1 - init_state
				end if

				! livecrootn_xfer  
				ptr => pptr%pns%livecrootn_xfer(p)
				init_state = ptr*wt_old
				change_state = ptr*dwt
				new_state = init_state+change_state
				if (wt_new /= 0._r8) then
					ptr = new_state/wt_new
					dwt_ptr1 = dwt_ptr1 + change_state
				else
					ptr = 0._r8
					dwt_ptr1 = dwt_ptr1 - init_state
				end if

				! deadcrootn 
				ptr => pptr%pns%deadcrootn(p)
				dwt_ptr0 => dwt_deadcrootn_to_litter(p)
				init_state = ptr*wt_old
				change_state = ptr*dwt
				new_state = init_state+change_state
				if (wt_new /= 0._r8) then
					ptr = new_state/wt_new
					dwt_ptr0 = dwt_ptr0 - change_state
				else
					ptr = 0._r8
					dwt_ptr0 = dwt_ptr0 + init_state
				end if

				! deadcrootn_storage  
				ptr => pptr%pns%deadcrootn_storage(p)
				init_state = ptr*wt_old
				change_state = ptr*dwt
				new_state = init_state+change_state
				if (wt_new /= 0._r8) then
					ptr = new_state/wt_new
					dwt_ptr1 = dwt_ptr1 + change_state
				else
					ptr = 0._r8
					dwt_ptr1 = dwt_ptr1 - init_state
				end if

				! deadcrootn_xfer  
				ptr => pptr%pns%deadcrootn_xfer(p)
				init_state = ptr*wt_old
				change_state = ptr*dwt
				new_state = init_state+change_state
				if (wt_new /= 0._r8) then
					ptr = new_state/wt_new
					dwt_ptr1 = dwt_ptr1 + change_state
				else
					ptr = 0._r8
					dwt_ptr1 = dwt_ptr1 - init_state
				end if

				! retransn  
				ptr => pptr%pns%retransn(p)
				init_state = ptr*wt_old
				change_state = ptr*dwt
				new_state = init_state+change_state
				if (wt_new /= 0._r8) then
					ptr = new_state/wt_new
					dwt_ptr1 = dwt_ptr1 + change_state
				else
					ptr = 0._r8
					dwt_ptr1 = dwt_ptr1 - init_state
				end if

				! npool  
				ptr => pptr%pns%npool(p)
				init_state = ptr*wt_old
				change_state = ptr*dwt
				new_state = init_state+change_state
				if (wt_new /= 0._r8) then
					ptr = new_state/wt_new
					dwt_ptr1 = dwt_ptr1 + change_state
				else
					ptr = 0._r8
					dwt_ptr1 = dwt_ptr1 - init_state
				end if
				
				! pft_ntrunc  
				ptr => pptr%pns%pft_ntrunc(p)
				init_state = ptr*wt_old
				change_state = ptr*dwt
				new_state = init_state+change_state
				if (wt_new /= 0._r8) then
					ptr = new_state/wt_new
					dwt_ptr1 = dwt_ptr1 + change_state
				else
					ptr = 0._r8
					dwt_ptr1 = dwt_ptr1 - init_state
				end if
				
			end if       ! weight decreasing
		end if           ! is soil
	end do               ! pft loop
    
	! calculate column-level seeding fluxes
!dir$ nointerchange
	do pi = 1,max_pft_per_col
!dir$ concurrent
!cdir nodep
		do c = begc, endc
			if ( pi <=  cptr%npfts(c) ) then
				p = cptr%pfti(c) + pi - 1
				
				! C fluxes
				cptr%ccf%dwt_seedc_to_leaf(c) = cptr%ccf%dwt_seedc_to_leaf(c) + dwt_leafc_seed(p)/dt
				cptr%ccf%dwt_seedc_to_deadstem(c) = cptr%ccf%dwt_seedc_to_deadstem(c) + dwt_deadstemc_seed(p)/dt
				
				! C13 fluxes
				cptr%cc13f%dwt_seedc_to_leaf(c) = cptr%cc13f%dwt_seedc_to_leaf(c) + dwt_leafc13_seed(p)/dt
				cptr%cc13f%dwt_seedc_to_deadstem(c) = cptr%cc13f%dwt_seedc_to_deadstem(c) + dwt_deadstemc13_seed(p)/dt
				
				! N fluxes
				cptr%cnf%dwt_seedn_to_leaf(c) = cptr%cnf%dwt_seedn_to_leaf(c) + dwt_leafn_seed(p)/dt
				cptr%cnf%dwt_seedn_to_deadstem(c) = cptr%cnf%dwt_seedn_to_deadstem(c) + dwt_deadstemn_seed(p)/dt
			end if
		end do
	end do


	! calculate pft-to-column for fluxes into litter and CWD pools
!dir$ nointerchange
	do pi = 1,max_pft_per_col
!dir$ concurrent
!cdir nodep
		do c = begc, endc
			if ( pi <=  cptr%npfts(c) ) then
				p = cptr%pfti(c) + pi - 1

				! fine root litter carbon fluxes
				cptr%ccf%dwt_frootc_to_litr1c(c) = cptr%ccf%dwt_frootc_to_litr1c(c) + (dwt_frootc_to_litter(p)*pftcon%fr_flab(pptr%itype(p)))/dt
				cptr%ccf%dwt_frootc_to_litr2c(c) = cptr%ccf%dwt_frootc_to_litr2c(c) + (dwt_frootc_to_litter(p)*pftcon%fr_fcel(pptr%itype(p)))/dt
				cptr%ccf%dwt_frootc_to_litr3c(c) = cptr%ccf%dwt_frootc_to_litr3c(c) + (dwt_frootc_to_litter(p)*pftcon%fr_flig(pptr%itype(p)))/dt

				! fine root litter C13 fluxes
				cptr%cc13f%dwt_frootc_to_litr1c(c) = cptr%cc13f%dwt_frootc_to_litr1c(c) + &
                                                            (dwt_frootc13_to_litter(p)*pftcon%fr_flab(pptr%itype(p)))/dt
				cptr%cc13f%dwt_frootc_to_litr2c(c) = cptr%cc13f%dwt_frootc_to_litr2c(c) + &
                                                            (dwt_frootc13_to_litter(p)*pftcon%fr_fcel(pptr%itype(p)))/dt
				cptr%cc13f%dwt_frootc_to_litr3c(c) = cptr%cc13f%dwt_frootc_to_litr3c(c) + &
                                                            (dwt_frootc13_to_litter(p)*pftcon%fr_flig(pptr%itype(p)))/dt

				! fine root litter nitrogen fluxes
				cptr%cnf%dwt_frootn_to_litr1n(c) = cptr%cnf%dwt_frootn_to_litr1n(c) + (dwt_frootn_to_litter(p)*pftcon%fr_flab(pptr%itype(p)))/dt
				cptr%cnf%dwt_frootn_to_litr2n(c) = cptr%cnf%dwt_frootn_to_litr2n(c) + (dwt_frootn_to_litter(p)*pftcon%fr_fcel(pptr%itype(p)))/dt
				cptr%cnf%dwt_frootn_to_litr3n(c) = cptr%cnf%dwt_frootn_to_litr3n(c) + (dwt_frootn_to_litter(p)*pftcon%fr_flig(pptr%itype(p)))/dt

				! livecroot fluxes to cwd
				cptr%ccf%dwt_livecrootc_to_cwdc(c) = cptr%ccf%dwt_livecrootc_to_cwdc(c) + (dwt_livecrootc_to_litter(p))/dt
				cptr%cc13f%dwt_livecrootc_to_cwdc(c) = cptr%cc13f%dwt_livecrootc_to_cwdc(c) + (dwt_livecrootc13_to_litter(p))/dt
				cptr%cnf%dwt_livecrootn_to_cwdn(c) = cptr%cnf%dwt_livecrootn_to_cwdn(c) + (dwt_livecrootn_to_litter(p))/dt

				! deadcroot fluxes to cwd
				cptr%ccf%dwt_deadcrootc_to_cwdc(c) = cptr%ccf%dwt_deadcrootc_to_cwdc(c) + (dwt_deadcrootc_to_litter(p))/dt
				cptr%cc13f%dwt_deadcrootc_to_cwdc(c) = cptr%cc13f%dwt_deadcrootc_to_cwdc(c) + (dwt_deadcrootc13_to_litter(p))/dt
				cptr%cnf%dwt_deadcrootn_to_cwdn(c) = cptr%cnf%dwt_deadcrootn_to_cwdn(c) + (dwt_deadcrootn_to_litter(p))/dt
			end if
		end do
	end do

	! calculate pft-to-column for fluxes into product pools and conversion flux
!dir$ nointerchange
	do pi = 1,max_pft_per_col
!dir$ concurrent
!cdir nodep
		do c = begc,endc
			if (pi <= cptr%npfts(c)) then
				p = cptr%pfti(c) + pi - 1

				! column-level fluxes are accumulated as positive fluxes.
				! column-level C flux updates
				cptr%ccf%dwt_conv_cflux(c) = cptr%ccf%dwt_conv_cflux(c) - conv_cflux(p)/dt
				cptr%ccf%dwt_prod10c_gain(c) = cptr%ccf%dwt_prod10c_gain(c) - prod10_cflux(p)/dt
				cptr%ccf%dwt_prod100c_gain(c) = cptr%ccf%dwt_prod100c_gain(c) - prod100_cflux(p)/dt

				! column-level C13 flux updates
				cptr%cc13f%dwt_conv_cflux(c) = cptr%cc13f%dwt_conv_cflux(c) - conv_c13flux(p)/dt
				cptr%cc13f%dwt_prod10c_gain(c) = cptr%cc13f%dwt_prod10c_gain(c) - prod10_c13flux(p)/dt
				cptr%cc13f%dwt_prod100c_gain(c) = cptr%cc13f%dwt_prod100c_gain(c) - prod100_c13flux(p)/dt

				! column-level N flux updates
				cptr%cnf%dwt_conv_nflux(c) = cptr%cnf%dwt_conv_nflux(c) - conv_nflux(p)/dt
				cptr%cnf%dwt_prod10n_gain(c) = cptr%cnf%dwt_prod10n_gain(c) - prod10_nflux(p)/dt
				cptr%cnf%dwt_prod100n_gain(c) = cptr%cnf%dwt_prod100n_gain(c) - prod100_nflux(p)/dt

			end if
		end do
	end do

	! calculate column-level losses from product pools
	! the following (1/s) rate constants result in ~90% loss of initial state over 10 and 100 years,
	! respectively, using a discrete-time fractional decay algorithm.
	kprod10 = 7.2e-9
	kprod100 = 7.2e-10
!dir$ concurrent
!cdir nodep
	do c = begc,endc
		! calculate fluxes (1/sec)
		cptr%ccf%dwt_prod10c_loss(c) = cptr%ccs%prod10c(c) * kprod10
		cptr%ccf%dwt_prod100c_loss(c) = cptr%ccs%prod100c(c) * kprod100
		cptr%cc13f%dwt_prod10c_loss(c) = cptr%cc13s%prod10c(c) * kprod10
		cptr%cc13f%dwt_prod100c_loss(c) = cptr%cc13s%prod100c(c) * kprod100
		cptr%cnf%dwt_prod10n_loss(c) = cptr%cns%prod10n(c) * kprod10
		cptr%cnf%dwt_prod100n_loss(c) = cptr%cns%prod100n(c) * kprod100
	end do

	! Deallocate pft-level flux arrays
        deallocate(dwt_leafc_seed)
        deallocate(dwt_leafn_seed)
        deallocate(dwt_leafc13_seed)
        deallocate(dwt_deadstemc_seed)
        deallocate(dwt_deadstemn_seed)
        deallocate(dwt_deadstemc13_seed)
	deallocate(dwt_frootc_to_litter)
	deallocate(dwt_livecrootc_to_litter)
	deallocate(dwt_deadcrootc_to_litter)
	deallocate(dwt_frootc13_to_litter)
	deallocate(dwt_livecrootc13_to_litter)
	deallocate(dwt_deadcrootc13_to_litter)
	deallocate(dwt_frootn_to_litter)
	deallocate(dwt_livecrootn_to_litter)
	deallocate(dwt_deadcrootn_to_litter)
	deallocate(conv_cflux)
	deallocate(prod10_cflux)
	deallocate(prod100_cflux)
	deallocate(conv_c13flux)
	deallocate(prod10_c13flux)
	deallocate(prod100_c13flux)
	deallocate(conv_nflux)
	deallocate(prod10_nflux)
	deallocate(prod100_nflux)
    
end subroutine pftdyn_cnbal
#endif

end module pftdynMod
