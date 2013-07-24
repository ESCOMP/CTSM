module pftdynMod

  !---------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Determine pft weights at current time using dynamic landuse datasets.
  ! ASSUMES that only have one dynamic landuse dataset.
  !
  ! !USES:
  use spmdMod
  use clmtype
  use decompMod   , only : bounds_type
  use clm_varsur  , only : pctspec
  use clm_varpar  , only : max_pft_per_col
  use clm_varctl  , only : iulog, use_c13, use_c14, use_cn
  use shr_sys_mod , only : shr_sys_flush
  use abortutils  , only : endrun
  use ncdio_pio
  !
  ! !PUBLIC TYPES:
  implicit none
  private
  save
  public :: pftdyn_init
  public :: pftdyn_interp
  public :: pftdyn_wbal_init
  public :: pftdyn_wbal
  public :: pftdyn_cnbal
  public :: pftwt_init
  public :: pftwt_interp
  public :: CNHarvest
  public :: CNHarvestPftToColumn
  !
  ! ! PRIVATE TYPES
  integer  , allocatable   :: yearspft (:) ! 
  real(r8) , allocatable   :: wtpft1 (:,:) ! weight of each pft relative to the natural veg landunit, time 1 
  real(r8) , allocatable   :: wtpft2 (:,:) ! weight of each pft relative to the natural veg landunit, time 2 
  real(r8) , allocatable   :: harvest (:)  ! 
  real(r8) , allocatable   :: wtcol_old (:)! 
  integer :: nt1
  integer :: nt2
  integer :: ntimes
  logical :: do_harvest
  type(file_desc_t)  :: ncid   ! netcdf id
  ! default multiplication factor for epsilon for error checks
  real(r8), private, parameter :: eps_fact = 2._r8
  !---------------------------------------------------------------------------

contains
  
  !-----------------------------------------------------------------------
  subroutine pftdyn_init(bounds)
    !
    ! !DESCRIPTION:
    ! Initialize dynamic landuse dataset (position it to the right time samples
    ! that bound the initial model date)
    !
    ! !USES:
    use clm_time_manager, only : get_curr_date
    use clm_varctl  , only : fpftdyn
    use clm_varpar  , only : numpft, maxpatch_pft, natpft_lb, natpft_ub, natpft_size
    use clm_varcon  , only : numurbl, istsoil, istcrop
    use clm_varsur  , only : wt_lunit
    use fileutils   , only : getfil
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds  ! bounds
    !
    ! !LOCAL VARIABLES:
    integer  :: i,j,m,n,g,nl                    ! indices
    real(r8) :: sumpct                          ! sum for error check
    integer  :: varid                           ! netcdf ids
    integer  :: year                            ! year (0, ...) for nstep+1
    integer  :: mon                             ! month (1, ..., 12) for nstep+1
    integer  :: day                             ! day of month (1, ..., 31) for nstep+1
    integer  :: sec                             ! seconds into current date for nstep+1
    integer  :: ier, ret                        ! error status
    logical  :: found                           ! true => input dataset bounding dates found
    logical  :: readvar	                        ! true => variable is on input dataset
    ! leave the following as pointers instead of changing to allocatable since the ncd_io
    ! interface expects a pointer type
    real(r8) , pointer :: pctnatveg (:)        ! percent of gcell is natural vegetation landunit 
    real(r8) , pointer :: pctcrop (:)          ! percent of gcell is crop landunit 
    real(r8) , pointer :: pctgla (:)           ! percent of gcell is glacier 
    real(r8) , pointer :: pctlak (:)           ! percent of gcell is lake 
    real(r8) , pointer :: pctwet (:)           ! percent of gcell is wetland 
    real(r8) , pointer :: pcturb (:,:)         ! percent of gcell is urbanized 
    real(r8) , pointer :: pcturb_tot (:)       ! percent of grid cell is urban (sum over density classes) 
    character(len=256) :: locfn                ! local file name
    character(len= 32) :: subname='pftdyn_init'! subroutine name
    !-----------------------------------------------------------------------

    ! Error check

    if ( maxpatch_pft /= numpft+1 )then
       call endrun( subname//' maxpatch_pft does NOT equal numpft+1 -- this is invalid for dynamic PFT case' )
    end if

    allocate(pctnatveg(bounds%begg:bounds%endg))
    allocate(pctcrop(bounds%begg:bounds%endg))
    allocate(pctgla(bounds%begg:bounds%endg))
    allocate(pctlak(bounds%begg:bounds%endg))
    allocate(pctwet(bounds%begg:bounds%endg))
    allocate(pcturb(bounds%begg:bounds%endg,numurbl))
    allocate(pcturb_tot(bounds%begg:bounds%endg))

    ! pctspec must be saved between time samples
    ! position to first time sample - assume that first time sample must match starting date
    ! check consistency -  special landunits, grid, frac and mask
    ! only do this once

    ! read data PCT_NAT_PFT corresponding to correct year

    allocate(wtpft1(bounds%begg:bounds%endg,natpft_lb:natpft_ub), &
             wtpft2(bounds%begg:bounds%endg,natpft_lb:natpft_ub), &
             stat=ier)
    if (ier /= 0) then
       call endrun( subname//' allocation error for wtpft1, wtpft2' )
    end if
    
    allocate(harvest(bounds%begg:bounds%endg),stat=ier)
    if (ier /= 0) then
       call endrun( subname//' allocation error for harvest')
    end if

    allocate(wtcol_old(bounds%begp:bounds%endp),stat=ier)
    if (ier /= 0) then
       call endrun( subname//' allocation error for wtcol_old' )
    end if

    if (masterproc) then
       write(iulog,*) 'Attempting to read pft dynamic landuse data .....'
    end if

    ! Obtain file
    call getfil (fpftdyn, locfn, 0)
    call ncd_pio_openfile (ncid, locfn, 0)

    ! Obtain pft years from dynamic landuse file
    
    call ncd_inqdid(ncid, 'time', varid)
    call ncd_inqdlen(ncid, varid, ntimes)

    ! Consistency check
    
    call check_dim(ncid, 'natpft', natpft_size)

    allocate (yearspft(ntimes), stat=ier)
    if (ier /= 0) then
       write(iulog,*) subname//' allocation error for yearspft'; call endrun()
    end if

    call ncd_io(ncid=ncid, varname='YEAR', flag='read', data=yearspft)

    call ncd_io(ncid=ncid, varname='PCT_NATVEG', flag='read', data=pctnatveg, &
         dim1name=grlnd, readvar=readvar)
    if (.not. readvar) call endrun( trim(subname)//' ERROR: PCT_NATVEG NOT on pftdyn file' )

    call ncd_io(ncid=ncid, varname='PCT_CROP', flag='read', data=pctcrop, &
         dim1name=grlnd, readvar=readvar)
    if (.not. readvar) call endrun( trim(subname)//' ERROR: PCT_CROP NOT on pftdyn file' )

    call ncd_io(ncid=ncid, varname='PCT_WETLAND', flag='read', data=pctwet, &
         dim1name=grlnd, readvar=readvar)
    if (.not. readvar) call endrun( trim(subname)//' ERROR: PCT_WETLAND NOT on pftdyn file' )

    call ncd_io(ncid=ncid, varname= 'PCT_LAKE', flag='read', data=pctlak, &
         dim1name=grlnd, readvar=readvar)
    if (.not. readvar) call endrun( trim(subname)//' ERROR: PCT_LAKE NOT on pftdyn file' )

    call ncd_io(ncid=ncid, varname= 'PCT_GLACIER', flag='read', data=pctgla, &
         dim1name=grlnd, readvar=readvar)
    if (.not. readvar) call endrun( trim(subname)//' ERROR: PCT_GLACIER NOT on pftdyn file' )

    call ncd_io(ncid=ncid, varname= 'PCT_URBAN'  , flag='read', data=pcturb, &
         dim1name=grlnd, readvar=readvar)
    if (.not. readvar) call endrun( trim(subname)//' ERROR: PCT_URBAN NOT on pftdyn file' )
    pcturb_tot(:) = 0._r8
    do n = 1, numurbl
       do nl = bounds%begg,bounds%endg
          pcturb_tot(nl) = pcturb_tot(nl) + pcturb(nl,n)
       enddo
    enddo

    ! Consistency check
    do g = bounds%begg,bounds%endg
       ! WJS (5-9-13): I am adding these pctnatveg and pctcrop consistency checks for now,
       ! although both these and the following pctspec consistency check (which was
       ! already here) may not be necessary now that pctpft is specified as % of landunit
       ! rather than % of grid cell. Furthermore, once we have dynamic landunits, both of
       ! these consistency checks may become problematic, and could be removed. If these
       ! consistency checks are removed, I think we could do quite a bit of cleanup:
       ! (1) I think that the time-invariant PCT fields no longer need to be on the pftdyn dataset
       ! (2) Similarly, much of this routine could be removed
       ! (3) I think that pctspec no longer needs to be saved (I think it's just saved
       !     for the sake of this consistency check)
       ! (4) wt_lunit no longer needs to be saved past the end of the initialize1 routine
       !     in clm_initializeMod (so we can move its deallocation from initialize2 to
       !     initialize1)
       if (abs(pctnatveg(g) - wt_lunit(g,istsoil)*100._r8) > 1.e-13_r8) then
          write(iulog,*) subname//'mismatch between input PCT_NATVEG = ', pctnatveg(g), &
               ' and that obtained from surface dataset ', wt_lunit(g,istsoil)*100._r8, &
               ' at g= ',g
          call endrun()
       end if
       if (abs(pctcrop(g) - wt_lunit(g,istcrop)*100._r8) > 1.e-13_r8) then
          write(iulog,*) subname//'mismatch between input PCT_CROP = ', pctcrop(g), &
               ' and that obtained from surface dataset ', wt_lunit(g,istcrop)*100._r8, &
               ' at g= ',g
          call endrun()
       end if

    !   This was causing a fail, even though values are the same to within 1e-15
    !   if (pctlak(g)+pctwet(g)+pcturb(g)+pctgla(g) /= pctspec(g)) then 
       if (abs((pctlak(g)+pctwet(g)+pcturb_tot(g)+pctgla(g))-pctspec(g)) > 1e-13_r8) then 
          write(iulog,*) subname//'mismatch between input pctspec = ',&
                     pctlak(g)+pctwet(g)+pcturb_tot(g)+pctgla(g),&
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
    ! PET: harvest - since harvest is specified as a rate for each year, this
    ! approach will not work. Instead, need to seta flag that indicates harvest is
    ! zero for the period before the beginning and after the end of the dynpft timeseries.

    call get_curr_date(year, mon, day, sec)

    if (year < yearspft(1)) then
       nt1 = 1
       nt2 = 1
       do_harvest = .false.
    else if (year >= yearspft(ntimes)) then
       nt1 = ntimes
       nt2 = ntimes
       do_harvest = .false.
    else
       found = .false.
       do n = 1,ntimes-1 
          if (year == yearspft(n)) then
             nt1 = n
             nt2 = nt1 + 1
             found = .true.
             do_harvest = .true.
          end if   
       end do
       if (.not. found) then
          write(iulog,*) subname//' error: model year not found in pftdyn timeseries'
          write(iulog,*)'model year = ',year
          call endrun()
       end if
    end if

    ! Get pctpft time samples bracketing the current time

    call pftdyn_getdata(bounds, nt1, wtpft1, natpft_lb,natpft_ub)
    call pftdyn_getdata(bounds, nt2, wtpft2, natpft_lb,natpft_ub)
    
    if (use_cn) then
       ! Get harvest rate at the nt1 time
       call pftdyn_getharvest(bounds, nt1)
    end if

    deallocate(pctnatveg, pctcrop, pctgla,pctlak,pctwet,pcturb,pcturb_tot)

  end subroutine pftdyn_init

  !-----------------------------------------------------------------------
  subroutine pftdyn_interp(bounds)
    !
    ! !DESCRIPTION:
    ! Time interpolate dynamic landuse data to get pft weights for model time
    ! Note that harvest data are stored as rates (not weights) and so time interpolation is 
    ! not necessary - the harvest rate is held constant through the year.  This is consistent with
    ! the treatment of changing PFT weights, where interpolation of the annual endpoint weights leads to 
    ! a constant rate of change in PFT weight through the year, with abrupt changes in the rate at
    ! annual boundaries. This routine is still used to get the next harvest time slice, when needed.
    ! This routine is also used to turn off the harvest switch when the model year runs past the end of
    ! the dynpft time series.
    !
    ! !USES:
    use clm_time_manager, only : get_curr_date, get_curr_calday, &
                                 get_days_per_year
    use clm_varcon      , only : istsoil
    use clm_varpar      , only : natpft_lb, natpft_ub
    implicit none
    !
    !
    ! !LOCAL VARIABLES:
    type(bounds_type), intent(in) :: bounds  ! bounds
    integer  :: i,j,m,p,l,g,c    ! indices
    integer  :: year             ! year (0, ...) for nstep+1
    integer  :: mon              ! month (1, ..., 12) for nstep+1
    integer  :: day              ! day of month (1, ..., 31) for nstep+1
    integer  :: sec              ! seconds into current date for nstep+1
    real(r8) :: cday             ! current calendar day (1.0 = 0Z on Jan 1)
    real(r8) :: days_per_year    ! days per year
    integer  :: ier              ! error status
    real(r8) :: wt1              ! time interpolation weights
    character(len=32) :: subname='pftdyn_interp' ! subroutine name
    !-----------------------------------------------------------------------

    ! Interpolate pctpft to current time step - output in pctpft_intp
    ! Map interpolated pctpft to subgrid weights
    ! assumes that maxpatch_pft = numpft + 1, that each landunit has only 1 column, 
    ! SCAM and CNDV have not been defined, and create_croplandunit = .false.

    ! If necessary, obtain new time sample

    ! Get current date

    call get_curr_date(year, mon, day, sec)

    ! Obtain new time sample if necessary.
    ! The first condition is the regular crossing of a year boundary
    ! when within the dynpft timeseries range. The second condition is
    ! the case of the first entry into the dynpft timeseries range from
    ! an earlier period of constant weights.

    if (year > yearspft(nt1) .or. (nt1 == 1 .and. nt2 == 1 .and. year == yearspft(1))) then

       if (year >= yearspft(ntimes)) then
          nt1 = ntimes
          nt2 = ntimes
       else
          nt1        = nt2
          nt2        = nt2 + 1
          do_harvest = .true.
       end if
       
       if (year > yearspft(ntimes)) then
          do_harvest = .false.
       endif
       
       if (nt2 > ntimes .and. masterproc) then
          write(iulog,*)subname,' error - current year is past input data boundary'
       end if
       
       do m = natpft_lb,natpft_ub
          do g = bounds%begg,bounds%endg
             wtpft1(g,m) = wtpft2(g,m)
          end do
       end do

       call pftdyn_getdata(bounds, nt2, wtpft2, natpft_lb,natpft_ub)

       if (use_cn) then
          call pftdyn_getharvest(bounds, nt1)
       end if

    end if  ! end of need new data if-block 

    ! Interpolate pft weight to current time

    cday          = get_curr_calday() 
    days_per_year = get_days_per_year()

    wt1 = ((days_per_year + 1._r8) - cday)/days_per_year

    do p = bounds%begp,bounds%endp
       c = pft%column(p)
       g = pft%gridcell(p)
       l = pft%landunit(p)

       ! Note that we only deal with the istsoil landunit here, NOT the istcrop landunit
       ! (if there is one)
       ! (However, currently [as of 5-9-13] the code won't let you run with transient
       ! PFTs combined with create_crop_landunit anyway, so it's a moot point.)
       if (lun%itype(l) == istsoil) then
          m = pft%itype(p)
          wtcol_old(p)      = pft%wtcol(p)

          ! Note that the following assignments assume that all PFTs share a single column

!         --- recoded for roundoff performance, tcraig 3/07 from k.lindsay
!         pft%wtcol(p)     = wtpft1(g,m)*wt1 + wtpft2(g,m)*wt2
          pft%wtcol(p)     = wtpft2(g,m) + wt1*(wtpft1(g,m)-wtpft2(g,m))
          pft%wtlunit(p)   = pft%wtcol(p)
          pft%wtgcell(p)   = pft%wtlunit(p) * lun%wtgcell(l)
       end if

    end do

  end subroutine pftdyn_interp

  !-----------------------------------------------------------------------
  subroutine pftdyn_getdata(bounds, ntime, wtpft, pft0, maxpft)
    !
    ! !DESCRIPTION:
    ! Obtain dynamic landuse data (wtpft) and make sure that
    ! percentage of PFTs sum to 100% cover for vegetated landunit
    !
    ! !USES:
    use surfrdMod   , only : surfrd_check_sums_equal_1
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds  ! bounds
    integer , intent(in)  :: ntime
    integer , intent(in)  :: pft0, maxpft
    real(r8), intent(out) :: wtpft(bounds%begg:bounds%endg,pft0:maxpft)  ! pft weights (sum to 1.0)
    !
    !
    ! !LOCAL VARIABLES:
    integer  :: i,j,m,n
    integer  :: err, ierr, ret
    real(r8) :: sumpct,sumerr                     ! temporary
    real(r8) , pointer :: arrayl (:,:)            ! temporary array 
    logical  :: readvar
    character(len=32) :: subname='pftdyn_getdata' ! subroutine name
    !-----------------------------------------------------------------------
    
    allocate(arrayl(bounds%begg:bounds%endg,pft0:maxpft))	
    call ncd_io(ncid=ncid, varname= 'PCT_NAT_PFT', flag='read', data=arrayl, &
         dim1name=grlnd, nt=ntime, readvar=readvar)
    wtpft(bounds%begg:bounds%endg,pft0:maxpft) = arrayl(bounds%begg:bounds%endg,pft0:maxpft) / 100._r8
    deallocate(arrayl)		
    if (.not. readvar) call endrun( trim(subname)//' ERROR: PCT_NAT_PFT NOT on pftdyn file' )
    
    call surfrd_check_sums_equal_1(wtpft, bounds%begg, 'PCT_NAT_PFT', subname)
    
  end subroutine pftdyn_getdata

  !-----------------------------------------------------------------------
  subroutine pftdyn_getharvest(bounds, ntime)
    !
    ! !DESCRIPTION:
    ! Obtain harvest data 
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds  ! bounds
    integer , intent(in)  :: ntime
    !
    ! !LOCAL VARIABLES:
    real(r8) , pointer :: arrayl (:)                    !  [real(r8) (:)]  temporary array 
    logical :: readvar 
    character(len=32) :: subname='pftdyn_getharvest' ! subroutine name
    !-----------------------------------------------------------------------
    
    allocate(arrayl(bounds%begg:bounds%endg))
    
    call ncd_io(ncid=ncid, varname= 'HARVEST_VH1', flag='read', data=arrayl, dim1name=grlnd, &
         nt=ntime, readvar=readvar)
    if (.not. readvar) call endrun( trim(subname)//' ERROR: HARVEST_VH1 not on pftdyn file' )
    harvest(bounds%begg:bounds%endg) = arrayl(bounds%begg:bounds%endg)
    
    call ncd_io(ncid=ncid, varname= 'HARVEST_VH2', flag='read', data=arrayl, dim1name=grlnd, &
         nt=ntime, readvar=readvar)
    if (.not. readvar) call endrun( trim(subname)//' ERROR: HARVEST_VH2 not on pftdyn file' )
    harvest(bounds%begg:bounds%endg) = harvest(bounds%begg:bounds%endg) + arrayl(bounds%begg:bounds%endg)
    
    call ncd_io(ncid=ncid, varname= 'HARVEST_SH1', flag='read', data=arrayl, dim1name=grlnd, &
         nt=ntime, readvar=readvar)
    if (.not. readvar) call endrun( trim(subname)//' ERROR: HARVEST_SH1 not on pftdyn file' )
    harvest(bounds%begg:bounds%endg) = harvest(bounds%begg:bounds%endg) + arrayl(bounds%begg:bounds%endg)
    
    call ncd_io(ncid=ncid, varname= 'HARVEST_SH2', flag='read', data=arrayl, dim1name=grlnd, &
         nt=ntime, readvar=readvar)
    if (.not. readvar) call endrun( trim(subname)//' ERROR: HARVEST_SH2 not on pftdyn file' )
    harvest(bounds%begg:bounds%endg) = harvest(bounds%begg:bounds%endg) + arrayl(bounds%begg:bounds%endg)
    
    call ncd_io(ncid=ncid, varname= 'HARVEST_SH3', flag='read', data=arrayl, dim1name=grlnd, &
         nt=ntime, readvar=readvar)
    if (.not. readvar) call endrun( trim(subname)//' ERROR: HARVEST_SH3 not on pftdyn file' )
    harvest(bounds%begg:bounds%endg) = harvest(bounds%begg:bounds%endg) + arrayl(bounds%begg:bounds%endg)

    deallocate(arrayl)

  end subroutine pftdyn_getharvest

  !-----------------------------------------------------------------------
  subroutine pftdyn_wbal_init(bounds)
    !
    ! !DESCRIPTION:
    ! initialize the column-level mass-balance correction term.
    ! Called in every timestep.
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds  ! bounds
    !
    ! !LOCAL VARIABLES:
    integer  :: c             ! indices
    !-----------------------------------------------------------------------

    ! set column-level canopy water mass balance correction flux
    ! term to 0 at the beginning of every timestep
    
    do c = bounds%begc,bounds%endc
       cwf%h2ocan_loss(c) = 0._r8
    end do
    
  end subroutine pftdyn_wbal_init

  !-----------------------------------------------------------------------
  subroutine pftdyn_wbal(bounds)
    !
    ! !DESCRIPTION:
    ! modify pft-level state and flux variables to maintain water balance with
    ! dynamic pft-weights.
    ! Canopy water balance does not need to consider harvest fluxes, since pft weights are
    ! not affected by harvest.
    !
    ! !USES:
    use clm_varcon  , only : istsoil
    use clm_varcon  , only : istcrop
    use clm_time_manager, only : get_step_size
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds  ! bounds
    !
    ! !LOCAL VARIABLES:
    integer  :: pi,p,c,l,g    ! indices
    integer  :: ier           ! error code
    real(r8) :: dtime         ! land model time step (sec)
    real(r8) :: dwt           ! change in pft weight (relative to column)
    real(r8) :: init_h2ocan   ! initial canopy water mass
    real(r8) :: new_h2ocan    ! canopy water mass after weight shift
    real(r8), allocatable :: loss_h2ocan(:) ! canopy water mass loss due to weight shift
    character(len=32) :: subname='pftdyn_wbal' ! subroutine name
    !-----------------------------------------------------------------------

    ! Allocate loss_h2ocan
    allocate(loss_h2ocan(bounds%begp:bounds%endp), stat=ier)
    if (ier /= 0) then
          write(iulog,*)subname,' allocation error for loss_h2ocan'; call endrun()
    end if

    ! Get time step

    dtime = get_step_size()

    ! set column-level canopy water mass balance correction flux
    ! term to 0 at the beginning of every weight-shifting timestep

    do c = bounds%begc,bounds%endc
       cwf%h2ocan_loss(c) = 0._r8 ! is this OR pftdyn_wbal_init redundant?
    end do

    do p = bounds%begp,bounds%endp
       l = pft%landunit(p)
       loss_h2ocan(p) = 0._r8

       if (lun%itype(l) == istsoil .or. lun%itype(l) == istcrop) then

          ! calculate the change in weight for the timestep
          dwt = pft%wtcol(p)-wtcol_old(p)
  
          if (dwt > 0._r8) then
          
             ! if the pft gained weight, then the 
             ! initial canopy water state is redistributed over the
             ! new (larger) area, conserving mass.

             pws%h2ocan(p) = pws%h2ocan(p) * (wtcol_old(p)/pft%wtcol(p))
          
          else if (dwt < 0._r8) then
          
             ! if the pft lost weight on the timestep, then the canopy water
             ! mass associated with the lost weight is directed to a 
             ! column-level flux term that gets added to the precip flux
             ! for every pft calculation in Hydrology1()
             
             init_h2ocan = pws%h2ocan(p) * wtcol_old(p)
             loss_h2ocan(p) = pws%h2ocan(p) * (-dwt)
             new_h2ocan = init_h2ocan - loss_h2ocan(p)
             if (abs(new_h2ocan) < 1e-8_r8) then
                new_h2ocan = 0._r8
                loss_h2ocan(p) = init_h2ocan
             end if
             if (pft%wtcol(p) /= 0._r8) then  
                pws%h2ocan(p) = new_h2ocan/pft%wtcol(p)
             else
                pws%h2ocan(p) = 0._r8
                loss_h2ocan(p) = init_h2ocan
             end if 
       

          end if

       end if
    end do

    do pi = 1,max_pft_per_col
       do c = bounds%begc,bounds%endc
          if (pi <= col%npfts(c)) then
             p = col%pfti(c) + pi - 1
             cwf%h2ocan_loss(c) = cwf%h2ocan_loss(c) + loss_h2ocan(p)/dtime
          end if
       end do
    end do

    ! Deallocate loss_h2ocan
    deallocate(loss_h2ocan)
    
  end subroutine pftdyn_wbal
  
  !-----------------------------------------------------------------------
  subroutine pftdyn_cnbal(bounds)
    !
    ! !DESCRIPTION:
    ! modify pft-level state and flux variables to maintain carbon and nitrogen balance with
    ! dynamic pft-weights.
    !
    ! !USES:
    use shr_const_mod,only : SHR_CONST_PDB
    use clm_varcon  , only : istsoil
    use clm_varpar  , only : numveg, nlevdecomp
    use clm_varcon  , only : istcrop
    use pftvarcon   , only : pconv, pprod10, pprod100
    use clm_varcon  , only : c13ratio, c14ratio
    use clm_time_manager, only : get_step_size
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds  ! bounds
    !
    ! !LOCAL VARIABLES:
    integer  :: pi,p,c,l,g,j    ! indices
    integer  :: ier           ! error code
    real(r8) :: dwt           ! change in pft weight (relative to column)
    real(r8) :: dt            ! land model time step (sec)
    real(r8) :: init_h2ocan   ! initial canopy water mass
    real(r8) :: new_h2ocan    ! canopy water mass after weight shift
    real(r8), allocatable :: dwt_leafc_seed(:)       ! pft-level mass gain due to seeding of new area
    real(r8), allocatable :: dwt_leafn_seed(:)       ! pft-level mass gain due to seeding of new area
    real(r8), allocatable :: dwt_deadstemc_seed(:)       ! pft-level mass gain due to seeding of new area
    real(r8), allocatable :: dwt_deadstemn_seed(:)       ! pft-level mass gain due to seeding of new area
    real(r8), allocatable :: dwt_frootc_to_litter(:)       ! pft-level mass loss due to weight shift
    real(r8), allocatable :: dwt_livecrootc_to_litter(:)   ! pft-level mass loss due to weight shift
    real(r8), allocatable :: dwt_deadcrootc_to_litter(:)   ! pft-level mass loss due to weight shift
    real(r8), allocatable, target :: dwt_frootn_to_litter(:)       ! pft-level mass loss due to weight shift
    real(r8), allocatable, target :: dwt_livecrootn_to_litter(:)   ! pft-level mass loss due to weight shift
    real(r8), allocatable, target :: dwt_deadcrootn_to_litter(:)   ! pft-level mass loss due to weight shift
    real(r8), allocatable :: conv_cflux(:)         ! pft-level mass loss due to weight shift
    real(r8), allocatable :: prod10_cflux(:)       ! pft-level mass loss due to weight shift
    real(r8), allocatable :: prod100_cflux(:)      ! pft-level mass loss due to weight shift
    real(r8), allocatable, target :: conv_nflux(:)         ! pft-level mass loss due to weight shift
    real(r8), allocatable, target :: prod10_nflux(:)       ! pft-level mass loss due to weight shift
    real(r8), allocatable, target :: prod100_nflux(:)      ! pft-level mass loss due to weight shift
    real(r8) :: t1,t2,wt_new,wt_old
    real(r8) :: init_state, change_state, new_state
    real(r8) :: tot_leaf, pleaf, pstor, pxfer
    real(r8) :: leafc_seed, leafn_seed
    real(r8) :: deadstemc_seed, deadstemn_seed
    real(r8), pointer :: dwt_ptr0, dwt_ptr1, dwt_ptr2, dwt_ptr3, ptr
    character(len=32) :: subname='pftdyn_cbal' ! subroutine name
    !! C13
    real(r8), allocatable :: dwt_leafc13_seed(:)     ! pft-level mass gain due to seeding of new area
    real(r8), allocatable :: dwt_deadstemc13_seed(:)     ! pft-level mass gain due to seeding of new area
    real(r8), allocatable, target :: dwt_frootc13_to_litter(:)     ! pft-level mass loss due to weight shift
    real(r8), allocatable, target :: dwt_livecrootc13_to_litter(:) ! pft-level mass loss due to weight shift
    real(r8), allocatable, target :: dwt_deadcrootc13_to_litter(:) ! pft-level mass loss due to weight shift
    real(r8), allocatable, target :: conv_c13flux(:)       ! pft-level mass loss due to weight shift
    real(r8), allocatable, target :: prod10_c13flux(:)     ! pft-level mass loss due to weight shift
    real(r8), allocatable, target :: prod100_c13flux(:)    ! pft-level mass loss due to weight shift
    real(r8) :: c3_del13c     ! typical del13C for C3 photosynthesis (permil, relative to PDB)
    real(r8) :: c4_del13c     ! typical del13C for C4 photosynthesis (permil, relative to PDB)
    real(r8) :: c3_r1_c13         ! isotope ratio (13c/12c) for C3 photosynthesis
    real(r8) :: c4_r1_c13         ! isotope ratio (13c/12c) for C4 photosynthesis
    real(r8) :: c3_r2_c13         ! isotope ratio (13c/[12c+13c]) for C3 photosynthesis
    real(r8) :: c4_r2_c13         ! isotope ratio (13c/[12c+13c]) for C4 photosynthesis
    real(r8) :: leafc13_seed, deadstemc13_seed
    !! C14
    real(r8), allocatable :: dwt_leafc14_seed(:)     ! pft-level mass gain due to seeding of new area
    real(r8), allocatable :: dwt_deadstemc14_seed(:)     ! pft-level mass gain due to seeding of new area
    real(r8), allocatable, target :: dwt_frootc14_to_litter(:)     ! pft-level mass loss due to weight shift
    real(r8), allocatable, target :: dwt_livecrootc14_to_litter(:) ! pft-level mass loss due to weight shift
    real(r8), allocatable, target :: dwt_deadcrootc14_to_litter(:) ! pft-level mass loss due to weight shift
    real(r8), allocatable, target :: conv_c14flux(:)       ! pft-level mass loss due to weight shift
    real(r8), allocatable, target :: prod10_c14flux(:)     ! pft-level mass loss due to weight shift
    real(r8), allocatable, target :: prod100_c14flux(:)    ! pft-level mass loss due to weight shift
    real(r8) :: c3_del14c     ! typical del14C for C3 photosynthesis (permil, relative to PDB)
    real(r8) :: c4_del14c     ! typical del14C for C4 photosynthesis (permil, relative to PDB)
    real(r8) :: c3_r1_c14         ! isotope ratio (14c/12c) for C3 photosynthesis
    real(r8) :: c4_r1_c14         ! isotope ratio (14c/12c) for C4 photosynthesis
    real(r8) :: c3_r2_c14         ! isotope ratio (14c/[12c+14c]) for C3 photosynthesis
    real(r8) :: c4_r2_c14         ! isotope ratio (14c/[12c+14c]) for C4 photosynthesis
    real(r8) :: leafc14_seed, deadstemc14_seed
    !-----------------------------------------------------------------------
    
   associate(& 
   lfpftd  =>  pps%lfpftd  & ! Output:  [real(r8) (:)] F. Li and S. Levis                                        
   )

    ! Allocate pft-level mass loss arrays
    allocate(dwt_leafc_seed(bounds%begp:bounds%endp), stat=ier)
    if (ier /= 0) then
          write(iulog,*)subname,' allocation error for dwt_leafc_seed'; call endrun()
    end if
    allocate(dwt_leafn_seed(bounds%begp:bounds%endp), stat=ier)
    if (ier /= 0) then
          write(iulog,*)subname,' allocation error for dwt_leafn_seed'; call endrun()
    end if
    allocate(dwt_deadstemc_seed(bounds%begp:bounds%endp), stat=ier)
    if (ier /= 0) then
          write(iulog,*)subname,' allocation error for dwt_deadstemc_seed'; call endrun()
    end if
    allocate(dwt_deadstemn_seed(bounds%begp:bounds%endp), stat=ier)
    if (ier /= 0) then
          write(iulog,*)subname,' allocation error for dwt_deadstemn_seed'; call endrun()
    end if
    allocate(dwt_frootc_to_litter(bounds%begp:bounds%endp), stat=ier)
    if (ier /= 0) then
          write(iulog,*)subname,' allocation error for dwt_frootc_to_litter'; call endrun()
    end if
    allocate(dwt_livecrootc_to_litter(bounds%begp:bounds%endp), stat=ier)
    if (ier /= 0) then
          write(iulog,*)subname,' allocation error for dwt_livecrootc_to_litter'; call endrun()
    end if
    allocate(dwt_deadcrootc_to_litter(bounds%begp:bounds%endp), stat=ier)
    if (ier /= 0) then
          write(iulog,*)subname,' allocation error for dwt_deadcrootc_to_litter'; call endrun()
    end if
    allocate(dwt_frootn_to_litter(bounds%begp:bounds%endp), stat=ier)
    if (ier /= 0) then
          write(iulog,*)subname,' allocation error for dwt_frootn_to_litter'; call endrun()
    end if
    allocate(dwt_livecrootn_to_litter(bounds%begp:bounds%endp), stat=ier)
    if (ier /= 0) then
          write(iulog,*)subname,' allocation error for dwt_livecrootn_to_litter'; call endrun()
    end if
    allocate(dwt_deadcrootn_to_litter(bounds%begp:bounds%endp), stat=ier)
    if (ier /= 0) then
          write(iulog,*)subname,' allocation error for dwt_deadcrootn_to_litter'; call endrun()
    end if
    allocate(conv_cflux(bounds%begp:bounds%endp), stat=ier)
    if (ier /= 0) then
          write(iulog,*)subname,' allocation error for conv_cflux'; call endrun()
    end if
    allocate(prod10_cflux(bounds%begp:bounds%endp), stat=ier)
    if (ier /= 0) then
          write(iulog,*)subname,' allocation error for prod10_cflux'; call endrun()
    end if
    allocate(prod100_cflux(bounds%begp:bounds%endp), stat=ier)
    if (ier /= 0) then
          write(iulog,*)subname,' allocation error for prod100_cflux'; call endrun()
    end if
    allocate(conv_nflux(bounds%begp:bounds%endp), stat=ier)
    if (ier /= 0) then
          write(iulog,*)subname,' allocation error for conv_nflux'; call endrun()
    end if
    allocate(prod10_nflux(bounds%begp:bounds%endp), stat=ier)
    if (ier /= 0) then
          write(iulog,*)subname,' allocation error for prod10_nflux'; call endrun()
    end if
    allocate(prod100_nflux(bounds%begp:bounds%endp), stat=ier)
    if (ier /= 0) then
          write(iulog,*)subname,' allocation error for prod100_nflux'; call endrun()
    end if

    if ( use_c13 ) then
       allocate(dwt_leafc13_seed(bounds%begp:bounds%endp), stat=ier)
       if (ier /= 0) then
          write(iulog,*)subname,' allocation error for dwt_leafc13_seed'; call endrun()
       end if
       allocate(dwt_deadstemc13_seed(bounds%begp:bounds%endp), stat=ier)
       if (ier /= 0) then
          write(iulog,*)subname,' allocation error for dwt_deadstemc13_seed'; call endrun()
       end if
       allocate(dwt_frootc13_to_litter(bounds%begp:bounds%endp), stat=ier)
       if (ier /= 0) then
          write(iulog,*)subname,' allocation error for dwt_frootc13_to_litter'; call endrun()
       end if
       allocate(dwt_livecrootc13_to_litter(bounds%begp:bounds%endp), stat=ier)
       if (ier /= 0) then
          write(iulog,*)subname,' allocation error for dwt_livecrootc13_to_litter'; call endrun()
       end if
       allocate(dwt_deadcrootc13_to_litter(bounds%begp:bounds%endp), stat=ier)
       if (ier /= 0) then
          write(iulog,*)subname,' allocation error for dwt_deadcrootc13_to_litter'; call endrun()
       end if
       allocate(conv_c13flux(bounds%begp:bounds%endp), stat=ier)
       if (ier /= 0) then
          write(iulog,*)subname,' allocation error for conv_c13flux'; call endrun()
       end if
       allocate(prod10_c13flux(bounds%begp:bounds%endp), stat=ier)
       if (ier /= 0) then
          write(iulog,*)subname,' allocation error for prod10_c13flux'; call endrun()
       end if
       allocate(prod100_c13flux(bounds%begp:bounds%endp), stat=ier)
       if (ier /= 0) then
          write(iulog,*)subname,' allocation error for prod100_c13flux'; call endrun()
       end if
    endif
    if ( use_c14 ) then
       allocate(dwt_leafc14_seed(bounds%begp:bounds%endp), stat=ier)
       if (ier /= 0) then
          write(iulog,*)subname,' allocation error for dwt_leafc14_seed'; call endrun()
       end if
       allocate(dwt_deadstemc14_seed(bounds%begp:bounds%endp), stat=ier)
       if (ier /= 0) then
          write(iulog,*)subname,' allocation error for dwt_deadstemc14_seed'; call endrun()
       end if
       allocate(dwt_frootc14_to_litter(bounds%begp:bounds%endp), stat=ier)
       if (ier /= 0) then
          write(iulog,*)subname,' allocation error for dwt_frootc14_to_litter'; call endrun()
       end if
       allocate(dwt_livecrootc14_to_litter(bounds%begp:bounds%endp), stat=ier)
       if (ier /= 0) then
          write(iulog,*)subname,' allocation error for dwt_livecrootc14_to_litter'; call endrun()
       end if
       allocate(dwt_deadcrootc14_to_litter(bounds%begp:bounds%endp), stat=ier)
       if (ier /= 0) then
          write(iulog,*)subname,' allocation error for dwt_deadcrootc14_to_litter'; call endrun()
       end if
       allocate(conv_c14flux(bounds%begp:bounds%endp), stat=ier)
       if (ier /= 0) then
          write(iulog,*)subname,' allocation error for conv_c14flux'; call endrun()
       end if
       allocate(prod10_c14flux(bounds%begp:bounds%endp), stat=ier)
       if (ier /= 0) then
          write(iulog,*)subname,' allocation error for prod10_c14flux'; call endrun()
       end if
       allocate(prod100_c14flux(bounds%begp:bounds%endp), stat=ier)
       if (ier /= 0) then
          write(iulog,*)subname,' allocation error for prod100_c14flux'; call endrun()
       end if
    endif
    
    ! Get time step
    dt = real( get_step_size(), r8 )
    
    do p = bounds%begp,bounds%endp
       c = pft%column(p)
       ! initialize all the pft-level local flux arrays
       dwt_leafc_seed(p) = 0._r8
       dwt_leafn_seed(p) = 0._r8
       dwt_deadstemc_seed(p) = 0._r8
       dwt_deadstemn_seed(p) = 0._r8
       dwt_frootc_to_litter(p) = 0._r8
       dwt_livecrootc_to_litter(p) = 0._r8
       dwt_deadcrootc_to_litter(p) = 0._r8
       dwt_frootn_to_litter(p) = 0._r8
       dwt_livecrootn_to_litter(p) = 0._r8
       dwt_deadcrootn_to_litter(p) = 0._r8
       conv_cflux(p) = 0._r8
       prod10_cflux(p) = 0._r8
       prod100_cflux(p) = 0._r8
       conv_nflux(p) = 0._r8
       prod10_nflux(p) = 0._r8
       prod100_nflux(p) = 0._r8
       
       if ( use_c13 ) then
          dwt_leafc13_seed(p) = 0._r8
          dwt_deadstemc13_seed(p) = 0._r8
          dwt_frootc13_to_litter(p) = 0._r8
          dwt_livecrootc13_to_litter(p) = 0._r8
          dwt_deadcrootc13_to_litter(p) = 0._r8
          conv_c13flux(p) = 0._r8
          prod10_c13flux(p) = 0._r8
          prod100_c13flux(p) = 0._r8
       endif
       
       if ( use_c14 ) then
          dwt_leafc14_seed(p) = 0._r8
          dwt_deadstemc14_seed(p) = 0._r8
          dwt_frootc14_to_litter(p) = 0._r8
          dwt_livecrootc14_to_litter(p) = 0._r8
          dwt_deadcrootc14_to_litter(p) = 0._r8
          conv_c14flux(p) = 0._r8
          prod10_c14flux(p) = 0._r8
          prod100_c14flux(p) = 0._r8
       endif
       
       l = pft%landunit(p)
       if (lun%itype(l) == istsoil .or. lun%itype(l) == istcrop) then
          
          ! calculate the change in weight for the timestep
          dwt = pft%wtcol(p)-wtcol_old(p)
            lfpftd(p)=-dwt
          ! PFTs for which weight increases on this timestep
          if (dwt > 0._r8) then
             
             ! first identify PFTs that are initiating on this timestep
             ! and set all the necessary state and flux variables
             if (wtcol_old(p) == 0._r8) then
                
                ! set initial conditions for PFT that is being initiated
                ! in this time step.  Based on the settings in cnIniTimeVar.
                
                ! pft-level carbon state variables
                pcs%leafc(p)              = 0._r8
                pcs%leafc_storage(p)      = 0._r8
                pcs%leafc_xfer(p)         = 0._r8
                pcs%frootc(p)             = 0._r8
                pcs%frootc_storage(p)     = 0._r8
                pcs%frootc_xfer(p)        = 0._r8
                pcs%livestemc(p)          = 0._r8
                pcs%livestemc_storage(p)  = 0._r8
                pcs%livestemc_xfer(p)     = 0._r8
                pcs%deadstemc(p)          = 0._r8
                pcs%deadstemc_storage(p)  = 0._r8
                pcs%deadstemc_xfer(p)     = 0._r8
                pcs%livecrootc(p)         = 0._r8
                pcs%livecrootc_storage(p) = 0._r8
                pcs%livecrootc_xfer(p)    = 0._r8
                pcs%deadcrootc(p)         = 0._r8
                pcs%deadcrootc_storage(p) = 0._r8
                pcs%deadcrootc_xfer(p)    = 0._r8
                pcs%gresp_storage(p)      = 0._r8
                pcs%gresp_xfer(p)         = 0._r8
                pcs%cpool(p)              = 0._r8
                pcs%xsmrpool(p)           = 0._r8
                pcs%pft_ctrunc(p)         = 0._r8
                pcs%dispvegc(p)           = 0._r8
                pcs%storvegc(p)           = 0._r8
                pcs%totvegc(p)            = 0._r8
                pcs%totpftc(p)            = 0._r8
                
                if ( use_c13 ) then
                   ! pft-level carbon-13 state variables
                   pc13s%leafc(p)              = 0._r8
                   pc13s%leafc_storage(p)      = 0._r8
                   pc13s%leafc_xfer(p)         = 0._r8
                   pc13s%frootc(p)             = 0._r8
                   pc13s%frootc_storage(p)     = 0._r8
                   pc13s%frootc_xfer(p)        = 0._r8
                   pc13s%livestemc(p)          = 0._r8
                   pc13s%livestemc_storage(p)  = 0._r8
                   pc13s%livestemc_xfer(p)     = 0._r8
                   pc13s%deadstemc(p)          = 0._r8
                   pc13s%deadstemc_storage(p)  = 0._r8
                   pc13s%deadstemc_xfer(p)     = 0._r8
                   pc13s%livecrootc(p)         = 0._r8
                   pc13s%livecrootc_storage(p) = 0._r8
                   pc13s%livecrootc_xfer(p)    = 0._r8
                   pc13s%deadcrootc(p)         = 0._r8
                   pc13s%deadcrootc_storage(p) = 0._r8
                   pc13s%deadcrootc_xfer(p)    = 0._r8
                   pc13s%gresp_storage(p)      = 0._r8
                   pc13s%gresp_xfer(p)         = 0._r8
                   pc13s%cpool(p)              = 0._r8
                   pc13s%xsmrpool(p)           = 0._r8
                   pc13s%pft_ctrunc(p)         = 0._r8
                   pc13s%dispvegc(p)           = 0._r8
                   pc13s%storvegc(p)           = 0._r8
                   pc13s%totvegc(p)            = 0._r8
                   pc13s%totpftc(p)            = 0._r8
                endif
                
                if ( use_c14 ) then
                   ! pft-level carbon-14 state variables
                   pc14s%leafc(p)              = 0._r8
                   pc14s%leafc_storage(p)      = 0._r8
                   pc14s%leafc_xfer(p)         = 0._r8
                   pc14s%frootc(p)             = 0._r8
                   pc14s%frootc_storage(p)     = 0._r8
                   pc14s%frootc_xfer(p)        = 0._r8
                   pc14s%livestemc(p)          = 0._r8
                   pc14s%livestemc_storage(p)  = 0._r8
                   pc14s%livestemc_xfer(p)     = 0._r8
                   pc14s%deadstemc(p)          = 0._r8
                   pc14s%deadstemc_storage(p)  = 0._r8
                   pc14s%deadstemc_xfer(p)     = 0._r8
                   pc14s%livecrootc(p)         = 0._r8
                   pc14s%livecrootc_storage(p) = 0._r8
                   pc14s%livecrootc_xfer(p)    = 0._r8
                   pc14s%deadcrootc(p)         = 0._r8
                   pc14s%deadcrootc_storage(p) = 0._r8
                   pc14s%deadcrootc_xfer(p)    = 0._r8
                   pc14s%gresp_storage(p)      = 0._r8
                   pc14s%gresp_xfer(p)         = 0._r8
                   pc14s%cpool(p)              = 0._r8
                   pc14s%xsmrpool(p)           = 0._r8
                   pc14s%pft_ctrunc(p)         = 0._r8
                   pc14s%dispvegc(p)           = 0._r8
                   pc14s%storvegc(p)           = 0._r8
                   pc14s%totvegc(p)            = 0._r8
                   pc14s%totpftc(p)            = 0._r8
                endif
                
                ! pft-level nitrogen state variables
                pns%leafn(p)	           = 0._r8
                pns%leafn_storage(p)      = 0._r8
                pns%leafn_xfer(p)         = 0._r8
                pns%frootn(p)	           = 0._r8
                pns%frootn_storage(p)     = 0._r8
                pns%frootn_xfer(p)        = 0._r8
                pns%livestemn(p)	       = 0._r8
                pns%livestemn_storage(p)  = 0._r8
                pns%livestemn_xfer(p)     = 0._r8
                pns%deadstemn(p)	       = 0._r8
                pns%deadstemn_storage(p)  = 0._r8
                pns%deadstemn_xfer(p)     = 0._r8
                pns%livecrootn(p)         = 0._r8
                pns%livecrootn_storage(p) = 0._r8
                pns%livecrootn_xfer(p)    = 0._r8
                pns%deadcrootn(p)         = 0._r8
                pns%deadcrootn_storage(p) = 0._r8
                pns%deadcrootn_xfer(p)    = 0._r8
                pns%retransn(p)	       = 0._r8
                pns%npool(p)	           = 0._r8
                pns%pft_ntrunc(p)         = 0._r8
                pns%dispvegn(p)           = 0._r8
                pns%storvegn(p)           = 0._r8
                pns%totvegn(p)            = 0._r8
                pns%totpftn (p)           = 0._r8
                
                ! initialize same flux and epv variables that are set
                ! in CNiniTimeVar
                pcf%psnsun(p) = 0._r8
                pcf%psnsha(p) = 0._r8
                pps%laisun(p) = 0._r8
                pps%laisha(p) = 0._r8
                
                pepv%dormant_flag(p) = 1._r8
                pepv%days_active(p) = 0._r8
                pepv%onset_flag(p) = 0._r8
                pepv%onset_counter(p) = 0._r8
                pepv%onset_gddflag(p) = 0._r8
                pepv%onset_fdd(p) = 0._r8
                pepv%onset_gdd(p) = 0._r8
                pepv%onset_swi(p) = 0.0_r8
                pepv%offset_flag(p) = 0._r8
                pepv%offset_counter(p) = 0._r8
                pepv%offset_fdd(p) = 0._r8
                pepv%offset_swi(p) = 0._r8
                pepv%lgsf(p) = 0._r8
                pepv%bglfr(p) = 0._r8
                pepv%bgtr(p) = 0._r8
                ! difference from CNiniTimeVar: using column-level
                ! information to initialize annavg_t2m.
                pepv%annavg_t2m(p) = cps%cannavg_t2m(c)
                pepv%tempavg_t2m(p) = 0._r8
                pepv%gpp(p) = 0._r8
                pepv%availc(p) = 0._r8
                pepv%xsmrpool_recover(p) = 0._r8
                pepv%alloc_pnow(p) = 1._r8
                pepv%c_allometry(p) = 0._r8
                pepv%n_allometry(p) = 0._r8
                pepv%plant_ndemand(p) = 0._r8
                pepv%tempsum_potential_gpp(p) = 0._r8
                pepv%annsum_potential_gpp(p) = 0._r8
                pepv%tempmax_retransn(p) = 0._r8
                pepv%annmax_retransn(p) = 0._r8
                pepv%avail_retransn(p) = 0._r8
                pepv%plant_nalloc(p) = 0._r8
                pepv%plant_calloc(p) = 0._r8
                pepv%excess_cflux(p) = 0._r8
                pepv%downreg(p) = 0._r8
                pepv%prev_leafc_to_litter(p) = 0._r8
                pepv%prev_frootc_to_litter(p) = 0._r8
                pepv%tempsum_npp(p) = 0._r8
                pepv%annsum_npp(p) = 0._r8
                
                if ( use_c13 ) then
                   pc13f%psnsun(p) = 0._r8
                   pc13f%psnsha(p) = 0._r8
                   
                   pps%alphapsnsun(p) = 0._r8
                   pps%alphapsnsha(p) = 0._r8
                   
                   pepv%xsmrpool_c13ratio(p) = c13ratio
                   
                   pepv%rc13_canair(p) = 0._r8
                   pepv%rc13_psnsun(p) = 0._r8
                   pepv%rc13_psnsha(p) = 0._r8
                endif
                
                if ( use_c14 ) then
                   pc14f%psnsun(p) = 0._r8
                   pc14f%psnsha(p) = 0._r8
                   pepv%rc14_atm(p) = c14ratio
                   pepv%rc14_atm(p) = 0._r8
                endif
                
             end if  ! end initialization of new pft
             
             ! (still in dwt > 0 block)
             
             ! set the seed sources for leaf and deadstem
             ! leaf source is split later between leaf, leaf_storage, leaf_xfer
             leafc_seed   = 0._r8
             leafn_seed   = 0._r8
             deadstemc_seed   = 0._r8
             deadstemn_seed   = 0._r8
             if ( use_c13 ) then
                leafc13_seed = 0._r8
                deadstemc13_seed = 0._r8
             endif
             if ( use_c14 ) then
                leafc14_seed = 0._r8
                deadstemc14_seed = 0._r8
             endif
             if (pft%itype(p) /= 0) then
                leafc_seed = 1._r8
                leafn_seed  = leafc_seed / pftcon%leafcn(pft%itype(p))
                if (pftcon%woody(pft%itype(p)) == 1._r8) then
                   deadstemc_seed = 0.1_r8
                   deadstemn_seed = deadstemc_seed / pftcon%deadwdcn(pft%itype(p))
                end if
                
                if ( use_c13 ) then
                   ! 13c state is initialized assuming del13c = -28 permil for C3, and -13 permil for C4.
                   ! That translates to ratios of (13c/(12c+13c)) of 0.01080455 for C3, and 0.01096945 for C4
                   ! based on the following formulae: 
                   ! r1 (13/12) = PDB + (del13c * PDB)/1000.0
                   ! r2 (13/(13+12)) = r1/(1+r1)
                   ! PDB = 0.0112372_R8  (ratio of 13C/12C in Pee Dee Belemnite, C isotope standard)
                   c3_del13c = -28._r8
                   c4_del13c = -13._r8
                   c3_r1_c13 = SHR_CONST_PDB + ((c3_del13c*SHR_CONST_PDB)/1000._r8)
                   c3_r2_c13 = c3_r1_c13/(1._r8 + c3_r1_c13)
                   c4_r1_c13 = SHR_CONST_PDB + ((c4_del13c*SHR_CONST_PDB)/1000._r8)
                   c4_r2_c13 = c4_r1_c13/(1._r8 + c4_r1_c13)
                   
                   if (pftcon%c3psn(pft%itype(p)) == 1._r8) then
                      leafc13_seed     = leafc_seed     * c3_r2_c13
                      deadstemc13_seed = deadstemc_seed * c3_r2_c13
                   else
                      leafc13_seed     = leafc_seed     * c4_r2_c13
                      deadstemc13_seed = deadstemc_seed * c4_r2_c13
                   end if
                endif
                
                if ( use_c14 ) then
                   ! 14c state is initialized assuming initial "modern" 14C of 1.e-12
                   if (pftcon%c3psn(pft%itype(p)) == 1._r8) then
                      leafc14_seed     = leafc_seed     * c14ratio
                      deadstemc14_seed = deadstemc_seed * c14ratio
                   else
                      leafc14_seed     = leafc_seed     * c14ratio
                      deadstemc14_seed = deadstemc_seed * c14ratio
                   end if
                endif
             end if
             
             ! When PFT area expands (dwt > 0), the pft-level mass density 
             ! is modified to conserve the original pft mass distributed
             ! over the new (larger) area, plus a term to account for the 
             ! introduction of new seed source for leaf and deadstem
             t1 = wtcol_old(p)/pft%wtcol(p)
             t2 = dwt/pft%wtcol(p)
             
             tot_leaf = pcs%leafc(p) + pcs%leafc_storage(p) + pcs%leafc_xfer(p)
             pleaf = 0._r8
             pstor = 0._r8
             pxfer = 0._r8
             if (tot_leaf /= 0._r8) then
                ! when adding seed source to non-zero leaf state, use current proportions
                pleaf = pcs%leafc(p)/tot_leaf
                pstor = pcs%leafc_storage(p)/tot_leaf
                pxfer = pcs%leafc_xfer(p)/tot_leaf
             else
                ! when initiating from zero leaf state, use evergreen flag to set proportions
                if (pftcon%evergreen(pft%itype(p)) == 1._r8) then
                   pleaf = 1._r8
                else
                   pstor = 1._r8
                end if
             end if
             pcs%leafc(p)         = pcs%leafc(p)*t1         + leafc_seed*pleaf*t2
             pcs%leafc_storage(p) = pcs%leafc_storage(p)*t1 + leafc_seed*pstor*t2
             pcs%leafc_xfer(p)    = pcs%leafc_xfer(p)*t1    + leafc_seed*pxfer*t2
             pcs%frootc(p)  		   = pcs%frootc(p) 			* t1
             pcs%frootc_storage(p)     = pcs%frootc_storage(p) 	* t1
             pcs%frootc_xfer(p) 	   = pcs%frootc_xfer(p)		* t1
             pcs%livestemc(p)		   = pcs%livestemc(p)  		* t1
             pcs%livestemc_storage(p)  = pcs%livestemc_storage(p)  * t1
             pcs%livestemc_xfer(p)     = pcs%livestemc_xfer(p) 	* t1
             pcs%deadstemc(p)     = pcs%deadstemc(p)*t1     + deadstemc_seed*t2
             pcs%deadstemc_storage(p)  = pcs%deadstemc_storage(p)  * t1
             pcs%deadstemc_xfer(p)     = pcs%deadstemc_xfer(p) 	* t1
             pcs%livecrootc(p)  	   = pcs%livecrootc(p) 		* t1
             pcs%livecrootc_storage(p) = pcs%livecrootc_storage(p) * t1
             pcs%livecrootc_xfer(p)    = pcs%livecrootc_xfer(p)	* t1
             pcs%deadcrootc(p)  	   = pcs%deadcrootc(p) 		* t1
             pcs%deadcrootc_storage(p) = pcs%deadcrootc_storage(p) * t1
             pcs%deadcrootc_xfer(p)    = pcs%deadcrootc_xfer(p)	* t1
             pcs%gresp_storage(p)	   = pcs%gresp_storage(p)  	* t1
             pcs%gresp_xfer(p)  	   = pcs%gresp_xfer(p) 		* t1
             pcs%cpool(p)			   = pcs%cpool(p)  			* t1
             pcs%xsmrpool(p)		   = pcs%xsmrpool(p)			* t1
             pcs%pft_ctrunc(p)  	   = pcs%pft_ctrunc(p) 		* t1
             pcs%dispvegc(p)		   = pcs%dispvegc(p)			* t1
             pcs%storvegc(p)		   = pcs%storvegc(p)			* t1
             pcs%totvegc(p) 		   = pcs%totvegc(p)			* t1
             pcs%totpftc(p) 		   = pcs%totpftc(p)			* t1
             
             if ( use_c13 ) then
                ! pft-level carbon-13 state variables 
                tot_leaf = pc13s%leafc(p) + pc13s%leafc_storage(p) + pc13s%leafc_xfer(p)
                pleaf = 0._r8
                pstor = 0._r8
                pxfer = 0._r8
                if (tot_leaf /= 0._r8) then
                   pleaf = pc13s%leafc(p)/tot_leaf
                   pstor = pc13s%leafc_storage(p)/tot_leaf
                   pxfer = pc13s%leafc_xfer(p)/tot_leaf
                else
                   ! when initiating from zero leaf state, use evergreen flag to set proportions
                   if (pftcon%evergreen(pft%itype(p)) == 1._r8) then
                      pleaf = 1._r8
                   else
                      pstor = 1._r8
                   end if
                end if
                pc13s%leafc(p)         = pc13s%leafc(p)*t1         + leafc13_seed*pleaf*t2
                pc13s%leafc_storage(p) = pc13s%leafc_storage(p)*t1 + leafc13_seed*pstor*t2
                pc13s%leafc_xfer(p)    = pc13s%leafc_xfer(p)*t1    + leafc13_seed*pxfer*t2
                pc13s%frootc(p)			 = pc13s%frootc(p) 		* t1
                pc13s%frootc_storage(p)	         = pc13s%frootc_storage(p) 	* t1
                pc13s%frootc_xfer(p)		 = pc13s%frootc_xfer(p)		* t1
                pc13s%livestemc(p) 		 = pc13s%livestemc(p)  		* t1
                pc13s%livestemc_storage(p)          = pc13s%livestemc_storage(p)      * t1
                pc13s%livestemc_xfer(p)	         = pc13s%livestemc_xfer(p) 	* t1
                pc13s%deadstemc(p)                  = pc13s%deadstemc(p)*t1     + deadstemc13_seed*t2
                pc13s%deadstemc_storage(p)          = pc13s%deadstemc_storage(p)      * t1
                pc13s%deadstemc_xfer(p)	         = pc13s%deadstemc_xfer(p) 	* t1
                pc13s%livecrootc(p)		 = pc13s%livecrootc(p) 		* t1
                pc13s%livecrootc_storage(p)         = pc13s%livecrootc_storage(p)     * t1
                pc13s%livecrootc_xfer(p)	         = pc13s%livecrootc_xfer(p)	* t1
                pc13s%deadcrootc(p)		 = pc13s%deadcrootc(p) 		* t1
                pc13s%deadcrootc_storage(p)         = pc13s%deadcrootc_storage(p)     * t1
                pc13s%deadcrootc_xfer(p)	         = pc13s%deadcrootc_xfer(p)	* t1
                pc13s%gresp_storage(p) 	         = pc13s%gresp_storage(p)  	* t1
                pc13s%gresp_xfer(p)		 = pc13s%gresp_xfer(p) 		* t1
                pc13s%cpool(p) 			 = pc13s%cpool(p)  		* t1
                pc13s%xsmrpool(p)  		 = pc13s%xsmrpool(p)		* t1
                pc13s%pft_ctrunc(p)		 = pc13s%pft_ctrunc(p) 		* t1
                pc13s%dispvegc(p)  		 = pc13s%dispvegc(p)		* t1
                pc13s%storvegc(p)  		 = pc13s%storvegc(p)		* t1
                pc13s%totvegc(p)			 = pc13s%totvegc(p)		* t1
                pc13s%totpftc(p)			 = pc13s%totpftc(p)		* t1
                
             endif
             
             if ( use_c14 ) then
                ! pft-level carbon-14 state variables 
                tot_leaf = pc14s%leafc(p) + pc14s%leafc_storage(p) + pc14s%leafc_xfer(p)
                pleaf = 0._r8
                pstor = 0._r8
                pxfer = 0._r8
                if (tot_leaf /= 0._r8) then
                   pleaf = pc14s%leafc(p)/tot_leaf
                   pstor = pc14s%leafc_storage(p)/tot_leaf
                   pxfer = pc14s%leafc_xfer(p)/tot_leaf
                else
                   ! when initiating from zero leaf state, use evergreen flag to set proportions
                   if (pftcon%evergreen(pft%itype(p)) == 1._r8) then
                      pleaf = 1._r8
                   else
                      pstor = 1._r8
                   end if
                end if
                pc14s%leafc(p)         = pc14s%leafc(p)*t1         + leafc14_seed*pleaf*t2
                pc14s%leafc_storage(p) = pc14s%leafc_storage(p)*t1 + leafc14_seed*pstor*t2
                pc14s%leafc_xfer(p)    = pc14s%leafc_xfer(p)*t1    + leafc14_seed*pxfer*t2
                pc14s%frootc(p)			 = pc14s%frootc(p) 		* t1
                pc14s%frootc_storage(p)	         = pc14s%frootc_storage(p) 	* t1
                pc14s%frootc_xfer(p)		 = pc14s%frootc_xfer(p)		* t1
                pc14s%livestemc(p) 		 = pc14s%livestemc(p)  		* t1
                pc14s%livestemc_storage(p)          = pc14s%livestemc_storage(p)      * t1
                pc14s%livestemc_xfer(p)	         = pc14s%livestemc_xfer(p) 	* t1
                pc14s%deadstemc(p)                  = pc14s%deadstemc(p)*t1     + deadstemc14_seed*t2
                pc14s%deadstemc_storage(p)          = pc14s%deadstemc_storage(p)      * t1
                pc14s%deadstemc_xfer(p)	         = pc14s%deadstemc_xfer(p) 	* t1
                pc14s%livecrootc(p)		 = pc14s%livecrootc(p) 		* t1
                pc14s%livecrootc_storage(p)         = pc14s%livecrootc_storage(p)     * t1
                pc14s%livecrootc_xfer(p)	         = pc14s%livecrootc_xfer(p)	* t1
                pc14s%deadcrootc(p)		 = pc14s%deadcrootc(p) 		* t1
                pc14s%deadcrootc_storage(p)         = pc14s%deadcrootc_storage(p)     * t1
                pc14s%deadcrootc_xfer(p)	         = pc14s%deadcrootc_xfer(p)	* t1
                pc14s%gresp_storage(p) 	         = pc14s%gresp_storage(p)  	* t1
                pc14s%gresp_xfer(p)		 = pc14s%gresp_xfer(p) 		* t1
                pc14s%cpool(p) 			 = pc14s%cpool(p)  		* t1
                pc14s%xsmrpool(p)  		 = pc14s%xsmrpool(p)		* t1
                pc14s%pft_ctrunc(p)		 = pc14s%pft_ctrunc(p) 		* t1
                pc14s%dispvegc(p)  		 = pc14s%dispvegc(p)		* t1
                pc14s%storvegc(p)  		 = pc14s%storvegc(p)		* t1
                pc14s%totvegc(p)			 = pc14s%totvegc(p)		* t1
                pc14s%totpftc(p)			 = pc14s%totpftc(p)		* t1
             endif
             
             
             tot_leaf = pns%leafn(p) + pns%leafn_storage(p) + pns%leafn_xfer(p)
             pleaf = 0._r8
             pstor = 0._r8
             pxfer = 0._r8
             if (tot_leaf /= 0._r8) then
                pleaf = pns%leafn(p)/tot_leaf
                pstor = pns%leafn_storage(p)/tot_leaf
                pxfer = pns%leafn_xfer(p)/tot_leaf
             else
                ! when initiating from zero leaf state, use evergreen flag to set proportions
                if (pftcon%evergreen(pft%itype(p)) == 1._r8) then
                   pleaf = 1._r8
                else
                   pstor = 1._r8
                end if
             end if
             ! pft-level nitrogen state variables
             pns%leafn(p)         = pns%leafn(p)*t1         + leafn_seed*pleaf*t2
             pns%leafn_storage(p) = pns%leafn_storage(p)*t1 + leafn_seed*pstor*t2
             pns%leafn_xfer(p)    = pns%leafn_xfer(p)*t1    + leafn_seed*pxfer*t2
             pns%frootn(p)  		   = pns%frootn(p) 		* t1
             pns%frootn_storage(p)         = pns%frootn_storage(p) 	* t1
             pns%frootn_xfer(p) 	   = pns%frootn_xfer(p)		* t1
             pns%livestemn(p)		   = pns%livestemn(p)  		* t1
             pns%livestemn_storage(p)      = pns%livestemn_storage(p)      * t1
             pns%livestemn_xfer(p)         = pns%livestemn_xfer(p) 	* t1
             pns%deadstemn(p)              = pns%deadstemn(p)*t1     + deadstemn_seed*t2
             pns%deadstemn_storage(p)      = pns%deadstemn_storage(p)      * t1
             pns%deadstemn_xfer(p)         = pns%deadstemn_xfer(p) 	* t1
             pns%livecrootn(p)  	   = pns%livecrootn(p) 		* t1
             pns%livecrootn_storage(p)     = pns%livecrootn_storage(p)     * t1
             pns%livecrootn_xfer(p)        = pns%livecrootn_xfer(p)	* t1
             pns%deadcrootn(p)  	   = pns%deadcrootn(p) 		* t1
             pns%deadcrootn_storage(p)     = pns%deadcrootn_storage(p)     * t1
             pns%deadcrootn_xfer(p)        = pns%deadcrootn_xfer(p)        * t1
             pns%retransn(p)		   = pns%retransn(p)		* t1
             pns%npool(p)		   = pns%npool(p)  		* t1
             pns%pft_ntrunc(p)  	   = pns%pft_ntrunc(p)        	* t1
             pns%dispvegn(p)		   = pns%dispvegn(p)		* t1
             pns%storvegn(p)		   = pns%storvegn(p)		* t1
             pns%totvegn(p) 		   = pns%totvegn(p)		* t1
             pns%totpftn(p) 		   = pns%totpftn(p)		* t1
             
             ! update temporary seed source arrays
             ! These are calculated in terms of the required contributions from
             ! column-level seed source
             dwt_leafc_seed(p)   = leafc_seed   * dwt
             if ( use_c13 ) then
                dwt_leafc13_seed(p) = leafc13_seed * dwt
                dwt_deadstemc13_seed(p) = deadstemc13_seed * dwt
             endif
             if ( use_c14 ) then
                dwt_leafc14_seed(p) = leafc14_seed * dwt
                dwt_deadstemc14_seed(p) = deadstemc14_seed * dwt
             endif
             dwt_leafn_seed(p)   = leafn_seed   * dwt
             dwt_deadstemc_seed(p)   = deadstemc_seed   * dwt
             dwt_deadstemn_seed(p)   = deadstemn_seed   * dwt
             
          else if (dwt < 0._r8) then
             
             ! if the pft lost weight on the timestep, then the carbon and nitrogen state
             ! variables are directed to litter, CWD, and wood product pools.
             
             ! N.B. : the conv_cflux, prod10_cflux, and prod100_cflux fluxes are accumulated
             ! as negative values, but the fluxes for pft-to-litter are accumulated as 
             ! positive values
             
             ! set local weight variables for this pft
             wt_new = pft%wtcol(p)
             wt_old = wtcol_old(p)
             
             !---------------
             ! C state update
             !---------------
             
             ! leafc 
             ptr => pcs%leafc(p)
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
             ptr => pcs%leafc_storage(p)
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
             ptr => pcs%leafc_xfer(p)
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
             ptr => pcs%frootc(p)
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
             ptr => pcs%frootc_storage(p)
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
             ptr => pcs%frootc_xfer(p)
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
             ptr => pcs%livestemc(p)
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
             ptr => pcs%livestemc_storage(p)
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
             ptr => pcs%livestemc_xfer(p)
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
             ptr => pcs%deadstemc(p)
             init_state = ptr*wt_old
             change_state = ptr*dwt
             new_state = init_state+change_state
             if (wt_new /= 0._r8) then
                ptr = new_state/wt_new
                conv_cflux(p) = conv_cflux(p) + change_state*pconv(pft%itype(p))
                prod10_cflux(p) = prod10_cflux(p) + change_state*pprod10(pft%itype(p))
                prod100_cflux(p) = prod100_cflux(p) + change_state*pprod100(pft%itype(p))
             else
                ptr = 0._r8
                conv_cflux(p) = conv_cflux(p) - init_state*pconv(pft%itype(p))
                prod10_cflux(p) = prod10_cflux(p) - init_state*pprod10(pft%itype(p))
                prod100_cflux(p) = prod100_cflux(p) - init_state*pprod100(pft%itype(p))
             end if
             
             ! deadstemc_storage 
             ptr => pcs%deadstemc_storage(p)
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
             ptr => pcs%deadstemc_xfer(p)
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
             ptr => pcs%livecrootc(p)
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
             ptr => pcs%livecrootc_storage(p)
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
             ptr => pcs%livecrootc_xfer(p)
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
             ptr => pcs%deadcrootc(p)
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
             ptr => pcs%deadcrootc_storage(p)
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
             ptr => pcs%deadcrootc_xfer(p)
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
             ptr => pcs%gresp_storage(p)
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
             ptr => pcs%gresp_xfer(p)
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
             ptr => pcs%cpool(p)
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
             ptr => pcs%xsmrpool(p)
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
             ptr => pcs%pft_ctrunc(p)
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

             if ( use_c13 ) then
                !-------------------
                ! C13 state update
                !-------------------
                
                ! set pointers to the conversion and product pool fluxes for this pft
                ! dwt_ptr0 is reserved for local assignment to dwt_xxx_to_litter fluxes
                dwt_ptr1 => conv_c13flux(p)
                dwt_ptr2 => prod10_c13flux(p)
                dwt_ptr3 => prod100_c13flux(p)
                
                ! leafc 
                ptr => pc13s%leafc(p)
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
                ptr => pc13s%leafc_storage(p)
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
                ptr => pc13s%leafc_xfer(p)
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
                ptr => pc13s%frootc(p)
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
                ptr => pc13s%frootc_storage(p)
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
                ptr => pc13s%frootc_xfer(p)
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
                ptr => pc13s%livestemc(p)
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
                ptr => pc13s%livestemc_storage(p)
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
                ptr => pc13s%livestemc_xfer(p)
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
                ptr => pc13s%deadstemc(p)
                init_state = ptr*wt_old
                change_state = ptr*dwt
                new_state = init_state+change_state
                if (wt_new /= 0._r8) then
                   ptr = new_state/wt_new
                   dwt_ptr1 = dwt_ptr1 + change_state*pconv(pft%itype(p))
                   dwt_ptr2 = dwt_ptr2 + change_state*pprod10(pft%itype(p))
                   dwt_ptr3 = dwt_ptr3 + change_state*pprod100(pft%itype(p))
                else
                   ptr = 0._r8
                   dwt_ptr1 = dwt_ptr1 - init_state*pconv(pft%itype(p))
                   dwt_ptr2 = dwt_ptr2 - init_state*pprod10(pft%itype(p))
                   dwt_ptr3 = dwt_ptr3 - init_state*pprod100(pft%itype(p))
                end if
                
                ! deadstemc_storage 
                ptr => pc13s%deadstemc_storage(p)
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
                ptr => pc13s%deadstemc_xfer(p)
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
                ptr => pc13s%livecrootc(p)
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
                ptr => pc13s%livecrootc_storage(p)
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
                ptr => pc13s%livecrootc_xfer(p)
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
                ptr => pc13s%deadcrootc(p)
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
                ptr => pc13s%deadcrootc_storage(p)
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
                ptr => pc13s%deadcrootc_xfer(p)
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
                ptr => pc13s%gresp_storage(p)
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
                ptr => pc13s%gresp_xfer(p)
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
                ptr => pc13s%cpool(p)
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
                ptr => pc13s%pft_ctrunc(p)
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
                
             endif

             if ( use_c14 ) then
                !-------------------
                ! C14 state update
                !-------------------
                
                ! set pointers to the conversion and product pool fluxes for this pft
                ! dwt_ptr0 is reserved for local assignment to dwt_xxx_to_litter fluxes
                dwt_ptr1 => conv_c14flux(p)
                dwt_ptr2 => prod10_c14flux(p)
                dwt_ptr3 => prod100_c14flux(p)
                
                ! leafc 
                ptr => pc14s%leafc(p)
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
                ptr => pc14s%leafc_storage(p)
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
                ptr => pc14s%leafc_xfer(p)
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
                ptr => pc14s%frootc(p)
                dwt_ptr0 => dwt_frootc14_to_litter(p)
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
                ptr => pc14s%frootc_storage(p)
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
                ptr => pc14s%frootc_xfer(p)
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
                ptr => pc14s%livestemc(p)
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
                ptr => pc14s%livestemc_storage(p)
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
                ptr => pc14s%livestemc_xfer(p)
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
                ptr => pc14s%deadstemc(p)
                init_state = ptr*wt_old
                change_state = ptr*dwt
                new_state = init_state+change_state
                if (wt_new /= 0._r8) then
                   ptr = new_state/wt_new
                   dwt_ptr1 = dwt_ptr1 + change_state*pconv(pft%itype(p))
                   dwt_ptr2 = dwt_ptr2 + change_state*pprod10(pft%itype(p))
                   dwt_ptr3 = dwt_ptr3 + change_state*pprod100(pft%itype(p))
                else
                   ptr = 0._r8
                   dwt_ptr1 = dwt_ptr1 - init_state*pconv(pft%itype(p))
                   dwt_ptr2 = dwt_ptr2 - init_state*pprod10(pft%itype(p))
                   dwt_ptr3 = dwt_ptr3 - init_state*pprod100(pft%itype(p))
                end if
                
                ! deadstemc_storage 
                ptr => pc14s%deadstemc_storage(p)
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
                ptr => pc14s%deadstemc_xfer(p)
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
                ptr => pc14s%livecrootc(p)
                dwt_ptr0 => dwt_livecrootc14_to_litter(p)
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
                ptr => pc14s%livecrootc_storage(p)
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
                ptr => pc14s%livecrootc_xfer(p)
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
                ptr => pc14s%deadcrootc(p)
                dwt_ptr0 => dwt_deadcrootc14_to_litter(p)
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
                ptr => pc14s%deadcrootc_storage(p)
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
                ptr => pc14s%deadcrootc_xfer(p)
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
                ptr => pc14s%gresp_storage(p)
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
                ptr => pc14s%gresp_xfer(p)
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
                ptr => pc14s%cpool(p)
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
                ptr => pc14s%pft_ctrunc(p)
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
             endif
             
             
             !---------------
             ! N state update
             !---------------
             
             ! set pointers to the conversion and product pool fluxes for this pft
             ! dwt_ptr0 is reserved for local assignment to dwt_xxx_to_litter fluxes
             dwt_ptr1 => conv_nflux(p)
             dwt_ptr2 => prod10_nflux(p)
             dwt_ptr3 => prod100_nflux(p)
             
             ! leafn 
             ptr => pns%leafn(p)
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
             ptr => pns%leafn_storage(p)
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
             ptr => pns%leafn_xfer(p)
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
             ptr => pns%frootn(p)
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
             ptr => pns%frootn_storage(p)
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
             ptr => pns%frootn_xfer(p)
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
             ptr => pns%livestemn(p)
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
             ptr => pns%livestemn_storage(p)
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
             ptr => pns%livestemn_xfer(p)
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
             ptr => pns%deadstemn(p)
             init_state = ptr*wt_old
             change_state = ptr*dwt
             new_state = init_state+change_state
             if (wt_new /= 0._r8) then
                ptr = new_state/wt_new
                dwt_ptr1 = dwt_ptr1 + change_state*pconv(pft%itype(p))
                dwt_ptr2 = dwt_ptr2 + change_state*pprod10(pft%itype(p))
                dwt_ptr3 = dwt_ptr3 + change_state*pprod100(pft%itype(p))
             else
                ptr = 0._r8
                dwt_ptr1 = dwt_ptr1 - init_state*pconv(pft%itype(p))
                dwt_ptr2 = dwt_ptr2 - init_state*pprod10(pft%itype(p))
                dwt_ptr3 = dwt_ptr3 - init_state*pprod100(pft%itype(p))
             end if
             
             ! deadstemn_storage 
             ptr => pns%deadstemn_storage(p)
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
             ptr => pns%deadstemn_xfer(p)
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
             ptr => pns%livecrootn(p)
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
             ptr => pns%livecrootn_storage(p)
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
             ptr => pns%livecrootn_xfer(p)
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
             ptr => pns%deadcrootn(p)
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
             ptr => pns%deadcrootn_storage(p)
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
             ptr => pns%deadcrootn_xfer(p)
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
             ptr => pns%retransn(p)
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
             ptr => pns%npool(p)
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
             ptr => pns%pft_ntrunc(p)
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
    do pi = 1,max_pft_per_col
       do c = bounds%begc, bounds%endc
          if ( pi <=  col%npfts(c) ) then
             p = col%pfti(c) + pi - 1
             
             ! C fluxes
             ccf%dwt_seedc_to_leaf(c) = ccf%dwt_seedc_to_leaf(c) + dwt_leafc_seed(p)/dt
             ccf%dwt_seedc_to_deadstem(c) = ccf%dwt_seedc_to_deadstem(c) &
                  + dwt_deadstemc_seed(p)/dt
             
             if ( use_c13 ) then
                cc13f%dwt_seedc_to_leaf(c) = cc13f%dwt_seedc_to_leaf(c) + dwt_leafc13_seed(p)/dt
                cc13f%dwt_seedc_to_deadstem(c) = cc13f%dwt_seedc_to_deadstem(c) &
                     + dwt_deadstemc13_seed(p)/dt
             endif
             
             if ( use_c14 ) then	
                cc14f%dwt_seedc_to_leaf(c) = cc14f%dwt_seedc_to_leaf(c) + dwt_leafc14_seed(p)/dt
                cc14f%dwt_seedc_to_deadstem(c) = cc14f%dwt_seedc_to_deadstem(c) &
                     + dwt_deadstemc14_seed(p)/dt
             endif
             
             ! N fluxes
             cnf%dwt_seedn_to_leaf(c) = cnf%dwt_seedn_to_leaf(c) + dwt_leafn_seed(p)/dt
             cnf%dwt_seedn_to_deadstem(c) = cnf%dwt_seedn_to_deadstem(c) &
                  + dwt_deadstemn_seed(p)/dt
          end if
       end do
    end do
    
    
    ! calculate pft-to-column for fluxes into litter and CWD pools
    do j = 1, nlevdecomp
       do pi = 1,max_pft_per_col
          do c = bounds%begc, bounds%endc
             if ( pi <=  col%npfts(c) ) then
                p = col%pfti(c) + pi - 1
                
                ! fine root litter carbon fluxes
                ccf%dwt_frootc_to_litr_met_c(c,j) = ccf%dwt_frootc_to_litr_met_c(c,j) + &
                     (dwt_frootc_to_litter(p)*pftcon%fr_flab(pft%itype(p)))/dt * pps%froot_prof(p,j)
                ccf%dwt_frootc_to_litr_cel_c(c,j) = ccf%dwt_frootc_to_litr_cel_c(c,j) + &
                     (dwt_frootc_to_litter(p)*pftcon%fr_fcel(pft%itype(p)))/dt * pps%froot_prof(p,j)
                ccf%dwt_frootc_to_litr_lig_c(c,j) = ccf%dwt_frootc_to_litr_lig_c(c,j) + &
                     (dwt_frootc_to_litter(p)*pftcon%fr_flig(pft%itype(p)))/dt * pps%froot_prof(p,j)
                
                
                ! fine root litter nitrogen fluxes
                cnf%dwt_frootn_to_litr_met_n(c,j) = cnf%dwt_frootn_to_litr_met_n(c,j) + &
                     (dwt_frootn_to_litter(p)*pftcon%fr_flab(pft%itype(p)))/dt * pps%froot_prof(p,j)
                cnf%dwt_frootn_to_litr_cel_n(c,j) = cnf%dwt_frootn_to_litr_cel_n(c,j) + &
                     (dwt_frootn_to_litter(p)*pftcon%fr_fcel(pft%itype(p)))/dt * pps%froot_prof(p,j)
                cnf%dwt_frootn_to_litr_lig_n(c,j) = cnf%dwt_frootn_to_litr_lig_n(c,j) + &
                     (dwt_frootn_to_litter(p)*pftcon%fr_flig(pft%itype(p)))/dt * pps%froot_prof(p,j)
                
                ! livecroot fluxes to cwd
                ccf%dwt_livecrootc_to_cwdc(c,j) = ccf%dwt_livecrootc_to_cwdc(c,j) + &
                     (dwt_livecrootc_to_litter(p))/dt * pps%croot_prof(p,j)
                cnf%dwt_livecrootn_to_cwdn(c,j) = cnf%dwt_livecrootn_to_cwdn(c,j) + &
                     (dwt_livecrootn_to_litter(p))/dt * pps%croot_prof(p,j)
                
                ! deadcroot fluxes to cwd
                ccf%dwt_deadcrootc_to_cwdc(c,j) = ccf%dwt_deadcrootc_to_cwdc(c,j) + &
                     (dwt_deadcrootc_to_litter(p))/dt * pps%croot_prof(p,j)
                cnf%dwt_deadcrootn_to_cwdn(c,j) = cnf%dwt_deadcrootn_to_cwdn(c,j) + &
                     (dwt_deadcrootn_to_litter(p))/dt * pps%croot_prof(p,j)
             
                if ( use_c13 ) then
                   ! C13 fine root litter fluxes
                   cc13f%dwt_frootc_to_litr_met_c(c,j) = cc13f%dwt_frootc_to_litr_met_c(c,j) + &
                        (dwt_frootc13_to_litter(p)*pftcon%fr_flab(pft%itype(p)))/dt * pps%froot_prof(p,j)
                   cc13f%dwt_frootc_to_litr_cel_c(c,j) = cc13f%dwt_frootc_to_litr_cel_c(c,j) + &
                        (dwt_frootc13_to_litter(p)*pftcon%fr_fcel(pft%itype(p)))/dt * pps%froot_prof(p,j)
                   cc13f%dwt_frootc_to_litr_lig_c(c,j) = cc13f%dwt_frootc_to_litr_lig_c(c,j) + &
                        (dwt_frootc13_to_litter(p)*pftcon%fr_flig(pft%itype(p)))/dt * pps%froot_prof(p,j)
                   ! livecroot fluxes to cwd
                   cc13f%dwt_livecrootc_to_cwdc(c,j) = cc13f%dwt_livecrootc_to_cwdc(c,j) + &
                        (dwt_livecrootc13_to_litter(p))/dt * pps%croot_prof(p,j)
                   ! deadcroot fluxes to cwd
                   cc13f%dwt_deadcrootc_to_cwdc(c,j) = cc13f%dwt_deadcrootc_to_cwdc(c,j) + &
                        (dwt_deadcrootc13_to_litter(p))/dt * pps%croot_prof(p,j)
                   
                endif
                
                if ( use_c14 ) then                   
                   ! C14 fine root litter fluxes
                   cc14f%dwt_frootc_to_litr_met_c(c,j) = cc14f%dwt_frootc_to_litr_met_c(c,j) + &
                        (dwt_frootc14_to_litter(p)*pftcon%fr_flab(pft%itype(p)))/dt * pps%froot_prof(p,j)
                   cc14f%dwt_frootc_to_litr_cel_c(c,j) = cc14f%dwt_frootc_to_litr_cel_c(c,j) + &
                        (dwt_frootc14_to_litter(p)*pftcon%fr_fcel(pft%itype(p)))/dt * pps%froot_prof(p,j)
                   cc14f%dwt_frootc_to_litr_lig_c(c,j) = cc14f%dwt_frootc_to_litr_lig_c(c,j) + &
                        (dwt_frootc14_to_litter(p)*pftcon%fr_flig(pft%itype(p)))/dt * pps%froot_prof(p,j)
                   ! livecroot fluxes to cwd
                   cc14f%dwt_livecrootc_to_cwdc(c,j) = cc14f%dwt_livecrootc_to_cwdc(c,j) + &
                        (dwt_livecrootc14_to_litter(p))/dt * pps%croot_prof(p,j)
                   ! deadcroot fluxes to cwd
                   cc14f%dwt_deadcrootc_to_cwdc(c,j) = cc14f%dwt_deadcrootc_to_cwdc(c,j) + &
                        (dwt_deadcrootc14_to_litter(p))/dt * pps%croot_prof(p,j)
                endif
                
             end if
          end do
       end do
    end do
    ! calculate pft-to-column for fluxes into product pools and conversion flux
    do pi = 1,max_pft_per_col
       do c = bounds%begc,bounds%endc
          if (pi <= col%npfts(c)) then
             p = col%pfti(c) + pi - 1
             
             ! column-level fluxes are accumulated as positive fluxes.
             ! column-level C flux updates
             ccf%dwt_conv_cflux(c) = ccf%dwt_conv_cflux(c) - conv_cflux(p)/dt
             ccf%dwt_prod10c_gain(c) = ccf%dwt_prod10c_gain(c) - prod10_cflux(p)/dt
             ccf%dwt_prod100c_gain(c) = ccf%dwt_prod100c_gain(c) - prod100_cflux(p)/dt

             ! These magic constants should be replaced with: nbrdlf_evr_trp_tree and nbrdlf_dcd_trp_tree
             if(pft%itype(p)==4.or.pft%itype(p)==6)then
                ccf%lf_conv_cflux(c) = ccf%lf_conv_cflux(c) - conv_cflux(p)/dt
             end if
             
             if ( use_c13 ) then
                ! C13 column-level flux updates
                cc13f%dwt_conv_cflux(c) = cc13f%dwt_conv_cflux(c) - conv_c13flux(p)/dt
                cc13f%dwt_prod10c_gain(c) = cc13f%dwt_prod10c_gain(c) - prod10_c13flux(p)/dt
                cc13f%dwt_prod100c_gain(c) = cc13f%dwt_prod100c_gain(c) - prod100_c13flux(p)/dt
             endif
             
             if ( use_c14 ) then
                ! C14 column-level flux updates
                cc14f%dwt_conv_cflux(c) = cc14f%dwt_conv_cflux(c) - conv_c14flux(p)/dt
                cc14f%dwt_prod10c_gain(c) = cc14f%dwt_prod10c_gain(c) - prod10_c14flux(p)/dt
                cc14f%dwt_prod100c_gain(c) = cc14f%dwt_prod100c_gain(c) - prod100_c14flux(p)/dt
             endif
             
             ! column-level N flux updates
             cnf%dwt_conv_nflux(c) = cnf%dwt_conv_nflux(c) - conv_nflux(p)/dt
             cnf%dwt_prod10n_gain(c) = cnf%dwt_prod10n_gain(c) - prod10_nflux(p)/dt
             cnf%dwt_prod100n_gain(c) = cnf%dwt_prod100n_gain(c) - prod100_nflux(p)/dt
             
          end if
       end do
    end do
    
    ! Deallocate pft-level flux arrays
    deallocate(dwt_leafc_seed)
    deallocate(dwt_leafn_seed)
    deallocate(dwt_deadstemc_seed)
    deallocate(dwt_deadstemn_seed)
    deallocate(dwt_frootc_to_litter)
    deallocate(dwt_livecrootc_to_litter)
    deallocate(dwt_deadcrootc_to_litter)
    deallocate(dwt_frootn_to_litter)
    deallocate(dwt_livecrootn_to_litter)
    deallocate(dwt_deadcrootn_to_litter)
    deallocate(conv_cflux)
    deallocate(prod10_cflux)
    deallocate(prod100_cflux)
    deallocate(conv_nflux)
    deallocate(prod10_nflux)
    deallocate(prod100_nflux)
             
    if ( use_c13 ) then
       deallocate(dwt_leafc13_seed)
       deallocate(dwt_deadstemc13_seed)
       deallocate(dwt_frootc13_to_litter)
       deallocate(dwt_livecrootc13_to_litter)
       deallocate(dwt_deadcrootc13_to_litter)
       deallocate(conv_c13flux)
       deallocate(prod10_c13flux)
       deallocate(prod100_c13flux)
    endif
             
    if ( use_c14 ) then
       deallocate(dwt_leafc14_seed)
       deallocate(dwt_deadstemc14_seed)
       deallocate(dwt_frootc14_to_litter)
       deallocate(dwt_livecrootc14_to_litter)
       deallocate(dwt_deadcrootc14_to_litter)
       deallocate(conv_c14flux)
       deallocate(prod10_c14flux)
       deallocate(prod100_c14flux)
    endif
    
    end associate 
   end subroutine pftdyn_cnbal

   !-----------------------------------------------------------------------
   subroutine pftwt_init(bounds)
     !
     ! !DESCRIPTION:
     ! Initialize time interpolation of cndv pft weights from annual to time step
     !
     ! !USES:
     use clm_varctl, only : nsrest, nsrStartup
     !
     ! !ARGUMENTS:
     implicit none
     type(bounds_type), intent(in) :: bounds  ! bounds
     !
     ! !LOCAL VARIABLES:
     integer  :: ier, p                        ! error status, do-loop index
     character(len=32) :: subname='pftwt_init' ! subroutine name
     !-----------------------------------------------------------------------

     allocate(wtcol_old(bounds%begp:bounds%endp),stat=ier)
     if (ier /= 0) then
        call endrun( subname//'::ERROR: pftwt_init allocation error for wtcol_old')
     end if

     if (nsrest == nsrStartup) then
        do p = bounds%begp,bounds%endp
           pdgvs%fpcgrid(p) = pft%wtcol(p)
           pdgvs%fpcgridold(p) = pft%wtcol(p)
           wtcol_old(p) = pft%wtcol(p)
        end do
     else
        do p = bounds%begp,bounds%endp
           wtcol_old(p) = pft%wtcol(p)
        end do
     end if

  end subroutine pftwt_init

  !-----------------------------------------------------------------------
  subroutine pftwt_interp( bounds )
    !
    ! !DESCRIPTION:
    ! Time interpolate cndv pft weights from annual to time step
    !
    ! !USES:
    use clm_time_manager, only : get_curr_calday, get_curr_date, get_days_per_year
    use clm_time_manager, only : get_step_size, get_nstep
    use clm_varcon      , only : istsoil ! CNDV incompatible with dynLU
    use clm_varctl      , only : finidat
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds  ! bounds
    !
    ! !LOCAL VARIABLES:
    integer  :: c,g,l,p            ! indices
    real(r8) :: cday               ! current calendar day (1.0 = 0Z on Jan 1)
    real(r8) :: wt1                ! time interpolation weights
    real(r8) :: dtime              ! model time step
    real(r8) :: days_per_year      ! days per year
    integer  :: nstep              ! time step number
    integer  :: year               ! year (0, ...) at nstep + 1
    integer  :: mon                ! month (1, ..., 12) at nstep + 1
    integer  :: day                ! day of month (1, ..., 31) at nstep + 1
    integer  :: sec                ! seconds into current date at nstep + 1
    character(len=32) :: subname='pftwt_interp' ! subroutine name
    !-----------------------------------------------------------------------

    ! Interpolate pft weight to current time step
    ! Map interpolated pctpft to subgrid weights
    ! assumes maxpatch_pft = numpft + 1, each landunit has 1 column, 
    ! SCAM not defined and create_croplandunit = .false.

    nstep         = get_nstep()
    dtime         = get_step_size()
    cday          = get_curr_calday(offset=-int(dtime))
    days_per_year = get_days_per_year()

    wt1 = ((days_per_year + 1._r8) - cday)/days_per_year

    call get_curr_date(year, mon, day, sec, offset=int(dtime))

    do p = bounds%begp,bounds%endp
       g = pft%gridcell(p)
       l = pft%landunit(p)

       if (lun%itype(l) == istsoil .and. lun%wtgcell(l) > 0._r8) then ! CNDV incompatible with dynLU
          wtcol_old(p)    = pft%wtcol(p)
          pft%wtcol(p)   = pdgvs%fpcgrid(p) + &
                     wt1 * (pdgvs%fpcgridold(p) - pdgvs%fpcgrid(p))
          pft%wtlunit(p) = pft%wtcol(p)
          pft%wtgcell(p) = pft%wtcol(p) * lun%wtgcell(l)

          if (mon==1 .and. day==1 .and. sec==dtime .and. nstep>0) then
             pdgvs%fpcgridold(p) = pdgvs%fpcgrid(p)
          end if
       end if
    end do

  end subroutine pftwt_interp

  !-----------------------------------------------------------------------
  subroutine CNHarvest (num_soilc, filter_soilc, num_soilp, filter_soilp)
    !
    ! !DESCRIPTION:
    ! Harvest mortality routine for coupled carbon-nitrogen code (CN)
    !
    ! !USES:
    use clmtype
    use pftvarcon       , only : noveg, nbrdlf_evr_shrub, pprodharv10
    use clm_varcon      , only : secspday
    use clm_time_manager, only : get_days_per_year
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: num_soilc       ! number of soil columns in filter
    integer, intent(in) :: filter_soilc(:) ! column filter for soil points
    integer, intent(in) :: num_soilp       ! number of soil pfts in filter
    integer, intent(in) :: filter_soilp(:) ! pft filter for soil points
    !
    ! !LOCAL VARIABLES:
    integer :: p                         ! pft index
    integer :: g                         ! gridcell index
    integer :: fp                        ! pft filter index
    real(r8):: am                        ! rate for fractional harvest mortality (1/yr)
    real(r8):: m                         ! rate for fractional harvest mortality (1/s)
    real(r8):: days_per_year             ! days per year
    !-----------------------------------------------------------------------

   associate(& 
   pgridcell                           =>   pft%gridcell                                 , & ! Input:  [integer (:)]  pft-level index into gridcell-level quantities     
   ivt                                 =>   pft%itype                                    , & ! Input:  [integer (:)]  pft vegetation type                                
   leafc                               =>    pcs%leafc                                   , & ! Input:  [real(r8) (:)]  (gC/m2) leaf C                                    
   frootc                              =>    pcs%frootc                                  , & ! Input:  [real(r8) (:)]  (gC/m2) fine root C                               
   livestemc                           =>    pcs%livestemc                               , & ! Input:  [real(r8) (:)]  (gC/m2) live stem C                               
   deadstemc                           =>    pcs%deadstemc                               , & ! Input:  [real(r8) (:)]  (gC/m2) dead stem C                               
   livecrootc                          =>    pcs%livecrootc                              , & ! Input:  [real(r8) (:)]  (gC/m2) live coarse root C                        
   deadcrootc                          =>    pcs%deadcrootc                              , & ! Input:  [real(r8) (:)]  (gC/m2) dead coarse root C                        
   xsmrpool                            =>    pcs%xsmrpool                                , & ! Input:  [real(r8) (:)]  (gC/m2) abstract C pool to meet excess MR demand  
   leafc_storage                       =>    pcs%leafc_storage                           , & ! Input:  [real(r8) (:)]  (gC/m2) leaf C storage                            
   frootc_storage                      =>    pcs%frootc_storage                          , & ! Input:  [real(r8) (:)]  (gC/m2) fine root C storage                       
   livestemc_storage                   =>    pcs%livestemc_storage                       , & ! Input:  [real(r8) (:)]  (gC/m2) live stem C storage                       
   deadstemc_storage                   =>    pcs%deadstemc_storage                       , & ! Input:  [real(r8) (:)]  (gC/m2) dead stem C storage                       
   livecrootc_storage                  =>    pcs%livecrootc_storage                      , & ! Input:  [real(r8) (:)]  (gC/m2) live coarse root C storage                
   deadcrootc_storage                  =>    pcs%deadcrootc_storage                      , & ! Input:  [real(r8) (:)]  (gC/m2) dead coarse root C storage                
   gresp_storage                       =>    pcs%gresp_storage                           , & ! Input:  [real(r8) (:)]  (gC/m2) growth respiration storage                
   leafc_xfer                          =>    pcs%leafc_xfer                              , & ! Input:  [real(r8) (:)]  (gC/m2) leaf C transfer                           
   frootc_xfer                         =>    pcs%frootc_xfer                             , & ! Input:  [real(r8) (:)]  (gC/m2) fine root C transfer                      
   livestemc_xfer                      =>    pcs%livestemc_xfer                          , & ! Input:  [real(r8) (:)]  (gC/m2) live stem C transfer                      
   deadstemc_xfer                      =>    pcs%deadstemc_xfer                          , & ! Input:  [real(r8) (:)]  (gC/m2) dead stem C transfer                      
   livecrootc_xfer                     =>    pcs%livecrootc_xfer                         , & ! Input:  [real(r8) (:)]  (gC/m2) live coarse root C transfer               
   deadcrootc_xfer                     =>    pcs%deadcrootc_xfer                         , & ! Input:  [real(r8) (:)]  (gC/m2) dead coarse root C transfer               
   gresp_xfer                          =>    pcs%gresp_xfer                              , & ! Input:  [real(r8) (:)]  (gC/m2) growth respiration transfer               
   leafn                               =>    pns%leafn                                   , & ! Input:  [real(r8) (:)]  (gN/m2) leaf N                                    
   frootn                              =>    pns%frootn                                  , & ! Input:  [real(r8) (:)]  (gN/m2) fine root N                               
   livestemn                           =>    pns%livestemn                               , & ! Input:  [real(r8) (:)]  (gN/m2) live stem N                               
   deadstemn                           =>    pns%deadstemn                               , & ! Input:  [real(r8) (:)]  (gN/m2) dead stem N                               
   livecrootn                          =>    pns%livecrootn                              , & ! Input:  [real(r8) (:)]  (gN/m2) live coarse root N                        
   deadcrootn                          =>    pns%deadcrootn                              , & ! Input:  [real(r8) (:)]  (gN/m2) dead coarse root N                        
   retransn                            =>    pns%retransn                                , & ! Input:  [real(r8) (:)]  (gN/m2) plant pool of retranslocated N            
   leafn_storage                       =>    pns%leafn_storage                           , & ! Input:  [real(r8) (:)]  (gN/m2) leaf N storage                            
   frootn_storage                      =>    pns%frootn_storage                          , & ! Input:  [real(r8) (:)]  (gN/m2) fine root N storage                       
   livestemn_storage                   =>    pns%livestemn_storage                       , & ! Input:  [real(r8) (:)]  (gN/m2) live stem N storage                       
   deadstemn_storage                   =>    pns%deadstemn_storage                       , & ! Input:  [real(r8) (:)]  (gN/m2) dead stem N storage                       
   livecrootn_storage                  =>    pns%livecrootn_storage                      , & ! Input:  [real(r8) (:)]  (gN/m2) live coarse root N storage                
   deadcrootn_storage                  =>    pns%deadcrootn_storage                      , & ! Input:  [real(r8) (:)]  (gN/m2) dead coarse root N storage                
   leafn_xfer                          =>    pns%leafn_xfer                              , & ! Input:  [real(r8) (:)]  (gN/m2) leaf N transfer                           
   frootn_xfer                         =>    pns%frootn_xfer                             , & ! Input:  [real(r8) (:)]  (gN/m2) fine root N transfer                      
   livestemn_xfer                      =>    pns%livestemn_xfer                          , & ! Input:  [real(r8) (:)]  (gN/m2) live stem N transfer                      
   deadstemn_xfer                      =>    pns%deadstemn_xfer                          , & ! Input:  [real(r8) (:)]  (gN/m2) dead stem N transfer                      
   livecrootn_xfer                     =>    pns%livecrootn_xfer                         , & ! Input:  [real(r8) (:)]  (gN/m2) live coarse root N transfer               
   deadcrootn_xfer                     =>    pns%deadcrootn_xfer                         , & ! Input:  [real(r8) (:)]  (gN/m2) dead coarse root N transfer               
   hrv_leafc_to_litter                 =>    pcf%hrv_leafc_to_litter                     , & ! Output: [real(r8) (:)]                                                    
   hrv_frootc_to_litter                =>    pcf%hrv_frootc_to_litter                    , & ! Output: [real(r8) (:)]                                                    
   hrv_livestemc_to_litter             =>    pcf%hrv_livestemc_to_litter                 , & ! Output: [real(r8) (:)]                                                    
   hrv_deadstemc_to_prod10c            =>    pcf%hrv_deadstemc_to_prod10c                , & ! Output: [real(r8) (:)]                                                    
   hrv_deadstemc_to_prod100c           =>    pcf%hrv_deadstemc_to_prod100c               , & ! Output: [real(r8) (:)]                                                    
   hrv_livecrootc_to_litter            =>    pcf%hrv_livecrootc_to_litter                , & ! Output: [real(r8) (:)]                                                    
   hrv_deadcrootc_to_litter            =>    pcf%hrv_deadcrootc_to_litter                , & ! Output: [real(r8) (:)]                                                    
   hrv_xsmrpool_to_atm                 =>    pcf%hrv_xsmrpool_to_atm                     , & ! Output: [real(r8) (:)]                                                    
   hrv_leafc_storage_to_litter         =>    pcf%hrv_leafc_storage_to_litter             , & ! Output: [real(r8) (:)]                                                    
   hrv_frootc_storage_to_litter        =>    pcf%hrv_frootc_storage_to_litter            , & ! Output: [real(r8) (:)]                                                    
   hrv_livestemc_storage_to_litter     =>    pcf%hrv_livestemc_storage_to_litter         , & ! Output: [real(r8) (:)]                                                    
   hrv_deadstemc_storage_to_litter     =>    pcf%hrv_deadstemc_storage_to_litter         , & ! Output: [real(r8) (:)]                                                    
   hrv_livecrootc_storage_to_litter    =>    pcf%hrv_livecrootc_storage_to_litter        , & ! Output: [real(r8) (:)]                                                    
   hrv_deadcrootc_storage_to_litter    =>    pcf%hrv_deadcrootc_storage_to_litter        , & ! Output: [real(r8) (:)]                                                    
   hrv_gresp_storage_to_litter         =>    pcf%hrv_gresp_storage_to_litter             , & ! Output: [real(r8) (:)]                                                    
   hrv_leafc_xfer_to_litter            =>    pcf%hrv_leafc_xfer_to_litter                , & ! Output: [real(r8) (:)]                                                    
   hrv_frootc_xfer_to_litter           =>    pcf%hrv_frootc_xfer_to_litter               , & ! Output: [real(r8) (:)]                                                    
   hrv_livestemc_xfer_to_litter        =>    pcf%hrv_livestemc_xfer_to_litter            , & ! Output: [real(r8) (:)]                                                    
   hrv_deadstemc_xfer_to_litter        =>    pcf%hrv_deadstemc_xfer_to_litter            , & ! Output: [real(r8) (:)]                                                    
   hrv_livecrootc_xfer_to_litter       =>    pcf%hrv_livecrootc_xfer_to_litter           , & ! Output: [real(r8) (:)]                                                    
   hrv_deadcrootc_xfer_to_litter       =>    pcf%hrv_deadcrootc_xfer_to_litter           , & ! Output: [real(r8) (:)]                                                    
   hrv_gresp_xfer_to_litter            =>    pcf%hrv_gresp_xfer_to_litter                , & ! Output: [real(r8) (:)]                                                    
   hrv_leafn_to_litter                 =>    pnf%hrv_leafn_to_litter                     , & ! Output: [real(r8) (:)]                                                    
   hrv_frootn_to_litter                =>    pnf%hrv_frootn_to_litter                    , & ! Output: [real(r8) (:)]                                                    
   hrv_livestemn_to_litter             =>    pnf%hrv_livestemn_to_litter                 , & ! Output: [real(r8) (:)]                                                    
   hrv_deadstemn_to_prod10n            =>    pnf%hrv_deadstemn_to_prod10n                , & ! Output: [real(r8) (:)]                                                    
   hrv_deadstemn_to_prod100n           =>    pnf%hrv_deadstemn_to_prod100n               , & ! Output: [real(r8) (:)]                                                    
   hrv_livecrootn_to_litter            =>    pnf%hrv_livecrootn_to_litter                , & ! Output: [real(r8) (:)]                                                    
   hrv_deadcrootn_to_litter            =>    pnf%hrv_deadcrootn_to_litter                , & ! Output: [real(r8) (:)]                                                    
   hrv_retransn_to_litter              =>    pnf%hrv_retransn_to_litter                  , & ! Output: [real(r8) (:)]                                                    
   hrv_leafn_storage_to_litter         =>    pnf%hrv_leafn_storage_to_litter             , & ! Output: [real(r8) (:)]                                                    
   hrv_frootn_storage_to_litter        =>    pnf%hrv_frootn_storage_to_litter            , & ! Output: [real(r8) (:)]                                                    
   hrv_livestemn_storage_to_litter     =>    pnf%hrv_livestemn_storage_to_litter         , & ! Output: [real(r8) (:)]                                                    
   hrv_deadstemn_storage_to_litter     =>    pnf%hrv_deadstemn_storage_to_litter         , & ! Output: [real(r8) (:)]                                                    
   hrv_livecrootn_storage_to_litter    =>    pnf%hrv_livecrootn_storage_to_litter        , & ! Output: [real(r8) (:)]                                                    
   hrv_deadcrootn_storage_to_litter    =>    pnf%hrv_deadcrootn_storage_to_litter        , & ! Output: [real(r8) (:)]                                                    
   hrv_leafn_xfer_to_litter            =>    pnf%hrv_leafn_xfer_to_litter                , & ! Output: [real(r8) (:)]                                                    
   hrv_frootn_xfer_to_litter           =>    pnf%hrv_frootn_xfer_to_litter               , & ! Output: [real(r8) (:)]                                                    
   hrv_livestemn_xfer_to_litter        =>    pnf%hrv_livestemn_xfer_to_litter            , & ! Output: [real(r8) (:)]                                                    
   hrv_deadstemn_xfer_to_litter        =>    pnf%hrv_deadstemn_xfer_to_litter            , & ! Output: [real(r8) (:)]                                                    
   hrv_livecrootn_xfer_to_litter       =>    pnf%hrv_livecrootn_xfer_to_litter           , & ! Output: [real(r8) (:)]                                                    
   hrv_deadcrootn_xfer_to_litter       =>    pnf%hrv_deadcrootn_xfer_to_litter             & ! Output: [real(r8) (:)]                                                    
   )


   days_per_year = get_days_per_year()

   ! pft loop
   do fp = 1,num_soilp
      p = filter_soilp(fp)
      g = pgridcell(p)
      
      ! If this is a tree pft, then
      ! get the annual harvest "mortality" rate (am) from harvest array
      ! and convert to rate per second
      if (ivt(p) > noveg .and. ivt(p) < nbrdlf_evr_shrub) then

         if (do_harvest) then
            am = harvest(g)
            m  = am/(days_per_year * secspday)
         else
            m = 0._r8
         end if   

         ! pft-level harvest carbon fluxes
         ! displayed pools
         hrv_leafc_to_litter(p)               = leafc(p)               * m
         hrv_frootc_to_litter(p)              = frootc(p)              * m
         hrv_livestemc_to_litter(p)           = livestemc(p)           * m
         hrv_deadstemc_to_prod10c(p)          = deadstemc(p)           * m * &
                                                pprodharv10(ivt(p))
         hrv_deadstemc_to_prod100c(p)         = deadstemc(p)           * m * &
                                                (1.0_r8 - pprodharv10(ivt(p)))
         hrv_livecrootc_to_litter(p)          = livecrootc(p)          * m
         hrv_deadcrootc_to_litter(p)          = deadcrootc(p)          * m
         hrv_xsmrpool_to_atm(p)               = xsmrpool(p)            * m

         ! storage pools
         hrv_leafc_storage_to_litter(p)       = leafc_storage(p)       * m
         hrv_frootc_storage_to_litter(p)      = frootc_storage(p)      * m
         hrv_livestemc_storage_to_litter(p)   = livestemc_storage(p)   * m
         hrv_deadstemc_storage_to_litter(p)   = deadstemc_storage(p)   * m
         hrv_livecrootc_storage_to_litter(p)  = livecrootc_storage(p)  * m
         hrv_deadcrootc_storage_to_litter(p)  = deadcrootc_storage(p)  * m
         hrv_gresp_storage_to_litter(p)       = gresp_storage(p)       * m

         ! transfer pools
         hrv_leafc_xfer_to_litter(p)          = leafc_xfer(p)          * m
         hrv_frootc_xfer_to_litter(p)         = frootc_xfer(p)         * m
         hrv_livestemc_xfer_to_litter(p)      = livestemc_xfer(p)      * m
         hrv_deadstemc_xfer_to_litter(p)      = deadstemc_xfer(p)      * m
         hrv_livecrootc_xfer_to_litter(p)     = livecrootc_xfer(p)     * m
         hrv_deadcrootc_xfer_to_litter(p)     = deadcrootc_xfer(p)     * m
         hrv_gresp_xfer_to_litter(p)          = gresp_xfer(p)          * m

         ! pft-level harvest mortality nitrogen fluxes
         ! displayed pools
         hrv_leafn_to_litter(p)               = leafn(p)               * m
         hrv_frootn_to_litter(p)              = frootn(p)              * m
         hrv_livestemn_to_litter(p)           = livestemn(p)           * m
         hrv_deadstemn_to_prod10n(p)          = deadstemn(p)           * m * &
                                                pprodharv10(ivt(p))
         hrv_deadstemn_to_prod100n(p)         = deadstemn(p)           * m * &
                                                (1.0_r8 - pprodharv10(ivt(p)))
         hrv_livecrootn_to_litter(p)          = livecrootn(p)          * m
         hrv_deadcrootn_to_litter(p)          = deadcrootn(p)          * m
         hrv_retransn_to_litter(p)            = retransn(p)            * m

         ! storage pools
         hrv_leafn_storage_to_litter(p)       = leafn_storage(p)       * m
         hrv_frootn_storage_to_litter(p)      = frootn_storage(p)      * m
         hrv_livestemn_storage_to_litter(p)   = livestemn_storage(p)   * m
         hrv_deadstemn_storage_to_litter(p)   = deadstemn_storage(p)   * m
         hrv_livecrootn_storage_to_litter(p)  = livecrootn_storage(p)  * m
         hrv_deadcrootn_storage_to_litter(p)  = deadcrootn_storage(p)  * m

         ! transfer pools
         hrv_leafn_xfer_to_litter(p)          = leafn_xfer(p)          * m
         hrv_frootn_xfer_to_litter(p)         = frootn_xfer(p)         * m
         hrv_livestemn_xfer_to_litter(p)      = livestemn_xfer(p)      * m
         hrv_deadstemn_xfer_to_litter(p)      = deadstemn_xfer(p)      * m
         hrv_livecrootn_xfer_to_litter(p)     = livecrootn_xfer(p)     * m
         hrv_deadcrootn_xfer_to_litter(p)     = deadcrootn_xfer(p)     * m
         
      end if  ! end tree block

   end do ! end of pft loop

   ! gather all pft-level litterfall fluxes from harvest to the column
   ! for litter C and N inputs

   call CNHarvestPftToColumn(num_soilc, filter_soilc)

    end associate 
 end subroutine CNHarvest

 !-----------------------------------------------------------------------
 subroutine CNHarvestPftToColumn (num_soilc, filter_soilc)
   !
   ! !DESCRIPTION:
   ! called at the end of CNHarvest to gather all pft-level harvest litterfall fluxes
   ! to the column level and assign them to the three litter pools
   !
   ! !USES:
   use clmtype
   use clm_varpar, only : max_pft_per_col, maxpatch_pft, nlevdecomp
   !
   ! !ARGUMENTS:
   implicit none
   integer, intent(in) :: num_soilc       ! number of soil columns in filter
   integer, intent(in) :: filter_soilc(:) ! soil column filter
   !
   ! !LOCAL VARIABLES:
   integer :: fc,c,pi,p,j               ! indices
   !-----------------------------------------------------------------------

   associate(& 
   lf_flab                             =>    pftcon%lf_flab                              , & ! Input:  [real(r8) (:)]  leaf litter labile fraction                       
   lf_fcel                             =>    pftcon%lf_fcel                              , & ! Input:  [real(r8) (:)]  leaf litter cellulose fraction                    
   lf_flig                             =>    pftcon%lf_flig                              , & ! Input:  [real(r8) (:)]  leaf litter lignin fraction                       
   fr_flab                             =>    pftcon%fr_flab                              , & ! Input:  [real(r8) (:)]  fine root litter labile fraction                  
   fr_fcel                             =>    pftcon%fr_fcel                              , & ! Input:  [real(r8) (:)]  fine root litter cellulose fraction               
   fr_flig                             =>    pftcon%fr_flig                              , & ! Input:  [real(r8) (:)]  fine root litter lignin fraction                  
   npfts                               =>   col%npfts                                    , & ! Input:  [integer (:)]  number of pfts for each column                     
   pfti                                =>   col%pfti                                     , & ! Input:  [integer (:)]  beginning pft index for each column                
   chrv_deadstemc_to_prod10c           =>    ccf%hrv_deadstemc_to_prod10c                , & ! InOut:  [real(r8) (:)]                                                    
   chrv_deadstemc_to_prod100c          =>    ccf%hrv_deadstemc_to_prod100c               , & ! InOut:  [real(r8) (:)]                                                    
   chrv_deadstemn_to_prod10n           =>    cnf%hrv_deadstemn_to_prod10n                , & ! InOut:  [real(r8) (:)]                                                    
   chrv_deadstemn_to_prod100n          =>    cnf%hrv_deadstemn_to_prod100n               , & ! InOut:  [real(r8) (:)]                                                    
   harvest_c_to_litr_met_c             =>    ccf%harvest_c_to_litr_met_c                 , & ! Input:  [real(r8) (:,:)]  C fluxes associated with harvest to litter metabolic pool (gC/m3/s)
   harvest_c_to_litr_cel_c             =>    ccf%harvest_c_to_litr_cel_c                 , & ! Input:  [real(r8) (:,:)]  C fluxes associated with harvest to litter cellulose pool (gC/m3/s)
   harvest_c_to_litr_lig_c             =>    ccf%harvest_c_to_litr_lig_c                 , & ! Input:  [real(r8) (:,:)]  C fluxes associated with harvest to litter lignin pool (gC/m3/s)
   harvest_c_to_cwdc                   =>    ccf%harvest_c_to_cwdc                       , & ! Input:  [real(r8) (:,:)]  C fluxes associated with harvest to CWD pool (gC/m3/s)
   harvest_n_to_litr_met_n             =>    cnf%harvest_n_to_litr_met_n                 , & ! Input:  [real(r8) (:,:)]  N fluxes associated with harvest to litter metabolic pool (gN/m3/s)
   harvest_n_to_litr_cel_n             =>    cnf%harvest_n_to_litr_cel_n                 , & ! Input:  [real(r8) (:,:)]  N fluxes associated with harvest to litter cellulose pool (gN/m3/s)
   harvest_n_to_litr_lig_n             =>    cnf%harvest_n_to_litr_lig_n                 , & ! Input:  [real(r8) (:,:)]  N fluxes associated with harvest to litter lignin pool (gN/m3/s)
   harvest_n_to_cwdn                   =>    cnf%harvest_n_to_cwdn                       , & ! Input:  [real(r8) (:,:)]  N fluxes associated with harvest to CWD pool (gN/m3/s)
   pactive                             =>    pft%active                                  , & ! Input:  [logical (:)]  true=>do computations on this pft (see reweightMod for details)
   ivt                                 =>   pft%itype                                    , & ! Input:  [integer (:)]  pft vegetation type                                
   wtcol                               =>   pft%wtcol                                    , & ! Input:  [real(r8) (:)]  pft weight relative to column (0-1)               
   hrv_leafc_to_litter                 =>    pcf%hrv_leafc_to_litter                     , & ! Input:  [real(r8) (:)]                                                    
   hrv_frootc_to_litter                =>    pcf%hrv_frootc_to_litter                    , & ! Input:  [real(r8) (:)]                                                    
   hrv_livestemc_to_litter             =>    pcf%hrv_livestemc_to_litter                 , & ! Input:  [real(r8) (:)]                                                    
   phrv_deadstemc_to_prod10c           =>    pcf%hrv_deadstemc_to_prod10c                , & ! Input:  [real(r8) (:)]                                                    
   phrv_deadstemc_to_prod100c          =>    pcf%hrv_deadstemc_to_prod100c               , & ! Input:  [real(r8) (:)]                                                    
   hrv_livecrootc_to_litter            =>    pcf%hrv_livecrootc_to_litter                , & ! Input:  [real(r8) (:)]                                                    
   hrv_deadcrootc_to_litter            =>    pcf%hrv_deadcrootc_to_litter                , & ! Input:  [real(r8) (:)]                                                    
   hrv_leafc_storage_to_litter         =>    pcf%hrv_leafc_storage_to_litter             , & ! Input:  [real(r8) (:)]                                                    
   hrv_frootc_storage_to_litter        =>    pcf%hrv_frootc_storage_to_litter            , & ! Input:  [real(r8) (:)]                                                    
   hrv_livestemc_storage_to_litter     =>    pcf%hrv_livestemc_storage_to_litter         , & ! Input:  [real(r8) (:)]                                                    
   hrv_deadstemc_storage_to_litter     =>    pcf%hrv_deadstemc_storage_to_litter         , & ! Input:  [real(r8) (:)]                                                    
   hrv_livecrootc_storage_to_litter    =>    pcf%hrv_livecrootc_storage_to_litter        , & ! Input:  [real(r8) (:)]                                                    
   hrv_deadcrootc_storage_to_litter    =>    pcf%hrv_deadcrootc_storage_to_litter        , & ! Input:  [real(r8) (:)]                                                    
   hrv_gresp_storage_to_litter         =>    pcf%hrv_gresp_storage_to_litter             , & ! Input:  [real(r8) (:)]                                                    
   hrv_leafc_xfer_to_litter            =>    pcf%hrv_leafc_xfer_to_litter                , & ! Input:  [real(r8) (:)]                                                    
   hrv_frootc_xfer_to_litter           =>    pcf%hrv_frootc_xfer_to_litter               , & ! Input:  [real(r8) (:)]                                                    
   hrv_livestemc_xfer_to_litter        =>    pcf%hrv_livestemc_xfer_to_litter            , & ! Input:  [real(r8) (:)]                                                    
   hrv_deadstemc_xfer_to_litter        =>    pcf%hrv_deadstemc_xfer_to_litter            , & ! Input:  [real(r8) (:)]                                                    
   hrv_livecrootc_xfer_to_litter       =>    pcf%hrv_livecrootc_xfer_to_litter           , & ! Input:  [real(r8) (:)]                                                    
   hrv_deadcrootc_xfer_to_litter       =>    pcf%hrv_deadcrootc_xfer_to_litter           , & ! Input:  [real(r8) (:)]                                                    
   hrv_gresp_xfer_to_litter            =>    pcf%hrv_gresp_xfer_to_litter                , & ! Input:  [real(r8) (:)]                                                    
   hrv_leafn_to_litter                 =>    pnf%hrv_leafn_to_litter                     , & ! Input:  [real(r8) (:)]                                                    
   hrv_frootn_to_litter                =>    pnf%hrv_frootn_to_litter                    , & ! Input:  [real(r8) (:)]                                                    
   hrv_livestemn_to_litter             =>    pnf%hrv_livestemn_to_litter                 , & ! Input:  [real(r8) (:)]                                                    
   phrv_deadstemn_to_prod10n           =>    pnf%hrv_deadstemn_to_prod10n                , & ! Input:  [real(r8) (:)]                                                    
   phrv_deadstemn_to_prod100n          =>    pnf%hrv_deadstemn_to_prod100n               , & ! Input:  [real(r8) (:)]                                                    
   hrv_livecrootn_to_litter            =>    pnf%hrv_livecrootn_to_litter                , & ! Input:  [real(r8) (:)]                                                    
   hrv_deadcrootn_to_litter            =>    pnf%hrv_deadcrootn_to_litter                , & ! Input:  [real(r8) (:)]                                                    
   hrv_retransn_to_litter              =>    pnf%hrv_retransn_to_litter                  , & ! Input:  [real(r8) (:)]                                                    
   hrv_leafn_storage_to_litter         =>    pnf%hrv_leafn_storage_to_litter             , & ! Input:  [real(r8) (:)]                                                    
   hrv_frootn_storage_to_litter        =>    pnf%hrv_frootn_storage_to_litter            , & ! Input:  [real(r8) (:)]                                                    
   hrv_livestemn_storage_to_litter     =>    pnf%hrv_livestemn_storage_to_litter         , & ! Input:  [real(r8) (:)]                                                    
   hrv_deadstemn_storage_to_litter     =>    pnf%hrv_deadstemn_storage_to_litter         , & ! Input:  [real(r8) (:)]                                                    
   hrv_livecrootn_storage_to_litter    =>    pnf%hrv_livecrootn_storage_to_litter        , & ! Input:  [real(r8) (:)]                                                    
   hrv_deadcrootn_storage_to_litter    =>    pnf%hrv_deadcrootn_storage_to_litter        , & ! Input:  [real(r8) (:)]                                                    
   hrv_leafn_xfer_to_litter            =>    pnf%hrv_leafn_xfer_to_litter                , & ! Input:  [real(r8) (:)]                                                    
   hrv_frootn_xfer_to_litter           =>    pnf%hrv_frootn_xfer_to_litter               , & ! Input:  [real(r8) (:)]                                                    
   hrv_livestemn_xfer_to_litter        =>    pnf%hrv_livestemn_xfer_to_litter            , & ! Input:  [real(r8) (:)]                                                    
   hrv_deadstemn_xfer_to_litter        =>    pnf%hrv_deadstemn_xfer_to_litter            , & ! Input:  [real(r8) (:)]                                                    
   hrv_livecrootn_xfer_to_litter       =>    pnf%hrv_livecrootn_xfer_to_litter           , & ! Input:  [real(r8) (:)]                                                    
   hrv_deadcrootn_xfer_to_litter       =>    pnf%hrv_deadcrootn_xfer_to_litter           , & ! Input:  [real(r8) (:)]                                                    
   leaf_prof                           =>    pps%leaf_prof                               , & ! InOut:  [real(r8) (:,:)]  (1/m) profile of leaves                         
   froot_prof                          =>    pps%froot_prof                              , & ! InOut:  [real(r8) (:,:)]  (1/m) profile of fine roots                     
   croot_prof                          =>    pps%croot_prof                              , & ! InOut:  [real(r8) (:,:)]  (1/m) profile of coarse roots                   
   stem_prof                           =>    pps%stem_prof                                 & ! InOut:  [real(r8) (:,:)]  (1/m) profile of stems                          
   )

   do j = 1, nlevdecomp
      do pi = 1,maxpatch_pft
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            
            if (pi <=  npfts(c)) then
               p = pfti(c) + pi - 1
               
               if (pactive(p)) then
                  
                  
                  ! leaf harvest mortality carbon fluxes
                  harvest_c_to_litr_met_c(c,j) = harvest_c_to_litr_met_c(c,j) + &
                       hrv_leafc_to_litter(p) * lf_flab(ivt(p)) * wtcol(p) * leaf_prof(p,j)
                  harvest_c_to_litr_cel_c(c,j) = harvest_c_to_litr_cel_c(c,j) + &
                       hrv_leafc_to_litter(p) * lf_fcel(ivt(p)) * wtcol(p) * leaf_prof(p,j)
                  harvest_c_to_litr_lig_c(c,j) = harvest_c_to_litr_lig_c(c,j) + &
                       hrv_leafc_to_litter(p) * lf_flig(ivt(p)) * wtcol(p) * leaf_prof(p,j)
                  
                  ! fine root harvest mortality carbon fluxes
                  harvest_c_to_litr_met_c(c,j) = harvest_c_to_litr_met_c(c,j) + &
                       hrv_frootc_to_litter(p) * fr_flab(ivt(p)) * wtcol(p) * froot_prof(p,j)
                  harvest_c_to_litr_cel_c(c,j) = harvest_c_to_litr_cel_c(c,j) + &
                       hrv_frootc_to_litter(p) * fr_fcel(ivt(p)) * wtcol(p) * froot_prof(p,j)
                  harvest_c_to_litr_lig_c(c,j) = harvest_c_to_litr_lig_c(c,j) + &
                       hrv_frootc_to_litter(p) * fr_flig(ivt(p)) * wtcol(p) * froot_prof(p,j)
                  
                  ! wood harvest mortality carbon fluxes
                  harvest_c_to_cwdc(c,j)  = harvest_c_to_cwdc(c,j)  + &
                       hrv_livestemc_to_litter(p)  * wtcol(p) * stem_prof(p,j) 
                  harvest_c_to_cwdc(c,j) = harvest_c_to_cwdc(c,j) + &
                       hrv_livecrootc_to_litter(p) * wtcol(p) * croot_prof(p,j)
                  harvest_c_to_cwdc(c,j) = harvest_c_to_cwdc(c,j) + &
                       hrv_deadcrootc_to_litter(p) * wtcol(p) * croot_prof(p,j) 
                  
                  ! storage harvest mortality carbon fluxes
                  harvest_c_to_litr_met_c(c,j)      = harvest_c_to_litr_met_c(c,j)      + &
                       hrv_leafc_storage_to_litter(p)      * wtcol(p) * leaf_prof(p,j)
                  harvest_c_to_litr_met_c(c,j)     = harvest_c_to_litr_met_c(c,j)     + &
                       hrv_frootc_storage_to_litter(p)     * wtcol(p) * froot_prof(p,j)
                  harvest_c_to_litr_met_c(c,j)  = harvest_c_to_litr_met_c(c,j)  + &
                       hrv_livestemc_storage_to_litter(p)  * wtcol(p) * stem_prof(p,j)
                  harvest_c_to_litr_met_c(c,j)  = harvest_c_to_litr_met_c(c,j)  + &
                       hrv_deadstemc_storage_to_litter(p)  * wtcol(p) * stem_prof(p,j)
                  harvest_c_to_litr_met_c(c,j) = harvest_c_to_litr_met_c(c,j) + &
                       hrv_livecrootc_storage_to_litter(p) * wtcol(p) * croot_prof(p,j)
                  harvest_c_to_litr_met_c(c,j) = harvest_c_to_litr_met_c(c,j) + &
                       hrv_deadcrootc_storage_to_litter(p) * wtcol(p) * croot_prof(p,j)
                  harvest_c_to_litr_met_c(c,j)      = harvest_c_to_litr_met_c(c,j)      + &
                       hrv_gresp_storage_to_litter(p)      * wtcol(p) * leaf_prof(p,j)
                  
                  ! transfer harvest mortality carbon fluxes
                  harvest_c_to_litr_met_c(c,j)      = harvest_c_to_litr_met_c(c,j)      + &
                       hrv_leafc_xfer_to_litter(p)      * wtcol(p) * leaf_prof(p,j)
                  harvest_c_to_litr_met_c(c,j)     = harvest_c_to_litr_met_c(c,j)     + &
                       hrv_frootc_xfer_to_litter(p)     * wtcol(p) * froot_prof(p,j)
                  harvest_c_to_litr_met_c(c,j)  = harvest_c_to_litr_met_c(c,j)  + &
                       hrv_livestemc_xfer_to_litter(p)  * wtcol(p) * stem_prof(p,j)
                  harvest_c_to_litr_met_c(c,j)  = harvest_c_to_litr_met_c(c,j)  + &
                       hrv_deadstemc_xfer_to_litter(p)  * wtcol(p) * stem_prof(p,j)
                  harvest_c_to_litr_met_c(c,j) = harvest_c_to_litr_met_c(c,j) + &
                       hrv_livecrootc_xfer_to_litter(p) * wtcol(p) * croot_prof(p,j)
                  harvest_c_to_litr_met_c(c,j) = harvest_c_to_litr_met_c(c,j) + &
                       hrv_deadcrootc_xfer_to_litter(p) * wtcol(p) * croot_prof(p,j)
                  harvest_c_to_litr_met_c(c,j)      = harvest_c_to_litr_met_c(c,j)      + &
                       hrv_gresp_xfer_to_litter(p)      * wtcol(p) * leaf_prof(p,j)
                  
                  ! leaf harvest mortality nitrogen fluxes
                  harvest_n_to_litr_met_n(c,j) = harvest_n_to_litr_met_n(c,j) + &
                       hrv_leafn_to_litter(p) * lf_flab(ivt(p)) * wtcol(p) * leaf_prof(p,j)
                  harvest_n_to_litr_cel_n(c,j) = harvest_n_to_litr_cel_n(c,j) + &
                       hrv_leafn_to_litter(p) * lf_fcel(ivt(p)) * wtcol(p) * leaf_prof(p,j)
                  harvest_n_to_litr_lig_n(c,j) = harvest_n_to_litr_lig_n(c,j) + &
                       hrv_leafn_to_litter(p) * lf_flig(ivt(p)) * wtcol(p) * leaf_prof(p,j)
                  
                  ! fine root litter nitrogen fluxes
                  harvest_n_to_litr_met_n(c,j) = harvest_n_to_litr_met_n(c,j) + &
                       hrv_frootn_to_litter(p) * fr_flab(ivt(p)) * wtcol(p) * froot_prof(p,j)
                  harvest_n_to_litr_cel_n(c,j) = harvest_n_to_litr_cel_n(c,j) + &
                       hrv_frootn_to_litter(p) * fr_fcel(ivt(p)) * wtcol(p) * froot_prof(p,j)
                  harvest_n_to_litr_lig_n(c,j) = harvest_n_to_litr_lig_n(c,j) + &
                       hrv_frootn_to_litter(p) * fr_flig(ivt(p)) * wtcol(p) * froot_prof(p,j)
                  
                  ! wood harvest mortality nitrogen fluxes
                  harvest_n_to_cwdn(c,j)  = harvest_n_to_cwdn(c,j)  + &
                       hrv_livestemn_to_litter(p)  * wtcol(p) * stem_prof(p,j)
                  harvest_n_to_cwdn(c,j) = harvest_n_to_cwdn(c,j) + &
                       hrv_livecrootn_to_litter(p) * wtcol(p) * croot_prof(p,j)
                  harvest_n_to_cwdn(c,j) = harvest_n_to_cwdn(c,j) + &
                       hrv_deadcrootn_to_litter(p) * wtcol(p) * croot_prof(p,j)
                  
                  ! retranslocated N pool harvest mortality fluxes
                  harvest_n_to_litr_met_n(c,j) = harvest_n_to_litr_met_n(c,j) + &
                       hrv_retransn_to_litter(p) * wtcol(p) * leaf_prof(p,j)
                  
                  ! storage harvest mortality nitrogen fluxes
                  harvest_n_to_litr_met_n(c,j)      = harvest_n_to_litr_met_n(c,j)      + &
                       hrv_leafn_storage_to_litter(p)      * wtcol(p) * leaf_prof(p,j)
                  harvest_n_to_litr_met_n(c,j)     = harvest_n_to_litr_met_n(c,j)     + &
                       hrv_frootn_storage_to_litter(p)     * wtcol(p) * froot_prof(p,j)
                  harvest_n_to_litr_met_n(c,j)  = harvest_n_to_litr_met_n(c,j)  + &
                       hrv_livestemn_storage_to_litter(p)  * wtcol(p) * stem_prof(p,j)
                  harvest_n_to_litr_met_n(c,j)  = harvest_n_to_litr_met_n(c,j)  + &
                       hrv_deadstemn_storage_to_litter(p)  * wtcol(p) * stem_prof(p,j)
                  harvest_n_to_litr_met_n(c,j) = harvest_n_to_litr_met_n(c,j) + &
                       hrv_livecrootn_storage_to_litter(p) * wtcol(p) * croot_prof(p,j)
                  harvest_n_to_litr_met_n(c,j) = harvest_n_to_litr_met_n(c,j) + &
                       hrv_deadcrootn_storage_to_litter(p) * wtcol(p) * croot_prof(p,j)
                  
                  ! transfer harvest mortality nitrogen fluxes
                  harvest_n_to_litr_met_n(c,j)      = harvest_n_to_litr_met_n(c,j)      + &
                       hrv_leafn_xfer_to_litter(p)      * wtcol(p) * leaf_prof(p,j)
                  harvest_n_to_litr_met_n(c,j)     = harvest_n_to_litr_met_n(c,j)     + &
                       hrv_frootn_xfer_to_litter(p)     * wtcol(p) * froot_prof(p,j)
                  harvest_n_to_litr_met_n(c,j)  = harvest_n_to_litr_met_n(c,j)  + &
                       hrv_livestemn_xfer_to_litter(p)  * wtcol(p) * stem_prof(p,j)
                  harvest_n_to_litr_met_n(c,j)  = harvest_n_to_litr_met_n(c,j)  + &
                       hrv_deadstemn_xfer_to_litter(p)  * wtcol(p) * stem_prof(p,j)
                  harvest_n_to_litr_met_n(c,j) = harvest_n_to_litr_met_n(c,j) + &
                       hrv_livecrootn_xfer_to_litter(p) * wtcol(p) * croot_prof(p,j)
                  harvest_n_to_litr_met_n(c,j) = harvest_n_to_litr_met_n(c,j) + &
                       hrv_deadcrootn_xfer_to_litter(p) * wtcol(p) * croot_prof(p,j)
                  
               end if
            end if
            
         end do
         
      end do
   end do
   
   do pi = 1,maxpatch_pft
      do fc = 1,num_soilc
         c = filter_soilc(fc)
         
         if (pi <=  npfts(c)) then
            p = pfti(c) + pi - 1
            
            if (pactive(p)) then
               
               
               ! wood harvest mortality carbon fluxes to product pools
               chrv_deadstemc_to_prod10c(c)  = chrv_deadstemc_to_prod10c(c)  + &
                    phrv_deadstemc_to_prod10c(p)  * wtcol(p)
               chrv_deadstemc_to_prod100c(c)  = chrv_deadstemc_to_prod100c(c)  + &
                    phrv_deadstemc_to_prod100c(p)  * wtcol(p)
               
               
               ! wood harvest mortality nitrogen fluxes to product pools
               chrv_deadstemn_to_prod10n(c)  = chrv_deadstemn_to_prod10n(c)  + &
                    phrv_deadstemn_to_prod10n(p)  * wtcol(p)
               chrv_deadstemn_to_prod100n(c)  = chrv_deadstemn_to_prod100n(c)  + &
                    phrv_deadstemn_to_prod100n(p)  * wtcol(p)
            end if
         end if
         
      end do
      
   end do
   
    end associate 
  end subroutine CNHarvestPftToColumn

end module pftdynMod
