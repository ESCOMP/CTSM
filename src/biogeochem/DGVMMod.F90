#include <misc.h>
#include <preproc.h>

module DGVMMod

#if (defined DGVM)

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: DGVMMod
!
! !DESCRIPTION:
! Module containing routines to drives the annual portion of lpj
! (called once per year), reset variables related to lpj,
! and initialize/Reset time invariant dgvm variables
!
! !USES:
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use clm_varpar     , only : maxpatch_pft, lsmlon, lsmlat, nlevsoi
  use abortutils     , only : endrun
!
! !PUBLIC TYPES:
  implicit none
  private
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public lpj                ! Drives the annual portion of lpj, called once
                            ! per year
  public lpjreset           ! Resets variables related to LPJ
  public resetTimeConstDGVM ! Initialize/Reset time invariant dgvm variables
  public resetWeightsDGVM   ! Reset DGVM subgrid weights and areas
  public histDGVM           ! Output DGVM history file
!
! !REVISION HISTORY:
! Module created by Mariana Vertenstein
!
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: lpj
!
! !INTERFACE:
  subroutine lpj(lbg, ubg, lbp, ubp, num_natvegp, filter_natvegp, kyr)
!
! !DESCRIPTION:
! Drives the annual portion of lpj, called once per year
!
! !USES:
    use clmtype
    use DGVMReproductionMod , only : Reproduction
    use DGVMTurnoverMod     , only : Turnover
    use DGVMAllocationMod   , only : Allocation
    use DGVMLightMod        , only : Light
    use DGVMMortalityMod    , only : Mortality
    use DGVMFireMod         , only : Fire
    use DGVMEstablishmentMod, only : Establishment
    use DGVMKillMod         , only : Kill
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: lbg, ubg       ! gridcell bounds
    integer, intent(in) :: lbp, ubp       ! pft bounds
    integer, intent(inout) :: num_natvegp ! number of naturally-vegetated
                                          ! pfts in filter
    integer, intent(inout) :: filter_natvegp(ubp-lbp+1) ! filter for
                                          ! naturally-vegetated pfts
    integer, intent(in) :: kyr            ! used in routine climate20 below
!
! !CALLED FROM:
!
! !REVISION HISTORY:
! Author: Sam Levis
!
! !LOCAL VARIABLES:
!
! local pointers to implicit in arguments
!
   integer , pointer :: mxy(:)         ! pft m index (for laixy(i,j,m),etc.)
   integer , pointer :: pgridcell(:)   ! gridcell of corresponding pft
   real(r8), pointer :: fpcgrid(:)     ! foliar projective cover on gridcell (fraction)
   real(r8), pointer :: agdd(:)        ! accumulated growing degree days above 5
   real(r8), pointer :: t_mo_min(:)    ! annual min of t_mo (Kelvin)
!
! local pointers to implicit inout arguments
!
   real(r8), pointer :: tmomin20(:)         ! 20-yr running mean of tmomin
   real(r8), pointer :: agdd20(:)           ! 20-yr running mean of agdd
   real(r8), pointer :: bm_inc(:)           ! biomass increment
   real(r8), pointer :: afmicr(:)           ! annual microbial respiration
   real(r8), pointer :: afirefrac_gcell(:)  ! fraction of gridcell affected by fire
   real(r8), pointer :: acfluxfire_gcell(:) ! gridcell C flux to atmosphere from biomass burning
   real(r8), pointer :: bmfm_gcell(:,:)     ! gridcell biomass
   real(r8), pointer :: afmicr_gcell(:,:)   ! gridcell microbial respiration
!
!EOP
!
! !LOCAL VARIABLES:
    integer  :: g,p                    ! indices
    real(r8) :: afirefrac(lbp:ubp)     ! for history write
    real(r8) :: acfluxfire(lbp:ubp)    ! for history write
!-----------------------------------------------------------------------

    ! Assign local pointers to derived type members (gridcell-level)

    afirefrac_gcell  => clm3%g%gdgvs%afirefrac
    acfluxfire_gcell => clm3%g%gdgvs%acfluxfire
    bmfm_gcell       => clm3%g%gdgvs%bmfm
    afmicr_gcell     => clm3%g%gdgvs%afmicr

    ! Assign local pointers to derived type members (pft-level)

    mxy       => clm3%g%l%c%p%mxy
    pgridcell => clm3%g%l%c%p%gridcell
    fpcgrid   => clm3%g%l%c%p%pdgvs%fpcgrid
    tmomin20  => clm3%g%l%c%p%pdgvs%tmomin20
    t_mo_min  => clm3%g%l%c%p%pdgvs%t_mo_min
    agdd      => clm3%g%l%c%p%pdgvs%agdd
    agdd20    => clm3%g%l%c%p%pdgvs%agdd20
    bm_inc    => clm3%g%l%c%p%pdgvs%bm_inc
    afmicr    => clm3%g%l%c%p%pdgvs%afmicr

    ! *************************************************************************
    ! S. Levis version of LPJ's routine climate20 - 'Returns' tmomin20 and agdd20
    ! for use in routine bioclim, which I have placed in routine Establishment
    ! Instead of 20-yr running mean of coldest monthly temperature,
    ! use 20-yr running mean of minimum 10-day running mean
    ! *************************************************************************

!dir$ concurrent
!cdir nodep
    do p = lbp,ubp
       if (kyr == 2) then
          tmomin20(p) = t_mo_min(p)
          agdd20(p) = agdd(p)
       end if
       tmomin20(p) = (19.0_r8 * tmomin20(p) + t_mo_min(p)) / 20.0_r8
       agdd20(p)   = (19.0_r8 * agdd20(p)   + agdd(p)    ) / 20.0_r8
    end do

    ! Determine grid values of npp and microbial respiration

    bmfm_gcell(lbg:ubg,1:maxpatch_pft) = 0._r8
!dir$ concurrent
!cdir nodep
    do p = lbp,ubp
       if (mxy(p) <= maxpatch_pft) then
          g = pgridcell(p)
          bmfm_gcell(g,mxy(p)) = bm_inc(p)  ! [gC/m2 patch] for output
       end if
       bm_inc(p) = bm_inc(p) * fpcgrid(p)   ! [gC/m2 cell vegetated area]
    end do

    afmicr_gcell(lbg:ubg,1:maxpatch_pft) = 0._r8
!dir$ concurrent
!cdir nodep
    do p = lbp,ubp
       if (mxy(p) <= maxpatch_pft) then
          g = pgridcell(p)
          afmicr_gcell(g,mxy(p)) = afmicr(p) * fpcgrid(p) ![gC/m2 cell veg'd area]
       end if
    end do

    ! Build filter of present natually-vegetated pfts

    call BuildNatVegFilter(lbp, ubp, num_natvegp, filter_natvegp)

    ! Returns updated bm_inc, litterag

    call Reproduction(lbp, ubp, num_natvegp, filter_natvegp)

    ! Returns turnover_ind and updated litterag,bg, l,s,h,rm_ind

    call Turnover(lbp, ubp, num_natvegp, filter_natvegp)

    ! Returns updated litterag, bg, and present

    call Kill(lbp, ubp, num_natvegp, filter_natvegp)

    ! Rebuild filter of present natually-vegetated pfts after Kill()

    call BuildNatVegFilter(lbp, ubp, num_natvegp, filter_natvegp)

    ! Returns lai_ind, lai_inc, updates crownarea, htop, l, s, h, rm_ind, litterag, litterbg

    call Allocation(lbp, ubp, num_natvegp, filter_natvegp)

    ! Returns lm,rm_ind, fpcgrid, nind, litterag,bg via modules
    ! reason for different set up (ie, no external patch loop):
    ! in this routine sub-grid patches (k) communicate at the grid cell level (i,j)

    call Light(lbg, ubg, lbp, ubp, num_natvegp, filter_natvegp)

    ! Obtain updated present, nind, litterag and bg

    call Mortality(lbp, ubp, num_natvegp, filter_natvegp)

    ! Returns updated litterag and nind

    call Fire(lbp, ubp, afirefrac, acfluxfire)

    afirefrac_gcell(lbg:ubg) = 0.0_r8
    acfluxfire_gcell(lbg:ubg) = 0.0_r8
    do p = lbp,ubp
       g = pgridcell(p)
       afirefrac_gcell(g) = afirefrac_gcell(g) + afirefrac(p)*fpcgrid(p)
       acfluxfire_gcell(g) = acfluxfire_gcell(g) + acfluxfire(p)
    end do

    ! Returns updated present, nind, *m_ind, crownarea, fpcgrid, htop,
    ! litter*g, and vegetation type
    ! reason for different set up (ie, no external patch loop):
    ! in this routine sub-grid patches (k) communicate at the grid cell level (i,j)

    call Establishment(lbp, ubp, lbg, ubg)

  end subroutine lpj

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: lpjreset
!
! !INTERFACE:
  subroutine lpjreset(lbg, ubg, lbc, ubc, lbp, ubp, &
                      num_nolakep, filter_nolakep, &
                      caldayp1, declinp1)
!
! !DESCRIPTION:
! Resets variables related to lpj!
!
! !USES:
    use clmtype
    use SurfaceAlbedoMod   , only : SurfaceAlbedo
    use DGVMEcosystemDynMod, only : DGVMEcosystemDyn
!
! !ARGUMENTS:
    implicit none
    integer , intent(in) :: lbg, ubg       ! gridcell bounds
    integer , intent(in) :: lbc, ubc       ! column bounds
    integer , intent(in) :: lbp, ubp       ! pft bounds
    integer , intent(in) :: num_nolakep    ! number of non-lake pfts in filter
    integer , intent(in) :: filter_nolakep(ubp-lbp+1) ! pft filter for non-lake points
    real(r8), intent(in) :: caldayp1       ! calendar day at Greenwich (1.00, ..., 365.99) for nstep+1
    real(r8), intent(in) :: declinp1       ! declination angle for next time step
!
! !CALLED FROM:
! subroutine driver() in driver.F90
!
! !REVISION HISTORY:
! Author: Sam Levis
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: p         ! indices
!-----------------------------------------------------------------------

    ! Reset a few variables here at the very end of the year
    ! First determine necessary per processor subgrid bounds

!dir$ concurrent
!cdir nodep
    do p = lbp,ubp
       clm3%g%l%c%p%pdgvs%annpsn(p)     = 0._r8
       clm3%g%l%c%p%pdgvs%annpsnpot(p)  = 0._r8
       clm3%g%l%c%p%pdgvs%bm_inc(p)     = 0._r8
       clm3%g%l%c%p%pdgvs%afmicr(p)     = 0._r8
       clm3%g%l%c%p%pdgvs%firelength(p) = 0._r8
       clm3%g%l%c%p%pdgvs%agddtw(p)     = 0._r8
       clm3%g%l%c%p%pdgvs%agdd(p)       = 0._r8
       clm3%g%l%c%p%pdgvs%t10min(p)     = 1.0e+36_r8
       clm3%g%l%c%p%pdgvs%t_mo_min(p)   = 1.0e+36_r8
    end do

    ! Call DGVMEcosystemDyn and SurfaceAlbedo because need information
    ! for first timestep of next year.

    call DGVMEcosystemDyn(lbp, ubp, num_nolakep, filter_nolakep, &
         doalb=.false., endofyr=.true.)

    call SurfaceAlbedo(lbg, ubg, lbc, ubc, lbp, ubp, &
         caldayp1, declinp1)

  end subroutine lpjreset

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: histDGVM
!
! !INTERFACE:
  subroutine histDGVM()
!
! !DESCRIPTION:
! Create DGVM history file
!
! !USES:
    use clmtype
    use ncdio
    use decompMod    , only : get_proc_bounds, get_proc_global
    use clm_varpar   , only : lsmlon, lsmlat, maxpatch_pft
    use domainMod    , only : ldomain
    use decompMod    , only : ldecomp
    use clm_varctl   , only : caseid, ctitle, finidat, fsurdat, fpftcon, &
                              frivinp_rtm, archive_dir, mss_wpass, mss_irt
    use clm_varcon   , only : spval
    use clm_time_manager , only : get_ref_date, get_nstep, get_curr_date, &
                              get_curr_time
    use fileutils    , only : set_filename, putfil, get_filename
    use shr_sys_mod  , only : shr_sys_getenv
    use spmdMod      , only : masterproc
    use spmdGathScatMod, only : gather_data_to_master
    use shr_const_mod, only : SHR_CONST_CDAY
!
! !ARGUMENTS:
    implicit none
!
! !CALLED FROM:
!
! !REVISION HISTORY:
! Author: Sam Levis
!
! !LOCAL VARIABLES:
!
! local pointers to implicit in arguments
!
   logical , pointer :: ifspecial(:)        ! true=>landunit is not vegetated (landunit-level)
   integer , pointer :: pgridcell(:)        ! gridcell index of corresponding pft (pft-level)
   integer , pointer :: plandunit(:)        ! landunit index of corresponding pft (pft-level)
   integer , pointer :: ivt(:)              ! pft vegetation (pft-level)
   integer , pointer :: mxy(:)              ! pft m index (for laixy(i,j,m),etc.)
   real(r8), pointer :: fpcgrid(:)          ! foliar projective cover on gridcell (fraction)
   real(r8), pointer :: lm_ind(:)           ! individual leaf mass
   real(r8), pointer :: sm_ind(:)           ! individual sapwood mass
   real(r8), pointer :: hm_ind(:)           ! individual heartwood mass
   real(r8), pointer :: rm_ind(:)           ! individual root mass
   real(r8), pointer :: nind(:)             ! number of individuals (#/m**2)
   real(r8), pointer :: afirefrac_gcell(:)  ! fraction of gridcell affected by fire
   real(r8), pointer :: acfluxfire_gcell(:) ! gridcell C flux to atmosphere from biomass burning
   real(r8), pointer :: bmfm_gcell(:,:)     ! gridcell biomass
   real(r8), pointer :: afmicr_gcell(:,:)   ! gridcell microbial respiration
   real(r8), pointer :: data(:)             ! temporary global array
!
!EOP
!
! !LOCAL VARIABLES:
    character(len=256) :: dgvm_fn      ! dgvm history filename
    integer :: ncid                    ! dgvm netcdf file id
    integer :: omode                   ! returned mode from netCDF call
    integer :: ncprec                  ! output precision
    integer :: g,p,l                   ! indices
    integer :: begp, endp              ! per-proc beginning and ending pft indices
    integer :: begc, endc              ! per-proc beginning and ending column indices
    integer :: begl, endl              ! per-proc beginning and ending landunit indices
    integer :: begg, endg              ! per-proc gridcell ending gridcell indices
    integer :: numg,numl,numc,nump     ! global glcp cells
    integer :: ier                     ! error status
    integer :: n,m                     ! index
    integer :: mdcur, mscur, mcdate    ! outputs from get_curr_time
    integer :: yr,mon,day,mcsec        ! outputs from get_curr_date
    integer :: hours,minutes,secs      ! hours,minutes,seconds of hh:mm:ss
    integer :: nstep                   ! time step
    integer :: nbsec                   ! seconds components of a date
    integer :: dimid                   ! dimension, variable id
    real(r8):: time                    ! current time
    real(r8),pointer :: lonvar(:)      ! only used for full grid
    real(r8),pointer :: latvar(:)      ! only used for full grid
    character(len=256) :: str          ! temporary string
    character(len=  8) :: curdate      ! current date
    character(len=  8) :: curtime      ! current time
    character(len= 10) :: basedate     ! base date (yyyymmdd)
    character(len=  8) :: basesec      ! base seconds
    character(len=256) :: rem_dir      ! remote (archive) directory
    character(len=256) :: rem_fn       ! remote (archive) filename
    real(r8), pointer :: rbuf2dg(:,:)  ! temporary
    integer , pointer :: ibuf2dg(:,:)  ! temporary
    character(len=32) :: subname='histDGVM'
!-----------------------------------------------------------------------

    ! Assign local pointers to derived type members (gridcell-level)

    afirefrac_gcell  => clm3%g%gdgvs%afirefrac
    acfluxfire_gcell => clm3%g%gdgvs%acfluxfire
    bmfm_gcell       => clm3%g%gdgvs%bmfm
    afmicr_gcell     => clm3%g%gdgvs%afmicr

    ! Assign local pointers to derived type members (landunit-level)

    ifspecial  => clm3%g%l%ifspecial

    ! Assign local pointers to derived subtypes components (pft-level)

    ivt       => clm3%g%l%c%p%itype
    mxy       => clm3%g%l%c%p%mxy
    pgridcell => clm3%g%l%c%p%gridcell
    plandunit => clm3%g%l%c%p%landunit
    fpcgrid   => clm3%g%l%c%p%pdgvs%fpcgrid
    lm_ind    => clm3%g%l%c%p%pdgvs%lm_ind
    rm_ind    => clm3%g%l%c%p%pdgvs%rm_ind
    sm_ind    => clm3%g%l%c%p%pdgvs%sm_ind
    hm_ind    => clm3%g%l%c%p%pdgvs%hm_ind
    nind      => clm3%g%l%c%p%pdgvs%nind

    ! Determine subgrid bounds for this processor and allocate dynamic memory

    call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)
    call get_proc_global(numg, numl, numc, nump)

    allocate(rbuf2dg(begg:endg,maxpatch_pft), ibuf2dg(begg:endg,maxpatch_pft), stat=ier)
    if (ier /= 0) call endrun('histDGVM: allocation error for rbuf2dg, ibuf2dg')

    ! Set output precision

    ncprec = nf_double

    if (masterproc) then

       ! -----------------------------------------------------------------------
       ! Create new netCDF file. File will be in define mode
       ! Set fill mode to "no fill" to optimize performance
       ! -----------------------------------------------------------------------

       dgvm_fn = set_dgvm_filename()
       call check_ret(nf_create(dgvm_fn, nf_clobber, ncid), subname)
       call check_ret(nf_set_fill(ncid, nf_nofill, omode), subname)

       ! -----------------------------------------------------------------------
       ! Create global attributes.
       ! -----------------------------------------------------------------------

       str = 'CF1.0'
       call check_ret(&
            nf_put_att_text (ncid, NF_GLOBAL, 'conventions', len_trim(str), trim(str)), subname)

       call getdatetime(curdate, curtime)
       str = 'created on ' // curdate // ' ' // curtime
       call check_ret(&
            nf_put_att_text(ncid, NF_GLOBAL,'history', len_trim(str), trim(str)), subname)

       call shr_sys_getenv('LOGNAME', str, ier)
       if (ier /= 0) call endrun('error: LOGNAME environment variable not defined')
       
       call check_ret(&
            nf_put_att_text (ncid, NF_GLOBAL, 'logname',len_trim(str), trim(str)), subname)
       
       call shr_sys_getenv('HOST', str, ier)
       call check_ret(&
            nf_put_att_text (ncid, NF_GLOBAL, 'host', len_trim(str), trim(str)), subname)
       
       str = 'Community Land Model: CLM3'
       call check_ret(&
            nf_put_att_text (ncid, NF_GLOBAL, 'source', len_trim(str), trim(str)), subname)
       
       str = '$Name$'
       call check_ret(&
            nf_put_att_text (ncid, NF_GLOBAL, 'version', len_trim(str), trim(str)), subname)
       
       str = '$Id$'
       call check_ret(&
            nf_put_att_text (ncid, NF_GLOBAL, 'revision_id', len_trim(str), trim(str)), subname)

       str = ctitle
       call check_ret(&
            nf_put_att_text (ncid, NF_GLOBAL, 'case_title', len_trim(str), trim(str)), subname)

       str = caseid
       call check_ret(&
            nf_put_att_text (ncid, NF_GLOBAL, 'case_id', len_trim(str), trim(str)), subname)
       
       str = get_filename(fsurdat)
       call check_ret(&
            nf_put_att_text(ncid, NF_GLOBAL, 'Surface_dataset', len_trim(str), trim(str)), subname)
       
       str = 'arbitrary initialization'
       if (finidat /= ' ') str = get_filename(finidat)
       call check_ret(&
            nf_put_att_text(ncid, NF_GLOBAL, 'Initial_conditions_dataset', len_trim(str), trim(str)), subname)
       
       str = get_filename(fpftcon)
       call check_ret(&
            nf_put_att_text(ncid, NF_GLOBAL, 'PFT_physiological_constants_dataset', len_trim(str), trim(str)), subname)
       
       if (frivinp_rtm /= ' ') then
          str = get_filename(frivinp_rtm)
          call check_ret(&
               nf_put_att_text(ncid, NF_GLOBAL, 'RTM_input_datset', len_trim(str), trim(str)), subname)
       end if

       ! -----------------------------------------------------------------------
       ! Define dimensions.
       ! -----------------------------------------------------------------------
       
       call check_ret(nf_def_dim (ncid, 'lon', lsmlon, dimid), subname)
       call check_ret(nf_def_dim (ncid, 'lat', lsmlat, dimid), subname)
       call check_ret(nf_def_dim (ncid, 'pft', maxpatch_pft, dimid), subname)
       call check_ret(nf_def_dim (ncid, 'time', nf_unlimited, dimid), subname)
       call check_ret(nf_def_dim (ncid, 'string_length', 80, dimid), subname)
       
       ! -----------------------------------------------------------------------
       ! Define variables
       ! -----------------------------------------------------------------------
       
       ! Define coordinate variables (including time)
       
       call ncd_defvar(ncid=ncid, varname='lon', xtype=ncprec, dim1name='lon',&
            long_name='coordinate longitude', units='degrees_east')
          
       call ncd_defvar(ncid=ncid, varname='lat', xtype=ncprec, dim1name='lat',&
            long_name='coordinate latitude', units='degrees_north')

       call get_curr_time(mdcur, mscur)
       call get_ref_date(yr, mon, day, nbsec)
       hours   = nbsec / 3600
       minutes = (nbsec - hours*3600) / 60
       secs    = (nbsec - hours*3600 - minutes*60)
       write(basedate,80) yr,mon,day
80     format(i4.4,'-',i2.2,'-',i2.2)
       write(basesec ,90) hours, minutes, secs
90     format(i2.2,':',i2.2,':',i2.2)
       str = 'days since ' // basedate // " " // basesec
       time = mdcur + mscur/SHR_CONST_CDAY
       
       call ncd_defvar(ncid=ncid, varname='time', xtype=nf_double, dim1name='time', &
            long_name='time', units=str)
       
       call ncd_defvar(ncid=ncid, varname='edgen', xtype=ncprec, &
            long_name='northern edge of surface grid', units='degrees_north')
          
       call ncd_defvar(ncid=ncid, varname='edgee', xtype=ncprec, &
            long_name='eastern edge of surface grid', units='degrees_east')
          
       call ncd_defvar(ncid=ncid, varname='edges', xtype=ncprec, &
            long_name='southern edge of surface grid', units='degrees_north')
          
       call ncd_defvar(ncid=ncid, varname='edgew', xtype=ncprec, &
            long_name='western edge of surface grid', units='degrees_east')

       ! Define surface grid (coordinate variables, latitude, longitude, surface type).
       
       call ncd_defvar(ncid=ncid, varname='longxy', xtype=ncprec, dim1name='lon', dim2name='lat', &
            long_name='longitude', units='degrees_east')
       
       call ncd_defvar(ncid=ncid, varname='latixy', xtype=ncprec, dim1name='lon', dim2name='lat', &
            long_name='latitude', units='degrees_north')
       
       call ncd_defvar(ncid=ncid, varname='landmask', xtype=nf_int, dim1name='lon', dim2name='lat', &
            long_name='land/ocean mask (0.=ocean and 1.=land)')
       
       ! Define time information
       
       call ncd_defvar(ncid=ncid, varname='mcdate', xtype=nf_int, dim1name='time',&
            long_name='current date (YYYYMMDD)')
       
       call ncd_defvar(ncid=ncid, varname='mcsec', xtype=nf_int, dim1name='time',&
            long_name='current seconds of current date', units='s')
       
       call ncd_defvar(ncid=ncid, varname='mdcur', xtype=nf_int, dim1name='time',&
            long_name='current day (from base day)')
       
       call ncd_defvar(ncid=ncid, varname='mscur', xtype=nf_int, dim1name='time',&
            long_name='current seconds of current day', units='s')

       call ncd_defvar(ncid=ncid, varname='nstep', xtype=nf_int, dim1name='time',&
            long_name='time step', units='s')

       ! Define time dependent variables
       
       call ncd_defvar(ncid=ncid, varname='BURN', xtype=ncprec, &
            dim1name='lon', dim2name='lat', dim3name='time', &
            long_name='fraction of vegetated area burned', &
            missing_value=spval, fill_value=spval)

       call ncd_defvar(ncid=ncid, varname='CFLUXFIRE', xtype=ncprec, &
            dim1name='lon', dim2name='lat', dim3name='time', &
            long_name='carbon flux to the atmosphere due to fire', units='grams C/m2 of vegetated area', &
            missing_value=spval, fill_value=spval)

       call ncd_defvar(ncid=ncid, varname='NPP', xtype=ncprec, &
            dim1name='lon', dim2name='lat', dim3name='pft', dim4name='time', &
            long_name='annual net primary production', units='grams C/m2 of patch', &
            missing_value=spval, fill_value=spval)

       call ncd_defvar(ncid=ncid, varname='Rh', xtype=ncprec, &
            dim1name='lon', dim2name='lat', dim3name='pft', dim4name='time', &
            long_name='annual heterotrophic respiration', units='grams C/m2 of vegetated area', &
            missing_value=spval, fill_value=spval)

       call ncd_defvar(ncid=ncid, varname='PFT', xtype=nf_int, &
            dim1name='lon', dim2name='lat', dim3name='pft', dim4name='time', &
            long_name='plant functional type', &
            imissing_value=9999, ifill_value=9999)

       call ncd_defvar(ncid=ncid, varname='FPCGRID', xtype=ncprec, &
            dim1name='lon', dim2name='lat', dim3name='pft', dim4name='time', &
            long_name='plant functional type cover', units='fraction of vegetated area', &
            missing_value=spval, fill_value=spval)

       call ncd_defvar(ncid=ncid, varname='LCIND', xtype=ncprec, &
            dim1name='lon', dim2name='lat', dim3name='pft', dim4name='time', &
            long_name='leaf carbon per individual', units='grams C', &
            missing_value=spval, fill_value=spval)

       call ncd_defvar(ncid=ncid, varname='RCIND', xtype=ncprec, &
            dim1name='lon', dim2name='lat', dim3name='pft', dim4name='time', &
            long_name='root carbon per individual', units='grams C', &
            missing_value=spval, fill_value=spval)

       call ncd_defvar(ncid=ncid, varname='SCIND', xtype=ncprec, &
            dim1name='lon', dim2name='lat', dim3name='pft', dim4name='time', &
            long_name='stem carbon per individual', units='grams C', &
            missing_value=spval, fill_value=spval)

       call ncd_defvar(ncid=ncid, varname='HCIND', xtype=ncprec, &
            dim1name='lon', dim2name='lat', dim3name='pft', dim4name='time', &
            long_name='heartwood carbon per individual', units='grams C', &
            missing_value=spval, fill_value=spval)

       call ncd_defvar(ncid=ncid, varname='NIND', xtype=ncprec, &
            dim1name='lon', dim2name='lat', dim3name='pft', dim4name='time', &
            long_name='number of individuals', units='individuals/m2 vegetated land', &
            missing_value=spval, fill_value=spval)
       
       call check_ret(nf_enddef(ncid), subname)

    end if   ! End of if-masterproc block

    ! -----------------------------------------------------------------------
    ! Write variables
    ! -----------------------------------------------------------------------

    call ncd_ioglobal(varname='edgen', data=ldomain%edges(1), ncid=ncid, flag='write')
    call ncd_ioglobal(varname='edgee', data=ldomain%edges(2), ncid=ncid, flag='write')
    call ncd_ioglobal(varname='edges', data=ldomain%edges(3), ncid=ncid, flag='write')
    call ncd_ioglobal(varname='edgew', data=ldomain%edges(4), ncid=ncid, flag='write')

    ! Write surface grid (coordinate variables, latitude, longitude, surface type).

    allocate(lonvar(lsmlon),latvar(lsmlat),data(numg))

    call gather_data_to_master (ldomain%lonc, data, clmlevel='gridcell')
    lonvar = spval
    do n = 1,lsmlon
    do m = 1,lsmlat
       g = ldecomp%glo2gdc((m-1)*lsmlon + n)
       if (g > 0 .and. g < lsmlon*lsmlat) lonvar(n) = data(g)
    enddo
    enddo

    call gather_data_to_master (ldomain%latc, data, clmlevel='gridcell')
    latvar = spval
    do n = 1,lsmlon
    do m = 1,lsmlat
       g = ldecomp%glo2gdc((m-1)*lsmlon + n)
       if (g > 0 .and. g < lsmlon*lsmlat) latvar(m) = data(g)
    enddo
    enddo

    if (masterproc) then
       call ncd_ioglobal(varname='lon', data=lonvar, ncid=ncid, flag='write')
       call ncd_ioglobal(varname='lat', data=latvar, ncid=ncid, flag='write')
    endif

    call ncd_iolocal(varname='longxy'  , data=ldomain%lonc, ncid=ncid, &
         flag='write', dim1name='gridcell', &
         nlonxy=ldomain%ni, nlatxy=ldomain%nj)
    call ncd_iolocal(varname='latixy'  , data=ldomain%latc, ncid=ncid, &
         flag='write', dim1name='gridcell', &
         nlonxy=ldomain%ni, nlatxy=ldomain%nj)
    call ncd_iolocal(varname='landmask', data=ldomain%mask, ncid=ncid, &
         flag='write', dim1name='gridcell', &
         nlonxy=ldomain%ni, nlatxy=ldomain%nj)

    deallocate(lonvar,latvar,data)

    ! Write current date, current seconds, current day, current nstep

    call get_curr_date(yr, mon, day, mcsec)
    mcdate = yr*10000 + mon*100 + day
    nstep = get_nstep()

    call ncd_ioglobal(varname='mcdate', data=mcdate, nt=1, ncid=ncid, flag='write')
    call ncd_ioglobal(varname='mcsec' , data=mcsec , nt=1, ncid=ncid, flag='write')
    call ncd_ioglobal(varname='mdcur' , data=mdcur , nt=1, ncid=ncid, flag='write')
    call ncd_ioglobal(varname='mscur' , data=mcsec , nt=1, ncid=ncid, flag='write')
    call ncd_ioglobal(varname='nstep' , data=nstep , nt=1, ncid=ncid, flag='write')
    call ncd_ioglobal(varname='time'  , data=time  , nt=1, ncid=ncid, flag='write')

    ! Write time dependent variables to DGVM history file

    call ncd_iolocal(varname='BURN', dim1name='gridcell', data=afirefrac_gcell, &
         nlonxy=lsmlon, nlatxy=lsmlat, nt=1, ncid=ncid, flag='write')

    call ncd_iolocal(varname='CFLUXFIRE', dim1name='gridcell', data=acfluxfire_gcell, &
         nlonxy=lsmlon, nlatxy=lsmlat, nt=1, ncid=ncid, flag='write')

    call ncd_iolocal(varname='NPP', dim1name='gridcell', dim2name='pft', data=bmfm_gcell, &
         nlonxy=lsmlon, nlatxy=lsmlat, nt=1, ncid=ncid, flag='write')

    call ncd_iolocal(varname='Rh', dim1name='gridcell', dim2name='pft', data=afmicr_gcell, &
         nlonxy=lsmlon, nlatxy=lsmlat, nt=1, ncid=ncid, flag='write')

    ! Note that checking on if-not ifspecial status below guarantees that the m index will
    ! always lie between 1 and maxpatch_pft

    ibuf2dg(:,:) = 0
!dir$ concurrent
!cdir nodep
    do p = begp,endp
       g = pgridcell(p)
       l = plandunit(p)
       if (.not. ifspecial(l)) ibuf2dg(g,mxy(p)) = ivt(p)
    end do
    call ncd_iolocal(varname='PFT', dim1name='gridcell', dim2name='pft', data=ibuf2dg, &
         nlonxy=lsmlon, nlatxy=lsmlat, nt=1, ncid=ncid, flag='write')

    rbuf2dg(:,:) = 0._r8
!dir$ concurrent
!cdir nodep
    do p = begp,endp
       g = pgridcell(p)
       l = plandunit(p)
       if (.not. ifspecial(l)) rbuf2dg(g,mxy(p)) = fpcgrid(p)*100._r8
    end do
    call ncd_iolocal(varname='FPCGRID', dim1name='gridcell', dim2name='pft', data=rbuf2dg, &
         nlonxy=lsmlon, nlatxy=lsmlat, nt=1, ncid=ncid, flag='write')

    rbuf2dg(:,:) = 0._r8
!dir$ concurrent
!cdir nodep
    do p = begp,endp
       g = pgridcell(p)
       l = plandunit(p)
       if (.not. ifspecial(l)) rbuf2dg(g,mxy(p)) = lm_ind(p)
    end do
    call ncd_iolocal(varname='LCIND', dim1name='gridcell', dim2name='pft', data=rbuf2dg, &
         nlonxy=lsmlon, nlatxy=lsmlat, nt=1, ncid=ncid, flag='write')

    rbuf2dg(:,:) = 0._r8
!dir$ concurrent
!cdir nodep
    do p = begp,endp
       g = pgridcell(p)
       l = plandunit(p)
       if (.not. ifspecial(l)) rbuf2dg(g,mxy(p)) = rm_ind(p)
    end do
    call ncd_iolocal(varname='RCIND', dim1name='gridcell', dim2name='pft', data=rbuf2dg, &
         nlonxy=lsmlon, nlatxy=lsmlat, nt=1, ncid=ncid, flag='write')

    rbuf2dg(:,:) = 0._r8
!dir$ concurrent
!cdir nodep
    do p = begp,endp
       g = pgridcell(p)
       l = plandunit(p)
       if (.not. ifspecial(l)) rbuf2dg(g,mxy(p)) = sm_ind(p)
    end do
    call ncd_iolocal(varname='SCIND', dim1name='gridcell', dim2name='pft', data=rbuf2dg, &
         nlonxy=lsmlon, nlatxy=lsmlat, nt=1, ncid=ncid, flag='write')

    rbuf2dg(:,:) = 0._r8
!dir$ concurrent
!cdir nodep
    do p = begp,endp
       g = pgridcell(p)
       l = plandunit(p)
       if (.not. ifspecial(l)) rbuf2dg(g,mxy(p)) = hm_ind(p)
    end do
    call ncd_iolocal(varname='HCIND', dim1name='gridcell', dim2name='pft', data=rbuf2dg, &
         nlonxy=lsmlon, nlatxy=lsmlat, nt=1, ncid=ncid, flag='write')

    rbuf2dg(:,:) = 0._r8
!dir$ concurrent
!cdir nodep
    do p = begp,endp
       g = pgridcell(p)
       l = plandunit(p)
       if (.not. ifspecial(l)) rbuf2dg(g,mxy(p)) = nind(p)
    end do
    call ncd_iolocal(varname='NIND',dim1name='gridcell', dim2name='pft', data=rbuf2dg, &
         nlonxy=lsmlon, nlatxy=lsmlat, nt=1, ncid=ncid, flag='write')

    ! Deallocate dynamic memory

    deallocate(rbuf2dg, ibuf2dg)

    !------------------------------------------------------------------
    ! Close and archive netcdf DGVM history file
    !------------------------------------------------------------------

    if (masterproc) then
       call check_ret(nf_close(ncid), subname)
       write(6,*)'(histDGVM): Finished writing clm DGVM history dataset ',&
            trim(dgvm_fn), 'at nstep = ',get_nstep()
       if (mss_irt > 0) then
          rem_dir = trim(archive_dir) // '/hist/'
          rem_fn = set_filename(rem_dir, dgvm_fn)
          call putfil (dgvm_fn, rem_fn, mss_wpass, mss_irt, .true.)
       end if
    end if

  end subroutine histDGVM

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: resetTimeConstDGVM
!
! !INTERFACE:
  subroutine resetTimeConstDGVM(lbp, ubp)
!
! !DESCRIPTION:
! Initialize/reset time invariant DGVM variables
!
! !USES:
    use clmtype
    use pftvarcon , only : roota_par, rootb_par, noveg
    use clm_varcon, only : spval
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: lbp, ubp       ! pft bounds
!
! !CALLED FROM:
! lpjreset1() in this module
! initialize() in initializeMod.F90
! iniTimeVar() in iniTimeVar.F90
!
! !REVISION HISTORY:
! Author: Gordon Bonan
!
! !LOCAL VARIABLES:
!
! local pointers to implicit in arguments
!
   real(r8), pointer :: zi(:,:)        ! interface level below a "z" level (m) (-nlevsno+0:nlevsoi)
   integer , pointer :: ivt(:)         ! pft vegetation
   integer , pointer :: pcolumn(:)     ! column of corresponding pft
   real(r8), pointer :: rootfr(:,:)    ! fraction of roots in each soil layer  (nlevsoi)
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: p,c,j        ! indices
!-----------------------------------------------------------------------

    ! Assign local pointers to derived subtypes components (column-level)

    zi      => clm3%g%l%c%cps%zi

    ! Assign local pointers to derived subtypes components (pft-level)

    ivt     => clm3%g%l%c%p%itype
    pcolumn => clm3%g%l%c%p%column
    rootfr  => clm3%g%l%c%p%pps%rootfr

    ! Initialize root fraction (computing from surface, d is depth in meter):
    ! Y = 1 -1/2 (exp(-ad)+exp(-bd) under the constraint that
    ! Y(d =0.1m) = 1-beta^(10 cm) and Y(d=d_obs)=0.99 with
    ! beta & d_obs given in Zeng et al. (1998).

!dir$ concurrent
!cdir nodep
    do p = lbp,ubp
       c = pcolumn(p)
       if (ivt(p) /= noveg) then
          do j = 1, nlevsoi-1
             rootfr(p,j) = .5_r8*( exp(-roota_par(ivt(p)) * zi(c,j-1))  &
                              + exp(-rootb_par(ivt(p)) * zi(c,j-1))  &
                              - exp(-roota_par(ivt(p)) * zi(c,j  ))  &
                              - exp(-rootb_par(ivt(p)) * zi(c,j  )) )
          end do
          rootfr(p,nlevsoi) = .5_r8*( exp(-roota_par(ivt(p)) * zi(c,nlevsoi-1))  &
                                 + exp(-rootb_par(ivt(p)) * zi(c,nlevsoi-1)) )
       else
          rootfr(p,1:nlevsoi) = spval
       end if
    end do

  end subroutine resetTimeConstDGVM

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: resetWeightsDGVM
!
! !INTERFACE:
  subroutine resetWeightsDGVM(lbg, ubg, lbc, ubc, lbp, ubp)
!
! !DESCRIPTION:
! Determine new subgrid weights and areas
! In CLM3 with satellite data, the number of veg pfts is determined once
! and is less than maxpatch_pft (4) in some cells.
! In LSM with LPJ, the number of veg patches could be dynamic. Until we
! implement it as such, we will make all grid cells have 10 veg patches.
!------------------------------------------------------------------

!
! !USES:
    use clmtype
    use decompMod, only : ldecomp
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: lbg, ubg       ! gridcell bounds
    integer, intent(in) :: lbc, ubc       ! column bounds
    integer, intent(in) :: lbp, ubp       ! pft bounds
!
! !CALLED FROM:
!  subroutine restart_dgvm in module DGVMRestMod: if the restart file is read
!  subroutine inicrd in module inicFileMod: if the initial file is read
!  subroutine mkarbinit in module iniTimeVar
!
! !REVISION HISTORY:
! Author: Gordon Bonan
!
! !LOCAL VARIABLES:
!
! local pointers to implicit in arguments
!
   integer , pointer :: ixy(:)            ! gridcell lon index (gridcell level)
   integer , pointer :: jxy(:)            ! gridcell lat index (gridcell level)
   integer , pointer :: ltype(:)          ! landunit type
   logical , pointer :: ifspecial(:)      ! true=>landunit is not vegetated
   real(r8), pointer :: lwtgcell(:)       ! weight (relative to gridcell) for this landunit
   real(r8), pointer :: fpcgrid(:)        ! weight of pft relative to vegetated landunit
!
! local pointers to implicit out arguments
!
   integer , pointer :: pcolumn(:)        ! index into column for each pft
   integer , pointer :: plandunit(:)      ! index into landunit for each pft
   integer , pointer :: pgridcell(:)      ! index into gridcell for each pft
   real(r8), pointer :: pwtcol(:)         ! weight (relative to column) for this pft (0-1)
   real(r8), pointer :: pwtlunit(:)       ! weight (relative to landunit) for this pft (0-1)
   real(r8), pointer :: pwtgcell(:)       ! weight (relative to gridcell) for this pft (0-1)
!
!EOP
!
! !LOCAL VARIABLES:
    integer  :: g,p,c,l             ! indices
    integer  :: fn,filterg(lbg:ubg) ! local gridcell filter for error check
    real(r8) :: sumwt(lbg:ubg)      ! consistency check
!-----------------------------------------------------------------------

    ! Assign local pointers to derived subtypes components (gridcell-level)

    ixy        => ldecomp%gdc2i
    jxy        => ldecomp%gdc2j

    ! Assign local pointers to derived subtypes components (landunit-level)

    ltype      => clm3%g%l%itype
    ifspecial  => clm3%g%l%ifspecial
    lwtgcell   => clm3%g%l%wtgcell

    ! Assign local pointers to derived subtypes components (pft-level)

    pgridcell  => clm3%g%l%c%p%gridcell
    plandunit  => clm3%g%l%c%p%landunit
    pcolumn    => clm3%g%l%c%p%column
    pwtcol     => clm3%g%l%c%p%wtcol
    pwtlunit   => clm3%g%l%c%p%wtlunit
    pwtgcell   => clm3%g%l%c%p%wtgcell
    fpcgrid    => clm3%g%l%c%p%pdgvs%fpcgrid

    ! Determine new pft properties

!dir$ concurrent
!cdir nodep
    do p = lbp,ubp
       g = pgridcell(p)
       l = plandunit(p)
       c = pcolumn(p)
       if (.not. ifspecial(l)) then

          ! Determine pft weight relative to column and relative to landunit
          ! One column per landunit - column and landunit areas are equal
          ! Weight relative to column and weight relative to landunit are equal
          pwtcol(p) = fpcgrid(p)
          pwtlunit(p) = fpcgrid(p)

          ! Determine new pft weight relative to grid cell
          pwtgcell(p) = pwtlunit(p) * lwtgcell(l)

       end if
    end do

    ! Consistency check - add up all the pft weights for a given gridcell
    ! and make sure they are not greater than one.

    sumwt(:) = 0._r8
    do p = lbp,ubp
       g = pgridcell(p)
       sumwt(g) = sumwt(g) + pwtgcell(p)
    end do
    fn = 0
    do g = lbg,ubg
       if (abs(sumwt(g) - 1.0_r8) > 1.0e-6_r8) then
          fn = fn + 1
          filterg(fn) = g
       end if
    end do
    if (fn > 0) then
       g = filterg(1)
       write(6,*) 'resetWeightsDGVM: sumwt of pfts for grid cell ',&
         'i,j = ',ixy(g),jxy(g),' not equal to 1'
       write(6,*) 'sum of pft weights for gridcell =',sumwt(g)
       call endrun
    end if

  end subroutine resetWeightsDGVM

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: set_dgvm_filename
!
! !INTERFACE:
  character(len=256) function set_dgvm_filename ()
!
! !DESCRIPTION:
! Determine initial dataset filenames
!
! !USES:
    use clm_varctl  , only : caseid
    use clm_time_manager, only : get_curr_date
!
! !ARGUMENTS:
    implicit none
!
! !CALLED FROM:
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
    character(len=256) :: cdate       !date char string
    integer :: day                    !day (1 -> 31)
    integer :: mon                    !month (1 -> 12)
    integer :: yr                     !year (0 -> ...)
    integer :: sec                    !seconds into current day
!-----------------------------------------------------------------------

    call get_curr_date (yr, mon, day, sec)
    write(cdate,'(i4.4,"-",i2.2,"-",i2.2,"-",i5.5)') yr,mon,day,sec
    set_dgvm_filename = "./"//trim(caseid)//".clm2.hv."//trim(cdate)//".nc"

  end function set_dgvm_filename

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: BuildNatVegFilter
!
! !INTERFACE:
  subroutine BuildNatVegFilter(lbp, ubp, num_natvegp, filter_natvegp)
!
! !DESCRIPTION:
! Reconstruct a filter of naturally-vegetated PFTs for use in DGVM
!
! !USES:
    use clmtype
    use pftvarcon , only : crop
!
! !ARGUMENTS:
    implicit none
    integer, intent(in)  :: lbp, ubp                   ! pft bounds
    integer, intent(out) :: num_natvegp                ! number of pfts in naturally-vegetated filter
    integer, intent(out) :: filter_natvegp(ubp-lbp+1)  ! pft filter for naturally-vegetated points
!
! !CALLED FROM:
! subroutine lpj in this module
!
! !REVISION HISTORY:
! Author: Forrest Hoffman
!
! !LOCAL VARIABLES:
! local pointers to implicit in arguments
    integer , pointer :: ivt(:)         ! pft vegetation (pft level)
    logical , pointer :: present(:)     ! whether this pft present in patch
!EOP
!
! !LOCAL VARIABLES:
    integer :: p
!-----------------------------------------------------------------------

    ! Assign local pointers to derived type members (pft-level)
    ivt       => clm3%g%l%c%p%itype
    present   => clm3%g%l%c%p%pdgvs%present

    num_natvegp = 0
    do p = lbp,ubp
       if (ivt(p) > 0 .and. present(p) .and. crop(ivt(p)) == 0._r8) then
          num_natvegp = num_natvegp + 1
          filter_natvegp(num_natvegp) = p
       end if
    end do

  end subroutine BuildNatVegFilter

#endif

end module DGVMMod
