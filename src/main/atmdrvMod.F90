#include <misc.h>
#include <preproc.h>

module atmdrvMod

#if (defined OFFLINE)

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: atmdrvMod
!
! !DESCRIPTION:
! Read and generate atmospheric grid data at model resolution
!
! !USES:
  use shr_kind_mod , only : r8 => shr_kind_r8
  use shr_const_mod, only : SHR_CONST_TKFRZ, SHR_CONST_PSTD
  use abortutils   , only : endrun
  use spmdMod      , only : masterproc, mpicom, comp_id, MPI_REAL8, MPI_INTEGER, iam
  use clm_mct_mod
  use decompMod    , only : gsMap_atm_gdc2glo, perm_atm_gdc2glo
  use perf_mod
!
! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: atmdrv_init  ! read atmospheric grid
  public :: atmdrv       ! read atmospheric data
!
! !REVISION HISTORY:
! Created by Gordon Bonan, Sam Levis and Mariana Vertenstein
! 2005.12.15 T Craig Updated
!
!EOP
!
! PRIVATE MEMBER FUNCTIONS:
  private :: atm_openfile  ! open atmospheric forcing netCDF file
  private :: atm_readdata  ! read atmospheric forcing data
  private :: interpa2s     ! area average fields from atmosphere to surface grid
!
! PRIVATE TYPES:
  private
!
! logical variables for file manipuation
!
  logical :: open_data=.true.          !true => open data file (first tstep &
                                       !        of the run or month)
  logical :: allocated_data=.false.    !true => allocate dynamic data
!
! datafile grid data
!
  integer  :: datlon                   !number of data longitudes
  integer  :: datlat                   !number of data latitudes
!
! primary field names
!
  integer,parameter:: aV_d2a_size=14
  character(len=8) :: aV_d2a_list(av_d2a_size)    ! local aV list for d2a
  data av_d2a_list( 1) /'f_txy   '/     ! bottom level temp (K)
  data av_d2a_list( 2) /'f_uxy   '/     ! bottom level zonal wind (m/s)
  data av_d2a_list( 3) /'f_vxy   '/     ! bottom level meridional wind (m/s)
  data av_d2a_list( 4) /'f_qxy   '/     ! bottom level specific humidity (kg/kg)
  data av_d2a_list( 5) /'zgcmxy  '/     ! bottom level height above surface (m)
  data av_d2a_list( 6) /'prcxy   '/     ! convective precip rate (mm H2O/s)
  data av_d2a_list( 7) /'prlxy   '/     ! large-scale precip rate (mm H2O/s)
  data av_d2a_list( 8) /'flwdsxy '/     ! downward longwave rad onto surface (W/m**2)
  data av_d2a_list( 9) /'f_sols  '/     ! vis direct beam solar rad onto srf (W/m**2)
  data av_d2a_list(10) /'f_soll  '/     ! nir direct beam solar rad onto srf (W/m**2)
  data av_d2a_list(11) /'f_solsd '/     ! vis diffuse solar rad onto srf (W/m**2)
  data av_d2a_list(12) /'f_solld '/     ! nir diffuse solar rad onto srf (W/m**2)
  data av_d2a_list(13) /'f_pbotxy'/     ! bottom level pressure (Pa)
  data av_d2a_list(14) /'f_psrfxy'/     ! surface pressure (Pa)
  integer :: if_txy, if_uxy, if_vxy, if_qxy, izgcmxy, iprcxy, iprlxy, iflwdsxy, &
             if_sols, if_soll, if_solsd, if_solld, if_pbotxy, if_psrfxy

  type(mct_gsMap) :: gsMap_drv_glo0
  type(mct_aVect) :: aV_drv_d2a, aV_atm_d2a
  type(mct_sMat)  :: sMat0_d2a
  type(mct_sMatP) :: sMatP_d2a
  logical         :: usevector=.false.

!
! atmospheric forcing variables on raw data grid
!
  real(r8), allocatable :: x(:,:,:)            !temp. array in which atm data is stored
!
! file netCDF id's
!
  integer :: ncid                !netCDF dataset id
  integer :: nvar                !number of variables in the data file
  integer :: nlon                !number of atm longitude points
  integer :: nlat                !number of atm latitude points
  integer :: ntim                !number of atm time slices per data file
  character(len=8) :: varnam(99) !variable names of atm. fields

!-----------------------------------------------------------------------

contains

!=======================================================================

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: atmdrv
!
! !INTERFACE:
  subroutine atmdrv(nstep)
!
! !DESCRIPTION:
! This code reads in atmospheric fields from an input file and generates
! the required atmospheric forcing. These data files have [atmmin] minute
! average data for each month. Input data files are named in month-year
! format (e.g., 09-0001 contains 240 3-hour time slices of data, 30*8, for
! September of year one). The model will cycle through however many full
! years of data are available [pyr]. At least one full year of data is
! necessary for cycling. The model may start on any base date, as long as
! this date corresponds to an existing data file. A run need not be an
! exact multiple of a year.
!
! ============================
! Possible atmospheric fields:
! ============================
! Name     Description                              Required/Optional
! -----------------------------------------------------------------------------
! TBOT     temperature (K)                          Required
! WIND     wind:sqrt(u**2+v**2) (m/s)               Required
! QBOT     specific humidity (kg/kg)                Required
! Tdew     dewpoint temperature (K)                 Alternative to Q
! RH       relative humidity (percent)              Alternative to Q
! ZBOT     reference height (m)                     optional
! PSRF     surface pressure (Pa)                    optional
! FSDS     total incident solar radiation (W/m**2)  Required
! FSDSdir  direct incident solar radiation (W/m**2) optional (replaces FSDS)
! FSDSdif  diffuse incident solar rad (W/m**2)      optional (replaces FSDS)
! FLDS     incident longwave radiation (W/m**2)     optional
! PRECTmms total precipitation (mm H2O / sec)       Required
! PRECCmms convective precipitation (mm H2O / sec)  optional (replaces PRECT)
! PRECLmms large-scale precipitation (mm H2O / sec) optional (replaces PRECT)
!
! ============
! Data format:
! ============
! Data format is netCDF with dimensions longitude x latitude
! for each time slice and field. Variable names can be as in above list
! or can be reset to desired names using [fldlst] in code below.
!
! ===============
! Namelist input:
! ===============
! character*256 offline_atmdir = directory for input atm data files (can be Mass Store)
!
! !USES:
    use nanMod
    use decompMod   , only : adecomp, get_proc_bounds_atm
    use clm_atmlnd  , only : clm_mapa2l, atm_a2l, clm_a2l
    use clm_varctl  , only : offline_atmdir, pertlim
    use clm_varcon  , only : rair, cpair, co2_ppmv_const, o2_molar_const, tcrit, c13ratio
    use clm_time_manager, only : get_step_size, get_curr_calday, get_curr_date
    use fileutils   , only : getfil
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: nstep    !current time step
!
! !REVISION HISTORY:
! Created by Sam Levis
!
!EOP
!
! LOCAL VARIABLES:
    integer :: i,j,n,k,g,g1           !indices
    integer :: itimlast               !last time index used in atmrd
    real(r8):: calday                 !calendar day at Greenwich (1.00 -> 365.99)
    integer :: kda                    !day (1 -> 31)
    integer :: kmo                    !month (1 -> 12)
    integer :: kyr                    !year (0 -> ...)
    integer :: ksec                   !current seconds of current date (0 -> 86400)
    integer :: mcdate                 !current date in integer format [yyyymmdd]
    integer :: mcsec                  !current time of day [seconds]
    integer :: dtime                  !time step size
    integer :: minpday = 1440         !minutes per day
    integer :: secpmin = 60           !seconds per minute
    integer, SAVE :: itim             !time index used in atmrd
    integer, SAVE :: atmmin           !temporal resolution of atm data (in minutes)
    character(len=256), SAVE :: locfn !full file name in case atmdir is in MSS
    real(r8):: coefb        ! Slope of "Alta" expression for dependence of flfall on temp
    real(r8):: coefa        ! Offset of  of "Alta" expression for dependence of flfall on temp
    integer :: begg, endg   ! per-proc gridcell ending gridcell indices
!------------------------------------------------------------------------

    ! Determine necessary indices

    call get_proc_bounds_atm(begg, endg)

    ! -----------------------------------------------------------------
    ! Open netcdf file and read data every [atmmin] minutes
    ! -----------------------------------------------------------------

    ! Calendar information for current [nstep]

    dtime = get_step_size()
    calday = get_curr_calday()
    call get_curr_date(kyr, kmo, kda, mcsec)
    mcdate = kyr*10000 + kmo*100 + kda

    ! If 1st tstep of the run or of the month, then enter next if-block
    ! Rest flag to open file and set flag to read data

    locfn = " "
    if (open_data) then
       call atm_openfile (kda, kmo, kyr, locfn, itim, atmmin)
    endif

    ! Calculate time index

    itimlast = itim
    itim = 1 + ((kda - 1)*minpday + (mcsec - dtime)/secpmin)/atmmin
    if (dtime == int(atmmin*secpmin)) then
       if (kda == 1 .and. mcsec == 0) itim = itimlast+1
    else
       if (kda == 1 .and. mcsec == 0) itim = itimlast
    endif

    ! Determine if new data is to be read

    if (open_data .or. mod(nstep-1,atmmin*secpmin/dtime) ==0 ) then

       ! Read data for current time slice

       call t_startf('atmdread')

       if (masterproc) then
          write (6,*)
          write (6,'(72a1)') ("-",i=1,60)
          write (6,*)'nstep= ',nstep,' date= ',mcdate,' sec= ',mcsec
          if ( len_trim(locfn) > 0 )then
             write (6,*)'ATMDRV: attempting to read data from ',trim(locfn)
          end if
          write (6,'(72a1)') ("-",i=1,60)
          write (6,*)
       endif
       call atm_readdata (locfn, kmo, itim)
       call t_stopf('atmdread')

       ! Map 2d atmospheric fields from atmospheric grid to land model surface grid.
       ! Area-average absolute value of winds (i.e., regardless of
       ! direction) since land model cares about magnitude not direction.
       ! Then need to adjust resultant stresses for direction of wind.

       call t_startf('atmdinterp')

       call interpa2s

       ! Map data fields to atm model: [datlon] x [datlat] grid ->
       ! [numland] vector of land points -> [numpatch] vector of subgrid patches

!$OMP PARALLEL DO PRIVATE (g,i,j,n)
#if !defined (USE_OMP)
!CSD$ PARALLEL DO PRIVATE (g,i,j,n)
#endif
!dir$ concurrent
!cdir nodep

       do g = begg, endg
          g1 = g - begg + 1
          i = adecomp%gdc2i(g)
          j = adecomp%gdc2j(g)
          n = adecomp%gdc2glo(g)

          !States
          atm_a2l%forc_t(g) = aV_atm_d2a%rAttr(if_txy,g1)
          atm_a2l%forc_u(g) = aV_atm_d2a%rAttr(if_uxy,g1)
          atm_a2l%forc_v(g) = aV_atm_d2a%rAttr(if_vxy,g1)
          atm_a2l%forc_wind(g) = sqrt(aV_atm_d2a%rAttr(if_uxy,g1)**2 + aV_atm_d2a%rAttr(if_vxy,g1)**2)
          atm_a2l%forc_q(g) = aV_atm_d2a%rAttr(if_qxy,g1)
          atm_a2l%forc_hgt(g) = aV_atm_d2a%rAttr(izgcmxy,g1)
          atm_a2l%forc_hgt_u(g) = aV_atm_d2a%rAttr(izgcmxy,g1) !observational height of wind [m]
          atm_a2l%forc_hgt_t(g) = aV_atm_d2a%rAttr(izgcmxy,g1) !observational height of temp [m]
          atm_a2l%forc_hgt_q(g) = aV_atm_d2a%rAttr(izgcmxy,g1) !observational height of humidity [m]
          atm_a2l%forc_pbot(g) = aV_atm_d2a%rAttr(if_pbotxy,g1)
          atm_a2l%forc_psrf(g) = aV_atm_d2a%rAttr(if_psrfxy,g1)
          atm_a2l%forc_th(g)  = atm_a2l%forc_t(g) * (atm_a2l%forc_psrf(g) &
               / atm_a2l%forc_pbot(g))**(rair/cpair)
          atm_a2l%forc_vp(g)  = atm_a2l%forc_q(g) * atm_a2l%forc_pbot(g) &
               / (0.622_r8 + 0.378_r8 * atm_a2l%forc_q(g))
          atm_a2l%forc_rho(g) = (atm_a2l%forc_pbot(g) - 0.378_r8 * atm_a2l%forc_vp(g)) &
               / (rair * atm_a2l%forc_t(g))

          !BGC tracers

          atm_a2l%forc_pco2(g) = co2_ppmv_const * 1.e-6_r8 * atm_a2l%forc_pbot(g)
          atm_a2l%forc_po2(g)  = o2_molar_const * atm_a2l%forc_pbot(g)
          ! 4/14/05: PET
          ! Adding isotope code
          atm_a2l%forc_pc13o2(g) = co2_ppmv_const * c13ratio * 1.e-6_r8 * atm_a2l%forc_pbot(g)

          !Fluxes

          atm_a2l%forc_lwrad(g) = aV_atm_d2a%rAttr(iflwdsxy,g1)
          atm_a2l%forc_solad(g,1) = aV_atm_d2a%rAttr(if_sols,g1)
          atm_a2l%forc_solad(g,2) = aV_atm_d2a%rAttr(if_soll,g1)
          atm_a2l%forc_solai(g,1) = aV_atm_d2a%rAttr(if_solsd,g1)
          atm_a2l%forc_solai(g,2) = aV_atm_d2a%rAttr(if_solld,g1)
          atm_a2l%forc_solar(g) = atm_a2l%forc_solad(g,1) + atm_a2l%forc_solad(g,2) &
               + atm_a2l%forc_solai(g,1) + atm_a2l%forc_solai(g,2)

          ! Snow and Rain
          ! Set upper limit of air temperature for snowfall at 275.65K.
          ! This cut-off was selected based on Fig. 1, Plate 3-1, of Snow
          ! Hydrology (1956).

          if (aV_atm_d2a%rAttr(iprcxy,g1) + aV_atm_d2a%rAttr(iprlxy,g1) > 0._r8) then
             if (atm_a2l%forc_t(g) > (SHR_CONST_TKFRZ + tcrit)) then
                atm_a2l%forc_rain(g) = aV_atm_d2a%rAttr(iprcxy,g1) + aV_atm_d2a%rAttr(iprlxy,g1)
                atm_a2l%forc_snow(g) = 0._r8
                atm_a2l%flfall(g) = 1._r8
             else
                atm_a2l%forc_rain(g) = 0._r8
                atm_a2l%forc_snow(g) = aV_atm_d2a%rAttr(iprcxy,g1) + aV_atm_d2a%rAttr(iprlxy,g1)
                if (atm_a2l%forc_t(g) <= SHR_CONST_TKFRZ) then
                   atm_a2l%flfall(g) = 0._r8
                else if (atm_a2l%forc_t(g) <= SHR_CONST_TKFRZ+2._r8) then
                   atm_a2l%flfall(g) = -0.2_r8*SHR_CONST_TKFRZ + 0.2_r8*atm_a2l%forc_t(g)
                else
                   atm_a2l%flfall(g) = 0.4_r8
                endif
             endif
          else
             atm_a2l%forc_rain(g) = 0._r8
             atm_a2l%forc_snow(g) = 0._r8
             atm_a2l%flfall(g) = 1._r8
          endif

       end do

       call clm_mapa2l(atm_a2l, clm_a2l)

       call t_stopf('atmdinterp')

    end if

    ! Reset open_data

    if (open_data) then
       open_data = .false.    !reset to false
    elseif (kda == 1 .and. mcsec == 0) then
       open_data = .true.     !for next time step
    endif

  end subroutine atmdrv

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: atmdrv_init
!
! !INTERFACE:
  subroutine atmdrv_init()
!
! !DESCRIPTION:
! Read atmospheric grid
!
! !USES:
    use nanMod
    use clm_varctl  , only : offline_atmdir
    use domainMod   , only : alatlon, latlon_type, latlon_check, latlon_clean
    use surfrdMod   , only : surfrd_get_latlon
    use decompMod   , only : adecomp
    use areaMod     , only : celledge, cellarea,map_setmapsAR
    use fileutils   , only : getfil
    use clm_time_manager, only : get_curr_date
    use ncdio
!
! !ARGUMENTS:
    implicit none
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!
! LOCAL VARIABLES:
    type(latlon_type)  :: dlatlon        ! data file domain
    integer :: kda                !day (1 -> 31)
    integer :: kmo                !month (1 -> 12)
    integer :: kyr                !year (0 -> ...)
    integer :: ksec               !current seconds of current date (0 -> 86400)
    integer :: mcsec              !current time of day [seconds]
    character(len=  7) :: ext     !month-year extension, e.g., 01-0005
    character(len=256) :: filenam !full file name, atmdir + ext
    character(len=256) :: locfn   !full file name in case atmdir is in MSS
    logical :: lexist             !true => file exists, used when looking for a file
    integer :: dimid              !netCDF dimension id
    integer :: varid              !netCDF variable id
    integer :: ier                !error status
    character(len=32) :: subname = 'atmdrv_init'
    integer :: atmlon,atmlat      !size of alatlon
    real(r8), allocatable :: mask_d(:)   !dummy field: atm grid mask
    real(r8), allocatable :: mask_a(:)   !dummy field: land model grid mask
    character(len=256) :: str     ! string
    integer :: n                  ! generic index
    integer :: ns                 ! size
    integer :: ngseg              ! gsmap size
    integer :: root               ! root pe number
    integer,allocatable :: start(:),length(:),pe_loc(:)  ! for gsmap
!------------------------------------------------------------------------

    atmlon = alatlon%ni
    atmlat = alatlon%nj

    ! ----------------------------------------------------------------------
    ! Read offline grid data and allocate dynamic memory
    ! ----------------------------------------------------------------------

    ! Build [month]-[year] extension for file name to be read
    ! append extension to path name to get full file name

    call get_curr_date(kyr, kmo, kda, mcsec)
    write (ext,'(i4.4,"-",i2.2)') kyr,kmo
    filenam = trim(offline_atmdir) // '/' // ext // '.nc'

    call surfrd_get_latlon(dlatlon, filenam)
    call latlon_check(dlatlon)

    datlon = dlatlon%ni
    datlat = dlatlon%nj

    ! Initialize gsmaps and attr vectors
    ngseg = 1
    root = 0
    allocate(start(ngseg),length(ngseg),pe_loc(ngseg))
    start = 1
    length = datlon*datlat
    pe_loc = root
    call mct_gsMap_init(gsMap_drv_glo0,ngseg,start,length,pe_loc,root,mpicom,comp_id)
    deallocate(start,length,pe_loc)

    str = trim(av_d2a_list(1))
    do n = 2,av_d2a_size
       str = trim(str)//':'//trim(av_d2a_list(n))
    enddo

    ns = mct_gsMap_lsize(gsMap_drv_glo0, mpicom)
    call mct_aVect_init(aV_drv_d2a,rlist=str,lsize=ns)
    ns = mct_gsMap_lsize(gsMap_atm_gdc2glo, mpicom)
    call mct_aVect_init(aV_atm_d2a,rlist=str,lsize=ns)

    if_txy    = mct_aVect_indexRA(aV_drv_d2a,'f_txy'   ,perrWith=subName)
    if_uxy    = mct_aVect_indexRA(aV_drv_d2a,'f_uxy'   ,perrWith=subName)
    if_vxy    = mct_aVect_indexRA(aV_drv_d2a,'f_vxy'   ,perrWith=subName)
    if_qxy    = mct_aVect_indexRA(aV_drv_d2a,'f_qxy'   ,perrWith=subName)
    izgcmxy   = mct_aVect_indexRA(aV_drv_d2a,'zgcmxy'  ,perrWith=subName)
    iprcxy    = mct_aVect_indexRA(aV_drv_d2a,'prcxy'   ,perrWith=subName)
    iprlxy    = mct_aVect_indexRA(aV_drv_d2a,'prlxy'   ,perrWith=subName)
    iflwdsxy  = mct_aVect_indexRA(aV_drv_d2a,'flwdsxy' ,perrWith=subName)
    if_sols   = mct_aVect_indexRA(aV_drv_d2a,'f_sols'  ,perrWith=subName)
    if_soll   = mct_aVect_indexRA(aV_drv_d2a,'f_soll'  ,perrWith=subName)
    if_solsd  = mct_aVect_indexRA(aV_drv_d2a,'f_solsd' ,perrWith=subName)
    if_solld  = mct_aVect_indexRA(aV_drv_d2a,'f_solld' ,perrWith=subName)
    if_pbotxy = mct_aVect_indexRA(aV_drv_d2a,'f_pbotxy',perrWith=subName)
    if_psrfxy = mct_aVect_indexRA(aV_drv_d2a,'f_psrfxy',perrWith=subName)
    
    allocate( x(datlon,datlat,14), stat=ier)
    if (ier /= 0) then
       write (6,*) 'atmdrv_init(): allocation error _d'
       call endrun
    end if

    ! Initialize gridmap_d2a

    allocate (mask_d(datlon*datlat),mask_a(atmlon*atmlat), stat=ier)
    if (ier /= 0) then
       write (6,*) 'mask_d, mask_a allocation error'
       call endrun
    end if

    mask_d = 1._r8
    mask_a = 0._r8
    do n = 1,atmlon*atmlat
       if (adecomp%glo2gdc(n) > 0) mask_a(n) = 1._r8
    enddo
    if (masterproc) then
       call map_setmapsAR(dlatlon, alatlon, sMat0_d2a, fracin=mask_d, fracout=mask_a)
    endif

    deallocate (mask_d, mask_a)

    call mct_sMatP_init(sMatP_d2a, sMat0_d2a, &
                        gsMap_drv_glo0, gsMap_atm_gdc2glo, &
                        'Xonly',0,mpicom,comp_id)

#ifdef CPP_VECTOR
    !--- initialize the vector parts of the sMat
    call mct_sMatP_Vecinit(sMatP_d2a)
#endif

    !--- clean up the root sMat0 datatypes

    if (masterproc) then
       call mct_sMat_clean(sMat0_d2a)
    endif

    if ( masterproc )then
       write (6,*) 'Successfully made atm -> srf interpolation'
       write (6,*) 'Successfully initialized area-averaging interpolation'
       write (6,*)
    end if

    call latlon_clean(dlatlon)

  end subroutine atmdrv_init

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: atm_openfile
!
! !INTERFACE:
  subroutine atm_openfile (kda, kmo, kyr, locfn, itim, atmmin)
!
! !DESCRIPTION:
! Open atmospheric forcing netCDF file
!
! !USES:
    use clm_varctl, only : offline_atmdir
    use fileutils , only : getfil
    use ncdio
!
! !ARGUMENTS:
    implicit none
    integer, intent(in)  :: kda              !day (1 -> 31)
    integer, intent(in)  :: kmo              !month (1 -> 12)
    integer, intent(in)  :: kyr              !year (0 -> ...)
    character(len=*), intent(inout) :: locfn !history file to open and read
    integer, intent(out) :: itim             !time index used in atmrd
    integer, intent(out) :: atmmin           !temporal resolution of atm data (in minutes)
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!
! LOCAL VARIABLES:
    integer :: i,j,k,n                !do loop indices
    integer :: dimid                  !netCDF dimension id
    integer :: status                 !netCDF error status
    integer :: cyc                    !current cycle
    integer :: nyr                    !year extension to data file name, e.g., 05
    integer :: pyr = 0                !complete years of data since basedate
    integer :: nmo = 0                !number of months past basedate
    integer :: mmo = -1               !number of months past end of data
    integer :: ier                    !error status
    character(len=256) :: filenam     !full file name, atmdir + ext
    character(len=256) :: locfnlast   !full file name saved from last call
    character(len=  7) :: ext         !month-year extension, e.g., 01-0005
    logical :: lexist                 !true => file exists, used when looking for a file
    integer :: ndaypm(12) =      &    !number of days per month
             (/ 31,28,31,30,31,30,31,31,30,31,30,31 /)
    integer :: minpday = 1440         !minutes per day
    character(len=32) :: subname = 'atm_openfil'
!------------------------------------------------------------------------

    if (masterproc) then

       ! Build [month]-[year] extension for file name to be read
       ! append extension to path name to get full file name
       ! first save old locfn, then obtain the file if it exists, and
       ! finally check if file to be opened exists (e.g. if are at the
       ! end of the set of files to be cycled, then the next file will
       ! not exists and will have to rewind to beginning of data set

       write (ext,'(i4.4,"-",i2.2)') kyr,kmo
       filenam = trim(offline_atmdir) //'/'// ext // '.nc'
       locfnlast = locfn
       call getfil(filenam, locfn, 1)
       inquire (file = locfn, exist = lexist)

       ! If file exists...
       !   makes sure that if run was restarted within the same month
       !   that the original 'initial' run was started, then nmo will
       !   not count that month twice
       ! If file doesn't exist...
       !   rewind to beginning of data set (by repeating the process of
       !   creating extension, appending to path name and looking for
       !   the file

       if (lexist) then
          if (locfn /= locfnlast) nmo = nmo + 1 !months past base date
          pyr = nmo / 12      !years past base date
       elseif (nmo == 0) then
          write(6,*) 'ATM_OPENFILE error: could not find initial atm datafile'
          call endrun
       else
          mmo = mmo + 1       !(months-1) past end of data
          cyc = mmo / (pyr*12) !cycles since end of data
          nyr = kyr - pyr * (cyc + 1) !rewind to beginning of dataset
          write (ext,'(i4.4,"-",i2.2)') kyr,kmo
          filenam = trim(offline_atmdir) //'/'// ext // '.nc'
          call getfil(filenam, locfn, 1)
          inquire (file = locfn, exist = lexist)
          if (.not. lexist) then
             if (nmo < 12) then
                write(6,*) 'ATM_OPENFILE error: You must supply at least a'
                write(6,*) 'year of input data if you wish to'
                write(6,*) 'cycle through it more than once'
             else
                write(6,*) 'ATM_OPENFILE error: Not finding ',locfn
             end if
             call endrun
          end if
       end if

       ! Open netCDF data file and get lengths of lat,lon,time dimensions
       ! Do this only at the first timestep of the run or of the month

       call check_ret(nf_open (locfn, nf_nowrite, ncid), subname)

       call check_ret(nf_inq_dimid  (ncid, 'lon', dimid), subname)
       call check_ret(nf_inq_dimlen (ncid, dimid, nlon), subname)

       call check_ret(nf_inq_dimid  (ncid, 'lat', dimid), subname)
       call check_ret(nf_inq_dimlen (ncid, dimid, nlat), subname)

       call check_ret(nf_inq_dimid  (ncid, 'time', dimid), subname)
       call check_ret(nf_inq_dimlen (ncid, dimid, ntim), subname)

       if (ntim == 0) then
          write (6,*) 'ATM_OPENFILE error: zero input time slices'
          call endrun
       end if
       if (nlon /= datlon) then
          write (6,*) 'ATM_OPENFILE error: nlon = ',nlon, &
               ' in data file not equal to datlon = ',datlon,' first read in'
          call endrun
       end if
       if (nlat /= datlat) then
          write (6,*) 'ATMRD error: nlat = ',nlat, &
               ' in data file not equal to datlat = ',datlat,' first read in'
          call endrun
       end if

       ! Get variable names

       status = nf_inq_nvars(ncid, nvar)
       if (status /= nf_noerr) then
          write (6,*) ' ATM_OPENFILE netcdf error = ',nf_strerror(status)
          call endrun
       end if
       do i = 1, nvar
          call check_ret(nf_inq_varname (ncid, i, varnam(i)), subname)
       end do

    endif     !end of if-masterproc block

    ! Calculate temporal resolution of the dataset in minutes and
    ! reset time index

    call mpi_bcast (ntim, 1, MPI_INTEGER, 0, mpicom, ier)
    atmmin = ndaypm(kmo) * minpday / ntim
    itim = 0

  end subroutine atm_openfile

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: atm_readdata
!
! !INTERFACE:
  subroutine atm_readdata (fname, kmo, itim)
!
! !DESCRIPTION:
! read atmospheric forcing single level fields from netCDF file
! this is done at every subsequent call to atmrd
! this includes a second call at the 1st tstep of the run or of the month
! Use nf_get_vara_real to read data for variable referenced by
! variable id = [i]
! nf_get_vara_real (ncid, i, beg, len, varval)
! integer beg - a vector of integers specifying the index in the
!               variable where the first data value is read.
!               must be dimensioned same as the variable's dimension
!               and a starting value must be given for each dimension
! integer len - a vector of integers specifying the number of data
!               values, for each of the variable's dimensions, to read.
!               must be dimensioned same as the variable's dimension
!               and a length must be given for each dimension
!
! !USES:
    use clm_varcon, only : sb
    use fileutils , only : getfil
    use ncdio
!
! !ARGUMENTS:
    implicit none
    character(len=*), intent(in) :: fname           !history file to open and read
    integer, intent(in)  :: kmo, itim               !current month and time index
!
! !REVISION HISTORY:
! Created by Sam Levis
!
!EOP
!
! LOCAL VARIABLES:
    integer :: i,j,k,n                  !do loop indices
    integer :: ier                      !error status
    integer :: varid                    !netCDF variable id
    integer :: status                   !netCDF error status
    integer :: beg3d(3)                 !netCDF 3-d start index (where to read first value)
    integer :: len3d(3)                 !netCDF 3-d count index (number of values to read)
    character(len=32) :: subname = 'atm_readdata'
!
! atm input field names
!
    real(r8) ea                    !atmospheric emissivity

    logical atmread_err

    ! use polynomials to calculate saturation vapor pressure and derivative with
    ! respect to temperature: over water when t > 0 c and over ice when t <= 0 c
    ! required to convert relative humidity to specific humidity

    real(r8) esatw                 !saturation vapor pressure over water (Pa)
    real(r8) esati                 !saturation vapor pressure over ice (Pa)
    real(r8) e                     !vapor pressure (Pa)
    real(r8) qsat                  !saturation specific humidity (kg/kg)
    real(r8) a0,a1,a2,a3,a4,a5,a6  !coefficients for esat over water
    real(r8) b0,b1,b2,b3,b4,b5,b6  !coefficients for esat over ice
    real(r8) tdc, t                !Kelvins to Celcius function and its input

    parameter (a0=6.107799961_r8    , a1=4.436518521e-01_r8, &
               a2=1.428945805e-02_r8, a3=2.650648471e-04_r8, &
               a4=3.031240396e-06_r8, a5=2.034080948e-08_r8, &
               a6=6.136820929e-11_r8)

    parameter (b0=6.109177956_r8    , b1=5.034698970e-01_r8, &
               b2=1.886013408e-02_r8, b3=4.176223716e-04_r8, &
               b4=5.824720280e-06_r8, b5=4.838803174e-08_r8, &
               b6=1.838826904e-10_r8)

    integer,parameter:: fldsize=14
    character(len=8) :: fldlst(fldsize)      !name of possible atm fields in input file
    data fldlst( 1) /'TBOT    '/
    data fldlst( 2) /'WIND    '/
    data fldlst( 3) /'QBOT    '/
    data fldlst( 4) /'Tdew    '/
    data fldlst( 5) /'RH      '/
    data fldlst( 6) /'ZBOT    '/
    data fldlst( 7) /'PSRF    '/
    data fldlst( 8) /'FSDS    '/
    data fldlst( 9) /'FSDSdir '/
    data fldlst(10) /'FSDSdif '/
    data fldlst(11) /'FLDS    '/
    data fldlst(12) /'PRECTmms'/
    data fldlst(13) /'PRECCmms'/
    data fldlst(14) /'PRECLmms'/
!
! function declarations
!
    tdc(t) = min( 50._r8, max(-50._r8,(t-SHR_CONST_TKFRZ)) )
    esatw(t) = 100._r8*(a0+t*(a1+t*(a2+t*(a3+t*(a4+t*(a5+t*a6))))))
    esati(t) = 100._r8*(b0+t*(b1+t*(b2+t*(b3+t*(b4+t*(b5+t*b6))))))
!------------------------------------------------------------------------



    ! Read single level fields

    if (masterproc) then

       ! initialize fields to the flag value

       x(:,:,:) = -1._r8

       ! read input data single-level fields

       beg3d(1) = 1     ;  len3d(1) = datlon
       beg3d(2) = 1     ;  len3d(2) = datlat
       beg3d(3) = itim  ;  len3d(3) = 1
       do k = 1, 14
          do n = 1, nvar
             if (varnam(n) == fldlst(k)) then
                call check_ret(nf_get_vara_double(ncid,n,beg3d,len3d,x(1,1,k)), subname)
             end if
          end do              !end loop of fields in input file
       end do                 !end loop of fields expected in input file

       ! Close file at the end of the month
       ! NOTE: as written will not close file if run ends mid-month

       if (itim == ntim) then
          call check_ret(nf_close (ncid), subname)
          write (6,*) '---------------------------------------'
          write (6,*) 'ATMRD: closing data for ',trim(fname)
          write (6,*) '---------------------------------------'
          write (6,*)
       end if

    endif     !end of if-masterproc block

    ! ----------------------------------------------------------------------
    ! Determine 2d atmospheric fields
    ! Follow order in fldlst(14) to determine what was read and what was not
    ! ----------------------------------------------------------------------

    ! Loop over atmospheric longitudes and latitudes

    if (masterproc) then

    atmread_err = .false.
!$OMP PARALLEL DO PRIVATE (i,j,e,ea,qsat,n)
#if !defined (USE_OMP)
!CSD$ PARALLEL DO PRIVATE (i,j,e,ea,qsat,n)
#endif
    do j = 1, datlat
       do i = 1, datlon
          n = (j-1)*datlon + i
          ! FORC_TXY
          if (nint(x(i,j,1)) == -1) then
             write(6,*)'ATM error: TBOT has not been read by atmrd'
             atmread_err = .true.
          else if (x(i,j,1) < 50._r8) then
             write(6,*)'ATM error: TBOT appears to be in deg C'
             write(6,*)'Converting to Kelvins now'
             aV_drv_d2a%rAttr(if_txy,n) = x(i,j,1) + SHR_CONST_TKFRZ
          else
             aV_drv_d2a%rAttr(if_txy,n) = x(i,j,1)
          end if

          ! FORC_UXY, FORC_VXY
          if (nint(x(i,j,2)) == -1) then
             write(6,*)'ATM error: WIND has not been read by atmrd'
             atmread_err = .true.
          else
             aV_drv_d2a%rAttr(if_uxy,n) = x(i,j,2) / sqrt(2._r8)
             aV_drv_d2a%rAttr(if_vxy,n) = x(i,j,2) / sqrt(2._r8)
          end if

          ! FORC_PSRFXY, FORC_PBOTXY
          if (nint(x(i,j,7)) == -1) then
             aV_drv_d2a%rAttr(if_psrfxy,n) = SHR_CONST_PSTD
          else
             aV_drv_d2a%rAttr(if_psrfxy,n) = x(i,j,7)
          end if
          aV_drv_d2a%rAttr(if_pbotxy,n)  = aV_drv_d2a%rAttr(if_psrfxy,n)

          !FORC_QXY
          if (nint(x(i,j,3)) == -1) then
             if (nint(x(i,j,4)) == -1) then
                if (nint(x(i,j,5)) == -1) then
                   write(6,*)'ATM error: Humidity has not been'
                   write(6,*)'read by atmrd'
                   atmread_err = .true.
                else          !using RH as %
                   if (aV_drv_d2a%rAttr(if_txy,n) > SHR_CONST_TKFRZ) then
                      e = x(i,j,5)/100._r8 * esatw(tdc(aV_drv_d2a%rAttr(if_txy,n)))
                   else
                      e = x(i,j,5)/100._r8 * esati(tdc(aV_drv_d2a%rAttr(if_txy,n)))
                   end if
                end if
                aV_drv_d2a%rAttr(if_qxy,n) = 0.622_r8*e / (aV_drv_d2a%rAttr(if_pbotxy,n) - 0.378_r8*e)
             else             !using Tdew
                if (x(i,j,4) < 50._r8) then
                   write(6,*)'ATM warning: Tdew appears to be in'
                   write(6,*)'deg C, so converting to Kelvin'
                   x(i,j,4) = x(i,j,4) + SHR_CONST_TKFRZ
                end if
                if (x(i,j,4) > aV_drv_d2a%rAttr(if_txy,n)) then
                   write(6,*)'ATM warning: Dewpt temp > temp!'
                end if
                if (x(i,j,4) > SHR_CONST_TKFRZ) then
                   e = esatw(tdc(x(i,j,4)))
                else
                   e = esati(tdc(x(i,j,4)))
                end if
                aV_drv_d2a%rAttr(if_qxy,n) = 0.622_r8*e / (aV_drv_d2a%rAttr(if_pbotxy,n) - 0.378_r8*e)
             end if
          else                !using QBOT in kg/kg
             if (aV_drv_d2a%rAttr(if_txy,n) > SHR_CONST_TKFRZ) then
                e = esatw(tdc(aV_drv_d2a%rAttr(if_txy,n)))
             else
                e = esati(tdc(aV_drv_d2a%rAttr(if_txy,n)))
             end if
             qsat = 0.622_r8*e / (aV_drv_d2a%rAttr(if_pbotxy,n) - 0.378_r8*e)
             if (qsat < x(i,j,3)) then
                aV_drv_d2a%rAttr(if_qxy,n) = qsat
!                 write(6,*)'ATM warning: qsat < q!'
             else
                aV_drv_d2a%rAttr(if_qxy,n) = x(i,j,3)
             end if
          end if

          ! ZGCMXY
          if (nint(x(i,j,6)) == -1) then
             aV_drv_d2a%rAttr(izgcmxy,n) = 30._r8
          else
             aV_drv_d2a%rAttr(izgcmxy,n) = x(i,j,6)
          end if

          ! FORC_SOLSXY, FORC_SOLLXY, FORC_SOLSDXY, FORC_SOLLDXY

          if (nint(x(i,j,9))==-1.or.nint(x(i,j,10))==-1) then
             if (nint(x(i,j,8)) /= -1) then
                aV_drv_d2a%rAttr(if_sols,n)  = 0.7_r8 * (0.5_r8 * x(i,j,8))
                aV_drv_d2a%rAttr(if_soll,n)  = aV_drv_d2a%rAttr(if_sols,n)
                aV_drv_d2a%rAttr(if_solsd,n) = 0.3_r8 * (0.5_r8 * x(i,j,8))
                aV_drv_d2a%rAttr(if_solld,n) = aV_drv_d2a%rAttr(if_solsd,n)
             else
                write(6,*)'ATM error: neither FSDSdir/dif nor'
                write(6,*)'       FSDS have been read in by atmrd'
                atmread_err = .true.
             end if
          else
             aV_drv_d2a%rAttr(if_sols,n)  = 0.5_r8 * x(i,j,9)
             aV_drv_d2a%rAttr(if_soll,n)  = aV_drv_d2a%rAttr(if_sols,n)
             aV_drv_d2a%rAttr(if_solsd,n) = 0.5_r8 * x(i,j,10)
             aV_drv_d2a%rAttr(if_solld,n) = aV_drv_d2a%rAttr(if_solsd,n)
          end if

          ! PRCXY, PRLXY

          if (nint(x(i,j,13))==-1.or.nint(x(i,j,14))==-1) then
             if (nint(x(i,j,12)).ne.-1) then
                aV_drv_d2a%rAttr(iprcxy,n) = 0.1_r8 * x(i,j,12)
                aV_drv_d2a%rAttr(iprlxy,n) = 0.9_r8 * x(i,j,12)
             else
                write(6,*)'ATM error: neither PRECC/L nor PRECT'
                write(6,*)'           have been read in by atmrd'
                atmread_err = .true.
             end if
          else
             aV_drv_d2a%rAttr(iprcxy,n) = x(i,j,13)
             aV_drv_d2a%rAttr(iprlxy,n) = x(i,j,14)
          end if

          ! FLWDSXY

          if (nint(x(i,j,11)) == -1) then
             e = aV_drv_d2a%rAttr(if_psrfxy,n) * aV_drv_d2a%rAttr(if_qxy,n) / (0.622_r8 + 0.378_r8 * aV_drv_d2a%rAttr(if_qxy,n))
             ea = 0.70_r8 + 5.95e-05_r8 * 0.01_r8*e * exp(1500.0_r8/aV_drv_d2a%rAttr(if_txy,n))
             aV_drv_d2a%rAttr(iflwdsxy,n) = ea * sb * aV_drv_d2a%rAttr(if_txy,n)**4
          else
             aV_drv_d2a%rAttr(iflwdsxy,n) = x(i,j,11)
          end if

       end do                 !end loop of latitudes
    end do                    !end loop of longitudes
#if !defined (USE_OMP)
!CSD$ END PARALLEL DO
#endif
!$OMP END PARALLEL DO

    if (atmread_err) then
       write(6,*) 'atm_readdata: error reading atm data'
       call endrun
    end if

    endif   ! masterproc

  end subroutine atm_readdata

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: interpa2s
!
! !INTERFACE:
  subroutine interpa2s ()
!
! !DESCRIPTION:
! Area average fields from atmosphere grid to surface grid
!
! !USES:
!
! !ARGUMENTS:
    implicit none
!
! !REVISION HISTORY:
! Created by Gordon Bonan
!
!EOP
!
! LOCAL VARIABLES:
    integer  :: i,j,n,g,k,g1              !longitude,latitude loop indices
    real(r8),target :: forc_u(datlon*datlat)  !dummy wind (u)
    real(r8),target :: forc_v(datlon*datlat)  !dummy wind (v)
    integer  :: cnt
    real(r8),pointer :: farray(:)
    integer  :: begg,endg
!------------------------------------------------------------------------

    ! area-average absolute value of winds (i.e., regardless of
    ! direction) since land model cares about magnitude not direction.
    ! then need to adjust resultant stresses for direction of wind.

    if (masterproc) then
       do j = 1, datlat
       do i = 1, datlon
          n = (j-1)*datlon + i
          av_drv_d2a%rAttr(if_uxy,n) = abs(av_drv_d2a%rAttr(if_uxy,n))
          av_drv_d2a%rAttr(if_vxy,n) = abs(av_drv_d2a%rAttr(if_vxy,n))
       end do
       end do
    endif

    call mct_Smat_AvMult(av_drv_d2a, sMatP_d2a, av_atm_d2a, vector=usevector)
    call mct_aVect_unpermute(av_atm_d2a, perm_atm_gdc2glo)

  end subroutine interpa2s

#endif

end module atmdrvMod

