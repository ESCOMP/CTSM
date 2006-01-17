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
  use areaMod      , only : gridmap_type
  use domainMod    , only : domain_type
  use abortutils   , only : endrun
#if (defined SPMD)
  use spmdMod      , only : masterproc, mpicom, MPI_REAL8, MPI_INTEGER
#else
  use spmdMod      , only : masterproc
#endif
  use shr_sys_mod  , only : shr_sys_flush
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
  private :: interpa2si    ! initialize atm->land model interpolation
  private :: interpa2s     ! area average fields from atmosphere to surface grid
!
! PRIVATE TYPES:
  private

  type(domain_type)  :: ddomain        ! data file domain
  type(gridmap_type) :: gridmap_d2a    ! local grid for datafile to atm grid
!
! logical variables for file manipuation
!
  logical :: open_data=.true.          !true => open data file (first tstep &
                                       !        of the run or month)
  logical :: allocated_data=.false.    !true => allocate dynamic data
!
! datafile, ddomain grid data
!
  integer  :: datlon                   !number of data longitudes
  integer  :: datlat                   !number of data latitudes
!
! atmospheric forcing variables on raw data grid
!
  real(r8), allocatable :: x(:,:,:)            !temp. array in which atm data is stored
  real(r8), allocatable :: forc_txy_d (:,:)    !atm bottom level temperature (Kelvin)
  real(r8), allocatable :: forc_uxy_d (:,:)    !atm bottom level zonal wind (m/s)
  real(r8), allocatable :: forc_vxy_d (:,:)    !atm bottom level meridional wind (m/s)
  real(r8), allocatable :: forc_qxy_d (:,:)    !atm bottom level specific humidity (kg/kg)
  real(r8), allocatable :: zgcmxy_d (:,:)      !atm bottom level height above surface (m)
  real(r8), allocatable :: prcxy_d  (:,:)      !convective precipitation rate (mm H2O/s)
  real(r8), allocatable :: prlxy_d  (:,:)      !large-scale precipitation rate (mm H2O/s)
  real(r8), allocatable :: flwdsxy_d(:,:)      !downward longwave rad onto surface (W/m**2)
  real(r8), allocatable :: forc_solsxy_d (:,:) !vis direct beam solar rad onto srf (W/m**2)
  real(r8), allocatable :: forc_sollxy_d (:,:) !nir direct beam solar rad onto srf (W/m**2)
  real(r8), allocatable :: forc_solsdxy_d(:,:) !vis diffuse solar rad onto srf (W/m**2)
  real(r8), allocatable :: forc_solldxy_d(:,:) !nir diffuse solar rad onto srf(W/m**2)
  real(r8), allocatable :: forc_pbotxy_d (:,:) !atm bottom level pressure (Pa)
  real(r8), allocatable :: forc_psrfxy_d (:,:) !atm surface pressure (Pa)
!
! atmospheric forcing variables on atmosphere model grid
!
  real(r8), allocatable :: forc_txy_a (:,:)         !atm bottom level temperature (Kelvin)
  real(r8), allocatable :: forc_uxy_a (:,:)         !atm bottom level zonal wind (m/s)
  real(r8), allocatable :: forc_vxy_a (:,:)         !atm bottom level meridional wind (m/s)
  real(r8), allocatable :: forc_qxy_a (:,:)         !atm bottom level specific humidity (kg/kg)
  real(r8), allocatable :: zgcmxy_a (:,:)           !atm bottom level height above surface (m)
  real(r8), allocatable :: prcxy_a  (:,:)           !convective precipitation rate (mm H2O/s)
  real(r8), allocatable :: prlxy_a  (:,:)           !large-scale precipitation rate (mm H2O/s)
  real(r8), allocatable :: flwdsxy_a(:,:)           !downward longwave rad onto surface (W/m**2)
  real(r8), allocatable :: forc_solsxy_a (:,:)      !vis direct beam solar rad onto srf (W/m**2)
  real(r8), allocatable :: forc_sollxy_a (:,:)      !nir direct beam solar rad onto srf (W/m**2)
  real(r8), allocatable :: forc_solsdxy_a(:,:)      !vis diffuse solar rad onto srf (W/m**2)
  real(r8), allocatable :: forc_solldxy_a(:,:)      !nir diffuse solar rad onto srf(W/m**2)
  real(r8), allocatable :: forc_pbotxy_a (:,:)      !atm bottom level pressure (Pa)
  real(r8), allocatable :: forc_psrfxy_a (:,:)      !atm surface pressure (Pa)
!
! file netCDF id's
!
  integer :: ncid                !netCDF dataset id
  integer :: nvar                !number of variables in the data file
  integer :: nlon                !number of atm longitude points
  integer :: nlat                !number of atm latitude points
  integer :: ntim                !number of atm time slices per data file
  character(len=8) :: varnam(99) !variable names of atm. fields

#if (defined PERGRO)
  logical :: do_perturb = .true. !if true, perturb tbot
#endif
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
    use clmtype
    use decompMod   , only : get_proc_bounds
    use clm_atmlnd  , only : clm_a2l, atm_a2l, gridmap_a2l, clm_mapa2l
    use clm_varctl  , only : offline_atmdir, pertlim
    use clm_varcon  , only : rair, cpair, co2_ppmv_const, o2_molar_const, tcrit, c13ratio
    use time_manager, only : get_step_size, get_curr_calday, get_curr_date
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
    integer :: i,j,k,g                !indices
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
#if (defined PERGRO)
    real(r8):: pertval
#endif
    real(r8):: coefb        ! Slope of "Alta" expression for dependence of flfall on temp
    real(r8):: coefa        ! Offset of  of "Alta" expression for dependence of flfall on temp
    integer :: begp, endp   ! per-proc beginning and ending pft indices
    integer :: begc, endc   ! per-proc beginning and ending column indices
    integer :: begl, endl   ! per-proc beginning and ending landunit indices
    integer :: begg, endg   ! per-proc gridcell ending gridcell indices
    type(gridcell_type), pointer :: gptr  ! pointer to gridcell derived subtype
!------------------------------------------------------------------------

    ! Determine necessary indices

    call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)

    ! Set pointers into derived type

    gptr => clm3%g

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

       if (masterproc) then
          write (6,*)
          write (6,'(72a1)') ("-",i=1,60)
          write (6,*)'nstep= ',nstep,' date= ',mcdate,' sec= ',mcsec
          if ( len_trim(locfn) > 0 )then
             write (6,*)'ATM: attempting to read data from ',trim(locfn)
          end if
          write (6,'(72a1)') ("-",i=1,60)
          write (6,*)
       endif
       call atm_readdata (locfn, kmo, itim)

       ! Map 2d atmospheric fields from atmospheric grid to land model surface grid.
       ! Area-average absolute value of winds (i.e., regardless of
       ! direction) since land model cares about magnitude not direction.
       ! Then need to adjust resultant stresses for direction of wind.

       call interpa2s

       ! Map data fields to atm model: [datlon] x [datlat] grid ->
       ! [numland] vector of land points -> [numpatch] vector of subgrid patches

#if (defined PERGRO)
       if (masterproc) then
          if (pertlim /= 0.0_r8 .and. do_perturb) then
             write(6,*)'ATMDRV: Adding random perturbation bounded by +/-', &
                  pertlim,' to initial temperature field'
          end if
       endif
#endif

#if (!defined PERGRO)
!$OMP PARALLEL DO PRIVATE (g,i,j)
#if !defined (USE_OMP)
!CSD$ PARALLEL DO PRIVATE (g,i,j)
#endif
#endif
!dir$ concurrent
!cdir nodep

       do g = begg, endg
          i = gptr%ixy(g)
          j = gptr%jxy(g)

          !States

          atm_a2l%forc_t(g) = forc_txy_a(i,j)
          atm_a2l%forc_u(g) = forc_uxy_a(i,j)
          atm_a2l%forc_v(g) = forc_vxy_a(i,j)
          atm_a2l%forc_wind(g) = sqrt(forc_uxy_a(i,j)**2 + forc_vxy_a(i,j)**2)
          atm_a2l%forc_q(g) = forc_qxy_a(i,j)
          atm_a2l%forc_hgt(g) = zgcmxy_a(i,j)
          atm_a2l%forc_hgt_u(g) = zgcmxy_a(i,j) !observational height of wind [m]
          atm_a2l%forc_hgt_t(g) = zgcmxy_a(i,j) !observational height of temp [m]
          atm_a2l%forc_hgt_q(g) = zgcmxy_a(i,j) !observational height of humidity [m]
          atm_a2l%forc_pbot(g) = forc_pbotxy_a(i,j)
          atm_a2l%forc_psrf(g) = forc_psrfxy_a(i,j)
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

          atm_a2l%forc_lwrad(g) = flwdsxy_a(i,j)
          atm_a2l%forc_solad(g,1) = forc_solsxy_a(i,j)
          atm_a2l%forc_solad(g,2) = forc_sollxy_a(i,j)
          atm_a2l%forc_solai(g,1) = forc_solsdxy_a(i,j)
          atm_a2l%forc_solai(g,2) = forc_solldxy_a(i,j)
          atm_a2l%forc_solar(g) = forc_solsxy_a(i,j) + forc_sollxy_a(i,j) &
               + forc_solsdxy_a(i,j) + forc_solldxy_a(i,j)

          ! Snow and Rain
          ! Set upper limit of air temperature for snowfall at 275.65K.
          ! This cut-off was selected based on Fig. 1, Plate 3-1, of Snow
          ! Hydrology (1956).

          if (prcxy_a(i,j) + prlxy_a(i,j) > 0._r8) then
             if (atm_a2l%forc_t(g) > (SHR_CONST_TKFRZ + tcrit)) then
                atm_a2l%forc_rain(g) = prcxy_a(i,j) + prlxy_a(i,j)
                atm_a2l%forc_snow(g) = 0._r8
                atm_a2l%flfall(g) = 1._r8
             else
                atm_a2l%forc_rain(g) = 0._r8
                atm_a2l%forc_snow(g) = prcxy_a(i,j) + prlxy_a(i,j)
#if (defined PERGRO)
                ! Note for cleanup: this PERGRO block has the same functional form as the
                ! non-PERGRO block - not sure why.

                coefb = 0.4_r8/2.0_r8
                coefa = -coefb*SHR_CONST_TKFRZ
                if (atm_a2l%forc_t(g) <= SHR_CONST_TKFRZ) then
                   atm_a2l%flfall = 0.0_r8
                else if (atm_a2l%forc_t(g) <= SHR_CONST_TKFRZ+2._r8) then
                   atm_a2l%flfall(g) = coefa + coefb * atm_a2l%forc_t(g)
                else
                   atm_a2l%flfall(g) = coefa + coefb * (SHR_CONST_TKFRZ+2._r8)
                endif
#else
                if (atm_a2l%forc_t(g) <= SHR_CONST_TKFRZ) then
                   atm_a2l%flfall(g) = 0._r8
                else if (atm_a2l%forc_t(g) <= SHR_CONST_TKFRZ+2._r8) then
                   atm_a2l%flfall(g) = -54.632_r8 + 0.2_r8 * atm_a2l%forc_t(g)
                else
                   atm_a2l%flfall(g) = 0.4_r8
                endif
#endif
             endif
          else
             atm_a2l%forc_rain(g) = 0._r8
             atm_a2l%forc_snow(g) = 0._r8
             atm_a2l%flfall(g) = 1._r8
          endif

#if (defined PERGRO)
          ! Add random perturbation to temperature if required

          if (pertlim /= 0.0_r8 .and. do_perturb) then
             call random_number (pertval)
             pertval = 2._r8*pertlim*(0.5_r8 - pertval)
             atm_a2l%forc_t(g) = (atm_a2l%forc_t(g))*(1._r8 + pertval)
          endif
#endif
       end do

#if (!defined PERGRO)
#if !defined (USE_OMP)
!CSD$ END PARALLEL DO
#endif
!$OMP END PARALLEL DO
#else
       if (do_perturb) then
          do_perturb = .false.
       endif
#endif

       call clm_mapa2l(atm_a2l,clm_a2l,gridmap_a2l)

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
    use domainMod   , only : domain_init, adomain
    use areaMod     , only : celledge, cellarea
    use fileutils   , only : getfil
    use time_manager, only : get_curr_date
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
    integer :: atmlon,atmlat      !size of adomain
!------------------------------------------------------------------------

    atmlon = adomain%ni
    atmlat = adomain%nj

    ! ----------------------------------------------------------------------
    ! Read offline grid data and allocate dynamic memory
    ! ----------------------------------------------------------------------

    if (masterproc) then

       ! Build [month]-[year] extension for file name to be read
       ! append extension to path name to get full file name

       call get_curr_date(kyr, kmo, kda, mcsec)
       write (ext,'(i4.4,"-",i2.2)') kyr,kmo
       filenam = trim(offline_atmdir) // '/' // ext // '.nc'
       call getfil(filenam, locfn, 1)
       inquire (file = locfn, exist = lexist)
       if (.not. lexist) then
          write(6,*) 'ATMDRV_INIT error: could not find initial atm datafile'
          call endrun
       endif

       ! Open netCDF data file and get lengths of lat,lon,time dimensions

       call check_ret(nf_open (locfn, nf_nowrite, ncid), subname)

       call check_ret(nf_inq_dimid  (ncid, 'lon', dimid), subname)
       call check_ret(nf_inq_dimlen (ncid, dimid, datlon), subname)

       call check_ret(nf_inq_dimid  (ncid, 'lat', dimid), subname)
       call check_ret(nf_inq_dimlen (ncid, dimid, datlat), subname)

       call check_ret(nf_inq_dimid  (ncid, 'time', dimid), subname)
       call check_ret(nf_inq_dimlen (ncid, dimid, ntim), subname)
       if (ntim == 0) then
          write (6,*) 'ATMDRV_INIT error: zero input time slices'
          call endrun
       end if

    endif  !end of if-masterproc block
#if (defined SPMD)
    call mpi_bcast (datlon, 1, MPI_INTEGER, 0, mpicom, ier)
    call mpi_bcast (datlat, 1, MPI_INTEGER, 0, mpicom, ier)
#endif

    ! Allocate space for dynamic variables

    if (.not. allocated_data) then
       if (masterproc) write(6,*)' ATM_GETRID: allocating dynamic space'
       call domain_init(ddomain,datlon,datlat)
       allocate( forc_txy_d(datlon,datlat), &
                 forc_uxy_d(datlon,datlat), &
                 forc_vxy_d(datlon,datlat), &
                 forc_qxy_d(datlon,datlat), &
                 zgcmxy_d(datlon,datlat), &
                 prcxy_d(datlon,datlat), &
                 prlxy_d(datlon,datlat), &
                 flwdsxy_d(datlon,datlat), &
                 forc_solsxy_d(datlon,datlat), &
                 forc_sollxy_d(datlon,datlat), &
                 forc_solsdxy_d(datlon,datlat), &
                 forc_solldxy_d(datlon,datlat), &
                 forc_pbotxy_d(datlon,datlat), &
                 forc_psrfxy_d(datlon,datlat), &
                 x(datlon,datlat,14), stat=ier)
       if (ier /= 0) then
          write (6,*) 'atmdrv_init(): allocation error _d'
          call endrun
       end if
       allocate( forc_txy_a(atmlon,atmlat), &
                 forc_uxy_a(atmlon,atmlat), &
                 forc_vxy_a(atmlon,atmlat), &
                 forc_qxy_a(atmlon,atmlat), &
                 zgcmxy_a(atmlon,atmlat), &
                 prcxy_a(atmlon,atmlat), &
                 prlxy_a(atmlon,atmlat), &
                 flwdsxy_a(atmlon,atmlat), &
                 forc_solsxy_a(atmlon,atmlat), &
                 forc_sollxy_a(atmlon,atmlat), &
                 forc_solsdxy_a(atmlon,atmlat), &
                 forc_solldxy_a(atmlon,atmlat), &
                 forc_pbotxy_a(atmlon,atmlat), &
                 forc_psrfxy_a(atmlon,atmlat), &
                 stat=ier)
       if (ier /= 0) then
          write (6,*) 'atmdrv_init(): allocation error _a'
          call endrun
       end if
       allocated_data = .true.
    endif

    ! Extract atmospheric data grid information and close file

    if (masterproc) then

       call check_ret(nf_inq_varid(ncid, 'EDGEN', varid), subname)
       call check_ret(nf_get_var_double(ncid, varid, ddomain%edges(1)), subname)

       call check_ret(nf_inq_varid(ncid, 'EDGEE', varid), subname)
       call check_ret(nf_get_var_double(ncid, varid, ddomain%edges(2)), subname)

       call check_ret(nf_inq_varid(ncid, 'EDGES', varid), subname)
       call check_ret(nf_get_var_double(ncid, varid, ddomain%edges(3)), subname)

       call check_ret(nf_inq_varid(ncid, 'EDGEW', varid), subname)
       call check_ret(nf_get_var_double(ncid, varid, ddomain%edges(4)), subname)

       call check_ret(nf_inq_varid(ncid, 'LONGXY', varid), subname)
       call check_ret(nf_get_var_double(ncid, varid, ddomain%lonc), subname)

       call check_ret(nf_inq_varid(ncid,'LATIXY', varid), subname)
       call check_ret(nf_get_var_double(ncid, varid, ddomain%latc), subname)

       call check_ret(nf_close (ncid), subname)
       write (6,*) 'ATMDRV_INIT: closing data for ',trim(locfn)

    end if !end of if-masterproc block

#if (defined SPMD)
    call mpi_bcast (ddomain%edges,size(ddomain%edges),MPI_REAL8, 0, mpicom, ier)
    call mpi_bcast (ddomain%lonc ,size(ddomain%lonc) ,MPI_REAL8, 0, mpicom, ier)
    call mpi_bcast (ddomain%latc ,size(ddomain%latc) ,MPI_REAL8, 0, mpicom, ier)
#endif

    call celledge (ddomain, &
                   ddomain%edges(1) , ddomain%edges(2) , &
                   ddomain%edges(3) , ddomain%edges(4)  )

    call cellarea (ddomain, &
                   ddomain%edges(1) , ddomain%edges(2) , &
                   ddomain%edges(3) , ddomain%edges(4)  )

    ! Initialize

    call interpa2si

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

#if (defined SPMD)
    call mpi_bcast (ntim, 1, MPI_INTEGER, 0, mpicom, ier)
#endif
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
    character(len=8) :: fldlst(14)      !name of possible atm fields in input file
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

#if (defined SPMD)
    call mpi_bcast (x, size(x), MPI_REAL8, 0, mpicom, ier)
#endif

    ! ----------------------------------------------------------------------
    ! Determine 2d atmospheric fields
    ! Follow order in fldlst(14) to determine what was read and what was not
    ! ----------------------------------------------------------------------

    ! Loop over atmospheric longitudes and latitudes

    atmread_err = .false.
!$OMP PARALLEL DO PRIVATE (i,j,e,ea,qsat)
#if !defined (USE_OMP)
!CSD$ PARALLEL DO PRIVATE (i,j,e,ea,qsat)
#endif
    do j = 1, datlat
       do i = 1, datlon
          ! FORC_TXY
          if (nint(x(i,j,1)) == -1) then
             write(6,*)'ATM error: TBOT has not been read by atmrd'
             atmread_err = .true.
          else if (x(i,j,1) < 50._r8) then
             write(6,*)'ATM error: TBOT appears to be in deg C'
             write(6,*)'Converting to Kelvins now'
             forc_txy_d(i,j) = x(i,j,1) + SHR_CONST_TKFRZ
          else
             forc_txy_d(i,j) = x(i,j,1)
          end if

          ! FORC_UXY, FORC_VXY
          if (nint(x(i,j,2)) == -1) then
             write(6,*)'ATM error: WIND has not been read by atmrd'
             atmread_err = .true.
          else
             forc_uxy_d(i,j) = x(i,j,2) / sqrt(2._r8)
             forc_vxy_d(i,j) = x(i,j,2) / sqrt(2._r8)
          end if

          ! FORC_PSRFXY, FORC_PBOTXY
          if (nint(x(i,j,7)) == -1) then
             forc_psrfxy_d(i,j) = SHR_CONST_PSTD
          else
             forc_psrfxy_d(i,j) = x(i,j,7)
          end if
          forc_pbotxy_d(i,j)  = forc_psrfxy_d(i,j)

          !FORC_QXY
          if (nint(x(i,j,3)) == -1) then
             if (nint(x(i,j,4)) == -1) then
                if (nint(x(i,j,5)) == -1) then
                   write(6,*)'ATM error: Humidity has not been'
                   write(6,*)'read by atmrd'
                   atmread_err = .true.
                else          !using RH as %
                   if (forc_txy_d(i,j) > SHR_CONST_TKFRZ) then
                      e = x(i,j,5)/100._r8 * esatw(tdc(forc_txy_d(i,j)))
                   else
                      e = x(i,j,5)/100._r8 * esati(tdc(forc_txy_d(i,j)))
                   end if
                end if
                forc_qxy_d(i,j) = 0.622_r8*e / (forc_pbotxy_d(i,j) - 0.378_r8*e)
             else             !using Tdew
                if (x(i,j,4) < 50._r8) then
                   write(6,*)'ATM warning: Tdew appears to be in'
                   write(6,*)'deg C, so converting to Kelvin'
                   x(i,j,4) = x(i,j,4) + SHR_CONST_TKFRZ
                end if
                if (x(i,j,4) > forc_txy_d(i,j)) then
                   write(6,*)'ATM warning: Dewpt temp > temp!'
                end if
                if (x(i,j,4) > SHR_CONST_TKFRZ) then
                   e = esatw(tdc(x(i,j,4)))
                else
                   e = esati(tdc(x(i,j,4)))
                end if
                forc_qxy_d(i,j) = 0.622_r8*e / (forc_pbotxy_d(i,j) - 0.378_r8*e)
             end if
          else                !using QBOT in kg/kg
             if (forc_txy_d(i,j) > SHR_CONST_TKFRZ) then
                e = esatw(tdc(forc_txy_d(i,j)))
             else
                e = esati(tdc(forc_txy_d(i,j)))
             end if
             qsat = 0.622_r8*e / (forc_pbotxy_d(i,j) - 0.378_r8*e)
             if (qsat < x(i,j,3)) then
                forc_qxy_d(i,j) = qsat
!                 write(6,*)'ATM warning: qsat < q!'
             else
                forc_qxy_d(i,j) = x(i,j,3)
             end if
          end if

          ! ZGCMXY
          if (nint(x(i,j,6)) == -1) then
             zgcmxy_d(i,j) = 30._r8
          else
             zgcmxy_d(i,j) = x(i,j,6)
          end if

          ! FORC_SOLSXY, FORC_SOLLXY, FORC_SOLSDXY, FORC_SOLLDXY

          if (nint(x(i,j,9))==-1.or.nint(x(i,j,10))==-1) then
             if (nint(x(i,j,8)) /= -1) then
                forc_solsxy_d(i,j)  = 0.7_r8 * (0.5_r8 * x(i,j,8))
                forc_sollxy_d(i,j)  = forc_solsxy_d(i,j)
                forc_solsdxy_d(i,j) = 0.3_r8 * (0.5_r8 * x(i,j,8))
                forc_solldxy_d(i,j) = forc_solsdxy_d(i,j)
             else
                write(6,*)'ATM error: neither FSDSdir/dif nor'
                write(6,*)'       FSDS have been read in by atmrd'
                atmread_err = .true.
             end if
          else
             forc_solsxy_d(i,j)  = 0.5_r8 * x(i,j,9)
             forc_sollxy_d(i,j)  = forc_solsxy_d(i,j)
             forc_solsdxy_d(i,j) = 0.5_r8 * x(i,j,10)
             forc_solldxy_d(i,j) = forc_solsdxy_d(i,j)
          end if

          ! PRCXY, PRLXY

          if (nint(x(i,j,13))==-1.or.nint(x(i,j,14))==-1) then
             if (nint(x(i,j,12)).ne.-1) then
                prcxy_d(i,j) = 0.1_r8 * x(i,j,12)
                prlxy_d(i,j) = 0.9_r8 * x(i,j,12)
             else
                write(6,*)'ATM error: neither PRECC/L nor PRECT'
                write(6,*)'           have been read in by atmrd'
                atmread_err = .true.
             end if
          else
             prcxy_d(i,j) = x(i,j,13)
             prlxy_d(i,j) = x(i,j,14)
          end if

          ! FLWDSXY

          if (nint(x(i,j,11)) == -1) then
             e = forc_psrfxy_d(i,j) * forc_qxy_d(i,j) / (0.622_r8 + 0.378_r8 * forc_qxy_d(i,j))
             ea = 0.70_r8 + 5.95e-05_r8 * 0.01_r8*e * exp(1500.0_r8/forc_txy_d(i,j))
             flwdsxy_d(i,j) = ea * sb * forc_txy_d(i,j)**4
          else
             flwdsxy_d(i,j) = x(i,j,11)
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

  end subroutine atm_readdata

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: interpa2si
!
! !INTERFACE:
  subroutine interpa2si()
!
! !DESCRIPTION:
! Initialize variables for atm->land model surface interp
!
! !USES:
    use areaMod  , only : areaini,gridmap_checkmap
    use domainMod, only : adomain

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
    integer :: i,j                       !indices
    integer :: ier                       !error status
    integer :: atmlon,atmlat             !size of adomain
    real(r8), allocatable :: mask_d(:,:) !dummy field: atm grid mask
    real(r8), allocatable :: mask_a(:,:) !dummy field: land model grid mask
!------------------------------------------------------------------------

    atmlon = adomain%ni
    atmlat = adomain%nj

    ! Dynamically allocate memory

    allocate (mask_d(datlon,datlat),mask_a(atmlon,atmlat), stat=ier)
    if (ier /= 0) then
       write (6,*) 'interpa2si(): allocation error'
       call endrun
    end if

    if ( masterproc )then
       write (6,*) 'Attempting to initialize atm->land model grid interpolation .....'
       write (6,*) 'Initializing atm -> srf interpolation .....'
    end if

    ! --------------------------------------------------------------------
    ! Map from atmosphere grid to surface grid
    ! --------------------------------------------------------------------

    ! [mask_d] = 1 means all grid cells on atm grid, regardless of whether
    ! land or ocean, will contribute to surface grid.

    do j = 1, datlat
       do i = 1, datlon
          mask_d(i,j) = 1._r8
       end do
    end do

    ! [mask_a] = 1 means all the surface grid is land. Used as dummy
    ! variable so code will not abort with false, non-valid error check

    do j = 1, atmlat
       do i = 1, atmlon
          mask_a(i,j) = 1._r8
       end do
    end do

    ! For each surface grid cell: get lat [jovr_a2s] and lon [iovr_a2s] indices
    ! and weights [wovr_a2s] of overlapping atm grid cells

    call areaini (ddomain, adomain, gridmap_d2a, &
                  fracin=mask_d, fracout=mask_a )

    call gridmap_checkmap(gridmap_d2a)

    deallocate (mask_d, mask_a)

    if ( masterproc )then
       write (6,*) 'Successfully made atm -> srf interpolation'
       write (6,*) 'Successfully initialized area-averaging interpolation'
       write (6,*)
    end if

  end subroutine interpa2si

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
    use areaMod
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
    integer  :: i,j                  !longitude,latitude loop indices
    real(r8) :: forc_u(datlon,datlat)  !dummy wind (u)
    real(r8) :: forc_v(datlon,datlat)  !dummy wind (v)
!------------------------------------------------------------------------

    ! area-average absolute value of winds (i.e., regardless of
    ! direction) since land model cares about magnitude not direction.
    ! then need to adjust resultant stresses for direction of wind.

!$OMP PARALLEL DO PRIVATE (j,i)
#if !defined (USE_OMP)
!CSD$ PARALLEL DO PRIVATE (j,i)
#endif
    do j = 1, datlat
       do i = 1, datlon
          forc_u(i,j) = abs(forc_uxy_d(i,j))
          forc_v(i,j) = abs(forc_vxy_d(i,j))
       end do
    end do
#if !defined (USE_OMP)
!CSD$ END PARALLEL DO
#endif
!$OMP END PARALLEL DO

    call areaave ( forc_txy_d    , forc_txy_a    , gridmap_d2a )
    call areaave ( forc_u        , forc_uxy_a    , gridmap_d2a )
    call areaave ( forc_v        , forc_vxy_a    , gridmap_d2a )
    call areaave ( forc_qxy_d    , forc_qxy_a    , gridmap_d2a )
    call areaave ( zgcmxy_d      , zgcmxy_a      , gridmap_d2a )
    call areaave ( prcxy_d       , prcxy_a       , gridmap_d2a )
    call areaave ( prlxy_d       , prlxy_a       , gridmap_d2a )
    call areaave ( flwdsxy_d     , flwdsxy_a     , gridmap_d2a )
    call areaave ( forc_solsxy_d , forc_solsxy_a , gridmap_d2a )
    call areaave ( forc_sollxy_d , forc_sollxy_a , gridmap_d2a )
    call areaave ( forc_solsdxy_d, forc_solsdxy_a, gridmap_d2a )
    call areaave ( forc_solldxy_d, forc_solldxy_a, gridmap_d2a )
    call areaave ( forc_pbotxy_d , forc_pbotxy_a , gridmap_d2a )
    call areaave ( forc_psrfxy_d , forc_psrfxy_a , gridmap_d2a )

  end subroutine interpa2s

#endif

end module atmdrvMod

