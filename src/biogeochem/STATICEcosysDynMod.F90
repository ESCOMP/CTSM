#include <misc.h>
#include <preproc.h>
#if ( defined SCAM )
#include <max.h>
#endif

module STATICEcosysdynMOD

#if (!defined DGVM)

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: STATICEcosysDynMod
!
! !DESCRIPTION:
! Static Ecosystem dynamics: phenology, vegetation.
!
! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
  use abortutils,   only: endrun
#if ( defined SCAM )
  use scamMod, only :initlonidx,initlatidx
#endif
!
! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: EcosystemDyn       ! Ecosystem dynamics: phenology, vegetation
  public :: EcosystemDynini    ! Dynamically allocate memory
  public :: interpMonthlyVeg   ! interpolate monthly vegetation data
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!
! !PRIVATE MEMBER FUNCTIONS:
  private :: readMonthlyVegetation   ! read monthly vegetation data for two months
!
! PRIVATE TYPES:
  integer , private :: InterpMonths1         ! saved month index
  real(r8), private :: timwt(2)              ! time weights for month 1 and month 2
  real(r8), private, allocatable :: mlai1(:) ! lai for interpolation (month 1)
  real(r8), private, allocatable :: mlai2(:) ! lai for interpolation (month 2)
  real(r8), private, allocatable :: msai1(:) ! sai for interpolation (month 1)
  real(r8), private, allocatable :: msai2(:) ! sai for interpolation (month 2)
  real(r8), private, allocatable :: mhvt1(:) ! top vegetation height for interpolation (month 1)
  real(r8), private, allocatable :: mhvt2(:) ! top vegetation height for interpolation (month 2)
  real(r8), private, allocatable :: mhvb1(:) ! bottom vegetation height for interpolation(month 1)
  real(r8), private, allocatable :: mhvb2(:) ! bottom vegetation height for interpolation(month 2)
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: EcosystemDynini
!
! !INTERFACE:
  subroutine EcosystemDynini ()
!
! !DESCRIPTION:
! Dynamically allocate memory and set to signaling NaN.
!
! !USES:
    use nanMod
    use decompMod, only : get_proc_global
!
! !ARGUMENTS:
    implicit none
!
! !REVISION HISTORY:
!
!EOP
!
! LOCAL VARIABLES:
    integer :: ier    ! error code
    integer :: numg   ! total number of gridcells across all processors
    integer :: numl   ! total number of landunits across all processors
    integer :: numc   ! total number of columns across all processors
    integer :: nump   ! total number of pfts across all processors
!-----------------------------------------------------------------------

    InterpMonths1 = -999  ! saved month index
    call get_proc_global(numg, numl, numc, nump) 

    ier = 0
    if(.not.allocated(mlai1))allocate (mlai1(nump), mlai2(nump), &
              msai1(nump), msai2(nump), &
              mhvt1(nump), mhvt2(nump), &
              mhvb1(nump), mhvb2(nump), stat=ier)
    if (ier /= 0) then
       write (6,*) 'EcosystemDynini allocation error'
       call endrun
    end if

    mlai1(:) = nan
    mlai2(:) = nan
    msai1(:) = nan
    msai2(:) = nan
    mhvt1(:) = nan
    mhvt2(:) = nan
    mhvb1(:) = nan
    mhvb2(:) = nan

  end subroutine EcosystemDynini

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: EcosystemDyn
!
! !INTERFACE:
  subroutine EcosystemDyn(lbp, ubp, num_nolakep, filter_nolakep, doalb)
!
! !DESCRIPTION:
! Ecosystem dynamics: phenology, vegetation
! Calculates leaf areas (tlai, elai),  stem areas (tsai, esai) and
! height (htop).
!
! !USES:
    use clmtype
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: lbp, ubp                    ! pft bounds
    integer, intent(in) :: num_nolakep                 ! number of column non-lake points in pft filter
    integer, intent(in) :: filter_nolakep(ubp-lbp+1)   ! pft filter for non-lake points
    logical, intent(in) :: doalb                       ! true = surface albedo calculation time step
!
! !CALLED FROM:
!
! !REVISION HISTORY:
! Author: Gordon Bonan
! 2/1/02, Peter Thornton: Migrated to new data structure.
!
! !LOCAL VARIABLES:
!
! local pointers to implicit in arguments
!
    integer , pointer :: pcolumn(:)  ! column index associated with each pft
    real(r8), pointer :: snowdp(:)   ! snow height (m)
#if (defined CASA)
    real(r8), pointer :: plai(:)     ! Prognostic lai
    integer , pointer :: ivt(:)      ! pft vegetation type
#endif
!
! local pointers to implicit out arguments
!
    real(r8), pointer :: tlai(:)     ! one-sided leaf area index, no burying by snow
    real(r8), pointer :: tsai(:)     ! one-sided stem area index, no burying by snow
    real(r8), pointer :: htop(:)     ! canopy top (m)
    real(r8), pointer :: hbot(:)     ! canopy bottom (m)
    real(r8), pointer :: elai(:)     ! one-sided leaf area index with burying by snow
    real(r8), pointer :: esai(:)     ! one-sided stem area index with burying by snow
    integer , pointer :: frac_veg_nosno_alb(:) ! frac of vegetation not covered by snow [-]
!
!EOP
!
! !OTHER LOCAL VARIABLES:
!
    integer  :: fp,p,c   ! indices
    real(r8) :: ol       ! thickness of canopy layer covered by snow (m)
    real(r8) :: fb       ! fraction of canopy layer covered by snow
!-----------------------------------------------------------------------

    if (doalb) then

       ! Assign local pointers to derived type scalar members (column-level)

       snowdp  => clm3%g%l%c%cps%snowdp

       ! Assign local pointers to derived type scalar members (pftlevel)

       pcolumn => clm3%g%l%c%p%column
       tlai    => clm3%g%l%c%p%pps%tlai
       tsai    => clm3%g%l%c%p%pps%tsai
       elai    => clm3%g%l%c%p%pps%elai
       esai    => clm3%g%l%c%p%pps%esai
       htop    => clm3%g%l%c%p%pps%htop
       hbot    => clm3%g%l%c%p%pps%hbot
       frac_veg_nosno_alb => clm3%g%l%c%p%pps%frac_veg_nosno_alb
#if (defined CASA)
       plai    => clm3%g%l%c%p%pps%plai
       ivt     => clm3%g%l%c%p%itype
#endif

!dir$ concurrent
!cdir nodep
       do fp = 1, num_nolakep
          p = filter_nolakep(fp)
          c = pcolumn(p)

          ! need to update elai and esai only every albedo time step so do not
          ! have any inconsistency in lai and sai between SurfaceAlbedo calls (i.e.,
          ! if albedos are not done every time step).
          ! leaf phenology
          ! Set leaf and stem areas based on day of year
          ! Interpolate leaf area index, stem area index, and vegetation heights
          ! between two monthly
          ! The weights below (timwt(1) and timwt(2)) were obtained by a call to
          ! routine InterpMonthlyVeg in subroutine NCARlsm.
          !                 Field   Monthly Values
          !                -------------------------
          ! leaf area index LAI  <- mlai1 and mlai2
          ! leaf area index SAI  <- msai1 and msai2
          ! top height      HTOP <- mhvt1 and mhvt2
          ! bottom height   HBOT <- mhvb1 and mhvb2

          tlai(p) = timwt(1)*mlai1(p) + timwt(2)*mlai2(p)
          tsai(p) = timwt(1)*msai1(p) + timwt(2)*msai2(p)
          htop(p) = timwt(1)*mhvt1(p) + timwt(2)*mhvt2(p)
          hbot(p) = timwt(1)*mhvb1(p) + timwt(2)*mhvb2(p)

#if (defined CASA)
! use PLAI from CASA
!! can also choose to use default tlai instead of plai - need option
! 03/03/11: don't reset lai of crops - use default lsm tlai/elai

          if (ivt(p) > 0 .and. ivt(p) < 15) tlai(p) = plai(p)
#endif

          ! adjust lai and sai for burying by snow. if exposed lai and sai
          ! are less than 0.05, set equal to zero to prevent numerical
          ! problems associated with very small lai and sai.

          ol = min( max(snowdp(c)-hbot(p), 0._r8), htop(p)-hbot(p))
          fb = 1._r8 - ol / max(1.e-06_r8, htop(p)-hbot(p))
          elai(p) = max(tlai(p)*fb, 0.0_r8)
          esai(p) = max(tsai(p)*fb, 0.0_r8)
          if (elai(p) < 0.05_r8) elai(p) = 0._r8
          if (esai(p) < 0.05_r8) esai(p) = 0._r8

          ! Fraction of vegetation free of snow

          if ((elai(p) + esai(p)) >= 0.05_r8) then
             frac_veg_nosno_alb(p) = 1
          else
             frac_veg_nosno_alb(p) = 0
          end if

       end do ! end of pft loop

    end if  !end of if-doalb block

  end subroutine EcosystemDyn

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: interpMonthlyVeg
!
! !INTERFACE:
  subroutine interpMonthlyVeg ()
!
! !DESCRIPTION:
! Determine if 2 new months of data are to be read.
!
! !USES:
    use clm_varctl  , only : fsurdat
#ifdef COUP_CAM
    use time_manager, only : get_curr_date, get_step_size, get_perp_date, is_perpetual
#else
    use time_manager, only : get_curr_date, get_step_size
#endif
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
    integer :: kyr         ! year (0, ...) for nstep+1
    integer :: kmo         ! month (1, ..., 12)
    integer :: kda         ! day of month (1, ..., 31)
    integer :: ksec        ! seconds into current date for nstep+1
    real(r8):: dtime       ! land model time step (sec)
    real(r8):: t           ! a fraction: kda/ndaypm
    integer :: it(2)       ! month 1 and month 2 (step 1)
    integer :: months(2)   ! months to be interpolated (1 to 12)
    integer, dimension(12) :: ndaypm= &
         (/31,28,31,30,31,30,31,31,30,31,30,31/) !days per month
!-----------------------------------------------------------------------

    dtime = get_step_size()

#ifdef COUP_CAM
    if ( is_perpetual() ) then
       call get_perp_date(kyr, kmo, kda, ksec, offset=int(dtime))
    else
       call get_curr_date(kyr, kmo, kda, ksec, offset=int(dtime))
    end if
#else
    call get_curr_date(kyr, kmo, kda, ksec, offset=int(dtime))
#endif

    t = (kda-0.5_r8) / ndaypm(kmo)
    it(1) = t + 0.5_r8
    it(2) = it(1) + 1
    months(1) = kmo + it(1) - 1
    months(2) = kmo + it(2) - 1
    if (months(1) <  1) months(1) = 12
    if (months(2) > 12) months(2) = 1
    timwt(1) = (it(1)+0.5_r8) - t
    timwt(2) = 1._r8-timwt(1)

    if (InterpMonths1 /= months(1)) then
       call readMonthlyVegetation (fsurdat, kmo, kda, months)
       InterpMonths1 = months(1)
    end if

  end subroutine interpMonthlyVeg

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: readMonthlyVegetation
!
! !INTERFACE:
  subroutine readMonthlyVegetation (fveg, kmo, kda, months)
!
! !DESCRIPTION:
! Read monthly vegetation data for two consec. months.
!
! !USES:
    use clmtype
    use decompMod   , only : get_proc_global
    use clm_varpar  , only : lsmlon, lsmlat, maxpatch_pft, maxpatch, npatch_crop, numpft
    use clm_varsur  , only : all_pfts_on_srfdat
    use pftvarcon   , only : noveg
    use fileutils   , only : getfil
#if (defined SPMD)
    use spmdMod     , only : masterproc, mpicom, MPI_REAL8
#else
    use spmdMod     , only : masterproc
#endif
    use time_manager, only : get_nstep
    use ncdio       , only : check_ret
#if ( defined SCAM )
    use getnetcdfdata
#endif

!
! !ARGUMENTS:
    implicit none
    include 'netcdf.inc'
    character(len=*), intent(in) :: fveg  ! file with monthly vegetation data
    integer, intent(in) :: kmo            ! month (1, ..., 12)
    integer, intent(in) :: kda            ! day of month (1, ..., 31)
    integer, intent(in) :: months(2)      ! months to be interpolated (1 to 12)
!
! !REVISION HISTORY:
! Created by Sam Levis
!
!EOP
!
! LOCAL VARIABLES:
    character(len=256) :: locfn           ! local file name
    integer :: i,j,k,l,m,p,ivt            ! indices
    integer :: ncid,dimid,varid           ! input netCDF id's
    integer :: beg4d(4),len4d(4)          ! netCDF variable edges
    integer :: ntim                       ! number of input data time samples
    integer :: nlon_i                     ! number of input data longitudes
    integer :: nlat_i                     ! number of input data latitudes
    integer :: npft_i                     ! number of input data pft types
    integer :: numg                       ! total number of gridcells across all processors
    integer :: numl                       ! total number of landunits across all processors
    integer :: numc                       ! total number of columns across all processors
    integer :: nump                       ! total number of pfts across all processors
    integer :: ier                        ! error code
    real(r8), allocatable :: mlai(:,:,:)  ! lai read from input files
    real(r8), allocatable :: msai(:,:,:)  ! sai read from input files
    real(r8), allocatable :: mhgtt(:,:,:) ! top vegetation height
    real(r8), allocatable :: mhgtb(:,:,:) ! bottom vegetation height
#if ( defined SCAM )
    real(r8), allocatable :: coldata(:) ! temporary for getncdata calls
#endif
    character(len=32) :: subname = 'readMonthlyVegetation'
!-----------------------------------------------------------------------

    if (all_pfts_on_srfdat) then
       allocate(mlai(lsmlon,lsmlat,0:numpft), &
                msai(lsmlon,lsmlat,0:numpft), &
                mhgtt(lsmlon,lsmlat,0:numpft), &
                mhgtb(lsmlon,lsmlat,0:numpft), stat=ier)
    else
       allocate(mlai(lsmlon,lsmlat,maxpatch_pft), &
                msai(lsmlon,lsmlat,maxpatch_pft), &
                mhgtt(lsmlon,lsmlat,maxpatch_pft), &
                mhgtb(lsmlon,lsmlat,maxpatch_pft), stat=ier)
    end if
    if (ier /= 0) then
       write(6,*)subname, 'allocation error '; call endrun()
    end if
#if ( defined SCAM )
    if (all_pfts_on_srfdat) then
       allocate(coldata(0:numpft), stat=ier)
    else
       allocate(coldata(maxpatch_pft), stat=ier)
    end if
    if (ier /= 0) then
       write(6,*)subname, 'allocation error '; call endrun()
    end if
#endif
    ! Determine necessary indices

    call get_proc_global(numg, numl, numc, nump)

    ! ----------------------------------------------------------------------
    ! Open monthly vegetation file
    ! Read data and convert from [lsmlon] x [lsmlat] grid to patch data
    ! ----------------------------------------------------------------------

    if (masterproc) then

       write (6,*) 'Attempting to read monthly vegetation data .....'
       write (6,*) 'nstep = ',get_nstep(),' month = ',kmo,' day = ',kda

       call getfil(fveg, locfn, 0)
       call check_ret(nf_open(locfn, 0, ncid), subname)

       call check_ret(nf_inq_dimid (ncid, 'lsmlon', dimid), subname)
       call check_ret(nf_inq_dimlen(ncid, dimid, nlon_i), subname)
#if ( !defined SCAM)
       if (nlon_i /= lsmlon) then
          write(6,*)subname,' parameter lsmlon= ',lsmlon,'does not equal input nlat_i= ',nlon_i
          call endrun()
       end if
#endif

       call check_ret(nf_inq_dimid(ncid, 'lsmlat', dimid), subname)
       call check_ret(nf_inq_dimlen(ncid, dimid, nlat_i), subname)

#if ( !defined SCAM )
       if (nlat_i /= lsmlat) then
          write(6,*)subname,' parameter lsmlat= ',lsmlat,'does not equal input nlat_i= ',nlat_i
          call endrun()
       end if
#endif

       call check_ret(nf_inq_dimid(ncid, 'lsmpft', dimid), subname)
       call check_ret(nf_inq_dimlen(ncid, dimid, npft_i), subname)

#if ( !defined SCAM )
       if (all_pfts_on_srfdat) then
          if (npft_i /= numpft+1) then
             write(6,*)subname,' parameter numpft+1 = ',numpft+1,'does not equal input npft_i= ',npft_i
             call endrun()
          end if
       else
          if (npft_i /= maxpatch_pft) then
             write(6,*)subname,' parameter maxpatch_pft = ',maxpatch_pft,'does not equal input npft_i= ',npft_i
             call endrun()
          end if
       end if
#endif

       call check_ret(nf_inq_dimid(ncid, 'time', dimid), subname)
       call check_ret(nf_inq_dimlen(ncid, dimid, ntim), subname)

       do k=1,2   !loop over months and read vegetated data

          beg4d(1) = 1         ; len4d(1) = nlon_i
          beg4d(2) = 1         ; len4d(2) = nlat_i
          beg4d(3) = 1         ; len4d(3) = npft_i
          beg4d(4) = months(k) ; len4d(4) = 1

#if ( defined SCAM )
          call getncdata (ncid, initLatIdx, initLonIdx, months(k),'MONTHLY_LAI', coldata, IER)
          mlai(1,1,:)=coldata
          
          call getncdata (ncid, initLatIdx, initLonIdx, months(k),'MONTHLY_SAI', coldata, IER)
          msai(1,1,:)=coldata

          call getncdata (ncid, initLatIdx, initLonIdx, months(k),'MONTHLY_HEIGHT_TOP', coldata, IER)
          mhgtt(1,1,:) = coldata

          call getncdata (ncid, initLatIdx, initLonIdx, months(k),'MONTHLY_HEIGHT_BOT', coldata, IER)
          mhgtb(1,1,:) = coldata
#else
          call check_ret(nf_inq_varid(ncid, 'MONTHLY_LAI', varid), subname)
          call check_ret(nf_get_vara_double(ncid, varid, beg4d, len4d, mlai), subname)

          call check_ret(nf_inq_varid(ncid, 'MONTHLY_SAI', varid), subname)
          call check_ret(nf_get_vara_double(ncid, varid, beg4d, len4d, msai), subname)

          call check_ret(nf_inq_varid(ncid, 'MONTHLY_HEIGHT_TOP', varid), subname)
          call check_ret(nf_get_vara_double(ncid, varid, beg4d, len4d, mhgtt), subname)

          call check_ret(nf_inq_varid(ncid, 'MONTHLY_HEIGHT_BOT', varid), subname)
          call check_ret(nf_get_vara_double(ncid, varid, beg4d, len4d, mhgtb), subname)
#endif

          ! store data directly in clmtype structure
          ! only vegetated pfts have nonzero values

          do p = 1, nump

             i = clm3%g%l%c%p%ixy(p)
             j = clm3%g%l%c%p%jxy(p)

             ! Assign lai/sai/hgtt/hgtb to the top [maxpatch_pft] pfts
             ! as determined in subroutine surfrd

             if (all_pfts_on_srfdat) then

                ivt = clm3%g%l%c%p%itype(p)
                if (ivt /= noveg) then     ! vegetated pft
                   do l = 0, numpft
                      if (l == ivt) then
                         if (k == 1) then
                            mlai1(p) = mlai(i,j,l)
                            msai1(p) = msai(i,j,l)
                            mhvt1(p) = mhgtt(i,j,l)
                            mhvb1(p) = mhgtb(i,j,l)
                         else !if (k == 2)
                            mlai2(p) = mlai(i,j,l)
                            msai2(p) = msai(i,j,l)
                            mhvt2(p) = mhgtt(i,j,l)
                            mhvb2(p) = mhgtb(i,j,l)
                         end if
                      end if
                   end do
                else                        ! non-vegetated pft
                   if (k == 1) then
                      mlai1(p) = 0._r8
                      msai1(p) = 0._r8
                      mhvt1(p) = 0._r8
                      mhvb1(p) = 0._r8
                   else   !if (k == 2)
                      mlai2(p) = 0._r8
                      msai2(p) = 0._r8
                      mhvt2(p) = 0._r8
                      mhvb2(p) = 0._r8
                   end if
                end if

             else

                m = clm3%g%l%c%p%mxy(p)
                if (m <= maxpatch_pft) then ! vegetated pft
                   if (k == 1) then
                      mlai1(p) = mlai(i,j,m)
                      msai1(p) = msai(i,j,m)
                      mhvt1(p) = mhgtt(i,j,m)
                      mhvb1(p) = mhgtb(i,j,m)
                   else !if (k == 2)
                      mlai2(p) = mlai(i,j,m)
                      msai2(p) = msai(i,j,m)
                      mhvt2(p) = mhgtt(i,j,m)
                      mhvb2(p) = mhgtb(i,j,m)
                   end if
                else                        ! non-vegetated pft
                   if (k == 1) then
                      mlai1(p) = 0._r8
                      msai1(p) = 0._r8
                      mhvt1(p) = 0._r8
                      mhvb1(p) = 0._r8
                   else !if (k == 2)
                      mlai2(p) = 0._r8
                      msai2(p) = 0._r8
                      mhvt2(p) = 0._r8
                      mhvb2(p) = 0._r8
                   end if
                end if

             end if

          end do   ! end of loop over pfts

       end do   ! end of loop over months

       call check_ret(nf_close(ncid), subname)

       write (6,*) 'Successfully read monthly vegetation data for'
       write (6,*) 'month ', months(1), ' and month ', months(2)
       write (6,*)

    end if ! end of if-masterproc if block

#if ( defined SPMD )
    ! pass surface data to all processors
    call mpi_bcast (mlai1, size(mlai1), MPI_REAL8, 0, mpicom, ier)
    call mpi_bcast (mlai2, size(mlai2), MPI_REAL8, 0, mpicom, ier)
    call mpi_bcast (msai1, size(msai1), MPI_REAL8, 0, mpicom, ier)
    call mpi_bcast (msai2, size(msai2), MPI_REAL8, 0, mpicom, ier)
    call mpi_bcast (mhvt1, size(mhvt1), MPI_REAL8, 0, mpicom, ier)
    call mpi_bcast (mhvt2, size(mhvt2), MPI_REAL8, 0, mpicom, ier)
    call mpi_bcast (mhvb1, size(mhvb1), MPI_REAL8, 0, mpicom, ier)
    call mpi_bcast (mhvb2, size(mhvb2), MPI_REAL8, 0, mpicom, ier)
#endif

    deallocate(mlai, msai, mhgtt, mhgtb)

  end subroutine readMonthlyVegetation

#endif

end module STATICEcosysDynMod
