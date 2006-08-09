#include <misc.h>
#include <preproc.h>

module RtmMod

#if (defined RTM)

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: RtmMod
!
! !DESCRIPTION:
! River Routing Model (U. of Texas River Transport
! Model)~\cite{Branstetter:2001}
!
! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use spmdMod     , only : masterproc
  use clm_varpar  , only : lsmlon, lsmlat, rtmlon, rtmlat
  use shr_sys_mod , only : shr_sys_flush
  use domainMod   , only : domain_type, domain_init
  use abortutils  , only : endrun
!
! !PUBLIC TYPES:
  implicit none
  save
  integer , parameter, public :: rtmloni = 1        ! per-proc beg lon index
  integer , parameter, public :: rtmlonf = rtmlon   ! per-proc end lon index
  integer , parameter, public :: rtmlati = 1        ! per-proc beg lat index
  integer , parameter, public :: rtmlatf = rtmlat   ! per-proc end lat index
  type (domain_type) , public :: rdomain            ! rtm grid 
  real(r8), allocatable, public :: volr(:,:)        ! cell h2o volume (m^3)
!
! !PUBLIC MEMBER FUNCTIONS:
  public Rtmini        ! Initialize RTM grid and land mask
  public Rtmfluxini    ! Initialize RTM fluxout
  public Rtmriverflux  ! Interface with RTM river routing model
  public RtmRest       ! Read/write RTM restart data (netcdf)
!
! !REVISION HISTORY:
! Author: Sam Levis
!
!EOP
!
! !PRIVATE MEMBER FUNCTIONS:
!
  private RtmUpdateInput   ! Update rtm inputs
  private Rtm           ! River routing model (based on U. Texas code)
!
! !PRIVATE TYPES:
  private
!
! RTM input
!
  real(r8), pointer :: rtmin_ave(:)          ! RTM averaging buffer for runoff
  real(r8), pointer :: rtmin_glob(:)         ! RTM global input
  real(r8), pointer :: rtmin(:)              ! RTM local input
  integer  :: ncount_rtm                     ! RTM time averaging = number of time samples to average over
  real(r8) :: delt_rtm                       ! RTM time step
!
! RTM 1/2 degree resolution variables
!
  real(r8), parameter   :: effvel = 0.35_r8  ! RTM effective velocity (m/s)
  integer , allocatable :: rdirc(:,:)        ! RTM input - rtm river flow direction (0-8)
  real(r8), allocatable :: ddist(:,:)        ! RTM input - downstream dist (m)
  real(r8), allocatable :: rivarea(:,:)      ! RTM input - cell area (m^2)
  real(r8), allocatable :: totrunin_r(:,:)   ! RTM input - surface roff (mm/s)
  real(r8), allocatable :: sfluxin(:,:)      ! RTM input - cell in h2o flux (m3/s)
  real(r8), allocatable :: fluxout(:,:)      ! RTM input/output - water flux out of cell (m^3/s)
  real(r8), allocatable :: flxlnd_r(:,:)     ! RTM output - river flux (m**3/s)
  real(r8), allocatable :: flxocn_r(:,:)     ! RTM output - river flux to the ocean (m**3/s)
  real(r8), allocatable :: dvolrdt_r(:,:)    ! RTM output - change in storage (mm/s)
  real(r8), allocatable :: dvolrdt_lnd_r(:,:)! RTM output - change in storage (mm/s)
  real(r8), allocatable :: dvolrdt_ocn_r(:,:)! RTM output - change in storage (mm/s)
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Rtmini
!
! !INTERFACE:
  subroutine Rtmini
!
! !DESCRIPTION:
! Initialize RTM grid and land mask.
!
! !USES:
    use shr_const_mod, only : SHR_CONST_PI
    use domainMod    , only : ldomain
    use areaMod      , only : celledge, cellarea
    use clm_varctl   , only : frivinp_rtm
    use clm_varcon   , only : re
    use decompMod    , only : get_proc_bounds, get_proc_global
    use clm_time_manager , only : get_curr_date
    use RunoffMod    , only : set_proc_rof_bounds, set_roflnd, set_rofocn
!
! !ARGUMENTS:
    implicit none
!
! !CALLED FROM:
! subroutine initialize in module initializeMod
!
! !REVISION HISTORY:
! Author: Sam Levis
!
!EOP
!
! !LOCAL VARIABLES:
    real(r8), dimension(4) :: rtmedge = (/ 90._r8, 180._r8, -90._r8, -180._r8 /)  !N,E,S,W edges of rtm grid
    integer  :: ioff(0:8) = (/0,0,1,1,1,0,-1,-1,-1/) !calc dist as in hydra
    integer  :: joff(0:8) = (/0,1,1,0,-1,-1,-1,0,1/) !of grid cell down stream
    integer  :: i,j,k,n,g,n2                  ! loop indices
    integer  :: i2,j2                         ! downstream i and j
    integer  :: is,js                         ! land model grid indices
    integer  :: ir,jr                         ! rtm grid indices
    real(r8) :: deg2rad                       ! pi/180
    real(r8) :: dx                            ! lon dist. betn grid cells (m)
    real(r8) :: dy                            ! lat dist. betn grid cells (m)
    real(r8) :: tempg(rtmlon,rtmlat)          ! temporary buffer
    integer  :: tempgp(0:rtmlon+1,0:rtmlat+1) ! temporary buffer
    integer  :: masktmp_r(rtmlon,rtmlat)      ! dummy mask
    integer  :: ier                           ! error code
    integer  :: pid                           ! processor id
    integer  :: mon                           ! month (1, ..., 12)
    integer  :: day                           ! day of month (1, ..., 31)
    integer  :: ncsec                         ! seconds of current date
    integer  :: numg                          ! tot num of gridcells on all pes
    integer  :: numl                          ! tot num of landunits on all pes
    integer  :: numc                          ! tot num of columns on all pes
    integer  :: nump                          ! tot num of pfts on all pes
    real(r8),allocatable :: ltmp(:,:),rtmp(:,:) ! temporaries
    character(len=16), dimension(50) :: river_name
    character(len=30), dimension(50) :: rivstat_name
    real(r8)         , dimension(50) :: rivstat_lon
    real(r8)         , dimension(50) :: rivstat_lat
    integer  :: nroflnd
    integer  :: nrofocn
!-----------------------------------------------------------------------

    ! Allocate rtm grid variables
    call domain_init(rdomain,rtmlon,rtmlat)

    ! Allocate rtm flux variables

    allocate (volr(rtmloni:rtmlonf,rtmlati:rtmlatf),&
              rdirc(rtmloni-1:rtmlonf+1,rtmlati-1:rtmlatf+1), &
              fluxout(rtmloni-1:rtmlonf+1,rtmlati-1:rtmlatf+1), &
              ddist(rtmloni:rtmlonf,rtmlati:rtmlatf), &
              rivarea(rtmloni:rtmlonf,rtmlati:rtmlatf), stat=ier)
    if (ier /= 0) then
       write(6,*)'Rtmgridini: Allocation error for ',&
            'rdirec, fluxout, ddist, rivarea'
       call endrun
    end if

    ! Allocate inputs and outputs to rtm at 1/2 degree resolution

    allocate (totrunin_r(rtmloni:rtmlonf,rtmlati:rtmlatf), &
              flxlnd_r(rtmloni:rtmlonf,rtmlati:rtmlatf), &
              flxocn_r(rtmloni:rtmlonf,rtmlati:rtmlatf), &
              dvolrdt_lnd_r(rtmloni:rtmlonf,rtmlati:rtmlatf), &
              dvolrdt_ocn_r(rtmloni:rtmlonf,rtmlati:rtmlatf), &
              dvolrdt_r(rtmloni:rtmlonf,rtmlati:rtmlatf), &
              sfluxin(rtmloni:rtmlonf,rtmlati:rtmlatf),  stat=ier)
    if (ier /= 0) then
       write(6,*)'Rtmgridini: Allocation error for ',&
            'totrunin_r, flxlnd_r, flxocn_r, dvolrdt_r, dvolrdt_lnd_r, dvolrdt_ocn_r, sfluxin'
       call endrun
    end if

    ! Useful constants and initial values

    deg2rad = SHR_CONST_PI / 180._r8
    volr(:,:) = 0._r8
    fluxout(:,:) = 0._r8
    flxocn_r(:,:) = 0._r8
    flxlnd_r(:,:) = 0._r8

    ! Open and read input data (river direction file)
    ! rtm operates from south to north and from the dateline
    ! River station data is currently not used by the model -
    ! it is only used by the diagnostic package
    ! If the river direction file is modified - the river station
    ! part must also be modified

    open (1,file=frivinp_rtm)
    do j = 1,rtmlat
       do i = 1,rtmlon
          n = (j-1)*rtmlon + i
          read(1,*) rdomain%latc(n),rdomain%lonc(n),tempg(i,j)
          tempgp(i,j) = nint(tempg(i,j))
       enddo
    enddo
    do n = 1,50
       read(1,10,iostat=ier) river_name(n), rivstat_lon(n), rivstat_lat(n), rivstat_name(n)
       if (ier /= 0) exit
10     format(1x,a16,f7.2,1x,f7.2,a30)
    end do
    close(1)

    if (masterproc) then
       write(6,*)'Columns in RTM = ',rtmlon
       write(6,*)'Rows in RTM    = ',rtmlat
       write(6,*)'read river direction data'
    end if

    ! Determine rtm mask, downstream distance and area
    ! The following assumes that there is no runoff 
    !   south of j=1 or north of j=rtmlat
    ! This is true for rdirc.05

    do i=1,rtmlon
       tempgp(i,0)        = 0._r8
       tempgp(i,rtmlat+1) = 0._r8
    enddo
    do j=0,rtmlat+1
       tempgp(0,j)        = tempgp(rtmlon,j)
       tempgp(rtmlon+1,j) = tempgp(1,j)
    enddo

    ! Determine rtm river flow direction (0-8)

    do j=rtmlati-1,rtmlatf+1
       do i=rtmloni-1,rtmlonf+1
          rdirc(i,j)=tempgp(i,j)
       enddo
    enddo

    ! Determine rtm ocn/land mask

    do j=rtmlati,rtmlatf
       do i=rtmloni,rtmlonf
          n = (j-1)*rtmlon + i
          if (rdirc(i,j) == 0) then
             rdomain%mask(n) = 0
          else
             rdomain%mask(n) = 1
          end if
       enddo
    enddo

    ! Determine RTM celledges and areas 

    call celledge (rdomain, &
                   rtmedge(1), rtmedge(2), rtmedge(3), rtmedge(4))

    call cellarea (rdomain, &
                   rtmedge(1), rtmedge(2), rtmedge(3), rtmedge(4))

    ! Determine downstream distance - instead of reading a distance file
    ! calculate the downstream distance

    do j=rtmlati,rtmlatf
       do i=rtmloni,rtmlonf
          n = (j-1)*rtmlon + i
          i2 = i + ioff(tempgp(i,j))
          j2 = j + joff(tempgp(i,j))
          if (i2 == 0) i2 = 2                 !avoids i2 out of bounds in the following
          if (i2 == rtmlon+1) i2 = rtmlon-1   !avoids i2 out of bounds in the following
          n2 = (j2-1)*rtmlon + i2
          dy = deg2rad * abs(rdomain%latc(n)-rdomain%latc(n2)) * re*1000._r8
          dx = deg2rad * abs(rdomain%lonc(n)-rdomain%lonc(n2)) * re*1000._r8 &
               *0.5_r8*(cos(rdomain%latc(n)*deg2rad)+cos(rdomain%latc(n2)*deg2rad))
          ddist(i,j) = sqrt(dx*dx + dy*dy)
          rivarea(i,j)=1.e6_r8 * rdomain%area(n)     !convert into m**2
       enddo
    enddo

    !---------------------------------------------------------------------
    ! Determine which ocean and land cells have runoff values
    !---------------------------------------------------------------------

    ! First loop over all ocean points and determine which are at the
    ! end of rivers by examining if any neighboring points are land and
    ! if that land neighbor points into this ocean point. Next loop over all
    ! ocean points and determine which overlap with at least one land cell.
    ! Allocate ocean runoff vector and indices and determine indices
    ! need to reset cpl runoff size to 0 and do the counting again because 
    ! need to first count to allocate vector and must now count to actually 
    ! determine indices

    allocate(ltmp(lsmlon,lsmlat),rtmp(rtmlon,rtmlat))
    do j = 1,lsmlat
    do i = 1,lsmlon
       n = (j-1)*lsmlon + i
       ltmp(i,j) = ldomain%frac(n)
    enddo
    enddo
!    call RtmMapClm2Rtm( ldomain%frac, rdomain%frac )
    call RtmMapClm2Rtm( ltmp, rtmp )
    do j = 1,rtmlat
    do i = 1,rtmlon
       n = (j-1)*rtmlon + i
       rdomain%frac(n) = rtmp(i,j)
    enddo
    enddo
    deallocate(ltmp,rtmp)

    nrofocn = 0
    masktmp_r(:,:) = 0
    do j=rtmlati,rtmlatf
       do i=rtmloni,rtmlonf
          n = (j-1)*rtmlon + i
          if (rdomain%mask(n) == 0) then
             if (rdirc(i  ,j-1)==1) masktmp_r(i,j) = 1
             if (rdirc(i-1,j-1)==2) masktmp_r(i,j) = 1
             if (rdirc(i-1,j  )==3) masktmp_r(i,j) = 1
             if (rdirc(i-1,j+1)==4) masktmp_r(i,j) = 1
             if (rdirc(i  ,j+1)==5) masktmp_r(i,j) = 1
             if (rdirc(i+1,j+1)==6) masktmp_r(i,j) = 1
             if (rdirc(i+1,j  )==7) masktmp_r(i,j) = 1
             if (rdirc(i+1,j-1)==8) masktmp_r(i,j) = 1
             if (masktmp_r(i,j) == 0 .and. rdomain%frac(n)>0._r8) masktmp_r(i,j) = 1
          endif
          if (masktmp_r(i,j) == 1) nrofocn = nrofocn +1
       enddo
    enddo

    call set_rofocn(rtmlon, rtmlat, masktmp_r, rdomain%area, nrofocn)

    ! Determine which land cells have runoff values

    nroflnd = 0
    masktmp_r(:,:) = 0
    do j=rtmlati,rtmlatf
       do i=rtmloni,rtmlonf
          n = (j-1)*rtmlon + i
          if (rdomain%mask(n) == 1) then
             masktmp_r(i,j) = 1
             nroflnd = nroflnd +1
          end if
       end do
    end do

    call set_roflnd(rtmlon, rtmlat, masktmp_r, rdomain%area, nroflnd)

    ! Determine per-processor runoff bounds

    call set_proc_rof_bounds()

    ! Allocate and initialize dynamic memory for rtm input

    call get_proc_global(numg, numl, numc, nump)
    allocate (rtmin(numg), rtmin_glob(numg), rtmin_ave(numg), stat=ier)
    if (ier /= 0) then
       write(6,*)'Rtmlandini: Allocation error for rtmin, rtmin_glob, rtmin_ave'
       call endrun
    end if
    rtmin(:)      = 0._r8
    rtmin_glob(:) = 0._r8
    rtmin_ave(:)  = 0._r8

  end subroutine Rtmini

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Rtmfluxini()
!
! !INTERFACE:
  subroutine Rtmfluxini()
!
! !DESCRIPTION:
! Initialize RTM fluxout for case of initial run when initial data is
! read in. For restart run, RTM fluxout is read from restart dataset.
!
! !USES:
    use clm_time_manager, only : get_step_size
    use clm_varctl  , only : rtm_nsteps
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
    integer :: i,j,n       !indices
    integer :: delt        !delt for rtm
!-----------------------------------------------------------------------

    delt = rtm_nsteps*get_step_size()
    do j = rtmlati,rtmlatf
       do i = rtmloni,rtmlonf
          n = (j-1)*rtmlon + i
          if (rdomain%mask(n)==1) then
             fluxout(i,j) = volr(i,j) * effvel/ddist(i,j)
             fluxout(i,j) = min(fluxout(i,j), volr(i,j) / delt)
          else
             fluxout(i,j) = 0._r8
          endif
       enddo
    enddo

  end subroutine Rtmfluxini

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Rtmriverflux
!
! !INTERFACE:
  subroutine Rtmriverflux()
!
! !DESCRIPTION:
! Interface with RTM river routing model.
!
! !USES:
    use RunoffMod   , only : UpdateRunoff
!
! !ARGUMENTS:
    implicit none
!
! !CALLED FROM:
! subroutine driver
!
! !REVISION HISTORY:
! Author: Sam Levis
!
!EOP
!
    logical  :: do_rtm                       ! true => perform rtm calculation
    real(r8),allocatable :: totruninxy(:,:)  ! surface runoff (mm H2O /s)
!-----------------------------------------------------------------------

    allocate(totruninxy(lsmlon,lsmlat))

    ! Determine RTM inputs on land model grid

    call RtmUpdateInput(do_rtm, totruninxy)

    if (do_rtm) then

       ! Map RTM inputs from land model grid to RTM grid (1/2 degree resolution)

       call RtmMapClm2Rtm( totruninxy, totrunin_r, landfrac_scale=.true. )

       ! Determine RTM runoff fluxes

       call t_startf('rtm_calc')
       call Rtm()
       call t_stopf('rtm_calc')

       ! Update coupler and history file info

       call t_startf('rtm_update')
       call UpdateRunoff(rtmloni, rtmlonf, rtmlati, rtmlatf, flxocn_r, flxlnd_r, &
                         dvolrdt_ocn_r, dvolrdt_lnd_r)
       call t_stopf('rtm_update')

    end if

    deallocate(totruninxy)

  end subroutine Rtmriverflux

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: RtmUpdateInput
!
! !INTERFACE:
  subroutine RtmUpdateInput(do_rtm, totruninxy)
!
! !DESCRIPTION:
! Update RTM inputs.
!
! !USES:
    use clmtype        , only : clm3,nameg
    use decompMod      , only : get_proc_bounds, get_proc_global, ldecomp
    use clm_varctl     , only : rtm_nsteps
    use clm_time_manager   , only : get_step_size, get_nstep
#if (defined SPMD)
    use spmdMod        , only : mpicom
    use spmdGathScatMod, only : scatter_data_from_master, gather_data_to_master, allgather_data
#endif
!
! !ARGUMENTS:
    implicit none
    logical , intent(out) :: do_rtm
    real(r8), intent(out) :: totruninxy(:,:) ! surface runoff (mm H2O /s)
!
! !CALLED FROM:
! subroutine driver
!
! !REVISION HISTORY:
! Author: Sam Levis
!
! !LOCAL VARIABLES:
!
! local pointers to implicit in arguments
!
    integer , pointer :: ixy(:)               ! gricell xy lon index
    integer , pointer :: jxy(:)               ! gricell xy lat index
    integer , pointer :: cgridcell(:)         ! corresponding gridcell index for each column
    real(r8), pointer :: wtgcell(:)           ! weight (relative to gridcell) for each column (0-1)
    real(r8), pointer :: qflx_qrgwl(:)        ! qflx_surf at glaciers, wetlands, lakes
    real(r8), pointer :: qflx_drain(:)        ! sub-surface runoff (mm H2O /s)
    real(r8), pointer :: qflx_evap_tot(:)     ! qflx_evap_soi + qflx_evap_veg + qflx_tran_veg
    real(r8), pointer :: qflx_surf(:)         ! surface runoff (mm H2O /s)
!
!EOP
!
! !OTHER LOCAL VARIABLES:
!
    integer :: i,j,k,n,g,l,c,p                ! indices
    integer :: io,jo,ir,jr,is,js              ! mapping indices
    integer :: begp, endp                     ! per-proc beginning and ending pft indices
    integer :: begc, endc                     ! per-proc beginning and ending column indices
    integer :: begl, endl                     ! per-proc beginning and ending landunit indices
    integer :: begg, endg                     ! per-proc gridcell ending gridcell indices
    integer :: numg                           ! total number of gridcells across all processors
    integer :: numl                           ! total number of landunits across all processors
    integer :: numc                           ! total number of columns across all processors
    integer :: nump                           ! total number of pfts across all processors
    integer :: ier                            ! error status
    integer :: nstep                          ! time step index
!-----------------------------------------------------------------------

#if (defined TIMING_BARRIERS)
    call t_startf ('sync_clmrtm')
    call mpi_barrier (mpicom, ier)
    call t_stopf ('sync_clmrtm')
#endif

   ! Assign local pointers to derived type members

    ixy           => ldecomp%gdc2i
    jxy           => ldecomp%gdc2j
    cgridcell     => clm3%g%l%c%gridcell
    wtgcell       => clm3%g%l%c%wtgcell
    qflx_qrgwl    => clm3%g%l%c%cwf%qflx_qrgwl
    qflx_drain    => clm3%g%l%c%cwf%qflx_drain
    qflx_evap_tot => clm3%g%l%c%cwf%pwf_a%qflx_evap_tot
    qflx_surf     => clm3%g%l%c%cwf%qflx_surf

    ! Determine subgrid bounds for this processor

    call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)
    call get_proc_global(numg, numl, numc, nump)

    ! Make gridded representation of runoff
    ! total surface runoff = surface runoff on soils + runoff on glaciers,  wetlands, lakes (P-E)

    rtmin(:) = 0._r8
    do c = begc, endc
       g = cgridcell(c)
       rtmin(g) = rtmin(g) + (qflx_surf(c) + qflx_qrgwl(c) + qflx_drain(c)) * wtgcell(c)
    end do

    ! RTM input averaging is not done

    if (rtm_nsteps <= 1) then
#if (defined SPMD)
       call allgather_data(rtmin, rtmin_glob, clmlevel=nameg)
#else
       rtmin_glob => rtmin
#endif
       totruninxy(:,:) = 0._r8
       do g = 1,numg
          i = ixy(g)
          j = jxy(g)
          totruninxy(i,j) = rtmin_glob(g)
       end do
       delt_rtm = get_step_size()
       do_rtm = .true.
    end if

    ! RTM input averaging is done - average fluxes for RTM calculation

    if (rtm_nsteps > 1) then
       do g = begg,endg
          rtmin_ave(g) = rtmin_ave(g) + rtmin(g)
       end do
       ncount_rtm = ncount_rtm + 1
       nstep = get_nstep()
       if ((mod(nstep,rtm_nsteps)==0) .and. (nstep>1)) then
#if (defined SPMD)
          call allgather_data(rtmin_ave, rtmin_glob, clmlevel=nameg)
#else
          rtmin_glob => rtmin_ave
#endif
          totruninxy(:,:) = 0._r8
          do g = 1,numg
             i = ixy(g)
             j = jxy(g)
             totruninxy(i,j) = rtmin_glob(g)/ncount_rtm
          end do
          delt_rtm = ncount_rtm*get_step_size()   !compute delt for rtm
          ncount_rtm = 0                          !reset counter to 0
          do g = begg,endg
             rtmin_ave(g) = 0._r8                    !reset averager
          end do
          do_rtm = .true.
       else
          do_rtm = .false.
       endif
    endif

  end subroutine RtmUpdateInput

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Rtm
!
! !INTERFACE:
  subroutine Rtm
!
! !DESCRIPTION:
! River routing model (based on U. Texas code).
! Input is totrunin\_r.
! Input/output is fluxout, volr.
! Outputs are dvolrdt\_r, dvolrdt\_lnd\_r, dvolrdt\_ocn\_r, flxocn\_r, flxlnd\_r.
!
! !USES:
!
! !ARGUMENTS:
    implicit none
!
! !CALLED FROM:
! subroutine Rtmriverflux in this module
!
! !REVISION HISTORY:
! Author: Sam Levis
!
!EOP
!
! !LOCAL VARIABLES:
    integer  :: i, j, n                     !loop indices
    real(r8) :: dvolrdt                     !change in storage (m3/s)
    real(r8) :: sumdvolr(rtmlat)            !global sum (m3/s)
    real(r8) :: sumrunof(rtmlat)            !global sum (m3/s)
    real(r8) :: sumdvolr_tot                !global sum (m3/s)
    real(r8) :: sumrunof_tot                !global sum (m3/s)
!-----------------------------------------------------------------------

    ! The following assumes that the river routing data has no communication across the
    ! poles - this is true for the river routing file rdirc.05.

    fluxout(0,:)        = fluxout(rtmlon,:)
    fluxout(rtmlon+1,:) = fluxout(1,:)

    ! Determine cell-to-cell transport - calculate sfluxin

!$OMP PARALLEL DO PRIVATE (i,j)
#if !defined (_OPENMP)
!CSD$ PARALLEL DO PRIVATE (i,j)
#endif
    do j=rtmlati,rtmlatf
       do i=rtmloni,rtmlonf
          sfluxin(i,j) = 0._r8
          if (rdirc(i  ,j-1)==1) sfluxin(i,j) = sfluxin(i,j) + fluxout(i  ,j-1)
          if (rdirc(i-1,j-1)==2) sfluxin(i,j) = sfluxin(i,j) + fluxout(i-1,j-1)
          if (rdirc(i-1,j  )==3) sfluxin(i,j) = sfluxin(i,j) + fluxout(i-1,j  )
          if (rdirc(i-1,j+1)==4) sfluxin(i,j) = sfluxin(i,j) + fluxout(i-1,j+1)
          if (rdirc(i  ,j+1)==5) sfluxin(i,j) = sfluxin(i,j) + fluxout(i  ,j+1)
          if (rdirc(i+1,j+1)==6) sfluxin(i,j) = sfluxin(i,j) + fluxout(i+1,j+1)
          if (rdirc(i+1,j  )==7) sfluxin(i,j) = sfluxin(i,j) + fluxout(i+1,j  )
          if (rdirc(i+1,j-1)==8) sfluxin(i,j) = sfluxin(i,j) + fluxout(i+1,j-1)
       enddo
    enddo
#if !defined (_OPENMP)
!CSD$ END PARALLEL DO
#endif
!$OMP END PARALLEL DO

    ! Loops above and below must remain separate because fluxout is updated below

    sumdvolr(:) = 0._r8
    sumrunof(:) = 0._r8
!$OMP PARALLEL DO PRIVATE (i,j,dvolrdt)
#if !defined (_OPENMP)
!CSD$ PARALLEL DO PRIVATE (i,j,dvolrdt)
#endif
    do j = rtmlati,rtmlatf
       do i = rtmloni,rtmlonf
          n = (j-1)*rtmlon + i

          ! calculate change in cell storage volume change units for
          ! totrunin from kg/m2s==mm/s -> m3/s

          dvolrdt = sfluxin(i,j) - fluxout(i,j) + 0.001_r8*totrunin_r(i,j)*rivarea(i,j)

          ! calculate flux out of a cell:
          ! land: do not permit change in cell storage volume greater than volume present
          ! make up for the difference with an extraction term (eg from aquifers)
          ! ocean: do not permit negative change in cell storage volume,
          ! because at ocean points cell storage volume equals zero
          ! water balance check (in mm/s), convert runinxy from mm/s to m/s (* 1.e-3)
          ! and land model area from km**2 to m**2 (* 1.e6)

          if (rdomain%mask(n) == 1) then         ! land points
             volr(i,j)     = volr(i,j) + dvolrdt*delt_rtm
             fluxout(i,j)  = volr(i,j) * effvel/ddist(i,j)
             fluxout(i,j)  = min(fluxout(i,j), volr(i,j) / delt_rtm)
             flxlnd_r(i,j) = fluxout(i,j)
             flxocn_r(i,j) = 0._r8
             dvolrdt_lnd_r(i,j) = 1000._r8*dvolrdt/rivarea(i,j)
             dvolrdt_ocn_r(i,j) = 0._r8
          else                               ! ocean points
             flxlnd_r(i,j) = 0._r8
             flxocn_r(i,j) = dvolrdt
             dvolrdt_lnd_r(i,j) = 0._r8
             dvolrdt_ocn_r(i,j) = 1000._r8*dvolrdt/rivarea(i,j)
          endif
          sumdvolr(j) = sumdvolr(j) + dvolrdt
          sumrunof(j) = sumrunof(j) + totrunin_r(i,j)*1000._r8*rdomain%area(n)
          dvolrdt_r(i,j) = 1000._r8*dvolrdt/rivarea(i,j)

       enddo
    enddo
#if !defined (_OPENMP)
!CSD$ END PARALLEL DO
#endif
!$OMP END PARALLEL DO

    ! Global water balance calculation and error check

    sumdvolr_tot = 0._r8
    sumrunof_tot = 0._r8
    do j = 1,rtmlat
       sumdvolr_tot = sumdvolr_tot + sumdvolr(j)
       sumrunof_tot = sumrunof_tot + sumrunof(j)
    end do
    if (abs((sumdvolr_tot-sumrunof_tot)/sumrunof_tot) > 0.01_r8) then
       write(6,*) 'RTM Error: sumdvolr= ',sumdvolr_tot,&
            ' not equal to sumrunof= ',sumrunof_tot
       call endrun
    end if

  end subroutine Rtm

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: RtmRest
!
! !INTERFACE:
  subroutine RtmRest(ncid, flag)
!
! !DESCRIPTION:
! Read/write RTM restart data.
!
! !USES:
    use runoffMod   , only : runoff
    use ncdio
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: ncid            ! netcdf id
    character(len=*), intent(in) :: flag   ! 'read' or 'write'
!
! !CALLED FROM:
! subroutine restart in module restFileMod
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
!
! !OTHER LOCAL VARIABLES:
    logical :: readvar          ! determine if variable is on initial file
!-----------------------------------------------------------------------

    ! water volume in cell (m^3)
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='RTMVOLR', xtype=nf_double,  &
            dim1name='rtmlon', dim2name='rtmlat', &
            long_name='water volumn in cell (volr)', units='m3')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_ioglobal(varname='RTMVOLR', data=volr, &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) then
             call endrun()
          else
             volr(:,:) = 0._r8
          end if
       end if
    end if

    ! water flux out of cell (m^3/s)
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='RTM_FLUXOUT', xtype=nf_double,  &
            dim1name='rtmlon', dim2name='rtmlat', &
            long_name='water volumn in cell (volr)', units='m3')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_ioglobal(varname='RTM_FLUXOUT', data=fluxout, &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) then
             call endrun()
          else
             call RTMfluxini()
          end if
       end if
    end if

    ! counter for rtm averaged input (on clm grid)
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='RTM_NCOUNT', xtype=nf_int,  &
            long_name='counter for RTM averaged input on CLM grid', units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_ioglobal(varname='RTM_NCOUNT', data=ncount_rtm, &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) then
             call endrun()
          end if
       end if
    end if

    ! rtm averaged input (on clm grid)
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='RTM_INPUT', xtype=nf_double,  &
            dim1name='gridcell', long_name='RTM averaged input on CLM grid', units='mm/s')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='RTM_INPUT', data=rtmin_ave, &
            dim1name='gridcell', ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    ! ocean river runoff - for history output 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='RTM_OCNROF', xtype=nf_double,  &
            dim1name='ocnrof', long_name='RTM ocean river runoff', units='m3/s')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_ioglobal(varname='RTM_OCNROF', data=runoff%ocn, &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    ! land river runoff - for history output 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='RTM_LNDROF', xtype=nf_double,  &
            dim1name='lndrof', long_name='RTM land river runoff', units='m3/s')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_ioglobal(varname='RTM_LNDROF', data=runoff%lnd, &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    ! ocean dvolrdt - for history output 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='RTM_OCN_DVOLRDT', xtype=nf_double, &
            dim1name='ocnrof', long_name='RTM river flow into ocean change', units='mm/s')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_ioglobal(varname='RTM_OCN_DVOLRDT', data=runoff%ocn_dvolrdt, &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    ! land dvolrdt - for history output 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='RTM_LND_DVOLRDT', xtype=nf_double, &
            dim1name='lndrof', long_name='RTM river flow into ocean change', units='mm/s')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_ioglobal(varname='RTM_LND_DVOLRDT', data=runoff%lnd_dvolrdt, &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

  end subroutine RtmRest

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: RtmMapClm2Rtm
!
! !INTERFACE:
  subroutine RtmMapClm2Rtm( field_s, field_r, landfrac_scale)
!
! !USES:
    use domainMod , only : ldomain
    use areaMod   , only : gridmap_type,areaini,areaave
!
! !ARGUMENTS:
    implicit none
    real(r8), intent(in) :: field_s(:,:)             ! input on clm grid
    real(r8), intent(out):: field_r(rtmlon, rtmlat)  ! output on rtm grid
    logical , intent(in), optional :: landfrac_scale ! scale field_s by landfrac
!
!EOP
!
! !PRIVATE TYPES:
!
    real(r8), dimension(4) :: rtmedge = (/ 90._r8, 180._r8, -90._r8, -180._r8 /)  ! RTM grid N,E,S,W edges 
    integer  :: i,j,ir,jr,is,js,n             ! indices
    real(r8),allocatable :: maskone_s(:)      ! global dummy field
    real(r8),allocatable :: maskone_r(:)      ! global dummy field
    real(r8) :: offset                        ! offset for interpolation from model->rtm grid
    integer  :: ier                           ! error code
    type (gridmap_type),save :: gridmap_l2r
    real(r8):: wt                             ! weight
    logical :: init = .true.
!-----------------------------------------------------------------------

    if (init) then 

       init = .false.

       allocate(maskone_s(lsmlon*lsmlat))
       allocate(maskone_r(rtmlon*rtmlat))

       ! --------------------------------------------------------------------
       ! Map weights from land model grid to rtm grid (global calculation)
       ! --------------------------------------------------------------------
       
       if (masterproc) then
          write(6,*) 'Initializing land model -> rtm interpolation .....'
       endif
    
       ! To find fraction of each land model grid cell that is land based on 
       ! rtm grid.
       ! For this purpose, want all rtm grid cells to contribute to grid cell
       ! average on land model grid, i.e., all cells used regardless of 
       ! whether land or ocean. Do this by setting [maskone_s] = 1
       ! [maskone_s]=1 means all grid cells on land model grid, regardless of 
       ! whether land or ocean, will contribute to rtm grid.

       maskone_s = 1._r8

       ! [maskone_r] = 1 means all the rtm grid is land. Used as dummy
       ! variable so code will not abort with false, non-valid error check

       maskone_r = 1._r8

       ! For each rtm grid cell: get lat [jovr_s2r] and lon [iovr_s2r] indices
       ! and weights [wovr_s2r] of overlapping atm grid cells

       call areaini (  ldomain, rdomain, gridmap_l2r, &
            fracin=maskone_s, fracout=maskone_r)

       deallocate(maskone_s,maskone_r)

       if (masterproc) then
          write(6,*) 'Successfully made land model -> rtm interpolation'
          write(6,*)
       endif

    end if

    ! Map RTM inputs from land model grid to RTM grid (1/2 degree resolution)

    call areaave(field_s, field_r, gridmap_l2r, scale_i=ldomain%frac)

  end subroutine RtmMapClm2Rtm

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: is_restart
!
! !INTERFACE:
  logical function is_restart( )
!
! !DESCRIPTION:
! Determine if restart run
!
! !USES:
    use clm_varctl, only : nsrest
!
! !ARGUMENTS:
    implicit none
!
! !CALLED FROM:
! subroutine initialize in this module
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!-----------------------------------------------------------------------

    if (nsrest == 1) then
       is_restart = .true.
    else
       is_restart = .false.
    end if

  end function is_restart

#endif

end module RtmMod
