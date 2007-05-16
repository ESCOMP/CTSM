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
  use spmdMod     , only : masterproc,npes,iam,mpicom,comp_id,MPI_REAL8
  use clm_varpar  , only : lsmlon, lsmlat, rtmlon, rtmlat
  use shr_sys_mod , only : shr_sys_flush
  use domainMod   , only : latlon_type, latlon_init, domain_clean
  use abortutils  , only : endrun
  use RunoffMod   , only : runoff
  use RunoffMod   , only : gsMap_rtm_gdc2glo,perm_rtm_gdc2glo,sMatP_l2r
  use clm_mct_mod
  use perf_mod
!
! !PUBLIC TYPES:
  implicit none
  private
!
! !PUBLIC MEMBER FUNCTIONS:
  public Rtmini        ! Initialize RTM grid and land mask
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
  private RtmUpdateInput  ! Update rtm inputs
  private Rtm             ! River routing model (based on U. Texas code)

! !PRIVATE TYPES:

! RTM input
  real(r8), pointer :: rtmin_acc(:)        ! RTM averaging buffer for runoff
  real(r8), pointer :: rtmin_avg(:)       ! RTM global input
  integer  :: ncount_rtm                   ! RTM time averaging = number of time samples to average over
  real(r8) :: delt_rtm                     ! RTM time step
  real(r8) :: delt_rtm_max                 ! RTM max timestep
  real(r8) :: cfl_scale = 0.1_r8           ! cfl scale factor, must be <= 1.0
  real(r8), parameter :: effvel = 0.35_r8  ! downstream velocity (m/s)

!glo
  integer , pointer :: dwnstrm_index(:)! downstream index

!gdc
  real(r8), pointer :: volr(:)         ! cell h2o volume (m^3)
  real(r8), pointer :: ddist(:)        ! downstream dist (m)
  real(r8), pointer :: evel(:)         ! effective velocity (m/s)
  real(r8), pointer :: sfluxin(:)      ! cell h2o influx (m3/s)
  real(r8), pointer :: fluxout(:)      ! cell h2o outlflux (m3/s)
  real(r8), pointer :: totrunin(:)     ! cell h2o lnd forcing on rtm grid (mm/s)

!map
  type(mct_sMat)     :: sMat0_l2r
  type(mct_sMat)     :: sMat0_l2r_d

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
! Initialize RTM grid, mask, decomp
!
! !USES:
    use shr_const_mod, only : SHR_CONST_PI
    use domainMod    , only : llatlon,ldomain,domain_check
    use areaMod      , only : celledge, cellarea, map_setmapsAR
    use clm_varctl   , only : frivinp_rtm
    use clm_varctl   , only : rtm_nsteps
    use clm_varcon   , only : re
    use decompMod    , only : get_proc_bounds, get_proc_global, ldecomp
    use decompMod    , only : gsMap_lnd_gdc2glo, perm_lnd_gdc2glo
    use clm_time_manager, only : get_curr_date
    use clm_time_manager, only : get_step_size
    use spmdMod
    use spmdGathScatMod
!
! !ARGUMENTS:
    implicit none
!
! !CALLED FROM:
! subroutine initialize in module initializeMod
!
! !REVISION HISTORY:
! Author: Sam Levis
! Update: T Craig, Dec 2006
!
!EOP
!
! !LOCAL VARIABLES:
    real(r8), dimension(4) :: rtmedge = (/ 90._r8, 180._r8, -90._r8, -180._r8 /)  !N,E,S,W edges of rtm grid
    integer  :: ioff(0:8) = (/0,0,1,1,1,0,-1,-1,-1/) !rdirc input to i
    integer  :: joff(0:8) = (/0,1,1,0,-1,-1,-1,0,1/) !rdirc input to j
    integer  :: i,j,k,n,g,n2                  ! loop indices
    integer  :: im1,ip1,jm1,jp1,ir,jr,nr      ! neighbor indices
    integer  :: i2,j2                         ! downstream i and j
    real(r8) :: deg2rad                       ! pi/180
    real(r8) :: dx,dx1,dx2,dx3                ! lon dist. betn grid cells (m)
    real(r8) :: dy                            ! lat dist. betn grid cells (m)
    real(r8),allocatable :: tempg(:)          ! temporary buffer
    integer ,allocatable :: rdirc(:)          ! temporary buffer
    integer ,allocatable :: iocn(:)           ! downstream ocean cell
    integer ,allocatable :: nocn(:)           ! number of rtm cells in basin
    integer ,allocatable :: pocn(:)           ! pe number assigned to basin
    integer ,allocatable :: nop(:)            ! number of rtm cells on a pe
    integer  :: pemin                         ! pe with min number of rtm cells
    integer  :: ier                           ! error code
    integer  :: mon                           ! month (1, ..., 12)
    integer  :: day                           ! day of month (1, ..., 31)
    integer  :: ncsec                         ! seconds of current date
    integer  :: begg,endg                     ! local start/end gridcell indices
    integer  :: numg                          ! tot num of gridcells on all pes
    integer  :: numl                          ! tot num of landunits on all pes
    integer  :: numc                          ! tot num of columns on all pes
    integer  :: nump                          ! tot num of pfts on all pes
    integer  :: begr,endr,numr                ! tot num of roff pts on all pes
    real(r8) :: dtover,dtovermax              ! ts calc temporaries
    character(len=16), dimension(50) :: river_name
    character(len=30), dimension(50) :: rivstat_name
    real(r8)         , dimension(50) :: rivstat_lon
    real(r8)         , dimension(50) :: rivstat_lat
    integer  :: nroflnd
    integer  :: nrofocn
    integer  :: pid,np,npmin,npmax,npint      ! log loop control
    integer,parameter  :: dbug = 1            ! 0 = none, 1=normal, 2=much, 3=max

    integer,pointer :: gindex(:)           ! index for permute
    integer lsize,gsize                    ! temporary for permute
    integer  :: na,nb,ns                   ! mct sizes
    integer  :: igrow,igcol,iwgt           ! mct field indices
    integer  :: ii,ji,ni,no,gi,go          ! tmps
    real(r8) :: wt                         ! mct wt
    real(r8),pointer :: lfield(:)          ! tmp lnd field
    real(r8),pointer :: rfield(:)          ! tmp rtm field
    real(r8),pointer :: glatc(:),glonc(:)  ! global lat/lon
    real(r8),pointer :: gfrac(:)           ! global frac
    real(r8),pointer :: lfrac(:)           ! global land frac
    integer ,pointer :: gmask(:)           ! global mask
    type(latlon_type):: rlatlon            ! rtm grid 

!-----------------------------------------------------------------------

    !--- Allocate rtm grid variables
    call latlon_init(rlatlon,rtmlon,rtmlat)

    !--- Allocate inputs and outputs to rtm at 1/2 degree resolution

    allocate (dwnstrm_index  (rtmlon*rtmlat), &
              runoff%rlat(rtmlat), runoff%rlon(rtmlon), &
              glatc(rtmlon*rtmlat), glonc(rtmlon*rtmlat), &
              gfrac(rtmlon*rtmlat), lfrac(lsmlon*lsmlat), &
              gmask(rtmlon*rtmlat), &
              stat=ier)
    if (ier /= 0) then
       write(6,*)'Rtmgridini: Allocation error for ',&
            'dwnstrm_index,rlat,rlon'
       call endrun
    end if

    !--- Allocate temporaries

    allocate(tempg(rtmlon*rtmlat),rdirc(rtmlon*rtmlat), &
         iocn(rtmlon*rtmlat),nocn(rtmlon*rtmlat),stat=ier)
    if (ier /= 0) then
       write(6,*)'Rtmgridini: Allocation error for ',&
            'tempg'
       call endrun
    end if

    !--- Useful constants and initial values

    deg2rad = SHR_CONST_PI / 180._r8

    ! Open and read input data (river direction file)
    ! rtm operates from south to north and from the dateline
    ! River station data is currently not used by the model -
    ! it is only used by the diagnostic package
    ! If the river direction file is modified - the river station
    ! part must also be modified

    rlatlon%edges(1:4) = rtmedge(1:4)
    open (1,file=frivinp_rtm)
    do n = 1,rtmlon*rtmlat
       read(1,*) glatc(n),glonc(n),tempg(n)
       rdirc(n) = nint(tempg(n))
    enddo
    do n = 1,50
       read(1,10,iostat=ier) river_name(n), rivstat_lon(n), rivstat_lat(n), &
                             rivstat_name(n)
       if (ier /= 0) exit
10     format(1x,a16,f7.2,1x,f7.2,a30)
    end do
    close(1)

    if (masterproc) then
       write(6,*)'Columns in RTM = ',rtmlon
       write(6,*)'Rows in RTM    = ',rtmlat
       write(6,*)'read river direction data'
    end if

    !--- set 1d lat/lon values

    do j=1,rtmlat
       n = (j-1)*rtmlon + 1
       runoff%rlat(j) = glatc(n)
       rlatlon%latc(j) = glatc(n)
    enddo
    do i=1,rtmlon
       n = i
       runoff%rlon(i) = glonc(n)
       rlatlon%lonc(i) = glonc(n)
    enddo

    !--- Set dwnstrm_index from rdirc values
    !--- The following assumes that there is no runoff 
    !---   south of j=1 or north of j=rtmlat
    !--- This is true for rdirc.05
    !--- Determine dwnstrmm_index from rtm river flow direction (0-8)

    dwnstrm_index = 0
    do j=1,rtmlat
    do i=1,rtmlon
       n = (j-1)*rtmlon + i
       if (rdirc(n) /= 0) then
          ir = i + ioff(rdirc(n))
          jr = j + joff(rdirc(n))
          if (ir < 1     ) ir = ir + rtmlon
          if (ir > rtmlon) ir = ir - rtmlon
          !--- check cross pole flow, etc
          if (jr < 1 .or. jr > rtmlat .or. ir < 1 .or. ir > rtmlon) then
             write(6,*) 'Rtmini ERROR ir jr bounds ',i,j,rdirc(n),ir,jr
             call endrun()
          endif
          nr = (jr-1)*rtmlon + ir
          if (n == nr) then
             write(6,*) 'Rtmini ERROR dwnstrm_index ',i,j,n,rdirc(n),ir,jr,nr
             call endrun()
          endif
          dwnstrm_index(n) = nr
       endif
    enddo
    enddo

    !--- Determine RTM celledges and areas 

    call celledge (rlatlon, &
                   rlatlon%edges(1), rlatlon%edges(2), &
                   rlatlon%edges(3), rlatlon%edges(4))


    !--- Set sMat0_l2r, full mapping weights for l2r, just on root pe
    !--- for now use lfield to "ignore" non-active land cells in sMat0_l2r
    !--- Later these will be "reduced" to just the useful weights
    !--- Compute rdomain%frac on root pe and bcast 

    call get_proc_bounds(begg, endg)
    call gather_data_to_master(ldomain%frac,lfrac,gsmap_lnd_gdc2glo,perm_lnd_gdc2glo,begg,endg)

    if (masterproc) then
       allocate(lfield(lsmlon*lsmlat),rfield(rtmlon*rtmlat))
       lfield = 0._r8
       do n = 1,lsmlon*lsmlat
          if (ldecomp%glo2gdc(n) > 0) lfield(n) = 1._r8
       enddo
       rfield = 1._r8

       call map_setmapsAR(llatlon,rlatlon,sMat0_l2r, &
          fracin=lfield, fracout=rfield)
       igrow = mct_sMat_indexIA(sMat0_l2r,'grow')
       igcol = mct_sMat_indexIA(sMat0_l2r,'gcol')
       iwgt  = mct_sMat_indexRA(sMat0_l2r,'weight')
       gfrac = 0._r8
       do n = 1,mct_sMat_lsize(sMat0_l2r)
          nr = sMat0_l2r%data%iAttr(igrow,n)
          ns = sMat0_l2r%data%iAttr(igcol,n)
          wt = sMat0_l2r%data%rAttr(iwgt ,n)
          gfrac(nr) = gfrac(nr) + lfrac(ns)
       enddo
       deallocate(lfield,rfield)
    endif
    call mpi_bcast(gfrac,size(gfrac),MPI_REAL8,0,mpicom,ier)

    !--- Determine rtm ocn/land mask, 0=none, 1=land, 2=ocean

    gmask = 0             ! assume neither land nor ocn

    do n=1,rtmlon*rtmlat         ! set downstream value first
       nr = dwnstrm_index(n)
       if (nr /= 0) then         ! assume downstream cell is ocn
          gmask(nr) = 2
       end if
    enddo

    do n=1,rtmlon*rtmlat         ! override downstream setting from local info
       nr = dwnstrm_index(n)
       if (nr /= 0) then         ! n is always land if dwnstrm_index exists
          gmask(n) = 1
       else                      ! n is ocn if no dwnstrm_index and some frac
          if (gfrac(n)>0._r8) gmask(n) = 2
       end if
    enddo

   !--- Compute river basins, actually compute ocean outlet gridcell
   !--- iocn = final downstream cell, index is global 1d ocean gridcell
   !--- nocn = number of source gridcells for ocean gridcell

    iocn = 0
    nocn = 0
    do nr=1,rtmlon*rtmlat
       n = nr
       if (gmask(n) == 1) then    ! land
          g = 0
          do while (gmask(n) == 1 .and. g < rtmlon*rtmlat)  ! follow downstream
             n = dwnstrm_index(n)
             g = g + 1
          end do
          if (gmask(n) == 2) then  ! found ocean outlet
             iocn(nr) = n                 ! set ocean outlet or nr to n
             nocn(n) = nocn(n) + 1        ! one more land cell for n
          elseif (gmask(n) == 1) then  ! no ocean outlet, warn user, ignore cell
             write(6,*) 'rtmini WARNING no downstream ocean cell - IGNORED', &
               g,nr,gmask(nr),dwnstrm_index(nr), &
               n,gmask(n),dwnstrm_index(n)
          else 
             write(6,*) 'rtmini ERROR downstream cell is non-ocean,non-land', &
               g,nr,gmask(nr),dwnstrm_index(nr), &
               n,gmask(n),dwnstrm_index(n)
             call endrun()
          endif
       elseif (gmask(n) == 2) then  ! ocean, give to self
          iocn(nr) = n
          nocn(n) = nocn(n) + 1
       endif
    enddo

    !--- Now allocate those basins to pes
    !--- pocn is the pe that gets the basin associated with ocean outlet nr
    !--- nop is a running count of the number of rtm cells/pe 

    allocate(pocn(rtmlon*rtmlat),nop(0:npes-1),runoff%glo2gdc(rtmlon*rtmlat),runoff%num_rtm(0:npes-1))
    nop = 0
    pocn = -99
    runoff%glo2gdc = 0
    do nr=1,rtmlon*rtmlat
       if (nocn(nr) /= 0) then
          pemin = 0
          do n = 1,npes-1
             if (nop(n) < nop(pemin)) pemin = n
          enddo
          if (pemin > npes-1 .or. pemin < 0) then
             write(6,*) 'error in decomp for rtm ',nr,npes,pemin
             call endrun()
          endif
          nop(pemin) = nop(pemin) + nocn(nr)
          pocn(nr) = pemin
       endif
    enddo

    !--- Count and distribute cells to runoff%glo2gdc

    runoff%numr   = 0
    runoff%numro  = 0
    runoff%numrl  = 0
    runoff%lnumr  = 0
    runoff%lnumro = 0
    runoff%lnumrl = 0
    runoff%num_rtm = 0
    g = 0
    do n = 0,npes-1
       if (iam == n) then
          runoff%begr  = runoff%numr  + 1
          runoff%begrl = runoff%numrl + 1
          runoff%begro = runoff%numro + 1
       endif
       do nr=1,rtmlon*rtmlat
          if (pocn(nr) == n .and. nocn(nr) /= 0) then
             runoff%num_rtm(n) = runoff%num_rtm(n) + nocn(nr)
             runoff%numr  = runoff%numr  + nocn(nr)
             runoff%numro = runoff%numro + 1
             runoff%numrl = runoff%numrl + nocn(nr) - 1
             k = g
             if (nocn(nr) == 1) then   ! avoid the double rtm nested loop
                n2 = nr
                g = g + 1
                runoff%glo2gdc(n2) = g
             else
                do n2 = 1,rtmlon*rtmlat
                   if (iocn(n2) == nr) then
                      g = g + 1
                      runoff%glo2gdc(n2) = g
                   endif
                enddo
             endif
             if ((g-k) /= nocn(nr)) then
                write(6,*) 'Rtmini ERROR rtm cell count ',n,nr,k,g,g-k,nocn(nr)
                call endrun()
             endif
             if (iam == n) then
                runoff%lnumr  = runoff%lnumr  + nocn(nr)
                runoff%lnumro = runoff%lnumro + 1
                runoff%lnumrl = runoff%lnumrl + nocn(nr) - 1
             endif
          endif
       enddo
       if (iam == n) then
          runoff%endr  = runoff%numr
          runoff%endro = runoff%begro + runoff%lnumro - 1
          runoff%endrl = runoff%begrl + runoff%lnumrl - 1
       endif
    enddo

    !--- set some local values

    nroflnd = runoff%numrl
    nrofocn = runoff%numro
    numr = nroflnd + nrofocn
    begr = runoff%begr
    endr = runoff%endr

    !--- Write per-processor runoff bounds depending on dbug level

#ifndef UNICOSMP
    call shr_sys_flush(6)
#endif
    if (masterproc) then
       write(6,*) 'total runoff cells numr = ',runoff%numr, &
          'numrl = ',runoff%numrl,'numro = ',runoff%numro
    endif
#ifndef UNICOSMP
    call shr_sys_flush(6)
#endif
    call mpi_barrier(mpicom,ier)
    npmin = 0
    npmax = npes-1
    npint = 1
    if (dbug == 0) then
       npmax = 0
    elseif (dbug == 1) then
       npmax = min(npes-1,4)
    elseif (dbug == 2) then
       npint = npes/8
    endif
    do np = npmin,npmax,npint
       pid = np
       if (dbug == 1) then
          if (np == 2) pid=npes/2-1
          if (np == 3) pid=npes-2
          if (np == 4) pid=npes-1
       endif
       pid = max(pid,0)
       pid = min(pid,npes-1)
       if (iam == pid) then
          write(6,*) 'rtm decomp info',' proc = ',iam, &
             ' begr = ',runoff%begr,' endr = ',runoff%endr, &
             ' numr = ',runoff%lnumr
          write(6,*) '               ',' proc = ',iam, &
             ' begrl= ',runoff%begrl,' endrl= ',runoff%endrl, &
             ' numrl= ',runoff%lnumrl
          write(6,*) '               ',' proc = ',iam, &
             ' begro= ',runoff%begro,' endro= ',runoff%endro, &
             ' numro= ',runoff%lnumro
       endif
#ifndef UNICOSMP
       call shr_sys_flush(6)
#endif
       call mpi_barrier(mpicom,ier)
    enddo

    !--- allocate runoff variables

    allocate(runoff%runoff(begr:endr),runoff%dvolrdt(begr:endr), &
             runoff%area(begr:endr), &
             runoff%lonc(begr:endr),  runoff%latc(begr:endr),  &
             runoff%dsi(begr:endr), stat=ier)
    if (ier /= 0) then
       write(6,*)'Rtmini ERROR allocation of runoff%runoff'
       call endrun
    end if

    allocate(runoff%gdc2glo(numr), runoff%mask(numr), &
             runoff%gdc2gsn(numr), &
             runoff%gdc2i(numr),runoff%gdc2j(numr),stat=ier)
    if (ier /= 0) then
       write(6,*)'Rtmini ERROR allocation of runoff%gcd2glo'
       call endrun
    end if

    !--- Allocate rtm flux variables

    allocate (volr    (begr:endr), &
              fluxout (begr:endr), &
              ddist   (begr:endr), &
              totrunin(begr:endr), &
              evel    (begr:endr), &
              sfluxin (begr:endr),  stat=ier)
    if (ier /= 0) then
       write(6,*)'Rtmgridini: Allocation error for ',&
            'volr, fluxout, ddist'
       call endrun
    end if
    volr = 0._r8
    fluxout = 0._r8
    ddist = 0._r8
    evel = effvel
    sfluxin = 0._r8

    !--- Initialize runoff data

    numr = 0
    do j = 1,rtmlat
    do i = 1,rtmlon
       n = (j-1)*rtmlon + i
       nr = runoff%glo2gdc(n)
       if (nr /= 0) then
          numr = numr + 1
          runoff%gdc2glo(nr) = n         
          runoff%gdc2i(nr) = i         
          runoff%gdc2j(nr) = j         
          runoff%gdc2gsn(nr) = numr
          runoff%mask(nr) = gmask(n)   ! global for hist file
       endif
    enddo
    enddo

    if (numr /= runoff%numr) then
       write(6,*) 'Rtmini ERROR numr numr ',numr,runoff%numr
       call endrun()
    endif

    do nr = begr,endr
       n = runoff%gdc2glo(nr)
       i = runoff%gdc2i(nr)
       j = runoff%gdc2j(nr)
       if (n <= 0 .or. n > rtmlon*rtmlat) then
          write(6,*) 'Rtmini ERROR gdc2glo ',nr,runoff%gdc2glo(nr)
          call endrun()
       endif
       runoff%runoff(nr) = 0._r8
       runoff%dvolrdt(nr) = 0._r8
       runoff%lonc(nr) = glonc(n)
       runoff%latc(nr) = glatc(n)

       dx = (rlatlon%lone(i) - rlatlon%lonw(i)) * deg2rad
       dy = sin(rlatlon%latn(j)*deg2rad) - sin(rlatlon%lats(j)*deg2rad)
       runoff%area(nr) = 1.e6_r8 * dx*dy*re*re
       if (dwnstrm_index(n) == 0) then
          runoff%dsi(nr) = 0
       else
          if (runoff%glo2gdc(dwnstrm_index(n)) == 0) then
             write(6,*) 'Rtmini ERROR glo2gdc dwnstrm ',nr,n,dwnstrm_index(n),runoff%glo2gdc(dwnstrm_index(n))
             call endrun()
          endif
          runoff%dsi(nr) = runoff%glo2gdc(dwnstrm_index(n))
       endif
    enddo

    !--- Determine downstream distance - instead of reading a distance file
    !--- calculate the downstream distance

    do nr=begr,endr
       g = runoff%dsi(nr)
       if (g == 0) then
          ddist(nr) = 0._r8
       elseif (g < begr .or. g > endr) then
          write(6,*) 'Rtmini: error in ddist calc ',nr,g,begr,endr
          call endrun
       else
          dy = deg2rad * abs(runoff%latc(nr)-runoff%latc(g)) * re*1000._r8
          dx = runoff%lonc(nr)-runoff%lonc(g)
          dx1 = abs(dx)
          dx2 = abs(dx+360._r8)
          dx3 = abs(dx-360._r8)
          dx = min(dx1,dx2,dx3)
          dx = deg2rad * dx * re*1000._r8 * &
               0.5_r8*(cos(runoff%latc(nr)*deg2rad)+ &
                       cos(runoff%latc(g)*deg2rad))
          ddist(nr) = sqrt(dx*dx + dy*dy)
       endif
    enddo

    !--- Compute timestep and subcycling number

    if (rtm_nsteps < 1) then
       write(6,*) 'rtm ERROR in rtm_nsteps',rtm_nsteps
       call endrun()
    endif
    delt_rtm = rtm_nsteps*get_step_size()

    dtover = 0._r8
    dtovermax = 0._r8
    do nr=begr,endr
       if (ddist(nr) /= 0._r8) then
          dtover = evel(nr)/ddist(nr)
       else
          dtover = 0._r8
       endif
       dtovermax = max(dtovermax,dtover)
    enddo
    dtover = dtovermax
    call mpi_allreduce(dtover,dtovermax,1,MPI_REAL8,MPI_MAX,mpicom,ier)
    if (dtovermax > 0._r8) then
       delt_rtm_max = (1.0_r8/dtovermax)*cfl_scale
    else
       write(6,*) 'rtmini error in delt_rtm_max ',delt_rtm_max,dtover
       call endrun
    endif
    if (masterproc) write(6,*) 'rtm max timestep = ',delt_rtm_max,' (sec) for cfl_scale = ',cfl_scale
    if (masterproc) write(6,*) 'rtm act timestep ~ ',delt_rtm

    !--- Allocate and initialize rtm input fields on clm decomp

    call get_proc_global(numg, numl, numc, nump)
    call get_proc_bounds(begg, endg)
    allocate (rtmin_avg(begg:endg), rtmin_acc(begg:endg), stat=ier)
    if (ier /= 0) then
       write(6,*)'Rtmlandini: Allocation error for rtmin, rtmin_avg, rtmin_acc'
       call endrun
    end if
    rtmin_avg(:) = 0._r8
    rtmin_acc(:) = 0._r8

    !--- clean up temporaries

    deallocate(tempg,rdirc,iocn,nocn)
    deallocate(pocn,nop,dwnstrm_index,glatc,glonc,gfrac,lfrac,gmask)

   !--- initialization rtm gsmap

    allocate(gindex(begr:endr))
    do n = begr,endr
       gindex(n) = runoff%gdc2glo(n)
    enddo
    lsize = endr-begr+1
    gsize = rtmlon * rtmlat
    allocate(perm_rtm_gdc2glo(lsize),stat=ier)
    call mct_indexset(perm_rtm_gdc2glo)
    call mct_indexsort(lsize,perm_rtm_gdc2glo,gindex)
    call mct_permute(gindex,perm_rtm_gdc2glo,lsize)
    call mct_gsMap_init( gsMap_rtm_gdc2glo, gindex, mpicom, comp_id, lsize, gsize )
    deallocate(gindex)

    !--- initialize sMat0_l2r_d, from sMat0_l2r - remove unused weights
    !--- root pe only

   if (masterproc) then
       na = llatlon%ni * llatlon%nj
       nb = rtmlon * rtmlat
       igrow = mct_sMat_indexIA(sMat0_l2r,'grow')
       igcol = mct_sMat_indexIA(sMat0_l2r,'gcol')
       iwgt  = mct_sMat_indexRA(sMat0_l2r,'weight')

       ns = 0
       do n = 1,mct_sMat_lsize(sMat0_l2r)
          ni = sMat0_l2r%data%iAttr(igcol,n)
          no = sMat0_l2r%data%iAttr(igrow,n)
          if (ldecomp%glo2gdc(ni) > 0 .and. runoff%glo2gdc(no) > 0) then
             ns = ns + 1
          endif
       enddo

       call mct_sMat_init(sMat0_l2r_d, nb, na, ns)

       ns = 0
       do n = 1,mct_sMat_lsize(sMat0_l2r)
          ni = sMat0_l2r%data%iAttr(igcol,n)
          no = sMat0_l2r%data%iAttr(igrow,n)
          if (ldecomp%glo2gdc(ni) > 0 .and. runoff%glo2gdc(no) > 0) then
             ns = ns + 1
             sMat0_l2r_d%data%iAttr(igcol,ns) = sMat0_l2r%data%iAttr(igcol,n)
             sMat0_l2r_d%data%iAttr(igrow,ns) = sMat0_l2r%data%iAttr(igrow,n)
             sMat0_l2r_d%data%rAttr(iwgt ,ns) = sMat0_l2r%data%rAttr(iwgt ,n)
          endif
       enddo
    endif   ! masterproc

    !--- initialize sMatP_l2r, scatter sMat0_l2r_d based on gsmaps
    
    call mct_sMatP_init(sMatP_l2r,  sMat0_l2r_d,  &
                        gsmap_lnd_gdc2glo, gsMap_rtm_gdc2glo, &
                       'Xonly',0,mpicom,comp_id)

#ifdef CPP_VECTOR
   !--- initialize the vector parts of the sMat
   call mct_sMatP_Vecinit(sMatP_l2r)
#endif

   !--- clean up the root sMat0 datatypes

   if (masterproc) then
      call mct_sMat_clean(sMat0_l2r)
      call mct_sMat_clean(sMat0_l2r_d)
      write(6,*) 'Rtmini complete'
   endif

  end subroutine Rtmini

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
    use decompMod      , only : get_proc_bounds, get_proc_global
    use decompMod      , only : gsMap_lnd_gdc2glo, perm_lnd_gdc2glo
    use domainMod      , only : ldomain
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
    integer  :: i,j,n,n2,nr,ns             ! indices
    logical  :: do_rtm                     ! true => perform rtm calculation
    integer  :: begg,endg
    logical  :: usevector = .false.
    type(mct_aVect)    :: aV_lndr,aV_rtmr
!-----------------------------------------------------------------------

    ! Determine RTM inputs on land model grid

    call RtmUpdateInput(do_rtm)

    if (do_rtm) then

       ! Map RTM inputs from land model grid to RTM grid (1/2 degree resolution)

       ns = mct_gsMap_lsize(gsmap_lnd_gdc2glo, mpicom)
       call mct_aVect_init(aV_lndr,rlist='rtminput',lsize=ns)
       ns = mct_gsMap_lsize(gsMap_rtm_gdc2glo, mpicom)
       call mct_aVect_init(aV_rtmr,rlist='rtminput',lsize=ns)

       call get_proc_bounds(begg, endg)
       nr = mct_aVect_indexRA(av_lndr,'rtminput',perrWith='Rtmriverflux')
       do n = begg,endg
          n2 = n-begg+1
          av_lndr%rAttr(nr,n2) = rtmin_avg(n)*ldomain%frac(n)
       enddo

       call mct_aVect_permute  (av_lndr,perm_lnd_gdc2glo)
       call mct_Smat_AvMult    (av_lndr,sMatP_l2r,av_rtmr,vector=usevector)
       call mct_aVect_unpermute(av_rtmr,perm_rtm_gdc2glo)
 
       nr = mct_aVect_indexRA(av_rtmr,'rtminput',perrWith='Rtmriverflux')
       do n = runoff%begr,runoff%endr
          n2 = n-runoff%begr+1
          totrunin(n) = av_rtmr%rAttr(nr,n2)
       enddo

       call mct_aVect_clean(aV_lndr)
       call mct_aVect_clean(aV_rtmr)

       ! Determine RTM runoff fluxes

       call t_startf('rtm_calc')
       call Rtm()
       call t_stopf('rtm_calc')

    end if

  end subroutine Rtmriverflux

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: RtmUpdateInput
!
! !INTERFACE:
  subroutine RtmUpdateInput(do_rtm)
!
! !DESCRIPTION:
! Update RTM inputs.
!
! !USES:
    use clmtype        , only : clm3,nameg
    use decompMod      , only : get_proc_bounds, get_proc_global, ldecomp
    use clm_varctl     , only : rtm_nsteps
    use clm_time_manager   , only : get_step_size, get_nstep
    use spmdGathScatMod, only : scatter_data_from_master, gather_data_to_master, allgather_data
!
! !ARGUMENTS:
    implicit none
    logical , intent(out) :: do_rtm
!
! !CALLED FROM:
! subroutine rtmriverflux
!
! !REVISION HISTORY:
! Author: Sam Levis
!
! !LOCAL VARIABLES:
!
! local pointers to implicit in arguments
!
    integer , pointer :: nxy(:)               ! gricell xy lon index
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

    call t_barrierf('sync_clmrtm', mpicom)

   ! Assign local pointers to derived type members

    nxy           => ldecomp%gdc2glo
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

    do c = begc, endc
       g = cgridcell(c)
       rtmin_acc(g) = rtmin_acc(g) + (qflx_surf(c) + qflx_qrgwl(c) + qflx_drain(c)) * wtgcell(c)
    end do
    ncount_rtm = ncount_rtm + 1

    nstep = get_nstep()
    if (mod(nstep,rtm_nsteps)==0 .and. nstep>1) then
       if (ncount_rtm*get_step_size() /= delt_rtm) then
          write(6,*) 'RtmUpdateInput timestep out of sync ',delt_rtm,ncount_rtm*get_step_size()
!          call endrun
          delt_rtm = ncount_rtm*get_step_size()
       endif
       do g = begg,endg
          rtmin_avg(g) = rtmin_acc(g)/ncount_rtm
          rtmin_acc(g) = 0._r8
       end do
       ncount_rtm = 0                          !reset counter to 0
       do_rtm = .true.
    else
       do_rtm = .false.
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
! Input is rtmin\_avg.
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
    integer  :: i, j, n, ns              !loop indices
    integer  :: ir,jr,nr                    !neighbor indices
    real(r8) :: dvolrdt                     !change in storage (m3/s)
    real(r8) :: sumdvolr_tot                !global sum (m3/s)
    real(r8) :: sumrunof_tot                !global sum (m3/s)
    integer  :: nsub                        !subcyling for cfl
    integer, save :: nsub_save              !previous nsub
    real(r8) :: delt                        !delt associated with subcycling
    real(r8), save :: delt_save             !previous delt
!-----------------------------------------------------------------------

    nsub = int(delt_rtm/delt_rtm_max) + 1
    delt = delt_rtm/float(nsub)

    if (delt /= delt_save) then
       if (masterproc) write(6,*) 'rtm delt update from/to',delt_save,delt,nsub_save,nsub
    endif

    nsub_save = nsub
    delt_save = delt

    sumdvolr_tot = 0._r8
    sumrunof_tot = 0._r8

    runoff%dvolrdt = 0._r8
    runoff%runoff = 0._r8

    !--- subcycling ---
    do ns = 1,nsub

       sfluxin = 0._r8
       do n = runoff%begr,runoff%endr
          nr = runoff%dsi(n)
          if (nr /= 0) then
             if (nr < runoff%begr .or. nr > runoff%endr) then
                write(6,*) 'Rtm ERROR: non-local communication ',n,nr
                call endrun()
             endif
             sfluxin(nr) = sfluxin(nr) + fluxout(n)
          endif
       enddo

       do n = runoff%begr,runoff%endr
          dvolrdt = sfluxin(n) - fluxout(n) + 0.001_r8*totrunin(n)*runoff%area(n)

          if (runoff%mask(n) == 1) then         ! land points
             volr(n)     = volr(n) + dvolrdt*delt
             fluxout(n)  = volr(n) * evel(n)/ddist(n)
!            --- old cfl constraint.  now use subcycling.  for reference only
!            fluxout(n)  = min(fluxout(n), volr(n)/delt)
!            --- this would stop upstream flow if volr/fluxout < 0
!            --- negative flow largely a function of negative forcing
!            --- can still have negative runoff where there is negative
!                forcing over a mask=2 (ocn) point since forcing is put onto
!                the ocean instantly at these points
!            --- also, want to allow negative flow so it doesn't build up
!            fluxout(n) = max(fluxout(n),0._r8)
          else
             volr(n) = 0._r8
             fluxout(n) = 0._r8
          endif

          if (runoff%mask(n) == 1) then
             runoff%runoff(n) = runoff%runoff(n) + fluxout(n)
          elseif (runoff%mask(n) == 2) then
             runoff%runoff(n) = runoff%runoff(n) + dvolrdt
          endif
          runoff%dvolrdt(n) = runoff%dvolrdt(n) + 1000._r8*dvolrdt/runoff%area(n)

          sumdvolr_tot = sumdvolr_tot + dvolrdt
          sumrunof_tot = sumrunof_tot + 0.001_r8*totrunin(n)*runoff%area(n)

       enddo

    enddo

    ! average fluxes over subcycling
    runoff%runoff = runoff%runoff / float(nsub)
    runoff%dvolrdt = runoff%dvolrdt / float(nsub)

    ! Global water balance calculation and error check

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
    integer :: n,g              ! indices
    real(r8),allocatable :: gtmp(:) ! temporary
!-----------------------------------------------------------------------

    allocate(gtmp(rtmlon*rtmlat))

    ! water volume in cell (m^3)
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='RTMVOLR', xtype=nf_double,  &
            dim1name='rtmlon', dim2name='rtmlat', &
            long_name='water volumn in cell (volr)', units='m3')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='RTMVOLR', data=volr, dim1name='allrof', &
            ncid=ncid, flag=flag, nlonxy=rtmlon,nlatxy=rtmlat,readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) then
             call endrun()
          else
             volr = 0._r8
          end if
       end if
    end if

    ! water flux out of cell (m^3/s)
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='RTMFLUXOUT', xtype=nf_double,  &
            dim1name='rtmlon', dim2name='rtmlat', &
            long_name='water fluxout in cell (fluxout)', units='m3/s')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='RTMFLUXOUT', data=fluxout, dim1name='allrof', &
            ncid=ncid, flag=flag, nlonxy=rtmlon,nlatxy=rtmlat,readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) then
             call endrun()
          else
             fluxout = 0._r8
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
          else
             ncount_rtm = 0
          end if
       end if
    end if

    ! rtm averaged input (on clm grid)
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='RTM_INPUT', xtype=nf_double,  &
            dim1name='gridcell', long_name='RTM averaged input on CLM grid', units='mm/s')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='RTM_INPUT', data=rtmin_acc, &
            dim1name='gridcell', ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) then
             call endrun()
          else
             rtmin_acc = 0._r8
          endif          
       end if
    end if

    ! runoff%runoff
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='RTM_RUNOFF', xtype=nf_double,  &
            dim1name='rtmlon', dim2name='rtmlat', &
            long_name='RTM runoff (runoff)', units='m3/s')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='RTM_RUNOFF', data=runoff%runoff, dim1name='allrof', &
            ncid=ncid, flag=flag, nlonxy=rtmlon,nlatxy=rtmlat,readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
             runoff%runoff = 0_r8
       end if
    end if

    ! runoff%dvolrdt
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='RTM_DVOLRDT', xtype=nf_double,  &
            dim1name='rtmlon', dim2name='rtmlat', &
            long_name='water volumn change in cell (dvolrdt)', units='mm/s')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='RTM_DVOLRDT', data=runoff%dvolrdt, dim1name='allrof', &
            ncid=ncid, flag=flag, nlonxy=rtmlon,nlatxy=rtmlat,readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
             runoff%dvolrdt = 0._r8
       end if
    end if

    deallocate(gtmp)

  end subroutine RtmRest

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
