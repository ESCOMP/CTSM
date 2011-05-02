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
  use spmdMod     , only : masterproc,npes,iam,mpicom,comp_id,MPI_REAL8,MPI_INTEGER
  use spmdMod     , only : MPI_MAX,MPI_SUM
  use clm_varpar  , only : lsmlon, lsmlat, rtmlon, rtmlat
  use clm_varcon  , only : re, spval
  use clm_varctl  , only : iulog
  use shr_sys_mod , only : shr_sys_flush
  use abortutils  , only : endrun
  use RunoffMod   , only : runoff, nt_rtm, rtm_tracers
  use RunoffMod   , only : gsMap_rtm_gdc2glo
  use clm_mct_mod
  use perf_mod
!
! !PUBLIC TYPES:
  implicit none
  private
!
! !PUBLIC MEMBER FUNCTIONS:
  public Rtmini          ! Initialize RTM grid and land mask
  public RtmInput        ! Update rtm inputs
  public RtmMapl2r       ! Maps input from clm to rtm grid
  public Rtm             ! River routing model (based on U. Texas code)
  public RtmRest         ! Read/write RTM restart data (netcdf)
!
! !REVISION HISTORY:
! Author: Sam Levis
!
!
! !PRIVATE MEMBER FUNCTIONS:
!

! !PRIVATE TYPES:

! RTM tracers
  character(len=256) :: rtm_trstr   ! tracer string

! RTM input
  real(r8), pointer :: rtmin_acc(:,:)      ! RTM averaging buffer for runoff
  real(r8), pointer :: rtmin_avg(:,:)      ! RTM global input
  integer  :: ncount_rtm = 0               ! RTM time averaging = number of time samples to average over
  real(r8) :: delt_rtm                     ! RTM time step
  real(r8) :: delt_rtm_max                 ! RTM max timestep
  real(r8) :: cfl_scale = 0.1_r8           ! cfl scale factor, must be <= 1.0
  real(r8), parameter :: effvel(nt_rtm) = 0.35_r8  ! downstream velocity (m/s)

!glo
  integer , pointer :: dwnstrm_index(:)! downstream index

!gdc
  real(r8), pointer :: ddist(:)        ! downstream dist (m)
  real(r8), pointer :: evel(:,:)       ! effective tracer velocity (m/s)
  real(r8), pointer :: sfluxin(:,:)    ! cell tracer influx (m3/s)
  real(r8), pointer :: fluxout(:,:)    ! cell tracer outlflux (m3/s)
  real(r8), pointer :: totrunin(:,:)   ! cell tracer lnd forcing on rtm grid (mm/s)

!map
  type(mct_sMat)  :: sMat0_l2r
  type(mct_sMat)  :: sMat0_l2r_d
  type(mct_sMatP) :: sMatP_l2r

!EOP
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
    use domainMod    , only : llatlon,ldomain,domain_check, &
                              latlon_type, latlon_init, latlon_clean
    use RtmMapMod    , only : celledge, map_setmapsAR
    use clm_varctl   , only : frivinp_rtm, rtm_nsteps
    use decompMod    , only : get_proc_bounds, ldecomp, gsMap_lnd_gdc2glo
    use clmtype      , only : grlnd
    use ncdio_pio    , only : ncd_pio_openfile, ncd_pio_closefile, ncd_io, &
                              ncd_inqdid, ncd_inqdlen, file_desc_t
    use fileutils    , only : getfil
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
!
! !LOCAL VARIABLES:
!EOP
    real(r8), dimension(4) :: rtmedge = (/ 90._r8, 180._r8, -90._r8, -180._r8 /)  !N,E,S,W edges of rtm grid
    integer  :: ioff(0:8) = (/0,0,1,1,1,0,-1,-1,-1/) !rdirc input to i
    integer  :: joff(0:8) = (/0,1,1,0,-1,-1,-1,0,1/) !rdirc input to j
    integer  :: i,j,k,n,g,n2,nt               ! loop indices
    integer  :: im1,ip1,jm1,jp1,ir,jr,nr      ! neighbor indices
    integer  :: i2,j2                         ! downstream i and j
    real(r8) :: deg2rad                       ! pi/180
    real(r8) :: dx,dx1,dx2,dx3                ! lon dist. betn grid cells (m)
    real(r8) :: dy                            ! lat dist. betn grid cells (m)
    real(r8),allocatable :: tempr(:,:)        ! temporary buffer
    integer ,allocatable :: rdirc(:)          ! temporary buffer
    integer ,allocatable :: iocn(:)           ! downstream ocean cell
    integer ,allocatable :: nocn(:)           ! number of rtm cells in basin
    integer ,allocatable :: pocn(:)           ! pe number assigned to basin
    integer ,allocatable :: nop(:)            ! number of rtm cells on a pe
    integer ,allocatable :: nba(:)            ! number of basins on each pe
    integer ,allocatable :: nrs(:)            ! begr on each pe
    integer ,allocatable :: baslc(:)          ! basin overlap land cell
    integer ,allocatable :: basnp(:)          ! basin pe number
    integer ,allocatable :: basin(:)          ! basin to rtm mapping
    integer  :: nbas                          ! number of basins
    integer  :: nrtm                          ! num of rtm points
    integer  :: baspe                         ! pe with min number of rtm cells
    integer  :: maxrtm                        ! max num of rtms per pe for decomp
    integer  :: minbas,maxbas                 ! used for decomp search
    integer  :: nl,nloops                     ! used for decomp search
    integer  :: ier                           ! error code
    integer  :: mon                           ! month (1, ..., 12)
    integer  :: day                           ! day of month (1, ..., 31)
    integer  :: ncsec                         ! seconds of current date
    integer  :: begg,endg                     ! local start/end gridcell indices
    integer  :: begr,endr,numr                ! tot num of roff pts on all pes
    real(r8) :: dtover,dtovermax              ! ts calc temporaries
    character(len=16), dimension(50) :: river_name
    character(len=30), dimension(50) :: rivstat_name
    real(r8)         , dimension(50) :: rivstat_lon
    real(r8)         , dimension(50) :: rivstat_lat
    type(file_desc_t)                :: ncid  ! netcdf id
    integer  :: dimid                         ! dimension identifier
    integer  :: nroflnd
    integer  :: nrofocn
    integer  :: pid,np,npmin,npmax,npint      ! log loop control
    integer,parameter  :: dbug = 1            ! 0 = none, 1=normal, 2=much, 3=max

    integer lsize,gsize                    ! sizes to initialize GsMap
    integer  :: na,nb,ns                   ! mct sizes
    integer  :: igrow,igcol,iwgt           ! mct field indices
    integer  :: ii,ji,ni,no,gi,go          ! tmps
    real(r8) :: wt                         ! mct wt
    integer ,pointer :: rgdc2glo(:)        ! temporary for initialization
    integer ,pointer :: rglo2gdc(:)        ! temporary for initialization
    real(r8),pointer :: lfield(:)          ! tmp lnd field
    real(r8),pointer :: rfield(:)          ! tmp rtm field
    real(r8),pointer :: glatc(:),glonc(:)  ! global lat/lon
    real(r8),pointer :: gfrac(:)           ! global frac
    real(r8),pointer :: lfrac(:)           ! global land frac
    integer ,pointer :: gmask(:)           ! global mask
    type(latlon_type):: rlatlon            ! rtm grid 
    logical          :: found              ! if variable found on rdirc file
    character(len=256):: locfn             ! local file name
    character(len=32) :: subname = 'Rtmini'! subroutine name

!-----------------------------------------------------------------------

    call t_startf('rtmi_grid')

    !--- Initialize rtm_trstr
    rtm_trstr = trim(rtm_tracers(1))
    do n = 2,nt_rtm
       rtm_trstr = trim(rtm_trstr)//':'//trim(rtm_tracers(n))
    enddo
    if (masterproc) then
       write(iulog,*)'rtm tracers = ',nt_rtm,trim(rtm_trstr)
    end if

    !--- Read in RTM file to get dimension sizes
    call getfil(frivinp_rtm, locfn, 0 )
    if (masterproc) then
       write(iulog,*)'Read in RTM file name: ',trim(frivinp_rtm)
       call shr_sys_flush(iulog)
    endif

    call ncd_pio_openfile (ncid, trim(locfn), 0)

    call ncd_inqdid(ncid,'ni',dimid)
    call ncd_inqdlen(ncid,dimid,rtmlon)
    call ncd_inqdid(ncid,'nj',dimid)
    call ncd_inqdlen(ncid,dimid,rtmlat)

    if (masterproc) then
       write(iulog,*) 'Values for rtmlon/rtmlat: ',rtmlon,rtmlat
       write(iulog,*) 'Successfully read RTM dimensions'
       call shr_sys_flush(iulog)
    endif

    !--- Allocate rtm grid variables
    call latlon_init(rlatlon,rtmlon,rtmlat)

    !--- Allocate inputs and outputs to rtm at 1/2 degree resolution

    allocate (glatc(rtmlon*rtmlat), glonc(rtmlon*rtmlat), &
              stat=ier)
    if (ier /= 0) then
       write(iulog,*)'Rtmgridini: Allocation error for ',&
            'glatc,glonc'
       call endrun
    end if

    allocate (runoff%rlat(rtmlat), runoff%rlon(rtmlon), &
              stat=ier)
    if (ier /= 0) then
       write(iulog,*)'Rtmgridini: Allocation error for ',&
            'rlat,rlon'
       call endrun
    end if

    !--- Allocate temporaries

    allocate(rdirc(rtmlon*rtmlat),stat=ier)
    if (ier /= 0) then
       write(iulog,*)'Rtmgridini: Allocation error for ',&
            'rdirc'
       call endrun
    end if

    !--- Useful constants and initial values

    deg2rad = SHR_CONST_PI / 180._r8

    ! Open and read input data (river direction file)
    rlatlon%edges(1:4) = rtmedge(1:4)
    allocate(tempr(rtmlon,rtmlat))

    call ncd_io(ncid=ncid, varname='RTM_FLOW_DIRECTION', flag='read', data=tempr, readvar=found)
    if ( .not. found ) call endrun( trim(subname)//' ERROR: RTM_FLOW_DIRECTION NOT on rdirc file' )
    do j=1,rtmlat
      do i=1,rtmlon
         n = (j-1)*rtmlon + i
         rdirc(n) = nint(tempr(i,j))
      enddo
    enddo
    call ncd_io(ncid=ncid, varname='xc', flag='read', data=tempr, readvar=found)
    if ( .not. found ) call endrun( trim(subname)//' ERROR: RTM longitudes NOT on rdirc file' )
    do j=1,rtmlat
       do i=1,rtmlon
          n = (j-1)*rtmlon + i
          glonc(n) = tempr(i,j)
       enddo
    enddo
    call ncd_io(ncid=ncid, varname='yc', flag='read', data=tempr, readvar=found)
    if ( .not. found ) call endrun( trim(subname)//' ERROR: RTM latitudes NOT on rdirc file' )
    do j=1,rtmlat
       do i=1,rtmlon
          n = (j-1)*rtmlon + i
          glatc(n) = tempr(i,j)
       enddo
    enddo

    call ncd_pio_closefile(ncid)
    deallocate(tempr)

    if (masterproc) then
       write(iulog,*)'RTM netcdf river direction file successfully read '
       call shr_sys_flush(iulog)
    endif

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

    deallocate(glatc,glonc)

    allocate (dwnstrm_index(rtmlon*rtmlat), &
              gfrac(rtmlon*rtmlat), lfrac(lsmlon*lsmlat), &
              gmask(rtmlon*rtmlat), &
              stat=ier)
    if (ier /= 0) then
       write(iulog,*)'Rtmgridini: Allocation error for ',&
            'dwnstrm_index,gfrac'
       call endrun
    end if

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
          ir = i + ioff(abs(rdirc(n)))
          jr = j + joff(abs(rdirc(n)))
          if (ir < 1     ) ir = ir + rtmlon
          if (ir > rtmlon) ir = ir - rtmlon
          !--- check cross pole flow, etc
          if (jr < 1 .or. jr > rtmlat .or. ir < 1 .or. ir > rtmlon) then
             write(iulog,*) 'Rtmini ERROR ir jr bounds ',i,j,rdirc(n),ir,jr
             call endrun()
          endif
          nr = (jr-1)*rtmlon + ir
          if (n == nr) then
             write(iulog,*) 'Rtmini ERROR dwnstrm_index ',i,j,n,rdirc(n),ir,jr,nr
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

    call t_stopf('rtmi_grid')

    !--- Set sMat0_l2r, full mapping weights for l2r, just on root pe
    !--- for now use lfield to "ignore" non-active land cells in sMat0_l2r
    !--- Later these will be "reduced" to just the useful weights
    !--- Compute gfrac on root pe and bcast 

    !--- Change 9/2007, compute on all pes for decomp based on overlap

    call t_startf('rtmi_setl2r')
    call get_proc_bounds(begg, endg)
    call gather_data_to_master(ldomain%frac,lfrac,grlnd)
    call mpi_bcast(lfrac,size(lfrac),MPI_REAL8,0,mpicom,ier)

    allocate(lfield(lsmlon*lsmlat),rfield(rtmlon*rtmlat))
    lfield = 0._r8
    do n = 1,lsmlon*lsmlat
       if (ldecomp%glo2gdc(n) > 0) lfield(n) = 1._r8
    enddo
    rfield = 1._r8

    call map_setmapsAR(llatlon,rlatlon,sMat0_l2r,fracin=lfield,fracout=rfield)

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

    call t_stopf('rtmi_setl2r')

    !--- Determine rtm ocn/land mask, 0=none, 1=land, 2=ocean outflow, -1=reroute over ocean to ocean outflow points

    call t_startf('rtmi_decomp')

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
          if (rdirc(n) > 0) then
             gmask(n) = 1
          else if (rdirc(n) < 0) then 
             gmask(n) = -1
          end if
       else                      ! n is ocn if no dwnstrm_index and some frac
          if (gfrac(n)>0._r8) gmask(n) = 2
       end if
    enddo

    deallocate(rdirc)
    deallocate(gfrac,lfrac)

   !--- Compute river basins, actually compute ocean outlet gridcell
   !--- iocn = final downstream cell, index is global 1d ocean gridcell
   !--- nocn = number of source gridcells for ocean gridcell

    allocate(iocn(rtmlon*rtmlat),nocn(rtmlon*rtmlat),stat=ier)
    if (ier /= 0) then
       write(iulog,*)'Rtmgridini: Allocation error for ',&
            'iocn,nocn'
       call endrun
    end if

    call t_startf('rtmi_dec_basins')
    iocn = 0
    nocn = 0
    do nr=1,rtmlon*rtmlat
       n = nr
       if (abs(gmask(n)) == 1) then    ! land
          g = 0
          do while (abs(gmask(n)) == 1 .and. g < rtmlon*rtmlat)  ! follow downstream
             n = dwnstrm_index(n)
             g = g + 1
          end do
          if (gmask(n) == 2) then  ! found ocean outlet
             iocn(nr) = n                 ! set ocean outlet or nr to n
             nocn(n) = nocn(n) + 1        ! one more land cell for n
          elseif (abs(gmask(n)) == 1) then  ! no ocean outlet, warn user, ignore cell
             write(iulog,*) 'rtmini WARNING no downstream ocean cell - IGNORED', &
               g,nr,gmask(nr),dwnstrm_index(nr), &
               n,gmask(n),dwnstrm_index(n)
          else 
             write(iulog,*) 'rtmini ERROR downstream cell is non-ocean,non-land', &
               g,nr,gmask(nr),dwnstrm_index(nr), &
               n,gmask(n),dwnstrm_index(n)
             call endrun()
          endif
       elseif (gmask(n) == 2) then  ! ocean, give to self
          iocn(nr) = n
          nocn(n) = nocn(n) + 1
       endif
    enddo
    call t_stopf('rtmi_dec_basins')

    !--- Now allocate those basins to pes
    !--- pocn is the pe that gets the basin associated with ocean outlet nr
    !--- nop is a running count of the number of rtm cells/pe 

    call t_startf('rtmi_dec_distr')

    nbas = 0
    nrtm = 0
    do nr=1,rtmlon*rtmlat
       if (nocn(nr) > 0) then
          nbas = nbas + 1
          nrtm = nrtm + nocn(nr)
       endif
    enddo

    allocate(pocn(rtmlon*rtmlat), &
             nop(0:npes-1),nba(0:npes-1),nrs(0:npes-1), &
             rglo2gdc(rtmlon*rtmlat),runoff%num_rtm(0:npes-1))

    nop = 0
    nba = 0
    nrs = 0
    pocn = -99
    rglo2gdc = 0
    baspe = 0
    maxrtm = int(float(nrtm)/float(npes)*0.445) + 1
    call shr_sys_flush(iulog)
    nloops = 3
    minbas = nrtm
    do nl=1,nloops
       maxbas = minbas - 1
       minbas = maxval(nocn)/(2**nl)
       if (nl == nloops) minbas = min(minbas,1)
    do nr=1,rtmlon*rtmlat
!       if (nocn(nr) /= 0) then
       if (nocn(nr) > 0 .and. nocn(nr) >= minbas .and. nocn(nr) <= maxbas) then
! Decomp options
!   find min pe (implemented but scales poorly)
!   use increasing thresholds (implemented, ok load balance for l2r or calc)
!   distribute basins using above methods but work from max to min basin size
!   distribute basins to minimize l2r time, basins put on pes associated 
!      with lnd forcing, need to know l2r map and lnd decomp
!
!--------------
! find min pe
!             baspe = 0
!             do n = 1,npes-1
!                if (nop(n) < nop(baspe)) baspe = n
!             enddo
!--------------
! find next pe below maxrtm threshhold and increment
             do while (nop(baspe) > maxrtm)
                baspe = baspe + 1
                if (baspe > npes-1) then
                   baspe = 0
                   maxrtm = max(maxrtm*1.5, maxrtm+1.0)   ! 3 loop, .445 and 1.5 chosen carefully
                endif
             enddo
!--------------
          if (baspe > npes-1 .or. baspe < 0) then
             write(iulog,*) 'error in decomp for rtm ',nr,npes,baspe
             call endrun()
          endif
          nop(baspe) = nop(baspe) + nocn(nr)
          nba(baspe) = nba(baspe) + 1
          pocn(nr) = baspe
       endif
    enddo ! nr
    enddo ! nl

    ! set pocn for land cells, was set for ocean above
    do nr=1,rtmlon*rtmlat
       if (iocn(nr) > 0) then
          pocn(nr) = pocn(iocn(nr))
          if (pocn(nr) < 0 .or. pocn(nr) > npes-1) then
             write(iulog,*) 'Rtmini ERROR pocn lnd setting ',nr,iocn(nr),iocn(iocn(nr)),pocn(iocn(nr)),pocn(nr),npes
             call endrun()
          endif
       endif
    enddo

    if (masterproc) write(iulog,*) 'rtm cells and basins total  = ',nrtm,nbas
    if (masterproc) write(iulog,*) 'rtm cells per basin avg/max = ',nrtm/nbas,maxval(nocn)
    if (masterproc) write(iulog,*) 'rtm cells per pe    min/max = ',minval(nop),maxval(nop)
    if (masterproc) write(iulog,*) 'basins    per pe    min/max = ',minval(nba),maxval(nba)

    !--- Count and distribute cells to rglo2gdc

    runoff%numr   = 0
    runoff%numro  = 0
    runoff%numrl  = 0
    runoff%lnumr  = 0
    runoff%lnumro = 0
    runoff%lnumrl = 0
    runoff%num_rtm = 0

    do n = 0,npes-1
       if (iam == n) then
          runoff%begr  = runoff%numr  + 1
          runoff%begrl = runoff%numrl + 1
          runoff%begro = runoff%numro + 1
       endif

       runoff%num_rtm(n) = runoff%num_rtm(n) + nop(n)
       runoff%numr  = runoff%numr  + nop(n)
       runoff%numro = runoff%numro + nba(n)
       runoff%numrl = runoff%numrl + nop(n) - nba(n)

       if (iam == n) then
          runoff%lnumr  = runoff%lnumr  + nop(n)
          runoff%lnumro = runoff%lnumro + nba(n)
          runoff%lnumrl = runoff%lnumrl + nop(n) - nba(n)
          runoff%endr  = runoff%begr  + runoff%lnumr  - 1
          runoff%endro = runoff%begro + runoff%lnumro - 1
          runoff%endrl = runoff%begrl + runoff%lnumrl - 1
       endif
    enddo

    ! nrs is begr on each pe
    nrs(0) = 1
    do n = 1,npes-1
       nrs(n) = nrs(n-1) + nop(n-1)
    enddo

    ! reuse nba for nop-like counter here
    ! pocn -99 is unused cell
    nba = 0
    do nr = 1,rtmlon*rtmlat
       if (pocn(nr) >= 0) then
          rglo2gdc(nr) = nrs(pocn(nr)) + nba(pocn(nr))
          nba(pocn(nr)) = nba(pocn(nr)) + 1          
       endif
    enddo
    do n = 0,npes-1
       if (nba(n) /= nop(n)) then
          write(iulog,*) 'Rtmini ERROR rtm cell count ',n,nba(n),nop(n)
          call endrun()
       endif
    enddo

    deallocate(nop,nba,nrs)
    deallocate(iocn,nocn)
    deallocate(pocn)
    call t_stopf('rtmi_dec_distr')

    !--- set some local values

    nroflnd = runoff%numrl
    nrofocn = runoff%numro
    numr = nroflnd + nrofocn
    begr = runoff%begr
    endr = runoff%endr

    call t_stopf('rtmi_decomp')

    !--- Write per-processor runoff bounds depending on dbug level

    call t_startf('rtmi_print')

    call shr_sys_flush(iulog)
    if (masterproc) then
       write(iulog,*) 'total runoff cells numr = ',runoff%numr, &
          'numrl = ',runoff%numrl,'numro = ',runoff%numro
    endif
    call shr_sys_flush(iulog)
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
          write(iulog,*) 'rtm decomp info',' proc = ',iam, &
             ' begr = ',runoff%begr,' endr = ',runoff%endr, &
             ' numr = ',runoff%lnumr
          write(iulog,*) '               ',' proc = ',iam, &
             ' begrl= ',runoff%begrl,' endrl= ',runoff%endrl, &
             ' numrl= ',runoff%lnumrl
          write(iulog,*) '               ',' proc = ',iam, &
             ' begro= ',runoff%begro,' endro= ',runoff%endro, &
             ' numro= ',runoff%lnumro
       endif
       call shr_sys_flush(iulog)
       call mpi_barrier(mpicom,ier)
    enddo

    call t_stopf('rtmi_print')

    !--- allocate runoff variables

    call t_startf('rtmi_vars')

    allocate(runoff%runoff(begr:endr,nt_rtm),runoff%dvolrdt(begr:endr,nt_rtm), &
             runoff%runofflnd(begr:endr,nt_rtm),runoff%dvolrdtlnd(begr:endr,nt_rtm), &
             runoff%runoffocn(begr:endr,nt_rtm),runoff%dvolrdtocn(begr:endr,nt_rtm), &
             runoff%area(begr:endr), &
             runoff%volr(begr:endr,nt_rtm), runoff%volrlnd(begr:endr,nt_rtm), &
             runoff%lonc(begr:endr),  runoff%latc(begr:endr),  &
             runoff%dsi(begr:endr), stat=ier)
    if (ier /= 0) then
       write(iulog,*)'Rtmini ERROR allocation of runoff%runoff'
       call endrun
    end if

    allocate(runoff%runofflnd_nt1(begr:endr),runoff%runofflnd_nt2(begr:endr), &
             runoff%runoffocn_nt1(begr:endr),runoff%runoffocn_nt2(begr:endr), &
             runoff%volr_nt1(begr:endr), runoff%volr_nt2(begr:endr), &
             runoff%dvolrdtlnd_nt1(begr:endr),runoff%dvolrdtlnd_nt2(begr:endr), &
             runoff%dvolrdtocn_nt1(begr:endr),runoff%dvolrdtocn_nt2(begr:endr), &
             stat=ier)
    if (ier /= 0) then
       write(iulog,*)'Rtmini ERROR allocation of runoff%runoff_nt'
       call endrun
    end if

    allocate(rgdc2glo(numr), runoff%mask(numr), stat=ier)
    if (ier /= 0) then
       write(iulog,*)'Rtmini ERROR allocation of runoff%gcd2glo'
       call endrun
    end if

    !--- Allocate rtm flux variables

    allocate (fluxout (begr:endr,nt_rtm), &
              ddist   (begr:endr), &
              totrunin(begr:endr,nt_rtm), &
              evel    (begr:endr,nt_rtm), &
              sfluxin (begr:endr,nt_rtm),  stat=ier)
    if (ier /= 0) then
       write(iulog,*)'Rtmgridini: Allocation error for ',&
            'volr, fluxout, ddist'
       call endrun
    end if
    fluxout = 0._r8
    ddist = 0._r8
    sfluxin = 0._r8
    do nt = 1,nt_rtm
    do nr = begr,endr
      evel(nr,nt) = effvel(nt)
    enddo
    enddo

    !--- Initialize runoff data

    numr = 0
    do j = 1,rtmlat
    do i = 1,rtmlon
       n = (j-1)*rtmlon + i
       nr = rglo2gdc(n)
       if (nr /= 0) then
          numr = numr + 1
          rgdc2glo(nr) = n         
          runoff%mask(nr) = gmask(n)   ! global for hist file
       endif
    enddo
    enddo

    if (numr /= runoff%numr) then
       write(iulog,*) 'Rtmini ERROR numr numr ',numr,runoff%numr
       call endrun()
    endif

    runoff%runoff = 0._r8
    runoff%runofflnd = spval
    runoff%runoffocn = spval
    runoff%dvolrdt = 0._r8
    runoff%dvolrdtlnd = spval
    runoff%dvolrdtocn = spval
    runoff%volr = 0._r8
    runoff%volrlnd = spval

    do nr = begr,endr
       n = rgdc2glo(nr)
       i = mod(n-1,rtmlon) + 1
       j = (n-1)/rtmlon + 1
       if (n <= 0 .or. n > rtmlon*rtmlat) then
          write(iulog,*) 'Rtmini ERROR gdc2glo ',nr,rgdc2glo(nr)
          call endrun()
       endif
       runoff%lonc(nr) = runoff%rlon(i)
       runoff%latc(nr) = runoff%rlat(j)

       if (runoff%mask(nr) == 1) then
          do nt = 1,nt_rtm
             runoff%runofflnd(nr,nt) = runoff%runoff(nr,nt)
             runoff%dvolrdtlnd(nr,nt)= runoff%dvolrdt(nr,nt)
             runoff%volrlnd(nr,nt)= runoff%volr(nr,nt)
          enddo
       elseif (runoff%mask(nr) == 2) then
          do nt = 1,nt_rtm
             runoff%runoffocn(nr,nt) = runoff%runoff(nr,nt)
             runoff%dvolrdtocn(nr,nt)= runoff%dvolrdt(nr,nt)
          enddo
       endif

       !tcx testing
       dx = (rlatlon%lone(i) - rlatlon%lonw(i)) * deg2rad
       dy = sin(rlatlon%latn(j)*deg2rad) - sin(rlatlon%lats(j)*deg2rad)
       runoff%area(nr) = 1.e6_r8 * dx*dy*re*re
       if (dwnstrm_index(n) == 0) then
          runoff%dsi(nr) = 0
       else
          if (rglo2gdc(dwnstrm_index(n)) == 0) then
             write(iulog,*) 'Rtmini ERROR glo2gdc dwnstrm ',nr,n,dwnstrm_index(n),rglo2gdc(dwnstrm_index(n))
             call endrun()
          endif
          runoff%dsi(nr) = rglo2gdc(dwnstrm_index(n))
       endif
    enddo

    !--- Determine downstream distance - instead of reading a distance file
    !--- calculate the downstream distance

    do nr=begr,endr
       g = runoff%dsi(nr)
       if (g == 0) then
          ddist(nr) = 0._r8
       elseif (g < begr .or. g > endr) then
          write(iulog,*) 'Rtmini: error in ddist calc ',nr,g,begr,endr
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

    dtover = 0._r8
    dtovermax = 0._r8
    do nt=1,nt_rtm
    do nr=begr,endr
       if (ddist(nr) /= 0._r8) then
          dtover = evel(nr,nt)/ddist(nr)
       else
          dtover = 0._r8
       endif
       dtovermax = max(dtovermax,dtover)
    enddo
    enddo
    dtover = dtovermax
    call mpi_allreduce(dtover,dtovermax,1,MPI_REAL8,MPI_MAX,mpicom,ier)
    if (dtovermax > 0._r8) then
       delt_rtm_max = (1.0_r8/dtovermax)*cfl_scale
    else
       write(iulog,*) 'rtmini error in delt_rtm_max ',delt_rtm_max,dtover
       call endrun
    endif
    if (masterproc) write(iulog,*) 'rtm max timestep = ',delt_rtm_max,' (sec) for cfl_scale = ',cfl_scale

    !--- Allocate and initialize rtm input fields on clm decomp

    call get_proc_bounds(begg, endg)
    allocate (rtmin_avg(begg:endg,nt_rtm), rtmin_acc(begg:endg,nt_rtm), stat=ier)
    if (ier /= 0) then
       write(iulog,*)'Rtmlandini: Allocation error for rtmin, rtmin_avg, rtmin_acc'
       call endrun
    end if
    rtmin_avg = 0._r8
    rtmin_acc = 0._r8

    !--- clean up temporaries

    deallocate(dwnstrm_index,gmask)
    call t_stopf('rtmi_vars')

   !--- initialization rtm gsmap

    call t_startf('rtmi_mctdata')

    allocate(runoff%gindex(begr:endr))
    do n = begr,endr
       runoff%gindex(n) = rgdc2glo(n)
    enddo
    lsize = endr-begr+1
    gsize = rtmlon * rtmlat
    call mct_gsMap_init( gsMap_rtm_gdc2glo, runoff%gindex, mpicom, comp_id, lsize, gsize )

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
          if (ldecomp%glo2gdc(ni) > 0 .and. rglo2gdc(no) > 0) then
             ns = ns + 1
          endif
       enddo

       call mct_sMat_init(sMat0_l2r_d, nb, na, ns)

       ns = 0
       do n = 1,mct_sMat_lsize(sMat0_l2r)
          ni = sMat0_l2r%data%iAttr(igcol,n)
          no = sMat0_l2r%data%iAttr(igrow,n)
          if (ldecomp%glo2gdc(ni) > 0 .and. rglo2gdc(no) > 0) then
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

   deallocate(rgdc2glo,rglo2gdc)
   call latlon_clean(rlatlon)

   !--- clean up the root sMat0 datatypes

   call mct_sMat_clean(sMat0_l2r)
   if (masterproc) then
      call mct_sMat_clean(sMat0_l2r_d)
      write(iulog,*) 'Rtmini complete'
   endif

   !--- update rtm history fields

   call rtm_sethist()

   call t_stopf('rtmi_mctdata')

  end subroutine Rtmini

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: RtmMap
!
! !INTERFACE:
  subroutine RtmMapl2r()
!
! !DESCRIPTION:
! Interface with RTM river routing model.
!
! !USES:
    use decompMod      , only : get_proc_bounds
    use decompMod      , only : gsMap_lnd_gdc2glo
    use domainMod      , only : ldomain
!
! !ARGUMENTS:
    implicit none
!
! !CALLED FROM:
! subroutine clm_driver2
!
! !REVISION HISTORY:
! Author: Sam Levis
!
!EOP
!
    integer  :: i,j,n,n2,nr,ns,nt          ! indices
    integer  :: begg,endg
    logical  :: usevector = .false.
    real(r8) :: suml(nt_rtm),sumr(nt_rtm),sumlt(nt_rtm),sumrt(nt_rtm)   ! water diagnostics
    integer  :: ier
    type(mct_aVect)   :: aV_lndr,aV_rtmr
    integer,parameter :: dbug = 1
!-----------------------------------------------------------------------

    ! Determine RTM inputs on land model grid

    call t_startf('clmrtm_l2r')

    ! Map RTM inputs from land model grid to RTM grid (1/2 degree resolution)

    ns = mct_gsMap_lsize(gsmap_lnd_gdc2glo, mpicom)
    call mct_aVect_init(aV_lndr,rlist=trim(rtm_trstr),lsize=ns)

    ns = mct_gsMap_lsize(gsMap_rtm_gdc2glo, mpicom)
    call mct_aVect_init(aV_rtmr,rlist=trim(rtm_trstr),lsize=ns)

    suml = 0._r8
    sumr = 0._r8
    call get_proc_bounds(begg, endg)
    do n = begg,endg
       do nt = 1,nt_rtm
          n2 = n-begg+1
          av_lndr%rAttr(nt,n2) = rtmin_avg(n,nt)*ldomain%frac(n)
          suml(nt) = suml(nt) + av_lndr%rAttr(nt,n2)*ldomain%area(n)
       enddo
    enddo
    
    call mct_Smat_AvMult(av_lndr, sMatP_l2r, av_rtmr, vector=usevector)
 
    do n = runoff%begr,runoff%endr
       do nt = 1,nt_rtm
          n2 = n-runoff%begr+1
          totrunin(n,nt) = av_rtmr%rAttr(nt,n2)
          sumr(nt) = sumr(nt) + totrunin(n,nt)*runoff%area(n)*1.0e-6_r8   ! area m2 to km2
       enddo
    enddo
    
    if (dbug > 1) then
       call mpi_reduce(suml, sumlt, nt_rtm, MPI_REAL8, MPI_SUM, 0, mpicom, ier)
       call mpi_reduce(sumr, sumrt, nt_rtm, MPI_REAL8, MPI_SUM, 0, mpicom, ier)
       if (masterproc) then
          do nt = 1,nt_rtm
             if (abs(sumlt(nt)+sumrt(nt)) > 0.0_r8) then
                if (abs(sumlt(nt) - sumrt(nt))/(sumlt(nt)+sumrt(nt)) > 1.0e-6) then
                   write(iulog,*) 'WARNING: l2r water not conserved ',nt,sumlt(nt),sumrt(nt)
                endif
             endif
          enddo
       endif
    endif
    
    call mct_aVect_clean(aV_lndr)
    call mct_aVect_clean(aV_rtmr)

    call t_stopf('clmrtm_l2r')
    
  end subroutine RtmMapl2r

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: RtmInput
!
! !INTERFACE:
  subroutine RtmInput(do_rtm, begg, endg, qflx_runoffg, qflx_snwcp_iceg)
!
! !DESCRIPTION:
! Update RTM inputs.
!
! !USES:
    use domainMod       , only : ldomain
    use clm_varctl      , only : rtm_nsteps
    use clm_time_manager, only : get_step_size, get_nstep
!
! !ARGUMENTS:
    implicit none
    logical , intent(out) :: do_rtm
    integer , intent(in)  :: begg, endg                 ! per-proc gridcell ending gridcell indices
    real(r8), intent(in)  :: qflx_runoffg(begg:endg)    ! total runoff (mm H2O /s)
    real(r8), intent(in)  :: qflx_snwcp_iceg(begg:endg) ! excess snowfall due to snow capping (mm H2O /s)
!
! !CALLED FROM:
! subroutine rtmMap
!
! !REVISION HISTORY:
! Author: Sam Levis
!
! !LOCAL VARIABLES:
!EOP
!
    integer  :: g,nt                           ! indices
    integer  :: nliq,nfrz                      ! field indices
    integer  :: nstep                          ! time step index
    logical  :: first_time = .true.
!-----------------------------------------------------------------------

    call t_barrierf('sync_clmrtm', mpicom)

    ! Make gridded representation of runoff from clm for tracers

    nliq = 0
    nfrz = 0
    do nt = 1,nt_rtm
       if (trim(rtm_tracers(nt)) == 'LIQ') then
          nliq = nt
       endif
       if (trim(rtm_tracers(nt)) == 'ICE') then
          nfrz = nt
       endif
    enddo
    if (nliq == 0 .or. nfrz == 0) then
       write(iulog,*)'RtmInput: ERROR in rtm_tracers LIQ ICE ',nliq,nfrz,nt_rtm,rtm_tracers
       call endrun()
    endif

    do g = begg, endg
       rtmin_acc(g,nliq) = rtmin_acc(g,nliq)+qflx_runoffg(g)
       rtmin_acc(g,nfrz) = rtmin_acc(g,nfrz)+qflx_snwcp_iceg(g)
    enddo

    ncount_rtm = ncount_rtm + 1

    if (first_time) then	
       if (rtm_nsteps < 1) then
          write(iulog,*) 'rtm ERROR in rtm_nsteps',rtm_nsteps
          call endrun()
       endif
       delt_rtm = rtm_nsteps * get_step_size()
       if (masterproc) write(iulog,*) 'rtm act timestep ~ ',delt_rtm
       first_time = .false.
    end if

    nstep = get_nstep()
    if (mod(nstep,rtm_nsteps)==0 .and. nstep>1) then
       if (ncount_rtm*get_step_size() /= delt_rtm) then
          if (masterproc) write(iulog,*) 'RtmInput timestep out of sync ',&
               delt_rtm,ncount_rtm*get_step_size()
          delt_rtm = ncount_rtm*get_step_size()
       endif
       do nt = 1,nt_rtm
       do g = begg,endg
          rtmin_avg(g,nt) = rtmin_acc(g,nt)*ldomain%asca(g)/(ncount_rtm*1.0_r8)
          rtmin_acc(g,nt) = 0._r8
       end do
       end do
       ncount_rtm = 0                          !reset counter to 0
       do_rtm = .true.
    else
       do_rtm = .false.
    endif

  end subroutine RtmInput

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
! subroutine RtmMap in this module
!
! !REVISION HISTORY:
! Author: Sam Levis
!
!
! !LOCAL VARIABLES:
!EOP
    integer  :: i, j, n, ns, nt             !loop indices
    integer  :: ir,jr,nr                    !neighbor indices
    real(r8) :: dvolrdt                     !change in storage in discharge units (m3/s)
    real(r8) :: sumfin(nt_rtm),sumfex(nt_rtm)
    real(r8) :: sumrin(nt_rtm),sumdvt(nt_rtm)
    real(r8) :: sum1,sum2
    integer  :: nsub                        !subcyling for cfl
    integer, save :: nsub_save              !previous nsub
    real(r8) :: delt                        !delt associated with subcycling
    real(r8), save :: delt_save             !previous delt
    integer,parameter :: dbug = 1           !local debug flag
!-----------------------------------------------------------------------

    nsub = int(delt_rtm/delt_rtm_max) + 1
    delt = delt_rtm/float(nsub)

    if (delt /= delt_save) then
       if (masterproc) write(iulog,*) 'rtm delt update from/to',delt_save,delt,nsub_save,nsub
    endif

    nsub_save = nsub
    delt_save = delt

    sumfin = 0._r8
    sumfex = 0._r8
    sumrin = 0._r8
    sumdvt = 0._r8

    runoff%runoff = 0._r8
    runoff%runofflnd = spval
    runoff%runoffocn = spval
    runoff%dvolrdt = 0._r8
    runoff%dvolrdtlnd = spval
    runoff%dvolrdtocn = spval

    !--- subcycling ---
    do ns = 1,nsub

       sfluxin = 0._r8
       do n = runoff%begr,runoff%endr
          nr = runoff%dsi(n)
          if (abs(runoff%mask(n)) == 1) then
             if (nr < runoff%begr .or. nr > runoff%endr) then
                write(iulog,*) 'Rtm ERROR: non-local communication ',n,nr
                call endrun()
             endif
             do nt = 1,nt_rtm
                sfluxin(nr,nt) = sfluxin(nr,nt) + fluxout(n,nt)
             enddo
          endif
       enddo

       if (dbug > 1) then
          do nt = 1,nt_rtm
             sum1 = 0._r8
             sum2 = 0._r8
             do n = runoff%begr,runoff%endr
                sum1 = sum1 + sfluxin(n,nt)
                sum2 = sum2 + fluxout(n,nt)
             enddo
             if (abs(sum1+sum2) > 0.0_r8) then
             if (abs(sum1-sum2)/(sum1+sum2) > 1.0e-12) then
                write(iulog,*) 'RTM Warning: fluxin = ',sum1,&
                     ' not equal to fluxout = ',sum2,' for tracer ',nt
             endif
             endif
          enddo
       endif

       do nt = 1,nt_rtm
       do n = runoff%begr,runoff%endr
          dvolrdt = sfluxin(n,nt) + 0.001_r8*totrunin(n,nt)*runoff%area(n) - fluxout(n,nt)

          if (dbug > 1) then
             sumfin(nt) = sumfin(nt) + sfluxin(n,nt)
             sumfex(nt) = sumfex(nt) + fluxout(n,nt)
             sumrin(nt) = sumrin(nt) + 0.001_r8*totrunin(n,nt)*runoff%area(n)
             sumdvt(nt) = sumdvt(nt) + dvolrdt
          endif

          if (abs(runoff%mask(n)) == 1) then         ! land points
             runoff%volr(n,nt)     = runoff%volr(n,nt) + dvolrdt*delt
             fluxout(n,nt)  = runoff%volr(n,nt) * evel(n,nt)/ddist(n)
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
             runoff%volr(n,nt)     = 0._r8
             fluxout(n,nt) = 0._r8
          endif

          if (abs(runoff%mask(n)) == 1) then
             runoff%runoff(n,nt) = runoff%runoff(n,nt) + fluxout(n,nt)
          elseif (runoff%mask(n) == 2) then
             runoff%runoff(n,nt) = runoff%runoff(n,nt) + dvolrdt
          endif
          ! Convert local dvolrdt (in m3/s) to output dvolrdt (in mm/s)
          runoff%dvolrdt(n,nt) = runoff%dvolrdt(n,nt) + 1000._r8*dvolrdt/runoff%area(n)

       enddo
       enddo

    enddo

    ! average fluxes over subcycling
    runoff%runoff = runoff%runoff / float(nsub)
    runoff%dvolrdt = runoff%dvolrdt / float(nsub)

    do n = runoff%begr,runoff%endr
       if (runoff%mask(n) == 1) then
          do nt = 1,nt_rtm
             runoff%volrlnd(n,nt)= runoff%volr(n,nt)
             runoff%runofflnd(n,nt) = runoff%runoff(n,nt)
             runoff%dvolrdtlnd(n,nt)= runoff%dvolrdt(n,nt)
          enddo
       elseif (runoff%mask(n) == 2) then
          do nt = 1,nt_rtm
             runoff%runoffocn(n,nt) = runoff%runoff(n,nt)
             runoff%dvolrdtocn(n,nt)= runoff%dvolrdt(n,nt)
          enddo
       endif
    enddo

    call rtm_sethist()

    ! Global water balance calculation and error check

    if (dbug > 1) then
       do nt = 1,nt_rtm
       if (abs(sumdvt(nt)+sumrin(nt)) > 0.0_r8) then
       if (abs((sumdvt(nt)-sumrin(nt))/(sumdvt(nt)+sumrin(nt))) > 1.0e-6) then
          write(iulog,*) 'RTM Warning: water balance nt,dvt,rin,fin,fex = ', &
             nt,sumdvt(nt),sumrin(nt),sumfin(nt),sumfex(nt)
!          call endrun
       endif
       endif
       enddo
    endif

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
    use ncdio_pio
    use decompMod       , only : get_proc_bounds
    use clm_time_manager, only : is_restart
!
! !ARGUMENTS:
    implicit none
    type(file_desc_t), intent(inout)  :: ncid ! netcdf id
    character(len=*) , intent(in) :: flag   ! 'read' or 'write'
!
! !CALLED FROM:
! subroutine restart in module restFileMod
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!
! !OTHER LOCAL VARIABLES:
!EOP
    logical :: readvar          ! determine if variable is on initial file
    integer :: nt,nv,n          ! indices
    integer :: begg,endg        ! start end indices
    real(r8) , pointer :: dfld(:) ! temporary array
    character(len=32)  :: vname,uname
    character(len=128) :: lname
!-----------------------------------------------------------------------

    do nv = 1,4
    do nt = 1,nt_rtm

       if (nv == 1) then
          vname = 'RTM_VOLR_'//trim(rtm_tracers(nt))
          lname = 'water volume in cell (volr)'
          uname = 'm3'
          dfld  => runoff%volr(:,nt)
       elseif (nv == 2) then
          vname = 'RTM_FLUXOUT_'//trim(rtm_tracers(nt))
          lname = 'water fluxout in cell (fluxout)'
          uname = 'm3/s'
          dfld  => fluxout(:,nt)
       elseif (nv == 3) then
          vname = 'RTM_RUNOFF_'//trim(rtm_tracers(nt))
          lname = 'runoff (runoff)'
          uname = 'm3/s'
          dfld  => runoff%runoff(:,nt)
       elseif (nv == 4) then
          vname = 'RTM_DVOLRDT_'//trim(rtm_tracers(nt))
          lname = 'water volume change in cell (dvolrdt)'
          uname = 'mm/s'
          dfld  => runoff%dvolrdt(:,nt)
       else
          write(iulog,*) 'Rtm ERROR: illegal nv value a ',nv
          call endrun()
       endif

       if (flag == 'define') then
          call ncd_defvar(ncid=ncid, varname=trim(vname), &
               xtype=ncd_double,  dim1name='rtmlon', dim2name='rtmlat', &
               long_name=trim(lname), units=trim(uname))
       else if (flag == 'read' .or. flag == 'write') then
          call ncd_io(varname=trim(vname), data=dfld, dim1name='allrof', &
               ncid=ncid, flag=flag, readvar=readvar)
          if (flag=='read' .and. .not. readvar) then
             if (is_restart()) then
                !tcx this is for backward compatability, will not be bfb, endrun should be used
                write(iulog,*) 'Rtm ERROR: data not found on restart, set to zero ',trim(vname)
                dfld = 0._r8
                ! call endrun()
             else
                dfld = 0._r8
             end if
          end if
       end if

    enddo
    enddo

    do nt = 1,nt_rtm
       vname = 'RTM_INPUT_'//trim(rtm_tracers(nt))
       lname = 'average input on clm grid (rtmin_acc)'
       uname = 'mm/s'
       dfld => rtmin_acc(:,nt)

       if (flag == 'define') then
          call ncd_defvar(ncid=ncid, varname=trim(vname), &
               xtype=ncd_double,  dim1name='gridcell', &
               long_name=trim(lname), units=trim(uname))
       else if (flag == 'read' .or. flag == 'write') then
          call ncd_io(varname=trim(vname), data=dfld, dim1name='gridcell', &
               ncid=ncid, flag=flag, readvar=readvar)
          if (flag=='read' .and. .not. readvar) then
             if (is_restart()) then
                !tcx this is for backward compatability, will not be bfb, endrun should be used
                write(iulog,*) 'Rtm ERROR: data not found on restart, set to zero ',trim(vname)
                dfld = 0._r8
                ! call endrun()
             else
                dfld = 0._r8
             end if
          end if
       end if

    enddo

    ! counter for rtm averaged input (on clm grid)
    vname = 'RTM_NCOUNT'
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname=trim(vname), xtype=ncd_int,  &
            long_name='counter for RTM averaged input on CLM grid', units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname=trim(vname), data=ncount_rtm, &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) then
             write(iulog,*) 'Rtm ERROR: data not found on restart RTM_NCOUNT'
             call endrun()
          else
             ncount_rtm = 0
          end if
       end if
    end if

    if (flag == 'read') then
       do n = runoff%begr,runoff%endr
          if (runoff%mask(n) == 1) then
             do nt = 1,nt_rtm
                runoff%runofflnd(n,nt) = runoff%runoff(n,nt)
                runoff%dvolrdtlnd(n,nt)= runoff%dvolrdt(n,nt)
                runoff%volrlnd(n,nt)= runoff%volr(n,nt)
             enddo
          elseif (runoff%mask(n) == 2) then
             do nt = 1,nt_rtm
                runoff%runoffocn(n,nt) = runoff%runoff(n,nt)
                runoff%dvolrdtocn(n,nt)= runoff%dvolrdt(n,nt)
             enddo
          endif
       enddo
       call rtm_sethist()
    endif

  end subroutine RtmRest

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: rtm_sethist
!
! !INTERFACE:
  subroutine rtm_sethist()
!
! !DESCRIPTION:
! set rtm history fields
!
! !USES:
!
! !ARGUMENTS:
    implicit none
!
! !REVISION HISTORY:
! Author: T Craig
!
!
! !OTHER LOCAL VARIABLES:
!EOP

    runoff%runofflnd_nt1(:) = runoff%runofflnd(:,1)
    runoff%runofflnd_nt2(:) = runoff%runofflnd(:,2)
    runoff%runoffocn_nt1(:) = runoff%runoffocn(:,1)
    runoff%runoffocn_nt2(:) = runoff%runoffocn(:,2)
    runoff%dvolrdtlnd_nt1(:) = runoff%dvolrdtlnd(:,1)
    runoff%dvolrdtlnd_nt2(:) = runoff%dvolrdtlnd(:,2)
    runoff%dvolrdtocn_nt1(:) = runoff%dvolrdtocn(:,1)
    runoff%dvolrdtocn_nt2(:) = runoff%dvolrdtocn(:,2)
    runoff%volr_nt1(:) = runoff%volrlnd(:,1)
    runoff%volr_nt2(:) = runoff%volrlnd(:,2)

  end subroutine rtm_sethist

#endif

end module RtmMod
