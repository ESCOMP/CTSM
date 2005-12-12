#include <misc.h>
#include <preproc.h>

module clm_camMod

#if (defined COUP_CAM)

!---------------------------------------------------------------------------
!BOP
!
! !MODULE: clm_camMod
!
! !DESCRIPTION:
! This module provides the interface between the atmosphere land modules.
! If running as part of cam, the land surface model must use the same
! grid as the cam. The land surface model calculates its own net solar
! radiation and net longwave radiation at the surface. The net longwave
! radiation at the surface will differ somewhat from that calculated in the
! atmospheric model because the atm model will use the upward longwave flux
! (or radiative temperature) from the previous time step whereas the land
! surface model uses the flux for the current time step. The net solar
! radiation should equal that calculated in the atmospheric model. If not,
! there is a problem in how the models are coupled.
!
! !USES:
  use shr_kind_mod , only : r8 => shr_kind_r8
  use ppgrid       , only : begchunk, endchunk
  use decompMod    , only : get_nclumps, get_clump_owner_id, &
                            get_clump_ncells_id, get_clump_coord_id, &
                            get_clump_gcell_info
#if (defined SPMD)
  use mpishorthand , only : mpir8, mpicom
  use spmd_utils   , only : npes, iam
#else
  use spmdMod      , only : npes, iam
#endif
  use abortutils   , only : endrun
!
! !PUBLIC TYPES:
  implicit none
  public clm_camInit                   ! Initialization
  public clm_camRun                    ! run method
  public clm_camFinal                  ! Finalization method
  SAVE
  private                              ! By default make data private
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
!
! !PRIVATe MEMBER FUNCTIONS:
  private clm_camCheckGrid               ! check consistency of cam/clm grid
  private lp_coupling_init               ! initialize clump<-->chunk mapping
  private lp_coupling_finalize           ! destroy clump<-->chunk mapping
  private lp_all2all_clump_to_chunk_init ! communicate fluxes from lnd to atm
  private lp_all2all_clump_to_chunk      ! communicate fluxes from lnd to atm
  private lp_all2all_chunk_to_clump      ! communicate fluxes from atm to lnd
!
! !PRIVATE TYPES:
   type clump2chunk
      integer :: lchunk
      integer :: col
   end type clump2chunk
   type(clump2chunk), dimension(:,:), allocatable, private :: clump2chunks

   type chunk2clump
      integer :: clumpid
      integer :: cell
   end type chunk2clump
   type(chunk2clump), dimension(:,:), allocatable, private :: chunk2clumps

   real(r8), dimension(:), target, allocatable :: lp_sendbuf ! lnd->phys send buf
   real(r8), dimension(:), target, allocatable :: lp_recvbuf ! lnd->phys receive buf
   real(r8), dimension(:), target, allocatable :: pl_sendbuf ! phys->lnd send buf
   real(r8), dimension(:), target, allocatable :: pl_recvbuf ! phys->lnd receive buf

   integer, dimension(:), allocatable :: lp_blkcnts ! l->p send/p->l recv blocks
   integer, dimension(:), allocatable :: lp_sndcnts ! lnd->phys send counts
   integer, dimension(:), allocatable :: lp_rcvcnts ! lnd->phys receive counts
   integer, dimension(:), allocatable :: lp_sdispls ! lnd->phys send dsplsmnt
   integer, dimension(:), allocatable :: lp_rdispls ! lnd->phys receive dsplsmnt

   integer, dimension(:), allocatable :: pl_blkcnts ! p->l send/l->p recv blocks
   integer, dimension(:), allocatable :: pl_sndcnts ! phys->lnd send counts
   integer, dimension(:), allocatable :: pl_rcvcnts ! phys->lnd receive counts
   integer, dimension(:), allocatable :: pl_sdispls ! phys->lnd send dsplsmnt
   integer, dimension(:), allocatable :: pl_rdispls ! phys->lnd receive dsplsmnt

   integer, parameter :: pl_nval = 16               ! phys->lnd flux values
   integer, parameter :: lp_nval = 13               ! lnd->phys flux values

   logical :: lpc_init_flag = .false.               ! set if initialized
   logical :: noland = .false.                      ! Flag if no land points here
!---------------------------------------------------------------------------

contains

!---------------------------------------------------------------------------
!BOP
!
! !IROUTINE: clm_camInit
!
! !INTERFACE:
  subroutine clm_camInit(srfflx_land, srfflx_merge)
!
! !DESCRIPTION:
! Initialize land surface model and obtain relevant atmospheric model arrays
! back from (i.e. albedos, surface temperature and snow cover over land).
!
! !USES:
    use shr_orb_mod     , only : shr_orb_decl
    use radiation       , only : radiation_get
    use ppgrid          , only : pcols
    use camsrfexch_types, only : srfflx_parm, surface_state, srfflx_state, srfcomp2hub_alloc
    use time_manager    , only : get_nstep, get_step_size, get_curr_calday
    use filenames       , only : mss_irt, caseid
    use history         , only : ctitle, inithist, nhtfrq, mfilt
    use clm_varctl      , only : cam_caseid, cam_ctitle, cam_irad, cam_nsrest, &
                                 cam_crtinic, cam_nhtfrq, cam_mfilt, cam_irt, finidat       
    use initializeMod   , only : initialize           
    use initSurfAlbMod  , only : initSurfAlb, do_initsurfalb 
    use lnd2atmMod      , only : lnd2atm
    use clm_varsur      , only : landfrac, numlon
    use clm_varpar      , only : lsmlat
#ifdef SCAM
#include <max.h>
    use scamMod         , only : switch, have_tg, isrestart
#endif
#include <comsol.h>
#include <comctl.h>
!
! !ARGUMENTS:
    type(srfflx_parm) , pointer    :: srfflx_land(:)
    type(srfflx_state), intent(in) :: srfflx_merge(begchunk:endchunk)
!
! !LOCAL VARIABLES:
    integer  :: i,j         ! indices
    real(r8) :: dtime       ! time step increment (sec)
    integer  :: nstep       ! model time step
    real(r8) :: calday      ! calendar day for nstep
    real(r8) :: caldaym1    ! calendar day for nstep-1
    real(r8) :: declin      ! solar declination angle in radians for nstep
    real(r8) :: declinm1    ! solar declination angle in radians for nstep-1
    real(r8) :: eccf        ! earth orbit eccentricity factor
    integer  :: ier         ! returned error code
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
!-----------------------------------------------------------------------

    call srfcomp2hub_alloc( srfflx_land )
    if ( adiabatic .or. ideal_phys .or. aqua_planet )then
        noland = .true.
        return
    end if
#if ( defined SCAM )
    if(switch(CRM_SW+1))then
       noland = .true.
       return
    end if
    if ( isrestart ) return
    if(.not. have_tg) return
#endif

    ! Set namelist variables that must be same

    call radiation_get(iradsw_out=cam_irad)
    cam_caseid  = caseid
    cam_ctitle  = ctitle
    cam_nsrest  = nsrest
    cam_crtinic = inithist
    cam_nhtfrq  = nhtfrq(1)
    cam_mfilt   = mfilt(1)
    cam_irt     = mss_irt

    ! Initialize clm

    call initialize( )

    ! Determine consistency with cam grid info - only need to do this on 
    ! master processor (note that cam latitudes and longitudes are
    ! computed each time by the cam model at startup 

    call clm_camCheckGrid(srfflx_land, srfflx_merge)

    ! Initialize albedos (correct pft filters are needed)

    if (nsrest == 0) then
       if (finidat == ' ' .or. do_initsurfalb) then
          calday = get_curr_calday()
          call shr_orb_decl( calday, eccen, mvelpp, lambm0, obliqr, declin, eccf )
          
          dtime = get_step_size()
          caldaym1 = get_curr_calday(offset=-int(dtime))
          call shr_orb_decl( caldaym1, eccen, mvelpp, lambm0, obliqr, declinm1, eccf )
          
          call initSurfAlb( calday, declin, declinm1 )
       end if
    end if
       
    ! Initialize mapping between the atmosphere physics chunks and the land clumps

    call lp_coupling_init()

    ! For initial run only 

    if (get_nstep() == 0) then

       ! Determine gridcell averaged properties to send to atm

       call lnd2atm(init=.true.)

       ! Get initial data back from land model in atmospheric data structures

#ifdef TIMING_BARRIERS
       call t_startf('sync_cl2ch_ini')
       call mpi_barrier(mpicom)
       call t_stopf('sync_cl2ch_ini')
#endif

       call t_startf('clump2chunkini')
       call lp_all2all_clump_to_chunk_init(srfflx_land)
       call t_stopf('clump2chunkini')
    endif

    ! Determine if have land point 

    noland=.true.
    do j = 1,lsmlat
       do i = 1,numlon(j)
          if (landfrac(i,j) > 0._r8) noland = .false.
       end do
    end do

  end subroutine clm_camInit

!---------------------------------------------------------------------------
!BOP
!
! !IROUTINE: clm_camCheckGrid
!
! !INTERFACE:
  subroutine clm_camCheckGrid(srfflx_land, srfflx_merge)
!
! !DESCRIPTION:
! Check that cam grid is consistent with clm grid read in from clm surface 
! dataset
!
! !USES:
    use shr_const_mod   , only : SHR_CONST_PI
    use commap          , only : clat, londeg
    use rgrid           , only : nlon
    use pmgrid          , only : plon, plat
    use ppgrid          , only : pcols
    use phys_grid       , only : gather_chunk_to_field
    use camsrfexch_types, only : srfflx_state, srfflx_parm
    use clm_varsur      , only : landfrac, landmask, latixy, longxy, numlon
    use spmdMod         , only : masterproc
!BFB_BEGIN
    use clmtype
    use decompMod       , only : get_proc_global
!BFB_END
!
! !ARGUMENTS:
    type(srfflx_parm) , pointer    :: srfflx_land(:) 
    type(srfflx_state), intent(in) :: srfflx_merge(begchunk:endchunk)
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein 2005-05-14
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: i,j                     ! indices
    integer :: ext_numlon(plat)        ! ext number of longitudes
    real(r8):: ext_longxy(plon,plat)   ! ext lon values
    real(r8):: ext_latixy(plon,plat)   ! ext lat values
    real(r8):: ext_landfrac(plon,plat) ! ext fractional land
    integer :: ext_landmask(plon,plat) ! ext land mask
    real(r8):: cam_landfrac(pcols,begchunk:endchunk) ! local landfrac
!BFB_BEGIN
    integer :: numg         ! total number of gridcells across all processors
    integer :: numl         ! total number of landunits across all processors
    integer :: numc         ! total number of columns across all processors
    integer :: nump         ! total number of pfts across all processors
    integer :: g
!BFB_END
!---------------------------------------------------------------------------

    ! Determine cam grid info

#ifdef TIMING_BARRIERS
    call t_startf('sync_gather_landfrac')
    call mpi_barrier(mpicom)
    call t_stopf('sync_gather_landfrac')
#endif

    do i = begchunk, endchunk
       srfflx_land(i)%areafrac(:pcols) = srfflx_merge(i)%landfrac(:pcols)
    end do

    ! Determine consistency of cam and clm grid - all processors

    do j = 1,plat
       ext_numlon(j) = nlon(j)	
       if (ext_numlon(j) /= numlon(j)) then
          write(6,*)'clm_camInit error: CAM numlon array different from surface dataset'
          call endrun()
       end if
       do i = 1,nlon(j)
          ext_longxy(i,j) = londeg(i,j)
          ext_latixy(i,j) = (180._r8/SHR_CONST_PI)*clat(j)
          if ( abs(ext_latixy(i,j)-latixy(i,j)) > 1.e-12_r8 ) then
             write(6,*)'clm_camInit error: CAM latitude ',ext_latixy(i,j),' and clm input latitude ', &
                  latixy(i,j),' has difference too large at i,j= ',i,j
             call endrun()
          end if
          if ( abs(ext_longxy(i,j)-longxy(i,j)) > 1.e-12_r8 ) then
             write(6,*)'clm_camInit error: CAM longitude ',ext_longxy(i,j),' and clm input longitude ', &
                  longxy(i,j),' has difference too large at i,j= ',i,j
             call endrun()
          end if
       end do
    end do

    ! Determine consistency of cam and clm landfrac/landmask - masterproc only
    ! First, gather global landfrac on master processor only

    do i = begchunk, endchunk
       cam_landfrac(:pcols,i) = srfflx_merge(i)%landfrac(:pcols)
    end do
    call gather_chunk_to_field(1, 1, 1, plon, cam_landfrac, ext_landfrac)

    if (masterproc) then
       do j = 1,plat
          do i = 1,nlon(j)
             if (ext_landfrac(i,j) > 0._r8) then
                ext_landmask(i,j) = 1
             else
                ext_landmask(i,j) = 0
             endif
             if (ext_landmask(i,j) /= landmask(i,j)) then
                write(6,*)'clm_camInit error: CAM land mask different from surface dataset at i,j= ',i,j
                call endrun()
             end if
             if (ext_landfrac(i,j) /= landfrac(i,j)) then
                write(6,*)'clm_camInit error: CAM fractional land differs from surface dataset at i,j= ',i,j
                call endrun()
             end if
          end do
       end do
    end if
       
!BFB_BEGIN
!    call get_proc_global(numg, numl, numc, nump)
!    do g = 1,numg
!       i = clm3%g%ixy(g)
!       j = clm3%g%jxy(g)
!       clm3%g%lat(g)    = ext_latixy(i,j) * SHR_CONST_PI/180._r8  
!       clm3%g%lon(g)    = ext_longxy(i,j) * SHR_CONST_PI/180._r8
!       clm3%g%latdeg(g) = ext_latixy(i,j) 
!       clm3%g%londeg(g) = ext_longxy(i,j) 
!    end do
!BFB_END

  end subroutine clm_camCheckGrid

!---------------------------------------------------------------------------
!BOP
!
! !IROUTINE: clm_camRun
!
! !INTERFACE:
  subroutine clm_camRun(srf_state, srfflx_land)
!
! !DESCRIPTION:
! Pack data to be sent to land model into a single array.  Send data to
! land model and call land model driver.  Receive data back from land
! model in a single array.  Unpack this data into component arrays.
! NOTE: component arrays are contained in module comsrf.  When coupling to
! an atmospheric model: solar radiation depends on surface albedos from
! the previous time step (based on current surface conditions and solar
! zenith angle for next time step).  Longwave radiation depends on upward
! longwave flux from previous time step.
!
! !USES:
    use shr_orb_mod     , only : shr_orb_decl
    use camsrfexch_types, only : srfflx_parm, surface_state
    use clm_varctl      , only : irad 
    use time_manager    , only : get_nstep, get_step_size, get_curr_calday
    use lnd2atmMod      , only : lnd2atm
    use driver          , only : driver1, driver2
#include <comsol.h>
!
! !ARGUMENTS:
    type(srfflx_parm)  , intent(inout) :: srfflx_land(begchunk:endchunk)
    type(surface_state), intent(inout) :: srf_state(begchunk:endchunk)
!
! !LOCAL VARIABLES:
    logical  :: doalb                 ! true if surface albedo calculation time step
    real(r8) :: dtime                 ! time step increment (sec)
    integer  :: nstep                 ! model time step
    real(r8) :: caldayp1              ! calendar day for nstep+1
    real(r8) :: declinp1              ! solar declination angle in radians for nstep+1
    real(r8) :: eccf                  ! earth orbit eccentricity factor
    integer  :: ier                   ! returned error code
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
!---------------------------------------------------------------------------

    if (noland) return

#ifdef TIMING_BARRIERS
    call t_startf ('sync_lnd')
    call mpibarrier (mpicom)
    call t_stopf ('sync_lnd')
#endif

    ! Fill in cam chunked structures

#ifdef TIMING_BARRIERS
    call t_startf('sync_chnk2clmp')
    call mpi_barrier(mpicom)
    call t_stopf('sync_chnk2clmp')
#endif

    call t_startf('chunk2clump')
    call lp_all2all_chunk_to_clump(srf_state)
    call t_stopf('chunk2clump')
    
    ! Determine doalb (true when the next time step is a radiation time step) 

    nstep = get_nstep()
    doalb = ((irad==1) .or. (mod(nstep,irad)==0 .and. nstep/=0))

    ! Determin declination angle for next time step
    
    dtime = get_step_size()
    caldayp1 = get_curr_calday( offset=int(dtime) )
    call shr_orb_decl( caldayp1, eccen, mvelpp, lambm0, obliqr, declinp1, eccf )

    ! Call land model driver1
    
    call t_startf('clm_driver1')
    call driver1(doalb, caldayp1, declinp1)
    call t_stopf('clm_driver1')

    ! Call land model driver2
    
    call t_startf('clm_driver2')
    call driver2(caldayp1, declinp1)
    call t_stopf('clm_driver2')

    ! Determine gridcell averaged properties to send to atm (l2as and l2af derived types)

    call lnd2atm()

    ! Fill in cam chunked structures

#ifdef TIMING_BARRIERS
    call t_startf('sync_clmp2chnk')
    call mpi_barrier(mpicom)
    call t_stopf('sync_clmp2chnk')
#endif
    
    call t_startf('clump2chunk')
    call lp_all2all_clump_to_chunk(srfflx_land)
    call t_stopf('clump2chunk')

#ifdef TIMING_BARRIERS
    call t_startf ('sync_after_lnd')
    call mpibarrier (mpicom)
    call t_stopf ('sync_after_lnd')
#endif

  end subroutine clm_camRun

!---------------------------------------------------------------------------
!BOP
!
! !IROUTINE: clm_camFinal
!
! !INTERFACE:
  subroutine clm_camFinal(srfflx_land)
!
! !DESCRIPTION:
! Finalize land surface model
!
! !USES:
    use camsrfexch_types, only : srfflx_parm
!------------------------------------------------------------------------------
!BOP
!
! !ARGUMENTS:
    type(srfflx_parm)  , pointer :: srfflx_land(:)
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
!---------------------------------------------------------------------------

    deallocate( srfflx_land )
  end subroutine clm_camFinal

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: lp_coupling_init
!
! !INTERFACE:
   subroutine lp_coupling_init()
!
! !DESCRIPTION:
! This subroutine initializes the mapping between the atmosphere physics
! chunks and the land clumps.  It may (and must) be called repeatedly to
! re-initialize the mapping if the decomposition of either the atmosphere
! physics or the land changes.  It allocates communication buffers
! constructs vectors of counts and displacements used for subsequent
! communication between MPI processes.
!
! !USES
  use phys_grid    , only : get_chunk_coord_owner_p
!
! !ARGUMENTS:
   implicit none
!
! !LOCAL VARIABLES:
   integer :: p, c, g                            ! loop indices
   integer :: nclumps                            ! number of clumps defined
   integer :: ncells                             ! number of clump cells
   integer :: clump_owner                        ! clump owner
   integer, dimension(:), allocatable :: lons    ! clump longitudes
   integer, dimension(:), allocatable :: lats    ! clump latitudes
   integer, dimension(:), allocatable :: lchnks  ! chunk ids
   integer, dimension(:), allocatable :: cols    ! chunk columns
   integer, dimension(:), allocatable :: chunk_owners  ! chunk owners
   integer :: max_gpc = 0                        ! max cells per clump
   integer :: ier                                ! error codes
!
! !REVISION HISTORY:
! 2002.09.11  Forrest Hoffman  Creation.
! 2002.11.18  Mariana Vertenstein and Forrest Hoffman Updated to clm2.1.
!
!EOP
!------------------------------------------------------------------------------

   ! If already initialized, then deallocate buffers and re-initialize everything

   call lp_coupling_finalize()

   allocate(lp_blkcnts(0:npes-1), lp_sndcnts(0:npes-1), lp_rcvcnts(0:npes-1), &
            pl_blkcnts(0:npes-1), pl_sndcnts(0:npes-1), pl_rcvcnts(0:npes-1), &
            stat=ier)
   if (ier /= 0) then
      write (6,*) 'lp_coupling_init(): allocation error for lp_blkcnts, ', &
         'lp_sndcnts, pl_blkcnts, and pl_sndcnts'
      call endrun
   end if
   lp_blkcnts(:) = 0
   lp_sndcnts(:) = 0
   lp_rcvcnts(:) = 0
   pl_blkcnts(:) = 0
   pl_sndcnts(:) = 0
   pl_rcvcnts(:) = 0

   ! Determine max_gpc and allocate dynamic memory

   nclumps = get_nclumps()
   do c = 1,nclumps
      ncells = get_clump_ncells_id(c)
      if (ncells > max_gpc) max_gpc = ncells
   end do
   allocate(lons(max_gpc), lats(max_gpc), lchnks(max_gpc), cols(max_gpc), &
        chunk_owners(max_gpc), stat=ier)
   if (ier /= 0) then
      write (6,*) 'lp_coupling_init(): allocation error for local lons, ', &
         'lats, lchnks, cols, and chunk_owners variables'
      call endrun
   end if

   ! Found above already
   ! nclumps = get_nclumps()
   ! lp_blkcnts: for each cam pid, determine the total number of sends that
   ! will be sent from clm
   ! pl_blkcnts: for each clm pid, determine the total number of sends that
   ! will be sent from cam

   do c = 1,nclumps
      clump_owner = get_clump_owner_id(c)
      ncells = get_clump_ncells_id(c)
      call get_clump_coord_id(c, ncells, lons, lats)
      call get_chunk_coord_owner_p(ncells, lons, lats, lchnks, cols, chunk_owners)
      do g = 1,ncells
         if (clump_owner == iam) then
            p = chunk_owners(g)
            lp_blkcnts(p) = lp_blkcnts(p) + 1
         endif
         if (chunk_owners(g) == iam) then
            p = clump_owner
            pl_blkcnts(p) = pl_blkcnts(p) + 1
         endif
      end do
   end do

   allocate(clump2chunks(0:npes-1, 1:maxval(pl_blkcnts)), &
            chunk2clumps(0:npes-1, 1:maxval(lp_blkcnts)), stat=ier)
   if (ier /= 0) then
      write (6,*) 'lp_coupling_init(): allocation error for clump2chunks ', &
         'and chunk2clumps'
      call endrun
   end if
   clump2chunks(:,:)%lchunk = 0
   clump2chunks(:,:)%col = 0
   chunk2clumps(:,:)%clumpid = 0
   chunk2clumps(:,:)%cell = 0


   ! Found above already
   ! nclumps = get_nclumps()

   lp_blkcnts(:) = 0
   pl_blkcnts(:) = 0
   do c = 1,nclumps
      clump_owner = get_clump_owner_id(c)
      ncells = get_clump_ncells_id(c)
      call get_clump_coord_id(c, ncells, lons, lats)
      call get_chunk_coord_owner_p(ncells,lons,lats,lchnks,cols,chunk_owners)
      do g = 1,ncells
         if (clump_owner == iam) then
            p = chunk_owners(g)
            lp_blkcnts(p) = lp_blkcnts(p) + 1
            chunk2clumps(p,lp_blkcnts(p))%clumpid=c
            chunk2clumps(p,lp_blkcnts(p))%cell = g
         end if
         if (chunk_owners(g) == iam) then
            p = clump_owner
            pl_blkcnts(p) = pl_blkcnts(p) + 1
            clump2chunks(p,pl_blkcnts(p))%lchunk = lchnks(g)
            clump2chunks(p,pl_blkcnts(p))%col = cols(g)
         end if
      end do
   end do

   deallocate(lons, lats, lchnks, cols, chunk_owners)

   pl_sndcnts(:) = pl_blkcnts(:) * pl_nval
   lp_rcvcnts(:) = pl_blkcnts(:) * lp_nval
   lp_sndcnts(:) = lp_blkcnts(:) * lp_nval
   pl_rcvcnts(:) = lp_blkcnts(:) * pl_nval

   allocate(pl_sendbuf(0:sum(pl_sndcnts)-1), lp_recvbuf(0:sum(lp_rcvcnts)-1), &
      stat=ier)
   if (ier /= 0) then
      write (6,*) 'lp_coupling_init(): allocation error for pl_sendbuf and ', &
         'lp_recvbuf'
      call endrun
   end if
   allocate(lp_sendbuf(0:sum(lp_sndcnts)-1), pl_recvbuf(0:sum(pl_rcvcnts)-1), &
      stat=ier)
   if (ier /= 0) then
      write (6,*) 'lp_coupling_init(): allocation error for lp_sendbuf and ', &
         'pl_recvbuf'
      call endrun
   end if

   allocate(lp_sdispls(0:npes-1), lp_rdispls(0:npes-1), pl_sdispls(0:npes-1), &
      pl_rdispls(0:npes-1), stat=ier)
   if (ier /= 0) then
      write (6,*) 'lp_coupling_init(): allocation error for lp_sdispls, ', &
         'lp_rdispls, pl_sdispls, and pl_rdispls'
      call endrun
   end if

   lp_sdispls(0) = 0
   lp_rdispls(0) = 0
   pl_sdispls(0) = 0
   pl_rdispls(0) = 0
   do p = 1,npes-1
      lp_sdispls(p) = lp_sdispls(p-1) + lp_sndcnts(p-1)
      lp_rdispls(p) = lp_rdispls(p-1) + lp_rcvcnts(p-1)
      pl_sdispls(p) = pl_sdispls(p-1) + pl_sndcnts(p-1)
      pl_rdispls(p) = pl_rdispls(p-1) + pl_rcvcnts(p-1)
   end do

   lpc_init_flag = .true.

   end subroutine lp_coupling_init

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: lp_coupling_finalize
!
! !INTERFACE:
   subroutine lp_coupling_finalize()
!
! !ARGUMENTS:
   implicit none
!
! !DESCRIPTION:
! This subroutine destroys the mapping between the atmsphere physics
! chunks and the land clumps if the lpc\_init\_flag flag is set.  It is
! called from lp\_coupling\_init() to ensure memory is recycled when a
! new mapping is to be created.
!
! !REVISION HISTORY:
! 2002.09.11  Forrest Hoffman  Creation.
! 2002.11.18  Mariana Vertenstein and Forrest Hoffman Updated to clm2.1.
!
!EOP
!------------------------------------------------------------------------------

   if (lpc_init_flag) then
      deallocate(clump2chunks, chunk2clumps)
      deallocate(lp_sendbuf, lp_recvbuf, pl_sendbuf, pl_recvbuf)
      deallocate(lp_blkcnts, pl_blkcnts)
      deallocate(lp_sndcnts, lp_rcvcnts, pl_sndcnts, pl_rcvcnts)
      deallocate(lp_sdispls, lp_rdispls, pl_sdispls, pl_rdispls)
      lpc_init_flag = .false.
   endif

   end subroutine lp_coupling_finalize

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: lp_all2all_clump_to_chunk_init
!
! !INTERFACE:
   subroutine lp_all2all_clump_to_chunk_init (srfflx_land)
!
! !DESCRIPTION:
! This subroutine performs the initial communication from the land model
! to the atmosphere physics (from clumps to chunks) based on the mapping
! constructed in lp\_coupling\_init().
!
! !USES:
   use ppgrid          , only : begchunk, endchunk
   use comsrf          , only : snowhland
   use camsrfexch_types, only : srfflx_parm
   use clm_varcon      , only : sb
   use clmtype
!
! !ARGUMENTS:
   implicit none
   type(srfflx_parm), intent(inout), dimension(begchunk:endchunk) :: srfflx_land
!
! !REVISION HISTORY:
! 2002.09.11  Forrest Hoffman  Creation.
! 2002.11.18  Mariana Vertenstein and Forrest Hoffman Updated to clm2.1.
! 2003.05.01  Mariana Vertenstein Updated to l2as data structures
!
!EOP
!
! !LOCAL VARIABLES:
   integer  :: p, n, k, m, g              ! indices
   integer  :: is                         ! send buffer index
   integer  :: lchnk, i                   ! local chunk and column
   integer  :: ier                        ! returned error code
   real(r8), pointer :: lp_rbufp(:)       ! recv buffer pointer
   type(gridcell_type), pointer :: gptr   ! pointer to gridcell derived subtype
!------------------------------------------------------------------------------

   ! Set pointers into derived type

   gptr => clm3%g

   ! Fill lnd->phys send buffer

   lp_sendbuf(:) = 0.0_r8

   do p = 0, npes-1
      do n = 1, lp_blkcnts(p)

         ! Determine clump 1d gridcell index
         call get_clump_gcell_info (chunk2clumps(p,n)%clumpid, chunk2clumps(p,n)%cell, g)

         ! Fill in send buffer
         is = lp_sdispls(p)+(n-1)*lp_nval + 0
         lp_sendbuf(is) = sqrt(sqrt(gptr%l2af%eflx_lwrad_out(g)/sb)) ! tsxy
         is = lp_sdispls(p)+(n-1)*lp_nval + 1
         lp_sendbuf(is) = gptr%l2as%albd(g,1)                ! asdir
         is = lp_sdispls(p)+(n-1)*lp_nval + 2
         lp_sendbuf(is) = gptr%l2as%albd(g,2)                ! aldir
         is = lp_sdispls(p)+(n-1)*lp_nval + 3
         lp_sendbuf(is) = gptr%l2as%albi(g,1)                ! asdif
         is = lp_sdispls(p)+(n-1)*lp_nval + 4
         lp_sendbuf(is) = gptr%l2as%albi(g,2)                ! aldif
         is = lp_sdispls(p)+(n-1)*lp_nval + 5
         lp_sendbuf(is) = gptr%l2as%h2osno(g)                ! snow [mm]->[m]
         is = lp_sdispls(p)+(n-1)*lp_nval + 6
         lp_sendbuf(is) = 1.e36_r8                              ! spval
         is = lp_sdispls(p)+(n-1)*lp_nval + 7
         lp_sendbuf(is) = 1.e36_r8                              ! spval
         is = lp_sdispls(p)+(n-1)*lp_nval + 8
         lp_sendbuf(is) = 1.e36_r8                              ! spval
         is = lp_sdispls(p)+(n-1)*lp_nval + 9
         lp_sendbuf(is) = 1.e36_r8                              ! spval
         is = lp_sdispls(p)+(n-1)*lp_nval + 10
         lp_sendbuf(is) = gptr%l2af%eflx_lwrad_out(g)        ! lwrad
         is = lp_sdispls(p)+(n-1)*lp_nval + 11
         lp_sendbuf(is) = 1.e36_r8                              ! spval
         is = lp_sdispls(p)+(n-1)*lp_nval + 12
         lp_sendbuf(is) = 1.e36_r8                              ! spval

      end do
   end do

#ifdef SPMD
   call mpi_alltoallv (lp_sendbuf, lp_sndcnts, lp_sdispls, mpir8, &
                       lp_recvbuf, lp_rcvcnts, lp_rdispls, mpir8, mpicom, ier)
   if (ier /= 0) then
      write (6,*) 'lp_all2all_clump_to_chunk_init(): MPI error ', ier, &
         ' from mpi_alltoallv()'
      call endrun
   endif
   lp_rbufp => lp_recvbuf
#else
   lp_rbufp => lp_sendbuf
#endif

   ! Extract lnd->phys receive buffer

   do p = 0,npes-1
      do n = 1,pl_blkcnts(p)
         lchnk = clump2chunks(p,n)%lchunk
         i     = clump2chunks(p,n)%col
         srfflx_land(lchnk)%ts(i)    = lp_rbufp(lp_rdispls(p)+(n-1)*lp_nval+0)
         srfflx_land(lchnk)%asdir(i) = lp_rbufp(lp_rdispls(p)+(n-1)*lp_nval+1)
         srfflx_land(lchnk)%aldir(i) = lp_rbufp(lp_rdispls(p)+(n-1)*lp_nval+2)
         srfflx_land(lchnk)%asdif(i) = lp_rbufp(lp_rdispls(p)+(n-1)*lp_nval+3)
         srfflx_land(lchnk)%aldif(i) = lp_rbufp(lp_rdispls(p)+(n-1)*lp_nval+4)
         snowhland(i,lchnk)       = lp_rbufp(lp_rdispls(p)+(n-1)*lp_nval+5)
         srfflx_land(lchnk)%lwup(i)  = lp_rbufp(lp_rdispls(p)+(n-1)*lp_nval+10)
      end do
   end do

   end subroutine lp_all2all_clump_to_chunk_init

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: lp_all2all_clump_to_chunk
!
! !INTERFACE:
   subroutine lp_all2all_clump_to_chunk (srfflx_land)
!
! !DESCRIPTION:
! This subroutine performs the communication from the land model to
! the atmosphere physics (from clumps to chunks) based on the mapping
! constructed in lp\_coupling\_init().
!
! !USES:
   use ppgrid          , only : begchunk, endchunk
   use comsrf          , only : snowhland
   use camsrfexch_types, only : srfflx_parm
   use constituents    , only : pcnst, pnats
   use clmtype
!
! !ARGUMENTS:
   implicit none
   type(srfflx_parm), intent(inout), dimension(begchunk:endchunk) :: srfflx_land
!
! !REVISION HISTORY:
! 2002.09.11  Forrest Hoffman  Creation.
! 2002.11.18  Mariana Vertenstein and Forrest Hoffman Updated to clm2.1.
! 2003.05.01  Mariana Vertenstein Updated to l2as data structures
!
!EOP
!
! !LOCAL VARIABLES:
   integer  :: p, n, k, m, g              ! loop indices
   integer  :: bpatch, npatch             ! patch id and count
   integer  :: lchnk, i                   ! local chunk and column
   integer  :: ier                        ! returned error code
   integer  :: is                         ! send buffer index
   real(r8) :: wt                         ! patch wt
   real(r8), pointer :: lp_rbufp(:)       ! recv buffer pointer
   type(gridcell_type), pointer :: gptr   ! pointer to gridcell derived subtype
!------------------------------------------------------------------------------
! lp_coupling.F90,v 1.1.6.12 2004/05/07 21:19:22 forrest Exp
! forrest
!------------------------------------------------------------------------------

    ! Set pointers into derived type

    gptr => clm3%g

   ! Fill lnd->phys send buffer

   lp_sendbuf(:) = 0.0_r8

   do p = 0,npes-1
      do n = 1,lp_blkcnts(p)

         ! Determine clump gridcell info
         call get_clump_gcell_info(chunk2clumps(p,n)%clumpid, chunk2clumps(p,n)%cell, g)

         ! Fill in send buffer
         is = lp_sdispls(p)+(n-1)*lp_nval+0   ! tsxy
         lp_sendbuf(is) = gptr%l2as%t_rad(g)
         is = lp_sdispls(p)+(n-1)*lp_nval+5   ! snow [mm]->[m]
         lp_sendbuf(is) = gptr%l2as%h2osno(g)
         is = lp_sdispls(p)+(n-1)*lp_nval+1   ! asdir
         lp_sendbuf(is) = gptr%l2as%albd(g,1)
         is = lp_sdispls(p)+(n-1)*lp_nval+2   ! aldir
         lp_sendbuf(is) = gptr%l2as%albd(g,2)
         is = lp_sdispls(p)+(n-1)*lp_nval+3   ! asdif
         lp_sendbuf(is) = gptr%l2as%albi(g,1)
         is = lp_sdispls(p)+(n-1)*lp_nval+4   ! aldif
         lp_sendbuf(is) = gptr%l2as%albi(g,2)
         is = lp_sdispls(p)+(n-1)*lp_nval+6   ! taux
         lp_sendbuf(is) = gptr%l2af%taux(g)
         is = lp_sdispls(p)+(n-1)*lp_nval+7   ! tauy
         lp_sendbuf(is) = gptr%l2af%tauy(g)
         is = lp_sdispls(p)+(n-1)*lp_nval+8   ! lhflx
         lp_sendbuf(is) = gptr%l2af%eflx_lh_tot(g)
         is = lp_sdispls(p)+(n-1)*lp_nval+9   ! shflx
         lp_sendbuf(is) = gptr%l2af%eflx_sh_tot(g)
         is = lp_sdispls(p)+(n-1)*lp_nval+10  ! lwrad
         lp_sendbuf(is) = gptr%l2af%eflx_lwrad_out(g)
         is = lp_sdispls(p)+(n-1)*lp_nval+11  ! qflx
         lp_sendbuf(is) = gptr%l2af%qflx_evap_tot(g)
         is = lp_sdispls(p)+(n-1)*lp_nval+12  ! tref
         lp_sendbuf(is) = gptr%l2as%t_ref2m(g)

      end do   ! end loop over destination processes
   end do   ! end loop over source processes

#ifdef SPMD
   call mpi_alltoallv (lp_sendbuf, lp_sndcnts, lp_sdispls, mpir8, &
                       lp_recvbuf, lp_rcvcnts, lp_rdispls, mpir8, mpicom, ier)
   if (ier /= 0) then
      write (6,*) 'alltoall_clump_to_chunk(): MPI error ', ier, &
         ' from mpi_alltoallv()'
      call endrun
   endif
   lp_rbufp => lp_recvbuf
#else
   lp_rbufp => lp_sendbuf
#endif

   ! Extract lnd->phys receive buffer

   do p = 0,npes-1
      do n = 1,pl_blkcnts(p)
         lchnk = clump2chunks(p,n)%lchunk
         i     = clump2chunks(p,n)%col
         srfflx_land(lchnk)%ts(i)     = lp_rbufp(lp_rdispls(p)+(n-1)*lp_nval+0)
         srfflx_land(lchnk)%asdir(i)  = lp_rbufp(lp_rdispls(p)+(n-1)*lp_nval+1)
         srfflx_land(lchnk)%aldir(i)  = lp_rbufp(lp_rdispls(p)+(n-1)*lp_nval+2)
         srfflx_land(lchnk)%asdif(i)  = lp_rbufp(lp_rdispls(p)+(n-1)*lp_nval+3)
         srfflx_land(lchnk)%aldif(i)  = lp_rbufp(lp_rdispls(p)+(n-1)*lp_nval+4)
         snowhland(i,lchnk)           = lp_rbufp(lp_rdispls(p)+(n-1)*lp_nval+5)
         srfflx_land(lchnk)%wsx(i)    = lp_rbufp(lp_rdispls(p)+(n-1)*lp_nval+6)
         srfflx_land(lchnk)%wsy(i)    = lp_rbufp(lp_rdispls(p)+(n-1)*lp_nval+7)
         srfflx_land(lchnk)%lhf(i)    = lp_rbufp(lp_rdispls(p)+(n-1)*lp_nval+8)
         srfflx_land(lchnk)%shf(i)    = lp_rbufp(lp_rdispls(p)+(n-1)*lp_nval+9)
         srfflx_land(lchnk)%lwup(i)   = lp_rbufp(lp_rdispls(p)+(n-1)*lp_nval+10)
         srfflx_land(lchnk)%cflx(i,1) = lp_rbufp(lp_rdispls(p)+(n-1)*lp_nval+11)
         srfflx_land(lchnk)%tref(i)   = lp_rbufp(lp_rdispls(p)+(n-1)*lp_nval+12)

         ! Reset all other constituent surface fluxes to zero over land
         do m = 2,pcnst+pnats
            srfflx_land(lchnk)%cflx(i,m) = 0.0_r8
         end do
      end do
   end do

   end subroutine lp_all2all_clump_to_chunk

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: lp_all2all_chunk_to_clump
!
! !INTERFACE:
   subroutine lp_all2all_chunk_to_clump (srf_state)
!
! !DESCRIPTION:
! This subroutine performs communication from the atmosphere physics to
! the land model (from chunks to clumps) based on the mapping constructed
! in lp\_coupling\_init().
!
! !USES:
   use ppgrid          , only: begchunk, endchunk
   use camsrfexch_types, only: surface_state
   use clmtype
   use clm_varcon      , only: rair, o2_molar_const, co2_ppmv_const, c13ratio
!
! !ARGUMENTS:
   implicit none
   type(surface_state), intent(in), dimension(begchunk:endchunk) :: srf_state
!
! !REVISION HISTORY:
! 2002.09.11  Forrest Hoffman  Creation.
! 2002.11.18  Mariana Vertenstein and Forrest Hoffman Updated to clm2.1.
!
!
!EOP
!
! !LOCAL VARIABLES:
   integer  :: g, p, n, k               ! indices
   integer  :: lchnk, i                 ! local chunk and column
   integer  :: ier                      ! returned error code
   real(r8) :: forc_rainc, forc_rainl   ! rainxy [mm/s]
   real(r8) :: forc_snowc, forc_snowl   ! snowfxy [mm/s]
   real(r8), pointer :: pl_rbufp(:)     ! recv buffer pointer
   type(gridcell_type), pointer :: gptr ! pointer to gridcell derived subtype
!------------------------------------------------------------------------------
! lp_coupling.F90,v 1.1.6.12 2004/05/07 21:19:22 forrest Exp
! forrest
!------------------------------------------------------------------------------

    ! Set pointers into derived type

    gptr => clm3%g

   ! Fill phys->lnd send buffer

   do p = 0,npes-1
      do n = 1,pl_blkcnts(p)
         lchnk = clump2chunks(p,n)%lchunk
         i = clump2chunks(p,n)%col

         ! Atmoshperic state variable [m]
         pl_sendbuf(pl_sdispls(p)+(n-1)*pl_nval+0) = srf_state(lchnk)%zbot(i)

         ! Atmoshperic state variable [m/s]
         pl_sendbuf(pl_sdispls(p)+(n-1)*pl_nval+1) = srf_state(lchnk)%ubot(i)

         ! Atmoshperic state variable [m/s]
         pl_sendbuf(pl_sdispls(p)+(n-1)*pl_nval+2) = srf_state(lchnk)%vbot(i)

         ! Atmoshperic state variable [K]
         pl_sendbuf(pl_sdispls(p)+(n-1)*pl_nval+3) = srf_state(lchnk)%thbot(i)

         ! Atmoshperic state variable [kg/kg]
         pl_sendbuf(pl_sdispls(p)+(n-1)*pl_nval+4) = srf_state(lchnk)%qbot(i,1)

         ! Atmoshperic state variable [Pa]
         pl_sendbuf(pl_sdispls(p)+(n-1)*pl_nval+5) = srf_state(lchnk)%pbot(i)

         ! Atmoshperic state variable [K]
         pl_sendbuf(pl_sdispls(p)+(n-1)*pl_nval+6) = srf_state(lchnk)%tbot(i)

         ! Atmoshperic flux [W/m^2]
         pl_sendbuf(pl_sdispls(p)+(n-1)*pl_nval+7) = srf_state(lchnk)%flwds(i)

         ! Atmoshperic flux [m/s] -> [mm/s]
         pl_sendbuf(pl_sdispls(p)+(n-1)*pl_nval+8) = &
              srf_state(lchnk)%precsc(i) * 1000._r8

         ! Atmoshperic flux [m/s] -> [mm/s]
         pl_sendbuf(pl_sdispls(p)+(n-1)*pl_nval+9) = &
            srf_state(lchnk)%precsl(i) * 1000._r8

         ! Atmoshperic flux [m/s] -> [mm/s]
         pl_sendbuf(pl_sdispls(p)+(n-1)*pl_nval+10) = &
            (srf_state(lchnk)%precc(i) - srf_state(lchnk)%precsc(i)) * 1000._r8

         ! Atmoshperic flux [m/s] -> [mm/s]
         pl_sendbuf(pl_sdispls(p)+(n-1)*pl_nval+11) = &
            (srf_state(lchnk)%precl(i) - srf_state(lchnk)%precsl(i)) * 1000._r8

         ! Atmoshperic flux [W/m^2]
         pl_sendbuf(pl_sdispls(p)+(n-1)*pl_nval+12) = srf_state(lchnk)%soll(i)

         ! Atmoshperic flux [W/m^2]
         pl_sendbuf(pl_sdispls(p)+(n-1)*pl_nval+13) = srf_state(lchnk)%sols(i)

         ! Atmoshperic flux [W/m^2]
         pl_sendbuf(pl_sdispls(p)+(n-1)*pl_nval+14) = srf_state(lchnk)%solld(i)

         ! Atmoshperic flux [W/m^2]
         pl_sendbuf(pl_sdispls(p)+(n-1)*pl_nval+15) = srf_state(lchnk)%solsd(i)

      end do
   end do

#ifdef SPMD
   call mpi_alltoallv (pl_sendbuf, pl_sndcnts, pl_sdispls, mpir8, &
                       pl_recvbuf, pl_rcvcnts, pl_rdispls, mpir8, mpicom, ier)
   if (ier /= 0) then
      write (6,*) 'lp_all2all_chunk_to_clump(): MPI error ', ier, &
         ' from mpi_alltoallv()'
      call endrun
   endif
   pl_rbufp => pl_recvbuf
#else
   pl_rbufp => pl_sendbuf
#endif

   ! Extract phys->lnd receive buffer

   do p = 0,npes-1
      do n = 1,lp_blkcnts(p)

         ! Determine clump gridcell info
         call get_clump_gcell_info(chunk2clumps(p,n)%clumpid, chunk2clumps(p,n)%cell, g)

         ! Fill in clmtype gridcell info
         gptr%a2ls%forc_hgt(g)     = pl_rbufp(pl_rdispls(p)+(n-1)*pl_nval+0)
         gptr%a2ls%forc_u(g)       = pl_rbufp(pl_rdispls(p)+(n-1)*pl_nval+1)
         gptr%a2ls%forc_v(g)       = pl_rbufp(pl_rdispls(p)+(n-1)*pl_nval+2)
         gptr%a2ls%forc_th(g)      = pl_rbufp(pl_rdispls(p)+(n-1)*pl_nval+3)
         gptr%a2ls%forc_q(g)       = pl_rbufp(pl_rdispls(p)+(n-1)*pl_nval+4)
         gptr%a2ls%forc_pbot(g)    = pl_rbufp(pl_rdispls(p)+(n-1)*pl_nval+5)
         gptr%a2ls%forc_t(g)       = pl_rbufp(pl_rdispls(p)+(n-1)*pl_nval+6)
         gptr%a2lf%forc_lwrad(g)   = pl_rbufp(pl_rdispls(p)+(n-1)*pl_nval+7)
         forc_snowc                = pl_rbufp(pl_rdispls(p)+(n-1)*pl_nval+8)
         forc_snowl                = pl_rbufp(pl_rdispls(p)+(n-1)*pl_nval+9)
         forc_rainc                = pl_rbufp(pl_rdispls(p)+(n-1)*pl_nval+10)
         forc_rainl                = pl_rbufp(pl_rdispls(p)+(n-1)*pl_nval+11)
         gptr%a2lf%forc_solad(g,2) = pl_rbufp(pl_rdispls(p)+(n-1)*pl_nval+12)
         gptr%a2lf%forc_solad(g,1) = pl_rbufp(pl_rdispls(p)+(n-1)*pl_nval+13)
         gptr%a2lf%forc_solai(g,2) = pl_rbufp(pl_rdispls(p)+(n-1)*pl_nval+14)
         gptr%a2lf%forc_solai(g,1) = pl_rbufp(pl_rdispls(p)+(n-1)*pl_nval+15)

         ! Determine derived quantities
         gptr%a2ls%forc_hgt_u(g) = gptr%a2ls%forc_hgt(g)  ! obs height of wind [m]
         gptr%a2ls%forc_hgt_t(g) = gptr%a2ls%forc_hgt(g)  ! obs height of temperature [m]
         gptr%a2ls%forc_hgt_q(g) = gptr%a2ls%forc_hgt(g)  ! obs height of humidity [m]
         gptr%a2ls%forc_vp(g)    = gptr%a2ls%forc_q(g) * gptr%a2ls%forc_pbot(g) &
                                   / (0.622_r8 + 0.378_r8 * gptr%a2ls%forc_q(g))
         gptr%a2ls%forc_rho(g)   = (gptr%a2ls%forc_pbot(g) - 0.378_r8 * gptr%a2ls%forc_vp(g)) &
                                   / (rair * gptr%a2ls%forc_t(g))
         gptr%a2ls%forc_pco2(g)  = co2_ppmv_const * 1.e-6_r8 * gptr%a2ls%forc_pbot(g)
         gptr%a2ls%forc_pc13o2(g) = co2_ppmv_const * c13ratio * 1.e-6_r8 * gptr%a2ls%forc_pbot(g)

         gptr%a2ls%forc_po2(g)   = o2_molar_const * gptr%a2ls%forc_pbot(g)
         gptr%a2ls%forc_wind(g)  = sqrt(gptr%a2ls%forc_u(g)**2 + gptr%a2ls%forc_v(g)**2)
         gptr%a2lf%forc_solar(g) = gptr%a2lf%forc_solad(g,1) + gptr%a2lf%forc_solai(g,1) + &
                                   gptr%a2lf%forc_solad(g,2) + gptr%a2lf%forc_solai(g,2)

         ! Determine precipitation needed by clm
#ifdef PERGRO
         ! For error growth only, allow rain, not snowfall
         forc_rainc           = forc_rainc + forc_snowc
         forc_rainl           = forc_rainl + forc_snowl
         forc_snowc           = 0.0_r8
         forc_snowl           = 0.0_r8
#endif
         gptr%a2lf%forc_rain(g) = forc_rainc + forc_rainl
         gptr%a2lf%forc_snow(g) = forc_snowc + forc_snowl

      end do   ! end loop over destination processes
   end do   ! end loop over source processes

   end subroutine lp_all2all_chunk_to_clump

#endif

end module clm_camMod
