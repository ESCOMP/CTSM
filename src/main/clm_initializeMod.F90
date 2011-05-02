module clm_initializeMod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: clm_initializeMod
!
! !DESCRIPTION:
! Performs land model initialization
!
! !USES:
  use shr_kind_mod    , only : r8 => shr_kind_r8
  use spmdMod         , only : masterproc
  use shr_sys_mod     , only : shr_sys_flush
  use abortutils      , only : endrun
  use clm_varctl      , only : nsrest, nsrStartup, nsrContinue, nsrBranch, downscale
  use clm_varctl      , only : iulog
  use clm_varctl      , only : create_glacier_mec_landunit
  use clm_varsur      , only : wtxy,vegxy
  use clm_varsur      , only : topoxy
  use clmtype         , only : gratm, grlnd, nameg, namel, namec, namep, allrof
  use perf_mod        , only : t_startf, t_stopf
  use ncdio_pio
  use mct_mod

! !PUBLIC TYPES:
  implicit none
  save

  private    ! By default everything is private

! !PUBLIC MEMBER FUNCTIONS:
  public :: initialize1  ! Phase one initialization
  public :: initialize2  ! Phase two initialization
!
! !REVISION HISTORY:
! Created by Gordon Bonan, Sam Levis and Mariana Vertenstein
!
!
! !PRIVATE MEMBER FUNCTIONS:
  private header         ! echo version numbers
  private do_restread    ! read a restart file
  private cellarea       ! area of grid cells (square km)
!-----------------------------------------------------------------------
! !PRIVATE DATA MEMBERS: None

!EOP
!-----------------------------------------------------------------------
contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: initialize1
!
! !INTERFACE:
  subroutine initialize1( )
!
! !DESCRIPTION:
! Land model initialization.
! o Initializes run control variables via the [clm_inparm] namelist.
! o Reads surface data on model grid.
! o Defines the multiple plant types and fraction areas for each surface type.
! o Builds the appropriate subgrid <-> grid mapping indices and weights.
! o Set up parallel processing.
! o Initializes time constant variables.
! o Reads restart data for a restart or branch run.
! o Reads initial data and initializes the time variant variables for an initial run.
! o Initializes history file output.
! o Initializes river routing model.
! o Initializes accumulation variables.
!
! !USES:
    use clm_varpar   , only : lsmlon, lsmlat, maxpatch, clm_varpar_init
    use pftvarcon    , only : pftconrd
    use decompInitMod, only : decompInit_atm, decompInit_lnd, decompInit_glcp
    use decompMod    , only : adecomp,ldecomp
    use decompMod    , only : get_proc_clumps, get_clump_bounds, &
                              get_proc_bounds, get_proc_bounds_atm
    use domainMod    , only : domain_check
    use domainMod    , only : adomain,ldomain
    use domainMod    , only : alatlon,llatlon,gatm,amask,pftm
    use domainMod    , only : latlon_check, latlon_setsame
    use surfrdMod    , only : surfrd,surfrd_get_grid,surfrd_get_frac,&
                              surfrd_get_topo, surfrd_get_latlon
    use clm_varctl   , only : fsurdat, fatmgrid, fatmlndfrc, &
                              fatmtopo, flndtopo, noland, downscale, fglcmask
    use controlMod   , only : control_init, control_print, nlfilename
    use UrbanInputMod, only : UrbanInput
    use ncdio_pio    , only : ncd_pio_init
    use downscaleMod , only : map_setgatmfm
!
! !ARGUMENTS:
!
! !REVISION HISTORY:
! Created by Gordon Bonan, Sam Levis and Mariana Vertenstein
!
!
! !LOCAL VARIABLES:
!EOP
    integer  :: ier                   ! error status
    integer  :: i,j,n1,n2,n,k,nr      ! loop indices
    integer  :: nl,nlg                ! gdc and glo lnd indices
    real(r8) :: rmaxlon,rmaxlat       ! local min/max vars
    integer  :: begg, endg            ! clump beg and ending gridcell indices
    integer  :: begg_atm, endg_atm    ! proc beg and ending gridcell indices
    character(len=32) :: subname = 'initialize1' ! subroutine name
!-----------------------------------------------------------------------

    ! ------------------------------------------------------------------------
    ! Initialize run control variables, timestep
    ! ------------------------------------------------------------------------

    call t_startf('init_control')
    call header()

    if (masterproc) then
       write(iulog,*) 'Attempting to initialize the land model .....'
       write(iulog,*)
       call shr_sys_flush(iulog)
    endif

    call clm_varpar_init()
    call control_init( )
    call ncd_pio_init( )

    if (masterproc) call control_print()

    call t_stopf('init_control')
    call t_startf('init_grids')

    ! ------------------------------------------------------------------------
    ! Read atm and land grids - global reads of data
    ! ------------------------------------------------------------------------

    if (masterproc) then
       write(iulog,*) 'Attempting to read alatlon from fatmgrid'
       call shr_sys_flush(iulog)
    endif

    call surfrd_get_latlon(alatlon, fatmgrid, amask, fatmlndfrc)
    call latlon_check(alatlon)
    if (masterproc) then
       write(iulog,*) 'amask size/min/max ',size(amask),minval(amask),maxval(amask)
       call shr_sys_flush(iulog)
    endif

    if (masterproc) then
       write(iulog,*) 'Attempting to read llatlon from fsurdat'
       call shr_sys_flush(iulog)
    endif

    call surfrd_get_latlon(llatlon, fsurdat, pftm, pftmflag=.true.)
    call latlon_check(llatlon)

    lsmlon = llatlon%ni
    lsmlat = llatlon%nj

    if (llatlon%ni < alatlon%ni .or. llatlon%nj < alatlon%nj) then
       if (masterproc) write(iulog,*) 'ERROR llatlon size > alatlon size: ', &
          n,llatlon%ni, llatlon%nj, alatlon%ni, alatlon%nj
       call endrun( trim(subname)//'ERROR: bad sizes' )
    endif

    !--- check if grids are "close", adjust, continue, or end  ---
    !--- set llatlon== alatlon if lats/lons < 0.001 different ---
    !--- exit if lats/lon > 0.001 and < 1.0 different          ---
    !--- continue if lat/lons > 1.0 different                  ---

    downscale = .true.
    if (alatlon%ni == llatlon%ni .and. alatlon%nj == llatlon%nj) then
       rmaxlon = 0.0_r8
       rmaxlat = 0.0_r8
       do n = 1,alatlon%ni
          rmaxlon = max(rmaxlon,abs(alatlon%lonc(n)-llatlon%lonc(n)))
       enddo
       do n = 1,alatlon%nj
          rmaxlat = max(rmaxlat,abs(alatlon%latc(n)-llatlon%latc(n)))
       enddo
       if (rmaxlon < 0.001_r8 .and. rmaxlat < 0.001_r8) then
          if (masterproc) write(iulog,*) trim(subname)//': set llatlon =~ alatlon', &
             ':continue',rmaxlon,rmaxlat
          call latlon_setsame(alatlon,llatlon)
          downscale = .false.
       elseif (rmaxlon < 1.0_r8 .and. rmaxlat < 1.0_r8) then
          if (masterproc) write(iulog,*) trim(subname)//': alatlon/llatlon mismatch', &
             ':error',rmaxlon,rmaxlat
          call endrun( trim(subname) )
       else
          if (masterproc) write(iulog,*) trim(subname)//': alatlon/llatlon different', &
              ':continue',rmaxlon,rmaxlat
       endif
    else
       if (masterproc) write(iulog,*) trim(subname)//': alatlon/llatlon different ', &
          'sizes:continue'
    endif

    ! Exit early if no valid land points
    if ( all(amask == 0) )then
       if (masterproc) write(iulog,*) trim(subname)//': no valid land points do NOT run clm'
       noland = .true.
       return
    end if

    ! ------------------------------------------------------------------------
    ! Initialize clump and processor decomposition
    ! ------------------------------------------------------------------------

    call decompInit_atm(alatlon,amask)

    if (downscale) then
       call map_setgatmfm(gatm,alatlon,llatlon,amask,pftm)
    else
       nlg = alatlon%ns
       allocate(gatm(nlg))
       gatm(:) = -1
       do n = 1,nlg
          if (amask(n) /= 0) then
             gatm(n)=n
          end if
       end do
    end if

    call decompInit_lnd(alatlon%ns, alatlon%ni, alatlon%nj, &
                        llatlon%ns, llatlon%ni, llatlon%nj)

    ! ------------------------------------------------------------------------
    ! Determine atm and lnd domains and surface properties
    ! ------------------------------------------------------------------------

    !--- Read atm grid -----------------------------------------------------

    call get_proc_bounds_atm(begg_atm, endg_atm)

    if (masterproc) then
       write(iulog,*) 'Attempting to read adomain from fatmgrid'
       call shr_sys_flush(iulog)
    endif
    call surfrd_get_grid(adomain, fatmgrid, begg_atm, endg_atm, gratm)

    if (masterproc) then
       write(iulog,*) 'Attempting to read atm landfrac from fatmlndfrc'
       call shr_sys_flush(iulog)
    endif

    if (create_glacier_mec_landunit) then
       call surfrd_get_frac(adomain, fatmlndfrc, fglcmask)

       ! Make sure the glc mask is a subset of the land mask
       do n = begg_atm, endg_atm
          if (adomain%glcmask(n)==1 .and. adomain%mask(n)==0) then
             write(iulog,*) 'landmask/glcmask mismatch'
             write(iulog,*) 'glc requires input where landmask = 0, gridcell index', n
             call shr_sys_flush(iulog)
             call endrun( trim(subname)//'ERROR: clm and land-ice masks are mismatched' )
          endif
       enddo
    else
       call surfrd_get_frac(adomain, fatmlndfrc)
    endif

    if (fatmtopo /= " ") then
       if (masterproc) then
          write(iulog,*) 'Attempting to read atm topo from fatmtopo'
          call shr_sys_flush(iulog)
       endif
       call surfrd_get_topo(adomain, fatmtopo)
    endif

    !--- compute area
    if (.not. adomain%areaset) then
       do nr = begg_atm,endg_atm
          n = adecomp%gdc2glo(nr)
          i = mod(n-1,adomain%ni) + 1
          j = (n-1)/adomain%ni + 1
          adomain%area(nr) = cellarea(alatlon,i,j)
       enddo
       adomain%areaset = .true.
    endif

    if (masterproc) then
       write(iulog,*) 'adomain status:'
       call domain_check(adomain)
    endif

    !--- Read lnd grid -----------------------------------------------------

    call get_proc_bounds(begg, endg)

    if (masterproc) then
       write(iulog,*) 'Attempting to read ldomain from fsurdat ',trim(fsurdat)
       call shr_sys_flush(iulog)
    endif
    ! Note this sets ldomain%pftm
    call surfrd_get_grid(ldomain, fsurdat, begg, endg, grlnd)

    if (flndtopo /= " ") then
       if (masterproc) then
          write(iulog,*) 'Attempting to read lnd topo from flndtopo ',trim(flndtopo)
          call shr_sys_flush(iulog)
       endif
       ! The following sets ldomain%topo
       call surfrd_get_topo(ldomain, flndtopo)
    endif

    !--- compute area
    if (.not. ldomain%areaset) then
       do nr = begg,endg
          n = ldecomp%gdc2glo(nr)
          i = mod(n-1,ldomain%ni) + 1
          j = (n-1)/ldomain%ni + 1
          ldomain%area(nr) = cellarea(alatlon,i,j)
       enddo
       ldomain%areaset = .true.
    endif

    if (masterproc) then
       call domain_check(ldomain)
    endif

    !--- overwrite ldomain if same grids -----------------------------------

    if (.not. downscale) then
       if (masterproc) write(iulog,*) trim(subname)//': downscale false, set ldomain =~ adomain'
       ! Don't copy [pftm, topo, mask, frac, nara, ntop]
       ! ldomain%pftm is set above in call to surfrd_get_grid
       ! ldomain%topo is set above in call to surfrd_get_topo
       ! ldomain%mask, ldomain%frac, ldomain%nara and ndomain%ntop are set in initialize2 

       ldomain%latc     = adomain%latc
       ldomain%lonc     = adomain%lonc
       ldomain%area     = adomain%area
       ldomain%set      = adomain%set
       ldomain%areaset  = adomain%areaset 
       ldomain%regional = adomain%regional
       ldomain%decomped = adomain%decomped
       ldomain%glcmask  = adomain%glcmask
    endif

    ! --------------------------------------------------------------------
    ! Initialize urban model input (initialize urbinp data structure)
    ! --------------------------------------------------------------------

    call UrbanInput(mode='initialize')

    ! Allocate surface grid dynamic memory (for wtxy and vegxy arrays)

    call t_stopf('init_grids')
    call t_startf('init_surdat')

    call get_proc_bounds    (begg    , endg)
    allocate (vegxy(begg:endg,maxpatch), &
              wtxy(begg:endg,maxpatch),  &
              stat=ier)   
    if (ier /= 0) then
       write(iulog,*)'initialize allocation error'; call endrun()
    endif

    ! Allocate additional dynamic memory for glacier_mec topo and thickness

    if (create_glacier_mec_landunit) then
       allocate (topoxy(begg:endg,maxpatch), stat=ier)
       if (ier /= 0) then
          write(iulog,*)'initialize allocation error'; call endrun()
       endif
    else
       allocate (topoxy(1,1), stat=ier)
       if (ier /= 0) then
          write(iulog,*)'initialize allocation error'; call endrun()
       endif
    endif

    ! --------------------------------------------------------------------
    ! Read list of PFTs and their corresponding parameter values
    ! Independent of model resolution
    ! Needs to stay before call surfrd
    ! --------------------------------------------------------------------

    call pftconrd()

    ! Read surface dataset and set up vegetation type [vegxy] and 
    ! weight [wtxy] arrays for [maxpatch] subgrid patches.
    
    call surfrd (fsurdat, ldomain)

    if (create_glacier_mec_landunit) then
       call decompInit_glcp (alatlon%ns,   alatlon%ni,   alatlon%nj, &
                             llatlon%ns,   llatlon%ni,   llatlon%nj, &
                             ldomain%glcmask)
    else
       call decompInit_glcp (alatlon%ns, alatlon%ni, alatlon%nj, &
                             llatlon%ns, llatlon%ni, llatlon%nj)
    endif

    call t_stopf('init_surdat')

  end subroutine initialize1

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: initialize2
!
! !INTERFACE:
  subroutine initialize2( )
!
! !DESCRIPTION:
! Land model initialization.
! o Initializes run control variables via the [clm_inparm] namelist.
! o Reads surface data on model grid.
! o Defines the multiple plant types and fraction areas for each surface type.
! o Builds the appropriate subgrid <-> grid mapping indices and weights.
! o Set up parallel processing.
! o Initializes time constant variables.
! o Reads restart data for a restart or branch run.
! o Reads initial data and initializes the time variant variables for an initial run.
! o Initializes history file output.
! o Initializes river routing model.
! o Initializes accumulation variables.
!
! !USES:
    use clm_atmlnd      , only : init_atm2lnd_type, init_lnd2atm_type, &
                                 clm_a2l, clm_l2a, atm_a2l, atm_l2a, &
	                         init_adiag_type
    use initGridCellsMod, only : initGridCells
    use clm_varctl      , only : finidat, fpftdyn
    use clmtypeInitMod  , only : initClmtype
    use domainMod       , only : gatm, ldomain, adomain, llatlon, alatlon
    use decompMod       , only : adecomp,ldecomp, &
                                 get_proc_clumps, get_clump_bounds, &
                                 get_proc_bounds, get_proc_bounds_atm
    use downscaleMod    , only : map_setmapsFM, map1dl_a2l, map1dl_l2a
    use filterMod       , only : allocFilters, setFilters
    use histFldsMod     , only : hist_initFlds
    use histFileMod     , only : hist_htapes_build
    use restFileMod     , only : restFile_getfile, &
                                 restFile_open, restFile_close, restFile_read 
    use accFldsMod      , only : initAccFlds, initAccClmtype
    use mkarbinitMod    , only : mkarbinit
    use pftdynMod       , only : pftdyn_init, pftdyn_interp
#ifdef CN
    use ndepStreamMod    , only : ndep_init, ndep_interp
    use CNEcosystemDynMod, only : CNEcosystemDynInit
#endif
#if (defined CNDV)
    use pftdynMod             , only : pftwt_init
    use CNDVEcosystemDyniniMod, only : CNDVEcosystemDynini
#endif
    use STATICEcosysDynMod , only : EcosystemDynini, readAnnualVegetation
    use STATICEcosysDynMod , only : interpMonthlyVeg
    use DustMod         , only : Dustini
#if (defined CASA)
    use CASAMod         , only : initCASA
    use CASAPhenologyMod, only : initCASAPhenology
    use CASAiniTimeVarMod,only : CASAiniTimeVar
#endif
#if (defined RTM) 
    use RtmMod          , only : Rtmini
#endif
    use clm_time_manager, only : get_curr_date, get_nstep, advance_timestep, &
                                 timemgr_init, timemgr_restart_io, timemgr_restart
    use clm_time_manager, only : get_step_size, get_curr_calday
    use fileutils       , only : getfil
    use UrbanMod        , only : UrbanClumpInit
    use UrbanInitMod    , only : UrbanInitTimeConst, UrbanInitTimeVar, UrbanInitAero 
    use UrbanInputMod   , only : UrbanInput
    use clm_glclnd      , only : init_glc2lnd_type, init_lnd2glc_type, &
                                 clm_x2s, clm_s2x, atm_x2s, atm_s2x
    use seq_drydep_mod  , only : n_drydep, drydep_method, DD_XLND
    use shr_orb_mod        , only : shr_orb_decl
    use initSurfAlbMod     , only : initSurfAlb, do_initsurfalb 
    use clm_atmlnd         , only : clm_map2gcell
    use clm_varorb         , only : eccen, mvelpp, lambm0, obliqr

! !Arguments    
    implicit none
!
! !REVISION HISTORY:
! Created by Gordon Bonan, Sam Levis and Mariana Vertenstein
!
!
! !LOCAL VARIABLES:
!EOP
    integer  :: nl,nlg,na,nag         ! indices
    integer  :: i,j,k,n1,n2           ! indices
    integer  :: yr                    ! current year (0, ...)
    integer  :: mon                   ! current month (1 -> 12)
    integer  :: day                   ! current day (1 -> 31)
    integer  :: ncsec                 ! current time of day [seconds]
    integer  :: nc                    ! clump index
    integer  :: nclumps               ! number of clumps on this processor
    integer  :: begp, endp            ! clump beg and ending pft indices
    integer  :: begc, endc            ! clump beg and ending column indices
    integer  :: begl, endl            ! clump beg and ending landunit indices
    integer  :: begg, endg            ! clump beg and ending gridcell indices
    integer  :: begg_atm, endg_atm    ! proc beg and ending gridcell indices
    character(len=256) :: fnamer      ! name of netcdf restart file 
    character(len=256) :: pnamer      ! full pathname of netcdf restart file
    type(file_desc_t)  :: ncid        ! netcdf id
    real(r8) :: dtime                 ! time step increment (sec)
    integer  :: nstep                 ! model time step
    real(r8) :: calday                ! calendar day for nstep
    real(r8) :: caldaym1              ! calendar day for nstep-1
    real(r8) :: declin                ! solar declination angle in radians for nstep
    real(r8) :: declinm1              ! solar declination angle in radians for nstep-1
    real(r8) :: eccf                  ! earth orbit eccentricity factor
    character(len=32) :: subname = 'initialize2' ! subroutine name
!----------------------------------------------------------------------

    call t_startf('clm_init2')

    ! ------------------------------------------------------------------------
    ! Set the a2l and l2a maps
    ! ------------------------------------------------------------------------

    call get_proc_bounds(begg=begg,endg=endg)

    if (downscale) then
       call t_startf('init_mapsFM')
       call map_setmapsFM(adomain, ldomain, gatm)
       call t_stopf('init_mapsFM')
       
       ! Set ldomain mask and frac based on adomain and mapping.
       ! Want ldomain%frac to match adomain%frac but scale by effective area.
       ! so the implied area of an ldomain cell is the actual area * scaled frac 
       ! which aggregated over all land cells under an atm cell,
       ! will match the area associated with the atm cell.
       
       call get_proc_bounds_atm(begg=begg_atm,endg=endg_atm)
       do nl = begg,endg
          nlg = ldecomp%gdc2glo(nl)
          nag = gatm(nlg)
          if (nag <= 0) then
             write(iulog,*) 'nag <= 0 ERROR ',nl,nlg,nag
             call endrun( trim(subname)//'ERROR: no grid cells' )
          endif
          na = adecomp%glo2gdc(nag)
          if (nlg < 1        .or. nlg > llatlon%ns .or. &
              nag < 1        .or. nag > alatlon%ns .or. &
              na  < begg_atm .or. na  > endg_atm) then
             write(iulog,*) 'index ERROR ',&
                  nl,nlg,nag,na,llatlon%ns,alatlon%ns,begg_atm,endg_atm
             call endrun( trim(subname)//'ERROR: no grid cells' )
          endif
          ldomain%mask(nl) = 1
          ldomain%frac(nl) = adomain%frac(na) * (ldomain%nara(nl)/ldomain%area(nl))
       enddo
    else
       do nl = begg,endg
          ldomain%mask(nl) = 1
          ldomain%frac(nl) = adomain%frac(nl) 
	  ldomain%nara(nl) = adomain%area(nl)
	  ldomain%ntop(nl) = adomain%topo(nl)
       end do
    end if

    ! ------------------------------------------------------------------------
    ! Allocate memory and initialize values of clmtype data structures
    ! ------------------------------------------------------------------------

    call t_startf('init_clmtype')
    call initClmtype()

    call get_proc_bounds    (begg    , endg    , begl, endl, &
                             begc    , endc    , begp, endp )
    call init_atm2lnd_type  (begg    , endg    , clm_a2l)
    call init_lnd2atm_type  (begg    , endg    , clm_l2a)

    call get_proc_bounds_atm(begg_atm, endg_atm)
    call init_atm2lnd_type  (begg_atm, endg_atm, atm_a2l)
    call init_lnd2atm_type  (begg_atm, endg_atm, atm_l2a)

    call init_adiag_type()

    if (create_glacier_mec_landunit) then
       call init_glc2lnd_type  (begg    , endg    , clm_x2s)
       call init_lnd2glc_type  (begg    , endg    , clm_s2x)

       call init_glc2lnd_type  (begg    , endg    , atm_x2s)
       call init_lnd2glc_type  (begg    , endg    , atm_s2x)
    endif

    ! ------------------------------------------------------------------------
    ! Build hierarchy and topological info for derived typees
    ! ------------------------------------------------------------------------

    call initGridCells()

    ! Deallocate surface grid dynamic memory (for wtxy and vegxy arrays)

    deallocate (vegxy, wtxy)
    deallocate (topoxy)

    call t_stopf('init_clmtype')

    ! ------------------------------------------------------------------------
    ! Initialize time constant variables 
    ! ------------------------------------------------------------------------

    ! Initialize Ecosystem Dynamics 

    call t_startf('init_ecosys')
#if (defined CNDV)
    call CNDVEcosystemDynini()
#elif (!defined CN)
    call EcosystemDynini()
#endif
#if (defined CN) || (defined CNDV)

    ! --------------------------------------------------------------
    ! Initialize CLMSP ecosystem dynamics when drydeposition is used
    ! so that estimates of monthly differences in LAI can be computed
    ! --------------------------------------------------------------

    if ( n_drydep > 0 .and. drydep_method == DD_XLND )then
       call EcosystemDynini()
    end if
#endif
    call t_stopf('init_ecosys')

    ! Initialize dust emissions model 

    call t_startf('init_dust')
    call Dustini()
    call t_stopf('init_dust')

    ! ------------------------------------------------------------------------
    ! Initialize time constant urban variables
    ! (this must be called before initCASA())
    ! ------------------------------------------------------------------------

    call t_startf('init_io1')
    call UrbanInitTimeConst()
    call iniTimeConst()

    ! ------------------------------------------------------------------------
    ! Obtain restart file if appropriate
    ! ------------------------------------------------------------------------

    if (do_restread()) then
       call restFile_getfile( file=fnamer, path=pnamer )
    end if

    ! ------------------------------------------------------------------------
    ! Initialize time manager
    ! ------------------------------------------------------------------------

    if (nsrest == nsrStartup) then  
       call timemgr_init()
    else
       call restFile_open( flag='read', file=fnamer, ncid=ncid )
       call timemgr_restart_io( ncid=ncid, flag='read' )
       call restFile_close( ncid=ncid )
       call timemgr_restart()
    end if
    call t_stopf('init_io1')

    ! ------------------------------------------------------------------------
    ! Initialize CN Ecosystem Dynamics (must be after time-manager initialization)
    ! ------------------------------------------------------------------------
#if (defined CN) || (defined CNDV)
    call CNEcosystemDynInit( begc, endc, begp, endp )
#endif

    ! --------------------------------------------------------------
    ! Initialize river routing model
    ! --------------------------------------------------------------

#if (defined RTM)
    if (masterproc) write(iulog,*)'Attempting to initialize RTM'
    call shr_sys_flush(iulog)
    call t_startf('init_rtm')
    call Rtmini()
    call t_stopf('init_rtm')
    if (masterproc) write(iulog,*)'Successfully initialized RTM'
#endif

    ! ------------------------------------------------------------------------
    ! Initialize accumulated fields
    ! ------------------------------------------------------------------------

    ! Initialize accumulator fields to be time accumulated for various purposes.
    ! The time manager needs to be initialized before this called is made, since
    ! the step size is needed.

    call t_startf('init_accflds')
    call initAccFlds()
    call t_stopf('init_accflds')

    ! ------------------------------------------------------------------------
    ! Set arbitrary initial conditions for time varying fields 
    ! used in coupled carbon-nitrogen code
    ! ------------------------------------------------------------------------
    
#if (defined CN)
    call t_startf('init_cninitim')
    if (nsrest == nsrStartup) then
       call CNiniTimeVar()
    end if
    call t_stopf('init_cninitim')
#elif (defined CASA)
    call t_startf('init_casainitim')
    call CASAiniTimeVar()
    call t_stopf('init_casainitim')
#endif

    ! ------------------------------------------------------------------------
    ! Initialization of dynamic pft weights
    ! ------------------------------------------------------------------------

    ! Determine correct pft weights (interpolate pftdyn dataset if initial run)
    ! Otherwise these are read in for a restart run

#if (defined CNDV)
    call pftwt_init()
#else
    if (fpftdyn /= ' ') then
       call t_startf('init_pftdyn')
       call pftdyn_init()
       call pftdyn_interp( )
       call t_stopf('init_pftdyn')
    end if
#endif


    ! ------------------------------------------------------------------------
    ! Read restart/initial info
    ! ------------------------------------------------------------------------

    ! No weight related information can be contained in the routines,  
    ! "mkarbinit, inicfile and restFile". 

    call t_startf('init_io2')
    if (do_restread()) then
       if (masterproc) write(iulog,*)'reading restart file ',fnamer
       call restFile_read( fnamer )
    else if (nsrest == nsrStartup .and. finidat == ' ') then
       call mkarbinit()
       call UrbanInitTimeVar( )
    end if
    call t_stopf('init_io2')

    ! ------------------------------------------------------------------------
    ! Initialize nitrogen deposition
    ! ------------------------------------------------------------------------

#ifdef CN
    call t_startf('init_ndep')
    call ndep_init()
    call ndep_interp()
    call t_stopf('init_ndep')
#endif
    
    ! ------------------------------------------------------------------------
    ! Initialization of model parameterizations that are needed after
    ! restart file is read in
    ! ------------------------------------------------------------------------

#if (defined CASA)
    call t_startf('init_casa')
    ! Initialize CASA 

    call initCASA()
    call initCASAPhenology()
    call t_stopf('init_casa')
#endif

    ! ------------------------------------------------------------------------
    ! Initialize history and accumator buffers
    ! ------------------------------------------------------------------------

    call t_startf('init_hist1')
    ! Initialize master history list. 

    call hist_initFlds()

    ! Initialize active history fields. This is only done if not a restart run. 
    ! If a restart run, then this information has already been obtained from the 
    ! restart data read above. Note that routine hist_htapes_build needs time manager 
    ! information, so this call must be made after the restart information has been read.

    if (nsrest == nsrStartup .or. nsrest == nsrBranch) call hist_htapes_build()

    ! Initialize clmtype variables that are obtained from accumulated fields.
    ! This routine is called in an initial run at nstep=0
    ! This routine is also always called for a restart run and must 
    ! therefore be called after the restart file is read in

    call initAccClmtype()

    call t_stopf('init_hist1')

    ! --------------------------------------------------------------
    ! Note - everything below this point needs updated weights
    ! --------------------------------------------------------------

    ! Initialize filters
    
    call t_startf('init_filters')

    call allocFilters()
    nclumps = get_proc_clumps()
!$OMP PARALLEL DO PRIVATE (nc)
    do nc = 1, nclumps
       call setFilters(nc)
    end do
!$OMP END PARALLEL DO

    call t_stopf('init_filters')

    ! Calculate urban "town" roughness length and displacement 
    ! height for urban landunits

    call UrbanInitAero()

    ! Initialize urban radiation model - this uses urbinp data structure

    call UrbanClumpInit()

    ! Finalize urban model initialization
    
    call UrbanInput(mode='finalize')

    !
    ! Even if CN is on, and dry-deposition is active, read CLMSP annual vegetation to get estimates of monthly LAI
    !
    if ( n_drydep > 0 .and. drydep_method == DD_XLND )then
       call readAnnualVegetation()
    end if

    ! End initialization

    call t_startf('init_wlog')
    if (masterproc) then
       write(iulog,*) 'Successfully initialized the land model'
       if (nsrest == nsrStartup) then
          write(iulog,*) 'begin initial run at: '
       else
          write(iulog,*) 'begin continuation run at:'
       end if
       call get_curr_date(yr, mon, day, ncsec)
       write(iulog,*) '   nstep= ',get_nstep(), ' year= ',yr,' month= ',mon,&
            ' day= ',day,' seconds= ',ncsec
       write(iulog,*)
       write(iulog,'(72a1)') ("*",i=1,60)
       write(iulog,*)
    endif
    call t_stopf('init_wlog')

    if (get_nstep() == 0 .or. nsrest == nsrStartup) then
       ! Initialize albedos (correct pft filters are needed)

       if (finidat == ' ' .or. do_initsurfalb) then
          call t_startf('init_orb')
          calday = get_curr_calday()
          call t_startf('init_orbd1')
          call shr_orb_decl( calday, eccen, mvelpp, lambm0, obliqr, declin, eccf )
          call t_stopf('init_orbd1')
          
          dtime = get_step_size()
          caldaym1 = get_curr_calday(offset=-int(dtime))
          call t_startf('init_orbd2')
          call shr_orb_decl( caldaym1, eccen, mvelpp, lambm0, obliqr, declinm1, eccf )
          call t_stopf('init_orbd2')
          
          call t_startf('init_orbSA')
          call initSurfAlb( calday, declin, declinm1 )
          call t_stopf('init_orbSA')
          call t_stopf('init_orb')
       else if ( n_drydep > 0 .and. drydep_method == DD_XLND )then
          ! Call interpMonthlyVeg for dry-deposition so that mlaidiff will be calculated
          ! This needs to be done even if CN or CNDV is on!
          call interpMonthlyVeg()
       end if

       ! Determine gridcell averaged properties to send to atm

       call t_startf('init_map2gc')
       call clm_map2gcell(init=.true.)
       call t_stopf('init_map2gc')

    end if
    call t_stopf('clm_init2')

  end subroutine initialize2

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: header
!
! !INTERFACE:
  subroutine header()
!
! !DESCRIPTION:
! Echo and save model version number
!
! !USES:
    use clm_varctl  , only : version
!
! !ARGUMENTS:
    implicit none
!
! !CALLED FROM:
! subroutine initialize in this module
!
! !REVISION HISTORY:
! Created by Gordon Bonan
!
!EOP
!-----------------------------------------------------------------------

    if ( masterproc )then
      write(iulog,*) trim(version)
      write(iulog,*)
    end if

  end subroutine header

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: do_restread
!
! !INTERFACE:
  logical function do_restread( )
!
! !DESCRIPTION:
! Determine if restart file will be read
!
! !USES:
    use clm_varctl, only : finidat
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

    do_restread = .false.
    if (nsrest == nsrStartup .and. finidat /= ' ') then
       do_restread = .true.
    end if
    if (nsrest == nsrContinue .or. nsrest == nsrBranch) then
       do_restread = .true.
    end if
  end function do_restread
  
!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cellarea
!
! !INTERFACE:
  real(r8) function cellarea(latlon,i,j)
!
! !DESCRIPTION:
! Comute area of grid cells (square kilometers)

! !USES:
    use domainMod , only : latlon_type
    use clm_varcon, only : re
!
! !ARGUMENTS:
    implicit none
    type(latlon_type), intent(in   ) :: latlon    ! latlon datatype
    integer          , intent(in)    :: i,j       ! latlon index to compute
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!
! !LOCAL VARIABLES:
!EOP
    real(r8) :: deg2rad ! pi/180
    real(r8) :: dx      ! cell width: E-W
    real(r8) :: dy      ! cell width: N-S
!------------------------------------------------------------------------

    deg2rad = SHR_CONST_PI / 180._r8
    dx = (latlon%lone(i) - latlon%lonw(i)) * deg2rad
    dy = sin(latlon%latn(j)*deg2rad) - sin(latlon%lats(j)*deg2rad)
    cellarea = dx*dy*re*re

  end function cellarea

end module clm_initializeMod
