module clm_initializeMod

  !-----------------------------------------------------------------------
  ! Performs land model initialization
  !
  use shr_kind_mod    , only : r8 => shr_kind_r8
  use spmdMod         , only : masterproc
  use shr_sys_mod     , only : shr_sys_flush, shr_sys_abort
  use abortutils      , only : endrun
  use clm_varctl      , only : nsrest, nsrStartup, nsrContinue, nsrBranch, &
                               create_glacier_mec_landunit, iulog, use_lch4, use_cn, &
                               use_cndv
  use clm_varsur      , only : wt_lunit, wt_nat_pft, wt_cft, wt_glc_mec, topo_glc_mec
  use perf_mod        , only : t_startf, t_stopf
  use readParamsMod   , only : readParameters
  use ncdio_pio
  use mct_mod

  implicit none
  save
  private    ! By default everything is private

  public :: initialize1  ! Phase one initialization
  public :: initialize2  ! Phase two initialization
  !
  private header         ! echo version numbers
  private do_restread    ! read a restart file
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
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
    use clmtypeInitMod   ,only: initClmtype
    use clm_varpar       ,only: clm_varpar_init, natpft_lb, natpft_ub, cft_lb, cft_ub, maxpatch_glcmec
    use clm_varcon       ,only: clm_varcon_init, max_lunit
    use clm_varctl       ,only: fsurdat, fatmlndfrc, flndtopo, fglcmask, noland 
    use pftvarcon        ,only: pftconrd
    use decompInitMod    ,only: decompInit_lnd, decompInit_clumps, decompInit_glcp
    use decompMod        ,only: bounds_type, get_proc_bounds 
    use domainMod        ,only: domain_check, ldomain, domain_init
    use surfrdMod        ,only: surfrd_get_globmask, surfrd_get_grid, surfrd_get_topo, surfrd_get_data 
    use controlMod       ,only: control_init, control_print, nlfilename
    use UrbanInputMod    ,only: UrbanInput
    use ncdio_pio        ,only: ncd_pio_init
    use clm_atmlnd       ,only: init_atm2lnd_type, init_lnd2atm_type, clm_a2l, clm_l2a
    use clm_glclnd       ,only: init_glc2lnd_type, init_lnd2glc_type, clm_x2s, clm_s2x
    use initGridCellsMod ,only: initGridCells
    use ch4varcon        ,only: ch4conrd
    !
    ! !LOCAL VARIABLES:
    integer  :: ier                   ! error status
    integer  :: i,j,n,k               ! loop indices
    integer  :: nl                    ! gdc and glo lnd indices
    integer  :: ns, ni, nj            ! global grid sizes
    integer  :: begg, endg            ! processor bounds
    type(bounds_type) :: bounds       ! bounds
    integer ,pointer  :: amask(:)     ! global land mask
    character(len=32) :: subname = 'initialize1' ! subroutine name
    !-----------------------------------------------------------------------

    call t_startf('clm_init1')

#if (defined _OPENMP)
    call shr_sys_abort( &
         'ERROR: OpenMP does NOT work on the CLM4.5 science branch. &
         Set the number of threads for LND to 1 and rerun.' )
#endif

    ! ------------------------------------------------------------------------
    ! Initialize run control variables, timestep
    ! ------------------------------------------------------------------------

    call header()

    if (masterproc) then
       write(iulog,*) 'Attempting to initialize the land model .....'
       write(iulog,*)
       call shr_sys_flush(iulog)
    endif

    call control_init()
    call clm_varpar_init()
    call clm_varcon_init()
    call ncd_pio_init()


    if (masterproc) call control_print()

    ! ------------------------------------------------------------------------
    ! Read in global land grid and land mask (amask)- needed to set decomposition
    ! ------------------------------------------------------------------------

    ! global memory for amask is allocate in surfrd_get_glomask - must be
    ! deallocated below
    if (masterproc) then
       write(iulog,*) 'Attempting to read global land mask from ',trim(fatmlndfrc)
       call shr_sys_flush(iulog)
    endif
    call surfrd_get_globmask(filename=fatmlndfrc, mask=amask, ni=ni, nj=nj)

    ! Exit early if no valid land points
    if ( all(amask == 0) )then
       if (masterproc) write(iulog,*) trim(subname)//': no valid land points do NOT run clm'
       noland = .true.
       return
    end if

    ! ------------------------------------------------------------------------
    ! Determine clm gridcell decomposition and processor bounds for gridcells
    ! ------------------------------------------------------------------------

    call decompInit_lnd(ni, nj, amask)
    deallocate(amask)

    ! *** Get JUST gridcell processor bounds ***
    ! Remaining bounds (landunits, columns, pfts) will be determined 
    ! after the call to decompInit_glcp - so get_proc_bounds is called
    ! twice and the gridcell information is just filled in twice

    call get_proc_bounds(begg, endg)

    ! ------------------------------------------------------------------------
    ! Get grid and land fraction (set ldomain)
    ! ------------------------------------------------------------------------

    if (masterproc) then
       write(iulog,*) 'Attempting to read ldomain from ',trim(fatmlndfrc)
       call shr_sys_flush(iulog)
    endif
    if (create_glacier_mec_landunit) then
       call surfrd_get_grid(begg, endg, ldomain, fatmlndfrc, fglcmask)
    else
       call surfrd_get_grid(begg, endg, ldomain, fatmlndfrc)
    endif
    if (masterproc) then
       call domain_check(ldomain)
    endif
    ldomain%mask = 1  !!! TODO - is this needed?

    ! Get topo if appropriate (set ldomain%topo)

    if (flndtopo /= " ") then
       if (masterproc) then
          write(iulog,*) 'Attempting to read atm topo from ',trim(flndtopo)
          call shr_sys_flush(iulog)
       endif
       call surfrd_get_topo(ldomain, flndtopo)  
    endif

    ! Initialize urban model input (initialize urbinp data structure)

    call UrbanInput(begg, endg, mode='initialize')

    ! Allocate surface grid dynamic memory (just gridcell bounds dependent)

    allocate (wt_lunit  (begg:endg, max_lunit))
    allocate (wt_nat_pft(begg:endg, natpft_lb:natpft_ub))
    allocate (wt_cft    (begg:endg, cft_lb:cft_ub))
    if (create_glacier_mec_landunit) then
       allocate (wt_glc_mec  (begg:endg, maxpatch_glcmec))
       allocate (topo_glc_mec(begg:endg, maxpatch_glcmec))
    else
       allocate (wt_glc_mec  (1,1))
       allocate (topo_glc_mec(1,1))
    endif

    ! Read list of PFTs and their corresponding parameter values
    ! Independent of model resolution, Needs to stay before surfrd_get_data

    call pftconrd()

    ! Read surface dataset and set up subgrid weight arrays
    
    call surfrd_get_data(begg, endg, ldomain, fsurdat)

    ! ------------------------------------------------------------------------
    ! Determine decomposition of subgrid scale landunits, columns, pfts
    ! ------------------------------------------------------------------------

    if (create_glacier_mec_landunit) then
       call decompInit_clumps(ns, ni, nj, ldomain%glcmask)
    else
       call decompInit_clumps(ns, ni, nj)
    endif

    ! *** Get ALL processor bounds - for gridcells, landunit, columns and pfts ***

    call get_proc_bounds(bounds)
    
    ! Allocate memory and initialize values of clmtype data structures
    ! This is needed here for the following call to initGridcells

    call initClmtype(bounds)

    ! Build hierarchy and topological info for derived types
    ! This is needed here for the following call to decompInit_glcp

    call initGridCells(bounds)

    ! Set global seg maps for gridcells, landlunits, columns and pfts

    if (create_glacier_mec_landunit) then
       call decompInit_glcp(ns, ni, nj, ldomain%glcmask)
    else
       call decompInit_glcp(ns, ni, nj)
    endif

    ! Initialize atm->lnd, lnd->atm, glc->lnd and lnd->glc data structures

    call init_atm2lnd_type(bounds, clm_a2l)
    call init_lnd2atm_type(bounds, clm_l2a)
    if (create_glacier_mec_landunit) then
       call init_glc2lnd_type(bounds, clm_x2s)
       call init_lnd2glc_type(bounds, clm_s2x)
    endif

    ! ------------------------------------------------------------------------
    ! Remainder of initialization1
    ! ------------------------------------------------------------------------

    ! Set CH4 Model Parameters from namelist.
    ! Need to do before iniTimeConst so that it knows whether to 
    ! look for several optional parameters on surfdata file.

    if (use_lch4) then
       call ch4conrd()
    end if

    ! Deallocate surface grid dynamic memory for variables that aren't needed elsewhere
    ! Note that wt_lunit is kept until the end of initialize2 so we can do some
    ! consistency checking. 

    deallocate (wt_nat_pft, wt_cft, wt_glc_mec, topo_glc_mec)

    call t_stopf('clm_init1')

  end subroutine initialize1


  !-----------------------------------------------------------------------
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
    use clm_atmlnd            , only : clm_map2gcell
    use clm_glclnd            , only : update_clm_s2x
    use clm_varctl            , only : finidat, fpftdyn
    use decompMod             , only : get_proc_clumps, get_proc_bounds, bounds_type
    use filterMod             , only : allocFilters
    use reweightMod           , only : reweightWrapup
    use histFldsMod           , only : hist_initFlds
    use histFileMod           , only : hist_htapes_build, htapes_fieldlist
    use restFileMod           , only : restFile_getfile, restFile_open, restFile_close, restFile_read 
    use accFldsMod            , only : initAccFlds, initAccClmtype
    use mkarbinitMod          , only : mkarbinit
    use pftdynMod             , only : pftdyn_init, pftdyn_interp
    use ndepStreamMod         , only : ndep_init, ndep_interp
    use CNEcosystemDynMod     , only : CNEcosystemDynInit
    use pftdynMod             , only : pftwt_init
    use CNDVEcosystemDyniniMod, only : CNDVEcosystemDynini

    use STATICEcosysDynMod    , only : EcosystemDynini, readAnnualVegetation, interpMonthlyVeg
    use DustMod               , only : Dustini
    use clm_time_manager      , only : get_step_size, get_curr_calday
    use fileutils             , only : getfil
    use UrbanMod              , only : UrbanParamInit
    use UrbanInitMod          , only : UrbanInitTimeConst, UrbanInitTimeVar, UrbanInitAero 
    use UrbanInputMod         , only : UrbanInput
    use initSLakeMod          , only : initSLake
    use initch4Mod            , only : initch4
    use clm_glclnd            , only : init_glc2lnd_type, init_lnd2glc_type, clm_x2s, clm_s2x
    use seq_drydep_mod        , only : n_drydep, drydep_method, DD_XLND
    use shr_orb_mod           , only : shr_orb_decl
    use initSurfAlbMod        , only : initSurfAlb, do_initsurfalb 
    use clm_varorb            , only : eccen, mvelpp, lambm0, obliqr
    use VOCEmissionMod        , only : VOCEmission_init
    use clm_time_manager      , only : get_curr_date, get_nstep, advance_timestep, &
                                       timemgr_init, timemgr_restart_io, timemgr_restart
    !
    ! !ARGUMENTS    
    implicit none
    !
    ! !LOCAL VARIABLES:
    integer  :: nl,na,nag             ! indices
    integer  :: i,j,k                 ! indices
    integer  :: yr                    ! current year (0, ...)
    integer  :: mon                   ! current month (1 -> 12)
    integer  :: day                   ! current day (1 -> 31)
    integer  :: ncsec                 ! current time of day [seconds]
    integer  :: nc                    ! clump index
    integer  :: nclumps               ! number of clumps on this processor
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
    type(bounds_type) :: bounds       ! bounds
    character(len=32) :: subname = 'initialize2' ! subroutine name
    logical           :: arbinit            ! Make arb init in initSLake
    !----------------------------------------------------------------------

    ! ------------------------------------------------------------------------
    ! Determine processor bounds and clumps for this processor
    ! ------------------------------------------------------------------------

    call get_proc_bounds(bounds)
    nclumps = get_proc_clumps()

    call t_startf('clm_init2')

    ! ------------------------------------------------------------------------
    ! Initialize time constant variables 
    ! ------------------------------------------------------------------------

    ! Initialize Ecosystem Dynamics 

    call t_startf('init_ecosys')
    if (use_cndv) then
       call CNDVEcosystemDynini(bounds)
    else if (.not. use_cn) then
       call EcosystemDynini(bounds)
    end if

    ! Initialize CLMSP ecosystem dynamics when drydeposition is used
    ! so that estimates of monthly differences in LAI can be computed
    if (use_cn .or. use_cndv) then
       if ( n_drydep > 0 .and. drydep_method == DD_XLND )then
          call EcosystemDynini(bounds)
       end if
    end if
    call t_stopf('init_ecosys')

    ! Initialize dust emissions model 

    call t_startf('init_dust')
    call Dustini(bounds)
    call t_stopf('init_dust')
    
    ! Initialize MEGAN emissions model 

    call VOCEmission_init()

    ! Initialize time constant urban variables

    call t_startf('init_io1')

    ! ------------------------------------------------------------------------
    ! Read in constants files
    ! ------------------------------------------------------------------------
    call readParameters()

    call UrbanInitTimeConst(bounds)
    call iniTimeConst(bounds)

    ! ------------------------------------------------------------------------
    ! Obtain restart file if appropriate
    ! ------------------------------------------------------------------------

    if (do_restread()) then
       call restFile_getfile(file=fnamer, path=pnamer)
    end if

    ! ------------------------------------------------------------------------
    ! Initialize master history list. 
    ! ------------------------------------------------------------------------
    call t_startf('hist_initFlds')

    call hist_initFlds()
    ! On restart process the history namelist. Later the namelist from the restart file
    ! will be used. But, this allows some basic checking to make sure you didn't
    ! try to change the history namelist on restart.
    if (nsrest == nsrContinue ) call htapes_fieldlist()

    call t_stopf('hist_initFlds')

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
    if (use_cn .or. use_cndv) then 
       call get_proc_bounds(bounds)
       call CNEcosystemDynInit(bounds)
    end if

    ! ------------------------------------------------------------------------
    ! Initialize accumulated fields
    ! ------------------------------------------------------------------------

    ! Initialize accumulator fields to be time accumulated for various purposes.
    ! The time manager needs to be initialized before this called is made, since
    ! the step size is needed.

    call t_startf('init_accflds')
    call initAccFlds(bounds)
    call t_stopf('init_accflds')

    ! ------------------------------------------------------------------------
    ! Set arbitrary initial conditions for time varying fields 
    ! used in coupled carbon-nitrogen code
    ! ------------------------------------------------------------------------
    
    if (use_cn) then
       call t_startf('init_cninitim')
       if (nsrest == nsrStartup) then
          call CNiniTimeVar(bounds)
       end if
       call t_stopf('init_cninitim')
    end if

    ! ------------------------------------------------------------------------
    ! Initialization of dynamic pft weights
    ! ------------------------------------------------------------------------

    ! Determine correct pft weights (interpolate pftdyn dataset if initial run)
    ! Otherwise these are read in for a restart run

    if (use_cndv) then
       call pftwt_init(bounds)
    else
       if (fpftdyn /= ' ') then
          call t_startf('init_pftdyn')
          call pftdyn_init(bounds)
          call pftdyn_interp(bounds)
          call t_stopf('init_pftdyn')
       end if
    end if

    ! ------------------------------------------------------------------------
    ! Read restart/initial info
    ! ------------------------------------------------------------------------

    ! No weight related information can be contained in the routines,  
    ! "mkarbinit, inicfile and restFile". 

    call t_startf('init_io2')
    if (do_restread()) then
       if (masterproc) write(iulog,*)'reading restart file ',fnamer
       call restFile_read(bounds, fnamer)

       arbinit = .false.
       call initSLake(bounds, arbinit)
       if (use_lch4) then
          arbinit = .false.
          call initch4(bounds, arbinit)
       end if
    else if (nsrest == nsrStartup .and. finidat == ' ') then
       call mkarbinit(bounds)
       call UrbanInitTimeVar(bounds)

       arbinit = .true.
       call initSLake(bounds, arbinit)
       if (use_lch4) then
          arbinit = .true.
          call initch4(bounds, arbinit)
       end if
    !!!!! Attn EK: The calls to initch4 and initSLake combine both setting of state vars, and constant + flux & diagnostic
    !          vars.  This is set up so that the submodels would be back-compatible with old restart files.
    !          It is intended to work and allow bfb restarts as is, but may not be consistent style.
    !          See these two routines for structure.  Feel free to modify this but be careful.
    !          You may want to keep at least the initch4 as is to allow CH4 to be run even if it wasn't spun up
    !          with CH4, as the CH4 sub-model comes to eq. within a month plus one year for annual mean variables.
    !          Of course if clm_varctl:anoxia is used or NITRIF_DENITRIF is defined, then CN would no longer be in eq.

    end if
    call t_stopf('init_io2')

    ! ------------------------------------------------------------------------
    ! Initialize nitrogen deposition
    ! ------------------------------------------------------------------------

    if (use_cn) then
       call t_startf('init_ndep')
       call ndep_init(bounds)
       call ndep_interp(bounds)
       call t_stopf('init_ndep')
    end if
    
    ! ------------------------------------------------------------------------
    ! Initialization of model parameterizations that are needed after
    ! restart file is read in
    ! ------------------------------------------------------------------------

    ! Initialize active history fields. This is only done if not a restart run. 
    ! If a restart run, then this information has already been obtained from the 
    ! restart data read above. Note that routine hist_htapes_build needs time manager 
    ! information, so this call must be made after the restart information has been read.

    if (nsrest == nsrStartup .or. nsrest == nsrBranch) call hist_htapes_build()

    ! Initialize clmtype variables that are obtained from accumulated fields.
    ! This routine is called in an initial run at nstep=0
    ! This routine is also always called for a restart run and must 
    ! therefore be called after the restart file is read in

    call initAccClmtype(bounds)

    call t_stopf('init_hist1')

    ! --------------------------------------------------------------
    ! Note - everything below this point needs updated weights
    ! --------------------------------------------------------------

    ! Initialize filters
    
    call t_startf('init_filters')

    call allocFilters()

    !$OMP PARALLEL DO PRIVATE (nc)
    do nc = 1, nclumps
       call reweightWrapup(nc, bounds)
    end do
    !$OMP END PARALLEL DO

    call t_stopf('init_filters')

    ! Calculate urban "town" roughness length and displacement 
    ! height for urban landunits

    call UrbanInitAero(bounds)

    ! Initialize urban radiation model - this uses urbinp data structure

    call UrbanParamInit(bounds)

    ! Finalize urban model initialization
    
    call UrbanInput(bounds%begg, bounds%endg, mode='finalize')

    ! Even if CN is on, and dry-deposition is active, read CLMSP annual vegetation to get estimates of monthly LAI

    if ( n_drydep > 0 .and. drydep_method == DD_XLND )then
       call readAnnualVegetation(bounds)
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
          call interpMonthlyVeg(bounds)
       end if

       ! Determine gridcell averaged properties to send to atm

       call t_startf('init_map2gc')
       call clm_map2gcell(bounds, init=.true.)
       call t_stopf('init_map2gc')

    end if

    ! Initialize sno export state
    if (create_glacier_mec_landunit) then
       call t_startf('init_create_s2x')
       call update_clm_s2x(bounds, init=.true.)
       call t_stopf('init_create_s2x')
    end if

    ! wt_lunit was allocated in initialize1, but needed to be kept around through
    ! initialize2 for some consistency checking; now it can be deallocated
    deallocate(wt_lunit)

    call t_stopf('clm_init2')

  end subroutine initialize2

  !-----------------------------------------------------------------------
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
    !-----------------------------------------------------------------------

    if ( masterproc )then
      write(iulog,*) trim(version)
      write(iulog,*)
    end if

  end subroutine header

  !-----------------------------------------------------------------------
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
    !-----------------------------------------------------------------------

    do_restread = .false.
    if (nsrest == nsrStartup .and. finidat /= ' ') then
       do_restread = .true.
    end if
    if (nsrest == nsrContinue .or. nsrest == nsrBranch) then
       do_restread = .true.
    end if
  end function do_restread
  
end module clm_initializeMod
