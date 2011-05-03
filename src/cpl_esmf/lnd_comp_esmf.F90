module lnd_comp_esmf
  
!---------------------------------------------------------------------------
!BOP
!
! !MODULE: lnd_comp_esmf
!
!  Interface of the active land model component of CESM the CLM (Community Land Model)
!  with the main CESM driver. This is a thin interface taking CESM driver information
!  in MCT (Model Coupling Toolkit) format and converting it to use by CLM and outputing
!  if in ESMF (Earth System Modelling Framework) format.
!
! !DESCRIPTION:
!
! !USES:
  use shr_kind_mod , only : r8 => shr_kind_r8, SHR_KIND_CL
  use abortutils   , only : endrun
  use esmf_mod
  use esmfshr_mod
!
! !PUBLIC MEMBER FUNCTIONS:
  implicit none
  SAVE
  private                              ! By default make data private
!
! !PUBLIC MEMBER FUNCTIONS:
!
  public :: lnd_register_esmf          ! register clm initial, run, final methods
  public :: lnd_init_esmf              ! clm initialization
  public :: lnd_run_esmf               ! clm run phase
  public :: lnd_final_esmf             ! clm finalization/cleanup
!
! ! PUBLIC DATA MEMBERS: None
!
! !REVISION HISTORY:
!             Author:                                      Mariana Vertenstein, Fei Liu
! Apr/27/2010 Update to lnd_comp_mct version and add documentation Erik Kluzek
!
! !PRIVATE MEMBER FUNCTIONS:
  private :: lnd_DistGrid_esmf        ! Distribute clm grid
  private :: lnd_domain_esmf          ! Set the land model domain information
  private :: lnd_export_esmf          ! export land data to CESM coupler
  private :: lnd_import_esmf          ! import data from the CESM coupler to the land model
#ifdef RTM
  private :: rof_DistGrid_esmf        ! Distribute river runoff grid
  private :: rof_domain_esmf          ! Set the river runoff model domain information
  private :: rof_export_esmf          ! Export the river runoff model data to the CESM coupler
#endif
  private :: sno_export_esmf
  private :: sno_import_esmf
!
! !PRIVATE DATA MEMBERS:
!
! Time averaged flux fields
!  
  type(ESMF_Array)  :: l2x_SNAP       ! Snapshot of land to coupler data on the land grid
  type(ESMF_Array)  :: l2x_SUM        ! Summation of land to coupler data on the land grid

  type(ESMF_Array)  :: s2x_s_SNAP     ! Snapshot of sno to coupler data on the land grid
  type(ESMF_Array)  :: s2x_s_SUM      ! Summation of sno to coupler data on the land grid
!
! Time averaged counter for flux fields
!
  integer :: avg_count                ! Number of times snapshots of above flux data summed together
  integer :: avg_count_sno
!
! Atmospheric mode  
!
  logical :: atm_prognostic           ! Flag if active atmosphere component or not
!
! Internal clm grid
!
  real(R8), pointer :: fptr_l2x_clm(:,:) ! land to cpl export state
  real(R8), pointer :: fptr_x2l_clm(:,:) ! cpl to land import state

!EOP
!===============================================================
contains
!===============================================================

!---------------------------------------------------------------------------
!BOP
!
! !IROUTINE: lnd_register_esmf
!
! !INTERFACE:
subroutine lnd_register_esmf(comp, rc)
! !DESCRIPTION:
! Register the clm initial, run, and final phase methods with ESMF.
!
! !USES:
    implicit none
!
! !ARGUMENTS:
    type(ESMF_GridComp)  :: comp  ! CLM grid component
    integer, intent(out) :: rc    ! return status
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
!-----------------------------------------------------------------------
    rc = ESMF_SUCCESS
    ! Register the callback routines.

    call ESMF_GridCompSetEntryPoint(comp, ESMF_SETINIT, &
      lnd_init_esmf, phase=ESMF_SINGLEPHASE, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
    call ESMF_GridCompSetEntryPoint(comp, ESMF_SETRUN, &
      lnd_run_esmf, phase=ESMF_SINGLEPHASE, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
    call ESMF_GridCompSetEntryPoint(comp, ESMF_SETFINAL, &
      lnd_final_esmf, phase=ESMF_SINGLEPHASE, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)

end subroutine

!---------------------------------------------------------------------------
!BOP
!
! !IROUTINE: lnd_init_esmf
!
! !INTERFACE:
  subroutine lnd_init_esmf(comp, import_state, export_state, EClock, rc)
!
! !DESCRIPTION:
! Initialize land surface model and obtain relevant atmospheric model arrays
! back from (i.e. albedos, surface temperature and snow cover over land).
!
! !USES:
    use abortutils       , only : endrun
    use downscaleMod     , only : map1dl_a2l, map1dl_l2a, map_maparrayl
    use clm_time_manager , only : get_nstep, get_step_size, set_timemgr_init, &
                                  set_nextsw_cday
    use clm_atmlnd       , only : clm_l2a
    use clm_initializeMod, only : initialize1, initialize2
    use clm_varctl       , only : finidat,single_column, set_clmvarctl, iulog, noland, &
                                  downscale
    use controlMod       , only : control_setNL
    use decompMod        , only : get_proc_bounds, get_proc_bounds_atm
    use domainMod        , only : adomain
    use clm_varpar       , only : rtmlon, rtmlat
    use clm_varorb       , only : eccen, obliqr, lambm0, mvelpp
    use abortutils       , only : endrun
    use esmf_mod
    use clm_varctl       , only : iulog, noland
    use shr_file_mod     , only : shr_file_setLogUnit, shr_file_setLogLevel, &
                                  shr_file_getLogUnit, shr_file_getLogLevel, &
                                  shr_file_getUnit, shr_file_setIO
    use spmdMod          , only : masterproc, spmd_init
    use seq_timemgr_mod  , only : seq_timemgr_EClockGetData
    use seq_infodata_mod , only : seq_infodata_start_type_cont, &
                                  seq_infodata_start_type_brnch, &
	                          seq_infodata_start_type_start
    use seq_flds_mod
    use clm_cpl_indices  , only : clm_cpl_indices_set, nflds_l2x, nflds_x2l
    use clm_glclnd       , only : clm_maps2x, clm_s2x, atm_s2x, create_clm_s2x
    use clm_varctl       , only : create_glacier_mec_landunit, nsrStartup, &
                                  nsrContinue, nsrBranch
    implicit none
!
! !ARGUMENTS:
   type(ESMF_GridComp)          :: comp            ! CLM gridded component
   type(ESMF_State)             :: import_state    ! CLM import state
   type(ESMF_State)             :: export_state    ! CLM export state
   type(ESMF_Clock)             :: EClock          ! ESMF synchronization clock
   integer, intent(out)         :: rc              ! Return code
   
!
! !LOCAL VARIABLES:
    integer              :: mpicom_lnd, mpicom_vm, gsize

    type(ESMF_DistGrid)  :: distgrid_l, distgrid_r, distgrid_s
    type(ESMF_Array)     :: dom_l, dom_r, dom_s, l2x, x2l, r2x, s2x, x2s
    type(ESMF_VM)        :: vm
    integer, allocatable :: gindex(:)
    integer  :: lsize                                ! size of attribute vector
    integer  :: g,i,j                                ! indices
    integer  :: dtime_sync                           ! coupling time-step from the input synchronization clock
    integer  :: dtime_clm                            ! clm time-step
    logical  :: exists                               ! true if file exists
    real(r8) :: scmlat                               ! single-column latitude
    real(r8) :: scmlon                               ! single-column longitude
    real(r8) :: nextsw_cday                          ! calday from clock of next radiation computation
    character(len=SHR_KIND_CL) :: caseid             ! case identifier name
    character(len=SHR_KIND_CL) :: ctitle             ! case description title
    character(len=SHR_KIND_CL) :: starttype          ! start-type (startup, continue, branch, hybrid)
    character(len=SHR_KIND_CL) :: calendar           ! calendar type name
    character(len=SHR_KIND_CL) :: hostname           ! hostname of machine running on
    character(len=SHR_KIND_CL) :: version            ! Model version
    character(len=SHR_KIND_CL) :: username           ! user running the model
    integer :: nsrest                                ! clm restart type
    integer :: perpetual_ymd                         ! perpetual date
    integer :: ref_ymd                               ! reference date (YYYYMMDD)
    integer :: ref_tod                               ! reference time of day (sec)
    integer :: start_ymd                             ! start date (YYYYMMDD)
    integer :: start_tod                             ! start time of day (sec)
    integer :: stop_ymd                              ! stop date (YYYYMMDD)
    integer :: stop_tod                              ! stop time of day (sec)
    logical :: brnch_retain_casename                 ! flag if should retain the case name on a branch start type
    logical :: perpetual_run                         ! flag if should cycle over a perpetual date or not
    logical :: samegrid_al                           ! true if atmosphere and land are on the same grid
    logical :: atm_aero                              ! Flag if aerosol data sent from atm model
    integer :: lbnum                                 ! input to memory diagnostic
    integer :: shrlogunit,shrloglev                  ! old values for log unit and log level
    integer :: begg_l, endg_l, begg_a, endg_a
    real(R8), pointer :: fptr(:, :)
    character(len=32), parameter :: sub = 'lnd_init_esmf'
    character(len=*),  parameter :: format = "('("//trim(sub)//") :',A)"
    character(ESMF_MAXSTR) :: convCIM, purpComp
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein, Fei Liu
!
!EOP
!-----------------------------------------------------------------------

    ! Determine indices

    call clm_cpl_indices_set()

    rc = ESMF_SUCCESS
 
    ! duplicate the mpi communicator from the current VM 
    call ESMF_VMGetCurrent(vm, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
    call ESMF_VMGet(vm, mpiCommunicator=mpicom_vm, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
    call MPI_Comm_dup(mpicom_vm, mpicom_lnd, rc)
    if(rc /= 0) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)

    ! Initialize clm MPI communicator

    call spmd_init( mpicom_lnd )

#if (defined _MEMTRACE)
    if(masterproc) then
       lbnum=1
       call memmon_dump_fort('memmon.out','lnd_init_esmf:start::',lbnum)
    endif
#endif                      

    ! Initialize io log unit

    call shr_file_getLogUnit (shrlogunit)
    if (masterproc) then
       inquire(file='lnd_modelio.nml',exist=exists)
       if (exists) then
          iulog = shr_file_getUnit()
          call shr_file_setIO('lnd_modelio.nml',iulog)
       end if
       write(iulog,format) "CLM land model initialization"
    else
       iulog = shrlogunit
    end if

    call shr_file_getLogLevel(shrloglev)
    call shr_file_setLogUnit (iulog)
    
    ! Use infodata to set orbital values

    call ESMF_AttributeGet(export_state, name="orb_eccen", value=eccen, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
    call ESMF_AttributeGet(export_state, name="orb_mvelpp", value=mvelpp, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
    call ESMF_AttributeGet(export_state, name="orb_lambm0", value=lambm0, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
    call ESMF_AttributeGet(export_state, name="orb_obliqr", value=obliqr, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)

    ! Consistency check on namelist filename	

    call control_setNL( 'lnd_in' )

    ! Initialize clm
    ! initialize1 reads namelist, grid and surface data
    ! initialize2 performs rest of initialization    

    call seq_timemgr_EClockGetData(EClock,                               &
                                   start_ymd=start_ymd,                  &
                                   start_tod=start_tod, ref_ymd=ref_ymd, &
                                   ref_tod=ref_tod, stop_ymd=stop_ymd,   &
                                   stop_tod=stop_tod,                    &
                                   calendar=calendar )

    call ESMF_AttributeGet(export_state, name="perpetual", value=perpetual_run, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
    call ESMF_AttributeGet(export_state, name="perpetual_ymd", value=perpetual_ymd, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
    call ESMF_AttributeGet(export_state, name="case_name", value=caseid, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
    call ESMF_AttributeGet(export_state, name="case_desc", value=ctitle, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
    call ESMF_AttributeGet(export_state, name="single_column", value=single_column, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
    call ESMF_AttributeGet(export_state, name="scmlat", value=scmlat, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
    call ESMF_AttributeGet(export_state, name="scmlon", value=scmlon, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
    call ESMF_AttributeGet(export_state, name="brnch_retain_casename", value=brnch_retain_casename, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
    call ESMF_AttributeGet(export_state, name="start_type", value=starttype, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
    call ESMF_AttributeGet(export_state, name="model_version", value=version, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
    call ESMF_AttributeGet(export_state, name="hostname", value=hostname, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
    call ESMF_AttributeGet(export_state, name="username", value=username, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
    call ESMF_AttributeGet(export_state, name="samegrid_al", value=samegrid_al, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
    
    call set_timemgr_init( calendar_in=calendar, start_ymd_in=start_ymd, start_tod_in=start_tod, &
                           ref_ymd_in=ref_ymd, ref_tod_in=ref_tod, stop_ymd_in=stop_ymd,         &
                           stop_tod_in=stop_tod,  perpetual_run_in=perpetual_run,                &
                           perpetual_ymd_in=perpetual_ymd )

    if (     trim(starttype) == trim(seq_infodata_start_type_start)) then
       nsrest = nsrStartup
    else if (trim(starttype) == trim(seq_infodata_start_type_cont) ) then
       nsrest = nsrContinue
    else if (trim(starttype) == trim(seq_infodata_start_type_brnch)) then
       nsrest = nsrBranch
    else
       call endrun( sub//' ERROR: unknown starttype' )
    end if

    call set_clmvarctl(    caseid_in=caseid, ctitle_in=ctitle,                     &
                           brnch_retain_casename_in=brnch_retain_casename,         &
                           single_column_in=single_column, scmlat_in=scmlat,       &
                           scmlon_in=scmlon, nsrest_in=nsrest, version_in=version, &
                           hostname_in=hostname, username_in=username )

    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)

    call initialize1( )

    ! If no land then exit out of initialization

    if ( noland) then
       call ESMF_AttributeSet(export_state, name="lnd_present", value=.false., rc=rc)
       if( rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
       call ESMF_AttributeSet(export_state, name="lnd_prognostic", value=.false., rc=rc)
       if( rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
       call ESMF_AttributeSet(export_state, name="rof_present", value=.false., rc=rc)
       if( rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
    end if

    ! If trigrid AND downscale -- abort as an error

    if ( (.not. samegrid_al) .and. downscale )then
       write(iulog,format) "Incompatible inputs, atmosphere and land are on different grids"
       write(iulog,format) "But, you are also trying to use the CLM fine-mesh to downscale to a finer grid"
       write(iulog,format) "When using CLM finemesh mode -- you must NOT run the atmosphere on a different grid"
       call endrun( sub//' ERROR: downscaling using finemesh AND running atmosphere on a different grid' )
    end if

    ! Determine if aerosol and dust deposition come from atmosphere component

    rc = ESMF_SUCCESS

    call ESMF_AttributeGet(export_state, name="atm_aero", value=atm_aero, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)

    if ( .not. atm_aero )then
       call endrun( sub//' ERROR: atmosphere model MUST send aerosols to CLM' )
    end if

    ! Initialize lnd distributed grid

    distgrid_l = lnd_DistGrid_esmf(gsize, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)

    call ESMF_AttributeSet(export_state, name="gsize_lnd", value=gsize, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)

    dom_l = mct2esmf_init(distgrid_l, attname=seq_flds_dom_fields, name="domain_l", rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)

    call lnd_domain_esmf( dom_l, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)

    ! Finish initializing clm

    call initialize2()

#ifdef RTM
    ! Initialize rof distgrid and domain

    distgrid_r = rof_DistGrid_esmf(gsize, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
    call ESMF_AttributeSet(export_state, name="gsize_rof", value=gsize, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)

    dom_r = mct2esmf_init(distgrid_r, attname=seq_flds_dom_fields, name="domain_r", rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)

    call rof_domain_esmf( dom_r, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
#endif

    ! Initialize sno distgrid and domain (currently same as lnd - ask Jon)

    distgrid_s = lnd_DistGrid_esmf(gsize, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
    call ESMF_AttributeSet(export_state, name="gsize_sno", value=gsize, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)

    dom_s = mct2esmf_init(distgrid_s, attname=seq_flds_dom_fields, name="domain_s", rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)

    call lnd_domain_esmf( dom_s, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
 
    ! Initialize lnd import and export states

    l2x = mct2esmf_init(distgrid_l, attname=seq_flds_l2x_fields, name="l2x", rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
    r2x = mct2esmf_init(distgrid_r, attname=seq_flds_r2x_fields, name="r2x", rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
    s2x = mct2esmf_init(distgrid_s, attname=seq_flds_s2x_fields, name="s2x", rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)

    x2l = mct2esmf_init(distgrid_l, attname=seq_flds_x2l_fields, name="x2l", rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
    x2s = mct2esmf_init(distgrid_s, attname=seq_flds_x2s_fields, name="x2s", rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)

    call ESMF_StateAdd(export_state, dom_l, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
    call ESMF_StateAdd(export_state, dom_r, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
    call ESMF_StateAdd(export_state, dom_s, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)

    call ESMF_StateAdd(export_state, l2x, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
    call ESMF_StateAdd(export_state, r2x, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
    call ESMF_StateAdd(export_state, s2x, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)

    call ESMF_StateAdd(import_state, x2l, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
    call ESMF_StateAdd(import_state, x2s, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)

    l2x_SNAP = mct2esmf_init(distgrid_l, attname=seq_flds_l2x_fluxes, &
        name="l2x_SNAP", rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
    l2x_SUM = mct2esmf_init(distgrid_l, attname=seq_flds_l2x_fluxes, &
        name="l2x_SUM", rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)

    if (masterproc) then
       write(iulog,format)'time averaging the following flux fields over the coupling interval'
       write(iulog,format) trim(seq_flds_l2x_fluxes)
    end if

    ! Check that clm internal dtime aligns with clm coupling interval

    call seq_timemgr_EClockGetData(EClock, dtime=dtime_sync )
    dtime_clm = get_step_size()
    if(masterproc) write(iulog,*)'dtime_sync= ',dtime_sync,' dtime_clm= ',dtime_clm,' mod = ',mod(dtime_sync,dtime_clm)
    if (mod(dtime_sync,dtime_clm) /= 0) then
       write(iulog,*)'clm dtime ',dtime_clm,' and Eclock dtime ',dtime_sync,' never align'
       call endrun( sub//' ERROR: time out of sync' )
    end if

    ! Create land export state then map this from the 
    ! clm internal grid to the clm driver (atm) grid 

    call get_proc_bounds_atm(begg_a, endg_a) ! clm grid - driver 
    call get_proc_bounds    (begg_l, endg_l) ! clm grid - internal (atm)

    allocate(fptr_x2l_clm(nflds_x2l, endg_l-begg_l+1))
    allocate(fptr_l2x_clm(nflds_l2x, endg_l-begg_l+1))

    call ESMF_ArrayGet(l2x, localDe=0, farrayPtr=fptr, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)

#ifdef RTM
    ! Initialize rof export state

    call rof_export_esmf( r2x , rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
#endif

    if (create_glacier_mec_landunit) then

       call sno_export_esmf( atm_s2x, s2x, rc=rc)
       if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)

    endif   ! create_glacier_mec_landunit

    ! Initialize averaging counter

    avg_count = 0

    ! Set land modes

    call ESMF_AttributeSet(export_state, name="lnd_prognostic", value=.true., rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
    call ESMF_AttributeSet(export_state, name="lnd_nx", value=adomain%ni, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
    call ESMF_AttributeSet(export_state, name="lnd_ny", value=adomain%nj, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)

#ifdef RTM
    call ESMF_AttributeSet(export_state, name="rof_present", value=.true., rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
    call ESMF_AttributeSet(export_state, name="rof_nx", value=rtmlon, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
    call ESMF_AttributeSet(export_state, name="rof_ny", value=rtmlat, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
#else
    call ESMF_AttributeSet(export_state, name="rof_present", value=.false., rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
#endif

    if (create_glacier_mec_landunit) then
       call ESMF_AttributeSet(export_state, name="sno_present", value=.true., rc=rc)
       if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
       call ESMF_AttributeSet(export_state, name="sno_prognostic", value=.true., rc=rc)
       if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
       call ESMF_AttributeSet(export_state, name="sno_nx", value=adomain%ni, rc=rc)
       if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
       call ESMF_AttributeSet(export_state, name="sno_ny", value=adomain%nj, rc=rc)
       if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
    else
       call ESMF_AttributeSet(export_state, name="sno_present", value=.false., rc=rc)
       if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
    endif

    call ESMF_AttributeGet(export_state, name="nextsw_cday", value=nextsw_cday, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
    call set_nextsw_cday( nextsw_cday )

    ! Create land export state 
    ! If downscale, map this from the clm internal grid (dst) to the clm driver grid (src)
    ! Reset landfrac on atmosphere grid to have the right domain

    if (downscale) then
       call lnd_export_esmf(clm_l2a, fptr_l2x_clm, begg_l, endg_l)
       call map_maparrayl(begg_l, endg_l, begg_a, endg_a, nflds_l2x, &
                          fptr_l2x_clm, fptr, map1dl_l2a, reverse_order=.true.)
    else
       call lnd_export_esmf(clm_l2a, fptr, begg_l, endg_l)
    end if
    
#ifdef RTM
    call rof_export_esmf( r2x , rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
#endif

    !call sno_export_esmf( s2x, rc=rc)
    !if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)

    ! Determine atmosphere modes

    call ESMF_AttributeGet(export_state, name="atm_prognostic", value=atm_prognostic, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)

    if (masterproc) then
       if ( atm_prognostic )then
          write(iulog,format) 'Atmospheric input is from a prognostic model'
       else
          write(iulog,format) 'Atmospheric input is from a data model'
       end if
    end if

    ! Reset shr logging to original values

    call shr_file_setLogUnit (shrlogunit)
    call shr_file_setLogLevel(shrloglev)

#if (defined _MEMTRACE)
    if(masterproc) then
       write(iulog,*) TRIM(Sub) // ':end::'
       lbnum=1
       call memmon_dump_fort('memmon.out','lnd_int_esmf:end::',lbnum)
       call memmon_reset_addr()
    endif
#endif

    convCIM  = "CIM 1.0"
    purpComp = "Model Component Simulation Description"

    call ESMF_AttributeAdd(comp,  &
                           convention=convCIM, purpose=purpComp, rc=rc)

    call ESMF_AttributeSet(comp, "ShortName", "CLM", &
                           convention=convCIM, purpose=purpComp, rc=rc)
    call ESMF_AttributeSet(comp, "LongName", &
                           "Community Land Model", &
                           convention=convCIM, purpose=purpComp, rc=rc)
    call ESMF_AttributeSet(comp, "Description", &
                  "The Community Land Model version 4.0 is " // &
                  "the land model used in the CESM1.0.  " // &
                  "More information on the CLM project " // &
                  "and access to previous CLM model versions and " // &
                  "documentation can be found via the CLM Web Page.", &
                           convention=convCIM, purpose=purpComp, rc=rc)
    call ESMF_AttributeSet(comp, "ReleaseDate", "2010", &
                           convention=convCIM, purpose=purpComp, rc=rc)
    call ESMF_AttributeSet(comp, "ModelType", "Land", &
                           convention=convCIM, purpose=purpComp, rc=rc)

    call ESMF_AttributeSet(comp, "Name", "Sam Levis", &
                           convention=convCIM, purpose=purpComp, rc=rc)
    call ESMF_AttributeSet(comp, "EmailAddress", &
                           "slevis@ucar.edu", &
                           convention=convCIM, purpose=purpComp, rc=rc)
    call ESMF_AttributeSet(comp, "ResponsiblePartyRole", "contact", &
                           convention=convCIM, purpose=purpComp, rc=rc)

  end subroutine lnd_init_esmf

!---------------------------------------------------------------------------
!BOP
!
! !IROUTINE: lnd_run_esmf
!
! !INTERFACE:
subroutine lnd_run_esmf(comp, import_state, export_state, EClock, rc)
!
! !DESCRIPTION:
! Run clm model
!
! !USES:
    use shr_kind_mod    ,only : r8 => shr_kind_r8
    use downscaleMod    ,only : map1dl_a2l, map1dl_l2a, map_maparrayl
    use clm_atmlnd      ,only : clm_l2a, atm_a2l, clm_a2l, clm_downscale_a2l
    use clm_driver      ,only : clm_drv
    use clm_varorb      ,only : eccen, obliqr, lambm0, mvelpp
    use clm_time_manager,only : get_curr_date, get_nstep, get_curr_calday, get_step_size, &
                                advance_timestep, set_nextsw_cday,update_rad_dtime
    use domainMod       ,only : adomain, ldomain
    use decompMod       ,only : get_proc_bounds_atm, get_proc_bounds
    use abortutils      ,only : endrun
    use esmf_mod
    use clm_varctl      ,only : iulog
    use shr_file_mod    ,only : shr_file_setLogUnit, shr_file_setLogLevel, &
                                shr_file_getLogUnit, shr_file_getLogLevel
    use seq_timemgr_mod ,only : seq_timemgr_EClockGetData, seq_timemgr_StopAlarmIsOn, &
                                seq_timemgr_RestartAlarmIsOn, seq_timemgr_EClockDateInSync
    use spmdMod         ,only : masterproc, mpicom
    use perf_mod        ,only : t_startf, t_stopf, t_barrierf
    use clm_varctl      ,only : create_glacier_mec_landunit, downscale
    use clm_glclnd      ,only : clm_maps2x, clm_mapx2s, clm_s2x, atm_s2x, atm_x2s, clm_x2s
    use clm_glclnd      ,only : create_clm_s2x, unpack_clm_x2s
    use shr_orb_mod     ,only : shr_orb_decl
    use clm_varorb      ,only : eccen, mvelpp, lambm0, obliqr
    use clm_cpl_indices ,only : nflds_l2x, nflds_x2l
    implicit none
!
! !ARGUMENTS:
   type(ESMF_GridComp)          :: comp            ! CLM gridded component
   type(ESMF_State)             :: import_state    ! CLM import state
   type(ESMF_State)             :: export_state    ! CLM export state
   type(ESMF_Clock)             :: EClock          ! ESMF synchronization clock
   integer, intent(out)         :: rc              ! Return code
!
! !LOCAL VARIABLES:
    type(ESMF_Array)  :: l2x, r2x, s2x, x2l, x2s, dom_l
    real(R8), pointer :: fptr(:, :)
    integer :: ymd_sync                   ! Sync date (YYYYMMDD)
    integer :: yr_sync                    ! Sync current year
    integer :: mon_sync                   ! Sync current month
    integer :: day_sync                   ! Sync current day
    integer :: tod_sync                   ! Sync current time of day (sec)
    integer :: ymd                        ! CLM current date (YYYYMMDD)
    integer :: yr                         ! CLM current year
    integer :: mon                        ! CLM current month
    integer :: day                        ! CLM current day
    integer :: tod                        ! CLM current time of day (sec)
    integer :: dtime                      ! time step increment (sec)
    integer :: nstep                      ! time step index
    logical :: rstwr_sync                 ! .true. ==> write restart file before returning
    logical :: rstwr                      ! .true. ==> write restart file before returning
    logical :: nlend_sync                 ! Flag signaling last time-step
    logical :: nlend                      ! .true. ==> last time-step
    logical :: dosend                     ! true => send data back to driver
    logical :: doalb                      ! .true. ==> do albedo calculation on this time step
    real(r8):: nextsw_cday                ! calday from clock of next radiation computation
    real(r8):: caldayp1                   ! clm calday plus dtime offset
    real(r8):: calday                     ! calendar day for nstep
    real(r8):: declin                     ! solar declination angle in radians for nstep
    real(r8):: declinp1                   ! solar declination angle in radians for nstep+1
    real(r8):: eccf                       ! earth orbit eccentricity factor
    integer :: shrlogunit,shrloglev       ! old values
    integer :: begg_l, endg_l, begg_a, endg_a ! Beginning and ending gridcell index numbers
    integer :: lbnum                          ! input to memory diagnostic
    integer :: g,i,ka                         ! counters
    logical,save :: first_call = .true.   ! first call work
    character(len=32)            :: rdate ! date char string for restart file names
    character(len=32), parameter :: sub = "lnd_run_esmf"
    logical :: glcrun_alarm    ! if true, sno data is averaged and sent to glc this step
    logical :: update_glc2sno_fields  ! if true, update glacier_mec fields
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein, Fei Liu
!
!EOP
!---------------------------------------------------------------------------

    call t_startf ('lc_lnd_run1')
    rc = ESMF_SUCCESS

#if (defined _MEMTRACE)
    if(masterproc) then
       lbnum=1
       call memmon_dump_fort('memmon.out','lnd_run_esmf:start::',lbnum)
    endif
#endif

    ! Reset shr logging to my log file

    call shr_file_getLogUnit (shrlogunit)
    call shr_file_getLogLevel(shrloglev)
    call shr_file_setLogUnit (iulog)

    ! Use infodata to set orbital values if it was updated at run time

    call ESMF_AttributeGet(export_state, name="orb_eccen", value=eccen, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
    call ESMF_AttributeGet(export_state, name="orb_mvelpp", value=mvelpp, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
    call ESMF_AttributeGet(export_state, name="orb_lambm0", value=lambm0, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
    call ESMF_AttributeGet(export_state, name="orb_obliqr", value=obliqr, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)

    ! Determine time of next atmospheric shortwave calculation

    call seq_timemgr_EClockGetData(EClock, &
         curr_ymd=ymd, curr_tod=tod_sync,  &
         curr_yr=yr_sync, curr_mon=mon_sync, curr_day=day_sync)
    call ESMF_AttributeGet(export_state, name="nextsw_cday", value=nextsw_cday, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)

    call set_nextsw_cday( nextsw_cday )
    dtime = get_step_size()

    write(rdate,'(i4.4,"-",i2.2,"-",i2.2,"-",i5.5)') yr_sync,mon_sync,day_sync,tod_sync
    nlend_sync = seq_timemgr_StopAlarmIsOn( EClock )
    rstwr_sync = seq_timemgr_RestartAlarmIsOn( EClock )

    call get_proc_bounds_atm(begg_a, endg_a)
    call get_proc_bounds    (begg_l, endg_l)

    if (first_call) then
       call ESMF_StateGet(export_state, itemName="domain_l", array=dom_l, rc=rc)
       if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
       call ESMF_ArrayGet(dom_l, localDe=0, farrayPtr=fptr, rc=rc)
       if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
       ka = esmfshr_util_ArrayGetIndex(dom_l,'ascale',rc=rc)
       if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
       do g = begg_a,endg_a
          i = 1 + (g - begg_a)
          adomain%asca(g) = fptr(ka, i)
          if ( .not. downscale )then
             ldomain%asca(g) = adomain%asca(g)
          end if
       end do
    endif
    call t_stopf ('lc_lnd_run1')
    
    ! Map ESMF to CLM data type
    ! Perform downscaling if appropriate

    call t_startf ('lc_lnd_import')
    call ESMF_StateGet(import_state, itemName="x2l", array=x2l, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
    call ESMF_ArrayGet(x2l, localDe=0, farrayPtr=fptr, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
    if (downscale) then
       call lnd_import_esmf( fptr, atm_a2l, begg_a, endg_a )
       call map_maparrayl(begg_a, endg_a, begg_l, endg_l, nflds_x2l, &
	                  fptr, fptr_x2l_clm, map1dl_a2l, reverse_order=.true.)
       call lnd_import_esmf( fptr_x2l_clm, clm_a2l, begg_l, endg_l )
       call clm_downscale_a2l(atm_a2l, clm_a2l)
    else
       call lnd_import_esmf( fptr, clm_a2l, begg_l, endg_l )
    end if
    
    ! Map to clm (only when state and/or fluxes need to be updated)

    if (create_glacier_mec_landunit) then
       update_glc2sno_fields  = .false.
       call ESMF_AttributeGet(export_state, name="glc_g2supdate", value=update_glc2sno_fields, rc=rc)
       if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)

       if (update_glc2sno_fields) then
          if (downscale) then
             call sno_import_esmf( x2s, atm_x2s, rc=rc )
             if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
             call clm_mapx2s(atm_x2s, clm_x2s)
          else
             call sno_import_esmf( x2s, clm_x2s, rc=rc )
             if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
          end if
          call unpack_clm_x2s(clm_x2s)
       endif   ! update_glc2sno
    endif   ! create_glacier_mec_landunit
    call t_stopf ('lc_lnd_import')

    ! Loop over time steps in coupling interval

    call t_startf ('lc_lnd_run2')

    dosend      = .false.
    do while(.not. dosend)

       ! Determine if dosend
       ! When time is not updated at the beginning of the loop - then return only if
       ! are in sync with clock before time is updated

       call get_curr_date( yr, mon, day, tod )
       ymd = yr*10000 + mon*100 + day
       tod = tod
       dosend = (seq_timemgr_EClockDateInSync( EClock, ymd, tod))

       ! Determine doalb based on nextsw_cday sent from atm model

       nstep = get_nstep()
       caldayp1 = get_curr_calday(offset=dtime)
       doalb = abs(nextsw_cday- caldayp1) < 1.e-10_r8
       if (nstep == 0) then
          doalb = .false. 	
       else if (nstep == 1) then 
          doalb = (abs(nextsw_cday- caldayp1) < 1.e-10_r8) 
       else
          doalb = (nextsw_cday >= -0.5_r8) 
       end if
       call update_rad_dtime(doalb)

       ! Determine if time to write cam restart and stop

       rstwr = .false.
       if (rstwr_sync .and. dosend) rstwr = .true.
       nlend = .false.
       if (nlend_sync .and. dosend) nlend = .true.

       ! Run clm 

       call t_barrierf('sync_clm_run', mpicom)
       call t_startf ('clm_run')
       calday = get_curr_calday()
       call shr_orb_decl( calday     , eccen, mvelpp, lambm0, obliqr, declin  , eccf )
       call shr_orb_decl( nextsw_cday, eccen, mvelpp, lambm0, obliqr, declinp1, eccf )
       call clm_drv(doalb, nextsw_cday, declinp1, declin, rstwr, nlend, rdate)
       call t_stopf ('clm_run')

       ! Map CLM data type to MCT
       ! Reset landfrac on atmosphere grid to have the right domain
       
       call t_startf ('lc_lnd_export')
       call ESMF_StateGet(export_state, itemName="l2x", array=l2x, rc=rc)
       if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
       call ESMF_ArrayGet(l2x, localDe=0, farrayPtr=fptr, rc=rc)
       if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)

       if (downscale) then
          call lnd_export_esmf( clm_l2a, fptr_l2x_clm, begg_l, endg_l )
          call map_maparrayl(begg_l, endg_l, begg_a, endg_a, nflds_l2x, &
	            	     fptr_l2x_clm, fptr, map1dl_l2a)
       else
          call lnd_export_esmf( clm_l2a, fptr, begg_l, endg_l )
       end if
       call t_stopf ('lc_lnd_export')
       
       ! Compute snapshot attribute vector for accumulation
       
       ! don't accumulate on first coupling freq ts0 and ts1
       ! for consistency with ccsm3 when flxave is off
 
       nstep = get_nstep()
       if (nstep <= 1) then
          call esmfshr_util_ArrayCopy(l2x, l2x_SUM, rc=rc)
          if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
          avg_count = 1
       else
          call esmfshr_util_ArrayCopy(l2x, l2x_SNAP, rc=rc)
          if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
          call esmfshr_util_ArraySum(l2x_SNAP, l2x_SUM, rc=rc)
          if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
          avg_count = avg_count + 1
       endif

       ! Map sno data type to MCT

       if (create_glacier_mec_landunit) then
          call create_clm_s2x(clm_s2x)
          call ESMF_StateGet(export_state, itemName="s2x", array=s2x, rc=rc)
          if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
          if (downscale) then
             call clm_maps2x(clm_s2x, atm_s2x)
             call sno_export_esmf( atm_s2x, s2x, rc=rc )
             if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
          else
             call sno_export_esmf( clm_s2x, s2x, rc=rc )
             if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
          end if

          if (nstep <= 1) then
             call esmfshr_util_ArrayCopy(s2x, s2x_s_SUM, rc=rc)
             if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
             avg_count_sno = 1
          else
             call esmfshr_util_ArrayCopy(s2x, s2x_s_SNAP, rc=rc)
             if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
             call esmfshr_util_ArraySum(s2x_s_SNAP, s2x_s_SUM, rc=rc)
             if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
             avg_count_sno = avg_count + 1
          endif
       endif    ! create_glacier_mec_landunit

       ! Advance clm time step
       
       call t_startf ('lc_clm2_adv_timestep')
       call advance_timestep()
       call t_stopf ('lc_clm2_adv_timestep')

    end do

    call t_stopf ('lc_lnd_run2')
    call t_startf ('lc_lnd_run3')

    ! Finish accumulation of attribute vector and average and zero out partial sum and counter
    call esmfshr_util_ArrayAvg(l2x_SUM, avg_count, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
    call esmfshr_util_ArrayCopy(l2x_SUM, l2x, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
    call esmfshr_util_ArrayZero(l2x_SUM, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
    avg_count = 0                   

    if (create_glacier_mec_landunit) then
       call ESMF_AttributeGet(export_state, name="glcrun_alarm", value=glcrun_alarm, rc=rc)
       if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
       if (glcrun_alarm) then
          call esmfshr_util_ArrayAvg(s2x_s_SUM, avg_count, rc=rc)
          if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
          call esmfshr_util_ArrayCopy(s2x_s_SUM, s2x, rc=rc)
          if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
          call esmfshr_util_ArrayZero(s2x_s_SUM, rc=rc)
          if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
          avg_count_sno = 0

       endif
    endif

#ifdef RTM
    ! Create river runoff output state

    call t_startf ('lc_rof_export')
    call ESMF_StateGet(export_state, itemName="r2x", array=r2x, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
    call rof_export_esmf(r2x, rc=rc )
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
    call t_stopf ('lc_rof_export')
#endif
       
    ! Check that internal clock is in sync with master clock

    call get_curr_date( yr, mon, day, tod, offset=-dtime )
    ymd = yr*10000 + mon*100 + day
    tod = tod
    if ( .not. seq_timemgr_EClockDateInSync( EClock, ymd, tod ) )then
       call seq_timemgr_EclockGetData( EClock, curr_ymd=ymd_sync, curr_tod=tod_sync )
       write(iulog,*)' clm ymd=',ymd     ,'  clm tod= ',tod
       write(iulog,*)'sync ymd=',ymd_sync,' sync tod= ',tod_sync
       call endrun( sub//":: CLM clock not in sync with Master Sync clock" )
    end if
    
    ! Reset shr logging to my original values

    call shr_file_setLogUnit (shrlogunit)
    call shr_file_setLogLevel(shrloglev)
  
#if (defined _MEMTRACE)
    if(masterproc) then
       lbnum=1
       call memmon_dump_fort('memmon.out','lnd_run_esmf:end::',lbnum)
       call memmon_reset_addr()
    endif
#endif

    first_call = .false.
    call t_stopf ('lc_lnd_run3')

  end subroutine lnd_run_esmf

!---------------------------------------------------------------------------
!BOP
!
! !IROUTINE: lnd_final_esmf
!
! !INTERFACE:
subroutine lnd_final_esmf(comp, import_state, export_state, EClock, rc)
!EOP
!
! !DESCRIPTION:
! Finalize land surface model
!
!------------------------------------------------------------------------------
!BOP
!
! !ARGUMENTS:
!
   implicit none
   type(ESMF_GridComp)          :: comp            ! CLM gridded component
   type(ESMF_State)             :: import_state    ! CLM import state
   type(ESMF_State)             :: export_state    ! CLM export state
   type(ESMF_Clock)             :: EClock          ! ESMF synchronization clock
   integer, intent(out)         :: rc              ! Return code

! !REVISION HISTORY:
! Author: Mariana Vertenstein, Fei Liu
!
!EOP
!---------------------------------------------------------------------------

   ! fill this in
    rc = ESMF_SUCCESS

    ! Destroy ESMF objects
    call esmfshr_util_StateArrayDestroy(export_state,'domain_l',rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
    call esmfshr_util_StateArrayDestroy(export_state,'domain_r',rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
    call esmfshr_util_StateArrayDestroy(export_state,'domain_s',rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
    call esmfshr_util_StateArrayDestroy(export_state,'l2x',rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
    call esmfshr_util_StateArrayDestroy(export_state,'r2x',rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
    call esmfshr_util_StateArrayDestroy(export_state,'s2x',rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)

    call esmfshr_util_StateArrayDestroy(import_state,'x2l',rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
    call esmfshr_util_StateArrayDestroy(import_state,'x2s',rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)

  end subroutine lnd_final_esmf

!=================================================================================

!---------------------------------------------------------------------------
!BOP
!
! !IROUTINE: lnd_DistGrid_esmf
!
! !INTERFACE:
  function lnd_DistGrid_esmf(gsize, rc)
!
! !DESCRIPTION:
! Setup distributed grid for CLM
!
! !USES:
    use shr_kind_mod , only : r8 => shr_kind_r8
    use decompMod    , only : get_proc_bounds_atm, adecomp
    use domainMod    , only : adomain
!
! !ARGUMENTS:
    integer, intent(out) :: gsize   ! grid size
    integer, intent(out) :: rc      ! return code
!
! RETURN:
    type(ESMF_DistGrid) :: lnd_DistGrid_esmf  ! Resulting distributed grid
!
! !LOCAL VARIABLES:
    integer,allocatable :: gindex(:)  ! grid indices
    integer :: i, j, n, gi            ! indices
    integer :: ier                    ! error code
    integer :: begg, endg             ! grid begining and ending indices
    integer :: lsize                  ! span of grid for this task
!
! !REVISION HISTORY:
! Author: Fei Liu
!
!EOP
!---------------------------------------------------------------------------
    ! Build the land grid numbering
    ! NOTE:  Numbering scheme is: West to East and South to North
    ! starting at south pole.  Should be the same as what's used in SCRIP
    
    rc = ESMF_SUCCESS

    call get_proc_bounds_atm(begg, endg)

    allocate(gindex(begg:endg),stat=ier)

    ! number the local grid

    do n = begg, endg
        gindex(n) = adecomp%gdc2glo(n)
    end do
    lsize = endg-begg+1
    gsize = adomain%ni*adomain%nj

    lnd_DistGrid_esmf = mct2esmf_init(gindex, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)

    deallocate(gindex)

  end function lnd_DistGrid_esmf

!=================================================================================

!---------------------------------------------------------------------------
!BOP
!
! !IROUTINE: lnd_export_esmf
!
! !INTERFACE:
  subroutine lnd_export_esmf( l2a, fptr, begg, endg )
!
! !DESCRIPTION:
!
! Convert the data to be sent from the clm model to the coupler from clm data types
! to ESMF data types.
!
!---------------------------------------------------------------------------
! !USES:
    use shr_kind_mod    , only : r8 => shr_kind_r8
    use clm_time_manager, only : get_nstep  
    use clm_atmlnd      , only : lnd2atm_type
    use domainMod       , only : adomain
    use seq_drydep_mod  , only : n_drydep
    use clm_cpl_indices
    implicit none
! !ARGUMENTS:
    type(lnd2atm_type), intent(inout) :: l2a         ! clm land to atmosphere exchange data type
    real(R8)          , pointer       :: fptr(:, :)
    integer           , intent(in)    :: begg, endg 
!
! !LOCAL VARIABLES:
    integer :: g,i           ! indices
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
!---------------------------------------------------------------------------

    ! cesm sign convention is that fluxes are positive downward

    fptr(:,:) = 0.0_r8

    do g = begg,endg
       i = 1 + (g-begg)
       fptr(index_l2x_Sl_t,i)        =  l2a%t_rad(g)
       fptr(index_l2x_Sl_snowh,i)    =  l2a%h2osno(g)
       fptr(index_l2x_Sl_avsdr,i)    =  l2a%albd(g,1)
       fptr(index_l2x_Sl_anidr,i)    =  l2a%albd(g,2)
       fptr(index_l2x_Sl_avsdf,i)    =  l2a%albi(g,1)
       fptr(index_l2x_Sl_anidf,i)    =  l2a%albi(g,2)
       fptr(index_l2x_Sl_tref,i)     =  l2a%t_ref2m(g)
       fptr(index_l2x_Sl_qref,i)     =  l2a%q_ref2m(g)
       fptr(index_l2x_Sl_u10,i)      =  l2a%u_ref10m(g)
       fptr(index_l2x_Fall_taux,i)   = -l2a%taux(g)
       fptr(index_l2x_Fall_tauy,i)   = -l2a%tauy(g)
       fptr(index_l2x_Fall_lat,i)    = -l2a%eflx_lh_tot(g)
       fptr(index_l2x_Fall_sen,i)    = -l2a%eflx_sh_tot(g)
       fptr(index_l2x_Fall_lwup,i)   = -l2a%eflx_lwrad_out(g)
       fptr(index_l2x_Fall_evap,i)   = -l2a%qflx_evap_tot(g)
       fptr(index_l2x_Fall_swnet,i)  =  l2a%fsa(g)
       if (index_l2x_Fall_fco2_lnd /= 0) then
          fptr(index_l2x_Fall_fco2_lnd,i) = -l2a%nee(g)  
       end if

       ! Additional fields for DUST, PROGSSLT, dry-deposition and VOC
       ! These are now standard fields, but the check on the index makes sure the driver handles them
       if (index_l2x_Sl_ram1      /= 0 )  fptr(index_l2x_Sl_ram1,i) = l2a%ram1(g)
       if (index_l2x_Sl_fv        /= 0 )  fptr(index_l2x_Sl_fv,i)   = l2a%fv(g)
       if (index_l2x_Fall_flxdst1 /= 0 )  fptr(index_l2x_Fall_flxdst1,i)= -l2a%flxdst(g,1)
       if (index_l2x_Fall_flxdst2 /= 0 )  fptr(index_l2x_Fall_flxdst2,i)= -l2a%flxdst(g,2)
       if (index_l2x_Fall_flxdst3 /= 0 )  fptr(index_l2x_Fall_flxdst3,i)= -l2a%flxdst(g,3)
       if (index_l2x_Fall_flxdst4 /= 0 )  fptr(index_l2x_Fall_flxdst4,i)= -l2a%flxdst(g,4)
       if (index_l2x_Sl_ddvel     /= 0 )  fptr(index_l2x_Sl_ddvel:index_l2x_Sl_ddvel+n_drydep-1,i) &
                                                                        = l2a%ddvel(g,:n_drydep)
       if (index_l2x_Fall_flxvoc1 /= 0 )  fptr(index_l2x_Fall_flxvoc1,i)= -l2a%flxvoc(g,1)
       if (index_l2x_Fall_flxvoc2 /= 0 )  fptr(index_l2x_Fall_flxvoc2,i)= -l2a%flxvoc(g,2)
       if (index_l2x_Fall_flxvoc3 /= 0 )  fptr(index_l2x_Fall_flxvoc3,i)= -l2a%flxvoc(g,3)
       if (index_l2x_Fall_flxvoc4 /= 0 )  fptr(index_l2x_Fall_flxvoc4,i)= -l2a%flxvoc(g,4)
       if (index_l2x_Fall_flxvoc5 /= 0 )  fptr(index_l2x_Fall_flxvoc5,i)= -l2a%flxvoc(g,5)
    end do

  end subroutine lnd_export_esmf

!====================================================================================

!---------------------------------------------------------------------------
!BOP
!
! !IROUTINE: lnd_import_esmf
!
! !INTERFACE:
  subroutine lnd_import_esmf( fptr, a2l, begg, endg)
!
! !DESCRIPTION:
!
! Convert the input data from the coupler to the land model from ESMF import state
! into internal clm data types.
!
!---------------------------------------------------------------------------
! !USES:
    use shr_kind_mod    , only: r8 => shr_kind_r8
    use clm_atmlnd      , only: atm2lnd_type
    use clm_varctl      , only: co2_type, co2_ppmv
    use clm_varcon      , only: rair, o2_molar_const
#if (defined C13)
    use clm_varcon      , only: c13ratio
#endif
    use shr_const_mod   , only: SHR_CONST_TKFRZ
    use decompMod       , only: get_proc_bounds_atm
    use clm_varctl      , only: iulog
    use clm_cpl_indices
    implicit none
! !ARGUMENTS:
    real(r8)          , pointer       :: fptr(:,:)
    type(atm2lnd_type), intent(inout) :: a2l          ! clm internal input data type
    integer           , intent(in)    :: begg, endg   ! beginning and ending gridcell indices
!
! !LOCAL VARIABLES:
    integer  :: g,i,nstep,ier        ! indices, number of steps, and error code
    real(r8) :: forc_rainc           ! rainxy Atm flux mm/s
    real(r8) :: e                    ! vapor pressure (Pa)
    real(r8) :: qsat                 ! saturation specific humidity (kg/kg)
    real(r8) :: forc_rainl           ! rainxy Atm flux mm/s
    real(r8) :: forc_snowc           ! snowfxy Atm flux  mm/s
    real(r8) :: forc_snowl           ! snowfxl Atm flux  mm/s
    real(r8) :: co2_ppmv_diag        ! temporary
    real(r8) :: co2_ppmv_prog        ! temporary
    real(r8) :: co2_ppmv_val         ! temporary
    integer  :: co2_type_idx         ! integer flag for co2_type options
    real(r8) :: esatw                ! saturation vapor pressure over water (Pa)
    real(r8) :: esati                ! saturation vapor pressure over ice (Pa)
    real(r8) :: a0,a1,a2,a3,a4,a5,a6 ! coefficients for esat over water
    real(r8) :: b0,b1,b2,b3,b4,b5,b6 ! coefficients for esat over ice
    real(r8) :: tdc, t               ! Kelvins to Celcius function and its input
    character(len=32), parameter :: sub = 'lnd_import_esmf'

    ! Constants to compute vapor pressure
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
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
!---------------------------------------------------------------------------
    co2_type_idx = 0
    if (co2_type == 'prognostic') then
       co2_type_idx = 1
    else if (co2_type == 'diagnostic') then
       co2_type_idx = 2
    end if
    if (co2_type == 'prognostic' .and. index_x2l_Sa_co2prog == 0) then
       call endrun( sub//' ERROR: must have nonzero index_x2l_Sa_co2prog for co2_type equal to prognostic' )
    else if (co2_type == 'diagnostic' .and. index_x2l_Sa_co2diag == 0) then
       call endrun( sub//' ERROR: must have nonzero index_x2l_Sa_co2diag for co2_type equal to diagnostic' )
    end if
    
    ! Note that the precipitation fluxes received  from the coupler
    ! are in units of kg/s/m^2. To convert these precipitation rates
    ! in units of mm/sec, one must divide by 1000 kg/m^3 and multiply
    ! by 1000 mm/m resulting in an overall factor of unity.
    ! Below the units are therefore given in mm/s.
    
    do g = begg,endg
        i = 1 + (g - begg)
       
        ! Determine required receive fields

        a2l%forc_hgt(g)     = fptr(index_x2l_Sa_z,i)         ! zgcmxy  Atm state m
        a2l%forc_u(g)       = fptr(index_x2l_Sa_u,i)         ! forc_uxy  Atm state m/s
        a2l%forc_v(g)       = fptr(index_x2l_Sa_v,i)         ! forc_vxy  Atm state m/s
        a2l%forc_th(g)      = fptr(index_x2l_Sa_ptem,i)      ! forc_thxy Atm state K
        a2l%forc_q(g)       = fptr(index_x2l_Sa_shum,i)      ! forc_qxy  Atm state kg/kg
        a2l%forc_pbot(g)    = fptr(index_x2l_Sa_pbot,i)      ! ptcmxy  Atm state Pa
        a2l%forc_t(g)       = fptr(index_x2l_Sa_tbot,i)      ! forc_txy  Atm state K
        a2l%forc_lwrad(g)   = fptr(index_x2l_Faxa_lwdn,i)    ! flwdsxy Atm flux  W/m^2
        forc_rainc          = fptr(index_x2l_Faxa_rainc,i)   ! mm/s
        forc_rainl          = fptr(index_x2l_Faxa_rainl,i)   ! mm/s
        forc_snowc          = fptr(index_x2l_Faxa_snowc,i)   ! mm/s
        forc_snowl          = fptr(index_x2l_Faxa_snowl,i)   ! mm/s
        a2l%forc_solad(g,2) = fptr(index_x2l_Faxa_swndr,i)   ! forc_sollxy  Atm flux  W/m^2
        a2l%forc_solad(g,1) = fptr(index_x2l_Faxa_swvdr,i)   ! forc_solsxy  Atm flux  W/m^2
        a2l%forc_solai(g,2) = fptr(index_x2l_Faxa_swndf,i)   ! forc_solldxy Atm flux  W/m^2
        a2l%forc_solai(g,1) = fptr(index_x2l_Faxa_swvdf,i)   ! forc_solsdxy Atm flux  W/m^2

        ! atmosphere coupling, for prognostic/prescribed aerosols
        a2l%forc_aer(g,1)  =  fptr(index_x2l_Faxa_bcphidry,i)
        a2l%forc_aer(g,2)  =  fptr(index_x2l_Faxa_bcphodry,i)
        a2l%forc_aer(g,3)  =  fptr(index_x2l_Faxa_bcphiwet,i)
        a2l%forc_aer(g,4)  =  fptr(index_x2l_Faxa_ocphidry,i)
        a2l%forc_aer(g,5)  =  fptr(index_x2l_Faxa_ocphodry,i)
        a2l%forc_aer(g,6)  =  fptr(index_x2l_Faxa_ocphiwet,i)
        a2l%forc_aer(g,7)  =  fptr(index_x2l_Faxa_dstwet1,i)
        a2l%forc_aer(g,8)  =  fptr(index_x2l_Faxa_dstdry1,i)
        a2l%forc_aer(g,9)  =  fptr(index_x2l_Faxa_dstwet2,i)
        a2l%forc_aer(g,10) =  fptr(index_x2l_Faxa_dstdry2,i)
        a2l%forc_aer(g,11) =  fptr(index_x2l_Faxa_dstwet3,i)
        a2l%forc_aer(g,12) =  fptr(index_x2l_Faxa_dstdry3,i)
        a2l%forc_aer(g,13) =  fptr(index_x2l_Faxa_dstwet4,i)
        a2l%forc_aer(g,14) =  fptr(index_x2l_Faxa_dstdry4,i)

        ! Determine optional receive fields

        if (index_x2l_Sa_co2prog /= 0) then
           co2_ppmv_prog = fptr(index_x2l_Sa_co2prog,i)   ! co2 atm state prognostic
        else
           co2_ppmv_prog = co2_ppmv
        end if
 
        if (index_x2l_Sa_co2diag /= 0) then
           co2_ppmv_diag = fptr(index_x2l_Sa_co2diag,i)   ! co2 atm state diagnostic
        else
           co2_ppmv_diag = co2_ppmv
        end if

        ! Determine derived quantities for required fields
        a2l%forc_hgt_u(g) = a2l%forc_hgt(g)    !observational height of wind [m]
        a2l%forc_hgt_t(g) = a2l%forc_hgt(g)    !observational height of temperature [m]
        a2l%forc_hgt_q(g) = a2l%forc_hgt(g)    !observational height of humidity [m]
        a2l%forc_vp(g)    = a2l%forc_q(g) * a2l%forc_pbot(g) &
                            / (0.622_r8 + 0.378_r8 * a2l%forc_q(g))
        a2l%forc_rho(g)   = (a2l%forc_pbot(g) - 0.378_r8 * a2l%forc_vp(g)) &
                            / (rair * a2l%forc_t(g))
        a2l%forc_po2(g)   = o2_molar_const * a2l%forc_pbot(g)
        a2l%forc_wind(g)  = sqrt(a2l%forc_u(g)**2 + a2l%forc_v(g)**2)
        a2l%forc_solar(g) = a2l%forc_solad(g,1) + a2l%forc_solai(g,1) + &
                            a2l%forc_solad(g,2) + a2l%forc_solai(g,2)
        a2l%forc_rain(g)  = forc_rainc + forc_rainl
        a2l%forc_snow(g)  = forc_snowc + forc_snowl
        a2l%rainf    (g)  = a2l%forc_rain(g) + a2l%forc_snow(g)

        if (a2l%forc_t(g) > SHR_CONST_TKFRZ) then
           e = esatw(tdc(a2l%forc_t(g)))
        else
           e = esati(tdc(a2l%forc_t(g)))
        end if
        qsat           = 0.622_r8*e / (a2l%forc_pbot(g) - 0.378_r8*e)
        a2l%forc_rh(g) = 100.0_r8*(a2l%forc_q(g) / qsat)
        ! Make sure relative humidity is properly bounded
        ! a2l%forc_rh(g) = min( 100.0_r8, a2l%forc_rh(g) )
        ! a2l%forc_rh(g) = max(   0.0_r8, a2l%forc_rh(g) )
        
        ! Determine derived quantities for optional fields
        ! Note that the following does unit conversions from ppmv to partial pressures (Pa)
        ! Note that forc_pbot is in Pa

        if (co2_type_idx == 1) then
           co2_ppmv_val = co2_ppmv_prog
        else if (co2_type_idx == 2) then
           co2_ppmv_val = co2_ppmv_diag 
        else
           co2_ppmv_val = co2_ppmv
        end if
        a2l%forc_pco2(g)   = co2_ppmv_val * 1.e-6_r8 * a2l%forc_pbot(g) 
#if (defined C13)
        a2l%forc_pc13o2(g) = co2_ppmv_val * c13ratio * 1.e-6_r8 * a2l%forc_pbot(g)
#endif

     end do

   end subroutine lnd_import_esmf

!===============================================================================

!---------------------------------------------------------------------------
!BOP
!
! !IROUTINE: lnd_domain_esmf
!
! !INTERFACE:
  subroutine lnd_domain_esmf( dom, rc )
!
! !DESCRIPTION:
!
! Send the land model domain information to the coupler
!
!---------------------------------------------------------------------------
! !USES:
    use shr_kind_mod, only : r8 => shr_kind_r8
    use clm_varcon  , only : re
    use domainMod   , only : adomain
    use decompMod   , only : get_proc_bounds_atm, adecomp
    implicit none
! !ARGUMENTS:
    type(ESMF_Array), intent(inout)     :: dom         ! CLM domain data
    integer, intent(out)                :: rc          ! return code
!
! !LOCAL VARIABLES:
    integer :: g,i                ! index
    integer :: begg, endg         ! beginning and ending gridcell indices
    integer :: klon,klat,karea,kmask,kfrac ! domain fields
    real(R8), pointer :: fptr (:,:)
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
!---------------------------------------------------------------------------
    !
    ! Initialize domain type
    ! lat/lon in degrees,  area in radians^2, mask is 1 (land), 0 (non-land)
    ! Note that in addition land carries around landfrac for the purposes of domain checking
    ! 
    rc = ESMF_SUCCESS

    call ESMF_ArrayGet(dom, localDe=0, farrayPtr=fptr, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)

    ! Fill in correct values for domain components
    klon  = esmfshr_util_ArrayGetIndex(dom,'lon ',rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
    klat  = esmfshr_util_ArrayGetIndex(dom,'lat ',rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
    karea = esmfshr_util_ArrayGetIndex(dom,'area',rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
    kmask = esmfshr_util_ArrayGetIndex(dom,'mask',rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
    kfrac = esmfshr_util_ArrayGetIndex(dom,'frac',rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)

    ! Determine bounds

    call get_proc_bounds_atm(begg, endg)

    ! Fill in correct values for domain components
    ! Note aream will be filled in in the atm-lnd mapper

    fptr(:,:) = -9999.0_R8
    fptr(kmask,:) = -0.0_R8
    do g = begg,endg
       i = 1 + (g - begg)
       fptr(klon, i) = adomain%lonc(g)
       fptr(klat, i) = adomain%latc(g)
       fptr(karea, i) = adomain%area(g)/(re*re)
       fptr(kmask, i) = real(adomain%mask(g), r8)
       fptr(kfrac, i) = real(adomain%frac(g), r8)
    end do

  end subroutine lnd_domain_esmf
    
!===============================================================================
    
#ifdef RTM
!---------------------------------------------------------------------------
!BOP
!
! !IROUTINE: rof_SetgsMap_esmf
!
! !INTERFACE:
  function rof_DistGrid_esmf(gsize, rc)
!
! !DESCRIPTION:
!
! Set the ESMF distributed grid for the runoff model
!
!---------------------------------------------------------------------------

! !USES:
    use shr_kind_mod, only : r8 => shr_kind_r8
    use clm_varpar  , only : rtmlon, rtmlat
    use RunoffMod   , only : runoff
    use abortutils  , only : endrun
    use clm_varctl  , only : iulog
    implicit none
! !RETURN:
    type(ESMF_DistGrid)                 :: rof_DistGrid_esmf   
! !ARGUMENTS:
    integer, intent(out)                :: gsize   ! Grid size
    integer, intent(out)                :: rc      ! Return code
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
!---------------------------------------------------------------------------
    !
    ! Local Variables
    !
    integer,allocatable :: gindex(:)         ! indexing for runoff grid cells
    integer :: n, ni                         ! indices
    integer :: lsize                         ! size of runoff data and number of grid cells
    integer :: ier                           ! error code
    character(len=32), parameter :: sub = 'rof_DistGrid_esmf'
    !-------------------------------------------------------------------

    ! Build the rof grid numbering
    ! NOTE:  Numbering scheme is: West to East and South to North
    ! starting at south pole.  Should be the same as what's used in SCRIP
    
    rc = ESMF_SUCCESS
    
    gsize = rtmlon*rtmlat
    lsize = runoff%lnumro
    allocate(gindex(lsize),stat=ier)

    ni = 0
    do n = runoff%begr,runoff%endr
       if (runoff%mask(n) == 2) then
          ni = ni + 1
          if (ni > runoff%lnumro) then
             write(iulog,*) sub, ' : ERROR runoff count',n,ni,runoff%lnumro
             call endrun( sub//' ERROR: runoff > expected' )
          endif
          gindex(ni) = runoff%gindex(n)
       endif
    end do
    if (ni /= runoff%lnumro) then
       write(iulog,*) sub, ' : ERROR runoff total count',ni,runoff%lnumro
       call endrun( sub//' ERROR: runoff not equal to expected' )
    endif

    rof_DistGrid_esmf = mct2esmf_init(gindex, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)

    deallocate(gindex)

  end function rof_DistGrid_esmf

!===============================================================================

!---------------------------------------------------------------------------
!BOP
!
! !IROUTINE: rof_domain_esmf
!
! !INTERFACE:
  subroutine rof_domain_esmf( dom, rc )
!
! !DESCRIPTION:
!
! Send the runoff model domain information to the coupler
!
!---------------------------------------------------------------------------

! !USES:
    use shr_kind_mod, only : r8 => shr_kind_r8
    use clm_varcon  , only : re
    use RunoffMod   , only : runoff
    use abortutils  , only : endrun
    use clm_varctl  , only : iulog
    use spmdMod     , only : iam
    implicit none
! !ARGUMENTS:
    type(ESMF_Array), intent(inout)     :: dom   ! Domain information
    integer, intent(out)                :: rc    ! Return code
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
!---------------------------------------------------------------------------
    !
    ! Local Variables
    !
    integer :: n, ni              ! index
    integer :: klon,klat,karea,kmask,kfrac ! domain fields
    real(R8),    pointer    :: fptr(:,:)
    character(len=32), parameter :: sub = 'rof_domain_esmf'
    !-------------------------------------------------------------------
    !
    ! Initialize domain type
    ! lat/lon in degrees,  area in radians^2, mask is 1 (land), 0 (non-land)
    ! Note that in addition land carries around landfrac for the purposes of domain checking
    ! 

    rc = ESMF_SUCCESS

    call ESMF_ArrayGet(dom, localDe=0, farrayPtr=fptr, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)

    ! Fill in correct values for domain components
    klon  = esmfshr_util_ArrayGetIndex(dom,'lon ',rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
    klat  = esmfshr_util_ArrayGetIndex(dom,'lat ',rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
    karea = esmfshr_util_ArrayGetIndex(dom,'area',rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
    kmask = esmfshr_util_ArrayGetIndex(dom,'mask',rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
    kfrac = esmfshr_util_ArrayGetIndex(dom,'frac',rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)

    ! Determine domain (numbering scheme is: West to East and South to North to South pole)
    ! Initialize attribute vector with special value

    ! Determine bounds numbering consistency

    ni = 0
    do n = runoff%begr,runoff%endr
       if (runoff%mask(n) == 2) then
          ni = ni + 1
          if (ni > runoff%lnumro) then
             write(iulog,*) sub, ' : ERROR runoff count',n,ni,runoff%lnumro
             call endrun( sub//' ERROR: runoff > expected' )
          endif
       end if
    end do
    if (ni /= runoff%lnumro) then
       write(iulog,*) sub, ' : ERROR runoff total count',ni,runoff%lnumro
       call endrun( sub//' ERROR: runoff not equal to expected' )
    endif

    ! Fill in correct values for domain components
    ! Note aream will be filled in in the atm-lnd mapper

    fptr(:,:) = -9999.0_R8
    ni = 0
    do n = runoff%begr,runoff%endr
       if (runoff%mask(n) == 2) then
          ni = ni + 1
          fptr(klon, ni) = runoff%lonc(n)
          fptr(klat, ni) = runoff%latc(n)
          fptr(karea, ni) = runoff%area(n)*1.0e-6_r8/(re*re)
          fptr(kmask, ni) = 1.0_r8
          fptr(kfrac, ni) = 1.0_r8
       end if
    end do

  end subroutine rof_domain_esmf

!====================================================================================

!---------------------------------------------------------------------------
!BOP
!
! !IROUTINE: rof_export_esmf
!
! !INTERFACE:
  subroutine rof_export_esmf( r2x_array, rc)
!
! !DESCRIPTION:
!
! Send the runoff model export state to the CESM coupler
!
!---------------------------------------------------------------------------
! !USES:
    use shr_kind_mod, only : r8 => shr_kind_r8
    use RunoffMod   , only : runoff, nt_rtm, rtm_tracers
    use abortutils  , only : endrun
    use clm_varctl  , only : iulog
#ifdef RTM
    use clm_varctl  , only : ice_runoff
#endif
    use clm_cpl_indices
    implicit none
! !ARGUMENTS:
    type(ESMF_Array), intent(inout)            :: r2x_array
    integer, intent(out)                       :: rc
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
!---------------------------------------------------------------------------
    !
    ! Local variables
    !
    integer :: ni, n, nt, nliq, nfrz
    real(R8), pointer :: fptr(:, :)
    character(len=*), parameter :: sub = 'rof_export_esmf'
    !-----------------------------------------------------
    
    rc = ESMF_SUCCESS

    call ESMF_ArrayGet(r2x_array, localDe=0, farrayPtr=fptr, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
    
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
       write(iulog,*)'RtmUpdateInput: ERROR in rtm_tracers LIQ ICE ',nliq,nfrz,rtm_tracers
       call endrun()
    endif

    ni = 0
    if ( ice_runoff )then
       do n = runoff%begr,runoff%endr
          if (runoff%mask(n) == 2) then
             ni = ni + 1
             ! liquid and ice runoff are treated separately
             fptr(index_r2x_Forr_roff,ni) = runoff%runoff(n,nliq)/(runoff%area(n)*1.0e-6_r8*1000._r8)
             fptr(index_r2x_Forr_ioff,ni) = runoff%runoff(n,nfrz)/(runoff%area(n)*1.0e-6_r8*1000._r8)
             if (ni > runoff%lnumro) then
                write(iulog,*) sub, ' : ERROR runoff count',n,ni
                call endrun( sub//' : ERROR runoff > expected' )
             endif
          endif
       end do
    else
       do n = runoff%begr,runoff%endr
          if (runoff%mask(n) == 2) then
             ni = ni + 1
             ! liquid and ice runoff are bundled together to liquid runoff, and then ice runoff set to zero
             fptr(index_r2x_Forr_roff,ni) =   &
               (runoff%runoff(n,nfrz)+runoff%runoff(n,nliq))/(runoff%area(n)*1.0e-6_r8*1000._r8)
             fptr(index_r2x_Forr_ioff,ni) = 0.0_r8
             if (ni > runoff%lnumro) then
                write(iulog,*) sub, ' : ERROR runoff count',n,ni
                call endrun( sub//' : ERROR runoff > expected' )
             endif
          endif
       end do
    end if
    if (ni /= runoff%lnumro) then
       write(iulog,*) sub, ' : ERROR runoff total count',ni,runoff%lnumro
       call endrun( sub//' : ERROR runoff not equal to expected' )
    endif

  end subroutine rof_export_esmf
#endif

!====================================================================================

  subroutine sno_export_esmf( s2x, s2x_array, rc )

    !-----------------------------------------------------
    use shr_kind_mod    , only : r8 => shr_kind_r8
    use clm_glclnd      , only : lnd2glc_type
    use decompMod       , only : get_proc_bounds_atm
    use clm_cpl_indices

    type(lnd2glc_type), intent(inout) :: s2x
    type(ESMF_Array)  , intent(inout) :: s2x_array
    integer, intent(out)              :: rc

    integer :: g,i
    integer :: begg, endg    ! beginning and ending gridcell indices
    real(r8), pointer :: fptr(:, :)
    !-----------------------------------------------------

    rc = ESMF_SUCCESS

    call get_proc_bounds_atm(begg, endg)

    call ESMF_ArrayGet(s2x_array, localDe=0, farrayPtr=fptr, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
   
    fptr(:,:) = 0.0_r8

    ! qice is positive if ice is growing, negative if melting

    ! For now, pass separate fields for each elevation class.
    ! Currently, the number of elevation classes must be 1, 3, 5, or 10.

    do g = begg,endg
       i = 1 + (g-begg)

#if (defined GLC_NEC_10 || defined GLC_NEC_5 || defined GLC_NEC_3 || defined GLC_NEC_1)
       fptr(index_s2x_Ss_tsrf01,i)   =  s2x%tsrf(g,1)
       fptr(index_s2x_Ss_topo01,i)   =  s2x%topo(g,1)
       fptr(index_s2x_Fgss_qice01,i) =  s2x%qice(g,1)
#endif
#if (defined GLC_NEC_10 || defined GLC_NEC_5 || defined GLC_NEC_3 )
       fptr(index_s2x_Ss_tsrf02,i)   =  s2x%tsrf(g,2)
       fptr(index_s2x_Ss_topo02,i)   =  s2x%topo(g,2)
       fptr(index_s2x_Fgss_qice02,i) =  s2x%qice(g,2)
       fptr(index_s2x_Ss_tsrf03,i)   =  s2x%tsrf(g,3)
       fptr(index_s2x_Ss_topo03,i)   =  s2x%topo(g,3)
       fptr(index_s2x_Fgss_qice03,i) =  s2x%qice(g,3)
#endif
#if (defined GLC_NEC_10 || defined GLC_NEC_5 )
       fptr(index_s2x_Ss_tsrf04,i)   =  s2x%tsrf(g,4)
       fptr(index_s2x_Ss_topo04,i)   =  s2x%topo(g,4)
       fptr(index_s2x_Fgss_qice04,i) =  s2x%qice(g,4)
       fptr(index_s2x_Ss_tsrf05,i)   =  s2x%tsrf(g,5)
       fptr(index_s2x_Ss_topo05,i)   =  s2x%topo(g,5)
       fptr(index_s2x_Fgss_qice05,i) =  s2x%qice(g,5)
#endif
#if (defined GLC_NEC_10 )
       fptr(index_s2x_Ss_tsrf06,i)   =  s2x%tsrf(g,6)
       fptr(index_s2x_Ss_topo06,i)   =  s2x%topo(g,6)
       fptr(index_s2x_Fgss_qice06,i) =  s2x%qice(g,6)
       fptr(index_s2x_Ss_tsrf07,i)   =  s2x%tsrf(g,7)
       fptr(index_s2x_Ss_topo07,i)   =  s2x%topo(g,7)
       fptr(index_s2x_Fgss_qice07,i) =  s2x%qice(g,7)
       fptr(index_s2x_Ss_tsrf08,i)   =  s2x%tsrf(g,8)
       fptr(index_s2x_Ss_topo08,i)   =  s2x%topo(g,8)
       fptr(index_s2x_Fgss_qice08,i) =  s2x%qice(g,8)
       fptr(index_s2x_Ss_tsrf09,i)   =  s2x%tsrf(g,9)
       fptr(index_s2x_Ss_topo09,i)   =  s2x%topo(g,9)
       fptr(index_s2x_Fgss_qice09,i) =  s2x%qice(g,9)
       fptr(index_s2x_Ss_tsrf10,i)   =  s2x%tsrf(g,10)
       fptr(index_s2x_Ss_topo10,i)   =  s2x%topo(g,10)
       fptr(index_s2x_Fgss_qice10,i) =  s2x%qice(g,10)
#endif

    end do   ! g

  end subroutine sno_export_esmf

!====================================================================================
  subroutine sno_import_esmf( x2s_array, x2s, rc )

    !-----------------------------------------------------
    use shr_kind_mod    , only: r8 => shr_kind_r8
    use clm_glclnd      , only: glc2lnd_type
    use decompMod       , only: get_proc_bounds_atm
    use abortutils      , only: endrun
    use clm_varctl      , only: iulog
    use clm_cpl_indices

    !
    ! Arguments
    !
    type(ESMF_Array)  , intent(inout) :: x2s_array
    type(glc2lnd_type), intent(inout) :: x2s
    integer, intent(out)              :: rc
    !
    ! Local Variables
    !
    integer  :: g,i
    integer  :: begg, endg           ! beginning and ending gridcell indices
    real(r8), pointer :: fptr(:, :)
    character(len=32), parameter :: sub = 'sno_import_esmf'

    !-----------------------------------------------------

    rc = ESMF_SUCCESS

    call ESMF_ArrayGet(x2s_array, localDe=0, farrayPtr=fptr, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)

    call get_proc_bounds_atm(begg, endg)

    do g = begg,endg
        i = 1 + (g - begg)

#if (defined GLC_NEC_10 || defined GLC_NEC_5 || defined GLC_NEC_3 || defined GLC_NEC_1)
        x2s%frac(g,1)  = fptr(index_x2s_Sg_frac01,i)
        x2s%topo(g,1)  = fptr(index_x2s_Sg_topo01,i)
        x2s%rofi(g,1)  = fptr(index_x2s_Fsgg_rofi01,i)
        x2s%rofl(g,1)  = fptr(index_x2s_Fsgg_rofl01,i)
        x2s%hflx(g,1)  = fptr(index_x2s_Fsgg_hflx01,i)
#endif
#if (defined GLC_NEC_10 || defined GLC_NEC_5 || defined GLC_NEC_3 )
        x2s%frac(g,2)  = fptr(index_x2s_Sg_frac02,i)
        x2s%topo(g,2)  = fptr(index_x2s_Sg_topo02,i)
        x2s%rofi(g,2)  = fptr(index_x2s_Fsgg_rofi02,i)
        x2s%rofl(g,2)  = fptr(index_x2s_Fsgg_rofl02,i)
        x2s%hflx(g,2)  = fptr(index_x2s_Fsgg_hflx02,i)
        x2s%frac(g,3)  = fptr(index_x2s_Sg_frac03,i)
        x2s%topo(g,3)  = fptr(index_x2s_Sg_topo03,i)
        x2s%rofi(g,3)  = fptr(index_x2s_Fsgg_rofi03,i)
        x2s%rofl(g,3)  = fptr(index_x2s_Fsgg_rofl03,i)
        x2s%hflx(g,3)  = fptr(index_x2s_Fsgg_hflx03,i)
#endif
#if (defined GLC_NEC_10 || defined GLC_NEC_5 )
        x2s%frac(g,4)  = fptr(index_x2s_Sg_frac04,i)
        x2s%topo(g,4)  = fptr(index_x2s_Sg_topo04,i)
        x2s%rofi(g,4)  = fptr(index_x2s_Fsgg_rofi04,i)
        x2s%rofl(g,4)  = fptr(index_x2s_Fsgg_rofl04,i)
        x2s%hflx(g,4)  = fptr(index_x2s_Fsgg_hflx04,i)
        x2s%frac(g,5)  = fptr(index_x2s_Sg_frac05,i)
        x2s%topo(g,5)  = fptr(index_x2s_Sg_topo05,i)
        x2s%rofi(g,5)  = fptr(index_x2s_Fsgg_rofi05,i)
        x2s%rofl(g,5)  = fptr(index_x2s_Fsgg_rofl05,i)
        x2s%hflx(g,5)  = fptr(index_x2s_Fsgg_hflx05,i)
#endif
#if (defined GLC_NEC_10 )
        x2s%frac(g,6)  = fptr(index_x2s_Sg_frac06,i)
        x2s%topo(g,6)  = fptr(index_x2s_Sg_topo06,i)
        x2s%rofi(g,6)  = fptr(index_x2s_Fsgg_rofi06,i)
        x2s%rofl(g,6)  = fptr(index_x2s_Fsgg_rofl06,i)
        x2s%hflx(g,6)  = fptr(index_x2s_Fsgg_hflx06,i)
        x2s%frac(g,7)  = fptr(index_x2s_Sg_frac07,i)
        x2s%topo(g,7)  = fptr(index_x2s_Sg_topo07,i)
        x2s%rofi(g,7)  = fptr(index_x2s_Fsgg_rofi07,i)
        x2s%rofl(g,7)  = fptr(index_x2s_Fsgg_rofl07,i)
        x2s%hflx(g,7)  = fptr(index_x2s_Fsgg_hflx07,i)
        x2s%frac(g,8)  = fptr(index_x2s_Sg_frac08,i)
        x2s%topo(g,8)  = fptr(index_x2s_Sg_topo08,i)
        x2s%rofi(g,8)  = fptr(index_x2s_Fsgg_rofi08,i)
        x2s%rofl(g,8)  = fptr(index_x2s_Fsgg_rofl08,i)
        x2s%hflx(g,8)  = fptr(index_x2s_Fsgg_hflx08,i)
        x2s%frac(g,9)  = fptr(index_x2s_Sg_frac09,i)
        x2s%topo(g,9)  = fptr(index_x2s_Sg_topo09,i)
        x2s%rofi(g,9)  = fptr(index_x2s_Fsgg_rofi09,i)
        x2s%rofl(g,9)  = fptr(index_x2s_Fsgg_rofl09,i)
        x2s%hflx(g,9)  = fptr(index_x2s_Fsgg_hflx09,i)
        x2s%frac(g,10) = fptr(index_x2s_Sg_frac10,i)
        x2s%topo(g,10) = fptr(index_x2s_Sg_topo10,i)
        x2s%rofi(g,10) = fptr(index_x2s_Fsgg_rofi10,i)
        x2s%rofl(g,10) = fptr(index_x2s_Fsgg_rofl10,i)
        x2s%hflx(g,10) = fptr(index_x2s_Fsgg_hflx10,i)
#endif

     end do    ! g


   end subroutine sno_import_esmf

!====================================================================================

end module lnd_comp_esmf
