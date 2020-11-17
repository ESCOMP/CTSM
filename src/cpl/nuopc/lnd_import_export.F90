module lnd_import_export

  use ESMF                  , only : ESMF_GridComp, ESMF_State, ESMF_Mesh, ESMF_StateGet
  use ESMF                  , only : ESMF_KIND_R8, ESMF_SUCCESS, ESMF_MAXSTR, ESMF_LOGMSG_INFO
  use ESMF                  , only : ESMF_LogWrite, ESMF_LOGMSG_ERROR, ESMF_LogFoundError
  use ESMF                  , only : ESMF_STATEITEM_NOTFOUND, ESMF_StateItem_Flag
  use ESMF                  , only : operator(/=), operator(==)
  use NUOPC                 , only : NUOPC_CompAttributeGet, NUOPC_Advertise, NUOPC_IsConnected
  use NUOPC_Model           , only : NUOPC_ModelGet
  use shr_kind_mod          , only : r8 => shr_kind_r8, cx=>shr_kind_cx, cxx=>shr_kind_cxx, cs=>shr_kind_cs
  use shr_sys_mod           , only : shr_sys_abort
  use clm_varctl            , only : iulog
  use clm_time_manager      , only : get_nstep
  use decompmod             , only : bounds_type
  use lnd2atmType           , only : lnd2atm_type
  use lnd2glcMod            , only : lnd2glc_type
  use atm2lndType           , only : atm2lnd_type
  use glc2lndMod            , only : glc2lnd_type
  use domainMod             , only : ldomain
  use spmdMod               , only : masterproc
  use seq_drydep_mod        , only : seq_drydep_readnl, n_drydep
  use shr_megan_mod         , only : shr_megan_readnl, shr_megan_mechcomps_n
  use shr_fire_emis_mod     , only : shr_fire_emis_readnl
  use shr_carma_mod         , only : shr_carma_readnl
  use shr_ndep_mod          , only : shr_ndep_readnl
  use nuopc_shr_methods     , only : chkerr
  use lnd_import_export_utils, only : derive_quantities, check_for_errors, check_for_nans

  implicit none
  private ! except

  public  :: advertise_fields
  public  :: realize_fields
  public  :: import_fields
  public  :: export_fields

  private :: fldlist_add
  private :: fldlist_realize
  private :: state_getimport
  private :: state_setexport
  private :: state_getfldptr

  type fld_list_type
     character(len=128) :: stdname
     integer :: ungridded_lbound = 0
     integer :: ungridded_ubound = 0
  end type fld_list_type

  integer, parameter     :: fldsMax = 100
  integer                :: fldsToLnd_num = 0
  integer                :: fldsFrLnd_num = 0
  type (fld_list_type)   :: fldsToLnd(fldsMax)
  type (fld_list_type)   :: fldsFrLnd(fldsMax)
  integer, parameter     :: gridTofieldMap = 2 ! ungridded dimension is innermost

  ! from atm->lnd
  integer                :: ndep_nflds       ! number  of nitrogen deposition fields from atm->lnd/ocn

  ! from lnd->atm
  character(len=cx)      :: carma_fields     ! List of CARMA fields from lnd->atm
  integer                :: drydep_nflds     ! number of dry deposition velocity fields lnd-> atm
  integer                :: megan_nflds      ! number of MEGAN voc fields from lnd-> atm
  integer                :: emis_nflds       ! number of fire emission fields from lnd-> atm

  logical                :: flds_co2a        ! use case
  logical                :: flds_co2b        ! use case
  logical                :: flds_co2c        ! use case
  integer                :: glc_nec          ! number of glc elevation classes
  integer, parameter     :: debug = 0        ! internal debug level

  character(*),parameter :: F01 = "('(lnd_import_export) ',a,i5,2x,i5,2x,d21.14)"
  character(*),parameter :: u_FILE_u = &
       __FILE__

!===============================================================================
contains
!===============================================================================

  subroutine advertise_fields(gcomp, flds_scalar_name, glc_present, cism_evolve, rof_prognostic, rc)

    use clm_varctl, only : ndep_from_cpl

    ! input/output variables
    type(ESMF_GridComp)            :: gcomp
    character(len=*) , intent(in)  :: flds_scalar_name
    logical          , intent(in)  :: glc_present
    logical          , intent(in)  :: cism_evolve
    logical          , intent(in)  :: rof_prognostic
    integer          , intent(out) :: rc

    ! local variables
    type(ESMF_State)       :: importState
    type(ESMF_State)       :: exportState
    character(ESMF_MAXSTR) :: stdname
    character(ESMF_MAXSTR) :: cvalue
    character(len=2)       :: nec_str
    integer                :: n, num
    character(len=128)     :: fldname
    character(len=*), parameter :: subname='(lnd_import_export:advertise_fields)'
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    call NUOPC_ModelGet(gcomp, importState=importState, exportState=exportState, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !--------------------------------
    ! determine necessary toggles for below
    !--------------------------------

    call NUOPC_CompAttributeGet(gcomp, name='flds_co2a', value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) flds_co2a
    call ESMF_LogWrite('flds_co2a = '// trim(cvalue), ESMF_LOGMSG_INFO)

    call NUOPC_CompAttributeGet(gcomp, name='flds_co2b', value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) flds_co2b
    call ESMF_LogWrite('flds_co2b = '// trim(cvalue), ESMF_LOGMSG_INFO)

    call NUOPC_CompAttributeGet(gcomp, name='flds_co2c', value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) flds_co2c
    call ESMF_LogWrite('flds_co2c = '// trim(cvalue), ESMF_LOGMSG_INFO)

    ! Determine  number of elevation classes
    call NUOPC_CompAttributeGet(gcomp, name='glc_nec', value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) glc_nec
    call ESMF_LogWrite('glc_nec = '// trim(cvalue), ESMF_LOGMSG_INFO)
    if (glc_nec < 1) then
       call shr_sys_abort('ERROR: In CLM4.5 and later, glc_nec must be at least 1.')
    end if

    !--------------------------------
    ! Advertise export fields
    !--------------------------------

    call fldlist_add(fldsFrLnd_num, fldsFrlnd, trim(flds_scalar_name))

    ! export land states
    call fldlist_add(fldsFrLnd_num, fldsFrlnd, 'Sl_lfrin'      )
    call fldlist_add(fldsFrLnd_num, fldsFrlnd, 'Sl_t'          )
    call fldlist_add(fldsFrLnd_num, fldsFrlnd, 'Sl_tref'       )
    call fldlist_add(fldsFrLnd_num, fldsFrlnd, 'Sl_qref'       )
    call fldlist_add(fldsFrLnd_num, fldsFrlnd, 'Sl_avsdr'      )
    call fldlist_add(fldsFrLnd_num, fldsFrlnd, 'Sl_anidr'      )
    call fldlist_add(fldsFrLnd_num, fldsFrlnd, 'Sl_avsdf'      )
    call fldlist_add(fldsFrLnd_num, fldsFrlnd, 'Sl_anidf'      )
    call fldlist_add(fldsFrLnd_num, fldsFrlnd, 'Sl_snowh'      )
    call fldlist_add(fldsFrLnd_num, fldsFrlnd, 'Sl_u10'        )
    call fldlist_add(fldsFrLnd_num, fldsFrlnd, 'Sl_fv'         )
    call fldlist_add(fldsFrLnd_num, fldsFrlnd, 'Sl_ram1'       )

    ! export fluxes to river
    if (rof_prognostic) then
       call fldlist_add(fldsFrLnd_num, fldsFrlnd, 'Flrl_rofsur'   )
       call fldlist_add(fldsFrLnd_num, fldsFrlnd, 'Flrl_rofgwl'   )
       call fldlist_add(fldsFrLnd_num, fldsFrlnd, 'Flrl_rofsub'   )
       call fldlist_add(fldsFrLnd_num, fldsFrlnd, 'Flrl_rofi'     )
       call fldlist_add(fldsFrLnd_num, fldsFrlnd, 'Flrl_irrig'    )
    end if

    ! export fluxes to atm
    call fldlist_add(fldsFrLnd_num, fldsFrlnd, 'Fall_taux'     )
    call fldlist_add(fldsFrLnd_num, fldsFrlnd, 'Fall_tauy'     )
    call fldlist_add(fldsFrLnd_num, fldsFrlnd, 'Fall_lat'      )
    call fldlist_add(fldsFrLnd_num, fldsFrlnd, 'Fall_sen'      )
    call fldlist_add(fldsFrLnd_num, fldsFrlnd, 'Fall_lwup'     )
    call fldlist_add(fldsFrLnd_num, fldsFrlnd, 'Fall_evap'     )
    call fldlist_add(fldsFrLnd_num, fldsFrlnd, 'Fall_swnet'    )

    ! call fldlist_add(fldsFrLnd_num, fldsFrlnd, 'Fall_methane'  )

    ! dust fluxes from land (4 sizes)
    call fldlist_add(fldsFrLnd_num, fldsFrLnd, 'Fall_flxdst', ungridded_lbound=1, ungridded_ubound=4)

    ! co2 fields from land
    if (flds_co2b .or. flds_co2c) then
       call fldlist_add(fldsFrLnd_num, fldsFrlnd, 'Fall_fco2_lnd' )
    end if

    ! Dry Deposition velocities from land - ALSO initialize drydep here
    call seq_drydep_readnl("drv_flds_in", drydep_nflds)
    if (drydep_nflds > 0) then
       call fldlist_add(fldsFrLnd_num, fldsFrLnd, 'Sl_ddvel', ungridded_lbound=1, ungridded_ubound=drydep_nflds)
    end if

    ! MEGAN VOC emissions fluxes from land
    call shr_megan_readnl('drv_flds_in', megan_nflds)
    if (shr_megan_mechcomps_n .ne. megan_nflds) call shr_sys_abort('ERROR: megan field count mismatch')
    if (shr_megan_mechcomps_n > 0) then
       call fldlist_add(fldsFrLnd_num, fldsFrLnd, 'Fall_voc', ungridded_lbound=1, ungridded_ubound=megan_nflds)
    end if

    ! Fire emissions fluxes from land
    call shr_fire_emis_readnl('drv_flds_in', emis_nflds)
    if (emis_nflds > 0) then
       call fldlist_add(fldsFrLnd_num, fldsFrLnd, 'Fall_fire', ungridded_lbound=1, ungridded_ubound=emis_nflds)
       call fldlist_add(fldsFrLnd_num, fldsFrLnd, 'Sl_fztop')
    end if
    ! CARMA volumetric soil water from land
    ! TODO: is the following correct - the CARMA field exchange is very confusing in mct
    call shr_carma_readnl('drv_flds_in', carma_fields)
    if (carma_fields /= ' ') then
       call fldlist_add(fldsFrLnd_num, fldsFrlnd, 'Sl_soilw') ! optional for carma
    end if

    if (glc_present .and. cism_evolve) then
       ! lnd->glc states from land all lnd->glc elevation classes (1:glc_nec) plus bare land (index 0).
       ! The following puts all of the elevation class fields as an
       ! undidstributed dimension in the export state field

       call fldlist_add(fldsFrLnd_num, fldsFrLnd, 'Sl_tsrf_elev'  , ungridded_lbound=1, ungridded_ubound=glc_nec+1)
       call fldlist_add(fldsFrLnd_num, fldsFrLnd, 'Sl_topo_elev'  , ungridded_lbound=1, ungridded_ubound=glc_nec+1)
       call fldlist_add(fldsFrLnd_num, fldsFrLnd, 'Flgl_qice_elev', ungridded_lbound=1, ungridded_ubound=glc_nec+1)
    end if

    ! Now advertise above export fields
    do n = 1,fldsFrLnd_num
       call NUOPC_Advertise(exportState, standardName=fldsFrLnd(n)%stdname, &
            TransferOfferGeomObject='will provide', rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    enddo

    !--------------------------------
    ! Advertise import fields
    !--------------------------------

    call fldlist_add(fldsToLnd_num, fldsToLnd, trim(flds_scalar_name))

    ! from atm - states
    call fldlist_add(fldsToLnd_num, fldsToLnd, 'Sa_z'         )
    call fldlist_add(fldsToLnd_num, fldsToLnd, 'Sa_topo'      )
    call fldlist_add(fldsToLnd_num, fldsToLnd, 'Sa_u'         )
    call fldlist_add(fldsToLnd_num, fldsToLnd, 'Sa_v'         )
    call fldlist_add(fldsToLnd_num, fldsToLnd, 'Sa_ptem'      )
    call fldlist_add(fldsToLnd_num, fldsToLnd, 'Sa_pbot'      )
    call fldlist_add(fldsToLnd_num, fldsToLnd, 'Sa_tbot'      )
    call fldlist_add(fldsToLnd_num, fldsToLnd, 'Sa_shum'      )
   !call fldlist_add(fldsToLnd_num, fldsToLnd, 'Sa_methane'   )

    ! from atm - fluxes
    call fldlist_add(fldsToLnd_num, fldsToLnd, 'Faxa_lwdn'    )
    call fldlist_add(fldsToLnd_num, fldsToLnd, 'Faxa_rainc'   )
    call fldlist_add(fldsToLnd_num, fldsToLnd, 'Faxa_rainl'   )
    call fldlist_add(fldsToLnd_num, fldsToLnd, 'Faxa_snowc'   )
    call fldlist_add(fldsToLnd_num, fldsToLnd, 'Faxa_snowl'   )
    call fldlist_add(fldsToLnd_num, fldsToLnd, 'Faxa_swndr'   )
    call fldlist_add(fldsToLnd_num, fldsToLnd, 'Faxa_swvdr'   )
    call fldlist_add(fldsToLnd_num, fldsToLnd, 'Faxa_swndf'   )
    call fldlist_add(fldsToLnd_num, fldsToLnd, 'Faxa_swvdf'   )

    ! from atm - black carbon deposition fluxes (3)
    ! (1) => bcphidry, (2) => bcphodry, (3) => bcphiwet
    call fldlist_add(fldsToLnd_num, fldsToLnd, 'Faxa_bcph',  ungridded_lbound=1, ungridded_ubound=3)

    ! from atm - organic carbon deposition fluxes (3)
    ! (1) => ocphidry, (2) => ocphodry, (3) => ocphiwet
    call fldlist_add(fldsToLnd_num, fldsToLnd, 'Faxa_ocph',  ungridded_lbound=1, ungridded_ubound=3)

    ! from atm - wet dust deposition frluxes (4 sizes)
    ! (1) => dstwet1, (2) => dstwet2, (3) => dstwet3, (4) => dstwet4
    call fldlist_add(fldsToLnd_num, fldsToLnd, 'Faxa_dstwet', ungridded_lbound=1, ungridded_ubound=4)

    ! from - atm dry dust deposition frluxes (4 sizes)
    call fldlist_add(fldsToLnd_num, fldsToLnd, 'Faxa_dstdry', ungridded_lbound=1, ungridded_ubound=4)

    ! from atm - nitrogen deposition
    call shr_ndep_readnl("drv_flds_in", ndep_nflds)
    if (ndep_nflds > 0) then
       call fldlist_add(fldsToLnd_num, fldsToLnd, 'Faxa_ndep', ungridded_lbound=1, ungridded_ubound=ndep_nflds)
       ! This sets a variable in clm_varctl
       ndep_from_cpl = .true.
    end if

    ! from atm - co2 exchange scenarios
    if (flds_co2a .or. flds_co2b .or. flds_co2c) then
       call fldlist_add(fldsToLnd_num, fldsToLnd, 'Sa_co2prog')
       call fldlist_add(fldsToLnd_num, fldsToLnd, 'Sa_co2diag')
    end if

    if (rof_prognostic) then
       ! from river
       call fldlist_add(fldsToLnd_num, fldsToLnd, 'Flrr_flood'   )
       call fldlist_add(fldsToLnd_num, fldsToLnd, 'Flrr_volr'    )
       call fldlist_add(fldsToLnd_num, fldsToLnd, 'Flrr_volrmch' )
    end if

    if (glc_present) then
       ! from land-ice (glc) - no elevation classes
       call fldlist_add(fldsToLnd_num, fldsToLnd, 'Sg_icemask'               ) ! mask of where cism is running
       call fldlist_add(fldsToLnd_num, fldsToLnd, 'Sg_icemask_coupled_fluxes') !

       ! from land-ice (glc) - fields for all glc->lnd elevation classes (1:glc_nec) plus bare land (index 0)
       call fldlist_add(fldsToLnd_num, fldsToLnd, 'Sg_ice_covered_elev', ungridded_lbound=1, ungridded_ubound=glc_nec+1)
       call fldlist_add(fldsToLnd_num, fldsToLnd, 'Sg_topo_elev'       , ungridded_lbound=1, ungridded_ubound=glc_nec+1)

       !current not used - but could be used in the future
       !call fldlist_add(fldsToLnd_num, fldsToLnd, 'Flgg_hflx_elev'     , ungridded_lbound=1, ungridded_ubound=glc_nec+1)
    end if

    ! Now advertise import fields
    do n = 1,fldsToLnd_num
       call NUOPC_Advertise(importState, standardName=fldsToLnd(n)%stdname, &
            TransferOfferGeomObject='will provide', rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    enddo

  end subroutine advertise_fields

  !===============================================================================

  subroutine realize_fields(gcomp, Emesh, flds_scalar_name, flds_scalar_num, rc)

    ! input/output variables
    type(ESMF_GridComp) , intent(inout) :: gcomp
    type(ESMF_Mesh)     , intent(in)    :: Emesh
    character(len=*)    , intent(in)    :: flds_scalar_name
    integer             , intent(in)    :: flds_scalar_num
    integer             , intent(out)   :: rc

    ! local variables
    type(ESMF_State)     :: importState
    type(ESMF_State)     :: exportState
    character(len=*), parameter :: subname='(lnd_import_export:realize_fields)'
    !---------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    call NUOPC_ModelGet(gcomp, importState=importState, exportState=exportState, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call fldlist_realize( &
         state=ExportState, &
         fldList=fldsFrLnd, &
         numflds=fldsFrLnd_num, &
         flds_scalar_name=flds_scalar_name, &
         flds_scalar_num=flds_scalar_num, &
         tag=subname//':clmExport',&
         mesh=Emesh, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call fldlist_realize( &
         state=importState, &
         fldList=fldsToLnd, &
         numflds=fldsToLnd_num, &
         flds_scalar_name=flds_scalar_name, &
         flds_scalar_num=flds_scalar_num, &
         tag=subname//':clmImport',&
         mesh=Emesh, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

  end subroutine realize_fields

  !===============================================================================

  subroutine import_fields( gcomp, bounds, glc_present, rof_prognostic, &
       atm2lnd_inst, glc2lnd_inst, wateratm2lndbulk_inst, rc)

    !---------------------------------------------------------------------------
    ! Convert the input data from the mediator to the land model
    !---------------------------------------------------------------------------

    use clm_varctl           , only: co2_type, co2_ppmv, use_c13, ndep_from_cpl
    use clm_varcon           , only: rair, o2_molar_const, c13ratio
    use shr_const_mod        , only: SHR_CONST_TKFRZ
    use Wateratm2lndBulkType , only: wateratm2lndbulk_type
    use QSatMod              , only: QSat

    ! input/output variabes
    type(ESMF_GridComp)                         :: gcomp
    type(bounds_type)           , intent(in)    :: bounds       ! bounds
    logical                     , intent(in)    :: glc_present    ! .true. => running with a non-stub GLC model
    logical                     , intent(in)    :: rof_prognostic ! .true. => running with a prognostic ROF model
    type(atm2lnd_type)          , intent(inout) :: atm2lnd_inst ! clm internal input data type
    type(glc2lnd_type)          , intent(inout) :: glc2lnd_inst ! clm internal input data type
    type(Wateratm2lndbulk_type) , intent(inout) :: wateratm2lndbulk_inst
    integer                     , intent(out)   :: rc

    ! local variables
    type(ESMF_State)          :: importState
    type(ESMF_StateItem_Flag) :: itemFlag
    real(r8), pointer         :: dataPtr(:)
    character(len=128)        :: fldname
    integer                   :: num
    integer                   :: begg, endg                          ! bounds
    integer                   :: g,i,k                               ! indices
    real(r8)                  :: qsat_kg_kg                          ! saturation specific humidity (kg/kg)
    real(r8)                  :: forc_pbot                           ! atmospheric pressure (Pa)
    real(r8)                  :: co2_ppmv_input(bounds%begg:bounds%endg)   ! temporary
    real(r8)                  :: forc_rainc(bounds%begg:bounds%endg) ! rainxy Atm flux mm/s
    real(r8)                  :: forc_rainl(bounds%begg:bounds%endg) ! rainxy Atm flux mm/s
    real(r8)                  :: forc_snowc(bounds%begg:bounds%endg) ! snowfxy Atm flux  mm/s
    real(r8)                  :: forc_snowl(bounds%begg:bounds%endg) ! snowfxl Atm flux  mm/s
    real(r8)                  :: forc_noy(bounds%begg:bounds%endg)
    real(r8)                  :: forc_nhx(bounds%begg:bounds%endg)
    real(r8)                  :: frac_grc(bounds%begg:bounds%endg, 0:glc_nec)
    real(r8)                  :: topo_grc(bounds%begg:bounds%endg, 0:glc_nec)
    real(r8)                  :: hflx_grc(bounds%begg:bounds%endg, 0:glc_nec)
    real(r8)                  :: icemask_grc(bounds%begg:bounds%endg)
    real(r8)                  :: icemask_coupled_fluxes_grc(bounds%begg:bounds%endg)
    character(len=*), parameter :: subname='(lnd_import_export:import_fields)'

    !---------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    ! Get import state
    call NUOPC_ModelGet(gcomp, importState=importState, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Set bounds
    begg = bounds%begg; endg=bounds%endg

    ! Note: precipitation fluxes received  from the coupler
    ! are in units of kg/s/m^2. To convert these precipitation rates
    ! in units of mm/sec, one must divide by 1000 kg/m^3 and multiply
    ! by 1000 mm/m resulting in an overall factor of unity.
    ! Below the units are therefore given in mm/s.

    !--------------------------
    ! Required atmosphere input fields
    !--------------------------

    call state_getimport(importState, 'Sa_z', bounds, output=atm2lnd_inst%forc_hgt_grc, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_getimport(importState, 'Sa_topo', bounds, output=atm2lnd_inst%forc_topo_grc, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_getimport(importState, 'Sa_u', bounds, output=atm2lnd_inst%forc_u_grc, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_getimport(importState, 'Sa_v', bounds, output=atm2lnd_inst%forc_v_grc, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_getimport(importState, 'Sa_ptem', bounds, output=atm2lnd_inst%forc_th_not_downscaled_grc, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_getimport(importState, 'Sa_shum', bounds, output=wateratm2lndbulk_inst%forc_q_not_downscaled_grc, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_getimport(importState, 'Sa_pbot', bounds, output=atm2lnd_inst%forc_pbot_not_downscaled_grc, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_getimport(importState, 'Sa_tbot', bounds, output=atm2lnd_inst%forc_t_not_downscaled_grc, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_getimport(importState, 'Faxa_rainc', bounds, output=forc_rainc, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_getimport(importState, 'Faxa_rainl', bounds, output=forc_rainl, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_getimport(importState, 'Faxa_snowc', bounds, output=forc_snowc, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_getimport(importState, 'Faxa_snowl', bounds, output=forc_snowl, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_getimport(importState, 'Faxa_lwdn', bounds, output=atm2lnd_inst%forc_lwrad_not_downscaled_grc, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_getimport(importState, 'Faxa_swvdr', bounds, output=atm2lnd_inst%forc_solad_grc(:,1), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_getimport(importState, 'Faxa_swndr', bounds, output=atm2lnd_inst%forc_solad_grc(:,2), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_getimport(importState, 'Faxa_swvdf', bounds, output=atm2lnd_inst%forc_solai_grc(:,1), rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_getimport(importState, 'Faxa_swndf', bounds, output=atm2lnd_inst%forc_solai_grc(:,2), rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Atmosphere prognostic/prescribed aerosol fields

    ! bcphidry
    call state_getimport(importState, 'Faxa_bcph', bounds, output=atm2lnd_inst%forc_aer_grc(:,1), &
         ungridded_index=1, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    ! bcphodry
    call state_getimport(importState, 'Faxa_bcph', bounds, output=atm2lnd_inst%forc_aer_grc(:,2), &
         ungridded_index=2, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    ! bcphiwet
    call state_getimport(importState, 'Faxa_bcph', bounds, output=atm2lnd_inst%forc_aer_grc(:,3), &
         ungridded_index=3, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! ocphidry
    call state_getimport(importState, 'Faxa_ocph', bounds, output=atm2lnd_inst%forc_aer_grc(:,4), &
         ungridded_index=1, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    ! ocphodry
    call state_getimport(importState, 'Faxa_ocph', bounds, output=atm2lnd_inst%forc_aer_grc(:,5), &
         ungridded_index=2, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    ! ocphiwet
    call state_getimport(importState, 'Faxa_ocph', bounds, output=atm2lnd_inst%forc_aer_grc(:,6), &
         ungridded_index=3, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_getimport(importState, 'Faxa_dstwet', bounds, output=atm2lnd_inst%forc_aer_grc(:,7), &
         ungridded_index=1, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport(importState, 'Faxa_dstdry', bounds, output=atm2lnd_inst%forc_aer_grc(:,8), &
         ungridded_index=1, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_getimport(importState, 'Faxa_dstwet', bounds, output=atm2lnd_inst%forc_aer_grc(:,9), &
         ungridded_index=2, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport(importState, 'Faxa_dstdry', bounds, output=atm2lnd_inst%forc_aer_grc(:,10), &
         ungridded_index=2, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_getimport(importState, 'Faxa_dstwet', bounds, output=atm2lnd_inst%forc_aer_grc(:,11), &
         ungridded_index=3, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport(importState, 'Faxa_dstdry', bounds, output=atm2lnd_inst%forc_aer_grc(:,12), &
         ungridded_index=3, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_getimport(importState, 'Faxa_dstwet', bounds, output=atm2lnd_inst%forc_aer_grc(:,13), &
         ungridded_index=4, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport(importState, 'Faxa_dstdry', bounds, output=atm2lnd_inst%forc_aer_grc(:,14), &
         ungridded_index=4, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_getimport(importState, 'Sa_methane', bounds, output=atm2lnd_inst%forc_pch4_grc, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! The mediator is sending ndep in units if kgN/m2/s - and ctsm uses units of gN/m2/sec
    ! so the following conversion needs to happen

    call state_getimport(importState, 'Faxa_ndep', bounds, output=forc_nhx, ungridded_index=1, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport(importState, 'Faxa_ndep', bounds, output=forc_noy, ungridded_index=2, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    do g = begg,endg
       atm2lnd_inst%forc_ndep_grc(g) = (forc_nhx(g) + forc_noy(g))*1000._r8
    end do

    !--------------------------
    ! Atmosphere co2
    !--------------------------

    ! Set default value to a constant and overwrite for prognostic and diagnostic
    do g = begg,endg
       co2_ppmv_input(g) = co2_ppmv
    end do
    if (co2_type == 'prognostic') then
       fldName = 'Sa_co2prog'
       call ESMF_StateGet(importState, trim(fldname), itemFlag, rc=rc)
       if ( ChkErr(rc,__LINE__,u_FILE_u)) return
       if (itemflag == ESMF_STATEITEM_NOTFOUND .and. co2_type == 'prognostic') then
          call shr_sys_abort( subname//' ERROR: must have Sa_co2prog in import state if co2_type is prognostic' )
       end if
       if (itemflag /= ESMF_STATEITEM_NOTFOUND) then
          call state_getfldptr(importState, trim(fldname), dataPtr, rc=rc )
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          do g = begg,endg
             co2_ppmv_input(g) = dataPtr(g-begg+1)   ! co2 atm prognostic
          end do
       end if
    else if (co2_type == 'diagnostic') then
       fldName = 'Sa_co2diag'
       call ESMF_StateGet(importState, trim(fldname), itemFlag, rc=rc)
       if ( ChkErr(rc,__LINE__,u_FILE_u)) return
       if (itemflag == ESMF_STATEITEM_NOTFOUND .and. co2_type == 'diagnostic') then
          call shr_sys_abort( subname//' ERROR: must have Sa_co2diag in import state if co2_type equal diagnostic')
       end if
       if (itemflag /= ESMF_STATEITEM_NOTFOUND) then
          call state_getfldptr(importState, trim(fldname), dataPtr, rc=rc )
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          do g = begg,endg
             co2_ppmv_input(g) = dataPtr(g-begg+1)   ! co2 atm diagnostic
          end do
       end if
    end if

    ! Note that the following does unit conversions from ppmv to partial pressures (Pa)
    ! Note that forc_pbot is in Pa
    do g = begg,endg
       forc_pbot = atm2lnd_inst%forc_pbot_not_downscaled_grc(g)
       atm2lnd_inst%forc_pco2_grc(g) = co2_ppmv_input(g) * 1.e-6_r8 * forc_pbot
       if (use_c13) then
          atm2lnd_inst%forc_pc13o2_grc(g) = co2_ppmv_input(g) * c13ratio * 1.e-6_r8 * forc_pbot
       end if
    end do

    !--------------------------
    ! Flooding back from river
    !--------------------------

    ! sign convention is positive downward and hierarchy is atm/glc/lnd/rof/ice/ocn.
    ! so water sent from rof to land is negative,
    ! change the sign to indicate addition of water to system.

    if (rof_prognostic) then
       call state_getimport(importState, 'Flrr_flood', bounds, output=wateratm2lndbulk_inst%forc_flood_grc, rc=rc )
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       wateratm2lndbulk_inst%forc_flood_grc(:) = -wateratm2lndbulk_inst%forc_flood_grc(:)
    else
       wateratm2lndbulk_inst%forc_flood_grc(:) = 0._r8
    end if

    if (rof_prognostic) then
       call state_getimport(importState, 'Flrr_volr', bounds, output=wateratm2lndbulk_inst%volr_grc, rc=rc )
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       wateratm2lndbulk_inst%volr_grc(:) = wateratm2lndbulk_inst%volr_grc(:) * (ldomain%area(:) * 1.e6_r8)
    else
       wateratm2lndbulk_inst%volr_grc(:) = 0._r8
    end if

    if (rof_prognostic) then
       call state_getimport(importState, 'Flrr_volrmch', bounds, output=wateratm2lndbulk_inst%volrmch_grc, rc=rc )
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       wateratm2lndbulk_inst%volrmch_grc(:) = wateratm2lndbulk_inst%volrmch_grc(:) * (ldomain%area(:) * 1.e6_r8)
    else
       wateratm2lndbulk_inst%volrmch_grc(:) = 0._r8
    end if

    !--------------------------
    ! Land-ice (glc) fields
    !--------------------------

    if (glc_present) then
       ! We could avoid setting these fields if glc_present is .false., if that would
       ! help with performance. (The downside would be that we wouldn't have these fields
       ! available for diagnostic purposes or to force a later T compset with dlnd.)

       do num = 0,glc_nec
          call state_getimport(importState, 'Sg_ice_covered_elev', bounds, frac_grc(:,num), ungridded_index=num+1, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          call state_getimport(importState, 'Sg_topo_elev'       , bounds, topo_grc(:,num), ungridded_index=num+1, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          call state_getimport(importState, 'Flgg_hflx_elev'     , bounds, hflx_grc(:,num), ungridded_index=num+1, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       end do
       call state_getimport(importState, 'Sg_icemask'               ,  bounds, icemask_grc, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call state_getimport(importState, 'Sg_icemask_coupled_fluxes',  bounds, icemask_coupled_fluxes_grc, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call glc2lnd_inst%set_glc2lnd_fields_nuopc( bounds, glc_present, &
            frac_grc, topo_grc, hflx_grc, icemask_grc, icemask_coupled_fluxes_grc )
    end if

    !--------------------------
    ! Derived quantities for required fields
    ! and corresponding error checks
    !--------------------------

    call derive_quantities(bounds, atm2lnd_inst, wateratm2lndbulk_inst, &
       forc_rainc, forc_rainl, forc_snowc, forc_snowl)

    call check_for_errors(bounds, atm2lnd_inst, wateratm2lndbulk_inst)

  end subroutine import_fields

  !===============================================================================

  subroutine export_fields( gcomp, bounds, glc_present, rof_prognostic, &
       waterlnd2atmbulk_inst, lnd2atm_inst, lnd2glc_inst, rc)

    !-------------------------------
    ! Pack the export state
    !-------------------------------

    use Waterlnd2atmBulkType , only: waterlnd2atmbulk_type

    ! input/output variables
    type(ESMF_GridComp)                         :: gcomp
    type(bounds_type)           , intent(in)    :: bounds       ! bounds
    logical                     , intent(in)    :: glc_present
    logical                     , intent(in)    :: rof_prognostic
    type(waterlnd2atmbulk_type) , intent(inout) :: waterlnd2atmbulk_inst
    type(lnd2atm_type)          , intent(inout) :: lnd2atm_inst ! land to atmosphere exchange data type
    type(lnd2glc_type)          , intent(inout) :: lnd2glc_inst ! land to atmosphere exchange data type
    integer                     , intent(out)   :: rc

    ! local variables
    type(ESMF_State)   :: exportState
    integer            :: i, g, num
    real(r8)           :: array(bounds%begg:bounds%endg)
    character(len=*), parameter :: subname='(lnd_import_export:export_fields)'
    !---------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    ! Get export state
    call NUOPC_ModelGet(gcomp, exportState=exportState, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! -----------------------
    ! output to mediator
    ! -----------------------

    call state_setexport(exportState, 'Sl_lfrin', bounds, input=ldomain%frac, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! -----------------------
    ! output to atm
    ! -----------------------

    call state_setexport(exportState, 'Sl_t', bounds, input=lnd2atm_inst%t_rad_grc, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_setexport(exportState, 'Sl_snowh', bounds, input=waterlnd2atmbulk_inst%h2osno_grc, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_setexport(exportState, 'Sl_avsdr', bounds, input=lnd2atm_inst%albd_grc(bounds%begg:,1), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_setexport(exportState, 'Sl_anidr', bounds, input=lnd2atm_inst%albd_grc(bounds%begg:,2), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_setexport(exportState, 'Sl_avsdf', bounds, input=lnd2atm_inst%albi_grc(bounds%begg:,1), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_setexport(exportState, 'Sl_anidf', bounds, input=lnd2atm_inst%albi_grc(bounds%begg:,2), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_setexport(exportState, 'Sl_tref', bounds, input=lnd2atm_inst%t_ref2m_grc, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_setexport(exportState, 'Sl_qref', bounds, input=waterlnd2atmbulk_inst%q_ref2m_grc, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_setexport(exportState, 'Sl_u10', bounds, input=lnd2atm_inst%u_ref10m_grc, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_setexport(exportState, 'Fall_taux', bounds, input=lnd2atm_inst%taux_grc, minus=.true., rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_setexport(exportState, 'Fall_tauy', bounds, input=lnd2atm_inst%tauy_grc, minus=.true., rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_setexport(exportState, 'Fall_lat', bounds, input=lnd2atm_inst%eflx_lh_tot_grc, minus=.true., rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_setexport(exportState, 'Fall_sen', bounds, input=lnd2atm_inst%eflx_sh_tot_grc, minus=.true., rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_setexport(exportState, 'Fall_lwup', bounds, input=lnd2atm_inst%eflx_lwrad_out_grc, minus=.true., rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_setexport(exportState, 'Fall_evap', bounds, input=waterlnd2atmbulk_inst%qflx_evap_tot_grc, minus=.true., rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_setexport(exportState, 'Fall_swnet', bounds, input=lnd2atm_inst%fsa_grc, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_setexport(exportState, 'Fall_flxdst', bounds, input=lnd2atm_inst%flxdst_grc(:,1), &
         minus=.true., ungridded_index=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_setexport(exportState, 'Fall_flxdst', bounds, input=lnd2atm_inst%flxdst_grc(:,2), &
         minus=.true., ungridded_index=2, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_setexport(exportState, 'Fall_flxdst', bounds, input=lnd2atm_inst%flxdst_grc(:,3), &
         minus=.true., ungridded_index=3, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_setexport(exportState, 'Fall_flxdst', bounds, input=lnd2atm_inst%flxdst_grc(:,4), &
         minus=.true., ungridded_index=4, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_setexport(exportState, 'Fall_methane', bounds, input=lnd2atm_inst%ch4_surf_flux_tot_grc, minus=.true., rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_setexport(exportState, 'Sl_ram1', bounds, input=lnd2atm_inst%ram1_grc, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_setexport(exportState, 'Sl_fv', bounds, input=lnd2atm_inst%fv_grc, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_setexport(exportState, 'Sl_soilw', bounds, &
         input=waterlnd2atmbulk_inst%h2osoi_vol_grc(:,1), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! co2 from land
    if (flds_co2b .or. flds_co2c) then
       call state_setexport(exportState, 'Fall_fco2_lnd', bounds, lnd2atm_inst%net_carbon_exchange_grc, minus=.true., rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    ! dry dep velocities
    do num = 1, drydep_nflds
       call state_setexport(exportState, 'Sl_ddvel', bounds, input=lnd2atm_inst%ddvel_grc(:,num), &
            ungridded_index=num, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end do

    ! MEGAN VOC emis fluxes
    do num = 1, shr_megan_mechcomps_n
       call state_setexport(exportState, 'Fall_voc', bounds, input=lnd2atm_inst%flxvoc_grc(:,num), minus=.true., &
            ungridded_index=num, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end do

    ! fire emis fluxes
    if (emis_nflds > 0) then
       do num = 1, emis_nflds
          call state_setexport(exportState, 'Fall_fire', bounds, input=lnd2atm_inst%fireflx_grc(:,num), minus=.true., &
               ungridded_index=num, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       end do
       call state_setexport(exportState, 'Sl_fztop', bounds, input=lnd2atm_inst%fireztop_grc, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    endif
    ! sign convention is positive downward with hierarchy of atm/glc/lnd/rof/ice/ocn.
    ! i.e. water sent from land to rof is positive

    ! -----------------------
    ! output to river
    ! -----------------------

    ! surface runoff is the sum of qflx_over, qflx_h2osfc_surf
    !    do g = bounds%begg,bounds%endg
    !       array(g) = waterlnd2atmbulk_inst%qflx_rofliq_qsur_grc(g) + waterlnd2atmbulk_inst%qflx_rofliq_h2osfc_grc(g)
    !    end do
    call state_setexport(exportState, 'Flrl_rofsur', bounds, input=waterlnd2atmbulk_inst%qflx_rofliq_qsur_grc, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! subsurface runoff is the sum of qflx_drain and qflx_perched_drain
    do g = bounds%begg,bounds%endg
       array(g) = waterlnd2atmbulk_inst%qflx_rofliq_qsub_grc(g) + waterlnd2atmbulk_inst%qflx_rofliq_drain_perched_grc(g)
    end do
    call state_setexport(exportState, 'Flrl_rofsub', bounds, input=array, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! qgwl sent individually to coupler
    call state_setexport(exportState, 'Flrl_rofgwl', bounds, input=waterlnd2atmbulk_inst%qflx_rofliq_qgwl_grc, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! ice  sent individually to coupler
    call state_setexport(exportState, 'Flrl_rofi', bounds, input=waterlnd2atmbulk_inst%qflx_rofice_grc, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! irrigation flux to be removed from main channel storage (negative)
    call state_setexport(exportState, 'Flrl_irrig', bounds, input=waterlnd2atmbulk_inst%qirrig_grc, minus=.true., rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! -----------------------
    ! output to glc
    ! -----------------------

    ! We could avoid setting these fields if glc_present is .false., if that would
    ! help with performance. (The downside would be that we wouldn't have these fields
    ! available for diagnostic purposes or to force a later T compset with dlnd.)

    do num = 0,glc_nec
       call state_setexport(exportState, 'Sl_tsrf_elev', bounds, input=lnd2glc_inst%tsrf_grc(:,num), &
            ungridded_index=num+1, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call state_setexport(exportState, 'Sl_topo_elev', bounds, input=lnd2glc_inst%topo_grc(:,num), &
            ungridded_index=num+1, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call state_setexport(exportState, 'Flgl_qice_elev', bounds, input=lnd2glc_inst%qice_grc(:,num), &
            ungridded_index=num+1, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end do

  end subroutine export_fields

  !===============================================================================

  subroutine fldlist_add(num, fldlist, stdname, ungridded_lbound, ungridded_ubound)

    ! input/output variables
    integer,                    intent(inout) :: num
    type(fld_list_type),        intent(inout) :: fldlist(:)
    character(len=*),           intent(in)    :: stdname
    integer,          optional, intent(in)    :: ungridded_lbound
    integer,          optional, intent(in)    :: ungridded_ubound

    ! local variables
    integer :: rc
    character(len=*), parameter :: subname='(lnd_import_export:fldlist_add)'
    !-------------------------------------------------------------------------------

    ! Set up a list of field information

    num = num + 1
    if (num > fldsMax) then
       call ESMF_LogWrite(trim(subname)//": ERROR num > fldsMax "//trim(stdname), &
            ESMF_LOGMSG_ERROR, line=__LINE__, file=__FILE__)
       call shr_sys_abort(trim(subname)//": ERROR: num > fldsMax")
    endif
    fldlist(num)%stdname = trim(stdname)

    if (present(ungridded_lbound) .and. present(ungridded_ubound)) then
       fldlist(num)%ungridded_lbound = ungridded_lbound
       fldlist(num)%ungridded_ubound = ungridded_ubound
    end if

  end subroutine fldlist_add

  !===============================================================================

  subroutine fldlist_realize(state, fldList, numflds, flds_scalar_name, flds_scalar_num, mesh, tag, rc)

    use NUOPC , only : NUOPC_IsConnected, NUOPC_Realize
    use ESMF  , only : ESMF_MeshLoc_Element, ESMF_FieldCreate, ESMF_TYPEKIND_R8
    use ESMF  , only : ESMF_MAXSTR, ESMF_Field, ESMF_State, ESMF_Mesh, ESMF_StateRemove
    use ESMF  , only : ESMF_LogFoundError, ESMF_LOGMSG_INFO, ESMF_SUCCESS
    use ESMF  , only : ESMF_LogWrite, ESMF_LOGMSG_ERROR, ESMF_LOGERR_PASSTHRU

    ! input/output variables
    type(ESMF_State)    , intent(inout) :: state
    type(fld_list_type) , intent(in)    :: fldList(:)
    integer             , intent(in)    :: numflds
    character(len=*)    , intent(in)    :: flds_scalar_name
    integer             , intent(in)    :: flds_scalar_num
    character(len=*)    , intent(in)    :: tag
    type(ESMF_Mesh)     , intent(in)    :: mesh
    integer             , intent(inout) :: rc

    ! local variables
    integer                :: n
    type(ESMF_Field)       :: field
    character(len=80)      :: stdname
    character(len=*),parameter  :: subname='(lnd_import_export:fldlist_realize)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    do n = 1, numflds
       stdname = fldList(n)%stdname
       if (NUOPC_IsConnected(state, fieldName=stdname)) then
          if (stdname == trim(flds_scalar_name)) then
             call ESMF_LogWrite(trim(subname)//trim(tag)//" Field = "//trim(stdname)//" is connected on root pe", &
                  ESMF_LOGMSG_INFO)
             ! Create the scalar field
             call SetScalarField(field, flds_scalar_name, flds_scalar_num, rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
          else
             ! Create the field
             if (fldlist(n)%ungridded_lbound > 0 .and. fldlist(n)%ungridded_ubound > 0) then
                field = ESMF_FieldCreate(mesh, ESMF_TYPEKIND_R8, name=stdname, meshloc=ESMF_MESHLOC_ELEMENT, &
                     ungriddedLbound=(/fldlist(n)%ungridded_lbound/), &
                     ungriddedUbound=(/fldlist(n)%ungridded_ubound/), &
                     gridToFieldMap=(/gridToFieldMap/), rc=rc)
                if (ChkErr(rc,__LINE__,u_FILE_u)) return
             else
                field = ESMF_FieldCreate(mesh, ESMF_TYPEKIND_R8, name=stdname, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
                if (ChkErr(rc,__LINE__,u_FILE_u)) return
             end if
             call ESMF_LogWrite(trim(subname)//trim(tag)//" Field = "//trim(stdname)//" is connected using mesh", &
                  ESMF_LOGMSG_INFO)
          endif

          ! NOW call NUOPC_Realize
          call NUOPC_Realize(state, field=field, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       else
          if (stdname /= trim(flds_scalar_name)) then
             call ESMF_LogWrite(subname // trim(tag) // " Field = "// trim(stdname) // " is not connected.", &
                  ESMF_LOGMSG_INFO)
             call ESMF_StateRemove(state, (/stdname/), rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
          end if
       end if
    end do

  contains  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    subroutine SetScalarField(field, flds_scalar_name, flds_scalar_num, rc)
      ! ----------------------------------------------
      ! create a field with scalar data on the root pe
      ! ----------------------------------------------
      use ESMF, only : ESMF_Field, ESMF_DistGrid, ESMF_Grid
      use ESMF, only : ESMF_DistGridCreate, ESMF_GridCreate, ESMF_LogFoundError, ESMF_LOGERR_PASSTHRU
      use ESMF, only : ESMF_FieldCreate, ESMF_GridCreate, ESMF_TYPEKIND_R8

      type(ESMF_Field) , intent(inout) :: field
      character(len=*) , intent(in)    :: flds_scalar_name
      integer          , intent(in)    :: flds_scalar_num
      integer          , intent(inout) :: rc

      ! local variables
      type(ESMF_Distgrid) :: distgrid
      type(ESMF_Grid)     :: grid
      character(len=*), parameter :: subname='(lnd_import_export:SetScalarField)'
      ! ----------------------------------------------

      rc = ESMF_SUCCESS

      ! create a DistGrid with a single index space element, which gets mapped onto DE 0.
      distgrid = ESMF_DistGridCreate(minIndex=(/1/), maxIndex=(/1/), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return

      grid = ESMF_GridCreate(distgrid, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return

      field = ESMF_FieldCreate(name=trim(flds_scalar_name), grid=grid, typekind=ESMF_TYPEKIND_R8, &
           ungriddedLBound=(/1/), ungriddedUBound=(/flds_scalar_num/), gridToFieldMap=(/2/), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return

    end subroutine SetScalarField

  end subroutine fldlist_realize

  !===============================================================================

  subroutine state_getimport(state, fldname, bounds, output, ungridded_index, rc)

    ! ----------------------------------------------
    ! Map import state field to output array
    ! ----------------------------------------------

    ! input/output variables
    type(ESMF_State)    , intent(in)    :: state
    type(bounds_type)   , intent(in)    :: bounds
    character(len=*)    , intent(in)    :: fldname
    real(r8)            , intent(out)   :: output(bounds%begg:bounds%endg)
    integer, optional   , intent(in)    :: ungridded_index
    integer             , intent(out)   :: rc

    ! local variables
    integer                     :: g, i,n
    real(R8), pointer           :: fldptr1d(:)
    real(R8), pointer           :: fldptr2d(:,:)
    type(ESMF_StateItem_Flag)   :: itemFlag
    character(len=cs)           :: cvalue
    character(len=*), parameter :: subname='(lnd_import_export:state_getimport)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    ! Determine if field with name fldname exists in state
    call ESMF_StateGet(state, trim(fldname), itemFlag, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! if field exists then create output array - else do nothing
    if (itemflag /= ESMF_STATEITEM_NOTFOUND) then

       ! get field pointer
       if (present(ungridded_index)) then
          write(cvalue,*) ungridded_index
          call ESMF_LogWrite(trim(subname)//": getting import for "//trim(fldname)//" index "//trim(cvalue), &
               ESMF_LOGMSG_INFO)
          call state_getfldptr(state, trim(fldname), fldptr2d=fldptr2d, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       else
          call ESMF_LogWrite(trim(subname)//": getting import for "//trim(fldname),ESMF_LOGMSG_INFO)
          call state_getfldptr(state, trim(fldname), fldptr1d=fldptr1d, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       end if

       ! determine output array
       if (present(ungridded_index)) then
          if (gridToFieldMap == 1) then
             do g = bounds%begg, bounds%endg
                n = g - bounds%begg + 1
                output(g) = fldptr2d(n,ungridded_index)
             end do
          else if (gridToFieldMap == 2) then
             do g = bounds%begg, bounds%endg
                n = g - bounds%begg + 1
                output(g) = fldptr2d(ungridded_index,n)
             end do
          end if
       else
          do g = bounds%begg, bounds%endg
             n = g - bounds%begg + 1
             output(g) = fldptr1d(n)
          end do
       end if

       ! write debug output if appropriate
       if (masterproc .and. debug > 0 .and. get_nstep() < 48) then
          do g = bounds%begg,bounds%endg
             i = 1 + g - bounds%begg
             write(iulog,F01)'import: nstep, n, '//trim(fldname)//' = ',get_nstep(),i,output(g)
          end do
       end if

       ! check for nans
       call check_for_nans(output, trim(fldname), bounds%begg)
    end if

  end subroutine state_getimport

  !===============================================================================

  subroutine state_setexport(state, fldname, bounds, input, minus, ungridded_index, rc)

    ! ----------------------------------------------
    ! Map input array to export state field
    ! ----------------------------------------------

    use shr_const_mod, only : fillvalue=>SHR_CONST_SPVAL

    ! input/output variables
    type(ESMF_State)    , intent(inout) :: state
    type(bounds_type)   , intent(in)    :: bounds
    character(len=*)    , intent(in)    :: fldname
    real(r8)            , intent(in)    :: input(bounds%begg:bounds%endg)
    logical, optional   , intent(in)    :: minus
    integer, optional   , intent(in)    :: ungridded_index
    integer             , intent(out)   :: rc

    ! local variables
    logical                     :: l_minus  ! local version of minus
    integer                     :: g, i, n
    real(R8), pointer           :: fldptr1d(:)
    real(R8), pointer           :: fldptr2d(:,:)
    character(len=cs)           :: cvalue
    type(ESMF_StateItem_Flag)   :: itemFlag
    character(len=*), parameter :: subname='(lnd_import_export:state_setexport)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    l_minus = .false.
    if (present(minus)) then
       l_minus = minus
    end if

    ! Determine if field with name fldname exists in state
    call ESMF_StateGet(state, trim(fldname), itemFlag, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! if field exists then create output array - else do nothing
    if (itemflag /= ESMF_STATEITEM_NOTFOUND) then

       ! get field pointer
       if (present(ungridded_index)) then
          write(cvalue,*) ungridded_index
          call ESMF_LogWrite(trim(subname)//": setting export for "//trim(fldname)//" index "//trim(cvalue), &
               ESMF_LOGMSG_INFO)
          call state_getfldptr(state, trim(fldname), fldptr2d=fldptr2d, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       else
          call ESMF_LogWrite(trim(subname)//": setting export for "//trim(fldname), ESMF_LOGMSG_INFO)
          call state_getfldptr(state, trim(fldname), fldptr1d=fldptr1d, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       end if

       ! TODO: if fillvalue = shr_const_spval the snowhl sent to the atm will have the spval over some points
       ! rather than 0 - this is very odd and needs to be understood
       !fldptr1d(:) = fillvalue

       ! determine output array
       if (present(ungridded_index)) then
          fldptr2d(ungridded_index,:) = 0._r8
          !fldptr2d(ungridded_index,:) = fillvalue
          do g = bounds%begg, bounds%endg
             n = g - bounds%begg + 1
             fldptr2d(ungridded_index,n) = input(g)
          end do
          if (l_minus) then
             fldptr2d(ungridded_index,:) = -fldptr2d(ungridded_index,:)
          end if
       else
          fldptr1d(:) = 0._r8
          !fldptr1d(:) = fillvalue
          do g = bounds%begg, bounds%endg
             n = g - bounds%begg + 1
             fldptr1d(n) = input(g)
          end do
          if (l_minus) then
             fldptr1d(:) = -fldptr1d(:)
          end if
       end if

       ! write debug output if appropriate
       if (masterproc .and. debug > 0 .and. get_nstep() < 48) then
          do g = bounds%begg,bounds%endg
             i = 1 + g - bounds%begg
             write(iulog,F01)'export: nstep, n, '//trim(fldname)//' = ',get_nstep(),i,input(g)
          end do
       end if

       ! check for nans
       call check_for_nans(input, trim(fldname), bounds%begg)
    end if

  end subroutine state_setexport

  !===============================================================================

  subroutine state_getfldptr(State, fldname, fldptr1d, fldptr2d, rc)

    ! ----------------------------------------------
    ! Get pointer to a state field
    ! ----------------------------------------------

    use ESMF , only : ESMF_State, ESMF_Field, ESMF_Mesh, ESMF_FieldStatus_Flag
    use ESMF , only : ESMF_StateGet, ESMF_FieldGet, ESMF_MeshGet
    use ESMF , only : ESMF_FIELDSTATUS_COMPLETE, ESMF_FAILURE

    ! input/output variables
    type(ESMF_State),             intent(in)    :: State
    character(len=*),             intent(in)    :: fldname
    real(R8), pointer, optional , intent(out)   :: fldptr1d(:)
    real(R8), pointer, optional , intent(out)   :: fldptr2d(:,:)
    integer,                      intent(out)   :: rc

    ! local variables
    type(ESMF_FieldStatus_Flag) :: status
    type(ESMF_Field)            :: lfield
    type(ESMF_Mesh)             :: lmesh
    integer                     :: nnodes, nelements
    character(len=*), parameter :: subname='(lnd_import_export:state_getfldptr)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    call ESMF_StateGet(State, itemName=trim(fldname), field=lfield, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_FieldGet(lfield, status=status, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (status /= ESMF_FIELDSTATUS_COMPLETE) then
       call ESMF_LogWrite(trim(subname)//": ERROR data not allocated ", ESMF_LOGMSG_INFO, rc=rc)
       rc = ESMF_FAILURE
       return
    else
       call ESMF_FieldGet(lfield, mesh=lmesh, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call ESMF_MeshGet(lmesh, numOwnedNodes=nnodes, numOwnedElements=nelements, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       if (nnodes == 0 .and. nelements == 0) then
          call ESMF_LogWrite(trim(subname)//": no local nodes or elements ", ESMF_LOGMSG_INFO)
          rc = ESMF_FAILURE
          return
       end if

       if (present(fldptr1d)) then
          call ESMF_FieldGet(lfield, farrayPtr=fldptr1d, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       else if (present(fldptr2d)) then
          call ESMF_FieldGet(lfield, farrayPtr=fldptr2d, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       else
          call shr_sys_abort("either fldptr1d or fldptr2d must be an input argument")
       end if
    endif  ! status

  end subroutine state_getfldptr

end module lnd_import_export
