module lnd_import_export

  use ESMF                  , only : ESMF_GridComp, ESMF_State, ESMF_Mesh, ESMF_StateGet
  use ESMF                  , only : ESMF_KIND_R8, ESMF_SUCCESS, ESMF_MAXSTR, ESMF_LOGMSG_INFO
  use ESMF                  , only : ESMF_LogWrite, ESMF_LOGMSG_ERROR, ESMF_LogFoundError
  use ESMF                  , only : ESMF_STATEITEM_NOTFOUND, ESMF_StateItem_Flag
  use ESMF                  , only : operator(/=), operator(==)
  use NUOPC                 , only : NUOPC_CompAttributeGet, NUOPC_Advertise, NUOPC_IsConnected
  use NUOPC_Model           , only : NUOPC_ModelGet
  use shr_kind_mod          , only : r8 => shr_kind_r8, cx=>shr_kind_cx, cxx=>shr_kind_cxx
  use shr_infnan_mod        , only : isnan => shr_infnan_isnan
  use shr_string_mod        , only : shr_string_listGetName, shr_string_listGetNum
  use shr_sys_mod           , only : shr_sys_abort
  use shr_nuopc_methods_mod , only : shr_nuopc_methods_chkerr
  use shr_nuopc_scalars_mod , only : flds_scalar_name, flds_scalar_num
  use clm_varctl            , only : iulog
  use clm_time_manager      , only : get_nstep
  use decompmod             , only : bounds_type
  use lnd2atmType           , only : lnd2atm_type
  use lnd2glcMod            , only : lnd2glc_type
  use atm2lndType           , only : atm2lnd_type
  use glc2lndMod            , only : glc2lnd_type
  use domainMod             , only : ldomain
  use spmdMod               , only : masterproc, mpicom, comp_id
  use seq_drydep_mod        , only : seq_drydep_readnl, seq_drydep_init, n_drydep
  use shr_megan_mod         , only : shr_megan_readnl, shr_megan_mechcomps_n
  use shr_fire_emis_mod     , only : shr_fire_emis_readnl, shr_fire_emis_mechcomps_n, shr_fire_emis_ztop_token
  use shr_carma_mod         , only : shr_carma_readnl
  use shr_ndep_mod          , only : shr_ndep_readnl
  use glc_elevclass_mod     , only : glc_elevclass_init, glc_get_num_elevation_classes, glc_elevclass_as_string

  implicit none
  private ! except

  public  :: advertise_fields
  public  :: realize_fields
  public  :: import_fields
  public  :: export_fields

  private :: fldlist_add
  private :: fldlist_realize
  private :: state_getfldptr
  private :: check_for_nans

  type fld_list_type
     character(len=128) :: stdname
  end type fld_list_type

  integer, parameter     :: fldsMax = 100
  integer                :: fldsToLnd_num = 0
  integer                :: fldsFrLnd_num = 0
  type (fld_list_type)   :: fldsToLnd(fldsMax)
  type (fld_list_type)   :: fldsFrLnd(fldsMax)

  character(len=cxx)     :: drydep_fields       ! List of dry-deposition fields
  character(len=cxx)     :: megan_voc_fields    ! List of MEGAN VOC emission fields
  character(len=cxx)     :: fire_emis_fields    ! List of fire emission fields
  character(len=cx)      :: carma_fields        ! List of CARMA fields from lnd->atm
  character(len=cx)      :: ndep_fields         ! List of nitrogen deposition fields from atm->lnd/ocn

  integer                :: glc_nec
  integer, parameter     :: dbug = 10
  integer, parameter     :: debug_import = 1 ! internal debug level
  integer, parameter     :: debug_export = 1 ! internal debug level
  character(*),parameter :: F01 = "('(lnd_import_export) ',a,i5,2x,i5,2x,d21.14)"
  character(*),parameter :: u_FILE_u = &
       __FILE__

!===============================================================================
contains
!===============================================================================

  subroutine advertise_fields(gcomp, rc)
    use clm_varctl, only : ndep_from_cpl

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_State)       :: importState
    type(ESMF_State)       :: exportState
    character(ESMF_MAXSTR) :: stdname
    character(ESMF_MAXSTR) :: cvalue
    character(len=2)       :: nec_str
    integer                :: dbrc
    integer                :: n, num
    logical                :: flds_co2a      ! use case
    logical                :: flds_co2b      ! use case
    logical                :: flds_co2c      ! use case
    logical                :: add_ndep_fields
    character(len=128)     :: fldname
    character(len=*), parameter :: subname='(lnd_import_export:advertise_fields)'
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    call NUOPC_ModelGet(gcomp, importState=importState, exportState=exportState, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !--------------------------------
    ! determine necessary toggles for below
    !--------------------------------

    call NUOPC_CompAttributeGet(gcomp, name='flds_co2a', value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) flds_co2a
    call ESMF_LogWrite('flds_co2a = '// trim(cvalue), ESMF_LOGMSG_INFO, rc=dbrc)

    call NUOPC_CompAttributeGet(gcomp, name='flds_co2b', value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) flds_co2b
    call ESMF_LogWrite('flds_co2b = '// trim(cvalue), ESMF_LOGMSG_INFO, rc=dbrc)

    call NUOPC_CompAttributeGet(gcomp, name='flds_co2c', value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) flds_co2c
    call ESMF_LogWrite('flds_co2c = '// trim(cvalue), ESMF_LOGMSG_INFO, rc=dbrc)

    ! Determine  number of elevation classes
    call NUOPC_CompAttributeGet(gcomp, name='glc_nec', value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) glc_nec
    call ESMF_LogWrite('glc_nec = '// trim(cvalue), ESMF_LOGMSG_INFO, rc=dbrc)
    if (glc_nec < 1) then
       call shr_sys_abort('ERROR: In CLM4.5 and later, glc_nec must be at least 1.')
    end if

    ! Initialize glc_elevclass module
    call glc_elevclass_init(glc_nec)

    !--------------------------------
    ! Advertise export fields
    !--------------------------------

    call fldlist_add(fldsFrLnd_num, fldsFrlnd, trim(flds_scalar_name))

    ! land states
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

    ! fluxes to river
    call fldlist_add(fldsFrLnd_num, fldsFrlnd, 'Flrl_rofsur'   )
    call fldlist_add(fldsFrLnd_num, fldsFrlnd, 'Flrl_rofgwl'   )
    call fldlist_add(fldsFrLnd_num, fldsFrlnd, 'Flrl_rofsub'   )
    call fldlist_add(fldsFrLnd_num, fldsFrlnd, 'Flrl_rofi'     )
    call fldlist_add(fldsFrLnd_num, fldsFrlnd, 'Flrl_irrig'    )

    ! fluxes to atm
    call fldlist_add(fldsFrLnd_num, fldsFrlnd, 'Fall_taux'     )
    call fldlist_add(fldsFrLnd_num, fldsFrlnd, 'Fall_tauy'     )
    call fldlist_add(fldsFrLnd_num, fldsFrlnd, 'Fall_lat'      )
    call fldlist_add(fldsFrLnd_num, fldsFrlnd, 'Fall_sen'      )
    call fldlist_add(fldsFrLnd_num, fldsFrlnd, 'Fall_lwup'     )
    call fldlist_add(fldsFrLnd_num, fldsFrlnd, 'Fall_evap'     )
    call fldlist_add(fldsFrLnd_num, fldsFrlnd, 'Fall_swnet'    )
    call fldlist_add(fldsFrLnd_num, fldsFrlnd, 'Fall_flxdst1'  )
    call fldlist_add(fldsFrLnd_num, fldsFrlnd, 'Fall_flxdst2'  )
    call fldlist_add(fldsFrLnd_num, fldsFrlnd, 'Fall_flxdst3'  )
    call fldlist_add(fldsFrLnd_num, fldsFrlnd, 'Fall_flxdst4'  )
    call fldlist_add(fldsFrLnd_num, fldsFrlnd, 'Fall_fco2_lnd' )
   !call fldlist_add(fldsFrLnd_num, fldsFrlnd, 'Fall_methane'  )

    ! Dry Deposition fluxes from land - ALSO initialize drydep here
    call seq_drydep_readnl("drv_flds_in", mpicom, masterproc, drydep_fields)
    do n = 1, n_drydep
       call shr_string_listGetName(drydep_fields, n, fldname)
       call fldlist_add(fldsFrLnd_num, fldsFrLnd, trim(fldname))
    end do
    call seq_drydep_init( )

    ! MEGAN VOC emissions fluxes from land
    call shr_megan_readnl('drv_flds_in', mpicom, masterproc, megan_voc_fields)
    do n = 1, shr_megan_mechcomps_n
       call shr_string_listGetName(megan_voc_fields, n, fldname)
       call fldlist_add(fldsFrLnd_num, fldsFrLnd, trim(fldname))
    end do

    ! Fire emissions fluxes from land
    call shr_fire_emis_readnl('drv_flds_in', mpicom, masterproc, fire_emis_fields)
    if (shr_fire_emis_mechcomps_n > 0) then
       do n = 1, shr_fire_emis_mechcomps_n
          call shr_string_listGetName(fire_emis_fields, n, fldname)
          call fldlist_add(fldsFrLnd_num, fldsFrLnd, trim(fldname))
       end do
       call fldlist_add(fldsFrLnd_num, fldsFrLnd, trim(shr_fire_emis_ztop_token))
    end if

    ! CARMA volumetric soil water from land
    ! TODO: is the following correct - the CARMA field exchange is very confusing in mct
    call shr_carma_readnl('drv_flds_in', mpicom, masterproc, carma_fields)
    if (carma_fields /= ' ') then
       call fldlist_add(fldsFrLnd_num, fldsFrlnd, 'Sl_soilw') ! optional for carma
    end if

    ! lnd->glc states from land all lnd->glc elevation classes (1:glc_nec) plus bare land (index 0).
    do num = 0,glc_nec
       nec_str = glc_elevclass_as_string(num)
       fldname = 'Sl_tsrf' // nec_str
       call fldlist_add(fldsFrLnd_num, fldsFrLnd, trim(fldname))
       fldname = 'Sl_topo' // nec_str
       call fldlist_add(fldsFrLnd_num, fldsFrLnd, trim(fldname))
       fldname = 'Flgl_qice' // nec_str
       call fldlist_add(fldsFrLnd_num, fldsFrLnd, trim(fldname))
    end do

    ! Now advertise above export fields
    do n = 1,fldsFrLnd_num
       call NUOPC_Advertise(exportState, standardName=fldsFrLnd(n)%stdname, &
            TransferOfferGeomObject='will provide', rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    enddo

    !--------------------------------
    ! Advertise import fields
    !--------------------------------

    call fldlist_add(fldsToLnd_num, fldsToLnd, trim(flds_scalar_name))

    ! from atm
    call fldlist_add(fldsToLnd_num, fldsToLnd, 'Sa_z'         )
    call fldlist_add(fldsToLnd_num, fldsToLnd, 'Sa_topo'      )
    call fldlist_add(fldsToLnd_num, fldsToLnd, 'Sa_u'         )
    call fldlist_add(fldsToLnd_num, fldsToLnd, 'Sa_v'         )
    call fldlist_add(fldsToLnd_num, fldsToLnd, 'Sa_ptem'      )
    call fldlist_add(fldsToLnd_num, fldsToLnd, 'Sa_pbot'      )
    call fldlist_add(fldsToLnd_num, fldsToLnd, 'Sa_tbot'      )
    call fldlist_add(fldsToLnd_num, fldsToLnd, 'Sa_shum'      )
   !call fldlist_add(fldsToLnd_num, fldsToLnd, 'Sa_methane'   )

    call fldlist_add(fldsToLnd_num, fldsToLnd, 'Flrr_volr'    )
    call fldlist_add(fldsToLnd_num, fldsToLnd, 'Flrr_volrmch' )
    call fldlist_add(fldsToLnd_num, fldsToLnd, 'Faxa_lwdn'    )
    call fldlist_add(fldsToLnd_num, fldsToLnd, 'Faxa_rainc'   )
    call fldlist_add(fldsToLnd_num, fldsToLnd, 'Faxa_rainl'   )
    call fldlist_add(fldsToLnd_num, fldsToLnd, 'Faxa_snowc'   )
    call fldlist_add(fldsToLnd_num, fldsToLnd, 'Faxa_snowl'   )
    call fldlist_add(fldsToLnd_num, fldsToLnd, 'Faxa_swndr'   )
    call fldlist_add(fldsToLnd_num, fldsToLnd, 'Faxa_swvdr'   )
    call fldlist_add(fldsToLnd_num, fldsToLnd, 'Faxa_swndf'   )
    call fldlist_add(fldsToLnd_num, fldsToLnd, 'Faxa_swvdf'   )
    call fldlist_add(fldsToLnd_num, fldsToLnd, 'Faxa_bcphidry')
    call fldlist_add(fldsToLnd_num, fldsToLnd, 'Faxa_bcphodry')
    call fldlist_add(fldsToLnd_num, fldsToLnd, 'Faxa_bcphiwet')
    call fldlist_add(fldsToLnd_num, fldsToLnd, 'Faxa_ocphidry')
    call fldlist_add(fldsToLnd_num, fldsToLnd, 'Faxa_ocphodry')
    call fldlist_add(fldsToLnd_num, fldsToLnd, 'Faxa_ocphiwet')
    call fldlist_add(fldsToLnd_num, fldsToLnd, 'Faxa_dstdry1' )
    call fldlist_add(fldsToLnd_num, fldsToLnd, 'Faxa_dstdry2' )
    call fldlist_add(fldsToLnd_num, fldsToLnd, 'Faxa_dstdry3' )
    call fldlist_add(fldsToLnd_num, fldsToLnd, 'Faxa_dstdry4' )
    call fldlist_add(fldsToLnd_num, fldsToLnd, 'Faxa_dstwet1' )
    call fldlist_add(fldsToLnd_num, fldsToLnd, 'Faxa_dstwet2' )
    call fldlist_add(fldsToLnd_num, fldsToLnd, 'Faxa_dstwet3' )
    call fldlist_add(fldsToLnd_num, fldsToLnd, 'Faxa_dstwet4' )

    ! from atm - co2 exchange scenarios
    if (flds_co2a .or. flds_co2b .or. flds_co2c) then
       call fldlist_add(fldsToLnd_num, fldsToLnd, 'Sa_co2prog')
       call fldlist_add(fldsToLnd_num, fldsToLnd, 'Sa_co2diag')
    end if

    ! from atm - nitrogen deposition
    call shr_ndep_readnl("drv_flds_in", mpicom, masterproc, ndep_fields, add_ndep_fields)
    if (add_ndep_fields) then
       do n = 1, shr_string_listGetNum(ndep_fields)
          call  shr_string_listGetName(ndep_fields, n, fldname)
          call fldlist_add(fldsToLnd_num, fldsToLnd, trim(fldname))
       end do
       ! This sets a variable in clm_varctl
       ndep_from_cpl = .true.
    end if

    ! from river
    call fldlist_add(fldsToLnd_num, fldsToLnd, 'Flrr_flood')

    ! from glc
    call fldlist_add(fldsToLnd_num, fldsToLnd, 'Sg_icemask'               )
    call fldlist_add(fldsToLnd_num, fldsToLnd, 'Sg_icemask_coupled_fluxes')
    do num = 0,glc_nec
       ! create coupling fields for all glc->lnd elevation classes (1:glc_nec) plus bare land (index 0).
       nec_str = glc_elevclass_as_string(num)
       fldname = 'Sg_ice_covered' // nec_str
       call fldlist_add(fldsToLnd_num, fldsToLnd, trim(fldname))
       fldname = 'Sg_topo' // nec_str
       call fldlist_add(fldsToLnd_num, fldsToLnd, trim(fldname))
       fldname = 'Flgg_hflx' // nec_str
       call fldlist_add(fldsToLnd_num, fldsToLnd, trim(fldname))
    end do

    ! Now advertise import fields
    do n = 1,fldsToLnd_num
       call NUOPC_Advertise(importState, standardName=fldsToLnd(n)%stdname, &
            TransferOfferGeomObject='will provide', rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    enddo

  end subroutine advertise_fields

  !===============================================================================

  subroutine realize_fields(gcomp, Emesh, rc)
    type(ESMF_GridComp)  :: gcomp
    type(ESMF_Mesh)      :: Emesh
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_State)     :: importState
    type(ESMF_State)     :: exportState
    character(len=*), parameter :: subname='(lnd_import_export:realize_fields)'
    !---------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    call NUOPC_ModelGet(gcomp, importState=importState, exportState=exportState, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call fldlist_realize( &
         state=ExportState, &
         fldList=fldsFrLnd, &
         numflds=fldsFrLnd_num, &
         flds_scalar_name=flds_scalar_name, &
         flds_scalar_num=flds_scalar_num, &
         tag=subname//':clmExport',&
         mesh=Emesh, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call fldlist_realize( &
         state=importState, &
         fldList=fldsToLnd, &
         numflds=fldsToLnd_num, &
         flds_scalar_name=flds_scalar_name, &
         flds_scalar_num=flds_scalar_num, &
         tag=subname//':clmImport',&
         mesh=Emesh, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

  end subroutine realize_fields

  !===============================================================================

  subroutine import_fields( gcomp, bounds, glc_present, atm2lnd_inst, glc2lnd_inst, rc)

    !---------------------------------------------------------------------------
    ! Convert the input data from the mediator to the land model
    !---------------------------------------------------------------------------

    use clm_varctl      , only: co2_type, co2_ppmv, use_c13, ndep_from_cpl
    use clm_varcon      , only: rair, o2_molar_const, c13ratio
    use shr_const_mod   , only: SHR_CONST_TKFRZ

    ! input/output variabes
    type(ESMF_GridComp)                :: gcomp
    type(bounds_type)  , intent(in)    :: bounds       ! bounds
    logical            , intent(in)    :: glc_present  ! .true. => running with a non-stub GLC model
    type(atm2lnd_type) , intent(inout) :: atm2lnd_inst ! clm internal input data type
    type(glc2lnd_type) , intent(inout) :: glc2lnd_inst ! clm internal input data type
    integer            , intent(out)   :: rc

    ! local variables
    type(ESMF_State)          :: importState
    type(ESMF_StateItem_Flag) :: itemFlag
    type(ESMF_StateItem_Flag) :: itemFlag1
    type(ESMF_StateItem_Flag) :: itemFlag2
    real(r8), pointer         :: dataPtr(:)
    real(r8), pointer         :: dataPtr1(:)
    real(r8), pointer         :: dataPtr2(:)
    character(len=128)        :: fldname
    character(len=128)        :: fldname1
    character(len=128)        :: fldname2
    integer                   :: dbrc
    character(len=2)          :: nec_str
    integer                   :: num
    integer                   :: begg, endg                             ! bounds
    integer                   :: g,i,k                                  ! indices
    real(r8)                  :: e                                      ! vapor pressure (Pa)
    real(r8)                  :: qsat                                   ! saturation specific humidity (kg/kg)
    real(r8)                  :: co2_ppmv_diag(bounds%begg:bounds%endg) ! temporary
    real(r8)                  :: co2_ppmv_prog(bounds%begg:bounds%endg) ! temporary
    real(r8)                  :: co2_ppmv_val                           ! temporary
    real(r8)                  :: esatw                                  ! saturation vapor pressure over water (Pa)
    real(r8)                  :: esati                                  ! saturation vapor pressure over ice (Pa)
    real(r8)                  :: a0,a1,a2,a3,a4,a5,a6                   ! coefficients for esat over water
    real(r8)                  :: b0,b1,b2,b3,b4,b5,b6                   ! coefficients for esat over ice
    real(r8)                  :: tdc, t                                 ! Kelvins to Celcius function and its input
    real(r8)                  :: forc_t                                 ! atmospheric temperature (Kelvin)
    real(r8)                  :: forc_q                                 ! atmospheric specific humidity (kg/kg)
    real(r8)                  :: forc_pbot                              ! atmospheric pressure (Pa)
    real(r8)                  :: forc_rainc(bounds%begg:bounds%endg)    ! rainxy Atm flux mm/s
    real(r8)                  :: forc_rainl(bounds%begg:bounds%endg)    ! rainxy Atm flux mm/s
    real(r8)                  :: forc_snowc(bounds%begg:bounds%endg)    ! snowfxy Atm flux  mm/s
    real(r8)                  :: forc_snowl(bounds%begg:bounds%endg)    ! snowfxl Atm flux  mm/s
    real(r8)                  :: frac_grc(bounds%begg:bounds%endg, 0:glc_nec)
    real(r8)                  :: topo_grc(bounds%begg:bounds%endg, 0:glc_nec)
    real(r8)                  :: hflx_grc(bounds%begg:bounds%endg, 0:glc_nec)
    real(r8)                  :: icemask_grc(bounds%begg:bounds%endg)
    real(r8)                  :: icemask_coupled_fluxes_grc(bounds%begg:bounds%endg)
    character(len=*), parameter :: subname='(lnd_import_export:import_fields)'

    ! Constants to compute vapor pressure
    parameter (a0=6.107799961_r8    , a1=4.436518521e-01_r8, &
         a2=1.428945805e-02_r8, a3=2.650648471e-04_r8, &
         a4=3.031240396e-06_r8, a5=2.034080948e-08_r8, &
         a6=6.136820929e-11_r8)

    parameter (b0=6.109177956_r8    , b1=5.034698970e-01_r8, &
         b2=1.886013408e-02_r8, b3=4.176223716e-04_r8, &
         b4=5.824720280e-06_r8, b5=4.838803174e-08_r8, &
         b6=1.838826904e-10_r8)

    ! function declarations
    tdc(t) = min( 50._r8, max(-50._r8,(t-SHR_CONST_TKFRZ)) )
    esatw(t) = 100._r8*(a0+t*(a1+t*(a2+t*(a3+t*(a4+t*(a5+t*a6))))))
    esati(t) = 100._r8*(b0+t*(b1+t*(b2+t*(b3+t*(b4+t*(b5+t*b6))))))
    !---------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    if (dbug > 5) call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO, rc=dbrc)

    ! Get import state
    call NUOPC_ModelGet(gcomp, importState=importState, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

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

    call state_getfldptr(importState, 'Sa_z', dataPtr, begg, rc=rc )
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    do g = begg,endg
       atm2lnd_inst%forc_hgt_grc(g) = dataPtr(g-begg+1)
    end do

    call state_getfldptr(importState, 'Sa_topo', dataPtr, begg, rc=rc )
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    do g = begg,endg
       atm2lnd_inst%forc_topo_grc(g) = dataPtr(g-begg+1)
    end do

    call state_getfldptr(importState, 'Sa_u', dataPtr, begg, rc=rc )
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    do g = begg,endg
       atm2lnd_inst%forc_u_grc(g) = dataPtr(g-begg+1)
    end do

    call state_getfldptr(importState, 'Sa_v', dataPtr, begg, rc=rc )
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    do g = begg,endg
       atm2lnd_inst%forc_v_grc(g) = dataPtr(g-begg+1)
    end do

    call state_getfldptr(importState, 'Sa_ptem', dataPtr, begg, rc=rc )
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    do g = begg,endg
       atm2lnd_inst%forc_th_not_downscaled_grc(g) = dataPtr(g-begg+1)
    end do

    call state_getfldptr(importState, 'Sa_shum', dataPtr, begg, rc=rc )
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    do g = begg,endg
       atm2lnd_inst%forc_q_not_downscaled_grc(g) = dataPtr(g-begg+1)
    end do

    call state_getfldptr(importState, 'Sa_pbot', dataPtr, begg, rc=rc )
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    do g = begg,endg
       atm2lnd_inst%forc_pbot_not_downscaled_grc(g) = dataPtr(g-begg+1)
    end do

    call state_getfldptr(importState, 'Sa_tbot', dataPtr, begg, rc=rc )
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    do g = begg,endg
       atm2lnd_inst%forc_t_not_downscaled_grc(g) = dataPtr(g-begg+1)
    end do

    call state_getfldptr(importState, 'Faxa_rainc', dataPtr, begg, rc=rc )
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    do g = begg,endg
       forc_rainc(g) = dataPtr(g-begg+1)
    end do

    call state_getfldptr(importState, 'Faxa_snowl', dataPtr, begg, rc=rc )
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    do g = begg,endg
       forc_snowl(g) = dataPtr(g-begg+1)
    end do

    call state_getfldptr(importState, 'Faxa_snowc', dataPtr, begg, rc=rc )
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    do g = begg,endg
       forc_snowc(g) = dataPtr(g-begg+1)
    end do

    call state_getfldptr(importState, 'Faxa_rainl', dataPtr, begg, rc=rc )
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    do g = begg,endg
       forc_rainl(g) = dataPtr(g-begg+1)
    end do

    call state_getfldptr(importState, 'Faxa_lwdn', dataPtr, begg, rc=rc )
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    do g = begg,endg
       atm2lnd_inst%forc_lwrad_not_downscaled_grc(g) = dataPtr(g-begg+1)
    end do

    call state_getfldptr(importState, 'Faxa_swndr', dataPtr, begg, rc=rc )
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    do g = begg,endg
       atm2lnd_inst%forc_solad_grc(g,2) = dataPtr(g-begg+1)
    end do

    call state_getfldptr(importState, 'Faxa_swvdr', dataPtr, begg, rc=rc )
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    do g = begg,endg
       atm2lnd_inst%forc_solad_grc(g,1) = dataPtr(g-begg+1)
    end do

    call state_getfldptr(importState, 'Faxa_swndf', dataPtr, begg, rc=rc )
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    do g = begg,endg
       atm2lnd_inst%forc_solai_grc(g,2) = dataPtr(g-begg+1)
    end do

    call state_getfldptr(importState, 'Faxa_swvdf', dataPtr, begg, rc=rc )
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    do g = begg,endg
       atm2lnd_inst%forc_solai_grc(g,1) = dataPtr(g-begg+1)
    end do

    ! Atmosphere prognostic/prescribed aerosol fields

    call state_getfldptr(importState, 'Faxa_bcphidry', dataPtr, begg, rc=rc )
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    do g = begg,endg
       atm2lnd_inst%forc_aer_grc(g,1) = dataPtr(g-begg+1)
    end do

    call state_getfldptr(importState, 'Faxa_bcphodry', dataPtr, begg, rc=rc )
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    do g = begg,endg
       atm2lnd_inst%forc_aer_grc(g,2) = dataPtr(g-begg+1)
    end do

    call state_getfldptr(importState, 'Faxa_bcphiwet', dataPtr, begg, rc=rc )
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    do g = begg,endg
       atm2lnd_inst%forc_aer_grc(g,3) = dataPtr(g-begg+1)
    end do

    call state_getfldptr(importState, 'Faxa_ocphidry', dataPtr, begg, rc=rc )
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    do g = begg,endg
       atm2lnd_inst%forc_aer_grc(g,4) = dataPtr(g-begg+1)
    end do

    call state_getfldptr(importState, 'Faxa_ocphodry', dataPtr, begg, rc=rc )
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    do g = begg,endg
       atm2lnd_inst%forc_aer_grc(g,5) = dataPtr(g-begg+1)
    end do

    call state_getfldptr(importState, 'Faxa_ocphiwet', dataPtr, begg, rc=rc )
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    do g = begg,endg
       atm2lnd_inst%forc_aer_grc(g,6) = dataPtr(g-begg+1)
    end do

    call state_getfldptr(importState, 'Faxa_dstwet1', dataPtr, begg, rc=rc )
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    do g = begg,endg
       atm2lnd_inst%forc_aer_grc(g,7) = dataPtr(g-begg+1)
    end do

    call state_getfldptr(importState, 'Faxa_dstdry1', dataPtr, begg, rc=rc )
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    do g = begg,endg
       atm2lnd_inst%forc_aer_grc(g,8) = dataPtr(g-begg+1)
    end do

    call state_getfldptr(importState, 'Faxa_dstwet2', dataPtr, begg, rc=rc )
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    do g = begg,endg
       atm2lnd_inst%forc_aer_grc(g,9) = dataPtr(g-begg+1)
    end do

    call state_getfldptr(importState, 'Faxa_dstdry2', dataPtr, begg, rc=rc )
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    do g = begg,endg
       atm2lnd_inst%forc_aer_grc(g,10) = dataPtr(g-begg+1)
    end do

    call state_getfldptr(importState, 'Faxa_dstwet3', dataPtr, begg, rc=rc )
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    do g = begg,endg
       atm2lnd_inst%forc_aer_grc(g,11) = dataPtr(g-begg+1)
    end do

    call state_getfldptr(importState, 'Faxa_dstdry3', dataPtr, begg, rc=rc )
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    do g = begg,endg
       atm2lnd_inst%forc_aer_grc(g,12) = dataPtr(g-begg+1)
    end do

    call state_getfldptr(importState, 'Faxa_dstwet4', dataPtr, begg, rc=rc )
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    do g = begg,endg
       atm2lnd_inst%forc_aer_grc(g,13) = dataPtr(g-begg+1)
    end do

    call state_getfldptr(importState, 'Faxa_dstdry4', dataPtr, begg, rc=rc )
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    do g = begg,endg
       atm2lnd_inst%forc_aer_grc(g,14) = dataPtr(g-begg+1)
    end do

    !--------------------------
    ! Atmosphere co2 
    !--------------------------

    fldName = 'Sa_co2prog'  
    call ESMF_StateGet(importState, trim(fldname), itemFlag, rc=rc)
    if ( shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    if (itemflag == ESMF_STATEITEM_NOTFOUND .and. co2_type == 'prognostic') then
       call shr_sys_abort( subname//' ERROR: must have nonzero Sa_co2prog for co2_type equal to prognostic' )
    end if
    call ESMF_StateGet(importState, trim(fldname), itemFlag, rc=rc)
    if ( shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    if (itemflag /= ESMF_STATEITEM_NOTFOUND) then 
       call state_getfldptr(importState, trim(fldname), dataPtr, begg, rc=rc )
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       do g = begg,endg
          co2_ppmv_prog(g) = dataPtr(g-begg+1)   ! co2 atm prognostic
       end do
    else
       do g = begg,endg
          co2_ppmv_prog(g) = co2_ppmv
       end do
    end if

    fldName = 'Sa_co2diag'
    call ESMF_StateGet(importState, trim(fldname), itemFlag, rc=rc)
    if ( shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    if (itemflag == ESMF_STATEITEM_NOTFOUND .and. co2_type == 'diagnostic') then
       call shr_sys_abort( subname//' ERROR: must have nonzero Sa_co2prog for co2_type equal to prognostic' )
    end if
    call ESMF_StateGet(importState, trim(fldname), itemFlag, rc=rc)
    if ( shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    if (itemflag /= ESMF_STATEITEM_NOTFOUND) then 
       call state_getfldptr(importState, trim(fldname), dataPtr, begg, rc=rc )
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       do g = begg,endg
          co2_ppmv_diag(g) = dataPtr(g-begg+1)   ! co2 atm diagnostic
       end do
    else
       do g = begg,endg
          co2_ppmv_diag(g) = co2_ppmv
       end do
    end if

    ! Note that the following does unit conversions from ppmv to partial pressures (Pa)
    ! Note that forc_pbot is in Pa
    do g = begg,endg
       if (co2_type == 'prognostic') then
          co2_ppmv_val = co2_ppmv_prog(g)
       else if (co2_type == 'diagnostic') then
          co2_ppmv_val = co2_ppmv_diag(g)
       else
          co2_ppmv_val = co2_ppmv
       end if
       forc_pbot = atm2lnd_inst%forc_pbot_not_downscaled_grc(g)
       atm2lnd_inst%forc_pco2_grc(g) = co2_ppmv_val * 1.e-6_r8 * forc_pbot
       if (use_c13) then
          atm2lnd_inst%forc_pc13o2_grc(g) = co2_ppmv_val * c13ratio * 1.e-6_r8 * forc_pbot
       end if
    end do

    !--------------------------
    ! Atmosphere methane
    !--------------------------

    call ESMF_StateGet(importState, 'Sa_methane', itemFlag, rc=rc)
    if ( shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    if (itemflag /= ESMF_STATEITEM_NOTFOUND) then 
       call state_getfldptr(importState, 'Sa_methane', dataPtr, begg, rc=rc )
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       do g = begg,endg
          atm2lnd_inst%forc_pch4_grc(g) = dataPtr(g-begg+1)
       end do
    end if

    !--------------------------
    ! Atmosphere nitrogen deposition
    !--------------------------

    ! The mediator is sending ndep in units if kgN/m2/s - and clm
    ! uses units of gN/m2/sec - so the following conversion needs to happen

    call ESMF_StateGet(importState, 'Faxa_nhx', itemFlag1, rc=rc)
    if ( shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_StateGet(importState, 'Faxa_noy', itemFlag2, rc=rc)
    if ( shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    if ( itemflag1 /= ESMF_STATEITEM_NOTFOUND .and. itemflag2 /= ESMF_STATEITEM_NOTFOUND) then
       call state_getfldptr(importState, 'Faxa_nhx', dataPtr1, begg, rc=rc )
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       call state_getfldptr(importState, 'Faxa_noy', dataPtr2, begg, rc=rc )
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       do g = begg,endg
          atm2lnd_inst%forc_ndep_grc(g) = (dataPtr1(g-begg+1) + dataPtr2(g-begg+1))*1000._r8
       end do
    end if

    !--------------------------
    ! Flooding back from river
    !--------------------------

    ! sign convention is positive downward and hierarchy is atm/glc/lnd/rof/ice/ocn.
    ! so water sent from rof to land is negative,
    ! change the sign to indicate addition of water to system.

    call state_getfldptr(importState, 'Flrr_flood', dataPtr, begg, rc=rc )
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    do g = begg,endg
       atm2lnd_inst%forc_flood_grc(g) = -dataPtr(g-begg+1)
    end do

    call state_getfldptr(importState, 'Flrr_volr', dataPtr, begg, rc=rc )
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    do g = begg,endg
       atm2lnd_inst%volr_grc(g) = dataPtr(g-begg+1) * (ldomain%area(g) * 1.e6_r8)
    end do

    call state_getfldptr(importState, 'Flrr_volrmch', dataPtr, begg, rc=rc )
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    do g = begg,endg
       atm2lnd_inst%volrmch_grc(g) = dataPtr(g-begg+1) * (ldomain%area(g) * 1.e6_r8)
    end do

    !--------------------------
    ! Land-ice (glc) fields
    !--------------------------

    ! We could avoid setting these fields if glc_present is .false., if that would
    ! help with performance. (The downside would be that we wouldn't have these fields
    ! available for diagnostic purposes or to force a later T compset with dlnd.)

    do num = 0,glc_nec
       nec_str = glc_elevclass_as_string(num)

       call state_getfldptr(importState, 'Sg_ice_covered'//nec_str, dataPtr, begg, rc=rc )
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       do g = begg,endg
          frac_grc(g,num) = dataPtr(g-begg+1)
       end do

       call state_getfldptr(importState, 'Sg_topo'//nec_str, dataPtr, begg, rc=rc )
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       do g = begg,endg
          topo_grc(g,num) = dataPtr(g-begg+1)
       end do

       call state_getfldptr(importState, 'Flgg_hflx'//nec_str, dataPtr, begg, rc=rc )
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       do g = begg,endg
          hflx_grc(g,num) = dataPtr(g-begg+1)
       end do
    end do

    call state_getfldptr(importState,'Sg_icemask', dataPtr, begg, rc=rc )
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    do g = begg,endg
       icemask_grc(g) = dataPtr(g-begg+1)
    end do

    call state_getfldptr(importState,'Sg_icemask_coupled_fluxes', dataPtr, begg, rc=rc )
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    do g = begg,endg
       icemask_coupled_fluxes_grc(g) = dataPtr(g-begg+1)
    end do

    call glc2lnd_inst%set_glc2lnd_fields_nuopc( bounds, glc_present, &
         frac_grc, topo_grc, hflx_grc, icemask_grc, icemask_coupled_fluxes_grc)

    !--------------------------
    ! Derived quantities 
    !--------------------------

    do g = begg, endg
       forc_t    = atm2lnd_inst%forc_t_not_downscaled_grc(g)
       forc_q    = atm2lnd_inst%forc_q_not_downscaled_grc(g)
       forc_pbot = atm2lnd_inst%forc_pbot_not_downscaled_grc(g)

       atm2lnd_inst%forc_hgt_u_grc(g) = atm2lnd_inst%forc_hgt_grc(g)    !observational height of wind [m]
       atm2lnd_inst%forc_hgt_t_grc(g) = atm2lnd_inst%forc_hgt_grc(g)    !observational height of temperature [m]
       atm2lnd_inst%forc_hgt_q_grc(g) = atm2lnd_inst%forc_hgt_grc(g)    !observational height of humidity [m]

       atm2lnd_inst%forc_vp_grc(g) = forc_q * forc_pbot  / (0.622_r8 + 0.378_r8 * forc_q)

       atm2lnd_inst%forc_rho_not_downscaled_grc(g) = &
            (forc_pbot - 0.378_r8 * atm2lnd_inst%forc_vp_grc(g)) / (rair * forc_t)

       atm2lnd_inst%forc_po2_grc(g) = o2_molar_const * forc_pbot

       atm2lnd_inst%forc_wind_grc(g) = sqrt(atm2lnd_inst%forc_u_grc(g)**2 + atm2lnd_inst%forc_v_grc(g)**2)

       atm2lnd_inst%forc_solar_grc(g) = atm2lnd_inst%forc_solad_grc(g,1) + atm2lnd_inst%forc_solai_grc(g,1) + &
                                        atm2lnd_inst%forc_solad_grc(g,2) + atm2lnd_inst%forc_solai_grc(g,2)

       atm2lnd_inst%forc_rain_not_downscaled_grc(g)  = forc_rainc(g) + forc_rainl(g)
       atm2lnd_inst%forc_snow_not_downscaled_grc(g)  = forc_snowc(g) + forc_snowl(g)

       if (forc_t > SHR_CONST_TKFRZ) then
          e = esatw(tdc(forc_t))
       else
          e = esati(tdc(forc_t))
       end if
       qsat = 0.622_r8*e / (forc_pbot - 0.378_r8*e)

       ! modify specific humidity if precip occurs
       if (1==2) then
          if ((forc_rainc(g)+forc_rainl(g)) > 0._r8) then
             forc_q = 0.95_r8*qsat
            !forc_q = qsat
             atm2lnd_inst%forc_q_not_downscaled_grc(g) = forc_q
          endif
       endif

       atm2lnd_inst%forc_rh_grc(g) = 100.0_r8*(forc_q / qsat)
    end do

    !--------------------------
    ! Error checks
    !--------------------------

    ! Check that solar, specific-humidity and LW downward aren't negative
    do g = begg,endg
       if ( atm2lnd_inst%forc_lwrad_not_downscaled_grc(g) <= 0.0_r8 ) then
          call shr_sys_abort( subname//&
               ' ERROR: Longwave down sent from the atmosphere model is negative or zero' )
       end if
       if ( (atm2lnd_inst%forc_solad_grc(g,1) < 0.0_r8) .or. &
            (atm2lnd_inst%forc_solad_grc(g,2) < 0.0_r8) .or. &
            (atm2lnd_inst%forc_solai_grc(g,1) < 0.0_r8) .or. &
            (atm2lnd_inst%forc_solai_grc(g,2) < 0.0_r8) ) then
          call shr_sys_abort( subname//&
               ' ERROR: One of the solar fields (indirect/diffuse, vis or near-IR)'// &
               ' from the atmosphere model is negative or zero' )
       end if
       if ( atm2lnd_inst%forc_q_not_downscaled_grc(g) < 0.0_r8 )then
          call shr_sys_abort( subname//&
               ' ERROR: Bottom layer specific humidty sent from the atmosphere model is less than zero' )
       end if
    end do

    ! Make sure relative humidity is properly bounded
    ! atm2lnd_inst%forc_rh_grc(g) = min( 100.0_r8, atm2lnd_inst%forc_rh_grc(g) )
    ! atm2lnd_inst%forc_rh_grc(g) = max(   0.0_r8, atm2lnd_inst%forc_rh_grc(g) )

    !--------------------------
    ! Debug output
    !--------------------------

    if (masterproc .and. debug_import > 0 .and. get_nstep() < 5) then
       do g = bounds%begg,bounds%endg
          i = 1 + (g-bounds%begg)
          write(iulog,F01)'import: nstep, n, Sa_z        = ',get_nstep(), i, atm2lnd_inst%forc_hgt_grc(g)
          write(iulog,F01)'import: nstep, n, Sa_topo     = ',get_nstep(), i, atm2lnd_inst%forc_topo_grc(g)
          write(iulog,F01)'import: nstep, n, Sa_u        = ',get_nstep(), i, atm2lnd_inst%forc_u_grc(g)
          write(iulog,F01)'import: nstep, n, Sa_v        = ',get_nstep(), i, atm2lnd_inst%forc_v_grc(g)
          write(iulog,F01)'import: nstep, n, Fa_swndr    = ',get_nstep(), i, atm2lnd_inst%forc_solad_grc(g,2)
          write(iulog,F01)'import: nstep, n, Fa_swvdr    = ',get_nstep(), i, atm2lnd_inst%forc_solad_grc(g,1)
          write(iulog,F01)'import: nstep, n, Fa_swndf    = ',get_nstep(), i, atm2lnd_inst%forc_solai_grc(g,2)
          write(iulog,F01)'import: nstep, n, Fa_swvdf    = ',get_nstep(), i, atm2lnd_inst%forc_solai_grc(g,1)
          write(iulog,F01)'import: nstep, n, Sa_ptem     = ',get_nstep(), i, atm2lnd_inst%forc_th_not_downscaled_grc(g)
          write(iulog,F01)'import: nstep, n, Sa_shum     = ',get_nstep(), i, atm2lnd_inst%forc_q_not_downscaled_grc(g)
          write(iulog,F01)'import: nstep, n, Sa_pbot     = ',get_nstep(), i, atm2lnd_inst%forc_pbot_not_downscaled_grc(g)
          write(iulog,F01)'import: nstep, n, Sa_tbot     = ',get_nstep(), i, atm2lnd_inst%forc_t_not_downscaled_grc(g)
          write(iulog,F01)'import: nstep, n, Fa_lwdn     = ',get_nstep(), i, atm2lnd_inst%forc_lwrad_not_downscaled_grc(g)
          write(iulog,F01)'import: nstep, n, Fa_rainc    = ',get_nstep(), i, forc_rainc(g)
          write(iulog,F01)'import: nstep, n, Fa_rainl    = ',get_nstep(), i, forc_rainl(g)
          write(iulog,F01)'import: nstep, n, Fa_snowc    = ',get_nstep(), i, forc_snowc(g)
          write(iulog,F01)'import: nstep, n, Fa_snowl    = ',get_nstep(), i, forc_snowl(g)
          write(iulog,F01)'import: nstep, n, Fa_bcphidry = ',get_nstep(), i, atm2lnd_inst%forc_aer_grc(g,1)
          write(iulog,F01)'import: nstep, n, Fa_fcphodry = ',get_nstep(), i, atm2lnd_inst%forc_aer_grc(g,2)
          write(iulog,F01)'import: nstep, n, Fa_bcphiwet = ',get_nstep(), i, atm2lnd_inst%forc_aer_grc(g,3)
          write(iulog,F01)'import: nstep, n, Fa_ocphidry = ',get_nstep(), i, atm2lnd_inst%forc_aer_grc(g,4)
          write(iulog,F01)'import: nstep, n, Fa_ocphodry = ',get_nstep(), i, atm2lnd_inst%forc_aer_grc(g,5)
          write(iulog,F01)'import: nstep, n, Fa_ocphiwet = ',get_nstep(), i, atm2lnd_inst%forc_aer_grc(g,6)
          write(iulog,F01)'import: nstep, n, Fa_dstwet1  = ',get_nstep(), i, atm2lnd_inst%forc_aer_grc(g,7)
          write(iulog,F01)'import: nstep, n, Fa_dstdry1  = ',get_nstep(), i, atm2lnd_inst%forc_aer_grc(g,8)
          write(iulog,F01)'import: nstep, n, Fa_dstwet2  = ',get_nstep(), i, atm2lnd_inst%forc_aer_grc(g,9)
          write(iulog,F01)'import: nstep, n, Fa_dstdry2  = ',get_nstep(), i, atm2lnd_inst%forc_aer_grc(g,10)
          write(iulog,F01)'import: nstep, n, Fa_dstwet3  = ',get_nstep(), i, atm2lnd_inst%forc_aer_grc(g,11)
          write(iulog,F01)'import: nstep, n, Fa_dstdry3  = ',get_nstep(), i, atm2lnd_inst%forc_aer_grc(g,12)
          write(iulog,F01)'import: nstep, n, Fa_dstwet4  = ',get_nstep(), i, atm2lnd_inst%forc_aer_grc(g,13)
          write(iulog,F01)'import: nstep, n, Fa_dstdry4  = ',get_nstep(), i, atm2lnd_inst%forc_aer_grc(g,14)
          write(iulog,F01)'import: nstep, n, Fa_co2      = ',get_nstep(), i, atm2lnd_inst%forc_pco2_grc(g) 
          write(iulog,F01)'import: nstep, n, Fa_ndep     = ',get_nstep(), i, atm2lnd_inst%forc_ndep_grc(g)
       end do
    end if

  end subroutine import_fields

  !===============================================================================

  subroutine export_fields( gcomp, bounds, lnd2atm_inst, lnd2glc_inst, rc)

    !-------------------------------
    ! Pack the export state
    !-------------------------------

    ! input/output variables
    type(ESMF_GridComp)               :: gcomp
    type(bounds_type) , intent(in)    :: bounds       ! bounds
    type(lnd2atm_type), intent(inout) :: lnd2atm_inst ! land to atmosphere exchange data type
    type(lnd2glc_type), intent(inout) :: lnd2glc_inst ! land to atmosphere exchange data type
    integer           , intent(out)   :: rc

    ! local variables
    type(ESMF_State)          :: exportState
    type(ESMF_StateItem_Flag) :: itemFlag
    character(len=128)        :: fldname
    character(len=2)          :: nec_str
    integer                   :: i, g, num
    integer                   :: begg, endg
    real(r8), pointer         :: dataPtr(:)
    character(len=*), parameter :: subname='(lnd_import_export:export_fields)'
    !---------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    ! Get export state
    call NUOPC_ModelGet(gcomp, exportState=exportState, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Set bounds
    begg = bounds%begg; endg=bounds%endg

    call state_getfldptr(exportState, 'Sl_lfrin', dataPtr, begg, rc=rc )
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    do g = begg,endg
       dataPtr(g-begg+1) = ldomain%frac(g)
    end do

    call state_getfldptr(exportState, 'Sl_t', dataPtr, begg, rc=rc )
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    do g = begg,endg
       dataPtr(g-begg+1) = lnd2atm_inst%t_rad_grc(g)
    end do

    call state_getfldptr(exportState, 'Sl_snowh', dataPtr, begg, rc=rc )
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    do g = begg,endg
       dataPtr(g-begg+1)= lnd2atm_inst%h2osno_grc(g)
    end do

    call state_getfldptr(exportState, 'Sl_avsdr', dataPtr, begg, rc=rc )
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    do g = begg,endg
       dataPtr(g-begg+1) = lnd2atm_inst%albd_grc(g,1)
    end do

    call state_getfldptr(exportState, 'Sl_anidr', dataPtr, begg, rc=rc )
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    do g = begg,endg
       dataPtr(g-begg+1) =  lnd2atm_inst%albd_grc(g,2)
    end do

    call state_getfldptr(exportState, 'Sl_avsdf', dataPtr, begg, rc=rc )
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    do g = begg,endg
       dataPtr(g-begg+1) =  lnd2atm_inst%albi_grc(g,1)
    end do

    call state_getfldptr(exportState, 'Sl_anidf', dataPtr, begg, rc=rc )
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    do g = begg,endg
       dataPtr(g-begg+1) =  lnd2atm_inst%albi_grc(g,2)
    end do

    call state_getfldptr(exportState, 'Sl_tref', dataPtr, begg, rc=rc )
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    do g = begg,endg
       dataPtr(g-begg+1) =  lnd2atm_inst%t_ref2m_grc(g)
    end do

    call state_getfldptr(exportState, 'Sl_qref', dataPtr, begg, rc=rc )
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    do g = begg,endg
       dataPtr(g-begg+1) =  lnd2atm_inst%q_ref2m_grc(g)
    end do

    call state_getfldptr(exportState, 'Sl_u10', dataPtr, begg, rc=rc )
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    do g = begg,endg
       dataPtr(g-begg+1) =  lnd2atm_inst%u_ref10m_grc(g)
    end do

    call state_getfldptr(exportState, 'Fall_taux', dataPtr, begg, rc=rc )
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    do g = begg,endg
       dataPtr(g-begg+1) = -lnd2atm_inst%taux_grc(g)
    end do

    call state_getfldptr(exportState, 'Fall_tauy', dataPtr, begg, rc=rc )
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    do g = begg,endg
       dataPtr(g-begg+1) = -lnd2atm_inst%tauy_grc(g)
    end do

    call state_getfldptr(exportState, 'Fall_lat', dataPtr, begg, rc=rc )
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    do g = begg,endg
       dataPtr(g-begg+1) = -lnd2atm_inst%eflx_lh_tot_grc(g)
    end do

    call state_getfldptr(exportState, 'Fall_sen', dataPtr, begg, rc=rc )
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    do g = begg,endg
       dataPtr(g-begg+1) = -lnd2atm_inst%eflx_sh_tot_grc(g)
    end do

    call state_getfldptr(exportState, 'Fall_lwup', dataPtr, begg, rc=rc )
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    do g = begg,endg
       dataPtr(g-begg+1) = -lnd2atm_inst%eflx_lwrad_out_grc(g)
    end do

    call state_getfldptr(exportState, 'Fall_evap', dataPtr, begg, rc=rc )
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    do g = begg,endg
       dataPtr(g-begg+1) = -lnd2atm_inst%qflx_evap_tot_grc(g)
    end do

    call state_getfldptr(exportState, 'Fall_swnet', dataPtr, begg, rc=rc )
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    do g = begg,endg
       dataPtr(g-begg+1) =  lnd2atm_inst%fsa_grc(g)
    end do

    if (NUOPC_IsConnected(exportState,  fieldName='Sl_ram1')) then
       call state_getfldptr(exportState,'Sl_ram1', dataPtr, begg, rc=rc )
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       do g = begg,endg
          dataPtr(g-begg+1) = lnd2atm_inst%ram1_grc(g)
       end do
    end if
    if (NUOPC_IsConnected(exportState, fieldName='Sl_fv')) then
       call state_getfldptr(exportState,'Sl_fv', dataPtr, begg, rc=rc )
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       do g = begg,endg
          dataPtr(g-begg+1) = lnd2atm_inst%fv_grc(g)
       end do
    end if

    call ESMF_StateGet(exportState, 'Sl_soilw', itemFlag, rc=rc)
    if ( shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    if (itemflag /= ESMF_STATEITEM_NOTFOUND) then 
       call state_getfldptr(exportState,'Sl_soilw', dataPtr, begg, rc=rc )
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       do g = begg,endg
          dataPtr(g-begg+1) = lnd2atm_inst%h2osoi_vol_grc(g,1)
       end do
    end if

    call ESMF_StateGet(exportState, 'Fall_flxdst1', itemFlag, rc=rc) 
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    if (itemflag /= ESMF_STATEITEM_NOTFOUND) then 
       call state_getfldptr(exportState,'Fall_flxdst1', dataPtr, begg, rc=rc )
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       do g = begg,endg
          dataPtr(g-begg+1) = -lnd2atm_inst%flxdst_grc(g,1)
       end do
    end if

    call ESMF_StateGet(exportState, 'Fall_flxdst2', itemFlag, rc=rc) 
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    if (itemflag /= ESMF_STATEITEM_NOTFOUND) then 
       call state_getfldptr(exportState,'Fall_flxdst2', dataPtr, begg, rc=rc )
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       do g = begg,endg
          dataPtr(g-begg+1) = -lnd2atm_inst%flxdst_grc(g,2)
       end do
    end if

    call ESMF_StateGet(exportState, 'Fall_flxdst3', itemFlag, rc=rc) 
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    if (itemflag /= ESMF_STATEITEM_NOTFOUND) then 
       call state_getfldptr(exportState,'Fall_flxdst3', dataPtr, begg, rc=rc )
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       do g = begg,endg
          dataPtr(g-begg+1) = -lnd2atm_inst%flxdst_grc(g,3)
       end do
    end if

    call ESMF_StateGet(exportState, 'Fall_flxdst4', itemFlag, rc=rc) 
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    if (itemflag /= ESMF_STATEITEM_NOTFOUND) then 
       call state_getfldptr(exportState,'Fall_flxdst4', dataPtr, begg, rc=rc )
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       do g = begg,endg
          dataPtr(g-begg+1) = -lnd2atm_inst%flxdst_grc(g,4)
       end do
    end if

    ! dry dep velocities
    do num = 1, shr_string_listGetNum(drydep_fields)
       call  shr_string_listGetName(drydep_fields, num, fldname)
       if (NUOPC_IsConnected(exportState,  fieldName=trim(fldname))) then
          call state_getfldptr(exportState, trim(fldname), dataPtr, begg, rc=rc )
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
          do g = begg,endg
             dataPtr(1+g-begg) = lnd2atm_inst%ddvel_grc(g,num)
          end do
       end if
    end do

    ! MEGAN VOC emis fluxes
    do num = 1, shr_megan_mechcomps_n
       call  shr_string_listGetName(megan_voc_fields, num, fldname)
       if (NUOPC_IsConnected(exportState,  fieldName=trim(fldname))) then
          call state_getfldptr(exportState, trim(fldname), dataPtr, begg, rc=rc )
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
          do g = begg,endg
             dataPtr(g-begg+1) = -lnd2atm_inst%flxvoc_grc(g,num)
          end do
       end if
    end do

    ! fire emis fluxes
    do num = 1, shr_string_listGetNum(fire_emis_fields)
       call  shr_string_listGetName(fire_emis_fields, num, fldname)
       if (NUOPC_IsConnected(exportState,  fieldName=trim(fldname))) then
          call state_getfldptr(exportState, trim(fldname), dataPtr, begg, rc=rc )
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
          do g = begg,endg
             dataPtr(g-begg+1) = -lnd2atm_inst%fireflx_grc(g,num)
          end do
       end if
    end do
    if (NUOPC_IsConnected(exportState,  fieldName=shr_fire_emis_ztop_token)) then
       call state_getfldptr(exportState, trim(shr_fire_emis_ztop_token), dataPtr, begg, rc=rc )
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       do g = begg,endg
          dataPtr(g-begg+1) = lnd2atm_inst%fireztop_grc(g)
       end do
    end if

    ! methane
    call ESMF_StateGet(exportState, 'Fall_methane', itemFlag, rc=rc) 
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    if (itemflag /= ESMF_STATEITEM_NOTFOUND) then 
       call state_getfldptr(exportState, 'Fall_methane', dataPtr, begg, rc=rc )
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       do g = begg,endg
          dataPtr(g-begg+1) = -lnd2atm_inst%flux_ch4_grc(g)
       end do
    end if

    ! sign convention is positive downward with hierarchy of atm/glc/lnd/rof/ice/ocn.
    ! i.e. water sent from land to rof is positive

    ! surface runoff is the sum of qflx_over, qflx_h2osfc_surf
    call state_getfldptr(exportState, 'Flrl_rofsur', dataPtr, begg, rc=rc )
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    do g = begg,endg
       dataPtr(g-begg+1) = lnd2atm_inst%qflx_rofliq_qsur_grc(g) + lnd2atm_inst%qflx_rofliq_h2osfc_grc(g)
    end do

    !  subsurface runoff is the sum of qflx_drain and qflx_perched_drain
    call state_getfldptr(exportState, 'Flrl_rofsub', dataPtr, begg, rc=rc )
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    do g = begg,endg
       dataPtr(g-begg+1) = lnd2atm_inst%qflx_rofliq_qsub_grc(g) + lnd2atm_inst%qflx_rofliq_drain_perched_grc(g)
    end do

    !  qgwl sent individually to coupler
    call state_getfldptr(exportState, 'Flrl_rofgwl', dataPtr, begg, rc=rc )
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    do g = begg,endg
       dataPtr(g-begg+1) = lnd2atm_inst%qflx_rofliq_qgwl_grc(g)
    end do

    ! ice  sent individually to coupler
    call state_getfldptr(exportState, 'Flrl_rofi', dataPtr, begg, rc=rc )
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    do g = begg,endg
       dataPtr(g-begg+1) = lnd2atm_inst%qflx_rofice_grc(g)
    end do

    ! irrigation flux to be removed from main channel storage (negative)
    call state_getfldptr(exportState, 'Flrl_irrig', dataPtr, begg, rc=rc )
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    do g = begg,endg
       dataPtr(g-begg+1) = - lnd2atm_inst%qirrig_grc(g)
    end do

    ! glc fields
    ! We could avoid setting these fields if glc_present is .false., if that would
    ! help with performance. (The downside would be that we wouldn't have these fields
    ! available for diagnostic purposes or to force a later T compset with dlnd.)

    do num = 0,glc_nec
       nec_str = glc_elevclass_as_string(num)

       fldname = 'Sl_tsrf' // nec_str
       call state_getfldptr(exportState, trim(fldname), dataPtr, begg, rc=rc )
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       do g = begg,endg
          dataPtr(g-begg+1) = lnd2glc_inst%tsrf_grc(g,num)
       end do

       fldname = 'Sl_topo' // nec_str
       call state_getfldptr(exportState, trim(fldname), dataPtr, begg, rc=rc )
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       do g = begg,endg
          dataPtr(g-begg+1) = lnd2glc_inst%topo_grc(g,num)
       end do

       fldname = 'Flgl_qice' // nec_str
       call state_getfldptr(exportState, trim(fldname), dataPtr, begg, rc=rc )
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       do g = begg,endg
          dataPtr(g-begg+1) = lnd2glc_inst%qice_grc(g,num)
       end do
    end do

    if (masterproc .and. debug_export > 0 .and. get_nstep() < 5) then
       do g = begg,endg
          i = 1 + (g-bounds%begg)
          write(iulog,F01)'export:  nstep, n, Sl_t       =  ', get_nstep(), i, lnd2atm_inst%t_rad_grc(g)
          write(iulog,F01)'export:  nstep, n, Sl_snowh   =  ', get_nstep(), i, lnd2atm_inst%h2osno_grc(g)
          write(iulog,F01)'export:  nstep, n, Sl_avsdr   =  ', get_nstep(), i, lnd2atm_inst%albd_grc(g,1)
          write(iulog,F01)'export:  nstep, n, Sl_anidr   =  ', get_nstep(), i, lnd2atm_inst%albd_grc(g,2)
          write(iulog,F01)'export:  nstep, n, Sl_avsdf   =  ', get_nstep(), i, lnd2atm_inst%albi_grc(g,1)
          write(iulog,F01)'export:  nstep, n, Sl_anidf   =  ', get_nstep(), i, lnd2atm_inst%albi_grc(g,2)
          write(iulog,F01)'export:  nstep, n, Sl_tref    =  ', get_nstep(), i, lnd2atm_inst%t_ref2m_grc(g)
          write(iulog,F01)'export:  nstep, n, Sl_qref    =  ', get_nstep(), i, lnd2atm_inst%q_ref2m_grc(g)
          write(iulog,F01)'export:  nstep, n, Sl_u10     =  ', get_nstep(), i, lnd2atm_inst%u_ref10m_grc(g)
          write(iulog,F01)'export:  nstep, n, Fl_taux    =  ', get_nstep(), i,-lnd2atm_inst%taux_grc(g)
          write(iulog,F01)'export:  nstep, n, Fl_tauy    =  ', get_nstep(), i,-lnd2atm_inst%tauy_grc(g)
          write(iulog,F01)'export:  nstep, n, Fl_lat     =  ', get_nstep(), i,-lnd2atm_inst%eflx_lh_tot_grc(g)
          write(iulog,F01)'export:  nstep, n, Fl_sen     =  ', get_nstep(), i,-lnd2atm_inst%eflx_sh_tot_grc(g)
          write(iulog,F01)'export:  nstep, n, Fl_lwup    =  ', get_nstep(), i,-lnd2atm_inst%eflx_lwrad_out_grc(g)
          write(iulog,F01)'export:  nstep, n, Fl_evap    =  ', get_nstep(), i,-lnd2atm_inst%qflx_evap_tot_grc(g)
          write(iulog,F01)'export:  nstep, n, Fl_swnet   =  ', get_nstep(), i, lnd2atm_inst%fsa_grc(g)
          write(iulog,F01)'export:  nstep, n, Sl_ram1    =  ', get_nstep(), i, lnd2atm_inst%ram1_grc(g)
          write(iulog,F01)'export:  nstep, n, Sl_fv      =  ', get_nstep(), i, lnd2atm_inst%fv_grc(g)
          write(iulog,F01)'export:  nstep, n, Fl_flxdst1 =  ', get_nstep(), i, lnd2atm_inst%flxdst_grc(g,1)
          write(iulog,F01)'export:  nstep, n, Fl_flxdst2 =  ', get_nstep(), i, lnd2atm_inst%flxdst_grc(g,2)
          write(iulog,F01)'export:  nstep, n, Fl_flxdst3 =  ', get_nstep(), i,-lnd2atm_inst%flxdst_grc(g,3)
          write(iulog,F01)'export:  nstep, n, Fl_flxdst4 =  ', get_nstep(), i,-lnd2atm_inst%flxdst_grc(g,4)
          write(iulog,F01)'export:  nstep, n, Sl_fv       = ', get_nstep(), i, lnd2atm_inst%fv_grc(g)
          write(iulog,F01)'export:  nstep, n, Flrl_rofsur = ', get_nstep(), i, lnd2atm_inst%qflx_rofliq_qsur_grc(g) + &
                                                                               lnd2atm_inst%qflx_rofliq_h2osfc_grc(g)
          write(iulog,F01)'export:  nstep, n, Flrl_rofsub = ', get_nstep(), i, lnd2atm_inst%qflx_rofliq_qsub_grc(g) + &
                                                                               lnd2atm_inst%qflx_rofliq_drain_perched_grc(g)
          write(iulog,F01)'export:  nstep, n, Flrl_rofgwl = ', get_nstep(), i, lnd2atm_inst%qflx_rofliq_qgwl_grc(g)
          write(iulog,F01)'export:  nstep, n, Flrl_rofi   = ', get_nstep(), i, lnd2atm_inst%qflx_rofice_grc(g)
          write(iulog,F01)'export:  nstep, n, Flrl_irrig  = ', get_nstep(), i,-lnd2atm_inst%qirrig_grc(g)
       end do
    end if

  end subroutine export_fields

  !===============================================================================

  subroutine fldlist_add(num, fldlist, stdname)
    integer,                    intent(inout) :: num
    type(fld_list_type),        intent(inout) :: fldlist(:)
    character(len=*),           intent(in)    :: stdname

    ! local variables
    integer :: rc
    integer :: dbrc
    character(len=*), parameter :: subname='(dshr_nuopc_mod:fldlist_add)'
    !-------------------------------------------------------------------------------

    ! Set up a list of field information

    num = num + 1
    if (num > fldsMax) then
       call ESMF_LogWrite(trim(subname)//": ERROR num > fldsMax "//trim(stdname), &
            ESMF_LOGMSG_ERROR, line=__LINE__, file=__FILE__, rc=dbrc)
       return
    endif
    fldlist(num)%stdname = trim(stdname)

  end subroutine fldlist_add

  !===============================================================================

  subroutine fldlist_realize(state, fldList, numflds, flds_scalar_name, flds_scalar_num, mesh, tag, rc)

    use NUOPC , only : NUOPC_IsConnected, NUOPC_Realize
    use ESMF  , only : ESMF_MeshLoc_Element, ESMF_FieldCreate, ESMF_TYPEKIND_R8
    use ESMF  , only : ESMF_MAXSTR, ESMF_Field, ESMF_State, ESMF_Mesh, ESMF_StateRemove
    use ESMF  , only : ESMF_LogFoundError, ESMF_LOGMSG_INFO, ESMF_SUCCESS
    use ESMF  , only : ESMF_LogWrite, ESMF_LOGMSG_ERROR, ESMF_LOGERR_PASSTHRU

    type(ESMF_State)    , intent(inout) :: state
    type(fld_list_type) , intent(in)    :: fldList(:)
    integer             , intent(in)    :: numflds
    character(len=*)    , intent(in)    :: flds_scalar_name
    integer             , intent(in)    :: flds_scalar_num
    character(len=*)    , intent(in)    :: tag
    type(ESMF_Mesh)     , intent(in)    :: mesh
    integer             , intent(inout) :: rc

    ! local variables
    integer                :: dbrc
    integer                :: n
    type(ESMF_Field)       :: field
    character(len=80)      :: stdname
    character(len=*),parameter  :: subname='(dshr_nuopc_mod:fldlist_realize)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    do n = 1, numflds
       stdname = fldList(n)%stdname
       if (NUOPC_IsConnected(state, fieldName=stdname)) then
          if (stdname == trim(flds_scalar_name)) then
             call ESMF_LogWrite(trim(subname)//trim(tag)//" Field = "//trim(stdname)//" is connected on root pe", &
                  ESMF_LOGMSG_INFO, rc=dbrc)
             ! Create the scalar field
             call SetScalarField(field, flds_scalar_name, flds_scalar_num, rc=rc)
             if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return
          else
             call ESMF_LogWrite(trim(subname)//trim(tag)//" Field = "//trim(stdname)//" is connected using mesh", &
                  ESMF_LOGMSG_INFO, rc=dbrc)
             ! Create the field
             field = ESMF_FieldCreate(mesh, ESMF_TYPEKIND_R8, name=stdname, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
             if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return
          endif

          ! NOW call NUOPC_Realize
          call NUOPC_Realize(state, field=field, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return
       else
          if (stdname /= trim(flds_scalar_name)) then
             call ESMF_LogWrite(subname // trim(tag) // " Field = "// trim(stdname) // " is not connected.", &
                  ESMF_LOGMSG_INFO, rc=dbrc)
             call ESMF_StateRemove(state, (/stdname/), rc=rc)
             if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return
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
      character(len=*), parameter :: subname='(dshr_nuopc_mod:SetScalarField)'
      ! ----------------------------------------------

      rc = ESMF_SUCCESS

      ! create a DistGrid with a single index space element, which gets mapped onto DE 0.
      distgrid = ESMF_DistGridCreate(minIndex=(/1/), maxIndex=(/1/), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return

      grid = ESMF_GridCreate(distgrid, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return

      field = ESMF_FieldCreate(name=trim(flds_scalar_name), grid=grid, typekind=ESMF_TYPEKIND_R8, &
           ungriddedLBound=(/1/), ungriddedUBound=(/flds_scalar_num/), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return

    end subroutine SetScalarField

  end subroutine fldlist_realize

  !===============================================================================

  subroutine state_set(State, fldname, bounds, array_import, array_export, rc)
    ! ----------------------------------------------
    ! Get pointer to a state field
    ! ----------------------------------------------
    use ESMF , only : ESMF_State, ESMF_Field, ESMF_Mesh, ESMF_FieldStatus_Flag
    use ESMF , only : ESMF_StateGet, ESMF_FieldGet, ESMF_MeshGet
    use ESMF , only : ESMF_FIELDSTATUS_COMPLETE, ESMF_FAILURE
    use ESMF , only : operator(/=)

    type(ESMF_State)    , intent(in)  :: State
    type(bounds_type)   , intent(in)  :: bounds       
    character(len=*)    , intent(in)  :: fldname
    real(r8), optional  , intent(out) :: array_import(bounds%begg:bounds%endg)
    real(r8), optional  , intent(in)  :: array_export(bounds%begg:bounds%endg)
    integer             , intent(out) :: rc

    ! local variables
    integer                     :: g 
    real(R8), pointer           :: fldptr(:)
    type(ESMF_FieldStatus_Flag) :: status
    type(ESMF_Field)            :: lfield
    type(ESMF_Mesh)             :: lmesh
    integer                     :: dbrc
    integer                     :: nnodes, nelements
    character(len=*), parameter :: subname='(lnd_import_export:State_GetFldPtr)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    if (dbug > 10) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

    call ESMF_StateGet(State, itemName=trim(fldname), field=lfield, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_FieldGet(lfield, status=status, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (status /= ESMF_FIELDSTATUS_COMPLETE) then
       call ESMF_LogWrite(trim(subname)//": ERROR data not allocated ", ESMF_LOGMSG_INFO, rc=rc)
       rc = ESMF_FAILURE
       return
    else
       call ESMF_FieldGet(lfield, mesh=lmesh, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

       call ESMF_MeshGet(lmesh, numOwnedNodes=nnodes, numOwnedElements=nelements, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

       if (nnodes == 0 .and. nelements == 0) then
          call ESMF_LogWrite(trim(subname)//": no local nodes or elements ", ESMF_LOGMSG_INFO, rc=dbrc)
          rc = ESMF_FAILURE
          return
       end if

       call ESMF_FieldGet(lfield, farrayPtr=fldptr, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    endif  ! status

    ! Now determine either array_in (for import state or array_out for export state)
    ! and check for nans

    if (present(array_import)) then
       do g = bounds%begg, bounds%endg
          array_import(g) = fldptr(g-bounds%begg+1)
       end do
    else if (present(array_export)) then
       do g = bounds%begg, bounds%endg
          fldptr(g-bounds%begg+1) = array_export(g)
       end do
    else
       call shr_sys_abort(trim(subname)//"either array_in or array_out must be present")
    end if
    call check_for_nans(fldptr, trim(fldname), bounds%begg)

    if (dbug > 10) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

  end subroutine state_set  

  !===============================================================================

  subroutine state_getfldptr(State, fldname, fldptr, begg, rc)
    ! ----------------------------------------------
    ! Get pointer to a state field
    ! ----------------------------------------------
    use ESMF , only : ESMF_State, ESMF_Field, ESMF_Mesh, ESMF_FieldStatus_Flag
    use ESMF , only : ESMF_StateGet, ESMF_FieldGet, ESMF_MeshGet
    use ESMF , only : ESMF_FIELDSTATUS_COMPLETE, ESMF_FAILURE

    type(ESMF_State),  intent(in)    :: State
    character(len=*),  intent(in)    :: fldname
    integer,           intent(in)    :: begg            
    real(R8), pointer, intent(out)   :: fldptr(:)
    integer,           intent(out)   :: rc

    ! local variables
    type(ESMF_FieldStatus_Flag) :: status
    type(ESMF_Field)            :: lfield
    type(ESMF_Mesh)             :: lmesh
    integer                     :: dbrc
    integer                     :: nnodes, nelements
    character(len=*), parameter :: subname='(lnd_import_export:state_getfldptr)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    if (dbug > 10) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

    call ESMF_StateGet(State, itemName=trim(fldname), field=lfield, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_FieldGet(lfield, status=status, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (status /= ESMF_FIELDSTATUS_COMPLETE) then
       call ESMF_LogWrite(trim(subname)//": ERROR data not allocated ", ESMF_LOGMSG_INFO, rc=rc)
       rc = ESMF_FAILURE
       return
    else
       call ESMF_FieldGet(lfield, mesh=lmesh, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

       call ESMF_MeshGet(lmesh, numOwnedNodes=nnodes, numOwnedElements=nelements, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

       if (nnodes == 0 .and. nelements == 0) then
          call ESMF_LogWrite(trim(subname)//": no local nodes or elements ", ESMF_LOGMSG_INFO, rc=dbrc)
          rc = ESMF_FAILURE
          return
       end if

       call ESMF_FieldGet(lfield, farrayPtr=fldptr, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    endif  ! status

    if (dbug > 10) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

  end subroutine state_getfldptr

  !===============================================================================

  subroutine check_for_nans(array, fname, begg)

    ! input/output variables
    real(r8), pointer             :: array(:)
    character(len=*) , intent(in) :: fname
    integer          , intent(in) :: begg

    ! local variables
    integer :: i
    !-------------------------------------------------------------------------------

    ! Check if any input from mediator or output to mediator is NaN

    if (any(isnan(array))) then
       write(iulog,*) '# of NaNs = ', count(isnan(array))
       write(iulog,*) 'Which are NaNs = ', isnan(array)
       do i = 1, size(array)
          if (isnan(array(i))) then
             write(iulog,*) "NaN found in field", trim(fname), ' at gridcell index ',begg+i-1 
          end if
       end do
       call shr_sys_abort(' ERROR: One or more of the output from CLM to the coupler are NaN ' )
    end if
  end subroutine check_for_nans

end module lnd_import_export
