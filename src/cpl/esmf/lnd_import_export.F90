module lnd_import_export

  ! TODO (mvertens, 2019-01-28): what to do about scalar fields?
  !use shr_nuopc_scalars_mod , only : flds_scalar_name, flds_scalar_num

  use ESMF                  , only : ESMF_GridComp, ESMF_State, ESMF_Mesh, ESMF_StateGet
  use ESMF                  , only : ESMF_KIND_R8, ESMF_SUCCESS, ESMF_MAXSTR, ESMF_LOGMSG_INFO
  use ESMF                  , only : ESMF_LogWrite, ESMF_LOGMSG_ERROR, ESMF_LogFoundError
  use ESMF                  , only : ESMF_STATEITEM_NOTFOUND, ESMF_StateItem_Flag
  use ESMF                  , only : operator(/=), operator(==)
  use shr_utils_mod         , only : ChkErr => shr_utils_ChkErr
  use shr_kind_mod          , only : r8 => shr_kind_r8, cx=>shr_kind_cx, cxx=>shr_kind_cxx
  use shr_infnan_mod        , only : isnan => shr_infnan_isnan
  use shr_string_mod        , only : shr_string_listGetName, shr_string_listGetNum
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
  use seq_drydep_mod        , only : seq_drydep_readnl, seq_drydep_init, n_drydep

  implicit none
  private ! except

  public  :: add_fields
  public  :: import_fields
  public  :: export_fields

  private :: fldlist_add
  private :: state_getimport
  private :: state_setexport
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

  character(len=cxx)     :: drydep_fields    ! List of dry-deposition fields
  character(len=cxx)     :: megan_voc_fields ! List of MEGAN VOC emission fields
  character(len=cxx)     :: fire_emis_fields ! List of fire emission fields
  character(len=cx)      :: carma_fields     ! List of CARMA fields from lnd->atm
  character(len=cx)      :: ndep_fields      ! List of nitrogen deposition fields from atm->lnd/ocn
  logical                :: flds_co2a        ! use case
  logical                :: flds_co2b        ! use case
  logical                :: flds_co2c        ! use case
  integer                :: glc_nec
  integer, parameter     :: debug = 0        ! internal debug level
  character(*),parameter :: F01 = "('(lnd_import_export) ',a,i5,2x,i5,2x,d21.14)"
  character(*),parameter :: u_FILE_u = &
       __FILE__

!===============================================================================
contains
!===============================================================================

  subroutine add_fields(gcomp, rc)

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
    integer                :: drydep_nflds, ndep_nflds, megan_nflds, emis_nflds
    character(len=128)     :: fldname
    character(len=*), parameter :: subname='(lnd_import_export:advertise_fields)'
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    call ESMF_GridCompGet(gcomp, importState=importState, exportState=exportState, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! TODO (mvertens, 2019-01-28): what to do about scalar fields?

    !--------------------------------
    ! Add default export fields
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

    ! Create fields in export state

    do n = 1,size(fldsFrLnd)
       name = trim(fldsFrLnd(n)%stdname)
       field = ESMF_FieldCreate(mesh, ESMF_TYPEKIND_R8, meshloc=ESMF_MESHLOC_ELEMENT, name=trim(name), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_StateAdd(ExportState, field, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end do

    !--------------------------------
    ! Add import fields
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

    ! from atm - nitrogen deposition
    ! TODO (mvertens, 2019-01-28): what to do here?
    ! call shr_ndep_readnl("drv_flds_in", ndep_fields, ndep_nflds)
    ! if (ndep_nflds>0) then
    !    do n = 1, ndep_nflds
    !       call  shr_string_listGetName(ndep_fields, n, fldname)
    !       call fldlist_add(fldsToLnd_num, fldsToLnd, trim(fldname))
    !    end do
    !    ! This sets a variable in clm_varctl
    !    ndep_from_cpl = .true.
    ! end if

    ! Create fields in import state

    do n = 1,size(fldsToLnd)
       name = trim(fldsToLnd(n)%stdname)
       field = ESMF_FieldCreate(mesh, ESMF_TYPEKIND_R8, meshloc=ESMF_MESHLOC_ELEMENT, name=trim(name), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_StateAdd(ImportState, field, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end do

  end subroutine advertise_fields

  !===============================================================================

  subroutine import_fields( gcomp, bounds, glc_present, atm2lnd_inst, glc2lnd_inst, wateratm2lndbulk_inst, rc)

    !---------------------------------------------------------------------------
    ! Convert the input data from the mediator to the land model
    !---------------------------------------------------------------------------

    use clm_varctl           , only: co2_type, co2_ppmv, use_c13, ndep_from_cpl
    use clm_varcon           , only: rair, o2_molar_const, c13ratio
    use shr_const_mod        , only: SHR_CONST_TKFRZ
    use Wateratm2lndBulkType , only: wateratm2lndbulk_type

    ! input/output variabes
    type(ESMF_GridComp)                :: gcomp
    type(bounds_type)  , intent(in)    :: bounds       ! bounds
    logical            , intent(in)    :: glc_present  ! .true. => running with a non-stub GLC model
    type(atm2lnd_type) , intent(inout) :: atm2lnd_inst ! clm internal input data type
    type(glc2lnd_type) , intent(inout) :: glc2lnd_inst ! clm internal input data type
    type(Wateratm2lndbulk_type), intent(inout) :: wateratm2lndbulk_inst
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
    real(r8)                  :: forc_noy(bounds%begg:bounds%endg)
    real(r8)                  :: forc_nhx(bounds%begg:bounds%endg)
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
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO, rc=dbrc)

    ! Get import state
    call ESMF_GridCompGet(gcomp, importState=importState, rc=rc)
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

    call state_getimport(importState, 'Faxa_bcphidry', bounds, output=atm2lnd_inst%forc_aer_grc(:,1), rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_getimport(importState, 'Faxa_bcphodry', bounds, output=atm2lnd_inst%forc_aer_grc(:,2), rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_getimport(importState, 'Faxa_bcphiwet', bounds, output=atm2lnd_inst%forc_aer_grc(:,3), rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_getimport(importState, 'Faxa_ocphidry', bounds, output=atm2lnd_inst%forc_aer_grc(:,4), rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_getimport(importState, 'Faxa_ocphodry', bounds, output=atm2lnd_inst%forc_aer_grc(:,5), rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_getimport(importState, 'Faxa_ocphiwet', bounds, output=atm2lnd_inst%forc_aer_grc(:,6), rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_getimport(importState, 'Faxa_dstwet1', bounds, output=atm2lnd_inst%forc_aer_grc(:,7), rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_getimport(importState, 'Faxa_dstdry1', bounds, output=atm2lnd_inst%forc_aer_grc(:,8), rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_getimport(importState, 'Faxa_dstwet2', bounds, output=atm2lnd_inst%forc_aer_grc(:,9), rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_getimport(importState, 'Faxa_dstdry2', bounds, output=atm2lnd_inst%forc_aer_grc(:,10), rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_getimport(importState, 'Faxa_dstwet3', bounds, output=atm2lnd_inst%forc_aer_grc(:,11), rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_getimport(importState, 'Faxa_dstdry3', bounds, output=atm2lnd_inst%forc_aer_grc(:,12), rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_getimport(importState, 'Faxa_dstwet4', bounds, output=atm2lnd_inst%forc_aer_grc(:,13), rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_getimport(importState, 'Faxa_dstdry4', bounds, output=atm2lnd_inst%forc_aer_grc(:,14), rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_getimport(importState, 'Sa_methane', bounds, output=atm2lnd_inst%forc_pch4_grc, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! The mediator is sending ndep in units if kgN/m2/s - and ctsm uses units of gN/m2/sec
    ! so the following conversion needs to happen
    call state_getimport(importState, 'Faxa_nhx', bounds, output=forc_nhx, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport(importState, 'Faxa_noy', bounds, output=forc_noy, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    do g = begg,endg
       atm2lnd_inst%forc_ndep_grc(g) = (forc_nhx(g) + forc_noy(g))*1000._r8
    end do

    !--------------------------
    ! Atmosphere co2
    !--------------------------

    fldName = 'Sa_co2prog'
    call ESMF_StateGet(importState, trim(fldname), itemFlag, rc=rc)
    if ( ChkErr(rc,__LINE__,u_FILE_u)) return
    if (itemflag == ESMF_STATEITEM_NOTFOUND .and. co2_type == 'prognostic') then
       call shr_sys_abort( subname//' ERROR: must have nonzero Sa_co2prog for co2_type equal to prognostic' )
    end if
    if (itemflag /= ESMF_STATEITEM_NOTFOUND) then
       call state_getfldptr(importState, trim(fldname), dataPtr, rc=rc )
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
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
    if ( ChkErr(rc,__LINE__,u_FILE_u)) return
    if (itemflag == ESMF_STATEITEM_NOTFOUND .and. co2_type == 'diagnostic') then
       call shr_sys_abort( subname//' ERROR: must have nonzero Sa_co2prog for co2_type equal to prognostic' )
    end if
    if (itemflag /= ESMF_STATEITEM_NOTFOUND) then
       call state_getfldptr(importState, trim(fldname), dataPtr, rc=rc )
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
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
    ! Derived quantities
    !--------------------------

    do g = begg, endg
       forc_t    = atm2lnd_inst%forc_t_not_downscaled_grc(g)
       forc_q    = wateratm2lndbulk_inst%forc_q_not_downscaled_grc(g)
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

       wateratm2lndbulk_inst%forc_rain_not_downscaled_grc(g)  = forc_rainc(g) + forc_rainl(g)
       wateratm2lndbulk_inst%forc_snow_not_downscaled_grc(g)  = forc_snowc(g) + forc_snowl(g)

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
             wateratm2lndbulk_inst%forc_q_not_downscaled_grc(g) = forc_q
          endif
       endif

       wateratm2lndbulk_inst%forc_rh_grc(g) = 100.0_r8*(forc_q / qsat)
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
       if ( wateratm2lndbulk_inst%forc_q_not_downscaled_grc(g) < 0.0_r8 )then
          call shr_sys_abort( subname//&
               ' ERROR: Bottom layer specific humidty sent from the atmosphere model is less than zero' )
       end if
    end do

    ! Make sure relative humidity is properly bounded
    ! atm2lnd_inst%forc_rh_grc(g) = min( 100.0_r8, atm2lnd_inst%forc_rh_grc(g) )
    ! atm2lnd_inst%forc_rh_grc(g) = max(   0.0_r8, atm2lnd_inst%forc_rh_grc(g) )

  end subroutine import_fields

  !===============================================================================

  subroutine export_fields( gcomp, bounds, waterlnd2atmbulk_inst, lnd2atm_inst, rc)

    use Waterlnd2atmBulkType , only: waterlnd2atmbulk_type

    !-------------------------------
    ! Pack the export state
    !-------------------------------

    ! input/output variables
    type(ESMF_GridComp)                         :: gcomp
    type(bounds_type)           , intent(in)    :: bounds       ! bounds
    type(waterlnd2atmbulk_type) , intent(inout) :: waterlnd2atmbulk_inst
    type(lnd2atm_type)          , intent(inout) :: lnd2atm_inst ! land to atmosphere exchange data type
    integer                     , intent(out)   :: rc

    ! local variables
    type(ESMF_State)   :: exportState
    character(len=128) :: fldname
    character(len=2)   :: nec_str
    integer            :: i, g, num
    real(r8)           :: array(bounds%begg:bounds%endg)
    character(len=*), parameter :: subname='(lnd_import_export:export_fields)'
    !---------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    ! Get export state
    call ESMF_GridCompGet(gcomp, exportState=exportState, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Set field pointers in export set

    call state_setexport(exportState, 'Sl_lfrin', bounds, input=ldomain%frac, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_setexport(exportState, 'Sl_t', bounds, input=lnd2atm_inst%t_rad_grc, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_setexport(exportState, 'Sl_snowh', bounds, &
         input=waterlnd2atmbulk_inst%h2osno_grc, rc=rc)
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

    call state_setexport(exportState, 'Sl_qref', bounds, &
         input=waterlnd2atmbulk_inst%q_ref2m_grc, rc=rc)
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

    call state_setexport(exportState, 'Fall_evap', bounds, &
         input=waterlnd2atmbulk_inst%qflx_evap_tot_grc, minus=.true., rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_setexport(exportState, 'Fall_swnet', bounds, input=lnd2atm_inst%fsa_grc, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_setexport(exportState, 'Fall_flxdst1', bounds, input=lnd2atm_inst%flxdst_grc(:,1), minus=.true., rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_setexport(exportState, 'Fall_flxdst2', bounds, input=lnd2atm_inst%flxdst_grc(:,2), minus=.true., rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_setexport(exportState, 'Fall_flxdst3', bounds, input=lnd2atm_inst%flxdst_grc(:,3), minus=.true., rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_setexport(exportState, 'Fall_flxdst4', bounds, input=lnd2atm_inst%flxdst_grc(:,4), minus=.true., rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_setexport(exportState, 'Fall_methane', bounds, input=lnd2atm_inst%flux_ch4_grc, minus=.true., rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_setexport(exportState, 'Sl_ram1', bounds, input=lnd2atm_inst%ram1_grc, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_setexport(exportState, 'Sl_fv', bounds, input=lnd2atm_inst%fv_grc, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_setexport(exportState, 'Sl_soilw', bounds, &
         input=waterlnd2atmbulk_inst%h2osoi_vol_grc(:,1), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

  end subroutine export_fields

  !===============================================================================

  subroutine fldlist_add(num, fldlist, stdname)
    integer,                    intent(inout) :: num
    type(fld_list_type),        intent(inout) :: fldlist(:)
    character(len=*),           intent(in)    :: stdname

    ! local variables
    integer :: rc
    integer :: dbrc
    character(len=*), parameter :: subname='(lnd_import_export:fldlist_add)'
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

  subroutine state_getimport(state, fldname, bounds, output, rc)

    ! ----------------------------------------------
    ! Map import state field to output array
    ! ----------------------------------------------

    ! input/output variables
    type(ESMF_State)    , intent(in)    :: state
    type(bounds_type)   , intent(in)    :: bounds
    character(len=*)    , intent(in)    :: fldname
    real(r8)            , intent(out)   :: output(bounds%begg:bounds%endg)
    integer             , intent(out)   :: rc

    ! local variables
    integer                     :: g, i
    real(R8), pointer           :: fldptr(:)
    type(ESMF_StateItem_Flag)   :: itemFlag
    integer                     :: dbrc
    character(len=*), parameter :: subname='(lnd_import_export:state_getimport)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    ! Determine if field with name fldname exists in state
    call ESMF_StateGet(state, trim(fldname), itemFlag, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! if field exists then create output array - else do nothing
    if (itemflag /= ESMF_STATEITEM_NOTFOUND) then

       ! get field pointer
       call state_getfldptr(state, trim(fldname), fldptr, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       ! determine output array
       do g = bounds%begg, bounds%endg
          output(g) = fldptr(g - bounds%begg + 1)
       end do

       ! write debug output if appropriate
       if (masterproc .and. debug > 0 .and. get_nstep() < 5) then
          do g = bounds%begg,bounds%endg
             i = 1 + g - bounds%begg
             write(iulog,F01)'import: nstep, n, '//trim(fldname)//' = ',get_nstep(),i,output(g)
          end do
       end if

       ! check for nans
       call check_for_nans(fldptr, trim(fldname), bounds%begg)
    end if

  end subroutine state_getimport

  !===============================================================================

  subroutine state_setexport(state, fldname, bounds, input, minus, rc)

    use shr_const_mod, only : fillvalue=>SHR_CONST_SPVAL

    ! ----------------------------------------------
    ! Map input array to export state field
    ! ----------------------------------------------

    ! input/output variables
    type(ESMF_State)    , intent(inout) :: state
    type(bounds_type)   , intent(in)    :: bounds
    character(len=*)    , intent(in)    :: fldname
    real(r8)            , intent(in)    :: input(bounds%begg:bounds%endg)
    logical, optional   , intent(in)    :: minus
    integer             , intent(out)   :: rc

    ! local variables
    integer                     :: g, i
    real(R8), pointer           :: fldptr(:)
    type(ESMF_StateItem_Flag)   :: itemFlag
    integer                     :: dbrc
    character(len=*), parameter :: subname='(lnd_import_export:state_setexport)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    ! Determine if field with name fldname exists in state
    call ESMF_StateGet(state, trim(fldname), itemFlag, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! if field exists then create output array - else do nothing
    if (itemflag /= ESMF_STATEITEM_NOTFOUND) then

       ! get field pointer
       call state_getfldptr(state, trim(fldname), fldptr, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       ! TODO: if fillvalue = shr_const_spval the snowhl sent to the atm will have the spval over some points
       ! rather than 0 - this is very odd and needs to be understood
       !fldptr(:) = fillvalue
       fldptr(:) = 0._r8

       ! set fldptr values to input array
       if (present(minus)) then
          do g = bounds%begg, bounds%endg
             fldptr(g-bounds%begg+1) = - input(g)
          end do
       else
          do g = bounds%begg, bounds%endg
             fldptr(g-bounds%begg+1) = input(g)
          end do
       end if

       ! write debug output if appropriate
       if (masterproc .and. debug > 0 .and. get_nstep() < 5) then
          do g = bounds%begg,bounds%endg
             i = 1 + g - bounds%begg
             write(iulog,F01)'export: nstep, n, '//trim(fldname)//' = ',get_nstep(),i,input(g)
          end do
       end if

       ! check for nans
       call check_for_nans(fldptr, trim(fldname), bounds%begg)
    end if

  end subroutine state_setexport

  !===============================================================================

  subroutine state_getfldptr(State, fldname, fldptr, rc)
    ! ----------------------------------------------
    ! Get pointer to a state field
    ! ----------------------------------------------
    use ESMF , only : ESMF_State, ESMF_Field, ESMF_Mesh, ESMF_FieldStatus_Flag
    use ESMF , only : ESMF_StateGet, ESMF_FieldGet, ESMF_MeshGet
    use ESMF , only : ESMF_FIELDSTATUS_COMPLETE, ESMF_FAILURE

    type(ESMF_State),  intent(in)    :: State
    character(len=*),  intent(in)    :: fldname
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
    call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)

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
          call ESMF_LogWrite(trim(subname)//": no local nodes or elements ", ESMF_LOGMSG_INFO, rc=dbrc)
          rc = ESMF_FAILURE
          return
       end if

       call ESMF_FieldGet(lfield, farrayPtr=fldptr, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    endif  ! status

    call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)

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
             write(iulog,*) "NaN found in field ", trim(fname), ' at gridcell index ',begg+i-1
          end if
       end do
       call shr_sys_abort(' ERROR: One or more of the output from CLM to the coupler are NaN ' )
    end if
  end subroutine check_for_nans

end module lnd_import_export
