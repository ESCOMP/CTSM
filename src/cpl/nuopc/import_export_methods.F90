module lnd_import_export

  use ESMF                  , only : ESMF_State, ESMF_KIND_R8, ESMF_SUCCESS, ESMF_MAXSTR
  use ESMF                  , only : ESMF_LogWrite, ESMF_LOGMSG_ERROR, ESMF_LogFoundError
  use NUOPC                 , only : NUOPC_CompAttributeGet, NUOPC_Advertise
  use NUOPC_Model           , only : NUOPC_ModelGet
  use shr_kind_mod          , only : r8 => shr_kind_r8, cl=>shr_kind_cl
  use shr_infnan_mod        , only : isnan => shr_infnan_isnan
  use shr_string_mod        , only : shr_string_listGetName, shr_string_listGetNum
  use shr_sys_mod           , only : shr_sys_abort
  use clm_varctl            , only : iulog
  use decompmod             , only : bounds_type
  use lnd2atmType           , only : lnd2atm_type
  use lnd2glcMod            , only : lnd2glc_type
  use atm2lndType           , only : atm2lnd_type
  use glc2lndMod            , only : glc2lnd_type
  use domainMod             , only : ldomain
  use seq_drydep_mod        , only : seq_drydep_init, seq_drydep_readnl, lnd_drydep
  use shr_megan_mod         , only : shr_megan_readnl, shr_megan_mechcomps_n
  use shr_fire_emis_mod     , only : shr_fire_emis_readnl, shr_fire_emis_mechcomps_n, shr_fire_emis_ztop_token
  use shr_carma_mod         , only : shr_carma_readnl
  use shr_ndep_mod          , only : shr_ndep_readnl
  use glc_elevclass_mod     , only : glc_elevclass_init, glc_get_num_elevation_classes, glc_elevclass_as_string
  use shr_nuopc_scalars_mod , only : flds_scalar_name
  use shr_nuopc_scalars_mod , only : flds_scalar_num

  implicit none
  private ! except

  public :: advertise_fields
  public :: realize_fields
  public :: import_fields
  public :: export_fields

  private :: fld_list_add
  private :: fld_list_realize
  private :: check_for_nans

  type fld_list_type
    character(len=128) :: stdname
  end type fld_list_type

  integer, parameter     :: fldsMax = 100
  integer                :: fldsToLnd_num = 0
  integer                :: fldsFrLnd_num = 0
  type (fld_list_type)   :: fldsToLnd(fldsMax)
  type (fld_list_type)   :: fldsFrLnd(fldsMax)

  integer                :: glc_nec
  integer, parameter     :: dbug = 10
  integer, parameter     :: debug_packing = 0 ! internal debug level
  character(*),parameter :: F01 = "('(lnd_import_export) ',a,i4,2x,i5,2x,i8,2x,d21.14)"

!===============================================================================
contains
!===============================================================================

  subroutine advertise_fields(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    character(ESMF_MAXSTR) :: stdname
    character(ESMF_MAXSTR) :: cvalue
    logical                :: flds_co2a      ! use case
    logical                :: flds_co2b      ! use case
    logical                :: flds_co2c      ! use case

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
    glc_nec = glc_get_num_elevation_classes()
    if (glc_nec < 1) then
       call shr_sys_abort('ERROR: In CLM4.5 and later, glc_nec must be at least 1.')
    end if
    call glc_elevclass_init(glc_nec)

    !--------------------------------
    ! Export fields
    !--------------------------------

    call fld_list_add(fldsFrLnd_num, fldsFrlnd, trim(flds_scalar_name))

    ! land states
    call fld_list_add(fldsFrLnd_num, fldsFrlnd, 'Sl_lfrin'      )
    call fld_list_add(fldsFrLnd_num, fldsFrlnd, 'Sl_t'          )
    call fld_list_add(fldsFrLnd_num, fldsFrlnd, 'Sl_tref'       )
    call fld_list_add(fldsFrLnd_num, fldsFrlnd, 'Sl_qref'       )
    call fld_list_add(fldsFrLnd_num, fldsFrlnd, 'Sl_avsdr'      )
    call fld_list_add(fldsFrLnd_num, fldsFrlnd, 'Sl_anidr'      )
    call fld_list_add(fldsFrLnd_num, fldsFrlnd, 'Sl_avsdf'      )
    call fld_list_add(fldsFrLnd_num, fldsFrlnd, 'Sl_anidf'      )
    call fld_list_add(fldsFrLnd_num, fldsFrlnd, 'Sl_snowh'      )
    call fld_list_add(fldsFrLnd_num, fldsFrlnd, 'Sl_u10'        )
    call fld_list_add(fldsFrLnd_num, fldsFrlnd, 'Sl_ddvel'      )
    call fld_list_add(fldsFrLnd_num, fldsFrlnd, 'Sl_fv'         )
    call fld_list_add(fldsFrLnd_num, fldsFrlnd, 'Sl_ram1'       )
    call fld_list_add(fldsFrLnd_num, fldsFrlnd, 'Sl_ztopfire'   ) ! optional
    call fld_list_add(fldsFrLnd_num, fldsFrlnd, 'Sl_soilw'      ) ! optional

    ! fluxes to river
    call fld_list_add(fldsFrLnd_num, fldsFrlnd, 'Flrl_rofsur'   )
    call fld_list_add(fldsFrLnd_num, fldsFrlnd, 'Flrl_rofgwl'   )
    call fld_list_add(fldsFrLnd_num, fldsFrlnd, 'Flrl_rofsub'   )
    call fld_list_add(fldsFrLnd_num, fldsFrlnd, 'Flrl_rofi'     )
    call fld_list_add(fldsFrLnd_num, fldsFrlnd, 'Flrl_irrig'    )

    ! fluxes to atm
    call fld_list_add(fldsFrLnd_num, fldsFrlnd, 'Fall_taux'     )
    call fld_list_add(fldsFrLnd_num, fldsFrlnd, 'Fall_tauy'     )
    call fld_list_add(fldsFrLnd_num, fldsFrlnd, 'Fall_lat'      )
    call fld_list_add(fldsFrLnd_num, fldsFrlnd, 'Fall_sen'      )
    call fld_list_add(fldsFrLnd_num, fldsFrlnd, 'Fall_lwup'     )
    call fld_list_add(fldsFrLnd_num, fldsFrlnd, 'Fall_evap'     )
    call fld_list_add(fldsFrLnd_num, fldsFrlnd, 'Fall_swnet'    )
    call fld_list_add(fldsFrLnd_num, fldsFrlnd, 'Fall_fco2_lnd' )
    call fld_list_add(fldsFrLnd_num, fldsFrlnd, 'Fall_flxdst1'  )
    call fld_list_add(fldsFrLnd_num, fldsFrlnd, 'Fall_flxdst2'  )
    call fld_list_add(fldsFrLnd_num, fldsFrlnd, 'Fall_flxdst3'  )
    call fld_list_add(fldsFrLnd_num, fldsFrlnd, 'Fall_flxdst4'  )
    call fld_list_add(fldsFrLnd_num, fldsFrlnd, 'Fall_flxvoc'   )
    call fld_list_add(fldsFrLnd_num, fldsFrlnd, 'Fall_flxfire'  )
    call fld_list_add(fldsFrLnd_num, fldsFrlnd, 'Fall_fco2_lnd' )
   !call fld_list_add(fldsFrLnd_num, fldsFrlnd, 'Fall_methane'  )

    ! Dry Deposition fluxes (from land)
    call seq_drydep_readnl("drv_flds_in", mpicom, mastertask, drydep_fields, drydep_fields_token)
    if (drydep_fields_token /= '') then
       do n = 1, shr_string_listGetNum(fire_emis_fields)
          call  shr_string_listGetName(fire_emis_fields, n, fldname)
          call fld_list_add(fldsFrLnd_num, fldsFrLnd, trim(fldname))
       end do
    end if

    ! MEGAN emissions fluxes
    call shr_megan_readnl('drv_flds_in', mpicom, mastertask, megan_voc_fields, megan_fields_token)
    if (megan_fields_token /= '') then
       do n = 1, shr_string_listGetNum(megan_voc_fields)
          call  shr_string_listGetName(megan_voc_fields, n, fldname)
          call fld_list_add(fldsFrLnd_num, fldsFrLnd, trim(fldname))
       end do
    end if

    ! Fire emissions fluxes
    call shr_fire_emis_readnl('drv_flds_in', mpicom, mastertask, fire_emis_fields, &
         fire_emis_fields_token=fire_emis_fields_token, fire_emis_ztop_token=fire_emis_ztop_token)
    if (fire_emis_fields_token /= '') then
       do n = 1, shr_string_listGetNum(fire_emis_fields)
          call  shr_string_listGetName(fire_emis_fields, n, fldname)
          call fld_list_add(fldsFrLnd_num, fldsFrLnd, trim(fldname))
       end do
    end if
    if (fire_emis_ztop_token /= '') then
       call fld_list_add(fldsFrLnd_num, fldsFrLnd, trim(fire_emis_ztop_token))
    end if

    ! Create coupling fields for all lnd->glc elevation classes (1:glc_nec) plus bare land (index 0).
    do num = 0,glc_nec
       nec_str = glc_elevclass_as_string(num)
       name = 'Sl_tsrf' // nec_str
       call fld_list_add(fldsFrLnd_num, fldsFrLnd, trim(name))
       name = 'Sl_topo' // nec_str
       call fld_list_add(fldsFrLnd_num, fldsFrLnd, trim(name))
       name = 'Flgl_qice' // nec_str
       call fld_list_add(fldsFrLnd_num, fldsFrLnd, trim(name))
    end do

    ! Now advertise above export fields
    do n = 1,fldsFrLnd_num
       call NUOPC_Advertise(exportState, standardName=fldsFrLnd(n)%stdname, &
            TransferOfferGeomObject='will provide', rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    enddo

    !--------------------------------
    ! Import
    !--------------------------------

    call fld_list_add(fldsToLnd_num, fldsToLnd, trim(flds_scalar_name))

    ! from atm
    call fld_list_add(fldsToLnd_num, fldsToLnd, 'Sa_z'         )
    call fld_list_add(fldsToLnd_num, fldsToLnd, 'Sa_topo'      )
    call fld_list_add(fldsToLnd_num, fldsToLnd, 'Sa_u'         )
    call fld_list_add(fldsToLnd_num, fldsToLnd, 'Sa_v'         )
    call fld_list_add(fldsToLnd_num, fldsToLnd, 'Sa_ptem'      )
    call fld_list_add(fldsToLnd_num, fldsToLnd, 'Sa_pbot'      )
    call fld_list_add(fldsToLnd_num, fldsToLnd, 'Sa_tbot'      )
    call fld_list_add(fldsToLnd_num, fldsToLnd, 'Sa_shum'      )
   !call fld_list_add(fldsToLnd_num, fldsToLnd, 'Sa_methane'   )

    call fld_list_add(fldsToLnd_num, fldsToLnd, 'Flrr_volr'    )
    call fld_list_add(fldsToLnd_num, fldsToLnd, 'Flrr_volrmch' )
    call fld_list_add(fldsToLnd_num, fldsToLnd, 'Faxa_lwdn'    )
    call fld_list_add(fldsToLnd_num, fldsToLnd, 'Faxa_rainc'   )
    call fld_list_add(fldsToLnd_num, fldsToLnd, 'Faxa_rainl'   )
    call fld_list_add(fldsToLnd_num, fldsToLnd, 'Faxa_snowc'   )
    call fld_list_add(fldsToLnd_num, fldsToLnd, 'Faxa_snowl'   )
    call fld_list_add(fldsToLnd_num, fldsToLnd, 'Faxa_swndr'   )
    call fld_list_add(fldsToLnd_num, fldsToLnd, 'Faxa_swvdr'   )
    call fld_list_add(fldsToLnd_num, fldsToLnd, 'Faxa_swndf'   )
    call fld_list_add(fldsToLnd_num, fldsToLnd, 'Faxa_swvdf'   )
    call fld_list_add(fldsToLnd_num, fldsToLnd, 'Faxa_bcphidry')
    call fld_list_add(fldsToLnd_num, fldsToLnd, 'Faxa_bcphodry')
    call fld_list_add(fldsToLnd_num, fldsToLnd, 'Faxa_bcphiwet')
    call fld_list_add(fldsToLnd_num, fldsToLnd, 'Faxa_ocphidry')
    call fld_list_add(fldsToLnd_num, fldsToLnd, 'Faxa_ocphodry')
    call fld_list_add(fldsToLnd_num, fldsToLnd, 'Faxa_ocphiwet')
    call fld_list_add(fldsToLnd_num, fldsToLnd, 'Faxa_dstdry1' )
    call fld_list_add(fldsToLnd_num, fldsToLnd, 'Faxa_dstdry2' )
    call fld_list_add(fldsToLnd_num, fldsToLnd, 'Faxa_dstdry3' )
    call fld_list_add(fldsToLnd_num, fldsToLnd, 'Faxa_dstdry4' )
    call fld_list_add(fldsToLnd_num, fldsToLnd, 'Faxa_dstwet1' )
    call fld_list_add(fldsToLnd_num, fldsToLnd, 'Faxa_dstwet2' )
    call fld_list_add(fldsToLnd_num, fldsToLnd, 'Faxa_dstwet3' )
    call fld_list_add(fldsToLnd_num, fldsToLnd, 'Faxa_dstwet4' )

    ! from atm - co2 exchange scenarios
    if (flds_co2a .or. flds_co2b .or. flds_co2c) then
       call fld_list_add(fldsToLnd_num, fldsToLnd, 'Sa_co2prog')
       call fld_list_add(fldsToLnd_num, fldsToLnd, 'Sa_co2diag')
    end if

    ! from atm - nitrogen deposition
    call shr_ndep_readnl(nlfilename="drv_flds_in", id=compid, &
         ndep_fields=ndep_fields, add_ndep_fields=add_ndep_fields)
    if (add_ndep_fields) then
       do n = 1, shr_string_listGetNum(ndep_fields)
          call  shr_string_listGetName(ndep_fields, n, fldname)
          call fld_list_add(fldsToLnd_num, fldsToLnd, trim(fldname))
       end do
       ! This sets a variable in clm_varctl
       ndep_from_cpl = .true.
    end if

    ! from river
    call fld_list_add(fldsToLnd_num, fldsToLnd, 'Flrr_flood')

    ! from glc
    call fld_list_add(fldsToLnd_num, fldsToLnd, 'Sg_icemask'               )
    call fld_list_add(fldsToLnd_num, fldsToLnd, 'Sg_icemask_coupled_fluxes')
    do num = 0,glc_nec
       ! create coupling fields for all glc->lnd elevation classes (1:glc_nec) plus bare land (index 0).
       nec_str = glc_elevclass_as_string(num)
       name = 'Sg_ice_covered' // nec_str
       call fld_list_add(fldsToLnd_num, fldsToLnd, trim(name))
       name = 'Sg_topo' // nec_str
       call fld_list_add(fldsToLnd_num, fldsToLnd, trim(name))
       name = 'Flgg_hflx' // nec_str
       call fld_list_add(fldsToLnd_num, fldsToLnd, trim(name))
    end do

    ! Now advertise import fields
    do n = 1,fldsToAtm_num
       call NUOPC_Advertise(importState, standardName=fldsToLnd(n)%stdname, &
            TransferOfferGeomObject='will provide', rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    enddo

  end subroutine advertise_fields

  !===============================================================================

  subroutine realize_fields(importState, exportState, Emesh, rc)
    type(ESMF_State)     :: importState
    type(ESMF_State)     :: exportState
    type(ESMF_Mesh)      :: Emesh
    integer, intent(out) :: rc

    rc = ESMF_SUCCESS
    !---------------------------------------------------------------------------

    call fld_list_realize( &
         state=ExportState, &
         fldList=fldsFrLnd, &
         numflds=fldsFrLnd_num, &
         flds_scalar_name=flds_scalar_name, &
         flds_scalar_num=flds_scalar_num, &
         tag=subname//':clmExport',&
         mesh=Emesh, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call fld_list_realize( &
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

  subroutine import_fields( gcomp, bounds, glc_present, atm2lnd_inst, glc2lnd_inst)

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
    type(ESMF_State)  :: importState, exportState
    real(r8), pointer :: dataPtr(:)
    real(r8), pointer :: dataPtr1(:)
    real(r8), pointer :: dataPtr2(:)
    integer           :: begg, endg                          ! bounds
    integer           :: g,i,k,nstep,ier                     ! indices, number of steps, and error code
    real(r8)          :: e                                   ! vapor pressure (Pa)
    real(r8)          :: qsat                                ! saturation specific humidity (kg/kg)
    real(r8)          :: co2_ppmv_diag                       ! temporary
    real(r8)          :: co2_ppmv_prog                       ! temporary
    real(r8)          :: co2_ppmv_val                        ! temporary
    real(r8)          :: esatw                               ! saturation vapor pressure over water (Pa)
    real(r8)          :: esati                               ! saturation vapor pressure over ice (Pa)
    real(r8)          :: a0,a1,a2,a3,a4,a5,a6                ! coefficients for esat over water
    real(r8)          :: b0,b1,b2,b3,b4,b5,b6                ! coefficients for esat over ice
    real(r8)          :: tdc, t                              ! Kelvins to Celcius function and its input
    real(r8)          :: forc_t                              ! atmospheric temperature (Kelvin)
    real(r8)          :: forc_q                              ! atmospheric specific humidity (kg/kg)
    real(r8)          :: forc_pbot                           ! atmospheric pressure (Pa)
    real(r8)          :: forc_rainc(bounds%begg:bounds:endg) ! rainxy Atm flux mm/s
    real(r8)          :: forc_rainl(bounds%begg:bounds:endg) ! rainxy Atm flux mm/s
    real(r8)          :: forc_snowc(bounds%begg:bounds:endg) ! snowfxy Atm flux  mm/s
    real(r8)          :: forc_snowl(bounds%begg:bounds:endg) ! snowfxl Atm flux  mm/s
    real(r8)          :: frac_grc(bounds%begg:bounds:endg, glc_nec)
    real(r8)          :: topo_grc(bounds%begg:bounds:endg, glc_nec)
    real(r8)          :: hflx_grc(bounds%begg:bounds:endg, glc_nec)
    real(r8)          :: glc2lnd_icemask(bounds%begg:bounds:endg)
    real(r8)          :: glc2lnd_icemask_coupled_fluxes(bounds%begg:bounds:endg)
    character(len=32), parameter :: sub = 'nuopc_import'

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
    call NUOPC_ModelGet(gcomp, modelClock=clock, importState=importState, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Error checks
    if (NUOPC_IsConnected(exportState, 'Sa_co2prog')) then
       if (co2_type == 'prognostic') then
          call shr_sys_abort( sub//' ERROR: must have nonzero Sa_co2prog for co2_type equal to prognostic' )
       end if
    end if
    if (NUOPC_IsConnected(exportState, 'Sa_co2diag')) then
       if (co2_type == 'diagnostic') then
          call shr_sys_abort( sub//' ERROR: must have nonzero Sa_co2diag for co2_type equal to diagnostic' )
       end if
    end if

    ! Set bounds
    begg = bounds%begg; endg=bounds%endg

    ! Note: precipitation fluxes received  from the coupler
    ! are in units of kg/s/m^2. To convert these precipitation rates
    ! in units of mm/sec, one must divide by 1000 kg/m^3 and multiply
    ! by 1000 mm/m resulting in an overall factor of unity.
    ! Below the units are therefore given in mm/s.

    ! Determine flooding input, sign convention is positive downward and
    ! hierarchy is atm/glc/lnd/rof/ice/ocn.  so water sent from rof to land is negative,
    ! change the sign to indicate addition of water to system.

    call State_getFldPtr(importState, 'Flrr_flood', dataPtr, rc=rc )
    do g = begg,endg
       i = 1 + (g-begg)
       atm2lnd_inst%forc_flood_grc(g) = -dataPtr(i)
    end do
    call check_for_nans(dataPtr, 'Flrr_flood', begg)

    call State_getFldPtr(importState, 'Flrr_volr', dataPtr, rc=rc )
    do g = begg,endg
       i = 1 + (g-begg)
       atm2lnd_inst%volr_grc(g) = dataPtr(i) * (ldomain%area(g) * 1.e6_r8)
    end do
    call check_for_nans(dataPtr, 'Flrr_volr', begg)

    call State_getFldPtr(importState, 'Flrr_volrmch', dataPtr, rc=rc )
    do g = begg,endg
       i = 1 + (g-begg)
       atm2lnd_inst%volrmch_grc(g) = dataPtr(i) * (ldomain%area(g) * 1.e6_r8)
    end do
    call check_for_nans(dataPtr, '', begg)

    ! Determine required input fields

    call State_getFldPtr(importState, 'Sa_z', dataPtr, rc=rc )
    do g = begg,endg
       i = 1 + (g-begg)
       atm2lnd_inst%forc_hgt_grc(g) = dataPtr(i)
    end do
    call check_for_nans(dataPtr, 'Sa_z', begg)

    call State_getFldPtr(importState, 'Sa_u', dataPtr, rc=rc )
    do g = begg,endg
       i = 1 + (g-begg)
       atm2lnd_inst%forc_u_grc(g) = dataPtr(i)
    end do
    call check_for_nans(dataPtr, 'Sa_u', begg)

    call State_getFldPtr(importState, 'Sa_v', dataPtr, rc=rc )
    do g = begg,endg
       i = 1 + (g-begg)
       atm2lnd_inst%forc_v_grc(g) = dataPtr(i)
    end do
    call check_for_nans(dataPtr, 'Sa_v', begg)

    call State_getFldPtr(importState, 'Sa_ptem', dataPtr, rc=rc )
    do g = begg,endg
       i = 1 + (g-begg)
       atm2lnd_inst%forc_th_not_downscaled_grc(g) = dataPtr(i)
    end do
    call check_for_nans(dataPtr, 'Sa_ptem', begg)

    call State_getFldPtr(importState, 'Sa_shum', dataPtr, rc=rc )
    do g = begg,endg
       i = 1 + (g-begg)
       atm2lnd_inst%forc_q_not_downscaled_grc(g) = dataPtr(i)
    end do
    call check_for_nans(dataPtr, 'Sa_shum', begg)

    call State_getFldPtr(importState, 'Sa_pbot', dataPtr, rc=rc )
    do g = begg,endg
       i = 1 + (g-begg)
       atm2lnd_inst%forc_pbot_not_downscaled_grc(g) = dataPtr(i)
    end do
    call check_for_nans(dataPtr, 'Sa_pbot', begg)

    call State_getFldPtr(importState, 'Sa_tbot', dataPtr, rc=rc )
    do g = begg,endg
       i = 1 + (g-begg)
       atm2lnd_inst%forc_t_not_downscaled_grc(g) = dataPtr(i)
    end do
    call check_for_nans(dataPtr, 'Sa_tbot', begg)

    call State_getFldPtr(importState, 'Faxa_rainc', dataPtr, rc=rc )
    do g = begg,endg
       i = 1 + (g-begg)
       forc_rainc(g) = dataPtr(i)
    end do
    call check_for_nans(dataPtr, 'Faxa_rainc', begg)

    call State_getFldPtr(importState, 'Faxa_snowc', dataPtr, rc=rc )
    do g = begg,endg
       i = 1 + (g-begg)
       forc_snowc(g) = dataPtr(i)
    end do
    call check_for_nans(dataPtr, 'Faxa_snowc', begg)

    call State_getFldPtr(importState, 'Faxa_rainl', dataPtr, rc=rc )
    do g = begg,endg
       i = 1 + (g-begg)
       forc_rainl(g) = dataPtr(i)
    end do
    call check_for_nans(dataPtr, 'Faxa_rainl', begg)

    call State_getFldPtr(importState, 'Faxa_snowl', dataPtr, rc=rc )
    do g = begg,endg
       i = 1 + (g-begg)
       forc_snowl(g) = dataPtr(i)
    end do
    call check_for_nans(dataPtr, 'Faxa_snowl', begg)

    call State_getFldPtr(importState, 'Faxa_lwdn', dataPtr, rc=rc )
    do g = begg,endg
       i = 1 + (g-begg)
       atm2lnd_inst%forc_lwrad_not_downscaled_grc(g) = dataPtr(i)
    end do
    call check_for_nans(dataPtr, 'Faxa_lwdn', begg)

    call State_getFldPtr(importState, 'Faxa_swndr', dataPtr, rc=rc )
    do g = begg,endg
       i = 1 + (g-begg)
       atm2lnd_inst%forc_solad_grc(g,2) = dataPtr(i)
    end do
    call check_for_nans(dataPtr, 'Faxa_swndr', begg)

    call State_getFldPtr(importState, 'Faxa_swvdr', dataPtr, rc=rc )
    do g = begg,endg
       i = 1 + (g-begg)
       atm2lnd_inst%forc_solad_grc(g,1) = dataPtr(i)
    end do
    call check_for_nans(dataPtr, 'Faxa_swvdr', begg)

    call State_getFldPtr(importState, 'Faxa_swndf', dataPtr, rc=rc )
    do g = begg,endg
       i = 1 + (g-begg)
       atm2lnd_inst%forc_solai_grc(g,2) = dataPtr(i)
    end do
    call check_for_nans(dataPtr, 'Faxa_swndf', begg)

    call State_getFldPtr(importState, 'Faxa_swvdf', dataPtr, rc=rc )
    do g = begg,endg
       i = 1 + (g-begg)
       atm2lnd_inst%forc_solai_grc(g,1) = dataPtr(i)
    end do
    call check_for_nans(dataPtr, 'Faxa_swvdf', begg)

    if (NUOPC_IsConnected(importState, fieldName='Faxa_bcphidry')) then
       call State_getFldPtr(exportState,'Faxa_bcphidry', dataPtr, rc=rc )
       do g = begg,endg
          i = 1 + (g-begg)
          atm2lnd_inst%forc_aer_grc(g,1) = dataPtr(i)
       end do
       call check_for_nans(dataPtr, 'Faxa_bcphidry', begg)
    end if

    ! atmosphere coupling, for prognostic/prescribed aerosols

    if (NUOPC_IsConnected(importState, fieldName='Faxa_bcphodry')) then
       call State_getFldPtr(exportState,'Faxa_bcphodry', dataPtr, rc=rc )
       do g = begg,endg
          i = 1 + (g-begg)
          atm2lnd_inst%forc_aer_grc(g,2) = dataPtr(i)
       end do
       call check_for_nans(dataPtr, 'Faxa_bcphodry', begg)
    end if

    if (NUOPC_IsConnected(importState, fieldName='Faxa_bcphiwet')) then
       call State_getFldPtr(exportState,'Faxa_bcphiwet', dataPtr, rc=rc )
       do g = begg,endg
          i = 1 + (g-begg)
          atm2lnd_inst%forc_aer_grc(g,3) = dataPtr(i)
       end do
       call check_for_nans(dataPtr, 'Faxa_bcphiwet', begg)
    end if

    if (NUOPC_IsConnected(importState, fieldName='Faxa_ocphidry')) then
       call State_getFldPtr(exportState,'Faxa_ocphidry', dataPtr, rc=rc )
       do g = begg,endg
          i = 1 + (g-begg)
          atm2lnd_inst%forc_aer_grc(g,4) = dataPtr(i)
       end do
       call check_for_nans(dataPtr, 'Faxa_ocphidry', begg)
    end if

    if (NUOPC_IsConnected(importState, fieldName='Faxa_ocphiodry')) then
       call State_getFldPtr(exportState,'Faxa_ocphiodry', dataPtr, rc=rc )
       do g = begg,endg
          i = 1 + (g-begg)
          atm2lnd_inst%forc_aer_grc(g,5) = dataPtr(i)
       end do
       call check_for_nans(dataPtr, 'Faxa_ocphiodry', begg)
    end if

    if (NUOPC_IsConnected(importState, fieldName='Faxa_ocphiwet')) then
       call State_getFldPtr(exportState,'Faxa_ocphiwet', dataPtr, rc=rc )
       do g = begg,endg
          i = 1 + (g-begg)
          atm2lnd_inst%forc_aer_grc(g,6) = dataPtr(i)
       end do
       call check_for_nans(dataPtr, 'Faxa_ocphiwet', begg)
    end if

    if (NUOPC_IsConnected(importState, fieldName='Faxa_wet1')) then
       call State_getFldPtr(exportState,'Faxa_wet1', dataPtr, rc=rc )
       do g = begg,endg
          i = 1 + (g-begg)
          atm2lnd_inst%forc_aer_grc(g,7) = dataPtr(i)
       end do
       call check_for_nans(dataPtr, 'Faxa_wet1', begg)
    end if

    if (NUOPC_IsConnected(importState, fieldName='Faxa_dry1')) then
       call State_getFldPtr(exportState,'Faxa_dry1', dataPtr, rc=rc )
       do g = begg,endg
          i = 1 + (g-begg)
          atm2lnd_inst%forc_aer_grc(g,8) = dataPtr(i)
       end do
       call check_for_nans(dataPtr, 'Faxa_dry1', begg)
    end if

    if (NUOPC_IsConnected(importState, fieldName='Faxa_wet2')) then
       call State_getFldPtr(exportState,'Faxa_wet2', dataPtr, rc=rc )
       do g = begg,endg
          i = 1 + (g-begg)
          atm2lnd_inst%forc_aer_grc(g,9) = dataPtr(i)
       end do
       call check_for_nans(dataPtr, 'Faxa_wet2', begg)
    end if

    if (NUOPC_IsConnected(importState, fieldName='Faxa_dry2')) then
       call State_getFldPtr(exportState,'Faxa_dry2', dataPtr, rc=rc )
       do g = begg,endg
          i = 1 + (g-begg)
          atm2lnd_inst%forc_aer_grc(g,10) = dataPtr(i)
       end do
       call check_for_nans(dataPtr, 'Faxa_dry2', begg)
    end if

    if (NUOPC_IsConnected(importState, fieldName='Faxa_wet3')) then
       call State_getFldPtr(exportState,'Faxa_wet3', dataPtr, rc=rc )
       do g = begg,endg
          i = 1 + (g-begg)
          atm2lnd_inst%forc_aer_grc(g,11) = dataPtr(i)
       end do
       call check_for_nans(dataPtr, 'Faxa_wet3', begg)
    end if

    if (NUOPC_IsConnected(importState, fieldName='Faxa_dry3')) then
       call State_getFldPtr(exportState,'Faxa_dry3', dataPtr, rc=rc )
       do g = begg,endg
          i = 1 + (g-begg)
          atm2lnd_inst%forc_aer_grc(g,12) = dataPtr(i)
       end do
       call check_for_nans(dataPtr, 'Faxa_dry3', begg)
    end if

    if (NUOPC_IsConnected(importState, fieldName='Faxa_wet4')) then
       call State_getFldPtr(exportState,'Faxa_wet4', dataPtr, rc=rc )
       do g = begg,endg
          i = 1 + (g-begg)
          atm2lnd_inst%forc_aer_grc(g,13) = dataPtr(i)
       end do
       call check_for_nans(dataPtr, 'Faxa_wet4', begg)
    end if

    if (NUOPC_IsConnected(importState, fieldName='Faxa_dry4')) then
       call State_getFldPtr(exportState,'Faxa_dry4', dataPtr, rc=rc )
       do g = begg,endg
          i = 1 + (g-begg)
          atm2lnd_inst%forc_aer_grc(g,14) = dataPtr(i)
       end do
       call check_for_nans(dataPtr, 'Faxa_dry4', begg)
    end if

    if ( NUOPC_IsConnected(importState, fieldName='Sa_co2prog')) then
       call State_getFldPtr(exportState,'Sa_co2prog', dataPtr, rc=rc )
       co2_ppmv_prog = dataPtr(1))   ! co2 atm state prognostic
    else
       co2_ppmv_prog = co2_ppmv
    end if

    if ( NUOPC_IsConnected(importState, fieldName='Sa_co2diag')) then
       call State_getFldPtr(exportState,'Sa_co2diag', dataPtr, rc=rc )
       co2_ppmv_diag = dataPtr(1))   ! co2 atm state diagnostic
    else
       co2_ppmv_diag = co2_ppmv
    end if

    if ( NUOPC_IsConnected(importState, fieldName='Sa_methane')) then
       call State_getFldPtr(exportState,'Sa_methane', dataPtr, rc=rc )
       do g = begg,endg
          i = 1 + (g-begg)
          atm2lnd_inst%forc_pch4_grc(g) = dataPtr(i)
       end do
    end if

     ! Determine derived quantities for required fields

     do g = begg, endg

       forc_t    = atm2lnd_inst%forc_t_not_downscaled_grc(g)
       forc_q    = atm2lnd_inst%forc_q_not_downscaled_grc(g)
       forc_pbot = atm2lnd_inst%forc_pbot_not_downscaled_grc(g)

       atm2lnd_inst%forc_hgt_u_grc(g) = atm2lnd_inst%forc_hgt_grc(g)    !observational height of wind [m]
       atm2lnd_inst%forc_hgt_t_grc(g) = atm2lnd_inst%forc_hgt_grc(g)    !observational height of temperature [m]
       atm2lnd_inst%forc_hgt_q_grc(g) = atm2lnd_inst%forc_hgt_grc(g)    !observational height of humidity [m]

       atm2lnd_inst%forc_vp_grc(g)    = forc_q * forc_pbot  / (0.622_r8 + 0.378_r8 * forc_q)

       atm2lnd_inst%forc_rho_not_downscaled_grc(g) = &
            (forc_pbot - 0.378_r8 * atm2lnd_inst%forc_vp_grc(g)) / (rair * forc_t)

       atm2lnd_inst%forc_po2_grc(g)   = o2_molar_const * forc_pbot

       atm2lnd_inst%forc_wind_grc(g)  = sqrt(atm2lnd_inst%forc_u_grc(g)**2 + atm2lnd_inst%forc_v_grc(g)**2)

       atm2lnd_inst%forc_solar_grc(g) = atm2lnd_inst%forc_solad_grc(g,1) + atm2lnd_inst%forc_solai_grc(g,1) + &
                                        atm2lnd_inst%forc_solad_grc(g,2) + atm2lnd_inst%forc_solai_grc(g,2)

       atm2lnd_inst%forc_rain_not_downscaled_grc(g)  = forc_rainc(g) + forc_rainl(g)

       atm2lnd_inst%forc_snow_not_downscaled_grc(g)  = forc_snowc(g) + forc_snowl(g)

       if (forc_t > SHR_CONST_TKFRZ) then
          e = esatw(tdc(forc_t))
       else
          e = esati(tdc(forc_t))
       end if
       qsat           = 0.622_r8*e / (forc_pbot - 0.378_r8*e)

       !modify specific humidity if precip occurs
       if (1==2) then
          if((forc_rainc+forc_rainl) > 0._r8) then
             forc_q = 0.95_r8*qsat
             !           forc_q = qsat
             atm2lnd_inst%forc_q_not_downscaled_grc(g) = forc_q
          endif
       endif

       atm2lnd_inst%forc_rh_grc(g) = 100.0_r8*(forc_q / qsat)
    end do

    ! Optional fields

    if ( NUOPC_IsConnected(importState, fieldName='Faxa_nhx') .and. &
         NUOPC_IsConnected(importState, fieldName='Faxa_noy')) then
       ! The mediator is sending ndep in units if kgN/m2/s - and clm
       ! uses units of gN/m2/sec - so the following conversion needs to happen
       call State_getFldPtr(importState, 'Faxa_nhx', dataPtr1, rc=rc )
       call State_getFldPtr(importState, 'Faxa_noy', dataPtr2, rc=rc )
       atm2lnd_inst%forc_ndep_grc(g) = (dataPtr1(i) + dataPtr2(i))*1000._r8
    end if

    ! Check that solar, specific-humidity and LW downward aren't negative
    do g = begg,endg
       if ( atm2lnd_inst%forc_lwrad_not_downscaled_grc(g) <= 0.0_r8 )then
          call shr_sys_abort( sub//' ERROR: Longwave down sent from the atmosphere model is negative or zero' )
       end if
       if ( (atm2lnd_inst%forc_solad_grc(g,1) < 0.0_r8) .or. &
            (atm2lnd_inst%forc_solad_grc(g,2) < 0.0_r8) .or. &
            (atm2lnd_inst%forc_solai_grc(g,1) < 0.0_r8) .or. &
            (atm2lnd_inst%forc_solai_grc(g,2) < 0.0_r8) ) then
          call shr_sys_abort( sub//&
               ' ERROR: One of the solar fields (indirect/diffuse, vis or near-IR)'// &
               ' from the atmosphere model is negative or zero' )
       end if
       if ( atm2lnd_inst%forc_q_not_downscaled_grc(g) < 0.0_r8 )then
          call shr_sys_abort( sub//&
               ' ERROR: Bottom layer specific humidty sent from the atmosphere model is less than zero' )
       end if
    end do

    ! Note that the following does unit conversions from ppmv to partial pressures (Pa)
    ! Note that forc_pbot is in Pa

    if (co2_type == 'prognostic') then
       co2_ppmv_val = co2_ppmv_prog
    else if (co2_type == 'diagnostic') then
       co2_ppmv_val = co2_ppmv_diag
    else
       co2_ppmv_val = co2_ppmv
    end if
    do begg, endg
       atm2lnd_inst%forc_pco2_grc(g)  = co2_ppmv_val * 1.e-6_r8 * forc_pbot
       if (use_c13) then
          atm2lnd_inst%forc_pc13o2_grc(g) = co2_ppmv_val * c13ratio * 1.e-6_r8 * forc_pbot
       end if
    end do

    ! Make sure relative humidity is properly bounded
    ! atm2lnd_inst%forc_rh_grc(g) = min( 100.0_r8, atm2lnd_inst%forc_rh_grc(g) )
    ! atm2lnd_inst%forc_rh_grc(g) = max(   0.0_r8, atm2lnd_inst%forc_rh_grc(g) )

    ! glc fields
    ! We could avoid setting these fields if glc_present is .false., if that would
    ! help with performance. (The downside would be that we wouldn't have these fields
    ! available for diagnostic purposes or to force a later T compset with dlnd.)
    do num = 0,glc_nec
       nec_str = glc_elevclass_as_string(num)

       name = 'Sg_ice_covered' // nec_str
       call State_getFldPtr(importState,trim(name), dataPtr, rc=rc )
       do g = begg,endg
          i = 1 + (g-begg)
          frac_grc(g,num) = dataPtr(i)
       end do
       call check_for_nans(dataPtr, trim(name), begg)

       name = 'Sg_topo' // nec_str
       call State_getFldPtr(importState,trim(name), dataPtr, rc=rc )
       do g = begg,endg
          i = 1 + (g-begg)
          topo_grc(g,num) = dataPtr(i)
       end do
       call check_for_nans(dataPtr, trim(name), begg)

       name = 'Flgg_hflx' // nec_str
       call State_getFldPtr(importState,trim(name), dataPtr, rc=rc )
       do g = begg,endg
          i = 1 + (g-begg)
          hflx_grc(g,num) = dataPtr(i)
       end do
       call check_for_nans(dataPtr, trim(name), begg)
    end do

    call State_getFldPtr(exportState,'Sg_icemask', dataPtr, rc=rc )
    do g = begg,endg
       i = 1 + (g-begg)
       glc2lnd_icemask(g) = dataPtr(i)
    end do
    call check_for_nans(dataPtr, 'Sg_icemask', begg)

    call State_getFldPtr(exportState,'Sg_icemask_coupled_fluxes', dataPtr, rc=rc )
    do g = begg,endg
       i = 1 + (g-begg)
       glc2lnd_icemask_coupled_fluxes(g) = dataPtr(i)
    end do
    call check_for_nans(dataPtr, 'Sg_icemask_coupled_fluxes', begg)

    call glc2lnd_inst%set_glc2lnd_fields( &
         bounds, glc_present, frac_grc, topo_grc, hflx_grc, &
         glc2lnd_ice_mask, glc2lnd_icemask_coupled_fluxes)

    ! Write debug output if appropriate
    if (masterproc .and. debug_import > 0 .and. get_nstep() < 5) then
       do g = bounds%begg,bounds%endg
          i = 1 + (g - bounds%begg)
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
          write(iulog,F01)'import: nstep, n, Fa_snowl    = ',get_nstep(), i, forc_snowc(l)
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
       end do
    end if

  end subroutine import_fields

  !===============================================================================

  subroutine export_fields( export_state, bounds, lnd2atm_inst, lnd2glc_inst, rc)

    !-------------------------------
    ! Pact the export state
    !-------------------------------

    ! input/output variabes
    type(ESMF_GridComp)               :: gcomp
    type(bounds_type) , intent(in)    :: bounds       ! bounds
    type(lnd2atm_type), intent(inout) :: lnd2atm_inst ! land to atmosphere exchange data type
    type(lnd2glc_type), intent(inout) :: lnd2glc_inst ! land to atmosphere exchange data type
    integer           , intent(out)   :: rc

    ! local variables
    type(ESMF_State)  :: export_state
    integer           :: i, g, num
    integer           :: begg, endg
    real(r8), pointer :: dataPtr(:)
    !---------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    begg = bounds%begg; endg=bounds%endg

    call State_getFldPtr(exportState, 'Sl_t', dataPtr, rc=rc )
    do g = begg,endg
       i = 1 + (g-begg)
       dataPtr(i) = lnd2atm_inst%t_rad_grc(g)
    end do
    call check_for_nans(dataPtr, 'Sl_t', begg)

    call State_getFldPtr(exportState, 'Sl_showh', dataPtr, rc=rc )
    do g = begg,endg
       i = 1 + (g-begg)
       dataPtr_Sl_snowh(i)= lnd2atm_inst%h2osno_grc(g)
    end do
    call check_for_nans(dataPtr, 'Sl_snowh', begg)

    call State_getFldPtr(exportState, 'Sl_avsdr', dataPtr, rc=rc )
    do g = begg,endg
       i = 1 + (g-begg)
       dataPtr(i) = lnd2atm_inst%albd_grc(g,1)
    end do
    call check_for_nans(dataPtr, 'Sl_avsdr', begg)

    call State_getFldPtr(exportState, 'Sl_anidr', dataPtr, rc=rc )
    do g = begg,endg
       i = 1 + (g-begg)
       dataPtr(i) =  lnd2atm_inst%albd_grc(g,2)
    end do
    call check_for_nans(dataPtr, 'Sl_anidr', begg)

    call State_getFldPtr(exportState, 'Sl_avsdf', dataPtr, rc=rc )
    do g = begg,endg
       i = 1 + (g-begg)
       dataPtr(i) =  lnd2atm_inst%albi_grc(g,1)
    end do
    call check_for_nans(dataPtr, 'Sl_avsdf', begg)

    call State_getFldPtr(exportState, 'Sl_anidf', dataPtr, rc=rc )
    do g = begg,endg
       i = 1 + (g-begg)
       dataPtr(i) =  lnd2atm_inst%albi_grc(g,2)
    end do
    call check_for_nans(dataPtr, 'Sl_anidf', begg)

    call State_getFldPtr(exportState, 'Sl_tref', dataPtr, rc=rc )
    do g = begg,endg
       i = 1 + (g-begg)
       dataPtr(i) =  lnd2atm_inst%t_ref2m_grc(g)
    end do
    call check_for_nans(dataPtr, 'Sl_tref', begg)

    call State_getFldPtr(exportState, 'Sl_qref', dataPtr, rc=rc )
    do g = begg,endg
       i = 1 + (g-begg)
       dataPtr(i) =  lnd2atm_inst%q_ref2m_grc(g)
    end do
    call check_for_nans(dataPtr, 'Sl_qref', begg)

    call State_getFldPtr(exportState, 'Sl_u10', dataPtr, rc=rc )
    do g = begg,endg
       i = 1 + (g-begg)
       dataPtr(i) =  lnd2atm_inst%u_ref10m_grc(g)
    end do
    call check_for_nans(dataPtr, 'Sl_u10', begg)

    call State_getFldPtr(exportState, 'Fall_taux', dataPtr, rc=rc )
    do g = begg,endg
       i = 1 + (g-begg)
       dataPtr(i) = -lnd2atm_inst%taux_grc(g)
    end do
    call check_for_nans(dataPtr, 'Fall_taux', begg)

    call State_getFldPtr(exportState, 'Fall_tauy', dataPtr, rc=rc )
    do g = begg,endg
       i = 1 + (g-begg)
       dataPtr(i) = -lnd2atm_inst%tauy_grc(g)
    end do
    call check_for_nans(dataPtr, 'Fall_tauy', begg)

    call State_getFldPtr(exportState, 'Fall_lat', dataPtr, rc=rc )
    do g = begg,endg
       i = 1 + (g-begg)
       dataPtr(i) = -lnd2atm_inst%eflx_lh_tot_grc(g)
    end do
    call check_for_nans(dataPtr, 'Fall_lat', begg)

    call State_getFldPtr(exportState, 'Fall_sen', dataPtr, rc=rc )
    do g = begg,endg
       i = 1 + (g-begg)
       dataPtr(i) = -lnd2atm_inst%eflx_sh_tot_grc(g)
    end do
    call check_for_nans(dataPtr, 'Fall_sen', begg)

    call State_getFldPtr(exportState, 'Fall_lwup', dataPtr, rc=rc )
    do g = begg,endg
       i = 1 + (g-begg)
       dataPtr(i) = -lnd2atm_inst%eflx_lwrad_out_grc(g)
    end do
    call check_for_nans(dataPtr, 'Fall_lwup', begg)

    call State_getFldPtr(exportState, 'Fall_evap', dataPtr, rc=rc )
    do g = begg,endg
       i = 1 + (g-begg)
       dataPtr(i) = -lnd2atm_inst%qflx_evap_tot_grc(g)
    end do
    call check_for_nans(dataPtr, 'Fall_evap', begg)

    call State_getFldPtr(exportState, 'Fall_swnet', dataPtr, rc=rc )
    do g = begg,endg
       i = 1 + (g-begg)
       dataPtr(i) =  lnd2atm_inst%fsa_grc(g)
    end do
    call check_for_nans(dataPtr, 'Fall_swnet', begg)

    if (NUOPC_IsConnected(exportState, fieldName='Sl_ram1')) then
       call State_getFldPtr(exportState,'Sl_ram1', dataPtr, rc=rc )
       do g = begg,endg
          i = 1 + (g-begg)
          dataPtr(i) = lnd2atm_inst%ram1_grc(g)
       end do
       call check_for_nans(dataPtr, 'Sl_ram1', begg)
    end if
    if (NUOPC_IsConnected(exportState, fieldName='Sl_fv')) then
       call State_getFldPtr(exportState,'Sl_fv', dataPtr, rc=rc )
       do g = begg,endg
          i = 1 + (g-begg)
          dataPtr(i) = lnd2atm_inst%fv_grc(g)
       end do
       call check_for_nans(dataPtr, 'Sl_fv', begg)
    end if
    if (NUOPC_IsConnected(exportState, fieldName='Sl_soilw')) then
       call State_getFldPtr(exportState,'Sl_soilw', dataPtr, rc=rc )
       do g = begg,endg
          i = 1 + (g-begg)
          dataPtr(i) = lnd2atm_inst%h2osoi_vol_grc(g,1)
       end do
       call check_for_nans(dataPtr, 'Sl_soilw', begg)
    end if
    if (NUOPC_IsConnected(exportState, fieldName='Fall_flxdst1')) then
       call State_getFldPtr(exportState,'Fall_flxdst1', dataPtr, rc=rc )
       do g = begg,endg
          i = 1 + (g-begg)
          dataPtr(i) = -lnd2atm_inst%flxdst_grc(g,1)
       end do
       call check_for_nans(dataPtr, 'Fall_flxdst1', begg)
    end if
    if (NUOPC_IsConnected(exportState, fieldName='Fall_flxdst2')) then
       call State_getFldPtr(exportState,'Fall_flxdst2', dataPtr, rc=rc )
       do g = begg,endg
          i = 1 + (g-begg)
          dataPtr(i) = -lnd2atm_inst%flxdst_grc(g,1)
       end do
       call check_for_nans(dataPtr, 'Fall_flxdst2', begg)
    end if
    if (NUOPC_IsConnected(exportState, fieldName='Fall_flxdst3')) then
       call State_getFldPtr(exportState,'Fall_flxdst3', dataPtr, rc=rc )
       do g = begg,endg
          i = 1 + (g-begg)
          dataPtr(i) = -lnd2atm_inst%flxdst_grc(g,1)
       end do
       call check_for_nans(dataPtr, 'Fall_flxdst3', begg)
    end if
    if (NUOPC_IsConnected(exportState, fieldName='Fall_flxdst4')) then
       call State_getFldPtr(exportState,'Fall_flxdst4', dataPtr, rc=rc )
       do g = begg,endg
          i = 1 + (g-begg)
          dataPtr(i) = -lnd2atm_inst%flxdst_grc(g,1)
       end do
       call check_for_nans(dataPtr, 'Fall_flxdst4', begg)
    end if

    ! dry dep velocities
    do num = 1, shr_string_listGetNum(drydep_fields)
       call  shr_string_listGetName(drydep_fields, num, name)
       if (NUOPC_IsConnected(exportState, fieldName=trim(name))) then
          call State_getFldPtr(exportState, trim(name), dataPtr, rc=rc )
          do g = begg,endg
             dataPtr(1+g-begg) = lnd2atm_inst%ddvel_grc(g,num)
          end do
          call check_for_nans(dataPtr, trim(name), begg)
       end if
    end do

    ! MEGAN VOC emis fluxes
    do num = 1, megan_fields_num
       call  shr_string_listGetName(megan_voc_fields, num, name)
       if (NUOPC_IsConnected(exportState, fieldName=trim(name))) then
          call State_getFldPtr(exportState, trim(name), dataPtr, rc=rc )
          do g = begg,endg
             i = 1 + (g-begg)
             dataPtr(i) = -lnd2atm_inst%flxvoc_grc(g,num)
          end do
          call check_for_nans(dataPtr, trim(name), begg)
       end if
    end do

    ! fire emis fluxes
    do num = 1, fire_emis_fields_num
       call  shr_string_listGetName(fire_emis_fields, num, name)
       if (NUOPC_IsConnected(exportState, fieldName=trim(name))) then
          call State_getFldPtr(exportState, trim(name), dataPtr, rc=rc )
          do g = begg,endg
             i = 1 + (g-begg)
             dataPtr(i) = -lnd2atm_inst%fireflx_grc(g,num)
          end do
          call check_for_nans(dataPtr, trim(name), begg)
       end if
    end do
    if (NUOPC_IsConnected(exportState, fire_ztop_name)) then
       call State_getFldPtr(exportState, trim(fire_ztop_name), dataPtr, rc=rc )
       do g = begg,endg
          i = 1 + (g-begg)
          dataPtr(i) = lnd2atm_inst%fireztop_grc(g)
       end do
       call check_for_nans(dataPtr, trim(fire_ztop_name), begg)
    end if

    ! methane
    if (NUOPC_IsConnected(exportState, 'Fall_methane')) then
       call State_getFldPtr(exportState, 'Fall_methane', dataPtr, rc=rc )
       do g = begg,endg
          i = 1 + (g-begg)
          dataPtr(i) = -lnd2atm_inst%flux_ch4_grc(g)
       end do
       call check_for_nans(dataPtr, 'Fall_methane', begg)
    end if

    ! sign convention is positive downward with hierarchy of atm/glc/lnd/rof/ice/ocn.
    ! i.e. water sent from land to rof is positive

    ! surface runoff is the sum of qflx_over, qflx_h2osfc_surf
    call State_getFldPtr(exportState, 'Flrl_rofsur', dataPtr, rc=rc )
    do g = begg,endg
       i = 1 + (g-begg)
       dataPtr(i) = lnd2atm_inst%qflx_rofliq_qsur_grc(g) + lnd2atm_inst%qflx_rofliq_h2osfc_grc(g)
    end do
    call check_for_nans(dataPtr, 'Flrl_rofsur', begg)

    !  subsurface runoff is the sum of qflx_drain and qflx_perched_drain
    call State_getFldPtr(exportState, 'Flrl_rofsub', dataPtr, rc=rc )
    do g = begg,endg
       i = 1 + (g-begg)
       dataPtr(i) = lnd2atm_inst%qflx_rofliq_qsub_grc(g) + lnd2atm_inst%qflx_rofliq_drain_perched_grc(g)
    end do
    call check_for_nans(dataPtr, 'Flrl_rofsub', begg)

    !  qgwl sent individually to coupler
    call State_getFldPtr(exportState, 'Flrl_rofgwl', dataPtr, rc=rc )
    do g = begg,endg
       i = 1 + (g-begg)
       dataPtr(i) = lnd2atm_inst%qflx_rofliq_qgwl_grc(g)
    end do
    call check_for_nans(dataPtr, 'Flrl_rofgwl', begg)

    ! ice  sent individually to coupler
    call State_getFldPtr(exportState, 'Flrl_rofice', dataPtr, rc=rc )
    do g = begg,endg
       i = 1 + (g-begg)
       dataPtr(i) = lnd2atm_inst%qflx_rofice_grc(g)
    end do
    call check_for_nans(dataPtr, 'Flrl_rofice', begg)

    ! irrigation flux to be removed from main channel storage (negative)
    call State_getFldPtr(exportState, 'Flrl_irrig', dataPtr, rc=rc )
    do g = begg,endg
       i = 1 + (g-begg)
       dataPtr(i) = - lnd2atm_inst%qirrig_grc(g)
    end do
    call check_for_nans(dataPtr, 'Flrl_irrig', begg)

    ! glc fields
    ! We could avoid setting these fields if glc_present is .false., if that would
    ! help with performance. (The downside would be that we wouldn't have these fields
    ! available for diagnostic purposes or to force a later T compset with dlnd.)
    do num = 0,glc_nec
       nec_str = glc_elevclass_as_string(num)
       name = 'Sl_tsrf' // nec_str
       call State_getFldPtr(exportState,trim(name), dataPtr, rc=rc )
       do g = begg,endg
          i = 1 + (g-begg)
          dataPtr(i) = lnd2glc_inst%tsrf_grc(g,num)
       end do
       call check_for_nans(dataPtr, trim(name), begg)

       name = 'Sl_topo' // nec_str
       call State_getFldPtr(exportState,trim(name), dataPtr, rc=rc )
       do g = begg,endg
          i = 1 + (g-begg)
          dataPtr(i) = lnd2glc_inst%topo_grc(g,num)
       end do
       call check_for_nans(dataPtr, trim(name), begg)

       name = 'Flgl_qice' // nec_str
       call State_getFldPtr(exportState,trim(name), dataPtr, rc=rc )
       do g = begg,endg
          i = 1 + (g-begg)
          dataPtr = lnd2glc_inst%qice_grc(g,num)
       end do
       call check_for_nans(dataPtr, trim(name), begg)
    end do

  end subroutine export_fields

  !===============================================================================

  subroutine check_for_nans(array, fname, begg)
    real(r8), pointer             :: array(:)
    character(len=*) , intent(in) :: fname
    integer          , intent(in) :: begg

    ! Check if any input from mediator or output to mediator is NaN
    if (any(isnan(array))) then
       write(iulog,*) '# of NaNs = ', count(isnan(array))
       write(iulog,*) 'Which are NaNs = ', isnan(array)
       do i = 1, size(array))
          if (isnan(array(i))) then
             write(iulog,*) trim(fname)
             write(iulog,*) 'gridcell index = ', begg+i-1
          end if
       end do
       call shr_sys_abort(' ERROR: One or more of the output from CLM to the coupler are NaN ' )
    end if
  end subroutine check_for_nan

  !===============================================================================

  subroutine fld_list_add(num, fldlist, stdname)
    integer,                    intent(inout) :: num
    type(fld_list_type),        intent(inout) :: fldlist(:)
    character(len=*),           intent(in)    :: stdname

    ! local variables
    integer :: rc
    integer :: dbrc
    character(len=*), parameter :: subname='(dshr_nuopc_mod:fld_list_add)'
    !-------------------------------------------------------------------------------

    ! Set up a list of field information

    num = num + 1
    if (num > fldsMax) then
      call ESMF_LogWrite(trim(subname)//": ERROR num > fldsMax "//trim(stdname), &
        ESMF_LOGMSG_ERROR, line=__LINE__, file=__FILE__, rc=dbrc)
      return
    endif
    fldlist(num)%stdname = trim(stdname)

  end subroutine fld_list_add

  !===============================================================================

  subroutine fld_list_realize(state, fldList, numflds, flds_scalar_name, flds_scalar_num, mesh, tag, rc)

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
    character(len=*),parameter  :: subname='(dshr_nuopc_mod:fld_list_realize)'
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

  end subroutine fld_list_realize

end module lnd_import_export
