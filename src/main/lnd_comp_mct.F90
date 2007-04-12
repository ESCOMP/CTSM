#include <misc.h>
#include <preproc.h>

module lnd_comp_mct
  
#if (defined SEQ_MCT)

!---------------------------------------------------------------------------
!BOP
!
! !MODULE: lnd_comp_mct
!
! !DESCRIPTION:
!
! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use shr_sys_mod , only : shr_sys_abort
  use mct_mod
  use seq_flds_mod
  use seq_flds_indices
  use seq_cdata_mod
  use seq_infobuf_mod
  use spmdMod
  use perf_mod
!
! !PUBLIC MEMBER FUNCTIONS:
  implicit none
  public :: lnd_init_mct
  public :: lnd_run_mct
  public :: lnd_final_mct
  SAVE
  private                              ! By default make data private
!
! ! PUBLIC DATA:
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
! !PRIVATE MEMBER FUNCTIONS:
  private :: lnd_SetgsMap_mct
  private :: lnd_domain_mct
  private :: lnd_export_mct
  private :: lnd_import_mct
!
! !PRIVATE VARIABLES
  integer, dimension(:), allocatable :: perm  ! permutation array to reorder points
!
! Time averaged flux fields
!  
  type(mct_aVect)   :: l2x_l_SNAP
  type(mct_aVect)   :: l2x_l_SUM
!
! Time averaged counter for flux fields
!
  integer :: avg_count

!===============================================================
contains
!===============================================================

!---------------------------------------------------------------------------
!BOP
!
! !IROUTINE: lnd_init_mct
!
! !INTERFACE:
  subroutine lnd_init_mct( cdata_l, x2l_l, l2x_l, &
                           cdata_r,        r2x_r, SyncClock, NLFilename )
!
! !DESCRIPTION:
! Initialize land surface model and obtain relevant atmospheric model arrays
! back from (i.e. albedos, surface temperature and snow cover over land).
!
! !USES:
    use clm_time_manager , only : get_nstep      
    use clm_atmlnd       , only : clm_mapl2a, clm_l2a, atm_l2a
    use domainMod        , only : adomain
    use clm_comp         , only : clm_init0, clm_init1, clm_init2
    use clm_varctl       , only : finidat,single_column
    use shr_inputinfo_mod, only : shr_inputInfo_initType,       &
                                  shr_inputInfo_initGetData
    use eshr_timemgr_mod , only : eshr_timemgr_clockType, eshr_timemgr_clockInfoType, &
                                  eshr_timemgr_clockGet
    use controlMod       , only : control_setNL
    use domainMod        , only : amask
!
! !ARGUMENTS:
    type(seq_cdata),              intent(inout) :: cdata_l
    type(mct_aVect),              intent(inout) :: x2l_l, l2x_l
    type(seq_cdata),              intent(inout) :: cdata_r
    type(mct_aVect),              intent(inout) ::        r2x_r
    type(eshr_timemgr_clockType), intent(in)    :: SyncClock
    character(len=*), optional,   intent(in)    :: NLFilename 
!
! !LOCAL VARIABLES:
    integer                                     :: LNDID	
    integer                                     :: mpicom_lnd       	
    type(mct_gsMap),              pointer       :: GSMap_lnd
    type(mct_gGrid),              pointer       :: dom_l
    type(shr_InputInfo_initType), pointer       :: CCSMInit
    integer  :: lsize                             ! size of attribute vector
    integer  :: i,j                               ! indices
    type(eshr_timemgr_clockInfoType) :: clockInfo ! Clock information including orbit
    character(len=32), parameter :: sub = 'lnd_init_mct'
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
!-----------------------------------------------------------------------

    ! set cdata data
    call seq_cdata_setptrs(cdata_l, ID=LNDID, mpicom=mpicom_lnd, &
         gsMap=GSMap_lnd, dom=dom_l, CCSMInit=CCSMInit)

    ! Initialize clm MPI communicator 

    call spmd_init( mpicom_lnd )
    
    ! Use CCSMInit to set orbital values

    call eshr_timemgr_clockGet( SyncClock, info=clockInfo )
    call lnd_setorb_mct( clockInfo )

    ! Consistency check on namelist filename	

    call control_setNL( 'lnd_in' )

    ! Initialize clm
    ! clm_init0 reads namelist, grid and surface data
    ! clm_init1 and clm_init2 performs rest of initialization	

    call clm_init0( CCSMInit )

    ! If in SCM mode and no land then exit out of initialization
    if ( single_column .and. amask(1)==0) return        

    call clm_init1( SyncClock )
    call clm_init2()

    ! Initialize MCT gsMap

    call lnd_SetgsMap_mct( mpicom_lnd, LNDID, gsMap_lnd ) 	
    lsize = mct_gsMap_lsize(gsMap_lnd, mpicom_lnd)

    ! Initialize MCT domain

    call lnd_domain_mct( lsize, gsMap_lnd, dom_l )

    ! Initialize MCT attribute vectors

    call mct_aVect_init(x2l_l, rList=seq_flds_x2l_fields, lsize=lsize)
    call mct_aVect_zero(x2l_l)

    call mct_aVect_init(l2x_l, rList=seq_flds_l2x_fields, lsize=lsize)
    call mct_aVect_zero(l2x_l)

    call mct_aVect_init(l2x_l_SNAP, rList=seq_flds_l2x_fluxes, lsize=lsize)
    call mct_aVect_zero(l2x_l_SNAP)

    call mct_aVect_init(l2x_l_SUM , rList=seq_flds_l2x_fluxes, lsize=lsize)
    call mct_aVect_zero(l2x_l_SUM )

    if (masterproc) then
       write(6,*)'LND_INIT_MCT: time averaging the following flux fields over the coupling interval'
       write(6,*)trim(seq_flds_l2x_fluxes)
    end if

    ! Map internal data structure into coupling data structure

    call clm_mapl2a(clm_l2a, atm_l2a)
    call lnd_export_mct( atm_l2a, l2x_l )

    ! Initialize averaging counter

    avg_count = 0

  end subroutine lnd_init_mct

!---------------------------------------------------------------------------
!BOP
!
! !IROUTINE: lnd_run_mct
!
! !INTERFACE:
  subroutine lnd_run_mct( cdata_l, x2l_l, l2x_l, cdata_r, r2x_r, SyncClock )
!
! !DESCRIPTION:
! Run clm model
!
! !USES:
    use clm_atmlnd      ,only : clm_mapl2a, clm_mapa2l
    use clm_atmlnd      ,only : clm_l2a, atm_l2a, atm_a2l, clm_a2l
    use clm_comp        ,only : clm_run1, clm_run2
    use eshr_timemgr_mod,only : eshr_timemgr_clockType,         &
                                eshr_timemgr_clockAlarmIsOnRes, &
                                eshr_timemgr_clockDateInSync
    use clm_time_manager,only : get_curr_date, advance_timestep
!
! !ARGUMENTS:
    type(seq_cdata)             , intent(in)    :: cdata_l
    type(mct_aVect)             , intent(inout) :: x2l_l
    type(mct_aVect)             , intent(inout) :: l2x_l
    type(seq_cdata)             , intent(in)    :: cdata_r
    type(mct_aVect)             , intent(inout) :: r2x_r
    type(eshr_timemgr_clockType), intent(in)    :: SyncClock
!
! !LOCAL VARIABLES:
    character(len=32), parameter :: sub = "lnd_run_mct"
    logical :: rstwr   ! If time to write restart file
    logical :: dosend  ! true => send data back to driver
    integer :: ymd     ! Current date (YYYYMMDD)
    integer :: yr      ! Current year
    integer :: mon     ! Current month
    integer :: day     ! Current day
    integer :: tod     ! Current time of day (sec)
    integer :: ier     ! MPI error return
    type(seq_infobuf), pointer :: infobuf_l  
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
!---------------------------------------------------------------------------

    ! Get cdata pointers

    call seq_cdata_setptrs(cdata_l, infobuf=infobuf_l)

    ! Check that internal clock in sync with master clock

    call get_curr_date(yr, mon, day, tod )
    ymd = yr*10000 + mon*100 + day
    if ( .not. eshr_timemgr_clockDateInSync( SyncClock, ymd, tod ) )then
       call shr_sys_abort( sub//":: Internal CLM clock not in sync with "// &
                    "Master Synchronization clock" )
    end if
    rstwr = eshr_timemgr_clockAlarmIsOnRes( SyncClock )

    ! Map MCT to land data type

    call t_startf ('lc_lnd_import')
    call lnd_import_mct( x2l_l, atm_a2l )
    call t_stopf ('lc_lnd_import')

    call t_startf ('lc_clm_mapa2l')
    call clm_mapa2l(atm_a2l, clm_a2l)
    call t_stopf ('lc_clm_mapa2l')
    
    ! Loop over time steps in coupling interval

    dosend = .false.
    do while(.not. dosend)

       ! Run clm 

       call t_barrierf('sync_clm_run1', mpicom)
       call t_startf ('lc_clm_run1')
       call clm_run1(infobuf_l%rbuf(rbuf_nextsw_cday), dosend)
       call t_stopf ('lc_clm_run1')

       call t_barrierf('sync_clm_run2', mpicom)
       call t_startf ('lc_clm_run2')
       call clm_run2( rstwr )
       call t_stopf ('lc_clm_run2')

       ! Map land data type to MCT
       
       call t_startf ('lc_clm_mapl2a')
       call clm_mapl2a(clm_l2a, atm_l2a)
       call t_stopf ('lc_clm_mapl2a')
       
       call t_startf ('lc_lnd_export')
       call lnd_export_mct( atm_l2a, l2x_l )
       call t_stopf ('lc_lnd_export')
       
       ! Compute snapshot attribute vector for accumulation
       
       call mct_aVect_copy( l2x_l, l2x_l_SNAP )
       call mct_aVect_accum( aVin=l2x_l_SNAP, aVout=l2x_l_SUM )
       avg_count = avg_count + 1
       
       ! Advance clm time step
       
       call t_startf ('lc_clm2_adv_timestep')
       call advance_timestep()
       call t_stopf ('lc_clm2_adv_timestep')
	
       ! ***For now - force data to be sent back every time step***
       dosend = .true.

    end do

    ! Finish accumulation of attribute vector and average and zero out partial sum and counter
    
    call mct_aVect_avg ( l2x_l_SUM, avg_count)
    call mct_aVect_copy( l2x_l_SUM, l2x_l )
    call mct_aVect_zero( l2x_l_SUM) 
    avg_count = 0                   

  end subroutine lnd_run_mct

!---------------------------------------------------------------------------
!BOP
!
! !IROUTINE: lnd_final_mct
!
! !INTERFACE:
  subroutine lnd_final_mct( )
!
! !DESCRIPTION:
! Finalize land surface model
!
!------------------------------------------------------------------------------
!BOP
!
! !ARGUMENTS:
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
!---------------------------------------------------------------------------

   ! fill this in
  end subroutine lnd_final_mct

!=================================================================================

  subroutine lnd_SetgsMap_mct( mpicom_lnd, LNDID, gsMap_lnd )

    !-------------------------------------------------------------------
    !
    ! Uses
    !
    use decompMod, only : get_proc_bounds_atm, adecomp
    use domainMod, only : adomain
    !
    ! Arguments
    !
    integer        , intent(in)  :: mpicom_lnd
    integer        , intent(in)  :: LNDID
    type(mct_gsMap), intent(out) :: gsMap_lnd
    !
    ! Local Variables
    !
    integer,allocatable :: gindex(:)
    integer :: i, j, n, gi
    integer :: lsize,gsize
    integer :: ier
    integer :: begg, endg    
    !-------------------------------------------------------------------

    ! Build the land grid numbering for MCT
    ! NOTE:  Numbering scheme is: West to East and South to North
    ! starting at south pole.  Should be the same as what's used in SCRIP
    
    call get_proc_bounds_atm(begg, endg)

    allocate(gindex(begg:endg),stat=ier)

    ! number the local grid

    do n = begg, endg
        gindex(n) = adecomp%gdc2glo(n)
    end do
    lsize = endg-begg+1
    gsize = adomain%ni*adomain%nj

    ! reorder gindex to be in ascending order, initialize a permutation array,
    ! derive a permutation that puts gindex in ascending order since the
    ! the default for IndexSort is ascending and finally sort gindex in-place

    allocate(perm(lsize),stat=ier)
    call mct_indexset(perm)
    call mct_indexsort(lsize,perm,gindex)
    call mct_permute(gindex,perm,lsize)
    call mct_gsMap_init( gsMap_lnd, gindex, mpicom_lnd, LNDID, lsize, gsize )

    deallocate(gindex)

  end subroutine lnd_SetgsMap_mct

!====================================================================================

  subroutine lnd_export_mct( l2a, l2x_l )   

    !-----------------------------------------------------
    use clm_time_manager, only : get_nstep  
    use clm_atmlnd  , only : lnd2atm_type
    use domainMod   , only : adomain
    use decompMod   , only : get_proc_bounds_atm, adecomp

    type(lnd2atm_type), intent(inout) :: l2a
    type(mct_aVect)   , intent(inout) :: l2x_l

    integer :: g,i
    integer :: begg, endg    ! beginning and ending gridcell indices
    !-----------------------------------------------------
    
    call get_proc_bounds_atm(begg, endg)

    l2x_l%rAttr(:,:) = 0.0_r8

    ! ccsm sign convention is that fluxes are positive downward

!dir$ concurrent
    do g = begg,endg
       i = 1 + (g-begg)
       l2x_l%rAttr(index_l2x_Sl_landfrac,i) =  adomain%frac(g)
       l2x_l%rAttr(index_l2x_Sl_t,i)        =  l2a%t_rad(g)
       l2x_l%rAttr(index_l2x_Sl_snowh,i)    =  l2a%h2osno(g)
       l2x_l%rAttr(index_l2x_Sl_avsdr,i)    =  l2a%albd(g,1)
       l2x_l%rAttr(index_l2x_Sl_anidr,i)    =  l2a%albd(g,2)
       l2x_l%rAttr(index_l2x_Sl_avsdf,i)    =  l2a%albi(g,1)
       l2x_l%rAttr(index_l2x_Sl_anidf,i)    =  l2a%albi(g,2)
       l2x_l%rAttr(index_l2x_Sl_tref,i)     =  l2a%t_ref2m(g)
       l2x_l%rAttr(index_l2x_Sl_qref,i)     =  l2a%q_ref2m(g)
       l2x_l%rAttr(index_l2x_Fall_taux,i)   = -l2a%taux(g)
       l2x_l%rAttr(index_l2x_Fall_tauy,i)   = -l2a%tauy(g)
       l2x_l%rAttr(index_l2x_Fall_lat,i)    = -l2a%eflx_lh_tot(g)
       l2x_l%rAttr(index_l2x_Fall_sen,i)    = -l2a%eflx_sh_tot(g)
       l2x_l%rAttr(index_l2x_Fall_lwup,i)   = -l2a%eflx_lwrad_out(g)
       l2x_l%rAttr(index_l2x_Fall_evap,i)   = -l2a%qflx_evap_tot(g)
       l2x_l%rAttr(index_l2x_Fall_swnet,i)  = -l2a%fsa(g)
       if (index_l2x_Fall_nee /= 0) then
          l2x_l%rAttr(index_l2x_Fall_nee,i) = -l2a%nee(g)  
       end if

       ! optional fields for dust.  The index = 0 is a good way to flag it,
       ! but I have set it up so that l2a doesn't have ram1,fv,flxdst[1-4] if
       ! progsslt or dust aren't running. 
#if ( defined DUST || defined PROGSSLT )
       if (index_l2x_Sl_ram1 /= 0 )  l2x_l%rAttr(index_l2x_Sl_ram1,i) = l2a%ram1(g)
       if (index_l2x_Sl_fv   /= 0 )  l2x_l%rAttr(index_l2x_Sl_fv,i)   = l2a%fv(g)
#endif
#if ( defined DUST )
       if (index_l2x_Fall_flxdst1 /= 0 )  l2x_l%rAttr(index_l2x_Fall_flxdst1,i)= -l2a%flxdst(g,1)
       if (index_l2x_Fall_flxdst2 /= 0 )  l2x_l%rAttr(index_l2x_Fall_flxdst2,i)= -l2a%flxdst(g,2)
       if (index_l2x_Fall_flxdst3 /= 0 )  l2x_l%rAttr(index_l2x_Fall_flxdst3,i)= -l2a%flxdst(g,3)
       if (index_l2x_Fall_flxdst4 /= 0 )  l2x_l%rAttr(index_l2x_Fall_flxdst4,i)= -l2a%flxdst(g,4)
#endif
    end do

    ! permute before using the Rearrange call.

    call mct_aVect_permute(l2x_l,perm)

  end subroutine lnd_export_mct

!====================================================================================

  subroutine lnd_import_mct( x2l_l, a2l )

    !-----------------------------------------------------
    use clm_atmlnd      , only: atm2lnd_type
    use clm_varctl      , only: co2_type
    use clm_varcon      , only: rair, o2_molar_const, co2_ppmv_const, c13ratio
    use decompMod       , only: get_proc_bounds_atm
    !
    ! Arguments
    !
    type(mct_aVect)   , intent(inout) :: x2l_l
    type(atm2lnd_type), intent(inout) :: a2l
    !
    ! Local Variables
    !
    integer  :: g,i,nstep,ier
    real(r8) :: forc_rainc    ! rainxy Atm flux mm/s
    real(r8) :: forc_rainl    ! rainxy Atm flux mm/s
    real(r8) :: forc_snowc    ! snowfxy Atm flux  mm/s
    real(r8) :: forc_snowl    ! snowfxl Atm flux  mm/s
    real(r8) :: co2_ppmv_diag ! temporary
    real(r8) :: co2_ppmv_prog ! temporary
    real(r8) :: co2_ppmv      ! temporary
    integer  :: begg, endg    ! beginning and ending gridcell indices
    integer  :: co2_type_idx  ! integer flag for co2_type options
    !-----------------------------------------------------

    ! unpermute after rearrange call and before copying into local arrays.

    call mct_aVect_unpermute(x2l_l, perm)

    call get_proc_bounds_atm(begg, endg)

    co2_type_idx = 0
    if (co2_type == 'prognostic') then
       co2_type_idx = 1
    else if (co2_type == 'diagnostic') then
       co2_type_idx = 2
    end if
    if (co2_type == 'prognostic' .and. index_x2l_Sa_co2prog == 0) then
       write(6,*)' must have nonzero index_x2l_Sa_co2prog for co2_type equal to prognostic'
       call shr_sys_abort()
    else if (co2_type == 'diagnostic' .and. index_x2l_Sa_co2diag == 0) then
       write(6,*)' must have nonzero index_x2l_Sa_co2diag for co2_type equal to diagnostic'
       call shr_sys_abort()
    end if
    
    ! Note that the precipitation fluxes received  from the coupler
    ! are in units of kg/s/m^2. To convert these precipitation rates
    ! in units of mm/sec, one must divide by 1000 kg/m^3 and multiply
    ! by 1000 mm/m resulting in an overall factor of unity.
    ! Below the units are therefore given in mm/s.
    
!dir$ concurrent
    do g = begg,endg
        i = 1 + (g - begg)
       
        ! Determine required receive fields

        a2l%forc_hgt(g)     = x2l_l%rAttr(index_x2l_Sa_z,i)         ! zgcmxy  Atm state m
        a2l%forc_u(g)       = x2l_l%rAttr(index_x2l_Sa_u,i)         ! forc_uxy  Atm state m/s
        a2l%forc_v(g)       = x2l_l%rAttr(index_x2l_Sa_v,i)         ! forc_vxy  Atm state m/s
        a2l%forc_th(g)      = x2l_l%rAttr(index_x2l_Sa_ptem,i)      ! forc_thxy Atm state K
        a2l%forc_q(g)       = x2l_l%rAttr(index_x2l_Sa_shum,i)      ! forc_qxy  Atm state kg/kg
        a2l%forc_pbot(g)    = x2l_l%rAttr(index_x2l_Sa_pbot,i)      ! ptcmxy  Atm state Pa
        a2l%forc_t(g)       = x2l_l%rAttr(index_x2l_Sa_tbot,i)      ! forc_txy  Atm state K
        a2l%forc_lwrad(g)   = x2l_l%rAttr(index_x2l_Faxa_lwdn,i)    ! flwdsxy Atm flux  W/m^2
        forc_rainc          = x2l_l%rAttr(index_x2l_Faxa_rainc,i)   ! mm/s
        forc_rainl          = x2l_l%rAttr(index_x2l_Faxa_rainl,i)   ! mm/s
        forc_snowc          = x2l_l%rAttr(index_x2l_Faxa_snowc,i)   ! mm/s
        forc_snowl          = x2l_l%rAttr(index_x2l_Faxa_snowl,i)   ! mm/s
        a2l%forc_solad(g,2) = x2l_l%rAttr(index_x2l_Faxa_swndr,i)   ! forc_sollxy  Atm flux  W/m^2
        a2l%forc_solad(g,1) = x2l_l%rAttr(index_x2l_Faxa_swvdr,i)   ! forc_solsxy  Atm flux  W/m^2
        a2l%forc_solai(g,2) = x2l_l%rAttr(index_x2l_Faxa_swndf,i)   ! forc_solldxy Atm flux  W/m^2
        a2l%forc_solai(g,1) = x2l_l%rAttr(index_x2l_Faxa_swvdf,i)   ! forc_solsdxy Atm flux  W/m^2

        ! Determine optional receive fields

        if (index_x2l_Sa_co2prog /= 0) then
           co2_ppmv_prog = x2l_l%rAttr(index_x2l_Sa_co2prog,i)   ! co2 atm state prognostic
        else
           co2_ppmv_prog = co2_ppmv_const
        end if
 
        if (index_x2l_Sa_co2diag /= 0) then
           co2_ppmv_diag = x2l_l%rAttr(index_x2l_Sa_co2diag,i)   ! co2 atm state diagnostic
        else
           co2_ppmv_diag = co2_ppmv_const
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
        
        ! Determine derived quantities for optional fields
        ! Note that the following does unit conversions from ppmv to partial pressures (Pa)
        ! Note that forc_pbot is in Pa

        if (co2_type_idx == 1) then
           co2_ppmv = co2_ppmv_prog
        else if (co2_type_idx == 2) then
           co2_ppmv = co2_ppmv_diag 
        else
           co2_ppmv = co2_ppmv_const      
        end if
        a2l%forc_pco2(g)   = co2_ppmv * 1.e-6_r8 * a2l%forc_pbot(g) 
        a2l%forc_pc13o2(g) = co2_ppmv * c13ratio * 1.e-6_r8 * a2l%forc_pbot(g)
	 
     end do

   end subroutine lnd_import_mct

!---------------------------------------------------------------------------
!BOP
!
! !IROUTINE: lnd_setorb_mct
!
! !INTERFACE:
  subroutine lnd_setorb_mct( ClockInfo )
!
! !DESCRIPTION:
! Determine clm orbital parameters
!
! !USES:
    use clm_varorb      , only : iyear_AD, eccen, obliq, nmvelp, obliqr, &
                                 lambm0, mvelpp
    use eshr_timemgr_mod, only : eshr_timemgr_clockInfoType, &
                                 eshr_timemgr_clockInfoGet
    use shr_orb_mod     , only : shr_orb_params
!
! !ARGUMENTS: 
    type(eshr_timemgr_clockInfoType), intent(in) :: ClockInfo

!
!EOP
!-----------------------------------------------------------------------

    ! Get orbital information from clockInfo object
    call eshr_timemgr_clockInfoGet( ClockInfo, orb_eccen=eccen, orb_mvelp=nmvelp, &
                                    orb_iyear_AD=iyear_AD, orb_obliq=obliq )
    ! Set orbital parameters based on the above values
    call shr_orb_params( iyear_AD, eccen, obliq, nmvelp,     &
                         obliqr, lambm0, mvelpp,   Log_Print=.false. )

  end subroutine lnd_setorb_mct

!===============================================================================

  subroutine lnd_domain_mct( lsize, gsMap_l, dom_l )

    !-------------------------------------------------------------------
    use clm_varcon, only : re
    use domainMod , only : adomain
    use decompMod , only : get_proc_bounds_atm, adecomp
    !
    ! Arguments
    !
    integer        , intent(in)    :: lsize
    type(mct_gsMap), intent(inout) :: gsMap_l
    type(mct_ggrid), intent(out)   :: dom_l      
    !
    ! Local Variables
    !
    integer :: g,i,j              ! index
    integer :: begg, endg         ! beginning and ending gridcell indices
    real(r8), pointer :: data(:)  ! temporary
    integer , pointer :: idata(:) ! temporary
    !-------------------------------------------------------------------
    !
    ! Initialize mct domain type
    !
    call mct_gGrid_init( GGrid=dom_l, CoordChars="lat:lon", OtherChars="area:mask:aream", lsize=lsize )
    !
    ! Allocate memory
    !
    allocate(data(lsize))
    !
    ! Determine global gridpoint number attribute, GlobGridNum, which is set automatically by MCT
    !
    call mct_gsMap_orderedPoints(gsMap_l, iam, idata)
    call mct_gGrid_importIAttr(dom_l,'GlobGridNum',idata,lsize)
    !
    ! Determine domain (numbering scheme is: West to East and South to North to South pole)
    ! Initialize attribute vector with special value
    !
    data(:) = -9999.0_R8 
    call mct_gGrid_importRAttr(dom_l,"lat"  ,data,lsize) 
    call mct_gGrid_importRAttr(dom_l,"lon"  ,data,lsize) 
    call mct_gGrid_importRAttr(dom_l,"area" ,data,lsize) 
    call mct_gGrid_importRAttr(dom_l,"aream",data,lsize) 
    data(:) = 0.0_R8     
    call mct_gGrid_importRAttr(dom_l,"mask" ,data,lsize) 
    !
    ! Determine bounds
    !
    call get_proc_bounds_atm(begg, endg)
    !
    ! Fill in correct values for domain components
    ! Note aream will be filled in in the atm-lnd mapper
    !
    do g = begg,endg
       i = 1 + (g - begg)
       data(i) = adomain%lonc(g)
    end do
    call mct_gGrid_importRattr(dom_l,"lon",data,lsize) 

    do g = begg,endg
       i = 1 + (g - begg)
       data(i) = adomain%latc(g)
    end do
    call mct_gGrid_importRattr(dom_l,"lat",data,lsize) 

    do g = begg,endg
       i = 1 + (g - begg)
       data(i) = adomain%area(g)/(re*re)
    end do
    call mct_gGrid_importRattr(dom_l,"area",data,lsize) 

    do g = begg,endg
       i = 1 + (g - begg)
       data(i) = real(adomain%mask(g), r8)
    end do
    call mct_gGrid_importRattr(dom_l,"mask",data,lsize) 

    ! Permute dom_l to have ascending order

    call mct_gGrid_permute(dom_l, perm)

    deallocate(data)
    deallocate(idata)

  end subroutine lnd_domain_mct
    
#endif

end module lnd_comp_mct
