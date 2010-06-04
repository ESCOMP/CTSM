module clm_comp
!---------------------------------------------------------------------------
!BOP
!
! !MODULE: clm_comp
!
! !DESCRIPTION: 
!
! High level initial and run phase component module for CLM. This module interfaces
! with the generic lnd_comp_<framework> module to send/receive land model information
! with the interface desired by CCSM for initialization phase and run time stepping
! phases.
!
!---------------------------------------------------------------------------

! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use shr_sys_mod , only : shr_sys_abort
  use perf_mod

  implicit none

  SAVE
  private                              ! By default make data private

! !PUBLIC MEMBER FUNCTIONS:

  public clm_init1   ! Phase 1 initialization, namelist, setup grid etc.
  public clm_init2   ! Phase 2 initialization, initialize parameterizations
  public clm_init3   ! Phase 3 initialization, initialize albedos
  public clm_run1    ! Phase 1 run
  public clm_run2    ! Phase 2 run
!EOP
!---------------------------------------------------------------------------

contains

!---------------------------------------------------------------------------
!BOP
!
! !IROUTINE: clm_init1
!
! !INTERFACE:
  subroutine clm_init1( )
!
! !DESCRIPTION:
! Initialize land surface model and obtain relevant atmospheric model arrays
! back from (i.e. albedos, surface temperature and snow cover over land).
!
! !USES:
    use clm_initializeMod,     only : initialize1
!
! !ARGUMENTS:
!
! !LOCAL VARIABLES:
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
!-----------------------------------------------------------------------

    call t_startf('clm_init1')
    call initialize1( )
    call t_stopf('clm_init1')

  end subroutine clm_init1

!---------------------------------------------------------------------------
!BOP
!
! !IROUTINE: clm_init2
!
! !INTERFACE:
  subroutine clm_init2( )
!
! !DESCRIPTION:
! Initialize land surface model and obtain relevant atmospheric model arrays
! back from (i.e. albedos, surface temperature and snow cover over land).
!
! !USES:
    use clm_initializeMod, only : initialize2
    use mct_mod
!
! !ARGUMENTS    
    implicit none
!
! !LOCAL VARIABLES:
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
!-----------------------------------------------------------------------

   call t_startf('clm_init2')
   call initialize2( )
   call t_stopf('clm_init2')

  end subroutine clm_init2

!---------------------------------------------------------------------------
!BOP
!
! !IROUTINE: clm_init3
!
! !INTERFACE:
  subroutine clm_init3( )
!
! !DESCRIPTION:
! Initialize land surface model and obtain relevant atmospheric model arrays
! back from (i.e. albedos, surface temperature and snow cover over land).
!
! !USES:
    use shr_orb_mod        , only : shr_orb_decl
    use clm_varctl         , only : finidat, nsrest
    use initSurfAlbMod     , only : initSurfAlb, do_initsurfalb 
    use clm_time_manager   , only : get_nstep, get_step_size, get_curr_calday
    use clm_atmlnd         , only : clm_map2gcell
    use clm_varorb         , only : eccen, mvelpp, lambm0, obliqr
    use seq_drydep_mod     , only : n_drydep, drydep_method, DD_XLND
    use STATICEcosysDynMod , only : interpMonthlyVeg
!
! !ARGUMENTS:
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
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
!-----------------------------------------------------------------------

    call t_startf('clm_init3')
    if (get_nstep() == 0 .or. nsrest == 0) then

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
    call t_stopf('clm_init3')

  end subroutine clm_init3

!---------------------------------------------------------------------------
!BOP
!
! !IROUTINE: clm_run1
!
! !INTERFACE:
  subroutine clm_run1( doalb, nextsw_cday )
!
! !DESCRIPTION:
! land model run1 phase
!
! !USES:
    use shr_orb_mod     , only : shr_orb_decl
    use clm_time_manager, only : get_nstep, get_step_size, get_curr_calday
    use clm_varorb      , only : eccen, mvelpp, lambm0, obliqr
    use clm_driver      , only : clm_driver1
    use clm_atmlnd      , only : clm_map2gcell
!
! !ARGUMENTS:
    logical , intent(IN) :: doalb       ! true if surface albedo calculation time step from atm
    real(r8), intent(IN) :: nextsw_cday ! calendar day for nstep+1
!
! !LOCAL VARIABLES:
    integer  :: dtime                 ! time step increment (sec)
    real(r8) :: calday                ! calendar day for nstep
    real(r8) :: declin                ! solar declination angle in radians for nstep
    real(r8) :: declinp1              ! solar declination angle in radians for nstep+1
    real(r8) :: eccf                  ! earth orbit eccentricity factor
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
!---------------------------------------------------------------------------

    ! Set default values first 

    dtime = get_step_size()
    calday   = get_curr_calday()

    ! Determine declination angle for current and next time step
    
    call shr_orb_decl( calday     , eccen, mvelpp, lambm0, obliqr, declin  , eccf )
    call shr_orb_decl( nextsw_cday, eccen, mvelpp, lambm0, obliqr, declinp1, eccf )

    ! Call land model driver1
    
    call t_startf('clm_driver1')
    call clm_driver1(doalb, nextsw_cday, declinp1, declin)
    call t_stopf('clm_driver1')

    ! Determine gridcell averaged properties to send to atm (l2as and l2af derived types)

    call t_startf('clm_map2gcell')
    call clm_map2gcell( )
    call t_stopf('clm_map2gcell')

  end subroutine clm_run1

!---------------------------------------------------------------------------
!BOP
!
! !IROUTINE: clm_run2
!
! !INTERFACE:
  subroutine clm_run2( nextsw_cday, rstwr, nlend, rdate )
!
! !DESCRIPTION:
! land model run2 phase
!
! !USES:
    use shr_orb_mod     , only : shr_orb_decl
    use clm_time_manager, only : get_nstep, get_step_size, get_curr_calday
    use clm_varorb      , only : eccen, mvelpp, lambm0, obliqr
    use clm_driver      , only : clm_driver2
!
! !ARGUMENTS:
    real(r8)                 ,intent(in) :: nextsw_cday ! calendar day for nextsw_cday
    logical                  ,intent(in) :: rstwr    ! true => write restart file this step
    logical                  ,intent(in) :: nlend    ! true => end run this step
    character(len=*)         ,intent(in) :: rdate    ! time stamp for restart file names
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!
! !LOCAL VARIABLES:
!EOP
    real(r8) :: eccf                  ! earth orbit eccentricity factor
    real(r8) :: declinp1              ! solar declination angle in radians for nstep+1
!---------------------------------------------------------------------------

    ! Call land model driver2
    
    call shr_orb_decl( nextsw_cday, eccen, mvelpp, lambm0, obliqr, declinp1, eccf )
    call clm_driver2(nextsw_cday, declinp1, rstwr, nlend, rdate)

  end subroutine clm_run2

end module clm_comp
