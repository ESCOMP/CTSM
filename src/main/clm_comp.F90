#include <misc.h>
#include <preproc.h>

module clm_comp

  use shr_kind_mod, only : r8 => shr_kind_r8
  use shr_sys_mod , only : shr_sys_abort
  use perf_mod

  implicit none

  SAVE
  private                              ! By default make data private

! !PUBLIC MEMBER FUNCTIONS:

  public clm_init0
  public clm_init1
  public clm_init2
  public clm_run1
  public clm_run2

contains

!---------------------------------------------------------------------------
!BOP
!
! !IROUTINE: clm_init0
!
! !INTERFACE:
  subroutine clm_init0( CCSMInit )
!
! !DESCRIPTION:
! Initialize land surface model and obtain relevant atmospheric model arrays
! back from (i.e. albedos, surface temperature and snow cover over land).
!
! !USES:
    use initializeMod,     only : initialize1
    use shr_InputInfo_mod, only : shr_InputInfo_initType
!
! !ARGUMENTS:
    type(shr_InputInfo_initType), intent(in), optional :: CCSMInit
!
! !LOCAL VARIABLES:
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
!-----------------------------------------------------------------------

    call t_startf('clm_init0')
    if ( present(CCSMInit) )then
       call initialize1( CCSMInit )
    else
       call initialize1( )
    end if
    call t_stopf('clm_init0')

  end subroutine clm_init0


!---------------------------------------------------------------------------
!BOP
!
! !IROUTINE: clm_init1
!
! !INTERFACE:
  subroutine clm_init1( SyncClock )
!
! !DESCRIPTION:
! Initialize land surface model and obtain relevant atmospheric model arrays
! back from (i.e. albedos, surface temperature and snow cover over land).
!
! !USES:
    use initializeMod,   only : initialize2
    use eshr_timemgr_mod,only : eshr_timemgr_clockType
!
! !ARGUMENTS:
   type(eshr_timeMgr_clockType), optional, intent(IN) :: SyncClock ! Synchronization clock
!
! !LOCAL VARIABLES:
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
!-----------------------------------------------------------------------

   call t_startf('clm_init1')
   if (present(SyncClock)) then	
      call initialize2( SyncClock )
   else
      call initialize2()
   end if
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
    use shr_orb_mod     , only : shr_orb_decl
    use clm_varctl      , only : finidat, nsrest
    use initSurfAlbMod  , only : initSurfAlb, do_initsurfalb 
    use clm_time_manager, only : get_nstep, get_step_size, get_curr_calday
    use clm_atmlnd      , only : clm_map2gcell
    use clm_varorb      , only : eccen, mvelpp, lambm0, obliqr
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

    call t_startf('clm_init2')
    if (get_nstep() == 0 .or. nsrest == 0) then

       ! Initialize albedos (correct pft filters are needed)

       if (finidat == ' ' .or. do_initsurfalb) then
          calday = get_curr_calday()
          call shr_orb_decl( calday, eccen, mvelpp, lambm0, obliqr, declin, eccf )
          
          dtime = get_step_size()
          caldaym1 = get_curr_calday(offset=-int(dtime))
          call shr_orb_decl( caldaym1, eccen, mvelpp, lambm0, obliqr, declinm1, eccf )
          
          call initSurfAlb( calday, declin, declinm1 )
       end if

       ! Determine gridcell averaged properties to send to atm

       call clm_map2gcell(init=.true.)

    end if
    call t_stopf('clm_init2')

  end subroutine clm_init2

!---------------------------------------------------------------------------
!BOP
!
! !IROUTINE: clm_run1
!
! !INTERFACE:
  subroutine clm_run1( nextsw_cday, dosend )
!
! !DESCRIPTION:
! land model run1 phase
!
! !USES:
    use shr_orb_mod     , only : shr_orb_decl
    use clm_varctl      , only : irad 
    use clm_time_manager, only : get_nstep, get_step_size, get_curr_calday
    use clm_varorb      , only : eccen, mvelpp, lambm0, obliqr
    use driver          , only : driver1
    use clm_atmlnd      , only : clm_map2gcell
!
! !ARGUMENTS:
    real(r8), intent(in) , optional :: nextsw_cday
    logical , intent(out), optional :: dosend
!
! !LOCAL VARIABLES:
    integer  :: nstep                 ! model time step
    integer  :: dtime                 ! time step increment (sec)
    logical  :: doalb                 ! true if surface albedo calculation time step from atm
    real(r8) :: caldayp1_clm          ! calendar day for nstep+1 for clm (debug need only)
    real(r8) :: caldayp1              ! calendar day for nstep+1
    real(r8) :: declinp1              ! solar declination angle in radians for nstep+1
    real(r8) :: eccf                  ! earth orbit eccentricity factor
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
!---------------------------------------------------------------------------

    call t_startf('clm_run1')

    ! Determine doalb (true when the next time step is a radiation time step) 

    nstep = get_nstep()
    doalb = ((irad==1) .or. (mod(nstep,irad)==0 .and. nstep/=0))
    dtime = get_step_size()

    if (present(nextsw_cday) .and. present(dosend)) then
       caldayp1 = nextsw_cday
       if (doalb) then
          caldayp1_clm = get_curr_calday( offset=dtime )
          if (caldayp1_clm /= caldayp1) then
             write(6,*)'caldayp1_atm= ',nextsw_cday,' caldayp1_atm= ',caldayp1_clm
             call shr_sys_abort('caldayp1 from atm and clm do not match')
          end if
       end if
       dosend = doalb
    else
       caldayp1 = get_curr_calday( offset=int(dtime) )
    endif

    ! Determine declination angle for next time step
    
    call shr_orb_decl( caldayp1, eccen, mvelpp, lambm0, obliqr, declinp1, eccf )

    ! Call land model driver1
    
    call t_startf('driver1')
    call driver1(doalb, caldayp1, declinp1)
    call t_stopf('driver1')

    ! Determine gridcell averaged properties to send to atm (l2as and l2af derived types)

    call t_startf('clm_map2gcell')
    call clm_map2gcell( )
    call t_stopf('clm_map2gcell')

    call t_stopf('clm_run1')

  end subroutine clm_run1

!---------------------------------------------------------------------------
!BOP
!
! !IROUTINE: clm_run2
!
! !INTERFACE:
  subroutine clm_run2(rstwr)
!
! !DESCRIPTION:
! land model run2 phase
!
! !USES:
    use shr_orb_mod     , only : shr_orb_decl
    use clm_time_manager, only : get_nstep, get_step_size, get_curr_calday
    use clm_varorb      , only : eccen, mvelpp, lambm0, obliqr
    use driver          , only : driver2
!
! !ARGUMENTS:
    logical,optional,intent(in) :: rstwr    ! true => write restart file this step
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
    real(r8) :: dtime                 ! time step increment (sec)
    real(r8) :: eccf                  ! earth orbit eccentricity factor
    real(r8) :: caldayp1              ! calendar day for nstep+1
    real(r8) :: declinp1              ! solar declination angle in radians for nstep+1
!---------------------------------------------------------------------------

    call t_startf('clm_run2')

    ! Call land model driver2
    
    dtime = get_step_size()
    caldayp1 = get_curr_calday( offset=int(dtime) )
    call shr_orb_decl( caldayp1, eccen, mvelpp, lambm0, obliqr, declinp1, eccf )
    call t_startf('driver2')
    if (present(rstwr)) then
       call driver2(caldayp1, declinp1, rstwr)
    else
       call driver2(caldayp1, declinp1)
    endif
    call t_stopf('driver2')

    call t_stopf('clm_run2')

  end subroutine clm_run2

end module clm_comp
