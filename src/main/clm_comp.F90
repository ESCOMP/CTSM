#include <misc.h>
#include <preproc.h>

module clm_comp

  use shr_kind_mod, only : r8 => shr_kind_r8
  implicit none

  SAVE
  private                              ! By default make data private

! !PUBLIC MEMBER FUNCTIONS:

  public clm_init1
  public clm_init2
  public clm_run1
  public clm_run2

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
    use initializeMod, only : initialize           
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

    call initialize()

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
    use time_manager    , only : get_nstep, get_step_size, get_curr_calday
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

  end subroutine clm_init2

!---------------------------------------------------------------------------
!BOP
!
! !IROUTINE: clm_run1
!
! !INTERFACE:
  subroutine clm_run1( )
!
! !DESCRIPTION:
! land model run1 phase
!
! !USES:
    use shr_orb_mod , only : shr_orb_decl
    use clm_varctl  , only : irad 
    use time_manager, only : get_nstep, get_step_size, get_curr_calday
    use clm_varorb  , only : eccen, mvelpp, lambm0, obliqr
    use driver      , only : driver1
    use clm_atmlnd  , only : clm_map2gcell
!
! !ARGUMENTS:
!
! !LOCAL VARIABLES:
    logical  :: doalb                 ! true if surface albedo calculation time step
    real(r8) :: dtime                 ! time step increment (sec)
    integer  :: nstep                 ! model time step
    real(r8) :: caldayp1              ! calendar day for nstep+1
    real(r8) :: declinp1              ! solar declination angle in radians for nstep+1
    real(r8) :: eccf                  ! earth orbit eccentricity factor
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
!---------------------------------------------------------------------------

    ! Determine doalb (true when the next time step is a radiation time step) 

    nstep = get_nstep()
    doalb = ((irad==1) .or. (mod(nstep,irad)==0 .and. nstep/=0))

    ! Determine declination angle for next time step
    
    dtime = get_step_size()
    caldayp1 = get_curr_calday( offset=int(dtime) )
    call shr_orb_decl( caldayp1, eccen, mvelpp, lambm0, obliqr, declinp1, eccf )

    ! Call land model driver1
    
    call t_startf('clm_driver1')
    call driver1(doalb, caldayp1, declinp1)
    call t_stopf('clm_driver1')

    ! Determine gridcell averaged properties to send to atm (l2as and l2af derived types)

    call clm_map2gcell( )

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
    use shr_orb_mod , only : shr_orb_decl
    use time_manager, only : get_nstep, get_step_size, get_curr_calday
    use clm_varorb  , only : eccen, mvelpp, lambm0, obliqr
    use driver      , only : driver2
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

    ! Call land model driver2
    
    call t_startf('clm_driver2')

    dtime = get_step_size()
    caldayp1 = get_curr_calday( offset=int(dtime) )
    call shr_orb_decl( caldayp1, eccen, mvelpp, lambm0, obliqr, declinp1, eccf )
    if (present(rstwr)) then
       call driver2(caldayp1, declinp1, rstwr)
    else
       call driver2(caldayp1, declinp1)
    endif

    call t_stopf('clm_driver2')

  end subroutine clm_run2

end module clm_comp
