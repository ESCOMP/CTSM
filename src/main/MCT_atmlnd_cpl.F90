#include <misc.h>
#include <preproc.h>

module MCT_atmlnd_cpl

#if (defined COUP_CAM)

!---------------------------------------------------------------------
!
! Purpose:
!
! Collect coupling routines for sequential coupling of LND-ATM.
!       
!
! Author: R. Jacob, M. Vertenstein
!
!---------------------------------------------------------------------

  use shr_sys_mod
  use m_AttrVect  , only: AttrVect
  use m_Rearranger, only: Rearranger
  implicit none

  private  ! except

!--------------------------------------------------------------------------
! Public interfaces
!--------------------------------------------------------------------------

  public :: MCT_lnd2atm_init
  public :: MCT_atm2lnd_init
  public :: MCT_atm2lnd
  public :: MCT_lnd2atm

!--------------------------------------------------------------------------
! Public data
!--------------------------------------------------------------------------

  type(Rearranger), public :: Re_lnd2atm
  type(Rearranger), public :: Re_atm2lnd

#if !(defined SPMD)
  integer, parameter :: mpicom = 1
#endif

  character(*),parameter :: subName = '(atmlnd_MCTCplMod) '

!=======================================================================
   contains
!=======================================================================

  subroutine MCT_atm2lnd_init()

    use m_Rearranger, only: MCT_Rearr_init => init
    use MCT_atm_comp, only: GSMap_atm
    use MCT_lnd_comp, only: GSMap_lnd
#if (defined SPMD)
    use mpishorthand , only : mpicom
#endif

    implicit none

    call MCT_Rearr_init(GSMap_atm, GSMap_lnd, mpicom, Re_atm2lnd)

  end subroutine MCT_atm2lnd_init

  !=======================================================================

  subroutine MCT_lnd2atm_init()

    use m_Rearranger, only: MCT_Rearr_init => init
    use MCT_atm_comp, only: GSMap_atm
    use MCT_lnd_comp, only: GSMap_lnd
#if (defined SPMD)
    use mpishorthand , only : mpicom
#endif

    implicit none

    call MCT_Rearr_init(GSMap_lnd, GSMap_atm, mpicom, Re_lnd2atm)

  end subroutine MCT_lnd2atm_init

  !=======================================================================

  subroutine MCT_atm2lnd( a2c_a, a2c_l )

    ! Perform redistribution form atmosphere decomposition to land decomposition

    use m_Rearranger, only: Rearranger, MCT_Rearrange => rearrange

    implicit none
    type(AttrVect), intent(in)  :: a2c_a
    type(AttrVect), intent(out) :: a2c_l

    ! Perform redistribution form atmosphere decomposition to land decomposition

    call MCT_Rearrange(a2c_a, a2c_l, Re_atm2lnd)

  end subroutine MCT_atm2lnd

  !=======================================================================

  subroutine MCT_lnd2atm( l2c_l, l2c_a)

    ! Perform redistribution form land decomposition to atmosphere decomposition

    use m_Rearranger, only: Rearranger, MCT_Rearrange => rearrange

    implicit none
    type(AttrVect), intent(in)  :: l2c_l
    type(AttrVect), intent(out) :: l2c_a

    call MCT_Rearrange(l2c_l, l2c_a, Re_lnd2atm)

  end subroutine MCT_lnd2atm     

#endif

end module MCT_atmlnd_cpl
