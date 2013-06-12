module initParametersMod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: initParametersMod
!
! !DESCRIPTION:
! Do some parameter initialization that doesn't fit well elsewhere
!
! !USES:
!
! !PUBLIC TYPES:
  implicit none
  private
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public initParameters      ! initialize parameters
!
! !REVISION HISTORY
!   Created by Bill Sacks
!
!EOP
!-----------------------------------------------------------------------

contains

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: clm_varpar_init
!
! !INTERFACE:
  subroutine initParameters()
!
! !DESCRIPTION:
! This subroutine does some parameter initialization that doesn't fit well elsewhere
!
! !USES:
    use clm_varpar, only : numpft, numcft, natpft_size, cft_size, &
                           natpft_lb, natpft_ub, cft_lb, cft_ub
    use clm_varctl, only : create_crop_landunit
!
! !ARGUMENTS:
    implicit none
!
! !REVISION HISTORY:
!   Created by Bill Sacks
! !LOCAL VARIABLES:
!
!EOP
!------------------------------------------------------------------------------
    
    ! The following initialization really belongs in clm_varpar_init, but it can't go
    ! there because it needs create_crop_landunit - and clm_varpar can't use clm_varctl's
    ! create_crop_landunit because that would create a circular dependency (because
    ! clm_varctl also uses clm_varpar)

    ! For arrays containing all PFTs (natural veg & crop), determine lower and upper bounds
    ! for (1) PFTs on the natural vegetation landunit (includes bare ground, and includes
    ! crops if create_crop_landunit=false), and (2) CFTs on the crop landunit (no elements
    ! if create_crop_landunit=false)
    if (create_crop_landunit) then
       natpft_size = (numpft + 1) - numcft    ! note that numpft doesn't include bare ground -- thus we add 1
       cft_size    = numcft
    else
       natpft_size = numpft + 1               ! note that numpft doesn't include bare ground -- thus we add 1
       cft_size    = 0
    end if

    natpft_lb = 0
    natpft_ub = natpft_lb + natpft_size - 1
    cft_lb = natpft_ub + 1
    cft_ub = cft_lb + cft_size - 1

    ! (Done initialization that belongs in clm_varpar_init)

  end subroutine initParameters
!------------------------------------------------------------------------------

end module initParametersMod
