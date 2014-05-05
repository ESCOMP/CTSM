module EDSetValuesMod
  !-----------------------------------------------------------------------
  ! !MODULE: EDsetValuesMod
  !
  ! !DESCRIPTION:
  ! set values in EDpcf instance (of type 
  !
  ! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8

  implicit none
  save
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: EDSetPcf
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine EDSetPcf( val, EDpcf )
    !
    ! !ARGUMENTS:
    use clm_varctl, only : use_ed  
    use EDclmType , only : EDpft_cflux_type

    implicit none

    real(r8), intent(in) :: val
    type (EDpft_cflux_type), intent(inout) :: EDpcf

    EDpcf%trimming(:)          = val
    EDpcf%canopy_spread(:)     = val
    EDpcf%PFTbiomass(:,:)      = val
    EDpcf%PFTleafbiomass(:,:)  = val
    EDpcf%PFTstorebiomass(:,:) = val
    EDpcf%PFTnindivs(:,:)      = val
    EDpcf%efpot(:)             = val
    EDpcf%rb(:)                = val

  end subroutine EDSetPcf

end module EDSetValuesMod
