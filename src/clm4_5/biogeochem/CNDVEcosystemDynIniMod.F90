module CNDVEcosystemDyniniMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! CNDV related initializations
  !
  ! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
  use decompMod   , only : bounds_type
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public  :: CNDVEcosystemDynini ! CNDV related initializations
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine CNDVEcosystemDynini(bounds)
    !
    ! !DESCRIPTION:
    ! CNDV related initializations
    !
    ! !USES:
    use clmtype
    use shr_const_mod, only : SHR_CONST_PI, SHR_CONST_TKFRZ
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds  ! bounds
    !
    ! !LOCAL VARIABLES:
    integer  :: p           ! pft index
    !-----------------------------------------------------------------------

    ! Some of the following came from LPJ subroutine initgrid

    do p = bounds%begp,bounds%endp
       pdgvs%present(p)   = .false.
       pdgvs%crownarea(p) = 0._r8
       pdgvs%nind(p)      = 0._r8
       pcs%leafcmax(p)    = 0._r8
       pdgvs%t_mo_min(p)  = 1.0e+36_r8
       pdgvs%agdd20(p)   = 0._r8
       pdgvs%tmomin20(p) = SHR_CONST_TKFRZ - 5._r8 !initialize this way for Phenology code
    end do

  end subroutine CNDVEcosystemDynini

end module CNDVEcosystemDyniniMod
