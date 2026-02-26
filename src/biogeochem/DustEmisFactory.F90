module DustEmisFactory
    !---------------------------------------------------------------------------
    !
    ! Factory to figure out whihc dust emission method to instantiate
    !
    !---------------------------------------------------------------------------
    use abortutils   , only : endrun
    use shr_log_mod  , only : errMsg => shr_log_errMsg
    use clm_varctl   , only : iulog

    implicit none
    save
    private
    !
    public :: create_dust_emissions  ! create an object of class dust_emis_base

    character(len=*), parameter, private :: sourcefile = &
        __FILE__

contains

    !---------------------------------------------------------------------------

    function create_dust_emissions(bounds, NLFilename) result(dust_emis)
        !---------------------------------------------------------------------------
        ! Create a dust_emission base class objecct
        ! The method implemented depends on namelist input
        !
        ! DESIGN NOTES:   Erik Kluzek 07/15/2024
        !     This implementation is different from for example the Fire Factory functions
        !     that use a direct namelist item with case statements to determine the method.
        !     Here we use logical functions from the shr_dust_emis_mod code. Because shr_dust_emis_mod
        !     is used by both CTSM and CAM I wanted it to be robust with neither CAM nor CTSM
        !     being able to change internal settings so a functional programming design was used
        !     (with function calls that can't change anything inside shr_dust_emis_mod). This is also
        !     why I added a unit-tester for the shr_dust_emis_mod code, so that both CTSM and CAM
        !     can rely on it's behavior.
        !---------------------------------------------------------------------------
        use DustEmisBase      , only : dust_emis_base_type
        use DustEmisZender2003, only : dust_emis_zender2003_type
        use DustEmisLeung2023 , only : dust_emis_leung2023_type
        use decompMod         , only : bounds_type
        use shr_dust_emis_mod , only : is_dust_emis_zender, is_dust_emis_leung
        implicit none
        ! Arguments
        class(dust_emis_base_type), allocatable :: dust_emis
        type(bounds_type), intent(in) :: bounds
        character(len=*),  intent(in) :: NLFilename

        if ( is_dust_emis_zender() )then
           allocate(dust_emis, source=dust_emis_zender2003_type() )

        else if ( is_dust_emis_leung() )then
           allocate(dust_emis, source=dust_emis_leung2023_type() )

        else
           write(iulog,*) 'ERROR: unknown dust_emis_method: ', &
                           errMsg(sourcefile, __LINE__)
           call endrun( "Unrecognized dust_emis_method"  )

        end if

        call dust_emis%Init(bounds, NLFilename)

    end function create_dust_emissions

    !---------------------------------------------------------------------------

end module DustEmisFactory
