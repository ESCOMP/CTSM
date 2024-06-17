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
        !---------------------------------------------------------------------------
        use DustEmisBase      , only : dust_emis_base_type
        use DustEmisZender2003, only : dust_emis_zender2003_type
        use clm_varctl        , only : dust_emis_method
        use decompMod         , only : bounds_type
        use shr_kind_mod      , only : CL => shr_kind_cl
        implicit none
        ! Arguments
        class(dust_emis_base_type), allocatable :: dust_emis
        type(bounds_type), intent(in) :: bounds
        character(len=*),  intent(in) :: NLFilename
        ! Local variables
        character(len=CL) :: method

        method = dust_emis_method

        select case ( trim(method) )

        case( "Zender_2003" )
           allocate(dust_emis, source=dust_emis_zender2003_type() )
        
        ! This will be added when the Leung2023 comes in
        !case( "Leung_2023" )
        !  allocate(dust_emis, source=dust_emis_zender2003_type() )
        case default
           write(iulog,*) 'ERROR: unknown dust_emis_method: ', method, &
                           errMsg(sourcefile, __LINE__)
           call endrun( "Unrecognized dust_emis_method"  )

        end select

        call dust_emis%Init(bounds, NLFilename)

    end function create_dust_emissions

    !---------------------------------------------------------------------------

end module DustEmisFactory
