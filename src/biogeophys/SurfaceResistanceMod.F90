module SurfaceResistanceMod

#include "shr_assert.h"

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module holding routines for calculation of surface resistances of the different tracers
  ! transported with BeTR. The surface here refers to water and soil, not including canopy
  !
  ! !USES:
  use shr_kind_mod  , only: r8 => shr_kind_r8
  use shr_const_mod , only: SHR_CONST_TKFRZ
  use clm_varctl    , only: iulog
  use SoilStateType , only: soilstate_type
  use WaterStateBulkType, only: waterstatebulk_type 
  use WaterDiagnosticBulkType, only: waterdiagnosticbulk_type 
   use TemperatureType   , only : temperature_type
  implicit none
  save
  private
  integer :: soil_resis_method   !choose the method for soil resistance calculation
  
  integer, parameter :: leepielke_1992 = 0 !
  integer, parameter :: sl_14 = 1 
  
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: calc_soilevap_resis
  public :: do_soilevap_beta, do_soil_resistance_sl14
!  public :: init_soil_resistance
  public :: soil_resistance_readNL

  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  !
  ! !REVISION HISTORY:
  ! 6/25/2013 Created by Jinyun Tang
  !-----------------------------------------------------------------------
  
contains

  !-----------------------------------------------------------------------
!!$  subroutine init_soil_resistance()
!!$   !
!!$   !DESCRIPTIONS
!!$   ! initialize method for soil resis calculation
!!$    !
!!$    ! !USES:
!!$    use abortutils      , only : endrun   
!!$    use fileutils       , only : getavu, relavu
!!$    use spmdMod         , only : mpicom, masterproc
!!$    use shr_mpi_mod     , only : shr_mpi_bcast
!!$    use clm_varctl      , only : iulog, use_bedrock
!!$    use controlMod      , only : NLFilename
!!$    use clm_nlUtilsMod  , only : find_nlgroup_name
!!$
!!$    ! !ARGUMENTS:
!!$    !------------------------------------------------------------------------------
!!$    implicit none
!!$    integer            :: nu_nml                     ! unit for namelist file
!!$    integer            :: nml_error                  ! namelist i/o error flag
!!$    character(*), parameter    :: subName = "('init_soil_resistance')"
!!$
!!$    !-----------------------------------------------------------------------
!!$
!!$! MUST agree with name in namelist and read statement
!!$    namelist /soil_resis_inparm/ soil_resis_method
!!$
!!$    ! Default values for namelist
!!$
!!$   soil_resis_method = sl_14
!!$
!!$    ! Read soil_resis namelist
!!$    if (masterproc) then
!!$       nu_nml = getavu()
!!$       open( nu_nml, file=trim(NLFilename), status='old', iostat=nml_error )
!!$       call find_nlgroup_name(nu_nml, 'soil_resis_inparm', status=nml_error)
!!$       if (nml_error == 0) then
!!$          read(nu_nml, nml=soil_resis_inparm,iostat=nml_error)
!!$          if (nml_error /= 0) then
!!$             call endrun(subname // ':: ERROR reading soil_resis namelist')
!!$          end if
!!$       end if
!!$       close(nu_nml)
!!$       call relavu( nu_nml )
!!$
!!$    endif
!!$
!!$    call shr_mpi_bcast(soil_resis_method, mpicom)
!!$
!!$    if (masterproc) then
!!$       write(iulog,*) ' '
!!$       write(iulog,*) 'soil_resis settings:'
!!$       write(iulog,*) '  soil_resis_method  = ',soil_resis_method
!!$    endif
!!$!scs   
!!$!   soil_resis_method = leepielke_1992
!!$!   soil_resis_method = sl_14
!!$!scs
!!$
!!$  end subroutine init_soil_resistance
   
  !-----------------------------------------------------------------------
  subroutine soil_resistance_readNL(NLFilename)
   !
   !DESCRIPTIONS
   ! Read the namelist for soil resistance method
    !
    ! !USES:
    use abortutils      , only : endrun   
    use fileutils       , only : getavu, relavu
    use spmdMod         , only : mpicom, masterproc
    use shr_mpi_mod     , only : shr_mpi_bcast
    use clm_varctl      , only : iulog
    use clm_nlUtilsMod  , only : find_nlgroup_name

    ! !ARGUMENTS:
    !------------------------------------------------------------------------------
    implicit none
    character(len=*), intent(IN) :: NLFilename ! Namelist filename
    integer                      :: nu_nml     ! unit for namelist file
    integer                      :: nml_error  ! namelist i/o error flag
    character(*), parameter      :: subName = "('init_soil_resistance')"

    !-----------------------------------------------------------------------

! MUST agree with name in namelist and read statement
    namelist /soil_resis_inparm/ soil_resis_method

    ! Default values for namelist

   soil_resis_method = sl_14

    ! Read soil_resis namelist
    if (masterproc) then
       nu_nml = getavu()
       open( nu_nml, file=trim(NLFilename), status='old', iostat=nml_error )
       call find_nlgroup_name(nu_nml, 'soil_resis_inparm', status=nml_error)
       if (nml_error == 0) then
          read(nu_nml, nml=soil_resis_inparm,iostat=nml_error)
          if (nml_error /= 0) then
             call endrun(subname // ':: ERROR reading soil_resis namelist')
          end if
       else
          call endrun(subname // ':: ERROR reading soil_resis namelist')
       end if
       close(nu_nml)
       call relavu( nu_nml )

    endif

    call shr_mpi_bcast(soil_resis_method, mpicom)

    if (masterproc) then
       write(iulog,*) ' '
       write(iulog,*) 'soil_resis settings:'
       write(iulog,*) '  soil_resis_method  = ',soil_resis_method
    endif

  end subroutine soil_resistance_readNL
   
   !------------------------------------------------------------------------------   
   subroutine calc_soilevap_resis(bounds, num_nolakec, filter_nolakec, &
        soilstate_inst, waterstatebulk_inst, waterdiagnosticbulk_inst, temperature_inst)
     !
     ! DESCRIPTIONS
     ! compute the resis factor for soil evaporation calculation
     !
     use shr_kind_mod  , only : r8 => shr_kind_r8     
     use shr_const_mod , only : SHR_CONST_PI  
     use decompMod     , only : bounds_type
     use ColumnType    , only : col
     use LandunitType  , only : lun
     use abortutils    , only : endrun      
     !
     ! !ARGUMENTS:
     implicit none
     type(bounds_type)     , intent(in)    :: bounds    ! bounds   
     integer               , intent(in)    :: num_nolakec
     integer               , intent(in)    :: filter_nolakec(:)
     type(soilstate_type)  , intent(inout) :: soilstate_inst
     type(waterstatebulk_type) , intent(in)    :: waterstatebulk_inst
     type(waterdiagnosticbulk_type) , intent(in)    :: waterdiagnosticbulk_inst
     type(temperature_type), intent(in)    :: temperature_inst
     character(len=32) :: subname = 'calc_soilevap_resis'  ! subroutine name
     associate(                &
          soilbeta =>  soilstate_inst%soilbeta_col  , & ! Output: [real(r8) (:)] factor that reduces ground evaporation
          dsl       =>  soilstate_inst%dsl_col      , & ! Output: [real(r8) (:)] soil dry surface layer thickness
          soilresis =>  soilstate_inst%soilresis_col  & ! Output: [real(r8) (:)] soil evaporative resistance
          )
   
       !select the right method and do the calculation
       select case (soil_resis_method)

       case (leepielke_1992)
          call calc_beta_leepielke1992(bounds, num_nolakec, filter_nolakec, &
               soilstate_inst, waterstatebulk_inst, waterdiagnosticbulk_inst, soilbeta(bounds%begc:bounds%endc))

               case (sl_14)
          call calc_soil_resistance_sl14(bounds, num_nolakec, filter_nolakec, &
               soilstate_inst, waterstatebulk_inst, temperature_inst, &
               dsl(bounds%begc:bounds%endc), soilresis(bounds%begc:bounds%endc))
               case default
          call endrun(subname // ':: a soilevap resis function must be specified!')     
       end select

     end associate

   end subroutine calc_soilevap_resis
   
   !------------------------------------------------------------------------------   
   subroutine calc_beta_leepielke1992(bounds, num_nolakec, filter_nolakec, &
        soilstate_inst, waterstatebulk_inst, waterdiagnosticbulk_inst, soilbeta)
     !
     ! DESCRIPTION
     ! compute the lee-pielke beta factor to scal actual soil evaporation from potential evaporation
     !
     ! USES
     use shr_kind_mod    , only : r8 => shr_kind_r8     
     use shr_const_mod   , only : SHR_CONST_PI
     use shr_log_mod     , only : errMsg => shr_log_errMsg   
     use shr_infnan_mod  , only : nan => shr_infnan_nan, assignment(=)
     use decompMod       , only : bounds_type
     use clm_varcon      , only : denh2o, denice
     use landunit_varcon , only : istice_mec, istwet, istsoil, istcrop
     use column_varcon   , only : icol_roof, icol_sunwall, icol_shadewall
     use column_varcon   , only : icol_road_imperv, icol_road_perv
     use ColumnType      , only : col
     use LandunitType    , only : lun
     !
     implicit none
     type(bounds_type)     , intent(in)    :: bounds    ! bounds   
     integer               , intent(in)    :: num_nolakec
     integer               , intent(in)    :: filter_nolakec(:)
     type(soilstate_type)  , intent(in)    :: soilstate_inst
     type(waterstatebulk_type) , intent(in)    :: waterstatebulk_inst
     type(waterdiagnosticbulk_type) , intent(in)    :: waterdiagnosticbulk_inst
     real(r8)              , intent(inout) :: soilbeta(bounds%begc:bounds%endc)

     !local variables
     real(r8) :: fac, fac_fc, wx      !temporary variables
     integer  :: c, l, fc     !indices

     SHR_ASSERT_ALL((ubound(soilbeta)    == (/bounds%endc/)), errMsg(sourcefile, __LINE__))

     associate(                                              &
          watsat      =>    soilstate_inst%watsat_col      , & ! Input:  [real(r8) (:,:)] volumetric soil water at saturation (porosity)
          watfc       =>    soilstate_inst%watfc_col       , & ! Input:  [real(r8) (:,:)] volumetric soil water at field capacity
          
          h2osoi_ice  =>    waterstatebulk_inst%h2osoi_ice_col , & ! Input:  [real(r8) (:,:)] ice lens (kg/m2)                       
          h2osoi_liq  =>    waterstatebulk_inst%h2osoi_liq_col , & ! Input:  [real(r8) (:,:)] liquid water (kg/m2)                   
          frac_sno    =>    waterdiagnosticbulk_inst%frac_sno_col   , & ! Input:  [real(r8) (:)] fraction of ground covered by snow (0 to 1)
          frac_h2osfc =>    waterdiagnosticbulk_inst%frac_h2osfc_col  & ! Input:  [real(r8) (:)]  fraction of ground covered by surface water (0 to 1)
          )

       do fc = 1,num_nolakec
          c = filter_nolakec(fc)
          l = col%landunit(c)   
          if (lun%itype(l)/=istwet .AND. lun%itype(l)/=istice_mec) then
             if (lun%itype(l) == istsoil .or. lun%itype(l) == istcrop) then
                wx   = (h2osoi_liq(c,1)/denh2o+h2osoi_ice(c,1)/denice)/col%dz(c,1)
                fac  = min(1._r8, wx/watsat(c,1))
                fac  = max( fac, 0.01_r8 )
                !! Lee and Pielke 1992 beta, added by K.Sakaguchi
                if (wx < watfc(c,1) ) then  !when water content of ths top layer is less than that at F.C.
                   fac_fc  = min(1._r8, wx/watfc(c,1))  !eqn5.66 but divided by theta at field capacity
                   fac_fc  = max( fac_fc, 0.01_r8 )
                   ! modify soil beta by snow cover. soilbeta for snow surface is one
                   soilbeta(c) = (1._r8-frac_sno(c)-frac_h2osfc(c)) &
                        *0.25_r8*(1._r8 - cos(SHR_CONST_PI*fac_fc))**2._r8 &
                        + frac_sno(c)+ frac_h2osfc(c)
                else   !when water content of ths top layer is more than that at F.C.
                   soilbeta(c) = 1._r8
                end if
             else if (col%itype(c) == icol_road_perv) then
                soilbeta(c) = 0._r8
             else if (col%itype(c) == icol_sunwall .or. col%itype(c) == icol_shadewall) then
                soilbeta(c) = 0._r8          
             else if (col%itype(c) == icol_roof .or. col%itype(c) == icol_road_imperv) then
                soilbeta(c) = 0._r8
             endif
          else
             soilbeta(c) =   1._r8
          endif
       enddo

     end associate

   end subroutine calc_beta_leepielke1992
   
   !------------------------------------------------------------------------------   
   function do_soilevap_beta()result(lres)
     !
     !DESCRIPTION
     ! return true if the moisture stress for soil evaporation is computed as beta factor
     ! otherwise false
     implicit none
     logical :: lres

     if(soil_resis_method==leepielke_1992)then
        lres=.true.
     else
        lres=.false.
     endif
     return

   end function do_soilevap_beta

  !------------------------------------------------------------------------------   
   subroutine calc_soil_resistance_sl14(bounds, num_nolakec, filter_nolakec, &
        soilstate_inst, waterstatebulk_inst, temperature_inst, dsl, soilresis)
     !
     ! DESCRIPTION
     ! compute the lee-pielke beta factor to scal actual soil evaporation from potential evaporation
     !
     ! USES
     use shr_kind_mod    , only : r8 => shr_kind_r8     
     use shr_const_mod   , only : SHR_CONST_PI
     use shr_log_mod     , only : errMsg => shr_log_errMsg   
     use shr_infnan_mod  , only : nan => shr_infnan_nan, assignment(=)
     use decompMod       , only : bounds_type
     use clm_varcon      , only : denh2o, denice
     use landunit_varcon , only : istice_mec, istwet, istsoil, istcrop
     use column_varcon   , only : icol_roof, icol_sunwall, icol_shadewall
     use column_varcon   , only : icol_road_imperv, icol_road_perv
     use ColumnType      , only : col
     use LandunitType    , only : lun
     !
     implicit none
     type(bounds_type)     , intent(in)    :: bounds    ! bounds   
     integer               , intent(in)    :: num_nolakec
     integer               , intent(in)    :: filter_nolakec(:)
     type(soilstate_type)  , intent(in)    :: soilstate_inst
     type(waterstatebulk_type) , intent(in)    :: waterstatebulk_inst
     type(temperature_type), intent(in)    :: temperature_inst
     real(r8)              , intent(inout) :: dsl(bounds%begc:bounds%endc)
     real(r8)              , intent(inout) :: soilresis(bounds%begc:bounds%endc)

   !local variables
     real(r8) :: aird, eps, dg, d0, vwc_liq
     real(r8) :: eff_por_top
     integer  :: c, l, fc     !indices
     
     SHR_ASSERT_ALL((ubound(dsl)    == (/bounds%endc/)), errMsg(sourcefile, __LINE__))
     SHR_ASSERT_ALL((ubound(soilresis)    == (/bounds%endc/)), errMsg(sourcefile, __LINE__))

     associate(                                              &
          dz                =>    col%dz                             , & ! Input:  [real(r8) (:,:) ]  layer thickness (m)                             
          watsat            =>    soilstate_inst%watsat_col      , & ! Input:  [real(r8) (:,:)] volumetric soil water at saturation (porosity)
          bsw               =>    soilstate_inst%bsw_col             , & ! Input:  [real(r8) (:,:) ]  Clapp and Hornberger "b"                        
          sucsat            =>    soilstate_inst%sucsat_col          , & ! Input:  [real(r8) (:,:) ]  minimum soil suction (mm)                       
!          eff_porosity      =>    soilstate_inst%eff_porosity_col    , & ! Input:  [real(r8) (:,:) ]  effective porosity = porosity - vol_ice         
          t_soisno          =>    temperature_inst%t_soisno_col      ,  & ! Input:  [real(r8) (:,:) ]  soil temperature (Kelvin)                       
         
          h2osoi_ice        =>    waterstatebulk_inst%h2osoi_ice_col , & ! Input:  [real(r8) (:,:)] ice lens (kg/m2)                       
          h2osoi_liq        =>    waterstatebulk_inst%h2osoi_liq_col  & ! Input:  [real(r8) (:,:)] liquid water (kg/m2)                   
          )

   do fc = 1,num_nolakec
      c = filter_nolakec(fc)
      l = col%landunit(c)  
      if (lun%itype(l)/=istwet .AND. lun%itype(l)/=istice_mec) then
         if (lun%itype(l) == istsoil .or. lun%itype(l) == istcrop) then
            vwc_liq = max(h2osoi_liq(c,1),1.0e-6_r8)/(dz(c,1)*denh2o)
! eff_porosity not calculated til SoilHydrology
             eff_por_top = max(0.01_r8,watsat(c,1)-min(watsat(c,1), h2osoi_ice(c,1)/(dz(c,1)*denice)))

! calculate diffusivity and air free pore space
            aird = watsat(c,1)*(sucsat(c,1)/1.e7_r8)**(1./bsw(c,1))
            d0 = 2.12e-5*(t_soisno(c,1)/273.15)**1.75 ![Bitelli et al., JH, 08]
            eps = watsat(c,1) - aird
            dg = eps*d0*(eps/watsat(c,1))**(3._r8/max(3._r8,bsw(c,1)))
            
!      dsl(c) = dzmm(c,1)*max(0.001_r8,(0.8*eff_porosity(c,1) - vwc_liq)) &
! try arbitrary scaling (not top layer thickness)
!            dsl(c) = 15._r8*max(0.001_r8,(0.8*eff_porosity(c,1) - vwc_liq)) &
            dsl(c) = 15._r8*max(0.001_r8,(0.8*eff_por_top - vwc_liq)) &
                 !           /max(0.001_r8,(watsat(c,1)- aird))
                 /max(0.001_r8,(0.8*watsat(c,1)- aird))
            
            dsl(c)=max(dsl(c),0._r8)
            dsl(c)=min(dsl(c),200._r8)
            
            soilresis(c) = dsl(c)/(dg*eps*1.e3) + 20._r8
            soilresis(c) = min(1.e6_r8,soilresis(c))

         else if (col%itype(c) == icol_road_perv) then
            soilresis(c) = 1.e6_r8
         else if (col%itype(c) == icol_sunwall .or. col%itype(c) == icol_shadewall) then
            soilresis(c) = 1.e6_r8          
         else if (col%itype(c) == icol_roof .or. col%itype(c) == icol_road_imperv) then
            soilresis(c) = 1.e6_r8
         endif   
      else
         soilresis(c) =   0._r8
      endif
   enddo   
   end associate
   end subroutine calc_soil_resistance_sl14

   !------------------------------------------------------------------------------   
   function do_soil_resistance_sl14()result(lres)
     !
     !DESCRIPTION
     ! return true if the soil evaporative resistance is computed using a DSL
     ! otherwise false
     implicit none
     logical :: lres

     if(soil_resis_method==sl_14)then
        lres=.true.
     else
        lres=.false.
     endif
     return

   end function do_soil_resistance_sl14

end module SurfaceResistanceMod
