module initch4Mod
#ifdef LCH4

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: initch4Mod
!
! !DESCRIPTION:
! Contains time constant (and flux / diagnostic vars) and time-varying (state vars) initialization code for CH4 scheme.
! 
!
! !PUBLIC TYPES:
  implicit none
  save
  private
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: initch4 ! driver
!
! !PRIVATE MEMBER FUNCTIONS:
  private :: initTimeConst_ch4        ! Set constant parameters.
  private :: makearbinit_ch4          ! Set time-variable parameters for spin up.
!
! !REVISION HISTORY:
! Created by Zack Subin, 2009.
!
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: initch4
!
! !INTERFACE:
  subroutine initch4( arbinit )
!
! !DESCRIPTION:
! Calls initTimeConst_ch4.
! Calls makearbinit_ch4 with logical arbinit. If arbinit == .true. OR if initial conditions file
! does not contain methane and oxygen concentrations, then initializes time varying values. This
! allows back-compatibility with initial condition files that have not been spun up with the new
! lake code. In future versions, this could be phased out.
!
! !USES:
  !use ch4varcon, only : ch4conrd
!
! !ARGUMENTS:
    implicit none
!
    logical, intent(in) ::  arbinit ! Whether mkarbinit has been called.
!
! !CALLED FROM:
! subroutine initialize2 in module initializeMod
!
! !REVISION HISTORY:
! Created by Zack Subin, 2009.
!
! !LOCAL VARIABLES:
!
!
!EOP
!

    call initTimeConst_ch4()
! Attn EK
! For now
    call makearbinit_ch4(arbinit)
! For future versions always using initial condition files spun up with the new ch4 code:
!    if (arbinit) call makearbinit_ch4(arbinit)

  end subroutine initch4

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: makearbinit_ch4
!
! !INTERFACE:
  subroutine makearbinit_ch4( arbinit )
!
! !DESCRIPTION:
! If arbinit == .true., or if methane & oxygen concentrations (or lake soil org matter) 
! have not been initialized, then sets time
! varying values.
! Initializes the following time varying variables:
! 
! conc_ch4_sat, conc_ch4_unsat, conc_o2_sat, conc_o2_unsat, lake_soilc, o2stress, finunduated
!
! !USES:
    use shr_kind_mod , only : r8 => shr_kind_r8
    use clmtype
    use clm_varpar   , only : nlevsoi, nlevgrnd
    use clm_varcon   , only : istsoil, istdlak, spval, istcrop
    use spmdMod      , only : masterproc
    use decompMod    , only : get_proc_bounds
    use clm_varctl   , only : iulog
!
! !ARGUMENTS:
    implicit none
    logical, intent(in) :: arbinit ! Whether mkarbinit has been called.
!
! !CALLED FROM:
! subroutine initialize in module initializeMod
!
! !REVISION HISTORY:
! Created by Zack Subin, 2009.
!
! !LOCAL VARIABLES:
!EOP
    integer :: j,l,c,p      ! indices
    integer :: begp, endp   ! per-proc beginning and ending pft indices
    integer :: begc, endc   ! per-proc beginning and ending column indices
    integer :: begl, endl   ! per-proc beginning and ending landunit indices
    integer :: begg, endg   ! per-proc gridcell ending gridcell indices
!-----------------------------------------------------------------------

    if ( masterproc ) write (iulog,*) 'Setting initial data to non-spun up values for CH4 Mod,', &
                                  'if no inicFile or no valid values for concentrations,', &
                                  'for CH4 Model.'
   associate(& 
   ltype                               =>    lun%itype                                   , & ! Input:  [integer (:)]  landunit type                                      
   clandunit                           =>   col%landunit                                 , & ! Input:  [integer (:)]  landunit index associated with each column         
   conc_ch4_sat                        =>    cch4%conc_ch4_sat                           , & ! Output: [real(r8) (:,:)]  CH4 conc in each soil layer (mol/m3) (nlevsoi)  
   conc_ch4_unsat                      =>    cch4%conc_ch4_unsat                         , & ! Output: [real(r8) (:,:)]  CH4 conc in each soil layer (mol/m3) (nlevsoi)  
   conc_o2_sat                         =>    cch4%conc_o2_sat                            , & ! Output: [real(r8) (:,:)]  O2 conc in each soil layer (mol/m3) (nlevsoi)   
   conc_o2_unsat                       =>    cch4%conc_o2_unsat                          , & ! Output: [real(r8) (:,:)]  O2 conc in each soil layer (mol/m3) (nlevsoi)   
   lake_soilc                          =>    cch4%lake_soilc                             , & ! Output: [real(r8) (:,:)]  total soil organic matter found in level (g C / m^3) (nlevsoi)
   cellorg                             =>    cps%cellorg                                 , & ! Input:  [real(r8) (:,:)]  column 3D organic matter (kg/m^3, 58% by mass carbon) (nlevsoi)
   qflx_surf_lag                       =>    cch4%qflx_surf_lag                          , & ! Output: [real(r8) (:)]  time-lagged surface runoff (mm H2O /s)            
   finundated_lag                      =>    cch4%finundated_lag                         , & ! Output: [real(r8) (:)]  time-lagged fractional inundated area             
   o2stress_sat                        =>    cch4%o2stress_sat                           , & ! Output: [real(r8) (:,:)]  Ratio of oxygen available to that demanded by roots, aerobes, & methanotrophs (nlevsoi)
   o2stress_unsat                      =>    cch4%o2stress_unsat                         , & ! Output: [real(r8) (:,:)]  Ratio of oxygen available to that demanded by roots, aerobes, & methanotrophs (nlevsoi)
   finundated                          =>    cws%finundated                              , & ! Output: [real(r8) (:)]  inundated gridcell fractional area (excluding dedicated wetland columns)
   fsat_bef                            =>    cch4%fsat_bef                               , & ! Input:  [real(r8) (:)]  finundated from previous timestep                 
   layer_sat_lag                       =>    cch4%layer_sat_lag                            & ! Output: [real(r8) (:,:)]  Lagged saturation status of soil layer in the unsaturated zone (1 = sat)
   )


    ! Determine subgrid bounds on this processor

    call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)

    do c = begc,endc

       l = clandunit(c)

       if (ltype(l) == istsoil .or. ltype(l) == istcrop) then
          do j=1,nlevsoi
             if (conc_ch4_sat(c,j) == spval .or. arbinit)   conc_ch4_sat(c,j)   = 0._r8
             if (conc_ch4_unsat(c,j) == spval .or. arbinit) conc_ch4_unsat(c,j) = 0._r8
             if (conc_o2_sat(c,j) == spval .or. arbinit)    conc_o2_sat(c,j)    = 0._r8
             if (conc_o2_unsat(c,j) == spval .or. arbinit)  conc_o2_unsat(c,j)  = 0._r8
             if (o2stress_sat(c,j) == spval .or. arbinit) o2stress_sat(c,j) = 1._r8
             if (o2stress_unsat(c,j) == spval .or. arbinit) o2stress_unsat(c,j) = 1._r8
             if (layer_sat_lag(c,j) == spval .or. arbinit) layer_sat_lag(c,j) = 1._r8
          end do
          if (qflx_surf_lag(c) == spval .or. arbinit) qflx_surf_lag(c) = 0._r8
          if (finundated_lag(c) == spval .or. arbinit) finundated_lag(c) = 0._r8
          ! finundated will be used to calculate soil decomposition if anoxia is used
          if (fsat_bef(c) == spval .or. arbinit) then
             finundated(c) = 0._r8
          else
             finundated(c) = fsat_bef(c)
          end if
       else if (ltype(l) == istdlak) then
          do j=1,nlevsoi
             if (conc_ch4_sat(c,j) == spval .or. arbinit)   conc_ch4_sat(c,j)   = 0._r8
             if (conc_o2_sat(c,j) == spval .or. arbinit)    conc_o2_sat(c,j)    = 0._r8
             if (lake_soilc(c,j) == spval .or. arbinit) lake_soilc(c,j) = 580._r8*cellorg(c,j)
             ! Need to convert from kg/m^3 organic matter to g C / m^3 (org matter is defined to be 58% C)
          end do
       end if 

       ! Set values for all columns equal to zero below nlevsoi
       do j=nlevsoi+1,nlevgrnd
          conc_ch4_sat(c,j) = 0._r8
          conc_ch4_unsat(c,j) = 0._r8
          conc_o2_sat(c,j) = 0._r8
          conc_o2_unsat(c,j) = 0._r8
          lake_soilc(c,j) = 0._r8
          o2stress_sat(c,j) = 1._r8
          o2stress_unsat(c,j) = 1._r8
          layer_sat_lag(c,j) = 1._r8
       end do
       
    end do


    end associate 
   end subroutine makearbinit_ch4


!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: initTimeConst_ch4
!
! !INTERFACE:
subroutine initTimeConst_ch4
!
! !DESCRIPTION:
! Initialize variables for ch4 code that will not be input from restart/inic file. Also set values for
! inactive CH4 columns to spval so that they will not be averaged in history file.
!
! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use clmtype
  use decompMod   , only : get_proc_bounds, get_proc_global
  use clm_varpar  , only : nlevsoi, ngases, nlevgrnd, nlevdecomp
  use clm_varcon  , only : istsoil, istdlak, spval, istcrop
  use clm_varctl  , only : iulog
  use spmdMod ,     only : masterproc
  use ch4varcon   , only : allowlakeprod
!
! !ARGUMENTS:
  implicit none
!
! !CALLED FROM:
! subroutine initialize2 in module initializeMod.
!
! !REVISION HISTORY:
! 9/09, Zack Subin
!
! !LOCAL VARIABLES:
!EOP
  integer  :: g,l,c,p,j        ! indices
  integer  :: begp, endp       ! per-proc beginning and ending pft indices
  integer  :: begc, endc       ! per-proc beginning and ending column indices
  integer  :: begl, endl       ! per-proc beginning and ending landunit indices
  integer  :: begg, endg       ! per-proc gridcell ending gridcell indices

   associate(& 
   watsat                              =>    cps%watsat                                  , & ! Input:  [real(r8) (:,:)]  volumetric soil water at saturation (porosity)  
   ch4_surf_diff_sat                   =>    cch4%ch4_surf_diff_sat                      , & ! Output: [real(r8) (:)]  CH4 surface flux (mol/m2/s)                       
   ch4_surf_diff_unsat                 =>    cch4%ch4_surf_diff_unsat                    , & ! Output: [real(r8) (:)]  CH4 surface flux (mol/m2/s)                       
   ch4_surf_diff_lake                  =>    cch4%ch4_surf_diff_lake                     , & ! Output: [real(r8) (:)]  CH4 surface flux (mol/m2/s)                       
   ch4_surf_aere_sat                   =>    cch4%ch4_surf_aere_sat                      , & ! Output: [real(r8) (:)]  Total column CH4 aerenchyma (mol/m2/s)            
   ch4_surf_aere_unsat                 =>    cch4%ch4_surf_aere_unsat                    , & ! Output: [real(r8) (:)]  Total column CH4 aerenchyma (mol/m2/s)            
   ch4_surf_ebul_sat                   =>    cch4%ch4_surf_ebul_sat                      , & ! Output: [real(r8) (:)]  CH4 ebullition to atmosphere (mol/m2/s)           
   ch4_surf_ebul_unsat                 =>    cch4%ch4_surf_ebul_unsat                    , & ! Output: [real(r8) (:)]  CH4 ebullition to atmosphere (mol/m2/s)           
   ch4_surf_ebul_lake                  =>    cch4%ch4_surf_ebul_lake                     , & ! Output: [real(r8) (:)]  CH4 ebullition to atmosphere (mol/m2/s)           
   ch4_oxid_depth_sat                  =>    cch4%ch4_oxid_depth_sat                     , & ! Output: [real(r8) (:,:)]  CH4 consumption rate via oxidation in each soil layer (mol/m3/s) (nlevsoi)
   ch4_oxid_depth_unsat                =>    cch4%ch4_oxid_depth_unsat                   , & ! Output: [real(r8) (:,:)]  CH4 consumption rate via oxidation in each soil layer (mol/m3/s) (nlevsoi)
   ch4_oxid_depth_lake                 =>    cch4%ch4_oxid_depth_lake                    , & ! Output: [real(r8) (:,:)]  CH4 consumption rate via oxidation in each soil layer (mol/m3/s) (nlevsoi)
   ch4_prod_depth_sat                  =>    cch4%ch4_prod_depth_sat                     , & ! Output: [real(r8) (:,:)]  production of CH4 in each soil layer (nlevsoi) (mol/m3/s)
   ch4_prod_depth_unsat                =>    cch4%ch4_prod_depth_unsat                   , & ! Output: [real(r8) (:,:)]  production of CH4 in each soil layer (nlevsoi) (mol/m3/s)
   ch4_prod_depth_lake                 =>    cch4%ch4_prod_depth_lake                    , & ! Output: [real(r8) (:,:)]  production of CH4 in each soil layer (nlevsoi) (mol/m3/s)
   ch4_ebul_total_sat                  =>    cch4%ch4_ebul_total_sat                     , & ! Output: [real(r8) (:)]  Total column CH4 ebullition (mol/m2/s)            
   ch4_ebul_total_unsat                =>    cch4%ch4_ebul_total_unsat                   , & ! Output: [real(r8) (:)]  Total column CH4 ebullition (mol/m2/s)            
   ch4_ebul_depth_sat                  =>    cch4%ch4_ebul_depth_sat                     , & ! Output: [real(r8) (:,:)]  CH4 loss rate via ebullition in each soil layer (mol/m3/s) (nlevsoi)
   ch4_ebul_depth_unsat                =>    cch4%ch4_ebul_depth_unsat                   , & ! Output: [real(r8) (:,:)]  CH4 loss rate via ebullition in each soil layer (mol/m3/s) (nlevsoi)
   ch4_aere_depth_sat                  =>    cch4%ch4_aere_depth_sat                     , & ! Output: [real(r8) (:,:)]  CH4 loss rate via aerenchyma in each soil layer (mol/m3/s) (nlevsoi)
   ch4_aere_depth_unsat                =>    cch4%ch4_aere_depth_unsat                   , & ! Output: [real(r8) (:,:)]  CH4 loss rate via aerenchyma in each soil layer (mol/m3/s) (nlevsoi)
   ch4_tran_depth_sat                  =>    cch4%ch4_tran_depth_sat                     , & ! Output: [real(r8) (:,:)]  CH4 loss rate via transpiration in each soil layer (mol/m3/s) (nlevsoi)
   ch4_tran_depth_unsat                =>    cch4%ch4_tran_depth_unsat                   , & ! Output: [real(r8) (:,:)]  CH4 loss rate via transpiration in each soil layer (mol/m3/s) (nlevsoi)
   co2_aere_depth_sat                  =>    cch4%co2_aere_depth_sat                     , & ! Output: [real(r8) (:,:)]  CO2 loss rate via aerenchyma in each soil layer (mol/m3/s) (nlevsoi)
   co2_aere_depth_unsat                =>    cch4%co2_aere_depth_unsat                   , & ! Output: [real(r8) (:,:)]  CO2 loss rate via aerenchyma in each soil layer (mol/m3/s) (nlevsoi)
   o2_oxid_depth_sat                   =>    cch4%o2_oxid_depth_sat                      , & ! Output: [real(r8) (:,:)]  O2 consumption rate via oxidation in each soil layer (mol/m3/s) (nlevsoi)
   o2_oxid_depth_unsat                 =>    cch4%o2_oxid_depth_unsat                    , & ! Output: [real(r8) (:,:)]  O2 consumption rate via oxidation in each soil layer (mol/m3/s) (nlevsoi)
   o2_decomp_depth_sat                 =>    cch4%o2_decomp_depth_sat                    , & ! Output: [real(r8) (:,:)]  O2 consumption during decomposition in each soil layer (nlevsoi) (mol/m3/s)
   o2_decomp_depth_unsat               =>    cch4%o2_decomp_depth_unsat                  , & ! Output: [real(r8) (:,:)]  O2 consumption during decomposition in each soil layer (nlevsoi) (mol/m3/s)
   o2_aere_depth_sat                   =>    cch4%o2_aere_depth_sat                      , & ! Output: [real(r8) (:,:)]  O2 gain rate via aerenchyma in each soil layer (mol/m3/s) (nlevsoi)
   o2_aere_depth_unsat                 =>    cch4%o2_aere_depth_unsat                    , & ! Output: [real(r8) (:,:)]  O2 gain rate via aerenchyma in each soil layer (mol/m3/s) (nlevsoi)
   co2_decomp_depth_sat                =>    cch4%co2_decomp_depth_sat                   , & ! Output: [real(r8) (:,:)]  CO2 production during decomposition in each soil layer (nlevsoi) (mol/m3/s)
   co2_decomp_depth_unsat              =>    cch4%co2_decomp_depth_unsat                 , & ! Output: [real(r8) (:,:)]  CO2 production during decomposition in each soil layer (nlevsoi) (mol/m3/s)
   co2_oxid_depth_sat                  =>    cch4%co2_oxid_depth_sat                     , & ! Output: [real(r8) (:,:)]  CO2 production rate via oxidation in each soil layer (mol/m3/s) (nlevsoi)
   co2_oxid_depth_unsat                =>    cch4%co2_oxid_depth_unsat                   , & ! Output: [real(r8) (:,:)]  CO2 production rate via oxidation in each soil layer (mol/m3/s) (nlevsoi)
   ch4_dfsat_flux                      =>    cch4%ch4_dfsat_flux                         , & ! Output: [real(r8) (:)]  CH4 flux to atm due to decreasing fsat (kg C/m^2/s) [+]
   zwt_ch4_unsat                       =>    cch4%zwt_ch4_unsat                          , & ! Output: [real(r8) (:)]  depth of water table for unsaturated fraction (m) 
   grnd_ch4_cond                       =>    pps_a%grnd_ch4_cond                         , & ! Output: [real(r8) (:)]  tracer conductance for boundary layer [m/s]       
   conc_ch4_lake                       =>    cch4%conc_ch4_lake                          , & ! Output: [real(r8) (:,:)]  CH4 conc in each soil layer (mol/m3) (nlevsoi)  
   conc_o2_lake                        =>    cch4%conc_o2_lake                           , & ! Output: [real(r8) (:,:)]  O2 conc in each soil layer (mol/m3) (nlevsoi)   
   finundated                          =>    cws%finundated                              , & ! Output: [real(r8) (:)]  inundated gridcell fractional area                
   fphr                                =>    cch4%fphr                                   , & ! Output: [real(r8) (:,:)]  fraction of potential HR                        
   sif                                 =>    cch4%sif                                    , & ! Output: [real(r8) (:)]  (unitless) ratio applied to sat. prod. to account for seasonal inundation
   rootfr                              =>    pps_a%rootfr                                , & ! Output: [real(r8) (:,:)]  column-averaged root fraction                   
   o2stress_unsat                      =>    cch4%o2stress_unsat                         , & ! Output: [real(r8) (:,:)]  Ratio of oxygen available to that demanded by roots, aerobes, & methanotrophs (nlevsoi)
   o2stress_sat                        =>    cch4%o2stress_sat                           , & ! Output: [real(r8) (:,:)]  Ratio of oxygen available to that demanded by roots, aerobes, & methanotrophs (nlevsoi)
   ch4stress_unsat                     =>    cch4%ch4stress_unsat                        , & ! Output: [real(r8) (:,:)]  Ratio of methane available to the total per-timestep methane sinks (nlevsoi)
   ch4stress_sat                       =>    cch4%ch4stress_sat                          , & ! Output: [real(r8) (:,:)]  Ratio of methane available to the total per-timestep methane sinks (nlevsoi)
   totcolch4                           =>    cch4%totcolch4                              , & ! Output: [real(r8) (:)]  total methane in soil column (g C / m^2)          
   conc_ch4_sat                        =>    cch4%conc_ch4_sat                           , & ! Output: [real(r8) (:,:)]  CH4 conc in each soil layer (mol/m3) (nlevsoi)  
   conc_o2_sat                         =>    cch4%conc_o2_sat                            , & ! Output: [real(r8) (:,:)]  O2 conc in each soil layer (mol/m3) (nlevsoi)   
   lake_soilc                          =>    cch4%lake_soilc                             , & ! Output: [real(r8) (:,:)]  total soil organic matter found in level (g C / m^3) (nlevsoi)
   cellorg                             =>    cps%cellorg                                 , & ! Input:  [real(r8) (:,:)]  column 3D organic matter (kg/m^3, 58% by mass carbon) (nlevsoi)
   ltype                               =>    lun%itype                                   , & ! Input:  [integer (:)]  landunit type index                                
   clandunit                           =>    col%landunit                                , & ! Input:  [integer (:)]  landunit index of column                           
   cgridcell                           =>    col%gridcell                                  & ! Input:  [integer (:)]  gridcell index of column                           
   )

  if (masterproc) write (iulog,*) 'Attempting to initialize non-state variables for CH4 Mod'

  call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)

  do c = begc,endc
     l = clandunit(c)

     ! First set levels from nlevsoi+1 to nlevgrnd = 0

     ch4_prod_depth_sat(c,nlevsoi+1:nlevgrnd) = 0._r8
     ch4_prod_depth_unsat(c,nlevsoi+1:nlevgrnd) = 0._r8
     ch4_prod_depth_lake(c,nlevsoi+1:nlevgrnd) = 0._r8
     ch4_oxid_depth_sat(c,nlevsoi+1:nlevgrnd) = 0._r8
     ch4_oxid_depth_unsat(c,nlevsoi+1:nlevgrnd) = 0._r8
     ch4_oxid_depth_lake(c,nlevsoi+1:nlevgrnd) = 0._r8
     o2_oxid_depth_sat(c,nlevsoi+1:nlevgrnd) = 0._r8
     o2_oxid_depth_unsat(c,nlevsoi+1:nlevgrnd) = 0._r8
     o2_decomp_depth_sat(c,nlevsoi+1:nlevgrnd) = 0._r8
     o2_decomp_depth_unsat(c,nlevsoi+1:nlevgrnd) = 0._r8
     o2_aere_depth_sat(c,nlevsoi+1:nlevgrnd) = 0._r8
     o2_aere_depth_unsat(c,nlevsoi+1:nlevgrnd) = 0._r8
     co2_decomp_depth_sat(c,nlevsoi+1:nlevgrnd) = 0._r8
     co2_decomp_depth_unsat(c,nlevsoi+1:nlevgrnd) = 0._r8
     co2_oxid_depth_sat(c,nlevsoi+1:nlevgrnd) = 0._r8
     co2_oxid_depth_unsat(c,nlevsoi+1:nlevgrnd) = 0._r8
     ch4_aere_depth_sat(c,nlevsoi+1:nlevgrnd) = 0._r8
     ch4_aere_depth_unsat(c,nlevsoi+1:nlevgrnd) = 0._r8
     ch4_tran_depth_sat(c,nlevsoi+1:nlevgrnd) = 0._r8
     ch4_tran_depth_unsat(c,nlevsoi+1:nlevgrnd) = 0._r8
     co2_aere_depth_sat(c,nlevsoi+1:nlevgrnd) = 0._r8
     co2_aere_depth_unsat(c,nlevsoi+1:nlevgrnd) = 0._r8
     ch4_ebul_depth_sat(c,nlevsoi+1:nlevgrnd) = 0._r8
     ch4_ebul_depth_unsat(c,nlevsoi+1:nlevgrnd) = 0._r8
     conc_ch4_lake(c,nlevsoi+1:nlevgrnd) = 0._r8
     conc_o2_lake(c,nlevsoi+1:nlevgrnd) = 0._r8
     rootfr(c,nlevsoi+1:nlevgrnd) = 0._r8
     ch4stress_unsat(c,nlevsoi+1:nlevgrnd) = 0._r8
     ch4stress_sat(c,nlevsoi+1:nlevgrnd) = 0._r8
     fphr(c,nlevdecomp+1:nlevgrnd) = 0._r8


     if (ltype(l) == istsoil .or. ltype(l) == istcrop) then

        conc_ch4_lake(c,:) = spval
        conc_o2_lake(c,:) = spval
        ch4_surf_diff_lake(c) = spval
        ch4_surf_ebul_lake(c) = spval
        ch4_prod_depth_lake(c,:) = spval
        ch4_oxid_depth_lake(c,:) = spval

     else if (ltype(l) == istdlak .and. allowlakeprod) then

        ch4_prod_depth_unsat(c,:) = spval
        ch4_oxid_depth_unsat(c,:) = spval
        o2_oxid_depth_unsat(c,:) = spval
        o2_decomp_depth_unsat(c,:) = spval
        o2_aere_depth_unsat(c,:) = spval
        co2_decomp_depth_unsat(c,:) = spval
        co2_oxid_depth_unsat(c,:) = spval
        ch4_aere_depth_unsat(c,:) = spval
        ch4_tran_depth_unsat(c,:) = spval
        co2_aere_depth_unsat(c,:) = spval
        ch4_surf_aere_unsat(c) = spval
        ch4_ebul_depth_unsat(c,:) = spval
        ch4_ebul_total_unsat(c) = spval
        ch4_surf_ebul_unsat(c) = spval
        ch4_surf_diff_unsat(c) = spval
        ch4_dfsat_flux(c) = spval
        zwt_ch4_unsat(c) = spval
        finundated(c) = spval
        fphr(c,:) = spval
        sif(c)       = spval
        rootfr(c,:) = spval
        o2stress_unsat(c,:) = spval
        ch4stress_unsat(c,:) = spval

     else  ! Inactive CH4 columns

        ch4_prod_depth_sat(c,:) = spval
        ch4_prod_depth_unsat(c,:) = spval
        ch4_prod_depth_lake(c,:) = spval
        ch4_oxid_depth_sat(c,:) = spval
        ch4_oxid_depth_unsat(c,:) = spval
        ch4_oxid_depth_lake(c,:) = spval
        o2_oxid_depth_sat(c,:) = spval
        o2_oxid_depth_unsat(c,:) = spval
        o2_decomp_depth_sat(c,:) = spval
        o2_decomp_depth_unsat(c,:) = spval
        o2_aere_depth_sat(c,:) = spval
        o2_aere_depth_unsat(c,:) = spval
        co2_decomp_depth_sat(c,:) = spval
        co2_decomp_depth_unsat(c,:) = spval
        co2_oxid_depth_sat(c,:) = spval
        co2_oxid_depth_unsat(c,:) = spval
        ch4_aere_depth_sat(c,:) = spval
        ch4_aere_depth_unsat(c,:) = spval
        ch4_tran_depth_sat(c,:) = spval
        ch4_tran_depth_unsat(c,:) = spval
        co2_aere_depth_sat(c,:) = spval
        co2_aere_depth_unsat(c,:) = spval
        ch4_surf_aere_sat(c) = spval
        ch4_surf_aere_unsat(c) = spval
        ch4_ebul_depth_sat(c,:) = spval
        ch4_ebul_depth_unsat(c,:) = spval
        ch4_ebul_total_sat(c) = spval
        ch4_ebul_total_unsat(c) = spval
        ch4_surf_ebul_sat(c) = spval
        ch4_surf_ebul_unsat(c) = spval
        ch4_surf_ebul_lake(c) = spval
        ch4_surf_diff_sat(c) = spval
        ch4_surf_diff_unsat(c) = spval
        ch4_surf_diff_lake(c) = spval
        ch4_dfsat_flux(c) = spval
        zwt_ch4_unsat(c) = spval
        grnd_ch4_cond(c) = spval
        conc_ch4_lake(c,:) = spval
        conc_o2_lake(c,:) = spval
        finundated(c) = spval
        fphr(c,:) = spval
        sif(c)  = spval
        rootfr(c,:) = spval
        o2stress_unsat(c,:) = spval
        o2stress_sat(c,:) = spval
        ch4stress_unsat(c,:) = spval
        ch4stress_sat(c,:) = spval
        totcolch4(c) = 0._r8  ! Set to zero for inactive columns so that this can be used
                              ! as an appropriate area-weighted gridcell average soil methane
                              ! content.


     end if

     if (ltype(l) == istdlak .and. .not. allowlakeprod) then
     ! To avoid rare pathologies with allowlakeprod switching between restarts
        conc_ch4_sat(c,:) = 0._r8
        conc_o2_sat(c,:) = 0._r8
        do j=1,nlevsoi
           lake_soilc(c,j) = 580._r8*cellorg(c,j)
        end do
     end if

  end do

  if (masterproc) write (iulog,*) 'Successfully initialized non-state variables for CH4 Mod'

    end associate 
 end subroutine initTimeConst_ch4

#endif
! Defined LCH4
end module initch4Mod
