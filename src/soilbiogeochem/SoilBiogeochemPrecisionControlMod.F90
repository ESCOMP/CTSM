module SoilBiogeochemPrecisionControlMod

  !----------------------------------------------------------------------- 
  ! !DESCRIPTION:
  ! controls on very low values in critical state variables 
  ! 
  ! !USES:
  use shr_kind_mod                    , only : r8 => shr_kind_r8
  use clm_varpar                      , only : ndecomp_pools
  use SoilBiogeochemCarbonStateType   , only : soilbiogeochem_carbonstate_type
  use SoilBiogeochemNitrogenStateType , only : soilbiogeochem_nitrogenstate_type
  use ColumnType                      , only : col
  !
  implicit none
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public:: SoilBiogeochemPrecisionControlInit ! Initialization
  public:: SoilBiogeochemPrecisionControl     ! Apply precision control to soil biogeochemistry carbon and nitrogen states

  ! !PUBLIC DATA:
  real(r8), public :: ccrit                   ! critical carbon state value for truncation (gC/m2)
  real(r8), public :: ncrit                   ! critical nitrogen state value for truncation (gN/m2)
  !----------------------------------------------------------------------- 

contains

  !-----------------------------------------------------------------------
  subroutine SoilBiogeochemPrecisionControlInit( soilbiogeochem_carbonstate_inst, c13_soilbiogeochem_carbonstate_inst, &
       c14_soilbiogeochem_carbonstate_inst, soilbiogeochem_nitrogenstate_inst)

    !
    ! !DESCRIPTION: 
    ! Initialization of soil biogeochemistry precision control
    !
    ! !USES:
    use clm_varctl , only : use_c13, use_c14
    !
    ! !ARGUMENTS:
    type(soilbiogeochem_carbonstate_type)   , intent(inout) :: soilbiogeochem_carbonstate_inst
    type(soilbiogeochem_carbonstate_type)   , intent(inout) :: c13_soilbiogeochem_carbonstate_inst
    type(soilbiogeochem_carbonstate_type)   , intent(inout) :: c14_soilbiogeochem_carbonstate_inst
    type(soilbiogeochem_nitrogenstate_type) , intent(inout) :: soilbiogeochem_nitrogenstate_inst
    !
    ! !LOCAL VARIABLES:
    real(r8), parameter :: totvegcthresh = 1.0_r8   ! Total vegetation carbon threshold to zero out decomposition pools
    !-----------------------------------------------------------------------
    ccrit    =  1.e-8_r8              ! critical carbon state value for truncation (gC/m2)
    ncrit    =  1.e-8_r8              ! critical nitrogen state value for truncation (gN/m2)

    call soilbiogeochem_carbonstate_inst%setTotVgCThresh( totvegcthresh )
    if ( use_c13 )then
        call c13_soilbiogeochem_carbonstate_inst%setTotVgCThresh( totvegcthresh )
    end if
    if ( use_c14 )then
        call c14_soilbiogeochem_carbonstate_inst%setTotVgCThresh( totvegcthresh )
    end if
    call soilbiogeochem_nitrogenstate_inst%setTotVgCThresh( totvegcthresh )

  end subroutine SoilBiogeochemPrecisionControlInit

  !-----------------------------------------------------------------------
  subroutine SoilBiogeochemPrecisionControl(num_bgc_soilc, filter_bgc_soilc, &
       soilbiogeochem_carbonstate_inst, c13_soilbiogeochem_carbonstate_inst, &
       c14_soilbiogeochem_carbonstate_inst, soilbiogeochem_nitrogenstate_inst)

    !
    ! !DESCRIPTION: 
    ! On the radiation time step, force leaf and deadstem c and n to 0 if
    ! they get too small.
    !
    ! !USES:
    use clm_varctl , only : iulog, use_c13, use_c14, use_nitrif_denitrif
    use clm_varpar , only : nlevdecomp
    use CNSharedParamsMod, only: use_fun
    !
    ! !ARGUMENTS:
    integer                                 , intent(in)    :: num_bgc_soilc       ! number of bgc soil columns in filter
    integer                                 , intent(in)    :: filter_bgc_soilc(:) ! filter for bgc soil columns
    type(soilbiogeochem_carbonstate_type)   , intent(inout) :: soilbiogeochem_carbonstate_inst
    type(soilbiogeochem_carbonstate_type)   , intent(inout) :: c13_soilbiogeochem_carbonstate_inst
    type(soilbiogeochem_carbonstate_type)   , intent(inout) :: c14_soilbiogeochem_carbonstate_inst
    type(soilbiogeochem_nitrogenstate_type) , intent(inout) :: soilbiogeochem_nitrogenstate_inst
    !
    ! !LOCAL VARIABLES:
    integer :: c,j,k ! indices
    integer :: fc    ! filter indices
    real(r8):: cc,cn ! truncation terms for column-level corrections
    real(r8):: cc13  ! truncation terms for column-level corrections
    real(r8):: cc14  ! truncation terms for column-level corrections
    !-----------------------------------------------------------------------

    ! soilbiogeochem_carbonstate_inst%ctrunc_vr_col          Output:  [real(r8) (:,:)   ]  (gC/m3) column-level sink for C truncation      
    ! soilbiogeochem_carbonstate_inst%decomp_cpools_vr_col   Output:  [real(r8) (:,:,:) ]  (gC/m3)  vertically-resolved decomposing (litter, cwd, soil) c pools

    ! soilbiogeochem_nitrogenstate_inst%ntrunc_vr_col        Output:  [real(r8) (:,:)   ]  (gN/m3) column-level sink for N truncation      
    ! soilbiogeochem_nitrogenstate_inst%decomp_npools_vr_col Output:  [real(r8) (:,:,:) ]  (gC/m3)  vertically-resolved decomposing (litter, cwd, soil) N pools
    ! soilbiogeochem_nitrogenstate_inst%smin_nh4_vr_col      Output:  [real(r8) (:,:)   ]  (gN/m3) soil mineral NH4                        
    ! soilbiogeochem_nitrogenstate_inst%smin_no3_vr_col      Output:  [real(r8) (:,:)   ]  (gN/m3) soil mineral NO3                        

    associate(&
         cs    => soilbiogeochem_carbonstate_inst     , &
         ns    => soilbiogeochem_nitrogenstate_inst   , & 
         c13cs => c13_soilbiogeochem_carbonstate_inst , &
         c14cs => c14_soilbiogeochem_carbonstate_inst   &
         )

      ! column loop
      do fc = 1,num_bgc_soilc
         c = filter_bgc_soilc(fc)

         do j = 1,nlevdecomp
            ! initialize the column-level C and N truncation terms
            cc = 0._r8
            if ( use_c13 ) cc13 = 0._r8
            if ( use_c14 ) cc14 = 0._r8
            cn = 0._r8

            ! do tests on state variables for precision control
            ! for linked C-N state variables, perform precision test on
            ! the C component, but truncate both C and N components


            ! all decomposing pools C and N
            do k = 1, ndecomp_pools

               if (abs(cs%decomp_cpools_vr_col(c,j,k)) < ccrit) then

                  cc = cc + cs%decomp_cpools_vr_col(c,j,k)
                  cs%decomp_cpools_vr_col(c,j,k) = 0._r8

                  cn = cn + ns%decomp_npools_vr_col(c,j,k)
                  ns%decomp_npools_vr_col(c,j,k) = 0._r8

                  if ( use_c13 ) then
                     cc13 = cc13 + c13cs%decomp_cpools_vr_col(c,j,k)
                     c13cs%decomp_cpools_vr_col(c,j,k) = 0._r8
                  endif
                  if ( use_c14 ) then
                     cc14 = cc14 + c14cs%decomp_cpools_vr_col(c,j,k)
                     c14cs%decomp_cpools_vr_col(c,j,k) = 0._r8
                  endif
               end if

            end do

            ! not doing precision control on soil mineral N, since it will
            ! be getting the N truncation flux anyway.

            cs%ctrunc_vr_col(c,j) = cs%ctrunc_vr_col(c,j) + cc

            ns%ntrunc_vr_col(c,j) = ns%ntrunc_vr_col(c,j) + cn

            if ( use_c13 ) then
               c13cs%ctrunc_vr_col(c,j) = c13cs%ctrunc_vr_col(c,j) + cc13
            endif
            if ( use_c14 ) then
               c14cs%ctrunc_vr_col(c,j) = c14cs%ctrunc_vr_col(c,j) + cc14
            endif
         end do

      end do   ! end of column loop

     if(.not.use_fun)then
      if (use_nitrif_denitrif) then
         ! remove small negative perturbations for stability purposes, if any should arise.
        
         do fc = 1,num_bgc_soilc
            c = filter_bgc_soilc(fc)
            do j = 1,nlevdecomp
               if (abs(ns%smin_no3_vr_col(c,j)) < ncrit/1e4_r8) then
                  if ( ns%smin_no3_vr_col(c,j)  < 0._r8 ) then
                     !write(iulog, *) '-10^-12 < smin_no3 < 0. resetting to zero.'
                     !write(iulog, *) 'smin_no3_vr_col(c,j), c, j: ', ns%smin_no3_vr_col(c,j), c, j
                     ns%smin_no3_vr_col(c,j) = 0._r8
                  endif
               end if
               if (abs(ns%smin_nh4_vr_col(c,j)) < ncrit/1e4_r8) then
                  if ( ns%smin_nh4_vr_col(c,j)  < 0._r8 ) then
                     !write(iulog, *) '-10^-12 < smin_nh4 < 0. resetting to zero.'
                     !write(iulog, *) 'smin_nh4_vr_col(c,j), c, j: ', ns%smin_nh4_vr_col(c,j), c, j
                     ns%smin_nh4_vr_col(c,j) = 0._r8
                  endif
               end if
            end do
         end do
      endif
     endif

    end associate

  end subroutine SoilBiogeochemPrecisionControl

end module SoilBiogeochemPrecisionControlMod
