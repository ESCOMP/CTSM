module CNPrecisionControlMod

#include "shr_assert.h"

  !----------------------------------------------------------------------- 
  ! !DESCRIPTION:
  ! controls on very low values in critical state variables 
  ! 
  ! !USES:
  use shr_kind_mod           , only : r8 => shr_kind_r8
  use CNVegCarbonStateType   , only : cnveg_carbonstate_type
  use CNVegNitrogenStateType , only : cnveg_nitrogenstate_type
  use PatchType              , only : patch
  use abortutils             , only : endrun
  !
  implicit none
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public:: CNPrecisionControlReadNML
  public:: CNPrecisionControl

  ! !PUBLIC DATA:
  real(r8), public :: ccrit    =  1.e-8_r8              ! critical carbon state value for truncation (gC/m2)
  real(r8), public :: cnegcrit = -6.e+1_r8              ! critical negative carbon state value for abort (gC/m2)
  real(r8), public :: ncrit    =  1.e-8_r8              ! critical nitrogen state value for truncation (gN/m2)
  real(r8), public :: nnegcrit = -7.e+0_r8              ! critical negative nitrogen state value for abort (gN/m2)
  real(r8), public, parameter :: n_min = 0.000000001_r8 ! Minimum Nitrogen value to use when calculating CN ratio (gN/m2)

  ! !PRIVATE DATA:
  logical, private :: prec_control_for_froot = .true.   ! If true do precision control for frootc/frootn
  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  !----------------------------------------------------------------------- 

contains

  !-----------------------------------------------------------------------
  subroutine CNPrecisionControlReadNML( NLFilename )
    !
    ! !DESCRIPTION:
    ! Read the namelist for CN Precision control
    !
    ! !USES:
    use fileutils      , only : getavu, relavu, opnfil
    use shr_nl_mod     , only : shr_nl_find_group_name
    use spmdMod        , only : masterproc, mpicom
    use shr_mpi_mod    , only : shr_mpi_bcast
    use clm_varctl     , only : iulog, use_nguardrail
    use shr_log_mod    , only : errMsg => shr_log_errMsg
    !
    ! !ARGUMENTS:
    character(len=*), intent(in) :: NLFilename ! Namelist filename
    !
    ! !LOCAL VARIABLES:
    integer :: ierr                 ! error code
    integer :: unitn                ! unit for namelist file

    character(len=*), parameter :: subname = 'CNPrecisionControlReadNML'
    character(len=*), parameter :: nmlname = 'cnprecision_inparm'
    !-----------------------------------------------------------------------
    namelist /cnprecision_inparm/ ncrit, ccrit, cnegcrit, nnegcrit

    if (masterproc) then
       unitn = getavu()
       write(iulog,*) 'Read in '//nmlname//'  namelist'
       call opnfil (NLFilename, unitn, 'F')
       call shr_nl_find_group_name(unitn, nmlname, status=ierr)
       if (ierr == 0) then
          read(unitn, nml=cnprecision_inparm, iostat=ierr)
          if (ierr /= 0) then
             call endrun(msg="ERROR reading "//nmlname//"namelist"//errmsg(sourcefile, __LINE__))
          end if
       else
          call endrun(msg="ERROR could NOT find "//nmlname//"namelist"//errmsg(sourcefile, __LINE__))
       end if
       call relavu( unitn )
    end if

    call shr_mpi_bcast (ncrit   , mpicom)
    call shr_mpi_bcast (ccrit   , mpicom)
    call shr_mpi_bcast (nnegcrit, mpicom)
    call shr_mpi_bcast (cnegcrit, mpicom)

    if (masterproc) then
       write(iulog,*) ' '
       write(iulog,*) nmlname//' settings:'
       write(iulog,nml=cnprecision_inparm)
       write(iulog,*) ' '
    end if

    ! Have precision control for froot be determined by use_nguardrail setting
    prec_control_for_froot = .not. use_nguardrail

  end subroutine CNPrecisionControlReadNML

  !-----------------------------------------------------------------------
  subroutine CNPrecisionControl(bounds, num_soilp, filter_soilp, &
       cnveg_carbonstate_inst, c13_cnveg_carbonstate_inst, c14_cnveg_carbonstate_inst, &
       cnveg_nitrogenstate_inst)
    !
    ! !DESCRIPTION: 
    ! Force leaf and deadstem c and n to 0 if they get too small.
    !
    ! !USES:
    use clm_varctl , only : iulog, use_c13, use_c14
    use clm_varctl , only : use_crop
    use pftconMod  , only : nc3crop
    use decompMod  , only : bounds_type
    !
    ! !ARGUMENTS:
    type(bounds_type)              , intent(in)    :: bounds          ! bounds
    integer                        , intent(in)    :: num_soilp       ! number of soil patchs in filter
    integer                        , intent(in)    :: filter_soilp(:) ! filter for soil patches
    type(cnveg_carbonstate_type)   , intent(inout) :: cnveg_carbonstate_inst
    type(cnveg_carbonstate_type)   , intent(inout) :: c13_cnveg_carbonstate_inst
    type(cnveg_carbonstate_type)   , intent(inout) :: c14_cnveg_carbonstate_inst
    type(cnveg_nitrogenstate_type) , intent(inout) :: cnveg_nitrogenstate_inst
    !
    ! !LOCAL VARIABLES:
    integer :: p,j,k                             ! indices
    integer :: fp                                ! filter indices
    real(r8):: pc(bounds%begp:bounds%endp)       ! truncation terms for patch-level corrections Carbon
    real(r8):: pn(bounds%begp:bounds%endp)       ! truncation terms for patch-level corrections nitrogen
    real(r8):: pc13(bounds%begp:bounds%endp)     ! truncation terms for patch-level corrections
    real(r8):: pc14(bounds%begp:bounds%endp)     ! truncation terms for patch-level corrections
    !-----------------------------------------------------------------------

    ! cnveg_carbonstate_inst%cpool_patch                     Output:  [real(r8) (:)     ]  (gC/m2) temporary photosynthate C pool            
    ! cnveg_carbonstate_inst%deadcrootc_patch                Output:  [real(r8) (:)     ]  (gC/m2) dead coarse root C                        
    ! cnveg_carbonstate_inst%deadcrootc_storage_patch        Output:  [real(r8) (:)     ]  (gC/m2) dead coarse root C storage                
    ! cnveg_carbonstate_inst%deadcrootc_xfer_patch           Output:  [real(r8) (:)     ]  (gC/m2) dead coarse root C transfer               
    ! cnveg_carbonstate_inst%deadstemc_patch                 Output:  [real(r8) (:)     ]  (gC/m2) dead stem C                               
    ! cnveg_carbonstate_inst%deadstemc_storage_patch         Output:  [real(r8) (:)     ]  (gC/m2) dead stem C storage                       
    ! cnveg_carbonstate_inst%deadstemc_xfer_patch            Output:  [real(r8) (:)     ]  (gC/m2) dead stem C transfer                      
    ! cnveg_carbonstate_inst%frootc_patch                    Output:  [real(r8) (:)     ]  (gC/m2) fine root C                               
    ! cnveg_carbonstate_inst%frootc_storage_patch            Output:  [real(r8) (:)     ]  (gC/m2) fine root C storage                       
    ! cnveg_carbonstate_inst%frootc_xfer_patch               Output:  [real(r8) (:)     ]  (gC/m2) fine root C transfer                      
    ! cnveg_carbonstate_inst%gresp_storage_patch             Output:  [real(r8) (:)     ]  (gC/m2) growth respiration storage                
    ! cnveg_carbonstate_inst%gresp_xfer_patch                Output:  [real(r8) (:)     ]  (gC/m2) growth respiration transfer               
    ! cnveg_carbonstate_inst%leafc_patch                     Output:  [real(r8) (:)     ]  (gC/m2) leaf C                                    
    ! cnveg_carbonstate_inst%leafc_storage_patch             Output:  [real(r8) (:)     ]  (gC/m2) leaf C storage                            
    ! cnveg_carbonstate_inst%leafc_xfer_patch                Output:  [real(r8) (:)     ]  (gC/m2) leaf C transfer                           
    ! cnveg_carbonstate_inst%livecrootc_patch                Output:  [real(r8) (:)     ]  (gC/m2) live coarse root C                        
    ! cnveg_carbonstate_inst%livecrootc_storage_patch        Output:  [real(r8) (:)     ]  (gC/m2) live coarse root C storage                
    ! cnveg_carbonstate_inst%livecrootc_xfer_patch           Output:  [real(r8) (:)     ]  (gC/m2) live coarse root C transfer               
    ! cnveg_carbonstate_inst%livestemc_patch                 Output:  [real(r8) (:)     ]  (gC/m2) live stem C                               
    ! cnveg_carbonstate_inst%livestemc_storage_patch         Output:  [real(r8) (:)     ]  (gC/m2) live stem C storage                       
    ! cnveg_carbonstate_inst%livestemc_xfer_patch            Output:  [real(r8) (:)     ]  (gC/m2) live stem C transfer                      
    ! cnveg_carbonstate_inst%ctrunc_patch                    Output:  [real(r8) (:)     ]  (gC/m2) patch-level sink for C truncation           
    ! cnveg_carbonstate_inst%xsmrpool_patch                  Output:  [real(r8) (:)     ]  (gC/m2) execss maint resp C pool                  
    ! cnveg_carbonstate_inst%grainc_patch                    Output:  [real(r8) (:)     ]  (gC/m2) grain C                                   
    ! cnveg_carbonstate_inst%grainc_storage_patch            Output:  [real(r8) (:)     ]  (gC/m2) grain C storage                           
    ! cnveg_carbonstate_inst%grainc_xfer_patch               Output:  [real(r8) (:)     ]  (gC/m2) grain C transfer                          
    
    ! cnveg_nitrogenstate_inst%deadcrootn_patch              Output:  [real(r8) (:)     ]  (gN/m2) dead coarse root N                        
    ! cnveg_nitrogenstate_inst%deadcrootn_storage_patch      Output:  [real(r8) (:)     ]  (gN/m2) dead coarse root N storage                
    ! cnveg_nitrogenstate_inst%deadcrootn_xfer_patch         Output:  [real(r8) (:)     ]  (gN/m2) dead coarse root N transfer               
    ! cnveg_nitrogenstate_inst%deadstemn_patch               Output:  [real(r8) (:)     ]  (gN/m2) dead stem N                               
    ! cnveg_nitrogenstate_inst%deadstemn_storage_patch       Output:  [real(r8) (:)     ]  (gN/m2) dead stem N storage                       
    ! cnveg_nitrogenstate_inst%deadstemn_xfer_patch          Output:  [real(r8) (:)     ]  (gN/m2) dead stem N transfer                      
    ! cnveg_nitrogenstate_inst%frootn_patch                  Output:  [real(r8) (:)     ]  (gN/m2) fine root N                               
    ! cnveg_nitrogenstate_inst%frootn_storage_patch          Output:  [real(r8) (:)     ]  (gN/m2) fine root N storage                       
    ! cnveg_nitrogenstate_inst%frootn_xfer_patch             Output:  [real(r8) (:)     ]  (gN/m2) fine root N transfer                      
    ! cnveg_nitrogenstate_inst%leafn_patch                   Output:  [real(r8) (:)     ]  (gN/m2) leaf N                                    
    ! cnveg_nitrogenstate_inst%leafn_storage_patch           Output:  [real(r8) (:)     ]  (gN/m2) leaf N storage                            
    ! cnveg_nitrogenstate_inst%leafn_xfer_patch              Output:  [real(r8) (:)     ]  (gN/m2) leaf N transfer                           
    ! cnveg_nitrogenstate_inst%livecrootn_patch              Output:  [real(r8) (:)     ]  (gN/m2) live coarse root N                        
    ! cnveg_nitrogenstate_inst%livecrootn_storage_patch      Output:  [real(r8) (:)     ]  (gN/m2) live coarse root N storage                
    ! cnveg_nitrogenstate_inst%livecrootn_xfer_patch         Output:  [real(r8) (:)     ]  (gN/m2) live coarse root N transfer               
    ! cnveg_nitrogenstate_inst%grainn_patch                  Output:  [real(r8) (:)     ]  (gC/m2) grain N                                   
    ! cnveg_nitrogenstate_inst%grainn_storage_patch          Output:  [real(r8) (:)     ]  (gC/m2) grain N storage                           
    ! cnveg_nitrogenstate_inst%grainn_xfer_patch             Output:  [real(r8) (:)     ]  (gC/m2) grain N transfer                          
    ! cnveg_nitrogenstate_inst%livestemn_patch               Output:  [real(r8) (:)     ]  (gN/m2) live stem N                               
    ! cnveg_nitrogenstate_inst%livestemn_storage_patch       Output:  [real(r8) (:)     ]  (gN/m2) live stem N storage                       
    ! cnveg_nitrogenstate_inst%livestemn_xfer_patch          Output:  [real(r8) (:)     ]  (gN/m2) live stem N transfer                      
    ! cnveg_nitrogenstate_inst%npool_patch                   Output:  [real(r8) (:)     ]  (gN/m2) temporary plant N pool                    
    ! cnveg_nitrogenstate_inst%ntrunc_patch                  Output:  [real(r8) (:)     ]  (gN/m2) patch-level sink for N truncation           
    ! cnveg_nitrogenstate_inst%retransn_patch                Output:  [real(r8) (:)     ]  (gN/m2) plant pool of retranslocated N            

    
    associate(                                           &
         cs     => cnveg_carbonstate_inst              , &
         ns     => cnveg_nitrogenstate_inst            , &
         c13cs  => c13_cnveg_carbonstate_inst          , &
         c14cs  => c14_cnveg_carbonstate_inst            &
         )

      ! patch loop
      do fp = 1,num_soilp
         p = filter_soilp(fp)

         ! initialize the patch-level C and N truncation terms
         pc(p) = 0._r8
         pn(p) = 0._r8
         if ( use_c13 ) pc13(p) = 0._r8
         if ( use_c14 ) pc14(p) = 0._r8
      end do

      ! do tests on state variables for precision control
      ! for linked C-N state variables, perform precision test on
      ! the C component, but truncate C, C13, and N components

      ! leaf C and N
      call TruncateCandNStates( bounds, filter_soilp, num_soilp, cs%leafc_patch(bounds%begp:bounds%endp), &
                                ns%leafn_patch(bounds%begp:bounds%endp), &
                                pc(bounds%begp:), pn(bounds%begp:), __LINE__, &
                                c13=c13cs%leafc_patch, c14=c14cs%leafc_patch, &
                                pc13=pc13(bounds%begp:), pc14=pc14(bounds%begp:) )

      ! leaf storage C and N
      call TruncateCandNStates( bounds, filter_soilp, num_soilp, cs%leafc_storage_patch(bounds%begp:bounds%endp), &
                                ns%leafn_storage_patch(bounds%begp:bounds%endp), pc(bounds%begp:), pn(bounds%begp:), __LINE__, &
                                c13=c13cs%leafc_storage_patch, c14=c14cs%leafc_storage_patch, &
                                pc13=pc13(bounds%begp:), pc14=pc14(bounds%begp:) )

      ! leaf transfer C and N
      call TruncateCandNStates( bounds, filter_soilp, num_soilp, cs%leafc_xfer_patch(bounds%begp:bounds%endp), &
                                ns%leafn_xfer_patch(bounds%begp:bounds%endp), pc(bounds%begp:), pn(bounds%begp:), __LINE__, &
                                c13=c13cs%leafc_xfer_patch, c14=c14cs%leafc_xfer_patch, &
                                pc13=pc13(bounds%begp:), pc14=pc14(bounds%begp:) )

      ! froot C and N
      ! EBK KO DML: For some reason frootc/frootn can go negative and allowing
      ! it to be negative is important for C4 crops (otherwise they die) Jun/3/2016
      if ( prec_control_for_froot ) then
         call TruncateCandNStates( bounds, filter_soilp, num_soilp, cs%frootc_patch(bounds%begp:bounds%endp),  &
                                   ns%frootn_patch(bounds%begp:bounds%endp), pc(bounds%begp:), pn(bounds%begp:), __LINE__, &
                                   c13=c13cs%frootc_patch, c14=c14cs%frootc_patch,  &
                                   pc13=pc13(bounds%begp:), pc14=pc14(bounds%begp:), allowneg=.true. )
      end if

      ! froot storage C and N
      call TruncateCandNStates( bounds, filter_soilp, num_soilp, cs%frootc_storage_patch(bounds%begp:bounds%endp), &
                      ns%frootn_storage_patch(bounds%begp:bounds%endp), pc(bounds%begp:), pn(bounds%begp:), &
                      __LINE__, c13=c13cs%frootc_storage_patch, c14=c14cs%frootc_storage_patch, &
                      pc13=pc13(bounds%begp:), pc14=pc14(bounds%begp:) )

      ! froot transfer C and N
      call TruncateCandNStates( bounds, filter_soilp, num_soilp, cs%frootc_xfer_patch(bounds%begp:bounds%endp), &
                                ns%frootn_xfer_patch(bounds%begp:bounds%endp), pc(bounds%begp:), pn(bounds%begp:), __LINE__, &
                                c13=c13cs%frootc_xfer_patch, c14=c14cs%frootc_xfer_patch, &
                                pc13=pc13(bounds%begp:), pc14=pc14(bounds%begp:) )

      if ( use_crop )then
         ! grain C and N
         call TruncateCandNStates( bounds, filter_soilp, num_soilp, cs%grainc_patch(bounds%begp:bounds%endp), &
                                   ns%grainn_patch(bounds%begp:bounds%endp), pc(bounds%begp:), pn(bounds%begp:), __LINE__, &
                                   c13=c13cs%grainc_patch, c14=c14cs%grainc_patch, &
                                   pc13=pc13(bounds%begp:), pc14=pc14(bounds%begp:), croponly=.true. )

         ! grain storage C and N
         call TruncateCandNStates( bounds, filter_soilp, num_soilp, cs%grainc_storage_patch(bounds%begp:bounds%endp), &
                                   ns%grainn_storage_patch(bounds%begp:bounds%endp), pc(bounds%begp:), pn(bounds%begp:), &
                                   __LINE__, c13=c13cs%grainc_storage_patch, c14=c14cs%grainc_storage_patch, &
                                   pc13=pc13(bounds%begp:), pc14=pc14(bounds%begp:), croponly=.true. )

         ! grain transfer C and N
         call TruncateCandNStates( bounds, filter_soilp, num_soilp, cs%grainc_xfer_patch(bounds%begp:bounds%endp), &
                                   ns%grainn_xfer_patch(bounds%begp:bounds%endp), pc(bounds%begp:), pn(bounds%begp:), __LINE__, &
                                   c13=c13cs%grainc_xfer_patch, c14=c14cs%grainc_xfer_patch, &
                                   pc13=pc13(bounds%begp:), pc14=pc14(bounds%begp:), croponly=.true. )

         ! grain transfer C and N
         call TruncateCandNStates( bounds, filter_soilp, num_soilp, cs%cropseedc_deficit_patch(bounds%begp:bounds%endp), &
                                   ns%cropseedn_deficit_patch(bounds%begp:bounds%endp), pc(bounds%begp:), &
                                   pn(bounds%begp:), __LINE__, &
                                   c13=c13cs%cropseedc_deficit_patch, c14=c14cs%cropseedc_deficit_patch, &
                                   pc13=pc13(bounds%begp:), pc14=pc14(bounds%begp:), allowneg=.true., croponly=.true. )

      end if

      ! livestem C and N
      call TruncateCandNStates( bounds, filter_soilp, num_soilp, cs%livestemc_patch(bounds%begp:bounds%endp), &
                                ns%livestemn_patch(bounds%begp:bounds%endp), pc(bounds%begp:), pn(bounds%begp:), __LINE__, &
                                c13=c13cs%livestemc_patch, c14=c14cs%livestemc_patch, &
                                pc13=pc13(bounds%begp:), pc14=pc14(bounds%begp:) )

      ! livestem storage C and N
      call TruncateCandNStates( bounds, filter_soilp, num_soilp, cs%livestemc_storage_patch(bounds%begp:bounds%endp), &
                      ns%livestemn_storage_patch(bounds%begp:bounds%endp), pc(bounds%begp:), pn(bounds%begp:), &
                      __LINE__, c13=c13cs%livestemc_storage_patch, c14=c14cs%livestemc_storage_patch, &
                      pc13=pc13(bounds%begp:), pc14=pc14(bounds%begp:) )

      ! livestem transfer C and N
      call TruncateCandNStates( bounds, filter_soilp, num_soilp, cs%livestemc_xfer_patch(bounds%begp:bounds%endp), &
                      ns%livestemn_xfer_patch(bounds%begp:bounds%endp), pc(bounds%begp:), pn(bounds%begp:), &
                      __LINE__, c13=c13cs%livestemc_xfer_patch, c14=c14cs%livestemc_xfer_patch, &
                      pc13=pc13(bounds%begp:), pc14=pc14(bounds%begp:) )

      ! deadstem C and N
      call TruncateCandNStates( bounds, filter_soilp, num_soilp, cs%deadstemc_patch(bounds%begp:bounds%endp), &
                                ns%deadstemn_patch(bounds%begp:bounds%endp), pc(bounds%begp:), pn(bounds%begp:), __LINE__, &
                                c13=c13cs%deadstemc_patch, c14=c14cs%deadstemc_patch, &
                                pc13=pc13(bounds%begp:), pc14=pc14(bounds%begp:) )
      ! deadstem storage C and N
      call TruncateCandNStates( bounds, filter_soilp, num_soilp, cs%deadstemc_storage_patch(bounds%begp:bounds%endp), &
                      ns%deadstemn_storage_patch(bounds%begp:bounds%endp), pc(bounds%begp:), pn(bounds%begp:), &
                      __LINE__, c13=c13cs%deadstemc_storage_patch, c14=c14cs%deadstemc_storage_patch, &
                      pc13=pc13(bounds%begp:), pc14=pc14(bounds%begp:) )

      ! deadstem transfer C and N
      call TruncateCandNStates( bounds, filter_soilp, num_soilp, cs%deadstemc_xfer_patch(bounds%begp:bounds%endp), &
                      ns%deadstemn_xfer_patch(bounds%begp:bounds%endp), pc(bounds%begp:), pn(bounds%begp:), &
                      __LINE__, c13=c13cs%deadstemc_xfer_patch, c14=c14cs%deadstemc_xfer_patch, &
                      pc13=pc13(bounds%begp:), pc14=pc14(bounds%begp:) )

      ! livecroot C and N
      call TruncateCandNStates( bounds, filter_soilp, num_soilp, cs%livecrootc_patch(bounds%begp:bounds%endp), &
                                ns%livecrootn_patch(bounds%begp:bounds%endp), pc(bounds%begp:), pn(bounds%begp:), __LINE__, &
                                c13=c13cs%livecrootc_patch, c14=c14cs%livecrootc_patch, &
                                pc13=pc13(bounds%begp:), pc14=pc14(bounds%begp:) )

      ! livecroot storage C and N
      call TruncateCandNStates( bounds, filter_soilp, num_soilp, cs%livecrootc_storage_patch(bounds%begp:bounds%endp), &
                      ns%livecrootn_storage_patch(bounds%begp:bounds%endp), pc(bounds%begp:), pn(bounds%begp:), &
                      __LINE__, c13=c13cs%livecrootc_storage_patch, c14=c14cs%livecrootc_storage_patch, &
                      pc13=pc13(bounds%begp:), pc14=pc14(bounds%begp:) )

      ! livecroot transfer C and N
      call TruncateCandNStates( bounds, filter_soilp, num_soilp, cs%livecrootc_xfer_patch(bounds%begp:bounds%endp), &
                      ns%livecrootn_xfer_patch(bounds%begp:bounds%endp), pc(bounds%begp:), pn(bounds%begp:), &
                      __LINE__, c13=c13cs%livecrootc_xfer_patch, c14=c14cs%livecrootc_xfer_patch, &
                      pc13=pc13(bounds%begp:), pc14=pc14(bounds%begp:) )

      ! deadcroot C and N
      call TruncateCandNStates( bounds, filter_soilp, num_soilp, cs%deadcrootc_patch(bounds%begp:bounds%endp), &
                                ns%deadcrootn_patch(bounds%begp:bounds%endp), pc(bounds%begp:), pn(bounds%begp:), __LINE__, &
                                c13=c13cs%deadcrootc_patch, c14=c14cs%deadcrootc_patch, &
                                pc13=pc13(bounds%begp:), pc14=pc14(bounds%begp:) )

      ! deadcroot storage C and N
      call TruncateCandNStates( bounds, filter_soilp, num_soilp, cs%deadcrootc_storage_patch(bounds%begp:bounds%endp), &
                      ns%deadcrootn_storage_patch(bounds%begp:bounds%endp), pc(bounds%begp:), pn(bounds%begp:), &
                      __LINE__, c13=c13cs%deadcrootc_storage_patch, c14=c14cs%deadcrootc_storage_patch, &
                      pc13=pc13(bounds%begp:), pc14=pc14(bounds%begp:) )

      ! deadcroot transfer C and N
      call TruncateCandNStates( bounds, filter_soilp, num_soilp, cs%deadcrootc_xfer_patch(bounds%begp:bounds%endp), &
                      ns%deadcrootn_xfer_patch(bounds%begp:bounds%endp), pc(bounds%begp:), pn(bounds%begp:), &
                      __LINE__, c13=c13cs%deadcrootc_xfer_patch, c14=c14cs%deadcrootc_xfer_patch, &
                      pc13=pc13(bounds%begp:), pc14=pc14(bounds%begp:) )

      ! gresp_storage (C only)
      call TruncateCStates( bounds, filter_soilp, num_soilp, cs%gresp_storage_patch(bounds%begp:bounds%endp), &
                            pc(bounds%begp:), __LINE__, &
                            c13=c13cs%gresp_storage_patch, c14=c14cs%gresp_storage_patch, &
                            pc13=pc13(bounds%begp:), pc14=pc14(bounds%begp:) )

      ! gresp_xfer(c only)
      call TruncateCStates( bounds, filter_soilp, num_soilp, cs%gresp_xfer_patch(bounds%begp:bounds%endp), &
                            pc(bounds%begp:), __LINE__, &
                            c13=c13cs%gresp_xfer_patch, c14=c14cs%gresp_xfer_patch, &
                            pc13=pc13(bounds%begp:), pc14=pc14(bounds%begp:) )

      ! cpool (C only)
      call TruncateCStates( bounds, filter_soilp, num_soilp, cs%cpool_patch(bounds%begp:bounds%endp), &
                            pc(bounds%begp:), __LINE__, &
                            c13=c13cs%cpool_patch, c14=c14cs%cpool_patch, &
                            pc13=pc13(bounds%begp:), pc14=pc14(bounds%begp:) )

      if ( use_crop )then
         ! xsmrpool (C only)
         ! xsmr is a pool to balance the budget and as such can be freely negative
         call TruncateCStates( bounds, filter_soilp, num_soilp, cs%xsmrpool_patch(bounds%begp:bounds%endp), &
                               pc(bounds%begp:), __LINE__, &
                               c13=c13cs%xsmrpool_patch, c14=c14cs%xsmrpool_patch, &
                               pc13=pc13(bounds%begp:), pc14=pc14(bounds%begp:), allowneg=.true., croponly=.true. )

      end if

      ! retransn (N only)
      call TruncateNStates( bounds, filter_soilp, num_soilp, ns%retransn_patch(bounds%begp:bounds%endp), pn(bounds%begp:), &
                            __LINE__ )

      ! npool (N only)
      call TruncateNStates( bounds, filter_soilp, num_soilp, ns%npool_patch(bounds%begp:bounds%endp), pn(bounds%begp:), &
                            __LINE__ )

      ! patch loop
      do fp = 1,num_soilp
         p = filter_soilp(fp)

         cs%ctrunc_patch(p) = cs%ctrunc_patch(p) + pc(p)

         ns%ntrunc_patch(p) = ns%ntrunc_patch(p) + pn(p)

         if ( use_c13 ) then
            c13cs%ctrunc_patch(p) = c13cs%ctrunc_patch(p) + pc13(p)
         endif
         if ( use_c14 ) then
            c14cs%ctrunc_patch(p) = c14cs%ctrunc_patch(p) + pc14(p)
         endif
       end do

    end associate

 end subroutine CNPrecisionControl

 subroutine TruncateCandNStates( bounds, filter_soilp, num_soilp, carbon_patch, nitrogen_patch, pc, pn, lineno, c13, c14, &
                                 pc13, pc14, croponly, allowneg )
    !
    ! !DESCRIPTION:
    ! Truncate paired Carbon and Nitrogen states. If a paired carbon and nitrogen state iare too small truncate 
    ! the pair of them to zero.
    !
    ! !USES:
    use shr_log_mod, only : errMsg => shr_log_errMsg
    use clm_varctl , only : use_c13, use_c14, use_nguardrail
    use clm_varctl , only : iulog
    use pftconMod  , only : nc3crop
    use decompMod  , only : bounds_type
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type)              , intent(in)    :: bounds          ! bounds
    integer                        , intent(in)    :: num_soilp       ! number of soil patchs in filter
    integer                        , intent(in)    :: filter_soilp(:) ! filter for soil patches
    real(r8), intent(inout) :: carbon_patch(bounds%begp:)
    real(r8), intent(inout) :: nitrogen_patch(bounds%begp:)
    real(r8), intent(inout) :: pc(bounds%begp:)
    real(r8), intent(inout) :: pn(bounds%begp:)
    integer,  intent(in)    :: lineno
    real(r8), intent(inout), optional, pointer :: c13(:)
    real(r8), intent(inout), optional, pointer :: c14(:)
    real(r8), intent(inout), optional :: pc13(bounds%begp:)
    real(r8), intent(inout), optional :: pc14(bounds%begp:)
    logical , intent(in)   , optional :: croponly
    logical , intent(in)   , optional :: allowneg

    logical :: lcroponly, lallowneg
    integer :: fp, p

    SHR_ASSERT_ALL((ubound(carbon_patch)   == (/bounds%endp/)), 'ubnd(carb)'//errMsg(sourcefile, lineno))
    SHR_ASSERT_ALL((ubound(nitrogen_patch) == (/bounds%endp/)), 'ubnd(nitro)'//errMsg(sourcefile, lineno))
    SHR_ASSERT_ALL((ubound(pc)             == (/bounds%endp/)), 'ubnd(pc)'//errMsg(sourcefile, lineno))
    SHR_ASSERT_ALL((ubound(pn)             == (/bounds%endp/)), 'ubnd(pn)'//errMsg(sourcefile, lineno))
#ifndef _OPENMP
    if ( present(c13) .and. use_c13 )then
       SHR_ASSERT_ALL((lbound(c13)         == (/bounds%begp/)), 'lbnd(c13)'//errMsg(sourcefile, lineno))
       SHR_ASSERT_ALL((ubound(c13)         == (/bounds%endp/)), 'ubnd(c13)'//errMsg(sourcefile, lineno))
    end if
    if ( present(c14) .and. use_c14 )then
       SHR_ASSERT_ALL((lbound(c14)         == (/bounds%begp/)), 'lbnd(c14)'//errMsg(sourcefile, lineno))
       SHR_ASSERT_ALL((ubound(c14)         == (/bounds%endp/)), 'ubnd(c14)'//errMsg(sourcefile, lineno))
    end if
#endif
    if ( present(pc13) )then
       SHR_ASSERT_ALL((ubound(pc13)        == (/bounds%endp/)), 'ubnd(pc13)'//errMsg(sourcefile, lineno))
    end if
    if ( present(pc14) )then
       SHR_ASSERT_ALL((ubound(pc14)        == (/bounds%endp/)), 'ubnd(pc14)'//errMsg(sourcefile, lineno))
    end if
    ! patch loop
    lcroponly = .false.
    if ( present(croponly) )then
      if ( croponly ) lcroponly = .true.
    end if
    lallowneg = .false.
    if ( present(allowneg) )then
      if (  allowneg ) lallowneg = .true.
    end if
    do fp = 1,num_soilp
       p = filter_soilp(fp)

       if ( .not. lcroponly .or. (patch%itype(p) >= nc3crop) ) then
          if ( .not. lallowneg .and. ((carbon_patch(p) < cnegcrit) .or. (nitrogen_patch(p) < nnegcrit)) ) then
             write(iulog,*) 'ERROR: Carbon or Nitrogen patch negative = ', carbon_patch(p), nitrogen_patch(p)
             write(iulog,*) 'ERROR: limits = ', cnegcrit, nnegcrit
             call endrun(msg='ERROR: carbon or nitrogen state critically negative '//errMsg(sourcefile, lineno))
          else if ( abs(carbon_patch(p)) < ccrit .or. (use_nguardrail .and. abs(nitrogen_patch(p)) < ncrit) ) then
             pc(p) = pc(p) + carbon_patch(p)
             carbon_patch(p) = 0._r8
      
             pn(p) = pn(p) + nitrogen_patch(p)
             nitrogen_patch(p) = 0._r8
   
             if ( use_c13 .and. present(c13) .and. present(pc13) ) then
                pc13(p) = pc13(p) + c13(p)
                c13(p) = 0._r8
             endif
             if ( use_c14 .and. present(c14) .and. present(pc14)) then
                pc14(p) = pc14(p) + c14(p)
                c14(p) = 0._r8
             endif
          end if
       end if
    end do
 end subroutine TruncateCandNStates

 subroutine TruncateCStates( bounds, filter_soilp, num_soilp, carbon_patch, pc, lineno, c13, c14, pc13, pc14, croponly, allowneg )
    !
    ! !DESCRIPTION:
    ! Truncate Carbon states. If a carbon state is too small truncate it to
    ! zero.
    !
    ! !USES:
    use abortutils , only : endrun
    use clm_varctl , only : iulog
    use shr_log_mod, only : errMsg => shr_log_errMsg
    use clm_varctl , only : use_c13, use_c14
    use pftconMod  , only : nc3crop
    use decompMod  , only : bounds_type
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in)    :: bounds          ! bounds
    integer          , intent(in)    :: num_soilp       ! number of soil patchs in filter
    integer          , intent(in)    :: filter_soilp(:) ! filter for soil patches
    real(r8)         , intent(inout) :: carbon_patch(bounds%begp:)
    real(r8)         , intent(inout) :: pc(bounds%begp:)
    integer          , intent(in)    :: lineno
    real(r8)         , intent(inout), optional, pointer :: c13(:)
    real(r8)         , intent(inout), optional, pointer :: c14(:)
    real(r8)         , intent(inout), optional :: pc13(bounds%begp:)
    real(r8)         , intent(inout), optional :: pc14(bounds%begp:)
    logical          , intent(in)   , optional :: croponly
    logical          , intent(in)   , optional :: allowneg

    logical :: lcroponly, lallowneg
    integer :: fp, p

    SHR_ASSERT_ALL((ubound(carbon_patch)   == (/bounds%endp/)), errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(pc)             == (/bounds%endp/)), errMsg(sourcefile, __LINE__))
#ifndef _OPENMP
    if ( present(c13) .and. use_c13 )then
       SHR_ASSERT_ALL((lbound(c13)         == (/bounds%begp/)), errMsg(sourcefile, __LINE__))
       SHR_ASSERT_ALL((ubound(c13)         == (/bounds%endp/)), errMsg(sourcefile, __LINE__))
    end if
    if ( present(c14) .and. use_c14 )then
       SHR_ASSERT_ALL((lbound(c14)         == (/bounds%begp/)), errMsg(sourcefile, __LINE__))
       SHR_ASSERT_ALL((ubound(c14)         == (/bounds%endp/)), errMsg(sourcefile, __LINE__))
    end if
#endif
    if ( present(pc13) )then
       SHR_ASSERT_ALL((ubound(pc13)        == (/bounds%endp/)), errMsg(sourcefile, __LINE__))
    end if
    if ( present(pc14) )then
       SHR_ASSERT_ALL((ubound(pc14)        == (/bounds%endp/)), errMsg(sourcefile, __LINE__))
    end if
    if ( -ccrit < cnegcrit )then
        call endrun(msg='ERROR: cnegcrit should be less than -ccrit: '//errMsg(sourcefile, lineno))
    end if
    lcroponly = .false.
    if ( present(croponly) )then
      if ( croponly ) lcroponly = .true.
    end if
    lallowneg = .false.
    if ( present(allowneg) )then
      if (  allowneg ) lallowneg = .true.
    end if
    do fp = 1,num_soilp
       p = filter_soilp(fp)

       if ( .not. lcroponly .or. (patch%itype(p) >= nc3crop) ) then
          if ( .not. lallowneg .and. (carbon_patch(p) < cnegcrit) ) then
             write(iulog,*) 'ERROR: Carbon patch negative = ', carbon_patch(p)
             write(iulog,*) 'ERROR: limit = ', cnegcrit
             call endrun(msg='ERROR: carbon state critically negative '//errMsg(sourcefile, lineno))
          else if ( abs(carbon_patch(p)) < ccrit) then
             pc(p) = pc(p) + carbon_patch(p)
             carbon_patch(p) = 0._r8
   
             if ( use_c13 .and. present(c13) .and. present(pc13) ) then
                pc13(p) = pc13(p) + c13(p)
                c13(p) = 0._r8
             endif
             if ( use_c14 .and. present(c14)  .and. present(pc14)) then
                pc14(p) = pc14(p) + c14(p)
                c14(p) = 0._r8
             endif
          end if
       end if
    end do
 end subroutine TruncateCStates

 subroutine TruncateNStates( bounds, filter_soilp, num_soilp, nitrogen_patch, pn, lineno )
    !
    ! !DESCRIPTION:
    ! Truncate Nitrogen states. If a nitrogen state is too small truncate it to
    ! zero.
    !
    ! !USES:
    use abortutils , only : endrun
    use shr_log_mod, only : errMsg => shr_log_errMsg
    use clm_varctl , only : iulog
    use decompMod  , only : bounds_type
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type)              , intent(in)    :: bounds          ! bounds
    integer                        , intent(in)    :: num_soilp       ! number of soil patchs in filter
    integer                        , intent(in)    :: filter_soilp(:) ! filter for soil patches
    real(r8), intent(inout) :: nitrogen_patch(bounds%begp:)
    real(r8), intent(inout) :: pn(bounds%begp:)
    integer,  intent(in)    :: lineno

    integer :: fp, p

    SHR_ASSERT_ALL((ubound(nitrogen_patch) == (/bounds%endp/)), errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(pn)             == (/bounds%endp/)), errMsg(sourcefile, __LINE__))
    do fp = 1,num_soilp
       p = filter_soilp(fp)
       if ( nitrogen_patch(p) < nnegcrit ) then
          !write(iulog,*) 'WARNING: Nitrogen patch negative = ', nitrogen_patch
          !call endrun(msg='ERROR: nitrogen state critically negative'//errMsg(sourcefile, lineno))
       else if ( abs(nitrogen_patch(p)) < ncrit) then
          pn(p) = pn(p) + nitrogen_patch(p)
          nitrogen_patch(p) = 0._r8

       end if
    end do
 end subroutine TruncateNStates

end module CNPrecisionControlMod
