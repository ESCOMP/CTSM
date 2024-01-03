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
  use CropReprPoolsMod           , only : nrepr
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
  real(r8), public :: cnegcrit =  -6.e+1_r8             ! critical negative carbon state value for abort (gC/m2)
  real(r8), public :: ncrit    =  1.e-8_r8              ! critical nitrogen state value for truncation (gN/m2)
  real(r8), public :: nnegcrit =  -7.e+0_r8             ! critical negative nitrogen state value for abort (gN/m2)
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
  subroutine CNPrecisionControl(bounds, num_bgc_vegp, filter_bgc_vegp, &
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
    integer                        , intent(in)    :: num_bgc_vegp       ! number of bgc veg patches in filter
    integer                        , intent(in)    :: filter_bgc_vegp(:) ! filter for bgc veg patches
    type(cnveg_carbonstate_type)   , intent(inout) :: cnveg_carbonstate_inst
    type(cnveg_carbonstate_type)   , intent(inout) :: c13_cnveg_carbonstate_inst
    type(cnveg_carbonstate_type)   , intent(inout) :: c14_cnveg_carbonstate_inst
    type(cnveg_nitrogenstate_type) , intent(inout) :: cnveg_nitrogenstate_inst
    !
    ! !LOCAL VARIABLES:
    integer :: p,j,k                             ! indices
    integer :: fp                                ! filter indices
    integer :: num_truncatep                     ! number of points in filter_truncatep
    integer :: filter_truncatep(bounds%endp-bounds%begp+1) ! filter for points that need truncation
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
    ! cnveg_carbonstate_inst%reproductivec_patch                    Output:  [real(r8) (:,:)     ]  (gC/m2) grain C
    ! cnveg_carbonstate_inst%reproductivec_storage_patch            Output:  [real(r8) (:,:)     ]  (gC/m2) grain C storage
    ! cnveg_carbonstate_inst%reproductivec_xfer_patch               Output:  [real(r8) (:,:)     ]  (gC/m2) grain C transfer
    
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
    ! cnveg_nitrogenstate_inst%reproductiven_patch                  Output:  [real(r8) (:,:)     ]  (gC/m2) grain N
    ! cnveg_nitrogenstate_inst%reproductiven_storage_patch          Output:  [real(r8) (:,:)     ]  (gC/m2) grain N storage
    ! cnveg_nitrogenstate_inst%reproductiven_xfer_patch             Output:  [real(r8) (:,:)     ]  (gC/m2) grain N transfer
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
      do fp = 1,num_bgc_vegp
         p = filter_bgc_vegp(fp)

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
      call TruncateCandNStates( bounds, filter_bgc_vegp, num_bgc_vegp, cs%leafc_patch(bounds%begp:bounds%endp), &
                                ns%leafn_patch(bounds%begp:bounds%endp), &
                                pc(bounds%begp:), pn(bounds%begp:), __LINE__, &
                                num_truncatep, filter_truncatep)

      if (use_c13) then
         call TruncateAdditional( bounds, num_truncatep, filter_truncatep, &
                                  c13cs%leafc_patch(bounds%begp:bounds%endp), pc13(bounds%begp:bounds%endp), &
                                  __LINE__)
                         end if
      if (use_c14) then
         call TruncateAdditional( bounds, num_truncatep, filter_truncatep, &
                                  c14cs%leafc_patch(bounds%begp:bounds%endp), pc14(bounds%begp:bounds%endp), &
                                  __LINE__)
      end if


      ! leaf storage C and N
      call TruncateCandNStates( bounds, filter_bgc_vegp, num_bgc_vegp, cs%leafc_storage_patch(bounds%begp:bounds%endp), &
                                ns%leafn_storage_patch(bounds%begp:bounds%endp), pc(bounds%begp:), pn(bounds%begp:), __LINE__, &
                                num_truncatep, filter_truncatep)
      if (use_c13) then
         call TruncateAdditional( bounds, num_truncatep, filter_truncatep, &
                                  c13cs%leafc_storage_patch(bounds%begp:bounds%endp), pc13(bounds%begp:bounds%endp), &
                                  __LINE__)
      end if
      if (use_c14) then
         call TruncateAdditional( bounds, num_truncatep, filter_truncatep, &
                                 c14cs%leafc_storage_patch(bounds%begp:bounds%endp), pc14(bounds%begp:bounds%endp), &
                                 __LINE__)
      end if

      ! leaf transfer C and N
      call TruncateCandNStates( bounds, filter_bgc_vegp, num_bgc_vegp, cs%leafc_xfer_patch(bounds%begp:bounds%endp), &
                                ns%leafn_xfer_patch(bounds%begp:bounds%endp), pc(bounds%begp:), pn(bounds%begp:), __LINE__, &
                                num_truncatep, filter_truncatep)
      if (use_c13) then
         call TruncateAdditional( bounds, num_truncatep, filter_truncatep, &
                                  c13cs%leafc_xfer_patch(bounds%begp:bounds%endp), pc13(bounds%begp:bounds%endp), &
                                  __LINE__)
      end if
      if (use_c14) then
         call TruncateAdditional( bounds, num_truncatep, filter_truncatep, &
                                  c14cs%leafc_xfer_patch(bounds%begp:bounds%endp), pc14(bounds%begp:bounds%endp), &
                                  __LINE__)
      end if

      ! froot C and N
      ! EBK KO DML: For some reason frootc/frootn can go negative and allowing
      ! it to be negative is important for C4 crops (otherwise they die) Jun/3/2016
      if ( prec_control_for_froot ) then
         call TruncateCandNStates( bounds, filter_bgc_vegp, num_bgc_vegp, cs%frootc_patch(bounds%begp:bounds%endp),  &
                                   ns%frootn_patch(bounds%begp:bounds%endp), pc(bounds%begp:), pn(bounds%begp:), __LINE__, &
                                   num_truncatep, filter_truncatep, allowneg=.true.)
          if (use_c13) then
             call TruncateAdditional( bounds, num_truncatep, filter_truncatep, &
                                      c13cs%frootc_patch(bounds%begp:bounds%endp), pc13(bounds%begp:bounds%endp), &
                                      __LINE__)
          end if
          if (use_c14) then
             call TruncateAdditional( bounds, num_truncatep, filter_truncatep, &
                                      c14cs%frootc_patch(bounds%begp:bounds%endp), pc14(bounds%begp:bounds%endp), &
                                      __LINE__)
          end if
      end if

      ! froot storage C and N
      call TruncateCandNStates( bounds, filter_bgc_vegp, num_bgc_vegp, cs%frootc_storage_patch(bounds%begp:bounds%endp), &
                                ns%frootn_storage_patch(bounds%begp:bounds%endp), pc(bounds%begp:), pn(bounds%begp:), &
                                __LINE__, num_truncatep, filter_truncatep)
      if (use_c13) then
         call TruncateAdditional( bounds, num_truncatep, filter_truncatep, &
                                  c13cs%frootc_storage_patch(bounds%begp:bounds%endp), pc13(bounds%begp:bounds%endp), &
                                  __LINE__)
      end if
      if (use_c14) then
         call TruncateAdditional( bounds, num_truncatep, filter_truncatep, &
                                  c14cs%frootc_storage_patch(bounds%begp:bounds%endp), pc14(bounds%begp:bounds%endp), &
                                  __LINE__)
      end if

      ! froot transfer C and N
      call TruncateCandNStates( bounds, filter_bgc_vegp, num_bgc_vegp, cs%frootc_xfer_patch(bounds%begp:bounds%endp), &
                                ns%frootn_xfer_patch(bounds%begp:bounds%endp), pc(bounds%begp:), pn(bounds%begp:), __LINE__, &
                                num_truncatep, filter_truncatep)
      if (use_c13) then
         call TruncateAdditional( bounds, num_truncatep, filter_truncatep, &
                                  c13cs%frootc_xfer_patch(bounds%begp:bounds%endp), pc13(bounds%begp:bounds%endp), &
                                  __LINE__)
      end if
      if (use_c14) then
         call TruncateAdditional( bounds, num_truncatep, filter_truncatep, &
                                  c14cs%frootc_xfer_patch(bounds%begp:bounds%endp), pc14(bounds%begp:bounds%endp), &
                                  __LINE__)
      end if

      if ( use_crop )then
         do k = 1, nrepr
            ! grain C and N
            call TruncateCandNStates( bounds, filter_bgc_vegp, num_bgc_vegp, cs%reproductivec_patch(bounds%begp:bounds%endp,k), &
                 ns%reproductiven_patch(bounds%begp:bounds%endp,k), pc(bounds%begp:), pn(bounds%begp:), __LINE__, &
                 num_truncatep, filter_truncatep, croponly=.true. )
            if (use_c13) then
               call TruncateAdditional( bounds, num_truncatep, filter_truncatep, &
                    c13cs%reproductivec_patch(bounds%begp:bounds%endp,k), pc13(bounds%begp:bounds%endp), &
                    __LINE__)
            end if
            if (use_c14) then
               call TruncateAdditional( bounds, num_truncatep, filter_truncatep, &
                    c14cs%reproductivec_patch(bounds%begp:bounds%endp,k), pc14(bounds%begp:bounds%endp), &
                    __LINE__)
            end if

            ! grain storage C and N
            call TruncateCandNStates( bounds, filter_bgc_vegp, num_bgc_vegp, &
                 cs%reproductivec_storage_patch(bounds%begp:bounds%endp,k), &
                 ns%reproductiven_storage_patch(bounds%begp:bounds%endp,k), pc(bounds%begp:), pn(bounds%begp:), &
                 __LINE__, num_truncatep, filter_truncatep, croponly=.true. )

            if (use_c13) then
               call TruncateAdditional( bounds, num_truncatep, filter_truncatep, &
                    c13cs%reproductivec_storage_patch(bounds%begp:bounds%endp,k), pc13(bounds%begp:bounds%endp), &
                    __LINE__)
            end if
            if (use_c14) then
               call TruncateAdditional( bounds, num_truncatep, filter_truncatep, &
                    c14cs%reproductivec_storage_patch(bounds%begp:bounds%endp,k), pc14(bounds%begp:bounds%endp), &
                    __LINE__)
            end if

            ! grain transfer C and N
            call TruncateCandNStates( bounds, filter_bgc_vegp, num_bgc_vegp, &
                 cs%reproductivec_xfer_patch(bounds%begp:bounds%endp,k), &
                 ns%reproductiven_xfer_patch(bounds%begp:bounds%endp,k), pc(bounds%begp:), pn(bounds%begp:), __LINE__, &
                 num_truncatep, filter_truncatep, croponly=.true.)
            if (use_c13) then
               call TruncateAdditional( bounds, num_truncatep, filter_truncatep, &
                    c13cs%reproductivec_xfer_patch(bounds%begp:bounds%endp,k), pc13(bounds%begp:bounds%endp), &
                    __LINE__)
            end if
            if (use_c14) then
               call TruncateAdditional( bounds, num_truncatep, filter_truncatep, &
                    c14cs%reproductivec_xfer_patch(bounds%begp:bounds%endp,k), pc14(bounds%begp:bounds%endp), &
                    __LINE__)
            end if
         end do
         ! grain transfer C and N
         call TruncateCandNStates( bounds, filter_bgc_vegp, num_bgc_vegp, cs%cropseedc_deficit_patch(bounds%begp:bounds%endp), &
                                   ns%cropseedn_deficit_patch(bounds%begp:bounds%endp), pc(bounds%begp:), &
                                   pn(bounds%begp:), __LINE__, &
                                   num_truncatep, filter_truncatep, &
                                   allowneg=.true., croponly=.true. )
          if (use_c13) then
             call TruncateAdditional( bounds, num_truncatep, filter_truncatep, &
                                      c13cs%cropseedc_deficit_patch(bounds%begp:bounds%endp), pc13(bounds%begp:bounds%endp), &
                                      __LINE__)
          end if
          if (use_c14) then
             call TruncateAdditional( bounds, num_truncatep, filter_truncatep, &
                                      c14cs%cropseedc_deficit_patch(bounds%begp:bounds%endp), pc14(bounds%begp:bounds%endp), &
                                      __LINE__)
          end if

      end if

      ! livestem C and N
      call TruncateCandNStates( bounds, filter_bgc_vegp, num_bgc_vegp, cs%livestemc_patch(bounds%begp:bounds%endp), &
                                ns%livestemn_patch(bounds%begp:bounds%endp), pc(bounds%begp:), pn(bounds%begp:), __LINE__, &
                                num_truncatep, filter_truncatep)
      if (use_c13) then
         call TruncateAdditional( bounds, num_truncatep, filter_truncatep, &
                                  c13cs%livestemc_patch(bounds%begp:bounds%endp), pc13(bounds%begp:bounds%endp), &
                                  __LINE__)
      end if
      if (use_c14) then
         call TruncateAdditional( bounds, num_truncatep, filter_truncatep, &
                                  c14cs%livestemc_patch(bounds%begp:bounds%endp), pc14(bounds%begp:bounds%endp), &
                                  __LINE__)
      end if

      ! livestem storage C and N
      call TruncateCandNStates( bounds, filter_bgc_vegp, num_bgc_vegp, cs%livestemc_storage_patch(bounds%begp:bounds%endp), &
                                ns%livestemn_storage_patch(bounds%begp:bounds%endp), pc(bounds%begp:), pn(bounds%begp:), &
                                __LINE__, num_truncatep, filter_truncatep)

      if (use_c13) then
         call TruncateAdditional( bounds, num_truncatep, filter_truncatep, &
                                  c13cs%livestemc_storage_patch(bounds%begp:bounds%endp), pc13(bounds%begp:bounds%endp), &
                                  __LINE__)
      end if
      if (use_c14) then
         call TruncateAdditional( bounds, num_truncatep, filter_truncatep, &
                                  c14cs%livestemc_storage_patch(bounds%begp:bounds%endp), pc14(bounds%begp:bounds%endp), &
                                  __LINE__)
      end if
      ! livestem transfer C and N
      call TruncateCandNStates( bounds, filter_bgc_vegp, num_bgc_vegp, cs%livestemc_xfer_patch(bounds%begp:bounds%endp), &
                                ns%livestemn_xfer_patch(bounds%begp:bounds%endp), pc(bounds%begp:), pn(bounds%begp:), &
                                __LINE__, num_truncatep, filter_truncatep)
      if (use_c13) then
         call TruncateAdditional( bounds, num_truncatep, filter_truncatep, &
                                  c13cs%livestemc_xfer_patch(bounds%begp:bounds%endp), pc13(bounds%begp:bounds%endp), &
                                  __LINE__)
      end if
      if (use_c14) then
         call TruncateAdditional( bounds, num_truncatep, filter_truncatep, &
                                  c14cs%livestemc_xfer_patch(bounds%begp:bounds%endp), pc14(bounds%begp:bounds%endp), &
                                  __LINE__)
      end if

      ! deadstem C and N
      call TruncateCandNStates( bounds, filter_bgc_vegp, num_bgc_vegp, cs%deadstemc_patch(bounds%begp:bounds%endp), &
                                ns%deadstemn_patch(bounds%begp:bounds%endp), pc(bounds%begp:), pn(bounds%begp:), __LINE__, &
                                num_truncatep, filter_truncatep)
      if (use_c13) then
         call TruncateAdditional( bounds, num_truncatep, filter_truncatep, &
                                  c13cs%deadstemc_patch(bounds%begp:bounds%endp), pc13(bounds%begp:bounds%endp), &
                                  __LINE__)
      end if
      if (use_c14) then
         call TruncateAdditional( bounds, num_truncatep, filter_truncatep, &
                                  c14cs%deadstemc_patch(bounds%begp:bounds%endp), pc14(bounds%begp:bounds%endp), &
                                  __LINE__)
      end if
      ! deadstem storage C and N
      call TruncateCandNStates( bounds, filter_bgc_vegp, num_bgc_vegp, cs%deadstemc_storage_patch(bounds%begp:bounds%endp), &
                                ns%deadstemn_storage_patch(bounds%begp:bounds%endp), pc(bounds%begp:), pn(bounds%begp:), &
                                __LINE__, num_truncatep, filter_truncatep)
      if (use_c13) then
         call TruncateAdditional( bounds, num_truncatep, filter_truncatep, &
                                  c13cs%deadstemc_storage_patch(bounds%begp:bounds%endp), pc13(bounds%begp:bounds%endp), &
                                  __LINE__)
      end if
      if (use_c14) then
         call TruncateAdditional( bounds, num_truncatep, filter_truncatep, &
                                  c14cs%deadstemc_storage_patch(bounds%begp:bounds%endp), pc14(bounds%begp:bounds%endp), &
                                  __LINE__)
      end if

      ! deadstem transfer C and N
      call TruncateCandNStates( bounds, filter_bgc_vegp, num_bgc_vegp, cs%deadstemc_xfer_patch(bounds%begp:bounds%endp), &
                                ns%deadstemn_xfer_patch(bounds%begp:bounds%endp), pc(bounds%begp:), pn(bounds%begp:), &
                                __LINE__, num_truncatep, filter_truncatep)
      if (use_c13) then
         call TruncateAdditional( bounds, num_truncatep, filter_truncatep, &
                                  c13cs%deadstemc_xfer_patch(bounds%begp:bounds%endp), pc13(bounds%begp:bounds%endp), &
                                  __LINE__)
      end if
      if (use_c14) then
         call TruncateAdditional( bounds, num_truncatep, filter_truncatep, &
                                  c14cs%deadstemc_xfer_patch(bounds%begp:bounds%endp), pc14(bounds%begp:bounds%endp), &
                                  __LINE__)
      end if

      ! livecroot C and N
      call TruncateCandNStates( bounds, filter_bgc_vegp, num_bgc_vegp, cs%livecrootc_patch(bounds%begp:bounds%endp), &
                                ns%livecrootn_patch(bounds%begp:bounds%endp), pc(bounds%begp:), pn(bounds%begp:), __LINE__, &
                                num_truncatep, filter_truncatep)
      if (use_c13) then
         call TruncateAdditional( bounds, num_truncatep, filter_truncatep, &
                                  c13cs%livecrootc_patch(bounds%begp:bounds%endp), pc13(bounds%begp:bounds%endp), &
                                  __LINE__)
      end if
      if (use_c14) then
         call TruncateAdditional( bounds, num_truncatep, filter_truncatep, &
                                  c14cs%livecrootc_patch(bounds%begp:bounds%endp), pc14(bounds%begp:bounds%endp), &
                                  __LINE__)
      end if

      ! livecroot storage C and N
      call TruncateCandNStates( bounds, filter_bgc_vegp, num_bgc_vegp, cs%livecrootc_storage_patch(bounds%begp:bounds%endp), &
                                ns%livecrootn_storage_patch(bounds%begp:bounds%endp), pc(bounds%begp:), pn(bounds%begp:), &
                                __LINE__, num_truncatep, filter_truncatep)
      if (use_c13) then
         call TruncateAdditional( bounds, num_truncatep, filter_truncatep, &
                                  c13cs%livecrootc_storage_patch(bounds%begp:bounds%endp), pc13(bounds%begp:bounds%endp), &
                                  __LINE__)
      end if
      if (use_c14) then
         call TruncateAdditional( bounds, num_truncatep, filter_truncatep, &
                                  c14cs%livecrootc_storage_patch(bounds%begp:bounds%endp), pc14(bounds%begp:bounds%endp), &
                                  __LINE__)
      end if

      ! livecroot transfer C and N
      call TruncateCandNStates( bounds, filter_bgc_vegp, num_bgc_vegp, cs%livecrootc_xfer_patch(bounds%begp:bounds%endp), &
                                ns%livecrootn_xfer_patch(bounds%begp:bounds%endp), pc(bounds%begp:), pn(bounds%begp:), &
                                __LINE__, num_truncatep, filter_truncatep)

      if (use_c13) then
         call TruncateAdditional( bounds, num_truncatep, filter_truncatep, &
                                  c13cs%livecrootc_xfer_patch(bounds%begp:bounds%endp), pc13(bounds%begp:bounds%endp), &
                                  __LINE__)
      end if
      if (use_c14) then
         call TruncateAdditional( bounds, num_truncatep, filter_truncatep, &
                                  c14cs%livecrootc_xfer_patch(bounds%begp:bounds%endp), pc14(bounds%begp:bounds%endp), &
                                  __LINE__)
      end if

      ! deadcroot C and N
      call TruncateCandNStates( bounds, filter_bgc_vegp, num_bgc_vegp, cs%deadcrootc_patch(bounds%begp:bounds%endp), &
                                ns%deadcrootn_patch(bounds%begp:bounds%endp), pc(bounds%begp:), pn(bounds%begp:), __LINE__, &
                                num_truncatep, filter_truncatep)
      if (use_c13) then
         call TruncateAdditional( bounds, num_truncatep, filter_truncatep, &
                                  c13cs%deadcrootc_patch(bounds%begp:bounds%endp), pc13(bounds%begp:bounds%endp), &
                                  __LINE__)
      end if
      if (use_c14) then
         call TruncateAdditional( bounds, num_truncatep, filter_truncatep, &
                                  c14cs%deadcrootc_patch(bounds%begp:bounds%endp), pc14(bounds%begp:bounds%endp), &
                                  __LINE__)
      end if

      ! deadcroot storage C and N
      call TruncateCandNStates( bounds, filter_bgc_vegp, num_bgc_vegp, cs%deadcrootc_storage_patch(bounds%begp:bounds%endp), &
                                ns%deadcrootn_storage_patch(bounds%begp:bounds%endp), pc(bounds%begp:), pn(bounds%begp:), &
                                __LINE__, num_truncatep, filter_truncatep)
      if (use_c13) then
         call TruncateAdditional( bounds, num_truncatep, filter_truncatep, &
                                  c13cs%deadcrootc_storage_patch(bounds%begp:bounds%endp), pc13(bounds%begp:bounds%endp), &
                                  __LINE__)
      end if
      if (use_c14) then
         call TruncateAdditional( bounds, num_truncatep, filter_truncatep, &
                                  c14cs%deadcrootc_storage_patch(bounds%begp:bounds%endp), pc14(bounds%begp:bounds%endp), &
                                  __LINE__)
      end if

      ! deadcroot transfer C and N
      call TruncateCandNStates( bounds, filter_bgc_vegp, num_bgc_vegp, cs%deadcrootc_xfer_patch(bounds%begp:bounds%endp), &
                                ns%deadcrootn_xfer_patch(bounds%begp:bounds%endp), pc(bounds%begp:), pn(bounds%begp:), &
                                __LINE__, num_truncatep, filter_truncatep)
      if (use_c13) then
         call TruncateAdditional( bounds, num_truncatep, filter_truncatep, &
                                  c13cs%deadcrootc_xfer_patch(bounds%begp:bounds%endp), pc13(bounds%begp:bounds%endp), &
                                  __LINE__)
      end if
      if (use_c14) then
         call TruncateAdditional( bounds, num_truncatep, filter_truncatep, &
                                  c14cs%deadcrootc_xfer_patch(bounds%begp:bounds%endp), pc14(bounds%begp:bounds%endp), &
                                  __LINE__)
      end if

      ! gresp_storage (C only)
      call TruncateCStates( bounds, filter_bgc_vegp, num_bgc_vegp, cs%gresp_storage_patch(bounds%begp:bounds%endp), &
                            pc(bounds%begp:), __LINE__, num_truncatep, filter_truncatep)
      if (use_c13) then
         call TruncateAdditional( bounds, num_truncatep, filter_truncatep, &
                                  c13cs%gresp_storage_patch(bounds%begp:bounds%endp), pc13(bounds%begp:bounds%endp), &
                                  __LINE__)
      end if
      if (use_c14) then
         call TruncateAdditional( bounds, num_truncatep, filter_truncatep, &
                                  c14cs%gresp_storage_patch(bounds%begp:bounds%endp), pc14(bounds%begp:bounds%endp), &
                                  __LINE__)
      end if

      ! gresp_xfer(c only)
      call TruncateCStates( bounds, filter_bgc_vegp, num_bgc_vegp, cs%gresp_xfer_patch(bounds%begp:bounds%endp), &
                            pc(bounds%begp:), __LINE__, num_truncatep, filter_truncatep)
      if (use_c13) then
         call TruncateAdditional( bounds, num_truncatep, filter_truncatep, &
                                  c13cs%gresp_xfer_patch(bounds%begp:bounds%endp), pc13(bounds%begp:bounds%endp), &
                                  __LINE__)
      end if
      if (use_c14) then
         call TruncateAdditional( bounds, num_truncatep, filter_truncatep, &
                                  c14cs%gresp_xfer_patch(bounds%begp:bounds%endp), pc14(bounds%begp:bounds%endp), &
                                  __LINE__)
      end if

      ! cpool (C only)
      call TruncateCStates( bounds, filter_bgc_vegp, num_bgc_vegp, cs%cpool_patch(bounds%begp:bounds%endp), &
                            pc(bounds%begp:), __LINE__, num_truncatep, filter_truncatep)
      if (use_c13) then
         call TruncateAdditional( bounds, num_truncatep, filter_truncatep, &
                                  c13cs%cpool_patch(bounds%begp:bounds%endp), pc13(bounds%begp:bounds%endp), &
                                  __LINE__)
      end if
      if (use_c14) then
         call TruncateAdditional( bounds, num_truncatep, filter_truncatep, &
                                  c14cs%cpool_patch(bounds%begp:bounds%endp), pc14(bounds%begp:bounds%endp), &
                                  __LINE__)
      end if

      if ( use_crop )then
         ! xsmrpool (C only)
         ! xsmr is a pool to balance the budget and as such can be freely negative
         call TruncateCStates( bounds, filter_bgc_vegp, num_bgc_vegp, cs%xsmrpool_patch(bounds%begp:bounds%endp), &
                                pc(bounds%begp:), __LINE__, num_truncatep, filter_truncatep, &
                                allowneg=.true., croponly=.true. )
         if (use_c13) then
             call TruncateAdditional( bounds, num_truncatep, filter_truncatep, &
                                      c13cs%xsmrpool_patch(bounds%begp:bounds%endp), pc13(bounds%begp:bounds%endp), &
                                      __LINE__)
         end if
         if (use_c14) then
             call TruncateAdditional( bounds, num_truncatep, filter_truncatep, &
                                      c14cs%xsmrpool_patch(bounds%begp:bounds%endp), pc14(bounds%begp:bounds%endp), &
                                      __LINE__)
         end if

      end if

      ! retransn (N only)
      call TruncateNStates( bounds, filter_bgc_vegp, num_bgc_vegp, ns%retransn_patch(bounds%begp:bounds%endp), pn(bounds%begp:), &
                            __LINE__ )

      ! npool (N only)
      call TruncateNStates( bounds, filter_bgc_vegp, num_bgc_vegp, ns%npool_patch(bounds%begp:bounds%endp), pn(bounds%begp:), &
                            __LINE__ )

      ! patch loop
      do fp = 1,num_bgc_vegp
         p = filter_bgc_vegp(fp)

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

 subroutine TruncateCandNStates( bounds, filter_bgc_vegp, num_bgc_vegp, carbon_patch, nitrogen_patch, pc, pn, lineno, &
                                 num_truncatep, filter_truncatep, croponly, allowneg )
    !
    ! !DESCRIPTION:
    ! Truncate paired Carbon and Nitrogen states. If a paired carbon and nitrogen state iare too small truncate 
    ! the pair of them to zero.
    !
    ! !USES:
    use shr_log_mod, only : errMsg => shr_log_errMsg
    use clm_varctl , only : use_c13, use_c14, use_nguardrail
    use CNSharedParamsMod, only : use_matrixcn
    use clm_varctl , only : iulog
    use pftconMod  , only : nc3crop
    use decompMod  , only : bounds_type, subgrid_level_patch
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type)              , intent(in)    :: bounds          ! bounds
    integer                        , intent(in)    :: num_bgc_vegp       ! number of bgc veg patches in filter
    integer                        , intent(in)    :: filter_bgc_vegp(:) ! filter for bgc veg patches
    real(r8), intent(inout) :: carbon_patch(bounds%begp:)
    real(r8), intent(inout) :: nitrogen_patch(bounds%begp:)
    real(r8), intent(inout) :: pc(bounds%begp:)
    real(r8), intent(inout) :: pn(bounds%begp:)
    integer,  intent(in)    :: lineno
    integer,  intent(out)   :: num_truncatep       ! number of points in filter_truncatep
    integer,  intent(out)   :: filter_truncatep(:) ! filter for points that need truncation
    logical , intent(in)   , optional :: croponly
    logical , intent(in)   , optional :: allowneg

    logical :: lcroponly, lallowneg
    integer :: fp, p

    SHR_ASSERT_ALL_FL((ubound(carbon_patch)   == (/bounds%endp/)), 'ubnd(carb)'//sourcefile, lineno)
    SHR_ASSERT_ALL_FL((ubound(nitrogen_patch) == (/bounds%endp/)), 'ubnd(nitro)'//sourcefile, lineno)
    SHR_ASSERT_ALL_FL((ubound(pc)             == (/bounds%endp/)), 'ubnd(pc)'//sourcefile, lineno)
    SHR_ASSERT_ALL_FL((ubound(pn)             == (/bounds%endp/)), 'ubnd(pn)'//sourcefile, lineno)

    ! patch loop
    lcroponly = .false.
    if ( present(croponly) )then
      if ( croponly ) lcroponly = .true.
    end if
    lallowneg = .false.
    if ( present(allowneg) )then
      if (  allowneg ) lallowneg = .true.
    end if

    num_truncatep = 0
    do fp = 1,num_bgc_vegp
       p = filter_bgc_vegp(fp)

       if ( .not. lcroponly .or. (patch%itype(p) >= nc3crop) ) then
          if ( .not. lallowneg .and. ((carbon_patch(p) < cnegcrit) .or. (nitrogen_patch(p) < nnegcrit)) ) then
             write(iulog,*) 'ERROR: Carbon or Nitrogen patch negative = ', carbon_patch(p), nitrogen_patch(p)
             write(iulog,*) 'ERROR: limits = ', cnegcrit, nnegcrit
             call endrun(subgrid_index=p, subgrid_level=subgrid_level_patch, &
                  msg='ERROR: carbon or nitrogen state critically negative '//errMsg(sourcefile, lineno))
          else 
             if (use_matrixcn)then
                ! The matrix code has a different check here
                if ( (carbon_patch(p) < ccrit .and. carbon_patch(p) > -ccrit * 1.e+6) .or. &
                     (use_nguardrail .and. nitrogen_patch(p) < ncrit .and. nitrogen_patch(p) > -ncrit*1.e+6) ) then
                   ! This part needs to be identical to the part for non-matrix
                   ! below
                   num_truncatep = num_truncatep + 1
                   filter_truncatep(num_truncatep) = p

                   pc(p) = pc(p) + carbon_patch(p)
                   carbon_patch(p) = 0._r8
          
                   pn(p) = pn(p) + nitrogen_patch(p)
                   nitrogen_patch(p) = 0._r8
    
                end if
             else
                if ( abs(carbon_patch(p)) < ccrit .or. (use_nguardrail .and. abs(nitrogen_patch(p)) < ncrit) ) then
                   ! This part needs to be identical to above
                   num_truncatep = num_truncatep + 1
                   filter_truncatep(num_truncatep) = p

                   pc(p) = pc(p) + carbon_patch(p)
                   carbon_patch(p) = 0._r8

                   pn(p) = pn(p) + nitrogen_patch(p)
                   nitrogen_patch(p) = 0._r8

                end if
             end if
          end if
       end if
    end do
 end subroutine TruncateCandNStates

 subroutine TruncateCStates( bounds, filter_bgc_vegp, num_bgc_vegp, carbon_patch, pc, lineno,  &
                             num_truncatep, filter_truncatep, croponly, allowneg )
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
    use decompMod  , only : bounds_type, subgrid_level_patch
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in)    :: bounds          ! bounds
    integer          , intent(in)    :: num_bgc_vegp       ! number of bgc veg patches in filter
    integer          , intent(in)    :: filter_bgc_vegp(:) ! filter for bgc veg patches
    real(r8)         , intent(inout) :: carbon_patch(bounds%begp:)
    real(r8)         , intent(inout) :: pc(bounds%begp:)
    integer          , intent(in)    :: lineno
    integer          , intent(out)   :: num_truncatep       ! number of points in filter_truncatep
    integer          , intent(out)   :: filter_truncatep(:) ! filter for points that need truncation
    logical          , intent(in)   , optional :: croponly
    logical          , intent(in)   , optional :: allowneg

    logical :: lcroponly, lallowneg
    integer :: fp, p

    SHR_ASSERT_ALL_FL((ubound(carbon_patch)   == (/bounds%endp/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(pc)             == (/bounds%endp/)), sourcefile, __LINE__)

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

    num_truncatep = 0
    do fp = 1,num_bgc_vegp
       p = filter_bgc_vegp(fp)

       if ( .not. lcroponly .or. (patch%itype(p) >= nc3crop) ) then
          if ( .not. lallowneg .and. (carbon_patch(p) < cnegcrit) ) then
             write(iulog,*) 'ERROR: Carbon patch negative = ', carbon_patch(p)
             write(iulog,*) 'ERROR: limit = ', cnegcrit
             call endrun(subgrid_index=p, subgrid_level=subgrid_level_patch, &
                  msg='ERROR: carbon state critically negative '//errMsg(sourcefile, lineno))
          else if ( abs(carbon_patch(p)) < ccrit) then

             num_truncatep = num_truncatep + 1
             filter_truncatep(num_truncatep) = p

             pc(p) = pc(p) + carbon_patch(p)
             carbon_patch(p) = 0._r8
          end if
       end if
    end do
 end subroutine TruncateCStates

 subroutine TruncateNStates( bounds, filter_bgc_vegp, num_bgc_vegp, nitrogen_patch, pn, lineno )
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
    integer                        , intent(in)    :: num_bgc_vegp       ! number of bgc veg patches in filter
    integer                        , intent(in)    :: filter_bgc_vegp(:) ! filter for bgc veg patches
    real(r8), intent(inout) :: nitrogen_patch(bounds%begp:)
    real(r8), intent(inout) :: pn(bounds%begp:)
    integer,  intent(in)    :: lineno

    integer :: fp, p

    SHR_ASSERT_ALL_FL((ubound(nitrogen_patch) == (/bounds%endp/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(pn)             == (/bounds%endp/)), sourcefile, __LINE__)
    do fp = 1,num_bgc_vegp
       p = filter_bgc_vegp(fp)
       if ( nitrogen_patch(p) < nnegcrit ) then
          ! write(iulog,*) 'WARNING: Nitrogen patch negative = ', nitrogen_patch
          ! call endrun(subgrid_index=p, subgrid_level=subgrid_level_patch, &
          !      msg='ERROR: nitrogen state critically negative'//errMsg(sourcefile, lineno))
       else if ( abs(nitrogen_patch(p)) < ncrit) then
          pn(p) = pn(p) + nitrogen_patch(p)
          nitrogen_patch(p) = 0._r8

       end if
    end do
 end subroutine TruncateNStates

 !-----------------------------------------------------------------------
 subroutine TruncateAdditional( bounds, num_truncatep, filter_truncatep, &
                                state_patch, truncation_patch, lineno)
   !
   ! !DESCRIPTION:
   ! Given a filter of points for which we have already determined that truncation should
   ! occur, do the truncation for the given patch-level state, putting the truncation
   ! amount in truncation_patch.
   !
   use decompMod  , only : bounds_type
   ! !ARGUMENTS:
   implicit none
   type(bounds_type) , intent (in)    :: bounds              ! bounds
   integer           , intent (in)    :: num_truncatep       ! number of points in filter_truncatep
   integer           , intent (in)    :: filter_truncatep(:) ! filter for points that need truncation
   real(r8)          , intent (inout) :: state_patch(bounds%begp: )
   real(r8)          , intent (inout) :: truncation_patch(bounds%begp: )
   integer           , intent (in)    :: lineno
   !
   ! !LOCAL VARIABLES:
   integer                     :: fp, p
   character(len=*), parameter :: subname = 'TruncateAdditional'
   !-----------------------------------------------------------------------

   SHR_ASSERT_FL((ubound(state_patch, 1)      == bounds%endp), 'state_patch '     //sourcefile, lineno)
   SHR_ASSERT_FL((ubound(truncation_patch, 1) == bounds%endp), 'truncation_patch '//sourcefile, lineno)

   do fp = 1, num_truncatep
      p = filter_truncatep(fp)
      truncation_patch(p) = truncation_patch(p) + state_patch(p)
      state_patch(p) = 0._r8
   end do

 end subroutine TruncateAdditional

end module CNPrecisionControlMod
