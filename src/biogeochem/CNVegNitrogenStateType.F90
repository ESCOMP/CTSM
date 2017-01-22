module CNVegNitrogenStateType

#include "shr_assert.h"

  use shr_kind_mod                       , only : r8 => shr_kind_r8
  use shr_infnan_mod                     , only : isnan => shr_infnan_isnan, nan => shr_infnan_nan, assignment(=)
  use shr_log_mod                        , only : errMsg => shr_log_errMsg
  use clm_varpar                         , only : ndecomp_cascade_transitions, ndecomp_pools, nlevcan
  use clm_varpar                         , only : nlevdecomp_full, nlevdecomp
  use clm_varcon                         , only : spval, ispval, dzsoi_decomp, zisoi
  use landunit_varcon                    , only : istcrop, istsoil 
  use clm_varctl                         , only : use_nitrif_denitrif, use_vertsoilc, use_century_decomp
  use clm_varctl                         , only : iulog, override_bgc_restart_mismatch_dump
  use clm_varctl                         , only : use_crop
  use CNSharedParamsMod                  , only : use_fun
  use decompMod                          , only : bounds_type
  use pftconMod                          , only : npcropmin, noveg, pftcon
  use SoilBiogeochemDecompCascadeConType , only : decomp_cascade_con
  use abortutils                         , only : endrun
  use spmdMod                            , only : masterproc 
  use LandunitType                       , only : lun                
  use ColumnType                         , only : col                
  use PatchType                          , only : patch                
  use dynPatchStateUpdaterMod, only : patch_state_updater_type
  use CNSpeciesMod   , only : CN_SPECIES_N
  use CNVegComputeSeedMod, only : ComputeSeedAmounts
  !
  ! !PUBLIC TYPES:
  implicit none

  private


  !
  type, public :: cnveg_nitrogenstate_type

     real(r8), pointer :: grainn_patch             (:) ! (gN/m2) grain N (crop)
     real(r8), pointer :: grainn_storage_patch     (:) ! (gN/m2) grain N storage (crop)
     real(r8), pointer :: grainn_xfer_patch        (:) ! (gN/m2) grain N transfer (crop)
     real(r8), pointer :: leafn_patch              (:) ! (gN/m2) leaf N 
     real(r8), pointer :: leafn_storage_patch      (:) ! (gN/m2) leaf N storage
     real(r8), pointer :: leafn_xfer_patch         (:) ! (gN/m2) leaf N transfer
     real(r8), pointer :: leafn_storage_xfer_acc_patch (:) ! (gN/m2) Accmulated leaf N transfer
     real(r8), pointer :: storage_ndemand_patch        (:) ! (gN/m2) N demand during the offset period 
     real(r8), pointer :: frootn_patch             (:) ! (gN/m2) fine root N
     real(r8), pointer :: frootn_storage_patch     (:) ! (gN/m2) fine root N storage
     real(r8), pointer :: frootn_xfer_patch        (:) ! (gN/m2) fine root N transfer
     real(r8), pointer :: livestemn_patch          (:) ! (gN/m2) live stem N
     real(r8), pointer :: livestemn_storage_patch  (:) ! (gN/m2) live stem N storage
     real(r8), pointer :: livestemn_xfer_patch     (:) ! (gN/m2) live stem N transfer
     real(r8), pointer :: deadstemn_patch          (:) ! (gN/m2) dead stem N
     real(r8), pointer :: deadstemn_storage_patch  (:) ! (gN/m2) dead stem N storage
     real(r8), pointer :: deadstemn_xfer_patch     (:) ! (gN/m2) dead stem N transfer
     real(r8), pointer :: livecrootn_patch         (:) ! (gN/m2) live coarse root N
     real(r8), pointer :: livecrootn_storage_patch (:) ! (gN/m2) live coarse root N storage
     real(r8), pointer :: livecrootn_xfer_patch    (:) ! (gN/m2) live coarse root N transfer
     real(r8), pointer :: deadcrootn_patch         (:) ! (gN/m2) dead coarse root N
     real(r8), pointer :: deadcrootn_storage_patch (:) ! (gN/m2) dead coarse root N storage
     real(r8), pointer :: deadcrootn_xfer_patch    (:) ! (gN/m2) dead coarse root N transfer
     real(r8), pointer :: retransn_patch           (:) ! (gN/m2) plant pool of retranslocated N
     real(r8), pointer :: npool_patch              (:) ! (gN/m2) temporary plant N pool
     real(r8), pointer :: ntrunc_patch             (:) ! (gN/m2) patch-level sink for N truncation
     real(r8), pointer :: cropseedn_deficit_patch  (:) ! (gN/m2) pool for seeding new crop growth; this is a NEGATIVE term, indicating the amount of seed usage that needs to be repaid
     real(r8), pointer :: seedn_grc                (:) ! (gN/m2) gridcell-level pool for seeding new pFTs via dynamic landcover

     ! summary (diagnostic) state variables, not involved in mass balance
     real(r8), pointer :: dispvegn_patch           (:) ! (gN/m2) displayed veg nitrogen, excluding storage
     real(r8), pointer :: storvegn_patch           (:) ! (gN/m2) stored vegetation nitrogen
     real(r8), pointer :: totvegn_patch            (:) ! (gN/m2) total vegetation nitrogen
     real(r8), pointer :: totvegn_col              (:) ! (gN/m2) total vegetation nitrogen (p2c)
     real(r8), pointer :: totn_patch               (:) ! (gN/m2) total patch-level nitrogen
     real(r8), pointer :: totn_p2c_col             (:) ! (gN/m2) totn_patch averaged to col
     real(r8), pointer :: totn_col                 (:) ! (gN/m2) total column nitrogen, incl veg
     real(r8), pointer :: totecosysn_col           (:) ! (gN/m2) total ecosystem nitrogen, incl veg  

   contains

     procedure , public  :: Init   
     procedure , public  :: Restart
     procedure , public  :: SetValues
     procedure , public  :: ZeroDWT
     procedure , public  :: Summary => Summary_nitrogenstate
     procedure , public  :: DynamicPatchAdjustments   ! adjust state variables when patch areas change
     procedure , private :: InitAllocate 
     procedure , private :: InitHistory  
     procedure , private :: InitCold     

  end type cnveg_nitrogenstate_type
  !------------------------------------------------------------------------

  ! !PRIVATE DATA:
  character(len=*), parameter :: sourcefile = &
       __FILE__

contains

  !------------------------------------------------------------------------
  subroutine Init(this, bounds,                           &
       leafc_patch, leafc_storage_patch, frootc_patch, frootc_storage_patch, deadstemc_patch)

    class(cnveg_nitrogenstate_type)   :: this
    type(bounds_type) , intent(in)    :: bounds  
    real(r8)          , intent(in)    :: leafc_patch          (bounds%begp:)
    real(r8)          , intent(in)    :: leafc_storage_patch  (bounds%begp:)
    real(r8)          , intent(in)    :: frootc_patch         (bounds%begp:)     
    real(r8)          , intent(in)    :: frootc_storage_patch (bounds%begp:)     
    real(r8)          , intent(in)    :: deadstemc_patch      (bounds%begp:)

    call this%InitAllocate (bounds )
    call this%InitHistory (bounds)
    call this%InitCold ( bounds, &
         leafc_patch, leafc_storage_patch, frootc_patch, frootc_storage_patch, deadstemc_patch)

  end subroutine Init

  !------------------------------------------------------------------------
  subroutine InitAllocate(this, bounds)
    !
    ! !ARGUMENTS:
    class (cnveg_nitrogenstate_type) :: this
    type(bounds_type) , intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer           :: begp,endp
    integer           :: begc,endc
    integer           :: begg,endg
    !------------------------------------------------------------------------

    begp = bounds%begp; endp = bounds%endp
    begc = bounds%begc; endc = bounds%endc
    begg = bounds%begg; endg = bounds%endg

    allocate(this%grainn_patch             (begp:endp)) ; this%grainn_patch             (:) = nan     
    allocate(this%grainn_storage_patch     (begp:endp)) ; this%grainn_storage_patch     (:) = nan
    allocate(this%grainn_xfer_patch        (begp:endp)) ; this%grainn_xfer_patch        (:) = nan     
    allocate(this%leafn_patch              (begp:endp)) ; this%leafn_patch              (:) = nan
    allocate(this%leafn_storage_patch      (begp:endp)) ; this%leafn_storage_patch      (:) = nan     
    allocate(this%leafn_xfer_patch         (begp:endp)) ; this%leafn_xfer_patch         (:) = nan     
    allocate(this%leafn_storage_xfer_acc_patch  (begp:endp)) ; this%leafn_storage_xfer_acc_patch         (:) = nan
    allocate(this%storage_ndemand_patch    (begp:endp)) ; this%storage_ndemand_patch    (:) = nan
    allocate(this%frootn_patch             (begp:endp)) ; this%frootn_patch             (:) = nan
    allocate(this%frootn_storage_patch     (begp:endp)) ; this%frootn_storage_patch     (:) = nan     
    allocate(this%frootn_xfer_patch        (begp:endp)) ; this%frootn_xfer_patch        (:) = nan     
    allocate(this%livestemn_patch          (begp:endp)) ; this%livestemn_patch          (:) = nan
    allocate(this%livestemn_storage_patch  (begp:endp)) ; this%livestemn_storage_patch  (:) = nan
    allocate(this%livestemn_xfer_patch     (begp:endp)) ; this%livestemn_xfer_patch     (:) = nan
    allocate(this%deadstemn_patch          (begp:endp)) ; this%deadstemn_patch          (:) = nan
    allocate(this%deadstemn_storage_patch  (begp:endp)) ; this%deadstemn_storage_patch  (:) = nan
    allocate(this%deadstemn_xfer_patch     (begp:endp)) ; this%deadstemn_xfer_patch     (:) = nan
    allocate(this%livecrootn_patch         (begp:endp)) ; this%livecrootn_patch         (:) = nan
    allocate(this%livecrootn_storage_patch (begp:endp)) ; this%livecrootn_storage_patch (:) = nan
    allocate(this%livecrootn_xfer_patch    (begp:endp)) ; this%livecrootn_xfer_patch    (:) = nan
    allocate(this%deadcrootn_patch         (begp:endp)) ; this%deadcrootn_patch         (:) = nan
    allocate(this%deadcrootn_storage_patch (begp:endp)) ; this%deadcrootn_storage_patch (:) = nan
    allocate(this%deadcrootn_xfer_patch    (begp:endp)) ; this%deadcrootn_xfer_patch    (:) = nan
    allocate(this%retransn_patch           (begp:endp)) ; this%retransn_patch           (:) = nan
    allocate(this%npool_patch              (begp:endp)) ; this%npool_patch              (:) = nan
    allocate(this%ntrunc_patch             (begp:endp)) ; this%ntrunc_patch             (:) = nan
    allocate(this%dispvegn_patch           (begp:endp)) ; this%dispvegn_patch           (:) = nan
    allocate(this%storvegn_patch           (begp:endp)) ; this%storvegn_patch           (:) = nan
    allocate(this%totvegn_patch            (begp:endp)) ; this%totvegn_patch            (:) = nan
    allocate(this%totn_patch               (begp:endp)) ; this%totn_patch               (:) = nan

    allocate(this%cropseedn_deficit_patch  (begp:endp)) ; this%cropseedn_deficit_patch  (:) = nan
    allocate(this%seedn_grc                (begg:endg)) ; this%seedn_grc                (:) = nan
    allocate(this%totvegn_col              (begc:endc)) ; this%totvegn_col              (:) = nan
    allocate(this%totn_p2c_col             (begc:endc)) ; this%totn_p2c_col             (:) = nan
    allocate(this%totn_col                 (begc:endc)) ; this%totn_col                 (:) = nan
    allocate(this%totecosysn_col           (begc:endc)) ; this%totecosysn_col           (:) = nan

  end subroutine InitAllocate

  !------------------------------------------------------------------------
  subroutine InitHistory(this, bounds)
    !
    ! !DESCRIPTION:
    ! add history fields for all CN variables, always set as default='inactive'
    !
    ! !USES:
    use histFileMod, only : hist_addfld1d
    !
    ! !ARGUMENTS:
    class(cnveg_nitrogenstate_type) :: this
    type(bounds_type)         , intent(in) :: bounds 
    !
    ! !LOCAL VARIABLES:
    integer           :: k,l,ii,jj 
    integer           :: begp,endp
    integer           :: begc,endc
    integer           :: begg,endg
    character(24)     :: fieldname
    character(100)    :: longname
    real(r8), pointer :: data1dptr(:)   ! temp. pointer for slicing larger arrays
    !---------------------------------------------------------------------

    begp = bounds%begp; endp = bounds%endp
    begc = bounds%begc; endc = bounds%endc
    begg = bounds%begg; endg = bounds%endg

    !-------------------------------
    ! patch state variables 
    !-------------------------------
    
    if (use_crop) then
       this%grainn_patch(begp:endp) = spval
       call hist_addfld1d (fname='GRAINN', units='gN/m^2', &
            avgflag='A', long_name='grain N', &
            ptr_patch=this%grainn_patch)
       call hist_addfld1d (fname='CROPSEEDN_DEFICIT', units='gN/m^2', &
            avgflag='A', long_name='N used for crop seed that needs to be repaid', &
            ptr_patch=this%cropseedn_deficit_patch)
    end if

    this%leafn_patch(begp:endp) = spval
    call hist_addfld1d (fname='LEAFN', units='gN/m^2', &
         avgflag='A', long_name='leaf N', &
         ptr_patch=this%leafn_patch)

    this%leafn_storage_patch(begp:endp) = spval
    call hist_addfld1d (fname='LEAFN_STORAGE', units='gN/m^2', &
         avgflag='A', long_name='leaf N storage', &
         ptr_patch=this%leafn_storage_patch)     

    this%leafn_xfer_patch(begp:endp) = spval
    call hist_addfld1d (fname='LEAFN_XFER', units='gN/m^2', &
         avgflag='A', long_name='leaf N transfer', &
         ptr_patch=this%leafn_xfer_patch)     

    if ( use_fun ) then
       this%leafn_storage_xfer_acc_patch(begp:endp) = spval
       call hist_addfld1d (fname='LEAFN_STORAGE_XFER_ACC', units='gN/m^2', &
            avgflag='A', long_name='Accmulated leaf N transfer', &
            ptr_patch=this%leafn_storage_xfer_acc_patch, default='inactive')

       this%storage_ndemand_patch(begp:endp)        = spval
       call hist_addfld1d (fname='STORAGE_NDEMAND', units='gN/m^2', &
            avgflag='A', long_name='N demand during the offset period', &
            ptr_patch=this%storage_ndemand_patch, default='inactive')
    end if

    this%frootn_patch(begp:endp) = spval
    call hist_addfld1d (fname='FROOTN', units='gN/m^2', &
         avgflag='A', long_name='fine root N', &
         ptr_patch=this%frootn_patch)

    this%frootn_storage_patch(begp:endp) = spval
    call hist_addfld1d (fname='FROOTN_STORAGE', units='gN/m^2', &
         avgflag='A', long_name='fine root N storage', &
         ptr_patch=this%frootn_storage_patch)     

    this%frootn_xfer_patch(begp:endp) = spval
    call hist_addfld1d (fname='FROOTN_XFER', units='gN/m^2', &
         avgflag='A', long_name='fine root N transfer', &
         ptr_patch=this%frootn_xfer_patch)     

    this%livestemn_patch(begp:endp) = spval
    call hist_addfld1d (fname='LIVESTEMN', units='gN/m^2', &
         avgflag='A', long_name='live stem N', &
         ptr_patch=this%livestemn_patch)

    this%livestemn_storage_patch(begp:endp) = spval
    call hist_addfld1d (fname='LIVESTEMN_STORAGE', units='gN/m^2', &
         avgflag='A', long_name='live stem N storage', &
         ptr_patch=this%livestemn_storage_patch)    

    this%livestemn_xfer_patch(begp:endp) = spval
    call hist_addfld1d (fname='LIVESTEMN_XFER', units='gN/m^2', &
         avgflag='A', long_name='live stem N transfer', &
         ptr_patch=this%livestemn_xfer_patch)     

    this%deadstemn_patch(begp:endp) = spval
    call hist_addfld1d (fname='DEADSTEMN', units='gN/m^2', &
         avgflag='A', long_name='dead stem N', &
         ptr_patch=this%deadstemn_patch)

    this%deadstemn_storage_patch(begp:endp) = spval
    call hist_addfld1d (fname='DEADSTEMN_STORAGE', units='gN/m^2', &
         avgflag='A', long_name='dead stem N storage', &
         ptr_patch=this%deadstemn_storage_patch)    

    this%deadstemn_xfer_patch(begp:endp) = spval
    call hist_addfld1d (fname='DEADSTEMN_XFER', units='gN/m^2', &
         avgflag='A', long_name='dead stem N transfer', &
         ptr_patch=this%deadstemn_xfer_patch)    

    this%livecrootn_patch(begp:endp) = spval
    call hist_addfld1d (fname='LIVECROOTN', units='gN/m^2', &
         avgflag='A', long_name='live coarse root N', &
         ptr_patch=this%livecrootn_patch)

    this%livecrootn_storage_patch(begp:endp) = spval
    call hist_addfld1d (fname='LIVECROOTN_STORAGE', units='gN/m^2', &
         avgflag='A', long_name='live coarse root N storage', &
         ptr_patch=this%livecrootn_storage_patch)    

    this%livecrootn_xfer_patch(begp:endp) = spval
    call hist_addfld1d (fname='LIVECROOTN_XFER', units='gN/m^2', &
         avgflag='A', long_name='live coarse root N transfer', &
         ptr_patch=this%livecrootn_xfer_patch)    

    this%deadcrootn_patch(begp:endp) = spval
    call hist_addfld1d (fname='DEADCROOTN', units='gN/m^2', &
         avgflag='A', long_name='dead coarse root N', &
         ptr_patch=this%deadcrootn_patch)

    this%deadcrootn_storage_patch(begp:endp) = spval
    call hist_addfld1d (fname='DEADCROOTN_STORAGE', units='gN/m^2', &
         avgflag='A', long_name='dead coarse root N storage', &
         ptr_patch=this%deadcrootn_storage_patch)    

    this%deadcrootn_xfer_patch(begp:endp) = spval
    call hist_addfld1d (fname='DEADCROOTN_XFER', units='gN/m^2', &
         avgflag='A', long_name='dead coarse root N transfer', &
         ptr_patch=this%deadcrootn_xfer_patch)    

    this%retransn_patch(begp:endp) = spval
    call hist_addfld1d (fname='RETRANSN', units='gN/m^2', &
         avgflag='A', long_name='plant pool of retranslocated N', &
         ptr_patch=this%retransn_patch)

    this%npool_patch(begp:endp) = spval
    call hist_addfld1d (fname='NPOOL', units='gN/m^2', &
         avgflag='A', long_name='temporary plant N pool', &
         ptr_patch=this%npool_patch)     

    this%ntrunc_patch(begp:endp) = spval
    call hist_addfld1d (fname='PFT_NTRUNC', units='gN/m^2', &
         avgflag='A', long_name='patch-level sink for N truncation', &
         ptr_patch=this%ntrunc_patch)

    this%dispvegn_patch(begp:endp) = spval
    call hist_addfld1d (fname='DISPVEGN', units='gN/m^2', &
         avgflag='A', long_name='displayed vegetation nitrogen', &
         ptr_patch=this%dispvegn_patch)

    this%storvegn_patch(begp:endp) = spval
    call hist_addfld1d (fname='STORVEGN', units='gN/m^2', &
         avgflag='A', long_name='stored vegetation nitrogen', &
         ptr_patch=this%storvegn_patch)

    this%totvegn_patch(begp:endp) = spval
    call hist_addfld1d (fname='TOTVEGN', units='gN/m^2', &
         avgflag='A', long_name='total vegetation nitrogen', &
         ptr_patch=this%totvegn_patch)

    this%totn_patch(begp:endp) = spval
    call hist_addfld1d (fname='TOTPFTN', units='gN/m^2', &
         avgflag='A', long_name='total patch-level nitrogen', &
         ptr_patch=this%totn_patch)

    !-------------------------------
    ! column state variables 
    !-------------------------------

    this%seedn_grc(begg:endg) = spval
    call hist_addfld1d (fname='SEEDN', units='gN/m^2', &
         avgflag='A', long_name='pool for seeding new PFTs via dynamic landcover', &
         ptr_gcell=this%seedn_grc)

    this%totecosysn_col(begc:endc) = spval
    call hist_addfld1d (fname='TOTECOSYSN', units='gN/m^2', &
         avgflag='A', long_name='total ecosystem N, excluding product pools', &
         ptr_col=this%totecosysn_col)

    this%totn_col(begc:endc) = spval
    call hist_addfld1d (fname='TOTCOLN', units='gN/m^2', &
         avgflag='A', long_name='total column-level N, excluding product pools', &
         ptr_col=this%totn_col)

  end subroutine InitHistory

  !-----------------------------------------------------------------------
  subroutine InitCold(this, bounds, &
       leafc_patch, leafc_storage_patch, frootc_patch, frootc_storage_patch, deadstemc_patch)
    !
    ! !DESCRIPTION:
    ! Initializes time varying variables used only in coupled carbon-nitrogen mode (CN):
    !
    use clm_varctl     , only : MM_Nuptake_opt   
    ! !ARGUMENTS:
    class(cnveg_nitrogenstate_type) :: this
    type(bounds_type) , intent(in) :: bounds  
    real(r8)          , intent(in) :: leafc_patch(bounds%begp:)
    real(r8)          , intent(in) :: leafc_storage_patch(bounds%begp:)
    real(r8)          , intent(in) :: frootc_patch(bounds%begp:)            
    real(r8)          , intent(in) :: frootc_storage_patch(bounds%begp:)    
    real(r8)          , intent(in) :: deadstemc_patch(bounds%begp:)
    !
    ! !LOCAL VARIABLES:
    integer :: fc,fp,g,l,c,p,j,k                       ! indices
    integer :: num_special_col                         ! number of good values in special_col filter
    integer :: num_special_patch                         ! number of good values in special_patch filter
    integer :: special_col   (bounds%endc-bounds%begc+1) ! special landunit filter - columns
    integer :: special_patch (bounds%endp-bounds%begp+1) ! special landunit filter - patches
    !------------------------------------------------------------------------

    SHR_ASSERT_ALL((ubound(leafc_patch)          == (/bounds%endp/)), errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(leafc_storage_patch)  == (/bounds%endp/)), errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(frootc_patch)         == (/bounds%endp/)), errMsg(sourcefile, __LINE__))   
    SHR_ASSERT_ALL((ubound(frootc_storage_patch) == (/bounds%endp/)), errMsg(sourcefile, __LINE__))   
    SHR_ASSERT_ALL((ubound(deadstemc_patch)      == (/bounds%endp/)), errMsg(sourcefile, __LINE__))

    ! Set column filters

    num_special_patch = 0
    do p = bounds%begp,bounds%endp
       l = patch%landunit(p)
       if (lun%ifspecial(l)) then
          num_special_patch = num_special_patch + 1
          special_patch(num_special_patch) = p
       end if
    end do

    ! Set patch filters

    num_special_col = 0
    do c = bounds%begc, bounds%endc
       l = col%landunit(c)
       if (lun%ifspecial(l)) then
          num_special_col = num_special_col + 1
          special_col(num_special_col) = c
       end if
    end do

    !-------------------------------------------
    ! initialize patch-level variables
    !-------------------------------------------

    do p = bounds%begp,bounds%endp

       l = patch%landunit(p)
       if (lun%itype(l) == istsoil .or. lun%itype(l) == istcrop) then

          if (patch%itype(p) == noveg) then
             this%leafn_patch(p) = 0._r8
             this%leafn_storage_patch(p) = 0._r8
             if (MM_Nuptake_opt .eqv. .true.) then   
                this%frootn_patch(p) = 0._r8            
                this%frootn_storage_patch(p) = 0._r8    
             end if 
          else
             this%leafn_patch(p)         = leafc_patch(p)         / pftcon%leafcn(patch%itype(p))
             this%leafn_storage_patch(p) = leafc_storage_patch(p) / pftcon%leafcn(patch%itype(p))
             if (MM_Nuptake_opt .eqv. .true.) then  
                this%frootn_patch(p) = frootc_patch(p) / pftcon%frootcn(patch%itype(p))           
                this%frootn_storage_patch(p) = frootc_storage_patch(p) / pftcon%frootcn(patch%itype(p))   
             end if 
          end if

          this%leafn_xfer_patch(p)        = 0._r8

          this%leafn_storage_xfer_acc_patch(p)        = 0._r8
          this%storage_ndemand_patch(p)   = 0._r8

          if ( use_crop )then
             this%grainn_patch(p)         = 0._r8
             this%grainn_storage_patch(p) = 0._r8
             this%grainn_xfer_patch(p)    = 0._r8
          end if
          if (MM_Nuptake_opt .eqv. .false.) then  ! if not running in floating CN ratio option 
             this%frootn_patch(p)            = 0._r8
             this%frootn_storage_patch(p)    = 0._r8
          end if 
          this%frootn_xfer_patch(p)       = 0._r8
          this%livestemn_patch(p)         = 0._r8
          this%livestemn_storage_patch(p) = 0._r8
          this%livestemn_xfer_patch(p)    = 0._r8

          ! tree types need to be initialized with some stem mass so that
          ! roughness length is not zero in canopy flux calculation

          if (pftcon%woody(patch%itype(p)) == 1._r8) then
             this%deadstemn_patch(p) = deadstemc_patch(p) / pftcon%deadwdcn(patch%itype(p))
          else
             this%deadstemn_patch(p) = 0._r8
          end if

          this%deadstemn_storage_patch(p)  = 0._r8
          this%deadstemn_xfer_patch(p)     = 0._r8
          this%livecrootn_patch(p)         = 0._r8
          this%livecrootn_storage_patch(p) = 0._r8
          this%livecrootn_xfer_patch(p)    = 0._r8
          this%deadcrootn_patch(p)         = 0._r8
          this%deadcrootn_storage_patch(p) = 0._r8
          this%deadcrootn_xfer_patch(p)    = 0._r8
          this%retransn_patch(p)           = 0._r8
          this%npool_patch(p)              = 0._r8
          this%ntrunc_patch(p)             = 0._r8
          this%dispvegn_patch(p)           = 0._r8
          this%storvegn_patch(p)           = 0._r8
          this%totvegn_patch(p)            = 0._r8
          this%totn_patch(p)               = 0._r8
          this%cropseedn_deficit_patch(p)  = 0._r8
       end if
    end do

    !-------------------------------------------
    ! initialize column-level variables
    !-------------------------------------------

    do c = bounds%begc, bounds%endc
       l = col%landunit(c)
       if (lun%itype(l) == istsoil .or. lun%itype(l) == istcrop) then
          ! total nitrogen pools
          this%totecosysn_col(c) = 0._r8
          this%totn_p2c_col(c)   = 0._r8
          this%totn_col(c)       = 0._r8
       end if
    end do


    do g = bounds%begg, bounds%endg
       this%seedn_grc(g) = 0._r8
    end do

    ! now loop through special filters and explicitly set the variables that
    ! have to be in place for biogeophysics

    ! initialize fields for special filters

    call this%SetValues (&
         num_patch=num_special_patch, filter_patch=special_patch, value_patch=0._r8, &
         num_column=num_special_col, filter_column=special_col, value_column=0._r8)

  end subroutine InitCold

  !-----------------------------------------------------------------------
  subroutine Restart ( this,  bounds, ncid, flag )
    !
    ! !DESCRIPTION: 
    ! Read/write restart data 
    !
    ! !USES:
    use restUtilMod
    use ncdio_pio
    use clm_varctl             , only : spinup_state, use_cndv
    use clm_time_manager       , only : get_nstep

    !
    ! !ARGUMENTS:
    class (cnveg_nitrogenstate_type) :: this
    type(bounds_type)          , intent(in)    :: bounds 
    type(file_desc_t)          , intent(inout) :: ncid   
    character(len=*)           , intent(in)    :: flag   !'read' or 'write' or 'define'
    !
    ! !LOCAL VARIABLES:
    integer            :: i,j,k,l,c
    logical            :: readvar
    real(r8), pointer  :: ptr1d(:)   ! temp. pointers for slicing larger arrays
    character(len=128) :: varname    ! temporary
    logical            :: exit_spinup  = .false.
    logical            :: enter_spinup = .false.
    integer            :: idata

    ! spinup state as read from restart file, for determining whether to enter or exit spinup mode.
    integer            :: restart_file_spinup_state

    !------------------------------------------------------------------------

    !--------------------------------
    ! patch nitrogen state variables
    !--------------------------------

    call restartvar(ncid=ncid, flag=flag, varname='leafn', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%leafn_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='leafn_storage', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%leafn_storage_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='leafn_xfer', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%leafn_xfer_patch) 

     if ( use_fun ) then
        call restartvar(ncid=ncid, flag=flag, varname='leafn_storage_xfer_acc', xtype=ncd_double,  &
             dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%leafn_storage_xfer_acc_patch)
    
        call restartvar(ncid=ncid, flag=flag, varname='storage_ndemand', xtype=ncd_double,  &
             dim1name='pft', long_name='', units='', &
             interpinic_flag='interp', readvar=readvar, data=this%storage_ndemand_patch)
     end if


    call restartvar(ncid=ncid, flag=flag, varname='frootn', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%frootn_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='frootn_storage', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%frootn_storage_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='frootn_xfer', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%frootn_xfer_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='livestemn', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%livestemn_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='livestemn_storage', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%livestemn_storage_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='livestemn_xfer', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%livestemn_xfer_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='deadstemn', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%deadstemn_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='deadstemn_storage', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%deadstemn_storage_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='deadstemn_xfer', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%deadstemn_xfer_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='livecrootn', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%livecrootn_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='livecrootn_storage', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%livecrootn_storage_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='livecrootn_xfer', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%livecrootn_xfer_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='deadcrootn', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%deadcrootn_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='deadcrootn_storage', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%deadcrootn_storage_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='deadcrootn_xfer', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%deadcrootn_xfer_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='retransn', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%retransn_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='npool', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%npool_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='pft_ntrunc', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%ntrunc_patch) 

    if (use_crop) then
       call restartvar(ncid=ncid, flag=flag,  varname='grainn', xtype=ncd_double,  &
            dim1name='pft',    long_name='grain N', units='gN/m2', &
            interpinic_flag='interp', readvar=readvar, data=this%grainn_patch)

       call restartvar(ncid=ncid, flag=flag,  varname='grainn_storage', xtype=ncd_double,  &
            dim1name='pft',    long_name='grain N storage', units='gN/m2', &
            interpinic_flag='interp', readvar=readvar, data=this%grainn_storage_patch)

       call restartvar(ncid=ncid, flag=flag,  varname='grainn_xfer', xtype=ncd_double,  &
            dim1name='pft',    long_name='grain N transfer', units='gN/m2', &
            interpinic_flag='interp', readvar=readvar, data=this%grainn_xfer_patch)

       call restartvar(ncid=ncid, flag=flag, varname='cropseedn_deficit', xtype=ncd_double,  &
            dim1name='pft', long_name='pool for seeding new crop growth', units='gN/m2', &
            interpinic_flag='interp', readvar=readvar, data=this%cropseedn_deficit_patch)
    end if

    !--------------------------------
    ! gridcell nitrogen state variables
    !--------------------------------

    ! BACKWARDS_COMPATIBILITY(wjs, 2017-01-12) Naming this with a _g suffix in order to
    ! distinguish it from the old column-level seedn restart variable
    call restartvar(ncid=ncid, flag=flag, varname='seedn_g', xtype=ncd_double,  &
         dim1name='gridcell', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%seedn_grc) 


    if (flag == 'read') then
       call restartvar(ncid=ncid, flag=flag, varname='spinup_state', xtype=ncd_int, &
         long_name='Spinup state of the model that wrote this restart file: ' &
         // ' 0 = normal model mode, 1 = AD spinup', units='', &
         interpinic_flag='copy', readvar=readvar,  data=idata)

       if (readvar) then
          restart_file_spinup_state = idata
       else
          restart_file_spinup_state = spinup_state
          if ( masterproc ) then
             write(iulog,*) ' CNRest: WARNING!  Restart file does not contain info ' &
                   // ' on spinup state used to generate the restart file. '
             write(iulog,*) '   Assuming the same as current setting: ', spinup_state
          end if
       end if
    end if

    if (flag == 'read' .and. spinup_state /= restart_file_spinup_state .and. .not. use_cndv) then
       if (spinup_state <= 1 .and. restart_file_spinup_state == 2 ) then
          if ( masterproc ) write(iulog,*) ' CNRest: taking Dead wood N pools out of AD spinup mode'
          exit_spinup = .true.
          if ( masterproc ) write(iulog, *) 'Multiplying stemn and crootn by 10 for exit spinup '
          do i = bounds%begp,bounds%endp
             this%deadstemn_patch(i) = this%deadstemn_patch(i) * 10._r8
             this%deadcrootn_patch(i) = this%deadcrootn_patch(i) * 10._r8
          end do
       else if (spinup_state == 2 .and. restart_file_spinup_state <= 1 ) then
          if ( masterproc ) write(iulog,*) ' CNRest: taking Dead wood N pools into AD spinup mode'
          enter_spinup = .true.
          if ( masterproc ) write(iulog, *) 'Dividing stemn and crootn by 10 for enter spinup '
          do i = bounds%begp,bounds%endp
             this%deadstemn_patch(i) = this%deadstemn_patch(i) / 10._r8
             this%deadcrootn_patch(i) = this%deadcrootn_patch(i) / 10._r8
          end do
       endif

    end if

  end subroutine Restart

  !-----------------------------------------------------------------------
  subroutine SetValues ( this, &
       num_patch, filter_patch, value_patch, &
       num_column, filter_column, value_column)
    !
    ! !DESCRIPTION:
    ! Set nitrogen state variables
    !
    ! !ARGUMENTS:
    class (cnveg_nitrogenstate_type) :: this
    integer , intent(in) :: num_patch
    integer , intent(in) :: filter_patch(:)
    real(r8), intent(in) :: value_patch
    integer , intent(in) :: num_column
    integer , intent(in) :: filter_column(:)
    real(r8), intent(in) :: value_column
    !
    ! !LOCAL VARIABLES:
    integer :: fi,i     ! loop index
    integer :: j,k      ! indices
    !------------------------------------------------------------------------

    do fi = 1,num_patch
       i = filter_patch(fi)

       this%leafn_patch(i)              = value_patch
       this%leafn_storage_patch(i)      = value_patch
       this%leafn_xfer_patch(i)         = value_patch
       this%leafn_storage_xfer_acc_patch(i) = value_patch
       this%frootn_patch(i)             = value_patch
       this%frootn_storage_patch(i)     = value_patch
       this%frootn_xfer_patch(i)        = value_patch
       this%livestemn_patch(i)          = value_patch
       this%livestemn_storage_patch(i)  = value_patch
       this%livestemn_xfer_patch(i)     = value_patch
       this%deadstemn_patch(i)          = value_patch
       this%deadstemn_storage_patch(i)  = value_patch
       this%deadstemn_xfer_patch(i)     = value_patch
       this%livecrootn_patch(i)         = value_patch
       this%livecrootn_storage_patch(i) = value_patch
       this%livecrootn_xfer_patch(i)    = value_patch
       this%deadcrootn_patch(i)         = value_patch
       this%deadcrootn_storage_patch(i) = value_patch
       this%deadcrootn_xfer_patch(i)    = value_patch
       this%retransn_patch(i)           = value_patch
       this%npool_patch(i)              = value_patch
       this%ntrunc_patch(i)             = value_patch
       this%dispvegn_patch(i)           = value_patch
       this%storvegn_patch(i)           = value_patch
       this%totvegn_patch(i)            = value_patch
       this%totn_patch(i)               = value_patch
       this%cropseedn_deficit_patch(i)  = value_patch
    end do

    if ( use_crop )then
       do fi = 1,num_patch
          i = filter_patch(fi)
          this%grainn_patch(i)          = value_patch
          this%grainn_storage_patch(i)  = value_patch
          this%grainn_xfer_patch(i)     = value_patch   
       end do
    end if

    do fi = 1,num_column
       i = filter_column(fi)

       this%totecosysn_col(i) = value_column
       this%totvegn_col(i)    = value_column
       this%totn_p2c_col(i)   = value_column
       this%totn_col(i)       = value_column
    end do

  end subroutine SetValues

  !-----------------------------------------------------------------------
  subroutine ZeroDwt( this, bounds )
    !
    ! !DESCRIPTION
    ! Initialize variables needed for dynamic land use.
    !
    ! !ARGUMENTS:
    class(cnveg_nitrogenstate_type) :: this
    type(bounds_type), intent(in)  :: bounds 
    !
    ! !LOCAL VARIABLES:
    integer  :: p          ! indices
    !-----------------------------------------------------------------------

    do p = bounds%begp,bounds%endp
       this%dispvegn_patch(p) = 0._r8
       this%storvegn_patch(p) = 0._r8
       this%totvegn_patch(p)  = 0._r8
       this%totn_patch(p)     = 0._r8
    end do

  end subroutine ZeroDwt

  !-----------------------------------------------------------------------
  subroutine Summary_nitrogenstate(this, bounds, num_allc, filter_allc, &
       num_soilc, filter_soilc, num_soilp, filter_soilp,&
       soilbiogeochem_nitrogenstate_inst)
    !
    ! !USES:
    use subgridAveMod, only : p2c
    use SoilBiogeochemNitrogenStateType, only : soilbiogeochem_nitrogenstate_type
    !
    ! !ARGUMENTS:
    class(cnveg_nitrogenstate_type)                      :: this
    type(bounds_type)                       , intent(in) :: bounds  
    integer                                 , intent(in) :: num_allc        ! number of columns in allc filter
    integer                                 , intent(in) :: filter_allc(:)  ! filter for all active columns
    integer                                 , intent(in) :: num_soilc       ! number of soil columns in filter
    integer                                 , intent(in) :: filter_soilc(:) ! filter for soil columns
    integer                                 , intent(in) :: num_soilp       ! number of soil patches in filter
    integer                                 , intent(in) :: filter_soilp(:) ! filter for soil patches
    type(soilbiogeochem_nitrogenstate_type) , intent(in) :: soilbiogeochem_nitrogenstate_inst
    !
    ! !LOCAL VARIABLES:
    integer  :: c,p,j,k,l ! indices
    integer  :: fp,fc       ! lake filter indices
    real(r8) :: maxdepth    ! depth to integrate soil variables
    !-----------------------------------------------------------------------

    ! --------------------------------------------
    ! patch level summary
    ! --------------------------------------------
    
    do fp = 1,num_soilp
       p = filter_soilp(fp)
         
	      
       ! displayed vegetation nitrogen, excluding storage (DISPVEGN)
       this%dispvegn_patch(p) = &
            this%leafn_patch(p)      + &
            this%frootn_patch(p)     + &
            this%livestemn_patch(p)  + &
            this%deadstemn_patch(p)  + &
            this%livecrootn_patch(p) + &
            this%deadcrootn_patch(p)

       ! stored vegetation nitrogen, including retranslocated N pool (STORVEGN)
       this%storvegn_patch(p) = &
            this%leafn_storage_patch(p)      + &
            this%frootn_storage_patch(p)     + &
            this%livestemn_storage_patch(p)  + &
            this%deadstemn_storage_patch(p)  + &
            this%livecrootn_storage_patch(p) + &
            this%deadcrootn_storage_patch(p) + &
            this%leafn_xfer_patch(p)         + &
            this%frootn_xfer_patch(p)        + &
            this%livestemn_xfer_patch(p)     + &
            this%deadstemn_xfer_patch(p)     + &
            this%livecrootn_xfer_patch(p)    + &
            this%deadcrootn_xfer_patch(p)    + &
            this%npool_patch(p)              + &
            this%retransn_patch(p)

       if ( use_crop .and. patch%itype(p) >= npcropmin )then
          this%dispvegn_patch(p) = &
               this%dispvegn_patch(p) + &
               this%grainn_patch(p)

          this%storvegn_patch(p) = &
               this%storvegn_patch(p) + &
               this%grainn_storage_patch(p)     + &
               this%grainn_xfer_patch(p) + &
               this%cropseedn_deficit_patch(p)
       end if

       ! total vegetation nitrogen (TOTVEGN)
       this%totvegn_patch(p) = &
            this%dispvegn_patch(p) + &
            this%storvegn_patch(p)

       ! total patch-level carbon (add ntrunc)
       this%totn_patch(p) = &
            this%totvegn_patch(p) + &
            this%ntrunc_patch(p)
            
    end do

    ! --------------------------------------------
    ! column level summary
    ! --------------------------------------------

    call p2c(bounds, num_soilc, filter_soilc, &
         this%totvegn_patch(bounds%begp:bounds%endp), &
         this%totvegn_col(bounds%begc:bounds%endc))

    call p2c(bounds, num_soilc, filter_soilc, &
         this%totn_patch(bounds%begp:bounds%endp), &
         this%totn_p2c_col(bounds%begc:bounds%endc))

    do fc = 1,num_allc
       c = filter_allc(fc)

       ! total ecosystem nitrogen, including veg (TOTECOSYSN)
       this%totecosysn_col(c) =    &
            soilbiogeochem_nitrogenstate_inst%cwdn_col(c)    + &
            soilbiogeochem_nitrogenstate_inst%totlitn_col(c) + &
            soilbiogeochem_nitrogenstate_inst%totsomn_col(c) + &
            soilbiogeochem_nitrogenstate_inst%sminn_col(c)   + &
            this%totvegn_col(c)                              

       ! total column nitrogen, including patch (TOTCOLN)

       this%totn_col(c) = this%totn_p2c_col(c)               + &
            soilbiogeochem_nitrogenstate_inst%cwdn_col(c)    + &
            soilbiogeochem_nitrogenstate_inst%totlitn_col(c) + &
            soilbiogeochem_nitrogenstate_inst%totsomn_col(c) + &
            soilbiogeochem_nitrogenstate_inst%sminn_col(c)   + &
            soilbiogeochem_nitrogenstate_inst%ntrunc_col(c)

    end do
    
    
    

  end subroutine Summary_nitrogenstate

  !-----------------------------------------------------------------------
  subroutine DynamicPatchAdjustments(this, bounds, &
       num_soilp_with_inactive, filter_soilp_with_inactive, &
       patch_state_updater, &
       leafc_seed, deadstemc_seed, &
       conv_nflux, product_nflux, &
       dwt_frootn_to_litter, &
       dwt_livecrootn_to_litter, &
       dwt_deadcrootn_to_litter, &
       dwt_leafn_seed, &
       dwt_deadstemn_seed)
    !
    ! !DESCRIPTION:
    ! Adjust state variables and compute associated fluxes when patch areas change due to
    ! dynamic landuse
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(cnveg_nitrogenstate_type) , intent(inout) :: this
    type(bounds_type)               , intent(in)    :: bounds
    integer                         , intent(in)    :: num_soilp_with_inactive ! number of points in filter
    integer                         , intent(in)    :: filter_soilp_with_inactive(:) ! soil patch filter that includes inactive points
    type(patch_state_updater_type)  , intent(in)    :: patch_state_updater
    real(r8)                        , intent(in)    :: leafc_seed  ! seed amount for leaf C
    real(r8)                        , intent(in)    :: deadstemc_seed ! seed amount for deadstem C
    real(r8)                        , intent(inout) :: conv_nflux( bounds%begp: )  ! patch-level conversion N flux to atm
    real(r8)                        , intent(inout) :: product_nflux( bounds%begp: ) ! patch-level product N flux
    real(r8)                        , intent(inout) :: dwt_frootn_to_litter( bounds%begp: ) ! patch-level fine root N to litter
    real(r8)                        , intent(inout) :: dwt_livecrootn_to_litter( bounds%begp: ) ! patch-level live coarse root N to litter
    real(r8)                        , intent(inout) :: dwt_deadcrootn_to_litter( bounds%begp: ) ! patch-level live coarse root N to litter
    real(r8)                        , intent(inout) :: dwt_leafn_seed( bounds%begp: ) ! patch-level mass gain due to seeding of new area: leaf N
    real(r8)                        , intent(inout) :: dwt_deadstemn_seed( bounds%begp: ) ! patch-level mass gain due to seeding of new area: deadstem N
    !
    ! !LOCAL VARIABLES:
    integer :: begp, endp

    logical  :: old_weight_was_zero(bounds%begp:bounds%endp)
    logical  :: patch_grew(bounds%begp:bounds%endp)

    ! The following are only set for growing patches:
    real(r8) :: seed_leafn_patch(bounds%begp:bounds%endp)
    real(r8) :: seed_leafn_storage_patch(bounds%begp:bounds%endp)
    real(r8) :: seed_leafn_xfer_patch(bounds%begp:bounds%endp)
    real(r8) :: seed_deadstemn_patch(bounds%begp:bounds%endp)

    character(len=*), parameter :: subname = 'DynamicPatchAdjustments'
    !-----------------------------------------------------------------------

    begp = bounds%begp
    endp = bounds%endp

    SHR_ASSERT_ALL((ubound(conv_nflux) == (/endp/)), errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(product_nflux) == (/endp/)), errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(dwt_frootn_to_litter) == (/endp/)), errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(dwt_livecrootn_to_litter) == (/endp/)), errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(dwt_deadcrootn_to_litter) == (/endp/)), errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(dwt_leafn_seed) == (/endp/)), errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(dwt_deadstemn_seed) == (/endp/)), errMsg(sourcefile, __LINE__))

    old_weight_was_zero = patch_state_updater%old_weight_was_zero(bounds)
    patch_grew = patch_state_updater%patch_grew(bounds)

    call ComputeSeedAmounts(bounds, &
         num_soilp_with_inactive, filter_soilp_with_inactive, &
         species = CN_SPECIES_N, &
         leafc_seed = leafc_seed, &
         deadstemc_seed = deadstemc_seed, &
         leaf_patch = this%leafn_patch(begp:endp), &
         leaf_storage_patch = this%leafn_storage_patch(begp:endp), &
         leaf_xfer_patch = this%leafn_xfer_patch(begp:endp), &

         ! Calculations only needed for patches that grew:
         compute_here_patch = patch_grew(begp:endp), &

         ! For patches that previously had zero area, ignore the current state for the
         ! sake of computing leaf proportions:
         ignore_current_state_patch = old_weight_was_zero(begp:endp), &

         seed_leaf_patch = seed_leafn_patch(begp:endp), &
         seed_leaf_storage_patch = seed_leafn_storage_patch(begp:endp), &
         seed_leaf_xfer_patch = seed_leafn_xfer_patch(begp:endp), &
         seed_deadstem_patch = seed_deadstemn_patch(begp:endp))

    call update_patch_state( &
         var = this%leafn_patch(begp:endp), &
         flux_out = conv_nflux(begp:endp), &
         seed = seed_leafn_patch(begp:endp), &
         seed_addition = dwt_leafn_seed(begp:endp))

    call update_patch_state( &
         var = this%leafn_storage_patch(begp:endp), &
         flux_out = conv_nflux(begp:endp), &
         seed = seed_leafn_storage_patch(begp:endp), &
         seed_addition = dwt_leafn_seed(begp:endp))

    call update_patch_state( &
         var = this%leafn_xfer_patch(begp:endp), &
         flux_out = conv_nflux(begp:endp), &
         seed = seed_leafn_xfer_patch(begp:endp), &
         seed_addition = dwt_leafn_seed(begp:endp))

    call update_patch_state( &
         var = this%frootn_patch(begp:endp), &
         flux_out = dwt_frootn_to_litter(begp:endp))

    call update_patch_state( &
         var = this%frootn_storage_patch(begp:endp), &
         flux_out = conv_nflux(begp:endp))

    call update_patch_state( &
         var = this%frootn_xfer_patch(begp:endp), &
         flux_out = conv_nflux(begp:endp))

    call update_patch_state( &
         var = this%livestemn_patch(begp:endp), &
         flux_out = conv_nflux(begp:endp))

    call update_patch_state( &
         var = this%livestemn_storage_patch(begp:endp), &
         flux_out = conv_nflux(begp:endp))

    call update_patch_state( &
         var = this%livestemn_xfer_patch(begp:endp), &
         flux_out = conv_nflux(begp:endp))

    call patch_state_updater%update_patch_state_partition_flux_by_type(bounds, &
         num_soilp_with_inactive, filter_soilp_with_inactive, &
         flux1_fraction_by_pft_type = pftcon%pconv, &
         var = this%deadstemn_patch(begp:endp), &
         flux1_out = conv_nflux(begp:endp), &
         flux2_out = product_nflux(begp:endp), &
         seed = seed_deadstemn_patch(begp:endp), &
         seed_addition = dwt_deadstemn_seed(begp:endp))

    call update_patch_state( &
         var = this%deadstemn_storage_patch(begp:endp), &
         flux_out = conv_nflux(begp:endp))

    call update_patch_state( &
         var = this%deadstemn_xfer_patch(begp:endp), &
         flux_out = conv_nflux(begp:endp))

    call update_patch_state( &
         var = this%livecrootn_patch(begp:endp), &
         flux_out = dwt_livecrootn_to_litter(begp:endp))

    call update_patch_state( &
         var = this%livecrootn_storage_patch(begp:endp), &
         flux_out = conv_nflux(begp:endp))

    call update_patch_state( &
         var = this%livecrootn_xfer_patch(begp:endp), &
         flux_out = conv_nflux(begp:endp))

    call update_patch_state( &
         var = this%deadcrootn_patch(begp:endp), &
         flux_out = dwt_deadcrootn_to_litter(begp:endp))

    call update_patch_state( &
         var = this%deadcrootn_storage_patch(begp:endp), &
         flux_out = conv_nflux(begp:endp))

    call update_patch_state( &
         var = this%deadcrootn_xfer_patch(begp:endp), &
         flux_out = conv_nflux(begp:endp))

    call update_patch_state( &
         var = this%retransn_patch(begp:endp), &
         flux_out = conv_nflux(begp:endp))

    call update_patch_state( &
         var = this%npool_patch(begp:endp), &
         flux_out = conv_nflux(begp:endp))

    call update_patch_state( &
         var = this%ntrunc_patch(begp:endp), &
         flux_out = conv_nflux(begp:endp))

  contains
    subroutine update_patch_state(var, flux_out, seed, seed_addition)
      ! Wraps call to update_patch_state, in order to remove duplication
      real(r8), intent(inout) :: var( bounds%begp: )
      real(r8), intent(inout), optional :: flux_out( bounds%begp: )
      real(r8), intent(in), optional :: seed( bounds%begp: )
      real(r8), intent(inout), optional :: seed_addition( bounds%begp: )
      
      call patch_state_updater%update_patch_state(bounds, &
         num_soilp_with_inactive, filter_soilp_with_inactive, &
         var = var, &
         flux_out = flux_out, &
         seed = seed, &
         seed_addition = seed_addition)
    end subroutine update_patch_state


  end subroutine DynamicPatchAdjustments

end module CNVegNitrogenStateType
