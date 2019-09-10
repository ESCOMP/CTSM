module CNVegCarbonStateType

#include "shr_assert.h"

  !-----------------------------------------------------------------------
  !
  ! !USES:
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
  use shr_const_mod  , only : SHR_CONST_PDB
  use shr_log_mod    , only : errMsg => shr_log_errMsg
  use pftconMod	     , only : noveg, npcropmin, pftcon
  use clm_varcon     , only : spval, c3_r2, c4_r2, c14ratio
  use clm_varctl     , only : iulog, use_cndv, use_crop
  use decompMod      , only : bounds_type
  use abortutils     , only : endrun
  use spmdMod        , only : masterproc 
  use LandunitType   , only : lun                
  use ColumnType     , only : col                
  use PatchType      , only : patch
  use CNSpeciesMod   , only : species_from_string, CN_SPECIES_C12
  use dynPatchStateUpdaterMod, only : patch_state_updater_type
  use CNVegComputeSeedMod, only : ComputeSeedAmounts
  ! 
  ! !PUBLIC TYPES:
  implicit none
  private
  !

  type, public :: cnveg_carbonstate_type

     integer :: species  ! c12, c13, c14

     real(r8), pointer :: grainc_patch             (:) ! (gC/m2) grain C (crop model)
     real(r8), pointer :: grainc_storage_patch     (:) ! (gC/m2) grain C storage (crop model)
     real(r8), pointer :: grainc_xfer_patch        (:) ! (gC/m2) grain C transfer (crop model)
     real(r8), pointer :: leafc_patch              (:) ! (gC/m2) leaf C
     real(r8), pointer :: leafc_storage_patch      (:) ! (gC/m2) leaf C storage
     real(r8), pointer :: leafc_xfer_patch         (:) ! (gC/m2) leaf C transfer
     real(r8), pointer :: leafc_storage_xfer_acc_patch   (:) ! (gC/m2) Accmulated leaf C transfer
     real(r8), pointer :: storage_cdemand_patch          (:) ! (gC/m2)       C use from the C storage pool 
     real(r8), pointer :: frootc_patch             (:) ! (gC/m2) fine root C
     real(r8), pointer :: frootc_storage_patch     (:) ! (gC/m2) fine root C storage
     real(r8), pointer :: frootc_xfer_patch        (:) ! (gC/m2) fine root C transfer
     real(r8), pointer :: livestemc_patch          (:) ! (gC/m2) live stem C
     real(r8), pointer :: livestemc_storage_patch  (:) ! (gC/m2) live stem C storage
     real(r8), pointer :: livestemc_xfer_patch     (:) ! (gC/m2) live stem C transfer
     real(r8), pointer :: deadstemc_patch          (:) ! (gC/m2) dead stem C
     real(r8), pointer :: deadstemc_storage_patch  (:) ! (gC/m2) dead stem C storage
     real(r8), pointer :: deadstemc_xfer_patch     (:) ! (gC/m2) dead stem C transfer
     real(r8), pointer :: livecrootc_patch         (:) ! (gC/m2) live coarse root C
     real(r8), pointer :: livecrootc_storage_patch (:) ! (gC/m2) live coarse root C storage
     real(r8), pointer :: livecrootc_xfer_patch    (:) ! (gC/m2) live coarse root C transfer
     real(r8), pointer :: deadcrootc_patch         (:) ! (gC/m2) dead coarse root C
     real(r8), pointer :: deadcrootc_storage_patch (:) ! (gC/m2) dead coarse root C storage
     real(r8), pointer :: deadcrootc_xfer_patch    (:) ! (gC/m2) dead coarse root C transfer
     real(r8), pointer :: gresp_storage_patch      (:) ! (gC/m2) growth respiration storage
     real(r8), pointer :: gresp_xfer_patch         (:) ! (gC/m2) growth respiration transfer
     real(r8), pointer :: cpool_patch              (:) ! (gC/m2) temporary photosynthate C pool
     real(r8), pointer :: xsmrpool_patch           (:) ! (gC/m2) abstract C pool to meet excess MR demand
     real(r8), pointer :: ctrunc_patch             (:) ! (gC/m2) patch-level sink for C truncation
     real(r8), pointer :: woodc_patch              (:) ! (gC/m2) wood C
     real(r8), pointer :: leafcmax_patch           (:) ! (gC/m2) ann max leaf C
     real(r8), pointer :: totc_patch               (:) ! (gC/m2) total patch-level carbon, including cpool
     real(r8), pointer :: rootc_col                (:) ! (gC/m2) root carbon at column level (fire)
     real(r8), pointer :: leafc_col                (:) ! (gC/m2) column-level leafc (fire)
     real(r8), pointer :: deadstemc_col            (:) ! (gC/m2) column-level deadstemc (fire)
     real(r8), pointer :: fuelc_col                (:) ! fuel load outside cropland
     real(r8), pointer :: fuelc_crop_col           (:) ! fuel load for cropland
     real(r8), pointer :: cropseedc_deficit_patch  (:) ! (gC/m2) pool for seeding new crop growth; this is a NEGATIVE term, indicating the amount of seed usage that needs to be repaid

     ! pools for dynamic landcover
     real(r8), pointer :: seedc_grc                (:) ! (gC/m2) gridcell-level pool for seeding new PFTs via dynamic landcover

     ! summary (diagnostic) state variables, not involved in mass balance
     real(r8), pointer :: dispvegc_patch           (:) ! (gC/m2) displayed veg carbon, excluding storage and cpool
     real(r8), pointer :: storvegc_patch           (:) ! (gC/m2) stored vegetation carbon, excluding cpool
     real(r8), pointer :: totvegc_patch            (:) ! (gC/m2) total vegetation carbon, excluding cpool
     real(r8), pointer :: totvegc_col              (:) ! (gC/m2) total vegetation carbon, excluding cpool averaged to column (p2c)

     ! Total C pools       
     real(r8), pointer :: totc_p2c_col             (:) ! (gC/m2) totc_patch averaged to col
     real(r8), pointer :: totc_col                 (:) ! (gC/m2) total column carbon, incl veg and cpool
     real(r8), pointer :: totecosysc_col           (:) ! (gC/m2) total ecosystem carbon, incl veg but excl cpool 

   contains

     procedure , public  :: Init   
     procedure , public  :: SetValues
     procedure , public  :: ZeroDWT
     procedure , public  :: Restart
     procedure , public  :: Summary => Summary_carbonstate
     procedure , public  :: DynamicPatchAdjustments   ! adjust state variables when patch areas change
     
     procedure , private :: InitAllocate    ! Allocate arrays
     procedure , private :: InitReadNML     ! Read in namelist
     procedure , private :: InitHistory     ! Initialize history
     procedure , private :: InitCold        ! Initialize arrays for a cold-start

  end type cnveg_carbonstate_type

  ! !PRIVATE DATA:

  type, private :: cnvegcarbonstate_const_type
      ! !PRIVATE MEMBER DATA:
      real(r8) :: initial_vegC = 20._r8    ! Initial vegetation carbon for leafc/frootc and storage
  end type
  type(cnvegcarbonstate_const_type), private :: cnvegcstate_const    ! Constants used here
  character(len=*), parameter :: sourcefile = &
       __FILE__

  !------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  subroutine Init(this, bounds, carbon_type, ratio, NLFilename, &
                  c12_cnveg_carbonstate_inst)

    class(cnveg_carbonstate_type)                       :: this
    type(bounds_type)            , intent(in)           :: bounds  
    real(r8)                     , intent(in)           :: ratio
    character(len=*)             , intent(in)           :: carbon_type                ! Carbon isotope type C12, C13 or C1
    character(len=*)             , intent(in)           :: NLFilename                 ! Namelist filename
    type(cnveg_carbonstate_type) , intent(in), optional :: c12_cnveg_carbonstate_inst ! cnveg_carbonstate for C12 (if C13 or C14)
    !-----------------------------------------------------------------------

    this%species = species_from_string(carbon_type)

    call this%InitAllocate ( bounds)
    call this%InitReadNML  ( NLFilename )
    call this%InitHistory ( bounds, carbon_type)
    if (present(c12_cnveg_carbonstate_inst)) then
       call this%InitCold  ( bounds, ratio, carbon_type, c12_cnveg_carbonstate_inst )
    else
       call this%InitCold  ( bounds, ratio, carbon_type )
    end if

  end subroutine Init

  !------------------------------------------------------------------------
  subroutine InitReadNML(this, NLFilename)
    !
    ! !DESCRIPTION:
    ! Read the namelist for CNVegCarbonState
    !
    !USES:
    use fileutils      , only : getavu, relavu, opnfil
    use shr_nl_mod     , only : shr_nl_find_group_name
    use spmdMod        , only : masterproc, mpicom
    use shr_mpi_mod    , only : shr_mpi_bcast
    use clm_varctl     , only : iulog
    !
    ! !ARGUMENTS:
    class(cnveg_carbonstate_type)                       :: this
    character(len=*)             , intent(in)           :: NLFilename                 ! Namelist filename
    !
    ! !LOCAL VARIABLES:
    integer :: ierr                 ! error code
    integer :: unitn                ! unit for namelist file

    character(len=*), parameter :: subname = 'InitReadNML'
    character(len=*), parameter :: nmlname = 'cnvegcarbonstate'   ! MUST match what is in namelist below
    !-----------------------------------------------------------------------
    real(r8) :: initial_vegC
    namelist /cnvegcarbonstate/ initial_vegC

    initial_vegC = cnvegcstate_const%initial_vegC

    if (masterproc) then
       unitn = getavu()
       write(iulog,*) 'Read in '//nmlname//'  namelist'
       call opnfil (NLFilename, unitn, 'F')
       call shr_nl_find_group_name(unitn, nmlname, status=ierr)
       if (ierr == 0) then
          read(unitn, nml=cnvegcarbonstate, iostat=ierr)
          if (ierr /= 0) then
             call endrun(msg="ERROR reading "//nmlname//"namelist"//errmsg(sourcefile, __LINE__))
          end if
       else
          call endrun(msg="ERROR could NOT find "//nmlname//"namelist"//errmsg(sourcefile, __LINE__))
       end if
       call relavu( unitn )
    end if

    call shr_mpi_bcast (initial_vegC            , mpicom)

    cnvegcstate_const%initial_vegC = initial_vegC

    if (masterproc) then
       write(iulog,*) ' '
       write(iulog,*) nmlname//' settings:'
       write(iulog,nml=cnvegcarbonstate)    ! Name here MUST be the same as in nmlname above!
       write(iulog,*) ' '
    end if

    !-----------------------------------------------------------------------

  end subroutine InitReadNML

  !------------------------------------------------------------------------
  subroutine InitAllocate(this, bounds)
    !
    ! !ARGUMENTS:
    class (cnveg_carbonstate_type) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer           :: begp,endp
    integer           :: begc,endc
    integer           :: begg,endg
    !------------------------------------------------------------------------

    begp = bounds%begp; endp = bounds%endp
    begc = bounds%begc; endc = bounds%endc
    begg = bounds%begg; endg = bounds%endg

    allocate(this%leafc_patch              (begp:endp)) ; this%leafc_patch              (:) = nan
    allocate(this%leafc_storage_patch      (begp:endp)) ; this%leafc_storage_patch      (:) = nan
    allocate(this%leafc_xfer_patch         (begp:endp)) ; this%leafc_xfer_patch         (:) = nan
    allocate(this%leafc_storage_xfer_acc_patch (begp:endp)) ; this%leafc_storage_xfer_acc_patch (:) = nan
    allocate(this%storage_cdemand_patch        (begp:endp)) ; this%storage_cdemand_patch        (:) = nan
    allocate(this%frootc_patch             (begp:endp)) ; this%frootc_patch             (:) = nan
    allocate(this%frootc_storage_patch     (begp:endp)) ; this%frootc_storage_patch     (:) = nan
    allocate(this%frootc_xfer_patch        (begp:endp)) ; this%frootc_xfer_patch        (:) = nan
    allocate(this%livestemc_patch          (begp:endp)) ; this%livestemc_patch          (:) = nan
    allocate(this%livestemc_storage_patch  (begp:endp)) ; this%livestemc_storage_patch  (:) = nan
    allocate(this%livestemc_xfer_patch     (begp:endp)) ; this%livestemc_xfer_patch     (:) = nan
    allocate(this%deadstemc_patch          (begp:endp)) ; this%deadstemc_patch          (:) = nan
    allocate(this%deadstemc_storage_patch  (begp:endp)) ; this%deadstemc_storage_patch  (:) = nan
    allocate(this%deadstemc_xfer_patch     (begp:endp)) ; this%deadstemc_xfer_patch     (:) = nan
    allocate(this%livecrootc_patch         (begp:endp)) ; this%livecrootc_patch         (:) = nan
    allocate(this%livecrootc_storage_patch (begp:endp)) ; this%livecrootc_storage_patch (:) = nan
    allocate(this%livecrootc_xfer_patch    (begp:endp)) ; this%livecrootc_xfer_patch    (:) = nan
    allocate(this%deadcrootc_patch         (begp:endp)) ; this%deadcrootc_patch         (:) = nan
    allocate(this%deadcrootc_storage_patch (begp:endp)) ; this%deadcrootc_storage_patch (:) = nan
    allocate(this%deadcrootc_xfer_patch    (begp:endp)) ; this%deadcrootc_xfer_patch    (:) = nan
    allocate(this%gresp_storage_patch      (begp:endp)) ; this%gresp_storage_patch      (:) = nan
    allocate(this%gresp_xfer_patch         (begp:endp)) ; this%gresp_xfer_patch         (:) = nan
    allocate(this%cpool_patch              (begp:endp)) ; this%cpool_patch              (:) = nan
    allocate(this%xsmrpool_patch           (begp:endp)) ; this%xsmrpool_patch           (:) = nan
    allocate(this%ctrunc_patch             (begp:endp)) ; this%ctrunc_patch             (:) = nan
    allocate(this%dispvegc_patch           (begp:endp)) ; this%dispvegc_patch           (:) = nan
    allocate(this%storvegc_patch           (begp:endp)) ; this%storvegc_patch           (:) = nan
    allocate(this%leafcmax_patch           (begp:endp)) ; this%leafcmax_patch           (:) = nan
    allocate(this%totc_patch               (begp:endp))  ; this%totc_patch               (:) = nan
    allocate(this%grainc_patch             (begp:endp)) ; this%grainc_patch             (:) = nan
    allocate(this%grainc_storage_patch     (begp:endp)) ; this%grainc_storage_patch     (:) = nan
    allocate(this%grainc_xfer_patch        (begp:endp)) ; this%grainc_xfer_patch        (:) = nan
    allocate(this%woodc_patch              (begp:endp)) ; this%woodc_patch              (:) = nan     

    allocate(this%cropseedc_deficit_patch  (begp:endp)) ; this%cropseedc_deficit_patch  (:) = nan
    allocate(this%seedc_grc                (begg:endg)) ; this%seedc_grc                (:) = nan
    allocate(this%rootc_col                (begc:endc)) ; this%rootc_col                (:) = nan
    allocate(this%leafc_col                (begc:endc)) ; this%leafc_col                (:) = nan
    allocate(this%deadstemc_col            (begc:endc)) ; this%deadstemc_col            (:) = nan
    allocate(this%fuelc_col                (begc:endc)) ; this%fuelc_col                (:) = nan
    allocate(this%fuelc_crop_col           (begc:endc)) ; this%fuelc_crop_col           (:) = nan

    allocate(this%totvegc_patch            (begp:endp)) ; this%totvegc_patch            (:) = nan
    allocate(this%totvegc_col              (begc:endc)) ; this%totvegc_col              (:) = nan

    allocate(this%totc_p2c_col             (begc:endc)) ; this%totc_p2c_col             (:) = nan
    allocate(this%totc_col                 (begc:endc)) ; this%totc_col                 (:) = nan
    allocate(this%totecosysc_col           (begc:endc)) ; this%totecosysc_col           (:) = nan

  end subroutine InitAllocate

  !------------------------------------------------------------------------
  subroutine InitHistory(this, bounds, carbon_type)
    !
    ! !DESCRIPTION:
    ! add history fields for all CN variables, always set as default='inactive'
    !
    ! !USES:
    use clm_varctl , only : use_c13, use_c14
    use histFileMod, only : hist_addfld1d, hist_addfld2d, hist_addfld_decomp 
    !
    ! !ARGUMENTS:
    class (cnveg_carbonstate_type) :: this
    type(bounds_type)         , intent(in) :: bounds 
    character(len=*)          , intent(in) :: carbon_type ! one of ['c12', c13','c14']
    !
    ! !LOCAL VARIABLES:
    integer           :: k,l,ii,jj 
    character(10)     :: active
    integer           :: begp,endp
    integer           :: begc,endc
    integer           :: begg,endg 
    character(24)     :: fieldname
    character(100)    :: longname
    real(r8), pointer :: data1dptr(:)   ! temp. pointer for slicing larger arrays
    real(r8), pointer :: data2dptr(:,:) ! temp. pointer for slicing larger arrays
    !---------------------------------------------------------------------

    begp = bounds%begp; endp = bounds%endp
    begc = bounds%begc; endc = bounds%endc
    begg = bounds%begg; endg = bounds%endg

    !-------------------------------
    ! C12 state variables
    !-------------------------------

    if (carbon_type == 'c12') then

       if (use_crop) then
          this%grainc_patch(begp:endp) = spval
          call hist_addfld1d (fname='GRAINC', units='gC/m^2', &
               avgflag='A', long_name='grain C (does not equal yield)', &
               ptr_patch=this%grainc_patch)
          this%cropseedc_deficit_patch(begp:endp) = spval
          call hist_addfld1d (fname='CROPSEEDC_DEFICIT', units='gC/m^2', &
               avgflag='A', long_name='C used for crop seed that needs to be repaid', &
               ptr_patch=this%cropseedc_deficit_patch)
       end if
       
       this%woodc_patch(begp:endp) = spval
       call hist_addfld1d (fname='WOODC', units='gC/m^2', &
            avgflag='A', long_name='wood C', &
            ptr_patch=this%woodc_patch)

       this%leafc_patch(begp:endp) = spval
       call hist_addfld1d (fname='LEAFC', units='gC/m^2', &
            avgflag='A', long_name='leaf C', &
            ptr_patch=this%leafc_patch)

       this%leafc_storage_patch(begp:endp) = spval
       call hist_addfld1d (fname='LEAFC_STORAGE', units='gC/m^2', &
            avgflag='A', long_name='leaf C storage', &
            ptr_patch=this%leafc_storage_patch, default='inactive')    

       this%leafc_xfer_patch(begp:endp) = spval
       call hist_addfld1d (fname='LEAFC_XFER', units='gC/m^2', &
            avgflag='A', long_name='leaf C transfer', &
            ptr_patch=this%leafc_xfer_patch, default='inactive')    

       this%leafc_storage_xfer_acc_patch(begp:endp) = spval
       call hist_addfld1d (fname='LEAFC_STORAGE_XFER_ACC', units='gC/m^2', &
            avgflag='A', long_name='Accumulated leaf C transfer', &
            ptr_patch=this%leafc_storage_xfer_acc_patch, default='inactive')

       this%storage_cdemand_patch(begp:endp) = spval
       call hist_addfld1d (fname='STORAGE_CDEMAND', units='gC/m^2', &
            avgflag='A', long_name='C use from the C storage pool', &
            ptr_patch=this%storage_cdemand_patch, default='inactive')

       this%frootc_patch(begp:endp) = spval
       call hist_addfld1d (fname='FROOTC', units='gC/m^2', &
            avgflag='A', long_name='fine root C', &
            ptr_patch=this%frootc_patch)

       this%frootc_storage_patch(begp:endp) = spval
       call hist_addfld1d (fname='FROOTC_STORAGE', units='gC/m^2', &
            avgflag='A', long_name='fine root C storage', &
            ptr_patch=this%frootc_storage_patch, default='inactive')   

       this%frootc_xfer_patch(begp:endp) = spval
       call hist_addfld1d (fname='FROOTC_XFER', units='gC/m^2', &
            avgflag='A', long_name='fine root C transfer', &
            ptr_patch=this%frootc_xfer_patch, default='inactive')    

       this%livestemc_patch(begp:endp) = spval
       call hist_addfld1d (fname='LIVESTEMC', units='gC/m^2', &
            avgflag='A', long_name='live stem C', &
            ptr_patch=this%livestemc_patch)

       this%livestemc_storage_patch(begp:endp) = spval
       call hist_addfld1d (fname='LIVESTEMC_STORAGE', units='gC/m^2', &
            avgflag='A', long_name='live stem C storage', &
            ptr_patch=this%livestemc_storage_patch, default='inactive')    

       this%livestemc_xfer_patch(begp:endp) = spval
       call hist_addfld1d (fname='LIVESTEMC_XFER', units='gC/m^2', &
            avgflag='A', long_name='live stem C transfer', &
            ptr_patch=this%livestemc_xfer_patch, default='inactive')     

       this%deadstemc_patch(begp:endp) = spval
       call hist_addfld1d (fname='DEADSTEMC', units='gC/m^2', &
            avgflag='A', long_name='dead stem C', &
            ptr_patch=this%deadstemc_patch)

       this%deadstemc_storage_patch(begp:endp) = spval
       call hist_addfld1d (fname='DEADSTEMC_STORAGE', units='gC/m^2', &
            avgflag='A', long_name='dead stem C storage', &
            ptr_patch=this%deadstemc_storage_patch, default='inactive')    

       this%deadstemc_xfer_patch(begp:endp) = spval
       call hist_addfld1d (fname='DEADSTEMC_XFER', units='gC/m^2', &
            avgflag='A', long_name='dead stem C transfer', &
            ptr_patch=this%deadstemc_xfer_patch, default='inactive')    

       this%livecrootc_patch(begp:endp) = spval
       call hist_addfld1d (fname='LIVECROOTC', units='gC/m^2', &
            avgflag='A', long_name='live coarse root C', &
            ptr_patch=this%livecrootc_patch)

       this%livecrootc_storage_patch(begp:endp) = spval
       call hist_addfld1d (fname='LIVECROOTC_STORAGE', units='gC/m^2', &
            avgflag='A', long_name='live coarse root C storage', &
            ptr_patch=this%livecrootc_storage_patch, default='inactive')     

       this%livecrootc_xfer_patch(begp:endp) = spval
       call hist_addfld1d (fname='LIVECROOTC_XFER', units='gC/m^2', &
            avgflag='A', long_name='live coarse root C transfer', &
            ptr_patch=this%livecrootc_xfer_patch, default='inactive')    

       this%deadcrootc_patch(begp:endp) = spval
       call hist_addfld1d (fname='DEADCROOTC', units='gC/m^2', &
            avgflag='A', long_name='dead coarse root C', &
            ptr_patch=this%deadcrootc_patch)

       this%deadcrootc_storage_patch(begp:endp) = spval
       call hist_addfld1d (fname='DEADCROOTC_STORAGE', units='gC/m^2', &
            avgflag='A', long_name='dead coarse root C storage', &
            ptr_patch=this%deadcrootc_storage_patch, default='inactive')   

       this%deadcrootc_xfer_patch(begp:endp) = spval
       call hist_addfld1d (fname='DEADCROOTC_XFER', units='gC/m^2', &
            avgflag='A', long_name='dead coarse root C transfer', &
            ptr_patch=this%deadcrootc_xfer_patch, default='inactive')   

       this%gresp_storage_patch(begp:endp) = spval
       call hist_addfld1d (fname='GRESP_STORAGE', units='gC/m^2', &
            avgflag='A', long_name='growth respiration storage', &
            ptr_patch=this%gresp_storage_patch, default='inactive')    

       this%gresp_xfer_patch(begp:endp) = spval
       call hist_addfld1d (fname='GRESP_XFER', units='gC/m^2', &
            avgflag='A', long_name='growth respiration transfer', &
            ptr_patch=this%gresp_xfer_patch, default='inactive')     

       this%cpool_patch(begp:endp) = spval
       call hist_addfld1d (fname='CPOOL', units='gC/m^2', &
            avgflag='A', long_name='temporary photosynthate C pool', &
            ptr_patch=this%cpool_patch)

       this%xsmrpool_patch(begp:endp) = spval
       call hist_addfld1d (fname='XSMRPOOL', units='gC/m^2', &
            avgflag='A', long_name='temporary photosynthate C pool', &
            ptr_patch=this%xsmrpool_patch)

       this%ctrunc_patch(begp:endp) = spval
       call hist_addfld1d (fname='PFT_CTRUNC', units='gC/m^2', &
            avgflag='A', long_name='patch-level sink for C truncation', &
            ptr_patch=this%ctrunc_patch, default='inactive')

       this%dispvegc_patch(begp:endp) = spval
       call hist_addfld1d (fname='DISPVEGC', units='gC/m^2', &
            avgflag='A', long_name='displayed veg carbon, excluding storage and cpool', &
            ptr_patch=this%dispvegc_patch)

       this%storvegc_patch(begp:endp) = spval
       call hist_addfld1d (fname='STORVEGC', units='gC/m^2', &
            avgflag='A', long_name='stored vegetation carbon, excluding cpool', &
            ptr_patch=this%storvegc_patch)

       this%totvegc_patch(begp:endp) = spval
       call hist_addfld1d (fname='TOTVEGC', units='gC/m^2', &
            avgflag='A', long_name='total vegetation carbon, excluding cpool', &
            ptr_patch=this%totvegc_patch)

       this%totc_patch(begp:endp) = spval
       call hist_addfld1d (fname='TOTPFTC', units='gC/m^2', &
            avgflag='A', long_name='total patch-level carbon, including cpool', &
            ptr_patch=this%totc_patch)

       this%seedc_grc(begg:endg) = spval
       call hist_addfld1d (fname='SEEDC', units='gC/m^2', &
            avgflag='A', long_name='pool for seeding new PFTs via dynamic landcover', &
            ptr_gcell=this%seedc_grc)

       this%fuelc_col(begc:endc) = spval
       call hist_addfld1d (fname='FUELC', units='gC/m^2', &
            avgflag='A', long_name='fuel load', &
            ptr_col=this%fuelc_col)

       this%totc_col(begc:endc) = spval
       call hist_addfld1d (fname='TOTCOLC', units='gC/m^2', &
            avgflag='A', long_name='total column carbon, incl veg and cpool but excl product pools', &
            ptr_col=this%totc_col)

       this%totecosysc_col(begc:endc) = spval
       call hist_addfld1d (fname='TOTECOSYSC', units='gC/m^2', &
            avgflag='A', long_name='total ecosystem carbon, incl veg but excl cpool and product pools', &
            ptr_col=this%totecosysc_col)

    end if

    !-------------------------------
    ! C13 state variables 
    !-------------------------------

    if ( carbon_type == 'c13' ) then

       this%leafc_patch(begp:endp) = spval
       call hist_addfld1d (fname='C13_LEAFC', units='gC13/m^2', &
            avgflag='A', long_name='C13 leaf C', &
            ptr_patch=this%leafc_patch, default='inactive')

       this%leafc_storage_patch(begp:endp) = spval
       call hist_addfld1d (fname='C13_LEAFC_STORAGE', units='gC13/m^2', &
            avgflag='A', long_name='C13 leaf C storage', &
            ptr_patch=this%leafc_storage_patch, default='inactive')

       this%leafc_xfer_patch(begp:endp) = spval
       call hist_addfld1d (fname='C13_LEAFC_XFER', units='gC13/m^2', &
            avgflag='A', long_name='C13 leaf C transfer', &
            ptr_patch=this%leafc_xfer_patch, default='inactive')

       this%leafc_storage_xfer_acc_patch(begp:endp) = spval
       call hist_addfld1d (fname='C13_LEAFC_STORAGE_XFER_ACC', units='gC13/m^2', &
            avgflag='A', long_name='Accumulated C13 leaf C transfer', &
            ptr_patch=this%leafc_storage_xfer_acc_patch, default='inactive')

       this%frootc_patch(begp:endp) = spval
       call hist_addfld1d (fname='C13_FROOTC', units='gC13/m^2', &
            avgflag='A', long_name='C13 fine root C', &
            ptr_patch=this%frootc_patch, default='inactive')

       this%frootc_storage_patch(begp:endp) = spval
       call hist_addfld1d (fname='C13_FROOTC_STORAGE', units='gC13/m^2', &
            avgflag='A', long_name='C13 fine root C storage', &
            ptr_patch=this%frootc_storage_patch, default='inactive')

       this%frootc_xfer_patch(begp:endp) = spval
       call hist_addfld1d (fname='C13_FROOTC_XFER', units='gC13/m^2', &
            avgflag='A', long_name='C13 fine root C transfer', &
            ptr_patch=this%frootc_xfer_patch, default='inactive')

       this%livestemc_patch(begp:endp) = spval
       call hist_addfld1d (fname='C13_LIVESTEMC', units='gC13/m^2', &
            avgflag='A', long_name='C13 live stem C', &
            ptr_patch=this%livestemc_patch, default='inactive')

       this%livestemc_storage_patch(begp:endp) = spval
       call hist_addfld1d (fname='C13_LIVESTEMC_STORAGE', units='gC13/m^2', &
            avgflag='A', long_name='C13 live stem C storage', &
            ptr_patch=this%livestemc_storage_patch, default='inactive')

       this%livestemc_xfer_patch(begp:endp) = spval
       call hist_addfld1d (fname='C13_LIVESTEMC_XFER', units='gC13/m^2', &
            avgflag='A', long_name='C13 live stem C transfer', &
            ptr_patch=this%livestemc_xfer_patch, default='inactive')

       this%deadstemc_patch(begp:endp) = spval
       call hist_addfld1d (fname='C13_DEADSTEMC', units='gC13/m^2', &
            avgflag='A', long_name='C13 dead stem C', &
            ptr_patch=this%deadstemc_patch, default='inactive')

       this%deadstemc_storage_patch(begp:endp) = spval
       call hist_addfld1d (fname='C13_DEADSTEMC_STORAGE', units='gC13/m^2', &
            avgflag='A', long_name='C13 dead stem C storage', &
            ptr_patch=this%deadstemc_storage_patch, default='inactive')

       this%deadstemc_xfer_patch(begp:endp) = spval
       call hist_addfld1d (fname='C13_DEADSTEMC_XFER', units='gC13/m^2', &
            avgflag='A', long_name='C13 dead stem C transfer', &
            ptr_patch=this%deadstemc_xfer_patch, default='inactive')

       this%livecrootc_patch(begp:endp) = spval
       call hist_addfld1d (fname='C13_LIVECROOTC', units='gC13/m^2', &
            avgflag='A', long_name='C13 live coarse root C', &
            ptr_patch=this%livecrootc_patch, default='inactive')

       this%livecrootc_storage_patch(begp:endp) = spval
       call hist_addfld1d (fname='C13_LIVECROOTC_STORAGE', units='gC13/m^2', &
            avgflag='A', long_name='C13 live coarse root C storage', &
            ptr_patch=this%livecrootc_storage_patch, default='inactive')

       this%livecrootc_xfer_patch(begp:endp) = spval
       call hist_addfld1d (fname='C13_LIVECROOTC_XFER', units='gC13/m^2', &
            avgflag='A', long_name='C13 live coarse root C transfer', &
            ptr_patch=this%livecrootc_xfer_patch, default='inactive')

       this%deadcrootc_patch(begp:endp) = spval
       call hist_addfld1d (fname='C13_DEADCROOTC', units='gC13/m^2', &
            avgflag='A', long_name='C13 dead coarse root C', &
            ptr_patch=this%deadcrootc_patch, default='inactive')

       this%deadcrootc_storage_patch(begp:endp) = spval
       call hist_addfld1d (fname='C13_DEADCROOTC_STORAGE', units='gC13/m^2', &
            avgflag='A', long_name='C13 dead coarse root C storage', &
            ptr_patch=this%deadcrootc_storage_patch,  default='inactive')

       this%deadcrootc_xfer_patch(begp:endp) = spval
       call hist_addfld1d (fname='C13_DEADCROOTC_XFER', units='gC13/m^2', &
            avgflag='A', long_name='C13 dead coarse root C transfer', &
            ptr_patch=this%deadcrootc_xfer_patch, default='inactive')

       this%gresp_storage_patch(begp:endp) = spval
       call hist_addfld1d (fname='C13_GRESP_STORAGE', units='gC13/m^2', &
            avgflag='A', long_name='C13 growth respiration storage', &
            ptr_patch=this%gresp_storage_patch, default='inactive')

       this%gresp_xfer_patch(begp:endp) = spval
       call hist_addfld1d (fname='C13_GRESP_XFER', units='gC13/m^2', &
            avgflag='A', long_name='C13 growth respiration transfer', &
            ptr_patch=this%gresp_xfer_patch, default='inactive')

       this%cpool_patch(begp:endp) = spval
       call hist_addfld1d (fname='C13_CPOOL', units='gC13/m^2', &
            avgflag='A', long_name='C13 temporary photosynthate C pool', &
            ptr_patch=this%cpool_patch, default='inactive')

       this%xsmrpool_patch(begp:endp) = spval
       call hist_addfld1d (fname='C13_XSMRPOOL', units='gC13/m^2', &
            avgflag='A', long_name='C13 temporary photosynthate C pool', &
            ptr_patch=this%xsmrpool_patch, default='inactive')

       this%ctrunc_patch(begp:endp) = spval
       call hist_addfld1d (fname='C13_PFT_CTRUNC', units='gC13/m^2', &
            avgflag='A', long_name='C13 patch-level sink for C truncation', &
            ptr_patch=this%ctrunc_patch, default='inactive')

       this%dispvegc_patch(begp:endp) = spval
       call hist_addfld1d (fname='C13_DISPVEGC', units='gC13/m^2', &
            avgflag='A', long_name='C13 displayed veg carbon, excluding storage and cpool', &
            ptr_patch=this%dispvegc_patch, default='inactive')

       this%storvegc_patch(begp:endp) = spval
       call hist_addfld1d (fname='C13_STORVEGC', units='gC13/m^2', &
            avgflag='A', long_name='C13 stored vegetation carbon, excluding cpool', &
            ptr_patch=this%storvegc_patch, default='inactive')

       this%totvegc_patch(begp:endp) = spval
       call hist_addfld1d (fname='C13_TOTVEGC', units='gC13/m^2', &
            avgflag='A', long_name='C13 total vegetation carbon, excluding cpool', &
            ptr_patch=this%totvegc_patch)

       this%totc_patch(begp:endp) = spval
       call hist_addfld1d (fname='C13_TOTPFTC', units='gC13/m^2', &
            avgflag='A', long_name='C13 total patch-level carbon, including cpool', &
            ptr_patch=this%totc_patch, default='inactive')

       this%seedc_grc(begg:endg) = spval
       call hist_addfld1d (fname='C13_SEEDC', units='gC13/m^2', &
            avgflag='A', long_name='C13 pool for seeding new PFTs via dynamic landcover', &
            ptr_gcell=this%seedc_grc, default='inactive')

       this%totc_col(begc:endc) = spval
       call hist_addfld1d (fname='C13_TOTCOLC', units='gC13/m^2', &
            avgflag='A', long_name='C13 total column carbon, incl veg and cpool but excl product pools', &
            ptr_col=this%totc_col, default='inactive')

       this%totecosysc_col(begc:endc) = spval
       call hist_addfld1d (fname='C13_TOTECOSYSC', units='gC13/m^2', &
            avgflag='A', long_name='C13 total ecosystem carbon, incl veg but excl cpool and product pools', &
            ptr_col=this%totecosysc_col)

       if (use_crop) then
          this%grainc_patch(begp:endp) = spval
          call hist_addfld1d (fname='C13_GRAINC', units='gC/m^2', &
               avgflag='A', long_name='C13 grain C (does not equal yield)', &
               ptr_patch=this%grainc_patch, default='inactive')
          this%cropseedc_deficit_patch(begp:endp) = spval
          call hist_addfld1d (fname='C13_CROPSEEDC_DEFICIT', units='gC/m^2', &
               avgflag='A', long_name='C13 C used for crop seed that needs to be repaid', &
               ptr_patch=this%cropseedc_deficit_patch, default='inactive')
       end if


    endif

    !-------------------------------
    ! C14 state variables 
    !-------------------------------

    if ( carbon_type == 'c14') then

       this%leafc_patch(begp:endp) = spval
       call hist_addfld1d (fname='C14_LEAFC', units='gC14/m^2', &
            avgflag='A', long_name='C14 leaf C', &
            ptr_patch=this%leafc_patch, default='inactive')

       this%leafc_storage_patch(begp:endp) = spval
       call hist_addfld1d (fname='C14_LEAFC_STORAGE', units='gC14/m^2', &
            avgflag='A', long_name='C14 leaf C storage', &
            ptr_patch=this%leafc_storage_patch, default='inactive')

       this%leafc_xfer_patch(begp:endp) = spval
       call hist_addfld1d (fname='C14_LEAFC_XFER', units='gC14/m^2', &
            avgflag='A', long_name='C14 leaf C transfer', &
            ptr_patch=this%leafc_xfer_patch, default='inactive')

        this%leafc_storage_xfer_acc_patch(begp:endp) = spval
        call hist_addfld1d (fname='C14_LEAFC_STORAGE_XFER_ACC', units='gC14/m^2', &
             avgflag='A', long_name='Accumulated C14 leaf C transfer', &
             ptr_patch=this%leafc_storage_xfer_acc_patch, default='inactive')

       this%frootc_patch(begp:endp) = spval
       call hist_addfld1d (fname='C14_FROOTC', units='gC14/m^2', &
            avgflag='A', long_name='C14 fine root C', &
            ptr_patch=this%frootc_patch, default='inactive')

       this%frootc_storage_patch(begp:endp) = spval
       call hist_addfld1d (fname='C14_FROOTC_STORAGE', units='gC14/m^2', &
            avgflag='A', long_name='C14 fine root C storage', &
            ptr_patch=this%frootc_storage_patch, default='inactive')

       this%frootc_xfer_patch(begp:endp) = spval
       call hist_addfld1d (fname='C14_FROOTC_XFER', units='gC14/m^2', &
            avgflag='A', long_name='C14 fine root C transfer', &
            ptr_patch=this%frootc_xfer_patch, default='inactive')

       this%livestemc_patch(begp:endp) = spval
       call hist_addfld1d (fname='C14_LIVESTEMC', units='gC14/m^2', &
            avgflag='A', long_name='C14 live stem C', &
            ptr_patch=this%livestemc_patch, default='inactive')

       this%livestemc_storage_patch(begp:endp) = spval
       call hist_addfld1d (fname='C14_LIVESTEMC_STORAGE', units='gC14/m^2', &
            avgflag='A', long_name='C14 live stem C storage', &
            ptr_patch=this%livestemc_storage_patch, default='inactive')

       this%livestemc_xfer_patch(begp:endp) = spval
       call hist_addfld1d (fname='C14_LIVESTEMC_XFER', units='gC14/m^2', &
            avgflag='A', long_name='C14 live stem C transfer', &
            ptr_patch=this%livestemc_xfer_patch, default='inactive')

       this%deadstemc_patch(begp:endp) = spval
       call hist_addfld1d (fname='C14_DEADSTEMC', units='gC14/m^2', &
            avgflag='A', long_name='C14 dead stem C', &
            ptr_patch=this%deadstemc_patch, default='inactive')

       this%deadstemc_storage_patch(begp:endp) = spval
       call hist_addfld1d (fname='C14_DEADSTEMC_STORAGE', units='gC14/m^2', &
            avgflag='A', long_name='C14 dead stem C storage', &
            ptr_patch=this%deadstemc_storage_patch, default='inactive')

       this%deadstemc_xfer_patch(begp:endp) = spval
       call hist_addfld1d (fname='C14_DEADSTEMC_XFER', units='gC14/m^2', &
            avgflag='A', long_name='C14 dead stem C transfer', &
            ptr_patch=this%deadstemc_xfer_patch, default='inactive')

       this%livecrootc_patch(begp:endp) = spval
       call hist_addfld1d (fname='C14_LIVECROOTC', units='gC14/m^2', &
            avgflag='A', long_name='C14 live coarse root C', &
            ptr_patch=this%livecrootc_patch, default='inactive')

       this%livecrootc_storage_patch(begp:endp) = spval
       call hist_addfld1d (fname='C14_LIVECROOTC_STORAGE', units='gC14/m^2', &
            avgflag='A', long_name='C14 live coarse root C storage', &
            ptr_patch=this%livecrootc_storage_patch, default='inactive')

       this%livecrootc_xfer_patch(begp:endp) = spval
       call hist_addfld1d (fname='C14_LIVECROOTC_XFER', units='gC14/m^2', &
            avgflag='A', long_name='C14 live coarse root C transfer', &
            ptr_patch=this%livecrootc_xfer_patch, default='inactive')

       this%deadcrootc_patch(begp:endp) = spval
       call hist_addfld1d (fname='C14_DEADCROOTC', units='gC14/m^2', &
            avgflag='A', long_name='C14 dead coarse root C', &
            ptr_patch=this%deadcrootc_patch, default='inactive')

       this%deadcrootc_storage_patch(begp:endp) = spval
       call hist_addfld1d (fname='C14_DEADCROOTC_STORAGE', units='gC14/m^2', &
            avgflag='A', long_name='C14 dead coarse root C storage', &
            ptr_patch=this%deadcrootc_storage_patch,  default='inactive')

       this%deadcrootc_xfer_patch(begp:endp) = spval
       call hist_addfld1d (fname='C14_DEADCROOTC_XFER', units='gC14/m^2', &
            avgflag='A', long_name='C14 dead coarse root C transfer', &
            ptr_patch=this%deadcrootc_xfer_patch, default='inactive')

       this%gresp_storage_patch(begp:endp) = spval
       call hist_addfld1d (fname='C14_GRESP_STORAGE', units='gC14/m^2', &
            avgflag='A', long_name='C14 growth respiration storage', &
            ptr_patch=this%gresp_storage_patch, default='inactive')

       this%gresp_xfer_patch(begp:endp) = spval
       call hist_addfld1d (fname='C14_GRESP_XFER', units='gC14/m^2', &
            avgflag='A', long_name='C14 growth respiration transfer', &
            ptr_patch=this%gresp_xfer_patch, default='inactive')

       this%cpool_patch(begp:endp) = spval
       call hist_addfld1d (fname='C14_CPOOL', units='gC14/m^2', &
            avgflag='A', long_name='C14 temporary photosynthate C pool', &
            ptr_patch=this%cpool_patch, default='inactive')

       this%xsmrpool_patch(begp:endp) = spval
       call hist_addfld1d (fname='C14_XSMRPOOL', units='gC14/m^2', &
            avgflag='A', long_name='C14 temporary photosynthate C pool', &
            ptr_patch=this%xsmrpool_patch, default='inactive')

       this%ctrunc_patch(begp:endp) = spval
       call hist_addfld1d (fname='C14_PFT_CTRUNC', units='gC14/m^2', &
            avgflag='A', long_name='C14 patch-level sink for C truncation', &
            ptr_patch=this%ctrunc_patch, default='inactive')

       this%dispvegc_patch(begp:endp) = spval
       call hist_addfld1d (fname='C14_DISPVEGC', units='gC14/m^2', &
            avgflag='A', long_name='C14 displayed veg carbon, excluding storage and cpool', &
            ptr_patch=this%dispvegc_patch, default='inactive')

       this%storvegc_patch(begp:endp) = spval
       call hist_addfld1d (fname='C14_STORVEGC', units='gC14/m^2', &
            avgflag='A', long_name='C14 stored vegetation carbon, excluding cpool', &
            ptr_patch=this%storvegc_patch, default='inactive')

       this%totvegc_patch(begp:endp) = spval
       call hist_addfld1d (fname='C14_TOTVEGC', units='gC14/m^2', &
            avgflag='A', long_name='C14 total vegetation carbon, excluding cpool', &
            ptr_patch=this%totvegc_patch)

       this%totc_patch(begp:endp) = spval
       call hist_addfld1d (fname='C14_TOTPFTC', units='gC14/m^2', &
            avgflag='A', long_name='C14 total patch-level carbon, including cpool', &
            ptr_patch=this%totc_patch, default='inactive')

       this%seedc_grc(begg:endg) = spval
       call hist_addfld1d (fname='C14_SEEDC', units='gC14/m^2', &
            avgflag='A', long_name='C14 pool for seeding new PFTs via dynamic landcover', &
            ptr_gcell=this%seedc_grc, default='inactive')

       this%totc_col(begc:endc) = spval
       call hist_addfld1d (fname='C14_TOTCOLC', units='gC14/m^2', &
            avgflag='A', long_name='C14 total column carbon, incl veg and cpool but excl product pools', &
            ptr_col=this%totc_col, default='inactive')

       this%totecosysc_col(begc:endc) = spval
       call hist_addfld1d (fname='C14_TOTECOSYSC', units='gC14/m^2', &
            avgflag='A', long_name='C14 total ecosystem carbon, incl veg but excl cpool and product pools', &
            ptr_col=this%totecosysc_col)

       if (use_crop) then
          this%grainc_patch(begp:endp) = spval
          call hist_addfld1d (fname='C14_GRAINC', units='gC/m^2', &
               avgflag='A', long_name='C14 grain C (does not equal yield)', &
               ptr_patch=this%grainc_patch, default='inactive')
          this%cropseedc_deficit_patch(begp:endp) = spval
          call hist_addfld1d (fname='C14_CROPSEEDC_DEFICIT', units='gC/m^2', &
               avgflag='A', long_name='C14 C used for crop seed that needs to be repaid', &
               ptr_patch=this%cropseedc_deficit_patch, default='inactive')
       end if


    endif

  end subroutine InitHistory

  !-----------------------------------------------------------------------
  subroutine InitCold(this, bounds, ratio, carbon_type, c12_cnveg_carbonstate_inst)
    !
    ! !DESCRIPTION:
    ! Initializes time varying variables used only in coupled carbon-nitrogen mode (CN):
    !
    ! !USES:
    use landunit_varcon	 , only : istsoil, istcrop 
    use clm_time_manager , only : is_restart, get_nstep
    use clm_varctl, only : MM_Nuptake_opt    
    !
    ! !ARGUMENTS:
    class(cnveg_carbonstate_type)                       :: this 
    type(bounds_type)            , intent(in)           :: bounds  
    real(r8)                     , intent(in)           :: ratio              ! Standard isotope ratio
    character(len=*)             , intent(in)           :: carbon_type        ! 'c12' or 'c13' or 'c14'
    type(cnveg_carbonstate_type) , optional, intent(in) :: c12_cnveg_carbonstate_inst
    !
    ! !LOCAL VARIABLES:
    integer  :: p,c,l,g,j,k,i
    integer  :: fc                                       ! filter index
    integer  :: num_special_col                          ! number of good values in special_col filter
    integer  :: num_special_patch                        ! number of good values in special_patch filter
    integer  :: special_col(bounds%endc-bounds%begc+1)   ! special landunit filter - columns
    integer  :: special_patch(bounds%endp-bounds%begp+1) ! special landunit filter - patches
    !-----------------------------------------------------------------------

    if (carbon_type == 'c13' .or. carbon_type == 'c14') then
       if (.not. present(c12_cnveg_carbonstate_inst)) then
          call endrun(msg=' ERROR: for C13 or C14 must pass in c12_cnveg_carbonstate_inst as argument' //&
               errMsg(sourcefile, __LINE__))
       end if
    end if

    ! Set column filters

    num_special_col = 0
    do c = bounds%begc, bounds%endc
       l = col%landunit(c)
       if (lun%ifspecial(l)) then
          num_special_col = num_special_col + 1
          special_col(num_special_col) = c
       end if
    end do

    ! Set patch filters

    num_special_patch = 0
    do p = bounds%begp,bounds%endp
       l = patch%landunit(p)
       if (lun%ifspecial(l)) then
          num_special_patch = num_special_patch + 1
          special_patch(num_special_patch) = p
       end if
    end do

    !-----------------------------------------------
    ! initialize patch-level carbon state variables
    !-----------------------------------------------

    do p = bounds%begp,bounds%endp

       this%leafcmax_patch(p) = 0._r8

       l = patch%landunit(p)
       if (lun%itype(l) == istsoil .or. lun%itype(l) == istcrop) then

          if (patch%itype(p) == noveg) then
             this%leafc_patch(p)          = 0._r8
             this%leafc_storage_patch(p)  = 0._r8
             this%frootc_patch(p)         = 0._r8            
             this%frootc_storage_patch(p) = 0._r8    
          else
             if (pftcon%evergreen(patch%itype(p)) == 1._r8) then
                this%leafc_patch(p)          = cnvegcstate_const%initial_vegC * ratio     
                this%leafc_storage_patch(p)  = 0._r8
                this%frootc_patch(p)         = cnvegcstate_const%initial_vegC * ratio           
                this%frootc_storage_patch(p) = 0._r8    
             else if (patch%itype(p) >= npcropmin) then ! prognostic crop types
                this%leafc_patch(p)          = 0._r8
                this%leafc_storage_patch(p)  = 0._r8
                this%frootc_patch(p)         = 0._r8            
                this%frootc_storage_patch(p) = 0._r8    
             else
                this%leafc_patch(p)          = 0._r8
                this%leafc_storage_patch(p)  = cnvegcstate_const%initial_vegC * ratio   
                this%frootc_patch(p)         = 0._r8            
                this%frootc_storage_patch(p) = cnvegcstate_const%initial_vegC * ratio   
             end if
          end if
          this%leafc_xfer_patch(p) = 0._r8
          this%leafc_storage_xfer_acc_patch(p)  = 0._r8
          this%storage_cdemand_patch(p)         = 0._r8

          if (MM_Nuptake_opt .eqv. .false.) then  ! if not running in floating CN ratio option 
             this%frootc_patch(p)            = 0._r8 
             this%frootc_storage_patch(p)    = 0._r8 
          end if     
          this%frootc_xfer_patch(p)       = 0._r8 

          this%livestemc_patch(p)         = 0._r8 
          this%livestemc_storage_patch(p) = 0._r8 
          this%livestemc_xfer_patch(p)    = 0._r8 

          if (pftcon%woody(patch%itype(p)) == 1._r8) then
             this%deadstemc_patch(p) = 0.1_r8 * ratio
          else
             this%deadstemc_patch(p) = 0._r8 
          end if
          this%deadstemc_storage_patch(p)  = 0._r8 
          this%deadstemc_xfer_patch(p)     = 0._r8 

          this%livecrootc_patch(p)         = 0._r8 
          this%livecrootc_storage_patch(p) = 0._r8 
          this%livecrootc_xfer_patch(p)    = 0._r8 

          this%deadcrootc_patch(p)         = 0._r8 
          this%deadcrootc_storage_patch(p) = 0._r8 
          this%deadcrootc_xfer_patch(p)    = 0._r8 

          this%gresp_storage_patch(p)      = 0._r8 
          this%gresp_xfer_patch(p)         = 0._r8 

          this%cpool_patch(p)              = 0._r8 
          this%xsmrpool_patch(p)           = 0._r8 
          this%ctrunc_patch(p)             = 0._r8 
          this%dispvegc_patch(p)           = 0._r8 
          this%storvegc_patch(p)           = 0._r8 
          this%woodc_patch(p)              = 0._r8
          this%totc_patch(p)               = 0._r8 

          if ( use_crop )then
             this%grainc_patch(p)         = 0._r8 
             this%grainc_storage_patch(p) = 0._r8 
             this%grainc_xfer_patch(p)    = 0._r8 
             this%cropseedc_deficit_patch(p)  = 0._r8
          end if

       endif

    end do

    ! -----------------------------------------------
    ! initialize column-level variables
    ! -----------------------------------------------

    do c = bounds%begc, bounds%endc
       l = col%landunit(c)
       if (lun%itype(l) == istsoil .or. lun%itype(l) == istcrop) then
!          this%totgrainc_col(c)  = 0._r8

          ! total carbon pools
          this%totecosysc_col(c) = 0._r8
          this%totc_p2c_col(c)   = 0._r8
          this%totc_col(c)       = 0._r8
       end if
    end do


    do g = bounds%begg, bounds%endg
       this%seedc_grc(g) = 0._r8
    end do

    if ( .not. is_restart() .and. get_nstep() == 1 ) then

       do p = bounds%begp,bounds%endp
          if (pftcon%c3psn(patch%itype(p)) == 1._r8) then
             this%grainc_patch(p)            = c12_cnveg_carbonstate_inst%grainc_patch(p)         * c3_r2
             this%grainc_storage_patch(p)    = c12_cnveg_carbonstate_inst%grainc_storage_patch(p) * c3_r2
             this%grainc_xfer_patch(p)       = c12_cnveg_carbonstate_inst%grainc_xfer_patch(p)    * c3_r2
             this%dispvegc_patch(p)          = c12_cnveg_carbonstate_inst%dispvegc_patch(p)       * c3_r2
             this%storvegc_patch(p)          = c12_cnveg_carbonstate_inst%storvegc_patch(p)       * c3_r2
             this%totvegc_patch(p)           = c12_cnveg_carbonstate_inst%totvegc_patch(p)        * c3_r2
             this%totc_patch(p)              = c12_cnveg_carbonstate_inst%totc_patch(p)           * c3_r2
             this%woodc_patch(p)             = c12_cnveg_carbonstate_inst%woodc_patch(p)          * c3_r2
          else
             this%grainc_patch(p)            = c12_cnveg_carbonstate_inst%grainc_patch(p)         * c4_r2
             this%grainc_storage_patch(p)    = c12_cnveg_carbonstate_inst%grainc_storage_patch(p) * c4_r2
             this%grainc_xfer_patch(p)       = c12_cnveg_carbonstate_inst%grainc_xfer_patch(p)    * c4_r2
             this%dispvegc_patch(p)          = c12_cnveg_carbonstate_inst%dispvegc_patch(p)       * c4_r2
             this%storvegc_patch(p)          = c12_cnveg_carbonstate_inst%storvegc_patch(p)       * c4_r2
             this%totvegc_patch(p)           = c12_cnveg_carbonstate_inst%totvegc_patch(p)        * c4_r2
             this%totc_patch(p)              = c12_cnveg_carbonstate_inst%totc_patch(p)           * c4_r2
             this%woodc_patch(p)             = c12_cnveg_carbonstate_inst%woodc_patch(p)          * c4_r2
          end if
       end do
    end if

    ! initialize fields for special filters

    call this%SetValues (&
         num_patch=num_special_patch, filter_patch=special_patch, value_patch=0._r8, &
         num_column=num_special_col, filter_column=special_col, value_column=0._r8)

  end subroutine InitCold

  !-----------------------------------------------------------------------
  subroutine Restart ( this,  bounds, ncid, flag, carbon_type, reseed_dead_plants, &
                       c12_cnveg_carbonstate_inst, filter_reseed_patch, &
                       num_reseed_patch)
    !
    ! !DESCRIPTION: 
    ! Read/write CN restart data for carbon state
    !
    ! !USES:
    use shr_infnan_mod   , only : isnan => shr_infnan_isnan, nan => shr_infnan_nan, assignment(=)
    use clm_varcon       , only : c13ratio, c14ratio
    use clm_varctl       , only : spinup_state, use_cndv, MM_Nuptake_opt
    use clm_time_manager , only : get_nstep, is_restart, get_nstep
    use landunit_varcon	 , only : istsoil, istcrop 
    use spmdMod          , only : mpicom
    use shr_mpi_mod      , only : shr_mpi_sum
    use restUtilMod
    use ncdio_pio
    !
    ! !ARGUMENTS:
    class (cnveg_carbonstate_type)                               :: this
    type(bounds_type)                     , intent(in)           :: bounds 
    type(file_desc_t)                     , intent(inout)        :: ncid   ! netcdf id
    character(len=*)                      , intent(in)           :: flag   !'read' or 'write'
    character(len=*)                      , intent(in)           :: carbon_type ! 'c12' or 'c13' or 'c14'
    logical                               , intent(in)           :: reseed_dead_plants
    type (cnveg_carbonstate_type)         , intent(in), optional :: c12_cnveg_carbonstate_inst 
    integer                               , intent(out), optional :: filter_reseed_patch(:)
    integer                               , intent(out), optional :: num_reseed_patch
    !
    ! !LOCAL VARIABLES:
    integer            :: i,j,k,l,c,p
    real(r8)           :: ratio
    character(len=128) :: varname   ! temporary
    logical            :: readvar
    integer            :: idata
    logical            :: exit_spinup  = .false.
    logical            :: enter_spinup = .false.
    ! flags for comparing the model and restart decomposition cascades
    integer            :: decomp_cascade_state, restart_file_decomp_cascade_state 
    ! spinup state as read from restart file, for determining whether to enter or exit spinup mode.
    integer            :: restart_file_spinup_state
    integer            :: total_num_reseed_patch      ! Total number of patches to reseed across all processors

    !------------------------------------------------------------------------

    if (carbon_type == 'c13' .or. carbon_type == 'c14') then
       if (.not. present(c12_cnveg_carbonstate_inst)) then
          call endrun(msg=' ERROR: for C14 must pass in c12_cnveg_carbonstate_inst as argument' //&
               errMsg(sourcefile, __LINE__))
       end if
    end if
    if (carbon_type == 'c12') then
       ratio = 1._r8
    else if (carbon_type == 'c13') then
       ratio = c13ratio
    else if (carbon_type == 'c14') then
       ratio = c14ratio
    end if

    if ( (      present(num_reseed_patch) .and. .not. present(filter_reseed_patch)) &
    .or. (.not. present(num_reseed_patch) .and.       present(filter_reseed_patch) ) )then
       call endrun(msg=' ERROR: filter_reseed_patch and num_reseed_patch both need to be entered ' //&
       errMsg(sourcefile, __LINE__))
    end if
    if ( present(num_reseed_patch) )then
       num_reseed_patch = 0
       filter_reseed_patch(:) = -1
    end if

    !--------------------------------
    ! patch carbon state variables (c12)
    !--------------------------------

    if (carbon_type == 'c12') then
       call restartvar(ncid=ncid, flag=flag, varname='leafc', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%leafc_patch) 

       call restartvar(ncid=ncid, flag=flag, varname='leafc_storage', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%leafc_storage_patch) 

       call restartvar(ncid=ncid, flag=flag, varname='leafc_xfer', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%leafc_xfer_patch) 

       call restartvar(ncid=ncid, flag=flag, varname='leafc_storage_xfer_acc', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%leafc_storage_xfer_acc_patch)
 
       call restartvar(ncid=ncid, flag=flag, varname='storage_cdemand', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%storage_cdemand_patch)

       call restartvar(ncid=ncid, flag=flag, varname='frootc', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%frootc_patch) 

       call restartvar(ncid=ncid, flag=flag, varname='frootc_storage', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%frootc_storage_patch) 

       call restartvar(ncid=ncid, flag=flag, varname='frootc_xfer', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%frootc_xfer_patch) 

       call restartvar(ncid=ncid, flag=flag, varname='livestemc', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%livestemc_patch) 

       call restartvar(ncid=ncid, flag=flag, varname='livestemc_storage', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%livestemc_storage_patch) 

       call restartvar(ncid=ncid, flag=flag, varname='livestemc_xfer', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%livestemc_xfer_patch) 

       call restartvar(ncid=ncid, flag=flag, varname='deadstemc', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%deadstemc_patch) 

       call restartvar(ncid=ncid, flag=flag, varname='deadstemc_storage', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%deadstemc_storage_patch) 

       call restartvar(ncid=ncid, flag=flag, varname='deadstemc_xfer', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%deadstemc_xfer_patch) 

       call restartvar(ncid=ncid, flag=flag, varname='livecrootc', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%livecrootc_patch) 

       call restartvar(ncid=ncid, flag=flag, varname='livecrootc_storage', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%livecrootc_storage_patch) 

       call restartvar(ncid=ncid, flag=flag, varname='livecrootc_xfer', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%livecrootc_xfer_patch) 

       call restartvar(ncid=ncid, flag=flag, varname='deadcrootc', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%deadcrootc_patch) 

       call restartvar(ncid=ncid, flag=flag, varname='deadcrootc_storage', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%deadcrootc_storage_patch) 

       call restartvar(ncid=ncid, flag=flag, varname='deadcrootc_xfer', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%deadcrootc_xfer_patch) 

       call restartvar(ncid=ncid, flag=flag, varname='gresp_storage', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%gresp_storage_patch) 

       call restartvar(ncid=ncid, flag=flag, varname='gresp_xfer', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%gresp_xfer_patch) 

       call restartvar(ncid=ncid, flag=flag, varname='cpool', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%cpool_patch) 

       call restartvar(ncid=ncid, flag=flag, varname='xsmrpool', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%xsmrpool_patch) 

       call restartvar(ncid=ncid, flag=flag, varname='pft_ctrunc', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%ctrunc_patch) 

       call restartvar(ncid=ncid, flag=flag, varname='leafcmax', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%leafcmax_patch)

       if (flag == 'read') then
          call restartvar(ncid=ncid, flag=flag, varname='spinup_state', xtype=ncd_int, &
            long_name='Spinup state of the model that wrote this restart file: ' &
            // ' 0 = normal model mode, 1 = AD spinup, 2 = AAD spinup', units='', &
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
          if ( masterproc ) write(iulog, *) 'exit_spinup ',exit_spinup,' restart_file_spinup_state ',restart_file_spinup_state
          if (spinup_state <= 1 .and. restart_file_spinup_state == 2 ) then
             if ( masterproc ) write(iulog,*) ' CNRest: taking Dead wood C pools out of AD spinup mode'
             exit_spinup = .true.
             if ( masterproc ) write(iulog, *) 'Multiplying stemc and crootc by 10 for exit spinup'
             do i = bounds%begp,bounds%endp
                this%deadstemc_patch(i) = this%deadstemc_patch(i) * 10._r8
                this%deadcrootc_patch(i) = this%deadcrootc_patch(i) * 10._r8
             end do
          else if (spinup_state == 2 .and. restart_file_spinup_state <= 1 )then
             if (spinup_state == 2 .and. restart_file_spinup_state <= 1 )then
                if ( masterproc ) write(iulog,*) ' CNRest: taking Dead wood C pools into AD spinup mode'
                enter_spinup = .true.
                if ( masterproc ) write(iulog, *) 'Dividing stemc and crootc by 10 for enter spinup '
                do i = bounds%begp,bounds%endp
                   this%deadstemc_patch(i) = this%deadstemc_patch(i) / 10._r8
                   this%deadcrootc_patch(i) = this%deadcrootc_patch(i) / 10._r8
                end do
             end if
          end if
       end if
       !--------------------------------
       ! C12 carbon state variables
       !--------------------------------

       if (carbon_type == 'c12') then
          call restartvar(ncid=ncid, flag=flag, varname='totvegc', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%totvegc_patch) 
          ! totvegc_col needed for resetting soil carbon stocks during AD spinup exit
          call restartvar(ncid=ncid, flag=flag, varname='totvegc_col', xtype=ncd_double,  &
               dim1name='column', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%totvegc_col)
       end if

       !--------------------------------
       ! C13 carbon state variables 
       !--------------------------------

       if ( carbon_type == 'c13')  then
          call restartvar(ncid=ncid, flag=flag, varname='totvegc_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%totvegc_patch) 
          if (flag=='read' .and. .not. readvar) then
             if ( masterproc ) write(iulog,*) 'initializing cnveg_carbonstate_inst%totvegc with atmospheric c13 value'
             do i = bounds%begp,bounds%endp
                if (pftcon%c3psn(patch%itype(i)) == 1._r8) then
                   this%totvegc_patch(i) = c12_cnveg_carbonstate_inst%totvegc_patch(i) * c3_r2
                else
                   this%totvegc_patch(i) = c12_cnveg_carbonstate_inst%totvegc_patch(i) * c4_r2
                endif
             end do
          end if

          call restartvar(ncid=ncid, flag=flag, varname='totvegc_col_13', xtype=ncd_double,  &
               dim1name='column', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%totvegc_col)
          if (flag=='read' .and. .not. readvar) then
             if ( masterproc ) write(iulog,*) 'initializing cnveg_carbonstate_inst%totvegc with atmospheric c13 value'
             do i = bounds%begc,bounds%endc
                if (pftcon%c3psn(patch%itype(i)) == 1._r8) then
                   this%totvegc_col(i) = c12_cnveg_carbonstate_inst%totvegc_col(i) * c3_r2
                else
                   this%totvegc_col(i) = c12_cnveg_carbonstate_inst%totvegc_col(i) * c4_r2
                endif
             end do
          end if

       end if

       !--------------------------------
       ! C14 patch carbon state variables 
       !--------------------------------

       if ( carbon_type == 'c14')  then
          call restartvar(ncid=ncid, flag=flag, varname='totvegc_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%totvegc_patch) 
          if (flag=='read' .and. .not. readvar) then
             if ( masterproc ) write(iulog,*) 'initializing this%totvegc_patch with atmospheric c14 value'
             do i = bounds%begp,bounds%endp
                if (this%totvegc_patch(i) /= spval .and. &
                    .not. isnan(this%totvegc_patch(i)) ) then
                   this%totvegc_patch(i) = c12_cnveg_carbonstate_inst%totvegc_patch(i) * c14ratio
                endif
             end do
          endif

          call restartvar(ncid=ncid, flag=flag, varname='totvegc_col_14', xtype=ncd_double,  &
               dim1name='column', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%totvegc_col)
          if (flag=='read' .and. .not. readvar) then
             if ( masterproc ) write(iulog,*) 'initializing cnveg_carbonstate_inst%totvegc with atmospheric c14 value'
             do i = bounds%begc,bounds%endc
                if (this%totvegc_col(i) /= spval .and. &
                    .not. isnan(this%totvegc_col(i)) ) then
                   this%totvegc_col(i) = c12_cnveg_carbonstate_inst%totvegc_col(i) * c14ratio
                endif
             end do
          end if
       end if


       if (  flag == 'read' .and. (enter_spinup .or. (reseed_dead_plants .and. .not. is_restart())) .and. .not. use_cndv) then
             if ( masterproc ) write(iulog, *) 'Reseeding dead plants for CNVegCarbonState'
             ! If a pft is dead (indicated by totvegc = 0) then we reseed that
             ! pft according to the cold start protocol in the InitCold subroutine.
             ! Thus, the variable totvegc is required to be read before here
             ! so that if it is zero for a given pft, the pft can be reseeded.
             do i = bounds%begp,bounds%endp
                if (this%totvegc_patch(i) .le. 0.0_r8) then
                   !-----------------------------------------------
                   ! initialize patch-level carbon state variables
                   !-----------------------------------------------

                   this%leafcmax_patch(i) = 0._r8

                   l = patch%landunit(i)
                   if (lun%itype(l) == istsoil )then
                      if ( present(num_reseed_patch) ) then
                         num_reseed_patch = num_reseed_patch + 1
                         filter_reseed_patch(num_reseed_patch) = i
                      end if

                      if (patch%itype(i) == noveg) then
                         this%leafc_patch(i)          = 0._r8
                         this%leafc_storage_patch(i)  = 0._r8
                         this%frootc_patch(i)         = 0._r8            
                         this%frootc_storage_patch(i) = 0._r8    
                      else
                         if (pftcon%evergreen(patch%itype(i)) == 1._r8) then
                            this%leafc_patch(i)          = cnvegcstate_const%initial_vegC * ratio     
                            this%leafc_storage_patch(i)  = 0._r8
                            this%frootc_patch(i)         = cnvegcstate_const%initial_vegC * ratio           
                            this%frootc_storage_patch(i) = 0._r8    
                         else
                            this%leafc_patch(i)          = 0._r8
                            this%leafc_storage_patch(i)  = cnvegcstate_const%initial_vegC * ratio   
                            this%frootc_patch(i)         = 0._r8            
                            this%frootc_storage_patch(i) = cnvegcstate_const%initial_vegC * ratio   
                         end if
                      end if
                      this%leafc_xfer_patch(i) = 0._r8
                      this%leafc_storage_xfer_acc_patch(i)  = 0._r8
                      this%storage_cdemand_patch(i)         = 0._r8

                      if (MM_Nuptake_opt .eqv. .false.) then  ! if not running in floating CN ratio option 
                         this%frootc_patch(i)            = 0._r8 
                         this%frootc_storage_patch(i)    = 0._r8 
                      end if     
                      this%frootc_xfer_patch(i)       = 0._r8 

                      this%livestemc_patch(i)         = 0._r8 
                      this%livestemc_storage_patch(i) = 0._r8 
                      this%livestemc_xfer_patch(i)    = 0._r8 

                      if (pftcon%woody(patch%itype(i)) == 1._r8) then
                         this%deadstemc_patch(i) = 0.1_r8 * ratio
                      else
                         this%deadstemc_patch(i) = 0._r8 
                      end if
                      this%deadstemc_storage_patch(i)  = 0._r8 
                      this%deadstemc_xfer_patch(i)     = 0._r8 

                      this%livecrootc_patch(i)         = 0._r8 
                      this%livecrootc_storage_patch(i) = 0._r8 
                      this%livecrootc_xfer_patch(i)    = 0._r8 

                      this%deadcrootc_patch(i)         = 0._r8 
                      this%deadcrootc_storage_patch(i) = 0._r8 
                      this%deadcrootc_xfer_patch(i)    = 0._r8 

                      this%gresp_storage_patch(i)      = 0._r8 
                      this%gresp_xfer_patch(i)         = 0._r8 

                      this%cpool_patch(i)              = 0._r8 
                      this%xsmrpool_patch(i)           = 0._r8 
                      this%ctrunc_patch(i)             = 0._r8 
                      this%dispvegc_patch(i)           = 0._r8 
                      this%storvegc_patch(i)           = 0._r8 
                      this%woodc_patch(i)              = 0._r8
                      this%totc_patch(i)               = 0._r8 

                      if ( use_crop )then
                         this%grainc_patch(i)         = 0._r8 
                         this%grainc_storage_patch(i) = 0._r8 
                         this%grainc_xfer_patch(i)    = 0._r8 
                         this%cropseedc_deficit_patch(i)  = 0._r8
                      end if

                      ! calculate totvegc explicitly so that it is available for the isotope 
                      ! code on the first time step.

                      this%totvegc_patch(i) = &
                           this%leafc_patch(i)              + &
                           this%leafc_storage_patch(i)      + &
                           this%leafc_xfer_patch(i)         + &
                           this%frootc_patch(i)             + &
                           this%frootc_storage_patch(i)     + &
                           this%frootc_xfer_patch(i)        + &
                           this%livestemc_patch(i)          + &
                           this%livestemc_storage_patch(i)  + &
                           this%livestemc_xfer_patch(i)     + &
                           this%deadstemc_patch(i)          + &
                           this%deadstemc_storage_patch(i)  + &
                           this%deadstemc_xfer_patch(i)     + &
                           this%livecrootc_patch(i)         + &
                           this%livecrootc_storage_patch(i) + &
                           this%livecrootc_xfer_patch(i)    + &
                           this%deadcrootc_patch(i)         + &
                           this%deadcrootc_storage_patch(i) + &
                           this%deadcrootc_xfer_patch(i)    + &
                           this%gresp_storage_patch(i)      + &
                           this%gresp_xfer_patch(i)         + &
                           this%cpool_patch(i)

                      if ( use_crop )then
                         this%totvegc_patch(i) =         &
                              this%totvegc_patch(i)    + &
                              this%grainc_patch(i)         + &
                              this%grainc_storage_patch(i) + &
                              this%grainc_xfer_patch(i)
                      end if

                   endif
                end if
             end do
             if ( .not. is_restart() .and. get_nstep() == 1 ) then

                do p = bounds%begp,bounds%endp
                  if (this%leafc_patch(p) .lt. 0.01_r8) then
                   if (pftcon%c3psn(patch%itype(p)) == 1._r8) then
                      this%grainc_patch(p)         = c12_cnveg_carbonstate_inst%grainc_patch(p)         * c3_r2
                      this%grainc_storage_patch(p) = c12_cnveg_carbonstate_inst%grainc_storage_patch(p) * c3_r2
                      this%grainc_xfer_patch(p)    = c12_cnveg_carbonstate_inst%grainc_xfer_patch(p)    * c3_r2
                      this%dispvegc_patch(p)       = c12_cnveg_carbonstate_inst%dispvegc_patch(p)       * c3_r2
                      this%storvegc_patch(p)       = c12_cnveg_carbonstate_inst%storvegc_patch(p)       * c3_r2
                      this%totvegc_patch(p)        = c12_cnveg_carbonstate_inst%totvegc_patch(p)        * c3_r2
                      this%totc_patch(p)           = c12_cnveg_carbonstate_inst%totc_patch(p)           * c3_r2
                      this%woodc_patch(p)          = c12_cnveg_carbonstate_inst%woodc_patch(p)          * c3_r2
                   else
                      this%grainc_patch(p)         = c12_cnveg_carbonstate_inst%grainc_patch(p)         * c4_r2
                      this%grainc_storage_patch(p) = c12_cnveg_carbonstate_inst%grainc_storage_patch(p) * c4_r2
                      this%grainc_xfer_patch(p)    = c12_cnveg_carbonstate_inst%grainc_xfer_patch(p)    * c4_r2
                      this%dispvegc_patch(p)       = c12_cnveg_carbonstate_inst%dispvegc_patch(p)       * c4_r2
                      this%storvegc_patch(p)       = c12_cnveg_carbonstate_inst%storvegc_patch(p)       * c4_r2
                      this%totvegc_patch(p)        = c12_cnveg_carbonstate_inst%totvegc_patch(p)        * c4_r2
                      this%totc_patch(p)           = c12_cnveg_carbonstate_inst%totc_patch(p)           * c4_r2
                      this%woodc_patch(p)          = c12_cnveg_carbonstate_inst%woodc_patch(p)          * c4_r2
                   end if
                  end if
                end do
             end if
             if ( present(num_reseed_patch) ) then
                call shr_mpi_sum( num_reseed_patch, total_num_reseed_patch, mpicom )
                if ( masterproc ) write(iulog,*) 'Total num_reseed, over all tasks = ', total_num_reseed_patch
             end if
       end if

    end if

    !--------------------------------
    ! C13 patch carbon state variables 
    !--------------------------------

    if ( carbon_type == 'c13')  then
       call restartvar(ncid=ncid, flag=flag, varname='leafc_13', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%leafc_patch)
       if (flag=='read' .and. .not. readvar) then
          if ( masterproc ) write(iulog,*) 'initializing this%leafc with atmospheric c13 value'
          do i = bounds%begp,bounds%endp
             if (pftcon%c3psn(patch%itype(i)) == 1._r8) then
                this%leafc_patch(i) = c12_cnveg_carbonstate_inst%leafc_patch(i) * c3_r2
             else
                this%leafc_patch(i) = c12_cnveg_carbonstate_inst%leafc_patch(i) * c4_r2
             endif
          end do
       end if

       call restartvar(ncid=ncid, flag=flag, varname='leafc_storage_13', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%leafc_storage_patch) 
       if (flag=='read' .and. .not. readvar) then
          if ( masterproc ) write(iulog,*) 'initializing this%leafc_storage with atmospheric c13 value'
          do i = bounds%begp,bounds%endp
             if (pftcon%c3psn(patch%itype(i)) == 1._r8) then
                this%leafc_storage_patch(i) = c12_cnveg_carbonstate_inst%leafc_storage_patch(i) * c3_r2
             else
                this%leafc_storage_patch(i) = c12_cnveg_carbonstate_inst%leafc_storage_patch(i) * c4_r2
                this%leafc_storage_patch(i) = c12_cnveg_carbonstate_inst%leafc_storage_patch(i) * c4_r2
             endif
          end do
       end if

       call restartvar(ncid=ncid, flag=flag, varname='leafc_xfer_13', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%leafc_xfer_patch) 
       if (flag=='read' .and. .not. readvar) then
          if ( masterproc ) write(iulog,*) 'initializing this%leafc_xfer with atmospheric c13 value'
          do i = bounds%begp,bounds%endp
             if (pftcon%c3psn(patch%itype(i)) == 1._r8) then
                this%leafc_xfer_patch(i) = c12_cnveg_carbonstate_inst%leafc_xfer_patch(i) * c3_r2
             else
                this%leafc_xfer_patch(i) = c12_cnveg_carbonstate_inst%leafc_xfer_patch(i) * c4_r2
             endif
          end do
       end if

       call restartvar(ncid=ncid, flag=flag, varname='frootc_13', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%frootc_patch) 
       if (flag=='read' .and. .not. readvar) then
          if ( masterproc ) write(iulog,*) 'initializing this%frootc with atmospheric c13 value'
          do i = bounds%begp,bounds%endp
             if (pftcon%c3psn(patch%itype(i)) == 1._r8) then
                this%frootc_patch(i) = c12_cnveg_carbonstate_inst%frootc_patch(i) * c3_r2
             else
                this%frootc_patch(i) = c12_cnveg_carbonstate_inst%frootc_patch(i) * c4_r2
             endif
          end do
       end if

       call restartvar(ncid=ncid, flag=flag, varname='frootc_storage_13', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%frootc_storage_patch) 
       if (flag=='read' .and. .not. readvar) then
          if ( masterproc ) write(iulog,*) 'initializing this%frootc_storage with atmospheric c13 value'
          do i = bounds%begp,bounds%endp
             if (pftcon%c3psn(patch%itype(i)) == 1._r8) then
                this%frootc_storage_patch(i) = c12_cnveg_carbonstate_inst%frootc_storage_patch(i) * c3_r2
             else
                this%frootc_storage_patch(i) = c12_cnveg_carbonstate_inst%frootc_storage_patch(i) * c4_r2
             endif
          end do
       end if

       call restartvar(ncid=ncid, flag=flag, varname='frootc_xfer_13', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%frootc_xfer_patch) 
       if (flag=='read' .and. .not. readvar) then
          if ( masterproc ) write(iulog,*) 'initializing this%frootc_xfer with atmospheric c13 value'
          do i = bounds%begp,bounds%endp
             if (pftcon%c3psn(patch%itype(i)) == 1._r8) then
                this%frootc_xfer_patch(i) = c12_cnveg_carbonstate_inst%frootc_xfer_patch(i) * c3_r2
             else
                this%frootc_xfer_patch(i) = c12_cnveg_carbonstate_inst%frootc_xfer_patch(i) * c4_r2
             endif
          end do
       end if

       call restartvar(ncid=ncid, flag=flag, varname='livestemc_13', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%livestemc_patch) 
       if (flag=='read' .and. .not. readvar) then
          if ( masterproc ) write(iulog,*) 'initializing this%livestemc with atmospheric c13 value'
          do i = bounds%begp,bounds%endp
             if (pftcon%c3psn(patch%itype(i)) == 1._r8) then
                this%livestemc_patch(i) = c12_cnveg_carbonstate_inst%livestemc_patch(i) * c3_r2
             else
                this%livestemc_patch(i) = c12_cnveg_carbonstate_inst%livestemc_patch(i) * c4_r2
             endif
          end do
       end if

       call restartvar(ncid=ncid, flag=flag, varname='livestemc_storage_13', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%livestemc_storage_patch) 
       if (flag=='read' .and. .not. readvar) then
          if ( masterproc ) write(iulog,*) 'initializing this%livestemc_storage with atmospheric c13 value'
          do i = bounds%begp,bounds%endp
             if (pftcon%c3psn(patch%itype(i)) == 1._r8) then
                this%livestemc_storage_patch(i) = c12_cnveg_carbonstate_inst%livestemc_storage_patch(i) * c3_r2
             else
                this%livestemc_storage_patch(i) = c12_cnveg_carbonstate_inst%livestemc_storage_patch(i) * c4_r2
             endif
          end do
       end if

       call restartvar(ncid=ncid, flag=flag, varname='livestemc_xfer_13', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%livestemc_xfer_patch) 
       if (flag=='read' .and. .not. readvar) then
          if ( masterproc ) write(iulog,*) 'initializing this%livestemc_xfer with atmospheric c13 value'
          do i = bounds%begp,bounds%endp
             if (pftcon%c3psn(patch%itype(i)) == 1._r8) then
                this%livestemc_xfer_patch(i) = c12_cnveg_carbonstate_inst%livestemc_xfer_patch(i) * c3_r2
             else
                this%livestemc_xfer_patch(i) = c12_cnveg_carbonstate_inst%livestemc_xfer_patch(i) * c4_r2
             endif
          end do
       end if

       call restartvar(ncid=ncid, flag=flag, varname='deadstemc_13', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%deadstemc_patch) 
       if (flag=='read' .and. .not. readvar) then
          if ( masterproc ) write(iulog,*) 'initializing this%deadstemc with atmospheric c13 value'
          do i = bounds%begp,bounds%endp
             if (pftcon%c3psn(patch%itype(i)) == 1._r8) then
                this%deadstemc_patch(i) = c12_cnveg_carbonstate_inst%deadstemc_patch(i) * c3_r2
             else
                this%deadstemc_patch(i) = c12_cnveg_carbonstate_inst%deadstemc_patch(i) * c4_r2
             endif
          end do
       end if

       call restartvar(ncid=ncid, flag=flag, varname='deadstemc_storage_13', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%deadstemc_storage_patch) 
       if (flag=='read' .and. .not. readvar) then
          if ( masterproc ) write(iulog,*) 'initializing this%deadstemc_storage with atmospheric c13 value'
          do i = bounds%begp,bounds%endp
             if (pftcon%c3psn(patch%itype(i)) == 1._r8) then
                this%deadstemc_storage_patch(i) = c12_cnveg_carbonstate_inst%deadstemc_storage_patch(i) * c3_r2
             else
                this%deadstemc_storage_patch(i) = c12_cnveg_carbonstate_inst%deadstemc_storage_patch(i) * c4_r2
             endif
          end do
       end if

       call restartvar(ncid=ncid, flag=flag, varname='deadstemc_xfer_13', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%deadstemc_xfer_patch) 
       if (flag=='read' .and. .not. readvar) then
          if ( masterproc ) write(iulog,*) 'initializing this%deadstemc_xfer with atmospheric c13 value'
          do i = bounds%begp,bounds%endp
             if (pftcon%c3psn(patch%itype(i)) == 1._r8) then
                this%deadstemc_xfer_patch(i) = c12_cnveg_carbonstate_inst%deadstemc_xfer_patch(i) * c3_r2
             else
                this%deadstemc_xfer_patch(i) = c12_cnveg_carbonstate_inst%deadstemc_xfer_patch(i) * c4_r2
             endif
          end do
       end if

       call restartvar(ncid=ncid, flag=flag, varname='livecrootc_13', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%livecrootc_patch) 
       if (flag=='read' .and. .not. readvar) then
          if ( masterproc ) write(iulog,*) 'initializing this%livecrootc with atmospheric c13 value'
          do i = bounds%begp,bounds%endp
             if (pftcon%c3psn(patch%itype(i)) == 1._r8) then
                this%livecrootc_patch(i) = c12_cnveg_carbonstate_inst%livecrootc_patch(i) * c3_r2
             else
                this%livecrootc_patch(i) = c12_cnveg_carbonstate_inst%livecrootc_patch(i) * c4_r2
             endif
          end do
       end if

       call restartvar(ncid=ncid, flag=flag, varname='livecrootc_storage_13', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%livecrootc_storage_patch) 
       if (flag=='read' .and. .not. readvar) then
          if ( masterproc ) write(iulog,*) 'initializing this%livecrootc_storage with atmospheric c13 value'
          do i = bounds%begp,bounds%endp
             if (pftcon%c3psn(patch%itype(i)) == 1._r8) then
                this%livecrootc_storage_patch(i) = c12_cnveg_carbonstate_inst%livecrootc_storage_patch(i) * c3_r2
             else
                this%livecrootc_storage_patch(i) = c12_cnveg_carbonstate_inst%livecrootc_storage_patch(i) * c4_r2
             endif
          end do
       end if

       call restartvar(ncid=ncid, flag=flag, varname='livecrootc_xfer_13', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%livecrootc_xfer_patch) 
       if (flag=='read' .and. .not. readvar) then
          if ( masterproc ) write(iulog,*) 'initializing this%livecrootc_xfer with atmospheric c13 value'
          do i = bounds%begp,bounds%endp
             if (pftcon%c3psn(patch%itype(i)) == 1._r8) then
                this%livecrootc_xfer_patch(i) = c12_cnveg_carbonstate_inst%livecrootc_xfer_patch(i) * c3_r2
             else
                this%livecrootc_xfer_patch(i) = c12_cnveg_carbonstate_inst%livecrootc_xfer_patch(i) * c4_r2
             endif
          end do
       end if

       call restartvar(ncid=ncid, flag=flag, varname='deadcrootc_13', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%deadcrootc_patch) 
       if (flag=='read' .and. .not. readvar) then
          if ( masterproc ) write(iulog,*) 'initializing this%deadcrootc with atmospheric c13 value'
          do i = bounds%begp,bounds%endp
             if (pftcon%c3psn(patch%itype(i)) == 1._r8) then
                this%deadcrootc_patch(i) = c12_cnveg_carbonstate_inst%deadcrootc_patch(i) * c3_r2
             else
                this%deadcrootc_patch(i) = c12_cnveg_carbonstate_inst%deadcrootc_patch(i) * c4_r2
             endif
          end do
       end if

       call restartvar(ncid=ncid, flag=flag, varname='deadcrootc_storage_13', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%deadcrootc_storage_patch) 
       if (flag=='read' .and. .not. readvar) then
          if ( masterproc ) write(iulog,*) 'initializing this%deadcrootc_storage with atmospheric c13 value'
          do i = bounds%begp,bounds%endp
             if (pftcon%c3psn(patch%itype(i)) == 1._r8) then
                this%deadcrootc_storage_patch(i) = c12_cnveg_carbonstate_inst%deadcrootc_storage_patch(i) * c3_r2
             else
                this%deadcrootc_storage_patch(i) = c12_cnveg_carbonstate_inst%deadcrootc_storage_patch(i) * c4_r2
             endif
          end do
       end if

       call restartvar(ncid=ncid, flag=flag, varname='deadcrootc_xfer_13', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%deadcrootc_xfer_patch) 
       if (flag=='read' .and. .not. readvar) then
          if ( masterproc ) write(iulog,*) 'initializing this%deadcrootc_xfer with atmospheric c13 value'
          do i = bounds%begp,bounds%endp
             if (pftcon%c3psn(patch%itype(i)) == 1._r8) then
                this%deadcrootc_xfer_patch(i) = c12_cnveg_carbonstate_inst%deadcrootc_xfer_patch(i) * c3_r2
             else
                this%deadcrootc_xfer_patch(i) = c12_cnveg_carbonstate_inst%deadcrootc_xfer_patch(i) * c4_r2
             endif
          end do
       end if

       call restartvar(ncid=ncid, flag=flag, varname='gresp_storage_13', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%gresp_storage_patch) 
       if (flag=='read' .and. .not. readvar) then
          if ( masterproc ) write(iulog,*) 'initializing this%gresp_storage with atmospheric c13 value'
          do i = bounds%begp,bounds%endp
             if (pftcon%c3psn(patch%itype(i)) == 1._r8) then
                this%gresp_storage_patch(i) = c12_cnveg_carbonstate_inst%gresp_storage_patch(i) * c3_r2
             else
                this%gresp_storage_patch(i) = c12_cnveg_carbonstate_inst%gresp_storage_patch(i) * c4_r2
             endif
          end do
       end if

       call restartvar(ncid=ncid, flag=flag, varname='gresp_xfer_13', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%gresp_xfer_patch) 
       if (flag=='read' .and. .not. readvar) then
          if ( masterproc ) write(iulog,*) 'initializing this%gresp_xfer with atmospheric c13 value'
          do i = bounds%begp,bounds%endp
             if (pftcon%c3psn(patch%itype(i)) == 1._r8) then
                this%gresp_xfer_patch(i) = c12_cnveg_carbonstate_inst%gresp_xfer_patch(i) * c3_r2
             else
                this%gresp_xfer_patch(i) = c12_cnveg_carbonstate_inst%gresp_xfer_patch(i) * c4_r2
             endif
          end do
       end if

       call restartvar(ncid=ncid, flag=flag, varname='cpool_13', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%cpool_patch) 
       if (flag=='read' .and. .not. readvar) then
          if ( masterproc ) write(iulog,*) 'initializing this%cpool with atmospheric c13 value'
          do i = bounds%begp,bounds%endp
             if (pftcon%c3psn(patch%itype(i)) == 1._r8) then
                this%cpool_patch(i) = c12_cnveg_carbonstate_inst%cpool_patch(i) * c3_r2
             else
                this%cpool_patch(i) = c12_cnveg_carbonstate_inst%cpool_patch(i) * c4_r2
             endif
          end do
       end if

       call restartvar(ncid=ncid, flag=flag, varname='xsmrpool_13', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%xsmrpool_patch) 
       if (flag=='read' .and. .not. readvar) then
          if ( masterproc ) write(iulog,*) 'initializing this%xsmrpool with atmospheric c13 value'
          do i = bounds%begp,bounds%endp
             if (pftcon%c3psn(patch%itype(i)) == 1._r8) then
                this%xsmrpool_patch(i) = c12_cnveg_carbonstate_inst%xsmrpool_patch(i) * c3_r2
             else
                this%xsmrpool_patch(i) = c12_cnveg_carbonstate_inst%xsmrpool_patch(i) * c4_r2
             endif
          end do
       end if

       call restartvar(ncid=ncid, flag=flag, varname='pft_ctrunc_13', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%ctrunc_patch) 
       if (flag=='read' .and. .not. readvar) then
          if ( masterproc ) write(iulog,*) 'initializing this%ctrunc with atmospheric c13 value'
          do i = bounds%begp,bounds%endp
             if (pftcon%c3psn(patch%itype(i)) == 1._r8) then
                this%ctrunc_patch(i) = c12_cnveg_carbonstate_inst%ctrunc_patch(i) * c3_r2
             else
                this%ctrunc_patch(i) = c12_cnveg_carbonstate_inst%ctrunc_patch(i) * c4_r2
             endif
          end do
       end if

    end if

    !--------------------------------
    ! C14 patch carbon state variables 
    !--------------------------------

    if ( carbon_type == 'c14')  then
       call restartvar(ncid=ncid, flag=flag, varname='leafc_14', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%leafc_patch) 
       if (flag=='read' .and. .not. readvar) then
          if ( masterproc ) write(iulog,*) 'initializing this%leafc_patch with atmospheric c14 value'
          do i = bounds%begp,bounds%endp
             if (this%leafc_patch(i) /= spval .and. &
                  .not. isnan(this%leafc_patch(i)) ) then
                this%leafc_patch(i) = c12_cnveg_carbonstate_inst%leafc_patch(i) * c14ratio
             endif
          end do
       end if

       call restartvar(ncid=ncid, flag=flag, varname='leafc_storage_14', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%leafc_storage_patch) 
       if (flag=='read' .and. .not. readvar) then
          if ( masterproc ) write(iulog,*) 'initializing this%leafc_storage_patch with atmospheric c14 value'
          do i = bounds%begp,bounds%endp
             if (this%leafc_storage_patch(i) /= spval .and. &
                  .not. isnan(this%leafc_storage_patch(i)) ) then
                this%leafc_storage_patch(i) = c12_cnveg_carbonstate_inst%leafc_storage_patch(i) * c14ratio
             endif
          end do
       end if

       call restartvar(ncid=ncid, flag=flag, varname='leafc_xfer_14', xtype=ncd_double,  &
            dim1name='pft',    long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%leafc_xfer_patch) 
       if (flag=='read' .and. .not. readvar) then
          if ( masterproc ) write(iulog,*) 'initializing this%leafc_xfer_patch with atmospheric c14 value'
          do i = bounds%begp,bounds%endp
             if (this%leafc_xfer_patch(i) /= spval .and. .not. isnan(this%leafc_xfer_patch(i)) ) then
                this%leafc_xfer_patch(i) = c12_cnveg_carbonstate_inst%leafc_xfer_patch(i) * c14ratio
             endif
          end do
       end if

       call restartvar(ncid=ncid, flag=flag, varname='frootc_14', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%frootc_patch) 
       if (flag=='read' .and. .not. readvar) then
          if ( masterproc ) write(iulog,*) 'initializing this%frootc_patch with atmospheric c14 value'
          do i = bounds%begp,bounds%endp
             if (this%frootc_patch(i) /= spval .and. &
                  .not. isnan(this%frootc_patch(i)) ) then
                this%frootc_patch(i) = c12_cnveg_carbonstate_inst%frootc_patch(i) * c14ratio
             endif
          end do
       end if

       call restartvar(ncid=ncid, flag=flag, varname='frootc_storage_14', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%frootc_storage_patch) 
       if (flag=='read' .and. .not. readvar) then
          if ( masterproc ) write(iulog,*) 'initializing this%frootc_storage_patch with atmospheric c14 value'
          do i = bounds%begp,bounds%endp
             if (this%frootc_storage_patch(i) /= spval .and. &
                  .not. isnan(this%frootc_storage_patch(i)) ) then
                this%frootc_storage_patch(i) = c12_cnveg_carbonstate_inst%frootc_storage_patch(i) * c14ratio
             endif
          end do
       end if

       call restartvar(ncid=ncid, flag=flag, varname='frootc_xfer_14', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%frootc_xfer_patch) 
       if (flag=='read' .and. .not. readvar) then
          if ( masterproc ) write(iulog,*) 'initializing this%frootc_xfer_patch with atmospheric c14 value'
          do i = bounds%begp,bounds%endp
             if (this%frootc_xfer_patch(i) /= spval .and. &
                  .not. isnan(this%frootc_xfer_patch(i)) ) then
                this%frootc_xfer_patch(i) = c12_cnveg_carbonstate_inst%frootc_xfer_patch(i) * c14ratio
             endif
          end do
       end if

       call restartvar(ncid=ncid, flag=flag, varname='livestemc_14', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%livestemc_patch) 
       if (flag=='read' .and. .not. readvar) then
          if ( masterproc ) write(iulog,*) 'initializing this%livestemc_patch with atmospheric c14 value'
          do i = bounds%begp,bounds%endp
             if (this%livestemc_patch(i) /= spval .and. .not. isnan(this%livestemc_patch(i)) ) then
                this%livestemc_patch(i) = c12_cnveg_carbonstate_inst%livestemc_patch(i) * c14ratio
             endif
          end do
       end if

       call restartvar(ncid=ncid, flag=flag, varname='livestemc_storage_14', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%livestemc_storage_patch) 
       if (flag=='read' .and. .not. readvar) then
          if ( masterproc ) write(iulog,*) 'initializing this%livestemc_storage_patch with atmospheric c14 value'
          do i = bounds%begp,bounds%endp
             if (this%livestemc_storage_patch(i) /= spval .and. .not. isnan(this%livestemc_storage_patch(i)) ) then
                this%livestemc_storage_patch(i) = c12_cnveg_carbonstate_inst%livestemc_storage_patch(i) * c14ratio
             endif
          end do
       end if

       call restartvar(ncid=ncid, flag=flag, varname='livestemc_xfer_14', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%livestemc_xfer_patch) 
       if (flag=='read' .and. .not. readvar) then
          if ( masterproc ) write(iulog,*) 'initializing this%livestemc_xfer_patch with atmospheric c14 value'
          do i = bounds%begp,bounds%endp
             if (this%livestemc_xfer_patch(i) /= spval .and. .not. isnan(this%livestemc_xfer_patch(i)) ) then
                this%livestemc_xfer_patch(i) = c12_cnveg_carbonstate_inst%livestemc_xfer_patch(i) * c14ratio
             endif
          end do
       end if

       call restartvar(ncid=ncid, flag=flag, varname='deadstemc_14', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%deadstemc_patch) 
       if (flag=='read' .and. .not. readvar) then
          if ( masterproc ) write(iulog,*) 'initializing this%deadstemc_patch with atmospheric c14 value'
          do i = bounds%begp,bounds%endp
             if (this%deadstemc_patch(i) /= spval .and. .not. isnan(this%deadstemc_patch(i)) ) then
                this%deadstemc_patch(i) = c12_cnveg_carbonstate_inst%deadstemc_patch(i) * c14ratio
             endif
          end do
       end if

       call restartvar(ncid=ncid, flag=flag, varname='deadstemc_storage_14', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%deadstemc_storage_patch) 
       if (flag=='read' .and. .not. readvar) then
          if ( masterproc ) write(iulog,*) 'initializing this%deadstemc_storage_patch with atmospheric c14 value'
          do i = bounds%begp,bounds%endp
             if (this%deadstemc_storage_patch(i) /= spval .and. .not. isnan(this%deadstemc_storage_patch(i)) ) then
                this%deadstemc_storage_patch(i) = c12_cnveg_carbonstate_inst%deadstemc_storage_patch(i) * c14ratio
             endif
          end do
       end if

       call restartvar(ncid=ncid, flag=flag, varname='deadstemc_xfer_14', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%deadstemc_xfer_patch) 
       if (flag=='read' .and. .not. readvar) then
          if ( masterproc ) write(iulog,*) 'initializing this%deadstemc_xfer_patch with atmospheric c14 value'
          do i = bounds%begp,bounds%endp
             if (this%deadstemc_xfer_patch(i) /= spval .and. .not. isnan(this%deadstemc_xfer_patch(i)) ) then
                this%deadstemc_xfer_patch(i) = c12_cnveg_carbonstate_inst%deadstemc_xfer_patch(i) * c14ratio
             endif
          end do
       end if

       call restartvar(ncid=ncid, flag=flag, varname='livecrootc_14', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%livecrootc_patch) 
       if (flag=='read' .and. .not. readvar) then
          if ( masterproc ) write(iulog,*) 'initializing this%livecrootc_patch with atmospheric c14 value'
          do i = bounds%begp,bounds%endp
             if (this%livecrootc_patch(i) /= spval .and. .not. isnan(this%livecrootc_patch(i)) ) then
                this%livecrootc_patch(i) = c12_cnveg_carbonstate_inst%livecrootc_patch(i) * c14ratio
             endif
          end do
       end if

       call restartvar(ncid=ncid, flag=flag, varname='livecrootc_storage_14', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%livecrootc_storage_patch) 
       if (flag=='read' .and. .not. readvar) then
          if ( masterproc ) write(iulog,*) 'initializing this%livecrootc_storage_patch with atmospheric c14 value'
          do i = bounds%begp,bounds%endp
             if (this%livecrootc_storage_patch(i) /= spval .and. .not. isnan(this%livecrootc_storage_patch(i)) ) then
                this%livecrootc_storage_patch(i) = c12_cnveg_carbonstate_inst%livecrootc_storage_patch(i) * c14ratio
             endif
          end do
       end if

       call restartvar(ncid=ncid, flag=flag, varname='livecrootc_xfer_14', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%livecrootc_xfer_patch) 
       if (flag=='read' .and. .not. readvar) then
          if ( masterproc ) write(iulog,*) 'initializing this%livecrootc_xfer_patch with atmospheric c14 value'
          do i = bounds%begp,bounds%endp
             if (this%livecrootc_xfer_patch(i) /= spval .and. .not. isnan(this%livecrootc_xfer_patch(i)) ) then
                this%livecrootc_xfer_patch(i) = c12_cnveg_carbonstate_inst%livecrootc_xfer_patch(i) * c14ratio
             endif
          end do
       end if

       call restartvar(ncid=ncid, flag=flag, varname='deadcrootc_14', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%deadcrootc_patch) 
       if (flag=='read' .and. .not. readvar) then
          if ( masterproc ) write(iulog,*) 'initializing this%deadcrootc_patch with atmospheric c14 value'
          do i = bounds%begp,bounds%endp
             if (this%deadcrootc_patch(i) /= spval .and. .not. isnan(this%deadcrootc_patch(i)) ) then
                this%deadcrootc_patch(i) = c12_cnveg_carbonstate_inst%deadcrootc_patch(i) * c14ratio
             endif
          end do
       end if

       call restartvar(ncid=ncid, flag=flag, varname='deadcrootc_storage_14', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%deadcrootc_storage_patch) 
       if (flag=='read' .and. .not. readvar) then
          if ( masterproc ) write(iulog,*) 'initializing this%deadcrootc_storage_patch with atmospheric c14 value'
          do i = bounds%begp,bounds%endp
             if (this%deadcrootc_storage_patch(i) /= spval .and. .not. isnan(this%deadcrootc_storage_patch(i)) ) then
                this%deadcrootc_storage_patch(i) = c12_cnveg_carbonstate_inst%deadcrootc_storage_patch(i) * c14ratio
             endif
          end do
       end if

       call restartvar(ncid=ncid, flag=flag, varname='deadcrootc_xfer_14', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%deadcrootc_xfer_patch) 
       if (flag=='read' .and. .not. readvar) then
          if ( masterproc ) write(iulog,*) 'initializing this%deadcrootc_xfer_patch with atmospheric c14 value'
          do i = bounds%begp,bounds%endp
             if (this%deadcrootc_xfer_patch(i) /= spval .and. .not. isnan(this%deadcrootc_xfer_patch(i)) ) then
                this%deadcrootc_xfer_patch(i) = c12_cnveg_carbonstate_inst%deadcrootc_xfer_patch(i) * c14ratio
             endif
          end do
       end if

       call restartvar(ncid=ncid, flag=flag, varname='gresp_storage_14', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%gresp_storage_patch) 
       if (flag=='read' .and. .not. readvar) then
          if ( masterproc ) write(iulog,*) 'initializing this%gresp_storage_patch with atmospheric c14 value'
          do i = bounds%begp,bounds%endp
             if (this%gresp_storage_patch(i) /= spval .and. .not. isnan(this%gresp_storage_patch(i)) ) then
                this%gresp_storage_patch(i) = c12_cnveg_carbonstate_inst%gresp_storage_patch(i) * c14ratio
             endif
          end do
       end if

       call restartvar(ncid=ncid, flag=flag, varname='gresp_xfer_14', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%gresp_xfer_patch) 
       if (flag=='read' .and. .not. readvar) then
          if ( masterproc ) write(iulog,*) 'initializing this%gresp_xfer_patch with atmospheric c14 value'
          do i = bounds%begp,bounds%endp
             if (this%gresp_xfer_patch(i) /= spval .and. .not. isnan(this%gresp_xfer_patch(i)) ) then
                this%gresp_xfer_patch(i) = c12_cnveg_carbonstate_inst%gresp_xfer_patch(i) * c14ratio
             endif
          end do
       end if

       call restartvar(ncid=ncid, flag=flag, varname='cpool_14', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%cpool_patch) 
       if (flag=='read' .and. .not. readvar) then
          if ( masterproc ) write(iulog,*) 'initializing this%cpool_patch with atmospheric c14 value'
          do i = bounds%begp,bounds%endp
             if (this%cpool_patch(i) /= spval .and. .not. isnan(this%cpool_patch(i)) ) then
                this%cpool_patch(i) = c12_cnveg_carbonstate_inst%cpool_patch(i) * c14ratio
             endif
          end do
       end if

       call restartvar(ncid=ncid, flag=flag, varname='xsmrpool_14', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%xsmrpool_patch) 
       if (flag=='read' .and. .not. readvar) then
          if ( masterproc ) write(iulog,*) 'initializing this%xsmrpool_patch with atmospheric c14 value'
          do i = bounds%begp,bounds%endp
             if (this%xsmrpool_patch(i) /= spval .and. .not. isnan(this%xsmrpool_patch(i)) ) then
                this%xsmrpool_patch(i) = c12_cnveg_carbonstate_inst%xsmrpool_patch(i) * c14ratio
             endif
          end do
       end if

       call restartvar(ncid=ncid, flag=flag, varname='pft_ctrunc_14', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%ctrunc_patch) 
       if (flag=='read' .and. .not. readvar) then
          if ( masterproc ) write(iulog,*) 'initializing this%ctrunc_patch with atmospheric c14 value'
          do i = bounds%begp,bounds%endp
             if (this%ctrunc_patch(i) /= spval .and. .not. isnan(this%ctrunc_patch(i)) ) then
                this%ctrunc_patch(i) = c12_cnveg_carbonstate_inst%ctrunc_patch(i) * c14ratio
             endif
          end do
       end if

    end if

    !--------------------------------
    ! patch prognostic crop variables
    !--------------------------------

    if (use_crop) then
       if (carbon_type == 'c12') then
          call restartvar(ncid=ncid, flag=flag,  varname='grainc', xtype=ncd_double,  &
               dim1name='pft', long_name='grain C', units='gC/m2', &
               interpinic_flag='interp', readvar=readvar, data=this%grainc_patch)

          call restartvar(ncid=ncid, flag=flag,  varname='grainc_storage', xtype=ncd_double,  &
               dim1name='pft', long_name='grain C storage', units='gC/m2', &
               interpinic_flag='interp', readvar=readvar, data=this%grainc_storage_patch)

          call restartvar(ncid=ncid, flag=flag,  varname='grainc_xfer', xtype=ncd_double,  &
               dim1name='pft', long_name='grain C transfer', units='gC/m2', &
               interpinic_flag='interp', readvar=readvar, data=this%grainc_xfer_patch)

          call restartvar(ncid=ncid, flag=flag, varname='cropseedc_deficit', xtype=ncd_double,  &
               dim1name='pft', long_name='pool for seeding new crop growth', units='gC/m2', &
               interpinic_flag='interp', readvar=readvar, data=this%cropseedc_deficit_patch)
       end if

       if (carbon_type == 'c13') then
          call restartvar(ncid=ncid, flag=flag, varname='grainc_13', xtype=ncd_double,  &
               dim1name='pft', long_name='c13 grain C', units='gC13/m2', &
               interpinic_flag='interp', readvar=readvar, data=this%grainc_patch)
          if (flag=='read' .and. .not. readvar) then
             call set_missing_from_template( &
                  my_var = this%grainc_patch, &
                  template_var = c12_cnveg_carbonstate_inst%grainc_patch, &
                  multiplier = c3_r2)
          end if

          call restartvar(ncid=ncid, flag=flag, varname='grainc_13_storage', xtype=ncd_double,  &
               dim1name='pft', long_name='c13 grain C storage', units='gC13/m2', &
               interpinic_flag='interp', readvar=readvar, data=this%grainc_storage_patch)
          if (flag=='read' .and. .not. readvar) then
             call set_missing_from_template( &
                  my_var = this%grainc_storage_patch, &
                  template_var = c12_cnveg_carbonstate_inst%grainc_storage_patch, &
                  multiplier = c3_r2)
          end if

          call restartvar(ncid=ncid, flag=flag, varname='grainc_13_xfer', xtype=ncd_double,  &
               dim1name='pft', long_name='c13 grain C transfer', units='gC13/m2', &
               interpinic_flag='interp', readvar=readvar, data=this%grainc_xfer_patch)
          if (flag=='read' .and. .not. readvar) then
             call set_missing_from_template( &
                  my_var = this%grainc_xfer_patch, &
                  template_var = c12_cnveg_carbonstate_inst%grainc_xfer_patch, &
                  multiplier = c3_r2)
          end if

          call restartvar(ncid=ncid, flag=flag, varname='cropseedc_13_deficit', xtype=ncd_double,  &
               dim1name='pft', long_name='pool for seeding new crop growth', units='gC13/m2', &
               interpinic_flag='interp', readvar=readvar, data=this%cropseedc_deficit_patch)
          if (flag=='read' .and. .not. readvar) then
             call set_missing_from_template( &
                  my_var = this%cropseedc_deficit_patch, &
                  template_var = c12_cnveg_carbonstate_inst%cropseedc_deficit_patch, &
                  multiplier = c3_r2)
          end if
       end if

       if ( carbon_type == 'c14' ) then

          call restartvar(ncid=ncid, flag=flag, varname='grainc_14', xtype=ncd_double,  &
               dim1name='pft', long_name='c14 grain C', units='gC14/m2', &
               interpinic_flag='interp', readvar=readvar, data=this%grainc_patch)
          if (flag=='read' .and. .not. readvar) then
             call set_missing_from_template( &
                  my_var = this%grainc_patch, &
                  template_var = c12_cnveg_carbonstate_inst%grainc_patch, &
                  multiplier = c3_r2)
          end if

          call restartvar(ncid=ncid, flag=flag, varname='grainc_14_storage', xtype=ncd_double,  &
               dim1name='pft', long_name='c14 grain C storage', units='gC14/m2', &
               interpinic_flag='interp', readvar=readvar, data=this%grainc_storage_patch)
          if (flag=='read' .and. .not. readvar) then
             call set_missing_from_template( &
                  my_var = this%grainc_storage_patch, &
                  template_var = c12_cnveg_carbonstate_inst%grainc_storage_patch, &
                  multiplier = c3_r2)
          end if

          call restartvar(ncid=ncid, flag=flag, varname='grainc_14_xfer', xtype=ncd_double,  &
               dim1name='pft', long_name='c14 grain C transfer', units='gC14/m2', &
               interpinic_flag='interp', readvar=readvar, data=this%grainc_xfer_patch)
          if (flag=='read' .and. .not. readvar) then
             call set_missing_from_template( &
                  my_var = this%grainc_xfer_patch, &
                  template_var = c12_cnveg_carbonstate_inst%grainc_xfer_patch, &
                  multiplier = c3_r2)
          end if

          call restartvar(ncid=ncid, flag=flag, varname='cropseedc_14_deficit', xtype=ncd_double,  &
               dim1name='pft', long_name='pool for seeding new crop growth', units='gC14/m2', &
               interpinic_flag='interp', readvar=readvar, data=this%cropseedc_deficit_patch)
          if (flag=='read' .and. .not. readvar) then
             if ( masterproc ) write(iulog,*) 'initializing this%cropseedc_deficit_patch with atmospheric c14 value'
             call set_missing_from_template( &
                  my_var = this%cropseedc_deficit_patch, &
                  template_var = c12_cnveg_carbonstate_inst%cropseedc_deficit_patch, &
                  multiplier = c14ratio)
          end if
       end if
    end if

    !--------------------------------
    ! gridcell carbon state variables
    !--------------------------------

    if (carbon_type == 'c12') then
       ! BACKWARDS_COMPATIBILITY(wjs, 2017-01-12) Naming this with a _g suffix in order
       ! to distinguish it from the old column-level seedc restart variable
       call restartvar(ncid=ncid, flag=flag, varname='seedc_g', xtype=ncd_double,  &
            dim1name='gridcell', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%seedc_grc) 
    end if

    !--------------------------------
    ! C13 gridcell carbon state variables
    !--------------------------------

    if (carbon_type == 'c13') then
       call restartvar(ncid=ncid, flag=flag, varname='seedc_13_g', xtype=ncd_double,  &
            dim1name='gridcell', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%seedc_grc) 
       if (flag=='read' .and. .not. readvar) then
          call set_missing_from_template( &
               my_var = this%seedc_grc, &
               template_var = c12_cnveg_carbonstate_inst%seedc_grc, &
               multiplier = c3_r2)
       end if
    end if

    !--------------------------------
    ! C14 column carbon state variables
    !--------------------------------

    if ( carbon_type == 'c14' ) then
       call restartvar(ncid=ncid, flag=flag, varname='seedc_14_g', xtype=ncd_double,  &
            dim1name='gridcell', &
            long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%seedc_grc) 
       if (flag=='read' .and. .not. readvar) then
          if ( masterproc ) write(iulog,*) 'initializing this%seedc_grc with atmospheric c14 value'
          call set_missing_from_template( &
               my_var = this%seedc_grc, &
               template_var = c12_cnveg_carbonstate_inst%seedc_grc, &
               multiplier = c14ratio)
       end if
    end if

  end subroutine Restart

  !-----------------------------------------------------------------------
  subroutine SetValues ( this, &
       num_patch, filter_patch, value_patch, &
       num_column, filter_column, value_column)
    !
    ! !DESCRIPTION:
    ! Set carbon state variables
    !
    ! !ARGUMENTS:
    class (cnveg_carbonstate_type) :: this
    integer , intent(in) :: num_patch
    integer , intent(in) :: filter_patch(:)
    real(r8), intent(in) :: value_patch
    integer , intent(in) :: num_column
    integer , intent(in) :: filter_column(:)
    real(r8), intent(in) :: value_column
    !
    ! !LOCAL VARIABLES:
    integer :: fi,i,j,k,l     ! loop index
    !------------------------------------------------------------------------

    do fi = 1,num_patch
       i  = filter_patch(fi)
       this%leafc_patch(i)              = value_patch
       this%leafc_storage_patch(i)      = value_patch
       this%leafc_xfer_patch(i)         = value_patch
       this%leafc_storage_xfer_acc_patch(i) = value_patch
       this%storage_cdemand_patch(i)        = value_patch        
       this%frootc_patch(i)             = value_patch
       this%frootc_storage_patch(i)     = value_patch
       this%frootc_xfer_patch(i)        = value_patch
       this%livestemc_patch(i)          = value_patch
       this%livestemc_storage_patch(i)  = value_patch
       this%livestemc_xfer_patch(i)     = value_patch
       this%deadstemc_patch(i)          = value_patch
       this%deadstemc_storage_patch(i)  = value_patch
       this%deadstemc_xfer_patch(i)     = value_patch
       this%livecrootc_patch(i)         = value_patch
       this%livecrootc_storage_patch(i) = value_patch
       this%livecrootc_xfer_patch(i)    = value_patch
       this%deadcrootc_patch(i)         = value_patch
       this%deadcrootc_storage_patch(i) = value_patch
       this%deadcrootc_xfer_patch(i)    = value_patch
       this%gresp_storage_patch(i)      = value_patch
       this%gresp_xfer_patch(i)         = value_patch
       this%cpool_patch(i)              = value_patch
       this%xsmrpool_patch(i)           = value_patch
       this%ctrunc_patch(i)             = value_patch
       this%dispvegc_patch(i)           = value_patch
       this%storvegc_patch(i)           = value_patch
       this%woodc_patch(i)              = value_patch
       this%totvegc_patch(i)            = value_patch
       this%totc_patch(i)               = value_patch
       if ( use_crop ) then
          this%grainc_patch(i)          = value_patch
          this%grainc_storage_patch(i)  = value_patch
          this%grainc_xfer_patch(i)     = value_patch
          this%cropseedc_deficit_patch(i)  = value_patch
       end if
    end do

    do fi = 1,num_column
       i  = filter_column(fi)
       this%rootc_col(i)                = value_column
       this%leafc_col(i)                = value_column
       this%deadstemc_col(i)            = value_column
       this%fuelc_col(i)                = value_column
       this%fuelc_crop_col(i)           = value_column
       this%totvegc_col(i)              = value_column
       this%totc_p2c_col(i)             = value_column
       this%totc_col(i)                 = value_column
       this%totecosysc_col(i)           = value_column
    end do

  end subroutine SetValues

  !-----------------------------------------------------------------------
  subroutine ZeroDwt( this, bounds )
    !
    ! !DESCRIPTION
    ! Initialize variables needed for dynamic land use.
    !
    ! !ARGUMENTS:
    class(cnveg_carbonstate_type) :: this
    type(bounds_type), intent(in)  :: bounds 
    !
    ! !LOCAL VARIABLES:
    integer  :: p          ! indices
    !-----------------------------------------------------------------------

    do p = bounds%begp,bounds%endp
       this%dispvegc_patch(p)   = 0._r8
       this%storvegc_patch(p)   = 0._r8
       this%totc_patch(p)       = 0._r8
    end do

  end subroutine ZeroDwt

  !-----------------------------------------------------------------------
  subroutine Summary_carbonstate(this, bounds, num_allc, filter_allc, &
       num_soilc, filter_soilc, num_soilp, filter_soilp, &
       soilbiogeochem_cwdc_col, soilbiogeochem_totlitc_col, soilbiogeochem_totsomc_col, &
       soilbiogeochem_ctrunc_col)
    !
    ! !USES:
    use subgridAveMod, only : p2c
    use clm_time_manager , only : get_nstep

    !
    ! !DESCRIPTION:
    ! Perform patch and column-level carbon summary calculations
    !
    ! !ARGUMENTS:
    class(cnveg_carbonstate_type)  :: this
    type(bounds_type) , intent(in) :: bounds          
    integer           , intent(in) :: num_allc        ! number of columns in allc filter
    integer           , intent(in) :: filter_allc(:)  ! filter for all active columns
    integer           , intent(in) :: num_soilc       ! number of soil columns in filter
    integer           , intent(in) :: filter_soilc(:) ! filter for soil columns
    integer           , intent(in) :: num_soilp       ! number of soil patches in filter
    integer           , intent(in) :: filter_soilp(:) ! filter for soil patches
    real(r8)          , intent(in) :: soilbiogeochem_cwdc_col(bounds%begc:)   
    real(r8)          , intent(in) :: soilbiogeochem_totlitc_col(bounds%begc:)
    real(r8)          , intent(in) :: soilbiogeochem_totsomc_col(bounds%begc:)
    real(r8)          , intent(in) :: soilbiogeochem_ctrunc_col(bounds%begc:)
    !
    ! !LOCAL VARIABLES:
    integer  :: c,p,j,k,l       ! indices
    integer  :: fp,fc           ! lake filter indices
    !-----------------------------------------------------------------------

    SHR_ASSERT_ALL((ubound(soilbiogeochem_cwdc_col)    == (/bounds%endc/)), errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(soilbiogeochem_totlitc_col) == (/bounds%endc/)), errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(soilbiogeochem_totsomc_col) == (/bounds%endc/)), errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(soilbiogeochem_ctrunc_col)  == (/bounds%endc/)), errMsg(sourcefile, __LINE__))

    ! calculate patch -level summary of carbon state

    do fp = 1,num_soilp
       p = filter_soilp(fp)

       ! displayed vegetation carbon, excluding storage and cpool (DISPVEGC)
       this%dispvegc_patch(p) =        &
            this%leafc_patch(p)      + &
            this%frootc_patch(p)     + &
            this%livestemc_patch(p)  + &
            this%deadstemc_patch(p)  + &
            this%livecrootc_patch(p) + &
            this%deadcrootc_patch(p)

       ! stored vegetation carbon, excluding cpool (STORVEGC)
       this%storvegc_patch(p) =                &
            this%cpool_patch(p)              + &
            this%leafc_storage_patch(p)      + &
            this%frootc_storage_patch(p)     + &
            this%livestemc_storage_patch(p)  + &
            this%deadstemc_storage_patch(p)  + &
            this%livecrootc_storage_patch(p) + &
            this%deadcrootc_storage_patch(p) + &
            this%leafc_xfer_patch(p)         + &
            this%frootc_xfer_patch(p)        + &
            this%livestemc_xfer_patch(p)     + &
            this%deadstemc_xfer_patch(p)     + &
            this%livecrootc_xfer_patch(p)    + &
            this%deadcrootc_xfer_patch(p)    + &
            this%gresp_storage_patch(p)      + &
            this%gresp_xfer_patch(p)

       if ( use_crop .and. patch%itype(p) >= npcropmin )then
          this%storvegc_patch(p) =            &
               this%storvegc_patch(p)       + &
               this%grainc_storage_patch(p) + &
               this%grainc_xfer_patch(p)

          this%dispvegc_patch(p) =            &
               this%dispvegc_patch(p)       + &
               this%grainc_patch(p)
       end if

       ! total vegetation carbon, excluding cpool (TOTVEGC)
       this%totvegc_patch(p) = &
            this%dispvegc_patch(p) + &
            this%storvegc_patch(p)

       ! total patch-level carbon, including xsmrpool, ctrunc
       this%totc_patch(p) = &
            this%totvegc_patch(p) + &
            this%xsmrpool_patch(p) + &
            this%ctrunc_patch(p)

       if (use_crop) then 
          this%totc_patch(p) = this%totc_patch(p) + this%cropseedc_deficit_patch(p)
       end if

       ! (WOODC) - wood C
       this%woodc_patch(p) = &
            this%deadstemc_patch(p)    + &
            this%livestemc_patch(p)    + &
            this%deadcrootc_patch(p)   + &
            this%livecrootc_patch(p)

    end do

    ! --------------------------------------------
    ! column level summary
    ! --------------------------------------------

    call p2c(bounds, num_soilc, filter_soilc, &
         this%totvegc_patch(bounds%begp:bounds%endp), &
         this%totvegc_col(bounds%begc:bounds%endc))

    call p2c(bounds, num_soilc, filter_soilc, &
         this%totc_patch(bounds%begp:bounds%endp), &
         this%totc_p2c_col(bounds%begc:bounds%endc))

    do fc = 1,num_allc
       c = filter_allc(fc)

       ! total ecosystem carbon, including veg but excluding cpool (TOTECOSYSC)
       this%totecosysc_col(c) =    &
            soilbiogeochem_cwdc_col(c)    + &
            soilbiogeochem_totlitc_col(c) + &
            soilbiogeochem_totsomc_col(c) + &
            this%totvegc_col(c)

       ! total column carbon, including veg and cpool (TOTCOLC)
       this%totc_col(c) =  this%totc_p2c_col(c) + &
            soilbiogeochem_cwdc_col(c)      + &
            soilbiogeochem_totlitc_col(c)   + &
            soilbiogeochem_totsomc_col(c)   + &
            soilbiogeochem_ctrunc_col(c)

    end do

  end subroutine Summary_carbonstate

  !-----------------------------------------------------------------------
  subroutine DynamicPatchAdjustments(this, bounds, &
       num_soilp_with_inactive, filter_soilp_with_inactive, &
       patch_state_updater, &
       leafc_seed, deadstemc_seed, &
       conv_cflux, wood_product_cflux, crop_product_cflux, &
       dwt_frootc_to_litter, &
       dwt_livecrootc_to_litter, &
       dwt_deadcrootc_to_litter, &
       dwt_leafc_seed, &
       dwt_deadstemc_seed)
    !
    ! !DESCRIPTION:
    ! Adjust state variables and compute associated fluxes when patch areas change due to
    ! dynamic landuse
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(cnveg_carbonstate_type)   , intent(inout) :: this
    type(bounds_type)               , intent(in)    :: bounds
    integer                         , intent(in)    :: num_soilp_with_inactive ! number of points in filter
    integer                         , intent(in)    :: filter_soilp_with_inactive(:) ! soil patch filter that includes inactive points
    type(patch_state_updater_type)  , intent(in)    :: patch_state_updater
    real(r8)                        , intent(in)    :: leafc_seed  ! seed amount for leaf C
    real(r8)                        , intent(in)    :: deadstemc_seed ! seed amount for deadstem C
    real(r8)                        , intent(inout) :: conv_cflux( bounds%begp: )  ! patch-level conversion C flux to atm (expressed per unit GRIDCELL area)
    real(r8)                        , intent(inout) :: wood_product_cflux( bounds%begp: ) ! patch-level product C flux (expressed per unit GRIDCELL area)
    real(r8)                        , intent(inout) :: crop_product_cflux( bounds%begp: ) ! patch-level crop product C flux (expressed per unit GRIDCELL area)
    real(r8)                        , intent(inout) :: dwt_frootc_to_litter( bounds%begp: ) ! patch-level fine root C to litter (expressed per unit COLUMN area)
    real(r8)                        , intent(inout) :: dwt_livecrootc_to_litter( bounds%begp: ) ! patch-level live coarse root C to litter (expressed per unit COLUMN area)
    real(r8)                        , intent(inout) :: dwt_deadcrootc_to_litter( bounds%begp: ) ! patch-level live coarse root C to litter (expressed per unit COLUMN area)
    real(r8)                        , intent(inout) :: dwt_leafc_seed( bounds%begp: ) ! patch-level mass gain due to seeding of new area: leaf C (expressed per unit GRIDCELL area)
    real(r8)                        , intent(inout) :: dwt_deadstemc_seed( bounds%begp: ) ! patch-level mass gain due to seeding of new area: deadstem C (expressed per unit GRIDCELL area)
    !
    ! !LOCAL VARIABLES:
    integer :: begp, endp

    logical  :: old_weight_was_zero(bounds%begp:bounds%endp)
    logical  :: patch_grew(bounds%begp:bounds%endp)

    ! The following are only set for growing patches:
    real(r8) :: seed_leafc_patch(bounds%begp:bounds%endp)
    real(r8) :: seed_leafc_storage_patch(bounds%begp:bounds%endp)
    real(r8) :: seed_leafc_xfer_patch(bounds%begp:bounds%endp)
    real(r8) :: seed_deadstemc_patch(bounds%begp:bounds%endp)

    character(len=*), parameter :: subname = 'DynamicPatchAdjustments'
    !-----------------------------------------------------------------------

    begp = bounds%begp
    endp = bounds%endp

    SHR_ASSERT_ALL((ubound(conv_cflux) == (/endp/)), errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(wood_product_cflux) == (/endp/)), errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(crop_product_cflux) == (/endp/)), errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(dwt_frootc_to_litter) == (/endp/)), errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(dwt_livecrootc_to_litter) == (/endp/)), errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(dwt_deadcrootc_to_litter) == (/endp/)), errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(dwt_leafc_seed) == (/endp/)), errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(dwt_deadstemc_seed) == (/endp/)), errMsg(sourcefile, __LINE__))

    old_weight_was_zero = patch_state_updater%old_weight_was_zero(bounds)
    patch_grew = patch_state_updater%patch_grew(bounds)

    call ComputeSeedAmounts(bounds, &
         num_soilp_with_inactive, filter_soilp_with_inactive, &
         species = this%species, &
         leafc_seed = leafc_seed, &
         deadstemc_seed = deadstemc_seed, &
         leaf_patch = this%leafc_patch(begp:endp), &
         leaf_storage_patch = this%leafc_storage_patch(begp:endp), &
         leaf_xfer_patch = this%leafc_xfer_patch(begp:endp), &

         ! Calculations only needed for patches that grew:
         compute_here_patch = patch_grew(begp:endp), &

         ! For patches that previously had zero area, ignore the current state for the
         ! sake of computing leaf proportions:
         ignore_current_state_patch = old_weight_was_zero(begp:endp), &

         seed_leaf_patch = seed_leafc_patch(begp:endp), &
         seed_leaf_storage_patch = seed_leafc_storage_patch(begp:endp), &
         seed_leaf_xfer_patch = seed_leafc_xfer_patch(begp:endp), &
         seed_deadstem_patch = seed_deadstemc_patch(begp:endp))

    call update_patch_state( &
         var = this%leafc_patch(begp:endp), &
         flux_out_grc_area = conv_cflux(begp:endp), &
         seed = seed_leafc_patch(begp:endp), &
         seed_addition = dwt_leafc_seed(begp:endp))

    call update_patch_state( &
         var = this%leafc_storage_patch(begp:endp), &
         flux_out_grc_area = conv_cflux(begp:endp), &
         seed = seed_leafc_storage_patch(begp:endp), &
         seed_addition = dwt_leafc_seed(begp:endp))

    call update_patch_state( &
         var = this%leafc_xfer_patch(begp:endp), &
         flux_out_grc_area = conv_cflux(begp:endp), &
         seed = seed_leafc_xfer_patch(begp:endp), &
         seed_addition = dwt_leafc_seed(begp:endp))

    call update_patch_state( &
         var = this%frootc_patch(begp:endp), &
         flux_out_col_area = dwt_frootc_to_litter(begp:endp))

    call update_patch_state( &
         var = this%frootc_storage_patch(begp:endp), &
         flux_out_grc_area = conv_cflux(begp:endp))

    call update_patch_state( &
         var = this%frootc_xfer_patch(begp:endp), &
         flux_out_grc_area = conv_cflux(begp:endp))

    call update_patch_state( &
         var = this%livestemc_patch(begp:endp), &
         flux_out_grc_area = conv_cflux(begp:endp))

    call update_patch_state( &
         var = this%livestemc_storage_patch(begp:endp), &
         flux_out_grc_area = conv_cflux(begp:endp))

    call update_patch_state( &
         var = this%livestemc_xfer_patch(begp:endp), &
         flux_out_grc_area = conv_cflux(begp:endp))

    call patch_state_updater%update_patch_state_partition_flux_by_type(bounds, &
         num_soilp_with_inactive, filter_soilp_with_inactive, &
         flux1_fraction_by_pft_type = pftcon%pconv, &
         var = this%deadstemc_patch(begp:endp), &
         flux1_out = conv_cflux(begp:endp), &
         flux2_out = wood_product_cflux(begp:endp), &
         seed = seed_deadstemc_patch(begp:endp), &
         seed_addition = dwt_deadstemc_seed(begp:endp))

    call update_patch_state( &
         var = this%deadstemc_storage_patch(begp:endp), &
         flux_out_grc_area = conv_cflux(begp:endp))

    call update_patch_state( &
         var = this%deadstemc_xfer_patch(begp:endp), &
         flux_out_grc_area = conv_cflux(begp:endp))

    call update_patch_state( &
         var = this%livecrootc_patch(begp:endp), &
         flux_out_col_area = dwt_livecrootc_to_litter(begp:endp))

    call update_patch_state( &
         var = this%livecrootc_storage_patch(begp:endp), &
         flux_out_grc_area = conv_cflux(begp:endp))

    call update_patch_state( &
         var = this%livecrootc_xfer_patch(begp:endp), &
         flux_out_grc_area = conv_cflux(begp:endp))

    call update_patch_state( &
         var = this%deadcrootc_patch(begp:endp), &
         flux_out_col_area = dwt_deadcrootc_to_litter(begp:endp))

    call update_patch_state( &
         var = this%deadcrootc_storage_patch(begp:endp), &
         flux_out_grc_area = conv_cflux(begp:endp))

    call update_patch_state( &
         var = this%deadcrootc_xfer_patch(begp:endp), &
         flux_out_grc_area = conv_cflux(begp:endp))

    call update_patch_state( &
         var = this%gresp_storage_patch(begp:endp), &
         flux_out_grc_area = conv_cflux(begp:endp))

    call update_patch_state( &
         var = this%gresp_xfer_patch(begp:endp), &
         flux_out_grc_area = conv_cflux(begp:endp))

    call update_patch_state( &
         var = this%cpool_patch(begp:endp), &
         flux_out_grc_area = conv_cflux(begp:endp))

    call update_patch_state( &
         var = this%xsmrpool_patch(begp:endp), &
         flux_out_grc_area = conv_cflux(begp:endp))

    call update_patch_state( &
         var = this%ctrunc_patch(begp:endp), &
         flux_out_grc_area = conv_cflux(begp:endp))

    if (use_crop) then
       call update_patch_state( &
            var = this%grainc_patch(begp:endp), &
            flux_out_grc_area = crop_product_cflux(begp:endp))

       call update_patch_state( &
            var = this%grainc_storage_patch(begp:endp), &
            flux_out_grc_area = conv_cflux(begp:endp))

       call update_patch_state( &
            var = this%grainc_xfer_patch(begp:endp), &
            flux_out_grc_area = conv_cflux(begp:endp))

       if (use_crop) then
          ! This is a negative pool. So any deficit that we haven't repaid gets sucked out
          ! of the atmosphere.
          call update_patch_state( &
               var = this%cropseedc_deficit_patch(begp:endp), &
               flux_out_grc_area = conv_cflux(begp:endp))
       end if
    end if

  contains
    subroutine update_patch_state(var, flux_out_col_area, flux_out_grc_area, &
         seed, seed_addition)
      ! Wraps call to update_patch_state, in order to remove duplication
      real(r8), intent(inout) :: var( bounds%begp: )
      real(r8), intent(inout), optional :: flux_out_col_area( bounds%begp: )
      real(r8), intent(inout), optional :: flux_out_grc_area( bounds%begp: )
      real(r8), intent(in), optional :: seed( bounds%begp: )
      real(r8), intent(inout), optional :: seed_addition( bounds%begp: )

      call patch_state_updater%update_patch_state(bounds, &
         num_soilp_with_inactive, filter_soilp_with_inactive, &
         var = var, &
         flux_out_col_area = flux_out_col_area, &
         flux_out_grc_area = flux_out_grc_area, &
         seed = seed, &
         seed_addition = seed_addition)
    end subroutine update_patch_state

  end subroutine DynamicPatchAdjustments

end module CNVegCarbonStateType
