module SectorWaterMod

    #include "shr_assert.h"
    use shr_kind_mod     , only : r8 => shr_kind_r8
    use decompMod        , only : bounds_type
    use shr_log_mod      , only : errMsg => shr_log_errMsg
    use abortutils       , only : endrun
    use clm_varctl       , only : iulog
    use clm_time_manager , only : get_curr_date
    use WaterType        , only : water_type
    use GridcellType     , only : grc           
    use clm_time_manager , only : get_step_size              
    use ncdio_pio
    
 
  implicit none
  private
 
  ! !PUBLIC TYPES:
 
  ! This type is public (and its components are public, too) to aid unit testing
  type, public :: sectorwater_params_type
       ! Threshold for how much of the current river water can be used to satisfy sectoral demand, 
       ! if limit_sectorwater is .true. (fraction of available river water). A threshold of 0
       ! means allow all current river water to be used; a threshold of 0.1 means allow 90% of the current
       ! river volume to be used; etc.
       ! This is done to protect against negative runoff, caused by extracting more water than currently available.
       real(r8) :: sectorwater_river_volume_threshold
 
       ! Whether sectorwater usage is limited based on river storage. This only applies if ROF is
       ! enabled (i.e., rof_prognostic is .true.) - otherwise we don't limit sectorwater usage,
       ! regardless of the value of this flag.
       logical :: limit_sectorwater_if_rof_enabled      
       
       ! The sectoral water usage is computed based on the provided input data
       ! path_sectorwater_input_data is the path to a .txt file containing the paths to each year .nc file with sectoral withdrawal and consumption data
       ! The right format of the .txt file is:
       ! (line 1) first year and last year provided as integers separated by coma (e.g. 1971,2010)
       ! (line 2) path_directory_with_the_input_sectorwater_data/name_file_for_year_1971.nc (here 1971 provided as example)
       ! (line n) path_directory_with_the_input_sectorwater_data/name_file_for_year_n.nc (n represent a year between 1971 and 2010)
       ! (line N) path_directory_with_the_input_sectorwater_data/name_file_for_year_2010.nc (N line corresponds to the last year input data path)
       ! N.B. Each new line represent the path to a monthly input dataset for the given year. It is the user responsability to make sure that the inputs path are given in the right consecutive order from first to last year.
       ! N.B. Make sure that there is no NaN values in your input files for the sector water usage variables. Instead of NaN values, use 0.
       character(len=256) :: path_sectorwater_input_data
 end type sectorwater_params_type
 
 
 type, public :: sectorwater_type
       ! private
       ! Private data members; set in initialization:
       
       type(sectorwater_params_type) :: params

       integer :: dtime                ! land model time step (sec)
 
       ! Private data members; time-varying:
       ! naming: dom = domestic, liv = livestock, elec = thermoelectric, mfc = manufacturing, min = mining
       ! naming: withd = withdrawal, cons = consumption, rf = return flow
 
       real(r8), pointer :: input_mon_dom_withd_grc   (:) ! input expected withdrawal for current month
       real(r8), pointer :: input_mon_dom_cons_grc    (:) ! input expected consumption for current month
       real(r8), pointer :: dom_withd_grc             (:) ! expected withdrawal flux for the day [mm/s]
       real(r8), pointer :: dom_cons_grc              (:) ! expected consumption flux for the day [mm/s]
       real(r8), pointer :: dom_withd_actual_grc      (:) ! actual withdrawal flux for the day [mm/s]
       real(r8), pointer :: dom_cons_actual_grc       (:) ! actual consumption flux for the day  [mm/s]
       real(r8), pointer :: dom_rf_actual_grc         (:) ! actual return flow flux for the day  [mm/s]
 
       real(r8), pointer :: input_mon_liv_withd_grc   (:) ! input expected withdrawal for current month
       real(r8), pointer :: input_mon_liv_cons_grc    (:) ! input expected consumption for current month
       real(r8), pointer :: liv_withd_grc             (:) ! expected withdrawal flux for the day [mm/s]
       real(r8), pointer :: liv_cons_grc              (:) ! expected consumption flux for the day [mm/s]
       real(r8), pointer :: liv_withd_actual_grc      (:) ! actual withdrawal flux for the day [mm/s]
       real(r8), pointer :: liv_cons_actual_grc       (:) ! actual consumption flux for the day  [mm/s]
       real(r8), pointer :: liv_rf_actual_grc         (:) ! actual return flow flux for the day  [mm/s]
 
       real(r8), pointer :: input_mon_elec_withd_grc   (:) ! input expected withdrawal for current month
       real(r8), pointer :: input_mon_elec_cons_grc    (:) ! input expected consumption for current month
       real(r8), pointer :: elec_withd_grc             (:) ! expected withdrawal flux for the day [mm/s]
       real(r8), pointer :: elec_cons_grc              (:) ! expected consumption flux for the day [mm/s]
       real(r8), pointer :: elec_withd_actual_grc      (:) ! actual withdrawal flux for the day [mm/s]
       real(r8), pointer :: elec_cons_actual_grc       (:) ! actual consumption flux for the day  [mm/s]
       real(r8), pointer :: elec_rf_actual_grc         (:) ! actual return flow flux for the day  [mm/s]
 
       real(r8), pointer :: input_mon_mfc_withd_grc   (:) ! input expected withdrawal for current month
       real(r8), pointer :: input_mon_mfc_cons_grc    (:) ! input expected consumption for current month
       real(r8), pointer :: mfc_withd_grc             (:) ! expected withdrawal flux for the day [mm/s]
       real(r8), pointer :: mfc_cons_grc              (:) ! expected consumption flux for the day [mm/s]
       real(r8), pointer :: mfc_withd_actual_grc      (:) ! actual withdrawal flux for the day [mm/s]
       real(r8), pointer :: mfc_cons_actual_grc       (:) ! actual consumption flux for the day  [mm/s]
       real(r8), pointer :: mfc_rf_actual_grc         (:) ! actual return flow flux for the day  [mm/s]
 
       real(r8), pointer :: input_mon_min_withd_grc   (:) ! input expected withdrawal for current month
       real(r8), pointer :: input_mon_min_cons_grc    (:) ! input expected consumption for current month
       real(r8), pointer :: min_withd_grc             (:) ! expected withdrawal flux for the day [mm/s]
       real(r8), pointer :: min_cons_grc              (:) ! expected consumption flux for the day [mm/s]
       real(r8), pointer :: min_withd_actual_grc      (:) ! actual withdrawal flux for the day [mm/s]
       real(r8), pointer :: min_cons_actual_grc       (:) ! actual consumption flux for the day  [mm/s]
       real(r8), pointer :: min_rf_actual_grc         (:) ! actual return flow flux for the day  [mm/s]

       real(r8), pointer :: sectorwater_total_actual_withd (:) ! total actual water withdrawal for all sectors (except irrigation) during current day [m3]

       
       
    contains
       ! Public routines
       procedure, public :: Init => SectorWaterInit
       procedure, public :: ReadSectorWaterData
       procedure, public :: CalcSectorWaterNeeded
       procedure, public :: Clean => SectorWaterClean ! deallocate memory

       ! Public routines to be added:
       ! procedure, public :: Restart
 
       ! Private routines
       procedure, private :: ReadNamelist
       procedure, private :: CheckNamelistValidity   ! Check for validity of input parameters
       procedure, private :: InitAllocate => SectorWaterInitAllocate
       procedure, private :: InitHistory => SectorWaterInitHistory
       procedure, private :: InitCold => SectorWaterInitCold
       procedure, private :: CalcSectorDemandVolrLimited   ! calculate actual sectoral abstractions limited by river volume
 end type sectorwater_type
 
 interface sectorwater_params_type
    module procedure sectorwater_params_constructor
 end interface sectorwater_params_type
 
 real(r8), parameter :: m3_over_km2_to_mm = 1.e-3_r8
 real(r8), parameter :: mm_to_m3_over_km2 = 1.0/m3_over_km2_to_mm

 character(len=*), parameter, private :: sourcefile = &
    __FILE__
 
 contains
 
    ! ========================================================================
    ! Constructors
    ! ========================================================================
 
    !-----------------------------------------------------------------------
     function sectorwater_params_constructor(sectorwater_river_volume_threshold, &
          limit_sectorwater_if_rof_enabled, path_sectorwater_input_data) &
          result(this)
       !
       ! !DESCRIPTION:
       ! Create an sectorwater_params instance
       !
       ! !USES:
       !
       ! !ARGUMENTS:
       type(sectorwater_params_type) :: this  ! function result
       real(r8), intent(in) :: sectorwater_river_volume_threshold
       logical , intent(in) :: limit_sectorwater_if_rof_enabled
       character(len=256), intent(in) :: path_sectorwater_input_data
       !
       ! !LOCAL VARIABLES:
       
       character(len=*), parameter :: subname = 'sectorwater_params_constructor'
       !-----------------------------------------------------------------------
       this%sectorwater_river_volume_threshold = sectorwater_river_volume_threshold
       this%limit_sectorwater_if_rof_enabled = limit_sectorwater_if_rof_enabled
       this%path_sectorwater_input_data = path_sectorwater_input_data

     end function sectorwater_params_constructor
 
   ! ========================================================================
   ! Infrastructure routines (initialization, restart, etc.)
   ! ========================================================================
 
   !------------------------------------------------------------------------
     subroutine SectorWaterInit(this, bounds, NLFilename)
          class(sectorwater_type) , intent(inout) :: this
          type(bounds_type)      , intent(in)     :: bounds
          character(len=*)       , intent(in)    :: NLFilename ! Namelist filename
          call this%ReadNamelist(NLFilename)
          call this%InitAllocate(bounds)
          call this%InitHistory(bounds)
          call this%InitCold(bounds)
     end subroutine SectorWaterInit
 
   !-----------------------------------------------------------------------
     subroutine ReadNamelist(this, NLFilename)
          !
          !!DESCRIPTION:
          ! Read the sectorwater namelist
          !
          ! !USES:
          use fileutils      , only : getavu, relavu, opnfil
          use shr_nl_mod     , only : shr_nl_find_group_name
          use spmdMod        , only : masterproc, mpicom
          use shr_mpi_mod    , only : shr_mpi_bcast
          use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
          !
          !!ARGUMENTS:
          class(sectorwater_type) , intent(inout) :: this
          character(len=*), intent(in)            :: NLFilename ! Namelist filename
          !
          !!LOCAL VARIABLES:
 
          ! temporary variables corresponding to the components of sectorwater_params_type
          real(r8) :: sectorwater_river_volume_threshold
          logical  :: limit_sectorwater_if_rof_enabled
          character(len=256) :: path_sectorwater_input_data
          integer  :: ierr                 ! error code
          integer  :: unitn                ! unit for namelist file

          character(len=*), parameter :: nmlname_sectorwater = 'sectorwater_inparm'
          character(len=*), parameter :: subname = 'ReadNamelist'
          !-----------------------------------------------------------------------
 
          namelist /sectorwater_inparm/ sectorwater_river_volume_threshold, limit_sectorwater_if_rof_enabled, path_sectorwater_input_data
    
          ! Initialize options to garbage defaults, forcing all to be specified explicitly in
          ! order to get reasonable results
          sectorwater_river_volume_threshold = nan
          limit_sectorwater_if_rof_enabled = .false.
          path_sectorwater_input_data = ' '
 
          if (masterproc) then
               unitn = getavu()
               write(iulog,*) 'Read in '//nmlname_sectorwater//'  namelist'
               call opnfil (NLFilename, unitn, 'F')
               call shr_nl_find_group_name(unitn, nmlname_sectorwater, status=ierr)
               
               if (ierr == 0) then
                    read(unitn, nml=sectorwater_inparm, iostat=ierr)
                    
                    if (ierr /= 0) then
                         call endrun(msg="ERROR reading "//nmlname_sectorwater//"namelist"//errmsg(sourcefile, __LINE__))
                    end if
               else
                    call endrun(msg="ERROR could NOT find "//nmlname_sectorwater//"namelist"//errmsg(sourcefile, __LINE__))
               end if

               call relavu( unitn )

          end if

          call shr_mpi_bcast(sectorwater_river_volume_threshold, mpicom)
          call shr_mpi_bcast(limit_sectorwater_if_rof_enabled, mpicom)
          call shr_mpi_bcast(path_sectorwater_input_data, mpicom)

          this%params = sectorwater_params_type(sectorwater_river_volume_threshold = sectorwater_river_volume_threshold, &
          limit_sectorwater_if_rof_enabled = limit_sectorwater_if_rof_enabled, &
          path_sectorwater_input_data = path_sectorwater_input_data)
          
          if (masterproc) then
               write(iulog,*) ' '
               write(iulog,*) nmlname_sectorwater//' settings:'
               ! Write settings one-by-one rather than with a nml write because
               ! sectorwater_river_volume_threshold may be NaN
               write(iulog,*) 'limit_sectorwater_if_rof_enabled = ', limit_sectorwater_if_rof_enabled
               if (limit_sectorwater_if_rof_enabled) then
                    write(iulog,*) 'sectorwater_river_volume_threshold = ', sectorwater_river_volume_threshold
               end if
               write(iulog,*) 'path_sectorwater_input_data = ', path_sectorwater_input_data
               write(iulog,*) ' '
          end if
                  
     end subroutine ReadNamelist
 
     !-----------------------------------------------------------------------
     subroutine CheckNamelistValidity(this)
          !
          ! !DESCRIPTION:
          ! Check for validity of input parameters.
          !
          ! Assumes that the inputs have already been packed into 'this%params'.
          !
          ! Only needs to be called by the master task, since parameters are the same for all
          ! tasks.
          !
          ! !ARGUMENTS:
          class(sectorwater_type), intent(in) :: this
          !
          ! !LOCAL VARIABLES:
          
          character(len=*), parameter :: subname = 'CheckNamelistValidity'
          !-----------------------------------------------------------------------
          
          associate( &
               sectorwater_river_volume_threshold => this%params%sectorwater_river_volume_threshold, &
               limit_sectorwater_if_rof_enabled => this%params%limit_sectorwater_if_rof_enabled)
          
          if (limit_sectorwater_if_rof_enabled) then
               if (sectorwater_river_volume_threshold < 0._r8 .or. sectorwater_river_volume_threshold > 1._r8) then
                    write(iulog,*) ' ERROR: sectorwater_river_volume_threshold must be between 0 and 1'
                    write(iulog,*) ' sectorwater_river_volume_threshold = ', sectorwater_river_volume_threshold
                    call endrun(msg=' ERROR: sectorwater_river_volume_threshold must be between 0 and 1 ' // &
                         errMsg(sourcefile, __LINE__))
               end if
          end if
          end associate
 
     end subroutine CheckNamelistValidity
  
     !-----------------------------------------------------------------------
     subroutine SectorWaterInitAllocate(this, bounds)
          !
          ! !DESCRIPTION:
          ! Initialize sector water data structure
          !
          ! !USES:
          use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
          !
          ! !ARGUMENTS:
          class(sectorwater_type) , intent(inout) :: this
          type(bounds_type)      , intent(in)    :: bounds
          !
          ! !LOCAL VARIABLES:
          integer :: begg, endg
          
          character(len=*), parameter :: subname = 'InitAllocate'
          !-----------------------------------------------------------------------
 
          begg = bounds%begg; endg= bounds%endg   
          
          allocate(this%input_mon_dom_withd_grc (begg:endg))            ; this%input_mon_dom_withd_grc        (:)   = 0
          allocate(this%input_mon_dom_cons_grc (begg:endg))             ; this%input_mon_dom_cons_grc         (:)   = 0
          allocate(this%dom_withd_grc (begg:endg))                      ; this%dom_withd_grc                  (:)   = 0
          allocate(this%dom_cons_grc (begg:endg))                       ; this%dom_cons_grc                   (:)   = 0
          allocate(this%dom_withd_actual_grc (begg:endg))               ; this%dom_withd_actual_grc           (:)   = 0
          allocate(this%dom_cons_actual_grc (begg:endg))                ; this%dom_cons_actual_grc            (:)   = 0
          allocate(this%dom_rf_actual_grc (begg:endg))                  ; this%dom_rf_actual_grc              (:)   = 0
 
          allocate(this%input_mon_liv_withd_grc (begg:endg))            ; this%input_mon_liv_withd_grc        (:)   = 0
          allocate(this%input_mon_liv_cons_grc (begg:endg))             ; this%input_mon_liv_cons_grc         (:)   = 0
          allocate(this%liv_withd_grc (begg:endg))                      ; this%liv_withd_grc                  (:)   = 0
          allocate(this%liv_cons_grc (begg:endg))                       ; this%liv_cons_grc                   (:)   = 0
          allocate(this%liv_withd_actual_grc (begg:endg))               ; this%liv_withd_actual_grc           (:)   = 0
          allocate(this%liv_cons_actual_grc (begg:endg))                ; this%liv_cons_actual_grc            (:)   = 0
          allocate(this%liv_rf_actual_grc (begg:endg))                  ; this%liv_rf_actual_grc              (:)   = 0
 
          allocate(this%input_mon_elec_withd_grc (begg:endg))           ; this%input_mon_elec_withd_grc       (:)   = 0
          allocate(this%input_mon_elec_cons_grc (begg:endg))            ; this%input_mon_elec_cons_grc        (:)   = 0
          allocate(this%elec_withd_grc (begg:endg))                     ; this%elec_withd_grc                 (:)   = 0
          allocate(this%elec_cons_grc (begg:endg))                      ; this%elec_cons_grc                  (:)   = 0
          allocate(this%elec_withd_actual_grc (begg:endg))              ; this%elec_withd_actual_grc          (:)   = 0
          allocate(this%elec_cons_actual_grc (begg:endg))               ; this%elec_cons_actual_grc           (:)   = 0
          allocate(this%elec_rf_actual_grc (begg:endg))                 ; this%elec_rf_actual_grc             (:)   = 0
 
          allocate(this%input_mon_mfc_withd_grc (begg:endg))            ; this%input_mon_mfc_withd_grc        (:)   = 0
          allocate(this%input_mon_mfc_cons_grc (begg:endg))             ; this%input_mon_mfc_cons_grc         (:)   = 0
          allocate(this%mfc_withd_grc (begg:endg))                      ; this%mfc_withd_grc                  (:)   = 0
          allocate(this%mfc_cons_grc (begg:endg))                       ; this%mfc_cons_grc                   (:)   = 0
          allocate(this%mfc_withd_actual_grc (begg:endg))               ; this%mfc_withd_actual_grc           (:)   = 0
          allocate(this%mfc_cons_actual_grc (begg:endg))                ; this%mfc_cons_actual_grc            (:)   = 0
          allocate(this%mfc_rf_actual_grc (begg:endg))                  ; this%mfc_rf_actual_grc              (:)   = 0
 
          allocate(this%input_mon_min_withd_grc (begg:endg))            ; this%input_mon_min_withd_grc        (:)   = 0
          allocate(this%input_mon_min_cons_grc (begg:endg))             ; this%input_mon_min_cons_grc         (:)   = 0
          allocate(this%min_withd_grc (begg:endg))                      ; this%min_withd_grc                  (:)   = 0
          allocate(this%min_cons_grc (begg:endg))                       ; this%min_cons_grc                   (:)   = 0
          allocate(this%min_withd_actual_grc (begg:endg))               ; this%min_withd_actual_grc           (:)   = 0
          allocate(this%min_cons_actual_grc (begg:endg))                ; this%min_cons_actual_grc            (:)   = 0
          allocate(this%min_rf_actual_grc (begg:endg))                  ; this%min_rf_actual_grc              (:)   = 0

          allocate(this%sectorwater_total_actual_withd (begg:endg))     ; this%sectorwater_total_actual_withd (:)   = 0

     end subroutine SectorWaterInitAllocate
 
 
     !-----------------------------------------------------------------------
     subroutine SectorWaterInitHistory(this, bounds)
          !
          ! !DESCRIPTION:
          ! Initialize sectoral water use history fields
          !
          ! !USES:
          use histFileMod  , only : hist_addfld1d
          !
          ! !ARGUMENTS:
          class(sectorWater_type) , intent(inout) :: this
          type(bounds_type)      , intent(in)    :: bounds
          !
          ! !LOCAL VARIABLES:
          integer :: begg, endg
    
          character(len=*), parameter :: subname = 'InitHistory'
      !-----------------------------------------------------------------------
 
          begg = bounds%begg; endg= bounds%endg

          ! Add output variables
          ! Domestic sector:
          this%input_mon_dom_withd_grc(begg:endg) = 0
          call hist_addfld1d (fname='INPUT_MON_DOM_WITHD', units='mm', &
               avgflag='A', long_name='input monthly domestic withdrawal', &
               ptr_gcell=this%input_mon_dom_withd_grc, default='inactive')
          
          this%input_mon_dom_cons_grc(begg:endg) = 0
          call hist_addfld1d (fname='INPUT_MON_DOM_CONS', units='mm', &
               avgflag='A', long_name='input monthly domestic consumption', &
               ptr_gcell=this%input_mon_dom_cons_grc, default='inactive')

          this%dom_cons_grc(begg:endg) = 0
          call hist_addfld1d (fname='DOM_EXPECTED_CONS', units='mm/s', &
               avgflag='A', long_name='domestic expected consumption flux', &
               ptr_gcell=this%dom_cons_grc, default='inactive')

          this%dom_cons_actual_grc(begg:endg) = 0
          call hist_addfld1d (fname='DOM_ACTUAL_CONS', units='mm/s', &
               avgflag='A', long_name='domestic actual consumption flux', &
               ptr_gcell=this%dom_cons_actual_grc, default='inactive')
 
          this%dom_withd_actual_grc(begg:endg) = 0
          call hist_addfld1d (fname='DOM_ACTUAL_WITHD', units='mm/s', &
               avgflag='A', long_name='domestic actual withdrawal flux', &
               ptr_gcell=this%dom_withd_actual_grc, default='inactive')

          this%dom_withd_grc(begg:endg) = 0
          call hist_addfld1d (fname='DOM_EXPECTED_WITHD', units='mm/s', &
               avgflag='A', long_name='domestic expected withdrawal flux', &
               ptr_gcell=this%dom_withd_grc, default='inactive')

          this%dom_rf_actual_grc(begg:endg) = 0
          call hist_addfld1d (fname='DOM_ACTUAL_RF', units='mm/s', &
               avgflag='A', long_name='domestic actual return flow flux', &
               ptr_gcell=this%dom_rf_actual_grc, default='inactive')
 
          ! Livestock
          this%input_mon_liv_withd_grc(begg:endg) = 0
          call hist_addfld1d (fname='INPUT_MON_LIV_WITHD', units='mm', &
               avgflag='A', long_name='input monthly livestock withdrawal', &
               ptr_gcell=this%input_mon_liv_withd_grc, default='inactive')
          
          this%input_mon_liv_cons_grc(begg:endg) = 0
          call hist_addfld1d (fname='INPUT_MON_LIV_CONS', units='mm', &
               avgflag='A', long_name='input monthly livestock consumption', &
               ptr_gcell=this%input_mon_liv_cons_grc, default='inactive')

          this%liv_cons_grc(begg:endg) = 0
          call hist_addfld1d (fname='LIV_EXPECTED_CONS', units='mm/s', &
               avgflag='A', long_name='livestock expected consumption flux', &
               ptr_gcell=this%liv_cons_grc, default='inactive')

          this%liv_cons_actual_grc(begg:endg) = 0
          call hist_addfld1d (fname='LIV_ACTUAL_CONS', units='mm/s', &
               avgflag='A', long_name='livestock actual consumption flux', &
               ptr_gcell=this%liv_cons_actual_grc, default='inactive')
          
          this%liv_withd_actual_grc(begg:endg) = 0
          call hist_addfld1d (fname='LIV_ACTUAL_WITHD', units='mm/s', &
               avgflag='A', long_name='livestock actual withdrawal flux', &
               ptr_gcell=this%liv_withd_actual_grc, default='inactive')

          this%liv_withd_grc(begg:endg) = 0
          call hist_addfld1d (fname='LIV_EXPECTED_WITHD', units='mm/s', &
               avgflag='A', long_name='livestock expected withdrawal flux', &
               ptr_gcell=this%liv_withd_grc, default='inactive')

          this%liv_rf_actual_grc(begg:endg) = 0
          call hist_addfld1d (fname='LIV_ACTUAL_RF', units='mm/s', &
               avgflag='A', long_name='livestock actual return flow flux', &
               ptr_gcell=this%liv_rf_actual_grc, default='inactive')

          ! Thermoelectric:
          this%input_mon_elec_withd_grc(begg:endg) = 0
          call hist_addfld1d (fname='INPUT_MON_ELEC_WITHD', units='mm', &
               avgflag='A', long_name='input monthly thermoelectric withdrawal', &
               ptr_gcell=this%input_mon_elec_withd_grc, default='inactive')
          
          this%input_mon_elec_cons_grc(begg:endg) = 0
          call hist_addfld1d (fname='INPUT_MON_ELEC_CONS', units='mm', &
               avgflag='A', long_name='input monthly thermoelectric consumption', &
               ptr_gcell=this%input_mon_elec_cons_grc, default='inactive')

          this%elec_cons_grc(begg:endg) = 0
          call hist_addfld1d (fname='ELEC_EXPECTED_CONS', units='mm/s', &
               avgflag='A', long_name='thermoelectric expected consumption flux', &
               ptr_gcell=this%elec_cons_grc, default='inactive')

          this%elec_cons_actual_grc(begg:endg) = 0
          call hist_addfld1d (fname='ELEC_ACTUAL_CONS', units='mm/s', &
               avgflag='A', long_name='thermoelectric actual consumption flux', &
               ptr_gcell=this%elec_cons_actual_grc, default='inactive')
          
          this%elec_withd_actual_grc(begg:endg) = 0
          call hist_addfld1d (fname='ELEC_ACTUAL_WITHD', units='mm/s', &
               avgflag='A', long_name='thermoelectric actual withdrawal flux', &
               ptr_gcell=this%elec_withd_actual_grc, default='inactive')

          this%elec_withd_grc(begg:endg) = 0
          call hist_addfld1d (fname='ELEC_EXPECTED_WITHD', units='mm/s', &
               avgflag='A', long_name='thermoelectric expected withdrawal flux', &
               ptr_gcell=this%elec_withd_grc, default='inactive')

          this%elec_rf_actual_grc(begg:endg) = 0
          call hist_addfld1d (fname='ELEC_ACTUAL_RF', units='mm/s', &
               avgflag='A', long_name='thermoelectric actual return flow flux', &
               ptr_gcell=this%elec_rf_actual_grc, default='inactive')

          ! Manufacturing:
          this%input_mon_mfc_withd_grc(begg:endg) = 0
          call hist_addfld1d (fname='INPUT_MON_MFC_WITHD', units='mm', &
               avgflag='A', long_name='input monthly manufacturing withdrawal', &
               ptr_gcell=this%input_mon_mfc_withd_grc, default='inactive')
          
          this%input_mon_mfc_cons_grc(begg:endg) = 0
          call hist_addfld1d (fname='INPUT_MON_MFC_CONS', units='mm', &
               avgflag='A', long_name='input monthly manufacturing consumption', &
               ptr_gcell=this%input_mon_mfc_cons_grc, default='inactive')

          this%mfc_cons_grc(begg:endg) = 0
          call hist_addfld1d (fname='MFC_EXPECTED_CONS', units='mm/s', &
               avgflag='A', long_name='manufacturing expected consumption flux', &
               ptr_gcell=this%mfc_cons_grc, default='inactive')

          this%mfc_cons_actual_grc(begg:endg) = 0
          call hist_addfld1d (fname='MFC_ACTUAL_CONS', units='mm/s', &
               avgflag='A', long_name='manufacturing actual consumption flux', &
               ptr_gcell=this%mfc_cons_actual_grc, default='inactive')
          
          this%mfc_withd_actual_grc(begg:endg) = 0
          call hist_addfld1d (fname='MFC_ACTUAL_WITHD', units='mm/s', &
               avgflag='A', long_name='manufacturing actual withdrawal flux', &
               ptr_gcell=this%mfc_withd_actual_grc, default='inactive')

          this%mfc_withd_grc(begg:endg) = 0
          call hist_addfld1d (fname='MFC_EXPECTED_WITHD', units='mm/s', &
               avgflag='A', long_name='manufacturing expected withdrawal flux', &
               ptr_gcell=this%mfc_withd_grc, default='inactive')

          this%mfc_rf_actual_grc(begg:endg) = 0
          call hist_addfld1d (fname='MFC_ACTUAL_RF', units='mm/s', &
               avgflag='A', long_name='manufacturing actual return flow flux', &
               ptr_gcell=this%mfc_rf_actual_grc, default='inactive')

          ! Mining:
          this%input_mon_min_withd_grc(begg:endg) = 0
          call hist_addfld1d (fname='INPUT_MON_MIN_WITHD', units='mm', &
               avgflag='A', long_name='input monthly mining withdrawal', &
               ptr_gcell=this%input_mon_min_withd_grc, default='inactive')
          
          this%input_mon_min_cons_grc(begg:endg) = 0
          call hist_addfld1d (fname='INPUT_MON_MIN_CONS', units='mm', &
               avgflag='A', long_name='input monthly mining consumption', &
               ptr_gcell=this%input_mon_min_cons_grc, default='inactive')

          this%min_cons_grc(begg:endg) = 0
          call hist_addfld1d (fname='MIN_EXPECTED_CONS', units='mm/s', &
               avgflag='A', long_name='mining expected consumption flux', &
               ptr_gcell=this%min_cons_grc, default='inactive')

          this%min_cons_actual_grc(begg:endg) = 0
          call hist_addfld1d (fname='MIN_ACTUAL_CONS', units='mm/s', &
               avgflag='A', long_name='mining actual consumption flux', &
               ptr_gcell=this%min_cons_actual_grc, default='inactive')
          
          this%min_withd_actual_grc(begg:endg) = 0
          call hist_addfld1d (fname='MIN_ACTUAL_WITHD', units='mm/s', &
               avgflag='A', long_name='mining actual withdrawal flux', &
               ptr_gcell=this%min_withd_actual_grc, default='inactive')

          this%min_withd_grc(begg:endg) = 0
          call hist_addfld1d (fname='MIN_EXPECTED_WITHD', units='mm/s', &
               avgflag='A', long_name='mining expected withdrawal flux', &
               ptr_gcell=this%min_withd_grc, default='inactive')

          this%min_rf_actual_grc(begg:endg) = 0
          call hist_addfld1d (fname='MIN_ACTUAL_RF', units='mm/s', &
               avgflag='A', long_name='mining actual return flow flux', &
               ptr_gcell=this%min_rf_actual_grc, default='inactive')
         
     end subroutine SectorWaterInitHistory
 
!-----------------------------------------------------------------------
     subroutine SectorWaterInitCold(this, bounds)
          !
          ! !DESCRIPTION:
          ! Do cold-start initialization for sector water abstractions data structure
          !
          ! !ARGUMENTS:
          class(sectorwater_type) , intent(inout) :: this
          type(bounds_type)      , intent(in)     :: bounds
          
          character(len=*), parameter :: subname = 'InitCold'
          !-----------------------------------------------------------------------  
          
          this%dtime = get_step_size()
          this%input_mon_dom_withd_grc(bounds%begg:bounds%endg)  = 0._r8
          this%input_mon_liv_withd_grc(bounds%begg:bounds%endg)  = 0._r8
          this%input_mon_elec_withd_grc(bounds%begg:bounds%endg) = 0._r8
          this%input_mon_mfc_withd_grc(bounds%begg:bounds%endg)  = 0._r8
          this%input_mon_min_withd_grc(bounds%begg:bounds%endg)  = 0._r8

          this%input_mon_dom_cons_grc(bounds%begg:bounds%endg)   = 0._r8
          this%input_mon_liv_cons_grc(bounds%begg:bounds%endg)   = 0._r8
          this%input_mon_elec_cons_grc(bounds%begg:bounds%endg)  = 0._r8
          this%input_mon_mfc_cons_grc(bounds%begg:bounds%endg)   = 0._r8
          this%input_mon_min_cons_grc(bounds%begg:bounds%endg)   = 0._r8

          this%dom_withd_grc(bounds%begg:bounds%endg)            = 0._r8
          this%liv_withd_grc(bounds%begg:bounds%endg)            = 0._r8
          this%elec_withd_grc(bounds%begg:bounds%endg)           = 0._r8
          this%mfc_withd_grc(bounds%begg:bounds%endg)            = 0._r8
          this%min_withd_grc(bounds%begg:bounds%endg)            = 0._r8

          this%dom_cons_grc(bounds%begg:bounds%endg)             = 0._r8
          this%liv_cons_grc(bounds%begg:bounds%endg)             = 0._r8
          this%elec_cons_grc(bounds%begg:bounds%endg)            = 0._r8
          this%mfc_cons_grc(bounds%begg:bounds%endg)             = 0._r8
          this%min_cons_grc(bounds%begg:bounds%endg)             = 0._r8

          this%dom_withd_actual_grc(bounds%begg:bounds%endg)     = 0._r8
          this%liv_withd_actual_grc(bounds%begg:bounds%endg)     = 0._r8
          this%elec_withd_actual_grc(bounds%begg:bounds%endg)    = 0._r8
          this%mfc_withd_actual_grc(bounds%begg:bounds%endg)     = 0._r8
          this%min_withd_actual_grc(bounds%begg:bounds%endg)     = 0._r8
          
          this%dom_cons_actual_grc(bounds%begg:bounds%endg)      = 0._r8
          this%liv_cons_actual_grc(bounds%begg:bounds%endg)      = 0._r8
          this%elec_cons_actual_grc(bounds%begg:bounds%endg)     = 0._r8
          this%mfc_cons_actual_grc(bounds%begg:bounds%endg)      = 0._r8
          this%min_cons_actual_grc(bounds%begg:bounds%endg)      = 0._r8
          
          this%dom_rf_actual_grc(bounds%begg:bounds%endg)        = 0._r8
          this%liv_rf_actual_grc(bounds%begg:bounds%endg)        = 0._r8
          this%elec_rf_actual_grc(bounds%begg:bounds%endg)       = 0._r8
          this%mfc_rf_actual_grc(bounds%begg:bounds%endg)        = 0._r8
          this%min_rf_actual_grc(bounds%begg:bounds%endg)        = 0._r8

          this%sectorwater_total_actual_withd(bounds%begg:bounds%endg)        = 0._r8

     
     end subroutine SectorWaterInitCold
 
     !-----------------------------------------------------------------------
     subroutine SectorWaterClean(this)
          !
          ! !DESCRIPTION:
          ! Deallocate memory
          !
          ! !ARGUMENTS:
          class(sectorwater_type), intent(inout) :: this
          !
          ! !LOCAL VARIABLES:
          
          character(len=*), parameter :: subname = 'Clean'
          !-----------------------------------------------------------------------
     
          deallocate(this%input_mon_dom_withd_grc)
          deallocate(this%input_mon_dom_cons_grc)
          deallocate(this%dom_withd_grc)
          deallocate(this%dom_cons_grc)
          deallocate(this%dom_withd_actual_grc)
          deallocate(this%dom_cons_actual_grc)
          deallocate(this%dom_rf_actual_grc)
          
          deallocate(this%input_mon_liv_withd_grc)
          deallocate(this%input_mon_liv_cons_grc)
          deallocate(this%liv_withd_grc)
          deallocate(this%liv_cons_grc)
          deallocate(this%liv_withd_actual_grc)
          deallocate(this%liv_cons_actual_grc)
          deallocate(this%liv_rf_actual_grc)
          
          deallocate(this%input_mon_elec_withd_grc)
          deallocate(this%input_mon_elec_cons_grc)
          deallocate(this%elec_withd_grc)
          deallocate(this%elec_cons_grc)
          deallocate(this%elec_withd_actual_grc)
          deallocate(this%elec_cons_actual_grc)
          deallocate(this%elec_rf_actual_grc)
          
          deallocate(this%input_mon_mfc_withd_grc)
          deallocate(this%input_mon_mfc_cons_grc)
          deallocate(this%mfc_withd_grc)
          deallocate(this%mfc_cons_grc)
          deallocate(this%mfc_withd_actual_grc)
          deallocate(this%mfc_cons_actual_grc)
          deallocate(this%mfc_rf_actual_grc)
          
          deallocate(this%input_mon_min_withd_grc)
          deallocate(this%input_mon_min_cons_grc)
          deallocate(this%min_withd_grc)
          deallocate(this%min_cons_grc)
          deallocate(this%min_withd_actual_grc)
          deallocate(this%min_cons_actual_grc)
          deallocate(this%min_rf_actual_grc)

          deallocate(this%sectorwater_total_actual_withd)


     end subroutine sectorWaterClean
 
 
     ! ========================================================================
     ! Science routines
     ! ========================================================================
     
     subroutine ReadSectorWaterData (this, bounds, year, mon)
          !
          ! !DESCRIPTION:
          ! read the input data, withdrawal and consumption, for all sectors and for the current month 
          !
          ! !USES:
          use fileutils       , only : getfil
          use clm_varcon      , only : nameg
          use ncdio_pio       , only : file_desc_t
          use spmdMod         , only : masterproc
          use netcdf
          
          ! !ARGUMENTS:
          class(sectorwater_type), intent(inout) :: this
          type(bounds_type)      , intent(in)    :: bounds
          integer                , intent(in)    :: year    ! current year (e.g. 2000)
          integer                , intent(in)    :: mon     ! month (1, ..., 12) for nstep+1
     
          !
          ! !LOCAL VARIABLES:
          type(file_desc_t) :: ncid               ! netcdf id
          integer :: start_year_input, end_year_input, current_line_number, i, read_status

          integer :: ier                          ! error code
          integer :: g                            ! gridcell index
          logical :: readvar
          real(r8), pointer :: mon_dom_withd(:)   ! monthly domestic withdrawal read from input files
          real(r8), pointer :: mon_dom_cons(:)    ! monthly domestic consumption read from input files
          real(r8), pointer :: mon_liv_withd(:)   ! monthly livestock withdrawal read from input files
          real(r8), pointer :: mon_liv_cons(:)    ! monthly livestock consumption read from input files
          real(r8), pointer :: mon_elec_withd(:)  ! monthly thermoelectric withdrawal read from input files
          real(r8), pointer :: mon_elec_cons(:)   ! monthly thermoelectric consumption read from input files
          real(r8), pointer :: mon_mfc_withd(:)   ! monthly manufacturing withdrawal read from input files
          real(r8), pointer :: mon_mfc_cons(:)    ! monthly manufacturing consumption read from input files
          real(r8), pointer :: mon_min_withd(:)   ! monthly mining withdrawal read from input files
          real(r8), pointer :: mon_min_cons(:)    ! monthly mining consumption read from input files
 
          character(len=256) :: current_line
          character(len=256) :: current_year_input_data ! path for the sectorwater input data for current year
          character(len=256) :: locfn             ! local file name
          character(len=20)  :: yearErrMessage
          character(len=20)  :: string_year
          character(len=32)  :: subname = 'ReadSectorWaterData'
          !-----------------------------------------------------------------------
 
          if (masterproc) then
               write (iulog,*) 'Attempting to read sectoral water usage data for current month .....'
          end if
          
 
          allocate(&
                    mon_dom_withd(bounds%begg:bounds%endg), &
                    mon_dom_cons(bounds%begg:bounds%endg), &
                    mon_liv_withd(bounds%begg:bounds%endg), &
                    mon_liv_cons(bounds%begg:bounds%endg), &
                    mon_elec_withd(bounds%begg:bounds%endg), &
                    mon_elec_cons(bounds%begg:bounds%endg), &
                    mon_mfc_withd(bounds%begg:bounds%endg), &
                    mon_mfc_cons(bounds%begg:bounds%endg), &
                    mon_min_withd(bounds%begg:bounds%endg), &
                    mon_min_cons(bounds%begg:bounds%endg), &
                    stat=ier)
                    
          ! Open the input .txt file:
          open(unit=10, file=this%params%path_sectorwater_input_data, status='old')
          ! Read first line to get start_year_input and end_year_input
          read(10,*,IOSTAT=i) start_year_input, end_year_input
          ! Check if current year is withing the input data limits
          if (year > end_year_input) then
               write(yearErrMessage, '(I0)') year
               yearErrMessage = trim(yearErrMessage)
               call endrun(msg='Error: there is no sector water demand data for current year. Please update the sector water input .txt file with the path to water demand inputs for current year '//yearErrMessage//errMsg(sourcefile, __LINE__)) 
          end if
          if (year < start_year_input) then
               write(yearErrMessage, '(I0)') year
               yearErrMessage = trim(yearErrMessage)
               call endrun(msg='Error: there is no sector water demand data for current year. Please update the sector water input .txt file with the path to water demand inputs for current year '//yearErrMessage//errMsg(sourcefile, __LINE__)) 
          end if

          ! Compute the current line number
          current_line_number = year - start_year_input + 2
          ! Rewind file
          rewind(10)

          i = 0
          do
               read(10,'(A)', ADVANCE='NO',IOSTAT=read_status) current_line
               i = i + 1
               if (i == current_line_number) exit
          end do
          ! convert year to string_year
          write(string_year, '(I0)') year
          ! Check if the current_line has the string_year in it
          ! If not, it may be that user didn't included the path for this year, or the paths are not provided in consecutive order by year
          if (index(current_line, trim(string_year)).eq.0) Then
               call endrun(msg=' ERROR: for current year, the sector water input path is '//current_line//". This path name do not contain the current year in its name (current year is "//trim(string_year)//")."//errMsg(sourcefile, __LINE__))
          end if
          ! Assign current line to path_current_year_input_data
          current_year_input_data = current_line
          ! Close the input .txt file
          Close(10)

          ! Determine necessary indices
          call getfil(current_year_input_data, locfn, 0)
          call ncd_pio_openfile (ncid, trim(locfn), 0)
          
          call ncd_io(ncid=ncid, varname='withd_dom', flag='read', data=mon_dom_withd, &
                    dim1name=nameg, nt=mon, readvar=readvar)
          if (.not. readvar) call endrun(msg=' ERROR: withd_dom NOT on surfdata file'//errMsg(sourcefile, __LINE__))   

          call ncd_io(ncid=ncid, varname='cons_dom', flag='read', data=mon_dom_cons, &
                    dim1name=nameg, nt=mon, readvar=readvar)
          if (.not. readvar) call endrun(msg=' ERROR: cons_dom NOT on surfdata file'//errMsg(sourcefile, __LINE__))           
          
          call ncd_io(ncid=ncid, varname='withd_liv', flag='read', data=mon_liv_withd, &
                    dim1name=nameg, nt=mon, readvar=readvar)
          if (.not. readvar) call endrun(msg=' ERROR: withd_liv NOT on surfdata file'//errMsg(sourcefile, __LINE__))   

          call ncd_io(ncid=ncid, varname='cons_liv', flag='read', data=mon_liv_cons, &
                    dim1name=nameg, nt=mon, readvar=readvar)
          if (.not. readvar) call endrun(msg=' ERROR: cons_liv NOT on surfdata file'//errMsg(sourcefile, __LINE__))   
          
          call ncd_io(ncid=ncid, varname='withd_elec', flag='read', data=mon_elec_withd, &
                    dim1name=nameg, nt=mon, readvar=readvar)
          if (.not. readvar) call endrun(msg=' ERROR: withd_elec NOT on surfdata file'//errMsg(sourcefile, __LINE__))           

          call ncd_io(ncid=ncid, varname='cons_elec', flag='read', data=mon_elec_cons, &
                    dim1name=nameg, nt=mon, readvar=readvar)     
          if (.not. readvar) call endrun(msg=' ERROR: cons_elec NOT on surfdata file'//errMsg(sourcefile, __LINE__))  

          call ncd_io(ncid=ncid, varname='withd_mfc', flag='read', data=mon_mfc_withd, &
                    dim1name=nameg, nt=mon, readvar=readvar)
          if (.not. readvar) call endrun(msg=' ERROR: withd_mfc NOT on surfdata file'//errMsg(sourcefile, __LINE__))   

          call ncd_io(ncid=ncid, varname='cons_mfc', flag='read', data=mon_mfc_cons, &
                    dim1name=nameg, nt=mon, readvar=readvar)
          if (.not. readvar) call endrun(msg=' ERROR: cons_mfc NOT on surfdata file'//errMsg(sourcefile, __LINE__))   

          call ncd_io(ncid=ncid, varname='withd_min', flag='read', data=mon_min_withd, &
                    dim1name=nameg, nt=mon, readvar=readvar)
          if (.not. readvar) call endrun(msg=' ERROR: withd_min NOT on surfdata file'//errMsg(sourcefile, __LINE__))   

          call ncd_io(ncid=ncid, varname='cons_min', flag='read', data=mon_min_cons, &
                    dim1name=nameg, nt=mon, readvar=readvar)
          if (.not. readvar) call endrun(msg=' ERROR: cons_min NOT on surfdata file'//errMsg(sourcefile, __LINE__))   

          call ncd_pio_closefile(ncid)
 
          do g = bounds%begg,bounds%endg
               this%input_mon_dom_withd_grc(g) = mon_dom_withd(g)
               this%input_mon_dom_cons_grc(g)  = mon_dom_cons(g)
          
               this%input_mon_liv_withd_grc(g) = mon_liv_withd(g)
               this%input_mon_liv_cons_grc(g)  = mon_liv_cons(g)
               
               this%input_mon_elec_withd_grc(g) = mon_elec_withd(g)
               this%input_mon_elec_cons_grc(g)  = mon_elec_cons(g)
               
               this%input_mon_mfc_withd_grc(g) = mon_mfc_withd(g)
               this%input_mon_mfc_cons_grc(g)  = mon_mfc_cons(g)
               
               this%input_mon_min_withd_grc(g) = mon_min_withd(g)
               this%input_mon_min_cons_grc(g)  = mon_min_cons(g)
               
          end do

          deallocate(mon_dom_withd, mon_dom_cons, mon_liv_withd, mon_liv_cons, mon_elec_withd, mon_elec_cons, mon_mfc_withd, mon_mfc_cons, mon_min_withd, mon_min_cons)

     endsubroutine ReadSectorWaterData
 
     subroutine CalcSectorWaterNeeded(this, bounds, volr, rof_prognostic)
 
          use shr_const_mod      , only : SHR_CONST_TKFRZ
          use clm_time_manager   , only : get_curr_date, is_end_curr_month, get_curr_days_per_year
          !
          ! !ARGUMENTS:
          class(sectorwater_type) , intent(inout) :: this
          type(bounds_type)       , intent(in)    :: bounds
          
          ! river water volume (m3) (ignored if rof_prognostic is .false.)
          real(r8), intent(in) :: volr(bounds%begg:bounds%endg)
          
          ! whether we're running with a prognostic ROF component; this is needed to determine
          ! whether we can limit demand based on river volume.
          logical, intent(in) :: rof_prognostic
          
          !
          ! !LOCAL VARIABLES:
          integer :: g       ! gridcell index
          integer :: year    ! year (0, ...) for nstep+1
          integer :: mon     ! month (1, ..., 12) for nstep+1
          integer :: day     ! day of month (1, ..., 31) for nstep+1
          integer :: sec     ! seconds into current date for nstep+1

          integer :: first_read ! variable to do first read, in future I may prefer to make a subroutine is_beg_curr_month to avoid exception for first reading

          real(r8) :: dayspyr                                    ! days per year
          real(r8) :: dayspm                                     ! days per month
          real(r8) :: secs_per_day                               ! seconds per day
          real(r8) :: flux_transform_from_monthly_to_second      ! factor to transform the demand from mm/month to mm/s for the given day

          real(r8) :: dom_demand(bounds%begg:bounds%endg)
          real(r8) :: dom_demand_volr_limited(bounds%begg:bounds%endg)
          
          real(r8) :: dom_consumption(bounds%begg:bounds%endg)
          real(r8) :: dom_consumption_volr_limited(bounds%begg:bounds%endg)
          
          real(r8) :: liv_demand(bounds%begg:bounds%endg)
          real(r8) :: liv_demand_volr_limited(bounds%begg:bounds%endg)
          
          real(r8) :: liv_consumption(bounds%begg:bounds%endg)
          real(r8) :: liv_consumption_volr_limited(bounds%begg:bounds%endg)
          
          real(r8) :: elec_demand(bounds%begg:bounds%endg)
          real(r8) :: elec_demand_volr_limited(bounds%begg:bounds%endg)
          
          real(r8) :: elec_consumption(bounds%begg:bounds%endg)
          real(r8) :: elec_consumption_volr_limited(bounds%begg:bounds%endg)
          
          real(r8) :: mfc_demand(bounds%begg:bounds%endg)
          real(r8) :: mfc_demand_volr_limited(bounds%begg:bounds%endg)
          
          real(r8) :: mfc_consumption(bounds%begg:bounds%endg)
          real(r8) :: mfc_consumption_volr_limited(bounds%begg:bounds%endg)
          
          real(r8) :: min_demand(bounds%begg:bounds%endg)
          real(r8) :: min_demand_volr_limited(bounds%begg:bounds%endg)
          
          real(r8) :: min_consumption(bounds%begg:bounds%endg)
          real(r8) :: min_consumption_volr_limited(bounds%begg:bounds%endg)
          
          ! Whether we should limit deficits by available volr
          logical :: limit_sectorwater
          
          character(len=*), parameter :: subname = 'CalcSectorWaterNeeded'
          !-----------------------------------------------------------------------
 
          ! Get current date
          call get_curr_date(year, mon, day, sec)
          dayspyr = get_curr_days_per_year()
          dayspm  = dayspyr/12._r8
          first_read = 1
          
          flux_transform_from_monthly_to_second = (1._r8/dayspm)/86400._r8
          
          
          ! For the first month of the simulation we need to read the data at the beginning of the month
          ! But for all the remaining months, the information is updated at the end of the month.
          if (first_read == 1) then 
               call this%ReadSectorWaterData(bounds, year, mon)
               first_read = 2
          endif

          ! Read input for the new month if end of current month
          if (is_end_curr_month()) then
               call this%ReadSectorWaterData(bounds, year, mon)
          endif
          
          ! Compute demand [mm]
          ! First initialize demand to 0 everywhere;
          dom_demand(bounds%begg:bounds%endg)       = 0._r8
          dom_consumption(bounds%begg:bounds%endg)  = 0._r8
          
          liv_demand(bounds%begg:bounds%endg)       = 0._r8
          liv_consumption(bounds%begg:bounds%endg)  = 0._r8
          
          elec_demand(bounds%begg:bounds%endg)      = 0._r8
          elec_consumption(bounds%begg:bounds%endg) = 0._r8
          
          mfc_demand(bounds%begg:bounds%endg)       = 0._r8
          mfc_consumption(bounds%begg:bounds%endg)  = 0._r8
          
          min_demand(bounds%begg:bounds%endg)       = 0._r8
          min_consumption(bounds%begg:bounds%endg)  = 0._r8
 
          do g = bounds%begg,bounds%endg
               dom_demand(g)       = this%input_mon_dom_withd_grc(g) *flux_transform_from_monthly_to_second
               dom_consumption(g)  = this%input_mon_dom_cons_grc(g)  *flux_transform_from_monthly_to_second
               
               liv_demand(g)       = this%input_mon_liv_withd_grc(g) *flux_transform_from_monthly_to_second
               liv_consumption(g)  = this%input_mon_liv_cons_grc(g)  *flux_transform_from_monthly_to_second
               
               elec_demand(g)      = this%input_mon_elec_withd_grc(g)*flux_transform_from_monthly_to_second
               elec_consumption(g) = this%input_mon_elec_cons_grc(g) *flux_transform_from_monthly_to_second
               
               mfc_demand(g)       = this%input_mon_mfc_withd_grc(g) *flux_transform_from_monthly_to_second
               mfc_consumption(g)  = this%input_mon_mfc_cons_grc(g)  *flux_transform_from_monthly_to_second
               
               min_demand(g)       = this%input_mon_min_withd_grc(g) *flux_transform_from_monthly_to_second
               min_consumption(g)  = this%input_mon_min_cons_grc(g)  *flux_transform_from_monthly_to_second
          
          end do ! end loop over gridcels
 
          ! Limit deficits by available volr, if desired. Note that we cannot do this limiting
          ! if running without a prognostic river model, since we need river volume to impose the limitation.
          limit_sectorwater = (this%params%limit_sectorwater_if_rof_enabled .and. rof_prognostic)
          if (limit_sectorwater) then
               call this%CalcSectorDemandVolrLimited( &
                    bounds = bounds, &
                    dom_demand = dom_demand(bounds%begg:bounds%endg), &
                    dom_consumption = dom_consumption(bounds%begg:bounds%endg), &
                    liv_demand = liv_demand(bounds%begg:bounds%endg), &
                    liv_consumption = liv_consumption(bounds%begg:bounds%endg), &
                    elec_demand = elec_demand(bounds%begg:bounds%endg), &
                    elec_consumption = elec_consumption(bounds%begg:bounds%endg), &
                    mfc_demand = mfc_demand(bounds%begg:bounds%endg), &
                    mfc_consumption = mfc_consumption(bounds%begg:bounds%endg), &
                    min_demand = min_demand(bounds%begg:bounds%endg), &
                    min_consumption = min_consumption(bounds%begg:bounds%endg), &
                    volr = volr(bounds%begg:bounds%endg), &
                    dom_demand_volr_limited = dom_demand_volr_limited(bounds%begg:bounds%endg), &
                    dom_consumption_volr_limited = dom_consumption_volr_limited(bounds%begg:bounds%endg), &
                    liv_demand_volr_limited = liv_demand_volr_limited(bounds%begg:bounds%endg), &
                    liv_consumption_volr_limited = liv_consumption_volr_limited(bounds%begg:bounds%endg), &
                    elec_demand_volr_limited = elec_demand_volr_limited(bounds%begg:bounds%endg), &
                    elec_consumption_volr_limited = elec_consumption_volr_limited(bounds%begg:bounds%endg), &
                    mfc_demand_volr_limited = mfc_demand_volr_limited(bounds%begg:bounds%endg), &
                    mfc_consumption_volr_limited = mfc_consumption_volr_limited(bounds%begg:bounds%endg), &
                    min_demand_volr_limited = min_demand_volr_limited(bounds%begg:bounds%endg), &
                    min_consumption_volr_limited = min_consumption_volr_limited(bounds%begg:bounds%endg))
          else
               dom_demand_volr_limited(bounds%begg:bounds%endg)       = dom_demand(bounds%begg:bounds%endg)
               dom_consumption_volr_limited(bounds%begg:bounds%endg)  = dom_consumption(bounds%begg:bounds%endg)
               
               liv_demand_volr_limited(bounds%begg:bounds%endg)       = liv_demand(bounds%begg:bounds%endg)
               liv_consumption_volr_limited(bounds%begg:bounds%endg)  = liv_consumption(bounds%begg:bounds%endg)
               
               elec_demand_volr_limited(bounds%begg:bounds%endg)      = elec_demand(bounds%begg:bounds%endg)
               elec_consumption_volr_limited(bounds%begg:bounds%endg) = elec_consumption(bounds%begg:bounds%endg)
               
               mfc_demand_volr_limited(bounds%begg:bounds%endg)       = mfc_demand(bounds%begg:bounds%endg)
               mfc_consumption_volr_limited(bounds%begg:bounds%endg)  = mfc_consumption(bounds%begg:bounds%endg)
               
               min_demand_volr_limited(bounds%begg:bounds%endg)       = min_demand(bounds%begg:bounds%endg)
               min_consumption_volr_limited(bounds%begg:bounds%endg)  = min_consumption(bounds%begg:bounds%endg)
          end if
 
          ! Convert demand to withdrawal rates [mm/s]
          ! Here it also seems like I could directly operate with the this%arrays instead of generating new ones (to check)
          do g = bounds%begg,bounds%endg
               ! Domestic
               this%dom_withd_grc(g)         = dom_demand(g)
               this%dom_withd_actual_grc(g)  = dom_demand_volr_limited(g)
               
               this%dom_cons_grc(g)          = dom_consumption(g)
               this%dom_cons_actual_grc(g)   = dom_consumption_volr_limited(g)
               
               this%dom_rf_actual_grc(g)     = this%dom_withd_actual_grc(g) - this%dom_cons_actual_grc(g)
               
               ! Livestock
               this%liv_withd_grc(g)         = liv_demand(g)
               this%liv_withd_actual_grc(g)  = liv_demand_volr_limited(g)
               
               this%liv_cons_grc(g)          = liv_consumption(g)
               this%liv_cons_actual_grc(g)   = liv_consumption_volr_limited(g)
               
               this%liv_rf_actual_grc(g)     = this%liv_withd_actual_grc(g) - this%liv_cons_actual_grc(g)
               
               ! Thermoelectric
               this%elec_withd_grc(g)        = elec_demand(g)
               this%elec_withd_actual_grc(g) = elec_demand_volr_limited(g)
               
               this%elec_cons_grc(g)         = elec_consumption(g)
               this%elec_cons_actual_grc(g)  = elec_consumption_volr_limited(g)
               
               this%elec_rf_actual_grc(g)    = this%elec_withd_actual_grc(g) - this%elec_cons_actual_grc(g)
               
               ! Manufacturing
               this%mfc_withd_grc(g)         = mfc_demand(g)
               this%mfc_withd_actual_grc(g)  = mfc_demand_volr_limited(g)
               
               this%mfc_cons_grc(g)          = mfc_consumption(g)
               this%mfc_cons_actual_grc(g)   = mfc_consumption_volr_limited(g)
               
               this%mfc_rf_actual_grc(g)     = this%mfc_withd_actual_grc(g) - this%mfc_cons_actual_grc(g)
               
               ! Mining
               this%min_withd_grc(g)         = min_demand(g)
               this%min_withd_actual_grc(g)  = min_demand_volr_limited(g)
               
               this%min_cons_grc(g)          = min_consumption(g)
               this%min_cons_actual_grc(g)   = min_consumption_volr_limited(g)
               
               this%min_rf_actual_grc(g)     = this%min_withd_actual_grc(g) - this%min_cons_actual_grc(g)

               ! Total actual withdrawal volume in m3 for all sectors during the current day
               ! Total actual sectoral withdrawal volume will be used to constrain how much water is available for irrigation taking into acount VOLR capacity (if limit on abstractions is active)
               this%sectorwater_total_actual_withd(g) =  (this%dom_withd_actual_grc(g) + this%liv_withd_actual_grc(g) + &
                                                          this%elec_withd_actual_grc(g) + this%mfc_withd_actual_grc(g) + &
                                                          this%min_withd_actual_grc(g)) * mm_to_m3_over_km2 * grc%area(g) * 86400._r8
          end do
          
     end subroutine CalcSectorWaterNeeded
 
 
     !-----------------------------------------------------------------------
     subroutine CalcSectorDemandVolrLimited(this, bounds, dom_demand, dom_consumption, liv_demand, liv_consumption, elec_demand, elec_consumption, &
     mfc_demand, mfc_consumption, min_demand, min_consumption, volr, dom_demand_volr_limited, dom_consumption_volr_limited, liv_demand_volr_limited, &
     liv_consumption_volr_limited, elec_demand_volr_limited, elec_consumption_volr_limited, mfc_demand_volr_limited, &
     mfc_consumption_volr_limited, min_demand_volr_limited, min_consumption_volr_limited)
          ! !USES:
          !
          ! !ARGUMENTS:
          class(sectorwater_type) , intent(in) :: this
          type(bounds_type)      , intent(in) :: bounds
          
          real(r8), intent(in)  :: dom_demand( bounds%begg:bounds%endg)
          real(r8), intent(in)  :: dom_consumption( bounds%begg:bounds%endg)
          
          real(r8), intent(in)  :: liv_demand( bounds%begg:bounds%endg)
          real(r8), intent(in)  :: liv_consumption( bounds%begg:bounds%endg)
          
          real(r8), intent(in)  :: elec_demand( bounds%begg:bounds%endg)
          real(r8), intent(in)  :: elec_consumption( bounds%begg:bounds%endg)
          
          real(r8), intent(in)  :: mfc_demand( bounds%begg:bounds%endg)
          real(r8), intent(in)  :: mfc_consumption( bounds%begg:bounds%endg)
          
          real(r8), intent(in)  :: min_demand( bounds%begg:bounds%endg)
          real(r8), intent(in)  :: min_consumption( bounds%begg:bounds%endg)
 
 
          real(r8), intent(in)  :: volr( bounds%begg:bounds%endg) ! river water volume [m3]
          
          real(r8), intent(out) :: dom_demand_volr_limited( bounds%begg:bounds%endg)
          real(r8), intent(out) :: dom_consumption_volr_limited( bounds%begg:bounds%endg)
          
          real(r8), intent(out) :: liv_demand_volr_limited( bounds%begg:bounds%endg)
          real(r8), intent(out) :: liv_consumption_volr_limited( bounds%begg:bounds%endg)
          
          real(r8), intent(out) :: elec_demand_volr_limited( bounds%begg:bounds%endg)
          real(r8), intent(out) :: elec_consumption_volr_limited( bounds%begg:bounds%endg)
          
          real(r8), intent(out) :: mfc_demand_volr_limited( bounds%begg:bounds%endg)
          real(r8), intent(out) :: mfc_consumption_volr_limited( bounds%begg:bounds%endg)
          
          real(r8), intent(out) :: min_demand_volr_limited( bounds%begg:bounds%endg)
          real(r8), intent(out) :: min_consumption_volr_limited( bounds%begg:bounds%endg)
 
          !
          ! !LOCAL VARIABLES:
          integer  :: g  ! gridcell index
          real(r8) :: available_volr ! volr available for withdrawal [m3]
          real(r8) :: max_demand_supported_by_volr ! [kg/m2] [i.e., mm]
          
          ! ratio of demand_volr_limited to demand for each grid cell
          real(r8) :: dom_demand_limited_ratio_grc(bounds%begg:bounds%endg)
          real(r8) :: liv_demand_limited_ratio_grc(bounds%begg:bounds%endg)
          real(r8) :: elec_demand_limited_ratio_grc(bounds%begg:bounds%endg)
          real(r8) :: mfc_demand_limited_ratio_grc(bounds%begg:bounds%endg)
          real(r8) :: min_demand_limited_ratio_grc(bounds%begg:bounds%endg)
 
 
          character(len=*), parameter :: subname = 'CalcSectorDemandVolrLimited'
          !-----------------------------------------------------------------------
          

          do g = bounds%begg, bounds%endg
               if (volr(g) > 0._r8) then
                    available_volr = volr(g) * (1._r8 - this%params%sectorwater_river_volume_threshold)
                    max_demand_supported_by_volr = (available_volr / grc%area(g) * m3_over_km2_to_mm)
               else
                    ! Ensure that negative volr is treated the same as 0 volr
                    max_demand_supported_by_volr = 0._r8
               end if
               
               ! I think the algorithm is potentially too conservative
               ! I will need to check river discharge when there is high amount of unsatisfied demand
               ! The reason why I am saying this is because we compare the expected demand for the entire day with the current volume available in the river...
               ! It would make sense to use such an algorithm under the condition that volr does not change much during a day for a given gridcell
               ! But if this is not the case we maybe underestimate the amount of water available for usage.
               ! If this would be done once a day, then no problem, but we do it at each time step
               ! This means that volr get 
               if (dom_demand(g) * 86400.0 > max_demand_supported_by_volr) then
                    ! inadequate river storage, adjust demand
                    dom_demand_limited_ratio_grc(g)  = max_demand_supported_by_volr / (dom_demand(g) * 86400.0)
                    liv_demand_limited_ratio_grc(g)  = 0._r8
                    elec_demand_limited_ratio_grc(g) = 0._r8
                    mfc_demand_limited_ratio_grc(g)  = 0._r8
                    min_demand_limited_ratio_grc(g)  = 0._r8
     
               else if (liv_demand(g) * 86400.0 > (max_demand_supported_by_volr - (dom_demand(g) * 86400.0))) then
                    dom_demand_limited_ratio_grc(g)  = 1._r8
                    liv_demand_limited_ratio_grc(g)  = (max_demand_supported_by_volr - (dom_demand(g) * 86400.0)) / (liv_demand(g) * 86400.0)
                    elec_demand_limited_ratio_grc(g) = 0._r8
                    mfc_demand_limited_ratio_grc(g)  = 0._r8
                    min_demand_limited_ratio_grc(g)  = 0._r8
               else if (elec_demand(g) * 86400.0 > (max_demand_supported_by_volr - (dom_demand(g) + liv_demand(g))* 86400.0 )) then
                    dom_demand_limited_ratio_grc(g)  = 1._r8
                    liv_demand_limited_ratio_grc(g)  = 1._r8
                    elec_demand_limited_ratio_grc(g) = (max_demand_supported_by_volr - (dom_demand(g) + liv_demand(g))* 86400.0)  / (elec_demand(g) * 86400.0)
                    mfc_demand_limited_ratio_grc(g)  = 0._r8
                    min_demand_limited_ratio_grc(g)  = 0._r8
               else if (mfc_demand(g) * 86400.0 > (max_demand_supported_by_volr - (dom_demand(g) + liv_demand(g) + elec_demand(g)) * 86400.0)) then
                    dom_demand_limited_ratio_grc(g)  = 1._r8
                    liv_demand_limited_ratio_grc(g)  = 1._r8
                    elec_demand_limited_ratio_grc(g) = 1._r8
                    mfc_demand_limited_ratio_grc(g)  = (max_demand_supported_by_volr - (dom_demand(g) + liv_demand(g) + elec_demand(g)) * 86400.0) / (mfc_demand(g) * 86400.0)
                    min_demand_limited_ratio_grc(g)  = 0._r8
               
               else if (min_demand(g) * 86400.0 > (max_demand_supported_by_volr - (dom_demand(g) + liv_demand(g) + elec_demand(g) + mfc_demand(g)) * 86400.0)) then
                    dom_demand_limited_ratio_grc(g)  = 1._r8
                    liv_demand_limited_ratio_grc(g)  = 1._r8
                    elec_demand_limited_ratio_grc(g) = 1._r8
                    mfc_demand_limited_ratio_grc(g)  = 1._r8
                    min_demand_limited_ratio_grc(g)  = (max_demand_supported_by_volr - (dom_demand(g) + liv_demand(g) + elec_demand(g) + mfc_demand(g)) * 86400.0)  / (min_demand(g) * 86400.0)
               
               else
                    dom_demand_limited_ratio_grc(g)  = 1._r8
                    liv_demand_limited_ratio_grc(g)  = 1._r8
                    elec_demand_limited_ratio_grc(g) = 1._r8
                    mfc_demand_limited_ratio_grc(g)  = 1._r8
                    min_demand_limited_ratio_grc(g)  = 1._r8
               end if
          end do
 
          dom_demand_volr_limited(bounds%begg:bounds%endg) = 0._r8
          dom_consumption_volr_limited(bounds%begg:bounds%endg) = 0._r8
          liv_demand_volr_limited(bounds%begg:bounds%endg) = 0._r8
          liv_consumption_volr_limited(bounds%begg:bounds%endg) = 0._r8
          elec_demand_volr_limited(bounds%begg:bounds%endg) = 0._r8
          elec_consumption_volr_limited(bounds%begg:bounds%endg) = 0._r8
          mfc_demand_volr_limited(bounds%begg:bounds%endg) = 0._r8
          mfc_consumption_volr_limited(bounds%begg:bounds%endg) = 0._r8
          min_demand_volr_limited(bounds%begg:bounds%endg) = 0._r8
          min_consumption_volr_limited(bounds%begg:bounds%endg) = 0._r8

          do g = bounds%begg, bounds%endg
               dom_demand_volr_limited(g) = dom_demand(g) * dom_demand_limited_ratio_grc(g)
               dom_consumption_volr_limited(g) = dom_consumption(g) * dom_demand_limited_ratio_grc(g)
               
               liv_demand_volr_limited(g) = liv_demand(g) * liv_demand_limited_ratio_grc(g)
               liv_consumption_volr_limited(g) = liv_consumption(g) * liv_demand_limited_ratio_grc(g)
               
               elec_demand_volr_limited(g) = elec_demand(g) * elec_demand_limited_ratio_grc(g)
               elec_consumption_volr_limited(g) = elec_consumption(g) * elec_demand_limited_ratio_grc(g)
               
               mfc_demand_volr_limited(g) = mfc_demand(g) * mfc_demand_limited_ratio_grc(g)
               mfc_consumption_volr_limited(g) = mfc_consumption(g) * mfc_demand_limited_ratio_grc(g)
               
               min_demand_volr_limited(g) = min_demand(g) * min_demand_limited_ratio_grc(g)
               min_consumption_volr_limited(g) = min_consumption(g) * min_demand_limited_ratio_grc(g)
          
          end do
     
     end subroutine CalcSectorDemandVolrLimited
 
 end module SectorWaterMod