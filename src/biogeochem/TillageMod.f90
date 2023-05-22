module TillageMod
  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module for soil tillage.
  !
  ! !USES:
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use shr_log_mod    , only : errMsg => shr_log_errMsg
  use abortutils     , only : endrun
  use clm_varctl     , only : iulog
  ! TODO: Confirm that i_soil1->i_act_som, i_soil2->i_slo_som, i_soil3->i_pas_som
  use clm_varpar     , only : i_litr2, i_litr3, ndecomp_pools
  !
  implicit none
  private
  ! !PUBLIC MEMBER PROCEDURES
  public :: tillage_init
  public :: get_do_tillage
  public :: get_apply_tillage_multipliers
  !
  ! !PRIVATE DATA MEMBERS
  logical  :: do_tillage_low   ! Do low-intensity tillage?
  logical  :: do_tillage_high  ! Do high-intensity tillage?
  logical  :: use_original_tillage ! Use get_tillage_multipliers_orig?
  real(r8), pointer :: tillage_mults(:,:)


!==============================================================================
contains
!==============================================================================

  subroutine tillage_init(bounds)
    !
    ! Read namelist parameters and allocate variables related to tillage
    !
    ! !USES:
    use spmdMod        , only : masterproc
    use controlMod     , only : NLFilename
    use clm_nlUtilsMod , only : find_nlgroup_name
    use shr_mpi_mod    , only : shr_mpi_bcast
    use decompMod      , only : bounds_type
    !
    ! !ARGUMENTS
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES
    integer                :: begp, endp   ! indices for allocating patch dimension
    integer                :: nu_nml       ! unit for namelist file
    integer                :: nml_error    ! namelist i/o error flag
    integer                :: mpicom       ! MPI communicator
    character(*), parameter :: subname = "('tillage_init')"

    namelist /tillage_nml/    &
        do_tillage_low,       &
        do_tillage_high,      &
        use_original_tillage

    ! Default values
    do_tillage_low = .false.
    do_tillage_high = .false.
    use_original_tillage = .false.

    ! Read tillage namelist
    if (masterproc) then
        open(newunit=nu_nml, file=trim(NLFilename), status='old', iostat=nml_error )
        call find_nlgroup_name(nu_nml, 'tillage_nml', status=nml_error)
        if (nml_error == 0) then
           read(nu_nml, nml=tillage_nml, iostat=nml_error)
           if (nml_error /= 0) then
              call endrun(subname // ':: ERROR reading tillage namelist')
           end if
        else
           call endrun(subname // ':: ERROR finding tillage namelist')
        end if
        close(nu_nml)
     endif
     call shr_mpi_bcast(do_tillage_low, mpicom)
     call shr_mpi_bcast(do_tillage_high , mpicom)
     call shr_mpi_bcast(use_original_tillage , mpicom)

     if (masterproc) then
        write(iulog,*) ' '
        write(iulog,*) 'tillage settings:'
        write(iulog,*) '  do_tillage_low  = ',do_tillage_low
        write(iulog,*) '  do_tillage_high   = ',do_tillage_high
        write(iulog,*) '  use_original_tillage   = ',use_original_tillage
     endif

     if (do_tillage_low .and. do_tillage_high) then
        call endrun(subname // ':: ERROR do_tillage_low and do_tillage_high are mutually exclusive')
     endif

     if (do_tillage_high) then
        call endrun(subname // ':: ERROR high-intensity tillage not yet set up')
     endif

     ! Allocate variables
     if (get_do_tillage()) then
        begp = bounds%begp
        endp = bounds%endp
        allocate(tillage_mults(begp:endp, ndecomp_pools)) ; tillage_mults(:,:) = 1.0_r8
     end if

  end subroutine tillage_init


  function get_do_tillage()
    logical :: get_do_tillage
    get_do_tillage = do_tillage_low .or. do_tillage_high
  end function get_do_tillage


  subroutine get_tillage_multipliers_orig(idop, num_soilp, filter_soilp, i_act_som, i_slo_som, i_pas_som)
    ! !DESCRIPTION:
    !
    !  Get the cultivation effective multiplier if prognostic crops are on and
    !  cultivation is turned on. Created by Sam Levis, modified by Michael Graham
    !  to use days past planting.
    !
    ! !USES:
    use clm_time_manager, only : get_curr_calday, get_curr_days_per_year
    use pftconMod       , only : npcropmin, ntmp_corn, nirrig_tmp_corn, ntmp_soybean, nirrig_tmp_soybean
    use PatchType       , only : patch
    ! !ARGUMENTS:
    integer          , intent(in) :: idop(:) ! patch day of planting
    integer          , intent(in) :: num_soilp          ! number of soil pfts in filter
    integer          , intent(in) :: filter_soilp(:)    ! filter for soil pfts
    integer          , intent(in) :: i_act_som, i_slo_som, i_pas_som  ! indices for soil organic matter pools
    !
    ! !LOCAL VARIABLES:
    !
    integer :: fp, p, c, g          ! Indices
    integer :: day                  ! julian day
    integer :: idpp                 ! days past planting
    real(r8) dayspyr                ! days per year
    !-----------------------------------------------------------------------
        
    !get info from externals
    day = get_curr_calday()
    dayspyr = get_curr_days_per_year()               !Add by MWG for IDPP-based routine

    do fp = 1,num_soilp
        p = filter_soilp(fp)
        c = patch%column(p)
        g = patch%gridcell(p)
    
        ! days past planting may determine harvest/tillage
        ! SSR: Unused!
        if (day >= idop(p)) then
            idpp = day - idop(p)
        else
            idpp = int(dayspyr) + day - idop(p)
        end if
    
        ! -----------------------------------------------------
        ! 3) assigning cultivation practices and mapping to the
        !    effect on soil C decomposition
        ! -----------------------------------------------------
        ! info from DAYCENT (Melannie Hartman CSU)
        ! temp. cereals: P 30 d bef, C 15 d bef, D on day of planting
        ! corn, soy    : P           C           D           & HW-7 30 d aftr

        if (day < idop(p)) then
            tillage_mults(c,:) = 1._r8
        else if (day >= idop(p) .and. day < idop(p)+15) then ! based on Point Chisel Tandem Disk multipliers
            tillage_mults(c,:) = 1._r8
            if (patch%itype(p) >= npcropmin) then
                tillage_mults(c,i_litr2) = 1.50_r8 !high 1.80,low 1.50
                tillage_mults(c,i_litr3) = 1.50_r8 !high 1.80,low 1.50
                tillage_mults(c,i_act_som) = 1.00_r8 !high 1.20,low 1.00
                tillage_mults(c,i_slo_som) = 3.00_r8 !high 4.80,low 3.00
                tillage_mults(c,i_pas_som) = 3.00_r8 !high 4.80,low 3.00
            end if
        else if (day >= idop(p)+15 .and. day < idop(p)+45) then ! based on Field and Row Cultivator multipliers
            tillage_mults(c,:) = 1._r8
            if (patch%itype(p) >= npcropmin) then
                tillage_mults(c,i_litr2) = 1.50_r8 !high 1.50,low 1.50
                tillage_mults(c,i_litr3) = 1.50_r8 !high 1.50,low 1.50
                tillage_mults(c,i_act_som) = 1.00_r8 !high 1.00,low 1.00
                tillage_mults(c,i_slo_som) = 1.60_r8 !high 3.50,low 1.60
                tillage_mults(c,i_pas_som) = 1.60_r8 !high 3.50,low 1.60
            end if
        else if (day >= idop(p)+45 .and. day <idop(p)+75) then ! based on Rod Weed Row Planter
            tillage_mults(c,:) = 1._r8
            if (patch%itype(p) >= npcropmin) then
                tillage_mults(c,i_litr2) = 1.10_r8 !high 1.10,low 1.10
                tillage_mults(c,i_litr3) = 1.10_r8 !high 1.10,low 1.10
                tillage_mults(c,i_act_som) = 1.00_r8 !high 1.00,low 1.00
                tillage_mults(c,i_slo_som) = 1.30_r8 !high 2.50,low 1.30
                tillage_mults(c,i_pas_som) = 1.30_r8 !high 2.50,low 1.30
            end if
        else if (day >= idop(p)+75 .and. day < idop(p)+80) then ! June 14
            tillage_mults(c,:) = 1._r8
            if (patch%itype(p) == ntmp_corn      .or. &
              patch%itype(p) == nirrig_tmp_corn .or. &
              patch%itype(p) == ntmp_soybean   .or. &
              patch%itype(p) == nirrig_tmp_soybean      ) then
                tillage_mults(c,i_litr2) = 1.00_r8
                tillage_mults(c,i_litr3) = 1.00_r8
                tillage_mults(c,i_act_som) = 1.00_r8
                tillage_mults(c,i_slo_som) = 1.00_r8
                tillage_mults(c,i_pas_som) = 1.00_r8
            end if
        else if (day >= idop(p)+80) then ! July 14
            tillage_mults(c,:) = 1._r8
        end if
    enddo
    
  end subroutine get_tillage_multipliers_orig


  subroutine get_tillage_multipliers_new(idop, num_soilp, filter_soilp, i_act_som, i_slo_som, i_pas_som)
    ! !DESCRIPTION:
    !
    !  Get the cultivation effective multiplier if prognostic crops are on and
    !  cultivation is turned on. Based on get_tillage_multipliers_orig().
    !  Changed by Sam Rabin:
    !  - Actually use idpp (corrected for crossing new year) instead of idop+N
    !
    ! !USES:
    use clm_time_manager, only : get_curr_calday, get_curr_days_per_year
    use pftconMod       , only : npcropmin, ntmp_corn, nirrig_tmp_corn, ntmp_soybean, nirrig_tmp_soybean
    use PatchType       , only : patch
    ! !ARGUMENTS:
    integer          , intent(in) :: idop(:) ! patch day of planting
    integer          , intent(in) :: num_soilp          ! number of soil pfts in filter
    integer          , intent(in) :: filter_soilp(:)    ! filter for soil pfts
    integer          , intent(in) :: i_act_som, i_slo_som, i_pas_som  ! indices for soil organic matter pools
    !
    ! !LOCAL VARIABLES:
    !
    integer :: fp, p, c, g          ! Indices
    integer :: day                  ! julian day
    integer :: idpp                 ! days past planting
    real(r8) dayspyr                ! days per year
    !-----------------------------------------------------------------------
    
    !get info from externals
    day = get_curr_calday()
    dayspyr = get_curr_days_per_year()               !Add by MWG for IDPP-based routine

    do fp = 1,num_soilp
        p = filter_soilp(fp)
        c = patch%column(p)
        g = patch%gridcell(p)
    
        ! days past planting may determine harvest/tillage
        if (day >= idop(p)) then
            idpp = day - idop(p)
        else
            idpp = int(dayspyr) + day - idop(p)
        end if
    
        ! -----------------------------------------------------
        ! 3) assigning cultivation practices and mapping to the
        !    effect on soil C decomposition
        ! -----------------------------------------------------
        ! info from DAYCENT (Melannie Hartman CSU)
        ! temp. cereals: P 30 d bef, C 15 d bef, D on day of planting
        ! corn, soy    : P           C           D           & HW-7 30 d aftr

        if (idpp < 0) then
            tillage_mults(c,:) = 1._r8
        else if (idpp < 15) then ! based on Point Chisel Tandem Disk multipliers
            tillage_mults(c,:) = 1._r8
            if (patch%itype(p) >= npcropmin) then
                tillage_mults(c,i_litr2) = 1.50_r8 !high 1.80,low 1.50
                tillage_mults(c,i_litr3) = 1.50_r8 !high 1.80,low 1.50
                tillage_mults(c,i_act_som) = 1.00_r8 !high 1.20,low 1.00
                tillage_mults(c,i_slo_som) = 3.00_r8 !high 4.80,low 3.00
                tillage_mults(c,i_pas_som) = 3.00_r8 !high 4.80,low 3.00
            end if
        else if (idpp < 45) then ! based on Field and Row Cultivator multipliers
            tillage_mults(c,:) = 1._r8
            if (patch%itype(p) >= npcropmin) then
                tillage_mults(c,i_litr2) = 1.50_r8 !high 1.50,low 1.50
                tillage_mults(c,i_litr3) = 1.50_r8 !high 1.50,low 1.50
                tillage_mults(c,i_act_som) = 1.00_r8 !high 1.00,low 1.00
                tillage_mults(c,i_slo_som) = 1.60_r8 !high 3.50,low 1.60
                tillage_mults(c,i_pas_som) = 1.60_r8 !high 3.50,low 1.60
            end if
        else if (idpp < 75) then ! based on Rod Weed Row Planter
            tillage_mults(c,:) = 1._r8
            if (patch%itype(p) >= npcropmin) then
                tillage_mults(c,i_litr2) = 1.10_r8 !high 1.10,low 1.10
                tillage_mults(c,i_litr3) = 1.10_r8 !high 1.10,low 1.10
                tillage_mults(c,i_act_som) = 1.00_r8 !high 1.00,low 1.00
                tillage_mults(c,i_slo_som) = 1.30_r8 !high 2.50,low 1.30
                tillage_mults(c,i_pas_som) = 1.30_r8 !high 2.50,low 1.30
            end if
        else if (idpp < 80) then ! June 14
            tillage_mults(c,:) = 1._r8
            if (patch%itype(p) == ntmp_corn     .or. &
              patch%itype(p) == nirrig_tmp_corn .or. &
              patch%itype(p) == ntmp_soybean    .or. &
              patch%itype(p) == nirrig_tmp_soybean      ) then
                tillage_mults(c,i_litr2) = 1.00_r8
                tillage_mults(c,i_litr3) = 1.00_r8
                tillage_mults(c,i_act_som) = 1.00_r8
                tillage_mults(c,i_slo_som) = 1.00_r8
                tillage_mults(c,i_pas_som) = 1.00_r8
            end if
        else ! July 14
            tillage_mults(c,:) = 1._r8
        end if
    enddo
    
  end subroutine get_tillage_multipliers_new


  ! Public interface to choose between and call either get_tillage_multipliers_orig() or get_tillage_multipliers_new()
  subroutine get_tillage_multipliers(idop, num_soilp, filter_soilp, i_act_som, i_slo_som, i_pas_som)
    ! !DESCRIPTION:
    !
    !  Public interface to choose between and call either original (buggy) or
    !  new (fixed), depending on use_original_tillage true or false.
    !
    ! !ARGUMENTS:
    integer          , intent(in) :: idop(:) ! patch day of planting
    integer          , intent(in) :: num_soilp          ! number of soil pfts in filter
    integer          , intent(in) :: filter_soilp(:)    ! filter for soil pfts
    integer          , intent(in) :: i_act_som, i_slo_som, i_pas_som  ! indices for soil organic matter pools

    if (use_original_tillage) then
        call get_tillage_multipliers_orig(idop, num_soilp, filter_soilp, i_act_som, i_slo_som, i_pas_som)
    else
        call get_tillage_multipliers_new(idop, num_soilp, filter_soilp, i_act_som, i_slo_som, i_pas_som)
    end if

  end subroutine get_tillage_multipliers


  subroutine get_apply_tillage_multipliers(idop, num_soilp, filter_soilp, num_soilc, filter_soilc, decomp_k, i_act_som, i_slo_som, i_pas_som)
    ! !DESCRIPTION:
    !
    ! Multiply decomposition rate constants by tillage coefficients.
    ! Written by Sam Rabin, based on original code by Michael Graham.
    !
    ! !USES
    !
    ! !ARGUMENTS:
    integer       , intent(in) :: idop(:) ! patch day of planting
    integer          , intent(in) :: num_soilp          ! number of soil pfts in filter
    integer          , intent(in) :: filter_soilp(:)    ! filter for soil pfts
    integer          , intent(in) :: num_soilc          ! number of soil columns in filter
    integer          , intent(in) :: filter_soilc(:)    ! filter for soil columns
    real(r8)         , intent(inout) :: decomp_k(:,:,:) ! Output: [real(r8) (:,:,:) ]  rate constant for decomposition (1./sec)
    integer          , intent(in) :: i_act_som, i_slo_som, i_pas_som  ! indices for soil organic matter pools
    !
    ! !LOCAL VARIABLES
    integer :: j, fc, c

    call get_tillage_multipliers(idop, num_soilp, filter_soilp, i_act_som, i_slo_som, i_pas_som)

    ! Top 5 layers (instead of all nlevdecomp) so that model only tills the top 26-40 cm
    ! of the soil surface, rather than whole soil - MWGraham
    do j = 1,5
        do fc = 1,num_soilc
           c = filter_soilc(fc)
           ! TODO: Loop through ALL pools, not just the ones that currently have non-1 values
           decomp_k(c,j,i_litr2) = decomp_k(c,j,i_litr2) * tillage_mults(c,i_litr2)
           decomp_k(c,j,i_litr3) = decomp_k(c,j,i_litr3) * tillage_mults(c,i_litr3) 
           decomp_k(c,j,i_act_som) = decomp_k(c,j,i_act_som) * tillage_mults(c,i_act_som)
           decomp_k(c,j,i_slo_som) = decomp_k(c,j,i_slo_som) * tillage_mults(c,i_slo_som)
           decomp_k(c,j,i_pas_som) = decomp_k(c,j,i_pas_som) * tillage_mults(c,i_pas_som)
        end do
    end do

  end subroutine get_apply_tillage_multipliers

end module TillageMod
