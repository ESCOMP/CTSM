module TillageMod
  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module for soil tillage.
  !
  ! !USES:
  use shr_kind_mod   , only : r8 => shr_kind_r8, CS => shr_kind_CS
  use shr_log_mod    , only : errMsg => shr_log_errMsg
  use abortutils     , only : endrun
  use clm_varctl     , only : iulog
  use clm_varpar     , only : ndecomp_pools
  use ColumnType     , only : col
  use PatchType      , only : patch
  !
  implicit none
  private
  ! !PUBLIC MEMBER PROCEDURES
  public :: readParams
  public :: get_do_tillage
  public :: get_apply_tillage_multipliers
  ! !PUBLIC DATA MEMBERS
  character(len=CS), public :: tillage_mode     ! off, low, high
  integer, parameter, public :: ntill_intensities_max = 2
  !
  ! !PRIVATE DATA MEMBERS
  integer             :: tillage_intensity
  integer, parameter  :: tillage_off = 0
  integer, parameter  :: tillage_low = 1
  integer, parameter  :: tillage_high = 2
  logical  :: use_original_tillage ! Use get_tillage_multipliers_orig?
  real(r8), pointer :: tillage_mults_allphases(:,:) ! (ndecomp_pools, ntill_stages_max)
  integer, parameter :: ntill_stages_max = 3 ! How many different tillage phases are there? (Not including all-1 phases.)

!==============================================================================
contains
!==============================================================================

  subroutine readParams_namelist()
    !
    ! Read namelist parameters related to tillage.
    !
    ! !USES:
    use spmdMod        , only : masterproc, mpicom
    use controlMod     , only : NLFilename
    use clm_nlUtilsMod , only : find_nlgroup_name
    use shr_mpi_mod    , only : shr_mpi_bcast
    !
    ! !LOCAL VARIABLES
    integer                :: nu_nml       ! unit for namelist file
    integer                :: nml_error    ! namelist i/o error flag
    character(*), parameter :: subname = "('readParams_namelist')"

    namelist /tillage_inparm/    &
        tillage_mode,            &
        use_original_tillage

    ! Default values
    tillage_mode = 'off'
    use_original_tillage = .false.

    ! Read tillage namelist
    if (masterproc) then
        open(newunit=nu_nml, file=trim(NLFilename), status='old', iostat=nml_error )
        call find_nlgroup_name(nu_nml, 'tillage_inparm', status=nml_error)
        if (nml_error == 0) then
           read(nu_nml, nml=tillage_inparm, iostat=nml_error)
           if (nml_error /= 0) then
              call endrun(subname // ':: ERROR reading tillage namelist')
           end if
        else
           call endrun(subname // ':: ERROR finding tillage namelist')
        end if
        close(nu_nml)
     endif
     call shr_mpi_bcast(tillage_mode, mpicom)
     call shr_mpi_bcast(use_original_tillage , mpicom)

     if (masterproc) then
        write(iulog,*) ' '
        write(iulog,*) 'tillage settings:'
        write(iulog,*) '  tillage_mode  = ',tillage_mode
        write(iulog,*) '  use_original_tillage   = ',use_original_tillage
     endif

     ! Assign these
     if (tillage_mode == "off") then
         tillage_intensity = tillage_off
     else if (tillage_mode == "low") then
         tillage_intensity = tillage_low
     else if (tillage_mode == "high") then
         tillage_intensity = tillage_high
     else
        call endrun(subname // ':: ERROR Unrecognized tillage_mode')
     end if

  end subroutine readParams_namelist


  subroutine readParams_netcdf(ncid)
    ! !DESCRIPTION:
    !
    ! Read paramfile parameters to be used in tillage.
    !
    ! !USES
    use ncdio_pio , only : file_desc_t, ncd_io
    use clm_varpar, only : ndecomp_pools_max
    use SoilBiogeochemDecompCascadeConType, only : no_soil_decomp, century_decomp, mimics_decomp, decomp_method
    !
    ! !ARGUMENTS:
    type(file_desc_t),intent(inout) :: ncid   ! pio netCDF file id
    !
    ! !LOCAL VARIABLES:
    character(len=32)  :: subname = 'readParams_netcdf'
    character(len=100) :: errCode = 'Error reading tillage params '
    logical            :: readv   ! has variable been read in or not
    real(r8), allocatable :: tempr(:,:,:)   ! temporary to read in constant
    character(len=100) :: tString ! temp. var for reading
    character(len=3)   :: decomp_method_str

    ! Initialize tillage multipliers as all 1, and exit if not tilling
    allocate(tillage_mults_allphases(ndecomp_pools, ntill_stages_max))
    tillage_mults_allphases(:,:) = 1.0_r8
    if (.not. get_do_tillage()) then
        return
    end if

    ! Handle decomposition method
    select case( decomp_method )
    case( no_soil_decomp ) 
       return
    case( century_decomp ) 
        tString = 'bgc_till_decompk_multipliers'
    case( mimics_decomp )
        tString = 'mimics_till_decompk_multipliers'
    case default
       write(decomp_method_str, '(I3)') decomp_method
       call endrun('Bad decomp_method = '//decomp_method_str )
    end select

    ! Read off of netcdf file
    allocate(tempr(2,ndecomp_pools_max,ntill_stages_max))
    call ncd_io(trim(tString), tempr, 'read', ncid, readvar = readv)
    if (.not. readv) then
        call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    end if

    ! Save
    tillage_mults_allphases = tempr(tillage_intensity,1:ndecomp_pools,:)

  end subroutine readParams_netcdf


  subroutine readParams(ncid)
    ! !USES
    use ncdio_pio , only : file_desc_t
    !
    ! !ARGUMENTS:
    type(file_desc_t),intent(inout) :: ncid   ! pio netCDF file id

    call readParams_namelist()
    call readParams_netcdf(ncid)

  end subroutine readParams


  function get_do_tillage()
    logical :: get_do_tillage
    get_do_tillage = tillage_intensity > tillage_off
  end function get_do_tillage


  subroutine get_tillage_multipliers(tillage_mults, idop)
    ! !DESCRIPTION:
    !
    !  Get the tillage effective multiplier if prognostic crops are on and
    !  tillage is turned on. Created by Sam Levis. Modified by Michael Graham
    !  to use days past planting (idpp).
    !
    !  Modified by Sam Rabin to fix a bug where idpp wasn't actually used, which
    !  would affect growing seasons that crossed over into a new calendar year.
    !  Also avoids day=1 in last timestep of year by using DaysPastPlanting(), which
    !  uses get_prev_calday() instead of get_curr_calday().
    !  Previous behavior can be requested with namelist variable use_original_tillage.
    !
    !  Original code had two versions depending on cell's GDP, but this seems to
    !  have been only an initial effort that was (a) never published and (b) not
    !  completely developed.
    !
    ! !USES:
    use clm_time_manager, only : get_curr_calday, get_curr_days_per_year
    use pftconMod       , only : ntmp_corn, nirrig_tmp_corn, ntmp_soybean, nirrig_tmp_soybean
    use CNPhenologyMod  , only : DaysPastPlanting
    ! !ARGUMENTS:
    real(r8)         , intent(inout) :: tillage_mults(:) ! tillage multipliers for this patch
    integer          , intent(in) :: idop    ! patch day of planting
    !
    ! !LOCAL VARIABLES:
    !
    integer :: day                  ! julian day
    integer :: idpp                 ! days past planting
    integer :: phase                ! which tillage phase are we in?
    real(r8) dayspyr                ! days per year
    !-----------------------------------------------------------------------
        
    ! info from DAYCENT (Melannie Hartman CSU)
    ! temp. cereals: P 30 d bef, C 15 d bef, D on day of planting
    ! corn, soy    : P           C           D           & HW-7 30 d aftr

    phase = 0

    if (use_original_tillage) then
        day = get_curr_calday()
        if (day >= idop .and. day < idop+15) then ! based on Point Chisel Tandem Disk multipliers
            phase = 1
        else if (day >= idop+15 .and. day < idop+45) then ! based on Field and Row Cultivator multipliers
            phase = 2
        else if (day >= idop+45 .and. day < idop+75) then ! based on Rod Weed Row Planter
            phase = 3
        end if
    else
        idpp = DaysPastPlanting(idop)
        if (idpp < 15) then ! based on Point Chisel Tandem Disk multipliers
            phase = 1
        else if (idpp < 45) then ! based on Field and Row Cultivator multipliers
            phase = 2
        else if (idpp < 75) then ! based on Rod Weed Row Planter
            phase = 3
        end if
    end if

    tillage_mults(:) = 1._r8
    if (phase > 0) then
        if (phase > ntill_stages_max) then
            call endrun(msg='Tillage phase > ntill_stages_max')
        end if
        tillage_mults = tillage_mults_allphases(:, phase)
    end if
    
  end subroutine get_tillage_multipliers


  subroutine get_apply_tillage_multipliers(idop, c, j, decomp_k)
    ! !DESCRIPTION:
    !
    ! Multiply decomposition rate constants by tillage coefficients.
    ! Written by Sam Rabin, based on original code by Michael Graham.
    !
    ! !USES
    use pftconMod , only : npcropmin
    use PatchType , only : patch
    !
    ! !ARGUMENTS:
    integer       , intent(in) :: idop(:) ! patch day of planting
    integer       , intent(in) :: c       ! index of column this is being called for
    integer       , intent(in) :: j       ! index of soil layer this is being called for
    real(r8), dimension(:), intent(inout) :: decomp_k ! Output: [real(r8) (:) ]  rate constant for decomposition (1./sec)
    !
    ! !LOCAL VARIABLES
    integer :: p, this_patch, n_noncrop
    real    :: sumwt ! sum of all patch weights, to check
    real(r8), dimension(ndecomp_pools) :: tillage_mults
    real(r8), dimension(ndecomp_pools) :: tillage_mults_1patch

    if (.not. col%active(c) .or. j > 5) then
        ! Top 5 layers (instead of all nlevdecomp) so that model only tills
        ! the top 26-40 cm of the soil surface, rather than whole soil - MWGraham
        return
    end if
    
    ! Initialize tillage multipliers to 0. We will loop through all patches in column,
    ! adding patch-weighted multipliers to this.
    tillage_mults(:) = 0.0_r8

    ! TODO: Figure out why adding ".and. col%lun_itype(c) == istcrop" to conditional
    !       controlling call of this subroutine didn't properly exclude non-crop columns.
    !       That working would allow some simplification here.
    this_patch = 0
    n_noncrop = 0
    sumwt = 0.0_r8
    do p = col%patchi(c),col%patchf(c)
        if (patch%active(p) .and. patch%wtcol(p) /= 0._r8) then
            if (patch%itype(p) >= npcropmin) then
                this_patch = p
                call get_tillage_multipliers(tillage_mults_1patch, idop(p))
                tillage_mults = tillage_mults + tillage_mults_1patch * patch%wtcol(p)
                sumwt = sumwt + patch%wtcol(p)
            else
                n_noncrop = n_noncrop + 1
            end if
        end if
    end do
    if (n_noncrop > 0) then
        if (this_patch > 0) then
            call endrun('ERROR Active, non-zero-weight crop AND non-crop patches found')
        end if
        return
    elseif (this_patch == 0) then
        call endrun('ERROR No active, non-zero-weight patches found (crop OR non-crop)')
    elseif (abs(1.0_r8 - sumwt) > 1.e-6_r8) then
        call endrun('ERROR Active crop patch weights does not sum to 1')
    end if

    ! Apply
    decomp_k = decomp_k * tillage_mults(:)

  end subroutine get_apply_tillage_multipliers

end module TillageMod
