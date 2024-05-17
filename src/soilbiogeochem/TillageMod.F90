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
  public :: get_fraction_tilled
  ! !PUBLIC DATA MEMBERS
  character(len=CS), public :: tillage_mode     ! off, low, high
  !
  ! !PRIVATE DATA MEMBERS
  integer             :: tillage_intensity
  integer, parameter  :: tillage_off = 0
  integer, parameter  :: tillage_low = 1
  integer, parameter  :: tillage_high = 2
  logical  :: use_original_tillage_phases ! Use buggy tillage phase determination?
  real(r8), pointer :: tillage_mults_allphases(:,:) ! (ndecomp_pools, ntill_stages_max)
  integer, parameter :: ntill_stages_max = 3 ! How many different tillage phases are there? (Not including all-1 phases.)
  integer, parameter :: ntill_intensities_max = 2 ! How many different tillage intensities are allowed (other than "off")?
  real(r8)           :: max_tillage_depth ! Maximum depth to till (m)

!==============================================================================
contains
!==============================================================================

  subroutine readParams_namelist(NLFilename)
    !
    ! Read namelist parameters related to tillage.
    !
    ! !USES:
    use spmdMod        , only : masterproc, mpicom
    use clm_nlUtilsMod , only : find_nlgroup_name
    use shr_mpi_mod    , only : shr_mpi_bcast
    ! !ARGUMENTS:
    character(len=*), intent(in) :: NLFilename ! Namelist filename
    !
    ! !LOCAL VARIABLES
    integer                :: nu_nml       ! unit for namelist file
    integer                :: nml_error    ! namelist i/o error flag
    character(*), parameter :: subname = "('readParams_namelist')"

    namelist /tillage_inparm/    &
        tillage_mode,            &
        use_original_tillage_phases,    &
        max_tillage_depth

    ! Default values
    tillage_mode = 'off'
    use_original_tillage_phases = .false.
    max_tillage_depth = 0.26_r8  ! Graham et al. (2021) unintentionally used 0.32

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
     call shr_mpi_bcast(use_original_tillage_phases , mpicom)
     call shr_mpi_bcast(max_tillage_depth, mpicom)

     if (masterproc) then
        write(iulog,*) ' '
        write(iulog,*) 'tillage settings:'
        write(iulog,*) '  tillage_mode  = ',tillage_mode
        write(iulog,*) '  use_original_tillage_phases   = ',use_original_tillage_phases
        write(iulog,*) '  max_tillage_depth = ',max_tillage_depth
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
    allocate(tempr(ntill_intensities_max,ndecomp_pools_max,ntill_stages_max))
    call ncd_io(trim(tString), tempr, 'read', ncid, readvar = readv, posNOTonfile = .true.)
    if (.not. readv) then
        call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    end if

    ! Save
    tillage_mults_allphases = tempr(tillage_intensity,1:ndecomp_pools,:)

  end subroutine readParams_netcdf


  subroutine readParams(ncid, NLFilename)
    ! !USES
    use ncdio_pio , only : file_desc_t
    !
    ! !ARGUMENTS:
    type(file_desc_t),intent(inout) :: ncid   ! pio netCDF file id
    character(len=*), intent(in) :: NLFilename ! Namelist filename

    call readParams_namelist(NLFilename)
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
    !  Previous behavior can be requested with namelist variable use_original_tillage_phases.
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

    if (use_original_tillage_phases) then
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


  function get_fraction_tilled(layer_bottom, layer_thickness, max_tillage_depth_gft) result(fraction_tilled)
    ! !ARGUMENTS
    real(r8), intent(in) :: layer_bottom    ! Soil layer interface (between j and j+1) depth (zisoi)
    real(r8), intent(in) :: layer_thickness ! Soil layer thickness (dzsoi_decomp)
    real(r8) :: max_tillage_depth_gft ! Maximum tillage depth
    ! !LOCAL VARIABLES
    real(r8) :: layer_top
    ! !RESULT
    real(r8) :: fraction_tilled ! Fraction of this layer that's within the tillage depth

    ! If the top of the layer is below the max tillage depth, do not till.
    layer_top = layer_bottom - layer_thickness
    if (layer_top > max_tillage_depth_gft) then
        fraction_tilled = 0._r8
        return
    end if

    ! Handle zero-thickness layers. This may not be necessary.
    if (layer_thickness == 0._r8) then
        if (layer_bottom <= max_tillage_depth_gft) then
            fraction_tilled = 1._r8
        else
            fraction_tilled = 0._r8
        end if
        return
    end if

    fraction_tilled = max(0._r8, min(1._r8, (max_tillage_depth_gft - layer_top) / layer_thickness))

  end function get_fraction_tilled


  subroutine get_apply_tillage_multipliers(idop, c, j, decomp_k)
    ! !DESCRIPTION:
    !
    ! Multiply decomposition rate constants by tillage coefficients.
    ! Written by Sam Rabin, based on original code by Michael Graham.
    !
    ! !USES
    use pftconMod , only : npcropmin
    use clm_varcon, only : zisoi, dzsoi_decomp
    use landunit_varcon , only : istcrop
    use PatchType , only : patch
    !
    ! !ARGUMENTS:
    integer       , intent(in) :: idop(:) ! patch day of planting
    integer       , intent(in) :: c       ! index of column this is being called for
    integer       , intent(in) :: j       ! index of soil layer this is being called for
    real(r8), dimension(:), intent(inout) :: decomp_k ! Output: [real(r8) (:) ]  rate constant for decomposition (1./sec)
    !
    ! !LOCAL VARIABLES
    integer :: p
    real    :: sumwt ! sum of all patch weights, to check
    real(r8), dimension(ndecomp_pools) :: tillage_mults
    real(r8), dimension(ndecomp_pools) :: tillage_mults_1patch
    real(r8) :: fraction_tilled ! Fraction of this layer that's within the tillage depth

    ! Skip tillage if column is inactive or this layer doesn't get tilled
    fraction_tilled = get_fraction_tilled(zisoi(j), dzsoi_decomp(j), max_tillage_depth)
    if (.not. col%active(c) .or. fraction_tilled == 0._r8 .or. col%lun_itype(c) /= istcrop) then
        return
    end if
    
    ! Initialize tillage multipliers to 0. We will loop through all patches in column,
    ! adding patch-weighted multipliers to this.
    tillage_mults(:) = 0.0_r8

    sumwt = 0.0_r8
    do p = col%patchi(c),col%patchf(c)
        if (patch%active(p) .and. patch%wtcol(p) /= 0._r8) then
            if (patch%itype(p) < npcropmin) then
                ! Do not till generic crops
                tillage_mults_1patch(:) = 1._r8
            else
                call get_tillage_multipliers(tillage_mults_1patch, idop(p))
            end if
            tillage_mults = tillage_mults + tillage_mults_1patch * patch%wtcol(p)
            sumwt = sumwt + patch%wtcol(p)
        end if
    end do
    if (abs(1.0_r8 - sumwt) > 1.e-6_r8) then
        call endrun('ERROR Active crop patch weights does not sum to 1')
    end if

    ! Adjust tillage_mults to consider fraction of this layer that's within tillage depth.
    tillage_mults = tillage_mults *          fraction_tilled &
                    + 1._r8       * (1._r8 - fraction_tilled)

    ! Apply
    decomp_k = decomp_k * tillage_mults(:)

  end subroutine get_apply_tillage_multipliers

end module TillageMod
