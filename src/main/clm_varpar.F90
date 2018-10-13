module clm_varpar

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module containing CLM parameters
  !
  ! !USES:
  use shr_kind_mod , only: r8 => shr_kind_r8
  use spmdMod      , only: masterproc
  use clm_varctl   , only: use_extralakelayers, use_vertsoilc
  use clm_varctl   , only: use_century_decomp, use_c13, use_c14
  use clm_varctl   , only: iulog, use_crop, create_crop_landunit, irrigate
  use clm_varctl   , only: use_vichydro, soil_layerstruct
  use clm_varctl   , only: use_fates

  !
  ! !PUBLIC TYPES:
  implicit none
  save
  private

  ! Note - model resolution is read in from the surface dataset

  integer, public, parameter :: nlev_equalspace   = 15
  integer, public, parameter :: toplev_equalspace =  6
  integer, public    :: nlevsoi               ! number of hydrologically active soil layers
  integer, public    :: nlevsoifl             ! number of soil layers on input file
  integer, public    :: nlevgrnd              ! number of ground layers 
                                              ! (includes lower layers that are hydrologically inactive)
  integer, public    :: nlevurb               ! number of urban layers
  integer, public    :: nlevlak               ! number of lake layers
  integer, public    :: nlevdecomp            ! number of biogeochemically active soil layers
  integer, public    :: nlevdecomp_full       ! number of biogeochemical layers 
                                              ! (includes lower layers that are biogeochemically inactive)
  integer, public    :: nlevsno     =  -1     ! maximum number of snow layers
  integer, public, parameter :: ngases      =   3     ! CH4, O2, & CO2
  integer, public, parameter :: nlevcan     =   1     ! number of leaf layers in canopy layer
  integer, public, parameter :: nvegwcs     =   4     ! number of vegetation water conductance segments
  !ED variables
  integer, public, parameter :: numwat      =   5     ! number of water types (soil, ice, 2 lakes, wetland)
  integer, public, parameter :: numrad      =   2     ! number of solar radiation bands: vis, nir
  integer, public, parameter :: ivis        =   1     ! index for visible band
  integer, public, parameter :: inir        =   2     ! index for near-infrared band
  integer, public, parameter :: numsolar    =   2     ! number of solar type bands: direct, diffuse
  integer, public, parameter :: ndst        =   4     ! number of dust size classes (BGC only)
  integer, public, parameter :: dst_src_nbr =   3     ! number of size distns in src soil (BGC only)
  integer, public, parameter :: sz_nbr      = 200     ! number of sub-grid bins in large bin of dust size distribution (BGC only)
  integer, public, parameter :: mxpft       =  78     ! maximum number of PFT's for any mode;
  ! FIX(RF,032414) might we set some of these automatically from reading pft-physiology?
  integer, public, parameter :: nlayer      =   3     ! number of VIC soil layer --Added by AWang
  integer, public    :: nlayert               ! number of VIC soil layer + 3 lower thermal layers
  integer, public, parameter :: nvariants   =   2     ! number of variants of PFT constants

  integer, public :: maxveg           ! # of pfts + cfts
  integer, public :: maxpatch_urb= 5       ! max number of urban patches (columns) in urban landunit

  integer, public :: maxpatch_pft     ! obsolete: max number of plant functional types in naturally vegetated landunit (namelist setting)
  integer, public :: maxsoil_patches  ! # of pfts + cfts + bare ground; replaces maxpatch_pft, which is obsolete

  ! constants for decomposition cascade

  integer, public, parameter :: i_met_lit  = 1
  integer, public, parameter :: i_cel_lit  = i_met_lit + 1
  integer, public, parameter :: i_lig_lit  = i_cel_lit + 1
  integer, public    :: i_cwd

  integer, public :: ndecomp_pools
  integer, public :: ndecomp_cascade_transitions

  ! Indices used in surface file read and set in clm_varpar_init

  integer, public :: natpft_lb          ! In PATCH arrays, lower bound of Patches on the natural veg landunit (i.e., bare ground index)
  integer, public :: natpft_ub          ! In PATCH arrays, upper bound of Patches on the natural veg landunit
  integer, public :: natpft_size        ! Number of Patches on natural veg landunit (including bare ground)

  ! The following variables pertain to arrays of all PFTs - e.g., those dimensioned (g,
  ! pft_index). These include unused CFTs that are merged into other CFTs. Thus, these
  ! variables do NOT give the actual number of CFTs on the crop landunit - that number
  ! will generally be less because CLM does not simulate all crop types (some crop types
  ! are merged into other types).
  integer, public :: cft_lb             ! In arrays of PFTs, lower bound of PFTs on the crop landunit
  integer, public :: cft_ub             ! In arrays of PFTs, upper bound of PFTs on the crop landunit
  integer, public :: cft_size           ! Number of PFTs on crop landunit in arrays of PFTs

  integer, public :: maxpatch_glcmec    ! max number of elevation classes
  integer, public :: max_patch_per_col
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public clm_varpar_init          ! set parameters
  !
  !-----------------------------------------------------------------------

contains

  !------------------------------------------------------------------------------
  subroutine clm_varpar_init(actual_maxsoil_patches, actual_numcft)
    !
    ! !DESCRIPTION:
    ! Initialize module variables 
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: actual_maxsoil_patches  ! value from surface dataset
    integer, intent(in) :: actual_numcft  ! Actual number of crops
    !
    ! !LOCAL VARIABLES:
    !
    character(len=32) :: subname = 'clm_varpar_init'  ! subroutine name
    !------------------------------------------------------------------------------

    ! actual_maxsoil_patches and actual_numcft were read directly from the
    ! surface dataset
    maxsoil_patches = actual_maxsoil_patches  ! # of patches with bare ground
    maxveg = maxsoil_patches - 1  ! # of patches without bare ground

    ! For arrays containing all Patches (natural veg & crop), determine lower and upper bounds
    ! for (1) Patches on the natural vegetation landunit (includes bare ground, and includes
    ! crops if create_crop_landunit=false), and (2) CFTs on the crop landunit (no elements
    ! if create_crop_landunit=false)

    if (create_crop_landunit) then
       natpft_size = maxsoil_patches - actual_numcft  ! includes bare ground
       cft_size    = actual_numcft
    else
       natpft_size = maxsoil_patches  ! includes bare ground
       cft_size    = 0
    end if

    natpft_lb = 0
    natpft_ub = natpft_lb + natpft_size - 1
    cft_lb = natpft_ub + 1
    cft_ub = cft_lb + cft_size - 1

    ! TODO(wjs, 2015-10-04, bugz 2227) Using actual_numcft in this 'max' gives a significant
    ! overestimate of max_patch_per_col when use_crop is true. This should be reworked -
    ! or, better, removed from the code entirely (because it is a maintenance problem, and
    ! I can't imagine that looping idioms that use it help performance that much, and
    ! likely they hurt performance.)
    max_patch_per_col= max(maxsoil_patches, actual_numcft, maxpatch_urb)

    nlevsoifl   =  10
    nlevurb     =  5
    if ( masterproc ) write(iulog, *) 'soil_layerstruct varpar ',soil_layerstruct
    if ( soil_layerstruct == '10SL_3.5m' ) then
       nlevsoi     =  nlevsoifl
       nlevgrnd    =  15
    else if ( soil_layerstruct == '23SL_3.5m' ) then 
       nlevsoi     =  8  + nlev_equalspace
       nlevgrnd    =  15 + nlev_equalspace
    else if ( soil_layerstruct == '49SL_10m' ) then
      nlevsoi     =  49 ! 10x10 + 9x100 + 30x300 = 1e4mm = 10m
!       nlevsoi     =  29 ! 10x10 + 9x100 + 10x300 = 4e3mm = 4m
       nlevgrnd    =  nlevsoi+5
    else if ( soil_layerstruct == '20SL_8.5m' ) then
      nlevsoi     =  20 
      nlevgrnd    =  nlevsoi+5
    endif
    if ( masterproc ) write(iulog, *) 'soil_layerstruct varpar ',soil_layerstruct,nlevsoi,nlevgrnd

    if (use_vichydro) then
       nlayert     =  nlayer + (nlevgrnd -nlevsoi)
    endif

    ! here is a switch to set the number of soil levels for the biogeochemistry calculations.
    ! currently it works on either a single level or on nlevsoi and nlevgrnd levels
    if (use_vertsoilc) then
       nlevdecomp      = nlevsoi
       nlevdecomp_full = nlevgrnd
    else
       nlevdecomp      = 1
       nlevdecomp_full = 1
    end if

    if (.not. use_extralakelayers) then
       nlevlak     =  10     ! number of lake layers
    else
       nlevlak     =  25     ! number of lake layers (Yields better results for site simulations)
    end if

    if ( masterproc )then
       write(iulog, *) 'CLM varpar subsurface discretization levels '
       write(iulog, '(a, i3)') '    nlevsoi = ', nlevsoi
       write(iulog, '(a, i3)') '    nlevgrnd = ', nlevgrnd
       write(iulog, '(a, i3)') '    nlevdecomp = ', nlevdecomp
       write(iulog, '(a, i3)') '    nlevdecomp_full = ', nlevdecomp_full
       write(iulog, '(a, i3)') '    nlevlak = ', nlevlak
       write(iulog, *)
    end if

    if ( use_fates ) then
       i_cwd = 0
       if (use_century_decomp) then
          ndecomp_pools = 6
          ndecomp_cascade_transitions = 8
       else
          ndecomp_pools = 7
          ndecomp_cascade_transitions = 7
       end if
    else
       i_cwd = 4
       if (use_century_decomp) then
          ndecomp_pools = 7
          ndecomp_cascade_transitions = 10
       else
          ndecomp_pools = 8
          ndecomp_cascade_transitions = 9
       end if
    endif

  end subroutine clm_varpar_init

end module clm_varpar
