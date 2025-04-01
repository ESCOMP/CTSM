module clm_varpar

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module containing CLM parameters
  !
  ! !USES:
  use shr_kind_mod , only: r8 => shr_kind_r8
  use shr_sys_mod  , only: shr_sys_abort
  use spmdMod      , only: masterproc
  use clm_varctl   , only: use_extralakelayers
  use clm_varctl   , only: use_c13, use_c14
  use clm_varctl   , only: iulog, use_crop, create_crop_landunit, irrigate
  use clm_varctl   , only: use_vichydro, rundef
  use clm_varctl   , only: soil_layerstruct_predefined
  use clm_varctl   , only: soil_layerstruct_userdefined
  use clm_varctl   , only: soil_layerstruct_userdefined_nlevsoi
  use clm_varctl   , only: use_fates, use_cn, use_fates_sp

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
  integer, public    :: nlevmaxurbgrnd        ! maximum of the number of ground and urban layers
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
  integer, public, parameter :: mxsowings   =   1     ! maximum number of crop growing seasons to begin in any year;
  integer, public            :: mxharvests            ! maximum number of crop harvests in any year
                                                      ! (allows for multiple harvests in a calendar year in case harvest occurs near
                                                      ! beginning/end of year);
  ! FIX(RF,032414) might we set some of these automatically from reading pft-physiology?
  integer, public, parameter :: nlayer      =   3     ! number of VIC soil layer --Added by AWang
  integer, public    :: nlayert               ! number of VIC soil layer + 3 lower thermal layers
  integer, public, parameter :: nvariants   =   2     ! number of variants of PFT constants

  ! CN Matrix solution sizes
  integer, public, parameter :: nvegpool_natveg = 18  ! number of vegetation matrix pool without crop
  integer, public, parameter :: nvegpool_crop   =  3  ! number of vegetation matrix pool with crop
  integer, public, parameter :: nveg_retransn   =  1  ! number of vegetation retranslocation pool
  integer, public :: nvegcpool                        ! number of vegetation C pools
  integer, public :: nvegnpool                        ! number of vegetation N pools

  integer, public :: maxveg                           ! # of pfts + cfts
  integer, public :: maxpatch_urb= 5                  ! max number of urban patches (columns) in urban landunit

  integer, public :: maxsoil_patches                  ! # of pfts + cfts + bare ground; replaces maxpatch_pft, which is obsolete

  ! constants for decomposition cascade

  integer, public, parameter :: i_litr1 = 1   ! TEMPORARY FOR CascadeCN TO BUILD
  integer, public            :: i_litr2 = -9  ! TEMPORARY FOR CascadeCN TO BUILD
  integer, public            :: i_litr3 = -9  ! TEMPORARY FOR CascadeCN TO BUILD
  ! The code currently expects i_litr_min = i_met_lit = 1 and
  !                            i_litr_max = 2 or 3
  integer, public            :: i_litr_min    = -9    ! min index of litter pools; overwritten in SoilBiogeochemDecompCascade*Mod
  integer, public            :: i_litr_max    = -9    ! max index of litter pools; overwritten in SoilBiogeochemDecompCascade*Mod
  integer, public            :: i_met_lit     = -9    ! index of metabolic litter pool; overwritten in SoilBiogeochemDecompCascade*Mod
  integer, public            :: i_str_lit     = -9    ! index of structural litter pool; overwritten in SoilBiogeochemDecompCascade*Mod
  integer, public            :: i_phys_som    = -9    ! index of physically protected Soil Organic Matter (SOM); overwritten in SoilBiogeochemDecompCascade*Mod
  integer, public            :: i_chem_som    = -9    ! index of chemically protected Soil Organic Matter (SOM); overwritten in SoilBiogeochemDecompCascade*Mod
  integer, public            :: i_cop_mic     = -9    ! index of copiotrophic microbial pool; overwritten in SoilBiogeochemDecompCascade*Mod
  integer, public            :: i_oli_mic     = -9    ! index of oligotrophic microbial pool; overwritten in SoilBiogeochemDecompCascade*Mod
  integer, public            :: i_cwd         = -9    ! index of cwd pool; overwritten in SoilBiogeochemDecompCascade*Mod
  integer, public            :: i_cwdl2       = -9    ! index of cwd to l2 transition; overwritten in SoilBiogeochemDecompCascade*Mod
  integer, public, parameter :: ileaf         = 1     ! leaf pool index
  integer, public, parameter :: ileaf_st      = 2     ! leaf storage pool index
  integer, public, parameter :: ileaf_xf      = 3     ! leaf transfer pool index
  integer, public, parameter :: ifroot        = 4     ! fine root pool index
  integer, public, parameter :: ifroot_st     = 5     ! fine root storage pool index
  integer, public, parameter :: ifroot_xf     = 6     ! fine root transfer pool index
  integer, public, parameter :: ilivestem     = 7     ! live stem pool index
  integer, public, parameter :: ilivestem_st  = 8     ! live stem storage pool index
  integer, public, parameter :: ilivestem_xf  = 9     ! live stem transfer pool index
  integer, public, parameter :: ideadstem     = 10    ! dead stem pool index
  integer, public, parameter :: ideadstem_st  = 11    ! dead stem storage pool index
  integer, public, parameter :: ideadstem_xf  = 12    ! dead stem transfer pool index
  integer, public, parameter :: ilivecroot    = 13    ! live coarse root pool index
  integer, public, parameter :: ilivecroot_st = 14    ! live coarse root storage pool index
  integer, public, parameter :: ilivecroot_xf = 15    ! live coarse root transfer pool index
  integer, public, parameter :: ideadcroot    = 16    ! dead coarse root pool index
  integer, public, parameter :: ideadcroot_st = 17    ! dead coarse root storage pool index
  integer, public, parameter :: ideadcroot_xf = 18    ! dead coarse root transfer pool index
  integer, public, parameter :: igrain        = 19    ! grain pool index
  integer, public, parameter :: igrain_st     = 20    ! grain storage pool index
  integer, public, parameter :: igrain_xf     = 21    ! grain transfer pool index

  integer, public :: ncphtrans      !maximum number of vegetation C transfers through phenology
  integer, public :: ncphouttrans   !maximum number of vegetation C transfers out of vegetation through phenology
  integer, public :: ncgmtrans      !maximum number of vegetation C transfers through gap mortality
  integer, public :: ncgmouttrans   !maximum number of vegetation C transfers out of vegetation through gap mortality
  integer, public :: ncfitrans      !maximum number of vegetation C transfers through fire
  integer, public :: ncfiouttrans   !maximum number of vegetation C transfers out of vegetation trhough fire
  integer, public :: nnphtrans      !maximum number of vegetation N transfers through phenology
  integer, public :: nnphouttrans   !maximum number of vegetation N transfers out of vegetation through phenology
  integer, public :: nngmtrans      !maximum number of vegetation N transfers through gap mortality
  integer, public :: nngmouttrans   !maximum number of vegetation N transfers out of vegetation through gap mortality
  integer, public :: nnfitrans      !maximum number of vegetation N transfers through fire
  integer, public :: nnfiouttrans   !maximum number of vegetation N transfers out of vegetation trhough fire

  integer, public :: iretransn      ! retranslocation pool index

  integer, public :: ioutc          ! external C pool index
  integer, public :: ioutn          ! external N pool index

  integer, public :: ndecomp_pools_max
  integer, public :: ndecomp_pools  ! total number of pools
  integer, public :: ndecomp_cascade_transitions
  integer, public :: ndecomp_cascade_outtransitions

  ! for soil matrix
  integer, public  :: ndecomp_pools_vr  ! ndecomp_pools * levels in the vertical

  ! Indices used in surface file read and set in clm_varpar_init

  integer, public :: natpft_lb          ! In PATCH arrays, lower bound of Patches on the natural veg landunit (i.e., bare ground index)
  integer, public :: natpft_ub          ! In PATCH arrays, upper bound of Patches on the natural veg landunit
  integer, public :: natpft_size        ! Number of Patches on natural veg landunit (including bare ground)

  integer, public :: surfpft_lb         ! Lower bound of PFTs in the surface file
                                        ! synonymous with natpft_lb for non-fates and fates-sp
  integer, public :: surfpft_ub         ! Upper bound of PFTs in the surface file
                                        ! synonymous with natpft_ub for non-fates and fates-sp

  
  ! The following variables pertain to arrays of all PFTs - e.g., those dimensioned (g,
  ! pft_index). These include unused CFTs that are merged into other CFTs. Thus, these
  ! variables do NOT give the actual number of CFTs on the crop landunit - that number
  ! will generally be less because CLM does not simulate all crop types (some crop types
  ! are merged into other types).
  integer, public :: cft_lb             ! In arrays of PFTs, lower bound of PFTs on the crop landunit
  integer, public :: cft_ub             ! In arrays of PFTs, upper bound of PFTs on the crop landunit
  integer, public :: cft_size           ! Number of PFTs on crop landunit in arrays of PFTs

  integer, public :: maxpatch_glc    ! max number of elevation classes
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public clm_varpar_init          ! set parameters
  !
  !-----------------------------------------------------------------------

contains

  !------------------------------------------------------------------------------
  subroutine clm_varpar_init(actual_maxsoil_patches, surf_numpft, surf_numcft, actual_nlevurb)
    !
    ! !DESCRIPTION:
    ! Initialize module variables 
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: actual_maxsoil_patches  ! Number of soil patches to allocate
                                                   ! This value comes either from the
                                                   ! surface dataset (non-fates) or 
                                                   ! from fates (via its parameter file)
    integer, intent(in) :: surf_numpft             ! Number of PFTs in the surf dataset
    integer, intent(in) :: surf_numcft             ! Number of CFTs in the surf dataset
    integer, intent(in) :: actual_nlevurb          ! nlevurb from surface dataset
    !
    ! !LOCAL VARIABLES:
    !
    integer :: j  ! loop index
    character(len=32) :: subname = 'clm_varpar_init'  ! subroutine name
    !------------------------------------------------------------------------------

    ! actual_maxsoil_patches is either the total number of cfts+pfts in the surface
    ! file (for non-fates), or the number of patches plus the bareground that
    ! fates requests.  If this is a fates-sp run, this value will also be the number 
    ! of cfts+pfts as in a nonfates run

    maxsoil_patches = actual_maxsoil_patches
    
    maxveg = maxsoil_patches - 1  ! # of patches without bare ground

    ! For arrays containing all Patches (natural veg & crop), determine lower and upper bounds
    ! for (1) Patches on the natural vegetation landunit (includes bare ground, and includes
    ! crops if create_crop_landunit=false), and (2) CFTs on the crop landunit (no elements
    ! if create_crop_landunit=false)
    ! As for when we don't have a crop LU, which is currently when FATES is on...
    ! These values are used to create the wt_nat_patch array that is used by fates_sp 
    ! and fixed biogeog. Also, the pft and cft vectors are concatenated into
    ! the natpft vector (wt_nat_patch), the wt_cft array is unused (size zero)
    ! The following values should not be used for allocating patch structures
    ! though.  That should be handled completely by maxoil_patches and maxveg

    if (create_crop_landunit) then
       
       natpft_size = surf_numpft    ! includes bare ground + natveg pfts
       cft_size    = surf_numcft
       natpft_lb   = 0
       natpft_ub   = natpft_lb + natpft_size - 1
       cft_lb      = natpft_ub + 1
       cft_ub      = cft_lb + cft_size - 1
       surfpft_lb  = natpft_lb
       surfpft_ub  = natpft_ub
       
    else ! only true when FATES is active
       
       natpft_size = maxsoil_patches
       cft_size    = 0
       natpft_lb   = 0
       natpft_ub   = natpft_lb + natpft_size - 1
       cft_lb      = 0
       cft_ub      = 0
       surfpft_lb  = 0
       surfpft_ub  = surf_numpft+surf_numcft-1
       
    end if
       
    if(use_fates_sp .and. (natpft_ub .ne. maxveg) ) then
       write(iulog,*) 'maxveg should match the upper bound for non-fates and fates-sp runs'
       write(iulog,*) 'the surface dataset PFT+CFT indices (ie lsmft), yours: ',natpft_ub,maxveg
       call shr_sys_abort(subname//' ERROR: conflict in maxveg and pft bounds')
    end if
    
    mxharvests = mxsowings + 1

    nlevsoifl   =  10
    nlevurb     =  actual_nlevurb

    if ( masterproc ) write(iulog, *) 'soil_layerstruct_predefined varpar ', soil_layerstruct_predefined
    if ( masterproc ) write(iulog, *) 'soil_layerstruct_userdefined varpar ', soil_layerstruct_userdefined

    if (soil_layerstruct_userdefined(1) /= rundef) then  ! user defined soil layers
       if (soil_layerstruct_predefined /= 'UNSET') then
          write(iulog,*) subname//' ERROR: Both soil_layerstruct_predefined and soil_layer_userdefined have values'
          call shr_sys_abort(subname//' ERROR: Cannot decide how to set the soil layer structure')
       else
          nlevgrnd = size(soil_layerstruct_userdefined)
          ! loops backwards until it hits the last valid user-defined value
          do j = nlevgrnd,1,-1
             if (soil_layerstruct_userdefined(j) /= rundef) then
                exit
             else
                nlevgrnd = nlevgrnd - 1
             end if
          end do
          nlevsoi = soil_layerstruct_userdefined_nlevsoi  ! read in namelist
          if (nlevsoi >= nlevgrnd) then
             write(iulog,*) subname//' ERROR: nlevsoi >= nlevgrnd; did you enter soil_layerstruct_userdefined_nlevsoi correctly in user_nl_clm?'
             call shr_sys_abort(subname//' ERROR: nlevsoi must be less than nlevgrnd')
          end if
       end if
    else  ! pre-defined soil structure options
       if ( soil_layerstruct_predefined == '10SL_3.5m' ) then
          nlevsoi     =  nlevsoifl
          nlevgrnd    =  15
       else if ( soil_layerstruct_predefined == '23SL_3.5m' ) then
          nlevsoi     =  8  + nlev_equalspace
          nlevgrnd    =  15 + nlev_equalspace
       else if ( soil_layerstruct_predefined == '49SL_10m' ) then
          nlevsoi     =  49 ! 10x10 + 9x100 + 30x300 = 1e4mm = 10m
!          nlevsoi     =  29 ! 10x10 + 9x100 + 10x300 = 4e3mm = 4m
          nlevgrnd    =  nlevsoi+5
       else if ( soil_layerstruct_predefined == '20SL_8.5m' ) then
         nlevsoi     =  20
         nlevgrnd    =  nlevsoi+5
       else if ( soil_layerstruct_predefined == '4SL_2m' ) then
          nlevsoi     =  4
          nlevgrnd    =  5
       else if (soil_layerstruct_predefined == 'UNSET') then
          write(iulog,*) subname//' ERROR: Both soil_layerstruct_predefined and soil_layer_userdefined currently undefined'
          call shr_sys_abort(subname//' ERROR: Cannot set the soil layer structure')
       else
          write(iulog,*) subname//' ERROR: Unrecognized pre-defined soil layer structure: ', trim(soil_layerstruct_predefined)
          call shr_sys_abort(subname//' ERROR: Unrecognized pre-defined soil layer structure')
       end if
    endif
    nlevmaxurbgrnd = max0(nlevurb,nlevgrnd)
    if ( masterproc ) write(iulog, *) 'nlevsoi, nlevgrnd varpar ', nlevsoi, nlevgrnd

    if (use_vichydro) then
       nlayert     =  nlayer + (nlevgrnd -nlevsoi)
    endif

    !
    ! Number of layers for soil decomposition
    !
    if ( use_cn .or. (use_fates .and. .not. use_fates_sp) ) then
       ! to set the number of soil levels for the biogeochemistry calculations.
       ! currently it works on nlevsoi and nlevgrnd levels
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

    ! CN Matrix settings
    if (use_crop)then
       nvegcpool = nvegpool_natveg + nvegpool_crop
       ncphtrans = 18
       nnphtrans = 37
       ncphouttrans = 4
       nnphouttrans = 5
    else
       nvegcpool = nvegpool_natveg
       ncphtrans = 17
       nnphtrans = 34
       ncphouttrans = 3
       nnphouttrans = 4
    end if
    ncgmtrans = 18
    ncgmouttrans = 18
    ncfitrans = 20
    ncfiouttrans = 18
    nngmtrans = 19
    nngmouttrans = 19
    nnfitrans = 21
    nnfiouttrans = 19
    nvegnpool = nvegcpool + 1
    iretransn = nvegnpool
    ioutc = nvegcpool + 1
    ioutn = nvegnpool + 1

  end subroutine clm_varpar_init

end module clm_varpar
