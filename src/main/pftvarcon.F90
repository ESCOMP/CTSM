module pftvarcon

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: pftvarcon
!
! !DESCRIPTION:
! Module containing vegetation constants and method to
! read and initialize vegetation (PFT) constants.
!
! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use abortutils  , only : endrun
  use clm_varpar  , only : numpft, numrad, ivis, inir
  use clm_varctl  , only : iulog
!
! !PUBLIC TYPES:
  implicit none
  save
!
! Vegetation type constants
!
  character(len=40) pftname(0:numpft) !PFT description

  integer :: noveg                  !value for not vegetated 
  integer :: ndllf_evr_tmp_tree     !value for Needleleaf evergreen temperate tree
  integer :: ndllf_evr_brl_tree     !value for Needleleaf evergreen boreal tree
  integer :: ndllf_dcd_brl_tree     !value for Needleleaf deciduous boreal tree
  integer :: nbrdlf_evr_trp_tree    !value for Broadleaf evergreen tropical tree
  integer :: nbrdlf_evr_tmp_tree    !value for Broadleaf evergreen temperate tree
  integer :: nbrdlf_dcd_trp_tree    !value for Broadleaf deciduous tropical tree
  integer :: nbrdlf_dcd_tmp_tree    !value for Broadleaf deciduous temperate tree
  integer :: nbrdlf_dcd_brl_tree    !value for Broadleaf deciduous boreal tree
  integer :: ntree                  !value for last type of tree
  integer :: nbrdlf_evr_shrub       !value for Broadleaf evergreen shrub
  integer :: nbrdlf_dcd_tmp_shrub   !value for Broadleaf deciduous temperate shrub
  integer :: nbrdlf_dcd_brl_shrub   !value for Broadleaf deciduous boreal shrub
  integer :: nc3_arctic_grass       !value for C3 arctic grass
  integer :: nc3_nonarctic_grass    !value for C3 non-arctic grass
  integer :: nc4_grass              !value for C4 grass
  integer :: nc3crop                !value for generic crop
  integer :: nirrig                 !value for irrigated generic crop

  real(r8):: dleaf(0:numpft)       !characteristic leaf dimension (m)
  real(r8):: c3psn(0:numpft)       !photosynthetic pathway: 0. = c4, 1. = c3
  real(r8):: vcmx25(0:numpft)      !max rate of carboxylation at 25C (umol CO2/m**2/s)
  real(r8):: mp(0:numpft)          !slope of conductance-to-photosynthesis relationship
  real(r8):: qe25(0:numpft)        !quantum efficiency at 25C (umol CO2 / umol photon)
  real(r8):: xl(0:numpft)          !leaf/stem orientation index
  real(r8):: rhol(0:numpft,numrad) !leaf reflectance: 1=vis, 2=nir
  real(r8):: rhos(0:numpft,numrad) !stem reflectance: 1=vis, 2=nir
  real(r8):: taul(0:numpft,numrad) !leaf transmittance: 1=vis, 2=nir
  real(r8):: taus(0:numpft,numrad) !stem transmittance: 1=vis, 2=nir
  real(r8):: z0mr(0:numpft)        !ratio of momentum roughness length to canopy top height (-)
  real(r8):: displar(0:numpft)     !ratio of displacement height to canopy top height (-)
  real(r8):: roota_par(0:numpft)   !CLM rooting distribution parameter [1/m]
  real(r8):: rootb_par(0:numpft)   !CLM rooting distribution parameter [1/m]
  real(r8):: crop(0:numpft)        ! crop pft: 0. = not crop, 1. = crop pft
  real(r8):: irrigated(0:numpft)   ! irrigated pft: 0. = not, 1. = irrigated
  real(r8):: smpso(0:numpft)       !soil water potential at full stomatal opening (mm)
  real(r8):: smpsc(0:numpft)       !soil water potential at full stomatal closure (mm)
  real(r8):: fnitr(0:numpft)       !foliage nitrogen limitation factor (-)
  ! begin new pft parameters for CN code
  real(r8):: slatop(0:numpft)      !SLA at top of canopy [m^2/gC]
  real(r8):: dsladlai(0:numpft)    !dSLA/dLAI [m^2/gC]
  real(r8):: leafcn(0:numpft)      !leaf C:N [gC/gN]
  real(r8):: flnr(0:numpft)        !fraction of leaf N in Rubisco [no units]
  real(r8):: woody(0:numpft)       !woody lifeform flag (0 or 1)
  real(r8):: lflitcn(0:numpft)      !leaf litter C:N (gC/gN)
  real(r8):: frootcn(0:numpft)      !fine root C:N (gC/gN)
  real(r8):: livewdcn(0:numpft)     !live wood (phloem and ray parenchyma) C:N (gC/gN)
  real(r8):: deadwdcn(0:numpft)     !dead wood (xylem and heartwood) C:N (gC/gN)

  real(r8):: froot_leaf(0:numpft)   !allocation parameter: new fine root C per new leaf C (gC/gC) 
  real(r8):: stem_leaf(0:numpft)    !allocation parameter: new stem c per new leaf C (gC/gC)
  real(r8):: croot_stem(0:numpft)   !allocation parameter: new coarse root C per new stem C (gC/gC)
  real(r8):: flivewd(0:numpft)      !allocation parameter: fraction of new wood that is live (phloem and ray parenchyma) (no units)
  real(r8):: fcur(0:numpft)         !allocation parameter: fraction of allocation that goes to currently displayed growth, remainder to storage
  real(r8):: fcurdv(0:numpft)       !alternate fcur for use with cndv
  real(r8):: lf_flab(0:numpft)      !leaf litter labile fraction
  real(r8):: lf_fcel(0:numpft)      !leaf litter cellulose fraction
  real(r8):: lf_flig(0:numpft)      !leaf litter lignin fraction
  real(r8):: fr_flab(0:numpft)      !fine root litter labile fraction
  real(r8):: fr_fcel(0:numpft)      !fine root litter cellulose fraction
  real(r8):: fr_flig(0:numpft)      !fine root litter lignin fraction
  real(r8):: dw_fcel(0:numpft)      !dead wood cellulose fraction
  real(r8):: dw_flig(0:numpft)      !dead wood lignin fraction
  real(r8):: leaf_long(0:numpft)    !leaf longevity (yrs)
  real(r8):: evergreen(0:numpft)    !binary flag for evergreen leaf habit (0 or 1)
  real(r8):: stress_decid(0:numpft) !binary flag for stress-deciduous leaf habit (0 or 1)
  real(r8):: season_decid(0:numpft) !binary flag for seasonal-deciduous leaf habit (0 or 1)

  ! new pft parameters for CN-fire code
  real(r8):: resist(0:numpft)       !resistance to fire (no units)

  ! pft parameters for CNDV code
  ! from LPJ subroutine pftparameters
  real(r8) pftpar20(0:numpft)       !tree maximum crown area (m2)
  real(r8) pftpar28(0:numpft)       !min coldest monthly mean temperature
  real(r8) pftpar29(0:numpft)       !max coldest monthly mean temperature
  real(r8) pftpar30(0:numpft)       !min growing degree days (>= 5 deg C)
  real(r8) pftpar31(0:numpft)       !upper limit of temperature of the warmest month (twmax)
  real(r8), parameter :: reinickerp = 1.6_r8 !parameter in allometric equation
  real(r8), parameter :: dwood  = 2.5e5_r8   !cn wood density (gC/m3); lpj:2.0e5
  real(r8), parameter :: allom1 = 100.0_r8   !parameters in
  real(r8), parameter :: allom2 =  40.0_r8   !...allometric
  real(r8), parameter :: allom3 =   0.5_r8   !...equations
  real(r8), parameter :: allom1s = 250.0_r8  !modified for shrubs by
  real(r8), parameter :: allom2s =   8.0_r8  !X.D.Z
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: pftconrd ! Read and initialize vegetation (PFT) constants
!
! !REVISION HISTORY:
! Created by Sam Levis (put into module form by Mariana Vertenstein)
! 10/21/03, Peter Thornton: Added new variables for CN code
! 06/24/09, Erik Kluzek: Add indices for all pft types, and add expected_pftnames array and comparision
! 09/17/10, David Lawrence: Modified code to read in netCDF pft physiology file
!
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: pftconrd
!
! !INTERFACE:
  subroutine pftconrd
!
! !DESCRIPTION:
! Read and initialize vegetation (PFT) constants
!
! !USES:
    use fileutils , only : getfil
    use ncdio_pio , only : ncd_io, ncd_pio_closefile, ncd_pio_openfile, file_desc_t, &
                           ncd_inqdid, ncd_inqdlen
    use clm_varctl, only : fpftcon
    use spmdMod   , only : masterproc, mpicom, MPI_REAL8, MPI_CHARACTER, MPI_INTEGER
    use nanMod    , only : inf
!
! !ARGUMENTS:
    implicit none
!
! !CALLED FROM:
! routine initialize in module initializeMod
!
! !REVISION HISTORY:
! Created by Gordon Bonan
!
!
! !LOCAL VARIABLES:
!EOP
    character(len=256) :: locfn ! local file name
    integer :: i,n              ! loop indices
    integer :: ier              ! error code
    type(file_desc_t) :: ncid   ! pio netCDF file id
    integer :: dimid            ! netCDF dimension id
    integer :: npft             ! number of pfts on pft-physiology file
    logical :: readv            ! read variable in or not
    character(len=32) :: subname = 'pftconrd'              ! subroutine name
    !
    ! Expected PFT names: The names expected on the fpftcon file and the order they are expected to be in.
    ! NOTE: similar types are assumed to be together, first trees (ending with broadleaf_deciduous_boreal_tree
    !       then shrubs, ending with broadleaf_deciduous_boreal_shrub, then grasses starting with c3_arctic_grass
    !       and finally crops, ending with wheat
    ! DO NOT CHANGE THE ORDER -- WITHOUT MODIFYING OTHER PARTS OF THE CODE WHERE THE ORDER MATTERS!
    !
    character(len=40), parameter :: expected_pftnames(0:numpft) = (/ &
                 'not_vegetated                      '  &
               , 'needleleaf_evergreen_temperate_tree'  &
               , 'needleleaf_evergreen_boreal_tree   '  &
               , 'needleleaf_deciduous_boreal_tree   '  &
               , 'broadleaf_evergreen_tropical_tree  '  &
               , 'broadleaf_evergreen_temperate_tree '  &
               , 'broadleaf_deciduous_tropical_tree  '  &
               , 'broadleaf_deciduous_temperate_tree '  &
               , 'broadleaf_deciduous_boreal_tree    '  &
               , 'broadleaf_evergreen_shrub          '  &
               , 'broadleaf_deciduous_temperate_shrub'  &
               , 'broadleaf_deciduous_boreal_shrub   '  &
               , 'c3_arctic_grass                    '  &
               , 'c3_non-arctic_grass                '  &
               , 'c4_grass                           '  &
               , 'c3_crop                            '  &
               , 'c3_irrigated                       '  &
    /)
!-----------------------------------------------------------------------

    ! Set specific vegetation type values

    if (masterproc) then
       write(iulog,*) 'Attempting to read PFT physiological data .....'
    end if
    call getfil (fpftcon, locfn, 0)
    call ncd_pio_openfile (ncid, trim(locfn), 0)
    call ncd_inqdid(ncid,'pft',dimid)
    call ncd_inqdlen(ncid,dimid,npft)

    call ncd_io('pftname',pftname, 'read', ncid, readvar=readv) 
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    call ncd_io('z0mr',z0mr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    call ncd_io('displar',displar, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    call ncd_io('dleaf',dleaf, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    call ncd_io('c3psn',c3psn, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    call ncd_io('vcmx25',vcmx25, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    call ncd_io('mp',mp, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    call ncd_io('qe25',qe25, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    call ncd_io('rholvis',rhol(:,ivis), 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    call ncd_io('rholnir',rhol(:,inir), 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    call ncd_io('rhosvis',rhos(:,ivis), 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    call ncd_io('rhosnir', rhos(:,inir), 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    call ncd_io('taulvis',taul(:,ivis), 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    call ncd_io('taulnir',taul(:,inir), 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    call ncd_io('tausvis',taus(:,ivis), 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    call ncd_io('tausnir',taus(:,inir), 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    call ncd_io('xl',xl, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    call ncd_io('roota_par',roota_par, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    call ncd_io('rootb_par',rootb_par, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    call ncd_io('slatop',slatop, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    call ncd_io('dsladlai',dsladlai, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    call ncd_io('leafcn',leafcn, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    call ncd_io('flnr',flnr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    call ncd_io('smpso',smpso, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    call ncd_io('smpsc',smpsc, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    call ncd_io('fnitr',fnitr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    call ncd_io('woody',woody, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    call ncd_io('lflitcn',lflitcn, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    call ncd_io('frootcn',frootcn, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    call ncd_io('livewdcn',livewdcn, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    call ncd_io('deadwdcn',deadwdcn, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    call ncd_io('froot_leaf',froot_leaf, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    call ncd_io('stem_leaf',stem_leaf, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    call ncd_io('croot_stem',croot_stem, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    call ncd_io('flivewd',flivewd, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    call ncd_io('fcur',fcur, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    call ncd_io('fcurdv',fcurdv, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    call ncd_io('lf_flab',lf_flab, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    call ncd_io('lf_fcel',lf_fcel, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    call ncd_io('lf_flig',lf_flig, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    call ncd_io('fr_flab',fr_flab, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    call ncd_io('fr_fcel',fr_fcel, 'read', ncid, readvar=readv)    
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    call ncd_io('fr_flig',fr_flig, 'read', ncid, readvar=readv)    
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    call ncd_io('dw_fcel',dw_fcel, 'read', ncid, readvar=readv)    
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    call ncd_io('dw_flig',dw_flig, 'read', ncid, readvar=readv)    
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    call ncd_io('leaf_long',leaf_long, 'read', ncid, readvar=readv)    
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    call ncd_io('evergreen',evergreen, 'read', ncid, readvar=readv)    
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    call ncd_io('stress_decid',stress_decid, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    call ncd_io('season_decid',season_decid, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    call ncd_io('resist',resist, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    call ncd_io('pftpar20',pftpar20, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    call ncd_io('pftpar28',pftpar28, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    call ncd_io('pftpar29',pftpar29, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    call ncd_io('pftpar30',pftpar30, 'read', ncid, readvar=readv)  
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    call ncd_io('pftpar31',pftpar31, 'read', ncid, readvar=readv)  
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    call ncd_pio_closefile(ncid)

    do i = 0, numpft
       if ( trim(adjustl(pftname(i))) /= trim(expected_pftnames(i)) )then
          write(iulog,*)'pftconrd: pftname is NOT what is expected, name = ', &
                        trim(pftname(i)), ', expected name = ', trim(expected_pftnames(i))
          call endrun( 'pftconrd: bad name for pft on fpftcon dataset' )
       end if
       if ( trim(pftname(i)) == 'not_vegetated'                       ) noveg               = i
       if ( trim(pftname(i)) == 'needleleaf_evergreen_temperate_tree' ) ndllf_evr_tmp_tree  = i
       if ( trim(pftname(i)) == 'needleleaf_evergreen_boreal_tree'    ) ndllf_evr_brl_tree  = i
       if ( trim(pftname(i)) == 'needleleaf_deciduous_boreal_tree'    ) ndllf_dcd_brl_tree  = i
       if ( trim(pftname(i)) == 'broadleaf_evergreen_tropical_tree'   ) nbrdlf_evr_trp_tree  = i
       if ( trim(pftname(i)) == 'broadleaf_evergreen_temperate_tree'  ) nbrdlf_evr_tmp_tree  = i
       if ( trim(pftname(i)) == 'broadleaf_deciduous_tropical_tree'   ) nbrdlf_dcd_trp_tree  = i
       if ( trim(pftname(i)) == 'broadleaf_deciduous_temperate_tree'  ) nbrdlf_dcd_tmp_tree  = i
       if ( trim(pftname(i)) == 'broadleaf_deciduous_boreal_tree'     ) nbrdlf_dcd_brl_tree  = i
       if ( trim(pftname(i)) == 'broadleaf_evergreen_shrub'           ) nbrdlf_evr_shrub     = i
       if ( trim(pftname(i)) == 'broadleaf_deciduous_temperate_shrub' ) nbrdlf_dcd_tmp_shrub = i
       if ( trim(pftname(i)) == 'broadleaf_deciduous_boreal_shrub'    ) nbrdlf_dcd_brl_shrub = i
       if ( trim(pftname(i)) == 'c3_arctic_grass'                     ) nc3_arctic_grass     = i
       if ( trim(pftname(i)) == 'c3_non-arctic_grass'                 ) nc3_nonarctic_grass  = i
       if ( trim(pftname(i)) == 'c4_grass'                            ) nc4_grass            = i
       if ( trim(pftname(i)) == 'c3_crop'                             ) nc3crop              = i
       if ( trim(pftname(i)) == 'c3_irrigated'                        ) nirrig               = i
    end do

    ntree                = nbrdlf_dcd_brl_tree  ! value for last type of tree

    ! Set some crop-related parameters explicitly here
    ! (in future will be on pft dataset)

    crop(:)              = 0
    crop(nc3crop:numpft) = 1
    irrigated(:)         = 0
    irrigated(nirrig)    = 1

#if (defined CNDV)
    fcur(:) = fcurdv(:)
#endif

    if (masterproc) then
       write(iulog,*) 'Successfully read PFT physiological data'
       write(iulog,*)
    end if

  end subroutine pftconrd

end module pftvarcon

