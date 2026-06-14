module SurfaceRadiationMod

  !------------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Calculate solar fluxes absorbed by vegetation and ground surface
  !
  ! !USES:
  use shr_kind_mod      , only : r8 => shr_kind_r8
  use shr_log_mod       , only : errMsg => shr_log_errMsg
  use clm_varctl        , only : use_snicar_frc, use_fates
  ! [PORTED by Hui Tang: nvp (moss/lichen) control switches for radiation]
  use clm_varctl        , only : use_nvp
  use decompMod         , only : bounds_type, subgrid_level_column
  use atm2lndType       , only : atm2lnd_type
  use WaterDiagnosticBulkType    , only : waterdiagnosticbulk_type
  use CanopyStateType   , only : canopystate_type
  use SurfaceAlbedoType , only : surfalb_type
  use SolarAbsorbedType , only : solarabs_type
  use GridcellType      , only : grc
  use LandunitType      , only : lun
  use ColumnType        , only : col
  use PatchType         , only : patch
  use landunit_varcon   , only : istdlak

  ! !PRIVATE TYPES:
  implicit none
  private

  logical, parameter :: local_debug = .false.  ! for debugging this module

  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: SurfaceRadiation         ! Solar fluxes absorbed by veg and ground surface
  public :: CanopySunShadeFracs      ! Sun/Shade fractions and some area indices computations

  !
  ! !PRIVATE DATA:
  type, public :: surfrad_type
     real(r8), pointer, private  :: sfc_frc_aer_patch     (:) ! patch surface forcing of snow with all aerosols (patch) [W/m2]
     real(r8), pointer, private  :: sfc_frc_bc_patch      (:) ! patch surface forcing of snow with BC (patch) [W/m2]
     real(r8), pointer, private  :: sfc_frc_oc_patch      (:) ! patch surface forcing of snow with OC (patch) [W/m2]
     real(r8), pointer, private  :: sfc_frc_dst_patch     (:) ! patch surface forcing of snow with dust (patch) [W/m2]
     real(r8), pointer, private  :: sfc_frc_aer_sno_patch (:) ! patch surface forcing of snow with all aerosols, averaged only when snow is present (patch) [W/m2]
     real(r8), pointer, private  :: sfc_frc_bc_sno_patch  (:) ! patch surface forcing of snow with BC, averaged only when snow is present (patch) [W/m2]
     real(r8), pointer, private  :: sfc_frc_oc_sno_patch  (:) ! patch surface forcing of snow with OC, averaged only when snow is present (patch) [W/m2]
     real(r8), pointer, private  :: sfc_frc_dst_sno_patch (:) ! patch surface forcing of snow with dust, averaged only when snow is present (patch) [W/m2]

     real(r8), pointer, private  :: parveg_ln_patch       (:) ! patch  absorbed par by vegetation at local noon (W/m**2)

     real(r8), pointer, private  :: fsr_sno_vd_patch      (:) ! patch reflected direct beam vis solar radiation from snow (W/m**2)
     real(r8), pointer, private  :: fsr_sno_nd_patch      (:) ! patch reflected direct beam NIR solar radiation from snow (W/m**2)
     real(r8), pointer, private  :: fsr_sno_vi_patch      (:) ! patch reflected diffuse vis solar radiation from snow (W/m**2)
     real(r8), pointer, private  :: fsr_sno_ni_patch      (:) ! patch reflected diffuse NIR solar radiation from snow (W/m**2)

     real(r8), pointer, private  :: fsr_vis_d_patch       (:) ! patch reflected direct beam vis solar radiation (W/m**2)
     real(r8), pointer, private  :: fsr_vis_i_patch       (:) ! patch reflected diffuse vis solar radiation (W/m**2)
     real(r8), pointer, private  :: fsr_vis_d_ln_patch    (:) ! patch reflected direct beam vis solar radiation at local noon (W/m**2)
     ! diagnostic fluxes:
     real(r8), pointer, private  :: fsrSF_vis_d_patch     (:) ! snow-free patch reflected direct beam vis solar radiation (W/m**2)
     real(r8), pointer, private  :: fsrSF_vis_i_patch     (:) ! snow-free patch reflected diffuse vis solar radiation (W/m**2)
     real(r8), pointer, private  :: fsrSF_vis_d_ln_patch  (:) ! snow-free patch reflected direct beam vis solar radiation at local noon (W/m**2)
     real(r8), pointer, private  :: ssre_fsr_vis_d_patch  (:) ! snow radiative effect
     real(r8), pointer, private  :: ssre_fsr_vis_i_patch  (:) ! snow radiative effect
     real(r8), pointer, private  :: ssre_fsr_vis_d_ln_patch(:)! snow radiative effect
     real(r8), pointer, private  :: fsds_sno_vd_patch     (:) ! patch incident visible, direct radiation on snow  (for history files)  [W/m2]
     real(r8), pointer, private  :: fsds_sno_nd_patch     (:) ! patch incident near-IR, direct radiation on snow  (for history files)  [W/m2]
     real(r8), pointer, private  :: fsds_sno_vi_patch     (:) ! patch incident visible, diffuse radiation on snow (for history files) [W/m2]
     real(r8), pointer, private  :: fsds_sno_ni_patch     (:) ! patch incident near-IR, diffuse radiation on snow (for history files) [W/m2]

     real(r8), pointer, private  :: fsds_vis_d_patch      (:) ! patch incident direct beam vis solar radiation (W/m**2)
     real(r8), pointer, private  :: fsds_vis_i_patch      (:) ! patch incident diffuse vis solar radiation (W/m**2)
     real(r8), pointer, private  :: fsds_vis_d_ln_patch   (:) ! patch incident direct beam vis solar radiation at local noon (W/m**2)
     real(r8), pointer, private  :: fsds_vis_i_ln_patch   (:) ! patch incident diffuse beam vis solar radiation at local noon (W/m**2)

   contains

     procedure, public  :: Init
     procedure, private :: InitAllocate
     procedure, private :: InitHistory
     procedure, private :: InitCold

  end type surfrad_type

  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  !-----------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  subroutine Init(this, bounds)

    class(surfrad_type) :: this
    type(bounds_type), intent(in) :: bounds

    call this%InitAllocate(bounds)
    call this%InitHistory(bounds)
    call this%InitCold(bounds)

  end subroutine Init

  !-----------------------------------------------------------------------
  subroutine InitAllocate(this, bounds)
    !
    ! !USES:
    use shr_infnan_mod, only : nan => shr_infnan_nan, assignment(=)
    !
    ! !ARGUMENTS:
    class(surfrad_type) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: begp, endp
    !---------------------------------------------------------------------

    begp = bounds%begp; endp = bounds%endp

    allocate(this%sfc_frc_aer_patch     (begp:endp))              ; this%sfc_frc_aer_patch     (:)   = nan
    allocate(this%sfc_frc_bc_patch      (begp:endp))              ; this%sfc_frc_bc_patch      (:)   = nan
    allocate(this%sfc_frc_oc_patch      (begp:endp))              ; this%sfc_frc_oc_patch      (:)   = nan
    allocate(this%sfc_frc_dst_patch     (begp:endp))              ; this%sfc_frc_dst_patch     (:)   = nan
    allocate(this%sfc_frc_aer_sno_patch (begp:endp))              ; this%sfc_frc_aer_sno_patch (:)   = nan
    allocate(this%sfc_frc_bc_sno_patch  (begp:endp))              ; this%sfc_frc_bc_sno_patch  (:)   = nan
    allocate(this%sfc_frc_oc_sno_patch  (begp:endp))              ; this%sfc_frc_oc_sno_patch  (:)   = nan
    allocate(this%sfc_frc_dst_sno_patch (begp:endp))              ; this%sfc_frc_dst_sno_patch (:)   = nan

    allocate(this%parveg_ln_patch       (begp:endp))              ; this%parveg_ln_patch       (:)   = nan

    allocate(this%fsr_vis_d_patch       (begp:endp))              ; this%fsr_vis_d_patch       (:)   = nan
    allocate(this%fsr_vis_d_ln_patch    (begp:endp))              ; this%fsr_vis_d_ln_patch    (:)   = nan
    allocate(this%fsr_vis_i_patch       (begp:endp))              ; this%fsr_vis_i_patch       (:)   = nan
    allocate(this%fsrSF_vis_d_patch     (begp:endp))              ; this%fsrSF_vis_d_patch     (:)   = nan
    allocate(this%fsrSF_vis_d_ln_patch  (begp:endp))              ; this%fsrSF_vis_d_ln_patch  (:)   = nan
    allocate(this%fsrSF_vis_i_patch     (begp:endp))              ; this%fsrSF_vis_i_patch     (:)   = nan
    allocate(this%ssre_fsr_vis_d_patch  (begp:endp))              ; this%ssre_fsr_vis_d_patch  (:)   = nan
    allocate(this%ssre_fsr_vis_d_ln_patch(begp:endp))             ; this%ssre_fsr_vis_d_ln_patch(:)  = nan
    allocate(this%ssre_fsr_vis_i_patch  (begp:endp))              ; this%ssre_fsr_vis_i_patch  (:)   = nan
    allocate(this%fsr_sno_vd_patch      (begp:endp))              ; this%fsr_sno_vd_patch      (:)   = nan
    allocate(this%fsr_sno_nd_patch      (begp:endp))              ; this%fsr_sno_nd_patch      (:)   = nan
    allocate(this%fsr_sno_vi_patch      (begp:endp))              ; this%fsr_sno_vi_patch      (:)   = nan
    allocate(this%fsr_sno_ni_patch      (begp:endp))              ; this%fsr_sno_ni_patch      (:)   = nan

    allocate(this%fsds_vis_d_patch      (begp:endp))              ; this%fsds_vis_d_patch      (:)   = nan
    allocate(this%fsds_vis_i_patch      (begp:endp))              ; this%fsds_vis_i_patch      (:)   = nan
    allocate(this%fsds_vis_d_ln_patch   (begp:endp))              ; this%fsds_vis_d_ln_patch   (:)   = nan
    allocate(this%fsds_vis_i_ln_patch   (begp:endp))              ; this%fsds_vis_i_ln_patch   (:)   = nan
    allocate(this%fsds_sno_vd_patch     (begp:endp))              ; this%fsds_sno_vd_patch     (:)   = nan
    allocate(this%fsds_sno_nd_patch     (begp:endp))              ; this%fsds_sno_nd_patch     (:)   = nan
    allocate(this%fsds_sno_vi_patch     (begp:endp))              ; this%fsds_sno_vi_patch     (:)   = nan
    allocate(this%fsds_sno_ni_patch     (begp:endp))              ; this%fsds_sno_ni_patch     (:)   = nan

  end subroutine InitAllocate

  !-----------------------------------------------------------------------
  subroutine InitHistory(this, bounds)
    !
    ! History fields initialization
    !
    ! !USES:
    use shr_infnan_mod, only : nan => shr_infnan_nan, assignment(=)
    use clm_varcon    , only : spval
    use histFileMod   , only : hist_addfld1d, hist_addfld2d
    use clm_varctl    , only : use_SSRE
    !
    ! !ARGUMENTS:
    class(surfrad_type) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: begp, endp
    integer :: begc, endc
    real(r8), pointer :: data2dptr(:,:) ! temp. pointers for slicing larger arrays
    !---------------------------------------------------------------------

    begp = bounds%begp; endp = bounds%endp
    begc = bounds%begc; endc = bounds%endc

    if (use_snicar_frc) then
       this%sfc_frc_aer_patch(begp:endp) = spval
       call hist_addfld1d (fname='SNOAERFRCL', units='W/m^2', &
            avgflag='A', long_name='surface forcing of all aerosols in snow (land) ', &
            ptr_patch=this%sfc_frc_aer_patch, set_urb=spval)

       this%sfc_frc_aer_sno_patch(begp:endp) = spval
       call hist_addfld1d (fname='SNOAERFRC2L', units='W/m^2', &
            avgflag='A', long_name='surface forcing of all aerosols in snow, averaged only when snow is present (land)', &
            ptr_patch=this%sfc_frc_aer_sno_patch, set_urb=spval)

       this%sfc_frc_bc_patch(begp:endp) = spval
       call hist_addfld1d (fname='SNOBCFRCL', units='W/m^2', &
            avgflag='A', long_name='surface forcing of BC in snow (land) ', &
            ptr_patch=this%sfc_frc_bc_patch, set_urb=spval)

       this%sfc_frc_bc_sno_patch(begp:endp) = spval
       call hist_addfld1d (fname='SNOBCFRC2L', units='W/m^2', &
            avgflag='A', long_name='surface forcing of BC in snow, averaged only when snow is present (land)', &
            ptr_patch=this%sfc_frc_bc_sno_patch, set_urb=spval)

       this%sfc_frc_oc_patch(begp:endp) = spval
       call hist_addfld1d (fname='SNOOCFRCL', units='W/m^2', &
            avgflag='A', long_name='surface forcing of OC in snow (land) ', &
            ptr_patch=this%sfc_frc_oc_patch, set_urb=spval)

       this%sfc_frc_oc_sno_patch(begp:endp) = spval
       call hist_addfld1d (fname='SNOOCFRC2L', units='W/m^2', &
            avgflag='A', long_name='surface forcing of OC in snow, averaged only when snow is present (land)', &
            ptr_patch=this%sfc_frc_oc_sno_patch, set_urb=spval)

       this%sfc_frc_dst_patch(begp:endp) = spval
       call hist_addfld1d (fname='SNODSTFRCL', units='W/m^2', &
            avgflag='A', long_name='surface forcing of dust in snow (land) ', &
            ptr_patch=this%sfc_frc_dst_patch, set_urb=spval)

       this%sfc_frc_dst_sno_patch(begp:endp) = spval
       call hist_addfld1d (fname='SNODSTFRC2L', units='W/m^2', &
            avgflag='A', long_name='surface forcing of dust in snow, averaged only when snow is present (land)', &
            ptr_patch=this%sfc_frc_dst_sno_patch, set_urb=spval)
    end if

    this%fsds_vis_d_patch(begp:endp) = spval
    call hist_addfld1d (fname='FSDSVD', units='W/m^2',  &
         avgflag='A', long_name='direct vis incident solar radiation', &
         ptr_patch=this%fsds_vis_d_patch)

    this%fsds_vis_i_patch(begp:endp) = spval
    call hist_addfld1d (fname='FSDSVI', units='W/m^2',  &
         avgflag='A', long_name='diffuse vis incident solar radiation', &
         ptr_patch=this%fsds_vis_i_patch)

    this%fsr_vis_d_patch(begp:endp) = spval
    call hist_addfld1d (fname='FSRVD', units='W/m^2',  &
         avgflag='A', long_name='direct vis reflected solar radiation', &
         ptr_patch=this%fsr_vis_d_patch, c2l_scale_type='urbanf')
    this%fsr_vis_i_patch(begp:endp) = spval
    call hist_addfld1d (fname='FSRVI', units='W/m^2',  &
         avgflag='A', long_name='diffuse vis reflected solar radiation', &
         ptr_patch=this%fsr_vis_i_patch, c2l_scale_type='urbanf')
    ! diagnostic fluxes
    if (use_SSRE) then
       this%fsrSF_vis_d_patch(begp:endp) = spval
       call hist_addfld1d (fname='FSRSFVD', units='W/m^2',  &
            avgflag='A', long_name='direct vis reflected solar radiation', &
            ptr_patch=this%fsrSF_vis_d_patch, c2l_scale_type='urbanf')
       this%fsrSF_vis_i_patch(begp:endp) = spval
       call hist_addfld1d (fname='FSRSFVI', units='W/m^2',  &
            avgflag='A', long_name='diffuse vis reflected solar radiation', &
            ptr_patch=this%fsrSF_vis_i_patch, c2l_scale_type='urbanf')

       this%ssre_fsr_vis_d_patch(begp:endp) = spval
       call hist_addfld1d (fname='SSRE_FSRVD', units='W/m^2',  &
            avgflag='A', long_name='surface snow radiatve effect on direct vis reflected solar radiation', &
            ptr_patch=this%ssre_fsr_vis_d_patch, c2l_scale_type='urbanf')
       this%ssre_fsr_vis_i_patch(begp:endp) = spval
       call hist_addfld1d (fname='SSRE_FSRVI', units='W/m^2',  &
            avgflag='A', long_name='surface snow radiatve effect on diffuse vis reflected solar radiation', &
            ptr_patch=this%ssre_fsr_vis_i_patch, c2l_scale_type='urbanf')
    end if
    this%fsds_vis_d_ln_patch(begp:endp) = spval
    call hist_addfld1d (fname='FSDSVDLN', units='W/m^2',  &
         avgflag='A', long_name='direct vis incident solar radiation at local noon', &
         ptr_patch=this%fsds_vis_d_ln_patch)

    this%fsds_vis_i_ln_patch(begp:endp) = spval
    call hist_addfld1d (fname='FSDSVILN', units='W/m^2',  &
         avgflag='A', long_name='diffuse vis incident solar radiation at local noon', &
         ptr_patch=this%fsds_vis_i_ln_patch)

    this%parveg_ln_patch(begp:endp) = spval
    call hist_addfld1d (fname='PARVEGLN', units='W/m^2',  &
         avgflag='A', long_name='absorbed par by vegetation at local noon', &
         ptr_patch=this%parveg_ln_patch)

    this%fsr_vis_d_ln_patch(begp:endp) = spval
    call hist_addfld1d (fname='FSRVDLN', units='W/m^2',  &
         avgflag='A', long_name='direct vis reflected solar radiation at local noon', &
         ptr_patch=this%fsr_vis_d_ln_patch, c2l_scale_type='urbanf')
    ! diagnostic flux
    if (use_SSRE) then
       this%fsrSF_vis_d_ln_patch(begp:endp) = spval
       call hist_addfld1d (fname='FSRSFVDLN', units='W/m^2',  &
            avgflag='A', long_name='direct vis reflected solar radiation at local noon', &
            ptr_patch=this%fsrSF_vis_d_ln_patch, c2l_scale_type='urbanf')
       this%ssre_fsr_vis_d_ln_patch(begp:endp) = spval
       call hist_addfld1d (fname='SSRE_FSRVDLN', units='W/m^2',  &
            avgflag='A', long_name='surface snow radiatve effect on direct vis reflected solar radiation at local noon', &
            ptr_patch=this%ssre_fsr_vis_d_ln_patch, c2l_scale_type='urbanf')
    end if
    this%fsds_sno_vd_patch(begp:endp) = spval
    call hist_addfld1d (fname='SNOFSDSVD', units='W/m^2',  &
         avgflag='A', long_name='direct vis incident solar radiation on snow', &
         ptr_patch=this%fsds_sno_vd_patch, default='inactive')

    this%fsds_sno_nd_patch(begp:endp) = spval
    call hist_addfld1d (fname='SNOFSDSND', units='W/m^2',  &
         avgflag='A', long_name='direct nir incident solar radiation on snow', &
         ptr_patch=this%fsds_sno_nd_patch, default='inactive')

    this%fsds_sno_vi_patch(begp:endp) = spval
    call hist_addfld1d (fname='SNOFSDSVI', units='W/m^2',  &
         avgflag='A', long_name='diffuse vis incident solar radiation on snow', &
         ptr_patch=this%fsds_sno_vi_patch, default='inactive')

    this%fsds_sno_ni_patch(begp:endp) = spval
    call hist_addfld1d (fname='SNOFSDSNI', units='W/m^2',  &
         avgflag='A', long_name='diffuse nir incident solar radiation on snow', &
         ptr_patch=this%fsds_sno_ni_patch, default='inactive')

    this%fsr_sno_vd_patch(begp:endp) = spval
    call hist_addfld1d (fname='SNOFSRVD', units='W/m^2',  &
         avgflag='A', long_name='direct vis reflected solar radiation from snow', &
         ptr_patch=this%fsr_sno_vd_patch)

    this%fsr_sno_nd_patch(begp:endp) = spval
    call hist_addfld1d (fname='SNOFSRND', units='W/m^2',  &
         avgflag='A', long_name='direct nir reflected solar radiation from snow', &
         ptr_patch=this%fsr_sno_nd_patch)

    this%fsr_sno_vi_patch(begp:endp) = spval
    call hist_addfld1d (fname='SNOFSRVI', units='W/m^2',  &
         avgflag='A', long_name='diffuse vis reflected solar radiation from snow', &
         ptr_patch=this%fsr_sno_vi_patch)

    this%fsr_sno_ni_patch(begp:endp) = spval
    call hist_addfld1d (fname='SNOFSRNI', units='W/m^2',  &
         avgflag='A', long_name='diffuse nir reflected solar radiation from snow', &
         ptr_patch=this%fsr_sno_ni_patch)


  end subroutine InitHistory

  !------------------------------------------------------------------------
  subroutine InitCold(this, bounds)
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(surfrad_type) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: p,l
    !-----------------------------------------------------------------------

    ! nothing for now

  end subroutine InitCold


  subroutine CanopySunShadeFracs(filter_nourbanp, num_nourbanp,  &
                                 atm2lnd_inst, surfalb_inst,     &
                                 canopystate_inst, solarabs_inst)

    ! ------------------------------------------------------------------------------------
    ! This subroutine calculates and returns patch vectors of
    !
    ! 1) absorbed PAR for sunlit leaves in canopy layer
    ! 2) absorbed PAR for shaded leaves in canopy layer
    ! 3) sunlit leaf area
    ! 4) shaded  leaf area
    ! 5) sunlit leaf area for canopy layer
    ! 6) shaded leaf area for canopy layer
    ! 7) sunlit fraction of canopy
    !
    ! This routine has a counterpart when the fates model is turned on.
    ! CLMEDInterf_CanopySunShadeFracs()
    ! If changes are applied to this routine, please take a moment to review that
    ! subroutine as well and consider if any new information related to these types of
    ! variables also needs to be augmented in that routine as well.
    ! ------------------------------------------------------------------------------------


    implicit none

    ! Arguments (in)

    integer, intent(in),dimension(:)      :: filter_nourbanp    ! patch filter for non-urban points
    integer, intent(in)                   :: num_nourbanp       ! size of the nonurban filter
    type(atm2lnd_type), intent(in)        :: atm2lnd_inst
    type(surfalb_type), intent(in)        :: surfalb_inst

    ! Arguments (inout)
    type(canopystate_type), intent(inout) :: canopystate_inst
    type(solarabs_type), intent(inout)    :: solarabs_inst

    ! local variables
    integer           :: fp                         ! non-urban filter patch index
    integer           :: p                          ! patch index
    integer           :: c                          ! column index
    integer           :: g                          ! gridcell index
    integer           :: iv                         ! canopy layer index
    integer,parameter :: ipar = 1                   ! The band index for PAR

    associate( tlai_z  => surfalb_inst%tlai_z_patch, &    ! tlai increment for canopy layer
          fsun_z      => surfalb_inst%fsun_z_patch, &     ! sunlit fraction of canopy layer
          elai        => canopystate_inst%elai_patch, &   ! one-sided leaf area index
          forc_solad_col  => atm2lnd_inst%forc_solad_downscaled_col, &   ! direct beam radiation, column (W/m**2)
          forc_solai  => atm2lnd_inst%forc_solai_grc, &   ! diffuse radiation (W/m**2)
          fabd_sun_z  => surfalb_inst%fabd_sun_z_patch, & ! absorbed sunlit leaf direct PAR
          fabd_sha_z  => surfalb_inst%fabd_sha_z_patch, & ! absorbed shaded leaf direct PAR
          fabi_sun_z  => surfalb_inst%fabi_sun_z_patch, & ! absorbed sunlit leaf diffuse PAR
          fabi_sha_z  => surfalb_inst%fabi_sha_z_patch, & ! absorbed shaded leaf diffuse PAR
          nrad        => surfalb_inst%nrad_patch, &       ! number of canopy layers
          parsun_z    => solarabs_inst%parsun_z_patch, &  ! absorbed PAR for sunlit leaves
          parsha_z    => solarabs_inst%parsha_z_patch, &  ! absorbed PAR for shaded leaves
          laisun      => canopystate_inst%laisun_patch, & ! sunlit leaf area
          laisha      => canopystate_inst%laisha_patch, & ! shaded  leaf area
          laisun_z    => canopystate_inst%laisun_z_patch, & ! sunlit leaf area for canopy layer
          laisha_z    => canopystate_inst%laisha_z_patch, & ! shaded leaf area for canopy layer
          fsun        => canopystate_inst%fsun_patch)       ! sunlit fraction of canopy

     do fp = 1,num_nourbanp

        p = filter_nourbanp(fp)

        do iv = 1, nrad(p)
           parsun_z(p,iv) = 0._r8
           parsha_z(p,iv) = 0._r8
           laisun_z(p,iv) = 0._r8
           laisha_z(p,iv) = 0._r8
        end do

        ! Loop over patches to calculate laisun_z and laisha_z for each layer.
        ! Derive canopy laisun, laisha, and fsun from layer sums.
        ! If sun/shade big leaf code, nrad=1 and fsun_z(p,1) and tlai_z(p,1) from
        ! SurfaceAlbedo is canopy integrated so that layer value equals canopy value.

        laisun(p) = 0._r8
        laisha(p) = 0._r8
        do iv = 1, nrad(p)
           laisun_z(p,iv) = tlai_z(p,iv) * fsun_z(p,iv)
           laisha_z(p,iv) = tlai_z(p,iv) * (1._r8 - fsun_z(p,iv))
           laisun(p) = laisun(p) + laisun_z(p,iv)
           laisha(p) = laisha(p) + laisha_z(p,iv)
        end do
        if (elai(p) > 0._r8) then
           fsun(p) = laisun(p) / elai(p)
        else
           fsun(p) = 0._r8
        end if

        ! Absorbed PAR profile through canopy
        ! If sun/shade big leaf code, nrad=1 and fluxes from SurfaceAlbedo
        ! are canopy integrated so that layer values equal big leaf values.

        g = patch%gridcell(p)
        c = patch%column(p)

        do iv = 1, nrad(p)
           parsun_z(p,iv) = forc_solad_col(c,ipar)*fabd_sun_z(p,iv) + forc_solai(g,ipar)*fabi_sun_z(p,iv)
           parsha_z(p,iv) = forc_solad_col(c,ipar)*fabd_sha_z(p,iv) + forc_solai(g,ipar)*fabi_sha_z(p,iv)
        end do

     end do ! end of fp = 1,num_nourbanp loop
   end associate
   return
 end subroutine CanopySunShadeFracs

  !------------------------------------------------------------------------------
  subroutine SurfaceRadiation(bounds, num_nourbanp, filter_nourbanp, &
       num_urbanp, filter_urbanp, num_urbanc, filter_urbanc, &
       atm2lnd_inst, waterdiagnosticbulk_inst, canopystate_inst, &
       surfalb_inst, solarabs_inst, surfrad_inst)
     !
     ! !DESCRIPTION:
     ! Solar fluxes absorbed by vegetation and ground surface
     ! Note possible problem when land is on different grid than atmosphere.
     ! Land may have sun above the horizon (coszen > 0) but atmosphere may
     ! have sun below the horizon (forc_solad = 0 and forc_solai = 0). This is okay
     ! because all fluxes (absorbed, reflected, transmitted) are multiplied
     ! by the incoming flux and all will equal zero.
     ! Atmosphere may have sun above horizon (forc_solad > 0 and forc_solai > 0) but
     ! land may have sun below horizon. This is okay because fabd, fabi,
     ! ftdd, ftid, and ftii all equal zero so that sabv=sabg=fsa=0. Also,
     ! albd and albi equal one so that fsr=forc_solad+forc_solai. In other words, all
     ! the radiation is reflected. NDVI should equal zero in this case.
     ! However, the way the code is currently implemented this is only true
     ! if (forc_solad+forc_solai)|vis = (forc_solad+forc_solai)|nir.
     ! Output variables are parsun,parsha,sabv,sabg,fsa,fsr,ndvi
     !
     ! !USES:
     use clm_varpar       , only : numrad, nlevsno
     use clm_varcon       , only : spval
     use landunit_varcon  , only : istsoil, istcrop 
     use clm_varctl       , only : use_subgrid_fluxes, use_snicar_frc, iulog, use_SSRE, do_sno_oc
     use clm_time_manager , only : get_step_size_real, is_near_local_noon, get_nstep  ! [PORTED by Hui Tang: get_nstep for Phase-4 diagnostic]
     use abortutils       , only : endrun
     !
     ! !ARGUMENTS:
     type(bounds_type)      , intent(in)            :: bounds
     integer                , intent(in)            :: num_nourbanp       ! number of patches in non-urban points in patch  filter
     integer                , intent(in)            :: filter_nourbanp(:) ! patch filter for non-urban points
     integer                , intent(in)            :: num_urbanp         ! number of patches in non-urban points in patch filter
     integer                , intent(in)            :: filter_urbanp(:)   ! patch filter for non-urban points
     integer                , intent(in)            :: num_urbanc         ! number of urban columns in clump
     integer                , intent(in)            :: filter_urbanc(:)   ! urban column filter
     type(atm2lnd_type)     , intent(in)            :: atm2lnd_inst
     type(waterdiagnosticbulk_type)  , intent(in)            :: waterdiagnosticbulk_inst
     type(surfalb_type)     , intent(in)            :: surfalb_inst
     type(canopystate_type) , intent(inout)         :: canopystate_inst
     type(solarabs_type)    , intent(inout)         :: solarabs_inst
     type(surfrad_type)     , intent(inout)         :: surfrad_inst
     !
     ! !LOCAL VARIABLES:
     integer , parameter :: nband = numrad           ! number of solar radiation waveband classes
     real(r8), parameter :: mpe = 1.e-06_r8          ! prevents overflow for division by zero
     integer  :: fp                                  ! non-urban filter patch index
     integer  :: p                                   ! patch index
     integer  :: c                                   ! column index
     integer  :: l                                   ! landunit index
     integer  :: g                                   ! grid cell index
     integer  :: ib                                  ! waveband number (1=vis, 2=nir)
     integer  :: iv                                  ! canopy layer
     real(r8) :: absrad                              ! absorbed solar radiation (W/m**2)
     integer  :: i                                   ! layer index [idx]
     real(r8) :: rnir                                ! reflected solar radiation [nir] (W/m**2)
     real(r8) :: rvis                                ! reflected solar radiation [vis] (W/m**2)
     real(r8) :: rnirSF                              ! snow-free reflected solar radiation [nir] (W/m**2)
     real(r8) :: rvisSF                              ! snow-free reflected solar radiation [vis] (W/m**2)
     real(r8) :: trd(bounds%begp:bounds%endp,numrad) ! transmitted solar radiation: direct (W/m**2)
     real(r8) :: tri(bounds%begp:bounds%endp,numrad) ! transmitted solar radiation: diffuse (W/m**2)
     real(r8) :: cad(bounds%begp:bounds%endp,numrad) ! direct beam absorbed by canopy (W/m**2)
     real(r8) :: cai(bounds%begp:bounds%endp,numrad) ! diffuse radiation absorbed by canopy (W/m**2)
     real(r8) :: dtime                               ! land model time step (sec)
     real(r8) :: sabg_snl_sum                        ! temporary, absorbed energy in all active snow layers [W/m2]
     ! [PORTED by Hui Tang: partial-snow NVP blend locals]
     real(r8) :: frac_nvp_eff_loc                    ! locally-computed exposed NVP area fraction
     real(r8) :: f_exp_loc                           ! fraction of NVP area that is exposed (not snow-buried)
     real(r8) :: sabg_nvp_beer                       ! Beer's law NVP absorption per unit column area [W/m2]
     real(r8) :: sabg_sum_chk                        ! [PORTED by Hui Tang: sum(sabg_lyr) (+ sabg_nvp for snl==0) for the SNICAR conservation check]
     real(r8) :: absrad_pur                          ! temp: absorbed solar radiation by pure snow [W/m2]
     real(r8) :: absrad_bc                           ! temp: absorbed solar radiation without BC [W/m2]
     real(r8) :: absrad_oc                           ! temp: absorbed solar radiation without OC [W/m2]
     real(r8) :: absrad_dst                          ! temp: absorbed solar radiation without dust [W/m2]
     real(r8) :: sabg_pur(bounds%begp:bounds%endp)   ! solar radiation absorbed by ground with pure snow [W/m2]
     real(r8) :: sabg_bc(bounds%begp:bounds%endp)    ! solar radiation absorbed by ground without BC [W/m2]
     real(r8) :: sabg_oc(bounds%begp:bounds%endp)    ! solar radiation absorbed by ground without OC [W/m2]
     real(r8) :: sabg_dst(bounds%begp:bounds%endp)   ! solar radiation absorbed by ground without dust [W/m2]
     real(r8) :: parveg(bounds%begp:bounds%endp)     ! absorbed par by vegetation (W/m**2)
     !
     !------------------------------------------------------------------------------

     associate(                                                     &
          snl             =>    col%snl                           , & ! Input:  [integer  (:)   ] negative number of snow layers [nbr]

          forc_solad_col  =>    atm2lnd_inst%forc_solad_downscaled_col     , & ! Input:  [real(r8) (:,:) ] direct beam radiation, column (W/m**2)
          forc_solai      =>    atm2lnd_inst%forc_solai_grc       , & ! Input:  [real(r8) (:,:) ] diffuse radiation (W/m**2)

          snow_depth      =>    waterdiagnosticbulk_inst%snow_depth_col    , & ! Input:  [real(r8) (:)   ] snow height (m)
          frac_sno        =>    waterdiagnosticbulk_inst%frac_sno_col      , & ! Input:  [real(r8) (:)   ] fraction of ground covered by snow (0 to 1)

          nrad            =>    surfalb_inst%nrad_patch           , & ! Input:  [integer  (:)   ] number of canopy layers, above snow for radiative transfer
          coszen          =>    surfalb_inst%coszen_col           , & ! Input:  [real(r8) (:)   ] column cosine of solar zenith angle
          albgrd          =>    surfalb_inst%albgrd_col           , & ! Input:  [real(r8) (:,:) ] ground albedo (direct)
          albgri          =>    surfalb_inst%albgri_col           , & ! Input:  [real(r8) (:,:) ] ground albedo (diffuse)
          albsod          =>    surfalb_inst%albsod_col           , & ! Input:  [real(r8) (:,:) ] direct-beam soil albedo (col,bnd) [frc]
          albgrd_oc       =>    surfalb_inst%albgrd_oc_col        , & ! Input:  [real(r8) (:,:) ] ground albedo without OC (direct) (col,bnd)
          albgri_oc       =>    surfalb_inst%albgri_oc_col        , & ! Input:  [real(r8) (:,:) ] ground albedo without OC (diffuse) (col,bnd)
          albgrd_dst      =>    surfalb_inst%albgrd_dst_col       , & ! Input:  [real(r8) (:,:) ] ground albedo without dust (direct) (col,bnd)
          albgri_dst      =>    surfalb_inst%albgri_dst_col       , & ! Input:  [real(r8) (:,:) ] ground albedo without dust (diffuse) (col,bnd)
          albsnd_hst      =>    surfalb_inst%albsnd_hst_col       , & ! Input:  [real(r8) (:,:) ] snow albedo, direct, for history files (col,bnd) [frc]
          albsni_hst      =>    surfalb_inst%albsni_hst_col       , & ! Input:  [real(r8) (:,:) ] snow ground albedo, diffuse, for history files (col,bnd
          flx_absdv       =>    surfalb_inst%flx_absdv_col        , & ! Input:  [real(r8) (:,:) ] direct flux absorption factor (col,lyr): VIS [frc]
          flx_absdn       =>    surfalb_inst%flx_absdn_col        , & ! Input:  [real(r8) (:,:) ] direct flux absorption factor (col,lyr): NIR [frc]
          flx_absiv       =>    surfalb_inst%flx_absiv_col        , & ! Input:  [real(r8) (:,:) ] diffuse flux absorption factor (col,lyr): VIS [frc]
          flx_absin       =>    surfalb_inst%flx_absin_col        , & ! Input:  [real(r8) (:,:) ] diffuse flux absorption factor (col,lyr): NIR [frc]
          albsoi          =>    surfalb_inst%albsoi_col           , & ! Input:  [real(r8) (:,:) ] diffuse soil albedo (col,bnd) [frc]
          albd            =>    surfalb_inst%albd_patch           , & ! Input:  [real(r8) (:,:) ] surface albedo (direct)
          albi            =>    surfalb_inst%albi_patch           , & ! Input:  [real(r8) (:,:) ] surface albedo (diffuse)
          albdSF          =>    surfalb_inst%albdSF_patch         , & ! Input:  [real(r8) (:,:) ] snow-free surface albedo (direct)
          albiSF          =>    surfalb_inst%albiSF_patch         , & ! Input:  [real(r8) (:,:) ] snow-free surface albedo (diffuse)
          fabd            =>    surfalb_inst%fabd_patch           , & ! Input:  [real(r8) (:,:) ] flux absorbed by canopy per unit direct flux
          fabd_sun        =>    surfalb_inst%fabd_sun_patch       , & ! Input:  [real(r8) (:,:) ] flux absorbed by sunlit canopy per unit direct flux
          fabd_sha        =>    surfalb_inst%fabd_sha_patch       , & ! Input:  [real(r8) (:,:) ] flux absorbed by shaded canopy per unit direct flux
          fabi            =>    surfalb_inst%fabi_patch           , & ! Input:  [real(r8) (:,:) ] flux absorbed by canopy per unit diffuse flux
          fabi_sun        =>    surfalb_inst%fabi_sun_patch       , & ! Input:  [real(r8) (:,:) ] flux absorbed by sunlit canopy per unit diffuse flux
          fabi_sha        =>    surfalb_inst%fabi_sha_patch       , & ! Input:  [real(r8) (:,:) ] flux absorbed by shaded canopy per unit diffuse flux
          ftdd            =>    surfalb_inst%ftdd_patch           , & ! Input:  [real(r8) (:,:) ] down direct flux below canopy per unit direct flux
          ftid            =>    surfalb_inst%ftid_patch           , & ! Input:  [real(r8) (:,:) ] down diffuse flux below canopy per unit direct flux
          ftii            =>    surfalb_inst%ftii_patch           , & ! Input:  [real(r8) (:,:) ] down diffuse flux below canopy per unit diffuse flux
          fabd_sun_z      =>    surfalb_inst%fabd_sun_z_patch     , & ! Input:  [real(r8) (:,:) ] absorbed sunlit leaf direct  PAR (per unit lai+sai) for each canopy layer
          fabd_sha_z      =>    surfalb_inst%fabd_sha_z_patch     , & ! Input:  [real(r8) (:,:) ] absorbed shaded leaf direct  PAR (per unit lai+sai) for each canopy layer
          fabi_sun_z      =>    surfalb_inst%fabi_sun_z_patch     , & ! Input:  [real(r8) (:,:) ] absorbed sunlit leaf diffuse PAR (per unit lai+sai) for each canopy layer
          fabi_sha_z      =>    surfalb_inst%fabi_sha_z_patch     , & ! Input:  [real(r8) (:,:) ] absorbed shaded leaf diffuse PAR (per unit lai+sai) for each canopy layer
          albgrd_pur      =>    surfalb_inst%albgrd_pur_col       , & ! Input:  [real(r8) (:,:) ] pure snow ground albedo (direct)
          albgri_pur      =>    surfalb_inst%albgri_pur_col       , & ! Input:  [real(r8) (:,:) ] pure snow ground albedo (diffuse)
          albgrd_bc       =>    surfalb_inst%albgrd_bc_col        , & ! Input:  [real(r8) (:,:) ] ground albedo without BC (direct) (col,bnd)
          albgri_bc       =>    surfalb_inst%albgri_bc_col        , & ! Input:  [real(r8) (:,:) ] ground albedo without BC (diffuse) (col,bnd)
          tlai            =>    canopystate_inst%tlai_patch       , & ! Input:  [real(r8) (:)   ] one-sided leaf area index
          elai            =>    canopystate_inst%elai_patch       , & ! Input:  [real(r8) (:)   ] one-sided leaf area index with burying by snow
          esai            =>    canopystate_inst%esai_patch       , & ! Input:  [real(r8) (:)   ] one-sided stem area index with burying by snow
          fsun            =>    canopystate_inst%fsun_patch       , & ! Output: [real(r8) (:)   ] sunlit fraction of canopy
          fsa             =>    solarabs_inst%fsa_patch           , & ! Output: [real(r8) (:)   ] solar radiation absorbed (total) (W/m**2)
          fsr             =>    solarabs_inst%fsr_patch           , & ! Output: [real(r8) (:)   ] solar radiation reflected (W/m**2)
          fsrSF           =>    solarabs_inst%fsrSF_patch         , & ! Output: [real(r8) (:)   ] diagnostic snow-free solar radiation reflected (W/m**2)
          ssre_fsr        =>    solarabs_inst%ssre_fsr_patch      , & ! Output: [real(r8) (:)   ] diagnostic snow-free solar radiation reflected (W/m**2)
          sabv            =>    solarabs_inst%sabv_patch          , & ! Output: [real(r8) (:)   ] solar radiation absorbed by vegetation (W/m**2)
          sabg            =>    solarabs_inst%sabg_patch          , & ! Output: [real(r8) (:)   ] solar radiation absorbed by ground (W/m**2)
          sabg_pen        =>    solarabs_inst%sabg_pen_patch      , & ! Output: [real(r8) (:)   ] solar (rural) radiation penetrating top soisno layer (W/m**2)
          sabg_soil       =>    solarabs_inst%sabg_soil_patch     , & ! Output: [real(r8) (:)   ] solar radiation absorbed by soil (W/m**2)
          sabg_nvp        =>    solarabs_inst%sabg_nvp_patch      , & ! [PORTED by Hui Tang: solar absorbed by exposed NVP moss surface (W/m**2)]
          sabg_soil_bandloop =>  solarabs_inst%sabg_soil_bandloop_patch , & ! [PORTED by Hui Tang: band-loop ground absorption snapshot]
          sabg_snow       =>    solarabs_inst%sabg_snow_patch     , & ! Output: [real(r8) (:)   ] solar radiation absorbed by snow (W/m**2)
          sabg_lyr        =>    solarabs_inst%sabg_lyr_patch      , & ! Output: [real(r8) (:,:) ] absorbed radiative flux (patch,lyr) [W/m2]
          fsr_nir_d       =>    solarabs_inst%fsr_nir_d_patch     , & ! Output: [real(r8) (:)   ] reflected direct beam nir solar radiation (W/m**2)
          fsr_nir_i       =>    solarabs_inst%fsr_nir_i_patch     , & ! Output: [real(r8) (:)   ] reflected diffuse nir solar radiation (W/m**2)
          fsr_nir_d_ln    =>    solarabs_inst%fsr_nir_d_ln_patch  , & ! Output: [real(r8) (:)   ] reflected direct beam nir solar rad at local noon (W/m**2)
          fsds_nir_d      =>    solarabs_inst%fsds_nir_d_patch    , & ! Output: [real(r8) (:)   ] incident direct beam nir solar radiation (W/m**2)
          fsds_nir_d_ln   =>    solarabs_inst%fsds_nir_d_ln_patch , & ! Output: [real(r8) (:)   ] incident direct beam nir solar rad at local noon (W/m**2)
          fsds_nir_i      =>    solarabs_inst%fsds_nir_i_patch    , & ! Output: [real(r8) (:)   ] incident diffuse nir solar radiation (W/m**2)
          fsrSF_nir_d     =>    solarabs_inst%fsrSF_nir_d_patch   , & ! Output: [real(r8) (:)   ] snow-free reflected direct beam nir solar radiation (W/m**2)
          fsrSF_nir_i     =>    solarabs_inst%fsrSF_nir_i_patch   , & ! Output: [real(r8) (:)   ] snow-free reflected diffuse nir solar radiation (W/m**2)
          fsrSF_nir_d_ln  =>    solarabs_inst%fsrSF_nir_d_ln_patch, & ! Output: [real(r8) (:)   ] snow-free reflected direct beam nir solar rad at local noon (W/m**2)
          ssre_fsr_nir_d  =>    solarabs_inst%ssre_fsr_nir_d_patch, & ! Output: [real(r8) (:)   ] snow-free reflected direct beam nir solar radiation (W/m**2)
          ssre_fsr_nir_i  =>    solarabs_inst%ssre_fsr_nir_i_patch, & ! Output: [real(r8) (:)   ] snow-free reflected diffuse nir solar radiation (W/m**2)
          ssre_fsr_nir_d_ln=>   solarabs_inst%ssre_fsr_nir_d_ln_patch,&!Output: [real(r8) (:)   ] snow-free reflected direct beam nir solar rad at local noon (W/m**2)
          fsa_r           =>    solarabs_inst%fsa_r_patch         , & ! Output: [real(r8) (:)   ] rural solar radiation absorbed (total) (W/m**2)
          sub_surf_abs_SW =>    solarabs_inst%sub_surf_abs_SW_patch,& ! Output: [real(r8) (:)   ] fraction of solar radiation absorbed below first snow layer (W/M**2)

          parveg_ln       =>    surfrad_inst%parveg_ln_patch      , & ! Output: [real(r8) (:)   ] absorbed par by vegetation at local noon (W/m**2)
          fsr_vis_d       =>    surfrad_inst%fsr_vis_d_patch      , & ! Output: [real(r8) (:)   ] reflected direct beam vis solar radiation (W/m**2)
          fsr_vis_i       =>    surfrad_inst%fsr_vis_i_patch      , & ! Output: [real(r8) (:)   ] reflected diffuse vis solar radiation (W/m**2)
          fsrSF_vis_d     =>    surfrad_inst%fsrSF_vis_d_patch    , & ! Output: [real(r8) (:)   ] snow-free reflected direct beam vis solar radiation (W/m**2)
          fsrSF_vis_i     =>    surfrad_inst%fsrSF_vis_i_patch    , & ! Output: [real(r8) (:)   ] snow-free reflected diffuse vis solar radiation (W/m**2)
          ssre_fsr_vis_d  =>    surfrad_inst%ssre_fsr_vis_d_patch , & ! Output: [real(r8) (:)   ] snow-free reflected direct beam vis solar radiation (W/m**2)
          ssre_fsr_vis_i  =>    surfrad_inst%ssre_fsr_vis_i_patch , & ! Output: [real(r8) (:)   ] snow-free reflected diffuse vis solar radiation (W/m**2)
          fsds_vis_i_ln   =>    surfrad_inst%fsds_vis_i_ln_patch  , & ! Output: [real(r8) (:)   ] incident diffuse beam vis solar rad at local noon (W/m**2)
          fsr_vis_d_ln    =>    surfrad_inst%fsr_vis_d_ln_patch   , & ! Output: [real(r8) (:)   ] reflected direct beam vis solar rad at local noon (W/m**2)
          fsrSF_vis_d_ln  =>    surfrad_inst%fsrSF_vis_d_ln_patch , & ! Output: [real(r8) (:)   ] snow-free reflected direct beam vis solar rad at local noon (W/m**2)
          fsds_vis_d      =>    surfrad_inst%fsds_vis_d_patch     , & ! Output: [real(r8) (:)   ] incident direct beam vis solar radiation (W/m**2)
          fsds_vis_i      =>    surfrad_inst%fsds_vis_i_patch     , & ! Output: [real(r8) (:)   ] incident diffuse vis solar radiation (W/m**2)
          fsds_vis_d_ln   =>    surfrad_inst%fsds_vis_d_ln_patch  , & ! Output: [real(r8) (:)   ] incident direct beam vis solar rad at local noon (W/m**2)
          sfc_frc_aer     =>    surfrad_inst%sfc_frc_aer_patch    , & ! Output: [real(r8) (:)   ] surface forcing of snow with all aerosols (patch) [W/m2]
          sfc_frc_aer_sno =>    surfrad_inst%sfc_frc_aer_sno_patch, & ! Output: [real(r8) (:)   ] surface forcing of snow with all aerosols, averaged only when snow is present (patch) [W/m2]
          sfc_frc_bc      =>    surfrad_inst%sfc_frc_bc_patch     , & ! Output: [real(r8) (:)   ] surface forcing of snow with BC (patch) [W/m2]
          sfc_frc_bc_sno  =>    surfrad_inst%sfc_frc_bc_sno_patch , & ! Output: [real(r8) (:)   ] surface forcing of snow with BC, averaged only when snow is present (patch) [W/m2]
          sfc_frc_oc      =>    surfrad_inst%sfc_frc_oc_patch     , & ! Output: [real(r8) (:)   ] surface forcing of snow with OC (patch) [W/m2]
          sfc_frc_oc_sno  =>    surfrad_inst%sfc_frc_oc_sno_patch , & ! Output: [real(r8) (:)   ] surface forcing of snow with OC, averaged only when snow is present (patch) [W/m2]
          sfc_frc_dst     =>    surfrad_inst%sfc_frc_dst_patch    , & ! Output: [real(r8) (:)   ] surface forcing of snow with dust (patch) [W/m2]
          sfc_frc_dst_sno =>    surfrad_inst%sfc_frc_dst_sno_patch, & ! Output: [real(r8) (:)   ] surface forcing of snow with dust, averaged only when snow is present (patch) [W/m2]
          fsr_sno_vd      =>    surfrad_inst%fsr_sno_vd_patch     , & ! Output: [real(r8) (:)   ] reflected visible, direct radiation from snow (for history files) (patch) [W/m2]
          fsr_sno_nd      =>    surfrad_inst%fsr_sno_nd_patch     , & ! Output: [real(r8) (:)   ] reflected near-IR, direct radiation from snow (for history files) (patch) [W/m2]
          fsr_sno_vi      =>    surfrad_inst%fsr_sno_vi_patch     , & ! Output: [real(r8) (:)   ] reflected visible, diffuse radiation from snow (for history files) (patch) [W/m2]
          fsr_sno_ni      =>    surfrad_inst%fsr_sno_ni_patch     , & ! Output: [real(r8) (:)   ] reflected near-IR, diffuse radiation from snow (for history files) (patch) [W/m2]
          fsds_sno_vd     =>    surfrad_inst%fsds_sno_vd_patch    , & ! Output: [real(r8) (:)   ] incident visible, direct radiation on snow (for history files) (patch) [W/m2]
          fsds_sno_nd     =>    surfrad_inst%fsds_sno_nd_patch    , & ! Output: [real(r8) (:)   ] incident near-IR, direct radiation on snow (for history files) (patch) [W/m2]
          fsds_sno_vi     =>    surfrad_inst%fsds_sno_vi_patch    , & ! Output: [real(r8) (:)   ] incident visible, diffuse radiation on snow (for history files) (patch) [W/m2]
          fsds_sno_ni     =>    surfrad_inst%fsds_sno_ni_patch    , & ! Output: [real(r8) (:)   ] incident near-IR, diffuse radiation on snow (for history files) (patch) [W/m2]
          frac_sno    => waterdiagnosticbulk_inst%frac_sno_col      & ! Input:  [real(r8)  (:)   ]  fraction of ground covered by snow (0 to 1)

          )

       ! Determine seconds off current time step
       dtime = get_step_size_real()
y
       ! Initialize fluxes

       do fp = 1,num_nourbanp
          p = filter_nourbanp(fp)
          l = patch%landunit(p)
          g = patch%gridcell(p)

          sabg_soil(p)  = 0._r8
          sabg_soil_bandloop(p) = 0._r8   ! [PORTED by Hui Tang]
          if (use_nvp) sabg_nvp(p) = 0._r8   ! [PORTED by Hui Tang: exposed-NVP surface solar]
          sabg_snow(p)  = 0._r8
          sabg(p)       = 0._r8
          sabv(p)       = 0._r8
          fsa(p)        = 0._r8
          if (lun%itype(l)==istsoil .or. lun%itype(l)==istcrop) then
             fsa_r(p) = 0._r8
          end if
          sabg_lyr(p,:) = 0._r8
          sabg_pur(p)   = 0._r8
          sabg_bc(p)    = 0._r8
          sabg_oc(p)    = 0._r8
          sabg_dst(p)   = 0._r8

       end do

       ! zero-out fsun for the urban patches
       ! the non-urban patches were set prior to this call
       ! and split into fates and non-fates specific functions
       do fp = 1,num_urbanp
          p = filter_urbanp(fp)
          fsun(p) = 0._r8
       end do

       ! Loop over nband wavebands
       do ib = 1, nband
          do fp = 1,num_nourbanp
             p = filter_nourbanp(fp)
             c = patch%column(p)
             l = patch%landunit(p)
             g = patch%gridcell(p)

             ! Absorbed by canopy

             cad(p,ib) = forc_solad_col(c,ib)*fabd(p,ib)
             cai(p,ib) = forc_solai(g,ib)*fabi(p,ib)
             sabv(p) = sabv(p) + cad(p,ib) + cai(p,ib)
             fsa(p)  = fsa(p)  + cad(p,ib) + cai(p,ib)
             if (ib == 1) then
                parveg(p) = cad(p,ib) + cai(p,ib)
             end if
             if (lun%itype(l)==istsoil .or. lun%itype(l)==istcrop) then
                fsa_r(p)  = fsa_r(p)  + cad(p,ib) + cai(p,ib)
             end if

             ! Transmitted = solar fluxes incident on ground

             trd(p,ib) = forc_solad_col(c,ib)*ftdd(p,ib)
             tri(p,ib) = forc_solad_col(c,ib)*ftid(p,ib) + forc_solai(g,ib)*ftii(p,ib)
             ! Solar radiation absorbed by ground surface
             ! calculate absorbed solar by soil/snow separately
             absrad  = trd(p,ib)*(1._r8-albsod(c,ib)) + tri(p,ib)*(1._r8-albsoi(c,ib))
             sabg_soil(p) = sabg_soil(p) + absrad
             absrad  = trd(p,ib)*(1._r8-albsnd_hst(c,ib)) + tri(p,ib)*(1._r8-albsni_hst(c,ib))
             sabg_snow(p) = sabg_snow(p) + absrad
             absrad  = trd(p,ib)*(1._r8-albgrd(c,ib)) + tri(p,ib)*(1._r8-albgri(c,ib))
             sabg(p) = sabg(p) + absrad
             fsa(p)  = fsa(p)  + absrad
             if (lun%itype(l)==istsoil .or. lun%itype(l)==istcrop) then
                fsa_r(p)  = fsa_r(p)  + absrad
             end if
             if (snl(c) == 0) then
                sabg_snow(p) = sabg(p)
                sabg_soil(p) = sabg(p)
             endif
             ! if no subgrid fluxes, make sure to set both components equal to weighted average
             if (.not. use_subgrid_fluxes .or. lun%itype(l) == istdlak) then
                sabg_snow(p) = sabg(p)
                sabg_soil(p) = sabg(p)
             endif

             if (use_snicar_frc) then
                ! Solar radiation absorbed by ground surface without BC
                absrad_bc = trd(p,ib)*(1._r8-albgrd_bc(c,ib)) + tri(p,ib)*(1._r8-albgri_bc(c,ib))
                sabg_bc(p) = sabg_bc(p) + absrad_bc

                ! Solar radiation absorbed by ground surface without OC
                absrad_oc = trd(p,ib)*(1._r8-albgrd_oc(c,ib)) + tri(p,ib)*(1._r8-albgri_oc(c,ib))
                sabg_oc(p) = sabg_oc(p) + absrad_oc

                ! Solar radiation absorbed by ground surface without dust
                absrad_dst = trd(p,ib)*(1._r8-albgrd_dst(c,ib)) + tri(p,ib)*(1._r8-albgri_dst(c,ib))
                sabg_dst(p) = sabg_dst(p) + absrad_dst

                ! Solar radiation absorbed by ground surface without any aerosols
                absrad_pur = trd(p,ib)*(1._r8-albgrd_pur(c,ib)) + tri(p,ib)*(1._r8-albgri_pur(c,ib))
                sabg_pur(p) = sabg_pur(p) + absrad_pur
             end if

          end do ! end of patch loop
       end do ! end nbands loop

       !   compute absorbed flux in each snow layer and top soil layer,
       !   based on flux factors computed in the radiative transfer portion of SNICAR.

       do fp = 1,num_nourbanp
          p = filter_nourbanp(fp)
          c = patch%column(p)
          l = patch%landunit(p)
          sabg_snl_sum = 0._r8

          sub_surf_abs_SW(p) = 0._r8

          ! [PORTED by Hui Tang: snapshot the band-loop (albsod) ground absorption BEFORE the NVP
          !  carve-out (CASE1, sabg_soil -= sabg_lyr0) and the SNICAR snow reassignment (CASE2,
          !  sabg_soil = sabg_lyr(p,1)) overwrite sabg_soil. This raw value lumps NVP+soil absorption
          !  for the exposed surface and is used only to build the diagnostic SABG tile in
          !  SoilTemperatureMod (true exposed-surface absorption during melt). Diagnostic-only.]
          sabg_soil_bandloop(p) = sabg_soil(p)

          ! CASE1: No snow layers: all energy is absorbed in top soil layer
          if (snl(c) == 0) then
             sabg_lyr(p,:) = 0._r8
             sabg_lyr(p,1) = sabg(p)
             sabg_snl_sum  = sabg_lyr(p,1)
             ! [PORTED by Hui Tang: no-snow - NVP layer (index 0) absorbs before soil]
             ! fabd_nvp_col/fabi_nvp_col are Beer's law absorptance fractions (dimensionless,
             ! per unit flux incident on NVP). trd/tri are below-canopy direct/diffuse fluxes.
             ! sabg(p) is unchanged (ground total = NVP + soil via modified albedo).
             ! sabg_soil is corrected because it was computed using soil-only albedo (albsod).
             ! [PORTED by Hui Tang: nest the NVP guard — Fortran does not short-circuit .and.,
             !  and the NVP arrays are only allocated when use_nvp=.true.; combining the
             !  use_nvp check with array access in one .and. dereferences a null pointer.]
             if (use_nvp) then
             if (col%nvp_layer_active(patch%column(p))) then
                ! [PORTED by Hui Tang: snl==0 — moss is fully exposed (no snow). The Beer's-law
                !  absorption is the moss SURFACE solar -> store as sabg_nvp (analogous to sabg_soil),
                !  carve it out of the soil layer/sabg_soil, and set the internal sabg_lyr(p,0)=0
                !  (there is no snow above the moss, so no buried/SNICAR internal absorption).]
                sabg_nvp(p) = 0._r8
                do ib = 1, nband
                   sabg_nvp(p) = sabg_nvp(p) + &
                        surfalb_inst%fabd_nvp_col(c,ib) * trd(p,ib) + &
                        surfalb_inst%fabi_nvp_col(c,ib) * tri(p,ib)
                end do
                sabg_nvp(p)   = max(0._r8, min(sabg_nvp(p), sabg_lyr(p,1)))  ! per-column

                ![PORTED by Hui Tang: Soil patches receive the same amount of radiation, no matter it is under moss or not (exposed), per-area]
                sabg_lyr(p,1) = sabg_lyr(p,1) - sabg_nvp(p)
                sabg_soil(p)  = sabg_soil(p)  - sabg_nvp(p)

                ! [PORTED by Hui Tang (2026-06-12): keep sabg_lyr(p,0) = sabg_nvp so the SNICAR
                !  energy-conservation guard (sum(sabg_lyr)==sabg_snow) is satisfied without special-
                !  casing. This does NOT re-inject internal solar: the j=0 RHS internal term is
                !  frac_sno_eff*sabg_lyr_col(c,0) = 0 for snl==0 (no snow), so the moss solar still
                !  flows ONLY through hs_nvp (= frac_nvp_eff*sabg_nvp, thermostatted). The accounting
                !  uses frac_nvp_eff*sabg_nvp (not sabg_lyr(p,0)).]
                sabg_lyr(p,0) = sabg_nvp(p)
             end if
             end if  ! [PORTED by Hui Tang: close outer use_nvp guard]

             ! CASE 2: Snow layers present: absorbed radiation is scaled according to
             ! flux factors computed by SNICAR
          else
             do i = -nlevsno+1,1,1
                sabg_lyr(p,i) = flx_absdv(c,i)*trd(p,1) + flx_absdn(c,i)*trd(p,2) + &
                     flx_absiv(c,i)*tri(p,1) + flx_absin(c,i)*tri(p,2)
                ! summed radiation in active snow layers:
                if (i >= snl(c)+1) then
                   sabg_snl_sum = sabg_snl_sum + sabg_lyr(p,i)
                endif
                if (i > snl(c)+1) then ! if snow layer is below surface snow layer
                   !accumulate subsurface flux as a diagnostic for history file
                   sub_surf_abs_SW(p) = sub_surf_abs_SW(p) + sabg_lyr(p,i)
                endif
             enddo

             ! Divide absorbed by total, to get fraction absorbed in subsurface
             if (sabg_snl_sum /= 0._r8) then
                sub_surf_abs_SW(p) = sub_surf_abs_SW(p)/sabg_snl_sum
             else
                sub_surf_abs_SW(p) = 0._r8
             endif

             ! Error handling: The situation below can occur when solar radiation is
             ! NOT computed every timestep.
             ! When the number of snow layers has changed in between computations of the
             ! absorbed solar energy in each layer, we must redistribute the absorbed energy
             ! to avoid physically unrealistic conditions. The assumptions made below are
             ! somewhat arbitrary, but this situation does not arise very frequently.
             ! This error handling is implemented to accomodate any value of the
             ! radiation frequency.
             ! change condition to match sabg_snow isntead of sabg
             if (abs(sabg_snl_sum-sabg_snow(p)) > 0.00001_r8) then
                if (snl(c) == 0) then
                   sabg_lyr(p,-nlevsno+1:0) = 0._r8
                   sabg_lyr(p,1) = sabg(p)
                elseif (snl(c) == -1) then
                   sabg_lyr(p,-nlevsno+1:-1) = 0._r8
                   sabg_lyr(p,0) = sabg_snow(p)*0.6_r8
                   sabg_lyr(p,1) = sabg_snow(p)*0.4_r8
                else
                   sabg_lyr(p,:) = 0._r8
                   sabg_lyr(p,snl(c)+1) = sabg_snow(p)*0.75_r8
                   sabg_lyr(p,snl(c)+2) = sabg_snow(p)*0.25_r8
                endif
             endif

             ! If shallow snow depth, all solar radiation absorbed in top or top two snow layers
             ! to prevent unrealistic timestep soil warming 
             if (.not. use_subgrid_fluxes .or. lun%itype(l) == istdlak) then 
                if (snow_depth(c) < 0.10_r8) then
                   if (snl(c) == 0) then
                      sabg_lyr(p,-nlevsno+1:0) = 0._r8
                      sabg_lyr(p,1) = sabg(p)
                   elseif (snl(c) == -1) then
                      sabg_lyr(p,-nlevsno+1:-1) = 0._r8
                      sabg_lyr(p,0) = sabg(p)
                      sabg_lyr(p,1) = 0._r8
                   else
                      sabg_lyr(p,:) = 0._r8
                      sabg_lyr(p,snl(c)+1) = sabg(p)*0.75_r8
                      sabg_lyr(p,snl(c)+2) = sabg(p)*0.25_r8
                   endif
                endif
             endif
             ! [PORTED by Hui Tang: partial-snow NVP blend — CASE2 (snl<0) with partially-exposed moss]
             ! When snow is partial (frac_sno_eff < frac_nvp), fraction f_exp of the NVP is exposed
             ! and receives unattenuated radiation via Beer's law (per column area, same formula as CASE1
             ! weighted by f_exp).  The buried fraction (1-f_exp) receives SNICAR-attenuated radiation
             ! that was set by the SNICAR loop above (per unit snow area); multiplying by frac_sno_eff
             ! converts it to per column area.
             ! Combined: sabg_lyr(p,0) = f_exp*beer_per_col + (1-f_exp)*frac_sno_eff*snicar_per_snow
             ! [PORTED by Hui Tang (2026-06-11): SPLIT the moss solar into two un-weighted quantities,
             !  exactly mirroring the soil pair (sabg_soil vs sabg_lyr(p,1)):
             !    sabg_nvp(p)   = Beer's-law absorption = EXPOSED-moss SURFACE solar (full). The
             !                    frac_nvp_eff exposure weighting is applied later in the thermal solve
             !                    (hs_nvp, via nvp_exp*wtcol) — so NO f_exp here (it would double-count).
             !    sabg_lyr(p,0) = SNICAR moss-layer absorption = BURIED-moss INTERNAL solar (left as set
             !                    by the SNICAR loop above). The fse weighting is applied in the solve
             !                    (fse*sabg_lyr_col(c,0)) — so NO fse pre-weighting here.
             !  This replaces the old blend sabg_lyr(p,0)=f_exp*beer+(1-f_exp)*fse*snicar, which both
             !  double-weighted the fractions AND put the exposed solar internally (no -dhsdT surface
             !  thermostat) -> moss overheating.]
             if (use_nvp) then
             if (col%nvp_layer_active(c)) then
                if (col%frac_nvp(c) > 0._r8) then
                   sabg_nvp_beer = 0._r8
                   do ib = 1, nband
                      sabg_nvp_beer = sabg_nvp_beer + &
                           surfalb_inst%fabd_nvp_col(c,ib) * trd(p,ib) + &
                           surfalb_inst%fabi_nvp_col(c,ib) * tri(p,ib)
                   end do
                   sabg_nvp(p)   = max(0._r8, min(sabg_nvp_beer, sabg(p)-sabg_snow(p)))  ! per-column

                  ! [PORTED by Hui Tang: snow - NVP layer-0 SNICAR]
                  ! When use_nvp and SNICAR NVP layer-0 is active, flx_absdv(c,0)/flx_absiv(c,0)
                  ! already hold NVP absorption (set by SNICAR_RT above), so sabg_lyr(p,0) is
                  ! correct from the SNICAR loop above. Correct sabg_soil to use SNICAR soil layer.
                  ! [PORTED by Hui Tang: nest the NVP guard — see line ~768 for rationale]
                  ! sabg_lyr(p,1) = SNICAR soil-layer absorption (excludes NVP); use it directly.
                   frac_nvp_eff_loc = min(1._r8 - frac_sno(c), max(0._r8, col%frac_nvp(c) - frac_sno(c)))

                   ! [PORTED by Hui Tang (2026-06-13): guard the exposed-soil back-out against full snow
                   !  cover. When frac_sno_eff==1 the denominator (1-frac_sno_eff)=0 -> sabg_soil=Inf/NaN,
                   !  which then contaminates sw_grnd via the 0*NaN trap in SoilFluxesMod (the zero soil
                   !  weight does NOT cancel a NaN). At full snow cover there is no exposed soil, so
                   !  sabg_soil must be a finite 0.]
                   if (frac_sno(c) < 1._r8) then
                      sabg_soil(p) = (sabg(p) - sabg_snow(p)*frac_sno(c) - (frac_nvp_eff_loc/col%frac_nvp(c))*sabg_nvp(p))/(1-frac_sno(c))
                   else
                      sabg_soil(p) = 0._r8
                   end if

                   ! [PORTED by Hui Tang (2026-06-13): Phase-4 diagnostic (snl<0 partial snow). Dumps the
                   !  raw ground-solar pieces so we can measure how sabg(p) (albgrd total, now incl. the
                   !  exposed moss via Phase 3) decomposes vs the FGR reconstruction and the SABG tile.
                   !  Offline: M_alb = sabg(p) - fse*sabg_snow - (1-fse)*sabg_soil_bandloop (opaque moss);
                   !  M_beer = nvp_exp*sabg_nvp; FGR_solar ≈ (1-fse)*sabg_lyr1 + fse*sabg_snow + M_beer;
                   !  SABG_tile = fse*sabg_snow + (1-fse)*bandloop + M_beer. Residual to close: sabg(p)-FGR
                   !  and the soil seam (1-fse)*(bandloop - sabg_lyr1). REMOVE after Phase 4 verified.]
                   if (frac_sno_eff(c) > 0._r8 .and. frac_sno_eff(c) < col%frac_nvp(c)) then
                      write(iulog,*) '[NVP P4 RAD] nstep,c,p=', get_nstep(), c, p, &
                           ' sabg=', sabg(p), ' sabg_snow=', sabg_snow(p), &
                           ' bandloop=', sabg_soil_bandloop(p), ' sabg_lyr0=', sabg_lyr(p,0), &
                           ' sabg_lyr1=', sabg_lyr(p,1), ' sabg_nvp=', sabg_nvp(p), &
                           ' fse=', frac_sno_eff(c), ' fsno=', frac_sno(c), &
                           ' frac_nvp=', col%frac_nvp(c), &
                           ' nvp_exp=', min(1._r8-frac_sno_eff(c), &
                                            max(0._r8,col%frac_nvp(c)-frac_sno_eff(c)))/col%frac_nvp(c), &
                           ' FGRrecon=', (1._r8-frac_sno_eff(c))*sabg_lyr(p,1) + frac_sno_eff(c)*sabg_snow(p) &
                                + (min(1._r8-frac_sno_eff(c), &
                                       max(0._r8,col%frac_nvp(c)-frac_sno_eff(c)))/col%frac_nvp(c))*sabg_nvp(p)
                   end if
                end if
             end if  ! [PORTED by Hui Tang: close nvp_layer_active guard]
             end if  ! [PORTED by Hui Tang: close use_nvp guard]
          endif

          ! This situation should not happen:
          ! [PORTED by Hui Tang: skip endrun for NVP partial-snow — after the Beer's-law blend above,
          !  sabg_lyr(p,0) is per column area, so sum(sabg_lyr) intentionally exceeds sabg_snow
          !  (per snow area) when snl<0 and frac_sno_eff < frac_nvp]
          ! [PORTED by Hui Tang (2026-06-12): sabg_lyr(p,0)=sabg_nvp for snl==0 (set above), so the
          !  moss is in sum(sabg_lyr) and this check passes unchanged — no NVP special-casing needed.]
          sabg_sum_chk = sum(sabg_lyr(p,:))
          if (abs(sabg_sum_chk-sabg_snow(p)) > 0.00001_r8 .and. &
               .not. (use_nvp .and. col%nvp_layer_active(c) .and. &
                      snl(c) < 0 .and. frac_sno_eff(c) < col%frac_nvp(c))) then
             write(iulog,*)"SNICAR ERROR: Absorbed ground radiation not equal to summed snow layer radiation"
             write(iulog,*)"Diff        = ",sabg_sum_chk-sabg_snow(p)
             write(iulog,*)"sabg_snow(p)= ",sabg_snow(p)
             write(iulog,*)"sabg_sum(p) = ",sabg_sum_chk
             write(iulog,*)"snl(c)      = ",snl(c)
             write(iulog,*)"flx_absdv1  = ",trd(p,1)*(1.-albgrd(c,1))
             write(iulog,*)"flx_absdv2  = ",sum(flx_absdv(c,:))*trd(p,1)
             write(iulog,*)"flx_absiv1  = ",tri(p,1)*(1.-albgri(c,1))
             write(iulog,*)"flx_absiv2  = ",sum(flx_absiv(c,:))*tri(p,1)
             write(iulog,*)"flx_absdn1  = ",trd(p,2)*(1.-albgrd(c,2))
             write(iulog,*)"flx_absdn2  = ",sum(flx_absdn(c,:))*trd(p,2)
             write(iulog,*)"flx_absin1  = ",tri(p,2)*(1.-albgri(c,2))
             write(iulog,*)"flx_absin2  = ",sum(flx_absin(c,:))*tri(p,2)
             write(iulog,*)"albgrd_nir  = ",albgrd(c,2)
             write(iulog,*)"coszen      = ",coszen(c)
             call endrun(subgrid_index=c, subgrid_level=subgrid_level_column, msg=errmsg(sourcefile, __LINE__))
          endif

          ! Diagnostic: shortwave penetrating ground (e.g. top layer)
          if (lun%itype(l) == istsoil .or. lun%itype(l) == istcrop) then
             sabg_pen(p) = sabg(p) - sabg_lyr(p, snl(c)+1)
          end if

          if (use_snicar_frc) then

             ! BC aerosol forcing (patch-level):
             sfc_frc_bc(p) = sabg(p) - sabg_bc(p)

             ! OC aerosol forcing (patch-level):
             if (do_sno_oc) then
                sfc_frc_oc(p) = sabg(p) - sabg_oc(p)
             else
                sfc_frc_oc(p) = 0._r8
             endif

             ! dust aerosol forcing (patch-level):
             sfc_frc_dst(p) = sabg(p) - sabg_dst(p)

             ! all-aerosol forcing (patch-level):
             sfc_frc_aer(p) = sabg(p) - sabg_pur(p)

             ! forcings averaged only over snow:
             if (frac_sno(c) > 0._r8) then
                sfc_frc_bc_sno(p)  = sfc_frc_bc(p)/frac_sno(c)
                sfc_frc_oc_sno(p)  = sfc_frc_oc(p)/frac_sno(c)
                sfc_frc_dst_sno(p) = sfc_frc_dst(p)/frac_sno(c)
                sfc_frc_aer_sno(p) = sfc_frc_aer(p)/frac_sno(c)
             else
                sfc_frc_bc_sno(p)  = spval
                sfc_frc_oc_sno(p)  = spval
                sfc_frc_dst_sno(p) = spval
                sfc_frc_aer_sno(p) = spval
             endif
          end if
       enddo

       ! Radiation diagnostics

       do fp = 1,num_nourbanp
          p = filter_nourbanp(fp)
          g = patch%gridcell(p)
          c = patch%column(p)

          ! NDVI and reflected solar radiation

          rvis = albd(p,1)*forc_solad_col(c,1) + albi(p,1)*forc_solai(g,1)
          rnir = albd(p,2)*forc_solad_col(c,2) + albi(p,2)*forc_solai(g,2)
          fsr(p) = rvis + rnir
          if (use_SSRE) then
             rvisSF = albdSF(p,1)*forc_solad_col(c,1) + albiSF(p,1)*forc_solai(g,1)
             rnirSF = albdSF(p,2)*forc_solad_col(c,2) + albiSF(p,2)*forc_solai(g,2)
             fsrSF(p) = rvisSF + rnirSF
             ssre_fsr(p) = fsr(p)-fsrSF(p)
          end if
          fsds_vis_d(p) = forc_solad_col(c,1)
          fsds_nir_d(p) = forc_solad_col(c,2)
          fsds_vis_i(p) = forc_solai(g,1)
          fsds_nir_i(p) = forc_solai(g,2)
          fsr_vis_d(p)  = albd(p,1)*forc_solad_col(c,1)
          fsr_nir_d(p)  = albd(p,2)*forc_solad_col(c,2)
          fsr_vis_i(p)  = albi(p,1)*forc_solai(g,1)
          fsr_nir_i(p)  = albi(p,2)*forc_solai(g,2)
          if (use_SSRE) then
             fsrSF_vis_d(p)  = albdSF(p,1)*forc_solad_col(c,1)
             fsrSF_nir_d(p)  = albdSF(p,2)*forc_solad_col(c,2)
             fsrSF_vis_i(p)  = albiSF(p,1)*forc_solai(g,1)
             fsrSF_nir_i(p)  = albiSF(p,2)*forc_solai(g,2)

             ssre_fsr_vis_d(p) = fsrSF_vis_d(p)-fsr_vis_d(p)
             ssre_fsr_nir_d(p) = fsrSF_nir_d(p)-fsr_nir_d(p)
             ssre_fsr_vis_i(p) = fsrSF_vis_i(p)-fsr_vis_i(p)
             ssre_fsr_nir_i(p) = fsrSF_nir_i(p)-fsr_nir_i(p)
          end if
          if ( is_near_local_noon( grc%londeg(g), deltasec=nint(dtime)/2 ) )then
             fsds_vis_d_ln(p) = forc_solad_col(c,1)
             fsds_nir_d_ln(p) = forc_solad_col(c,2)
             fsr_vis_d_ln(p) = albd(p,1)*forc_solad_col(c,1)
             fsr_nir_d_ln(p) = albd(p,2)*forc_solad_col(c,2)
             fsds_vis_i_ln(p) = forc_solai(g,1)
             parveg_ln(p)     = parveg(p)
          else
             fsds_vis_d_ln(p) = spval
             fsds_nir_d_ln(p) = spval
             fsr_vis_d_ln(p) = spval
             fsr_nir_d_ln(p) = spval
             fsds_vis_i_ln(p) = spval
             parveg_ln(p)     = spval
          end if
          if (use_SSRE) then
             if ( is_near_local_noon( grc%londeg(g), deltasec=nint(dtime)/2 ) )then
                fsrSF_vis_d_ln(p) = albdSF(p,1)*forc_solad_col(c,1)
                fsrSF_nir_d_ln(p) = albdSF(p,2)*forc_solad_col(c,2)
             else
                fsrSF_vis_d_ln(p) = spval
                fsrSF_nir_d_ln(p) = spval
             end if
          end if
          ! diagnostic variables (downwelling and absorbed radiation partitioning) for history files
          ! (OPTIONAL)
          c = patch%column(p)
          if (snl(c) < 0) then
             fsds_sno_vd(p) = forc_solad_col(c,1)
             fsds_sno_nd(p) = forc_solad_col(c,2)
             fsds_sno_vi(p) = forc_solai(g,1)
             fsds_sno_ni(p) = forc_solai(g,2)

             fsr_sno_vd(p) = fsds_vis_d(p)*albsnd_hst(c,1)
             fsr_sno_nd(p) = fsds_nir_d(p)*albsnd_hst(c,2)
             fsr_sno_vi(p) = fsds_vis_i(p)*albsni_hst(c,1)
             fsr_sno_ni(p) = fsds_nir_i(p)*albsni_hst(c,2)
          else
             fsds_sno_vd(p) = spval
             fsds_sno_nd(p) = spval
             fsds_sno_vi(p) = spval
             fsds_sno_ni(p) = spval

             fsr_sno_vd(p) = spval
             fsr_sno_nd(p) = spval
             fsr_sno_vi(p) = spval
             fsr_sno_ni(p) = spval
          endif
       end do

       ! TODO: urban snow-free albedos:
       do fp = 1,num_urbanp
          p = filter_urbanp(fp)
          g = patch%gridcell(p)
          c = patch%column(p)

          if(elai(p)==0.0_r8.and.fabd(p,1)>0._r8)then
             if ( local_debug ) write(iulog,*) 'absorption without LAI',elai(p),tlai(p),fabd(p,1),p
          endif

          ! Solar incident

          fsds_vis_d(p) = forc_solad_col(c,1)
          fsds_nir_d(p) = forc_solad_col(c,2)
          fsds_vis_i(p) = forc_solai(g,1)
          fsds_nir_i(p) = forc_solai(g,2)

          ! Determine local noon incident solar
          if ( is_near_local_noon( grc%londeg(g), deltasec=nint(dtime)/2 ) )then
             fsds_vis_d_ln(p) = forc_solad_col(c,1)
             fsds_nir_d_ln(p) = forc_solad_col(c,2)
             fsds_vis_i_ln(p) = forc_solai(g,1)
             parveg_ln(p)     = 0._r8
          else
             fsds_vis_d_ln(p) = spval
             fsds_nir_d_ln(p) = spval
             fsds_vis_i_ln(p) = spval
             parveg_ln(p)     = spval
          endif

          ! Solar reflected
          ! per unit ground area (roof, road) and per unit wall area (sunwall, shadewall)

          fsr_vis_d(p) = albd(p,1) * forc_solad_col(c,1)
          fsr_nir_d(p) = albd(p,2) * forc_solad_col(c,2)
          fsr_vis_i(p) = albi(p,1) * forc_solai(g,1)
          fsr_nir_i(p) = albi(p,2) * forc_solai(g,2)

          ! Determine local noon reflected solar
          if ( is_near_local_noon( grc%londeg(g), deltasec=nint(dtime)/2 ) )then
             fsr_vis_d_ln(p) = fsr_vis_d(p)
             fsr_nir_d_ln(p) = fsr_nir_d(p)
          else
             fsr_vis_d_ln(p) = spval
             fsr_nir_d_ln(p) = spval
          endif
          fsr(p) = fsr_vis_d(p) + fsr_nir_d(p) + fsr_vis_i(p) + fsr_nir_i(p)
       end do

     end associate

   end subroutine SurfaceRadiation

end module SurfaceRadiationMod
