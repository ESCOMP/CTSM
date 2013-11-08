module initSurfalbMod
  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Computes initial surface albedo calculation - 
  ! Initialization of ecosystem dynamics is needed for this
  !
  ! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use abortutils,   only : endrun
  use clm_varctl,   only : iulog, use_cn
  !
  ! !PUBLIC TYPES:
  implicit none
  logical, public :: do_initsurfalb
  ! save
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: InitSurfAlb
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine initSurfalb( calday, declin )
    !
    ! !DESCRIPTION:
    ! The variable, h2osoi_vol, is needed by the soil albedo routine - this is not needed
    ! on restart since it is computed before the soil albedo computation is called.
    ! The remaining variables are initialized by calls to ecosystem dynamics and
    ! albedo subroutines.
    !
    ! !USES:
    use shr_orb_mod         , only : shr_orb_decl
    use shr_const_mod       , only : SHR_CONST_PI
    use clmtype
    use spmdMod             , only : masterproc,iam
    use decompMod           , only : get_proc_clumps, get_clump_bounds, get_proc_bounds, bounds_type
    use filterMod           , only : filter, filter_inactive_and_active
    use clm_varpar          , only : nlevsoi, nlevsno, nlevlak, nlevgrnd
    use clm_varcon          , only : zlnd, istsoil, denice, denh2o, &
                                     icol_roof, icol_road_imperv, icol_road_perv
    use clm_varcon          , only : istcrop
    use clm_time_manager    , only : get_step_size
    use FracWetMod          , only : FracWet
    use SurfaceAlbedoMod    , only : SurfaceAlbedo
    use CNEcosystemDynMod   , only : CNEcosystemDynNoLeaching, CNEcosystemDynLeaching
    use CNVegStructUpdateMod, only : CNVegStructUpdate
    use STATICEcosysDynMod  , only : EcosystemDyn, interpMonthlyVeg
    use UrbanMod            , only : UrbanAlbedo
    use abortutils          , only : endrun
    !
    ! !ARGUMENTS:
    implicit none
    real(r8), intent(in) :: calday               ! calendar day for declin
    real(r8), intent(in) :: declin               ! declination angle (radians) for calday
    !
    ! !LOCAL VARIABLES:
    integer :: nc,j,l,c,p,fc ! indices
    integer :: nclumps       ! number of clumps on this processor
    integer :: ier           ! MPI return code
    real(r8):: snowbd        ! temporary calculation of snow bulk density (kg/m3)
    real(r8):: fmelt         ! snowbd/100
    type(bounds_type) :: bounds  ! bounds
    !-----------------------------------------------------------------------

   associate(& 
   dz                    =>    cps%dz                  , & ! Input:  [real(r8) (:,:)]  layer thickness depth (m)                       
   h2osoi_ice            =>    cws%h2osoi_ice          , & ! Input:  [real(r8) (:,:)]  ice lens (kg/m2)                                
   h2osoi_liq            =>    cws%h2osoi_liq          , & ! Input:  [real(r8) (:,:)]  liquid water (kg/m2)                            
   h2osoi_vol            =>    cws%h2osoi_vol          , & ! Output: [real(r8) (:,:)]  volumetric soil water (0<=h2osoi_vol<=watsat) [m3/m3]
   snow_depth            =>    cps%snow_depth          , & ! Output: [real(r8) (:)]  snow height (m)                                   
   h2osno                =>    cws%h2osno              , & ! Input:  [real(r8) (:)]  snow water (mm H2O)                               
   frac_sno              =>    cps%frac_sno            , & ! Output: [real(r8) (:)]  fraction of ground covered by snow (0 to 1)       
   soilpsi               =>    cps%soilpsi             , & ! Input:  [real(r8) (:,:)]  soil water potential in each soil layer (MPa)   
   frac_veg_nosno_alb    =>    pps%frac_veg_nosno_alb  , & ! Input:  [integer (:)]  fraction of vegetation not covered by snow (0 OR 1) [-]
   frac_veg_nosno        =>    pps%frac_veg_nosno      , & ! Output: [integer (:)]  fraction of vegetation not covered by snow (0 OR 1) [-]
   fwet                  =>    pps%fwet                , & ! Output: [real(r8) (:)]  fraction of canopy that is wet (0 to 1) (pft-level)
   htop                  =>    pps%htop                , & ! Output: [real(r8) (:)]  canopy top (m)                                    
   hbot                  =>    pps%hbot                , & ! Output: [real(r8) (:)]  canopy bottom (m)                                 
   tlai                  =>    pps%tlai                , & ! Output: [real(r8) (:)]  one-sided leaf area index, no burying by snow     
   tsai                  =>    pps%tsai                , & ! Output: [real(r8) (:)]  one-sided stem area index, no burying by snow     
   elai                  =>    pps%elai                , & ! Output: [real(r8) (:)]  one-sided leaf area index with burying by snow    
   esai                  =>    pps%esai                , & ! Output: [real(r8) (:)]  one-sided stem area index with burying by snow    
   fdry                  =>    pps%fdry                  & ! Output: [real(r8) (:)]  fraction of foliage that is green and dry [-] (new)
   )

    ! ========================================================================
    ! Determine surface albedo - initialized by calls to ecosystem dynamics and
    ! albedo subroutines. Note: elai, esai, frac_veg_nosno_alb are computed in
    ! Ecosysdyn and needed by routines FracWet and SurfaceAlbedo and 
    ! frac_veg_nosno is needed by FracWet
    ! fwet is needed in routine TwoStream (called by SurfaceAlbedo)
    ! frac_sno is needed by SoilAlbedo (called by SurfaceAlbedo)
    ! ========================================================================

    if (.not. use_cn) then
       ! the default mode uses prescribed vegetation structure
       ! Read monthly vegetation data for interpolation to daily values
       
       call get_proc_bounds(bounds)
       call interpMonthlyVeg(bounds)
    end if

    ! Determine clump bounds for this processor

    nclumps = get_proc_clumps()

    ! Loop over clumps on this processor
    !$OMP PARALLEL DO PRIVATE (nc,p,j,l,c,fc,bounds,snowbd,fmelt)
    do nc = 1,nclumps

       ! Determine clump bounds

       call get_clump_bounds(nc, bounds)

       ! Determine variables needed by SurfaceAlbedo for lake points

       do p = bounds%begp,bounds%endp
          l = pft%landunit(p)
          if (lun%lakpoi(l)) then
             fwet(p) = 0._r8
             fdry(p) = 0._r8
             elai(p) = 0._r8
             esai(p) = 0._r8
             htop(p) = 0._r8
             hbot(p) = 0._r8
             tlai(p) = 0._r8
             tsai(p) = 0._r8
             frac_veg_nosno_alb(p) = 0._r8
             frac_veg_nosno(p) = 0._r8
          end if
       end do

       ! ============================================================================
       ! Ecosystem dynamics: Uses CN, or static parameterizations
       ! ============================================================================

       if (use_cn) then
          do j = 1, nlevgrnd
             do fc = 1, filter(nc)%num_soilc
                c = filter(nc)%soilc(fc)
                soilpsi(c,j) = -15.0_r8
             end do
          end do
       end if

       ! Determine variables needed for SurfaceAlbedo for non-lake points

       do c = bounds%begc, bounds%endc
          l = col%landunit(c)
          if (lun%urbpoi(l)) then
             ! From Bonan 1996 (LSM technical note)
             frac_sno(c) = min( snow_depth(c)/0.05_r8, 1._r8)
          else
              frac_sno(c) = 0._r8
             ! snow cover fraction as in Niu and Yang 2007
             if(snow_depth(c) .gt. 0.0)  then
                snowbd   = min(400._r8,h2osno(c)/snow_depth(c)) !bulk density of snow (kg/m3)
                fmelt    = (snowbd/100.)**1.
                ! 100 is the assumed fresh snow density; 1 is a melting factor that could be
                ! reconsidered, optimal value of 1.5 in Niu et al., 2007
                frac_sno(c) = tanh( snow_depth(c) /(2.5 * zlnd * fmelt) )
             endif
          end if
       end do

       if (use_cn) then
          ! CN initialization is done only on the soil landunits.
          
          call CNEcosystemDynNoLeaching(bounds, &
               filter(nc)%num_soilc, filter(nc)%soilc, &
               filter(nc)%num_soilp, filter(nc)%soilp, &
               filter(nc)%num_pcropp, filter(nc)%pcropp, doalb=.true.)
          
          call CNEcosystemDynLeaching(bounds, &
               filter(nc)%num_soilc, filter(nc)%soilc, &
               filter(nc)%num_soilp, filter(nc)%soilp, &
               filter(nc)%num_pcropp, filter(nc)%pcropp, doalb=.true.)
          
       else
          ! this is the default call if CN not set
          
          call EcosystemDyn(bounds, &
               filter(nc)%num_nolakep, filter(nc)%nolakep, doalb=.true.)
       end if

       do p = bounds%begp, bounds%endp
          l = pft%landunit(p)
          if (.not. lun%lakpoi(l)) then
             frac_veg_nosno(p) = frac_veg_nosno_alb(p)
             fwet(p) = 0._r8
          end if
       end do
       
       call FracWet(filter(nc)%num_nolakep, filter(nc)%nolakep)
       
       ! Compute Surface Albedo - all land points (including lake) other than urban
       ! Needs as input fracion of soil covered by snow (Z.-L. Yang U. Texas)

       ! Note that in both the SurfaceAlbedo call and the UrbanAlbedo call, we use the
       ! version of the filters that includes inactive as well as active points, because
       ! some variables computed there are needed over points that later become active
       ! due to landuse change.

       call SurfaceAlbedo(bounds, &
            filter_inactive_and_active(nc)%num_nourbanc, &
            filter_inactive_and_active(nc)%nourbanc, &
            filter_inactive_and_active(nc)%num_nourbanp, &
            filter_inactive_and_active(nc)%nourbanp, &
            calday, declin)
       

       ! Determine albedos for urban landunits

       if (filter_inactive_and_active(nc)%num_urbanl > 0) then

          call UrbanAlbedo(bounds, &
               filter_inactive_and_active(nc)%num_urbanl, &
               filter_inactive_and_active(nc)%urbanl, &
               filter_inactive_and_active(nc)%num_urbanc, &
               filter_inactive_and_active(nc)%urbanc, &
               filter_inactive_and_active(nc)%num_urbanp, &
               filter_inactive_and_active(nc)%urbanp )

       end if

    end do   ! end of loop over clumps
    !$OMP END PARALLEL DO

    end associate 
   end subroutine initSurfalb

end module initSurfalbMod
