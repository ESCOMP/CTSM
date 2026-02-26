module DustEmisZender2003

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Routines in this module calculate Dust mobilization and dry deposition for dust.
  ! Simulates dust mobilization due to wind from the surface into the
  ! lowest atmospheric layer. On output flx_mss_vrt_dst(ndst) is the surface dust
  ! emission (kg/m**2/s) [ + = to atm].
  ! Calculates the turbulent component of dust dry deposition, (the turbulent deposition
  ! velocity through the lowest atmospheric layer). CAM will calculate the settling
  ! velocity through the whole atmospheric column. The two calculations will determine
  ! the dust dry deposition flux to the surface.
  !
  ! !USES:
  use shr_kind_mod         , only : r8 => shr_kind_r8
  use shr_log_mod          , only : errMsg => shr_log_errMsg
  use shr_infnan_mod       , only : nan => shr_infnan_nan, assignment(=)
  use clm_varpar           , only : dst_src_nbr, ndst
  use clm_varcon           , only : grav, spval
  use landunit_varcon      , only : istcrop, istsoil
  use clm_varctl           , only : iulog
  use abortutils           , only : endrun
  use decompMod            , only : bounds_type, subgrid_level_landunit
  use atm2lndType          , only : atm2lnd_type
  use SoilStateType        , only : soilstate_type
  use CanopyStateType      , only : canopystate_type
  use WaterStateBulkType       , only : waterstatebulk_type
  use WaterDiagnosticBulkType       , only : waterdiagnosticbulk_type
  use FrictionVelocityMod  , only : frictionvel_type
  use LandunitType         , only : lun
  use PatchType            , only : patch
  use ZenderSoilErodStreamType,  only : soil_erod_stream_type
  use DustEmisBase         , only : dust_emis_base_type
  !
  ! !PUBLIC TYPES
  implicit none
  private
  !
  ! !PRIVATE DATA:
  !
  !
  ! !PUBLIC DATA TYPES:
  !
  type, public, extends(dust_emis_base_type) :: dust_emis_zender2003_type

     real(r8), pointer, private :: mbl_bsn_fct_col           (:)   ! [dimensionless] basin factor, or soil rodibility, time-constant
     type(soil_erod_stream_type), private :: soil_erod_stream      ! Zender soil erodibility stream data

   contains

     procedure , public  :: Init => InitZender2003
     procedure , public  :: DustEmission    ! Dust mobilization
     procedure , public  :: Clean => CleanZender2003
     procedure , private :: InitAllocate    ! Allocate data
     procedure , private :: InitHistory     ! History initialization
     procedure , private :: InitCold

  end type dust_emis_zender2003_type

  interface dust_emis_zender2003_type
     ! initialize a new dust emission object
      module procedure constructor
  end interface dust_emis_zender2003_type
  !------------------------------------------------------------------------

  character(len=*), parameter, private :: sourcefile = &
       __FILE__

contains

  !-----------------------------------------------------------------------
  type(dust_emis_zender2003_type) function constructor()
  !
  ! Creates a dust emission object for Zender-2003 type
  ! For now this is just a placeholder
  !-----------------------------------------------------------------------

  end function constructor

  !------------------------------------------------------------------------

  subroutine InitZender2003(this, bounds, NLFilename)

   ! Initialization for this extended class, calling base level initiation and adding to it
    class(dust_emis_zender2003_type) :: this
    type(bounds_type), intent(in) :: bounds
    character(len=*),  intent(in) :: NLFilename

    call this%soil_erod_stream%Init( bounds, NLFilename )
    call this%InitBase(bounds, NLFilename)
    call this%InitAllocate (bounds)
    call this%InitHistory  (bounds)
    call this%InitCold     (bounds)

  end subroutine InitZender2003

  !------------------------------------------------------------------------

  subroutine InitAllocate(this, bounds)
   !
   ! !ARGUMENTS:
   class (dust_emis_zender2003_type) :: this
   type(bounds_type), intent(in) :: bounds
   !
   ! !LOCAL VARIABLES:
   integer :: begc,endc
   !------------------------------------------------------------------------

   begc = bounds%begc ; endc = bounds%endc

   allocate(this%mbl_bsn_fct_col           (begc:endc))        ; this%mbl_bsn_fct_col           (:)   = nan

  end subroutine InitAllocate

  !------------------------------------------------------------------------

  subroutine CleanZender2003(this)
    !
    ! Deallocation for this extended class, calling base level deallocation and adding to it
    ! !ARGUMENTS:
    class (dust_emis_zender2003_type) :: this
    !
    ! !LOCAL VARIABLES:
    !------------------------------------------------------------------------

    call this%CleanBase()
    deallocate(this%mbl_bsn_fct_col)

  end subroutine CleanZender2003

  !------------------------------------------------------------------------

  subroutine InitHistory(this, bounds)
    !
    ! !USES:
    use histFileMod, only : hist_addfld1d
    !
    !
    ! !ARGUMENTS:
    class (dust_emis_zender2003_type) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: begc,endc
    !------------------------------------------------------------------------

    begc = bounds%begc; endc = bounds%endc

    if ( this%soil_erod_stream%UseStreams() )then
       this%mbl_bsn_fct_col(begc:endc) = spval
       call hist_addfld1d (fname='LND_MBL', units='fraction',  &
               avgflag='A', long_name='Soil erodibility factor', &
               ptr_col=this%mbl_bsn_fct_col, default='inactive')
    end if

  end subroutine InitHistory

  !-----------------------------------------------------------------------

  subroutine InitCold(this, bounds)
    !
    ! Initialize values from a cold start
    ! !ARGUMENTS:
    class (dust_emis_zender2003_type) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: c,l
    !-----------------------------------------------------------------------

    if ( this%soil_erod_stream%UseStreams() )then
       call this%soil_erod_stream%CalcDustSource( bounds, &
                               this%mbl_bsn_fct_col(bounds%begc:bounds%endc) )
    else
       this%mbl_bsn_fct_col(:) = 1.0_r8
    end if

  end subroutine InitCold

  !------------------------------------------------------------------------

  subroutine DustEmission (this, bounds, &
       num_nolakep, filter_nolakep, &
       atm2lnd_inst, soilstate_inst, canopystate_inst, waterstatebulk_inst, waterdiagnosticbulk_inst, &
       frictionvel_inst)
    !
    ! !DESCRIPTION:
    ! Dust mobilization. This code simulates dust mobilization due to wind
    ! from the surface into the lowest atmospheric layer
    ! On output flx_mss_vrt_dst(ndst) is the surface dust emission
    ! (kg/m**2/s) [ + = to atm]
    ! Source: C. Zender's dust model
    !
    ! !USES
    use shr_const_mod, only : SHR_CONST_RHOFW
    use subgridaveMod, only : p2g
    !
    ! !ARGUMENTS:
    class (dust_emis_zender2003_type)      :: this
    type(bounds_type)      , intent(in)    :: bounds
    integer                , intent(in)    :: num_nolakep                 ! number of column non-lake points in patch filter
    integer                , intent(in)    :: filter_nolakep(num_nolakep) ! patch filter for non-lake points
    type(atm2lnd_type)     , intent(in)    :: atm2lnd_inst
    type(soilstate_type)   , intent(in)    :: soilstate_inst
    type(canopystate_type) , intent(in)    :: canopystate_inst
    type(waterstatebulk_type)  , intent(in)    :: waterstatebulk_inst
    type(waterdiagnosticbulk_type)  , intent(in)    :: waterdiagnosticbulk_inst
    type(frictionvel_type) , intent(in)    :: frictionvel_inst

    !
    ! !LOCAL VARIABLES
    integer  :: fp,p,c,l,g,m,n      ! indices
    real(r8) :: liqfrac             ! fraction of total water that is liquid
    real(r8) :: wnd_frc_rat         ! [frc] Wind friction threshold over wind friction
    real(r8) :: wnd_frc_slt_dlt     ! [m s-1] Friction velocity increase from saltatn
    real(r8) :: wnd_rfr_dlt         ! [m s-1] Reference windspeed excess over threshld
    real(r8) :: dst_slt_flx_rat_ttl
    real(r8) :: flx_mss_hrz_slt_ttl
    real(r8) :: flx_mss_vrt_dst_ttl(bounds%begp:bounds%endp)
    real(r8) :: frc_thr_wet_fct
    real(r8) :: frc_thr_rgh_fct
    real(r8) :: wnd_frc_thr_slt
    real(r8) :: wnd_rfr_thr_slt
    real(r8) :: wnd_frc_slt
    real(r8) :: lnd_frc_mbl(bounds%begp:bounds%endp)
    real(r8) :: bd
    real(r8) :: gwc_sfc
    real(r8) :: ttlai(bounds%begp:bounds%endp)
    real(r8) :: tlai_lu(bounds%begl:bounds%endl)
    real(r8) :: sumwt(bounds%begl:bounds%endl) ! sum of weights
    logical  :: found                          ! temporary for error check
    integer  :: index
    !
    ! constants
    !
    real(r8), parameter :: cst_slt = 2.61_r8           ! [frc] Saltation constant
    real(r8), parameter :: flx_mss_fdg_fct = 5.0e-4_r8 ! [frc] Empir. mass flx tuning eflx_lh_vegt
    real(r8), parameter :: vai_mbl_thr = 0.3_r8        ! [m2 m-2] VAI threshold quenching dust mobilization
    character(len=*),parameter :: subname = 'DUSTEmission'
    !------------------------------------------------------------------------

    associate(                                                         &
         forc_rho            => atm2lnd_inst%forc_rho_downscaled_col , & ! Input:  [real(r8) (:)   ]  downscaled density (kg/m**3)

         gwc_thr             => soilstate_inst%gwc_thr_col           , & ! Input:  [real(r8) (:)   ]  threshold gravimetric soil moisture based on clay content
         mss_frc_cly_vld     => soilstate_inst%mss_frc_cly_vld_col   , & ! Input:  [real(r8) (:)   ]  [frc] Mass fraction clay limited to 0.20
         watsat              => soilstate_inst%watsat_col            , & ! Input:  [real(r8) (:,:) ]  saturated volumetric soil water

         tlai                => canopystate_inst%tlai_patch          , & ! Input:  [real(r8) (:)   ]  one-sided leaf area index, no burying by snow
         tsai                => canopystate_inst%tsai_patch          , & ! Input:  [real(r8) (:)   ]  one-sided stem area index, no burying by snow

         frac_sno            => waterdiagnosticbulk_inst%frac_sno_col         , & ! Input:  [real(r8) (:)   ]  fraction of ground covered by snow (0 to 1)
         h2osoi_vol          => waterstatebulk_inst%h2osoi_vol_col       , & ! Input:  [real(r8) (:,:) ]  volumetric soil water (0<=h2osoi_vol<=watsat)
         h2osoi_liq          => waterstatebulk_inst%h2osoi_liq_col       , & ! Input:  [real(r8) (:,:) ]  liquid soil water (kg/m2)
         h2osoi_ice          => waterstatebulk_inst%h2osoi_ice_col       , & ! Input:  [real(r8) (:,:) ]  frozen soil water (kg/m2)

         fv                  => frictionvel_inst%fv_patch            , & ! Input:  [real(r8) (:)   ]  friction velocity (m/s) (for dust model)
         u10                 => frictionvel_inst%u10_patch           , & ! Input:  [real(r8) (:)   ]  10-m wind (m/s) (created for dust model)

         mbl_bsn_fct         => this%mbl_bsn_fct_col                 , & ! Input:  [real(r8) (:)   ]  basin factor
         flx_mss_vrt_dst     => this%flx_mss_vrt_dst_patch           , & ! Output: [real(r8) (:,:) ]  surface dust emission (kg/m**2/s)
         flx_mss_vrt_dst_tot => this%flx_mss_vrt_dst_tot_patch         & ! Output: [real(r8) (:)   ]  total dust flux back to atmosphere (pft)
         )

      ttlai(bounds%begp : bounds%endp) = 0._r8
      ! make lai average at landunit level
      do fp = 1,num_nolakep
         p = filter_nolakep(fp)
         ttlai(p) = tlai(p)+tsai(p)
      enddo

      tlai_lu(bounds%begl : bounds%endl) = spval
      sumwt(bounds%begl : bounds%endl) = 0._r8
      do p = bounds%begp,bounds%endp
         if (ttlai(p) /= spval .and. patch%active(p) .and. patch%wtlunit(p) /= 0._r8) then
            c = patch%column(p)
            l = patch%landunit(p)
            if (sumwt(l) == 0._r8) tlai_lu(l) = 0._r8
            tlai_lu(l) = tlai_lu(l) + ttlai(p) * patch%wtlunit(p)
            sumwt(l) = sumwt(l) + patch%wtlunit(p)
         end if
      end do
      found = .false.
      do l = bounds%begl,bounds%endl
         if (sumwt(l) > 1.0_r8 + 1.e-6_r8) then
            found = .true.
            index = l
            exit
         else if (sumwt(l) /= 0._r8) then
            tlai_lu(l) = tlai_lu(l)/sumwt(l)
         end if
      end do
      if (found) then
         write(iulog,*) subname//':: error: sumwt is greater than 1.0 at l= ',index
         call endrun(subgrid_index=index, subgrid_level=subgrid_level_landunit, msg=errMsg(sourcefile, __LINE__))
      end if

      ! Loop through patches

      ! initialize variables which get passed to the atmosphere
      flx_mss_vrt_dst(bounds%begp:bounds%endp,:)=0._r8

      do fp = 1,num_nolakep
         p = filter_nolakep(fp)
         c = patch%column(p)
         l = patch%landunit(p)

         ! the following code from subr. lnd_frc_mbl_get was adapted for lsm use
         ! purpose: return fraction of each gridcell suitable for dust mobilization

         ! the "bare ground" fraction of the current sub-gridscale cell decreases
         ! linearly from 1 to 0 as VAI(=tlai+tsai) increases from 0 to vai_mbl_thr
         ! if ice sheet, wetland, or lake, no dust allowed

         if (lun%itype(l) == istsoil .or. lun%itype(l) == istcrop) then
            if (tlai_lu(l) < vai_mbl_thr) then
               lnd_frc_mbl(p) = 1.0_r8 - (tlai_lu(l))/vai_mbl_thr
            else
               lnd_frc_mbl(p) = 0.0_r8
            endif
            lnd_frc_mbl(p) = lnd_frc_mbl(p) * (1.0_r8 - frac_sno(c))
         else
            lnd_frc_mbl(p) = 0.0_r8
         end if
      end do

      do fp = 1,num_nolakep
         p = filter_nolakep(fp)
         if (lnd_frc_mbl(p)>1.0_r8 .or. lnd_frc_mbl(p)<0.0_r8) then
            write(iulog,*)'Error dstmbl: pft= ',p,' lnd_frc_mbl(p)= ',lnd_frc_mbl(p), &
                           errMsg(sourcefile, __LINE__)
            call endrun("Bad value for dust mobilization fraction")
            return
         end if
      end do

      ! reset history output variables before next if-statement to avoid output = inf

      do fp = 1,num_nolakep
         p = filter_nolakep(fp)
         flx_mss_vrt_dst_tot(p) = 0.0_r8
      end do
      do n = 1, ndst
         do fp = 1,num_nolakep
            p = filter_nolakep(fp)
            flx_mss_vrt_dst(p,n) = 0.0_r8
         end do
      end do

      do fp = 1,num_nolakep
         p = filter_nolakep(fp)
         c = patch%column(p)
         l = patch%landunit(p)
         g = patch%gridcell(p)

         ! only perform the following calculations if lnd_frc_mbl is non-zero

         if (lnd_frc_mbl(p) > 0.0_r8) then

            ! the following comes from subr. frc_thr_rgh_fct_get
            ! purpose: compute factor by which surface roughness increases threshold
            !          friction velocity (currently a constant)

            frc_thr_rgh_fct = 1.0_r8

            ! the following comes from subr. frc_thr_wet_fct_get
            ! purpose: compute factor by which soil moisture increases threshold friction velocity
            ! adjust threshold velocity for inhibition by moisture
            ! modified 4/5/2002 (slevis) to use gravimetric instead of volumetric
            ! water content

            bd = (1._r8-watsat(c,1))*2.7e3_r8      ![kg m-3] Bulk density of dry surface soil
            gwc_sfc = h2osoi_vol(c,1)*SHR_CONST_RHOFW/bd    ![kg kg-1] Gravimetric H2O cont
            if (gwc_sfc > gwc_thr(c)) then
               frc_thr_wet_fct = sqrt(1.0_r8 + 1.21_r8 * (100.0_r8*(gwc_sfc - gwc_thr(c)))**0.68_r8)
            else
               frc_thr_wet_fct = 1.0_r8
            end if

            ! slevis: adding liqfrac here, because related to effects from soil water

            liqfrac = max( 0.0_r8, min( 1.0_r8, h2osoi_liq(c,1) / (h2osoi_ice(c,1)+h2osoi_liq(c,1)+1.0e-6_r8) ) )

            ! the following lines come from subr. dst_mbl
            ! purpose: adjust threshold friction velocity to acct for moisture and
            !          roughness. The ratio saltation_factor / sqrt(forc_rho) comes from
            !          subr. wnd_frc_thr_slt_get which computes dry threshold
            !          friction velocity for saltation

            wnd_frc_thr_slt = this%saltation_factor / sqrt(forc_rho(c)) * frc_thr_wet_fct * frc_thr_rgh_fct

            ! reset these variables which will be updated in the following if-block

            wnd_frc_slt = fv(p)
            flx_mss_hrz_slt_ttl = 0.0_r8
            flx_mss_vrt_dst_ttl(p) = 0.0_r8

            ! the following line comes from subr. dst_mbl
            ! purpose: threshold saltation wind speed

            wnd_rfr_thr_slt = u10(p) * wnd_frc_thr_slt / fv(p)

            ! the following if-block comes from subr. wnd_frc_slt_get
            ! purpose: compute the saltating friction velocity
            ! theory: saltation roughens the boundary layer, AKA "Owen's effect"

            if (u10(p) >= wnd_rfr_thr_slt) then
               wnd_rfr_dlt = u10(p) - wnd_rfr_thr_slt
               wnd_frc_slt_dlt = 0.003_r8 * wnd_rfr_dlt * wnd_rfr_dlt
               wnd_frc_slt = fv(p) + wnd_frc_slt_dlt
            end if

            ! the following comes from subr. flx_mss_hrz_slt_ttl_Whi79_get
            ! purpose: compute vertically integrated streamwise mass flux of particles

            if (wnd_frc_slt > wnd_frc_thr_slt) then
               wnd_frc_rat = wnd_frc_thr_slt / wnd_frc_slt
               flx_mss_hrz_slt_ttl = cst_slt * forc_rho(c) * (wnd_frc_slt**3.0_r8) * &
                    (1.0_r8 - wnd_frc_rat) * (1.0_r8 + wnd_frc_rat) * (1.0_r8 + wnd_frc_rat) / grav

               ! the following loop originates from subr. dst_mbl
               ! purpose: apply land sfc and veg limitations and global tuning factor
               ! slevis: multiply flx_mss_hrz_slt_ttl by liqfrac to incude the effect
               ! of frozen soil

               flx_mss_hrz_slt_ttl = flx_mss_hrz_slt_ttl * lnd_frc_mbl(p) * mbl_bsn_fct(c) * &
                    flx_mss_fdg_fct * liqfrac
            end if

            ! the following comes from subr. flx_mss_vrt_dst_ttl_MaB95_get
            ! purpose: diagnose total vertical mass flux of dust from vertically
            !          integrated streamwise mass flux

            dst_slt_flx_rat_ttl = 100.0_r8 * exp( log(10.0_r8) * (13.4_r8 * mss_frc_cly_vld(c) - 6.0_r8) )
            flx_mss_vrt_dst_ttl(p) = flx_mss_hrz_slt_ttl * dst_slt_flx_rat_ttl

         end if   ! lnd_frc_mbl > 0.0

      end do

      ! the following comes from subr. flx_mss_vrt_dst_prt in C. Zender's code
      ! purpose: partition total vertical mass flux of dust into transport bins

      do n = 1, ndst
         do m = 1, dst_src_nbr
            do fp = 1,num_nolakep
               p = filter_nolakep(fp)
               if (lnd_frc_mbl(p) > 0.0_r8) then
                  flx_mss_vrt_dst(p,n) = flx_mss_vrt_dst(p,n) +  this%ovr_src_snk_mss(m,n) * flx_mss_vrt_dst_ttl(p)
               end if
            end do
         end do
      end do

      do n = 1, ndst
         do fp = 1,num_nolakep
            p = filter_nolakep(fp)
            if (lnd_frc_mbl(p) > 0.0_r8) then
               flx_mss_vrt_dst_tot(p) = flx_mss_vrt_dst_tot(p) + flx_mss_vrt_dst(p,n)
            end if
         end do
      end do

    end associate

  end subroutine DustEmission

  !------------------------------------------------------------------------

end module DustEmisZender2003