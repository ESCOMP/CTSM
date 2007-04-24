#include <misc.h>
#include <preproc.h>

module SurfaceAlbedoMod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: SurfaceAlbedoMod
!
! !DESCRIPTION:
! Performs surface albedo calculations
!
! !PUBLIC TYPES:
  use clm_varcon , only : istsoil
  use spmdMod    , only : masterproc
  use shr_kind_mod, only : r8 => shr_kind_r8
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: SurfaceAlbedo  ! Surface albedo and two-stream fluxes
  public :: SnowAge        ! Update snow age
!
! !PRIVATE MEMBER FUNCTIONS:
  private :: SnowAlbedo    ! Determine snow albedos
  private :: SoilAlbedo    ! Determine ground surface albedo
  private :: TwoStream     ! Two-stream fluxes for canopy radiative transfer
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: SurfaceAlbedo
!
! !INTERFACE:
  subroutine SurfaceAlbedo(lbg, ubg, lbc, ubc, lbp, ubp, &
       caldayp1, declinp1)
!
! !DESCRIPTION:
! Surface albedo and two-stream fluxes
! Surface albedos. Also fluxes (per unit incoming direct and diffuse
! radiation) reflected, transmitted, and absorbed by vegetation.
! Also sunlit fraction of the canopy.
! The calling sequence is:
! -> SurfaceAlbedo:   albedos for next time step
!    -> SnowAlbedo:   snow albedos: direct beam
!    -> SnowAlbedo:   snow albedos: diffuse
!    -> SoilAlbedo:   soil/lake/glacier/wetland albedos
!    -> TwoStream:    absorbed, reflected, transmitted solar fluxes (vis dir,vis dif, nir dir, nir dif)
!
! !USES:
    use clmtype
    use clm_varpar      , only : numrad
    use shr_orb_mod
    use clm_time_manager, only : get_nstep
!
! !ARGUMENTS:
    implicit none
    integer , intent(in) :: lbg, ubg ! gridcell bounds
    integer , intent(in) :: lbc, ubc ! column bounds
    integer , intent(in) :: lbp, ubp ! pft bounds
    real(r8), intent(in) :: caldayp1 ! calendar day at Greenwich (1.00, ..., 365.99)
    real(r8), intent(in) :: declinp1 ! declination angle (radians) for next time step
!
! !CALLED FROM:
! subroutine lpjreset1 in module DGVMMod (only applicable when cpp token DGVM is defined)
! subroutine driver
! subroutine iniTimeVar
!
! !REVISION HISTORY:
! Author: Gordon Bonan
! 2/1/02, Peter Thornton: Migrate to new data structures
! 8/20/03, Mariana Vertenstein: Vectorized routine
! 11/3/03, Peter Thornton: added decl(c) output for use in CN code.
!
! !LOCAL VARIABLES:
!
! local pointers to implicit in arguments
!
    integer , pointer :: pgridcell(:) ! gridcell of corresponding pft
    integer , pointer :: plandunit(:) ! index into landunit level quantities
    integer , pointer :: itypelun(:)  ! landunit type
    integer , pointer :: pcolumn(:)   ! column of corresponding pft
    integer , pointer :: cgridcell(:) ! gridcell of corresponding column
    real(r8), pointer :: pwtgcell(:)  ! weight of pft wrt corresponding gridcell
    real(r8), pointer :: lat(:)       ! gridcell latitude (radians)
    real(r8), pointer :: lon(:)       ! gridcell longitude (radians)
    real(r8), pointer :: elai(:)      ! one-sided leaf area index with burying by snow
    real(r8), pointer :: esai(:)      ! one-sided stem area index with burying by snow
    real(r8), pointer :: h2osno(:)    ! snow water (mm H2O)
    real(r8), pointer :: snowage(:)   ! non dimensional snow age [-]
    real(r8), pointer :: rhol(:,:)    ! leaf reflectance: 1=vis, 2=nir
    real(r8), pointer :: rhos(:,:)    ! stem reflectance: 1=vis, 2=nir
    real(r8), pointer :: taul(:,:)    ! leaf transmittance: 1=vis, 2=nir
    real(r8), pointer :: taus(:,:)    ! stem transmittance: 1=vis, 2=nir
    integer , pointer :: ivt(:)       ! pft vegetation type
!
! local pointers toimplicit out arguments
!
    real(r8), pointer :: coszen(:)	  ! cosine of solar zenith angle
    real(r8), pointer :: fsun(:)      ! sunlit fraction of canopy
    real(r8), pointer :: albgrd(:,:)  ! ground albedo (direct)
    real(r8), pointer :: albgri(:,:)  ! ground albedo (diffuse)
    real(r8), pointer :: albd(:,:)    ! surface albedo (direct)
    real(r8), pointer :: albi(:,:)    ! surface albedo (diffuse)
    real(r8), pointer :: fabd(:,:)    ! flux absorbed by veg per unit direct flux
    real(r8), pointer :: fabi(:,:)    ! flux absorbed by veg per unit diffuse flux
    real(r8), pointer :: ftdd(:,:)    ! down direct flux below veg per unit dir flx
    real(r8), pointer :: ftid(:,:)    ! down diffuse flux below veg per unit dir flx
    real(r8), pointer :: ftii(:,:)    ! down diffuse flux below veg per unit dif flx
    real(r8), pointer :: decl(:)      ! solar declination angle (radians)
    real(r8), pointer :: gdir(:)      ! leaf projection in solar direction (0 to 1)
    real(r8), pointer :: omega(:,:)   ! fraction of intercepted radiation that is scattered (0 to 1)
!
!EOP
!
! !OTHER LOCAL VARIABLES:
!
    real(r8), parameter :: mpe = 1.e-06_r8    ! prevents overflow for division by zero
    integer  :: fp,g,c,p                   ! indices
    integer  :: ib                         ! band index
    integer  :: ic                         ! 0=unit incoming direct; 1=unit incoming diffuse
    real(r8) :: wl(lbp:ubp)                ! fraction of LAI+SAI that is LAI
    real(r8) :: ws(lbp:ubp)                ! fraction of LAI+SAI that is SAI
    real(r8) :: vai(lbp:ubp)               ! elai+esai
    real(r8) :: rho(lbp:ubp,numrad)        ! leaf/stem refl weighted by fraction LAI and SAI
    real(r8) :: tau(lbp:ubp,numrad)        ! leaf/stem tran weighted by fraction LAI and SAI
    real(r8) :: ftdi(lbp:ubp,numrad)       ! down direct flux below veg per unit dif flux = 0
    real(r8) :: albsnd(lbc:ubc,numrad)     ! snow albedo (direct)
    real(r8) :: albsni(lbc:ubc,numrad)     ! snow albedo (diffuse)
    real(r8) :: ext(lbp:ubp)               ! optical depth direct beam per unit LAI+SAI
    real(r8) :: coszen_gcell(lbg:ubg)      ! cosine solar zenith angle for next time step (gridcell level)
    real(r8) :: coszen_col(lbc:ubc)        ! cosine solar zenith angle for next time step (pft level)
    real(r8) :: coszen_pft(lbp:ubp)        ! cosine solar zenith angle for next time step (pft level)
    integer  :: num_vegsol                 ! number of vegetated pfts where coszen>0
    integer  :: filter_vegsol(ubp-lbp+1)   ! pft filter where vegetated and coszen>0
    integer  :: num_novegsol               ! number of vegetated pfts where coszen>0
    integer  :: filter_novegsol(ubp-lbp+1) ! pft filter where vegetated and coszen>0
    integer  :: num_solar                  ! number of gridcells where coszen>0
  !-----------------------------------------------------------------------

    ! Assign local pointers to derived subtypes components (gridcell-level)

    lat       => clm3%g%lat_a
    lon       => clm3%g%lon_a

    ! Assign local pointers to derived subtypes components (landunit level)

    itypelun       => clm3%g%l%itype

    ! Assign local pointers to derived subtypes components (column-level)

    cgridcell => clm3%g%l%c%gridcell
    h2osno    => clm3%g%l%c%cws%h2osno
    snowage   => clm3%g%l%c%cps%snowage
    albgrd    => clm3%g%l%c%cps%albgrd
    albgri    => clm3%g%l%c%cps%albgri
    decl      => clm3%g%l%c%cps%decl 
    coszen    => clm3%g%l%c%cps%coszen 

    ! Assign local pointers to derived subtypes components (pft-level)

    plandunit => clm3%g%l%c%p%landunit
    pgridcell => clm3%g%l%c%p%gridcell
    pcolumn   => clm3%g%l%c%p%column
    pwtgcell  => clm3%g%l%c%p%wtgcell
    albd      => clm3%g%l%c%p%pps%albd
    albi      => clm3%g%l%c%p%pps%albi
    fabd      => clm3%g%l%c%p%pps%fabd
    fabi      => clm3%g%l%c%p%pps%fabi
    ftdd      => clm3%g%l%c%p%pps%ftdd
    ftid      => clm3%g%l%c%p%pps%ftid
    ftii      => clm3%g%l%c%p%pps%ftii
    fsun      => clm3%g%l%c%p%pps%fsun
    elai      => clm3%g%l%c%p%pps%elai
    esai      => clm3%g%l%c%p%pps%esai
    gdir      => clm3%g%l%c%p%pps%gdir
    omega     => clm3%g%l%c%p%pps%omega
    ivt       => clm3%g%l%c%p%itype
    rhol      => pftcon%rhol
    rhos      => pftcon%rhos
    taul      => pftcon%taul
    taus      => pftcon%taus

    ! Cosine solar zenith angle for next time step

!dir$ concurrent
!cdir nodep
    do g = lbg, ubg
       coszen_gcell(g) = shr_orb_cosz (caldayp1, lat(g), lon(g), declinp1)
    end do

    ! Save coszen and declination values to  clm3 data structures for
    ! use in other places in the CN code

!dir$ concurrent
!cdir nodep
    do c = lbc, ubc
       g = cgridcell(c)
       coszen_col(c) = coszen_gcell(g)
       coszen(c) = coszen_col(c)
       decl(c) = declinp1
    end do

!dir$ concurrent
!cdir nodep
    do p = lbp, ubp
       g = pgridcell(p)
       coszen_pft(p) = coszen_gcell(g)
    end do

    ! Initialize output because solar radiation only done if coszen > 0

    do ib = 1, numrad
!dir$ concurrent
!cdir nodep
       do c = lbc,ubc
          albgrd(c,ib) = 0._r8
          albgri(c,ib) = 0._r8
       end do
!dir$ concurrent
!cdir nodep
       do p = lbp,ubp
          albd(p,ib) = 1._r8
          albi(p,ib) = 1._r8
          fabd(p,ib) = 0._r8
          fabi(p,ib) = 0._r8
          ftdd(p,ib) = 0._r8
          ftid(p,ib) = 0._r8
          ftii(p,ib) = 0._r8
          omega(p,ib)= 0._r8
          if (ib==1) then
             gdir(p) = 0._r8
          end if
       end do
    end do

    ! Index points with positive coszen for subsequent calculations
    ! Return if all coszen are not positive

    num_solar = 0
!dir$ concurrent
!cdir nodep
    do g = lbg,ubg
       if (coszen_gcell(g) > 0._r8) num_solar = num_solar + 1
    end do
    if (num_solar <= 0._r8) return

    ! Snow albedos
    ! Note that snow albedo routine will only compute nonzero snow albedos
    ! where h2osno> 0 and coszen > 0

    ic = 0; call SnowAlbedo(lbc, ubc, coszen_col, ic, albsnd)
    ic = 1; call SnowAlbedo(lbc, ubc, coszen_col, ic, albsni)

    ! Ground surface albedos
    ! Note that ground albedo routine will only compute nonzero snow albedos
    ! where coszen > 0

    call SoilAlbedo(lbc, ubc, coszen_col, albsnd, albsni)

    ! Creat solar-vegetated filter for the following calculations

    num_vegsol = 0
    num_novegsol = 0
    do p = lbp,ubp
       if (pwtgcell(p)>0._r8) then
          if (coszen_pft(p) > 0._r8) then
             if (itypelun(plandunit(p)) == istsoil .and. (elai(p) + esai(p)) > 0._r8) then                       
                num_vegsol = num_vegsol + 1
                filter_vegsol(num_vegsol) = p
             else if (itypelun(plandunit(p)) /= istsoil .or. &
                     (itypelun(plandunit(p)) == istsoil .and. (elai(p) + esai(p)) == 0._r8)) then
                num_novegsol = num_novegsol + 1
                filter_novegsol(num_novegsol) = p
             end if
          end if
       end if
    end do

    ! Weight reflectance/transmittance by lai and sai
    ! Only perform on vegetated pfts where coszen > 0

!dir$ concurrent
!cdir nodep
    do fp = 1,num_vegsol
       p = filter_vegsol(fp)
       vai(p) = elai(p) + esai(p)
       wl(p) = elai(p) / max( vai(p), mpe )
       ws(p) = esai(p) / max( vai(p), mpe )
    end do

    do ib = 1, numrad
!dir$ concurrent
!cdir nodep
       do fp = 1,num_vegsol
          p = filter_vegsol(fp)
          rho(p,ib) = max( rhol(ivt(p),ib)*wl(p) + rhos(ivt(p),ib)*ws(p), mpe )
          tau(p,ib) = max( taul(ivt(p),ib)*wl(p) + taus(ivt(p),ib)*ws(p), mpe )
       end do
    end do

    ! Calculate surface albedos and fluxes 
    ! Only perform on vegetated pfts where coszen > 0

    call TwoStream (lbc, ubc, lbp, ubp, filter_vegsol, num_vegsol, &
                    coszen_pft, vai, rho, tau)
       
    ! Determine values for non-vegetated pfts where coszen > 0

    do ib = 1,numrad
!dir$ concurrent
!cdir nodep
       do fp = 1,num_novegsol
          p = filter_novegsol(fp)
          c = pcolumn(p)
          fabd(p,ib) = 0._r8
          fabi(p,ib) = 0._r8
          ftdd(p,ib) = 1._r8
          ftid(p,ib) = 0._r8
          ftii(p,ib) = 1._r8
          albd(p,ib) = albgrd(c,ib)
          albi(p,ib) = albgri(c,ib)
          gdir(p) = 0._r8
       end do
    end do

  end subroutine SurfaceAlbedo

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: SnowAlbedo
!
! !INTERFACE:
  subroutine SnowAlbedo (lbc, ubc, coszen, ind, alb)
!
! !DESCRIPTION:
! Determine snow albedos
!
! !USES:
    use clmtype
!
! !ARGUMENTS:
    implicit none
    integer , intent(in) :: lbc, ubc                 ! column bounds
    real(r8), intent(in) :: coszen(lbc:ubc)          ! cosine solar zenith angle for next time step
    integer , intent(in) :: ind                      ! 0=direct beam, 1=diffuse radiation
    real(r8), intent(out):: alb(lbc:ubc,2)           ! snow albedo by waveband (assume 2 wavebands)
!
! !CALLED FROM:
! subroutine SurfaceAlbedo in this module
!
! !REVISION HISTORY:
! Author: Gordon Bonan
! 2/5/02, Peter Thornton: Migrated to new data structures. Eliminated
! reference to derived types in this subroutine, and made consistent use
! of the assumption that numrad = 2, with index values: 1=visible,2=NIR
! 8/20/03, Mariana Vertenstein: Vectorized routine
!
!
! !LOCAL VARIABLES:
!
! local pointers to implicit in arguments
!
    real(r8), pointer :: h2osno(:)    ! snow water (mm H2O)
    real(r8), pointer :: snowage(:)   ! non dimensional snow age [-]
!
!EOP
!
! !OTHER LOCAL VARIABLES:
!
! variables and constants for snow albedo calculation
!
    real(r8), parameter :: snal0 = 0.95_r8 ! vis albedo of new snow for sza<60
    real(r8), parameter :: snal1 = 0.65_r8 ! nir albedo of new snow for sza<60
    real(r8), parameter :: cons  = 0.2_r8  ! constant for visible snow albedo calculation [-]
    real(r8), parameter :: conn  = 0.5_r8  ! constant for nir snow albedo calculation [-]
    real(r8), parameter :: sl    = 2.0_r8  ! factor that helps control alb zenith dependence [-]
    integer  :: c                       ! index
    real(r8) :: age                     ! factor to reduce visible snow alb due to snow age [-]
    real(r8) :: albs                    ! temporary vis snow albedo
    real(r8) :: albl                    ! temporary nir snow albedo
    real(r8) :: cff                     ! snow alb correction factor for zenith angle > 60 [-]
    real(r8) :: czf                     ! solar zenith correction for new snow albedo [-]
!-----------------------------------------------------------------------

    ! Assign local pointers to derived subtypes components (column-level)

    h2osno  => clm3%g%l%c%cws%h2osno
    snowage => clm3%g%l%c%cps%snowage

    ! this code assumes that numrad = 2 , with the following
    ! index values: 1 = visible, 2 = NIR

    ! Albedo for snow cover.
    ! Snow albedo depends on snow-age, zenith angle, and thickness of snow,
    ! age gives reduction of visible radiation
    ! age below is correction for snow age
    ! czf below corrects albedo of new snow for solar zenith

!dir$ concurrent
!cdir nodep
    do c = lbc, ubc
       if (coszen(c) > 0._r8 .and. h2osno(c) > 0._r8) then
          age = 1._r8-1._r8/(1._r8+snowage(c))
          albs = snal0*(1._r8-cons*age)
          albl = snal1*(1._r8-conn*age)
          if (ind == 0) then
             cff  = ((1._r8+1._r8/sl)/(1._r8+max(0.001_r8, coszen(c))*2._r8*sl )- 1._r8/sl)
             cff  = max(cff, 0._r8)
             czf  = 0.4_r8 * cff * (1._r8-albs)
             albs = albs + czf
             czf  = 0.4_r8 *cff * (1._r8-albl)
             albl = albl + czf
          end if
          alb(c,1) = albs
          alb(c,2) = albl
       else
          alb(c,1) = 0._r8
          alb(c,2) = 0._r8
       end if
    end do

  end subroutine SnowAlbedo

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: SoilAlbedo
!
! !INTERFACE:
  subroutine SoilAlbedo (lbc, ubc, coszen, albsnd, albsni)
!
! !DESCRIPTION:
! Determine ground surface albedo, accounting for snow
!
! !USES:
    use clmtype
    use clm_varpar, only : numrad
    use clm_varcon, only : albsat, albdry, alblak, albice, tfrz, istice, istsoil
!
! !ARGUMENTS:
    implicit none
    integer , intent(in) :: lbc, ubc                ! column bounds
    real(r8), intent(in) :: coszen(lbc:ubc)         ! cos solar zenith angle next time step (column-level)
    real(r8), intent(in) :: albsnd(lbc:ubc,numrad)  ! snow albedo (direct)
    real(r8), intent(in) :: albsni(lbc:ubc,numrad)  ! snow albedo (diffuse)
!
! !CALLED FROM:
! subroutine SurfaceAlbedo in this module
!
! !REVISION HISTORY:
! Author: Gordon Bonan
! 2/5/02, Peter Thornton: Migrated to new data structures.
! 8/20/03, Mariana Vertenstein: Vectorized routine
!
! !LOCAL VARIABLES:
!
! local pointers to original implicit in arguments
!
    integer , pointer :: clandunit(:)    ! landunit of corresponding column
    integer , pointer :: ltype(:)        ! landunit type
    integer , pointer :: isoicol(:)      ! soil color class
    real(r8), pointer :: t_grnd(:)       ! ground temperature (Kelvin)
    real(r8), pointer :: frac_sno(:)     ! fraction of ground covered by snow (0 to 1)
    real(r8), pointer :: h2osoi_vol(:,:) ! volumetric soil water [m3/m3]
!
! local pointers to original implicit out arguments
!
    real(r8), pointer:: albgrd(:,:)      ! ground albedo (direct)
    real(r8), pointer:: albgri(:,:)      ! ground albedo (diffuse)
!
!EOP
!
! !OTHER LOCAL VARIABLES:
!
    integer, parameter :: nband =numrad ! number of solar radiation waveband classes
    integer  :: c,l           ! indices
    integer  :: ib            ! waveband number (1=vis, 2=nir)
    real(r8) :: inc           ! soil water correction factor for soil albedo
    real(r8) :: albsod        ! soil albedo (direct)
    real(r8) :: albsoi        ! soil albedo (diffuse)
    integer  :: soilcol       ! soilcolor
!-----------------------------------------------------------------------
!dir$ inlinenever SoilAlbedo

    ! Assign local pointers to derived subtypes components (column-level)

    clandunit  => clm3%g%l%c%landunit
    isoicol    => clm3%g%l%c%cps%isoicol
    t_grnd     => clm3%g%l%c%ces%t_grnd
    frac_sno   => clm3%g%l%c%cps%frac_sno
    h2osoi_vol => clm3%g%l%c%cws%h2osoi_vol
    albgrd     => clm3%g%l%c%cps%albgrd
    albgri     => clm3%g%l%c%cps%albgri

    ! Assign local pointers to derived subtypes components (landunit-level)

    ltype      => clm3%g%l%itype

    ! Compute soil albedos

    do ib = 1, nband
!dir$ concurrent
!cdir nodep
       do c = lbc, ubc
          if (coszen(c) > 0._r8) then
             l = clandunit(c)

             if (ltype(l) == istsoil)  then              ! soil
                inc    = max(0.11_r8-0.40_r8*h2osoi_vol(c,1), 0._r8)
                soilcol = isoicol(c)
                albsod = min(albsat(soilcol,ib)+inc, albdry(soilcol,ib))
                albsoi = albsod
             else if (ltype(l) == istice)  then          ! land ice
                albsod = albice(ib)
                albsoi = albsod
             else if (t_grnd(c) > tfrz) then             ! unfrozen lake, wetland
                albsod = 0.05_r8/(max(0.001_r8,coszen(c)) + 0.15_r8)
                albsoi = albsod
             else                                        ! frozen lake, wetland
                albsod = alblak(ib)
                albsoi = albsod
             end if

             albgrd(c,ib) = albsod*(1._r8-frac_sno(c)) + albsnd(c,ib)*frac_sno(c)
             albgri(c,ib) = albsoi*(1._r8-frac_sno(c)) + albsni(c,ib)*frac_sno(c)
          end if
       end do
    end do

  end subroutine SoilAlbedo

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: TwoStream
!
! !INTERFACE:
  subroutine TwoStream (lbc, ubc, lbp, ubp, filter_vegsol, num_vegsol, &
                        coszen, vai, rho, tau)
!
! !DESCRIPTION:
! Two-stream fluxes for canopy radiative transfer
! Use two-stream approximation of Dickinson (1983) Adv Geophysics
! 25:305-353 and Sellers (1985) Int J Remote Sensing 6:1335-1372
! to calculate fluxes absorbed by vegetation, reflected by vegetation,
! and transmitted through vegetation for unit incoming direct or diffuse
! flux given an underlying surface with known albedo.
!
! !USES:
    use clmtype
    use clm_varpar, only : numrad
    use clm_varcon, only : omegas, tfrz, betads, betais
!
! !ARGUMENTS:
    implicit none
    integer , intent(in)  :: lbc, ubc                 ! column bounds
    integer , intent(in)  :: lbp, ubp                 ! pft bounds
    integer , intent(in)  :: filter_vegsol(ubp-lbp+1) ! filter for vegetated pfts with coszen>0
    integer , intent(in)  :: num_vegsol               ! number of vegetated pfts where coszen>0
    real(r8), intent(in)  :: coszen(lbp:ubp)          ! cosine solar zenith angle for next time step
    real(r8), intent(in)  :: vai(lbp:ubp)             ! elai+esai
    real(r8), intent(in)  :: rho(lbp:ubp,numrad)      ! leaf/stem refl weighted by fraction LAI and SAI
    real(r8), intent(in)  :: tau(lbp:ubp,numrad)      ! leaf/stem tran weighted by fraction LAI and SAI
!
! !CALLED FROM:
! subroutine SurfaceAlbedo in this module
!
! !REVISION HISTORY:
! Author: Gordon Bonan
! Modified for speedup: Mariana Vertenstein, 8/26/02
! Vectorized routine: Mariana Vertenstein:  8/20/03
!
! !LOCAL VARIABLES:
!
! local pointers to implicit in scalars
!
    integer , pointer :: pcolumn(:)    ! column of corresponding pft
    real(r8), pointer :: albgrd(:,:)   ! ground albedo (direct) (column-level)
    real(r8), pointer :: albgri(:,:)   ! ground albedo (diffuse)(column-level)
    real(r8), pointer :: t_veg(:)      ! vegetation temperature (Kelvin)
    real(r8), pointer :: fwet(:)       ! fraction of canopy that is wet (0 to 1)
    integer , pointer :: ivt(:)        ! pft vegetation type
    real(r8), pointer :: xl(:)         ! ecophys const - leaf/stem orientation index
!
! local pointers to implicit out scalars
!
    real(r8), pointer :: albd(:,:)     ! surface albedo (direct)
    real(r8), pointer :: albi(:,:)     ! surface albedo (diffuse)
    real(r8), pointer :: fabd(:,:)     ! flux absorbed by veg per unit direct flux
    real(r8), pointer :: fabi(:,:)     ! flux absorbed by veg per unit diffuse flux
    real(r8), pointer :: ftdd(:,:)     ! down direct flux below veg per unit dir flx
    real(r8), pointer :: ftid(:,:)     ! down diffuse flux below veg per unit dir flx
    real(r8), pointer :: ftii(:,:)     ! down diffuse flux below veg per unit dif flx
    real(r8), pointer :: gdir(:)		   ! leaf projection in solar direction (0 to 1)
	 real(r8), pointer :: omega(:,:)    ! fraction of intercepted radiation that is scattered (0 to 1)
!
!EOP
!
! !OTHER LOCAL VARIABLES:
!
    integer  :: fp,p,c           ! array indices
    !integer  :: ic               ! 0=unit incoming direct; 1=unit incoming diffuse
    integer  :: ib               ! waveband number
    real(r8) :: cosz             ! 0.001 <= coszen <= 1.000
    real(r8) :: asu              ! single scattering albedo
    real(r8) :: chil(lbp:ubp)    ! -0.4 <= xl <= 0.6
    real(r8) :: twostext(lbp:ubp)! optical depth of direct beam per unit leaf area
    real(r8) :: avmu(lbp:ubp)    ! average diffuse optical depth
    real(r8) :: omegal           ! omega for leaves
    real(r8) :: betai            ! upscatter parameter for diffuse radiation
    real(r8) :: betail           ! betai for leaves
    real(r8) :: betad            ! upscatter parameter for direct beam radiation
    real(r8) :: betadl           ! betad for leaves
    real(r8) :: tmp0,tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8,tmp9 ! temporary
    real(r8) :: p1,p2,p3,p4,s1,s2,u1,u2,u3                        ! temporary
    real(r8) :: b,c1,d,d1,d2,f,h,h1,h2,h3,h4,h5,h6,h7,h8,h9,h10   ! temporary
    real(r8) :: phi1,phi2,sigma                                   ! temporary
    real(r8) :: temp0(lbp:ubp),temp1,temp2(lbp:ubp)               ! temporary
    real(r8) :: t1
!-----------------------------------------------------------------------

    ! Assign local pointers to derived subtypes components (column-level)

    albgrd  => clm3%g%l%c%cps%albgrd
    albgri  => clm3%g%l%c%cps%albgri

    ! Assign local pointers to derived subtypes components (pft-level)

    pcolumn => clm3%g%l%c%p%column
    fwet    => clm3%g%l%c%p%pps%fwet
    t_veg   => clm3%g%l%c%p%pes%t_veg
    ivt     => clm3%g%l%c%p%itype
    albd    => clm3%g%l%c%p%pps%albd
    albi    => clm3%g%l%c%p%pps%albi
    fabd    => clm3%g%l%c%p%pps%fabd
    fabi    => clm3%g%l%c%p%pps%fabi
    ftdd    => clm3%g%l%c%p%pps%ftdd
    ftid    => clm3%g%l%c%p%pps%ftid
    ftii    => clm3%g%l%c%p%pps%ftii
    gdir    => clm3%g%l%c%p%pps%gdir
    omega   => clm3%g%l%c%p%pps%omega
    xl      => pftcon%xl

    ! Calculate two-stream parameters omega, betad, betai, avmu, gdir, twostext.
    ! Omega, betad, betai are adjusted for snow. Values for omega*betad
    ! and omega*betai are calculated and then divided by the new omega
    ! because the product omega*betai, omega*betad is used in solution.
    ! Also, the transmittances and reflectances (tau, rho) are linear
    ! weights of leaf and stem values.

!dir$ concurrent
!cdir nodep
    do fp = 1,num_vegsol
       p = filter_vegsol(fp)
       
       ! note that the following limit only acts on cosz values > 0 and less than 
       ! 0.001, not on values cosz = 0, since these zero have already been filtered
       ! out in filter_vegsol
       cosz = max(0.001_r8, coszen(p))
       
       chil(p) = min( max(xl(ivt(p)), -0.4_r8), 0.6_r8 )
       if (abs(chil(p)) <= 0.01_r8) chil(p) = 0.01_r8
       phi1 = 0.5_r8 - 0.633_r8*chil(p) - 0.330_r8*chil(p)*chil(p)
       phi2 = 0.877_r8 * (1._r8-2._r8*phi1)
       gdir(p) = phi1 + phi2*cosz
       twostext(p) = gdir(p)/cosz
       avmu(p) = ( 1._r8 - phi1/phi2 * log((phi1+phi2)/phi1) ) / phi2
       temp0(p) = gdir(p) + phi2*cosz
       temp1 = phi1*cosz
       temp2(p) = ( 1._r8 - temp1/temp0(p) * log((temp1+temp0(p))/temp1) )
    end do

    do ib = 1, numrad
!dir$ concurrent
!cdir nodep
       do fp = 1,num_vegsol
          p = filter_vegsol(fp)
          c = pcolumn(p)

          omegal = rho(p,ib) + tau(p,ib)
          asu = 0.5_r8*omegal*gdir(p)/temp0(p) *temp2(p)
          betadl = (1._r8+avmu(p)*twostext(p))/(omegal*avmu(p)*twostext(p))*asu
          betail = 0.5_r8 * ((rho(p,ib)+tau(p,ib)) + (rho(p,ib)-tau(p,ib)) &
               * ((1._r8+chil(p))/2._r8)**2) / omegal

          ! Adjust omega, betad, and betai for intercepted snow

          if (t_veg(p) > tfrz) then                             !no snow
             tmp0 = omegal
             tmp1 = betadl
             tmp2 = betail
          else
             tmp0 =   (1._r8-fwet(p))*omegal        + fwet(p)*omegas(ib)
             tmp1 = ( (1._r8-fwet(p))*omegal*betadl + fwet(p)*omegas(ib)*betads ) / tmp0
             tmp2 = ( (1._r8-fwet(p))*omegal*betail + fwet(p)*omegas(ib)*betais ) / tmp0
          end if
          omega(p,ib) = tmp0           
          betad = tmp1 
          betai = tmp2  

          ! Absorbed, reflected, transmitted fluxes per unit incoming radiation

          b = 1._r8 - omega(p,ib) + omega(p,ib)*betai
          c1 = omega(p,ib)*betai
          tmp0 = avmu(p)*twostext(p)
          d = tmp0 * omega(p,ib)*betad
          f = tmp0 * omega(p,ib)*(1._r8-betad)
          tmp1 = b*b - c1*c1
          h = sqrt(tmp1) / avmu(p)
          sigma = tmp0*tmp0 - tmp1
          p1 = b + avmu(p)*h
          p2 = b - avmu(p)*h
          p3 = b + tmp0
          p4 = b - tmp0
          
#if (defined CN)          
          ! PET, 3/1/04: added this test to avoid floating point errors in exp()
          t1 = min(h*vai(p), 40._r8)
          s1 = exp(-t1)
          t1 = min(twostext(p)*vai(p), 40._r8)
          s2 = exp(-t1)
#else
          s1 = exp(-h*vai(p))
          s2 = exp(-twostext(p)*vai(p))
#endif
          
          ! Determine fluxes for vegetated pft for unit incoming direct 
          ! Loop over incoming direct and incoming diffuse
          ! 0=unit incoming direct; 1=unit incoming diffuse

          ! ic = 0 unit incoming direct flux
          ! ========================================

          u1 = b - c1/albgrd(c,ib)
          u2 = b - c1*albgrd(c,ib)
          u3 = f + c1*albgrd(c,ib)

          tmp2 = u1 - avmu(p)*h
          tmp3 = u1 + avmu(p)*h
          d1 = p1*tmp2/s1 - p2*tmp3*s1
          tmp4 = u2 + avmu(p)*h
          tmp5 = u2 - avmu(p)*h
          d2 = tmp4/s1 - tmp5*s1
          h1 = -d*p4 - c1*f
          tmp6 = d - h1*p3/sigma
          tmp7 = ( d - c1 - h1/sigma*(u1+tmp0) ) * s2
          h2 = ( tmp6*tmp2/s1 - p2*tmp7 ) / d1
          h3 = - ( tmp6*tmp3*s1 - p1*tmp7 ) / d1
          h4 = -f*p3 - c1*d
          tmp8 = h4/sigma
          tmp9 = ( u3 - tmp8*(u2-tmp0) ) * s2
          h5 = - ( tmp8*tmp4/s1 + tmp9 ) / d2
          h6 = ( tmp8*tmp5*s1 + tmp9 ) / d2
          h7 = (c1*tmp2) / (d1*s1)
          h8 = (-c1*tmp3*s1) / d1
          h9 = tmp4 / (d2*s1)
          h10 = (-tmp5*s1) / d2

          ! Downward direct and diffuse fluxes below vegetation (ic = 0)

          ftdd(p,ib) = s2
          ftid(p,ib) = h4*s2/sigma + h5*s1 + h6/s1

          ! Flux reflected by vegetation (ic = 0)

          albd(p,ib) = h1/sigma + h2 + h3

          ! Flux absorbed by vegetation (ic = 0)

          fabd(p,ib) = 1._r8 - albd(p,ib) &
               - (1._r8-albgrd(c,ib))*ftdd(p,ib) - (1._r8-albgri(c,ib))*ftid(p,ib)

          ! ic = 1 unit incoming diffuse
          ! ========================================

          u1 = b - c1/albgri(c,ib)
          u2 = b - c1*albgri(c,ib)
          u3 = f + c1*albgri(c,ib)

          tmp2 = u1 - avmu(p)*h
          tmp3 = u1 + avmu(p)*h
          d1 = p1*tmp2/s1 - p2*tmp3*s1
          tmp4 = u2 + avmu(p)*h
          tmp5 = u2 - avmu(p)*h
          d2 = tmp4/s1 - tmp5*s1
          h1 = -d*p4 - c1*f
          tmp6 = d - h1*p3/sigma
          tmp7 = ( d - c1 - h1/sigma*(u1+tmp0) ) * s2
          h2 = ( tmp6*tmp2/s1 - p2*tmp7 ) / d1
          h3 = - ( tmp6*tmp3*s1 - p1*tmp7 ) / d1
          h4 = -f*p3 - c1*d
          tmp8 = h4/sigma
          tmp9 = ( u3 - tmp8*(u2-tmp0) ) * s2
          h5 = - ( tmp8*tmp4/s1 + tmp9 ) / d2
          h6 = ( tmp8*tmp5*s1 + tmp9 ) / d2
          h7 = (c1*tmp2) / (d1*s1)
          h8 = (-c1*tmp3*s1) / d1
          h9 = tmp4 / (d2*s1)
          h10 = (-tmp5*s1) / d2

          ! Downward direct and diffuse fluxes below vegetation

          ftii(p,ib) = h9*s1 + h10/s1

          ! Flux reflected by vegetation

          albi(p,ib) = h7 + h8

          ! Flux absorbed by vegetation

          fabi(p,ib) = 1._r8 - albi(p,ib) - (1._r8-albgri(c,ib))*ftii(p,ib)

       end do   ! end of pft loop
    end do   ! end of radiation band loop

  end subroutine TwoStream

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: SnowAge
!
! !INTERFACE:
  subroutine SnowAge (lbc, ubc)
!
! !DESCRIPTION:
! Updates snow age Based on BATS code.
!
! !USES:
    use clmtype
    use clm_varcon, only : tfrz
    use clm_time_manager, only : get_step_size
!
! !ARGUMENTS:
    implicit none
    integer , intent(in) :: lbc, ubc ! column bounds
!
! !CALLED FROM:
!
! !REVISION HISTORY:
! Original Code:  Robert Dickinson
! 15 September 1999: Yongjiu Dai; Integration of code into CLM
! 15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision
! 3/4/02, Peter Thornton: Migrated to new data structures.
! 8/20/03, Mariana Vertenstein: Vectorized routine
!
! !LOCAL VARIABLES:
!
! local pointers to implicit in arguments
!
    real(r8), pointer :: t_grnd(:)     ! ground temperature (Kelvin)
    real(r8), pointer :: h2osno(:)     ! snow water (mm H2O)
    real(r8), pointer :: h2osno_old(:) ! snow mass for previous time step (kg/m2)
!
! local pointers to implicit inout arguments
!
    real(r8), pointer :: snow_age(:)   ! non dimensional snow age [-]
!
!EOP
!
! !OTHER LOCAL VARIABLES:
!
    integer  :: c     ! index
    real(r8) :: age1  ! snow aging factor due to crystal growth [-]
    real(r8) :: age2  ! snow aging factor due to surface growth [-]
    real(r8) :: age3  ! snow aging factor due to accum of other particles [-]
    real(r8) :: arg   ! temporary variable used in snow age calculation [-]
    real(r8) :: arg2  ! temporary variable used in snow age calculation [-]
    real(r8) :: dela  ! temporary variable used in snow age calculation [-]
    real(r8) :: dels  ! temporary variable used in snow age calculation [-]
    real(r8) :: sge   ! temporary variable used in snow age calculation [-]
    real(r8) :: dtime ! land model time step (sec)
!-----------------------------------------------------------------------

   ! Assign local pointers to derived type members (column-level)

    snow_age   => clm3%g%l%c%cps%snowage
    t_grnd     => clm3%g%l%c%ces%t_grnd
    h2osno     => clm3%g%l%c%cws%h2osno
    h2osno_old => clm3%g%l%c%cws%h2osno_old

    ! Determine snow age

    dtime = get_step_size()

!dir$ concurrent
!cdir nodep
    do c = lbc, ubc
       if (h2osno(c) <= 0._r8) then
          snow_age(c) = 0._r8
       else if (h2osno(c) > 800._r8) then       ! Over Antarctica
          snow_age(c) = 0._r8
       else                                  ! Away from Antarctica
          age3 = 0.3_r8
          arg  = 5.e3_r8*(1._r8/tfrz-1._r8/t_grnd(c))
          arg2 = min(0._r8, 10._r8*arg)
          age2 = exp(arg2)
          age1 = exp(arg)
          dela = 1.e-6_r8 * dtime * (age1+age2+age3)
          dels = 0.1_r8*max(0.0_r8,  h2osno(c)-h2osno_old(c))
          sge  = (snow_age(c)+dela) * (1.0_r8-dels)
          snow_age(c) = max(0.0_r8, sge)
       end if
    end do

  end subroutine SnowAge

end module SurfaceAlbedoMod
