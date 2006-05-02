#include <misc.h>
#include <preproc.h>

module lnd_comp_mct

#if (defined COUP_CAM)

!---------------------------------------------------------------------------
!BOP
!
! !MODULE: lnd_comp_mct
!
! !DESCRIPTION:
! This interface 
!
! !USES:
  use seq_mct_mod
  use seq_mct_init
  use seq_flds_mod
  use seq_flds_indices
  use shr_kind_mod,   only : r8 => shr_kind_r8
  use abortutils,     only : endrun
#if (defined SPMD)
  use mpishorthand,   only : mpicom
#endif
!
! !PUBLIC MEMBER FUNCTIONS:
  implicit none
  public :: lnd_init_mct
  public :: lnd_run_mct
  public :: lnd_final_mct
  SAVE
  private                              ! By default make data private
!
! ! PUBLIC DATA:
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
! !PRIVATE MEMBER FUNCTIONS:
  private :: lnd_export_mct
  private :: lnd_exportinit_mct
  private :: lnd_import_mct
  private :: lnd_SetGSMap_mct
  private :: lnd_CheckGrid_mct         ! check consistency of cam/clm grid
!
! !PRIVATE VARIABLES
  integer, dimension(:), allocatable :: perm  ! permutation array to reorder points
#if !(defined SPMD)
    integer, parameter :: mpicom = 1
#endif
!---------------------------------------------------------------------------

!===============================================================
contains
!===============================================================

!---------------------------------------------------------------------------
!BOP
!
! !IROUTINE: lnd_init_mct
!
! !INTERFACE:
  subroutine lnd_init_mct( gsMap_lnd, x2l_l, l2x_l, land_present )
!
! !DESCRIPTION:
! Initialize land surface model and obtain relevant atmospheric model arrays
! back from (i.e. albedos, surface temperature and snow cover over land).
!
! !USES:
    use radiation       , only : radiation_get      !(cam use)
    use filenames       , only : mss_irt, caseid    !(cam use) 
    use history         , only : ctitle, inithist   !(cam use)
    use time_manager    , only : get_nstep      
    use clm_atmlnd      , only : clm_l2a, atm_l2a
    use clm_atmlnd      , only : gridmap_l2a, clm_mapl2a
    use domainMod       , only : adomain
    use clm_comp        , only : clm_init0, clm_init1, clm_init2
    use clm_varctl      , only : cam_caseid, cam_ctitle, cam_irad, cam_nsrest, &
                                 cam_crtinic, cam_irt, finidat       
#include <comctl.h>
#include <comsol.h>
!
! !ARGUMENTS:
    type(mct_gsMap), intent(inout) :: GSMap_lnd
    type(mct_aVect), intent(inout) :: x2l_l, l2x_l
    logical        , intent(inout) :: land_present
!
! !LOCAL VARIABLES:
    integer  :: i,j         ! indices
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
!-----------------------------------------------------------------------

    !=============================================================
    ! Determine if must return now
    !=============================================================

    if (adiabatic .or. ideal_phys .or. aqua_planet) then
       land_present = .false.
       return
    end if

    !=============================================================
    ! Preset clm namelist variables and orbital params same as cam
    !=============================================================

    call radiation_get(iradsw_out=cam_irad)
    cam_caseid  = caseid
    cam_ctitle  = ctitle
    cam_nsrest  = nsrest
    cam_crtinic = inithist
    cam_irt     = mss_irt
    call lnd_setorb_mct( eccen, obliqr, lambm0, mvelpp )

    !=============================================================
    ! Initialize clm phase 1 - read namelist, grid and surface data
    !=============================================================

    call clm_init0()

#if ( defined SCAM )
    if (adomain%frac(1,1)==0) then
       land_present = .false.
       return
    end if
#endif

    ! Determine consistency with cam grid info 
    ! only need to do this on  master processor (note that cam latitudes
    ! and longitudes are computed each time by the cam model at startup)

    call lnd_CheckGrid_mct( )

    !=============================================================
    ! Initialize clm phase 2 - rest of initialization
    !=============================================================

    call clm_init1()
    call clm_init2()

    !=============================================================
    ! Initialize MCT attribute vectors and indices
    !=============================================================

    call lnd_SetGSMap_mct( GSMap_lnd ) 	

    call mct_aVect_init(x2l_l, rList=seq_flds_x2l_fields,    &
         lsize=MCT_GSMap_lsize(GSMap_lnd, mpicom))
    call mct_aVect_zero(x2l_l)

    call mct_aVect_init(l2x_l, rList=seq_flds_l2x_fields,    &
         lsize=MCT_GSMap_lsize(GSMap_lnd, mpicom))
    call mct_aVect_zero(l2x_l)

    !=============================================================
    ! Map internal data structure into coupling data structure
    !=============================================================

    if (get_nstep() == 0) then
       call clm_mapl2a( clm_l2a, atm_l2a, gridmap_l2a )
       call lnd_exportinit_mct( atm_l2a, l2x_l)
    endif

  end subroutine lnd_init_mct

!---------------------------------------------------------------------------
!BOP
!
! !IROUTINE: lnd_run_mct
!
! !INTERFACE:
  subroutine lnd_run_mct( x2l_l, l2x_l, rstwr )
!
! !DESCRIPTION:
! Run clm model
!
! !USES:
    use clm_atmlnd, only : clm_a2l, clm_l2a, atm_a2l, atm_l2a, &
                           gridmap_l2a, clm_mapl2a, gridmap_a2l, clm_mapa2l
    use clm_comp  , only : clm_run1, clm_run2
!
! !ARGUMENTS:
    type(mct_aVect), intent(inout) :: x2l_l
    type(mct_aVect), intent(inout) :: l2x_l
    logical        , intent(in)    :: rstwr    ! true => write restart file this step
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
!---------------------------------------------------------------------------

    ! Map MCT to land data type

    call lnd_import_mct( x2l_l, atm_a2l )
    call clm_mapa2l( atm_a2l, clm_a2l, gridmap_a2l )
    
    ! Run clm

    call clm_run1( )
    call clm_run2( rstwr )

    ! Map land data type to MCT

    call clm_mapl2a( clm_l2a, atm_l2a, gridmap_l2a )
    call lnd_export_mct( atm_l2a, l2x_l )

  end subroutine lnd_run_mct

!---------------------------------------------------------------------------
!BOP
!
! !IROUTINE: lnd_final_mct
!
! !INTERFACE:
  subroutine lnd_final_mct( )
!
! !DESCRIPTION:
! Finalize land surface model
!
!------------------------------------------------------------------------------
!BOP
!
! !ARGUMENTS:
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
!---------------------------------------------------------------------------

   ! fill this in
  end subroutine lnd_final_mct

!=================================================================================

  subroutine lnd_SetGSMap_mct( gsMap_lnd )

    !-------------------------------------------------------------------
    use decompMod     , only : get_proc_bounds, adecomp
    use domainMod     , only : adomain
    use clmtype

    implicit none
    type(mct_gsMap), intent(inout) :: gsMap_lnd

    integer,allocatable :: gindex(:)
    integer :: i, j, n, gi
    integer :: lsize,gsize
    integer :: ier
    integer :: begp, endp    ! beginning and ending pft indices
    integer :: begc, endc    ! beginning and ending column indices
    integer :: begl, endl    ! beginning and ending landunit indices
    integer :: begg, endg    ! beginning and ending gridcell indices
    !-------------------------------------------------------------------

    ! Build the land grid numbering for MCT
    ! NOTE:  Numbering scheme is: West to East and South to North
    ! starting at south pole.  Should be the same as what's used
    ! in SCRIP
    
    call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)

    allocate(gindex(begg:endg),stat=ier)

    ! number the local grid

    do n = begg, endg
        i = adecomp%gdc2i(n)
        j = adecomp%gdc2j(n)
        gindex(n) = (j-1)*adomain%ni + i
    end do

    lsize = endg-begg+1
    gsize = adomain%ni*adomain%nj

    allocate(perm(lsize),stat=ier)

    ! reorder gindex to be in ascending order.

    ! initialize a permutation array

    call mct_indexset(perm)

    ! derive a permutation that puts gindex in ascending order
    ! the default for IndexSort is Ascending.

    call mct_indexsort(lsize,perm,gindex)

    ! Sort gindex in-place

    call mct_permute(gindex,perm,lsize)

    call seq_mct_init_SetgsMap( lsize, gsize, gindex, mpicom, LNDID, gsMap_lnd )

    deallocate(gindex)

  end subroutine lnd_SetGSMap_mct

!=================================================================================

  subroutine lnd_exportinit_mct( l2a, l2x_l )

    !-----------------------------------------------------
    use clm_atmlnd, only : lnd2atm_type
    use domainMod , only : adomain
    use decompMod , only : get_proc_bounds_atm, adecomp
    use clm_varcon, only : sb

    implicit none

    type(lnd2atm_type), intent(inout) :: l2a
    type(mct_aVect)   , intent(inout) :: l2x_l

    integer :: g,i
    integer :: begg, endg    ! beginning and ending gridcell indices
    !-----------------------------------------------------

    call get_proc_bounds_atm(begg, endg)

!dir$ concurrent
    do g = begg,endg
       i = 1 + (g - begg)
       l2x_l%rAttr(index_l2x_Sl_landfrac,i) =  adomain%frac(adecomp%gdc2i(g), adecomp%gdc2j(g))
       l2x_l%rAttr(index_l2x_Sl_t,i)        =  sqrt(sqrt(l2a%eflx_lwrad_out(g)/sb))
       l2x_l%rAttr(index_l2x_Sl_snowh,i)    =  l2a%h2osno(g)
       l2x_l%rAttr(index_l2x_Sl_avsdr,i)    =  l2a%albd(g,1)
       l2x_l%rAttr(index_l2x_Sl_anidr,i)    =  l2a%albd(g,2)
       l2x_l%rAttr(index_l2x_Sl_avsdf,i)    =  l2a%albi(g,1)
       l2x_l%rAttr(index_l2x_Sl_anidf,i)    =  l2a%albi(g,2)
       l2x_l%rAttr(index_l2x_Fall_lwup,i)   =  l2a%eflx_lwrad_out(g)
    end do

    ! permute before using the Rearrange call

    call mct_Avect_permute(l2x_l,perm)
    
  end subroutine lnd_exportinit_mct

!====================================================================================

  subroutine lnd_export_mct( l2a, l2x_l )   

    !-----------------------------------------------------
    use time_manager, only : get_nstep
    use clm_atmlnd  , only : lnd2atm_type
    use domainMod   , only : adomain
    use decompMod   , only : get_proc_bounds_atm, adecomp
    implicit none

    type(lnd2atm_type), intent(inout) :: l2a
    type(mct_aVect)   , intent(inout) :: l2x_l

    integer :: g,i
    integer :: begg, endg    ! beginning and ending gridcell indices
    !-----------------------------------------------------
    
    call get_proc_bounds_atm(begg, endg)

    l2x_l%rAttr(:,:) = 0.0_r8

!dir$ concurrent
    do g = begg,endg
       i = 1 + (g-begg)
       l2x_l%rAttr(index_l2x_Sl_landfrac,i) = adomain%frac(adecomp%gdc2i(g), adecomp%gdc2j(g))
       l2x_l%rAttr(index_l2x_Sl_t,i)        = l2a%t_rad(g)
       l2x_l%rAttr(index_l2x_Sl_snowh,i)    = l2a%h2osno(g)
       l2x_l%rAttr(index_l2x_Sl_avsdr,i)    = l2a%albd(g,1)
       l2x_l%rAttr(index_l2x_Sl_anidr,i)    = l2a%albd(g,2)
       l2x_l%rAttr(index_l2x_Sl_avsdf,i)    = l2a%albi(g,1)
       l2x_l%rAttr(index_l2x_Sl_anidf,i)    = l2a%albi(g,2)
       l2x_l%rAttr(index_l2x_Sl_tref,i)     = l2a%t_ref2m(g)
       l2x_l%rAttr(index_l2x_Sl_qref,i)     = l2a%q_ref2m(g)
       l2x_l%rAttr(index_l2x_Fall_taux,i)   = l2a%taux(g)
       l2x_l%rAttr(index_l2x_Fall_tauy,i)   = l2a%tauy(g)
       l2x_l%rAttr(index_l2x_Fall_lat,i)    = l2a%eflx_lh_tot(g)
       l2x_l%rAttr(index_l2x_Fall_sen,i)    = l2a%eflx_sh_tot(g)
       l2x_l%rAttr(index_l2x_Fall_lwup,i)   = l2a%eflx_lwrad_out(g)
       l2x_l%rAttr(index_l2x_Fall_evap,i)   = l2a%qflx_evap_tot(g)
       l2x_l%rAttr(index_l2x_Fall_swnet,i)  = l2a%fsa(g)
       if (index_l2x_Fall_nee /= 0) then
          l2x_l%rAttr(index_l2x_Fall_nee,i) = l2a%nee(g)
       end if
    end do

    ! permute before using the Rearrange call.

    call mct_aVect_permute(l2x_l,perm)

  end subroutine lnd_export_mct

!====================================================================================

  subroutine lnd_import_mct( x2l_l, a2l )

    !-----------------------------------------------------
    use clm_atmlnd      , only: atm2lnd_type
    use clm_varctl      , only: co2_type
    use clm_varcon      , only: rair, o2_molar_const, co2_ppmv_const, c13ratio
    use decompMod       , only: get_proc_bounds_atm

    implicit none
    type(mct_aVect)   , intent(inout) :: x2l_l
    type(atm2lnd_type), intent(inout) :: a2l

    integer  :: g,i,nstep,ier
    real(r8) :: forc_rainc    ! rainxy Atm flux mm/s
    real(r8) :: forc_rainl    ! rainxy Atm flux mm/s
    real(r8) :: forc_snowc    ! snowfxy Atm flux  mm/s
    real(r8) :: forc_snowl    ! snowfxl Atm flux  mm/s
    real(r8) :: co2_ppmv_diag ! temporary
    real(r8) :: co2_ppmv_prog ! temporary
    real(r8) :: co2_ppmv      ! temporary
    integer  :: begg, endg    ! beginning and ending gridcell indices
    integer  :: co2_type_idx  ! integer flag for co2_type options
    !-----------------------------------------------------

    ! unpermute after rearrange call and before copying into local arrays.

    call mct_aVect_unpermute(x2l_l, perm)

    call get_proc_bounds_atm(begg, endg)

    co2_type_idx = 0
    if (co2_type == 'prognostic') then
       co2_type_idx = 1
    else if (co2_type == 'diagnostic') then
       co2_type_idx = 2
    end if
    if (co2_type == 'prognostic' .and. index_x2l_Sa_co2prog == 0) then
       write(6,*)' must have nonzero index_x2l_Sa_co2prog for co2_type equal to prognostic'
       call endrun()
    else if (co2_type == 'diagnostic' .and. index_x2l_Sa_co2diag == 0) then
       write(6,*)' must have nonzero index_x2l_Sa_co2diag for co2_type equal to diagnostic'
       call endrun()
    end if
    
    ! Note that the precipitation fluxes received  from the coupler
    ! are in units of kg/s/m^2. To convert these precipitation rates
    ! in units of mm/sec, one must divide by 1000 kg/m^3 and multiply
    ! by 1000 mm/m resulting in an overall factor of unity.
    ! Below the units are therefore given in mm/s.
    
!dir$ concurrent
    do g = begg,endg
        i = 1 + (g - begg)
       
        ! Determine required receive fields

        a2l%forc_hgt(g)     = x2l_l%rAttr(index_x2l_Sa_z,i)         ! zgcmxy  Atm state m
        a2l%forc_u(g)       = x2l_l%rAttr(index_x2l_Sa_u,i)         ! forc_uxy  Atm state m/s
        a2l%forc_v(g)       = x2l_l%rAttr(index_x2l_Sa_v,i)         ! forc_vxy  Atm state m/s
        a2l%forc_th(g)      = x2l_l%rAttr(index_x2l_Sa_ptem,i)      ! forc_thxy Atm state K
        a2l%forc_q(g)       = x2l_l%rAttr(index_x2l_Sa_shum,i)      ! forc_qxy  Atm state kg/kg
        a2l%forc_pbot(g)    = x2l_l%rAttr(index_x2l_Sa_pbot,i)      ! ptcmxy  Atm state Pa
        a2l%forc_t(g)       = x2l_l%rAttr(index_x2l_Sa_tbot,i)      ! forc_txy  Atm state K
        a2l%forc_lwrad(g)   = x2l_l%rAttr(index_x2l_Faxa_lwdn,i)    ! flwdsxy Atm flux  W/m^2
        forc_rainc          = x2l_l%rAttr(index_x2l_Faxa_rainc,i)   ! mm/s
        forc_rainl          = x2l_l%rAttr(index_x2l_Faxa_rainl,i)   ! mm/s
        forc_snowc          = x2l_l%rAttr(index_x2l_Faxa_snowc,i)   ! mm/s
        forc_snowl          = x2l_l%rAttr(index_x2l_Faxa_snowl,i)   ! mm/s
        a2l%forc_solad(g,2) = x2l_l%rAttr(index_x2l_Faxa_swndr,i)   ! forc_sollxy  Atm flux  W/m^2
        a2l%forc_solad(g,1) = x2l_l%rAttr(index_x2l_Faxa_swvdr,i)   ! forc_solsxy  Atm flux  W/m^2
        a2l%forc_solai(g,2) = x2l_l%rAttr(index_x2l_Faxa_swndf,i)   ! forc_solldxy Atm flux  W/m^2
        a2l%forc_solai(g,1) = x2l_l%rAttr(index_x2l_Faxa_swvdf,i)   ! forc_solsdxy Atm flux  W/m^2

        ! Determine optional receive fields

        if (index_x2l_Sa_co2prog /= 0) then
           co2_ppmv_prog = x2l_l%rAttr(index_x2l_Sa_co2prog,i)   ! co2 atm state prognostic
        else
           co2_ppmv_prog = co2_ppmv_const
        end if
 
        if (index_x2l_Sa_co2diag /= 0) then
           co2_ppmv_diag = x2l_l%rAttr(index_x2l_Sa_co2diag,i)   ! co2 atm state diagnostic
        else
           co2_ppmv_diag = co2_ppmv_const
        end if

        ! Determine derived quantities for required fields

        a2l%forc_hgt_u(g) = a2l%forc_hgt(g)    !observational height of wind [m]
        a2l%forc_hgt_t(g) = a2l%forc_hgt(g)    !observational height of temperature [m]
        a2l%forc_hgt_q(g) = a2l%forc_hgt(g)    !observational height of humidity [m]
        a2l%forc_vp(g)    = a2l%forc_q(g) * a2l%forc_pbot(g) &
                            / (0.622_r8 + 0.378_r8 * a2l%forc_q(g))
        a2l%forc_rho(g)   = (a2l%forc_pbot(g) - 0.378_r8 * a2l%forc_vp(g)) &
                            / (rair * a2l%forc_t(g))
        a2l%forc_po2(g)   = o2_molar_const * a2l%forc_pbot(g)
        a2l%forc_wind(g)  = sqrt(a2l%forc_u(g)**2 + a2l%forc_v(g)**2)
        a2l%forc_solar(g) = a2l%forc_solad(g,1) + a2l%forc_solai(g,1) + &
                            a2l%forc_solad(g,2) + a2l%forc_solai(g,2)
        a2l%forc_rain(g)  = forc_rainc + forc_rainl
        a2l%forc_snow(g)  = forc_snowc + forc_snowl
        
        ! Determine derived quantities for optional fields
        ! Note that the following does unit conversions from ppmv to partial pressures (Pa)
        ! Note that forc_pbot is in Pa

        if (co2_type_idx == 1) then
           co2_ppmv = co2_ppmv_prog
        else if (co2_type_idx == 2) then
           co2_ppmv = co2_ppmv_diag 
        else
           co2_ppmv = co2_ppmv_const      
        end if
        a2l%forc_pco2(g)   = co2_ppmv * 1.e-6_r8 * a2l%forc_pbot(g) 
        a2l%forc_pc13o2(g) = co2_ppmv * c13ratio * 1.e-6_r8 * a2l%forc_pbot(g)
	 
     end do

   end subroutine lnd_import_mct

!---------------------------------------------------------------------------
!BOP
!
! !IROUTINE: lnd_CheckGrid_mct
!
! !INTERFACE:
  subroutine lnd_CheckGrid_mct( )
!
! !DESCRIPTION:
! Check that cam grid is consistent with clm grid read in from clm surface 
! dataset
!
! !USES:
    use shr_const_mod   , only : SHR_CONST_PI
    use commap          , only : clat, londeg              !(cam use)
    use comsrf          , only : ext_frac => landfrac_glob !(cam use)
    use domainMod       , only : adomain 
    use spmdMod         , only : masterproc
!
! !ARGUMENTS:
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein 2005-05-14
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: i,j        ! indices
    real(r8):: ext_lonc   ! ext lon values
    real(r8):: ext_latc   ! ext lat values
    integer :: ext_mask   ! ext land mask values	
!---------------------------------------------------------------------------

    do j = 1,adomain%nj
       do i = 1,adomain%ni

          ! Determine consistency of cam and clm grid - all processors

          ext_latc = (180._r8/SHR_CONST_PI)*clat(j)
          if ( abs(ext_latc - adomain%latc(i,j)) > 1.e-12_r8 ) then
             write(6,*)'MCT_lnd_checkgrid error: CAM latitude ',ext_latc,' and clm input latitude ', &
                  adomain%latc(i,j),' has difference too large at i,j= ',i,j
             call endrun()
          end if
          ext_lonc = londeg(i,j)
          if ( abs(ext_lonc - adomain%lonc(i,j)) > 1.e-12_r8 ) then
             write(6,*)'MCT_lnd_checkgrid error: CAM longitude ',ext_lonc,' and clm input longitude ', &
                  adomain%lonc(i,j),' has difference too large at i,j= ',i,j
             call endrun()
          end if

          ! Determine consistency of cam and clm landfrac/landmask - masterproc only

          if (masterproc) then
             if (ext_frac(i,j) > 0._r8) then
                ext_mask = 1
             else
                ext_mask = 0
             endif
             if (ext_mask /= adomain%mask(i,j)) then
                write(6,*)'MCT_lnd_checkgrid error: CAM land mask different from surface dataset at i,j= ',i,j
                call endrun()
             end if
             if (ext_frac(i,j) /= adomain%frac(i,j)) then
                write(6,*)'MCT_lnd_checkgrid error: CAM fractional land differs from surface dataset at i,j= ',i,j
                call endrun()
             end if
          end if
    
       end do
    end do
       
  end subroutine lnd_CheckGrid_mct

!---------------------------------------------------------------------------
!BOP
!
! !IROUTINE: lnd_setorb_mct_mct
!
! !INTERFACE:
  subroutine lnd_setorb_mct(clm_eccen, clm_obliqr, clm_lambm0, clm_mvelpp)
!
! !DESCRIPTION:
! Determine clm orbital parameters
!
! !USES:
    use clm_varorb, only : eccen, obliqr, lambm0, mvelpp
!
! !ARGUMENTS: 
    implicit none
    real(r8), intent(in) :: clm_eccen
    real(r8), intent(in) :: clm_obliqr 
    real(r8), intent(in) :: clm_lambm0
    real(r8), intent(in) :: clm_mvelpp
!
!EOP
!-----------------------------------------------------------------------

    eccen  = clm_eccen  
    obliqr = clm_obliqr 
    lambm0 = clm_lambm0
    mvelpp = clm_mvelpp

  end subroutine lnd_setorb_mct

#endif 

end module lnd_comp_mct

