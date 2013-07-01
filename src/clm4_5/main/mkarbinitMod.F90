module mkarbinitMod
!---------------------------------------------------------------------------
!BOP
!
! !MODULE: mkarbinitMod
!
! !DESCRIPTION:
!
!
!---------------------------------------------------------------------------

! !USES:
    use shr_kind_mod , only : r8 => shr_kind_r8
    use clm_varctl   , only : iulog
    use shr_sys_mod  , only : shr_sys_flush
    use spmdMod      , only : masterproc

    implicit none

    SAVE
    private                              ! By default make data private

! !PUBLIC MEMBER FUNCTIONS:

    public mkarbinit   ! Make arbitrary initial conditions

!EOP
!-----------------------------------------------------------------------
contains
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: mkarbinit
!
! !INTERFACE:
  subroutine mkarbinit()
!
! !DESCRIPTION:
! Initializes the following time varying variables:
! water      : h2osno, h2ocan, h2osoi_liq, h2osoi_ice, h2osoi_vol
! snow       : snow_depth, snl, dz, z, zi
! temperature: t_soisno, t_veg, t_grnd
!
! !USES:
    use shr_const_mod, only : SHR_CONST_TKFRZ
    use clmtype
    use clm_varpar   , only : nlevsoi, nlevgrnd, nlevsno, nlevlak, nlevurb
    use clm_varcon   , only : bdsno, istice, istwet, istsoil, &
                              denice, denh2o, spval, sb, icol_road_perv, &
                              icol_road_imperv, icol_roof, icol_sunwall, &
                              icol_shadewall
    use clm_varcon   , only : istcrop
    use clm_varcon   , only : istice_mec, h2osno_max
    use clm_varctl   , only : iulog
    use spmdMod      , only : masterproc
    use decompMod    , only : get_proc_bounds
    use SNICARMod    , only : snw_rds_min

!
! !ARGUMENTS:
    implicit none
!
! !CALLED FROM:
! subroutine initialize in module initializeMod
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
! 3/07/08 Keith Oleson: initialize h2osoi_vol for all soil layers to 0.3
! 3/18/08 David Lawrence, initialize deep layers
! 03/28/08 Mark Flanner, initialize snow aerosols and grain size
!
! !LOCAL VARIABLES:
!
!
!
!

!
!
!EOP
    integer :: j,l,c,p      ! indices
    integer :: nlevs        ! number of levels
    integer :: begp, endp   ! per-proc beginning and ending pft indices
    integer :: begc, endc   ! per-proc beginning and ending column indices
    integer :: begl, endl   ! per-proc beginning and ending landunit indices
    integer :: begg, endg   ! per-proc gridcell ending gridcell indices
    real(r8):: vwc,psi      ! for calculating soilpsi
!-----------------------------------------------------------------------

    if ( masterproc )then
        write(iulog,*) 'Setting initial data to non-spun up values'
    end if


   associate(& 
   h2osfc                              =>    cws%h2osfc                                  , & ! Input:  [real(r8) (:)]  surface water (mm)                                
   t_h2osfc                            =>    ces%t_h2osfc                                , & ! Input:  [real(r8) (:)]  surface water temperature                         
   frac_h2osfc                         =>    cps%frac_h2osfc                             , & ! Input:  [real(r8) (:)]  fraction of ground covered by surface water (0 to 1)
   qflx_h2osfc_surf                    =>    cwf%qflx_h2osfc_surf                        , & ! Input:  [real(r8) (:)] surface water runoff (mm/s)                        
   qflx_snow_melt                      =>    cwf%qflx_snow_melt                          , & ! Input:  [real(r8) (:)]  snow melt (net)                                   
   frost_table                         =>    cws%frost_table                             , & ! Input:  [real(r8) (:)]  frost table depth (m)                             
   zwt_perched                         =>    cws%zwt_perched                             , & ! Input:  [real(r8) (:)]  perched water table depth (m)                     
   int_snow                            =>    cws%int_snow                                , & ! Input:  [real(r8) (:)]  integrated snowfall                               
   ltype                               =>    lun%itype                                   , & ! Input:  [integer (:)]  landunit type                                      
   lakpoi                              =>   lun%lakpoi                                   , & ! Input:  [logical (:)]  true => landunit is a lake point                   
   ifspecial                           =>   lun%ifspecial                                , & ! Input:  [logical (:)]  true => landunit is not vegetated                  
   urbpoi                              =>   lun%urbpoi                                   , & ! Input:  [logical (:)]  true => landunit is an urban point                 
   ctype                               =>    col%itype                                   , & ! Input:  [integer (:)]  column type                                        
   clandunit                           =>   col%landunit                                 , & ! Input:  [integer (:)]  landunit index associated with each column         
   snl                                 =>    cps%snl                                     , & ! Output: [integer (:)]  number of snow layers                              
   dz                                  =>    cps%dz                                      , & ! Input:  [real(r8) (:,:)]  layer thickness depth (m)                       
   watsat                              =>    cps%watsat                                  , & ! Input:  [real(r8) (:,:)]  volumetric soil water at saturation (porosity)  
   sucsat                              =>    cps%sucsat                                  , & ! Input:  [real(r8) (:,:)]  minimum soil suction (mm)                       
   bsw                                 =>    cps%bsw                                     , & ! Input:  [real(r8) (:,:)]  Clapp and Hornberger "b"                        
   soilpsi                             =>    cps%soilpsi                                 , & ! Output: [real(r8) (:,:)]  soil water potential in each soil layer (MPa)   
   h2osoi_ice                          =>    cws%h2osoi_ice                              , & ! Input:  [real(r8) (:,:)]  ice lens (kg/m2)                                
   h2osoi_liq                          =>    cws%h2osoi_liq                              , & ! Input:  [real(r8) (:,:)]  liquid water (kg/m2)                            
   h2osoi_vol                          =>    cws%h2osoi_vol                              , & ! Output: [real(r8) (:,:)]  volumetric soil water (0<=h2osoi_vol<=watsat) [m3/m3]
   h2ocan_col                          =>    pws_a%h2ocan                                , & ! Output: [real(r8) (:)]  canopy water (mm H2O) (column-level)              
   qflx_irrig                          =>    cwf%qflx_irrig                              , & ! Output: [real(r8) (:)]  irrigation flux (mm H2O/s)                        
   snow_depth                          =>    cps%snow_depth                              , & ! Output: [real(r8) (:)]  snow height (m)                                   
   h2osno                              =>    cws%h2osno                                  , & ! Output: [real(r8) (:)]  snow water (mm H2O)                               
   t_soisno                            =>    ces%t_soisno                                , & ! Output: [real(r8) (:,:)]  soil temperature (Kelvin)  (-nlevsno+1:nlevgrnd)
   t_lake                              =>    ces%t_lake                                  , & ! Output: [real(r8) (:,:)]  lake temperature (Kelvin)  (1:nlevlak)          
   t_grnd                              =>    ces%t_grnd                                  , & ! Output: [real(r8) (:)]  ground temperature (Kelvin)                       
   tsoi17                              =>    ces%tsoi17                                  , & ! Output: [real(r8) (:)]  soil T for top 0.17 m                             
   zi                                  =>    cps%zi                                      , & ! Input:  [real(r8) (:,:)]  interface level below a "z" level (m)           
   wa                                  =>    cws%wa                                      , & ! Input:  [real(r8) (:)]  water in the unconfined aquifer (mm)              
   zwt                                 =>    cws%zwt                                     , & ! Input:  [real(r8) (:)]  water table depth (m)                             
   fsat                                =>    cws%fsat                                    , & ! Output: [real(r8) (:)] fractional area with water table at surface        
   snw_rds                             =>    cps%snw_rds                                 , & ! Output: [real(r8) (:,:)]  effective snow grain radius (col,lyr) [microns, m^-6]
   snw_rds_top                         =>    cps%snw_rds_top                             , & ! Output: [real(r8) (:)]  snow grain size, top (col) [microns]              
   sno_liq_top                         =>    cps%sno_liq_top                             , & ! Output: [real(r8) (:)]  liquid water fraction (mass) in top snow layer (col) [frc]
   mss_bcpho                           =>    cps%mss_bcpho                               , & ! Output: [real(r8) (:,:)]  mass of hydrophobic BC in snow (col,lyr) [kg]   
   mss_bcphi                           =>    cps%mss_bcphi                               , & ! Output: [real(r8) (:,:)]  mass of hydrophillic BC in snow (col,lyr) [kg]  
   mss_bctot                           =>    cps%mss_bctot                               , & ! Output: [real(r8) (:,:)]  total mass of BC (pho+phi) (col,lyr) [kg]       
   mss_bc_col                          =>    cps%mss_bc_col                              , & ! Output: [real(r8) (:)]  total mass of BC in snow column (col) [kg]        
   mss_bc_top                          =>    cps%mss_bc_top                              , & ! Output: [real(r8) (:)]  total mass of BC in top snow layer (col) [kg]     
   mss_cnc_bcphi                       =>    cps%mss_cnc_bcphi                           , & ! Output: [real(r8) (:,:)]  mass concentration of BC species 1 (col,lyr) [kg/kg]
   mss_cnc_bcpho                       =>    cps%mss_cnc_bcpho                           , & ! Output: [real(r8) (:,:)]  mass concentration of BC species 2 (col,lyr) [kg/kg]
   mss_ocpho                           =>    cps%mss_ocpho                               , & ! Output: [real(r8) (:,:)]  mass of hydrophobic OC in snow (col,lyr) [kg]   
   mss_ocphi                           =>    cps%mss_ocphi                               , & ! Output: [real(r8) (:,:)]  mass of hydrophillic OC in snow (col,lyr) [kg]  
   mss_octot                           =>    cps%mss_octot                               , & ! Output: [real(r8) (:,:)]  total mass of OC (pho+phi) (col,lyr) [kg]       
   mss_oc_col                          =>    cps%mss_oc_col                              , & ! Output: [real(r8) (:)]  total mass of OC in snow column (col) [kg]        
   mss_oc_top                          =>    cps%mss_oc_top                              , & ! Output: [real(r8) (:)]  total mass of OC in top snow layer (col) [kg]     
   mss_cnc_ocphi                       =>    cps%mss_cnc_ocphi                           , & ! Output: [real(r8) (:,:)]  mass concentration of OC species 1 (col,lyr) [kg/kg]
   mss_cnc_ocpho                       =>    cps%mss_cnc_ocpho                           , & ! Output: [real(r8) (:,:)]  mass concentration of OC species 2 (col,lyr) [kg/kg]
   mss_dst1                            =>    cps%mss_dst1                                , & ! Output: [real(r8) (:,:)]  mass of dust species 1 in snow (col,lyr) [kg]   
   mss_dst2                            =>    cps%mss_dst2                                , & ! Output: [real(r8) (:,:)]  mass of dust species 2 in snow (col,lyr) [kg]   
   mss_dst3                            =>    cps%mss_dst3                                , & ! Output: [real(r8) (:,:)]  mass of dust species 3 in snow (col,lyr) [kg]   
   mss_dst4                            =>    cps%mss_dst4                                , & ! Output: [real(r8) (:,:)]  mass of dust species 4 in snow (col,lyr) [kg]   
   mss_dsttot                          =>    cps%mss_dsttot                              , & ! Output: [real(r8) (:,:)]  total mass of dust in snow (col,lyr) [kg]       
   mss_dst_col                         =>    cps%mss_dst_col                             , & ! Output: [real(r8) (:)]  total mass of dust in snow column (col) [kg]      
   mss_dst_top                         =>    cps%mss_dst_top                             , & ! Output: [real(r8) (:)]  total mass of dust in top snow layer (col) [kg]   
   mss_cnc_dst1                        =>    cps%mss_cnc_dst1                            , & ! Output: [real(r8) (:,:)]  mass concentration of dust species 1 (col,lyr) [kg/kg]
   mss_cnc_dst2                        =>    cps%mss_cnc_dst2                            , & ! Output: [real(r8) (:,:)]  mass concentration of dust species 2 (col,lyr) [kg/kg]
   mss_cnc_dst3                        =>    cps%mss_cnc_dst3                            , & ! Output: [real(r8) (:,:)]  mass concentration of dust species 3 (col,lyr) [kg/kg]
   mss_cnc_dst4                        =>    cps%mss_cnc_dst4                            , & ! Output: [real(r8) (:,:)]  mass concentration of dust species 4 (col,lyr) [kg/kg]
   n_irrig_steps_left                  =>    cps%n_irrig_steps_left                      , & ! Output: [integer (:)]  number of time steps for which we still need to irrigate today (if 0, ignore irrig_rate)
   irrig_rate                          =>    cps%irrig_rate                              , & ! Output: [real(r8) (:)]  current irrigation rate [mm/s]                    
   pcolumn                             =>   pft%column                                   , & ! Input:  [integer (:)]  column index associated with each pft              
   h2ocan_pft                          =>    pws%h2ocan                                  , & ! Output: [real(r8) (:)]  canopy water (mm H2O) (pft-level)                 
   t_veg                               =>    pes%t_veg                                   , & ! Output: [real(r8) (:)]  vegetation temperature (Kelvin)                   
   t_ref2m                             =>    pes%t_ref2m                                 , & ! Output: [real(r8) (:)]  2 m height surface air temperature (Kelvin)       
   t_ref2m_u                           =>    pes%t_ref2m_u                               , & ! Output: [real(r8) (:)]  Urban 2 m height surface air temperature (Kelvin) 
   t_ref2m_r                           =>    pes%t_ref2m_r                               , & ! Output: [real(r8) (:)]  Rural 2 m height surface air temperature (Kelvin) 
   plandunit                           =>   pft%landunit                                 , & ! Input:  [integer (:)]  landunit index associated with each pft            
   eflx_lwrad_out                      =>    pef%eflx_lwrad_out                            & ! Output: [real(r8) (:)]  emitted infrared (longwave) radiation (W/m**2)    
   )

    ! Determine subgrid bounds on this processor

    call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)

    ! NOTE: h2ocan, h2osno, and snow_depth has valid values everywhere
    ! canopy water (pft level)

    do p = begp, endp
       h2ocan_pft(p) = 0._r8
    end do

    ! initialize h2osfc, frac_h2osfc, t_h2osfc, qflx_snow_melt
    do c = begc,endc
       h2osfc(c)           = 0._r8
       frac_h2osfc(c)      = 0._r8
       !       t_h2osfc(c) = spval
       t_h2osfc(c)         = 274._r8
       qflx_h2osfc_surf(c) = 0._r8
       qflx_snow_melt(c)   = 0._r8
    enddo

    do c = begc,endc

       ! canopy water (column level)

       h2ocan_col(c) = 0._r8

       ! snow water

       l = clandunit(c)

       ! Note: Glacier_mec columns are initialized with half the maximum snow cover.
       ! This gives more realistic values of qflx_glcice sooner in the simulation
       !  for columns with net ablation, at the cost of delaying ice formation
       !  in columns with net accumulation.
       if (ltype(l)==istice) then
          h2osno(c) = h2osno_max
       elseif (ltype(l)==istice_mec) then
          h2osno(c) = 0.5_r8 * h2osno_max   ! 50 cm if h2osno_max = 1 m
       else
          h2osno(c) = 0._r8
       endif

       ! initialize int_snow, int_melt
       int_snow(c) = h2osno(c)
       ! snow depth

       snow_depth(c)  = h2osno(c) / bdsno

       ! Initialize Irrigation to zero
       if (ltype(l)==istsoil) then
          n_irrig_steps_left(c) = 0
          irrig_rate(c)         = 0.0_r8
       end if

    end do

    ! Set snow layer number, depth and thickiness

    call snow_depth2lev(begc, endc)

    ! Set snow/soil temperature, note:
    ! t_soisno only has valid values over non-lake
    ! t_lake   only has valid values over lake
    ! t_grnd has valid values over all land
    ! t_veg  has valid values over all land

    ! NOTE: THESE MEMORY COPIES ARE INEFFICIENT -- SINCE nlev LOOP IS NESTED FIRST!!!!
    do c = begc,endc

       t_soisno(c,-nlevsno+1:nlevgrnd) = spval
       t_lake(c,1:nlevlak) = spval

       l = clandunit(c)
       if (.not. lakpoi(l)) then  !not lake
          t_soisno(c,-nlevsno+1:0) = spval
          if (snl(c) < 0) then    !snow layer temperatures
             do j = snl(c)+1, 0
                t_soisno(c,j) = 250._r8
             enddo
          endif
          if (ltype(l)==istice .or. ltype(l)==istice_mec) then
             do j = 1, nlevgrnd
                t_soisno(c,j) = 250._r8
             end do
          else if (ltype(l) == istwet) then
             do j = 1, nlevgrnd
                t_soisno(c,j) = 277._r8
             end do
          else if (urbpoi(l)) then
#if (defined VANCOUVER)
             if (ctype(c) == icol_road_perv .or. ctype(c) == icol_road_imperv) then 
               ! Set road top layer to initial air temperature and interpolate other
               ! layers down to 20C in bottom layer
               do j = 1, nlevgrnd
                  t_soisno(c,j) = 297.56 - (j-1) * ((297.56-293.16)/(nlevgrnd-1)) 
               end do
             ! Set wall and roof layers to initial air temperature
             else if (ctype(c) == icol_sunwall .or. ctype(c) == icol_shadewall .or. ctype(c) == icol_roof) then
               do j = 1, nlevurb
                  t_soisno(c,j) = 297.56
               end do
             else
               do j = 1, nlevgrnd
                  t_soisno(c,j) = 283._r8
               end do
             end if
#elif (defined MEXICOCITY)
             if (ctype(c) == icol_road_perv .or. ctype(c) == icol_road_imperv) then 
               ! Set road top layer to initial air temperature and interpolate other
               ! layers down to 22C in bottom layer
               do j = 1, nlevgrnd
                  t_soisno(c,j) = 289.46 - (j-1) * ((289.46-295.16)/(nlevgrnd-1)) 
               end do
             ! Set wall and roof layers to initial air temperature
             else if (ctype(c) == icol_sunwall .or. ctype(c) == icol_shadewall .or. ctype(c) == icol_roof) then
               do j = 1, nlevurb
                  t_soisno(c,j) = 289.46
               end do
             else
               do j = 1, nlevgrnd
                  t_soisno(c,j) = 283._r8
               end do
             end if
#else
             if (ctype(c) == icol_road_perv .or. ctype(c) == icol_road_imperv) then 
               do j = 1, nlevgrnd
                 t_soisno(c,j) = 274._r8
               end do
             ! Set sunwall, shadewall, roof to fairly high temperature to avoid initialization
             ! shock from large heating/air conditioning flux
             else if (ctype(c) == icol_sunwall .or. ctype(c) == icol_shadewall &
                      .or. ctype(c) == icol_roof) then
               do j = 1, nlevurb
                 t_soisno(c,j) = 292._r8
               end do
             end if
#endif
          else
             do j = 1, nlevgrnd
                t_soisno(c,j) = 274._r8
             end do
          endif
          t_grnd(c) = t_soisno(c,snl(c)+1)
       else                     !lake
          t_lake(c,1:nlevlak) = 277._r8
          t_grnd(c) = t_lake(c,1)
       endif
       tsoi17(c) = t_grnd(c)

    end do

    do p = begp, endp
       c = pcolumn(p)
       l = plandunit(p)

       ! Initialize Irrigation to zero
       if (ltype(l)==istsoil) then
          qflx_irrig(c)      = 0.0_r8
       end if

#if (defined VANCOUVER)
       t_veg(p) = 297.56
       t_ref2m(p) = 297.56
       if (urbpoi(l)) then
         t_ref2m_u(p) = 297.56
       else
         t_ref2m_u(p) = spval
       end if
       if (ifspecial(l)) then
         t_ref2m_r(p) = spval
       else
         t_ref2m_r(p) = 297.56
       end if 
#elif (defined MEXICOCITY)
       t_veg(p) = 289.46
       t_ref2m(p) = 289.46
       if (urbpoi(l)) then
         t_ref2m_u(p) = 289.46
       else
         t_ref2m_u(p) = spval
       end if
       if (ifspecial(l)) then
         t_ref2m_r(p) = spval
       else
         t_ref2m_r(p) = 289.46
       end if 
#else
       t_veg(p) = 283._r8
       t_ref2m(p) = 283._r8
       if (urbpoi(l)) then
         t_ref2m_u(p) = 283._r8
       else
         t_ref2m_u(p) = spval
       end if
       if (ifspecial(l)) then
         t_ref2m_r(p) = spval
       else
         t_ref2m_r(p) = 283._r8
       end if 
#endif
       eflx_lwrad_out(p) = sb * (t_grnd(c))**4
    end do

    ! Set snow/soil ice and liquid mass

    ! volumetric water is set first and liquid content and ice lens are obtained
    ! NOTE: h2osoi_vol, h2osoi_liq and h2osoi_ice only have valid values over soil
    ! and urban pervious road (other urban columns have zero soil water)

    h2osoi_vol(begc:endc,         1:) = spval
    h2osoi_liq(begc:endc,-nlevsno+1:) = spval
    h2osoi_ice(begc:endc,-nlevsno+1:) = spval

    wa(begc:endc)  = 5000._r8
    zwt(begc:endc) = 0._r8

    do c = begc,endc
       l = clandunit(c)
       if (.not. lakpoi(l)) then  !not lake
          if (urbpoi(l)) then
             if (ctype(c) == icol_road_perv) then
                wa(c)  = 4800._r8
                zwt(c) = (25._r8 + zi(c,nlevsoi)) - wa(c)/0.2_r8 /1000._r8  ! One meter below soil column
             else
                wa(c)  = spval
                zwt(c) = spval
             end if
             ! initialize frost_table, zwt_perched
             zwt_perched(c) = spval
             frost_table(c) = spval
          else
             wa(c)  = 4000._r8
             zwt(c) = (25._r8 + zi(c,nlevsoi)) - wa(c)/0.2_r8 /1000._r8  ! One meter below soil column
             ! initialize frost_table, zwt_perched to bottom of soil column
             zwt_perched(c) = zi(c,nlevsoi)
             frost_table(c) = zi(c,nlevsoi)
          end if
       end if
    end do

    do c = begc,endc
       l = clandunit(c)
       if (.not. lakpoi(l)) then  !not lake

          ! volumetric water
          if (ltype(l) == istsoil .or. ltype(l) == istcrop) then
             nlevs = nlevgrnd
             do j = 1, nlevs
                if (j > nlevsoi) then
                   h2osoi_vol(c,j) = 0.0_r8
                else
                   h2osoi_vol(c,j) = 0.15_r8
                endif
             end do
          else if (urbpoi(l)) then
             if (ctype(c) == icol_road_perv) then
               nlevs = nlevgrnd
               do j = 1, nlevs
                  if (j <= nlevsoi) then
                     h2osoi_vol(c,j) = 0.3_r8
                  else
                     h2osoi_vol(c,j) = 0.0_r8
                  end if
               end do
             else if (ctype(c) == icol_road_imperv) then
               nlevs = nlevgrnd
               do j = 1, nlevs
                  h2osoi_vol(c,j) = 0.0_r8
               end do
             else
               nlevs = nlevurb
               do j = 1, nlevs
                  h2osoi_vol(c,j) = 0.0_r8
               end do
             end if
          else if (ltype(l) == istwet) then
             nlevs = nlevgrnd
             do j = 1, nlevs
                if (j > nlevsoi) then
                   h2osoi_vol(c,j) = 0.0_r8
                else
                   h2osoi_vol(c,j) = 1.0_r8
                endif
             end do
          else if (ltype(l) == istice .or. ltype(l) == istice_mec) then
             nlevs = nlevgrnd 
             do j = 1, nlevs
                h2osoi_vol(c,j) = 1.0_r8
             end do
          endif
          do j = 1, nlevs
             h2osoi_vol(c,j) = min(h2osoi_vol(c,j),watsat(c,j))
        
             ! soil layers
             if (t_soisno(c,j) <= SHR_CONST_TKFRZ) then
                h2osoi_ice(c,j)  = dz(c,j)*denice*h2osoi_vol(c,j)
                h2osoi_liq(c,j) = 0._r8
             else
                h2osoi_ice(c,j) = 0._r8
                h2osoi_liq(c,j) = dz(c,j)*denh2o*h2osoi_vol(c,j)
             endif
          end do

#if (defined CN) 
          ! soil water potential (added 10/21/03, PET)
          ! required for CN code
          if (ltype(l) == istsoil .or. ltype(l) == istcrop) then
             nlevs = nlevgrnd
             do j = 1, nlevs
                if (h2osoi_liq(c,j) > 0._r8) then
                   vwc = h2osoi_liq(c,j)/(dz(c,j)*denh2o)
                   psi = sucsat(c,j) * (-9.8e-6_r8) * (vwc/watsat(c,j))**(-bsw(c,j))  ! Mpa
                   soilpsi(c,j) = max(psi, -15.0_r8)
                   soilpsi(c,j) = min(soilpsi(c,j),0.0_r8)
                end if
             end do
          end if
          fsat(c)   = 0.0_r8
#endif
       end if

    end do

    ! Set snow

    do j = -nlevsno+1, 0
       do c = begc,endc
          l = clandunit(c)
          if (.not. lakpoi(l)) then  !not lake
             if (j > snl(c)) then
                h2osoi_ice(c,j) = dz(c,j)*250._r8
                h2osoi_liq(c,j) = 0._r8
             end if
          end if
       end do
    end do


    ! initialize SNICAR fields:
    do c = begc,endc
       mss_bctot(c,:) = 0._r8
       mss_bcpho(c,:) = 0._r8
       mss_bcphi(c,:) = 0._r8
       mss_cnc_bcphi(c,:)=0._r8
       mss_cnc_bcpho(c,:)=0._r8

       mss_octot(c,:) = 0._r8
       mss_ocpho(c,:) = 0._r8
       mss_ocphi(c,:) = 0._r8
       mss_cnc_ocphi(c,:)=0._r8
       mss_cnc_ocpho(c,:)=0._r8
       
       mss_dst1(c,:) = 0._r8
       mss_dst2(c,:) = 0._r8
       mss_dst3(c,:) = 0._r8
       mss_dst4(c,:) = 0._r8
       mss_dsttot(c,:) = 0._r8
       mss_cnc_dst1(c,:)=0._r8
       mss_cnc_dst2(c,:)=0._r8
       mss_cnc_dst3(c,:)=0._r8
       mss_cnc_dst4(c,:)=0._r8
       
       if (snl(c) < 0) then
          snw_rds(c,snl(c)+1:0)        = snw_rds_min
          snw_rds(c,-nlevsno+1:snl(c)) = 0._r8
          snw_rds_top(c)               = snw_rds_min
          sno_liq_top(c) = h2osoi_liq(c,snl(c)+1) / (h2osoi_liq(c,snl(c)+1)+h2osoi_ice(c,snl(c)+1))
       elseif (h2osno(c) > 0._r8) then
          snw_rds(c,0)             = snw_rds_min
          snw_rds(c,-nlevsno+1:-1) = 0._r8
          snw_rds_top(c)           = spval
          sno_liq_top(c)           = spval
       else
          snw_rds(c,:)   = 0._r8
          snw_rds_top(c) = spval
          sno_liq_top(c) = spval
       endif
    enddo


    end associate 
   end subroutine mkarbinit

!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: snow_depth2lev
!
! !INTERFACE:
  subroutine snow_depth2lev(lbc, ubc)
!
! !DESCRIPTION:
! Create snow layers and interfaces given snow depth.
! Note that cps%zi(0) is set in routine iniTimeConst.
!
! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use clmtype
  use clm_varpar  , only : nlevsno
!
! !ARGUMENTS:
  implicit none
  integer, intent(in) :: lbc, ubc                    ! column bounds
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
! !LOCAL VARIABLES:
!
!
!
!
!
!
! !LOCAL VARIABLES:
!EOP
  integer :: c,l,j      !indices
!-----------------------------------------------------------------------
 

   associate(& 
   lakpoi                              =>   lun%lakpoi                                   , & ! Input:  [logical (:)]  true => landunit is a lake point                   
   clandunit                           =>   col%landunit                                 , & ! Input:  [integer (:)]  landunit index associated with each column         
   snow_depth                          =>    cps%snow_depth                              , & ! Input:  [real(r8) (:)]  snow height (m)                                   
   snl                                 =>    cps%snl                                     , & ! Output: [integer (:)]  number of snow layers                              
   zi                                  =>    cps%zi                                      , & ! Output: [real(r8) (:,:)]  interface depth (m) over snow only              
   dz                                  =>    cps%dz                                      , & ! Output: [real(r8) (:,:)]  layer thickness depth (m) over snow only        
   z                                   =>    cps%z                                         & ! Output: [real(r8) (:,:)]  layer depth  (m) over snow only                 
   )

  ! Initialize snow levels and interfaces (lake and non-lake points)

  do c = lbc, ubc
     dz(c,-nlevsno+1: 0) = 1.e36_r8
     z (c,-nlevsno+1: 0) = 1.e36_r8
     zi(c,-nlevsno  :-1) = 1.e36_r8
  end do

  ! Determine snow levels and interfaces for non-lake points

  do c = lbc,ubc
     l = clandunit(c)
     if (.not. lakpoi(l)) then
        if (snow_depth(c) < 0.01_r8) then
           snl(c) = 0
           dz(c,-nlevsno+1:0) = 0._r8
           z (c,-nlevsno+1:0) = 0._r8
           zi(c,-nlevsno+0:0) = 0._r8
        else
           if ((snow_depth(c) >= 0.01_r8) .and. (snow_depth(c) <= 0.03_r8)) then
              snl(c) = -1
              dz(c,0)  = snow_depth(c)
           else if ((snow_depth(c) > 0.03_r8) .and. (snow_depth(c) <= 0.04_r8)) then
              snl(c) = -2
              dz(c,-1) = snow_depth(c)/2._r8
              dz(c, 0) = dz(c,-1)
           else if ((snow_depth(c) > 0.04_r8) .and. (snow_depth(c) <= 0.07_r8)) then
              snl(c) = -2
              dz(c,-1) = 0.02_r8
              dz(c, 0) = snow_depth(c) - dz(c,-1)
           else if ((snow_depth(c) > 0.07_r8) .and. (snow_depth(c) <= 0.12_r8)) then
              snl(c) = -3
              dz(c,-2) = 0.02_r8
              dz(c,-1) = (snow_depth(c) - 0.02_r8)/2._r8
              dz(c, 0) = dz(c,-1)
           else if ((snow_depth(c) > 0.12_r8) .and. (snow_depth(c) <= 0.18_r8)) then
              snl(c) = -3
              dz(c,-2) = 0.02_r8
              dz(c,-1) = 0.05_r8
              dz(c, 0) = snow_depth(c) - dz(c,-2) - dz(c,-1)
           else if ((snow_depth(c) > 0.18_r8) .and. (snow_depth(c) <= 0.29_r8)) then
              snl(c) = -4
              dz(c,-3) = 0.02_r8
              dz(c,-2) = 0.05_r8
              dz(c,-1) = (snow_depth(c) - dz(c,-3) - dz(c,-2))/2._r8
              dz(c, 0) = dz(c,-1)
           else if ((snow_depth(c) > 0.29_r8) .and. (snow_depth(c) <= 0.41_r8)) then
              snl(c) = -4
              dz(c,-3) = 0.02_r8
              dz(c,-2) = 0.05_r8
              dz(c,-1) = 0.11_r8
              dz(c, 0) = snow_depth(c) - dz(c,-3) - dz(c,-2) - dz(c,-1)
           else if ((snow_depth(c) > 0.41_r8) .and. (snow_depth(c) <= 0.64_r8)) then
              snl(c) = -5
              dz(c,-4) = 0.02_r8
              dz(c,-3) = 0.05_r8
              dz(c,-2) = 0.11_r8
              dz(c,-1) = (snow_depth(c) - dz(c,-4) - dz(c,-3) - dz(c,-2))/2._r8
              dz(c, 0) = dz(c,-1)
           else if (snow_depth(c) > 0.64_r8) then
              snl(c) = -5
              dz(c,-4) = 0.02_r8
              dz(c,-3) = 0.05_r8
              dz(c,-2) = 0.11_r8
              dz(c,-1) = 0.23_r8
              dz(c, 0)=snow_depth(c)-dz(c,-4)-dz(c,-3)-dz(c,-2)-dz(c,-1)
           endif
        end if
     end if
  end do

  ! The following loop is currently not vectorized

  do c = lbc,ubc
     l = clandunit(c)
     if (.not. lakpoi(l)) then
        do j = 0, snl(c)+1, -1
           z(c,j)    = zi(c,j) - 0.5_r8*dz(c,j)
           zi(c,j-1) = zi(c,j) - dz(c,j)
        end do
     end if
  end do

  ! Determine snow levels and interfaces for lake points

  do c = lbc,ubc
     l = clandunit(c)
     if (lakpoi(l)) then
        snl(c) = 0
        dz(c,-nlevsno+1:0) = 0._r8
        z (c,-nlevsno+1:0) = 0._r8
        zi(c,-nlevsno+0:0) = 0._r8
     end if
  end do

    end associate 
   end subroutine snow_depth2lev


!-----------------------------------------------------------------------

end module mkarbinitMod
