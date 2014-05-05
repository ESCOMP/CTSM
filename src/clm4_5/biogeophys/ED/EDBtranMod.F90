 module EDBtranMod

  !------------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! This routine accumulates NPP, GPP and respiration of each cohort over the course of each 24 hour period. 
  ! The fluxes are stored per cohort, and the npp_clm (etc) fluxes are calcualted in EDPhotosynthesis
  ! This routine cannot be in EDPhotosynthesis because EDPhotosynthesis is a loop and therefore would
  ! erroneously add these things up multiple times. 
  ! Rosie Fisher. March 2014. 
  
  
  ! !USES:
  use EDtypesMod, only : patch, cohort, gridCellEdState,numpft_ed
  use clmtype,    only : gps, ces, cps, cws, pft, pps, pftcon

  implicit none
  save
  private
  
  public :: BTRAN_ED
  
  type(cohort), pointer  :: currentCohort ! current cohort
  type(patch) , pointer  :: currentPatch ! current patch
   
  contains 
  !------------------------------------------------------------------------------
  subroutine BTRAN_ED(p)
      
   use shr_kind_mod       ,  only : r8 => shr_kind_r8
   use clm_varpar         ,  only : nlevgrnd
   use clm_varctl         ,  only : iulog
   use clm_varcon         ,  only : tfrz, denice, denh2o
   use shr_const_mod      ,  only : shr_const_pi  
   use EDclmType          ,  only : EDpft, EDpftcon
   
   implicit none

   integer, intent(in) :: p  !patch/'p'

   integer :: iv !leaf layer
   integer :: g  !gridcell
   integer :: c  !column
   integer :: j  !soil layer
   integer :: ft ! plant functional type index
 
    ! Inputs to model from CLM. To be read in through an input file for the purposes of this program.      
    integer, parameter    :: nv = 5        ! Number of canopy layers
    real(r8) :: xksat                      ! maximum hydraulic conductivity of soil [mm/s]
    real(r8) :: s1                         ! HC intermediate
    real(r8) :: swp_mpa(nlevgrnd)          ! matrix potential - MPa
    real(r8) :: hk(nlevgrnd)               ! hydraulic conductivity [mm h2o/s] 
    real(r8) :: rootxsecarea               ! root X-sectional area (m2)
    real(r8) :: rootmass(nlevgrnd)         ! root mass in each layer (g) 
    real(r8) :: rootlength(nlevgrnd)       ! root length in each layer (m) 
    real(r8) :: soilr1(nlevgrnd)           ! soil-to-root resistance in each layer (MPa s m2 mmol-1)
    real(r8) :: soilr2(nlevgrnd)           ! internal root resistance in each layer (MPa s m2 mmol-1)
    real(r8) :: rs                         ! intermediate variable
    real(r8) :: soilr_z(nlevgrnd)          ! soil-to-xylem resistance in each layer (MPa s m2 mmol-1) 
    real(r8) :: lsoil(nlevgrnd)            ! hydraulic conductivity in each soil layer  

    real(r8) :: estevap(nlevgrnd)          ! potential suction from each soil layer (mmol m-2 s-1)
    real(r8) :: totestevap                 ! potential suction from each soil layer (mmol m-2 s-1)
    real(r8) :: fraction_uptake(nlevgrnd)  ! Uptake of water from each soil layer (-) 
    real(r8) :: maxevap(nlevgrnd)          ! potential suction from each soil layer (mmol m-2 s-1)
    real(r8) :: totmaxevap                 ! potential suction from each soil layer (mmol m-2 s-1)
    real(r8) :: fleaf                      ! fraction of leaves in each canopy layer

    ! Model parameters
    real(r8) :: head         = 0.009807_r8    ! head of pressure  (MPa/m)   
    real(r8) :: rootdens     = 0.5e6_r8     ! root density, g biomass m-3 root
    real(r8) :: pi           = shr_const_pi
    real(r8) :: vol_ice               ! partial volume of ice lens in layer
    real(r8) :: eff_porosity          ! effective porosity in layer
    real(r8) :: vol_liq               ! partial volume of liquid water in layer
    real(r8) :: s_node                ! vol_liq/eff_porosity
    real(r8) :: smp_node              ! matrix potential

    ! To be read in from pft file ultimately.       
    real(r8) :: minlwp  =  -2.5_r8             ! minimum leaf water potential in MPa
    real(r8) :: rootrad  =  0.001_r8           ! root radius in metres

    ! Outputs to CLM_SPA
    real(r8) :: weighted_SWP                   ! weighted apparent soil water potential: MPa.
    real(r8) :: canopy_soil_resistance(nv)     ! Resistance experienced by each canopy layer: MPa s m2 mmol-1

    ! SPA Pointers from CLM type. 
    logical, parameter :: SPA_soil=.false.     ! Is the BTRAN model SPA or CLM? FIX(SPM,032414) ed - make this a namelist var

    real(r8) :: rresis_ft(numpft_ed,nlevgrnd)  ! resistance to water uptake per pft and soil layer.
    real(r8) :: pftgs(numpft_ed)               !pft weighted stomatal conductance s/m
    real(r8) :: temprootr                   
  
     ! SPA Pointers from CLM type. 

    !------------------------------------------------------------------------------

   associate(& 
   max_dayl    => gps%max_dayl    , & ! Input:  [real(r8) (:)]  maximum daylength for this grid cell (s)
   t_soisno    => ces%t_soisno    , & ! Input:  [real(r8) (:,:)]  soil temperature (Kelvin)
   watsat      => cps%watsat      , & ! Input:  [real(r8) (:,:)]  volumetric soil water at saturation (porosity)
   watdry      => cps%watdry      , & ! Input:  [real(r8) (:,:)]  btran parameter for btran=0
   watopt      => cps%watopt      , & ! Input:  [real(r8) (:,:)]  btran parameter for btran = 1
   h2osoi_ice  => cws%h2osoi_ice  , & ! Input:  [real(r8) (:,:)]  ice lens (kg/m2)
   h2osoi_vol  => cws%h2osoi_vol  , & ! Input:  [real(r8) (:,:)]  
                                      ! volumetric soil water (0<=h2osoi_vol<=watsat) [m3/m3] by F. Li and S. Levis
   dz          => cps%dz          , & ! Input:  [real(r8) (:,:)]  layer depth (m)
   h2osoi_liq  => cws%h2osoi_liq  , & ! Input:  [real(r8) (:,:)]  liquid water (kg/m2) 
   sucsat      => cps%sucsat      , & ! Input:  [real(r8) (:,:)]  minimum soil suction (mm) 
   t_grnd      => ces%t_grnd      , & ! Input:  [real(r8) (:)]  ground surface temperature [K]
   soilbeta    => cws%soilbeta    , & ! Input:  [real(r8) (:)]  soil wetness relative to field capacity 
   ivt         => pft%itype       , & ! Input:  [integer (:)]  pft vegetation type
   pcolumn     => pft%column      , & ! Input:  [integer (:)]  pft's column index  
   plandunit   => pft%landunit    , & ! Input:  [integer (:)]  pft's landunit index
   pgridcell   => pft%gridcell    , & ! Input:  [integer (:)]  pft's gridcell index 
   rootfr      => pps%rootfr      , & ! Input:  [real(r8) (:,:)]  fraction of roots in each soil layer 
   elai        => pps%elai        , & ! Input:  [real(r8) (:)]  one-sided leaf area index with burying by snow 
   esai        => pps%esai        , & ! Input:  [real(r8) (:)]  one-sided stem area index with burying by snow
   fdry        => pps%fdry        , & ! Input:  [real(r8) (:)]  fraction of foliage that is green and dry [-]
   bsw         => cps%bsw         , & ! Input:  [real(r8) (:,:)]  Clapp and Hornberger "b" 
   laisun      => pps%laisun      , & ! Input:  [real(r8) (:)]  sunlit leaf area 
   laisha      => pps%laisha      , & ! Input:  [real(r8) (:)]  shaded leaf area  
   smpso       => pftcon%smpso    , & ! Input:  [real(r8) (:)]  soil water potential at full stomatal opening (mm)
   smpsc       => pftcon%smpsc    , & ! Input:  [real(r8) (:)]  soil water potential at full stomatal closure (mm)
   rb1         => pps%rb1         , & ! Output: [real(r8) (:)]  boundary layer resistance (s/m)
   btran       => pps%btran       , & ! Output: [real(r8) (:)]  transpiration wetness factor (0 to 1)
   btran2      => pps%btran2      , & ! Output: [real(r8) (:)] F. Li and S. Levis 
   rresis      => pps%rresis      , & ! Output: [real(r8) (:,:)]  root resistance by layer (0-1)  (nlevgrnd) 
   canopy_cond => pps%canopy_cond , & ! Output: [real(r8) (:)] tracer conductance for canopy [m/s]
   sand        => pps%sandfrac    , & ! Input:  [real(r8) (:)]   % sand of soil 
   rootr       => pps%rootr       , & ! Output: [real(r8) (:,:)] Fraction of water uptake in each layer
   ED_patch    => EDpft%ED_patch  , &
   rscanopy    => pps%rscanopy      & ! Output: [real(r8) (:,:)]    !  canopy resistance s/m  
   )
   
          if(ED_patch(p)==1)then
             g = pgridcell(p)
             currentPatch => gridCellEdState(g)%spnt%oldest_patch   
             do while(p /= currentPatch%clm_pno)
                currentPatch => currentPatch%younger
             enddo

             c = pcolumn(p)
             do FT = 1,numpft_ed
                currentPatch%btran_ft(FT) = 0.0_r8
                do j = 1,nlevgrnd
                   !Root resistance factors
                   vol_ice = min(watsat(c,j), h2osoi_ice(c,j)/(dz(c,j)*denice))
                   eff_porosity = watsat(c,j)-vol_ice
                   vol_liq = min(eff_porosity, h2osoi_liq(c,j)/(dz(c,j)*denh2o))
                   if (vol_liq  <=  0._r8 .or. t_soisno(c,j)  <=  tfrz-2._r8) then
                      currentPatch%rootr_ft(FT,j) = 0._r8
                   else
                      s_node = max(vol_liq/eff_porosity,0.01_r8)
                      smp_node = max(smpsc(FT), -sucsat(c,j)*s_node**(-bsw(c,j)))
                      !FIX(RF,032414) for junipers
                      rresis_ft(FT,j) = min( (eff_porosity/watsat(c,j))* &
                           (smp_node - smpsc(FT)) / (smpso(FT) - smpsc(FT)), 1._r8)

                      currentPatch%rootr_ft(FT,j) = currentPatch%rootfr_ft(FT,j)*rresis_FT(FT,j)
                      ! root water uptake is not linearly proportional to root density,
                      ! to allow proper deep root funciton. Replace with equations from SPA/Newman. FIX(RF,032414)
                      ! currentPatch%rootr_ft(FT,j) = currentPatch%rootfr_ft(FT,j)**0.3*rresis_FT(FT,j)/ &
                      ! sum(currentPatch%rootfr_ft(FT,1:nlevgrnd)**0.3)
                      currentPatch%btran_ft(FT) = currentPatch%btran_ft(FT) + currentPatch%rootr_ft(FT,j)
                   end if
                end do !j
                btran(p) = currentPatch%btran_ft(1) !FIX(RF,032414) for TRF where is this used?
                ! Normalize root resistances to get layer contribution to ET 
                do j = 1,nlevgrnd    
                   if (currentPatch%btran_ft(FT)  >  0.0_r8) then
                      currentPatch%rootr_ft(FT,j) = currentPatch%rootr_ft(FT,j)/currentPatch%btran_ft(FT)
                   else
                      currentPatch%rootr_ft(FT,j) = 0._r8
                   end if
                end do

             end do !PFT

             ! PFT-averaged point level root fraction for extraction purposese.
             ! This probably needs to be weighted by actual transpiration from each pft. FIX(RF,032414).
             pftgs(:) = 0._r8
             currentCohort => currentPatch%tallest
             do while(associated(currentCohort))
                pftgs(currentCohort%pft) = pftgs(currentCohort%pft) + currentCohort%gscan * currentCohort%n    
                currentCohort => currentCohort%shorter
             enddo

             do j = 1,nlevgrnd       
                rootr(p,j) = 0._r8
                btran(p) = 0.0_r8
                do FT = 1,numpft_ed
                   if(sum(pftgs) > 0._r8)then !prevent problem with the first timestep - might fail
                      !bit-retart test as a result? FIX(RF,032414)  
                      rootr(p,j) = rootr(p,j) + currentPatch%rootr_ft(FT,j) * pftgs(ft)/sum(pftgs)
                   else
                      rootr(p,j) = rootr(p,j) + currentPatch%rootr_ft(FT,j) * 1./numpft_ed
                   end if
                enddo
             enddo


           !---------------------------------------------------------------------------------------
             ! SPA based recalculation of BTRAN and water uptake. 
             !---------------------------------------------------------------------------------------

             if (SPA_soil) then   ! normal case don't run this.
                rootr(p,:) = 0._r8
                do FT = 1,numpft_ed 
  
                   ! Soil Physics  
                   do j = 1,nlevgrnd
                      ! CLM water retention curve. Clapp and Hornberger equation.    
                      s1 = max(h2osoi_vol(c,j)/watsat(c,j), 0.01_r8)
                      s1 = min(1.0_r8,s1)  
                      smp_node = -sucsat(c,j)*s1**(-bsw(c,j))
                      swp_mpa(j)  = smp_node *10_r8/1000000_r8  !convert from mm to Mpa 

                      ! CLM hydraulic conductivity curve. 
                      ! As opposed to the Richard's equation solution in SoilHydrology.Mod 
                      ! the conductivity here is defined in the middle of the layer in question, not at the edge... 
                      xksat   = 0.0070556_r8 * (10._r8**(-0.884_r8+0.0153_r8*sand(p)) ) 
                      hk(j)    =   xksat*s1**(2._r8*bsw(c,j)+2._r8)  !removed the ice from here to avoid 1st ts crashing        
                   enddo

                   ! Root resistance
                   rootxsecarea=3.14159*rootrad**2
                   do j = 1,nlevgrnd      
                      rootmass(j) =  EDpftcon%soilbeta(FT) * currentPatch%rootfr_ft(FT,j)
                      rootlength(j) = rootmass(j)/(rootdens*rootxsecarea)   !m m-3 soil
                      Lsoil(j)     = hk(j)/1000/head !converts from mms-1 to ms-1 and then to m2 s-1 MPa-1    
                      if(Lsoil(j) < 1e-35_r8.or.currentPatch%rootfr_ft(ft,j) <= 0.0_r8)then   !prevent floating point error
                         soilr_z(j) = 1e35_r8
                         soilr2(j)  = 1e35_r8
                      else 
                         ! Soil-to-root water uptake from Newman (1969). 
                         rs = sqrt (1._r8  / (rootlength(j) * pi)) 
                         soilr1(j) = log(rs/rootrad) / (2.0_r8 * pi * rootlength(j) * Lsoil(j) * dz(c,j)) 
                         ! convert from MPa s m2 m-3 to MPa s m2 mmol-1     
                         soilr1(j) = soilr1(j) * 1E-6_r8 * 18_r8 * 0.001_r8  
                         ! second component of below ground resistance is related to root hydraulics
                         soilr2(j) = EDpftcon%rootresist(FT)/(rootmass(j)*dz(c,j))
                         soilr_z(j) = soilr1(j)+soilr2(j)
                      end if
                   enddo

                   ! Aggregate soil layers
                   totestevap=0._r8
                   weighted_SWP=0._r8
                   estevap=0._r8
                   fraction_uptake=0._r8
                   canopy_soil_resistance=0._r8  !Reset Counters
                   totmaxevap = 0._r8

                   ! Estimated max transpiration from LWP gradient / soil resistance
                   do j = 1,nlevgrnd    
                      estevap(j) = (swp_mpa(j) - minlwp)/(soilr_z(j))
                      estevap(j) = max(0._r8,estevap(j))         ! no negative uptake 
                      maxevap(j) = (0.0_r8 - minlwp)/(soilr2(j)) 
                   enddo
                   totestevap = sum(estevap)
                   totmaxevap = sum(maxevap)  
  
                   ! Weighted soil water potential
                   do j = 1,nlevgrnd 
                      if(totestevap > 0._r8)then
                         fraction_uptake(j) = estevap(j)/totestevap   !Fraction of total ET taken from this soil layer 
                      else
                         estevap(j) = 0._r8
                         fraction_uptake(j)=1._r8/nlevgrnd
                      end if
                      weighted_SWP = weighted_SWP + swp_mpa(j) * estevap(j)    
                   enddo


                   if(totestevap > 0._r8)then
                      weighted_swp = weighted_swp/totestevap 
                      ! weight SWP for the total evaporation 
                   else   
                      write(iulog,*) 'empty soil', totestevap
                      ! error check
                      weighted_swp = minlwp
                   end if

                   ! Weighted soil-root resistance. Aggregate the conductances (1/soilR) for each soil layer
                   do iv = 1,nv !leaf layers
                      fleaf = 1.0_r8/nv
                      do j = 1,nlevgrnd !root layers
                         ! Soil resistance for each canopy layer is related to leaf area
                         ! The conductance of the root system to the 
                         ! whole canopy is reduced by the fraction of leaves in this layer...
                         canopy_soil_resistance(iv) = canopy_soil_resistance(iv)+fleaf * 1.0_r8/(soilr_z(j))         
                      enddo
                      ! Turn aggregated conductance back into resistance. mmol MPa-1 s-1 m-2  to  MPa s m2 mmol-1
                      canopy_soil_resistance(iv) = 1./canopy_soil_resistance(iv)
                   enddo

                   currentPatch%btran_ft(FT) =  totestevap/totmaxevap         
                   do j = 1,nlevgrnd       
                     if(sum(pftgs) > 0._r8)then !prevent problem with the first timestep - might fail
                         !bit-retart test as a result? FIX(RF,032414)   
                         rootr(p,j) = rootr(p,j) + fraction_uptake(j) * pftgs(ft)/sum(pftgs)
                      else
                         rootr(p,j) = rootr(p,j) + fraction_uptake(j) * 1./numpft_ed
                      end if
                   enddo

                enddo !pft loop

             end if !
                !---------------------------------------------------------------------------------------
             ! end of SPA based recalculation of BTRAN and water uptake. 
             !---------------------------------------------------------------------------------------

             !weight patch level output BTRAN for the
             btran(p) = 0.0_r8
             do FT = 1,numpft_ed
                if(sum(pftgs) > 0._r8)then !prevent problem with the first timestep - might fail
                   !bit-retart test as a result? FIX(RF,032414)   
                   btran(p)   = btran(p) + currentPatch%btran_ft(FT)  * pftgs(ft)/sum(pftgs)
                else
                   btran(p)   = btran(p) + currentPatch%btran_ft(FT) * 1./numpft_ed
                end if
             enddo

             temprootr = sum(rootr(p,:))
             if(temprootr /= 1.0_r8)then
                !write(iulog,*) 'error with rootr in canopy fluxes',sum(rootr(p,:))
                if(temprootr > 0._r8)then
                   do j = 1,nlevgrnd
                      rootr(p,j) = rootr(p,j) / temprootr
                   enddo
                end if
             end if

          else ! edpatch
             currentPatch%btran_ft(1:numpft_ed) = 1._r8
          end if ! edpatch
     
     end associate
     
end subroutine BTRAN_ED

end module EDBtranMod
