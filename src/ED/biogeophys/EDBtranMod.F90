module EDBtranMod
   
   !-------------------------------------------------------------------------------------
   ! Description:
   ! 
   ! ------------------------------------------------------------------------------------
   
   use pftconMod         , only : pftcon
   use clm_varcon        , only : tfrz
   use EDTypesMod        , only : ed_site_type,       &
                                  ed_patch_type,      &
                                  ed_cohort_type,     &
                                  numpft_ed,          &
                                  cp_numlevgrnd
   use shr_kind_mod      , only : r8 => shr_kind_r8
   use FatesInterfaceMod , only : bc_in_type, &
                                  bc_out_type
   use clm_varctl        , only : iulog   !INTERF-TODO: THIS SHOULD BE MOVED

   !
   implicit none
   private
   
   public :: btran_ed
   public :: get_active_suction_layers
   
contains 
   
  ! ====================================================================================

  logical function check_layer_water(h2o_liq_vol, tempk)
    
    implicit none
    ! Arguments
    real(r8),intent(in) :: h2o_liq_vol
    real(r8),intent(in) :: tempk
    
    check_layer_water = .false.

    if ( h2o_liq_vol .gt. 0._r8 ) then
       if ( tempk .gt. tfrz-2._r8) then
          check_layer_water = .true.
       end if
    end if
    return
  end function check_layer_water

  ! =====================================================================================
  
  subroutine get_active_suction_layers(nsites, sites, bc_in, bc_out)
    
    ! Arguments
    
    integer,intent(in)                      :: nsites
    type(ed_site_type),intent(inout),target :: sites(nsites)
    type(bc_in_type),intent(in)             :: bc_in(nsites)
    type(bc_out_type),intent(inout)         :: bc_out(nsites)
    
    ! !LOCAL VARIABLES:
    integer  :: s                 ! site
    integer  :: j                 ! soil layer
    !------------------------------------------------------------------------------
    
      do s = 1,nsites
         if (bc_in(s)%filter_btran) then
            do j = 1,cp_numlevgrnd
               bc_out(s)%active_suction_gl(j) = check_layer_water( bc_in(s)%h2o_liqvol_gl(j),bc_in(s)%tempk_gl(j) )
            end do
         else
            bc_out(s)%active_suction_gl(:) = .false.
         end if
      end do

  end subroutine get_active_suction_layers
  
  ! =====================================================================================

  subroutine btran_ed( nsites, sites, bc_in, bc_out)
      
      ! ---------------------------------------------------------------------------------
      ! Calculate the transpiration wetness function (BTRAN) and the root uptake
      ! distribution (ROOTR).
      ! Boundary conditions in: bc_in(s)%eff_porosity_gl(j)    unfrozen porosity
      !                         bc_in(s)%watsat_gl(j)          porosity
      !                         bc_in(s)%active_uptake_gl(j)   frozen/not frozen
      !                         bc_in(s)%smp_gl(j)             suction
      ! Boundary conditions out: bc_out(s)%rootr_pagl          root uptake distribution
      !                          bc_out(s)%btran_pa            wetness factor
      ! ---------------------------------------------------------------------------------
      
      ! Arguments
      
      integer,intent(in)                      :: nsites
      type(ed_site_type),intent(inout),target :: sites(nsites)
      type(bc_in_type),intent(in)             :: bc_in(nsites)
      type(bc_out_type),intent(inout)         :: bc_out(nsites)
      
      !
      ! !LOCAL VARIABLES:
      type(ed_patch_type),pointer             :: cpatch ! Current Patch Pointer
      type(ed_cohort_type),pointer            :: ccohort ! Current cohort pointer
      integer  :: s                 ! site
      integer  :: j                 ! soil layer
      integer  :: ifp               ! patch vector index for the site
      integer  :: ft                ! plant functional type index
      real(r8) :: smp_node          ! matrix potential
      real(r8) :: rresis            ! suction limitation to transpiration independent
                                    ! of root density
      real(r8) :: pftgs(numpft_ed)  ! pft weighted stomatal conductance s/m
      real(r8) :: temprootr                   
      !------------------------------------------------------------------------------
      
      associate(                                 &
            smpsc     => pftcon%smpsc          , &  ! INTERF-TODO: THESE SHOULD BE FATES PARAMETERS
            smpso     => pftcon%smpso            &  ! INTERF-TODO: THESE SHOULD BE FATES PARAMETERS
            )
        
        do s = 1,nsites

           bc_out(s)%rootr_pagl(:,:) = 0._r8

           ifp = 0
           cpatch => sites(s)%oldest_patch
           do while (associated(cpatch))                 
              ifp=ifp+1
              
              ! THIS SHOULD REALLY BE A COHORT LOOP ONCE WE HAVE rootfr_ft FOR COHORTS (RGK)
              
              do ft = 1,numpft_ed
                 cpatch%btran_ft(ft) = 0.0_r8
                 do j = 1,cp_numlevgrnd
                    
                    ! Calculations are only relevant where liquid water exists
                    ! see clm_fates%wrap_btran for calculation with CLM/ALM
                    
                    if ( check_layer_water(bc_in(s)%h2o_liqvol_gl(j),bc_in(s)%tempk_gl(j)) )  then
                       
                       smp_node = max(smpsc(ft), bc_in(s)%smp_gl(j))
                       
                       rresis  = min( (bc_in(s)%eff_porosity_gl(j)/bc_in(s)%watsat_gl(j))*               &
                            (smp_node - smpsc(ft)) / (smpso(ft) - smpsc(ft)), 1._r8)
                       
                       cpatch%rootr_ft(ft,j) = cpatch%rootfr_ft(ft,j)*rresis
                       
                       ! root water uptake is not linearly proportional to root density,
                       ! to allow proper deep root funciton. Replace with equations from SPA/Newman. FIX(RF,032414)
                       ! cpatch%rootr_ft(ft,j) = cpatch%rootfr_ft(ft,j)**0.3*rresis_ft(ft,j)/ &
                       ! sum(cpatch%rootfr_ft(ft,1:nlevgrnd)**0.3)
                       cpatch%btran_ft(ft) = cpatch%btran_ft(ft) + cpatch%rootr_ft(ft,j)
                       
                    else
                       cpatch%rootr_ft(ft,j) = 0._r8
                    end if
                    
                 end do !j
                 
                 ! Normalize root resistances to get layer contribution to ET
                 do j = 1,cp_numlevgrnd    
                    if (cpatch%btran_ft(ft)  >  0.0_r8) then
                       cpatch%rootr_ft(ft,j) = cpatch%rootr_ft(ft,j)/cpatch%btran_ft(ft)
                    else
                       cpatch%rootr_ft(ft,j) = 0._r8
                    end if
                 end do
                 
              end do !PFT
              
              ! PFT-averaged point level root fraction for extraction purposese.
              ! This probably needs to be weighted by actual transpiration from each pft. FIX(RF,032414).
              pftgs(:) = 0._r8
              ccohort => cpatch%tallest
              do while(associated(ccohort))
                 pftgs(ccohort%pft) = pftgs(ccohort%pft) + ccohort%gscan * ccohort%n    
                 ccohort => ccohort%shorter
              enddo
              
              ! Process the boundary output, this is necessary for calculating the soil-moisture
              ! sink term across the different layers in driver/host.  Photosynthesis will
              ! pass the host a total transpiration for the patch.  This needs rootr to be
              ! distributed over the soil layers.
              
              do j = 1,cp_numlevgrnd
                 bc_out(s)%rootr_pagl(ifp,j) = 0._r8
                 do ft = 1,numpft_ed
                    if(sum(pftgs) > 0._r8)then !prevent problem with the first timestep - might fail
                       !bit-retart test as a result? FIX(RF,032414)  
                       bc_out(s)%rootr_pagl(ifp,j) = bc_out(s)%rootr_pagl(ifp,j) + &
                            cpatch%rootr_ft(ft,j) * pftgs(ft)/sum(pftgs)
                    else
                       bc_out(s)%rootr_pagl(ifp,j) = bc_out(s)%rootr_pagl(ifp,j) + &
                            cpatch%rootr_ft(ft,j) * 1./numpft_ed
                    end if
                 enddo
              enddo
              
              !weight patch level output BTRAN for the
              bc_out(s)%btran_pa(ifp) = 0.0_r8
              do ft = 1,numpft_ed
                 if(sum(pftgs) > 0._r8)then !prevent problem with the first timestep - might fail
                    !bit-retart test as a result? FIX(RF,032414)   
                    bc_out(s)%btran_pa(ifp)   = bc_out(s)%btran_pa(ifp) + cpatch%btran_ft(ft)  * pftgs(ft)/sum(pftgs)
                 else
                    bc_out(s)%btran_pa(ifp)   = bc_out(s)%btran_pa(ifp) + cpatch%btran_ft(ft) * 1./numpft_ed
                 end if
              enddo
              
              temprootr = sum(bc_out(s)%rootr_pagl(ifp,:))
              if(abs(1.0_r8-temprootr) > 1.0e-10_r8 .and. temprootr > 1.0e-10_r8)then
                 write(iulog,*) 'error with rootr in canopy fluxes',temprootr,sum(pftgs),sum(cpatch%rootr_ft(1:2,:),dim=2)
                 do j = 1,cp_numlevgrnd
                    bc_out(s)%rootr_pagl(ifp,j) = bc_out(s)%rootr_pagl(ifp,j)/temprootr
                 enddo
              end if
              
              cpatch => cpatch%younger
           end do
        

        end do
        
      end associate
      
    end subroutine btran_ed

   ! =========================================================================================

         !---------------------------------------------------------------------------------------
         ! SPA based recalculation of BTRAN and water uptake. 
         !---------------------------------------------------------------------------------------

!         if (SPA_soil) then   ! normal case don't run this.
!            rootr(p,:) = 0._r8
!            do FT = 1,numpft_ed 

!               ! Soil Physics  
!               do j = 1,nlevgrnd
!                  ! CLM water retention curve. Clapp and Hornberger equation.    
!                  s1 = max(h2osoi_vol(c,j)/watsat(c,j), 0.01_r8)
!                  s1 = min(1.0_r8,s1)  
!                  smp_node = -sucsat(c,j)*s1**(-bsw(c,j))
!                  swp_mpa(j)  = smp_node *10.0_r8/1000000.0_r8  !convert from mm to Mpa 

!                  ! CLM hydraulic conductivity curve. 
!                  ! As opposed to the Richard's equation solution in SoilHydrology.Mod 
!                  ! the conductivity here is defined in the middle of the layer in question, not at the edge... 
!                  xksat   = 0.0070556_r8 * (10._r8**(-0.884_r8+0.0153_r8*sand(p)) ) 
!                  hk(j)    =   xksat*s1**(2._r8*bsw(c,j)+2._r8)  !removed the ice from here to avoid 1st ts crashing        
!               enddo

!               ! Root resistance
!               rootxsecarea=3.14159*rootrad**2
!               do j = 1,nlevgrnd      
!                  rootmass(j) =  EDecophyscon%soilbeta(FT) * cpatch%rootfr_ft(FT,j)
!                  rootlength(j) = rootmass(j)/(rootdens*rootxsecarea)   !m m-3 soil
!                  Lsoil(j)     = hk(j)/1000/head !converts from mms-1 to ms-1 and then to m2 s-1 MPa-1    
!                  if(Lsoil(j) < 1e-35_r8.or.cpatch%rootfr_ft(ft,j) <= 0.0_r8)then   !prevent floating point error
!                     soilr_z(j) = 1e35_r8
!                     soilr2(j)  = 1e35_r8
!                  else 
!                     ! Soil-to-root water uptake from Newman (1969). 
!                     rs = sqrt (1._r8  / (rootlength(j) * pi)) 
!                     soilr1(j) = log(rs/rootrad) / (2.0_r8 * pi * rootlength(j) * Lsoil(j) * dz(c,j)) 
!                     ! convert from MPa s m2 m-3 to MPa s m2 mmol-1     
!                     soilr1(j) = soilr1(j) * 1E-6_r8 * 18_r8 * 0.001_r8  
!                     ! second component of below ground resistance is related to root hydraulics
!                     soilr2(j) = EDecophyscon%rootresist(FT)/(rootmass(j)*dz(c,j))
!                     soilr_z(j) = soilr1(j)+soilr2(j)
!                  end if
!               enddo

               ! Aggregate soil layers
!               totestevap=0._r8
!               weighted_SWP=0._r8
!               estevap=0._r8
!               fraction_uptake=0._r8
!               canopy_soil_resistance=0._r8  !Reset Counters
!               totmaxevap = 0._r8

               ! Estimated max transpiration from LWP gradient / soil resistance
!               do j = 1,nlevgrnd    
!                  estevap(j) = (swp_mpa(j) - minlwp)/(soilr_z(j))
!                  estevap(j) = max(0._r8,estevap(j))         ! no negative uptake 
!                  maxevap(j) = (0.0_r8 - minlwp)/(soilr2(j)) 
!               enddo
!               totestevap = sum(estevap)
!               totmaxevap = sum(maxevap)  

               ! Weighted soil water potential
!               do j = 1,nlevgrnd 
!                  if(totestevap > 0._r8)then
!                     fraction_uptake(j) = estevap(j)/totestevap   !Fraction of total ET taken from this soil layer 
!                  else
!                     estevap(j) = 0._r8
!                     fraction_uptake(j)=1._r8/nlevgrnd
!                  end if
!                  weighted_SWP = weighted_SWP + swp_mpa(j) * estevap(j)    
!               enddo

!               if(totestevap > 0._r8)then
!                  weighted_swp = weighted_swp/totestevap 
!                  ! weight SWP for the total evaporation 
!               else   
!                  write(iulog,*) 'empty soil', totestevap
!                  ! error check
!                  weighted_swp = minlwp
!               end if

               ! Weighted soil-root resistance. Aggregate the conductances (1/soilR) for each soil layer
!               do iv = 1,nv !leaf layers
!                  fleaf = 1.0_r8/nv
!                  do j = 1,nlevgrnd !root layers
!                     ! Soil resistance for each canopy layer is related to leaf area
!                     ! The conductance of the root system to the 
!                     ! whole canopy is reduced by the fraction of leaves in this layer...
!                     canopy_soil_resistance(iv) = canopy_soil_resistance(iv)+fleaf * 1.0_r8/(soilr_z(j))         
!                  enddo
!                  ! Turn aggregated conductance back into resistance. mmol MPa-1 s-1 m-2  to  MPa s m2 mmol-1
!                  canopy_soil_resistance(iv) = 1./canopy_soil_resistance(iv)
!               enddo
!
!               cpatch%btran_ft(FT) =  totestevap/totmaxevap         
!               do j = 1,nlevgrnd       
!                  if(sum(pftgs) > 0._r8)then !prevent problem with the first timestep - might fail
!                     !bit-retart test as a result? FIX(RF,032414)   
!                     rootr(p,j) = rootr(p,j) + fraction_uptake(j) * pftgs(ft)/sum(pftgs)
!                  else
!                     rootr(p,j) = rootr(p,j) + fraction_uptake(j) * 1./numpft_ed
!                  end if
!               enddo
!            enddo !pft loop
!         end if !
         !---------------------------------------------------------------------------------------
         ! end of SPA based recalculation of BTRAN and water uptake. 
         !---------------------------------------------------------------------------------------

  

end module EDBtranMod
