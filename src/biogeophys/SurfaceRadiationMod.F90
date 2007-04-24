#include <misc.h>
#include <preproc.h>

module SurfaceRadiationMod

!------------------------------------------------------------------------------
!BOP
!
! !MODULE: SurfaceRadiationMod
!
! !DESCRIPTION:
! Calculate solar fluxes absorbed by vegetation and ground surface
!
! !USES:
   use shr_kind_mod, only: r8 => shr_kind_r8
!
! !PUBLIC TYPES:
   implicit none
   save
!
! !PUBLIC MEMBER FUNCTIONS:
   public :: SurfaceRadiation ! Solar fluxes absorbed by veg and ground surface
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
! 11/26/03, Peter Thornton: Added new routine for improved treatment of
!    sunlit/shaded canopy radiation.
! 4/26/05, Peter Thornton: Adopted the sun/shade algorithm as the default,
!    removed the old SurfaceRadiation(), and renamed SurfaceRadiationSunShade()
!    as SurfaceRadiation().
!
!EOP
!------------------------------------------------------------------------------

contains

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: SurfaceRadiation
!
! !INTERFACE:
   subroutine SurfaceRadiation(lbp, ubp)
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
     use clmtype
     use clm_atmlnd  , only : clm_a2l
     use clm_varpar  , only : numrad
     use clm_varcon  , only : spval
     use clm_time_manager, only : get_curr_date, get_step_size
!
! !ARGUMENTS:
     implicit none
     integer, intent(in) :: lbp, ubp      ! pft upper and lower bounds
!
! !CALLED FROM:
! subroutine Biogeophysics1 in module Biogeophysics1Mod
! subroutine BiogeophysicsLake in module BiogeophysicsLakeMod
!
! !REVISION HISTORY:
! Author: Gordon Bonan
! 2/18/02, Peter Thornton: Migrated to new data structures. Added a pft loop.
! 6/05/03, Peter Thornton: Modified sunlit/shaded canopy treatment. Original code
! had all radiation being absorbed in the sunlit canopy, and now the sunlit and shaded
! canopies are each given the appropriate fluxes.  There was also an inconsistency in
! the original code, where parsun was not being scaled by leaf area, and so represented
! the entire canopy flux.  This goes into Stomata (in CanopyFluxes) where it is assumed
! to be a flux per unit leaf area. In addition, the fpsn flux coming out of Stomata was
! being scaled back up to the canopy by multiplying by lai, but the input radiation flux was
! for the entire canopy to begin with.  Corrected this inconsistency in this version, so that
! the parsun and parsha fluxes going into canopy fluxes are per unit lai in the sunlit and
! shaded canopies.
! 6/9/03, Peter Thornton: Moved coszen from g%gps to c%cps to avoid problem
! with OpenMP threading over columns, where different columns hit the radiation
! time step at different times during execution.
! 6/10/03, Peter Thornton: Added constraint on negative tot_aid, instead of
! exiting with error. Appears to be happening only at roundoff level.
! 6/11/03, Peter Thornton: Moved calculation of ext inside if (coszen),
! and added check on laisun = 0 and laisha = 0 in calculation of sun_aperlai
! and sha_aperlai.
! 11/26/03, Peter Thornton: During migration to new vector code, created 
!   this as a new routine to handle sunlit/shaded canopy calculations.
!
! !LOCAL VARIABLES:
!
! local pointers to original implicit in arguments
!
     integer , pointer :: ivt(:)           ! pft vegetation type
     integer , pointer :: pcolumn(:)       ! pft's column index
     integer , pointer :: pgridcell(:)     ! pft's gridcell index
     real(r8), pointer :: pwtgcell(:)      ! pft's weight relative to corresponding gridcell
     real(r8), pointer :: elai(:)          ! one-sided leaf area index with burying by snow
     real(r8), pointer :: esai(:)          ! one-sided stem area index with burying by snow
     real(r8), pointer :: londeg(:)        ! longitude (degrees)
     real(r8), pointer :: latdeg(:)        ! latitude (degrees)
     real(r8), pointer :: slasun(:)        ! specific leaf area for sunlit canopy, projected area basis (m^2/gC)
     real(r8), pointer :: slasha(:)        ! specific leaf area for shaded canopy, projected area basis (m^2/gC)
     real(r8), pointer :: gdir(:)	   ! leaf projection in solar direction (0 to 1)
     real(r8), pointer :: omega(:,:)       ! fraction of intercepted radiation that is scattered (0 to 1)
     real(r8), pointer :: coszen(:)	   ! cosine of solar zenith angle
     real(r8), pointer :: forc_solad(:,:)  ! direct beam radiation (W/m**2)
     real(r8), pointer :: forc_solai(:,:)  ! diffuse radiation (W/m**2)
     real(r8), pointer :: fabd(:,:)        ! flux absorbed by veg per unit direct flux
     real(r8), pointer :: fabi(:,:)        ! flux absorbed by veg per unit diffuse flux
     real(r8), pointer :: ftdd(:,:)        ! down direct flux below veg per unit dir flx
     real(r8), pointer :: ftid(:,:)        ! down diffuse flux below veg per unit dir flx
     real(r8), pointer :: ftii(:,:)        ! down diffuse flux below veg per unit dif flx
     real(r8), pointer :: albgrd(:,:)      ! ground albedo (direct)
     real(r8), pointer :: albgri(:,:)      ! ground albedo (diffuse)
     real(r8), pointer :: albd(:,:)        ! surface albedo (direct)
     real(r8), pointer :: albi(:,:)        ! surface albedo (diffuse)
     real(r8), pointer :: slatop(:)        ! specific leaf area at top of canopy, projected area basis [m^2/gC]
     real(r8), pointer :: dsladlai(:)      ! dSLA/dLAI, projected area basis [m^2/gC]
!
! local pointers to original implicit out arguments
!
     real(r8), pointer :: fsun(:)          ! sunlit fraction of canopy
     real(r8), pointer :: laisun(:)        ! sunlit leaf area
     real(r8), pointer :: laisha(:)        ! shaded leaf area
     real(r8), pointer :: sabg(:)          ! solar radiation absorbed by ground (W/m**2)
     real(r8), pointer :: sabv(:)          ! solar radiation absorbed by vegetation (W/m**2)
     real(r8), pointer :: fsa(:)           ! solar radiation absorbed (total) (W/m**2)
     real(r8), pointer :: parsun(:)        ! average absorbed PAR for sunlit leaves (W/m**2)
     real(r8), pointer :: parsha(:)        ! average absorbed PAR for shaded leaves (W/m**2)
     real(r8), pointer :: fsr(:)           ! solar radiation reflected (W/m**2)
     real(r8), pointer :: fsds_vis_d(:)    ! incident direct beam vis solar radiation (W/m**2)
     real(r8), pointer :: fsds_nir_d(:)    ! incident direct beam nir solar radiation (W/m**2)
     real(r8), pointer :: fsds_vis_i(:)    ! incident diffuse vis solar radiation (W/m**2)
     real(r8), pointer :: fsds_nir_i(:)    ! incident diffuse nir solar radiation (W/m**2)
     real(r8), pointer :: fsr_vis_d(:)     ! reflected direct beam vis solar radiation (W/m**2)
     real(r8), pointer :: fsr_nir_d(:)     ! reflected direct beam nir solar radiation (W/m**2)
     real(r8), pointer :: fsr_vis_i(:)     ! reflected diffuse vis solar radiation (W/m**2)
     real(r8), pointer :: fsr_nir_i(:)     ! reflected diffuse nir solar radiation (W/m**2)
     real(r8), pointer :: fsds_vis_d_ln(:) ! incident direct beam vis solar rad at local noon (W/m**2)
     real(r8), pointer :: fsds_nir_d_ln(:) ! incident direct beam nir solar rad at local noon (W/m**2)
     real(r8), pointer :: fsr_vis_d_ln(:)  ! reflected direct beam vis solar rad at local noon (W/m**2)
     real(r8), pointer :: fsr_nir_d_ln(:)  ! reflected direct beam nir solar rad at local noon (W/m**2)
     real(r8), pointer :: eff_kid(:,:)     ! effective extinction coefficient for indirect from direct
     real(r8), pointer :: eff_kii(:,:)     ! effective extinction coefficient for indirect from indirect
     real(r8), pointer :: sun_faid(:,:)    ! fraction sun canopy absorbed indirect from direct
     real(r8), pointer :: sun_faii(:,:)    ! fraction sun canopy absorbed indirect from indirect
     real(r8), pointer :: sha_faid(:,:)    ! fraction shade canopy absorbed indirect from direct
     real(r8), pointer :: sha_faii(:,:)    ! fraction shade canopy absorbed indirect from indirect
     real(r8), pointer :: sun_add(:,:)     ! sun canopy absorbed direct from direct (W/m**2)
     real(r8), pointer :: tot_aid(:,:)     ! total canopy absorbed indirect from direct (W/m**2)
     real(r8), pointer :: sun_aid(:,:)     ! sun canopy absorbed indirect from direct (W/m**2)
     real(r8), pointer :: sun_aii(:,:)     ! sun canopy absorbed indirect from indirect (W/m**2)
     real(r8), pointer :: sha_aid(:,:)     ! shade canopy absorbed indirect from direct (W/m**2)
     real(r8), pointer :: sha_aii(:,:)     ! shade canopy absorbed indirect from indirect (W/m**2)
     real(r8), pointer :: sun_atot(:,:)    ! sun canopy total absorbed (W/m**2)
     real(r8), pointer :: sha_atot(:,:)    ! shade canopy total absorbed (W/m**2)
     real(r8), pointer :: sun_alf(:,:)     ! sun canopy total absorbed by leaves (W/m**2)
     real(r8), pointer :: sha_alf(:,:)     ! shade canopy total absored by leaves (W/m**2)
     real(r8), pointer :: sun_aperlai(:,:) ! sun canopy total absorbed per unit LAI (W/m**2)
     real(r8), pointer :: sha_aperlai(:,:) ! shade canopy total absorbed per unit LAI (W/m**2)
!
!EOP
!
! !OTHER LOCAL VARIABLES:
!
     integer , parameter :: nband = numrad ! number of solar radiation waveband classes
     real(r8), parameter :: mpe = 1.e-06_r8   ! prevents overflow for division by zero
     integer  :: p                   ! pft index
     integer  :: c                   ! column index
     integer  :: g                   ! grid cell index
     integer  :: ib                  ! waveband number (1=vis, 2=nir)
     real(r8) :: abs                 ! absorbed solar radiation (W/m**2)
     real(r8) :: rnir                ! reflected solar radiation [nir] (W/m**2)
     real(r8) :: rvis                ! reflected solar radiation [vis] (W/m**2)
     real(r8) :: laifra              ! leaf area fraction of canopy
     real(r8) :: trd                 ! transmitted solar radiation: direct (W/m**2)
     real(r8) :: tri                 ! transmitted solar radiation: diffuse (W/m**2)
     real(r8) :: cad(lbp:ubp,numrad) ! direct beam absorbed by canopy (W/m**2)
     real(r8) :: cai(lbp:ubp,numrad) ! diffuse radiation absorbed by canopy (W/m**2)
     real(r8) :: vai(lbp:ubp)        ! total leaf area index + stem area index, one sided
     real(r8) :: ext                 ! optical depth direct beam per unit LAI+SAI
     real(r8) :: t1, t2              ! temporary variables
     real(r8) :: cosz
     integer  :: local_secp1         ! seconds into current date in local time
     real(r8) :: dtime               ! land model time step (sec)
     integer  :: year,month,day,secs !  calendar info for current time step
!------------------------------------------------------------------------------

     ! Assign local pointers to multi-level derived type members (gridcell level)

     londeg        => clm3%g%londeg
     latdeg        => clm3%g%latdeg
     forc_solad    => clm_a2l%forc_solad
     forc_solai    => clm_a2l%forc_solai

     ! Assign local pointers to multi-level derived type members (column level)

     albgrd        => clm3%g%l%c%cps%albgrd
     albgri        => clm3%g%l%c%cps%albgri
     coszen        => clm3%g%l%c%cps%coszen

     ! Assign local pointers to derived type members (pft-level)

     ivt           => clm3%g%l%c%p%itype
     pcolumn       => clm3%g%l%c%p%column
     pgridcell     => clm3%g%l%c%p%gridcell
     pwtgcell      => clm3%g%l%c%p%wtgcell
     elai          => clm3%g%l%c%p%pps%elai
     esai          => clm3%g%l%c%p%pps%esai
     slasun        => clm3%g%l%c%p%pps%slasun
     slasha        => clm3%g%l%c%p%pps%slasha
     gdir          => clm3%g%l%c%p%pps%gdir
     omega         => clm3%g%l%c%p%pps%omega
     laisun        => clm3%g%l%c%p%pps%laisun
     laisha        => clm3%g%l%c%p%pps%laisha
     fabd          => clm3%g%l%c%p%pps%fabd
     fabi          => clm3%g%l%c%p%pps%fabi
     ftdd          => clm3%g%l%c%p%pps%ftdd
     ftid          => clm3%g%l%c%p%pps%ftid
     ftii          => clm3%g%l%c%p%pps%ftii
     albd          => clm3%g%l%c%p%pps%albd
     albi          => clm3%g%l%c%p%pps%albi
     fsun          => clm3%g%l%c%p%pps%fsun
     sabg          => clm3%g%l%c%p%pef%sabg
     sabv          => clm3%g%l%c%p%pef%sabv
     fsa           => clm3%g%l%c%p%pef%fsa
     fsr           => clm3%g%l%c%p%pef%fsr
     parsun        => clm3%g%l%c%p%pef%parsun
     parsha        => clm3%g%l%c%p%pef%parsha
     fsds_vis_d    => clm3%g%l%c%p%pef%fsds_vis_d
     fsds_nir_d    => clm3%g%l%c%p%pef%fsds_nir_d
     fsds_vis_i    => clm3%g%l%c%p%pef%fsds_vis_i
     fsds_nir_i    => clm3%g%l%c%p%pef%fsds_nir_i
     fsr_vis_d     => clm3%g%l%c%p%pef%fsr_vis_d
     fsr_nir_d     => clm3%g%l%c%p%pef%fsr_nir_d
     fsr_vis_i     => clm3%g%l%c%p%pef%fsr_vis_i
     fsr_nir_i     => clm3%g%l%c%p%pef%fsr_nir_i
     fsds_vis_d_ln => clm3%g%l%c%p%pef%fsds_vis_d_ln
     fsds_nir_d_ln => clm3%g%l%c%p%pef%fsds_nir_d_ln
     fsr_vis_d_ln  => clm3%g%l%c%p%pef%fsr_vis_d_ln
     fsr_nir_d_ln  => clm3%g%l%c%p%pef%fsr_nir_d_ln
     eff_kid =>       clm3%g%l%c%p%pps%eff_kid
     eff_kii =>       clm3%g%l%c%p%pps%eff_kii
     sun_faid =>      clm3%g%l%c%p%pps%sun_faid
     sun_faii =>      clm3%g%l%c%p%pps%sun_faii
     sha_faid =>      clm3%g%l%c%p%pps%sha_faid
     sha_faii =>      clm3%g%l%c%p%pps%sha_faii
     sun_add =>       clm3%g%l%c%p%pef%sun_add
     tot_aid =>       clm3%g%l%c%p%pef%tot_aid
     sun_aid =>       clm3%g%l%c%p%pef%sun_aid
     sun_aii =>       clm3%g%l%c%p%pef%sun_aii
     sha_aid =>       clm3%g%l%c%p%pef%sha_aid
     sha_aii =>       clm3%g%l%c%p%pef%sha_aii
     sun_atot =>      clm3%g%l%c%p%pef%sun_atot
     sha_atot =>      clm3%g%l%c%p%pef%sha_atot
     sun_alf =>       clm3%g%l%c%p%pef%sun_alf
     sha_alf =>       clm3%g%l%c%p%pef%sha_alf
     sun_aperlai =>   clm3%g%l%c%p%pef%sun_aperlai
     sha_aperlai =>   clm3%g%l%c%p%pef%sha_aperlai
     
     ! Assign local pointers to derived type members (ecophysiological)

     slatop        => pftcon%slatop
     dsladlai      => pftcon%dsladlai
     
     ! Determine seconds off current time step
     
     dtime = get_step_size()
     call get_curr_date (year, month, day, secs)

     ! Determine fluxes

!dir$ concurrent
!cdir nodep
     do p = lbp,ubp
        if (pwtgcell(p)>0._r8) then
           sabg(p) = 0._r8
           sabv(p) = 0._r8
           fsa(p)  = 0._r8
        end if
     end do 

     ! Loop over pfts to calculate fsun, etc
!dir$ concurrent
!cdir nodep
     do p = lbp,ubp
        if (pwtgcell(p)>0._r8) then
           c = pcolumn(p)
           g = pgridcell(p)
        
           vai(p) = elai(p) + esai(p)
           if (coszen(c) > 0._r8 .and. elai(p) > 0._r8 .and. gdir(p) > 0._r8) then
              cosz = max(0.001_r8, coszen(c))
              ext = gdir(p)/cosz
              t1 = min(ext*elai(p), 40.0_r8)
              t2 = exp(-t1)
              fsun(p) = (1._r8-t2)/t1
              
              ! new control on low lai, to avoid numerical problems in
              ! calculation of slasun, slasha
              ! PET: 2/29/04
              
              if (elai(p) > 0.01_r8) then
                 laisun(p) = elai(p)*fsun(p)
                 laisha(p) = elai(p)*(1._r8-fsun(p))
                 
                 ! calculate the average specific leaf area for sunlit and shaded
                 ! canopies, when effective LAI > 0
                 slasun(p) = (t2*dsladlai(ivt(p))*ext*elai(p) + &
                              t2*dsladlai(ivt(p)) + &
                              t2*slatop(ivt(p))*ext - &
                              dsladlai(ivt(p)) - &
                              slatop(ivt(p))*ext) / &
                              (ext*(t2-1._r8))
                 slasha(p) = ((slatop(ivt(p)) + &
                             (dsladlai(ivt(p)) * elai(p)/2.0_r8)) * elai(p) - &
                             laisun(p)*slasun(p)) / laisha(p)
              else
                 ! special case for low elai
                 fsun(p) = 1._r8
                 laisun(p) = elai(p)
                 laisha(p) = 0._r8
                 slasun(p) = slatop(ivt(p))
                 slasha(p) = 0._r8
              end if
           else
              fsun(p)   = 0._r8
              laisun(p) = 0._r8
              laisha(p) = elai(p)
              slasun(p) = 0._r8
              slasha(p) = 0._r8
           end if
        end if
     end do
        
     ! Loop over nband wavebands
     do ib = 1, nband
!dir$ concurrent
!cdir nodep
        do p = lbp,ubp
           if (pwtgcell(p)>0._r8) then
              c = pcolumn(p)
              g = pgridcell(p)
              
              ! Absorbed by canopy
              
              cad(p,ib) = forc_solad(g,ib)*fabd(p,ib)
              cai(p,ib) = forc_solai(g,ib)*fabi(p,ib)
              sabv(p) = sabv(p) + cad(p,ib) + cai(p,ib)
              fsa(p)  = fsa(p)  + cad(p,ib) + cai(p,ib)
              
              ! Transmitted = solar fluxes incident on ground
              
              trd = forc_solad(g,ib)*ftdd(p,ib)
              tri = forc_solad(g,ib)*ftid(p,ib) + forc_solai(g,ib)*ftii(p,ib)
              
              ! Solar radiation absorbed by ground surface
              
              abs = trd*(1._r8-albgrd(c,ib)) + tri*(1._r8-albgri(c,ib))
              sabg(p) = sabg(p) + abs
              fsa(p)  = fsa(p)  + abs
              
              ! New sunlit.shaded canopy algorithm
              
              if (coszen(c) > 0._r8 .and. elai(p) > 0._r8) then
                 
                 ! 1. calculate flux of direct beam radiation absorbed in the 
                 ! sunlit canopy as direct (sun_add), and the flux of direct
                 ! beam radiation absorbed in the total canopy as indirect
                 
                 sun_add(p,ib) = forc_solad(g,ib) * (1._r8-ftdd(p,ib)) * (1._r8-omega(p,ib))
                 tot_aid(p,ib) = (forc_solad(g,ib) * fabd(p,ib)) - sun_add(p,ib)
                 
                 ! the following constraint set to catch round-off level errors
                 ! that can cause negative tot_aid
                 
                 tot_aid(p,ib) = max(tot_aid(p,ib), 0._r8)
                 
                 ! 2. calculate the effective extinction coefficients for indirect
                 ! transmission originating from direct and indirect streams,
                 ! using ftid and ftii
                 
                 !eff_kid(p,ib) = -(log(ftid(p,ib)))/vai(p)
                 !eff_kii(p,ib) = -(log(ftii(p,ib)))/vai(p)
                 
                 ! 3. calculate the fraction of indirect radiation being absorbed 
                 ! in the sunlit and shaded canopy fraction. Some of this indirect originates in
                 ! the direct beam and some originates in the indirect beam.

                 !sun_faid(p,ib) = 1.-exp(-eff_kid(p,ib) * vaisun(p))
                 !sun_faii(p,ib) = 1.-exp(-eff_kii(p,ib) * vaisun(p))
                 sun_faid(p,ib) = fsun(p)
                 sun_faii(p,ib) = fsun(p)
                 sha_faid(p,ib) = 1._r8-sun_faid(p,ib)
                 sha_faii(p,ib) = 1._r8-sun_faii(p,ib)

                 ! 4. calculate the total indirect flux absorbed by the sunlit
                 ! and shaded canopy based on these fractions and the fabd and
                 ! fabi from surface albedo calculations

                 sun_aid(p,ib) = tot_aid(p,ib) * sun_faid(p,ib)
                 sun_aii(p,ib) = forc_solai(g,ib)*fabi(p,ib)*sun_faii(p,ib)
                 sha_aid(p,ib) = tot_aid(p,ib) * sha_faid(p,ib)
                 sha_aii(p,ib) = forc_solai(g,ib)*fabi(p,ib)*sha_faii(p,ib)
                 
                 ! 5. calculate the total flux absorbed in the sunlit and shaded
                 ! canopy as the sum of these terms
                 
                 sun_atot(p,ib) = sun_add(p,ib) + sun_aid(p,ib) + sun_aii(p,ib)
                 sha_atot(p,ib) = sha_aid(p,ib) + sha_aii(p,ib)
                 
                 ! 6. calculate the total flux absorbed by leaves in the sunlit
                 ! and shaded canopies
                 
                 laifra = elai(p)/vai(p)
                 sun_alf(p,ib) = sun_atot(p,ib) * laifra
                 sha_alf(p,ib) = sha_atot(p,ib) * laifra
                 
                 ! 7. calculate the fluxes per unit lai in the sunlit and shaded
                 ! canopies
                 
                 if (laisun(p) > 0._r8) then
                    sun_aperlai(p,ib) = sun_alf(p,ib)/laisun(p)
                 else
                    sun_aperlai(p,ib) = 0._r8
                 endif
                 if (laisha(p) > 0._r8) then
                    sha_aperlai(p,ib) = sha_alf(p,ib)/laisha(p)
                 else
                    sha_aperlai(p,ib) = 0._r8
                 endif
             
              else   ! coszen = 0 or elai = 0
                 
                 sun_add(p,ib)     = 0._r8
                 tot_aid(p,ib)     = 0._r8
                 eff_kid(p,ib)     = 0._r8
                 eff_kii(p,ib)     = 0._r8
                 sun_faid(p,ib)    = 0._r8
                 sun_faii(p,ib)    = 0._r8
                 sha_faid(p,ib)    = 0._r8
                 sha_faii(p,ib)    = 0._r8
                 sun_aid(p,ib)     = 0._r8
                 sun_aii(p,ib)     = 0._r8
                 sha_aid(p,ib)     = 0._r8
                 sha_aii(p,ib)     = 0._r8
                 sun_atot(p,ib)    = 0._r8
                 sha_atot(p,ib)    = 0._r8
                 sun_alf(p,ib)     = 0._r8
                 sha_alf(p,ib)     = 0._r8
                 sun_aperlai(p,ib) = 0._r8
                 sha_aperlai(p,ib) = 0._r8
                 
              end if
           end if
        end do ! end of pft loop
     end do ! end nbands loop   

!dir$ concurrent
!cdir nodep
     do p = lbp,ubp
        if (pwtgcell(p)>0._r8) then
           g = pgridcell(p)
        
           ! Final step of new sunlit/shaded canopy algorithm
           ! 8. calculate the total and per-unit-lai fluxes for PAR in the
           ! sunlit and shaded canopy leaf fractions
           
           parsun(p) = sun_aperlai(p,1)
           parsha(p) = sha_aperlai(p,1)
           
           ! The following code is duplicated from SurfaceRadiation
           ! NDVI and reflected solar radiation
           
           rvis = albd(p,1)*forc_solad(g,1) + albi(p,1)*forc_solai(g,1)
           rnir = albd(p,2)*forc_solad(g,2) + albi(p,2)*forc_solai(g,2)
           fsr(p) = rvis + rnir
           
           fsds_vis_d(p) = forc_solad(g,1)
           fsds_nir_d(p) = forc_solad(g,2)
           fsds_vis_i(p) = forc_solai(g,1)
           fsds_nir_i(p) = forc_solai(g,2)
           fsr_vis_d(p)  = albd(p,1)*forc_solad(g,1)
           fsr_nir_d(p)  = albd(p,2)*forc_solad(g,2)
           fsr_vis_i(p)  = albi(p,1)*forc_solai(g,1)
           fsr_nir_i(p)  = albi(p,2)*forc_solai(g,2)
           
           local_secp1 = secs + nint((londeg(g)/15._r8*3600._r8)/dtime)*dtime
           local_secp1 = mod(local_secp1,86400)
           if (local_secp1 == 43200) then
              fsds_vis_d_ln(p) = forc_solad(g,1)
              fsds_nir_d_ln(p) = forc_solad(g,2)
              fsr_vis_d_ln(p) = albd(p,1)*forc_solad(g,1)
              fsr_nir_d_ln(p) = albd(p,2)*forc_solad(g,2)
           else
              fsds_vis_d_ln(p) = spval
              fsds_nir_d_ln(p) = spval
              fsr_vis_d_ln(p) = spval
              fsr_nir_d_ln(p) = spval
           end if
        end if
     end do 

   end subroutine SurfaceRadiation

end module SurfaceRadiationMod
