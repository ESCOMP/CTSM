#include <misc.h>
#include <preproc.h>

module clm_camMod

#if (defined COUP_CAM)

!---------------------------------------------------------------------------
!BOP
!
! !MODULE: clm_camMod
!
! !DESCRIPTION:
! This module provides the interface between the atmosphere land modules.
! If running as part of cam, the land surface model must use the same
! grid as the cam. The land surface model calculates its own net solar
! radiation and net longwave radiation at the surface. The net longwave
! radiation at the surface will differ somewhat from that calculated in the
! atmospheric model because the atm model will use the upward longwave flux
! (or radiative temperature) from the previous time step whereas the land
! surface model uses the flux for the current time step. The net solar
! radiation should equal that calculated in the atmospheric model. If not,
! there is a problem in how the models are coupled.
!
! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
#if (defined SPMD)
  use mpishorthand, only : mpicom
#endif
  use m_AttrVect  , only : AttrVect
!
! !PUBLIC MEMBER FUNCTIONS:
  implicit none
  public clm_camInit                   ! Initialization
  public clm_camRun                    ! run method
  public clm_camFinal                  ! Finalization method
  SAVE
  private                              ! By default make data private
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
!
! !PRIVATE MEMBER FUNCTIONS:
  private clm_camCheckGrid             ! check consistency of cam/clm grid
!
! !PRIVATE VARIABLES
  type(AttrVect) :: a2c_a, l2c_a
  type(AttrVect) :: a2c_l, l2c_l 
  type(AttrVect) :: c2a_a, i2c_a, o2c_a
  type(AttrVect) :: c2l_l 
  logical :: noland = .false.          ! Flag if no land points here
!---------------------------------------------------------------------------

contains

!---------------------------------------------------------------------------
!BOP
!
! !IROUTINE: clm_camInit
!
! !INTERFACE:
  subroutine clm_camInit( lnd_out, atm_in )
!
! !DESCRIPTION:
! Initialize land surface model and obtain relevant atmospheric model arrays
! back from (i.e. albedos, surface temperature and snow cover over land).
!
! !USES:
    ! cam uses
    use ppgrid          , only : begchunk, endchunk
    use radiation       , only : radiation_get
    use camsrfexch_types, only : srfflx_parm, srfflx_state, srfcomp2hub_alloc
    use time_manager    , only : get_nstep
    use filenames       , only : mss_irt, caseid
    use history         , only : ctitle, inithist
    ! clm uses
    use clm_atmlnd      , only : clm_l2a, atm_l2a
    use clm_atmlnd      , only : gridmap_l2a, clm_mapl2a
    use domainMod       , only : adomain
    use clm_comp        , only : clm_init1, clm_init2
    use clm_varctl      , only : cam_caseid, cam_ctitle, cam_irad, cam_nsrest, &
                                 cam_crtinic, cam_irt, finidat       
    ! mct uses
    use MCT_lnd_comp
    use MCT_atm_comp
    use MCT_atmlnd_cpl
#ifdef SCAM
#include <max.h>
    use scamMod         , only : switch, have_tg, isrestart
#endif
#include <comctl.h>
#include <comsol.h>
!
! !ARGUMENTS:
    type(srfflx_parm) , pointer    :: lnd_out(:)
    type(srfflx_state), intent(in) :: atm_in(begchunk:endchunk)
!
! !LOCAL VARIABLES:
    integer  :: i,j         ! indices
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
!-----------------------------------------------------------------------

    ! Initialize memeory for lnd_out (note that this currently cannot
    ! be done in MCT_atmhub_lndImportInit since that is only called at
    ! nstep = 0

    call srfcomp2hub_alloc( lnd_out )

    !=============================================================
    ! Determine if must return now
    !=============================================================

    if (adiabatic .or. ideal_phys .or. aqua_planet) then
       noland = .true.
       return
    end if

    !=============================================================
    ! Set clm namelist variables that must be same as cam
    !=============================================================

    call radiation_get(iradsw_out=cam_irad)
    cam_caseid  = caseid
    cam_ctitle  = ctitle
    cam_nsrest  = nsrest
    cam_crtinic = inithist
    cam_irt     = mss_irt

    !=============================================================
    ! Determine clm orbital parameters to those of cam
    !=============================================================

    call clm_setorb( eccen, obliqr, lambm0, mvelpp )

    !=============================================================
    ! Initialize clm
    !=============================================================

    call clm_init1()

#if ( defined SCAM )
    if (switch(CRM_SW+1)) noland = .true.
    if (adomain%frac(1,1)==0) noland = .true.
    if (noland) return
#endif

    call clm_init2()

    !=============================================================
    ! Initialize MCT
    !=============================================================

    ! Initialize atmospere attribute vectors (this needs to be moved)
    ! Not all of these attribute vectors are currently being used 
    
    call MCT_atm_AvInit( a2c_a, c2a_a, l2c_a, i2c_a, o2c_a )

    ! Initialize land attribute vectors
    ! Not all of these attribute vectors are currently being used 

    call MCT_lnd_AvInit( a2c_l, c2l_l, l2c_l )

    ! Initialize atm->land and land->atm couplers

    call MCT_atm2lnd_init()
    call MCT_lnd2atm_init()

    !=============================================================
    ! Send info back to application driver at time step zero
    !=============================================================

    if (get_nstep() == 0) then
#ifdef TIMING_BARRIERS
       call t_startf('sync_cl2ch_ini')
       call mpi_barrier(mpicom)
       call t_stopf('sync_cl2ch_ini')
#endif
       call t_startf('clump2chunkini')
       call clm_mapl2a(clm_l2a, atm_l2a, gridmap_l2a)
       call MCT_lnd_ExportInit(atm_l2a, l2c_l)
!tcxz       call MCT_lnd_ExportInit(clm_l2a, l2c_l)
       call MCT_lnd2atm( l2c_l, l2c_a ) 
       call MCT_atmhub_lndImportInit( l2c_a, lnd_out ) 
       call t_stopf('clump2chunkini')
    endif

    !=============================================================
    ! Determine consistency with cam grid info 
    !=============================================================

    ! only need to do this on  master processor (note that cam latitudes
    ! and longitudes are computed each time by the cam model at startup 

    call clm_camCheckGrid( lnd_out, atm_in )

    !=============================================================
    ! Determine if have land point 
    !=============================================================

    noland=.true.
    do j = 1,adomain%nj
       do i = 1,adomain%ni
          if (adomain%frac(i,j) > 0._r8) noland = .false.
       end do
    end do

  end subroutine clm_camInit

!---------------------------------------------------------------------------
!BOP
!
! !IROUTINE: clm_camRun
!
! !INTERFACE:
  subroutine clm_camRun( atm_out, lnd_out, rstwr )
!
! !DESCRIPTION:
! Pack data to be sent to land model into a single array.  Send data to
! land model and call land model driver.  Receive data back from land
! model in a single array.  Unpack this data into component arrays.
! NOTE: component arrays are contained in module comsrf.  When coupling to
! an atmospheric model: solar radiation depends on surface albedos from
! the previous time step (based on current surface conditions and solar
! zenith angle for next time step).  Longwave radiation depends on upward
! longwave flux from previous time step.
!
! !USES:
    use ppgrid          , only : begchunk, endchunk
    use camsrfexch_types, only : srfflx_parm, surface_state
    use clm_atmlnd      , only : clm_a2l, clm_l2a, atm_a2l, atm_l2a
    use clm_atmlnd      , only : gridmap_l2a, clm_mapl2a
    use clm_atmlnd      , only : gridmap_a2l, clm_mapa2l
    use clm_comp        , only : clm_run1, clm_run2
    use MCT_lnd_comp 
    use MCT_atm_comp
    use MCT_atmlnd_cpl
!
! !ARGUMENTS:
    type(surface_state), intent(inout) :: atm_out(begchunk:endchunk)
    type(srfflx_parm)  , intent(inout) :: lnd_out(begchunk:endchunk)
    logical,             intent(in)    :: rstwr    ! true => write restart file this step
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
!---------------------------------------------------------------------------

    if (noland) return

    ! a->l coupling

#ifdef TIMING_BARRIERS
    call t_startf('sync_chnk2clmp')
    call mpi_barrier(mpicom)
    call t_stopf('sync_chnk2clmp')
#endif
    call t_startf('chunk2clump')
    call MCT_atm_Export( atm_out, a2c_a )       
    call MCT_atm2lnd( a2c_a, a2c_l )
!tcxz    call MCT_lnd_Import( clm_a2l, a2c_l )
    call MCT_lnd_Import( atm_a2l, a2c_l )
    call clm_mapa2l(atm_a2l, clm_a2l, gridmap_a2l)
    call t_stopf('chunk2clump')
    
    ! Run clm

    call clm_run1( )
    call clm_run2(rstwr)

    ! l->a coupling

#ifdef TIMING_BARRIERS
    call t_startf('sync_clmp2chnk')
    call mpi_barrier(mpicom)
    call t_stopf('sync_clmp2chnk')
#endif
    call t_startf('clump2chunk')
    call clm_mapl2a(clm_l2a, atm_l2a, gridmap_l2a)
    call MCT_lnd_Export(atm_l2a, l2c_l)
!tcxz    call MCT_lnd_Export(clm_l2a, l2c_l)
    call MCT_lnd2atm( l2c_l, l2c_a )
    call MCT_atmhub_lndImport( l2c_a, lnd_out)
    call t_stopf('clump2chunk')
#ifdef TIMING_BARRIERS
    call t_startf ('sync_after_lnd')
    call mpibarrier (mpicom)
    call t_stopf ('sync_after_lnd')
#endif

  end subroutine clm_camRun

!---------------------------------------------------------------------------
!BOP
!
! !IROUTINE: clm_camFinal
!
! !INTERFACE:
  subroutine clm_camFinal(lnd_out)
!
! !DESCRIPTION:
! Finalize land surface model
!
! !USES:
    use camsrfexch_types, only : srfflx_parm
!------------------------------------------------------------------------------
!BOP
!
! !ARGUMENTS:
    type(srfflx_parm), pointer :: lnd_out(:)
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
!---------------------------------------------------------------------------

    deallocate( lnd_out )

  end subroutine clm_camFinal

!---------------------------------------------------------------------------
!BOP
!
! !IROUTINE: clm_camCheckGrid
!
! !INTERFACE:
  subroutine clm_camCheckGrid( lnd_out, atm_in )
!
! !DESCRIPTION:
! Check that cam grid is consistent with clm grid read in from clm surface 
! dataset
!
! !USES:
    use ppgrid          , only : begchunk, endchunk
    use shr_const_mod   , only : SHR_CONST_PI
    use commap          , only : clat, londeg
    use domainMod       , only : adomain
    use rgrid           , only : nlon
    use pmgrid          , only : plon, plat
    use ppgrid          , only : pcols
    use phys_grid       , only : gather_chunk_to_field
    use camsrfexch_types, only : srfflx_state, srfflx_parm
    use spmdMod         , only : masterproc
    use abortutils      , only : endrun
!
! !ARGUMENTS:
    type(srfflx_parm) , pointer    :: lnd_out(:) 
    type(srfflx_state), intent(in) :: atm_in(begchunk:endchunk)
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein 2005-05-14
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: i,j                     ! indices
    integer :: ext_numlon(plat)        ! ext number of longitudes
    real(r8):: ext_lonc(plon,plat)   ! ext lon values
    real(r8):: ext_latc(plon,plat)   ! ext lat values
    real(r8):: ext_landfrac(plon,plat) ! ext fractional land
    integer :: ext_landmask(plon,plat) ! ext land mask
    real(r8):: cam_landfrac(pcols,begchunk:endchunk) ! local landfrac
!---------------------------------------------------------------------------

    ! Determine cam grid info

#ifdef TIMING_BARRIERS
    call t_startf('sync_gather_landfrac')
    call mpi_barrier(mpicom)
    call t_stopf('sync_gather_landfrac')
#endif

    do i = begchunk, endchunk
       lnd_out(i)%areafrac(:pcols) = atm_in(i)%landfrac(:pcols)
    end do

    ! Determine consistency of cam and clm grid - all processors

    do j = 1,plat
       ext_numlon(j) = nlon(j)	
       if (ext_numlon(j) /= adomain%ni) then
          write(6,*)'clm_camInit error: CAM numlon array not consistent'
          call endrun()
       end if
       do i = 1,nlon(j)
          ext_lonc(i,j) = londeg(i,j)
          ext_latc(i,j) = (180._r8/SHR_CONST_PI)*clat(j)
          if ( abs(ext_latc(i,j)-adomain%latc(i,j)) > 1.e-12_r8 ) then
             write(6,*)'clm_camInit error: CAM latitude ',ext_latc(i,j),' and clm input latitude ', &
                  adomain%latc(i,j),' has difference too large at i,j= ',i,j
             call endrun()
          end if
          if ( abs(ext_lonc(i,j)-adomain%lonc(i,j)) > 1.e-12_r8 ) then
             write(6,*)'clm_camInit error: CAM longitude ',ext_lonc(i,j),' and clm input longitude ', &
                  adomain%lonc(i,j),' has difference too large at i,j= ',i,j
             call endrun()
          end if
       end do
    end do

    ! Determine consistency of cam and clm landfrac/landmask - masterproc only
    ! First, gather global landfrac on master processor only

    do i = begchunk, endchunk
       cam_landfrac(:pcols,i) = atm_in(i)%landfrac(:pcols)
    end do
    call gather_chunk_to_field(1, 1, 1, plon, cam_landfrac, ext_landfrac)

    if (masterproc) then
       do j = 1,plat
          do i = 1,nlon(j)
             if (ext_landfrac(i,j) > 0._r8) then
                ext_landmask(i,j) = 1
             else
                ext_landmask(i,j) = 0
             endif
             if (ext_landmask(i,j) /= adomain%mask(i,j)) then
                write(6,*)'clm_camInit error: CAM land mask different from surface dataset at i,j= ',i,j
                call endrun()
             end if
             if (ext_landfrac(i,j) /= adomain%frac(i,j)) then
                write(6,*)'clm_camInit error: CAM fractional land differs from surface dataset at i,j= ',i,j
                call endrun()
             end if
          end do
       end do
    end if
       
  end subroutine clm_camCheckGrid

!---------------------------------------------------------------------------
!BOP
!
! !IROUTINE: clm_setorb
!
! !INTERFACE:
  subroutine clm_setorb(clm_eccen, clm_obliqr, clm_lambm0, clm_mvelpp)
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

  end subroutine clm_setorb

#endif

end module clm_camMod
