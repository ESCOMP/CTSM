module SatellitePhenologyMod

#include "shr_assert.h"

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! CLM Satelitte Phenology model (SP) ecosystem dynamics (phenology, vegetation).
  ! Allow some subroutines to be used by the CLM Carbon Nitrogen model (CLMCN)
  ! so that DryDeposition code can get estimates of LAI differences between months.
  !
  ! !USES:
  use shr_kind_mod , only : r8 => shr_kind_r8, CS => shr_kind_CS
  use shr_log_mod  , only : errMsg => shr_log_errMsg
  use decompMod    , only : bounds_type
  use abortutils   , only : endrun
  use clm_varctl   , only : iulog, use_lai_streams
  use perf_mod     , only : t_startf, t_stopf
  use spmdMod      , only : masterproc, mpicom, iam
  use laiStreamMod , only : lai_init, lai_advance, lai_interp
  use ncdio_pio
  !
  ! !PUBLIC TYPES:
  implicit none
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: CalcSatellitePhenologyTimeInterp ! put the data into the correct format
  public :: UpdateSatellitePhenologyCanopy ! CLM(BGC)-SP phenology and vegetation
  public :: SatellitePhenologyInit ! Dynamically allocate memory
  public :: interpMonthlyVeg       ! interpolate monthly vegetation data
  public :: readAnnualVegetation   ! Read in annual vegetation (needed for Dry-deposition)
  !
  ! !PRIVATE MEMBER FUNCTIONS:
  private :: readMonthlyVegetation   ! read monthly vegetation data for two months
  !
  ! !PRIVATE MEMBER DATA:
  integer                     :: InterpMonths1      ! saved month index
  real(r8)                    :: timwt(2)           ! time weights for month 1 and month 2
  real(r8), allocatable       :: mlai2t(:,:)        ! lai for interpolation (2 months)
  real(r8), allocatable       :: msai2t(:,:)        ! sai for interpolation (2 months)
  real(r8), allocatable       :: mhvt2t(:,:)        ! top vegetation height for interpolation (2 months)
  real(r8), allocatable       :: mhvb2t(:,:)        ! bottom vegetation height for interpolation(2 months)

  character(len=*), parameter :: sourcefile = &
       __FILE__

!==============================================================================
contains
!==============================================================================

  subroutine SatellitePhenologyInit (bounds)
    !
    ! !DESCRIPTION:
    ! Dynamically allocate memory and set to signaling NaN.
    !
    ! !USES:
    use shr_infnan_mod, only : nan => shr_infnan_nan, assignment(=)
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: ier    ! error code
    !-----------------------------------------------------------------------

    InterpMonths1 = -999  ! saved month index

    ier = 0
    if (.not.allocated(mlai2t)) then
       allocate (mlai2t(bounds%begp:bounds%endp,2), &
                 msai2t(bounds%begp:bounds%endp,2), &
                 mhvt2t(bounds%begp:bounds%endp,2), &
                 mhvb2t(bounds%begp:bounds%endp,2), stat=ier)
    end if
    if (ier /= 0) then
       write(iulog,*) 'EcosystemDynini allocation error'
       call endrun(msg=errMsg(sourcefile, __LINE__))
    end if

    mlai2t(bounds%begp : bounds%endp, :) = nan
    msai2t(bounds%begp : bounds%endp, :) = nan
    mhvt2t(bounds%begp : bounds%endp, :) = nan
    mhvb2t(bounds%begp : bounds%endp, :) = nan

    if (use_lai_streams) then
       call lai_init(bounds)
    endif

  end subroutine SatellitePhenologyInit

  !================================================================
  subroutine CalcSatellitePhenologyTimeInterp(bounds, num_filter, filter, canopystate_inst)
    !
    ! !DESCRIPTION:
    ! Ecosystem dynamics: phenology, vegetation
    ! Calculates leaf areas (tlai, elai),  stem areas (tsai, esai) and height (htop).
    !
    ! !USES:
    use pftconMod               , only : noveg, nbrdlf_dcd_brl_shrub
    use WaterDiagnosticBulkType , only : waterdiagnosticbulk_type
    use CanopyStateType         , only : canopystate_type
    use PatchType               , only : patch

    !
    ! !ARGUMENTS:
    type(bounds_type)              , intent(in)    :: bounds
    integer                        , intent(in)    :: num_filter                        ! number of column points in patch filter
    integer                        , intent(in)    :: filter(bounds%endp-bounds%begp+1) ! patch filter points
    type(canopystate_type)         , intent(inout) :: canopystate_inst
    !
    ! !LOCAL VARIABLES:
    integer  :: fp,p,c                            ! indices
    !-----------------------------------------------------------------------

    associate(                                                           &
         tlai_tinterp        => canopystate_inst%tlai_input_patch    ,   & ! Output: [real(r8) (:) ] one-sided leaf area index, no burying by snow
         tsai_tinterp        => canopystate_inst%tsai_input_patch    ,   & ! Output: [real(r8) (:) ] one-sided stem area index, no burying by snow
         htop_tinterp        => canopystate_inst%htop_input_patch    ,   & ! Output: [real(r8) (:) ] canopy top (m)
         hbot_tinterp        => canopystate_inst%hbot_input_patch        & ! Output: [real(r8) (:) ] canopy bottom (m)
         )

      if (use_lai_streams) then
         call lai_interp(bounds, canopystate_inst)
      endif

      do fp = 1, num_filter
         p = filter(fp)
         c = patch%column(p)

         ! need to update elai and esai only every albedo time step so do not
         ! have any inconsistency in lai and sai between SurfaceAlbedo calls (i.e.,
         ! if albedos are not done every time step).
         ! leaf phenology
         ! Set leaf and stem areas based on day of year
         ! Interpolate leaf area index, stem area index, and vegetation heights
         ! between two monthly
         ! The weights below (timwt(1) and timwt(2)) were obtained by a call to
         ! routine InterpMonthlyVeg in subroutine NCARlsm.
         !                 Field   Monthly Values
         !                -------------------------
         ! leaf area index LAI  <- mlai1 and mlai2
         ! leaf area index SAI  <- msai1 and msai2
         ! top height      HTOP <- mhvt1 and mhvt2
         ! bottom height   HBOT <- mhvb1 and mhvb2

         if (.not. use_lai_streams) then
            tlai_tinterp(p) = timwt(1)*mlai2t(p,1) + timwt(2)*mlai2t(p,2)
         endif

         tsai_tinterp(p) = timwt(1)*msai2t(p,1) + timwt(2)*msai2t(p,2)
         htop_tinterp(p) = timwt(1)*mhvt2t(p,1) + timwt(2)*mhvt2t(p,2)
         hbot_tinterp(p) = timwt(1)*mhvb2t(p,1) + timwt(2)*mhvb2t(p,2)
         
      end do
      end associate
      
   end subroutine CalcSatellitePhenologyTimeInterp
   
   !==============================================================================
   
   subroutine UpdateSatellitePhenologyCanopy(bounds, num_filter, filter, &
    waterdiagnosticbulk_inst, canopystate_inst)
    !
    ! !DESCRIPTION:
    ! Ecosystem dynamics: phenology, vegetation
    ! Sets the canopystate_inst% data structure for non-FATES runs
    !
    ! !USES:
    use pftconMod               , only : noveg, nbrdlf_dcd_brl_shrub
    use WaterDiagnosticBulkType , only : waterdiagnosticbulk_type
    use CanopyStateType         , only : canopystate_type
    use PatchType               , only : patch
    use clm_varctl              , only : use_fates
    
    !
    ! !ARGUMENTS:
    type(bounds_type)              , intent(in)    :: bounds
    integer                        , intent(in)    :: num_filter                        ! number of column points in patch filter
    integer                        , intent(in)    :: filter(bounds%endp-bounds%begp+1) ! patch filter points
    type(waterdiagnosticbulk_type) , intent(in)    :: waterdiagnosticbulk_inst
    type(canopystate_type)         , intent(inout) :: canopystate_inst
    !
    ! !LOCAL VARIABLES:
    integer  :: fp,p,c                            ! indices
    real(r8) :: ol                                ! thickness of canopy layer covered by snow (m)
    real(r8) :: fb                                ! fraction of canopy layer covered by snow
    
    if (use_fates) then 
      write(iulog,*) 'Should not be calling this method when use_fates == .true.'
     call endrun(msg=errMsg(sourcefile, __LINE__))
    end if 
    
    associate(                                                           &
      frac_sno           => waterdiagnosticbulk_inst%frac_sno_col   , & ! Input:  [real(r8) (:) ] fraction of ground covered by snow (0 to 1)
      snow_depth         => waterdiagnosticbulk_inst%snow_depth_col , & ! Input:  [real(r8) (:) ] snow height (m)
      tlai_driver        => canopystate_inst%tlai_input_patch    ,    & ! Input: [real(r8) (:) ] SP driver data for one-sided leaf area index, no burying by snow
      tsai_driver        => canopystate_inst%tsai_input_patch    ,    & ! Input: [real(r8) (:) ] SP driver data for one-sided stem area index, no burying by snow
      tlai               => canopystate_inst%tlai_patch,              & ! Output: [real(r8) (:)]  one-sided leaf area index, no burying by snow
      tsai               => canopystate_inst%tsai_patch,              & ! Output: [real(r8) (:)]  one-sided stem area index, no burying by snow
      elai               => canopystate_inst%elai_patch    ,          & ! Output: [real(r8) (:) ] one-sided leaf area index with burying by snow
      esai               => canopystate_inst%esai_patch    ,          & ! Output: [real(r8) (:) ] one-sided stem area index with burying by snow
      htop_driver        => canopystate_inst%htop_input_patch    ,    & ! Input: [real(r8) (:) ] SP driver data for canopy top (m)
      hbot_driver        => canopystate_inst%hbot_input_patch    ,    & ! Input: [real(r8) (:) ] SP driver data for canopy bottom (m)
      htop               => canopystate_inst%htop_patch    ,          & ! Output: [real(r8) (:) ] canopy top (m)
      hbot               => canopystate_inst%hbot_patch    ,          & ! Output: [real(r8) (:) ] canopy bottom (m)
      frac_veg_nosno_alb => canopystate_inst%frac_veg_nosno_alb_patch & ! Output: [integer  (:) ] fraction of vegetation not covered by snow (0 OR 1) [-]
      )
      
    do fp = 1, num_filter
      
      p = filter(fp)
      c = patch%column(p)
      
      ! for regular CLM (non-FATES), this is just a 1:1 mapping
      tlai(p) = tlai_driver(p)
      tsai(p) = tsai_driver(p)
      htop(p) = htop_driver(p)
      hbot(p) = hbot_driver(p)
      
      ! adjust lai and sai for burying by snow. if exposed lai and sai
      ! are less than 0.05, set equal to zero to prevent numerical
      ! problems associated with very small lai and sai.

      ! snow burial fraction for short vegetation (e.g. grasses, crops) changes with vegetation height
      ! accounts for a 20% bending factor, as used in Lombardozzi et al. (2018) GRL 45(18), 9889-9897

      ! NOTE: The following snow burial code is duplicated in CNVegStructUpdateMod.
      ! Changes in one place should be accompanied by similar changes in the other.

      if (patch%itype(p) > noveg .and. patch%itype(p) <= nbrdlf_dcd_brl_shrub) then
        ol = min(max(snow_depth(c) - hbot(p), 0.0_r8), htop(p) - hbot(p))
        fb = 1._r8 - ol / max(1.e-06_r8, htop(p)-hbot(p))
      else
        fb = 1._r8 - (max(min(snow_depth(c),max(0.05,htop(p)*0.8_r8)),0._r8)/(max(0.05,htop(p)*0.8_r8)))
      endif

      elai(p) = max(tlai(p)*(1.0_r8 - frac_sno(c)) + tlai(p)*fb*frac_sno(c), 0.0_r8)
      esai(p) = max(tsai(p)*(1.0_r8 - frac_sno(c)) + tsai(p)*fb*frac_sno(c), 0.0_r8)
      if (elai(p) < 0.05_r8) elai(p) = 0._r8
      if (esai(p) < 0.05_r8) esai(p) = 0._r8

      ! Fraction of vegetation free of snow
      if ((elai(p) + esai(p)) >= 0.05_r8) then
        frac_veg_nosno_alb(p) = 1
      else
        frac_veg_nosno_alb(p) = 0
      end if

    end do ! end of patch loop

    end associate

  end subroutine UpdateSatellitePhenologyCanopy

  !==============================================================================
  subroutine interpMonthlyVeg (bounds, canopystate_inst)
    !
    ! !DESCRIPTION:
    ! Determine if 2 new months of data are to be read.
    !
    ! !USES:
    use clm_varctl       , only : fsurdat
    use clm_time_manager , only : get_curr_date, get_step_size_real, get_nstep
    use CanopyStateType  , only : canopystate_type
    !
    ! !ARGUMENTS:
    type(bounds_type)      , intent(in)    :: bounds
    type(canopystate_type) , intent(inout) :: canopystate_inst
    !
    ! !LOCAL VARIABLES:
    integer :: kyr         ! year (0, ...) for nstep+1
    integer :: kmo         ! month (1, ..., 12)
    integer :: kda         ! day of month (1, ..., 31)
    integer :: ksec        ! seconds into current date for nstep+1
    real(r8):: dtime       ! land model time step (sec)
    real(r8):: t           ! a fraction: kda/ndaypm
    integer :: it(2)       ! month 1 and month 2 (step 1)
    integer :: months(2)   ! months to be interpolated (1 to 12)
    integer, dimension(12) :: ndaypm= &
         (/31,28,31,30,31,30,31,31,30,31,30,31/) !days per month
    !-----------------------------------------------------------------------

    dtime = get_step_size_real()

    call get_curr_date(kyr, kmo, kda, ksec, offset=int(dtime))

    t = (kda-0.5_r8) / ndaypm(kmo)
    it(1) = t + 0.5_r8
    it(2) = it(1) + 1
    months(1) = kmo + it(1) - 1
    months(2) = kmo + it(2) - 1
    if (months(1) <  1) months(1) = 12
    if (months(2) > 12) months(2) = 1
    timwt(1) = (it(1)+0.5_r8) - t
    timwt(2) = 1._r8-timwt(1)

    if (InterpMonths1 /= months(1)) then
       if (masterproc) then
          write(iulog,*) 'Attempting to read monthly vegetation data .....'
          write(iulog,*) 'nstep = ',get_nstep(),' month = ',kmo,' day = ',kda
       end if
       call t_startf('readMonthlyVeg')
       call readMonthlyVegetation (bounds, fsurdat, months, canopystate_inst)
       InterpMonths1 = months(1)
       call t_stopf('readMonthlyVeg')
    end if

  end subroutine interpMonthlyVeg

  !==============================================================================
  subroutine readAnnualVegetation (bounds, canopystate_inst)
    !
    ! !DESCRIPTION:
    ! read 12 months of veg data for dry deposition
    !
    ! !USES:
    use clm_varpar      , only : maxveg, maxsoil_patches
    use pftconMod       , only : noveg
    use fileutils       , only : getfil
    use clm_varctl      , only : fsurdat
    use domainMod       , only : ldomain
    use clm_varcon      , only : grlnd
    use PatchType       , only : patch
    use CanopyStateType , only : canopystate_type
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds
    type(canopystate_type), intent(inout) :: canopystate_inst
    !
    ! !LOCAL VARIABLES:
    type(file_desc_t) :: ncid             ! netcdf id
    real(r8), pointer :: annlai(:,:)      ! 12 months of monthly lai from input data set
    real(r8), pointer :: mlai(:,:)        ! lai read from input files
    integer :: ier                        ! error code
    integer :: g,k,l,m,n,p                ! indices
    integer :: ni,nj,ns                   ! indices
    integer :: dimid,varid                ! input netCDF id's
    integer :: ntim                       ! number of input data time samples
    integer :: nlon_i                     ! number of input data longitudes
    integer :: nlat_i                     ! number of input data latitudes
    integer :: npft_i                     ! number of input data patch types
    logical :: isgrid2d                   ! true => file is 2d
    character(len=256) :: locfn           ! local file name
    character(len=32) :: subname = 'readAnnualVegetation'
    !-----------------------------------------------------------------------

    annlai    => canopystate_inst%annlai_patch

    ! Determine necessary indices

    allocate(mlai(bounds%begg:bounds%endg,0:maxveg), stat=ier)
    if (ier /= 0) then
       write(iulog,*)subname, 'allocation error '
       call endrun(msg=errMsg(sourcefile, __LINE__))
    end if

    if (masterproc) then
       write (iulog,*) 'Attempting to read annual vegetation data .....'
    end if

    call getfil(fsurdat, locfn, 0)
    call ncd_pio_openfile (ncid, trim(locfn), 0)
    call ncd_inqfdims (ncid, isgrid2d, ni, nj, ns)

    if (ldomain%ns /= ns .or. ldomain%ni /= ni .or. ldomain%nj /= nj) then
       write(iulog,*)trim(subname), 'ldomain and input file do not match dims '
       write(iulog,*)trim(subname), 'ldomain%ni,ni,= ',ldomain%ni,ni
       write(iulog,*)trim(subname), 'ldomain%nj,nj,= ',ldomain%nj,nj
       write(iulog,*)trim(subname), 'ldomain%ns,ns,= ',ldomain%ns,ns
       call endrun(msg=errMsg(sourcefile, __LINE__))
    end if
    call check_dim_size(ncid, 'lsmpft', maxsoil_patches)

    do k=1,12   !! loop over months and read vegetated data

       call ncd_io(ncid=ncid, varname='MONTHLY_LAI', flag='read', data=mlai, &
            dim1name=grlnd, nt=k)

       !! only vegetated patches have nonzero values
       !! Assign lai/sai/hgtt/hgtb to the top [maxsoil_patches] patches
       !! as determined in subroutine surfrd

       do p = bounds%begp,bounds%endp
          g =patch%gridcell(p)
          if (patch%itype(p) /= noveg) then     !! vegetated pft
             do l = 0, maxveg
                if (l == patch%itype(p)) then
                   annlai(k,p) = mlai(g,l)
                end if
             end do
          else                       !! non-vegetated pft
             annlai(k,p) = 0._r8
          end if
       end do   ! end of loop over patches

    enddo ! months loop

    call ncd_pio_closefile(ncid)

    deallocate(mlai)

  endsubroutine readAnnualVegetation

  !==============================================================================
  subroutine readMonthlyVegetation (bounds, fveg, months, canopystate_inst)
    !
    ! !DESCRIPTION:
    ! Read monthly vegetation data for two consec. months.
    !
    ! !USES:
    use clm_varpar       , only : maxveg
    use pftconMod        , only : noveg
    use fileutils        , only : getfil
    use spmdMod          , only : masterproc, mpicom, MPI_REAL8, MPI_INTEGER
    use clm_time_manager , only : get_nstep
    use CanopyStateType  , only : canopystate_type
    use PatchType        , only : patch
    use clm_varcon       , only : grlnd
    use netcdf
    !
    ! !ARGUMENTS:
    type(bounds_type) , intent(in) :: bounds
    character(len=*)  , intent(in) :: fveg      ! file with monthly vegetation data
    integer           , intent(in) :: months(2) ! months to be interpolated (1 to 12)
    type(canopystate_type), intent(inout) :: canopystate_inst
    !
    ! !LOCAL VARIABLES:
    character(len=256) :: locfn           ! local file name
    type(file_desc_t)  :: ncid            ! netcdf id
    integer :: g,n,k,l,m,p,ni,nj,ns       ! indices
    integer :: dimid,varid                ! input netCDF id's
    integer :: ntim                       ! number of input data time samples
    integer :: nlon_i                     ! number of input data longitudes
    integer :: nlat_i                     ! number of input data latitudes
    integer :: npft_i                     ! number of input data patch types
    integer :: ier                        ! error code
    logical :: readvar
    real(r8), pointer :: mlai(:,:)        ! lai read from input files
    real(r8), pointer :: msai(:,:)        ! sai read from input files
    real(r8), pointer :: mhgtt(:,:)       ! top vegetation height
    real(r8), pointer :: mhgtb(:,:)       ! bottom vegetation height
    character(len=32) :: subname = 'readMonthlyVegetation'
    !-----------------------------------------------------------------------

    ! Determine necessary indices

    allocate(&
         mlai(bounds%begg:bounds%endg,0:maxveg), &
         msai(bounds%begg:bounds%endg,0:maxveg), &
         mhgtt(bounds%begg:bounds%endg,0:maxveg), &
         mhgtb(bounds%begg:bounds%endg,0:maxveg), &
         stat=ier)
    if (ier /= 0) then
       write(iulog,*)subname, 'allocation big error '
       call endrun(msg=errMsg(sourcefile, __LINE__))
    end if

    ! ----------------------------------------------------------------------
    ! Open monthly vegetation file
    ! Read data and convert from gridcell to patch data
    ! ----------------------------------------------------------------------

    call getfil(fveg, locfn, 0)
    call ncd_pio_openfile (ncid, trim(locfn), 0)

    do k=1,2   !loop over months and read vegetated data

       call ncd_io(ncid=ncid, varname='MONTHLY_LAI', flag='read', data=mlai, dim1name=grlnd, &
            nt=months(k), readvar=readvar)
       if (.not. readvar) call endrun(msg=' ERROR: MONTHLY_LAI NOT on fveg file'//errMsg(sourcefile, __LINE__))

       call ncd_io(ncid=ncid, varname='MONTHLY_SAI', flag='read', data=msai, dim1name=grlnd, &
            nt=months(k), readvar=readvar)
       if (.not. readvar) call endrun(msg=' ERROR: MONTHLY_SAI NOT on fveg file'//errMsg(sourcefile, __LINE__))

       call ncd_io(ncid=ncid, varname='MONTHLY_HEIGHT_TOP', flag='read', data=mhgtt, dim1name=grlnd, &
            nt=months(k), readvar=readvar)
       if (.not. readvar) call endrun(msg=' ERROR: MONTHLY_HEIGHT_TOP NOT on fveg file'//errMsg(sourcefile, __LINE__))

       call ncd_io(ncid=ncid, varname='MONTHLY_HEIGHT_BOT', flag='read', data=mhgtb, dim1name=grlnd, &
            nt=months(k), readvar=readvar)
       if (.not. readvar) call endrun(msg=' ERROR: MONTHLY_HEIGHT_TOP NOT on fveg file'//errMsg(sourcefile, __LINE__))

       ! Only vegetated patches have nonzero values
       ! Assign lai/sai/hgtt/hgtb to the top [maxsoil_patches] patches
       ! as determined in subroutine surfrd

       do p = bounds%begp,bounds%endp
          g =patch%gridcell(p)
          if (patch%itype(p) /= noveg) then     ! vegetated pft
             do l = 0, maxveg
                if (l == patch%itype(p)) then
                   mlai2t(p,k) = mlai(g,l)
                   msai2t(p,k) = msai(g,l)
                   mhvt2t(p,k) = mhgtt(g,l)
                   mhvb2t(p,k) = mhgtb(g,l)
                end if
             end do
          else                        ! non-vegetated pft
             mlai2t(p,k) = 0._r8
             msai2t(p,k) = 0._r8
             mhvt2t(p,k) = 0._r8
             mhvb2t(p,k) = 0._r8
          end if
       end do   ! end of loop over patches

    end do   ! end of loop over months

    call ncd_pio_closefile(ncid)

    if (masterproc) then
       k = 2
       write(iulog,*) 'Successfully read monthly vegetation data for'
       write(iulog,*) 'month ', months(k)
       write(iulog,*)
    end if

    deallocate(mlai, msai, mhgtt, mhgtb)

    do p = bounds%begp,bounds%endp
       canopystate_inst%mlaidiff_patch(p) = mlai2t(p,1)-mlai2t(p,2)
    enddo

  end subroutine readMonthlyVegetation

end module SatellitePhenologyMod
