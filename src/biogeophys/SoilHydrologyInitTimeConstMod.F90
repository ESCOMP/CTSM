module SoilHydrologyInitTimeConstMod

  !------------------------------------------------------------------------------
  ! DESCRIPTION:
  ! Initialize time constant variables for SoilHydrologyType
  !
  ! !USES
  use shr_kind_mod      , only : r8 => shr_kind_r8
  use shr_log_mod       , only : errMsg => shr_log_errMsg
  use decompMod         , only : bounds_type
  use SoilStateType     , only : soilstate_type
  use SoilHydrologyType , only : soilhydrology_type
  use WaterStateBulkType, only : waterstatebulk_type
  use LandunitType      , only : lun                
  use ColumnType        , only : col                
  use SoilStateInitTimeConstMod, only: organic_max
  !
  implicit none
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: SoilHydrologyInitTimeConst
  public :: readParams
  !
  ! !PRIVATE MEMBER FUNCTIONS:
  private :: initSoilParVIC    ! Convert default CLM soil properties to VIC parameters
  private :: initCLMVICMap     ! Initialize map from VIC to CLM layers
  private :: linear_interp     ! function for linear interperation 
  type, private :: params_type
     real(r8) :: pc            ! Threshold probability for surface water (unitless)
     real(r8) :: om_frac_sf    ! Scale factor for organic matter fraction (unitless)
  end type params_type
  type(params_type), private ::  params_inst

  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  !-----------------------------------------------------------------------
  !
contains

  !-----------------------------------------------------------------------
  subroutine readParams( ncid )
    !
    ! !USES:
    use ncdio_pio, only: file_desc_t
    use paramUtilMod, only: readNcdioScalar
    !
    ! !ARGUMENTS:
    implicit none
    type(file_desc_t),intent(inout) :: ncid   ! pio netCDF file id
    !
    ! !LOCAL VARIABLES:
    character(len=*), parameter :: subname = 'readParams_SoilHydrologyInitTimeConst'
    !--------------------------------------------------------------------

    ! Threshold probability for surface water (unitless)
    call readNcdioScalar(ncid, 'pc', subname, params_inst%pc)
    ! Scale factor for om_frac (unitless)
    call readNcdioScalar(ncid, 'om_frac_sf', subname, params_inst%om_frac_sf)

  end subroutine readParams

  !-----------------------------------------------------------------------
  subroutine SoilHydrologyInitTimeConst(bounds, soilhydrology_inst, soilstate_inst)
    !
    ! !USES:
    use shr_const_mod   , only : shr_const_pi
    use shr_spfn_mod    , only : shr_spfn_erf
    use abortutils      , only : endrun
    use spmdMod         , only : masterproc
    use clm_varctl      , only : fsurdat, iulog, use_vichydro
    use clm_varpar      , only : toplev_equalspace
    use clm_varpar      , only : nlevsoi, nlevgrnd, nlayer, nlayert
    use clm_varcon      , only : dzsoi, spval, nlvic, dzvic, grlnd
    use clm_varcon      , only : aquifer_water_baseline
    use landunit_varcon , only : istwet, istdlak, istice_mec
    use column_varcon   , only : icol_shadewall, icol_road_perv, icol_road_imperv, icol_roof, icol_sunwall
    use fileutils       , only : getfil
    use ncdio_pio       , only : file_desc_t, ncd_io, ncd_pio_openfile, ncd_pio_closefile
    !
    ! !ARGUMENTS:
    type(bounds_type)        , intent(in)    :: bounds                                    
    type(soilhydrology_type) , intent(inout) :: soilhydrology_inst
    type(soilstate_type)     , intent(in)    :: soilstate_inst
    !
    ! !LOCAL VARIABLES:
    integer            :: p,c,j,l,g,lev,nlevs 
    integer            :: ivic,ivicstrt,ivicend   
    real(r8)           :: maxslope, slopemax, minslope
    real(r8)           :: d, fd, dfdd, slope0,slopebeta
    real(r8) ,pointer  :: tslope(:)  
    logical            :: readvar 
    type(file_desc_t)  :: ncid        
    character(len=256) :: locfn       
    real(r8) ,pointer  :: b2d        (:)   ! read in - VIC b  
    real(r8) ,pointer  :: ds2d       (:)   ! read in - VIC Ds 
    real(r8) ,pointer  :: dsmax2d    (:)   ! read in - VIC Dsmax 
    real(r8) ,pointer  :: ws2d       (:)   ! read in - VIC Ws 
    real(r8), pointer  :: sandcol    (:,:) ! column level sand fraction for calculating VIC parameters
    real(r8), pointer  :: claycol    (:,:) ! column level clay fraction for calculating VIC parameters
    real(r8), pointer  :: om_fraccol (:,:) ! column level organic matter fraction for calculating VIC parameters
    !-----------------------------------------------------------------------
    ! Initialize VIC variables

    if (use_vichydro) then

       allocate(b2d        (bounds%begg:bounds%endg))
       allocate(ds2d       (bounds%begg:bounds%endg))
       allocate(dsmax2d    (bounds%begg:bounds%endg))
       allocate(ws2d       (bounds%begg:bounds%endg))
       allocate(sandcol    (bounds%begc:bounds%endc,1:nlevgrnd ))
       allocate(claycol    (bounds%begc:bounds%endc,1:nlevgrnd ))
       allocate(om_fraccol (bounds%begc:bounds%endc,1:nlevgrnd ))

       call getfil (fsurdat, locfn, 0)
       call ncd_pio_openfile (ncid, locfn, 0)
       call ncd_io(ncid=ncid, varname='binfl', flag='read', data=b2d, dim1name=grlnd, readvar=readvar)
       if (.not. readvar) then
          call endrun(msg=' ERROR: binfl NOT on surfdata file'//errMsg(sourcefile, __LINE__))
       end if
       call ncd_io(ncid=ncid, varname='Ds', flag='read', data=ds2d, dim1name=grlnd, readvar=readvar)
       if (.not. readvar) then
          call endrun(msg=' ERROR: Ds NOT on surfdata file'//errMsg(sourcefile, __LINE__))
       end if
       call ncd_io(ncid=ncid, varname='Dsmax', flag='read', data=dsmax2d, dim1name=grlnd, readvar=readvar)
       if (.not. readvar) then
          call endrun(msg=' ERROR: Dsmax NOT on surfdata file'//errMsg(sourcefile, __LINE__))
       end if
       call ncd_io(ncid=ncid, varname='Ws', flag='read', data=ws2d, dim1name=grlnd, readvar=readvar)
       if (.not. readvar) then
          call endrun(msg=' ERROR: Ws NOT on surfdata file'//errMsg(sourcefile, __LINE__))
       end if
       call ncd_pio_closefile(ncid)

       !define the depth of VIC soil layers here
       nlvic(1) = 3
       nlvic(2) = toplev_equalspace - nlvic(1)
       nlvic(3) = nlevsoi - (nlvic(1) + nlvic(2))

       dzvic(:) = 0._r8
       ivicstrt = 1

       do ivic = 1,nlayer
          ivicend = ivicstrt+nlvic(ivic)-1
          do j = ivicstrt,ivicend
             dzvic(ivic) = dzvic(ivic)+dzsoi(j)
          end do
          ivicstrt = ivicend+1
       end do

       do c = bounds%begc, bounds%endc
          g = col%gridcell(c)
          soilhydrology_inst%b_infil_col(c) = b2d(g)
          soilhydrology_inst%ds_col(c)      = ds2d(g)
          soilhydrology_inst%dsmax_col(c)   = dsmax2d(g)
          soilhydrology_inst%Wsvic_col(c)   = ws2d(g)
       end do

       do c = bounds%begc, bounds%endc
          do lev = 1, nlayer
             soilhydrology_inst%ice_col(c,lev)       = spval
             soilhydrology_inst%moist_col(c,lev)     = spval
             soilhydrology_inst%moist_vol_col(c,lev) = spval
             soilhydrology_inst%max_moist_col(c,lev) = spval
             soilhydrology_inst%porosity_col(c,lev)  = spval
             soilhydrology_inst%expt_col(c,lev)      = spval
             soilhydrology_inst%ksat_col(c,lev)      = spval
             soilhydrology_inst%phi_s_col(c,lev)     = spval
             soilhydrology_inst%depth_col(c,lev)     = spval
             sandcol(c,lev)            = spval
             claycol(c,lev)            = spval
             om_fraccol(c,lev)         = spval
          end do
       end do

       do c = bounds%begc, bounds%endc
          g = col%gridcell(c)
          l = col%landunit(c)

          if (lun%itype(l) /= istdlak) then  ! soil columns of both urban and non-urban types
             if (lun%itype(l)==istwet .or. lun%itype(l)==istice_mec) then
                ! do nothing
             else if (lun%urbpoi(l) .and. (col%itype(c) /= icol_road_perv) .and. (col%itype(c) /= icol_road_imperv) )then
                ! do nothing
             else
                do lev = 1,nlevgrnd
                   if ( lev <= nlevsoi )then
                      claycol(c,lev)    = soilstate_inst%cellclay_col(c,lev)
                      sandcol(c,lev)    = soilstate_inst%cellsand_col(c,lev)
                      om_fraccol(c,lev) = min(params_inst%om_frac_sf*soilstate_inst%cellorg_col(c,lev) / organic_max, 1._r8)
                   else
                      claycol(c,lev)    = soilstate_inst%cellclay_col(c,nlevsoi)
                      sandcol(c,lev)    = soilstate_inst%cellsand_col(c,nlevsoi)
                      om_fraccol(c,lev) = 0.0_r8
                   end if

                end do
                soilhydrology_inst%depth_col(c, 1:nlayer) = dzvic
                soilhydrology_inst%depth_col(c, nlayer+1:nlayert) = col%dz(c, nlevsoi+1:nlevgrnd)

                ! create weights to map soil moisture profiles (10 layer) to 3 layers for VIC hydrology, M.Huang
                call initCLMVICMap(c, soilhydrology_inst)
                call initSoilParVIC(c, claycol, sandcol, om_fraccol, soilhydrology_inst)
             end if
          end if ! end of if not lake
       end do ! end of loop over columns

      deallocate(b2d, ds2d, dsmax2d, ws2d)
      deallocate(sandcol, claycol, om_fraccol)

    end if ! end of if use_vichydro

    associate(micro_sigma => col%micro_sigma)
      do c = bounds%begc, bounds%endc
         
         ! determine h2osfc threshold ("fill & spill" concept)
         ! set to zero for no h2osfc (w/frac_infclust =large)
         
         soilhydrology_inst%h2osfc_thresh_col(c) = 0._r8
         if (micro_sigma(c) > 1.e-6_r8 .and. (soilhydrology_inst%h2osfcflag /= 0)) then
            d = 0.0
            do p = 1,4
               fd   = 0.5*(1.0_r8+shr_spfn_erf(d/(micro_sigma(c)*sqrt(2.0)))) - params_inst%pc
               dfdd = exp(-d**2/(2.0*micro_sigma(c)**2))/(micro_sigma(c)*sqrt(2.0*shr_const_pi))
               d    = d - fd/dfdd
            enddo
            soilhydrology_inst%h2osfc_thresh_col(c) = 0.5*d*(1.0_r8+shr_spfn_erf(d/(micro_sigma(c)*sqrt(2.0)))) + &
                 micro_sigma(c)/sqrt(2.0*shr_const_pi)*exp(-d**2/(2.0*micro_sigma(c)**2))         
            soilhydrology_inst%h2osfc_thresh_col(c) = 1.e3_r8 * soilhydrology_inst%h2osfc_thresh_col(c) !convert to mm from meters
         else
            soilhydrology_inst%h2osfc_thresh_col(c) = 0._r8
         endif

         if (soilhydrology_inst%h2osfcflag == 0) then 
            soilhydrology_inst%h2osfc_thresh_col(c) = 0._r8    ! set to zero for no h2osfc (w/frac_infclust =large)
         endif

         ! set decay factor
         soilhydrology_inst%hkdepth_col(c) = 1._r8/2.5_r8

      end do
    end associate

  end subroutine SoilhydrologyInitTimeConst

  !-----------------------------------------------------------------------
  subroutine initSoilParVIC(c, claycol, sandcol, om_fraccol, soilhydrology_inst)
    !
    ! !DESCRIPTION:
    ! Convert default CLM soil properties to VIC parameters
    ! to be used for runoff simulations (added by M. Huang)
    !
    ! !USES:
    use clm_varpar, only : nlevsoi, nlayert, nlayer
    !
    ! !ARGUMENTS:
    integer                  , intent(in)    :: c               ! column index
    real(r8)                 , pointer       :: sandcol(:,:)    ! read in - soil texture: percent sand
    real(r8)                 , pointer       :: claycol(:,:)    ! read in - soil texture: percent clay
    real(r8)                 , pointer       :: om_fraccol(:,:) ! read in - organic matter: kg/m3
    type(soilhydrology_type) , intent(inout) :: soilhydrology_inst

    ! !LOCAL VARIABLES:
    real(r8) :: om_watsat    = 0.9_r8             ! porosity of organic soil
    real(r8) :: om_hksat     = 0.1_r8             ! saturated hydraulic conductivity of organic soil [mm/s]
    real(r8) :: om_tkm       = 0.25_r8            ! thermal conductivity of organic soil (Farouki, 1986) [W/m/K]
    real(r8) :: om_sucsat    = 10.3_r8            ! saturated suction for organic matter (Letts, 2000)
    real(r8) :: om_csol      = 2.5_r8             ! heat capacity of peat soil *10^6 (J/K m3) (Farouki, 1986)
    real(r8) :: om_tkd       = 0.05_r8            ! thermal conductivity of dry organic soil (Farouki, 1981)
    real(r8) :: om_b         = 2.7_r8             ! Clapp Hornberger paramater for oragnic soil (Letts, 2000)
    real(r8) :: om_expt      = 3._r8+2._r8*2.7_r8 ! soil expt for VIC        
    real(r8) :: csol_bedrock = 2.0e6_r8           ! vol. heat capacity of granite/sandstone  J/(m3 K)(Shabbir, 2000)
    real(r8) :: pc           = 0.5_r8             ! percolation threshold
    real(r8) :: pcbeta       = 0.139_r8           ! percolation exponent
    real(r8) :: xksat                             ! maximum hydraulic conductivity of soil [mm/s]
    real(r8) :: perc_frac                         ! "percolating" fraction of organic soil
    real(r8) :: perc_norm                         ! normalize to 1 when 100% organic soil
    real(r8) :: uncon_hksat                       ! series conductivity of mineral/organic soil
    real(r8) :: uncon_frac                        ! fraction of "unconnected" soil
    real(r8) :: temp_sum_frac                     ! sum of node fractions in each VIC layer
    real(r8) :: sandvic(1:nlayert)                ! temporary, weighted averaged sand% for VIC layers
    real(r8) :: clayvic(1:nlayert)                ! temporary, weighted averaged clay% for VIC layers
    real(r8) :: om_fracvic(1:nlayert)             ! temporary, weighted averaged organic matter fract for VIC layers
    integer  :: i, j                              ! indices
    !-------------------------------------------------------------------------------------------

    ! soilhydrology_inst%depth_col(:,:)           Output: layer depth of upper layer(m) 
    ! soilhydrology_inst%vic_clm_fract_col(:,:,:) Output: fraction of VIC layers in CLM layers
    ! soilhydrology_inst%c_param_col(:)           Output: baseflow exponent (Qb)
    ! soilhydrology_inst%expt_col(:,:)            Output: pore-size distribution related paramter(Q12)
    ! soilhydrology_inst%ksat_col(:,:)            Output: Saturated hydrologic conductivity (mm/s)
    ! soilhydrology_inst%phi_s_col(:,:)           Output: soil moisture dissusion parameter
    ! soilhydrology_inst%porosity_col(:,:)        Output: soil porosity
    ! soilhydrology_inst%max_moist_col(:,:)       Output: maximum soil moisture (ice + liq)

    !  map parameters between VIC layers and CLM layers
    soilhydrology_inst%c_param_col(c) = 2._r8

    ! map the CLM layers to VIC layers 
    do i = 1, nlayer      

       sandvic(i)    = 0._r8
       clayvic(i)    = 0._r8   
       om_fracvic(i) = 0._r8  
       temp_sum_frac = 0._r8     
       do j = 1, nlevsoi
          sandvic(i)    = sandvic(i)    + sandcol(c,j)    * soilhydrology_inst%vic_clm_fract_col(c,i,j)
          clayvic(i)    = clayvic(i)    + claycol(c,j)    * soilhydrology_inst%vic_clm_fract_col(c,i,j)
          om_fracvic(i) = om_fracvic(i) + om_fraccol(c,j) * soilhydrology_inst%vic_clm_fract_col(c,i,j) 
          temp_sum_frac = temp_sum_frac +                   soilhydrology_inst%vic_clm_fract_col(c,i,j)
       end do

       !average soil properties, M.Huang, 08/11/2010
       sandvic(i) = sandvic(i)/temp_sum_frac
       clayvic(i) = clayvic(i)/temp_sum_frac
       om_fracvic(i) = om_fracvic(i)/temp_sum_frac

       !make sure sand, clay and om fractions are between 0 and 100% 
       sandvic(i)    = min(100._r8 , sandvic(i))
       clayvic(i)    = min(100._r8 , clayvic(i))
       om_fracvic(i) = min(100._r8 , om_fracvic(i))
       sandvic(i)    = max(0._r8   , sandvic(i))
       clayvic(i)    = max(0._r8   , clayvic(i))
       om_fracvic(i) = max(0._r8   , om_fracvic(i))

       !calculate other parameters based on teh percentages
       soilhydrology_inst%porosity_col(c, i) = 0.489_r8 - 0.00126_r8*sandvic(i)
       soilhydrology_inst%expt_col(c, i)     = 3._r8+ 2._r8*(2.91_r8 + 0.159_r8*clayvic(i))
       xksat = 0.0070556 *( 10.**(-0.884+0.0153*sandvic(i)) )

       !consider organic matter, M.Huang 
       soilhydrology_inst%expt_col(c, i)    = &
            (1._r8 - om_fracvic(i))*soilhydrology_inst%expt_col(c, i)    + om_fracvic(i)*om_expt 
       soilhydrology_inst%porosity_col(c,i) = &
            (1._r8 - om_fracvic(i))*soilhydrology_inst%porosity_col(c,i) + om_watsat*om_fracvic(i) 

       ! perc_frac is zero unless perf_frac greater than percolation threshold
       if (om_fracvic(i) > pc) then
          perc_norm=(1._r8 - pc)**(-pcbeta)
          perc_frac=perc_norm*(om_fracvic(i) - pc)**pcbeta
       else
          perc_frac=0._r8
       endif
       ! uncon_frac is fraction of mineral soil plus fraction of "nonpercolating" organic soil
       uncon_frac=(1._r8-om_fracvic(i))+(1._r8-perc_frac)*om_fracvic(i)

       ! uncon_hksat is series addition of mineral/organic conductivites
       if (om_fracvic(i) < 1._r8) then
          uncon_hksat=uncon_frac/((1._r8-om_fracvic(i))/xksat &
               +((1._r8-perc_frac)*om_fracvic(i))/om_hksat)
       else
          uncon_hksat = 0._r8
       end if

       soilhydrology_inst%ksat_col(c,i)  = &
            uncon_frac*uncon_hksat + (perc_frac*om_fracvic(i))*om_hksat

       soilhydrology_inst%max_moist_col(c,i) = &
            soilhydrology_inst%porosity_col(c,i) * soilhydrology_inst%depth_col(c,i) * 1000._r8 !in mm!

       soilhydrology_inst%phi_s_col(c,i) = &
            -(exp((1.54_r8 - 0.0095_r8*sandvic(i) + &
            0.0063_r8*(100.0_r8-sandvic(i)-clayvic(i)))*log(10.0_r8))*9.8e-5_r8)

    end do ! end of loop over layers

  end subroutine initSoilParVIC

   !-----------------------------------------------------------------------
   subroutine initCLMVICMap(c, soilhydrology_inst)
     !
     ! !DESCRIPTION:
     ! Calculates mapping between CLM and VIC layers
     ! added by AWang, modified by M.Huang for CLM4 
     ! NOTE: in CLM h2osoil_liq unit is kg/m2, in VIC moist is mm
     ! h2osoi_ice is actually water equavlent ice content.
     !
     ! !USES:
     use clm_varpar  , only : nlevsoi, nlayer
     !
     ! !ARGUMENTS:
     integer , intent(in)  :: c
     type(soilhydrology_type), intent(inout) :: soilhydrology_inst
     !
     ! !REVISION HISTORY:
     ! Created by Maoyi Huang
     ! 11/13/2012, Maoyi Huang: rewrite the mapping modules in CLM4VIC 
     !
     ! !LOCAL VARIABLES
     real(r8) :: sum_frac(1:nlayer)                  ! sum of fraction for each layer
     real(r8) :: deltal(1:nlayer+1)                  ! temporary
     real(r8) :: zsum                                ! temporary
     real(r8) :: lsum                                ! temporary
     real(r8) :: temp                                ! temporary
     integer :: i, j, fc
     !-----------------------------------------------------------------------

     associate(                                                    & 
          dz            =>    col%dz    ,                          & ! Input:  [real(r8) (:,:)   ]  layer depth (m)                       
          zi            =>    col%zi    ,                          & ! Input:  [real(r8) (:,:)   ]  interface level below a "z" level (m) 
          z             =>    col%z     ,                          & ! Input:  [real(r8) (:,:)   ]  layer thickness (m)                   

          depth         =>    soilhydrology_inst%depth_col ,       & ! Input:  [real(r8) (:,:)   ]  layer depth of VIC (m)                
          vic_clm_fract =>    soilhydrology_inst%vic_clm_fract_col & ! Output: [real(r8) (:,:,:) ]  fraction of VIC layers in clm layers
          )

       !  set fraction of VIC layer in each CLM layer

       lsum = 0._r8
       do i = 1, nlayer
          deltal(i) = depth(c,i)
       end do
       do i = 1, nlayer
          zsum = 0._r8
          sum_frac(i) = 0._r8
          do j = 1, nlevsoi
             if( (zsum < lsum) .and. (zsum + dz(c,j) >= lsum ))  then
                call linear_interp(lsum, temp, zsum, zsum + dz(c,j), 0._r8, 1._r8)
                vic_clm_fract(c,i,j) = 1._r8 - temp
                if(lsum + deltal(i) < zsum + dz(c,j)) then
                   call linear_interp(lsum + deltal(i), temp, zsum, zsum + dz(c,j), 1._r8, 0._r8)
                   vic_clm_fract(c,i,j) = vic_clm_fract(c,i,j) - temp
                end if
             else if( (zsum < lsum + deltal(i)) .and. (zsum + dz(c,j) >= lsum + deltal(i)) ) then
                call linear_interp(lsum + deltal(i), temp, zsum, zsum + dz(c,j), 0._r8, 1._r8)
                vic_clm_fract(c,i,j) = temp
                if(zsum<=lsum) then
                   call linear_interp(lsum, temp, zsum, zsum + dz(c,j), 0._r8, 1._r8)
                   vic_clm_fract(c,i,j) = vic_clm_fract(c,i,j) - temp
                end if
             else if( (zsum >= lsum .and. zsum + dz(c,j) <= lsum + deltal(i)) )  then
                vic_clm_fract(c,i,j) = 1._r8
             else
                vic_clm_fract(c,i,j) = 0._r8
             end if
             zsum = zsum + dz(c,j)
             sum_frac(i) = sum_frac(i) + vic_clm_fract(c,i,j)
          end do                           ! end CLM layer calculation
          lsum = lsum + deltal(i)
       end do                             ! end VIC layer calcultion 

     end associate 

   end subroutine initCLMVICMap

   !-------------------------------------------------------------------
   subroutine linear_interp(x,y, x0, x1, y0, y1)
     !
     ! !DESCRIPTION:
     ! Provides linear interpolation
     !
     ! !ARGUMENTS:
     real(r8), intent(in)  :: x, x0, y0, x1, y1
     real(r8), intent(out) :: y
     !-------------------------------------------------------------------

     y = y0 + (x - x0) * (y1 - y0) / (x1 - x0)

   end subroutine linear_interp

end module SoilHydrologyInitTimeConstMod
