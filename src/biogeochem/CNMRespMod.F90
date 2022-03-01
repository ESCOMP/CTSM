module CNMRespMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module holding maintenance respiration routines for coupled carbon
  ! nitrogen code.
  !
  ! !USES:
  use shr_kind_mod           , only : r8 => shr_kind_r8
  use shr_const_mod          , only : SHR_CONST_TKFRZ
  use clm_varpar             , only : nlevgrnd
  use clm_varcon             , only : spval
  use decompMod              , only : bounds_type
  use abortutils             , only : endrun
  use shr_log_mod            , only : errMsg => shr_log_errMsg
  use pftconMod              , only : npcropmin, pftcon
  use SoilStateType          , only : soilstate_type
  use CanopyStateType        , only : canopystate_type
  use TemperatureType        , only : temperature_type
  use PhotosynthesisMod      , only : photosyns_type
  use CNVegcarbonfluxType    , only : cnveg_carbonflux_type
  use CNVegnitrogenstateType , only : cnveg_nitrogenstate_type
  use CNSharedParamsMod      , only : CNParamsShareInst
  use CropReprPoolsMod           , only : nrepr
  use PatchType              , only : patch                
  !
  implicit none
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: CNMRespReadNML       ! Read in namelist (CALL FIRST!)
  public :: readParams           ! Read in parameters from file
  public :: CNMResp              ! Apply maintenance respiration

  type, private :: params_type
     real(r8) :: br      = spval ! base rate for maintenance respiration (gC/gN/s)
     real(r8) :: br_root = spval ! base rate for maintenance respiration for roots (gC/gN/s)
  end type params_type

  type(params_type), private :: params_inst

  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine CNMRespReadNML( NLFilename )
    !
    ! !DESCRIPTION:
    ! Read the namelist for CNMResp (MUST BE CALLED BEFORE readParams!!!)
    !
    ! !USES:
    use fileutils      , only : getavu, relavu, opnfil
    use shr_nl_mod     , only : shr_nl_find_group_name
    use spmdMod        , only : masterproc, mpicom
    use shr_mpi_mod    , only : shr_mpi_bcast
    use clm_varctl     , only : iulog
    !
    ! !ARGUMENTS:
    character(len=*), intent(in) :: NLFilename ! Namelist filename
    !
    ! !LOCAL VARIABLES:
    integer :: ierr                 ! error code
    integer :: unitn                ! unit for namelist file

    character(len=*), parameter :: subname = 'CNMRespReadNML'
    character(len=*), parameter :: nmlname = 'cnmresp_inparm'
    real(r8) :: br_root = spval ! base rate for maintenance respiration for roots (gC/gN/s)
    !-----------------------------------------------------------------------

    namelist /cnmresp_inparm/ br_root

    ! Initialize options to default values, in case they are not specified in
    ! the namelist

    if (masterproc) then
       unitn = getavu()
       write(iulog,*) 'Read in '//nmlname//'  namelist'
       call opnfil (NLFilename, unitn, 'F')
       call shr_nl_find_group_name(unitn, nmlname, status=ierr)
       if (ierr == 0) then
          read(unitn, nml=cnmresp_inparm, iostat=ierr)
          if (ierr /= 0) then
             call endrun(msg="ERROR reading "//nmlname//"namelist"//errmsg(sourcefile, __LINE__))
          end if
       else
          call endrun(msg="ERROR finding "//nmlname//"namelist"//errmsg(sourcefile, __LINE__))
       end if
       call relavu( unitn )
    end if

    call shr_mpi_bcast (br_root, mpicom)

    params_inst%br_root = br_root

    if (masterproc) then
       write(iulog,*) ' '
       write(iulog,*) nmlname//' settings:'
       write(iulog,nml=cnmresp_inparm)
       write(iulog,*) ' '
    end if

  end subroutine CNMRespReadNML
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  subroutine readParams ( ncid )
    !
    ! !DESCRIPTION:
    ! Read parameters (call AFTER CNMRespReadNML!)
    !
    ! !USES:
    use ncdio_pio , only : file_desc_t,ncd_io
    !
    ! !ARGUMENTS:
    implicit none
    type(file_desc_t),intent(inout) :: ncid   ! pio netCDF file id
    !
    ! !LOCAL VARIABLES:
    character(len=32)  :: subname = 'CNMRespParamsType'
    character(len=100) :: errCode = '-Error reading in parameters file:'
    logical            :: readv ! has variable been read in or not
    real(r8)           :: tempr ! temporary to read in constant
    character(len=100) :: tString ! temp. var for reading
    !-----------------------------------------------------------------------

    tString='br_mr'
    call ncd_io(varname=trim(tString),data=tempr, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%br=tempr

    if ( params_inst%br_root == spval ) then
       params_inst%br_root = params_inst%br
    end if

  end subroutine readParams

  !-----------------------------------------------------------------------
  ! FIX(SPM,032414) this shouldn't even be called with fates on.
  !
  subroutine CNMResp(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, &
       canopystate_inst, soilstate_inst, temperature_inst, photosyns_inst, &
       cnveg_carbonflux_inst, cnveg_nitrogenstate_inst)
    !
    ! !DESCRIPTION:
    !
    ! !ARGUMENTS:
    use clm_varcon  , only : tfrz
    
    type(bounds_type)              , intent(in)    :: bounds          
    integer                        , intent(in)    :: num_soilc       ! number of soil points in column filter
    integer                        , intent(in)    :: filter_soilc(:) ! column filter for soil points
    integer                        , intent(in)    :: num_soilp       ! number of soil points in patch filter
    integer                        , intent(in)    :: filter_soilp(:) ! patch filter for soil points
    type(canopystate_type)         , intent(in)    :: canopystate_inst
    type(soilstate_type)           , intent(in)    :: soilstate_inst
    type(temperature_type)         , intent(in)    :: temperature_inst
    type(photosyns_type)           , intent(in)    :: photosyns_inst
    type(cnveg_carbonflux_type)    , intent(inout) :: cnveg_carbonflux_inst
    type(cnveg_nitrogenstate_type) , intent(in)    :: cnveg_nitrogenstate_inst
    !
    ! !LOCAL VARIABLES:
    integer :: c,p,j,k ! indices
    integer :: fp      ! soil filter patch index
    integer :: fc      ! soil filter column index
    real(r8):: br      ! base rate (gC/gN/s)
    real(r8):: br_root ! root base rate (gC/gN/s)
    real(r8):: q10     ! temperature dependence

    real(r8):: tc      ! temperature correction, 2m air temp (unitless)
    real(r8):: tcsoi(bounds%begc:bounds%endc,nlevgrnd) ! temperature correction by soil layer (unitless)
    !-----------------------------------------------------------------------

    associate(                                                            &    
         ivt            =>    patch%itype                                 , & ! Input:  [integer  (:)   ]  patch vegetation type                                

         woody          =>    pftcon%woody                              , & ! Input:  binary flag for woody lifeform (1=woody, 0=not woody)

         frac_veg_nosno =>    canopystate_inst%frac_veg_nosno_patch     , & ! Input:  [integer  (:)   ]  fraction of vegetation not covered by snow (0 OR 1) [-]
         laisun         =>    canopystate_inst%laisun_patch             , & ! Input:  [real(r8) (:)   ]  sunlit projected leaf area index                  
         laisha         =>    canopystate_inst%laisha_patch             , & ! Input:  [real(r8) (:)   ]  shaded projected leaf area index                  

         crootfr        =>    soilstate_inst%crootfr_patch              , & ! Input:  [real(r8) (:,:) ]  fraction of roots for carbon in each soil layer  (nlevgrnd)

         t_soisno       =>    temperature_inst%t_soisno_col             , & ! Input:  [real(r8) (:,:) ]  soil temperature (Kelvin)  (-nlevsno+1:nlevgrnd)
         t_ref2m        =>    temperature_inst%t_ref2m_patch            , & ! Input:  [real(r8) (:)   ]  2 m height surface air temperature (Kelvin)       

         t10            => temperature_inst%t_a10_patch                 , & ! Input:  [real(r8) (:)   ]  10-day running mean of the 2 m temperature (K)
 
         lmrsun         =>    photosyns_inst%lmrsun_patch               , & ! Input:  [real(r8) (:)   ]  sunlit leaf maintenance respiration rate (umol CO2/m**2/s)
         lmrsha         =>    photosyns_inst%lmrsha_patch               , & ! Input:  [real(r8) (:)   ]  shaded leaf maintenance respiration rate (umol CO2/m**2/s)
         rootstem_acc   =>    photosyns_inst%rootstem_acc               , & ! Input:  [logical        ]  root and stem acclimation switch

         frootn         =>    cnveg_nitrogenstate_inst%frootn_patch     , & ! Input:  [real(r8) (:)   ]  (gN/m2) fine root N                               
         livestemn      =>    cnveg_nitrogenstate_inst%livestemn_patch  , & ! Input:  [real(r8) (:)   ]  (gN/m2) live stem N                               
         livecrootn     =>    cnveg_nitrogenstate_inst%livecrootn_patch , & ! Input:  [real(r8) (:)   ]  (gN/m2) live coarse root N                        
         reproductiven         =>    cnveg_nitrogenstate_inst%reproductiven_patch     , & ! Input:  [real(r8) (:,:)   ]  (kgN/m2) grain N

         leaf_mr        =>    cnveg_carbonflux_inst%leaf_mr_patch       , & ! Output: [real(r8) (:)   ]                                                    
         froot_mr       =>    cnveg_carbonflux_inst%froot_mr_patch      , & ! Output: [real(r8) (:)   ]                                                    
         livestem_mr    =>    cnveg_carbonflux_inst%livestem_mr_patch   , & ! Output: [real(r8) (:)   ]                                                    
         livecroot_mr   =>    cnveg_carbonflux_inst%livecroot_mr_patch  , & ! Output: [real(r8) (:)   ]                                                    
         reproductive_mr       =>    cnveg_carbonflux_inst%reproductive_mr_patch        & ! Output: [real(r8) (:,:)   ]

         )

      ! base rate for maintenance respiration is from:
      ! M. Ryan, 1991. Effects of climate change on plant respiration.
      ! Ecological Applications, 1(2), 157-167.
      ! Original expression is br = 0.0106 molC/(molN h)
      ! Conversion by molecular weights of C and N gives 2.525e-6 gC/(gN s)
      ! set constants
      br      = params_inst%br
      br_root = params_inst%br_root
     
      ! Peter Thornton: 3/13/09 
      ! Q10 was originally set to 2.0, an arbitrary choice, but reduced to 1.5 as part of the tuning
      ! to improve seasonal cycle of atmospheric CO2 concentration in global
      ! simulatoins
      Q10 = CNParamsShareInst%Q10

      ! column loop to calculate temperature factors in each soil layer
      do j=1,nlevgrnd
         do fc = 1, num_soilc
            c = filter_soilc(fc)

            ! calculate temperature corrections for each soil layer, for use in
            ! estimating fine root maintenance respiration with depth
            tcsoi(c,j) = Q10**((t_soisno(c,j)-SHR_CONST_TKFRZ - 20.0_r8)/10.0_r8)
         end do
      end do

      ! patch loop for leaves and live wood
      do fp = 1, num_soilp
         p = filter_soilp(fp)

         ! calculate maintenance respiration fluxes in
         ! gC/m2/s for each of the live plant tissues.
         ! Leaf and live wood MR

         tc = Q10**((t_ref2m(p)-SHR_CONST_TKFRZ - 20.0_r8)/10.0_r8)
         
         !RF: acclimation of root and stem respiration fluxes
         ! n.b. we do not yet know if this is defensible scientifically (awaiting data analysis)
         ! turning this on will increase R and decrease productivity in boreal forests, A LOT. :)

         if(rootstem_acc)then 
           br      = br      * 10._r8**(-0.00794_r8*((t10(p)-tfrz)-25._r8))
           br_root = br_root * 10._r8**(-0.00794_r8*((t10(p)-tfrz)-25._r8))
         end if

         if (frac_veg_nosno(p) == 1) then

            leaf_mr(p) = lmrsun(p) * laisun(p) * 12.011e-6_r8 + &
                         lmrsha(p) * laisha(p) * 12.011e-6_r8

         else !nosno

            leaf_mr(p) = 0._r8

         end if

         if (woody(ivt(p)) == 1) then
            livestem_mr(p) = livestemn(p)*br*tc
            livecroot_mr(p) = livecrootn(p)*br_root*tc
         else if (ivt(p) >= npcropmin) then
            livestem_mr(p) = livestemn(p)*br*tc
            do k = 1, nrepr
               reproductive_mr(p,k) = reproductiven(p,k)*br*tc
            end do
         end if
      end do

      ! soil and patch loop for fine root

      do j = 1,nlevgrnd
         do fp = 1,num_soilp
            p = filter_soilp(fp)
            c = patch%column(p)

            ! Fine root MR
            ! crootfr(j) sums to 1.0 over all soil layers, and
            ! describes the fraction of root mass for carbon that is in each
            ! layer.  This is used with the layer temperature correction
            ! to estimate the total fine root maintenance respiration as a
            ! function of temperature and N content.
            if(rootstem_acc)then
               br_root = br_root * 10._r8**(-0.00794_r8*((t10(p)-tfrz)-25._r8))
            end if
            froot_mr(p) = froot_mr(p) + frootn(p)*br_root*tcsoi(c,j)*crootfr(p,j)
            
         end do
      end do

    end associate

  end subroutine CNMResp

end module CNMRespMod
