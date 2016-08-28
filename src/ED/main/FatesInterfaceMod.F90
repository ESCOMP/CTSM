module FatesInterfaceMod

   ! ------------------------------------------------------------------------------------
   ! This is the FATES public API
   ! A host land model has defined and allocated a structure "fates" as
   ! defined by fates_interface_type
   !
   ! It is also likely/possible that this type is defined as a vector
   ! which is allocated by thread
   ! ------------------------------------------------------------------------------------

   ! ------------------------------------------------------------------------------------
   ! Used CLM Modules
   ! INTERF-TODO:  NO CLM MODULES SHOULD BE ACCESSIBLE BY THE FATES
   ! PUBLIC API!!!!
   ! ------------------------------------------------------------------------------------

   use EDtypesMod            , only : ed_site_type,      &
                                      numPatchesPerCol,  &
                                      cp_nclmax,         &
                                      cp_numSWb,         &
                                      cp_numlevgrnd,     &
                                      cp_maxSWb,         &
                                      cp_numlevdecomp,   &
                                      cp_numlevdecomp_full 

   use shr_kind_mod          , only : r8 => shr_kind_r8  ! INTERF-TODO: REMOVE THIS

   

   ! ------------------------------------------------------------------------------------
   ! Notes on types
   ! For floating point arrays, it is sometimes the convention to define the arrays as
   ! POINTER instead of ALLOCATABLE.  This usually achieves the same result with subtle
   ! differences.  POINTER arrays can point to scalar values, discontinuous array slices
   ! or alias other variables, ALLOCATABLES cannnot.  According to S. Lionel 
   ! (Intel-Forum Post), ALLOCATABLES are better perfomance wise as long as they point 
   ! to contiguous memory spaces and do not alias other variables, the case here.
   ! Naming conventions:   _gl  means ground layer dimensions
   !                       _pa  means patch dimensions
   !                       _rb  means radiation band
   ! ------------------------------------------------------------------------------------
   
   type, public :: bc_in_type

      ! The actual number of FATES' ED patches
      integer :: npatches

      ! Radiation variables for calculating sun/shade fractions
      ! ---------------------------------------------------------------------------------

      ! Downwelling direct beam radiation (patch,radiation-band) [W/m2]
      real(r8), allocatable :: solad_parb(:,:)  

      ! Downwelling diffuse (I-ndirect) radiation (patch,radiation-band) [W/m2]
      real(r8), allocatable :: solai_parb(:,:)

      ! Hydrology variables for BTRAN
      ! ---------------------------------------------------------------------------------

      ! Soil suction potential of layers in each site, negative, [mm]
      real(r8), allocatable :: smp_gl(:)

      ! Effective porosity = porosity - vol_ic, of layers in each site [-]
      real(r8), allocatable :: eff_porosity_gl(:)

      ! volumetric soil water at saturation (porosity)
      real(r8), allocatable :: watsat_gl(:)

      ! Temperature of ground layers [K]
      real(r8), allocatable :: tempk_gl(:)

      ! Liquid volume in ground layer
      real(r8), allocatable :: h2o_liqvol_gl(:)

      ! Site level filter for uptake response functions
      logical               :: filter_btran

      ! Photosynthesis variables
      ! ---------------------------------------------------------------------------------

      ! Patch level filter flag for photosynthesis calculations
      ! has a short memory, flags:
      ! 1 = patch has not been called
      ! 2 = patch is currently marked for photosynthesis
      ! 3 = patch has been called for photosynthesis at least once
      integer, allocatable  :: filter_photo_pa(:)

      ! atmospheric pressure (Pa)
      real(r8)              :: forc_pbot             

      ! daylength scaling factor (0-1)
      real(r8), allocatable :: dayl_factor_pa(:)
      
      ! saturation vapor pressure at t_veg (Pa)
      real(r8), allocatable :: esat_tv_pa(:)

      ! vapor pressure of canopy air (Pa)
      real(r8), allocatable :: eair_pa(:)

      ! Atmospheric O2 partial pressure (Pa)
      real(r8), allocatable :: oair_pa(:)

      ! Atmospheric CO2 partial pressure (Pa)
      real(r8), allocatable :: cair_pa(:)

      ! boundary layer resistance (s/m)
      real(r8), allocatable :: rb_pa(:)

      ! vegetation temperature (Kelvin)
      real(r8), allocatable :: t_veg_pa(:)
             
      ! air temperature at agcm reference height (kelvin)
      real(r8), allocatable :: tgcm_pa(:)

      ! soil temperature (Kelvin)
      real(r8), allocatable :: t_soisno_gl(:)

      ! Canopy Radiation Boundaries
      ! ---------------------------------------------------------------------------------
      
      ! Filter for vegetation patches with a positive zenith angle (daylight)
      logical, allocatable :: filter_vegzen_pa(:)

      ! Cosine of the zenith angle (0-1), by patch
      ! Note RGK: It does not seem like the code would currently generate
      !           different zenith angles for different patches (nor should it)
      !           I am leaving it at this scale for simplicity.  Patches should
      !           have no spacially variable information
      real(r8), allocatable :: coszen_pa(:)
      
      ! Abledo of the ground for direct radiation, by site broadband (0-1)
      real(r8), allocatable :: albgr_dir_rb(:)

      ! Albedo of the ground for diffuse radiation, by site broadband (0-1)
      real(r8), allocatable :: albgr_dif_rb(:)
      
      ! LitterFlux Boundaries
      ! the index of the deepest model soil level where roots may be
      ! due to permafrost or bedrock constraints
      integer  :: max_rooting_depth_index_col


   end type bc_in_type


   type, public :: bc_out_type

      ! Sunlit fraction of the canopy for this patch [0-1]
      real(r8),allocatable :: fsun_pa(:)

      ! Sunlit canopy LAI
      real(r8),allocatable :: laisun_pa(:)
      
      ! Shaded canopy LAI
      real(r8),allocatable :: laisha_pa(:)
      
      ! Logical stating whether a ground layer can have water uptake by plants
      ! The only condition right now is that liquid water exists
      ! The name (suction) is used to indicate that soil suction should be calculated
      logical, allocatable :: active_suction_gl(:)

      ! Effective fraction of roots in each soil layer 
      real(r8), allocatable :: rootr_pagl(:,:)

      ! Integrated (vertically) transpiration wetness factor (0 to 1) 
      ! (diagnostic, should not be used by HLM)
      real(r8), allocatable :: btran_pa(:)

      ! Sunlit canopy resistance [s/m]
      real(r8), allocatable :: rssun_pa(:)

      ! Shaded canopy resistance [s/m]
      real(r8), allocatable :: rssha_pa(:)

      ! Canopy conductance [mmol m-2 s-1]
      real(r8), allocatable :: gccanopy_pa(:)

      ! patch sunlit leaf photosynthesis (umol CO2 /m**2/ s)
      real(r8), allocatable :: psncanopy_pa(:)

      ! patch sunlit leaf maintenance respiration rate (umol CO2/m**2/s) 
      real(r8), allocatable :: lmrcanopy_pa(:)

      ! Canopy Radiation Boundaries
      ! ---------------------------------------------------------------------------------
      
      ! Surface albedo (direct) (HLMs use this for atm coupling and balance checks)
      real(r8), allocatable :: albd_parb(:,:)
      
      ! Surface albedo (diffuse) (HLMs use this for atm coupling and balance checks)
      real(r8), allocatable :: albi_parb(:,:)                 
      
      ! Flux absorbed by canopy per unit direct flux (HLMs use this for balance checks)
      real(r8), allocatable :: fabd_parb(:,:) 
      
      ! Flux absorbed by canopy per unit diffuse flux (HLMs use this for balance checks)
      real(r8), allocatable :: fabi_parb(:,:)

      ! Down direct flux below canopy per unit direct flx (HLMs use this for balance checks)
      real(r8), allocatable :: ftdd_parb(:,:)

      ! Down diffuse flux below canopy per unit direct flx (HLMs use this for balance checks)
      real(r8), allocatable :: ftid_parb(:,:)
      
      ! Down diffuse flux below canopy per unit diffuse flx (HLMs use this for balance checks)
      real(r8), allocatable :: ftii_parb(:,:)


      ! litterfall fluxes of C from FATES patches to BGC columns

      ! total labile    litter coming from ED. gC/m3/s
      real(r8), allocatable :: FATES_c_to_litr_lab_c_col(:)      

      !total cellulose litter coming from ED. gC/m3/s
      real(r8), allocatable :: FATES_c_to_litr_cel_c_col(:)      
      
      !total lignin    litter coming from ED. gC/m3/s
      real(r8), allocatable :: FATES_c_to_litr_lig_c_col(:)      

      
   end type bc_out_type


   type, public :: fates_interface_type
      
      ! This is the root of the ED/FATES hierarchy of instantaneous state variables
      ! ie the root of the linked lists. Each path list is currently associated with a 
      ! grid-cell, this is intended to be migrated to columns 

      integer                         :: nsites

      type(ed_site_type), pointer :: sites(:)

      ! These are boundary conditions that the FATES models are required to be filled.  
      ! These values are filled by the driver or HLM.  Once filled, these have an 
      ! intent(in) status.  Each site has a derived type structure, which may include 
      ! a scalar for site level data, a patch vector, potentially cohort vectors (but 
      ! not yet atm) and other dimensions such as soil-depth or pft.  These vectors 
      ! are initialized by maximums, and the allocations are static in time to avoid
      ! having to allocate/de-allocate memory

      type(bc_in_type), allocatable   :: bc_in(:)

      ! These are the boundary conditions that the FATES model returns to its HLM or 
      ! driver. It has the same allocation strategy and similar vector types.
      
      type(bc_out_type), allocatable  :: bc_out(:)

   contains
      
      procedure, public :: zero_bcs

   end type fates_interface_type

   public :: FatesInterfaceInit
   public :: set_fates_ctrlparms


contains

   ! ====================================================================================
  subroutine FatesInterfaceInit(log_unit, global_verbose)

    use FatesGlobals, only : FatesGlobalsInit

    implicit none

    integer, intent(in) :: log_unit
    logical, intent(in) :: global_verbose

    call FatesGlobalsInit(log_unit, global_verbose)

  end subroutine FatesInterfaceInit

   ! ====================================================================================

   ! INTERF-TODO: THIS IS A PLACE-HOLDER ROUTINE, NOT CALLED YET...
   subroutine fates_clean(this)
      
      implicit none
      
      ! Input Arguments
      class(fates_interface_type), intent(inout) :: this
      
      ! Incrementally walk through linked list and deallocate
      
      
      
      ! Deallocate the site list
      deallocate (this%sites)
      
      return
   end subroutine fates_clean


   ! ====================================================================================
   

   subroutine allocate_bcin(bc_in)
      
      ! ---------------------------------------------------------------------------------
      ! Allocate and Initialze the FATES boundary condition vectors
      ! ---------------------------------------------------------------------------------
      
      implicit none
      type(bc_in_type), intent(inout) :: bc_in
      
      ! Allocate input boundaries
      
      ! Radiation
      allocate(bc_in%solad_parb(numPatchesPerCol,cp_numSWb))
      allocate(bc_in%solai_parb(numPatchesPerCol,cp_numSWb))
      
      ! Hydrology
      allocate(bc_in%smp_gl(cp_numlevgrnd))
      allocate(bc_in%eff_porosity_gl(cp_numlevgrnd))
      allocate(bc_in%watsat_gl(cp_numlevgrnd))
      allocate(bc_in%tempk_gl(cp_numlevgrnd))
      allocate(bc_in%h2o_liqvol_gl(cp_numlevgrnd))

      ! Photosynthesis
      allocate(bc_in%filter_photo_pa(numPatchesPerCol))
      allocate(bc_in%dayl_factor_pa(numPatchesPerCol))
      allocate(bc_in%esat_tv_pa(numPatchesPerCol))
      allocate(bc_in%eair_pa(numPatchesPerCol))
      allocate(bc_in%oair_pa(numPatchesPerCol))
      allocate(bc_in%cair_pa(numPatchesPerCol))
      allocate(bc_in%rb_pa(numPatchesPerCol))
      allocate(bc_in%t_veg_pa(numPatchesPerCol))
      allocate(bc_in%tgcm_pa(numPatchesPerCol))
      allocate(bc_in%t_soisno_gl(cp_numlevgrnd))

      ! Canopy Radiation
      allocate(bc_in%filter_vegzen_pa(numPatchesPerCol))
      allocate(bc_in%coszen_pa(numPatchesPerCol))
      allocate(bc_in%albgr_dir_rb(cp_numSWb))
      allocate(bc_in%albgr_dif_rb(cp_numSWb))

      return
   end subroutine allocate_bcin
   
   subroutine allocate_bcout(bc_out)

      ! ---------------------------------------------------------------------------------
      ! Allocate and Initialze the FATES boundary condition vectors
      ! ---------------------------------------------------------------------------------
      
      implicit none
      type(bc_out_type), intent(inout) :: bc_out
      
      
      ! Radiation
      allocate(bc_out%fsun_pa(numPatchesPerCol))
      allocate(bc_out%laisun_pa(numPatchesPerCol))
      allocate(bc_out%laisha_pa(numPatchesPerCol))
      
      ! Hydrology
      allocate(bc_out%active_suction_gl(cp_numlevgrnd))
      allocate(bc_out%rootr_pagl(numPatchesPerCol,cp_numlevgrnd))
      allocate(bc_out%btran_pa(numPatchesPerCol))
      
      ! Photosynthesis
      allocate(bc_out%rssun_pa(numPatchesPerCol))
      allocate(bc_out%rssha_pa(numPatchesPerCol))
      allocate(bc_out%gccanopy_pa(numPatchesPerCol))
      allocate(bc_out%lmrcanopy_pa(numPatchesPerCol))
      allocate(bc_out%psncanopy_pa(numPatchesPerCol))
      
      ! Canopy Radiation
      allocate(bc_out%albd_parb(numPatchesPerCol,cp_numSWb))
      allocate(bc_out%albi_parb(numPatchesPerCol,cp_numSWb))
      allocate(bc_out%fabd_parb(numPatchesPerCol,cp_numSWb))
      allocate(bc_out%fabi_parb(numPatchesPerCol,cp_numSWb))
      allocate(bc_out%ftdd_parb(numPatchesPerCol,cp_numSWb))
      allocate(bc_out%ftid_parb(numPatchesPerCol,cp_numSWb))
      allocate(bc_out%ftii_parb(numPatchesPerCol,cp_numSWb))

      ! biogeochemistry
      allocate(bc_out%FATES_c_to_litr_lab_c_col(cp_numlevdecomp_full))        
      allocate(bc_out%FATES_c_to_litr_cel_c_col(cp_numlevdecomp_full))
      allocate(bc_out%FATES_c_to_litr_lig_c_col(cp_numlevdecomp_full))

      return
   end subroutine allocate_bcout

   ! ====================================================================================

   subroutine zero_bcs(this,s)

      implicit none
      class(fates_interface_type), intent(inout) :: this
      integer, intent(in) :: s

      ! Input boundaries

      this%bc_in(s)%solad_parb(:,:)     = 0.0_r8
      this%bc_in(s)%solai_parb(:,:)     = 0.0_r8
      this%bc_in(s)%smp_gl(:)           = 0.0_r8
      this%bc_in(s)%eff_porosity_gl(:)  = 0.0_r8
      this%bc_in(s)%watsat_gl(:)        = 0.0_r8
      this%bc_in(s)%tempk_gl(:)         = 0.0_r8
      this%bc_in(s)%h2o_liqvol_gl(:)    = 0.0_r8
      this%bc_in(s)%filter_vegzen_pa(:) = .false.
      this%bc_in(s)%coszen_pa(:)        = 0.0_r8
      this%bc_in(s)%albgr_dir_rb(:)     = 0.0_r8
      this%bc_in(s)%albgr_dif_rb(:)     = 0.0_r8
      this%bc_in(s)%max_rooting_depth_index_col = 0

      
      ! Output boundaries
      this%bc_out(s)%active_suction_gl(:) = .false.
      this%bc_out(s)%fsun_pa(:)      = 0.0_r8
      this%bc_out(s)%laisun_pa(:)    = 0.0_r8
      this%bc_out(s)%laisha_pa(:)    = 0.0_r8
      this%bc_out(s)%rootr_pagl(:,:) = 0.0_r8
      this%bc_out(s)%btran_pa(:)     = 0.0_r8

      this%bc_out(s)%FATES_c_to_litr_lab_c_col(:) = 0.0_r8
      this%bc_out(s)%FATES_c_to_litr_cel_c_col(:) = 0.0_r8
      this%bc_out(s)%FATES_c_to_litr_lig_c_col(:) = 0.0_r8

      this%bc_out(s)%rssun_pa(:)     = 0.0_r8
      this%bc_out(s)%rssha_pa(:)     = 0.0_r8
      this%bc_out(s)%gccanopy_pa(:)  = 0.0_r8
      this%bc_out(s)%psncanopy_pa(:) = 0.0_r8
      this%bc_out(s)%lmrcanopy_pa(:) = 0.0_r8

      this%bc_out(s)%albd_parb(:,:) = 0.0_r8
      this%bc_out(s)%albi_parb(:,:) = 0.0_r8
      this%bc_out(s)%fabd_parb(:,:) = 0.0_r8
      this%bc_out(s)%fabi_parb(:,:) = 0.0_r8
      this%bc_out(s)%ftdd_parb(:,:) = 0.0_r8
      this%bc_out(s)%ftid_parb(:,:) = 0.0_r8
      this%bc_out(s)%ftii_parb(:,:) = 0.0_r8
      
      return
    end subroutine zero_bcs
   
   ! ==================================================================================== 

   subroutine set_fates_ctrlparms(tag,dimval)
      
      ! ---------------------------------------------------------------------------------
      ! INTERF-TODO:  NEED ALLOWANCES FOR REAL AND CHARACTER ARGS..
      ! Certain model control parameters and dimensions used by FATES are dictated by 
      ! the the driver or the host mode. To see which parameters should be filled here
      ! please also look at the ctrl_parms_type in FATESTYpeMod, in the section listing
      ! components dictated by the host model.
      !
      ! Some important points:
      ! 1. Calls to this function are likely from the clm_fates module in the HLM.
      ! 2. The calls should be preceeded by a flush function.
      ! 3. All values in ctrl_parm (FATESTypesMod.F90) that are classified as 
      !    'dictated by the HLM' must be listed in this subroutine
      ! 4. Should look like this:
      ! 
      ! call set_fates_ctrlparms('flush_to_unset')
      ! call set_fates_ctrlparms('num_sw_bbands',numrad)  ! or other variable
      ! ...
      ! call set_fates_ctrlparms('num_lev_ground',nlevgrnd)   ! or other variable
      ! call set_fates_ctrlparms('check_allset') 
      !
      ! RGK-2016
      ! ---------------------------------------------------------------------------------
     use FatesGlobals, only : fates_log, fates_global_verbose

      ! Arguments
      integer, optional, intent(in)  :: dimval
      character(len=*),intent(in)    :: tag
      
      ! local variables
      logical              :: all_set
      integer,  parameter  :: unset_int = -999
      real(r8), parameter  :: unset_double = -999.9
      
      select case (trim(tag))
      case('flush_to_unset')
         if (fates_global_verbose()) then
            write(fates_log(), *) 'Flushing FATES control parameters prior to transfer from host'
         end if
         cp_numSwb     = unset_int
         cp_numlevgrnd = unset_int
         cp_numlevdecomp_full = unset_int
         cp_numlevdecomp      = unset_int


      case('check_allset')
         
         if(cp_numSWb .eq. unset_int) then
            if (fates_global_verbose()) then
               write(fates_log(), *) 'FATES dimension/parameter unset: num_sw_rad_bbands'
            end if
            ! INTERF-TODO: FATES NEEDS INTERNAL end_run
            ! end_run('MESSAGE')
         end if

         if(cp_numSWb > cp_maxSWb) then
            if (fates_global_verbose()) then
               write(fates_log(), *) 'FATES sets a maximum number of shortwave bands'
               write(fates_log(), *) 'for some scratch-space, cp_maxSWb'
               write(fates_log(), *) 'it defaults to 2, but can be increased as needed'
               write(fates_log(), *) 'your driver or host model is intending to drive'
               write(fates_log(), *) 'FATES with:',cp_numSWb,' bands.'
               write(fates_log(), *) 'please increase cp_maxSWb in EDTypes to match'
               write(fates_log(), *) 'or exceed this value'
            end if
            ! end_run('MESSAGE')
         end if

         if(cp_numlevgrnd .eq. unset_int) then
            if (fates_global_verbose()) then
               write(fates_log(), *) 'FATES dimension/parameter unset: numlevground'
            end if
            ! INTERF-TODO: FATES NEEDS INTERNAL end_run
            ! end_run('MESSAGE')
         end if

         if(cp_numlevdecomp_full .eq. unset_int) then
            if (fates_global_verbose()) then
               write(fates_log(), *) 'FATES dimension/parameter unset: numlevdecomp_full'
            end if
            ! INTERF-TODO: FATES NEEDS INTERNAL end_run
            ! end_run('MESSAGE')
         end if

         if(cp_numlevdecomp .eq. unset_int) then
            if (fates_global_verbose()) then
               write(fates_log(), *) 'FATES dimension/parameter unset: numlevdecomp'
            end if
            ! INTERF-TODO: FATES NEEDS INTERNAL end_run
            ! end_run('MESSAGE')
         end if

         if (fates_global_verbose()) then
            write(fates_log(), *) 'Checked. All control parameters sent to FATES.'
         end if
         
      case default

         if(present(dimval))then
            select case (trim(tag))
               
            case('num_sw_bbands')
               
               cp_numSwb = dimval
               if (fates_global_verbose()) then
                  write(fates_log(), *) 'Transfering num_sw_bbands = ',dimval,' to FATES'
               end if
               
            case('num_lev_ground')
               
               cp_numlevgrnd = dimval
               if (fates_global_verbose()) then

                  write(fates_log(), *) 'Transfering num_lev_ground = ',dimval,' to FATES'
               end if
               
            case('num_levdecomp_full')
               
               cp_numlevdecomp_full = dimval
               if (fates_global_verbose()) then
                  write(fates_log(), *) 'Transfering num_levdecomp_full = ',dimval,' to FATES'
               end if
            
            case('num_levdecomp')
               
               cp_numlevdecomp = dimval
               if (fates_global_verbose()) then
                  write(fates_log(), *) 'Transfering num_levdecomp = ',dimval,' to FATES'
               end if

            case default
               if (fates_global_verbose()) then
                  write(fates_log(), *) 'tag not recognized:',trim(tag)
               end if
               ! end_run
            end select
         else
            if (fates_global_verbose()) then
               write(fates_log(), *) 'no value was provided for the tag'
            end if
         end if
         
      end select
      
      
      return
   end subroutine set_fates_ctrlparms



end module FatesInterfaceMod
