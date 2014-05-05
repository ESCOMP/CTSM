module EDhistFldsMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module containing initialization of ED history fields and files
  ! This is the module that the user must modify in order to add new
  ! history fields or modify defaults associated with existing history
  ! fields.
  !
  ! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8

  implicit none
  save
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public EDhist_initFlds ! Builds ED variables into master field list of 
                         ! all possible history file fields
  !------------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine EDhist_initFlds()
    !
    ! !DESCRIPTION:
    ! Build master field list of all possible ED fields in a history file.
    ! Each field has associated with it a ``long\_name'' netcdf attribute that
    ! describes what the field is, and a ``units'' attribute. A subroutine is
    ! called to add each field to the masterlist.
    !
    ! !USES:
    use EDClmtype,      only : EDpcf
    use clmtype,        only : pps
    use histFileMod,    only : hist_addfld1d, hist_addfld2d
    !
    ! !ARGUMENTS:
    implicit none
    !
    ! !LOCAL VARIABLES:
    !-----------------------------------------------------------------------

          call hist_addfld1d (fname='TRIMMING', units='none',  &
                avgflag='A', long_name='Degree to which canopy expansion is limited by leaf economics', &
                ptr_pft=EDpcf%trimming, set_lake=0._r8, set_urb=0._r8)  
       
           call hist_addfld1d (fname='AREA_PLANT', units='m2',  &
                avgflag='A', long_name='area occupied by all plants', &
                ptr_pft=EDpcf%area_plant, set_lake=0._r8, set_urb=0._r8)
       
           call hist_addfld1d (fname='AREA_TREES', units='m2',  &
                avgflag='A', long_name='area occupied by woody plants', &
                ptr_pft=EDpcf%area_trees, set_lake=0._r8, set_urb=0._r8)
       
           call hist_addfld1d (fname='CANOPY_SPREAD', units='none',  &   
                avgflag='A', long_name='Scaling factor between tree basal area and canopy area', &
                ptr_pft=EDpcf%canopy_spread, set_lake=0._r8, set_urb=0._r8)   
           
           call hist_addfld1d (fname='GCCANOPY', units='none',  &
                avgflag='A', long_name='Canopy Conductance: mmol m-2 s-1', &
                ptr_pft=EDpcf%GCcanopy, set_lake=0._r8, set_urb=0._r8)  
            
           call hist_addfld2d (fname='PFTbiomass',  units='kgC/m2', type2d='levgrnd', &
                avgflag='A', long_name='total PFT level biomass', &
                ptr_pft=EDpcf%PFTbiomass)
           
           call hist_addfld2d (fname='PFTleafbiomass',  units='kgC/m2', type2d='levgrnd', &
                avgflag='A', long_name='total PFT level biomass', &
                ptr_pft=EDpcf%PFTleafbiomass)
               
           call hist_addfld2d (fname='PFTstorebiomass',  units='kgC/m2', type2d='levgrnd', &
                avgflag='A', long_name='total PFT level biomass', &
                ptr_pft=EDpcf%PFTstorebiomass)
           
           call hist_addfld2d (fname='PFTnindivs',  units='kgC/m2', type2d='levgrnd', &
                avgflag='A', long_name='total PFT level biomass', &
                ptr_pft=EDpcf%PFTnindivs)
       
           call hist_addfld1d (fname='FIRE_NESTEROV_INDEX', units='none',  &
                avgflag='A', long_name='nesterov_fire_danger index', &
                ptr_pft=EDpcf%nesterov_fire_danger, set_lake=0._r8, set_urb=0._r8)
           
           call hist_addfld1d (fname='FIRE_ROS', units='m/min',  &
                avgflag='A', long_name='fire rate of spread m/min', &
                ptr_pft=EDpcf%spitfire_ROS, set_lake=0._r8, set_urb=0._r8)
 
           call hist_addfld1d (fname='EFFECT_WSPEED', units='none',  &
                avgflag='A', long_name='effective windspeed for fire spread', &
                ptr_pft=EDpcf%effect_wspeed, set_lake=0._r8, set_urb=0._r8)
           
           call hist_addfld1d (fname='FIRE_TFC_ROS', units='none',  &
                avgflag='A', long_name='total fuel consumed', &
                ptr_pft=EDpcf%TFC_ROS, set_lake=0._r8, set_urb=0._r8)
           
           call hist_addfld1d (fname='FIRE_INTENSITY', units='kJ/m/s',  &
                avgflag='A', long_name='spitfire fire intensity: kJ/m/s', &
                ptr_pft=EDpcf%fire_intensity, set_lake=0._r8, set_urb=0._r8)
              
           call hist_addfld1d (fname='FIRE_AREA', units='fraction',  &
                avgflag='A', long_name='spitfire fire area:m2', &
                ptr_pft=EDpcf%fire_area, set_lake=0._r8, set_urb=0._r8)
       
           call hist_addfld1d (fname='SCORCH_HEIGHT', units='m',  &
                avgflag='A', long_name='spitfire fire area:m2', &
                ptr_pft=EDpcf%scorch_height, set_lake=0._r8, set_urb=0._r8)
  
           call hist_addfld1d (fname='fire_fuel_mef', units='m',  &
                avgflag='A', long_name='spitfire fuel moisture', &
                ptr_pft=EDpcf%fire_fuel_mef, set_lake=0._r8, set_urb=0._r8)

           call hist_addfld1d (fname='fire_fuel_bulkd', units='m',  &
                avgflag='A', long_name='spitfire fuel bulk density', &
                ptr_pft=EDpcf%fire_fuel_bulkd, set_lake=0._r8, set_urb=0._r8)

           call hist_addfld1d (fname='fire_fuel_eff_moist', units='m',  &
                avgflag='A', long_name='spitfire fuel moisture', &
                ptr_pft=EDpcf%fire_fuel_eff_moist, set_lake=0._r8, set_urb=0._r8)
       
           call hist_addfld1d (fname='fire_fuel_sav', units='m',  &
                avgflag='A', long_name='spitfire fuel surface/volume ', &
                ptr_pft=EDpcf%fire_fuel_sav, set_lake=0._r8, set_urb=0._r8)

           call hist_addfld1d (fname='TFC_ROS', units='m',  &
                avgflag='A', long_name='spitfire fuel surface/volume ', &
                ptr_pft=EDpcf%TFC_ROS, set_lake=0._r8, set_urb=0._r8)

            call hist_addfld1d (fname='SUM_FUEL', units=' KgC m-2 y-1',  &
                avgflag='A', long_name='Litter flux in leaves', &
                ptr_pft=EDpcf%sum_fuel, set_lake=0._r8, set_urb=0._r8)
    

           call hist_addfld1d (fname='LITTER_IN', units=' KgC m-2 y-1',  &
                avgflag='A', long_name='Litter flux in leaves', &
                ptr_pft=EDpcf%litter_in, set_lake=0._r8, set_urb=0._r8)
    
           call hist_addfld1d (fname='LITTER_OUT', units=' KgC m-2 y-1',  &
                avgflag='A', long_name='Litter flux out leaves', &
                ptr_pft=EDpcf%litter_out, set_lake=0._r8, set_urb=0._r8)

           call hist_addfld1d (fname='SEED_BANK', units=' KgC m-2',  &
                avgflag='A', long_name='Total Seed Mass of all PFTs', &
                ptr_pft=EDpcf%seed_bank, set_lake=0._r8, set_urb=0._r8)

           call hist_addfld1d (fname='SEEDS_IN', units=' KgC m-2 y-1',  &
                avgflag='A', long_name='Seed Production Rate', &
                ptr_pft=EDpcf%seeds_in, set_lake=0._r8, set_urb=0._r8)
      
           call hist_addfld1d (fname='SEED_GERMINATION', units=' KgC m-2 y-1',  &
                avgflag='A', long_name='Seed mass converted into new cohorts', &
                ptr_pft=EDpcf%seed_germination, set_lake=0._r8, set_urb=0._r8)
      
           call hist_addfld1d (fname='SEED_DECAY', units=' KgC m-2 y-1',  &
                avgflag='A', long_name='Seed mass decay', &
                ptr_pft=EDpcf%seed_decay, set_lake=0._r8, set_urb=0._r8)              

           call hist_addfld1d (fname='RB', units=' s m-1',  &
                avgflag='A', long_name='leaf boundary resistance', &
                ptr_pft=EDpcf%rb, set_lake=0._r8, set_urb=0._r8)

           call hist_addfld1d (fname='EFPOT', units='',  &
                avgflag='A', long_name='potential evap', &
                ptr_pft=EDpcf%efpot, set_lake=0._r8, set_urb=0._r8)

           ! FIX(SPM, 041614) Should rscanopy be ed specific?  right now it is in clmtype
           call hist_addfld1d (fname='RSCANOPY', units=' s m-1',  &
                avgflag='A', long_name='canopy resistance', &
                ptr_pft=pps%rscanopy, set_lake=0._r8, set_urb=0._r8)

  end subroutine EDhist_initFlds

end module EDhistFldsMod
