module GroundwaterMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Calculates prognostic groundwater and Water Table Depth (WTD) considering pumping
  ! and groundwater lateral flow
  !
  ! Usage:
  !
  ! Design notes:
  !
  !
  !   Assumptions: The head across a given radius is constant (i.e., cone of depression is circular)
  !
  !
  ! HISTORY:
  !
  !   First Version: Farshid Felfelani 2019-03-01
  !                  Developed by Farshid Felfelani, first results are presented in the WRR paper:
  !                  Felfelani et al. (2021) "Representing Intercell Lateral Groundwater Flow and 
  !                                           Aquifer Pumping in the Community Land Model"
  !                  https://doi.org/10.1029/2020WR027531
  !
  !   Farshid Felfelani 2022-05-23
  !   There seems to be a water balance issue, to solve it:
  !   I added residual lateral flow (when >0) to the sub-surface runoff
  !
  ! !USES:
#include "shr_assert.h"
  use shr_kind_mod      , only : r8 => shr_kind_r8
  use decompMod         , only : bounds_type, get_proc_global
  use shr_log_mod       , only : errMsg => shr_log_errMsg
  use abortutils        , only : endrun
  use clm_varctl        , only : iulog, use_pumping
  use clm_varcon        , only : isecspday, degpsec, denh2o, spval, namec, rpi, aquifer_water_baseline, watmin, pondmx
  use clm_varpar        , only : nlevsoi, nlevgrnd
  use GridcellType      , only : grc
  use LandunitType      , only : lun    
  use ColumnType        , only : col
  use PatchType         , only : patch                
  use subgridAveMod     , only : p2c, c2g
  use filterColMod      , only : filter_col_type, col_filter_from_logical_array
  use SoilHydrologyType , only : soilhydrology_type  
  use SoilStateType     , only : soilstate_type
  use WaterfluxType     , only : waterflux_type
  use WaterstateType    , only : waterstate_type
  use IrrigationMod     , only : irrigation_type
  use spmdMod           , only : iam, masterproc  ! FFelfelani: to get processor number
  use decompMod         , only : get_proc_global, get_proc_bounds, get_clump_bounds,get_proc_clumps  ! FFelfelani: to get number of gridcells
  use SoilWaterRetentionCurveMod, only : soil_water_retention_curve_type

  ! !PUBLIC TYPES:
  implicit none
  private
  
   type, public :: groundwater_type
     private
    contains
   
     ! Public routines
     procedure, public :: UpdateGWFanLatPump
     procedure, public :: UpdateGWDefaultPump

     ! Private routines
     procedure, private :: TransmissivityFromFan
   end type groundwater_type 

  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  
contains

  ! ======================================================================== 
  ! Infrastructure routines (initialization, restart, etc.)
  ! ========================================================================
  
  !------------------------------------------------------------------------
  subroutine UpdateGWFanLatPump(this, bounds, num_hydrologyc, filter_hydrologyc, &
        soilhydrology_inst, soilstate_inst,waterstate_inst, irrigation_inst, waterflux_inst)

    ! !DESCRIPTION:
    !   In principle, Theim equation is applied between 
    !   r = r_e and r = dx (center to center of four
    !   surrounding cells) to take into account the GW pumping.
    !
    !   Assumptions: The head across a given radius is constant 
    !  (i.e., cone of depression is circular) 

    ! !USES:
    use spmdMod         , only : MPI_REAL8, MPI_SUM, mpicom, MPI_INTEGER
    use decompMod       , only : ldecomp, get_proc_global
    use shr_const_mod   , only : SHR_CONST_PI
    use GridcellType    , only : grc
    use clm_time_manager, only : get_step_size, get_curr_date, get_nstep
    use landunit_varcon , only : istwet, istsoil, istice_mec, istcrop

    ! !ARGUMENTS:
    class(groundwater_type)  , intent(inout) :: this
    type(bounds_type)        , intent(in)    :: bounds  
    integer                  , intent(in)    :: num_hydrologyc       ! number of column soil points in column filter
    integer                  , intent(in)    :: filter_hydrologyc(:) ! column filter for soil points
    type(soilhydrology_type) , intent(inout) :: soilhydrology_inst
    type(soilstate_type)     , intent(in)    :: soilstate_inst
    type(waterstate_type)    , intent(inout) :: waterstate_inst
    type(irrigation_type)    , intent(in)    :: irrigation_inst
    type(waterflux_type)     , intent(inout) :: waterflux_inst

    ! !LOCAL VARIABLES:	
 
    character(len=32) :: subname = 'GroundwaterMod' ! subroutine name 
    real(r8), pointer :: zwt_long(:)        , zwt_glob(:)            ! GW_out array for all grid cells; complete grid cell array of GW_out
    real(r8), pointer :: Qn_glob(:)                                  ! Lateral flow in (mm3)
    real(r8), pointer :: g_totCweight(:)
    real(r8), pointer :: g_watsat_long(:)
    real(r8), pointer :: g_sucst_long(:)
    real(r8), pointer :: g_bsw_long(:)
    real(r8), pointer :: g_cellarea_long(:) , g_cellarea_glob(:)
    real(r8), pointer :: qirrig_long(:)
    real(r8), pointer :: GW_ratio_long(:)
    real(r8), pointer :: AqTransmiss_long(:), AqTransmiss_glob(:)
    real(r8), pointer :: qcharge_long(:)

    integer,  pointer :: ZeroHydroCell(:)          ! cells that have no hydrological columns
    integer,  pointer :: ZeroHydroCell_glob(:)     ! global cells that have no hydrological columns

    real(r8) :: rous, aRatio, colArea              ! aquifer yield (-); area ratio of center/neighbor
    real(r8) :: AqTransmissMean                    ! mean aquifer transmissivity of two adjucent cells
    real(r8) :: widMean, lenMean, deltaxMean       ! mean contact width and lenght of two adjucent cells
    real(r8) :: s_y, dummysum, dummysum2
    real(r8) :: TransmissMean
    real(r8) :: pump_tot, pump_layer
    real(r8) :: Qgw_lateral_tot, Qgw_lateral_layer
    integer  :: jwt(bounds%begc:bounds%endc)       ! index of the soil layer right above the water table (-)
    real(r8) :: dtime                              ! land model time step (sec)
	
    integer :: ng, nl, nc, np, nCohorts            ! total number of grid cells,landunits,columns,patches
    integer :: g, c                                ! patch, gridcell, column indices
    integer :: ier                                 ! error code
    integer :: j,fc,i

    integer :: begg, endg       ! beginning and ending gridcell index of current proc
    integer :: begl, endl       ! beginning and ending landunit index of current proc
    integer :: begc, endc       ! beginning and ending column index of current proc
    integer :: begp, endp       ! beginning and ending pft index of current proc

    integer :: year       ! year (0, ...) for nstep
    integer :: month      ! month (1, ..., 12) for nstep
    integer :: day        ! day of month (1, ..., 31) for nstep
    integer :: secs       ! seconds into current date for nstep
    integer :: nstep

    ! Conversion factors

    real(r8), parameter :: km_to_mm    = 1.e6_r8
    real(r8), parameter :: km2_to_mm2  = 1.e12_r8
    real(r8), parameter :: mm_to_m     = 1.e-3_r8
    real(r8), parameter :: m_to_mm     = 1.e3_r8
    real(r8), parameter :: HydroThresh = 0.1_r8
    !-----------------------------------------------------------------------
     associate(&                               
          zi                 =>    col%zi                                , & ! Input:  [real(r8) (:,:) ]  interface level below a "z" level (m)           
          GW_ratio           =>    col%GW_ratio                          , & ! Input:  [real(r8) (:)   ]  USGS GW ratio as irrigation source                                                
                                              
          qflx_irrig         =>    irrigation_inst%qflx_irrig_col        , & ! irrigation flux (mm H2O /s)

          bsw                =>    soilstate_inst%bsw_col                , & ! Input:  [real(r8) (:,:) ]  Clapp and Hornberger "b"                        
          hksat              =>    soilstate_inst%hksat_col              , & ! Input:  [real(r8) (:,:) ]  hydraulic conductivity at saturation (mm H2O /s)
          sucsat             =>    soilstate_inst%sucsat_col             , & ! Input:  [real(r8) (:,:) ]  minimum soil suction (mm)                       
          watsat             =>    soilstate_inst%watsat_col             , & ! Input:  [real(r8) (:,:) ]  volumetric soil water at saturation (porosity)

          zwt                =>    soilhydrology_inst%zwt_col            , & ! Input and Output: [real(r8) (:)   ]  water table depth (m)                                        
          wa                 =>    soilhydrology_inst%wa_col             , & ! Output: [real(r8) (:)   ]  water in the unconfined aquifer (mm)              
          qcharge            =>    soilhydrology_inst%qcharge_col        , & ! Input:  [real(r8) (:)   ]  aquifer recharge rate (mm/s)
          Qgw_lateral        =>    soilhydrology_inst%Qgw_lateral_col    , & ! Output: [real(r8) (:)   ]  GW lateral flow (mm)
          AqTransmiss        =>    soilhydrology_inst%AqTransmiss_col    , & ! Output: [real(r8) (:)   ]  Aquifer Transmissivity(mm2/s)
          Pump_wa            =>    soilhydrology_inst%Pump_wa_col        , & ! Output: [real(r8) (:)   ]  Pumped Water from the aquifer(mm)

          qflx_drain         =>    waterflux_inst%qflx_drain_col         , & ! Input and Output: [real(r8) (:)   ] sub-surface runoff (mm H2O /s)                    

          h2osoi_liq         =>    waterstate_inst%h2osoi_liq_col        & ! Output: [real(r8) (:,:) ] liquid water (kg/m2)
          )
       !-----------------------------------------------
       dtime = get_step_size()
       nstep = get_nstep()
       call get_curr_date (year, month, day, secs)
       ! if (masterproc) then
       !if (iam == 200) then
       !   write(*,*) 'year, month, day, secs, dtime, ', year, month, day, secs, dtime
       !end if    		  
 
       ! Initialize for the mpi_allreduce located between the out and 
       ! in loops
       call get_proc_global(ng=ng, nl=nl, nc=nc, np=np, nCohorts=nCohorts)
       ! Variables to gather from all PEs, while between the out and  
       ! in loops 

       allocate(Qn_glob(ng))
       allocate(ZeroHydroCell(ng))
       allocate(ZeroHydroCell_glob(ng))
       allocate(zwt_long(ng))
       allocate(zwt_glob(ng))
       allocate(g_totCweight(ng))
       allocate(g_watsat_long(ng))
       allocate(g_sucst_long(ng))
       allocate(g_bsw_long(ng))
       allocate(g_cellarea_long(ng))
       allocate(g_cellarea_glob(ng))
       allocate(qirrig_long(ng))
       allocate(GW_ratio_long(ng))
       allocate(AqTransmiss_long(ng))
       allocate(AqTransmiss_glob(ng))
       allocate(qcharge_long(ng))


       ! Initialize to 0 so as to MPI_SUM zeros in all PEs but one per grid cell
       ! FFELFELANI: the length is the total number of gridcell in the entire domain
       ! but only those cells taken care of each processor get value and the rest
       ! remain zero, finally mpi_allreduce would sum them up
       zwt_long(:)                = 0._r8
       Qn_glob(:)                 = 0._r8
       g_totCweight(:)            = 0._r8
       g_watsat_long(:)           = 0._r8
       g_sucst_long(:)            = 0._r8
       g_bsw_long(:)              = 0._r8
       ZeroHydroCell(:)           = 0
       ZeroHydroCell_glob(:)      = 10000
       g_cellarea_long(:)         = 0._r8
       qirrig_long(:)             = 0._r8
       GW_ratio_long(:)           = 0._r8
       AqTransmiss_long(:)        = 0._r8
       qcharge_long(:)            = 0._r8

       !Initialize to 1e36 to help make errors stand out
       zwt_glob(:)           = 1.e36_r8
       g_cellarea_glob(:)    = 1.e36_r8
       AqTransmiss_glob(:)   = 1.e36_r8

       AqTransmiss(bounds%begc:bounds%endc) = this%TransmissivityFromFan(bounds, num_hydrologyc, filter_hydrologyc, soilstate_inst, soilhydrology_inst)

       ! The layer index of the first unsaturated layer, i.e., the layer right above
       ! the water table
       do fc = 1, num_hydrologyc
          c = filter_hydrologyc(fc) 
          jwt(c) = nlevsoi
          ! allow jwt to equal zero when zwt is in top layer  
          do j = 1,nlevsoi
             if(zwt(c) <= zi(c,j)) then
                jwt(c) = j-1 
                exit
             end if
          enddo
       end do

      ! Removing the pumped water from the soil column	   
       if (use_pumping == .true.) then
          do fc = 1, num_hydrologyc
             c = filter_hydrologyc(fc)
             g = col%gridcell(c)
                 !!! use analytical expression for aquifer specific yield
                 rous = watsat(c,nlevsoi) &
                      * ( 1. - (1.+1.e3*zwt(c)/sucsat(c,nlevsoi))**(-1./bsw(c,nlevsoi)))
                 rous=max(rous,0.02_r8)

                 pump_tot = - GW_ratio(c) * qflx_irrig(c) * dtime
                 Pump_wa(c) = GW_ratio(c) * qflx_irrig(c)
                 ! if (nstep > 400 .and. c == 838456) write(*,*) 'c,', c, pump_tot, Pump_wa(c) 
                 !!!--  water table is below the soil column  --------------------------------------
                 if(jwt(c) == nlevsoi) then             
                    wa(c)  = wa(c) + pump_tot
                    zwt(c) = zwt(c) - pump_tot/1000._r8/rous
 
                 else                                
                    !!!-- water table within soil layers 1-9  --------------------------------------
                    !!!============================== pump_tot ========================================= 
                    !!!--  Now remove water via pump_tot

                    !!!should never be positive... but include for completeness
                    if(pump_tot > 0.) then !rising water table

                       call endrun(msg="pump_tot IS POSITIVE in Groundwater!"//errmsg(sourcefile, __LINE__))

                    else ! deepening water table
                       do j = jwt(c)+1, nlevsoi
                          !!! use analytical expression for specific yield
                          s_y = watsat(c,j) &
                               * ( 1. - (1.+1.e3*zwt(c)/sucsat(c,j))**(-1./bsw(c,j)))
                          s_y=max(s_y,0.02_r8)

                          pump_layer=max(pump_tot,-(s_y*(zi(c,j) - zwt(c))*1.e3))
                          pump_layer=min(pump_layer,0._r8)
                          h2osoi_liq(c,j) = h2osoi_liq(c,j) + pump_layer

                          pump_tot = pump_tot - pump_layer

                          if (pump_tot >= 0.) then 
                             zwt(c) = zwt(c) - pump_layer/s_y/1000._r8
                             exit
                          else
                             zwt(c) = zi(c,j)
                          endif
                       enddo

                       !!!--  remove residual pump_tot  ---------------------------------------------
                       zwt(c) = zwt(c) - pump_tot/1000._r8/rous
                       wa(c)  = wa(c) + pump_tot
                    endif

                    !!!-- recompute jwt  ---------------------------------------------------------
                    !!! allow jwt to equal zero when zwt is in top layer
                    jwt(c) = nlevsoi
                    do j = 1,nlevsoi
                       if(zwt(c) <= zi(c,j)) then
                          jwt(c) = j-1
                          exit
                       end if
                    enddo
                 end if! end of jwt if construct

                 zwt(c) = max(0.0_r8,zwt(c))
                 zwt(c) = min(80._r8,zwt(c))    
          end do
       end if
       ! To get the weighted average ZWT accross all columns within each cell  
       ! loop over the columns that each processor handles
       ! do c = bounds%begc,bounds%endc
       do fc = 1, num_hydrologyc
          c = filter_hydrologyc(fc)
          g = col%gridcell(c)
         ! Saving for mpi_allreduce coming up
          g_totCweight(g)    = g_totCweight(g)     + col%wtgcell(c)
          g_cellarea_long(g) = g_cellarea_long(g)  + col%wtgcell(c) * grc%area(g)

          zwt_long(g)          = zwt_long(g)         + col%wtgcell(c) * zwt(c)  ! GW_od goes to history
          AqTransmiss_long(g)  = AqTransmiss_long(g) + col%wtgcell(c) * AqTransmiss(c)

          qirrig_long(g)   = qirrig_long(g)   + col%wtgcell(c) * qflx_irrig(c)
          GW_ratio_long(g) = GW_ratio_long(g) + col%wtgcell(c) * GW_ratio(c)
          qcharge_long(g)  = qcharge_long(g)  + col%wtgcell(c) * qcharge(c)

          g_watsat_long(g)  = g_watsat_long(g) + col%wtgcell(c) * watsat(c,nlevsoi) 
          g_sucst_long(g)   = g_sucst_long(g)  + col%wtgcell(c) * sucsat(c,nlevsoi) 
          g_bsw_long(g)     = g_bsw_long(g)    + col%wtgcell(c) * bsw(c,nlevsoi)
       end do  

       do g = bounds%begg,bounds%endg
          if (g_totCweight(g) <= 1.000001_r8 .and. g_totCweight(g) > HydroThresh) then

              zwt_long(g)         = zwt_long(g) / g_totCweight(g)
              AqTransmiss_long(g) = AqTransmiss_long(g)/ g_totCweight(g)
              qirrig_long(g)      = qirrig_long(g) / g_totCweight(g)
              GW_ratio_long(g)    = GW_ratio_long(g) / g_totCweight(g)
              qcharge_long(g)     = qcharge_long(g)/ g_totCweight(g)
              g_watsat_long(g)    = g_watsat_long(g) / g_totCweight(g)
              g_sucst_long(g)     = g_sucst_long(g) / g_totCweight(g)
              g_bsw_long(g)       = g_bsw_long(g) / g_totCweight(g)
          else
              ZeroHydroCell(g) = 1
          end if
       end do


       ! Need action: what happens to the zwt_long of cells g_totCweight(g) < 0.01???????????????
	
       ! ------------------------------------------------
       ! Between the out and in loops
       ! all processors share all the relevant data needed to come up with GW_id
       ! Could the same be done in if (iam == 0) followed by call mpi_bcast?
       ! ------------------------------------------------

       ! Gather zwt_glob from zwt_long, ultimately
       !                  from GW_od

       ! I do MPI_SUM because each processor only has values for the patch/clump its dealing  
       ! with and the rest of the array is zero   
       call mpi_allreduce(zwt_long, zwt_glob, ng, &
                          MPI_REAL8, MPI_SUM, mpicom, ier) 

       call mpi_allreduce(AqTransmiss_long, AqTransmiss_glob, ng, &
                          MPI_REAL8, MPI_SUM, mpicom, ier) 

       call mpi_allreduce(g_cellarea_long, g_cellarea_glob, ng, &
                          MPI_REAL8, MPI_SUM, mpicom, ier)

       ! we need to exclude those cells that have no hydrologically active columns  
       call mpi_allreduce(ZeroHydroCell, ZeroHydroCell_glob, ng, &
                          MPI_INTEGER, MPI_SUM, mpicom, ier)
						  
       call mpi_barrier(mpicom,ier)

       ! gathering the information from the neigboring cells.
       do  g = 1, ng
          ! The GW lateral flow is ruled by Darcy's
          if (ZeroHydroCell_glob(g)== 0) then

             if (ldecomp%gtoplft(g) <= ng .and. ldecomp%gtoplft(g) >= 1 .and. ZeroHydroCell_glob(ldecomp%gtoplft(g)) == 0) then
                 AqTransmissMean = (AqTransmiss_glob(g) + AqTransmiss_glob(ldecomp%gtoplft(g)))/2._r8
                 lenMean = (sqrt(g_cellarea_glob(g)) + sqrt(g_cellarea_glob(ldecomp%gtoplft(g)))) * km_to_mm * sqrt(2._r8) / 2._r8
                 deltaxMean = (sqrt(g_cellarea_glob(g)) + sqrt(g_cellarea_glob(ldecomp%gtoplft(g)))) * km_to_mm / 2._r8
                 widMean = deltaxMean * sqrt(0.5_r8 * tan(rpi/8._r8))
                 Qn_glob(g) = Qn_glob(g) + widMean * AqTransmissMean * (zwt_glob(g) - zwt_glob(ldecomp%gtoplft(g))) * m_to_mm * dtime / lenMean

             end if

             if (ldecomp%gtop(g) <= ng .and. ldecomp%gtop(g) >= 1 .and. ZeroHydroCell_glob(ldecomp%gtop(g)) == 0) then
                 AqTransmissMean = (AqTransmiss_glob(g) + AqTransmiss_glob(ldecomp%gtop(g)))/2._r8
                 lenMean = (sqrt(g_cellarea_glob(g)) + sqrt(g_cellarea_glob(ldecomp%gtop(g)))) * km_to_mm / 2._r8
                 widMean = lenMean * sqrt(0.5_r8 * tan(rpi/8._r8))
                 Qn_glob(g) = Qn_glob(g) + widMean * AqTransmissMean * (zwt_glob(g) - zwt_glob(ldecomp%gtop(g))) * m_to_mm * dtime / lenMean

             end if

             if (ldecomp%gtoprgt(g) <= ng .and. ldecomp%gtoprgt(g) >= 1 .and. ZeroHydroCell_glob(ldecomp%gtoprgt(g)) == 0) then
                 AqTransmissMean = (AqTransmiss_glob(g) + AqTransmiss_glob(ldecomp%gtoprgt(g)))/2._r8
                 lenMean = (sqrt(g_cellarea_glob(g)) + sqrt(g_cellarea_glob(ldecomp%gtoprgt(g)))) * km_to_mm * sqrt(2._r8) / 2._r8
                 deltaxMean = (sqrt(g_cellarea_glob(g)) + sqrt(g_cellarea_glob(ldecomp%gtoprgt(g)))) * km_to_mm / 2._r8
                 widMean = deltaxMean * sqrt(0.5_r8 * tan(rpi/8._r8))
                 Qn_glob(g) = Qn_glob(g) + widMean * AqTransmissMean * (zwt_glob(g) - zwt_glob(ldecomp%gtoprgt(g))) * m_to_mm * dtime / lenMean

             end if

             if (ldecomp%grgt(g) <= ng .and. ldecomp%grgt(g) >= 1 .and. ZeroHydroCell_glob(ldecomp%grgt(g)) == 0) then
                 AqTransmissMean = (AqTransmiss_glob(g) + AqTransmiss_glob(ldecomp%grgt(g)))/2._r8
                 lenMean = (sqrt(g_cellarea_glob(g)) + sqrt(g_cellarea_glob(ldecomp%grgt(g)))) * km_to_mm / 2._r8
                 widMean = lenMean * sqrt(0.5_r8 * tan(rpi/8._r8))
                 Qn_glob(g) = Qn_glob(g) + widMean * AqTransmissMean * (zwt_glob(g) - zwt_glob(ldecomp%grgt(g))) * m_to_mm * dtime / lenMean

             end if	

             if (ldecomp%gbotrgt(g) <= ng .and. ldecomp%gbotrgt(g) >= 1 .and. ZeroHydroCell_glob(ldecomp%gbotrgt(g)) == 0) then
                 AqTransmissMean = (AqTransmiss_glob(g) + AqTransmiss_glob(ldecomp%gbotrgt(g)))/2._r8
                 lenMean = (sqrt(g_cellarea_glob(g)) + sqrt(g_cellarea_glob(ldecomp%gbotrgt(g)))) * km_to_mm * sqrt(2._r8) / 2._r8
                 deltaxMean = (sqrt(g_cellarea_glob(g)) + sqrt(g_cellarea_glob(ldecomp%gbotrgt(g)))) * km_to_mm / 2._r8
                 widMean = deltaxMean * sqrt(0.5_r8 * tan(rpi/8._r8))
                 Qn_glob(g) = Qn_glob(g) + widMean * AqTransmissMean * (zwt_glob(g) - zwt_glob(ldecomp%gbotrgt(g))) * m_to_mm * dtime / lenMean

             end if

             if (ldecomp%gbot(g) <= ng .and. ldecomp%gbot(g) >= 1 .and. ZeroHydroCell_glob(ldecomp%gbot(g)) == 0) then
                 AqTransmissMean = (AqTransmiss_glob(g) + AqTransmiss_glob(ldecomp%gbot(g)))/2._r8
                 lenMean = (sqrt(g_cellarea_glob(g)) + sqrt(g_cellarea_glob(ldecomp%gbot(g)))) * km_to_mm / 2._r8
                 widMean = lenMean * sqrt(0.5_r8 * tan(rpi/8._r8))
                 Qn_glob(g) = Qn_glob(g) + widMean * AqTransmissMean * (zwt_glob(g) - zwt_glob(ldecomp%gbot(g))) * m_to_mm * dtime / lenMean

             end if

             if (ldecomp%gbotlft(g) <= ng .and. ldecomp%gbotlft(g) >= 1 .and. ZeroHydroCell_glob(ldecomp%gbotlft(g)) == 0) then
                 AqTransmissMean = (AqTransmiss_glob(g) + AqTransmiss_glob(ldecomp%gbotlft(g)))/2._r8
                 lenMean = (sqrt(g_cellarea_glob(g)) + sqrt(g_cellarea_glob(ldecomp%gbotlft(g)))) * km_to_mm * sqrt(2._r8) / 2._r8
                 deltaxMean = (sqrt(g_cellarea_glob(g)) + sqrt(g_cellarea_glob(ldecomp%gbotlft(g)))) * km_to_mm / 2._r8
                 widMean = deltaxMean * sqrt(0.5_r8 * tan(rpi/8._r8))
                 Qn_glob(g) = Qn_glob(g) + widMean * AqTransmissMean * (zwt_glob(g) - zwt_glob(ldecomp%gbotlft(g))) * m_to_mm * dtime / lenMean

             end if

             if (ldecomp%glft(g) <= ng .and. ldecomp%glft(g) >= 1 .and. ZeroHydroCell_glob(ldecomp%glft(g)) == 0) then
                 AqTransmissMean = (AqTransmiss_glob(g) + AqTransmiss_glob(ldecomp%glft(g)))/2._r8
                 lenMean = (sqrt(g_cellarea_glob(g)) + sqrt(g_cellarea_glob(ldecomp%glft(g)))) * km_to_mm / 2._r8
                 widMean = lenMean * sqrt(0.5_r8 * tan(rpi/8._r8))
                 Qn_glob(g) = Qn_glob(g) + widMean * AqTransmissMean * (zwt_glob(g) - zwt_glob(ldecomp%glft(g))) * m_to_mm * dtime / lenMean

             end if

          ! IF there is pumping, the GW lateral flow is ruled by Combination of Darcy's and Theim
          ! else if (ZeroHydroCell_glob(g)== 0 .and. GW_ratio_long(g) * qirrig_long(g) > 0._r8) then

          end if
       end do
   
       ! Checking the lateral water balance (in terms of volume) 
       ! if (iam == 200) then
          ! dummysum = 0._r8
          ! do  g = 1, ng
              ! dummysum = dummysum + Qn_glob(g)
          ! end do
          ! write(*,*) 'Felfelani: this is the sum of the Latera GW; ', dummysum
       ! end if
       ! convert m to mm water


       do fc = 1, num_hydrologyc
          c = filter_hydrologyc(fc)
          g = col%gridcell(c)

          aRatio  = col%wtgcell(c) / g_totCweight(g)
          colArea = col%wtgcell(c) * grc%area(g) * km2_to_mm2
          Qgw_lateral(c) = Qn_glob(g) * aRatio / colArea  !unit is mm 

          rous = watsat(c,nlevsoi) &
               * ( 1. - (1.+1.e3*zwt(c)/sucsat(c,nlevsoi))**(-1./bsw(c,nlevsoi)))
          rous=max(rous,0.02_r8)
          ! if (nstep  > 400 .and. c == 838456) write(*,*) 'c,', c, Qgw_lateral(c)
          if (Qgw_lateral(c) > 0._r8) then

                Qgw_lateral_tot = Qgw_lateral(c) * dtime * 1._r8
                if(jwt(c) == nlevsoi) then             
                   wa(c)  = wa(c) + Qgw_lateral_tot
                   zwt(c) = zwt(c) - Qgw_lateral_tot/1000._r8/rous
                else   
                   do j = jwt(c)+1, 1,-1
                       !! use analytical expression for specific yield
                       s_y = watsat(c,j) &
                            * ( 1. -  (1.+1.e3*zwt(c)/sucsat(c,j))**(-1./bsw(c,j)))
                       s_y=max(s_y,0.02_r8)

                       Qgw_lateral_layer=min(Qgw_lateral_tot,(s_y*(zwt(c) - zi(c,j-1))*1.e3))
                       Qgw_lateral_layer=max(Qgw_lateral_layer,0._r8)

                       h2osoi_liq(c,j) = h2osoi_liq(c,j) + Qgw_lateral_layer
				
                       if(s_y > 0._r8) zwt(c) = zwt(c) - Qgw_lateral_layer/s_y/1000._r8

                       Qgw_lateral_tot = Qgw_lateral_tot - Qgw_lateral_layer
                       if (Qgw_lateral_tot <= 0.) exit
                   enddo

                   ! add residual to the sub-surface runoff (both lateral flow and runoff are positive here)
				   qflx_drain(c) = qflx_drain(c) + Qgw_lateral_tot / dtime


                end if

          else if (Qgw_lateral(c) < 0._r8) then

              Qgw_lateral_tot = Qgw_lateral(c) * dtime * 1._r8
              !! --  water table is below the soil column  -------------------------------------- 
              if(jwt(c) == nlevsoi) then             
                 wa(c)  = wa(c) + Qgw_lateral_tot
                 zwt(c) = zwt(c) - Qgw_lateral_tot/1000._r8/rous
              else                                
                 !! -- water table within soil layers 1-9  -------------------------------------
                 !! ============================== Qgw_lateral_tot ========================================= 
                 !! --  Now remove water via Qgw_lateral_tot

                 !! should never be positive... but include for completeness 
                 if(Qgw_lateral_tot > 0.) then !rising water table

                    call endrun(msg="Qgw_lateral_tot IS POSITIVE in Groundwater!"//errmsg(sourcefile, __LINE__))

                 else ! deepening water table
                    do j = jwt(c)+1, nlevsoi
                       !! use analytical expression for specific yield
                       s_y = watsat(c,j) &
                            * ( 1. - (1.+1.e3*zwt(c)/sucsat(c,j))**(-1./bsw(c,j)))
                       s_y=max(s_y,0.02_r8)

                       Qgw_lateral_layer=max(Qgw_lateral_tot,-(s_y*(zi(c,j) - zwt(c))*1.e3))
                       Qgw_lateral_layer=min(Qgw_lateral_layer,0._r8)
                       h2osoi_liq(c,j) = h2osoi_liq(c,j) + Qgw_lateral_layer

                       Qgw_lateral_tot = Qgw_lateral_tot - Qgw_lateral_layer

                       if (Qgw_lateral_tot >= 0.) then 
                         zwt(c) = zwt(c) - Qgw_lateral_layer/s_y/1000._r8
                          exit
                       else
                          zwt(c) = zi(c,j)
                       endif
                    enddo

                    !! --  remove residual Qgw_lateral_tot  ---------------------------------------------
                    zwt(c) = zwt(c) - Qgw_lateral_tot/1000._r8/rous
                    wa(c) = wa(c) + Qgw_lateral_tot
                 endif

                 !! -- recompute jwt  ---------------------------------------------------------
                 !! allow jwt to equal zero when zwt is in top layer
                 jwt(c) = nlevsoi
                 do j = 1,nlevsoi
                    if(zwt(c) <= zi(c,j)) then
                       jwt(c) = j-1
                       exit
                    end if
                 enddo
              end if! end of jwt if construct

              zwt(c) = max(0.0_r8,zwt(c))
              zwt(c) = min(80._r8,zwt(c))
          end if
       end do

       deallocate(Qn_glob)
       deallocate(ZeroHydroCell)
       deallocate(ZeroHydroCell_glob)
       deallocate(zwt_long)
       deallocate(zwt_glob)
       deallocate(g_totCweight)
       deallocate(g_watsat_long)
       deallocate(g_sucst_long)
       deallocate(g_bsw_long)
       deallocate(g_cellarea_long)
       deallocate(g_cellarea_glob)
       deallocate(qirrig_long)
       deallocate(GW_ratio_long)
       deallocate(AqTransmiss_long)
       deallocate(AqTransmiss_glob)
       deallocate(qcharge_long)  
     end associate
  end subroutine UpdateGWFanLatPump
  !-----------------------------------------------------------------------
  subroutine UpdateGWDefaultPump(this, bounds, num_hydrologyc, filter_hydrologyc, &
        soilhydrology_inst, soilstate_inst,waterstate_inst, irrigation_inst)

    ! !DESCRIPTION:
    !   In principle, Theim equation is applied between
    !   r = r_e and r = dx (center to center of four
    !   surrounding cells) to take into account the GW pumping.
    !
    !   Assumptions: The head across a given radius is constant
    !  (i.e., cone of depression is circular)

    ! !USES:
    use spmdMod         , only : MPI_REAL8, MPI_SUM, mpicom, MPI_INTEGER
    use decompMod       , only : ldecomp, get_proc_global
    use shr_const_mod   , only : SHR_CONST_PI
    use GridcellType    , only : grc
    use clm_time_manager, only : get_step_size, get_curr_date, get_nstep
    use landunit_varcon , only : istwet, istsoil, istice_mec, istcrop

    ! !ARGUMENTS:
    class(groundwater_type)  , intent(inout) :: this
    type(bounds_type)        , intent(in)    :: bounds  
    integer                  , intent(in)    :: num_hydrologyc       ! number of column soil points in column filter
    integer                  , intent(in)    :: filter_hydrologyc(:) ! column filter for soil points
    type(soilhydrology_type) , intent(inout) :: soilhydrology_inst
    type(soilstate_type)     , intent(in)    :: soilstate_inst
    type(waterstate_type)    , intent(inout) :: waterstate_inst
    type(irrigation_type)    , intent(in)    :: irrigation_inst

    ! !LOCAL VARIABLES:
 
    character(len=32) :: subname = 'GroundwaterMod' ! subroutine name

    real(r8) :: rous, aRatio, colArea              ! aquifer yield (-); area ratio of center/neighbor
    real(r8) :: s_y, dummysum, dummysum2
    real(r8) :: pump_tot, pump_layer
    integer  :: jwt(bounds%begc:bounds%endc)       ! index of the soil layer right above the water table (-)
    real(r8) :: dtime                              ! land model time step (sec)


    integer :: ng, nl, nc, np, nCohorts            ! total number of grid cells,landunits,columns,patches
    integer :: g, c                                ! patch, gridcell, column indices
    integer :: ier                                 ! error code
    integer :: j,fc,i

    integer :: begg, endg       ! beginning and ending gridcell index of current proc
    integer :: begl, endl       ! beginning and ending landunit index of current proc
    integer :: begc, endc       ! beginning and ending column index of current proc
    integer :: begp, endp       ! beginning and ending pft index of current proc

    integer :: year       ! year (0, ...) for nstep
    integer :: month      ! month (1, ..., 12) for nstep
    integer :: day        ! day of month (1, ..., 31) for nstep
    integer :: secs       ! seconds into current date for nstep
    integer :: nstep

    ! Conversion factors

    real(r8), parameter :: km_to_mm    = 1.e6_r8
    real(r8), parameter :: km2_to_mm2  = 1.e12_r8
    real(r8), parameter :: mm_to_m     = 1.e-3_r8
    real(r8), parameter :: m_to_mm     = 1.e3_r8
    real(r8), parameter :: HydroThresh = 0.1_r8
    !-----------------------------------------------------------------------
     associate(&                              
          zi                 =>    col%zi                                , & ! Input:  [real(r8) (:,:) ]  interface level below a "z" level (m)          
          GW_ratio           =>    col%GW_ratio                          , & ! Input:  [real(r8) (:)   ]  USGS GW ratio as irrigation source                                                
                                             
          qflx_irrig         =>    irrigation_inst%qflx_irrig_col        , & ! irrigation flux (mm H2O /s)

          bsw                =>    soilstate_inst%bsw_col                , & ! Input:  [real(r8) (:,:) ]  Clapp and Hornberger "b"                        
          hksat              =>    soilstate_inst%hksat_col              , & ! Input:  [real(r8) (:,:) ]  hydraulic conductivity at saturation (mm H2O /s)
          sucsat             =>    soilstate_inst%sucsat_col             , & ! Input:  [real(r8) (:,:) ]  minimum soil suction (mm)                      
          watsat             =>    soilstate_inst%watsat_col             , & ! Input:  [real(r8) (:,:) ]  volumetric soil water at saturation (porosity)

          zwt                =>    soilhydrology_inst%zwt_col                , & ! Input and Output: [real(r8) (:)   ]  water table depth (m)
          wa                 =>    soilhydrology_inst%wa_col                 , & ! Output: [real(r8) (:)   ]  water in the unconfined aquifer (mm)
          qcharge            =>    soilhydrology_inst%qcharge_col            , & ! Input:  [real(r8) (:)   ]  aquifer recharge rate (mm/s)
          Qgw_lateral        =>    soilhydrology_inst%Qgw_lateral_col        , & ! Output: [real(r8) (:)   ]  GW lateral flow (mm/s)
          AqTransmiss        =>    soilhydrology_inst%AqTransmiss_col        , & ! Output: [real(r8) (:)   ]  Aquifer Transmissivity(mm2/s)
          Pump_wa            =>    soilhydrology_inst%Pump_wa_col            , & ! Output: [real(r8) (:)   ]  Pumped Water from the aquifer(mm/s)
          QlatField_north    =>    soilhydrology_inst%QlatField_northing_grc , & !  Output: [real(r8) (:)   ] Northward lateral GW flow
          QlatField_east     =>    soilhydrology_inst%QlatField_easting_grc  , & !  Output: [real(r8) (:)   ] Eastward lateral GW flow

          h2osoi_liq         =>    waterstate_inst%h2osoi_liq_col              & ! Output: [real(r8) (:,:) ] liquid water (kg/m2)
          )
       !-----------------------------------------------
       dtime = get_step_size()
       nstep = get_nstep()
       call get_curr_date (year, month, day, secs)
       call get_proc_global(ng=ng, nl=nl, nc=nc, np=np, nCohorts=nCohorts)


       ! The layer index of the first unsaturated layer, i.e., the layer right above
       ! the water table
       do fc = 1, num_hydrologyc
          c = filter_hydrologyc(fc)
          jwt(c) = nlevsoi
          ! allow jwt to equal zero when zwt is in top layer  
          do j = 1,nlevsoi
             if(zwt(c) <= zi(c,j)) then
                jwt(c) = j-1
                exit
             end if
          enddo
       end do

      ! Removing the pumped water from the soil column  
       if (use_pumping == .true.) then
          do fc = 1, num_hydrologyc
             c = filter_hydrologyc(fc)
             g = col%gridcell(c)
                 !!! use analytical expression for aquifer specific yield
                 rous = watsat(c,nlevsoi) &
                      * ( 1. - (1.+1.e3*zwt(c)/sucsat(c,nlevsoi))**(-1./bsw(c,nlevsoi)))
                 rous=max(rous,0.02_r8)

                 pump_tot = - GW_ratio(c) * qflx_irrig(c) * dtime
                 Pump_wa(c) = GW_ratio(c) * qflx_irrig(c)
                 ! if (nstep > 400 .and. c == 838456) write(*,*) 'c,', c, pump_tot, Pump_wa(c)
                 !!!--  water table is below the soil column  --------------------------------------
                 if(jwt(c) == nlevsoi) then            
                    wa(c)  = wa(c) + pump_tot
                    zwt(c) = zwt(c) - pump_tot/1000._r8/rous
 
                 else                                
                    !!!-- water table within soil layers 1-9  --------------------------------------
                    !!!============================== pump_tot =========================================
                    !!!--  Now remove water via pump_tot

                    !!!should never be positive... but include for completeness
                    if(pump_tot > 0.) then !rising water table

                       call endrun(msg="pump_tot IS POSITIVE in Groundwater!"//errmsg(sourcefile, __LINE__))

                    else ! deepening water table
                       do j = jwt(c)+1, nlevsoi
                          !!! use analytical expression for specific yield
                          s_y = watsat(c,j) &
                               * ( 1. - (1.+1.e3*zwt(c)/sucsat(c,j))**(-1./bsw(c,j)))
                          s_y=max(s_y,0.02_r8)

                          pump_layer=max(pump_tot,-(s_y*(zi(c,j) - zwt(c))*1.e3))
                          pump_layer=min(pump_layer,0._r8)
                          h2osoi_liq(c,j) = h2osoi_liq(c,j) + pump_layer

                          pump_tot = pump_tot - pump_layer

                          if (pump_tot >= 0.) then
                             zwt(c) = zwt(c) - pump_layer/s_y/1000._r8
                             exit
                          else
                             zwt(c) = zi(c,j)
                          endif
                       enddo

                       !!!--  remove residual pump_tot  ---------------------------------------------
                       zwt(c) = zwt(c) - pump_tot/1000._r8/rous
                       wa(c)  = wa(c) + pump_tot
                    endif

                    !!!-- recompute jwt  ---------------------------------------------------------
                    !!! allow jwt to equal zero when zwt is in top layer
                    jwt(c) = nlevsoi
                    do j = 1,nlevsoi
                       if(zwt(c) <= zi(c,j)) then
                          jwt(c) = j-1
                          exit
                       end if
                    enddo
                 end if! end of jwt if construct

                 zwt(c) = max(0.0_r8,zwt(c))
                 ! zwt(c) = min(80._r8,zwt(c))    
          end do
       end if

     end associate
  end subroutine UpdateGWDefaultPump
  !-----------------------------------------------------------------------

  function TransmissivityFromFan(this, bounds, num_hydrologyc, filter_hydrologyc, &
                                 soilstate_inst, soilhydrology_inst) &
    result(Transmiss)

    ! !DESCRIPTION:
    !  Calculating the transmissivity based on the Fan et al. (2007)
    !  and Zeng et al. (2016)

    ! !ARGUMENTS:
    class(groundwater_type)  , intent(in)    :: this
    type(bounds_type)        , intent(in)    :: bounds
    integer                  , intent(in)    :: num_hydrologyc       ! number of column soil points in column filter
    integer                  , intent(in)    :: filter_hydrologyc(:) ! column filter for soil points
    type(soilstate_type)     , intent(in)    :: soilstate_inst
    type(soilhydrology_type) , intent(in)    :: soilhydrology_inst

    ! !LOCAL VARIABLES:
    integer  :: c,j,fc,i,g                                   ! indices
    real(r8) :: e_folding_length
    real(r8) :: Transmiss(bounds%begc:bounds%endc)           ! Transmissivity of the aquifer (mm^2/s)
    real(r8) :: beta_rad                                     ! terrain slope (rad)
    integer  :: jwt2(bounds%begc:bounds%endc)                ! index of the soil layer right above the water table (-)

    real(r8), parameter :: m_to_mm = 1.e3_r8	     

    associate(                                                             &                              
          zi                 =>    col%zi                                , & ! Input:  [real(r8) (:,:) ] interface level below a "z" level (m)           
          dz                 =>    col%dz                                , & ! Input:  [real(r8) (:,:) ] layer depth (m)  

          hksat              =>    soilstate_inst%hksat_col              , & ! Input:  [real(r8) (:,:) ]  hydraulic conductivity at saturation (mm H2O /s)                
          clayP              =>    soilstate_inst%cellclay_col           , & ! Input:  [real(r8) (:,:) ] percent clay (0 < ~ < 100)-----> It is not a fraction!!!
          zwt                =>    soilhydrology_inst%zwt_col              & ! Input: [real(r8) (:)   ]  water table depth (m)
          )
    !-----------------------------------------------------------------------

    ! The layer index of the first unsaturated layer, i.e., the layer right above
    ! the water table
    do fc = 1, num_hydrologyc
       c = filter_hydrologyc(fc)
       jwt2(c) = 100 ! arbitrarily 100 for water table below the soil column 
       ! allow jwt2 to equal zero when zwt is in top layer
       do j = 1,nlevsoi
          if(zwt(c) <= zi(c,j)) then
             jwt2(c) = j
             exit
          end if
       enddo
    end do

    Transmiss(:) = 0._r8
    do fc = 1, num_hydrologyc
       c = filter_hydrologyc(fc)
       g = col%gridcell(c)

       ! Zeng et al.(2018): e_folding_length (m) calculations for bedrock	 
       !beta_rad = (rpi/180.) * col%topo_slope(c)
       !if (beta_rad <= 0.16) then
       !    e_folding_length = 20._r8/(1._r8 + 125._r8*beta_rad)
       !else if (beta_rad > 0.16) then
       !    e_folding_length = 1._r8
       !end if

       ! Fan et al.(2007): e_folding_length (m) calculations for regolith
       beta_rad = (rpi/180.) * col%topo_slope(c)
       if (beta_rad <= 0.16) then
           e_folding_length = 120._r8/(1._r8 + 150._r8*beta_rad)
       else if (beta_rad > 0.16) then
           e_folding_length = 5._r8
       end if



       if (jwt2(c) < nlevsoi) then
          Transmiss(c) = clayP(c,jwt2(c)) * hksat(c,jwt2(c)) * (zi(c,jwt2(c))-zwt(c)) * m_to_mm
          do j = jwt2(c)+1,nlevsoi
             Transmiss(c) = Transmiss(c) + clayP(c,j) * hksat(c,j) * dz(c,j) * m_to_mm
          end do
          Transmiss(c) = Transmiss(c) + clayP(c,nlevsoi) * hksat(c,nlevsoi) * e_folding_length * m_to_mm

       else if (jwt2(c) .eq. nlevsoi) then
          Transmiss(c) = clayP(c,nlevsoi) * hksat(c,nlevsoi) * (zi(c,nlevsoi)-zwt(c)) * m_to_mm &
                         + clayP(c,nlevsoi) * hksat(c,nlevsoi) * e_folding_length * m_to_mm

       else if (jwt2(c) .eq. 100) then
          Transmiss(c) = clayP(c,nlevsoi) * hksat(c,nlevsoi) * e_folding_length * m_to_mm &
                         * exp((zi(c,nlevsoi)-zwt(c))/e_folding_length)
       end if

       !if (grc%latdeg(g) < 22.0 .and. grc%latdeg(g) > 21.0 .and. grc%londeg(g) < 243.0 .and. grc%londeg(g) > 242.0) then

       ! if ((c .eq. 1662431) .or. (c .eq. 1486770) .or. (c .eq. 1398102)) then
           ! write(*,*)  '--------------------------------------------------------------' 
           ! write(*,*)  'Felfelani Transmiss: Processor Num, g, c, lat(g), lon(g), grc%area(g)' 
           ! write(*,*)   iam, g, c, grc%latdeg(g), grc%londeg(g), grc%area(g)
           ! write(*,*)  'Transmiss(c), clayP(c,1), clayP(c,nlevsoi), hksat(c,1), hksat(c,nlevsoi)' 
           ! write(*,*)   Transmiss(c), clayP(c,1), clayP(c,nlevsoi), hksat(c,1), hksat(c,nlevsoi) 
           ! write(*,*)  'e_folding_length,zi(c,nlevsoi),zwt(c),beta_rad,col%topo_slope(c)'
           ! write(*,*)   e_folding_length,zi(c,nlevsoi),zwt(c),beta_rad,col%topo_slope(c)
           ! write(*,*)  '--------------------------------------------------------------'
       ! end if 

    end do
    end associate
  end function TransmissivityFromFan

end module GroundwaterMod