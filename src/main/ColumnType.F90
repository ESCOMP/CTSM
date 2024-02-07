module ColumnType

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Column data type allocation and initialization
  ! -------------------------------------------------------- 
  ! column types can have values of
  ! -------------------------------------------------------- 
  !   1  => (istsoil)          soil (vegetated or bare soil)
  !   2  => (istcrop)          crop (only for crop configuration)
  !   3  => (UNUSED)           (formerly non-multiple elevation class land ice; currently unused)
  !   4  => (istice)           land ice
  !   5  => (istdlak)          deep lake
  !   6  => (istwet)           wetland
  !   71 => (icol_roof)        urban roof
  !   72 => (icol_sunwall)     urban sunwall
  !   73 => (icol_shadewall)   urban shadewall
  !   74 => (icol_road_imperv) urban impervious road
  !   75 => (icol_road_perv)   urban pervious road
  !
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
  use clm_varpar     , only : nlevsno, nlevgrnd, nlevlak, nlevmaxurbgrnd
  use clm_varcon     , only : spval, ispval
  use shr_sys_mod    , only : shr_sys_abort
  use clm_varctl     , only : iulog
  use column_varcon  , only : is_hydrologically_active
  use LandunitType   , only : lun
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  private
  !
  type, public :: column_type
     ! g/l/c/p hierarchy, local g/l/c/p cells only
     integer , pointer :: landunit             (:)   ! index into landunit level quantities
     real(r8), pointer :: wtlunit              (:)   ! weight (relative to landunit)
     integer , pointer :: gridcell             (:)   ! index into gridcell level quantities
     real(r8), pointer :: wtgcell              (:)   ! weight (relative to gridcell)
     integer , pointer :: patchi               (:)   ! beginning patch index for each column
     integer , pointer :: patchf               (:)   ! ending patch index for each column
     integer , pointer :: npatches             (:)   ! number of patches for each column

     ! topological mapping functionality
     integer , pointer :: itype                (:)   ! column type (after init, should only be modified via update_itype routine)
     integer , pointer :: lun_itype            (:)   ! landunit type (col%lun_itype(ci) is the
                                                     ! same as lun%itype(col%landunit(ci)), but is often a more convenient way to access this type
     logical , pointer :: active               (:)   ! true=>do computations on this column
     logical , pointer :: type_is_dynamic      (:)   ! true=>itype can change throughout the run
     
     logical , pointer :: is_fates             (:)   ! .true. -> this is a fates column
                                                     ! .false. -> this is NOT a fates column
     
     ! topography
     ! TODO(wjs, 2016-04-05) Probably move these things into topoMod
     real(r8), pointer :: micro_sigma          (:)   ! microtopography pdf sigma (m)
     real(r8), pointer :: topo_slope           (:)   ! gridcell topographic slope
     real(r8), pointer :: topo_std             (:)   ! gridcell elevation standard deviation

     ! vertical levels
     integer , pointer :: snl                  (:)   ! number of snow layers
     real(r8), pointer :: dz                   (:,:) ! layer thickness (m)  (-nlevsno+1:nlevgrnd) 
     real(r8), pointer :: z                    (:,:) ! layer depth (m) (-nlevsno+1:nlevgrnd) 
     real(r8), pointer :: zi                   (:,:) ! interface level below a "z" level (m) (-nlevsno+0:nlevgrnd) 
     real(r8), pointer :: zii                  (:)   ! convective boundary height [m]
     real(r8), pointer :: dz_lake              (:,:) ! lake layer thickness (m)  (1:nlevlak)
     real(r8), pointer :: z_lake               (:,:) ! layer depth for lake (m)
     real(r8), pointer :: lakedepth            (:)   ! variable lake depth (m)                             
     integer , pointer :: nbedrock             (:)   ! variable depth to bedrock index
     ! hillslope hydrology variables
     integer,  pointer :: col_ndx              (:)   ! column index of column (hillslope hydrology)
     integer,  pointer :: colu                 (:)   ! column index of uphill column (hillslope hydrology)
     integer,  pointer :: cold                 (:)   ! column index of downhill column (hillslope hydrology)
     integer,  pointer :: hillslope_ndx        (:)   ! hillslope identifier
     real(r8), pointer :: hill_elev            (:)   ! mean elevation of column relative to stream channel (m)
     real(r8), pointer :: hill_slope           (:)   ! mean along-hill slope (m/m)
     real(r8), pointer :: hill_area            (:)   ! mean surface area (m2)
     real(r8), pointer :: hill_width           (:)   ! across-hill width of bottom boundary of column (m)
     real(r8), pointer :: hill_distance        (:)   ! along-hill distance of column from bottom of hillslope (m)
     real(r8), pointer :: hill_aspect          (:)   ! azimuth angle of column wrt to north, positive to east (radians)

     ! other column characteristics
     logical , pointer :: is_hillslope_column(:)     ! true if this column is a hillslope element
     logical , pointer :: hydrologically_active(:)   ! true if this column is a hydrologically active type
     logical , pointer :: urbpoi               (:)   ! true=>urban point

     ! levgrnd_class gives the class in which each layer falls. This is relevant for
     ! columns where there are 2 or more fundamentally different layer types. For
     ! example, this distinguishes between soil and bedrock layers. The particular value
     ! assigned to each class is irrelevant; the important thing is that different
     ! classes (e.g., soil vs. bedrock) have different values of levgrnd_class.
     !
     ! levgrnd_class = ispval indicates that the given layer is completely unused for
     ! this column (i.e., this column doesn't use the full nlevgrnd layers).
     integer , pointer :: levgrnd_class        (:,:) ! class in which each layer falls (1:nlevgrnd)
   contains

     procedure, public :: Init
     procedure, public :: Clean

     ! Update the column type for one column. Any updates to col%itype after
     ! initialization should be made via this routine.
     procedure, public :: update_itype

  end type column_type

  type(column_type), public, target :: col !column data structure (soil/snow/canopy columns)
  !------------------------------------------------------------------------

contains
  
  !------------------------------------------------------------------------
  subroutine Init(this, begc, endc)
    !
    ! !ARGUMENTS:
    class(column_type)  :: this
    integer, intent(in) :: begc,endc
    !------------------------------------------------------------------------

    ! The following is set in initGridCellsMod
    allocate(this%gridcell    (begc:endc))                     ; this%gridcell    (:)   = ispval
    allocate(this%wtgcell     (begc:endc))                     ; this%wtgcell     (:)   = nan
    allocate(this%landunit    (begc:endc))                     ; this%landunit    (:)   = ispval
    allocate(this%wtlunit     (begc:endc))                     ; this%wtlunit     (:)   = nan
    allocate(this%patchi      (begc:endc))                     ; this%patchi      (:)   = ispval
    allocate(this%patchf      (begc:endc))                     ; this%patchf      (:)   = ispval
    allocate(this%npatches     (begc:endc))                    ; this%npatches     (:)   = ispval
    allocate(this%itype       (begc:endc))                     ; this%itype       (:)   = ispval
    allocate(this%lun_itype   (begc:endc))                     ; this%lun_itype   (:)   = ispval
    allocate(this%active      (begc:endc))                     ; this%active      (:)   = .false.
    allocate(this%type_is_dynamic(begc:endc))                  ; this%type_is_dynamic(:) = .false.

    allocate(this%is_fates(begc:endc))                         ; this%is_fates(:) = .false.
    
    ! The following is set in initVerticalMod
    allocate(this%snl         (begc:endc))                     ; this%snl         (:)   = ispval  !* cannot be averaged up
    allocate(this%dz          (begc:endc,-nlevsno+1:nlevmaxurbgrnd)) ; this%dz          (:,:) = nan
    allocate(this%z           (begc:endc,-nlevsno+1:nlevmaxurbgrnd)) ; this%z           (:,:) = nan
    allocate(this%zi          (begc:endc,-nlevsno+0:nlevmaxurbgrnd)) ; this%zi          (:,:) = nan
    allocate(this%zii         (begc:endc))                     ; this%zii         (:)   = nan
    allocate(this%lakedepth   (begc:endc))                     ; this%lakedepth   (:)   = spval  
    allocate(this%dz_lake     (begc:endc,nlevlak))             ; this%dz_lake     (:,:) = nan
    allocate(this%z_lake      (begc:endc,nlevlak))             ; this%z_lake      (:,:) = nan
    allocate(this%col_ndx    (begc:endc))                      ; this%col_ndx(:) = ispval
    allocate(this%colu       (begc:endc))                      ; this%colu   (:) = ispval
    allocate(this%cold       (begc:endc))                      ; this%cold   (:) = ispval
    allocate(this%hillslope_ndx(begc:endc))                    ; this%hillslope_ndx (:) = ispval
    allocate(this%hill_elev(begc:endc))                        ; this%hill_elev     (:) = spval
    allocate(this%hill_slope(begc:endc))                       ; this%hill_slope    (:) = spval
    allocate(this%hill_area(begc:endc))                        ; this%hill_area     (:) = spval
    allocate(this%hill_width(begc:endc))                       ; this%hill_width    (:) = spval
    allocate(this%hill_distance(begc:endc))                    ; this%hill_distance (:) = spval
    allocate(this%hill_aspect(begc:endc))                      ; this%hill_aspect (:) = spval
    allocate(this%nbedrock   (begc:endc))                      ; this%nbedrock   (:)   = ispval  
    allocate(this%levgrnd_class(begc:endc,nlevmaxurbgrnd))     ; this%levgrnd_class(:,:) = ispval
    allocate(this%micro_sigma (begc:endc))                     ; this%micro_sigma (:)   = nan
    allocate(this%topo_slope  (begc:endc))                     ; this%topo_slope  (:)   = nan
    allocate(this%topo_std    (begc:endc))                     ; this%topo_std    (:)   = nan
    allocate(this%is_hillslope_column(begc:endc))              ; this%is_hillslope_column(:) = .false.
    allocate(this%hydrologically_active(begc:endc))            ; this%hydrologically_active(:) = .false.
    allocate(this%urbpoi      (begc:endc))                     ; this%urbpoi      (:)   = .false.

  end subroutine Init

  !------------------------------------------------------------------------
  subroutine Clean(this)
    !
    ! !ARGUMENTS:
    class(column_type) :: this
    !------------------------------------------------------------------------

    deallocate(this%gridcell   )
    deallocate(this%wtgcell    )
    deallocate(this%landunit   )
    deallocate(this%wtlunit    )
    deallocate(this%patchi     )
    deallocate(this%patchf     )
    deallocate(this%npatches    )
    deallocate(this%itype      )
    deallocate(this%lun_itype  )
    deallocate(this%active     )
    deallocate(this%is_fates   )
    deallocate(this%type_is_dynamic)
    deallocate(this%snl        )
    deallocate(this%dz         )
    deallocate(this%z          )
    deallocate(this%zi         )
    deallocate(this%zii        )
    deallocate(this%lakedepth  )
    deallocate(this%dz_lake    )
    deallocate(this%z_lake     )
    deallocate(this%micro_sigma)
    deallocate(this%topo_slope )
    deallocate(this%topo_std   )
    deallocate(this%nbedrock   )
    deallocate(this%levgrnd_class)
    deallocate(this%is_hillslope_column)
    deallocate(this%hydrologically_active)
    deallocate(this%col_ndx    )
    deallocate(this%colu       )
    deallocate(this%cold       )
    deallocate(this%hillslope_ndx)
    deallocate(this%hill_elev    )
    deallocate(this%hill_slope   )
    deallocate(this%hill_area    )
    deallocate(this%hill_width   )
    deallocate(this%hill_distance)
    deallocate(this%hill_aspect  )
    deallocate(this%urbpoi       )
  end subroutine Clean

  !-----------------------------------------------------------------------
  subroutine update_itype(this, c, itype)
    !
    ! !DESCRIPTION:
    ! Update the column type for one column. Any updates to col%itype after
    ! initialization should be made via this routine.
    !
    ! This can NOT be used to change the landunit type: it can only be used to change the
    ! column type within a fixed landunit.
    !
    ! !ARGUMENTS:
    class(column_type), intent(inout) :: this
    integer, intent(in) :: c
    integer, intent(in) :: itype
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'update_itype'
    !-----------------------------------------------------------------------

    if (col%type_is_dynamic(c)) then
       col%itype(c) = itype
       col%hydrologically_active(c) = is_hydrologically_active( &
            col_itype = itype, &
            lun_itype = col%lun_itype(c))
       ! Properties that are tied to the landunit's properties (like urbpoi) are assumed
       ! not to change here.
    else
       write(iulog,*) subname//' ERROR: attempt to update itype when type_is_dynamic is false'
       write(iulog,*) 'c, col%itype(c), itype = ', c, col%itype(c), itype
       ! Need to use shr_sys_abort rather than endrun, because using endrun would cause
       ! circular dependencies
       call shr_sys_abort(subname//' ERROR: attempt to update itype when type_is_dynamic is false')
    end if
  end subroutine update_itype



end module ColumnType
