module PatchType

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Patch data type allocation 
  ! -------------------------------------------------------- 
  ! patch types can have values of
  ! -------------------------------------------------------- 
  !   0  => not_vegetated
  !   1  => needleleaf_evergreen_temperate_tree
  !   2  => needleleaf_evergreen_boreal_tree
  !   3  => needleleaf_deciduous_boreal_tree
  !   4  => broadleaf_evergreen_tropical_tree
  !   5  => broadleaf_evergreen_temperate_tree
  !   6  => broadleaf_deciduous_tropical_tree
  !   7  => broadleaf_deciduous_temperate_tree
  !   8  => broadleaf_deciduous_boreal_tree
  !   9  => broadleaf_evergreen_shrub
  !   10 => broadleaf_deciduous_temperate_shrub
  !   11 => broadleaf_deciduous_boreal_shrub
  !   12 => c3_arctic_grass
  !   13 => c3_non-arctic_grass
  !   14 => c4_grass
  !   15 => c3_crop
  !   16 => c3_irrigated
  !   17 => temperate_corn
  !   18 => irrigated_temperate_corn
  !   19 => spring_wheat
  !   20 => irrigated_spring_wheat
  !   21 => winter_wheat
  !   22 => irrigated_winter_wheat
  !   23 => temperate_soybean
  !   24 => irrigated_temperate_soybean
  !   25 => barley
  !   26 => irrigated_barley
  !   27 => winter_barley
  !   28 => irrigated_winter_barley
  !   29 => rye
  !   30 => irrigated_rye
  !   31 => winter_rye
  !   32 => irrigated_winter_rye
  !   33 => cassava
  !   34 => irrigated_cassava
  !   35 => citrus
  !   36 => irrigated_citrus
  !   37 => cocoa
  !   38 => irrigated_cocoa
  !   39 => coffee
  !   40 => irrigated_coffee
  !   41 => cotton
  !   42 => irrigated_cotton
  !   43 => datepalm
  !   44 => irrigated_datepalm
  !   45 => foddergrass
  !   46 => irrigated_foddergrass
  !   47 => grapes
  !   48 => irrigated_grapes
  !   49 => groundnuts
  !   50 => irrigated_groundnuts
  !   51 => millet
  !   52 => irrigated_millet
  !   53 => oilpalm
  !   54 => irrigated_oilpalm
  !   55 => potatoes
  !   56 => irrigated_potatoes
  !   57 => pulses
  !   58 => irrigated_pulses
  !   59 => rapeseed
  !   60 => irrigated_rapeseed
  !   61 => rice
  !   62 => irrigated_rice
  !   63 => sorghum
  !   64 => irrigated_sorghum
  !   65 => sugarbeet
  !   66 => irrigated_sugarbeet
  !   67 => sugarcane
  !   68 => irrigated_sugarcane
  !   69 => sunflower
  !   70 => irrigated_sunflower
  !   71 => miscanthus
  !   72 => irrigated_miscanthus
  !   73 => switchgrass
  !   74 => irrigated_switchgrass
  !   75 => tropical_corn
  !   76 => irrigated_tropical_corn
  !   77 => tropical_soybean
  !   78 => irrigated_tropical_soybean
  ! -------------------------------------------------------- 
  !
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
  use clm_varcon     , only : ispval
  use clm_varctl     , only : use_ed 
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  private
  !
  type, public :: patch_type

     ! g/l/c/p hierarchy, local g/l/c/p cells only
     integer , pointer :: column   (:) ! index into column level quantities
     real(r8), pointer :: wtcol    (:) ! weight (relative to column) 
     integer , pointer :: landunit (:) ! index into landunit level quantities
     real(r8), pointer :: wtlunit  (:) ! weight (relative to landunit) 
     integer , pointer :: gridcell (:) ! index into gridcell level quantities
     real(r8), pointer :: wtgcell  (:) ! weight (relative to gridcell) 

     ! Non-ED only 
     integer , pointer :: itype    (:) ! patch vegetation 
     integer , pointer :: mxy      (:) ! m index for laixy(i,j,m),etc. (undefined for special landunits)
     logical , pointer :: active   (:) ! true=>do computations on this patch

     ! ED only
     logical , pointer :: is_veg   (:)
     logical , pointer :: is_bareground  (:)
     real(r8), pointer :: wt_ed       (:) !TODO mv ? can this be removed

   contains

     procedure, public :: Init
     procedure, public :: Clean
     
  end type patch_type
  type(patch_type), public, target :: patch  ! patch type data structure
  !------------------------------------------------------------------------

contains
  
  !------------------------------------------------------------------------
  subroutine Init(this, begp, endp)
    !
    ! !ARGUMENTS:
    class(patch_type)   :: this
    integer, intent(in) :: begp,endp
    !
    ! LOCAL VARAIBLES:
    !------------------------------------------------------------------------

    ! The following is set in InitGridCells

    allocate(this%gridcell      (begp:endp)); this%gridcell   (:) = ispval
    allocate(this%wtgcell       (begp:endp)); this%wtgcell    (:) = nan

    allocate(this%landunit      (begp:endp)); this%landunit   (:) = ispval
    allocate(this%wtlunit       (begp:endp)); this%wtlunit    (:) = nan

    allocate(this%column        (begp:endp)); this%column     (:) = ispval
    allocate(this%wtcol         (begp:endp)); this%wtcol      (:) = nan

    allocate(this%mxy           (begp:endp)); this%mxy        (:) = ispval
    allocate(this%active        (begp:endp)); this%active     (:) = .false.

    ! TODO (MV, 10-17-14): The following must be commented out because
    ! currently the logic checking if patch%itype(p) is not equal to noveg
    ! is used in RootBiogeophysMod in zeng2001_rootfr- a filter is not used
    ! in that routine - which would elimate this problem

    !    if (.not. use_ed) then
    allocate(this%itype      (begp:endp)); this%itype      (:) = ispval
    !    end if

    if (use_ed) then
       allocate(this%is_veg  (begp:endp)); this%is_veg  (:) = .false.
       allocate(this%is_bareground (begp:endp)); this%is_bareground (:) = .false.
       allocate(this%wt_ed      (begp:endp)); this%wt_ed      (:) = nan 
    end if

  end subroutine Init

  !------------------------------------------------------------------------
  subroutine Clean(this)
    !
    ! !ARGUMENTS:
    class(patch_type) :: this
    !------------------------------------------------------------------------

    deallocate(this%gridcell)
    deallocate(this%wtgcell )
    deallocate(this%landunit)
    deallocate(this%wtlunit )
    deallocate(this%column  )
    deallocate(this%wtcol   )
    deallocate(this%itype   )
    deallocate(this%mxy     )
    deallocate(this%active  )
    
    if (use_ed) then
       deallocate(this%is_veg)
       deallocate(this%is_bareground)
       deallocate(this%wt_ed)
    end if

  end subroutine Clean

end module PatchType
