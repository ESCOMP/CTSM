module mkpftMod

  use ESMF
  use pio
  use shr_kind_mod   , only : r8 => shr_kind_r8, r4=>shr_kind_r4
  use shr_sys_mod    , only : shr_sys_abort
  use mkpioMod       , only : mkpio_get_rawdata, mkpio_get_dimlengths
  use mkpioMod       , only : pio_iotype, pio_ioformat, pio_iosystem
  use mkpioMod       , only : mkpio_iodesc_rawdata, mkpio_get_rawdata_level
  use mkesmfMod      , only : regrid_rawdata, create_routehandle_r8, get_meshareas
  use mkutilsMod     , only : chkerr
  use mkvarctl       , only : numpft, root_task, ndiag
  use mkvarpar       , only : numstdpft, numstdcft, noveg
  use mkpftConstantsMod

  implicit none
  private           ! By default make data private

#include <mpif.h>

  public :: mkpftInit          ! Initialization
  public :: mkpft              ! Set PFT
  public :: mkpft_parse_oride  ! Parse the string with PFT fraction/index info to override

  private :: mkpft_check_oride  ! Check the pft_frc and pft_idx values for correctness

  ! When pft_idx and pft_frc are set, they must be set together, and they will cause the
  ! entire area to be covered with vegetation and zero out other landunits.
  ! The sum of pft_frc must = 100%, and each pft_idx point in the array corresponds to
  ! the fraction in pft_frc. Only the first few points are used until pft_frc = 0.0.

  integer :: m                     ! index

  ! PFT vegetation index to override with
  integer, public :: pft_idx(0:maxpft) = (/ ( -1,  m = 0, maxpft ) /)

  ! PFT vegetation fraction to override with
  real(r8), public :: pft_frc(0:maxpft) = (/ ( 0.0_r8, m = 0, maxpft ) /)

  ! private data members:
  logical, public, protected :: use_input_pft = .false. ! Flag to override PFT with input values
  logical, public, protected :: presc_cover   = .false. ! Flag to prescribe vegetation coverage
  integer, private           :: nzero                   ! index of first zero fraction

  type, public :: pft_oride              ! Public only for unit testing
     real(r8) :: crop                    ! Percent covered by crops
     real(r8) :: natveg                  ! Percent covered by natural vegetation
     real(r8), allocatable :: natpft(:)  ! Percent of each natural PFT within the natural veg landunit
     real(r8), allocatable :: cft(:)     ! Percent of each crop CFT within the crop landunit
   contains
     procedure, public :: InitZeroOut      ! Initialize the PFT override object to zero out all vegetation
     procedure, public :: InitAllPFTIndex  ! Initialize the PFT override object with PFT indeces for all veg and crop types
     procedure, public :: Clean            ! Clean up a PFT Override object
  end type pft_oride

  interface pft_oride
     module procedure :: constructor  ! PFT Overide object constructor
  end interface pft_oride

  character(len=35) :: veg(0:maxpft) ! vegetation types

  ! Module instance of PFT override object
  ! Used for both zeroing out PFT's as well as setting specified PFT's over the gridcell
  type(pft_oride), private :: pft_override

  character(len=*) , parameter :: u_FILE_u = &
       __FILE__

!===============================================================
contains
!===============================================================

  subroutine mkpftInit( zero_out_l, all_veg_l )
    !
    ! Initialize of Make PFT data
    !
    ! input/output variables
    logical, intent(in) :: zero_out_l ! If veg should be zero'ed out
    logical, intent(in) :: all_veg_l  ! If should zero out other fractions so that all land-cover is vegetation

    ! local variables:
    logical             :: error_happened    ! If an error was triggered so should return
    real(r8), parameter :: hndrd = 100.0_r8  ! A hundred percent
    character(len=*), parameter :: subname = ' (mkpftInit) '
    !-----------------------------------------------------------------------

    if ( maxpft < numpft ) then
       write(6,*) subname//'number PFT is > max allowed!'
       call shr_sys_abort()
    end if
    nzero = -1
    call mkpft_check_oride( error_happened )
    if ( error_happened )then
       write(6,*) subname//'Problem setting pft override settings'
    end if
    if ( zero_out_l .and. use_input_pft )then
       write(6,*) subname//"trying to both zero out all PFT's as well as set them to specific values"
       call shr_sys_abort()
    end if

    ! If zeroing out, set use_input_pft to true so the pft_override will be used
    if( zero_out_l )then
       nzero = 0
       pft_frc(0) = 0.0_r8
       pft_idx(0) = noveg
       use_input_pft = .true.
    end if
    if ( use_input_pft ) then
       write(6,*) 'Set PFT fraction to : ', pft_frc(0:nzero)
       write(6,*) 'With PFT index      : ', pft_idx(0:nzero)
    end if
    if ( all_veg_l .and. .not. use_input_pft )then
       write(6,*) subname//'if all_veg is set to true then specified PFT indices must be provided (i.e. pft_frc and pft_idx)'
       call shr_sys_abort()
    end if

    if ( zero_out_l .and. all_veg_l )then
       write(6,*) subname//'zeroing out vegetation and setting vegetation to 100% is a contradiction!'
       call shr_sys_abort()
    end if

    ! Determine number of PFTs on the natural vegetation landunit, and number of CFTs on
    ! the crop landunit.
    !
    ! For the sake of dynamic PFTs and dynamic landunits, it helps for the structure of the
    ! surface dataset to reflect the subgrid structure that will be used by CLM. Currently
    ! generic crops will always go on the crop landunit, regardless of whether or not we're
    ! using the extra specific crops (so we always run CLM with create_crop_landunit=.true.).
    ! When we create a surface dataset WITH the extra specific crops, all crops
    ! (including the generic crops) again go on the crop landunit.

    num_natpft = numstdpft - numstdcft
    num_cft    = numpft - num_natpft

    ! Determine array bounds for arrays of just natural pfts and just crops. Note that
    ! these are set up so that they always span 0:numpft, so that there is a 1:1
    ! correspondence between an element in a full 0:numpft array and an element with the
    ! same index in either a natpft array or a cft array.
    natpft_lb = noveg
    natpft_ub = num_natpft
    cft_lb    = num_natpft+1
    cft_ub    = cft_lb + num_cft - 1

    ! Make sure the array indices have been set up properly, to ensure the 1:1
    ! correspondence mentioned above
    if (cft_ub /= numpft) then
       write(6,*) 'CFT_UB set up incorrectly: cft_ub, numpft = ', cft_ub, numpft
       call shr_sys_abort()
    end if

    ! Set the PFT override values if applicable
    pft_override = pft_oride()
    presc_cover = .false.
    if( zero_out_l )then
       call pft_override%InitZeroOut()
       presc_cover = .true.
    else if ( use_input_pft ) then
       call pft_override%InitAllPFTIndex()
       if ( .not. all_veg_l )then
          if ( pft_override%crop <= 0.0 )then
             write(6,*) "Warning: PFT/CFT's are being overridden, but no crop type is being asked for"
          end if
          if ( pft_override%natveg <= 0.0 )then
             write(6,*) "Warning: PFT/CFT's are being overridden, but no natural vegetation type is being asked for"
          end if
          presc_cover = .false.
       else
          presc_cover = .true.
       end if
    end if

    ! -----------------------------------------------------------------
    ! Set the vegetation types
    ! -----------------------------------------------------------------

    if ( numpft >= numstdpft ) then
       veg(0:maxpft) = (/                          &
            'not vegetated                      ', &
            'needleleaf evergreen temperate tree', &
            'needleleaf evergreen boreal tree   ', &
            'needleleaf deciduous boreal tree   ', &
            'broadleaf evergreen tropical tree  ', &
            'broadleaf evergreen temperate tree ', &
            'broadleaf deciduous tropical tree  ', &
            'broadleaf deciduous temperate tree ', &
            'broadleaf deciduous boreal tree    ', &
            'broadleaf evergreen shrub          ', &
            'broadleaf deciduous temperate shrub', &
            'broadleaf deciduous boreal shrub   ', &
            'c3 arctic grass                    ', &
            'c3 non-arctic grass                ', &
            'c4 grass                           ', &
            'c3_crop                            ', &
            'c3_irrigated                       ', &
            'temperate_corn                     ', &
            'irrigated_temperate_corn           ', &
            'spring_wheat                       ', &
            'irrigated_spring_wheat             ', &
            'winter_wheat                       ', &
            'irrigated_winter_wheat             ', &
            'temperate_soybean                  ', &
            'irrigated_temperate_soybean        ', &
            'barley                             ', &
            'irrigated_barley                   ', &
            'winter_barley                      ', &
            'irrigated_winter_barley            ', &
            'rye                                ', &
            'irrigated_rye                      ', &
            'winter_rye                         ', &
            'irrigated_winter_rye               ', &
            'cassava                            ', &
            'irrigated_cassava                  ', &
            'citrus                             ', &
            'irrigated citrus                   ', &
            'cocoa                              ', &
            'irrigated_cocoa                    ', &
            'coffee                             ', &
            'irrigated_coffee                   ', &
            'cotton                             ', &
            'irrigated_cotton                   ', &
            'datepalm                           ', &
            'irrigated_datepalm                 ', &
            'foddergrass                        ', &
            'irrigated_foddergrass              ', &
            'grapes                             ', &
            'irrigated_grapes                   ', &
            'groundnuts                         ', &
            'irrigated_groundnuts               ', &
            'millet                             ', &
            'irrigated_millet                   ', &
            'oilpalm                            ', &
            'irrigated_oilpalm                  ', &
            'potatoes                           ', &
            'irrigated_potatoes                 ', &
            'pulses                             ', &
            'irrigated_pulses                   ', &
            'rapeseed                           ', &
            'irrigated_rapeseed                 ', &
            'rice                               ', &
            'irrigated_rice                     ', &
            'sorghum                            ', &
            'irrigated_sorghum                  ', &
            'sugarbeet                          ', &
            'irrigated_sugarbeet                ', &
            'sugarcane                          ', &
            'irrigated_sugarcane                ', &
            'sunflower                          ', &
            'irrigated_sunflower                ', &
            'miscanthus                         ', &
            'irrigated_miscanthus               ', &
            'switchgrass                        ', &
            'irrigated_switchgrass              ', &
            'tropical_corn                      ', &
            'irrigated_tropical_corn            ', &
            'tropical_soybean                   ', &
            'irrigated_tropical_soybean         ' /)
    end if
    if (numpft == numstdpft )then
       if (root_task) then
          write(ndiag, '(a,i8)')'Creating surface datasets with the standard # of PFTs =', numpft
       end if
    else if ( numpft > numstdpft )then
       if (root_task) then
          write(ndiag,'(a,i8)')'Creating surface datasets with extra types for crops; total pfts =', numpft
       end if
    else
       write(6,*) subname//': parameter numpft is NOT set to a known value (should be 16 or more) =',numpft
       call shr_sys_abort()
    end if

  end subroutine mkpftInit

  !===============================================================
  subroutine mkpft(file_mesh_i, file_data_i, mesh_o, pctlnd_o, pctnatpft_o, pctcft_o, rc)
    !
    ! Make PFT data
    !
    ! This dataset consists of the %cover of the [numpft]+1 PFTs used by
    ! the model. The input %cover pertains to the "vegetated" portion of the
    ! grid cell and sums to 100. The real portion of each grid cell
    ! covered by each PFT is the PFT cover times the fraction of the
    ! grid cell that is land. This is the quantity preserved when
    ! area-averaging from the input grid to the models grid.
    !
    ! Upon return from this routine, the % cover of the natural veg + crop landunits is
    ! generally 100% everywhere; this will be normalized later to account for special landunits.
    !
    use mkpctPftTypeMod,   only : pct_pft_type
    use mkpftConstantsMod, only : natpft_lb, natpft_ub, num_cft, cft_lb, cft_ub
    !
    ! input/output variables
    character(len=*)  , intent(in)    :: file_mesh_i    ! input mesh file name
    character(len=*)  , intent(in)    :: file_data_i    ! input data file name
    type(ESMF_Mesh)   , intent(in)    :: mesh_o         ! model mesh
    real(r8)          , intent(inout) :: pctlnd_o(:)    ! output grid:%land/gridcell
    type(pct_pft_type), intent(inout) :: pctnatpft_o(:) ! natural PFT cover
    type(pct_pft_type), intent(inout) :: pctcft_o(:)    ! crop (CFT) cover
    integer           , intent(out)   :: rc
    !
    ! local variables:
    type(ESMF_RouteHandle)          :: routehandle
    type(ESMF_Mesh)                 :: mesh_i
    type(file_desc_t)               :: pioid
    integer                         :: dimid
    integer                         :: ndims              ! number of dimensions for a variable on the file
    integer                         :: dimlens(3)         ! dimension lengths for a variable on the file
    integer                         :: ns_i, ns_o         ! input/output bounds
    integer                         :: ni,no              ! indices
    integer                         :: k,n,m              ! indices
    type(pct_pft_type), allocatable :: pctnatpft_i(:)     ! input grid: natural PFT cover
    type(pct_pft_type), allocatable :: pctcft_i(:)        ! input grid: crop (CFT) cover
    real(r8), allocatable           :: pct_cft_i(:,:)     ! input grid: CFT (Crop Functional Type) percent (% of landunit cell)
    real(r8), allocatable           :: pct_cft_o(:,:)     ! output grid: CFT (Crop Functional Type) percent (% of landunit cell)
    real(r8), allocatable           :: pct_nat_pft_i(:,:) ! input grid: natural PFT percent (% of landunit cell)
    real(r8), allocatable           :: pct_nat_pft_o(:,:) ! output grid: natural PFT percent (% of landunit cell)
    real(r8), allocatable           :: output_pct_nat_pft_o(:,:)
    real(r8), allocatable           :: output_pct_cft_o(:,:)
    integer , allocatable           :: mask_i(:)
    real(r8), allocatable           :: frac_i(:)
    real(r8), allocatable           :: frac_o(:)
    real(r8), allocatable           :: area_i(:)
    real(r8), allocatable           :: area_o(:)
    real(r8), allocatable           :: pctnatveg_i(:)     ! input grid: natural veg percent (% of grid cell)
    real(r8), allocatable           :: pctnatveg_o(:)     ! output grid: natural veg percent (% of grid cell)
    real(r8), allocatable           :: pctcrop_i(:)       ! input grid: all crop percent (% of grid cell)
    real(r8), allocatable           :: pctcrop_o(:)       ! output grid: all crop percent (% of grid cell)
    real(r8), allocatable           :: pctpft_i(:,:)      ! input grid: PFT percent (for error checks)
    real(r8), allocatable           :: pctpft_o(:,:)      ! output grid: PFT percent (% of grid cell)
    real(r8), allocatable           :: temp_i(:,:)        ! input grid: temporary 2D variable to read in
    integer                         :: numpft_i           ! num of plant types input data
    integer                         :: natpft_i           ! num of natural plant types input data
    integer                         :: ncft_i             ! num of crop types input data
    real(r8)                        :: sum_fldo           ! global sum of dummy output fld
    real(r8)                        :: sum_fldi           ! global sum of dummy input fld
    real(r8)                        :: wst_sum            ! sum of %pft
    real(r8), allocatable           :: gpft_o(:)          ! output grid: global area pfts
    real(r8)                        :: garea_o            ! output grid: global area
    real(r8), allocatable           :: gpft_i(:)          ! input grid: global area pfts
    real(r8)                        :: garea_i            ! input grid: global area
    integer                         :: ier, rcode         ! error status
    real(r8)                        :: relerr = 0.0001_r8 ! max error: sum overlap wts ne 1
    logical                         :: error_happened     ! If an error was triggered so should return
    character(len=*), parameter :: subname = 'mkpf'
    !-----------------------------------------------------------------------

    if (root_task) then
       write(ndiag,*)
       write(ndiag,'(1x,80a1)') ('=',k=1,80)
       write(ndiag,*)
       write(ndiag,'(a)') 'Attempting to make PFTs .....'
       write(ndiag,'(a)') ' Input file is '//trim(file_data_i)
       write(ndiag,'(a)') ' Input mesh file is '//trim(file_mesh_i)
    end if

    ! Open input pft file
    rcode = pio_openfile(pio_iosystem, pioid, pio_iotype, trim(file_data_i), pio_nowrite)
    call ESMF_VMLogMemInfo("After pio_openfile "//trim(file_data_i))

    ! Read in input mesh
    mesh_i = ESMF_MeshCreate(filename=trim(file_mesh_i), fileformat=ESMF_FILEFORMAT_ESMFMESH, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMLogMemInfo("After create mesh_i in "//trim(subname))

    ! Determine ns_i
    call ESMF_MeshGet(mesh_i, numOwnedElements=ns_i, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Determine ns_o
    call ESMF_MeshGet(mesh_o, numOwnedElements=ns_o, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Get the landmask from the file and reset the mesh mask based on that
    allocate(frac_i(ns_i), stat=ier)
    if (ier/=0) call shr_sys_abort()
    allocate(mask_i(ns_i), stat=ier)
    if (ier/=0) call shr_sys_abort()
    call mkpio_get_rawdata(pioid, 'LANDMASK', mesh_i, frac_i, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    do ni = 1,ns_i
       if (frac_i(ni) > 0._r4) then
          mask_i(ni) = 1
       else
          mask_i(ni) = 0
       end if
    end do
    call ESMF_MeshSet(mesh_i, elementMask=mask_i, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! Create Error checks if presc_cover is false
    if ( .not. presc_cover ) then
       rcode = pio_inq_dimid(pioid, 'natpft', dimid)
       rcode = pio_inq_dimlen(pioid, dimid, natpft_i)
       rcode = pio_inq_dimid(pioid, 'cft', dimid)
       rcode = pio_inq_dimlen(pioid, dimid, ncft_i)
       numpft_i = natpft_i + ncft_i

       ! Check if the number of pfts on the input matches the expected number. A mismatch
       ! is okay if the input raw dataset has prognostic crops and the output does not.
       if (numpft_i /= numpft+1) then
          if (numpft_i == numstdpft+1) then
             write(6,*) subname//' ERROR: trying to use non-crop input file'
             write(6,*) 'for a surface dataset with crops.'
             call shr_sys_abort()
          else if (numpft_i > numstdpft+1 .and. numpft_i == maxpft+1) then
             write(6,*) subname//' WARNING: using a crop input raw dataset for a non-crop output surface dataset'
          else
             write(6,*) subname//': parameter numpft+1= ',numpft+1, &
                  'does not equal input dataset numpft= ',numpft_i
             call shr_sys_abort()
          end if
       endif
    end if

    ! ----------------------------------------
    ! Create a route handle between the input and output mesh and get frac_o
    ! ----------------------------------------
    if (.not. presc_cover) then
       allocate(frac_o(ns_o),stat=ier)
       if (ier/=0) call shr_sys_abort()
       call create_routehandle_r8(mesh_i, mesh_o, routehandle, frac_o=frac_o, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call ESMF_VMLogMemInfo("After create routehandle in "//trim(subname))
    end if

    ! ----------------------------------------
    ! Determine pctlnd_o(:) (output argument)
    ! ----------------------------------------
    if ( use_input_pft .and. presc_cover ) then
       pctlnd_o(:) = 100._r8
    else
       pctlnd_o(:) = frac_o(:) * 100._r8
    end if

    ! DEBUG
    RETURN
    !DEBUG

    ! ----------------------------------------
    ! Determine pctnatveg_o(:)
    ! ----------------------------------------
    if ( use_input_pft .and. presc_cover ) then
       do no = 1,ns_o
          pctnatveg_o(no) = pft_override%natveg
       end do
    else
       allocate(pctnatveg_i(ns_i), stat=ier)
       if (ier/=0) call shr_sys_abort()
       call mkpio_get_rawdata(pioid, 'PCT_NATVEG', mesh_i, pctnatveg_i, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call regrid_rawdata(mesh_i, mesh_o, routehandle, pctnatveg_i, pctnatveg_o, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    ! ----------------------------------------
    ! Determine pctcrop_o(:)
    ! ----------------------------------------
    if ( use_input_pft .and. presc_cover ) then
       do no = 1,ns_o
          pctcrop_o(no) = pft_override%crop
       end do
    else
       allocate(pctcrop_i(ns_i), stat=ier)
       if (ier/=0) call shr_sys_abort()
       call mkpio_get_rawdata(pioid, 'PCT_CROP', mesh_i, pctcrop_i, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call regrid_rawdata(mesh_i, mesh_o, routehandle, pctcrop_i, pctcrop_o, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    ! ----------------------------------------
    ! Determine pct_nat_pft_o(:,:)
    ! ----------------------------------------
    if ( .not. use_input_pft )then
       allocate(pct_nat_pft_i(0:num_natpft,ns_i))
       if (ier/=0) call shr_sys_abort()
       allocate(pct_nat_pft_o(0:num_natpft,ns_o))
       if (ier/=0) call shr_sys_abort()
       call mkpio_get_rawdata(pioid, 'PCT_NAT_PFT', mesh_i, pct_nat_pft_i, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       do m = 0,num_natpft
          do ni = 1,ns_i
             pct_nat_pft_i(m,ni) = pct_nat_pft_i(m,ni) * pctnatveg_i(ni) * 0.01_r8 * mask_i(ni)
          end do
       end do
       call regrid_rawdata(mesh_i, mesh_o, routehandle, pct_nat_pft_i, pct_nat_pft_o, 0, num_natpft, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       do m = 0,num_natpft
          do no = 1,ns_o
             pct_nat_pft_o(m,no) = pct_nat_pft_o(m,no) / pctnatveg_o(no) * 0.01_r8
             if (pctlnd_o(no) < 1.0e-6 .or. pctnatveg_o(no) < 1.0e-6) then
                if (m == 0) then
                   pct_nat_pft_o(m,no) = 100._r8
                else
                   pct_nat_pft_o(m,no) = 0._r8
                endif
             end if
          end do
       end do
    end if

    ! ----------------------------------------
    ! Determine pct_cft_o(:,:)
    ! ----------------------------------------
    allocate(pct_cft_i(1:num_cft,ns_i), stat=ier)
    if (ier/=0) call shr_sys_abort()
    allocate(pct_cft_o(1:num_cft,ns_o), stat=ier)
    if (ier/=0) call shr_sys_abort()

    if ( .not. use_input_pft )then
       call mkpio_get_dimlengths(pioid, 'PCT_CFT', ndims, dimlens(:))
       if (ndims == 3 .and. dimlens(1)*dimlens(2) == ns_i .and. dimlens(3) == num_cft )then
          call mkpio_get_rawdata(pioid, 'PCT_CFT', mesh_i, pct_cft_i, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
       else if ( ndims == 3 .and. dimlens(1)*dimlens(2) == ns_i .and. dimlens(3) > num_cft )then
          ! Read in the whole array: then sum the rainfed and irrigated seperately
          allocate(temp_i(dimlens(3),ns_i))
          call mkpio_get_rawdata(pioid, 'PCT_CFT', mesh_i, temp_i, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          do n = 1, num_cft
             pct_cft_i(n,:) = 0.0_r8
             do m = n, dimlens(3), 2
                pct_cft_i(n,:) = pct_cft_i(n,:) + temp_i(m,:)
             end do
          end do
          deallocate(temp_i)
       else
          call shr_sys_abort(subname//' error: dimensions for PCT_CROP are NOT what is expected')
       end if

       do m = 0,num_cft
          do ni = 1,ns_i
             pct_cft_i(m,ni) = pct_cft_i(m,ni) * pctcrop_i(ni) * 0.01_r8 * mask_i(ni)
          end do
       end do
       call regrid_rawdata(mesh_i, mesh_o, routehandle, pct_cft_i, pct_cft_o, 1, num_cft, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       do no = 1,ns_o
          pct_cft_o(m,no) = pct_cft_o(m,no) / pctcrop_o(no) * 0.01_r8
          if (pctlnd_o(no) < 1.0e-6 .or. pctcrop_o(no) < 1.0e-6) then
             if (m == 1) then
                pct_cft_o(no,m) = 100._r8
             else
                pct_cft_o(no,m) = 0._r8
             endif
          end if
       enddo
    end if

    ! ----------------------------------------
    ! If specific PFT/CFT's are prescribed set them directly
    ! But first do error checking to make sure specific veg
    ! types are given where nat-veg and crop is assigned
    ! ----------------------------------------
    if (use_input_pft) then
       do no = 1,ns_o
          if (pctlnd_o(no) > 1.0e-6 .and. pctnatveg_o(no) > 1.0e-6) then
             if ( pft_override%natveg <= 0.0_r8 )then
                write(6,*) subname//': ERROR: no natural vegetation PFTs are being prescribed '//&
                     ' but there are natural vegetation areas: provide at least one natural veg PFT'
                call shr_sys_abort()
             end if
          end if
          if (pctlnd_o(no) > 1.0e-6 .and. pctcrop_o(no) > 1.0e-6) then
             if ( pft_override%crop <= 0.0_r8 )then
                write(6,*) subname//': ERROR: no crop CFTs are being prescribed '&
                     //' but there are crop areas: provide at least one CFT'
                call shr_sys_abort()
             end if
          end if
       end do

       do no = 1,ns_o
          if (pctlnd_o(no) > 1.0e-6 .and. pctnatveg_o(no) > 1.0e-6) then
             pct_nat_pft_o(noveg:num_natpft,no) = pft_override%natpft(noveg:num_natpft)
          else
             pct_nat_pft_o(noveg,no)    = 100._r8
             pct_nat_pft_o(noveg+1:,no) = 0._r8
          end if
          if (pctlnd_o(no) > 1.0e-6 .and. pctcrop_o(no) > 1.0e-6) then
             pct_cft_o(1:num_cft,no) = pft_override%cft(1:num_cft)
          else
             pct_cft_o(1,no)  = 100._r8
             pct_cft_o(2:,no) = 0._r8
          end if
          pctpft_o(natpft_lb:natpft_ub,no)   = pct_nat_pft_o(0:num_natpft,no)
          pctpft_o(cft_lb:cft_ub,no)         = pct_cft_o(1:num_cft,no)
       end do
    end if

    ! ----------------------------------------
    ! Correct sums so that if they differ slightly from 100, they are
    ! corrected to equal 100 more exactly.
    ! ----------------------------------------
    do no = 1,ns_o
       wst_sum = 0.
       do m = 0, num_natpft
          wst_sum = wst_sum + pct_nat_pft_o(m,no)
       enddo
       ! Error check: percents should sum to 100 for land grid cells, within roundoff
       if (abs(wst_sum-100._r8) > relerr) then
          write (6,*) subname//'error: nat pft = ', (pct_nat_pft_o(m,no), m = 0, num_natpft), &
               ' do not sum to 100. at no = ',no,' but to ', wst_sum
          call shr_sys_abort()
       end if

       ! Correct sum so that if it differs slightly from 100, it is corrected to equal
       ! 100 more exactly
       do m = 1, num_natpft
          pct_nat_pft_o(m,no) = pct_nat_pft_o(m,no) * 100._r8 / wst_sum
       end do

       wst_sum = 0.
       do m = 1, num_cft
          wst_sum = wst_sum + pct_cft_o(m,no)
       enddo
       ! Error check: percents should sum to 100 for land grid cells, within roundoff
       if (abs(wst_sum-100._r8) > relerr) then
          write (6,*) subname//'error: crop cft = ',(pct_cft_o(no,m), m = 1, num_cft), &
               ' do not sum to 100. at no = ',no,' but to ', wst_sum
          call shr_sys_abort()
       end if

       ! Correct sum so that if it differs slightly from 100, it is corrected to equal
       ! 100 more exactly
       do m = 1, num_cft
          pct_cft_o(m,no) = pct_cft_o(m,no) * 100._r8 / wst_sum
       end do
    end do

    ! ----------------------------------------
    ! Convert % pft as % of grid cell to % pft on the landunit and % of landunit on the grid cell
    ! ----------------------------------------
    do no = 1,ns_o
       output_pct_nat_pft_o(no,:) = pct_nat_pft_o(:,no)
       pctnatpft_o(no) = pct_pft_type( output_pct_nat_pft_o(no,:), pctnatveg_o(no), first_pft_index=natpft_lb )

       output_pct_cft_o(no,:) = pct_cft_o(:,no)
       pctcft_o(no)    = pct_pft_type( output_pct_cft_o(no,:),     pctcrop_o(no),   first_pft_index=cft_lb    )
    end do

    ! -----------------------------------------------------------------
    ! Error check
    ! Compare global areas on input and output grids
    ! Only when you aren't prescribing the vegetation coverage everywhere
    ! If use_input_pft is set this will compare the global coverage of
    ! the prescribed vegetation to the coverage of PFT/CFT's on the input
    ! datasets.
    ! -----------------------------------------------------------------

    !     if (presc_cover) then
    !        ns_i = 1
    !        numpft_i = numpft+1
    !     end if
    
    !     if ( .not. presc_cover ) then
    !        ! Derived types
    !        allocate(pctnatpft_i(ns_i), stat=ier)
    !        if (ier/=0) call shr_sys_abort()
    !        allocate(pctcft_i(ns_i), stat=ier)
    !        if (ier/=0) call shr_sys_abort()
    
    !        do ni = 1, ns_i
    !           pctnatpft_i(ni) = pct_pft_type( pct_nat_pft_i(ni,:), pctnatveg_i(ni), first_pft_index=natpft_lb )
    !           pctcft_i(ni)    = pct_pft_type( pct_cft_i(ni,:),     pctcrop_i(ni),   first_pft_index=cft_lb    )
    !        end do
    
    !        do no = 1,ns_o
    !           pctpft_o(no,natpft_lb:natpft_ub) = pctnatpft_o(no)%get_pct_p2g()
    !           pctpft_o(no,cft_lb:cft_ub)       = pctcft_o(no)%get_pct_p2g()
    !        end do
    
    !        allocate(gpft_i(0:numpft_i-1))
    !        allocate(gpft_o(0:numpft_i-1))
    
    !        ! input grid
    !        allocate(pctpft_i(ns_i,0:(numpft_i-1)), stat=ier)
    !        if (ier/=0) call shr_sys_abort()
    !        allocate(pctpft_o(ns_o,0:(numpft_i-1)), stat=ier)
    !        if (ier/=0) call shr_sys_abort()
    
    !        gpft_i(:) = 0.
    !        garea_i   = 0.
    !        do ni = 1,ns_i
    !           garea_i = garea_i + area_src(ni)*re**2
    !           do m = 0, numpft_i - 1
    !              gpft_i(m) = gpft_i(m) + pctpft_i(ni,m) * area_src(ni) * mask(ni)*re**2
    !           end do
    !        end do
    !        if ( allocated(pctpft_i) ) deallocate (pctpft_i)
    
    !        ! output grid
    !        gpft_o(:) = 0.
    !        garea_o   = 0.
    !        do no = 1,ns_o
    !           garea_o = garea_o + area_dst(no)*re**2
    !           do m = 0, numpft_i - 1
    !              gpft_o(m) = gpft_o(m) + pctpft_o(no,m)*area_dst(no)*frac_o(no)*re**2
    !           end do
    !        end do
    
    !        ! comparison
    !        write (ndiag,*)
    !        write (ndiag,'(1x,70a1)') ('=',k=1,70)
    !        write (ndiag,*) 'PFTs Output'
    !        write (ndiag,'(1x,70a1)') ('=',k=1,70)
    
    !        write (ndiag,*)
    !        write (ndiag,'(1x,70a1)') ('.',k=1,70)
    !        write (ndiag,1001)
    ! 1001   format (1x,'plant type     ',20x,' input grid area',' output grid area',/ &
    !             1x,33x,'     10**6 km**2','      10**6 km**2')
    !        write (ndiag,'(1x,70a1)') ('.',k=1,70)
    !        write (ndiag,*)
    !        do m = 0, numpft_i - 1
    !           write (ndiag,1002) veg(m), gpft_i(m)*1.e-06/100.,gpft_o(m)*1.e-06/100.
    !        end do
    ! 1002   format (1x,a35,f16.3,f17.3)
    
    !        deallocate(gpft_i, gpft_o, frac_dst)
    
    !     end if
    
    ! Clean up memory
    call ESMF_RouteHandleDestroy(routehandle, nogarbage = .true., rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort()
    call ESMF_MeshDestroy(mesh_i, nogarbage = .true., rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort()

    if (root_task) then
       write (ndiag,'(a)') 'Successfully made PFTs'
       write (ndiag,*)
    end if

  end subroutine mkpft

  !-----------------------------------------------------------------------
  subroutine mkpft_parse_oride( string )
    !
    ! Parse the string with pft fraction and index information on it, to override
    ! the file with this information rather than reading from a file.
    !

    use shr_string_mod, only: shr_string_betweenTags, shr_string_countChar

    ! !ARGUMENTS:
    character(len=256), intent(IN) :: string  ! String to parse with PFT fraction  and index data
    !
    ! !LOCAL VARIABLES:
    integer :: rc                         ! error return code
    integer :: num_elms                   ! number of elements
    character(len=256) :: substring       ! string between tags
    character(len=*), parameter :: frc_start = "<pft_f>"
    character(len=*), parameter :: frc_end   = "</pft_f>"
    character(len=*), parameter :: idx_start = "<pft_i>"
    character(len=*), parameter :: idx_end   = "</pft_i>"
    character(len=*), parameter :: subname = 'mkpft_parse_oride'
    !-----------------------------------------------------------------------

    ! NOTE(bja, 2015-02) pft_frc and pft_index can be reset multiple
    ! times by calls to this function. If the number of elements being
    ! set is different each time, then we are working with out of date
    ! information, and the sums may not sum to 100%.
    pft_frc = 0.0_r8
    pft_idx = -1

    call shr_string_betweenTags( string, frc_start, frc_end, substring, rc )
    if ( rc /= 0 )then
       write(6,*) subname//'Trouble finding pft_frac start end tags'
       call shr_sys_abort()
       return
    end if
    num_elms = shr_string_countChar( substring, ",", rc )
    read(substring,*) pft_frc(0:num_elms)
    call shr_string_betweenTags( string, idx_start, idx_end, substring, rc )
    if ( rc /= 0 )then
       write(6,*) subname//'Trouble finding pft_index start end tags'
       call shr_sys_abort()
       return
    end if
    if ( num_elms /= shr_string_countChar( substring, ",", rc ) )then
       write(6,*) subname//'number of elements different between frc and idx fields'
       call shr_sys_abort()
       return
    end if
    read(substring,*) pft_idx(0:num_elms)

  end subroutine mkpft_parse_oride

  !-----------------------------------------------------------------------
  subroutine  mkpft_check_oride( error_happened )
    !
    ! Check that the pft override values are valid
    !
    ! input/output variables
    logical, intent(out) :: error_happened ! Result, true if there was a problem
    !
    ! local variables:
    integer  :: i, j                         ! indices
    real(r8) :: sumpft                       ! Sum of pft_frc
    real(r8), parameter :: hndrd = 100.0_r8  ! A hundred percent
    character(len=*), parameter :: subname = 'mkpftMod::mkpft_check_oride() '
    !-----------------------------------------------------------------------

    error_happened = .false.
    sumpft = sum(pft_frc)
    if (          sumpft == 0.0 )then
       ! PFT fraction is NOT used
       use_input_pft = .false.
    else if ( abs(sumpft - hndrd) > 1.e-6 )then
       write(6, '(a, a, f15.12)') subname, 'Sum of PFT fraction is NOT equal to 100% =', sumpft
       write(6,*) 'Set PFT fraction to : ', pft_frc(0:nzero)
       write(6,*) 'With PFT index      : ', pft_idx(0:nzero)
       error_happened = .true.
       call shr_sys_abort()
    else
       use_input_pft = .true.
       nzero = numpft
       do i = 0, numpft
          if ( pft_frc(i) == 0.0_r8 )then
             nzero = i-1
             exit
          end if
       end do
       ! PFT fraction IS used, and sum is OK, now check details
       do i = 0, nzero
          if ( pft_frc(i) < 0.0_r8 .or. pft_frc(i) > hndrd )then
             write(6,*) subname//'PFT fraction is out of range: pft_frc=', pft_frc(i)
             error_happened = .true.
             call shr_sys_abort()
          else if ( pft_frc(i) > 0.0_r8 .and. pft_idx(i) == -1 )then
             write(6,*) subname//'PFT fraction > zero, but index NOT set: pft_idx=', pft_idx(i)
             error_happened = .true.
             call shr_sys_abort()
          end if
          ! PFT index out of range
          if ( pft_idx(i) < 0 .or. pft_idx(i) > numpft )then
             write(6,*) subname//'PFT index is out of range: ', pft_idx(i)
             error_happened = .true.
             call shr_sys_abort()
          end if
          ! Make sure index values NOT used twice
          do j = 0, i-1
             if ( pft_idx(i) == pft_idx(j) )then
                write(6,*) subname//'Same PFT index is used twice: ', pft_idx(i)
                error_happened = .true.
                call shr_sys_abort()
             end if
          end do
       end do
       ! Make sure the rest of the fraction is zero and index are not set as well
       do i = nzero+1, numpft
          if ( pft_frc(i) /= 0.0_r8 .or. pft_idx(i) /= -1 )then
             write(6,*) subname//'After PFT fraction is zeroed out, fraction is non zero, or index set'
             error_happened = .true.
             call shr_sys_abort()
          end if
       end do
    end if

  end subroutine mkpft_check_oride

  !-----------------------------------------------------------------------
  function constructor( ) result(this)
    !
    ! Construct a new PFT override object
    !
    ! input/output variables
    type(pft_oride) :: this
    character(len=*), parameter :: subname = 'mkpftMod::constructor() '

    this%crop   = -1.0_r8
    this%natveg = -1.0_r8
    if ( num_natpft < 0 )then
       write(6,*) subname//'num_natpft is NOT set = ', num_natpft
       call shr_sys_abort()
       return
    end if
    if ( num_cft < 0 )then
       write(6,*) subname//'num_cft is NOT set = ', num_cft
       call shr_sys_abort()
       return
    end if
    allocate( this%natpft(noveg:num_natpft) )
    allocate( this%cft(1:num_cft) )
    this%natpft(:) = -1.0_r8
    this%cft(:)    = -1.0_r8
    call this%InitZeroOut()

  end function constructor

  !-----------------------------------------------------------------------
  subroutine InitZeroOut( this )
    !
    ! Initialize a pft_oride object with vegetation that's zeroed out
    !
    ! input/output variables
    class(pft_oride), intent(inout) :: this

    this%crop          = 0.0_r8
    this%natveg        = 0.0_r8

    this%natpft        = 0.0_r8
    this%natpft(noveg) = 100.0_r8
    this%cft           = 0.0_r8
    this%cft(1)        = 100.0_r8

  end subroutine InitZeroOut

  !-----------------------------------------------------------------------
  subroutine InitAllPFTIndex( this )
    !
    ! Initialize a pft_oride object with vegetation that's zeroed out
    !
    ! input/otuput variables
    class(pft_oride), intent(inout) :: this

    ! local variables
    integer :: m, i                  ! Indices
    real(r8) :: croptot              ! Total of crop
    real(r8) :: natvegtot            ! Total of natural vegetation
    character(len=*), parameter :: subname = 'mkpftMod::coInitAllPFTIndex() '

    croptot     = 0.0_r8
    natvegtot   = 0.0_r8
    this%natpft = 0.0_r8
    this%cft    = 0.0_r8
    do m = noveg, nzero
       i = pft_idx(m)
       if ( (i < noveg) .or. (i > numpft) )then
          write(6,*)  subname//'PFT index is out of valid range'
          call shr_sys_abort()
          return
       else if ( i <= num_natpft )then
          this%natpft(i) = pft_frc(m)
          natvegtot = natvegtot + pft_frc(m)
       else
          this%cft(i-num_natpft) = pft_frc(m)
          croptot = croptot + pft_frc(m)
       end if
    end do
    this%crop   = croptot
    this%natveg = natvegtot
    ! Renormalize
    if ( natvegtot > 0.0_r8 )then
       this%natpft = 100.0_r8 * this%natpft / natvegtot
    else
       this%natpft(noveg) = 100.0_r8
    end if
    if (croptot > 0.0_r8 )then
       this%cft = 100.0_r8 * this%cft / croptot
    else
       this%cft(1) = 100.0_r8
    end if

  end subroutine InitAllPFTIndex

  !-----------------------------------------------------------------------
  subroutine Clean( this )
    !
    ! Clean up a PFT Oride object
    !
    ! !ARGUMENTS:
    class(pft_oride), intent(inout) :: this

    this%crop   = -1.0_r8
    this%natveg = -1.0_r8
    deallocate( this%natpft )
    deallocate( this%cft    )

  end subroutine Clean

end module mkpftMod
