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
  use mkvarctl       , only : numpft, root_task, ndiag, mpicom
  use mkvarpar       , only : numstdpft, numstdcft, noveg
  use mkpftConstantsMod

  implicit none
  private           ! By default make data private

#include <mpif.h>

  public :: mkpftInit ! Initialization
  public :: mkpft     ! Set PFT

  integer :: m ! index

  character(len=35) :: veg(0:maxpft) ! vegetation types
  real(r8), allocatable :: frac_o_nonorm(:)
  type(ESMF_RouteHandle) :: routehandle_nonorm

  character(len=*) , parameter :: u_FILE_u = &
       __FILE__

!===============================================================
contains
!===============================================================

  subroutine mkpftInit( )
    !
    ! Initialize of PFT data
    !
    ! local variables:
    character(len=*), parameter :: subname = ' (mkpftInit) '
    !-----------------------------------------------------------------------

    if ( maxpft < numpft ) then
       write(6,*) subname//'number PFT is > max allowed!'
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

    if (root_task) then
       write(ndiag, '(a, i8)') subname//' num_natpft = ',num_natpft
       write(ndiag, '(a, i8)') subname//' natpft_lb  = ',natpft_lb
       write(ndiag, '(a, i8)') subname//' natpft_ub  = ',natpft_ub
       write(ndiag, '(a, i8)') subname//' num_cft    = ',num_cft
       write(ndiag, '(a, i8)') subname//' cft_lb     = ',cft_lb
       write(ndiag, '(a, i8)') subname//' cft_ub     = ',cft_ub
    end if

    ! Make sure the array indices have been set up properly, to ensure the 1:1
    ! correspondence mentioned above
    if (cft_ub /= numpft) then
       write(6,*) 'CFT_UB set up incorrectly: cft_ub, numpft = ', cft_ub, numpft
       call shr_sys_abort()
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
  subroutine mkpft(file_mesh_i, file_data_i, mesh_o, pctlnd_o, pctnatpft_o, &
                   pctcft_o, rc)
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
    use mkinputMod,        only : mksrf_fdynuse
    !
    ! input/output variables
    character(len=*)  , intent(in)    :: file_mesh_i    ! input mesh file name
    character(len=*)  , intent(in)    :: file_data_i    ! input data file name
    type(ESMF_Mesh)   , intent(in)    :: mesh_o         ! model mesh
    real(r8)          , intent(out)   :: pctlnd_o(:)    ! output grid:%land/gridcell
    type(pct_pft_type), intent(inout) :: pctnatpft_o(:) ! natural PFT cover
    type(pct_pft_type), intent(inout) :: pctcft_o(:)    ! crop (CFT) cover

    integer           , intent(out)   :: rc
    !
    ! local variables:
    type(ESMF_Mesh)                 :: mesh_i
    type(file_desc_t)               :: pioid
    integer                         :: dimid
    integer                         :: ndims              ! number of dimensions for a variable on the file
    integer, allocatable            :: dimlens(:)         ! dimension lengths for a variable on the file
    integer                         :: ns_i, ns_o         ! input/output bounds
    integer                         :: ni,no              ! indices
    integer                         :: k,n,m              ! indices
    type(pct_pft_type), allocatable :: pctnatpft_i(:)     ! input natural PFT cover
    type(pct_pft_type), allocatable :: pctcft_i(:)        ! input crop (CFT) cover
    real(r8), allocatable           :: pct_cft_i(:,:)     ! input  CFT (Crop Functional Type) percent (% of landunit cell)
    real(r8), allocatable           :: pct_cft_o(:,:)     ! output CFT (Crop Functional Type) percent (% of landunit cell)
    real(r8), allocatable           :: pct_nat_pft_i(:,:) ! input  natural PFT percent (% of landunit cell)
    real(r8), allocatable           :: pct_nat_pft_o(:,:) ! output natural PFT percent (% of landunit cell)
    real(r8), allocatable           :: output_pct_nat_pft_o(:,:)
    real(r8), allocatable           :: output_pct_cft_o(:,:)
    integer , allocatable           :: mask_i(:)
    real(r8), allocatable           :: frac_i(:)
    real(r8), allocatable           :: pctlnd_i(:)        ! input  land fraction
    real(r8), allocatable           :: pctnatveg_i(:)     ! input  natural veg percent (% of grid cell)
    real(r8), allocatable           :: pctnatveg_o(:)     ! output natural veg percent (% of grid cell)
    real(r8), allocatable           :: pctcrop_i(:)       ! input  all crop percent (% of grid cell)
    real(r8), allocatable           :: pctcrop_o(:)       ! output all crop percent (% of grid cell)
    real(r8), allocatable           :: pctpft_i(:,:)      ! input  PFT percent (for error checks)
    real(r8), allocatable           :: pctpft_o(:,:)      ! output PFT percent (% of grid cell)
    real(r8), allocatable           :: temp_i(:,:)        ! input  temporary 2D variable to read in
    integer                         :: numpft_i           ! num of plant types input data
    integer                         :: natpft_i           ! num of natural plant types input data
    integer                         :: ncft_i             ! num of crop types input data
    real(r8)                        :: wst_sum            ! sum of %pft
    real(r8), allocatable           :: area_o(:)
    real(r8), allocatable           :: loc_gpft_o(:)     ! output global area pfts
    real(r8), allocatable           :: glob_gpft_o(:)     ! output global area pfts
    integer                         :: ier, rcode         ! error status
    real(r8)                        :: relerr = 0.0001_r8 ! max error: sum overlap wts ne 1
    character(len=*), parameter :: subname = 'mkpft'
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

    ! Create error checks
    rcode = pio_inq_dimid(pioid, 'natpft', dimid)
    rcode = pio_inq_dimlen(pioid, dimid, natpft_i)
    rcode = pio_inq_dimid(pioid, 'cft', dimid)
    rcode = pio_inq_dimlen(pioid, dimid, ncft_i)
    numpft_i = natpft_i + ncft_i

    ! Check if the number of pfts on the input matches the expected number. A mismatch
    ! is okay if the input raw dataset has prognostic crops and the output does not.
    if (numpft_i /= numpft+1) then
       if (numpft_i == numstdpft+1) then
          if (root_task) then
             write(ndiag,*) subname//' ERROR: trying to use non-crop input file'
             write(ndiag,*) 'for a surface dataset with crops.'
          end if
          call shr_sys_abort()
       else if (numpft_i > numstdpft+1 .and. numpft_i == maxpft+1) then
          if (root_task) then
             write(ndiag,*) subname//' WARNING: using a crop input raw dataset for a non-crop output surface dataset'
          end if
       else
          if (root_task) then
             write(ndiag,*) subname//': parameter numpft+1= ',numpft+1, &
                  'does not equal input dataset numpft= ',numpft_i
          end if
          call shr_sys_abort()
       end if
    endif

    ! ----------------------------------------
    ! Create a route handle between the input and output mesh and get frac_o_nonorm
    ! ----------------------------------------
    if (.not. ESMF_RouteHandleIsCreated(routehandle_nonorm)) then
       allocate(frac_o_nonorm(ns_o),stat=ier)
       if (ier/=0) call shr_sys_abort()
       ! Note that norm_by_fracs is false in the following because this routehandle is
       ! used to map fields that are expressed in terms of % of the grid cell.
       call create_routehandle_r8(mesh_i=mesh_i, mesh_o=mesh_o, norm_by_fracs=.false., &
            routehandle=routehandle_nonorm, frac_o=frac_o_nonorm, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call ESMF_VMLogMemInfo("After create routehandle in "//trim(subname))
    end if

    ! ----------------------------------------
    ! Determine pctlnd_o(:)
    ! ----------------------------------------
    allocate(pctlnd_i(ns_i), stat=ier)
    if (ier/=0) call shr_sys_abort('error allocating pctlnd_i')

    call mkpio_get_rawdata(pioid, 'LANDFRAC', mesh_i, pctlnd_i, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call regrid_rawdata(mesh_i, mesh_o, routehandle_nonorm, pctlnd_i, pctlnd_o, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    ! Convert from fraction to percent:
    pctlnd_o(:) = pctlnd_o(:) * 100._r8

    ! ----------------------------------------
    ! Determine pct_nat_pft_o(:,:)
    ! ----------------------------------------
    allocate(pctnatveg_o(ns_o), stat=ier)
    if (ier/=0) call shr_sys_abort('error in allocating pctnatveg_o')

    ! First determine pctnatveg_o(:)
    ! Read in pctnatveg_i
    allocate(pctnatveg_i(ns_i), stat=ier)
    if (ier/=0) call shr_sys_abort('error allocating pctnatveg_i')

    call mkpio_get_rawdata(pioid, 'PCT_NATVEG', mesh_i, pctnatveg_i, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! Regrid to determine pctnatveg_o
    call regrid_rawdata(mesh_i, mesh_o, routehandle_nonorm, pctnatveg_i, pctnatveg_o, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    allocate(pct_nat_pft_i(0:num_natpft,ns_i))
    if (ier/=0) call shr_sys_abort()
    allocate(pct_nat_pft_o(0:num_natpft,ns_o))
    if (ier/=0) call shr_sys_abort()

    ! Read in pct_nat_pft_i
    call mkpio_get_rawdata(pioid, 'PCT_NAT_PFT', mesh_i, pct_nat_pft_i, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    do ni = 1,ns_i
       do m = 0,num_natpft
          pct_nat_pft_i(m,ni) = pct_nat_pft_i(m,ni) * (pctnatveg_i(ni) * 0.01_r8 * mask_i(ni))
       end do
    end do

    ! Readgrid to determine pct_nat_pft_o
    call regrid_rawdata(mesh_i, mesh_o, routehandle_nonorm, pct_nat_pft_i, pct_nat_pft_o, 0, num_natpft, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Rescale pct_nat_pft_o, and set tiny pctnatveg to 0
    do no = 1,ns_o
       if (pctnatveg_o(no) >= 1.0e-6_r8) then
          do m = 0,num_natpft
             pct_nat_pft_o(m,no) = pct_nat_pft_o(m,no) / (pctnatveg_o(no) * 0.01_r8)
          end do
       else
          pctnatveg_o(no) = 0._r8
          pct_nat_pft_o(0,no) = 100._r8
          pct_nat_pft_o(1:num_natpft,no) = 0._r8
       end if

       ! Correct sums so that if they differ slightly from 100, they are
       ! corrected to equal 100 more exactly.
       ! Error check: percents should sum to 100 for land grid cells, within roundoff
       wst_sum = 0.
       do m = 0, num_natpft
          wst_sum = wst_sum + pct_nat_pft_o(m,no)
       enddo
       if (abs(wst_sum - 100._r8) > relerr) then
          write (6,*) subname//'error: nat pft = ', (pct_nat_pft_o(m,no), m = 0, num_natpft), &
               ' do not sum to 100. at no = ',no,' but to ', wst_sum
          call shr_sys_abort()
       end if
       do m = 1, num_natpft
          pct_nat_pft_o(m,no) = pct_nat_pft_o(m,no) * 100._r8 / wst_sum
       end do
    end do

    ! ----------------------------------------
    ! Determine pct_cft_o(:,:)
    ! ----------------------------------------

    ! First Determine pctcrop_o(:)
    allocate(pctcrop_o(ns_o), stat=ier)
    if (ier/=0) call shr_sys_abort('error allocating pctcrop_o_o')
    allocate(pctcrop_i(ns_i), stat=ier)
    if (ier/=0) call shr_sys_abort()
    call mkpio_get_rawdata(pioid, 'PCT_CROP', mesh_i, pctcrop_i, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call regrid_rawdata(mesh_i, mesh_o, routehandle_nonorm, pctcrop_i, pctcrop_o, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    allocate(pct_cft_i(1:num_cft,ns_i), stat=ier)
    if (ier/=0) call shr_sys_abort()
    allocate(pct_cft_o(1:num_cft,ns_o), stat=ier)
    if (ier/=0) call shr_sys_abort()

    ! Get dimensions for PCT_CFT
    allocate(dimlens(3))
    call mkpio_get_dimlengths(pioid, 'PCT_CFT', ndims, dimlens(:))
    if (root_task) then
       do n = 1,ndims
          write(ndiag,'(a,i8,i8)')'   dimid, length= ',n,dimlens(n)
       end do
       write(ndiag,'(a,i8)')'   num_cft = ',num_cft
    end if

    ! Read in pct_cft_i
    if (dimlens(ndims) == num_cft )then
       call mkpio_get_rawdata(pioid, 'PCT_CFT', mesh_i, pct_cft_i, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    else if (dimlens(ndims) > num_cft )then
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
    do ni = 1,ns_i
       do m = 1,num_cft
          pct_cft_i(m,ni) = pct_cft_i(m,ni) * (pctcrop_i(ni) * 0.01_r8 * mask_i(ni))
       end do
    end do

    ! Readgrid pct_cft_i to determine pct_cft_o
    call regrid_rawdata(mesh_i, mesh_o, routehandle_nonorm, pct_cft_i, pct_cft_o, 1, num_cft, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Rescale pct_cft_o, and set tiny pctcrop to 0
    do no = 1,ns_o
       if (pctcrop_o(no) >= 1.0e-6_r8) then
          do m = 1,num_cft
             pct_cft_o(m,no) = pct_cft_o(m,no) / (pctcrop_o(no) * 0.01_r8)
          end do
       else
          pctcrop_o(no) = 0._r8
          pct_cft_o(1,no) = 100._r8
          pct_cft_o(2:num_cft,no) = 0._r8
       end if

       ! Correct sums so that if they differ slightly from 100, they are
       ! corrected to equal 100 more exactly.
       ! Error check: percents should sum to 100 for land grid cells, within roundoff
       wst_sum = 0.
       do m = 1, num_cft
          wst_sum = wst_sum + pct_cft_o(m,no)
       enddo
       if (abs(wst_sum-100._r8) > relerr) then
          write (6,*) subname//'error: crop cft = ',(pct_cft_o(no,m), m = 1, num_cft), &
               ' do not sum to 100. at no = ',no,' but to ', wst_sum
          call shr_sys_abort()
       end if
       do m = 1, num_cft
          pct_cft_o(m,no) = pct_cft_o(m,no) * 100._r8 / wst_sum
       end do
    enddo
    ! ----------------------------------------
    ! Convert % pft as % of grid cell to % pft on the landunit and % of landunit on the grid cell
    ! *** NOTE*** pctnatpft_o and pctcft_o are output arguments
    ! ----------------------------------------
    allocate(output_pct_nat_pft_o(ns_o, 0:num_natpft), stat=ier)
    if (ier/=0) call shr_sys_abort('error in allocating output_pct_nat_pft_o')
    allocate(output_pct_cft_o(ns_o, 1:num_cft), stat=ier)
    if (ier/=0) call shr_sys_abort('error in allocating output_pct_cft_o')

    do no = 1,ns_o
       output_pct_nat_pft_o(no,:) = pct_nat_pft_o(:,no)
       pctnatpft_o(no) = pct_pft_type( output_pct_nat_pft_o(no,:), pctnatveg_o(no), first_pft_index=natpft_lb )

       output_pct_cft_o(no,:) = pct_cft_o(:,no)
       pctcft_o(no) = pct_pft_type( output_pct_cft_o(no,:), pctcrop_o(no), first_pft_index=cft_lb)
    end do

    deallocate(output_pct_nat_pft_o)
    deallocate(output_pct_cft_o)

    ! -----------------------------------------------------------------
    ! Error check
    ! Output global sums on output grid
    ! -----------------------------------------------------------------

    allocate(area_o(ns_o))
    call get_meshareas(mesh_o, area_o, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    allocate(pctpft_o(ns_o,0:(numpft_i-1)), stat=ier)
    if (ier/=0) call shr_sys_abort()
    do no = 1,ns_o
       pctpft_o(no,natpft_lb:natpft_ub) = pctnatpft_o(no)%get_pct_p2g()
       pctpft_o(no,cft_lb:cft_ub)       = pctcft_o(no)%get_pct_p2g()
    end do

    allocate(loc_gpft_o(0:numpft_i-1))
    allocate(glob_gpft_o(0:numpft_i-1))
    loc_gpft_o(:) = 0.
    do no = 1,ns_o
       do m = 0, numpft_i-1
          loc_gpft_o(m) = loc_gpft_o(m) + pctpft_o(no,m) * area_o(no) * frac_o_nonorm(no)
       end do
    end do
    do m = 0,numpft_i-1
       call mpi_reduce(loc_gpft_o(m), glob_gpft_o(m), 1, MPI_REAL8, MPI_SUM, 0, mpicom, ier)
    end do

    if (root_task) then
       write (ndiag,*)
       write (ndiag,'(1x,70a1)') ('.',k=1,70)
       write (ndiag,*) 'PFTs Output'
       write (ndiag,101)
101    format (1x,'plant type     ',20x,' output grid area',/ &
               1x,33x,'      10**6 km**2')
       write (ndiag,'(1x,70a1)') ('.',k=1,70)
       write (ndiag,*)
       do m = 0, numpft_i - 1
          write (ndiag,102) veg(m), glob_gpft_o(m)*1.e-06/100.
       end do
102    format (1x,a35,f17.3)
    end if

    ! Clean up memory
    if (mksrf_fdynuse == ' ') then  ! ...else we will reuse it
       deallocate(frac_o_nonorm)
       call ESMF_RouteHandleDestroy(routehandle_nonorm, nogarbage = .true., rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort()
    end if
    call ESMF_MeshDestroy(mesh_i, nogarbage = .true., rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort()

    if (root_task) then
       write (ndiag,'(a)') 'Successfully made PFTs'
       write (ndiag,*)
    end if

  end subroutine mkpft

end module mkpftMod
