module mkpftMod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: mkpft
!
! !DESCRIPTION:
! Make PFT data
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!-----------------------------------------------------------------------
!!USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use shr_sys_mod , only : shr_sys_flush
  use mkvarctl    , only : numpft
  use mkdomainMod , only : domain_checksame
  use mkpftConstantsMod

  implicit none

  private           ! By default make data private
!
! !PUBLIC MEMBER FUNCTIONS:
!
  public mkpftInit          ! Initialization
  public mkpft              ! Set PFT
  public mkpft_parse_oride  ! Parse the string with PFT fraction/index info to override
  public mkpftAtt           ! Write out attributes to output file on pft
!
! !PUBLIC DATA MEMBERS: 
!

  !
  ! When pft_idx and pft_frc are set, they must be set together, and they will cause the
  ! entire area to be covered with vegetation and zero out other landunits.
  ! The sum of pft_frc must = 100%, and each pft_idx point in the array corresponds to
  ! the fraction in pft_frc. Only the first few points are used until pft_frc = 0.0.
  !
  integer            :: m                     ! index
  integer, public    :: pft_idx(0:maxpft) = & ! PFT vegetation index to override with
                             (/ ( -1,  m = 0, maxpft )   /)
  real(r8), public   :: pft_frc(0:maxpft) = & ! PFT vegetation fraction to override with
                             (/ ( 0.0_r8, m = 0, maxpft ) /)
!
! !PRIVATE DATA MEMBERS:
!
  logical, private :: zero_out      = .false. ! Flag to zero out PFT
  logical, private :: use_input_pft = .false. ! Flag to override PFT with input values
  integer, private :: nzero                   ! index of first zero fraction
!
! !PRIVATE MEMBER FUNCTIONS:
!
  private :: mkpft_check_oride  ! Check the pft_frc and pft_idx values for correctness
!EOP
!===============================================================
contains
!===============================================================

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: mkpftInit
!
! !INTERFACE:
subroutine mkpftInit( zero_out_l, all_veg )
!
! !DESCRIPTION:
! Initialize of Make PFT data
! !USES:
  use mkvarpar, only : numstdpft, numstdcft
!
! !ARGUMENTS:
  implicit none
  logical, intent(IN)  :: zero_out_l ! If veg should be zero'ed out
  logical, intent(OUT) :: all_veg    ! If should zero out other fractions so that
                                     ! all land-cover is vegetation
!
! !CALLED FROM:
! subroutine mksrfdat in module mksrfdatMod
!
! !REVISION HISTORY:
! Author: Erik Kluzek
!
!
! !LOCAL VARIABLES:
!EOP
  real(r8), parameter :: hndrd = 100.0_r8  ! A hundred percent
  character(len=32) :: subname = 'mkpftMod::mkpftInit() '
!-----------------------------------------------------------------------
  write (6, '(a, a, a)') "In ", trim(subname), "..."
  call mkpft_check_oride( )
  if ( use_input_pft ) then
     if ( maxpft < numpft ) then
        write(6,*) subname//'number PFT is > max allowed!'
        call abort()
     end if
     write(6,*) 'Set PFT fraction to : ', pft_frc(0:nzero-1)
     write(6,*) 'With PFT index      : ', pft_idx(0:nzero-1)
  end if

  all_veg = use_input_pft

  if ( zero_out_l .and. all_veg )then
     write(6,*) subname//'zeroing out vegetation and setting vegetation to 100% is a contradiction!'
     call abort()
  end if

  ! Copy local zero out to module data version
  zero_out = zero_out_l

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
  natpft_lb = 0
  natpft_ub = num_natpft
  cft_lb    = num_natpft+1
  cft_ub    = cft_lb + num_cft - 1

  ! Make sure the array indices have been set up properly, to ensure the 1:1
  ! correspondence mentioned above
  if (cft_ub /= numpft) then
     write(6,*) 'CFT_UB set up incorrectly: cft_ub, numpft = ', cft_ub, numpft
     call abort()
  end if

end subroutine mkpftInit

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: mkpft
!
! !INTERFACE:
subroutine mkpft(ldomain, mapfname, fpft, ndiag, allow_no_crops, &
     pctlnd_o, pctnatpft_o, pctcft_o, &
     pctcft_o_saved)
!
! !DESCRIPTION:
! Make PFT data
!
! This dataset consists of the %cover of the [numpft]+1 PFTs used by
! the model. The input %cover pertains to the "vegetated" portion of the
! grid cell and sums to 100. The real portion of each grid cell
! covered by each PFT is the PFT cover times the fraction of the
! grid cell that is land. This is the quantity preserved when
! area-averaging from the input (1/2 degree) grid to the models grid.
!
! Upon return from this routine, the % cover of the natural veg + crop landunits is
! generally 100% everywhere; this will be normalized later to account for special landunits.
!
! If allow_no_crops is true, then we allow the input dataset to have no prognostic crop
! information (i.e., only contain information about the "standard" PFTs). In this case,
! pctcft_o_saved MUST be given. If the input dataset is found to not have information
! about the prognostic crops, then we take the generic c3 crop cover from the input
! dataset to specify the crop landunit area, and we take the individual crop breakdown
! from pctcft_o_saved (which will generally have come from some other input dataset that
! DID contain prognostic crop information).
!
! !USES:
  use mkdomainMod, only : domain_type, domain_clean, domain_read
  use mkgridmapMod
  use mkvarpar	
  use mkvarctl    
  use mkncdio
  use mkpctPftTypeMod,   only : pct_pft_type
  use mkpftUtilsMod,     only : convert_from_p2g
  use mkpftConstantsMod, only : natpft_lb, natpft_ub, num_cft, cft_lb, cft_ub
!
! !ARGUMENTS:
  implicit none
  type(domain_type), intent(inout) :: ldomain
  character(len=*)  , intent(in) :: mapfname              ! input mapping file name
  character(len=*)  , intent(in) :: fpft                  ! input pft dataset file name
  integer           , intent(in) :: ndiag                 ! unit number for diag out
  logical           , intent(in) :: allow_no_crops        ! if it's okay to not have prognostic crops in the input file
  real(r8)          , intent(out):: pctlnd_o(:)           ! output grid:%land/gridcell
  type(pct_pft_type), intent(out):: pctnatpft_o(:)        ! natural PFT cover
  type(pct_pft_type), intent(out):: pctcft_o(:)           ! crop (CFT) cover

! saved crop cover information, in case the input dataset does not contain information about prognostic crops
  type(pct_pft_type), intent(in), optional :: pctcft_o_saved(:)
!
! !CALLED FROM:
! subroutine mksrfdat in module mksrfdatMod
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!
! !LOCAL VARIABLES:
!EOP
  type(pct_pft_type), allocatable:: pctnatpft_i(:)         ! input grid: natural PFT cover
  type(pct_pft_type), allocatable:: pctcft_i(:)            ! input grid: crop (CFT) cover
  type(domain_type)    :: tdomain             ! local domain
  type(gridmap_type)    :: tgridmap           ! local gridmap
  real(r8), allocatable :: pctpft_i(:,:)      ! input grid: PFT percent
  real(r8), allocatable :: pctpft_o(:,:)      ! output grid: PFT percent (% of grid cell)
  real(r8), allocatable :: pctnatveg_i(:)     ! input grid: natural veg percent (% of grid cell)
  real(r8), allocatable :: pctnatveg_o(:)     ! output grid: natural veg percent (% of grid cell)
  real(r8), allocatable :: pctcrop_i(:)       ! input grid: all crop percent (% of grid cell)
  real(r8), allocatable :: pctcrop_o(:)       ! output grid: all crop percent (% of grid cell)
  real(r8), allocatable :: pct_cft_i(:,:)     ! input grid: CFT (Crop Functional Type) percent (% of landunit cell)
  real(r8), allocatable :: temp_i(:,:)        ! input grid: temporary 2D variable to read in
  real(r8), allocatable :: pct_cft_o(:,:)     ! output grid: CFT (Crop Functional Type) percent (% of landunit cell)
  real(r8), allocatable :: pct_nat_pft_i(:,:) ! input grid: natural PFT percent (% of landunit cell)
  real(r8), allocatable :: pct_nat_pft_o(:,:) ! output grid: natural PFT percent (% of landunit cell)
  integer  :: numpft_i                        ! num of plant types input data
  integer  :: natpft_i                        ! num of natural plant types input data
  integer  :: ncft_i                          ! num of crop types input data
  real(r8) :: sum_fldo                        ! global sum of dummy output fld
  real(r8) :: sum_fldi                        ! global sum of dummy input fld
  real(r8) :: wst_sum                         ! sum of %pft
  real(r8), allocatable :: gpft_o(:)          ! output grid: global area pfts
  real(r8) :: garea_o                         ! output grid: global area
  real(r8), allocatable :: gpft_i(:)          ! input grid: global area pfts
  real(r8) :: garea_i                         ! input grid: global area
  integer  :: k,n,m,ni,no,ns_i,ns_o           ! indices
  integer  :: ncid,dimid,varid                ! input netCDF id's
  integer  :: ndims                           ! number of dimensions for a variable on the file
  integer  :: dimlens(3)                      ! dimension lengths for a variable on the file
  integer  :: ier                             ! error status
  logical  :: missing_crops                   ! if we need prognostic crop info, but the input dataset is missing this crop info
  real(r8) :: relerr = 0.0001_r8              ! max error: sum overlap wts ne 1
  logical  :: oldformat                       ! if input file is in the old format or not (based on what variables exist)

  character(len=35)  veg(0:maxpft)            ! vegetation types
  character(len=32) :: subname = 'mkpftMod::mkpft()'
!-----------------------------------------------------------------------

  write (6,*)
  write (6, '(a, a, a)') "In ", trim(subname), "..."
  write (6,*) 'Attempting to make PFTs .....'
  call shr_sys_flush(6)

  if (allow_no_crops) then
     if (.not. present(pctcft_o_saved)) then
        write(6,*) subname, ' ERROR: when allow_no_crops is true, pctcft_o_saved must be given'
        call abort()
     end if
  end if

  ! Start by assuming the input dataset is NOT missing crop info
  missing_crops = .false.

  ! -----------------------------------------------------------------
  ! Set the vegetation types
  ! -----------------------------------------------------------------
  if ( numpft >= numstdpft )then
     veg(0:maxpft) = (/                                   &
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
  if (      numpft == numstdpft )then
     write(6,*)'Creating surface datasets with the standard # of PFTs =', numpft
  else if ( numpft > numstdpft )then
     write(6,*)'Creating surface datasets with extra types for crops; total pfts =', numpft
  else
     write(6,*) subname//': parameter numpft is NOT set to a known value (should be 16 or more) =',numpft
     call abort()
  end if

  ! -----------------------------------------------------------------
  ! Read input PFT file
  ! -----------------------------------------------------------------

  ns_o = ldomain%ns

  if ( .not. use_input_pft ) then
     ! Obtain input grid info, read PCT_PFT

     call domain_read(tdomain,fpft)
     ns_i = tdomain%ns

     write (6,*) 'Open PFT file: ', trim(fpft)
     call check_ret(nf_open(fpft, 0, ncid), subname)

     ! Check what variables exist to determine what format the file is in
     call check_ret(nf_inq_varid (ncid, 'PCT_PFT', varid), subname, varexists=oldformat)

     if ( .not. oldformat ) then
        call check_ret(nf_inq_dimid  (ncid, 'natpft', dimid), subname)
        call check_ret(nf_inq_dimlen (ncid, dimid, natpft_i), subname)
        call check_ret(nf_inq_dimid  (ncid, 'cft', dimid), subname)
        call check_ret(nf_inq_dimlen (ncid, dimid, ncft_i), subname)
        numpft_i = natpft_i + ncft_i
     else
        call check_ret(nf_inq_dimid  (ncid, 'pft', dimid), subname)
        call check_ret(nf_inq_dimlen (ncid, dimid, numpft_i), subname)
     end if

     ! Check if the number of pfts on the input matches the expected number. A mismatch
     ! is okay in the case that the input has the standard number of pfts (i.e., no
     ! prognostic crop info), if allow_no_crops is true. Otherwise, a mismatch is an error.
     if (numpft_i .ne. numpft+1) then
        if (numpft_i .eq. numstdpft+1) then
           if (allow_no_crops) then
              write(6,*) subname//': using non-crop input file for a surface dataset with crops'
              write(6,*) "(this is okay: we'll use the saved crop breakdown from the non-transient input file)"
              missing_crops = .true.
           else
              write(6,*) subname//' ERROR: trying to use non-crop input file'
              write(6,*) 'for a surface dataset with crops, but allow_no_crops is false'
              write(6,*) "(This can happen if you're trying to use a non-crop input file"
              write(6,*) "for the surface dataset itself: a non-crop input file is only"
              write(6,*) "allowed for the transient PFT information.)"
              call abort()
           end if
        else if (numpft_i > numstdpft+1 .and. numpft_i == maxpft+1 .and. .not. oldformat) then
           write(6,*) subname//' WARNING: using a crop input raw dataset for a non-crop output surface dataset'
        else
           write(6,*) subname//': parameter numpft+1= ',numpft+1, &
                'does not equal input dataset numpft= ',numpft_i
           call abort()
        end if
     endif


     ! If file is in the new format, expect the following variables: 
     !      PCT_NATVEG, PCT_CROP, PCT_NAT_PFT, PCT_CFT
     if ( .not. oldformat )then
        allocate(pctnatveg_i(ns_i), &
                 pctnatveg_o(ns_o), &
                 pctcrop_i(ns_i),   &
                 pctcrop_o(ns_o),   &
                 pct_cft_i(ns_i,1:num_cft), &
                 pct_cft_o(ns_o,1:num_cft), &
                 pct_nat_pft_i(ns_i,0:num_natpft), &
                 pct_nat_pft_o(ns_o,0:num_natpft), &
                 stat=ier)
        if (ier/=0) call abort()

        call check_ret(nf_inq_varid (ncid, 'PCT_NATVEG', varid), subname)
        call check_ret(nf_get_var_double (ncid, varid, pctnatveg_i), subname)
        call check_ret(nf_inq_varid (ncid, 'PCT_CROP', varid), subname)
        call check_ret(nf_get_var_double (ncid, varid, pctcrop_i), subname)
        call check_ret(nf_inq_varid (ncid, 'PCT_CFT', varid), subname)
        call get_dim_lengths(ncid, 'PCT_CFT', ndims, dimlens(:) )
        if (      ndims == 3 .and. dimlens(1)*dimlens(2) == ns_i .and. dimlens(3) == num_cft )then
           call check_ret(nf_get_var_double (ncid, varid, pct_cft_i), subname)
        else if ( ndims == 3 .and. dimlens(1)*dimlens(2) == ns_i .and. dimlens(3) > num_cft )then
           ! Read in the whole array: then sum the rainfed and irrigated
           ! seperately
           allocate( temp_i(ns_i,dimlens(3)) )
           call check_ret(nf_get_var_double (ncid, varid, temp_i), subname)
           do n = 1, num_cft
              pct_cft_i(:,n) = 0.0_r8
              do m = n, dimlens(3), 2
                 pct_cft_i(:,n) = pct_cft_i(:,n) + temp_i(:,m)
              end do
           end do
           deallocate( temp_i )
        else
           write(6,*) subname//': ERROR: dimensions for PCT_CROP are NOT what is expected'
           call abort()
        end if
        call check_ret(nf_inq_varid (ncid, 'PCT_NAT_PFT', varid), subname)
        call check_ret(nf_get_var_double (ncid, varid, pct_nat_pft_i), subname)

     ! Read in from the old format with PCT_PFT alone
     else
        allocate(pctpft_i(ns_i,0:(numpft_i-1)), &
                 pctpft_o(ns_o,0:(numpft_i-1)), &
                 stat=ier)
        if (ier/=0) call abort()

        call check_ret(nf_inq_varid (ncid, 'PCT_PFT', varid), subname)
        call check_ret(nf_get_var_double (ncid, varid, pctpft_i), subname)
     end if

     call check_ret(nf_close(ncid), subname)

  else
     oldformat = .true.
     ns_i = 1
     numpft_i = numpft+1
     allocate(pctpft_o(ns_o,0:numpft), stat=ier)
     if (ier/=0) call abort()
  end if

  ! Determine pctpft_o on output grid

  if ( zero_out ) then

     pctpft_o(:,:)  = 0._r8
     pctlnd_o(:)    = 100._r8
     pctnatveg_o(:) = 0._r8
     pctcrop_o(:)   = 0._r8
     pct_nat_pft_o(:,:) =   0._r8
     pct_nat_pft_o(:,0) = 100._r8
     pct_cft_o(:,:)     =   0._r8
     pct_cft_o(:,1)     = 100._r8

  else if ( use_input_pft ) then

     call mkpft_check_oride( )

     ! set PFT based on input pft_frc and pft_idx
     pctpft_o(:,:) = 0._r8
     pctlnd_o(:)   = 100._r8
     do m = 0, numpft
        ! Once reach a PFT where fraction goes to zero -- exit
        if ( pft_frc(m) .eq. 0.0_r8 ) exit
        do no = 1,ns_o
           pctpft_o(no,pft_idx(m)) = pft_frc(m)
        end do
     end do

  else

     ! Compute pctlnd_o, pctpft_o

     call gridmap_mapread(tgridmap, mapfname)

     ! Error checks for domain and map consistencies

     call domain_checksame( tdomain, ldomain, tgridmap )

     ! Area-average percent cover on input grid [pctpft_i] to output grid 
     ! [pctpft_o] and correct [pctpft_o] according to land landmask
     ! Note that percent cover is in terms of total grid area.
  
     do no = 1,ns_o
        pctlnd_o(no)     = tgridmap%frac_dst(no) * 100._r8
        ldomain%frac(no) = tgridmap%frac_dst(no) 
     end do

     ! New format with extra variables on input
     if ( .not. oldformat ) then
        call gridmap_areaave(tgridmap, pctnatveg_i, pctnatveg_o, nodata=0._r8)
        call gridmap_areaave(tgridmap, pctcrop_i,   pctcrop_o,   nodata=0._r8)

        do m = 0, num_natpft
           call gridmap_areaave_scs(tgridmap, pct_nat_pft_i(:,m), pct_nat_pft_o(:,m), nodata=0._r8,src_wt=pctnatveg_i*0.01_r8,dst_wt=pctnatveg_o*0.01_r8)
           do no = 1,ns_o
              if (pctlnd_o(no) < 1.0e-6 .or. pctnatveg_o(no) < 1.0e-6) then
                 if (m == 0) then
                    pct_nat_pft_o(no,m) = 100._r8
                 else
                    pct_nat_pft_o(no,m) = 0._r8
                 endif
              end if
           enddo
        end do
        do m = 1, num_cft
           call gridmap_areaave_scs(tgridmap, pct_cft_i(:,m), pct_cft_o(:,m), nodata=0._r8,src_wt=pctcrop_i*0.01_r8,dst_wt=pctcrop_o*0.01_r8)
           do no = 1,ns_o
              if (pctlnd_o(no) < 1.0e-6 .or. pctcrop_o(no) < 1.0e-6) then
                 if (m == 1) then
                    pct_cft_o(no,m) = 100._r8
                 else
                    pct_cft_o(no,m) = 0._r8
                 endif
              end if
           enddo
        end do
     ! Old format with just PCTPFT
     else
        do m = 0, numpft_i - 1
           call gridmap_areaave(tgridmap, pctpft_i(:,m), pctpft_o(:,m), nodata=0._r8)
           do no = 1,ns_o
              if (pctlnd_o(no) < 1.0e-6) then
                 if (m == 0) then
                    pctpft_o(no,m) = 100._r8
                 else
                    pctpft_o(no,m) = 0._r8
                 endif
              end if
           enddo
        enddo
     end if

  end if

  ! Error check: percents should sum to 100 for land grid cells, within roundoff
  ! Also correct sums so that if they differ slightly from 100, they are corrected to
  ! equal 100 more exactly.

  if ( (.not. zero_out) .and. oldformat) then
     do no = 1,ns_o
        wst_sum = 0.
        do m = 0, numpft_i - 1
           wst_sum = wst_sum + pctpft_o(no,m)
        enddo
        if (abs(wst_sum-100._r8) > relerr) then
           write (6,*) subname//'error: pft = ', &
                (pctpft_o(no,m), m = 0, numpft_i-1), &
                ' do not sum to 100. at no = ',no,' but to ', wst_sum
           stop
        end if

        ! Correct sum so that if it differs slightly from 100, it is corrected to equal
        ! 100 more exactly
        do m = 0, numpft_i - 1
           pctpft_o(no,m) = pctpft_o(no,m) * 100._r8 / wst_sum
        end do

     end do
  else if ( (.not. zero_out) .and. (.not. oldformat) ) then
     do no = 1,ns_o
        wst_sum = 0.
        do m = 0, num_natpft
           wst_sum = wst_sum + pct_nat_pft_o(no,m)
        enddo
        if (abs(wst_sum-100._r8) > relerr) then
           write (6,*) subname//'error: nat pft = ', &
                (pct_nat_pft_o(no,m), m = 0, num_natpft), &
                ' do not sum to 100. at no = ',no,' but to ', wst_sum
           stop
        end if

        ! Correct sum so that if it differs slightly from 100, it is corrected to equal
        ! 100 more exactly
        do m = 1, num_natpft
           pct_nat_pft_o(no,m) = pct_nat_pft_o(no,m) * 100._r8 / wst_sum
        end do

        wst_sum = 0.
        do m = 1, num_cft
           wst_sum = wst_sum + pct_cft_o(no,m)
        enddo
        if (abs(wst_sum-100._r8) > relerr) then
           write (6,*) subname//'error: crop cft = ', &
                (pct_cft_o(no,m), m = 1, num_cft), &
                ' do not sum to 100. at no = ',no,' but to ', wst_sum
           stop
        end if

        ! Correct sum so that if it differs slightly from 100, it is corrected to equal
        ! 100 more exactly
        do m = 1, num_cft
           pct_cft_o(no,m) = pct_cft_o(no,m) * 100._r8 / wst_sum
        end do

     end do
  end if

  ! Convert % pft as % of grid cell to % pft on the landunit and % of landunit on the
  ! grid cell
  if (missing_crops) then
     do no = 1,ns_o
        call convert_from_p2g(pct_p2g=pctpft_o(no,:), pctcft_saved=pctcft_o_saved(no), &
             pctnatpft=pctnatpft_o(no), pctcft=pctcft_o(no))
     end do
  else if ( .not. oldformat ) then
     do no = 1,ns_o
        pctnatpft_o(no) = pct_pft_type( pct_nat_pft_o(no,:), pctnatveg_o(no), first_pft_index=natpft_lb )
        pctcft_o(no)    = pct_pft_type( pct_cft_o(no,:),     pctcrop_o(no),   first_pft_index=cft_lb    )
     end do
  else
     do no = 1,ns_o
        call convert_from_p2g(pct_p2g=pctpft_o(no,:), &
             pctnatpft=pctnatpft_o(no), pctcft=pctcft_o(no))
     end do
  end if

  ! -----------------------------------------------------------------
  ! Error check
  ! Compare global areas on input and output grids
  ! -----------------------------------------------------------------

  if ( .not. (zero_out .or. use_input_pft) ) then

     ! Convert to pctpft over grid if using new format
     if ( .not. oldformat ) then
        allocate(pctpft_i(ns_i,0:(numpft_i-1)), &
                 pctpft_o(ns_o,0:(numpft_i-1)), &
                 pctnatpft_i(ns_i),             &
                 pctcft_i(ns_i),                &
                 stat=ier)
        if (ier/=0) call abort()
        do ni = 1, ns_i
           pctnatpft_i(ni) = pct_pft_type( pct_nat_pft_i(ni,:), pctnatveg_i(ni), first_pft_index=natpft_lb )
           pctcft_i(ni)    = pct_pft_type( pct_cft_i(ni,:),     pctcrop_i(ni),   first_pft_index=cft_lb    )
        end do

        do no = 1,ns_o
           pctpft_o(no,natpft_lb:natpft_ub) = pctnatpft_o(no)%get_pct_p2g()
           pctpft_o(no,cft_lb:cft_ub)       = pctcft_o(no)%get_pct_p2g()
        end do
     end if
     allocate(gpft_i(0:numpft_i-1))
     allocate(gpft_o(0:numpft_i-1))

     ! input grid

     gpft_i(:) = 0.
     garea_i   = 0.
     do ni = 1,ns_i
        garea_i = garea_i + tgridmap%area_src(ni)*re**2
        do m = 0, numpft_i - 1
           gpft_i(m) = gpft_i(m) + pctpft_i(ni,m)*tgridmap%area_src(ni)*&
                                                  tgridmap%frac_src(ni)*re**2
        end do
     end do
     if ( allocated(pctpft_i) ) deallocate (pctpft_i)

     ! output grid

     gpft_o(:) = 0.
     garea_o   = 0.
     do no = 1,ns_o
        garea_o = garea_o + tgridmap%area_dst(no)*re**2
        do m = 0, numpft_i - 1
           gpft_o(m) = gpft_o(m) + pctpft_o(no,m)*tgridmap%area_dst(no)*&
                                                  tgridmap%frac_dst(no)*re**2
        end do
     end do

     ! comparison

     write (ndiag,*)
     write (ndiag,'(1x,70a1)') ('=',k=1,70)
     write (ndiag,*) 'PFTs Output'
     write (ndiag,'(1x,70a1)') ('=',k=1,70)

     write (ndiag,*)
     write (ndiag,'(1x,70a1)') ('.',k=1,70)
     write (ndiag,1001)
1001 format (1x,'plant type     ',20x,' input grid area',' output grid area',/ &
             1x,33x,'     10**6 km**2','      10**6 km**2')
     write (ndiag,'(1x,70a1)') ('.',k=1,70)
     write (ndiag,*)
     do m = 0, numpft_i - 1
        write (ndiag,1002) veg(m), gpft_i(m)*1.e-06/100.,gpft_o(m)*1.e-06/100.
     end do
1002 format (1x,a35,f16.3,f17.3)
     call shr_sys_flush(ndiag)

     deallocate(gpft_i, gpft_o)
     if ( .not. oldformat ) then
        deallocate(pctpft_o)
     end if

  end if


  ! Deallocate dynamic memory

  if ( .not. oldformat ) then
     deallocate(pctnatveg_i)
     deallocate(pctnatveg_o)
     deallocate(pctcrop_i)
     deallocate(pctcrop_o)
     deallocate(pct_cft_i)
     deallocate(pct_cft_o)
     deallocate(pct_nat_pft_i)
     deallocate(pct_nat_pft_o)
  else
     deallocate(pctpft_o)
  end if
  call domain_clean(tdomain) 
  if ( .not. zero_out .and. .not. use_input_pft ) then
     call gridmap_clean(tgridmap)
  end if

  write (6,*) 'Successfully made PFTs'
  write (6,*)


end subroutine mkpft

!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: mkpft_parse_oride
!
! !INTERFACE:
subroutine mkpft_parse_oride( string )
!
! !DESCRIPTION:
! Parse the string with pft fraction and index information on it, to override
! the file with this information rather than reading from a file.
!
! !USES:
   use shr_string_mod, only: shr_string_betweenTags, shr_string_countChar
! !ARGUMENTS:
   character(len=256), intent(IN) :: string  ! String to parse with PFT fraction 
                                             ! and index data
!
! !CALLED FROM:
! subroutine mksrfdat in module mksrfdatMod
!
! !REVISION HISTORY:
! Author: Erik Kluzek
!
!
! !LOCAL VARIABLES:
!EOP
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
     call abort()
  end if
  num_elms = shr_string_countChar( substring, ",", rc )
  read(substring,*) pft_frc(0:num_elms)
  call shr_string_betweenTags( string, idx_start, idx_end, substring, rc )
  if ( rc /= 0 )then
     write(6,*) subname//'Trouble finding pft_index start end tags'
     call abort()
  end if
  if ( num_elms /= shr_string_countChar( substring, ",", rc ) )then
     write(6,*) subname//'number of elements different between frc and idx fields'
     call abort()
  end if
  read(substring,*) pft_idx(0:num_elms)
!-----------------------------------------------------------------------

end subroutine mkpft_parse_oride

!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: mkpft_check_oride
!
! !INTERFACE:
subroutine mkpft_check_oride( )
!
! !DESCRIPTION:
! Check that the pft override values are valid
! !USES:
!
! !ARGUMENTS:
  implicit none
!
! !REVISION HISTORY:
! Author: Erik Kluzek
!
!
! !LOCAL VARIABLES:
!EOP
  integer  :: i, j                         ! indices
  real(r8) :: sumpft                       ! Sum of pft_frc
  real(r8), parameter :: hndrd = 100.0_r8  ! A hundred percent
  character(len=32) :: subname = 'mkpftMod::mkpft_check_oride() '
!-----------------------------------------------------------------------

  sumpft = sum(pft_frc)
  if (          sumpft == 0.0 )then
    ! PFT fraction is NOT used
    use_input_pft = .false.
  else if ( abs(sumpft - hndrd) > 1.e-6 )then
    write(6, '(a, a, f15.12)') trim(subname), 'Sum of PFT fraction is NOT equal to 100% =', sumpft
    write(6,*) 'Set PFT fraction to : ', pft_frc(0:nzero-1)
    write(6,*) 'With PFT index      : ', pft_idx(0:nzero-1)
    call abort()
  else
    use_input_pft = .true.
    nzero = 0
    do i = 0, numpft
       if ( pft_frc(i) == 0.0_r8 )then
          nzero = i
          exit
       end if
    end do
    ! PFT fraction IS used, and sum is OK, now check details
    do i = 0, nzero -1
      if ( pft_frc(i) < 0.0_r8 .or. pft_frc(i) > hndrd )then
         write(6,*) subname//'PFT fraction is out of range: pft_frc=', pft_frc(i)
         call abort()
      else if ( pft_frc(i) > 0.0_r8 .and. pft_idx(i) == -1 )then
         write(6,*) subname//'PFT fraction > zero, but index NOT set: pft_idx=', pft_idx(i)
         call abort()
      end if
      ! PFT index out of range
      if ( pft_idx(i) < 0 .or. pft_idx(i) > numpft )then
         write(6,*) subname//'PFT index is out of range: ', pft_idx(i)
         call abort()
      end if
      ! Make sure index values NOT used twice
      do j = 0, i-1
         if ( pft_idx(i) == pft_idx(j) )then
            write(6,*) subname//'Same PFT index is used twice: ', pft_idx(i)
            call abort()
         end if
      end do
    end do
    ! Make sure the rest of the fraction is zero and index are not set as well
    do i = nzero, numpft
      if ( pft_frc(i) /= 0.0_r8 .or. pft_idx(i) /= -1 )then
         write(6,*) subname//'After PFT fraction is zeroed out, fraction is non zero, or index set'
         call abort()
      end if
    end do
  end if

end subroutine mkpft_check_oride

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: mkpftAtt
!
! !INTERFACE:
subroutine mkpftAtt( ncid, dynlanduse, xtype )
!
! !DESCRIPTION:
! make PFT attributes on the output file
!
  use mkncdio    , only : check_ret, ncd_defvar, ncd_def_spatial_var
  use fileutils  , only : get_filename
  use mkvarctl   , only : mksrf_fvegtyp, mksrf_flai
  use mkvarpar   

! !ARGUMENTS:
  implicit none
  include 'netcdf.inc'
  integer, intent(in) :: ncid         ! NetCDF file ID to write out to
  logical, intent(in) :: dynlanduse   ! if dynamic land-use file
  integer, intent(in) :: xtype        ! external type to output real data as
!
! !CALLED FROM:
! subroutine mkfile in module mkfileMod
!
! !REVISION HISTORY:
! Author: Erik Kluzek
!
!
! !LOCAL VARIABLES:
!EOP
  integer :: pftsize              ! size of lsmpft dimension
  integer :: natpftsize           ! size of natpft dimension
  integer :: dimid                ! input netCDF id's
  character(len=256) :: str       ! global attribute string
  character(len=32) :: subname = 'mkpftAtt'

  ! Define dimensions
  call check_ret(nf_def_dim (ncid, 'time'   , nf_unlimited, dimid), subname)

  if (.not. dynlanduse) then
     pftsize = numpft + 1
     call check_ret(nf_def_dim (ncid, 'lsmpft' , pftsize     , dimid), subname)
  end if

  natpftsize = num_natpft + 1
  call check_ret(nf_def_dim (ncid, 'natpft' , natpftsize  , dimid), subname)

  ! zero-size dimensions can cause problems, so we only include the cft dimension if num_cft > 0
  ! Note that this implies that we can only include PCT_CFT on the dataset if num_cft > 0
  if (num_cft > 0) then
     call check_ret(nf_def_dim (ncid, 'cft'    , num_cft     , dimid), subname)
  end if

  ! Add global attributes

  if (.not. dynlanduse) then
     str = get_filename(mksrf_flai)
     call check_ret(nf_put_att_text(ncid, NF_GLOBAL, &
            'Lai_raw_data_file_name', len_trim(str), trim(str)), subname)
  end if

  if ( use_input_pft ) then
     str = 'TRUE'
     call check_ret(nf_put_att_text (ncid, NF_GLOBAL, &
          'pft_override', len_trim(str), trim(str)), subname)
  else if ( zero_out )then
     str = 'TRUE'
     call check_ret(nf_put_att_text (ncid, NF_GLOBAL, &
          'zero_out_pft_override', len_trim(str), trim(str)), subname)
  else
     str = get_filename(mksrf_fvegtyp)
     call check_ret(nf_put_att_text(ncid, NF_GLOBAL, &
       'Vegetation_type_raw_data_filename', len_trim(str), trim(str)), subname)
  end if

  ! Define variables

  ! Coordinate variable for indices of natural PFTs
  call ncd_defvar(ncid=ncid, varname='natpft', xtype=nf_int, &
       dim1name='natpft', long_name='indices of natural PFTs', units='index')

  ! Coordinate variable for indices of CFTs
  if (num_cft > 0) then
     call ncd_defvar(ncid=ncid, varname='cft', xtype=nf_int, &
          dim1name='cft', long_name='indices of CFTs', units='index')
  end if

  call ncd_def_spatial_var(ncid=ncid, varname='LANDFRAC_PFT', xtype=nf_double, &
       long_name='land fraction from pft dataset', units='unitless')
  
  call ncd_def_spatial_var(ncid=ncid, varname='PFTDATA_MASK', xtype=nf_int, &
       long_name='land mask from pft dataset, indicative of real/fake points', units='unitless')
  
  if (.not. dynlanduse) then
     call ncd_def_spatial_var(ncid=ncid, varname='PCT_NATVEG', xtype=xtype, &
          long_name='total percent natural vegetation landunit', units='unitless')
  end  if

  ! PCT_CROP
  if (.not. dynlanduse) then
     call ncd_def_spatial_var(ncid=ncid, varname='PCT_CROP', xtype=xtype, &
          long_name='total percent crop landunit', units='unitless')
  else
     call ncd_def_spatial_var(ncid=ncid, varname='PCT_CROP', xtype=xtype, &
          lev1name='time', &
          long_name='total percent crop landunit', units='unitless')
     call ncd_def_spatial_var(ncid=ncid, varname='PCT_CROP_MAX', xtype=xtype, &
          long_name='maximum total percent crop landunit during time period', units='unitless')
  end if

  ! PCT_NAT_PFT
  if (.not. dynlanduse) then
     call ncd_def_spatial_var(ncid=ncid, varname='PCT_NAT_PFT', xtype=xtype, &
          lev1name='natpft', &
          long_name='percent plant functional type on the natural veg landunit (% of landunit)', units='unitless')
  else
     call ncd_def_spatial_var(ncid=ncid, varname='PCT_NAT_PFT', xtype=xtype, &
          lev1name='natpft', lev2name='time', &
          long_name='percent plant functional type on the natural veg landunit (% of landunit)', units='unitless')
     call ncd_def_spatial_var(ncid=ncid, varname='PCT_NAT_PFT_MAX', xtype=xtype, &
          lev1name='natpft', &
          long_name='maximum percent plant functional type during time period (% of landunit)', units='unitless')
  end if

  ! PCT_CFT
  if (num_cft > 0) then
     if (.not. dynlanduse) then
        call ncd_def_spatial_var(ncid=ncid, varname='PCT_CFT', xtype=xtype, &
             lev1name='cft', &
             long_name='percent crop functional type on the crop landunit (% of landunit)', units='unitless')
     else
        call ncd_def_spatial_var(ncid=ncid, varname='PCT_CFT', xtype=xtype, &
             lev1name='cft', lev2name='time', &
             long_name='percent crop functional type on the crop landunit (% of landunit)', units='unitless')
        call ncd_def_spatial_var(ncid=ncid, varname='PCT_CFT_MAX', xtype=xtype, &
             lev1name='cft', &
             long_name='maximum percent crop functional type during time period (% of landunit)', units='unitless')
     end if
  end if

  ! LAI,SAI,HTOP,HBOT
  if (.not. dynlanduse) then
     call ncd_def_spatial_var(ncid=ncid, varname='MONTHLY_LAI', xtype=xtype,  &
          lev1name='lsmpft', lev2name='time', &
          long_name='monthly leaf area index', units='unitless')

     call ncd_def_spatial_var(ncid=ncid, varname='MONTHLY_SAI', xtype=xtype,  &
          lev1name='lsmpft', lev2name='time', &
          long_name='monthly stem area index', units='unitless')

     call ncd_def_spatial_var(ncid=ncid, varname='MONTHLY_HEIGHT_TOP', xtype=xtype,  &
          lev1name='lsmpft', lev2name='time', &
          long_name='monthly height top', units='meters')

     call ncd_def_spatial_var(ncid=ncid, varname='MONTHLY_HEIGHT_BOT', xtype=xtype,  &
          lev1name='lsmpft', lev2name='time', &
          long_name='monthly height bottom', units='meters')
  end if

  ! OTHER
  if (dynlanduse) then
     call ncd_defvar(ncid=ncid, varname='YEAR', xtype=nf_int,  &
            dim1name='time', &
            long_name='Year of PFT data', units='unitless')
     call ncd_defvar(ncid=ncid, varname='time', xtype=nf_int,  &
            dim1name='time', &
            long_name='year', units='unitless')
     call ncd_defvar(ncid=ncid, varname='input_pftdata_filename', xtype=nf_char,  &
            dim1name='nchar', &
            dim2name='time',  &
            long_name='Input filepath for PFT values for this year', units='unitless')
  else
     call ncd_defvar(ncid=ncid, varname='time', xtype=nf_int,  &
            dim1name='time', &
            long_name='Calendar month', units='month')
  end if

end subroutine mkpftAtt

!-----------------------------------------------------------------------

end module mkpftMod
