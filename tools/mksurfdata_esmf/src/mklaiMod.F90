module mklaiMod

  !-----------------------------------------------------------------------
  ! Make LAI/SAI/height data
  !-----------------------------------------------------------------------

  use ESMF
  use pio
  use shr_kind_mod   , only : r8 => shr_kind_r8, r4=>shr_kind_r4
  use shr_sys_mod    , only : shr_sys_abort
  use mkpioMod       , only : mkpio_get_rawdata, mkpio_get_dimlengths
  use mkpioMod       , only : pio_iotype, pio_ioformat, pio_iosystem
  use mkesmfMod      , only : regrid_rawdata, create_routehandle_r8, get_meshareas
  use mkutilsMod     , only : chkerr
  use mkvarctl


  implicit none
  private

  public  :: mklai
  private :: pft_laicheck

!=================================================================================
contains
!=================================================================================

  subroutine mklai(ldomain, mapfname, datfname, ndiag, ncido)
    !
    ! Make LAI/SAI/height data
    ! Portions of this code could be moved out of the month loop
    ! for improved efficiency
    !
    use mkvarpar, only : re
    use mkpftConstantsMod, only : c3cropindex, c3irrcropindex
    !
    ! input/output variables
    character(len=*)  , intent(in)  :: file_mesh_i ! input mesh file name
    character(len=*)  , intent(in)  :: file_data_i ! input data file name
    type(ESMF_Mesh)   , intent(in)  :: mesh_o      ! output mesh
    integer           , intent(out) :: glacier_region_o(:) ! glacier region
    integer           , intent(out) :: rc
    !
    ! local variables
    integer                :: numpft_i          ! number of plant types on input
    type(ESMF_RouteHandle) :: routehandle       ! nearest neighbor routehandle
    type(ESMF_Mesh)        :: mesh_i
    type(file_desc_t)      :: pioid
    integer                :: ni,no,l 
    integer                :: ns_i, ns_o
    integer , allocatable  :: mask_i(:)
    real(r8), allocatable  :: frac_i(:)
    real(r8), allocatable  :: frac_o(:)
    real(r8), allocatable  :: data_i(:,:)
    real(r8), allocatable  :: data_o(:,:)
    real(r8), allocatable  :: mlai_o(:,:)       ! monthly lai
    real(r8), allocatable  :: msai_o(:,:)       ! monthly sai
    real(r8), allocatable  :: mhgtt_o(:,:)      ! monthly height (top)
    real(r8), allocatable  :: mhgtb_o(:,:)      ! monthly height (bottom)
    real(r8), allocatable  :: mlai_max(:,:)     ! monthly lai
    real(r8), allocatable  :: msai_max(:,:)     ! monthly sai
    real(r8), allocatable  :: mhgtt_max(:,:)    ! monthly height (top)
    real(r8), allocatable  :: mhgtb_max(:,:)    ! monthly height (bottom)
    real(r8), allocatable  :: mlai_i(:,:)       ! monthly lai in
    real(r8), allocatable  :: msai_i(:,:)       ! monthly sai in
    real(r8), allocatable  :: mhgtt_i(:,:)      ! monthly height (top) in
    real(r8), allocatable  :: mhgtb_i(:,:)      ! monthly height (bottom) in
    real(r8), allocatable  :: frac_dst(:)       ! output fractions: same as frac_dst
    integer,  allocatable  :: laimask(:,:)      ! lai+sai output mask for each plant function type
    real(r8)               :: glai_o(0:numpft)  ! output grid: global area pfts
    real(r8)               :: gsai_o(0:numpft)  ! output grid: global area pfts
    real(r8)               :: ghgtt_o(0:numpft) ! output grid: global area pfts
    real(r8)               :: ghgtb_o(0:numpft) ! output grid: global area pfts
    real(r8)               :: glai_i(0:numpft)  ! input grid: global area pfts
    real(r8)               :: gsai_i(0:numpft)  ! input grid: global area pfts
    real(r8)               :: ghgtt_i(0:numpft) ! input grid: global area pfts
    real(r8)               :: ghgtb_i(0:numpft) ! input grid: global area pfts
    real(r8)               :: garea_i           ! input  grid: global area
    real(r8)               :: garea_o           ! output grid: global area
    integer                :: mwts              ! number of weights
    integer                :: ni,no,ns_i,ns_o   ! indices
    integer                :: k,l,n,m           ! indices
    integer                :: dimids(4)         ! netCDF dimension ids
    integer                :: bego(4),leno(4)   ! netCDF bounds
    integer                :: begi(4),leni(4)   ! netCDF bounds 
    integer                :: ntim              ! number of input time samples
    integer                :: ier               ! error status
    real(r8)               :: relerr = 0.00001  ! max error: sum overlap wts ne 1
    character(len=256)     :: name              ! name of attribute
    character(len=256)     :: unit              ! units of attribute
    character(len= 32) :: subname = 'mklai'
    !-----------------------------------------------------------------------

    rc = ESMF_SUCCESS

    if (root_task) then
       write(ndiag,*)
       write (ndiag,'(a)') 'Attempting to make LAIs/SAIs/heights .....'
       write(ndiag,'(a)') ' Input file is '//trim(file_data_i)
    end if


    ! Open input data file
    call ESMF_VMLogMemInfo("Before pio_openfile for "//trim(file_data_i))
    rcode = pio_openfile(pio_iosystem, pioid, pio_iotype, trim(file_data_i), pio_nowrite)
    call ESMF_VMLogMemInfo("After pio_openfile "//trim(file_data_i))

    ! Read in input mesh
    call ESMF_VMLogMemInfo("Before create mesh_i in "//trim(subname))
    mesh_i = ESMF_MeshCreate(filename=trim(file_mesh_i), fileformat=ESMF_FILEFORMAT_ESMFMESH, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMLogMemInfo("After create mesh_i in "//trim(subname))

    ! Determine ns_i
    call ESMF_MeshGet(mesh_i, numOwnedElements=ns_i, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Determine ns_o
    call ESMF_MeshGet(mesh_o, numOwnedElements=ns_o, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Get the landmask from the input data file and reset the mesh mask based on that
    allocate(frac_i(ns_i), stat=ier)
    if (ier/=0) call shr_sys_abort()
    allocate(mask_i(ns_i), stat=ier)
    if (ier/=0) call shr_sys_abort()
    call mkpio_get_rawdata(pioid, 'LANDMASK', mesh_i, frac_i, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    do ni = 1,ns_i
       if (frac_i(ni) > 0._r8) then
          mask_i(ni) = 1
       else
          mask_i(ni) = 0
       end if
    end do
    call ESMF_MeshSet(mesh_i, elementMask=mask_i, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! Create a route handle between the input and output mesh
    allocate(frac_o(ns_o))
    call create_routehandle_r8(mesh_i, mesh_o, routehandle, frac_o=frac_o, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMLogMemInfo("After create routehandle in "//trim(subname))

    rcode = pio_inq_dimid(pioid 'pft', dimid)
    rcode = pio_inq_dimlen(pioid, dimid, numpft_i)
    rcode = pio_inq_dimid(pioid, 'time', dimid)
    rcode = pio_inq_dimlen(pioid, dimid, ntim)

    if (numpft_i /= numpft+1) then
       write(6,*) 'WARNING: ' // trim(subname) // '(): parameter numpft+1 = ', numpft+1, &
            'does not equal input dataset numpft = ', numpft_i
       write(6,*)'This inconsistency used to stop the program. Now we allow it '
       write(6,*)'because crop pfts 17-last are assumed to never use satellite lai data.'
       if (numpft_i > numpft + 1) then
          ! NOTE(bja, 2015-01) If this error check is determined to be
          ! invalid, all the loop bounds over output data in this
          ! routine will need to be double checked!
          write(6, *) "ERROR:" // trim(subname) // "(): input numpft must be less than or equal to output numpft+1."
          call shr_sys_abort()
       end if
    endif
    if (ntim /= 12) then
       write(6,*)'MKLAI: must have 12 time samples on input data'
       call shr_sys_abort()
    endif

    ! Dynamic allocation of variables of size 0:numpft_i
    allocate(mlai_i(ns_i,0:numpft_i),  &
             msai_i(ns_i,0:numpft_i),  &
             mhgtt_i(ns_i,0:numpft_i), &
             mhgtb_i(ns_i,0:numpft_i), stat=ier)
    if (ier /= 0) then
       write(6,*)'mklai allocation error'; call shr_sys_abort()
    end if

    ! Dynamic allocation of variables of size 0:numpft
    allocate(mlai_o(ns_o,0:numpft),  &
             msai_o(ns_o,0:numpft),  &
             mhgtt_o(ns_o,0:numpft), &
             mhgtb_o(ns_o,0:numpft), &
             laimask(ns_i,0:numpft), stat=ier )
    if (ier /= 0) then
       write(6,*)'mklai allocation error'; call shr_sys_abort()
    end if

    allocate(data_i(0:numpft_i,ns_i),stat=ier)
    if (ier/=0) call shr_sys_abort()

    ! Determine number of dimensions in output by querying MONTHLY_LAI


    ! Loop over months 
    mlai_o(:,:)  = 0.
    msai_o(:,:)  = 0.
    mhgtt_o(:,:) = 0.
    mhgtb_o(:,:) = 0.

    do nt = 1, ntim
       call mkpio_get_rawdata(pioid, 'MONTHLY_LAI', mesh_i, data_i, nt=nt, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call regrid_rawdata(mesh_i, mesh_o, routehandle, data_i, data_o, 0, numpft, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       do l = 0,numpft_i-1
          do no = 1,ns_o
             mlai_o(no,l) = data_o(l,no)
          end do
       end do

       call mkpio_get_rawdata(pioid, 'MONTHLY_SAI', mesh_i, data_i, nt=nt, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call regrid_rawdata(mesh_i, mesh_o, routehandle, data_i, data_o, 0, numpft, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       do l = 0,numpft_i-1
          do no = 1,ns_o
             lai_o(no,l) = data_o(l,no)
          end do
       end do

       call mkpio_get_rawdata(pioid, 'MONTHLY_HEIGHT_TOP', mesh_i, data_i, nt=nt, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call regrid_rawdata(mesh_i, mesh_o, routehandle, data_i, data_o, 0, numpft, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       do l = 0,numpft_i-1
          do no = 1,ns_o
             mhgtt_o(no,l) = data_o(l,no)
          end do
       end do

       call mkpio_get_rawdata(pioid, 'MONTHLY_HEIGHT_BOT', mesh_i, data_i, nt=nt, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call regrid_rawdata(mesh_i, mesh_o, routehandle, data_i, data_o, 0, numpft, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       do l = 0,numpft_i-1
          do no = 1,ns_o
             mhgtb_o(no,l) = data_o(l,no)
          end do
       end do

       ! copy LAI, SAI, & heights from the C3 crop (pft15)
       ! to the irrigated (pft16) whether crop is on or off
       mlai_o(:,c3irrcropindex)  = mlai_o(:,c3cropindex)
       msai_o(:,c3irrcropindex)  = msai_o(:,c3cropindex)
       mhgtt_o(:,c3irrcropindex) = mhgtt_o(:,c3cropindex)
       mhgtb_o(:,c3irrcropindex) = mhgtb_o(:,c3cropindex)

       ! Determine laimask
       laimask(:,:) = 0

       ! -----------------------------------------------------------------
       ! Output model resolution LAI/SAI/HEIGHT data
       ! -----------------------------------------------------------------

       if ( outnc_double ) then
          xtype = PIO_DOUBLE
       else
          xtype = PIO_REAL
       end if

       rcode = pio_redef(pioid_o)

       call mkpio_def_spatial_var(pioid_o, varname='MONTHLY_LAI', xtype=xtype,  &
          lev1name='lsmpft', lev2name='time', long_name='monthly leaf area index', units='unitless')

       call mkpio_def_spatial_var(pioid_o, varname='MONTHLY_SAI', xtype=xtype,  &
            lev1name='lsmpft', lev2name='time', long_name='monthly stem area index', units='unitless')

       call mkpio_def_spatial_var(pioid_o, varname='MONTHLY_HEIGHT_TOP', xtype=xtype,  &
          lev1name='lsmpft', lev2name='time',long_name='monthly height top', units='meters')

       call mkpio_def_spatial_var(pioid_o, varname='MONTHLY_HEIGHT_BOT', xtype=xtype,  &
            lev1name='lsmpft', lev2name='time', long_name='monthly height bottom', units='meters')

       ! End define model
       rcode = pio_enddef(pioid_o)

       ! Only need to define 1 PIO descriptor here
       call mkpio_iodesc_output(pioid_o, mesh_o, 'MONTHLY_LAI', pio_iodesc, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in generating an iodesc for MONTHLY_LAI')

       rcode = pio_inq_varid(pioid_o, 'MONTHLY_LAI', pio_varid)
       call pio_write_darray(pioid_o, pio_varid, pio_iodesc, mlai_o, rcode)

       rcode = pio_inq_varid(pioid_o, 'MONTHLY_SAI', pio_varid)
       call pio_write_darray(pioid_o, pio_varid, pio_iodesc, msai_o, rcode)

       rcode = pio_inq_varid(pioid_o, 'MONTHLY_HEIGHT_TOP', pio_varid)
       call pio_write_darray(pioid_o, pio_varid, pio_iodesc, mhgtt_o, rcode)

       rcode = pio_inq_varid(pioid_o, 'MONTHLY_HEIGHT_BOT', pio_varid)
       call pio_write_darray(pioid_o, pio_varid, pio_iodesc, mhgtt_o, rcode)

       ! TODO - output time - write out the wholE time array at once
       !rcode = pio_inq_varid(pioid_o, 'time', pio_varid)
       !rcode = pio_put_var(pioid_o, pio_varid, bego(ndimso), leno(ndimso), m)

       ! -----------------------------------------------------------------
       ! Compare global areas on input and output grids
       ! -----------------------------------------------------------------

       ! Input grid global area

       garea_i    = 0.
       do ni = 1,ns_i
          garea_i = garea_i + tgridmap%area_src(ni)
       end do

       glai_i(:)  = 0.
       gsai_i(:)  = 0.
       ghgtt_i(:) = 0.
       ghgtb_i(:) = 0.
       do l = 0, numpft_i - 1
          do ni = 1, ns_i
             glai_i(l)  = glai_i(l) + mlai_i(ni,l) *tgridmap%area_src(ni)*&
                  tdomain%mask(ni)*re**2
             gsai_i(l)  = gsai_i(l) + msai_i(ni,l) *tgridmap%area_src(ni)*&
                  tdomain%mask(ni)*re**2
             ghgtt_i(l) = ghgtt_i(l)+ mhgtt_i(ni,l)*tgridmap%area_src(ni)*&
                  tdomain%mask(ni)*re**2
             ghgtb_i(l) = ghgtb_i(l)+ mhgtb_i(ni,l)*tgridmap%area_src(ni)*&
                  tdomain%mask(ni)*re**2
          end do
       end do

       ! Output grid global area

       garea_o    = 0.
       do no = 1,ns_o
          garea_o = garea_o + tgridmap%area_dst(no)
       end do

       glai_o(:)  = 0.
       gsai_o(:)  = 0.
       ghgtt_o(:) = 0.
       ghgtb_o(:) = 0.
       do l = 0, numpft_i - 1
          do no = 1,ns_o
             glai_o(l)  = glai_o(l) + mlai_o(no,l)*tgridmap%area_dst(no)* &
                  frac_dst(no)*re**2
             gsai_o(l)  = gsai_o(l) + msai_o(no,l)*tgridmap%area_dst(no)* &
                  frac_dst(no)*re**2
             ghgtt_o(l) = ghgtt_o(l)+ mhgtt_o(no,l)*tgridmap%area_dst(no)* &
                  frac_dst(no)*re**2
             ghgtb_o(l) = ghgtb_o(l)+ mhgtb_o(no,l)*tgridmap%area_dst(no)* &
                  frac_dst(no)*re**2
          end do
       end do

       ! Comparison

       write (ndiag,*)
       write (ndiag,'(1x,70a1)') ('=',k=1,70)
       write (ndiag,*) 'LAI Output for month ',m
       write (ndiag,'(1x,70a1)') ('=',k=1,70)

       write (ndiag,*)
       write (ndiag,'(1x,70a1)') ('.',k=1,70)
       write (ndiag,1001)
1001   format (1x,'PFT input grid area output grid area',/ &
            1x,3x,'     10**6 km**2','      10**6 km**2')
       write (ndiag,'(1x,70a1)') ('.',k=1,70)
       write (ndiag,*)
       do l = 0, numpft
          write (ndiag,1002) l, glai_i(l)*1.e-06*1.e-02,glai_o(l)*1.e-06*1.e-02
1002      format (1x,i3,f16.3,f17.3)
       end do

       write (6,*) 'Successfully made LAIs/SAIs/heights for month ', m

    enddo
    write (6,*)

    ! Close input file
    call check_ret(nf_close(ncidi), subname)

    ! consistency check that PFT and LAI+SAI make sense
    !call pft_laicheck( ni_s, pft_i, laimask )

    ! Deallocate dynamic memory
    deallocate(mlai_i)
    deallocate(msai_i)
    deallocate(mhgtt_i)
    deallocate(mhgtb_i)
    deallocate(mlai_o)
    deallocate(msai_o)
    deallocate(mhgtt_o)
    deallocate(mhgtb_o)
    deallocate(laimask)
    deallocate(frac_dst)

    call gridmap_clean(tgridmap)
    call domain_clean(tdomain) 

  end subroutine mklai

  !-----------------------------------------------------------------------
  subroutine pft_laicheck( ni_s, pctpft_i, laimask )

    ! !USES:
    !
    ! !DESCRIPTION:
    !
    ! consistency check that PFT and LAI+SAI make sense
    !
    ! !ARGUMENTS:
    integer , intent(in) :: ni_s          ! input PFT grid resolution
    real(r8), intent(in) :: pctpft_i(:,:)  ! % plant function types
    integer,  intent(in) :: laimask(:,:)   ! mask where LAI+SAI > 0

    ! local variables
    character(len=*), parameter :: subName="pft_laicheck"
    integer :: ni,l,n,nc      ! Indices
    !-----------------------------------------------------------------------

    do l  = 0, numpft
       n  = 0
       nc = 0
       do ni = 1,ni_s
          if ( pctpft_i(ni,l) > 0.0_r8 ) nc = nc + 1
          if ( (pctpft_i(ni,l) > 0.0_r8) .and. (laimask(ni,l) /= 1) )then
             write (6,*) subName//' :: warning: pft and LAI+SAI mask not consistent!'
             write (6,*) 'ni,l   = ', ni, l
             write (6,*) 'pctpft_i  = ',pctpft_i(ni,l)
             write (6,*) 'laimask   = ', laimask(ni,l)
             n = n + 1
          end if
       end do
       if ( n > max(4,nc/4) ) then
          write (6,*) subName//' :: pft/LAI+SAI inconsistency over more than 25% land-cover'
          write (6,*) '# inconsistent points, total PFT pts, total LAI+SAI pts = ', &
               n, nc, sum(laimask(:,l))
          call abort()
       end if
    end do

  end subroutine pft_laicheck

end module mklaiMod
