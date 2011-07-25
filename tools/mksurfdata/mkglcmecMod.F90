module mkglcmecMod
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: mkglcmecMod
!
! !DESCRIPTION:
! Make glacier multi-elevation class  data
!
! !REVISION HISTORY:
! Author: Erik Kluzek
!
!-----------------------------------------------------------------------
!!USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use shr_sys_mod , only : shr_sys_flush
  implicit none

  private           ! By default make data private
!
! !PUBLIC MEMBER FUNCTIONS:
!
  public mkglcmecInit  ! Initialization
  public mkglcmec      ! Set glacier multi-elevation class
  public mkglacier     ! Set percent glacier
!
! !PUBLIC DATA MEMBERS: 
!
  integer, public       :: nglcec         =  0   ! number of elevation classes for glaciers
  real(r8), pointer     :: elevclass(:)          ! elevation classes
!EOP
!===============================================================
contains
!===============================================================

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: mkglcmecInit
!
! !INTERFACE:
subroutine mkglcmecInit( elevclass_o )
!
! !DESCRIPTION:
! Initialize of Make glacier multi-elevation class data
! !USES:
!
! !ARGUMENTS:
  implicit none
  real(r8), intent(OUT) :: elevclass_o(:)          ! elevation classes
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
  character(len=32) :: subname = 'mkglcmecInit:: '
!-----------------------------------------------------------------------
  allocate( elevclass(nglcec+1) )

  ! -----------------------------------------------------------------
  ! Define elevation classes, represents lower boundary of each class
  ! -----------------------------------------------------------------

  if (      nglcec == 10 )then
     elevclass(1)  =     0.
     elevclass(2)  =   200.
     elevclass(3)  =   400.
     elevclass(4)  =   700.
     elevclass(5)  =  1000.
     elevclass(6)  =  1300.
     elevclass(7)  =  1600.
     elevclass(8)  =  2000.
     elevclass(9)  =  2500.
     elevclass(10) =  3000.
     elevclass(11) = 10000.
  else if ( nglcec == 5  )then
     elevclass(1)  =     0.
     elevclass(2)  =   500.
     elevclass(3)  =  1000.
     elevclass(4)  =  1500.
     elevclass(5)  =  2000.
     elevclass(6)  = 10000.
  else if ( nglcec == 3  )then
     elevclass(1)  =     0.
     elevclass(2)  =  1000.
     elevclass(3)  =  2000.
     elevclass(4)  = 10000.
  else if ( nglcec == 1  )then
     elevclass(1)  =     0.
     elevclass(2)  = 10000.
  else if ( nglcec == 0  )then
     elevclass(1)  = 10000.
  else
     write(6,*) subname//"ERROR:: nglcec must be 0, 1, 3, 5, or 10 to work with CLM: "
     call abort()
  end if

  elevclass_o(:) = elevclass(:)

end subroutine mkglcmecInit

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: mkglcmec
!
! !INTERFACE:
subroutine mkglcmec(lsmlon, lsmlat, fname1, fname2, fname3, ndiag, pctglac_o, pctglcmec_o, &
                    topoglcmec_o, thckglcmec_o )
!
! !DESCRIPTION:
! make percent glacier on multiple elevation classes and mean elevation for each elevation class
!
! !USES:
  use domainMod   , only : domain_type,domain_clean,domain_setptrs
  use creategridMod, only : read_domain
  use mkvarsur    , only : ldomain
  use mkvarctl    
  use areaMod     , only : gridmap_type,gridmap_clean,gridmap_setptrs
  use areaMod     , only : areaini,areaave
  use ncdio
!
! !ARGUMENTS:
  implicit none
  integer , intent(in) :: lsmlon, lsmlat          ! clm grid resolution
  character(len=*), intent(in) :: fname1           ! input dataset file name
  character(len=*), intent(in) :: fname2           ! input dataset file name
  character(len=*), intent(in) :: fname3           ! input dataset file name
  integer , intent(in) :: ndiag                   ! unit number for diag out
  real(r8), intent(in) :: pctglac_o(lsmlon,lsmlat)      ! % glac (output grid)
  real(r8), intent(out) :: pctglcmec_o(lsmlon,lsmlat,nglcec)  ! 
  real(r8), intent(out) :: topoglcmec_o(lsmlon,lsmlat,nglcec) ! 
  real(r8), intent(out) :: thckglcmec_o(lsmlon,lsmlat,nglcec) ! 
!
! !CALLED FROM:
! subroutine mksrfdat in module mksrfdatMod
!
! !REVISION HISTORY:
! Author: David Lawrence
!
!
! !LOCAL VARIABLES:
!EOP
  integer  :: nlong_i                        ! input grid glacier: lon points
  integer  :: nlatg_i                        ! input grid glacier: lat points
  integer  :: nlont_i                        ! input grid topography: lon points
  integer  :: nlatt_i                        ! input grid topography: lat points

  type(domain_type)     :: tdomain1          ! local domain: topo_ice , topo_bedrock
  type(domain_type)     :: tdomain2          ! local domain: fracdata
  type(domain_type)     :: tdomain3          ! local domain: pct_glacier

  type(gridmap_type)    :: tgridmap         ! local gridmap
  integer          :: mxovr                 ! total num of overlapping cells
  integer ,pointer :: novr(:,:)             ! number of overlapping input cells
  integer ,pointer :: iovr(:,:,:)           ! lon index of overlap input cell
  integer ,pointer :: jovr(:,:,:)           ! lat index of overlap input cell
  real(r8),pointer :: wovr(:,:,:)           ! weight    of overlap input cell

  real(r8), allocatable :: pctglcmec_i(:,:,:)  ! % for each elevation class on input glacier grid
  real(r8), allocatable :: topoglcmec_i(:,:,:) ! mean elevation for each elevation classs on input glacier grid
  real(r8), allocatable :: thckglcmec_i(:,:,:) ! mean ice thickness for each elevation classs on input glacier grid
  real(r8), allocatable :: topoice_i(:,:)      ! topo of ice surface
  real(r8), allocatable :: topobedrock_i(:,:)  ! topo of bedrock surface
  real(r8), allocatable :: landmask_i(:,:)  ! input landmask for topography
  real(r8), allocatable :: glac_i(:,:)      ! input glacier fraction
  real(r8), allocatable :: maskt_i(:,:)     ! input grid topo: mask (0, 1)
  real(r8), allocatable :: maskg_i(:,:)     ! input grid glacier: mask (0, 1) 
  real(r8), allocatable :: mask_o(:,:)      ! output grid: mask (0, 1)
  integer  :: nptsec(nglcec)                ! number of points in an elevation class 
  real(r8) :: topoec(nglcec)                ! sum of all elevations, weighted by area in elevation class
  real(r8) :: thckec(nglcec)                ! sum of all elevations, weighted by area in elevation class
  real(r8) :: areaec(nglcec)                ! total area for all points within elevation class
  real(r8) :: areatot                       ! total area for glacier grid cell
  real(r8) :: pctareaec                     ! accumulator for pct area for elevation class
  real(r8) :: pctectot                      ! accumulator for total percentage area of elevation classes
  real(r8) :: pcttot                        ! accumulator for total percentage area in grid 
  real(r8) :: pcttotec(nglcec)              ! accumulator for total percentage area in grid 
  integer  :: nptsectot                     ! total number of points within glacier grid cell that are not ocean
  integer  :: po                            ! temporary flag

  integer  :: ii,ji                         ! indices
  integer  :: io,jo                         ! indices
  integer  :: k,l,n,m                       ! indices
  integer  :: ncid,dimid,varid              ! input netCDF id's
  integer  :: ier                           ! error status
  real(r8), parameter :: minglac = 1.e-6_r8 ! Minimum glacier amount
  character(len=32) :: subname = 'mkglcmec'
!-----------------------------------------------------------------------

  ! Initialize output to zero

  pctglcmec_o = 0.
  topoglcmec_o = 0.
  thckglcmec_o = 0.

  ! -----------------------------------------------------------------
  ! Exit early, if no glaciers exist or if nglcec = 0
  ! -----------------------------------------------------------------
  if ( all(pctglac_o < minglac ) )then
     write (6,*) 'No glaciers exist, set glcmec to zero as well'
     call shr_sys_flush(6)
     return
  end if

  if ( nglcec == 0 )then
     write (6,*) 'Number of glacier elevation classes is zero -- set glcmec to zero as well'
     call shr_sys_flush(6)
     return
  end if

  write (6,*) 'Attempting to make percent elevation class and mean elevation for glaciers .....'
  call shr_sys_flush(6)

  ! -----------------------------------------------------------------
  ! Read input file
  ! -----------------------------------------------------------------

  ! Obtain input grid info, read local fields

  call read_domain(tdomain1,fname1)
  call domain_setptrs(tdomain1,ni=nlont_i,nj=nlatt_i)

  call check_ret(nf_open(fname1, 0, ncid), subname)

  allocate(topoice_i(nlont_i, nlatt_i), topobedrock_i(nlont_i, nlatt_i), stat=ier)
  if (ier/=0) call abort()

  call check_ret(nf_inq_varid (ncid, 'TOPO_ICE', varid), subname)
  call check_ret(nf_get_var_double (ncid, varid, topoice_i), subname)
  call check_ret(nf_inq_varid (ncid, 'TOPO_BEDROCK', varid), subname)
  call check_ret(nf_get_var_double (ncid, varid, topobedrock_i), subname)

  call check_ret(nf_close(ncid), subname)

  call read_domain(tdomain2,fname2)
  call domain_setptrs(tdomain2,ni=nlont_i,nj=nlatt_i)

  call check_ret(nf_open(fname2, 0, ncid), subname)

  allocate(landmask_i(nlont_i, nlatt_i), stat=ier)                         
  if (ier/=0) call abort()

  call check_ret(nf_inq_varid (ncid, 'LANDMASK', varid), subname)    
  call check_ret(nf_get_var_double (ncid, varid, landmask_i), subname)   

  call check_ret(nf_close(ncid), subname)

  ! Get pct_glacier data

  call read_domain(tdomain3,fname3)
  call domain_setptrs(tdomain3,ni=nlong_i,nj=nlatg_i)

  call check_ret(nf_open(fname3, 0, ncid), subname)

  allocate(glac_i(nlong_i, nlatg_i), stat=ier)
  if (ier/=0) call abort()

  call check_ret(nf_inq_varid (ncid, 'PCT_GLACIER', varid), subname)
  call check_ret(nf_get_var_double (ncid, varid, glac_i), subname)

  call check_ret(nf_close(ncid), subname)

  ! First calculate elevation classes on input glacier grid

  allocate(maskt_i(nlont_i,nlatt_i),maskg_i(nlong_i,nlatg_i),mask_o(lsmlon,lsmlat))
  allocate(pctglcmec_i(nlong_i,nlatg_i,nglcec),topoglcmec_i(nlong_i,nlatg_i,nglcec),&
           thckglcmec_i(nlong_i,nlatg_i,nglcec))
  maskg_i = 1.0_r8
  maskt_i = 1.0_r8
  mask_o  = 1.0_r8

  topoglcmec_i = 0.
  thckglcmec_i = 0.
  pctglcmec_i = 0.

  call areaini(tdomain2,tdomain3,tgridmap,fracin=maskt_i,fracout=maskg_i)

  call gridmap_setptrs(tgridmap,mx_ovr=mxovr,n_ovr=novr,i_ovr=iovr,j_ovr=jovr,w_ovr=wovr)

  do jo = 1, nlatg_i
  do io = 1, nlong_i

     ! Initialize local variables
     nptsec = 0
     topoec = 0.
     thckec = 0.
     areaec = 0.
     areatot = 0.
     pctareaec = 0.
     pctectot = 0.
     po = 0
     nptsectot = 0
     
     if (glac_i(io,jo) > minglac) then 

        do n = 1, novr(io,jo)
           ii = iovr(io,jo,n)
           ji = jovr(io,jo,n)

           areatot = areatot + tdomain2%area(ii,ji)

           if (tdomain2%mask(ii,ji) .eq. 1) then
              do m = 1,nglcec
                 if (topoice_i(ii,ji) .ge. elevclass(m) .and. topoice_i(ii,ji) .lt. elevclass(m+1)) then
                    nptsec(m) = nptsec(m)+1
                    nptsectot = nptsectot + nptsec(m)
                    topoec(m) = topoec(m) + wovr(io,jo,n)*tdomain2%area(ii,ji)*topoice_i(ii,ji)

                    ! bedrock cannot be below mean sea level; required to avoid overly thick ice sheets 
                    ! in ice shelf terrain (?)
                    if ( (topoice_i(ii,ji) - topobedrock_i(ii,ji)) .gt. elevclass(m+1) ) then
                       topobedrock_i(ii,ji) = 0
                    endif
                    thckec(m) = thckec(m) + wovr(io,jo,n)*tdomain2%area(ii,ji)*&
                                (topoice_i(ii,ji)-topobedrock_i(ii,ji))
                    areaec(m) = areaec(m) + wovr(io,jo,n)*tdomain2%area(ii,ji)
                 endif
              enddo
           endif
 
        enddo

        do m = nglcec,1,-1
           pctareaec = pctareaec + 100.*areaec(m)/areatot
           if (pctareaec .le. glac_i(io,jo) .and. areaec(m) .gt. 0) then 
              pctglcmec_i(io,jo,m) = 100.*areaec(m)/areatot
              topoglcmec_i(io,jo,m) = topoec(m)/areaec(m)
              thckglcmec_i(io,jo,m) = thckec(m)/areaec(m)
!              if (thckglcmec_i(io,jo,m) .gt. elevclass(m+1)+400) then 
!                 write(6,*) 'TOPO ',io,jo,m,topoglcmec_i(io,jo,m),thckglcmec_i(io,jo,m)
!              endif
           elseif (pctareaec .gt. glac_i(io,jo) .and. po .eq. 0) then 
              pctglcmec_i(io,jo,m) = 100.*areaec(m)/areatot
              topoglcmec_i(io,jo,m) = topoec(m)/areaec(m)
              thckglcmec_i(io,jo,m) = thckec(m)/areaec(m)
              po = 1
           endif

           ! if all topo points within a glacier grid point are zero, then glacier is ice 
           ! shelf with an assumed elevation of 5m
           if (nptsectot .eq. 0) then
              pctglcmec_i(io,jo,1) = 100.
              topoglcmec_i(io,jo,1) = 5
              thckglcmec_i(io,jo,1) = 5
           endif
           pctectot = pctectot + pctglcmec_i(io,jo,m)
        enddo

        ! Error check: are all elevations within elevation class range

        do m = 1,nglcec
           if ( (topoglcmec_i(io,jo,m) .lt. elevclass(m) .or. topoglcmec_i(io,jo,m) .gt. elevclass(m+1)) &
              .and. topoglcmec_i(io,jo,m) .ne. 0) then
              write(6,*) 'Warning: mean elevation does not fall within elevation class '
              write(6,*) elevclass(m),elevclass(m+1),topoglcmec_i(io,jo,m),m,io,jo  
           endif
        enddo

        ! normalize so that sum of pctglcmec_i adds up to one
        do m = 1,nglcec
           if ( pctectot /= 0.0_r8 ) &
           pctglcmec_i(io,jo,m) = (100./pctectot)*pctglcmec_i(io,jo,m)
        enddo
     endif
     
  enddo 
  enddo

!! Average from input pct_glacier to output grid

  call areaini(tdomain3,ldomain,tgridmap,fracin=maskg_i,fracout=mask_o)

  call gridmap_setptrs(tgridmap,mx_ovr=mxovr,n_ovr=novr,i_ovr=iovr,j_ovr=jovr,w_ovr=wovr)

  do jo = 1, ldomain%nj
  do io = 1, ldomain%numlon(jo)

     pcttot = 0.
     pcttotec = 0.
     if (pctglac_o(io,jo) .gt. minglac) then 

        do n = 1, novr(io,jo)
           ii = iovr(io,jo,n)
           ji = jovr(io,jo,n)
           pcttot = pcttot + wovr(io,jo,n)*glac_i(ii,ji)                                                 
           do m = 1,nglcec
              pcttotec(m) = pcttotec(m) + wovr(io,jo,n)*glac_i(ii,ji)*pctglcmec_i(ii,ji,m)
           enddo 
        enddo

        do n = 1, novr(io,jo)
           ii = iovr(io,jo,n)
           ji = jovr(io,jo,n)

           do m = 1,nglcec
              pctglcmec_o(io,jo,m) = pctglcmec_o(io,jo,m) + &
                                     wovr(io,jo,n)*glac_i(ii,ji)*pctglcmec_i(ii,ji,m)/pcttot
              if (pcttotec(m) .ne. 0) then
                 topoglcmec_o(io,jo,m) = topoglcmec_o(io,jo,m) + &
                    wovr(io,jo,n)*pctglcmec_i(ii,ji,m)*glac_i(ii,ji)*topoglcmec_i(ii,ji,m)/pcttotec(m)
                 thckglcmec_o(io,jo,m) = thckglcmec_o(io,jo,m) + &
                    wovr(io,jo,n)*pctglcmec_i(ii,ji,m)*glac_i(ii,ji)*thckglcmec_i(ii,ji,m)/pcttotec(m)
              endif
           enddo
        enddo

        ! Scale according to grid cell pct_glacier
        do m = 1,nglcec
           pctglcmec_o(io,jo,m) = pctglac_o(io,jo)*pctglcmec_o(io,jo,m)/100._r8
        enddo

        ! Error check: are all elevations within elevation class range
        do m = 1,nglcec
           if ( (topoglcmec_o(io,jo,m) .lt. elevclass(m) .or. topoglcmec_o(io,jo,m) .gt. elevclass(m+1)) &
                 .and. topoglcmec_o(io,jo,m) .ne. 0) then
              write(6,*) 'Warning: mean elevation does not fall within elevation class '
              write(6,*) elevclass(m),elevclass(m+1),topoglcmec_o(io,jo,m),m,io,jo
           endif
        enddo
     endif
  enddo
  enddo

  write (6,*) 'Successfully made percent elevation class and mean elevation for glaciers'
  write (6,*)
  call shr_sys_flush(6)

  ! Deallocate dynamic memory

  call domain_clean(tdomain1)
  call domain_clean(tdomain2)
  call domain_clean(tdomain3)
  call gridmap_clean(tgridmap)
  deallocate (topoice_i, topobedrock_i, landmask_i, glac_i)
  deallocate (maskt_i, maskg_i, mask_o)
  deallocate (pctglcmec_i, topoglcmec_i, thckglcmec_i)

end subroutine mkglcmec

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: mkglacier
!
! !INTERFACE:
subroutine mkglacier(lsmlon, lsmlat, fname, ndiag, zero_out, glac_o )
!
! !DESCRIPTION:
! make percent glacier
!
! !USES:
  use domainMod   , only : domain_type,domain_clean,domain_setptrs
  use creategridMod, only : read_domain
  use mkvarpar	
  use mkvarsur    , only : ldomain
  use mkvarctl    
  use areaMod     , only : areaini,areaave,gridmap_type,gridmap_clean  
  use ncdio
!
! !ARGUMENTS:
  implicit none
  integer , intent(in) :: lsmlon, lsmlat          ! clm grid resolution
  character(len=*), intent(in) :: fname           ! input dataset file name
  logical , intent(in) :: zero_out                ! if should zero glacier out
  integer , intent(in) :: ndiag                   ! unit number for diag out
  real(r8), intent(out):: glac_o(lsmlon,lsmlat)   ! output grid: %glacier
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
  integer  :: nlon_i                          ! input grid : lon points
  integer  :: nlat_i                          ! input grid : lat points

  type(domain_type)     :: tdomain            ! local domain
  type(gridmap_type)    :: tgridmap           ! local gridmap

  real(r8), allocatable :: glac_i(:,:)        ! input grid: percent glac
  real(r8), allocatable :: mask_i(:,:)        ! input grid: mask (0, 1)
  real(r8), allocatable :: mask_o(:,:)        ! output grid: mask (0, 1)
  real(r8), allocatable :: fld_i(:,:)         ! input grid: dummy field
  real(r8), allocatable :: fld_o(:,:)         ! output grid: dummy field
  real(r8) :: sum_fldi                        ! global sum of dummy input fld
  real(r8) :: sum_fldo                        ! global sum of dummy output fld
  real(r8) :: gglac_i                          ! input  grid: global glac
  real(r8) :: garea_i                         ! input  grid: global area
  real(r8) :: gglac_o                          ! output grid: global glac
  real(r8) :: garea_o                         ! output grid: global area

  integer  :: ii,ji                           ! indices
  integer  :: io,jo                           ! indices
  integer  :: k,n,m                           ! indices
  integer  :: ncid,dimid,varid                ! input netCDF id's
  integer  :: ier                             ! error status
  real(r8) :: relerr = 0.00001                ! max error: sum overlap wts ne 1
  character(len=32) :: subname = 'mkglacier'
!-----------------------------------------------------------------------

  write (6,*) 'Attempting to make %glacier .....'
  call shr_sys_flush(6)

  ! -----------------------------------------------------------------
  ! Read input file
  ! -----------------------------------------------------------------

  ! Obtain input grid info, read local fields

  call read_domain(tdomain,fname)
  call domain_setptrs(tdomain,ni=nlon_i,nj=nlat_i)

  call check_ret(nf_open(fname, 0, ncid), subname)

  allocate(glac_i(nlon_i,nlat_i), stat=ier)
  if (ier/=0) call abort()

  call check_ret(nf_inq_varid (ncid, 'PCT_GLACIER', varid), subname)
  call check_ret(nf_get_var_double (ncid, varid, glac_i), subname)

  call check_ret(nf_close(ncid), subname)

  ! Compute local fields _o

  allocate(mask_i(nlon_i,nlat_i),mask_o(lsmlon,lsmlat))
  allocate( fld_i(nlon_i,nlat_i), fld_o(lsmlon,lsmlat))

  mask_i = 1.0_r8
  mask_o = 1.0_r8
  call areaini(tdomain,ldomain,tgridmap,fracin=mask_i,fracout=mask_o)

  if ( zero_out )then

     do jo = 1, ldomain%nj
     do io = 1, ldomain%numlon(jo)
        glac_o(io,jo) = 0.
     enddo
     enddo

  else

     mask_i = float(tdomain%mask(:,:))
     call areaave(mask_i,mask_o,tgridmap)

     call gridmap_clean(tgridmap)
     call areaini(tdomain,ldomain,tgridmap,fracin=mask_i,fracout=mask_o)

     ! Area-average percent cover on input grid to output grid 
     ! and correct according to land landmask
     ! Note that percent cover is in terms of total grid area.

     call areaave(glac_i,glac_o,tgridmap)

     do jo = 1, ldomain%nj
     do io = 1, ldomain%numlon(jo)
        if (glac_o(io,jo) < 1.) glac_o(io,jo) = 0.
     enddo
     enddo

  end if

  ! Check for conservation

  do jo = 1, ldomain%nj
  do io = 1, ldomain%numlon(jo)
     if ((glac_o(io,jo)) > 100.000001_r8) then
        write (6,*) 'MKGLACIER error: glacier = ',glac_o(io,jo), &
                ' greater than 100.000001 for column, row = ',io,jo
        call shr_sys_flush(6)
        stop
     end if
  enddo
  enddo

  ! Some error checking and writing of global values before and after the regrid

  if ( .not. zero_out )then

     ! Global sum of output field -- must multiply by fraction of
     ! output grid that is land as determined by input grid

     sum_fldi = 0.0_r8
     do ji = 1, nlat_i
     do ii = 1, nlon_i
       fld_i(ii,ji) = ((ji-1)*nlon_i + ii) * tdomain%mask(ii,ji)
       sum_fldi = sum_fldi + tdomain%area(ii,ji) * fld_i(ii,ji)
     enddo
     enddo

     call areaave(fld_i,fld_o,tgridmap)

     sum_fldo = 0.
     do jo = 1, ldomain%nj
     do io = 1, ldomain%numlon(jo)
        fld_o(io,jo) = fld_o(io,jo)*mask_o(io,jo)
        sum_fldo = sum_fldo + ldomain%area(io,jo) * fld_o(io,jo)
     end do
     end do

     ! -----------------------------------------------------------------
     ! Error check1
     ! Compare global sum fld_o to global sum fld_i.
     ! -----------------------------------------------------------------

     if ( trim(mksrf_gridtype) == 'global') then
        if ( abs(sum_fldo/sum_fldi-1.) > relerr ) then
           write (6,*) 'MKGLACIER error: input field not conserved'
           write (6,'(a30,e20.10)') 'global sum output field = ',sum_fldo
           write (6,'(a30,e20.10)') 'global sum input  field = ',sum_fldi
           stop
        end if
     end if

     ! -----------------------------------------------------------------
     ! Error check2
     ! Compare global areas on input and output grids
     ! -----------------------------------------------------------------

     ! Input grid

     gglac_i = 0.
     garea_i = 0.

     do ji = 1, nlat_i
     do ii = 1, nlon_i
        garea_i = garea_i + tdomain%area(ii,ji)
        gglac_i = gglac_i + glac_i(ii,ji)*(tdomain%area(ii,ji)/100.) * &
                            tdomain%frac(ii,ji)
     end do
     end do

     ! Output grid

     gglac_o = 0.
     garea_o = 0.

     do jo = 1, ldomain%nj
     do io = 1, ldomain%numlon(jo)
        garea_o = garea_o + ldomain%area(io,jo)
        gglac_o = gglac_o + glac_o(io,jo)*(ldomain%area(io,jo)/100.) * &
                            ldomain%frac(io,jo)
     end do
     end do

     ! Diagnostic output

     write (ndiag,*)
     write (ndiag,'(1x,70a1)') ('=',k=1,70)
     write (ndiag,*) 'Glacier Output'
     write (ndiag,'(1x,70a1)') ('=',k=1,70)

     write (ndiag,*)
     write (ndiag,'(1x,70a1)') ('.',k=1,70)
     write (ndiag,2001)
2001 format (1x,'surface type   input grid area  output grid area'/ &
             1x,'                 10**6 km**2      10**6 km**2   ')
     write (ndiag,'(1x,70a1)') ('.',k=1,70)
     write (ndiag,*)
     write (ndiag,2002) gglac_i*1.e-06,gglac_o*1.e-06
     write (ndiag,2004) garea_i*1.e-06,garea_o*1.e-06
2002 format (1x,'glaciers    ',f14.3,f17.3)
2004 format (1x,'all surface ',f14.3,f17.3)

     if (lsmlat > 1) then
        k = lsmlat/2
        write (ndiag,*)
        write (ndiag,*) 'For reference the area on the output grid of a cell near the equator is: '
        write (ndiag,'(f10.3,a14)')ldomain%area(1,k)*1.e-06,' x 10**6 km**2'
        write (ndiag,*)
     endif
     call shr_sys_flush(ndiag)

  end if

  write (6,*) 'Successfully made %glacier'
  write (6,*)
  call shr_sys_flush(6)

  ! Deallocate dynamic memory

  call domain_clean(tdomain)
  call gridmap_clean(tgridmap)
  deallocate (glac_i)
  deallocate (mask_i,mask_o)
  deallocate ( fld_i, fld_o)

end subroutine mkglacier

!-----------------------------------------------------------------------

end module mkglcmecMod
