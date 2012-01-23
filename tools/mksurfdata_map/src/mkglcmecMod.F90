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
! Author: Erik Kluzek, Mariana Vertenstein
!
!-----------------------------------------------------------------------
!!USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use shr_sys_mod , only : shr_sys_flush
  use mkdomainMod , only : domain_checksame
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
  integer, public       :: nglcec         = 10   ! number of elevation classes for glaciers
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
     write(6,*) subname//"ERROR:: nglcec must be 0, 1, 3, 5, or 10",&
          " to work with CLM: "
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
subroutine mkglcmec(ldomain, mapfname, &
                    datfname_fglctopo, datfname_fglacier, ndiag, &
                    pctglac_o, pctglac_o_uncorrected, &
                    pctglcmec_o, topoglcmec_o, thckglcmec_o, &
                    pctglc_gic_o, pctglc_icesheet_o, pctglc_float_o, pctglc_ground_o)
!
! !DESCRIPTION:
! make percent glacier on multiple elevation classes and mean elevation for each elevation class
!
! Note: the input glacier and topo datasets (datfname_fglctopo, datfname_fglacier) must
! agree in resolution, and also must have identical longitudes and latitudes (e.g.,
! agreement in ordering and convention of longitudes). However, this requirement is
! removed if nglcec==0, in which case this subroutine doesn't do anything. (So, for
! example, it is okay for pct_glacier to be given at 1/2 deg and topo at 10' resolution if
! nglcec==0.)
!
! Also note that the topo dataset is regridded here using the same mapping file as is used
! for regridding the input glacier data, so that topo_ice & pct_glacier data are regridded
! in the same way. This mapping file may differ (e.g., in its use of a landmask) from the
! mapping file used to regrid the topo data elsewhere in mksurfdat.
!
! !USES:
  use mkdomainMod, only : domain_type, domain_clean, domain_read
  use mkgridmapMod
  use mkvarpar	
  use mkncdio
!
! !ARGUMENTS:
  implicit none
  type(domain_type) , intent(in) :: ldomain
  character(len=*)  , intent(in) :: mapfname           ! input mapping file name
  character(len=*)  , intent(in) :: datfname_fglctopo  ! raw glc topo data
  character(len=*)  , intent(in) :: datfname_fglacier  ! raw glacier data (same resolution as fglctopo)
  integer           , intent(in) :: ndiag              ! unit number for diag out
  real(r8)          , intent(in) :: pctglac_o(:)       ! % glac on output glacier grid
  real(r8)          , intent(in) :: pctglac_o_uncorrected(:) ! % glac on output glacier grid, before any corrections were done
  real(r8)          , intent(out):: pctglcmec_o (:,:)  ! % for each elevation class on output glacier
  real(r8)          , intent(out):: topoglcmec_o(:,:)  ! mean elevation for each elevation classs on output glacier 
  real(r8)          , intent(out):: thckglcmec_o(:,:)  ! mean ice thickness for each elevation classs on input glacier 
  real(r8)          , intent(out):: pctglc_gic_o(:)       ! % glc gic on output grid
  real(r8)          , intent(out):: pctglc_icesheet_o(:)  ! % glc ice sheet on output grid
  real(r8)          , intent(out):: pctglc_float_o(:)     ! % glc float on output grid
  real(r8)          , intent(out):: pctglc_ground_o(:)    ! % glc ground on output grid
!
! !CALLED FROM:
! subroutine mksrfdat in module mksrfdatMod
!
! !REVISION HISTORY:
! Author: David Lawrence
! 7/12/11: Bill Sacks: substantial rewrite to use input topo and % glacier at same resolution
!
!
! !LOCAL VARIABLES:
!EOP
  type(domain_type)     :: tdomain_topo        ! local domain: topo_ice , topo_bedrock
  type(gridmap_type)    :: tgridmap            ! local gridmap
  real(r8), allocatable :: pctglac_g(:)        ! input glacier percentage (on input glacier grid)
  real(r8), allocatable :: topoice_g(:)        ! topo of ice surface (on input topo grid - same as input glacier grid)
  real(r8), allocatable :: pctglc_gic_g(:)     ! % glc gic on input topo/glacier grid
  real(r8), allocatable :: pctglc_icesheet_g(:)! % glc ice sheet on input topo/glacier grid
  real(r8), allocatable :: pctglc_ground_g(:)  ! % glc floating ice on input topo/glacier grid
  real(r8), allocatable :: pctglc_float_g(:)   ! % glc grounded ice on input topo/glacier grid
  real(r8), allocatable :: topoglcmec_unnorm_o(:,:) ! same as topoglcmec_o, but unnormalized
  real(r8) :: wt, frac                         ! weighting factors for remapping
  integer  :: ni,no,ns_o,nst                   ! indices
  integer  :: m,n                              ! indices
  integer  :: ncid,dimid,varid                 ! input netCDF id's
  real(r8) :: sum                              ! temporary
  integer  :: ier                              ! error status
  logical  :: errors                           ! error status
  real(r8), parameter :: minglac = 1.e-6_r8    ! Minimum glacier amount
  character(len=32) :: subname = 'mkglcmec'
!-----------------------------------------------------------------------

  ! Initialize output to zero

  pctglcmec_o  = 0.
  topoglcmec_o = 0.
  thckglcmec_o = 1.e36

  ns_o = ldomain%ns

  ! -----------------------------------------------------------------
  ! Exit early, if no glaciers exist
  ! -----------------------------------------------------------------
  if ( all(pctglac_o < minglac ) )then
     write (6,*) 'No glaciers exist, set glcmec to zero as well'
     return
  end if

  if ( nglcec == 0 )then
     write (6,*) 'Number of glacier elevation classes is zero ',&
          '-- set glcmec to zero as well'
     call shr_sys_flush(6)
     return
  end if

  write (6,*) 'Attempting to make percent elevation class ',&
       'and mean elevation for glaciers .....'
  call shr_sys_flush(6)

  ! -----------------------------------------------------------------
  ! Read input files
  ! -----------------------------------------------------------------

  ! Get raw topo data 

  call domain_read(tdomain_topo,datfname_fglctopo)
  nst = tdomain_topo%ns

  allocate(topoice_g(nst), stat=ier)
  if (ier/=0) call abort()

  allocate(pctglc_gic_g(nst), pctglc_icesheet_g(nst), &
       pctglc_ground_g(nst), pctglc_float_g(nst), stat=ier)
  if (ier/=0) call abort()

  write (6,*) 'Open glacier topo file: ', trim(datfname_fglctopo)
  call check_ret(nf_open(datfname_fglctopo, 0, ncid), subname)
  call check_ret(nf_inq_varid (ncid, 'TOPO_ICE', varid), subname)
  call check_ret(nf_get_var_double (ncid, varid, topoice_g), subname)

  ier = nf_inq_varid (ncid, 'pct_glc_gic', varid)
  if (ier == NF_NOERR) then
     write(6,*)'reading pct_glc_gic data'
     call check_ret(nf_get_var_double (ncid, varid, pctglc_gic_g), subname)
  else
     pctglc_gic_g(:) = 0
  end if

  ier = nf_inq_varid (ncid, 'pct_glc_icesheet', varid)
  if (ier == NF_NOERR) then
     write(6,*)'reading pct_glc_icesheet data'
     call check_ret(nf_get_var_double (ncid, varid, pctglc_icesheet_g), subname)
  else
     pctglc_icesheet_g(:) = 0
  end if

  ier = nf_inq_varid (ncid, 'pct_glc_ground', varid)
  if (ier == NF_NOERR) then
     write(6,*)'reading pct_glc_ground data'
     call check_ret(nf_get_var_double (ncid, varid, pctglc_ground_g), subname)
  else
     pctglc_ground_g(:) = 0
  end if

  ier = nf_inq_varid (ncid, 'pct_glc_float', varid)
  if (ier == NF_NOERR) then
     write(6,*)'reading pct_glc_float data'
     call check_ret(nf_get_var_double (ncid, varid, pctglc_float_g), subname)
  else
     pctglc_float_g(:) = 0
  end if

  call check_ret(nf_close(ncid), subname)

  ! Get raw glacier data 

  call domain_read(tdomain_topo, datfname_fglacier)

  allocate(pctglac_g(nst), stat=ier)
  if (ier/=0) call abort()

  write (6,*) 'Open glacier file: ', trim(datfname_fglacier)
  call check_ret(nf_open(datfname_fglacier, 0, ncid), subname)
  call check_ret(nf_inq_varid (ncid, 'PCT_GLACIER', varid), subname)
  call check_ret(nf_get_var_double (ncid, varid, pctglac_g), subname)
  call check_ret(nf_close(ncid), subname)

  ! Mapping for raw glacier or topo -> model output grid

  call gridmap_mapread(tgridmap, mapfname )

  ! Error checks for domain and map consistencies: ensure that both the topo and glacier
  ! domains are consistent with the given mapping file
  call domain_checksame( tdomain_topo, ldomain, tgridmap )

  ! -------------------------------------------------------------------- 
  ! Compute fields on the output grid
  ! -------------------------------------------------------------------- 

  allocate(topoglcmec_unnorm_o(ns_o,nglcec), stat=ier)
  if (ier/=0) call abort()

  pctglcmec_o(:,:)  = 0.
  topoglcmec_unnorm_o(:,:) = 0.
  topoglcmec_o(:,:) = 0.
     
  do n = 1,tgridmap%ns
     ni = tgridmap%src_indx(n)
     no = tgridmap%dst_indx(n)
     wt = tgridmap%wovr(n)

     ! fraction of this destination cell that is covered by source cells that are within the source landmask
     frac = tgridmap%frac_dst(no)

     ! We don't bother with the following if pctglac_o(no) <= minglac, because the
     ! pctglcmec_o values will end up being <= minglac.
     ! Also, if frac == 0, then we can't do this, to avoid divide by 0. In this case, the
     ! outputs remain equal to 0 (their initialized value).
     if (pctglac_o(no) > minglac .and. frac > 0) then
        m = get_elevclass(topoice_g(ni))
        if (m < 1 .or. m > nglcec) then 
           call abort()
        end if

        pctglcmec_o(no,m) = pctglcmec_o(no,m) + wt*pctglac_g(ni) / frac

        ! note that, by weighting the following by pctglac_g, we are getting something
        ! like the average topographic height over glaciated areas - NOT the average
        ! topographic height of the entire grid cell
        topoglcmec_unnorm_o(no,m) = topoglcmec_unnorm_o(no,m) + wt*pctglac_g(ni)*topoice_g(ni) / frac
     end if
  end do

  ! Now pctglcmec_o is already normalized appropriately, because for destarea
  ! normalization, the sum of the wovr values for a given destination grid cell should
  ! equal frac_dst, so sum(wt/frac) = 1 for each destination grid cell.
  !
  ! However, we need to normalize topoglcmec_o. To do this, note that pctglcmec_o(n,m) is
  ! equal to the sum of the weights used in doing the weighted average of topoice_g
  ! (weight = wt*pctglac_g(ni)/frac); hence pctglcmec_o(n,m) is the correct normalization
  ! factor
  do no = 1,ns_o
     do m = 1,nglcec
        if (pctglcmec_o(no,m) > 0) then
           topoglcmec_o(no,m) = topoglcmec_unnorm_o(no,m) / pctglcmec_o(no,m)
        end if
     end do
  end do

  ! Map other fields

  call gridmap_areaave(tgridmap, pctglc_gic_g     , pctglc_gic_o)
  call gridmap_areaave(tgridmap, pctglc_icesheet_g, pctglc_icesheet_o)
  call gridmap_areaave(tgridmap, pctglc_float_g   , pctglc_float_o)
  call gridmap_areaave(tgridmap, pctglc_ground_g  , pctglc_ground_o)

  ! We want the sum of pctglcmec_o across elevation classes to equal pctglac_o for each
  ! grid cell. This should be approximately true at this point, but there have been a
  ! number of adjustments to pctglac_o since its original remapping (e.g., setting values
  ! < 1 to 0, and rescaling all percentages so that the total landcover is equal to
  ! 100%). We need to essentially apply these same adjustments to pctglcmec_o.
  ! 
  ! To do this, we use the stored, uncorrected values of pctglac_o. In theory, we could do
  ! this renormalization without reference to these uncorrected values, since
  ! pctglac_o_uncorrected should equal the sum of pctglcmec_o across elevation classes for
  ! each grid cell. However, pctglac_o_uncorrected is helpful for error checking: it
  ! allows us to verify that the pctglcmec sums are correct (this check is done below).
  ! 
  ! Here we also rescale the pctglc_gic, pctglc_icesheet, etc. values similarly, to
  ! achieve the same thing (i.e., we want sums to add up to pct_glacier)
  do no = 1,ns_o
     if (pctglac_o_uncorrected(no) > 0) then
        pctglcmec_o(no,:) = pctglcmec_o(no,:) * (pctglac_o(no) / pctglac_o_uncorrected(no))
        pctglc_gic_o(no) = pctglc_gic_o(no) * (pctglac_o(no) / pctglac_o_uncorrected(no))
        pctglc_icesheet_o(no) = pctglc_icesheet_o(no) * (pctglac_o(no) / pctglac_o_uncorrected(no))
        pctglc_float_o(no) = pctglc_float_o(no) * (pctglac_o(no) / pctglac_o_uncorrected(no))
        pctglc_ground_o(no) = pctglc_ground_o(no) * (pctglac_o(no) / pctglac_o_uncorrected(no))

     else if (pctglac_o_uncorrected(no) == 0 .and. pctglac_o(no) > 0) then
        ! There isn't a clear way to handle the case when pctglac_o_uncorrected==0 (and
        ! hence pctglcmec_o==0 for all elevation classes) but pctglac_o==0. Fortunately,
        ! this should never happen: the only place where there is a potential for this is
        ! over the Ross ice shelf and the south pole, where pctglac_o values are set to
        ! 100%. However that isn't a problem as long as ice shelves are handled in the
        ! input data by (1) being called land by the landmask, (2) having some non-zero
        ! pct_glacier value, and (3) having some real topo value (i.e., not just 0 or
        ! missing value everywhere); these conditions are satisfied as of 7-12-11. In that
        ! case, ice shelves should be treated properly - i.e., treated like other
        ! glaciers. However, if those conditions weren't satisfied, then this check could
        ! fail for some points, and we would have to handle ice shelves specially (e.g.,
        ! in the old code, points made up entirely of ocean were treated as ice shelves
        ! with topo = 5m).

        write(6,*) 'ERROR in ', subname
        write(6,*) 'pctglac_o_uncorrected==0 but pctglac_o > 0'
        write(6,*) 'no = ', no
        write(6,*) 'lon = ', tgridmap%xc_dst(no)
        write(6,*) 'lat = ', tgridmap%yc_dst(no)
        call abort()
     end if

     ! The only other possibility is that (pctglac_o_uncorrected(no) == 0 .and.
     ! pctglac_o(no) == 0); but in this case, all pctglcmec_o(no,:) values should be 0
     ! and can remain 0 (and similarly for pctglc_gic, pctglc_icesheet, etc.)
  end do

  errors = .false.

  ! Check that sum over pctglcmec_o (from 1 to nglcec) is equal to pctglac_o(no)
  do no = 1,ns_o
     sum = 0.
     do m = 1,nglcec
        sum = sum + pctglcmec_o(no,m)
     end do
     if (abs(sum - pctglac_o(no)) > 1.e-6) then
        write(6,*)'no,pctglc,pctglac= ',no,sum,pctglac_o(no)
        errors = .true.
     end if
  end do

  ! Check that sums of various pctglc_* variables equal pctglac_o
  ! Use threshold of 1e-5, because 1e-6 is occasionally violated (I think because
  ! pct_glacier is given as a float whereas the pctglc_* values are given as doubles on
  ! the input datasets)
  do no = 1,ns_o
     if (abs((pctglc_gic_o(no) + pctglc_icesheet_o(no)) - pctglac_o(no)) > 1.e-5) then
        write(6,*)'no,pctglc_gic,pctglc_icesheet,pctglac,pctglac_uncorrected,lon,lat=',no,pctglc_gic_o(no),&
             pctglc_icesheet_o(no),pctglac_o(no),pctglac_o_uncorrected(no),tgridmap%&
             xc_dst(no),tgridmap%yc_dst(no)
        errors = .true.
     end if

     if (abs((pctglc_ground_o(no) + pctglc_float_o(no)) - pctglac_o(no)) > 1.e-5) then
        write(6,*)'no,pctglc_ground,pctglc_float,pctglac,pctglac_uncorrected,lon,lat=',no,pctglc_ground_o(no),&
             pctglc_float_o(no),pctglac_o(no),pctglac_o_uncorrected(no),tgridmap%&
             xc_dst(no),tgridmap%yc_dst(no)
        errors = .true.
     end if
  end do

  ! Error check: are all elevations within elevation class range
  do no = 1,ns_o
     if (pctglac_o(no) .gt. minglac) then 
        do m = 1,nglcec
           ! WJS (7-12-11): the original check was ((topo outside range) and (topo .ne.
           ! 0)). I think that the second condition is meant to essentially check whether
           ! pctglcmec_o(no,m) > 0. So I think we could just check whether ((topo outside
           ! range) and (pctglcmec_o(no,m) > 0). However, I am keeping a warning message
           ! if ((topo outside range) and (topo .ne. 0)) as well, because I don't want to
           ! get rid of error checks.
           if ( (topoglcmec_o(no,m) .lt. elevclass(m) .or. topoglcmec_o(no,m) .gt. elevclass(m+1)) &
                .and. (pctglcmec_o(no,m) .gt. 0 .or. topoglcmec_o(no,m) .ne. 0)) then
              write(6,*) 'Warning: mean elevation does not fall within elevation class '
              write(6,*) elevclass(m),elevclass(m+1),topoglcmec_o(no,m),pctglcmec_o(no,m),m,no
           endif
        end do
     end if
  end do

  if (errors) then
     call abort()
  end if

  ! Deallocate dynamic memory

  call domain_clean(tdomain_topo)
  call gridmap_clean(tgridmap)
  deallocate (topoice_g, pctglc_gic_g, pctglc_icesheet_g, pctglc_ground_g, pctglc_float_g,&
       pctglac_g)
  deallocate(topoglcmec_unnorm_o)

  write (6,*) 'Successfully made percent elevation class and mean elevation for glaciers'
  write (6,*)
  call shr_sys_flush(6)

end subroutine mkglcmec

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: mkglacier
!
! !INTERFACE:
subroutine mkglacier(ldomain, mapfname, datfname, ndiag, zero_out, glac_o, glac_uncorrected)
!
! !DESCRIPTION:
! make percent glacier
!
! !USES:
  use mkdomainMod , only : domain_type, domain_clean, domain_read
  use mkgridmapMod
  use mkvarpar	
  use mkvarctl    
  use mkncdio
!
! !ARGUMENTS:
  implicit none
  type(domain_type), intent(in) :: ldomain
  character(len=*)  , intent(in) :: mapfname  ! input mapping file name
  character(len=*)  , intent(in) :: datfname  ! input data file name
  integer           , intent(in) :: ndiag     ! unit number for diag out
  logical           , intent(in) :: zero_out  ! if should zero glacier out
  real(r8)          , intent(out):: glac_o(:) ! output grid: %glacier
  real(r8)          , intent(out):: glac_uncorrected(:) ! output grid: %glacier before any
                                                        ! corrections are done
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
  type(gridmap_type)   :: tgridmap
  type(domain_type)    :: tdomain            ! local domain
  real(r8), allocatable :: glac_i(:)          ! input grid: percent glac
  real(r8) :: sum_fldi                        ! global sum of dummy input fld
  real(r8) :: sum_fldo                        ! global sum of dummy output fld
  real(r8) :: gglac_i                         ! input  grid: global glac
  real(r8) :: garea_i                         ! input  grid: global area
  real(r8) :: gglac_o                         ! output grid: global glac
  real(r8) :: garea_o                         ! output grid: global area
  integer  :: ni,no,k,n,m,ns                  ! indices
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

  call domain_read(tdomain,datfname)
  ns = tdomain%ns
  allocate(glac_i(ns), stat=ier)
  if (ier/=0) call abort()

  write (6,*) 'Open glacier file: ', trim(datfname)
  call check_ret(nf_open(datfname, 0, ncid), subname)
  call check_ret(nf_inq_varid (ncid, 'PCT_GLACIER', varid), subname)
  call check_ret(nf_get_var_double (ncid, varid, glac_i), subname)
  call check_ret(nf_close(ncid), subname)

  ! Area-average percent cover on input grid to output grid 
  ! and correct according to land landmask
  ! Note that percent cover is in terms of total grid area.

  if ( zero_out )then

     do no = 1, ldomain%ns
        glac_o(no) = 0.
     enddo

  else

     call gridmap_mapread(tgridmap, mapfname )

     ! Error checks for domain and map consistencies
     call domain_checksame( tdomain, ldomain, tgridmap )
     
     ! Determine glac_o on output grid

     call gridmap_areaave(tgridmap, glac_i, glac_o)
     
     ! Save a copy of glac_o before any corrections are done. This is needed for
     ! normalization in mkglcmec
     glac_uncorrected(:) = glac_o(:)

     do no = 1, ldomain%ns
        if (glac_o(no) < 1.) glac_o(no) = 0.
     enddo
  end if

  ! Check for conservation

  do no = 1, ldomain%ns
     if ((glac_o(no)) > 100.000001_r8) then
        write (6,*) 'MKGLACIER error: glacier = ',glac_o(no), &
                ' greater than 100.000001 for column, row = ',no
        call shr_sys_flush(6)
        stop
     end if
  enddo

  ! Some error checking and writing of global values before and after the regrid

  if ( .not. zero_out )then

     ! Global sum of output field -- must multiply by fraction of
     ! output grid that is land as determined by input grid

     sum_fldi = 0.0_r8
     do ni = 1, tdomain%ns
        sum_fldi = sum_fldi + tgridmap%area_src(ni) * tgridmap%frac_src(ni)
     enddo

     sum_fldo = 0.
     do no = 1, ldomain%ns
        sum_fldo = sum_fldo + tgridmap%area_dst(no) * tgridmap%frac_dst(no)
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
     do ni = 1, tdomain%ns
        garea_i = garea_i + tgridmap%area_src(ni)*re**2
        gglac_i = gglac_i + glac_i(ni)*(tgridmap%area_src(ni)/100.)*&
                                        tgridmap%frac_src(ni)*re**2
     end do

     ! Output grid

     gglac_o = 0.
     garea_o = 0.
     do no = 1, ldomain%ns
        garea_o = garea_o + tgridmap%area_dst(no)*re**2
        gglac_o = gglac_o + glac_o(no)*(tgridmap%area_dst(no)/100.)*&
                                        tgridmap%frac_dst(no)*re**2
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

  end if

  ! Deallocate dynamic memory

  call domain_clean(tdomain) 
  if ( .not. zero_out )then
     call gridmap_clean(tgridmap)
     deallocate (glac_i)
  end if

  write (6,*) 'Successfully made %glacier'
  write (6,*)
  call shr_sys_flush(6)

end subroutine mkglacier

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_elevclass
!
! !INTERFACE:
integer function get_elevclass(topo, writewarn)
!
! !DESCRIPTION:
! Returns elevation class index (1..nglcec) given the topographic height.
! If topo is lower than the lowest elevation class, returns 0.
! If topo is higher than the highest elevation class, returns (nglcec+1).
! In either of the two latter cases, the function also writes a warning message, unless
! writewarn is present and false.
!
! !ARGUMENTS:
   implicit none
   real(r8), intent(in) :: topo  ! topographic height (m)
   logical, intent(in), optional :: writewarn  ! should warning messages be written? (default: true)
!
! !REVISION HISTORY:
! Author: Bill Sacks
!
! !LOCAL VARIABLES:
!EOP
   integer :: m
   logical :: my_writewarn
   character(len=32) :: subname = 'mkglcmec'
!-----------------------------------------------------------------------
   
   if (present(writewarn)) then
      my_writewarn = writewarn
   else
      my_writewarn = .true.
   end if

   if (topo < elevclass(1)) then
      if (my_writewarn) then
         write(6,*) 'WARNING in ', trim(subname)
         write(6,*) 'topo out of bounds'
         write(6,*) 'topo = ', topo
         write(6,*) 'elevclass(1) = ', elevclass(1)
      end if
      get_elevclass = 0
      return
   end if
   
   do m = 1, nglcec
      if (topo < elevclass(m+1)) then
         ! note that we already know that topo >= elevclass(m), otherwise we would have
         ! returned earlier
         get_elevclass = m
         return
      end if
   end do

   if (my_writewarn) then
      write(6,*) 'WARNING in ', trim(subname)
      write(6,*) 'topo out of bounds'
      write(6,*) 'topo = ', topo
      write(6,*) 'elevclass(nglcec+1) = ', elevclass(nglcec+1)
   end if
   get_elevclass = nglcec+1

end function get_elevclass
         
!-----------------------------------------------------------------------

end module mkglcmecMod
