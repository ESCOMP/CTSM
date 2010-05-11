!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: mkglcmec
!
! !INTERFACE:
subroutine mkglcmec(lsmlon, lsmlat, fname1, fname2, fname3, ndiag, pctglac_o, pctglcmec_o, &
                    topoglcmec_o, thckglcmec_o, elevclass)
!
! !DESCRIPTION:
! make percent glacier on multiple elevation classes and mean elevation for each elevation class
!
! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use shr_sys_mod , only : shr_sys_flush
  use fileutils   , only : getfil
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
  real(r8), intent(out) :: elevclass(nglcec+1)                ! elevation classes
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
  character(len=256) locfn                  ! local dataset file name
  character(len=32) :: subname = 'mkglcmec'
!-----------------------------------------------------------------------

  write (6,*) 'Attempting to make percent elevation class and mean elevation for glaciers .....'
  call shr_sys_flush(6)

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
  else
     write(6,*) "ERROR:: nglcec must be 1, 3, 5, or 10 to work with CLM: "
     call abort()
  end if

  ! -----------------------------------------------------------------
  ! Read input file
  ! -----------------------------------------------------------------

  ! Obtain input grid info, read local fields

  call getfil (fname1, locfn, 0)

  call read_domain(tdomain1,locfn)
  call domain_setptrs(tdomain1,ni=nlont_i,nj=nlatt_i)

  call check_ret(nf_open(locfn, 0, ncid), subname)

  allocate(topoice_i(nlont_i, nlatt_i), topobedrock_i(nlont_i, nlatt_i), stat=ier)
  if (ier/=0) call abort()

  call check_ret(nf_inq_varid (ncid, 'TOPO_ICE', varid), subname)
  call check_ret(nf_get_var_double (ncid, varid, topoice_i), subname)
  call check_ret(nf_inq_varid (ncid, 'TOPO_BEDROCK', varid), subname)
  call check_ret(nf_get_var_double (ncid, varid, topobedrock_i), subname)

  call check_ret(nf_close(ncid), subname)

  call getfil (fname2, locfn, 0)

  call read_domain(tdomain2,locfn)
  call domain_setptrs(tdomain2,ni=nlont_i,nj=nlatt_i)

  call check_ret(nf_open(locfn, 0, ncid), subname)

  allocate(landmask_i(nlont_i, nlatt_i), stat=ier)                         
  if (ier/=0) call abort()

  call check_ret(nf_inq_varid (ncid, 'LANDMASK', varid), subname)    
  call check_ret(nf_get_var_double (ncid, varid, landmask_i), subname)   

  call check_ret(nf_close(ncid), subname)

  ! Get pct_glacier data


  call getfil (fname3, locfn, 0)

  call read_domain(tdomain3,locfn)
  call domain_setptrs(tdomain3,ni=nlong_i,nj=nlatg_i)

  call check_ret(nf_open(locfn, 0, ncid), subname)

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
     
     if (glac_i(io,jo) > 1.e-6) then 

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

  pctglcmec_o = 0.
  topoglcmec_o = 0.
  thckglcmec_o = 0.

  call areaini(tdomain3,ldomain,tgridmap,fracin=maskg_i,fracout=mask_o)

  call gridmap_setptrs(tgridmap,mx_ovr=mxovr,n_ovr=novr,i_ovr=iovr,j_ovr=jovr,w_ovr=wovr)

  do jo = 1, ldomain%nj
  do io = 1, ldomain%numlon(jo)

     pcttot = 0.
     pcttotec = 0.
     if (pctglac_o(io,jo) .gt. 1.e-6) then 

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

