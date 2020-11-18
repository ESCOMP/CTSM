module HillslopeHydrologyMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Read geomorphological parameters for hillslope columns 
  !
  ! !USES:
#include "shr_assert.h"
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use shr_log_mod    , only : errMsg => shr_log_errMsg
  use spmdMod        , only : masterproc
  use abortutils     , only : endrun
  use clm_varctl     , only : iulog
  use decompMod      , only : bounds_type
  use clm_varcon     , only : rpi

  ! !PUBLIC TYPES:
  implicit none
  private   
  save

  ! !PUBLIC MEMBER FUNCTIONS:
  public InitHillslope
  public HillslopeSoilThicknessProfile
  public HillslopeSetLowlandUplandPfts
  public HillslopeDominantPft
  public HillslopeDominantLowlandPft
  public HillslopePftFromFile

  ! PRIVATE 
  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  integer, private, parameter :: soil_profile_set_lowland_upland    = 0 
  integer, private, parameter :: soil_profile_linear                = 1


  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------

  subroutine InitHillslope(bounds,fsurdat)
    !
    ! !DESCRIPTION:
    ! Initialize hillslope geomorphology from input dataset
    !
    ! !USES:
    use LandunitType    , only : lun                
    use GridcellType    , only : grc                
    use ColumnType      , only : col                
    use clm_varctl      , only : nhillslope, nmax_col_per_hill
    use clm_varcon      , only : zmin_bedrock, zisoi
    use clm_varpar      , only : nlevsoi
    use spmdMod         , only : masterproc
    use fileutils       , only : getfil
    use clm_varcon      , only : spval, ispval, grlnd 
    use landunit_varcon , only : istsoil
    use ncdio_pio

    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds
    character(len=*) , intent(in) :: fsurdat    ! surface data file name
    integer,  pointer     :: ihillslope_in(:,:) ! read in - integer
    integer,  allocatable :: hill_ndx(:,:)      ! hillslope index
    integer,  allocatable :: col_ndx(:,:)       ! column index
    integer,  allocatable :: col_dndx(:,:)      ! downhill column index
    integer,  allocatable :: hill_pftndx(:,:)   ! hillslope pft index []
    real(r8), pointer     :: fhillslope_in(:,:) ! read in - float
    real(r8), allocatable :: pct_hillslope(:,:) ! percent of landunit occupied by hillslope
    real(r8), allocatable :: hill_slope(:,:)    ! hillslope slope  [m/m]
    real(r8), allocatable :: hill_aspect(:,:)   ! hillslope azimuth [radians]
    real(r8), allocatable :: hill_area(:,:)     ! hillslope area   [m2]
    real(r8), allocatable :: hill_length(:,:)   ! hillslope length [m]
    real(r8), allocatable :: hill_width(:,:)    ! hillslope width  [m]
    real(r8), allocatable :: hill_height(:,:)   ! hillslope height [m]
    real(r8), allocatable :: hill_bedrock(:,:)  ! hillslope bedrock depth [m]

    type(file_desc_t)     :: ncid                 ! netcdf id
    logical               :: readvar              ! check whether variable on file    
    character(len=256)    :: locfn                ! local filename
    integer               :: ierr                 ! error code
    integer               :: c, l, g, i, j, ci, nh       ! indices

    real(r8)              :: hillslope_area       ! total area of hillslope

    character(len=*), parameter :: subname = 'InitHillslope'

    !-----------------------------------------------------------------------

    ! Open surface dataset to read in data below 

    call getfil (fsurdat, locfn, 0)
    call ncd_pio_openfile (ncid, locfn, 0)
    
    allocate( &
         pct_hillslope(bounds%begl:bounds%endl,nhillslope),  &
         hill_ndx     (bounds%begl:bounds%endl,nmax_col_per_hill), &
         col_ndx      (bounds%begl:bounds%endl,nmax_col_per_hill), &
         col_dndx     (bounds%begl:bounds%endl,nmax_col_per_hill), &
         hill_slope   (bounds%begl:bounds%endl,nmax_col_per_hill), &
         hill_aspect  (bounds%begl:bounds%endl,nmax_col_per_hill), &
         hill_area    (bounds%begl:bounds%endl,nmax_col_per_hill), &
         hill_length  (bounds%begl:bounds%endl,nmax_col_per_hill), &
         hill_width   (bounds%begl:bounds%endl,nmax_col_per_hill), &
         hill_height  (bounds%begl:bounds%endl,nmax_col_per_hill), &
         stat=ierr)
    
    allocate(fhillslope_in(bounds%begg:bounds%endg,nhillslope))
    
    call ncd_io(ncid=ncid, varname='pct_hillslope', flag='read', data=fhillslope_in, dim1name=grlnd, readvar=readvar)
    if (.not. readvar) then
       if (masterproc) then
          call endrun( 'ERROR:: pct_hillslope not found on surface data set.'//errmsg(sourcefile, __LINE__) )
       end if
    end if
    do l = bounds%begl,bounds%endl
       g = lun%gridcell(l)
       pct_hillslope(l,:) = fhillslope_in(g,:)
    enddo
    deallocate(fhillslope_in)
    
    allocate(ihillslope_in(bounds%begg:bounds%endg,nmax_col_per_hill))
    
    call ncd_io(ncid=ncid, varname='hillslope_index', flag='read', data=ihillslope_in, dim1name=grlnd, readvar=readvar)
    if (.not. readvar) then
       if (masterproc) then
          call endrun( 'ERROR:: hillslope_index not found on surface data set.'//errmsg(sourcefile, __LINE__) )
       end if
    end if
    do l = bounds%begl,bounds%endl
       g = lun%gridcell(l)
       hill_ndx(l,:) = ihillslope_in(g,:)
    enddo
    
    call ncd_io(ncid=ncid, varname='column_index', flag='read', data=ihillslope_in, dim1name=grlnd, readvar=readvar)
    if (.not. readvar) then
       if (masterproc) then
          call endrun( 'ERROR:: column_index not found on surface data set.'//errmsg(sourcefile, __LINE__) )
       end if
    end if
    do l = bounds%begl,bounds%endl
       g = lun%gridcell(l)
       col_ndx(l,:) = ihillslope_in(g,:)
    enddo
    
    call ncd_io(ncid=ncid, varname='downhill_column_index', flag='read', data=ihillslope_in, dim1name=grlnd, readvar=readvar)
    if (.not. readvar) then
       if (masterproc) then
          call endrun( 'ERROR:: downhill_column_index not found on surface data set.'//errmsg(sourcefile, __LINE__) )
       end if
    end if
    do l = bounds%begl,bounds%endl
       g = lun%gridcell(l)
       col_dndx(l,:) = ihillslope_in(g,:)
    enddo
    deallocate(ihillslope_in)
    
    allocate(fhillslope_in(bounds%begg:bounds%endg,nmax_col_per_hill))
    call ncd_io(ncid=ncid, varname='h_slope', flag='read', data=fhillslope_in, dim1name=grlnd, readvar=readvar)
    if (.not. readvar) then
       if (masterproc) then
          call endrun( 'ERROR:: h_slope not found on surface data set.'//errmsg(sourcefile, __LINE__) )
       end if
    end if
    
    do l = bounds%begl,bounds%endl
       g = lun%gridcell(l)
       hill_slope(l,:) = fhillslope_in(g,:)
    enddo
    
    call ncd_io(ncid=ncid, varname='h_aspect', flag='read', data=fhillslope_in, dim1name=grlnd, readvar=readvar)
    if (.not. readvar) then
       if (masterproc) then
          call endrun( 'ERROR:: h_aspect not found on surface data set.'//errmsg(sourcefile, __LINE__) )
       end if
    end if
    
    do l = bounds%begl,bounds%endl
       g = lun%gridcell(l)
       hill_aspect(l,:) = fhillslope_in(g,:)
    enddo
    
    call ncd_io(ncid=ncid, varname='h_area', flag='read', data=fhillslope_in, dim1name=grlnd, readvar=readvar)
    if (.not. readvar) then
       if (masterproc) then
          call endrun( 'ERROR:: h_area not found on surface data set.'//errmsg(sourcefile, __LINE__) )
       end if
    end if
    do l = bounds%begl,bounds%endl
       g = lun%gridcell(l)
       hill_area(l,:) = fhillslope_in(g,:)
    enddo
    call ncd_io(ncid=ncid, varname='h_length', flag='read', data=fhillslope_in, dim1name=grlnd, readvar=readvar)
    if (.not. readvar) then
       if (masterproc) then
          call endrun( 'ERROR:: h_length not found on surface data set.'//errmsg(sourcefile, __LINE__) )
       end if
    end if
    
    do l = bounds%begl,bounds%endl
       g = lun%gridcell(l)
       hill_length(l,:) = fhillslope_in(g,:)
    enddo
    
    call ncd_io(ncid=ncid, varname='h_width', flag='read', data=fhillslope_in, dim1name=grlnd, readvar=readvar)
    if (.not. readvar) then
       if (masterproc) then
          call endrun( 'ERROR:: h_width not found on surface data set.'//errmsg(sourcefile, __LINE__) )
       end if
    end if
    do l = bounds%begl,bounds%endl
       g = lun%gridcell(l)
       hill_width(l,:) = fhillslope_in(g,:)
    enddo
    
    call ncd_io(ncid=ncid, varname='h_height', flag='read', data=fhillslope_in, dim1name=grlnd, readvar=readvar)
    if (.not. readvar) then
       if (masterproc) then
          call endrun( 'ERROR:: h_height not found on surface data set.'//errmsg(sourcefile, __LINE__) )
       end if
    end if
    do l = bounds%begl,bounds%endl
       g = lun%gridcell(l)
       hill_height(l,:) = fhillslope_in(g,:)
    enddo

    call ncd_io(ncid=ncid, varname='h_bedrock', flag='read', data=fhillslope_in, dim1name=grlnd, readvar=readvar)
    if (readvar) then
       allocate(hill_bedrock (bounds%begl:bounds%endl,nmax_col_per_hill), stat=ierr)
       if (masterproc) then
          write(iulog,*) 'h_bedrock found on surface data set'
       end if
    end if
    do l = bounds%begl,bounds%endl
       g = lun%gridcell(l)
       hill_bedrock(l,:) = fhillslope_in(g,:)
    enddo

    deallocate(fhillslope_in)
    
    allocate(ihillslope_in(bounds%begg:bounds%endg,nmax_col_per_hill))
    call ncd_io(ncid=ncid, varname='h_pftndx', flag='read', data=ihillslope_in, dim1name=grlnd, readvar=readvar)
    if (readvar) then
       allocate(hill_pftndx (bounds%begl:bounds%endl,nmax_col_per_hill), stat=ierr)
       if (masterproc) then
          write(iulog,*) 'h_pftndx found on surface data set'
       end if
       do l = bounds%begl,bounds%endl
          g = lun%gridcell(l)
          hill_pftndx(l,:) = ihillslope_in(g,:)
       enddo
    end if

    deallocate(ihillslope_in)

    !  Set hillslope hydrology column level variables
    !  This needs to match how columns set up in subgridMod
    do l = bounds%begl,bounds%endl
       g = lun%gridcell(l)
       if(lun%itype(l) == istsoil) then
          ! map external column index to internal column index
          do c = lun%coli(l), lun%colf(l)
             ! ci should span [1:nhillcolumns(l)]
             ci = c-lun%coli(l)+1

             if (col_dndx(l,ci) <= -999) then
                ! lowermost column of hillslope has no downstream neighbor
                col%cold(c) = ispval
             else
                ! relative separation should be the same
                col%cold(c) = c + (col_dndx(l,ci) - col_ndx(l,ci))
             endif

          enddo
          
          do c = lun%coli(l), lun%colf(l)
             ci = c-lun%coli(l)+1
             
             col%hillslope_ndx(c) = hill_ndx(l,ci)

             ! Find uphill neighbors (this may not actually be useful...)
             col%colu(c) = ispval
             do i = lun%coli(l), lun%colf(l)
                if(c == col%cold(i)) then
                   col%colu(c) = i
                endif
             enddo
             
             ! distance of lower edge of column from hillslope bottom
             col%hill_distance(c) = hill_length(l,ci)
             ! width of lower edge of column 
             col%hill_width(c) = hill_width(l,ci)
             ! mean elevation of column relative to gridcell mean elevation
             col%hill_elev(c) = hill_height(l,ci)
             ! mean along-hill slope of column
             col%hill_slope(c) = hill_slope(l,ci)
             ! area of column
             col%hill_area(c) = hill_area(l,ci)
             ! azimuth of column
             col%hill_aspect(c) = hill_aspect(l,ci)
             ! pft index of column
             if ( allocated(hill_bedrock) ) then
                do j = 1,nlevsoi
                   if(zisoi(j-1) > zmin_bedrock) then
                      if (zisoi(j-1) < hill_bedrock(l,ci) .and. zisoi(j) >= hill_bedrock(l,ci)) then
                         col%nbedrock(c) = j
                      end if
                   endif
                enddo
             endif
             if ( allocated(hill_pftndx) ) &
                  col%hill_pftndx(c) = hill_pftndx(l,ci)

          enddo

          ! Calculate total (representative) hillslope area on landunit
          ! weighted relative to one another via pct_hillslope
          hillslope_area = 0._r8
          do c = lun%coli(l), lun%colf(l)
             nh = col%hillslope_ndx(c)
             if (nh > 0) then
                hillslope_area = hillslope_area &
                     + col%hill_area(c)*(pct_hillslope(l,nh)*0.01_r8)
             endif
          enddo
          
          ! if missing hillslope information on surface dataset, fill data
          ! and recalculate hillslope_area
          if (hillslope_area == 0._r8) then
             do c = lun%coli(l), lun%colf(l)
                col%hill_area(c) = (grc%area(g)/real(lun%ncolumns(l),r8))*1.e6 ! km2 to m2
                col%hill_distance(c) = sqrt(col%hill_area(c))
                col%hill_width(c)    = sqrt(col%hill_area(c))
                col%hill_elev(c)     = col%topo_std(c)
                col%hill_slope(c)    = tan((rpi/180.)*col%topo_slope(c))
                col%hill_aspect(c)   = (rpi/2.) ! east (arbitrarily chosen)
                nh = col%hillslope_ndx(c)
                pct_hillslope(l,nh)  = 100/nhillslope
                hillslope_area = hillslope_area &
                     + col%hill_area(c)*(pct_hillslope(l,nh)*0.01_r8)
             enddo
          endif
          
          ! Recalculate column weights using input areas
          do c = lun%coli(l), lun%colf(l)
             nh = col%hillslope_ndx(c)
             if (nh > 0) then
                col%wtlunit(c) = col%hill_area(c) &
                     * (pct_hillslope(l,nh)*0.01_r8)/hillslope_area
             else
                col%wtlunit(c) = 0._r8
             endif
          enddo

       endif
    enddo
    
    deallocate(pct_hillslope,hill_ndx,col_ndx,col_dndx, &
         hill_slope,hill_area,hill_length, &
         hill_width,hill_height,hill_aspect)

    if ( allocated(hill_bedrock) ) then
       deallocate(hill_bedrock)
    else
       ! Modify hillslope soil thickness profile
       call HillslopeSoilThicknessProfile(bounds,&
            soil_profile_method=soil_profile_set_lowland_upland,&
            soil_depth_lowland_in=8.0_r8,soil_depth_upland_in=8.0_r8)

    endif
    if ( allocated(hill_pftndx) ) then
       deallocate(hill_pftndx)
       call HillslopePftFromFile()
    else
       ! Modify pft distributions
       ! this may require modifying subgridMod/natveg_patch_exists
       ! to ensure patch exists in every gridcell

       call HillslopeDominantPft()
    
       !upland_ivt  = 13 ! c3 non-arctic grass
       !lowland_ivt = 7  ! broadleaf deciduous tree 
       !call HillslopeSetLowlandUplandPfts(lowland_ivt=7,upland_ivt=13)
    endif

    call ncd_pio_closefile(ncid)
    
  end subroutine InitHillslope

  !-----------------------------------------------------------------------
  subroutine HillslopeSoilThicknessProfile(bounds,&
       soil_profile_method,soil_depth_lowland_in,soil_depth_upland_in)
    !
    ! !DESCRIPTION:
    ! Modify soil thickness across hillslope by changing
    ! col%nbedrock
    !
    ! !USES:
    use LandunitType    , only : lun                
    use GridcellType    , only : grc                
    use ColumnType      , only : col                
    use clm_varcon      , only : zmin_bedrock, zisoi
    use clm_varpar      , only : nlevsoi
    use spmdMod         , only : masterproc
    use fileutils       , only : getfil
    use clm_varcon      , only : spval, ispval, grlnd 
    use landunit_varcon , only : istsoil
    use ncdio_pio

    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds
    integer,  intent(in)  :: soil_profile_method
    real(r8), intent(in), optional  :: soil_depth_lowland_in
    real(r8), intent(in), optional  :: soil_depth_upland_in

    integer               :: c, l, g, i, j
    real(r8)              :: min_hill_dist, max_hill_dist
    real(r8)              :: m, b           ! linear soil thickness slope/intercept
    real(r8)              :: soil_depth_col
    real(r8)              :: soil_depth_lowland
    real(r8)              :: soil_depth_upland
    real(r8), parameter   :: soil_depth_lowland_default = 8.0
    real(r8), parameter   :: soil_depth_upland_default  = 8.0

    character(len=*), parameter :: subname = 'HillslopeSoilThicknessProfile'

    !-----------------------------------------------------------------------

    if(present(soil_depth_lowland_in)) then
       soil_depth_lowland = soil_depth_lowland_in
    else
       soil_depth_lowland = soil_depth_lowland_default
    endif
    
    if(present(soil_depth_upland_in)) then
       soil_depth_upland = soil_depth_upland_in
    else
       soil_depth_upland = soil_depth_upland_default
    endif

    do l = bounds%begl,bounds%endl
       if(lun%itype(l) == istsoil) then
          ! Specify lowland/upland soil thicknesses separately
          if(soil_profile_method == soil_profile_set_lowland_upland) then
             do c =  lun%coli(l), lun%colf(l)
                if(col%cold(c) /= ispval) then 
                   do j = 1,nlevsoi
                      if(zisoi(j-1) > zmin_bedrock) then
                         if (zisoi(j-1) < soil_depth_upland .and. zisoi(j) >= soil_depth_upland) then
                            col%nbedrock(c) = j
                         end if
                      end if
                   enddo
                else 
                   do j = 1,nlevsoi 
                      if(zisoi(j-1) > zmin_bedrock) then
                         if (zisoi(j-1) < soil_depth_lowland .and. zisoi(j) >= soil_depth_lowland) then
                            col%nbedrock(c) = j
                         end if
                      end if
                   enddo
                endif
             end do
          endif
          
          ! Linear soil thickness profile
          if(soil_profile_method == soil_profile_linear) then

             min_hill_dist = minval(col%hill_distance(lun%coli(l):lun%colf(l)))
             max_hill_dist = maxval(col%hill_distance(lun%coli(l):lun%colf(l)))
             m = (soil_depth_lowland - soil_depth_upland)/ &
                  (max_hill_dist - min_hill_dist)
             b = soil_depth_upland
             
             do c =  lun%coli(l), lun%colf(l)
                soil_depth_col = m*(max_hill_dist - col%hill_distance(c)) + b

                do j = 1,nlevsoi
                   if ((zisoi(j-1) <  soil_depth_col) .and. (zisoi(j) >= soil_depth_col)) then
                      col%nbedrock(c) = j
                   end if
                enddo
             enddo
          endif

       endif ! end of istsoil
    enddo    ! end of loop over landunits
       
  end subroutine HillslopeSoilThicknessProfile

  !------------------------------------------------------------------------
  subroutine HillslopeSetLowlandUplandPfts(lowland_ivt,upland_ivt)
    !
    ! !DESCRIPTION: 
    ! Reassign patch weights such that each column has a single pft.
    ! upland_ivt/lowland_ivt patches must be allocated
    ! in natveg_patch_exists (subgridMod) even if zero weight on fsurdat

    !
    ! !USES
    use LandunitType    , only : lun                
    use ColumnType      , only : col                
    use decompMod       , only : get_clump_bounds, get_proc_clumps
    use clm_varcon      , only : ispval
    use landunit_varcon , only : istsoil, istcrop
    use PatchType       , only : patch
    !
    ! !ARGUMENTS:
    !
    ! !LOCAL VARIABLES:
    integer :: n,nc,p,pu,pl,l,c    ! indices
    integer :: nclumps             ! number of clumps on this processor
    
    integer, intent(in) :: upland_ivt
    integer, intent(in) :: lowland_ivt
    real(r8) :: sum_wtcol, sum_wtlun, sum_wtgrc
    type(bounds_type) :: bounds_proc
    type(bounds_type) :: bounds_clump

    !------------------------------------------------------------------------

    nclumps = get_proc_clumps()

    !$OMP PARALLEL DO PRIVATE (nc, bounds_clump, l, nh, n, c)
    do nc = 1, nclumps

       call get_clump_bounds(nc, bounds_clump)

       do l = bounds_clump%begl, bounds_clump%endl

          if (lun%itype(l) == istsoil) then
             do c = lun%coli(l), lun%colf(l)

                sum_wtcol = sum(patch%wtcol(col%patchi(c):col%patchf(c)))
                sum_wtlun = sum(patch%wtlunit(col%patchi(c):col%patchf(c)))
                sum_wtgrc = sum(patch%wtgcell(col%patchi(c):col%patchf(c)))
                pl = ispval
                pu = ispval
                do p = col%patchi(c), col%patchf(c)
                   if(patch%itype(p) == lowland_ivt) pl = p
                   if(patch%itype(p) == upland_ivt)  pu = p
                enddo

                ! only reweight if pfts exist within column
                if (pl /= ispval .and. pu /= ispval) then
                   patch%wtcol(col%patchi(c):col%patchf(c)) = 0._r8
                   patch%wtlunit(col%patchi(c):col%patchf(c)) = 0._r8
                   patch%wtgcell(col%patchi(c):col%patchf(c)) = 0._r8

                ! hillbottom
                   if(col%cold(c) == ispval) then
                      patch%wtcol(pl) = sum_wtcol
                      patch%wtlunit(pl) = sum_wtlun
                      patch%wtgcell(pl) = sum_wtgrc
                   else
                      patch%wtcol(pu) = sum_wtcol
                      patch%wtlunit(pu) = sum_wtlun
                      patch%wtgcell(pu) = sum_wtgrc
                   endif
                endif
             enddo    ! end loop c
          endif
       enddo ! end loop l
    enddo    ! end loop nc
    !$OMP END PARALLEL DO

  end subroutine HillslopeSetLowlandUplandPfts

  !------------------------------------------------------------------------
  subroutine HillslopeDominantPft()
    !
    ! !DESCRIPTION: 
    ! Reassign patch weights such that each column has a single pft  
    ! determined by each column's most dominant pft on input dataset.  

    !
    ! !USES
    use LandunitType    , only : lun                
    use ColumnType      , only : col                
    use decompMod       , only : get_clump_bounds, get_proc_clumps
    use clm_varcon      , only : ispval
    use landunit_varcon , only : istsoil
    use PatchType       , only : patch
    !
    ! !ARGUMENTS:
    !
    ! !LOCAL VARIABLES:
    integer :: n,nc,p,pu,pl,l,c    ! indices
    integer :: nclumps             ! number of clumps on this processor
    integer :: pdom(1)
    real(r8) :: sum_wtcol, sum_wtlun, sum_wtgrc
    type(bounds_type) :: bounds_proc
    type(bounds_type) :: bounds_clump

    !------------------------------------------------------------------------

    nclumps = get_proc_clumps()

    !$OMP PARALLEL DO PRIVATE (nc, bounds_clump, l, nh, n, c)
    do nc = 1, nclumps

       call get_clump_bounds(nc, bounds_clump)

       do l = bounds_clump%begl, bounds_clump%endl

          if (lun%itype(l) == istsoil) then
             do c = lun%coli(l), lun%colf(l)

                pdom = maxloc(patch%wtcol(col%patchi(c):col%patchf(c)))
                pdom = pdom + (col%patchi(c) - 1)

                sum_wtcol = sum(patch%wtcol(col%patchi(c):col%patchf(c)))
                sum_wtlun = sum(patch%wtlunit(col%patchi(c):col%patchf(c)))
                sum_wtgrc = sum(patch%wtgcell(col%patchi(c):col%patchf(c)))

                patch%wtcol(col%patchi(c):col%patchf(c)) = 0._r8
                patch%wtlunit(col%patchi(c):col%patchf(c)) = 0._r8
                patch%wtgcell(col%patchi(c):col%patchf(c)) = 0._r8

                patch%wtcol(pdom(1)) = sum_wtcol
                patch%wtlunit(pdom(1)) = sum_wtlun
                patch%wtgcell(pdom(1)) = sum_wtgrc
             enddo    ! end loop c
          endif
       enddo ! end loop l
    enddo    ! end loop nc
    !$OMP END PARALLEL DO

  end subroutine HillslopeDominantPft

  !------------------------------------------------------------------------
  subroutine HillslopeDominantLowlandPft()
    !
    ! !DESCRIPTION: 
    ! Reassign patch weights such that each column has a single, 
    ! dominant pft.  Use largest weight for lowland, 2nd largest
    ! weight for uplands

    !
    ! !USES
    use LandunitType    , only : lun                
    use ColumnType      , only : col                
    use decompMod       , only : get_clump_bounds, get_proc_clumps
    use clm_varcon      , only : ispval
    use landunit_varcon , only : istsoil
    use PatchType       , only : patch
    !
    ! !ARGUMENTS:
    !
    ! !LOCAL VARIABLES:
    integer :: n,nc,p,pu,pl,l,c    ! indices
    integer :: nclumps             ! number of clumps on this processor
    integer :: pdom(1),psubdom(1)
    integer :: plow, phigh
    real(r8) :: sum_wtcol, sum_wtlun, sum_wtgrc
    real(r8),allocatable :: mask(:)
    type(bounds_type) :: bounds_proc
    type(bounds_type) :: bounds_clump

    !------------------------------------------------------------------------

    nclumps = get_proc_clumps()

    !$OMP PARALLEL DO PRIVATE (nc, bounds_clump, l, nh, n, c)
    do nc = 1, nclumps

       call get_clump_bounds(nc, bounds_clump)

       do l = bounds_clump%begl, bounds_clump%endl

          if (lun%itype(l) == istsoil) then
             do c = lun%coli(l), lun%colf(l)

                pdom = maxloc(patch%wtcol(col%patchi(c):col%patchf(c)))
                ! create mask to exclude pdom
                allocate(mask(col%npatches(c)))
                mask(:) = 1.
                mask(pdom(1)) = 0.

                pdom = pdom + (col%patchi(c) - 1)

                psubdom = maxloc(mask*patch%wtcol(col%patchi(c):col%patchf(c)))
                psubdom = psubdom + (col%patchi(c) - 1)
                deallocate(mask)

                sum_wtcol = sum(patch%wtcol(col%patchi(c):col%patchf(c)))
                sum_wtlun = sum(patch%wtlunit(col%patchi(c):col%patchf(c)))
                sum_wtgrc = sum(patch%wtgcell(col%patchi(c):col%patchf(c)))

                patch%wtcol(col%patchi(c):col%patchf(c)) = 0._r8
                patch%wtlunit(col%patchi(c):col%patchf(c)) = 0._r8
                patch%wtgcell(col%patchi(c):col%patchf(c)) = 0._r8

                ! assumes trees are 1-8, shrubs 9-11, and grasses 12-14
                ! and puts the lowest ivt on the lowland column
                if ((patch%itype(pdom(1)) > patch%itype(psubdom(1)) &
                     .and. patch%itype(psubdom(1)) > 0) .or. &
                     (patch%itype(pdom(1)) == 0)) then
                   plow = psubdom(1)
                   phigh = pdom(1)
                else
                   plow = pdom(1)
                   phigh = psubdom(1)
                endif

                ! Special cases (subjective)
                
                ! if NET/BDT assign BDT to lowland
                if ((patch%itype(pdom(1)) == 1) .and. (patch%itype(psubdom(1)) < 9)) then
                   plow = psubdom(1)        
                   phigh = pdom(1)
                endif
                ! if C3/C4 assign C4 to lowland
                if ((patch%itype(pdom(1)) == 14) .and. (patch%itype(psubdom(1)) == 13)) then
                   plow = pdom(1)        
                   phigh = psubdom(1)
                endif
                if ((patch%itype(pdom(1)) == 13) .and. (patch%itype(psubdom(1)) == 14)) then
                   plow = psubdom(1)        
                   phigh = pdom(1)
                endif
                
                if(col%cold(c) == ispval) then
                   ! lowland column
                   patch%wtcol(plow)   = sum_wtcol
                   patch%wtlunit(plow) = sum_wtlun
                   patch%wtgcell(plow) = sum_wtgrc
                else
                   ! upland columns
                   patch%wtcol(phigh)   = sum_wtcol
                   patch%wtlunit(phigh) = sum_wtlun
                   patch%wtgcell(phigh) = sum_wtgrc
                endif
             enddo    ! end loop c
          endif
       enddo ! end loop l
    enddo    ! end loop nc
    !$OMP END PARALLEL DO

  end subroutine HillslopeDominantLowlandPft

  !------------------------------------------------------------------------
  subroutine HillslopePftFromFile()
    !
    ! !DESCRIPTION: 
    ! Reassign patch weights using indices from surface data file
    !
    ! !USES
    use LandunitType    , only : lun                
    use ColumnType      , only : col                
    use decompMod       , only : get_clump_bounds, get_proc_clumps
    use clm_varcon      , only : ispval
    use landunit_varcon , only : istsoil
    use PatchType       , only : patch

    use GridcellType    , only : grc
    !
    ! !ARGUMENTS:
    !
    ! !LOCAL VARIABLES:
    integer :: n,nc,p,pc,l,c    ! indices
    integer :: nclumps             ! number of clumps on this processor
    real(r8) :: sum_wtcol, sum_wtlun, sum_wtgrc
    type(bounds_type) :: bounds_proc
    type(bounds_type) :: bounds_clump

    !------------------------------------------------------------------------

    nclumps = get_proc_clumps()

    !$OMP PARALLEL DO PRIVATE (nc, bounds_clump, l, nh, n, c)
    do nc = 1, nclumps

       call get_clump_bounds(nc, bounds_clump)

       do l = bounds_clump%begl, bounds_clump%endl

          if (lun%itype(l) == istsoil) then
             do c = lun%coli(l), lun%colf(l)
                ! this may require modifying
                ! subgridMod/natveg_patch_exists to ensure that
                ! a patch exists on each column
                
                ! find patch index of specified vegetation type
                pc = ispval
                do p = col%patchi(c), col%patchf(c)
                   if(patch%itype(p) == col%hill_pftndx(c)) pc = p
!scs                   if(patch%itype(p) == 1) pc = p
                enddo
                
                ! only reweight if pft exist within column
                if (pc /= ispval) then
                   sum_wtcol = sum(patch%wtcol(col%patchi(c):col%patchf(c)))
                   sum_wtlun = sum(patch%wtlunit(col%patchi(c):col%patchf(c)))
                   sum_wtgrc = sum(patch%wtgcell(col%patchi(c):col%patchf(c)))
                   
                   patch%wtcol(col%patchi(c):col%patchf(c)) = 0._r8
                   patch%wtlunit(col%patchi(c):col%patchf(c)) = 0._r8
                   patch%wtgcell(col%patchi(c):col%patchf(c)) = 0._r8

                   patch%wtcol(pc)   = sum_wtcol
                   patch%wtlunit(pc) = sum_wtlun
                   patch%wtgcell(pc) = sum_wtgrc

                else
                   write(iulog,*) 'no pft in column ',c, col%hill_pftndx(c)
                   write(iulog,*) 'pfts ',c,patch%itype(col%patchi(c):col%patchf(c))
                   write(iulog,*) 'weights ',c,patch%wtcol(col%patchi(c):col%patchf(c))
                   write(iulog,*) 'location ',c,grc%londeg(col%gridcell(c)),grc%latdeg(col%gridcell(c))
                   
                endif
             enddo    ! end loop c
          endif
       enddo ! end loop l
    enddo    ! end loop nc
    !$OMP END PARALLEL DO

  end subroutine HillslopePftFromFile

end module HillslopeHydrologyMod
