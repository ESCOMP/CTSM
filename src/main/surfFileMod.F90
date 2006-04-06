#include <misc.h>
#include <preproc.h>
#if ( defined SCAM )
#include <max.h>
#endif

module surfFileMod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: surfFileMod
!
! !DESCRIPTION:
! Contains methods for reading in surface data file and determining
! two-dimensional subgrid weights as well as writing out new surface
! dataset. When reading in the surface dataset, determines array
! which sets the PFT for each of the [maxpatch] patches and
! array which sets the relative abundance of the PFT.
! Also fills in the PFTs for vegetated portion of each grid cell.
! Fractional areas for these points pertain to "vegetated"
! area not to total grid area. Need to adjust them for fraction of grid
! that is vegetated. Also fills in urban, lake, wetland, and glacier patches.
!
! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use shr_sys_mod , only : shr_sys_flush
  use abortutils  , only : endrun
#if ( defined SCAM )
  use scamMod     , only : initlonidx,initlatidx
#endif
  use clm_varpar  , only : lsmlon, lsmlat, nlevsoi, numpft, &
	                   maxpatch_pft, maxpatch_cft, maxpatch, &
	                   npatch_urban, npatch_lake, npatch_wet, npatch_glacier
  use pftvarcon   , only : crop, noveg
  use ncdio       
  use spmdMod                         
#if ( defined SCAM )
  use getnetcdfdata
#endif
!
! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: surfrd  ! Read surface dataset and determine subgrid weights
  public :: surfrd_get_grid  ! Read surface dataset into domain
  public :: surfrd_get_frac  ! Read surface dataset into domain
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!
! !PRIVATE MEMBER FUNCTIONS:
  private :: surfrd_wtxy_special
  private :: surfrd_wtxy_veg_rank
  private :: surfrd_wtxy_veg_all
  private :: surfrd_wtxy_veg_dgvm
  private :: mkrank
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: surfrd
!
! !INTERFACE:
  subroutine surfrd(veg, wtxy, lfsurdat)
!
! !DESCRIPTION:
! Read the surface dataset and create subgrid weights.
! The model's surface dataset recognizes 5 basic land cover types within
! a grid cell: lake, wetland, urban, glacier, and vegetated. The vegetated
! portion of the grid cell is comprised of up to [maxpatch_pft] PFTs. These
! subgrid patches are read in explicitly for each grid cell. This is in
! contrast to LSMv1, where the PFTs were built implicitly from biome types.
!    o real edges of grid
!    o integer  number of longitudes per latitude
!    o real latitude  of grid cell (degrees)
!    o real longitude of grid cell (degrees)
!    o integer surface type: 0 = ocean or 1 = land
!    o integer soil color (1 to 9) for use with soil albedos
!    o real soil texture, %sand, for thermal and hydraulic properties
!    o real soil texture, %clay, for thermal and hydraulic properties
!    o real % of cell covered by lake    for use as subgrid patch
!    o real % of cell covered by wetland for use as subgrid patch
!    o real % of cell that is urban      for use as subgrid patch
!    o real % of cell that is glacier    for use as subgrid patch
!    o integer PFTs
!    o real % abundance PFTs (as a percent of vegetated area)
!
! !USES:
    use clm_varsur  , only : all_pfts_on_srfdat, pctspec
    use clm_varctl  , only : allocate_all_vegpfts
    use fileutils   , only : getfil
!
! !ARGUMENTS:
    implicit none
    include 'netcdf.inc'
    integer , intent(out) :: veg(lsmlon,lsmlat,maxpatch)   ! PFT
    real(r8), intent(out) :: wtxy(lsmlon,lsmlat,maxpatch)  ! subgrid weights
    character(len=*), intent(in) :: lfsurdat               ! surf filename
!
! !CALLED FROM:
! subroutine initialize in module initializeMod
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein, Sam Levis and Gordon Bonan
!
!EOP
!
! !LOCAL VARIABLES:
    character(len=256) :: locfn                          ! local file name
    integer  :: ncid,dimid,varid                         ! netCDF id's
    logical  :: found                                    ! temporary for error check
    integer  :: iindx, jindx                             ! temporary for error check
    integer  :: ier                                      ! error status
    character(len=32) :: subname = 'surfrd'              ! subroutine name
!-----------------------------------------------------------------------

    veg(:,:,:)   = noveg
    wtxy(:,:,:)  = 0._r8
    pctspec(:,:) = 0._r8

    if (masterproc) then
       write (6,*) 'Attempting to read surface boundary data .....'
       if (lfsurdat == ' ') then
          write(6,*)'lfsurdat must be specified'; call endrun()
       endif
    endif

    ! Read surface data

    if (masterproc) then
       call getfil( lfsurdat, locfn, 0 )
       call check_ret( nf_open(locfn, 0, ncid), subname )
    end if

    ! If PFT variable is on surface dataset then not all pfts are
    ! on the surface dataset and old surface dataset form is assumed

    if (masterproc) then
       all_pfts_on_srfdat = .true.
       ier = nf_inq_varid (ncid, 'PFT', varid)
       if (ier == NF_NOERR) all_pfts_on_srfdat = .false.
    end if
#if (defined SPMD)
    call mpi_bcast (all_pfts_on_srfdat, 1, MPI_LOGICAL, 0, mpicom, ier)
#endif
     
    ! Obtain surface dataset special landunit info

    call surfrd_wtxy_special(ncid, pctspec, veg, wtxy)

    ! Obtain surface dataset vegetated landunit info

#if (! defined DGVM)     
    if (allocate_all_vegpfts) then
       call surfrd_wtxy_veg_all(ncid, pctspec, veg, wtxy)
    else
       call surfrd_wtxy_veg_rank(ncid, pctspec, veg, wtxy)
    end if
#else
    call surfrd_wtxy_veg_dgvm(pctspec, veg, wtxy)
#endif

    if ( masterproc )then
       call check_ret(nf_close(ncid), subname)
       write (6,*) 'Successfully read surface boundary data'
       write (6,*)
    end if

  end subroutine surfrd

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: surfrd_get_grid
!
! !INTERFACE:
  subroutine surfrd_get_grid(domain,filename)
!
! !DESCRIPTION:
! Read the surface dataset grid related information:
! o real edges of grid
! o integer  number of longitudes per latitude
! o real latitude  of grid cell (degrees)
! o real longitude of grid cell (degrees)
! For offline mode only the the grid does not have to be global.
! If grid is read in from dataset, grid is assumed to be global but
! does not have to be regular.
! If grid is generated by model, grid does not have to be global but 
! must then define the north, east, south, and west edges:
!
! o edges(1)    = northern edge of grid (degrees): >  -90 and <= 90
! o edges(2)    = eastern edge of grid (degrees) : see following notes
! o edges(3)    = southern edge of grid (degrees): >= -90 and <  90
! o edges(4)    = western edge of grid (degrees) : see following notes
!
!   For partial grids, northern and southern edges are any latitude
!   between 90 (North Pole) and -90 (South Pole). Western and eastern
!   edges are any longitude between -180 and 180, with longitudes
!   west of Greenwich negative. That is, western edge >= -180 and < 180;
!   eastern edge > western edge and <= 180.
!
!   For global grids, northern and southern edges are 90 (North Pole)
!   and -90 (South Pole). The western and eastern edges depend on
!   whether the grid starts at Dateline or Greenwich. Regardless,
!   these edges must span 360 degrees. Examples:
!
!                              West edge    East edge
!                            --------------------------------------------------
!  (1) Dateline            :        -180 to 180       (negative W of Greenwich)
!  (2) Greenwich (centered):    0 - dx/2 to 360 - dx/2
!
!    Grid 1 is the grid for offline mode
!    Grid 2 is the grid for cam and csm mode since the NCAR CAM
!    starts at Greenwich, centered on Greenwich
!
! !USES:
    use clm_varcon, only : spval
    use domainMod , only : domain_type,domain_init
    use areaMod   , only : celledge, cellarea                      
    use fileutils , only : getfil
!
! !ARGUMENTS:
    implicit none
    include 'netcdf.inc'
    type(domain_type),intent(inout) :: domain   ! domain to init
    character(len=*) ,intent(in)    :: filename ! grid filename
!
! !CALLED FROM:
! subroutine surfrd in this module
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: i,j                 ! indices
    integer :: ni,nj               ! size of grid on file
    integer :: ncid,dimid,varid    ! netCDF id's
    integer :: ier                 ! error status
    integer,allocatable :: numlon(:)  ! local numlon from file
    character(len=256)  :: locfn   ! local file name
#if ( defined SCAM )
    integer :: ret, time_index
#endif
    logical :: AREAset             ! true if area read from grid file
    logical :: NSEWset             ! true if lat/lon NSEW read from grid file
    logical :: EDGEset             ! true if EDGE NSEW read from grid file
    character(len=32) :: subname = 'surfrd_get_grid'     ! subroutine name
!-----------------------------------------------------------------------

    AREAset = .false.
    NSEWset = .false.
    EDGEset = .false.

    if (masterproc) then

       if (filename == ' ') then
          write(6,*) trim(subname),' ERROR: filename must be specified '
          call endrun()
       endif

       call getfil( filename, locfn, 0 )
       call check_ret( nf_open(locfn, 0, ncid), subname )

#if ( ! defined SCAM )
       call check_ret(nf_inq_dimid (ncid, 'lsmlon', dimid), subname)
       call check_ret(nf_inq_dimlen(ncid, dimid, ni), subname)
       call check_ret(nf_inq_dimid (ncid, 'lsmlat', dimid), subname)
       call check_ret(nf_inq_dimlen(ncid, dimid, nj), subname)
#else
       ni = lsmlon
       nj = lsmlat
#endif

    endif

#if (defined SPMD)
    call mpi_bcast (ni, 1, MPI_INTEGER, 0, mpicom, ier)
    call mpi_bcast (nj, 1, MPI_INTEGER, 0, mpicom, ier)
#endif

    call domain_init(domain,ni,nj)

    if (masterproc) then

       domain%edges(:) = spval
       ier = nf_inq_varid (ncid, 'EDGEN', varid)
       if (ier == NF_NOERR) then
          EDGEset = .true.
          call check_ret(nf_inq_varid(ncid, 'EDGEN', varid), subname)
          call check_ret(nf_get_var_double(ncid, varid, domain%edges(1)), subname)
          call check_ret(nf_inq_varid(ncid, 'EDGEE', varid), subname)
          call check_ret(nf_get_var_double(ncid, varid, domain%edges(2)), subname)
          call check_ret(nf_inq_varid(ncid, 'EDGES', varid), subname)
          call check_ret(nf_get_var_double(ncid, varid, domain%edges(3)), subname)
          call check_ret(nf_inq_varid (ncid, 'EDGEW', varid), subname)
          call check_ret(nf_get_var_double(ncid, varid, domain%edges(4)), subname)
          if (maxval(domain%edges) > 1.0e35) EDGEset = .false. !read garbage
       endif

       ier = nf_inq_varid (ncid, 'LATN', varid)
       if (ier == NF_NOERR) then
          NSEWset = .true.
          call check_ret(nf_inq_varid(ncid, 'LATN', varid), subname)
          call check_ret(nf_get_var_double(ncid, varid, domain%latn), subname)
          call check_ret(nf_inq_varid(ncid, 'LONE', varid), subname)
          call check_ret(nf_get_var_double(ncid, varid, domain%lone), subname)
          call check_ret(nf_inq_varid(ncid, 'LATS', varid), subname)
          call check_ret(nf_get_var_double(ncid, varid, domain%lats), subname)
          call check_ret(nf_inq_varid(ncid, 'LONW', varid), subname)
          call check_ret(nf_get_var_double(ncid, varid, domain%lonw), subname)
       endif

       ier = nf_inq_varid (ncid, 'AREA', varid)
       if (ier == NF_NOERR) then
          AREAset = .true.
          call check_ret(nf_inq_varid(ncid, 'AREA', varid), subname)
          call check_ret(nf_get_var_double(ncid, varid, domain%area), subname)
       endif

#if ( defined SCAM )
       time_index=1
       call getncdata (ncid, initLatIdx, initLonIdx, time_index,'LONGXY'  , domain%lonc, RET)
       call getncdata (ncid, initLatIdx, initLonIdx, time_index,'LATIXY'  , domain%latc, RET)
       domain%mask = 1
       call getncdata (ncid, initLatIdx, initLonIdx, time_index,'LANDMASK', domain%mask, RET)
       domain%pftm = domain%mask
       where (domain%mask <= 0)
          domain%pftm = -1
       endwhere
       call getncdata (ncid, initLatIdx, initLonIdx, time_index,'LANDFRAC', domain%frac, RET)
#else

       ! if NUMLON is on file, make sure it's not a reduced grid.
       ier = nf_inq_varid (ncid, 'NUMLON', varid)
       if (ier == NF_NOERR) then
          allocate(numlon(domain%nj))
          call check_ret(nf_get_var_int(ncid, varid, numlon), subname)
          if (minval(numlon) /= maxval(numlon)) then
             write(6,*) 'ERROR surfFileMod: NUMLON no longer supported ', &
                         minval(numlon),maxval(numlon)
             call endrun
          endif
          deallocate(numlon)
       endif

       call check_ret(nf_inq_varid(ncid, 'LONGXY' , varid), subname)
       call check_ret(nf_get_var_double(ncid, varid, domain%lonc), subname)

       call check_ret(nf_inq_varid(ncid, 'LATIXY', varid), subname)
       call check_ret(nf_get_var_double(ncid, varid, domain%latc), subname)

! set mask to 1 everywhere by default, override if LANDMASK exists
! if landmask exists, use it to set pftm (for older datasets)
! pftm should be overwritten below for newer datasets
       domain%mask = 1
       ier = nf_inq_varid(ncid, 'LANDMASK', varid)
       if (ier == NF_NOERR) then
          call check_ret(nf_get_var_int(ncid, varid, domain%mask), subname)
          domain%pftm = domain%mask
          where (domain%mask <= 0)
             domain%pftm = -1
          endwhere
       endif

       ier = nf_inq_varid (ncid, 'PFTDATA_MASK', varid)
       if (ier == NF_NOERR) then
          call check_ret(nf_get_var_int(ncid, varid, domain%pftm), subname)
       endif

#endif

!tcx fix, this or a test/abort should be added so overlaps can be computed
!tcx fix, this is demonstrated not bfb in cam bl311 test.
!tcx fix, see also lat_o_local in areaMod.F90
#if (1 == 0)
       ! Check lat limited to -90,90
       if (minval(domain%latc) < -90.0_r8 .or. &
           maxval(domain%latc) >  90.0_r8) then
           write(6,*) trim(subname),' Limiting lat/lon to [-90/90] from ', &
              minval(domain%latc),maxval(domain%latc)
           do j = 1,domain%nj
           do i = 1,domain%ni
             domain%latc(i,j) = min(domain%latc(i,j),  90.0_r8)
             domain%latc(i,j) = max(domain%latc(i,j), -90.0_r8)
           enddo
           enddo
       endif
#endif

       ! -------------------------------------------------------------------
       ! Define edges and area of grid cells
       ! -------------------------------------------------------------------

#if (defined COUP_CAM) || (defined COUP_CSM)
       if (.not.NSEWset) call celledge (domain)
       if (.not.AREAset) call cellarea (domain)
#endif

#if (defined OFFLINE)
       if (.not.NSEWset) then    
          if (.not.EDGEset) then           ! global grid without use of edges
            call celledge (domain)
          else                             ! regional regular grid 
            call celledge (domain, &
                           domain%edges(1), domain%edges(2), &
                           domain%edges(3), domain%edges(4))
          end if	
       endif
       if (.not.AREAset) then
          if (.not.EDGEset) then           ! global grid without use of edges
             call cellarea (domain)
          else                             ! regional regular grid 
             call cellarea (domain, &
                            domain%edges(1), domain%edges(2), &
                            domain%edges(3), domain%edges(4))
          endif
       end if	
#endif

       call check_ret(nf_close(ncid), subname)

    end if   ! end of if-masterproc block

#if (defined SPMD)
    call mpi_bcast (domain%latn , size(domain%latn) , MPI_REAL8  , 0, mpicom, ier)
    call mpi_bcast (domain%lats , size(domain%lats) , MPI_REAL8  , 0, mpicom, ier)
    call mpi_bcast (domain%lonw , size(domain%lonw) , MPI_REAL8  , 0, mpicom, ier)
    call mpi_bcast (domain%lone , size(domain%lone) , MPI_REAL8  , 0, mpicom, ier)
    call mpi_bcast (domain%area , size(domain%area) , MPI_REAL8  , 0, mpicom, ier)
    call mpi_bcast (domain%latc , size(domain%latc) , MPI_REAL8  , 0, mpicom, ier)
    call mpi_bcast (domain%lonc , size(domain%lonc) , MPI_REAL8  , 0, mpicom, ier)
    call mpi_bcast (domain%pftm , size(domain%pftm) , MPI_INTEGER, 0, mpicom, ier)
    call mpi_bcast (domain%edges, size(domain%edges), MPI_REAL8  , 0, mpicom, ier)
#endif

  end subroutine surfrd_get_grid

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: surfrd_get_frac
!
! !INTERFACE:
  subroutine surfrd_get_frac(domain,filename)
!
! !DESCRIPTION:
! Read the landfrac dataset grid related information:
! Assume domain has already been initialized and read
!
! !USES:
    use domainMod , only : domain_type
    use fileutils , only : getfil
!
! !ARGUMENTS:
    implicit none
    include 'netcdf.inc'
    type(domain_type),intent(inout) :: domain   ! domain to init
    character(len=*) ,intent(in)    :: filename ! grid filename
!
! !CALLED FROM:
! subroutine initialize
!
! !REVISION HISTORY:
! Created by T Craig
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: i,j                 ! indices
    integer :: ni,nj               ! size of grid on file
    integer :: ncid,dimid,varid    ! netCDF id's
    integer :: ier                 ! error status
    real(r8):: eps = 1.0e-12_r8    ! lat/lon error tolerance
    integer,allocatable :: numlon(:)  ! local numlon from file
    real(r8),allocatable:: lonc(:,:),latc(:,:)  ! local lat/lon
    character(len=256)  :: locfn   ! local file name
#if ( defined SCAM )
    integer :: ret, time_index
#endif
    character(len=32) :: subname = 'surfrd_get_frac'     ! subroutine name
!-----------------------------------------------------------------------

    if (masterproc) then

       if (filename == ' ') then
          write(6,*) trim(subname),' ERROR: filename must be specified '
          call endrun()
       endif

       call getfil( filename, locfn, 0 )
       call check_ret( nf_open(locfn, 0, ncid), subname )

#if ( ! defined SCAM )
       call check_ret(nf_inq_dimid (ncid, 'lsmlon', dimid), subname)
       call check_ret(nf_inq_dimlen(ncid, dimid, ni), subname)
       call check_ret(nf_inq_dimid (ncid, 'lsmlat', dimid), subname)
       call check_ret(nf_inq_dimlen(ncid, dimid, nj), subname)
#else
       ni = lsmlon
       nj = lsmlat
#endif

       if (domain%ni /= ni .or. domain%nj /= nj) then
          write(6,*) trim(subname),' ERROR: landfrac file mismatch ni,nj'
          call endrun()
       endif

       allocate(latc(ni,nj),lonc(ni,nj))

#if ( defined SCAM )
       time_index=1
       call getncdata (ncid, initLatIdx, initLonIdx, time_index,'LONGXY'  , lonc, RET)
       call getncdata (ncid, initLatIdx, initLonIdx, time_index,'LATIXY'  , latc, RET)
       do j=1,nj
       do i=1,ni
          if (abs(latc(i,j)-domain%latc(i,j)) > eps .or. &
              abs(lonc(i,j)-domain%lonc(i,j)) > eps) then
             write(6,*) trim(subname),' ERROR: landfrac file mismatch lat,lon'
             call endrun()
          endif
       enddo
       enddo
       call getncdata (ncid, initLatIdx, initLonIdx, time_index,'LANDMASK', domain%mask, RET)
       call getncdata (ncid, initLatIdx, initLonIdx, time_index,'LANDFRAC', domain%frac, RET)
#else

       ! if NUMLON is on file, make sure it's not a reduced grid.
       ier = nf_inq_varid (ncid, 'NUMLON', varid)
       if (ier == NF_NOERR) then
          allocate(numlon(domain%nj))
          call check_ret(nf_get_var_int(ncid, varid, numlon), subname)
          if (minval(numlon) /= maxval(numlon)) then
             write(6,*) 'ERROR surfFileMod: NUMLON no longer supported ', &
                         minval(numlon),maxval(numlon)
             call endrun
          endif
          deallocate(numlon)
       endif

       call check_ret(nf_inq_varid(ncid, 'LONGXY' , varid), subname)
       call check_ret(nf_get_var_double(ncid, varid, lonc), subname)

       call check_ret(nf_inq_varid(ncid, 'LATIXY', varid), subname)
       call check_ret(nf_get_var_double(ncid, varid, latc), subname)

       do j=1,nj
       do i=1,ni
          if (abs(latc(i,j)-domain%latc(i,j)) > eps .or. &
              abs(lonc(i,j)-domain%lonc(i,j)) > eps) then
             write(6,*) trim(subname),' ERROR: landfrac file mismatch lat,lon'
             call endrun()
          endif
       enddo
       enddo

       call check_ret(nf_inq_varid(ncid, 'LANDMASK', varid), subname)
       call check_ret(nf_get_var_int(ncid, varid, domain%mask), subname)

       call check_ret(nf_inq_varid(ncid, 'LANDFRAC', varid), subname)
       call check_ret(nf_get_var_double(ncid, varid, domain%frac), subname)

#endif

       deallocate(latc,lonc)

       call check_ret(nf_close(ncid), subname)

    end if   ! end of if-masterproc block

#if (defined SPMD)
    call mpi_bcast (domain%mask , size(domain%mask) , MPI_INTEGER, 0, mpicom, ier)
    call mpi_bcast (domain%frac , size(domain%frac) , MPI_REAL8  , 0, mpicom, ier)
#endif

  end subroutine surfrd_get_frac

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: surfrd_wtxy_special 
!
! !INTERFACE:
  subroutine surfrd_wtxy_special(ncid, pctspec, veg, wtxy)
!
! !DESCRIPTION:
! Determine weight with respect to gridcell of all special "pfts" as well
! as soil color and percent sand and clay
!
! !USES:
!
! !ARGUMENTS:
    implicit none
    include 'netcdf.inc'
    integer , intent(in)    :: ncid                         ! netcdf file id 
    real(r8), intent(inout) :: pctspec(lsmlon, lsmlat)      ! percent wrt gridcell of special landuntis
    integer , intent(inout) :: veg(lsmlon,lsmlat,maxpatch)  ! PFT
    real(r8), intent(inout) :: wtxy(lsmlon,lsmlat,maxpatch) ! subgrid weights
!
! !CALLED FROM:
! subroutine surfrd in this module
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein, Sam Levis and Gordon Bonan
!
!EOP
!
! !LOCAL VARIABLES:
    integer  :: i,j                        ! indices
    integer  :: dimid,varid                ! netCDF id's
#if ( defined SCAM )
    integer  :: ret, time_index
    real(r8) :: nlevsoidata(nlevsoi)
#endif
    logical  :: found                      ! temporary for error check
    integer  :: iindx, jindx               ! temporary for error check
    integer  :: ier                        ! error status
    real(r8) :: pctgla(lsmlon,lsmlat)      ! percent of grid cell that is glacier
    real(r8) :: pctlak(lsmlon,lsmlat)      ! percent of grid cell that is lake
    real(r8) :: pctwet(lsmlon,lsmlat)      ! percent of grid cell that is wetland
    real(r8) :: pcturb(lsmlon,lsmlat)      ! percent of grid cell that is urbanized
    character(len=32) :: subname = 'surfrd_wtxy_special'  ! subroutine name
!!-----------------------------------------------------------------------

    if (masterproc) then

       call check_dim(ncid, 'nlevsoi', nlevsoi)

       ! Obtain non-grid surface properties of surface dataset other than percent pft

#if ( defined SCAM )
       time_index=1
       call getncdata (ncid, initLatIdx, initLonIdx, time_index, 'PCT_WETLAND', pctwet     , ret)
       call getncdata (ncid, initLatIdx, initLonIdx, time_index, 'PCT_LAKE'   , pctlak     , ret)
       call getncdata (ncid, initLatIdx, initLonIdx, time_index, 'PCT_GLACIER', pctgla     , ret)
       call getncdata (ncid, initLatIdx, initLonIdx, time_index, 'PCT_URBAN'  , pcturb     , ret)
#else
       call check_ret(nf_inq_varid(ncid, 'PCT_WETLAND', varid), subname)
       call check_ret(nf_get_var_double(ncid, varid, pctwet), subname)

       call check_ret(nf_inq_varid(ncid, 'PCT_LAKE', varid), subname)
       call check_ret(nf_get_var_double(ncid, varid, pctlak), subname)

       call check_ret(nf_inq_varid(ncid, 'PCT_GLACIER', varid), subname)
       call check_ret(nf_get_var_double(ncid, varid, pctgla), subname)

       call check_ret(nf_inq_varid(ncid, 'PCT_URBAN', varid), subname)
       call check_ret(nf_get_var_double(ncid, varid, pcturb), subname)
#endif 

       ! Error check: glacier, lake, wetland, urban sum must be less than 100

       found = .false.
       do j = 1,lsmlat
          do i = 1,lsmlon
             pctspec(i,j) = pctlak(i,j)+pctwet(i,j)+pcturb(i,j)+pctgla(i,j)
             if (pctspec(i,j) > 100._r8+1.e-04_r8) then
                found = .true.
                iindx = i
                jindx = j
                exit
             end if
          end do
          if (found) exit
       end do
       if ( found ) then
          write(6,*)'surfrd error: PFT cover>100 for i,j=',iindx,jindx; call endrun()
       end if

       ! Error check that urban parameterization is not yet finished

       found = .false.
       do j = 1,lsmlat
          do i = 1,lsmlon
             if (pcturb(i,j) /= 0._r8) then
                found = .true.
                iindx = i
                jindx = j
                exit
             end if
          end do
          if (found) exit
       end do
       if ( found ) then
          write (6,*)'surfrd error: urban parameterization not implemented at i,j= ',iindx,jindx,pcturb(iindx,jindx)
          call endrun()
       end if

    end if   ! end of if-masterproc

#if (defined SPMD)
    call mpi_bcast (pctwet  , size(pctwet)  , MPI_REAL8  , 0, mpicom, ier)
    call mpi_bcast (pctlak  , size(pctlak)  , MPI_REAL8  , 0, mpicom, ier)
    call mpi_bcast (pctgla  , size(pctgla)  , MPI_REAL8  , 0, mpicom, ier)
    call mpi_bcast (pcturb  , size(pcturb)  , MPI_REAL8  , 0, mpicom, ier)
    call mpi_bcast (pctspec , size(pctspec) , MPI_REAL8  , 0, mpicom, ier)
#endif

    ! Determine veg and wtxy for special landunits

    do j = 1,lsmlat
       do i = 1,lsmlon

          veg(i,j,npatch_urban)    = noveg
          wtxy(i,j,npatch_urban)   = pcturb(i,j)/100._r8

          veg(i,j,npatch_lake)     = noveg
          wtxy(i,j,npatch_lake)    = pctlak(i,j)/100._r8

          veg(i,j,npatch_wet)      = noveg
          wtxy(i,j,npatch_wet)     = pctwet(i,j)/100._r8

          veg(i,j,npatch_glacier)  = noveg
          wtxy(i,j,npatch_glacier) = pctgla(i,j)/100._r8

       end do
    end do

  end subroutine surfrd_wtxy_special

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: surfrd_wtxy_veg_rank
!
! !INTERFACE:
  subroutine surfrd_wtxy_veg_rank(ncid, pctspec, veg, wtxy)
!
! !DESCRIPTION:
! Determine wtxy and veg arrays for non-dynamic landuse mode
!
! !USES:
    use clm_varsur, only : all_pfts_on_srfdat	
    use clm_varctl, only : create_crop_landunit
    use domainMod , only : ldomain
!
! !ARGUMENTS:
    implicit none
    include 'netcdf.inc'
    integer , intent(in)    :: ncid                         ! netcdf file id 
    real(r8), intent(in)    :: pctspec(lsmlon, lsmlat)      ! percent wrt gridcell of special landuntis
    integer , intent(inout) :: veg(lsmlon,lsmlat,maxpatch)  ! PFT
    real(r8), intent(inout) :: wtxy(lsmlon,lsmlat,maxpatch) ! subgrid weights
!
! !CALLED FROM:
! subroutine surfrd in this module
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein, Sam Levis and Gordon Bonan
!
!EOP
!
! !LOCAL VARIABLES:
    integer  :: i,j,k,m,k1,k2                            ! indices
    integer  :: dimid,varid                              ! netCDF id's
    integer  :: cropcount                                ! temporary counter
    real(r8) :: sumvec(lsmlon,lsmlat)                    ! temporary vector sum
    logical  :: found                                    ! temporary for error check
    integer  :: iindx, jindx                             ! temporary for error check
    integer  :: miss = 99999                             ! missing data indicator
    real(r8) :: wst(0:numpft)                            ! as pft at specific i, j
    integer  :: wsti(maxpatch_pft)                       ! ranked indices of largest values in wst
    real(r8) :: wst_sum                                  ! sum of %pft
    real(r8) :: sumpct                                   ! sum of %pft over maxpatch_pft
    real(r8) :: diff                                     ! the difference (wst_sum - sumpct)
    real(r8) :: rmax                                     ! maximum patch cover
    integer  :: pft(lsmlon,lsmlat,maxpatch_pft)          ! PFT
    integer  :: cft(lsmlon,lsmlat,maxpatch_cft)          ! CFT
    real(r8) :: pctcft_lunit(lsmlon,lsmlat,maxpatch_cft) ! % of crop landunit area for CFTs
    real(r8) :: pctpft_lunit(lsmlon,lsmlat,maxpatch_pft) ! % of vegetated landunit area for PFTs
    integer  :: ier                                      ! error status
    real(r8) :: pctpft(lsmlon,lsmlat,0:numpft)           ! percent of vegetated gridcell area for PFTs
#if ( defined SCAM )
    integer  :: ret, time_index
    real(r8) :: rmaxpatchdata(maxpatch_pft)
    integer  :: imaxpatchdata(maxpatch_pft)
    real(r8) :: numpftp1data(0:numpft)         
#endif
    character(len=32) :: subname = 'surfrd_wtxy_veg_rank'  ! subroutine name
!-----------------------------------------------------------------------

    if (masterproc) then

       !-----------------------------------------------------------------
       ! Old style surface dataset
       !-----------------------------------------------------------------

       if (.not. all_pfts_on_srfdat) then

          if (create_crop_landunit) then
             write(6,*)'surfrd error: create_crop_landunit is only valid when surface dataset has all pfts'
             call endrun()
          end if

          call check_dim(ncid, 'lsmpft', maxpatch_pft)
#if ( defined SCAM )
          call getncdata (ncid, initLatIdx, initLonIdx, time_index,'PFT', imaxpatchdata, ret)
          pft(1,1,:)=imaxpatchdata(:)
          call getncdata (ncid, initLatIdx, initLonIdx, time_index,'PCT_PFT', rmaxpatchdata, ret)
          pctpft_lunit(1,1,:)=rmaxpatchdata(:)
#else
          call check_ret(nf_inq_varid(ncid, 'PFT', varid), subname)
          call check_ret(nf_get_var_int(ncid, varid, pft), subname)

          call check_ret(nf_inq_varid(ncid, 'PCT_PFT', varid), subname)
          call check_ret(nf_get_var_double(ncid, varid, pctpft_lunit), subname)
#endif

          ! Error check: valid PFTs and sum of cover must equal 100

          sumvec(:,:) = abs(sum(pctpft_lunit,dim=3)-100._r8)
          do j = 1,lsmlat
             do i = 1,lsmlon
                do m = 1, maxpatch_pft
                   if (pft(i,j,m)<0 .or. pft(i,j,m)>numpft) then
                      write(6,*)'surfrd error: invalid PFT for i,j,m=',i,j,m,pft(i,j,m)
                      call endrun()
                   end if
                end do
                if (sumvec(i,j) > 1.e-04_r8 .and. ldomain%pftm(i,j) >= 0) then
                   write(6,*)'surfrd error: PFT cover not equal to 100 for i,j=',i,j
                   do m=1,maxpatch_pft
                      write(6,*)'m= ',m,' pft= ',pft(i,j,m)
                   end do
                   write(6,*)'sumvec= ',sumvec(i,j)
                   call endrun()
                end if
             end do
          end do

       end if

       !-----------------------------------------------------------------
       ! New format surface dataset
       !-----------------------------------------------------------------

       if (all_pfts_on_srfdat) then

          call check_dim(ncid, 'lsmpft', numpft+1)
#if ( defined SCAM )
          call getncdata (ncid, initLatIdx, initLonIdx, time_index,'PCT_PFT', numpftp1data, ret)
          pctpft(1,1,:) = numpftp1data(:)
#else
          call check_ret(nf_inq_varid(ncid, 'PCT_PFT', varid), subname)
          call check_ret(nf_get_var_double(ncid, varid, pctpft), subname)
#endif

          ! 1. pctpft must go back to %vegetated landunit instead of %gridcell
          ! 2. pctpft bare = 100 when landmask = 1 and 100% special landunit
          ! NB: (1) and (2) do not apply to crops.
          ! For now keep all cfts (< 4 anyway) instead of 4 most dominant cfts

          do j = 1,lsmlat
             do i = 1,lsmlon
                cft(i,j,:) = 0
                pctcft_lunit(i,j,:) = 0._r8
                cropcount = 0

                if (pctspec(i,j) < 100._r8) then

                   do m = 0, numpft
                      if (create_crop_landunit) then
                         ! Separate crop landunit is to be created

                         if (crop(m) == 1._r8 .and. pctpft(i,j,m) > 0._r8) then
                            cropcount = cropcount + 1
                            if (cropcount > maxpatch_cft) then
                               write(6,*) 'ERROR surfFileMod: cropcount>maxpatch_cft'
                               call endrun()
                            end if
                            cft(i,j,cropcount) = m
                            pctcft_lunit(i,j,cropcount) = pctpft(i,j,m) * 100._r8/(100._r8-pctspec(i,j))
                            pctpft(i,j,m) = 0.0_r8
                         else if (crop(m) == 0._r8) then
                            pctpft(i,j,m) = pctpft(i,j,m) * 100._r8/(100._r8-pctspec(i,j))
                         end if

                      else
                         ! Separate crop landunit is not created

                         pctpft(i,j,m) = pctpft(i,j,m) * 100._r8/(100._r8-pctspec(i,j))
                      end if
                   end do

                else if (pctspec(i,j) == 100._r8) then

                   pctpft(i,j,0)        = 100._r8
                   pctpft(i,j,1:numpft) =   0._r8

                else

                   write(6,*)subname, 'error: pcturb+pctgla+pctlak+pctwet = ',pctspec(i,j), &
                        ' must be less than or equal to 100'
                   call endrun()

                end if
             end do
          end do

          ! Find pft and pct arrays
          ! Save percent cover by PFT [wst] and total percent cover [wst_sum]

          do j=1,lsmlat
             do i=1,lsmlon

                wst_sum = 0._r8
                sumpct = 0
                do m = 0, numpft
                   wst(m) = pctpft(i,j,m)
                   wst_sum = wst_sum + pctpft(i,j,m)
                end do

                if (ldomain%pftm(i,j) >= 0) then

                   ! Rank [wst] in ascendg order to obtain the top [maxpatch_pft] PFTs

                   call mkrank (numpft, wst, miss, wsti, maxpatch_pft)

                   ! Fill in [pft] and [pctpft] with data for top [maxpatch_pft] PFTs.
                   ! If land model grid cell is ocean, set to no PFTs.
                   ! If land model grid cell is land then:
                   !  1. If [pctlnd_o] = 0, there is no PFT data from the input grid.
                   !     Since need land data, use bare ground.
                   !  2. If [pctlnd_o] > 0, there is PFT data from the input grid but:
                   !     a. use the chosen PFT so long as it is not a missing value
                   !     b. missing value means no more PFTs with cover > 0
                   
                   do m = 1, maxpatch_pft
                      if (wsti(m) /=  miss) then
                         pft(i,j,m) = wsti(m)
                         pctpft_lunit(i,j,m) = wst(wsti(m))
                      else
                         pft(i,j,m) = noveg
                         pctpft_lunit(i,j,m) = 0._r8
                      end if
                      sumpct = sumpct + pctpft_lunit(i,j,m)
                   end do

                else                               ! model grid wants ocean

                   do m = 1, maxpatch_pft
                      pft(i,j,m) = 0
                      pctpft_lunit(i,j,m) = 0._r8
                   end do

                end if

                ! Correct for the case of more than [maxpatch_pft] PFTs present
                
                if (sumpct < wst_sum) then
                   diff  = wst_sum - sumpct
                   sumpct = 0._r8
                   do m = 1, maxpatch_pft
                      pctpft_lunit(i,j,m) = pctpft_lunit(i,j,m) + diff/maxpatch_pft
                      sumpct = sumpct + pctpft_lunit(i,j,m)
                   end do
                end if

                ! Error check: make sure have a valid PFT

                do m = 1,maxpatch_pft
                   if (pft(i,j,m) < 0 .or. pft(i,j,m) > numpft) then
                      write (6,*)'surfrd error: invalid PFT at gridcell i,j=',i,j,pft(i,j,m)
                      call endrun()
                   end if
                end do

                ! As done in mksrfdatMod.F90 for other percentages, truncate pctpft to
                ! ensure that weight relative to landunit is not nonzero
                ! (i.e. a very small number such as 1e-16) where it really should be zero
                ! The following if-block is here to preserve roundoff level differences
                ! between the call to surfrd_wtxy_veg_all and surfrd_wtxy_veg_rank

                if (maxpatch_pft < numpft+1) then
                   do m=1,maxpatch_pft
                      pctpft_lunit(i,j,m) = float(nint(pctpft_lunit(i,j,m)))
                   end do
                   do m=1,maxpatch_cft
                      pctcft_lunit(i,j,m) = float(nint(pctcft_lunit(i,j,m)))
                   end do
                end if
                   
                ! Make sure sum of PFT cover == 100 for land points. If not,
                ! subtract excess from most dominant PFT.

                rmax = -9999._r8
                k1 = -9999
                k2 = -9999
                sumpct = 0._r8
                do m = 1, maxpatch_pft
                   sumpct = sumpct + pctpft_lunit(i,j,m)
                   if (pctpft_lunit(i,j,m) > rmax) then
                      k1 = m
                      rmax = pctpft_lunit(i,j,m)
                   end if
                end do
                do m = 1, maxpatch_cft
                   sumpct = sumpct + pctcft_lunit(i,j,m)
                   if (pctcft_lunit(i,j,m) > rmax) then
                      k2 = m
                      rmax = pctcft_lunit(i,j,m)
                   end if
                end do
                if (k1 == -9999 .and. k2 == -9999) then
                   write(6,*)'surfrd error: largest PFT patch not found'
                   call endrun()
                else if (ldomain%pftm(i,j) >= 0) then
                   if (sumpct < 95 .or. sumpct > 105._r8) then
                      write(6,*)'surfrd error: sum of PFT cover =',sumpct,' at i,j=',i,j
                      call endrun()
                   else if (sumpct /= 100._r8 .and. k2 /= -9999) then
                      pctcft_lunit(i,j,k2) = pctcft_lunit(i,j,k2) - (sumpct-100._r8)
                   else if (sumpct /= 100._r8) then
                      pctpft_lunit(i,j,k1) = pctpft_lunit(i,j,k1) - (sumpct-100._r8)
                   end if
                end if

                ! Error check: make sure PFTs sum to 100% cover

                sumpct = 0._r8
                do m = 1, maxpatch_pft
                   sumpct = sumpct + pctpft_lunit(i,j,m)
                end do
                do m = 1, maxpatch_cft
                   sumpct = sumpct + pctcft_lunit(i,j,m)
                end do
                if (ldomain%pftm(i,j) >= 0) then
                   if (abs(sumpct - 100._r8) > 0.000001_r8) then
                      write(6,*)'surfFileMod error: sum(pct) over maxpatch_pft is not = 100.'
                      write(6,*)sumpct, i,j
                      call endrun()
                   end if
                   if (sumpct < -0.000001_r8) then
                      write(6,*)'surfFileMod error: sum(pct) over maxpatch_pft is < 0.'
                      write(6,*)sumpct, i,j
                      call endrun()
                   end if
                end if

             end do   ! end of longitude loop
          end do   ! end of latitude loop

       end if

    endif   ! end of if-masterproc block

#if (defined SPMD)
    call mpi_bcast (pft         , size(pft)         , MPI_INTEGER, 0, mpicom, ier)
    call mpi_bcast (pctpft_lunit, size(pctpft_lunit), MPI_REAL8  , 0, mpicom, ier)
    call mpi_bcast (cft         , size(cft)         , MPI_INTEGER, 0, mpicom, ier)
    call mpi_bcast (pctcft_lunit, size(pctcft_lunit), MPI_REAL8  , 0, mpicom, ier)
#endif

    ! Determine array [veg], which sets the PFT type for each of the [maxpatch]
    ! patches and array [wtxy], which sets the relative abundance of the PFT.
    ! Fill in PFTs for vegetated portion of grid cell. Fractional areas for
    ! these points [pctpft] pertain to "vegetated" area not to total grid area.
    ! So need to adjust them for fraction of grid that is vegetated.
    ! Next, fill in urban, lake, wetland, and glacier patches.

    do j = 1,lsmlat
       do i = 1,lsmlon
          if (ldomain%pftm(i,j) >= 0) then

             ! Naturally vegetated landunit

             do m = 1, maxpatch_pft
                veg(i,j,m)  = pft(i,j,m)
                wtxy(i,j,m) = pctpft_lunit(i,j,m) * (100._r8-pctspec(i,j))/10000._r8
#if (defined CN)
                ! the following test prevents the assignment of temperate deciduous
                ! vegetation types in the tropics
                ! 1. broadleaf deciduous temperate tree -> broadleaf deciduous tropical tree

                if (veg(i,j,m) == 7 .and. abs(ldomain%latc(i,j)) < 23.5_r8) veg(i,j,m) = 6

                ! 2. broadleaf deciduous temperate shrub -> broadleaf deciduous tropical tree
                ! this reassignment from shrub to tree is necessary because there is currently no
                ! tropical deciduous broadleaf shrub type defined.

                 if (veg(i,j,m) == 10 .and. abs(ldomain%latc(i,j)) < 23.5_r8) veg(i,j,m) = 6
#endif
             end do

             ! Crop landunit

             if (create_crop_landunit) then
                do m = 1,maxpatch_cft
                   veg(i,j,npatch_glacier+m)  = cft(i,j,m)
                   wtxy(i,j,npatch_glacier+m) = pctcft_lunit(i,j,m) * (100._r8-pctspec(i,j))/10000._r8
                end do
             end if

          end if
       end do
    end do

    ! Error check

    found = .false.
    sumvec(:,:) = abs(sum(wtxy,dim=3)-1._r8)
    do j=1,lsmlat
       do i=1,lsmlon
          if (sumvec(i,j) > 1.e-06_r8 .and. ldomain%pftm(i,j)>=0) then
             found = .true.
             iindx = i
             jindx = j
             exit
          endif
       end do
       if (found) exit
    end do
    if ( found ) then
       write (6,*)'surfrd error: WTXY > 1 occurs at i,j= ',iindx,jindx; call endrun()
    end if

  end subroutine surfrd_wtxy_veg_rank

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: surfrd_wtxy_veg_all
!
! !INTERFACE:
  subroutine surfrd_wtxy_veg_all(ncid, pctspec, veg, wtxy)
!
! !DESCRIPTION:
! Determine wtxy and veg arrays for non-dynamic landuse mode
!
! !USES:
    use domainMod, only : ldomain
!
! !ARGUMENTS:
    implicit none
    include 'netcdf.inc'
    integer , intent(in)    :: ncid                         ! netcdf file id 
    real(r8), intent(in)    :: pctspec(lsmlon, lsmlat)      ! percent wrt gridcell of special landuntis
    integer , intent(inout) :: veg(lsmlon,lsmlat,maxpatch)  ! PFT
    real(r8), intent(inout) :: wtxy(lsmlon,lsmlat,maxpatch) ! subgrid weights
!
! !CALLED FROM:
! subroutine surfrd in this module
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein, Sam Levis and Gordon Bonan
!
!EOP
!
! !LOCAL VARIABLES:
    integer  :: i,j,m,mp7,mp8,mp11             ! indices
    integer  :: dimid,varid                    ! netCDF id's
    integer  :: ier                            ! error status	
    real(r8) :: sumpct                         ! sum of %pft over maxpatch_pft
    real(r8) :: pctpft(lsmlon,lsmlat,0:numpft) ! percent of vegetated gridcell area for PFTs
#if (defined SCAM)
    integer  :: ret, time_index
    integer  :: imaxpatchdata(maxpatch_pft)
    real(r8) :: rmaxpatchdata(maxpatch_pft)
    real(r8) :: numpftp1data(0:numpft)         
#endif
    character(len=32) :: subname = 'surfrd_wtxy_veg_all'  ! subroutine name
!-----------------------------------------------------------------------

    if (masterproc) then
       call check_dim(ncid, 'lsmpft', numpft+1)
#if ( defined SCAM )
       call getncdata (ncid, initLatIdx, initLonIdx, time_index,'PCT_PFT', numpftp1data, ret)
       pctpft(1,1,:) = numpftp1data(:)
#else
       call check_ret(nf_inq_varid(ncid, 'PCT_PFT', varid), subname)
       call check_ret(nf_get_var_double(ncid, varid, pctpft), subname)
#endif
    end if

#if (defined SPMD)
    call mpi_bcast (pctpft, size(pctpft), MPI_REAL8, 0, mpicom, ier)
#endif

    do j = 1,lsmlat
       do i = 1,lsmlon
          if (ldomain%pftm(i,j) >= 0) then

             ! Error check: make sure PFTs sum to 100% cover for vegetated landunit 
             ! (convert pctpft from percent with respect to gridcel to percent with 
             ! respect to vegetated landunit)

             if (pctspec(i,j) < 100._r8) then
                sumpct = 0._r8
                do m = 0,numpft
                   sumpct = sumpct + pctpft(i,j,m) * 100._r8/(100._r8-pctspec(i,j))
                end do
                 if (abs(sumpct - 100._r8) > 0.1e-4_r8) then
                   write(6,*)'surfFileMod error: sum(pct) over numpft+1 is not = 100.'
                   write(6,*) sumpct, sumpct-100._r8, i, j
                   call endrun()
                end if
                if (sumpct < -0.000001_r8) then
                   write(6,*)'surfFileMod error: sum(pct) over numpft+1 is < 0.'
                   write(6,*) sumpct, i,j
                   call endrun()
                end if
             end if

             ! Set weight of each pft wrt gridcell (note that maxpatch_pft = numpft+1 here)

             do m = 1,numpft+1
                veg(i,j,m)  = m-1
                wtxy(i,j,m) = pctpft(i,j,m-1) / 100._r8
             end do

#if (defined CN)
             ! the following test prevents the assignment of temperate deciduous
             ! vegetation types in the tropics
             ! 1. broadleaf deciduous temperate tree (type7) -> broadleaf deciduous tropical tree (type6)
			 ! N.B. the veg and wtxy arrays start at 1, so index 1 corresponds to
             ! veg type 0.  So in this case I want to trap on veg types 7 and 10, 
             ! which are indices 8 and 11. Moving to vegtype6, or index 7
             mp7 = 7
             mp8 = 8
             mp11 = 11
             if (abs(ldomain%latc(i,j)) < 23.5_r8 .and. wtxy(i,j,mp8) > 0._r8) then
             	 if (masterproc) then
                   write(6,*)'surfFileMod warning: reassigning temperate tree -> tropical tree'
                   write(6,*)'i,j,lat,veg7wt,veg6wt,type'
                   write(6,*) i,j,ldomain%latc(i,j),wtxy(i,j,mp8),wtxy(i,j,mp7),veg(i,j,mp8)
                end if
             	 wtxy(i,j,mp7) = wtxy(i,j,mp7) + wtxy(i,j,mp8)
                wtxy(i,j,mp8) = 0._r8
             	 if (masterproc) then
                   write(6,*) i,j,ldomain%latc(i,j),wtxy(i,j,mp8),wtxy(i,j,mp7)
                end if
             end if

             ! 2. broadleaf deciduous temperate shrub (type10) -> broadleaf deciduous tropical tree (type6)
             ! this reassignment from shrub to tree is necessary because there is currently no
             ! tropical deciduous broadleaf shrub type defined.

             if (abs(ldomain%latc(i,j)) < 23.5_r8 .and. wtxy(i,j,mp11) > 0._r8) then
             	 if (masterproc) then
                   write(6,*)'surfFileMod warning: reassigning temperate shrub -> tropical tree'
                   write(6,*)'i,j,lat,veg10wt,veg6wt,type'
                   write(6,*) i,j,ldomain%latc(i,j),wtxy(i,j,mp11),wtxy(i,j,mp7),veg(i,j,mp11)
                end if
             	 wtxy(i,j,mp7) = wtxy(i,j,mp7) + wtxy(i,j,mp11)
                wtxy(i,j,mp11) = 0._r8
             	 if (masterproc) then
                   write(6,*) i,j,ldomain%latc(i,j),wtxy(i,j,mp11),wtxy(i,j,mp7)
                end if
             end if
#endif
          end if
       end do
    end do

  end subroutine surfrd_wtxy_veg_all

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: surfrd_wtxy_veg_dgvm
!
! !INTERFACE:
  subroutine surfrd_wtxy_veg_dgvm(pctspec, veg, wtxy)
!
! !DESCRIPTION:
! Determine wtxy and vegxy for DGVM mode.
!
! !USES:
!
! !ARGUMENTS:
    implicit none
    real(r8), intent(in)    :: pctspec(lsmlon,lsmlat)        ! percent wrt gridcell of  special landunits
    integer , intent(inout) :: veg(lsmlon,lsmlat,maxpatch)   ! PFT
    real(r8), intent(inout) :: wtxy(lsmlon,lsmlat,maxpatch)  ! subgrid weights
!
! !CALLED FROM:
! subroutine surfrd in this module
!
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein 12/04
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: i,j,m  ! indices
!-----------------------------------------------------------------------

    do j = 1,lsmlat
       do i = 1,lsmlon
          do m = 1, maxpatch_pft
             veg(i,j,m)  = noveg 
             wtxy(i,j,m) = 1.0_r8/maxpatch_pft * (100._r8-pctspec(i,j))/100._r8
          end do
       end do
    end do

  end subroutine surfrd_wtxy_veg_dgvm
   
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: mkrank
!
! !INTERFACE:
  subroutine mkrank (n, a, miss, iv, num)
!
! !DESCRIPTION:
! Return indices of largest [num] values in array [a]
!
! !USES:
    use shr_kind_mod, only: r8 => shr_kind_r8
    use abortutils, only : endrun
!
! !ARGUMENTS:
    implicit none
    integer , intent(in) :: n        ! array length
    real(r8), intent(in) :: a(0:n)   ! array to be ranked
    integer , intent(in) :: miss     ! missing data value
    integer , intent(in) :: num      ! number of largest values requested
    integer , intent(out):: iv(num)  ! index to [num] largest values in array [a]
!
! !CALLED FROM:
! ! subroutine surfrd_wtxy_veg_rank in this module
!
! !REVISION HISTORY:
! Author: Gordon Bonan
!
!EOP
!
! !LOCAL VARIABLES:
    real(r8) :: a_max       ! maximum value in array
    integer  :: i           ! array index
    real(r8) :: delmax      ! tolerance for finding if larger value
    integer  :: m           ! do loop index
    integer  :: k           ! do loop index
    logical  :: exclude     ! true if data value has already been chosen
!-----------------------------------------------------------------------

    delmax = 1.e-06_r8

    ! Find index of largest non-zero number
    
    iv(1) = miss
    a_max = -9999._r8

    do i = 0, n
       if (a(i)>0._r8 .and. (a(i)-a_max)>delmax) then
          a_max = a(i)
          iv(1)  = i
       end if
    end do

    ! iv(1) = miss indicates no values > 0. this is an error

    if (iv(1) == miss) then
       write (6,*) 'MKRANK error: iv(1) = missing'
       call endrun
    end if

    ! Find indices of the next [num]-1 largest non-zero number.
    ! iv(m) = miss if there are no more values > 0

    do m = 2, num
       iv(m) = miss
       a_max = -9999._r8
       do i = 0, n

          ! exclude if data value has already been chosen

          exclude = .false.
          do k = 1, m-1
             if (i == iv(k)) exclude = .true.
          end do

          ! if not already chosen, see if it is the largest of
          ! the remaining values

          if (.not. exclude) then
             if (a(i)>0._r8 .and. (a(i)-a_max)>delmax) then
                a_max = a(i)
                iv(m)  = i
             end if
          end if
       end do
    end do

  end subroutine mkrank

end module surfFileMod
