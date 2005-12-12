#include <misc.h>
#include <preproc.h>

module decompMod
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: decompMod
!
! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
#ifdef SPMD
#ifdef COUP_CAM
  use spmd_utils  , only : masterproc, iam, npes
#else
  use spmdMod     , only : masterproc, iam, npes
#endif
#else
  use spmdMod     , only : masterproc, iam, npes
#endif
  use shr_sys_mod , only : shr_sys_flush
  use abortutils  , only : endrun
!
! !PUBLIC TYPES:
  implicit none
  integer, public :: clump_pproc ! number of clumps per MPI process
!
! !PUBLIC MEMBER FUNCTIONS:
  public initDecomp              ! initializes land surface decomposition
                                 ! into clumps and processors
  public get_nclumps             ! returns the number of clumps defined
  public get_clump_cell_id_coord ! returns clump/cell ids based on lon/lat
  public get_clump_owner_id      ! returns clump owner based on clump id
  public get_clump_ncells_proc   ! returns number of cells for process
  public get_clump_ncells_id     ! returns number of cells in clump
  public get_clump_coord_id      ! returns lon/lat coordinates based on id
  public get_clump_gcell_info    ! returns 1d gridcell index
  public get_clump_bounds        ! beg and end gridcell, landunit, column,
                                 ! pft indices for clump
  public get_proc_clumps         ! number of clumps for this processor
  public get_proc_bounds         ! beg and end gridcell, landunit, column,
                                 ! pft indices for this processor
  public get_proc_total          ! total number of gridcells, landunits,
                                 ! columns and pfts for any processor
  public get_proc_global         ! total gridcells, landunits, columns, pfts
                                 ! across all processors
  save
!
! !DESCRIPTION:
! Module provides a descomposition into a clumped data structure which can
! be mapped back to atmosphere physics chunks.
!
! !REVISION HISTORY:
! 2002.09.11  Forrest Hoffman  Creation.
!
!EOP
!
! !PRIVATE TYPES:
  private

  integer :: nclumps             ! total number of clumps across all processors

  type processor
     integer :: nclumps          ! number of clumps for processor iam
     integer, pointer :: cid(:)  ! array of clump id's for processor iam
     integer :: begg, endg       ! beginning and ending gridcell index on processor iam
     integer :: begl, endl       ! beginning and ending landunit index on processor iam
     integer :: begc, endc       ! beginning and ending column index on processor iam
     integer :: begp, endp       ! beginning and ending pft index on processor iam
     integer :: ncells           ! total number of gridcells on processor iam
     integer :: nlunits          ! total number of landunits on processor iam
     integer :: ncols            ! total number of columns on processor iam
     integer :: npfts            ! total number of pfts on processor iam
     integer :: numg             ! total number of gridcells across all processors
     integer :: numl             ! total number of landunits across all processors
     integer :: numc             ! total number of columns across all processors
     integer :: nump             ! total number of pfts across all processors
  end type processor
  type(processor), allocatable :: procs(:)

  type clump
     integer :: owner            ! process id owning clump
     integer :: ncells           ! number of gridcells in clump
     integer :: nlunits          ! number of landunits in clump
     integer :: ncols            ! number of columns in clump
     integer :: npfts            ! number of pfts in clump
     integer :: begg, endg       ! beginning and ending gridcell index in clump
     integer :: begl, endl       ! beginning and ending landunit index in clump
     integer :: begc, endc       ! beginning and ending column index in clump
     integer :: begp, endp       ! beginning and ending pft index in clump
     integer, pointer :: ixy(:)  ! cell longitude indices
     integer, pointer :: jxy(:)  ! cell latitude indices
     integer, pointer :: gi(:)   ! global 1d grid cell index
  end type clump
  type(clump), allocatable :: clumps(:)

  type pmulc
     integer :: clumpid          ! clump id for (lon, lat)
     integer :: cell             ! matching cell id
  end type pmulc
  type(pmulc), allocatable :: pmulcs(:,:)

  type gcell_decomp
     integer :: gsn     ! corresponding cell index in south->north gridcell array
     integer :: ixy     ! cell longitude index
     integer :: jxy     ! cell latitude index
     integer :: li      ! beginning landunit index
     integer :: lf      ! ending landunit index
     integer :: ci      ! beginning column index
     integer :: cf      ! ending column index
     integer :: pi      ! beginning pft index
     integer :: pf      ! ending pft index
  end type gcell_decomp
  public gcell_decomp	
  type(gcell_decomp), public, allocatable :: gcelldc(:)

  type gcell_south_north
     integer :: gdc     ! corresponding cell index in decomposition gridcell array
     integer :: ixy     ! cell longitude index
     integer :: jxy     ! cell latitude index
     integer :: li      ! beginning landunit index
     integer :: lf      ! ending landunit index
     integer :: ci      ! beginning column index
     integer :: cf      ! ending column index
     integer :: pi      ! beginning pft index
     integer :: pf      ! ending pft index
  end type gcell_south_north
  public gcell_south_north
  type(gcell_south_north), public, allocatable :: gcellsn(:)

  interface map_dc2sn
     module procedure map_dc2sn_sl_real
     module procedure map_dc2sn_sl_int
     module procedure map_dc2sn_ml1_real
     module procedure map_dc2sn_ml1_int
  end interface
  public map_dc2sn

  interface map_sn2dc
     module procedure map_sn2dc_sl_real
     module procedure map_sn2dc_sl_int
     module procedure map_sn2dc_ml1_real
     module procedure map_sn2dc_ml1_int
  end interface
  public map_sn2dc
!------------------------------------------------------------------------------
! $Id$
! $Author$
!------------------------------------------------------------------------------

contains

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: initDecomp
!
! !INTERFACE:
  subroutine initDecomp(wtxy)
!
! !DESCRIPTION:
! This subroutine initializes the land surface decomposition into a clump
! data structure.
!
! !USES:
    use clmtype
    use clm_varsur    , only : numlon, landmask
    use clm_varpar    , only : lsmlon, lsmlat, maxpatch
    use initSubgridMod, only : get_gcell_info
!
! !ARGUMENTS:
    implicit none
    real(r8), intent(in) :: wtxy(lsmlon, lsmlat, maxpatch) ! subgrid patch
                                                           ! weights
!
! !LOCAL VARIABLES:
    integer :: ppc                    ! min number of pfts per clump
    integer :: lpc                    ! min number of landunits per clump
    integer :: ppclump                ! min pfts per clump
    integer :: i,j,cid,pid            ! indices
    integer :: gi,li,ci,pi            ! indices
    integer :: gf,lf,cf,pf            ! indices
    integer :: g,l,c,p,n,m            ! indices
    integer :: gdc,gsn                ! indices
    integer :: nzero                  ! first clump with zero gridcells
    integer :: ncells                 ! total gridcells
    integer :: nlunits                ! total landunits
    integer :: ncols                  ! total columns
    integer :: npfts                  ! total pfts
    integer :: nveg                   ! number of pfts in vegetated landunit
    integer :: numg                   ! total number of gridcells across all
                                      ! processors
    integer :: numl                   ! total number of landunits across all
                                      ! processors
    integer :: numc                   ! total number of columns across all
                                      ! processors
    integer :: nump                   ! total number of pfts across all
                                      ! processors
    logical, pointer :: clumpfull(:)  ! true => clump is full
    logical :: validclump             ! temporary for finding full clump
    logical :: error = .false.        ! temporary for finding full clump
    integer :: clumpcount             ! temporary for finding full clump
    integer :: ilunits, icols, ipfts  ! temporaries
    integer :: ng                     ! temporaries
    integer :: nl                     ! temporaries
    integer :: nc                     ! temporaries
    integer :: np                     ! temporaries
    integer :: ier                    ! error code
!
! !CALLED FROM:
! subroutine initialize
!
! !REVISION HISTORY:
! 2002.09.11  Forrest Hoffman  Creation.
!
!EOP
!------------------------------------------------------------------------------

    ! Dynamic memory allocation for procs

    ier = 0
    if( .not. allocated(procs)) allocate(procs(0:npes-1), stat=ier)
    if (ier /= 0) then
       write (6,*) 'initDecomp(): allocation error for procs'
       call endrun
    end if

    ! Find total global number of grid cells, landunits, columns and pfts

    ncells = 0
    nlunits = 0
    ncols = 0
    npfts = 0
    do j = 1, lsmlat
       do i = 1, numlon(j)
          if (landmask(i,j) == 1) then
             call get_gcell_info (i, j, wtxy, nlunits=ilunits, ncols=icols, npfts=ipfts)
             ncells  = ncells  + 1
             nlunits = nlunits + ilunits
             ncols   = ncols   + icols
             npfts   = npfts   + ipfts
          end if
       end do
    end do

    ! Allocate dynamic memory for gridcell derived types

    if( .not. allocated(gcellsn)) allocate(gcellsn(ncells), gcelldc(ncells), stat=ier)
    if (ier /= 0) then
       write (6,*) 'initDecomp(): allocation error for gcellsn and gcelldc'
       call endrun()
    end if

    ! Error check on total number of gridcells

    if (npes > ncells) then
       write (6,*) 'initDecomp(): Number of processes exceeds number ', &
            'of land grid cells'
       call endrun
    end if

    ! Diagnostic output

    if (masterproc) then
       write (6,*)' Surface Grid Characteristics'
       write (6,*)'   longitude points          = ',lsmlon
       write (6,*)'   latitude points           = ',lsmlat
       write (6,*)'   minimum number of longitude points per latitude = ',minval(numlon)
       write (6,*)'   maximum number of longitude points per latitude = ',maxval(numlon)
       write (6,*)'   total number of gridcells = ',ncells
       write (6,*)'   total number of landunits = ',nlunits
       write (6,*)'   total number of columns   = ',ncols
       write (6,*)'   total number of pfts      = ',npfts
       write (6,*)' Decomposition Characteristics'
       write (6,*)'   clumps per process        = ',clump_pproc
       write (6,*)
    end if

    ! Determine number of gridcell clumps and allocate dynamic memory.
    ! Decompose by a fixed number of clumps per MPI process.

    if (clump_pproc > 0) then
       nclumps = clump_pproc * npes
       if (nclumps < npes) then
          write (6,*) 'initDecomp(): Number of gridcell clumps= ',nclumps, &
               ' is less than the number of processes = ', npes
          call endrun()
       end if
       ppc = npfts / nclumps
       if (ppc * nclumps < npfts) ppc = ppc + 1
    else
       write(6,*)'clump_pproc= ',clump_pproc,'  must be greater than 0'
       call endrun()
    end if

    if ( .not. allocated(clumps)) allocate(clumps(1:nclumps), stat=ier)
    if (ier /= 0) then
       write (6,*) 'initDecomp(): allocation error for clumps'
       call endrun()
    end if
    allocate(clumpfull(1:nclumps), stat=ier)
    if (ier /= 0) then
       write (6,*) 'initDecomp(): allocation error for temporary clump arrays'
       call endrun()
    end if
    if( .not. allocated(pmulcs)) allocate(pmulcs(1:lsmlon, 1:lsmlat), stat=ier)
    if (ier /= 0) then
       write (6,*) 'initDecomp(): allocation error for pmulcs'
       call endrun()
    end if

    ! Determine number of grid cells in each clump - assign gridcells to clumps
    ! in a round robin fashion

    clumps(:)%ncells  = 0
    clumps(:)%nlunits = 0
    clumps(:)%ncols   = 0
    clumps(:)%npfts   = 0
    clumpfull(:) = .false.
    cid = 0
    do j = 1,lsmlat
       do i = 1,numlon(j)
          if (landmask(i,j) == 1) then

             ! Set clump id for gridcell - if clump is full, then determine
             ! first non-full clump

             cid = cid + 1
             if (cid == nclumps+1) cid = 1
             if (clumpfull(cid)) then
                clumpcount = 0
                validclump = .false.
                do while (.not. validclump .and. clumpcount <= nclumps)
                   clumpcount = clumpcount + 1
                   cid = cid + 1
                   if (cid == nclumps+1) cid = 1
                   if (.not. clumpfull(cid)) validclump = .true.
                end do
                if (.not. validclump) error = .true.
                if (error) exit
             end if

            ! Determine grid cell info

             call get_gcell_info (i, j, wtxy, nlunits=ilunits, ncols=icols, npfts=ipfts)

             clumps(cid)%ncells  = clumps(cid)%ncells  + 1
             clumps(cid)%nlunits = clumps(cid)%nlunits + ilunits
             clumps(cid)%ncols   = clumps(cid)%ncols   + icols
             clumps(cid)%npfts   = clumps(cid)%npfts   + ipfts

             if (clumps(cid)%npfts >= ppc .and. cid < nclumps) then
                clumpfull(cid) = .true.
             end if
          end if
       end do
       if (error) exit
    end do
    if (error) then
       write(6,*)'initDecomp: error encountered while trying to find an unfull clump'
       call endrun()
    end if

    do cid = 1,nclumps
       if (clumps(cid)%ncells == 0   .or. clumps(cid)%nlunits == 0 .or. &
           clumps(cid)%ncols == 0 .or. clumps(cid)%npfts == 0) then
          write(6,*)'Invalid clump at clump index ',cid
          write(6,*)'clumps(cid)%ncells  = ',clumps(cid)%ncells
          write(6,*)'clumps(cid)%nlunits = ',clumps(cid)%nlunits
          write(6,*)'clumps(cid)%ncols   = ',clumps(cid)%ncols
          write(6,*)'clumps(cid)%npfts   = ',clumps(cid)%npfts
          call endrun
       end if
    end do

    ! Allocate dynamic memory for clumps

    do cid = 1,nclumps
       ncells  = clumps(cid)%ncells
       nlunits = clumps(cid)%nlunits
       ncols   = clumps(cid)%ncols
       npfts   = clumps(cid)%npfts
       allocate(clumps(cid)%ixy(ncells), clumps(cid)%jxy(ncells), clumps(cid)%gi(ncells), stat=ier)
       if (ier /= 0) then
          write (6,*) 'initDecomp(): allocation errors for clump indices'
          call endrun()
       end if
    end do

    ! Redo the above calculation for gridcells distribution of clumps to determine
    ! the xy i and j indices for each clump gridcell

    clumps(:)%ncells  = 0
    clumps(:)%nlunits = 0
    clumps(:)%ncols   = 0
    clumps(:)%npfts   = 0
    clumpfull(:) = .false.
    cid = 0
    do j = 1,lsmlat
       do i = 1,numlon(j)
          if (landmask(i,j) == 1) then

             ! Set clump id for gridcell - if clump is full, then determine
             ! first non-full clump

             cid = cid + 1
             if (cid == nclumps+1) cid = 1
             if (clumpfull(cid)) then
                clumpcount = 0
                validclump = .false.
                do while (.not. validclump .and. clumpcount <= nclumps)
                   clumpcount = clumpcount + 1
                   cid = cid + 1
                   if (cid == nclumps+1) cid = 1
                   if (.not. clumpfull(cid)) validclump = .true.
                end do
                if (.not. validclump) error = .true.
                if (error) exit
             end if

            ! Determine grid cell info

             call get_gcell_info (i, j, wtxy, nlunits=ilunits, ncols=icols, npfts=ipfts)

             clumps(cid)%ncells  = clumps(cid)%ncells  + 1
             clumps(cid)%nlunits = clumps(cid)%nlunits + ilunits
             clumps(cid)%ncols   = clumps(cid)%ncols   + icols
             clumps(cid)%npfts   = clumps(cid)%npfts   + ipfts

             if (clumps(cid)%npfts >= ppc .and. cid < nclumps) then
                clumpfull(cid) = .true.
             end if

             ! Now determine the ixy and jxy indices since memory has been
             ! allocated for these arrays

             n = clumps(cid)%ncells
             clumps(cid)%ixy(n) = i
             clumps(cid)%jxy(n) = j

          end if
       end do
       if (error) exit
    end do
    if (error) then
       write(6,*)'initDecomp: error encountered while trying to find an unfull clump'
       call endrun()
    end if

    ! Assign clumps to processors in a round robin fashion.

    procs(:)%ncells  = 0
    procs(:)%nlunits = 0
    procs(:)%ncols   = 0
    procs(:)%npfts   = 0
    procs(:)%nclumps = 0
    pid = 0
    do cid = 1,nclumps
       clumps(cid)%owner  = pid
       procs(pid)%nclumps = procs(pid)%nclumps + 1
       procs(pid)%ncells  = procs(pid)%ncells  + clumps(cid)%ncells
       procs(pid)%nlunits = procs(pid)%nlunits + clumps(cid)%nlunits
       procs(pid)%ncols   = procs(pid)%ncols   + clumps(cid)%ncols
       procs(pid)%npfts   = procs(pid)%npfts   + clumps(cid)%npfts
       pid = pid + 1
       if (pid == npes) pid = 0
    end do

    ! Determine clump ids for each clump on this processor

    do pid = 0,npes-1
       n = procs(pid)%nclumps
       allocate(procs(pid)%cid(n), stat=ier)
       if (ier /= 0) then
          write (6,*) 'initProc(): allocation error for procs(pid)%cid'
          call endrun()
       end if
       n = 1
       do cid = 1,nclumps
          if (clumps(cid)%owner == pid) then
             procs(pid)%cid(n) = cid
             n = n + 1
          end if
       end do
    end do

    ! Determine per processor indices

    procs(0)%begg = 1
    procs(0)%begl = 1
    procs(0)%begc = 1
    procs(0)%begp = 1
    procs(0)%endg = procs(0)%ncells
    procs(0)%endl = procs(0)%nlunits
    procs(0)%endc = procs(0)%ncols
    procs(0)%endp = procs(0)%npfts
    do pid = 1,npes-1
       procs(pid)%begg = procs(pid-1)%endg + 1
       procs(pid)%begl = procs(pid-1)%endl + 1
       procs(pid)%begc = procs(pid-1)%endc + 1
       procs(pid)%begp = procs(pid-1)%endp + 1
       procs(pid)%endg = procs(pid-1)%endg + procs(pid)%ncells
       procs(pid)%endl = procs(pid-1)%endl + procs(pid)%nlunits
       procs(pid)%endc = procs(pid-1)%endc + procs(pid)%ncols
       procs(pid)%endp = procs(pid-1)%endp + procs(pid)%npfts
    end do

    ! Determine per clump indices

    gi = 1
    li = 1
    ci = 1
    pi = 1
    gf = clumps(1)%ncells
    lf = clumps(1)%nlunits
    cf = clumps(1)%ncols
    pf = clumps(1)%npfts
    do pid = 0,npes-1
       do n = 1,procs(pid)%nclumps
          cid = procs(pid)%cid(n)
          if (pid == 0 .and. n == 0 .and. cid /= 1) then
             write(6,*)'initProc error: clump 1 must always be the first clump ',&
                  ' on the first processor'
             call endrun()
          end if
          if (cid /= 1) then
             gi = gf + 1
             li = lf + 1
             ci = cf + 1
             pi = pf + 1
             gf = gi + clumps(cid)%ncells  - 1
             lf = li + clumps(cid)%nlunits - 1
             cf = ci + clumps(cid)%ncols   - 1
             pf = pi + clumps(cid)%npfts   - 1
          end if
          clumps(cid)%begg = gi
          clumps(cid)%begl = li
          clumps(cid)%begc = ci
          clumps(cid)%begp = pi
          clumps(cid)%endg = gf
          clumps(cid)%endl = lf
          clumps(cid)%endc = cf
          clumps(cid)%endp = pf
       end do
    end do

    ! Determine information across all processors

    numg = procs(0)%ncells
    numl = procs(0)%nlunits
    numc = procs(0)%ncols
    nump = procs(0)%npfts
    do pid = 1,npes-1
       numg = numg + procs(pid)%ncells
       numl = numl + procs(pid)%nlunits
       numc = numc + procs(pid)%ncols
       nump = nump + procs(pid)%npfts
    end do
    procs(0:npes-1)%numg = numg
    procs(0:npes-1)%numl = numl
    procs(0:npes-1)%numc = numc
    procs(0:npes-1)%nump = nump

    ! Determine pmulcs components and index into decomposition gridcells for each
    ! clump gridcell

    pmulcs(:,:)%clumpid = 0
    pmulcs(:,:)%cell = 0
    g = 0
    do pid = 0,npes-1
       do nc = 1,procs(pid)%nclumps
          cid = procs(pid)%cid(nc)
          do n = 1,clumps(cid)%ncells
             g = g + 1
             i = clumps(cid)%ixy(n)
             j = clumps(cid)%jxy(n)
             pmulcs(i,j)%clumpid = cid
             pmulcs(i,j)%cell = n
             clumps(cid)%gi(n) = g
          end do
       end do
    end do

    ! Determine derived type components for south->north gridcell array
    ! and decomposition gridcell array

    gcellsn(1)%li = 1
    gcellsn(1)%ci = 1
    gcellsn(1)%pi = 1
    g = 0
    do j = 1,lsmlat
       do i = 1,numlon(j)
          if (landmask(i,j) == 1) then
             call get_gcell_info (i, j, wtxy, nlunits=ilunits, ncols=icols, npfts=ipfts)
             g = g+1
             gcellsn(g)%ixy = i
             gcellsn(g)%jxy = j
             gcellsn(g)%lf = gcellsn(g)%li + ilunits - 1
             gcellsn(g)%cf = gcellsn(g)%ci + icols   - 1
             gcellsn(g)%pf = gcellsn(g)%pi + ipfts   - 1
             if (g <= numg-1) then
                gcellsn(g+1)%li = gcellsn(g)%lf + 1
                gcellsn(g+1)%ci = gcellsn(g)%cf + 1
                gcellsn(g+1)%pi = gcellsn(g)%pf + 1
             end if
          end if
       end do
    end do

    gcelldc(1)%li = 1
    gcelldc(1)%ci = 1
    gcelldc(1)%pi = 1
    g = 0
    do pid = 0,npes-1
       do nc = 1,procs(pid)%nclumps
          cid = procs(pid)%cid(nc)
          do n = 1,clumps(cid)%ncells
             i = clumps(cid)%ixy(n)
             j = clumps(cid)%jxy(n)
             call get_gcell_info (i, j, wtxy, nlunits=ilunits, ncols=icols, npfts=ipfts)
             g = g + 1
             gcelldc(g)%ixy = i
             gcelldc(g)%jxy = j
             gcelldc(g)%lf = gcelldc(g)%li + ilunits - 1
             gcelldc(g)%cf = gcelldc(g)%ci + icols   - 1
             gcelldc(g)%pf = gcelldc(g)%pi + ipfts   - 1
             if (g <= numg-1) then
                gcelldc(g+1)%li = gcelldc(g)%lf + 1
                gcelldc(g+1)%ci = gcelldc(g)%cf + 1
                gcelldc(g+1)%pi = gcelldc(g)%pf + 1
             end if
          end do
       end do
    end do

    ! Find corresponding south->north gridcell for each decomposition gridcell
    ! and corresponding decomposition gridcell for each south->north gridcell

    do gdc = 1,numg
       i = gcelldc(gdc)%ixy
       j = gcelldc(gdc)%jxy
       do gsn = 1,numg
          if (gcellsn(gsn)%ixy == i .and. gcellsn(gsn)%jxy == j) then
             gcellsn(gsn)%gdc = gdc
             gcelldc(gdc)%gsn = gsn
          end if
       end do
    end do

    ! Write out clump and proc info

    if (masterproc) then
       do pid = 0,npes-1
          write(6,*)
          write(6,*)'proc= ',pid,' beg gridcell= ',procs(pid)%begg,' end gridcell= ',procs(pid)%endg, &
               ' total gridcells per proc= ',procs(pid)%ncells
          write(6,*)'proc= ',pid,' beg landunit= ',procs(pid)%begl,' end landunit= ',procs(pid)%endl, &
               ' total landunits per proc = ',procs(pid)%nlunits
          write(6,*)'proc= ',pid,' beg column  = ',procs(pid)%begc,' end column  = ',procs(pid)%endc, &
               ' total columns per proc  = ',procs(pid)%ncols
          write(6,*)'proc= ',pid,' beg pft     = ',procs(pid)%begp,' end pft     = ',procs(pid)%endp, &
               ' total pfts per proc     = ',procs(pid)%npfts
          do n = 1,procs(pid)%nclumps
             cid = procs(pid)%cid(n)
             write(6,*)'proc= ',pid,' clump no = ',n,' clump id= ',procs(pid)%cid(n), &
                  ' beg gridcell= ',clumps(cid)%begg,' end gridcell= ',clumps(cid)%endg, &
                  ' total gridcells per clump= ',clumps(cid)%ncells
             write(6,*)'proc= ',pid,' clump no = ',n,' clump id= ',procs(pid)%cid(n), &
                  ' beg landunit= ',clumps(cid)%begl,' end landunit= ',clumps(cid)%endl, &
                  ' total landunits per clump = ',clumps(cid)%nlunits
             write(6,*)'proc= ',pid,' clump no = ',n,' clump id= ',procs(pid)%cid(n), &
                  ' beg column  = ',clumps(cid)%begc,' end column  = ',clumps(cid)%endc, &
                  ' total columns per clump  = ',clumps(cid)%ncols
             write(6,*)'proc= ',pid,' clump no = ',n,' clump id= ',procs(pid)%cid(n), &
                  ' beg pft     = ',clumps(cid)%begp,' end pft     = ',clumps(cid)%endp, &
                  ' total pfts per clump     = ',clumps(cid)%npfts
          end do
       end do
       write(6,*)' total gridcells over all procs = ', procs(0)%numg
       write(6,*)' total landunts  over all procs = ', procs(0)%numl
       write(6,*)' total columns   over all procs = ', procs(0)%numc
       write(6,*)' total pfts      over all procs = ', procs(0)%nump
       write(6,*)
       call shr_sys_flush(6)
    end if

    deallocate(clumpfull)

  end subroutine initDecomp

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_nclumps
!
! !INTERFACE:
   integer function get_nclumps()
!
! !DESCRIPTION:
! This function returns the number of clumps on the land model grid.
!
! !ARGUMENTS:
     implicit none
!
! !CALLED FROM:
! subroutine lp_coupling_init() in module lp_coupling (lp_coupling.F90)
!
! !REVISION HISTORY:
! 2002.09.11  Forrest Hoffman  Creation.
!
!EOP
!------------------------------------------------------------------------------

     get_nclumps = nclumps

   end function get_nclumps

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_clump_cell_id_coord
!
! !INTERFACE:
   subroutine get_clump_cell_id_coord(lon, lat, cid, cell)
!
! !DESCRIPTION:
! This subroutine returns the id of the clump and cell corresponding to
! the longitude/latitude indices provided.
!
! !ARGUMENTS:
     implicit none
     integer, intent(in)  :: lon, lat      ! longitude/latitude indices
     integer, intent(out) :: cid, cell     ! clump and cell id
!
! !CALLED FROM:
! Unused.
!
! !REVISION HISTORY:
! 2002.09.11  Forrest Hoffman  Creation.
!
!EOP
!------------------------------------------------------------------------------

     cid  = pmulcs(lon,lat)%clumpid
     cell = pmulcs(lon,lat)%cell

   end subroutine get_clump_cell_id_coord

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_clump_owner_id
!
! !INTERFACE:
   integer function get_clump_owner_id(cid)
!
! !DESCRIPTION:
! This function returns the MPI process id (rank) responsible for the clump
! identified by the clump id cid.
!
! !ARGUMENTS:
     implicit none
     integer, intent(in) :: cid                    ! clump id
!
! !CALLED FROM:
! subroutine lp_coupling_init() in module lp_coupling (lp_coupling.F90)
!
! !REVISION HISTORY:
! 2002.09.11  Forrest Hoffman  Creation.
!
!EOP
!------------------------------------------------------------------------------

     get_clump_owner_id = clumps(cid)%owner

   end function get_clump_owner_id

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_clump_ncells_proc
!
! !INTERFACE:
   integer function get_clump_ncells_proc(pid)
!
! !DESCRIPTION:
! This function returns the number of cells contained within clumps owned
! by process pid.
!
! !ARGUMENTS:
     implicit none
     integer, intent(in) :: pid                     ! process id
!
! !LOCAL VARIABLES:
     integer :: cid                                 ! loop indices
     integer :: ncells                              ! cell counter
!
! !CALLED FROM:
! subroutine clm_map() (clm_map.F90)
!
! !REVISION HISTORY:
! 2002.09.11  Forrest Hoffman  Creation.
!
!EOP
!------------------------------------------------------------------------------

     ncells = 0
     do cid = 1,nclumps
        if (clumps(cid)%owner == pid) ncells = ncells + clumps(cid)%ncells
     end do
     get_clump_ncells_proc = ncells

   end function get_clump_ncells_proc

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_clump_ncells_id
!
! !INTERFACE:
   integer function get_clump_ncells_id(cid)
!
! !REVISION HISTORY:
! 2002.09.11  Forrest Hoffman  Creation.
!
! !ARGUMENTS:
     implicit none
     integer, intent(in) :: cid                    ! clump id
!
! !CALLED FROM:
! subroutine lp_coupling_init() in module lp_coupling (lp_coupling.F90)
!
! !DESCRIPTION:
! This function returns the number of cells contained within the clump
! identified by the clump id cid.
!
!EOP
!------------------------------------------------------------------------------

     get_clump_ncells_id = clumps(cid)%ncells

   end function get_clump_ncells_id

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_clump_coord_id
!
! !INTERFACE:
   subroutine get_clump_coord_id(cid, ncells, lons, lats)
!
! !DESCRIPTION:
! This subroutine returns the first ncells longitutde/latitude indices for
! the clump identified by the clump id cid.  In practice, this ncells is
! always the number of cells in clump cid: clumps(cid)%ncells.
!
! !ARGUMENTS:
     implicit none
     integer, intent(in)  :: cid                    ! clump id
     integer, intent(in)  :: ncells                 ! number of grid cells
     integer, intent(out) :: lons(ncells)           ! longitude indices
     integer, intent(out) :: lats(ncells)           ! latitude indices
!
! !LOCAL VARIABLES:
     integer :: i                                   ! loop index
!
! !CALLED FROM:
! subroutine lp_coupling_init() in module lp_coupling (lp_coupling.F90)
!
! !REVISION HISTORY:
! 2002.09.11  Forrest Hoffman  Creation.
!
!EOP
!------------------------------------------------------------------------------

     do i = 1,ncells
        lons(i) = clumps(cid)%ixy(i)
        lats(i) = clumps(cid)%jxy(i)
     end do

   end subroutine get_clump_coord_id

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_clump_gcell_info
!
! !INTERFACE:
   subroutine get_clump_gcell_info(cid, cell, gi)
!
! !DESCRIPTION:
! Determine beginning grid cell index
!
! !ARGUMENTS:
     implicit none
     integer, intent(in)  :: cid                    ! clump id
     integer, intent(in)  :: cell                   ! clump cell id
     integer, intent(out) :: gi                     ! 1d gridcell index
!
! !CALLED FROM:
! subroutines alltoall_clump_to_chunk_init(), alltoall_clump_to_chunk(), and
! alltoall_chunk_to_clump() in module lp_coupling (lp_coupling.F90)
!
! !REVISION HISTORY:
! 2002.11.17  Mariana Vertenstein  Creation.
!
!EOP
!------------------------------------------------------------------------------

     gi = clumps(cid)%gi(cell)

   end subroutine get_clump_gcell_info

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_clump_bounds
!
! !INTERFACE:
   subroutine get_clump_bounds (n, begg, endg, begl, endl, begc, endc, &
                                begp, endp)
!
! !ARGUMENTS:
     implicit none
     integer, intent(in)  :: n           ! proc clump index
     integer, intent(out) :: begp, endp  ! clump beginning and ending
                                         ! pft indices
     integer, intent(out) :: begc, endc  ! clump beginning and ending
                                         ! column indices
     integer, intent(out) :: begl, endl  ! clump beginning and ending
                                         ! landunit indices
     integer, intent(out) :: begg, endg  ! clump beginning and ending
                                         ! gridcell indices
!
! !DESCRIPTION:
! Determine clump beginning and ending pft, column, landunit and
! gridcell indices.
!
! !REVISION HISTORY:
! 2003.09.12  Mariana Vertenstein  Creation.
!
!EOP
!
! !LOCAL VARIABLES:
     integer :: cid                        ! clump id
!------------------------------------------------------------------------------

     cid  = procs(iam)%cid(n)
     begp = clumps(cid)%begp
     endp = clumps(cid)%endp
     begc = clumps(cid)%begc
     endc = clumps(cid)%endc
     begl = clumps(cid)%begl
     endl = clumps(cid)%endl
     begg = clumps(cid)%begg
     endg = clumps(cid)%endg

   end subroutine get_clump_bounds

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_proc_bounds
!
! !INTERFACE:
   subroutine get_proc_bounds (begg, endg, begl, endl, begc, endc, &
                               begp, endp)
!
! !ARGUMENTS:
     implicit none
     integer, intent(out) :: begp, endp  ! proc beginning and ending
                                         ! pft indices
     integer, intent(out) :: begc, endc  ! proc beginning and ending
                                         ! column indices
     integer, intent(out) :: begl, endl  ! proc beginning and ending
                                         ! landunit indices
     integer, intent(out) :: begg, endg  ! proc beginning and ending
                                         ! gridcell indices
! !DESCRIPTION:
! Retrieve gridcell, landunit, column, and pft bounds for process.
!
! !REVISION HISTORY:
! 2003.09.12  Mariana Vertenstein  Creation.
!
!EOP
!------------------------------------------------------------------------------

     begp = procs(iam)%begp
     endp = procs(iam)%endp
     begc = procs(iam)%begc
     endc = procs(iam)%endc
     begl = procs(iam)%begl
     endl = procs(iam)%endl
     begg = procs(iam)%begg
     endg = procs(iam)%endg

   end subroutine get_proc_bounds

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_proc_total
!
! !INTERFACE:
   subroutine get_proc_total(pid, ncells, nlunits, ncols, npfts)
!
! !DESCRIPTION:
! Count up gridcells, landunits, columns, and pfts on process.
!
! !ARGUMENTS:
     implicit none
     integer, intent(in)  :: pid     ! proc id
     integer, intent(out) :: ncells  ! total number of gridcells
                                     ! on the processor
     integer, intent(out) :: nlunits ! total number of landunits
                                     ! on the processor
     integer, intent(out) :: ncols   ! total number of columns
                                     ! on the processor
     integer, intent(out) :: npfts   ! total number of pfts
                                     ! on the processor
!
! !REVISION HISTORY:
! 2003.09.12  Mariana Vertenstein  Creation.
!
!EOP
!
! !LOCAL VARIABLES:
   integer :: cid       ! clump index
!------------------------------------------------------------------------------

     npfts   = 0
     nlunits = 0
     ncols   = 0
     ncells  = 0
     do cid = 1,nclumps
        if (clumps(cid)%owner == pid) then
           ncells  = ncells  + clumps(cid)%ncells
           nlunits = nlunits + clumps(cid)%nlunits
           ncols   = ncols   + clumps(cid)%ncols
           npfts   = npfts   + clumps(cid)%npfts
        end if
     end do

   end subroutine get_proc_total

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_proc_global
!
! !INTERFACE:
   subroutine get_proc_global(numg, numl, numc, nump)
!
! !DESCRIPTION:
! Return number of gridcells, landunits, columns, and pfts across all
! processes.
!
! !ARGUMENTS:
     implicit none
     integer, intent(out) :: numg  ! total number of gridcells
                                   ! across all processors
     integer, intent(out) :: numl  ! total number of landunits
                                   ! across all processors
     integer, intent(out) :: numc  ! total number of columns
                                   ! across all processors
     integer, intent(out) :: nump  ! total number of pfts
                                   ! across all processors
! !REVISION HISTORY:
! 2003.09.12  Mariana Vertenstein  Creation.
!
!EOP
!------------------------------------------------------------------------------

     nump = procs(iam)%nump
     numc = procs(iam)%numc
     numl = procs(iam)%numl
     numg = procs(iam)%numg

   end subroutine get_proc_global

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_proc_clumps
!
! !INTERFACE:
   integer function get_proc_clumps()
!
! !DESCRIPTION:
! Return the number of clumps.
!
! !ARGUMENTS:
     implicit none
!
! !RETURN VALUE:
!    integer :: get_proc_clumps
!
! !REVISION HISTORY:
! 2003.09.12  Mariana Vertenstein  Creation.
!
!EOP
!------------------------------------------------------------------------------

     get_proc_clumps = procs(iam)%nclumps

   end function get_proc_clumps

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: map_dc2sn_sl_real
!
! !INTERFACE:
   subroutine map_dc2sn_sl_real(arraydc, arraysn, type1d)
!
! !DESCRIPTION:
!  Maps a decomposition single level real array into a
!  south->north single level real array
!
! !USES:
    use clmtype, only : nameg, namel, namec, namep
!
! !ARGUMENTS:
     implicit none
     real(r8), pointer :: arraydc(:)
     real(r8), pointer :: arraysn(:)
     character(len=*), intent(in) :: type1d
!
! !REVISION HISTORY:
! 2003.12.01  Mariana Vertenstein  Creation.
!
!EOP
!
! !LOCAL VARIABLES:
     integer :: dc,sn,gdc,gsn
     integer :: dci,dcf,sni,snf
     integer :: n,l,c,p
     logical :: error
!------------------------------------------------------------------------------

     error = .false.
     do gdc = 1,procs(iam)%numg
        gsn = gcelldc(gdc)%gsn
        if (type1d == nameg) then
           arraysn(gsn) = arraydc(gdc)
        else if (type1d == namel) then
           dci = gcelldc(gdc)%li
           dcf = gcelldc(gdc)%lf
           sni = gcellsn(gsn)%li
           snf = gcellsn(gsn)%lf
           n = dcf - dci + 1
           if (n /= snf - sni + 1) then
              error = .true.
              exit
           end if
           do l = 1,n
              sn = gcellsn(gsn)%li + l - 1
              dc = gcelldc(gdc)%li + l - 1
              arraysn(sn) = arraydc(dc)
           end do
        else if (type1d == namec) then
           dci = gcelldc(gdc)%ci
           dcf = gcelldc(gdc)%cf
           sni = gcellsn(gsn)%ci
           snf = gcellsn(gsn)%cf
           n = dcf - dci + 1
           if (n /= snf - sni + 1) then
              error = .true.
              exit
           end if
           do c = 1,n
              sn = gcellsn(gsn)%ci + c - 1
              dc = gcelldc(gdc)%ci + c - 1
              arraysn(sn) = arraydc(dc)
           end do
        else if (type1d == namep) then
           dci = gcelldc(gdc)%pi
           dcf = gcelldc(gdc)%pf
           sni = gcellsn(gsn)%pi
           snf = gcellsn(gsn)%pf
           n = dcf - dci + 1
           if (n /= snf - sni + 1) then
              error = .true.
              exit
           end if
           do p = 1,n
              sn = gcellsn(gsn)%pi + p - 1
              dc = gcelldc(gdc)%pi + p - 1
              arraysn(sn) = arraydc(dc)
           end do
        end if
     end do
     if (error)  then
        write(6,*)'error in  map_dc2n_sl_real '
        call endrun()
     end if

   end subroutine map_dc2sn_sl_real

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: map_dc2sn_sl_int
!
! !INTERFACE:
   subroutine map_dc2sn_sl_int(arraydc, arraysn, type1d)
!
! !DESCRIPTION:
!  Maps a decomposition single level integer array into a
!  south->north single level real array
!
! !USES:
    use clmtype, only : nameg, namel, namec, namep
!
! !ARGUMENTS:
     implicit none
     integer, pointer :: arraydc(:)
     integer, pointer :: arraysn(:)
     character(len=*), intent(in) :: type1d
!
! !REVISION HISTORY:
! 2003.12.01  Mariana Vertenstein  Creation.
!
!EOP
!
! !LOCAL VARIABLES:
     integer :: dc,sn,gdc,gsn
     integer :: dci,dcf,sni,snf
     integer :: n,l,c,p
     logical :: error
!------------------------------------------------------------------------------

     error = .false.
     do gdc = 1,procs(iam)%numg
        gsn = gcelldc(gdc)%gsn
        if (type1d == nameg) then
           arraysn(gsn) = arraydc(gdc)
        else if (type1d == namel) then
           dci = gcelldc(gdc)%li
           dcf = gcelldc(gdc)%lf
           sni = gcellsn(gsn)%li
           snf = gcellsn(gsn)%lf
           n = dcf - dci + 1
           if (n /= snf - sni + 1) then
              error = .true.
              exit
           end if
           do l = 1,n
              sn = gcellsn(gsn)%li + l - 1
              dc = gcelldc(gdc)%li + l - 1
              arraysn(sn) = arraydc(dc)
           end do
        else if (type1d == namec) then
           dci = gcelldc(gdc)%ci
           dcf = gcelldc(gdc)%cf
           sni = gcellsn(gsn)%ci
           snf = gcellsn(gsn)%cf
           n = dcf - dci + 1
           if (n /= snf - sni + 1) then
              error = .true.
              exit
           end if
           do c = 1,n
              sn = gcellsn(gsn)%ci + c - 1
              dc = gcelldc(gdc)%ci + c - 1
              arraysn(sn) = arraydc(dc)
           end do
        else if (type1d == namep) then
           dci = gcelldc(gdc)%pi
           dcf = gcelldc(gdc)%pf
           sni = gcellsn(gsn)%pi
           snf = gcellsn(gsn)%pf
           n = dcf - dci + 1
           if (n /= snf - sni + 1) then
              error = .true.
              exit
           end if
           do p = 1,n
              sn = gcellsn(gsn)%pi + p - 1
              dc = gcelldc(gdc)%pi + p - 1
              arraysn(sn) = arraydc(dc)
           end do
        end if
     end do
     if (error)  then
        write(6,*)'error in  map_dc2n_sl '
        call endrun()
     end if

   end subroutine map_dc2sn_sl_int

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: map_dc2sn_ml1_real
!
! !INTERFACE:
   subroutine map_dc2sn_ml1_real(arraydc, arraysn, type1d, lb1, ub1, revord)
!
! !DESCRIPTION:
!  Maps a decomposition (dc) multilevel real array into a
!  south->north (sn) multilevel level real array
!
! !USES:
    use clmtype, only : nameg, namel, namec, namep
!
! !ARGUMENTS:
     implicit none
     real(r8), pointer :: arraydc(:,:)
     real(r8), pointer :: arraysn(:,:)
     character(len=*), intent(in) :: type1d
     integer, intent(in) :: lb1
     integer, intent(in) :: ub1
     logical, intent(in), optional :: revord
!
! !REVISION HISTORY:
! 2003.12.01  Mariana Vertenstein  Creation.
!
!EOP
!
! !LOCAL VARIABLES:
     integer :: dc,sn,gdc,gsn,lev
     integer :: dci,dcf,sni,snf
     integer :: n,l,c,p
     logical :: error
!------------------------------------------------------------------------------

     error = .false.
     do gdc = 1,procs(iam)%numg
        gsn = gcelldc(gdc)%gsn
        if (type1d == nameg) then
           do lev = lb1, ub1
              sn = gsn
              dc = gdc
              if (present(revord)) then
                 if (revord) then
                    arraysn(sn,lev) = arraydc(dc,lev)
                 end if
              else
                 arraysn(lev,sn) = arraydc(lev,dc)
              end if
           end do
        else if (type1d == namel) then
           dci = gcelldc(gdc)%li
           dcf = gcelldc(gdc)%lf
           sni = gcellsn(gsn)%li
           snf = gcellsn(gsn)%lf
           n = dcf - dci + 1
           if (n /= snf - sni + 1) then
              error = .true.
              exit
           end if
           do l = 1,n
              sn = gcellsn(gsn)%li + l - 1
              dc = gcelldc(gdc)%li + l - 1
              do lev = lb1, ub1
                 if (present(revord)) then
                    if (revord) then
                       arraysn(sn,lev) = arraydc(dc,lev)
                    end if
                 else
                    arraysn(lev,sn) = arraydc(lev,dc)
                 end if
              end do
           end do
        else if (type1d == namec) then
           dci = gcelldc(gdc)%ci
           dcf = gcelldc(gdc)%cf
           sni = gcellsn(gsn)%ci
           snf = gcellsn(gsn)%cf
           n = dcf - dci + 1
           if (n /= snf - sni + 1) then
              error = .true.
              exit
           end if
           do c = 1,n
              sn = gcellsn(gsn)%ci + c - 1
              dc = gcelldc(gdc)%ci + c - 1
              do lev = lb1, ub1
                 if (present(revord)) then
                    if (revord) then
                       arraysn(sn,lev) = arraydc(dc,lev)
                    end if
                 else
                    arraysn(lev,sn) = arraydc(lev,dc)
                 end if
              end do
           end do
        else if (type1d == namep) then
           dci = gcelldc(gdc)%pi
           dcf = gcelldc(gdc)%pf
           sni = gcellsn(gsn)%pi
           snf = gcellsn(gsn)%pf
           n = dcf - dci + 1
           if (n /= snf - sni + 1) then
              error = .true.
              exit
           end if
           do p = 1,n
              sn = gcellsn(gsn)%pi + p - 1
              dc = gcelldc(gdc)%pi + p - 1
              do lev = lb1, ub1
                 if (present(revord)) then
                    if (revord) then
                       arraysn(sn,lev) = arraydc(dc,lev)
                    end if
                 else
                    arraysn(lev,sn) = arraydc(lev,dc)
                 end if
              end do
           end do
        end if
     end do
     if (error)  then
        write(6,*)'ndc = ',n,' nsn= ',snf - sni + 1
        write(6,*)'snf= ',snf,' sni= ',sni
        write(6,*)'dcf= ',dcf,' dci= ',dci
        write(6,*)'error in  map_dc2sn_ml1 '
        call endrun()
     end if

   end subroutine map_dc2sn_ml1_real

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: map_dc2sn_ml1_int
!
! !INTERFACE:
   subroutine map_dc2sn_ml1_int(arraydc, arraysn, type1d, lb1, ub1, revord)
!
! !DESCRIPTION:
!  Maps a decomposition (dc) multilevel integer array into a
!  south->north (sn) multilevel level integer array
!
! !USES:
    use clmtype, only : nameg, namel, namec, namep
!
! !ARGUMENTS:
     implicit none
     integer, pointer :: arraydc(:,:)
     integer, pointer :: arraysn(:,:)
     character(len=*), intent(in) :: type1d
     integer, intent(in) :: lb1
     integer, intent(in) :: ub1
     logical, intent(in), optional :: revord
!
! !REVISION HISTORY:
! 2003.12.01  Mariana Vertenstein  Creation.
!
!EOP
!
! !LOCAL VARIABLES:
     integer :: dc,sn,gdc,gsn,lev
     integer :: dci,dcf,sni,snf
     integer :: n,l,c,p
     logical :: error
!------------------------------------------------------------------------------

     error = .false.
     do gdc = 1,procs(iam)%numg
        gsn = gcelldc(gdc)%gsn
        if (type1d == nameg) then
           sn = gsn
           dc = gdc
           do lev = lb1, ub1
              if (present(revord)) then
                 if (revord) then
                    arraysn(sn,lev) = arraydc(dc,lev)
                 end if
              else
                 arraysn(lev,sn) = arraydc(lev,dc)
              end if
           end do
        else if (type1d == namel) then
           dci = gcelldc(gdc)%li
           dcf = gcelldc(gdc)%lf
           sni = gcellsn(gsn)%li
           snf = gcellsn(gsn)%lf
           n = dcf - dci + 1
           if (n /= snf - sni + 1) then
              error = .true.
              exit
           end if
           do l = 1,n
              sn = gcellsn(gsn)%li + l - 1
              dc = gcelldc(gdc)%li + l - 1
              do lev = lb1, ub1
                 if (present(revord)) then
                    if (revord) then
                       arraysn(sn,lev) = arraydc(dc,lev)
                    end if
                 else
                    arraysn(lev,sn) = arraydc(lev,dc)
                 end if
              end do
           end do
        else if (type1d == namec) then
           dci = gcelldc(gdc)%ci
           dcf = gcelldc(gdc)%cf
           sni = gcellsn(gsn)%ci
           snf = gcellsn(gsn)%cf
           n = dcf - dci + 1
           if (n /= snf - sni + 1) then
              error = .true.
              exit
           end if
           do c = 1,n
              sn = gcellsn(gsn)%ci + c - 1
              dc = gcelldc(gdc)%ci + c - 1
              do lev = lb1,ub1
                 if (present(revord)) then
                    if (revord) then
                       arraysn(sn,lev) = arraydc(dc,lev)
                    end if
                 else
                    arraysn(lev,sn) = arraydc(lev,dc)
                 end if
              end do
           end do
        else if (type1d == namep) then
           dci = gcelldc(gdc)%pi
           dcf = gcelldc(gdc)%pf
           sni = gcellsn(gsn)%pi
           snf = gcellsn(gsn)%pf
           n = dcf - dci + 1
           if (n /= snf - sni + 1) then
              error = .true.
              exit
           end if
           do p = 1,n
              sn = gcellsn(gsn)%pi + p - 1
              dc = gcelldc(gdc)%pi + p - 1
              do lev = lb1,ub1
                 if (present(revord)) then
                    if (revord) then
                       arraysn(sn,lev) = arraydc(dc,lev)
                    end if
                 else
                    arraysn(lev,sn) = arraydc(lev,dc)
                 end if
              end do
           end do
        end if
     end do
     if (error)  then
        write(6,*)'ndc = ',n,' nsn= ',snf - sni + 1
        write(6,*)'snf= ',snf,' sni= ',sni
        write(6,*)'dcf= ',dcf,' dci= ',dci
        write(6,*)'error in map_dc2sn_ml1 '
        call endrun()
     end if

   end subroutine map_dc2sn_ml1_int

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: map_sn2dc_sl_real
!
! !INTERFACE:
   subroutine map_sn2dc_sl_real(arraysn, arraydc, type1d)
!
! !DESCRIPTION:
!  Maps a south->north single level real array into a
!  decomposition single level real array
!
! !USES:
    use clmtype, only : nameg, namel, namec, namep
!
! !ARGUMENTS:
     implicit none
     real(r8), pointer :: arraysn(:)
     real(r8), pointer :: arraydc(:)
     character(len=*), intent(in) :: type1d
!
! !REVISION HISTORY:
! 2003.12.01  Mariana Vertenstein  Creation.
!
!EOP
!
! !LOCAL VARIABLES:
     integer :: dc,sn,gdc,gsn
     integer :: dci,dcf,sni,snf
     integer :: n,l,c,p
     logical :: error
!------------------------------------------------------------------------------

     error = .false.
     do gsn = 1,procs(iam)%numg
        gdc = gcellsn(gsn)%gdc
        if (type1d == nameg) then
           sn = gsn
           dc = gdc
           arraydc(dc) = arraysn(sn)
        else if (type1d == namel) then
           sni = gcellsn(gsn)%li
           snf = gcellsn(gsn)%lf
           dci = gcelldc(gdc)%li
           dcf = gcelldc(gdc)%lf
           n = snf - sni + 1
           if (n /= dcf - dci + 1) then
              error = .true.
              exit
           end if
           do l = 1,n
              dc = gcelldc(gdc)%li + l - 1
              sn = gcellsn(gsn)%li + l - 1
              arraydc(dc) = arraysn(sn)
           end do
        else if (type1d == namec) then
           sni = gcellsn(gsn)%ci
           snf = gcellsn(gsn)%cf
           dci = gcelldc(gdc)%ci
           dcf = gcelldc(gdc)%cf
           n = snf - sni + 1
           if (n /= dcf - dci + 1) then
              error = .true.
              exit
           end if
           do c = 1,n
              dc = gcelldc(gdc)%ci + c - 1
              sn = gcellsn(gsn)%ci + c - 1
              arraydc(dc) = arraysn(sn)
           end do
        else if (type1d == namep) then
           sni = gcellsn(gsn)%pi
           snf = gcellsn(gsn)%pf
           dci = gcelldc(gdc)%pi
           dcf = gcelldc(gdc)%pf
           n = snf - sni + 1
           if (n /= dcf - dci + 1) then
              error = .true.
              exit
           end if
           do p = 1,n
              dc = gcelldc(gdc)%pi + p - 1
              sn = gcellsn(gsn)%pi + p - 1
              arraydc(dc) = arraysn(sn)
           end do
        end if
     end do
     if (error)  then
        write(6,*)'error in  map_sn2dc_sl '
        call endrun()
     end if

   end subroutine map_sn2dc_sl_real

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: map_sn2dc_sl_int
!
! !INTERFACE:
   subroutine map_sn2dc_sl_int(arraysn, arraydc, type1d)
!
! !DESCRIPTION:
!  Maps a south->north single level integer array into a
!  decomposition single level integer array
!
! !USES:
    use clmtype, only : nameg, namel, namec, namep
!
! !ARGUMENTS:
     implicit none
     integer, pointer :: arraysn(:)
     integer, pointer :: arraydc(:)
     character(len=*), intent(in) :: type1d
!
! !REVISION HISTORY:
! 2003.12.01  Mariana Vertenstein  Creation.
!
!EOP
!
! !LOCAL VARIABLES:
     integer :: dc,sn,gdc,gsn
     integer :: dci,dcf,sni,snf
     integer :: n,l,c,p
     logical :: error
!------------------------------------------------------------------------------

     error = .false.
     do gsn = 1,procs(iam)%numg
        gdc = gcellsn(gsn)%gdc
        if (type1d == nameg) then
           sn = gsn
           dc = gdc
           arraydc(dc) = arraysn(sn)
        else if (type1d == namel) then
           sni = gcellsn(gsn)%li
           snf = gcellsn(gsn)%lf
           dci = gcelldc(gdc)%li
           dcf = gcelldc(gdc)%lf
           n = snf - sni + 1
           if (n /= dcf - dci + 1) then
              error = .true.
              exit
           end if
           do l = 1,n
              dc = gcelldc(gdc)%li + l - 1
              sn = gcellsn(gsn)%li + l - 1
              arraydc(dc) = arraysn(sn)
           end do
        else if (type1d == namec) then
           sni = gcellsn(gsn)%ci
           snf = gcellsn(gsn)%cf
           dci = gcelldc(gdc)%ci
           dcf = gcelldc(gdc)%cf
           n = snf - sni + 1
           if (n /= dcf - dci + 1) then
              error = .true.
              exit
           end if
           do c = 1,n
              dc = gcelldc(gdc)%ci + c - 1
              sn = gcellsn(gsn)%ci + c - 1
              arraydc(dc) = arraysn(sn)
           end do
        else if (type1d == namep) then
           sni = gcellsn(gsn)%pi
           snf = gcellsn(gsn)%pf
           dci = gcelldc(gdc)%pi
           dcf = gcelldc(gdc)%pf
           n = snf - sni + 1
           if (n /= dcf - dci + 1) then
              error = .true.
              exit
           end if
           do p = 1,n
              dc = gcelldc(gdc)%pi + p - 1
              sn = gcellsn(gsn)%pi + p - 1
              arraydc(dc) = arraysn(sn)
           end do
        end if
     end do
     if (error)  then
        write(6,*)'error in  map_sn2dc_sl '
        call endrun()
     end if

   end subroutine map_sn2dc_sl_int

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: map_sn2dc_ml1_real
!
! !INTERFACE:
   subroutine map_sn2dc_ml1_real(arraysn, arraydc, type1d, lb1, ub1)
!
! !DESCRIPTION:
!  Maps a south->north multi level real array into a
!  decomposition multi level real array
!
! !USES:
    use clmtype, only : nameg, namel, namec, namep
!
! !ARGUMENTS:
     implicit none
     real(r8), pointer :: arraysn(:,:)
     real(r8), pointer :: arraydc(:,:)
     character(len=*), intent(in) :: type1d
     integer, intent(in)  :: lb1
     integer, intent(in)  :: ub1
!
! !REVISION HISTORY:
! 2003.12.01  Mariana Vertenstein  Creation.
!
!EOP
!
! !LOCAL VARIABLES:
     integer :: dc,sn,gdc,gsn
     integer :: dci,dcf,sni,snf
     integer :: n,l,c,p,lev
     logical :: error

!------------------------------------------------------------------------------

     error = .false.
     do gsn = 1,procs(iam)%numg
        gdc = gcellsn(gsn)%gdc
        if (type1d == nameg) then
           sn = gsn
           dc = gdc
           do lev = lb1, ub1
              arraydc(lev,dc) = arraysn(lev,sn)
           end do
        else if (type1d == namel) then
           sni = gcellsn(gsn)%li
           snf = gcellsn(gsn)%lf
           dci = gcelldc(gdc)%li
           dcf = gcelldc(gdc)%lf
           n = snf - sni + 1
           if (n /= dcf - dci + 1) then
              error = .true.
              exit
           end if
           do l = 1,n
              dc = gcelldc(gdc)%li + l - 1
              sn = gcellsn(gsn)%li + l - 1
              do lev = lb1, ub1
                 arraydc(lev,dc) = arraysn(lev,sn)
              end do
           end do
        else if (type1d == namec) then
           sni = gcellsn(gsn)%ci
           snf = gcellsn(gsn)%cf
           dci = gcelldc(gdc)%ci
           dcf = gcelldc(gdc)%cf
           n = snf - sni + 1
           if (n /= dcf - dci + 1) then
              error = .true.
              exit
           end if
           do c = 1,n
              dc = gcelldc(gdc)%ci + c - 1
              sn = gcellsn(gsn)%ci + c - 1
              do lev = lb1,ub1
                 arraydc(lev,dc) = arraysn(lev,sn)
              end do
           end do
        else if (type1d == namep) then
           sni = gcellsn(gsn)%pi
           snf = gcellsn(gsn)%pf
           dci = gcelldc(gdc)%pi
           dcf = gcelldc(gdc)%pf
           n = snf - sni + 1
           if (n /= dcf - dci + 1) then
              error = .true.
              exit
           end if
           do p = 1,n
              dc = gcelldc(gdc)%pi + p - 1
              sn = gcellsn(gsn)%pi + p - 1
              do lev = lb1,ub1
                 arraydc(lev,dc) = arraysn(lev,sn)
              end do
           end do
        end if
     end do
     if (error)  then
        write(6,*)'error in  map_sn2n_ml1 '
        call endrun()
     end if

   end subroutine map_sn2dc_ml1_real

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: map_sn2dc_ml1_int
!
! !INTERFACE:
   subroutine map_sn2dc_ml1_int(arraysn, arraydc, type1d, lb1, ub1)
!
! !DESCRIPTION:
!  Maps a south->north multi level integer array into a
!  decomposition multi level integer array
!
! !USES:
    use clmtype, only : nameg, namel, namec, namep
!
! !ARGUMENTS:
     implicit none
     integer, pointer :: arraysn(:,:)
     integer, pointer :: arraydc(:,:)
     character(len=*), intent(in) :: type1d
     integer, intent(in)  :: lb1
     integer, intent(in)  :: ub1
!
! !REVISION HISTORY:
! 2003.12.01  Mariana Vertenstein  Creation.
!
!EOP
!
! !LOCAL VARIABLES:
     integer :: dc,sn,gdc,gsn
     integer :: dci,dcf,sni,snf
     integer :: n,l,c,p,lev
     logical :: error
!------------------------------------------------------------------------------

     error = .false.
     do gsn = 1,procs(iam)%numg
        gdc = gcellsn(gsn)%gdc
        if (type1d == nameg) then
           sn = gsn
           dc = gdc
           do lev = lb1, ub1
              arraydc(lev,dc) = arraysn(lev,sn)
           end do
        else if (type1d == namel) then
           sni = gcellsn(gsn)%li
           snf = gcellsn(gsn)%lf
           dci = gcelldc(gdc)%li
           dcf = gcelldc(gdc)%lf
           n = snf - sni + 1
           if (n /= dcf - dci + 1) then
              error = .true.
              exit
           end if
           do l = 1,n
              dc = gcelldc(gdc)%li + l - 1
              sn = gcellsn(gsn)%li + l - 1
              do lev = lb1, ub1
                 arraydc(lev,dc) = arraysn(lev,sn)
              end do
           end do
        else if (type1d == namec) then
           sni = gcellsn(gsn)%ci
           snf = gcellsn(gsn)%cf
           dci = gcelldc(gdc)%ci
           dcf = gcelldc(gdc)%cf
           n = snf - sni + 1
           if (n /= dcf - dci + 1) then
              error = .true.
              exit
           end if
           do c = 1,n
              dc = gcelldc(gdc)%ci + c - 1
              sn = gcellsn(gsn)%ci + c - 1
              do lev = lb1,ub1
                 arraydc(lev,dc) = arraysn(lev,sn)
              end do
           end do
        else if (type1d == namep) then
           sni = gcellsn(gsn)%pi
           snf = gcellsn(gsn)%pf
           dci = gcelldc(gdc)%pi
           dcf = gcelldc(gdc)%pf
           n = snf - sni + 1
           if (n /= dcf - dci + 1) then
              error = .true.
              exit
           end if
           do p = 1,n
              dc = gcelldc(gdc)%pi + p - 1
              sn = gcellsn(gsn)%pi + p - 1
              do lev = lb1,ub1
                 arraydc(lev,dc) = arraysn(lev,sn)
              end do
           end do
        end if
     end do
     if (error)  then
        write(6,*)'error in  map_sn2n_ml1 '
        call endrun()
     end if

   end subroutine map_sn2dc_ml1_int

end module decompMod
