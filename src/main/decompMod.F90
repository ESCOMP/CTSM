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
!                                 ! into clumps and processors
  public get_nclumps             ! returns the number of clumps defined
  public get_clump_owner_id      ! returns clump owner based on clump id
  public get_clump_ncells_proc   ! returns number of cells for process
  public get_clump_ncells_id     ! returns number of cells in clump
  public get_clump_gcell_info    ! returns 1d gridcell index
  public get_clump_bounds        ! beg and end gridcell, landunit, column,
                                 ! pft indices for clump
  public get_clump_coord_id      ! i,j indices for clump gridcells
  public get_proc_clumps         ! number of clumps for this processor
  public get_proc_bounds_atm     ! beg and end gridcell for atm
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
! 2005.11.01  T Craig  Rewrite
!
!EOP
!
! !PRIVATE TYPES:
  private

  integer :: nclumps     ! total number of clumps across all processors

  !---global information on each pe
  type processor_type
     integer :: nclumps          ! number of clumps for processor_type iam
     integer,pointer :: cid(:)   ! clump indices
     integer :: ncells           ! number of gridcells in proc
     integer :: nlunits          ! number of landunits in proc
     integer :: ncols            ! number of columns in proc
     integer :: npfts            ! number of pfts in proc
     integer :: begg, endg       ! beginning and ending gridcell index
     integer :: begl, endl       ! beginning and ending landunit index
     integer :: begc, endc       ! beginning and ending column index
     integer :: begp, endp       ! beginning and ending pft index
     integer :: numg             ! total number of gridcells on all procs
     integer :: numl             ! total number of landunits on all procs
     integer :: numc             ! total number of columns on all procs
     integer :: nump             ! total number of pfts on all procs
  end type processor_type
  type(processor_type), allocatable :: procs(:)

  !---global information on each pe
  type clump_type
     integer :: owner            ! process id owning clump
     integer :: ncells           ! number of gridcells in clump
     integer :: nlunits          ! number of landunits in clump
     integer :: ncols            ! number of columns in clump
     integer :: npfts            ! number of pfts in clump
     integer :: begg, endg       ! beginning and ending gridcell index
     integer :: begl, endl       ! beginning and ending landunit index
     integer :: begc, endc       ! beginning and ending column index
     integer :: begp, endp       ! beginning and ending pft index
  end type clump_type
  type(clump_type),  allocatable :: clumps(:)

  !---global information on each pe
  type decomp_type
     integer,pointer :: gsn2gdc(:)    ! 1d gsn to 1d gdc
     integer,pointer :: gdc2gsn(:)    ! 1d gdc to 1d gsn
     integer,pointer :: gsn2i(:)      ! 1d gsn to 2d sn i index
     integer,pointer :: gsn2j(:)      ! 1d gsn to 2d sn j index
     integer,pointer :: gdc2i(:)      ! 1d gdc to 2d sn i index
     integer,pointer :: gdc2j(:)      ! 1d gdc to 2d sn j index
     integer,pointer :: ij2gsn(:,:)   ! 2d sn i,j index to 1d gsn
     integer,pointer :: ij2gdc(:,:)   ! 2d sn i,j index to 1d gdc
  end type decomp_type
  public decomp_type
  type(decomp_type),public,target :: ldecomp
  type(decomp_type),public,target :: adecomp

  integer,pointer,save,public     :: abegg(:)
  integer,pointer,save,public     :: aendg(:)

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
  subroutine initDecomp(wtxy,n_ovr,i_ovr,j_ovr)
!
! !DESCRIPTION:
! This subroutine initializes the land surface decomposition into a clump
! data structure.
!
! !USES:
    use domainMod     , only : ldomain,adomain
    use initSubgridMod, only : get_gcell_info
!
! !ARGUMENTS:
    implicit none
    real(r8), intent(in) :: wtxy(:,:,:)   ! subgrid patch weights
    integer , pointer :: n_ovr(:,:)    ! num of cells for each dst
    integer , pointer :: i_ovr(:,:,:)  ! i index, map input cell
    integer , pointer :: j_ovr(:,:,:)  ! j index, map input cell

!
! !LOCAL VARIABLES:
    integer :: ppc                    ! min number of pfts per clump
    integer :: lpc                    ! min number of landunits per clump
    integer :: ppclump                ! min pfts per clump
    integer :: lni,lnj                ! land domain global size
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
    logical :: error                  ! temporary for finding full clump
    integer :: clumpcount             ! temporary for finding full clump
    integer, pointer :: cidi(:,:)     ! temporary for setting decomp
    integer, pointer :: cidj(:,:)     ! temporary for setting decomp
    logical, pointer :: lcell(:,:)    ! temporary for tracking used land cells
    integer :: ilunits, icols, ipfts  ! temporaries
    integer :: ng                     ! temporaries
    integer :: nl                     ! temporaries
    integer :: nc                     ! temporaries
    integer :: np                     ! temporaries
    integer :: ier                    ! error code

    integer :: ani,anj                ! atm domain global size
    integer :: ai,aj,ag               ! indices
    integer :: ancells                ! atm ncells
    integer, pointer :: acidi(:,:)    ! temporary for setting adecomp
    integer, pointer :: acidj(:,:)    ! temporary for setting adecomp
    integer, pointer :: clumps_ancells(:) ! tmp for number of atm cells/clump
!
! !CALLED FROM:
! subroutine initialize
!
! !REVISION HISTORY:
! 2002.09.11  Forrest Hoffman  Creation.
! 2005.12.15  T Craig Updated for finemesh
!
!EOP
!------------------------------------------------------------------------------

    lni = ldomain%ni
    lnj = ldomain%nj
    ani = adomain%ni
    anj = adomain%nj

    ier = 0
    allocate(lcell(lni,lnj),stat=ier)
    if (ier /= 0) then
       write (6,*) 'initDecomp(): allocation error for lcell'
       call endrun
    end if
    lcell = .false.

    ! Dynamic memory allocation for procs

    if( .not. allocated(procs)) allocate(procs(0:npes-1), stat=ier)
    if (ier /= 0) then
       write (6,*) 'initDecomp(): allocation error for procs'
       call endrun
    end if

    allocate(abegg(0:npes-1),aendg(0:npes-1), stat=ier)
    if (ier /= 0) then
       write (6,*) 'initDecomp(): allocation error for abegg,aendg'
       call endrun
    end if

    ! Find total global number of grid cells, landunits, columns and pfts

    ancells = 0
    ncells  = 0
    nlunits = 0
    ncols   = 0
    npfts   = 0
    do aj = 1, anj
    do ai = 1, ani
       if (adomain%mask(ai,aj) == 1) then
          ancells  = ancells  + 1
          !--- check again that cell in adomain has ldomain valid underneath
          if (n_ovr(ai,aj) < 1) then
             write (6,*) 'initDecomp(): map overlap error at ',ai,aj,adomain%mask(ai,aj),n_ovr(ai,aj)
             call endrun()
          endif
          do n = 1,n_ovr(ai,aj)
             i = i_ovr(ai,aj,n)
             j = j_ovr(ai,aj,n)
             call get_gcell_info (i, j, wtxy, nlunits=ilunits, &
                                  ncols=icols, npfts=ipfts)
             ncells  = ncells  + 1
             nlunits = nlunits + ilunits
             ncols   = ncols   + icols
             npfts   = npfts   + ipfts
             lcell(i,j) = .true.
          enddo
       end if
    end do
    end do

    ! Allocate dynamic memory for adecomp, ldecomp derived type

    allocate(adecomp%gdc2gsn(ancells), adecomp%gsn2gdc(ancells), &
             adecomp%gdc2i  (ancells), adecomp%gdc2j  (ancells), &
             adecomp%gsn2i  (ancells), adecomp%gsn2j  (ancells), &
             stat=ier)
    if (ier /= 0) then
       write (6,*) 'initDecomp(): allocation error1 for adecomp'
       call endrun()
    end if

    allocate(adecomp%ij2gsn(ani,anj), adecomp%ij2gdc(ani,anj), &
             stat=ier)
    if (ier /= 0) then
       write (6,*) 'initDecomp(): allocation error2 for adecomp'
       call endrun()
    end if

    adecomp%gdc2gsn(:)  = 0
    adecomp%gsn2gdc(:)  = 0
    adecomp%gsn2i(:)    = 0
    adecomp%gsn2j(:)    = 0
    adecomp%gdc2i(:)    = 0
    adecomp%gdc2j(:)    = 0
    adecomp%ij2gsn(:,:) = 0
    adecomp%ij2gdc(:,:) = 0

    allocate(ldecomp%gdc2gsn(ncells), ldecomp%gsn2gdc(ncells), &
             ldecomp%gdc2i  (ncells), ldecomp%gdc2j  (ncells), &
             ldecomp%gsn2i  (ncells), ldecomp%gsn2j  (ncells), &
             stat=ier)
    if (ier /= 0) then
       write (6,*) 'initDecomp(): allocation error1 for ldecomp'
       call endrun()
    end if

    allocate(ldecomp%ij2gsn(lni,lnj), ldecomp%ij2gdc(lni,lnj), &
             stat=ier)
    if (ier /= 0) then
       write (6,*) 'initDecomp(): allocation error2 for ldecomp'
       call endrun()
    end if

    ldecomp%gdc2gsn(:)  = 0
    ldecomp%gsn2gdc(:)  = 0
    ldecomp%gsn2i(:)    = 0
    ldecomp%gsn2j(:)    = 0
    ldecomp%gdc2i(:)    = 0
    ldecomp%gdc2j(:)    = 0
    ldecomp%ij2gsn(:,:) = 0
    ldecomp%ij2gdc(:,:) = 0

    ! Error check on total number of gridcells

    if (npes > ncells) then
       write (6,*) 'initDecomp(): Number of processes exceeds number ', &
            'of land grid cells'
       call endrun
    end if

    ! Diagnostic output

    if (masterproc) then
       write (6,*)' Atm Grid Characteristics'
       write (6,*)'   longitude points          = ',ani
       write (6,*)'   latitude points           = ',anj
       write (6,*)'   total number of gridcells = ',ancells
       write (6,*)' Surface Grid Characteristics'
       write (6,*)'   longitude points          = ',lni
       write (6,*)'   latitude points           = ',lnj
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
    allocate(clumps_ancells(1:nclumps), stat=ier)
    if (ier /= 0) then
       write (6,*) 'initDecomp(): allocation error for clumps_ancells'
       call endrun()
    end if
    allocate(clumpfull(1:nclumps), stat=ier)
    if (ier /= 0) then
       write (6,*) 'initDecomp(): allocation error for temporary clump arrays'
       call endrun()
    end if
    allocate(acidi(1:nclumps,ancells),acidj(1:nclumps,ancells),stat=ier)
    if (ier /= 0) then
       write (6,*) 'initDecomp(): allocation error for acidi, acidj'
       call endrun()
    end if
    allocate(cidi(1:nclumps,ncells),cidj(1:nclumps,ncells),stat=ier)
    if (ier /= 0) then
       write (6,*) 'initDecomp(): allocation error for cidi, cidj'
       call endrun()
    end if


    ! Determine number of grid cells in each clump - assign gridcells to clumps
    ! in a round robin fashion

    acidi(:,:) = 0
    acidj(:,:) = 0
    clumps_ancells(:) = 0
    clumps(:)%ncells  = 0
    clumps(:)%nlunits = 0
    clumps(:)%ncols   = 0
    clumps(:)%npfts   = 0
    clumpfull(:) = .false.
    cidi(:,:) = 0
    cidj(:,:) = 0
    cid = 0

    ! Set ldecomp sn indexing based on cells to be used and i,j order
    g  = 0
    do j = 1,lnj
    do i = 1,lni
       if (lcell(i,j)) then
          g = g + 1
          ldecomp%ij2gsn(i,j) = g
          ldecomp%gsn2i(g) = i
          ldecomp%gsn2j(g) = j
       endif
    enddo
    enddo

    ag = 0
    error = .false.
    do aj = 1,anj
    do ai = 1,ani
       if (adomain%mask(ai,aj) == 1) then

          ag = ag + 1
          adecomp%ij2gsn(ai,aj) = ag
          adecomp%gsn2i(ag) = ai
          adecomp%gsn2j(ag) = aj

          ! Set clump id for gridcell - if clump is full, then determine
          ! first non-full clump

          clumpcount = 0
          validclump = .false.
          do while (.not. validclump .and. clumpcount <= nclumps)
             clumpcount = clumpcount + 1
             cid = cid + 1
             if (cid == nclumps+1) cid = 1
             if (.not. clumpfull(cid)) validclump = .true.
          end do
          if (.not.validclump)  error = .true.
          if (n_ovr(ai,aj) < 1) error = .true.
          if (error) exit

          clumps_ancells(cid) = clumps_ancells(cid) + 1
          acidi(cid,clumps_ancells(cid)) = ai
          acidj(cid,clumps_ancells(cid)) = aj

          do n = 1,n_ovr(ai,aj)
             i = i_ovr(ai,aj,n)
             j = j_ovr(ai,aj,n)

             ! Determine grid cell info
             call get_gcell_info (i, j, wtxy, nlunits=ilunits, &
                                  ncols=icols, npfts=ipfts)

             clumps(cid)%ncells  = clumps(cid)%ncells  + 1
             clumps(cid)%nlunits = clumps(cid)%nlunits + ilunits
             clumps(cid)%ncols   = clumps(cid)%ncols   + icols
             clumps(cid)%npfts   = clumps(cid)%npfts   + ipfts

             ! set tmps, cid,ncell <-> i,j

             cidi(cid,clumps(cid)%ncells) = i
             cidj(cid,clumps(cid)%ncells) = j
             if (clumps(cid)%npfts >= ppc .and. cid < nclumps) then
                clumpfull(cid) = .true.
             end if
          end do

          if (error) exit
       endif
    end do
    end do
    if (error) then
       write(6,*)'initDecomp: error encountered trying to find an unfull clump'
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
    gf = 0
    lf = 0
    cf = 0
    pf = 0
    do pid = 0,npes-1
       do n = 1,procs(pid)%nclumps
          cid = procs(pid)%cid(n)
          gi = gf + 1
          gf = gi + clumps(cid)%ncells  - 1
          li = lf + 1
          lf = li + clumps(cid)%nlunits - 1
          ci = cf + 1
          cf = ci + clumps(cid)%ncols   - 1
          pi = pf + 1
          pf = pi + clumps(cid)%npfts   - 1
          clumps(cid)%begg = gi
          clumps(cid)%endg = gf
          clumps(cid)%begl = li
          clumps(cid)%endl = lf
          clumps(cid)%begc = ci
          clumps(cid)%endc = cf
          clumps(cid)%begp = pi
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

    ! Set ldecomp and adecomp data

    ag = 0
    g = 0
    do pid = 0,npes-1
       abegg(pid) = ag+1
       do nc=1,procs(pid)%nclumps
          cid = procs(pid)%cid(nc)

          do n = 1,clumps_ancells(cid)
             ag=ag+1
             i = acidi(cid,n)
             j = acidj(cid,n)
             if (i==0.or.j==0) then
                write(6,*)'initDecomp ERROR: acidi,acidj ',pid,nc,cid,n,ag,i,j
                call endrun()
             endif
             adecomp%gdc2i(ag) = i
             adecomp%gdc2j(ag) = j
             adecomp%gdc2gsn(ag) = adecomp%ij2gsn(i,j)
             adecomp%ij2gdc(i,j) = ag
             adecomp%gsn2gdc(adecomp%ij2gsn(i,j)) = ag
          enddo

          do n = 1,clumps(cid)%ncells
             g=g+1
             i = cidi(cid,n)
             j = cidj(cid,n)
             if (i==0.or.j==0) then
                write(6,*)'initDecomp ERROR: cidi,cidj ',pid,nc,cid,n,g,i,j
                call endrun()
             endif
             ldecomp%gdc2i(g) = i
             ldecomp%gdc2j(g) = j
             ldecomp%gdc2gsn(g) = ldecomp%ij2gsn(i,j)
             ldecomp%ij2gdc(i,j) = g
             ldecomp%gsn2gdc(ldecomp%ij2gsn(i,j)) = g
          enddo

       enddo
       aendg(pid) = ag
    enddo

    ! Write out clump and proc info

    if (masterproc) then
       do pid = 0,npes-1
          write(6,*)
          write(6,*)'proc= ',pid,' beg atmcell = ',abegg(pid),      &
               ' end atmcell = ',aendg(pid),                        &
               ' total atmcells per proc = ',aendg(pid)-abegg(pid)+1
          write(6,*)'proc= ',pid,' beg gridcell= ',procs(pid)%begg, &
               ' end gridcell= ',procs(pid)%endg,                   &
               ' total gridcells per proc= ',procs(pid)%ncells
          write(6,*)'proc= ',pid,' beg landunit= ',procs(pid)%begl, &
               ' end landunit= ',procs(pid)%endl,                   &
               ' total landunits per proc = ',procs(pid)%nlunits
          write(6,*)'proc= ',pid,' beg column  = ',procs(pid)%begc, &
               ' end column  = ',procs(pid)%endc,                   &
               ' total columns per proc  = ',procs(pid)%ncols
          write(6,*)'proc= ',pid,' beg pft     = ',procs(pid)%begp, &
               ' end pft     = ',procs(pid)%endp,                   &
               ' total pfts per proc     = ',procs(pid)%npfts
          do n = 1,procs(pid)%nclumps
             cid = procs(pid)%cid(n)
             write(6,*)'proc= ',pid,' clump no = ',n, &
                  ' clump id= ',procs(pid)%cid(n),    &
                  ' beg gridcell= ',clumps(cid)%begg, &
                  ' end gridcell= ',clumps(cid)%endg, &
                  ' total gridcells per clump= ',clumps(cid)%ncells
             write(6,*)'proc= ',pid,' clump no = ',n, &
                  ' clump id= ',procs(pid)%cid(n),    &
                  ' beg landunit= ',clumps(cid)%begl, &
                  ' end landunit= ',clumps(cid)%endl, &
                  ' total landunits per clump = ',clumps(cid)%nlunits
             write(6,*)'proc= ',pid,' clump no = ',n, &
                  ' clump id= ',procs(pid)%cid(n),    &
                  ' beg column  = ',clumps(cid)%begc, &
                  ' end column  = ',clumps(cid)%endc, &
                  ' total columns per clump  = ',clumps(cid)%ncols
             write(6,*)'proc= ',pid,' clump no = ',n, &
                  ' clump id= ',procs(pid)%cid(n),    &
                  ' beg pft     = ',clumps(cid)%begp, &
                  ' end pft     = ',clumps(cid)%endp, &
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

    deallocate(clumps_ancells)
    deallocate(clumpfull)
    deallocate(acidi,acidj)
    deallocate(cidi,cidj)
    deallocate(lcell)

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
! !IROUTINE: get_clump_owner_id
!
! !INTERFACE:
   integer function get_clump_owner_id(cid)
!
! !DESCRIPTION:
! This function returns the MPI process id (rank) responsible for the clump
! identified by the clump id cid.
!
! !USES:
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
! !USES:
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
! !USES:
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
! !USES:
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
     integer :: gdc                                 ! global dc index
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
        gdc = clumps(cid)%begg + i - 1
        lons(i) = ldecomp%gdc2i(gdc)
        lats(i) = ldecomp%gdc2j(gdc)
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
! !USES:
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

     gi = clumps(cid)%begg+cell-1

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
! !USES:
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
! !USES:
!
! !ARGUMENTS:
     implicit none
     integer, optional, intent(out) :: begp, endp  ! proc beginning and ending
                                                   ! pft indices
     integer, optional, intent(out) :: begc, endc  ! proc beginning and ending
                                                   ! column indices
     integer, optional, intent(out) :: begl, endl  ! proc beginning and ending
                                                   ! landunit indices
     integer, optional, intent(out) :: begg, endg  ! proc beginning and ending
                                                   ! gridcell indices
! !DESCRIPTION:
! Retrieve gridcell, landunit, column, and pft bounds for process.
!
! !REVISION HISTORY:
! 2003.09.12  Mariana Vertenstein  Creation.
!
!EOP
!------------------------------------------------------------------------------

     if (present(begp)) then
        begp = procs(iam)%begp
     endif
     if (present(endp)) then
        endp = procs(iam)%endp
     endif
     if (present(begc)) then
        begc = procs(iam)%begc
     endif
     if (present(endc)) then
        endc = procs(iam)%endc
     endif
     if (present(begl)) then
        begl = procs(iam)%begl
     endif
     if (present(endl)) then
        endl = procs(iam)%endl
     endif
     if (present(begg)) then
        begg = procs(iam)%begg
     endif
     if (present(endg)) then
        endg = procs(iam)%endg
     endif

   end subroutine get_proc_bounds

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_proc_bounds_atm
!
! !INTERFACE:
   subroutine get_proc_bounds_atm (begg, endg)
!
! !USES:
!
! !ARGUMENTS:
     implicit none
     integer, intent(out) :: begg, endg  ! proc beginning and ending
                                         ! gridcell indices for atm grid
! !DESCRIPTION:
! Retrieve gridcell begg, endg for atm decomp
!
! !REVISION HISTORY:
! 2005.12.15  T Craig Added
!
!EOP
!------------------------------------------------------------------------------

   begg = abegg(iam)
   endg = aendg(iam)

   end subroutine get_proc_bounds_atm

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
! !USES:
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
! !USES:
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
! !USES:
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
     integer :: n,l
!------------------------------------------------------------------------------
     do gdc = 1,procs(iam)%numg
        gsn = ldecomp%gdc2gsn(gdc)
        call map_indexes(gdc,gsn,type1d,dci,dcf,sni,snf)
        n = dcf - dci + 1
        do l = 1,n
           sn = sni + l - 1
           dc = dci + l - 1
           arraysn(sn) = arraydc(dc)
        end do
     end do

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
     integer :: n,l
!------------------------------------------------------------------------------

     do gdc = 1,procs(iam)%numg
        gsn = ldecomp%gdc2gsn(gdc)
        call map_indexes(gdc,gsn,type1d,dci,dcf,sni,snf)
        n = dcf - dci + 1
        do l = 1,n
           sn = sni + l - 1
           dc = dci + l - 1
           arraysn(sn) = arraydc(dc)
        end do
     end do

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
     integer :: n,l
!------------------------------------------------------------------------------

     do gdc = 1,procs(iam)%numg
        gsn = ldecomp%gdc2gsn(gdc)
        call map_indexes(gdc,gsn,type1d,dci,dcf,sni,snf)
        n = dcf - dci + 1
        do l = 1,n
           sn = sni + l - 1
           dc = dci + l - 1
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
     end do

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
     integer :: n,l
!------------------------------------------------------------------------------

     do gdc = 1,procs(iam)%numg
        gsn = ldecomp%gdc2gsn(gdc)
        call map_indexes(gdc,gsn,type1d,dci,dcf,sni,snf)
        n = dcf - dci + 1
        do l = 1,n
           sn = sni + l - 1
           dc = dci + l - 1
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
     end do

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
     integer :: n,l
!------------------------------------------------------------------------------

     do gsn = 1,procs(iam)%numg
        gdc = ldecomp%gsn2gdc(gsn)
        call map_indexes(gdc,gsn,type1d,dci,dcf,sni,snf)
        n = dcf - dci + 1
        do l = 1,n
           sn = sni + l - 1
           dc = dci + l - 1
           arraydc(dc) = arraysn(sn)
        end do
     end do

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
     integer :: n,l
!------------------------------------------------------------------------------

     do gsn = 1,procs(iam)%numg
        gdc = ldecomp%gsn2gdc(gsn)
        call map_indexes(gdc,gsn,type1d,dci,dcf,sni,snf)
        n = dcf - dci + 1
        do l = 1,n
           sn = sni + l - 1
           dc = dci + l - 1
           arraydc(dc) = arraysn(sn)
        end do
     end do

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
     integer :: n,l,lev

!------------------------------------------------------------------------------

     do gsn = 1,procs(iam)%numg
        gdc = ldecomp%gsn2gdc(gsn)
        call map_indexes(gdc,gsn,type1d,dci,dcf,sni,snf)
        n = dcf - dci + 1
        do l = 1,n
           sn = sni + l - 1
           dc = dci + l - 1
           do lev = lb1, ub1
              arraydc(lev,dc) = arraysn(lev,sn)
           end do
        end do
     end do

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
     integer :: n,l,lev
!------------------------------------------------------------------------------

     do gsn = 1,procs(iam)%numg
        gdc = ldecomp%gsn2gdc(gsn)
        call map_indexes(gdc,gsn,type1d,dci,dcf,sni,snf)
        n = dcf - dci + 1
        do l = 1,n
           sn = sni + l - 1
           dc = dci + l - 1
           do lev = lb1, ub1
              arraydc(lev,dc) = arraysn(lev,sn)
           end do
        end do
     end do

   end subroutine map_sn2dc_ml1_int

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: map_indexes
!
! !INTERFACE:
   subroutine map_indexes(gdc,gsn,type1d,dci,dcf,sni,snf)
!
! !DESCRIPTION:
!  Gets indices for dc2sn mapping routines
!
! !USES:
   use clmtype, only : nameg, namel, namec, namep
   use initSubgridMod, only : gcellsn,gcelldc
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: gdc,gsn
    character(len=*), intent(in) :: type1d
    integer, intent(out) :: dci,dcf,sni,snf
!
! !REVISION HISTORY:
! 2005.11.15  T Craig Extracted from map_*_* subroutines
!
!EOP
!
! !LOCAL VARIABLES:
    logical :: error
!------------------------------------------------------------------------------

     error = .false.
     if (type1d == nameg) then
        dci = gdc
        dcf = gdc
        sni = gsn
        snf = gsn
     else if (type1d == namel) then
        dci = gcelldc%g_li(gdc)
        dcf = gcelldc%g_lf(gdc)
        sni = gcellsn%g_li(gsn)
        snf = gcellsn%g_lf(gsn)
     else if (type1d == namec) then
        dci = gcelldc%g_ci(gdc)
        dcf = gcelldc%g_cf(gdc)
        sni = gcellsn%g_ci(gsn)
        snf = gcellsn%g_cf(gsn)
     else if (type1d == namep) then
        dci = gcelldc%g_pi(gdc)
        dcf = gcelldc%g_pf(gdc)
        sni = gcellsn%g_pi(gsn)
        snf = gcellsn%g_pf(gsn)
     end if
     if (dcf - dci /= snf - sni) then
        write(6,*)'error in map_indexes ',gsn,gdc,trim(type1d),dci,dcf,sni,snf
        call endrun()
     end if

end subroutine map_indexes

!------------------
end module decompMod
