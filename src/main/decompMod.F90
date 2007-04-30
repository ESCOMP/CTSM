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
  use spmdMod     , only : masterproc, iam, npes, mpicom, comp_id
  use clm_mct_mod
  use shr_sys_mod , only : shr_sys_flush
  use abortutils  , only : endrun
!
! !PUBLIC TYPES:
  implicit none
  integer, public :: clump_pproc ! number of clumps per MPI process
!
! !PUBLIC MEMBER FUNCTIONS:
!  public decomp_init             ! initializes land surface decomposition
!                                 ! into clumps and processors
  public decomp_atm_init         ! initializes atm grid decomposition
                                 ! into clumps and processors
  public decomp_lnd_init         ! initializes lnd grid decomposition
                                 ! into clumps and processors
  public decomp_glcp_init        ! initializes g,l,c,p decomp info
  public decomp_domg2l           ! create local domain from global domain
  public get_clump_bounds        ! beg and end gridcell, landunit, column,
                                 ! pft indices for clump
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
! 2006.06.06  T Craig  Reduce memory, cleanup
!
!EOP
!
! !PRIVATE TYPES:
  private

  integer :: nclumps     ! total number of clumps across all processors
  integer :: numg        ! total number of gridcells on all procs
  integer :: numl        ! total number of landunits on all procs
  integer :: numc        ! total number of columns on all procs
  integer :: nump        ! total number of pfts on all procs
  integer :: numa        ! total number of atm gridcells on all procs
  integer, pointer :: acid(:)       ! temporary for setting adecomp/ldecomp

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
     integer :: abegg,aendg      ! beginning and ending atm gridcell index
  end type processor_type
  type(processor_type) :: procinfo

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
  !--- i,j = 2d global
  !--- glo = 1d global sn ordered
  !--- gsn = 1d global sn ordered compressed
  !--- gdc = 1d global dc ordered compressed
  type decomp_type
     integer,pointer :: gsn2gdc(:)    ! 1d gsn to 1d gdc
     integer,pointer :: gdc2gsn(:)    ! 1d gdc to 1d gsn
     integer,pointer :: glo2gdc(:)    ! 1d glo to 1d gdc
     integer,pointer :: gdc2glo(:)    ! 1d gdc to 1d glo
     integer,pointer :: gdc2i(:)      ! 1d gdc to 2d sn i index
     integer,pointer :: gdc2j(:)      ! 1d gdc to 2d sn j index
  end type decomp_type
  public decomp_type
  type(decomp_type),public,target :: ldecomp
  type(decomp_type),public,target :: adecomp

  type(mct_gsMap)  ,public :: gsMap_lnd_gdc2glo
  integer,pointer  ,public ::  perm_lnd_gdc2glo(:)
  type(mct_gsMap)  ,public :: gsMap_atm_gdc2glo
  integer,pointer  ,public ::  perm_atm_gdc2glo(:)

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
! 
! 
!------------------------------------------------------------------------------

contains

!------------------------------------------------------------------------------
#if (1 == 0)
! tcx DO NOT DELETE THIS YET
!BOP
!
! !IROUTINE: decomp_init
!
! !INTERFACE:
  subroutine decomp_init(wtxy)
!
! !DESCRIPTION:
! This subroutine initializes the land surface decomposition into a clump
! data structure.  This assumes each pe has the same number of clumps
! set by clump_pproc
!
! !USES:
    use domainMod , only : ldomain,adomain
    use subgridMod, only : subgrid_get_gcellinfo
    use clm_varctl, only : nsegspc
!
! !ARGUMENTS:
    implicit none
    real(r8), intent(in) :: wtxy(:,:)   ! subgrid patch weights
!
! !LOCAL VARIABLES:
    integer :: lni,lnj                ! land domain global size
    integer :: ani,anj                ! atm domain global size
    integer :: lns,lg,ln,li,lj        ! indices
    integer :: ans,ag,an,ai,aj        ! indices
    integer :: anumg                  ! atm num gridcells
    integer :: anumg_tot              ! precompute of anumg
    real(r8):: rnsegspc               ! real value associated with nsegspc
!    integer :: anumpc                 ! min atm gridcells/clump
!    integer :: anumxtra               ! extra atm cells
    integer :: cid,pid                ! indices
    integer, pointer :: lcid(:)       ! temporary for setting adecomp
    integer, pointer :: acid(:)       ! temporary for setting adecomp
    integer :: n,m,np                 ! indices
    integer :: ilunits, icols, ipfts  ! temporaries
    integer :: ier                    ! error code
    integer :: cnt                    ! local counter
    integer, parameter :: dbug=3      ! 0 = min, 1=normal, 2=much, 3=max
    integer :: npmin,npmax,npint      ! do loop values for printing
    integer :: clmin,clmax,clint      ! do loop values for printing
    integer :: beg,end,lsize,gsize    ! used for gsmap init
    integer, pointer :: gindex(:)     ! global index for gsmap init

    integer, pointer :: lncnt(:)      ! lnd cell count per atm cell
    integer, pointer :: lnoff(:)      ! atm cell offset in lnmap
    integer, pointer :: lnmap(:)      ! map from atm cell to lnd cells
    integer, pointer :: lglo2gsn(:)   ! map from glo 2 gsn temporary
    integer, pointer :: aglo2gsn(:)   ! map from glo 2 gsn temporary
    integer :: lnidx

! !CALLED FROM:
! subroutine initialize
!
! !REVISION HISTORY:
! 2002.09.11  Forrest Hoffman  Creation.
! 2005.12.15  T Craig  Updated for finemesh
! 2006.08.18  P Worley Performance optimizations
!
!EOP
!------------------------------------------------------------------------------

    lni = ldomain%ni
    lnj = ldomain%nj
    ani = adomain%ni
    anj = adomain%nj
    lns = ldomain%ns
    ans = adomain%ns

    !--- set and verify nclumps ---
    if (clump_pproc > 0) then
       nclumps = clump_pproc * npes
       if (nclumps < npes) then
          write (6,*) 'decomp_init(): Number of gridcell clumps= ',nclumps, &
               ' is less than the number of processes = ', npes
          call endrun()
       end if
    else
       write(6,*)'clump_pproc= ',clump_pproc,'  must be greater than 0'
       call endrun()
    end if

    !--- allocate and initialize procinfo and clumps ---
    !--- beg and end indices initialized for simple addition of cells later ---

    allocate(procinfo%cid(clump_pproc), stat=ier)
    if (ier /= 0) then
       write (6,*) 'decomp_init(): allocation error for procinfo%cid'
       call endrun()
    endif

    procinfo%nclumps = clump_pproc
    procinfo%cid(:)  = -1
    procinfo%ncells  = 0
    procinfo%nlunits = 0
    procinfo%ncols   = 0
    procinfo%npfts   = 0
    procinfo%begg    = 1
    procinfo%begl    = 1
    procinfo%begc    = 1
    procinfo%begp    = 1
    procinfo%endg    = 0
    procinfo%endl    = 0
    procinfo%endc    = 0
    procinfo%endp    = 0
    procinfo%abegg   = 1
    procinfo%aendg   = 0

    allocate(clumps(nclumps), stat=ier)
    if (ier /= 0) then
       write (6,*) 'decomp_init(): allocation error for clumps'
       call endrun()
    end if
    clumps(:)%owner   = -1
    clumps(:)%ncells  = 0
    clumps(:)%nlunits = 0
    clumps(:)%ncols   = 0
    clumps(:)%npfts   = 0
    clumps(:)%begg    = 1
    clumps(:)%begl    = 1
    clumps(:)%begc    = 1
    clumps(:)%begp    = 1
    clumps(:)%endg    = 0
    clumps(:)%endl    = 0
    clumps(:)%endc    = 0
    clumps(:)%endp    = 0

    !--- assign clumps to proc round robin ---
    cid = 0
    do n = 1,nclumps
       pid = mod(n-1,npes)
       if (pid < 0 .or. pid > npes-1) then
          write (6,*) 'decomp_init(): round robin pid error ',n,pid,npes
          call endrun()
       endif
       clumps(n)%owner = pid
       if (iam == pid) then
          cid = cid + 1
          if (cid < 1 .or. cid > clump_pproc) then
             write (6,*) 'decomp_init(): round robin pid error ',n,pid,npes
             call endrun()
          endif
          procinfo%cid(cid) = n
       endif
    enddo

    allocate(lncnt(ans),lnoff(ans),lnmap(lns))

    lncnt = 0
    do ln = 1,lns
       an = ldomain%gatm(ln)
       if ((an > 0) .and. (an .le. ans)) then
          lncnt(an) = lncnt(an) + 1
       endif
    enddo

    lnoff(1) = 1
    do an = 2,ans
       lnoff(an) = lnoff(an-1) + lncnt(an-1)
    enddo

    lncnt = 0
    lnmap = -1
    do ln = 1,lns
       an = ldomain%gatm(ln)
       if ((an > 0) .and. (an .le. ans)) then
         lnmap(lnoff(an)+lncnt(an)) = ln
         lncnt(an) = lncnt(an) + 1
       endif
    enddo

!    !--- count total pfts, like loop below, in case you need them "early"
!    nump  = 0
    anumg_tot = 0
    do an = 1,ans
       if (adomain%mask(an) == 1) then
          anumg_tot  = anumg_tot  + 1
!          do lnidx = 0,lncnt(an)-1
!             ln = lnmap(lnoff(an)+lnidx)          
!             call subgrid_get_gcellinfo (ln, wtxy, nlunits=ilunits, &
!                                  ncols=icols, npfts=ipfts)
!             nump = nump + ipfts
!          enddo
       endif
    enddo
!    nsegspc  = 1000000                  ! number of segments/clump
    rnsegspc = min(float(nsegspc),float(anumg_tot)/float(nclumps))
!    anumpc   = anumg_tot/nclumps
!    anumxtra = mod(anumg_tot,nclumps)

!    if (masterproc) write(6,*) 'precompute total pfts ',nump
!    if (masterproc) write(6,*) 'precompute total anumg ',anumg_tot,anumpc,anumxtra
    if (masterproc) write(6,*) 'precompute total anumg ',anumg_tot,rnsegspc

    !--- assign gridcells to clumps (and thus pes) ---
    allocate(lcid(lns),acid(ans))
    lcid = 0
    acid = 0
    anumg = 0
    numg  = 0
    numl  = 0
    numc  = 0
    nump  = 0
    do an = 1,ans
       if (adomain%mask(an) == 1) then
          anumg  = anumg  + 1

#if (1 == 0) 
          !--- find clump with fewest pfts ---
          cid = 1
          do n = 2,nclumps
             if (clumps(n)%npfts < clumps(cid)%npfts) then
                cid = n
             endif
          enddo
#endif
#if (1 == 0)
          !--- give to clumps in simple round robin
          cid = mod((anumg-1),nclumps) + 1
#endif
#if (1 == 1)
          !--- give to clumps in order based on nsegspc
          cid = int(rnsegspc*float(nclumps*(anumg-1))/float(anumg_tot))
          cid = mod(cid,nclumps) + 1
#endif
          acid(an) = cid

          !--- give atm cell to pe that owns cid ---
          if (iam >  clumps(cid)%owner) then
             procinfo%abegg = procinfo%abegg + 1
          endif
          if (iam >= clumps(cid)%owner) then
             procinfo%aendg = procinfo%aendg + 1
          endif

          cnt = 0
          do lnidx = 0,lncnt(an)-1
             ln = lnmap(lnoff(an)+lnidx)          
             cnt = cnt + 1
             call subgrid_get_gcellinfo (ln, wtxy, nlunits=ilunits, &
                                  ncols=icols, npfts=ipfts)
             lcid(ln) = cid

             !--- overall total ---
             numg = numg + 1
             numl = numl + ilunits
             numc = numc + icols
             nump = nump + ipfts

             !--- give gridcell to cid ---
             !--- increment the beg and end indices ---
             clumps(cid)%ncells  = clumps(cid)%ncells  + 1
             clumps(cid)%nlunits = clumps(cid)%nlunits + ilunits
             clumps(cid)%ncols   = clumps(cid)%ncols   + icols
             clumps(cid)%npfts   = clumps(cid)%npfts   + ipfts

             do m = 1,nclumps
                if ((clumps(m)%owner >  clumps(cid)%owner) .or. &
                    (clumps(m)%owner == clumps(cid)%owner .and. m > cid)) then
                   clumps(m)%begg = clumps(m)%begg + 1
                   clumps(m)%begl = clumps(m)%begl + ilunits
                   clumps(m)%begc = clumps(m)%begc + icols
                   clumps(m)%begp = clumps(m)%begp + ipfts
                endif

                if ((clumps(m)%owner >  clumps(cid)%owner) .or. &
                    (clumps(m)%owner == clumps(cid)%owner .and. m >= cid)) then
                   clumps(m)%endg = clumps(m)%endg + 1
                   clumps(m)%endl = clumps(m)%endl + ilunits
                   clumps(m)%endc = clumps(m)%endc + icols
                   clumps(m)%endp = clumps(m)%endp + ipfts
                endif
             enddo

             !--- give gridcell to the proc that owns the cid ---
             !--- increment the beg and end indices ---
             if (iam == clumps(cid)%owner) then
                procinfo%ncells  = procinfo%ncells  + 1
                procinfo%nlunits = procinfo%nlunits + ilunits
                procinfo%ncols   = procinfo%ncols   + icols
                procinfo%npfts   = procinfo%npfts   + ipfts
             endif

             if (iam >  clumps(cid)%owner) then
                procinfo%begg = procinfo%begg + 1
                procinfo%begl = procinfo%begl + ilunits
                procinfo%begc = procinfo%begc + icols
                procinfo%begp = procinfo%begp + ipfts
             endif

             if (iam >= clumps(cid)%owner) then
                procinfo%endg = procinfo%endg + 1
                procinfo%endl = procinfo%endl + ilunits
                procinfo%endc = procinfo%endc + icols
                procinfo%endp = procinfo%endp + ipfts
             endif
          enddo
          !--- check that atm cell has at least 1 lnd grid cell
          if (cnt < 1) then
             write (6,*) 'decomp_init(): map overlap error at ',an, &
                adomain%mask(an),cnt
             call endrun()
          endif
       end if
    enddo

    ! Error check on total number of gridcells

    if (npes > anumg) then
       write (6,*) 'decomp_init(): Number of processes exceeds number ', &
            'of atm grid cells'
       call endrun()
    end if

    ! Allocate dynamic memory for adecomp, ldecomp derived type

    allocate(adecomp%gdc2gsn(anumg), adecomp%gsn2gdc(anumg), &
             adecomp%gdc2glo(anumg), adecomp%glo2gdc(ani*anj), &
             adecomp%gdc2i  (anumg), adecomp%gdc2j  (anumg), &
             stat=ier)
    if (ier /= 0) then
       write (6,*) 'decomp_init(): allocation error1 for adecomp'
       call endrun()
    end if

    adecomp%gdc2gsn(:)  = 0
    adecomp%gsn2gdc(:)  = 0
    adecomp%gdc2glo(:)  = 0
    adecomp%glo2gdc(:)  = 0
    adecomp%gdc2i(:)    = 0
    adecomp%gdc2j(:)    = 0

    allocate(ldecomp%gdc2gsn(numg), ldecomp%gsn2gdc(numg), &
             ldecomp%gdc2glo(numg), ldecomp%glo2gdc(lni*lnj), &
             ldecomp%gdc2i  (numg), ldecomp%gdc2j  (numg), &
             stat=ier)
    if (ier /= 0) then
       write (6,*) 'decomp_init(): allocation error1 for ldecomp'
       call endrun()
    end if

    ldecomp%gdc2gsn(:)  = 0
    ldecomp%gsn2gdc(:)  = 0
    ldecomp%gdc2glo(:)  = 0
    ldecomp%glo2gdc(:)  = 0
    ldecomp%gdc2i(:)    = 0
    ldecomp%gdc2j(:)    = 0

    !--- temporaries for decomp mappings
    allocate(aglo2gsn(ani*anj),lglo2gsn(lni*lnj),stat=ier)
    if (ier /= 0) then
       write (6,*) 'decomp_init(): allocation error for al-glo2gsn'
       call endrun()
    end if
    aglo2gsn(:) = 0
    lglo2gsn(:) = 0

    ag = 0
    do aj = 1,anj
    do ai = 1,ani
       an = (aj-1)*ani + ai
       if (acid(an) > 0) then
          ag  = ag  + 1
          aglo2gsn(an) = ag
       endif
    enddo
    enddo

    ! Set ldecomp sn indexing based on cells to be used and i,j order
    lg  = 0
    do lj = 1,lnj
    do li = 1,lni
       ln = (lj-1)*lni + li
       if (lcid(ln) > 0) then
          lg = lg + 1
          lglo2gsn(ln) = lg
       endif
    enddo
    enddo

    ! Set ldecomp and adecomp data
    ag = 0
    lg = 0
    do pid = 0,npes-1
    do cid = 1,nclumps
       if (clumps(cid)%owner == pid) then

          do aj = 1,anj
          do ai = 1,ani
             an = (aj-1)*ani + ai
             if (acid(an) == cid) then
                ag = ag + 1
                adecomp%gdc2i(ag) = ai
                adecomp%gdc2j(ag) = aj
                adecomp%gdc2gsn(ag) = aglo2gsn(an)
                adecomp%gdc2glo(ag) = an
                adecomp%gsn2gdc(aglo2gsn(an)) = ag
                adecomp%glo2gdc(an) = ag
             endif
          enddo
          enddo

          do lj = 1,lnj
          do li = 1,lni
             ln = (lj-1)*lni + li
             if (lcid(ln) == cid) then
                lg = lg + 1
                ldecomp%gdc2i(lg) = li
                ldecomp%gdc2j(lg) = lj
                ldecomp%gdc2gsn(lg) = lglo2gsn(ln)
                ldecomp%gdc2glo(lg) = ln
                ldecomp%gsn2gdc(lglo2gsn(ln)) = lg
                ldecomp%glo2gdc(ln) = lg
             endif
          enddo
          enddo
       endif
    enddo
    enddo

    deallocate(aglo2gsn,lglo2gsn)
    deallocate(acid,lcid)
    deallocate(lncnt,lnoff,lnmap)

    ! set gsMap_lnd_gdc2glo, perm_lnd_gdc2glo
    call get_proc_bounds(beg, end)
    allocate(gindex(beg:end))
    do n = beg,end
       gindex(n) = ldecomp%gdc2glo(n)
    enddo
    lsize = end-beg+1
    gsize = lni * lnj
    allocate(perm_lnd_gdc2glo(lsize),stat=ier)
    call mct_indexset(perm_lnd_gdc2glo)
    call mct_indexsort(lsize,perm_lnd_gdc2glo,gindex)
    call mct_permute(gindex,perm_lnd_gdc2glo,lsize)
    call mct_gsMap_init(gsMap_lnd_gdc2glo, gindex, mpicom, comp_id, lsize, gsize )
    deallocate(gindex)

    ! set gsMap_atm_gdc2glo, perm_atm_gdc2glo
    call get_proc_bounds_atm(beg, end)
    allocate(gindex(beg:end))
    do n = beg,end
       gindex(n) = adecomp%gdc2glo(n)
    enddo
    lsize = end-beg+1
    gsize = ani * anj
    allocate(perm_atm_gdc2glo(lsize),stat=ier)
    call mct_indexset(perm_atm_gdc2glo)
    call mct_indexsort(lsize,perm_atm_gdc2glo,gindex)
    call mct_permute(gindex,perm_atm_gdc2glo,lsize)
    call mct_gsMap_init(gsMap_atm_gdc2glo, gindex, mpicom, comp_id, lsize, gsize )
    deallocate(gindex)

    ! Diagnostic output

    if (masterproc) then
       write (6,*)' Atm Grid Characteristics'
       write (6,*)'   longitude points          = ',ani
       write (6,*)'   latitude points           = ',anj
       write (6,*)'   total number of gridcells = ',anumg
       write (6,*)' Surface Grid Characteristics'
       write (6,*)'   longitude points          = ',lni
       write (6,*)'   latitude points           = ',lnj
       write (6,*)'   total number of gridcells = ',numg
       write (6,*)'   total number of landunits = ',numl
       write (6,*)'   total number of columns   = ',numc
       write (6,*)'   total number of pfts      = ',nump
       write (6,*)' Decomposition Characteristics'
       write (6,*)'   clumps per process        = ',clump_pproc
       write (6,*)' gsMap Characteristics'
       write (6,*) '  lnd gsmap glo num of segs = ',mct_gsMap_ngseg(gsMap_lnd_gdc2glo)
       write (6,*) '  atm gsmap glo num of segs = ',mct_gsMap_ngseg(gsMap_atm_gdc2glo)
       write (6,*)
    end if

    ! Write out clump and proc info, one pe at a time, 
    ! barrier to control pes overwriting each other on stdout

#ifndef UNICOSMP
    call shr_sys_flush(6)
#endif
    call mpi_barrier(mpicom,ier)
     npmin = 0
     npmax = npes-1
     npint = 1
     if (dbug == 0) then
        npmax = 0
     elseif (dbug == 1) then
        npmax = min(npes-1,4)
     elseif (dbug == 2) then
        npint = npes/8
     endif
     do np = npmin,npmax,npint
       pid = np
       if (dbug == 1) then
          if (np == 2) pid=npes/2-1
          if (np == 3) pid=npes-2
          if (np == 4) pid=npes-1
       endif
       pid = max(pid,0)
       pid = min(pid,npes-1)

       if (iam == pid) then
          write(6,*)
          write(6,*)'proc= ',pid,' beg atmcell = ',procinfo%abegg, &
               ' end atmcell = ',procinfo%aendg,                   &
               ' total atmcells per proc = ',procinfo%aendg-procinfo%abegg+1
          write(6,*)'proc= ',pid,' beg gridcell= ',procinfo%begg, &
               ' end gridcell= ',procinfo%endg,                   &
               ' total gridcells per proc= ',procinfo%ncells
          write(6,*)'proc= ',pid,' beg landunit= ',procinfo%begl, &
               ' end landunit= ',procinfo%endl,                   &
               ' total landunits per proc= ',procinfo%nlunits
          write(6,*)'proc= ',pid,' beg column  = ',procinfo%begc, &
               ' end column  = ',procinfo%endc,                   &
               ' total columns per proc  = ',procinfo%ncols
          write(6,*)'proc= ',pid,' beg pft     = ',procinfo%begp, &
               ' end pft     = ',procinfo%endp,                   &
               ' total pfts per proc     = ',procinfo%npfts
          write(6,*)'proc= ',pid,' lnd ngseg   = ',mct_gsMap_ngseg(gsMap_lnd_gdc2glo), &
               ' lnd nlseg   = ',mct_gsMap_nlseg(gsMap_lnd_gdc2glo,iam)
          write(6,*)'proc= ',pid,' atm ngseg   = ',mct_gsMap_ngseg(gsMap_atm_gdc2glo), &
               ' atm nlseg   = ',mct_gsMap_nlseg(gsMap_atm_gdc2glo,iam)
          write(6,*)'proc= ',pid,' nclumps = ',procinfo%nclumps

          clmin = 1
          clmax = procinfo%nclumps
          if (dbug == 1) then
            clmax = 1
          elseif (dbug == 0) then
            clmax = -1
          endif
          do n = clmin,clmax
             cid = procinfo%cid(n)
             write(6,*)'proc= ',pid,' clump no = ',n, &
                  ' clump id= ',procinfo%cid(n),    &
                  ' beg gridcell= ',clumps(cid)%begg, &
                  ' end gridcell= ',clumps(cid)%endg, &
                  ' total gridcells per clump= ',clumps(cid)%ncells
             write(6,*)'proc= ',pid,' clump no = ',n, &
                  ' clump id= ',procinfo%cid(n),    &
                  ' beg landunit= ',clumps(cid)%begl, &
                  ' end landunit= ',clumps(cid)%endl, &
                  ' total landunits per clump = ',clumps(cid)%nlunits
             write(6,*)'proc= ',pid,' clump no = ',n, &
                  ' clump id= ',procinfo%cid(n),    &
                  ' beg column  = ',clumps(cid)%begc, &
                  ' end column  = ',clumps(cid)%endc, &
                  ' total columns per clump  = ',clumps(cid)%ncols
             write(6,*)'proc= ',pid,' clump no = ',n, &
                  ' clump id= ',procinfo%cid(n),    &
                  ' beg pft     = ',clumps(cid)%begp, &
                  ' end pft     = ',clumps(cid)%endp, &
                  ' total pfts per clump     = ',clumps(cid)%npfts
          end do
       end if
#ifndef UNICOSMP
       call shr_sys_flush(6)
#endif
       call mpi_barrier(mpicom,ier)
    end do
#ifndef UNICOSMP
    call shr_sys_flush(6)
#endif

  end subroutine decomp_init
#endif
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: decomp_lnd_init
!
! !INTERFACE:
  subroutine decomp_lnd_init(ans,ani,anj,lns,lni,lnj)
!
! !DESCRIPTION:
! This subroutine initializes the land surface decomposition into a clump
! data structure.  This assumes each pe has the same number of clumps
! set by clump_pproc
!
! !USES:
    use domainMod , only : gatm
!
! !ARGUMENTS:
    implicit none
    integer , intent(in) :: lns,lni,lnj ! land domain global size
    integer , intent(in) :: ans,ani,anj ! atm domain global size
!
! !LOCAL VARIABLES:
    integer :: lg,ln,li,lj        ! indices
    integer :: ag,an,ai,aj        ! indices
    integer :: anumg                  ! atm num gridcells
    integer :: cid,pid                ! indices
    integer, pointer :: lcid(:)       ! temporary for setting adecomp
    integer :: n,m,np                 ! indices
    integer :: ier                    ! error code
    integer :: cnt                    ! local counter
    integer :: beg,end,lsize,gsize    ! used for gsmap init
    integer, pointer :: gindex(:)     ! global index for gsmap init

    integer, pointer :: lncnt(:)      ! lnd cell count per atm cell
    integer, pointer :: lnoff(:)      ! atm cell offset in lnmap
    integer, pointer :: lnmap(:)      ! map from atm cell to lnd cells
    integer, pointer :: lglo2gsn(:)   ! map from glo 2 gsn temporary
    integer :: lnidx

! !CALLED FROM:
! subroutine initialize
!
! !REVISION HISTORY:
! 2002.09.11  Forrest Hoffman  Creation.
! 2005.12.15  T Craig  Updated for finemesh
! 2006.08.18  P Worley Performance optimizations
!
!EOP
!------------------------------------------------------------------------------

    allocate(lncnt(ans),lnoff(ans),lnmap(lns))

    lncnt = 0
    do ln = 1,lns
       an = gatm(ln)
       if ((an > 0) .and. (an .le. ans)) then
          lncnt(an) = lncnt(an) + 1
       endif
    enddo

    lnoff(1) = 1
    do an = 2,ans
       lnoff(an) = lnoff(an-1) + lncnt(an-1)
    enddo

    lncnt = 0
    lnmap = -1
    do ln = 1,lns
       an = gatm(ln)
       if ((an > 0) .and. (an .le. ans)) then
         lnmap(lnoff(an)+lncnt(an)) = ln
         lncnt(an) = lncnt(an) + 1
       endif
    enddo

    !--- assign gridcells to clumps (and thus pes) ---
    allocate(lcid(lns))
    lcid = 0
    numg = 0
    do anumg = 1,numa
          an = adecomp%gdc2glo(anumg)
          cid = acid(an)

          cnt = 0
          do lnidx = 0,lncnt(an)-1
             ln = lnmap(lnoff(an)+lnidx)          
             cnt = cnt + 1
             lcid(ln) = cid

             !--- overall total ---
             numg = numg + 1

             !--- give gridcell to cid ---
             !--- increment the beg and end indices ---
             clumps(cid)%ncells  = clumps(cid)%ncells  + 1

             do m = 1,nclumps
                if ((clumps(m)%owner >  clumps(cid)%owner) .or. &
                    (clumps(m)%owner == clumps(cid)%owner .and. m > cid)) then
                   clumps(m)%begg = clumps(m)%begg + 1
                endif

                if ((clumps(m)%owner >  clumps(cid)%owner) .or. &
                    (clumps(m)%owner == clumps(cid)%owner .and. m >= cid)) then
                   clumps(m)%endg = clumps(m)%endg + 1
                endif
             enddo

             !--- give gridcell to the proc that owns the cid ---
             !--- increment the beg and end indices ---
             if (iam == clumps(cid)%owner) then
                procinfo%ncells  = procinfo%ncells  + 1
             endif

             if (iam >  clumps(cid)%owner) then
                procinfo%begg = procinfo%begg + 1
             endif

             if (iam >= clumps(cid)%owner) then
                procinfo%endg = procinfo%endg + 1
             endif
          enddo
          !--- check that atm cell has at least 1 lnd grid cell
          if (cnt < 1) then
             write (6,*) 'decomp_lnd_init(): map overlap error at ',an,cnt
             call endrun()
          endif
    enddo

    allocate(ldecomp%gdc2gsn(numg), ldecomp%gsn2gdc(numg), &
             ldecomp%gdc2glo(numg), ldecomp%glo2gdc(lni*lnj), &
             ldecomp%gdc2i  (numg), ldecomp%gdc2j  (numg), &
             stat=ier)
    if (ier /= 0) then
       write (6,*) 'decomp_lnd_init(): allocation error1 for ldecomp'
       call endrun()
    end if

    ldecomp%gdc2gsn(:)  = 0
    ldecomp%gsn2gdc(:)  = 0
    ldecomp%gdc2glo(:)  = 0
    ldecomp%glo2gdc(:)  = 0
    ldecomp%gdc2i(:)    = 0
    ldecomp%gdc2j(:)    = 0

    !--- temporaries for decomp mappings
    allocate(lglo2gsn(lni*lnj),stat=ier)
    if (ier /= 0) then
       write (6,*) 'decomp_lnd_init(): allocation error for lglo2gsn'
       call endrun()
    end if
    lglo2gsn(:) = 0

    ! Set ldecomp sn indexing based on cells to be used and i,j order
    lg  = 0
    do lj = 1,lnj
    do li = 1,lni
       ln = (lj-1)*lni + li
       if (lcid(ln) > 0) then
          lg = lg + 1
          lglo2gsn(ln) = lg
       endif
    enddo
    enddo

    ! Set ldecomp and adecomp data
    lg = 0
    do pid = 0,npes-1
    do cid = 1,nclumps
       if (clumps(cid)%owner == pid) then

          do lj = 1,lnj
          do li = 1,lni
             ln = (lj-1)*lni + li
             if (lcid(ln) == cid) then
                lg = lg + 1
                ldecomp%gdc2i(lg) = li
                ldecomp%gdc2j(lg) = lj
                ldecomp%gdc2gsn(lg) = lglo2gsn(ln)
                ldecomp%gdc2glo(lg) = ln
                ldecomp%gsn2gdc(lglo2gsn(ln)) = lg
                ldecomp%glo2gdc(ln) = lg
             endif
          enddo
          enddo
       endif
    enddo
    enddo

    deallocate(lglo2gsn)
    deallocate(lcid)
    deallocate(lncnt,lnoff,lnmap)

    ! set gsMap_lnd_gdc2glo, perm_lnd_gdc2glo
    call get_proc_bounds(beg, end)
    allocate(gindex(beg:end))
    do n = beg,end
       gindex(n) = ldecomp%gdc2glo(n)
    enddo
    lsize = end-beg+1
    gsize = lni * lnj
    allocate(perm_lnd_gdc2glo(lsize),stat=ier)
    call mct_indexset(perm_lnd_gdc2glo)
    call mct_indexsort(lsize,perm_lnd_gdc2glo,gindex)
    call mct_permute(gindex,perm_lnd_gdc2glo,lsize)
    call mct_gsMap_init(gsMap_lnd_gdc2glo, gindex, mpicom, comp_id, lsize, gsize )
    deallocate(gindex)

    ! Diagnostic output

    if (masterproc) then
       write (6,*)' Atm Grid Characteristics'
       write (6,*)'   longitude points          = ',ani
       write (6,*)'   latitude points           = ',anj
       write (6,*)'   total number of gridcells = ',numa
       write (6,*)' Surface Grid Characteristics'
       write (6,*)'   longitude points          = ',lni
       write (6,*)'   latitude points           = ',lnj
       write (6,*)'   total number of gridcells = ',numg
       write (6,*)' Decomposition Characteristics'
       write (6,*)'   clumps per process        = ',clump_pproc
       write (6,*)' gsMap Characteristics'
       write (6,*) '  lnd gsmap glo num of segs = ',mct_gsMap_ngseg(gsMap_lnd_gdc2glo)
       write (6,*) '  atm gsmap glo num of segs = ',mct_gsMap_ngseg(gsMap_atm_gdc2glo)
       write (6,*)
    end if

    call shr_sys_flush(6)

  end subroutine decomp_lnd_init

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: decomp_glcp_init
!
! !INTERFACE:
  subroutine decomp_glcp_init(ans,ani,anj,lns,lni,lnj)
!
! !DESCRIPTION:
! This subroutine initializes the land surface decomposition into a clump
! data structure.  This assumes each pe has the same number of clumps
! set by clump_pproc
!
! !USES:
    use spmdMod
    use domainMod , only : gatm
    use subgridMod, only : subgrid_get_gcellinfo
!
! !ARGUMENTS:
    implicit none
    integer , intent(in) :: lns,lni,lnj ! land domain global size
    integer , intent(in) :: ans,ani,anj ! atm domain global size
!
! !LOCAL VARIABLES:
    integer :: lg,ln,li,lj        ! indices
    integer :: ag,an,ai,aj        ! indices
    integer :: abegg,aendg,anumg  ! atm num gridcells
    integer :: begg,endg          ! lnd num gridcells
    integer :: cid,pid                ! indices
    integer :: n,m,np                 ! indices
    integer :: icells, ilunits, icols, ipfts  ! temporaries
    integer :: ier                    ! error code
    integer :: cnt                    ! local counter
    integer, parameter :: dbug=3      ! 0 = min, 1=normal, 2=much, 3=max
    integer :: npmin,npmax,npint      ! do loop values for printing
    integer :: clmin,clmax,clint      ! do loop values for printing
    integer :: beg,end,lsize,gsize    ! used for gsmap init
    integer, pointer :: gindex(:)     ! global index for gsmap init

    integer, pointer :: lncnt(:)      ! lnd cell count per atm cell
    integer, pointer :: lnoff(:)      ! atm cell offset in lnmap
    integer, pointer :: lnmap(:)      ! map from atm cell to lnd cells
    integer, allocatable :: allvecg(:,:)  ! temporary vector "global"
    integer, allocatable :: allvecl(:,:)  ! temporary vector "local"
    integer, pointer :: lglo2gsn(:)   ! map from glo 2 gsn temporary
    integer :: lnidx

! !CALLED FROM:
! subroutine initialize
!
! !REVISION HISTORY:
! 2002.09.11  Forrest Hoffman  Creation.
! 2005.12.15  T Craig  Updated for finemesh
! 2006.08.18  P Worley Performance optimizations
!
!EOP
!------------------------------------------------------------------------------

    allocate(lncnt(ans),lnoff(ans),lnmap(lns))

    lncnt = 0
    do ln = 1,lns
       an = gatm(ln)
       if ((an > 0) .and. (an .le. ans)) then
          lncnt(an) = lncnt(an) + 1
       endif
    enddo

    lnoff(1) = 1
    do an = 2,ans
       lnoff(an) = lnoff(an-1) + lncnt(an-1)
    enddo

    lncnt = 0
    lnmap = -1
    do ln = 1,lns
       an = gatm(ln)
       if ((an > 0) .and. (an .le. ans)) then
         lnmap(lnoff(an)+lncnt(an)) = ln
         lncnt(an) = lncnt(an) + 1
       endif
    enddo

    !--- assign gridcells to clumps (and thus pes) ---
    call get_proc_bounds_atm(abegg, aendg)
    call get_proc_bounds(begg, endg)
    allocate(allvecg(nclumps,4),allvecl(nclumps,4))   ! 3 = gcells,lunit,cols,pfts
    allvecg  = 0
    allvecl = 0
    do anumg = abegg,aendg
          an = adecomp%gdc2glo(anumg)
          cid = acid(an)
          do lnidx = 0,lncnt(an)-1
!             ln = lnmap(lnoff(an)+lnidx)
             ln = ldecomp%glo2gdc(lnmap(lnoff(an)+lnidx))
             call subgrid_get_gcellinfo (ln, nlunits=ilunits, &
                                  ncols=icols, npfts=ipfts)
             allvecl(cid,1) = allvecl(cid,1) + 1
             allvecl(cid,2) = allvecl(cid,2) + ilunits
             allvecl(cid,3) = allvecl(cid,3) + icols
             allvecl(cid,4) = allvecl(cid,4) + ipfts
          enddo
    enddo
    call mpi_allreduce(allvecl,allvecg,size(allvecg),MPI_INTEGER,MPI_SUM,mpicom,ier)

    numg  = 0
    numl  = 0
    numc  = 0
    nump  = 0
    do cid = 1,nclumps
             icells  = allvecg(cid,1)
             ilunits = allvecg(cid,2)
             icols   = allvecg(cid,3)
             ipfts   = allvecg(cid,4)

             !--- overall total ---
             numg = numg + icells
             numl = numl + ilunits
             numc = numc + icols
             nump = nump + ipfts

             !--- give gridcell to cid ---
             !--- increment the beg and end indices ---
!             clumps(cid)%ncells  = clumps(cid)%ncells  + icells
             clumps(cid)%nlunits = clumps(cid)%nlunits + ilunits
             clumps(cid)%ncols   = clumps(cid)%ncols   + icols
             clumps(cid)%npfts   = clumps(cid)%npfts   + ipfts

             do m = 1,nclumps
                if ((clumps(m)%owner >  clumps(cid)%owner) .or. &
                    (clumps(m)%owner == clumps(cid)%owner .and. m > cid)) then
!                   clumps(m)%begg = clumps(m)%begg + icells
                   clumps(m)%begl = clumps(m)%begl + ilunits
                   clumps(m)%begc = clumps(m)%begc + icols
                   clumps(m)%begp = clumps(m)%begp + ipfts
                endif

                if ((clumps(m)%owner >  clumps(cid)%owner) .or. &
                    (clumps(m)%owner == clumps(cid)%owner .and. m >= cid)) then
!                   clumps(m)%endg = clumps(m)%endg + icells
                   clumps(m)%endl = clumps(m)%endl + ilunits
                   clumps(m)%endc = clumps(m)%endc + icols
                   clumps(m)%endp = clumps(m)%endp + ipfts
                endif
             enddo

             !--- give gridcell to the proc that owns the cid ---
             !--- increment the beg and end indices ---
             if (iam == clumps(cid)%owner) then
!                procinfo%ncells  = procinfo%ncells  + icells
                procinfo%nlunits = procinfo%nlunits + ilunits
                procinfo%ncols   = procinfo%ncols   + icols
                procinfo%npfts   = procinfo%npfts   + ipfts
             endif

             if (iam >  clumps(cid)%owner) then
!                procinfo%begg = procinfo%begg + icells
                procinfo%begl = procinfo%begl + ilunits
                procinfo%begc = procinfo%begc + icols
                procinfo%begp = procinfo%begp + ipfts
             endif

             if (iam >= clumps(cid)%owner) then
!                procinfo%endg = procinfo%endg + icells
                procinfo%endl = procinfo%endl + ilunits
                procinfo%endc = procinfo%endc + icols
                procinfo%endp = procinfo%endp + ipfts
             endif
    enddo

#if (1 == 0)
    numg  = 0
    numl  = 0
    numc  = 0
    nump  = 0
    do anumg = 1,numa
          an = adecomp%gdc2glo(anumg)
          cid = acid(an)

          cnt = 0
          do lnidx = 0,lncnt(an)-1
             ln = lnmap(lnoff(an)+lnidx)          
             cnt = cnt + 1
!             call subgrid_get_gcellinfo (ln, wtxy, nlunits=ilunits, &
             call subgrid_get_gcellinfo (ln, nlunits=ilunits, &
                                  ncols=icols, npfts=ipfts, global=.true.)

             !--- overall total ---
             numg = numg + 1
             numl = numl + ilunits
             numc = numc + icols
             nump = nump + ipfts

             !--- give gridcell to cid ---
             !--- increment the beg and end indices ---
!             clumps(cid)%ncells  = clumps(cid)%ncells  + 1
             clumps(cid)%nlunits = clumps(cid)%nlunits + ilunits
             clumps(cid)%ncols   = clumps(cid)%ncols   + icols
             clumps(cid)%npfts   = clumps(cid)%npfts   + ipfts

             do m = 1,nclumps
                if ((clumps(m)%owner >  clumps(cid)%owner) .or. &
                    (clumps(m)%owner == clumps(cid)%owner .and. m > cid)) then
!                   clumps(m)%begg = clumps(m)%begg + 1
                   clumps(m)%begl = clumps(m)%begl + ilunits
                   clumps(m)%begc = clumps(m)%begc + icols
                   clumps(m)%begp = clumps(m)%begp + ipfts
                endif

                if ((clumps(m)%owner >  clumps(cid)%owner) .or. &
                    (clumps(m)%owner == clumps(cid)%owner .and. m >= cid)) then
!                   clumps(m)%endg = clumps(m)%endg + 1
                   clumps(m)%endl = clumps(m)%endl + ilunits
                   clumps(m)%endc = clumps(m)%endc + icols
                   clumps(m)%endp = clumps(m)%endp + ipfts
                endif
             enddo

             !--- give gridcell to the proc that owns the cid ---
             !--- increment the beg and end indices ---
             if (iam == clumps(cid)%owner) then
!                procinfo%ncells  = procinfo%ncells  + 1
                procinfo%nlunits = procinfo%nlunits + ilunits
                procinfo%ncols   = procinfo%ncols   + icols
                procinfo%npfts   = procinfo%npfts   + ipfts
             endif

             if (iam >  clumps(cid)%owner) then
!                procinfo%begg = procinfo%begg + 1
                procinfo%begl = procinfo%begl + ilunits
                procinfo%begc = procinfo%begc + icols
                procinfo%begp = procinfo%begp + ipfts
             endif

             if (iam >= clumps(cid)%owner) then
!                procinfo%endg = procinfo%endg + 1
                procinfo%endl = procinfo%endl + ilunits
                procinfo%endc = procinfo%endc + icols
                procinfo%endp = procinfo%endp + ipfts
             endif
          enddo
          !--- check that atm cell has at least 1 lnd grid cell
          if (cnt < 1) then
             write (6,*) 'decomp_glcp_init(): map overlap error at ',an,cnt
             call endrun()
          endif
    enddo
#endif

    do n = 1,nclumps
       if (clumps(n)%ncells  /= allvecg(n,1) .or. &
           clumps(n)%nlunits /= allvecg(n,2) .or. &
           clumps(n)%ncols   /= allvecg(n,3) .or. &
           clumps(n)%npfts   /= allvecg(n,4)) then
          write(6,*) 'decomp_glcp_init(): allvecg error ncells ',iam,n,clumps(n)%ncells ,allvecg(n,1)
          write(6,*) 'decomp_glcp_init(): allvecg error lunits ',iam,n,clumps(n)%nlunits,allvecg(n,2)
          write(6,*) 'decomp_glcp_init(): allvecg error ncols  ',iam,n,clumps(n)%ncols  ,allvecg(n,3)
          write(6,*) 'decomp_glcp_init(): allvecg error pfts   ',iam,n,clumps(n)%npfts  ,allvecg(n,4)
          call shr_sys_flush(6)
          call shr_sys_abort()
       endif
    enddo

    deallocate(allvecg,allvecl)
    deallocate(acid)
    deallocate(lncnt,lnoff,lnmap)

    ! Diagnostic output

    if (masterproc) then
       write (6,*)' Atm Grid Characteristics'
       write (6,*)'   longitude points          = ',ani
       write (6,*)'   latitude points           = ',anj
       write (6,*)'   total number of gridcells = ',numa
       write (6,*)' Surface Grid Characteristics'
       write (6,*)'   longitude points          = ',lni
       write (6,*)'   latitude points           = ',lnj
       write (6,*)'   total number of gridcells = ',numg
       write (6,*)'   total number of landunits = ',numl
       write (6,*)'   total number of columns   = ',numc
       write (6,*)'   total number of pfts      = ',nump
       write (6,*)' Decomposition Characteristics'
       write (6,*)'   clumps per process        = ',clump_pproc
       write (6,*)' gsMap Characteristics'
       write (6,*) '  lnd gsmap glo num of segs = ',mct_gsMap_ngseg(gsMap_lnd_gdc2glo)
       write (6,*) '  atm gsmap glo num of segs = ',mct_gsMap_ngseg(gsMap_atm_gdc2glo)
       write (6,*)
    end if

    ! Write out clump and proc info, one pe at a time, 
    ! barrier to control pes overwriting each other on stdout

    call shr_sys_flush(6)
    call mpi_barrier(mpicom,ier)
     npmin = 0
     npmax = npes-1
     npint = 1
     if (dbug == 0) then
        npmax = 0
     elseif (dbug == 1) then
        npmax = min(npes-1,4)
     elseif (dbug == 2) then
        npint = npes/8
     endif
     do np = npmin,npmax,npint
       pid = np
       if (dbug == 1) then
          if (np == 2) pid=npes/2-1
          if (np == 3) pid=npes-2
          if (np == 4) pid=npes-1
       endif
       pid = max(pid,0)
       pid = min(pid,npes-1)

       if (iam == pid) then
          write(6,*)
          write(6,*)'proc= ',pid,' beg atmcell = ',procinfo%abegg, &
               ' end atmcell = ',procinfo%aendg,                   &
               ' total atmcells per proc = ',procinfo%aendg-procinfo%abegg+1
          write(6,*)'proc= ',pid,' beg gridcell= ',procinfo%begg, &
               ' end gridcell= ',procinfo%endg,                   &
               ' total gridcells per proc= ',procinfo%ncells
          write(6,*)'proc= ',pid,' beg landunit= ',procinfo%begl, &
               ' end landunit= ',procinfo%endl,                   &
               ' total landunits per proc= ',procinfo%nlunits
          write(6,*)'proc= ',pid,' beg column  = ',procinfo%begc, &
               ' end column  = ',procinfo%endc,                   &
               ' total columns per proc  = ',procinfo%ncols
          write(6,*)'proc= ',pid,' beg pft     = ',procinfo%begp, &
               ' end pft     = ',procinfo%endp,                   &
               ' total pfts per proc     = ',procinfo%npfts
          write(6,*)'proc= ',pid,' lnd ngseg   = ',mct_gsMap_ngseg(gsMap_lnd_gdc2glo), &
               ' lnd nlseg   = ',mct_gsMap_nlseg(gsMap_lnd_gdc2glo,iam)
          write(6,*)'proc= ',pid,' atm ngseg   = ',mct_gsMap_ngseg(gsMap_atm_gdc2glo), &
               ' atm nlseg   = ',mct_gsMap_nlseg(gsMap_atm_gdc2glo,iam)
          write(6,*)'proc= ',pid,' nclumps = ',procinfo%nclumps

          clmin = 1
          clmax = procinfo%nclumps
          if (dbug == 1) then
            clmax = 1
          elseif (dbug == 0) then
            clmax = -1
          endif
          do n = clmin,clmax
             cid = procinfo%cid(n)
             write(6,*)'proc= ',pid,' clump no = ',n, &
                  ' clump id= ',procinfo%cid(n),    &
                  ' beg gridcell= ',clumps(cid)%begg, &
                  ' end gridcell= ',clumps(cid)%endg, &
                  ' total gridcells per clump= ',clumps(cid)%ncells
             write(6,*)'proc= ',pid,' clump no = ',n, &
                  ' clump id= ',procinfo%cid(n),    &
                  ' beg landunit= ',clumps(cid)%begl, &
                  ' end landunit= ',clumps(cid)%endl, &
                  ' total landunits per clump = ',clumps(cid)%nlunits
             write(6,*)'proc= ',pid,' clump no = ',n, &
                  ' clump id= ',procinfo%cid(n),    &
                  ' beg column  = ',clumps(cid)%begc, &
                  ' end column  = ',clumps(cid)%endc, &
                  ' total columns per clump  = ',clumps(cid)%ncols
             write(6,*)'proc= ',pid,' clump no = ',n, &
                  ' clump id= ',procinfo%cid(n),    &
                  ' beg pft     = ',clumps(cid)%begp, &
                  ' end pft     = ',clumps(cid)%endp, &
                  ' total pfts per clump     = ',clumps(cid)%npfts
          end do
       end if
#ifndef UNICOSMP
       call shr_sys_flush(6)
#endif
       call mpi_barrier(mpicom,ier)
    end do
    call shr_sys_flush(6)

  end subroutine decomp_glcp_init

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: decomp_atm_init
!
! !INTERFACE:
  subroutine decomp_atm_init(alatlon,amask)
!
! !DESCRIPTION:
! This subroutine initializes the land surface decomposition into a clump
! data structure.  This assumes each pe has the same number of clumps
! set by clump_pproc
!
! !USES:
    use clm_varctl, only : nsegspc
    use domainMod , only : latlon_type
!
! !ARGUMENTS:
    implicit none
    type(latlon_type),intent(in) :: alatlon
    integer          ,intent(in) :: amask(:)
!
! !LOCAL VARIABLES:
    integer :: ani,anj                ! atm domain global size
    integer :: ans,ag,an,ai,aj        ! indices
    integer :: anumg                  ! atm num gridcells
    integer :: anumg_tot              ! precompute of anumg
    real(r8):: rnsegspc               ! real value associated with nsegspc
    integer :: cid,pid                ! indices
    integer :: n,m,np                 ! indices
    integer :: ier                    ! error code
    integer, parameter :: dbug=3      ! 0 = min, 1=normal, 2=much, 3=max
    integer :: npmin,npmax,npint      ! do loop values for printing
    integer :: clmin,clmax,clint      ! do loop values for printing
    integer :: beg,end,lsize,gsize    ! used for gsmap init
    integer, pointer :: gindex(:)     ! global index for gsmap init

    integer, pointer :: aglo2gsn(:)   ! map from glo 2 gsn temporary

! !CALLED FROM:
! subroutine initialize
!
! !REVISION HISTORY:
! 2002.09.11  Forrest Hoffman  Creation.
! 2005.12.15  T Craig  Updated for finemesh
! 2006.08.18  P Worley Performance optimizations
! 2007.01.24  T Craig  Created decomp_atm_init from decomp_init
!
!EOP
!------------------------------------------------------------------------------

    ani = alatlon%ni
    anj = alatlon%nj
    ans = alatlon%ns

    !--- set and verify nclumps ---
    if (clump_pproc > 0) then
       nclumps = clump_pproc * npes
       if (nclumps < npes) then
          write (6,*) 'decomp_atm_init(): Number of gridcell clumps= ',nclumps, &
               ' is less than the number of processes = ', npes
          call endrun()
       end if
    else
       write(6,*)'clump_pproc= ',clump_pproc,'  must be greater than 0'
       call endrun()
    end if

    !--- allocate and initialize procinfo and clumps ---
    !--- beg and end indices initialized for simple addition of cells later ---

    allocate(procinfo%cid(clump_pproc), stat=ier)
    if (ier /= 0) then
       write (6,*) 'decomp_atm_init(): allocation error for procinfo%cid'
       call endrun()
    endif

    procinfo%nclumps = clump_pproc
    procinfo%cid(:)  = -1
    procinfo%ncells  = 0
    procinfo%nlunits = 0
    procinfo%ncols   = 0
    procinfo%npfts   = 0
    procinfo%begg    = 1
    procinfo%begl    = 1
    procinfo%begc    = 1
    procinfo%begp    = 1
    procinfo%endg    = 0
    procinfo%endl    = 0
    procinfo%endc    = 0
    procinfo%endp    = 0
    procinfo%abegg   = 1
    procinfo%aendg   = 0

    allocate(clumps(nclumps), stat=ier)
    if (ier /= 0) then
       write (6,*) 'decomp_atm_init(): allocation error for clumps'
       call endrun()
    end if
    clumps(:)%owner   = -1
    clumps(:)%ncells  = 0
    clumps(:)%nlunits = 0
    clumps(:)%ncols   = 0
    clumps(:)%npfts   = 0
    clumps(:)%begg    = 1
    clumps(:)%begl    = 1
    clumps(:)%begc    = 1
    clumps(:)%begp    = 1
    clumps(:)%endg    = 0
    clumps(:)%endl    = 0
    clumps(:)%endc    = 0
    clumps(:)%endp    = 0

    !--- assign clumps to proc round robin ---
    cid = 0
    do n = 1,nclumps
       pid = mod(n-1,npes)
       if (pid < 0 .or. pid > npes-1) then
          write (6,*) 'decomp_atm_init(): round robin pid error ',n,pid,npes
          call endrun()
       endif
       clumps(n)%owner = pid
       if (iam == pid) then
          cid = cid + 1
          if (cid < 1 .or. cid > clump_pproc) then
             write (6,*) 'decomp_atm_init(): round robin pid error ',n,pid,npes
             call endrun()
          endif
          procinfo%cid(cid) = n
       endif
    enddo

    !--- count total atm gridcells
    anumg_tot = 0
    do an = 1,ans
       if (amask(an) == 1) then
          anumg_tot  = anumg_tot  + 1
       endif
    enddo
    numa = anumg_tot

    rnsegspc = min(float(nsegspc),float(anumg_tot)/float(nclumps))
    if (masterproc) write(6,*) 'precompute total anumg ',anumg_tot,numa,rnsegspc

    !--- assign gridcells to clumps (and thus pes) ---
    allocate(acid(ans))
    acid = 0
    anumg = 0
    do an = 1,ans
       if (amask(an) == 1) then
          anumg  = anumg  + 1

#if (1 == 0) 
          !--- find clump with fewest pfts ---
          cid = 1
          do n = 2,nclumps
             if (clumps(n)%npfts < clumps(cid)%npfts) then
                cid = n
             endif
          enddo
#endif
#if (1 == 0)
          !--- give to clumps in simple round robin
          cid = mod((anumg-1),nclumps) + 1
#endif
#if (1 == 1)
          !--- give to clumps in order based on nsegspc
          cid = int(rnsegspc*float(nclumps*(anumg-1))/float(anumg_tot))
          cid = mod(cid,nclumps) + 1
#endif
          acid(an) = cid

          !--- give atm cell to pe that owns cid ---
          if (iam >  clumps(cid)%owner) then
             procinfo%abegg = procinfo%abegg + 1
          endif
          if (iam >= clumps(cid)%owner) then
             procinfo%aendg = procinfo%aendg + 1
          endif

       end if
    enddo

    ! Error check on total number of gridcells

    if (anumg /= numa) then
       write (6,*) 'decomp_atm_init(): Number of atm gridcells inconsistent',anumg,numa
       call endrun()
    end if

    if (npes > anumg) then
       write (6,*) 'decomp_atm_init(): Number of processes exceeds number ', &
            'of atm grid cells',npes,anumg
       call endrun()
    end if

    ! Allocate dynamic memory for adecomp, ldecomp derived type

    allocate(adecomp%gdc2gsn(anumg), adecomp%gsn2gdc(anumg), &
             adecomp%gdc2glo(anumg), adecomp%glo2gdc(ani*anj), &
             adecomp%gdc2i  (anumg), adecomp%gdc2j  (anumg), &
             stat=ier)
    if (ier /= 0) then
       write (6,*) 'decomp_atm_init(): allocation error1 for adecomp'
       call endrun()
    end if

    adecomp%gdc2gsn(:)  = 0
    adecomp%gsn2gdc(:)  = 0
    adecomp%gdc2glo(:)  = 0
    adecomp%glo2gdc(:)  = 0
    adecomp%gdc2i(:)    = 0
    adecomp%gdc2j(:)    = 0

    !--- temporaries for decomp mappings
    allocate(aglo2gsn(ani*anj),stat=ier)
    if (ier /= 0) then
       write (6,*) 'decomp_atm_init(): allocation error for aglo2gsn'
       call endrun()
    end if
    aglo2gsn(:) = 0

    ag = 0
    do aj = 1,anj
    do ai = 1,ani
       an = (aj-1)*ani + ai
       if (acid(an) > 0) then
          ag  = ag  + 1
          aglo2gsn(an) = ag
       endif
    enddo
    enddo

    ! Set ldecomp and adecomp data
    ag = 0
    do pid = 0,npes-1
    do cid = 1,nclumps
       if (clumps(cid)%owner == pid) then

          do aj = 1,anj
          do ai = 1,ani
             an = (aj-1)*ani + ai
             if (acid(an) == cid) then
                ag = ag + 1
                adecomp%gdc2i(ag) = ai
                adecomp%gdc2j(ag) = aj
                adecomp%gdc2gsn(ag) = aglo2gsn(an)
                adecomp%gdc2glo(ag) = an
                adecomp%gsn2gdc(aglo2gsn(an)) = ag
                adecomp%glo2gdc(an) = ag
             endif
          enddo
          enddo

       endif
    enddo
    enddo

    deallocate(aglo2gsn)

    ! set gsMap_atm_gdc2glo, perm_atm_gdc2glo
    call get_proc_bounds_atm(beg, end)
    allocate(gindex(beg:end))
    do n = beg,end
       gindex(n) = adecomp%gdc2glo(n)
    enddo
    lsize = end-beg+1
    gsize = ani * anj
    allocate(perm_atm_gdc2glo(lsize),stat=ier)
    call mct_indexset(perm_atm_gdc2glo)
    call mct_indexsort(lsize,perm_atm_gdc2glo,gindex)
    call mct_permute(gindex,perm_atm_gdc2glo,lsize)
    call mct_gsMap_init(gsMap_atm_gdc2glo, gindex, mpicom, comp_id, lsize, gsize )
    deallocate(gindex)

    ! Diagnostic output

    if (masterproc) then
       write (6,*)' Atm Grid Characteristics'
       write (6,*)'   longitude points          = ',ani
       write (6,*)'   latitude points           = ',anj
       write (6,*)'   total number of gridcells = ',anumg
       write (6,*)' Decomposition Characteristics'
       write (6,*)'   clumps per process        = ',clump_pproc
       write (6,*)' gsMap Characteristics'
       write (6,*) '  atm gsmap glo num of segs = ',mct_gsMap_ngseg(gsMap_atm_gdc2glo)
       write (6,*)
    end if

    ! Write out clump and proc info, one pe at a time, 
    ! barrier to control pes overwriting each other on stdout

    call shr_sys_flush(6)
    call mpi_barrier(mpicom,ier)
     npmin = 0
     npmax = npes-1
     npint = 1
     if (dbug == 0) then
        npmax = 0
     elseif (dbug == 1) then
        npmax = min(npes-1,4)
     elseif (dbug == 2) then
        npint = npes/8
     endif
     do np = npmin,npmax,npint
       pid = np
       if (dbug == 1) then
          if (np == 2) pid=npes/2-1
          if (np == 3) pid=npes-2
          if (np == 4) pid=npes-1
       endif
       pid = max(pid,0)
       pid = min(pid,npes-1)

       if (iam == pid) then
          write(6,*)
          write(6,*)'proc= ',pid,' beg atmcell = ',procinfo%abegg, &
               ' end atmcell = ',procinfo%aendg,                   &
               ' total atmcells per proc = ',procinfo%aendg-procinfo%abegg+1
          write(6,*)'proc= ',pid,' atm ngseg   = ',mct_gsMap_ngseg(gsMap_atm_gdc2glo), &
               ' atm nlseg   = ',mct_gsMap_nlseg(gsMap_atm_gdc2glo,iam)
          write(6,*)'proc= ',pid,' nclumps = ',procinfo%nclumps
       end if
       call shr_sys_flush(6)
       call mpi_barrier(mpicom,ier)
    end do
#ifndef UNICOSMP
    call shr_sys_flush(6)
#endif

  end subroutine decomp_atm_init

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: decomp_domg2l
!
! !INTERFACE:
  subroutine decomp_domg2l(gdomain,ldomain,decomp,nbeg,nend)
!
! !DESCRIPTION:
! This subroutine allocates and nans the domain type
!
! !USES:
    use domainMod, only : domain_type, domain_init
!
! !ARGUMENTS:
    implicit none
    type(domain_type) :: gdomain       ! global domain
    type(domain_type) :: ldomain       ! local domain
    type(decomp_type) :: decomp        ! decomp type
    integer           :: nbeg,nend     ! local beg/end
!
! !REVISION HISTORY:
!   Created by T Craig
!
!EOP
!
! LOCAL VARIABLES:
    integer ier,nl,ng
!
!------------------------------------------------------------------------------

    call domain_init(ldomain,gdomain%ni,gdomain%nj,nbeg,nend)

    ldomain%edges   = gdomain%edges
    do nl = nbeg,nend
       ng = decomp%gdc2glo(nl)
       ldomain%mask(nl) = gdomain%mask(ng)
       ldomain%frac(nl) = gdomain%frac(ng)
       ldomain%topo(nl) = gdomain%topo(ng)
       ldomain%latc(nl) = gdomain%latc(ng)
       ldomain%lonc(nl) = gdomain%lonc(ng)
       ldomain%area(nl) = gdomain%area(ng)
       ldomain%lats(nl) = gdomain%lats(ng)
       ldomain%latn(nl) = gdomain%latn(ng)
       ldomain%lonw(nl) = gdomain%lonw(ng)
       ldomain%lone(nl) = gdomain%lone(ng)

       ldomain%pftm(nl) = gdomain%pftm(ng)
       ldomain%nara(nl) = gdomain%nara(ng)
       ldomain%ntop(nl) = gdomain%ntop(ng)
!       ldomain%gatm(nl) = gdomain%gatm(ng)
    enddo

end subroutine decomp_domg2l
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

     cid  = procinfo%cid(n)
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
        begp = procinfo%begp
     endif
     if (present(endp)) then
        endp = procinfo%endp
     endif
     if (present(begc)) then
        begc = procinfo%begc
     endif
     if (present(endc)) then
        endc = procinfo%endc
     endif
     if (present(begl)) then
        begl = procinfo%begl
     endif
     if (present(endl)) then
        endl = procinfo%endl
     endif
     if (present(begg)) then
        begg = procinfo%begg
     endif
     if (present(endg)) then
        endg = procinfo%endg
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

   begg = procinfo%abegg
   endg = procinfo%aendg

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
   subroutine get_proc_global(ng, nl, nc, np)
!
! !DESCRIPTION:
! Return number of gridcells, landunits, columns, and pfts across all
! processes.
!
! !USES:
!
! !ARGUMENTS:
     implicit none
     integer, intent(out) :: ng  ! total number of gridcells
                                 ! across all processors
     integer, intent(out) :: nl  ! total number of landunits
                                 ! across all processors
     integer, intent(out) :: nc  ! total number of columns
                                 ! across all processors
     integer, intent(out) :: np  ! total number of pfts
                                 ! across all processors
! !REVISION HISTORY:
! 2003.09.12  Mariana Vertenstein  Creation.
!
!EOP
!------------------------------------------------------------------------------

     np = nump
     nc = numc
     nl = numl
     ng = numg

   end subroutine get_proc_global

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_proc_global_atm
!
! !INTERFACE:
   subroutine get_proc_global_atm(na)
!
! !DESCRIPTION:
! Return number of gridcells, landunits, columns, and pfts across all
! processes.
!
! !USES:
!
! !ARGUMENTS:
     implicit none
     integer, intent(out) :: na  ! total number of atm gridcells
                                 ! across all processors
! !REVISION HISTORY:
! 2003.09.12  Mariana Vertenstein  Creation.
!
!EOP
!------------------------------------------------------------------------------

     na = numa

   end subroutine get_proc_global_atm


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

     get_proc_clumps = procinfo%nclumps

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
   use subgridMod, only : subgrid_get_indexes
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
     do gdc = 1,numg
        gsn = ldecomp%gdc2gsn(gdc)
        call subgrid_get_indexes(gdc,gsn,type1d,dci,dcf,sni,snf)
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
   use subgridMod, only : subgrid_get_indexes
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

     do gdc = 1,numg
        gsn = ldecomp%gdc2gsn(gdc)
        call subgrid_get_indexes(gdc,gsn,type1d,dci,dcf,sni,snf)
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
   use subgridMod, only : subgrid_get_indexes
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

     do gdc = 1,numg
        gsn = ldecomp%gdc2gsn(gdc)
        call subgrid_get_indexes(gdc,gsn,type1d,dci,dcf,sni,snf)
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
   use clmtype, only : nameg, namel, namec, namep!
   use subgridMod, only : subgrid_get_indexes

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

     do gdc = 1,numg
        gsn = ldecomp%gdc2gsn(gdc)
        call subgrid_get_indexes(gdc,gsn,type1d,dci,dcf,sni,snf)
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
   use subgridMod, only : subgrid_get_indexes
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

     do gsn = 1,numg
        gdc = ldecomp%gsn2gdc(gsn)
        call subgrid_get_indexes(gdc,gsn,type1d,dci,dcf,sni,snf)
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
   use subgridMod, only : subgrid_get_indexes
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

     do gsn = 1,numg
        gdc = ldecomp%gsn2gdc(gsn)
        call subgrid_get_indexes(gdc,gsn,type1d,dci,dcf,sni,snf)
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
   use subgridMod, only : subgrid_get_indexes
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

     do gsn = 1,numg
        gdc = ldecomp%gsn2gdc(gsn)
        call subgrid_get_indexes(gdc,gsn,type1d,dci,dcf,sni,snf)
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
   use subgridMod, only : subgrid_get_indexes
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

     do gsn = 1,numg
        gdc = ldecomp%gsn2gdc(gsn)
        call subgrid_get_indexes(gdc,gsn,type1d,dci,dcf,sni,snf)
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

end module decompMod
