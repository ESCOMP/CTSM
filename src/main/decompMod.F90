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
  use spmdMod     , only : masterproc, iam, npes, mpicom
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
  public decomp_init             ! initializes land surface decomposition
                                 ! into clumps and processors
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
     integer,pointer :: glo2gsn(:)    ! 1d glo to 1d gsn
     integer,pointer :: gsn2glo(:)    ! 1d gsn to 1d glo
     integer,pointer :: glo2gdc(:)    ! 1d glo to 1d gdc
     integer,pointer :: gdc2glo(:)    ! 1d gdc to 1d glo
     integer,pointer :: gsn2i(:)      ! 1d gsn to 2d sn i index
     integer,pointer :: gsn2j(:)      ! 1d gsn to 2d sn j index
     integer,pointer :: gdc2i(:)      ! 1d gdc to 2d sn i index
     integer,pointer :: gdc2j(:)      ! 1d gdc to 2d sn j index
     integer,pointer :: glo2i(:)      ! 1d glo to 2d sn j index
     integer,pointer :: glo2j(:)      ! 1d glo to 2d sn j index
     integer,pointer :: ij2gsn(:,:)   ! 2d sn i,j index to 1d gsn
     integer,pointer :: ij2gdc(:,:)   ! 2d sn i,j index to 1d gdc
     integer,pointer :: ij2glo(:,:)   ! 2d sn i,j index to 1d glo
  end type decomp_type
  public decomp_type
  type(decomp_type),public,target :: ldecomp
  type(decomp_type),public,target :: adecomp

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
    integer :: cid,pid                ! indices
    integer, pointer :: lcid(:)       ! temporary for setting adecomp
    integer, pointer :: acid(:)       ! temporary for setting adecomp
    integer :: n,m                    ! indices
    integer :: ilunits, icols, ipfts  ! temporaries
    integer :: ier                    ! error code
    integer :: cnt                    ! local counter


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

          !--- find clump with fewest pfts ---
          cid = 1
          do n = 2,nclumps
             if (clumps(n)%npfts < clumps(cid)%npfts) then
                cid = n
             endif
          enddo
          acid(an) = cid

          !--- give atm cell to pe that owns cid ---
          if (iam >  clumps(cid)%owner) then
             procinfo%abegg = procinfo%abegg + 1
          endif
          if (iam >= clumps(cid)%owner) then
             procinfo%aendg = procinfo%aendg + 1
          endif

          cnt = 0
          do ln = 1,lns
          if (ldomain%gatm(ln) == an) then         
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
          endif  ! ldomain%gatm == an
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
             adecomp%gdc2glo(anumg), adecomp%gsn2glo(anumg), &
             adecomp%gdc2i  (anumg), adecomp%gsn2i  (anumg), &
             adecomp%gdc2j  (anumg), adecomp%gsn2j  (anumg), &
             stat=ier)
    if (ier /= 0) then
       write (6,*) 'decomp_init(): allocation error1 for adecomp'
       call endrun()
    end if

    allocate(adecomp%ij2gsn (ani,anj), adecomp%ij2gdc (ani,anj), &
             adecomp%ij2glo (ani,anj),                           &
             adecomp%glo2gsn(ani*anj), adecomp%glo2gdc(ani*anj), &
             adecomp%glo2i  (ani*anj), adecomp%glo2j  (ani*anj), &
             stat=ier)
    if (ier /= 0) then
       write (6,*) 'decomp_init(): allocation error2 for adecomp'
       call endrun()
    end if

    adecomp%gdc2gsn(:)  = 0
    adecomp%gdc2glo(:)  = 0
    adecomp%gdc2i(:)    = 0
    adecomp%gdc2j(:)    = 0
    adecomp%gsn2gdc(:)  = 0
    adecomp%gsn2glo(:)  = 0
    adecomp%gsn2i(:)    = 0
    adecomp%gsn2j(:)    = 0
    adecomp%glo2gsn(:)  = 0
    adecomp%glo2gdc(:)  = 0
    adecomp%glo2i(:)    = 0
    adecomp%glo2j(:)    = 0
    adecomp%ij2gsn(:,:) = 0
    adecomp%ij2gdc(:,:) = 0
    adecomp%ij2glo(:,:) = 0

    allocate(ldecomp%gdc2gsn(numg), ldecomp%gsn2gdc(numg), &
             ldecomp%gdc2glo(numg), ldecomp%gsn2glo(numg), &
             ldecomp%gdc2i  (numg), ldecomp%gsn2i  (numg), &
             ldecomp%gdc2j  (numg), ldecomp%gsn2j  (numg), &
             stat=ier)
    if (ier /= 0) then
       write (6,*) 'decomp_init(): allocation error1 for ldecomp'
       call endrun()
    end if

    allocate(ldecomp%ij2gsn (lni,lnj), ldecomp%ij2gdc (lni,lnj), &
             ldecomp%ij2glo (lni,lnj),                           &
             ldecomp%glo2gsn(lni*lnj), ldecomp%glo2gdc(lni*lnj), &
             ldecomp%glo2i  (lni*lnj), ldecomp%glo2j  (lni*lnj), &
             stat=ier)
    if (ier /= 0) then
       write (6,*) 'decomp_init(): allocation error2 for ldecomp'
       call endrun()
    end if

    ldecomp%gdc2gsn(:)  = 0
    ldecomp%gdc2glo(:)  = 0
    ldecomp%gdc2i(:)    = 0
    ldecomp%gdc2j(:)    = 0
    ldecomp%gsn2gdc(:)  = 0
    ldecomp%gsn2glo(:)  = 0
    ldecomp%gsn2i(:)    = 0
    ldecomp%gsn2j(:)    = 0
    ldecomp%glo2gsn(:)  = 0
    ldecomp%glo2gdc(:)  = 0
    ldecomp%glo2i(:)    = 0
    ldecomp%glo2j(:)    = 0
    ldecomp%ij2gsn(:,:) = 0
    ldecomp%ij2gdc(:,:) = 0
    ldecomp%ij2glo(:,:) = 0

    ag = 0
    do aj = 1,anj
    do ai = 1,ani
       an = (aj-1)*ani + ai
       adecomp%ij2glo(ai,aj) = an
       adecomp%glo2i(an) = ai
       adecomp%glo2j(an) = aj
       if (acid(an) > 0) then
          ag  = ag  + 1
          adecomp%ij2gsn(ai,aj) = ag
          adecomp%gsn2i(ag) = ai
          adecomp%gsn2j(ag) = aj
          adecomp%glo2gsn(an) = ag
          adecomp%gsn2glo(ag) = an
       endif
    enddo
    enddo

    ! Set ldecomp sn indexing based on cells to be used and i,j order
    lg  = 0
    do lj = 1,lnj
    do li = 1,lni
       ln = (lj-1)*lni + li
       ldecomp%ij2glo(li,lj) = ln
       ldecomp%glo2i(ln) = li
       ldecomp%glo2j(ln) = lj
       if (lcid(ln) > 0) then
          lg = lg + 1
          ldecomp%ij2gsn(li,lj) = lg
          ldecomp%gsn2i(lg) = li
          ldecomp%gsn2j(lg) = lj
          ldecomp%glo2gsn(ln) = lg
          ldecomp%gsn2glo(lg) = ln
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
                adecomp%gdc2gsn(ag) = adecomp%ij2gsn(ai,aj)
                adecomp%gdc2glo(ag) = adecomp%ij2glo(ai,aj)
                adecomp%ij2gdc(ai,aj) = ag
                adecomp%gsn2gdc(adecomp%ij2gsn(ai,aj)) = ag
                adecomp%glo2gdc(adecomp%ij2glo(ai,aj)) = ag
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
                ldecomp%gdc2gsn(lg) = ldecomp%ij2gsn(li,lj)
                ldecomp%gdc2glo(lg) = ldecomp%ij2glo(li,lj)
                ldecomp%ij2gdc(li,lj) = lg
                ldecomp%gsn2gdc(ldecomp%ij2gsn(li,lj)) = lg
                ldecomp%glo2gdc(ldecomp%ij2glo(li,lj)) = lg
             endif
          enddo
          enddo
       endif
    enddo
    enddo

    deallocate(acid,lcid)


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
       write (6,*)
    end if


    ! Write out clump and proc info, one pe at a time, 
    ! barrier to control pes overwriting each other on stdout

    call shr_sys_flush(6)
#if (defined SPMD)
    call mpi_barrier(mpicom,ier)
#endif
    do pid = 0,npes-1
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
               ' total landunits per proc = ',procinfo%nlunits
          write(6,*)'proc= ',pid,' beg column  = ',procinfo%begc, &
               ' end column  = ',procinfo%endc,                   &
               ' total columns per proc  = ',procinfo%ncols
          write(6,*)'proc= ',pid,' beg pft     = ',procinfo%begp, &
               ' end pft     = ',procinfo%endp,                   &
               ' total pfts per proc     = ',procinfo%npfts
          do n = 1,procinfo%nclumps
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
       call shr_sys_flush(6)
#if (defined SPMD)
       call mpi_barrier(mpicom,ier)
#endif
    end do

    call shr_sys_flush(6)


  end subroutine decomp_init

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
