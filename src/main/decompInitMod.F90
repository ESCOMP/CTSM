module decompInitMod

!------------------------------------------------------------------------------
!BOP
!
! !MODULE: decompInitMod
!
! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use spmdMod     , only : masterproc, iam, npes, mpicom, comp_id
  use clm_mct_mod
  use shr_sys_mod , only : shr_sys_flush
  use abortutils  , only : endrun
  use clm_varctl  , only : iulog
  use decompMod
!
! !PUBLIC TYPES:
  implicit none
!
! !PUBLIC MEMBER FUNCTIONS:
  public decompInit_atm         ! initializes atm grid decomposition
                                 ! into clumps and processors
  public decompInit_lnd         ! initializes lnd grid decomposition
                                 ! into clumps and processors
  public decompInit_glcp        ! initializes g,l,c,p decomp info

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
!
! !PRIVATE TYPES:
  private

  integer, pointer :: acid(:)       ! temporary for setting adecomp/ldecomp

!EOP
!------------------------------------------------------------------------------

contains

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: decompInit_atm
!
! !INTERFACE:
  subroutine decompInit_atm(alatlon,amask)
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
    logical :: seglen1                ! is segment length one
    real(r8):: seglen                 ! average segment length
    real(r8):: rcid                   ! real value of cid
    integer :: cid,pid                ! indices
    integer :: n,m,np                 ! indices
    integer :: ier                    ! error code
    integer, parameter :: dbug=1      ! 0 = min, 1=normal, 2=much, 3=max
    integer :: npmin,npmax,npint      ! do loop values for printing
    integer :: clmin,clmax,clint      ! do loop values for printing
    integer :: beg,end,lsize,gsize    ! used for gsmap init
    integer, pointer :: gindex(:)     ! global index for gsmap init


! !CALLED FROM:
! subroutine initialize
!
! !REVISION HISTORY:
! 2002.09.11  Forrest Hoffman  Creation.
! 2005.12.15  T Craig  Updated for finemesh
! 2006.08.18  P Worley Performance optimizations
! 2007.01.24  T Craig  Created decompInit_atm from decomp_init
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
          write(iulog,*) 'decompInit_atm(): Number of gridcell clumps= ',nclumps, &
               ' is less than the number of processes = ', npes
          call endrun()
       end if
    else
       write(iulog,*)'clump_pproc= ',clump_pproc,'  must be greater than 0'
       call endrun()
    end if

    !--- allocate and initialize procinfo and clumps ---
    !--- beg and end indices initialized for simple addition of cells later ---

    allocate(procinfo%cid(clump_pproc), stat=ier)
    if (ier /= 0) then
       write(iulog,*) 'decompInit_atm(): allocation error for procinfo%cid'
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
       write(iulog,*) 'decompInit_atm(): allocation error for clumps'
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
          write(iulog,*) 'decompInit_atm(): round robin pid error ',n,pid,npes
          call endrun()
       endif
       clumps(n)%owner = pid
       if (iam == pid) then
          cid = cid + 1
          if (cid < 1 .or. cid > clump_pproc) then
             write(iulog,*) 'decompInit_atm(): round robin pid error ',n,pid,npes
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

    if (npes > numa) then
       write(iulog,*) 'decompInit_atm(): Number of processes exceeds number ', &
            'of atm grid cells',npes,numa
       call endrun()
    end if

    if (nclumps > numa) then
       write(iulog,*) 'decompInit_atm(): Number of clumps exceeds number ', &
            'of atm grid cells',nclumps,numa
       call endrun()
    end if

    if (float(anumg_tot)/float(nclumps) < float(nsegspc)) then
       seglen1 = .true.
       seglen = 1.0_r8
    else
       seglen1 = .false.
       seglen = dble(anumg_tot)/(dble(nsegspc)*dble(nclumps))
    endif

    if (masterproc) write(iulog,*) ' atm decomp precompute anumg,nclumps,seglen1,avg_seglen,nsegspc=', &
      numa,nclumps,seglen1,sngl(seglen),sngl(dble(anumg_tot)/(seglen*dble(nclumps)))

    !--- assign gridcells to clumps (and thus pes) ---
    allocate(acid(ans))
    acid = 0
    anumg = 0
    do an = 1,ans
       if (amask(an) == 1) then
          anumg  = anumg  + 1

          !--- give to clumps in order based on nsegspc
          if (seglen1) then
              cid = mod(anumg-1,nclumps) + 1
          else
              rcid = (dble(anumg-1)/dble(anumg_tot))*dble(nsegspc)*dble(nclumps)
              cid = mod(int(rcid),nclumps) + 1
          endif
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
       write(iulog,*) 'decompInit_atm(): Number of atm gridcells inconsistent',anumg,numa
       call endrun()
    end if

    ! Allocate dynamic memory for adecomp, ldecomp derived type

    allocate(adecomp%gdc2glo(anumg), adecomp%glo2gdc(ani*anj), &
             stat=ier)
    if (ier /= 0) then
       write(iulog,*) 'decompInit_atm(): allocation error1 for adecomp'
       call endrun()
    end if

    adecomp%gdc2glo(:)  = 0
    adecomp%glo2gdc(:)  = 0

    ! Set adecomp
    ag = 0
    do pid = 0,npes-1
    do cid = 1,nclumps
       if (clumps(cid)%owner == pid) then

          do aj = 1,anj
          do ai = 1,ani
             an = (aj-1)*ani + ai
             if (acid(an) == cid) then
                ag = ag + 1
                adecomp%gdc2glo(ag) = an
                adecomp%glo2gdc(an) = ag
             endif
          enddo
          enddo

       endif
    enddo
    enddo

    ! set gsMap_atm_gdc2glo
    call get_proc_bounds_atm(beg, end)
    allocate(gindex(beg:end))
    do n = beg,end
       gindex(n) = adecomp%gdc2glo(n)
    enddo
    lsize = end-beg+1
    gsize = ani * anj
    call mct_gsMap_init(gsMap_atm_gdc2glo, gindex, mpicom, comp_id, lsize, gsize )
    deallocate(gindex)

    ! Diagnostic output

    if (masterproc) then
       write(iulog,*)' Atm Grid Characteristics'
       write(iulog,*)'   longitude points          = ',ani
       write(iulog,*)'   latitude points           = ',anj
       write(iulog,*)'   total number of gridcells = ',anumg
       write(iulog,*)' Decomposition Characteristics'
       write(iulog,*)'   clumps per process        = ',clump_pproc
       write(iulog,*)' gsMap Characteristics'
       write(iulog,*) '  atm gsmap glo num of segs = ',mct_gsMap_ngseg(gsMap_atm_gdc2glo)
       write(iulog,*)
    end if

    ! Write out clump and proc info, one pe at a time, 
    ! barrier to control pes overwriting each other on stdout

    call shr_sys_flush(iulog)
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
          write(iulog,*)
          write(iulog,*)'proc= ',pid,' beg atmcell = ',procinfo%abegg, &
               ' end atmcell = ',procinfo%aendg,                   &
               ' total atmcells per proc = ',procinfo%aendg-procinfo%abegg+1
          write(iulog,*)'proc= ',pid,' atm ngseg   = ',mct_gsMap_ngseg(gsMap_atm_gdc2glo), &
               ' atm nlseg   = ',mct_gsMap_nlseg(gsMap_atm_gdc2glo,iam)
          write(iulog,*)'proc= ',pid,' nclumps = ',procinfo%nclumps
       end if
       call shr_sys_flush(iulog)
       call mpi_barrier(mpicom,ier)
    end do
    call shr_sys_flush(iulog)

  end subroutine decompInit_atm

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: decompInit_lnd
!
! !INTERFACE:
  subroutine decompInit_lnd(ans,ani,anj,lns,lni,lnj)
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
             write(iulog,*) 'decompInit_lnd(): map overlap error at ',an,cnt
             call endrun()
          endif
    enddo

    allocate(ldecomp%gdc2glo(numg), ldecomp%glo2gdc(lni*lnj), &
             stat=ier)
    if (ier /= 0) then
       write(iulog,*) 'decompInit_lnd(): allocation error1 for ldecomp'
       call endrun()
    end if

    ldecomp%gdc2glo(:)  = 0
    ldecomp%glo2gdc(:)  = 0

    ! Set ldecomp
    lg = 0
    do pid = 0,npes-1
    do cid = 1,nclumps
       if (clumps(cid)%owner == pid) then

          do lj = 1,lnj
          do li = 1,lni
             ln = (lj-1)*lni + li
             if (lcid(ln) == cid) then
                lg = lg + 1
                ldecomp%gdc2glo(lg) = ln
                ldecomp%glo2gdc(ln) = lg
             endif
          enddo
          enddo
       endif
    enddo
    enddo

    deallocate(lcid)
    deallocate(lncnt,lnoff,lnmap)

    ! set gsMap_lnd_gdc2glo
    call get_proc_bounds(beg, end)
    allocate(gindex(beg:end))
    do n = beg,end
       gindex(n) = ldecomp%gdc2glo(n)
    enddo
    lsize = end-beg+1
    gsize = lni * lnj
    call mct_gsMap_init(gsMap_lnd_gdc2glo, gindex, mpicom, comp_id, lsize, gsize )
    deallocate(gindex)

    ! Diagnostic output

    if (masterproc) then
       write(iulog,*)' Atm Grid Characteristics'
       write(iulog,*)'   longitude points          = ',ani
       write(iulog,*)'   latitude points           = ',anj
       write(iulog,*)'   total number of gridcells = ',numa
       write(iulog,*)' Surface Grid Characteristics'
       write(iulog,*)'   longitude points          = ',lni
       write(iulog,*)'   latitude points           = ',lnj
       write(iulog,*)'   total number of gridcells = ',numg
       write(iulog,*)' Decomposition Characteristics'
       write(iulog,*)'   clumps per process        = ',clump_pproc
       write(iulog,*)' gsMap Characteristics'
       write(iulog,*) '  lnd gsmap glo num of segs = ',mct_gsMap_ngseg(gsMap_lnd_gdc2glo)
       write(iulog,*) '  atm gsmap glo num of segs = ',mct_gsMap_ngseg(gsMap_atm_gdc2glo)
       write(iulog,*)
    end if

    call shr_sys_flush(iulog)

  end subroutine decompInit_lnd

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: decompInit_glcp
!
! !INTERFACE:
  subroutine decompInit_glcp(ans,ani,anj,lns,lni,lnj,glcmask)
!
! !DESCRIPTION:
! This subroutine initializes the land surface decomposition into a clump
! data structure.  This assumes each pe has the same number of clumps
! set by clump_pproc
!
! !USES:
    use clmtype   , only : grlnd, nameg, namel, namec, namep
    use spmdMod
    use spmdGathScatMod
    use domainMod , only : gatm
    use subgridMod, only : subgrid_get_gcellinfo
!
! !ARGUMENTS:
    implicit none
    integer , intent(in) :: lns,lni,lnj ! land domain global size
    integer , intent(in) :: ans,ani,anj ! atm domain global size
    integer , pointer, optional :: glcmask(:)  ! glc mask
!
! !LOCAL VARIABLES:
    integer :: lg,ln,li,lj        ! indices
    integer :: ag,an,ai,aj        ! indices
    integer :: abegg,aendg,anumg  ! atm num gridcells
    integer :: begg,endg          ! lnd num gridcells
    integer :: begl,endl          ! lnd num gridcells
    integer :: begc,endc          ! lnd num gridcells
    integer :: begp,endp          ! lnd num gridcells
    integer :: cid,pid                ! indices
    integer :: n,m,np                 ! indices
    integer :: icells, ilunits, icols, ipfts  ! temporaries
    integer :: ier                    ! error code
    integer :: cnt                    ! local counter
    integer, parameter :: dbug=1      ! 0 = min, 1=normal, 2=much, 3=max
    integer :: npmin,npmax,npint      ! do loop values for printing
    integer :: clmin,clmax,clint      ! do loop values for printing
    integer :: lsize,gsize            ! used for gsmap init
    integer :: ng                     ! number of gridcells in gsmap
    integer, pointer :: gindex(:)     ! global index for gsmap init

    integer, pointer :: lncnt(:)      ! lnd cell count per atm cell
    integer, pointer :: lnoff(:)      ! atm cell offset in lnmap
    integer, pointer :: lnmap(:)      ! map from atm cell to lnd cells
    integer, allocatable :: allvecg(:,:)  ! temporary vector "global"
    integer, allocatable :: allvecl(:,:)  ! temporary vector "local"
    integer :: lnidx

    integer, pointer :: arrayg(:)
    integer :: val1, val2
    integer :: i,g,l,c,p,k
    integer, pointer :: gstart(:),gcount(:)
    integer, pointer :: lstart(:),lcount(:)
    integer, pointer :: cstart(:),ccount(:)
    integer, pointer :: pstart(:),pcount(:)
    integer :: beg,end,num
    type(mct_gsmap),pointer :: gsmap
    integer, pointer :: start(:),count(:)
    integer, pointer :: tarr1(:),tarr2(:)
    integer :: ntest
    character(len=8) :: clmlevel
    character(len=32), parameter :: subname = 'decompInit_glcp'

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
    allocate(gstart(begg:endg),lstart(begg:endg),cstart(begg:endg),pstart(begg:endg))
    allocate(gcount(begg:endg),lcount(begg:endg),ccount(begg:endg),pcount(begg:endg))
    allocate(allvecg(nclumps,4),allvecl(nclumps,4))   ! 3 = gcells,lunit,cols,pfts
    allvecg  = 0
    allvecl = 0
    gcount = 0
    lcount = 0
    ccount = 0
    pcount = 0
    do anumg = abegg,aendg
          an = adecomp%gdc2glo(anumg)
          cid = acid(an)
          do lnidx = 0,lncnt(an)-1
             ln = ldecomp%glo2gdc(lnmap(lnoff(an)+lnidx))
             if (present(glcmask)) then
                call subgrid_get_gcellinfo (ln, nlunits=ilunits, &
                                  ncols=icols, npfts=ipfts, glcmask=glcmask(ln))
             else
                call subgrid_get_gcellinfo (ln, nlunits=ilunits, &
                                  ncols=icols, npfts=ipfts)
             endif
             allvecl(cid,1) = allvecl(cid,1) + 1
             allvecl(cid,2) = allvecl(cid,2) + ilunits
             allvecl(cid,3) = allvecl(cid,3) + icols
             allvecl(cid,4) = allvecl(cid,4) + ipfts
             gcount(ln) = 1
             lcount(ln) = ilunits
             ccount(ln) = icols
             pcount(ln) = ipfts
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

    do n = 1,nclumps
       if (clumps(n)%ncells  /= allvecg(n,1) .or. &
           clumps(n)%nlunits /= allvecg(n,2) .or. &
           clumps(n)%ncols   /= allvecg(n,3) .or. &
           clumps(n)%npfts   /= allvecg(n,4)) then
          write(iulog,*) 'decompInit_glcp(): allvecg error ncells ',iam,n,clumps(n)%ncells ,allvecg(n,1)
          write(iulog,*) 'decompInit_glcp(): allvecg error lunits ',iam,n,clumps(n)%nlunits,allvecg(n,2)
          write(iulog,*) 'decompInit_glcp(): allvecg error ncols  ',iam,n,clumps(n)%ncols  ,allvecg(n,3)
          write(iulog,*) 'decompInit_glcp(): allvecg error pfts   ',iam,n,clumps(n)%npfts  ,allvecg(n,4)
          call endrun()
       endif
    enddo

    deallocate(allvecg,allvecl)
    deallocate(acid)
    deallocate(lncnt,lnoff,lnmap)

    ! set gsMaps, perms for lun, col, pft

    ! this was just "set" above in procinfo, be careful not to move it up
    call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)

    ng = mct_gsmap_gsize(gsmap_lnd_gdc2glo)
    allocate(arrayg(ng))

    ! for each subgrid gsmap (l, c, p)
    ! gather the gdc subgrid counts to masterproc in glo order
    ! compute glo ordered start indices from the counts
    ! scatter the subgrid start indices back out to the gdc gridcells
    ! set the local gindex array for the subgrid from the subgrid start and count arrays

    do k = 1,4
       if (k == 1) then
          clmlevel = nameg
          beg = begg
          end = endg
          num = numg
          gsmap => gsmap_gce_gdc2glo
          start => gstart
          count => gcount
       elseif (k == 2) then
          clmlevel = namel
          beg = begl
          end = endl
          num = numl
          gsmap => gsmap_lun_gdc2glo
          start => lstart
          count => lcount
       elseif (k == 3) then
          clmlevel = namec
          beg = begc
          end = endc
          num = numc
          gsmap => gsmap_col_gdc2glo
          start => cstart
          count => ccount
       elseif (k == 4) then
          clmlevel = namep
          beg = begp
          end = endp
          num = nump
          gsmap => gsmap_pft_gdc2glo
          start => pstart
          count => pcount
       else
          write(iulog,*) 'decompInit_glcp error in k ',k
          call endrun()
       endif

       arrayg = 0
       call gather_data_to_master(count,arrayg,grlnd)
       if (masterproc) then
          gsize = arrayg(1)
          val1 = arrayg(1)
          arrayg(1) = 1
          do n = 2,ng
             gsize = gsize + arrayg(n)
             val2 = arrayg(n)
             arrayg(n) = arrayg(n-1) + val1
             val1 = val2
          enddo
       endif
       call scatter_data_from_master(start,arrayg,grlnd)

       allocate(gindex(beg:end))
       i = beg-1
       do g = begg,endg
          if (count(g) <  1) then
             write(iulog,*) 'decompInit_glcp warning count g ',k,iam,g,count(g)
          endif
          do l = 1,count(g)
             i = i + 1
             if (i < beg .or. i > end) then
                write(iulog,*) 'decompInit_glcp error i ',i,beg,end
                call endrun()
             endif
             gindex(i) = start(g) + l-1
          enddo
       enddo
      if (i /= end) then
         write(iulog,*) 'decompInit_glcp error size ',i,beg,end
         call endrun()
      endif

      lsize = end-beg+1
      gsize = num

      call mct_gsMap_init(gsMap, gindex, mpicom, comp_id, lsize, gsize )

      !--- test gsmap ---
      ntest = mct_gsMap_gsize(gsMap)
      allocate(tarr1(ntest),tarr2(beg:end))
      call gather_data_to_master(gindex,tarr1,clmlevel)
      call scatter_data_from_master(tarr2,tarr1,clmlevel)
      !--- verify gather/scatter produces same result
      do l = beg,end
         if (tarr2(l) /= gindex(l)) then
            write(iulog,*) 'decompInit_glcp error tarr2 ',k,l,gindex(l),tarr2(l)
            call endrun()
         endif
      enddo
      !--- verify gather of gindex on new gsmap produces ordered indices
      if (masterproc) then
         if (tarr1(1) /= 1) then
            write(iulog,*) 'decompInit_glcp error tarr1 ',k,1,tarr1(1)
            call endrun()
         endif
         do l = 2,ntest
            if (tarr1(l)-tarr1(l-1) /= 1) then
               write(iulog,*) 'decompInit_glcp error tarr1 ',k,l,tarr1(l-1),tarr1(l)
               call endrun()
            endif
         enddo
      endif
      deallocate(tarr1,tarr2)
      if (masterproc) then
         write(iulog,*) 'decompInit_glcp gsmap [l,c,p] test passes for ',k
      endif
      !--- end test section      

      deallocate(gindex)

    enddo

    deallocate(gstart,lstart,cstart,pstart)
    deallocate(gcount,lcount,ccount,pcount)

     ! Diagnostic output

    if (masterproc) then
       write(iulog,*)' Atm Grid Characteristics'
       write(iulog,*)'   longitude points          = ',ani
       write(iulog,*)'   latitude points           = ',anj
       write(iulog,*)'   total number of gridcells = ',numa
       write(iulog,*)' Surface Grid Characteristics'
       write(iulog,*)'   longitude points          = ',lni
       write(iulog,*)'   latitude points           = ',lnj
       write(iulog,*)'   total number of gridcells = ',numg
       write(iulog,*)'   total number of landunits = ',numl
       write(iulog,*)'   total number of columns   = ',numc
       write(iulog,*)'   total number of pfts      = ',nump
       write(iulog,*)' Decomposition Characteristics'
       write(iulog,*)'   clumps per process        = ',clump_pproc
       write(iulog,*)' gsMap Characteristics'
       write(iulog,*) '  atm gsmap glo num of segs = ',mct_gsMap_ngseg(gsMap_atm_gdc2glo)
       write(iulog,*) '  lnd gsmap glo num of segs = ',mct_gsMap_ngseg(gsMap_lnd_gdc2glo)
       write(iulog,*) '  gce gsmap glo num of segs = ',mct_gsMap_ngseg(gsMap_gce_gdc2glo)
       write(iulog,*) '  lun gsmap glo num of segs = ',mct_gsMap_ngseg(gsMap_lun_gdc2glo)
       write(iulog,*) '  col gsmap glo num of segs = ',mct_gsMap_ngseg(gsMap_col_gdc2glo)
       write(iulog,*) '  pft gsmap glo num of segs = ',mct_gsMap_ngseg(gsMap_pft_gdc2glo)
       write(iulog,*)
    end if

    ! Write out clump and proc info, one pe at a time, 
    ! barrier to control pes overwriting each other on stdout

    call shr_sys_flush(iulog)
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
          write(iulog,*)
          write(iulog,*)'proc= ',pid,' beg atmcell = ',procinfo%abegg, &
               ' end atmcell = ',procinfo%aendg,                   &
               ' total atmcells per proc = ',procinfo%aendg-procinfo%abegg+1
          write(iulog,*)'proc= ',pid,' beg gridcell= ',procinfo%begg, &
               ' end gridcell= ',procinfo%endg,                   &
               ' total gridcells per proc= ',procinfo%ncells
          write(iulog,*)'proc= ',pid,' beg landunit= ',procinfo%begl, &
               ' end landunit= ',procinfo%endl,                   &
               ' total landunits per proc= ',procinfo%nlunits
          write(iulog,*)'proc= ',pid,' beg column  = ',procinfo%begc, &
               ' end column  = ',procinfo%endc,                   &
               ' total columns per proc  = ',procinfo%ncols
          write(iulog,*)'proc= ',pid,' beg pft     = ',procinfo%begp, &
               ' end pft     = ',procinfo%endp,                   &
               ' total pfts per proc     = ',procinfo%npfts
          write(iulog,*)'proc= ',pid,' atm ngseg   = ',mct_gsMap_ngseg(gsMap_atm_gdc2glo), &
               ' atm nlseg   = ',mct_gsMap_nlseg(gsMap_atm_gdc2glo,iam)
          write(iulog,*)'proc= ',pid,' lnd ngseg   = ',mct_gsMap_ngseg(gsMap_lnd_gdc2glo), &
               ' lnd nlseg   = ',mct_gsMap_nlseg(gsMap_lnd_gdc2glo,iam)
          write(iulog,*)'proc= ',pid,' gce ngseg   = ',mct_gsMap_ngseg(gsMap_gce_gdc2glo), &
               ' gce nlseg   = ',mct_gsMap_nlseg(gsMap_gce_gdc2glo,iam)
          write(iulog,*)'proc= ',pid,' lun ngseg   = ',mct_gsMap_ngseg(gsMap_lun_gdc2glo), &
               ' lun nlseg   = ',mct_gsMap_nlseg(gsMap_lun_gdc2glo,iam)
          write(iulog,*)'proc= ',pid,' col ngseg   = ',mct_gsMap_ngseg(gsMap_col_gdc2glo), &
               ' col nlseg   = ',mct_gsMap_nlseg(gsMap_col_gdc2glo,iam)
          write(iulog,*)'proc= ',pid,' pft ngseg   = ',mct_gsMap_ngseg(gsMap_pft_gdc2glo), &
               ' pft nlseg   = ',mct_gsMap_nlseg(gsMap_pft_gdc2glo,iam)
          write(iulog,*)'proc= ',pid,' nclumps = ',procinfo%nclumps

          clmin = 1
          clmax = procinfo%nclumps
          if (dbug == 1) then
            clmax = 1
          elseif (dbug == 0) then
            clmax = -1
          endif
          do n = clmin,clmax
             cid = procinfo%cid(n)
             write(iulog,*)'proc= ',pid,' clump no = ',n, &
                  ' clump id= ',procinfo%cid(n),    &
                  ' beg gridcell= ',clumps(cid)%begg, &
                  ' end gridcell= ',clumps(cid)%endg, &
                  ' total gridcells per clump= ',clumps(cid)%ncells
             write(iulog,*)'proc= ',pid,' clump no = ',n, &
                  ' clump id= ',procinfo%cid(n),    &
                  ' beg landunit= ',clumps(cid)%begl, &
                  ' end landunit= ',clumps(cid)%endl, &
                  ' total landunits per clump = ',clumps(cid)%nlunits
             write(iulog,*)'proc= ',pid,' clump no = ',n, &
                  ' clump id= ',procinfo%cid(n),    &
                  ' beg column  = ',clumps(cid)%begc, &
                  ' end column  = ',clumps(cid)%endc, &
                  ' total columns per clump  = ',clumps(cid)%ncols
             write(iulog,*)'proc= ',pid,' clump no = ',n, &
                  ' clump id= ',procinfo%cid(n),    &
                  ' beg pft     = ',clumps(cid)%begp, &
                  ' end pft     = ',clumps(cid)%endp, &
                  ' total pfts per clump     = ',clumps(cid)%npfts
          end do
       end if
       call shr_sys_flush(iulog)
       call mpi_barrier(mpicom,ier)
    end do
    call shr_sys_flush(iulog)

  end subroutine decompInit_glcp

!------------------------------------------------------------------------------

end module decompInitMod
