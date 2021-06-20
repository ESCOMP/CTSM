module decompInitMod

  !------------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module provides a descomposition into a clumped data structure which can
  ! be mapped back to atmosphere physics chunks.
  !
  ! !USES:
  use shr_kind_mod           , only : r8 => shr_kind_r8
  use shr_sys_mod            , only : shr_sys_flush
  use shr_log_mod            , only : errMsg => shr_log_errMsg
  use spmdMod                , only : masterproc, iam, npes, mpicom
  use abortutils             , only : endrun
  use clm_varctl             , only : iulog, use_fates
  use LandunitType           , only : lun
  use ColumnType             , only : col
  use PatchType              , only : patch
  use glcBehaviorMod         , only : glc_behavior_type
  use FatesInterfaceTypesMod , only : fates_maxElementsPerSite
  use decompMod
  !
  implicit none
  private
  !
  ! !PUBLIC TYPES:
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public decompInit_lnd    ! initializes lnd grid decomposition into clumps and processors
  public decompInit_clumps ! initializes atm grid decomposition into clumps
  public decompInit_glcp   ! initializes g,l,c,p decomp info
  !
  ! !PRIVATE MEMBER FUNCTIONS:
  private :: set_subgrid_start
  private :: set_gsmap_global
  !
  ! !PRIVATE TYPES:
  integer, pointer :: lcid(:)       ! temporary for setting decomposition
  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  !------------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------------
  subroutine decompInit_lnd(lni,lnj,amask)
    !
    ! !DESCRIPTION:
    ! This subroutine initializes the land surface decomposition into a clump
    ! data structure.  This assumes each pe has the same number of clumps
    ! set by clump_pproc
    !
    ! !USES:
    use clm_varctl      , only : nsegspc
    use decompMod       , only : nglob_x, nglob_y, gindex_global
    !
    ! !ARGUMENTS:
    implicit none
    integer , intent(in) :: amask(:)
    integer , intent(in) :: lni,lnj   ! domain global size
    !
    ! !LOCAL VARIABLES:
    integer :: lns                    ! global domain size
    integer :: ln                     ! indices
    integer :: ag,an,ai,aj            ! indices
    integer :: numg                   ! number of land gridcells
    logical :: seglen1                ! is segment length one
    real(r8):: seglen                 ! average segment length
    real(r8):: rcid                   ! real value of cid
    integer :: cid,pid                ! indices
    integer :: n,m,ng                 ! indices
    integer :: ier                    ! error code
    integer :: begg, endg             ! beg and end gridcells
    integer, pointer  :: clumpcnt(:)  ! clump index counter
    integer, allocatable :: gdc2glo(:)! used to create gindex_global
    type(bounds_type) :: bounds       ! contains subgrid bounds data 
    !------------------------------------------------------------------------------

    lns = lni * lnj

    !--- set and verify nclumps ---
    if (clump_pproc > 0) then
       nclumps = clump_pproc * npes
       if (nclumps < npes) then
          write(iulog,*) 'decompInit_lnd(): Number of gridcell clumps= ',nclumps, &
               ' is less than the number of processes = ', npes
          call endrun(msg=errMsg(sourcefile, __LINE__))
       end if
    else
       write(iulog,*)'clump_pproc= ',clump_pproc,'  must be greater than 0'
       call endrun(msg=errMsg(sourcefile, __LINE__))
    end if

    ! allocate and initialize procinfo and clumps
    ! beg and end indices initialized for simple addition of cells later

    allocate(procinfo%cid(clump_pproc), stat=ier)
    if (ier /= 0) then
       write(iulog,*) 'decompInit_lnd(): allocation error for procinfo%cid'
       call endrun(msg=errMsg(sourcefile, __LINE__))
    endif
    procinfo%nclumps   = clump_pproc
    procinfo%cid(:)    = -1
    procinfo%ncells    = 0
    procinfo%nlunits   = 0
    procinfo%ncols     = 0
    procinfo%npatches  = 0
    procinfo%nCohorts  = 0
    procinfo%begg      = 1
    procinfo%begl      = 1
    procinfo%begc      = 1
    procinfo%begp      = 1
    procinfo%begCohort = 1
    procinfo%endg      = 0
    procinfo%endl      = 0
    procinfo%endc      = 0
    procinfo%endp      = 0
    procinfo%endCohort = 0

    allocate(clumps(nclumps), stat=ier)
    if (ier /= 0) then
       write(iulog,*) 'decompInit_lnd(): allocation error for clumps'
       call endrun(msg=errMsg(sourcefile, __LINE__))
    end if
    clumps(:)%owner     = -1
    clumps(:)%ncells    = 0
    clumps(:)%nlunits   = 0
    clumps(:)%ncols     = 0
    clumps(:)%npatches  = 0
    clumps(:)%nCohorts  = 0
    clumps(:)%begg      = 1
    clumps(:)%begl      = 1
    clumps(:)%begc      = 1
    clumps(:)%begp      = 1
    clumps(:)%begCohort = 1
    clumps(:)%endg      = 0
    clumps(:)%endl      = 0
    clumps(:)%endc      = 0
    clumps(:)%endp      = 0
    clumps(:)%endCohort = 0

    ! assign clumps to proc round robin
    cid = 0
    do n = 1,nclumps
       pid = mod(n-1,npes)
       if (pid < 0 .or. pid > npes-1) then
          write(iulog,*) 'decompInit_lnd(): round robin pid error ',n,pid,npes
          call endrun(msg=errMsg(sourcefile, __LINE__))
       endif
       clumps(n)%owner = pid
       if (iam == pid) then
          cid = cid + 1
          if (cid < 1 .or. cid > clump_pproc) then
             write(iulog,*) 'decompInit_lnd(): round robin pid error ',n,pid,npes
             call endrun(msg=errMsg(sourcefile, __LINE__))
          endif
          procinfo%cid(cid) = n
       endif
    enddo

    ! count total land gridcells
    numg = 0
    do ln = 1,lns
       if (amask(ln) == 1) then
          numg = numg + 1
       endif
    enddo

    if (npes > numg) then
       write(iulog,*) 'decompInit_lnd(): Number of processes exceeds number ', &
            'of land grid cells',npes,numg
       call endrun(msg=errMsg(sourcefile, __LINE__))
    end if
    if (nclumps > numg) then
       write(iulog,*) 'decompInit_lnd(): Number of clumps exceeds number ', &
            'of land grid cells',nclumps,numg
       call endrun(msg=errMsg(sourcefile, __LINE__))
    end if

    if (float(numg)/float(nclumps) < float(nsegspc)) then
       seglen1 = .true.
       seglen = 1.0_r8
    else
       seglen1 = .false.
       seglen = dble(numg)/(dble(nsegspc)*dble(nclumps))
    endif

    if (masterproc) then
       write(iulog,*) ' decomp precompute numg,nclumps,seglen1,avg_seglen,nsegspc=', &
            numg,nclumps,seglen1,&
            sngl(seglen),sngl(dble(numg)/(seglen*dble(nclumps)))
    end if

    ! Assign gridcells to clumps (and thus pes) ---

    allocate(lcid(lns))
    lcid(:) = 0
    ng = 0
    do ln = 1,lns
       if (amask(ln) == 1) then
          ng = ng  + 1

          !--- give to clumps in order based on nsegspc
          if (seglen1) then
             cid = mod(ng-1,nclumps) + 1
          else
             rcid = (dble(ng-1)/dble(numg))*dble(nsegspc)*dble(nclumps)
             cid = mod(int(rcid),nclumps) + 1
          endif
          lcid(ln) = cid

          !--- give gridcell cell to pe that owns cid ---
          !--- this needs to be done to subsequently use function
          !--- get_proc_bounds(begg,endg)
          if (iam == clumps(cid)%owner) then
             procinfo%ncells  = procinfo%ncells  + 1
          endif
          if (iam >  clumps(cid)%owner) then
             procinfo%begg = procinfo%begg + 1
          endif
          if (iam >= clumps(cid)%owner) then
             procinfo%endg = procinfo%endg + 1
          endif

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

       end if
    enddo

    ! Set gindex_global

    allocate(gdc2glo(numg), stat=ier)  
    if (ier /= 0) then
       write(iulog,*) 'decompInit_lnd(): allocation error1 for gdc2glo , etc'
       call endrun(msg=errMsg(sourcefile, __LINE__))
    end if
    gdc2glo(:) = 0
    allocate(clumpcnt(nclumps),stat=ier)
    if (ier /= 0) then
       write(iulog,*) 'decompInit_lnd(): allocation error1 for clumpcnt'
       call endrun(msg=errMsg(sourcefile, __LINE__))
    end if
    
    ! clumpcnt is the start gdc index of each clump

    ag = 0
    clumpcnt = 0
    ag = 1
    do pid = 0,npes-1
    do cid = 1,nclumps
       if (clumps(cid)%owner == pid) then
         clumpcnt(cid) = ag
         ag = ag + clumps(cid)%ncells
       endif
    enddo
    enddo

    ! now go through gridcells one at a time and increment clumpcnt
    ! in order to set gdc2glo

    do aj = 1,lnj
    do ai = 1,lni
       an = (aj-1)*lni + ai
       cid = lcid(an)
       if (cid > 0) then
          ag = clumpcnt(cid)
          gdc2glo(ag) = an
          clumpcnt(cid) = clumpcnt(cid) + 1
       end if
    end do
    end do

    ! Initialize global gindex (non-compressed, includes ocean points)
    ! Note that gindex_global goes from (1:endg-begg_1)
    nglob_x = lni !  decompMod module variables
    nglob_y = lnj !  decompMod module variables
    call get_proc_bounds(bounds)
    begg = bounds%begg; endg = bounds%endg
    allocate(gindex_global(endg-begg+1))
    do n = begg,endg
       gindex_global(n-begg+1) = gdc2glo(n)
    enddo

    deallocate(clumpcnt)
    deallocate(gdc2glo)

    ! Initialize gsmap_global
    call set_gsmap_global(gindex_global)

    ! Diagnostic output
    if (masterproc) then
       write(iulog,*)' Surface Grid Characteristics'
       write(iulog,*)'   longitude points               = ',lni
       write(iulog,*)'   latitude points                = ',lnj
       write(iulog,*)'   total number of land gridcells = ',numg
       write(iulog,*)' Decomposition Characteristics'
       write(iulog,*)'   clumps per process             = ',clump_pproc
       write(iulog,*)
    end if
    call shr_sys_flush(iulog)

  end subroutine decompInit_lnd

  !------------------------------------------------------------------------------
  subroutine decompInit_clumps(lni,lnj,glc_behavior)
    !
    ! !DESCRIPTION:
    ! This subroutine initializes the land surface decomposition into a clump
    ! data structure.  This assumes each pe has the same number of clumps
    ! set by clump_pproc
    !
    ! !USES:
    use subgridMod, only : subgrid_get_gcellinfo
    use spmdMod
    !
    ! !ARGUMENTS:
    implicit none
    integer , intent(in) :: lni,lnj ! land domain global size
    type(glc_behavior_type), intent(in) :: glc_behavior
    !
    ! !LOCAL VARIABLES:
    integer :: ln,an              ! indices
    integer :: i,g,l,k            ! indices
    integer :: cid                ! indices
    integer :: n,m,np             ! indices
    integer :: anumg              ! lnd num gridcells
    integer :: icells             ! temporary
    integer :: begg, endg         ! temporary
    integer :: ilunits            ! temporary
    integer :: icols              ! temporary
    integer :: ipatches           ! temporary
    integer :: icohorts           ! temporary
    integer :: ier                ! error code
    integer, allocatable :: allvecg(:,:)  ! temporary vector "global"
    integer, allocatable :: allvecl(:,:)  ! temporary vector "local"
    integer :: ntest
    type(bounds_type) :: bounds
    character(len=32), parameter :: subname = 'decompInit_clumps'
    !------------------------------------------------------------------------------

    !--- assign gridcells to clumps (and thus pes) ---
    call get_proc_bounds(bounds)
    begg = bounds%begg; endg = bounds%endg    

    allocate(allvecl(nclumps,5))   ! local  clumps [gcells,lunit,cols,patches,coh]
    allocate(allvecg(nclumps,5))   ! global clumps [gcells,lunit,cols,patches,coh]

    ! Determine the number of gridcells, landunits, columns, and patches, cohorts
    ! on this processor
    ! Determine number of landunits, columns and patches for each global
    ! gridcell index (an) that is associated with the local gridcell index (ln)

    ilunits=0
    icols=0
    ipatches=0
    icohorts=0

    allvecg= 0
    allvecl= 0
    do anumg = begg,endg
       an  = gindex_global(anumg - begg + 1)
       cid = lcid(an)
       ln  = anumg
       call subgrid_get_gcellinfo (ln, nlunits=ilunits, ncols=icols, npatches=ipatches, &
            ncohorts=icohorts, glc_behavior=glc_behavior)
       allvecl(cid,1) = allvecl(cid,1) + 1
       allvecl(cid,2) = allvecl(cid,2) + ilunits  ! number of landunits for local clump cid
       allvecl(cid,3) = allvecl(cid,3) + icols    ! number of columns for local clump cid
       allvecl(cid,4) = allvecl(cid,4) + ipatches ! number of patches for local clump cid
       allvecl(cid,5) = allvecl(cid,5) + icohorts ! number of cohorts for local clump cid
    enddo
    call mpi_allreduce(allvecl,allvecg,size(allvecg),MPI_INTEGER,MPI_SUM,mpicom,ier)

    ! Determine overall  total gridcells, landunits, columns and patches and distribute
    ! gridcells over clumps

    numg = 0
    numl = 0
    numc = 0
    nump = 0
    numCohort = 0

    do cid = 1,nclumps
       icells   = allvecg(cid,1)  ! number of all clump cid gridcells (over all processors)
       ilunits  = allvecg(cid,2)  ! number of all clump cid landunits (over all processors)
       icols    = allvecg(cid,3)  ! number of all clump cid columns (over all processors)
       ipatches = allvecg(cid,4)  ! number of all clump cid patches (over all processors)
       icohorts = allvecg(cid,5)  ! number of all clump cid cohorts (over all processors)

       !--- overall total ---
       numg = numg + icells             ! total number of gridcells
       numl = numl + ilunits            ! total number of landunits
       numc = numc + icols              ! total number of columns
       nump = nump + ipatches           ! total number of patches
       numCohort = numCohort + icohorts ! total number of cohorts

       !--- give gridcell to cid ---
       !--- increment the beg and end indices ---
       clumps(cid)%nlunits  = clumps(cid)%nlunits  + ilunits
       clumps(cid)%ncols    = clumps(cid)%ncols    + icols
       clumps(cid)%npatches = clumps(cid)%npatches    + ipatches
       clumps(cid)%nCohorts = clumps(cid)%nCohorts + icohorts

       do m = 1,nclumps
          if ((clumps(m)%owner >  clumps(cid)%owner) .or. &
              (clumps(m)%owner == clumps(cid)%owner .and. m > cid)) then
             clumps(m)%begl = clumps(m)%begl + ilunits
             clumps(m)%begc = clumps(m)%begc + icols
             clumps(m)%begp = clumps(m)%begp + ipatches
             clumps(m)%begCohort = clumps(m)%begCohort + icohorts
          endif

          if ((clumps(m)%owner >  clumps(cid)%owner) .or. &
              (clumps(m)%owner == clumps(cid)%owner .and. m >= cid)) then
             clumps(m)%endl = clumps(m)%endl + ilunits
             clumps(m)%endc = clumps(m)%endc + icols
             clumps(m)%endp = clumps(m)%endp + ipatches
             clumps(m)%endCohort = clumps(m)%endCohort + icohorts
          endif
       enddo

       !--- give gridcell to the proc that owns the cid ---
       !--- increment the beg and end indices ---
       if (iam == clumps(cid)%owner) then
          procinfo%nlunits  = procinfo%nlunits  + ilunits
          procinfo%ncols    = procinfo%ncols    + icols
          procinfo%npatches = procinfo%npatches + ipatches
          procinfo%nCohorts = procinfo%nCohorts + icohorts
       endif

       if (iam >  clumps(cid)%owner) then
          procinfo%begl = procinfo%begl + ilunits
          procinfo%begc = procinfo%begc + icols
          procinfo%begp = procinfo%begp + ipatches
          procinfo%begCohort = procinfo%begCohort + icohorts
       endif

       if (iam >= clumps(cid)%owner) then
          procinfo%endl = procinfo%endl + ilunits
          procinfo%endc = procinfo%endc + icols
          procinfo%endp = procinfo%endp + ipatches
          procinfo%endCohort = procinfo%endCohort + icohorts
       endif
    enddo

    do n = 1,nclumps
       if (clumps(n)%ncells   /= allvecg(n,1) .or. &
           clumps(n)%nlunits  /= allvecg(n,2) .or. &
           clumps(n)%ncols    /= allvecg(n,3) .or. &
           clumps(n)%npatches /= allvecg(n,4) .or. &
           clumps(n)%nCohorts /= allvecg(n,5)) then

          write(iulog ,*) 'decompInit_glcp(): allvecg error ncells ',iam,n,clumps(n)%ncells   ,allvecg(n,1)
          write(iulog ,*) 'decompInit_glcp(): allvecg error lunits ',iam,n,clumps(n)%nlunits  ,allvecg(n,2)
          write(iulog ,*) 'decompInit_glcp(): allvecg error ncols  ',iam,n,clumps(n)%ncols    ,allvecg(n,3)
          write(iulog ,*) 'decompInit_glcp(): allvecg error patches',iam,n,clumps(n)%npatches ,allvecg(n,4)
          write(iulog ,*) 'decompInit_glcp(): allvecg error cohorts',iam,n,clumps(n)%nCohorts ,allvecg(n,5)

          call endrun(msg=errMsg(sourcefile, __LINE__))
       endif
    enddo

    deallocate(allvecg,allvecl)
    deallocate(lcid)

  end subroutine decompInit_clumps

  !------------------------------------------------------------------------------
  subroutine decompInit_glcp(lni,lnj,glc_behavior)
    !
    ! !DESCRIPTION:
    ! Determine gindex for landunits, columns, patches and cohorts
    !
    ! !USES:
    use spmdMod
    use subgridMod, only : subgrid_get_gcellinfo
    !
    ! !ARGUMENTS:
    integer                 , intent(in) :: lni,lnj ! land domain global size
    type(glc_behavior_type) , intent(in) :: glc_behavior
    !
    ! !LOCAL VARIABLES:
    integer :: gi,li,ci,pi,coi           ! indices
    integer :: i,l,n,np                  ! indices
    integer :: cid,pid                   ! indices
    integer :: begg,endg                 ! beg,end gridcells
    integer :: begl,endl                 ! beg,end landunits
    integer :: begc,endc                 ! beg,end columns
    integer :: begp,endp                 ! beg,end patches
    integer :: begCohort,endCohort       ! beg,end cohorts
    integer :: numg                      ! total number of land gridcells across all processors
    integer :: numl                      ! total number of landunits across all processors
    integer :: numc                      ! total number of columns across all processors
    integer :: nump                      ! total number of patches across all processors
    integer :: numCohort                 ! fates cohorts
    integer :: ilunits                   ! temporary
    integer :: icols                     ! temporary
    integer :: ipatches                  ! temporary
    integer :: icohorts                  ! temporary
    integer :: ier                       ! error code
    integer :: npmin,npmax,npint         ! do loop values for printing
    integer :: clmin,clmax               ! do loop values for printing
    integer :: ng                        ! number of global gridcells
    integer, pointer :: gindex_global(:) ! global index
    integer, pointer :: array_glob(:)    ! temporaroy
    integer, pointer :: gstart(:), gcount(:)
    integer, pointer :: lstart(:), lcount(:)
    integer, pointer :: cstart(:), ccount(:)
    integer, pointer :: pstart(:), pcount(:)
    integer, pointer :: coStart(:), coCount(:)
    integer, pointer :: ioff(:)
    type(bounds_type):: bounds
    integer, parameter :: dbug=1         ! 0 = min, 1=normal, 2=much, 3=max
    character(len=32), parameter :: subname = 'decompInit_glcp'
    !------------------------------------------------------------------------------

    ! Get processor bounds
    
    call get_proc_bounds(bounds)
    begg = bounds%begg; endg = bounds%endg
    begl = bounds%begl; endl = bounds%endl
    begc = bounds%begc; endc = bounds%endc
    begp = bounds%begp; endp = bounds%endp
    begCohort = bounds%begCoHort; endCohort = bounds%endCoHort

    call get_proc_global(ng=numg, nl=numl, nc=numc, np=nump, nCohorts=numCohort)

    ! Allocate start and count for determining subgrid level global index space
    allocate(gcount(begg:endg)) ; gcount(:) = 0
    allocate(lcount(begg:endg)) ; lcount(:) = 0
    allocate(ccount(begg:endg)) ; ccount(:) = 0
    allocate(pcount(begg:endg)) ; pcount(:) = 0
    allocate(coCount(begg:endg)); coCount(:) = 0

    ! Determine gcount, lcount, ccount and pcount
    do gi = begg,endg
       call subgrid_get_gcellinfo (gi, nlunits=ilunits, ncols=icols, npatches=ipatches, &
            ncohorts=icohorts, glc_behavior=glc_behavior)
       gcount(gi)  = 1         ! number of gridcells for local gridcell index gi
       lcount(gi)  = ilunits   ! number of landunits for local gridcell index gi
       ccount(gi)  = icols     ! number of columns for local gridcell index gi
       pcount(gi)  = ipatches  ! number of patches for local gridcell index gi
       coCount(gi) = icohorts  ! number of fates cohorts for local gricell index gi
    enddo

    ! Determine gstart, lstart, cstart, pstart, coStart for the OUTPUT 1d data structures
    ! gather the gdc subgrid counts to masterproc in glo order
    ! compute glo ordered start indices from the counts
    ! scatter the subgrid start indices back out to the gdc gridcells
    ! set the local gindex array for the subgrid from the subgrid start and count arrays

    ! ---------------------------------------
    ! Determine total number of global gridcells (including ocean)
    ! ---------------------------------------

    ng = nglob_x * nglob_y
    allocate(array_glob(ng))
    allocate(ioff(begg:endg)); ioff(:) = 0

    ! ---------------------------------------
    ! Gridcell gindex (compressed, no ocean points)
    ! ---------------------------------------

    allocate(gstart(begg:endg)) ; gstart(:) = 0
    call set_subgrid_start(gsmap_global, array_glob, gcount, gstart)
    allocate(gindex_grc(endg-begg+1))
    i = 0
    do gi = begg,endg
       if (gcount(gi) <  1) then
          write(iulog,*) 'decompInit_glcp warning count g ',iam,gi,gcount(gi)
       endif
       do l = 1,gcount(gi)
          i = i + 1
          gindex_grc(i) = gstart(gi) + l - 1
       enddo
    enddo
    deallocate(gstart)
    deallocate(gcount)

    ! ---------------------------------------
    ! Landunit gindex
    ! ---------------------------------------

    ! lstart for gridcell (n) is the total number of the landunits over gridcells 1->n-1
    allocate(lstart(begg:endg)) ; lstart(:) = 0
    call set_subgrid_start(gsmap_global, array_glob, lcount, lstart)
    allocate(gindex_lun(endl-begl+1))
    ioff(:) = 0
    do li = begl,endl
       !this is determined internally from how landunits are spread
       !out in memory
       gi = lun%gridcell(li)
       ! the output gindex is ALWAYS the same regardless of how
       ! landuntis are spread out in memory
       gindex_lun(li-begl+1) = lstart(gi) + ioff(gi)
       ioff(gi)  = ioff(gi) + 1
       ! check that this is less than [lstart(gi) + lcount(gi)]
    enddo
    deallocate(lstart)
    deallocate(lcount)

    ! ---------------------------------------
    ! Column gindex
    ! ---------------------------------------

    allocate(cstart(begg:endg)) ; cstart(:) = 0
    call set_subgrid_start(gsmap_global, array_glob, ccount, cstart)
    allocate(gindex_col(endc-begc+1))
    ioff(:) = 0
    do ci = begc,endc
       gi = col%gridcell(ci)
       gindex_col(ci-begc+1) = cstart(gi) + ioff(gi)
       ioff(gi) = ioff(gi) + 1
       ! check that this is less than [cstart(gi) + ccount(gi)]
    enddo
    deallocate(cstart)
    deallocate(ccount)

    ! ---------------------------------------
    ! PATCH gindex
    ! ---------------------------------------

    allocate(pstart(begg:endg)) ; pstart(:) = 0
    call set_subgrid_start(gsmap_global, array_glob, pcount, pstart)
    allocate(gindex_patch(endp-begp+1))
    ioff(:) = 0
    do pi = begp,endp
       gi = patch%gridcell(pi)
       gindex_patch(pi-begp+1) = pstart(gi) + ioff(gi)
       ioff(gi) = ioff(gi) + 1
       ! check that this is less than [pstart(gi) + pcount(gi)]
    enddo
    deallocate(pcount)
    deallocate(pstart)

    ! ---------------------------------------
    ! FATES gindex for the cohort/element vector
    ! ---------------------------------------

    if ( use_fates ) then
       allocate(coStart(begg:endg)); coStart(:) = 0
       call set_subgrid_start(gsmap_global, array_glob, coCount, coStart)
       allocate(gindex_cohort(endCohort-begCohort+1))
       ioff(:) = 0
       gi = begg
       do coi = begCohort,endCohort
          gindex_cohort(coi-begCohort+1) = coStart(gi) + ioff(gi)
          ioff(gi) = ioff(gi) + 1
          if ( mod(coi, fates_maxElementsPerSite ) == 0 ) gi = gi + 1
       enddo
       deallocate(coStart)
       deallocate(coCount)
    endif

    ! ---------------------------------------
    ! Deallocate memory and diagnostic output
    ! ---------------------------------------

    ! Deallocate start/count arrays
    deallocate(array_glob)
    deallocate(ioff)

    ! Diagnostic output
    if (masterproc) then
       write(iulog,*)' Surface Grid Characteristics'
       write(iulog,*)'   longitude points          = ',lni
       write(iulog,*)'   latitude points           = ',lnj
       write(iulog,*)'   total number of gridcells = ',numg
       write(iulog,*)'   total number of landunits = ',numl
       write(iulog,*)'   total number of columns   = ',numc
       write(iulog,*)'   total number of patches   = ',nump
       write(iulog,*)'   total number of cohorts   = ',numCohort
       write(iulog,*)' Decomposition Characteristics'
       write(iulog,*)'   clumps per process        = ',clump_pproc
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
          write(iulog,*)'proc= ',pid,&
               ' beg gridcell= ',procinfo%begg, &
               ' end gridcell= ',procinfo%endg,                   &
               ' total gridcells per proc= ',procinfo%ncells
          write(iulog,*)'proc= ',pid,&
               ' beg landunit= ',procinfo%begl, &
               ' end landunit= ',procinfo%endl,                   &
               ' total landunits per proc= ',procinfo%nlunits
          write(iulog,*)'proc= ',pid,&
               ' beg column  = ',procinfo%begc, &
               ' end column  = ',procinfo%endc,                   &
               ' total columns per proc  = ',procinfo%ncols
          write(iulog,*)'proc= ',pid,&
               ' beg patch     = ',procinfo%begp, &
               ' end patch     = ',procinfo%endp,                   &
               ' total patches per proc = ',procinfo%npatches
          write(iulog,*)'proc= ',pid,&
               ' beg coh     = ',procinfo%begCohort, &
               ' end coh     = ',procinfo%endCohort,                   &
               ' total coh per proc     = ',procinfo%nCohorts
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
                  ' beg patch     = ',clumps(cid)%begp, &
                  ' end patch     = ',clumps(cid)%endp, &
                  ' total patches per clump = ',clumps(cid)%npatches
             write(iulog,*)'proc= ',pid,' clump no = ',n, &
                  ' clump id= ',procinfo%cid(n),    &
                  ' beg cohort     = ',clumps(cid)%begCohort, &
                  ' end cohort     = ',clumps(cid)%endCohort, &
                  ' total cohorts per clump     = ',clumps(cid)%nCohorts
          end do
       end if
       call shr_sys_flush(iulog)
       call mpi_barrier(mpicom,ier)
    end do
    call shr_sys_flush(iulog)

  end subroutine decompInit_glcp

  !------------------------------------------------------------------------------
  subroutine set_subgrid_start(gsmap, array_glob, count, start)
    ! !USES:
    !
    use mct_mod , only : mct_aVect, mct_gsMap
    use mct_mod , only : mct_aVect_init, mct_aVect_importIattr, mct_aVect_scatter
    use mct_mod , only : mct_aVect_gather, mct_aVect_exportIattr, mct_aVect_clean
    use mct_mod , only : mct_aVect_exportRattr, mct_aVect_importRattr

    ! !ARGUMENTS:
    type(mct_gsmap)  :: gsmap         ! global gsmap
    integer, pointer :: array_glob(:) ! input
    integer, pointer :: count(:)      ! input
    integer, pointer :: start(:)      ! output

    ! !LOCAL VARIABLES:
    integer            :: n,lb,ub    ! indices
    integer            :: lsize      ! size of local array
    type(mct_aVect)    :: AVi, AVo   ! attribute vectors
    integer ,pointer   :: adata(:)   ! temporary data array
    integer            :: val1, val2 ! temporaries
    !-----------------------------------------------------------------------

    ! Initialize array_glob
    array_glob(:) = 0

    ! Gather count to master and place the data in arrray_global
    lsize = size(count, dim=1)
    lb = lbound(count, dim=1); ub = ubound(count, dim=1)
    call mct_aVect_init(AVi, rList="", iList='f1', lsize=lsize)
    allocate(adata(lsize))
    do n = lb,ub
       adata(n-lb+1) = count(n)
    enddo
    call mct_aVect_importIattr(AVi, 'f1', adata, lsize)
    deallocate(adata)
    call mct_aVect_gather(AVi, AVo, gsmap, 0, mpicom)
    if (masterproc) then
       lsize = size(array_glob,dim=1)
       call mct_aVect_exportIattr(AVo, 'f1', array_glob, lsize)
       call mct_aVect_clean(AVo)
    endif
    call mct_aVect_clean(AVi)

    ! Create global start array in array_glob
    if (masterproc) then
       val1 = array_glob(1)
       array_glob(1) = 1
       do n = 2,size(array_glob, dim=1)
          val2 = array_glob(n)
          array_glob(n) = array_glob(n-1) + val1
          val1 = val2
       enddo
    endif

    ! Now scatter start array (i.e. array_glob from master)
    if (masterproc) then
       call mct_aVect_init(AVi, rList="", iList="f1", lsize=lsize)
       call mct_aVect_importIattr(AVi, 'f1', array_glob, size(array_glob,dim=1))
    endif
    call mct_aVect_scatter(AVi, AVo, gsmap, 0, mpicom)
    lsize = size(start, dim=1)
    allocate(adata(lsize))
    call mct_aVect_exportIattr(AVo, 'f1', adata, lsize)
    lb = lbound(start, dim=1); ub = ubound(start, dim=1)
    do n = lb,ub
       start(n) = adata(n-lb+1)
    enddo
    deallocate(adata)
    if (masterproc) then
       call mct_aVect_clean(AVi)
    endif
    call mct_aVect_clean(AVo)

  end subroutine set_subgrid_start

  !------------------------------------------------------------------------------
  subroutine set_gsmap_global(gindex_global)
    !
    ! !USES:
    use spmdMod   , only : mpicom, comp_id
    use decompMod , only : nglob_x, nglob_y
    use mct_mod   , only : mct_gsMap_init
    !
    ! !ARGUMENTS:
    integer, intent(in) :: gindex_global(:)
    !
    ! !LOCAL VARIABLES:
    integer :: lsize, gsize
    !-----------------------------------------------------------------------

    lsize = size(gindex_global)
    gsize = nglob_x * nglob_y

    call mct_gsMap_init(gsmap_global, gindex_global, mpicom, comp_id, lsize, gsize)

  end subroutine set_gsmap_global

end module decompInitMod
