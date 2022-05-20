module decompInitMod

  !------------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module provides a descomposition into a clumped data structure which can
  ! be mapped back to atmosphere physics chunks.
  !
  ! !USES:
  use shr_kind_mod , only : r8 => shr_kind_r8
  use shr_sys_mod  , only : shr_sys_flush
  use shr_log_mod  , only : errMsg => shr_log_errMsg
  use spmdMod      , only : masterproc, iam, npes, mpicom
  use abortutils   , only : endrun
  use clm_varctl   , only : iulog
  !
  implicit none
  private
  !
  ! !PUBLIC TYPES:
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: decompInit_lnd    ! initializes lnd grid decomposition into clumps and processors
  public :: decompInit_clumps ! initializes atm grid decomposition into clumps
  public :: decompInit_glcp   ! initializes g,l,c,p decomp info
  !
  ! !PRIVATE MEMBER FUNCTIONS:
  !
  ! PUBLIC TYPES:
  integer, public :: clump_pproc ! number of clumps per MPI process
  !
  ! !PRIVATE TYPES:
  integer, pointer   :: lcid(:)          ! temporary for setting decomposition
  integer            :: nglob_x, nglob_y ! global sizes
  integer, parameter :: dbug=0           ! 0 = min, 1=normal, 2=much, 3=max
  character(len=*), parameter :: sourcefile = &
       __FILE__

#include <mpif.h>         ! mpi library include file
  !------------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------------
  subroutine decompInit_lnd(lni, lnj, amask)
    !
    ! !DESCRIPTION:
    ! This subroutine initializes the land surface decomposition into a clump
    ! data structure.  This assumes each pe has the same number of clumps
    ! set by clump_pproc
    !
    ! !USES:
    use clm_varctl , only : nsegspc
    use decompMod  , only : gindex_global, nclumps, clumps
    use decompMod  , only : bounds_type, get_proc_bounds, procinfo
    !
    ! !ARGUMENTS:
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
    ! Note that gindex_global goes from (1:endg)
    nglob_x = lni !  decompMod module variables
    nglob_y = lnj !  decompMod module variables
    call get_proc_bounds(bounds)
    allocate(gindex_global(1:bounds%endg))
    do n = procinfo%begg,procinfo%endg
       gindex_global(n-procinfo%begg+1) = gdc2glo(n)
    enddo

    deallocate(clumpcnt)
    deallocate(gdc2glo)

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
    use subgridMod     , only : subgrid_get_gcellinfo
    use decompMod      , only : bounds_type, clumps, nclumps, procinfo
    use decompMod      , only : get_proc_global, get_proc_bounds
    use decompMod      , only : numg, numl, numc, nump, numCohort
    use decompMod      , only : gindex_global
    use glcBehaviorMod , only : glc_behavior_type
    !
    ! !ARGUMENTS:
    integer                 , intent(in) :: lni,lnj ! land domain global size
    type(glc_behavior_type) , intent(in) :: glc_behavior
    !
    ! !LOCAL VARIABLES:
    integer           :: ln,an             ! indices
    integer           :: i,g,l,k           ! indices
    integer           :: cid,pid           ! indices
    integer           :: n,m,np            ! indices
    integer           :: anumg             ! lnd num gridcells
    integer           :: icells            ! temporary
    integer           :: begg, endg        ! temporary
    integer           :: ilunits           ! temporary
    integer           :: icols             ! temporary
    integer           :: ipatches          ! temporary
    integer           :: icohorts          ! temporary
    integer           :: ier               ! error code
    integer           :: npmin,npmax,npint ! do loop values for printing
    integer           :: clmin,clmax       ! do loop values for printing
    type(bounds_type) :: bounds            ! bounds
    integer, allocatable :: allvecg(:,:)   ! temporary vector "global"
    integer, allocatable :: allvecl(:,:)  ! temporary vector "local"
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

    ! Diagnostic output

    call get_proc_global(ng=numg, nl=numl, nc=numc, np=nump, nCohorts=numCohort)
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
          write(iulog,'(4(a,2x,i10))')'proc = ',pid,                                        &
               ' beg gridcell= ',procinfo%begg,' end gridcell= ',procinfo%endg,             &
               ' gridcells per proc = ',procinfo%ncells
          write(iulog,'(4(a,2x,i10))')'proc = ',pid,                                        &
               ' beg landunit= ',procinfo%begl,' end landunit= ',procinfo%endl,             &
               ' landunits per proc = ',procinfo%nlunits
          write(iulog,'(4(a,2x,i10))')'proc = ',pid,                                        &
               ' beg column  = ',procinfo%begc,' end column  = ',procinfo%endc,           &
               ' columns per proc = ',procinfo%ncols
          write(iulog,'(4(a,2x,i10))')'proc = ',pid,                                        &
               ' beg patch   = ',procinfo%begp,' end patch   = ',procinfo%endp,           &
               ' patches per proc = ',procinfo%npatches
          write(iulog,'(4(a,2x,i10))')'proc = ',pid,                                        &
               ' beg cohort  = ',procinfo%begCohort,' end cohort  = ',procinfo%endCohort, &
               ' coh per proc = ',procinfo%nCohorts
          write(iulog,'(2(a,2x,i10))')'proc = ',pid,' nclumps = ',procinfo%nclumps
          if (dbug == 0) then
             clmax = -1
          else
             clmax = procinfo%nclumps 
          endif
          do n = 1,clmax
             cid = procinfo%cid(n)
             write(iulog,'(6(a,2x,i10))')'proc = ',pid,' clump no = ',n,                             &
                  ' clump id= ',procinfo%cid(n),                                                     &
                  ' beg gridcell= ',clumps(cid)%begg,' end gridcell= ',clumps(cid)%endg,             &
                  ' gridcells per clump= ',clumps(cid)%ncells
             write(iulog,'(6(a,2x,i10))')'proc = ',pid,' clump no = ',n,                             &
                  ' clump id= ',procinfo%cid(n),                                                     &
                  ' beg landunit= ',clumps(cid)%begl,' end landunit= ',clumps(cid)%endl,             &
                  ' landunits per clump = ',clumps(cid)%nlunits
             write(iulog,'(6(a,2x,i10))')'proc = ',pid,' clump no = ',n,                             &
                  ' clump id= ',procinfo%cid(n),                                                     &
                  ' beg column  = ',clumps(cid)%begc,' end column  = ',clumps(cid)%endc,           &
                  ' columns per clump = ',clumps(cid)%ncols
             write(iulog,'(6(a,2x,i10))')'proc = ',pid,' clump no = ',n,                             &
                  ' clump id= ',procinfo%cid(n),                                                     &
                  ' beg patch   = ',clumps(cid)%begp,' end patch   = ',clumps(cid)%endp,           &
                  ' patches per clump = ',clumps(cid)%npatches
             write(iulog,'(6(a,2x,i10))')'proc = ',pid,' clump no = ',n,                             &
                  ' clump id= ',procinfo%cid(n),                                                     &
                  ' beg cohort  = ',clumps(cid)%begCohort,' end cohort  = ',clumps(cid)%endCohort, &
                  ' cohorts per clump = ',clumps(cid)%nCohorts

          end do
       end if
       call shr_sys_flush(iulog)
       call mpi_barrier(mpicom,ier)
    end do

  end subroutine decompInit_clumps

  !------------------------------------------------------------------------------
  subroutine decompInit_glcp(lni,lnj,glc_behavior)
    !
    ! !DESCRIPTION:
    ! Determine gindex for landunits, columns, patches and cohorts
    !
    ! !USES:
    use clm_varctl             , only : use_fates
    use subgridMod             , only : subgrid_get_gcellinfo
    use decompMod              , only : bounds_type, get_proc_global, get_proc_bounds
    use decompMod              , only : gindex_global
    use decompMod              , only : gindex_grc, gindex_lun, gindex_col, gindex_patch, gindex_Cohort
    use decompMod              , only : procinfo, clump_type, clumps, get_proc_global
    use LandunitType           , only : lun
    use ColumnType             , only : col
    use PatchType              , only : patch
    use FatesInterfaceTypesMod , only : fates_maxElementsPerSite
    use glcBehaviorMod         , only : glc_behavior_type
    !
    ! !ARGUMENTS:
    integer                 , intent(in) :: lni,lnj ! land domain global size
    type(glc_behavior_type) , intent(in) :: glc_behavior
    !
    ! !LOCAL VARIABLES:
    integer              :: gi,li,ci,pi,coi     ! indices
    integer              :: i,l,n,np            ! indices
    integer              :: cid,pid             ! indices
    integer              :: numg                ! total number of land gridcells across all processors
    integer              :: numl                ! total number of landunits across all processors
    integer              :: numc                ! total number of columns across all processors
    integer              :: nump                ! total number of patches across all processors
    integer              :: numCohort           ! fates cohorts
    integer              :: ilunits             ! temporary
    integer              :: icols               ! temporary
    integer              :: ipatches            ! temporary
    integer              :: icohorts            ! temporary
    integer              :: ier                 ! error code
    integer, pointer     :: gcount(:)
    integer, pointer     :: lcount(:)
    integer, pointer     :: ccount(:)
    integer, pointer     :: pcount(:)
    integer, pointer     :: coCount(:)
    type(bounds_type)    :: bounds
    integer, allocatable :: ioff(:)
    integer, allocatable :: gridcells_per_pe(:) ! needed for gindex at all levels
    integer, allocatable :: gridcell_offsets(:) ! needed for gindex at all levels
    integer, allocatable :: index_gridcells(:)  ! needed for gindex at all levels
    integer, allocatable :: start_global(:)
    integer, allocatable :: start(:) 
    integer, allocatable :: index_lndgridcells(:)
    integer              :: count
    integer              :: temp
    integer              :: lsize_g, lsize_l, lsize_c, lsize_p, lsize_cohort
    integer              :: gsize
    Character(len=32), parameter :: subname = 'decompInit_glcp'
    !------------------------------------------------------------------------------

    ! Get processor bounds

    call get_proc_bounds(bounds)
    call get_proc_global(ng=numg, nl=numl, nc=numc, np=nump, nCohorts=numCohort)

    lsize_g = bounds%endg
    lsize_l = bounds%endl
    lsize_c = bounds%endc
    lsize_p = bounds%endp
    lsize_cohort = bounds%endCohort
    gsize = nglob_x * nglob_y

    ! allocate module variables in decompMod.F90
    allocate(gindex_grc(lsize_g))
    allocate(gindex_lun(lsize_l))
    allocate(gindex_col(lsize_c))
    allocate(gindex_patch(lsize_p))
    allocate(gindex_cohort(lsize_cohort))

    ! Determine counts
    allocate(gcount(lsize_g))  ; gcount(:) = 0
    allocate(lcount(lsize_g))  ; lcount(:) = 0
    allocate(ccount(lsize_g))  ; ccount(:) = 0
    allocate(pcount(lsize_g))  ; pcount(:) = 0
    allocate(coCount(lsize_g)) ; coCount(:) = 0
    do gi = 1,lsize_g
       call subgrid_get_gcellinfo (gi, nlunits=ilunits, ncols=icols, npatches=ipatches, &
            ncohorts=icohorts, glc_behavior=glc_behavior)
       gcount(gi)  = 1         ! number of gridcells for local gridcell index gi
       lcount(gi)  = ilunits   ! number of landunits for local gridcell index gi
       ccount(gi)  = icols     ! number of columns for local gridcell index gi
       pcount(gi)  = ipatches  ! number of patches for local gridcell index gi
       coCount(gi) = icohorts  ! number of fates cohorts for local gricell index gi
    enddo

    ! ---------------------------------------
    ! Arrays needed to determine gindex_xxx(:)
    ! ---------------------------------------

    allocate(ioff(lsize_g))

    if (masterproc) then
       allocate (gridcells_per_pe(0:npes-1))
    else
       allocate(gridcells_per_pe(0))
    endif
    call mpi_gather(lsize_g, 1, MPI_INTEGER, gridcells_per_pe, 1, MPI_INTEGER, 0, mpicom, ier)

    if (masterproc) then
       allocate(gridcell_offsets(0:npes-1))
       gridcell_offsets(0) = 0
       do n = 1 ,npes-1
          gridcell_offsets(n) = gridcell_offsets(n-1) + gridcells_per_pe(n-1)
       end do
    else
       allocate(gridcell_offsets(0))
    end if

    if (masterproc) then
       allocate(start_global(numg)) ! number of landunits in a gridcell
    else
       allocate(start_global(0))
    end if

    allocate(start(lsize_g))

    ! ---------------------------------------
    ! Gridcell gindex (compressed, no ocean points)
    ! ---------------------------------------

    ! gstart_global the global index of all of the land points in task order 
    call mpi_gatherv(gindex_global, lsize_g, MPI_INTEGER, start_global, &
         gridcells_per_pe, gridcell_offsets, MPI_INTEGER, 0, mpicom, ier)

    if (masterproc) then
       ! Create a global size index_gridcells that will have 0 for all ocean points
       ! Fill the location of each land point with the gatherv location of that land point
       allocate(index_gridcells(gsize))
       index_gridcells(:) = 0
       do n = 1,numg
          ! if n = 3, start_global(3)=100, index_gridcells(100)=3
          ! n is the task order location - so for global index 100 - the task order location is 3
          index_gridcells(start_global(n)) = n
       end do

       ! Create a land-only global index based on the original global index ordering
       ! Count is the running global land index
       allocate(index_lndgridcells(numg))
       count = 0
       do n = 1,gsize
          if (index_gridcells(n) > 0) then
             count = count + 1
             ! e.g. n=20, count=4 and index_gridcells(20)=100, then start_global(100)=4
             start_global(index_gridcells(n)) = count  
             index_lndgridcells(count) = index_gridcells(n)
          end if
       end do
       deallocate(index_gridcells)
    end if

    ! Determine gindex_grc
    call mpi_scatterv(start_global, gridcells_per_pe, gridcell_offsets, MPI_INTEGER, gindex_grc, &
         lsize_g, MPI_INTEGER, 0, mpicom, ier)
    deallocate(gcount)

    ! ---------------------------------------
    ! Landunit gindex
    ! ---------------------------------------

    start(:) = 0 
    call mpi_gatherv(lcount, lsize_g, MPI_INTEGER, start_global, &
         gridcells_per_pe, gridcell_offsets, MPI_INTEGER, 0, mpicom, ier)
    if (masterproc) then
       count = 1
       do n = 1,numg
          temp = start_global(index_lndgridcells(n))
          start_global(index_lndgridcells(n)) = count
          count = count + temp
       end do
    endif
    call mpi_scatterv(start_global, gridcells_per_pe, gridcell_offsets, MPI_INTEGER, start, &
         lsize_g, MPI_INTEGER, 0, mpicom, ier)

    ioff(:) = 0
    do li = 1,lsize_l
       gi = lun%gridcell(li)
       gindex_lun(li) = start(gi) + ioff(gi)
       ioff(gi)  = ioff(gi) + 1
    enddo
    deallocate(lcount)

    ! ---------------------------------------
    ! Column gindex
    ! ---------------------------------------

    start(:) = 0
    call mpi_gatherv(ccount, lsize_g, MPI_INTEGER, start_global, &
         gridcells_per_pe, gridcell_offsets, MPI_INTEGER, 0, mpicom, ier)
    if (masterproc) then
       count = 1
       do n = 1,numg
          temp = start_global(index_lndgridcells(n))
          start_global(index_lndgridcells(n)) = count
          count = count + temp
       end do
    endif
    call mpi_scatterv(start_global, gridcells_per_pe, gridcell_offsets, MPI_INTEGER, start, &
         lsize_g, MPI_INTEGER, 0, mpicom, ier)

    ioff(:) = 0
    do ci = 1,lsize_c
       gi = col%gridcell(ci)
       gindex_col(ci) = start(gi) + ioff(gi)
       ioff(gi) = ioff(gi) + 1
    enddo
    deallocate(ccount)

    ! ---------------------------------------
    ! PATCH gindex
    ! ---------------------------------------

    start(:) = 0
    call mpi_gatherv(pcount, lsize_g, MPI_INTEGER, start_global, &
         gridcells_per_pe, gridcell_offsets, MPI_INTEGER, 0, mpicom, ier)
    if (masterproc) then
       count = 1
       do n = 1,numg
          temp = start_global(index_lndgridcells(n))
          start_global(index_lndgridcells(n)) = count
          count = count + temp
       end do
    endif
    call mpi_scatterv(start_global, gridcells_per_pe, gridcell_offsets, MPI_INTEGER, start, &
         lsize_g, MPI_INTEGER, 0, mpicom, ier)

    ioff(:) = 0
    do pi = 1,lsize_p
       gi = patch%gridcell(pi)
       gindex_patch(pi) = start(gi) + ioff(gi)
       ioff(gi) = ioff(gi) + 1
    enddo
    deallocate(pcount)

    ! ---------------------------------------
    ! FATES gindex for the cohort/element vector
    ! ---------------------------------------

    if ( use_fates ) then
       start(:) = 0
       call mpi_gatherv(coCount, lsize_g, MPI_INTEGER, start_global, &
            gridcells_per_pe, gridcell_offsets, MPI_INTEGER, 0, mpicom, ier)
       if (masterproc) then
          count = 1
          do n = 1,numg
             temp = start_global(index_lndgridcells(n))
             start_global(index_lndgridcells(n)) = count
             count = count + temp
          end do
       endif
       call mpi_scatterv(start_global, gridcells_per_pe, gridcell_offsets, MPI_INTEGER, start, &
            lsize_g, MPI_INTEGER, 0, mpicom, ier)

       ioff(:) = 0
       gi = 1
       do coi = 1, lsize_cohort
          gindex_cohort(coi) = start(gi) + ioff(gi)
          ioff(gi) = ioff(gi) + 1
          if ( mod(coi, fates_maxElementsPerSite ) == 0 ) then
             gi = gi + 1
          end if
       enddo
       deallocate(coCount)
    endif

    ! ---------------------------------------
    ! Deallocate memory
    ! ---------------------------------------

    deallocate(ioff)
    deallocate(gridcells_per_pe)
    deallocate(gridcell_offsets)
    deallocate(start)
    deallocate(start_global)
    if (allocated(index_lndgridcells)) deallocate(index_lndgridcells)

  end subroutine decompInit_glcp

end module decompInitMod
