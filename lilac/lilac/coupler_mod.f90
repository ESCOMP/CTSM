    module CouplerMod

    use ESMF

    implicit none

    private

    ! Public entry point
    public Coupler_register

    contains


!------------------------------------------------------------------------------
!BOPI
! !IROUTINE: Coupler_register - public SetServices entry point

! !INTERFACE:
     subroutine Coupler_register(comp, rc)
!
! !ARGUMENTS:
     type(ESMF_CplComp)   :: comp
     integer, intent(out) :: rc
!
! !DESCRIPTION:
!     User-supplied setservices routine.
!
!     The arguments are:
!     \begin{description}
!     \item[comp]
!          Component.
!     \item[rc]
!          Return code; equals {\tt ESMF\_SUCCESS} if there are no errors,
!          otherwise {\tt ESMF\_FAILURE}.
!     \end{description}
!
!EOPI

      ! because none of the arguments to this subroutine will ever be optional,
      ! go ahead and set rc to an initial return code before using it below.
      ! (this makes some eager error-checking compilers happy.)
      rc = ESMF_FAILURE

      ! Register the callback routines.

      call ESMF_CplCompSetEntryPoint(comp, ESMF_METHOD_INITIALIZE, userRoutine=coupler_init, rc=rc)
      if(rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT, rc=rc)
      call ESMF_CplCompSetEntryPoint(comp, ESMF_METHOD_RUN, userRoutine=coupler_run, rc=rc)
      if(rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT, rc=rc)
      call ESMF_CplCompSetEntryPoint(comp, ESMF_METHOD_FINALIZE, userRoutine=coupler_final, rc=rc)
      if(rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT, rc=rc)

      print *, "CouplerMod: Registered Initialize, Run, and Finalize routines"

    end subroutine


!------------------------------------------------------------------------------
!BOPI
! !IROUTINE: coupler_init - coupler init routine

! !INTERFACE:
      subroutine coupler_init(comp, importState, exportState, clock, rc)

!
! !ARGUMENTS:
      type(ESMF_CplComp)   :: comp
      type(ESMF_State)     :: importState, exportState
      type(ESMF_Clock)     :: clock
      integer, intent(out) :: rc
!
! !DESCRIPTION:
!     User-supplied init routine.
!
!     The arguments are:
!     \begin{description}
!     \item[comp]
!          Component.
!     \item[importState]
!          Nested state object containing import data.
!     \item[exportState]
!          Nested state object containing export data.
!     \item[clock]
!          External clock.
!     \item[rc]
!          Return code; equals {\tt ESMF\_SUCCESS} if there are no errors,
!          otherwise {\tt ESMF\_FAILURE}.
!     \end{description}
!
!EOPI

!   ! Local variables
    type(ESMF_Field) :: src_field, dst_field
    type(ESMF_VM) :: vm
    character(ESMF_MAXSTR) :: statename

    print *, "Coupler Init starting"

    ! because none of the arguments to this subroutine will ever be optional,
    ! go ahead and set rc to an initial return code before using it below.
    ! (this makes some eager error-checking compilers happy.)
    rc = ESMF_FAILURE

    ! Get VM from coupler component to use in computing redistribution
    call ESMF_CplCompGet(comp, vm=vm, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT, rc=rc)

    ! Use placeholder SIE
    call ESMF_StateGet(importState, name=statename, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT, rc=rc)
    call ESMF_StateGet(importState, "SIE", src_field, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT, rc=rc)
    call ESMF_StateGet(exportState, "SIE", dst_field, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT, rc=rc)

    ! Compute routehandle
    ! Since state items are needed by default, mark Fields not needed during coupling
    if (trim(statename) .eq. "FlowSolver Feedback") then
      call setFieldNeeded(importState, "U", .false., rc=rc)
      if(rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT, rc=rc)
      call setFieldNeeded(importState, "P", .false., rc=rc)
      if(rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT, rc=rc)
      call setFieldNeeded(importState, "Q", .false., rc=rc)
      if(rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT, rc=rc)

      call ESMF_FieldRedistStore(src_field, dst_field, &
                                 routehandle=fromFlow_rh, rc=rc)
      if(rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT, rc=rc)

    endif

    if (trim(statename) .eq. "Injection Feedback") then
      call setFieldNeeded(importState, "U", .false., rc=rc)
      if(rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT, rc=rc)
      call setFieldNeeded(importState, "P", .false., rc=rc)
      if(rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT, rc=rc)
      call setFieldNeeded(importState, "Q", .false., rc=rc)
      if(rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT, rc=rc)
      call setFieldNeeded(importState, "FLAG", .false., rc=rc)
      if(rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT, rc=rc)

      call ESMF_FieldRedistStore(src_field, dst_field, &
                                 routehandle=fromInject_rh, rc=rc)
      if(rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT, rc=rc)

    endif

    print *, "Coupler Init returning"

    end subroutine coupler_init


!------------------------------------------------------------------------------
!BOPI
! !IROUTINE: coupler_run - coupler run routine

! !INTERFACE:
      subroutine coupler_run(comp, importState, exportState, clock, rc)

!
! !ARGUMENTS:
     type(ESMF_CplComp)   :: comp
     type(ESMF_State)     :: importState, exportState
     type(ESMF_Clock)     :: clock
     integer, intent(out) :: rc
!
! !DESCRIPTION:
!     User-supplied run routine.
!
!     The arguments are:
!     \begin{description}
!     \item[comp]
!          Component.
!     \item[importState]
!          Nested state object containing import data.
!     \item[exportState]
!          Nested state object containing export data.
!     \item[clock]
!          External clock.
!     \item[rc]
!          Return code; equals {\tt ESMF\_SUCCESS} if there are no errors,
!          otherwise {\tt ESMF\_FAILURE}.
!     \end{description}
!
!EOPI

      ! Local variables
        type(ESMF_Field) :: srcfield, dstfield
        type(ESMF_RouteHandle) :: rh

        character(len=ESMF_MAXSTR) :: statename

        integer :: i, datacount
        character(len=ESMF_MAXSTR), dimension(7) :: datanames

        ! none of the arguments to this subroutine will ever be optional, so
        ! go ahead and set rc to an initial return code before using it below.
        ! (this makes some eager error-checking compilers happy.)
        rc = ESMF_FAILURE

        datacount = 7
        datanames(1) = "SIE"
        datanames(2) = "U"
        datanames(3) = "V"
        datanames(4) = "RHO"
        datanames(5) = "P"
        datanames(6) = "Q"
        datanames(7) = "FLAG"

        ! In this case, the coupling is symmetric - you call redist going
        ! both ways - so we only care about the coupling direction in order
        ! to get the right routehandle selected.
        call ESMF_StateGet(importState, name=statename, rc=rc)
        if(rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT, rc=rc)
        if (trim(statename) .eq. "FlowSolver Feedback") then
            rh = fromFlow_rh
        else
            rh = fromInject_rh
        endif

        do i=1, datacount

           ! check isneeded flag here
           if (.not. isFieldNeeded(importState, datanames(i), rc=rc)) then
               !print *, "skipping field ", trim(datanames(i)), " not needed"
               cycle
           endif

           !print *, "processing field ", trim(datanames(i)), " as needed"
!BOE
! !DESCRIPTION:
! \subsubsection{Example of Redist Usage}
!
!   The following piece of code provides an example of calling the data
!   redistribution routine  between two Fields in the Coupler Component.
!   Unlike regrid, which translates between
!   different Grids, redist translates between different DELayouts on
!   the same Grid.   The first two lines get the Fields from the
!   States, each corresponding to a different subcomponent.  One is
!   an Export State and the other is an Import State.
!
!BOC
           call ESMF_StateGet(importState, datanames(i), srcfield, rc=rc)
           if(rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT, rc=rc)
           call ESMF_StateGet(exportState, datanames(i), dstfield, rc=rc)
           if(rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT, rc=rc)
!EOC
!
!   The redist routine uses information contained in the Fields and the
!   Coupler VM object to call the communication routines to move the data.
!   Because many Fields may share the same Grid association, the same
!   routing information may be needed repeatedly.  Route information is
!   saved so the precomputed information can be retained.  The following
!   is an example of a Field redist call:
!
!BOC
           call ESMF_FieldRedist(srcfield, dstfield, routehandle=rh, rc=rc)
           if(rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT, rc=rc)

!EOC
!EOE

        enddo

        ! rc has the last error code already

    end subroutine coupler_run


!------------------------------------------------------------------------------
!BOPI
! !IROUTINE:  coupler_final - finalization routine

! !INTERFACE:
      subroutine coupler_final(comp, importState, exportState, clock, rc)

!
! !ARGUMENTS:
      type(ESMF_CplComp)   :: comp
      type(ESMF_State)     :: importState, exportState
      type(ESMF_Clock)     :: clock
      integer, intent(out) :: rc
!
! !DESCRIPTION:
!     User-supplied finalize routine.
!
!     The arguments are:
!     \begin{description}
!     \item[comp]
!          Component.
!     \item[importState]
!          Nested state object containing import data.
!     \item[exportState]
!          Nested state object containing export data.
!     \item[clock]
!          External clock.
!     \item[rc]
!          Return code; equals {\tt ESMF\_SUCCESS} if there are no errors,
!          otherwise {\tt ESMF\_FAILURE}.
!     \end{description}
!
!EOPI

        print *, "Coupler Final starting"

        ! none of the arguments to this subroutine will ever be optional, so
        ! go ahead and set rc to an initial return code before using it below.
        ! (this makes some eager error-checking compilers happy.)
        rc = ESMF_FAILURE

        ! Only thing to do here is release redist and route handles
        call ESMF_FieldRedistRelease(fromFlow_rh, rc=rc)
        if(rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT, rc=rc)

        call ESMF_FieldRedistRelease(fromInject_rh, rc=rc)
        if(rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT, rc=rc)

        rc = ESMF_SUCCESS

        print *, "Coupler Final returning"

    end subroutine coupler_final


    end module CouplerMod
