module cpl_mod

    !-----------------------------------------------------------------------
    ! !DESCRIPTION:
    ! Module containing all routines for both couplers
    ! 1- coupler 1 : atm ---> lnd (cpl_atm2lnd)
    ! 2- coupler 2 : lnd ---> atm (cpl_lnd2atm)
    !-----------------------------------------------------------------------
    ! !USES
    use ESMF
     use clm_varctl    ,  only : iulog
    implicit none

    include 'mpif.h'

    private


    public                       :: cpl_atm2lnd_register
    public                       :: cpl_lnd2atm_register

    character(*), parameter      :: modname =  "  cpl_mod"
    type(ESMF_RouteHandle), save :: rh_atm2lnd, rh_lnd2atm


    integer :: mpierror, numprocs
    integer :: i, myid
    integer status(MPI_STATUS_SIZE)

    character(len=128)                                    :: fldname
    integer, parameter     :: begc = 1   !-- internal debug level
    integer, parameter     :: endc = 3312/4/2/2   !-- internal debug level
    character(*),parameter :: F01 = "('[cpl_mod] ',a,i5,2x,i5,2x,d21.14)"
    character(*),parameter :: F02 =  "('[cpl_mod]',a,i5,2x,d26.19)"
    integer, parameter                              :: debug = 1   !-- internaldebug level
    !======================================================================
     contains
    !======================================================================

    subroutine cpl_atm2lnd_register(cplcomp, rc)
        type(ESMF_CplComp   )             :: cplcomp
        integer, intent(out )             :: rc
        character(len=*     ) , parameter :: subname=trim(modname ) //' : [cpl_atm2lnd_register] '

        rc = ESMF_SUCCESS
        print *, "in cpl_atm2lnd_register routine"

        ! Register the callback routines.
        ! Set the entry points for coupler ESMF Component methods
        call ESMF_CplCompSetEntryPoint(cplcomp, ESMF_METHOD_INITIALIZE, userRoutine= cpl_atm2lnd_init, rc=rc)
        if(rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT, rc=rc)

        call ESMF_CplCompSetEntryPoint(cplcomp, ESMF_METHOD_RUN       , userRoutine=cpl_atm2lnd_run  , rc=rc)
        if(rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT, rc=rc)

        call ESMF_CplCompSetEntryPoint(cplcomp, ESMF_METHOD_FINALIZE  , userRoutine=cpl_atm2lnd_final, rc=rc)
        if(rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT, rc=rc)
    end subroutine cpl_atm2lnd_register

    subroutine cpl_lnd2atm_register(cplcomp, rc)
        type(ESMF_CplComp   )             :: cplcomp
        integer, intent(out )             :: rc
        character(len=*     ) , parameter :: subname=trim(modname ) //' : [cpl_lnd2atm_register] '


        rc = ESMF_SUCCESS
        print *, "in cpl_lnd2atm_register routine"

        ! Register the callback routines.
        ! Set the entry points for coupler ESMF Component methods
        call ESMF_CplCompSetEntryPoint(cplcomp, ESMF_METHOD_INITIALIZE, cpl_lnd2atm_init, rc=rc)
        if(rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT, rc=rc)

        call ESMF_CplCompSetEntryPoint(cplcomp, ESMF_METHOD_RUN       , userRoutine=cpl_lnd2atm_run  , rc=rc)
        if(rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT, rc=rc)

        call ESMF_CplCompSetEntryPoint(cplcomp, ESMF_METHOD_FINALIZE  , userRoutine=cpl_lnd2atm_final, rc=rc)
        if(rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT, rc=rc)
    end subroutine cpl_lnd2atm_register

    !--------------------------------------------------------------------------
    ! couplers init....
    !--------------------------------------------------------------------------

    subroutine cpl_atm2lnd_init(cplcomp, importState, exportState, clock, rc)

        type (ESMF_CplComp    )             :: cplcomp
        type (ESMF_State      )             :: importState
        type (ESMF_State      )             :: exportState
        type (ESMF_Clock      )             :: clock
        type (ESMF_FieldBundle)             :: import_fieldbundle, export_fieldbundle
        integer, intent(out   )             :: rc
        character(len=*       ) , parameter :: subname=trim(modname) //': [cpl_atm2lnd_init] '

        rc = ESMF_SUCCESS
        print *, "Coupler for atmosphere to land initialize routine called"
        call ESMF_LogWrite(subname//"-----------------!", ESMF_LOGMSG_INFO)

        call MPI_Comm_size(MPI_COMM_WORLD, numprocs, mpierror)
        call MPI_Comm_rank(MPI_COMM_WORLD, myid, mpierror)



        call ESMF_StateGet(importState, trim("a2c_fb"), import_fieldbundle, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out

        call ESMF_StateGet(exportState, trim("c2l_fb"), export_fieldbundle, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out


        if (myid == 0) then
            print *, "PRINTING FIELDBUNDLES"
            call ESMF_FieldBundlePrint (import_fieldbundle, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
            call ESMF_FieldBundlePrint (export_fieldbundle, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
        end if


        call ESMF_FieldBundleRedistStore(import_fieldbundle, export_fieldbundle, routehandle=rh_atm2lnd, rc=rc)
        !call ESMF_FieldBundleRegridStore(import_fieldbundle, export_fieldbundle, routehandle=rh_atm2lnd, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
        call ESMF_LogWrite(subname//"cpl init finished!", ESMF_LOGMSG_INFO)
    end subroutine cpl_atm2lnd_init

    subroutine cpl_lnd2atm_init(cplcomp, importState, exportState, clock, rc)

        type (ESMF_CplComp     )             :: cplcomp
        type (ESMF_State       )             :: importState
        type (ESMF_State       )             :: exportState
        type (ESMF_Clock       )             :: clock
        type (ESMF_FieldBundle )             :: import_fieldbundle, export_fieldbundle
        integer, intent(out    )             :: rc
        character(len=*        ) , parameter :: subname=trim(modname ) //': [cpl_lnd2atm_init] '

        rc = ESMF_SUCCESS
        print *, "Coupler for land to atmosphere initialize routine called"
        call ESMF_LogWrite(subname//"-----------------!", ESMF_LOGMSG_INFO)

        call ESMF_StateGet(importState, "l2c_fb", import_fieldbundle, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out

        call ESMF_StateGet(exportState, "c2a_fb", export_fieldbundle, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out

        call ESMF_FieldBundleRedistStore(import_fieldbundle, export_fieldbundle, routehandle=rh_lnd2atm, rc=rc)
        !call ESMF_FieldBundleRegridStore(import_fieldbundle, export_fieldbundle, routehandle=rh_lnd2atm, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
        call ESMF_LogWrite(subname//"cpl init finished!", ESMF_LOGMSG_INFO)
    end subroutine cpl_lnd2atm_init

    !--------------------------------------------------------------------------
    ! Couplers Run phase 
    !--------------------------------------------------------------------------

    subroutine cpl_atm2lnd_run(cplcomp, importState, exportState, clock, rc)

        type(ESMF_CplComp      )             :: cplcomp
        type(ESMF_State        )             :: importState
        type(ESMF_State        )             :: exportState
        type(ESMF_Clock        )             :: clock
        integer, intent(out    )             :: rc
        type (ESMF_FieldBundle )             :: import_fieldbundle, export_fieldbundle
        character(len=*        ) , parameter :: subname=trim(modname       ) //': [cpl_atm2lnd_run] '

        real, pointer           :: fldptr1d(:)

        rc = ESMF_SUCCESS
        print *, "Running cpl_atm2lnd_run"
        call ESMF_LogWrite(subname//"-----------------!", ESMF_LOGMSG_INFO)

        call ESMF_StateGet(importState, trim("a2c_fb"), import_fieldbundle, rc=rc)
        !call ESMF_StateGet(importState, itemName=trim("a2c_fb"), item=import_fieldbundle, rc=rc)     ! this syntax was not working???
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
        call ESMF_LogWrite(subname//" got a2c fieldbundle!", ESMF_LOGMSG_INFO)

        call ESMF_StateGet(exportState, trim("c2l_fb"), export_fieldbundle, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
        call ESMF_LogWrite(subname//" got c2l fieldbundle!", ESMF_LOGMSG_INFO)

        !fldname = 'Sa_topo'
        !call state_getfldptr(exportState, trim(fldname), fldptr1d=fldptr1d, rc=rc)

        !call ESMF_FieldBundleRegrid(import_fieldbundle, export_fieldbundle, routehandle=rh_atm2lnd, rc=rc)
        call ESMF_FieldBundleRedist(import_fieldbundle, export_fieldbundle, routehandle=rh_atm2lnd, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
        call ESMF_LogWrite(subname//" regridding fieldbundles from atmos to land!", ESMF_LOGMSG_INFO)

    end subroutine cpl_atm2lnd_run


    subroutine cpl_lnd2atm_run(cplcomp, importState, exportState, clock, rc)

        type(ESMF_CplComp      )             :: cplcomp
        type(ESMF_State        )             :: importState
        type(ESMF_State        )             :: exportState
        type(ESMF_Clock        )             :: clock
        integer, intent(out    )             :: rc
        type (ESMF_FieldBundle )             :: import_fieldbundle, export_fieldbundle
        character(len=*        ) , parameter :: subname=trim(modname       ) //': [cpl_lnd2atm_run] '

        rc = ESMF_SUCCESS
        print *, "Running cpl_lnd2atm_run"
        call ESMF_LogWrite(subname//"-----------------!", ESMF_LOGMSG_INFO)

        call ESMF_StateGet(importState, "l2c_fb", import_fieldbundle, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out

        call ESMF_StateGet(exportState, "c2a_fb", export_fieldbundle, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out

        call ESMF_FieldBundleRedist(import_fieldbundle, export_fieldbundle, routehandle=rh_lnd2atm, rc=rc)
        !call ESMF_FieldBundleRegrid(import_fieldbundle, export_fieldbundle, routehandle=rh_lnd2atm, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
        call ESMF_LogWrite(subname//" regridding fieldbundles  from land to atmos!", ESMF_LOGMSG_INFO)

    end subroutine cpl_lnd2atm_run

    !--------------------------------------------------------------------------
    ! couplers final phase
    !--------------------------------------------------------------------------

    subroutine cpl_atm2lnd_final(cplcomp, importState, exportState, clock, rc)

        type (ESMF_CplComp     )             :: cplcomp
        type (ESMF_State       )             :: importState
        type (ESMF_State       )             :: exportState
        type (ESMF_Clock       )             :: clock
        type (ESMF_FieldBundle )             :: import_fieldbundle, export_fieldbundle
        integer, intent(out    )             :: rc
        character(len=*        ) , parameter :: subname=trim(modname       ) //': [cpl_atm2lnd_final] '

        rc = ESMF_SUCCESS

        call ESMF_LogWrite(subname//"---------------------------------!", ESMF_LOGMSG_INFO)

        ! Only thing to do here is release redist (or regrid) and route handles
        call ESMF_FieldBundleRegridRelease (routehandle=rh_atm2lnd, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out

        call ESMF_LogWrite(subname//" rh_atm2lnd route handle released!", ESMF_LOGMSG_INFO)

    end subroutine cpl_atm2lnd_final

    subroutine cpl_lnd2atm_final(cplcomp, importState, exportState, clock, rc)

        type (ESMF_CplComp     )             :: cplcomp
        type (ESMF_State       )             :: importState
        type (ESMF_State       )             :: exportState
        type (ESMF_Clock       )             :: clock
        type (ESMF_FieldBundle )             :: import_fieldbundle, export_fieldbundle
        integer, intent(out    )             :: rc
        character(len=*        ) , parameter :: subname=trim(modname) //': [cpl_lnd2atm_final] '

        rc = ESMF_SUCCESS

        call ESMF_LogWrite(subname//"---------------------------------!", ESMF_LOGMSG_INFO)
        ! Only thing to do here is release redist (or regrid) and route handles
        call ESMF_FieldBundleRegridRelease (routehandle=rh_lnd2atm , rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out

        call ESMF_LogWrite(subname//" rh_lnd2atm route handle released!", ESMF_LOGMSG_INFO)

    end subroutine cpl_lnd2atm_final


     !===============================================================================

  subroutine state_getfldptr(State, fldname, fldptr1d, fldptr2d, rc)

    ! ----------------------------------------------
    ! Get pointer to a state field
    ! ----------------------------------------------

    use ESMF , only : ESMF_State, ESMF_Field, ESMF_Mesh, ESMF_FieldStatus_Flag
    use ESMF , only : ESMF_FieldBundle
    use ESMF , only : ESMF_StateGet, ESMF_FieldGet, ESMF_MeshGet
    use ESMF , only : ESMF_FIELDSTATUS_COMPLETE, ESMF_FAILURE
    use ESMF                  , only : ESMF_FieldBundleGet

    ! input/output variables
    type(ESMF_State),             intent(in)    :: State
    character(len=*),             intent(in)    :: fldname
    real , pointer, optional , intent(out)   :: fldptr1d(:)
    real , pointer, optional , intent(out)   :: fldptr2d(:,:)
    integer,                      intent(out)   :: rc

    ! local variables
    type(ESMF_FieldStatus_Flag) :: status
    type(ESMF_Field)            :: lfield
    type(ESMF_Mesh)             :: lmesh
    integer                     :: nnodes, nelements
    character(len=*), parameter :: subname='(lnd_import_export:state_getfldptr)'

    type(ESMF_StateItem_Flag)   :: itemFlag
    type(ESMF_FieldBundle)      :: fieldBundle
    logical                       :: isPresent
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    ! Determine if this field bundle exist....
    ! TODO: combine the error checks....


    call ESMF_StateGet(state, "c2l_fb", itemFlag, rc=rc)
    !call ESMF_StateGet(State, itemName=trim(fldname), field=lfield, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out

    ! Get the fieldbundle from state...
    call ESMF_StateGet(state, "c2l_fb", fieldBundle, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out


    call ESMF_FieldBundleGet(fieldBundle,fieldName=trim(fldname), field=lfield,  isPresent=isPresent, rc=rc)
    !call ESMF_FieldBundleGet(fieldBundle,trim(fldname), lfield,  isPresent, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out

    call ESMF_FieldGet(lfield, status=status, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out

    if (status /= ESMF_FIELDSTATUS_COMPLETE) then
       call ESMF_LogWrite(trim(subname)//": ERROR data not allocated ", ESMF_LOGMSG_INFO, rc=rc)
       rc = ESMF_FAILURE
       return
    else
       call ESMF_FieldGet(lfield, mesh=lmesh, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out

       call ESMF_MeshGet(lmesh, numOwnedNodes=nnodes, numOwnedElements=nelements, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out

       if (nnodes == 0 .and. nelements == 0) then
          call ESMF_LogWrite(trim(subname)//": no local nodes or elements ", ESMF_LOGMSG_INFO)
          rc = ESMF_FAILURE
          return
       end if

       if (present(fldptr1d)) then
          call ESMF_FieldGet(lfield, farrayPtr=fldptr1d, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
          if ( debug > 0)  then
             write(iulog,F01)' in '//trim(subname)//'fldptr1d for '//trim(fldname)//' is  '
          end if
          !print *, "FLDPTR1D is"
          !print *, FLDPTR1d
       else if (present(fldptr2d)) then
          call ESMF_FieldGet(lfield, farrayPtr=fldptr2d, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
       else
          !call shr_sys_abort("either fldptr1d or fldptr2d must be an input argument")
       end if
     endif  ! status


  end subroutine state_getfldptr




end module cpl_mod

