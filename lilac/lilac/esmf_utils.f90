module esmf_utils

  ! Wrappers and derived types exposing ESMF components to LILAC


#include "ESMC.h"
  use ESMF

  implicit none
  private

  ! Consider renaming ESMFInfoType (add lilac to name)
  type, public :: ESMFInfoType
     private
     character(len=MAXFILELENGTH) :: name

     type(ESMF_VM)                :: vm
     type(ESMF_State)             :: land_import
     type(ESMF_State)             :: land_export
     type(ESMF_State)             :: atmos_import
     type(E SMF_State)             :: atmos_export
     type(ESMF_GridComp)          :: atmos_comp
     type(ESMF_GridComp)          :: land_comp
     type(ESMF_CplComp)           :: cpl_comp

   contains
     procedure, public   :: init      => init
     procedure, public   :: run       => run
     procedure, public   :: final     => final

     procedure, private :: atmos_register => atmos_register
     procedure, private :: land_register  => land_register
     procedure, private :: cpl_register => cpl_register
  end type ESMFInfoType

contains

  subroutine init(self, name)
    implicit none
    class(ESMFInfoType), intent(inout)  :: self
    character(len=MAXVARLENGTH), intent(in) :: name

    ! TODO define subroutines: https://stackoverflow.com/questions/32809769/how-to-pass-subroutine-names-as-arguments-in-fortran

    ! Local variables
    integer :: localPet, petCount, localrc, rc=ESMF_SUCCESS, userrc=ESMF_SUCCESS
    character(len=ESMF_MAXSTR) :: cname1, cname2

    print *, "esmf_info%init()"

    self%name = name

    ! Create section
    !-------------------------------------------------------------------------

    ! Initialize framework and get back default global VM

    ! only run if not esmf_isintialized()
    call ESMF_Initialize(vm=self%vm, defaultlogfilename="lilac.log", logkindflag=ESMF_LOGKIND_MULTI, rc=localrc)
    call check(localrc, rc)

    ! Get number of PETs we are running with
    call ESMF_VMGet(self%vm, petCount=petCount, localPet=localPet, rc=localrc)
    call check(localrc, rc)

    ! Create the 2 model components and a coupler
    cname1 = "land"
    ! use petList to define land on all PET
    self%land_grid = ESMF_GridCompCreate(name=cname1, rc=localrc)
    print *, "Created component ", trim(cname1), "rc =", localrc
    call check(localrc, rc)

    cname2 = "atmosphere"
    ! use petList to define atmosphere on all PET
    self%atmos_comp = ESMF_GridCompCreate(name=cname2, rc=localrc)
    print *, "Created component ", trim(cname2), "rc =", localrc
    call check(localrc, rc)

    cplname = "lilac coupler"
    ! no petList means that coupler component runs on all PETs
    self%cpl_comp = ESMF_CplCompCreate(name=cplname, rc=localrc)
    print *, "Created component ", trim(cplname), ", rc =", localrc
    call check(localrc, rc)

    print *, "Comp Creates finished"

    ! Register section
    !-------------------------------------------------------------------------
    call ESMF_GridCompSetServices(self%atmos_comp, userRoutine=atmos_register, userRc=userrc, rc=localrc)
    print *, "atmos SetServices finished, rc= ", localrc
    call check(localrc, rc)
    call check(userrc, rc)

    call ESMF_GridCompSetServices(self%land_comp, userRoutine=land_register, userRc=userrc, rc=localrc)
    print *, "land SetServices finished, rc= ", localrc
    call check(localrc, rc)
    call check(userrc, rc)

    call ESMF_CplCompSetServices(self%cpl_comp, userRoutine=cpl_register, userRc=userrc, rc=localrc)
    print *, "Cpl SetServices finished, rc= ", localrc
    call check(localrc, rc)
    call check(userrc, rc)

    ! Init section
    !-------------------------------------------------------------------------
    ! land import/export states
    self%land_import = ESMF_StateCreate(name="land import", stateintent=ESMF_STATEINTENT_IMPORT, rc=localrc)
    call check(localrc, rc)
    self%land_export = ESMF_StateCreate(name="land export", stateintent=ESMF_STATEINTENT_EXPORT, rc=localrc)
    call check(localrc, rc)
    call ESMF_GridCompInitialize(land, importState=self%land_import, exportState=self%land_export, userRc=userrc, rc=localrc)
    call check(localrc, rc)
    call check(userrc, rc)
    print *, "Land Initialize finished, rc =", localrc

    ! atmosphere import/export state
    self%atmos_import = ESMF_StateCreate(name="atmos import",  &
         stateintent=ESMF_STATEINTENT_IMPORT, rc=localrc)
    call check(localrc, rc)

    self%atmos_export = ESMF_StateCreate(name="atmos export",  &
         stateintent=ESMF_STATEINTENT_EXPORT, rc=localrc)
    call check(localrc, rc)
    call ESMF_GridCompInitialize(self%atmos_comp, exportState=self%atmos_export, userRc=userrc, rc=localrc)
    print *, "Atmosphere Initialize finished, rc =", localrc
    call check(localrc, rc)
    call check(userrc, rc)

    ! call ESMF_CPLCompInitialize twice (once for each grid comp)

  end subroutine init

  subroutine run(self)
    implicit none
    integer :: localrc, rc=ESMF_SUCCESS, userrc=ESMF_SUCCESS
    print *, "esmf_info%run()"

    ! TODO: need some help on order of imports/exports/runs and whether the land/atm both need import/export states

    ! atmosphere run
    ! copy the atmos state and put it into atmos export
    call ESMF_GridCompRun(self%atmos_comp, exportState=self%atmos_export, phase=1, userRc=userrc, rc=localrc)
    print *, "Atmosphere Run returned, rc =", localrc
    call check(localrc, rc)
    call check(userrc, rc)

    ! coupler run
    call ESMF_CplCompRun(self%cpl_comp, importState=self%atoms_export, exportState=self%land_import, &
         userRc=userrc, rc=localrc)
    print *, "Coupler Run returned, rc =", localrc
    call check(localrc, rc)
    call check(userrc, rc)

    ! land run
    call ESMF_GridCompRun(self%land_comp, importState=self%land_import, exportState=self%land_export, userRc=userrc, rc=localrc)
    print *, "Land Run returned, rc =", localrc
    call check(localrc, rc)
    call check(userrc, rc)

    ! coupler run
    call ESMF_CplCompRun(self%cpl_comp, importState=self%land_export, exportState=self%atmos_import, &
         userRc=userrc, rc=localrc)
    print *, "Coupler Run returned, rc =", localrc
    call check(localrc, rc)
    call check(userrc, rc)

    call ESMF_GridCompRun(self%atmos_comp, importState%atmos_import, phase=2, userRc=userrc, rc=localrc)
    print *, "Atmosphere Run returned, rc =", localrc
    call check(localrc, rc)
    call check(userrc, rc)

  end subroutine run

  subroutine final(self)
    implicit none
    class(ESMFInfoType), intent(inout)  :: self
    integer :: localrc, rc=ESMF_SUCCESS

    print *, "esmf_info%final()"

    ! Destroy section
    call ESMF_GridCompDestroy(self%atmos_comp, rc=localrc)
    check(localrc, rc)
    call ESMF_GridCompDestroy(self%land_comp, rc=localrc)
    check(localrc, rc)
    call ESMF_CplCompDestroy(self%cpl_comp, rc=localrc)
    check(localrc, rc)

    call ESMF_StateDestroy(self%land_export, rc=localrc)
    call ESMF_StateDestroy(self%land_import, rc=localrc)
    check(localrc, rc)
    call ESMF_StateDestroy(self%atmos_export, rc=localrc)
    call ESMF_StateDestroy(self%atmos_import, rc=localrc)
    ! do this everywhere
    if return_error(localrc, rc) return

    print *, "All Destroy routines done"

  end subroutine final

  subroutine atoms_register(comp, rc)
    type(ESMF_GridComp)   :: comp   ! must not be optional
    integer, intent(out)  :: rc     ! must not be optional

    ! Set the entry points for standard ESMF Component methods
    call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_INITIALIZE, &
         userRoutine=atoms_init, rc=rc)
    call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_RUN, &
         userRoutine=atoms_copy_atm_to_lilac, phase=1, rc=rc)
    call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_RUN, &
         userRoutine=atoms_copy_lilac_to_atm, phase=2, rc=rc)
    call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_FINALIZE, &
         userRoutine=atoms_final, rc=rc)
    ! TODO: check rcs

    rc = ESMF_SUCCESS

  end subroutine atoms_register

  subroutine land_register(comp, rc)
    type(ESMF_GridComp)   :: comp   ! must not be optional
    integer, intent(out)  :: rc     ! must not be optional

    ! land_* comes from ctsm esmf cap

    ! Set the entry points for standard ESMF Component methods
    call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_INITIALIZE, &
         userRoutine=land_init, rc=rc)
    call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_RUN, &
         userRoutine=land_run, rc=rc)
    call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_FINALIZE, &
         userRoutine=land_final, rc=rc)
    ! TODO: check rcs

    rc = ESMF_SUCCESS

  end subroutine land_register

  subroutine cpl_register(comp, rc)
    type(ESMF_GridComp)   :: comp   ! must not be optional
    integer, intent(out)  :: rc     ! must not be optional

    rc = ESMF_FAILURE

    ! Register the callback routines.

    call ESMF_CplCompSetEntryPoint(comp, ESMF_METHOD_INITIALIZE, userRoutine=coupler_init, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT, rc=rc)
    call ESMF_CplCompSetEntryPoint(comp, ESMF_METHOD_RUN, userRoutine=coupler_run, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT, rc=rc)
    call ESMF_CplCompSetEntryPoint(comp, ESMF_METHOD_FINALIZE, userRoutine=coupler_final, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT, rc=rc)

    print *, "CouplerMod: Registered Initialize, Run, and Finalize routines"

    rc = ESMF_SUCCESS

  end subroutine cpl_register

  function return_error(rc, returnrc)
    ! fight with this later
    integer, intent(in) :: rc, returnrc
    if (ESMF_LogFoundError(rc, ESMF_ERR_PASSTHRU, ESMF_CONTEXT, rcToReturn=returnrc)) then
       return_error = .true.
    else
       return_error = .false.
    endif

  end function return_error

end module esmf_utils
