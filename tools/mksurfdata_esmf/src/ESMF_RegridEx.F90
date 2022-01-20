program ESMF_RegridEx

  !==============================================================================
  ! This program shows examples of using Regrid on Field data
  !-----------------------------------------------------------------------------
  use ESMF

  implicit none

  ! Local variables to be used in the Regrid method calls.
  ! The code creating and filling these variables is not included in the
  ! example documentation because those interfaces are not specific to
  ! Regrid.
  type(ESMF_Field)                             :: field1, field2
  type(ESMF_IGrid)                             :: srcigrid, dstigrid
  type(ESMF_RouteHandle)                       :: regrid_rh
  type(ESMF_DELayout)                          :: layout1, layout2
  integer                                      :: rc
  integer                                      :: finalrc, npets
  integer                                      :: i, j, lb(2), ub(2), halo
  type(ESMF_ArraySpec)                         :: arrayspec
  type(ESMF_VM)                                :: vm
  real (ESMF_KIND_R8), dimension(:,:), pointer :: f90ptr1
  real (ESMF_KIND_R8), dimension(2)            :: mincoords, maxcoords

  finalrc = ESMF_SUCCESS

  !-------------------------------------------------------------------------
  !   ! Setup:
  !   !
  !   !  Create a source and destination igrid with data on it, to use
  !   !  in the Regrid calls below.

  call ESMF_Initialize(rc=rc)
  call ESMF_VMGetGlobal(vm, rc=rc)
  call ESMF_VMGet(vm, petCount=npets, rc=rc)

  layout1 = ESMF_DELayoutCreate(vm, (/ 1, npets /), rc=rc)
  layout2 = ESMF_DELayoutCreate(vm, (/ npets, 1 /), rc=rc)

  mincoords = (/  0.0,  0.0 /)
  maxcoords = (/ 20.0, 30.0 /)
  srcigrid = ESMF_IGridCreateHorzXYUni((/ 90, 180 /), mincoords, maxcoords, &
       horzStagger=ESMF_IGRID_HORZ_STAGGER_A, name="srcigrid", rc=rc)
  call ESMF_IGridDistribute(srcigrid, delayout=layout1, rc=rc)

  ! same igrid coordinates, but different layout
  dstigrid = ESMF_IGridCreateHorzXYUni((/ 90, 180 /), mincoords, maxcoords, horzStagger=ESMF_IGRID_HORZ_STAGGER_A, name="srcigrid", rc=rc)
  call ESMF_IGridDistribute(dstigrid, delayout=layout2, rc=rc)
  call ESMF_ArraySpecSet(arrayspec, 2, ESMF_TYPEKIND_R8, rc)

  ! allow for a halo width of 3, let field create data space
  halo = 3
  field1 = ESMF_FieldCreate(srcigrid, arrayspec, horzRelloc=ESMF_CELL_CENTER, haloWidth=3, name="src pressure", rc=rc)

  ! get a fortran pointer to the data spacd
  call ESMF_FieldGetDataPointer(field1, f90ptr1, ESMF_DATACOPY_REFERENCE, rc=rc)

  lb(:) = lbound(f90ptr1)
  ub(:) = ubound(f90ptr1)

  f90ptr1(:,:) = 0.0
  do j=lb(2)+halo, ub(2)-halo
     do i=lb(1)+halo, ub(1)-halo
        f90ptr1(i, j) = i*1000 + j
     enddo
  enddo

  field2 = ESMF_FieldCreate(dstigrid, arrayspec, horzRelloc=ESMF_CELL_CENTER, name="dst pressure", rc=rc)

  ! fields all ready to go
  !  The user has already created an {\tt ESMF\_IGrid}, an
  !  {\tt ESMF\_Array} with data, and put them together in an {\tt ESMF\_Field}.
  !  An {\tt ESMF\_RouteHandle} is created by the regrid store call 
  !  and the data movement needed to
  !  execute the regrid is stored with that handle by the store method. 
  !  To actually execute the operation, the source and destination data
  !  objects must be supplied, along with the same {\tt ESMF\_RouteHandle}.

  call ESMF_FieldRegridStore(field1, field2, vm, routehandle=regrid_rh, regridmethod=ESMF_REGRIDMETHOD_BILINEAR, rc=rc)
  call ESMF_FieldRegrid(field1, field2, regrid_rh, rc=rc)
  call ESMF_FieldRegridRelease(regrid_rh, rc=rc)

  call ESMF_FieldDestroy(field1, rc=rc)
  call ESMF_FieldDestroy(field2, rc=rc)
  call ESMF_IGridDestroy(srcigrid, rc=rc)
  call ESMF_IGridDestroy(dstigrid, rc=rc)
  call ESMF_Finalize(rc=rc)

end program ESMF_RegridEx

    
    
