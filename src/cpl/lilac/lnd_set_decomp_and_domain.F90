module lnd_set_decomp_and_domain

  use ESMF
  use shr_kind_mod      , only : r8 => shr_kind_r8, cl=>shr_kind_cl
  use spmdMod           , only : masterproc
  use clm_varctl        , only : iulog
  use perf_mod          , only : t_startf, t_stopf, t_barrierf

  implicit none
  private ! except

  ! Module public routines
  public :: lnd_set_decomp_and_domain_from_meshinfo

  ! Module private routines
  private :: chkerr

  character(len=*) , parameter :: u_FILE_u = &
       __FILE__
  character(len=*), parameter, private :: sourcefile = &
       __FILE__

!===============================================================================
contains
!===============================================================================

  subroutine lnd_set_decomp_and_domain_from_meshinfo(model_meshfile, mesh_ctsm, ni, nj, rc)

    use decompInitMod , only : decompInit_ocn, decompInit_lnd, decompInit_lnd3D
    use domainMod     , only : ldomain, domain_init, lon1d, lat1d
    use decompMod     , only : ldecomp, bounds_type, get_proc_bounds
    use clm_varpar    , only : nlevsoi
    use clm_varctl    , only : fatmlndfrc, fsurdat, use_soil_moisture_streams, single_column
    use clm_varcon    , only : re
    use ncdio_pio     , only : ncd_io, file_desc_t, ncd_pio_openfile, ncd_pio_closefile, ncd_inqdlen
    use abortutils    , only : endrun
    use shr_log_mod   , only : errMsg => shr_log_errMsg
    use fileutils     , only : getfil

    ! input/output variables
    character(len=*)    , intent(in)    :: model_meshfile
    type(ESMF_Mesh)     , intent(out)   :: mesh_ctsm
    integer             , intent(out)   :: ni,nj  ! global sizes of dimensions
    integer             , intent(out)   :: rc

    ! local variables
    type(ESMF_VM)          :: vm
    type(ESMF_Mesh)        :: mesh_input
    type(ESMF_DistGrid)    :: distgrid_ctsm
    type(ESMF_DistGrid)    :: distgrid_input
    character(CL)          :: cvalue         ! config data
    integer                :: nlnd, nocn     ! local size ofarrays
    integer                :: g,n            ! indices
    type(bounds_type)      :: bounds         ! bounds
    integer                :: begg,endg
    integer  , pointer     :: gindex_lnd(:)  ! global index space for just land points
    integer  , pointer     :: gindex_ocn(:)  ! global index space for just ocean points
    integer  , pointer     :: gindex(:)      ! global index space for land and ocean points
    integer  , pointer     :: gindex_temp(:) ! temporary global index space
    integer  , pointer     :: mask(:)        ! local land/ocean mask
    integer  , pointer     :: lndmask_loc(:)
    real(r8) , pointer     :: lndfrac_loc(:)
    integer  , pointer     :: lndmask_glob(:)
    real(r8) , pointer     :: lndfrac_glob(:)
    real(r8) , pointer     :: lndlats_glob(:)
    real(r8) , pointer     :: lndlons_glob(:)
    real(r8) , pointer     :: rtemp_glob(:)
    integer  , pointer     :: itemp_glob(:)
    real(r8) , pointer     :: dataptr1d(:)
    integer                :: lsize,gsize
    logical                :: isgrid2d
    integer                :: numownedelements
    real(R8) , pointer     :: ownedElemCoords(:)
    integer                :: spatialDim
    type(ESMF_Field)       :: areaField
    type(ESMF_Array)       :: elemMaskArray
    character(len=CL)      :: locfn 
    type(file_desc_t)      :: ncid  ! netcdf file id
    integer                :: dimid ! netCDF dimension id
    logical                :: readvar  ! read variable in or not
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    ! Get current vm
    call ESMF_VMGetCurrent(vm, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Determine global 2d sizes from read of dimensions of surface dataset
    if (masterproc) then
       write(iulog,*) 'Attempting to global dimensions from surface dataset'
       if (fsurdat == ' ') then
          write(iulog,*)'fsurdat must be specified'
          call endrun(msg=errMsg(sourcefile, __LINE__))
       endif
    endif
    call getfil(fsurdat, locfn, 0 )
    call ncd_pio_openfile (ncid, trim(locfn), 0)
    call ncd_inqdlen(ncid, dimid, ni, 'lsmlon')
    call ncd_inqdlen(ncid, dimid, nj, 'lsmlat')
    call ncd_pio_closefile(ncid)
    gsize = ni*nj  
    if (single_column) then
       isgrid2d = .true.
    else if (nj == 1) then
       isgrid2d = .false.
    else
       isgrid2d = .true.
    end if
    if (masterproc) then
       write(iulog,'(a,2(i8,2x))') 'global ni,nj = ',ni,nj
       if (isgrid2d) then
          write(iulog,'(a)') 'model grid is 2-dimensional'
       else
          write(iulog,'(a)') 'model grid is not 2-dimensional'
       end if
    end if

    ! read in the land mesh from the file
    mesh_input = ESMF_MeshCreate(filename=trim(model_meshfile), fileformat=ESMF_FILEFORMAT_ESMFMESH, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (masterproc) then
       write(iulog,'(a)')'land mesh file ',trim(model_meshfile)
    end if

    ! Obtain global land amsk
    if (trim(fatmlndfrc) /= 'null') then
       if (masterproc) then
          write(iulog,*) 'Generating ctsm decomposition from ',trim(fatmlndfrc)
       endif
    else
       if (masterproc) then
          write(iulog,*) 'Generating ctsm decomposition from ',trim(model_meshfile)
       endif
    end if

    allocate(lndmask_glob(ni*nj)); lndmask_glob(:) = 0
    allocate(rtemp_glob(gsize))
    
    if (trim(fatmlndfrc) /= 'null') then

       ! Read in global land mask and land fraction from fatmlndfrc
       call getfil( trim(fatmlndfrc), locfn, 0 )
       call ncd_pio_openfile (ncid, trim(locfn), 0)
       call ncd_io(ncid=ncid, varname='mask', data=lndmask_glob, flag='read', readvar=readvar)
       if (.not. readvar) call endrun( msg=' ERROR: variable mask not on fatmlndfrc file'//errMsg(sourcefile, __LINE__))
       allocate(lndfrac_glob(ni*nj)); lndfrac_glob(:) = 0._r8
       call ncd_io(ncid=ncid, varname='frac', data=lndfrac_glob, flag='read', readvar=readvar)
       if (.not. readvar) call endrun( msg=' ERROR: variable frac not on fatmlndfrc file'//errMsg(sourcefile, __LINE__))
       call ncd_pio_closefile(ncid)

    else

       ! Obtain land mask from land mesh file - ASSUME THAT LAND FRAC IS IDENTICAL TO LAND MASK
       call ESMF_MeshGet(mesh_input, elementdistGrid=distgrid_input, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_DistGridGet(distgrid_input, localDe=0, elementCount=lsize, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       ! Determine lndmask_loc
       allocate(lndmask_loc(lsize))
       elemMaskArray = ESMF_ArrayCreate(distgrid_input, lndmask_loc, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       ! The following calls fills in the values of lndmask_loc
       call ESMF_MeshGet(mesh_input, elemMaskArray=elemMaskArray, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       ! Determine lndfrac_loc
       ! ASSUME that land fraction is identical to land mask in this case
       allocate(lndfrac_loc(lsize))
       lndfrac_loc(:) = lndmask_loc(:)
       
       ! determine global landmask_glob - needed to determine the ctsm decomposition
       ! land frac, lats, lons and areas will be done below
       allocate(gindex_temp(lsize))
       call ESMF_DistGridGet(distgrid_input, 0, seqIndexList=gindex_temp, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       allocate(lndmask_glob(gsize)); lndmask_glob(:) = 0
       do n = 1,lsize
          lndmask_glob(gindex(n)) = lndmask_loc(n)
       end do
       allocate(itemp_glob(gsize))
       call ESMF_VMAllReduce(vm, sendData=lndmask_glob, recvData=itemp_glob, count=gsize, &
            reduceflag=ESMF_REDUCE_SUM, rc=rc)
       lndmask_glob(:) = int(itemp_glob(:))
       deallocate(itemp_glob)

    end if

    ! Determine lnd decomposition that will be used by ctsm
    call decompInit_lnd(lni=ni, lnj=nj, amask=lndmask_glob)
    if (use_soil_moisture_streams) then
       call decompInit_lnd3D(lni=ni, lnj=nj, lnk=nlevsoi)
    end if

    ! *** Get JUST gridcell processor bounds ***
    ! Remaining bounds (landunits, columns, patches) will be set after calling decompInit_glcp
    ! so get_proc_bounds is called twice and the gridcell information is just filled in twice
    call get_proc_bounds(bounds)
    begg = bounds%begg
    endg = bounds%endg

    ! Create gindex_lnd 
    nlnd = endg - begg + 1
    allocate(gindex_lnd(nlnd))
    do g = begg, endg
       n = 1 + (g - begg)
       gindex_lnd(n) = ldecomp%gdc2glo(g)
    end do

    ! Create gindex_ocn
    ! Need this decomposition to create the full mesh
    ! Note that the memory for gindex_ocn will be allocated in the following call
    call decompInit_ocn(ni=ni, nj=nj, amask=lndmask_glob, gindex_ocn=gindex_ocn)

    ! Create a global index that includes both land and ocean points
    nocn = size(gindex_ocn)
    allocate(gindex(nlnd + nocn))
    do n = 1,nlnd+nocn
       if (n <= nlnd) then
          gindex(n) = gindex_lnd(n)
       else
          gindex(n) = gindex_ocn(n-nlnd)
       end if
    end do

    ! Generate a new distgrid based on gindex
    distgrid_ctsm = ESMF_DistGridCreate(arbSeqIndexList=gindex, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Generate the ctsm mesh on the gindex decomposition
    mesh_ctsm = ESMF_MeshCreate(mesh_input, elementDistGrid=distgrid_ctsm, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Initialize domain data structure
    call domain_init(domain=ldomain, isgrid2d=isgrid2d, ni=ni, nj=nj, nbeg=begg, nend=endg)

    ! Determine ldomain%mask
    do g = begg, endg
       n = gindex(g-begg+1)
       ldomain%mask(g) = lndmask_glob(n)
    end do
    deallocate(lndmask_glob)

    ! Determine ldomain%frac 
    ! note that lndfrac_glob was read in from fatmlndfrc above if it was not set to null
    if (trim(fatmlndfrc) == 'null') then
       allocate(lndfrac_glob(gsize))
       do n = 1,nlnd
          lndfrac_glob(gindex_lnd(n)) = lndfrac_loc(n)
       end do
       call ESMF_VMAllReduce(vm, sendData=lndfrac_glob, recvData=rtemp_glob, count=gsize, reduceflag=ESMF_REDUCE_SUM, rc=rc)
       do g = begg, endg
          n = gindex(g-begg+1)
          ldomain%frac(g) = rtemp_glob(n)
       end do
       deallocate(lndfrac_glob)
    else
       do g = begg, endg
          n = gindex(g-begg+1)
          ldomain%frac(g) = lndfrac_glob(n)
       end do
       deallocate(lndfrac_glob)
    end if

    ! Determine ldoman%latc and ldomain%lonc
    call ESMF_MeshGet(mesh_ctsm, spatialDim=spatialDim, numOwnedElements=numOwnedElements, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    allocate(ownedElemCoords(spatialDim*numownedelements))
    call ESMF_MeshGet(mesh_ctsm, ownedElemCoords=ownedElemCoords)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_MeshGet(mesh_ctsm, ownedElemCoords=ownedElemCoords, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    do g = begg,endg
       n = g - begg + 1
       ldomain%lonc(g) = ownedElemCoords(2*n-1)
       if (ldomain%lonc(g) == 360._r8) ldomain%lonc(g) = 0._r8 ! TODO: why the difference?
       ldomain%latc(g) = ownedElemCoords(2*n)
    end do

    ! Determine ldomain%area by querying the mesh on the ctsm decomposition
    areaField = ESMF_FieldCreate(mesh_ctsm, ESMF_TYPEKIND_R8, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldRegridGetArea(areaField, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(areaField, farrayPtr=dataptr1d, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    do g = begg, endg
       ldomain%area(g) = dataptr1d(g-begg+1) * (re*re)
    end do
    call ESMF_FieldDestroy(areaField)

    ! If grid is 2d, determine lon1d and lat1d
    if (isgrid2d) then
       ! Determine lon1d
       allocate(lndlons_glob(gsize))
       lndlons_glob(:) = 0._r8
       do n = 1,numownedelements
          if (ownedElemCoords(2*n-1) == 360._r8) then ! TODO: why is this needed?
             lndlons_glob(gindex(n)) = 0._r8
          else
             lndlons_glob(gindex(n)) = ownedElemCoords(2*n-1)
          end if
       end do
       call ESMF_VMAllReduce(vm, sendData=lndlons_glob, recvData=rtemp_glob, count=gsize, reduceflag=ESMF_REDUCE_SUM, rc=rc)
       deallocate(lndlons_glob)
       allocate(lon1d(ni))
       do n = 1,ni
          lon1d(n) = rtemp_glob(n)
       end do

       ! Determine lat1d
       allocate(lndlats_glob(gsize))
       lndlats_glob(:) = 0._r8
       do n = 1,numownedelements
          lndlats_glob(gindex(n)) = ownedElemCoords(2*n)
       end do
       call ESMF_VMAllReduce(vm, sendData=lndlats_glob, recvData=rtemp_glob, count=gsize, reduceflag=ESMF_REDUCE_SUM, rc=rc)
       deallocate(lndlats_glob)
       allocate(lat1d(nj))
       do n = 1,nj
          lat1d(n) = rtemp_glob((n-1)*ni + 1)
       end do
    end if

    deallocate(ownedElemCoords)
    deallocate(rtemp_glob)
    deallocate(gindex)

  end subroutine lnd_set_decomp_and_domain_from_meshinfo

  !===============================================================================
  logical function chkerr(rc, line, file)
    integer, intent(in) :: rc
    integer, intent(in) :: line
    character(len=*), intent(in) :: file
    integer :: lrc
    chkerr = .false.
    lrc = rc
    if (ESMF_LogFoundError(rcToCheck=lrc, msg=ESMF_LOGERR_PASSTHRU, line=line, file=file)) then
       chkerr = .true.
    endif
  end function chkerr

end module lnd_set_decomp_and_domain
