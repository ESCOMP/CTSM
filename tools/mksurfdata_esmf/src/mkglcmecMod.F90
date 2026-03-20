module mkglcmecMod

  !-----------------------------------------------------------------------
  ! Make glacier multi-elevation class  data
  !-----------------------------------------------------------------------

  use ESMF
  use shr_kind_mod   , only : r8 => shr_kind_r8, r4=>shr_kind_r4
  use shr_sys_mod    , only : shr_sys_abort
  use pio            , only : file_desc_t, pio_openfile, pio_closefile, pio_nowrite
  use pio            , only : var_desc_t, io_desc_t, Pio_Offset_Kind, pio_setframe
  use pio            , only : pio_inq_dimid, pio_inq_dimlen, pio_inq_varid
  use pio            , only : pio_put_var, pio_get_var
  use mkpioMod       , only : mkpio_get_rawdata, pio_iotype, pio_iosystem
  use mkpioMod       , only : mkpio_iodesc_rawdata, mkpio_get_rawdata_level
  use mkesmfMod      , only : regrid_rawdata, create_routehandle_r8
  use mkvarctl       , only : ndiag, root_task, outnc_3dglc, spval, nglcec
  use mkutilsMod     , only : chkerr
  use mkutilsMod     , only : slightly_below, slightly_above
  use mkfileMod      , only : mkfile_output
  use mkdiagnosticsMod , only : output_diagnostics_area

  implicit none
  private           ! By default make data private

  public  :: mkglcmecInit      ! Initialization
  public  :: mkglcmec          ! Set glacier multi-elevation class
  public  :: mkglacier         ! Set percent glacier

  private :: get_elevclass     ! get elevation class index
  private :: mean_elevation_vc ! get the elevation of a virtual column

  real(r8), allocatable :: elevclass(:)          ! elevation classes

  character(len=*) , parameter :: u_FILE_u = &
       __FILE__

!=================================================================================
contains
!=================================================================================

  subroutine mkglcmecInit( pioid_o )
    !
    ! Initialize of Make glacier multi-elevation class data
    !
    ! input/output variables
    type(file_desc_t) , intent(inout) :: pioid_o

    ! local variables:
    type(var_desc_t)      :: pio_varid
    type(io_desc_t)       :: pio_iodesc
    integer               :: rcode
    character(len=*), parameter :: subname = 'mkglcmecInit'
    !-----------------------------------------------------------------------

    allocate( elevclass(nglcec+1) )

    ! Define elevation classes, represents lower boundary of each class

    if ( nglcec == 36 )then
       elevclass(:) = (/ 0.,     200.,   400.,   600.,   800.,  &
                        1000.,  1200.,  1400.,  1600.,  1800.,  &
                        2000.,  2200.,  2400.,  2600.,  2800.,  &
                        3000.,  3200.,  3400.,  3600.,  3800.,  &
                        4000.,  4200.,  4400.,  4600.,  4800.,  &
                        5000.,  5200.,  5400.,  5600.,  5800.,  &
                        6000.,  6200.,  6400.,  6600.,  6800.,  &
                        7000., 10000./)
    else if ( nglcec == 10 )then
       elevclass(1)  =     0.
       elevclass(2)  =   200.
       elevclass(3)  =   400.
       elevclass(4)  =   700.
       elevclass(5)  =  1000.
       elevclass(6)  =  1300.
       elevclass(7)  =  1600.
       elevclass(8)  =  2000.
       elevclass(9)  =  2500.
       elevclass(10) =  3000.
       elevclass(11) = 10000.
    else if ( nglcec == 5  )then
       elevclass(1)  =     0.
       elevclass(2)  =   500.
       elevclass(3)  =  1000.
       elevclass(4)  =  1500.
       elevclass(5)  =  2000.
       elevclass(6)  = 10000.
    else if ( nglcec == 3  )then
       elevclass(1)  =     0.
       elevclass(2)  =  1000.
       elevclass(3)  =  2000.
       elevclass(4)  = 10000.
    else if ( nglcec == 1  )then
       elevclass(1)  =     0.
       elevclass(2)  = 10000.
    else
       write(6,*) subname//"error:: nglcec must be 1, 3, 5, 10 or 36",&
            " to work with CLM: "
       call shr_sys_abort()
    end if

    if (root_task)  write(ndiag, '(a)') trim(subname)//" writing out GLC_MEC"
    rcode = pio_inq_varid(pioid_o, 'GLC_MEC', pio_varid)
    rcode = pio_put_var(pioid_o, pio_varid, elevclass)

  end subroutine mkglcmecInit

  !=================================================================================
  subroutine mkglcmec(file_mesh_i, file_data_i, mesh_o,  pioid_o, rc)
    !
    ! make percent glacier on multiple elevation classes, mean elevation for each
    ! elevation class, and associated fields
    !
    ! Note that the raw glacier data are specified by level, and thus implicitly include the
    ! necessary topo data for breaking pct glacier into elevation classes. Each level in the
    ! input data is assigned to an elevation (given by BIN_CENTERS in the input data). Thus,
    ! all of the input glacier in level 1 is treated as being at the same elevation, and
    ! likewise for each other level. These elevations are then used in assigning pct_glacier
    ! to the appropriate elevation class in the output data, as well as determining the mean
    ! topographic height of each elevation class in the output data.
    !
    ! Note that the various percentages computed here are given as % of the glc_mec landunit.
    ! If the input glacier area is 0 for a given grid cell, this requires setting these %
    ! variables in an arbitrary way.
    !
    ! input/output variables
    character(len=*)  , intent(in)    :: file_mesh_i ! input mesh file name
    character(len=*)  , intent(in)    :: file_data_i ! input data file name
    type(ESMF_Mesh)   , intent(in)    :: mesh_o      ! output mesh
    type(file_desc_t) , intent(inout) :: pioid_o
    integer           , intent(out)   :: rc

    ! local variables:
    type(ESMF_RouteHandle) :: routehandle_nonorm
    type(ESMF_Mesh)        :: mesh_i
    type(file_desc_t)      :: pioid_i
    type(var_desc_t)       :: pio_varid
    type(io_desc_t)        :: pio_iodesc
    integer                :: pio_vartype
    integer                :: ni,no,lev
    integer                :: ns_i, ns_o
    integer                :: n,m,k                     ! indices
    integer                :: dimid                     ! dimension ids
    integer                :: ndims                     ! number of dimensions in input variables
    integer                :: nlev                      ! number of levels in input file
    integer , allocatable  :: mask_i(:)
    real(r8), allocatable  :: frac_i(:)
    real(r8), allocatable  :: frac_o_nonorm(:)
    real(r8), allocatable  :: area_i(:)
    real(r8), allocatable  :: area_o(:)
    real(r8), allocatable  :: data_pctglc_i(:)
    real(r8), allocatable  :: data_pctglc_o(:)
    real(r8), allocatable  :: data_pctglc_icesheet_i(:) ! input icesheet percentage
    real(r8), allocatable  :: data_pctglc_icesheet_o(:) ! input icesheet percentage
    real(r8), allocatable  :: data_pctglc_gic_i(:)
    real(r8), allocatable  :: data_pctglc_gic_o(:)
    real(r8), allocatable  :: topoglcmec_unnorm_o(:,:)  ! same as topoglcmec_o, but unnormalized
    real(r8), allocatable  :: pctglc_tot_o(:)           ! total glacier cover for the grid cell
    real(r4), allocatable  :: topoice_i(:)              ! topographic height of this level
    real(r8), allocatable  :: pctglcmec_o (:,:)         ! % for each elevation class on output glacier grid (% of landunit)
    real(r8), allocatable  :: topoglcmec_o(:,:)         ! mean elevation for each elevation classs on output glacier grid
    real(r8), allocatable  :: pctglcmec_gic_o(:,:)      ! % glc gic on output grid, by elevation class (% of landunit)
    real(r8), allocatable  :: pctglcmec_icesheet_o(:,:) ! % glc ice sheet on output grid, by elevation class (% of landunit)
    real(r8), allocatable  :: pctglc_gic_o(:)           ! % glc gic on output grid, summed across elevation classes (% of landunit)
    real(r8), allocatable  :: pctglc_icesheet_o(:)      ! % glc ice sheet on output grid, summed across elevation classes (% of landunit)
    integer                :: unlimited_index           ! z level
    real(r8)               :: glc_sum                   ! temporary
    integer                :: ier, rcode                ! error status
    logical                :: errors                    ! error status
    real(r8), parameter :: eps = 2.e-5_r8              ! epsilon for error checks (note that we use a large-ish value
                                                       ! because data are stored as single-precision floats in the raw dataset)
    real(r8), parameter :: eps_small = 1.e-12_r8       ! epsilon for error checks that expect close match
    character(len=*), parameter :: subname = 'mkglcmec'
    !-----------------------------------------------------------------------

    rc = ESMF_SUCCESS

    if (root_task) then
       write(ndiag,*)
       write(ndiag,'(1x,80a1)') ('=',k=1,80)
       write(ndiag,*)
       write(ndiag,'(a)') 'Attempting to make percent elevation class ',&
            'and mean elevation for glaciers .....'
       write(ndiag,'(a)') ' Input file is '//trim(file_data_i)
       write(ndiag,'(a)') ' Input mesh file is '//trim(file_mesh_i)
    end if

    ! Open input data file
    call ESMF_VMLogMemInfo("Before pio_openfile for "//trim(file_data_i))
    rcode = pio_openfile(pio_iosystem, pioid_i, pio_iotype, trim(file_data_i), pio_nowrite)
    call ESMF_VMLogMemInfo("After pio_openfile "//trim(file_data_i))

    ! Read in input mesh
    call ESMF_VMLogMemInfo("Before create mesh_i in "//trim(subname))
    mesh_i = ESMF_MeshCreate(filename=trim(file_mesh_i), fileformat=ESMF_FILEFORMAT_ESMFMESH, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMLogMemInfo("After create mesh_i in "//trim(subname))

    ! Determine ns_i
    call ESMF_MeshGet(mesh_i, numOwnedElements=ns_i, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Determine ns_o and allocate output data
    call ESMF_MeshGet(mesh_o, numOwnedElements=ns_o, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    allocate(pctglcmec_o(ns_o,nglcec))
    pctglcmec_o(:,:) = 0.

    allocate(topoglcmec_o(ns_o,nglcec))
    topoglcmec_o(:,:) = 0.
    if ( outnc_3dglc )then
       allocate(pctglcmec_gic_o(ns_o,nglcec))
       pctglcmec_gic_o(:,:) = 0.

       allocate(pctglcmec_icesheet_o(ns_o,nglcec))
       pctglcmec_icesheet_o(:,:) = 0.

       allocate(pctglc_gic_o(ns_o))
       pctglc_gic_o(:) = 0.

       allocate(pctglc_icesheet_o(ns_o))
       pctglc_icesheet_o(:) = 0.
    end if

    ! Get the landmask from the file and reset the mesh mask based on that
    allocate(frac_i(ns_i), stat=ier)
    if (ier/=0) call shr_sys_abort(subname//" error in allocating frac_i")
    allocate(mask_i(ns_i), stat=ier)
    if (ier/=0) call shr_sys_abort(subname//" error in allocating mask_i")
    call mkpio_get_rawdata(pioid_i, 'LANDMASK', mesh_i, frac_i, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    do ni = 1,ns_i
       if (frac_i(ni) > 0._r8) then
          mask_i(ni) = 1
       else
          mask_i(ni) = 0
       end if
    end do
    call ESMF_MeshSet(mesh_i, elementMask=mask_i, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! Read in BIN_CENTERS on all tasks and check validity
    rcode = pio_inq_dimid (pioid_i, 'z', dimid)
    rcode = pio_inq_dimlen(pioid_i, dimid, nlev)
    ! TODO: hard-wiring topoice to be r4 - this needs to be generalized to query the variable type
    ! on the netcdf file
    allocate(topoice_i(nlev))
    rcode = pio_inq_varid (pioid_i, 'BIN_CENTERS', pio_varid)
    rcode = pio_get_var(pioid_i, pio_varid, topoice_i)
    do lev = 1,nlev
       m = get_elevclass(topoice_i(lev))
       if (m < 1 .or. m > nglcec) then
          call shr_sys_abort(subname//" error m<1 or m > nglcec")
       end if
    end do

    ! Allocate memory for reading in one level at a time
    allocate(data_pctglc_i(ns_i), stat=ier)
    if (ier/=0) call shr_sys_abort(subname//" error in allocating data_pctglc_i")
    allocate(data_pctglc_o(ns_o), stat=ier)
    if (ier/=0) call shr_sys_abort(subname//" error in allocating data_pctglc_o")

    allocate(data_pctglc_gic_i(ns_i), stat=ier)
    if (ier/=0) call shr_sys_abort(subname//" error in allocating data_pctglc_gic_i")
    allocate(data_pctglc_gic_o(ns_o), stat=ier)
    if (ier/=0) call shr_sys_abort(subname//" error in allocating data_pctglc_gic_o")

    allocate(data_pctglc_icesheet_i(ns_i), stat=ier)
    if (ier/=0) call shr_sys_abort(subname//" error in allocating data_pctglc_icesheet_i")
    allocate(data_pctglc_icesheet_o(ns_o), stat=ier)
    if (ier/=0) call shr_sys_abort(subname//" error in allocating data_pctglc_icesheet_o")

    allocate(topoglcmec_unnorm_o(ns_o,nglcec), stat=ier)
    if (ier/=0) call shr_sys_abort(subname//" error in allocating topoglcmec_unnorm_o")
    topoglcmec_unnorm_o(:,:) = 0.

    ! Create iodescriptor for a single level of the input data
    call mkpio_iodesc_rawdata(mesh_i, 'PCT_GLC_GIC', pioid_i, pio_varid, pio_vartype, pio_iodesc, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! Create a route handle between the input and output mesh and get frac_o_nonorm
    allocate(frac_o_nonorm(ns_o),stat=ier)
    if (ier/=0) call shr_sys_abort(subname//" error in allocating frac_o_nonorm")
    ! Note that norm_by_fracs is false in the following because this routehandle is
    ! used to map fields that are expressed in terms of % of the grid cell.
    call create_routehandle_r8(mesh_i=mesh_i, mesh_o=mesh_o, norm_by_fracs=.false., &
         routehandle=routehandle_nonorm, frac_o=frac_o_nonorm, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMLogMemInfo("After create routehandle in "//trim(subname))

    ! Compute pctglcmec_gic_o, pctglcmec_gic_o, pctglcmec_icesheet_o and topoglcmec_unnorm_o
    ! Note that topoglcmec_unnorm_o is the average topographic height over glaciated areas -
    ! NOT the average topographic height of the entire grid cell
    do lev = 1, nlev
       write(6,'(i4)',advance='no') lev

       ! Read in one level of data
       rcode = pio_inq_varid(pioid_i, 'PCT_GLC_GIC', pio_varid)
       call pio_setframe(pioid_i, pio_varid, int(lev,kind=Pio_Offset_Kind))
       call mkpio_get_rawdata_level(pioid_i, pio_iodesc, lev, 'PCT_GLC_GIC', data_pctglc_gic_i)

       rcode = pio_inq_varid(pioid_i, 'PCT_GLC_ICESHEET', pio_varid)
       call pio_setframe(pioid_i, pio_varid, int(lev,kind=Pio_Offset_Kind))
       call mkpio_get_rawdata_level(pioid_i, pio_iodesc, lev, 'PCT_GLC_ICESHEET', data_pctglc_icesheet_i)

       ! Compute derived input data
       data_pctglc_i(:) = data_pctglc_gic_i(:) + data_pctglc_icesheet_i(:)

       ! Map level of data to output grid
       call regrid_rawdata(mesh_i, mesh_o, routehandle_nonorm, data_pctglc_i, data_pctglc_o, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call regrid_rawdata(mesh_i, mesh_o, routehandle_nonorm, data_pctglc_gic_i, data_pctglc_gic_o, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call regrid_rawdata(mesh_i, mesh_o, routehandle_nonorm, data_pctglc_icesheet_i, data_pctglc_icesheet_o, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       ! Compute output variables
       m = get_elevclass(topoice_i(lev))
       do no = 1,ns_o
          pctglcmec_o(no,m) = pctglcmec_o(no,m) + data_pctglc_o(no)
          if (outnc_3dglc) then
             pctglcmec_gic_o(no,m) = pctglcmec_gic_o(no,m)  + data_pctglc_gic_o(no)
             pctglcmec_icesheet_o(no,m) = pctglcmec_icesheet_o(no,m) + data_pctglc_icesheet_o(no)
          end if
          topoglcmec_unnorm_o(no,m) = topoglcmec_unnorm_o(no,m) + data_pctglc_o(no)*topoice_i(lev)
       end do
    end do

    ! Close glacier input file
    call pio_closefile(pioid_i)
    call ESMF_VMLogMemInfo("After pio_closefile in "//trim(subname))

    ! At this point, the various percentages are given as % of grid cell.
    ! We now renormalize these to be given as % of landunit.

    ! Normalize topoglcmec_o. To do this, note that pctglcmec_o(n,m) is equal to the sum of
    ! the weights used in doing the weighted average of topoice_i (weight = wt*pctglc_i/frac);
    ! hence pctglcmec_o(n,m) is the correct normalization factor
    do no = 1,ns_o
       do m = 1,nglcec
          if (pctglcmec_o(no,m) > 0) then
             topoglcmec_o(no,m) = topoglcmec_unnorm_o(no,m) / pctglcmec_o(no,m)
          else
             topoglcmec_o(no,m) = mean_elevation_vc(m)
          end if

          ! Correct for rounding errors that put topoglcmec_o(no,m) slightly outside the
          ! allowed bounds for this elevation class
          if (slightly_below(topoglcmec_o(no,m), elevclass(m))) then
             write(6,*) 'Warning: topoglcmec_o was slightly lower than lower bound; setting equal&
                  & to lower bound; for: ', no, m, topoglcmec_o(no,m), elevclass(m)
             write(6,*) '(this is informational only, and probably just indicates rounding error)'
             topoglcmec_o(no,m) = elevclass(m)
          else if (slightly_above(topoglcmec_o(no,m), elevclass(m+1))) then
             write(6,*) 'Warning: topoglcmec_o was slightly higher than upper bound; setting equal&
                  & to upper bound; for: ', no, m, topoglcmec_o(no,m), elevclass(m+1)
             write(6,*) '(this is informational only, and probably just indicates rounding error)'
             topoglcmec_o(no,m) = elevclass(m+1)
          end if
       end do
    end do

    ! Renormalize percentages to be given as % of landunit rather than % of grid cell.
    allocate(pctglc_tot_o(ns_o), stat=ier)
    if (ier/=0) call shr_sys_abort(subname//" error in allocating pctglc_tot")
    do no = 1,ns_o
       pctglc_tot_o(no) = sum(pctglcmec_o(no,:))
       if (pctglc_tot_o(no) > 0._r8) then
          pctglcmec_o(no,:) = pctglcmec_o(no,:) / pctglc_tot_o(no) * 100._r8
          if ( outnc_3dglc )then
             pctglcmec_gic_o(no,:) = pctglcmec_gic_o(no,:) / pctglc_tot_o(no) * 100._r8
             pctglcmec_icesheet_o(no,:) = pctglcmec_icesheet_o(no,:) / pctglc_tot_o(no) * 100._r8
          end if
       else
          ! Division of landunit is ambiguous. Apply the rule that all area is assigned to
          ! the lowest elevation class, and all GIC.
          pctglcmec_o(no,1) = 100._r8
          if ( outnc_3dglc ) then
             pctglcmec_gic_o(no,1) = 100._r8
          end if
       end if
    end do

    ! Set pctglc_gic_o to sum of pctglcmec_gic_o across elevation classes, and similarly for pctglc_icesheet_o
    if ( outnc_3dglc )then
       pctglc_gic_o      = sum(pctglcmec_gic_o, dim=2)
       pctglc_icesheet_o = sum(pctglcmec_icesheet_o, dim=2)
    end if

    ! --------------------------------------------------------------------
    ! Perform various sanity checks
    ! --------------------------------------------------------------------

    ! Confirm that the sum over pctglcmec_o (from 1 to nglcec) is 100%
    errors = .false.
    do no = 1,ns_o
       glc_sum = sum(pctglcmec_o(no,:))
       if (abs(glc_sum - 100._r8) > eps_small) then
          write(6,*)'glc_sum differs from 100% at no,pctglc= ',no,glc_sum
          errors = .true.
       end if
    end do

    ! Confirm that GIC + ICESHEET = 100%
    if ( outnc_3dglc )then
       do no = 1,ns_o
          if (abs((pctglc_gic_o(no) + pctglc_icesheet_o(no)) - 100._r8) > eps) then
             write(6,*)'GIC + ICESHEET differs from 100% at no,pctglc_gic,pctglc_icesheet =', &
                  no,pctglc_gic_o(no),pctglc_icesheet_o(no)
             errors = .true.
             ! TODO: output lat and lon out above
          end if
       end do

       ! Check that GIC + ICESHEET = total glacier at each elevation class
       do m = 1, nglcec
          do no = 1,ns_o
             if (abs((pctglcmec_gic_o(no,m) + pctglcmec_icesheet_o(no,m)) - &
                  pctglcmec_o(no,m)) > eps) then
                write(6,*)'GIC + ICESHEET differs from total GLC '
                write(6,*)'at no,m,pctglcmec,pctglcmec_gic,pctglcmec_icesheet = '
                write(6,*) no,m,pctglcmec_o(no,m),pctglcmec_gic_o(no,m),pctglcmec_icesheet_o(no,m)
                errors = .true.
             end if
          end do
       end do
    end if

    ! Error check: are all elevations within elevation class range
    do no = 1,ns_o
       do m = 1,nglcec
          if (topoglcmec_o(no,m) < elevclass(m) .or. topoglcmec_o(no,m) > elevclass(m+1)) then
             write(6,*) 'Error: mean elevation does not fall within elevation class '
             write(6,*) elevclass(m),elevclass(m+1),topoglcmec_o(no,m),m,no
             errors = .true.
          endif
       end do
    end do

    if (errors) then
       call shr_sys_abort(subname//" error in error checks")
    end if

    ! Output data tp fo;e
    if (root_task)  write(ndiag, '(a)') trim(subname)//" writing out pct_glc_mec"
    call mkfile_output(pioid_o, mesh_o, 'PCT_GLC_MEC', pctglcmec_o, lev1name='nglcec', rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in calling mkfile_output')

    if (root_task)  write(ndiag, '(a)') trim(subname)//" writing out topo_glc_mec"
    call mkfile_output(pioid_o, mesh_o, 'TOPO_GLC_MEC', topoglcmec_o, lev1name='nglcec', rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in calling mkfile_output')

    if (outnc_3dglc ) then
       if (root_task)  write(ndiag, '(a)') trim(subname)//" writing out pct_glc_mec_gic"
       call mkfile_output(pioid_o, mesh_o, 'PCT_GLC_MEC_GIC', pctglcmec_gic_o, lev1name='nglcec', rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in calling mkfile_output')

       if (root_task)  write(ndiag, '(a)') trim(subname)//" writing out pct_glc_mec_icesheet"
       call mkfile_output(pioid_o, mesh_o, 'PCT_GLC_MEC_ICESHEET', pctglcmec_icesheet_o, lev1name='nglcec', rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in calling mkfile_output')

       if (root_task)  write(ndiag, '(a)') trim(subname)//" writing out pct_glc_gic"
       call mkfile_output(pioid_o, mesh_o, 'PCT_GLC_GIC', pctglc_gic_o, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in calling mkfile_output')

       if (root_task)  write(ndiag, '(a)') trim(subname)//" writing out pct_icesheet"
       call mkfile_output(pioid_o, mesh_o, 'PCT_GLC_ICESHEET', pctglc_icesheet_o, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in calling mkfile_output')
    end if

    ! Deallocate dynamic memory
    call ESMF_RouteHandleDestroy(routehandle_nonorm, nogarbage = .true., rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort()
    call ESMF_MeshDestroy(mesh_i, nogarbage = .true., rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort()
    call ESMF_VMLogMemInfo("After destroy operations in "//trim(subname))

    if (root_task) then
       write (ndiag,'(a)') 'Successfully made percent elevation class and mean elevation for glaciers'
    end if

  end subroutine mkglcmec

  !=================================================================================
  subroutine mkglacier(file_mesh_i, file_data_i, mesh_o, glac_o, rc)
    !
    ! make percent glacier
    !
    ! In contrast to mkglcmec, this uses a "flat" PCT_GLACIER field (not separated by
    ! elevation class, and not separated into icesheet vs GIC).
    !
    ! This simpler routine is sufficient for cases when we run without multiple elevation
    ! classes. This routine is also used when running with multiple elevation classes: we
    ! first regrid the flat PCT_GLACIER field, then later create the multiple elevation class
    ! data. This multi-step process makes it easier to do corrections on the total
    ! PCT_GLACIER, and make sure these corrections apply appropriately to the multi-level
    ! output. The assumption is that PCT_GLACIER is the sum of both PCT_GLC_GIC and
    ! PCT_GLC_ICESHEET across all elevation bins.
    !
    ! input/output variables
    character(len=*) , intent(in)  :: file_mesh_i ! input mesh file name
    character(len=*) , intent(in)  :: file_data_i ! input data file name
    type(ESMF_Mesh)  , intent(in)  :: mesh_o      ! output mesh
    real(r8)         , intent(out) :: glac_o(:)   ! output grid: %glacier
    integer          , intent(out) :: rc
    !
    ! local variables
    type(ESMF_RouteHandle) :: routehandle_nonorm
    type(ESMF_Mesh)        :: mesh_i
    type(file_desc_t)      :: pioid
    integer                :: ni,no,k
    integer                :: ns_i, ns_o
    integer , allocatable  :: mask_i(:)
    real(r8), allocatable  :: frac_i(:)
    real(r8), allocatable  :: frac_o_nonorm(:)
    real(r8), allocatable  :: area_i(:)
    real(r8), allocatable  :: area_o(:)
    real(r8), allocatable  :: glac_i(:)            ! input grid: percent glac
    real(r8)               :: sum_fldi             ! global sum of dummy input fld
    real(r8)               :: sum_fldo             ! global sum of dummy output fld
    real(r8)               :: gglac_i              ! input  grid: global glac
    real(r8)               :: garea_i              ! input  grid: global area
    real(r8)               :: gglac_o              ! output grid: global glac
    real(r8)               :: garea_o              ! output grid: global area
    integer                :: ier, rcode           ! error status
    real(r8)               :: relerr = 0.00001     ! max error: sum overlap wts ne 1
    character(len=*), parameter :: subname = 'mkglacier'
    !-----------------------------------------------------------------------

    if (root_task) then
       write(ndiag,*)
       write(ndiag,'(1x,80a1)') ('=',k=1,80)
       write(ndiag,*)
       write(ndiag,'(a)') 'Attempting to make %glacier .....'
       write(ndiag,'(a)') ' Input file is '//trim(file_data_i)
       write(ndiag,'(a)') ' Input mesh file is '//trim(file_mesh_i)
    end if

    ! Open input data file
    rcode = pio_openfile(pio_iosystem, pioid, pio_iotype, trim(file_data_i), pio_nowrite)
    call ESMF_VMLogMemInfo("After pio_openfile "//trim(file_data_i))

    ! Read in input mesh
    mesh_i = ESMF_MeshCreate(filename=trim(file_mesh_i), fileformat=ESMF_FILEFORMAT_ESMFMESH, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMLogMemInfo("After create mesh_i in "//trim(subname))

    ! Determine ns_i
    call ESMF_MeshGet(mesh_i, numOwnedElements=ns_i, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Determine ns_o
    call ESMF_MeshGet(mesh_o, numOwnedElements=ns_o, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Get the landmask from the file and reset the mesh mask based on that
    allocate(frac_i(ns_i), stat=ier)
    if (ier/=0) call shr_sys_abort(subname//" ERROR in allocating frac_i")
    allocate(mask_i(ns_i), stat=ier)
    if (ier/=0) call shr_sys_abort(subname//" ERROR in allocating mask_i")
    call mkpio_get_rawdata(pioid, 'LANDMASK', mesh_i, frac_i, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    do ni = 1,ns_i
       if (frac_i(ni) > 0._r4) then
          mask_i(ni) = 1
       else
          mask_i(ni) = 0
       end if
    end do
    call ESMF_MeshSet(mesh_i, elementMask=mask_i, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! Read in glac_i
    allocate(glac_i(ns_i), stat=ier)
    if (ier/=0) call shr_sys_abort(subname//" error in allocating glac_i")
    call mkpio_get_rawdata(pioid, 'PCT_GLACIER', mesh_i, glac_i, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMLogMemInfo("After mkpio_getrawdata in "//trim(subname))

    ! Create a route handle between the input and output mesh and get frac_o_nonorm
    allocate(frac_o_nonorm(ns_o),stat=ier)
    if (ier/=0) call shr_sys_abort()
    ! Note that norm_by_fracs is false in the following because this routehandle is
    ! used to map fields that are expressed in terms of % of the grid cell.
    call create_routehandle_r8(mesh_i=mesh_i, mesh_o=mesh_o, norm_by_fracs=.false., &
         routehandle=routehandle_nonorm, frac_o=frac_o_nonorm, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMLogMemInfo("After create routehandle in "//trim(subname))

    ! Area-average percent cover on input grid to output grid (with correction for landmask)
    ! Note that percent cover is in terms of total grid area.
    ! Regrid glac_i to glac_o
    call regrid_rawdata(mesh_i, mesh_o, routehandle_nonorm, glac_i, glac_o, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    do no = 1,ns_o
       if (glac_o(no) < 1._r8) then
          glac_o(no) = 0._r8
       else if ((glac_o(no)) > 101._r8) then
          write(6,*) 'MKGLACIER error: glacier = ', glac_o(no), &
               ' > 101 for no = ', no
          call shr_sys_abort()
       else if ((glac_o(no)) > 100._r8) then
          if ((glac_o(no)) > 100.000001_r8) then
             write(6,*) 'MKGLACIER warning: glacier = ', glac_o(no), &
                ' > 100.000001 for no = ', no, ' Changing glacier > 100 to 100.'
          end if
          glac_o(no) = 100._r8
       end if
    enddo

    ! Check global areas
    call output_diagnostics_area(mesh_i, mesh_o, mask_i, frac_o_nonorm, &
         glac_i, glac_o, "pct glacier", percent=.true., ndiag=ndiag, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Deallocate dynamic memory
    deallocate (glac_i, frac_o_nonorm, mask_i)
    call ESMF_RouteHandleDestroy(routehandle_nonorm, nogarbage = .true., rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort()
    call ESMF_MeshDestroy(mesh_i, nogarbage = .true., rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort()

    if (root_task) then
       write (ndiag,'(a)') 'Successfully made %glacier'
    end if

  end subroutine mkglacier

  !=================================================================================
  integer function get_elevclass(topo, writewarn)
    !
    ! Returns elevation class index (1..nglcec) given the topographic height.
    ! If topo is lower than the lowest elevation class, returns 0.
    ! If topo is higher than the highest elevation class, returns (nglcec+1).
    ! In either of the two latter cases, the function also writes a warning message, unless
    ! writewarn is present and false.
    !
    ! input/output variables
    real(r4), intent(in) :: topo  ! topographic height (m)
    logical, intent(in), optional :: writewarn  ! should warning messages be written? (default: true)
    !
    ! local variables
    integer :: m
    logical :: my_writewarn
    character(len=*), parameter :: subname = 'get_elevclass'
    !-----------------------------------------------------------------------

    if (present(writewarn)) then
       my_writewarn = writewarn
    else
       my_writewarn = .true.
    end if

    if (topo < elevclass(1)) then
       if (my_writewarn) then
          write(6,*) 'WARNING in ', trim(subname)
          write(6,*) 'topo out of bounds'
          write(6,*) 'topo = ', topo
          write(6,*) 'elevclass(1) = ', elevclass(1)
       end if
       get_elevclass = 0
       return
    end if

    do m = 1, nglcec
       if (topo < elevclass(m+1)) then
          ! note that we already know that topo >= elevclass(m), otherwise we would have
          ! returned earlier
          get_elevclass = m
          return
       end if
    end do

    if (my_writewarn) then
       write(6,*) 'WARNING in ', trim(subname)
       write(6,*) 'topo out of bounds'
       write(6,*) 'topo = ', topo
       write(6,*) 'elevclass(nglcec+1) = ', elevclass(nglcec+1)
    end if
    get_elevclass = nglcec+1

  end function get_elevclass

  !=================================================================================
  real(r8) function mean_elevation_vc(class)
    !
    ! For a virtual column (thus, a column that has no true elevation data), return the
    ! "mean" elevation of the given elevation class.
    !
    ! input/output variables
    integer, intent(in) :: class  ! elevation class
    !
    ! local variables
    character(len=*), parameter :: subname = 'mean_elevation_vc'
    !-----------------------------------------------------------------------

    if (class < nglcec) then
       mean_elevation_vc = 0.5_r8 * (elevclass(class) + elevclass(class+1))
    else if (class == nglcec) then
       ! In the top elevation class; in this case, assignment of a "mean" elevation is
       ! somewhat arbitrary

       if (nglcec > 1) then
          mean_elevation_vc = 2.0_r8*elevclass(class) - elevclass(class-1)
       else
          ! entirely arbitrary
          mean_elevation_vc = 1000._r8
       end if
    else
       write(6,*) 'ERROR in ', trim(subname), ': class out of bounds= ', class
       call shr_sys_abort()
    end if

  end function mean_elevation_vc

end module mkglcmecMod
