module mkurbanparCommonMod
  
  !-----------------------------------------------------------------------
  ! Common routines for making urban parameter data, independent of the method used for
  ! making the urban parameters (e.g., averages, dominant type, etc.)
  !
  ! (WJS 4-18-12: In the past, this contained routines shared between mkurbanparDomMod and
  ! mkurbanparAvgMod; now there is just a single module, mkurbanparMod, but I am keeping the
  ! separate mkurbanparCommonMod in case a similar split comes back in the future. However,
  ! if such a split seems unlikely in the future, these routines could be moved back into
  ! mkurbanparMod.)
  !-----------------------------------------------------------------------

  use ESMF
  use pio
  use shr_kind_mod , only : r8 => shr_kind_r8, r4 => shr_kind_r4, cs => shr_kind_cs
  use shr_sys_mod  , only : shr_sys_abort
  use mkpioMod     , only : mkpio_get_rawdata, pio_iotype, pio_ioformat, pio_iosystem
  use mkpioMod     , only : mkpio_iodesc_output, mkpio_def_spatial_var, mkpio_wopen
  use mkpioMod     , only : mkpio_get_dimlengths, mkpio_get_rawdata
  use mkpioMod     , only : pio_iotype, pio_ioformat, pio_iosystem
  use mkesmfMod    , only : regrid_rawdata, create_routehandle_r8
  use mkutilsMod   , only : chkerr
  use mkvarctl     , only : ndiag

  implicit none
  private

  ! !public member functions:
  public :: mkurban_pct_diagnostics        ! print diagnostics related to pct urban
#ifdef TODO
  public :: mkurban_topo                   ! Get elevation to reduce urban for high elevation areas
#endif

  ! !public data members:
  real(r8), parameter :: MIN_DENS = 0.1_r8 ! minimum urban density (% of grid cell) - below this value, urban % is set to 0

  public :: MIN_DENS

  character(len=*) , parameter :: u_FILE_u = &
       __FILE__

!===============================================================
contains
!===============================================================

  subroutine mkurban_pct_diagnostics(area_i, area_o, mask_i, frac_o, urbn_i, urbn_o, dens_class)
    !
    ! print diagnostics related to pct urban
    ! Compare global areas on input and output grids
    !
    ! This is intended to be called after mkurban_pct, but is split out into a separate
    ! routine so that modifications to urbn_o can be made in between the two calls (e.g.,
    ! setting urbn_o to 0 wherever it is less than a certain threshold; the rules for doing
    ! this can't always be applied inline in mkurban_pct).
    !
    use mkvarpar

    ! input/output variables
    real(r8)          , intent(in) :: area_i(:)
    real(r8)          , intent(in) :: area_o(:)
    integer           , intent(in) :: mask_i(:)
    real(r8)          , intent(in) :: frac_o(:) 
    real(r8)          , intent(in) :: urbn_i(:)  ! input grid: percent urban
    real(r8)          , intent(in) :: urbn_o(:)  ! output grid: percent urban
    integer , intent(in), optional :: dens_class ! density class

    ! local variables:
    real(r8) :: gurbn_i ! input  grid: global urbn
    real(r8) :: garea_i ! input  grid: global area
    real(r8) :: gurbn_o ! output grid: global urbn
    real(r8) :: garea_o ! output grid: global area
    integer  :: ni,no,k ! indices
    character(len=*), parameter :: subname = 'mkurban_pct_diagnostics'
    !-----------------------------------------------------------------------

    ! Input grid
    gurbn_i = 0._r8
    garea_i = 0._r8
    do ni = 1, size(area_i)
       garea_i = garea_i + area_i(ni)*re**2
       gurbn_i = gurbn_i + urbn_i(ni)*(area_i(ni)/100._r8)* mask_i(ni)*re**2
    end do

    ! Output grid
    gurbn_o = 0._r8
    garea_o = 0._r8
    do no = 1, size(area_o)
       garea_o = garea_o + area_o(no)*re**2
       gurbn_o = gurbn_o + urbn_o(no)* (area_o(no)/100._r8)*frac_o(no)*re**2
    end do

    ! Diagnostic output
    write (ndiag,*)
    write (ndiag,'(1x,70a1)') ('=',k=1,70)
    if (present(dens_class)) then
       write (ndiag,'(1x,a,i0)') 'Urban Output -- class ', dens_class
    else
       write (ndiag,'(1x,a)') 'Urban Output'
    end if
    write (ndiag,'(1x,70a1)') ('=',k=1,70)
    write (ndiag,*)
    write (ndiag,'(1x,70a1)') ('.',k=1,70)
    write (ndiag,2001)
2001 format (1x,'surface type   input grid area  output grid area'/&
             1x,'                 10**6 km**2      10**6 km**2   ')
    write (ndiag,'(1x,70a1)') ('.',k=1,70)
    write (ndiag,*)
    write (ndiag,2003) gurbn_i*1.e-06,gurbn_o*1.e-06
    write (ndiag,2004) garea_i*1.e-06,garea_o*1.e-06
2002 format (1x,'urban       ',f14.3,f17.3)
2003 format (1x,'urban       ',f14.3,f22.8)
2004 format (1x,'all surface ',f14.3,f17.3)

  end subroutine mkurban_pct_diagnostics

#ifdef TODO
  !===============================================================
  subroutine mkurban_topo(mesh_i, mesh_o, file_data_i, varname, elev_o)
    !
    ! Make elevation data
    !
    use mkdiagnosticsMod, only : output_diagnostics_continuous
    !
    ! !ARGUMENTS:
    character(len=*)  , intent(in) :: datfname  ! input data file name
    integer           , intent(in) :: ndiag     ! unit number for diag out
    character(len=*)  , intent(in) :: varname   ! topo variable name
    real(r8)          , intent(out):: elev_o(:) ! output elevation data

    ! !LOCAL VARIABLES:
    real(r8), allocatable :: elev_i(:)  ! canyon_height to width ratio in
    real(r8), allocatable :: frac_o(:)  ! output fractions
    integer               :: ns_i,ns_o  ! indices
    integer               :: k,l,n,m,ni ! indices
    integer               :: ier        ! error status
    character(len=CS)     :: name       ! name of attribute
    character(len=CS)     :: unit       ! units of attribute
    character(len=*), parameter :: subname = 'mkelev'
    !-----------------------------------------------------------------------

    write (6,*) 'Attempting to make elevation .....'

    ! Query local mesh sizes
    call ESMF_MeshGet(mesh_i, numOwnedElements=ns_i, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_MeshGet(mesh_o, numOwnedElements=ns_o, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Read topo elev dataset with unit mask everywhere
    write (6,*) 'Open elevation file: ', trim(datfname)
    allocate(elev_i(ns_i), stat=ier)
    if (ier/=0) call shr_sys_abort()
    rcode = pio_openfile(pio_iosystem, pioid, pio_iotype, trim(file_data_i), pio_nowrite)
    call mkpio_get_rawdata(pioid, trim(varname), mesh_i, elev_i, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call pio_closefile(pioid)

    ! Create a route handle between the input and output mesh
    allocate(frac_o(ns_o), stat=ier)
    if (ier/=0) call shr_sys_abort()
    call create_routehandle_r8(mesh_i, mesh_o, routehandle, frac_o=frac_o, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMLogMemInfo("After create routehandle in "//trim(subname))

    ! Regrid input data to model resolution - determine elev_o on output grid
    elev_o(:) = 0.
    if (ier/=0) call shr_sys_abort()
    call regrid_rawdata(mesh_i, mesh_o, routehandle, elev_i, elev_o, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call output_diagnostics_continuous(elev_i, elev_o, tgridmap, &
         "Urban elev variable", "m", ndiag, mask_i, frac_o)

    ! Deallocate dynamic memory
    deallocate (elev_i)
    deallocate (frac_o)

    call ESMF_RouteHandleDestroy(routehandle, nogarbage = .true., rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort()

    if (root_task) then
       write (ndiag,'(a)') 'Successfully made elevation' 
       write (ndiag,'(a)')
    end if

  end subroutine mkurban_topo
#endif

end module mkurbanparCommonMod

