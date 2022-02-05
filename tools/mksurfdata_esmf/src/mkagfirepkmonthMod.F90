module mkagfirepkmonthMod

  !-----------------------------------------------------------------------
  ! Make agricultural fire peak month data
  !-----------------------------------------------------------------------

  use ESMF
  use pio
  use shr_kind_mod   , only : r8 => shr_kind_r8, r4=>shr_kind_r4
  use shr_sys_mod    , only : shr_sys_abort
  use mkpioMod       , only : mkpio_get_rawdata, mkpio_get_dimlengths
  use mkpioMod       , only : pio_iotype, pio_ioformat, pio_iosystem
  use mkpioMod       , only : mkpio_iodesc_rawdata, mkpio_get_rawdata_level
  use mkesmfMod      , only : regrid_rawdata, create_routehandle_r4, get_meshareas
  use mkutilsMod     , only : chkerr
  use mkvarctl       , only : ndiag, root_task
  use mkvarpar       , only : re
  use mkchecksMod    , only : min_bad, max_bad

  implicit none
  private           ! By default make data private

#include <mpif.h>

  public  :: mkagfirepkmon       ! Set agricultural fire peak month
  private :: define_months       ! define month strings

  integer , parameter :: min_valid_value = 1
  integer , parameter :: max_valid_value = 12
  integer , parameter :: unsetmon        = 13   ! flag to indicate agricultural fire peak month NOT set

  character(len=*) , parameter :: u_FILE_u = &
       __FILE__

!===============================================================
contains
!===============================================================

  subroutine mkagfirepkmon(file_mesh_i, file_data_i, mesh_o, agfirepkmon_o, rc)
    !
    ! Make agricultural fire peak month data from higher resolution data
    !
    ! input/output variables
    character(len=*) , intent(in)    :: file_mesh_i      ! input mesh file name
    character(len=*) , intent(in)    :: file_data_i      ! input data file name
    type(ESMF_Mesh)  , intent(in)    :: mesh_o           ! output mesh
    integer          , intent(inout) :: agfirepkmon_o(:) ! agricultural fire peak month
    integer          , intent(out)   :: rc
    !
    ! local variables:
    type(ESMF_RouteHandle)         :: routehandle
    type(ESMF_Mesh)                :: mesh_i
    type(file_desc_t)              :: pioid
    integer                        :: ni,no,l
    integer                        :: ns_i, ns_o
    integer , allocatable          :: mask_i(:)
    real(r4), allocatable          :: frac_i(:)
    real(r4), allocatable          :: frac_o(:)
    real(r8), allocatable          :: area_i(:)
    real(r8), allocatable          :: area_o(:)
    real(r4), allocatable          :: data_i(:,:)
    real(r4), allocatable          :: data_o(:,:)
    integer                        :: max_index(1)
    integer , allocatable          :: agfirepkmon_i(:)  ! input grid: agricultural fire peak month
    integer                        :: nagfirepkmon      ! number of peak months
    character(len=35), allocatable :: month(:)          ! name of each month
    integer                        :: rcode, ier        ! error status
    integer, parameter             :: unsetmon = 13     ! flag to indicate agricultural fire peak month NOT set
    integer, parameter             :: miss = unsetmon   ! missing data indicator
    integer, parameter             :: min_valid = 1     ! minimum valid value
    integer, parameter             :: max_valid = 13    ! maximum valid value
    real(r8), allocatable          :: loc_gast_i(:)     ! local global area, by surface type
    real(r8), allocatable          :: loc_gast_o(:)     ! local global area, by surface type
    real(r8), allocatable          :: gast_i(:)         ! global area, by surface type
    real(r8), allocatable          :: gast_o(:)         ! global area, by surface type
    character(len=*), parameter    :: subname = 'mkagfirepkmon'
    !-----------------------------------------------------------------------

    if (root_task) then
       write(ndiag,*)
       write(ndiag,'(a)') 'Attempting to make agricultural fire peak month data .....'
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
       if (frac_i(ni) > 0.) then
          mask_i(ni) = 1
       else
          mask_i(ni) = 0
       end if
    end do
    call ESMF_MeshSet(mesh_i, elementMask=mask_i, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! Read in agfirepkmon_i
    allocate(agfirepkmon_i(ns_i), stat=ier)
    if (ier/=0) call shr_sys_abort(subname//" error in allocating agfirepkmon_i")
    call mkpio_get_rawdata(pioid, 'abm', mesh_i, agfirepkmon_i, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMLogMemInfo("After mkpio_getrawdata in "//trim(subname))

    ! Create a route handle between the input and output mesh and get frac_o
    allocate(frac_o(ns_o),stat=ier)
    if (ier/=0) call shr_sys_abort()
    call create_routehandle_r4(mesh_i, mesh_o, routehandle, frac_o=frac_o, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMLogMemInfo("After create routehandle in "//trim(subname))

    ! Now determine data_i as a real 2d array - for every possible soil color create a global
    ! field with gridcells equal to 1 for that soil color and zero elsewhere
    allocate(data_i(min_valid_value:max_valid,ns_i),stat=ier)
    if (ier/=0) call shr_sys_abort(subname//' error in allocating data_i')
    data_i(:,:) = 0._r4
    do l = min_valid,max_valid
       do ni = 1,ns_i
          if (int(agfirepkmon_i(ni)) == l .and. int(agfirepkmon_i(ni)) /= miss) then
             data_i(l,ni) = 1._r4 * mask_i(ni)
          end if
       end do
    end do

    ! Regrid data_i to data_o
    ! Note that any input point that is outside the range [min_valid_value,max_valid_value]
    ! will be ignored; this ignores input points with value of unsetmon
    allocate(data_o(min_valid:max_valid, ns_o),stat=ier)
    if (ier/=0) call shr_sys_abort(subname//' error in allocating data_o')
    call regrid_rawdata(mesh_i, mesh_o, routehandle, data_i, data_o, min_valid, max_valid, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_LogWrite(subname//'after regrid rawdata in '//trim(subname))

    do no = 1,ns_o
       max_index = maxloc(data_o(:,no))
       write(6,*)'DEBUG: no,max_index = ',no,max_index
       agfirepkmon_o(no) = max_index(1)
    end do

    ! Check validity of output data
    if (min_bad(agfirepkmon_o, min_valid, 'agfirepkmon') .or. &
        max_bad(agfirepkmon_o, max_valid, 'agfirepkmon')) then
       call shr_sys_abort()
    end if

    ! -----------------------------------------------------------------
    ! Output diagnostics comparing global area of each peak month on input and output grids
    !
    ! WJS (3-4-13): I am trying to generally put these diagnostics in mkdiagnosticsMod, but
    ! so far there isn't a general diagnostics routine for categorical data
    !
    ! TODO(wjs, 2016-01-22) Now there is a routine for this: output_diagnostics_index.
    ! However, it currently doesn't provide the capability for named months. Either add
    ! that capability or decide it's not important, then delete the below code, instead
    ! calling output_diagnostics_index.
    ! -----------------------------------------------------------------

    !     nagfirepkmon = maxval(agfirepkmon_i)
    !     allocate(gast_i(1:nagfirepkmon),gast_o(1:nagfirepkmon),month(1:nagfirepkmon))
    !     call define_months(nagfirepkmon, month)
    
    !     gast_i(:) = 0.0_r8
    !     do ni = 1,ns_i
    !        k = agfirepkmon_i(ni)
    !        gast_i(k) = gast_i(k) + area_src(ni)*mask(ni)*re**2
    !     end do
    !     gast_o(:) = 0.0_r8
    !     do no = 1,ns_o
    !        k = agfirepkmon_o(no)
    !        gast_o(k) = gast_o(k) + area_dst(no)*frac_dst(no)*re**2
    !     end do
    
    !     area comparison
    
    !     write (ndiag,*)
    !     write (ndiag,'(1x,70a1)') ('=',k=1,70)
    !     write (ndiag,*) 'Agricultural fire peak month Output'
    !     write (ndiag,'(1x,70a1)') ('=',k=1,70)
    
    !     write (ndiag,*)
    !     write (ndiag,'(1x,70a1)') ('.',k=1,70)
    !     write (ndiag,1001)
    ! 1001 format (1x,'peak month',20x,' input grid area output grid area',/ &
    !          1x,33x,'     10**6 km**2','      10**6 km**2')
    !     write (ndiag,'(1x,70a1)') ('.',k=1,70)
    !     write (ndiag,*)
    
    !     do k = 1, nagfirepkmon
    !        write (ndiag,1002) month(k),gast_i(k)*1.e-6,gast_o(k)*1.e-6
    ! 1002   format (1x,a35,f16.3,f17.3)
    !     end do
    
    ! Close the file 
    call pio_closefile(pioid)

    ! Release memory
    call ESMF_RouteHandleDestroy(routehandle, nogarbage = .true., rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort()
    call ESMF_MeshDestroy(mesh_i, nogarbage = .true., rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort()
    !call ESMF_VMLogMemInfo("After destroy operations in "//trim(subname))

    if (root_task) then
       write (ndiag,'(a)') 'Successfully made Agricultural fire peak month'
    end if
    write(6,*)'finished mkagfirepkmon'

  end subroutine mkagfirepkmon

  !===============================================================
  subroutine define_months(nagfirepkmon, month)
    !
    ! Define month strings
    !
    ! input/output variables
    integer         , intent(in) :: nagfirepkmon    ! max input value (including the 'unset' special value)
    character(len=*), intent(out):: month(:)        ! name of each month value
    !-----------------------------------------------------------------------

    if (nagfirepkmon == unsetmon) then
       if (size(month) < 13) then
          write(6,*) 'month array too small: ', size(month), ' < 13'
          call shr_sys_abort()
       end if
       month(1)  = 'January                             '
       month(2)  = 'February                            '
       month(3)  = 'March                               '
       month(4)  = 'April                               '
       month(5)  = 'May                                 '
       month(6)  = 'June                                '
       month(7)  = 'July                                '
       month(8)  = 'August                              '
       month(9)  = 'September                           '
       month(10) = 'October                             '
       month(11) = 'November                            '
       month(12) = 'December                            '
       month(13) = 'no agricultural fire peak month data'
    else
       write(6,*)'nagfirepkmon value of ',nagfirepkmon,' not supported'
       call shr_sys_abort()
    end if

  end subroutine define_months

end module mkagfirepkmonthMod
