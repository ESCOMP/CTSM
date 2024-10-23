module mklaiMod

  !-----------------------------------------------------------------------
  ! Make LAI/SAI/height data
  !-----------------------------------------------------------------------

  use ESMF
  use pio
  use shr_kind_mod      , only : r8 => shr_kind_r8, r4=>shr_kind_r4
  use shr_sys_mod       , only : shr_sys_abort
  use mkpioMod          , only : pio_iotype, pio_ioformat, pio_iosystem
  use mkpioMod          , only : mkpio_get_rawdata, mkpio_get_rawdata_level
  use mkpioMod          , only : mkpio_iodesc_output, mkpio_iodesc_rawdata, mkpio_put_time_slice
  use mkesmfMod         , only : regrid_rawdata, create_routehandle_r8, get_meshareas
  use mkutilsMod        , only : chkerr
  use mkpftConstantsMod , only : c3cropindex, c3irrcropindex
  use mkvarctl          , only : root_task, ndiag, outnc_double, numpft, mpicom

  implicit none
  private

#include <mpif.h>

  public  :: mklai
  private :: check_global_sums

  character(len=*) , parameter :: u_FILE_u = &
       __FILE__

!=================================================================================
contains
!=================================================================================

  subroutine mklai(file_mesh_i, file_data_i, mesh_o, pioid_o, rc)
    !
    ! Make LAI/SAI/height data
    !
    ! input/output variables
    character(len=*)  , intent(in)    :: file_mesh_i ! input mesh file name
    character(len=*)  , intent(in)    :: file_data_i ! input data file name
    type(ESMF_Mesh)   , intent(in)    :: mesh_o      ! output mesh
    type(file_desc_t) , intent(inout) :: pioid_o
    integer           , intent(out)   :: rc
    !
    ! local variables
    type(ESMF_RouteHandle) :: routehandle           ! nearest neighbor routehandle
    type(ESMF_Mesh)        :: mesh_i
    type(file_desc_t)      :: pioid_i
    type(io_desc_t)        :: pio_iodesc_i
    type(io_desc_t)        :: pio_iotype_i
    type(var_desc_t)       :: pio_varid_i
    integer                :: pio_vartype_i
    type(io_desc_t)        :: pio_iodesc_o
    type(var_desc_t)       :: pio_varid_o
    integer                :: dimid
    integer                :: ni,no
    integer                :: ns_i, ns_o
    integer                :: k,l,m,nt              ! indices
    integer                :: numpft_i              ! number of plant types on input
    integer                :: ntime                  ! number of input time samples
    integer , allocatable  :: mask_i(:)
    real(r8), allocatable  :: frac_i(:)
    real(r8), allocatable  :: frac_o(:)
    real(r8), allocatable  :: data_i(:,:)
    real(r8), allocatable  :: data_o(:,:)
    real(r8), allocatable  :: mlai_o(:,:)           ! monthly lai
    real(r8), allocatable  :: msai_o(:,:)           ! monthly sai
    real(r8), allocatable  :: mhgtt_o(:,:)          ! monthly height (top)
    real(r8), allocatable  :: mhgtb_o(:,:)          ! monthly height (bottom)
    integer,  allocatable  :: laimask(:,:)          ! lai+sai output mask for each plant function type
    real(r8), allocatable  :: area_i(:)
    real(r8), allocatable  :: area_o(:)
    integer                :: ier, rcode            ! error status
    integer                :: xtype                 ! external type
    character(len=*), parameter :: subname = 'mklai'
    !-----------------------------------------------------------------------

    rc = ESMF_SUCCESS

    if (root_task) then
       write(ndiag,*)
       write(ndiag,'(a)') 'Attempting to make LAIs/SAIs/heights .....'
       write(ndiag,'(a)') ' Input file is '//trim(file_data_i)
       write(ndiag,'(a)') ' Input mesh file is '//trim(file_mesh_i)
    end if

    ! Open input data file
    rcode = pio_openfile(pio_iosystem, pioid_i, pio_iotype, trim(file_data_i), pio_nowrite)

    ! Read in input mesh
    call ESMF_VMLogMemInfo("Before create mesh_i in "//trim(subname))
    mesh_i = ESMF_MeshCreate(filename=trim(file_mesh_i), fileformat=ESMF_FILEFORMAT_ESMFMESH, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMLogMemInfo("After create mesh_i in "//trim(subname))

    ! Determine ns_i
    call ESMF_MeshGet(mesh_i, numOwnedElements=ns_i, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Determine ns_o
    call ESMF_MeshGet(mesh_o, numOwnedElements=ns_o, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Get the landmask from the input data file and reset the mesh mask based on that
    allocate(frac_i(ns_i), stat=ier)
    if (ier/=0) call shr_sys_abort()
    allocate(mask_i(ns_i), stat=ier)
    if (ier/=0) call shr_sys_abort()
    call mkpio_get_rawdata(pioid_i, 'LANDMASK', mesh_i, frac_i, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    do ni = 1,ns_i
       if (frac_i(ni) > 0._r8) then
          mask_i(ni) = 1
       else
          mask_i(ni) = 0
       end if
    end do
    deallocate(frac_i)
    call ESMF_MeshSet(mesh_i, elementMask=mask_i, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! Create a route handle between the input and output mesh
    allocate(frac_o(ns_o))
    call create_routehandle_r8(mesh_i=mesh_i, mesh_o=mesh_o, norm_by_fracs=.true., &
         routehandle=routehandle, frac_o=frac_o, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMLogMemInfo("After create routehandle in "//trim(subname))

    rcode = pio_inq_dimid(pioid_i, 'pft', dimid)
    rcode = pio_inq_dimlen(pioid_i, dimid, numpft_i)
    rcode = pio_inq_dimid(pioid_i, 'time', dimid)
    rcode = pio_inq_dimlen(pioid_i, dimid, ntime)

    if (numpft_i /= numpft+1) then
       if (root_task) then
          write(ndiag,*) 'WARNING: ' // trim(subname) // '(): parameter numpft+1 = ', numpft+1, &
               'does not equal input dataset numpft = ', numpft_i
          write(ndiag,*)'This inconsistency used to stop the program. Now we allow it '
          write(ndiag,*)'because crop pfts 17-last are assumed to never use satellite lai data.'
       end if
       if (numpft_i > numpft + 1) then
          ! NOTE(bja, 2015-01) If this error check is determined to be
          ! invalid, all the loop bounds over output data in this
          ! routine will need to be double checked!
          if (root_task) then
             write(ndiag,*) (subname) //' error input numpft must be less than or equal to output numpft+1.'
          end if
          call shr_sys_abort()
       end if
    endif
    if (ntime /= 12) then
       if (root_task) then
          write(ndiag,*) subname // 'error must have 12 time samples on input data'
       end if
       call shr_sys_abort()
    endif

    ! Dynamic allocation of variables of size 0:numpft
    allocate(mlai_o(ns_o,0:numpft),  &
             msai_o(ns_o,0:numpft),  &
             mhgtt_o(ns_o,0:numpft), &
             mhgtb_o(ns_o,0:numpft), &
             laimask(ns_i,0:numpft), stat=ier )
    if (ier /= 0) then
       call shr_sys_abort(subname //' mklai allocation error ')
    end if

    ! Create iodescriptor for a single level of the input data
    call mkpio_iodesc_rawdata(mesh_i, 'MONTHLY_LAI', pioid_i, pio_varid_i, pio_vartype_i, pio_iodesc_i, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! Create iodescriptor for a single level of the output data
    call mkpio_iodesc_output(pioid_o, mesh_o, 'MONTHLY_LAI', pio_iodesc_o, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in generating an iodesc for MONTHLY_LAI')

    ! Allocate memory that will be used in time loop below
    allocate(data_i(0:numpft_i,ns_i),stat=ier)
    if (ier/=0) call shr_sys_abort()
    allocate(data_o(0:numpft_i,ns_o),stat=ier)
    if (ier/=0) call shr_sys_abort()

    ! The following is needed for the global check
    allocate(area_i(ns_i), stat=ier)
    if (ier/=0) call shr_sys_abort()
    allocate(area_o(ns_o), stat=ier)
    if (ier/=0) call shr_sys_abort()
    call get_meshareas(mesh_i, area_i, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call get_meshareas(mesh_o, area_o, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! Loop over the 12 months and write out data
    mlai_o(:,:)  = 0.
    msai_o(:,:)  = 0.
    mhgtt_o(:,:) = 0.
    mhgtb_o(:,:) = 0.

    do nt = 1, ntime

       ! time is months for LAI, SAI, and pft heights
       rcode = pio_inq_varid(pioid_o, 'time', pio_varid_o)
       rcode = pio_put_var(pioid_o, pio_varid_o, (/nt/), nt)

       ! Below - copy LAI, SAI, & heights from the C3 crop (pft15)
       ! to the irrigated (pft16) whether crop is on or off
       ! Hence loop to numpft_i - 1 for other pfts

       ! Read in one time slice of data for mlai, regrid and write out
       rcode = pio_inq_varid(pioid_i, 'MONTHLY_LAI', pio_varid_i)
       call pio_setframe(pioid_i, pio_varid_i, int(nt, kind=Pio_Offset_Kind))
       call mkpio_get_rawdata_level(pioid_i, pio_iodesc_i, nt, 'MONTHLY_LAI', data_i)
       call regrid_rawdata(mesh_i, mesh_o, routehandle, data_i, data_o, 0, numpft_i, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       do l = 0,numpft_i-1
          do no = 1,ns_o
             mlai_o(no,l) = data_o(l,no)
          end do
       end do
       mlai_o(:,c3irrcropindex)  = mlai_o(:,c3cropindex)
       rcode = pio_inq_varid(pioid_o, 'MONTHLY_LAI', pio_varid_o)
       call mkpio_put_time_slice(pioid_o, pio_varid_o, pio_iodesc_o, nt, mlai_o)
       call check_global_sums('LAI', ns_i, ns_o, numpft_i, nt, &
            data_i, data_o, area_i, area_o, mask_i, frac_o)

       ! Read in one time slice of data for msai, regrid and write out
       rcode = pio_inq_varid(pioid_i, 'MONTHLY_SAI', pio_varid_i)
       call pio_setframe(pioid_i, pio_varid_i, int(nt, kind=Pio_Offset_Kind))
       call mkpio_get_rawdata_level(pioid_i, pio_iodesc_i, nt, 'MONTHLY_SAI', data_i)
       call regrid_rawdata(mesh_i, mesh_o, routehandle, data_i, data_o, 0, numpft_i, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       do l = 0,numpft_i-1
          do no = 1,ns_o
             msai_o(no,l) = data_o(l,no)
          end do
       end do
       msai_o(:,c3irrcropindex)  = msai_o(:,c3cropindex)
       rcode = pio_inq_varid(pioid_o, 'MONTHLY_SAI', pio_varid_o)
       call mkpio_put_time_slice(pioid_o, pio_varid_o, pio_iodesc_o, nt, msai_o)

       ! Read in one time slice of data for msai, regrid and write out
       rcode = pio_inq_varid(pioid_i, 'MONTHLY_HEIGHT_TOP', pio_varid_i)
       call pio_setframe(pioid_i, pio_varid_i, int(nt, kind=Pio_Offset_Kind))
       call mkpio_get_rawdata_level(pioid_i, pio_iodesc_i, nt, 'MONTHLY_HEIGHT_TOP', data_i)
       call regrid_rawdata(mesh_i, mesh_o, routehandle, data_i, data_o, 0, numpft_i, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       do l = 0,numpft_i-1
          do no = 1,ns_o
             mhgtt_o(no,l) = data_o(l,no)
          end do
       end do
       mhgtt_o(:,c3irrcropindex) = mhgtt_o(:,c3cropindex)
       rcode = pio_inq_varid(pioid_o, 'MONTHLY_HEIGHT_TOP', pio_varid_o)
       call mkpio_put_time_slice(pioid_o, pio_varid_o, pio_iodesc_o, nt, mhgtt_o)

       ! Read in one time slice of data for msai, regrid and write out
       rcode = pio_inq_varid(pioid_i, 'MONTHLY_HEIGHT_BOT', pio_varid_i)
       call pio_setframe(pioid_i, pio_varid_i, int(nt, kind=Pio_Offset_Kind))
       call mkpio_get_rawdata_level(pioid_i, pio_iodesc_i, nt, 'MONTHLY_HEIGHT_BOT', data_i)
       call regrid_rawdata(mesh_i, mesh_o, routehandle, data_i, data_o, 0, numpft_i, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       do l = 0,numpft_i-1
          do no = 1,ns_o
             mhgtb_o(no,l) = data_o(l,no)
          end do
       end do
       mhgtb_o(:,c3irrcropindex) = mhgtb_o(:,c3cropindex)
       rcode = pio_inq_varid(pioid_o, 'MONTHLY_HEIGHT_BOT', pio_varid_o)
       call mkpio_put_time_slice(pioid_o, pio_varid_o, pio_iodesc_o, nt, mhgtb_o)

       if (root_task) then
          write (ndiag,*) 'Successfully made LAIs/SAIs/heights for month ', nt
       end if

    end do  ! end loop over months
    call pio_syncfile(pioid_o)

    ! Free the decomps and close the file
    call pio_freedecomp(pioid_o, pio_iodesc_o)
    call pio_freedecomp(pioid_i, pio_iodesc_i)
    call pio_closefile(pioid_i)
    call ESMF_VMLogMemInfo("After pio_closefile for input in "//trim(subname))

    ! Release memory
    call ESMF_RouteHandleDestroy(routehandle, nogarbage = .true., rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort()
    call ESMF_MeshDestroy(mesh_i, nogarbage = .true., rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort()
    call ESMF_VMLogMemInfo("After destroy operations in "//trim(subname))

  end subroutine mklai

  !=================================================================================
  subroutine check_global_sums(name, ns_i, ns_o, numpft_i, nt, &
       data_i, data_o, area_i, area_o, mask_i, frac_o)

    ! Compare global areas on input and output grids
    ! NB. data_i and data_o started at 0 outside this subroutine but start
    ! at 1 within the subroutine, so the loops in the subroutine start at 1

    ! input/otuput variables
    character(len=*) , intent(in) :: name
    integer          , intent(in) :: ns_i
    integer          , intent(in) :: ns_o
    integer          , intent(in) :: nt
    integer          , intent(in) :: numpft_i
    real(r8)         , intent(in) :: data_i(:,:)
    real(r8)         , intent(in) :: data_o(:,:)
    real(r8)         , intent(in) :: area_i(:)
    real(r8)         , intent(in) :: area_o(:)
    integer          , intent(in) :: mask_i(:)
    real(r8)         , intent(in) :: frac_o(:)

    ! local variables
    integer  :: ni, no, l, k
    integer  :: ier
    real(r8) :: local_i(numpft_i)  ! local global area, by surface type
    real(r8) :: local_o(numpft_i)  ! local global area, by surface type
    real(r8) :: global_i(numpft_i)   ! input grid: global area pfts
    real(r8) :: global_o(numpft_i)   ! output grid: global area pfts
    !-----------------------------------------------------------------------

    ! Input grid global area
    local_i(:) = 0.
    do l = 1, numpft_i
       do ni = 1, ns_i
          local_i(l)  = local_i(l) + data_i(l,ni) *area_i(ni)*mask_i(ni)
       end do
    end do
    call mpi_reduce(local_i, global_i , numpft_i, MPI_REAL8, MPI_SUM, 0, mpicom, ier)

    ! Output grid global area
    local_o(:) = 0.
    do l = 1, numpft_i
       do no = 1, ns_o
          local_o(l) = local_o(l) + data_o(l,no) *area_o(no)*frac_o(no)
       end do
    end do
    call mpi_reduce(local_o, global_o , numpft_i, MPI_REAL8, MPI_SUM, 0, mpicom, ier)

    ! Comparison
    if (root_task) then
       write (ndiag,*)
       write (ndiag,*) trim(name)//' Output for month ',nt
       write (ndiag,'(1x,70a1)') ('.',k=1,70)
       write (ndiag,101)
101    format (1x,'PFT input grid area output grid area',/ &
               1x,3x,'     10**6 km**2','      10**6 km**2')
       write (ndiag,'(1x,70a1)') ('.',k=1,70)
       write (ndiag,*)
       do l = 1, numpft_i
          write (ndiag,102) l-1, global_i(l)*1.e-06*1.e-02, global_o(l)*1.e-06*1.e-02
102       format (1x,i3,f16.3,f17.3)
       end do
    end if

  end subroutine check_global_sums

end module mklaiMod
