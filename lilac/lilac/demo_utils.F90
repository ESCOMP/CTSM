module demo_utils
!----------------------------------------------------------------------------
    use mpi           ,  only : MPI_COMM_WORLD, MPI_COMM_NULL, MPI_Init, MPI_FINALIZE, MPI_SUCCESS
    use spmdMod       ,  only : masterproc
    implicit none
    private
    public :: demo_init
    public :: read_netcdf_mesh
    integer                            :: ierr
    integer                                         :: COMP_COMM
    integer                            :: npts   ! domain global size
    integer                            :: num_local
!----------------------------------------------------------------------------
contains
!----------------------------------------------------------------------------
    subroutine demo_init(gindex_atm)
        !! TODO: IS THE INTENT CORRECT FOR GINDEX_ATM
        integer , allocatable, intent(inout)  :: gindex_atm(:)
        integer                :: ntasks
        integer                :: mytask
        !-----------------------------------------------------------------------------
        ! Initiallize MPI
        !-----------------------------------------------------------------------------

        npts = 3312
        ! this is coming from
        ! /glade/work/mvertens/ctsm.nuopc/cime/src/drivers/nuopc/drivers/cime/esmApp.F90
        call MPI_init(ierr)
        COMP_COMM = MPI_COMM_WORLD

        !https://github.com/yudong-tian/LIS-CLM4.5SP/blob/8cec515a628325c73058cfa466db63210cd562ac/xlis-bld/xlis_main.F90
        if (ierr .ne. MPI_SUCCESS) then
            print *,'Error starting MPI program. Terminating.'
            call MPI_ABORT(MPI_COMM_WORLD, ierr)
        end if

        !

        call MPI_COMM_RANK(COMP_COMM, mytask, ierr)
        call MPI_COMM_SIZE(COMP_COMM, ntasks, ierr)

        if (masterproc) then
            print *, "MPI initialization done ..., ntasks=", ntasks
        end if
        
        call decompInit_atm( ntasks, mytask, gindex_atm)
        print *, "gindex_atm for ", mytask,"is: ", gindex_atm
        print *, "size gindex_atm for ", mytask,"is: ", size(gindex_atm)
    end subroutine demo_init

    subroutine decompInit_atm( ntasks, mytask,  gindex_atm)

        ! !DESCRIPTION:

        ! !USES:

        ! !ARGUMENTS:
        integer , intent(in)               :: ntasks
        integer , intent(in)               :: mytask
        integer , allocatable, intent(out) :: gindex_atm(:) ! this variable is allocated here, and is assumed to start unallocated
        ! !LOCAL VARIABLES:
        integer                            :: my_start
        integer                            :: my_end
        integer                            :: i_local
        integer                            :: i_global
        !------------------------------------------------------------------------------
        ! create the a global index array for ocean points

        num_local = npts / ntasks

        my_start = num_local*mytask + min(mytask, mod(npts, ntasks)) + 1
        ! The first mod(npts,ntasks) of ntasks are the ones that have an extra point
        if (mytask < mod(npts, ntasks)) then
           num_local = num_local + 1
        end if
        my_end = my_start + num_local - 1

        allocate(gindex_atm(num_local))

        i_global = my_start
        do i_local = 1, num_local
            gindex_atm(i_local) = i_global
            i_global = i_global +1
        end do 

    end subroutine decompInit_atm

    subroutine read_netcdf_mesh(filename)

        use netcdf
        implicit none

    !
    !  Parameters
    !

    !
    !  Arguments | Global Variables 
    !
        character(*)                     , intent(in)  :: filename


    !
    !  Local Variables
    !

        integer :: idfile

        integer :: ierror
        integer :: dimid_node
        integer :: dimid_elem
        integer :: dimid_maxnodepe
        integer :: dimid_coordDim

        integer :: iddim_node
        integer :: iddim_elem
        integer :: iddim_maxnodepe
        integer :: iddim_coordDim

        integer :: idvar_nodeCoords
        integer :: idvar_CenterCoords

        character (len=100) :: string

        integer :: nnode
        integer :: nelem
        integer :: maxnodePE
        integer :: coordDim
        !-----------------------------------------------------------------------------
        ! Open mesh file and get the idfile
        ierror  = nf90_open ( filename, NF90_NOWRITE, idfile); call nc_check_err(ierror, "opening file", filename)

        ! Get the dimid of  dimensions 
        ierror  = nf90_inq_dimid(idfile, 'nodeCount'        , dimid_node ); call nc_check_err(ierror, "inq_dimid nodeCount", filename)
        ierror  = nf90_inq_dimid(idfile, 'elementCount'     , dimid_elem ); call nc_check_err(ierror, "inq_dimid elementCount", filename)
        ierror  = nf90_inq_dimid(idfile, 'maxNodePElement'  , dimid_maxnodepe ); call nc_check_err(ierror, "inq_dimid maxNodePElement", filename)
        ierror  = nf90_inq_dimid(idfile, 'coordDim'         , dimid_coordDim  ); call nc_check_err(ierror, "coordDim", filename)

        ! Inquire dimensions based on their dimeid(s)
        ierror = nf90_inquire_dimension(idfile, iddim_node        , string, nnode     ); call nc_check_err(ierror, "inq_dim nodeCount", filename)
        ierror = nf90_inquire_dimension(idfile, iddim_elem        , string, nelem     ); call nc_check_err(ierror, "inq_dim elementCount", filename)
        ierror = nf90_inquire_dimension(idfile, iddim_maxnodepe   , string, maxnodePE ); call nc_check_err(ierror, "inq_dim maxNodePElement", filename)
        ierror = nf90_inquire_dimension(idfile, iddim_coordDim    , string, coordDim  ); call nc_check_err(ierror, "inq_dim coordDim", filename)


        ! Get variable IDs (varid)
        ierror = nf90_inq_varid(idfile, 'nodeCoords'   , idvar_nodeCoords      ); call nc_check_err(ierror, "inq_varid nodeCoords", filename)
        ierror = nf90_inq_varid(idfile, 'CenterCoords' , idvar_CenterCoords   ); call nc_check_err(ierror, "inq_varid CenterCoords", filename)

        ! Get variables values from varids
        !ierror = nf90_get_var(idfile, idvar_nodeCoords       , nodeCoords       , start=(/ 1,1/)   , count=(/ nnode, coordDim /)     ); call nc_check_err(ierror,"get_var nodeCoords", filename)
        !ierror = nf90_get_var(idfile, idvar_CenterCoords     , CenterCoords     , start=(/ 1,1/)   , count=(/ nelem, coordDim /)     ); call nc_check_err(ierror,"get_var CenterCoords", filename)




    end subroutine read_netcdf_mesh

end module demo_utils

