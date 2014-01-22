module dynpftFileMod

  !---------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Handle reading of the pftdyn dataset, which specifies transient areas of natural PFTs
  !
  ! !USES:
  use clmtype
  use shr_kind_mod          , only : r8 => shr_kind_r8
  use decompMod             , only : bounds_type, BOUNDS_LEVEL_PROC
  use dynFileMod            , only : dyn_file_type
  use dynVarTimeInterpMod   , only : dyn_var_time_interp_type
  use clm_varctl            , only : iulog
  use abortutils            , only : endrun
  use spmdMod               , only : masterproc
  use shr_assert_mod        , only : shr_assert
  use shr_log_mod           , only : errMsg => shr_log_errMsg
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  implicit none
  private
  save
  public :: dynpft_init     ! initialize information read from pftdyn dataset
  public :: dynpft_interp   ! interpolate pftdyn information to current time step
  !
  ! ! PRIVATE TYPES
  type(dyn_file_type), target    :: dynpft_file   ! information for the pftdyn file
  type(dyn_var_time_interp_type) :: wtpft         ! weight of each pft relative to the natural veg landunit
  !---------------------------------------------------------------------------

contains


  !-----------------------------------------------------------------------
  subroutine dynpft_init(bounds)
    !
    ! !DESCRIPTION:
    ! Initialize dynamic pft dataset (position it to the right time samples
    ! that bound the initial model date)
    !
    ! This also calls dynpft_interp for the initial time
    !
    ! !USES:
    use clm_varctl  , only : fpftdyn
    use clm_varpar  , only : numpft, maxpatch_pft, natpft_lb, natpft_ub, natpft_size
    use clm_varcon  , only : numurbl, istsoil, istcrop
    use clm_varsur  , only : pctspec, wt_lunit
    use ncdio_pio
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds  ! proc-level bounds
    !
    ! !LOCAL VARIABLES:
    integer  :: n,g,nl                          ! indices
    integer  :: ier                             ! error status
    logical  :: readvar	                        ! true => variable is on input dataset
    integer  :: wtpft_shape(2)                  ! shape of the wtpft data
    ! leave the following as pointers instead of changing to allocatable since the ncd_io
    ! interface expects a pointer type
    real(r8) , pointer     :: pctnatveg (:)        ! percent of gcell is natural vegetation landunit 
    real(r8) , pointer     :: pctcrop (:)          ! percent of gcell is crop landunit 
    real(r8) , pointer     :: pctgla (:)           ! percent of gcell is glacier 
    real(r8) , pointer     :: pctlak (:)           ! percent of gcell is lake 
    real(r8) , pointer     :: pctwet (:)           ! percent of gcell is wetland 
    real(r8) , pointer     :: pcturb (:,:)         ! percent of gcell is urbanized 
    real(r8) , pointer     :: pcturb_tot (:)       ! percent of grid cell is urban (sum over density classes) 
    character(len= 32)     :: subname='dynpft_init'! subroutine name
    !-----------------------------------------------------------------------

    call shr_assert(bounds%level == BOUNDS_LEVEL_PROC, subname // ': argument must be PROC-level bounds')

    ! Error check

    if ( maxpatch_pft /= numpft+1 )then
       call endrun( subname//' maxpatch_pft does NOT equal numpft+1 -- this is invalid for dynamic PFT case' )
    end if

    allocate(pctnatveg(bounds%begg:bounds%endg))
    allocate(pctcrop(bounds%begg:bounds%endg))
    allocate(pctgla(bounds%begg:bounds%endg))
    allocate(pctlak(bounds%begg:bounds%endg))
    allocate(pctwet(bounds%begg:bounds%endg))
    allocate(pcturb(bounds%begg:bounds%endg,numurbl))
    allocate(pcturb_tot(bounds%begg:bounds%endg))

    ! pctspec must be saved between time samples
    ! position to first time sample - assume that first time sample must match starting date
    ! check consistency -  special landunits, grid, frac and mask
    ! only do this once

    if (masterproc) then
       write(iulog,*) 'Attempting to read pft dynamic landuse data .....'
    end if

    dynpft_file = dyn_file_type(fpftdyn)

    ! Consistency check
    
    call check_dim(dynpft_file, 'natpft', natpft_size)

    call ncd_io(ncid=dynpft_file, varname='PCT_NATVEG', flag='read', data=pctnatveg, &
         dim1name=grlnd, readvar=readvar)
    if (.not. readvar) call endrun( trim(subname)//' ERROR: PCT_NATVEG NOT on pftdyn file' )

    call ncd_io(ncid=dynpft_file, varname='PCT_CROP', flag='read', data=pctcrop, &
         dim1name=grlnd, readvar=readvar)
    if (.not. readvar) call endrun( trim(subname)//' ERROR: PCT_CROP NOT on pftdyn file' )

    call ncd_io(ncid=dynpft_file, varname='PCT_WETLAND', flag='read', data=pctwet, &
         dim1name=grlnd, readvar=readvar)
    if (.not. readvar) call endrun( trim(subname)//' ERROR: PCT_WETLAND NOT on pftdyn file' )

    call ncd_io(ncid=dynpft_file, varname= 'PCT_LAKE', flag='read', data=pctlak, &
         dim1name=grlnd, readvar=readvar)
    if (.not. readvar) call endrun( trim(subname)//' ERROR: PCT_LAKE NOT on pftdyn file' )

    call ncd_io(ncid=dynpft_file, varname= 'PCT_GLACIER', flag='read', data=pctgla, &
         dim1name=grlnd, readvar=readvar)
    if (.not. readvar) call endrun( trim(subname)//' ERROR: PCT_GLACIER NOT on pftdyn file' )

    call ncd_io(ncid=dynpft_file, varname= 'PCT_URBAN'  , flag='read', data=pcturb, &
         dim1name=grlnd, readvar=readvar)
    if (.not. readvar) call endrun( trim(subname)//' ERROR: PCT_URBAN NOT on pftdyn file' )
    pcturb_tot(bounds%begg : bounds%endg) = 0._r8
    do n = 1, numurbl
       do nl = bounds%begg,bounds%endg
          pcturb_tot(nl) = pcturb_tot(nl) + pcturb(nl,n)
       enddo
    enddo

    ! Consistency check
    do g = bounds%begg,bounds%endg
       ! WJS (5-9-13): I am adding these pctnatveg and pctcrop consistency checks for now,
       ! although both these and the following pctspec consistency check (which was
       ! already here) may not be necessary now that pctpft is specified as % of landunit
       ! rather than % of grid cell. Furthermore, once we have dynamic landunits, both of
       ! these consistency checks may become problematic, and could be removed. If these
       ! consistency checks are removed, I think we could do quite a bit of cleanup:
       ! (1) I think that the time-invariant PCT fields no longer need to be on the pftdyn dataset
       ! (2) Similarly, much of this routine could be removed
       ! (3) I think that pctspec no longer needs to be saved (I think it's just saved
       !     for the sake of this consistency check)
       ! (4) wt_lunit no longer needs to be saved past the end of the initialize1 routine
       !     in clm_initializeMod (so we can move its deallocation from initialize2 to
       !     initialize1)
       if (abs(pctnatveg(g) - wt_lunit(g,istsoil)*100._r8) > 1.e-13_r8) then
          write(iulog,*) subname//'mismatch between input PCT_NATVEG = ', pctnatveg(g), &
               ' and that obtained from surface dataset ', wt_lunit(g,istsoil)*100._r8, &
               ' at g= ',g
          call endrun()
       end if
       if (abs(pctcrop(g) - wt_lunit(g,istcrop)*100._r8) > 1.e-13_r8) then
          write(iulog,*) subname//'mismatch between input PCT_CROP = ', pctcrop(g), &
               ' and that obtained from surface dataset ', wt_lunit(g,istcrop)*100._r8, &
               ' at g= ',g
          call endrun()
       end if

    !   This was causing a fail, even though values are the same to within 1e-15
    !   if (pctlak(g)+pctwet(g)+pcturb(g)+pctgla(g) /= pctspec(g)) then 
       if (abs((pctlak(g)+pctwet(g)+pcturb_tot(g)+pctgla(g))-pctspec(g)) > 1e-13_r8) then 
          write(iulog,*) subname//'mismatch between input pctspec = ',&
                     pctlak(g)+pctwet(g)+pcturb_tot(g)+pctgla(g),&
                    ' and that obtained from surface dataset ', pctspec(g),' at g= ',g
           call endrun()
       end if
    end do

    ! read data PCT_NAT_PFT corresponding to correct year
    !
    ! Note: if you want to change PCT_NAT_PFT so that it is NOT interpolated, but instead
    ! jumps to each year's value on Jan 1 of that year, simply change wtpft to be of type
    ! dyn_var_time_uninterp_type (rather than dyn_var_time_interp_type), and change the
    ! following constructor to construct a variable of dyn_var_time_uninterp_type. That's
    ! all you need to do.

    wtpft_shape = [(bounds%endg-bounds%begg+1), natpft_size]
    wtpft = dyn_var_time_interp_type( &
         dyn_file=dynpft_file, varname='PCT_NAT_PFT', &
         dim1name=grlnd, conversion_factor=100._r8, &
         do_check_sums_equal_1=.true., data_shape=wtpft_shape)

    call dynpft_interp(bounds)

    deallocate(pctnatveg, pctcrop, pctgla,pctlak,pctwet,pcturb,pcturb_tot)

  end subroutine dynpft_init

  !-----------------------------------------------------------------------
  subroutine dynpft_interp(bounds)
    !
    ! !DESCRIPTION:
    ! Time interpolate dynamic pft data to get pft weights for model time
    !
    ! !USES:
    use clm_varcon      , only : istsoil
    use clm_varpar      , only : natpft_lb, natpft_ub
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds  ! proc-level bounds
    !
    ! !LOCAL VARIABLES:
    integer               :: m,p,l,g          ! indices
    real(r8), allocatable :: wtpft_cur(:,:)   ! current pft weights
    character(len=32) :: subname='dynpft_interp' ! subroutine name
    !-----------------------------------------------------------------------

    ! Interpolate pctpft to current time step - output in pctpft_intp
    ! Map interpolated pctpft to subgrid weights
    ! assumes that maxpatch_pft = numpft + 1, that each landunit has only 1 column, 
    ! SCAM and CNDV have not been defined, and create_croplandunit = .false.

    call shr_assert(bounds%level == BOUNDS_LEVEL_PROC, subname // ': argument must be PROC-level bounds')

    ! Get pft weights for this time step

    call dynpft_file%update_time_info()

    allocate(wtpft_cur(bounds%begg:bounds%endg, natpft_lb:natpft_ub))
    call wtpft%get_current_data(wtpft_cur)

    do p = bounds%begp,bounds%endp
       g = pft%gridcell(p)
       l = pft%landunit(p)

       ! Note that we only deal with the istsoil landunit here, NOT the istcrop landunit
       ! (if there is one)
       ! (However, currently [as of 5-9-13] the code won't let you run with transient
       ! PFTs combined with create_crop_landunit anyway, so it's a moot point.)
       if (lun%itype(l) == istsoil) then
          m = pft%itype(p)

          ! Note that the following assignment assumes that all PFTs share a single column
          pft%wtcol(p) = wtpft_cur(g,m)
       end if

    end do

    deallocate(wtpft_cur)

  end subroutine dynpft_interp

end module dynpftFileMod
