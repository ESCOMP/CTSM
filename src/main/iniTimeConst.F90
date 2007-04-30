#include <misc.h>
#include <preproc.h>

!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: iniTimeConst
!
! !INTERFACE:
subroutine iniTimeConst
!
! !DESCRIPTION:
! Initialize time invariant clm variables
! 1) removed references to shallow lake - since it is not used
! 2) ***Make c%z, c%zi and c%dz allocatable depending on if you
!    have lake or soil
! 3) rootfr only initialized for soil points
!
! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use nanMod
  use clmtype
  use decompMod   , only : get_proc_bounds, get_proc_global
  use decompMod   , only : gsMap_lnd_gdc2glo, perm_lnd_gdc2glo
  use clm_atmlnd  , only : clm_a2l
  use clm_varpar  , only : nlevsoi, nlevlak, lsmlon, lsmlat, numpft, numrad
  use clm_varcon  , only : istice, istdlak, istwet, isturb, &
                           zlak, dzlak, zsoi, dzsoi, zisoi, spval, &
                           albsat, albdry
  use clm_varctl  , only : nsrest, fsurdat,scmlon,scmlat,single_column
  use pftvarcon   , only : ncorn, nwheat, noveg, ntree, roota_par, rootb_par,  &
                           smpso, smpsc, fnitr, &
                           z0mr, displar, dleaf, rhol, rhos, taul, taus, xl, &
                           qe25, vcmx25, mp, c3psn, slatop, dsladlai, leafcn, flnr, woody, &
                           lflitcn, frootcn, livewdcn, deadwdcn, froot_leaf, stem_leaf, croot_stem, &
                           flivewd, fcur, lf_flab, lf_fcel, lf_flig, fr_flab, fr_fcel, fr_flig, &
                           dw_fcel, dw_flig, leaf_long, evergreen, stress_decid, season_decid, &
                           resist, &
                           pftpar , tree   , summergreen, raingreen  , sla     , &
                           lm_sapl, sm_sapl, hm_sapl    , rm_sapl    , latosa  , &
                           allom1 , allom2 , allom3     , reinickerp , wooddens
  use clm_time_manager, only : get_step_size
  use abortutils  , only : endrun
  use fileutils   , only : getfil
  use ndepFileMod , only : ndeprd
  use pftvarcon   , only : pftconrd
  use ncdio
  use spmdMod
!
! !ARGUMENTS:
  implicit none
!
! !CALLED FROM:
! subroutine initialize in module initializeMod.
!
! !REVISION HISTORY:
! Created by Gordon Bonan.
! Updated to clm2.1 data structrues by Mariana Vertenstein
! 4/26/05, Peter Thornton: Eliminated exponential decrease in saturated hydraulic
!   conductivity (hksat) with depth. 
!
! !LOCAL VARIABLES:
!
! local pointers to implicit in arguments
!
  integer , pointer :: ivt(:)             !  vegetation type index
  integer , pointer :: pcolumn(:)         ! column index of corresponding pft
  integer , pointer :: pgridcell(:)       ! gridcell index of corresponding pft
  integer , pointer :: clandunit(:)       ! landunit index of column
  integer , pointer :: cgridcell(:)       ! gridcell index of column
  integer , pointer :: ltype(:)           ! landunit type index
!
! local pointers to implicit out arguments
!
  real(r8), pointer :: z(:,:)             ! layer depth (m)
  real(r8), pointer :: zi(:,:)            ! interface level below a "z" level (m)
  real(r8), pointer :: dz(:,:)            ! layer thickness depth (m)
  real(r8), pointer :: rootfr(:,:)        ! fraction of roots in each soil layer
  real(r8), pointer :: rresis(:,:)        !root resistance by layer (0-1)  (nlevsoi)	
  real(r8), pointer :: dewmx(:)           ! maximum allowed dew [mm]
  real(r8), pointer :: bsw(:,:)           ! Clapp and Hornberger "b" (nlevsoi)  
  real(r8), pointer :: bsw2(:,:)          ! Clapp and Hornberger "b" for CN code
  real(r8), pointer :: psisat(:,:)        ! soil water potential at saturation for CN code (MPa)
  real(r8), pointer :: vwcsat(:,:)        ! volumetric water content at saturation for CN code (m3/m3)
  real(r8), pointer :: watsat(:,:)        ! volumetric soil water at saturation (porosity) (nlevsoi) 
  real(r8), pointer :: watdry(:,:)        ! btran parameter for btran=0
  real(r8), pointer :: watopt(:,:)        ! btran parameter for btran = 1
  real(r8), pointer :: hksat(:,:)         ! hydraulic conductivity at saturation (mm H2O /s) (nlevsoi) 
  real(r8), pointer :: sucsat(:,:)        ! minimum soil suction (mm) (nlevsoi) 
  real(r8), pointer :: csol(:,:)          ! heat capacity, soil solids (J/m**3/Kelvin) (nlevsoi) 
  real(r8), pointer :: tkmg(:,:)          ! thermal conductivity, soil minerals  [W/m-K] (new) (nlevsoi) 
  real(r8), pointer :: tkdry(:,:)         ! thermal conductivity, dry soil (W/m/Kelvin) (nlevsoi) 
  real(r8), pointer :: tksatu(:,:)        ! thermal conductivity, saturated soil [W/m-K] (new) (nlevsoi) 
  real(r8), pointer :: wtfact(:)          ! maximum saturated fraction for a gridcell
  real(r8), pointer :: smpmin(:)          ! restriction for min of soil potential (mm) (new)
  real(r8), pointer :: hkdepth(:)         ! decay factor (m)
  integer , pointer :: isoicol(:)         ! soil color class
  real(r8), pointer :: gwc_thr(:)         ! threshold soil moisture based on clay content
  real(r8), pointer :: mss_frc_cly_vld(:) ! [frc] Mass fraction clay limited to 0.20
  real(r8), pointer :: forc_ndep(:)       ! nitrogen deposition rate (gN/m2/s)
#if (defined CASA)
  real(r8), pointer :: sandfrac(:)
  real(r8), pointer :: clayfrac(:)
#endif
!
!EOP
!
! !OTHER LOCAL VARIABLES:
  integer  :: ncid             ! netCDF file id 
  integer  :: n,i,j,ib,lev,bottom      ! indices
  integer  :: g,l,c,p          ! indices
  integer  :: m                ! vegetation type index
  real(r8) :: bd               ! bulk density of dry soil material [kg/m^3]
  real(r8) :: tkm              ! mineral conductivity
  real(r8) :: xksat            ! maximum hydraulic conductivity of soil [mm/s]
  real(r8) :: scalez = 0.025_r8   ! Soil layer thickness discretization (m)
  real(r8) :: clay,sand        ! temporaries
  real(r8) :: slope,intercept        ! temporary, for rooting distribution
  integer  :: begp, endp       ! per-proc beginning and ending pft indices
  integer  :: begc, endc       ! per-proc beginning and ending column indices
  integer  :: begl, endl       ! per-proc beginning and ending landunit indices
  integer  :: begg, endg       ! per-proc gridcell ending gridcell indices
  integer  :: numg             ! total number of gridcells across all processors
  integer  :: numl             ! total number of landunits across all processors
  integer  :: numc             ! total number of columns across all processors
  integer  :: nump             ! total number of pfts across all processors
  real(r8),pointer :: arrayl(:)   ! generic global array
  integer ,pointer :: irrayg(:)   ! generic global array
  integer ,pointer :: soic2d(:)   ! read in - soil color
  real(r8),pointer :: sand3d(:,:) ! read in - soil texture: percent sand
  real(r8),pointer :: clay3d(:,:) ! read in - soil texture: percent clay
  real(r8),pointer :: ndep(:)     ! read in - annual nitrogen deposition rate (gN/m2/yr)
  real(r8),pointer :: gti(:)      ! read in - fmax
  integer  :: start(3),count(3)   ! netcdf start/count arrays
  integer  :: dimid,varid      ! netCDF id's
  integer  :: ret, time_index
  real(r8) :: nlevsoidata(nlevsoi)
  integer  :: ier                                ! error status
  character(len=256) :: locfn                    ! local filename
  character(len= 32) :: subname = 'iniTimeConst' ! subroutine name
  integer :: mxsoil_color                        ! maximum number of soil color classes

  integer :: closelatidx,closelonidx
  real(r8):: closelat,closelon

!------------------------------------------------------------------------

  if (masterproc) write (6,*) 'Attempting to initialize time invariant variables'

  call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)
  call get_proc_global(numg, numl, numc, nump)

  allocate(soic2d(begg:endg),ndep(begg:endg), gti(begg:endg))
  allocate(sand3d(begg:endg,nlevsoi),clay3d(begg:endg,nlevsoi))

  ! Assign local pointers to derived subtypes components (landunit-level)

  ltype           => clm3%g%l%itype

  ! Assign local pointers to derived subtypes components (column-level)

  clandunit       => clm3%g%l%c%landunit
  cgridcell       => clm3%g%l%c%gridcell
  z               => clm3%g%l%c%cps%z
  dz              => clm3%g%l%c%cps%dz
  zi              => clm3%g%l%c%cps%zi
  bsw             => clm3%g%l%c%cps%bsw
  bsw2            => clm3%g%l%c%cps%bsw2
  psisat          => clm3%g%l%c%cps%psisat
  vwcsat          => clm3%g%l%c%cps%vwcsat
  watsat          => clm3%g%l%c%cps%watsat
  watdry          => clm3%g%l%c%cps%watdry  
  watopt          => clm3%g%l%c%cps%watopt  
  hksat           => clm3%g%l%c%cps%hksat
  sucsat          => clm3%g%l%c%cps%sucsat
  tkmg            => clm3%g%l%c%cps%tkmg
  tksatu          => clm3%g%l%c%cps%tksatu
  tkdry           => clm3%g%l%c%cps%tkdry
  csol            => clm3%g%l%c%cps%csol
  smpmin          => clm3%g%l%c%cps%smpmin
  hkdepth         => clm3%g%l%c%cps%hkdepth
  wtfact          => clm3%g%l%c%cps%wtfact
  isoicol         => clm3%g%l%c%cps%isoicol
  gwc_thr         => clm3%g%l%c%cps%gwc_thr
  mss_frc_cly_vld => clm3%g%l%c%cps%mss_frc_cly_vld
  forc_ndep       => clm_a2l%forc_ndep

  ! Assign local pointers to derived subtypes components (pft-level)

  ivt             => clm3%g%l%c%p%itype
  pgridcell       => clm3%g%l%c%p%gridcell
  pcolumn         => clm3%g%l%c%p%column
  dewmx           => clm3%g%l%c%p%pps%dewmx
  rootfr          => clm3%g%l%c%p%pps%rootfr
  rresis          => clm3%g%l%c%p%pps%rresis
#if (defined CASA)
  sandfrac        => clm3%g%l%c%p%pps%sandfrac
  clayfrac        => clm3%g%l%c%p%pps%clayfrac
#endif

!  ! Determine necessary subgrid bounds
!
!  call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)
!  call get_proc_global(numg, numl, numc, nump)

  ! --------------------------------------------------------------------
  ! Read soil color, sand and clay from surface dataset 
  ! --------------------------------------------------------------------

  if (masterproc) then
     write (6,*) 'Attempting to read soil color, sand and clay boundary data .....'
     call getfil (fsurdat, locfn, 0)
     call check_ret(nf_open(locfn, 0, ncid), subname)

     ! Determine number of soil color classes - if number of soil color classes is not
     ! on input dataset set it to 8

     ier = nf_inq_varid(ncid, 'mxsoil_color', varid)
     if (ier == NF_NOERR) then
        call check_ret(nf_inq_varid(ncid, 'mxsoil_color', varid), subname)
        call check_ret(nf_get_var_int(ncid, varid, mxsoil_color), subname)
     else
        mxsoil_color = 8  
     end if
  endif
  call mpi_bcast( mxsoil_color,        1  , MPI_INTEGER, 0, mpicom, ier )
  count(1) = lsmlon
  count(2) = lsmlat
  if (single_column) then
     call scam_setlatlonidx(ncid,scmlat,scmlon,closelat,closelon,closelatidx,closelonidx)
     start(1) = closelonidx
     start(2) = closelatidx
  else
     start(1) = 1
     start(2) = 1
  end if
  start(3) = 1
  count(3) = 1

  ! Read fmax
  call ncd_iolocal(ncid, 'FMAX', 'read', gti, begg, endg, gsMap_lnd_gdc2glo, perm_lnd_gdc2glo,start(:2),count(:2))

  ! Read in soil color, sand and clay fraction

  call ncd_iolocal(ncid, 'SOIL_COLOR', 'read', soic2d, begg, endg, gsMap_lnd_gdc2glo, perm_lnd_gdc2glo,start(:2),count(:2))

  allocate(arrayl(begg:endg))
  do n = 1,nlevsoi
     start(3) = n
     call ncd_iolocal(ncid,'PCT_SAND','read',arrayl,begg,endg,gsMap_lnd_gdc2glo,perm_lnd_gdc2glo,start,count)
     sand3d(begg:endg,n) = arrayl(begg:endg)
     call ncd_iolocal(ncid,'PCT_CLAY','read',arrayl,begg,endg,gsMap_lnd_gdc2glo,perm_lnd_gdc2glo,start,count)
     clay3d(begg:endg,n) = arrayl(begg:endg)
  enddo
  deallocate(arrayl)

  if (masterproc) then
     call check_ret(nf_close(ncid), subname)
     write (6,*) 'Successfully read fmax, soil color, sand and clay boundary data'
     write (6,*)
  endif
  ! Determine saturated and dry soil albedos for n color classes and 
  ! numrad wavebands (1=vis, 2=nir)

  allocate(albsat(mxsoil_color,numrad), albdry(mxsoil_color,numrad), stat=ier)
  if (ier /= 0) then
     write (6,*)'iniTimeConst: allocation error for albsat, albdry'
     call endrun()
  end if

  if (mxsoil_color == 8) then
     albsat(1:8,1) = (/0.12_r8,0.11_r8,0.10_r8,0.09_r8,0.08_r8,0.07_r8,0.06_r8,0.05_r8/)
     albsat(1:8,2) = (/0.24_r8,0.22_r8,0.20_r8,0.18_r8,0.16_r8,0.14_r8,0.12_r8,0.10_r8/)
     albdry(1:8,1) = (/0.24_r8,0.22_r8,0.20_r8,0.18_r8,0.16_r8,0.14_r8,0.12_r8,0.10_r8/)
     albdry(1:8,2) = (/0.48_r8,0.44_r8,0.40_r8,0.36_r8,0.32_r8,0.28_r8,0.24_r8,0.20_r8/)
  else if (mxsoil_color == 20) then
     albsat(1:20,1) = (/0.25_r8,0.23_r8,0.21_r8,0.20_r8,0.19_r8,0.18_r8,0.17_r8,0.16_r8,&
                        0.15_r8,0.14_r8,0.13_r8,0.12_r8,0.11_r8,0.10_r8,0.09_r8,0.08_r8,0.07_r8,0.06_r8,0.05_r8,0.04_r8/)
     albsat(1:20,2) = (/0.50_r8,0.46_r8,0.42_r8,0.40_r8,0.38_r8,0.36_r8,0.34_r8,0.32_r8,&
                        0.30_r8,0.28_r8,0.26_r8,0.24_r8,0.22_r8,0.20_r8,0.18_r8,0.16_r8,0.14_r8,0.12_r8,0.10_r8,0.08_r8/)
     albdry(1:20,1) = (/0.36_r8,0.34_r8,0.32_r8,0.31_r8,0.30_r8,0.29_r8,0.28_r8,0.27_r8,&
                        0.26_r8,0.25_r8,0.24_r8,0.23_r8,0.22_r8,0.20_r8,0.18_r8,0.16_r8,0.14_r8,0.12_r8,0.10_r8,0.08_r8/)
     albdry(1:20,2) = (/0.61_r8,0.57_r8,0.53_r8,0.51_r8,0.49_r8,0.48_r8,0.45_r8,0.43_r8,&
                        0.41_r8,0.39_r8,0.37_r8,0.35_r8,0.33_r8,0.31_r8,0.29_r8,0.27_r8,0.25_r8,0.23_r8,0.21_r8,0.16_r8/)
  else
     write(6,*)'maximum color class = ',mxsoil_color,' is not supported'
     call endrun
  end if
  
#if (defined CASA)
!dir$ concurrent
!cdir nodep
  do p = begp,endp
     g = pgridcell(p)
     sandfrac(p) = sand3d(g,1)/100.0_r8
     clayfrac(p) = clay3d(g,1)/100.0_r8
  end do
#endif

  ! --------------------------------------------------------------------
  ! Read list of PFTs and their corresponding parameter values
  ! This is independent of the model resolution
  ! --------------------------------------------------------------------

  call pftconrd()

  ! --------------------------------------------------------------------
  ! If a nitrogen deposition dataset has been specified, read it
  ! --------------------------------------------------------------------
  
  call ndeprd(ndep)

  ! --------------------------------------------------------------------
  ! Initialize time constant arrays of ecophysiological constants and
  ! arrays of dgvm ecophysiological constants
  ! --------------------------------------------------------------------

!dir$ concurrent
!cdir nodep
   do m = 0,numpft
      pftcon%ncorn(m) = ncorn
      pftcon%nwheat(m) = nwheat
      pftcon%noveg(m) = noveg
      pftcon%ntree(m) = ntree
      pftcon%z0mr(m) = z0mr(m)
      pftcon%displar(m) = displar(m)
      pftcon%dleaf(m) = dleaf(m)
      pftcon%xl(m) = xl(m)
      do ib = 1,numrad
         pftcon%rhol(m,ib) = rhol(m,ib)
         pftcon%rhos(m,ib) = rhos(m,ib)
         pftcon%taul(m,ib) = taul(m,ib)
         pftcon%taus(m,ib) = taus(m,ib)
      end do
      pftcon%qe25(m) = qe25(m)
      pftcon%vcmx25(m) = vcmx25(m)
      pftcon%mp(m) = mp(m)
      pftcon%c3psn(m) = c3psn(m)
      pftcon%sla(m) = sla(m)
      pftcon%slatop(m) = slatop(m)
      pftcon%dsladlai(m) = dsladlai(m)
      pftcon%leafcn(m) = leafcn(m)
      pftcon%flnr(m) = flnr(m)
      pftcon%smpso(m) = smpso(m)
      pftcon%smpsc(m) = smpsc(m)
      pftcon%fnitr(m) = fnitr(m)
      pftcon%woody(m) = woody(m)
      pftcon%lflitcn(m) = lflitcn(m)
      pftcon%frootcn(m) = frootcn(m)
      pftcon%livewdcn(m) = livewdcn(m)
      pftcon%deadwdcn(m) = deadwdcn(m)
      pftcon%froot_leaf(m) = froot_leaf(m)
      pftcon%stem_leaf(m) = stem_leaf(m)
      pftcon%croot_stem(m) = croot_stem(m)
      pftcon%flivewd(m) = flivewd(m)
      pftcon%fcur(m) = fcur(m)
      pftcon%lf_flab(m) = lf_flab(m)
      pftcon%lf_fcel(m) = lf_fcel(m)
      pftcon%lf_flig(m) = lf_flig(m)
      pftcon%fr_flab(m) = fr_flab(m)
      pftcon%fr_fcel(m) = fr_fcel(m)
      pftcon%fr_flig(m) = fr_flig(m)
      pftcon%dw_fcel(m) = dw_fcel(m)
      pftcon%dw_flig(m) = dw_flig(m)
      pftcon%leaf_long(m) = leaf_long(m)
      pftcon%evergreen(m) = evergreen(m)
      pftcon%stress_decid(m) = stress_decid(m)
      pftcon%season_decid(m) = season_decid(m)
      pftcon%resist(m) = resist(m)
   end do

!dir$ concurrent
!cdir nodep
   do m = 0,numpft
      dgv_pftcon%respcoeff(m) = pftpar(m,5)
      dgv_pftcon%flam(m) = pftpar(m,6)
      dgv_pftcon%resist(m) = pftpar(m,8)
      dgv_pftcon%l_turn(m) = pftpar(m,9)
      dgv_pftcon%l_long(m) = pftpar(m,10)
      dgv_pftcon%s_turn(m) = pftpar(m,11)
      dgv_pftcon%r_turn(m) = pftpar(m,12)
      dgv_pftcon%l_cton(m) = pftpar(m,13)
      dgv_pftcon%s_cton(m) = pftpar(m,14)
      dgv_pftcon%r_cton(m) = pftpar(m,15)
      dgv_pftcon%l_morph(m) = pftpar(m,16)
      dgv_pftcon%l_phen(m) = pftpar(m,17)
      dgv_pftcon%lmtorm(m) = pftpar(m,18)
      dgv_pftcon%crownarea_max(m) = pftpar(m,20)
      dgv_pftcon%init_lai(m) = pftpar(m,21)
      dgv_pftcon%x(m) = pftpar(m,22)
      dgv_pftcon%tcmin(m) = pftpar(m,28)
      dgv_pftcon%tcmax(m) = pftpar(m,29)
      dgv_pftcon%gddmin(m) = pftpar(m,30)
      dgv_pftcon%twmax(m) = pftpar(m,31)
      dgv_pftcon%lm_sapl(m) = lm_sapl(m)
      dgv_pftcon%sm_sapl(m) = sm_sapl(m)
      dgv_pftcon%hm_sapl(m) = hm_sapl(m)
      dgv_pftcon%rm_sapl(m) = rm_sapl(m)
      dgv_pftcon%tree(m) = tree(m)
      dgv_pftcon%summergreen(m) = summergreen(m)
      dgv_pftcon%raingreen(m) = raingreen(m)
      dgv_pftcon%reinickerp(m) = reinickerp
      dgv_pftcon%wooddens(m) = wooddens
      dgv_pftcon%latosa(m) = latosa
      dgv_pftcon%allom1(m) = allom1
      dgv_pftcon%allom2(m) = allom2
      dgv_pftcon%allom3(m) = allom3
   end do

   ! --------------------------------------------------------------------
   ! Define layer structure for soil and lakes
   ! Vertical profile of snow is initialized in routine iniTimeVar
   ! --------------------------------------------------------------------

   ! check that lake and soil levels are the same for now

   if (nlevlak /= nlevsoi) then
      write(6,*)'number of soil levels and number of lake levels must be the same'
      write(6,*)'nlevsoi= ',nlevsoi,' nlevlak= ',nlevlak
      call endrun
   endif

   ! Lake layers (assumed same for all lake patches)

   dzlak(1) = 0.1_r8
   dzlak(2) = 1._r8
   dzlak(3) = 2._r8
   dzlak(4) = 3._r8
   dzlak(5) = 4._r8
   dzlak(6) = 5._r8
   dzlak(7) = 7._r8
   dzlak(8) = 7._r8
   dzlak(9) = 10.45_r8
   dzlak(10)= 10.45_r8

   zlak(1) =  0.05_r8
   zlak(2) =  0.6_r8
   zlak(3) =  2.1_r8
   zlak(4) =  4.6_r8
   zlak(5) =  8.1_r8
   zlak(6) = 12.6_r8
   zlak(7) = 18.6_r8
   zlak(8) = 25.6_r8
   zlak(9) = 34.325_r8
   zlak(10)= 44.775_r8

   ! Soil layers and interfaces (assumed same for all non-lake patches)
   ! "0" refers to soil surface and "nlevsoi" refers to the bottom of model soil

   do j = 1, nlevsoi
      zsoi(j) = scalez*(exp(0.5_r8*(j-0.5_r8))-1._r8)    !node depths
   enddo

   dzsoi(1) = 0.5_r8*(zsoi(1)+zsoi(2))             !thickness b/n two interfaces
   do j = 2,nlevsoi-1
      dzsoi(j)= 0.5_r8*(zsoi(j+1)-zsoi(j-1))
   enddo
   dzsoi(nlevsoi) = zsoi(nlevsoi)-zsoi(nlevsoi-1)

   zisoi(0) = 0._r8
   do j = 1, nlevsoi-1
      zisoi(j) = 0.5_r8*(zsoi(j)+zsoi(j+1))         !interface depths
   enddo
   zisoi(nlevsoi) = zsoi(nlevsoi) + 0.5_r8*dzsoi(nlevsoi)

   ! --------------------------------------------------------------------
   ! Initialize nitrogen deposition values 
   ! for now these are constants by gridcell, eventually they
   ! will be variables from the atmosphere, and at some point in between
   ! they will be specified time varying fields.
   ! --------------------------------------------------------------------

   ! Grid level initialization
   do g = begg, endg

      ! nitrogen deposition (forcing flux from atmosphere)
      ! convert rate from 1/yr -> 1/s
      
      forc_ndep(g) = ndep(g)/(86400._r8 * 365._r8)
      
   end do

   ! --------------------------------------------------------------------
   ! Initialize soil and lake levels
   ! Initialize soil color, thermal and hydraulic properties
   ! --------------------------------------------------------------------

   ! Column level initialization
!dir$ concurrent
!cdir nodep
   do c = begc, endc

      ! Set gridcell and landunit indices
      g = cgridcell(c)
      l = clandunit(c)

      ! Initialize restriction for min of soil potential (mm)
      smpmin(c) = -1.e8_r8

      ! Decay factor (m)
      hkdepth(c) = 1._r8/2.5_r8

      ! Maximum saturated fraction
      wtfact(c) = gti(g)

      ! Soil color
      isoicol(c) = soic2d(g)

      ! Soil hydraulic and thermal properties
      if (ltype(l)==istdlak .or. ltype(l)==istwet .or. &
          ltype(l)==istice .or. ltype(l)==isturb ) then
         do lev = 1,nlevsoi
            bsw(c,lev) = spval
            bsw2(c,lev) = spval
            psisat(c,lev) = spval
            vwcsat(c,lev) = spval
            watsat(c,lev) = spval
            hksat(c,lev) = spval
            sucsat(c,lev) = spval
            tkmg(c,lev) = spval
            tksatu(c,lev) = spval
            tkdry(c,lev) = spval
            csol(c,lev) = spval
            watdry(c,lev) = spval 
            watopt(c,lev) = spval 
         end do
      else
         do lev = 1,nlevsoi
            clay = clay3d(g,lev)
            sand = sand3d(g,lev)
            watsat(c,lev) = 0.489_r8 - 0.00126_r8*sand
            bd = (1._r8-watsat(c,lev))*2.7e3_r8
            xksat = 0.0070556_r8 *( 10._r8**(-0.884_r8+0.0153_r8*sand) ) ! mm/s
            tkm = (8.80_r8*sand+2.92_r8*clay)/(sand+clay)          ! W/(m K)

            bsw(c,lev) = 2.91_r8 + 0.159_r8*clay
            bsw2(c,lev) = -(3.10_r8 + 0.157_r8*clay - 0.003_r8*sand)
            psisat(c,lev) = -(exp((1.54_r8 - 0.0095_r8*sand + 0.0063_r8*(100.0_r8-sand-clay))*log(10.0_r8))*9.8e-5_r8)
            vwcsat(c,lev) = (50.5_r8 - 0.142_r8*sand - 0.037_r8*clay)/100.0_r8
            hksat(c,lev) = xksat
            sucsat(c,lev) = 10._r8 * ( 10._r8**(1.88_r8-0.0131_r8*sand) )
            tkmg(c,lev) = tkm ** (1._r8- watsat(c,lev))
            tksatu(c,lev) = tkmg(c,lev)*0.57_r8**watsat(c,lev)
            tkdry(c,lev) = (0.135_r8*bd + 64.7_r8) / (2.7e3_r8 - 0.947_r8*bd)
            csol(c,lev) = (2.128_r8*sand+2.385_r8*clay) / (sand+clay)*1.e6_r8  ! J/(m3 K)
            watdry(c,lev) = watsat(c,lev) * (316230._r8/sucsat(c,lev)) ** (-1._r8/bsw(c,lev)) 
            watopt(c,lev) = watsat(c,lev) * (158490._r8/sucsat(c,lev)) ** (-1._r8/bsw(c,lev)) 
         end do
      endif

      ! Define lake or non-lake levels layers
      if (ltype(l) == istdlak) then
         z(c,1:nlevlak) = zlak(1:nlevlak)
         dz(c,1:nlevlak) = dzlak(1:nlevlak)
      else
         z(c,1:nlevsoi) = zsoi(1:nlevsoi)
         zi(c,0:nlevsoi) = zisoi(0:nlevsoi)
         dz(c,1:nlevsoi) = dzsoi(1:nlevsoi)
      end if

      ! Initialize terms needed for dust model
      clay = clay3d(g,1)
      gwc_thr(c) = 0.17_r8 + 0.14_r8*clay*0.01_r8
      mss_frc_cly_vld(c) = min(clay*0.01_r8, 0.20_r8)

   end do

   ! pft level initialization
!dir$ concurrent
!cdir nodep
   do p = begp, endp

      ! Initialize maximum allowed dew

      dewmx(p)  = 0.1_r8

      ! Initialize root fraction (computing from surface, d is depth in meter):
      ! Y = 1 -1/2 (exp(-ad)+exp(-bd) under the constraint that
      ! Y(d =0.1m) = 1-beta^(10 cm) and Y(d=d_obs)=0.99 with
      ! beta & d_obs given in Zeng et al. (1998).

      c = pcolumn(p)
      if (ivt(p) /= noveg) then
         do lev = 1, nlevsoi-1
            rootfr(p,lev) = .5_r8*( exp(-roota_par(ivt(p)) * zi(c,lev-1))  &
                               + exp(-rootb_par(ivt(p)) * zi(c,lev-1))  &
                               - exp(-roota_par(ivt(p)) * zi(c,lev  ))  &
                               - exp(-rootb_par(ivt(p)) * zi(c,lev  )) )
         end do
         rootfr(p,nlevsoi) = .5_r8*( exp(-roota_par(ivt(p)) * zi(c,nlevsoi-1))  &
          
                                + exp(-rootb_par(ivt(p)) * zi(c,nlevsoi-1)) )

#if (defined CN)
         ! replacing the exponential rooting distribution
         ! with a linear decrease, going to zero at the bottom of the lowest
         ! soil layer for woody pfts, but going to zero at the bottom of
         ! layer 8 for non-woody pfts.  This corresponds to 3.43 m for woody
         ! bottom, vs 1.38 m for non-woody bottom.
         if (woody(ivt(p)) == 1) then
            bottom = nlevsoi
            slope = -2._r8/(zi(c,bottom)*zi(c,bottom))
            intercept   = 2._r8/zi(c,bottom)
            do lev = 1, bottom
               rootfr(p,lev) = dz(c,lev) * 0.5_r8 * ((intercept+slope*zi(c,lev-1)) + (intercept+slope*zi(c,lev)))
            end do
            if (bottom < nlevsoi) then
               do lev=bottom+1,nlevsoi
                  rootfr(p,lev) = 0._r8
               end do
            end if
         else
            bottom = 8
            slope = -2._r8/(zi(c,bottom)*zi(c,bottom))
            intercept   = 2._r8/zi(c,bottom)
            do lev=1,bottom
               rootfr(p,lev) = dz(c,lev) * 0.5_r8 * ((intercept+slope*zi(c,lev-1)) + (intercept+slope*zi(c,lev)))
            end do
            if (bottom < nlevsoi) then
               do lev=bottom+1,nlevsoi
                  rootfr(p,lev) = 0._r8
               end do
            end if
         end if
#endif
      else
         rootfr(p,1:nlevsoi) = 0._r8
      endif
      
      ! initialize rresis, for use in ecosystemdyn
      do lev = 1,nlevsoi
         rresis(p,lev) = 0._r8
      end do

   end do ! end pft level initialization
   
#if (defined CN)
   ! initialize the CN variables for special landunits, including lake points
   call CNiniSpecial()
#endif

   deallocate(soic2d,ndep,sand3d,clay3d,gti)

   if (masterproc) write (6,*) 'Successfully initialized time invariant variables'

end subroutine iniTimeConst
