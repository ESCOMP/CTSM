#include <misc.h>
#include <preproc.h>

module pftvarcon

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: pftvarcon
!
! !DESCRIPTION:
! Module containing vegetation constants and method to
! eads and initialize vegetation (PFT) constants.
!
! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use abortutils  , only : endrun
  use clm_varpar
!
! !PUBLIC TYPES:
  implicit none
  save
!
! Vegetation type constants
!
  character(len=40) pftname(0:numpft) !PFT description

  integer ncorn                  !value for corn
  integer nwheat                 !value for wheat
  integer noveg                  !value for not vegetated 
  integer ntree                  !value for last type of tree

  real(r8):: dleaf(0:numpft)       !characteristic leaf dimension (m)
  real(r8):: c3psn(0:numpft)       !photosynthetic pathway: 0. = c4, 1. = c3
  real(r8):: vcmx25(0:numpft)      !max rate of carboxylation at 25C (umol CO2/m**2/s)
  real(r8):: mp(0:numpft)          !slope of conductance-to-photosynthesis relationship
  real(r8):: qe25(0:numpft)        !quantum efficiency at 25C (umol CO2 / umol photon)
  real(r8):: xl(0:numpft)          !leaf/stem orientation index
  real(r8):: rhol(0:numpft,numrad) !leaf reflectance: 1=vis, 2=nir
  real(r8):: rhos(0:numpft,numrad) !stem reflectance: 1=vis, 2=nir
  real(r8):: taul(0:numpft,numrad) !leaf transmittance: 1=vis, 2=nir
  real(r8):: taus(0:numpft,numrad) !stem transmittance: 1=vis, 2=nir
  real(r8):: z0mr(0:numpft)        !ratio of momentum roughness length to canopy top height (-)
  real(r8):: displar(0:numpft)     !ratio of displacement height to canopy top height (-)
  real(r8):: roota_par(0:numpft)   !CLM rooting distribution parameter [1/m]
  real(r8):: rootb_par(0:numpft)   !CLM rooting distribution parameter [1/m]
  real(r8):: crop(0:numpft)        ! crop pft: 0. = not crop, 1. = crop pft
  real(r8):: sla(0:numpft)         ! specific leaf area [m2 leaf g-1 carbon]
  real(r8):: smpso(0:numpft)       !soil water potential at full stomatal opening (mm)
  real(r8):: smpsc(0:numpft)       !soil water potential at full stomatal closure (mm)
  real(r8):: fnitr(0:numpft)       !foliage nitrogen limitation factor (-)
  ! begin new pft parameters for CN code
  real(r8):: slatop(0:numpft)      !SLA at top of canopy [m^2/gC]
  real(r8):: dsladlai(0:numpft)    !dSLA/dLAI [m^2/gC]
  real(r8):: leafcn(0:numpft)      !leaf C:N [gC/gN]
  real(r8):: flnr(0:numpft)        !fraction of leaf N in Rubisco [no units]
  real(r8):: woody(0:numpft)       !woody lifeform flag (0 or 1)
  real(r8):: lflitcn(0:numpft)      !leaf litter C:N (gC/gN)
  real(r8):: frootcn(0:numpft)      !fine root C:N (gC/gN)
  real(r8):: livewdcn(0:numpft)     !live wood (phloem and ray parenchyma) C:N (gC/gN)
  real(r8):: deadwdcn(0:numpft)     !dead wood (xylem and heartwood) C:N (gC/gN)
  real(r8):: froot_leaf(0:numpft)   !allocation parameter: new fine root C per new leaf C (gC/gC) 
  real(r8):: stem_leaf(0:numpft)    !allocation parameter: new stem c per new leaf C (gC/gC)
  real(r8):: croot_stem(0:numpft)   !allocation parameter: new coarse root C per new stem C (gC/gC)
  real(r8):: flivewd(0:numpft)      !allocation parameter: fraction of new wood that is live (phloem and ray parenchyma) (no units)
  real(r8):: fcur(0:numpft)         !allocation parameter: fraction of allocation that goes to currently displayed growth, remainder to storage
  real(r8):: lf_flab(0:numpft)      !leaf litter labile fraction
  real(r8):: lf_fcel(0:numpft)      !leaf litter cellulose fraction
  real(r8):: lf_flig(0:numpft)      !leaf litter lignin fraction
  real(r8):: fr_flab(0:numpft)      !fine root litter labile fraction
  real(r8):: fr_fcel(0:numpft)      !fine root litter cellulose fraction
  real(r8):: fr_flig(0:numpft)      !fine root litter lignin fraction
  real(r8):: dw_fcel(0:numpft)      !dead wood cellulose fraction
  real(r8):: dw_flig(0:numpft)      !dead wood lignin fraction
  real(r8):: leaf_long(0:numpft)    !leaf longevity (yrs)
  real(r8):: evergreen(0:numpft)    !binary flag for evergreen leaf habit (0 or 1)
  real(r8):: stress_decid(0:numpft) !binary flag for stress-deciduous leaf habit (0 or 1)
  real(r8):: season_decid(0:numpft) !binary flag for seasonal-deciduous leaf habit (0 or 1)
  ! begin new pft parameters for CN-fire code
  real(r8):: resist(0:numpft)       !resistance to fire (no units)
  ! end new pft parameters for CN code
  ! begin pft parameters for DGVM code
  real(r8) pftpar(0:numpft,1:npftpar) ! the rest for use with DGVM
  real(r8) lm_sapl(0:numpft)
  real(r8) sm_sapl(0:numpft)
  real(r8) hm_sapl(0:numpft)
  real(r8) rm_sapl(0:numpft)
  logical  tree(0:numpft)
  logical  summergreen(0:numpft)
  logical  raingreen(0:numpft)
  ! end pft parameters for DGVM code

  real(r8), parameter :: reinickerp = 1.6_r8 !parameter in allometric equation
  real(r8), parameter :: wooddens = 2.0e5_r8 !wood density (gC/m3)
  real(r8), parameter :: latosa = 8.0e3_r8   !ratio of leaf area to sapwood cross-sectional
                                          !area (Shinozaki et al 1964a,b)
  real(r8), parameter :: allom1 = 100.0_r8   !parameters in allometric
  real(r8), parameter :: allom2 =  40.0_r8
  real(r8), parameter :: allom3 =   0.5_r8
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: pftconrd ! Read and initialize vegetation (PFT) constants
!
! !REVISION HISTORY:
! Created by Sam Levis (put into module form by Mariana Vertenstein)
! 10/21/03, Peter Thornton: Added new variables for CN code
!
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: pftconrd
!
! !INTERFACE:
  subroutine pftconrd
!
! !DESCRIPTION:
! Read and initialize vegetation (PFT) constants
!
! !USES:
    use fileutils , only : opnfil, getfil, relavu, getavu
    use clm_varctl, only : fpftcon
    use spmdMod   , only : masterproc, mpicom, MPI_REAL8
!
! !ARGUMENTS:
    implicit none
!
! !CALLED FROM:
! routine initialize in module initializeMod
!
! !REVISION HISTORY:
! Created by Gordon Bonan
!
!EOP
!
! !LOCAL VARIABLES:
    character(len=256) :: locfn ! local file name
    integer :: i,n              ! loop indices
    integer :: ier              ! error code
!-----------------------------------------------------------------------

    ! Set specific vegetation type values

    ncorn  = 15
    nwheat = 16
    ntree  =  8  ! value for last type of tree
    noveg  =  0  ! value for non-vegetated

    ! Assign unit number to file. Get local file.
    ! Open file and read PFT's.
    ! Close and release file.

    if (masterproc) then
       write (6,*) 'Attempting to read PFT physiological data .....'
       n = getavu()
       call getfil (fpftcon, locfn, 0)
       call opnfil (locfn, n, 'f')
       do i = 1, numpft
#if (defined CN)
          read (n,*,iostat=ier)  pftname(i),              &
               z0mr(i)   , displar(i), dleaf(i)  , c3psn(i)  , &
               vcmx25(i) , mp(i)     , qe25(i)   , rhol(i,1) , &
               rhol(i,2) , rhos(i,1) , rhos(i,2) , taul(i,1) , &
               taul(i,2) , taus(i,1) , taus(i,2) , xl(i)     , &
               roota_par(i), rootb_par(i), slatop(i), dsladlai(i), &
					leafcn(i), flnr(i), &
                                        smpso(i), smpsc(i), fnitr(i), &
               woody(i), lflitcn(i), frootcn(i), livewdcn(i), &
               deadwdcn(i), froot_leaf(i), stem_leaf(i), croot_stem(i), &
               flivewd(i), fcur(i), lf_flab(i), lf_fcel(i), lf_flig(i), &
               fr_flab(i), fr_fcel(i), fr_flig(i), dw_fcel(i), dw_flig(i), &
               leaf_long(i), evergreen(i), stress_decid(i), season_decid(i), &
               resist(i)
#else
          read (n,*,iostat=ier)  pftname(i), &
               z0mr(i)   , displar(i), dleaf(i)  , c3psn(i)  , &
               vcmx25(i) , mp(i)     , qe25(i)   , rhol(i,1) , &
               rhol(i,2) , rhos(i,1) , rhos(i,2) , taul(i,1) , &
               taul(i,2) , taus(i,1) , taus(i,2) , xl(i)     , &
               roota_par(i), rootb_par(i), slatop(i), dsladlai(i), &
					leafcn(i), flnr(i), &
                                        smpso(i), smpsc(i), fnitr(i)
#endif 
          if (ier /= 0) then
             write(6,*)'pftconrd: error in reading in pft data'
             call endrun()
          end if
       end do
       call relavu (n)
    end if

    ! Set crop explicitly here (in future will be on pft dataset)

    if (masterproc) then
       crop(:)  = 0
       crop(15) = 1
       crop(16) = 1
    end if

    ! Define PFT zero to be bare ground

    pftname(noveg) = 'not_vegetated'
    z0mr(noveg) = 0._r8
    displar(noveg) = 0._r8
    dleaf(noveg) = 0._r8
    c3psn(noveg) = 1._r8
    vcmx25(noveg) = 0._r8
    mp(noveg) = 9._r8
    qe25(noveg) = 0._r8
    rhol(noveg,1) = 0._r8
    rhol(noveg,2) = 0._r8
    rhos(noveg,1) = 0._r8
    rhos(noveg,2) = 0._r8
    taul(noveg,1) = 0._r8
    taul(noveg,2) = 0._r8
    taus(noveg,1) = 0._r8
    taus(noveg,2) = 0._r8
    xl(noveg) = 0._r8
    roota_par(noveg) = 0._r8
    rootb_par(noveg) = 0._r8
    crop(noveg) = 0._r8
    smpso(noveg) = 0._r8
    smpsc(noveg) = 0._r8
    fnitr(noveg) = 0._r8
    slatop(noveg) = 0._r8
    dsladlai(noveg) = 0._r8
    leafcn(noveg) = 1._r8
    flnr(noveg) = 0._r8
    ! begin variables used only for CN code
    woody(noveg) = 0._r8
    lflitcn(noveg) = 1._r8
    frootcn(noveg) = 1._r8
    livewdcn(noveg) = 1._r8
    deadwdcn(noveg) = 1._r8
    froot_leaf(noveg) = 0._r8
    stem_leaf(noveg) = 0._r8
    croot_stem(noveg) = 0._r8
    flivewd(noveg) = 0._r8
    fcur(noveg) = 0._r8
    lf_flab(noveg) = 0._r8
    lf_fcel(noveg) = 0._r8
    lf_flig(noveg) = 0._r8
    fr_flab(noveg) = 0._r8
    fr_fcel(noveg) = 0._r8
    fr_flig(noveg) = 0._r8
    dw_fcel(noveg) = 0._r8
    dw_flig(noveg) = 0._r8
    leaf_long(noveg) = 0._r8
    evergreen(noveg) = 0._r8
    stress_decid(noveg) = 0._r8
    season_decid(noveg) = 0._r8
    resist(noveg) = 1._r8
    ! end variables used only for CN code

    call mpi_bcast (z0mr, size(z0mr), MPI_REAL8, 0, mpicom, ier)
    call mpi_bcast (displar, size(displar), MPI_REAL8, 0, mpicom, ier)
    call mpi_bcast (dleaf, size(dleaf), MPI_REAL8, 0, mpicom, ier)
    call mpi_bcast (c3psn, size(c3psn), MPI_REAL8, 0, mpicom, ier)
    call mpi_bcast (vcmx25, size(vcmx25), MPI_REAL8, 0, mpicom, ier)
    call mpi_bcast (mp, size(mp), MPI_REAL8, 0, mpicom, ier)
    call mpi_bcast (qe25, size(qe25), MPI_REAL8, 0, mpicom, ier)
    call mpi_bcast (rhol, size(rhol), MPI_REAL8, 0, mpicom, ier)
    call mpi_bcast (rhos, size(rhos), MPI_REAL8, 0, mpicom, ier)
    call mpi_bcast (taul, size(taul), MPI_REAL8, 0, mpicom, ier)
    call mpi_bcast (taus, size(taus), MPI_REAL8, 0, mpicom, ier)
    call mpi_bcast (xl, size(xl), MPI_REAL8, 0, mpicom, ier)
    call mpi_bcast (roota_par, size(roota_par), MPI_REAL8, 0, mpicom, ier)
    call mpi_bcast (rootb_par, size(rootb_par), MPI_REAL8, 0, mpicom, ier)
    call mpi_bcast (crop, size(crop), MPI_REAL8, 0, mpicom, ier)
    call mpi_bcast (smpso, size(smpso), MPI_REAL8, 0, mpicom, ier)
    call mpi_bcast (smpsc, size(smpsc), MPI_REAL8, 0, mpicom, ier)
    call mpi_bcast (fnitr, size(fnitr), MPI_REAL8, 0, mpicom, ier)

    ! begin variables used only for CN code
    call mpi_bcast (slatop, size(slatop), MPI_REAL8, 0, mpicom, ier)
    call mpi_bcast (dsladlai, size(dsladlai), MPI_REAL8, 0, mpicom, ier)
    call mpi_bcast (leafcn, size(leafcn), MPI_REAL8, 0, mpicom, ier)
    call mpi_bcast (flnr, size(flnr), MPI_REAL8, 0, mpicom, ier)
    call mpi_bcast (woody, size(woody), MPI_REAL8, 0, mpicom, ier)
    call mpi_bcast (lflitcn, size(lflitcn), MPI_REAL8, 0, mpicom, ier)
    call mpi_bcast (frootcn, size(frootcn), MPI_REAL8, 0, mpicom, ier)
    call mpi_bcast (livewdcn, size(livewdcn), MPI_REAL8, 0, mpicom, ier)
    call mpi_bcast (deadwdcn, size(deadwdcn), MPI_REAL8, 0, mpicom, ier)
    call mpi_bcast (froot_leaf, size(froot_leaf), MPI_REAL8, 0, mpicom, ier)
    call mpi_bcast (stem_leaf, size(stem_leaf), MPI_REAL8, 0, mpicom, ier)
    call mpi_bcast (croot_stem, size(croot_stem), MPI_REAL8, 0, mpicom, ier)
    call mpi_bcast (flivewd, size(flivewd), MPI_REAL8, 0, mpicom, ier)
    call mpi_bcast (fcur, size(fcur), MPI_REAL8, 0, mpicom, ier)
    call mpi_bcast (lf_flab, size(lf_flab), MPI_REAL8, 0, mpicom, ier)
    call mpi_bcast (lf_fcel, size(lf_fcel), MPI_REAL8, 0, mpicom, ier)
    call mpi_bcast (lf_flig, size(lf_flig), MPI_REAL8, 0, mpicom, ier)
    call mpi_bcast (fr_flab, size(fr_flab), MPI_REAL8, 0, mpicom, ier)
    call mpi_bcast (fr_fcel, size(fr_fcel), MPI_REAL8, 0, mpicom, ier)
    call mpi_bcast (fr_flig, size(fr_flig), MPI_REAL8, 0, mpicom, ier)
    call mpi_bcast (dw_fcel, size(dw_fcel), MPI_REAL8, 0, mpicom, ier)
    call mpi_bcast (dw_flig, size(dw_flig), MPI_REAL8, 0, mpicom, ier)
    call mpi_bcast (leaf_long, size(leaf_long), MPI_REAL8, 0, mpicom, ier)
    call mpi_bcast (evergreen, size(evergreen), MPI_REAL8, 0, mpicom, ier)
    call mpi_bcast (stress_decid, size(stress_decid), MPI_REAL8, 0, mpicom, ier)
    call mpi_bcast (season_decid, size(season_decid), MPI_REAL8, 0, mpicom, ier)
    call mpi_bcast (resist, size(resist), MPI_REAL8, 0, mpicom, ier)
    ! end variables used only for CN code

    if (masterproc) then
       write (6,*) 'Successfully read PFT physiological data'
       write (6,*)
    end if

  end subroutine pftconrd

end module pftvarcon

