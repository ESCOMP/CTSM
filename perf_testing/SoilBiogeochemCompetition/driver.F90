program SoilBiogeochemCompetition_driver

  !-----------------------------------------------------------------------
  ! Standalone timing harness for SoilBiogeochemCompetition.
  !
  ! Allocates synthetic column- and column-by-level arrays, fills them with
  ! deterministic non-trivial values, calls the routine niters times,
  ! computes a multi-output checksum, writes last_run.txt, and compares
  ! against baseline_checksum.txt when present.
  !
  ! Built-in timing (system_clock around the call loop, plus printed
  ! 'elapsed (s)' / 'per call (s)' lines) is gated by the cpp macro
  ! PERF_TIMING. Build with -DPERF_TIMING (the Makefile default) to keep
  ! it; build without (e.g. 'make TIMING=0') to leave timing to an
  ! external profiler.
  !
  ! Usage:
  !   ./driver [ncol [nlevdecomp [ndct [numfc [niters [use_nitrif_denitrif [carbon_only]]]]]]]
  !
  ! Defaults: 8000 10 8 -1 1 .true. .false.   (numfc=-1 means use ncol)
  !
  ! NOTE: outputs are zero-initialised once before the call loop. With
  ! niters > 1 the routine accumulates into some outputs (actual_immob,
  ! potential_immob in the non-nitrif branch), so the checksum after
  ! niters calls is not equal to niters * (checksum after 1 call).
  ! Canonical baseline uses niters=1.
  !-----------------------------------------------------------------------

  use SoilBiogeochemCompetition_mod, only : r8, SoilBiogeochemCompetition

  implicit none

  ! Run parameters (override-able via command-line args)
  integer :: ncol           = 8000
  integer :: nlevdecomp     = 10
  integer :: ndct           = 8       ! ndecomp_cascade_transitions
  integer :: numfc          = -1      ! -1 sentinel -> default to ncol
  integer :: niters         = 1
  logical :: use_nitrif_denitrif = .true.
  logical :: carbon_only         = .false.

  ! Fixed-for-now config (could be promoted to args later)
  integer , parameter :: decomp_method  = 2     ! mimics_decomp = 2
  integer , parameter :: mimics_decomp  = 2
  integer , parameter :: i_cop_mic      = 3
  integer , parameter :: i_oli_mic      = 4
  real(r8), parameter :: dt   = 1800.0_r8       ! 30-min decomp timestep
  real(r8), parameter :: bdnr = 1.0e-4_r8       ! bulk denit rate (1/s scaled)

  ! Competition coefficients (would come from params_inst in CTSM)
  real(r8), parameter :: compet_plant_no3  = 1.0_r8
  real(r8), parameter :: compet_plant_nh4  = 0.8_r8
  real(r8), parameter :: compet_decomp_no3 = 0.5_r8
  real(r8), parameter :: compet_decomp_nh4 = 1.0_r8
  real(r8), parameter :: compet_denit      = 0.3_r8
  real(r8), parameter :: compet_nit        = 0.5_r8

  integer :: begc, endc

  ! Pointer arrays (must be pointer to match the routine's signature)
  integer , pointer :: cascade_receiver_pool(:) => null()
  integer , pointer :: landunit(:)              => null()
  real(r8), pointer :: fpg(:)                   => null()
  real(r8), pointer :: fpi(:)                   => null()
  real(r8), pointer :: fpi_vr(:,:)              => null()
  real(r8), pointer :: nfixation_prof(:,:)      => null()
  real(r8), pointer :: plant_ndemand(:)         => null()
  real(r8), pointer :: sminn_vr(:,:)            => null()
  real(r8), pointer :: smin_nh4_vr(:,:)         => null()
  real(r8), pointer :: smin_no3_vr(:,:)         => null()
  real(r8), pointer :: c_overflow_vr(:,:,:)     => null()
  real(r8), pointer :: pot_f_nit_vr(:,:)        => null()
  real(r8), pointer :: pot_f_denit_vr(:,:)      => null()
  real(r8), pointer :: f_nit_vr(:,:)            => null()
  real(r8), pointer :: f_denit_vr(:,:)          => null()
  real(r8), pointer :: potential_immob(:)       => null()
  real(r8), pointer :: actual_immob(:)          => null()
  real(r8), pointer :: sminn_to_plant(:)        => null()
  real(r8), pointer :: sminn_to_denit_excess_vr(:,:)  => null()
  real(r8), pointer :: actual_immob_no3_vr(:,:) => null()
  real(r8), pointer :: actual_immob_nh4_vr(:,:) => null()
  real(r8), pointer :: smin_no3_to_plant_vr(:,:)=> null()
  real(r8), pointer :: smin_nh4_to_plant_vr(:,:)=> null()
  real(r8), pointer :: n2_n2o_ratio_denit_vr(:,:) => null()
  real(r8), pointer :: f_n2o_denit_vr(:,:)      => null()
  real(r8), pointer :: f_n2o_nit_vr(:,:)        => null()
  real(r8), pointer :: supplement_to_sminn_vr(:,:) => null()
  real(r8), pointer :: sminn_to_plant_vr(:,:)   => null()
  real(r8), pointer :: potential_immob_vr(:,:)  => null()
  real(r8), pointer :: actual_immob_vr(:,:)     => null()

  ! Allocatable / non-pointer args
  integer , allocatable :: filter_bgc_soilc(:)
  real(r8), allocatable :: dzsoi_decomp(:)
  real(r8), allocatable :: pmnf_decomp_cascade(:,:,:)
  real(r8), allocatable :: p_decomp_cn_gain(:,:,:)

  integer  :: iter
  real(r8) :: checksum
#ifdef PERF_TIMING
  integer(kind=8) :: t_start, t_end, t_rate
  real(r8) :: elapsed_s, per_call_s
#endif

  call parse_args()
  if (numfc < 0) numfc = ncol
  begc = 1
  endc = ncol

  call allocate_arrays()
  call fill_inputs()

#ifdef PERF_TIMING
  call system_clock(count_rate=t_rate)
  call system_clock(t_start)
#endif
  do iter = 1, niters
     call SoilBiogeochemCompetition( &
          begc, endc, nlevdecomp, ndct, &
          numfc, filter_bgc_soilc, &
          dt, bdnr, &
          use_nitrif_denitrif, carbon_only, &
          decomp_method, mimics_decomp, i_cop_mic, i_oli_mic, &
          compet_plant_no3, compet_plant_nh4, &
          compet_decomp_no3, compet_decomp_nh4, &
          compet_denit, compet_nit, &
          dzsoi_decomp, cascade_receiver_pool, landunit, &
          fpg, fpi, fpi_vr, nfixation_prof, plant_ndemand, &
          sminn_vr, smin_nh4_vr, smin_no3_vr, &
          c_overflow_vr, &
          pot_f_nit_vr, pot_f_denit_vr, f_nit_vr, f_denit_vr, &
          potential_immob, actual_immob, sminn_to_plant, &
          sminn_to_denit_excess_vr, &
          actual_immob_no3_vr, actual_immob_nh4_vr, &
          smin_no3_to_plant_vr, smin_nh4_to_plant_vr, &
          n2_n2o_ratio_denit_vr, f_n2o_denit_vr, f_n2o_nit_vr, &
          supplement_to_sminn_vr, sminn_to_plant_vr, &
          potential_immob_vr, actual_immob_vr, &
          pmnf_decomp_cascade, p_decomp_cn_gain)
  end do
#ifdef PERF_TIMING
  call system_clock(t_end)
  elapsed_s  = real(t_end - t_start, r8) / real(t_rate, r8)
  per_call_s = elapsed_s / real(niters, r8)
#endif

  call compute_checksum(checksum)
  call report(checksum)
  call write_last_run(checksum)
  call compare_to_baseline(checksum)

contains

  !---------------------------------------------------------------------
  subroutine parse_args()
    integer :: nargs
    character(len=64) :: arg

    nargs = command_argument_count()
    if (nargs >= 1) then; call get_command_argument(1, arg); read(arg,*) ncol;        end if
    if (nargs >= 2) then; call get_command_argument(2, arg); read(arg,*) nlevdecomp;  end if
    if (nargs >= 3) then; call get_command_argument(3, arg); read(arg,*) ndct;        end if
    if (nargs >= 4) then; call get_command_argument(4, arg); read(arg,*) numfc;       end if
    if (nargs >= 5) then; call get_command_argument(5, arg); read(arg,*) niters;      end if
    if (nargs >= 6) then; call get_command_argument(6, arg); read(arg,*) use_nitrif_denitrif; end if
    if (nargs >= 7) then; call get_command_argument(7, arg); read(arg,*) carbon_only;       end if
  end subroutine parse_args

  !---------------------------------------------------------------------
  subroutine allocate_arrays()
    allocate(filter_bgc_soilc(numfc))
    allocate(dzsoi_decomp(nlevdecomp))
    allocate(cascade_receiver_pool(ndct))
    allocate(landunit(begc:endc))

    allocate(fpg(begc:endc))
    allocate(fpi(begc:endc))
    allocate(plant_ndemand(begc:endc))
    allocate(potential_immob(begc:endc))
    allocate(actual_immob(begc:endc))
    allocate(sminn_to_plant(begc:endc))

    allocate(fpi_vr(begc:endc, 1:nlevdecomp))
    allocate(nfixation_prof(begc:endc, 1:nlevdecomp))
    allocate(sminn_vr(begc:endc, 1:nlevdecomp))
    allocate(smin_nh4_vr(begc:endc, 1:nlevdecomp))
    allocate(smin_no3_vr(begc:endc, 1:nlevdecomp))
    allocate(pot_f_nit_vr(begc:endc, 1:nlevdecomp))
    allocate(pot_f_denit_vr(begc:endc, 1:nlevdecomp))
    allocate(f_nit_vr(begc:endc, 1:nlevdecomp))
    allocate(f_denit_vr(begc:endc, 1:nlevdecomp))
    allocate(sminn_to_denit_excess_vr(begc:endc, 1:nlevdecomp))
    allocate(actual_immob_no3_vr(begc:endc, 1:nlevdecomp))
    allocate(actual_immob_nh4_vr(begc:endc, 1:nlevdecomp))
    allocate(smin_no3_to_plant_vr(begc:endc, 1:nlevdecomp))
    allocate(smin_nh4_to_plant_vr(begc:endc, 1:nlevdecomp))
    allocate(n2_n2o_ratio_denit_vr(begc:endc, 1:nlevdecomp))
    allocate(f_n2o_denit_vr(begc:endc, 1:nlevdecomp))
    allocate(f_n2o_nit_vr(begc:endc, 1:nlevdecomp))
    allocate(supplement_to_sminn_vr(begc:endc, 1:nlevdecomp))
    allocate(sminn_to_plant_vr(begc:endc, 1:nlevdecomp))
    allocate(potential_immob_vr(begc:endc, 1:nlevdecomp))
    allocate(actual_immob_vr(begc:endc, 1:nlevdecomp))

    allocate(c_overflow_vr(begc:endc, 1:nlevdecomp, 1:ndct))
    allocate(pmnf_decomp_cascade(begc:endc, 1:nlevdecomp, 1:ndct))
    allocate(p_decomp_cn_gain(begc:endc, 1:nlevdecomp, 1:ndct))
  end subroutine allocate_arrays

  !---------------------------------------------------------------------
  subroutine fill_inputs()
    integer  :: c, j, k, fc
    real(r8) :: frac, depthf

    ! Filter selects the first numfc columns
    do fc = 1, numfc
       filter_bgc_soilc(fc) = fc
    end do

    ! Per-layer thickness: 0.10, 0.15, ... — same shape as CTSM dzsoi_decomp
    do j = 1, nlevdecomp
       dzsoi_decomp(j) = 0.05_r8 + 0.05_r8 * real(j, r8)
    end do

    ! cascade_receiver_pool cycles through 1..5 so the if-test against
    ! i_cop_mic=3 and i_oli_mic=4 in the MIMICS branch hits sometimes.
    do k = 1, ndct
       cascade_receiver_pool(k) = mod(k - 1, 5) + 1
    end do

    ! landunit is read but never used downstream — value irrelevant.
    landunit(:) = 1

    ! Inputs that vary by column / level so different cells take
    ! different branches inside the routine.
    do c = begc, endc
       frac = real(c - begc + 1, r8) / real(endc - begc + 1, r8)   ! (0,1]
       plant_ndemand(c) = 1.0e-5_r8 * (1.0_r8 + 9.0_r8 * (1.0_r8 - frac))

       do j = 1, nlevdecomp
          depthf = real(j, r8) / real(nlevdecomp, r8)
          ! Range of sminn_vr: 0.05 .. ~2 g/m3 — exercises both
          ! 'supply > demand' and 'supply < demand' branches.
          sminn_vr(c, j)            = 0.05_r8 + 1.95_r8 * frac
          smin_nh4_vr(c, j)         = 0.4_r8 * sminn_vr(c, j)
          smin_no3_vr(c, j)         = 0.6_r8 * sminn_vr(c, j)
          nfixation_prof(c, j)      = 1.0_r8 / real(nlevdecomp, r8)
          potential_immob_vr(c, j)  = 5.0e-6_r8 * depthf
          pot_f_nit_vr(c, j)        = 1.0e-7_r8 * depthf
          pot_f_denit_vr(c, j)      = 2.0e-7_r8 * depthf
          n2_n2o_ratio_denit_vr(c, j) = 5.0_r8

          do k = 1, ndct
             pmnf_decomp_cascade(c, j, k) = 1.0e-6_r8 * frac * depthf * &
                                            real(k, r8) / real(ndct, r8)
             p_decomp_cn_gain(c, j, k)    = 8.0_r8 + real(mod(k, 4), r8)
          end do
       end do
    end do

    ! Zero all output / inout arrays. Routine accumulates into some of
    ! these (actual_immob, potential_immob in the non-nitrif branch),
    ! so they MUST start at zero for the first call.
    fpg(:)                       = 0.0_r8
    fpi(:)                       = 0.0_r8
    potential_immob(:)           = 0.0_r8
    actual_immob(:)              = 0.0_r8
    sminn_to_plant(:)            = 0.0_r8
    fpi_vr(:,:)                  = 0.0_r8
    f_nit_vr(:,:)                = 0.0_r8
    f_denit_vr(:,:)              = 0.0_r8
    sminn_to_denit_excess_vr(:,:)= 0.0_r8
    actual_immob_no3_vr(:,:)     = 0.0_r8
    actual_immob_nh4_vr(:,:)     = 0.0_r8
    smin_no3_to_plant_vr(:,:)    = 0.0_r8
    smin_nh4_to_plant_vr(:,:)    = 0.0_r8
    f_n2o_denit_vr(:,:)          = 0.0_r8
    f_n2o_nit_vr(:,:)            = 0.0_r8
    supplement_to_sminn_vr(:,:)  = 0.0_r8
    sminn_to_plant_vr(:,:)       = 0.0_r8
    actual_immob_vr(:,:)         = 0.0_r8
    c_overflow_vr(:,:,:)         = 0.0_r8
  end subroutine fill_inputs

  !---------------------------------------------------------------------
  subroutine compute_checksum(cs)
    real(r8), intent(out) :: cs
    cs = sum(sminn_to_plant) + sum(actual_immob) + sum(potential_immob) + &
         sum(fpg) + sum(fpi) + &
         sum(sminn_to_plant_vr) + sum(actual_immob_vr) + sum(fpi_vr) + &
         sum(f_nit_vr) + sum(f_denit_vr) + &
         sum(f_n2o_nit_vr) + sum(f_n2o_denit_vr) + &
         sum(sminn_to_denit_excess_vr) + sum(supplement_to_sminn_vr) + &
         sum(actual_immob_no3_vr) + sum(actual_immob_nh4_vr) + &
         sum(smin_no3_to_plant_vr) + sum(smin_nh4_to_plant_vr) + &
         sum(c_overflow_vr)
  end subroutine compute_checksum

  !---------------------------------------------------------------------
  subroutine report(checksum)
    real(r8), intent(in) :: checksum

    write(*,'(a)') '=== SoilBiogeochemCompetition standalone driver ==='
    write(*,'(a,i0)')      '  ncol                 = ', ncol
    write(*,'(a,i0)')      '  nlevdecomp           = ', nlevdecomp
    write(*,'(a,i0)')      '  ndct                 = ', ndct
    write(*,'(a,i0)')      '  numfc                = ', numfc
    write(*,'(a,i0)')      '  niters               = ', niters
    write(*,'(a,l1)')      '  use_nitrif_denitrif  = ', use_nitrif_denitrif
    write(*,'(a,l1)')      '  carbon_only          = ', carbon_only
    write(*,'(a,i0)')      '  decomp_method        = ', decomp_method
#ifdef PERF_TIMING
    write(*,'(a,es14.6)')  '  elapsed (s)          = ', elapsed_s
    write(*,'(a,es14.6)')  '  per call (s)         = ', per_call_s
#endif
    write(*,'(a,es24.16)') '  checksum             = ', checksum
  end subroutine report

  !---------------------------------------------------------------------
  subroutine write_last_run(checksum)
    real(r8), intent(in) :: checksum
    integer :: u
    open(newunit=u, file='last_run.txt', status='replace', action='write')
    write(u,'(a,i0)')      'ncol                 ', ncol
    write(u,'(a,i0)')      'nlevdecomp           ', nlevdecomp
    write(u,'(a,i0)')      'ndct                 ', ndct
    write(u,'(a,i0)')      'numfc                ', numfc
    write(u,'(a,i0)')      'niters               ', niters
    write(u,'(a,l1)')      'use_nitrif_denitrif  ', use_nitrif_denitrif
    write(u,'(a,l1)')      'carbon_only          ', carbon_only
    write(u,'(a,i0)')      'decomp_method        ', decomp_method
    write(u,'(a,es24.16)') 'checksum             ', checksum
    close(u)
  end subroutine write_last_run

  !---------------------------------------------------------------------
  subroutine compare_to_baseline(checksum)
    real(r8), intent(in) :: checksum
    integer  :: u, ios
    logical  :: exists
    character(len=64) :: key
    integer  :: b_ncol, b_nlevdecomp, b_ndct, b_numfc, b_niters, b_decomp_method
    logical  :: b_use_nitrif_denitrif, b_carbon_only
    real(r8) :: b_checksum, tol, diff

    inquire(file='baseline_checksum.txt', exist=exists)
    if (.not. exists) then
       write(*,'(a)') '  baseline             = (no baseline_checksum.txt found; skipping compare)'
       return
    end if

    open(newunit=u, file='baseline_checksum.txt', status='old', action='read', iostat=ios)
    if (ios /= 0) then
       write(*,'(a)') '  baseline             = (could not open baseline_checksum.txt)'
       return
    end if
    read(u,*,iostat=ios) key, b_ncol;                if (ios /= 0) goto 99
    read(u,*,iostat=ios) key, b_nlevdecomp;          if (ios /= 0) goto 99
    read(u,*,iostat=ios) key, b_ndct;                if (ios /= 0) goto 99
    read(u,*,iostat=ios) key, b_numfc;               if (ios /= 0) goto 99
    read(u,*,iostat=ios) key, b_niters;              if (ios /= 0) goto 99
    read(u,*,iostat=ios) key, b_use_nitrif_denitrif; if (ios /= 0) goto 99
    read(u,*,iostat=ios) key, b_carbon_only;         if (ios /= 0) goto 99
    read(u,*,iostat=ios) key, b_decomp_method;       if (ios /= 0) goto 99
    read(u,*,iostat=ios) key, b_checksum;            if (ios /= 0) goto 99
    close(u)

    if (b_ncol /= ncol .or. b_nlevdecomp /= nlevdecomp .or. &
        b_ndct /= ndct .or. b_numfc /= numfc .or. &
        b_niters /= niters .or. b_decomp_method /= decomp_method .or. &
        b_use_nitrif_denitrif .neqv. use_nitrif_denitrif .or. &
        b_carbon_only .neqv. carbon_only) then
       write(*,'(a)') '  baseline             = (param set differs; skipping compare)'
       return
    end if

    tol  = 1.0e-10_r8 * max(abs(b_checksum), 1.0_r8)
    diff = abs(checksum - b_checksum)
    if (diff <= tol) then
       write(*,'(a,es14.6,a)') '  baseline             = MATCH (|diff| = ', diff, ')'
    else
       write(*,'(a,es24.16)')  '  baseline value       = ', b_checksum
       write(*,'(a,es14.6,a,es14.6,a)') &
           '  baseline             = MISMATCH (|diff| = ', diff, ', tol = ', tol, ')'
    end if
    return

 99 continue
    close(u)
    write(*,'(a)') '  baseline             = (parse error in baseline_checksum.txt; skipping compare)'
  end subroutine compare_to_baseline

end program SoilBiogeochemCompetition_driver
