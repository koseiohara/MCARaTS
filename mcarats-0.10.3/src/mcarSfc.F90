!-*-f90-*-
 
! @LICENSE_HEADER_START@
!
!   This file is part of MCARaTS.
!   
!   --
!   MCARaTS: Monte Carlo Atmospheric Radiative Transfer Simulator
!   
!   Copyright (C) 2006-2016 Hironobu Iwabuchi.
!   
!   MCARaTS is free software: you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation, either version 3 of the License, or
!   (at your option) any later version.
!   
!   MCARaTS is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!   GNU General Public License for more details.
!   
!   You should have received a copy of the GNU General Public License
!   along with MCARaTS. If not, see <http://www.gnu.org/licenses/>.
!
! @LICENSE_HEADER_END@
 
 


!+
! Module for surfaces
!-
module mcarSfc 

  use globals
#if UseMPI == 1
  use globals_MPI
#endif
  use hparx
  use mcarUtl
  use mcarDsm
  use mcarRpv
  use mcarLsrt
  implicit none
  private

  ! Public
  public :: mcarSfc__user_init
  public :: mcarSfc__init_set
  public :: mcarSfc__init
  public :: mcarSfc__tune_jobs
  public :: mcarSfc__user_job
  public :: mcarSfc__prep_job
  public :: mcarSfc__final
  public :: mcarSfc__angDistr_R
  public :: mcarSfc__bsAlbedo_R
  public :: mcarSfc__emittance_R
  public :: mcarSfc__fMC_emit_1R
  public :: mcarSfc__intpOptProp
  public :: mcarSfc__prep_temp
  public :: mcarSfc__prep_optProp
  public :: mcarSfc__newDirec

  ! (Partly-public, read-only) User variables in namelist (1) for initialization
  integer, parameter :: Sfc_NPAR = 5
  character(256), save :: Sfc_inpfile = ' ' ! file name for input of surface properties
  integer,  save :: Sfc_mfmt = 1 ! flag for input file format (0=text, 1=binary)
  integer,  save :: Sfc_mbrdf(4) = (/1, 1, 1, 1/) ! flags of on/off status for BRDF models
  integer,  save :: Sfc_nxb = 1  ! # of X grid points
  integer,  save :: Sfc_nyb = 1  ! # of Y grid points
  integer,  save :: Sfc_nsco = 60 ! # of uz0 for coefficient table
  integer,  save :: Sfc_nsuz = 200 ! # of uz0 grids for albedo LUTs
  namelist /mcarSfc_nml_init/Sfc_inpfile, Sfc_mfmt, Sfc_mbrdf, Sfc_nxb, Sfc_nyb, Sfc_nsco, Sfc_nsuz

  ! (Partly-public, read-only) User variables in namelist (2) for core calculations
  integer,  save :: Sfc_idread = 0  ! index of data to be read (0: no data to be read)
  real(R_), save :: Sfc_tmp = 300.0_R_ ! temperature
  integer,  save :: Sfc_mtype = 1   ! surface BRDF type
  real(R_), save :: Sfc_param(Sfc_NPAR) = (/1.0_R_, 0.0_R_, 0.0_R_, 0.0_R_, 0.0_R_/) ! BRDF parameters
  !integer,  save :: Sfc_mtype  = 2 ! surface BRDF type (2:DSM)
  !real(R_), save :: Sfc_param(Sfc_NPAR) = (/0.2_R_, 0.1_R_, 1.33_R_, 2.0e-9_R_, 0.04_R_/)
  integer,  save :: Sfc_nudsm  = 14 ! # of table grid points for DSM model
  integer,  save :: Sfc_nurpv  = 8  ! # of table grid points for RPV model
  integer,  save :: Sfc_nulsrt = 14 ! # of table grid points for LSRT model
  integer,  save :: Sfc_nqpot  = 24 ! # of quadrature points for preprocess
  real(R_), save :: Sfc_rrmax  = 5.0_R_ ! max factor for relative BRDF used for random directions
  real(R_), save :: Sfc_rrexp  = 0.5_R_ ! scaling exponent for relative BRDF used for random directions
  namelist /mcarSfc_nml_job/Sfc_idread, Sfc_tmp, Sfc_mtype, Sfc_param, Sfc_nqpot, Sfc_rrmax, Sfc_rrexp, &
       & Sfc_nudsm, Sfc_nurpv, Sfc_nulsrt

  ! (Public, read-only) Surface model parameters
  real(R_), save :: Sfc_uz0min = 0.01_R_
  real(R_), save :: Sfc_r2min = 0.001_R_
  real(R_), save :: Sfc_facx, Sfc_facy
  integer,  save, allocatable :: Sfc_jsfc2d(:,:)   !(nxb,nyb)
  real(R_), save, allocatable :: Sfc_tmps2d(:,:)   !(nxb,nyb)
  real(R_), save, allocatable :: Sfc_psfc2d(:,:,:) !(nxb,nyb,Sfc_NPAR)
  real(R_), save, allocatable :: Sfc_scolut(:,:,:) !(nsco,nxb,nyb)
  real(R_), save, allocatable :: Sfc_salbs(:,:,:)  !(2,nxb,nyb)
  real(R_), save, allocatable :: Sfc_salut1(:)     !(nsuz)
  real(R_), save, allocatable :: Sfc_salut2(:)     !(nsuz)
  public :: Sfc_nxb, Sfc_nyb
  public :: Sfc_nudsm, Sfc_nurpv, Sfc_nulsrt
  public :: Sfc_uz0min, Sfc_r2min, Sfc_facx, Sfc_facy
  public :: Sfc_tmps2d, Sfc_psfc2d, Sfc_scolut, Sfc_salbs, Sfc_jsfc2d
  public :: Sfc_salut1, Sfc_salut2

  ! Private variables
  integer,  save :: Sfc_iud = -1     ! unit index of input data file (-1: not opened)
  integer,  save :: Sfc_iroot = 0    ! index of root PE
  integer,  save :: Sfc_npsfc(4) = (/1, 5, 4, 3/) ! # of BRDF parameters

  ! Common usage of parameter packet, pac(1:6)
  !   isfc=1 : Lambertian
  !     pac(1) : albedo
  !     pac(6) : albedo
  !   isfc=2 : DSM (diffuse-specular mixture)
  !     pac(1) : diffuse albedo
  !     pac(2) : fraction of diffuse reflection
  !     pac(3) : real part of refractive index
  !     pac(4) : imaginary part of refractive index
  !     pac(5) : variance of facet slope
  !     pac(6) : albedo
  !   isfc=3 : RPV (Rahmen-Pinty-Verstraete)
  !     pac(1) : "rho", the first parameter
  !     pac(2) : "k", exponent parameter
  !     pac(3) : "THETA", asymmetry parameter
  !     pac(4) : "delta", sharpness of the hot spot
  !     pac(5) : function A
  !     pac(6) : albedo
  !   isfc=4 : LSRT (Li-Sparse-Ross-Thick)
  !     pac(1) : "k_L", weight for Lambertian reflection
  !     pac(2) : "k_g", weight for geometrical optics
  !     pac(3) : "k_v", weight for volume scattering
  !     pac(6) : albedo
  !// Some elements are not used in some procedures.

contains

  !+
  ! Read in namelist variables
  !-
  subroutine mcarSfc__user_init(iu, filepath) 

    integer, intent(in) :: iu   ! file unit index for namelist input
    character(*), intent(in) :: filepath ! path name for I/O files
    integer :: ios
    read (iu, nml=mcarSfc_nml_init, iostat=ios)
    call err_read(ios, iu, 'mcarSfc__user_init: Invalid namelist input for mcarSfc_nml_init.')
    if (len_trim(Sfc_inpfile) > 0) Sfc_inpfile = trim(filepath)//'/'//Sfc_inpfile

  end subroutine mcarSfc__user_init


  !+
  ! Set parameters for initialization
  !-
  subroutine mcarSfc__init_set(mbrdf, nxb, nyb, nsco) 

    integer,  intent(in) :: mbrdf(:), nxb, nyb, nsco
    Sfc_mbrdf(1:4) = mbrdf(1:4)
    Sfc_nxb = nxb
    Sfc_nyb = nyb
    Sfc_nsco = nsco
    Sfc_nsuz = 200

  end subroutine mcarSfc__init_set


  !+
  ! Initialize this module
  !-
  subroutine mcarSfc__init(irank, iroot) 

    integer, intent(in) :: irank, iroot
    integer :: nrec

    ! Allocate
    if (allocated(Sfc_jsfc2d)) call mcarSfc__final()
    allocate (Sfc_jsfc2d(Sfc_nxb, Sfc_nyb))
    allocate (Sfc_tmps2d(Sfc_nxb, Sfc_nyb))
    allocate (Sfc_psfc2d(Sfc_nxb, Sfc_nyb, Sfc_NPAR))
    allocate (Sfc_scolut(Sfc_nsco, Sfc_nxb, Sfc_nyb))
    allocate (Sfc_salbs(2, Sfc_nxb, Sfc_nyb))
    allocate (Sfc_salut1(Sfc_nsuz))
    allocate (Sfc_salut2(Sfc_nsuz))

    ! Open input data file
    Sfc_iroot = iroot
    if (len_trim(Sfc_inpfile) > 0 .and. irank == iroot) then ! root PE only
       Sfc_iud = freeUnit_I(10)
       if (Sfc_mfmt == 0) then
          call open_seq(Sfc_iud, Sfc_inpfile, 'old')
       else
          nrec = recordLen_I(1, Sfc_nxb * Sfc_nyb, 4) ! Apr. 27, 2012
          call open_dir(Sfc_iud, Sfc_inpfile, nrec, 'old')
       end if
    else
       Sfc_iud = -1
    end if

  end subroutine mcarSfc__init


  !+
  ! Finalize this module
  !-
  subroutine mcarSfc__final() 

    if (allocated(Sfc_jsfc2d)) then
       deallocate (Sfc_jsfc2d, Sfc_tmps2d, Sfc_psfc2d, Sfc_scolut, Sfc_salbs)
       deallocate (Sfc_salut1, Sfc_salut2)
       if (Sfc_iud >= 1) then ! root PE only
          close (Sfc_iud)
          Sfc_iud = -1
       end if
    end if

  end subroutine mcarSfc__final


  !+
  ! Tune a preset for technical parameters
  !-
  subroutine mcarSfc__tune_jobs(moptim) 

    integer, intent(in) :: moptim ! optimization flag [-2,3]

    if (moptim <= -2) then ! as defaults
       return
    else if (moptim == -1) then ! no optimization
       call mcarSfc__set_tech(64, 20.0_R_, 0.5_R_, 30, 20, 30)
    else if (moptim == 0) then ! unbiased optimization
       call mcarSfc__set_tech(64, 10.0_R_, 0.5_R_, 30, 20, 30)
    else if (moptim == 1) then ! convervative
       call mcarSfc__set_tech(40, 6.0_R_, 0.5_R_, 16, 9, 16)
    else if (moptim == 2) then ! standard
       call mcarsfc__set_tech(20, 4.0_R_, 0.5_R_, 12, 7, 12)
    else ! quick-and-dirty
       call mcarSfc__set_tech(10, 2.0_R_, 0.5_R_, 7, 4, 7)
    end if
    
  end subroutine mcarSfc__tune_jobs


  !+
  ! Set technical parameters
  !-
  subroutine mcarSfc__set_tech(nqpot, rrmax, rrexp, nudsm, nurpv, nulsrt)

    real(R_), intent(in) :: rrmax, rrexp
    integer,  intent(in) :: nqpot, nudsm, nurpv, nulsrt
    Sfc_nqpot  = nqpot
    Sfc_rrmax  = rrmax
    Sfc_rrexp  = rrexp
    Sfc_nudsm  = min(nudsm,  Sfc_nsco/2)
    Sfc_nurpv  = min(nurpv,  Sfc_nsco/3)
    Sfc_nulsrt = min(nulsrt, Sfc_nsco/2)

  end subroutine mcarSfc__set_tech


  !+
  ! Read in namelist variables
  !-
  subroutine mcarSfc__user_job(iu, idread0) 

    integer, intent(in) :: iu   ! file unit index for namelist input
    integer, intent(in) :: idread0 ! default value for Sfc_idread
    integer :: ios
    Sfc_idread = idread0
    read (iu, nml=mcarSfc_nml_job, iostat=ios)
    call err_read(ios, iu, 'mcarSfc__user_job: Invalid namelist input for mcarSfc_nml_job.')

  end subroutine mcarSfc__user_job


  !+
  ! Initialize this module
  !-
  subroutine mcarSfc__prep_job(irank, mbswap, xmax, ymax, tmpmin, tmpmax)

    integer,  intent(in) :: irank, mbswap
    real(R_), intent(in) :: xmax, ymax
    real(R_), intent(out) :: tmpmin, tmpmax
    integer :: i

    ! Initialize
    Sfc_facx = real(Sfc_nxb, R_) / xmax
    Sfc_facy = real(Sfc_nyb, R_) / ymax
    !// Note there is a possibility that ixb = int(x * facx) + 1 > nxb. Use ixb = min(nxb, ixb).
    Sfc_uz0min = 0.02_R_
    Sfc_r2min = 1.0e-4_R_

    ! New dataset
    if (Sfc_idread >= 1) then ! only when new dataset is requested

       ! Initialize data by domain-wide constant properties
       Sfc_tmps2d(:,:) = Sfc_tmp
       Sfc_jsfc2d(:,:) = Sfc_mtype
       do i = 1, Sfc_NPAR
          Sfc_psfc2d(:,:,i) = Sfc_param(i)
       end do

       ! Import data
       if (irank == Sfc_iroot .and. Sfc_iud >= 1) call mcarSfc__read(mbswap) ! root PE only
#if UseMPI == 1
       call mcarSfc__MPI_bcast_read(Sfc_iroot)
#endif

       ! Preprocess data
       call mcarSfc__prep_optProp()
       call mcarSfc__prep_temp(tmpmin, tmpmax)
    end if

  end subroutine mcarSfc__prep_job


  !+
  ! Import data from file
  !-
  subroutine mcarSfc__read(mbswap) 

    integer, intent(in) :: mbswap
    integer :: nmax, i

    ! Initialize
    nmax = 0
    do i = 1, 4
       if (Sfc_mbrdf(i) == 1) nmax = max(nmax, Sfc_npsfc(i))
    end do
       
    ! Read in
    if (Sfc_mfmt == 0) then ! text
       call mcarUtl__txtSfc_read(Sfc_iud, Sfc_nxb, Sfc_nyb, nmax, Sfc_tmps2d, Sfc_psfc2d, Sfc_jsfc2d)
    else
       call mcarUtl__binSfc_read(Sfc_iud, Sfc_nxb, Sfc_nyb, nmax, Sfc_tmps2d, Sfc_psfc2d, &
            & Sfc_jsfc2d, mbswap, Sfc_idread)
    end if
    
  end subroutine mcarSfc__read


#if UseMPI == 1
  !+
  ! MPI broadcast data read in mcarSfc__read
  !-
  subroutine mcarSfc__MPI_bcast_read(iroot)

    include 'inc_mpi.f90'
    integer, intent(in) :: iroot ! root PE index
    integer :: ierr, nxyb, ip

    nxyb = Sfc_nxb * Sfc_nyb
    call MPI_Bcast(Sfc_tmps2d(:,:), nxyb, MPI_R_,      iroot, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(Sfc_jsfc2d(:,:), nxyb, MPI_INTEGER, iroot, MPI_COMM_WORLD, ierr)
    do ip = 1, Sfc_NPAR
       call MPI_Bcast(Sfc_psfc2d(:,:,ip), nxyb, MPI_R_, iroot, MPI_COMM_WORLD, ierr)
    end do

  end subroutine mcarSfc__MPI_bcast_read
#endif


  !+
  ! Setup surface parameters: make LUT for albedo etc
  !-
  subroutine mcarSfc__prep_optProp()

    real(R_),  parameter :: DLTMIN = 0.02_R_ ! min of delta(uz0) for searching albmin
    real(R_)  :: pac(10), pmin(4, 10), pmax(4, 10)
    real(R_)  :: fawrk(Sfc_nsco/2), w1wrk(Sfc_nsco/2), w2wrk(Sfc_nsco/2) !(AUTO)
    real(R_)  :: agwrk(Sfc_nsco/2), avwrk(Sfc_nsco/2)
    real(R_)  :: uz1wrk(Sfc_nqpot), fgwrk(Sfc_nqpot, Sfc_nsco/2), fvwrk(Sfc_nqpot, Sfc_nsco/2)
    real(R_)  :: ag, albbar, albmin, av, cosq0, cosq1, delta, duz0, duz1, fg, fv
    integer   :: i, iam, iqsfc, isfc, isfcnew, isco, isuz, iuz0, ixb, iyb
    integer   :: nphi, nuz0
    real(R_)  :: pacnew, sinq0, sinq1

    ! Initialize
    nphi = Sfc_nqpot / 2 + 1      ! # of quadrature points for phi
    pmin(1, 1) = 0.0_R_          ! Lambertian
    pmax(1, 1) = 1.0_R_
    pmin(2, 1) = 0.0_R_          ! DSM
    pmin(2, 2) = 0.0_R_
    pmin(2, 3) = 0.0_R_
    pmin(2, 4) = 0.0_R_
    pmin(2, 5) = 0.0_R_
    pmax(2, 1) = 1.0_R_
    pmax(2, 2) = 1.0_R_
    pmax(2, 3) = 100.0_R_
    pmax(2, 4) = 100.0_R_
    pmax(2, 5) = 100.0_R_
    pmin(3, 1) = 0.0_R_          ! RPV
    pmin(3, 2) = 0.0_R_
    pmin(3, 3) = -1.0_R_
    pmin(3, 4) = 0.0_R_
    pmax(3, 1) = 100.0_R_
    pmax(3, 2) = 1.0_R_
    pmax(3, 3) = 1.0_R_
    pmax(3, 4) = RLRG_
    pmin(4, 1) = -1.0_R_         ! LSRT
    pmin(4, 2) = -1.0_R_
    pmin(4, 3) = -1.0_R_
    pmax(4, 1) = 100.0_R_
    pmax(4, 2) = 100.0_R_
    pmax(4, 3) = 100.0_R_
    do iyb = 1, Sfc_nyb
       do ixb = 1, Sfc_nxb
          isfc = Sfc_jsfc2d(ixb, iyb)
          isfc = max(0, min(4, isfc))
          if (isfc >= 1) then
             do i = 1, Sfc_npsfc(isfc)
                pac(i) = Sfc_psfc2d(ixb, iyb, i)
                call check_iR('mcarSfc__prep_optProp: Sfc_psfc2d', pac(i), pmin(isfc, i), pmax(isfc, i))
             end do
             if (isfc == 1 .and. pac(1) <= 0.0_R_) then ! simplifications
                isfc = 0         ! black
             else if (isfc == 2 .and. pac(2) > 0.99999_R_) then
                isfc = 1         ! Lambertian
             end if
          end if
          Sfc_jsfc2d(ixb, iyb) = isfc
       end do
    end do

    ! LUTs for LSRT model
    duz0 = (1.0_R_ - Sfc_uz0min) / (Sfc_nsuz - 1)
    do isuz = 1, Sfc_nsuz
       cosq0 = -min(1.0_R_, Sfc_uz0min + duz0 * (isuz - 1))
       sinq0 = sqrt(1.0_R_ - cosq0 * cosq0)
       call mcarLsrt__albKernel(sinq0, cosq0, Sfc_nqpot, nphi, ag, av)
       Sfc_salut1(isuz) = ag
       Sfc_salut2(isuz) = av
    end do
    duz1 = 1.0_R_ / real(Sfc_nqpot, R_)
    do iqsfc = 1, Sfc_nqpot
       uz1wrk(iqsfc) = real(iqsfc - 0.5_R_) * duz1
    end do
    nuz0 = Sfc_nulsrt
    duz0 = (1.0_R_ - Sfc_uz0min) / (nuz0 - 1)
    do iuz0 = 1, nuz0
       cosq0 = -min(1.0_R_, Sfc_uz0min + duz0 * (iuz0 - 1))
       sinq0 = sqrt(1.0_R_ - cosq0 * cosq0)
       call mcarLsrt__albKernel(sinq0, cosq0, Sfc_nqpot, nphi, ag, av)
       agwrk(iuz0) = ag
       avwrk(iuz0) = av
       do iqsfc = 1, Sfc_nqpot
          cosq1 = uz1wrk(iqsfc)
          sinq1 = sqrt(1.0_R_ - cosq1 * cosq1)
          call mcarLsrt__BRFKernel(sinq0, cosq0, sinq1, cosq1, nphi, fg, fv)
          fgwrk(iqsfc, iuz0) = fg
          fvwrk(iqsfc, iuz0) = fv
       end do
    end do

    ! Loop for pixels
    do iyb = 1, Sfc_nyb
       do ixb = 1, Sfc_nxb
          Sfc_salbs(1, ixb, iyb) = 0.0_R_
          Sfc_salbs(2, ixb, iyb) = 0.0_R_

          isfcnew = Sfc_jsfc2d(ixb, iyb)
          delta = abs(isfc - isfcnew)
          isfc = isfcnew
          if (isfc == 0) cycle ! black
          do i = 1, Sfc_npsfc(isfc) ! parameters setting
             pacnew = Sfc_psfc2d(ixb, iyb, i)
             delta = delta + abs(pac(i) - pacnew)
             pac(i) = pacnew
          end do

          !     Duplications
          if (ixb >= 2 .and. delta < RSML_) then
             Sfc_salbs(1, ixb, iyb) = Sfc_salbs(1, ixb-1, iyb)
             Sfc_salbs(2, ixb, iyb) = Sfc_salbs(2, ixb-1, iyb)
             if (isfc >= 2) then
                do isco = 1, Sfc_nsco
                   Sfc_scolut(isco, ixb, iyb) = Sfc_scolut(isco, ixb-1, iyb)
                end do
             end if
             cycle
          end if

          !     Make LUTs for albedo & coefficients
          if (isfc <= 1) then ! Lambertian
             Sfc_salbs(1, ixb, iyb) = pac(1)
             Sfc_salbs(2, ixb, iyb) = pac(1)
          else
             if (isfc == 2) then ! DSM
                nuz0 = Sfc_nudsm
                call mcarDsm__prep_makeTab(pac, Sfc_r2min, Sfc_rrmax, Sfc_rrexp, Sfc_uz0min, &
                     &                 nuz0, Sfc_nqpot, nphi, fawrk, w1wrk, albmin, iam)
                call mcarSfc__albStat(isfc, pac, fawrk, avwrk, Sfc_uz0min, &
                     &                 nuz0, Sfc_nqpot, DLTMIN, albbar, albmin, iam)
                do iuz0 = 1, nuz0
                   Sfc_scolut(iuz0,        ixb, iyb) = log(fawrk(iuz0)) ! log(A)
                   Sfc_scolut(iuz0 + nuz0, ixb, iyb) = log(w1wrk(iuz0)) ! log(b')
                end do
             else if (isfc == 3) then ! RPV
                nuz0 = Sfc_nurpv
                call mcarRpv__prep_makeTab(pac, Sfc_rrmax, Sfc_rrexp, Sfc_uz0min, nuz0, &
                     &                 Sfc_nqpot, nphi, fawrk, w2wrk, w1wrk, albmin, iam)
                call mcarSfc__albStat(isfc, pac, fawrk, avwrk, Sfc_uz0min, &
                     &                 nuz0, Sfc_nqpot, DLTMIN, albbar, albmin, iam)
                do iuz0 = 1, nuz0
                   Sfc_scolut(iuz0,          ixb, iyb) = fawrk(iuz0) ! A
                   Sfc_scolut(iuz0 +   nuz0, ixb, iyb) = w1wrk(iuz0) ! b'
                   Sfc_scolut(iuz0 + 2*nuz0, ixb, iyb) = w2wrk(iuz0) ! a1
                end do
             else             ! LSRT
                nuz0 = Sfc_nulsrt
                call mcarLsrt__prep_makeTab(pac, Sfc_rrmax, Sfc_rrexp, Sfc_uz0min, nuz0, &
                     &                 Sfc_nqpot, uz1wrk, agwrk, avwrk, fgwrk, &
                     &                 fvwrk, w2wrk, w1wrk, albmin, iam)
                call mcarSfc__albStat(isfc, pac, agwrk, avwrk, Sfc_uz0min, &
                     &                 nuz0, Sfc_nqpot, DLTMIN, albbar, albmin, iam)
                do iuz0 = 1, nuz0
                   Sfc_scolut(iuz0,        ixb, iyb) = w1wrk(iuz0) ! b'
                   Sfc_scolut(iuz0 + nuz0, ixb, iyb) = w2wrk(iuz0) ! a1
                end do
             end if
             Sfc_salbs(1, ixb, iyb) = albbar
             Sfc_salbs(2, ixb, iyb) = albmin
          end if

       end do
    end do

  end subroutine mcarSfc__prep_optProp


  !+
  ! Get temperature min & max
  !-
  subroutine mcarSfc__prep_temp(tmpmin, tmpmax) 

    real(R_), intent(out) :: tmpmin, tmpmax
    tmpmin = minval(Sfc_tmps2d)
    tmpmax = maxval(Sfc_tmps2d)

  end subroutine mcarSfc__prep_temp


  !+
  ! Angular distribution function (ADF) to be used for local estimate of radiance
  !// Note: ADF is defined as the PDF (1/steradian) for angular (re)distribution
  !   -thermal emission at surface : ADF = |cos(q1)| * (Ems / Ems_bar) / pi
  !   -reflection at surface       : ADF = |cos(q1)| * BRDF(u0,f0,u1,f1) / BS_albedo
  !-
  function mcarSfc__angDistr_R(imod, isub, pac, dir1, dir2, ix, iy) result(pdf)

    integer,  intent(in) :: imod
    !// imod = 2 : thermal source emission at surface, 4 : surface reflection
    integer,  intent(in) :: isub
    !// isub = isfc (0:black, 1:Lambertian, 2:DSM, 3:RPV, 4:LSRT)
    real(R_), intent(in) :: pac(:)
    real(R_), intent(in) :: dir1(:)
    real(R_), intent(in) :: dir2(:)
    integer,  intent(in) :: ix, iy
    real(R_)  :: pdf ! result, PDF
    real(R_),  parameter :: API = 1.0_R_ / PI_
    real(R_),  parameter :: PDFMAX = 1.0e+5_R_
    real(R_)  :: r, alb, yfa

    ! Source emission from surface
    if (imod == 2) then
       pdf = abs(dir2(3)) * mcarSfc__emittance_R(dir2(3), ix, iy, isub, pac)   

       ! Reflection at surface
    else if (dir2(3) <= 0.0_R_ .or. isub <= 0) then ! black
       pdf = 0.0_R_
    else if (isub == 1) then ! Lambertian
       pdf = abs(dir2(3)) * API
    else
       if (isub == 2) then ! DSM
          !nco = 1
          !pac5o = pac(5) ! save the original
          !if (Pho_chi < 0.7_R_) pac(5) = max(0.03_R_, pac5o) ! smooth BRDF (temporalily)
          !call mcarSfc__intpOptProp(nco, isub, ix, iy, -dir1(3), yfa1, yfa2, yco)
          !alb = mcarSfc__bsAlbedo_R(isub, pac, yfa1, yfa2, -dir1(3)) ! albedo (can be > 1)
          !// Albedo should be modified for better energy conservation if the surface
          !   roughness is modified. This requires a complicated calculation of the albedo.
          alb = pac(6)
          r = mcarDsm__BRDFpA_R(pac, alb, dir1, dir2)
          !pac(5) = pac5o ! recover the original
       else if (isub == 3) then ! RPV
          yfa = pac(5)
          r = mcarRpv__BRDFpA_R(pac, yfa, dir1, dir2)
       else                     ! LSRT
          r = mcarLsrt__BRDF_R(pac, dir1, dir2) / pac(6)
       end if
       pdf = min(PDFMAX, abs(dir2(3)) * r)
    end if

  end function mcarSfc__angDistr_R


  !+
  ! Black-sky albedo of DSM (isfc=2), RPV (3), or LSRT (4) model
  !-
  function mcarSfc__bsAlbedo_R(isfc, pac, fa1, fa2, uz0) result(y)

    real(R_), intent(in) :: pac(:) ! pac(1-3) : parameter packet (definitions are different by isfc
    real(R_), intent(in) :: uz0 ! absolute value of cos(theta0), where theta0 is the zenith angle
    real(R_), intent(in) :: fa1, fa2 ! function A, 1st & 2nd kinds
    integer,  intent(in) :: isfc ! surface model ID (2 for DSM, 3 for RPV, 4 for LSRT)
    real(R_)  :: y ! result, black-sky albedo
    ! Note: The output albedo can be larger than unity!
    real(R_),  parameter :: ALBEPS = 1.0e-14_R_ ! minimum albedo for 1/albedo calculation

    ! DSM
    if (isfc == 2) then
       y = pac(2) * pac(1) + (1.0_R_ - pac(2)) * (fa1 / uz0)

       ! RPV
    else if (isfc == 3) then
       y = max(ALBEPS, fa1 * uz0**(pac(2) - 1.0_R_))

       ! LSRT
    else if (isfc == 4) then
       y = max(ALBEPS, pac(1) + pac(2) * fa1 + pac(3) * fa2)

       ! Lambertian
    else if (isfc == 1) then
       y = pac(1)

       ! Black
    else
       y = 0.0_R_
    end if

  end function mcarSfc__bsAlbedo_R


  !+
  ! Get albedo average & minimum
  !-
  subroutine mcarSfc__albStat(isfc, pac, falut1, falut2, uz0min, nuz0, &
       &     nquad, dltmin, albbar, albmin, iam)

    integer,  intent(in) :: isfc ! surface model ID (2 for DSM, 3 for RPV, 4 for LSRT)
    real(R_), intent(in) :: pac(:) ! pac(1-5): parameter packet (definitions are different by isfc
    real(R_), intent(in) :: falut1(:) ! LUT of function A (1st kind) for iuz0=1-nuz0
    real(R_), intent(in) :: falut2(:) ! LUT of function A (2nd kind) for iuz0=1-nuz0
    integer,  intent(in) :: nquad ! number of quadrature points for albedo avraging
    integer,  intent(in) :: nuz0  ! number of cosines of incoming zenith angle (should be > 3)
    real(R_), intent(in) :: uz0min ! minimum of cosine of incident zenith angle
    real(R_), intent(in) :: dltmin ! minimum of delta(uz0) for searching albmin (usualy 0.02_R_)
    real(R_), intent(out) :: albbar ! albedo average
    real(R_), intent(inout) :: albmin ! albedo minimum
    integer,  intent(inout) :: iam ! iuz0 corresponding to the minimum albedo
    real(R_),  parameter :: GOLD = 0.61803399_R_, GOLDC = 1.0_R_ - gold
    integer,   parameter :: KNQUAD = 200
    real(R_), save :: xq(KNQUAD), wq(KNQUAD)
    integer,  save :: nq = 0
    real(R_)  :: alb, buz0, fa1, fa2
    integer   :: iq, iuz0
    integer   :: nqtmp
    real(R_)  :: tsum, uz0, x, x0, x1, x2, x3, y1, y2

    ! Initialization
    nqtmp = min(KNQUAD, nquad)
    if (nq /= nqtmp) then
       nq = nqtmp
       call gaussLegen(nq, xq, wq)
       do iq = 1, nq
          xq(iq) = 0.5_R_ * (1.0_R_ + xq(iq))
       end do                 ! xq is set as 0 to 1
    end if

    ! Average albedo: integration of cos(q0)*albedo for the hemisphere
    tsum = 0.0_R_
    do iq = 1, nq
       x = xq(iq)
       !fa1 = intplutxeq_R(1, 3, falut1, 1, nuz0, uz0min, 1.0_R_, x)
       !if (isfc == 4) fa2 = intplutxeq_R(0, 3, falut2, 1, nuz0, uz0min, 1.0_R_, x)
       fa1 = intp_tabCub_eqX_R(falut1, nuz0, uz0min, 1.0_R_, x)
       if (isfc == 4) fa2 = intp_tabCub_eqX_R(falut2, nuz0, uz0min, 1.0_R_, x)
       alb = mcarSfc__bsAlbedo_R(isfc, pac, fa1, fa2, x)
       tsum = tsum + wq(iq) * x * min(1.0_R_, alb)
    end do
    albbar = min(1.0_R_, 2.0_R_ * tsum) ! average albedo

    ! Minimum albedo (first guess)
    buz0 = (1.0_R_ - uz0min) / (nuz0 - 1)
    if (iam == 0) then
       albmin = 1.0_R_
       iam = 1
       do iuz0 = 1, nuz0
          uz0 = uz0min + buz0 * (iuz0 - 1)
          fa1 = falut1(iuz0)
          fa2 = falut2(iuz0)
          alb = mcarSfc__bsAlbedo_R(isfc, pac, fa1, fa2, uz0)
          if (alb < albmin) then
             albmin = alb
             iam = iuz0
          end if
       end do
    end if

    ! Minimum albedo using the golden section search
    if (iam > 1 .and. iam < nuz0) then
       x1 = uz0min + buz0 * (iam - 1)
       x0 = x1 - buz0
       x3 = x1 + buz0
       x2 = x1 + GOLDC * (x3 - x1)
       !fa1 = intplutxeq_R(1, 3, falut1, 1, nuz0, uz0min, 1.0_R_, x2)
       !if (isfc == 4) fa2 = intplutxeq_R(0, 3, falut2, 1, nuz0, uz0min, 1.0_R_, x2)
       fa1 = intp_tabCub_eqX_R(falut1, nuz0, uz0min, 1.0_R_, x2)
       if (isfc == 4) fa2 = intp_tabCub_eqX_R(falut2, nuz0, uz0min, 1.0_R_, x2)
       y2 = mcarSfc__bsAlbedo_R(isfc, pac, fa1, fa2, x2)
       y1 = albmin
       do
          if (x3 - x0 < dltmin) exit
          if (y2 < y1) then
             x0 = x1
             x1 = x2
             x2 = GOLD * x1 + GOLDC * x3
             y1 = y2
             !fa1 = intplutxeq_R(1, 3, falut1, 1, nuz0, uz0min, 1.0_R_, x2)
             !if (isfc == 4) fa2 = intplutxeq_R(0, 3, falut2, 1, nuz0, uz0min, 1.0_R_, x2)
             fa1 = intp_tabCub_eqX_R(falut1, nuz0, uz0min, 1.0_R_, x2)
             if (isfc == 4) fa2 = intp_tabCub_eqX_R(falut2, nuz0, uz0min, 1.0_R_, x2)
             y2 = mcarSfc__bsAlbedo_R(isfc, pac, fa1, fa2, x2)
          else
             x3 = x2
             x2 = x1
             x1 = GOLD * x2 + GOLDC * x0
             y2 = y1
             !fa1 = intplutxeq_R(1, 3, falut1, 1, nuz0, uz0min, 1.0_R_, x1)
             !if (isfc == 4) fa2 = intplutxeq_R(0, 3, falut2, 1, nuz0, uz0min, 1.0_R_, x1)
             fa1 = intp_tabCub_eqX_R(falut1, nuz0, uz0min, 1.0_R_, x1)
             if (isfc == 4) fa2 = intp_tabCub_eqX_R(falut2, nuz0, uz0min, 1.0_R_, x1)
             y1 = mcarSfc__bsAlbedo_R(isfc, pac, fa1, fa2, x1)
          end if
       end do
       albmin = min(y1, y2)
    end if

  end subroutine mcarSfc__albStat


  !+
  ! Normalized angular PDF of surface emission
  !-
  function mcarSfc__emittance_R(uz2, ixb, iyb, isub, pac) result(adf)

    integer,  intent(in) :: isub
    integer,  intent(in) :: ixb, iyb
    real(R_), intent(in) :: uz2
    real(R_), intent(in) :: pac(:)
    real(R_)  :: adf
    real(R_),  parameter :: API = 1.0_R_ / PI_
    real(R_)  :: yco(2)
    real(R_)  :: alb, ems, emsbar, yfa1, yfa2

    ! Downward
    if (uz2 < 0.0_R_) then
       adf = 0.0_R_

       ! Black or Lambertian
    else if (isub <= 1) then
       adf = API

       ! DSM, RPV or LSRT
    else
       call mcarSfc__intpOptProp(0, isub, ixb, iyb, uz2, yfa1, yfa2, yco)
       alb = mcarSfc__bsAlbedo_R(isub, pac, yfa1, yfa2, uz2) ! albedo (can be > 1)
       ems = 1.0_R_ - min(1.0_R_, alb)
       emsbar = max(REPS_, 1.0_R_ - Sfc_salbs(1, ixb, iyb))
       adf = API * (ems / emsbar)
    end if

  end function mcarSfc__emittance_R


  !+
  ! Thermal emission from horizontal plane: determination of source direction
  !-
  function mcarSfc__fMC_emit_1R(ixb, iyb, isfc, pac) result(dir)

    integer,  intent(in) :: isfc     ! surface model ID
    !// (0: black; 1: Lambertian; 2: DSM; 3: RPV, 4: LSRT)
    integer,  intent(in) :: ixb, iyb ! pixel indexes
    real(R_), intent(in) :: pac(:)   ! parameter packet (definitions are different by isfc
    real(R_)  :: dir(3) ! result, outgoing direction vector (upward)
    real(R_)  :: yco(2), alb, cosf, ems, emsmax, uz, r2, sinf, sinq, yfa1, yfa2
    integer   :: itrial

    ! Random azimuth
    call rand_point_circle(1.0e-12_R_, r2, sinf, cosf)

    ! Lambertian emission
    if (isfc <= 1) then
       dir(:) = geo3D_UVec_LambertH_1R(r2, sinf, cosf)

       ! Non-Lambertian
    else
       emsmax = (1.0_R_ - Sfc_salbs(2, ixb, iyb)) * 1.1_R_
       !// correction because min(albedo) was roughly estimated

       ! Get a random direction
       do itrial = 1, 100
          r2 = mseq_rand_R()
          uz = sqrt(r2) ! Lambertian (tentative)
          call mcarSfc__intpOptProp(0, isfc, ixb, iyb, uz, yfa1, yfa2, yco)
          alb = mcarSfc__bsAlbedo_R(isfc, pac, yfa1, yfa2, uz) ! albedo
          ems = 1.0_R_ - alb
          if (mseq_rand_R() * emsmax < ems) exit ! emissivity scaling
       end do

       ! Vector components
       sinq = sqrt(1.0_R_ - r2)
       dir(1) = sinq * cosf
       dir(2) = sinq * sinf
       dir(3) = uz
    end if

  end function mcarSfc__fMC_emit_1R


  !+
  ! Interpolation of surface parameters by linear/3rd-order Lagrangian polynomial using LUT
  !-
  subroutine mcarSfc__intpOptProp(nco, isfc, ixb, iyb, uint, yfa1, yfa2, yco)

    integer,  intent(in) :: isfc
    integer,  intent(in) :: ixb, iyb
    integer,  intent(in) :: nco
    !   nco=0: retrieve yfa1 (and yfa2 if isfc=4)
    !   nco=1: retrieve yfa* and yco(1)
    !   nco=2: retrieve yfa*, yco(1), yco(2)
    real(R_), intent(in) :: uint
    real(R_), intent(out) :: yfa1, yfa2
    real(R_), intent(out) :: yco(:)
    integer   :: ico, ilut, iu, nu
    real(R_)  :: rat, rat3, unorm, zu

    ! Initializations
    unorm = (uint - Sfc_uz0min) / (1.0_R_ - Sfc_uz0min)
    unorm = max(0.0_R_, min(1.0_R_, unorm)) ! 0 to 1

    ! DSM or RPV
    if (isfc <= 3) then
       if (isfc == 2) then
          nu = Sfc_nudsm
       else
          nu = Sfc_nurpv
       end if
       zu = unorm * real(nu - 1, R_)
       iu = min(nu - 1, int(zu) + 1) ! 1 to nu-1
       if (iu >= 2 .and. iu <= nu - 2) then
          ilut = iu - 1
          rat3 = zu - iu + 2.0_R_
       else if (iu <= 1) then
          ilut = 1
          rat3 = zu - iu + 1.0_R_
       else
          ilut = nu - 3
          rat3 = zu - iu + 3.0_R_
       end if
       yfa1 = intp_cubic_eqX_R(rat3 - 1.0_R_, Sfc_scolut(ilut:ilut+3, ixb, iyb))
       if (isfc == 2) yfa1 = exp(yfa1)
       ilut = ilut + nu
       do ico = 1, nco
          yco(ico) = intp_cubic_eqX_R(rat3 - 1.0_R_, Sfc_scolut(ilut:ilut+3, ixb, iyb))
          if (isfc == 2) yco(ico) = exp(yco(ico))
          ilut = ilut + nu
       end do

       ! LSRT
    else
       zu = unorm * real(Sfc_nsuz - 1, R_)
       iu = min(Sfc_nsuz - 1, int(zu) + 1)
       rat = zu - real(iu - 1, R_)
       yfa1 = (1.0_R_ - rat) * Sfc_salut1(iu) + rat * Sfc_salut1(iu + 1)
       yfa2 = (1.0_R_ - rat) * Sfc_salut2(iu) + rat * Sfc_salut2(iu + 1)
       if (nco >= 1) then
          nu = Sfc_nulsrt
          zu = unorm * real(nu - 1, R_)
          iu = min(nu - 1, int(zu) + 1) ! 1 to nu-1
          if (iu >= 2 .and. iu <= nu - 2) then
             ilut = iu - 1
             rat3 = zu - iu + 1.0_R_
          else if (iu <= 1) then
             ilut = 1
             rat3 = zu - iu
          else
             ilut = nu - 3
             rat3 = zu - iu + 2.0_R_
          end if
          do ico = 1, nco
             yco(ico) = intp_cubic_eqX_R(rat3, Sfc_scolut(ilut:ilut+3, ixb, iyb))
             ilut = ilut + nu
          end do
       end if
    end if

  end subroutine mcarSfc__intpOptProp


  !+
  ! A new random reflection direction and angular PDF
  !-
  subroutine mcarSfc__newDirec(isfc, pac, yco, odir, rdir, pps)

    integer,  intent(in) :: isfc     ! surface BRDF model index
    real(R_), intent(in) :: pac(:)   ! surface parameter package
    real(R_), intent(in) :: yco(:)   ! coefficients for scattering
    real(R_), intent(in) :: odir(:)  ! incoming, original direction vector
    real(R_), intent(out) :: rdir(:) ! outgoing, reflection direction vector
    real(R_), intent(out) :: pps     ! angular PDF (/steradian)
    real(R_),  parameter :: API = 1.0_R_ / PI_
    real(R_),  parameter :: PPSMAX = 1.0e+5_R_
    real(R_) :: a

    if (isfc <= 1) then ! Lambertian
       call mcarDsm__newDir_diff(rdir, pps)
    else if (isfc == 2) then ! DSM
       a = pac(6) ! albedo
       if (mseq_rand_R() * a < pac(2) * pac(1)) then ! diffuse reflection
          call mcarDsm__newDir_diff(rdir, pps)
       else ! specular reflection
          call mcarDsm__newDir_spec(pac(3), pac(4), pac(5), yco(1), Sfc_r2min, a, odir, rdir, pps)
       end if
    else if (isfc == 3) then ! RPV
       a = pac(5) ! function A
       call mcarRpv__newDir(odir, pac, a, yco(2), yco(1), rdir, pps)
    else ! LSRT
       a = pac(6) ! albedo
       call mcarLsrt__newDir(odir, pac, a, yco(2), yco(1), rdir, pps)
    end if
    pps = min(PPSMAX, pps)
       
  end subroutine mcarSfc__newDirec

end module mcarSfc
