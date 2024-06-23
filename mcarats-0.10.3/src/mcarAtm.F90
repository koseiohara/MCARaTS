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
! Module for atmosphere
!-
module mcarAtm 

  use globals
#if UseMPI == 1
  use globals_MPI
#endif
  use hparx
  use mcarUtl
  use mcarSca
  use mcarPho
  implicit none
  private

  ! Public
  public :: mcarAtm__user_init
  public :: mcarAtm__init
  public :: mcarAtm__user_job
  public :: mcarAtm__prep_job
  public :: mcarAtm__check_opt3DStats
  public :: mcarAtm__final
  public :: mcarAtm__tune_jobs
  public :: mcarAtm__set_tech
  public :: mcarAtm__set_extMin
  public :: mcarAtm__prep_opt3D
  public :: mcarAtm__prep_optTot
  public :: mcarAtm__prep_temp
  public :: mcarAtm__prep_extMax
  public :: mcarAtm__angDistr_R
  public :: mcarAtm__rt_escape
  public :: mcarAtm__rt_bound1D
  public :: mcarAtm__rt_bound3D
  public :: mcarAtm__rt_renewLocXY
  public :: mcarAtm__rt_path1D_R
  public :: mcarAtm__rt_path3D
  public :: mcarAtm__rt_pathSup

  ! (Partly-Public) User variables in namelist (1) for initialization
  character(256), save :: Atm_inpfile = ' ' ! file name for input of 3D otpical properties
  integer,  save :: Atm_mfmt = 1   ! flag for a format of input file (0=text, 1=binary)
  integer,  save :: Atm_nx   = 1   ! # of X grid points
  integer,  save :: Atm_ny   = 1   ! # of Y grid points
  integer,  save :: Atm_nz   = 1   ! # of Z grid points (layers)
  integer,  save :: Atm_iz3l = 1   ! layer index for the lowest 3-D layer
  integer,  save :: Atm_nz3  = 0   ! # of 3-D inhomogeneous layers
  integer,  save :: Atm_nkd  = 1   ! # of k-distribution terms
  integer,  save :: Atm_np1d = 1   ! # of scattering components in 1-D medum
  integer,  save :: Atm_np3d = 1   ! # of scattering components in 3-D medum
  integer,  save :: Atm_mtprof = 0 ! flag for temperature profile (0=const, 1=linear)
  integer,  save :: Atm_nqlay = 10 ! # of Gaussian quadrature points per layer
  namelist /mcarAtm_nml_init/Atm_inpfile, Atm_mfmt, Atm_iz3l, Atm_nx, Atm_ny, Atm_nz, &
       & Atm_nz3, Atm_nkd, Atm_np1d, Atm_np3d, Atm_mtprof, Atm_nqlay

  ! (Private) User variables in namelist (2) for core calculations
  integer :: i
  integer, parameter :: KNKD  = 100
  integer, parameter :: KNP1D = 100
  integer, parameter :: KNP3D = 100 
  integer, parameter :: KNZ   = 3000
  integer,  save :: Atm_idread = 0     ! location of data to be read (0: no data to be read)
  real(R_), save :: Atm_dx = 1.0e+4_R_ ! X size of cell
  real(R_), save :: Atm_dy = 1.0e+4_R_ ! Y size of cell
  real(R_), save :: Atm_wkd0(KNKD) = (/1.0_R_, (0.0_R_, i = 1, KNKD - 1)/) ! k-distribution weights
  real(R_), save :: Atm_zgrd0(KNZ+1) = (/(10.0_R_*i, i = 1, KNZ+1)/) ! Z (height) at layer interfaces
  real(R_), save :: Atm_tmp1d(KNZ+1) = (/(300.0_R_,  i = 1, KNZ+1)/) ! temperatures
  real(R_), save :: Atm_ext1d(KNZ, KNP1D) ! extinction coefficients
  real(R_), save :: Atm_omg1d(KNZ, KNP1D) ! single scattering albedos
  real(R_), save :: Atm_apf1d(KNZ, KNP1D) ! phase function specification parameters
  real(R_), save :: Atm_abs1d(KNZ, KNKD)  ! absorption coefficients
  real(R_), save :: Atm_fext1d(KNP1D) = (/(1.0_R_, i = 1, KNP1D)/) ! scaling factor for Atm_ext1d
  real(R_), save :: Atm_fext3d(KNP3D) = (/(1.0_R_, i = 1, KNP3D)/) ! scaling factor for Atm_ext3d
  real(R_), save :: Atm_fabs1d = 1.0_R_ ! scaling factor for Atm_abs1d
  real(R_), save :: Atm_fabs3d = 1.0_R_ ! scaling factor for Atm_abs3d
  real(R_), save :: Atm_fsupg  = 2.0_R_ ! parameter for adaptive super-cell formation
  real(R_), save :: Atm_fsupi  = 1.3_R_ ! parameter for adaptive super-cell formation
  integer,  save :: Atm_nzzsup = 50     ! max # of layers in a single super-cell
  real(R_), save :: Atm_taumin = 2.0_R_ ! min column optical thickness for collision forcing 
  namelist /mcarAtm_nml_job/Atm_idread, Atm_dx, Atm_dy, Atm_zgrd0, Atm_tmp1d, Atm_ext1d, Atm_omg1d, &
       & Atm_apf1d, Atm_wkd0, Atm_abs1d, Atm_fext1d, Atm_fext3d, Atm_fabs1d, Atm_fabs3d, &
       & Atm_taumin, Atm_fsupg, Atm_fsupi, Atm_nzzsup

  ! (Public, read only)
  integer,  save :: Atm_iz3u = 0
  real(R_), save :: Atm_facx, Atm_facy
  real(R_), save :: Atm_xmax, Atm_ymax
  integer,  save, allocatable :: Atm_mlay(:) !(0:nz+1)
  !// Atm_mlay(iz) : layer index
  !    =0 : underlying surface
  !    =1 : horizontally-homogeneous layer (absorbing)
  !    =2 : horizontally-homogeneous layer (absorption + scattering)
  !    =3 : inhomogeneous layer
  !    =4 : above the domain
  real(R_), save, allocatable :: Atm_xgrd(:) !(0:nx)
  real(R_), save, allocatable :: Atm_ygrd(:) !(0:ny)
  real(R_), save, allocatable :: Atm_zgrd(:) !(0:nz+1) Z-grid values, heights (m) of layer interfaces
  real(R_), save, allocatable :: Atm_zdep(:) !(nz+1) layer depth (m)
  real(R_), save, allocatable :: Atm_wkd(:)  !(nkd)
  public :: Atm_nkd, Atm_nx, Atm_ny, Atm_nz, Atm_nz3, Atm_iz3l, Atm_iz3u
  public :: Atm_facx, Atm_facy, Atm_xmax, Atm_ymax, Atm_wkd
  public :: Atm_xgrd, Atm_ygrd, Atm_zgrd, Atm_zdep, Atm_mlay
  integer,  save :: Atm_nppomax = 4
  integer,  save, allocatable :: Atm_nppo(:)         !(nz)
  real(R_), save, allocatable :: Atm_extvar(:)       !(nz) normalized variance of extinctin
  real(R_), save, allocatable :: Atm_extmin(:)       !(nz)
  real(R_), save, allocatable :: Atm_extscl(:,:)     !(nchi,nkd)
  real(R_), save, allocatable :: Atm_apfp3d(:,:,:,:) !(nx,ny,nz,nppo)
  real(R_), save, allocatable :: Atm_scap3d(:,:,:,:) !(nx,ny,nz,nppo) ! component scattering coeff.
  real(R_), save, allocatable :: Atm_scat3d(:,:,:,:) !(nx,ny,nz,nchi) ! total scattering coeff.
  real(R_), save, allocatable :: Atm_abst3d(:,:,:,:) !(nx,ny,nz,nkd)  ! total absorption coeff.
  real(R_), save, allocatable :: Atm_tmpa3d(:,:,:)   !(nx,ny,nz+1)
  public :: Atm_mtprof, Atm_nppo, Atm_apfp3d, Atm_tmpa3d, Atm_extvar, Atm_extmin, Atm_extscl
  public :: Atm_scap3d, Atm_scat3d, Atm_abst3d
  integer,  save :: Atm_nzs
  integer,  save, allocatable :: Atm_jzstab(:) !(0:nz+1)
  real(R_), save, allocatable :: Atm_xsup(:,:) !(0:nxs,nzs)
  real(R_), save, allocatable :: Atm_ysup(:,:) !(0:nys,nzs)
  real(R_), save, allocatable :: Atm_zsup(:)   !(0:nzs)
  integer,  save, allocatable :: Atm_jsup(:)   !(0:nzs+1)
  integer,  save, allocatable :: Atm_jzupr(:)  !(0:nzs)
  integer,  save, allocatable :: Atm_nxs(:), Atm_nys(:) !(nzs)
  integer,  save, allocatable :: Atm_nxx(:), Atm_nyy(:) !(nzs)
  real(R_), save, allocatable :: Atm_extmax(:,:,:,:) !(nxs,nys,nzs,nchi)
  public :: Atm_nzs, Atm_jzstab
  public :: Atm_xsup, Atm_ysup, Atm_zsup, Atm_jsup, Atm_jzupr, Atm_nxs, Atm_nys
  public :: Atm_nxx, Atm_nyy, Atm_extmax

  ! (Private)
  integer,  save :: Atm_iud = -1     ! unit index of input data file (-1: not opened)
  integer,  save :: Atm_iroot = 0    ! index of root PE
  integer,  save :: Atm_nchi = 3     ! # of truncation orders
  integer,  save :: Atm_ngsmax = 10  ! max factor to limit ngsray
  integer,  save :: Atm_nxsmax, Atm_nysmax, Atm_nzsmax
  integer,  save, allocatable :: Atm_ngsray(:)   !(nz+1) # of quadrature points
  real(R_), save, allocatable :: Atm_xgsray(:,:) !(ngsray,nz+1) quadrature abscissas
  real(R_), save, allocatable :: Atm_wgsray(:,:) !(ngsray,nz+1) quadrature weights
  real(R_), save, allocatable :: Atm_dgsray(:,:) !(2,nz+1) directional cosine threshold
  !// for Gaussian quadrature

contains

  !+
  ! Read in namelist variables (1) for initialization
  !-
  subroutine mcarAtm__user_init(iu, filepath)

    integer, intent(in) :: iu   ! file unit index for namelist input
    character(*), intent(in) :: filepath ! path name for I/O files
    integer :: ios

    ! Read in
    read (iu, nml=mcarAtm_nml_init, iostat=ios)
    call err_read(ios, iu, 'mcarAtm__user_init: Invalid namelist input for mcarAtm_nml_init.')
    if (len_trim(Atm_inpfile) > 0) Atm_inpfile = trim(filepath)//'/'//Atm_inpfile
    call mcarAtm__init_check()

    ! Defaults for core calculations
    Atm_ext1d(:,:) = 0.0_R_
    Atm_omg1d(:,:) = 1.0_R_
    Atm_apf1d(:,:) = 0.0_R_
    Atm_abs1d(:,:) = 0.0_R_

  end subroutine mcarAtm__user_init


  !+
  ! Check for initialization
  !-
  subroutine mcarAtm__init_check() 

    call check_iI('mcarAtm__init_check: Atm_iz3l', Atm_iz3l, 1, Atm_nz)
    call check_iI('mcarAtm__init_check: Atm_nz3',  Atm_nz3,  0, Atm_nz)

  end subroutine mcarAtm__init_check


  !+
  ! Initialize this module
  !-
  subroutine mcarAtm__init(irank, iroot, nchi) 

    integer, intent(in) :: irank, iroot, nchi
    integer :: nrec

    ! Assign
    Atm_nchi = nchi
    Atm_nppomax = Atm_np1d + Atm_np3d
    Atm_nxsmax = Atm_nx ! bug fix, Apr. 27, 2012
    Atm_nysmax = Atm_ny
    Atm_nzsmax = Atm_nz ! should be = nz
    Atm_iz3u = Atm_iz3l + Atm_nz3 - 1
    Atm_iroot = iroot

    ! Allocate
    if (allocated(Atm_mlay)) call mcarAtm__final()
    allocate (Atm_mlay(0:Atm_nz+1))
    allocate (Atm_xgrd(0:Atm_nx), Atm_ygrd(0:Atm_ny), Atm_zgrd(0:Atm_nz+1), Atm_zdep(Atm_nz+1))
    allocate (Atm_wkd(Atm_nkd))
    allocate (Atm_nppo(Atm_nz))
    allocate (Atm_tmpa3d(Atm_nx, Atm_ny, Atm_nz + 1))
    allocate (Atm_apfp3d(Atm_nx, Atm_ny, Atm_nz, Atm_nppomax))
    allocate (Atm_scap3d(Atm_nx, Atm_ny, Atm_nz, Atm_nppomax))
    allocate (Atm_scat3d(Atm_nx, Atm_ny, Atm_nz, Atm_nchi))
    allocate (Atm_abst3d(Atm_nx, Atm_ny, Atm_nz, Atm_nkd))
    allocate (Atm_extvar(Atm_nz))
    allocate (Atm_extmin(Atm_nz))
    allocate (Atm_extscl(Atm_nchi, Atm_nkd))
    allocate (Atm_jzstab(0 : Atm_nz + 1))
    allocate (Atm_xsup(0 : Atm_nxsmax, Atm_nzsmax))
    allocate (Atm_ysup(0 : Atm_nysmax, Atm_nzsmax))
    allocate (Atm_zsup(0 : Atm_nzsmax))
    allocate (Atm_jsup(0 : Atm_nzsmax + 1))
    allocate (Atm_jzupr(0 : Atm_nzsmax))
    allocate (Atm_nxs(Atm_nzsmax), Atm_nys(Atm_nzsmax))
    allocate (Atm_nxx(Atm_nzsmax), Atm_nyy(Atm_nzsmax))
    allocate (Atm_extmax(Atm_nxsmax, Atm_nysmax, Atm_nzsmax, Atm_nchi))
    allocate (Atm_ngsray(Atm_nz+1))
    allocate (Atm_xgsray(Atm_ngsmax * Atm_nqlay, Atm_nz+1))
    allocate (Atm_wgsray(Atm_ngsmax * Atm_nqlay, Atm_nz+1))
    allocate (Atm_dgsray(2, Atm_nz+1))

    ! Open input data file
    if (len_trim(Atm_inpfile) > 0 .and. irank == iroot) then ! root PE only
       Atm_iud = freeUnit_I(10)
       if (Atm_mfmt == 0) then
          call open_seq(Atm_iud, Atm_inpfile, 'old')
       else
          nrec = recordLen_I(1, Atm_nx * Atm_ny, 4) ! Apr. 27, 2012
          call open_dir(Atm_iud, Atm_inpfile, nrec, 'old')
       end if
    else
       Atm_iud = -1
    end if

  end subroutine mcarAtm__init


  !+
  ! Finalize this module
  !-
  subroutine mcarAtm__final() 

    if (allocated(Atm_mlay)) then
       deallocate (Atm_mlay, Atm_xgrd, Atm_ygrd, Atm_zgrd, Atm_zdep, Atm_wkd)
       deallocate (Atm_nppo, Atm_apfp3d, Atm_tmpa3d, Atm_extvar, Atm_extmin, Atm_extscl)
       deallocate (Atm_scap3d, Atm_scat3d, Atm_abst3d)
       deallocate (Atm_jzstab, Atm_xsup, Atm_ysup, Atm_zsup, Atm_jsup, Atm_jzupr)
       deallocate (Atm_nxs, Atm_nys, Atm_nxx, Atm_nyy, Atm_extmax)
       deallocate (Atm_ngsray, Atm_xgsray, Atm_wgsray, Atm_dgsray)
       if (Atm_iud >= 1) then ! root PE only
          close (Atm_iud)
          Atm_iud = -1
       end if
    end if

  end subroutine mcarAtm__final


  !+
  ! Tune a preset for technical parameters
  !-
  subroutine mcarAtm__tune_jobs(moptim) 

    integer, intent(in) :: moptim ! optimization flag [-2,3]

    if (moptim <= -2) then ! no preset
       return
    else if (moptim == -1) then ! no optimization
       call mcarAtm__set_tech(2.0_R_, 1.3_R_, 0.0_R_)
    else if (moptim == 0) then ! unbiased optimization
       call mcarAtm__set_tech(2.0_R_, 1.3_R_, 1.0_R_)
    else if (moptim == 1) then ! conservative
       call mcarAtm__set_tech(2.0_R_, 1.3_R_, 1.0_R_)
    else if (moptim == 2) then ! standard
       call mcarAtm__set_tech(2.0_R_, 1.3_R_, 2.0_R_)
    else ! quick-and-dirty
       call mcarAtm__set_tech(2.0_R_, 1.3_R_, 3.0_R_)
    end if

  end subroutine mcarAtm__tune_jobs

 
  !+
  ! Set technical parameters
  !-
  subroutine mcarAtm__set_tech(fsupg, fsupi, taumin) 

    real(R_), intent(in) :: fsupg, fsupi, taumin
    Atm_fsupg  = fsupg
    Atm_fsupi  = fsupi
    Atm_taumin = taumin

  end subroutine mcarAtm__set_tech


  !+
  ! Read in namelist variables (2) for core calculations
  !-
  subroutine mcarAtm__user_job(iu, idread0)

    integer, intent(in) :: iu   ! file unit index for namelist input
    integer, intent(in) :: idread0 ! default value for Atm_idread
    integer :: ios

    Atm_idread = idread0
    read (iu, nml=mcarAtm_nml_job, iostat=ios)
    call err_read(ios, iu, 'mcarAtm__user_job: Invalid namelist input for mcarAtm_nml_job.')
    call mcarAtm__core_check()

  end subroutine mcarAtm__user_job


  !+
  ! Check for core calculations
  !-
  subroutine mcarAtm__core_check() 

    integer :: iz
    do iz = 1, Atm_nz
       call check_iR('mcarAtm__core_check: Atm_zgrd0', Atm_zgrd0(iz+1), Atm_zgrd0(iz))
    end do ! zgrd should be increasing
    call check_iR('mcarAtm__core_check: Atm_taumin',  Atm_taumin,  0.0_R_,  1000.0_R_)

  end subroutine mcarAtm__core_check


  !+
  ! Set mininum extinction coefficient for collision-forcing
  !-
  subroutine mcarAtm__set_extMin(extmin) 

    real(R_), intent(in) :: extmin(:)
    Atm_extmin(1:Atm_nz) = extmin(1:Atm_nz)

  end subroutine mcarAtm__set_extMin

    
  !+
  ! Prepare for core calculations
  !-
  subroutine mcarAtm__prep_job(irank, mbswap, npf, tmpmin, tmpmax)

    integer,  intent(in) :: irank, mbswap
    integer,  intent(in) :: npf
    real(R_), intent(out) :: tmpmin, tmpmax
    integer :: iz
    real(R_) :: s

    ! Initialize
    Atm_xmax = Atm_dx * Atm_nx
    Atm_ymax = Atm_dy * Atm_ny
    Atm_facx = Atm_nx / Atm_xmax * 0.9999995_R_
    Atm_facy = Atm_ny / Atm_ymax * 0.9999995_R_
    Atm_xgrd(0:Atm_nx) = gridValues_1R(Atm_nx + 1, 0.0_R_, Atm_xmax)
    Atm_ygrd(0:Atm_ny) = gridValues_1R(Atm_ny + 1, 0.0_R_, Atm_ymax)
    Atm_zgrd(0:Atm_nz) = Atm_zgrd0(1 : Atm_nz + 1)
    Atm_zgrd(Atm_nz + 1) = Atm_zgrd0(Atm_nz + 1) + (Atm_zgrd0(Atm_nz + 1) - Atm_zgrd0(1))
    Atm_zdep(1 : Atm_nz + 1) = Atm_zgrd(1 : Atm_nz + 1) - Atm_zgrd(0 : Atm_nz)
    Atm_wkd(1:Atm_nkd) = Atm_wkd0(1:Atm_nkd)
    s = sum(Atm_wkd(1:Atm_nkd))
    Atm_wkd(:) = Atm_wkd(:) / s ! normalized CKD weights

    ! Import new dataset
    if (Atm_idread >= 1) then ! only when new dataset is requested
       call mcarAtm__prep_optInit()
       if (irank == Atm_iroot .and. Atm_iud >= 1) call mcarAtm__read(mbswap) ! root PE only
#if UseMPI == 1
       call mcarAtm__MPI_bcast_read(Atm_iroot)
#endif
       call mcarAtm__prep_opt3D(npf)
       call mcarAtm__prep_optTot()
       call mcarAtm__prep_temp(tmpmin, tmpmax)
       !// These are not parallelized. This could be improved in the future.
    end if

    ! Quadrature variables
    Atm_ngsray(1:Atm_nz) = Atm_nqlay * min(Atm_ngsmax, int(1.0_R_ + Atm_extvar(1:Atm_nz)))
    Atm_ngsray(Atm_nz+1) = Atm_nqlay
    do iz = 1, Atm_nz + 1
       call gaussLegen(Atm_ngsray(iz), Atm_xgsray(:, iz), Atm_wgsray(:, iz))
    end do
    Atm_xgsray(:, :) = 0.5_R_ * (1.0_R_ + Atm_xgsray(:, :)) ! [0,1]
    Atm_dgsray(1, :) = 2.0_R_ * Atm_ngsray(:) * Atm_dx / Atm_zdep(:) ! tanQx = dx/dz
    Atm_dgsray(2, :) = 2.0_R_ * Atm_ngsray(:) * Atm_dy / Atm_zdep(:) ! tanQy = dy/dz

  end subroutine mcarAtm__prep_job


  !+
  ! Initialize all optical properties
  !-
  subroutine mcarAtm__prep_optInit() 

    Atm_scap3d(:,:,:,:) =  0.0_R_
    Atm_abst3d(:,:,:,:) =  0.0_R_
    Atm_apfp3d(:,:,:,:) = -1.0_R_ ! default, Rayleigh
    Atm_tmpa3d(:,:,:) = 0.0_R_ ! temperature perturbation

  end subroutine mcarAtm__prep_optInit



  !+
  ! Check and diagnose optical properties
  !-
  subroutine mcarAtm__check_opt3DStats() 

    integer :: n, iz, ikd, ippo
    real(R_) :: aave, amin, amax

    n = Atm_nx * Atm_ny
    write (*,*) '# Absorption layer optical depth (Atm_abst3d*dZ)'
    write (*,*) '## iz, ikd, mean, min, max'
    do ikd = 1, Atm_nkd
       do iz = 1, Atm_nz
          aave = sum(Atm_abst3d(:,:,iz,ikd) * Atm_zdep(iz)) / n
          amin = minval(Atm_abst3d(:,:,iz,ikd)) * Atm_zdep(iz)
          amax = minval(Atm_abst3d(:,:,iz,ikd)) * Atm_zdep(iz)
          write (*,'(2i7, 3es13.5)') iz, ikd, aave, amin, amax
       end do
    end do

    write (*,*) '# Scattering layer optical depth (Atm_scap3d*dZ)'
    write (*,*) '## iz, ippo, mean, min, max'
    do ippo = 1, Atm_nppomax
       do iz = 1, Atm_nz
          aave = sum(Atm_scap3d(:,:,iz,ippo) * Atm_zdep(iz)) / n
          amin = minval(Atm_scap3d(:,:,iz,ippo)) * Atm_zdep(iz)
          amax = minval(Atm_scap3d(:,:,iz,ippo)) * Atm_zdep(iz)
          write (*,'(2i7, 3es13.5)') iz, ippo, aave, amin, amax
       end do
    end do

    write (*,*) '# Scattering phase function idex (Atm_apfp3d)'
    write (*,*) '## iz, ippo, mean, min, max'
    do ippo = 1, Atm_nppomax
       do iz = 1, Atm_nz
          aave = sum(Atm_apfp3d(:,:,iz,ippo)) / n
          amin = minval(Atm_apfp3d(:,:,iz,ippo))
          amax = minval(Atm_apfp3d(:,:,iz,ippo))
          write (*,'(2i7, 3es13.5)') iz, ippo, aave, amin, amax
       end do
    end do

    write (*,*) '# Atmospheric temperature (Atm_tmpa3d)'
    write (*,*) '## iz, mean, min, max'
    do iz = 1, Atm_nz
       aave = sum(Atm_tmpa3d(:,:,iz)) / n
       amin = minval(Atm_tmpa3d(:,:,iz))
       amax = minval(Atm_tmpa3d(:,:,iz))
       write (*,'(1i7, 3es13.5)') iz, aave, amin, amax
    end do

  end subroutine mcarAtm__check_opt3DStats


  !+
  ! Import data from file
  !// Note this procedure touches only layers iz = [iz3l:iz3u].
  !-
  subroutine mcarAtm__read(mbswap) 

    integer,  intent(in) :: mbswap
    integer :: ip, ikd, iz, iz3
    real(R_), allocatable :: omgp3d(:,:,:,:)

    ! Initialize
    allocate (omgp3d(Atm_nx, Atm_ny, Atm_nz3, Atm_np3d))

    ! Open & read in
    if (Atm_mfmt == 0) then ! text
       call mcarUtl__txtAtm_read(Atm_iud, Atm_nx, Atm_ny, Atm_nz3, Atm_nkd, Atm_np3d, &
            & Atm_mtprof, Atm_tmpa3d(:,:,Atm_iz3l:), Atm_abst3d(:,:,Atm_iz3l:,:), &
            & Atm_apfp3d(:,:,Atm_iz3l:,:), Atm_scap3d(:,:,Atm_iz3l:,:), omgp3d)
    else ! binary
       call mcarUtl__binAtm_read(Atm_iud, Atm_nx, Atm_ny, Atm_nz3, Atm_nkd, Atm_np3d, &
            & Atm_mtprof, Atm_tmpa3d(:,:,Atm_iz3l:), Atm_abst3d(:,:,Atm_iz3l:,:), &
            & Atm_apfp3d(:,:,Atm_iz3l:,:), Atm_scap3d(:,:,Atm_iz3l:,:), omgp3d, mbswap, Atm_idread)
    end if
    !// Atm_scap3d() represents tentatively the extinction coefficients, here.
    !// Atm_abst3d() represents tentatively the gaseous absorption coefficients, here.

    ! Scale gaseous absorption coefficients
    if (Atm_fabs3d /= 1.0_R_) Atm_abst3d(:,:, Atm_iz3l:Atm_iz3u, :) = &
         & Atm_abst3d(:,:, Atm_iz3l:Atm_iz3u, :) * Atm_fabs3d

    ! Bext & Omega --> Bsca & Babs
    do ip = 1, Atm_np3d
       do iz = Atm_iz3l, Atm_iz3u
          iz3 = iz - Atm_iz3l + 1
          do ikd = 1, Atm_nkd
             Atm_abst3d(:, :, iz, ikd) = Atm_abst3d(:, :, iz, ikd) + Atm_scap3d(:, :, iz, ip) &
                  & * (1.0_R_ - omgp3d(:, :, iz3, ip)) * Atm_fext3d(ip)
          end do
          Atm_scap3d(:, :, iz, ip) = Atm_scap3d(:, :, iz, ip) * omgp3d(:, :, iz3, ip) * Atm_fext3d(ip)
          !// Now, Atm_abst3d() represents total absorption coefficients.
          !// Now, Atm_scap3d() represents scattering coefficients.
       end do
    end do
    deallocate (omgp3d)

  end subroutine mcarAtm__read


#if UseMPI == 1
  !+
  ! MPI broadcast data read in mcarAtm__read
  !-
  subroutine mcarAtm__MPI_bcast_read(iroot)

    include 'inc_mpi.f90'
    integer, intent(in) :: iroot ! root PE index
    integer :: ierr, nxy, iz, ip, ikd

    nxy = Atm_nx * Atm_ny
    do iz = Atm_iz3l, Atm_iz3u
       call MPI_Bcast(Atm_tmpa3d(:, :, iz), nxy, MPI_R_, iroot, MPI_COMM_WORLD, ierr)
       do ikd = 1, Atm_nkd
          call MPI_Bcast(Atm_abst3d(:, :, iz, ikd), nxy, MPI_R_, iroot, MPI_COMM_WORLD, ierr)
       end do
       do ip = 1, Atm_np3d
          call MPI_Bcast(Atm_scap3d(:, :, iz, ip), nxy, MPI_R_, iroot, MPI_COMM_WORLD, ierr)
          call MPI_Bcast(Atm_apfp3d(:, :, iz, ip), nxy, MPI_R_, iroot, MPI_COMM_WORLD, ierr)
       end do
    end do
    iz = Atm_iz3u + 1
    call MPI_Bcast(Atm_tmpa3d(:, :, iz), nxy, MPI_R_, iroot, MPI_COMM_WORLD, ierr)

  end subroutine mcarAtm__MPI_bcast_read
#endif


  !+
  ! Synthesize an entire 3-D field of optical properties, from 3-D and 1-D distributions
  !-
  subroutine mcarAtm__prep_opt3D(npf)

    integer,  intent(in) :: npf
    real(R_) :: tamin, tsmin, amin, smin, a, s, e, apfmax
    integer  :: ikd, ippo, iz, ip1d, ip3d
    real(R_), parameter :: TAUSMIN = 1.0e-8_R_, TAUAMIN = 1.0e-8_R_

    ! 3-D scattering coefficients
    tsmin = TAUSMIN / (Atm_nz * Atm_nppomax) ! min scattering optical thickness per layer
    do ip3d = 1, Atm_np3d
       do iz = Atm_iz3l, Atm_iz3u
          smin = tsmin / Atm_zdep(iz)
          Atm_scap3d(:, :, iz, ip3d) = max(smin, Atm_scap3d(:, :, iz, ip3d))
       end do
    end do

    ! Add 1-D distributions of particle properties
    apfmax = max(real(npf,R_), 0.999995_R_)
    do iz = 1, Atm_nz
       smin = tsmin / Atm_zdep(iz)
       if (iz >= Atm_iz3l .and. iz <= Atm_iz3u) then
          ippo = Atm_np3d
       else
          ippo = 0
          Atm_abst3d(:, :, iz, :) = 0.0_R_
       end if
       do ip1d = 1, Atm_np1d
          e = Atm_ext1d(iz, ip1d) * Atm_fext1d(ip1d)
          a = e * (1.0_R_ - Atm_omg1d(iz, ip1d))
          s = e * Atm_omg1d(iz, ip1d)
          Atm_abst3d(:, :, iz, :) =  Atm_abst3d(:, :, iz, :) + a
          if (s >= smin .or. ippo == 0) then ! register active components only
             ippo = ippo + 1
             Atm_apfp3d(:, :, iz, ippo) = Atm_apf1d(iz, ip1d)
             Atm_scap3d(:, :, iz, ippo) = s
          end if
       end do
       Atm_nppo(iz) = ippo ! should be > 0 (i.e. at least one component should be registered)
    end do
    Atm_apfp3d(:,:,:,:) = max(-2.99_R_, min(apfmax, Atm_apfp3d(:,:,:,:)))

    !print *, 'nppo=', Atm_nppo(:)
    !print *, 'apfp3d=', Atm_apfp3d(:,:,:,:) ! DEBUG

    ! Add 1-D distribution of gaseous absorption coefficients
    tamin = TAUAMIN / Atm_nz
    do ikd = 1, Atm_nkd
       do iz = 1, Atm_nz
          amin = tamin / Atm_zdep(iz)
          a = Atm_abs1d(iz, ikd) * Atm_fabs1d
          Atm_abst3d(:, :, iz, ikd) = max(amin, a + Atm_abst3d(:, :, iz, ikd))
       end do
    end do

  end subroutine mcarAtm__prep_opt3D


  !+
  ! Setup total optical properties & super-cells
  !// This procedure uses a mcarSca procedure, so that tables in mcarSca should be
  !   initialized before using this procedure.
  !-
  subroutine mcarAtm__prep_optTot()

    real(R_) :: scaa(Atm_nz, Atm_nchi), absa(Atm_nz, Atm_nkd), extxy(Atm_nx, Atm_ny) !(AUTO)
    real(R_) :: tau, ave, var
    integer  :: ichi, ikd, ippo, iy, iz

    ! Total scaled scattering coefficients
    Atm_scat3d(:, :, :, :) = 0.0_R_
    do ichi = 1, Atm_nchi
       do iz = 1, Atm_nz
          do ippo = 1, Atm_nppo(iz)
             do iy = 1, Atm_ny
                Atm_scat3d(:, iy, iz, ichi) = Atm_scat3d(:, iy, iz, ichi) + Atm_scap3d(:, iy, iz, ippo) &
                     & * mcarSca__fracDelC_1R(Atm_apfp3d(:, iy, iz, ippo), ichi)
             end do
          end do
       end do
    end do

    ! Kinds of layer
    ichi = 1
    Atm_mlay(0) = 0 ! surface
    Atm_mlay(Atm_nz + 1) = 4 ! space
    do iz = 1, Atm_nz
       if (iz >= Atm_iz3l .and. iz <= Atm_iz3u) then
          Atm_mlay(iz) = 3 ! 3-D
       else if (Atm_scat3d(1, 1, iz, ichi) < REPS_ * minval(Atm_abst3d(1, 1, iz, :))) then
          Atm_mlay(iz) = 1 ! 1-D, absorbing
       else
          Atm_mlay(iz) = 2 ! 1-D, with scattering
       end if
    end do

    ! Average scattering coefficients
    do ichi = 1, Atm_nchi
       scaa(:, ichi) = Atm_scat3d(1, 1, :, ichi) ! tentative
       do iz = Atm_iz3l, Atm_iz3u
          scaa(iz, ichi) = sum(Atm_scat3d(:, :, iz, ichi)) / (Atm_nx * Atm_ny) ! overwrite
          !print *, iz, ichi, extp(iz, ichi), Atm_mlay(iz)
       end do
    end do

    ! Average absorption coefficients
    do ikd = 1, Atm_nkd
       absa(:, ikd) = Atm_abst3d(1, 1, :, ikd) ! tentative
       do iz = Atm_iz3l, Atm_iz3u
          absa(iz, ikd) = sum(Atm_abst3d(:, :, iz, ikd)) / (Atm_nx * Atm_ny) ! overwrite
          !print *, iz, ikd, absg(iz, ikd), Atm_mlay(iz)
       end do
    end do

    ! Scaling parameters for collision-forcing
    Atm_extscl(:,:) = 1.0_R_ ! tentative
    do ikd = 1, Atm_nkd
       do ichi = 1, Atm_nchi
          tau = 0.0_R_
          do iz = 1, Atm_nz
             tau = tau + (scaa(iz, ichi) + absa(iz, ikd)) * Atm_zdep(iz)
          end do
          if (tau < Atm_taumin) Atm_extscl(ichi, ikd) = Atm_taumin / tau
          !print *, ichi, ikd, tau, taumin, Atm_extscl(ichi, ikd)
       end do
    end do

    ! Profile of normazlied variance of extinction
    ichi = 1
    ikd  = 1
    Atm_extvar(:) = 0.0_R_
    do iz = Atm_iz3l, Atm_iz3u
       extxy(:,:) = Atm_scat3d(:, :, iz, ichi) + Atm_abst3d(:, :, iz, ikd)
       ave = sum(extxy) / (Atm_nx * Atm_ny)
       var = sum((extxy(:,:) - ave)**2) / (Atm_nx * Atm_ny)
       Atm_extvar(iz) = var / ave**2
       !print *, iz, Atm_extvar(iz), ave, var
    end do

    ! Generate super-cells
    call mcarAtm__prep_supCell(scaa, absa)

  end subroutine mcarAtm__prep_optTot


  !+
  ! Generate super-cells
  !-
  subroutine mcarAtm__prep_supCell(scaa, absa)

    real(R_), intent(in) :: scaa(:,:)
    real(R_), intent(in) :: absa(:,:)
    real(R_)  :: ewrk(Atm_nz), wwrk(Atm_nx, Atm_ny), awrk(Atm_nx, Atm_ny) !(AUTO)
    real(R_)  :: dxs, dys, dzs, w, a
    integer   :: ichi, ikd, ixs, iys, iz, izs

    ! Initializations
    ikd = 1
    ichi = 1
    ewrk(:) = scaa(:, ichi) + absa(:, ikd)

    ! Super-layer generation
    Atm_jsup(0) = 0
    Atm_jzupr(0) = 0
    Atm_zsup(0) = Atm_zgrd(0)
    iz = 0
    do izs = 1, Atm_nzsmax
       iz = iz + 1
       Atm_jsup(izs) = Atm_mlay(iz)
       w = ewrk(iz)
       a = ewrk(iz) * Atm_zdep(iz)
       do
          if (Atm_jsup(izs) /= Atm_mlay(iz + 1)) exit
          w = max(w, ewrk(iz + 1))
          a = a + ewrk(iz + 1) * Atm_zdep(iz + 1)
          dzs = Atm_zgrd(iz + 1) - Atm_zsup(izs - 1)
          if (w * Atm_fsupg * dzs > 1.0_R_ .and. w * dzs > Atm_fsupi * a) exit
          if (iz + 1 - Atm_jzupr(izs - 1) > Atm_nzzsup) exit
          iz = iz + 1
       end do
       Atm_jzupr(izs) = iz
       Atm_zsup(izs) = Atm_zgrd(iz)
       !write (*,*) izs, iz, Atm_zsup(izs)
       if (iz >= Atm_nz) exit
    end do
    Atm_nzs = izs
    Atm_jsup(Atm_nzs + 1) = Atm_mlay(Atm_nz + 1) ! free space above TOA

    ! Inverse reference table
    Atm_jzstab(0) = 0
    Atm_jzstab(Atm_nz + 1) = Atm_nzs + 1
    do izs = 1, Atm_nzs
       do iz = Atm_jzupr(izs - 1) + 1, Atm_jzupr(izs)
          Atm_jzstab(iz) = izs
       end do
    end do

    ! Super-voxel generation
    do izs = 1, Atm_nzs
       if (Atm_jsup(izs) <= 2) then ! 1-D super-layer
          Atm_nxx(izs) = 1
          Atm_nyy(izs) = 1
          Atm_nxs(izs) = 1
          Atm_nys(izs) = 1
       else                   ! 3-D super-layer
          wwrk(:,:) = 0.0_R_
          do iz = Atm_jzupr(izs - 1) + 1, Atm_jzupr(izs)
             wwrk(:,:) = max(wwrk(:,:), Atm_scat3d(:, :, iz, ichi) + Atm_abst3d(:, :, iz, ikd))
          end do
          awrk(:,:) = wwrk(:,:)
          call mcarAtm__set_supVox(wwrk, awrk, Atm_nx, Atm_ny, Atm_xmax, Atm_ymax, &
               & Atm_fsupg, Atm_fsupi, Atm_nxx(izs), Atm_nyy(izs), Atm_nxs(izs), Atm_nys(izs))
       end if
       dxs = Atm_xmax / real(Atm_nx, R_) * real(Atm_nxx(izs), R_)
       dys = Atm_ymax / real(Atm_ny, R_) * real(Atm_nyy(izs), R_)
       do ixs = 0, Atm_nxs(izs) - 1
          Atm_xsup(ixs, izs) = dxs * real(ixs, R_)
       end do
       Atm_xsup(Atm_nxs(izs), izs) = Atm_xmax
       do iys = 0, Atm_nys(izs) - 1
          Atm_ysup(iys, izs) = dys * real(iys, R_)
       end do
       Atm_ysup(Atm_nys(izs), izs) = Atm_ymax
       !write (*,*) izs, Atm_jzupr(izs), Atm_nxx(izs), Atm_nyy(izs), Atm_nxs(izs), Atm_nys(izs)
    end do

  end subroutine mcarAtm__prep_supCell


  !+
  ! Set temperatures and get min & max
  !-
  subroutine mcarAtm__prep_temp(tmpmin, tmpmax)

    real(R_), intent(out) :: tmpmin, tmpmax
    real(R_) :: zmid(Atm_nz), tmid(Atm_nz), alev(Atm_nz + 1), tmin, tmax
    integer :: ix, iy, iz, izm(Atm_nz + 1)

    ! Merge 1-D and 3-D data
    if (Atm_iz3u < Atm_iz3l .or. Atm_nz3 <= 0) then ! all 1-D layers
       do iz = 1, Atm_nz + 1
          Atm_tmpa3d(:,:, iz) = Atm_tmp1d(iz)
       end do
    else ! 1-D & 3-D layers
       do iz = 1, Atm_iz3l - 1
          Atm_tmpa3d(:,:, iz) = Atm_tmp1d(iz)
       end do
       do iz = Atm_iz3l, Atm_iz3u + Atm_mtprof
          Atm_tmpa3d(:,:, iz) = Atm_tmpa3d(:,:, iz) + Atm_tmp1d(iz)
       end do
       do iz = Atm_iz3u + 1 + Atm_mtprof, Atm_nz + Atm_mtprof
          Atm_tmpa3d(:,:, iz) = Atm_tmp1d(iz)
       end do
    end if

    ! Interpolation for (nz+1) boundaries (if needed)
    if (Atm_mtprof == 0) then
       if (Atm_nz == 1) then
          Atm_tmpa3d(:,:,2) = Atm_tmpa3d(:,:,1)
       else
          zmid(:) = (Atm_zgrd(0:Atm_nz-1) + Atm_zgrd(1:Atm_nz)) * 0.5_R_ ! layer mid level
          if (Atm_nz <= 3) then ! linear interpolation
             do iz = 1, Atm_nz + 1
                izm(iz) = max(1, min(Atm_nz - 1, iz - 1))
                alev(iz) = (Atm_zgrd(iz - 1) - zmid(izm(iz))) / (zmid(izm(iz) + 1) - zmid(izm(iz)))
             end do
          else ! cubic interpolation
             do iz = 1, Atm_nz + 1
                izm(iz) = max(2, min(Atm_nz - 2, iz - 1)) - 1
             end do
          end if
          do iy = 1, Atm_ny
             do ix = 1, Atm_nx
                tmid(1:Atm_nz) = Atm_tmpa3d(ix, iy, 1:Atm_nz) ! save original data
                tmin = minval(tmid) * 0.5_R_
                tmax = maxval(tmid) * 2.0_R_
                do iz = 1, Atm_nz + 1
                   if (Atm_nz <= 3) then ! linear
                      Atm_tmpa3d(ix, iy, iz) = (1.0_R_ - alev(iz)) * tmid(izm(iz)) &
                           & + alev(iz) * tmid(izm(iz) + 1)
                   else ! cubic
                      Atm_tmpa3d(ix, iy, iz) = intp_cubic_R(Atm_zgrd(iz - 1), &
                           & zmid(izm(iz) : izm(iz)+3), tmid(izm(iz) : izm(iz)+3))
                   end if
                   Atm_tmpa3d(ix, iy, iz) = max(tmin, min(tmax, Atm_tmpa3d(ix, iy, iz)))
                end do
             end do
          end do
       end if
    end if

    ! Min & max
    tmpmin = minval(Atm_tmpa3d)
    tmpmax = maxval(Atm_tmpa3d)

  end subroutine mcarAtm__prep_temp


  !+
  ! Prepare the maximum extinction coefficients in super-grids
  !-
  subroutine mcarAtm__prep_extMax(ikd)

    integer,  intent(in) :: ikd
    real(R_), parameter :: TAUEPS = 1.0e-7_R_ ! very small column optical tchikness
    real(R_)  :: e, eeps, emin, dtaueps
    integer   :: ichi, ixmax, ixmin, ixs, iymax, iymin, iys, izmax, izmin, izs

    dtaueps = TAUEPS / Atm_nz ! optical thickness threshold per layer
    !// This is used to detect clear super-voxels with no active component.

    ! Loop for super-layers
    do izs = 1, Atm_nzs
       izmin = Atm_jzupr(izs - 1) + 1
       izmax = Atm_jzupr(izs)
       eeps = dtaueps / (Atm_zgrd(izmax) - Atm_zgrd(izmin - 1))
       emin = maxval(Atm_extmin(izmin:izmax))

       ! 1-D super-layers
       if (Atm_jsup(izs) <= 2) then
          do ichi = 1, Atm_nchi
             e = maxval(Atm_scat3d(1, 1, izmin:izmax, ichi) + Atm_abst3d(1, 1, izmin:izmax, ikd))
             if (e < eeps) then ! completely clear
                Atm_extmax(1, 1, izs, ichi) = 0.0_R_
             else ! with at least one active cell (so that CF is used)
                Atm_extmax(1, 1, izs, ichi) = Atm_extscl(ichi, ikd) * max(e, emin)
             end if
          end do

          ! 3-D super-layers
       else
          do iys = 1, Atm_nys(izs)
             iymin = Atm_nyy(izs) * (iys - 1) + 1
             iymax = Atm_nyy(izs) * iys
             if (iys == Atm_nys(izs)) iymax = Atm_ny
             do ixs = 1, Atm_nxs(izs)
                ixmin = Atm_nxx(izs) * (ixs - 1) + 1
                ixmax = Atm_nxx(izs) * ixs
                if (ixs == Atm_nxs(izs)) ixmax = Atm_nx
                do ichi = 1, Atm_nchi ! for all truncation regimes
                   e = maxval(Atm_scat3d(ixmin:ixmax, iymin:iymax, izmin:izmax, ichi) &
                        &   + Atm_abst3d(ixmin:ixmax, iymin:iymax, izmin:izmax, ikd))
                   if (e < eeps) then ! completely clear
                      Atm_extmax(ixs, iys, izs, ichi) = 0.0_R_
                   else ! with at least one active cell (so that CF is used)
                      Atm_extmax(ixs, iys, izs, ichi) = Atm_extscl(ichi, ikd) * max(e, emin)
                   end if
                end do
             end do
          end do
       end if
    end do

  end subroutine mcarAtm__prep_extMax


  !+
  ! PDF of angular distribution (/steradian)
  !// Note: PDF (1/steradian) for angular (re)distribution
  !    PDF = 1/(4*pi) * P"(Q)
  !-
  function mcarAtm__angDistr_R(dir1, dir2, apf, ichi) result(res)

    real(R_), intent(in) :: dir1(:)
    real(R_), intent(in) :: dir2(:)
    real(R_), intent(in) :: apf
    integer,  intent(in) :: ichi
    real(R_)  :: res
    real(R_),  parameter :: API4 = 0.25_R_ / PI_
    real(R_),  parameter :: PHSFMAX = 1.0e+5_R_
    real(R_)  :: cosq, q, rpf
    integer   :: ipf
    
    ! Isotropic
    ipf = int(apf)
    if (ipf <= -2) then
       res = API4
       
       ! Anisotropic
    else
       !call geo3D_twoUVec(dir1, dir2, cosq, sinq)
       cosq = geo3D_scaProd_R(dir1, dir2)
       cosq = max(-1.0_R_, min(1.0_R_, cosq))
       if (ipf == -1) then ! Rayleigh
          res = API4 * rayl_scatPF_R(cosq)
       else if (ipf == 0) then ! H-G
          res = API4 * scatPF_HG_R(apf, cosq)
       else ! tabulated phase function
          !q = atan2(sinq, cosq) ! atan2 is better than acos when cosq ~= +-1
          q = acos(cosq)
          rpf = apf - ipf ! index interpolation factor
          res = API4 * min(PHSFMAX, mcarSca__phaseFunc_R(q, ipf, rpf, ichi))
       end if
    end if

  end function mcarAtm__angDistr_R


  !+
  ! Escape from the event point to the detector point
  !-
  subroutine mcarAtm__rt_escape(tau, taumax, plen, loc, dir, itr, ix, iy, iz, ikd, ichi, zdest, mv3D)

    real(R_), intent(inout) :: plen(0:)
    integer,  intent(in) :: ichi
    integer,  intent(in) :: ikd
    real(R_), intent(inout) :: tau
    real(R_), intent(in) :: taumax
    real(R_), intent(inout) :: loc(:)
    real(R_), intent(in) :: dir(:)
    integer,  intent(in) :: itr(:)
    real(R_), intent(in) :: zdest
    logical,  intent(in) :: mv3D
    integer,  intent(inout) :: ix
    integer,  intent(inout) :: iy
    integer,  intent(inout) :: iz
    !// Important note: iz should be in [0,nz+1]. Otherwise, a fatal error occurs.
    real(R_)  :: extt, path
    integer   :: icase, iq, icell, ncell
    real(R_)  :: locnew(3), dpath, auz

    ! No possibility to escape
    if ((zdest - loc(3)) * dir(3) <= 0.0_R_) then
       tau = RLRG_
       return
    end if
    auz = 1.0_R_ / dir(3)

    ! From TOA/BOA
    if (iz <= 0 .or. iz >= Atm_nz + 1) then ! would be local source or solar source
       if ((iz <= 0 .and. zdest > Atm_zgrd(0)) .or. &
            & (iz >= Atm_nz + 1 .and. zdest < Atm_zgrd(Atm_nz))) then ! entering the atmosphere
          path = (Atm_zgrd(iz - 1 + itr(3)) - loc(3)) * auz
          plen(iz) = plen(iz) + path
          call mcarAtm__rt_bound1D(itr(3), iz, loc(3)) ! to the TOA/BOA boundary
          if (mv3D) call mcarAtm__rt_renewLocXY(dir, path, loc) 
          ix = 0
       end if
    end if

    ! Loop for atmospheric layers
    loop_domain: do

       !  1-D layer
       if (Atm_mlay(iz) == 1 .or. Atm_mlay(iz) == 2) then
          loop_1d: do ! loop for 1-D layers
             extt = Atm_abst3d(1, 1, iz, ikd) + Atm_scat3d(1, 1, iz, ichi)
             locnew(3) = Atm_zgrd(iz - 1 + itr(3))
             if ((zdest - locnew(3)) * dir(3) <= 0.0_R_) exit loop_domain ! destination
             path = (locnew(3) - loc(3)) * auz
             tau = tau + path * extt
             if (tau > taumax) return
             plen(iz) = plen(iz) + path
             call mcarAtm__rt_bound1D(itr(3), iz, loc(3))
             if (mv3D) call mcarAtm__rt_renewLocXY(dir, path, loc) 
             if (Atm_mlay(iz) < 1 .or. Atm_mlay(iz) > 2) exit loop_1d
          end do loop_1d
          ix = 0

          !  3-D layer
       else if (Atm_mlay(iz) == 3) then
          if (ix <= 0) then
             ix = gridIdx_unif0_I(Atm_xgrd, loc(1), 1, Atm_nx)
             iy = gridIdx_unif0_I(Atm_ygrd, loc(2), 1, Atm_ny)
          else
             call gridIdx_fix0(Atm_xgrd, loc(1), Atm_nx, ix)
             call gridIdx_fix0(Atm_ygrd, loc(2), Atm_ny, iy)
          end if

          loop_3d: do ! loop for 3-D layers
             if (Atm_mlay(iz) /= 3) exit loop_3d

             ! Vertical transfer only
             if (.not. mv3D) then
                extt = Atm_abst3d(ix, iy, iz, ikd) + Atm_scat3d(ix, iy, iz, ichi)
                locnew(3) = Atm_zgrd(iz - 1 + itr(3))
                if ((zdest - locnew(3)) * dir(3) <= 0.0_R_) exit loop_domain ! destination
                path = (locnew(3) - loc(3)) * auz
                tau = tau + path * extt
                if (tau > taumax) return
                plen(iz) = plen(iz) + path
                call mcarAtm__rt_bound1D(itr(3), iz, loc(3))

                ! 3-D transfer
             else
                if (   abs(dir(1)) > Atm_dgsray(1, iz) * abs(dir(3)) .or. &
                     & abs(dir(2)) > Atm_dgsray(2, iz) * abs(dir(3))) then
                   ncell = Atm_ngsray(iz)
                else
                   ncell = IMAX_
                end if

                ! Cell-by-cell integration (exact)
                do icell = 1, ncell ! loop for cells
                   extt = Atm_abst3d(ix, iy, iz, ikd) + Atm_scat3d(ix, iy, iz, ichi)
                   call mcarAtm__rt_path3D(loc, ix, iy, iz, dir, itr, icase, path)
                   locnew(3) = loc(3) + path * dir(3)
                   if ((zdest - locnew(3)) * dir(3) <= 0.0_R_) exit loop_domain ! destination
                   tau = tau + path * extt
                   if (tau > taumax) return
                   plen(iz) = plen(iz) + path
                   call mcarAtm__rt_bound3D(icase, path, dir, itr, loc, ix, iy, iz)
                   if (icase == 3) cycle loop_3d ! move to the next layer
                end do

                ! Gaussian quadrature (for nearly-horizontal transfer)
                locnew(3) = Atm_zgrd(iz - 1 + itr(3))
                if ((zdest - locnew(3)) * dir(3) <= 0.0_R_) locnew(3) = zdest ! destination
                path = (locnew(3) - loc(3)) * auz
                do iq = 1, Atm_ngsray(iz) ! loop for quadrature points
                   dpath = path * Atm_xgsray(iq, iz)
                   locnew(1:2) = loc(1:2)
                   call mcarAtm__rt_renewLocXY(dir, dpath, locnew) 
                   ix = min(Atm_nx, int(locnew(1) * Atm_facx) + 1)
                   iy = min(Atm_ny, int(locnew(2) * Atm_facy) + 1)
                   extt = Atm_abst3d(ix, iy, iz, ikd) + Atm_scat3d(ix, iy, iz, ichi)
                   tau = tau + extt * path * Atm_wgsray(iq, iz)
                   if (tau > taumax) return
                end do
                if (locnew(3) == zdest) then
                   extt = 0.0_R_
                   exit loop_domain
                end if
                plen(iz) = plen(iz) + path
                call mcarAtm__rt_bound1D(itr(3), iz, loc(3))
                call mcarAtm__rt_renewLocXY(dir, path, loc)
                ix = min(Atm_nx, int(loc(1) * Atm_facx) + 1)
                iy = min(Atm_ny, int(loc(2) * Atm_facy) + 1)
             end if
          end do loop_3d

          ! Exit from the atmosphere (iz = 0 or nz + 1)
       else !// The destination is outside of the atmosphere.
          ix = 0
          extt = 0.0_R_
          exit loop_domain
       end if
    end do loop_domain

    ! The last cell
    path = (zdest - loc(3)) * auz !/ dir(3) ! should be > 0
    if (mv3D) call mcarAtm__rt_renewLocXY(dir, path, loc)
    loc(3) = zdest
    tau = tau + extt * path
    plen(iz) = plen(iz) + path

  end subroutine mcarAtm__rt_escape


  !+
  ! Motion to the 1-D layer boundary : modify Z
  !-
  subroutine mcarAtm__rt_bound1D(itr, iz, z) 

    integer,  intent(in) :: itr ! 0 for downward, 1 for upward
    integer,  intent(inout) :: iz
    real(R_), intent(out) :: z
    z = Atm_zgrd(iz - 1 + itr)
    if (itr == 0) then ! downward
       iz = iz - 1
    else               ! upward
       iz = iz + 1
    end if

  end subroutine mcarAtm__rt_bound1D


  !+
  ! Renew the X & Y location
  !-
  subroutine mcarAtm__rt_renewLocXY(dir, path, loc) 

    real(R_), intent(in) :: dir(:)
    real(R_), intent(in) :: path
    real(R_), intent(inout) :: loc(:)
    loc(1:2) = loc(1:2) + path * dir(1:2) ! move
    if (loc(1) < 0.0_R_ .or. loc(1) >= Atm_xmax) then ! cycle the X location
       loc(1) = mod(loc(1), Atm_xmax)
       if (loc(1) < 0.0_R_) loc(1) = loc(1) + Atm_xmax
    end if
    if (loc(2) < 0.0_R_ .or. loc(2) >= Atm_ymax) then ! cycle the Y location
       loc(2) = mod(loc(2), Atm_ymax)
       if (loc(2) < 0.0_R_) loc(2) = loc(2) + Atm_ymax
    end if
    !call mcarUtl__horiShift(loc(1), dir(1), path, Atm_xmax)
    !call mcarUtl__horiShift(loc(2), dir(2), path, Atm_ymax)

  end subroutine mcarAtm__rt_renewLocXY


  !+
  ! Move to the boundary of the box
  !-
  subroutine mcarAtm__rt_bound3D(icase, path, dir, itr, loc, ix, iy, iz)

    integer,  intent(in) :: icase
    integer,  intent(inout) :: ix
    integer,  intent(inout) :: iy
    integer,  intent(inout) :: iz
    real(R_), intent(in) :: path
    real(R_), intent(inout) :: loc(:)
    real(R_), intent(in) :: dir(:)
    integer, intent(in) :: itr(:)

    ! Z
    if (icase == 3) then
       loc(1:2) = loc(1:2) + path * dir(1:2) ! x,y
       if (itr(3) == 0) then
          loc(3) = Atm_zgrd(iz - 1)
          iz = iz - 1
       else
          loc(3) = Atm_zgrd(iz)
          iz = iz + 1
       end if

       ! Y
    else if (icase == 2) then
       loc(1) = loc(1) + path * dir(1) ! x
       loc(3) = loc(3) + path * dir(3) ! z
       if (itr(2) == 0) then
          loc(2) = Atm_ygrd(iy - 1)
          iy = iy - 1
          if (iy <= 0) then
             loc(2) = Atm_ymax
             iy = Atm_ny
          end if
       else
          loc(2) = Atm_ygrd(iy)
          iy = iy + 1
          if (iy > Atm_ny) then
             loc(2) = 0.0_R_
             iy = 1
          end if
       end if

       ! X
    else
       loc(2:3) = loc(2:3) + path * dir(2:3) ! y,z
       if (itr(1) == 0) then
          loc(1) = Atm_xgrd(ix - 1)
          ix = ix - 1
          if (ix <= 0) then
             loc(1) = Atm_xmax
             ix = Atm_nx
          end if
       else
          loc(1) = Atm_xgrd(ix)
          ix = ix + 1
          if (ix > Atm_nx) then
             loc(1) = 0.0_R_
             ix = 1
          end if
       end if
    end if

  end subroutine mcarAtm__rt_bound3D


  !+
  ! Path length in 1-D layer
  !-
  function mcarAtm__rt_path1D_R(z, iz, uz, itr) result(path) 

    real(R_), intent(in) :: z, uz
    integer,  intent(in) :: iz
    integer,  intent(in) :: itr
    real(R_)  :: path
    path = (Atm_zgrd(iz - 1 + itr) - z) / uz

  end function mcarAtm__rt_path1D_R


  !+
  ! Path to a boundary of the box
  !-
  subroutine mcarAtm__rt_path3D(loc, ix, iy, iz, dir, itr, icase, path)

    real(R_), intent(in) :: loc(:) ! location vector
    real(R_), intent(in) :: dir(:) ! direction vector
    integer,  intent(in) :: itr(:) ! even parity flags
    integer,  intent(in) :: ix
    integer,  intent(in) :: iy
    integer,  intent(in) :: iz
    integer,  intent(out) :: icase
    real(R_), intent(out) :: path
    real(R_) :: dc(3)

    ! Find the closest plane
    dc(3) = Atm_zgrd(iz - 1 + itr(3)) - loc(3)
    dc(2) = Atm_ygrd(iy - 1 + itr(2)) - loc(2)
    dc(1) = Atm_xgrd(ix - 1 + itr(1)) - loc(1)
    icase = 3
    if (abs(dc(2) * dir(icase)) < abs(dc(icase) * dir(2))) icase = 2
    if (abs(dc(1) * dir(icase)) < abs(dc(icase) * dir(1))) icase = 1

    ! Pathlength
    path = dc(icase) / dir(icase) ! should be [0,+Inf), not be +Inf

  end subroutine mcarAtm__rt_path3D


  !+
  ! Path length to the super-grid box boundary
  !-
  subroutine mcarAtm__rt_pathSup(ixs, iys, izs, isup, path, icase)

    integer,  intent(out) :: icase
    integer,  intent(in) :: isup
    integer,  intent(in) :: ixs
    integer,  intent(in) :: iys
    integer,  intent(in) :: izs
    real(R_), intent(out) :: path
    real(R_)  :: dc(3)

    ! Find the closest plane
    dc(3) = Atm_zsup(izs - 1 + Pho_itr(3)) - Pho_loc(3)
    icase = 3
    if (isup == 3 .and. Pho_mv3D) then
       dc(2) = Atm_ysup(iys - 1 + Pho_itr(2), izs) - Pho_loc(2)
       dc(1) = Atm_xsup(ixs - 1 + Pho_itr(1), izs) - Pho_loc(1)
       if (abs(dc(2) * Pho_dir(icase)) < abs(dc(icase) * Pho_dir(2))) icase = 2
       if (abs(dc(1) * Pho_dir(icase)) < abs(dc(icase) * Pho_dir(1))) icase = 1
    end if

    ! Pathlength in the box
    path = dc(icase) / Pho_dir(icase) ! should be [0,+Inf), not be +Inf
    !call check_valid('mcarAtm__rt_pathSup: path', path) !TEST

    ! Unexpected case: Photon leaked from the super-voxel
    if (path < 0.0_R_) then
       if (Pho_loc(3) < Atm_zsup(izs - 1)) then
          Pho_loc(3) = Atm_zsup(izs - 1)
          icase = 3
       else if (Pho_loc(3) > Atm_zsup(izs)) then
          Pho_loc(3) = Atm_zsup(izs)
          icase = 3
       end if
       if (isup == 3 .and. Pho_mv3D) then
          if (Pho_loc(2) < Atm_ysup(iys - 1, izs)) then
             Pho_loc(2) = Atm_ysup(iys - 1, izs)
             icase = 2
          else if (Pho_loc(2) > Atm_ysup(iys, izs)) then
             Pho_loc(2) = Atm_ysup(iys, izs)
             icase = 2
          end if
          if (Pho_loc(1) < Atm_xsup(ixs - 1, izs)) then
             Pho_loc(1) = Atm_xsup(ixs - 1, izs)
             icase = 1
          else if (Pho_loc(1) > Atm_xsup(ixs, izs)) then
             Pho_loc(1) = Atm_xsup(ixs, izs)
             icase = 1
          end if
       end if
       path = 0.0_R_
    end if

  end subroutine mcarAtm__rt_pathSup


  !+
  ! Determine horizontal gridding of super-voxels
  !-
  subroutine mcarAtm__set_supVox(wwrk, awrk, nx, ny, xmax, ymax, fsupg, fsupi, nxx, nyy, nxs, nys)


    real(R_), intent(in) :: fsupg
    real(R_), intent(in) :: fsupi
    integer,  intent(in) :: nx
    integer,  intent(out) :: nxs
    integer,  intent(out) :: nxx
    integer,  intent(in) :: ny
    integer,  intent(out) :: nys
    integer,  intent(out) :: nyy
    real(R_), intent(in) :: xmax
    real(R_), intent(in) :: ymax
    real(R_), intent(inout) :: wwrk(:,:), awrk(:,:) ! max & average values
    real(R_),  parameter :: EXTEPS = RSML_
    real(R_)  :: a, detg, deti
    integer   :: ixend, ixr, ixs, iyend, iyr, iys, nxsnew, nxxnew, nysnew, nyynew
    real(R_)  :: w, xbin, ybin

    ! Initializations
    xbin = xmax / real(nx, R_)
    ybin = ymax / real(ny, R_)
    nxx = 1
    nyy = 1
    nxs = nx
    nys = ny
    ixend = 0
    iyend = 0
    if (nx <= 1) ixend = 1
    if (ny <= 1) iyend = 1

    ! Loop
    do
       if (ixend == 1 .and. iyend == 1) exit

       ! Merge along X-axis
       if (ixend == 0) then
          nxsnew = (nxs + 1) / 2
          nxxnew = nxx * 2
          if (nxsnew == 1) ixend = 1
          detg = 0.0_R_
          deti = 0.0_R_
          do iys = 1, nys
             do ixs = 1, nxsnew
                ixr = 2 * ixs - 1
                if (ixs < nxsnew .or. nxsnew * 2 == nxs) then
                   w = max(wwrk(ixr, iys), wwrk(ixr + 1, iys))
                   a = 0.5_R_ * (awrk(ixr, iys) + awrk(ixr + 1, iys))
                else
                   w = wwrk(ixr, iys)
                   a = awrk(ixr, iys)
                end if
                wwrk(ixs, iys) = w
                awrk(ixs, iys) = a
                detg = detg + w
                if (w > EXTEPS) then
                   deti = deti + a / w
                else
                   deti = deti + 1.0_R_
                end if
             end do
          end do
          detg = detg / real(nxsnew * nys, R_)
          deti = deti / real(nxsnew * nys, R_)
          if (detg * fsupg * xbin * real(nxxnew, R_) > 1.0_R_ .and. deti * fsupi < 1.0_R_) then
             ixend = 1
          else
             nxs = nxsnew
             nxx = nxxnew
          end if
       end if

       ! Merge along Y-axis
       if (iyend == 0) then
          nysnew = (nys + 1) / 2
          nyynew = nyy * 2
          if (nysnew == 1) iyend = 1
          detg = 0.0_R_
          deti = 0.0_R_
          do ixs = 1, nxs
             do iys = 1, nysnew
                iyr = 2 * iys - 1
                if (iys < nysnew .or. nysnew * 2 == nys) then
                   w = max(wwrk(ixs, iyr), wwrk(ixs, iyr + 1))
                   a = 0.5_R_ * (awrk(ixs, iyr) + awrk(ixs, iyr + 1))
                else
                   w = wwrk(ixs, iyr)
                   a = awrk(ixs, iyr)
                end if
                wwrk(ixs, iys) = w
                awrk(ixs, iys) = a
                detg = detg + w
                if (w > EXTEPS) then
                   deti = deti + a / w
                else
                   deti = deti + 1.0_R_
                end if
             end do
          end do
          detg = detg / real(nxs * nysnew, R_)
          deti = deti / real(nxs * nysnew, R_)
          if (detg * fsupg * ybin * real(nyynew, R_) > 1.0_R_ .and. deti * fsupi < 1.0_R_) then
             iyend = 1
          else
             nys = nysnew
             nyy = nyynew
          end if
       end if
    end do

  end subroutine mcarAtm__set_supVox

end module mcarAtm
