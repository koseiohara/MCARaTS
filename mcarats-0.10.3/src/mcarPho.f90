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
! Module for model photon packet
!-
module mcarPho 

  use globals
  use hparx
  implicit none
  private

  ! Public
  public :: mcarPho__tune_init
  public :: mcarPho__user_init
  public :: mcarPho__init      !* Initialize this module *
  public :: mcarPho__final     !* Finalize this module *
  public :: mcarPho__set_tech
  public :: mcarPho__set_iso_SS
  public :: mcarPho__register
  public :: mcarPho__rt_newPhoton
  public :: mcarPho__rt_update
  public :: mcarPho__rt_newDirecV
  public :: mcarPho__rt_ftau_new
  public :: mcarPho__rt_ftau_update
  public :: mcarPho__rt_hist_rec
  public :: mcarPho__rt_newStatus
  public :: mcarPho__rt_scatOrders
  public :: mcarPho__rt_scaleWgt
  public :: mcarPho__rt_killSmall
  public :: mcarPho__rt_duplicate
  public :: mcarPho__rt_ichi_next

  ! Namelist for Pho
  integer,  save :: Pho_iso_SS  = 1 ! max scattering order for 3-D transfer
  integer,  save :: Pho_iso_tru = 0 ! min scattering order for truncation approximations
  integer,  save :: Pho_iso_max = 10000000 ! max scattering order for simulation
  real(R_), save :: Pho_wsml   = 0.05_R_ ! Russian roulette weight limit for local estimates
  real(R_), save :: Pho_wmin   = 0.2_R_  ! min photon weight
  real(R_), save :: Pho_wmax   = 3.0_R_  ! max photon weight
  real(R_), save :: Pho_wfac   = 1.0_R_  ! factor for ideal photon weight
  real(R_), save :: Pho_pfpeak = 30.0_R_ ! phase function peak threshold
  namelist /mcarPho_nml_init/Pho_iso_SS, Pho_iso_tru, Pho_iso_max, Pho_wsml, Pho_wmin, Pho_wmax, &
       & Pho_wfac, Pho_pfpeak

  ! Public (read only) : Common status of photons
  integer,  save :: Pho_ikd    ! index of k-distribution term
  real(R_), save :: Pho_eone   ! energy power
  real(R_), save :: Pho_hdwx, Pho_hdwy
  public :: Pho_ikd, Pho_eone, Pho_hdwx, Pho_hdwy

  ! Public (readable/writable) : Photon packet variables for real trajectory
  logical,  save :: Pho_dead   ! true if the photon is dead (wgt = 0, iso > iso_max)
  real(R_), save :: Pho_wrr    ! RR weight limit for photon trancing
  logical,  save :: Pho_mv3D   ! true if the photon can move in the 3D (iso <= iso_SS)
  real(R_), save :: Pho_ftau   ! free optical thickness
  real(R_), save :: Pho_sdif   ! spread of photon packet (numerical diffusion length)
  real(R_), save :: Pho_loc(3) !# location vector
  integer,  save :: Pho_iz     !# Z index
  integer,  save :: Pho_ichi   ! index of chi, used for truncation approximations
  real(R_), save :: Pho_dir(3) ! direction vector
  integer,  save :: Pho_itr(3) ! Cartesian, transport direction flags
  real(R_), save :: Pho_wgt    ! weight of the power
  real(R_), save :: Pho_chi    ! directionality parameter, chi
  integer,  save :: Pho_iso    ! next scattering order (1: direct beam)
  real(R_), save, allocatable :: Pho_plen(:) ! path-length distribution
  public :: Pho_dead, Pho_ftau, Pho_sdif, Pho_wrr
  public :: Pho_plen, Pho_wgt, Pho_loc, Pho_dir, Pho_itr, Pho_iz
  public :: Pho_iso, Pho_ichi, Pho_chi, Pho_mv3D

  ! Public (readable/writable) : Photon packet variables for virtual trajectory
  logical,  save :: PhoV_mv3D   ! true if the photon can move in the 3D (iso <= iso_SS)
  real(R_), save :: PhoV_loc(3) ! location vector
  integer,  save :: PhoV_iz     ! Z index
  integer,  save :: PhoV_ichi   ! index of chi, used for truncation approximations
  real(R_), save :: PhoV_dir(3) ! direction vector
  integer,  save :: PhoV_itr(3) ! Cartesian, transport direction flags
  real(R_), save :: PhoV_wgt    ! weight of the power
  real(R_), save :: PhoV_chi    ! directionality parameter, chi
  integer,  save :: PhoV_iso    ! next scattering order (1: direct beam)
  real(R_), save, allocatable :: PhoV_plen(:) ! path-length distribution
  public :: PhoV_plen, PhoV_wgt, PhoV_loc, PhoV_dir, PhoV_itr, PhoV_iz
  public :: PhoV_iso, PhoV_ichi, PhoV_chi, PhoV_mv3D

  ! Private
  integer,  save :: Pho_nz
  integer,  save :: Pho_nchi
  real(R_), save :: Pho_ppeak = 30.0_R_ / (4.0_R_ * PI_) ! threshold for adaptive truncation
  !// Truncation order depends on P, probability per solid angle (steradian)
  !   For P >= ppeak, ichinew = ichi (the same truncation order)
  !   For P <  ppeak, ichinew = ichi + 1 (the truncation order will be updated)

contains

  !+
  ! Tune a preset for technical parameters
  !-
  subroutine mcarPho__tune_init(mtarget) 

    integer, intent(in) :: mtarget ! target flag [0,2]
    !// = 0: no optimazation (deactivate all default optimization)
    !     1: optimization for flux & HR
    !     2: optimization for radiance

    if (mtarget <= 0) then ! no optimization
       call mcarPho__set_tech(2.0_R_, 0.01_R_, 1.0_R_, 0.0_R_)
    else if (mtarget == 1) then ! for fluxes & HRs
       call mcarPho__set_tech(3.0_R_, 0.2_R_, 1.01_R_, 300.0_R_)
    else ! for radiances
       call mcarPho__set_tech(3.0_R_, 0.2_R_, 1.01_R_, 30.0_R_)
    end if

  end subroutine mcarPho__tune_init
 
  !+
  ! Set technical parameters
  !-
  subroutine mcarPho__set_tech(wmax, wmin, wfac, pfpeak) 

    real(R_), intent(in) :: wmax, wmin, wfac, pfpeak
    Pho_wmax = wmax
    Pho_wmin = wmin
    Pho_wfac = wfac
    Pho_pfpeak = pfpeak

  end subroutine mcarPho__set_tech

  !+
  ! Read in namelist variables
  !-
  subroutine mcarPho__user_init(iu) 

    integer, intent(in) :: iu   ! file unit index for namelist input
    integer :: ios
    read (iu, nml=mcarPho_nml_init, iostat=ios)
    call err_read(ios, iu, 'mcarPho__user_init: Invalid namelist input for mcarPho_nml_init.')
    call check_iI('mcarPho__user_init: Pho_iso_max', Pho_iso_max, 1)
    call check_iR('mcarPho__user_init: Pho_wmin',    Pho_wmin,    0.0_R_,  1.0_R_)
    call check_iR('mcarPho__user_init: Pho_wmax',    Pho_wmax,    1.0_R_,  100.0_R_)
    call check_iR('mcarPho__user_init: Pho_wfac',    Pho_wfac,    Pho_wmin, Pho_wmax)

  end subroutine mcarPho__user_init

  !+
  ! Initialize this module
  !-
  subroutine mcarPho__init(nz, nchi) 

    integer, intent(in) :: nz, nchi
    Pho_nz = nz
    Pho_nchi = nchi
    Pho_ppeak = Pho_pfpeak / (4.0_R_ * PI_)
    if (allocated(Pho_plen)) call mcarPho__final()
    allocate (Pho_plen(0:nz+1))
    allocate (PhoV_plen(0:nz+1))

  end subroutine mcarPho__init

  !+
  ! Finalize this module
  !-
  subroutine mcarPho__final() 

    if (allocated(Pho_plen)) then
       deallocate (Pho_plen, PhoV_plen)
    end if

  end subroutine mcarPho__final

  !+
  ! Set the 3D-to-1D scheme-shift order
  !-
  subroutine mcarPho__set_iso_SS(isol) 

    integer, intent(in) :: isol ! solver flag
    if (isol <= 0) then      ! F3D
       Pho_iso_SS = 1000000
    else if (isol >= 2) then ! ICA
       Pho_iso_SS = 0
    end if

  end subroutine mcarPho__set_iso_SS

  !+
  ! Register common photon status
  !-
  subroutine mcarPho__register(ikd, eone, hdwx, hdwy) 

    integer,  intent(in) :: ikd
    real(R_), intent(in) :: eone, hdwx, hdwy
    Pho_ikd = ikd
    Pho_eone = eone
    Pho_hdwx = hdwx
    Pho_hdwy = hdwy

  end subroutine mcarPho__register

  !+
  ! Generate a new photon packet
  !-
  subroutine mcarPho__rt_newPhoton(wgt0, dir0, loc0, iz0, chi0)

    real(R_), intent(in) :: wgt0, dir0(:), loc0(:), chi0
    integer,  intent(in) :: iz0
    real(R_), parameter :: EPSU = REPS_

    Pho_dead = .false.
    Pho_iso = 0
    Pho_chi = chi0
    Pho_ichi = 1
    Pho_wgt = wgt0
    Pho_wrr = max(Pho_wmin, min(Pho_wmax, Pho_wgt))
    Pho_loc(:) = loc0(:)
    Pho_iz = iz0
    Pho_sdif = 0.0_R_
    Pho_plen(:) = 0.0_R_
    Pho_ftau = 0.0_R_
    Pho_dir(:) = dir0(:)
    Pho_dir(3) = nonZero_R(Pho_dir(3), EPSU)
    Pho_itr(:) = mcarPho__rt_itr_1R(Pho_dir)

  end subroutine mcarPho__rt_newPhoton


  !+
  ! Even parity flags
  !-
  function mcarPho__rt_itr_1R(dir) result(itr) 

    real(R_), intent(in) :: dir(:) ! a new direction
    integer :: itr(size(dir))

    if (dir(1) < 0.0_R_) then
       itr(1) = 0
    else
       itr(1) = 1
    end if
    if (dir(2) < 0.0_R_) then
       itr(2) = 0
    else
       itr(2) = 1
    end if
    if (dir(3) < 0.0_R_) then
       itr(3) = 0
    else
       itr(3) = 1
    end if

  end function mcarPho__rt_itr_1R


  !+
  ! Update direction, free optical thickness & truncation order (real trajectory)
  !-
  subroutine mcarPho__rt_update(dir, pps)

    real(R_), intent(in) :: dir(:)
    real(R_), intent(in) :: pps  ! probability per solid angle for the new direction    
    real(R_), parameter :: EPSU = REPS_

    Pho_dir(:) = dir(:)
    Pho_dir(3) = nonZero_R(Pho_dir(3), EPSU)
    Pho_itr(:) = mcarPho__rt_itr_1R(Pho_dir)

    Pho_ftau = rand_exp_R() ! new random optical thickness

    if (Pho_iso <= Pho_iso_tru) then ! discard truncation approx.
       Pho_ichi = 1
    else if (pps < Pho_ppeak) then ! increse only for non-peak scattering
       Pho_ichi = min(Pho_nchi, Pho_ichi + 1)
    end if

  end subroutine mcarPho__rt_update


  !+
  ! Register a new direction (virtual trajectory)
  !-
  subroutine mcarPho__rt_newDirecV(dir) 

    real(R_), intent(in) :: dir(:)
    real(R_), parameter :: EPSU = REPS_

    PhoV_dir(:) = dir(:)
    PhoV_dir(3) = nonZero_R(PhoV_dir(3), EPSU)
    PhoV_itr(:) = mcarPho__rt_itr_1R(PhoV_dir)

  end subroutine mcarPho__rt_newDirecV


  !+
  ! Update truncation order for virtual trajectory
  !-
  subroutine mcarPho__rt_ichi_next(pps) 

    real(R_), intent(in) :: pps  ! probability per solid angle for the new direction    
    if (PhoV_iso <= Pho_iso_tru) then ! discard truncation approx.
       PhoV_ichi = 1
    else if (pps < Pho_ppeak) then ! increse only for non-peak scattering
       PhoV_ichi = min(Pho_nchi, Pho_ichi + 1)
    else
       PhoV_ichi = Pho_ichi
    end if

  end subroutine mcarPho__rt_ichi_next


  !+
  ! New free optical thickness
  !-
  subroutine mcarPho__rt_ftau_new(mtau) 

    integer,  intent(in) :: mtau ! method [0,1]
    if (mtau == 0) then
       Pho_ftau = rand_exp_R()
    else
       Pho_ftau = Pho_ftau + rand_exp_R()
    end if
    
  end subroutine mcarPho__rt_ftau_new


  !+
  ! Update free optical thickness
  !-
  subroutine mcarPho__rt_ftau_update(dtau) 

    real(R_), intent(in) :: dtau ! optical thickness the photon travelled
    Pho_ftau = max(0.0_R_, Pho_ftau - dtau) ! to avoid a fatal error with ftau < 0
    
  end subroutine mcarPho__rt_ftau_update


  !+
  ! Record history of the trajectry
  !-
  subroutine mcarPho__rt_hist_rec(iz, path) 

    integer,  intent(in) :: iz   ! layer index the history is recorded   
    real(R_), intent(in) :: path ! pathlength segment
    Pho_sdif = Pho_sdif + (1.0_R_ - Pho_chi) * path
    Pho_plen(iz) = Pho_plen(iz) + path

  end subroutine mcarPho__rt_hist_rec


  !+
  ! Randomly kill the photon packet with small weight
  !-
  subroutine mcarPho__rt_killSmall() 

    real(R_), parameter :: WACT = 0.7_R_
    if (Pho_wgt < Pho_wrr * WACT) then
       if (Pho_wgt > Pho_wrr * mseq_rand_R()) then
          Pho_wgt = Pho_wrr ! survive
       else
          Pho_wgt = 0.0_R_ ! killed
          Pho_dead = .true.
       end if
    end if

  end subroutine mcarPho__rt_killSmall


  !+
  ! Scale the weight and randomly kill the photon packet with too smalle weight
  !-
  subroutine mcarPho__rt_scaleWgt(fac)

    real(R_), intent(in) :: fac
    real(R_), parameter  :: WACT = 0.7_R_

    Pho_wgt = Pho_wgt * fac
    if (Pho_wgt < Pho_wsml * WACT) then
       if (Pho_wgt > Pho_wsml * mseq_rand_R()) then
          Pho_wgt = Pho_wsml ! survive
       else
          Pho_wgt = 0.0_R_ ! killed
          Pho_dead = .true.
       end if
    end if

  end subroutine mcarPho__rt_scaleWgt


  !+
  ! New scattering order & anisotropy parameters
  !-
  subroutine mcarPho__rt_newStatus(gtc)

    real(R_), intent(in) :: gtc

    PhoV_chi = Pho_chi * gtc ! anisotropy parameter
    PhoV_iso = Pho_iso + 1 ! scattering order
    if (PhoV_iso <= Pho_iso_SS) then ! tentative
       PhoV_mv3D = .true.
    else
       PhoV_mv3D = .false.
    end if

  end subroutine mcarPho__rt_newStatus


  !+
  ! New scattering order & anisotropy parameters
  !-
  subroutine mcarPho__rt_scatOrders()

    ! Update
    Pho_chi  = PhoV_chi  ! anisotropy parameter
    Pho_iso  = PhoV_iso  ! scattering order

    ! Status by scattering order
    if (Pho_iso > Pho_iso_max) then ! tentative
       Pho_dead = .true.
       return
    else
       Pho_dead = .false.
    end if
    if (Pho_iso <= Pho_iso_SS) then ! tentative
       Pho_mv3D = .true.
    else
       Pho_mv3D = .false.
    end if

    ! Change Wrr, weight for Russian roulette
    Pho_wrr = max(Pho_wmin, min(Pho_wmax, Pho_wrr * Pho_wfac))

  end subroutine mcarPho__rt_scatOrders


  !+
  ! Duplicate a photon packet (for local estimates)
  !-
  subroutine mcarPho__rt_duplicate() 

    PhoV_wgt     = Pho_wgt
    PhoV_iz      = Pho_iz
    PhoV_loc(:)  = Pho_loc(:)
    PhoV_plen(:) = Pho_plen(:)

  end subroutine mcarPho__rt_duplicate

end module mcarPho
