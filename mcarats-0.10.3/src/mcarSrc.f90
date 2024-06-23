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
! Module for radiation sources
!-
module mcarSrc 

  use globals
  use hparx
  use mcarAtm
  use mcarSfc
  implicit none
  private

  ! Public
  public :: mcarSrc__user_init
  public :: mcarSrc__init
  public :: mcarSrc__final
  public :: mcarSrc__user_job
  public :: mcarSrc__prep_isrc

  ! (Public, read-only) User variables in namelist (1) for initialization
  integer,  save :: Src_nsrc = 1 ! # of sources
  integer,  save :: Src_mtherm = 1 ! flag for thermal source (1=yes, 0=no)
  namelist /mcarSrc_nml_init/Src_nsrc, Src_mtherm
  public :: Src_nsrc

  ! (Public, read-only) User variables in namelist (2) for core calculations
  integer :: i
  integer, parameter :: KNSRC = 3000
  integer,  save :: Src_mtype(KNSRC) = (/(1, i = 1, KNSRC)/) ! flag for source type (0,1,2,3)
  integer,  save :: Src_mphi(KNSRC)  = (/(0, i = 1, KNSRC)/) ! flag for random azimuth
  real(R_), save :: Src_wlen(KNSRC)  = (/(10.0_R_,  i = 1, KNSRC)/) ! wavelength for thermal sorce
  real(R_), save :: Src_dwlen(KNSRC) = (/(0.1_R_,   i = 1, KNSRC)/) ! wevelength width
  real(R_), save :: Src_flx(KNSRC)   = (/(1.0_R_,   i = 1, KNSRC)/) ! source flux density or power
  real(R_), save :: Src_qmax(KNSRC)  = (/(0.0_R_,   i = 1, KNSRC)/) ! full cone angle
  real(R_), save :: Src_apsize(KNSRC) = (/(0.0_R_,  i = 1, KNSRC)/) ! aperture size
  real(R_), save :: Src_the(KNSRC)   = (/(120.0_R_, i = 1, KNSRC)/) ! zenith angle
  real(R_), save :: Src_phi(KNSRC)   = (/(0.0_R_,   i = 1, KNSRC)/) ! azimuth angle
  real(R_), save :: Src_xpos(KNSRC)  = (/(0.5_R_,   i = 1, KNSRC)/) ! X relative position
  real(R_), save :: Src_ypos(KNSRC)  = (/(0.5_R_,   i = 1, KNSRC)/) ! Y relative position
  real(R_), save :: Src_zloc(KNSRC)  = (/(0.0_R_,   i = 1, KNSRC)/) ! Z location
  namelist /mcarSrc_nml_job/Src_wlen, Src_dwlen, Src_flx, Src_the, Src_phi, Src_qmax, Src_xpos, &
       & Src_ypos, Src_zloc, Src_apsize, Src_mtype, Src_mphi
  public :: Src_wlen, Src_dwlen, Src_flx, Src_the, Src_phi, Src_qmax, Src_xpos, &
       & Src_ypos, Src_zloc, Src_apsize, Src_mtype, Src_mphi

  ! (Public, read-only) Module variables
  real(R_), save :: Src_esfc  ! total surface source power (W/micron)
  real(R_), save, allocatable :: Src_fsfc(:,:)   ! surface source distribution
  real(R_), save, allocatable :: Src_batm(:,:,:) ! atmospheric Planck functions
  public :: Src_esfc, Src_fsfc, Src_batm

contains

  !+
  ! Read in namelist variables (1) for initialization
  !-
  subroutine mcarSrc__user_init(iu) 

    integer, intent(in) :: iu   ! file unit index for namelist input
    integer :: ios
    read (iu, nml=mcarSrc_nml_init, iostat=ios)
    call err_read(ios, iu, 'mcarSrc__user_init: Invalid namelist input for mcarSrc_nml_init.')

  end subroutine mcarSrc__user_init


  !+
  ! Initialize this module
  !-
  subroutine mcarSrc__init() 

    if (allocated(Src_fsfc)) call mcarSrc__final()
    if (Src_mtherm == 1) then
       allocate (Src_fsfc(Sfc_nxb, Sfc_nyb))
       allocate (Src_batm(Atm_nx, Atm_ny, Atm_nz + 1))
    end if

  end subroutine mcarSrc__init


  !+
  ! Finalize this module
  !-
  subroutine mcarSrc__final() 

    if (allocated(Src_fsfc)) then
       deallocate (Src_fsfc, Src_batm)    
    end if

  end subroutine mcarSrc__final


  !+
  ! Read in namelist variables (2) for core calculations
  !-
  subroutine mcarSrc__user_job(iu) 

    integer, intent(in) :: iu   ! file unit index for namelist input
    integer :: ios, isrc

    ! Read in
    read (iu, nml=mcarSrc_nml_job, iostat=ios)
    call err_read(ios, iu, 'mcarSrc__user_job: Invalid namelist input for mcarSrc_nml_job.')

    ! Fix data
    do isrc = 1, Src_nsrc
       if (Src_the(isrc) > 90.0_R_) then
          Src_the(isrc) = max(90.01_R_, Src_the(isrc))
       else
          Src_the(isrc) = min(89.99_R_, Src_the(isrc))
       end if
       if (Src_mtype(isrc) >= 2 .and. Src_mtherm /= 1) then
          Src_mtherm = 1 
          call mcarSrc__init()
       end if
    end do

    ! Check values
    do isrc = 1, Src_nsrc
       call check_iR('mcarSrc__user_job: Src_the',   Src_the(isrc),   0.0_R_, 180.0_R_)
       call check_iR('mcarSrc__user_job: Src_qmax',  Src_qmax(isrc),  0.0_R_, 180.0_R_)
    end do

  end subroutine mcarSrc__user_job


  !+
  ! Prepare variables for a single source
  !-
  subroutine mcarSrc__prep_isrc(isrc, tmpmin, tmpmax, sdir, q0, f0, dq0)

    integer,  intent(in) :: isrc
    real(R_), intent(in) :: tmpmin, tmpmax
    real(R_), intent(out) :: sdir(:)
    real(R_), intent(out) :: q0, f0, dq0
    integer,  parameter :: NPLKTAB = 100000, NQMAX = 100
    integer  :: ix, ixb, iy, iyb, iz
    real(R_) :: wlen1(1), dwlen1(1), emt, sumf

    ! Local or solar source
    q0  = Src_the(isrc)  * DTOR_
    f0  = Src_phi(isrc)  * DTOR_
    dq0 = Src_qmax(isrc) * DTOR_ * 0.5_R_ ! a half cone angle
    sdir(:) = geo3D_aVec_1R(1.0_R_, q0, f0)

    ! Thermal source
    if (Src_mtype(isrc) >= 2) then ! Planck function table
       wlen1(1)  = Src_wlen(isrc)
       dwlen1(1) = Src_dwlen(isrc)
       call plkTab_init(tmpmin, tmpmax, NPLKTAB, wlen1, dwlen1, NQMAX) ! Planck function LUT

       ! Atmosphere Planck functions
       do iy = 1, Atm_ny
          do ix = 1, Atm_nx
             do iz = 1, Atm_nz + 1
                Src_batm(ix, iy, iz) = plkTab_intp_R(Atm_tmpa3d(ix, iy, iz), 1)
             end do
          end do
       end do

       ! Surface source powers
       do iyb = 1, Sfc_nyb
          do ixb = 1, Sfc_nxb
             emt = 1.0_R_ - Sfc_salbs(1, ixb, iyb) ! surface emittance
             Src_fsfc(ixb, iyb) = PI_ * emt * plkTab_intp_R(Sfc_tmps2d(ixb, iyb), 1) ! fluxes
          end do
       end do
       sumf = sum(Src_fsfc)
       Src_esfc = sumf * Atm_xmax * Atm_ymax / (Sfc_nxb * Sfc_nyb) ! surface source power
       Src_fsfc(:,:) = Src_fsfc(:,:) / sumf ! normalized distribution
    end if

  end subroutine mcarSrc__prep_isrc

end module mcarSrc
