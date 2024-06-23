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
! Module for Rahman-Pinty-Varstlaete (RPV) BRDF model
! Reference:
!   Rahman, H., B. Pinty, & M. M. Verstraete, 1993: Coupled surface-atmosphere
!   reflectance (CSAR) model. 2.0_R_ Semiempirical surface model usable with NOAA
!   Advanced High Resolution Radiometer data. J. Geophys. Res., 98, D11,
!   20791-20801
!-
module mcarRpv

  use globals
  use hparx
  use mcarUtl
  implicit none
  private

  ! Public
  public :: mcarRpv__BRDFpA_R
  public :: mcarRpv__fitUz1
  public :: mcarRpv__funcB_R
  public :: mcarRpv__funcC_R
  public :: mcarRpv__prep_makeTab
  public :: mcarRpv__newDir

contains

  !+
  ! BRDF/albedo of RPV model
  !-
  function mcarRpv__BRDFpA_R(pac, fa, vec0, vec1) result(res)

    real(R_), intent(in) :: fa      ! function A
    real(R_), intent(in) :: vec0(:) ! incoming direction vector (downward motion to the plane)
    real(R_), intent(in) :: vec1(:) ! outgoing direction vector (upward motion from the plane)
    real(R_), intent(in) :: pac(:)  ! pac(1-4) : parameter packet
    !    pac(1) : "rho0", the first parameter
    !    pac(2) : "k", anisotroy integer,   parameter :: if isotropic then k=1
    !    pac(3) : "THETA", parameter for the H-G type function
    !    pac(4) : "delta", sharpness parameter for hot spot
    real(R_)  :: res
    real(R_)  :: cosa, deno, f, h, rnum, s, tt, uzuz

    uzuz = -vec0(3) * vec1(3)
    cosa = -geo3D_scaProd_R(vec0, vec1)
    tt = pac(3) * pac(3)
    s = 1.0_R_ + tt + 2.0_R_ * pac(3) * cosa
    f = (1.0_R_ - tt) / (s * sqrt(s))
    h = (1.0_R_ - pac(1)) * uzuz / (pac(4) * uzuz &
         & + sqrt(max(0.0_R_, vec0(3)**2 + vec1(3)**2 - 2.0_R_ * uzuz * cosa)))

    rnum = pac(1) * f * (1.0_R_ + h)
    deno = PI_ * fa * (vec1(3) * (vec1(3) - vec0(3)))**(1.0_R_ - pac(2))

    if (deno > REPS_ * rnum) then
       res = rnum / deno
    else
       res = 1.0_R_ / REPS_
    end if

  end function mcarRpv__BRDFpA_R


  !+
  ! Get the function A of the RPV BRDF model
  !  and simultaneously fit a linear function to the function B
  !-
  subroutine mcarRpv__fitUz1(pac, uz0, nuz1, nphi, funca, a1, a0, v1, v0)

    real(R_), intent(out) :: a0, a1 ! regression coeffficients for y = a1 * x + a0
    real(R_), intent(out) :: funca ! the function A
    integer,  intent(in) :: nphi ! number of quadrature points for phi
    integer,  intent(in) :: nuz1 ! number of quadrature points for uz1
    real(R_), intent(in) :: uz0 ! cosine of incoming zenith angle (downward motion to the plane)
    !          Note uz0 should be < or = 0
    real(R_), intent(out) :: v0, v1 ! variance of the coefficients of a1 & a0, respectively
    real(R_), intent(in) :: pac(:) ! pac(1-4) : parameter packet
    !     pac(1) : "rho0", the first parameter
    !     pac(2) : "k", anisotroy integer,   parameter :: if isotropic then k=1
    !     pac(3) : "THETA", parameter for the H-G type function
    !     pac(4) : "delta", sharpness parameter for hot spot
    real(R_),  parameter :: FPI = 2.0_R_ / PI_
    integer,   parameter :: KNQUAD = 200
    real(R_), save :: x(KNQUAD), w(KNQUAD), stt, sx
    integer,  save :: nx = 0
    integer   :: ix, n
    real(R_)  :: dx, sa, sty, sy, wy

    ! Initializations
    n = min(KNQUAD, nuz1)
    if (nx /= n) then
       nx = n
       call gaussLegen(nx, x, w)
       sx = 0.0_R_
       do ix = 1, nx
          x(ix) = 0.5_R_ * (1.0_R_ + x(ix))
          sx = sx + w(ix) * x(ix)
       end do
       stt = 0.0_R_
       do ix = 1, nx
          dx = x(ix) - sx
          stt = stt + w(ix) * (dx * dx)
       end do
    end if

    ! Integration of B
    sa  = 0.0_R_
    sty = 0.0_R_
    sy  = 0.0_R_
    do ix = 1, nx
       wy = w(ix) * mcarRpv__funcB_R(pac, uz0, x(ix), nphi)
       sa  = sa  + wy
       sty = sty + wy * (x(ix) - sx)
       sy  = sy  + wy
    end do

    ! Function A, regression coeffs. & uncertainties
    funca = sa * FPI * pac(1)
    a1 = sty / stt
    a0 = sy - sx * a1
    v1 = 1.0_R_ / stt
    v0 = 1.0_R_ + sx * sx * v1

  end subroutine mcarRpv__fitUz1


  !+
  ! Function B of the RPV BRDF model
  !-
  function mcarRpv__funcB_R(pac, uz0, uz1, nphi) result(res)

    integer,  intent(in) :: nphi ! number of quadrature points for phi
    real(R_), intent(in) :: uz0  ! cosine of incoming zenith angle (downward motion to the plane)
    real(R_), intent(in) :: uz1  ! cosine of outgoing zenith angle (upward motion from the plane)
    real(R_), intent(in) :: pac(:) ! pac(1-4) : parameter packet
    !     pac(1) : "rho0", the first parameter
    !     pac(2) : "k", anisotroy integer,   parameter :: if isotropic then k=1
    !     pac(3) : "THETA", parameter for the H-G type function
    !     pac(4) : "delta", sharpness parameter for hot spot
    real(R_)  :: res
    integer,   parameter :: KNQUAD = 200
    real(R_), save :: x(KNQUAD), w(KNQUAD)
    integer,  save :: nx = 0
    integer   :: ix, n
    real(R_)  :: tsum

    ! Initializations
    n = min(KNQUAD, nphi)
    if (nx /= n) then
       nx = n
       call gaussLegen(nx, x, w)
       do ix = 1, nx
          x(ix) = cos(0.5_R_ * PI_ * (1.0_R_ + x(ix)))
       end do
    end if

    ! Integration of C
    tsum = w(1) * mcarRpv__funcC_R(1, pac, uz0, uz1, x(1))
    do ix = 2, nx
       tsum = tsum + w(ix) * mcarRpv__funcC_R(0, pac, uz0, uz1, x(ix))
    end do
    tsum = PI_ * tsum
    res = tsum * (uz1 * (uz1 - uz0))**pac(2) / (uz1 - uz0) ! function B

  end function mcarRpv__funcB_R


  !+
  ! Function C in the RPV BRDF model
  !-
  function mcarRpv__funcC_R(init, pac, uz0, uz1, cosphi) result(res)

    real(R_), intent(in) :: cosphi ! cosine of relative azimuth angle
    integer,  intent(in) :: init ! If 1 then initialization, 
    !          else then the initialization is skipped.
    !          Give 0 if pac, uz0, & uz1 are the same as the previous call.
    real(R_), intent(in) :: uz0 ! cosine of incoming zenith angle (downward motion to the plane)
    real(R_), intent(in) :: uz1 ! cosine of outgoing zenith angle (upward motion from the plane)
    real(R_), intent(in) :: pac(:) ! pac(1-4) : parameter packet
    !    pac(1) : "rho0", the first parameter
    !    pac(2) : "k", anisotroy integer,   parameter :: if isotropic then k=1
    !    pac(3) : "THETA", parameter for the H-G type function
    !    pac(4) : "delta", sharpness parameter for hot spot
    real(R_)  :: res
    real(R_)  :: dcosg, f, h, s, tt
    real(R_), save :: d1 = 0.0_R_, d2 = 0.0_R_, d3 = 0.0_R_, d4 = 0.0_R_, d5 = 0.0_R_

    ! Initialization
    if (init == 1) then
       tt = pac(3) * pac(3)
       d1 = 1.0_R_ - tt
       d2 = 1.0_R_ + tt
       d3 = -uz0 * uz1
       d4 = sqrt(max(0.0_R_, (1.0_R_ - uz0*uz0) * (1.0_R_ - uz1*uz1)))
       d5 = uz0 * uz0 + uz1 * uz1
    end if

    ! The function C
    dcosg = 2.0_R_ * (d3 + d4 * cosphi)
    s = d2 + pac(3) * dcosg
    f = d1 / (s * sqrt(s))
    h = (1.0_R_ - pac(1)) * d3 / (pac(4) * d3 + sqrt(max(0.0_R_, d5 - d3 * dcosg)))
    res = f * (1.0_R_ + h)

  end function mcarRpv__funcC_R


  !+
  ! Get albedo and coefficients a and b of the RPV BRDF model
  !-
  subroutine mcarRpv__prep_makeTab(pac, bmin, bsex, uz0min, nuz0, nuz1, nphi, &
       &     falut, calut, cblut, albmin, iam)

    real(R_), intent(in) :: bmin ! minimum threshold for factor b' (about 3 is recommended)
    real(R_), intent(in) :: bsex ! scaling exponent for b' (should be < or = 1)
    !            Note: If s=1, higher accuracy of angular distribution of
    !            reflection will be achieved; about 0.5_R_ is recommended.
    integer,  intent(in) :: nphi ! number of quadrature points for phi
    integer,  intent(in) :: nuz1 ! number of quadrature points for theta
    !            This is also used as number of grid points for finding "b"
    real(R_), intent(in) :: uz0min ! min cosines of incoming zenith angle
    real(R_), intent(in) :: pac(:) ! pac(1-4) : parameter packet
    !     pac(1) : "rho0", the first parameter
    !     pac(2) : "k", anisotroy integer,   parameter :: if isotropic then k=1
    !     pac(3) : "THETA", parameter for the H-G type function
    !     pac(4) : "delta", sharpness parameter for hot spot
    integer,  intent(inout) :: nuz0 ! number of cosines of incoming zenith angle (should be > 3)
    real(R_), intent(out) :: falut(:) ! function A for iuz0=1-nuz0
    real(R_), intent(out) :: calut(:) ! coeffficient a1
    real(R_), intent(out) :: cblut(:) ! coeffficient b'
    real(R_), intent(out) :: albmin ! minimum albedo in the table
    integer,  intent(out) :: iam ! iuz0 corresponding to the minimum albedo
    real(R_),  parameter :: A1MIN = -1.0_R_ / PI_ * 0.5_R_ 
    real(R_),  parameter :: A1MAX =  1.0_R_ / PI_ * 0.99_R_
    !// thresholds for slope of the linear comparison function
    real(R_)  :: a0, a1, aa0, aa1, alb, b, bb, buz0, buz1, fac, funca
    integer   :: igrid, iuz0
    integer   :: ngrid
    real(R_)  :: pwr, vec0(3), vec1(3), v0, v1

    ! Macros
    nuz0 = max(4, nuz0)
    buz0 = (1.0_R_ - uz0min) / (nuz0 - 1)
    ngrid = max(nuz1, 10)     ! # of points to find "b"
    pwr = pac(2) - 1.0_R_

    ! Loop for various uz0
    albmin = 1.0_R_
    iam = 1
    do iuz0 = 1, nuz0
       vec0(3) = -min(1.0_R_, uz0min + buz0 * (iuz0 - 1))
       vec0(1) = sqrt(1.0_R_ - vec0(3)**2)
       vec0(2) = 0.0_R_

       ! Function A & Linear regression
       call mcarRpv__fitUz1(pac, vec0(3), nuz1, nphi, funca, aa1, aa0, v1, v0)
       a1 = aa1 / (PI_ * (aa1 + 2.0_R_ * aa0))
       a1 = max(A1MIN, min(A1MAX, a1))
       a0 = 0.5_R_ * (1.0_R_ / PI_ - a1)

       ! Find the factor b for the comparison function
       vec1(2) = 0.0_R_
       vec1(3) = -vec0(3)    ! hot spot
       vec1(1) = -sqrt(1.0_R_ - vec1(3)**2)
       b = vec1(3) / (a1 * vec1(3) + a0) * mcarRpv__BRDFpA_R(pac, funca, vec0, vec1)
       vec1(3) = 1.0_R_     ! zenith
       vec1(1) = REPS_
       b = max(b, mcarRpv__BRDFpA_R(pac, funca, vec0, vec1) / (a1 + a0))
       buz1 = 1.0_R_ / ngrid
       do igrid = 0, ngrid - 1      ! other discrete points
          vec1(3) = 1.0_R_ - buz1 * igrid
          vec1(1) = -sqrt(1.0_R_ - vec1(3)**2)
          fac = vec1(3) / (a1 * vec1(3) + a0)
          b = max(b, fac * mcarRpv__BRDFpA_R(pac, funca, vec0, vec1))
          vec1(1) = -vec1(1)
          b = max(b, fac * mcarRpv__BRDFpA_R(pac, funca, vec0, vec1))
       end do
       b = b + 0.01_R_ ! add torrelence
       if (b < bmin) then
          bb = b
       else
          bb = bmin - 1.0_R_ + (b - bmin + 1.0_R_)**bsex
       end if

       ! Results
       falut(iuz0) = funca
       calut(iuz0) = a1
       cblut(iuz0) = bb

       ! Min albedo
       alb = funca * (-vec0(3))**pwr
       if (alb < albmin) then
          albmin = alb
          iam = iuz0
       end if
    end do

  end subroutine mcarRpv__prep_makeTab


  !+
  ! Scattering of RPV model : determination of a new direction and angular PDF
  ! Algorithm
  !   - The rejection method is used with a comparison function of linear function.
  !-
  subroutine mcarRpv__newDir(vec0, pac, fa, a1, bb, vec1, pps)

    real(R_), intent(in) :: vec0(:) ! incoming direction vector (downward motion to the plane)
    real(R_), intent(in) :: a1      ! coefficients a1
    real(R_), intent(in) :: bb      ! coefficients b'
    real(R_), intent(in) :: fa      ! function A
    real(R_), intent(in) :: pac(:)  ! parameter packet
    !     pac(1) : "rho0", the first parameter
    !     pac(2) : "k", anisotroy integer,   parameter :: if isotropic then k=1
    !     pac(3) : "THETA", parameter for the H-G type function
    !     pac(4) : "delta", sharpness parameter for hot spot
    real(R_), intent(out) :: vec1(:) ! outgoing direction vector (upward motion from the plane)
    real(R_), intent(out) :: pps     ! angular PDF (/steradian)
    real(R_),  parameter :: API = 1.0_R_ / PI_
    real(R_),  parameter :: R2MIN = REPS_**2
    real(R_)  :: a0, cosf, d, deno, r2, sinf, vz
    integer   :: itrial

    ! Rejection method
    a0 = 0.5_R_ * (API - a1)
    do itrial = 1, 200

       ! Tentative direction
       call rand_point_circle(R2MIN, r2, sinf, cosf)
       if (abs(a1) > REPS_) then
          d = a0 * a0 + a1 * r2 * API
          vec1(3) = min(1.0_R_, (sqrt(d) - a0) / a1)
       else
          vec1(3) = r2
       end if
       vz = sqrt(1.0_R_ - vec1(3)**2)
       vec1(1) = vz * cosf
       vec1(2) = vz * sinf

       !  Rejection?
       pps = vec1(3) * mcarRpv__BRDFpA_R(pac, fa, vec0, vec1)
       deno = bb * (a1 * vec1(3) + a0)
       if (pps > deno) exit ! unexpected case
       if (mseq_rand_R() * deno <= pps) exit
    end do
    call geo3D_UVec_renorm(vec1) ! correction of the unit vector

  end subroutine mcarRpv__newDir

end module mcarRpv
