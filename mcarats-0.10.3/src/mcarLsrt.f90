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
! Module for Li-Sparse-Ross-Thick (LSRT) BRDF model
! Reference:
!  Lucht, W., C.B. Schaaf, A.H. Strahler: An algorithm for the retrieval
!  of albedo from space using semiempirical BRDF models. IEEE Trans. Geosci.
!  Remote Sens., vol. 38, 977-998 (2000).
!-
module mcarLsrt

  use globals
  use hparx
  use mcarUtl
  implicit none
  private

  ! Public
  public :: mcarLsrt__BRDF_R
  public :: mcarLsrt__albKernel
  public :: mcarLsrt__BRFKernel
  public :: mcarLsrt__funcFg_R
  public :: mcarLsrt__funcFv_R
  public :: mcarLsrt__prep_makeTab
  public :: mcarLsrt__newDir

contains

  !+
  ! BRDF of LSRT model
  !-
  function mcarLsrt__BRDF_R(pac, vec0, vec1) result(res)

    real(R_), intent(in) :: pac(:)  !  pac(1-3) : parameter packet
    !    pac(1) : "k_L", weight of Lambertian reflection
    !    pac(2) : "k_g", weight of geometrical optics reflection
    !    pac(3) : "k_v", weight of volume scattering
    real(R_), intent(in) :: vec0(:) ! incoming direction vector (downward motion to the plane)
    real(R_), intent(in) :: vec1(:) ! outgoing direction vector (upward motion from the plane)
    real(R_)  :: res
    real(R_)  :: fg, fv
    real(R_),  parameter :: API = 1.0_R_ / PI_

    ! The kernel functions f_g & f_v
    fg = mcarLsrt__funcFg_R(vec0, vec1)
    fv = mcarLsrt__funcFv_R(vec0, vec1)

    ! BRDF
    res = max(0.0_R_, pac(1) + pac(2) * fg + pac(3) * fv) * API

  end function mcarLsrt__BRDF_R


  !+
  ! Functions a_g & a_v of LSRT BRDF model
  !-
  subroutine mcarLsrt__albKernel(sinq0, cosq0, nthe, nphi, ag, av)

    real(R_), intent(in) :: sinq0, cosq0 ! sine & cosine of incoming zenith angle (downward motion to the plane)
    integer,  intent(in) :: nthe ! number of quadrature points for theta1
    integer,  intent(in) :: nphi ! number of quadrature points for phi
    real(R_), intent(out) :: ag, av ! the function a_g & a_v
    real(R_), parameter :: FAC = 2.0_R_ / PI_
    real(R_)  :: sumg, sumv, fbarg, fbarv
    integer   :: nqtmp, iq
    integer,  parameter :: KNQUAD = 200
    real(R_), save :: xq(KNQUAD), wq(KNQUAD), cosq1(KNQUAD), sinq1(KNQUAD)
    integer,  save :: nq = 0

    ! Initializations
    nqtmp = min(KNQUAD, nthe)
    if (nq /= nqtmp) then
       nq = nqtmp
       call gaussLegen(nq, xq, wq)
       do iq = 1, nq
          cosq1(iq) = 0.5_R_ * (1.0_R_ + xq(iq))
          sinq1(iq) = sqrt(1.0_R_ - cosq1(iq)**2)
       end do
    end if

    ! Integrations of f_bar_g & f_bar_v
    sumg = 0.0_R_
    sumv = 0.0_R_
    do iq = 1, nq
       call mcarLsrt__BRFKernel(sinq0, cosq0, sinq1(iq), cosq1(iq), nphi, fbarg, fbarv)
       sumg = sumg + wq(iq) * fbarg * cosq1(iq)
       sumv = sumv + wq(iq) * fbarv * cosq1(iq)
    end do
    ag = FAC * sumg
    av = FAC * sumv

  end subroutine mcarLsrt__albKernel


  !+
  ! Functions f_bar_g & f_bar_v of LSRT BRDF model
  !-
  subroutine mcarLsrt__BRFKernel(sinq0, cosq0, sinq1, cosq1, nphi, fbarg, fbarv)

    real(R_), intent(in) :: sinq0, cosq0 ! sine & cosine of incoming zenith angle (downward motion to the plane)
    real(R_), intent(in) :: sinq1, cosq1 ! sine &  cosine of outgoing zenith angle (upward motion from the plane)
    integer,  intent(in) :: nphi ! number of quadrature points for phi
    real(R_), intent(out) :: fbarg, fbarv ! the function f_bar_g & f_bar_v
    real(R_)  :: f1, sumg, sumv, vec0(3), vec1(3)
    integer   :: nqtmp, iq
    integer,  parameter :: KNQUAD = 200
    real(R_), save :: xq(KNQUAD), wq(KNQUAD), cosf1(KNQUAD), sinf1(KNQUAD)
    integer,  save :: nq = 0

    ! Initializations
    nqtmp = min(KNQUAD, nphi)
    if (nq /= nqtmp) then
       nq = nqtmp
       call gaussLegen(nq, xq, wq)
       do iq = 1, nq
          f1 = 0.5_R_ * PI_ * (1.0_R_ + xq(iq))
          cosf1(iq) = cos(f1)
          sinf1(iq) = sin(f1)
       end do
    end if

    ! Vector components
    vec0(1) = sinq0
    vec0(2) = 0.0_R_
    vec0(3) = cosq0

    ! Integrations of f_g & f_v
    sumg = 0.0_R_
    sumv = 0.0_R_
    do iq = 1, nq
       vec1(1) = sinq1 * cosf1(iq)
       vec1(2) = sinq1 * sinf1(iq)
       vec1(3) = cosq1
       sumg = sumg + wq(iq) * mcarLsrt__funcFg_R(vec0, vec1)
       sumv = sumv + wq(iq) * mcarLsrt__funcFv_R(vec0, vec1)
    end do
    fbarg = PI_ * sumg
    fbarv = PI_ * sumv

  end subroutine mcarLsrt__BRFKernel


  !+
  ! The function f_g of LSRT model
  !-
  function mcarLsrt__funcFg_R(vec0, vec1) result(res)

    real(R_), intent(in) :: vec0(:) ! incoming direction vector (downward motion to the plane)
    real(R_), intent(in) :: vec1(:) ! outgoing direction vector (upward motion from the plane)
    real(R_)  :: res
    real(R_),  parameter :: BR = 1.0_R_, HB = 2.0_R_
    real(R_),  parameter :: SBR = BR * BR
    real(R_),  parameter :: API = 1.0_R_ / PI_
    real(R_)  :: fnorm, vec0p(size(vec0)), vec1p(size(vec1)), cosgp, rnum, deno
    real(R_)  :: cost, sint, t

    ! Prime vectors
    if (BR == 1.0_R_) then ! this is the case in the current implementation
       vec0p(3) = vec0(3)
       vec1p(3) = vec1(3)
       cosgp = -geo3D_scaProd_R(vec0, vec1)
    else
       fnorm = 1.0_R_ / sqrt(SBR + (1.0_R_ - SBR) * vec0(3)**2)
       vec0p(1:2) = vec0(1:2) * fnorm * BR
       vec0p(3)   = vec0(3)   * fnorm
       fnorm = 1.0_R_ / sqrt(SBR + (1.0_R_ - SBR) * vec1(3)**2)
       vec1p(1:2) = vec1(1:2) * fnorm * BR
       vec1p(3)   = vec1(3)   * fnorm
       cosgp = geo3D_scaProd_R(vec0p, vec1p)
    end if
    cosgp = max(-1.0_R_, min(1.0_R_, cosgp)) ! cos(gamma')

    ! cost, sint, & t
    rnum = HB * sqrt(1.0_R_ - cosgp**2)
    deno = vec1p(3) - vec0p(3)   ! > or = 0
    if (rnum >= deno) then
       cost = 1.0_R_
       sint = 0.0_R_
       t = 0.0_R_
    else 
       cost = rnum / deno
       sint = sqrt(deno * deno - rnum * rnum) / deno
       t = acos(max(-1.0_R_, min(1.0_R_, cost))) ! mtab_acosX_R(cost)
    end if

    ! The function f_g
    rnum = (1.0_R_ - (t - sint * cost) * API) * deno - 0.5_R_ * (1.0_R_ + cosgp)
    deno = vec1p(3) * vec0p(3) ! overwritten
    if (abs(deno) > abs(1.0e-3_R_ * rnum)) then
       res = rnum / deno
    else if (deno * rnum >= 0.0_R_) then
       res =  1.0e+3_R_
    else
       res = -1.0e+3_R_
    end if

  end function mcarLsrt__funcFg_R


  !+
  ! The function f_v of LSRT model
  !-
  function mcarLsrt__funcFv_R(vec0, vec1) result(res)

    real(R_), intent(in) :: vec0(:) ! incoming direction vector (downward motion to the plane)
    real(R_), intent(in) :: vec1(:) ! outgoing direction vector (upward motion from the plane)
    real(R_)  :: res
    real(R_)  :: cosg, sing, g, rnum, deno

    ! cosg, sing, & g
    cosg = geo3D_scaProd_R(vec0, vec1)
    cosg = max(-1.0_R_, min(1.0_R_, -cosg))
    sing = sqrt(1.0_R_ - cosg**2)
    g = acos(cosg) !mtab_acosX_R(cosg)

    ! The function f_v
    rnum = (0.5_R_ * PI_ - g) * cosg + sing
    deno = vec1(3) - vec0(3)       ! > or = 0
    if (deno > 1.0e-3_R_ * rnum) then
       res = rnum / deno - 0.25_R_ * PI_
    else if (deno * rnum >= 0.0_R_) then
       res =  1.0e+3_R_
    else
       res = -1.0e+3_R_
    end if

  end function mcarLsrt__funcFv_R


  !+
  ! Get albedo and coefficients a and b of the RPV BRDF model
  !-
  subroutine mcarLsrt__prep_makeTab(pac, bmin, bsex, uz0min, nuz0, nuz1, &
       &     uz1wrk, agwrk, avwrk, fgwrk, fvwrk, cawrk, cbwrk, albmin, iam)

    real(R_), intent(in) :: pac(:) ! pac(1-3) : parameter packet
    !     pac(1) : "k_L", weight of Lambertian reflection
    !     pac(2) : "k_g", weight of geometrical optics reflection
    !     pac(3) : "k_v", weight of volume scattering
    real(R_), intent(in) :: uz1wrk(:)
    real(R_), intent(in) :: agwrk(:), avwrk(:)
    real(R_), intent(in) :: fgwrk(:, :), fvwrk(:, :)
    real(R_), intent(in) :: bmin ! minimum threshold for factor b' (about 3 is recommended)
    real(R_), intent(in) :: bsex ! scaling exponent for b' (should be < or = 1)
    !            Note: If s=1, higher accuracy of angular distribution of
    !            reflection will be achieved; about 0.5_R_ is recommended.
    integer,  intent(in) :: nuz0 ! number of cosines of incoming zenith angle (should be > 3)
    integer,  intent(in) :: nuz1 ! number of grid points for theta for fitting
    !            This is also used as number of grid points for finding "b"
    real(R_), intent(in) :: uz0min ! min cosines of incoming zenith angle
    real(R_), intent(out) :: cawrk(:) ! coeffficient a1
    real(R_), intent(out) :: cbwrk(:) ! coeffficient b'
    real(R_), intent(out) :: albmin ! minimum albedo in the table
    integer,  intent(out) :: iam ! iuz0 corresponding to the minimum albedo
    real(R_),  parameter :: A1MIN = -1.0_R_ / PI_ * 0.5_R_  ! min of the slope
    real(R_),  parameter :: A1MAX =  1.0_R_ / PI_ * 0.99_R_ ! max of the slope
    real(R_)  :: a0, a1, aa0, aa1, alb, b, bb, buz1, duz0, fnorm
    real(R_)  :: r, stt, sty, sx, sy, vec0(3), vec1(3), x, y
    integer   :: igrid, iuz0, iuz1, ngrid

    ! Macros
    duz0 = (1.0_R_ - uz0min) / (nuz0 - 1)
    ngrid = max(nuz1, 20)     ! # of points to find "b"

    ! Loop for various uz0
    albmin = 1.0_R_
    iam = 1
    do iuz0 = 1, nuz0

       ! Linear regression
       fnorm = 1.0_R_ / real(nuz1, R_)
       sx  = 0.0_R_
       do iuz1 = 1, nuz1
          sx  = sx  + uz1wrk(iuz1)
       end do
       sx = sx * fnorm
       sty = 0.0_R_
       sy  = 0.0_R_
       stt = 0.0_R_
       do iuz1 = 1, nuz1
          x = uz1wrk(iuz1)
          y = max(0.0_R_, pac(1) * PI_ + pac(2) * fgwrk(iuz1, iuz0) + pac(3) * fvwrk(iuz1, iuz0)) * x
          sy  = sy  + y
          sty = sty + y * (x - sx)
          stt = stt + (x - sx)**2
       end do
       sy  = sy  * fnorm
       sty = sty * fnorm
       stt = stt * fnorm
       aa1 = sty / stt
       aa0 = sy - sx * aa1
       !v1 = 1.0_R_ / stt          ! uncertainty of aa1
       !v0 = 1.0_R_ + sx * sx * v1 ! uncertainty of aa0
       a1 = aa1 / (PI_ * (aa1 + 2.0_R_ * aa0))
       a1 = max(A1MIN, min(A1MAX, a1))
       a0 = 0.5_R_ * (1.0_R_ / PI_ - a1)

       ! Find the factor b
       alb = pac(1) + pac(2) * agwrk(iuz0) + pac(3) * avwrk(iuz0)
       alb = max(1.0e-14_R_, alb) ! reject negative albedo
       vec0(3) = -min(1.0_R_, uz0min + duz0 * (iuz0 - 1))
       vec0(2) = 0.0_R_
       vec0(1) = sqrt(1.0_R_ - vec0(3)**2)
       vec1(:) = -vec0(:) ! hot spot
       b = vec1(3) / ((a1 * vec1(3) + a0) * alb) * mcarLsrt__BRDF_R(pac, vec0, vec1)
       vec1(3) = 1.0_R_   ! zenith
       vec1(1) = REPS_
       b = max(b, 1.0_R_ / ((a1 + a0) * alb) * mcarLsrt__BRDF_R(pac, vec0, vec1))
       buz1 = 1.0_R_ / ngrid
       do igrid = 1, ngrid - 1      ! other discrete points
          vec1(3) = 1.0_R_ - buz1 * igrid
          vec1(1) = -sqrt(1.0_R_ - vec1(3)**2)
          r = mcarLsrt__BRDF_R(pac, vec0, vec1)
          vec1(1) = -vec1(1)
          r = max(r, mcarLsrt__BRDF_R(pac, vec0, vec1))
          b = max(b, r * vec1(3) / ((a1 * vec1(3) + a0) * alb))
       end do
       b = b + 0.01_R_
       if (b < bmin) then
          bb = b
       else
          bb = bmin - 1.0_R_ + (b - bmin + 1.0_R_)**bsex
       end if

       ! Results
       cawrk(iuz0) = a1
       cbwrk(iuz0) = bb
       if (alb < albmin) then
          albmin = alb
          iam = iuz0
       end if
    end do

  end subroutine mcarLsrt__prep_makeTab


  !+
  ! Scattering of LSRT model : determination of a new direction
  ! Algorithm
  !   - The rejection method is used with a comparison function of linear function.
  !-
  subroutine mcarLsrt__newDir(vec0, pac, alb, a1, bb, vec1, pps)

    real(R_), intent(in) :: vec0(:) ! incoming direction vector (downward motion to the plane)
    real(R_), intent(in) :: alb     ! black-sky albedo (can be > 1)
    real(R_), intent(in) :: a1, bb  ! coefficients a1 & b'
    real(R_), intent(in) :: pac(:)  !  parameter packet
    !    pac(1) : "k_L", weight of Lambertian reflection
    !    pac(2) : "k_g", weight of geometrical optics reflection
    !    pac(3) : "k_v", weight of volume scattering
    real(R_), intent(out) :: vec1(:) ! outgoing direction vector (upward motion from the plane)
    real(R_), intent(out) :: pps     ! angular PDF (/steradian)
    real(R_),  parameter :: API = 1.0_R_ / PI_
    real(R_),  parameter :: R2MIN = REPS_**2
    real(R_)  :: a0, r2, sinf, cosf, d, vz, rnum, deno
    integer   :: itrial

    ! The rejection method
    a0 = 0.5_R_ * (API - a1)
    do itrial = 1, 100

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

       ! Rejection?
       rnum = vec1(3) * mcarLsrt__BRDF_R(pac, vec0, vec1)
       deno = bb * (a1 * vec1(3) + a0) * alb
       if (rnum > deno) exit
       if (mseq_rand_R() * deno <= rnum) exit
    end do

    ! A new direction & angular PDF
    call geo3D_UVec_renorm(vec1) ! correction of the unit vector
    pps = rnum / alb

  end subroutine mcarLsrt__newDir

end module mcarLsrt
