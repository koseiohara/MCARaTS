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
! Module for Diffuse-Specular-Mixture (DSM) BRDF model
! Reference:
!   Nakajima, T., & M. Tanaka, 1983: Effect of wind-generated waves on
!   the transfer of solar radiation in the atmosphere-ocean system.
!   J. Quant. Spectrosc. Radiat. Transfer, 29, 521-537.
!-
module mcarDsm

  use globals
  use hparx
  use mcarUtl
  implicit none
  private

  ! Public
  public :: mcarDsm__BRDFpA_R
  public :: mcarDsm__fshad_R
  public :: mcarDsm__funcA_R
  public :: mcarDsm__funcB_R
  public :: mcarDsm__prep_makeTab
  public :: mcarDsm__newDir_diff
  public :: mcarDsm__newDir_spec
  public :: mcarDsm__specVec

  ! Private
  real(R_), parameter :: Dsm_VARMIN = 1.0e-4_R_

contains

  !+
  ! BRDF/albedo of DSM model
  !-
  function mcarDsm__BRDFpA_R(pac, alb, vec0, vec1) result(res)

    real(R_), intent(in) :: pac(:) ! parameter packet
    !//  pac(1) : "alpha_d", albedo for the diffuse reflection
    !    pac(2) : "f_d", fraction of the diffuse reflection
    !    pac(3) : "r_r", real part of refractive index
    !    pac(4) : "r_i", imaginary part of refractive index
    !    pac(5) : "sigma^2", surface roughness, variance of tangent
    real(R_), intent(in) :: alb     ! average albedo
    real(R_), intent(in) :: vec0(:) ! incoming direction vector (downward motion to the plane)
    real(R_), intent(in) :: vec1(:) ! outgoing direction vector (upward motion from the plane)
    real(R_)  :: res ! result, BRDF/albedo
    real(R_),  parameter :: API = 1.0_R_ / PI_
    real(R_),  parameter :: PMAX = 1.0e+8_R_
    real(R_)  :: deno, fshad, fshad0, funcs
    real(R_)  :: rho, rnum, rpad, rpas, uzn2, vecn(size(vec0)), v !(AUTO)

    ! Purely diffusive model
    if (pac(1) >= 1.0_R_) then
       res = API

       ! Invalid directions
    else if (vec0(3) * vec1(3) > -REPS_) then
       res = pac(2) * pac(1) * API / alb

       ! Mixed model
    else
       v = max(Dsm_VARMIN, pac(5))
       fshad0 = mcarDsm__fshad_R(-vec0(3), v)
       call mcarDsm__specVec(pac, vec0, vec1, fshad0, vecn, rho, fshad, funcs)
       uzn2 = vecn(3)**2
       rnum = funcs * API / v * exp( -(1.0_R_ - uzn2) / (v * uzn2) )
       deno = 4.0_R_ * (-vec0(3)) * vec1(3) * alb
       if (deno * PMAX > rnum) then ! to avoid too large values
          rpas = rnum / deno
       else
          rpas = PMAX
       end if
       rpad = pac(1) * API / alb ! diffuse component
       res = pac(2) * rpad + (1.0_R_ - pac(2)) * rpas
    end if

  end function mcarDsm__BRDFpA_R


  !+
  ! Function F for the shadowing factor of rough micro-facet model
  !-
  function mcarDsm__fshad_R(uz, vari) result(res)

    real(R_), intent(in) :: uz
    real(R_), intent(in) :: vari
    real(R_)  :: res
    real(R_),  parameter :: FPI = 0.56418958_R_
    real(R_)  :: t, v

    v = abs(uz) / max(REPS_, sqrt(vari * max(0.0_R_, 1.0_R_ - uz * uz)))
    v = max(REPS_, v)

    if (v < 10.0_R_) then
       !t = 0.5_R_ * (exp(-v * v) / v**2 * fpi - erfc_R(v)) ! incorrect eq.
       t = 0.5_R_ * (exp(-v * v) / v * FPI - erfc_R(v)) ! non-standard Fortran (erfc)
       !// 01/11/2008, debug
       res = max(0.0_R_, t)
    else
       res = 0.0_R_
    end if

  end function mcarDsm__fshad_R


  !+
  ! Function A of the DSM BRDF model
  !-
  function mcarDsm__funcA_R(pac, fshad0, sinq0, cosq0, nthe, nphi) result(res)

    real(R_), intent(in) :: cosq0, sinq0 ! sine/cosine of incoming zenith angle (both positive)
    real(R_), intent(in) :: fshad0 ! shadowing function F for incoming direction
    integer,  intent(in) :: nphi   ! number of quadrature points for phi
    integer,  intent(in) :: nthe   ! number of quadrature points for theta
    real(R_), intent(in) :: pac(:) ! pac(1-5) : parameter packet (only pac(3-5) are used)
    !    pac(3) : "r_r", real(R_)  :: part of refractive index
    !    pac(4) : "r_i", imaginary part of refractive index
    !    pac(5) : "sigma^2", surface roughness real(R_),  parameter :: variance of tangent
    real(R_)  :: res

    real(R_),  parameter :: CMAX = 3.0_R_
    integer,   parameter :: KNQUAD = 500
    real(R_), save :: x(KNQUAD), w(KNQUAD)
    integer,  save :: nx = 0

    real(R_)  :: cosqq, cv, d, dlt
    integer   :: ix
    integer   :: n
    real(R_)  :: rho, rho1, rho2, sinqq, tsum, uzt

    ! Initializations
    n = min(KNQUAD, nthe)
    if (nx /= n) then
       nx = n
       call gaussLegen(nx, x, w)
       do ix = 1, nx
          x(ix) = 0.5_R_ * (1.0_R_ + x(ix))
       end do                 ! x is set as 0 to 1
    end if

    ! Flat surface
    if (pac(5) < Dsm_VARMIN) then
       call fresnelRef1(pac(3), pac(4), cosq0, uzt, rho1, rho2, rho)
       res = rho * cosq0

       ! Rough surface
    else
       cv = CMAX**2 * pac(5)
       dlt = 2.0_R_ * cv / (1.0_R_ + cv)
       tsum = 0.0_R_
       do ix = nx, 1, -1      ! reverse order to avoid rounding error
          d = dlt * x(ix)
          sinqq = sqrt(d * (2.0_R_ - d))
          cosqq = 1.0_R_ - d
          tsum = tsum + w(ix) * mcarDsm__funcB_R(pac, fshad0, sinq0, cosq0, sinqq, cosqq, nphi)
       end do
       res = tsum * dlt
    end if

  end function mcarDsm__funcA_R


  !+
  ! Function B of the DSM BRDF model
  !-
  function mcarDsm__funcB_R(pac, fshad0, sinq0, cosq0, sinqq, cosqq, nphi) result(res)

    real(R_), intent(in) :: cosq0, sinq0 ! sine/cosine of incoming zenith angle (both positive)
    real(R_), intent(in) :: cosqq, sinqq ! sine/cosine of angle from the specular direction
    real(R_), intent(in) :: fshad0 ! shadowing function F for incoming direction
    integer,  intent(in) :: nphi   ! number of quadrature points for phi
    real(R_), intent(in) :: pac(:) ! parameter packet (only pac(3-5) are used)
    !    pac(3) : "r_r", real(R_)  :: part of refractive index
    !    pac(4) : "r_i", imaginary part of refractive index
    !    pac(5) : "sigma^2", surface roughness real(R_),  parameter :: variance of tangent
    !              the variance should be > 0
    real(R_)  :: res
    integer,   parameter :: KNQUAD = 500
    real(R_), save :: w(KNQUAD), cosff(KNQUAD), sinff(KNQUAD)
    integer,  save :: nx = 0
    real(R_)  :: cc, cosffmin, ff, fs
    integer   :: ix, ixs
    integer   :: n
    real(R_)  :: pmod, rho, rho1, rho2, squzb, squzn, ss, tsum, tmp
    real(R_)  :: vec0(3), vec1(3), uzb, uzt
    real(R_)  :: y

    ! Initializations
    n = min(KNQUAD, nphi)
    if (nx /= n) then
       nx = n
       call gaussLegen(nx, cosff, w) ! cosff is used as X
       do ix = 1, nx
          ff = 0.5_R_ * PI_ * (1.0_R_ + cosff(ix)) ! 0 to pi
          cosff(ix) = cos(ff) ! decreasing
          sinff(ix) = sin(ff)
       end do
    end if
    vec0(1) = sinq0
    vec0(2) = 0.0_R_
    vec0(3) = cosq0

    ! Min rotational azimuth angle
    cc = cosq0 * cosqq
    ss = sinq0 * sinqq
    if (cc >= ss) then
       cosffmin = 1.0_R_
    else if (cc <= -ss) then ! unexpected case?
       res = 0.0_R_
       return
    else
       cosffmin = cc / ss
    end if

    ! Initial point (to reject the downward hemisphere)
    ixs = 1 + gridIdx_bin_I(cosff, 1, nx, cosffmin)

    ! Integration over the rotational azimuth
    tsum = 0.0_R_
    do ix = ixs, nx
       vec1(1:2) = -vec0(1:2)
       vec1(3)   =  vec0(3)
       vec1(:) = geo3D_rotateU_1R(vec1, sinqq, cosqq, sinff(ix), cosff(ix))
       squzb = 0.5_R_ * (1.0_R_ + geo3D_scaProd_R(vec0, vec1)) ! = (1+cosX)/2
       uzb = sqrt(squzb)
       call fresnelRef1(pac(3), pac(4), uzb, uzt, rho1, rho2, rho)
       fs = 1.0_R_ / (1.0_R_ + fshad0 + mcarDsm__fshad_R(vec1(3), pac(5)))
       squzn = (vec0(3) + vec1(3))**2 / squzb * 0.25_R_
       tmp = (1.0_R_ - squzn) / (pac(5)*squzn)
       pmod = mtab_expNX_R(tmp)
       y = rho * fs * pmod / (squzn * squzn)
       tsum = tsum + w(ix) * y
    end do
    res = tsum / (2.0_R_ * pac(5))
    !// Note: This algorithm is not very good for accuracy.
    !   The integration should be done for the range cos(phi) = [-1,cosffmin], 
    !   and the quadrature points should be set depending on (q & q0).
    !   It would be needed to compute sin(phi) & cos(phi) every time this subroutine is called.

  end function mcarDsm__funcB_R


  !+
  ! Get albedo and coefficients a and b of the DSM BRDF model
  !-
  subroutine mcarDsm__prep_makeTab(pac, r2min, srmin, srexp, uz0min, nuz0, &
       &     nthe, nphi, falut, sxlut, albmin, iam)

    integer,  intent(in) :: nphi   ! number of quadrature points for phi
    integer,  intent(in) :: nthe   ! number of quadrature points for theta
    !          This is also used as number of grid points for finding Smax
    real(R_), intent(in) :: r2min  ! minimum threshold of r^2 to determine the max tangent
    real(R_), intent(in) :: srexp  ! scaling exponent for "S" ratio (should be < or = 1)
    !          Note: If this is 1, higher accuracy of angular distribution of
    !          reflection will be achieved; about 0.5_R_ is recommended.
    real(R_), intent(in) :: srmin  ! minimum threshold for "S" ratio (about 4 is recommended)
    real(R_), intent(in) :: uz0min ! minimum of cosine of incident zenith angle
    real(R_), intent(in) :: pac(:) ! pac(1-5) : parameter packet
    !    pac(1) : "alpha_d", albedo for the diffuse reflection
    !    pac(2) : "f_d", fraction of the diffuse reflection
    !    pac(3) : "r_r", real(R_)  :: part of refractive index
    !    pac(4) : "r_i", imaginary part of refractive index
    !    pac(5) : "sigma^2", surface roughness real(R_),  parameter :: variance of tangent
    integer,  intent(inout) :: nuz0   ! number of cosines of incoming zenith angle (should be > 3)
    real(R_), intent(out) :: falut(:) ! LUT of function A for iuz0=1-nuz0
    real(R_), intent(out) :: sxlut(:) ! LUT of coefficient Smax for iuz0=1-nuz0
    real(R_), intent(out) :: albmin ! minimum albedo in the table
    integer,  intent(out) :: iam    ! iuz0 corresponding to the minimum albedo
    real(R_),  parameter :: UEPS = 0.001_R_
    integer   :: ithe, iuz0, ngrid
    real(R_)  :: alb, albd, buz0, cc, cc2, cosq0, cost, fs, fshad, fshad0, funca, funcs
    real(R_)  :: rho, rho1, rho2, rr, s0, sinq0, sint, smax, ss, ss2, u, ub, ue, us
    real(R_)  :: vec0(3), vec1(3), vecn(3), uzt

    ! Macros
    nuz0 = max(4, nuz0)
    buz0 = (1.0_R_ - uz0min) / (nuz0 - 1)
    ngrid = max(nthe, 8)      ! # of points to find Smax
    rr = -pac(5) * log(r2min) ! tan(q_max)
    cc = 1.0_R_ / sqrt(1.0_R_ + rr) ! cos(q_max)
    ss = cc * sqrt(rr)        ! sin(q_max)
    cc2 = 1.0_R_ - 2.0_R_ * ss * ss ! cos(2 * q_max)
    ss2 = 2.0_R_ * ss * cc       ! sin(2 * q_max)
    albd = pac(2) * pac(1)
    fs = 1.0_R_ - pac(2)

    ! Loop for various uz0
    albmin = 1.0_R_
    iam = 1      
    do iuz0 = 1, nuz0
       cosq0 = min(1.0_R_, uz0min + buz0 * (iuz0 - 1))
       sinq0 = sqrt(1.0_R_ - cosq0 * cosq0)
       vec0(1) = -sinq0 ! incident vector (downward)
       vec0(2) = 0.0_R_
       vec0(3) = -cosq0

       ! Flat case
       if (pac(5) < Dsm_VARMIN) then
          call fresnelRef1(pac(3), pac(4), cosq0, uzt, rho1, rho2, rho)
          funca = rho * cosq0
          smax = RSML_

          ! Rough case
       else
          ! Function A
          fshad0 = mcarDsm__fshad_R(cosq0, pac(5))
          funca = mcarDsm__funcA_R(pac, fshad0, sinq0, cosq0, nthe, nphi)
          funca = max(1.0e-20_R_, funca)
          vec1(1) =  vec0(1) ! to get S value when facet normal is zenith
          vec1(2) =  vec0(2)
          vec1(3) = -vec0(3)
          call mcarDsm__specVec(pac, vec0, vec1, fshad0, vecn, rho, fshad, funcs)
          s0 = max(1.0e-25_R_, funcs)
          smax = s0 ! initial estimate of max S

          ! Gridding
          sint = sinq0 * cc2 - cosq0 * ss2 ! = sin(q0 - 2*qmax)
          cost = cosq0 * cc2 + sinq0 * ss2 ! = cos(q0 - 2*qmax)
          cost = max(UEPS, cost)
          if (sint < 0.0_R_) then
             us = 2.0_R_ - cost
          else
             us = cost
          end if
          ue = max(UEPS, cosq0 * cc2 - sinq0 * ss2)
          ub = (ue - us) / real(ngrid, R_)

          ! Find Smax
          do ithe = 0, ngrid - 1 ! for discrete points
             u = us + ub * (ithe + 0.5_R_)
             if (u > 1.0_R_) then
                vec1(3) = 2.0_R_ - u
                vec1(1) = sqrt(1.0_R_ - vec1(3)**2)
             else
                vec1(3) = u
                vec1(1) = -sqrt(1.0_R_ - vec1(3)**2)
             end if
             call mcarDsm__specVec(pac, vec0, vec1, fshad0, vecn, rho, fshad, funcs)
             smax = max(smax, funcs)
          end do
          if (smax > s0 * srmin) smax = s0 * (srmin - 1.0_R_ + (smax/s0 - srmin + 1.0_R_)**srexp)
       end if

       ! Results
       falut(iuz0) = funca
       sxlut(iuz0) = smax * 1.2_R_ ! because smax is a rough estimate

       ! Min black-sky albedo
       alb = albd + fs * (funca / cosq0)
       if (alb < albmin) then
          albmin = alb
          iam = iuz0
       end if
    end do

  end subroutine mcarDsm__prep_makeTab


  !+
  ! New direction by diffusive Lambertian reflection
  !-
  subroutine mcarDsm__newDir_diff(vec1, pps)

    real(R_), intent(out) :: vec1(:) ! outgoing upward direction vector
    real(R_), intent(out) :: pps     ! angular PDF (/steradian)
    real(R_)  :: r2, sinf, cosf
    real(R_),  parameter :: API = 1.0_R_ / PI_

    call rand_point_circle(1.0e-12_R_, r2, sinf, cosf)
    vec1(:) = geo3D_UVec_LambertH_1R(r2, sinf, cosf)
    pps = abs(vec1(3)) * API

  end subroutine mcarDsm__newDir_diff


  !+
  ! Scattering of specular reflection : determination of a new direction
  !-
  subroutine mcarDsm__newDir_spec(rr, ri, var, smax, r2min, alb, vec0, vec1, pps)

    real(R_), intent(in) :: rr    ! real part of refractive index
    real(R_), intent(in) :: ri    ! imaginary part of refractive index
    real(R_), intent(in) :: smax  ! coefficient Smax (max of function S)
    real(R_), intent(in) :: var   ! "sigma^2", surface roughness parameter, variance of tangent
    real(R_), intent(in) :: r2min ! minimum threshold of r^2 to determine the max tangent
    real(R_), intent(in) :: alb   ! black-sky albedo
    real(R_), intent(in) :: vec0(:)  ! incoming downward direction vector
    real(R_), intent(out) :: vec1(:) ! outgoing upward direction vector
    real(R_), intent(out) :: pps     ! angular PDF (/steradian)
    real(R_)  :: cosb, cosf, fs, fshad0, funcs, rnum, deno, rpas, uzn2
    real(R_)  :: r2, rho, rho1, rho2, sinf, t2, vecn(3), uzt, vzn
    integer   :: itrial
    real(R_),  parameter :: PMAX = 1.0e+8_R_
    real(R_),  parameter :: API = 1.0_R_ / PI_

    ! Flat surface
    if (var < Dsm_VARMIN) then
       vec1(1:2) = vec0(1:2)
       vec1(3) = -vec0(3)
       pps = abs(vec1(3)) * PMAX
       return
    end if

    ! Rough surface: Rejection method
    fshad0 = mcarDsm__fshad_R(-vec0(3), var) ! shadowning function for the incoming direction
    do itrial = 1, 30 ! will be iterated only one or two times usually
       
       ! A new candidate vector
       do
          call rand_point_circle(r2min, r2, sinf, cosf)
          t2 = -var * log(r2) ! tan(q)**2, should be >= 0
          vecn(3) = 1.0_R_ / sqrt(1.0_R_ + t2) ! normal vector
          vzn = vecn(3) * sqrt(t2)
          vecn(1) = vzn * cosf
          vecn(2) = vzn * sinf
          cosb = -geo3D_scaProd_R(vec0, vecn) ! is sometimes negative
          vec1(3) = vec0(3) + 2.0_R_ * cosb * vecn(3) ! reflection vector
          if (cosb > 0.0_R_ .and. vec1(3) > REPS_) exit
       end do
       vec1(1:2) = vec0(1:2) + 2.0_R_ * cosb * vecn(1:2)

       ! Accept or reject
       call fresnelRef1(rr, ri, cosb, uzt, rho1, rho2, rho)
       fs = 1.0_R_ / (1.0_R_ + fshad0 + mcarDsm__fshad_R(vec1(3), var))
       uzn2 = vecn(3)**2
       funcs = rho * fs / uzn2**2 ! function S
       if (mseq_rand_R() * smax <= funcs) exit ! accept the candidate?
    end do

    ! A new direction and PDF
    call geo3D_UVec_renorm(vec1)
    rnum = funcs * API / var * exp( -(1.0_R_ - uzn2) / (var * uzn2) )
    deno = 4.0_R_ * (-vec0(3)) * vec1(3) * alb
    if (deno * PMAX > rnum) then ! to avoid too large values
       rpas = rnum / deno
    else
       rpas = PMAX
    end if
    pps = abs(vec1(3)) * rpas

  end subroutine mcarDsm__newDir_spec


  !+
  !* Specular reflection parameters: facet vector, Fresnell reflectance,
  !  & function S of DSM model from specified incoming/outgoing vectors *
  !// Caution! This subprogram does not accept vec0(3)=0, vec1(3)=0 or pac(5)=0
  !-
  subroutine mcarDsm__specVec(pac, vec0, vec1, fshad0, vecn, rho, fshad, funcs)

    real(R_), intent(in)  :: pac(:)  ! parameter packet (only pac(3-5) are used)
    !//  pac(3) : "r_r", real part of refractive index
    !    pac(4) : "r_i", imaginary part of refractive index
    !    pac(5) : "sigma^2", surface roughness, variance of tangent
    real(R_), intent(in)  :: vec0(:) ! incoming direction vector (downward motion to the plane)
    real(R_), intent(in)  :: vec1(:) ! outgoing direction vector (upward motion from the plane)
    real(R_), intent(in)  :: fshad0  ! shadowing function F for incoming direction
    real(R_), intent(out) :: vecn(:) ! facet normal vector (upward)
    real(R_), intent(out) :: rho     ! Fresnell reflectance
    real(R_), intent(out) :: fshad   ! shadowing factor for the bidirectional geometry
    real(R_), intent(out) :: funcs   ! function S
    real(R_)  :: rho1, rho2, uza, uzt

    ! Facet normal vector
    vecn(1:3) = vec1(1:3) - vec0(1:3)
    vecn(:) = geo3D_UVec_1R(vecn)

    ! Other functions
    uza = geo3D_scaProd_R(vecn, vec0)
    call fresnelRef1(pac(3), pac(4), -uza, uzt, rho1, rho2, rho) ! Fresnel reflectance
    fshad = 1.0_R_ / (1.0_R_ + fshad0 + mcarDsm__fshad_R(vec1(3), pac(5))) ! shadowing factor
    funcs = rho * fshad / vecn(3)**4 ! function S

  end subroutine mcarDsm__specVec

end module mcarDsm
