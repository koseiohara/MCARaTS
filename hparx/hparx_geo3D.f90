!-*-f90-*-
 
! @LICENSE_HEADER_START@
!
!   This file is part of HPARX.
!   
!   --
!   HPARX: Fortran code library for atmospheric sciences
!   
!   Copyright (C) 2006-2016
!   Hironobu Iwabuchi, Souichiro Hioki, and Rintaro Okamura
!   
!   HPARX is free software: you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation, either version 3 of the License, or
!   (at your option) any later version.
!   
!   HPARX is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!   GNU General Public License for more details.
!   
!   You should have received a copy of the GNU General Public License
!   along with HPARX. If not, see <http://www.gnu.org/licenses/>.
!
! @LICENSE_HEADER_END@
 
 
 
 
 
 


!+
! Library of utilities for 3-D geometrical calculations
!-
module hparx_geo3D

  use globals, only : R_, RD_, REPS_, RSML_, RLRG_, PI_
  use hparx_math
  use hparx_rand
  implicit none
  private

  ! Public
  public :: geo3D_angles
  public :: geo3D_aVec_1R
  public :: geo3D_norm_R
  public :: geo3D_randDir_cone_1R
  public :: geo3D_randLoc_circle_1R
  public :: geo3D_rayNPlane
  public :: geo3D_rayCone_L
  public :: geo3D_rotate_polar_1R
  public :: geo3D_rotate_axis_1R
  public :: geo3D_rotate_coord_1R
  public :: geo3D_rotateU_1R
  public :: geo3D_rotateU_perp_1R
  public :: geo3D_rotateU_ring_1R
  public :: geo3D_rotateU_orient_1R
  public :: geo3D_rotateU_orientH_1R
  public :: geo3D_rotMat_ZYZ_2R
  public :: geo3D_scaProd_R
  public :: geo3D_twoUVec
  public :: geo3D_UVec_1R
  public :: geo3D_UVec_isotropic_1R
  public :: geo3D_UVec_Lambert_1R
  public :: geo3D_UVec_LambertH_1R
  public :: geo3D_UVec_renorm
  public :: geo3D_vecProd_1R

contains

  !+
  ! Polar and azimuth angles of a direction
  !-
  subroutine geo3D_angles(dir, the, phi) 

    real(R_), intent(in) :: dir(:) ! unit direction vector
    real(R_), intent(out) :: the ! polar   angle (radian) [0, PI]
    real(R_), intent(out) :: phi ! azimuth angle (radian) [-PI, PI]
    the = atan2(sqrt(dir(1)**2 + dir(2)**2), dir(3))
    phi = atan2(dir(2), dir(1))
    
  end subroutine geo3D_angles


  !+
  ! A vector in Cartesian coordinates (right-hand system)
  !-
  function geo3D_aVec_1R(r, the, phi) result(vec)

    real(R_), intent(in) :: r   ! radius
    real(R_), intent(in) :: the ! polar angle
    real(R_), intent(in) :: phi ! azimuth angle
    real(R_) :: vec(3)
    real(R_) :: cosq, sinq

    cosq = cos(the)
    sinq = sin(the)
    vec(1) = r * sinq * cos(phi)
    vec(2) = r * sinq * sin(phi)
    vec(3) = r * cosq

  end function geo3D_aVec_1R


  !+
  ! Norm of the vector
  !-
  function geo3D_norm_R(a) result(s) 

    real(R_), intent(in)  :: a(:)    ! a 3-D vector A
    real(R_) :: s  ! the norm of A
    s = sqrt(sum(a(1:3)**2)) !a(1)**2 + a(2)**2 + a(3)**2)

  end function geo3D_norm_R


  !+
  ! A random direction within a conical angular region
  !-
  function geo3D_randDir_cone_1R(adir, cosdc, mproj) result(bdir)

    real(R_), intent(in) :: adir(:) ! conical center direction vector
    real(R_), intent(in) :: cosdc   ! 1 - cos(qc), where qc is the half cone angle
    integer,  intent(in) :: mproj   ! method for projection
    !// (0:isotropic, 1:projection to horizontal plane)
    real(R_) :: bdir(size(adir)) ! a new direction
    real(R_) :: cosf, cosf0, cosqc, cosq, u1, u2, umax, r2, sinf, sinf0, sinq, sinq0, cosd, sind
    integer  :: itrial

    ! Collimated emission
    if (cosdc < 1.0e-12_R_) then
       bdir(:) = adir(:)

       ! Random direction
    else
       sinq0 = sqrt(adir(1)**2 + adir(2)**2)

       ! Isotropic source
       if (mproj <= 0) then
          call rand_point_circle(1e-12_R_, r2, sinf, cosf)
          cosqc = r2 * cosdc
          cosq = 1.0_R_ - cosqc
          sinq = sqrt(cosqc * (2.0_R_ - cosqc))
          bdir(3) = cosq * adir(3) - sinq * cosf * sinq0
          
          ! Projection to horizontal plane
       else
          cosd = 1.0_R_ - cosdc
          sind = sqrt(cosdc * (2.0_R_ - cosdc))
          if (cosd < abs(adir(3))) then
             umax = 1.0_R_
          else
             u1 = adir(3) * cosd - sinq0 * sind ! cos(q0+dq)
             u2 = adir(3) * cosd + sinq0 * sind ! cos(q0-dq)
             umax = max(abs(u1), abs(u2))
          end if
          do itrial = 1, 100 ! loop for random directions (average # of interations is <= 2)
             call rand_point_circle(1e-12_R_, r2, sinf, cosf) ! a random direction
             cosqc = r2 * cosdc
             cosq = 1.0_R_ - cosqc
             sinq = sqrt(cosqc * (2.0_R_ - cosqc))
             bdir(3) = cosq * adir(3) - sinq * cosf * sinq0
             if (mseq_rand_R() * umax <= abs(bdir(3))) exit ! accept or reject
          end do
       end if
       
       ! Normalized new vector
       if (sinq0 > 1.0e-12_R_) then
          cosf0 = adir(1) / sinq0
          sinf0 = adir(2) / sinq0
          bdir(1) = cosq * adir(1) + sinq * (cosf * adir(3) * cosf0 - sinf * sinf0)
          bdir(2) = cosq * adir(2) + sinq * (cosf * adir(3) * sinf0 + sinf * cosf0)
       else
          bdir(1) = sinq * cosf
          bdir(2) = sinq * sinf
       end if
       bdir(:) = geo3D_UVec_1R(bdir) ! normalize
    end if

  end function geo3D_randDir_cone_1R


  !+
  ! A random location within an arbitrary oriented circle
  !-
  function geo3D_randLoc_circle_1R(aloc, eta, rc) result(bloc)

    real(R_), intent(in) :: aloc(:) ! the center point of the circle
    real(R_), intent(in) :: eta(:)  ! the normal direction of the circle
    real(R_), intent(in) :: rc      ! radius of the circle
    real(R_) :: bloc(size(aloc)) ! a new location
    real(R_) :: rn, sinf, cosf

    if (rc > 0.0_R_) then
       call rand_point_circle(1.0e-12_R_, rn, sinf, cosf)
       bloc(:) = aloc(:) + rc * sqrt(rn) * geo3D_rotateU_perp_1R(eta, sinf, cosf)
    else
       bloc(:) = aloc(:)
    end if

  end function geo3D_randLoc_circle_1R


  !+
  ! Test whether the ray has an intersection point with a coordinate-normal plane
  !-
  subroutine geo3D_rayNPlane(pmin, pmax, oloc1, cdir1, x, oloc2, cdir2, &
       & ymin, ymax, oloc3, cdir3, zmin, zmax, path, hit)

    real(R_), intent(in) :: pmin, pmax
    real(R_), intent(in) :: oloc1, cdir1
    real(R_), intent(in) :: x
    real(R_), intent(in) :: oloc2, cdir2
    real(R_), intent(in) :: ymin, ymax
    real(R_), intent(in) :: oloc3, cdir3
    real(R_), intent(in) :: zmin, zmax
    real(R_), intent(out) :: path
    logical,  intent(out) :: hit ! true if the ray hits the plane
    real(R_) :: vec2, vec3

    hit = .false.
    path = (x - oloc1) / cdir1
    if (path > 0.0_R_) then
       vec2 = oloc2 + path * cdir2
       vec3 = oloc3 + path * cdir3
       if ((ymin - vec2) * (ymax - vec2) <= 0.0_R_ .and. &
            & (zmin - vec3) * (zmax - vec3) <= 0.0_R_ .and. &
            & (pmin - path) * (pmax - path) <= 0.0_R_) then
          hit = .true.
       end if
    end if

  end subroutine geo3D_rayNPlane


  !+
  ! Test whether a line segment is within a conical region
  !-
  function geo3D_rayCone_L(rdir, dloc, smax, cdir, cosqmax) result(sc)

    real(R_), intent(in) :: rdir(:) ! ray direction
    real(R_), intent(in) :: dloc(:) ! location of the start point (C to R vector)
    real(R_), intent(in) :: smax    ! max transport distance
    real(R_), intent(in) :: cdir(:) ! cone center direction
    real(R_), intent(in) :: cosqmax ! cos(half cone angle)
    logical :: sc ! true if the ray has intersection
    real(R_) :: ar(3), ac(3), ur(3), uc(3), wrk(3)
    real(R_) :: a, b, c, s1, s2
    
    ! Test the start point
    if (geo3D_scaProd_R(cdir, dloc) >= geo3D_norm_R(dloc) * cosqmax) then
       sc = .true.
       return
    end if

    ! Coefficient, a
    ar(1) = rdir(1) * rdir(2)
    ar(2) = rdir(2) * rdir(3)
    ar(3) = rdir(3) * rdir(1)
    ac(1) = cdir(1) * cdir(2)
    ac(2) = cdir(2) * cdir(3)
    ac(3) = cdir(3) * cdir(1)
    ur(1:3) = rdir(1:3)**2
    uc(1:3) = cdir(1:3)**2
    a = geo3D_scaProd_R(ur, uc) + 2.0_R_ * geo3D_scaProd_R(ar, ac) - cosqmax**2

    ! Coefficients, b & c
    wrk(:) = dloc(:) * rdir(:)
    b = geo3D_scaProd_R(wrk, uc(:) - cosqmax**2)
    wrk(1) = dloc(1) * rdir(2) + dloc(2) * rdir(1)
    wrk(2) = dloc(2) * rdir(3) + dloc(3) * rdir(2)
    wrk(3) = dloc(3) * rdir(1) + dloc(1) * rdir(3)
    b = 2.0_R_ * (b + geo3D_scaProd_R(wrk, ac))
    wrk(1) = dloc(1) * dloc(2)
    wrk(2) = dloc(2) * dloc(3)
    wrk(3) = dloc(3) * dloc(1)
    c = geo3D_scaProd_R(dloc(:)**2, uc(:) - cosqmax**2) + 2.0_R_ * geo3D_scaProd_R(wrk, ac)

    ! Find roots
    call root_quad(a, b, c, s1, s2)
    sc = .false.
    if (s1 > 0.0_R_ .and. s1 < smax .and. geo3D_scaProd_R(cdir, dloc(:) + s1 * rdir(:)) > 0.0_R_) then
       sc = .true.
    else if (s2 > 0.0_R_ .and. s2 < smax .and. geo3D_scaProd_R(cdir, dloc(:) + s2 * rdir(:)) > 0.0_R_) then
       sc = .true.
    end if

  end function geo3D_rayCone_L


  !+
  ! Rotate a vector V from a polar vector P, rotating V about (P x V)
  !-
  function geo3D_rotate_polar_1R(vec, eta, sind, cosd) result(vec1)

    real(R_), intent(in) :: vec(:)  ! a vector V
    real(R_), intent(in) :: eta(:)  ! the polar vector P
    real(R_), intent(in) :: sind    ! sin(rotation angle)
    real(R_), intent(in) :: cosd    ! cos(rotation angle)
    real(R_) :: vec1(size(vec, 1))  ! rotated vector
    real(R_) :: a(3), aa, fac, a1s, a2s, a3s, aa1

    ! Vector product, P x V
    a(:) = geo3D_vecProd_1R(eta, vec)
    aa = geo3D_norm_R(a)

    ! Rotate V
    if (aa > RSML_) then
       fac = sind / aa
       a1s = a(1) * fac
       a2s = a(2) * fac
       a3s = a(3) * fac
       vec1(1) = cosd * vec(1) -  a3s * vec(2) +  a2s * vec(3)
       vec1(2) =  a3s * vec(1) + cosd * vec(2) -  a1s * vec(3)
       vec1(3) = -a2s * vec(1) +  a1s * vec(2) + cosd * vec(3)
       aa1 = geo3D_norm_R(vec1)
       if (aa1 > RSML_) vec1(1:3) = vec1(1:3) * geo3D_norm_R(vec) / aa1
    else
       vec1(1:3) = vec(1:3)
    end if

  end function geo3D_rotate_polar_1R


  !+
  ! Rotate a vector V about an axis vector A, where A is normalized
  !-
  function geo3D_rotate_axis_1R(vec, axis, sind, cosd) result(vec1)

    real(R_), intent(in) :: vec(:)  ! a vector V
    real(R_), intent(in) :: axis(:) ! the axis vector A (should be normalized!)
    real(R_), intent(in) :: sind    ! sin(rotation angle)
    real(R_), intent(in) :: cosd    ! cos(rotation angle)
    real(R_) :: vec1(size(vec, 1))  ! rotated vector
    real(R_) :: rot(3, 3), cosd_c, c1, c2, c3, aa1

    ! 1 - cos
    if (1.0_R_ + cosd < REPS_) then
       cosd_c = 2.0_R_
    else
       cosd_c = sind**2 / (1.0_R_ + cosd)
    end if

    ! Matrix
    c1 = axis(1) * axis(2) * cosd_c
    c2 = axis(2) * axis(3) * cosd_c
    c3 = axis(3) * axis(1) * cosd_c
    rot(1, 1) = cosd + cosd_c * axis(1)**2
    rot(2, 2) = cosd + cosd_c * axis(2)**2
    rot(3, 3) = cosd + cosd_c * axis(3)**2
    rot(2, 1) = c1 - axis(3) * sind
    rot(1, 2) = c1 + axis(3) * sind
    rot(3, 2) = c2 - axis(1) * sind
    rot(2, 3) = c2 + axis(1) * sind
    rot(1, 3) = c3 - axis(2) * sind
    rot(3, 1) = c3 + axis(2) * sind

    ! New vector
    vec1(1) = rot(1, 1) * vec(1) + rot(2, 1) * vec(2) + rot(3, 1) * vec(3)
    vec1(2) = rot(1, 2) * vec(1) + rot(2, 2) * vec(2) + rot(3, 2) * vec(3)
    vec1(3) = rot(1, 3) * vec(1) + rot(2, 3) * vec(2) + rot(3, 3) * vec(3)
    aa1 = geo3D_norm_R(vec1)
    if (aa1 > RSML_) vec1(1:3) = vec1(1:3) * geo3D_norm_R(vec) / aa1

  end function geo3D_rotate_axis_1R


  !+
  ! Rotate a vector V about a coordinate axis X/Y/Z
  !-
  function geo3D_rotate_coord_1R(vec, maxis, sind, cosd) result(vec1)

    real(R_), intent(in) :: vec(:)  ! a vector V
    integer,  intent(in) :: maxis   ! axis, 1:X, 2:Y, 3:Z
    real(R_), intent(in) :: sind    ! sin(rotation angle)
    real(R_), intent(in) :: cosd    ! cos(rotation angle)
    real(R_) :: vec1(size(vec, 1))  ! rotated vector

    ! About X-axis
    if (maxis == 1) then
       vec1(1) = vec(1)
       vec1(2) = cosd * vec(2) - sind * vec(3)
       vec1(3) = sind * vec(2) + cosd * vec(3)

       ! About Y-axis
    else if (maxis == 2) then
       vec1(1) = sind * vec(3) + cosd * vec(1)
       vec1(2) = vec(2)
       vec1(3) = cosd * vec(3) - sind * vec(1)

       ! About Z-axis
    else if (maxis == 3) then
       vec1(1) = cosd * vec(1) - sind * vec(2)
       vec1(2) = sind * vec(1) + cosd * vec(2)
       vec1(3) = vec(3)
    end if

  end function geo3D_rotate_coord_1R


  !+
  ! Rotate a unit vector
  !-
  function geo3D_rotateU_1R(uvec, sinq, cosq, sinf, cosf) result(uvec1)

    real(R_), intent(in) :: uvec(:)     ! unit vector
    real(R_), intent(in) :: sinq, cosq  ! sin & cos of theta
    real(R_), intent(in) :: sinf, cosf  ! sin & cos of phi
    real(R_) :: uvec1(size(uvec, 1))    ! rotated vector
    real(R_) :: sinq02, sinq0, cosf0, sinf0

    sinq02 = uvec(1)**2 + uvec(2)**2

    if (sinq02 > REPS_**2) then
       sinq0 = sqrt(sinq02)
       cosf0 = uvec(1) / sinq0
       sinf0 = uvec(2) / sinq0
       uvec1(1) = cosq * uvec(1) + sinq * (cosf * uvec(3) * cosf0 - sinf * sinf0)
       uvec1(2) = cosq * uvec(2) + sinq * (cosf * uvec(3) * sinf0 + sinf * cosf0)
       uvec1(3) = cosq * uvec(3) - sinq * cosf * sinq0
    else if (uvec(3) > 0.0_R_) then
       uvec1(1) = sinq * cosf
       uvec1(2) = sinq * sinf
       uvec1(3) = cosq
    else
       uvec1(1) = -sinq * cosf
       uvec1(2) = -sinq * sinf
       uvec1(3) = -cosq
    end if

    call geo3D_UVec_renorm(uvec1) ! renormalize

  end function geo3D_rotateU_1R


  !+
  ! Rotate a unit vector to a perpendicular direction to the original
  !-
  function geo3D_rotateU_perp_1R(uvec, sinf, cosf) result(uvec1)

    real(R_), intent(in) :: uvec(:)     ! unit vector
    real(R_), intent(in) :: sinf, cosf  ! sin & cos of phi
    real(R_) :: uvec1(size(uvec, 1))    ! rotated vector
    real(R_) :: sinq02, sinq0

    sinq02 = uvec(1)**2 + uvec(2)**2
    if (sinq02 > REPS_**2) then
       sinq0 = sqrt(sinq02)
       uvec1(1) = (cosf * uvec(3) * uvec(1) - sinf * uvec(2)) / sinq0
       uvec1(2) = (cosf * uvec(3) * uvec(2) + sinf * uvec(1)) / sinq0
       uvec1(3) = -cosf * sinq0
    else if (uvec(3) > 0.0_R_) then
       uvec1(1) = cosf
       uvec1(2) = sinf
       uvec1(3) = 0.0_R_
    else
       uvec1(1) = -cosf
       uvec1(2) = -sinf
       uvec1(3) = 0.0_R_
    end if
    call geo3D_UVec_renorm(uvec1) ! renormalize

  end function geo3D_rotateU_perp_1R


  !+
  ! Rotate a unit vector with a triangular distribution for rotation polar angle
  !-
  function geo3D_rotateU_orient_1R(vec0, xi, sinf, cosf, cosqs_c) result(vec1)

    real(R_), intent(in) :: vec0(:)    ! an input vector
    real(R_), intent(in) :: xi         ! a random number to determine the polar angle
    real(R_), intent(in) :: sinf, cosf ! sin & cos of azimuth angle
    real(R_), intent(in) :: cosqs_c    ! 1-cosQs, half cone size (Qs < pi/2)
    real(R_) :: vec1(size(vec0)) ! result, output vector
    real(R_) :: cosq1_c, sinq1, cosq1

    cosq1_c = cosqs_c * (1.0_R_ - sqrt(xi))   ! triangular
    cosq1 = 1.0_R_ - cosq1_c
    sinq1 = sqrt(cosq1_c * (2.0_R_ - cosq1_c))      
    vec1(:) = geo3D_rotateU_1R(vec0, sinq1, cosq1, sinf, cosf)

  end function geo3D_rotateU_orient_1R


  !+
  !* Rotate a unit vector with a triangular distribution for rotation polar angle,
  !  modified version with a output vector in a specific hemisphere *
  !-
  function geo3D_rotateU_orientH_1R(vec0, xi, sinf, cosf, cosqs_c, eta) result(vec1)

    real(R_), intent(in) :: vec0(:)    ! an input vector
    real(R_), intent(in) :: xi         ! a random number to determine the polar angle
    real(R_), intent(in) :: sinf, cosf ! sin & cos of azimuth angle
    real(R_), intent(in) :: cosqs_c    ! 1-cosQs, half cone size
    real(R_), intent(in) :: eta(:)     ! a normal vector that specify the hemisphere
    real(R_) :: vec1(size(vec0)) ! result, output vector 
    real(R_) :: sinq, cosq, sind, cosd, cosq1, cosq1_c, sinq1

    ! Angles
    cosq1_c = cosqs_c * (1.0_R_ - sqrt(xi))   ! triangular
    cosq1 = 1.0_R_ - cosq1_c
    sinq1 = sqrt(cosq1_c * (2.0_R_ - cosq1_c))      

    ! Rotate
    vec1(:) = vec0(:)
    call geo3D_twoUVec(eta, vec1, cosq, sinq)
    if (cosq < 0.0_R_) vec1(:) = geo3D_rotate_polar_1R(vec1, eta, cosq, sinq) ! sinD=cosq, cosD=sinq
    vec1(:) = geo3D_rotateU_1R(vec1, sinq1, cosq1, sinf, cosf)

    ! If in ill hemisphere, then fix
    call geo3D_twoUVec(eta, vec1, cosq, sinq)
    if (cosq < 0.0_R_) then
       sind = 2.0_R_ * sinq * cosq  
       cosd = 1.0_R_ - 2.0_R_ * cosq**2
       vec1(:) = geo3D_rotate_polar_1R(vec1, eta, sind, cosd)
    end if

  end function geo3D_rotateU_orientH_1R


  !+
  ! Determine a scattering direction in a Equatorial ring around the local normal vector
  !-
  function geo3D_rotateU_ring_1R(vec0, xi, sinf, cosf, cosqe, reps, mhem) result(vec1)

    real(R_), intent(in) :: vec0(:)    ! an input vector
    real(R_), intent(in) :: xi         ! a random number to determine the polar angle
    real(R_), intent(in) :: sinf, cosf ! sin & cos of azimuth angle
    real(R_), intent(in) :: cosqe      ! cosqe, half width of the ring
    real(R_), intent(in) :: reps       ! epsilon to avoid complete horizontal directions
    integer,  intent(in) :: mhem       ! method (0:both hemispheres, 1:forward hemisphere)
    real(R_) :: vec1(size(vec0)) ! result, output vector 
    real(R_) :: cosq1, sinq1

    ! Angles
    if (mhem == 0) then ! both hemispheres
       cosq1 = cosqe * max(reps, 1.0_R_ - sqrt(abs(2.0_R_ * xi - 1.0_R_)))
       cosq1 = cosq1 * sign(1.0_R_, xi - 0.5_R_)
       sinq1 = sqrt(1.0_R_ - cosq1**2)
    else ! forward hemisphere only
       cosq1 = cosqe * max(reps, 1.0_R_ - sqrt(xi))
       sinq1 = sqrt(1.0_R_ - cosq1**2)
    end if

    ! Rotate
    vec1(:) = geo3D_rotateU_1R(vec0, sinq1, cosq1, sinf, cosf)

  end function geo3D_rotateU_ring_1R


  !+
  ! Matrix for Z-Y-Z rotations
  !-
  function geo3D_rotMat_ZYZ_2R(phi, the, psi) result(res) 

    real(R_), intent(in) :: phi ! angle around Z0
    real(R_), intent(in) :: the ! angle around Y1
    real(R_), intent(in) :: psi ! angle around Z2
    real(R_) :: res(3,3)
    real(R_) :: tmp(3,3), s, c, z, u

    z = 0.0_R_
    u = 1.0_R_
    s = sin(phi)
    c = cos(phi)
    res(1:3, 1) = (/ c, s, z/) ! matrix A
    res(1:3, 2) = (/-s, c, z/)
    res(1:3, 3) = (/ z, z, u/)
    s = sin(the)
    c = cos(the)
    tmp(1:3, 1) = (/c, z, -s/) ! matrix B
    tmp(1:3, 2) = (/z, u,  z/)
    tmp(1:3, 3) = (/s, z,  c/)
    res(:,:) = matmul(tmp, res) ! B x A
    s = sin(psi)
    c = cos(psi)
    tmp(1:3, 1) = (/ c, s, z/) ! matrix C
    tmp(1:3, 2) = (/-s, c, z/)
    tmp(1:3, 3) = (/ z, z, u/)
    res(:,:) = matmul(tmp, res) ! C x B x A

  end function geo3D_rotMat_ZYZ_2R


  !+
  ! Scalar product (A.B)
  !-
  function geo3D_scaProd_R(a, b) result(s) 

    real(R_), intent(in)  :: a(:)    ! a 3-D vector A
    real(R_), intent(in)  :: b(:)    ! a 3-D vector B
    real(R_) :: s  ! scalar product, A x B
    s = a(1) * b(1) + a(2) * b(2) + a(3) * b(3)

  end function geo3D_scaProd_R


  !+
  ! Geometry of two unit vectors
  !-
  subroutine geo3D_twoUVec(vec1, vec2, cosq, sinq) 

    real(R_), intent(in)  :: vec1(:) ! a 3-D unit vector A
    real(R_), intent(in)  :: vec2(:) ! a 3-D unit vector B
    real(R_), intent(out) :: cosq, sinq ! cos & sin of angle between A and B
    cosq = geo3D_scaProd_R(vec1, vec2)
    sinq = geo3D_norm_R(geo3D_vecProd_1R(vec1, vec2))

  end subroutine geo3D_twoUVec


  !+
  ! Unit vector, normalized by the norm of vector A
  !-
  function geo3D_UVec_1R(a) result(u)

    real(R_), intent(in)  :: a(:) ! a 3-D vector A
    real(R_) :: u(size(a, 1))     ! the unit vector
    real(R_) :: s

    s = geo3D_norm_R(a)
    if (abs(s) > RSML_) then
       u(1:3) = a(1:3) / s
    else
       u(1:3) = 0.0_R_
    end if

  end function geo3D_UVec_1R


  !+
  ! New direction with isotropic distribution
  !-
  function geo3D_UVec_isotropic_1R(xi, sinf, cosf) result(vec1)

    real(R_), intent(in) :: xi         ! a random number to determine the polar angle
    real(R_), intent(in) :: sinf, cosf ! sin & cos of azimuth angle
    real(R_) :: vec1(3) ! result, output vector 
    real(R_) :: cosq, sinq

    cosq = 2.0_R_ * xi - 1.0_R_
    sinq = 2.0_R_ * sqrt(xi * (1.0_R_ - xi))
    vec1(1) = sinq * cosf
    vec1(2) = sinq * sinf
    vec1(3) = cosq
    call geo3D_UVec_renorm(vec1)

  end function geo3D_UVec_isotropic_1R


  !+
  ! New direction by Lambert law, horizontal surface
  !-
  function geo3D_UVec_LambertH_1R(xi, sinf, cosf) result(vec1)

    real(R_), intent(in) :: xi         ! a random number to determine the polar angle
    real(R_), intent(in) :: sinf, cosf ! sin & cos of azimuth angle
    real(R_) :: vec1(3) ! result, output vector 
    real(R_) :: cosq, sinq

    cosq = sqrt(xi)
    sinq = sqrt(1.0_R_ - xi)
    vec1(1) = sinq * cosf
    vec1(2) = sinq * sinf
    vec1(3) = cosq
    call geo3D_UVec_renorm(vec1)

  end function geo3D_UVec_LambertH_1R


  !+
  ! New direction by Lambert law, general surface with arbitrary orientation
  !-
  function geo3D_UVec_Lambert_1R(eta, xi, sinf, cosf) result(vec1)

    real(R_), intent(in) :: eta(:)     ! a normal vector that specify the hemisphere
    real(R_), intent(in) :: xi         ! a random number to determine the polar angle
    real(R_), intent(in) :: sinf, cosf ! sin & cos of azimuth angle
    real(R_) :: vec1(3) ! result, output vector 
    real(R_) :: cosq, sinq

    cosq = sqrt(xi)
    sinq = sqrt(1.0_R_ - xi)
    vec1(:) = geo3D_rotateU_1R(eta, sinq, cosq, sinf, cosf)

  end function geo3D_UVec_Lambert_1R


  !+
  ! Renormalize a unit vector by using the Newtonian approximation
  !  A given unit vector must satisfy that ux*ux + uy*uy + uz*uz is nearly equal to 1
  !-
  subroutine geo3D_UVec_renorm(uvec)

    real(R_), intent(inout) :: uvec(:) ! unit vector
    real(R_) :: fone
    fone = 0.5_R_ * (3.0_R_ - (uvec(1)**2 + uvec(2)**2 + uvec(3)**2)) ! 1st-order Newton method
    uvec(1:3) = uvec(1:3) * fone

  end subroutine geo3D_UVec_renorm


  !+
  ! Vector product (A x B)
  !-
  function geo3D_vecProd_1R(a, b) result(axb) 

    real(R_), intent(in)  :: a(:), b(:)    ! 3-D vectors A & B
    real(R_) :: axb(size(a, 1))      ! vector product, A x B
    axb(1) = -a(3) * b(2) + a(2) * b(3)
    axb(2) = -a(1) * b(3) + a(3) * b(1)
    axb(3) = -a(2) * b(1) + a(1) * b(2)

  end function geo3D_vecProd_1R

end module hparx_geo3D
