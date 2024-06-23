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
! Library of table-based mathmatical functions
! Explicit initialization is required to use this module: call mtab_init().
! Important note: Functions in this module were made solely for quick and dirty calculations.
!  They may be faster than intrinsic functions, though less accurate. Therefore, they are 
!  useful if high accuracy is not needed. 
!-
module hparx_mtab 

  use globals, only : R_, RD_, PI_, PI2_, PIH_, FRAC16_, RHUGE_
  implicit none
  private

  ! Public
  public :: mtab_init        !* Initialize this module *
  public :: mtab_final       !* Finalize this module *
  public :: mtab_acosX_init  !* Initialize table for acos(x) *
  public :: mtab_acosX_R     !* acos(x) *
  public :: mtab_atanYX_R    !* atan2(y,x) *
  public :: mtab_expX_init   !* Initialize tables for exp(X) *
  public :: mtab_expX_1R     !* Vector of y = exp(x) *
  public :: mtab_expX_R      !* exp(x) *
  public :: mtab_expXC_1R    !* Vector of y = exp(x) - 1 *
  public :: mtab_expXC_R     !* exp(x)-1 *
  public :: mtab_expNX_1R    !* Vector of y = exp(-x) *
  public :: mtab_expNX_R     !* exp(-x) *
  public :: mtab_expNXC_1R   !* Vector of y = 1 - exp(-x) *
  public :: mtab_expNXC_R    !* 1-exp(-x) *
  public :: mtab_expNX_NXC   !* Get both of exp(-x) & 1-exp(-x) *
  public :: mtab_expX_b_1R   !* Better accuracy version, Vector of y = exp(x) *
  public :: mtab_expX_b_R    !* Better accuracy version, exp(x) *
  public :: mtab_expXC_b_1R  !* Better accuracy version, Vector of y = exp(x) - 1 *
  public :: mtab_expXC_b_R   !* Better accuracy version, exp(x)-1 *
  public :: mtab_expNX_b_1R  !* Better accuracy version, Vector of y = exp(-x) *
  public :: mtab_expNX_b_R   !* Better accuracy version, exp(-x) *
  public :: mtab_expNXC_b_1R !* Better accuracy version, Vector of y = 1 - exp(-x) *
  public :: mtab_expNXC_b_R  !* Better accuracy version, 1-exp(-x) *
  public :: mtab_sincos_init   !* Initialize tables for sin(X) & cos(X) *
  public :: mtab_sincos_cycleX !* Get sin(x) & cos(x) using LUTs *
  public :: mtab_sincos        !* Get sin(x) & cos(x) using LUTs, for X in [0,2*pi] *

  ! Tables for exp(x) and exp(-x)
  integer,  save :: Exp_ntab = 100000
  real(R_), save :: Exp_xmin = 0.005_R_
  real(R_), save :: Exp_xmax = 82.0_R_
  real(R_), save :: Exp_xsml = 0.01_R_
  real(R_), save :: Exp_xlrg = 37.0_R_  ! 1 - exp(-xlrg) ~= 1
  real(R_), save :: Exp_xinf = 82.0_R_  ! log(RHUGE_)
  real(R_), save :: Exp_facx
  real(R_), save :: Exp_xbin
  real(R_), save, allocatable :: Exp_etab_p(:), Exp_etab_n(:)

  ! Tables for sin(X) & cos(X)
  integer,  save :: Sctab_ntab = 50000
  real(R_), save :: Sctab_facx
  real(R_), save :: Sctab_xbin
  real(R_), save, allocatable :: Sctab_sin(:), Sctab_cos(:)

  ! Table for acos(x)
  real(R_), parameter :: Acos_ZMAX = 0.99_R_
  integer,  save :: Acos_ntab = 100000
  real(R_), save :: Acos_zbin
  real(R_), save :: Acos_znum
  real(R_), save, allocatable :: Acos_actab(:)

contains

  !+
  ! Initialize this module
  !-
  subroutine mtab_init(ntab_exp, ntab_sc, ntab_ac) 

    integer, intent(in), optional :: ntab_exp, ntab_sc, ntab_ac ! # of table grid points
    !// Note: Numbers of grid points are very important for accuracy.

    call mtab_expX_init(ntab_exp)
    call mtab_sincos_init(ntab_sc)
    call mtab_acosX_init(ntab_ac)

  end subroutine mtab_init


  !+
  ! Finalize this module
  !-
  subroutine mtab_final() 

    if (allocated(Exp_etab_p)) deallocate(Exp_etab_p)
    if (allocated(Exp_etab_n)) deallocate(Exp_etab_n)
    if (allocated(Sctab_sin) ) deallocate(Sctab_sin)
    if (allocated(Sctab_cos) ) deallocate(Sctab_cos)
    if (allocated(Acos_actab)) deallocate(Acos_actab)

  end subroutine mtab_final


  !+
  ! Initialize table for acos(x)
  !-
  subroutine mtab_acosX_init(ntab) 

    integer, intent(in), optional :: ntab ! # of table grid points
    integer :: itab

    ! Asign & allocate
    if (present(ntab)) Acos_ntab = ntab
    Acos_zbin = Acos_ZMAX / Acos_ntab
    Acos_znum = Acos_ntab / Acos_ZMAX
    if (allocated(Acos_actab)) deallocate(Acos_actab)
    allocate (Acos_actab(0:Acos_ntab)) ! the last element is dummy

    ! Make table
    do itab = 0, Acos_ntab - 1
       Acos_actab(itab) = acos(Acos_zbin * (real(itab, R_) + 0.5_R_))
    end do

  end subroutine mtab_acosX_init


  !+
  ! acos(x)
  !-
  function mtab_acosX_R(x) result(res) 

    real(R_), intent(in) :: x ! should be in the range [-1, 1]
    real(R_) :: res, z, y, a, dz

    ! Z in the range from 0 to 1
    if (x >= 0.0_R_) then
       z = x
    else
       z = -x
    end if

    ! Small Z: use LUT
    if (z < Acos_ZMAX) then
       y = Acos_actab(int(z * Acos_znum))

       ! Large Z: 1st-order Newtonian approximation
    else if (z < 1.0_R_) then
       dz = 1.0_R_ - z
       a = 30.0_R_ + dz * dz
       y = sqrt(2.0_R_ * dz) * ((a - 7.5_R_ * dz) / (a - 10.0_R_ * dz))

       ! Too large Z (invalid input!)
    else
       y = 0.0_R_
    end if

    ! Result
    if (x >= 0.0_R_) then
       res = y
    else
       res = PI_ - y
    end if

  end function mtab_acosX_R


  !+
  ! atan2(y,x)
  !// This may be faster than intrinsic atan2(y,x) function by a factor of 2
  !   (assuming x**2 + y**2 = 1).
  !-
  function mtab_atanYX_R(y, x) result(res) 

    real(R_), intent(in) :: y, x ! (Y, X) = (sinA, cosA), shoud be normalized as x**2 + y**2 = 1
    real(R_) :: res ! result, the angle A in [0,pi*2]

    if (abs(x) < 0.8_R_) then ! use x mainly
       res = sign(mtab_acosX_R(x), y)
    else ! use y mainly
       res = sign(mtab_acosX_R(y), -x) + PIH_
       if (res > PI_) res = res - PI2_
    end if

  end function mtab_atanYX_R


  !+
  ! Initialize tables for exp(X) & exp(-X)
  !-
  subroutine mtab_expX_init(ntab, xmin, xmax, xsml, xlrg)

    integer,  intent(in), optional :: ntab ! # of table grid points
    real(R_), intent(in), optional :: xmin, xmax ! min & max of X
    real(R_), intent(in), optional :: xsml, xlrg ! small & large threshold for X
    !// For abs(x) < xsml, 1-exp(+-x) is approximated by using the Taylor expansion
    !   For x > xlrg, 1-exp(-x) = 1
    integer  :: ilut

    ! Asign & allocate
    if (present(ntab)) Exp_ntab = ntab
    if (present(xmin)) Exp_xmin = xmin
    if (present(xmax)) Exp_xmax = xmax
    if (present(xsml)) Exp_xsml = xsml
    if (present(xlrg)) Exp_xlrg = xlrg
    Exp_facx = Exp_ntab / Exp_xmax
    Exp_xbin = Exp_xmax / Exp_ntab
    Exp_xinf = log(RHUGE_)
    if (allocated(Exp_etab_p)) deallocate(Exp_etab_p)
    if (allocated(Exp_etab_n)) deallocate(Exp_etab_n)
    allocate(Exp_etab_p(0:Exp_ntab)) ! the last element is dummy
    allocate(Exp_etab_n(0:Exp_ntab))

    ! Generate a LUT
    do ilut = 0, Exp_ntab
       Exp_etab_p(ilut) = exp( (ilut + 0.5_R_) * Exp_xbin)
       Exp_etab_n(ilut) = exp(-(ilut + 0.5_R_) * Exp_xbin)
    end do

  end subroutine mtab_expX_init


  !+
  ! Vector of y = exp(x)
  !-
  function mtab_expX_1R(x) result(y) 

    real(R_), intent(in) :: x(:)  ! input vector x
    real(R_) :: y(size(x))
    integer  :: ix
    do ix = 1, size(x)
       y(ix) = mtab_expX_R(x(ix))
    end do

  end function mtab_expX_1R


  !+
  ! Vector of y = exp(x) - 1
  !-
  function mtab_expXC_1R(x) result(y) 

    real(R_), intent(in) :: x(:)  ! input vector x
    real(R_) :: y(size(x))
    integer  :: ix
    do ix = 1, size(x)
       y(ix) = mtab_expXC_R(x(ix))
    end do

  end function mtab_expXC_1R


  !+
  ! exp(x)
  !-
  function mtab_expX_R(x) result(res) 

    real(R_), intent(in) :: x  ! argument x for exp(x), in [-Inf, +Inf]
    real(R_) :: res

    ! Small X
    if (x > -Exp_xmin .and. x < Exp_xmin) then
       res = 1.0_R_ + x * (1.0_R_ + 0.5_R_ * x)

       ! Positive X
    else if (x > 0.0_R_) then
       if (x < Exp_xmax) then
          res = Exp_etab_p(int(x * Exp_facx))
       else if (x > Exp_xinf) then
          res = RHUGE_
       else
          res = exp(x)
       end if

       ! Negative X
    else
       if (x > -Exp_xmax) then
          res = Exp_etab_n(int(-x * Exp_facx))
       else
          res = exp(x) ! nearly zero
       end if
    end if

  end function mtab_expX_R


  !+
  ! exp(x)-1
  !-
  function mtab_expXC_R(x) result(res) 

    real(R_), intent(in) :: x  ! argument X, in [-Inf, +Inf]
    real(R_) :: res ! result

    if (x > -Exp_xsml .and. x < Exp_xsml) then ! samall X
       res = x * (1.0_R_ + x * (0.5_R_ + FRAC16_ * x)) ! eps < 1e-4 %
    else if (x < 0.0_R_) then ! negative X
       if (x > -Exp_xlrg) then ! usual case
          res = Exp_etab_n(int(-x * Exp_facx)) - 1.0_R_
       else ! large negative
          res = -1.0_R_
       end if
    else ! positive X
       if (x < Exp_xmax) then
          res = Exp_etab_p(int(x * Exp_facx)) - 1.0_R_
       else if (x > Exp_xinf) then
          res = RHUGE_
       else
          res = exp(x) - 1.0_R_
       end if
    end if

  end function mtab_expXC_R


  !+
  ! Vector of y = exp(-x)
  !-
  function mtab_expNX_1R(x) result(y) 

    real(R_), intent(in) :: x(:)  ! input vector x
    real(R_) :: y(size(x))
    integer  :: ix
    do ix = 1, size(x)
       y(ix) = mtab_expNX_R(x(ix))
    end do

  end function mtab_expNX_1R


  !+
  ! Vector of y = 1 - exp(-x)
  !-
  function mtab_expNXC_1R(x) result(y) 

    real(R_), intent(in) :: x(:)  ! input vector x
    real(R_) :: y(size(x))
    integer  :: ix
    do ix = 1, size(x)
       y(ix) = mtab_expNXC_R(x(ix))
    end do

  end function mtab_expNXC_1R


  !+
  ! exp(-x)
  !-
  function mtab_expNX_R(x) result(res) 

    real(R_), intent(in) :: x ! X, positive, in [0, +Inf]
    real(R_) :: res

    if (x < Exp_xmin) then ! small X
       res = 1.0_R_ - x * (1.0_R_ - 0.5_R_ * x)
    else if (x < Exp_xmax) then ! usual case
       res = Exp_etab_n(int(x * Exp_facx))
    else ! large X
       res = exp(-x)
    end if

  end function mtab_expNX_R


  !+
  ! 1-exp(-x)
  !-
  function mtab_expNXC_R(x) result(res) 

    real(R_), intent(in) :: x  ! X, positive, in [0, +Inf]
    real(R_) :: res ! result

    if (x < Exp_xsml) then ! samall X
       res = x * (1.0_R_ - x * (0.5_R_ - FRAC16_ * x)) ! eps < 1e-4 %
    else if (x < Exp_xlrg) then ! usual case
       res = 1.0_R_ - Exp_etab_n(int(x * Exp_facx))
    else ! large X
       res = 1.0_R_
    end if

  end function mtab_expNXC_R


  !+
  ! Get both of exp(-x) & 1-exp(-x)
  !-
  subroutine mtab_expNX_NXC(x, e, ec) 

    real(R_), intent(in) :: x ! X, positive, in [0, +Inf]
    real(R_), intent(out) :: e, ec ! exp(-x), 1-exp(-x)

    if (x < Exp_xsml) then ! samall X
       ec = x * (1.0_R_ - x * (0.5_R_ - FRAC16_ * x)) ! eps < 1e-4 %
       e = 1.0_R_ - ec
    else if (x < Exp_xmax) then ! usual case
       e = Exp_etab_n(int(x * Exp_facx))
       ec = 1.0_R_ - e
    else ! large X
       e = exp(-x)
       ec = 1.0_R_
    end if

  end subroutine mtab_expNX_NXC


  !+
  ! Vector of y = exp(x), better accuracy
  !-
  function mtab_expX_b_1R(x) result(y) 

    real(R_), intent(in) :: x(:)  ! input vector x
    real(R_) :: y(size(x))
    integer  :: ix
    do ix = 1, size(x)
       y(ix) = mtab_expX_b_R(x(ix))
    end do

  end function mtab_expX_b_1R


  !+
  ! Vector of y = exp(x) - 1, better accuracy
  !-
  function mtab_expXC_b_1R(x) result(y) 

    real(R_), intent(in) :: x(:)  ! input vector x
    real(R_) :: y(size(x))
    integer  :: ix
    do ix = 1, size(x)
       y(ix) = mtab_expXC_b_R(x(ix))
    end do

  end function mtab_expXC_b_1R


  !+
  ! exp(x), better accuracy
  !-
  function mtab_expX_b_R(x) result(res) 

    real(R_), intent(in) :: x  ! argument x for exp(x), in [-Inf, +Inf]
    real(R_) :: res
    integer :: itab

    ! Small X
    if (x > -Exp_xmin .and. x < Exp_xmin) then
       res = 1.0_R_ + x * (1.0_R_ + 0.5_R_ * x)

       ! Positive X
    else if (x > 0.0_R_) then
       if (x < Exp_xmax) then
          itab = int(x * Exp_facx)
          res = Exp_etab_p(itab) * (1.0_R_ + (x - (itab + 0.5_R_) * Exp_xbin))
       else if (x > Exp_xinf) then
          res = RHUGE_
       else
          res = exp(x)
       end if

       ! Negative X
    else
       if (x > -Exp_xmax) then
          itab = int(-x * Exp_facx)
          res = Exp_etab_n(itab) * (1.0_R_ + (x + (itab + 0.5_R_) * Exp_xbin))
       else
          res = exp(x) ! nearly zero
       end if
    end if

  end function mtab_expX_b_R


  !+
  ! exp(x)-1, better accuracy
  !-
  function mtab_expXC_b_R(x) result(res) 

    real(R_), intent(in) :: x  ! argument X, in [-Inf, +Inf]
    real(R_) :: res ! result
    integer :: itab

    if (x > -Exp_xsml .and. x < Exp_xsml) then ! samall X
       res = x * (1.0_R_ + x * (0.5_R_ + FRAC16_ * x)) ! eps < 1e-4_R_%
    else if (x < 0.0_R_) then ! negative X
       if (x > -Exp_xlrg) then ! usual case
          itab = int(-x * Exp_facx)
          res = Exp_etab_n(itab) * (1.0_R_ + (x + (itab + 0.5_R_) * Exp_xbin)) - 1.0_R_
       else ! large negative
          res = -1.0_R_
       end if
    else ! positive X
       if (x < Exp_xmax) then
          itab = int(x * Exp_facx)
          res = Exp_etab_p(itab) * (1.0_R_ + (x - (itab + 0.5_R_) * Exp_xbin)) - 1.0_R_
       else if (x > Exp_xinf) then
          res = RHUGE_
       else
          res = exp(x) - 1.0_R_
       end if
    end if

  end function mtab_expXC_b_R


  !+
  ! Vector of y = exp(-x), better accuracy
  !-
  function mtab_expNX_b_1R(x) result(y) 

    real(R_), intent(in) :: x(:)  ! input vector x
    real(R_) :: y(size(x))
    integer  :: ix
    do ix = 1, size(x)
       y(ix) = mtab_expNX_b_R(x(ix))
    end do

  end function mtab_expNX_b_1R


  !+
  ! Vector of y = 1 - exp(-x), better accuracy
  !-
  function mtab_expNXC_b_1R(x) result(y) 

    real(R_), intent(in) :: x(:)  ! input vector x
    real(R_) :: y(size(x))
    integer  :: ix
    do ix = 1, size(x)
       y(ix) = mtab_expNXC_b_R(x(ix))
    end do

  end function mtab_expNXC_b_1R


  !+
  ! exp(-x), better accuracy
  !-
  function mtab_expNX_b_R(x) result(res) 

    real(R_), intent(in) :: x ! X, positive, in [0, +Inf]
    real(R_) :: res
    integer :: itab

    if (x < Exp_xmin) then ! small X
       res = 1.0_R_ - x * (1.0_R_ - 0.5_R_ * x)
    else if (x < Exp_xmax) then ! usual case
       itab = int(x * Exp_facx)
       res = Exp_etab_n(itab) * (1.0_R_ - (x - (itab + 0.5_R_) * Exp_xbin))
    else ! large X
       res = exp(-x)
    end if

  end function mtab_expNX_b_R


  !+
  ! 1-exp(-x), better accuracy
  !-
  function mtab_expNXC_b_R(x) result(res) 

    real(R_), intent(in) :: x  ! X, positive, in [0, +Inf]
    real(R_) :: res ! result
    integer :: itab

    if (x < Exp_xsml) then ! samall X
       res = x * (1.0_R_ - x * (0.5_R_ - FRAC16_ * x)) ! eps < 1e-4 %
    else if (x < Exp_xlrg) then ! usual case
       itab = int(x * Exp_facx)
       res = 1.0_R_ - Exp_etab_n(itab) * (1.0_R_ - (x - (itab + 0.5_R_) * Exp_xbin))
    else ! large X
       res = 1.0_R_
    end if

  end function mtab_expNXC_b_R


  !+
  ! Initialize tables for sin(X) & cos(X)
  !-
  subroutine mtab_sincos_init(ntab) 

    integer, intent(in), optional :: ntab ! # of table grid points
    integer  :: ilut

    ! Asign & allocate
    if (present(ntab)) Sctab_ntab = ntab
    Sctab_facx = Sctab_ntab / PI2_
    Sctab_xbin = PI2_ / Sctab_ntab
    if (allocated(Sctab_sin)) deallocate(Sctab_sin)
    if (allocated(Sctab_cos)) deallocate(Sctab_cos)
    allocate (Sctab_sin(0:Sctab_ntab+1))
    allocate (Sctab_cos(0:Sctab_ntab+1))

    ! Make tables
    do ilut = 0, Sctab_ntab
       Sctab_sin(ilut) = sin(ilut * Sctab_xbin)
       Sctab_cos(ilut) = cos(ilut * Sctab_xbin)
    end do

  end subroutine mtab_sincos_init


  !+
  ! Get sin(x) & cos(x) using LUTs
  !-
  subroutine mtab_sincos_cycleX(x, sinx, cosx) 

    real(R_), intent(in) :: x  ! argument x
    real(R_), intent(out) :: sinx, cosx ! sin(x), cos(x)
    real(R_) :: x1

    if (x >= 0.0_R_ .and. x <= PI2_) then
       call mtab_sincos(x, sinx, cosx)
    else
       x1 = mod(x, PI2_)
       if (x1 < 0.0_R_) x1 = x1 + PI2_
       call mtab_sincos(x1, sinx, cosx)
    end if

  end subroutine mtab_sincos_cycleX


  !+
  ! Get sin(x) & cos(x) using LUTs, for X in [0,2*pi]
  !-
  subroutine mtab_sincos(x, sinx, cosx) 

    real(R_), intent(in) :: x  ! argument x in [0, 2*pi]
    !// Caution! If x is < 0 or > 2*pi, then SIGSEG can appear!
    real(R_), intent(out) :: sinx, cosx ! sin(x), cos(x)
    integer  :: ilut
    real(R_) :: dx

    ilut = int(x * Sctab_facx + 0.5_R_)
    dx = x - ilut * Sctab_xbin
    sinx = Sctab_sin(ilut) + dx * Sctab_cos(ilut)
    cosx = Sctab_cos(ilut) - dx * Sctab_sin(ilut)

  end subroutine mtab_sincos

end module hparx_mtab
