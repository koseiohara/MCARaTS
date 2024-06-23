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
! Library of statistics procedures
!-
module hparx_stat

  use globals
  use hparx_math
  use hparx_lina
  implicit none
  private

  ! Public
  public :: prof_Gamma_1R
  public :: stat_basic
  public :: stat_conf_95p_R
  public :: stat_coVec
  public :: stat_covar
  public :: stat_Gamma_frac_R
  public :: stat_Gamma_kernels_2R
  public :: stat_Gamma_params
  public :: stat_Gamma_stat
  public :: stat_Gauss_pdf_R
  public :: stat_Gauss_R       ! Gaussian distribution function
  public :: stat_hist0_1R
  public :: stat_lognorm_kernels_2R
  public :: stat_lognorm_params
  public :: stat_lognorm_pdf_R
  public :: stat_lognorm_R     ! Lognormal Distribution Function
  public :: stat_lognorm_stat
  public :: stat_meanSdev
  public :: stat_moments_1R
  public :: stat_multiCorCov
  public :: stat_probDens_1R
  public :: stat_prof_Gamma_1R
  public :: stat_histStat
  public :: stat_sReg
  public :: stat_Tdist_5p_R

contains


  !+
  ! Returns a normalized profile weights synthesized by multiple Gamma distributions
  !-
  function prof_Gamma_1R(frac, xmin, xmax, ialp, xave, xgrd, nx) result(prof)

    real(R_), intent(in)  :: frac(:) ! (nk) fractions of each kernel
    real(R_), intent(in)  :: xmin(:) ! (nk) minimum of X
    real(R_), intent(in)  :: xmax(:) ! (nk) maximum of X
    integer,  intent(in)  :: ialp(:) ! (nk) alpha (scaling parameter) (should be < ~10)
    real(R_), intent(in)  :: xave(:) ! (nk) X average of Gamma distribution
    real(R_), intent(in)  :: xgrd(:) ! (nx+1) X grid values
    integer,  intent(in)  :: nx      ! # of X intervals
    real(R_) :: prof(nx) ! normalized profile (weighting function)
    integer  :: ix, iker
    real(R_) :: b, xb, xt, fac

    prof(:) = 0.0_R_
    do iker = 1, size(ialp)
       b = xave(iker) / real(ialp(iker), R_)
       fac = frac(iker) / (sum(frac) * stat_Gamma_frac_R(xmin(iker), xmax(iker), ialp(iker), b))
       do ix = 1, nx
          xb = max(xgrd(ix), xmin(iker))
          xt = min(xgrd(ix + 1), xmax(iker))
          if (xb < xt) prof(ix) = prof(ix) + fac * stat_Gamma_frac_R(xb, xt, ialp(iker), b)
       end do
    end do
    b = 1.0_R_ / sum(prof)
    prof(:) = prof(:) * b

  end function prof_Gamma_1R


  !+
  ! Basic statistics
  !-
  subroutine stat_basic(xvec, xmin, xmax, xmn, xsd, xsk, xkr, mask, xfrc)

    real(R_), intent(in)  :: xvec(:)    ! vector data
    real(R_), intent(out) :: xmin, xmax ! min & max
    real(R_), intent(out) :: xmn, xsd   ! mean & standard deviation
    real(R_), intent(out) :: xsk, xkr   ! skewness & kurtosis 
    logical,  intent(in),  optional :: mask(:) ! mask for xvec
    real(R_), intent(out), optional :: xfrc    ! fraction of active data
    real(RD_) :: sum2, sum3, sum4
    real(R_)  :: d, d2
    integer   :: nx, ix, na

    ! Initialize
    nx = size(xvec)
    sum2 = 0.0_RD_
    sum3 = 0.0_RD_
    sum4 = 0.0_RD_

    ! If mask is available
    if (present(mask) .and. present(xfrc)) then
       na = count(mask)
       if (na == 0) then
          xmin = 0.0_R_
          xmax = 0.0_R_
       else
          xmin = minval(xvec, mask)
          xmax = maxval(xvec, mask)
       end if
       xmn  = sum(xvec, mask) / max(1, na)
       do ix = 1, nx
          if (mask(ix)) then
             d = xvec(ix) - xmn
             d2 = d * d
             sum2 = sum2 + d2
             sum3 = sum3 + d2 * d
             sum4 = sum4 + d2 * d2
          end if
       end do
       xfrc = real(na, R_) / real(nx, R_)

       ! Else
    else
       na = nx
       xmin = minval(xvec)
       xmax = maxval(xvec)
       xmn  = sum(xvec) / na
       do ix = 1, nx
          d = xvec(ix) - xmn
          d2 = d * d
          sum2 = sum2 + d2
          sum3 = sum3 + d2 * d
          sum4 = sum4 + d2 * d2
       end do
    end if

    ! Higher moments
    sum2 = sum2 / max(1, na)
    sum3 = sum3 / max(1, na)
    sum4 = sum4 / max(1, na)
    xsd  = sqrt(sum2)
    if (sum2 > 0.0_RD_) then
       xsk = sum3 / (sum2 * xsd)
       xkr = sum4 / (sum2 * sum2)
    else
       xsk = 0.0_R_
       xkr = 0.0_R_
    endif

  end subroutine stat_basic


  !+
  ! Confidency width of mean at 95% significant level
  !-
  function stat_conf_95p_R(drms, n) result(res)

    real(R_), intent(in) :: drms ! root-mean-square
    integer,  intent(in) :: n    ! # of samples
    real(R_)  :: res
    res = stat_Tdist_5p_R(n - 1) * drms / sqrt(real(n, R_))

  end function stat_conf_95p_R


  !+
  ! Covariance of two vectors
  !-
  subroutine stat_covar(vec1, vec2, covar, ave1, ave2)

    real(R_), intent(in)  :: vec1(:), vec2(:) ! vectars
    real(R_), intent(out) :: covar ! covariance
    real(R_), optional, intent(out) :: ave1, ave2 ! mean values
    real(R_)  :: a1, a2
    integer   :: n

    n = size(vec1)
    a1 = sum(vec1) / n
    a2 = sum(vec2) / n
    covar = sum(vec1 * vec2) / n - a1 * a2
    if (present(ave1)) ave1 = a1
    if (present(ave2)) ave2 = a2

  end subroutine stat_covar


  !+
  ! Calc. mean, stdev, covariance, correlation coef. of given two vectors
  !-
  subroutine stat_coVec(dat1, dat2, cor, cov, dave1, dave2, dsig1, dsig2)

    real(R_), intent(in)  :: dat1(:), dat2(:) ! data vectors
    real(R_), intent(out) :: cor, cov     ! correlation & covariance
    real(R_), intent(out) :: dave1, dave2 ! averages
    real(R_), intent(out) :: dsig1, dsig2 ! standard deviation
    integer   :: ndat, idat
    real(RD_) :: sum11, sum22, sum12
    real(R_)  :: dev1, dev2

    ! Too few data
    ndat = size(dat1)
    if (ndat <= 1) then
       dave1 = 0.0_R_
       dave2 = 0.0_R_
       dsig1 = 0.0_R_
       dsig2 = 0.0_R_
       cov = 0.0_R_
       cor = 0.0_R_
       return
    endif

    ! Analyze
    dave1 = sum(dat1) / ndat
    dave2 = sum(dat2) / ndat
    sum11 = 0.0_RD_
    sum22 = 0.0_RD_
    sum12 = 0.0_RD_
    do idat = 1, ndat
       dev1 = dat1(idat) - dave1
       dev2 = dat2(idat) - dave2
       sum11 = sum11 + dev1**2
       sum22 = sum22 + dev2**2
       sum12 = sum12 + dev1 * dev2
    end do
    dsig1 = sqrt(sum11 / ndat)
    dsig2 = sqrt(sum22 / ndat)
    cov = sum12 / ndat
    cor = 0.0_R_
    if (dsig1 * dsig2 > 0.0_R_) cor = cov / (dsig1 * dsig2)

  end subroutine stat_coVec


  !+
  ! Fraction of a range of Gamma distribution with an integer scale parameter (alpha)
  !-
  function stat_Gamma_frac_R(xmin, xmax, ialp, beta) result(f)

    real(R_), intent(in) :: xmin
    real(R_), intent(in) :: xmax
    integer,  intent(in) :: ialp
    real(R_), intent(in) :: beta
    real(R_) :: f
    real(RD_) :: smin, smax, zmin, zmax
    integer  :: i, ii

    zmin = real(xmin / beta, RD_)
    zmax = real(xmax / beta, RD_)
    smin = 1.0_RD_
    smax = 1.0_RD_
    ii = 1
    do i = 1, ialp - 1
       ii = ii * i
       smin = smin + zmin**i / ii
       smax = smax + zmax**i / ii
    end do
    f = exp(-zmin) * smin - exp(-zmax) * smax

  end function stat_Gamma_frac_R


  !+
  ! Synthesize integration kernels for the Gamma distribution
  !-
  function stat_Gamma_kernels_2R(alp, bet, xvec, ndim) result(fknl)

    real(R_), intent(in) :: alp     ! alpha, the shape parameter (> 0)
    real(R_), intent(in) :: bet     ! beta, the scaling parameter (> 0)
    real(R_), intent(in) :: xvec(:) ! x vector
    integer,  intent(in) :: ndim    ! the max order (>= 0) for kernels
    real(R_) :: fknl(size(xvec)-1, 0:ndim) ! normalized kernels
    !// fknl(ibin, idim) = integral(x=xvec(ibin),xvec(ibin+1)){PDF(x) * x**idim}
    real(R_) :: xlo, xhi, x, w
    integer :: ibin, idim, nbin, ig
    !real(R_), parameter :: ALPMIN = 0.001_R_
    integer,  parameter :: NG = 20
    real(R_), save :: xg(NG), wg(NG)
    integer,  save :: init = 0

    ! Initialize
    nbin = size(xvec) - 1
    if (init == 0) then
       call gaussLegen(NG, xg, wg)
       init = 1
    end if

    ! Loop for size bins
    do ibin = 1, nbin
       xlo = xvec(ibin)     / bet
       xhi = xvec(ibin + 1) / bet

       ! Dummy bin
       if (xlo >= xhi) then
          fknl(ibin, :) = 0.0_R_

          ! Positive alpha
          !else if (alp > ALPMIN) then ! use analytical formulas
          !   w = 1.0_R_
          !   fknl(ibin, 0) = gammp_s(alp, xhi) - gammp_s(alp, xlo)
          !   do idim = 1, ndim
          !      w = w * (alp + idim - 1) * bet
          !      fknl(ibin, idim) = w * (gammp_s(alp + idim, xhi) - gammp_s(alp + idim, xlo))
          !   end do

          ! Nearly-zero or negative
       else ! use numerical integration (Gaussian quadrature)
          fknl(ibin, :) = 0.0_R_
          do ig = 1, NG
             x = xlo + (xhi - xlo) * (xg(ig) + 1.0_R_) * 0.5_R_
             w = wg(ig) * (xhi - xlo) * (x**(alp - 1.0_R_) * exp(-x))
             fknl(ibin, 0) = fknl(ibin, 0) + w
             do idim = 1, ndim
                w = w * x * bet
                fknl(ibin, idim) = fknl(ibin, idim) + w
             end do
          end do
       end if
    end do

    ! Normalize
    w = sum(fknl(1:nbin, 0))
    fknl(:,:) = fknl(:,:) / w

  end function stat_Gamma_kernels_2R


  !+
  ! Moments of Gamma Distribution (GD)
  !  mean, stdev, variance, skewness, kurtsis
  !-
  subroutine stat_Gamma_stat(alpha, beta, xave, xsig, xvar, xskw, xkrt)

    real(R_), intent(in) :: alpha
    real(R_), intent(in) :: beta
    real(R_), intent(out) :: xave, xsig, xvar
    real(R_), intent(out) :: xskw, xkrt

    xave = alpha * beta
    xvar = xave * beta
    xsig = sqrt(xvar)
    xskw = 2.0_R_ / sqrt(alpha)
    xkrt = 3.0_R_ + 6.0_R_ / alpha

  end subroutine stat_Gamma_stat


  !+
  ! Parameter estimation of Gamma Distribution (GD) by Moment method
  !-
  subroutine stat_Gamma_params(xave, xvar, alpha, beta)

    real(R_), intent(in)  :: xave
    real(R_), intent(in)  :: xvar
    real(R_), intent(out) :: alpha
    real(R_), intent(out) :: beta
    alpha = xave**2 / xvar
    beta = xvar / xave

  end subroutine stat_Gamma_params


  !+
  ! Gaussian distribution PDF
  !-
  function stat_Gauss_pdf_R(xmyu, xsgm, x) result(res)

    real(R_)  :: xmyu ! average
    real(R_)  :: xsgm ! sigma
    real(R_)  :: x    ! X
    real(R_)  :: res, z

    if (xsgm > 0.0_R_) then
       z = (x - xmyu) / xsgm
       res = gauss_R(z)
    else
       res = 0.0_R_
    endif

  end function stat_Gauss_pdf_R


  !+
  ! Gaussian distribution function
  !-
  function stat_Gauss_R(xmyu, xsgm, x) result(res) 

    real(R_), intent(in)  :: xmyu ! average
    real(R_), intent(in)  :: xsgm ! sigma
    real(R_), intent(in)  :: x    ! X
    real(R_)  :: res, z

    if (xsgm > 0.0_R_) then
       z = (x - xmyu) / xsgm
       if (z >= 0.0_R_) then
          res = 0.5_R_ * (1.0_R_ + erf_R(z))
       else
          res = 0.5_R_ * (1.0_R_ - erf_R(-z))
       endif
    else
       res = 0.0_R_
    endif

  end function stat_Gauss_R


  !+
  ! Calculate histgram from vector data
  !-
  function stat_hist0_1R(xdat, xmin, xmax, nbin, wgt) result (pdf)

    real(R_), intent(in) :: xdat(:)    ! X data
    real(R_), intent(in) :: xmin, xmax ! X min & max for the histgram
    integer,  intent(in) :: nbin       ! # of bins
    real(R_), intent(in), optional :: wgt(:) ! weights for the data
    real(R_) :: pdf(0:nbin+1) ! histgram (pdf(0) for xdat < xmin, pdf(nbin+1) for xdat >= xmax)
    integer  :: ndat, idat, ibin, nbin1
    real(R_) :: f, slp

    ! Initialize
    ndat = size(xdat)
    nbin1 = nbin + 1
    slp = nbin / (xmax - xmin)

    ! Analyze
    pdf(:) = 0.0_R_
    do idat = 1, ndat
       ibin = int(slp * (xdat(idat) - xmin) + 1.0_R_)
       ibin = max(0, min(nbin + 1, ibin))
       if (present(wgt)) then
          pdf(ibin) = pdf(ibin) + wgt(idat)
       else
          pdf(ibin) = pdf(ibin) + 1.0_R_
       end if
    end do

    ! Results
    if (present(wgt)) then
       f = 1.0_R_ / sum(wgt)
    else
       f = 1.0_R_ / real(ndat, R_)
    end if
    pdf(:) = pdf(:) * f ! normalized PDF

  end function stat_hist0_1R


  !+
  ! Histogram to basic statistics
  !-
  subroutine stat_histStat(x, pdf, sump, av, sd)

    real(R_), intent(in) :: x(:), pdf(:) ! bin-center X & frequency at each X bin
    real(R_), intent(out) :: sump   ! total of the pdf
    real(R_), intent(out) :: av, sd ! mean & sigma of X
    real(RD_) :: sum0, sum1, sum2
    integer :: ix

    sum0 = 0.0_RD_
    sum1 = 0.0_RD_
    sum2 = 0.0_RD_
    do ix = 1, size(pdf)
       sum0 = sum0 + pdf(ix)
       sum1 = sum1 + pdf(ix) * x(ix)
       sum2 = sum2 + pdf(ix) * x(ix)**2
    end do
    sump = sum0
    if (sum0 > RDSML_) then
       av = sum1 / sum0
       sd = sqrt(max(0.0_RD_, sum2 / sum0 - av**2))
    else
       av = 0.0_R_
       sd = 0.0_R_
    end if

  end subroutine stat_histStat


  !+
  ! Synthesize integration kernels for the lognormal distribution
  !-
  function stat_lognorm_kernels_2R(xmod, sig, xvec, ndim) result(fknl)

    real(R_), intent(in) :: xmod    ! mode (> 0)
    real(R_), intent(in) :: sig     ! standard deviation (sigma = ln(s))
    real(R_), intent(in) :: xvec(:) ! x vector
    integer,  intent(in) :: ndim    ! the max order (>= 0) for kernels
    real(R_) :: fknl(size(xvec)-1, 0:ndim) ! normalized kernels
    !// fknl(ibin, idim) = integral(x=xvec(ibin),xvec(ibin+1)){PDF(x) * x**idim}
    real(R_) :: xlo, xhi, x, w
    integer :: ibin, idim, nbin, ig
    integer,  parameter :: NG = 10
    real(R_), save :: xg(NG), wg(NG)
    integer,  save :: init = 0

    ! Initialize
    nbin = size(xvec) - 1
    if (init == 0) then
       call gaussLegen(NG, xg, wg)
       init = 1
    end if

    ! Loop for size bins
    do ibin = 1, nbin
       xlo = xvec(ibin)     / xmod
       xhi = xvec(ibin + 1) / xmod

       ! Dummy bin
       if (xlo >= xhi) then
          fknl(ibin, :) = 0.0_R_

          ! Numerical integration (Gaussian quadrature)
       else
          fknl(ibin, :) = 0.0_R_
          do ig = 1, NG
             x = xlo + (xhi - xlo) * (xg(ig) + 1.0_R_) * 0.5_R_
             w = wg(ig) * (xhi - xlo) * exp(-0.5 * (log(x) / sig)**2) / (x * xmod * sig)
             fknl(ibin, 0) = fknl(ibin, 0) + w
             do idim = 1, ndim
                w = w * x * xmod
                fknl(ibin, idim) = fknl(ibin, idim) + w
             end do
          end do
       end if
    end do

    ! Normalize
    w = sum(fknl(1:nbin, 0))
    fknl(:,:) = fknl(:,:) / w

  end function stat_lognorm_kernels_2R


  !+
  ! Lognormal Distribution Function
  !-
  function stat_lognorm_R(xmyu, xsgm, x) result(res) 

    real(R_), intent(in) :: xmyu ! average of log(x)
    real(R_), intent(in) :: xsgm ! standard deviation of log(x)
    real(R_), intent(in) :: x
    real(R_)  :: res, z

    if (xsgm > 0.0_R_ .and. x > 0.0_R_) then
       z = (log(x) - xmyu) / xsgm
       if (z >= 0.0_R_) then
          res = 0.5_R_ * (1.0_R_ + erf_R(z))
       else
          res = 0.5_R_ * (1.0_R_ - erf_R(-z))
       endif
    else
       res = 0.0_R_
    endif

  end function stat_lognorm_R


  !+
  ! Moments of Lognormal Distribution (LND)
  !  mean, stdev, variance, skewness, kurtsis.
  !-
  subroutine stat_lognorm_stat(xmyu, xsgm, xave, xsig, xvar, xskw, xkrt)

    real(R_), intent(in)  :: xmyu, xsgm ! average & standard deviation of log(x)
    real(R_), intent(out) :: xave, xsig, xvar ! average, standard deviation, and variance of x
    real(R_), intent(out) :: xskw, xkrt ! skewness and kurtosis of x
    real(R_)  :: expxmyu, expxsgm2, xsgm2

    xsgm2 = xsgm * xsgm
    expxmyu  = exp(xmyu)
    expxsgm2 = exp(xsgm2)
    xave = expxmyu * sqrt(expxsgm2)
    xvar = expxmyu * expxmyu * expxsgm2 * (expxsgm2 - 1.0_R_)
    xsig = sqrt(max(0.0_R_, xvar))
    xskw = (expxsgm2 + 2.0_R_) * sqrt(max(0.0_R_, expxsgm2 - 1.0_R_))
    xkrt = expxsgm2 * expxsgm2 * (3.0_R_ + expxsgm2 * (2.0_R_ + expxsgm2)) - 3.0_R_

  end subroutine stat_lognorm_stat


  !+
  ! Lognormal Distribution Probability Density Function, dP/dlnX
  !-
  function stat_lognorm_pdf_R(xmyu, xsgm, x) result(res)

    real(R_), intent(in) :: xmyu, xsgm ! average and standard deviation of log(x)
    real(R_), intent(in) :: x ! x data at the target point
    real(R_)  :: res, z

    if (xsgm > 0.0_R_ .and. x > 0.0_R_) then
       z = (log(x) - xmyu) / xsgm
       res = stat_Gauss_pdf_R(0.0_R_, 1.0_R_, z)
    else
       res = 0.0_R_
    endif

  end function stat_lognorm_pdf_R


  !+
  ! Parameter estimation of Lognormal Distribution (LND) by Moment method
  !-
  subroutine stat_lognorm_params(xmean, xvar, xmyu, xsgm, xsgm2)

    real(R_), intent(in)  :: xmean, xvar ! average & variance of x
    real(R_), intent(out) :: xmyu, xsgm, xsgm2 ! average, standard deviation, and variance of log(x)
    real(R_)  :: xmean2

    xmean2 = xmean * xmean
    xmyu  = log(xmean2 / sqrt(xmean2 + xvar))
    xsgm2 = log((xmean2 + xvar) / xmean2)
    xsgm  = sqrt(max(0.0_R_, xsgm2))

  end subroutine stat_lognorm_params


  !+
  ! Mean and standard deviation of a single vector data
  !-
  subroutine stat_meanSdev(vec, ave, dev)

    real(R_), intent(in)  :: vec(:)   ! vector data
    real(R_), intent(out) :: ave, dev ! mean & standard deviation
    integer   :: n

    n = size(vec)
    ave = sum(vec) / n
    dev = sqrt(sum((vec - ave)**2) / n)

  end subroutine stat_meanSdev


  !+
  ! Calculate the moments of arbitrary orders
  !-
  function stat_moments_1R(dat, nm) result(res) 

    real(R_), intent(in)  :: dat(:) ! data vector
    integer,  intent(in)  :: nm     ! max order of moments
    real(R_) :: res(nm) ! moments of orders 1-nm
    real(R_) :: wrk(size(dat))
    integer  :: n, im

    n = size(dat)
    wrk(:) = dat(:)
    res(1) = sum(wrk(:)) / n
    do im = 2, nm
       wrk(:) = wrk(:) * dat(:)
       res(im) = sum(wrk(:)) / n
    end do

    !real(R_)  :: d, dd(nm), wrk(size(dat))
    !real(RD_) :: wsum(nm)
    !wsum(:) = 0.0_RD_
    !do id = 1, size(dat)
    !   d = dat(id)
    !   dd(1) = d
    !   do im = 2, nm
    !      dd(im) = dd(im - 1) * d
    !   end do
    !   wsum(:) = wsum(:) + dd(:)
    !end do
    !mom(:) = wsum(:) / size(dat)

  end function stat_moments_1R


  !+
  ! Correlation matrix & covariance matrix for multiple variables
  !-
  subroutine stat_multiCorCov(dat, ave, var, cor, cov)

    real(R_), intent(in)  :: dat(:, :) !(nv,ndat) data
    real(R_), intent(out) :: ave(:), var(:) !(nv) average & variance
    real(R_), intent(out) :: cor(:, :) !(nv,nv) correlation matrix
    real(R_), intent(out) :: cov(:, :) !(nv,nv) covariance matrix
    integer   :: nv, ndat, i, idat, iv, j
    real(RD_) :: sum2

    ! Average
    nv   = size(dat, 1)
    ndat = size(dat, 2)
    do iv = 1, nv
       ave(iv) = sum(dat(iv, :)) / real(ndat, R_)
    end do

    ! Variance & covariance matrix
    do i = 1, nv
       do j = i, nv
          sum2 = 0.0_RD_
          do idat = 1, ndat
             sum2 = sum2 + (dat(i, idat) - ave(i)) * (dat(j, idat) - ave(j))
          end do
          cov(i, j) = sum2 / real(ndat - 1, R_)
          cov(j, i) = cov(i, j)
       end do
       var(i) = cov(i, i)
    end do

    ! Correlation matrix
    do i = 1, nv
       cor(i, i:nv) = cov(i, i:nv) / sqrt(var(i) * var(i:nv))
       cor(i:nv, i) = cor(i, i:nv)
    end do

  end subroutine stat_multiCorCov


  !+
  ! Probability density function, p(x)=dP/dx normalized
  !-
  function stat_probDens_1R(md, pd, x) result(res)

    integer,  intent(in) :: md    ! flag for a kind of distribution
    real(R_), intent(in) :: pd(:) ! parameter packet (see below)
    real(R_), intent(in) :: x(:)  ! X data vector
    ! -   md=1: Power-law
    !             pd(1)=X0, pd(2)=X1, pd(3)=X2, pd(4)=a1,
    !             pd(5)=a2, pd(6)=b1, pd(7)=b2, pd(8)=f
    !     md=2: Gamma
    !             pd(1)=beta, pd(2)=alpha
    !     md=3: Modified Gamma
    !             pd(1)=beta, pd(2)=alpha, pd(3)=gamma
    !     md=4: Lognormal dX/dr
    !             pd(1)=Xm, pd(2)=s>1
    !     md=5: Lognormal volume dY/dr, Y=X**3
    !             pd(1)=Xm, pd(2)=s>1
    real(R_)  :: res(size(x)), a, b
    integer  :: ix

    ! Power-law
    if (md == 1) then
       do ix = 1, size(x)
          if (x(ix) < pd(1) .or. x(ix) > pd(3)) then
             res(ix) = 0.0_R_
          else if (x(ix) < pd(2)) then
             res(ix) = pd(8) * pd(4) * x(ix)**pd(6)
          else
             res(ix) = pd(8) * pd(5) * x(ix)**pd(7)
          end if
       end do

       ! Gamma
    else if (md == 2) then
       b = gammaLog_R(pd(2))
       res(:) = x(:) / pd(1)
       res(:) = exp(pd(2) * log(res(:)) - b - res(:)) / x(:)

       ! Modified Gamma
    else if (md == 3) then
       b = gammaLog_R(pd(2))
       res(:) = (x(:) / pd(1))**pd(3)
       res(:) = exp(pd(2) * log(res(:)) - b - res(:)) / x(:)

       ! Lognormal dX/dx
    else if (md == 4) then
       a = 1.0_R_ / log(pd(2))
       res(:) = x(:) / pd(1)
       res(:) = exp(-0.5_R_ * (log(res(:)) * a)**2) * a / (sqrt(PI2_) * x(:))

       ! Lognormal dY/dr, Y=X**3
    else if (md == 5) then
       a = 1.0_R_ / log(pd(2))
       b = -4.5_R_ * log(pd(2))**2
       res(:) = x(:) / pd(1)
       res(:) = exp(-0.5_R_ * (log(res(:)) * a)**2 + b) * a / (sqrt(PI2_) * x(:) * res(:)**3)

       ! Others
    else
       res = 0.0_R_
    end if

  end function stat_probDens_1R


  !+
  ! Returns a normalized profile weights synthesized by multiple Gamma distributions
  !-
  function stat_prof_Gamma_1R(frac, alp, xave, xgrd) result(prof)

    real(R_), intent(in) :: frac(:) ! (nk) fractions of each kernel
    real(R_), intent(in) ::  alp(:) ! (nk) alpha (scaling parameter)
    real(R_), intent(in) :: xave(:) ! (nk) X mean of Gamma distribution
    real(R_), intent(in) :: xgrd(:) ! (nx+1) X grid values
    real(R_) :: prof(size(xgrd)-1)  ! (nx) normalized profile (weighting function)
    real(R_) :: fknl(size(xgrd)-1, 0:0) ! normalized kernel, integral(x=xgrd(ibin),xgrd(ibin+1)){PDF(x)}
    real(R_) :: bet
    integer  :: iker

    ! Sum up for the kernels
    prof(:) = 0.0_R_
    do iker = 1, size(alp)
       bet = xave(iker) / real(alp(iker), R_)
       fknl(:,:) = stat_Gamma_kernels_2R(alp(iker), bet, xgrd, 0)
       prof(:) = prof(:) + frac(iker) / sum(frac) * fknl(:,0)
    end do

    ! Normalize
    bet = 1.0_R_ / sum(prof)
    prof(:) = prof(:) * bet

  end function stat_prof_Gamma_1R


  !+
  ! Single-parameter regression
  !  modeling by y = a * x + b
  !-
  subroutine stat_sReg(x, y, a, b, cor, cov, xrmse, yrmse, xmean, ymean, xsig, ysig)

    real(R_), intent(in)  :: x(:), y(:) ! tabulated (x,y) data
    real(R_), intent(out) :: a, b       ! regression coefficients
    real(R_), intent(out), optional :: cor, cov     ! correlation & covariance
    real(R_), intent(out), optional :: xrmse, yrmse ! rms errors of x & y
    real(R_), intent(out), optional :: xmean, ymean ! averages of x & y
    real(R_), intent(out), optional :: xsig,  ysig  ! standard deviation of x & y
    integer  :: ndat
    real(R_) :: xmean1, ymean1, xsig1, ysig1, xrmse1, yrmse1, cov1, cor1

    ! Analyze
    call stat_coVec(x, y, cor1, cov1, xmean1, ymean1, xsig1, ysig1)
    ndat = size(y)
    a = 0.0_R_
    if (xsig1**2 > 0.0_R_) a = cov1 / xsig1**2
    b = ymean1 - a * xmean1
    xrmse1 = 0.0_R_
    yrmse1 = 0.0_R_
    if (ndat >= 3) yrmse1 = ysig1 * sqrt((1.0_R_ - cor1**2) * (ndat - 1) / (ndat - 2))
    if (ysig1 > 0.0_R_) xrmse1 = yrmse1 / ysig1 * xsig1

    ! Copy output values (if needed)
    if (present(cor)) cor = cor1
    if (present(cov)) cov = cov1
    if (present(xrmse)) xrmse = xrmse1
    if (present(yrmse)) yrmse = yrmse1
    if (present(xmean)) xmean = xmean1
    if (present(ymean)) ymean = ymean1
    if (present(xsig) ) xsig  = xsig1
    if (present(ysig) ) ysig  = ysig1

  end subroutine stat_sReg


  !+
  ! Approximate 5% point of T-distributoin
  !-
  function stat_Tdist_5p_R(nfree) result(res) 

    integer, intent(in)  :: nfree
    real(R_) :: res
    res = (1.4_R_ + 1.0_R_ / real(nfree, R_))**2

  end function stat_Tdist_5p_R

end module hparx_stat
