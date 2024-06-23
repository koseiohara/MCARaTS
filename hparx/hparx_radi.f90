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
! Library for utilities of the atmospheric radiation
!-
module hparx_radi

  use globals, only : R_, RD_, RSML_, RDSML_, PI_, FRAC13_, FRAC23_, GASCON_, AVOGAD_, &
       CSPEED_, PLANCK_, BOLTZ_
  use hparx_math, only : gaussLegen, effAve_exp_R, intp_Akima_coef_1R, intp_cubic_R, &
       intp_Hermite3_coef
  use hparx_vecmat, only : gridIdx_bin_I, gridIdx_loc
  use hparx_stat, only : stat_basic
  use hparx_rand, only : mseq_rand_R
  implicit none
  private

  ! Public procedures
  public :: autoMaxIntens_R
  public :: bcm_params
  public :: bcm_gen1D_1R
  public :: bcm_gen2D_2R
  public :: fresnelRef0
  public :: fresnelRef1
  public :: planck_R
  public :: planck_tvec_1R
  public :: planck_vecs
  public :: plkTab_init
  public :: plkTab_final
  public :: plkTab_intp_R
  public :: plkTab_intp_1R
  public :: randScat_HG
  public :: scatOpt_MieOpt_read
  public :: scatOpt_phys2opt
  public :: scatPF_biAsym
  public :: scatPF_asym_R
  public :: scatPF_distFunc
  public :: scatPF_HG_R
  public :: scatPF_HG_asym_1R
  public :: scatPF_HG_ang_1R
  public :: scatPF_intp
  public :: scatPF_Legendre
  public :: scatPF_normal
  public :: scatPF_normal2
  public :: scatPF_reconst
  public :: scatPF_triAsym
  public :: scatPF_truncDFlat_HG
  public :: scatPF_truncDFlat
  public :: scatPF_truncVert

  ! Constants for the Planck function
  real(R_), parameter :: Plk_C1 = CSPEED_**2 * (PLANCK_ * 2.0e24_R_)      ! = 1.191042759e8
  real(R_), parameter :: Plk_C2 = CSPEED_ * (PLANCK_ / BOLTZ_ * 1.0e6_R_) ! = 1.4387751359e4

  ! Private module variables: Band-averaged Planck function table
  real(R_), save, allocatable :: Plk_ptab(:,:) ! Planck functions (W/m^2/sr/micron)
  real(R_), save :: Plk_tmin
  real(R_), save :: Plk_tmax
  real(R_), save :: Plk_tfac
  integer,  save :: Plk_ntmp
  integer,  save :: Plk_nwav

contains

  !+
  ! Automatically-determined max intensity value, with excluding hot spots
  !-
  function autoMaxIntens_R(vec, flim, fmax) result(res)

    real(R_), intent(in) :: vec(:) ! vector data (only non-negative data are used)
    real(R_), intent(in), optional :: flim ! factor for threshold, used for excluding hot spots
    real(R_), intent(in), optional :: fmax ! factor for the max
    real(R_) :: res ! max value of vec(:) excluding hot spots
    real(R_) :: xmin, xmax, xmn, xsd, xsk, xkr, xfrc, xlim
    real(R_), parameter :: FLIM0 = 2.0_R_
    real(R_), parameter :: FMAX0 = 1.0_R_
    real(R_), parameter :: FACMIN = 1.3_R_
    
    call stat_basic(vec, xmin, xmax, xmn, xsd, xsk, xkr, (vec(:) >= 0.0_R_), xfrc)
    if (present(flim)) then
       xlim = xmn + xsd * flim
    else
       xlim = xmn + xsd * FLIM0
    end if
    call stat_basic(vec, xmin, xmax, xmn, xsd, xsk, xkr, (vec(:) >= 0.0_R_ .and. vec(:) < xlim), xfrc)
    if (present(fmax)) then
       res = max(FACMIN * xmn, xmn + xsd * fmax)
    else
       res = max(FACMIN * xmn, xmn + xsd * FMAX0)
    end if
    
  end function autoMaxIntens_R


  !+
  ! Bounded Cascade Model (BCM), 1-D version
  !-
  function bcm_gen1D_1R(ave, f, c, n) result(dat)

    real(R_), intent(in) :: ave ! average
    real(R_), intent(in) :: f   ! fractal parameter
    real(R_), intent(in) :: c   ! C in Cahalan et al. (1994, JAS)
    integer,  intent(in) :: n   ! # of data points (should be 2**P)
    real(R_) :: dat(n), dat1, frac, tmp0, tmp1
    integer  :: ntrn, i0, i1, ibin, ip2, itrn, np2

    ntrn = int(log(n + 0.1_R_) / log(2.0_R_))
    dat(1) = ave
    do itrn = 1, ntrn
       frac = f * c**real(itrn - 1, R_)
       np2  = 2**(itrn - 1)
       ibin = 2**(ntrn - itrn)
       do ip2 = 1, np2
          i0 = 2 * ibin * (ip2 - 1) + 1
          i1 = i0 + ibin
          dat1 = dat(i0)
          if (mseq_rand_R() < 0.5_R_) then
             tmp0 = (1.0_R_ + frac) * dat1
             tmp1 = (1.0_R_ - frac) * dat1
          else
             tmp0 = (1.0_R_ + frac) * dat1
             tmp1 = (1.0_R_ - frac) * dat1
          end if
          dat(i0) = tmp0
          dat(i1) = tmp1
       end do
    end do

  end function bcm_gen1D_1R


  !+
  ! Bounded Cascade Model (BCM), 2-D version
  !-
  function bcm_gen2D_2R(ave, f, c, n) result(dat)

    real(R_), intent(in) :: ave ! average
    real(R_), intent(in) :: f   ! fractal parameter
    real(R_), intent(in) :: c   ! C in Cahalan et al. (1994, JAS)
    integer,  intent(in) :: n   ! # of data points
    real(R_) :: dat(n, n), dat1, frac, tmp0, tmp1
    integer  :: ntrn, ip2x, ip2y, i0x, i0y, i1x, i1y, ibin, itrn, np2

    ntrn = int(log(n + 0.1_R_) / log(2.0_R_)) ! # of transfer processes
    dat(1, 1) = ave
    do itrn = 1, ntrn
       frac = f * c**real(itrn - 1, R_)
       np2  = 2**(itrn - 1)
       ibin = 2**(ntrn - itrn)
       do ip2y = 1, np2
          i0y = 2 * ibin * (ip2y - 1) + 1
          i1y = i0y + ibin
          do ip2x = 1, np2
             i0x = 2 * ibin * (ip2x - 1) + 1
             i1x = i0x + ibin
             dat1 = dat(i0x, i0y)
             if (mseq_rand_R() < 0.5_R_) then
                tmp0 = (1.0_R_ + frac) * dat1
                tmp1 = (1.0_R_ - frac) * dat1
             else
                tmp0 = (1.0_R_ + frac) * dat1
                tmp1 = (1.0_R_ - frac) * dat1
             end if
             dat(i0x, i0y) = tmp0
             dat(i1x, i0y) = tmp1
             dat(i0x, i1y) = tmp1
             dat(i1x, i1y) = tmp0
          end do
       end do
    end do

  end function bcm_gen2D_2R


  !+
  ! Get parameters of Bounded Cascade Model (BCM)
  !-
  subroutine bcm_params(xlogm, xlogs, c, f, xi, ave, del)

    real(R_), intent(in)  :: xlogm, xlogs ! average & sigma of logX
    real(R_), intent(in)  :: c   ! C in Cahalan et al. (1994)
    real(R_), intent(out) :: f   ! fractal parameter
    real(R_), intent(out) :: xi  ! reduction factor
    real(R_), intent(out) :: ave ! average
    real(R_), intent(out) :: del ! log(ave) - xlogm
    real(R_), parameter :: RATTH = 0.01_R_, FACTH = 0.99_R_
    integer,  parameter :: NMAX = 100
    integer,  parameter :: NTAB = 100
    real(R_), save :: flut(NTAB), slut(NTAB) ! LUT of S=mean(log(X)) wrt f
    integer,  save :: init = 0
    real(R_) :: f0, f1, fbin, fctr, fac10, ratio, spai, sum, tmp
    integer  :: itab, i

    ! Make LUT of S(f)
    fac10 = 1.0_R_ / log(10.0_R_)
    if (init == 0) then
       fbin = 0.9_R_ / real(NTAB, R_)
       flut(1) = 0.0_R_
       slut(1) = 0.0_R_
       do itab = 2, NTAB
          f0 = fbin * (itab - 1)
          flut(itab) = f0 
          sum = 0.0_R_
          ratio = 1.0_R_
          i = 0
          do
             if (i >= NMAX .or. ratio <= RATTH) exit
             f1 = f0 * c**i
             tmp = log((1.0_R_ + f1) / (1.0_R_ - f1)) * fac10 * 0.5_R_
             sum = sum + tmp * tmp
             ratio = tmp * tmp / sum
             i = i + 1
          end do
          slut(itab) = sqrt(sum)
       end do
       init = 1
    end if

    ! S ===> f, xi
    itab = gridIdx_bin_I(slut, 1, NTAB, xlogs)
    f = (xlogs - slut(itab)) / (slut(itab + 1) - slut(itab)) 
    f = flut(itab) + (flut(itab + 1) - flut(itab)) * f
    spai = 1.0_R_
    fctr = 0.0_R_
    i = 0
    do
       if (i >= NMAX .or. fctr >= FACTH) exit
       f1 = f * c**i
       fctr = 1.0_R_ - f1 * f1
       spai = spai * fctr
       i = i + 1
    end do
    xi = sqrt(spai)
    del = -log(xi) * fac10
    ave = 10.0_R_**(del + xlogm)

  end subroutine bcm_params


  !+
  ! Fresnel reflection for a real refractive index
  !-
  subroutine fresnelRef0(rr, uzi, uzt, rhov, rhoh, rho)

    real(R_), intent(in)  :: rr   ! real part of refractive index
    real(R_), intent(in)  :: uzi  ! cosine of incident   vector nadir angle (positive)
    real(R_), intent(out) :: uzt  ! cosine of refraction vector nadir angle (positive)
    real(R_), intent(out) :: rhoh ! V-component reflectance
    real(R_), intent(out) :: rhov ! H-component reflectance
    real(R_), intent(out) :: rho  ! average reflectance
    real(R_)  :: g, g2

    ! Ill input
    if (uzi <= 0.0_R_ .or. rr <= 0.0_R_) then
       uzt  =  0.0_R_
       rhov = -1.0_R_
       rhoh = -1.0_R_
       rho  = -1.0_R_
    else
       g2 = rr * rr + uzi * uzi - 1.0_R_

       ! Partial reflection
       if (g2 > RSML_) then
          g = sqrt(g2)
          uzt = g / rr
          rhov = (g - uzi) / (g + uzi)
          rhoh = rhov * (uzi * (g + uzi) - 1.0_R_) / (uzi * (g - uzi) + 1.0_R_)
          rhov = rhov * rhov
          rhoh = rhoh * rhoh
          rho  = 0.5_R_ * (rhov + rhoh)

          ! 100% reflection
       else
          uzt  = 0.0_R_
          rhov = 1.0_R_
          rhoh = 1.0_R_
          rho  = 1.0_R_
       end if
    end if

  end subroutine fresnelRef0


  !+
  ! Fresnel reflection for a complex refractive index
  !-
  subroutine fresnelRef1(rr, ri, uzi, uzt, rhov, rhoh, rho)

    real(R_), intent(in)  :: rr   ! real      part of refractive index
    real(R_), intent(in)  :: ri   ! imaginary part of refractive index (should be >= 0)
    real(R_), intent(in)  :: uzi  ! cosine of incident vector   nadir angle (positive)
    real(R_), intent(out) :: uzt  ! cosine of refraction vector nadir angle (positive)
    real(R_), intent(out) :: rhoh ! V-component reflectance
    real(R_), intent(out) :: rhov ! H-component reflectance
    real(R_), intent(out) :: rho  ! average reflectance
    real(R_)  :: g2, ri2, rr2, u, uzi2, v, w1, w2, w3, wa

    ! Ill input
    if (uzi <= 0.0_R_ .or. rr <= 0.0_R_) then
       uzt  =  0.0_R_
       rhov = -1.0_R_
       rhoh = -1.0_R_
       rho  = -1.0_R_

    else
       uzi2 = uzi * uzi
       rr2 = rr * rr
       ri2 = ri * ri
       g2 = rr2 + uzi2 - 1.0_R_

       ! Partial reflection
       if (g2 > RSML_) then
          uzt = sqrt(g2) / rr
          w1 = rr2 - ri2
          w2 = 2.0_R_ * rr * abs(ri)
          w3 = g2 - ri2
          wa = sqrt(w3 * w3 + w2 * w2)
          u = sqrt(0.5_R_ * abs(wa + w3))
          v = sqrt(0.5_R_ * abs(wa - w3))
          rhov = ((uzi - u)**2 + v*v) / ((uzi + u)**2 + v*v)
          rhoh =   ((w1 * uzi - u)**2 + (w2 * uzi - v)**2) &
               & / ((w1 * uzi + u)**2 + (w2 * uzi + v)**2)
          rho = 0.5_R_ * (rhov + rhoh)

          ! 100% reflection
       else
          uzt  = 0.0_R_
          rhov = 1.0_R_
          rhoh = 1.0_R_
          rho  = 1.0_R_
       end if
    end if

  end subroutine fresnelRef1


  !+
  ! Planck function
  !-
  function planck_R(tmp, wav) result(res) 

    real(R_),  intent(in) :: tmp  ! temperature (K)
    real(R_),  intent(in) :: wav  ! wavelength (micron)
    real(R_)  :: res ! result, Plank function (W/m2/sr/micron)
    real(R_), parameter :: AX = 15.0_R_
    real(R_), parameter :: ENAX = 3.0590232050182e-7_R_ ! = exp(-AX)
    real(R_) :: a
    
    a = Plk_C2 / (wav * tmp)
    if (a < 1.0e-4_R_) then  ! long wavelength approx.
       res = Plk_C1 * tmp / (Plk_C2 * wav**4)
    else if (a > AX) then    ! short wavelength approx.
       res = Plk_C1 / wav**5 * ENAX * exp(-a + AX)
    else                     ! usual case
       res = Plk_C1 / (wav**5 * (exp(a) - 1.0_R_))
    end if

  end function planck_R


  !+
  ! Planck functions for a temperature vector
  !-
  function planck_tvec_1R(tmp, wav) result(res) 

    real(R_),  intent(in) :: tmp(:) ! temperature (K)
    real(R_),  intent(in) :: wav    ! wavelength (micron)
    real(R_)  :: res(size(tmp))     ! result, Plank function (W/m2/sr/micron)
    real(R_), parameter :: AX = 15.0_R_
    real(R_), parameter :: ENAX = 3.0590232050182e-7_R_ ! = exp(-AX)
    real(R_) :: amin, amax
    
    amin = Plk_C2 / (wav * maxval(tmp))
    amax = Plk_C2 / (wav * minval(tmp))
    if (amax < 1.0e-4_R_) then  ! long wavelength approx.
       res(:) = Plk_C1 * tmp(:) / (Plk_C2 * wav**4)
    else if (amin > AX) then    ! short wavelength approx.
       res(:) = Plk_C1 / wav**5 * ENAX * exp(-(Plk_C2 / (wav * tmp(:))) + AX)
    else                        ! usual case
       res(:) = Plk_C1 / (wav**5 * (exp(Plk_C2 / (wav * tmp(:))) - 1.0_R_))
    end if

  end function planck_tvec_1R


  !+
  ! Calculate vectors of Planck functions & their derivatives wrt temperature
  !-
  subroutine planck_vecs(wav, tmp, bplk, dplk)

    real(R_), intent(in)  :: wav(:)  ! wavelength (micron)
    real(R_), intent(in)  :: tmp(:)  ! temperature T (K)
    real(R_), intent(out) :: bplk(:) ! Planck function, B (W/m2/sr/micron)
    real(R_), intent(out) :: dplk(:) ! derivative, dlnB/dT = dB/dT/B (/K)
    real(R_), parameter :: AX = 15.0_R_
    real(R_), parameter :: ENAX = 3.0590232050182e-7_R_ ! = exp(-AX)
    real(R_) :: a
    integer  :: i

    do i = 1, size(wav, 1)
       a = Plk_C2 / (wav(i) * tmp(i))
       if (a < 1.0e-4_R_) then  ! long wavelength approx.
          bplk(i) = Plk_C1 * tmp(i) / (Plk_C2 * wav(i)**4)
          dplk(i) = (1.0_R_ + a) / tmp(i)
       else if (a > AX) then    ! short wavelength approx.
          bplk(i) = Plk_C1 / wav(i)**5 * ENAX * exp(-a + AX)
          dplk(i) = a / tmp(i)
       else                     ! usual case
          bplk(i) = Plk_C1 / (wav(i)**5 * (exp(a) - 1.0_R_))
          dplk(i) = Plk_C2 / (wav(i) * tmp(i)**2 * (1.0_R_ - exp(-a)))
       end if
    end do

  end subroutine planck_vecs


  !+
  ! Initialize Plank function table for monochromatic wavelengths or spectral bands
  !-
  subroutine plkTab_init(tmin, tmax, ntmp, wav, dwav, nqmax)

    real(R_), intent(in) :: tmin, tmax
    integer,  intent(in) :: ntmp
    real(R_), intent(in) :: wav(:)
    real(R_), intent(in) :: dwav(:)
    integer,  intent(in) :: nqmax
    real(R_)  :: xq(nqmax), wq(nqmax)
    real(R_) :: dgam, w, t, dwavs, sumb, w0
    integer  :: nqold, nq, ns, iq, is, itmp, iwav

    ! Allocate
    if (allocated(Plk_ptab)) call plkTab_final()
    Plk_tmin = tmin
    Plk_tmax = tmax
    Plk_tfac = ntmp / (tmax - tmin)
    Plk_ntmp = ntmp
    Plk_nwav = size(wav)
    allocate (Plk_ptab(0:Plk_ntmp, Plk_nwav))

    ! Make table
    nqold = 0
    do iwav = 1, Plk_nwav
       dgam = dwav(iwav) / wav(iwav)
       nq = int(sqrt(dgam) * nqmax)
       if (nq <= 1) then
          w = wav(iwav) + 0.5_R_ * dwav(iwav)
          do itmp = 0, ntmp
             t = tmin + (tmax - tmin) * itmp / ntmp
             Plk_ptab(itmp, iwav) = planck_R(t, w)
          end do
       else
          if (dgam < 1.0_R_) then
             ns = 1
             nq = int(sqrt(dgam) * nqmax)
          else
             ns = int(sqrt(dgam))
             nq = nqmax
          end if
          if (nq /= nqold) then
             call gaussLegen(nq, xq, wq)
             nqold = nq
          end if
          dwavs = dwav(iwav) / ns
          xq(1:nq) = dwavs * 0.5_R_ * (1.0_R_ + xq(1:nq)) ! [0, dwavs]
          do itmp = 0, ntmp
             t = tmin + (tmax - tmin) * itmp / ntmp
             sumb = 0.0_R_
             do is = 1, ns
                w0 = wav(iwav) + dwavs * (is - 1)
                do iq = 1, nq
                   w = w0 + xq(iq)
                   sumb = sumb + wq(iq) * planck_R(t, w)
                end do
             end do
             Plk_ptab(itmp, iwav) = sumb / ns
          end do
       end if
    end do

  end subroutine plkTab_init


  !+
  ! Finalize Planck function table, Plk_*
  !-
  subroutine plkTab_final() 

    if (allocated(Plk_ptab)) deallocate (Plk_ptab)

  end subroutine plkTab_final


  !+
  ! Interpolate Planck function at a wavelength/band for a temperature value
  !-
  function plkTab_intp_R(t, iwav) result(b)

    real(R_), intent(in) :: t    ! temperature (K)
    integer,  intent(in) :: iwav ! wavelength index
    real(R_)  :: b ! result, interpolated Planck function
    integer   :: itmp
    real(R_)  :: x, a

    if (t <= Plk_tmin) then ! too cold
       b = Plk_ptab(0, iwav)
    else if (t >= Plk_tmax) then ! too hot
       b = Plk_ptab(Plk_ntmp, iwav)
    else ! normal case
       x = (t - Plk_tmin) * Plk_tfac
       itmp = min(Plk_ntmp - 1, int(x))
       a = x - itmp
       b = (1.0_R_ - a) * Plk_ptab(itmp, iwav) + a * Plk_ptab(itmp + 1, iwav)
    end if

  end function plkTab_intp_R


  !+
  ! Interpolate Planck functions at multiple wavelengths/bands for a temperature value
  !-
  function plkTab_intp_1R(t, iwmin, iwmax) result(b)

    real(R_), intent(in) :: t     ! temperature (K)
    integer,  intent(in) :: iwmin ! min of wavelength index
    integer,  intent(in) :: iwmax ! max of wavelength index
    real(R_)  :: b(iwmax - iwmin + 1) ! result, interpolated Planck functions
    integer   :: itmp
    real(R_)  :: x, a

    if (t <= Plk_tmin) then ! too cold
       b(:) = Plk_ptab(0, iwmin:iwmax)
    else if (t >= Plk_tmax) then ! too hot
       b(:) = Plk_ptab(Plk_ntmp, iwmin:iwmax)
    else ! normal case
       x = (t - Plk_tmin) * Plk_tfac
       itmp = min(Plk_ntmp - 1, int(x))
       a = x - itmp
       b(:) = (1.0_R_ - a) * Plk_ptab(itmp, iwmin:iwmax) + a * Plk_ptab(itmp + 1, iwmin:iwmax)
    end if

  end function plkTab_intp_1R


  !+
  ! Determine a random scattering angle for Henyey-Greenstein phase function
  !-
  subroutine randScat_HG(g, xi, cosq, sinq)

    real(R_), intent(in)  :: g          ! asymmetry factor
    real(R_), intent(in)  :: xi         ! a random number
    real(R_), intent(out) :: sinq, cosq ! sin/cos of scattering angle
    !// For xi = 0, cosq = +1 (forward scattering)
    !   For xi = 1, cosq = -1 (backscattering)
    !// Aug. 13, 2016, HI : Reversed for consistency with other randScat_* routines
    real(R_) :: xi1, gg

    xi1 = 1.0_R_ - 2.0_R_ * xi
    gg = g
    if (abs(gg) < 1.0e-6_R_) gg = sign(1.0e-6_R_, g)
    !cosq = (1.0_R_ + gg**2 - ((1.0_R_ - gg**2) / (1.0_R_ - gg * xi1))**2) / (2.0_R_ * gg)
    cosq = (1.0_R_ + gg**2 - ((1.0_R_ - gg**2) / (1.0_R_ + gg * xi1))**2) / (2.0_R_ * gg)
    cosq = max(-1.0_R_, min(1.0_R_, cosq))
    sinq = sqrt(1.0_R_ - cosq**2)

  end subroutine randScat_HG


  !+
  ! Import Mie optical properties from the database file
  !-
  subroutine scatOpt_MieOpt_read(iu, ndisp, rvdry, redry, rhodry, rvwet, rewet, rhowet, &
       & qe0, qe, qa, omg, asym, nang, npol, ang, qq)

    integer,  intent(in)  :: iu
    integer,  intent(out) :: ndisp
    real(R_), intent(out) :: rvdry(:), redry(:), rhodry(:)
    real(R_), intent(out) :: rvwet(:), rewet(:), rhowet(:)
    real(R_), intent(out) :: qe0(:), qe(:), qa(:), omg(:), asym(:)
    integer,  intent(out) :: nang, npol
    real(R_), intent(out) :: ang(:), qq(:,:,:)
    integer  :: idisp, iang, ipol

    ! Read in file
    read (iu, *)
    read (iu, *) 
    read (iu, *) ndisp
    read (iu, *)
    do idisp = 1, ndisp
       read (iu, *) rvdry(idisp), redry(idisp), rhodry(idisp), &
            &  rvwet(idisp), rewet(idisp), rhowet(idisp)
    end do
    read (iu, *)
    do idisp = 1, ndisp
       read (iu, *) qe0(idisp), qe(idisp), qa(idisp), omg(idisp), asym(idisp)
    end do
    read (iu, *)
    read (iu, *) nang, npol
    do idisp = 1, ndisp
       read (iu, *)
       do iang = 1, nang
          read (iu, *) ang(iang), (qq(iang, ipol, idisp), ipol = 1, npol)
       end do
    end do

  end subroutine scatOpt_MieOpt_read


  !+
  ! Convert volume and effictive radius of particles to extinction and absorption coefficients
  !-
  subroutine scatOpt_phys2opt(cola, effa, re, qe, qa, nd, be, ba, id, rd)

    real(R_), intent(in)  :: cola  ! some column amount (m3/m2 or m2/m2) for particle volume
    real(R_), intent(in)  :: effa  ! effective radius (micron)
    real(R_), intent(in)  :: re(:) ! table of effective radii (micron)
    real(R_), intent(in)  :: qe(:) ! table of extinction efficiency (m2/m2)
    real(R_), intent(in)  :: qa(:) ! table of absorption efficiency (m2/m2)
    integer,  intent(in)  :: nd    ! # of table grids
    real(R_), intent(out) :: be    ! extinction coefficient (m2/m3) or optical thickness
    real(R_), intent(out) :: ba    ! absorption coefficient (m2/m3) or optical thickness
    integer,  intent(out) :: id    ! table grid index found
    real(R_), intent(out) :: rd    ! table interpolation factor
    real(R_) :: aqe, aqa
    real(R_), parameter :: FEXT = 3.0_R_ / 4.0_R_ * 1.0e+6_R_

    if (effa * 1.0e+30_R_ > cola) then
       call gridIdx_loc(re, effa, id, rd)
       if (effa > re(1) .and. effa < re(nd)) rd = rd * (re(id + 1) / effa) ! non-linear correction
       aqe = (1.0_R_ - rd) * qe(id) + rd * qe(id + 1)
       aqa = (1.0_R_ - rd) * qa(id) + rd * qa(id + 1)
       be = FEXT * cola / effa * aqe
       ba = FEXT * cola / effa * aqa
    else
       be = 0.0_R_
       ba = 0.0_R_
       id = 1
       rd = 0.0_R_
    end if

  end subroutine scatOpt_phys2opt


  !+
  ! Asymmetry factor of scattering phase function
  !-
  function scatPF_asym_R(cosa, sina, phsf) result(res)

    real(R_), intent(in) :: cosa(:), sina(:) ! cos & sin of angles
    real(R_), intent(in) :: phsf(:) ! phase functions
    real(R_) :: res, dx
    integer  :: nang, iang
    real(RD_) :: sump, sumg

    ! Initialize
    nang = size(phsf)
    sump = 0.0_RD_
    sumg = 0.0_RD_

    ! Integrate (trapezoidal)
    dx = abs(sina(nang - 1)**2 / (1.0_R_ - cosa(nang - 1)) - sina(nang)**2 / (1.0_R_ - cosa(nang)))
    sump = sump + dx * phsf(nang)
    sumg = sumg + dx * phsf(nang) * cosa(nang)
    do iang = nang - 1, 2, -1
       if (cosa(iang)**2 < 0.9_R_) then
          dx = abs(cosa(iang - 1) - cosa(iang + 1))
       else if (cosa(iang) > 0.0_R_) then
          dx = abs(sina(iang - 1)**2 / (1.0_R_ + cosa(iang - 1)) &
               & - sina(iang + 1)**2 / (1.0_R_ + cosa(iang + 1)))
       else
          dx = abs(sina(iang - 1)**2 / (1.0_R_ - cosa(iang - 1)) &
               & - sina(iang + 1)**2 / (1.0_R_ - cosa(iang + 1)))
       end if
       sump = sump + dx * phsf(iang)
       sumg = sumg + dx * phsf(iang) * cosa(iang)
    end do
    dx = abs(sina(1)**2 / (1.0_R_ + cosa(1)) - sina(2)**2 / (1.0_R_ + cosa(2)))
    sump = sump + dx * phsf(1)
    sumg = sumg + dx * phsf(1) * cosa(1)
    res = sumg / sump

  end function scatPF_asym_R


  !+
  ! Calculate asymmetry factors of scattering phase function for bi-sections
  !  forward and backward regions
  !-
  subroutine scatPF_biAsym(cosa, sina, phs, coss, g0fwd, g0bwd, &
       & g1fwd, g1bwd, g1, g2fwd, g2bwd, g2)

    real(R_), intent(in)  :: cosa(:), sina(:) ! sin & cos
    real(R_), intent(in)  :: phs(:)           ! phase function
    real(R_), intent(in)  :: coss             ! cos at the section separation point
    real(R_), intent(out) :: g0fwd, g0bwd     ! fractions of forward/backward parts
    real(R_), intent(out) :: g1fwd, g1bwd, g1 ! first moments (asymmetry factors)
    real(R_), intent(out) :: g2fwd, g2bwd, g2 ! second moments
    real(RD_) :: sum0, sum1, sum2
    integer   :: nang, iang, iangs
    real(R_)  :: dcosa, phss, rat, w

    ! Separation point & P interpolation
    nang = size(cosa)
    iangs = gridIdx_bin_I(cosa, 1, nang, coss)
    dcosa = cosa(iangs + 1) - cosa(iangs)
    if (dcosa > RSML_) then
       phss = phs(iangs)
    else
       rat = (coss - cosa(iangs)) / dcosa
       phss = phs(iangs) * (1.0_R_ - rat) + phs(iangs + 1) * rat
    end if

    ! Forward region
    sum0 = 0.0_R_
    sum1 = 0.0_R_
    sum2 = 0.0_R_
    do iang = 1, iangs - 1
       w = phs(iang) + phs(iang+1)
       if (cosa(iang) * cosa(iang+1) > 0.99_R_) then
          if (cosa(iang) > 0.0_R_) then
             w = w * abs(sina(iang+1)**2 / (1.0_R_ + cosa(iang+1)) &
                  & - sina(iang)**2 / (1.0_R_ + cosa(iang)))
          else
             w = w * abs(sina(iang)**2 / (1.0_R_ - cosa(iang)) &
                  & - sina(iang+1)**2 / (1.0_R_ - cosa(iang+1)))
          end if
       else
          w = w * abs(cosa(iang) - cosa(iang+1))
       end if
       sum0 = sum0 + w
       sum1 = sum1 + w * 0.5_R_ * (cosa(iang)    + cosa(iang+1))
       sum2 = sum2 + w * 0.5_R_ * (cosa(iang)**2 + cosa(iang+1)**2)
    end do
    w = phs(iangs) + phss
    w = w * abs(cosa(iang) - coss)
    sum0 = sum0 + w
    sum1 = sum1 + w * 0.5_R_ * (cosa(iang)    + coss)
    sum2 = sum2 + w * 0.5_R_ * (cosa(iang)**2 + coss**2)
    g0fwd = sum0
    if (sum0 > RDSML_) then
       g1fwd = sum1 / sum0
       g2fwd = sum2 / sum0
    else
       g1fwd = 0.0_R_
       g2fwd = 0.0_R_
    end if

    ! Backward region
    sum0 = 0.0_R_
    sum1 = 0.0_R_
    sum2 = 0.0_R_
    do iang = nang - 1, iangs + 1, -1
       w = phs(iang) + phs(iang+1)
       if (cosa(iang) * cosa(iang+1) > 0.99_R_) then
          if (cosa(iang) > 0.0_R_) then
             w = w * abs(sina(iang+1)**2 / (1.0_R_ + cosa(iang+1)) &
                  & - sina(iang)**2 / (1.0_R_ + cosa(iang)))
          else
             w = w * abs(sina(iang)**2 / (1.0_R_ - cosa(iang)) &
                  & - sina(iang+1)**2 / (1.0_R_ - cosa(iang+1)))
          end if
       else
          w = w * abs(cosa(iang) - cosa(iang+1))
       end if
       sum0 = sum0 + w
       sum1 = sum1 + w * 0.5_R_ * (cosa(iang)    + cosa(iang+1))
       sum2 = sum2 + w * 0.5_R_ * (cosa(iang)**2 + cosa(iang+1)**2)
    end do
    w = phs(iangs + 1) + phss
    w = w * abs(coss - cosa(iang+1))
    sum0 = sum0 + w
    sum1 = sum1 + w * 0.5_R_ * (cosa(iang+1)    + coss)
    sum2 = sum2 + w * 0.5_R_ * (cosa(iang+1)**2 + coss**2)
    g0bwd = sum0
    if (sum0 > RDSML_) then
       g1bwd = sum1 / sum0
       g2bwd = sum2 / sum0
    else
       g1bwd = 0.0_R_
       g2bwd = 0.0_R_
    end if

    ! Moments
    sum0 = g0fwd + g0bwd
    g0fwd = g0fwd / sum0
    g0bwd = g0bwd / sum0
    g1 = g0fwd * g1fwd + g0bwd * g1bwd
    g2 = g0fwd * g2fwd + g0bwd * g2bwd

  end subroutine scatPF_biAsym


  !+
  ! Probability Density Function (PDF) & Cumulative Distribution Function (CDF)
  !  of the scattering angle
  !-
  subroutine scatPF_distFunc(cosa, cosa_c, phs, pdf, cdf)

    real(R_), intent(in)    :: cosa(:)   ! cosA
    real(R_), intent(in)    :: cosa_c(:) ! 1 - cosA for cosA > 0, 1 + cosA for cosA < 0
    real(R_), intent(inout) :: phs(:)    ! phase function (output will be normalized)
    real(R_), intent(out)   :: pdf(:)    ! PDF, pdf(iang) = cdf(iang) - cdf(iang+1)
    real(R_), intent(out)   :: cdf(:)    ! CDF, decrease with increasing iang
    ! - Note: cdf is integrated in the backward order (pi to 0)
    !    cdf(1) = 1.0_R_, cdf(nang) = 0.0_R_, pdf(nang) = 0.0_R_
    integer  :: iang, nang
    real(R_) :: c1, c0, cc0, cc1, y1, y0, d
    real(RD_) :: sump

    ! Integrate P (in the backward order)
    nang = size(phs)
    sump = 0.0_RD_
    pdf(nang) = 0.0_R_
    cdf(nang) = 0.0_R_
    y1  = phs(nang)
    c1  = cosa(nang)
    cc1 = cosa_c(nang)
    do iang = nang - 1, 1, -1
       y0  = phs(iang)
       c0  = cosa(iang)
       cc0 = cosa_c(iang)
       if (c0 * c1 < 0.99_R_) then
          d = (y0 + y1) * abs(c0 - c1)
       else
          d = (y0 + y1) * abs(cc1 - cc0)
       end if
       sump = sump + d
       pdf(iang) = d
       cdf(iang) = sump
       y1  = y0
       c1  = c0
       cc1 = cc0
    end do

    ! Normalizations: (2*pi)/(4*pi)*sum{Q=0,pi;(P(Q)*sinQ*dQ)} = 1
    d = 1.0_R_ / real(sump, R_)
    phs(:) = phs(:) * d * 4.0_R_
    pdf(:) = pdf(:) * d
    cdf(:) = cdf(:) * d
    cdf(1) = 1.0_R_
    cdf(nang) = 0.0_R_

  end subroutine scatPF_distFunc


  !+
  ! A value of Henyey-Greenstein phase function at a specific direction
  !-
  function scatPF_HG_R(g, cosq) result(pfun)

    real(R_), intent(in) :: g    ! asymmetry factor
    real(R_), intent(in) :: cosq ! cos(scattering angle)
    real(R_) :: pfun

    pfun = 1.0_R_ + g**2 - 2.0_R_ * g * cosq
    pfun = (1.0_R_ - g**2) / (pfun * sqrt(pfun))

  end function scatPF_HG_R


  !+
  ! Vector of Henyey-Greenstein phase function values for different asymmetry factors
  !-
  function scatPF_HG_asym_1R(g, cosq) result(pfun) 

    real(R_), intent(in)  :: g(:)    ! asymmetry factors
    real(R_), intent(in)  :: cosq    ! cos(scattering angle)
    real(R_) :: pfun(size(g, 1))     ! phase function

    pfun(:) = 1.0_R_ + g(:)**2 - 2.0_R_ * g(:) * cosq
    pfun(:) = (1.0_R_ - g(:)**2) / (pfun(:) * sqrt(pfun(:)))

  end function scatPF_HG_asym_1R


  !+
  ! Get Henyey-Greenstein phase function at specified angles
  !-
  function scatPF_HG_ang_1R(g, cosq) result(pfun)

    real(R_), intent(in)  :: g       ! asymmetry factor
    real(R_), intent(in)  :: cosq(:) ! cos(scattering angle)
    real(R_) :: pfun(size(cosq, 1))  ! phase function (the same array size as cosq(:))

    pfun(:) = 1.0_R_ + g**2 - 2.0_R_ * g * cosq(:)
    pfun(:) = (1.0_R_ - g**2) / (pfun(:) * sqrt(pfun(:)))

  end function scatPF_HG_ang_1R


  !+
  ! Interpolate phase functions by cubic polynomial
  !-
  subroutine scatPF_intp(ang, phs, nang, angi, nangi, phsi)

    real(R_), intent(in) :: ang(:), phs(:)
    real(R_), intent(in) :: angi(:)
    integer,  intent(in) :: nang, nangi
    real(R_), intent(out) :: phsi(:)
    real(R_) :: x, y
    integer  :: iangi, iang
    real(R_), parameter :: PHSMIN = 1.0e-35_R_

    do iangi = 1, nangi
       x = angi(iangi)
       iang = gridIdx_bin_I(ang, 1, nang, x) ! sequential search is better
       iang = max(2, min(nang - 2, iang))
       !y = intp_cubic_R(x, ang(iang-1:iang+2), log(max(PHSMIN, phs(iang-1:iang+2))))
       !phsi(iangi) = exp(y)
       y = intp_cubic_R(x, ang(iang-1:iang+2), phs(iang-1:iang+2))
       phsi(iangi) = max(PHSMIN, y)
    end do

  end subroutine scatPF_intp


  !+
  ! Legendre polynomial expansion:
  !  Compute the expansion coefficients of a phase function
  !  using the Gaussian quadrature for each angle interval
  !-
  subroutine scatPF_Legendre(ang, phs, gmom, nquad, angmin, angmax, gmom0)

    real(R_), intent(in)  :: ang(:)  ! angle (radian)
    real(R_), intent(in)  :: phs(:)  ! phase function (not normalized)
    real(R_), intent(out) :: gmom(:) ! moments of phase function (normalized)
    integer,  intent(in),  optional :: nquad  ! # of quadrature points in each interval 
    real(R_), intent(in),  optional :: angmin ! min angles (radian)
    real(R_), intent(in),  optional :: angmax ! max angles (radian)
    real(R_), intent(out), optional :: gmom0  ! 0th moment not normalized
    real(RD_) :: sumg(0:size(gmom)), sumq(0:size(gmom))
    real(R_)  :: c0(size(gmom)), c1(size(gmom))
    real(R_)  :: p0, p1, p2, angmin1, angmax1, xmin, xmax, dx, x, y, ax, c, d, cosx, wy
    integer   :: nmom, iangmin, iangmax, nang, iang, imom, ii, nq, iq
    integer,  parameter :: NQDEF = 10, NQMAX = 32
    real(R_), allocatable :: wquad(:), xquad(:), dydx(:)

    ! Initialize
    nang = size(ang)
    nmom = size(gmom)
    nq = NQDEF
    if (present(nquad)) nq = min(NQMAX, nquad)
    allocate (wquad(nq), xquad(nq), dydx(nang))
    call gaussLegen(nq, xquad, wquad)
    xquad(:) = (xquad(:) + 1.0_R_) * 0.5_R_
    dydx(:) = intp_Akima_coef_1R(ang, phs) ! coefficients for the Akima interpolation
    ii = -1
    do imom = 1, nmom
       ii = ii + 2
       c0(imom) = real(imom-1, R_) / real(imom, R_)
       c1(imom) = real(    ii, R_) / real(imom, R_)
    end do

    ! Angle range
    if (present(angmin)) then
       angmin1 = angmin
       iangmin = gridIdx_bin_I(ang, 1, nang, angmin)
    else
       angmin1 = ang(1)
       iangmin = 1
    end if
    if (present(angmax)) then
       angmax1 = angmax
       iangmax = gridIdx_bin_I(ang, 1, nang, angmax)
    else
       angmax1 = ang(nang)
       iangmax = nang - 1
    end if

    ! Loop for integration over angle intervals
    sumg(:) = 0.0_RD_
    do iang = iangmax, iangmin, -1  ! reverse order due to strong forward scattering

       ! Interval
       xmin = ang(iang)
       xmax = ang(iang+1)
       dx = max(RSML_, xmax - xmin) ! original interval
       if (iang == iangmin) xmin = angmin1
       if (iang == iangmax) xmax = angmax1

       ! Loop for Gaussian points
       sumq(:) = 0.0_RD_
       do iq = nq, 1, -1    ! integrate by Gaussian quadrature
          x = xmin + (xmax - xmin) * xquad(iq)
          cosx = cos(x)
          ax = (x - ang(iang)) / dx
          call intp_Hermite3_coef(phs(iang), phs(iang+1), dydx(iang)*dx, dydx(iang+1)*dx, c, d)
          y = phs(iang) + ax * (dydx(iang)*dx + ax * (c + ax * d)) ! Akima interpolation
          wy = 0.5_R_ * y * sin(x) * wquad(iq) * (xmax - xmin)
          p0 = 1.0_R_
          p1 = cosx
          sumq(0) = sumq(0) + wy
          sumq(1) = sumq(1) + wy * p1
          do imom = 2, nmom
             p2 = c1(imom) * cosx * p1 - c0(imom) * p0
             sumq(imom) = sumq(imom) + wy * p2
             p0 = p1
             p1 = p2
          end do
       end do
       sumg(:) = sumg(:) + sumq(:)
    end do

    ! Results
    gmom(1:nmom) = sumg(1:nmom) / sumg(0)  ! normalized moments
    if (present(gmom0)) gmom0 = sumg(0) ! non-normalized 0th moment
    deallocate (wquad, xquad, dydx)

  end subroutine scatPF_Legendre


  !+
  ! Caution! This is not good! Use scatPF_normal2 instead!
  !
  ! Normalize a phase function, and calculate the asymmetry parameter
  !-
  subroutine scatPF_normal(cosa, cosa_c, nang, phsf, asym)

    real(R_), intent(in)    :: cosa(:)   ! cos(angle)
    real(R_), intent(in)    :: cosa_c(:) ! (1 +/- cos)
    integer,  intent(in)    :: nang      ! # of angles
    real(R_), intent(inout) :: phsf(:)   ! phase function
    real(R_), intent(out)   :: asym      ! asymmetry parameter
    integer  :: iang
    real(R_) :: dx
    real(RD_) :: sump, sumg

    ! Zeros
    sump = 0.0_RD_
    sumg = 0.0_RD_

    ! Integrate (trapezoidal)
    dx = abs(cosa_c(nang - 1) - cosa_c(nang))
    sump = sump + dx * phsf(nang)
    sumg = sumg + dx * phsf(nang) * cosa(nang)
    do iang = nang - 1, 2, -1
       if (cosa(iang)**2 < 0.9_R_) then
          dx = abs(cosa(iang - 1) - cosa(iang + 1))
       else
          dx = abs(cosa_c(iang - 1) - cosa_c(iang + 1))
       end if
       sump = sump + dx * phsf(iang)
       sumg = sumg + dx * phsf(iang) * cosa(iang)
    end do
    dx = abs(cosa_c(1) - cosa_c(2))
    sump = sump + dx * phsf(1)
    sumg = sumg + dx * phsf(1) * cosa(1)

    ! Normalize
    phsf(:) = phsf(:) * (4.0_R_ / real(sump, R_))
    asym = sumg / sump

  end subroutine scatPF_normal


  !+
  ! Normalize a phase function, and calculate the asymmetry parameter
  !-
  subroutine scatPF_normal2(ang, sina, phsf, cosa, asym)

    real(R_), intent(in)    :: ang(:)  ! angles (radian)
    real(R_), intent(in)    :: sina(:) ! sin(angle)
    real(R_), intent(inout) :: phsf(:) ! phase function
    real(R_), intent(in),  optional :: cosa(:) ! cos(angle)
    real(R_), intent(out), optional :: asym    ! asymmetry parameter
    integer  :: iang, nang
    real(R_) :: dxy, dxy2
    real(RD_) :: sump, sumg

    ! Zeros
    nang = size(phsf)
    sump = 0.0_RD_
    sumg = 0.0_RD_

    ! Integrate (trapezoidal) for the body part
    do iang = nang - 1, 2, -1
       dxy = abs(ang(iang + 1) - ang(iang - 1)) * sina(iang) * phsf(iang)
       sump = sump + dxy
       if (present(cosa)) sumg = sumg + dxy * cosa(iang)
    end do

    ! Integrate for the tips
    dxy  = abs(ang(nang) - ang(nang - 1)) * sina(nang) * phsf(nang)
    dxy2 = abs(ang(2) - ang(1)) * sina(1) * phsf(1)
    sump = sump + (dxy + dxy2)
    if (present(cosa)) sumg = sumg + (dxy * cosa(nang) + dxy2 * cosa(1))

    ! Normalize
    phsf(:) = phsf(:) * (4.0_R_ / real(sump, R_))
    if (present(asym)) asym = sumg / sump

  end subroutine scatPF_normal2


  !+
  ! Reconstruction of phase function from Legendre expansion coefficients
  !  using the delta-M or delta-Legendre approximation
  !-
  subroutine scatPF_reconst(gmo, ngmo, cosang, fd, nang, phs, ngmo0)

    integer, intent(in) :: nang
    integer, intent(in) :: ngmo
    real(R_), intent(in) :: fd
    real(R_), intent(in) :: cosang(:) ! cos(angle)
    real(R_), intent(in) :: gmo(:)    ! moments
    real(R_), intent(out) :: phs(:)    ! phase function
    integer, intent(out) :: ngmo0
    real(RD_) :: sump
    real(R_) :: p0, p1, p2, x
    integer :: igmo, iang, i

    ! Max order
    do igmo = 1, ngmo
       if (gmo(igmo) < fd) exit
    end do
    ngmo0 = igmo - 1         ! max order of truncated series

    ! Construct
    do iang = 1, nang
       x = cosang(iang)
       p2 = 1.0_R_
       sump = (1.0_R_ - fd) * p2
       p1 = x
       i = 3
       if (ngmo0 >= 1) sump = sump + real(i, R_) * (gmo(1) - fd) * p1
       do igmo = 2, ngmo0
          p0 = (real(i, R_) / real(igmo, R_)) * x * p1 - (real(igmo-1, R_) / real(igmo, R_)) * p2
          i = i + 2
          sump = sump + real(i, R_) * (gmo(igmo) - fd) * p0
          p2 = p1
          p1 = p0
       end do
       phs(iang) = sump / (1.0_R_ - fd)
    end do

  end subroutine scatPF_reconst


  !+
  ! Calculate asymmetry factors of scattering phase function for three sections
  !-
  subroutine scatPF_triAsym(ang, phs, nang, angf, angb, g0fwd, g0mid, g0bwd, &
       & g1fwd, g1mid, g1bwd, g1, g2fwd, g2mid, g2bwd, g2)

    real(R_), intent(in)  :: ang(:)     ! angles (radian)
    real(R_), intent(in)  :: phs(:)     ! phase function
    real(R_), intent(in)  :: angf, angb ! angles (radian) at the section separation points
    real(R_), intent(out) :: g0fwd, g0mid, g0bwd     ! fractions of forward/backward parts
    real(R_), intent(out) :: g1fwd, g1mid, g1bwd, g1 ! first moments (asymmetry factors)
    real(R_), intent(out) :: g2fwd, g2mid, g2bwd, g2 ! second moments
    integer   :: nang, iang, iangf, iangb
    real(RD_) :: sum0, sum1, sum2
    real(R_)  :: cosa0, cosa1, dang, phsb, phsf, rat, sina0, sina1, w

    ! Separation points & P interpolation
    nang = size(ang)
    iangf = gridIdx_bin_I(ang, 1, nang, angf)
    iangb = gridIdx_bin_I(ang, 1, nang, angb)
    dang = ang(iangf + 1) - ang(iangf)
    if (dang > RSML_) then
       phsf = phs(iangf)
    else
       rat = (angf - ang(iangf)) / dang
       phsf = phs(iangf) * (1.0_R_ - rat) + phs(iangf + 1) * rat
    end if
    dang = ang(iangb + 1) - ang(iangb)
    if (dang > RSML_) then
       phsb = phs(iangb)
    else
       rat = (angb - ang(iangb)) / dang
       phsb = phs(iangb) * (1.0_R_ - rat) + phs(iangb + 1) * rat
    end if

    ! Forward region
    sum0 = 0.0_R_
    sum1 = 0.0_R_
    sum2 = 0.0_R_
    sina0 = sin(ang(1))
    cosa0 = cos(ang(1))
    do iang = 1, iangf - 1
       sina1 = sin(ang(iang + 1))
       cosa1 = cos(ang(iang + 1))
       w = (ang(iang + 1) - ang(iang)) * (phs(iang + 1) * sina1 + phs(iang) * sina0)
       sum0 = sum0 + w
       sum1 = sum1 + w * 0.5_R_ * (cosa0    + cosa1)
       sum2 = sum2 + w * 0.5_R_ * (cosa0**2 + cosa1**2)
       sina0 = sina1
       cosa0 = cosa1
    end do
    sina1 = sin(angf)
    cosa1 = cos(angf)
    w = (angf - ang(iangf)) * (phsf * sina1 + phs(iangf) * sina0)
    sum0 = sum0 + w
    sum1 = sum1 + w * 0.5_R_ * (cosa0    + cosa1)
    sum2 = sum2 + w * 0.5_R_ * (cosa0**2 + cosa1**2)
    g0fwd = sum0
    if (sum0 > RSML_) then
       g1fwd = sum1 / sum0
       g2fwd = sum2 / sum0
    else
       g1fwd = 0.0_R_
       g2fwd = 0.0_R_
    end if

    ! Backward region
    sum0 = 0.0_R_
    sum1 = 0.0_R_
    sum2 = 0.0_R_
    sina1 = sin(ang(nang))
    cosa1 = cos(ang(nang))
    do iang = nang - 1, iangb + 1, -1
       sina0 = sin(ang(iang))
       cosa0 = cos(ang(iang))
       w = (ang(iang+1) - ang(iang)) * (phs(iang + 1) * sina1 + phs(iang) * sina0)
       sum0 = sum0 + w
       sum1 = sum1 + w * 0.5_R_ * (cosa0    + cosa1)
       sum2 = sum2 + w * 0.5_R_ * (cosa0**2 + cosa1**2)
       sina1 = sina0
       cosa1 = cosa0
    end do
    sina0 = sin(angb)
    cosa0 = cos(angb)
    w = (ang(iangb+1) - angb) * (phs(iangb + 1) * sina1 + phsb * sina0)
    sum0 = sum0 + w
    sum1 = sum1 + w * 0.5_R_ * (cosa0    + cosa1)
    sum2 = sum2 + w * 0.5_R_ * (cosa0**2 + cosa1**2)
    g0bwd = sum0
    if (sum0 > RSML_) then
       g1bwd = sum1 / sum0
       g2bwd = sum2 / sum0
    else
       g1bwd = 0.0_R_
       g2bwd = 0.0_R_
    end if

    ! Mid region
    sum0 = 0.0_R_
    sum1 = 0.0_R_
    sum2 = 0.0_R_
    sina0 = sin(ang(iangb))
    cosa0 = cos(ang(iangb))
    sina1 = sin(angb)
    cosa1 = cos(angb)
    w = (angb - ang(iangb)) * (phsb * sina1 + phs(iangb) * sina0)
    sum0 = sum0 + w
    sum1 = sum1 + w * 0.5_R_ * (cosa0    + cosa1)
    sum2 = sum2 + w * 0.5_R_ * (cosa0**2 + cosa1**2)
    sina1 = sina0
    cosa1 = cosa0
    do iang = iangb - 1, iangf + 1, -1
       sina0 = sin(ang(iang))
       cosa0 = cos(ang(iang))
       w = (ang(iang+1) - ang(iang)) * (phs(iang + 1) * sina1 + phs(iang) * sina0)
       sum0 = sum0 + w
       sum1 = sum1 + w * 0.5_R_ * (cosa0    + cosa1)
       sum2 = sum2 + w * 0.5_R_ * (cosa0**2 + cosa1**2)
       sina1 = sina0
       cosa1 = cosa0
    end do
    sina0 = sin(angf)
    cosa0 = cos(angf)
    sina1 = sin(ang(iangf + 1))
    cosa1 = cos(ang(iangf + 1))
    w = (ang(iangf + 1) - angf) * (phs(iangf + 1) * sina1 + phsf * sina0)
    sum0 = sum0 + w
    sum1 = sum1 + w * 0.5_R_ * (cosa0    + cosa1)
    sum2 = sum2 + w * 0.5_R_ * (cosa0**2 + cosa1**2)
    g0mid = sum0
    if (sum0 > RSML_) then
       g1mid = sum1 / sum0
       g2mid = sum2 / sum0
    else
       g1mid = 0.0_R_
       g2mid = 0.0_R_
    end if

    ! Moments
    g0fwd = g0fwd / (g0fwd + g0mid + g0bwd)
    g0mid = g0mid / (g0fwd + g0mid + g0bwd)
    g0bwd = g0bwd / (g0fwd + g0mid + g0bwd)
    g1 = g0fwd * g1fwd + g0mid * g1mid + g0bwd * g1bwd
    g2 = g0fwd * g2fwd + g0mid * g2mid + g0bwd * g2bwd

  end subroutine scatPF_triAsym


  !+
  ! Calculate parameters for a truncation approximation,
  !  for Henyey-Greenstein phase function
  !-
  subroutine scatPF_truncDFlat_HG(g, cosqf_c, fdec, fbod, pfwd)

    real(R_), intent(in)  :: g       ! asymmetry factor
    real(R_), intent(in)  :: cosqf_c ! 1 - cos(scattering angle) at truncation point
    real(R_), intent(out) :: fdec    ! 1 - (delta fraction)
    real(R_), intent(out) :: fbod    ! fraction of body part
    real(R_), intent(out) :: pfwd    ! phase function at flat forward part
    real(R_) :: b, c, r, s0, s1

    b = 1.0_R_ + g**2
    r = sqrt(b - 2.0_R_ * g * (1.0_R_ - cosqf_c))
    c = (g + r - 1.0_R_) * (1.0_R_ - g**2) / (2.0_R_ * g)
    s0 = c / (r * (1.0_R_ - g))
    s1 = (b * s0 - c) / (2.0_R_ * g)

    fbod = 1.0_R_ - s0
    pfwd = (s0 - s1) * 4.0_R_ / cosqf_c**2
    fdec = 1.0_R_ - (2.0_R_ * s1 - (2.0_R_ - cosqf_c) * s0) / cosqf_c

  end subroutine scatPF_truncDFlat_HG


  !+
  ! Calculate parameters for a truncation approximation
  !-
  subroutine scatPF_truncDFlat(phs, cosa, cosa_c, iangf, fdec, fbod, pfwd)

    real(R_), intent(in)  :: phs(:)    ! phase function
    real(R_), intent(in)  :: cosa(:)   ! cos(angle)
    real(R_), intent(in)  :: cosa_c(:) ! (1 +/- cos)
    integer,  intent(in)  :: iangf     ! truncation angle index
    real(R_), intent(out) :: fdec      ! 1 - (delta fraction)
    real(R_), intent(out) :: fbod      ! fraction of body part
    real(R_), intent(out) :: pfwd      ! phase function at flat forward part
    integer  :: iang
    real(R_) :: w
    real(RD_) :: sum0, sum1, fdel

    ! No truncation
    if (iangf <= 1) then
       fdec = 1.0_R_
       fbod = 1.0_R_
       pfwd = 1.0_R_

       ! Delta and flat forward part
    else
       sum0 = 0.0_RD_
       sum1 = 0.0_RD_
       do iang = 1, iangf - 1
          if (abs(cosa(iang)) < 0.99_R_) then
             w = (phs(iang) + phs(iang + 1)) * abs(cosa(iang) - cosa(iang + 1))
          else
             w = (phs(iang) + phs(iang + 1)) * abs(cosa_c(iang) - cosa_c(iang + 1))
          end if
          sum0 = sum0 + w
          sum1 = sum1 + w * 0.5_R_ * (cosa(iang) + cosa(iang + 1))
       end do
       sum0 = sum0 * 0.25_RD_
       sum1 = sum1 * 0.25_RD_
       fbod = 1.0_RD_ - sum0
       pfwd = real(sum0 - sum1, R_) * 4.0_R_ / cosa_c(iangf)**2
       fdel = (2.0_RD_ * sum1 - (1.0_RD_ + cosa(iangf)) * sum0) / cosa_c(iangf)
       fdec = 1.0_RD_ - fdel
    end if

  end subroutine scatPF_truncDFlat


  !+
  ! Forward-end truncation approximation (FTA)
  !   Truncation of forward peak of phase function with a correction
  !   for conservativation of asymmetry factor.
  ! Note: cump() is in the reverse order (1 to 0)
  ! The renormalized phase function is
  !    P = 0                     for ang <= angtc
  !        phsf(ang) / (1 - ftc) for ang  > angtc
  !-
  subroutine scatPF_truncVert(cosa, cump, pdf, ang, g, ftmax, ft, ftc, gtc, iangtc, rat, angtc)

    real(R_), intent(in)  :: cosa(:)
    real(R_), intent(in)  :: cump(:)
    real(R_), intent(in)  :: pdf(:)
    real(R_), intent(in)  :: ang(:)
    real(R_), intent(in)  :: g
    real(R_), intent(in)  :: ftmax
    real(R_), intent(out) :: ft
    real(R_), intent(out) :: ftc
    real(R_), intent(out) :: gtc
    integer,  intent(out) :: iangtc
    real(R_), intent(out) :: rat
    real(R_), intent(out) :: angtc
    integer :: nang, iang, iangt
    real(RD_) :: sum0, sum1, gtc0, gtc1, gtc2, dft, ratc
    real(R_)  :: cumpmin

    ! No truncation
    if (ftmax <= 0.00001_R_) then
       ft  = 0.0_R_
       ftc = 0.0_R_
       gtc = g
       iangtc = 1
       rat  = 0.0_R_
       angtc = 0.0_R_
       return
    end if

    ! Find the delta fraction
    nang = size(cosa)
    cumpmin = 1.0_R_ - ftmax
    iangt = gridIdx_bin_I(cump, 1, nang, cumpmin)
    dft = min(0.99999_RD_, 1.0_RD_ - cump(iangt))

    ! P(Q) truncation point
    gtc2 = (g - dft) / (1.0_RD_ - dft)
    gtc2 = max(-1.0_RD_, min(1.0_RD_, gtc2)) ! ideal G value
    sum0 = 0.0_RD_
    sum1 = 0.0_RD_
    do iang = 1, iangt - 1
       sum0 = sum0 + pdf(iang)
       sum1 = sum1 + pdf(iang) * 0.5_R_ * (cosa(iang) + cosa(iang + 1))
    end do
    gtc0 = (g - sum1) / (1.0_RD_ - sum0)
    gtc1 = gtc0
    do iang = iangt, nang - 1
       sum0 = sum0 + pdf(iang)
       sum1 = sum1 + pdf(iang) * 0.5_R_ * (cosa(iang) + cosa(iang + 1))
       gtc1 = (g - sum1) / (1.0_RD_ - sum0)
       if (gtc1 < gtc2) exit
       gtc0 = gtc1
    end do
    iangtc = iang
    ratc = 1.0_RD_
    if (abs(gtc1 - gtc0) > 1.0e-7_RD_) ratc = (gtc1 - gtc2) / (gtc1 - gtc0)
    angtc = ang(iangtc) * ratc + ang(iangtc + 1) * (1.0_RD_ - ratc)
    ftc = 1.0_RD_ - cump(iangtc + 1) - pdf(iangtc) * ratc ! double precision
    ft = dft   ! double --> single
    gtc = gtc2
    rat = 1.0_RD_ - ratc

  end subroutine scatPF_truncVert

end module hparx_radi
