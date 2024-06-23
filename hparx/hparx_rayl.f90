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
! Library of Rayleigh scattering utilities
!-
module hparx_rayl

  use globals, only : R_, PI_, AVOGAD_, GASCON_, FRAC13_, REPS_
  use hparx_rand, only : mseq_rand_R
  implicit none
  private

  ! Public
  public :: rayl_kingFac_R
  public :: rayl_randScat
  public :: rayl_randScat2
  public :: rayl_refAnom_R
  public :: rayl_refAnom_acc_R
  public :: rayl_scatPF_R
  public :: rayl_scatPF_1R
  public :: rayl_xsec_R

  ! Depolarization coefficients
  real(R_), parameter :: Rayl_DEPOL = 0.0279_R_ ! depolarization factor, after Young (1980 AO)
  real(R_), parameter :: Rayl_PCOA = 3.0_R_ * (1.0_R_ + Rayl_DEPOL) / (4.0_R_ + 2.0_R_ * Rayl_DEPOL)
  real(R_), parameter :: Rayl_PCOB = 3.0_R_ * (1.0_R_ - Rayl_DEPOL) / (4.0_R_ + 2.0_R_ * Rayl_DEPOL)

  ! Normalization factors for air refraction
  real(R_), parameter :: Rayl_PSTD = 101325.0_R_ ! standard pressure (Pa)
  real(R_), parameter :: Rayl_FACS = 1.0e-8_R_ * GASCON_ * 288.15_R_ / Rayl_PSTD ! a factor (10^-8 m^3/mol)
  real(R_), parameter :: Rayl_FACL = 1.0e-8_R_ * GASCON_ * 273.15_R_ / Rayl_PSTD ! a factor (10^-8 m^3/mol)
  real(R_), parameter :: Rayl_FACB = 1.0e-8_R_ * GASCON_ * 292.75_R_ / Rayl_PSTD ! a factor (10^-8 m^3/mol)

  !+ Air refraction formulas
  !  mform = 1 : Ed53, Edlen (1953, JOSA)
  !  mform = 2 : Ed66, Edlen (1966, Metrologia)
  !  mform = 3 : PR72, Peck and Reeder (1972, JOSA)
  !  mform = 4 : Ba84, Bates (1984, Planet. Space Sci.)
  !  mform = 5 : BD94, Birch & Downs (1993, 1994, Metrologia)
  !  Ed53 : T = 288.15 K, p = 101325 Pa, CO2 = 300 ppmV
  !  Ed66 : T = 288.15 K, p = 101325 Pa, CO2 = 300 ppmV, WL = 0.3889-1.529 micron
  !  PR72 : T = 288.15 K, p = 101325 Pa, CO2 = 300 (330?) ppmV, WL = 0.23-1.69 micron
  !  Ba84 : T = 273.15 K, p = 101325 Pa, CO2 = ??? ppmV
  !  BD94 : T = 292.75 K, p = 101325 Pa, CO2 = 450 ppmV, wavelength dependence is from Ed66
  !-

  ! Coefficients for refraction of standard air
  real(R_), parameter :: Rayl_E5_C0 = 6.4328e+3_R_,  Rayl_E5_C1 = 2.94981e+6_R_,  Rayl_E5_C2 = 2.554e+4_R_
  real(R_), parameter :: Rayl_ED_C0 = 8.34213e+3_R_, Rayl_ED_C1 = 2.40603e+6_R_,  Rayl_ED_C2 = 1.5997e+4_R_
  real(R_), parameter :: Rayl_PR_C0 = 8.06051e+3_R_, Rayl_PR_C1 = 2.48099e+6_R_,  Rayl_PR_C2 = 1.74557e+4_R_
  real(R_), parameter :: Rayl_BA_C0 = 7.041e+3_R_,   Rayl_BA_C1 = 3.159e+6_R_,    Rayl_BA_C2 = 8.4127e+4_R_
  real(R_), parameter :: Rayl_BD_C0 = 8.34254e+3_R_, Rayl_BD_C1 = 2.406147e+6_R_, Rayl_BD_C2 = 1.5998e+4_R_
  real(R_), parameter :: Rayl_E5_C3 = 146.0_R_,   Rayl_E5_C4 = 41.0_R_
  real(R_), parameter :: Rayl_ED_C3 = 130.0_R_,   Rayl_ED_C4 = 38.9_R_
  real(R_), parameter :: Rayl_PR_C3 = 132.274_R_, Rayl_PR_C4 = 39.32957_R_
  real(R_), parameter :: Rayl_BA_C3 = 157.39_R_,  Rayl_BA_C4 = 50.429_R_
  real(R_), parameter :: Rayl_BD_C3 = 130.0_R_,   Rayl_BD_C4 = 38.9_R_

  ! Coefficients for H2O correction
  real(R_), parameter :: Rayl_ED_H2O_C0 = 4.29185e-2_R_, Rayl_ED_H2O_C1 = 3.42778e-4_R_
  real(R_), parameter :: Rayl_BD_H2O_C0 = 3.7345e-2_R_,  Rayl_BD_H2O_C1 = 4.01e-4_R_

  ! Coefficients for CO2 correction
  real(R_), parameter :: Rayl_ED_CO2_C0 = 0.54_R_,  Rayl_ED_CO2_C1 = 0.0003_R_
  real(R_), parameter :: Rayl_BD_CO2_C0 = 0.536_R_, Rayl_BD_CO2_C1 = 0.00045_R_

  ! Coefficients for dependence on p and T
  real(R_), parameter :: Rayl_ED_P0 = 96095.43_R_
  real(R_), parameter :: Rayl_ED_T0 = 0.6128e-8_R_, Rayl_ED_T1 = 9.97582e-11_R_, Rayl_ED_T2 = 0.003661_R_
  real(R_), parameter :: Rayl_BD_T0 = 0.601e-8_R_,  Rayl_BD_T1 = 9.72e-11_R_,    Rayl_BD_T2 = 0.003661_R_

contains

  !+
  ! King factor of dry air, according to Bodhaine et al. (1999)
  !-
  function rayl_kingFac_R(wav, xco2) result(res)

    real(R_), intent(in), optional :: wav  ! wavelength (micron)
    real(R_), intent(in), optional :: xco2 ! volume mixing ratio of CO2 (360e-6_R_ for 360 ppmv)
    real(R_) :: res, aww, fn2, fo2, p
    ! Bodhaine et al. (1999, JAOT), based on Bates (1984)
    real(R_), parameter :: C0_N2 = 1.034_R_, C1_N2 = 3.17e-4_R_ ! coefficients for N2
    real(R_), parameter :: C0_O2 = 1.096_R_, C1_O2 = 1.385e-3_R_, C2_O2 = 1.448e-4_R_ ! coefficients for O2
    real(R_), parameter :: F_AR  = 1.0_R_  ! King factor for Ar
    real(R_), parameter :: F_CO2 = 1.15_R_ ! King factor for CO2
    real(R_), parameter :: P_N2 = 78.084e-2_R_, P_O2 = 20.946e-2_R_, P_AR = 0.934e-2_R_, P_CO2 = 360e-6_R_
    real(R_), parameter :: P_SUM = P_N2 + P_O2 + P_AR

    ! Bodhaine et al. (1999) with wavelength dependence
    if (present(wav)) then
       aww = 1.0_R_ / max(0.1_R_, min(3.0_R_, wav))**2
       fn2 = C0_N2 + C1_N2 * aww
       fo2 = C0_O2 + C1_O2 * aww + C2_O2 * aww**2
       p = P_CO2 ! default CO2
       if (present(xco2)) p = xco2
       res = (P_N2 * fn2 + P_O2 * fo2 + P_AR * F_AR + p * F_CO2) / (P_SUM + p)

       ! Young (1980, AO)
    else
       res = (6.0_R_ + 3.0_R_ * Rayl_DEPOL) / (6.0_R_ - 7.0_R_ * Rayl_DEPOL)
    end if

  end function rayl_kingFac_R


  !+
  ! Determines a random scattering angle of Rayleigh scattering event
  !  using the rejection method
  !-
  subroutine rayl_randScat(sinq, cosq)

    real(R_), intent(out) :: sinq, cosq ! sine & cosine of the scattering angles
    integer  :: iloop
    real(R_) :: rn1, rn2, gam, c
    integer, parameter :: NLOOP_R = 20
    integer, parameter :: NLOOP_N = 2

    ! Loop due to the rejection method, without depolarization
    do iloop = 1, NLOOP_R
       rn1 = mseq_rand_R()
       rn2 = mseq_rand_R()
       gam = rn1 * (1.0_R_ - rn1)
       if (rn2 >= 2.0_R_ * gam) exit
    end do
    cosq = 1.0_R_ - 2.0_R_ * rn1

    ! Newtonian iteration, taking into account the depolarization
    c = 3.0_R_ * (1.0_R_ + cosq - Rayl_PCOA) - Rayl_PCOB
    do iloop = 1, NLOOP_N
       cosq = (2.0_R_ * Rayl_PCOB * cosq**3 + c) / (Rayl_PCOB * cosq**2 + Rayl_PCOA) * FRAC13_
    end do
    if (cosq < -1.0_R_) cosq = -1.0_R_
    if (cosq >  1.0_R_) cosq =  1.0_R_
    sinq = sqrt(1.0_R_ - cosq**2)

  end subroutine rayl_randScat


  !+
  ! Determines a random scattering angle of Rayleigh scattering event
  !  using a deterministic, Newtonian method
  !-
  subroutine rayl_randScat2(xi, sinq, cosq)

    real(R_), intent(in)  :: xi         ! a uniform random number
    real(R_), intent(out) :: sinq, cosq ! sine & cosine of the scattering angles
    integer,  parameter :: NLOOP = 6    ! # of revisions for better accuracy
    !// For xi = 0, cosq = +1 (forward scattering)
    !   For xi = 1, cosq = -1 (backscattering)
    integer  :: iloop
    real(R_) :: c, cosq1, delta

    ! Initial estimate for isotropic scattering
    cosq = 1.0_R_ - 2.0_R_ * xi

    ! Newtonian iteration
    c = 3.0_R_ * (1.0_R_ + cosq - Rayl_PCOA) - Rayl_PCOB
    do iloop = 1, NLOOP
       cosq1 = (2.0_R_ * Rayl_PCOB * cosq**3 + c) / (Rayl_PCOB * cosq**2 + Rayl_PCOA) * FRAC13_       
       delta = abs(cosq1 - cosq)
       cosq = cosq1
       if (delta < REPS_) exit
    end do
    if (cosq < -1.0_R_) cosq = -1.0_R_
    if (cosq >  1.0_R_) cosq =  1.0_R_
    sinq = sqrt(1.0_R_ - cosq**2)

  end subroutine rayl_randScat2


  !+
  ! Refractive index anomaly (from 1) of air per unit number density (m^3/mol) using simple formulae
  !  calculated taking into account CO2 and H2O mixing ratios
  ! Note: This function returns n-1 not taking p and T into account, which implicitly assumes an
  !  air near the standard state.
  !-
  function rayl_refAnom_R(w, mform, xco2, xh2o) result (res)

    real(R_), intent(in) :: w      ! wavelength (micometer)
    !// Note used formula is accurate for a wavelength range of 0.25-1.7 micron.
    !   The uncertainty rapidly increases for shorter wavelendths than 0.25 micron.
    integer,  intent(in), optional :: mform ! flag (1/2/3/4/5) for an original paper
    !// mform = 1 : Ed53, Edlen (1953, JOSA)
    !// mform = 2 : Ed66, Edlen (1966, Metrologia)
    !// mform = 3 : PR72, Peck and Reeder (1972, JOSA)
    !// mform = 4 : Ba84, Bates (1984, Planet. Space Sci.)
    !// mform = 5 : BD94, Birch & Downs (1993, 1994, Metrologia)
    real(R_), intent(in), optional :: xco2  ! volume mixing ratio of CO2 ( 380e-6_R_ for  380 ppmV)
    real(R_), intent(in), optional :: xh2o  ! volume mixing ratio of H2O (1500e-6_R_ for 1500 ppmV)
    !// By default, mform = 5, xco2 = 380 ppmV, xh2o = 0 ppmV
    real(R_) :: res   ! refractive index anomaly (from 1) per unit number density of air (m^3/mol)
    real(R_) :: aww, xxco2
    integer :: mf
    real(R_), save :: fnorm(5) = (/ Rayl_FACS, Rayl_FACS, Rayl_FACS, Rayl_FACL, Rayl_FACS /) ! normalization factors

    ! Initialize
    mf = 5
    if (present(mform)) mf = max(1, min(5, mform))
    xxco2 = 380.0_R_
    if (present(xco2)) xxco2 = xco2
    aww = 1.0_R_ / w**2

    ! Refractive index anomaly
    if (mf == 1) then ! Ed53
       res = Rayl_E5_C0 + Rayl_E5_C1 / (Rayl_E5_C3 - aww) + Rayl_E5_C2 / (Rayl_E5_C4 - aww)
    else if (mf == 2) then ! Ed66
       res = Rayl_ED_C0 + Rayl_ED_C1 / (Rayl_ED_C3 - aww) + Rayl_ED_C2 / (Rayl_ED_C4 - aww)
    else if (mf == 3) then ! PR72
       res = Rayl_PR_C0 + Rayl_PR_C1 / (Rayl_PR_C3 - aww) + Rayl_PR_C2 / (Rayl_PR_C4 - aww)
    else if (mf == 4) then ! Ba84
       res = Rayl_BA_C0 + Rayl_BA_C1 / (Rayl_BA_C3 - aww) + Rayl_BA_C2 / (Rayl_BA_C4 - aww)
    else ! BD94
       res = Rayl_BD_C0 + Rayl_BD_C1 / (Rayl_BD_C3 - aww) + Rayl_BD_C2 / (Rayl_BD_C4 - aww)
    end if

    ! Corrections
    if (mf <= 4) then ! Ed53, Ed66, PR72, or Ba84
       if (present(xh2o)) res = res - xh2o * Rayl_PSTD * (Rayl_ED_H2O_C0 - Rayl_ED_H2O_C1 * aww)
       res = res * fnorm(mf) ! normalize by air number density
       res = res * (1.0_R_ + Rayl_ED_CO2_C0 * (xxco2 - Rayl_ED_CO2_C1))
    else ! BD94
       if (present(xh2o)) res = res - xh2o * Rayl_PSTD * (Rayl_BD_H2O_C0 - Rayl_BD_H2O_C1 * aww)
       res = res * fnorm(mf) ! normalize by air number density
       res = res * (1.0_R_ + Rayl_BD_CO2_C0 * (xxco2 - Rayl_BD_CO2_C1))
    end if

  end function rayl_refAnom_R


  !+
  ! Accurate refractive index anomaly (from 1) of air at a specific pressure and temperature,
  !  calculated taking into account CO2 and H2O mixing ratios
  !-
  function rayl_refAnom_acc_R(w, p, t, mform, xco2, xh2o) result (res)

    real(R_), intent(in) :: w ! wavelength (micometer)
    !// Note used formula is accurate for a wavelength range of 0.25-1.7 micron.
    !   The uncertainty rapidly increases for shorter wavelendths than 0.25 micron.
    real(R_), intent(in) :: p ! pressure (Pa)
    real(R_), intent(in) :: t ! temperature (K)
    integer,  intent(in), optional :: mform ! flag (1/2/3/4/5) for an original paper
    !// mform = 1 : Ed53, Edlen (1953, JOSA)
    !// mform = 2 : Ed66, Edlen (1966, Metrologia)
    !// mform = 3 : PR72, Peck and Reeder (1972, JOSA)
    !// mform = 4 : Ba84, Bates (1984, Planet. Space Sci.)
    !// mform = 5 : BD94, Birch & Downs (1993, 1994, Metrologia)
    real(R_), intent(in), optional :: xco2  ! volume mixing ratio of CO2 ( 380e-6_R_ for  380 ppmV)
    real(R_), intent(in), optional :: xh2o  ! volume mixing ratio of H2O (1500e-6_R_ for 1500 ppmV)
    !// By default, mform = 5, xco2 = 380 ppmV, xh2o = 0 ppmV
    real(R_) :: res   ! refractive index anomaly (from 1)
    real(R_) :: aww, xxco2, tt
    integer :: mf

    ! Initialize
    mf = 5
    if (present(mform)) mf = max(1, min(5, mform))
    xxco2 = 380.0_R_
    if (present(xco2)) xxco2 = xco2
    aww = 1.0_R_ / w**2
    tt = t - 273.15_R_

    ! Refractive index anomaly
    if (mf == 1) then ! Ed53
       res = Rayl_E5_C0 + Rayl_E5_C1 / (Rayl_E5_C3 - aww) + Rayl_E5_C2 / (Rayl_E5_C4 - aww)
    else if (mf == 2) then ! Ed66
       res = Rayl_ED_C0 + Rayl_ED_C1 / (Rayl_ED_C3 - aww) + Rayl_ED_C2 / (Rayl_ED_C4 - aww)
    else if (mf == 3) then ! PR72
       res = Rayl_PR_C0 + Rayl_PR_C1 / (Rayl_PR_C3 - aww) + Rayl_PR_C2 / (Rayl_PR_C4 - aww)
    else if (mf == 4) then ! Ba84
       res = Rayl_BA_C0 + Rayl_BA_C1 / (Rayl_BA_C3 - aww) + Rayl_BA_C2 / (Rayl_BA_C4 - aww)
       res = res * (273.15_R_ / 288.15_R_) ! correction due to standard air temperature difference
    else ! BD94
       res = Rayl_BD_C0 + Rayl_BD_C1 / (Rayl_BD_C3 - aww) + Rayl_BD_C2 / (Rayl_BD_C4 - aww)
    end if

    ! Corrections
    if (mf <= 4) then ! Ed53, Ed66, PR72, or Ba84
       res = res * (p / Rayl_ED_P0) * (1.0_R_ + (Rayl_ED_T0 - Rayl_ED_T1 * tt)) / (1.0_R_ + Rayl_ED_T2 * tt)
       if (present(xh2o)) res = res - xh2o * p * (288.15_R_ / t) * (Rayl_ED_H2O_C0 - Rayl_ED_H2O_C1 * aww)
       res = res * (1.0_R_ + Rayl_ED_CO2_C0 * (xxco2 - Rayl_ED_CO2_C1))
    else ! BD94
       res = res * (p / Rayl_ED_P0) * (1.0_R_ + (Rayl_BD_T0 - Rayl_BD_T1 * tt)) / (1.0_R_ + Rayl_BD_T2 * tt)
       if (present(xh2o)) res = res - xh2o * p * (292.75_R_ / t) * (Rayl_BD_H2O_C0 - Rayl_BD_H2O_C1 * aww)
       res = res * (1.0_R_ + Rayl_BD_CO2_C0 * (xxco2 - Rayl_BD_CO2_C1))
    end if
    res = res * 1.0e-8_R_
    !// Note that we use a slight modification of the water vapor term, introducing temperature dependence.
    !   This modification can significantly improve the accuracy of the water vapor term, so as to account
    !   for number density change by temperature.

  end function rayl_refAnom_acc_R


  !+
  ! Rayleigh scattering phase function
  !-
  function rayl_scatPF_R(cosa) result(phs) 

    real(R_), intent(in) :: cosa ! cosine of angle
    real(R_) :: phs ! phase function
    phs = Rayl_PCOA + Rayl_PCOB * cosa**2

  end function rayl_scatPF_R


  !+
  ! Rayleigh scattering phase function
  !-
  function rayl_scatPF_1R(cosa) result(phs) 

    real(R_), intent(in) :: cosa(:) ! cosines of angles
    real(R_) :: phs(size(cosa)) ! phase functions
    phs(:) = Rayl_PCOA + Rayl_PCOB * cosa(:)**2

  end function rayl_scatPF_1R


  !+
  ! Rayleigh scattering cross section of air per mole (m^2/mol)
  !-
  function rayl_xsec_R(wl, mform, xco2, xh2o) result(res)

    real(R_), intent(in) :: wl ! wavelength (micron)
    !// Note used formula is accurate for a wavelength range of 0.25-1.6 micron.
    integer,  intent(in), optional :: mform ! flag (-1/0/1/2/3/4/5) for an original paper
    !// mform = -1 (default): Bodhaine et al. (1999), 0: Callan (1984)
    !// mform >= 1: See rayl_refAnom_R().
    real(R_), intent(in), optional :: xco2  ! volume mixing ratio of CO2 ( 360e-6_R_ for  360 ppmV)
    real(R_), intent(in), optional :: xh2o  ! volume mixing ratio of H2O (1500e-6_R_ for 1500 ppmV)
    real(R_) :: res ! Rayleigh scattering cross section of air per mole (m^2/mol)
    real(R_) :: ww, aww, drefr
    integer  :: mf
    ! Bodhaine et al. (1999), Eq. (29), 360 ppmv CO2
    !  Fitting error is 0.01% over the 0.25-0.85 micron range, and better than 0.05% over the 0.25-1.0 micron
    real(R_), parameter :: BOD_A0 = 1.0455996_R_, BOD_A1 = -341.29061_R_, BOD_A2 = -0.9023085_R_
    real(R_), parameter :: BOD_B1 = 2.7059889e-3_R_, BOD_B2 = -85.968563_R_
    ! Callan (1984), ??? ppmv CO2
    real(R_), parameter :: CAL_A0 = 3.9729066_R_,    CAL_A1 = 4.6547659e-2_R_
    real(R_), parameter :: CAL_A2 = 4.5055995e-4_R_, CAL_A3 = 2.3229848e-5_R_
    ! Edlen (1966) for corrections due to CO2
    real(R_), parameter :: ED_CO2_C1 = 0.00036_R_
    !// Note Edlen's data is for 300 ppmV but Bodhaine et al. for 360 ppmV. We assume an influence due
    !   to the difference is small and the correction coefficient is the same for 360 ppmV as for 300 ppmV.
    ! Other parameters
    real(R_), parameter :: FACC = AVOGAD_ * 1.0e-32_R_ ! (1e-32 /mol)
    real(R_), parameter :: FACR = 24.0_R_ * PI_**3

    ! Initialize
    mf = -1 ! default
    if (present(mform)) mf = mform

    ! Empirical formulas
    if (mf <= 0) then

       ! Standard state
       ww = max(0.25_R_, min(1.0_R_, wl))**2
       aww = 1.0_R_ / ww
       if (mf == -1) then ! Bodhaine et al. (1999)
          res = (BOD_A0 + BOD_A1 * aww + BOD_A2 * ww) / (1.0_R_ + BOD_B1 * aww + BOD_B2 * ww) * FACC
       else ! Callan (1984)
          res = (CAL_A0 + (CAL_A1 + (CAL_A2 + CAL_A3 * aww) * aww) * aww) * aww**2 * FACC
       end if

       ! Extrapolation for too short/long wavelength
       if (wl > 1.0_R_) then
          res = res * (1.0_R_ / wl)**4
       else if (wl < 0.25_R_) then
          res = res * (0.25_R_ / wl)**4
       end if

       ! Corrections due to CO2 & H2O (Edlen, 1953, 1966)
       if (present(xco2)) then
          ww = Rayl_ED_CO2_C0 * (xco2 - ED_CO2_C1) 
          res = res * (ww * (ww + 2.0_R_) + 1.0_R_) ! = (1 + ww)**2 because sigma ~ (n - 1)**2
       end if
       if (present(xh2o)) then
          ww = -xh2o * Rayl_PSTD * (Rayl_ED_H2O_C0 - Rayl_ED_H2O_C1 * aww) / &
               (Rayl_ED_C0 + Rayl_ED_C1 / (Rayl_ED_C3 - aww) + Rayl_ED_C2 / (Rayl_ED_C4 - aww))
          res = res * (ww * (ww + 2.0_R_) + 1.0_R_)
       end if

       ! Formulas according to refractive index
    else
       ww = max(0.25_R_, wl)
       drefr = rayl_refAnom_R(ww, mform, xco2, xh2o)
       res = drefr * (2.0_R_ + drefr) / (3.0_R_ + drefr * (2.0_R_ + drefr)) ! = (n^2 - 1) / (n^2 + 2)
       res = FACR * (1.0e+24_R_ / AVOGAD_) * res**2 / wl**4 * rayl_kingFac_R(ww, xco2)
    end if

  end function rayl_xsec_R

end module hparx_rayl
