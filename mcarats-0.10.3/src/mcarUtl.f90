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
! Utilities for MCARaTS world
!-
module mcarUtl 

  use globals
  use hparx
  implicit none
  private

  ! Public
  public :: mcarUtl__coneProj
  public :: mcarUtl__getGtab1
  public :: mcarUtl__getGtab2
  public :: mcarUtl__head_read  ! obsolete
  public :: mcarUtl__head_write ! obsolete
  public :: mcarUtl__horiShift
  public :: mcarUtl__binAtm_read
  public :: mcarUtl__binSfc_read
  public :: mcarUtl__binSca_read
  public :: mcarUtl__txtAtm_read
  public :: mcarUtl__txtSfc_read
  public :: mcarUtl__txtSca_read
  public :: mcarUtl__scaleADF
  public :: mcarUtl__imgSolAng_pol_2R
  public :: mcarUtl__imgSolAng_rec_2R

contains

  !+
  ! Get projection parameters for a conical source emission
  !-
  subroutine mcarUtl__coneProj(nquad, q0, dq, cosdc, ppcone, fcone)

    integer,  intent(in) :: nquad ! # of quadrature points
    real(R_), intent(in) :: q0 ! zenith angle (radian) of the center of the cone
    real(R_), intent(in) :: dq ! half cone angle (radian)
    !            If too small, then will be modified to a prescribed value
    real(R_), intent(out) :: cosdc ! 1 - cos(dq), with high accuracy
    real(R_), intent(out) :: fcone ! factor for horizontal-plane flux
    !// (cos(q0) for collimated source with dq=0)
    real(R_), intent(out) :: ppcone ! a projection parameter
    !  fcone = ppcone / pscone, and
    !  horizontal-plane flux = f_source * fcone,
    ! where f_source is source irradiance on the plane normal to the center direction
    ! of the cone.
    !  pscone = 2 * pi * integral(for s=cos(dq),1){abs(s)*ds}
    !  ppcone = integral(for p=0,2*pi){integral(for s=cos(dq),1){abs(p)ds}}
    !  p = cos(q), q is zenith angle for a direction in the cone
    integer,  parameter :: KNQ = 200
    real(R_), save :: xq(KNQ), wq(KNQ)
    integer,  save :: nq = 0
    real(RD_) :: tsum
    real(R_)  :: dq1, abscosq, cosd, cosq, gam
    real(R_)  :: ppa, ppb, pscone, sind, sind2, sinq, sinq2, tmp, x, y
    integer   :: iq, n

    ! Constants
    dq1 = max(REPS_, dq) ! limit minimum because of single precision
    sinq = sin(q0)
    cosq = cos(q0)
    sind = sin(dq1)
    cosd = cos(dq1) ! 1 for collimated source, -1 for isotropic source
    sind2 = sind * sind
    abscosq = abs(cosq)

    ! 1 - cos(Delta)
    if (cosd > 0.0_R_) then
       cosdc = sind2 / (1.0_R_ + cosd)
    else
       cosdc = 1.0_R_ - cosd
    end if

    ! Narrow cone
    if (cosd >= sinq * 0.999999_R_) then
       ppcone = PI_ * sind2 * abscosq
       fcone = abscosq ! cos(q0)

       ! Very wide cone
    else if (cosd <= -sinq * 0.999999_R_) then
       ppcone = PI_ * (2.0_R_ - sind2 * abscosq)
       fcone = (2.0_R_ - sind2 * abscosq) / (2.0_R_ - sind2)

       ! Wide cone
    else
       n = min(KNQ, nquad)
       if (n /= nq) then
          nq = n
          call gaussLegen(nq, xq, wq)
          do iq = 1, nq
             xq(iq) = 0.5_R_ * (1.0_R_ + xq(iq))
          end do              ! xq is set as 0 to 1            
       end if
       gam = cosd / sinq
       sinq2 = sinq * sinq
       tsum = 0.0_R_
       do iq = 1, nq
          x = (1.0_R_ - xq(iq)) * gam + xq(iq)
          tmp = cosq * x / sqrt(1.0_R_ - sinq2 * x*x)
          if (tmp > -1.0_R_ .and. tmp < 1.0_R_) then
             y = x * acos(tmp)
          else if (tmp >= 1.0_R_) then
             y = 0.0_R_
          else
             y = x * PI_
          end if
          tsum = tsum + wq(iq) * y
       end do
       tsum = tsum * (1.0_R_ - gam)
       tmp = PI_ * (1.0_R_ - gam * gam) - 4.0_R_ * tsum
       ppb = sinq * sinq * (PI_ - 2.0_R_ * asin(gam) + cosq * tmp) &
            &        - 2.0_R_ * cosd * sinq * sqrt(1.0_R_ - gam * gam)
       ppa = PI_ * abscosq**3
       ppcone = ppa + ppb
       pscone = PI_ * (1.0_R_ - cosd * abs(cosd))
       fcone = ppcone / pscone
    end if

  end subroutine mcarUtl__coneProj


  !+
  ! Make 1st kind of x or y grid number LUT for numerical diffusion
  !-
  subroutine mcarUtl__getGtab1(ixlut, knd, nxdlo, nxdhi, nx, ix)

    integer,  intent(in) :: ix
    integer,  intent(in) :: nx
    integer,  intent(in) :: knd
    integer,  intent(in) :: nxdhi
    integer,  intent(in) :: nxdlo
    integer,  intent(out) :: ixlut(-knd:knd)
    integer   :: ixd, ixt

    do ixd = -nxdlo, nxdhi
       ixt = ix + ixd
       if (ixt > nx) then
          ixt = ixt - nx
       else if (ixt < 1) then
          ixt = ixt + nx
       end if
       ixlut(ixd) = ixt
    end do

  end subroutine mcarUtl__getGtab1


  !+
  ! Make 2nd kind of x or y grid number LUT for numerical diffusion
  !-
  subroutine mcarUtl__getGtab2(ixlut, knd, nxdlo, nxdhi, x, dxr, fx, nx)

    real(R_), intent(in) :: dxr
    real(R_), intent(in) :: fx
    integer,  intent(in) :: nx
    integer,  intent(in) :: knd
    integer,  intent(in) :: nxdhi
    integer,  intent(in) :: nxdlo
    real(R_), intent(in) :: x
    integer,  intent(out) :: ixlut(-knd:knd)
    integer   :: ixd, ixt
    real(R_)  :: x1

    do ixd = -nxdlo, nxdhi
       x1 = x + dxr * real(ixd, R_)
       if (x1 < 0.0_R_) then
          ixt = int(x1 * fx) + nx
       else
          ixt = int(x1 * fx) + 1
          if (ixt > nx) ixt = ixt - nx
       end if
       ixlut(ixd) = ixt
    end do

  end subroutine mcarUtl__getGtab2


  !+
  ! Horizontal shift
  !-
  subroutine mcarUtl__horiShift(x, ux, path, xmax)

    real(R_), intent(in) :: path
    real(R_), intent(in) :: ux
    real(R_), intent(inout) :: x
    real(R_), intent(in) :: xmax

    x = x + path * ux
    if (x < 0.0_R_ .or. x >= xmax) then
       x = mod(x, xmax)
       if (x < 0.0_R_) x = x + xmax
    end if

  end subroutine mcarUtl__horiShift


  !+
  ! Get & put MCARaTS haeder information
  !-
  subroutine mcarUtl__head_io(iui, iuo, ierr)

    integer,  intent(out) :: ierr
    integer,  intent(in) :: iui
    integer,  intent(in) :: iuo
    integer   :: nx, ny, nz1, nz2, nvar1, nvar2
    real(R_)  :: ptot

    call mcarUtl__head_read(iui, ptot, nx, ny, nz1, nz2, nvar1, nvar2, ierr)
    if (ierr /= 0) return
    call mcarUtl__head_write(iuo, ptot, nx, ny, nz1, nz2, nvar1, nvar2)

  end subroutine mcarUtl__head_io


  !+
  ! Get MCARaTS header information
  !-
  subroutine mcarUtl__head_read(iu, ptot, nx, ny, nz1, nz2, nvar1, nvar2, ierr)

    integer,  intent(out) :: ierr
    integer,  intent(in) :: iu
    integer,  intent(out) :: nx, ny, nz1, nz2, nvar1, nvar2
    real(R_), intent(out) :: ptot
    character(100), save :: signs(1) = (/'%mcarats-0.10'/)

    call skipToSigns(iu, signs, ierr)
    if (ierr /= 0) then
       ierr = 1
       return
    end if

    read (iu, *, err = 1) ptot
    read (iu, *, err = 1) nx, ny, nz1, nz2
    read (iu, *, err = 1) nvar1, nvar2
    ierr = 0
    return
1   call err_read(1, iu, 'mcarUtl__head_read')

  end subroutine mcarUtl__head_read


  !+
  ! Put MCARaTS header information
  !-
  subroutine mcarUtl__head_write(iu, ptot, nx, ny, nz1, nz2, nvar1, nvar2)

    integer,  intent(in) :: iu
    integer,  intent(in) :: nx, ny, nz1, nz2, nvar1, nvar2
    real(R_), intent(in) :: ptot

    write (iu, '(a)',           err=1) '%mcarats-0.10'
    write (iu, '(es12.5e2, a)', err=1) ptot, ' : ptot'
    write (iu, '(4(1x, i5), a)', err=1) nx, ny, nz1, nz2, ' : nx, ny, nz1, nz2'
    write (iu, '(2(1x, i5), a)', err=1) nvar1, nvar2, ' : nvar1, nvar2'
    return
1   call err_write(1, iu, 'mcarUtl__head_write')

  end subroutine mcarUtl__head_write


  !+
  ! Read in 3-D atmosphere model parameters from a text file
  !-
  subroutine mcarUtl__txtAtm_read(iu, nx, ny, nz3, nkd, np3d, mtprof, &
       & tmpa3d, absg3d, apfp3d, extp3d, omgp3d)

    integer,  intent(in) :: iu
    integer,  intent(in) :: mtprof
    integer,  intent(in) :: nkd
    integer,  intent(in) :: np3d
    integer,  intent(in) :: nx, ny
    integer,  intent(in) :: nz3
    real(R_), intent(out) :: tmpa3d(:, :, :), absg3d(:, :, :, :)
    real(R_), intent(out) :: extp3d(:, :, :, :), omgp3d(:, :, :, :)
    real(R_), intent(out) :: apfp3d(:, :, :, :)
    integer   :: ierr, ikd, ip3d, ix, iy, iz

    ! Initialize
    if (nz3 <= 0) return
    call skipToSign(iu, '%mdla3d', ierr)
    call err_read(ierr, iu, 'mcarUtl__txtAtm_read: %mdla3d is not found.')

    ! Temperature
    read (iu, *, err=1, end=9)
    do iz = 1, nz3 + mtprof
       read (iu, *, err=1) ((tmpa3d(ix, iy, iz), ix = 1, nx), iy = 1, ny)
    end do

    ! Gaseous absorption coefficient
    do ikd = 1, nkd
       read (iu, *, err=2, end=9)
       do iz = 1, nz3
          read (iu, *, err=2) ((absg3d(ix, iy, iz, ikd), ix = 1, nx), iy = 1, ny)
       end do
    end do

    ! Cloud & aerosol parameters
    do ip3d = 1, np3d
       read (iu, *, err=3, end=9)
       do iz = 1, nz3 ! extinction coefficients
          read (iu, *, err=3) ((extp3d(ix, iy, iz, ip3d), ix = 1, nx), iy = 1, ny)
       end do
       read (iu, *, err=4, end=9)
       do iz = 1, nz3 ! single scattering albedo
          read (iu, *, err=4) ((omgp3d(ix, iy, iz, ip3d), ix = 1, nx), iy = 1, ny)
       end do
       read (iu, *, err=5, end=9)
       do iz = 1, nz3 ! phase function specification parameter
          read (iu, *, err=5) ((apfp3d(ix, iy, iz, ip3d), ix = 1, nx), iy = 1, ny)
       end do
    end do

    ! Exit
    return
1   call err_read(1, iu, 'mcarUtl__txtAtm_read: tmpa3d at iz ='//num2str_AN(iz))
2   call err_read(1, iu, 'mcarUtl__txtAtm_read: absg3d at iz ='//num2str_AN(iz))
3   call err_read(1, iu, 'mcarUtl__txtAtm_read: extp3d at iz ='//num2str_AN(iz))
4   call err_read(1, iu, 'mcarUtl__txtAtm_read: omgp3d at iz ='//num2str_AN(iz))
5   call err_read(1, iu, 'mcarUtl__txtAtm_read: apfp3d at iz ='//num2str_AN(iz))
9   call err_read(-1, iu, 'mcarUtl__txtAtm_read: Unexpected EOF. Array sizes would be invalid.')

  end subroutine mcarUtl__txtAtm_read


  !+
  ! Read in 3-D atmosphere model parameters from a binary file
  !-
  subroutine mcarUtl__binAtm_read(iu, nx, ny, nz3, nkd, np3d, mtprof, &
       & tmpa3d, absg3d, apfp3d, extp3d, omgp3d, mbswap, idloc)

    integer,  intent(in) :: iu
    integer,  intent(in) :: mtprof ! should be 0 or 1
    integer,  intent(in) :: nkd
    integer,  intent(in) :: np3d
    integer,  intent(in) :: nx, ny
    integer,  intent(in) :: nz3
    real(R_), intent(out) :: tmpa3d(:, :, :), absg3d(:, :, :, :)
    real(R_), intent(out) :: extp3d(:, :, :, :), omgp3d(:, :, :, :)
    real(R_), intent(out) :: apfp3d(:, :, :, :)
    integer, intent(in), optional :: mbswap ! flag for byte swapping (0=no, 1=yes)
    integer, intent(in), optional :: idloc  ! dataset location index
    real(R4_) :: wrk(nx, ny) !(AUTO)
    integer   :: ikd, ip3d, iz, mswap, irec

    ! Initialize
    if (nz3 <= 0) return
    mswap = 0
    if (present(mbswap) .and. mbswap == 1) mswap = 4
    if (present(idloc)) then
       irec = (nz3 + mtprof + nz3 * (nkd + 3*np3d)) * (idloc - 1) ! = (record index to be read) - 1
       call fileRec_set(iu, irec)
    end if

    ! Read in
    do iz = 1, nz3 + mtprof ! temperature
       call bin_read_o2R4(iu, wrk, mswap=mswap)
       tmpa3d(1:nx, 1:ny, iz) = wrk(1:nx, 1:ny)
    end do
    do ikd = 1, nkd ! gaseous absorption coefficient
       do iz = 1, nz3
          call bin_read_o2R4(iu, wrk, mswap=mswap)
          absg3d(1:nx, 1:ny, iz, ikd) = wrk(1:nx, 1:ny)
       end do
    end do
    do ip3d = 1, np3d
       do iz = 1, nz3 ! extinction coefficients
          call bin_read_o2R4(iu, wrk, mswap=mswap)
          extp3d(1:nx, 1:ny, iz, ip3d) = wrk(1:nx, 1:ny)
       end do
       do iz = 1, nz3 ! single scattering albedo
          call bin_read_o2R4(iu, wrk, mswap=mswap)
          omgp3d(1:nx, 1:ny, iz, ip3d) = wrk(1:nx, 1:ny)
       end do
       do iz = 1, nz3 ! phase function specification parameter
          call bin_read_o2R4(iu, wrk, mswap=mswap)
          apfp3d(1:nx, 1:ny, iz, ip3d) = wrk(1:nx, 1:ny)
       end do
    end do

  end subroutine mcarUtl__binAtm_read


  !+
  ! Get phase functions from a text-format file
  !-
  subroutine mcarUtl__txtSca_read(iu, nanci, na, ang, phs)

    integer,  intent(in)  :: iu       ! file unit index
    integer,  intent(in)  :: nanci    ! # of ancillary data blocks
    integer,  intent(out) :: na(:)    ! # of angles
    real(R_), intent(out) :: ang(:,:) ! angles (degrees)
    real(R_), intent(out) :: phs(:,:) ! phase functions
    !// Note: Number of ancillary data lines should be nanci*(np + 1).
    integer :: ia, ipf, nangi, npf, np, ip, naa, ianci, ierr

    ! Header
    nangi = size(ang, 1)
    npf   = size(ang, 2)

    ! Data
    ipf = 1
    loop_trial: do
       call skipToSign(iu, '%mdlphs', ierr)
       call err_read(ierr, iu, 'mcarUtl__txtSca_read: %mdlphs is not found.')
       read (iu, *, err=1, end=1)
       read (iu, *, err=1, end=1) naa
       call check_iI('mcarUtl__txtSca_read: naa', naa, 1, nangi)
       read (iu, *, err=1, end=1) np
       do ianci = 1, nanci
          do ip = 1, np + 1
             read (iu, *, err=1, end=1)
          end do
       end do
       do ip = 1, np
          read (iu, *, err=1, end=1)
          na(ipf) = naa
          do ia = 1, na(ipf)
             read (iu, *, err=1, end=1) ang(ia, ipf), phs(ia, ipf)
          end do
          ipf = ipf + 1
          if (ipf > npf) exit loop_trial
       end do
    end do loop_trial

    return
1   call err_read(1, iu, 'mcarUtl__txtSca_read: ang & phs at ipf='//num2str_AN(ipf))

  end subroutine mcarUtl__txtSca_read


  !+
  ! Read in phase functions from a binary file
  !-
  subroutine mcarUtl__binSca_read(iu, nangi, nanci, nphs, ang, phs, nskip, mbswap)

    integer,  intent(in) :: iu    ! file unit index
    integer,  intent(in) :: nangi ! # of angles
    integer,  intent(in) :: nanci ! # of ancillary data
    integer,  intent(in) :: nphs  ! # of phase functions
    real(R_), intent(out) :: ang(:) ! angles (degrees)
    real(R_), intent(out) :: phs(:,:) ! phase functions
    integer, intent(in), optional :: nskip  ! # of record lines to be skipped
    integer, intent(in), optional :: mbswap ! flag for byte swapping (0=no, 1=yes)
    real(R4_) :: wrk(nangi + nanci) !(AUTO)
    integer :: iphs, mswap, irec

    ! Initialize
    irec = 1
    mswap = 0
    if (present(nskip)) irec = nskip + 1
    if (present(mbswap) .and. mbswap == 1) mswap = 4

    ! Read in
    call bin_read_o1R4(iu, wrk, irec=irec, mswap=mswap)
    ang(1:nangi) = wrk(1:nangi)
    do iphs = 1, nphs
       irec = irec + 1
       call bin_read_o1R4(iu, wrk, irec=irec, mswap=mswap)
       phs(1:nangi, iphs) = wrk(1:nangi)
    end do

  end subroutine mcarUtl__binSca_read


  !+
  ! Read in surface model parameters from a text file
  !-
  subroutine mcarUtl__txtSfc_read(iu, nxb, nyb, npsfc, tmps2d, psfc2d, jsfc2d)

    integer,  intent(in) :: iu
    integer,  intent(in) :: nxb
    integer,  intent(in) :: nyb
    integer,  intent(in) :: npsfc
    integer,  intent(out) :: jsfc2d(:, :)
    real(R_), intent(out) :: tmps2d(:, :), psfc2d(:, :, :)
    integer   :: ierr, ipsfc, ixb, iyb

    ! Signature
    call skipToSign(iu, '%mdlsfc', ierr)
    call err_read(ierr, iu, 'mcarUtl__txtSfc_read: %mdlsfc is not found.')

    ! Surface model
    ipsfc = 0
    read (iu, *, err=1, end=2)
    read (iu, *, err=1, end=2) ((tmps2d(ixb, iyb), ixb = 1, nxb), iyb = 1, nyb)
    read (iu, *, err=1, end=2)
    read (iu, *, err=1, end=2) ((jsfc2d(ixb, iyb), ixb = 1, nxb), iyb = 1, nyb)
    read (iu, *, err=1, end=2)
    do ipsfc = 1, npsfc
       read (iu, *, err=1, end=2) ((psfc2d(ixb, iyb, ipsfc), ixb = 1, nxb), iyb = 1, nyb)
    end do

    return
1   call err_read( 1, iu, 'mcarUtl__txtSfc_read: ipsfc='//num2str_AN(ipsfc))
2   call err_read(-1, iu, 'mcarUtl__txtSfc_read: ipsfc='//num2str_AN(ipsfc))

  end subroutine mcarUtl__txtSfc_read


  !+
  ! Read in surface model parameters from a binary file
  !-
  subroutine mcarUtl__binSfc_read(iu, nxb, nyb, npsfc, tmps2d, psfc2d, jsfc2d, mbswap, idloc)

    integer,  intent(in) :: iu
    integer,  intent(in) :: nxb
    integer,  intent(in) :: nyb
    integer,  intent(in) :: npsfc
    integer,  intent(out) :: jsfc2d(:, :)
    real(R_), intent(out) :: tmps2d(:, :), psfc2d(:, :, :)
    integer, intent(in), optional :: mbswap ! flag for byte swapping (0=no, 1=yes)
    integer, intent(in), optional :: idloc  ! dataset location index
    real(R4_) :: wrk(nxb, nyb) !(AUTO)
    integer   :: ipsfc, mswap, irec

    ! Initialize
    mswap = 0
    if (present(mbswap) .and. mbswap == 1) mswap = 4
    if (present(idloc)) then
       irec = (2 + npsfc) * (idloc - 1) ! = (record index to be read) - 1
       call fileRec_set(iu, irec)
    end if

    ! Read in
    call bin_read_o2R4(iu, wrk, mswap=mswap)
    tmps2d(1:nxb, 1:nyb) = wrk(1:nxb, 1:nyb)
    call bin_read_o2R4(iu, wrk, mswap=mswap)
    jsfc2d(1:nxb, 1:nyb) = nint(wrk(1:nxb, 1:nyb))
    do ipsfc = 1, npsfc
       call bin_read_o2R4(iu, wrk, mswap=mswap)
       psfc2d(1:nxb, 1:nyb, ipsfc) = wrk(1:nxb, 1:nyb)
    end do

  end subroutine mcarUtl__binSfc_read


  !+
  ! Scale ADF (angular distribution function)
  !-
  subroutine mcarUtl__scaleADF(adf, adfmin, taulim)

    real(R_), intent(inout) :: adf
    real(R_), intent(in) :: adfmin
    real(R_), intent(out) :: taulim

    ! ADF is small
    if (adf <= adfmin) then
       if (mseq_rand_R() * adfmin >= adf) then
          adf = 0.0_R_
       else
          adf = adfmin
       end if
       taulim = 0.0_R_

       ! ADF is large
    else
       taulim = -log(adfmin / adf)
    end if

  end subroutine mcarUtl__scaleADF


  !+
  ! Integrals over solid angles of every pixels for a (distorted) rectangular image
  !-
  function mcarUtl__imgSolAng_rec_2R(nxr, nyr, mfunc, hc, vc) result(ds)

    integer,  intent(in) :: nxr, nyr ! # of image pixels (horizontal & vertical)
    integer,  intent(in) :: mfunc    ! 0 for w=1, 1 for w=cosQ
    real(R_), intent(in) :: hc, vc   ! half FOV angles (horizontal & vertical)
    real(R_) :: ds(nxr, nyr) ! result, integrals over solid angles
    !// If mfunc=0, ds = integral(for a solid angle){1} = solid angle
    !   If mfunc=1, ds = integral(for a solid angle){cosQ}, Q = angle from the image center direction
    real(RD_) :: sum1
    real(R_)  :: auc, avc, dh, dhmax, dv, dvmax, h, hmin, r, v, vmin, d, f
    integer   :: ih, iv, ixr, ixr2, iyr, iyr2, nh, nv, nxr2, nyr2

    ! Initialize
    d = 0.1_R_ * DTOR_
    auc = 2.0_R_ * hc
    avc = 2.0_R_ * vc
    dhmax = auc / real(nxr, R_)
    dvmax = avc / real(nyr, R_)
    nh = min(10, int(dhmax / d) + 1)
    nv = min(10, int(dvmax / d) + 1)
    dh = dhmax / real(nh, R_)
    dv = dvmax / real(nv, R_)
    nxr2 = (nxr + 1) / 2
    nyr2 = (nyr + 1) / 2

    ! Loop for a quarter iamge
    do iyr2 = 1, nyr2
       vmin = ((iyr2 - 1.0_R_) / real(nyr, R_) - 0.5_R_) * avc
       do ixr2 = 1, nxr2
          hmin = ((ixr2 - 1.0_R_) / real(nxr, R_) - 0.5_R_) * auc
          sum1 = 0.0_RD_
          do ih = 1, nh
             h = hmin + dh * (ih - 0.5_R_)
             do iv = 1, nv
                v = vmin + dv * (iv - 0.5_R_)
                r = max(RSML_, sqrt(h * h + v * v))
                f = sin(r) / r
                if (mfunc == 1) f = f * max(REPS_, cos(r))
                sum1 = sum1 + f
             end do
          end do
          ds(ixr2, iyr2) = sum1 * dh * dv ! solid angle [steradian]
       end do
    end do

    ! Full image
    do iyr = 1, nyr
       iyr2 = iyr
       if (iyr > nyr2) iyr2 = nyr - iyr + 1
       do ixr = 1, nxr
          ixr2 = ixr
          if (ixr > nxr2) ixr2 = nxr - ixr + 1
          ds(ixr, iyr) = ds(ixr2, iyr2)
       end do
    end do

  end function mcarUtl__imgSolAng_rec_2R


  !+
  ! Integrals over solid angles of every pixels for a polar image
  !-
  function mcarUtl__imgSolAng_pol_2R(nxr, nyr, mfunc, hc, vc) result(ds)

    integer,  intent(in) :: nxr, nyr ! # of image pixels (theta & phi)
    integer,  intent(in) :: mfunc    ! 0 for w=1, 1 for w=cosQ
    real(R_), intent(in) :: hc, vc   ! half FOV angles (the max polar and azimuth angles)
    real(R_) :: ds(nxr, nyr) ! result, solid angles
    !// If mfunc=0, ds = integral(for a solid angle){1} = solid angle
    !   If mfunc=1, ds = integral(for a solid angle){cosQ}, Q = angle from the image center direction
    real(R_) :: dx, dy, xmin, xmax, q
    integer  :: ixr

    ! Initialize
    dx = hc / real(nxr, R_)
    dy = vc / real(nyr, R_)

    ! Loop for pixels
    do ixr = 1, nxr
       xmax = dx * ixr
       xmin = xmax - dx
       if (xmin > 25.0_R_ * DTOR_) then
          if (mfunc == 0) then
             q = dy * abs(cos(xmin) - cos(xmax))
          else
             q = dy * abs(cos(xmin)**2 - cos(xmax)**2)
          end if
       else
          if (mfunc == 0) then
             q = dy * abs(sin(xmin)**2 / (1.0_R_ + cos(xmin)) - sin(xmax)**2 / (1.0_R_ + cos(xmax)))
          else
             q = dy * abs(sin(xmin)**2 - sin(xmax)**2)
          end if
       end if !// difference of (1-cos), instead of one of cos
       ds(ixr, :) = q
    end do

  end function mcarUtl__imgSolAng_pol_2R

end module mcarUtl
