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
! Module for flux samplers
!-
module mcarFlx 

  use globals
#if UseMPI == 1
  use globals_MPI
#endif
  use hparx
  use mcarUtl
  use mcarAtm
  use mcarSfc
  use mcarPho  ! Pho_* is readable & writable
  implicit none
  private

  ! Public
  public :: mcarFlx__tune_init
  public :: mcarFlx__user_init
  public :: mcarFlx__init
  public :: mcarFlx__final
  public :: mcarFlx__qry_flags
  public :: mcarFlx__qry_dims
  public :: mcarFlx__zero
  public :: mcarFlx__sampSrc
  public :: mcarFlx__sampFlux
  public :: mcarFlx__sampFluxSfc
  public :: mcarFlx__sampHeat
  public :: mcarFlx__sampHeatD
#if UseMPI == 1
  public :: mcarFlx__reduce
#endif
  public :: mcarFlx__normal
  public :: mcarFlx__writeBin
  public :: mcarFlx__writeTxt
  public :: mcarFlx__writeCtl

  ! (Private) User variables in namelist
  integer,  save :: Flx_mflx = 1 ! flag for flux density calculations (0=off, 1=on)
  integer,  save :: Flx_mhrt = 1 ! flag for heating rate calculations (0=off, 1=on)
  real(R_), save :: Flx_diff0 = 1.0_R_ ! numerical diffusion parameter
  real(R_), save :: Flx_diff1 = 0.1_R_ ! numerical diffusion parameter
  namelist /mcarFlx_nml_init/Flx_mflx, Flx_mhrt, Flx_diff0, Flx_diff1

  ! (Private) Flux & heating rate integrator
  integer,   save :: Flx_nflx = 3 ! should be = 3
  integer,   save :: Flx_nx
  integer,   save :: Flx_ny
  integer,   save :: Flx_nz
  integer,   save :: Flx_nzf
  integer,   save :: Flx_nzh
  real(R_),  save :: Flx_anxf, Flx_anyf
  integer,   save :: Flx_ndset = 0 ! # of output datasets
  real(RD_), save :: Flx_psrc ! # of source photons used
  real(RD_), save :: Flx_esrc ! sampler for average total energy
  integer,   save, allocatable :: Flx_iizf(:)       !(0:nz)
  real(RD_), save, allocatable :: Flx_pflx(:,:,:,:) !(nx,ny,nzf,nflx)
  real(RD_), save, allocatable :: Flx_phrt(:,:,:)   !(nx,ny,nzh)
  integer,   parameter :: Flx_KND = 30 ! max # of pixels for diffusion

contains

  !+
  ! Tune a preset for technical parameters
  !-
  subroutine mcarFlx__tune_init(moptim) 

    integer, intent(in) :: moptim ! optimization flag [-2,3]

    if (moptim <= -2) then ! no preset
       return
    else if (moptim == -1) then ! no optimization
       call mcarFlx__set_tech(0.0_R_, 0.0_R_)
    else if (moptim == 0) then ! unbiased optimization
       call mcarFlx__set_tech(0.0_R_, 0.0_R_)
    else if (moptim == 1) then ! conservative
       call mcarFlx__set_tech(1.0_R_, 0.01_R_)
    else if (moptim == 2) then ! standard
       call mcarFlx__set_tech(2.0_R_, 0.1_R_)
    else ! quick-and-dirty
       call mcarFlx__set_tech(4.0_R_, 0.4_R_)
    end if

  end subroutine mcarFlx__tune_init
 
  !+
  ! Set technical parameters
  !-
  subroutine mcarFlx__set_tech(diff0, diff1) 

    real(R_), intent(in) :: diff0, diff1
    Flx_diff0 = diff0
    Flx_diff1 = diff1

  end subroutine mcarFlx__set_tech

  !+
  ! Read in namelist variables
  !-
  subroutine mcarFlx__user_init(iu) 

    integer, intent(in) :: iu   ! file unit index for namelist input
    integer :: ios
    read (iu, nml=mcarFlx_nml_init, iostat=ios)
    call err_read(ios, iu, 'mcarFlx__user_init: Invalid namelist input for mcarFlx_nml_init.')

  end subroutine mcarFlx__user_init


  !+
  ! Setups for fluxes & heating rates
  !-
  subroutine mcarFlx__init(nx, ny, nz, iz3l, iz3u)

    integer,  intent(in) :: nx, ny, nz, iz3l, iz3u
    integer :: iz, izf

    ! Constants
    Flx_nx = nx
    Flx_ny = ny
    Flx_nz = nz
    Flx_anxf = min(10.0_R_, 0.5_R_ * real(nx, R_))
    Flx_anyf = min(10.0_R_, 0.5_R_ * real(ny, R_))

    ! Activate samplers
    Flx_ndset = 0 ! initialize # of output datasets
    Flx_nzh = 0
    if (Flx_mhrt >= 1) Flx_nzh = nz
    if (Flx_mflx <= 0) then
       Flx_nzf = 0
    else if (Flx_mflx == 1) then
       Flx_nzf = 2
    else if (Flx_mflx == 2 .and. iz3l < iz3u) then
       Flx_nzf = 0
       do iz = 0, nz
          if (iz < iz3l .or. iz >= iz3u) then
             Flx_nzf = Flx_nzf + 1
          end if
       end do
    else
       Flx_nzf = nz + 1
    end if

    ! Allocate
    if (allocated(Flx_iizf)) call mcarFlx__final()
    allocate (Flx_iizf(0:Flx_nz))
    allocate (Flx_pflx(Flx_nx, Flx_ny, Flx_nzf, Flx_nflx))
    allocate (Flx_phrt(Flx_nx, Flx_ny, Flx_nzh))

    ! Flux index table
    Flx_iizf(:) = 0
    if (Flx_mflx <= 0) then
    else if (Flx_mflx == 1) then
       Flx_iizf(0)  = 1
       Flx_iizf(nz) = 2
    else if (Flx_mflx == 2 .and. iz3l < iz3u) then
       izf = 0
       do iz = 0, nz
          if (iz < iz3l .or. iz >= iz3u) then
             izf = izf + 1
             Flx_iizf(iz) = izf
          end if
       end do
    else
       Flx_iizf(0:nz) = gridValues_1I(nz+1, 1, 1) ! = [1:nz+1]
       !// BUGFIX, 7/28/2012, <-- gridValues_1I(nz,...
    end if

  end subroutine mcarFlx__init


  !+
  ! Finalize this module
  !-
  subroutine mcarFlx__final() 

    if (allocated(Flx_iizf)) then
       deallocate (Flx_iizf, Flx_pflx, Flx_phrt)
    end if

  end subroutine mcarFlx__final

  !+
  ! Query about calculation flags
  !-
  subroutine mcarFlx__qry_flags(mflx, mhrt) 

    integer, intent(out), optional :: mflx, mhrt
    if (present(mflx)) mflx = Flx_mflx
    if (present(mhrt)) mhrt = Flx_mhrt

  end subroutine mcarFlx__qry_flags

  !+
  ! Query about dimension sizes
  !-
  subroutine mcarFlx__qry_dims(nx, ny, nzf, nzh) 

    integer, intent(out), optional :: nx, ny, nzf, nzh
    if (present(nx )) nx  = Flx_nx
    if (present(ny )) ny  = Flx_ny
    if (present(nzf)) nzf = Flx_nzf
    if (present(nzh)) nzh = Flx_nzh

  end subroutine mcarFlx__qry_dims


  !+
  ! Zero samplers
  !-
  subroutine mcarFlx__zero() 

    Flx_psrc = 0.0_R_
    Flx_esrc = 0.0_R_
    Flx_pflx(:,:,:,:) = 0.0_RD_
    Flx_phrt(:,:,:)   = 0.0_RD_

  end subroutine mcarFlx__zero


  !+
  ! Sample # of source photons & energy
  !-
  subroutine mcarFlx__sampSrc(np) 

    integer, intent(in) :: np
    Flx_psrc = Flx_psrc + np
    Flx_esrc = Flx_esrc + np * Pho_eone

  end subroutine mcarFlx__sampSrc


  !+
  ! Sample flux contributions around a point
  !-
  subroutine mcarFlx__sampFlux(x, y, izn)

    real(R_), intent(in) :: x, y
    integer,  intent(in) :: izn
    real(R_) :: p, pwgt
    integer  :: ixf, iyf, izf, ixt, iyt, ixd, iyd, nxyd, nxdlo, nxdhi, nydlo, nydhi, iflx
    integer  :: ixlut(-Flx_KND:Flx_KND), iylut(-Flx_KND:Flx_KND)

    ! Initialize
    izf = Flx_iizf(izn)
    if (izf <= 0) return
    iflx = 2 + Pho_itr(3) ! 2 for downward, 3 for upward

    ! Sample contributions
    p = Pho_eone * Pho_wgt
    call mcarFlx__diffParams(x, y, Pho_wgt, ixf, iyf, nxyd, nxdlo, nxdhi, nydlo, nydhi)
    if (nxyd == 1) then
       Flx_pflx(ixf, iyf, izf, iflx) = Flx_pflx(ixf, iyf, izf, iflx) + p
       if (Pho_iso == 1) Flx_pflx(ixf, iyf, izf, 1) = Flx_pflx(ixf, iyf, izf, 1) + p
    else
       pwgt = p / real(nxyd, R_)
       call mcarUtl__getGtab1(ixlut, Flx_KND, nxdlo, nxdhi, Flx_nx, ixf)
       call mcarUtl__getGtab1(iylut, Flx_KND, nydlo, nydhi, Flx_ny, iyf)
       do iyd = -nydlo, nydhi
          iyt = iylut(iyd)
          do ixd = -nxdlo, nxdhi
             ixt = ixlut(ixd)
             Flx_pflx(ixt, iyt, izf, iflx) = Flx_pflx(ixt, iyt, izf, iflx) + pwgt
             if (Pho_iso == 1) Flx_pflx(ixt, iyt, izf, 1) = Flx_pflx(ixt, iyt, izf, 1) + pwgt
          end do
       end do
    end if

  end subroutine mcarFlx__sampFlux


  !+
  ! Sample contributions for upward transport fluxes from the bottom surface
  !-
  subroutine mcarFlx__sampFluxSfc()

    integer   :: ixlut(-Flx_KND:Flx_KND), iylut(-Flx_KND:Flx_KND)
    integer   :: iflx, ixf, iyf, ixd, ixt, iyd, iyt, izf, nxyd, nxdhi, nxdlo, nydhi, nydlo
    real(R_)  :: p, pwgt

    ! Initialize
    izf = Flx_iizf(0)
    if (izf <= 0) return
    iflx = 3
    call mcarFlx__diffParams(Pho_loc(1), Pho_loc(2), Pho_wgt, ixf, iyf, nxyd, nxdlo, nxdhi, nydlo, nydhi)
    p = Pho_eone * Pho_wgt

    ! Sample contributions
    if (nxyd == 1) then
       Flx_pflx(ixf, iyf, izf, iflx) = Flx_pflx(ixf, iyf, izf, iflx) + p
    else
       call mcarUtl__getGtab1(ixlut, Flx_KND, nxdlo, nxdhi, Flx_nx, ixf)
       call mcarUtl__getGtab1(iylut, Flx_KND, nydlo, nydhi, Flx_ny, iyf)
       pwgt = p / real(nxyd, R_)
       do iyd = -nydlo, nydhi
          iyt = iylut(iyd)
          do ixd = -nxdlo, nxdhi
             ixt = ixlut(ixd)
             Flx_pflx(ixt, iyt, izf, iflx) = Flx_pflx(ixt, iyt, izf, iflx) + pwgt
          end do
       end do
    end if

  end subroutine mcarFlx__sampFluxSfc


  !+
  ! Calculate numerical diffusion parameters
  !-
  subroutine mcarFlx__diffParams(x, y, wgt, ixf, iyf, nxyd, nxdlo, nxdhi, nydlo, nydhi)

    real(R_), intent(in) :: x, y
    real(R_), intent(in) :: wgt
    integer,  intent(out) :: ixf, iyf
    integer,  intent(out) :: nxyd
    integer,  intent(out) :: nxdlo, nxdhi
    integer,  intent(out) :: nydlo, nydhi
    real(R_) :: rixf, riyf, rnx, rny, rnx1, rny1, a

    ! Pixel indexes
    rixf = Atm_facx * x
    riyf = Atm_facy * y
    ixf = int(rixf) + 1
    iyf = int(riyf) + 1
    
    ! Numerical diffusion
    if (Pho_mv3d .and. (Flx_diff0 > 0.0_R_ .or. Flx_diff1 > 0.0_R_)) then
       a = sqrt(wgt) * (Flx_diff0 + Flx_diff1 * Pho_sdif)
       rnx = Atm_facx * Pho_hdwx * a
       rny = Atm_facy * Pho_hdwy * a
       if (rnx < 0.5_R_ .and. rny < 0.5_R_) then
          nxyd = 1
       else
          rnx = min(Flx_anxf, rnx)
          rny = min(Flx_anyf, rny)
          rnx1 = real(ixf, R_) - rixf - 0.5_R_
          rny1 = real(iyf, R_) - riyf + 0.5_R_
          nxdlo = int(rnx + rnx1)
          nxdhi = int(rnx - rnx1)
          nydlo = int(rny + rny1)
          nydhi = int(rny - rny1)
          nxyd = (nxdlo + nxdhi + 1) * (nydlo + nydhi + 1)
       end if
    else ! no diffusion
       nxyd = 1
    end if

  end subroutine mcarFlx__diffParams


  !+
  ! Sample heating energy power
  !-
  subroutine mcarFlx__sampHeat(ix, iy, iz, e) 

    integer,  intent(in) :: ix, iy, iz ! grid indexes
    real(R_), intent(in) :: e ! energy sample
    Flx_phrt(ix, iy, iz) = Flx_phrt(ix, iy, iz) + e

  end subroutine mcarFlx__sampHeat


  !+
  ! Integration of heating rate, possibly with numerical diffusion
  !-
  subroutine mcarFlx__sampHeatD(fdltc, coas, isup)

    real(R_), intent(in) :: coas
    real(R_), intent(in) :: fdltc
    integer,  intent(in) :: isup
    integer   :: ixlut(-Flx_KND:Flx_KND), iylut(-Flx_KND:Flx_KND)
    real(R_)  :: wrk(-Flx_KND:Flx_KND, -Flx_KND:Flx_KND)
    integer   :: ixf, iyf, ixd, ixt, iyd, iyt, nxyd, nxdhi, nxdlo, nydhi, nydlo
    real(R_)  :: pwgt, tsum, pabs

    ! Initialization
    if (Flx_nzh <= 0 .or. coas <= 0.0_R_) return
    pabs = Pho_eone * Pho_wgt * coas
    call mcarFlx__diffParams(Pho_loc(1), Pho_loc(2), Pho_wgt * fdltc, ixf, iyf, nxyd, &
         & nxdlo, nxdhi, nydlo, nydhi)
    !if (iyf > Flx_ny .or. iyf < 1) print *, iyf, Pho_loc(2)

    ! Sample with no diffusion
    if (nxyd == 1) then
       Flx_phrt(ixf, iyf, Pho_iz) = Flx_phrt(ixf, iyf, Pho_iz) + pabs
       return
    end if

    ! Sample with diffusion
    call mcarUtl__getGtab1(ixlut, Flx_KND, nxdlo, nxdhi, Flx_nx, ixf)
    call mcarUtl__getGtab1(iylut, Flx_KND, nydlo, nydhi, Flx_ny, iyf)
    if (isup <= 2) then ! 1-D super-layer
       pwgt = pabs / real(nxyd, R_)
       do iyd = -nydlo, nydhi
          iyt = iylut(iyd)
          do ixd = -nxdlo, nxdhi
             ixt = ixlut(ixd)
             Flx_phrt(ixt, iyt, Pho_iz) = Flx_phrt(ixt, iyt, Pho_iz) + pwgt
          end do
       end do
    else ! 3-D super-layer
       do iyd = -nydlo, nydhi
          iyt = iylut(iyd)
          do ixd = -nxdlo, nxdhi
             ixt = ixlut(ixd)
             wrk(ixd, iyd) = sqrt(Atm_abst3d(ixt, iyt, Pho_iz, Pho_ikd))
          end do
       end do
       tsum = sum(wrk(-nxdlo:nxdhi, -nydlo:nydhi))
       if (tsum < RSML_) then ! unexpected case
          Flx_phrt(ixf, iyf, Pho_iz) = Flx_phrt(ixf, iyf, Pho_iz) + pabs
       else
          pwgt = pabs / tsum
          do iyd = -nydlo, nydhi
             iyt = iylut(iyd)
             do ixd = -nxdlo, nxdhi
                ixt = ixlut(ixd)
                Flx_phrt(ixt, iyt, Pho_iz) = Flx_phrt(ixt, iyt, Pho_iz) + pwgt * wrk(ixd, iyd)
             end do
          end do
       end if
    end if
 
  end subroutine mcarFlx__sampHeatD


#if UseMPI == 1
  !+
  ! Reduce the paralell computed integrals
  !-
  subroutine mcarFlx__reduce(irank, iroot)

    include 'inc_mpi.f90'
    integer, intent(in) :: irank
    integer, intent(in) :: iroot
    real(R_)  :: wrk1(Flx_nx * Flx_ny), wrk0(Flx_nx * Flx_ny), s1(2), s2(2)
    integer   :: ierr, iflx, izf, izh, na(1), nb(2)

    ! Array sizes
    nb(1) = Flx_nx
    nb(2) = Flx_ny
    na(1) = nb(1) * nb(2)

    ! Fluxes
    do iflx = 1, Flx_nflx
       do izf = 1, Flx_nzf
          wrk0(1:na(1)) = reshape(Flx_pflx(:, :, izf, iflx), na)
          call MPI_Reduce(wrk0, wrk1, na(1), MPI_R_, MPI_SUM, iroot, MPI_COMM_WORLD, ierr)
          if (irank == iroot) Flx_pflx(:, :, izf, iflx) = reshape(wrk1(1:na(1)), nb)
       end do
    end do

    ! Heating rates
    do izh = 1, Flx_nzh
       wrk0(1:na(1)) = reshape(Flx_phrt(:, :, izh), na)
       call MPI_Reduce(wrk0, wrk1, na(1), MPI_R_, MPI_SUM, iroot, MPI_COMM_WORLD, ierr)
       if (irank == iroot) Flx_phrt(:, :, izh) = reshape(wrk1(1:na(1)), nb)
    end do

    ! Sources
    s1(1) = Flx_psrc
    s1(2) = Flx_esrc
    call MPI_Reduce(s1, s2, 2, MPI_R_, MPI_SUM, iroot, MPI_COMM_WORLD, ierr)
    if (irank == iroot) then
       Flx_psrc = s2(1)
       Flx_esrc = s2(2)
    end if

  end subroutine mcarFlx__reduce
#endif


  !+
  ! Normalize results
  !-
  subroutine mcarFlx__normal(atot) 

    real(R_),  intent(in) :: atot   ! total horizontal area (m^2) of the domain
    real(RD_) :: fnorm
    integer   :: izh

    ! Normalize fluxes
    Flx_esrc = Flx_esrc / Flx_psrc ! average source power [W]
    if (Flx_nzf > 0) then
       fnorm = real(Flx_nx * Flx_ny, RD_) / (atot * Flx_psrc) ! [/m^2]
       Flx_pflx(:,:,:,:) = fnorm * Flx_pflx(:,:,:,:)
    end if

    ! Normalize heating rates
    do izh = 1, Flx_nzh
       fnorm = real(Flx_nx * Flx_ny, RD_) / (atot * Flx_psrc * Atm_zdep(izh)) ! [/m^3]
       Flx_phrt(:,:,izh) = fnorm * Flx_phrt(:,:,izh)
    end do

  end subroutine mcarFlx__normal


  !+
  ! Write out results to a binary fil
  !-
  subroutine mcarFlx__writeBin(iuo) 

    integer,   intent(in) :: iuo ! file unit index for output
    real(R4_), allocatable  :: wrk(:,:)
    integer   :: izf, izh, iflx

    ! Write out
    allocate (wrk(Flx_nx, Flx_ny))
    do iflx = 1, Flx_nflx ! fluxes
       do izf = 1, Flx_nzf
          wrk(:, :) = Flx_pflx(:, :, izf, iflx)
          call bin_write_i2R4(iuo, wrk)
       end do
    end do
    do izh = 1, Flx_nzh ! heating rates
       wrk(:, :) =  Flx_phrt(:, :, izh)
       call bin_write_i2R4(iuo, wrk)
    end do
    deallocate (wrk)
    Flx_ndset = Flx_ndset + 1 ! count # of output datasets

  end subroutine mcarFlx__writeBin


  !+
  ! Write out results to a text fil
  !-
  subroutine mcarFlx__writeTxt(iuo) 

    integer,   intent(in) :: iuo    ! file unit index for output
    real(R4_), allocatable  :: wrk(:,:)
    integer   :: izf, izh, iflx, ix, iy

    ! Write out flux
    allocate (wrk(Flx_nx, Flx_ny))
    if (Flx_nzf > 0) then
       do iflx = 1, Flx_nflx
          write (iuo, '(a, i4)', err=1) '# iflx= ', iflx
          do izf = 1, Flx_nzf
             wrk(:, :) = Flx_pflx(:, :, izf, iflx)
             write (iuo, *, err=1) ((wrk(ix, iy), ix = 1, Flx_nx), iy = 1, Flx_ny)
          end do
       end do
    end if

    ! Write out heating rate
    if (Flx_nzh > 0) then
       write (iuo, '(a, i4)', err=1) '# heating rate'
       do izh = 1, Flx_nzh
          wrk(:, :) =  Flx_phrt(:, :, izh)
          write (iuo, *, err=1) ((wrk(ix, iy), ix = 1, Flx_nx), iy = 1, Flx_ny)
       end do
    end if
    deallocate (wrk)

    return
1   call err_write(1, iuo, 'mcarFlx__writeTxt')

  end subroutine mcarFlx__writeTxt


  !+
  ! Make a (GrADS) control file
  !-
  subroutine mcarFlx__writeCtl(iuc, outfile) 

    integer, intent(in) :: iuc
    character(*), intent(in) :: outfile
    integer :: nz, nvar, i
    character(256) :: varstr(4)

    nz = max(Flx_nzf, Flx_nzh)
    nvar = 0
    if (Flx_nzf >= 1) then
       do i = 1, 3
          nvar = nvar + 1
          varstr(nvar) = 'a'//trim(num2str_AN(i))//' '//trim(num2str_AN(Flx_nzf))//' 99 Flux Density'
       end do
    end if
    if (Flx_nzh >= 1) then
       nvar = nvar + 1
       varstr(nvar) = 'b'//trim(num2str_AN(1))//' '//trim(num2str_AN(Flx_nzh))//' 99 Heating Rate'
    end if
    call gradsCtl_write(iuc, fileName_AN(outfile), Flx_nx, Flx_ny, nz, Flx_ndset, nvar, varstr)

  end subroutine mcarFlx__writeCtl

end module mcarFlx
