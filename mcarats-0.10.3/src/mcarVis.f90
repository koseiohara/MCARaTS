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
! Module for visualization of the atmosphre-surface system
!-
module mcarVis 

  use globals
  use hparx
  use mcarAtm
  use mcarSfc
  use mcarSca
  implicit none
  private

  ! Public
  public :: mcarVis__user_init
  public :: mcarVis__init
  public :: mcarVis__final
  public :: mcarVis__prep_job
  public :: mcarVis__prep_prop
  public :: mcarVis__prep_srcj
  public :: mcarVis__render

  ! (Private) User variables in namelist
  integer,  save :: Vis_mrend = 2 ! method for rendering (0/1/2/3/4)
  !// = 0 : integration of density along the ray (e.g. for calculating optical thickness)
  !     1 : as 1, but with attenuation (assuming that the source function is uniformly 1)
  !     2 : RTE-based solver, assuming horizontally-uniform source functions
  !     3 : as 2, but taking into accout 3-D distribution of 1st-order source functions (J1)
  !     4 : as 3, but with TMS correction for anisotropic scattering
  real(R_), save :: Vis_epserr = 1.0e-4_R_ ! convergence criterion for multiple scattering components
  real(R_), save :: Vis_fpsmth = 0.5_R_ ! phase function smoothing fraction for relaxed TMS correction
  real(R_), save :: Vis_fatten = 1.0_R_ ! attenuation factor (1 for physics-based rendering)
  integer,  save :: Vis_nqhem = 1       ! # of Gaussian quadrature points in a hemisphere
  integer,  save :: Vis_nqlay = 10      ! # of Gaussian quadrature points per layer
  namelist /mcarVis_nml_init/Vis_mrend, Vis_epserr, Vis_fpsmth, Vis_nqhem, Vis_nqlay, Vis_fatten

  ! Private
  integer,  save :: Vis_nx = 1
  integer,  save :: Vis_ny = 1
  integer,  save :: Vis_nz = 1
  integer,  save :: Vis_nxb = 1
  integer,  save :: Vis_nyb = 1
  integer,  save :: Vis_ngsmax = 10  ! max factor to limit ngsray
  real(R_), save :: Vis_facx, Vis_facy
  real(R_), save :: Vis_xmax, Vis_ymax
  real(R_), save :: Vis_sdir(3)         ! solar source transport direction vector
  real(R_), save, allocatable :: Vis_xgrd(:)     !(0:nx)
  real(R_), save, allocatable :: Vis_ygrd(:)     !(0:ny)
  real(R_), save, allocatable :: Vis_zgrd(:)     !(0:nz) Z-node values
  real(R_), save, allocatable :: Vis_zdep(:)     !(nz) layer depth (m)
  real(R_), save, allocatable :: Vis_dens(:,:,:) !(nx,ny,nz) some density
  real(R_), save, allocatable :: Vis_aalb(:,:,:) !(nx,ny,nz) single scattering albedo
  real(R_), save, allocatable :: Vis_srcj(:,:,:) !(nx,ny,nz) source function
  real(R_), save, allocatable :: Vis_acor(:,:,:) !(nx,ny,nz) function for TMS correction
  !// For TMS correction, acor = F0 / (4*pi) * exp(-tau0#) / Bext#, where # means scaled property.
  real(R_), save, allocatable :: Vis_splk(:,:)   !(nxb,nyb) Plank function at the surface for TR(1)
  real(R_), save, allocatable :: Vis_sfs1(:,:)   !(nx,ny) downward irradiance at the surface for SR(1)
  real(R_), save, allocatable :: Vis_sf2p(:,:)   !(nx,ny) downward irradiance at the surface for R(2+)
  integer,  save, allocatable :: Vis_ngsray(:)   !(nz) # of quadrature points
  real(R_), save, allocatable :: Vis_xgsray(:,:) !(ngsray,nz) quadrature abscissas
  real(R_), save, allocatable :: Vis_wgsray(:,:) !(ngsray,nz) quadrature weights
  real(R_), save, allocatable :: Vis_dgsray(:,:) !(2,nz) directional cosine threshold
  !// for Gaussian quadrature

contains

  !+
  ! Read in namelist variables
  !-
  subroutine mcarVis__user_init(iu) 

    integer, intent(in) :: iu   ! file unit index for namelist input
    integer :: ios
    read (iu, nml=mcarVis_nml_init, iostat=ios)
    call err_read(ios, iu, 'mcarVis__user_init: Invalid namelist input for mcarVis_nml_init.')

  end subroutine mcarVis__user_init


  !+
  ! Initialize this module
  !-
  subroutine mcarVis__init() 

    ! Assign
    Vis_nx = Atm_nx
    Vis_ny = Atm_ny
    Vis_nz = Atm_nz
    Vis_nxb = Sfc_nxb
    Vis_nyb = Sfc_nyb

    ! Allocate
    if (allocated(Vis_zgrd)) call mcarVis__final()
    allocate (Vis_xgrd(0:Vis_nx), Vis_ygrd(0:Vis_ny))
    allocate (Vis_zgrd(0:Vis_nz), Vis_zdep(Vis_nz))
    allocate (Vis_dens(Vis_nx, Vis_ny, Vis_nz))
    allocate (Vis_aalb(Vis_nx, Vis_ny, Vis_nz))
    allocate (Vis_srcj(Vis_nx, Vis_ny, Vis_nz))
    allocate (Vis_splk(Vis_nxb, Vis_nyb))
    allocate (Vis_acor(Vis_nx, Vis_ny, Vis_nz))
    allocate (Vis_sfs1(Vis_nx, Vis_ny))
    allocate (Vis_sf2p(Vis_nx, Vis_ny))
    allocate (Vis_ngsray(Vis_nz))
    allocate (Vis_xgsray(Vis_ngsmax * Vis_nqlay, Vis_nz))
    allocate (Vis_wgsray(Vis_ngsmax * Vis_nqlay, Vis_nz))
    allocate (Vis_dgsray(2, Vis_nz))

  end subroutine mcarVis__init
    

  !+
  ! Finalize this module
  !-
  subroutine mcarVis__final() 

    if (allocated(Vis_xgrd)) then
       deallocate (Vis_xgrd, Vis_ygrd, Vis_zgrd, Vis_zdep, Vis_dens)
       deallocate (Vis_aalb, Vis_srcj, Vis_splk, Vis_acor, Vis_sfs1, Vis_sf2p)
       deallocate (Vis_ngsray, Vis_xgsray, Vis_wgsray, Vis_dgsray)
    end if

  end subroutine mcarVis__final


  !+
  ! Prepare for core calculations
  !-
  subroutine mcarVis__prep_job() 

    integer :: iz

    ! Initialize
    Vis_xmax = Atm_xmax
    Vis_ymax = Atm_ymax
    Vis_facx = Vis_nx / Vis_xmax * 0.9999995_R_
    Vis_facy = Vis_ny / Vis_ymax * 0.9999995_R_
    Vis_xgrd(0:Vis_nx) = Atm_xgrd(0:Vis_nx)
    Vis_ygrd(0:Vis_ny) = Atm_ygrd(0:Vis_ny)
    Vis_zgrd(0:Vis_nz) = Atm_zgrd(0:Vis_nz)
    Vis_zdep(1:Vis_nz) = Atm_zdep(1:Vis_nz)

    ! Quadrature variables
    do iz = 1, Vis_nz
       Vis_ngsray(iz) = Vis_nqlay * min(Vis_ngsmax, int(1.0_R_ + Atm_extvar(iz)))
       call gaussLegen(Vis_ngsray(iz), Vis_xgsray(:, iz), Vis_wgsray(:, iz))
    end do
    Vis_xgsray(:, :) = 0.5_R_ * (1.0_R_ + Vis_xgsray(:, :)) ! [0,1]
    Vis_dgsray(1, :) = 2.0_R_ * Vis_ngsray(:) * Vis_xmax / Vis_nx / Vis_zdep(:) ! tanQx = dx/dz
    Vis_dgsray(2, :) = 2.0_R_ * Vis_ngsray(:) * Vis_ymax / Vis_ny / Vis_zdep(:) ! tanQy = dy/dz

  end subroutine mcarVis__prep_job


  !+
  ! Synthesize a density distribution
  !-
  subroutine mcarVis__prep_prop(mdens, ikd)

    integer,  intent(in) :: mdens ! flag of a kind of density
    integer,  intent(in) :: ikd   ! index of k-distribution term
    !// = 1 : dens = extinction coefficient scaled
    !     2 : dens = scattering coefficient scaled, aalb = 1
    !     3 : dens = extinction coefficient (non-scaled)
    !     4 : dens = scattering coefficient (non-scaled), aalb = 1
    !     5 : dens = absorption coefficient, aalb = 0
    integer  :: ippo, ix, iy, iz

    ! Scattering & absorption coefficients
    if (mdens <= 2) then
       Vis_aalb(:,:,:) = 0.0_R_
       do iz = 1, Vis_nz
          do ippo = 1, Atm_nppo(iz)
             do iy = 1, Vis_ny
                do ix = 1, Vis_nx
                   Vis_aalb(ix, iy, iz) = Vis_aalb(ix, iy, iz) + Atm_scap3d(ix, iy, iz, ippo) &
                        & * (1.0_R_ - mcarSca__asymStar_R(Atm_apfp3d(ix, iy, iz, ippo), 1))
                end do
             end do
          end do
       end do
    else
       Vis_aalb(:,:,:) = Atm_scat3d(:,:,:,1)
    end if
    Vis_dens(:,:,:) = Vis_aalb(:,:,:) + Atm_abst3d(:,:,:,ikd)
    !// Now, aalb & dens are respectively scattering & extinction coefficients.

    ! Results
    if (mdens == 1 .or. mdens == 3) then
       Vis_aalb(:,:,:) = Vis_aalb(:,:,:) / Vis_dens(:,:,:) ! single scattering albedo
    else if (mdens == 2 .or. mdens == 4) then
       Vis_dens(:,:,:) = Vis_aalb(:,:,:)
       Vis_aalb(:,:,:) = 1.0_R_
    else
       Vis_dens(:,:,:) = Vis_dens(:,:,:) - Vis_aalb(:,:,:)
       Vis_aalb(:,:,:) = 0.0_R_
    end if

  end subroutine mcarVis__prep_prop


  !+
  ! Prepare source functions
  !-
  subroutine mcarVis__prep_srcj(mstype, sdir, fsrc) 

    integer,  intent(in) :: mstype  ! source type (1/2/3)
    real(R_), intent(in) :: sdir(:) ! source transport direction vector (downward)
    real(R_), intent(in) :: fsrc    ! solar incident irradiance
    real(R_), parameter :: API4 = 1.0_R_ / (4.0_R_ * PI_)
    real(R_), parameter :: API  = 1.0_R_ / PI_
    integer,  parameter :: NPASS_MAX = 100 ! max # of iterations for scattering source functions
    real(R_), allocatable :: tau0(:,:,:)
    real(R_)  :: dtau(Vis_nz), alb(0:Vis_nz), sold(0:Vis_nz), snew(0:Vis_nz) !(AUTO)
    real(R_)  :: r0, t, ta, tb, tc, td, ax0, ay0, dx, dy, path, x0, y0, pfsrc
    integer   :: ix0, iy0, ix0s, iy0s, ix, iy, iz, ixb, iyb

    ! Uniform distribution
    if (Vis_mrend <= 1) then
       Vis_srcj(:,:,:) = 1.0_R_
       Vis_splk(:,:) = 0.0_R_
       return
    end if

    ! Initialization
    Vis_sdir(:) = sdir(:) ! save the solar transport direction
    Vis_sdir(3) = nonZero_R(Vis_sdir(3), REPS_)
    pfsrc = fsrc * abs(Vis_sdir(3))
    alb(0) = sum(Sfc_salbs(1,:,:)) / (Vis_nxb * Vis_nyb)! average surface albedo
    do iz = 1, Vis_nz
       dtau(iz) = sum(Vis_dens(:, :, iz)) / (Vis_nx * Vis_ny) * Vis_zdep(iz) ! average dTAU
       alb(iz)  = sum(Vis_aalb(:, :, iz)) / (Vis_nx * Vis_ny) ! average albedo
    end do

    ! SR(0) & SR(1)
    Vis_srcj(:,:,:) = 0.0_R_
    if (mstype <= 2 .and. Vis_sdir(3) < 0.0_R_) then
       if (Vis_mrend == 2) then ! 1-D approximation
          t = 0.0_R_
          do iz = Vis_nz, 1, -1
             r0 = fsrc * API4 * exp(-(t + dtau(iz) * 0.5_R_) / abs(Vis_sdir(3)))
             Vis_srcj(:, :, iz) = Vis_srcj(:, :, iz) + Vis_aalb(:, :, iz) * r0
             t = t + dtau(iz)
          end do
          Vis_sfs1(:,:) = pfsrc * exp(-t / abs(Vis_sdir(3))) ! SR(0) irradiance
       else ! explicit calculations
          allocate (tau0(Vis_nx + 1, Vis_ny + 1, Vis_nz + 1))
          call mcarVis__dBeams(Vis_sdir, tau0)

          ! Atmospheric SR(1) source, 3-D
          do iz = 1, Vis_nz
             path = ((Vis_zgrd(iz - 1) + Vis_zgrd(iz)) * 0.5_R_ - Vis_zgrd(Vis_nz)) / Vis_sdir(3)
             dx = path * Vis_sdir(1)
             dy = path * Vis_sdir(2)
             x0 = Vis_xmax / Vis_nx * 0.5_R_ - dx ! x0 for (ix,iy)=(1,1)
             y0 = Vis_ymax / Vis_ny * 0.5_R_ - dy ! y0 for (ix,iy)=(1,1)
             if (x0 < 0.0_R_ .or. x0 >= Vis_xmax) then ! cycle the X location
                x0 = mod(x0, Vis_xmax)
                if (x0 < 0.0_R_) x0 = x0 + Vis_xmax
             end if
             if (y0 < 0.0_R_ .or. y0 >= Vis_ymax) then ! cycle the Y location
                y0 = mod(y0, Vis_ymax)
                if (y0 < 0.0_R_) y0 = y0 + Vis_ymax
             end if
             ax0 = x0 * Vis_facx + 0.5_R_
             ay0 = y0 * Vis_facy + 0.5_R_
             ix0s = min(Vis_nx, int(ax0)) ! ix0 for (ix,iy)=(1,1)
             iy0s = min(Vis_ny, int(ay0)) ! iy0 for (ix,iy)=(1,1)
             ax0 = ax0 - ix0s ! interpolation factor for X
             ay0 = ay0 - iy0s ! interpolation factor for Y
             if (ix0s == 0) ix0s = Vis_nx ! cycle for a side margin
             if (iy0s == 0) iy0s = Vis_ny ! cycle for a side margin
             do iy = 1, Vis_ny
                iy0 = iy0s + iy - 1
                if (iy0 > Vis_ny) iy0 = iy0 - Vis_ny
                do ix = 1, Vis_nx
                   ix0 = ix0s + ix - 1
                   if (ix0 > Vis_nx) ix0 = ix0 - Vis_nx
                   ta = (1.0_R_ - ax0) * tau0(ix0, iy0,   iz)   + ax0 * tau0(ix0+1, iy0,   iz)
                   tb = (1.0_R_ - ax0) * tau0(ix0, iy0+1, iz)   + ax0 * tau0(ix0+1, iy0+1, iz)
                   tc = (1.0_R_ - ax0) * tau0(ix0, iy0,   iz+1) + ax0 * tau0(ix0+1, iy0,   iz+1)
                   td = (1.0_R_ - ax0) * tau0(ix0, iy0+1, iz+1) + ax0 * tau0(ix0+1, iy0+1, iz+1)
                   t = ((1.0_R_ - ay0) * (ta + tc) + ay0 * (tb + td)) * 0.5_R_ ! trilinear
                   r0 = fsrc * API4 * exp(-t)
                   t = Vis_aalb(ix, iy, iz) * r0
                   Vis_srcj(ix, iy, iz) = Vis_srcj(ix, iy, iz) + t
                   if (Vis_mrend == 4) Vis_acor(ix, iy, iz) = t
                end do
             end do
          end do
          
          ! Surface irradiance of direct beams, SR(0), 2-D
          path = (Vis_zgrd(0) - Vis_zgrd(Vis_nz)) / Vis_sdir(3)
          dx = path * Vis_sdir(1)
          dy = path * Vis_sdir(2)
          x0 = Vis_xmax / Vis_nx * 0.5_R_ - dx
          y0 = Vis_ymax / Vis_ny * 0.5_R_ - dy
          if (x0 < 0.0_R_ .or. x0 >= Vis_xmax) then ! cycle the X location
             x0 = mod(x0, Vis_xmax)
             if (x0 < 0.0_R_) x0 = x0 + Vis_xmax
          end if
          if (y0 < 0.0_R_ .or. y0 >= Vis_ymax) then ! cycle the Y location
             y0 = mod(y0, Vis_ymax)
             if (y0 < 0.0_R_) y0 = y0 + Vis_ymax
          end if
          ax0 = x0 * Vis_facx + 1.0_R_
          ay0 = y0 * Vis_facy + 1.0_R_
          ix0s = min(Vis_nx, int(ax0)) ! ix0 for (ix,iy)=(1,1)
          iy0s = min(Vis_ny, int(ay0)) ! iy0 for (ix,iy)=(1,1)
          ax0 = ax0 - ix0s ! interpolation factor for X
          ay0 = ay0 - iy0s ! interpolation factor for Y
          do iy = 1, Vis_ny
             iy0 = iy0s + iy - 1
             if (iy0 > Vis_ny) iy0 = iy0 - Vis_ny
             do ix = 1, Vis_nx
                ix0 = ix0s + ix - 1
                if (ix0 > Vis_nx) ix0 = ix0 - Vis_nx
                iz = 1
                ta = (1.0_R_ - ax0) * tau0(ix0, iy0,   iz) + ax0 * tau0(ix0+1, iy0,   iz)
                tb = (1.0_R_ - ax0) * tau0(ix0, iy0+1, iz) + ax0 * tau0(ix0+1, iy0+1, iz)
                t = (1.0_R_ - ay0) * ta + ay0 * tb ! bilinear
                Vis_sfs1(ix, iy) = pfsrc * exp(-t) ! SR(0) downward irradiance for SR(1)
             end do
          end do
          deallocate (tau0)
       end if
    end if

    ! Sources of TR(1)
    if (mstype >= 2) then
       do iz = 1, Vis_nz ! atmosphere, 3-D
          do iy = 1, Vis_ny
             do ix = 1, Vis_nx
                t = (Atm_tmpa3d(ix, iy, iz) + Atm_tmpa3d(ix, iy, iz + 1)) * 0.5_R_ ! average T
                Vis_srcj(ix, iy, iz) = Vis_srcj(ix, iy, iz) &
                     & + (1.0_R_ - Vis_aalb(ix, iy, iz)) * plkTab_intp_R(t, 1)
             end do
          end do
       end do
       do iyb = 1, Vis_nyb ! surface, 2-D
          do ixb = 1, Vis_nxb
             Vis_splk(ixb, iyb) = plkTab_intp_R(Sfc_tmps2d(ixb, iyb), 1)
          end do
       end do
    end if

    ! Solve a 1-D RTE
    sold(0) = alb(0) * sum(Vis_sfs1(:,:)) / (Vis_nx * Vis_ny) ! average of SR(1)
    sold(0) = sold(0) + sum(Vis_splk(:,:) * (1.0_R_ - Sfc_salbs(1,:,:))) / (Vis_nxb * Vis_nyb)
    !// add average of TR(1)
    do iz = 1, Vis_nz
       sold(iz) = sum(Vis_srcj(:,:,iz)) / (Vis_nx * Vis_ny)
    end do !print *, 'Initial source: ', sold(:)
    call mcarVis__solveRT_PPH(dtau, alb, Vis_nqhem, Vis_epserr, NPASS_MAX, sold, snew)
    !print *, 'Final, fraction of multiple scattering components', (snew(:) - sold(:)) / sold(:)

    ! Add parameterized multiple scattering components: SR(2+) + TR(2+)
    r0 = PI_ * (snew(0) - sold(0)) / alb(0) ! downward irradiance for SR(2+) + TR(2+)
    if (Vis_mrend == 2) then ! 1-D approximation
       Vis_sf2p(:,:) = r0
    else ! approximated 3-D
       t = sold(0) / snew(0) * 0.5_R_ ! a power exponent in (0,1]/2, modified on Nov 7, 2010
       Vis_sf2p(:,:) = r0 * ((1.0_R_ - t) &
            * (0.25_R_ + 0.75_R_ * sqrt(Vis_dens(:,:,1) * Vis_zdep(1) / dtau(1))) &
            + t * Vis_srcj(:,:,1) / sold(1))**(t**2)
       !// depend on the bottom layer
    end if
    do iz = 1, Vis_nz
       r0 = (snew(iz) - sold(iz)) / alb(iz) ! mean radiance of SR(2+) + TR(2+)
       if (Vis_mrend == 2) then ! 1-D approximation
          Vis_srcj(:,:,iz) = Vis_srcj(:,:,iz) + Vis_aalb(:,:,iz) * r0
       else ! approximated 3-D
          t = sold(iz) / snew(iz) * 0.5_R_ ! a power exponent in (0,1]/2, modified on Nov 7, 2010
          Vis_srcj(:,:,iz) = Vis_srcj(:,:,iz) + Vis_aalb(:,:,iz) * r0 &
               * ((1.0_R_ - t) * (0.25_R_ + 0.75_R_ * sqrt(Vis_dens(:,:,iz) * Vis_zdep(iz) / dtau(iz))) &
               + t * Vis_srcj(:,:,iz) / sold(iz))**(t**2)
       end if
    end do

    ! Prepare for the TMS correction
    if (mstype <= 2 .and. Vis_sdir(3) < 0.0_R_ .and. Vis_mrend == 4) then
       Vis_srcj(:,:,:) = Vis_srcj(:,:,:) - Vis_acor(:,:,:) ! subtract the SR(1) component
       Vis_acor(:,:,:) = Vis_acor(:,:,:) / (Vis_aalb(:,:,:) * Vis_dens(:,:,:))
       !// divide by scaled scattering coefficient
    end if

  end subroutine mcarVis__prep_srcj


  !+
  ! Solve a RTE under the transport approximation, by the Picard iteration
  !-
  subroutine mcarVis__solveRT_PPH(dtau, alb, nqh, epse, npass, sold, snew)

    real(R_), intent(in) :: dtau(:) ! layer optical thickness
    real(R_), intent(in) :: alb(0:) ! single scattering albedo
    integer,  intent(in) :: nqh     ! # of quadrature points in a hemisphere
    real(R_), intent(in) :: epse    ! convergence criterion (RMS difference of source function)
    integer,  intent(in) :: npass   ! max # of calculation passes
    real(R_), intent(in) :: sold(0:) ! old source functions
    real(R_), intent(out) :: snew(0:) ! new source functions
    real(R_) :: xq(nqh*2), wq(nqh*2), stau, tau, rold, rnew, e, rsfc !(AUTO)
    real(R_) :: wrk(0:size(sold)-1) ! work vector used for mean radiances or source functions
    integer  :: iz, nz, iqh, iqb, ipass

    ! Initialize
    nz = size(dtau)
    call gaussLegen(nqh*2, xq, wq)
    snew(:) = sold(:)

    ! Loop for multiple passes
    do ipass = 1, npass

       ! New mean radiances
       rsfc = 0.0_R_
       wrk(:) = 0.0_R_ ! Here, wrk(0:nz) are tentatively mean radiances
       do iqh = 1, nqh
          iqb = nqh*2 - iqh + 1
          rnew = 0.0_R_
          wrk(nz) = wrk(nz) + wq(iqh) * rnew
          rold = rnew ! downward at TOA
          tau = 0.0_R_
          do iz = nz, 1, -1 ! downward
             stau = dtau(iz) / abs(xq(iqh))
             rnew = rold * exp(-tau) + snew(iz) * mtab_expNXC_R(stau)
             wrk(iz - 1) = wrk(iz - 1) + wq(iqh) * rnew
             rold = rnew
             tau = tau + stau
          end do
          rnew = rold * alb(0) + sold(0) ! upward at BOA
          rsfc = rsfc + 2.0_R_ * wq(iqb) * rnew
          wrk(0) = wrk(0) + wq(iqb) * rnew
          rold = rnew
          tau = 0.0_R_
          do iz = 1, nz ! upward
             stau = dtau(iz) / abs(xq(iqb))
             rnew = rold * exp(-tau) + snew(iz) * mtab_expNXC_R(stau)
             wrk(iz) = wrk(iz) + wq(iqb) * rnew
             rold = rnew
             tau = tau + stau
          end do
       end do
    
       ! New source functions & test convergence
       wrk(0) = rsfc ! upward mean radiance at BOA
       wrk(1:nz) = sold(1:nz) + (wrk(0:nz-1) + wrk(1:nz)) * 0.5_R_ * alb(1:nz)
       !// Here, wrk(1:nz) have been changed to source functions
       e = sqrt(sum(((wrk(:) - snew(:)) / wrk(:))**2) / (nz + 1)) ! RMS error
       snew(:) = wrk(:)
       !print *, 'Pass', ipass, ', error =', e
       if (e < epse) exit
    end do

  end subroutine mcarVis__solveRT_PPH


  !+
  ! Calculate optical thicknesses for solar direct beams
  !-
  subroutine mcarVis__dBeams(dir, tau0)

    real(R_), intent(in) :: dir(:)
    real(R_), intent(out) :: tau0(:,:,:)
    integer  :: ix, iy, iz, itr(3), ix0, iy0
    real(R_) :: loc(3), tau, rad
    real(R_), parameter :: TAUMAX = 80.0_R_

    ! Initialize
    itr(:) = 0
    if (dir(1) >= 0.0_R_) itr(1) = 1
    if (dir(2) >= 0.0_R_) itr(2) = 1
    if (dir(3) >= 0.0_R_) itr(3) = 1
    if (itr(3) == 1) then
       tau0(:,:,:) = RLRG_ ! no solar source
       return
    end if

    ! Trace rays
    do iy0 = 1, Vis_ny
       do ix0 = 1, Vis_nx
          loc(1) = Vis_xmax / Vis_nx * (ix0 - 0.5_R_)
          loc(2) = Vis_ymax / Vis_ny * (iy0 - 0.5_R_)
          loc(3) = Vis_zgrd(Vis_nz)
          !// Initial locations are at the center of each (x,y) pixel.
          ix = 0
          iz = Vis_nz
          tau = 0.0_R_
          tau0(ix0, iy0, iz + 1) = tau
          do ! loop for layers (top to bottom)
             call mcarVis__trace(0, dir, itr, 0.0_R_, 1.0_R_, 0.0_R_, loc, ix, iy, iz, tau, rad)
             tau0(ix0, iy0, iz + 1) = min(TAUMAX, tau)
             if (tau < -1.0_R_) call err_issue(1, 'mcarVis__dBeams: Negative TAU')
             if (iz <= 0) exit
          end do
       end do
    end do

    ! Copy margins under cyclic boundary conditions
    tau0(Vis_nx + 1, :, :) = tau0(1, :, :)
    tau0(:, Vis_ny + 1, :) = tau0(:, 1, :)

  end subroutine mcarVis__dBeams


  !+
  ! Core of volume rendering : Calculate (pseudo-) radiance for a single ray
  !-
  subroutine mcarVis__render(taumax, dir, loc, tau, rad)

    real(R_), intent(in) :: taumax
    real(R_), intent(in) :: dir(:)  ! should not be 0
    real(R_), intent(inout) :: loc(:)
    real(R_), intent(out) :: tau ! optical thickness
    real(R_), intent(out) :: rad ! (pseudo-)radiance
    integer  :: ix, iy, iz, ixb, iyb, itr(3), isfc
    real(R_) :: csca, qsca, path, w, pac(10), yco(2), yfa1, yfa2

    ! Initialize
    tau = 0.0_R_
    rad = 0.0_R_
    iz = gridIdx_bin0_I(Vis_zgrd, loc(3), 0, Vis_nz+1) ! [0, nz+1]
    itr(:) = 0
    if (dir(1) >= 0.0_R_) itr(1) = 1
    if (dir(2) >= 0.0_R_) itr(2) = 1
    if (dir(3) >= 0.0_R_) itr(3) = 1
    if (Vis_mrend == 4) then
       csca = geo3D_scaProd_R(Vis_sdir, -dir) ! cos(scattering angle)
       csca = max(-1.0_R_, min(1.0_R_, csca))
       qsca = acos(csca)
    end if

    ! Entering the atmosphere from outside
    if (iz <= 0 .or. iz >= Vis_nz + 1) then
       if ((iz <= 0 .and. itr(3) == 0) .or. (iz >= Vis_nz + 1 .and. itr(3) == 1)) return
       !// viewing space or underground
       path = (Vis_zgrd(iz - 1 + itr(3)) - loc(3)) / dir(3)
       call mcarVis__rt_bound1D(itr(3), iz, loc(3)) ! to the TOA/BOA boundary
       call mcarVis__rt_renewLocXY(dir, path, loc) 
    end if

    ! Trace the ray
    ix = 0
    do
       call mcarVis__trace(Vis_mrend, dir, itr, qsca, csca, Vis_fpsmth, loc, ix, iy, iz, tau, rad)
       if (tau > taumax .or. iz <= 0 .or. iz >= Vis_nz + 1) exit
    end do

    ! Result, adding surface contribution
    if (Vis_mrend <= 0) then
       rad = tau
    else if (Vis_mrend >= 2 .and. tau <= taumax .and. iz <= 0) then ! add surface contributions
       ixb = min(Vis_nxb, int(loc(1) / Vis_xmax * Vis_nxb) + 1)
       iyb = min(Vis_nyb, int(loc(2) / Vis_ymax * Vis_nyb) + 1)
       isfc = Sfc_jsfc2d(ixb, iyb)
       if (isfc >= 1) then
          if (ix == 0) then
             ix = min(Vis_nx, int(loc(1) * Vis_facx) + 1)
             iy = min(Vis_nx, int(loc(2) * Vis_facy) + 1)
          end if
          pac(1:5) = Sfc_psfc2d(ixb, iyb, 1:5)
          if (isfc == 1) then ! Lambertian
             pac(6) = pac(1)
          else ! DSM, RPV, or LSRT
             call mcarSfc__intpOptProp(0, isfc, ixb, iyb, -dir(3), yfa1, yfa2, yco)
             if (isfc == 3) pac(5) = yfa1
             pac(6) = mcarSfc__bsAlbedo_R(isfc, pac, yfa1, yfa2, -dir(3)) ! albedo (can be > 1)
          end if
          w = exp(-tau)
          rad = rad + ((1.0_R_ - pac(6)) * Vis_splk(ixb, iyb) + pac(6) * Vis_sf2p(ix, iy) / PI_) * w
          !// add SR(2+) + TR(1+)
          if (Vis_sdir(3) < 0.0_R_) then ! BRDF for SR(1)
             if (isfc == 1) then ! Lambertian
                pac(6) = pac(1)
             else ! DSM, RPV, or LSRT
                call mcarSfc__intpOptProp(0, isfc, ixb, iyb, -Vis_sdir(3), yfa1, yfa2, yco)
                if (isfc == 3) pac(5) = yfa1
                pac(6) = mcarSfc__bsAlbedo_R(isfc, pac, yfa1, yfa2, -Vis_sdir(3)) ! albedo (can be > 1)
             end if
             w = w * mcarSfc__angDistr_R(4, isfc, pac, Vis_sdir, -dir, ixb, iyb)
             rad = rad + Vis_sfs1(ix, iy) * w * pac(6) / abs(dir(3)) ! add SR(1)
          end if
       end if
    end if

  end subroutine mcarVis__render


  !###
  !+
  ! Trace a ray in a single atmospheric layer
  !-
  subroutine mcarVis__trace(mrend, dir, itr, qsca, csca, fpsmth, loc, ix, iy, iz, tau, rad)

    integer,  intent(in) :: mrend ! method for rendering (0/1/2/3/4)
    real(R_), intent(in) :: dir(:)
    integer,  intent(in) :: itr(:)
    real(R_), intent(in) :: qsca, csca ! scattering angle (radian) and cosine, used only when mrend=4
    real(R_), intent(in) :: fpsmth ! phase function smoothing fraction for relaxed TMS correction
    real(R_), intent(inout) :: loc(:)
    integer,  intent(inout) :: ix, iy, iz
    real(R_), intent(out) :: tau ! optical thickness
    real(R_), intent(out) :: rad ! radiance
    integer  :: icase, iq, icell, ncell
    real(R_) :: path, locnew(3), dpath, dtau, w, scamin, adir3

    ! Initialize
    adir3 = abs(dir(3))
    if (abs(dir(1)) > Vis_dgsray(1, iz) * adir3 .or. abs(dir(2)) > Vis_dgsray(2, iz) * adir3) then
       ncell = Vis_ngsray(iz)
    else
       ncell = IMAX_
    end if
    if (ix == 0) then
       ix = min(Vis_nx, int(loc(1) * Vis_facx) + 1)
       iy = min(Vis_ny, int(loc(2) * Vis_facy) + 1)
    end if

    ! Cell-by-cell integration (exact)
    do icell = 1, ncell ! loop for cells
       call mcarVis__rt_path3D(loc, ix, iy, iz, dir, itr, icase, path)
       dtau = path * Vis_dens(ix, iy, iz)
       if (mrend >= 1) then
          w = mtab_expNXC_R(dtau * Vis_fatten) * mtab_expNX_R(tau * Vis_fatten)
          rad = rad + w * Vis_srcj(ix, iy, iz)
          if (mrend == 4 .and. Vis_sdir(3) < 0.0_R_) then ! TMS correction for SR(1)
             scamin = Vis_dens(ix, iy, iz) * REPS_
             w = w * mcarVis__angCoef(ix, iy, iz, csca, qsca, scamin, fpsmth) ! a factor
             rad = rad + w * Vis_acor(ix, iy, iz) ! add corrected SR(1) source function
          end if
       end if
       tau = tau + dtau
       call mcarVis__rt_bound3D(icase, path, dir, itr, loc, ix, iy, iz)
       if (icase == 3) return ! escape from this layer
    end do

    ! Gaussian quadrature (for nearly-horizontal transfer)
    path = (Vis_zgrd(iz - 1 + itr(3)) - loc(3)) / dir(3)
    do iq = 1, Vis_ngsray(iz) ! loop for quadrature points
       dpath = path * Vis_xgsray(iq, iz)
       locnew(1:2) = loc(1:2)
       call mcarVis__rt_renewLocXY(dir, dpath, locnew) ! location
       ix = min(Vis_nx, int(locnew(1) * Vis_facx) + 1) ! cell index
       iy = min(Vis_ny, int(locnew(2) * Vis_facy) + 1)
       dtau = Vis_dens(ix, iy, iz) * path * Vis_wgsray(iq, iz) ! tau segment
       if (mrend >= 1) then
          w = mtab_expNXC_R(dtau * Vis_fatten) * mtab_expNX_R(tau * Vis_fatten)
          rad = rad + w * Vis_srcj(ix, iy, iz) ! add radiance contribution
          if (mrend == 4 .and. Vis_sdir(3) < 0.0_R_) then ! TMS correction for SR(1)
             scamin = Vis_dens(ix, iy, iz) * REPS_
             w = w * mcarVis__angCoef(ix, iy, iz, csca, qsca, scamin, fpsmth) ! a factor
             rad = rad + w * Vis_acor(ix, iy, iz) ! add corrected SR(1) source function
          end if
       end if
       tau = tau + dtau
    end do
    call mcarVis__rt_bound1D(itr(3), iz, loc(3)) ! update loc(3)
    call mcarVis__rt_renewLocXY(dir, path, loc) ! update loc(1:2)
    ix = 0 ! cell index is unknown now

  end subroutine mcarVis__trace


  !+
  ! Angular scattering coefficient (Bsca * Psca), non-scaled original value
  !-
  function mcarVis__angCoef(ix, iy, iz, csca, qsca, scamin, fpsmth) result(res)

    integer,  intent(in) :: ix, iy, iz
    real(R_), intent(in) :: csca, qsca
    real(R_), intent(in) :: scamin
    real(R_), intent(in) :: fpsmth ! artificial smoothing factor
    real(R_) :: res, sca, apf, rpf, phs
    integer :: ippo, ipf

    res = 0.0_R_
    do ippo = 1, Atm_nppo(iz)
       sca = Atm_scap3d(ix, iy, iz, ippo)
       if (sca > scamin) then ! to check active components only
          apf = Atm_apfp3d(ix, iy, iz, ippo)
          ipf = int(apf)
          if (ipf >= 1) then ! tabulated phase function
             rpf = apf - ipf ! index interpolation factor
             phs = mcarSca__phaseFunc_R(qsca, ipf, rpf, 1)
          else if (ipf == 0) then ! H-G
             phs = scatPF_HG_R(apf, csca)
          else if (ipf == -1) then ! Rayleigh
             phs = rayl_scatPF_R(csca)
          else ! isotropic
             phs = 1.0_R_
          end if
          res = res + sca * (fpsmth + (1.0_R_ - fpsmth) * phs)
       end if
    end do

  end function mcarVis__angCoef


  !+
  ! Motion to the 1-D layer boundary : modify Z
  !-
  subroutine mcarVis__rt_bound1D(itr, iz, z) 

    integer,  intent(in) :: itr ! 0 for downward, 1 for upward
    integer,  intent(inout) :: iz
    real(R_), intent(out) :: z

    z = Vis_zgrd(iz - 1 + itr)
    if (itr == 0) then ! downward
       iz = iz - 1
    else               ! upward
       iz = iz + 1
    end if

  end subroutine mcarVis__rt_bound1D


  !+
  ! Renew the X & Y location
  !-
  subroutine mcarVis__rt_renewLocXY(dir, path, loc) 

    real(R_), intent(in) :: dir(:)
    real(R_), intent(in) :: path
    real(R_), intent(inout) :: loc(:)

    loc(1:2) = loc(1:2) + path * dir(1:2) ! move
    if (loc(1) < 0.0_R_ .or. loc(1) >= Vis_xmax) then ! cycle the X location
       loc(1) = mod(loc(1), Vis_xmax)
       if (loc(1) < 0.0_R_) loc(1) = loc(1) + Vis_xmax
    end if
    if (loc(2) < 0.0_R_ .or. loc(2) >= Vis_ymax) then ! cycle the Y location
       loc(2) = mod(loc(2), Vis_ymax)
       if (loc(2) < 0.0_R_) loc(2) = loc(2) + Vis_ymax
    end if

  end subroutine mcarVis__rt_renewLocXY


  !+
  ! Move to the boundary of the box
  !-
  subroutine mcarVis__rt_bound3D(icase, path, dir, itr, loc, ix, iy, iz)

    integer,  intent(in) :: icase
    real(R_), intent(in) :: path
    real(R_), intent(in) :: dir(:)
    integer, intent(in) :: itr(:)
    integer,  intent(inout) :: ix, iy, iz
    real(R_), intent(inout) :: loc(:)

    ! Z
    if (icase == 3) then
       loc(1:2) = loc(1:2) + path * dir(1:2) ! x,y
       if (itr(3) == 0) then
          loc(3) = Vis_zgrd(iz - 1)
          iz = iz - 1
       else
          loc(3) = Vis_zgrd(iz)
          iz = iz + 1
       end if

       ! Y
    else if (icase == 2) then
       loc(1) = loc(1) + path * dir(1) ! x
       loc(3) = loc(3) + path * dir(3) ! z
       if (itr(2) == 0) then
          loc(2) = Vis_ygrd(iy - 1)
          iy = iy - 1
          if (iy <= 0) then
             loc(2) = Vis_ymax
             iy = Vis_ny
          end if
       else
          loc(2) = Vis_ygrd(iy)
          iy = iy + 1
          if (iy > Vis_ny) then
             loc(2) = 0.0_R_
             iy = 1
          end if
       end if

       ! X
    else
       loc(2:3) = loc(2:3) + path * dir(2:3) ! y,z
       if (itr(1) == 0) then
          loc(1) = Vis_xgrd(ix - 1)
          ix = ix - 1
          if (ix <= 0) then
             loc(1) = Vis_xmax
             ix = Vis_nx
          end if
       else
          loc(1) = Vis_xgrd(ix)
          ix = ix + 1
          if (ix > Vis_nx) then
             loc(1) = 0.0_R_
             ix = 1
          end if
       end if
    end if

  end subroutine mcarVis__rt_bound3D


  !+
  ! Path to a boundary of the box
  !-
  subroutine mcarVis__rt_path3D(loc, ix, iy, iz, dir, itr, icase, path)

    real(R_), intent(in) :: loc(:) ! location vector
    real(R_), intent(in) :: dir(:) ! direction vector
    integer,  intent(in) :: itr(:) ! even parity flags
    integer,  intent(in) :: ix
    integer,  intent(in) :: iy
    integer,  intent(in) :: iz
    integer,  intent(out) :: icase
    real(R_), intent(out) :: path
    real(R_) :: dc(3)

    ! Find the closest plane
    dc(3) = Vis_zgrd(iz - 1 + itr(3)) - loc(3)
    dc(2) = Vis_ygrd(iy - 1 + itr(2)) - loc(2)
    dc(1) = Vis_xgrd(ix - 1 + itr(1)) - loc(1)
    icase = 3
    if (abs(dc(2) * dir(icase)) < abs(dc(icase) * dir(2))) icase = 2
    if (abs(dc(1) * dir(icase)) < abs(dc(icase) * dir(1))) icase = 1

    ! Pathlength
    path = dc(icase) / dir(icase) ! should be [0,+Inf), not be +Inf

  end subroutine mcarVis__rt_path3D

end module mcarVis
