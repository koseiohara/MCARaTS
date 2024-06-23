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
! Module for radiance samplers
!-
module mcarRad 

  use globals
#if UseMPI == 1
  use globals_MPI
#endif
  use hparx
  use mcarUtl
  use mcarAtm
  use mcarSfc
  use mcarPho  ! Pho_* & PhoV_* are readable & writable
  use mcarVis
  implicit none
  private

  ! Public
  public :: mcarRad__user_init
  public :: mcarRad__init
  public :: mcarRad__final
  public :: mcarRad__tune_jobs
  public :: mcarRad__user_job
  public :: mcarRad__prep_job
  public :: mcarRad__qry_dims
  public :: mcarRad__qry_flags
  public :: mcarRad__extMin
  public :: mcarRad__zero
  public :: mcarRad__sampSrc
  public :: mcarRad__samp
  public :: mcarRad__samp3
  public :: mcarRad__render
#if UseMPI == 1
  public :: mcarRad__reduce
#endif
  public :: mcarRad__normal
  public :: mcarRad__writeBin
  public :: mcarRad__writeTxt
  public :: mcarRad__writeCtl

  ! (Private) User variables in namelist (1) for initialization
  integer,   save :: Rad_mrkind = 2 ! a kind of radiance
  !// = 0 : no radiance calculation
  !     1 : 1st kind, local radiance averaged over solid angle
  !     2 : 2nd kind, averaged over pixel (horizontal cross section of atmospheric column)
  !     3 : 3rd kind, averaged over solid angle and horizontal plane
  integer,   save :: Rad_mpmap = 1 ! method for pixel mapping
  !// = 1 : rectangular (U = theta * cos(phi), V = theta * sin(phi))
  !     2 : polar (U = theta, V = phi)
  integer,   save :: Rad_mplen = 0 ! method of calculation of pathlength statistics
  !// = 0 : (npr = 0)  no calculation of pathlength statistics
  !     1 : (npr = nz) layer-by-layer pathlength distribution
  !     2 : (npr = nwf) average pathlengths with weighting functions
  !     3 : (npr = ntp) histogram of total, integrated pathelength
  integer,  save :: Rad_nrad = 0  ! # of radiances
  integer,  save :: Rad_nxr = 1   ! # of X grid points
  integer,  save :: Rad_nyr = 1   ! # of Y grid points
  integer,  save :: Rad_nwf = 1   ! # of weighting functions
  integer,  save :: Rad_ntp = 100 ! # of total pathlength bins
  real(R_), save :: Rad_tpmin = 0.0_R_    ! min of total pathlength
  real(R_), save :: Rad_tpmax = 1.0e+5_R_ ! max of total pathlength
  namelist /mcarRad_nml_init/Rad_mrkind, Rad_mpmap, Rad_mplen, Rad_nrad, Rad_nxr, &
       & Rad_nyr, Rad_nwf, Rad_ntp, Rad_tpmin, Rad_tpmax

  ! (Private) User variable in namelist (2) for core calculations
  integer, parameter :: KNZ   = 1000
  integer, parameter :: KNWF  = 30
  integer, parameter :: KNRAD = 3000
  integer, private :: i
  integer,  save :: Rad_mrproj = 0 ! flag for angular weighting (for mrkind = 1 or 3)
  !// = 0 : w = 1, results are radiances simply averaged over solid angle
  !   = 1 : w = cosQ for Q = angle from the camera center direction, results are weighted average
  !   Eamples: When FOV = hemisphere (nxr = 1, nyr = 1, mpmap = 2, umax = 90 deg.),
  !    mrproj = 0 for hemispherical-mean radiance (actinic flux density)
  !    mrproj = 1 for irradiance (flux density)
  real(R_), save :: Rad_difr0 = 30.0_R_  ! numerical diffusion parameter
  real(R_), save :: Rad_difr1 = 0.01_R_  ! numerical diffusion parameter
  real(R_), save :: Rad_frmod = 0.3_R_   ! factor for mode distance (from camera)
  real(R_), save :: Rad_zetamin = 0.3_R_ ! threshold for radiance contribution function
  real(R_), save :: Rad_wfunc0(KNZ, KNWF) ! weighting functions used when Rad_mplen=2
  real(R_), save :: Rad_rmin0(KNRAD) = (/(1.0e-17_R_, i = 1, KNRAD)/) ! min distance (from camera)
  real(R_), save :: Rad_rmax0(KNRAD) = (/(1.0e+17_R_, i = 1, KNRAD)/) ! max distance (from camera)
  real(R_), save :: Rad_phi(KNRAD) = (/(0.0_R_,   i = 1, KNRAD)/) ! phi,   angle around Z0
  real(R_), save :: Rad_the(KNRAD) = (/(0.0_R_,   i = 1, KNRAD)/) ! theta, angle around Y1
  real(R_), save :: Rad_psi(KNRAD) = (/(270.0_R_, i = 1, KNRAD)/) ! psi,   angle around Z2
  !// Camera coordinates: Z-Y-Z, three rotations
  !    1: rotation about Z0 (original Z-axis in world coordinates) by phi
  !    2: rotation about Y1 (Y-axis in (X1,Y1,Z1) coordinates) by theta
  !    3: rotation about Z2 (Z-axis in (X2,Y2,Z2) coordinates) by psi
  real(R_), save :: Rad_umax(KNRAD) = (/(180.0_R_, i = 1, KNRAD)/) ! max angle along U-direction
  real(R_), save :: Rad_vmax(KNRAD) = (/(180.0_R_, i = 1, KNRAD)/) ! max angle along V-direction
  real(R_), save :: Rad_qmax(KNRAD) = (/(180.0_R_, i = 1, KNRAD)/) ! max angle of FOV cone
  real(R_), save :: Rad_xpos(KNRAD) = (/(  0.5_R_, i = 1, KNRAD)/) ! X relative position
  real(R_), save :: Rad_ypos(KNRAD) = (/(  0.5_R_, i = 1, KNRAD)/) ! Y relative position
  real(R_), save :: Rad_zloc(KNRAD) = (/(  0.0_R_, i = 1, KNRAD)/) ! Z location
  real(R_), save :: Rad_apsize(KNRAD) = (/(0.0_R_, i = 1, KNRAD)/) ! aperture size
  real(R_), save :: Rad_zref(KNRAD) = (/(  0.0_R_, i = 1, KNRAD)/) ! Z location of the reference level height
  namelist /mcarRad_nml_job/Rad_mrproj, Rad_difr0, Rad_difr1, Rad_zetamin, Rad_frmod, Rad_xpos, Rad_ypos, &
       & Rad_wfunc0, Rad_zloc, Rad_rmin0, Rad_rmax0, Rad_the, Rad_phi, Rad_psi, Rad_umax, &
       & Rad_vmax, Rad_qmax, Rad_apsize, Rad_zref

  ! Private : Radiance samplers
  integer,   parameter :: Rad_NBSIZE = 200000 ! fixed memory size for domain index table
  integer,   save :: Rad_nz     ! should be = nz
  integer,   save :: Rad_npr    ! # of grid points for pathlength statistics
  real(R_),  save :: Rad_tpfac  ! width of total pathlength grids
  real(R_),  save :: Rad_adfmin !
  real(R_),  save :: Rad_anxr, Rad_anyr
  real(R_),  save :: Rad_facx, Rad_facy
  real(RD_), save :: Rad_psrc ! # of source photons used
  real(RD_), save :: Rad_esrc ! sampler for average total energy
  integer,   save, allocatable :: Rad_mznode(:)   !(0:nz)
  integer,   save, allocatable :: Rad_iizn(:)     !(nrad)
  real(R_),  save, allocatable :: Rad_wfunc(:,:)  !(nz,nwf)
  integer,   save, allocatable :: Rad_ndxy(:)     !(nrad) # of domains in FOV
  integer,   save, allocatable :: Rad_iidx(:,:)   !(kndxy, nrad) domain index for X
  integer,   save, allocatable :: Rad_iidy(:,:)   !(kndxy, nrad) domain index for Y
  real(R_),  save, allocatable :: Rad_cloc(:,:)   !(3,nrad)
  real(R_),  save, allocatable :: Rad_crot(:,:,:) !(3,3,nrad) world-to-camera rotation matrix
  real(R_),  save, allocatable :: Rad_cdir(:,:)   !(3,nrad)
  real(R_),  save, allocatable :: Rad_pixp(:,:)   !(6,nrad)
  real(R_),  save, allocatable :: Rad_disp(:,:)   !(4,nrad)
  real(RD_), save, allocatable :: Rad_smpr(:,:,:) !(nxr,nyr,nrad)
  real(RD_), save, allocatable :: Rad_smpp(:,:,:,:) !(nxr,nyr,npr,nrad)
  real(RD_), save, allocatable :: Rad_smpr1(:)    !(nrad)
  real(RD_), save, allocatable :: Rad_smpp1(:,:)  !(npr,nrad)
  integer,   save :: Rad_ndset = 0 ! # of output datasets 

contains

  !+
  ! Read in namelist variables (1) for initialization
  !-
  subroutine mcarRad__user_init(iu) 

    integer, intent(in) :: iu   ! file unit index for namelist input
    integer :: ios
    read (iu, nml=mcarRad_nml_init, iostat=ios)
    call err_read(ios, iu, 'mcarRad__user_init: Invalid namelist input for mcarRad_nml_init.')

  end subroutine mcarRad__user_init

  !+
  ! Initialize this module
  !-
  subroutine mcarRad__init(nz) 

    integer,  intent(in) :: nz
    integer   :: kndxy

    ! Discard use of this module
    if (Rad_mrkind <= 0) then
       Rad_nrad = 0
       return
    end if

    ! Assign
    Rad_nz = nz
    Rad_tpfac = Rad_ntp / (Rad_tpmax - Rad_tpmin)
    if (Rad_mplen <= 0) then
       Rad_npr = 0
    else if (Rad_mplen == 1) then
       Rad_npr = nz
    else if (Rad_mplen == 2) then
       Rad_npr = Rad_nwf
    else if (Rad_mplen == 3) then
       Rad_npr = Rad_ntp
    end if
    Rad_ndset = 0 ! initialize # of datasets

    ! Allocate
    if (allocated(Rad_mznode)) call mcarRad__final()
    kndxy = Rad_NBSIZE / Rad_nrad
    allocate (Rad_mznode(0:Rad_nz), Rad_iizn(Rad_nrad))
    allocate (Rad_wfunc(Rad_nz, Rad_nwf)) !, Rad_ndom(4, Rad_nrad)
    allocate (Rad_cloc(3, Rad_nrad), Rad_cdir(3, Rad_nrad))
    allocate (Rad_crot(3, 3, Rad_nrad), Rad_pixp(6, Rad_nrad), Rad_disp(4, Rad_nrad))
    allocate (Rad_ndxy(Rad_nrad), Rad_iidx(kndxy, Rad_nrad), Rad_iidy(kndxy, Rad_nrad))
    allocate (Rad_smpr(Rad_nxr, Rad_nyr, Rad_nrad), Rad_smpr1(Rad_nrad))
    allocate (Rad_smpp(Rad_nxr, Rad_nyr, Rad_npr, Rad_nrad), Rad_smpp1(Rad_npr, Rad_nrad))

  end subroutine mcarRad__init

  !+
  ! Finalize this module
  !-
  subroutine mcarRad__final() 

    if (allocated(Rad_mznode)) then
       deallocate (Rad_mznode, Rad_iizn) !, Rad_ndom
       deallocate (Rad_wfunc, Rad_cloc, Rad_cdir, Rad_crot, Rad_pixp, Rad_disp)
       deallocate (Rad_ndxy, Rad_iidx, Rad_iidy)
       deallocate (Rad_smpr, Rad_smpr1, Rad_smpp, Rad_smpp1)
    end if

  end subroutine mcarRad__final

  !+
  ! Tune a preset for technical parameters
  !-
  subroutine mcarRad__tune_jobs(moptim) 

    integer, intent(in) :: moptim ! optimization flag [-2,3]

    if (moptim <= -2) then ! no preset
       return
    else if (moptim == -1) then ! no optimization
       call mcarRad__set_tech(0.0_R_, 0.0_R_, 0.0_R_)
    else if (moptim == 0) then ! unbiased optimization
       call mcarRad__set_tech(0.0_R_, 0.0_R_, 0.3_R_)
    else if (moptim == 1) then ! conservative
       call mcarRad__set_tech(3.0_R_, 0.001_R_, 0.3_R_)
    else if (moptim == 2) then ! standard
       call mcarRad__set_tech(30.0_R_, 0.01_R_, 0.3_R_)
    else ! quick-and-dirty
       call mcarRad__set_tech(100.0_R_, 0.03_R_, 0.3_R_)
    end if

  end subroutine mcarRad__tune_jobs
 
  !+
  ! Set technical parameters
  !-
  subroutine mcarRad__set_tech(difr0, difr1, zetamin) 

    real(R_), intent(in) :: difr0, difr1, zetamin
    Rad_difr0 = difr0
    Rad_difr1 = difr1
    Rad_zetamin = zetamin

  end subroutine mcarRad__set_tech

  !+
  ! Read in namelist variables (2) for core calculations
  !-
  subroutine mcarRad__user_job(iu) 

    integer, intent(in) :: iu   ! file unit index for namelist input
    integer :: ios

    Rad_wfunc0(:,:) = 1.0_R_
    read (iu, nml=mcarRad_nml_job, iostat=ios)
    call err_read(ios, iu, 'mcarRad__user_job: Invalid namelist input for mcarRad_nml_job.')
    call mcarUsr__core_check()

  end subroutine mcarRad__user_job

  !+
  ! Check user variables (2) for core calculations
  !-
  subroutine mcarUsr__core_check() 

    integer :: irad
    call check_iR('mcarRad__check: Rad_zetamin', Rad_zetamin, 0.0_R_)
    call check_iR('mcarRad__check: Rad_frmod',   Rad_frmod,   0.0_R_,  1.0_R_)
    do irad = 1, Rad_nrad
       call check_iR('mcarRad__check: Rad_umax',  Rad_umax(irad),  0.1_R_, 180.0_R_)
       call check_iR('mcarRad__check: Rad_vmax',  Rad_vmax(irad),  0.1_R_, 180.0_R_)
       call check_iR('mcarRad__check: Rad_qmax',  Rad_qmax(irad),  0.1_R_, 180.0_R_)
       call check_iR('mcarRad__check: Rad_the',   Rad_the(irad),   0.0_R_, 180.0_R_)
       call check_iR('mcarRad__check: Rad_rmin0',  Rad_rmin0(irad),  0.0_R_)
       call check_iR('mcarRad__check: Rad_rmax0',  Rad_rmax0(irad),  Rad_rmin0(irad))
    end do

  end subroutine mcarUsr__core_check

  !+
  ! Prepare core calculations
  !-
  subroutine mcarRad__prep_job() 

    real(R_), parameter :: EPSU = REPS_
    real(R_), parameter :: BUZMIN = 0.1_R_
    real(R_) :: z, rmin1, rmax1, auc, avc, amc, athe, aphi, apsi
    integer  :: irad, izn

    ! Initialize
    Rad_adfmin = max(RSML_, Rad_zetamin / PI_)
    Rad_anxr = min(15.0_R_, 0.5_R_ * real(Rad_nxr, R_))
    Rad_anyr = min(15.0_R_, 0.5_R_ * real(Rad_nyr, R_))
    Rad_facx = real(Rad_nxr, R_) / Atm_xmax * 0.9999995_R_
    Rad_facy = real(Rad_nyr, R_) / Atm_ymax * 0.9999995_R_
    !// Note: A factor slightly smaller than 1 is needed to avoid ixr(iyr) > nxr(nyr)
    Rad_wfunc(1:Rad_nz, 1:Rad_nwf) = Rad_wfunc0(1:Rad_nz, 1:Rad_nwf)

    ! For radiances of the 3rd kind
    Rad_mznode(:) = 0
    Rad_iizn(:) = -1
    if (Rad_mrkind == 3) then
       do irad = 1, Rad_nrad
          z = Rad_zloc(irad)
          if (z >= (Atm_zgrd(Atm_nz - 1) + Atm_zdep(Atm_nz)) * 0.5_R_) then
             izn = Atm_nz
          else
             do izn = 0, Rad_nz - 1
                if (z < (Atm_zgrd(izn) + Atm_zgrd(izn + 1)) * 0.5_R_) exit
             end do
          end if
          Rad_iizn(irad) = izn ! table of layer interface index
          Rad_mznode(izn) = 1 ! set on the detector on the layer interface
       end do
    end if

    ! Setup for every cameras
    do irad = 1, Rad_nrad

       ! Direction
       aphi = Rad_phi(irad) * DTOR_
       !athe = max(0.01_R_, min(179.99_R_, Rad_the(irad))) * DTOR_
       athe = max(0.0001_R_, min(179.9999_R_, Rad_the(irad))) * DTOR_ ! 06/25/2014, fixed by HI
       apsi = Rad_psi(irad) * DTOR_
       !// Note: Cameras should not be pointed at nadir/zenith, because azimuth can
       !   not be specified by x, y, z-components of a unit vector in such a case.
       Rad_cdir(:, irad) = geo3D_aVec_1R(1.0_R_, athe, aphi)
       Rad_cdir(3, irad) = nonZero_R(Rad_cdir(3, irad), EPSU) ! to avoid horizontal directions
       Rad_crot(:, :, irad) = geo3D_rotMat_ZYZ_2R(-aphi, -athe, -apsi) ! world-to-camera

       ! Location
       Rad_cloc(1, irad) = Rad_xpos(irad) * Atm_xmax
       Rad_cloc(2, irad) = Rad_ypos(irad) * Atm_ymax
       Rad_cloc(3, irad) = Rad_zloc(irad)
       call mcarUtl__horiShift(Rad_cloc(1, irad), 1.0_R_, 0.0_R_, Atm_xmax) ! cyclic condition only
       call mcarUtl__horiShift(Rad_cloc(2, irad), 1.0_R_, 0.0_R_, Atm_ymax) ! in the current codes

       ! Pixel parameters
       auc = Rad_umax(irad) * DTOR_ * 0.5_R_
       avc = Rad_vmax(irad) * DTOR_ * 0.5_R_
       amc = min(179.97_R_, Rad_qmax(irad)) * DTOR_ * 0.5_R_
       if (Rad_mpmap == 1) then ! rectangular
          amc = min(amc, sqrt(auc**2 + avc**2))
          Rad_pixp(4, irad) = real(Rad_nxr, R_) / auc * 0.5_R_ * 0.9999995_R_
          Rad_pixp(5, irad) = real(Rad_nyr, R_) / avc * 0.5_R_ * 0.9999995_R_
       else ! polar
          amc = min(amc, auc)
          Rad_pixp(4, irad) = real(Rad_nxr, R_) / auc  * 0.9999995_R_
          Rad_pixp(5, irad) = real(Rad_nyr, R_) / PI2_ * 0.9999995_R_
       end if !// Note: A factor slightly smaller than 1 is needed to avoid ixr(iyr) > nxr(nyr)
       Rad_pixp(1, irad) = cos(amc)
       Rad_pixp(2, irad) = auc
       Rad_pixp(3, irad) = avc
       Rad_pixp(6, irad) = auc * avc / real(Rad_nxr * Rad_nyr, R_) ! 1/4 (radian^2) of a pixel

       ! View region
       rmin1 = Rad_rmin0(irad)
       rmax1 = Rad_rmax0(irad)
       if (Rad_mrkind == 1) then
          call mcarRad__prep_viewRegion(Rad_cloc(:, irad), Rad_cdir(:, irad), Rad_apsize(irad), &
               & sin(amc), cos(amc), BUZMIN, Atm_xmax, Atm_ymax, Atm_zgrd(0), Atm_zgrd(Atm_nz), &
               & rmin1, rmax1, Rad_ndxy(irad), Rad_iidx(:, irad), Rad_iidy(:, irad))
       else
          Rad_ndxy(irad) = 0
       end if
       Rad_disp(1, irad) = rmin1
       Rad_disp(2, irad) = ((1.0_R_ - Rad_frmod) * sqrt(rmin1) + Rad_frmod * sqrt(rmax1))**2
       Rad_disp(3, irad) = rmax1
       Rad_disp(4, irad) = Rad_apsize(irad) * 0.5_R_ ! lens radius
    end do

  end subroutine mcarRad__prep_job

  !+
  ! Query dimension sizes
  !-
  subroutine mcarRad__qry_dims(nxr, nyr, npr, nrad) 

    integer, intent(out), optional :: nxr, nyr
    integer, intent(out), optional :: npr
    integer, intent(out), optional :: nrad
    if (present(nxr))  nxr  = Rad_nxr
    if (present(nyr))  nyr  = Rad_nyr
    if (present(npr))  npr  = Rad_npr
    if (present(nrad)) nrad = Rad_nrad

  end subroutine mcarRad__qry_dims

  !+
  ! Query about calculation flags
  !-
  subroutine mcarRad__qry_flags(mrkind, mpmap, mplen) 

    integer, intent(out), optional :: mrkind, mpmap, mplen
    if (present(mrkind)) mrkind = Rad_mrkind
    if (present(mpmap )) mpmap  = Rad_mpmap
    if (present(mplen )) mplen  = Rad_mplen

  end subroutine mcarRad__qry_flags


  !+
  ! Prepare a view region for a camera
  !-
  subroutine mcarRad__prep_viewRegion(cloc, cdir, asize, sinqmax, cosqmax, buzmin, xmax, ymax, &
       & zmin, zmax, rmin, rmax, ndxy, iidx, iidy)

    real(R_), intent(in) :: cloc(:), cdir(:)
    real(R_), intent(in) :: asize
    real(R_), intent(in) :: sinqmax, cosqmax
    real(R_), intent(in) :: buzmin
    real(R_), intent(in) :: xmax, ymax
    real(R_), intent(in) :: zmin, zmax
    real(R_), intent(inout) :: rmin, rmax
    integer,  intent(out) :: ndxy
    integer,  intent(out) :: iidx(:), iidy(:)
    real(R_) :: oloc(3), doc, doa, pmin, pmax, apmin, opmax
    real(R_) :: xvmin, xvmax, yvmin, yvmax
    integer :: i, idx, idy, idxs, idxe, idys, idye
    integer, parameter :: NDMAX = 1000 ! tentative

    ! Initialize
    doa = asize * 0.5_R_ / sinqmax  ! distance O-Amax
    doc = doa * cosqmax ! distance O-C
    oloc(:) = cloc(:) - doc * cdir(:) ! origin of the view region

    ! Modified Rmin & Rmax
    call mcarRad__prep_PMinMax(oloc(3), cdir, sinqmax, cosqmax, doa, buzmin, zmin, zmax, apmin, opmax)
    rmin = max(rmin, apmin) ! modify min distance A-P (input rmin could be 0)
    rmax = min(rmax, opmax - doc) ! modify max distance A-P (input rmax could be very large)
    pmin = rmin + doc ! O to rmin-sphere
    pmax = rmax + doa ! O to rmax-sphere
    !// Note the difference of doc/doa is due to avoid the danger to clip too large region.
    !print *, 'pmin, pmax, rmin, rmax', pmin, pmax, rmin, rmax

    ! Min & max of X/Y coordinates
    call mcarRad__prep_viewXY(oloc, cdir, sinqmax, cosqmax, zmin, zmax, pmin, pmax, &
         & xvmin, xvmax, yvmin, yvmax)
    idxs = int(xvmin / xmax) !- 1
    idxe = int(xvmax / xmax) !+ 1
    idys = int(yvmin / ymax) !- 1
    idye = int(yvmax / ymax) !+ 1
    if (xvmin < 0.0_R_) idxs = idxs - 1
    if (xvmax < 0.0_R_) idxe = idxe - 1
    if (yvmin < 0.0_R_) idys = idys - 1
    if (yvmax < 0.0_R_) idye = idye - 1
    idx = (idxs + idxe) / 2
    idy = (idys + idye) / 2
    idxs = max(idx - NDMAX, idxs)
    idxe = max(idx + NDMAX, idxe)
    idys = max(idy - NDMAX, idys)
    idye = max(idy + NDMAX, idye)
    ndxy = (idxe - idxs + 1) * (idye - idys + 1)
    !print *, 'xvmin, xvmax, yvmin, yvmax: ', xvmin, xvmax, yvmin, yvmax
    !print *, 'idxs, idxe, idys, idye: ', idxs, idxe, idys, idye

    ! Check intersection of the FOV and each column
    ndxy = 0
    loop_y: do idy = idys, idye
       yvmin = ymax * idy
       yvmax = ymax * (idy + 1)
       loop_x: do idx = idxs, idxe
          xvmin = xmax * idx
          xvmax = xmax * (idx + 1)
          i = mcarRad__prep_viewBox_I(oloc, cdir, cosqmax, pmin, pmax, xvmin, xvmax, &
               & yvmin, yvmax, zmin, zmax)
          if (i >= 1) then
             ndxy = ndxy + 1
             if (ndxy > size(iidx)) then
                print *, 'mcarRad__prep_view: Too many duplicated domains in the view region.'
                exit loop_y
             end if
             iidx(ndxy) = idx
             iidy(ndxy) = idy
          end if
       end do loop_x
    end do loop_y
    !print *, 'view domain fraction = ', real(ndxy) / real((idxe - idxs + 1) * (idye - idys + 1))

  end subroutine mcarRad__prep_viewRegion


  !+
  ! Automatically set max view distance
  !-
  subroutine mcarRad__prep_PMinMax(oriz, cdir, sinqmax, cosqmax, doa, buzmin, zmin, zmax, &
       & apmin, opmax)

    real(R_), intent(in) :: oriz
    real(R_), intent(in) :: cdir(:)
    real(R_), intent(in) :: sinqmax, cosqmax
    real(R_), intent(in) :: doa
    real(R_), intent(in) :: buzmin
    real(R_), intent(in) :: zmin, zmax
    real(R_), intent(out) :: apmin, opmax ! min & max of OP distance
    real(R_) :: vzc, uz1, uz2, uzz, path1, path2, path, sdir(3)

    ! The most oblique direction
    vzc = sqrt(cdir(1)**2 + cdir(2)**2)
    uz1 = cosqmax * cdir(3) - sinqmax * vzc ! cosQ for edge direction 1
    uz2 = cosqmax * cdir(3) + sinqmax * vzc ! cosQ for edge direction 2
    apmin = RLRG_
    opmax = RSML_

    ! Max allowed distance
    sdir(1:2) = cdir(1:2) / vzc * sqrt(1.0_R_ - buzmin**2)
    sdir(3) = buzmin
    if (geo3D_scaProd_R(cdir, sdir) > cosqmax) then ! nearly-horizontal direction in FOV
       path1 = (zmax - oriz) / buzmin
       path2 = (zmin - oriz) / buzmin
       opmax = max(opmax, path1, path2)
    end if
    !print *, 'opmax: ', opmax, path1, path2
    sdir(3) = -buzmin
    if (geo3D_scaProd_R(cdir, sdir) > cosqmax) then ! nearly-horizontal direction in FOV
       path1 = -(zmax - oriz) / buzmin
       path2 = -(zmin - oriz) / buzmin
       opmax = max(opmax, path1, path2)
    end if
    !print *, 'opmax: ', opmax, path1, path2

    ! Pathlength at cone edge directions
    path1 = (zmin - oriz) / uz1 ! distance O-Zmin1
    path2 = (zmax - oriz) / uz1 ! distance O-Zmax1
    path = max(path1, path2)
    if (path1 >= doa) apmin = min(apmin, path1 - doa)
    if (path2 >= doa) apmin = min(apmin, path2 - doa)
    if ((path1 - doa) * (path2 - doa) <= 0.0_R_) apmin = 0.0_R_ ! the aperture is in the system
    !print *, 'path1, doa, path2, apmin = ', path1, doa, path2, apmin

    path1 = (zmin - oriz) / uz2 ! distance O-Zmin2
    path2 = (zmax - oriz) / uz2 ! distance O-Zmax2
    path = max(path, path1, path2)
    if (path1 >= doa) apmin = min(apmin, path1 - doa)
    if (path2 >= doa) apmin = min(apmin, path2 - doa)
    if ((path1 - doa) * (path2 - doa) <= 0.0_R_) apmin = 0.0_R_ ! the aperture is in the system
    !print *, 'path1, doa, path2, apmin = ', path1, doa, path2, apmin

    ! Final max & min OP
    opmax = max(0.0_R_, opmax, path) ! max distance O-P
    !// opmax = 0 if the camera is completely directed to the space or underground
    if (abs(cdir(3)) >= cosqmax) then ! FOV includes the zenith/nadir
       uzz = sign(1.0_R_, cdir(3))
       path1 = (zmin - oriz) * uzz ! distance O-Zmin0
       path2 = (zmax - oriz) * uzz ! distance O-Zmax0
       if (path1 >= doa) apmin = min(apmin, path1 - doa)
       if (path2 >= doa) apmin = min(apmin, path2 - doa)
       if ((path1 - doa) * (path2 - doa) <= 0.0_R_) apmin = 0.0_R_ ! the aperture is in the system
    end if
    !print *, 'Final, opmax, apmin = ', opmax, apmin

  end subroutine mcarRad__prep_PMinMax


  !+
  ! Calculate min & max of X/Y coordinates in the view region
  !-
  subroutine mcarRad__prep_viewXY(oloc, cdir, sinqmax, cosqmax, zmin, zmax, pmin, pmax, &
       & xvmin, xvmax, yvmin, yvmax)

    real(R_), intent(in) :: oloc(:), cdir(:)
    real(R_), intent(in) :: sinqmax, cosqmax
    real(R_), intent(in) :: zmin, zmax
    real(R_), intent(in) :: pmin, pmax
    real(R_), intent(out) :: xvmin, xvmax
    real(R_), intent(out) :: yvmin, yvmax
    real(R_) :: sdir(3), path1, path2
    real(R_), parameter :: EPSU = 0.1_R_

    ! X min
    sdir(:) = (/-1.0_R_, 0.0_R_, 0.0_R_/)
    sdir(:) = mcarRad__prep_coneNearDirec_1R(cdir, sinqmax, cosqmax, sdir)
    xvmin = oloc(1) + pmin * sdir(1)
    if (abs(sdir(3)) < EPSU) then
       xvmin = min(xvmin, oloc(1) + pmax * sdir(1))
    else
       path1 = (zmin - oloc(3)) / sdir(3) ! distance O-Xmin1
       path2 = (zmax - oloc(3)) / sdir(3) ! distance O-Xmin2
       if (path1 > 0.0_R_) xvmin = min(xvmin, oloc(1) + min(pmax, path1) * sdir(1))
       if (path2 > 0.0_R_) xvmin = min(xvmin, oloc(1) + min(pmax, path2) * sdir(1))
    end if

    ! X max
    sdir(:) = (/1.0_R_, 0.0_R_, 0.0_R_/)
    sdir(:) = mcarRad__prep_coneNearDirec_1R(cdir, sinqmax, cosqmax, sdir)
    xvmax = oloc(1) + pmin * sdir(1)
    if (abs(sdir(3)) < EPSU) then
       xvmax = max(xvmax, oloc(1) + pmax * sdir(1))
    else
       path1 = (zmin - oloc(3)) / sdir(3) ! distance O-Xmax1
       path2 = (zmax - oloc(3)) / sdir(3) ! distance O-Xmax2
       if (path1 > 0.0_R_) xvmax = max(xvmax, oloc(1) + min(pmax, path1) * sdir(1))
       if (path2 > 0.0_R_) xvmax = max(xvmax, oloc(1) + min(pmax, path2) * sdir(1))
    end if

    ! Y min
    sdir(:) = (/0.0_R_, -1.0_R_, 0.0_R_/)
    sdir(:) = mcarRad__prep_coneNearDirec_1R(cdir, sinqmax, cosqmax, sdir)
    yvmin = oloc(2) + pmin * sdir(2)
    if (abs(sdir(3)) < EPSU) then
       yvmin = min(yvmin, oloc(2) + pmax * sdir(2))
    else
       path1 = (zmin - oloc(3)) / sdir(3) ! distance O-Ymin1
       path2 = (zmax - oloc(3)) / sdir(3) ! distance O-Ymin2
       if (path1 > 0.0_R_) yvmin = min(yvmin, oloc(2) + min(pmax, path1) * sdir(2))
       if (path2 > 0.0_R_) yvmin = min(yvmin, oloc(2) + min(pmax, path2) * sdir(2))
    end if

    ! Y max
    sdir(:) = (/0.0_R_, 1.0_R_, 0.0_R_/)
    sdir(:) = mcarRad__prep_coneNearDirec_1R(cdir, sinqmax, cosqmax, sdir)
    yvmax = oloc(2) + pmin * sdir(2)
    if (abs(sdir(3)) < EPSU) then
       yvmax = max(yvmax, oloc(2) + pmax * sdir(2))
    else
       path1 = (zmin - oloc(3)) / sdir(3) ! distance O-Ymax1
       path2 = (zmax - oloc(3)) / sdir(3) ! distance O-Ymax2
       if (path1 > 0.0_R_) yvmax = max(yvmax, oloc(2) + min(pmax, path1) * sdir(2))
       if (path2 > 0.0_R_) yvmax = max(yvmax, oloc(2) + min(pmax, path2) * sdir(2))
    end if

  end subroutine mcarRad__prep_viewXY


  !+
  ! Roughly determine whether the box is in the view region
  !-
  function mcarRad__prep_viewBox_I(oloc, cdir, cosqmax, pmin, pmax, xmin, xmax, &
       & ymin, ymax, zmin, zmax) result(m)

    real(R_), intent(in) :: oloc(:), cdir(:)
    real(R_), intent(in) :: cosqmax
    real(R_), intent(in) :: pmin, pmax
    real(R_), intent(in) :: xmin, xmax
    real(R_), intent(in) :: ymin, ymax
    real(R_), intent(in) :: zmin, zmax
    integer  :: m ! result, 0 for no intersection, 1 or larger for possible intersection
    integer  :: i
    logical  :: hit
    real(R_) :: dloc(3, 4), vvec(3), dov, smax, path

    ! The center point of the box
    m = 0
    vvec(1) = (xmin + xmax) * 0.5_R_
    vvec(2) = (ymin + ymax) * 0.5_R_
    vvec(3) = (zmin + zmax) * 0.5_R_
    vvec(:) = vvec(:) - oloc(:) ! OV vector
    dov = geo3D_norm_R(vvec) ! OV distance
    if (geo3D_scaProd_R(cdir, vvec) >= dov * cosqmax .and. &
         & (pmin - dov) * (pmax - dov) <= 0.0_R_) then
       m = 1
       return
    end if

    ! The center direction vs the box : 6 planes
    call geo3D_rayNPlane(pmin, pmax, oloc(1), cdir(1), xmin, oloc(2), cdir(2), &
         & ymin, ymax, oloc(3), cdir(3), zmin, zmax, path, hit)
    if (hit) then
       m = 11
       return
    end if
    call geo3D_rayNPlane(pmin, pmax, oloc(1), cdir(1), xmax, oloc(2), cdir(2), &
         & ymin, ymax, oloc(3), cdir(3), zmin, zmax, path, hit)
    if (hit) then
       m = 12
       return
    end if
    call geo3D_rayNPlane(pmin, pmax, oloc(2), cdir(2), ymin, oloc(3), cdir(3), &
         & zmin, zmax, oloc(1), cdir(1), xmin, xmax, path, hit)
    if (hit) then
       m = 13
       return
    end if
    call geo3D_rayNPlane(pmin, pmax, oloc(2), cdir(2), ymax, oloc(3), cdir(3), &
         & zmin, zmax, oloc(1), cdir(1), xmin, xmax, path, hit)
    if (hit) then
       m = 14
       return
    end if
    call geo3D_rayNPlane(pmin, pmax, oloc(3), cdir(3), zmin, oloc(1), cdir(1), &
         & xmin, xmax, oloc(2), cdir(2), ymin, ymax, path, hit)
    if (hit) then
       m = 15
       return
    end if
    call geo3D_rayNPlane(pmin, pmax, oloc(3), cdir(3), zmax, oloc(1), cdir(1), &
         & xmin, xmax, oloc(2), cdir(2), ymin, ymax, path, hit)
    if (hit) then
       m = 16
       return
    end if

    ! The cone vs box edges : 12 line segments
    smax = zmax - zmin
    vvec(:) = (/0.0_R_, 0.0_R_, 1.0_R_/)
    dloc(:, 1) = (/xmin, ymin, zmin/) - oloc(:)
    dloc(:, 2) = (/xmax, ymin, zmin/) - oloc(:)
    dloc(:, 3) = (/xmin, ymax, zmin/) - oloc(:)
    dloc(:, 4) = (/xmax, ymax, zmin/) - oloc(:)
    m = 4
    do i = 1, 4
       if (dloc(1, i)**2 + dloc(2, i)**2 > pmax**2) m = m - 1
    end do
    if (m == 0) return ! no hit
    do i = 1, 4
       if (geo3D_rayCone_L(vvec, dloc(:, i), smax, cdir, cosqmax)) then
          m = 130 + i
          return
       end if
    end do
    smax = xmax - xmin
    vvec(:) = (/1.0_R_, 0.0_R_, 0.0_R_/)
    dloc(:, 1) = (/xmin, ymin, zmin/) - oloc(:)
    dloc(:, 2) = (/xmin, ymax, zmin/) - oloc(:)
    dloc(:, 3) = (/xmin, ymin, zmax/) - oloc(:)
    dloc(:, 4) = (/xmin, ymax, zmax/) - oloc(:)
    do i = 1, 4
       if (geo3D_rayCone_L(vvec, dloc(:, i), smax, cdir, cosqmax)) then
          m = 110 + i
          return
       end if
    end do
    smax = ymax - ymin
    vvec(:) = (/0.0_R_, 1.0_R_, 0.0_R_/)
    dloc(:, 1) = (/xmin, ymin, zmin/) - oloc(:)
    dloc(:, 2) = (/xmax, ymin, zmin/) - oloc(:)
    dloc(:, 3) = (/xmin, ymin, zmax/) - oloc(:)
    dloc(:, 4) = (/xmax, ymin, zmax/) - oloc(:)
    do i = 1, 4
       if (geo3D_rayCone_L(vvec, dloc(:, i), smax, cdir, cosqmax)) then
          m = 120 + i
          return
       end if
    end do
    m = 0 ! finally, no intersection with the box

  end function mcarRad__prep_viewBox_I


  !+
  ! The nearest direction to a specific direction, within a conical region
  !-
  function mcarRad__prep_coneNearDirec_1R(cdir, sinqmax, cosqmax, adir) result(ndir)

    real(R_), intent(in) :: sinqmax, cosqmax
    real(R_), intent(in) :: cdir(:)
    real(R_), intent(in) :: adir(:)
    real(R_)  :: ndir(size(cdir)) ! the nearest direction
    real(R_)  :: cosq, sinq
    
    ! Thee cone vs the direction A
    call geo3D_twoUVec(cdir, adir, cosq, sinq)
    if (cosq >= 0.0_R_ .and. sinq < sinqmax) then ! use sinQ mainly (cosQ could be inaccurate)
       ndir(:) = adir(:)
    else
       ndir(:) = geo3D_rotate_polar_1R(cdir, adir, -sinqmax, cosqmax)
       call geo3D_twoUVec(cdir, ndir, cosq, sinq)
    end if

  end function mcarRad__prep_coneNearDirec_1R


  !+
  ! Set mininum extinction coefficient for LCF
  !-
  subroutine mcarRad__extMin(dtaumin, nlcf, rlcf, extmin)

    real(R_), intent(in)  :: dtaumin
    integer,  intent(in)  :: nlcf
    real(R_), intent(in)  :: rlcf
    real(R_), intent(out) :: extmin(:)
    real(R_)  :: dz, ext, za, zc, vec(3), eta(3), sinavc, cosavc, uzmin, uzmax
    integer   :: irad, iz, izc

    ! Initialize
    eta(:) = (/0.0_R_, 0.0_R_, 1.0_R_/)
    extmin(:) = 0.0_R_

    ! Loop for cameras
    do irad = 1, Rad_nrad
       zc = Rad_cloc(3, irad)
       izc = gridIdx_bin0_I(Atm_zgrd, zc, 1, Atm_nz)
       cosavc = Rad_pixp(1, irad)
       sinavc = sqrt(1.0_R_ - cosavc**2)
       vec(:) = geo3D_rotate_polar_1R(Rad_cdir(:, irad), eta,  sinavc, cosavc)
       uzmin = vec(3)
       uzmax = vec(3)
       vec(:) = geo3D_rotate_polar_1R(Rad_cdir(:, irad), eta, -sinavc, cosavc)
       uzmin = min(uzmin, vec(3))
       uzmax = max(uzmax, vec(3))

       ! When having upward looking pixels
       if (uzmax > 0.01_R_) then
          do iz = izc + 1, Atm_nz
             dz = Atm_zdep(iz)
             za = 0.5_R_ * (Atm_zgrd(iz - 1) + Atm_zgrd(iz))
             ext = dtaumin / dz
             if (za - zc > rlcf) then
                ext = ext * (rlcf / (za - zc))**nlcf
             end if
             extmin(iz) = max(extmin(iz), ext)
          end do
       end if
       
       ! When having downward looking pixels
       if (uzmin < -0.01_R_) then
          do iz = 1, izc
             dz = Atm_zdep(iz)
             za = 0.5_R_ * (Atm_zgrd(iz - 1) + Atm_zgrd(iz))
             ext = dtaumin / dz
             if (zc - za > rlcf) then
                ext = ext * (rlcf / (zc - za))**nlcf
             end if
             extmin(iz) = max(extmin(iz), ext)
          end do
       end if
    end do

  end subroutine mcarRad__extMin


  !+
  ! Zero samplers
  !-
  subroutine mcarRad__zero() 

    Rad_psrc = 0.0_RD_
    Rad_esrc = 0.0_RD_
    Rad_smpr1(:)      = 0.0_RD_
    Rad_smpp1(:,:)    = 0.0_RD_
    Rad_smpr(:,:,:)   = 0.0_RD_
    Rad_smpp(:,:,:,:) = 0.0_RD_

  end subroutine mcarRad__zero


  !+
  ! Sample # of source photons & energy
  !-
  subroutine mcarRad__sampSrc(np) 

    integer, intent(in) :: np
    Rad_psrc = Rad_psrc + np
    Rad_esrc = Rad_esrc + np * Pho_eone

  end subroutine mcarRad__sampSrc


  !+
  ! Sample radiance contributions by the local estimation method
  !-
  subroutine mcarRad__samp(ix, iy, imod, isub, pac)

    integer,  intent(in) :: imod ! [-1,4]
    integer,  intent(in) :: isub
    integer,  intent(in) :: ix, iy
    real(R_), intent(in) :: pac(:)

    if (Rad_mrkind == 1) then
       call mcarRad__samp1(ix, iy, imod, isub, pac)
    else if (Rad_mrkind == 2) then
       call mcarRad__samp2(ix, iy, imod, isub, pac)
    end if

  end subroutine mcarRad__samp


  !+
  ! Sample contributions for radiances of the first kind
  !// Defs:
  !   Radiance [W/m^2/str]
  !        L = sum{E/dS*P*exp(-tau)}
  !   Solid angle dSa = dS*abs(u1)/r^2, where r is distance of the camera from
  !   the source. Note this program computes dSa*L [W/m2]:
  !        dSa*L = sum{E*abs(u1)/r^2*P*exp(-tau)}
  !// S: source, point  of scattering or source emission at original domain
  !   T: target, points of scattering or source emission at duplicated domains
  !   C: camera, point  of random location within the lens surface
  !-
  subroutine mcarRad__samp1(ix, iy, imod, isub, pac)

    integer,  intent(in) :: imod
    integer,  intent(in) :: isub
    integer,  intent(in) :: ix, iy
    real(R_), intent(in) :: pac(:)
    real(R_)  :: cloc(3), t2c(3), s2c(3), cdir(3), crot(3,3), vdir(3)
    real(R_)  :: adf, auc, alpha, amc, avc, cosq, sinq, hsrc, wfac
    real(R_)  :: dsa, dsac, f, facdsa, fdif, fuc, fvc, hnum, zdest
    real(R_)  :: p, path, q, r, rapt, rmax, rmid, rmin, sfac, sina, tau, taulim, taumax
    integer   :: irad, idet, idx, idx0, idy, idy0 !, idxmax, idxmin, idymax, idymin
    integer   :: ix2, ixr, iy2, iyr, ndxy, idxy, idxy1, ndxy1 !, n
    integer   :: ndx = 1, ndy = 1
    integer, parameter :: NRND = 500

    ! Loop for cameras
    do irad = 1, Rad_nrad

       ! Initialize
       crot(:,:) = Rad_crot(:,:,irad)
       cdir(:) = Rad_cdir(:, irad)
       rmin = Rad_disp(1, irad)
       rmid = Rad_disp(2, irad)
       rmax = Rad_disp(3, irad)
       rapt = Rad_disp(4, irad)
       amc  = Rad_pixp(1, irad)
       auc  = Rad_pixp(2, irad)
       avc  = Rad_pixp(3, irad)
       fuc  = Rad_pixp(4, irad)
       fvc  = Rad_pixp(5, irad)
       dsac = Rad_pixp(6, irad)
       fdif = Rad_difr0 + Rad_difr1 * Pho_sdif
       cloc(:) = geo3D_randLoc_circle_1R(Rad_cloc(:, irad), cdir, rapt) ! a lens surface point
       !// This is possibly random.
       s2c(:) = cloc(:) - Pho_loc(:) ! source to camera vector
       if (abs(s2c(3)) < RSML_) cycle ! too close to the source location, skip all below

       ! Domain shift for solar or point source
       facdsa = 1.0_R_ ! usually 1
       if (imod <= 0) then
          path = abs(s2c(3) / Pho_dir(3))
          !if (path > rmax) then ! move the source if it is too far
          !   facdsa = (rmax / path)**2 ! scaling factor depending on the distance
          !   path = rmax
          !   s2c(3) = path * Pho_dir(3) ! the source moves to the point of rmax
          !end if
          !idx0 = int(min(1.0e+4_R_, path / Atm_xmax * Pho_dir(1)))
          !idy0 = int(min(1.0e+4_R_, path / Atm_ymax * Pho_dir(2)))
          idx0 = int(path / Atm_xmax * Pho_dir(1))
          idy0 = int(path / Atm_ymax * Pho_dir(2))
          if (Pho_itr(1) == 0) idx0 = idx0 - 1
          if (Pho_itr(2) == 0) idy0 = idy0 - 1
          s2c(1) = s2c(1) + idx0 * Atm_xmax ! shift domain
          s2c(2) = s2c(2) + idy0 * Atm_ymax
          hsrc = abs(s2c(3)) * tan(min(PIH_ * 0.9_R_, 1.3_R_ * pac(1))) / Pho_dir(3)**2
          ndx = int(hsrc / Atm_xmax + 1.0_R_) ! add tolerence of 1
          ndy = int(hsrc / Atm_ymax + 1.0_R_)
          ndxy = (2 * ndx + 1) * (2 * ndy + 1) ! # of shifted domains
       else
          ndxy = Rad_ndxy(irad)
       end if
       ndxy1 = min(ndxy, NRND) ! # of randomly chosen domains

       ! Loop for horizontally-shifted domains
       do idxy1 = 1, ndxy1

          ! Domain index
          wfac = facdsa
          idxy = idxy1
          if (ndxy1 < ndxy) then
             idxy = min(ndxy, int(mseq_rand_R() * real(ndxy, R_)) + 1)
             wfac = wfac * real(ndxy, R_) / real(ndxy1, R_) ! scale weight
          end if
          if (imod <= 0) then
             idy = (idxy - 1) / (2 * ndx + 1) + 1
             idx = idxy - (2 * ndx + 1) * (idy - 1)
             idx = idx - ndx - 1
             idy = idy - ndy - 1
          else
             idx = Rad_iidx(idxy, irad)
             idy = Rad_iidy(idxy, irad)
          end if

          ! View vector & distance
          t2c(1) = s2c(1) - real(idx, R_) * Atm_xmax
          t2c(2) = s2c(2) - real(idy, R_) * Atm_ymax
          t2c(3) = s2c(3)
          !// Note: t2c(:) corresponds to emission-to-camera direction
          r = geo3D_norm_R(t2c) ! distance
          if (r < rmin) cycle
          if (imod >= 1 .and. r > rmax) cycle
          call mcarPho__rt_newDirecV(t2c(:) / r) ! view vector (world coordinates)
          vdir(:) = matmul(crot, -PhoV_dir) ! view vector (camera coordinates)
          if (vdir(3) < amc) cycle ! out of the FOV cone
          if (abs(PhoV_dir(3)) < REPS_) cycle ! reject horizontal directions

          ! Solar/point direct beam
          if (imod <= 0) then
             call geo3D_twoUVec(Pho_dir, PhoV_dir, cosq, sinq)
             q = atan2(sinq, cosq) ! [0,pi] a better formula for small q
             if (q > pac(1)) cycle ! out of the emission cone
          end if

          ! Random cutoff of too far samples
          if (imod >= 1 .and. r > rmid) then
             f = max(0.03_R_, rmid / r)
             if (mseq_rand_R() > f) cycle ! Russian roulette
             wfac = wfac / f ! scaling due to RR
          end if

          ! Image geometry
          call mcarRad__pixel(vdir, auc, avc, fuc, fvc, Rad_nxr, Rad_nyr, sina, alpha, &
               & ixr, iyr, idet)
          if (idet /= 1) cycle

          ! Duplicate a photon packet
          call mcarPho__rt_duplicate()
          ix2 = ix
          iy2 = iy

          ! PDF of angular distribution
          adf = mcarRad__angDistr_R(imod, isub, pac, ix, iy)
          call mcarPho__rt_ichi_next(adf) ! update
          PhoV_wgt = PhoV_wgt * adf
          call mcarUtl__scaleADF(PhoV_wgt, Rad_adfmin, taulim)
          if (PhoV_wgt < RSML_) cycle

          ! Trace the virtual trajectory
          zdest = cloc(3) - PhoV_dir(3) * rmin ! move the destination point to the rmin-sphere
          tau = 0.0_R_
          taumax = taulim + rand_exp_R()
          call mcarAtm__rt_escape(tau, taumax, PhoV_plen, PhoV_loc, PhoV_dir, PhoV_itr, &
               & ix2, iy2, PhoV_iz, Pho_ikd, PhoV_ichi, zdest, PhoV_mv3D)
          if (tau >= taumax) cycle
          !// Now, the virtual photon locates on the surface of a sphere centered at C
          !   with a radius of rmin. Transmittance for r = [0, rmin] is artificially set at 1.
          tau = min(tau, taulim)

          ! Contribution function
          PhoV_wgt = PhoV_wgt * mtab_expNX_R(tau)
          dsa = wfac / r**2 ! = 1 / r^2, possibly scaled
          p = Pho_eone * PhoV_wgt * dsa
          if (imod == 3) p = p * pac(2)

          ! Integration
          if (PhoV_iso >= 2 .and. PhoV_mv3d) then ! scattering component and 3DRT
             sfac = fdif * sqrt(PI_ * PhoV_wgt)
             hnum = sfac * sqrt(Pho_hdwx * Pho_hdwy * abs(PhoV_dir(3)) &
                  & * (dsa * alpha) / (dsac * sina))
          else
             sfac = 0.0_R_
             hnum = 0.0_R_
          end if
          call mcarRad__samp1_diff(ixr, iyr, irad, p, tau, hnum, sfac, Pho_loc(1), Pho_loc(2), &
               & t2c, PhoV_dir(3), vdir(3), rmax, crot, amc, auc, avc, fuc, fvc)
       end do
    end do

  end subroutine mcarRad__samp1


  !+
  ! Integrate camera radiance
  !-
  subroutine mcarRad__samp1_diff(ixr, iyr, irad, p, tau, hnum, sfac, x, y, t2c, uze, cosq, rmax, crot, &
       & amc, auc, avc, fuc, fvc)

    real(R_), intent(in) :: auc, avc
    real(R_), intent(in) :: amc
    real(R_), intent(in) :: t2c(:)
    real(R_), intent(in) :: fuc, fvc
    real(R_), intent(inout) :: hnum
    integer,  intent(in) :: irad
    integer,  intent(inout) :: ixr, iyr
    real(R_), intent(inout) :: p
    real(R_), intent(in) :: rmax
    real(R_), intent(in) :: sfac
    real(R_), intent(in) :: tau
    real(R_), intent(in) :: crot(:,:)
    real(R_), intent(in) :: uze
    real(R_), intent(in) :: cosq
    real(R_), intent(in) :: x, y
    real(R_),  parameter :: HNUMMAX = 16.0_R_
    integer,   parameter :: KNWRK = (2*18 + 1)**2
    real(R_)  :: wgtwrk(KNWRK)
    integer   :: ixrwrk(KNWRK), iyrwrk(KNWRK)
    integer   :: idet, ixd, iyd, n, nxdhi, nxdlo, nydhi, nydlo, iprmin, iprmax, nfull, i
    real(R_)  :: alpha, fcor, pdeg, r1, rixr, riyr, rnxd1, rnyd1, sina, wgt, xrbin, yrbin
    real(R_)  :: p1, edir(3), vdir(3), pp(Rad_npr) !(AUTO)

    ! Pathlength statistics
    call mcarRad__pathLen(PhoV_plen, pp, iprmin, iprmax)

    ! No diffusion
    if (hnum < 0.5_R_) then
       if (Rad_mrproj == 1) p = p * cosq
       Rad_smpr(ixr, iyr, irad) = Rad_smpr(ixr, iyr, irad) + p
       Rad_smpp(ixr, iyr, iprmin:iprmax, irad) = Rad_smpp(ixr, iyr, iprmin:iprmax, irad) &
            & + p * pp(iprmin:iprmax)

       ! Diffusion
    else
       ! Truncate extremes
       if (hnum > HNUMMAX) then
          pdeg = HNUMMAX / hnum
          hnum = HNUMMAX
          p1 = p * (1.0_R_ - pdeg)
          if (Rad_mrproj == 1) p1 = p1 * cosq
          Rad_smpr1(irad) = Rad_smpr1(irad) + p1
          Rad_smpp1(iprmin:iprmax, irad) = Rad_smpp1(iprmin:iprmax, irad) + pp(iprmin:iprmax) * p1
          p = p * pdeg
       end if
       !// The excess should be collected and redistributed later.
       !   At present, it is truncated simply, biasing the domain average.

       ! Position parameters
       xrbin = min(0.5_R_ * Atm_xmax, sfac * Pho_hdwx) / hnum
       yrbin = min(0.5_R_ * Atm_ymax, sfac * Pho_hdwy) / hnum
       rixr = x / xrbin
       riyr = y / yrbin
       rnxd1 = real(int(rixr), R_) - rixr + 0.5_R_
       rnyd1 = real(int(riyr), R_) - riyr + 0.5_R_
       nxdlo = int(hnum + rnxd1)
       nxdhi = int(hnum - rnxd1)
       nydlo = int(hnum + rnyd1)
       nydhi = int(hnum - rnyd1)
       !call mcarUtl__getGtab2(ixlut, KND, nxdlo, nxdhi, x, xrbin, fx, nx)
       !call mcarUtl__getGtab2(iylut, KND, nydlo, nydhi, y, yrbin, fy, ny)

       ! Sample
       n = 0
       nfull = ((nxdhi + nxdlo + 1) * (nydhi + nydlo + 1))
       p = p / nfull
       do iyd = -nydlo, nydhi
          do ixd = -nxdlo, nxdhi
             edir(1) = t2c(1) - xrbin * (ixd - 0.5_R_ + mseq_rand_R()) ! stratified MC
             edir(2) = t2c(2) - yrbin * (iyd - 0.5_R_ + mseq_rand_R())
             edir(3) = t2c(3)
             r1 = geo3D_norm_R(edir)
             if (r1 > rmax .or. r1 < RSML_) cycle
             !if (Atm_mlay(iz) == 3) then ! 3-D layer
             !   ixt = ixlut(ixd)
             !   iyt = iylut(iyd)
             !   if (Atm_scat3d(ixt, iyt, iz, ichi) <= exteps) cycle
             !end if
             edir(:) = edir(:) / r1 ! view vector (world coordinates)
             vdir(:) = matmul(crot, -edir) ! view vector (camera coordinates)
             if (vdir(3) < amc) cycle ! out of the FOV cone
             call mcarRad__pixel(vdir, auc, avc, fuc, fvc, Rad_nxr, Rad_nyr, &
                  & sina, alpha, ixr, iyr, idet)
             if (idet /= 1) cycle
             fcor = uze / edir(3) ! correction due to zenith angle difference
             wgt = mtab_expX_R(tau * (1.0_R_ - fcor)) * fcor ! corection
             n = n + 1
             if (Rad_mrproj == 1) wgt = wgt * vdir(3)
             Rad_smpr(ixr, iyr, irad) = Rad_smpr(ixr, iyr, irad) + p * wgt
             Rad_smpp(ixr, iyr, iprmin:iprmax, irad) = Rad_smpp(ixr, iyr, iprmin:iprmax, irad) &
                  & + p * wgt * pp(iprmin:iprmax)
             wgtwrk(n) = wgt
             ixrwrk(n) = ixr
             iyrwrk(n) = iyr
          end do
       end do

       ! Recover trimed margins
       if (n < nfull .and. n > 0) then
          p = p * (nfull - n) / n
          do i = 1, n
             wgt = wgtwrk(i)
             ixr = ixrwrk(i)
             iyr = iyrwrk(i)
             Rad_smpr(ixr, iyr, irad) = Rad_smpr(ixr, iyr, irad) + p * wgt
             Rad_smpp(ixr, iyr, iprmin:iprmax, irad) = Rad_smpp(ixr, iyr, iprmin:iprmax, irad) &
                  & + p * wgt * pp(iprmin:iprmax)
          end do
       end if
    end if

  end subroutine mcarRad__samp1_diff


  !+
  ! Sample contributions for radiances of the second kind
  ! Def: Radiance [W/m2/str], L = sum{E/dS*P*exp(-tau)}
  !   Note this program computes cos(q)*dS*L [W/str].
  !-
  subroutine mcarRad__samp2(ix, iy, imod, isub, pac)

    integer,  intent(in) :: imod ! [-1,4]
    integer,  intent(in) :: isub
    integer,  intent(in) :: ix, iy
    real(R_), intent(in) :: pac(:)
    real(R_)  :: adf, cosq, sinq, fs, zdest
    real(R_)  :: p, q, rnx, rny, tau, taulim, taumax
    integer   :: irad, ix2, ixr, iy2, iyr

    ! Loop for directions
    do irad = 1, Rad_nrad

       ! Initialize
       call mcarPho__rt_newDirecV(-Rad_cdir(1:3, irad)) ! view vector (incoming to camera)
       zdest = Rad_cloc(3, irad) ! index of radiance sampling level
       if (PhoV_dir(3) * (zdest - Pho_loc(3)) <= 0.0_R_) cycle

       ! Solar/point source
       if (imod <= 0) then
          call geo3D_twoUVec(PhoV_dir, Pho_dir, cosq, sinq)
          q = atan2(sinq, cosq) ! [0,pi] better formula for small q
          if (q > pac(1)) cycle
       end if

       ! Duplicate a photon packet
       call mcarPho__rt_duplicate()
       ix2 = ix
       iy2 = iy
       adf = mcarRad__angDistr_R(imod, isub, pac, ix, iy)
       call mcarPho__rt_ichi_next(adf) ! update
       PhoV_wgt = PhoV_wgt * adf
       call mcarUtl__scaleADF(PhoV_wgt, Rad_adfmin, taulim)
       if (PhoV_wgt < RSML_) cycle

       ! Escaping energy
       tau = 0.0_R_
       taumax = taulim + rand_exp_R()
       call mcarAtm__rt_escape(tau, taumax, PhoV_plen, PhoV_loc, PhoV_dir, PhoV_itr, &
            & ix2, iy2, PhoV_iz, Pho_ikd, PhoV_ichi, zdest, PhoV_mv3D)
       if (tau >= taumax) cycle
       if (tau >= taulim) tau = taulim
       if (imod == 3) PhoV_wgt = PhoV_wgt * pac(2)
       PhoV_wgt = PhoV_wgt * mtab_expNX_R(tau)
       p = Pho_eone * PhoV_wgt

       ! Integration
       ixr = int(PhoV_loc(1) * Rad_facx) + 1
       iyr = int(PhoV_loc(2) * Rad_facy) + 1
       if (PhoV_iso >= 2 .and. PhoV_mv3d) then ! scattering component and 3DRT
          fs = (Rad_difr0 + Rad_difr1 * Pho_sdif) * sqrt(PI_ * PhoV_wgt)
          rnx = Rad_facx * Pho_hdwx * fs
          rny = Rad_facy * Pho_hdwy * fs
       else
          rnx = 0.0_R_
          rny = 0.0_R_
       end if
       call mcarRad__samp2_diff(ixr, iyr, irad, p, rnx, rny, PhoV_loc(1), PhoV_loc(2), &
            & Pho_loc(1), Pho_loc(2), Pho_iz, PhoV_ichi)
    end do

  end subroutine mcarRad__samp2


  !+
  ! Integrate radiance
  !-
  subroutine mcarRad__samp2_diff(ixr, iyr, irad, p, rnx, rny, x2, y2, x, y, iz, ichi)

    integer,  intent(in) :: ichi
    integer,  intent(in) :: irad
    integer,  intent(in) :: ixr, iyr
    integer,  intent(in) :: iz
    real(R_), intent(inout) :: p
    real(R_), intent(inout) :: rnx, rny
    real(R_), intent(in) :: x, y
    real(R_), intent(in) :: x2, y2
    real(R_), parameter :: RNMAX = 30.0_R_ ! max half number of pixels for diffusion
    integer,  parameter :: KND = 40, KNWRK = (2 * KND + 1)**2
    integer  :: ixlut(-KND:KND), iylut(-KND:KND), ixrlut(-KND:KND), iyrlut(-KND:KND)
    integer  :: iixr2(KNWRK), iiyr2(KNWRK)
    real(R_) :: pp(Rad_npr) !(AUTO)
    real(R_) :: dxr, dyr, pdeg, rnx1, rny1
    integer  :: ixd, ixt, iyd, iyt, ixu, iyu, n, nxdhi, nxdlo, nydhi, nydlo
    integer  :: iprmin, iprmax, i, nfull

    ! Pathlength statistics
    call mcarRad__pathLen(PhoV_plen, pp, iprmin, iprmax)
       
    ! No diffusion
    if (rnx < 0.5_R_ .and. rny < 0.5_R_) then
       Rad_smpr(ixr, iyr, irad) = Rad_smpr(ixr, iyr, irad) + p
       Rad_smpp(ixr, iyr, iprmin:iprmax, irad) = Rad_smpp(ixr, iyr, iprmin:iprmax, irad) &
            & + p * pp(iprmin:iprmax)

       ! Diffusion
    else
       ! Truncate extreme values
       pdeg = 1.0_R_
       if (rnx > RNMAX) then 
          pdeg = pdeg * RNMAX / rnx
          rnx = RNMAX
       end if
       if (rny > RNMAX) then
          pdeg = pdeg * RNMAX / rny
          rny = RNMAX
       end if
       if (pdeg < 1.0_R_) then
          Rad_smpr1(irad) = Rad_smpr1(irad) + p * (1.0_R_ - pdeg)
          Rad_smpp1(iprmin:iprmax, irad) = Rad_smpp1(iprmin:iprmax, irad) &
               & + pp(iprmin:iprmax) * p * (1.0_R_ - pdeg)
          p = p * pdeg
       end if
       !// The excess should be collected and redistributed later.
       !   At present, it is truncated simply, biasing the domain average.

       ! Mapping parameters
       rnx = min(Rad_anxr, rnx)
       rny = min(Rad_anyr, rny)
       rnx1 = real(ixr, R_) - Rad_facx * x2 - 0.5_R_
       rny1 = real(iyr, R_) - Rad_facy * y2 - 0.5_R_
       nxdlo = int(rnx + rnx1)
       nxdhi = int(rnx - rnx1)
       nydlo = int(rny + rny1)
       nydhi = int(rny - rny1)
       call mcarUtl__getGtab1(ixrlut, KND, nxdlo, nxdhi, Rad_nxr, ixr)
       call mcarUtl__getGtab1(iyrlut, KND, nydlo, nydhi, Rad_nyr, iyr)

       ! Distribute and sample
       nfull = (nxdlo + nxdhi + 1) * (nydlo + nydhi + 1)
       p = p / nfull
       if (Atm_mlay(iz) == 3) then ! 3-D layer
          dxr = Atm_xmax / real(Rad_nxr, R_)
          dyr = Atm_ymax / real(Rad_nyr, R_)
          call mcarUtl__getGtab2(ixlut, KND, nxdlo, nxdhi, x, dxr, Atm_facx, Atm_nx)
          call mcarUtl__getGtab2(iylut, KND, nydlo, nydhi, y, dyr, Atm_facy, Atm_ny)
          n = 0
          do iyd = -nydlo, nydhi
             iyt =  iylut(iyd)
             iyu = iyrlut(iyd)
             do ixd = -nxdlo, nxdhi
                ixt =  ixlut(ixd)
                ixu = ixrlut(ixd)
                if (Atm_scat3d(ixt, iyt, iz, ichi) > 0.0_R_) then
                   n = n + 1
                   Rad_smpr(ixu, iyu, irad) = Rad_smpr(ixu, iyu, irad) + p
                   Rad_smpp(ixu, iyu, iprmin:iprmax, irad) = Rad_smpp(ixu, iyu, iprmin:iprmax, irad) &
                        & + p * pp(iprmin:iprmax)
                   iixr2(n) = ixu
                   iiyr2(n) = iyu
                end if
             end do
          end do
          if (n < nfull .and. n > 0) then
             p = p * (nfull - n) / n
             do i = 1, n
                ixu = iixr2(i)
                iyu = iiyr2(i)
                Rad_smpr(ixu, iyu, irad) = Rad_smpr(ixu, iyu, irad) + p
                Rad_smpp(ixu, iyu, iprmin:iprmax, irad) = Rad_smpp(ixu, iyu, iprmin:iprmax, irad) &
                  & + p * pp(iprmin:iprmax)
             end do
          end if
       else ! 1D, surface, space, or underground
          do iyd = -nydlo, nydhi
             iyu = iyrlut(iyd)
             do ixd = -nxdlo, nxdhi
                ixu = ixrlut(ixd)
                Rad_smpr(ixu, iyu, irad) = Rad_smpr(ixu, iyu, irad) + p
                Rad_smpp(ixu, iyu, iprmin:iprmax, irad) = Rad_smpp(ixu, iyu, iprmin:iprmax, irad) &
                     & + p * pp(iprmin:iprmax)
             end do
          end do
       end if
    end if

  end subroutine mcarRad__samp2_diff


  !+
  ! Angular PDF for one geometry
  !-
  function mcarRad__angDistr_R(imod, isub, pac, ix, iy) result(res)

    integer,  intent(in) :: imod
    !// imod = -1 : local source, 0 : solar source, 1 : thermal source in atmosphere
    integer,  intent(in) :: isub
    real(R_), intent(in) :: pac(:)
    integer,  intent(in) :: ix, iy
    real(R_)  :: res ! result, PDF (1/steradian) for angular (re)distribution
    real(R_),  parameter :: API4 = 0.25_R_ / PI_

    ! Source emission (not surface source)
    if (imod <= 1) then
       if (imod <= -1) then ! local source
          res = pac(2)
       else if (imod == 0) then ! solar source
          res = pac(2) * abs(PhoV_dir(3))
       else ! thermal source emission (atmosphere)
          res = API4
       end if
       !// direct local source beam : ADF = 1/dW, where dW = solid angle of the beam
       !   solar direct beam        : ADF = 1/dWs * |cos(q1)|,
       !     where dWs = integral(for Q=0,dQ){integral(for F=0,2*pi){|cos(q2)|}}
       !   thermal source in atmosphere : ADF = 1/(4*pi)

       ! Other cases
    else if (imod == 2 .or. imod == 4) then ! surface event
       res = mcarSfc__angDistr_R(imod, isub, pac, Pho_dir, PhoV_dir, ix, iy)
    else ! scattering
       res = mcarAtm__angDistr_R(Pho_dir, PhoV_dir, pac(1), Pho_ichi)
    end if

  end function mcarRad__angDistr_R


  !+
  ! Get pathlength statistics
  !-
  subroutine mcarRad__pathLen(plen, pp, iprmin, iprmax)

    real(R_), intent(in) :: plen(0:)
    integer,  intent(out) :: iprmin, iprmax
    real(R_), intent(out) :: pp(:)
    real(R_) :: tlen
    integer :: ipr

    ! Pathlength statistics
    iprmin = 1
    iprmax = 0
    if (Rad_mplen <= 0) then
    else if (Rad_mplen == 1) then ! layer average
       iprmax = Rad_npr
       pp(iprmin:iprmax) = plen(iprmin:iprmax)
    else if (Rad_mplen == 2) then ! slant column
       iprmax = Rad_npr
       do ipr = 1, Rad_npr
          pp(ipr) = sum(Rad_wfunc(1:Rad_nz, ipr) * plen(1:Rad_nz)) ! in the atmosphere
       end do
    else ! total pathlength
       tlen = sum(plen(0 : Rad_nz + 1))
       if (tlen >= Rad_tpmin) then
          ipr = int((tlen - Rad_tpmin) * Rad_tpfac) + 1
          if (ipr <= Rad_npr) then
             iprmin = ipr
             iprmax = ipr
             pp(ipr) = 1.0_R_
          end if
       end if
    end if

  end subroutine mcarRad__pathLen


  !+
  ! Camera image pixel & geometry
  !-
  subroutine mcarRad__pixel(vdir, auc, avc, fuc, fvc, nxr, nyr, sina, alpha, ixr, iyr, idet)

    real(R_), intent(in) :: vdir(:) ! view vector in the camera coordinates
    real(R_), intent(in) :: auc, avc
    real(R_), intent(in) :: fuc, fvc
    integer,  intent(in) :: nxr, nyr
    integer,  intent(out) :: idet
    integer,  intent(out) :: ixr, iyr
    real(R_), intent(out) :: alpha
    real(R_), intent(out) :: sina
    real(R_)  :: cosf, u, sinf, v

    ! Transform
    idet = 0
    sina = sqrt(vdir(1)**2 + vdir(2)**2)
    alpha = atan2(sina, vdir(3)) ! theta [0,pi/2]

    ! Image pixel indexes
    if (sina < 1.0e-17_R_) then ! center of camera iamge
       return
    else if (Rad_mpmap == 1) then ! rectangular mapping
       cosf = vdir(1) / sina
       sinf = vdir(2) / sina
       u = auc + alpha * cosf ! U-axis --> image pixel index
       v = avc + alpha * sinf ! V-axis --> image line index
       ixr = int(u * fuc) + 1
       iyr = int(v * fvc) + 1
       if (u < 0.0_R_ .or. ixr > nxr) return
       if (v < 0.0_R_ .or. iyr > nyr) return
    else ! polar mapping
       u = alpha
       v = atan2(-vdir(2), -vdir(1)) + PI_ ! azimuth [0,2*pi] from U-axis
       ixr = int(u * fuc) + 1
       iyr = int(v * fvc) + 1
       if (ixr < 1 .or. ixr > nxr) return !print *, ixr, u
       if (iyr < 1 .or. iyr > nyr) return !print *, iyr, v
    end if
    idet = 1

  end subroutine mcarRad__pixel


  !+
  ! Sample contributions for the radiances averaged over area and solid angle
  ! Def: Radiance [W/m2/str], L = sum{E/dS/(dA*|cos(q)|)}
  !-
  subroutine mcarRad__samp3(izn)

    integer, intent(in) :: izn ! layer node index [0,nz]
    real(R_), parameter :: UZMIN = 1.0e-4_R_ ! min of cos(q)
    real(R_) :: auc, avc, fuc, fvc, vdir(3), sina, alpha, p
    real(R_) :: pp(Rad_npr) !(AUTO)
    integer  :: irad, ixr, iyr, idet, iprmin, iprmax

    ! Loop for cameras located at the level
    if (Rad_mrkind /= 3 .or. Rad_mznode(izn) == 0) return
    do irad = 1, Rad_nrad
       if (Rad_iizn(irad) /= izn) cycle ! node level is different

       ! Camera coordinates
       vdir(:) = matmul(Rad_crot(:,:,irad), -Pho_dir) ! view vector (camera coordinates)
       if (vdir(3) < Rad_pixp(1, irad)) cycle ! out of the FOV cone

       ! Image mapping
       auc = Rad_pixp(2, irad)
       avc = Rad_pixp(3, irad)
       fuc = Rad_pixp(4, irad)
       fvc = Rad_pixp(5, irad)
       call mcarRad__pixel(vdir, auc, avc, fuc, fvc, Rad_nxr, Rad_nyr, sina, alpha, ixr, iyr, idet)
       if (idet /= 1) cycle

       ! Sample
       p = Pho_eone * Pho_wgt / nonZero_R(abs(Pho_dir(3)), UZMIN)
       if (Rad_mrproj == 1) p = p * vdir(3)
       call mcarRad__pathLen(Pho_plen, pp, iprmin, iprmax)
       Rad_smpr(ixr, iyr, irad) = Rad_smpr(ixr, iyr, irad) + p
       Rad_smpp(ixr, iyr, iprmin:iprmax, irad) = Rad_smpp(ixr, iyr, iprmin:iprmax, irad) &
            & + p * pp(iprmin:iprmax)
    end do

  end subroutine mcarRad__samp3


  !+
  ! Volume rendering for visualization or approximate radiance calculations
  !-
  subroutine mcarRad__render(nrank, irank, nray, wgt)

    integer,  intent(in) :: nrank ! # of PEs
    integer,  intent(in) :: irank ! my PE index
    integer,  intent(in) :: nray  ! # of rays per angular bin
    real(R_), intent(in) :: wgt   ! some weight for radiances
    real(R_), parameter :: TAUMAX = 16.0_R_
    real(R_) :: vdir(3), vloc(3), p, qc, uc, vc, w, sump, sumw
    real(R_) :: cosq0, cosq1, cosq, sinq, phi0, phi1, phi, tau, rn1, rn2
    integer  :: irad, ixr, iyr, iray, itask, ntask, nxyr, ixyr

    ! Task sharing among PEs
    nxyr = Rad_nxr * Rad_nyr
    ntask = nxyr * Rad_nrad

    ! Loops for tasks (for cameras and pixels)
    do itask = irank + 1, ntask, nrank ! MPI parallel
       irad = (itask - 1) / nxyr + 1
       ixyr = itask - nxyr * (irad - 1)
       iyr = (ixyr - 1) / Rad_nxr + 1
       ixr = ixyr - Rad_nxr * (iyr - 1)

       ! Loop for rays
       sump = 0.0_R_
       sumw = 0.0_R_
       do iray = 1, nray
          if (nray == 1) then
             rn1 = 0.5_R_
             rn2 = 0.5_R_
          else
             rn1 = mseq_rand_R()
             rn2 = mseq_rand_R()
          end if
          
          ! View location & direction
          if (Rad_mrkind == 1) then ! 1st-kind of radiance
             if (Rad_mpmap == 1) then ! rectanguler mapping
                uc = (ixr - 1 + rn1) / Rad_pixp(4, irad) - Rad_pixp(2, irad) ! MC sampling
                vc = (iyr - 1 + rn2) / Rad_pixp(5, irad) - Rad_pixp(3, irad)
                qc = sqrt(uc**2 + vc**2)
                cosq = cos(qc)
                sinq = sin(qc)
                if (cosq < Rad_pixp(1, irad)) cycle
                vdir(1) = sinq * uc / qc ! view vector (camera coordinates)
                vdir(2) = sinq * vc / qc
                vdir(3) = cosq
             else ! polar mapping
                cosq0 = cos((ixr - 1) / Rad_pixp(4, irad))
                cosq1 = cos(ixr / Rad_pixp(4, irad)) ! [0,1]
                phi0 = (iyr - 1) / Rad_pixp(5, irad)
                phi1 = iyr / Rad_pixp(5, irad)
                cosq = cosq0 + rn1 * (cosq1 - cosq0) ! MC sampling
                phi  = phi0  + rn2 * (phi1  - phi0)
                sinq = sqrt(max(0.0_R_, 1.0_R_ - cosq**2))
                vdir(1) = sinq * cos(phi) ! view vector (camera coordinates)
                vdir(2) = sinq * sin(phi)
                vdir(3) = cosq
             end if
             vloc(:) = Rad_cloc(:, irad)
             vdir(:) = matmul(transpose(Rad_crot(:,:,irad)), vdir) ! camera to world
             vdir(3) = nonZero_R(vdir(3), 1.0e-3_R_) ! to avoid completely horizontal directions
             w = 1.0_R_
             if (Rad_mrproj == 0) w = cosq
          else ! 2nd-kind of radiance
             vloc(1) = Atm_xmax / Rad_nxr * (ixr - 1 + rn1)
             vloc(2) = Atm_ymax / Rad_nyr * (iyr - 1 + rn2)
             vloc(3) = Rad_cloc(3, irad)
             vdir(:) = Rad_cdir(:, irad)
             w = 1.0_R_
          end if

          ! Sample contribution
          call mcarVis__render(TAUMAX, vdir, vloc, tau, p)
          sumw = sumw + w
          sump = sump + w * p
       end do
       if (sumw > 0.0_R_) Rad_smpr(ixr, iyr, irad) = Rad_smpr(ixr, iyr, irad) + wgt * sump / sumw
    end do

  end subroutine mcarRad__render
  

  !+
  ! Sample BRDF, just for validation of surface reflection modeling
  !// This procedure is simply for a test.
  !-
  subroutine mcarRad__sampBRDF(ix, iy, imod, isub, pac)

    integer,  intent(in) :: imod
    integer,  intent(in) :: isub
    integer,  intent(in) :: ix, iy
    real(R_), intent(in) :: pac(:)
    real(R_) :: fuc, fvc, vdir(3), p, dzmin, pps, ds
    real(R_) :: rn1, rn2, cosq0, cosq1, cosq, sinq, phi0, phi1, phi
    integer  :: irad, ixr, iyr

    ! Only at the BOA
    if (Rad_mrkind /= -1) return
    dzmin = (Atm_zgrd(Atm_nz) - Atm_zgrd(0)) * RSML_
    if (abs(Pho_loc(3) - Atm_zgrd(0)) > dzmin) return

    ! Loops for cameras and pixels
    do irad = 1, Rad_nrad
       do iyr = 1, Rad_nyr
          do ixr = 1, Rad_nxr

             ! The reflection vector
             fuc = Rad_pixp(4, irad)
             fvc = Rad_pixp(5, irad)
             if (Rad_mpmap == 2) then ! polar mapping
                cosq1 = cos(ixr / fuc) ! [0,1]
                cosq0 = cos((ixr - 1) / fuc)
                rn1 = mseq_rand_R()
                rn2 = mseq_rand_R()
                cosq = (1.0_R_ - rn1) * cosq0 + rn1 * cosq1 ! MC sampling
                sinq = sqrt(max(0.0_R_, 1.0_R_ - cosq**2))
                phi1 = iyr / fvc
                phi0 = (iyr - 1) / fvc
                phi =  (1.0_R_ - rn2) * phi0 + rn2 * phi1 ! MC sampling
                ds = abs((cosq0 - cosq1) * (phi1 - phi0)) ! solid angle (steradian)
                vdir(1) = sinq * cos(phi)
                vdir(2) = sinq * sin(phi)
                vdir(3) = cosq
                vdir(:) = matmul(transpose(Rad_crot(:,:,irad)), vdir) ! camera to world
                vdir(:) = -vdir(:) ! reverse vector
             else
                cycle
             end if

             ! Sample contribution
             pps = mcarSfc__angDistr_R(imod, isub, pac, Pho_dir, vdir, ix, iy)
             !// PDF = |cos(q1)| * BRDF / (albedo)
             p = Pho_eone * Pho_wgt / nonZero_R(abs(vdir(3)), REPS_) * pps * ds
             if (Rad_mrproj == 1) p = p * vdir(3)
             Rad_smpr(ixr, iyr, irad) = Rad_smpr(ixr, iyr, irad) + p
          end do
       end do
    end do

  end subroutine mcarRad__sampBRDF


#if UseMPI == 1
  !+
  ! Reduce the paralell computed integrals
  !-
  subroutine mcarRad__reduce(irank, iroot)

    include 'inc_mpi.f90'
    integer, intent(in) :: irank
    integer, intent(in) :: iroot
    real(R_)  :: wrk0(Rad_nxr * Rad_nyr), wrk1(Rad_nxr * Rad_nyr) !(AUTO)
    real(R_)  :: wrk2(Rad_nrad), wrk3(Rad_nrad), s1(2), s2(2)
    integer   :: ierr, irad, ipr, na(1), nb(2)

    ! Array sizes
    nb(1) = size(Rad_smpr, 1)
    nb(2) = size(Rad_smpr, 2)
    na(1) = nb(1) * nb(2)

    ! Radiances & pathlength statistics
    do irad = 1, Rad_nrad
       wrk0(1:na(1)) = reshape(Rad_smpr(:, :, irad), na)
       call MPI_Reduce(wrk0, wrk1, na(1), MPI_R_, MPI_SUM, iroot, MPI_COMM_WORLD, ierr)
       if (irank == iroot) Rad_smpr(:, :, irad) = reshape(wrk1(1:na(1)), nb)
       do ipr = 1, Rad_npr
          wrk0(1:na(1)) = reshape(Rad_smpp(:, :, ipr, irad), na)
          call MPI_Reduce(wrk0, wrk1, na(1), MPI_R_, MPI_SUM, iroot, MPI_COMM_WORLD, ierr)
          if (irank == iroot) Rad_smpp(:, :, ipr, irad) = reshape(wrk1(1:na(1)), nb)
       end do
    end do

    ! Truncated extrems
    wrk2(:) = Rad_smpr1(:) 
    call MPI_Reduce(wrk2, wrk3, Rad_nrad, MPI_R_, MPI_SUM, iroot, MPI_COMM_WORLD, ierr)
    if (irank == iroot) Rad_smpr1(:) = wrk3(:)
    do ipr = 1, Rad_npr
       wrk2(:) = Rad_smpp1(ipr, :) 
       call MPI_Reduce(wrk2, wrk3, Rad_nrad, MPI_R_, MPI_SUM, iroot, MPI_COMM_WORLD, ierr)
       if (irank == iroot) Rad_smpp1(ipr, :) = wrk3(:)
    end do

    ! Sources
    s1(1) = Rad_psrc
    s1(2) = Rad_esrc
    call MPI_Reduce(s1, s2, 2, MPI_R_, MPI_SUM, iroot, MPI_COMM_WORLD, ierr)
    if (irank == iroot) then
       Rad_psrc = s2(1)
       Rad_esrc = s2(2)
    end if

  end subroutine mcarRad__reduce
#endif


  !+
  ! Normalize results
  !-
  subroutine mcarRad__normal(atot) 

    real(R_),  intent(in) :: atot ! total horizontal area (m^2) of the domain
    real(R_), allocatable :: wrk(:,:)
    real(RD_) :: sum1, fac, fnorm
    real(R_)  :: pnorm
    integer   :: irad, ipr

    ! Initialize
    allocate (wrk(Rad_nxr, Rad_nyr))
    fnorm = 1.0_RD_ / Rad_psrc ! divided by # of source photons
    if (Rad_mrkind == 2) then
       fnorm = fnorm * real(Rad_nxr * Rad_nyr, RD_) / atot ! [m^-2] divided by the pixel area
    else if (Rad_mrkind == 3) then
       fnorm = fnorm / atot ! [m^-2] divided by the area of horizontal plane
    end if
    Rad_esrc = Rad_esrc / Rad_psrc ! average source power [W]

    ! Loop for cameras
    do irad = 1, Rad_nrad

       ! Redistribution of truncated radiance contributions
       sum1 = sum(Rad_smpr(:,:,irad))
       if (sum1 > RDSML_) then
          fac = (sum1 + Rad_smpr1(irad)) / sum1 ! redistribution factor
          Rad_smpr(:,:,irad) = fac * Rad_smpr(:,:,irad) ! corrected radiance samples
       end if
       !// Note: Redistribution over different anguler bins is inherently problematic.
       !   It is unknown how different the method should be by weighting function of the integral.
       !   That is, different methods should be used for mrproj = 1 and 2.
       !   An impovement could be future work.

       ! Normalize pathlength statistics
       do ipr = 1, Rad_npr
          sum1 = sum(Rad_smpp(:, :, ipr, irad))
          fac = 1.0_RD_
          if (sum1 > RDSML_) fac = (sum1 + Rad_smpp1(ipr, irad)) / sum1 ! redistribution factor
          if (Rad_mplen == 1) then
             pnorm = 1.0_R_ / Atm_zdep(ipr)
          else if (Rad_mplen == 2) then
             pnorm = 1.0_R_ / sum(Rad_wfunc(1:Rad_nz, ipr) * Atm_zdep(1:Rad_nz)) ! 1/VCD
          else
             pnorm = 1.0_R_
          end if
          Rad_smpp(:,:, ipr, irad) = pnorm * fac * Rad_smpp(:,:, ipr, irad) &
               & / max(RDSML_, Rad_smpr(:,:, irad))
       end do

       ! Normalize radiances
       if (Rad_mrkind == 2) then
          Rad_smpr(:,:,irad) = fnorm * Rad_smpr(:,:,irad) / abs(Rad_cdir(3, irad)) ! normalize
       else
          if (Rad_mpmap == 1) then
             wrk(:,:) = mcarUtl__imgSolAng_rec_2R(Rad_nxr, Rad_nyr, Rad_mrproj, Rad_pixp(2, irad), &
                  & Rad_pixp(3, irad))
          else
             wrk(:,:) = mcarUtl__imgSolAng_pol_2R(Rad_nxr, Rad_nyr, Rad_mrproj, Rad_pixp(2, irad), &
                  & Rad_pixp(3, irad))
          end if
          Rad_smpr(:,:,irad) = fnorm * Rad_smpr(:,:,irad) / wrk(:,:)
       end if

       ! Correct parallax
       if (Rad_mrkind == 2) then
          call shiftPixsMrkind2(Rad_smpr(:,:,irad), irad)
          do ipr = 1, Rad_npr
             call shiftPixsMrkind2(Rad_smpp(:,:,ipr,irad), irad)
          end do
       end if
    end do

  end subroutine mcarRad__normal


  !+
  ! Write out results to a binary file
  !-
  subroutine mcarRad__writeBin(iuo) 

    integer,   intent(in) :: iuo
    real(R4_), allocatable :: wrk(:,:)
    integer   :: irad, ipr

    ! Write out
    allocate (wrk(Rad_nxr, Rad_nyr))
    do irad = 1, Rad_nrad
       wrk(:,:) = Rad_smpr(:, :, irad)
       call bin_write_i2R4(iuo, wrk)
    end do
    do irad = 1, Rad_nrad
       do ipr = 1, Rad_npr
          wrk(:,:) = Rad_smpp(:, :, ipr, irad)
          call bin_write_i2R4(iuo, wrk)
       end do
    end do
    deallocate (wrk)
    Rad_ndset = Rad_ndset + 1 ! count # of datasets

  end subroutine mcarRad__writeBin


  !+
  ! Write out results to a text file
  !-
  subroutine mcarRad__writeTxt(iuo) 

    integer,   intent(in) :: iuo
    real(R4_), allocatable :: wrk(:,:)
    integer   :: irad, ipr, ixr, iyr

    ! Write out
    allocate (wrk(Rad_nxr, Rad_nyr))
    do irad = 1, Rad_nrad
       wrk(:,:) = Rad_smpr(:, :, irad)
       write (iuo, '(a, 1i4)', err=1) '# Radiances, irad =', irad
       write (iuo, *, err=1) ((wrk(ixr, iyr), ixr = 1, Rad_nxr), iyr = 1, Rad_nyr)
    end do
    do irad = 1, Rad_nrad
       write (iuo, '(a, 2i4)', err=1) '# Pathlength statistics, irad = ', irad
       do ipr = 1, Rad_npr
          wrk(:,:) = Rad_smpp(:, :, ipr, irad)
          write (iuo, *, err=1) ((wrk(ixr, iyr), ixr = 1, Rad_nxr), iyr = 1, Rad_nyr)
       end do
    end do
    deallocate (wrk)

    return
1   call err_write(1, iuo, 'mcarRad__writeTxt')

  end subroutine mcarRad__writeTxt


  !+
  ! Make a (GrADS) control file
  !-
  subroutine mcarRad__writeCtl(iuc, outfile) 

    integer, intent(in) :: iuc
    character(*), intent(in) :: outfile
    integer :: nvar, irad, nz
    character(80) :: varstr(Rad_nrad*2) !(AUTO)

    nz = max(1, Rad_npr)
    nvar = 0
    do irad = 1, Rad_nrad
       nvar = nvar + 1
       varstr(nvar) = 'a'//trim(num2str_AN(irad))//' '//trim(num2str_AN(1))//' 99 Radiance'
    end do
    if (Rad_npr >= 1) then
       do irad = 1, Rad_nrad
          nvar = nvar + 1
          varstr(nvar) = 'b'//trim(num2str_AN(irad))//' '//trim(num2str_AN(Rad_npr)) &
               & //' 99 Pathlength Statistics'
       end do
    end if
    call gradsCtl_write(iuc, fileName_AN(outfile), Rad_nxr, Rad_nyr, nz, Rad_ndset, nvar, varstr)

  end subroutine mcarRad__writeCtl

  !+
  ! Fix the position of the output when mrkind=2
  !+
  subroutine shiftPixsMrkind2(wrk, irad)

    real(R_), intent(inout) :: wrk(:,:)
    integer, intent(in) :: irad
    real(R_)  :: zdiff
    real(R_)  :: dif_x, dif_y
    real(R_)  :: sin_the, cos_the, sin_phi, cos_phi, tan_the

    sin_the =   sin(Rad_the(irad) / 180.0_R_ * PI_)
    cos_the = - cos(Rad_the(irad) / 180.0_R_ * PI_)
    sin_phi =   sin(Rad_phi(irad) / 180.0_R_ * PI_)
    cos_phi =   cos(Rad_phi(irad) / 180.0_R_ * PI_)
    tan_the = sin_the / cos_the
    zdiff = Rad_zloc(irad) - Rad_zref(irad) ! height difference between zloc and zref
    dif_x = zdiff * tan_the * cos_phi ! difference along x axis (m)
    dif_y = zdiff * tan_the * sin_phi ! difference along y axis (m)
    wrk = cshift(wrk, - nint(mod(dif_x, Atm_xmax) / Atm_xmax * real(Rad_nxr, R_)), 1) ! shift in dimension x
    wrk = cshift(wrk, - nint(mod(dif_y, Atm_ymax) / Atm_ymax * real(Rad_nyr, R_)), 2) ! shift in dimension y

  end subroutine shiftPixsMrkind2

end module mcarRad
