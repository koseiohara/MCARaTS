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
! MCARaTS world module
!-
module mcarWld 

  use globals
#if UseMPI == 1
  use globals_MPI
#endif
  use hparx
  use mcarPho  ! Pho_*, readable & writable
  use mcarAtm  ! Atm_*, read-only
  use mcarSfc  ! Sfc_*, read-only
  use mcarSrc  ! Src_*, read-only
  use mcarSca
  use mcarFlx
  use mcarRad
  use mcarVis
  use mcarUtl
  implicit none
  private

  ! Public
  public :: mcarWld__set_env  !* Set environment of this module *
  public :: mcarWld__init     !* Initialize this module with user configuration *
  public :: mcarWld__jobs     !* Do jobs for user configuration *
  public :: mcarWld__final    !* Finalize this module *

  ! (Private) User variables in namelist (1) for initialization
  integer,  save :: Wld_mverb = 0   ! verbose mode (0=quiet, 1=yes, 2=more, 3=most)
  integer,  save :: Wld_jseed = 0   ! seed for random number generator (0 for automatic)
  integer,  save :: Wld_mbswap = 0  ! flag for byte swapping for binary input files (0=no, 1=yes)
  integer,  save :: Wld_mtarget = 1 ! flag for target quantities
  !// = 1 : fluxes and heating rates, calculated by MC
  !     2 : radiances, calculated by MC
  !     3 : quasi-radiances (or some signals), calculated by volume rendering
  integer,  save :: Wld_moptim = 2  ! flag for optimization of calculation techniques
  !// = -2 : no tuning (use default optimization)
  !     -1 : no optimization (deactivate all optimization)
  !      0 : unbiased optimization (deactivate all biasing optimizations)
  !      1 : conservative optimizations (possibly for smaller biases)
  !      2 : standard optimizations (recommended)
  !      3 : quick-and-dirty optimizations (for speed when small biases are acceptable)
  integer,  save :: Wld_njob = 1    ! # of jobs in a single experiment
  !// An experiment corresponds to one namelist input file and one output file.
  namelist /mcarWld_nml_init/Wld_mverb, Wld_jseed, Wld_mbswap, Wld_mtarget, Wld_moptim, Wld_njob

  ! (Private) User variables in namelist (2) for core calculations
  integer,  save :: Wld_nplcf  = 0  ! power exponent for local collision forcing (LCF)
  real(R_), save :: Wld_rmlcf  = 100.0_R_ ! mode distance for LCF
  real(R_), save :: Wld_dtaumin = 0.0_R_  ! min layer optical thickness for collision forcing
  namelist /mcarWld_nml_job/Wld_rmlcf, Wld_nplcf, Wld_dtaumin

  ! Private
  integer,  save :: Wld_iroot = 0  ! root PE rank index
  integer,  save :: Wld_irank = 0  ! my PE rank index
  integer,  save :: Wld_nrank = 1  ! total # of PEs
  integer,  save :: Wld_iuo = 10   ! file unit index for output file
  integer,  save :: Wld_iuc = 11   ! file unit index for control file of the output file
  integer,  save :: Wld_mdet(4) = (/1, 1, 0, 0/)
  !// mdet(1) = 1 : to activate sampling of flux
  !   mdet(2) = 1 : to activate sampling of heating rate
  !   mdet(3) = 1 : to activate sampling of radiance of the 1st/2nd kind
  !   mdet(4) = 1 : to activate sampling of radiance of the 3rd kind
  character(256), save :: Wld_outfile = 'out' ! name of output file

  !- Common usage of parameter packet, pac(1:6)
  ! This array is used in this module, mcarRad, and mcarSfc.
  !  imod=-1 : Local source
  !    pac(1) = qmaxsrc : a half cone angle (radian) of the source emission
  !    pac(2) : angular PDF
  !  imod=0 : Solar source
  !    pac(1) = qmaxsrc : a half cone angle (radian) of the source emission
  !    pac(2) : angular PDF = adfsrc * abs(cosQ0)
  !  imod=1 : Thermal source in the atmosphere
  !  imod=2 : Thermal source from the surface
  !    pac(1:3) = Sfc_psfc2d(ixb, iyb, 1:3)
  !  imod=3 : Scattering in the atmosphere
  !    pac(1) = apf   : phase function parameter
  !    pac(2) = pdltc : probability of physical scattering
  !  imod=4 : Reflection at the surface
  !   isfc=0 : black
  !   isfc=1 : Lambertian
  !     pac(1) : albedo
  !     pac(6) : albedo
  !   isfc=2 : DSM (diffuse-specular mixture)
  !     pac(1) : diffuse albedo
  !     pac(2) : fraction of diffuse reflection
  !     pac(3) : real part of refractive index
  !     pac(4) : imaginary part of refractive index
  !     pac(5) : variance of facet slope
  !     pac(6) : albedo
  !   isfc=3 : RPV (Rahmen-Pinty-Verstraete)
  !     pac(1) : "rho", the first parameter
  !     pac(2) : "k", exponent parameter
  !     pac(3) : "THETA", asymmetry parameter
  !     pac(4) : "delta", sharpness of the hot spot
  !     pac(5) : function A
  !     pac(6) : albedo
  !   isfc=4 : LSRT (Li-Sparse-Ross-Thick)
  !     pac(1) : "k_L", weight for Lambertian reflection
  !     pac(2) : "k_g", weight for geometrical optics
  !     pac(3) : "k_v", weight for volume scattering
  !     pac(6) : albedo
  !// Some elements are not used in some procedures.

contains

  !+
  ! Set environment for this module
  !-
  subroutine mcarWld__set_env(iroot, irank, nrank) 

    integer, intent(in), optional :: iroot ! index of root PE
    integer, intent(in), optional :: irank ! index of my PE
    integer, intent(in), optional :: nrank ! # of PEs
    Wld_iroot = iroot
    Wld_irank = irank
    Wld_nrank = nrank
    !print *, iroot, irank, nrank !TEST

  end subroutine mcarWld__set_env

  !+
  ! Set technical parameters
  !-
  subroutine mcarWld__set_tech(dtaumin, nplcf) 

    real(R_), intent(in), optional :: dtaumin !, rmlcf
    integer,  intent(in), optional :: nplcf
    Wld_dtaumin = dtaumin
    Wld_nplcf  = nplcf
    
  end subroutine mcarWld__set_tech

  !+
  ! Initialize the MCARaTS world (this module and all sub-modules)
  !-
  subroutine mcarWld__init(iui, filepath, outfile, isol)

    integer, intent(in) :: iui ! file unit index for the namelist input
    character(*), intent(in) :: filepath
    character(*), intent(in) :: outfile
    integer, intent(in) :: isol ! solver flag (0=F3D, 1=P3D, 2=1D(ICA))
    integer :: jseed, nrec, nchi, mflx, mhrt, nxr, nyr, mrkind

    ! Read in namelists (1) for initialization
    call mcarWld__user_init(iui)
    call mcarSca__tune_init(Wld_moptim) ! tune technical parameters
    call mcarFlx__tune_init(Wld_moptim)
    call mcarPho__tune_init(Wld_mtarget)
    call mcarSca__user_init(iui, filepath)
    call mcarAtm__user_init(iui, filepath)
    call mcarSfc__user_init(iui, filepath)
    call mcarSrc__user_init(iui)
    call mcarFlx__user_init(iui)
    call mcarRad__user_init(iui)
    call mcarVis__user_init(iui)
    call mcarPho__user_init(iui)
    call mcarRad__qry_flags(mrkind)
    call mcarWld__tune_jobs(mrkind) ! tune technical parameters
    call mcarAtm__tune_jobs(Wld_moptim)
    call mcarSfc__tune_jobs(Wld_moptim)
    call mcarRad__tune_jobs(Wld_moptim)

    ! Initialize utilities
    if (Wld_irank /= Wld_iroot) Wld_mverb = 0 ! Only the root PE can report.
    call mtab_init(200000, 200000, 50000)
    jseed = Wld_jseed
    if (jseed <= 0) jseed = rand_newSeed_I() ! a time-dependent seed
    if (Wld_mverb >= 1) write (*,*)
    if (Wld_mverb >= 1) write (*,*) '* Random number seed =', jseed, 'for PE', Wld_irank
    jseed = rand_seed_I(Wld_irank, jseed) ! possibly modifiy the seed
    call mseq_init(jseed) ! initialize the random number generator

    ! Initialize sub-modules
    call mcarSca__qry_dims(nchi=nchi)
    call mcarSca__init()
    call mcarAtm__init(Wld_irank, Wld_iroot, nchi)
    call mcarSfc__init(Wld_irank, Wld_iroot)
    call mcarSrc__init() ! using mcarAtm variables
    if (Wld_mtarget == 1) call mcarFlx__init(Atm_nx, Atm_ny, Atm_nz, Atm_iz3l, Atm_iz3u)
    if (Wld_mtarget >= 2) call mcarRad__init(Atm_nz)
    if (Wld_mtarget == 3) call mcarVis__init()
    call mcarPho__init(Atm_nz, nchi)
    call mcarPho__set_iso_SS(isol) ! setup solver
    call mcarSca__prep(Wld_irank, Wld_nrank, Wld_mbswap) ! prepare scattering tables

    ! Setup MC samplers
    Wld_mdet(:) = 0
    if (Wld_mtarget == 1) then
       call mcarFlx__qry_flags(mflx, mhrt)
       if (mflx >= 1) Wld_mdet(1) = 1
       if (mhrt >= 1) Wld_mdet(2) = 1
    else if (Wld_mtarget == 2) then
       call mcarRad__qry_flags(mrkind)
       if (mrkind == 3) Wld_mdet(3) = 1
       if (mrkind /= 3) Wld_mdet(4) = 1
    end if

    ! Open output file
    if (Wld_mtarget == 1) then
       nrec = recordLen_I(1, Atm_nx * Atm_ny, 4) ! Apr. 27, 2012
    else
       call mcarRad__qry_dims(nxr, nyr)
       nrec = recordLen_I(1, nxr * nyr, 4) ! Apr. 27, 2012
    end if
    Wld_iuo = freeUnit_I(Wld_iuo)
    Wld_outfile = outfile
    if (Wld_irank == Wld_iroot) call open_dir(Wld_iuo, outfile, nrec, 'unknown') ! binary
    !if (Wld_irank == Wld_iroot) call open_seq(Wld_iuo, outfile, 'unknown') ! text
    if (Wld_mverb >= 3) write (*,*) 'mcarWld__init: finish'

  end subroutine mcarWld__init

  !+
  ! Read in namelist variables (1) for initialization
  !-
  subroutine mcarWld__user_init(iu) 

    integer, intent(in) :: iu   ! file unit index for namelist input
    integer :: ios
    read (iu, nml=mcarWld_nml_init, iostat=ios)
    call err_read(ios, iu, 'mcarWld__user_init: Invalid namelist input for mcarWld_nml_init.')

  end subroutine mcarWld__user_init

  !+
  ! Finalize this module
  !-
  subroutine mcarWld__final() 

    integer :: iuc

    ! Complete the output
    if (Wld_irank == Wld_iroot) then ! make a control file
       iuc = freeUnit_I(12)
       call open_seq(iuc, trim(Wld_outfile)//'.ctl', 'unknown')
       if (Wld_mtarget == 1) call mcarFlx__writeCtl(iuc, Wld_outfile)
       if (Wld_mtarget >= 2) call mcarRad__writeCtl(iuc, Wld_outfile)
       close (iuc)
    end if
    close (Wld_iuo)

    ! Finalize all
    call mcarSca__final()
    call mcarAtm__final()
    call mcarSfc__final()
    call mcarSrc__final()
    call mcarFlx__final()
    call mcarRad__final()
    call mcarVis__final()
    call mcarPho__final()

  end subroutine mcarWld__final

  !+
  ! Tune a preset for technical parameters
  !// Note: This procedure should be called after mcarAtm and mcarFlx are initialized.
  !-
  subroutine mcarWld__tune_jobs(mrkind) 

    integer, intent(in) :: mrkind
    integer :: mflx, mhrt
    real(R_) :: dtaumin

    if (Wld_moptim <= -2) then ! as defaults
    else if (Wld_moptim == -1) then ! no optimization
       call mcarWld__set_tech(0.0_R_, 0)
    else
       if (Wld_mtarget == 1) then ! flux & HR
          call mcarFlx__qry_flags(mflx, mhrt)
          if (mhrt == 1) then
             dtaumin = 0.5_R_ / Atm_nz * (Wld_moptim + 1)
             call mcarWld__set_tech(dtaumin, 0) ! GCF
          else
             call mcarWld__set_tech(0.0_R_, 0) ! no CF
          end if
       else if (mrkind == 1) then ! radiance of the 1st kind
          dtaumin = 0.1_R_ * (Wld_moptim + 1)
          call mcarWld__set_tech(dtaumin, 1) ! LCF
       else if (mrkind == 2 .or. mrkind == 3) then ! radiance of the 2nd/3rd kind
          call mcarWld__set_tech(0.0_R_, 0)
       end if
    end if
    
  end subroutine mcarWld__tune_jobs

  !+
  ! Do jobs for user configuration
  !-
  subroutine mcarWld__jobs(iui, ptot) 

    integer,  intent(in) :: iui  ! file unit index for the namelist input
    real(R_), intent(in) :: ptot ! # of used photons
    real(R_) :: tmpmin, tmpmax, tmpmin1, tmpmax1, extmin(Atm_nz) !(AUTO)
    real(R_) :: q0, f0, dq0, sdir(3), t0, t1, sloc(3)
    integer :: isrc, npf, ikd, ijob, idread0

    ! Loop for jobs
    do ijob = 1, Wld_njob
       if (Wld_mverb >= 1) call cpu_time(t0) ! for timing study

       ! Read in namelists (2) for this job
       idread0 = 1
       if (ijob >= 2) idread0 = 0
       call mcarWld__user_job(iui)
       call mcarAtm__user_job(iui, idread0)
       call mcarSfc__user_job(iui, idread0)
       call mcarSrc__user_job(iui)
       call mcarRad__user_job(iui)

       ! Prepare for core calculations
       call mcarSca__qry_dims(npf)
       call mcarAtm__prep_job(Wld_irank, Wld_mbswap, npf, tmpmin, tmpmax)
       call mcarSfc__prep_job(Wld_irank, Wld_mbswap, Atm_xmax, Atm_ymax, tmpmin1, tmpmax1)
       !// using Atm_* set in mcarAtm
       if (Wld_mtarget >= 2) call mcarRad__prep_job()
       if (Wld_mtarget == 3) call mcarVis__prep_job()
       tmpmin = min(tmpmin, tmpmin1)
       tmpmax = max(tmpmax, tmpmax1)
       if (Wld_mverb >= 2) call mcarAtm__check_opt3DStats() 

       ! Min extinction coefficients (for collision-forcing)
       if (Wld_nplcf == 0 .or. Wld_mtarget == 1) then ! nonlocal collision forcing
          extmin(1:Atm_nz) = Wld_dtaumin / Atm_zdep(1:Atm_nz)
       else ! local collision forcing (force collision near cameras)
          call mcarRad__extMin(Wld_dtaumin, Wld_nplcf, Wld_rmlcf, extmin)
       end if
       call mcarAtm__set_extMin(extmin)
       if (Wld_mverb >= 3) write (*,*) 'mcarWld__jobs: Preprocess finish'

       ! Loop for sources
       do isrc = 1, Src_nsrc
          if (Wld_mverb >= 2) write (*,*) 'mcarWld__jobs: start, isrc =', isrc
          if (Wld_mtarget == 1) call mcarFlx__zero() ! flux and/or HR
          if (Wld_mtarget >= 2) call mcarRad__zero() ! radiance
          call mcarSrc__prep_isrc(isrc, tmpmin, tmpmax, sdir, q0, f0, dq0)
          if (Wld_mtarget == 3) then ! volume rendering
             do ikd = 1, Atm_nkd
                call mcarVis__prep_prop(1, ikd) !mdens)
                call mcarVis__prep_srcj(Src_mtype(isrc), sdir, Src_flx(isrc))
                call mcarRad__render(Wld_nrank, Wld_irank, 1, Atm_wkd(ikd)) ! nray=1
             end do
          else ! forward MC
             sloc(1) = Src_xpos(isrc) * Atm_xmax
             sloc(2) = Src_ypos(isrc) * Atm_ymax
             sloc(3) = Src_zloc(isrc)
             call mcarUtl__horiShift(sloc(1), 1.0_R_, 0.0_R_, Atm_xmax) ! cyclic condition only
             call mcarUtl__horiShift(sloc(2), 1.0_R_, 0.0_R_, Atm_ymax) ! in the current codes
             call mcarWld__fMC(ptot, q0, f0, dq0, Src_mphi(isrc), Src_mtype(isrc), &
                  & Src_flx(isrc), sloc, Src_apsize(isrc)) !, Src_dwlen(isrc))
          end if
          call mcarWld__job_writeData()
          if (Wld_mverb >= 3) write (*,*) 'mcarWld__jobs: end, isrc =', isrc
       end do
       if (Wld_mverb >= 2) write (*,*) 'mcarWld__jobs: finished a job', ijob
       if (Wld_mverb >= 1) then ! report performance
          call cpu_time(t1)
          write (*, *) 'CPU time (s) =', t1 - t0
          write (*, *) 'Time per simulation (s) =', (t1 - t0) / Src_nsrc
       end if
    end do

  end subroutine mcarWld__jobs


  !+
  ! Read in namelist variables (2) for core calculations
  !-
  subroutine mcarWld__user_job(iu) 

    integer, intent(in) :: iu   ! file unit index for namelist input
    integer :: ios
    read (iu, nml=mcarWld_nml_job, iostat=ios)
    call err_read(ios, iu, 'mcarWld__user_job: Invalid namelist input for mcarWld_nml_job.')
    call check_iI('mcarWld__user_job: Wld_nplcf',   Wld_nplcf, 0, 3)
    call check_iR('mcarWld__user_job: Wld_dtaumin', Wld_dtaumin, 0.0_R_,  100.0_R_)
    call check_iR('mcarWld__user_job: Wld_rmlcf',   Wld_rmlcf,   0.0_R_)

  end subroutine mcarWld__user_job


  !+
  ! Write out results to a file
  !-
  subroutine mcarWld__job_writeData() 

    real(R_) :: atot

    ! Reduce MPI results
#if UseMPI == 1
    if (Wld_nrank > 1) then
       if (Wld_mtarget == 1) call mcarFlx__reduce(Wld_irank, Wld_iroot)
       if (Wld_mtarget >= 2) call mcarRad__reduce(Wld_irank, Wld_iroot)
    end if
#endif

    ! Write out
    if (Wld_irank == Wld_iroot) then
       atot = Atm_xmax * Atm_ymax ! horizontal cross section of the domain
       if (Wld_mtarget == 1) then
          call mcarFlx__normal(atot)
          call mcarFlx__writeBin(Wld_iuo)
       else if (Wld_mtarget == 2) then
          call mcarRad__normal(atot)
          call mcarRad__writeBin(Wld_iuo)
       else if (Wld_mtarget == 3) then
          call mcarRad__writeBin(Wld_iuo)
       end if
    end if

  end subroutine mcarWld__job_writeData


  !+
  !, dwlen)
  ! Simulations for a source condition
  !-
  subroutine mcarWld__fMC(ptot, thesrc, phisrc, qmaxsrc, mphi, mstype, fsrc, sloc, asize) 

    integer,  intent(in) :: mstype
    integer,  intent(in) :: mphi
    real(R_), intent(in) :: asize
    real(R_), intent(in) :: qmaxsrc ! a half cone angle (radian)
    !real(R_), intent(in) :: dwlen
    real(R_), intent(in) :: fsrc
    real(R_), intent(in) :: phisrc
    real(R_), intent(in) :: ptot
    real(R_), intent(in) :: thesrc
    real(R_), intent(in) :: sloc(:)
    integer, allocatable :: ixalut(:), iyalut(:), ixblut(:), iyblut(:)
    integer   :: iwrk, ip, ix0, iy0, iz0, np, nquad, ikd, nxa, nya, nxb, nyb, nwrk, nwrkmin
    real(R_)  :: pdfsrc(Atm_nz + 1), coemt(Atm_nz, 2), hsrc(Atm_nx, Atm_ny) !(AUTO)
    real(R_)  :: adfsrc, anp, anp_atm, anp_sfc, chi0, cosdc, dx0, dy0
    real(R_)  :: ebsrc, fbsrc, fcone, etot, fcol, dir0(3), fsatm, fssfc
    real(R_)  :: ppcone, sindq, x0min, y0min, atot, hdwx, hdwy !, dscol

    ! Randomly distribute works
    nwrkmin = min(20 * Wld_nrank, 20000 + Atm_nx * Atm_ny * Atm_nz / 3) ! complicated!
    call mcarWld__fMC_numWorks(Atm_nx, Atm_ny, nwrkmin, ptot, nxa, nya)
    allocate (ixalut(nxa * nya), iyalut(nxa * nya))
    call mcarWld__fMC_randWorks(nxa, nya, ixalut, iyalut)
    nwrk = nxa * nya
    if (Wld_mverb >= 3) write (*,*) 'set source x-y grids, nxa, nya =', nxa, nya
    if (mstype >= 2) then ! with thermal source
       nwrkmin = min(20 * Wld_nrank, 20000 + Sfc_nxb * Sfc_nyb) ! complicated!
       call mcarWld__fMC_numWorks(Sfc_nxb, Sfc_nyb, nwrkmin, ptot, nxb, nyb)
       allocate (ixblut(nxb * nyb), iyblut(nxb * nyb))
       call mcarWld__fMC_randWorks(nxb, nyb, ixblut, iyblut)
       nwrk = nwrk + nxb * nyb
       if (Wld_mverb >= 3) write (*,*) 'set source x-y grids, nxb, nyb =', nxb, nyb
    end if
    !// Used later: nwrk, ixalut, iyalut, nxa, nya, (ixblut, iyblut, nxb, nyb)

    ! Source parameters
    hdwx = 0.5_R_ * Atm_xmax / sqrt(ptot)
    hdwy = 0.5_R_ * Atm_ymax / sqrt(ptot)
    atot = Atm_xmax * Atm_ymax ! horizontal cross section of the domain
    dir0(:) = geo3D_aVec_1R(1.0_R_, thesrc, phisrc)
    if (mstype <= 0) then  ! point source
       sindq = sin(max(REPS_, qmaxsrc))
       if (qmaxsrc < PIH_) then
          cosdc = sindq**2 / (1.0_R_ + cos(qmaxsrc))
       else
          cosdc = 1.0_R_ - cos(qmaxsrc)
       end if
       adfsrc = 1.0_R_ / (PI2_ * cosdc) ! 1/(solid angle)
       iz0 = gridIdx_bin0_I(Atm_zgrd, sloc(3), 0, Atm_nz+1) ! = [0, nz+1]
       chi0 = 1.0_R_ - 0.5_R_ * cosdc ! source directionality
       ebsrc = fsrc ! spectral source energy (W/micron)
    else if (mstype <= 2) then ! solar source
       nquad = 64
       call mcarUtl__coneProj(nquad, thesrc, qmaxsrc, cosdc, ppcone, fcone)
       adfsrc = 1.0_R_ / ppcone
       chi0 = 1.0_R_ - 0.5_R_ * cosdc
       ebsrc = fsrc * fcone * atot
    else ! thermal emission alone
       ebsrc = 0.0_R_
    end if
    !if (dwlen > RSML_) ebsrc = ebsrc * dwlen ! (W)
    !// Now, no wavelendth width is multiplied, so that output results of this code will be 
    !   always spectrally-averaged quantities.
    fbsrc = ebsrc / atot ! averaged local/solar source flux
    !dscol = atot / (Atm_nx * Atm_ny) ! column cross section

    ! Loop for the k-distribution
    do ikd = 1, Atm_nkd
       if (Wld_mverb >= 2) write (*,*) '** Simulating for ikd=', ikd

       ! Preprocess for this k-term
       !   These have not yet been parallelized with MPI.
       call mcarAtm__prep_extMax(ikd)
       anp = ptot * Atm_wkd(ikd) / real(nxa * nya, R_) ! average # of photons per column
       if (mstype <= 1) then
          etot = ebsrc
          anp_atm = anp
          anp_sfc = 0.0_R_
       else
          call mcarWld__fMC_srcDist(ebsrc, ikd, ptot, hsrc, fsatm, fssfc, etot)
          anp_atm = ptot * Atm_wkd(ikd) / real(nxa * nya, R_) * fsatm
          anp_sfc = ptot * Atm_wkd(ikd) / real(nxb * nyb, R_) * fssfc
          if (Wld_mverb >= 2) write (*,*) ' fsatm, fssfc, etot :', fsatm, fssfc, etot
       end if

       ! Loop for works (possibly with different source columns)
       do iwrk = Wld_irank + 1, nwrk, Wld_nrank
          if (Wld_mverb >= 3) write (*,*) ' iwrk :', iwrk

          ! Configure photons for this work
          if (mstype <= 0) then
          else if (iwrk <= nxa * nya) then ! source from sun/atmosphere, mstype = 1/2/3
             dx0 = Atm_xmax / nxa
             dy0 = Atm_ymax / nya
             x0min = dx0 * real(ixalut(iwrk) - 1, R_)
             y0min = dy0 * real(iyalut(iwrk) - 1, R_)
             ix0 = gridIdx_unif0_I(Atm_xgrd, (x0min + 0.5_R_ * dx0), 1, Atm_nx)
             iy0 = gridIdx_unif0_I(Atm_ygrd, (y0min + 0.5_R_ * dy0), 1, Atm_ny)
             if (mstype >= 2) then ! solar + thermal
                anp = anp_atm * hsrc(ix0, iy0) * (Atm_nx * Atm_ny)
                call mcarWld__fMC_srcProf(ix0, iy0, fbsrc, ikd, pdfsrc, coemt, fcol)
             end if
             !if (Wld_mverb >= 3) write (*,*) ' ix0, iy0 :', ix0, iy0
          else ! source from surface, mstype = 2 or 3
             dx0 = Atm_xmax / nxb
             dy0 = Atm_ymax / nyb
             x0min = dx0 * real(ixblut(iwrk - nxa * nya) - 1, R_)
             y0min = dy0 * real(iyblut(iwrk - nxa * nya) - 1, R_)
             ix0 = min(Sfc_nxb, int((x0min + 0.5_R_ * dx0) * Sfc_facx) + 1)
             iy0 = min(Sfc_nyb, int((y0min + 0.5_R_ * dy0) * Sfc_facy) + 1)
             anp = anp_sfc * Src_fsfc(ix0, iy0) * (Sfc_nxb * Sfc_nyb)
          end if

          ! Start
          call mcarPho__register(ikd, etot, hdwx, hdwy) ! register common photon status
          np = int(anp)
          if (mseq_rand_R() < anp - real(np, R_)) np = np + 1
          if (Wld_mtarget == 1) call mcarFlx__sampSrc(np)
          if (Wld_mtarget == 2) call mcarRad__sampSrc(np)
          if (Wld_mverb >= 3) write (*,*) '  np :', np

          ! Monte Carlo core
          do ip = 1, np ! loop for photons
             if (mstype <= 0) then ! localized source
                call mcarWld__fMC_srcLoc(sloc, dir0, iz0, chi0, qmaxsrc, cosdc, adfsrc, asize)
             else if (iwrk <= nxa * nya) then ! solar or atmospheric source
                call mcarWld__fMC_srcAtm(mstype, mphi, pdfsrc, coemt, x0min, y0min, &
                     & dx0, dy0, ix0, iy0, dir0, chi0, qmaxsrc, cosdc, adfsrc)
             else ! surface source
                call mcarWld__fMC_srcSfc(x0min, y0min, dx0, dy0, ix0, iy0)
             end if
             do
                if (Pho_dead .or. Pho_iz <= 0 .or. Pho_iz >= Atm_nz + 1) exit
                call mcarWld__rt_atmos()
                if (Pho_iz <= 0) call mcarWld__rt_surf()
             end do
          end do

       end do
    end do
    deallocate (ixalut, iyalut)
    if (mstype >= 2) deallocate (ixblut, iyblut)

  end subroutine mcarWld__fMC


  !+
  ! Determines numbers of works, taking into account several factors
  !-
  subroutine mcarWld__fMC_numWorks(nx, ny, nwrkmin, ptot, nxw, nyw)

#if UseMPI == 1
    include 'inc_mpi.f90'
#endif
   integer,  intent(in) :: nx, ny
    integer,  intent(in) :: nwrkmin
    real(R_), intent(in) :: ptot
    integer,  intent(out) :: nxw, nyw
    integer :: i
    real(R_) :: anp, fanp

    ! Determine
    if (Wld_irank == Wld_iroot) then ! the root PE only
       anp = ptot / real(Atm_nkd * nx * ny, R_) ! average # of photons per column
       nxw = nx
       nyw = ny
       fanp = 1.0_R_
       if (Wld_nrank > 1) then
          do i = 1, 1000
             if (nxw * nyw >= nwrkmin .or. fanp * anp <= 8.0_R_) exit ! could be simplified
             nxw = 2 * nxw
             nyw = 2 * nyw
             fanp = fanp * 0.25_R_
          end do
       end if
    end if

    ! MPI broadcast
#if UseMPI == 1
    call MPI_Bcast(nxw, 1, MPI_INTEGER, Wld_iroot, MPI_COMM_WORLD, i)
    call MPI_Bcast(nyw, 1, MPI_INTEGER, Wld_iroot, MPI_COMM_WORLD, i)
#endif

  end subroutine mcarWld__fMC_numWorks


  !+
  ! Randomly distribute the source packets, generating LUTs for x and y grid mapping
  !-
  subroutine mcarWld__fMC_randWorks(nxa, nya, ixalut, iyalut)

#if UseMPI == 1
    include 'inc_mpi.f90'
#endif
    integer,  intent(in) :: nxa, nya
    integer,  intent(out) :: ixalut(:), iyalut(:)
    integer, allocatable :: iixy(:)
    integer :: iwrk, ires, ires1, ixa, i, iya, nres, nwrk

    ! Generate a LUT
    nwrk = nxa * nya
    if (Wld_irank == Wld_iroot) then ! the root PE only
       allocate (iixy(nwrk))
       do iwrk = 1, nwrk
          iixy(iwrk) = iwrk
       end do
       nres = nwrk
       do iwrk = 1, nwrk
          ires = min(nres, int(nres * mseq_rand_R() + 1.0_R_)) ! random grid index
          i = iixy(ires)
          iya = (i - 1) / nxa + 1
          ixa = i - nxa * (iya - 1)
          ixalut(iwrk) = ixa
          iyalut(iwrk) = iya
          nres = nres - 1
          do ires1 = ires, nres
             iixy(ires1) = iixy(ires1 + 1) ! remaining grid list
          end do
       end do
       deallocate (iixy)
    end if

    ! MPI broadcast
#if UseMPI == 1
    call MPI_Bcast(ixalut, nwrk, MPI_INTEGER, Wld_iroot, MPI_COMM_WORLD, i)
    call MPI_Bcast(iyalut, nwrk, MPI_INTEGER, Wld_iroot, MPI_COMM_WORLD, i)
#endif

  end subroutine mcarWld__fMC_randWorks


  !+
  ! Compute 2-D spatial distribution of emitted sources
  !-
  subroutine mcarWld__fMC_srcDist(esol, ikd, ptot, hsrc, fsatm, fssfc, etot)

    real(R_), intent(in) :: esol ! TOA solar source power
    integer,  intent(in) :: ikd
    real(R_), intent(in) :: ptot
    real(R_), intent(out) :: hsrc(:, :) ! probability of source power for each column
    real(R_), intent(out) :: fsatm, fssfc ! fractions of atmosphere and surface sources
    real(R_), intent(out) :: etot ! total source power
    integer   ::ix, iy, iz
    real(R_)  :: dscola, b, bb, bt, rat, etoa, eatm, ecell, stot

    ! Constants
    dscola = 4.0_R_ * PI_ * Atm_xmax * Atm_ymax / (Atm_nx * Atm_ny) ! coloumn cross section
    etoa = esol / (Atm_nx * Atm_ny) ! TOA source power per column

    ! Integrate atmospheric source power
    eatm = 0.0_R_
    do iy = 1, Atm_ny
       do ix = 1, Atm_nx
          stot = etoa
          do iz = 1, Atm_nz
             bb = Src_batm(ix, iy, iz)
             bt = Src_batm(ix, iy, iz + 1)
             rat = bt / bb
             if (abs(rat - 1.0_R_) < 0.0001_R_) then
                b = (bb + bt) * 0.5_R_
             else
                b = (bt - bb) / log(rat) ! assuming exponential profile
             end if
             ecell = b * dscola * Atm_zdep(iz) * Atm_abst3d(ix, iy, iz, ikd) ! source power of the cell (W)
             stot = stot + ecell
             if (Wld_mdet(2) == 1 .and. Wld_irank == Wld_iroot) then ! root PE only
                ecell = -ecell * Atm_wkd(ikd) * ptot ! cooling power * (# of photons) (W)
                call mcarFlx__sampHeat(ix, iy, iz, ecell) ! IR cooling rate
             end if
          end do
          hsrc(ix, iy) = stot
       end do
       eatm = eatm + sum(hsrc(1:Atm_nx, iy))
    end do

    ! Normalize
    etot = eatm + Src_esfc ! total source power
    fsatm = eatm / etot
    fssfc = Src_esfc / etot
    hsrc(:,:) = hsrc(:,:) / eatm ! probability of source power for each column

  end subroutine mcarWld__fMC_srcDist


  !+
  ! Compute atmosphere emitted source profile possibly with solar source
  !-
  subroutine mcarWld__fMC_srcProf(ix, iy, ftoa, ikd, pdfsrc, coemt, ftot)

    real(R_), intent(in) :: ftoa
    integer,  intent(in) :: ikd
    integer,  intent(in) :: ix, iy
    real(R_), intent(out) :: ftot
    real(R_), intent(out) :: pdfsrc(:) ! CDF (top to bottom)
    real(R_), intent(out) :: coemt(:,:)
    integer   :: iz
    real(R_)  :: bb, bt, rat, b

    ! Profile of source fluxes
    pdfsrc(Atm_nz + 1) = ftoa ! solar source
    do iz = Atm_nz, 1, -1 ! top to bottom
       bb = Src_batm(ix, iy, iz)
       bt = Src_batm(ix, iy, iz + 1)
       rat = bt / bb
       if (abs(rat - 1.0_R_) < 0.0001_R_) then
          b = (bb + bt) * 0.5_R_
       else
          b = (bt - bb) / log(rat) ! assuming exponential profile
       end if
       coemt(iz, 1) = rat
       coemt(iz, 2) = 1.0_R_ / log(rat)       
       pdfsrc(iz) = pdfsrc(iz + 1) + 4.0_R_ * PI_ * b * Atm_abst3d(ix, iy, iz, ikd) * Atm_zdep(iz) ! flux
    end do
    ftot = pdfsrc(1)
    pdfsrc(:) = pdfsrc(:) / ftot ! normalize
    pdfsrc(1) = 1.0_R_

  end subroutine mcarWld__fMC_srcProf


  !+
  ! Conical emission from a point source, initialization of a photon packet.
  !-
  subroutine mcarWld__fMC_srcLoc(sloc, dir0, iz0, chi0, qmaxsrc, cosdc, adfsrc, asize)

    real(R_), intent(in) :: adfsrc
    real(R_), intent(in) :: chi0
    real(R_), intent(in) :: cosdc
    real(R_), intent(in) :: asize
    real(R_), intent(in) :: qmaxsrc
    integer,  intent(in) :: iz0
    real(R_), intent(in) :: dir0(:)
    real(R_), intent(in) :: sloc(:)
    real(R_)  :: pac(10), path, znew, adir(3)
    integer   :: ix, iy, izold, izn

    ! General initializations
    call mcarPho__rt_newPhoton(1.0_R_, dir0, sloc, iz0, 1.0_R_)
    call mcarPho__rt_newStatus(chi0) ! Now, iso = 1 (direct beam)

    ! A random location for emission
    Pho_loc(:) = geo3D_randLoc_circle_1R(Pho_loc(:), Pho_dir, asize)
    Pho_iz = gridIdx_seq0_I(Atm_zgrd, Pho_loc(3), Pho_iz, 0, Atm_nz + 1) !TEST

    ! Local estimates of direct beam
    if (Wld_mdet(4) == 1) then
       ix = 0
       iy = 0
       pac(1) = qmaxsrc
       pac(2) = adfsrc
       call mcarRad__samp(ix, iy, -1, 0, pac)
    end if
    call mcarPho__rt_scatOrders() ! update status
    if (Pho_dead) return

    ! Emission vector
    adir(:) = geo3D_randDir_cone_1R(Pho_dir, cosdc, 0) ! mproj = 0 (isotropic)
    call mcarPho__rt_update(adir, adfsrc) ! update status

    ! Transport to the domain boundary
    if (Pho_iz >= 1 .and. Pho_iz <= Atm_nz) then
    else if (real(1 - Pho_iz, R_) * Pho_dir(3) <= 0.0_R_) then ! escape to space or ground
       !call mcarPho__rt_scaleWgt(0.0_R_)
       !Pho_dead = .true.
    else ! enter the domain
       izold = Pho_iz ! should be = 0 or nz+1, otherwise fatal error
       if (Pho_iz >= Atm_nz + 1) then ! from TOA
          znew = Atm_zgrd(Atm_nz)
          Pho_iz = Atm_nz
       else ! form BOA
          znew = Atm_zgrd(0)
          Pho_iz = 1
       end if
       path = (znew - Pho_loc(3)) / Pho_dir(3)
       Pho_loc(3) = znew
       if (Pho_mv3D) call mcarAtm__rt_renewLocXY(Pho_dir, path, Pho_loc) 
       call mcarPho__rt_hist_rec(izold, path)
       !Pho_sdif = Pho_sdif + (1.0_R_ - Pho_chi) * path
       !Pho_plen(izold) = Pho_plen(izold) + path
       izn = Pho_iz - Pho_itr(3)
       if (Wld_mdet(1) == 1) call mcarFlx__sampFlux(Pho_loc(1), Pho_loc(2), izn)
       if (Wld_mdet(3) == 1) call mcarRad__samp3(izn)
    end if

  end subroutine mcarWld__fMC_srcLoc


  !+
  ! Source: initialization of a photon packet
  !-
  subroutine mcarWld__fMC_srcAtm(mstype, mphi, pdfsrc, coemt, &
       & x0min, y0min, dx0, dy0, ix0, iy0, dir0, chi0, qmaxsrc, cosdc, adfsrc)

    real(R_), intent(in) :: adfsrc ! for solar source, angular PDF = adfsrc * abs(cosQ0)
    real(R_), intent(in) :: cosdc
    real(R_), intent(in) :: dx0, dy0
    integer,  intent(in) :: ix0, iy0
    integer,  intent(in) :: mstype
    !// = 1, external source from TOA (e.g., solar sources)
    !     2, external source plus internal source
    !     3, internal source (thermal emission)
    integer,  intent(in) :: mphi
    real(R_), intent(in) :: x0min, y0min
    real(R_), intent(in) :: pdfsrc(:)
    real(R_), intent(in) :: coemt(:,:)
    real(R_), intent(in) :: qmaxsrc
    real(R_), intent(inout) :: chi0
    real(R_), intent(inout) :: dir0(:)
    real(R_),  parameter :: PPSISO = 0.25_R_ / PI_
    real(R_)  :: pac(10), cosf, path, r2, rn, sinf, vz, loc0(3), adir(3), wgt0, pps
    integer   :: iz0, imod, ix, iy

    ! General initializations
    loc0(1) = x0min + dx0 * mseq_rand_R()
    loc0(2) = y0min + dy0 * mseq_rand_R()
    ix = ix0
    iy = iy0
    wgt0 = 1.0_R_

    ! Layer determination
    if (mstype == 1) then
       iz0 = Atm_nz + 1
    else
       rn = mseq_rand_R()
       do iz0 = Atm_nz + 1, 1, -1
          if (rn < pdfsrc(iz0)) exit
       end do !// should be replaced by the binary search
       iz0 = max(iz0, 1) ! because iz0 = 0 if rn = 1.0
    end if

    ! External source from top
    if (iz0 > Atm_nz) then
       if (mphi >= 1) then ! random azimuth
          vz = sqrt(dir0(1)**2 + dir0(2)**2)
          call rand_point_circle(1.0e-12_R_, r2, sinf, cosf)
          dir0(1) = vz * cosf
          dir0(2) = vz * sinf
       end if
       loc0(3) = Atm_zgrd(iz0)
       call mcarPho__rt_newPhoton(wgt0, dir0, loc0, iz0, 1.0_R_)
       call mcarPho__rt_newStatus(chi0)
       if (Wld_mdet(4) == 1) then
          imod = 0
          pac(1) = qmaxsrc
          pac(2) = adfsrc
          call mcarRad__samp(ix, iy, imod, 0, pac)
       end if
       call mcarPho__rt_scatOrders() ! update status
       if (Pho_dead) return
       adir(:) = geo3D_randDir_cone_1R(Pho_dir, cosdc, 1) ! mproj = 1 (projection)
       pps = adfsrc * abs(adir(3))
       call mcarPho__rt_update(adir, pps) ! update status
       if (Pho_itr(3) == 1) then ! exit to space
          Pho_iz = Pho_iz + 1
          !Pho_dead = .true.
          return
       end if
       path = (Atm_zgrd(Pho_iz - 1 + Pho_itr(3)) - Pho_loc(3)) / Pho_dir(3)
       call mcarPho__rt_hist_rec(Pho_iz, path)
       call mcarAtm__rt_bound1D(Pho_itr(3), Pho_iz, Pho_loc(3))
       if (Pho_mv3D) call mcarAtm__rt_renewLocXY(Pho_dir, path, Pho_loc) 
       if (Wld_mdet(1) == 1) call mcarFlx__sampFlux(Pho_loc(1), Pho_loc(2), Pho_iz)
       if (Wld_mdet(3) == 1) call mcarRad__samp3(Pho_iz)

       ! Thermal emission from atmosphere
    else
       chi0 = 0.5_R_
       rn = mseq_rand_R()
       if (abs(coemt(iz0, 1) - 1.0_R_) > 0.0001_R_) &
            & rn = coemt(iz0, 2) * log(1.0_R_ - rn + rn * coemt(iz0, 1)) ! assuming exponential profile
       loc0(3) = (1.0_R_ - rn) * Atm_zgrd(iz0 - 1) + rn * Atm_zgrd(iz0)
       call mcarPho__rt_newPhoton(wgt0, dir0, loc0, iz0, 1.0_R_)
       call mcarPho__rt_newStatus(chi0)
       if (Wld_mdet(4) >= 1) then
          imod = 1
          call mcarRad__samp(ix, iy, imod, 0, pac)
       end if
       call mcarPho__rt_scatOrders() ! update status
       if (Pho_dead) return
       call rand_point_sphere(1.0e-13_R_, rn, adir(1), adir(2), adir(3)) ! isotropic
       call mcarPho__rt_update(adir, PPSISO) ! update status
       call mcarWld__fMC_netEmit(ix, iy)
    end if

  end subroutine mcarWld__fMC_srcAtm


  !+
  ! Surface source: initialization of a photon packet
  !-
  subroutine mcarWld__fMC_srcSfc(x0min, y0min, dx0, dy0, ixb, iyb)

    real(R_), intent(in) :: x0min, y0min
    real(R_), intent(in) :: dx0, dy0
    integer,  intent(in) :: ixb, iyb
    real(R_)  :: chi0, pac(10), loc0(3), dir0(3), adir(3), wgt0, pps
    integer   :: iz0, imod, isfc

    ! General initializations
    iz0 = 0
    loc0(1) = x0min + dx0 * mseq_rand_R()
    loc0(2) = y0min + dy0 * mseq_rand_R()
    loc0(3) = Atm_zgrd(0)
    dir0(1:3) = (/0.0_R_, 0.0_R_, 1.0_R_/)
    wgt0 = 1.0_R_
    chi0 = 0.5_R_

    ! New photon
    isfc = Sfc_jsfc2d(ixb, iyb)
    pac(1:3) = Sfc_psfc2d(ixb, iyb, 1:3)
    call mcarPho__rt_newPhoton(wgt0, dir0, loc0, iz0, 1.0_R_)
    call mcarPho__rt_newStatus(chi0)
    if (Wld_mdet(4) == 1) then
       imod = 2
       call mcarRad__samp(ixb, iyb, imod, isfc, pac)
    end if
    call mcarPho__rt_scatOrders() ! update status
    if (Pho_dead) return

    ! New status
    adir(:) = mcarSfc__fMC_emit_1R(ixb, iyb, isfc, pac)
    pps = abs(adir(3)) * mcarSfc__emittance_R(adir(3), ixb, iyb, isfc, pac)
    call mcarPho__rt_update(adir, pps) ! update status
    Pho_iz = 1
    if (Wld_mdet(1) == 1) call mcarFlx__sampFlux(Pho_loc(1), Pho_loc(2), Pho_iz - 1)
    if (Wld_mdet(3) == 1) call mcarRad__samp3(Pho_iz - 1)

  end subroutine mcarWld__fMC_srcSfc


  !+
  ! Net emission from the initial grid, transmission or 1st-order scattering
  !-
  subroutine mcarWld__fMC_netEmit(ix, iy)

    integer,  intent(inout) :: ix, iy
    real(R_)  :: pac(10), apf, abst, fpath, sca, scat, path, scaran, tsum, w0
    real(R_)  :: pps, adir(3), gstar
    integer   :: icase, ippo, izn

    ! Layer specifications
    abst = Atm_abst3d(ix, iy, Pho_iz, Pho_ikd)
    w0 = Pho_wgt

    ! Loop for scattering events
    do
       ! Initial path
       scat = Atm_scat3d(ix, iy, Pho_iz, Pho_ichi) ! total scattering coefficient
       if (Pho_mv3D) then
          call mcarAtm__rt_path3D(Pho_loc, ix, iy, Pho_iz, Pho_dir, Pho_itr, icase, path)
       else
          path = (Atm_zgrd(Pho_iz - 1 + Pho_itr(3)) - Pho_loc(3)) * Pho_dir(3)
          icase = 3
       end if
       if (Pho_ftau >= scat * path) then ! escape from the box
          call mcarPho__rt_scaleWgt(mtab_expNX_R(abst * path)) ! scale weight
          exit
       end if

       ! Move to a scattering point and scale weight
       fpath = Pho_ftau / scat ! Note: fatal error if scat = 0 and Pho_ftau < 0
       call mcarPho__rt_scaleWgt(mtab_expNX_R(abst * fpath))
       if (Pho_dead) exit
       Pho_loc(3) = Pho_loc(3) + fpath * Pho_dir(3) ! scattering point
       if (Pho_mv3D) call mcarAtm__rt_renewLocXY(Pho_dir, fpath, Pho_loc)
       call mcarPho__rt_hist_rec(Pho_iz, fpath)
       !Pho_plen(Pho_iz) = Pho_plen(Pho_iz) + fpath
       !Pho_sdif = Pho_sdif + (1.0_R_ - Pho_chi) * fpath

       ! Determine a scattering component
       scaran = scat * mseq_rand_R() ! random scattering coeff.
       tsum = 0.0_R_
       do ippo = 1, Atm_nppo(Pho_iz)
          sca = Atm_scap3d(ix, iy, Pho_iz, ippo)
          apf = Atm_apfp3d(ix, iy, Pho_iz, ippo)
          if (sca > scat * REPS_) then ! to check active components only
             !// Note: scat is total SCALED scattering coefficient. 
             !   This is used here as an approximate guide of NON-SCALED scattering coefficient.
             tsum = tsum + sca * mcarSca__fracDelC_R(apf, Pho_ichi)
             if (scaran <= tsum) exit
          end if
       end do
       pac(1) = apf
       pac(2) = 1.0_R_ ! collision forcing is not used
       gstar = abs(mcarSca__asymStar_R(apf, Pho_ichi))

       ! Scattering & local estimates
       call mcarPho__rt_newStatus(gstar) ! new status
       if (Wld_mdet(4) == 1) call mcarRad__samp(ix, iy, 3, 0, pac)
       call mcarPho__rt_scatOrders() ! update status
       if (Pho_dead) exit
       adir(:) = Pho_dir(:)
       call mcarSca__newDirec(apf, Pho_ichi, adir, pps) ! a new direction
       call mcarPho__rt_update(adir, pps) ! update status
    end do

    ! Integrate HR
    if (Wld_mdet(2) == 1) call mcarFlx__sampHeat(ix, iy, Pho_iz, Pho_eone * (w0 - Pho_wgt))
    if (Pho_dead) return

    ! Motion to the boundary
    call mcarPho__rt_ftau_update(scat * path)
    call mcarPho__rt_hist_rec(Pho_iz, path)
    !Pho_plen(Pho_iz) = Pho_plen(Pho_iz) + path
    !Pho_sdif = Pho_sdif + (1.0_R_ - Pho_chi) * path
    if (Pho_mv3D) then
       call mcarAtm__rt_bound3D(icase, path, Pho_dir, Pho_itr, Pho_loc, ix, iy, Pho_iz)
    else
       call mcarAtm__rt_bound1D(Pho_itr(3), Pho_iz, Pho_loc(3))
    end if

    ! Integrate flux
    if (icase == 3) then
       izn = Pho_iz - Pho_itr(3)
       if (Wld_mdet(1) == 1) call mcarFlx__sampFlux(Pho_loc(1), Pho_loc(2), izn)
       if (Wld_mdet(3) == 1) call mcarRad__samp3(izn)
    end if
    call mcarPho__rt_killSmall()

  end subroutine mcarWld__fMC_netEmit


  !+
  ! Traces a trajectory in model atmosphere
  !// Algorithm
  ! - Atm_mlay=1: The photon penetrates the layer straightly. The photon moves
  !    layer-by-layer, and heating rate is integrated at each layer.
  ! - Atm_mlay=2 or 3: Maximum cross section method is used (the photon moves freely in 
  !    each super-voxel). Heating rate is integrated at each collision point.
  !-
  subroutine mcarWld__rt_atmos()

    integer  :: nloop1, nloop3
    real(R_) :: pac(10) ! parameter packet
    !// pac(1): apf,   phase function specification parameter
    !   pac(2): pdltc, probability of physical scattering
    real(R_) :: coas, extm, abst, fpath, path, patha, pps, adir(3), gstar
    integer  :: icase, iloop, isup, ix, ixs, iy, iys, izs, ikind, izn
    !// ikind : a kind of scattering event
    !    = -3 : mathmatical scattering at clear region
    !    = -2 : mathmatical scattering at non-clear region
    !    =  1 : physical scattering

    ! Super-layer indices
    izs = Atm_jzstab(Pho_iz)
    isup = Atm_jsup(izs)

    ! 1-D absorbing super-layers
    if (isup == 1) then
       nloop1 = 1000 + 2 * Atm_nz ! should be modified
       do iloop = 1, nloop1
          path = (Atm_zgrd(Pho_iz - 1 + Pho_itr(3)) - Pho_loc(3)) / Pho_dir(3)
          abst = Atm_abst3d(1, 1, Pho_iz, Pho_ikd)
          if (abst > 0.0_R_) then
             coas = mtab_expNXC_R(abst * path)
             if (Wld_mdet(2) == 1) then ! HR sampling
                patha = -log(max(REPS_, 1.0_R_ - mseq_rand_R() * coas)) / abst
                Pho_loc(3) = Pho_loc(3) + patha * Pho_dir(3)
                if (Pho_mv3D) call mcarAtm__rt_renewLocXY(Pho_dir, patha, Pho_loc) 
                call mcarPho__rt_hist_rec(Pho_iz, patha)
                path = path - patha

                !if (Pho_loc(2) >= 0.0 .and. Pho_loc(2) <= 5000.0) then
                !else
                !   print *, 'abs, ', Pho_loc(:)
                !endif

                call mcarFlx__sampHeatD(1.0_R_, coas, isup)
             end if
             call mcarPho__rt_scaleWgt(1.0_R_ - coas) ! scale weight
             if (Pho_dead) return
          end if
          call mcarPho__rt_hist_rec(Pho_iz, path)
          call mcarAtm__rt_bound1D(Pho_itr(3), Pho_iz, Pho_loc(3))
          if (Pho_mv3D) call mcarAtm__rt_renewLocXY(Pho_dir, path, Pho_loc) ! motion to the boundary
          izn = Pho_iz - Pho_itr(3)
          if (Wld_mdet(1) == 1) call mcarFlx__sampFlux(Pho_loc(1), Pho_loc(2), izn)
          if (Wld_mdet(3) == 1) call mcarRad__samp3(izn)
          if (Atm_jsup(Atm_jzstab(Pho_iz)) /= 1) return
       end do

       ! 1-D/3-D scattering super-layers
    else
       if (isup == 2) then
          ixs = 1
          iys = 1
       else
          ixs = gridIdx_unif0_I(Atm_xsup(:, izs), Pho_loc(1), 1, Atm_nxs(izs))
          iys = gridIdx_unif0_I(Atm_ysup(:, izs), Pho_loc(2), 1, Atm_nys(izs))
       end if
       nloop3 = 1000 + Atm_nx * Atm_ny * Atm_nzs ! should be modified

       ! Loop for super-voxels
       do iloop = 1, nloop3
          call mcarAtm__rt_pathSup(ixs, iys, izs, isup, path, icase)
          extm = Atm_extmax(ixs, iys, izs, Pho_ichi)
          do ! loop for collisions
             if (Pho_ftau >= extm * path) exit ! exit if extm = 0 because ftau should be >= 0
             fpath = Pho_ftau / extm

             call mcarWld__rt_collision(extm, fpath, isup, ix, iy, pac, ikind) ! a collision
             if (Pho_dead) return

             ! Real collision or virtual scattering due to CF
             if (ikind == 1) then
                gstar = abs(mcarSca__asymStar_R(pac(1), Pho_ichi)) ! abs(g*)
                call mcarPho__rt_newStatus(gstar)
                if (Wld_mdet(4) == 1) call mcarRad__samp(ix, iy, 3, 0, pac) ! local estimates
                call mcarPho__rt_killSmall()
                if (Pho_dead) return
                if (rand_isRare_L(pac(2))) then ! real scattering
                   call mcarPho__rt_scatOrders() ! update
                   if (Pho_dead) return
                   adir(:) = Pho_dir(:)
                   
                   !if (abs(pac(1)) < 1.0e-6_R_) print *, pac(1), ix, iy, Pho_iz

                   call mcarSca__newDirec(pac(1), Pho_ichi, adir, pps) ! a new direction
                   call mcarPho__rt_update(adir, pps) ! update status
                   call mcarAtm__rt_pathSup(ixs, iys, izs, isup, path, icase)
                   extm = Atm_extmax(ixs, iys, izs, Pho_ichi)
                else       ! virtual scattering due to CF
                   path = path - fpath
                   call mcarPho__rt_ftau_new(0)
                end if
             else ! virtual scattering due to MCS
                call mcarPho__rt_ftau_new(1)
             end if
          end do
          call mcarPho__rt_ftau_update(extm * path)
          call mcalWld__rt_moveSup(ixs, iys, izs, isup, path, icase)
          if (icase <= 0) return ! exit upward/downward
       end do
    end if

    ! Trap: too many loops
    call mcarPho__rt_scaleWgt(0.0_R_)
    
  end subroutine mcarWld__rt_atmos


  !+
  ! Next position, grid numbers, and optical properties in a 3D-layer
  !-
  subroutine mcarWld__rt_collision(extm, fpath, isup, ix, iy, pac, ikind)

    real(R_), intent(in) :: extm
    real(R_), intent(in) :: fpath
    integer,  intent(out) :: ikind
    integer,  intent(in) :: isup
    integer,  intent(out) :: ix, iy
    real(R_), intent(out) :: pac(:)
    real(R_)  :: abst, coas, extt, extt1, extran, scat, sca, fdltc
    real(R_)  :: scaran, tsum, apf, locnew(3)
    integer   :: ippo, iznew

    ! New position
    locnew(3) = Pho_loc(3) + fpath * Pho_dir(3)
    iznew = gridIdx_seq0_I(Atm_zgrd, locnew(3), Pho_iz, 1, Atm_nz)
    locnew(1:2) = Pho_loc(1:2)
    if (isup == 2) then ! 1-D layer
       if (Pho_mv3D) call mcarAtm__rt_renewLocXY(Pho_dir, fpath, locnew) 
       ix = 1
       iy = 1
    else ! 3-D layer
       if (Pho_mv3D) locnew(1:2) = locnew(1:2) + fpath * Pho_dir(1:2)
       ix = gridIdx_unif0_I(Atm_xgrd, locnew(1), 1, Atm_nx)
       iy = gridIdx_unif0_I(Atm_ygrd, locnew(2), 1, Atm_ny)
    end if

    ! New properties
    scat = Atm_scat3d(ix, iy, iznew, Pho_ichi)
    abst = Atm_abst3d(ix, iy, iznew, Pho_ikd)
    extt = abst + scat
    if (extt <= 0.0_R_) then ! clear voxel
       ikind = -3
       return
    end if
    extt1 = Atm_extscl(Pho_ichi, Pho_ikd) * max(extt, Atm_extmin(iznew)) ! scaling due to CF
    extran = extm * mseq_rand_R() ! random extinction coefficient
    if (extran > extt1) then ! virtual collision due to MCS
       ikind = -2
       return
    end if
    fdltc = extt / extt1 ! probability of real collision
    !// Now, the collision is real collision or virtual scattering due to CF

    ! Motion to the collision point
    call mcarWld__rt_sampSup(locnew(3), iznew, fpath, isup)
    Pho_loc(:) = locnew(:)
    Pho_iz = iznew

    ! Energy deposition (scaled)
    coas = abst / extt1

    !if (Pho_loc(2) >= 0.0 .and. Pho_loc(2) <= 5000.0) then
    !else
    !   print *, 'coll, ', Pho_loc(:), Pho_dir(:)
    !endif

    if (Wld_mdet(2) == 1) call mcarFlx__sampHeatD(fdltc, coas, isup)
    call mcarPho__rt_scaleWgt(1.0_R_ - coas) ! scale weight
    if (Pho_dead) return ! killed

    ! Determine a scattering component
    scaran = scat * (extran / extt1) ! random scattering coeff.
    tsum = 0.0_R_
    do ippo = 1, Atm_nppo(Pho_iz)
       sca = Atm_scap3d(ix, iy, Pho_iz, ippo)
       apf = Atm_apfp3d(ix, iy, Pho_iz, ippo)
       if (sca > scat * REPS_) then ! to check active scattering components only
          !// Note: scat is total SCALED scattering coefficient. 
          !   This is used here as an approximate guide of NON-SCALED scattering coefficient.
          tsum = tsum + sca * mcarSca__fracDelC_R(apf, Pho_ichi) ! sum scaled scattering coefficients
          if (scaran <= tsum) exit
       end if
    end do

    ! Save parameters
    pac(1) = apf
    pac(2) = 1.0_R_  ! probability of physical scattering
    if (extt < extt1) pac(2) = (extt - abst) / (extt1 - abst)
    ikind = 1

  end subroutine mcarWld__rt_collision


  !+
  ! Surface reflection
  !-
  subroutine mcarWld__rt_surf() 

    real(R_)  :: pac(10)              ! surface parameter package
    real(R_)  :: yco(2)               ! coefficients for scattering
    real(R_)  :: alb, gtc, pac5o, yfa1, yfa2, rdir(3), pps
    integer   :: isfc, ixb, iyb
    integer   :: nco

    ! Initializations
    ixb = min(Sfc_nxb, int(Pho_loc(1) * Sfc_facx) + 1)
    iyb = min(Sfc_nyb, int(Pho_loc(2) * Sfc_facy) + 1)
    isfc = Sfc_jsfc2d(ixb, iyb)
    if (isfc <= 0) then ! black
       alb = 0.0_R_
    else
       pac(1:5) = Sfc_psfc2d(ixb, iyb, 1:5)
       if (isfc == 1) then ! Lambertian
          alb = pac(1)
       else
          nco = 2
          if (isfc == 2) nco = 1
          call mcarSfc__intpOptProp(nco, isfc, ixb, iyb, -Pho_dir(3), yfa1, yfa2, yco)
          !yco(1) = max(1.0_R_, yco(1)) * 1.1_R_ ! coefficient b' <-- BUGFIX, Feb. 2009
          ! yco is Smax for DSM. This could be > 1. The above bug may affect angular distribution.
          yco(1) = yco(1) * 1.1_R_ ! Smax or coefficient b'
          if (isfc == 3) pac(5) = yfa1
          alb = mcarSfc__bsAlbedo_R(isfc, pac, yfa1, yfa2, -Pho_dir(3)) ! albedo (can be > 1)
       end if
    end if
    pac(6) = alb

    ! Scale the weight
    call mcarPho__rt_scaleWgt(min(1.0_R_, alb)) ! scale weight
    if (Pho_dead) return

    ! New order of scattering
    if (isfc == 2) then
       gtc = 0.9_R_ * (1.0_R_ - pac(2)) + 0.5_R_ * pac(2)
    else
       gtc = 0.5_R_
    end if

    ! Local estimates
    call mcarPho__rt_newStatus(gtc) ! new status
    if (Wld_mdet(4) == 1) then
       pac5o = pac(5)
       if (isfc == 2) pac(5) = max(0.1_R_ * (1.0_R_ - Pho_chi), pac5o) ! smooth BRDF (temporalily)
       call mcarRad__samp(ixb, iyb, 4, isfc, pac)
       if (isfc == 2) pac(5) = pac5o ! recover the original
    end if

    ! Reflection
    call mcarPho__rt_scatOrders() ! update status
    if (Pho_dead) return
    call mcarSfc__newDirec(isfc, pac, yco, Pho_dir, rdir, pps)
    call mcarPho__rt_update(rdir, pps) ! update status
    Pho_iz = 1
    if (Wld_mdet(1) == 1) call mcarFlx__sampFluxSfc() ! upward flux
    if (Wld_mdet(3) == 1) call mcarRad__samp3(Pho_iz - 1)
    call mcarPho__rt_killSmall() ! Russian roulette
    if (Pho_dead) return

  end subroutine mcarWld__rt_surf


  !+
  ! Movement to the super-grid boundary in 3D layer
  !-
  subroutine mcalWld__rt_moveSup(ixs, iys, izs, isup, path, icase)

    integer,  intent(inout) :: icase
    integer,  intent(inout) :: isup
    integer,  intent(inout) :: ixs, iys, izs
    real(R_), intent(in) :: path
    integer   :: iznew
    real(R_)  :: znew

    ! Z-plane
    if (icase == 3) then
       if (Pho_itr(3) == 0) then
          znew = Atm_zsup(izs - 1)
          iznew = Atm_jzupr(izs - 1)
       else
          znew = Atm_zsup(izs)
          iznew = Atm_jzupr(izs) + 1
       end if
       call mcarWld__rt_sampSup(znew, iznew, path, isup)
       Pho_loc(3) = znew
       Pho_iz = iznew
       if (Pho_mv3D) then
          if (isup == 2) then ! 1-D super-layer
             call mcarAtm__rt_renewLocXY(Pho_dir, path, Pho_loc) 
          else                ! 3-D super-layer
             Pho_loc(1:2) = Pho_loc(1:2) + path * Pho_dir(1:2)
          end if
       end if
       izs = Atm_jzstab(Pho_iz)
       isup = Atm_jsup(izs)
       if (isup < 2 .or. isup > 3) then ! exit from the super-layers
          icase = -icase
          return
       end if
       if (isup == 2) then
          ixs = 1
          iys = 1
       else
          ixs = gridIdx_unif0_I(Atm_xsup(:, izs), Pho_loc(1), 1, Atm_nxs(izs))
          iys = gridIdx_unif0_I(Atm_ysup(:, izs), Pho_loc(2), 1, Atm_nys(izs))
       end if

       ! Else
    else
       znew = Pho_loc(3) + path * Pho_dir(3)
       iznew = gridIdx_seq0_I(Atm_zgrd, znew, Pho_iz, 1, Atm_nz)
       !if (znew < Atm_zgrd(iznew-1) .or. znew > Atm_zgrd(iznew)) &
       !     & call err_issue(1, 'mcarWld__rt_boundSup: iznew ='//num2str_AN(iznew)) !TEST
       call mcarWld__rt_sampSup(znew, iznew, path, isup)
       Pho_loc(3) = znew
       Pho_iz = iznew

       ! X-plane
       if (icase == 1) then
          Pho_loc(2) = Pho_loc(2) + path * Pho_dir(2)
          if (Pho_itr(1) == 0) then
             ixs = ixs - 1
             if (ixs <= 0) ixs = Atm_nxs(izs)
             Pho_loc(1) = Atm_xsup(ixs, izs)
          else
             ixs = ixs + 1
             if (ixs > Atm_nxs(izs)) ixs = 1
             Pho_loc(1) = Atm_xsup(ixs - 1, izs)
          end if

          ! Y-plane
       else
          Pho_loc(1) = Pho_loc(1) + path * Pho_dir(1)
          if (Pho_itr(2) == 0) then
             iys = iys - 1
             if (iys <= 0) iys = Atm_nys(izs)
             Pho_loc(2) = Atm_ysup(iys, izs)
          else
             iys = iys + 1
             if (iys > Atm_nys(izs)) iys = 1
             Pho_loc(2) = Atm_ysup(iys - 1, izs)
          end if
       end if
    end if

  end subroutine mcalWld__rt_moveSup


  !+
  ! Integrate path lengths and fluxes when moving in super-grid box
  !-
  subroutine mcarWld__rt_sampSup(znew, iznew, fpath, isup)

    real(R_), intent(in) :: fpath
    integer,  intent(in) :: isup
    integer,  intent(in) :: iznew
    real(R_), intent(in) :: znew
    integer   :: izl, izlb, izle, izls, izd
    real(R_)  :: auz, path, aloc(3)

    ! No flux sampling
    if (iznew == Pho_iz) then   ! the same layer
       call mcarPho__rt_hist_rec(Pho_iz, fpath)
       return
    end if
    auz = 1.0_R_ / Pho_dir(3)
    if (Wld_mdet(1) /= 1 .and. Wld_mdet(3) /= 1) then ! flux & radiance samplers are not active
       aloc(3) = Pho_loc(3)
       if (Pho_itr(3) == 0) then ! downward
          do izl = Pho_iz, iznew + 1, -1
             path = (Atm_zgrd(izl - 1) - aloc(3)) * auz
             call mcarPho__rt_hist_rec(izl, path)
             aloc(3) = Atm_zgrd(izl - 1)
          end do
       else ! upward
          do izl = Pho_iz, iznew - 1
             path = (Atm_zgrd(izl) - aloc(3)) * auz
             call mcarPho__rt_hist_rec(izl, path)
             aloc(3) = Atm_zgrd(izl)
          end do
       end if
       path = (znew - aloc(3)) * auz
       call mcarPho__rt_hist_rec(iznew, path)
       return
    end if

    ! Index parameters for layer nodes
    izd = 1 - Pho_itr(3) ! 1 for downward, 0 for upward
    izls = Pho_iz - izd
    izle = iznew - Pho_itr(3)
    if (Pho_itr(3) == 0) then
       izlb = -1
    else
       izlb = 1
    end if
    aloc(1:3) = Pho_loc(1:3)

    ! Move and sample at layer node points
    do izl = izls, izle, izlb
       path = (Atm_zgrd(izl) - aloc(3)) * auz
       call mcarPho__rt_hist_rec(izl + izd, path)
       aloc(3) = Atm_zgrd(izl)
       if (Pho_mv3D) then ! 3-D transfer?
          if (isup == 2) then ! 1-D super-layer
             call mcarAtm__rt_renewLocXY(Pho_dir, path, aloc) 
          else                ! 3-D super-layer
             aloc(1:2) = aloc(1:2) + path * Pho_dir(1:2)
          end if
       end if
       if (Wld_mdet(1) == 1) call mcarFlx__sampFlux(aloc(1), aloc(2), izl)
       if (Wld_mdet(3) == 1) call mcarRad__samp3(izl)
    end do

    ! The last path
    path = (znew - aloc(3)) * auz
    call mcarPho__rt_hist_rec(iznew, path)

  end subroutine mcarWld__rt_sampSup

end module mcarWld
