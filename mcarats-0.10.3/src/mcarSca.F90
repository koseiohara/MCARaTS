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
! Module for scattering
!-
module mcarSca 

  use globals
#if UseMPI == 1
  use globals_MPI
#endif
  use hparx
  use mcarUtl
  implicit none
  private

  ! Public procedures
  public :: mcarSca__tune_init
  public :: mcarSca__user_init
  public :: mcarSca__init
  public :: mcarSca__final
  public :: mcarSca__prep
  public :: mcarSca__qry_dims
  public :: mcarSca__prep_makeTab
  public :: mcarSca__asymStar_R
  public :: mcarSca__fracDelC_1R
  public :: mcarSca__fracDelC_R
  public :: mcarSca__newDirec
  public :: mcarSca__phaseFunc_R

  ! (Private) User variables in namelist
  character(256), save :: Sca_inpfile = ' ' ! file name for input of scattering properties
  integer,  save :: Sca_mfmt = 1     ! flag for input file format (0=text, 1=binary)
  integer,  save :: Sca_npf  = 1     ! # of tabulated phase functions
  integer,  save :: Sca_nskip = 0    ! # of data record lines to be skipped
  integer,  save :: Sca_nanci = 0    ! # of ancillary data
  integer,  save :: Sca_nangi = 1000 ! max # of angles for the input table
  integer,  save :: Sca_nchi = 4     ! # of orders for truncation approximation (>= 2)
  integer,  save :: Sca_ntg  = 20000 ! # of table grids for angles & probabilities
  real(R_), save :: Sca_qtfmax = 20.0_R_ ! geometrical truncation angle (deg.)
  namelist /mcarSca_nml_init/Sca_inpfile, Sca_mfmt, Sca_npf, Sca_nangi, Sca_nanci, &
       & Sca_nskip, Sca_nchi, Sca_qtfmax, Sca_ntg
  !// Note nchi should be >= 2.
  !   ichi=1 : no truncation approx; ichi=[2,nchi] : truncation approx.

  ! (Private)
  real(R_), save :: Sca_dsang, Sca_fsang
  real(R_), save, allocatable :: Sca_wtfTab(:)   !(nchi)
  real(R_), save, allocatable :: Sca_qtfTab(:)   !(nchi)
  real(R_), save, allocatable :: Sca_phsTab(:,:) !(ntg+1,npf)
  !// phsTab(1) for Q=0, phsTab(ntg) for Q=pi
  real(R_), save, allocatable :: Sca_scaTab(:,:) !(ntg+1,npf)
  !// scaTab(1) = 0, scatTab(ntg) = pi, increasing with itg = [1,ntg]
  real(R_), save, allocatable :: Sca_fdcTab(:,:) !(nchi,npf)
  real(R_), save, allocatable :: Sca_ftbTab(:,:) !(nchi,npf)
  real(R_), save, allocatable :: Sca_ptfTab(:,:) !(nchi,npf)
  real(R_), save, allocatable :: Sca_gtcTab(:,:) !(nchi,npf)
  real(R_), save, allocatable :: Sca_ang(:,:)    !(knai,npf)
  real(R_), save, allocatable :: Sca_phs(:,:)    !(knai,npf)
  integer,  save, allocatable :: Sca_naVec(:)   !(npf)
  !// Phase function parameter definition, apf
  !   ipf = int(apf), rpf = apf - ipf
  !   ipf = -2, apf = (-3,-2] : isotropic
  !   ipf = -1, apf = (-2,-1] : Rayleigh
  !   ipf =  0, apf = (-1,+1) : Henyey-Greenstein phase function, apf = g (asymmetry factor)
  !                             (Note apf should never be = 1)
  !   ipf >= 1, apf >= 1 : user-specified, tabulated phase functions, rpf = interpolation factor
  !                      P = (1 - rpf) * P_tab(ipf) + rpf * P_tab(ipf+1)

contains

  !+
  ! Tune a preset for technical parameters
  !-
  subroutine mcarSca__tune_init(moptim) 

    integer, intent(in) :: moptim ! optimization flag [-2,3]

    if (moptim <= -2) then ! no preset
       return
    else if (moptim == -1) then ! no optimization
       call mcarSca__set_tech(0.0_R_, 2, 40000)
    else if (moptim == 0) then ! unbiased optimization
       call mcarSca__set_tech(0.0_R_, 2, 40000)
    else if (moptim == 1) then ! conservative
       call mcarSca__set_tech(15.0_R_, 4, 40000)
    else if (moptim == 2) then ! standard
       call mcarSca__set_tech(30.0_R_, 4, 20000)
    else ! quick-and-dirty
       call mcarSca__set_tech(45.0_R_, 3, 2000)
    end if

  end subroutine mcarSca__tune_init


  !+
  ! Set technical parameters
  !-
  subroutine mcarSca__set_tech(qtfmax, nchi, ntg) 

    real(R_), intent(in) :: qtfmax
    integer,  intent(in) :: nchi, ntg
    Sca_qtfmax = qtfmax
    Sca_nchi = nchi
    Sca_ntg  = ntg
    
  end subroutine mcarSca__set_tech
  

  !+
  ! Read in namelist variables
  !-
  subroutine mcarSca__user_init(iu, filepath) 

    integer, intent(in) :: iu   ! file unit index for namelist input
    character(*), intent(in) :: filepath ! path name for I/O files
    integer :: ios

    read (iu, nml=mcarSca_nml_init, iostat=ios)
    call err_read(ios, iu, 'mcarSca_nml_init_user: Invalid namelist input for mcarSca_nml_init.')
    if (len_trim(Sca_inpfile) > 0) Sca_inpfile = trim(filepath)//'/'//Sca_inpfile

    call check_iI('mcarSca__user_init: Sca_nchi',   Sca_nchi, 2, 100)
    call check_iR('mcarSca__user_init: Sca_qtfmax', Sca_qtfmax, 0.0_R_, 80.0_R_)

  end subroutine mcarSca__user_init


  !+
  ! Initialize this module
  !-
  subroutine mcarSca__init() 

    Sca_dsang = PI_ / (Sca_ntg - 1)
    Sca_fsang = (Sca_ntg - 1) / PI_
    if (allocated(Sca_qtfTab)) call mcarSca__final()
    allocate (Sca_qtfTab(Sca_nchi))
    allocate (Sca_wtfTab(Sca_nchi))
    allocate (Sca_phsTab(Sca_ntg+1, Sca_npf+1))
    allocate (Sca_scaTab(Sca_ntg+1, Sca_npf+1))
    allocate (Sca_fdcTab(Sca_nchi, Sca_npf+1))
    allocate (Sca_ftbTab(Sca_nchi, Sca_npf+1))
    allocate (Sca_ptfTab(Sca_nchi, Sca_npf+1))
    allocate (Sca_gtcTab(Sca_nchi, Sca_npf+1))
    allocate (Sca_ang(Sca_nangi, Sca_npf))
    allocate (Sca_phs(Sca_nangi, Sca_npf))
    allocate (Sca_naVec(Sca_npf))
    !// Elements of Sca_*Tab arrays at ipf=Sca_npf+1 are margins, which are dummy.
    !   They are nessesary, however, because ocasional accesses due to interpolation w.r.t. ipf.

  end subroutine mcarSca__init


  !+
  ! Finalize this module
  !-
  subroutine mcarSca__final() 

    if (allocated(Sca_qtfTab)) then
       deallocate (Sca_qtfTab, Sca_wtfTab)
       deallocate (Sca_phsTab, Sca_scaTab, Sca_fdcTab, Sca_ftbTab, Sca_ptfTab)
       deallocate (Sca_gtcTab, Sca_ang, Sca_phs, Sca_naVec)
    end if

  end subroutine mcarSca__final


  !+
  ! Query about dimension sizes
  !-
  subroutine mcarSca__qry_dims(npf, nchi, ntg) 

    integer, intent(out), optional :: npf, nchi, ntg
    if (present(npf))  npf  = Sca_npf
    if (present(nchi)) nchi = Sca_nchi
    if (present(ntg))  ntg  = Sca_ntg

  end subroutine mcarSca__qry_dims


  !+
  ! Prepare the tables for scattering
  !-
  subroutine mcarSca__prep(irank, nrank, mbswap) 

    integer, intent(in) :: irank, nrank, mbswap
    integer :: iroot

    ! Read in
    iroot = 0
    if (irank == iroot) call mcarSca__read(mbswap)
#if UseMPI == 1
    call mcarSca__MPI_bcast_read(iroot)
#endif

    ! Make tables
    call mcarSca__prep_makeTab(irank, nrank)
#if UseMPI == 1
    call mcarSca__tab_bcast_MPI(irank, nrank)
#endif

  end subroutine mcarSca__prep


  !+
  ! Import data from file
  !-
  subroutine mcarSca__read(mbswap) 

    integer, intent(in) :: mbswap
    integer :: ipf, iai, ia, iud, nrec

    ! Import
    if (len_trim(Sca_inpfile) > 0) then
       iud = freeUnit_I(10)
       if (Sca_mfmt == 0) then ! text format
          call open_seq(iud, Sca_inpfile, 'old')
          call mcarUtl__txtSca_read(iud, Sca_nanci, Sca_naVec, Sca_ang, Sca_phs)
       else ! binary format
          nrec = recordLen_I(1, Sca_nangi + Sca_nanci, 4) ! Apr. 27, 2012
          call open_dir(iud, Sca_inpfile, nrec, 'old')
          call mcarUtl__binSca_read(iud, Sca_nangi, Sca_nanci, Sca_npf, Sca_ang(:,1), Sca_phs, &
               & Sca_nskip, mbswap)
          Sca_naVec(:) = Sca_nangi
          do ipf = 2, Sca_npf
             Sca_ang(:, ipf) = Sca_ang(:, 1) ! copy
          end do
       end if
       close (iud)
    end if

    ! Fix data
    do ipf = 1, Sca_npf
       Sca_ang(1:Sca_naVec(ipf), ipf) = Sca_ang(1:Sca_naVec(ipf), ipf) * DTOR_ ! degree to radian
       ia = 1
       do iai = 2, Sca_naVec(ipf)
          if (Sca_ang(iai, ipf) /= Sca_ang(iai - 1, ipf)) then ! collect unique data only
             ia = ia + 1
             Sca_ang(ia, ipf) = Sca_ang(iai, ipf)
             Sca_phs(ia, ipf) = Sca_phs(iai, ipf)
          end if
       end do
       Sca_naVec(ipf) = ia
       if (Sca_ang(1, ipf) > Sca_ang(ia, ipf)) then ! reverse table
          Sca_ang(1:ia, ipf) = Sca_ang(ia:1:-1, ipf)
          Sca_phs(1:ia, ipf) = Sca_phs(ia:1:-1, ipf)
       end if
    end do

  end subroutine mcarSca__read


#if UseMPI == 1
  !+
  ! MPI broadcast data read in mcarSca__read
  !-
  subroutine mcarSca__MPI_bcast_read(iroot) 

    include 'inc_mpi.f90'
    integer, intent(in) :: iroot ! root PE index
    integer :: ierr, n
    n = size(Sca_phs, 1) * Sca_npf
    call MPI_Bcast(Sca_naVec(:), Sca_npf, MPI_INTEGER, iroot, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(Sca_ang(:,:), n, MPI_R_, iroot, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(Sca_phs(:,:), n, MPI_R_, iroot, MPI_COMM_WORLD, ierr)

  end subroutine mcarSca__MPI_bcast_read
#endif


  !+
  ! Make LUTs for phase function & scattering angle
  !-
  subroutine mcarSca__prep_makeTab(irank, nrank)

    integer, intent(in) :: irank, nrank
    real(R_)  :: wrkphs(Sca_ntg), wrkang(Sca_ntg) !(AUTO)
    real(R_)  :: wrkcum(Sca_ntg), wrkpdf(Sca_ntg) 
    real(R_)  :: wrkcos(Sca_ntg), wrksin(Sca_ntg)
    real(R_)  :: wrkcosC(Sca_ntg), wrksca(Sca_ntg)
    integer   :: iitf(Sca_nchi)
    real(R_)  :: fbod, fdec, pfwd, utf
    real(R_)  :: g0bh, g0fh, g1a, g1bh, g1fh, g1t, g2a, g2bh, g2fh
    integer   :: ichi, ipf, itg

    ! Angle & cosine
    wrkang(:) = gridValues_1R(Sca_ntg, 0.0_R_, PI_)
    wrkcos(:) = cos(wrkang(:))
    wrksin(:) = sin(wrkang(:))
    do itg = 1, Sca_ntg
       if (wrkcos(itg) > 0.0_R_) then
          wrkcosC(itg) = wrksin(itg)**2 / (1.0_R_ + wrkcos(itg)) ! 1-cos
       else
          wrkcosC(itg) = wrksin(itg)**2 / (1.0_R_ - wrkcos(itg)) ! 1+cos
       end if
    end do

    ! Truncation angles
    Sca_qtfTab(:) = gridValues_1R(Sca_nchi, 0.0_R_, Sca_qtfmax * DTOR_)
    do ichi = 1, Sca_nchi
       iitf(ichi) = gridIdx_bin_I(wrkang, 1, Sca_ntg, Sca_qtfTab(ichi))
    end do
    Sca_qtfTab(:) = wrkang(iitf(:))
    Sca_wtfTab(:) = 1.0_R_ - cos(Sca_qtfTab(:))

    ! Loop for phase functions
    do ipf = irank + 1, Sca_npf, nrank

       ! Make tables
       call scatPF_intp(Sca_ang(:, ipf), Sca_phs(:, ipf), Sca_naVec(ipf), wrkang, Sca_ntg, wrkphs)
       call scatPF_distFunc(wrkcos, wrkcosC, wrkphs, wrkpdf, wrkcum)
       call mcarSca__invertTab(wrkang, wrkcum, wrksca)
       Sca_phsTab(1:Sca_ntg, ipf) = wrkphs(1:Sca_ntg)
       Sca_scaTab(1:Sca_ntg, ipf) = wrksca(1:Sca_ntg)
       Sca_phsTab(Sca_ntg+1, ipf) = Sca_phsTab(Sca_ntg, ipf)
       Sca_scaTab(Sca_ntg+1, ipf) = Sca_scaTab(Sca_ntg, ipf)

       ! Truncation
       do ichi = 1, Sca_nchi
          if (Sca_qtfTab(ichi) <= 0.0_R_ .or. ichi == 1) then
             fdec = 1.0_R_
             fbod = 1.0_R_
             pfwd = 1.0_R_
             g1t = scatPF_asym_R(wrkcos, wrksin, wrkphs) 
          else
             utf = cos(Sca_qtfTab(ichi))
             call scatPF_truncDFlat(wrkphs, wrkcos, wrkcosC, iitf(ichi), fdec, fbod, pfwd)
             call scatPF_biAsym(wrkcos, wrksin, wrkphs, utf, g0fh, g0bh, g1fh, g1bh, g1a, &
                  & g2fh, g2bh, g2a)
             g1t = ((fdec - fbod) * (1.0_R_ + utf) * 0.5_R_ + fbod * g1bh) / fdec
          end if
          Sca_fdcTab(ichi, ipf) = fdec
          Sca_ftbTab(ichi, ipf) = fbod
          Sca_ptfTab(ichi, ipf) = pfwd
          Sca_gtcTab(ichi, ipf) = g1t
          !if (ichi == Sca_nchi) print *, ipf, fdec, fbod, pfwd !TEST
       end do
    end do

  end subroutine mcarSca__prep_makeTab


#if UseMPI == 1
  !+
  ! Make LUTs for phase function & scattering angle
  !-
  subroutine mcarSca__tab_bcast_MPI(irank, nrank)

    include 'inc_mpi.f90'
    integer, intent(in) :: irank, nrank
    integer   :: ichi, ipf, ierr, iroot, n
    real(R_)  :: rbuf(Sca_nchi*4), sbuf(Sca_ntg + 1), pbuf(Sca_ntg + 1) !(AUTO)

    ! MPI broadcast
    do ipf = 1, Sca_npf
       iroot = mod(ipf - 1, nrank)
       if (irank == iroot) then ! work data
          do ichi = 1, Sca_nchi
             rbuf(ichi)                = Sca_fdcTab(ichi, ipf)
             rbuf(ichi + Sca_nchi)     = Sca_ftbTab(ichi, ipf)
             rbuf(ichi + Sca_nchi * 2) = Sca_ptfTab(ichi, ipf)
             rbuf(ichi + Sca_nchi * 3) = Sca_gtcTab(ichi, ipf)
          end do
          pbuf(1:Sca_ntg+1) = Sca_phsTab(1:Sca_ntg+1, ipf)
          sbuf(1:Sca_ntg+1) = Sca_scaTab(1:Sca_ntg+1, ipf)
       end if
       n = Sca_nchi * 4
       call MPI_Bcast(rbuf, n, MPI_R_, iroot, MPI_COMM_WORLD, ierr)
       n = Sca_ntg + 1
       call MPI_Bcast(pbuf, n, MPI_R_, iroot, MPI_COMM_WORLD, ierr)
       call MPI_Bcast(sbuf, n, MPI_R_, iroot, MPI_COMM_WORLD, ierr)
       do ichi = 1, Sca_nchi
          Sca_fdcTab(ichi, ipf) = rbuf(ichi)
          Sca_ftbTab(ichi, ipf) = rbuf(ichi + Sca_nchi)
          Sca_ptfTab(ichi, ipf) = rbuf(ichi + Sca_nchi * 2)
          Sca_gtcTab(ichi, ipf) = rbuf(ichi + Sca_nchi * 3)
       end do
       Sca_phsTab(1:n, ipf) = pbuf(1:n)
       Sca_scaTab(1:n, ipf) = sbuf(1:n)
    end do

  end subroutine mcarSca__tab_bcast_MPI
#endif
  

  !+
  ! Make LUTs of scattering angle for uniform distributions
  !-
  subroutine mcarSca__invertTab(wrkang, wrkcum, wrksca)

    real(R_), intent(in)  :: wrkang(:) ! angles (radian)
    real(R_), intent(in)  :: wrkcum(:) ! CDF, 1 to 0, decreasing
    real(R_), intent(out) :: wrksca(:) ! angles for equally-spaced probabilities
    !// wrksca(1) = 0, wrksca(ntg) = pi, increasing with itg = [1,ntg]
    integer   :: ntg, itg, iwrk
    real(R_)  :: p, pdlt, rat, dx

    ! Scattering angle inversion
    ntg = size(wrksca)
    pdlt = 1.0_R_ / (ntg - 1)
    iwrk = 1
    do itg = 2, ntg - 1
       p = pdlt * real(ntg - itg, R_) ! decrease with increaing itg
       iwrk = gridIdx_seq_I(wrkcum, p, iwrk, iwrk, size(wrkcum))
       dx = wrkcum(iwrk) - wrkcum(iwrk + 1)
       if (dx > REPS_) then
          rat = (p - wrkcum(iwrk + 1)) / dx
          wrksca(itg) = rat * wrkang(iwrk) + (1.0_R_ - rat) * wrkang(iwrk + 1)
       else
          wrksca(itg) = wrkang(iwrk)
       end if
    end do
    wrksca(1)   = 0.0_R_
    wrksca(ntg) = PI_

  end subroutine mcarSca__invertTab


  !+
  ! Asymmetry factor, g*, for truncated phase function
  !-
  function mcarSca__asymStar_R(apf, ichi) result(g)

    integer,  intent(in) :: ichi
    real(R_), intent(in) :: apf
    real(R_) :: g
    real(R_) :: rpf, fdec, wtf, b, c, r, s0, s1
    integer  :: ipf
    
    if (apf >= 1.0_R_) then ! anisotropic, tabulated
       ipf = int(apf)
       rpf = apf - ipf
       g = Sca_gtcTab(ichi, ipf)
       if (rpf > REPS_) g = g + rpf * (Sca_gtcTab(ichi, ipf+1) - g)
    else if (apf > -1.0_R_) then ! H-G
       wtf = Sca_wtfTab(ichi)
       g = apf
       if (apf > 0.0_R_ .and. wtf > 0.0_R_) then
          b = 1.0_R_ + g**2
          r = sqrt(b - 2.0_R_ * g * (1.0_R_ - wtf))
          c = (g + r - 1.0_R_) * (1.0_R_ - g**2) / (2.0_R_ * g)
          s0 = c / (r * (1.0_R_ - g)) ! g should not be = 1
          s1 = (b * s0 - c) / (2.0_R_ * g)
          fdec = 1.0_R_ - (2.0_R_ * s1 - (2.0_R_ - wtf) * s0) / wtf
          g = 1.0_R_ - (1.0_R_ - g) / fdec ! g*
       end if
    else ! Rayleigh or isotropic
       g = 0.0_R_
    end if

  end function mcarSca__asymStar_R


  !+
  ! Complementary fraction (1 - f) of delta phase function, a vector version
  !-
  function mcarSca__fracDelC_1R(apf, ichi) result(fdec)

    integer,  intent(in) :: ichi
    real(R_), intent(in) :: apf(:)
    real(R_) :: fdec(size(apf))
    real(R_) :: r, g, b, c, s0, s1, wtf
    integer  :: i, ipf

    wtf = Sca_wtfTab(ichi)
    if (wtf <= 0.0_R_) then ! no truncation
       fdec(:) = 1.0_R_
    else
       do i = 1, size(apf)
          if (apf(i) <= 0.0_R_) then ! no truncation
             fdec(i) = 1.0_R_
          else if (apf(i) >= 1.0_R_) then ! tabulated phase function
             ipf = int(apf(i))
             r = apf(i) - ipf
             fdec(i) = (1.0_R_ - r) * Sca_fdcTab(ichi, ipf) + r * Sca_fdcTab(ichi, ipf + 1)
          else ! H-G
             g = apf(i)
             b = 1.0_R_ + g**2
             r = sqrt(b - 2.0_R_ * g * (1.0_R_ - wtf))
             c = (g + r - 1.0_R_) * (1.0_R_ - g**2) / (2.0_R_ * g)
             s0 = c / (r * (1.0_R_ - g)) ! g should not be = 1
             s1 = (b * s0 - c) / (2.0_R_ * g)
             fdec(i) = 1.0_R_ - (2.0_R_ * s1 - (2.0_R_ - wtf) * s0) / wtf
          end if
       end do
    end if

  end function mcarSca__fracDelC_1R


  !+
  ! Complementary fraction (1 - f) of delta phase function
  !-
  function mcarSca__fracDelC_R(apf, ichi) result(fdec)

    integer,  intent(in) :: ichi
    real(R_), intent(in) :: apf
    real(R_) :: fdec
    real(R_) :: r, g, b, c, s0, s1, wtf
    integer  :: ipf

    wtf = Sca_wtfTab(ichi)
    if (apf <= 0.0_R_ .or. wtf <= 0.0_R_) then ! no truncation
       fdec = 1.0_R_
    else if (apf >= 1.0_R_) then ! tabulated phase function
       ipf = int(apf)
       r = apf - ipf
       fdec = (1.0_R_ - r) * Sca_fdcTab(ichi, ipf) + r * Sca_fdcTab(ichi, ipf + 1)
    else ! H-G
       g = apf
       b = 1.0_R_ + g**2
       r = sqrt(b - 2.0_R_ * g * (1.0_R_ - wtf))
       c = (g + r - 1.0_R_) * (1.0_R_ - g**2) / (2.0_R_ * g)
       s0 = c / (r * (1.0_R_ - g)) ! g should not be = 1
       s1 = (b * s0 - c) / (2.0_R_ * g)
       fdec = 1.0_R_ - (2.0_R_ * s1 - (2.0_R_ - wtf) * s0) / wtf
    end if

  end function mcarSca__fracDelC_R


  !+
  ! New direction due to scattering by an atmospheric particle
  !-
  subroutine mcarSca__newDirec(apf, ichi, dir, pps)

    real(R_), intent(in) :: apf  ! phase function specification parameter
    integer,  intent(in) :: ichi ! truncation order
    real(R_), intent(inout) :: dir(:) ! direction vector
    real(R_), intent(out) :: pps  ! probability density function per solid angle
    real(R_),  parameter :: API4 = 0.25_R_ / PI_
    real(R_)  :: rpf, ftb, fdc, cosf, r2, sinf, wtf, p, cosqC, cosq, sinq
    integer   :: ipf

    ! Taburated, anisotropic scattering
    ipf = int(apf)
    if (ipf >= 1) then
       call rand_point_circle(1.0e-12_R_, r2, sinf, cosf) ! a random point in a unit circle
       wtf = Sca_wtfTab(ichi)
       rpf = apf - ipf ! index interpolation factor
       if (wtf <= 0.0_R_) then ! no truncation
          call mcarSca__randScat_uTab(r2, ipf, rpf, sinq, cosq, p)
          pps = API4 * p
       else ! truncated
          ftb = Sca_ftbTab(ichi, ipf)
          fdc = Sca_fdcTab(ichi, ipf)
          if (rpf > REPS_) then
             ftb = ftb + rpf * (Sca_ftbTab(ichi, ipf + 1) - ftb)
             fdc = fdc + rpf * (Sca_fdcTab(ichi, ipf + 1) - fdc)
          end if
          if (r2 * fdc <= ftb) then ! body part
             r2 = 1.0_R_ - r2 * fdc ! renormalized in (1-ftb, 1]
             call mcarSca__randScat_uTab(r2, ipf, rpf, sinq, cosq, p)
          else ! flat, forward part
             cosqC = Sca_wtfTab(ichi) * mseq_rand_R() ! 1 - cos, accurate
             cosq = 1.0_R_ - cosqC ! not accurate
             sinq = sqrt(cosqC * (2 - cosqC)) ! accurate
             p = (1.0_R_ - rpf) * Sca_ptfTab(ichi, ipf) + rpf * Sca_ptfTab(ichi, ipf + 1)
          end if
          pps = API4 * p / fdc
       end if
       dir(:) = geo3D_rotateU_1R(dir, sinq, cosq, sinf, cosf)

       ! H-G
    else if (ipf == 0) then
       call rand_point_circle(1.0e-12_R_, r2, sinf, cosf) ! a random point in a unit circle
       wtf = Sca_wtfTab(ichi)
       if (wtf <= 0.0_R_ .or. apf <= 0.0_R_) then ! no truncation
          call randScat_HG(apf, r2, cosq, sinq)
          pps = API4 * scatPF_HG_R(apf, cosq)
       else ! truncated
          call scatPF_truncDFlat_HG(apf, wtf, fdc, ftb, p)
          !if (r2 * fdc <= ftb) then ! body part
          !   r2 = r2 * fdc ! renormalized in [0,ftb)
          if (r2 * fdc <= fdc - ftb) then ! body part
             r2 = 1.0_R_ - fdc + r2 * fdc ! renormalized in (1-ftb,1]
             call randScat_HG(apf, r2, cosq, sinq)
             p = scatPF_HG_R(apf, cosq)
          else ! flat, forward part
             cosqC = wtf * mseq_rand_R() ! 1 - cos, accurate
             cosq = 1.0_R_ - cosqC ! not accurate
             sinq = sqrt(cosqC * (2 - cosqC)) ! accurate
          end if
          pps = API4 * p / fdc
       end if

       !if (cosq >= -1.0_R_ .and. cosq <= 1.0_R_) then
       !else
       !   print *, 'Sca:', sinq, cosq, sinf, cosf
       !   print *, r2, apf, wtf, ichi
       !endif

       dir(:) = geo3D_rotateU_1R(dir, sinq, cosq, sinf, cosf)

       ! Rayleigh scattering
    else if (ipf == -1) then
       call rand_point_circle(1.0e-12_R_, r2, sinf, cosf) ! a random point in a unit circle
       call rayl_randScat(sinq, cosq)
       pps = API4 * rayl_scatPF_R(cosq)
       dir(:) = geo3D_rotateU_1R(dir, sinq, cosq, sinf, cosf)

       ! Isotropic scattering
    else
       call rand_point_sphere(1.0e-12_R_, r2, dir(1), dir(2), dir(3))
       pps = API4
    end if
    dir(3) = nonZero_R(dir(3), REPS_) ! non-zero to avoid horizontal directions

  end subroutine mcarSca__newDirec


  !+
  ! Determine a random scattering angle
  !-
  subroutine mcarSca__randScat_uTab(xi, ipf, rpf, sinq, cosq, p)

    real(R_), intent(in) :: xi ! a random number
    !// For xi = 0, Q = 0  and cosq = +1 (forward scattering)
    !   For xi = 1, Q = pi and cosq = -1 (backscattering)
    integer,  intent(in) :: ipf
    real(R_), intent(in) :: rpf
    real(R_), intent(out) :: sinq, cosq
    real(R_), intent(out) :: p
    real(R_) :: ritg, rtg, q0, q1, q, p0, p1
    integer :: itg

    ritg = (Sca_ntg - 1) * xi + 1.0_R_
    itg = int(ritg)
    rtg = ritg - itg
    q0 = Sca_scaTab(itg,   ipf)
    q1 = Sca_scaTab(itg+1, ipf)
    p0 = Sca_phsTab(itg,   ipf)
    p1 = Sca_phsTab(itg+1, ipf)
    if (rpf > REPS_) then
       q0 = q0 + rpf * (Sca_scaTab(itg,   ipf+1) - q0)
       q1 = q1 + rpf * (Sca_scaTab(itg+1, ipf+1) - q1)
       p0 = p0 + rpf * (Sca_phsTab(itg,   ipf+1) - p0)
       p1 = p1 + rpf * (Sca_phsTab(itg+1, ipf+1) - p1)
    end if
    q = (1.0_R_ - rtg) * q0 + rtg * q1 ! radian
    p = (1.0_R_ - rtg) * p0 + rtg * p1
    call mtab_sincos(q, sinq, cosq)

  end subroutine mcarSca__randScat_uTab


  !+
  ! Phase function of atmospheric particles
  !-
  function mcarSca__phaseFunc_R(q, ipf, rpf, ichi) result(res)

    real(R_), intent(in) :: q    ! scattering angle (radian)
    integer,  intent(in) :: ipf  ! phase function index (should be >= 1)
    real(R_), intent(in) :: rpf  ! phase function interpolation factor
    integer,  intent(in) :: ichi
    real(R_)  :: res
    integer   :: itg
    real(R_)  :: rtg, ritg, p0, p1, fdc

    ! Flat, forward part
    if (q < Sca_qtfTab(ichi)) then
       res = Sca_ptfTab(ichi, ipf)
       fdc = Sca_fdcTab(ichi, ipf)
       if (rpf > REPS_) then
          res = res + rpf * (Sca_ptfTab(ichi, ipf+1) - res)
          fdc = fdc + rpf * (Sca_fdcTab(ichi, ipf+1) - fdc)
       end if
       res = res / fdc ! normalized truncated phase function

       ! Body part of tabulated phase function
    else
       ritg = Sca_fsang * q + 1.0_R_
       itg = int(ritg)
       rtg = ritg - itg
       p0  = Sca_phsTab(itg,   ipf)
       p1  = Sca_phsTab(itg+1, ipf)
       fdc = Sca_fdcTab(ichi,  ipf)
       if (rpf > REPS_) then
          p0  = p0  + rpf * (Sca_phsTab(itg,   ipf+1) - p0)
          p1  = p1  + rpf * (Sca_phsTab(itg+1, ipf+1) - p1)
          fdc = fdc + rpf * (Sca_fdcTab(ichi,  ipf+1) - fdc)
       end if
       res = (1.0_R_ - rtg) * p0 + rtg * p1
       res = res / fdc ! normalized truncated phase function
    end if

  end function mcarSca__phaseFunc_R

end module mcarSca
