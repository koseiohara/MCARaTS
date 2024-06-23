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
! Library of random number utilities
!-
module hparx_rand

  use globals, only : R_, RD_, I4_, I8_, REPS_, RSML_, IMAX_
  use hparx_base, only : shuffle_m1I, shuffle_m1A, cyclicAdd_I
  implicit none
  private

  ! Public
  public :: crypt_m1I
  public :: crypt_m1A
  public :: mseq_seed_I
  public :: mseq_init
  public :: mseq_shuffle
  public :: mseq_rand_1R
  public :: mseq_rand_1RD
  public :: mseq_rand_R
  public :: mseq_rand_RD
  public :: mseq_rand_I
  public :: rand_exp_R
  public :: rand_Gauss_1R
  public :: rand_isRare_L
  public :: rand_point_circle
  public :: rand_point_sphere
  public :: rand_newSeed_I
  public :: rand_seed_I
  public :: randStk_alloc
  public :: randStk_dealloc
  public :: randStk_setAuto
  public :: randStk_setUser
  public :: randStk_rand_R
  public :: randStk_rand_1R
  public :: russianRoulette
  public :: stat_approxStat

  ! M-sequence random number generators
  integer,  parameter :: Mseq_NBUF = 521      ! table size
  integer(I4_), save :: Mseq_jbuf(Mseq_NBUF)  ! the table
  integer,  save :: Mseq_ibuf = Mseq_NBUF + 1 ! pointer
  integer,  save :: Mseq_istat = 0            ! status of initialization (0 = not yet)

  ! Random number stack
  integer,  save :: Stk_nummax   ! max # of elements in chache
  integer,  save :: Stk_num   ! # of active values
  real(R_), save, allocatable :: Stk_rand(:)  ! chache for random numbers

contains

  !+
  ! Encrypt/Decrypt an integer vector data
  !-
  subroutine crypt_m1I(jbuf, iseed, menc) 

    integer,  intent(inout) :: jbuf(:)
    integer,  intent(in)    :: iseed
    integer,  intent(in)    :: menc     ! 1:enc, -1:dec
    integer,  allocatable :: jran(:)
    integer  :: itot, ntot

    ! Initialize
    ntot = size(jbuf, 1)
    allocate (jran(ntot))
    call mseq_init(iseed)
    do itot = 1, ntot
       jran(itot) = mseq_rand_I(1, ntot) ! store random numbers
    end do

    ! Encrypt/Decrypt
    if (menc == 1) then ! encrypt
       call shuffle_m1I(jbuf, jran, 1, ntot, 1)
       do itot = 1, ntot
          jbuf(itot) = cyclicAdd_I(jbuf(itot), mseq_rand_I(1, IMAX_))
       end do
    else                ! decrypt
       do itot = 1, ntot
          jbuf(itot) = cyclicAdd_I(jbuf(itot), -mseq_rand_I(1, IMAX_))
       end do
       call shuffle_m1I(jbuf, jran, ntot, 1, -1)
    end if
    deallocate (jran)

  end subroutine crypt_m1I


  !+
  ! Encrypt/Decrypt a character vector data
  !-
  subroutine crypt_m1A(abuf, iseed, menc)

    character(1), intent(inout) :: abuf(:)
    integer,  intent(in)  :: iseed
    integer,  intent(in)  :: menc     ! 1:enc, -1:dec
    integer,  allocatable :: jwrk(:)
    integer  :: ntot, nc, ni, ii

    ! Initialize
    ntot = size(abuf, 1)
    ni = ntot / 4 ! = size(transfer(abuf, jwrk))
    nc = ni * 4
    allocate (jwrk(ni))
    call mseq_init(iseed)
    do ii = 1, ni
       jwrk(ii) = mseq_rand_I(1, ntot) ! store random numbers in [1,ntot]
    end do

    ! Encrypt
    if (menc == 1) then
       call shuffle_m1A(abuf, jwrk, 4, nc, 4)  ! forward order
       jwrk(1:ni) = transfer(abuf(1:nc), jwrk)   ! A to I
       do ii = 1, ni
          jwrk(ii) = cyclicAdd_I(jwrk(ii), mseq_rand_I(1, IMAX_)) ! add positive
       end do
       abuf(1:nc) = transfer(jwrk(1:ni), abuf)   ! I to A

       ! Decrypt
    else
       jwrk(1:ni) = transfer(abuf(1:nc), jwrk)   ! A to I
       do ii = 1, ni
          jwrk(ii) = cyclicAdd_I(jwrk(ii), -mseq_rand_I(1, IMAX_)) ! add negative
       end do
       abuf(1:nc) = transfer(jwrk(1:ni), abuf)   ! I to A
       call mseq_init(iseed)                     ! initialize, again!
       do ii = 1, ni
          jwrk(ii) = mseq_rand_I(1, ntot)
       end do
       call shuffle_m1A(abuf, jwrk, nc, 4, -4) ! reverse order
    end if
    deallocate (jwrk)

  end subroutine crypt_m1A


  function mseq_seed_I(iseq) result(iseed)

    integer,  intent(in), optional :: iseq   ! my sequence index in [0,?] (not so large)
    integer  :: iseed  ! a seed generated in [1,IMAX_]
    write (*,*) 'mseq_seed_I: This is obsolete. Use rand_seed_I instead.'
    iseed = rand_seed_I(iseq)

  end function mseq_seed_I


  !+
  ! Initialize MLRS random number generator
  !
  ! References
  !   Fushimi, M., 1989: Random number. Tokyo Daigaku Shuppan Kai.
  !   Press et al., Numerical Recipes in Fortran 77._R_ Second Edition.
  !      Cambridge University Press.
  !-
  subroutine mseq_init(iseed)

    integer,  intent(in) :: iseed    ! a seed [1, 2**31-1]
    integer,  parameter :: JA = 16807, JM = IMAX_, JH = JM / 2 + 1
    integer,  parameter :: JQ = 127773, JR = 2836
    integer  :: ia(Mseq_NBUF), i, i1, i2, ibuf, ipass, npass
    integer(I4_) :: mtmp

    ! Validate the seed
    if (iseed >= 0) then ! 0 or positive
       i = iseed
    else                 ! negative
       i = iseed + JM + 1
    end if
    i = max(1, i) ! now corrected in valid range [1, 2**31-1]

    ! Warm-up the sequence
    do ipass = 1, 100      ! i=mod(JA*i,JM) without overflows by Schrage's method
       mtmp = i / JQ
       i = JA * (i - JQ * mtmp) - JR * mtmp
       if (i < 0) i = i + JM
    end do
    npass = mod(i, 100)    ! semi-randomly chosen # of passes
    do ipass = 1, npass
       mtmp = i / JQ
       i = JA * (i - JQ * mtmp) - JR * mtmp
       if (i < 0) i = i + JM
    end do

    ! Initial 521 bits
    do ibuf = 1, Mseq_NBUF
       mtmp = i / JQ
       i = JA * (i - JQ * mtmp) - JR * mtmp
       if (i < 0) i = i + JM
       if (i >= JH) then    ! 0 or 1 semi-randomly
          ia(ibuf) = 1
       else
          ia(ibuf) = 0
       end if
    end do

    ! Buffer array with random integer values
    i1 = 1
    i2 = 490
    do ibuf = 1, Mseq_NBUF
       mtmp = 0
       do i = 0, 31
          if (i <= 30) mtmp = mtmp * 2 + ia(i1)
          ia(i1) = mod(ia(i1) + ia(i2), 2)
          i1 = i1 + 1
          i2 = i2 + 1
          if (i1 > Mseq_NBUF) i1 = 1
          if (i2 > Mseq_NBUF) i2 = 1
       end do
       Mseq_jbuf(ibuf) = mtmp
    end do
    Mseq_istat = 1 ! initialization finished

    ! A new sequence (could be multiple passes)
    call mseq_newSeq()

  end subroutine mseq_init


  !+
  ! Check status of M-sequence and possibly generate a new sequence
  !-
  subroutine mseq_checkNewSeq() 

    if (Mseq_ibuf > Mseq_NBUF) then
       if (Mseq_istat == 0) call mseq_init(rand_seed_I()) ! safety net
       call mseq_newSeq()
    end if

  end subroutine mseq_checkNewSeq


  !+
  ! Generate a new sequence of uniform random numbers between [0,1)
  !  by Maximum-length Linearly-Recurring Sequence (MLRS) method.
  !
  ! Algorithm: Fushimi, M. (1989): The random number. 
  !   UP Ohyo Sugaku Sensho 12, Tokyo University Press (in Japanese).
  ! Note: Call mseq_init before use of this subroutine!
  !-
  subroutine mseq_newSeq()

    integer  :: ibuf
    do ibuf = 1, 32
       Mseq_jbuf(ibuf) = ieor(Mseq_jbuf(ibuf), Mseq_jbuf(ibuf - 32 + Mseq_NBUF))
    end do
    do ibuf = 33, Mseq_NBUF
       Mseq_jbuf(ibuf) = ieor(Mseq_jbuf(ibuf), Mseq_jbuf(ibuf - 32))
    end do
    Mseq_ibuf = 1  ! reset the pointer

  end subroutine mseq_newSeq


  !+
  ! Suffle the present sequence
  !-
  subroutine mseq_shuffle(ns)

    integer,  intent(in) :: ns
    integer  :: is, ibuf1, ibuf2
    integer(I4_) :: jbuf1
    do is = 1, ns
       ibuf1 = min(Mseq_NBUF, int(Mseq_NBUF * mseq_rand_R()) + 1)
       ibuf2 = min(Mseq_NBUF, int(Mseq_NBUF * mseq_rand_R()) + 1)
       jbuf1 = Mseq_jbuf(ibuf1)
       Mseq_jbuf(ibuf1) = Mseq_jbuf(ibuf2)
       Mseq_jbuf(ibuf2) = jbuf1
    end do

  end subroutine mseq_shuffle


  !+
  ! Generate uniform random numbers
  !-
  function mseq_rand_1R(n) result(rn)

    integer, intent(in) :: n ! # of random numbers
    real(R_) :: rn(n) ! uniform random numbers in [0.0_R_, 1.0_R_] 
    real(RD_), parameter :: FNORM = 1.0_RD_ / IMAX_ ! 1/2.147483647e+9_RD_
    integer :: i

    if (Mseq_ibuf + n - 1 > Mseq_NBUF) then
       do i = 1, n
          rn(i) = mseq_rand_RD()
       end do
    else
       rn(1:n) = Mseq_jbuf(Mseq_ibuf : Mseq_ibuf + n - 1) * FNORM
       Mseq_ibuf = Mseq_ibuf + n
    end if

  end function mseq_rand_1R


  !+
  ! Generate uniform random numbers
  !-
  function mseq_rand_1RD(n) result(rn)

    integer, intent(in) :: n ! # of random numbers
    real(RD_) :: rn(n) ! uniform random numbers in [0.0_RD_, 1.0_RD_] 
    real(RD_), parameter :: FNORM = 1.0_RD_ / IMAX_ ! 1/2.147483647e+9_RD_
    integer :: i

    if (Mseq_ibuf + n - 1 > Mseq_NBUF) then
       do i = 1, n
          rn(i) = mseq_rand_RD()
       end do
    else
       rn(1:n) = Mseq_jbuf(Mseq_ibuf : Mseq_ibuf + n - 1) * FNORM
       Mseq_ibuf = Mseq_ibuf + n
    end if

  end function mseq_rand_1RD


  !+
  ! Generate a uniform random number
  !-
  function mseq_rand_R() result(rn)

    real(R_) :: rn
    real(RD_), parameter :: FNORM = 1.0_RD_ / IMAX_ ! 1/2.147483647e+9_RD_
    call mseq_checkNewSeq()
    rn = Mseq_jbuf(Mseq_ibuf) * FNORM  ! a random number
    Mseq_ibuf = Mseq_ibuf + 1  ! next pointer

  end function mseq_rand_R


  !+
  ! Generate a uniform random number
  !-
  function mseq_rand_RD() result(rn)

    real(RD_) :: rn
    real(RD_), parameter :: FNORM = 1.0_RD_ / IMAX_ ! 1/2.147483647e+9_RD_
    call mseq_checkNewSeq()
    rn = Mseq_jbuf(Mseq_ibuf) * FNORM  ! a random number
    Mseq_ibuf = Mseq_ibuf + 1  ! next pointer

  end function mseq_rand_RD


  !+
  ! Generate an integer, uniform random number
  !-
  function mseq_rand_I(jmin, jmax) result(j)

    integer,  intent(in) :: jmin, jmax
    integer  :: j ! random number between [jmin, jmax]
    call mseq_checkNewSeq()
    j = int( int(Mseq_jbuf(Mseq_ibuf), I8_) * int(jmax - jmin + 1, I8_) / (int(IMAX_, I8_) + 1), I4_)
    j = jmin + j
    Mseq_ibuf = Mseq_ibuf + 1  ! next pointer

  end function mseq_rand_I


  !+
  ! Randomly-determined variable that obeys to the standard exponetial distribution
  !-
  function rand_exp_R() result(res)

    real(R_) :: res
    real(R_) :: rn

    rn = mseq_rand_R()
    if (rn > RSML_) then
       res = -log(rn)
    else
       res = -log(RSML_)
    end if

  end function rand_exp_R


  !+ Generate a random Gaussian vector data with given mean & sigma values
  !  with/without corrections for exact mean & sigma
  !-
  function rand_Gauss_1R(n, mcor, ave, sig) result(dat)

    integer,  intent(in) :: n    ! # of requested data points
    integer,  intent(in) :: mcor ! flag for correction of mean & sigma (1:yes) 
    real(R_), intent(in), optional :: ave, sig ! average & sigma (possibly not 0 & 1, respectively)
    real(R_) :: dat(n) ! random numbers that obeys Gaussian
    real(R_) :: r, x, y, rn, ave1, sig1
    integer  :: ii

    ! Random X value
    do ii = 1, (n + 1) / 2
       call rand_point_circle(REPS_**2, rn, x, y)
       r = sqrt(-2.0_R_ * log(rn))
       dat(2 * ii - 1) = r * x
       if (2 * ii <= n) dat(2 * ii) = r * y
    end do

    ! Correct mean & sigma
    if (mcor == 1) then
       ave1 = sum(dat) / n ! should be =~ 0
       sig1 = sqrt(max(0.0_R_, sum(dat**2) / n - ave1**2)) ! should be =~ 1
       if (sig1**2 > RSML_) then ! usual case
          dat(:) = (dat(:) - ave1) / sig1
       else ! unexpected case
          dat(:) = dat(:) - ave1
       end if
    end if

    ! Set mean & sigma
    if (present(ave) .and. present(sig)) dat(:) = dat(:) * sig + ave

  end function rand_Gauss_1R


  !+
  ! Determines whether a rare event occurs
  !-
  function rand_isRare_L(psmall) result(is)

    real(R_), intent(in) :: psmall ! a small probability for rare event
    logical :: is ! result

    ! Not very small psmall
    if (psmall > 1.0e-6_R_) then
       if (mseq_rand_R() < psmall) then
          is = .true.
       else
          is = .false.
       end if

       ! Very small psmall
    else
       if (mseq_rand_R() >= 1.0e-6_R_) then
          is = .false.
       else if (mseq_rand_R() < psmall * 1.0e+6_R_) then
          is = .true.
       else
          is = .false.
       end if
    end if

  end function rand_isRare_L


  !+
  ! Automatically generate a seed for the random number generator
  !-
  function rand_newSeed_I() result(iseed)

    integer  :: iseed  ! a seed generated in [1,IMAX_]
    integer  :: i, idt(8)
    call date_and_time(values=idt)
    i = 24 * (60 * idt(7) + idt(6)) + idt(5)   ! [0,86400]
    iseed = 29 * (12 * (31 * i + idt(3)) + idt(2)) ! 29*[13,32140831]

  end function rand_newSeed_I


  !+
  ! Initialize a seed for the random number generator
  !   Returns a seed possibly for one of several sequences for parallel computing
  !-
  function rand_seed_I(iseq, iseed0) result(iseed)

    integer,  intent(in), optional :: iseq   ! my sequence index in [0,?] (not so large)
    integer,  intent(in), optional :: iseed0 ! a candidate seed
    integer  :: iseed  ! a seed generated in [1,IMAX_]
    integer  :: i

    if (present(iseed0)) then
       i = iseed0
       if (i <= 0) i = i + IMAX_
       if (i <= 0) i = i + IMAX_       
    else
       i = rand_newSeed_I()
    end if

    if (present(iseq)) then
       i = cyclicAdd_I(i, iseq)
       if (i <= 0) i = i + IMAX_
       if (i <= 0) i = i + IMAX_
    end if
    iseed = i

  end function rand_seed_I


  !+
  ! Determines a random point on the unit circle
  !  This can be used to get a random unit vector (x,y) and a uniform random number
  !-
  subroutine rand_point_circle(r2min, rn, x, y)

    real(R_), intent(in)  :: r2min  ! (min radius)^2
    real(R_), intent(out) :: rn     ! a uniform random number   
    real(R_), intent(out) :: x, y   ! radius, x/y coordinate
    real(R_) :: r
    integer  :: i

    ! Rejection method
    do i = 1, 100
       x = 1.0_R_ - 2.0_R_ * mseq_rand_R()
       y = 1.0_R_ - 2.0_R_ * mseq_rand_R()
       rn = x**2 + y**2
       if (rn > r2min .and. rn <= 1.0_R_) exit
    end do

    ! The point
    if (i > 100) then ! unexpected case
       x = 1.0_R_
       y = 0.0_R_
       rn = 1.0_R_
    else              ! usual case
       r = 1.0_R_ / sqrt(rn)
       x = x * r
       y = y * r
    end if

  end subroutine rand_point_circle


  !+
  ! Determines a random point on the surface of the unit sphere
  !  This can be used to get a random unit vector (x,y,z) and 
  !   a uniform random number (=r*r2)
  !-
  subroutine rand_point_sphere(r2min, rn, x, y, z)

    real(R_),  intent(in)  :: r2min   ! (min radius)^2
    real(R_),  intent(out) :: rn      ! a random number
    real(R_),  intent(out) :: x, y, z ! x/y/z coordinates
    real(R_) :: r, r2
    integer  :: i

    !- Rejection method
    do i = 1, 200
       x = 1.0_R_ - 2.0_R_ * mseq_rand_R()
       y = 1.0_R_ - 2.0_R_ * mseq_rand_R()
       z = 1.0_R_ - 2.0_R_ * mseq_rand_R()
       r2 = x**2 + y**2 + z**2
       if (r2 > r2min .and. r2 <= 1.0_R_) exit
    end do

    !- The point
    if (i > 200) then ! unexpected case
       x = 0.0_R_
       y = 0.0_R_
       z = 1.0_R_
       rn = 1.0_R_
    else              ! usual case
       r = sqrt(r2)
       rn = r * r2
       r = 1.0_R_ / r
       x = x * r
       y = y * r
       z = z * r
    end if

  end subroutine rand_point_sphere


  !+
  ! Initialize the stack
  !-
  subroutine randStk_alloc(nstk) 

    integer,  intent(in) :: nstk
    Stk_nummax = nstk
    Stk_num = 0
    allocate (Stk_rand(nstk))

  end subroutine randStk_alloc


  !+
  ! Finalize the stack
  !-
  subroutine randStk_dealloc() 

    deallocate (Stk_rand)

  end subroutine randStk_dealloc


  !+
  ! Set automatically values to the stack
  !-
  subroutine randStk_setAuto(nval)

    integer,  intent(in) :: nval
    integer  :: ival, n

    if (nval <= 0) then
       n = Stk_nummax
    else
       n = min(Stk_nummax, nval)
    end if

    do ival = 1, n
       Stk_rand(ival) = mseq_rand_R()
    end do
    Stk_num = n

  end subroutine randStk_setAuto


  !+
  ! Fill in the stack with user-specified values
  !-
  subroutine randStk_setUser(nval, vals) 

    integer,  intent(in) :: nval
    real(R_), intent(in) :: vals(:)
    integer  :: n
    n = min(Stk_nummax, nval)
    Stk_rand(1:n) = vals(1:n)
    Stk_num = n

  end subroutine randStk_setUser


  !+
  ! A random number from the stack
  !-
  function randStk_rand_R() result(xi)

    real(R_) :: xi    ! a random number

    if (Stk_num >= 1) then  ! if the stack is not empty
       xi = Stk_rand(Stk_num)
       Stk_num = Stk_num - 1
    else                      ! if empty
       xi = mseq_rand_R()
    end if

  end function randStk_rand_R


  !+
  ! A random number from the stack
  !-
  function randStk_rand_1R(nval) result(xi)

    integer,  intent(in) :: nval
    real(R_) :: xi(nval)    ! random numbers
    integer  :: ival

    ! If enough values are available
    if (Stk_num >= nval) then
       xi(1:nval) = Stk_rand(Stk_num - nval + 1 : Stk_num)
       Stk_num = Stk_num - nval

       ! If the stack is not empty
    else if (Stk_num >= 1) then
       xi(1:Stk_num) = Stk_rand(1:Stk_num)
       do ival = Stk_num + 1, nval
          xi(ival) = mseq_rand_R()
       end do
       Stk_num = 0

       ! If empty
    else
       do ival = 1, nval
          xi(ival) = mseq_rand_R()
       end do
    end if

  end function randStk_rand_1R


  !+
  ! Russian roulette: survive or killed?
  !-
  subroutine russianRoulette(w, wact, wrr)

    real(R_), intent(inout) :: w     ! weight
    real(R_), intent(in)    :: wact  ! threshold for activation
    real(R_), intent(in)    :: wrr   ! cutoff threshold (should be >= wact)

    if (w < wact) then
       if (w > wrr * mseq_rand_R()) then
          w = wrr             ! survive
       else
          w = 0.0_R_             ! killed
       end if
    end if

  end subroutine russianRoulette


  !+
  ! Compute average & standard deviation from limited samples
  !-
  subroutine stat_approxStat(buf, ndat, nave_max, aver, sdev)

    real(R_), intent(in)  :: buf(:)    ! data buffer
    integer,  intent(in)  :: ndat      ! # of data
    integer,  intent(in)  :: nave_max  ! max # of samples for calculating statistics
    real(R_), intent(out) :: aver      ! average
    real(R_), intent(out) :: sdev      ! standard deviation
    real(R_) :: sum1, sum2
    integer  :: iave, idat, nave

    nave = min(nave_max, ndat)
    if (ndat <= nave_max) then                          ! regular sampling
       sum1 = sum(buf(1:ndat))
       sum2 = sum(buf(1:ndat)**2)
    else
       sum1 = 0.0_R_
       sum2 = 0.0_R_
       do iave = 1, nave
          idat = min(ndat, int(mseq_rand_R() * ndat) + 1) ! random sampling
          sum1 = sum1 + buf(idat)
          sum2 = sum2 + buf(idat)**2
       end do
    end if
    aver = sum1 / nave
    sdev = sqrt(max(0.0_R_, sum2 / nave - aver**2))

  end subroutine stat_approxStat

end module hparx_rand
