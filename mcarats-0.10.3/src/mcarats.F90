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
! MCARaTS top interface
!
! -MCARaTS version 0.10.3
! -Date: Sep. 6, 2016
! -Purpose: To simulate 3-D radiative transfer in atmosphere-ocean-land system
! -Copyright (C) 2006-2012 Hironobu Iwabuchi
! -Contact: Hiro Iwabuchi (hiroiwa at m.tohoku.ac.jp (" at " = @))
!
! Note: Only four interfaces are provided at present. However, other interfaces could be easily
!  made because this module is very simple.
!-
module mcarats 

  use globals
#if UseMPI == 1
  use globals_MPI
#endif
  use hparx
  use mcarWld
  implicit none
  private

  ! Public
  public :: mcarats__init      !* Initialize this module *
  public :: mcarats__final     !* Finalize this module *
  public :: mcarats__exec      !* Execute the MCARaTS core *
  public :: mcarats__exec_args !* Get command-line arguments and execute the core *
  !// Usage : Execute the followings in a correct order.
  !   1)  call mcarats__init()
  !   2a) call mcarats__exec_args(marg)  ! or use the next
  !   2b) call mcarats__exec(isol, ptot, inpfile, outfile)
  !   3)  call mcarats__final()

  ! Private
  integer, save :: Mcar_iroot = 0 ! rank index for the root PE
  integer, save :: Mcar_irank = 0 ! rank index for my PE
  integer, save :: Mcar_nrank = 1 ! total # of PEs
  
contains

  !+
  ! Initialize this module
  !-
  subroutine mcarats__init() 

#if UseMPI == 1
    call mcarats__MPI_init() ! MPI
#else
    Mcar_iroot = 0
    Mcar_irank = 0
    Mcar_nrank = 1
#endif
    call mcarWld__set_env(Mcar_iroot, Mcar_irank, Mcar_nrank) ! set simulation environment

  end subroutine mcarats__init


  !+
  ! Finalize this module
  !-
  subroutine mcarats__final 

#if UseMPI == 1
    call mcarats__MPI_final()
#endif

  end subroutine mcarats__final


  !+
  ! Get command-line arguments and execute the core
  !-
  subroutine mcarats__exec_args(marg) 

    integer, intent(in), optional :: marg ! method for obtaining the arguments (default = 1)
    !// marg = 0 for reading from standard input
    !          1 for reading from the command line (NON-STANDARD Fortran)
    character(256) :: inpfile, outfile
    real(R_) :: ptot
    integer  :: isol, marg1

    ! Get arguments
    if (Mcar_irank == Mcar_iroot) then
       marg1 = 1 ! default
       if (present(marg)) marg1 = marg
       call mcarats__args(marg1, isol, ptot, inpfile, outfile) ! arguments
    end if
#if UseMPI == 1
    call mcarats__MPI_bcast_args(isol, ptot, inpfile, outfile) ! broadcast
#endif

    ! Execute the core
    call mcarats__exec(isol, ptot, inpfile, outfile)

  end subroutine mcarats__exec_args


  !+
  ! Execute the MCARaTS core
  !-
  subroutine mcarats__exec(isol, ptot, inpfile, outfile) 

    integer,  intent(in) :: isol ! solver flag (0=F3D, 1=P3D, 2=1D(ICA))
    real(R_), intent(in) :: ptot ! total # of used photons
    character(*), intent(in) :: inpfile ! filename for namelist input
    character(*), intent(in) :: outfile ! filename for result output
    character(256) :: filepath
    integer :: iui

    filepath = filePath_AN(inpfile, '/') ! assuming UNIX-compatible
    iui = freeUnit_I(10) ! unit index for the namelist input file
    call open_seq(iui, inpfile, 'old')

    call mcarWld__init(iui, filepath, outfile, isol) ! initialize all
    call mcarWld__jobs(iui, ptot) ! do jobs
    call mcarWld__final() ! finalize all

    !rewind (iui) !DEBUG, to check correct deallocation
    !call mcarWld__init(iui, filepath, outfile, isol) ! reinitialize
    !call mcarWld__jobs(iui, ptot) ! do it again
    !call mcarWld__final()

    close (iui)

  end subroutine mcarats__exec


  !+
  ! Get arguments from the standard input or the command line
  !-
  subroutine mcarats__args(marg, isol, ptot, inpfile, outfile)

    integer,  intent(in) :: marg ! method for obtaining the arguments
    !// marg = 0 for reading from standard input
    !          1 for reading from the command line (NON-STANDARD Fortran)
    integer,  intent(out) :: isol ! solver flag
    real(R_), intent(out) :: ptot ! total # of used photons
    character(*), intent(out) :: inpfile  ! file name for input
    character(*), intent(out) :: outfile ! file name for output
    character(256) :: argv(10), usage
    character(512) :: argmsg
    integer :: ios, narg, i

    ! Initialize
    narg = 4
    usage = 'usage : mcarats ptot isol inpfile outfile\n\n'
    argmsg = '  ptot    : number of source photons (e.g., 1e6)\n' &
         & //'  isol    : solver flag (0:F3D, 1:P3D, 2:ICA)\n' &
         & //'  inpfile : filename for namelist input\n' &
         & //'  outfile : filename for result output'

    ! Get arguments from the user
    if (marg == 0) then ! reading from standard input
       write (*,*) 'mcarats : Put the following arguments (one per line).'
       call writeStr(6, argmsg)
       do i = 1, narg
          read (*, '(a)', iostat=ios) argv(i)
          call err_issue(ios, argmsg)
       end do
    else ! reading from the command line (NON-STANDARD)
       argmsg = trim(usage)//argmsg
       call getCmdArgs(narg, argv, argmsg) ! using NON-STANDARD getarg() and iargc()
    end if

    ! Get values
    read (argv(1), *, iostat=ios) ptot
    call err_issue(ios, argmsg)
    read (argv(2), *, iostat=ios) isol
    call err_issue(ios, argmsg)
    inpfile = argv(3)
    outfile = argv(4)
    call check_iR('mcarats__args: ptot', ptot, 0.99_R_)

  end subroutine mcarats__args


  ! The followings are used only when MPI is active
#if UseMPI == 1

  !+
  ! Initialize environment for mcarats (mainly for MPI)
  !-
  subroutine mcarats__MPI_init()

    include 'inc_mpi.f90'
    integer :: ierr
    call MPI_Init(ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, Mcar_irank, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, Mcar_nrank, ierr)
    Mcar_iroot = 0

  end subroutine mcarats__MPI_init


  !+
  ! Finalize MPI environment for mcarats
  !-
  subroutine mcarats__MPI_final() 

    include 'inc_mpi.f90'
    integer :: ierr
    call MPI_Finalize(ierr)

  end subroutine mcarats__MPI_final


  !+
  ! MPI broadcast data taken in mcarats_args
  !-
  subroutine mcarats__MPI_bcast_args(isol, ptot, inpfile, outfile)

    include 'inc_mpi.f90'
    integer,  intent(inout) :: isol ! solver flag
    real(R_), intent(inout) :: ptot ! total # of used photons
    character(*), intent(inout) :: inpfile  ! file name for input
    character(*), intent(inout) :: outfile ! file name for output
    integer :: ierr, n
    call MPI_Bcast(isol, 1, MPI_INTEGER, Mcar_iroot, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(ptot, 1, MPI_R_,    Mcar_iroot, MPI_COMM_WORLD, ierr)
    n = len(inpfile)
    call MPI_Bcast(inpfile, n, MPI_CHARACTER, Mcar_iroot, MPI_COMM_WORLD, ierr)
    n = len(outfile)
    call MPI_Bcast(outfile, n, MPI_CHARACTER, Mcar_iroot, MPI_COMM_WORLD, ierr)

  end subroutine mcarats__MPI_bcast_args
#endif

end module mcarats
