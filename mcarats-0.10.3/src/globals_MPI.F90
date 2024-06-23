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
! Global constants for MPI
!-
module globals_MPI 

#if UseMPI == 1

  include 'inc_mpi.f90'
  implicit none
  private

  ! Kinds of real or complex variables
  !integer,  parameter :: MPI_R_  = MPI_REAL
  integer,  parameter :: MPI_R_  = MPI_DOUBLE_PRECISION
  integer,  parameter :: MPI_RD_ = MPI_DOUBLE_PRECISION
  !// This module is provided for power users who know well MPI and how to encapsulate 
  !    types of Fortran floating point variables.
  !   The above two parameters should be consistent with types defined in hparx/globals.f90.
  !   MPI procedure call statements should all use encapsulated types MPI_R_ or MPI_RD_.
  !   Never use MPI_REAL or MPI_DOUBLE_PRECISION directly!

  ! Public
  public :: MPI_R_, MPI_RD_
  public :: globals_MPI__print !* Print values of some global constants *

contains

  !+
  ! Print values of some global constants
  !-
  subroutine globals_MPI__print() 

    write (*,*) 'MPI_R_, MPI_RD_ : ', MPI_R_, MPI_RD_

  end subroutine globals_MPI__print
  
#endif

end module globals_MPI


!program main
!
!  use globals_MPI
!
!#if UseMPI == 1
!  call globals_MPI__print()
!#endif
!  print *, 'Done.'
!
!end program main
