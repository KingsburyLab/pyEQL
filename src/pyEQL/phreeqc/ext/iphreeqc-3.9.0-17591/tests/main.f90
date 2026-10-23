!  main.f90 
!
!  FUNCTIONS:
!  main - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: test_f90
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************

    program test_f90

    implicit none

    INTEGER(KIND=4),PARAMETER :: EXIT_SUCCESS = 0
    INTEGER(KIND=4),PARAMETER :: EXIT_FAILURE = 1

    integer(KIND=4) F_MAIN
    integer(KIND=4) I

    I = F_MAIN()

    if (I .NE. EXIT_SUCCESS) then
      STOP EXIT_FAILURE
    endif

    end program test_f90
