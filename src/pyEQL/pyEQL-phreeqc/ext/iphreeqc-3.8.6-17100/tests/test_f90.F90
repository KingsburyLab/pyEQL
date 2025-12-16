FUNCTION F_MAIN()
#ifndef IPHREEQC_NO_FORTRAN_MODULE
  USE IPhreeqc
#endif
  IMPLICIT NONE
#ifdef IPHREEQC_NO_FORTRAN_MODULE
  INCLUDE 'IPhreeqc.f90.inc'
#endif
  INTEGER(KIND=4) id
  
  INTEGER(KIND=4)   r
  INTEGER(KIND=4)   c
  INTEGER(KIND=4)   t
  REAL(KIND=8)      d
  CHARACTER(LEN=80) s
  
  INTEGER(KIND=4) F_MAIN
  INTEGER(KIND=4) TestGetSet
  INTEGER(KIND=4) TestGetSetInitOn
  INTEGER(KIND=4) TestGetSetName
  
  INTEGER(KIND=4),PARAMETER :: EXIT_SUCCESS = 0
  INTEGER(KIND=4),PARAMETER :: EXIT_FAILURE = 1

  id = CreateIPhreeqc()
  IF (id.LT.0) THEN
     F_MAIN = EXIT_FAILURE
     RETURN
  END IF
  
  IF (AccumulateLine(id, "SOLUTION 1").NE.IPQ_OK) THEN
     CALL OutputErrorString(id)
     F_MAIN = EXIT_FAILURE
     RETURN
  END IF 
  
  IF (ClearAccumulatedLines(id).NE.IPQ_OK) THEN
     CALL OutputErrorString(id)
     F_MAIN = EXIT_FAILURE
     RETURN
  END IF
  
  ! Dump
  IF (TestGetSet(id,GetDumpFileOn,SetDumpFileOn).NE.0) THEN
     F_MAIN = EXIT_FAILURE
     RETURN
  END IF
  
  ! Dump string
  IF (TestGetSet(id,GetDumpStringOn,SetDumpStringOn).NE.0) THEN
     F_MAIN = EXIT_FAILURE
     RETURN
  END IF
  
  ! Dump filename
  IF (TestGetSetName(id,GetDumpFileName,SetDumpFileName).NE.0) THEN
     F_MAIN = EXIT_FAILURE
     RETURN
  END IF
  
  ! Error file
  IF (TestGetSet(id,GetErrorFileOn,SetErrorFileOn).NE.0) THEN
     F_MAIN = EXIT_FAILURE
     RETURN
  END IF

  ! Error
  IF (TestGetSetInitOn(id,GetErrorOn,SetErrorOn).NE.0) THEN
     F_MAIN = EXIT_FAILURE
     RETURN
  END IF

  ! Log
  IF (TestGetSet(id,GetLogFileOn,SetLogFileOn).NE.0) THEN
     F_MAIN = EXIT_FAILURE
     RETURN
  END IF
  
  ! Output
  IF (TestGetSet(id,GetOutputFileOn,SetOutputFileOn).NE.0) THEN
     F_MAIN = EXIT_FAILURE
     RETURN
  END IF
  
  ! Output filename
  IF (TestGetSetName(id,GetOutputFileName,SetOutputFileName).NE.0) THEN
     F_MAIN = EXIT_FAILURE
     RETURN
  END IF  
  
  ! Selected output
  IF (TestGetSet(id,GetSelectedOutputFileOn,SetSelectedOutputFileOn).NE.0) THEN
     F_MAIN = EXIT_FAILURE
     RETURN
  END IF
  
  IF (LoadDatabase(id, "phreeqc.dat").NE.0) THEN
     CALL OutputErrorString(id)
     F_MAIN = EXIT_FAILURE
     RETURN
  END IF
  
  IF (SetOutputStringOn(id, .TRUE.).NE.IPQ_OK) THEN
     CALL OutputErrorString(id)
     F_MAIN = EXIT_FAILURE
     RETURN
  END IF  
  
  IF (RunFile(id, "ex2").NE.0) THEN
     CALL OutputErrorString(id)
     F_MAIN = EXIT_FAILURE
     RETURN
  END IF
  
  DO r=0,GetSelectedOutputRowCount(id)
     DO c=1,GetSelectedOutputColumnCount(id)
        IF (GetSelectedOutputValue(id,r,c,t,d,s).NE.IPQ_OK) THEN
           CALL OutputErrorString(id)
           F_MAIN = EXIT_FAILURE
           RETURN
        END IF
     END DO
  END DO
  
  DO r=1,GetOutputStringLineCount(id)
     CALL GetOutputStringLine(id, r, s)
  END DO 
    
  IF (DestroyIPhreeqc(id).NE.0) THEN
     CALL OutputErrorString(id)
     F_MAIN = EXIT_FAILURE
     RETURN
  END IF
  
  F_MAIN = EXIT_SUCCESS
  RETURN
  
END FUNCTION F_MAIN


FUNCTION TestGetSet(id,getFunc,setFunc)
#ifndef IPHREEQC_NO_FORTRAN_MODULE
  USE IPhreeqc
#endif

  IMPLICIT NONE
#ifdef IPHREEQC_NO_FORTRAN_MODULE
  INCLUDE 'IPhreeqc.f90.inc'
#endif
  INTEGER(KIND=4) id
  INTEGER(KIND=4) TestGetSet
  INTERFACE
     FUNCTION getFunc(id)
       INTEGER(KIND=4), INTENT(in) :: id
       LOGICAL(KIND=4) getFunc
     END FUNCTION getFunc
  END INTERFACE
  INTERFACE
     FUNCTION setFunc(id,flag)
       INTEGER(KIND=4), INTENT(in) :: id
       LOGICAL(KIND=4), INTENT(in) :: flag
       INTEGER(KIND=4) setFunc
     END FUNCTION setFunc
  END INTERFACE
  INTEGER(KIND=4),PARAMETER :: EXIT_SUCCESS = 0
  INTEGER(KIND=4),PARAMETER :: EXIT_FAILURE = 1

  IF (getFunc(id)) THEN
     TestGetSet = EXIT_FAILURE
     WRITE(*,*) "FAILURE"
     RETURN
  END IF

  IF (setFunc(id,.TRUE.).NE.IPQ_OK) THEN
     TestGetSet = EXIT_FAILURE
     WRITE(*,*) "FAILURE"
     RETURN
  END IF

  IF (.NOT.getFunc(id)) THEN
     TestGetSet = EXIT_FAILURE
     WRITE(*,*) "FAILURE"
     RETURN
  END IF

  IF (setFunc(id,.FALSE.).NE.IPQ_OK) THEN
     TestGetSet = EXIT_FAILURE
     WRITE(*,*) "FAILURE"
     RETURN
  END IF

  IF (getFunc(id)) THEN
     TestGetSet = EXIT_FAILURE
     WRITE(*,*) "FAILURE"
     RETURN
  END IF

  TestGetSet = EXIT_SUCCESS
  RETURN

END FUNCTION TestGetSet


FUNCTION TestGetSetInitOn(id,getFunc,setFunc)
#ifndef IPHREEQC_NO_FORTRAN_MODULE
  USE IPhreeqc
#endif

  IMPLICIT NONE
#ifdef IPHREEQC_NO_FORTRAN_MODULE
  INCLUDE 'IPhreeqc.f90.inc'
#endif
  INTEGER(KIND=4) id
  INTEGER(KIND=4) TestGetSetInitOn
  INTERFACE
     FUNCTION getFunc(id)
       INTEGER(KIND=4), INTENT(in) :: id
       LOGICAL(KIND=4) getFunc
     END FUNCTION getFunc
  END INTERFACE
  INTERFACE
     FUNCTION setFunc(id,flag)
       INTEGER(KIND=4), INTENT(in) :: id
       LOGICAL(KIND=4), INTENT(in) :: flag
       INTEGER(KIND=4) setFunc
     END FUNCTION setFunc
  END INTERFACE
  INTEGER(KIND=4),PARAMETER :: EXIT_SUCCESS = 0
  INTEGER(KIND=4),PARAMETER :: EXIT_FAILURE = 1

  IF (.NOT.(getFunc(id))) THEN
     TestGetSetInitOn = EXIT_FAILURE
     WRITE(*,*) "FAILURE" 
     RETURN
  END IF

  IF (setFunc(id,.TRUE.).NE.IPQ_OK) THEN
     TestGetSetInitOn = EXIT_FAILURE
     WRITE(*,*) "FAILURE" 
     RETURN
  END IF

  IF (.NOT.getFunc(id)) THEN
     TestGetSetInitOn = EXIT_FAILURE
     WRITE(*,*) "FAILURE" 
     RETURN
  END IF

  IF (setFunc(id,.FALSE.).NE.IPQ_OK) THEN
     TestGetSetInitOn = EXIT_FAILURE
     WRITE(*,*) "FAILURE" 
     RETURN
  END IF

  IF (getFunc(id)) THEN
     TestGetSetInitOn = EXIT_FAILURE
     WRITE(*,*) "FAILURE" 
     RETURN
  END IF

  TestGetSetInitOn = EXIT_SUCCESS
  RETURN

END FUNCTION TestGetSetInitOn


FUNCTION TestGetSetName(id,getFuncName,setFuncName)
#ifndef IPHREEQC_NO_FORTRAN_MODULE
  USE IPhreeqc
#endif
  
  IMPLICIT NONE
#ifdef IPHREEQC_NO_FORTRAN_MODULE
  INCLUDE 'IPhreeqc.f90.inc'
#endif
  INTEGER(KIND=4) id
  INTEGER(KIND=4) TestGetSetName
  INTERFACE
     SUBROUTINE getFuncName(id,fname)
       INTEGER(KIND=4),  INTENT(in)  :: id
       CHARACTER(LEN=*), INTENT(out) :: fname
     END SUBROUTINE getFuncName
  END INTERFACE
  INTERFACE
     FUNCTION setFuncName(id,fname)
       INTEGER(KIND=4),  INTENT(in) :: id
       CHARACTER(LEN=*), INTENT(in) :: fname
       INTEGER(KIND=4) setFuncName
     END FUNCTION setFuncName
  END INTERFACE
  INTEGER(KIND=4),PARAMETER :: EXIT_SUCCESS = 0
  INTEGER(KIND=4),PARAMETER :: EXIT_FAILURE = 1
  CHARACTER(LEN=80) FILEN
  
  CALL getFuncName(id,FILEN)
  
  FILEN = 'ABCDEFG'
  
  IF (setFuncName(id,FILEN).NE.IPQ_OK) THEN
     TestGetSetName = EXIT_FAILURE
     WRITE(*,*) "FAILURE" 
     RETURN
  END IF
  
  CALL getFuncName(id,FILEN)
  IF (.NOT.LLE('ABCDEFG', FILEN)) THEN
     TestGetSetName = EXIT_FAILURE
     WRITE(*,*) "FAILURE" 
     RETURN
  END IF
  
  IF (setFuncName(id,'XYZ').NE.IPQ_OK) THEN
     TestGetSetName = EXIT_FAILURE
     WRITE(*,*) "FAILURE" 
     RETURN
  END IF
  
  CALL getFuncName(id,FILEN)
  IF (.NOT.LLE('XYZ', FILEN)) THEN
     TestGetSetName = EXIT_FAILURE
     WRITE(*,*) "FAILURE" 
     RETURN
  END IF  
  
  TestGetSetName = EXIT_SUCCESS
  RETURN
  
END FUNCTION TestGetSetName
