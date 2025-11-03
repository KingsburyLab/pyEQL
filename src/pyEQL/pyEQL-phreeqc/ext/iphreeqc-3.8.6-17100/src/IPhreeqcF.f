!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION AccumulateLine(ID,LINE)
        IMPLICIT NONE
        INTEGER(KIND=4)  :: ID
        CHARACTER(LEN=*) :: LINE
        INTEGER(KIND=4)  :: AccumulateLine
        INTEGER(KIND=4)  :: AccumulateLineF
        AccumulateLine = AccumulateLineF(ID,LINE)
      END FUNCTION AccumulateLine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION AddError(ID,ERROR_MSG)
        IMPLICIT NONE
        INTEGER(KIND=4)  :: ID
        CHARACTER(LEN=*) :: ERROR_MSG
        INTEGER(KIND=4)  :: AddError
        INTEGER(KIND=4)  :: AddErrorF
        AddError = AddErrorF(ID,ERROR_MSG)
      END FUNCTION AddError
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION AddWarning(ID,WARN_MSG)
        IMPLICIT NONE
        INTEGER(KIND=4)  :: ID
        CHARACTER(LEN=*) :: WARN_MSG
        INTEGER(KIND=4)  :: AddWarning
        INTEGER(KIND=4)  :: AddWarningF
        AddWarning = AddWarningF(ID,WARN_MSG)
      END FUNCTION AddWarning
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION ClearAccumulatedLines(ID)
        IMPLICIT NONE
        INTEGER(KIND=4)  :: ID
        INTEGER(KIND=4)  :: ClearAccumulatedLines
        INTEGER(KIND=4)  :: ClearAccumulatedLinesF
        ClearAccumulatedLines = ClearAccumulatedLinesF(ID)
      END FUNCTION ClearAccumulatedLines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION CreateIPhreeqc()
        IMPLICIT NONE
        INTEGER(KIND=4)  :: CreateIPhreeqc
        INTEGER(KIND=4)  :: CreateIPhreeqcF
        CreateIPhreeqc = CreateIPhreeqcF()
      END FUNCTION CreateIPhreeqc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION DestroyIPhreeqc(ID)
        IMPLICIT NONE
        INTEGER(KIND=4)  :: ID
        INTEGER(KIND=4)  :: DestroyIPhreeqc
        INTEGER(KIND=4)  :: DestroyIPhreeqcF
        DestroyIPhreeqc = DestroyIPhreeqcF(ID)
      END FUNCTION DestroyIPhreeqc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION GetComponent(ID,N,COMP)
        IMPLICIT NONE
        INTEGER(KIND=4)  :: ID
        INTEGER(KIND=4)  :: N
        CHARACTER(LEN=*) :: COMP
        INTEGER(KIND=4)  :: GetComponent
        INTEGER(KIND=4)  :: GetComponentF
        GetComponent = GetComponentF(ID,N,COMP)
      END FUNCTION GetComponent
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION GetComponentCount(ID)
        IMPLICIT NONE
        INTEGER(KIND=4) :: ID
        INTEGER(KIND=4) :: GetComponentCount
        INTEGER(KIND=4) :: GetComponentCountF
        GetComponentCount = GetComponentCountF(ID)
      END FUNCTION GetComponentCount
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION GetCurrentSelectedOutputUserNumber(ID)
        IMPLICIT NONE
        INTEGER(KIND=4) :: ID
        INTEGER(KIND=4) :: GetCurrentSelectedOutputUserNumber
        INTEGER(KIND=4) :: GetCurrentSelectedOutputUserNumberF
        GetCurrentSelectedOutputUserNumber = 
     &                     GetCurrentSelectedOutputUserNumberF(ID)
      END FUNCTION GetCurrentSelectedOutputUserNumber
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE GetDumpFileName(ID,FNAME)
        IMPLICIT NONE
        INTEGER(KIND=4)  :: ID
        CHARACTER(LEN=*) :: FNAME
        CALL GetDumpFileNameF(ID,FNAME)
      END SUBROUTINE GetDumpFileName
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION GetDumpFileOn(ID)
        IMPLICIT NONE
        INTEGER(KIND=4) :: ID
        LOGICAL(KIND=4) :: GetDumpFileOn
        INTEGER(KIND=4) :: GetDumpFileOnF
        IF (GetDumpFileOnF(ID).EQ.0) THEN
          GetDumpFileOn = .FALSE.
        ELSE
          GetDumpFileOn = .TRUE.
        ENDIF
      END FUNCTION GetDumpFileOn
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! GetDumpString      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION GetDumpStringLine(ID,N,LINE)
        IMPLICIT NONE
        INTEGER(KIND=4)  :: ID
        INTEGER(KIND=4)  :: N
        CHARACTER(LEN=*) :: LINE
        INTEGER(KIND=4)  :: GetDumpStringLine
        INTEGER(KIND=4)  :: GetDumpStringLineF
        GetDumpStringLine = GetDumpStringLineF(ID,N,LINE)
      END FUNCTION GetDumpStringLine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION GetDumpStringLineCount(ID)
        IMPLICIT NONE
        INTEGER(KIND=4) :: ID
        INTEGER(KIND=4) :: GetDumpStringLineCount
        INTEGER(KIND=4) :: GetDumpStringLineCountF
        GetDumpStringLineCount = GetDumpStringLineCountF(ID)
      END FUNCTION GetDumpStringLineCount
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION GetDumpStringOn(ID)
        IMPLICIT NONE
        INTEGER(KIND=4) :: ID
        LOGICAL(KIND=4) :: GetDumpStringOn
        INTEGER(KIND=4) :: GetDumpStringOnF
        IF (GetDumpStringOnF(ID).EQ.0) THEN
          GetDumpStringOn = .FALSE.
        ELSE
          GetDumpStringOn = .TRUE.
        ENDIF
      END FUNCTION GetDumpStringOn
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE GetErrorFileName(ID,FNAME)
        IMPLICIT NONE
        INTEGER(KIND=4)  :: ID
        CHARACTER(LEN=*) :: FNAME
        CALL GetErrorFileNameF(ID,FNAME)
      END SUBROUTINE GetErrorFileName
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION GetErrorFileOn(ID)
        IMPLICIT NONE
        INTEGER(KIND=4) :: ID
        LOGICAL(KIND=4) :: GetErrorFileOn
        INTEGER(KIND=4) :: GetErrorFileOnF
        IF (GetErrorFileOnF(ID).EQ.0) THEN
          GetErrorFileOn = .FALSE.
        ELSE
          GetErrorFileOn = .TRUE.
        ENDIF
      END FUNCTION GetErrorFileOn
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION GetErrorOn(ID)
        IMPLICIT NONE
        INTEGER(KIND=4) :: ID
        LOGICAL(KIND=4) :: GetErrorOn
        INTEGER(KIND=4) :: GetErrorOnF
        IF (GetErrorOnF(ID).EQ.0) THEN
          GetErrorOn = .FALSE.
        ELSE
          GetErrorOn = .TRUE.
        ENDIF
      END FUNCTION GetErrorOn
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! GetErrorString      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION GetErrorStringLine(ID,N,LINE)
        IMPLICIT NONE
        INTEGER(KIND=4)  :: ID
        INTEGER(KIND=4)  :: N
        CHARACTER(LEN=*) :: LINE
        INTEGER(KIND=4)  :: GetErrorStringLine
        INTEGER(KIND=4)  :: GetErrorStringLineF
        GetErrorStringLine = GetErrorStringLineF(ID,N,LINE)
      END FUNCTION GetErrorStringLine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION GetErrorStringLineCount(ID)
        IMPLICIT NONE
        INTEGER(KIND=4) :: ID
        INTEGER(KIND=4) :: GetErrorStringLineCount
        INTEGER(KIND=4) :: GetErrorStringLineCountF
        GetErrorStringLineCount = GetErrorStringLineCountF(ID)
      END FUNCTION GetErrorStringLineCount
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION GetErrorStringOn(ID)
        IMPLICIT NONE
        INTEGER(KIND=4) :: ID
        LOGICAL(KIND=4) :: GetErrorStringOn
        INTEGER(KIND=4) :: GetErrorStringOnF
        IF (GetErrorStringOnF(ID).EQ.0) THEN
          GetErrorStringOn = .FALSE.
        ELSE
          GetErrorStringOn = .TRUE.
        ENDIF
      END FUNCTION GetErrorStringOn
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE GetLogFileName(ID,FNAME)
        IMPLICIT NONE
        INTEGER(KIND=4)  :: ID
        CHARACTER(LEN=*) :: FNAME
        CALL GetLogFileNameF(ID,FNAME)
      END SUBROUTINE GetLogFileName
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION GetLogFileOn(ID)
        IMPLICIT NONE
        INTEGER(KIND=4) :: ID
        LOGICAL(KIND=4) :: GetLogFileOn
        INTEGER(KIND=4) :: GetLogFileOnF
        IF (GetLogFileOnF(ID).EQ.0) THEN
          GetLogFileOn = .FALSE.
        ELSE
          GetLogFileOn = .TRUE.
        ENDIF
      END FUNCTION GetLogFileOn
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! GetLogString      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION GetLogStringLine(ID,N,LINE)
        IMPLICIT NONE
        INTEGER(KIND=4)  :: ID
        INTEGER(KIND=4)  :: N
        CHARACTER(LEN=*) :: LINE
        INTEGER(KIND=4)  :: GetLogStringLine
        INTEGER(KIND=4)  :: GetLogStringLineF
        GetLogStringLine = GetLogStringLineF(ID,N,LINE)
      END FUNCTION GetLogStringLine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION GetLogStringLineCount(ID)
        IMPLICIT NONE
        INTEGER(KIND=4) :: ID
        INTEGER(KIND=4) :: GetLogStringLineCount
        INTEGER(KIND=4) :: GetLogStringLineCountF
        GetLogStringLineCount = GetLogStringLineCountF(ID)
      END FUNCTION GetLogStringLineCount
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION GetNthSelectedOutputUserNumber(ID,N)
        IMPLICIT NONE
        INTEGER(KIND=4) :: ID
        INTEGER(KIND=4) :: N
        INTEGER(KIND=4) :: GetNthSelectedOutputUserNumber
        INTEGER(KIND=4) :: GetNthSelectedOutputUserNumberF
        GetNthSelectedOutputUserNumber = 
     &                     GetNthSelectedOutputUserNumberF(ID,N)
      END FUNCTION GetNthSelectedOutputUserNumber
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION GetLogStringOn(ID)
        IMPLICIT NONE
        INTEGER(KIND=4) :: ID
        LOGICAL(KIND=4) :: GetLogStringOn
        INTEGER(KIND=4) :: GetLogStringOnF
        IF (GetLogStringOnF(ID).EQ.0) THEN
          GetLogStringOn = .FALSE.
        ELSE
          GetLogStringOn = .TRUE.
        ENDIF
      END FUNCTION GetLogStringOn
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE GetOutputFileName(ID,FNAME)
        IMPLICIT NONE
        INTEGER(KIND=4)  :: ID
        CHARACTER(LEN=*) :: FNAME
        CALL GetOutputFileNameF(ID,FNAME)
      END SUBROUTINE GetOutputFileName
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION GetOutputFileOn(ID)
        IMPLICIT NONE
        INTEGER(KIND=4) :: ID
        LOGICAL(KIND=4) :: GetOutputFileOn
        INTEGER(KIND=4) :: GetOutputFileOnF
        IF (GetOutputFileOnF(ID).EQ.0) THEN
          GetOutputFileOn = .FALSE.
        ELSE
          GetOutputFileOn = .TRUE.
        ENDIF
      END FUNCTION GetOutputFileOn
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! GetOutputString      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION GetOutputStringLine(ID,N,LINE)
        IMPLICIT NONE
        INTEGER(KIND=4)  :: ID
        INTEGER(KIND=4)  :: N
        CHARACTER(LEN=*) :: LINE
        INTEGER(KIND=4)  :: GetOutputStringLine
        INTEGER(KIND=4)  :: GetOutputStringLineF
        GetOutputStringLine = GetOutputStringLineF(ID,N,LINE)
      END FUNCTION GetOutputStringLine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION GetOutputStringLineCount(ID)
        IMPLICIT NONE
        INTEGER(KIND=4) :: ID
        INTEGER(KIND=4) :: GetOutputStringLineCount
        INTEGER(KIND=4) :: GetOutputStringLineCountF
        GetOutputStringLineCount = GetOutputStringLineCountF(ID)
      END FUNCTION GetOutputStringLineCount
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION GetOutputStringOn(ID)
        IMPLICIT NONE
        INTEGER(KIND=4) :: ID
        LOGICAL(KIND=4) :: GetOutputStringOn
        INTEGER(KIND=4) :: GetOutputStringOnF
        IF (GetOutputStringOnF(ID).EQ.0) THEN
          GetOutputStringOn = .FALSE.
        ELSE
          GetOutputStringOn = .TRUE.
        ENDIF
      END FUNCTION GetOutputStringOn
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION GetSelectedOutputColumnCount(ID)
        IMPLICIT NONE
        INTEGER(KIND=4) :: ID
        INTEGER(KIND=4) :: GetSelectedOutputColumnCount
        INTEGER(KIND=4) :: GetSelectedOutputColumnCountF
        GetSelectedOutputColumnCount = GetSelectedOutputColumnCountF(ID)
      END FUNCTION GetSelectedOutputColumnCount
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION GetSelectedOutputCount(ID)
        IMPLICIT NONE
        INTEGER(KIND=4) :: ID
        INTEGER(KIND=4) :: GetSelectedOutputCount
        INTEGER(KIND=4) :: GetSelectedOutputCountF
        GetSelectedOutputCount = GetSelectedOutputCountF(ID)
      END FUNCTION GetSelectedOutputCount
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE GetSelectedOutputFileName(ID,FNAME)
        IMPLICIT NONE
        INTEGER(KIND=4)  :: ID
        CHARACTER(LEN=*) :: FNAME
        CALL GetSelectedOutputFileNameF(ID,FNAME)
      END SUBROUTINE GetSelectedOutputFileName
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION GetSelectedOutputFileOn(ID)
        IMPLICIT NONE
        INTEGER(KIND=4) :: ID
        LOGICAL(KIND=4) :: GetSelectedOutputFileOn
        INTEGER(KIND=4) :: GetSelectedOutputFileOnF
        IF (GetSelectedOutputFileOnF(ID).EQ.0) THEN
          GetSelectedOutputFileOn = .FALSE.
        ELSE
          GetSelectedOutputFileOn = .TRUE.
        ENDIF
      END FUNCTION GetSelectedOutputFileOn
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION GetSelectedOutputRowCount(ID)
        IMPLICIT NONE
        INTEGER(KIND=4) :: ID
        INTEGER(KIND=4) :: GetSelectedOutputRowCount
        INTEGER(KIND=4) :: GetSelectedOutputRowCountF
        GetSelectedOutputRowCount = GetSelectedOutputRowCountF(ID)
      END FUNCTION GetSelectedOutputRowCount
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION GetSelectedOutputStringLine(ID,N,LINE)
        IMPLICIT NONE
        INTEGER(KIND=4)  :: ID
        INTEGER(KIND=4)  :: N
        CHARACTER(LEN=*) :: LINE        
        INTEGER(KIND=4)  :: GetSelectedOutputStringLine
        INTEGER(KIND=4)  :: GetSelectedOutputStringLineF
        GetSelectedOutputStringLine =
     &                      GetSelectedOutputStringLineF(ID,N,LINE)
      END FUNCTION GetSelectedOutputStringLine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION GetSelectedOutputStringLineCount(ID)
        IMPLICIT NONE
        INTEGER(KIND=4) :: ID
        INTEGER(KIND=4) :: GetSelectedOutputStringLineCount
        INTEGER(KIND=4) :: GetSelectedOutputStringLineCountF
        GetSelectedOutputStringLineCount = 
     &                     GetSelectedOutputStringLineCountF(ID)
      END FUNCTION GetSelectedOutputStringLineCount
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION GetSelectedOutputStringOn(ID)
        IMPLICIT NONE
        INTEGER(KIND=4) :: ID
        LOGICAL(KIND=4) :: GetSelectedOutputStringOn
        INTEGER(KIND=4) :: GetSelectedOutputStringOnF
        IF (GetSelectedOutputStringOnF(ID).EQ.0) THEN
          GetSelectedOutputStringOn = .FALSE.
        ELSE
          GetSelectedOutputStringOn = .TRUE.
        ENDIF
      END FUNCTION GetSelectedOutputStringOn
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION GetSelectedOutputValue(ID,ROW,COL,VTYPE,DVALUE,SVALUE)
        IMPLICIT NONE
        INTEGER(KIND=4)  :: ID
        INTEGER(KIND=4)  :: ROW
        INTEGER(KIND=4)  :: COL
        INTEGER(KIND=4)  :: VTYPE
        REAL(KIND=8)     :: DVALUE
        CHARACTER(LEN=*) :: SVALUE
        INTEGER(KIND=4)  :: GetSelectedOutputValue
        INTEGER(KIND=4)  :: GetSelectedOutputValueF
        GetSelectedOutputValue = GetSelectedOutputValueF(ID,ROW,
     &                     COL,VTYPE,DVALUE,SVALUE)
      END FUNCTION GetSelectedOutputValue
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE GetVersionString(VERSION)
        IMPLICIT NONE
        CHARACTER(LEN=*) :: VERSION
        CALL GetVersionStringF(VERSION)
      END SUBROUTINE GetVersionString
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! GetWarningString
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION GetWarningStringLine(ID,N,LINE)
        IMPLICIT NONE
        INTEGER(KIND=4)  :: ID
        INTEGER(KIND=4)  :: N
        CHARACTER(LEN=*) :: LINE        
        INTEGER(KIND=4)  :: GetWarningStringLine
        INTEGER(KIND=4)  :: GetWarningStringLineF
        GetWarningStringLine = GetWarningStringLineF(ID,N,LINE)
      END FUNCTION GetWarningStringLine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION GetWarningStringLineCount(ID)
        IMPLICIT NONE
        INTEGER(KIND=4) :: ID
        INTEGER(KIND=4) :: GetWarningStringLineCount
        INTEGER(KIND=4) :: GetWarningStringLineCountF
        GetWarningStringLineCount = GetWarningStringLineCountF(ID)
      END FUNCTION GetWarningStringLineCount
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION LoadDatabase(ID,FILENAME)
        IMPLICIT NONE
        INTEGER (KIND=4) :: ID
        CHARACTER(LEN=*) :: FILENAME
        INTEGER(KIND=4)  :: LoadDatabase
        INTEGER(KIND=4)  :: LoadDatabaseF
        LoadDatabase = LoadDatabaseF(ID,FILENAME)
      END FUNCTION LoadDatabase
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION LoadDatabaseString(ID,INPUT)
        IMPLICIT NONE
        INTEGER(KIND=4)  :: ID
        CHARACTER(LEN=*) :: INPUT
        INTEGER(KIND=4)  :: LoadDatabaseString
        INTEGER(KIND=4)  :: LoadDatabaseStringF
        LoadDatabaseString = LoadDatabaseStringF(ID,INPUT)
      END FUNCTION LoadDatabaseString
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE OutputAccumulatedLines(ID)
        IMPLICIT NONE
        INTEGER(KIND=4) :: ID
        CALL OutputAccumulatedLinesF(ID)
      END SUBROUTINE OutputAccumulatedLines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE OutputErrorString(ID)
        IMPLICIT NONE
        INTEGER(KIND=4) :: ID
        CALL OutputErrorStringF(ID)
      END SUBROUTINE OutputErrorString
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE OutputWarningString(ID)
        IMPLICIT NONE
        INTEGER(KIND=4) :: ID
        CALL OutputWarningStringF(ID)
      END SUBROUTINE OutputWarningString
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION RunAccumulated(ID)
        IMPLICIT NONE
        INTEGER(KIND=4)  :: ID
        INTEGER(KIND=4)  :: RunAccumulated
        INTEGER(KIND=4)  :: RunAccumulatedF
        RunAccumulated = RunAccumulatedF(ID)
      END FUNCTION RunAccumulated
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION RunFile(ID,FILENAME)
        IMPLICIT NONE
        INTEGER(KIND=4)  :: ID
        CHARACTER(LEN=*) :: FILENAME
        INTEGER(KIND=4)  :: RunFile
        INTEGER(KIND=4)  :: RunFileF
        RunFile = RunFileF(ID,FILENAME)
      END FUNCTION RunFile
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION RunString(ID,INPUT)
        IMPLICIT NONE
        INTEGER(KIND=4)  :: ID
        CHARACTER(LEN=*) :: INPUT
        INTEGER(KIND=4)  :: RunString
        INTEGER(KIND=4)  :: RunStringF
        RunString = RunStringF(ID,INPUT)
      END FUNCTION RunString
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION SetBasicFortranCallback(ID,COOKIE)
        IMPLICIT NONE
        INTEGER(KIND=4) :: ID
        INTERFACE 
          DOUBLE PRECISION FUNCTION cookie(x1, x2, str)
            DOUBLE PRECISION, INTENT(in) :: x1
            DOUBLE PRECISION, INTENT(in) :: x2
            CHARACTER(*), INTENT(in)     :: str
          END FUNCTION 
        END INTERFACE 
        INTEGER(KIND=4) :: SetBasicFortranCallback
        INTEGER(KIND=4) :: SetBasicFortranCallbackF
        SetBasicFortranCallback = SetBasicFortranCallbackF(ID,COOKIE)
      END FUNCTION SetBasicFortranCallback          
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION SetCurrentSelectedOutputUserNumber(ID,N)
        IMPLICIT NONE
        INTEGER(KIND=4)  :: ID
        INTEGER(KIND=4)  :: N
        INTEGER(KIND=4)  :: SetCurrentSelectedOutputUserNumber
        INTEGER(KIND=4)  :: SetCurrentSelectedOutputUserNumberF
        SetCurrentSelectedOutputUserNumber = 
     &                      SetCurrentSelectedOutputUserNumberF(ID,N)
      END FUNCTION SetCurrentSelectedOutputUserNumber
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION SetDumpFileName(ID,FNAME)
        IMPLICIT NONE
        INTEGER(KIND=4)  :: ID
        CHARACTER(LEN=*) :: FNAME
        INTEGER(KIND=4)  :: SetDumpFileName
        INTEGER(KIND=4)  :: SetDumpFileNameF
        SetDumpFileName = SetDumpFileNameF(ID,FNAME)
      END FUNCTION SetDumpFileName
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION SetDumpFileOn(ID,DUMP_ON)
        IMPLICIT NONE
        INTEGER(KIND=4) :: ID
        LOGICAL(KIND=4) :: DUMP_ON
        INTEGER(KIND=4) :: SetDumpFileOn
        INTEGER(KIND=4) :: SetDumpFileOnF
        SetDumpFileOn = SetDumpFileOnF(ID,DUMP_ON)
      END FUNCTION SetDumpFileOn
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION SetDumpStringOn(ID,DUMP_STRING_ON)
        IMPLICIT NONE
        INTEGER(KIND=4) :: ID
        LOGICAL(KIND=4) :: DUMP_STRING_ON
        INTEGER(KIND=4) :: SetDumpStringOn
        INTEGER(KIND=4) :: SetDumpStringOnF
        SetDumpStringOn = SetDumpStringOnF(ID,DUMP_STRING_ON)
      END FUNCTION SetDumpStringOn
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION SetErrorFileName(ID,FNAME)
        IMPLICIT NONE
        INTEGER(KIND=4)  :: ID
        CHARACTER(LEN=*) :: FNAME
        INTEGER(KIND=4)  :: SetErrorFileName
        INTEGER(KIND=4)  :: SetErrorFileNameF
        SetErrorFileName = SetErrorFileNameF(ID,FNAME)
      END FUNCTION SetErrorFileName
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION SetErrorFileOn(ID,ERROR_ON)
        IMPLICIT NONE
        INTEGER(KIND=4) :: ID
        LOGICAL(KIND=4) :: ERROR_ON
        INTEGER(KIND=4) :: SetErrorFileOn
        INTEGER(KIND=4) :: SetErrorFileOnF
        SetErrorFileOn = SetErrorFileOnF(ID,ERROR_ON)
      END FUNCTION SetErrorFileOn
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION SetErrorOn(ID,ERROR_ON)
        IMPLICIT NONE
        INTEGER(KIND=4) :: ID
        LOGICAL(KIND=4) :: ERROR_ON
        INTEGER(KIND=4) :: SetErrorOn
        INTEGER(KIND=4) :: SetErrorOnF
        SetErrorOn = SetErrorOnF(ID,ERROR_ON)
      END FUNCTION SetErrorOn
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION SetErrorStringOn(ID,ERROR_STRING_ON)
        IMPLICIT NONE
        INTEGER(KIND=4) :: ID
        LOGICAL(KIND=4) :: ERROR_STRING_ON
        INTEGER(KIND=4) :: SetErrorStringOn
        INTEGER(KIND=4) :: SetErrorStringOnF
        SetErrorStringOn = SetErrorStringOnF(ID,ERROR_STRING_ON)
      END FUNCTION SetErrorStringOn
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION SetLogFileName(ID,FNAME)
        IMPLICIT NONE
        INTEGER(KIND=4)  :: ID
        CHARACTER(LEN=*) :: FNAME
        INTEGER(KIND=4)  :: SetLogFileName
        INTEGER(KIND=4)  :: SetLogFileNameF
        SetLogFileName = SetLogFileNameF(ID,FNAME)
      END FUNCTION SetLogFileName
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION SetLogFileOn(ID,LOG_ON)
        IMPLICIT NONE
        INTEGER(KIND=4) :: ID
        LOGICAL(KIND=4) :: LOG_ON
        INTEGER(KIND=4) :: SetLogFileOn
        INTEGER(KIND=4) :: SetLogFileOnF
        SetLogFileOn = SetLogFileOnF(ID,LOG_ON)
      END FUNCTION SetLogFileOn
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION SetLogStringOn(ID,LOG_STRING_ON)
        IMPLICIT NONE
        INTEGER(KIND=4) :: ID
        LOGICAL(KIND=4) :: LOG_STRING_ON
        INTEGER(KIND=4) :: SetLogStringOn
        INTEGER(KIND=4) :: SetLogStringOnF
        SetLogStringOn = SetLogStringOnF(ID,LOG_STRING_ON)
      END FUNCTION SetLogStringOn
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION SetOutputFileName(ID,FNAME)
        IMPLICIT NONE
        INTEGER(KIND=4)  :: ID
        CHARACTER(LEN=*) :: FNAME
        INTEGER(KIND=4)  :: SetOutputFileName
        INTEGER(KIND=4)  :: SetOutputFileNameF
        SetOutputFileName = SetOutputFileNameF(ID,FNAME)
      END FUNCTION SetOutputFileName
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION SetOutputFileOn(ID,OUTPUT_FILE_ON)
        IMPLICIT NONE
        INTEGER(KIND=4) :: ID
        LOGICAL(KIND=4) :: OUTPUT_FILE_ON
        INTEGER(KIND=4) :: SetOutputFileOn
        INTEGER(KIND=4) :: SetOutputFileOnF
        SetOutputFileOn = SetOutputFileOnF(ID,OUTPUT_FILE_ON)
      END FUNCTION SetOutputFileOn
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION SetOutputStringOn(ID,OUTPUT_STRING_ON)
        IMPLICIT NONE
        INTEGER(KIND=4) :: ID
        LOGICAL(KIND=4) :: OUTPUT_STRING_ON
        INTEGER(KIND=4) :: SetOutputStringOn
        INTEGER(KIND=4) :: SetOutputStringOnF
        SetOutputStringOn = SetOutputStringOnF(ID,OUTPUT_STRING_ON)
      END FUNCTION SetOutputStringOn
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION SetSelectedOutputFileName(ID,FNAME)
        IMPLICIT NONE
        INTEGER(KIND=4)  :: ID
        CHARACTER(LEN=*) :: FNAME
        INTEGER(KIND=4)  :: SetSelectedOutputFileName
        INTEGER(KIND=4)  :: SetSelectedOutputFileNameF
        SetSelectedOutputFileName = SetSelectedOutputFileNameF(ID,FNAME)
      END FUNCTION SetSelectedOutputFileName
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION SetSelectedOutputFileOn(ID,SELOUT_FILE_ON)
        IMPLICIT NONE
        INTEGER(KIND=4) :: ID
        LOGICAL(KIND=4) :: SELOUT_FILE_ON
        INTEGER(KIND=4) :: SetSelectedOutputFileOn
        INTEGER(KIND=4) :: SetSelectedOutputFileOnF
        SetSelectedOutputFileOn = SetSelectedOutputFileOnF(ID,
     &                     SELOUT_FILE_ON)
      END FUNCTION SetSelectedOutputFileOn
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION SetSelectedOutputStringOn(ID,SELOUT_STRING_ON)
        IMPLICIT NONE
        INTEGER(KIND=4) :: ID
        LOGICAL(KIND=4) :: SELOUT_STRING_ON
        INTEGER(KIND=4) :: SetSelectedOutputStringOn
        INTEGER(KIND=4) :: SetSelectedOutputStringOnF
        SetSelectedOutputStringOn = SetSelectedOutputStringOnF(ID,
     &                     SELOUT_STRING_ON)
      END FUNCTION SetSelectedOutputStringOn
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

