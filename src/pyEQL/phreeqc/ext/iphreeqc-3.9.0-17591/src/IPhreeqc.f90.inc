       ! GetSelectedOutputValue TYPES
       INTEGER(KIND=4),PARAMETER :: TT_EMPTY  = 0
       INTEGER(KIND=4),PARAMETER :: TT_ERROR  = 1
       INTEGER(KIND=4),PARAMETER :: TT_DOUBLE = 3
       INTEGER(KIND=4),PARAMETER :: TT_STRING = 4


       ! ERROR RETURN TYPES
       INTEGER(KIND=4),PARAMETER :: IPQ_OK           =  0
       INTEGER(KIND=4),PARAMETER :: IPQ_OUTOFMEMORY  = -1
       INTEGER(KIND=4),PARAMETER :: IPQ_BADVARTYPE   = -2
       INTEGER(KIND=4),PARAMETER :: IPQ_INVALIDARG   = -3
       INTEGER(KIND=4),PARAMETER :: IPQ_INVALIDROW   = -4
       INTEGER(KIND=4),PARAMETER :: IPQ_INVALIDCOL   = -5
       INTEGER(KIND=4),PARAMETER :: IPQ_BADINSTANCE  = -6


       INTERFACE
        FUNCTION AccumulateLine(ID,LINE)
         INTEGER(KIND=4),  INTENT(IN) :: ID
         CHARACTER(LEN=*), INTENT(IN) :: LINE
         INTEGER(KIND=4)              :: AccumulateLine
        END FUNCTION AccumulateLine
       END INTERFACE


       INTERFACE
        FUNCTION AddError(ID,ERROR_MSG)
         INTEGER(KIND=4),  INTENT(IN) :: ID
         CHARACTER(LEN=*), INTENT(IN) :: ERROR_MSG
         INTEGER(KIND=4)              :: AddError
        END FUNCTION AddError
       END INTERFACE


       INTERFACE
        FUNCTION AddWarning(ID,WARN_MSG)
         INTEGER(KIND=4),  INTENT(IN) :: ID
         CHARACTER(LEN=*), INTENT(IN) :: WARN_MSG
         INTEGER(KIND=4)              :: AddWarning
        END FUNCTION AddWarning
       END INTERFACE
       
       
       INTERFACE
        FUNCTION ClearAccumulatedLines(ID)
         INTEGER(KIND=4), INTENT(IN) :: ID
         INTEGER(KIND=4)             :: ClearAccumulatedLines
        END FUNCTION ClearAccumulatedLines
       END INTERFACE


       INTERFACE
        FUNCTION CreateIPhreeqc()
         INTEGER(KIND=4) :: CreateIPhreeqc
        END FUNCTION CreateIPhreeqc
       END INTERFACE


       INTERFACE
        FUNCTION DestroyIPhreeqc(ID)
         INTEGER(KIND=4), INTENT(IN) :: ID
         INTEGER(KIND=4)             :: DestroyIPhreeqc
        END FUNCTION DestroyIPhreeqc
       END INTERFACE


       INTERFACE
        SUBROUTINE GetComponent(ID,N,COMP)
         INTEGER(KIND=4),  INTENT(IN)  :: ID
         INTEGER(KIND=4),  INTENT(IN)  :: N
         CHARACTER(LEN=*), INTENT(OUT) :: COMP
        END SUBROUTINE GetComponent
       END INTERFACE


       INTERFACE
        FUNCTION GetComponentCount(ID)
         INTEGER(KIND=4),  INTENT(IN)  :: ID
         INTEGER(KIND=4)               :: GetComponentCount
        END FUNCTION GetComponentCount
       END INTERFACE
       

       INTERFACE
        FUNCTION GetCurrentSelectedOutputUserNumber(ID)
         INTEGER(KIND=4),  INTENT(IN)  :: ID
         INTEGER(KIND=4)               :: GetCurrentSelectedOutputUserNumber
        END FUNCTION GetCurrentSelectedOutputUserNumber
       END INTERFACE


       INTERFACE
        SUBROUTINE GetDumpFileName(ID,FNAME)
         INTEGER(KIND=4),  INTENT(IN)  :: ID
         CHARACTER(LEN=*), INTENT(OUT) :: FNAME
        END SUBROUTINE GetDumpFileName
       END INTERFACE


       INTERFACE
        FUNCTION GetDumpFileOn(ID)
         INTEGER(KIND=4),  INTENT(IN)  :: ID
         LOGICAL(KIND=4)               :: GetDumpFileOn
        END FUNCTION GetDumpFileOn
       END INTERFACE


       INTERFACE
        SUBROUTINE GetDumpStringLine(ID,N,LINE)
         INTEGER(KIND=4),   INTENT(IN)   :: ID
         INTEGER(KIND=4),   INTENT(IN)   :: N
         CHARACTER(LEN=*),  INTENT(OUT)  :: LINE
        END SUBROUTINE
       END INTERFACE


       INTERFACE
        FUNCTION GetDumpStringLineCount(ID)
         INTEGER(KIND=4),  INTENT(IN)  :: ID
         INTEGER(KIND=4)               :: GetDumpStringLineCount
        END FUNCTION GetDumpStringLineCount
       END INTERFACE


       INTERFACE
        FUNCTION GetDumpStringOn(ID)
         INTEGER(KIND=4),  INTENT(IN)  :: ID
         LOGICAL(KIND=4)               :: GetDumpStringOn
        END FUNCTION GetDumpStringOn
       END INTERFACE


       INTERFACE
        SUBROUTINE GetErrorFileName(ID,FNAME)
         INTEGER(KIND=4),  INTENT(IN)  :: ID
         CHARACTER(LEN=*), INTENT(OUT) :: FNAME
        END SUBROUTINE GetErrorFileName
       END INTERFACE


       INTERFACE
        FUNCTION GetErrorFileOn(ID)
         INTEGER(KIND=4),  INTENT(IN)  :: ID
         LOGICAL(KIND=4)               :: GetErrorFileOn
        END FUNCTION GetErrorFileOn
       END INTERFACE


       INTERFACE
        FUNCTION GetErrorOn(ID)
         INTEGER(KIND=4),  INTENT(IN)  :: ID
         LOGICAL(KIND=4)               :: GetErrorOn
        END FUNCTION GetErrorOn
       END INTERFACE


       INTERFACE
        SUBROUTINE GetErrorStringLine(ID,N,LINE)
         INTEGER(KIND=4),  INTENT(IN)  :: ID
         INTEGER(KIND=4),  INTENT(IN)  :: N
         CHARACTER(LEN=*), INTENT(OUT) :: LINE
        END SUBROUTINE
       END INTERFACE

       INTERFACE
        FUNCTION GetErrorStringLineCount(ID)
         INTEGER(KIND=4),  INTENT(IN) :: ID
         INTEGER(KIND=4)              :: GetErrorStringLineCount
        END FUNCTION GetErrorStringLineCount
       END INTERFACE


       INTERFACE
        FUNCTION GetErrorStringOn(ID)
         INTEGER(KIND=4),  INTENT(IN)  :: ID
         LOGICAL(KIND=4)               :: GetErrorStringOn
        END FUNCTION GetErrorStringOn
       END INTERFACE


       INTERFACE
        SUBROUTINE GetLogFileName(ID,FNAME)
         INTEGER(KIND=4),  INTENT(IN)  :: ID
         CHARACTER(LEN=*), INTENT(OUT) :: FNAME
        END SUBROUTINE GetLogFileName
       END INTERFACE


       INTERFACE
        FUNCTION GetLogFileOn(ID)
         INTEGER(KIND=4),  INTENT(IN)  :: ID
         LOGICAL(KIND=4)               :: GetLogFileOn
        END FUNCTION GetLogFileOn
       END INTERFACE


       INTERFACE
        SUBROUTINE GetLogStringLine(ID,N,LINE)
         INTEGER(KIND=4),   INTENT(IN)   :: ID
         INTEGER(KIND=4),   INTENT(IN)   :: N
         CHARACTER(LEN=*),  INTENT(OUT)  :: LINE
        END SUBROUTINE
       END INTERFACE


       INTERFACE
        FUNCTION GetLogStringLineCount(ID)
         INTEGER(KIND=4),  INTENT(IN)  :: ID
         INTEGER(KIND=4)               :: GetLogStringLineCount
        END FUNCTION GetLogStringLineCount
       END INTERFACE


       INTERFACE
        FUNCTION GetLogStringOn(ID)
         INTEGER(KIND=4),  INTENT(IN) :: ID
         LOGICAL(KIND=4)              :: GetLogStringOn
        END FUNCTION GetLogStringOn
       END INTERFACE


       INTERFACE
        FUNCTION GetNthSelectedOutputUserNumber(ID,N)
         INTEGER(KIND=4),  INTENT(IN)   :: ID
         INTEGER(KIND=4),  INTENT(IN)   :: N
         INTEGER(KIND=4)                :: GetNthSelectedOutputUserNumber
        END FUNCTION GetNthSelectedOutputUserNumber
       END INTERFACE


       INTERFACE
        SUBROUTINE GetOutputFileName(ID,FNAME)
         INTEGER(KIND=4),  INTENT(IN)  :: ID
         CHARACTER(LEN=*), INTENT(OUT) :: FNAME
        END SUBROUTINE GetOutputFileName
       END INTERFACE


       INTERFACE
        FUNCTION GetOutputFileOn(ID)
         INTEGER(KIND=4),  INTENT(IN)  :: ID
         LOGICAL(KIND=4)               :: GetOutputFileOn
        END FUNCTION GetOutputFileOn
       END INTERFACE


       INTERFACE
        SUBROUTINE GetOutputStringLine(ID,N,LINE)
         INTEGER(KIND=4),   INTENT(IN)   :: ID
         INTEGER(KIND=4),   INTENT(IN)   :: N
         CHARACTER(LEN=*),  INTENT(OUT)  :: LINE
        END SUBROUTINE
       END INTERFACE


       INTERFACE
        FUNCTION GetOutputStringLineCount(ID)
         INTEGER(KIND=4),  INTENT(IN)  :: ID
         INTEGER(KIND=4)               :: GetOutputStringLineCount
        END FUNCTION GetOutputStringLineCount
       END INTERFACE
       
       
       INTERFACE
        FUNCTION GetOutputStringOn(ID)
         INTEGER(KIND=4),  INTENT(IN) :: ID
         LOGICAL(KIND=4)              :: GetOutputStringOn
        END FUNCTION GetOutputStringOn
       END INTERFACE


       INTERFACE
        FUNCTION GetSelectedOutputColumnCount(ID)
         INTEGER(KIND=4),  INTENT(IN) :: ID
         INTEGER(KIND=4)              :: GetSelectedOutputColumnCount
        END FUNCTION GetSelectedOutputColumnCount
       END INTERFACE


       INTERFACE
        FUNCTION GetSelectedOutputCount(ID)
         INTEGER(KIND=4),  INTENT(IN) :: ID
         INTEGER(KIND=4)              :: GetSelectedOutputCount
        END FUNCTION GetSelectedOutputCount
       END INTERFACE


       INTERFACE
        SUBROUTINE GetSelectedOutputFileName(ID,FNAME)
         INTEGER(KIND=4),  INTENT(IN)  :: ID
         CHARACTER(LEN=*), INTENT(OUT) :: FNAME
        END SUBROUTINE GetSelectedOutputFileName
       END INTERFACE


       INTERFACE
        FUNCTION GetSelectedOutputFileOn(ID)
         INTEGER(KIND=4),  INTENT(IN)  :: ID
         LOGICAL(KIND=4)               :: GetSelectedOutputFileOn
        END FUNCTION GetSelectedOutputFileOn
       END INTERFACE


       INTERFACE
        FUNCTION GetSelectedOutputRowCount(ID)
         INTEGER(KIND=4),  INTENT(IN) :: ID
         INTEGER(KIND=4)              :: GetSelectedOutputRowCount
        END FUNCTION GetSelectedOutputRowCount
       END INTERFACE


       INTERFACE
        SUBROUTINE GetSelectedOutputStringLine(ID,N,LINE)
         INTEGER(KIND=4),   INTENT(IN)   :: ID
         INTEGER(KIND=4),   INTENT(IN)   :: N
         CHARACTER(LEN=*),  INTENT(OUT)  :: LINE
        END SUBROUTINE
       END INTERFACE

       INTERFACE
        FUNCTION GetSelectedOutputStringLineCount(ID)
         INTEGER(KIND=4),  INTENT(IN)  :: ID
         INTEGER(KIND=4)  :: GetSelectedOutputStringLineCount
        END FUNCTION GetSelectedOutputStringLineCount
       END INTERFACE


       INTERFACE
        FUNCTION GetSelectedOutputStringOn(ID)
         INTEGER(KIND=4),  INTENT(IN) :: ID
         LOGICAL(KIND=4)              :: GetSelectedOutputStringOn
        END FUNCTION GetSelectedOutputStringOn
       END INTERFACE


       INTERFACE
        FUNCTION GetSelectedOutputValue(ID,ROW,COL,VTYPE,DVALUE,SVALUE)
         INTEGER(KIND=4),  INTENT(IN)  :: ID
         INTEGER(KIND=4),  INTENT(IN)  :: ROW
         INTEGER(KIND=4),  INTENT(IN)  :: COL
         INTEGER(KIND=4),  INTENT(OUT) :: VTYPE
         REAL(KIND=8),     INTENT(OUT) :: DVALUE
         CHARACTER(LEN=*), INTENT(OUT) :: SVALUE
         INTEGER(KIND=4)               :: GetSelectedOutputValue
        END FUNCTION GetSelectedOutputValue
       END INTERFACE


       INTERFACE
        SUBROUTINE GetVersionString(VER)
         CHARACTER(LEN=*), INTENT(OUT) :: VER
        END SUBROUTINE
       END INTERFACE


       INTERFACE
        SUBROUTINE GetWarningStringLine(ID,N,LINE)
         INTEGER(KIND=4),  INTENT(IN)  :: ID
         INTEGER(KIND=4),  INTENT(IN)  :: N
         CHARACTER(LEN=*), INTENT(OUT) :: LINE
        END SUBROUTINE
       END INTERFACE


       INTERFACE
        FUNCTION GetWarningStringLineCount(ID)
         INTEGER(KIND=4),  INTENT(IN) :: ID
         INTEGER(KIND=4)              :: GetWarningStringLineCount
        END FUNCTION GetWarningStringLineCount
       END INTERFACE


       INTERFACE
        FUNCTION LoadDatabase(ID,FILENAME)
         INTEGER(KIND=4),  INTENT(IN) :: ID
         CHARACTER(LEN=*), INTENT(IN) :: FILENAME
         INTEGER(KIND=4)              :: LoadDatabase
        END FUNCTION LoadDatabase
       END INTERFACE


       INTERFACE
        FUNCTION LoadDatabaseString(ID,INPUT)
         INTEGER(KIND=4),   INTENT(IN) :: ID
         CHARACTER(LEN=*),  INTENT(IN) :: INPUT
         INTEGER(KIND=4)               :: LoadDatabaseString
        END FUNCTION LoadDatabaseString
       END INTERFACE


       INTERFACE
        SUBROUTINE OutputAccumulatedLines(ID)
         INTEGER(KIND=4),  INTENT(IN) :: ID
        END SUBROUTINE OutputAccumulatedLines
       END INTERFACE


       INTERFACE 
        SUBROUTINE OutputErrorString(ID)
         INTEGER(KIND=4),  INTENT(IN) :: ID
        END SUBROUTINE OutputErrorString
       END INTERFACE


       INTERFACE 
        SUBROUTINE OutputWarningString(ID)
         INTEGER(KIND=4),  INTENT(IN) :: ID
        END SUBROUTINE OutputWarningString
       END INTERFACE


       INTERFACE
        FUNCTION RunAccumulated(ID)
         INTEGER(KIND=4),  INTENT(IN) :: ID
         INTEGER(KIND=4)              :: RunAccumulated
        END FUNCTION RunAccumulated
       END INTERFACE


       INTERFACE
        FUNCTION RunFile(ID,FNAME)
         INTEGER(KIND=4),  INTENT(IN) :: ID
         CHARACTER(LEN=*), INTENT(IN) :: FNAME
         INTEGER(KIND=4)              :: RunFile
        END FUNCTION RunFile
       END INTERFACE


       INTERFACE
        FUNCTION RunString(ID,INPUT)
         INTEGER(KIND=4),   INTENT(IN) :: ID
         CHARACTER(LEN=*),  INTENT(IN) :: INPUT
         INTEGER(KIND=4)               :: RunString
        END FUNCTION RunString
       END INTERFACE


       INTERFACE
        FUNCTION SetBasicFortranCallback(ID,COOKIE)
         INTEGER(KIND=4),  INTENT(IN) :: ID
        INTERFACE
          DOUBLE PRECISION FUNCTION cookie(x1, x2, str)
            DOUBLE PRECISION, INTENT(in) :: x1
            DOUBLE PRECISION, INTENT(in) :: x2
            CHARACTER(*), INTENT(in)                   :: str
          END FUNCTION 
        END INTERFACE
         INTEGER(KIND=4)              :: SetBasicFortranCallback
        END FUNCTION SetBasicFortranCallback
       END INTERFACE

       INTERFACE
        FUNCTION SetCurrentSelectedOutputUserNumber(ID,N)
         INTEGER(KIND=4), INTENT(IN)  :: ID
         INTEGER(KIND=4), INTENT(IN)  :: N
         INTEGER(KIND=4)              :: SetCurrentSelectedOutputUserNumber
        END FUNCTION SetCurrentSelectedOutputUserNumber
       END INTERFACE


       INTERFACE
        FUNCTION SetDumpFileName(ID,FNAME)
         INTEGER(KIND=4),  INTENT(IN) :: ID
         CHARACTER(LEN=*), INTENT(IN) :: FNAME
         INTEGER(KIND=4)              :: SetDumpFileName
        END FUNCTION SetDumpFileName
       END INTERFACE


       INTERFACE
        FUNCTION SetDumpFileOn(ID,DUMP_ON)
         INTEGER(KIND=4),  INTENT(IN) :: ID
         LOGICAL(KIND=4),  INTENT(IN) :: DUMP_ON
         INTEGER(KIND=4)              :: SetDumpFileOn
        END FUNCTION SetDumpFileOn
       END INTERFACE


       INTERFACE
        FUNCTION SetDumpStringOn(ID,DUMP_STRING_ON)
         INTEGER(KIND=4),  INTENT(IN) :: ID
         LOGICAL(KIND=4),  INTENT(IN) :: DUMP_STRING_ON
         INTEGER(KIND=4)              :: SetDumpStringOn
        END FUNCTION SetDumpStringOn
       END INTERFACE


       INTERFACE
        FUNCTION SetErrorFileNameF(ID,FNAME)
         INTEGER(KIND=4),  INTENT(IN) :: ID
         CHARACTER(LEN=*), INTENT(IN) :: FNAME
         INTEGER(KIND=4)              :: SetErrorFileNameF
        END FUNCTION SetErrorFileNameF
       END INTERFACE


       INTERFACE
        FUNCTION SetErrorFileOn(ID,ERR_FILE_ON)
         INTEGER(KIND=4),  INTENT(IN) :: ID
         LOGICAL(KIND=4),  INTENT(IN) :: ERR_FILE_ON
         INTEGER(KIND=4)              :: SetErrorFileOn
        END FUNCTION SetErrorFileOn
       END INTERFACE


       INTERFACE
        FUNCTION SetErrorOn(ID,ERR_ON)
         INTEGER(KIND=4),  INTENT(IN) :: ID
         LOGICAL(KIND=4),  INTENT(IN) :: ERR_ON
         INTEGER(KIND=4)              :: SetErrorOn
        END FUNCTION SetErrorOn
       END INTERFACE


       INTERFACE
        FUNCTION SetErrorStringOn(ID,ERR_STRING_ON)
         INTEGER(KIND=4),  INTENT(IN) :: ID
         LOGICAL(KIND=4),  INTENT(IN) :: ERR_STRING_ON
         INTEGER(KIND=4)              :: SetErrorStringOn
        END FUNCTION SetErrorStringOn
       END INTERFACE


       INTERFACE
        FUNCTION SetLogFileName(ID,FNAME)
         INTEGER(KIND=4),  INTENT(IN) :: ID
         CHARACTER(LEN=*), INTENT(IN) :: FNAME
         INTEGER(KIND=4)              :: SetLogFileName
        END FUNCTION SetLogFileName
       END INTERFACE


       INTERFACE
        FUNCTION SetLogFileOn(ID,LOG_ON)
         INTEGER(KIND=4),  INTENT(IN) :: ID
         LOGICAL(KIND=4),  INTENT(IN) :: LOG_ON
         INTEGER(KIND=4)              :: SetLogFileOn
        END FUNCTION SetLogFileOn
       END INTERFACE


       INTERFACE
        FUNCTION SetLogStringOn(ID,LOG_STRING_ON)
         INTEGER(KIND=4),  INTENT(IN) :: ID
         LOGICAL(KIND=4),  INTENT(IN) :: LOG_STRING_ON
         INTEGER(KIND=4)              :: SetLogStringOn
        END FUNCTION SetLogStringOn
       END INTERFACE


       INTERFACE
        FUNCTION SetOutputFileName(ID,FNAME)
         INTEGER(KIND=4),  INTENT(IN) :: ID
         CHARACTER(LEN=*), INTENT(IN) :: FNAME
         INTEGER(KIND=4)              :: SetOutputFileName
        END FUNCTION SetOutputFileName
       END INTERFACE


       INTERFACE
        FUNCTION SetOutputFileOn(ID,OUT_ON)
         INTEGER(KIND=4),  INTENT(IN) :: ID
         LOGICAL(KIND=4),  INTENT(IN) :: OUT_ON
         INTEGER(KIND=4)              :: SetOutputFileOn
        END FUNCTION SetOutputFileOn
       END INTERFACE


       INTERFACE
        FUNCTION SetOutputStringOn(ID,OUT_STRING_ON)
         INTEGER(KIND=4),  INTENT(IN) :: ID
         LOGICAL(KIND=4),  INTENT(IN) :: OUT_STRING_ON
         INTEGER(KIND=4)              :: SetOutputStringOn
        END FUNCTION SetOutputStringOn
       END INTERFACE


       INTERFACE
        FUNCTION SetSelectedOutputFileName(ID,FNAME)
         INTEGER(KIND=4),  INTENT(IN) :: ID
         CHARACTER(LEN=*), INTENT(IN) :: FNAME
         INTEGER(KIND=4)              :: SetSelectedOutputFileName
        END FUNCTION SetSelectedOutputFileName
       END INTERFACE


       INTERFACE
        FUNCTION SetSelectedOutputFileOn(ID,SEL_ON)
         INTEGER(KIND=4),  INTENT(IN) :: ID
         LOGICAL(KIND=4),  INTENT(IN) :: SEL_ON
         INTEGER(KIND=4)              :: SetSelectedOutputFileOn
        END FUNCTION SetSelectedOutputFileOn
       END INTERFACE


       INTERFACE
        FUNCTION SetSelectedOutputStringOn(ID,SEL_ON)
         INTEGER(KIND=4),  INTENT(IN) :: ID
         LOGICAL(KIND=4),  INTENT(IN) :: SEL_ON
         INTEGER(KIND=4)              :: SetSelectedOutputStringOn
        END FUNCTION SetSelectedOutputStringOn
       END INTERFACE
