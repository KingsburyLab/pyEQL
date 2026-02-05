!
!     GetSelectedOutputValue TYPES
!      
      INTEGER :: TT_EMPTY
      PARAMETER (TT_EMPTY=0)
          
      INTEGER :: TT_ERROR
      PARAMETER (TT_ERROR=1)
          
      INTEGER :: TT_DOUBLE
      PARAMETER (TT_DOUBLE=3)
          
      INTEGER :: TT_STRING
      PARAMETER (TT_STRING=4)
          
!
!     IPQ_RESULT TYPES
!      
          
      INTEGER :: IPQ_OK
      PARAMETER (IPQ_OK=0)
          
      INTEGER :: IPQ_OUTOFMEMORY
      PARAMETER (IPQ_OUTOFMEMORY=-1)
         
      INTEGER :: IPQ_BADVARTYPE
      PARAMETER (IPQ_BADVARTYPE=-2)
          
      INTEGER :: IPQ_INVALIDARG
      PARAMETER (IPQ_INVALIDARG=-3)
          
      INTEGER :: IPQ_INVALIDROW
      PARAMETER (IPQ_INVALIDROW=-4)
          
      INTEGER :: IPQ_INVALIDCOL
      PARAMETER (IPQ_INVALIDCOL=-5)
          
      INTEGER :: IPQ_BADINSTANCE
      PARAMETER (IPQ_BADINSTANCE=-6)
          
!
!     function prototypes
!      
          
      INTEGER(KIND=4) AccumulateLine
      INTEGER(KIND=4) AddError
      INTEGER(KIND=4) AddWarning
      INTEGER(KIND=4) ClearAccumulatedLines
      INTEGER(KIND=4) CreateIPhreeqc
      INTEGER(KIND=4) DestroyIPhreeqc
      INTEGER(KIND=4) GetComponentCount
      INTEGER(KIND=4) GetCurrentSelectedOutputUserNumber
      LOGICAL(KIND=4) GetDumpFileOn
      INTEGER(KIND=4) GetDumpStringLineCount
      LOGICAL(KIND=4) GetDumpStringOn
      LOGICAL(KIND=4) GetErrorFileOn
      LOGICAL(KIND=4) GetErrorOn
      INTEGER(KIND=4) GetErrorStringLine
      INTEGER(KIND=4) GetErrorStringLineCount
      LOGICAL(KIND=4) GetLogFileOn      
      LOGICAL(KIND=4) GetLogStringLineCount
      LOGICAL(KIND=4) GetLogStringOn
      LOGICAL(KIND=4) GetOutputFileOn
      INTEGER(KIND=4) GetOutputStringLineCount
      LOGICAL(KIND=4) GetOutputStringOn
      INTEGER(KIND=4) GetSelectedOutputColumnCount
      LOGICAL(KIND=4) GetSelectedOutputFileOn
      INTEGER(KIND=4) GetSelectedOutputRowCount
      INTEGER(KIND=4) GetSelectedOutputStringLineCount
      INTEGER(KIND=4) GetSelectedOutputValue
      INTEGER(KIND=4) GetWarningStringLine
      INTEGER(KIND=4) GetWarningStringLineCount
      INTEGER(KIND=4) LoadDatabase
      INTEGER(KIND=4) LoadDatabaseString
      INTEGER(KIND=4) RunAccumulated
      INTEGER(KIND=4) RunFile
      INTEGER(KIND=4) RunString
      INTEGER(KIND=4) SetDumpFileName
      INTEGER(KIND=4) SetDumpFileOn
      INTEGER(KIND=4) SetDumpStringOn
      INTEGER(KIND=4) SetErrorFileOn
      INTEGER(KIND=4) SetErrorOn
      INTEGER(KIND=4) SetErrorStringOn
      INTEGER(KIND=4) SetLogFileName
      INTEGER(KIND=4) SetLogFileOn
      INTEGER(KIND=4) SetLogStringOn
      INTEGER(KIND=4) SetOutputFileName
      INTEGER(KIND=4) SetOutputFileOn
      INTEGER(KIND=4) SetOutputStringOn
      INTEGER(KIND=4) SetSelectedOutputFileOn
      INTEGER(KIND=4) SetSelectedOutputStringOn
