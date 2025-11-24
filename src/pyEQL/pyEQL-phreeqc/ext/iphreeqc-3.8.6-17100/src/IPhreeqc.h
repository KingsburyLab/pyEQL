/*! @file IPhreeqc.h
	@brief C/Fortran Documentation
*/
#ifndef INC_IPHREEQC_H
#define INC_IPHREEQC_H

#include "Var.h"

#ifdef IPHREEQC_NO_FORTRAN_MODULE
#include <stddef.h>
#endif

/**
 * @mainpage IPhreeqc Library Documentation (3.8.6-17100)
 *
 *  @htmlonly
 *  <table>
 *   <tr><td class="indexkey"><a class="el" href="IPhreeqc_8h.html">IPhreeqc.h</a> </td><td class="indexvalue">C/Fortran Documentation </td></tr>
 *   <tr><td class="indexkey"><a class="el" href="IPhreeqc_8hpp.html">IPhreeqc.hpp</a> </td><td class="indexvalue">C++ Documentation </td></tr>
 *   <tr><td class="indexkey"><a class="el" href="Var_8h.html">Var.h</a></td><td class="indexvalue">IPhreeqc VARIANT Documentation </td></tr>
 *  </table>
 *  @endhtmlonly
 */

/*! @brief Enumeration used to return error codes.
*/
typedef enum {
	IPQ_OK            =  0,  /*!< Success */
	IPQ_OUTOFMEMORY   = -1,  /*!< Failure, Out of memory */
	IPQ_BADVARTYPE    = -2,  /*!< Failure, Invalid VAR type */
	IPQ_INVALIDARG    = -3,  /*!< Failure, Invalid argument */
	IPQ_INVALIDROW    = -4,  /*!< Failure, Invalid row */
	IPQ_INVALIDCOL    = -5,  /*!< Failure, Invalid column */
	IPQ_BADINSTANCE   = -6   /*!< Failure, Invalid instance id */
} IPQ_RESULT;


#if defined(__cplusplus)
extern "C" {
#endif

/**
 *  Accumlulate line(s) for input to phreeqc.
 *  @param id            The instance id returned from @ref CreateIPhreeqc.
 *  @param line          The line(s) to add for input to phreeqc.
 *  @retval IPQ_OK Success
 *  @retval IPQ_OUTOFMEMORY Out of memory
 *  @see                 ClearAccumulatedLines, OutputAccumulatedLines, RunAccumulated
 *  @par Fortran90 Interface:
 *  @htmlonly
 *  <CODE>
 *  <PRE>
 *  FUNCTION AccumulateLine(ID,LINE)
 *    INTEGER(KIND=4),   INTENT(IN)  :: ID
 *    CHARACTER(LEN=*),  INTENT(IN)  :: LINE
 *    INTEGER(KIND=4)                :: AccumulateLine
 *  END FUNCTION AccumulateLine
 *  </PRE>
 *  </CODE>
 *  @endhtmlonly
 *
 *  @par C Example:
 *  @include AccumulateLine.c
 *
 *  @par Fortran90 Example:
 *  see @ref GetDumpStringLine_f90 "GetDumpStringLine"
 */
	IPQ_DLL_EXPORT IPQ_RESULT  AccumulateLine(int id, const char *line);


/**
 *  Appends the given error message and increments the error count.
 *  Internally used to create an error condition.
 *  @param id            The instance id returned from @ref CreateIPhreeqc.
 *  @param error_msg     The error message to display.
 *  @return              The current error count if successful; otherwise a negative value indicates an error occurred (see @ref IPQ_RESULT).
 *  @see                 GetErrorString, GetErrorStringLine, GetErrorStringLineCount, OutputErrorString
 *  @par Fortran90 Interface:
 *  @htmlonly
 *  <CODE>
 *  <PRE>
 *  FUNCTION AddError(ID,ERROR_MSG)
 *    INTEGER(KIND=4),  INTENT(IN) :: ID
 *    CHARACTER(LEN=*), INTENT(IN) :: ERROR_MSG
 *    INTEGER(KIND=4)              :: AddError
 *  END FUNCTION AddError
 *  </PRE>
 *  </CODE>
 *  @endhtmlonly
 */
	IPQ_DLL_EXPORT int         AddError(int id, const char* error_msg);


/**
 *  Appends the given warning message and increments the warning count.
 *  Internally used to create a warning condition.
 *  @param id            The instance id returned from @ref CreateIPhreeqc.
 *  @param warn_msg      The warning message to display.
 *  @return              The current warning count if successful; otherwise a negative value indicates an error occurred (see @ref IPQ_RESULT).
 *  @see                 GetWarningString, GetWarningStringLine, GetWarningStringLineCount, OutputWarningString
 *  @par Fortran90 Interface:
 *  @htmlonly
 *  <CODE>
 *  <PRE>
 *  FUNCTION AddWarning(ID,WARN_MSG)
 *    INTEGER(KIND=4),  INTENT(IN) :: ID
 *    CHARACTER(LEN=*), INTENT(IN) :: WARN_MSG
 *    INTEGER(KIND=4)              :: AddWarning
 *  END FUNCTION AddWarning
 *  </PRE>
 *  </CODE>
 *  @endhtmlonly
 */
	IPQ_DLL_EXPORT int         AddWarning(int id, const char* warn_msg);



/**
 *  Clears the accumulated input buffer.  Input buffer is accumulated from calls to @ref AccumulateLine.
 *  @retval IPQ_OK           Success.
 *  @retval IPQ_BADINSTANCE  The given id is invalid.
 *  @see                     AccumulateLine, OutputAccumulatedLines, RunAccumulated
 *  @par Fortran90 Interface:
 *  @htmlonly
 *  <CODE>
 *  <PRE>
 *  FUNCTION ClearAccumulatedLines(ID)
 *    INTEGER(KIND=4), INTENT(IN) :: ID
 *    INTEGER(KIND=4)             :: ClearAccumulatedLines
 *  END FUNCTION ClearAccumulatedLines
 *  </PRE>
 *  </CODE>
 *  @endhtmlonly
 */
	IPQ_DLL_EXPORT IPQ_RESULT  ClearAccumulatedLines(int id);


/**
 *  Create a new IPhreeqc instance.
 *  @return      A non-negative value if successful; otherwise a negative value indicates an error occurred (see @ref IPQ_RESULT).
 *  @see         DestroyIPhreeqc
 *  @par Fortran90 Interface:
 *  @htmlonly
 *  <CODE>
 *  <PRE>
 *  FUNCTION CreateIPhreeqc()
 *    INTEGER(KIND=4)  :: CreateIPhreeqc
 *  END FUNCTION CreateIPhreeqc
 *  </PRE>
 *  </CODE>
 *  @endhtmlonly
 *
 *  @anchor CreateIPhreeqc_c
 *  @par C Example:
 *  @include CreateIPhreeqc.c
 *
 *  @anchor CreateIPhreeqc_f90
 *  @par Fortran90 Example:
 *  @include F90CreateIPhreeqc.f90
 */
	IPQ_DLL_EXPORT int         CreateIPhreeqc(void);


/**
 *  Release an IPhreeqc instance from memory.
 *  @param id            The instance id returned from @ref CreateIPhreeqc.
 *  @retval IPQ_OK Success
 *  @retval IPQ_BADINSTANCE The given id is invalid.
 *  @see                 CreateIPhreeqc
 *  @par Fortran90 Interface:
 *  @htmlonly
 *  <CODE>
 *  <PRE>
 *  FUNCTION DestroyIPhreeqc(ID)
 *    INTEGER(KIND=4),  INTENT(IN)  :: ID
 *    INTEGER(KIND=4)               :: DestroyIPhreeqc
 *  END FUNCTION DestroyIPhreeqc
 *  </PRE>
 *  </CODE>
 *  @endhtmlonly
 *
 *  @par C Example:
 *  see @ref CreateIPhreeqc_c "CreateIPhreeqc"
 *
 *  @par Fortran90 Example:
 *  see @ref CreateIPhreeqc_f90 "CreateIPhreeqc"
 */
	IPQ_DLL_EXPORT IPQ_RESULT  DestroyIPhreeqc(int id);


/**
 *  Retrieves the given component.
 *  @param id            The instance id returned from @ref CreateIPhreeqc.
 *  @param n             The zero-based index of the component to retrieve.
 *  @return              A null terminated string containing the given component.
 *                       Returns an empty string if n is out of range.
 *  @see                 GetComponentCount
 *  @par Fortran90 Interface:
 *  (Note: N is one-based for the Fortran interface)
 *  @htmlonly
 *  <CODE>
 *  <PRE>
 *  SUBROUTINE GetComponent(ID,N,COMP)
 *    INTEGER(KIND=4),   INTENT(IN)   :: ID
 *    INTEGER(KIND=4),   INTENT(IN)   :: N
 *    CHARACTER(LEN=*),  INTENT(OUT)  :: COMP
 *  END SUBROUTINE GetComponent
 *  </PRE>
 *  </CODE>
 *  @endhtmlonly
 *
 *  @anchor GetComponent_c
 *  @par  C Example:
 *  @include GetComponent.c
 *
 *  @anchor GetComponent_f90
 *  @par  Fortran90 Example:
 *  @include F90GetComponent.f90
 */
	IPQ_DLL_EXPORT const char* GetComponent(int id, int n);


/**
 *  Retrieves the number of components in the current component list.
 *  @param id            The instance id returned from @ref CreateIPhreeqc.
 *  @return              The current count of components.
 *                       A negative value indicates an error occurred (see @ref IPQ_RESULT).
 *  @see                 GetComponent
 *  @par Fortran90 Interface:
 *  @htmlonly
 *  <CODE>
 *  <PRE>
 *  FUNCTION GetComponentCount(ID)
 *    INTEGER(KIND=4),  INTENT(IN)  :: ID
 *    INTEGER(KIND=4)               :: GetComponentCount
 *  END FUNCTION GetComponentCount
 *  </PRE>
 *  </CODE>
 *  @endhtmlonly
 *
 *  @par C Example:
 *  see @ref GetComponent_c "GetComponent"
 *
 *  @par Fortran90 Example:
 *  see @ref GetComponent_f90 "GetComponent"
 */
	IPQ_DLL_EXPORT int         GetComponentCount(int id);

/**
 *  Retrieves the current <b>SELECTED_OUTPUT</b> user number for use in subsequent calls to (@ref GetSelectedOutputColumnCount,
 *  GetSelectedOutputFileName, GetSelectedOutputRowCount, GetSelectedOutputString, GetSelectedOutputStringLine,
 *  GetSelectedOutputStringLineCount, GetSelectedOutputValue, GetSelectedOutputValue2) routines.
 *  The initial setting after calling @ref CreateIPhreeqc is 1.
 *  @param id                   The instance id returned from @ref CreateIPhreeqc.
 *  @return                     The current <b>SELECTED_OUTPUT</b> user number.
 *  @see                        GetNthSelectedOutputUserNumber, GetSelectedOutputCount, SetCurrentSelectedOutputUserNumber
 *  @par Fortran90 Interface:
 *  @htmlonly
 *  <CODE>
 *  <PRE>
 *  FUNCTION GetCurrentSelectedOutputUserNumber(ID)
 *    INTEGER(KIND=4),  INTENT(IN)  :: ID
 *    INTEGER(KIND=4)               :: GetCurrentSelectedOutputUserNumber
 *  END FUNCTION GetCurrentSelectedOutputUserNumber
 *  </PRE>
 *  </CODE>
 *  @endhtmlonly
 */
	IPQ_DLL_EXPORT int         GetCurrentSelectedOutputUserNumber(int id);

/**
 *  Retrieves the name of the dump file.  This file name is used if not specified within <B>DUMP</B> input.
 *  The default value is <B><I>dump.id.out</I></B>.
 *  @param id               The instance id returned from @ref CreateIPhreeqc.
 *  @return filename        The name of the file to write <B>DUMP</B> output to.
 *  @see                    GetDumpFileOn, GetDumpString, GetDumpStringOn, GetDumpStringLine, GetDumpStringLineCount, SetDumpFileName, SetDumpFileOn, SetDumpStringOn
 *  @par Fortran90 Interface:
 *  @htmlonly
 *  <CODE>
 *  <PRE>
 *  SUBROUTINE GetDumpFileName(ID,FILENAME)
 *    INTEGER(KIND=4),   INTENT(IN)   :: ID
 *    CHARACTER(LEN=*),  INTENT(OUT)  :: FILENAME
 *  END SUBROUTINE GetDumpFileName
 *  </PRE>
 *  </CODE>
 *  @endhtmlonly
 */
	IPQ_DLL_EXPORT const char* GetDumpFileName(int id);


/**
 *  Retrieves the current value of the dump file switch.
 *  @param id            The instance id returned from @ref CreateIPhreeqc.
 *  @return              Non-zero if output is written to the <B>DUMP</B> (<B><I>dump.id.out</I></B> if unspecified) file, 0 (zero) otherwise.
 *  @see                 GetDumpString, GetDumpStringLine, GetDumpStringLineCount, GetDumpStringOn, SetDumpFileOn, SetDumpStringOn
 *  @par Fortran90 Interface:
 *  @htmlonly
 *  <CODE>
 *  <PRE>
 *  FUNCTION GetDumpFileOn(ID)
 *    INTEGER(KIND=4),  INTENT(IN)  :: ID
 *    LOGICAL(KIND=4)               :: GetDumpFileOn
 *  END FUNCTION GetDumpFileOn
 *  </PRE>
 *  </CODE>
 *  @endhtmlonly
 */
	IPQ_DLL_EXPORT int         GetDumpFileOn(int id);


/**
 *  Retrieves the string buffer containing <b>DUMP</b> output.
 *  @param id            The instance id returned from @ref CreateIPhreeqc.
 *  @return              A null terminated string containing <b>DUMP</b> output.
 *  @pre                 @ref SetDumpStringOn must have been set to true (non-zero) in order to receive <b>DUMP</b> output.
 *  @see                 GetDumpFileOn, GetDumpStringLine, GetDumpStringLineCount, SetDumpFileOn, GetDumpStringOn, SetDumpStringOn
 *  @par Fortran90 Interface:
 *  Not implemented. (see @ref GetDumpStringLineCount, @ref GetDumpStringLine)
 *
 *  @anchor GetDumpString_c
 *  @par  C Example:
 *  @include GetDumpString.c
 */
	IPQ_DLL_EXPORT const char* GetDumpString(int id);


/**
 *  Retrieves the given dump line.
 *  @param id            The instance id returned from @ref CreateIPhreeqc.
 *  @param n             The zero-based index of the line to retrieve.
 *  @return              A null terminated string containing the given line.
 *                       Returns an empty string if n is out of range.
 *  @pre                 @ref SetDumpStringOn must have been set to true (non-zero).
 *  @see                 GetDumpFileOn, GetDumpString, GetDumpStringLineCount, GetDumpStringOn, SetDumpFileOn, SetDumpStringOn
 *  @par Fortran90 Interface:
 *  @htmlonly
 *  (Note: N is one-based for the Fortran interface.)
 *  <CODE>
 *  <PRE>
 *  SUBROUTINE GetDumpStringLine(ID,N,LINE)
 *    INTEGER(KIND=4),   INTENT(IN)   :: ID
 *    INTEGER(KIND=4),   INTENT(IN)   :: N
 *    CHARACTER(LEN=*),  INTENT(OUT)  :: LINE
 *  END SUBROUTINE GetDumpStringLine
 *  </PRE>
 *  </CODE>
 *  @endhtmlonly
 *
 *  @anchor GetDumpStringLine_f90
 *  @par  Fortran90 Example:
 *  @include F90GetDumpStringLine.f90
 */
	IPQ_DLL_EXPORT const char* GetDumpStringLine(int id, int n);


/**
 *  Retrieves the number of lines in the current dump string buffer.
 *  @param id            The instance id returned from @ref CreateIPhreeqc.
 *  @return              The number of lines.
 *  @pre                 @ref SetDumpStringOn must have been set to true (non-zero).
 *  @see                 GetDumpFileOn, GetDumpString, GetDumpStringLine, GetDumpStringOn, SetDumpFileOn, SetDumpStringOn
 *  @par Fortran90 Interface:
 *  @htmlonly
 *  <CODE>
 *  <PRE>
 *  FUNCTION GetDumpStringLineCount(ID)
 *    INTEGER(KIND=4),  INTENT(IN)  :: ID
 *    INTEGER(KIND=4)               :: GetDumpStringLineCount
 *  END FUNCTION GetDumpStringLineCount
 *  </PRE>
 *  </CODE>
 *  @endhtmlonly
 *
 *  @par Fortran90 Example:
 *  see @ref GetDumpStringLine_f90 "GetDumpStringLine"
 */
	IPQ_DLL_EXPORT int         GetDumpStringLineCount(int id);


/**
 *  Retrieves the current value of the dump string switch.
 *  @param id            The instance id returned from @ref CreateIPhreeqc.
 *  @return              Non-zero if output defined by the <B>DUMP</B> keyword is stored, 0 (zero) otherwise.
 *  @see                 GetDumpFileOn, GetDumpString, GetDumpStringLine, GetDumpStringLineCount, SetDumpFileOn, SetDumpStringOn
 *  @par Fortran90 Interface:
 *  @htmlonly
 *  <CODE>
 *  <PRE>
 *  FUNCTION GetDumpStringOn(ID)
 *    INTEGER(KIND=4),  INTENT(IN)  :: ID
 *    LOGICAL(KIND=4)               :: GetDumpStringOn
 *  END FUNCTION GetDumpStringOn
 *  </PRE>
 *  </CODE>
 *  @endhtmlonly
 */
	IPQ_DLL_EXPORT int         GetDumpStringOn(int id);


/**
 *  Retrieves the name of the error file.  The default name is <B><I>phreeqc.id.err</I></B>.
 *  @param id               The instance id returned from @ref CreateIPhreeqc.
 *  @return filename        The name of the error file.
 *  @see                    GetErrorFileOn, GetErrorString, GetErrorStringOn, GetErrorStringLine, GetErrorStringLineCount, SetErrorFileName, SetErrorFileOn, SetErrorStringOn
 *  @par Fortran90 Interface:
 *  @htmlonly
 *  <CODE>
 *  <PRE>
 *  SUBROUTINE GetErrorFileName(ID,FILENAME)
 *    INTEGER(KIND=4),   INTENT(IN)   :: ID
 *    CHARACTER(LEN=*),  INTENT(OUT)  :: FILENAME
 *  END SUBROUTINE GetErrorFileName
 *  </PRE>
 *  </CODE>
 *  @endhtmlonly
 */
	IPQ_DLL_EXPORT const char* GetErrorFileName(int id);


/**
 *  Retrieves the current value of the error file switch.
 *  @param id            The instance id returned from @ref CreateIPhreeqc.
 *  @return              Non-zero if errors are written to the <B><I>phreeqc.id.err</I></B> file, 0 (zero) otherwise.
 *  @see                 SetErrorFileOn
 *  @par Fortran90 Interface:
 *  @htmlonly
 *  <CODE>
 *  <PRE>
 *  FUNCTION GetErrorFileOn(ID)
 *    INTEGER(KIND=4),  INTENT(IN)  :: ID
 *    LOGICAL(KIND=4)               :: GetErrorFileOn
 *  END FUNCTION GetErrorFileOn
 *  </PRE>
 *  </CODE>
 *  @endhtmlonly
 */
	IPQ_DLL_EXPORT int         GetErrorFileOn(int id);

/**
 *  Retrieves the current value of the error on switch.
 *  @param id            The instance id returned from @ref CreateIPhreeqc.
 *  @return              Non-zero if errors are generated, 0 (zero) otherwise.
 *  @see                 SetErrorOn
 *  @par Fortran90 Interface:
 *  @htmlonly
 *  <CODE>
 *  <PRE>
 *  FUNCTION GetErrorOn(ID)
 *    INTEGER(KIND=4),  INTENT(IN)  :: ID
 *    LOGICAL(KIND=4)               :: GetErrorOn
 *  END FUNCTION GetErrorOn
 *  </PRE>
 *  </CODE>
 *  @endhtmlonly
 */
	IPQ_DLL_EXPORT int         GetErrorOn(int id);


/**
 *  Retrieves the error messages from the last call to @ref RunAccumulated, @ref RunFile, @ref RunString, @ref LoadDatabase, or @ref LoadDatabaseString.
 *  @param id            The instance id returned from @ref CreateIPhreeqc.
 *  @return              A null terminated string containing error messages.
 *  @see                 GetErrorFileOn, GetErrorStringLine, GetErrorStringLineCount, OutputErrorString, SetErrorFileOn
 *  @par Fortran90 Interface:
 *  Not implemented. (see @ref GetErrorStringLineCount, @ref GetErrorStringLine, @ref OutputErrorString)
 *
 *  @anchor GetErrorString_c
 *  @par  C Example:
 *  @include GetErrorString.c
 */
	IPQ_DLL_EXPORT const char* GetErrorString(int id);


/**
 *  Retrieves the given error line.
 *  @param id            The instance id returned from @ref CreateIPhreeqc.
 *  @param n             The zero-based index of the line to retrieve.
 *  @return              A null terminated string containing the given line of the error string buffer.
 *  @see                 GetErrorFileOn, GetErrorString, GetErrorStringLineCount, OutputErrorString, SetErrorFileOn
 *  @par Fortran90 Interface:
 *  (Note: N is one-based for the Fortran interface.)
 *  @htmlonly
 *  <CODE>
 *  <PRE>
 *  SUBROUTINE GetErrorStringLine(ID,N,LINE)
 *    INTEGER(KIND=4),   INTENT(IN)   :: ID
 *    INTEGER(KIND=4),   INTENT(IN)   :: N
 *    CHARACTER(LEN=*),  INTENT(OUT)  :: LINE
 *  END SUBROUTINE GetErrorStringLine
 *  </PRE>
 *  </CODE>
 *  @endhtmlonly
 *
 *  @anchor GetErrorStringLine_f90
 *  @par  Fortran90 Example:
 *  @include F90GetErrorStringLine.f90
 */
	IPQ_DLL_EXPORT const char* GetErrorStringLine(int id, int n);


/**
 *  Retrieves the number of lines in the current error string buffer.
 *  @param id            The instance id returned from @ref CreateIPhreeqc.
 *  @return              The number of lines.
 *  @see                 GetErrorFileOn, GetErrorString, GetErrorStringLine, OutputErrorString, SetErrorFileOn
 *  @par Fortran90 Interface:
 *  @htmlonly
 *  <CODE>
 *  <PRE>
 *  FUNCTION GetErrorStringLineCount(ID)
 *    INTEGER(KIND=4),  INTENT(IN)  :: ID
 *    INTEGER(KIND=4)               :: GetErrorStringLineCount
 *  END FUNCTION GetErrorStringLineCount
 *  </PRE>
 *  </CODE>
 *  @endhtmlonly
 */
	IPQ_DLL_EXPORT int         GetErrorStringLineCount(int id);

/**
 *  Retrieves the current value of the error string switch.
 *  @param id            The instance id returned from @ref CreateIPhreeqc.
 *  @return              Non-zero if output is stored, 0 (zero) otherwise.
 *  @see                 GetErrorFileOn, GetErrorString, GetErrorStringLine, GetErrorStringLineCount, SetErrorFileOn, SetErrorStringOn
 *  @par Fortran90 Interface:
 *  @htmlonly
 *  <CODE>
 *  <PRE>
 *  FUNCTION GetErrorStringOn(ID)
 *    INTEGER(KIND=4),  INTENT(IN)  :: ID
 *    LOGICAL(KIND=4)               :: GetErrorStringOn
 *  END FUNCTION GetErrorStringOn
 *  </PRE>
 *  </CODE>
 *  @endhtmlonly
 */
	IPQ_DLL_EXPORT int         GetErrorStringOn(int id);

/**
 *  Retrieves the name of the log file.  The default name is <B><I>phreeqc.id.log</I></B>.
 *  @param id               The instance id returned from @ref CreateIPhreeqc.
 *  @return filename        The name of the log file.
 *  @see                    GetLogFileOn, GetLogString, GetLogStringOn, GetLogStringLine, GetLogStringLineCount, SetLogFileName, SetLogFileOn, SetLogStringOn
 *  @par Fortran90 Interface:
 *  @htmlonly
 *  <CODE>
 *  <PRE>
 *  SUBROUTINE GetLogFileName(ID,FILENAME)
 *    INTEGER(KIND=4),   INTENT(IN)   :: ID
 *    CHARACTER(LEN=*),  INTENT(OUT)  :: FILENAME
 *  END SUBROUTINE GetLogFileName
 *  </PRE>
 *  </CODE>
 *  @endhtmlonly
 */
	IPQ_DLL_EXPORT const char* GetLogFileName(int id);


/**
 *  Retrieves the current value of the log file switch.
 *  @param id            The instance id returned from @ref CreateIPhreeqc.
 *  @return              Non-zero if log messages are written to the <B><I>phreeqc.id.log</I></B> file, 0 (zero) otherwise.
 *  @remarks             Logging must be enabled through the use of the KNOBS -logfile option in order to receive any log messages.
 *  @see                 SetLogFileOn
 *  @par Fortran90 Interface:
 *  @htmlonly
 *  <CODE>
 *  <PRE>
 *  FUNCTION GetLogFileOn(ID)
 *    INTEGER(KIND=4),  INTENT(IN)  :: ID
 *    LOGICAL(KIND=4)               :: GetLogFileOn
 *  END FUNCTION GetLogFileOn
 *  </PRE>
 *  </CODE>
 *  @endhtmlonly
 */
	IPQ_DLL_EXPORT int         GetLogFileOn(int id);


/**
 *  Retrieves the string buffer containing log output.
 *  @param id            The instance id returned from @ref CreateIPhreeqc.
 *  @return              A null terminated string containing log output.
 *  @remarks             Logging must be enabled through the use of the KNOBS -logfile option in order to receive any log messages.
 *  @pre                 @ref SetLogStringOn must have been set to true (non-zero) in order to receive log output.
 *  @see                 GetLogFileOn, GetLogStringLine, GetLogStringLineCount, SetLogFileOn, GetLogStringOn, SetLogStringOn
 *  @par Fortran90 Interface:
 *  Not implemented. (see @ref GetLogStringLineCount, @ref GetLogStringLine)
 *
 *  @anchor GetLogString_c
 *  @par  C Example:
 *  @include GetLogString.c
 */
	IPQ_DLL_EXPORT const char* GetLogString(int id);


/**
 *  Retrieves the given log line.
 *  @param id            The instance id returned from @ref CreateIPhreeqc.
 *  @param n             The zero-based index of the line to retrieve.
 *  @return              A null terminated string containing the given line.
 *                       Returns an empty string if n is out of range.
 *  @pre                 @ref SetLogStringOn must have been set to true (non-zero).
 *  @see                 GetLogFileOn, GetLogString, GetLogStringLineCount, GetLogStringOn, SetLogFileOn, SetLogStringOn
 *  @par Fortran90 Interface:
 *  @htmlonly
 *  (Note: N is one-based for the Fortran interface.)
 *  <CODE>
 *  <PRE>
 *  SUBROUTINE GetLogStringLine(ID,N,LINE)
 *    INTEGER(KIND=4),   INTENT(IN)   :: ID
 *    INTEGER(KIND=4),   INTENT(IN)   :: N
 *    CHARACTER(LEN=*),  INTENT(OUT)  :: LINE
 *  END SUBROUTINE GetLogStringLine
 *  </PRE>
 *  </CODE>
 *  @endhtmlonly
 *
 *  @anchor GetLogStringLine_f90
 *  @par  Fortran90 Example:
 *  @include F90GetLogStringLine.f90
 */
	IPQ_DLL_EXPORT const char* GetLogStringLine(int id, int n);

/**
 *  Retrieves the number of lines in the current log string buffer.
 *  @param id            The instance id returned from @ref CreateIPhreeqc.
 *  @return              The number of lines.
 *  @pre                 @ref SetLogStringOn must have been set to true (non-zero).
 *  @see                 GetLogFileOn, GetLogString, GetLogStringLine, GetLogStringOn, SetLogFileOn, SetLogStringOn
 *  @par Fortran90 Interface:
 *  @htmlonly
 *  <CODE>
 *  <PRE>
 *  FUNCTION GetLogStringLineCount(ID)
 *    INTEGER(KIND=4),  INTENT(IN)  :: ID
 *    INTEGER(KIND=4)               :: GetLogStringLineCount
 *  END FUNCTION GetLogStringLineCount
 *  </PRE>
 *  </CODE>
 *  @endhtmlonly
 *
 *  @par Fortran90 Example:
 *  see @ref GetLogStringLine_f90 "GetLogStringLine"
 */
	IPQ_DLL_EXPORT int         GetLogStringLineCount(int id);


/**
 *  Retrieves the current value of the log string switch.
 *  @param id            The instance id returned from @ref CreateIPhreeqc.
 *  @return              Non-zero if output is stored, 0 (zero) otherwise.
 *  @see                 GetLogFileOn, GetLogString, GetLogStringLine, GetLogStringLineCount, SetLogFileOn, SetLogStringOn
 *  @par Fortran90 Interface:
 *  @htmlonly
 *  <CODE>
 *  <PRE>
 *  FUNCTION GetLogStringOn(ID)
 *    INTEGER(KIND=4),  INTENT(IN)  :: ID
 *    LOGICAL(KIND=4)               :: GetLogStringOn
 *  END FUNCTION GetLogStringOn
 *  </PRE>
 *  </CODE>
 *  @endhtmlonly
 */
	IPQ_DLL_EXPORT int         GetLogStringOn(int id);


/**
 *  Retrieves the nth user number of the currently defined <B>SELECTED_OUTPUT</B> keyword blocks.
 *  @param id            The instance id returned from @ref CreateIPhreeqc.
 *  @param n             The zero-based index of the <B>SELECTED_OUTPUT</B> user number to retrieve.
 *  @return              The nth defined user number; a negative value indicates an error occurred.
 *  @see                 GetCurrentSelectedOutputUserNumber, GetSelectedOutputCount, SetCurrentSelectedOutputUserNumber
 *  @pre @ref RunAccumulated, @ref RunFile, @ref RunString must have been called and returned 0 (zero) errors.
 *  @par Fortran90 Interface:
 *  @htmlonly
 *  (Note: N is one-based for the Fortran interface.)
 *  <CODE>
 *  <PRE>
 *  FUNCTION GetNthSelectedOutputUserNumber(ID,N)
 *    INTEGER(KIND=4),  INTENT(IN)  :: ID
 *    INTEGER(KIND=4),  INTENT(IN)  :: N
 *    INTEGER(KIND=4)               :: GetNthSelectedOutputUserNumber
 *  END FUNCTION GetNthSelectedOutputUserNumber
 *  </PRE>
 *  </CODE>
 *  @endhtmlonly
 *
 *  @par C Example:
 *  see @ref SetCurrentSelectedOutputUserNumber_c "SetCurrentSelectedOutputUserNumber"
 */
	IPQ_DLL_EXPORT int         GetNthSelectedOutputUserNumber(int id, int n);

/**
 *  Retrieves the name of the output file.  The default name is <B><I>phreeqc.id.out</I></B>.
 *  @param id               The instance id returned from @ref CreateIPhreeqc.
 *  @return filename        The name of the output file.
 *  @see                    GetOutputFileOn, GetOutputString, GetOutputStringOn, GetOutputStringLine, GetOutputStringLineCount, SetOutputFileName, SetOutputFileOn, SetOutputStringOn
 *  @par Fortran90 Interface:
 *  @htmlonly
 *  <CODE>
 *  <PRE>
 *  SUBROUTINE GetOutputFileName(ID,FILENAME)
 *    INTEGER(KIND=4),   INTENT(IN)   :: ID
 *    CHARACTER(LEN=*),  INTENT(OUT)  :: FILENAME
 *  END SUBROUTINE GetOutputFileName
 *  </PRE>
 *  </CODE>
 *  @endhtmlonly
 */
	IPQ_DLL_EXPORT const char* GetOutputFileName(int id);


/**
 *  Retrieves the current value of the output file switch.
 *  @param id            The instance id returned from @ref CreateIPhreeqc.
 *  @return              Non-zero if output is written to the <B><I>phreeqc.id.out</I></B> file, 0 (zero) otherwise.
 *  @see                 SetOutputFileOn
 *  @par Fortran90 Interface:
 *  @htmlonly
 *  <CODE>
 *  <PRE>
 *  FUNCTION GetOutputFileOn(ID)
 *    INTEGER(KIND=4),  INTENT(IN)  :: ID
 *    LOGICAL(KIND=4)               :: GetOutputFileOn
 *  END FUNCTION GetOutputFileOn
 *  </PRE>
 *  </CODE>
 *  @endhtmlonly
 */
	IPQ_DLL_EXPORT int         GetOutputFileOn(int id);

/**
 *  Retrieves the string buffer containing phreeqc output.
 *  @param id            The instance id returned from @ref CreateIPhreeqc.
 *  @return              A null terminated string containing phreeqc output.
 *  @pre                 @ref SetOutputStringOn must have been set to true (non-zero) in order to receive phreeqc output.
 *  @see                 GetOutputFileOn, GetOutputStringLine, GetOutputStringLineCount, SetOutputFileOn, GetOutputStringOn, SetOutputStringOn
 *  @par Fortran90 Interface:
 *  Not implemented. (see @ref GetOutputStringLineCount, @ref GetOutputStringLine)
 *
 *  @anchor GetOutputString_c
 *  @par  C Example:
 *  @include GetOutputString.c
 */
	IPQ_DLL_EXPORT const char* GetOutputString(int id);

/**
 *  Retrieves the given output line.
 *  @param id            The instance id returned from @ref CreateIPhreeqc.
 *  @param n             The zero-based index of the line to retrieve.
 *  @return              A null terminated string containing the given line.
 *                       Returns an empty string if n is out of range.
 *  @pre                 @ref SetOutputStringOn must have been set to true (non-zero).
 *  @see                 GetOutputFileOn, GetOutputString, GetOutputStringLineCount, GetOutputStringOn, SetOutputFileOn, SetOutputStringOn
 *  @par Fortran90 Interface:
 *  @htmlonly
 *  (Note: N is one-based for the Fortran interface.)
 *  <CODE>
 *  <PRE>
 *  SUBROUTINE GetOutputStringLine(ID,N,LINE)
 *    INTEGER(KIND=4),   INTENT(IN)   :: ID
 *    INTEGER(KIND=4),   INTENT(IN)   :: N
 *    CHARACTER(LEN=*),  INTENT(OUT)  :: LINE
 *  END SUBROUTINE GetOutputStringLine
 *  </PRE>
 *  </CODE>
 *  @endhtmlonly
 *
 *  @anchor GetOutputStringLine_f90
 *  @par  Fortran90 Example:
 *  @include F90GetOutputStringLine.f90
 */
	IPQ_DLL_EXPORT const char* GetOutputStringLine(int id, int n);

/**
 *  Retrieves the number of lines in the current output string buffer.
 *  @param id            The instance id returned from @ref CreateIPhreeqc.
 *  @return              The number of lines.
 *  @pre                 @ref SetOutputStringOn must have been set to true (non-zero).
 *  @see                 GetOutputFileOn, GetOutputString, GetOutputStringLine, GetOutputStringOn, SetOutputFileOn, SetOutputStringOn
 *  @par Fortran90 Interface:
 *  @htmlonly
 *  <CODE>
 *  <PRE>
 *  FUNCTION GetOutputStringLineCount(ID)
 *    INTEGER(KIND=4),  INTENT(IN)  :: ID
 *    INTEGER(KIND=4)               :: GetOutputStringLineCount
 *  END FUNCTION GetOutputStringLineCount
 *  </PRE>
 *  </CODE>
 *  @endhtmlonly
 *
 *  @par Fortran90 Example:
 *  see @ref GetOutputStringLine_f90 "GetOutputStringLine"
 */
	IPQ_DLL_EXPORT int         GetOutputStringLineCount(int id);

/**
 *  Retrieves the current value of the output string switch.
 *  @param id            The instance id returned from @ref CreateIPhreeqc.
 *  @return              Non-zero if output is stored, 0 (zero) otherwise.
 *  @see                 GetOutputFileOn, GetOutputString, GetOutputStringLine, GetOutputStringLineCount, SetOutputFileOn, SetOutputStringOn
 *  @par Fortran90 Interface:
 *  @htmlonly
 *  <CODE>
 *  <PRE>
 *  FUNCTION GetOutputStringOn(ID)
 *    INTEGER(KIND=4),  INTENT(IN)  :: ID
 *    LOGICAL(KIND=4)               :: GetOutputStringOn
 *  END FUNCTION GetOutputStringOn
 *  </PRE>
 *  </CODE>
 *  @endhtmlonly
 */
	IPQ_DLL_EXPORT int         GetOutputStringOn(int id);


/**
 *  Retrieves the number of columns in the selected-output buffer.
 *  @param id            The instance id returned from @ref CreateIPhreeqc.
 *  @return              The number of columns.
 *  @see                 GetSelectedOutputFileOn, GetSelectedOutputRowCount, GetSelectedOutputValue, SetSelectedOutputFileOn
 *  @par Fortran90 Interface:
 *  @htmlonly
 *  <CODE>
 *  <PRE>
 *  FUNCTION GetSelectedOutputColumnCount(ID)
 *    INTEGER(KIND=4),  INTENT(IN)  :: ID
 *    INTEGER(KIND=4)               :: GetSelectedOutputColumnCount
 *  END FUNCTION GetSelectedOutputColumnCount
 *  </PRE>
 *  </CODE>
 *  @endhtmlonly
 */
	IPQ_DLL_EXPORT int         GetSelectedOutputColumnCount(int id);

/**
 *  Retrieves the count of <B>SELECTED_OUTPUT</B> blocks that are currently defined.
 *  @param id            The instance id returned from @ref CreateIPhreeqc.
 *  @return              The number of <B>SELECTED_OUTPUT</B> blocks.
 *  @see                 GetCurrentSelectedOutputUserNumber, GetNthSelectedOutputUserNumber, SetCurrentSelectedOutputUserNumber
 *  @pre                 (@ref RunAccumulated, @ref RunFile, @ref RunString) must have been called and returned 0 (zero) errors.
 *  @par Fortran90 Interface:
 *  @htmlonly
 *  <CODE>
 *  <PRE>
 *  FUNCTION GetSelectedOutputCount(ID)
 *    INTEGER(KIND=4),  INTENT(IN)  :: ID
 *    INTEGER(KIND=4)               :: GetSelectedOutputCount
 *  END FUNCTION GetSelectedOutputCount
 *  </PRE>
 *  </CODE>
 *  @endhtmlonly
 *
 *  @par C Example:
 *  see @ref SetCurrentSelectedOutputUserNumber_c "SetCurrentSelectedOutputUserNumber"
 */
	IPQ_DLL_EXPORT int         GetSelectedOutputCount(int id);


/**
 *  Retrieves the name of the current selected output file (see @ref SetCurrentSelectedOutputUserNumber).  This file name is used if not specified within <B>SELECTED_OUTPUT</B> input.
 *  The default value is <B><I>selected_n.id.out</I></B>.
 *  @param id               The instance id returned from @ref CreateIPhreeqc.
 *  @return filename        The name of the file to write <B>SELECTED_OUTPUT</B> output to.
 *  @see                    GetSelectedOutputFileOn, GetSelectedOutputString, GetSelectedOutputStringOn, GetSelectedOutputStringLine, GetSelectedOutputStringLineCount, SetCurrentSelectedOutputUserNumber, SetSelectedOutputFileName, SetSelectedOutputFileOn, SetSelectedOutputStringOn
 *  @par Fortran90 Interface:
 *  @htmlonly
 *  <CODE>
 *  <PRE>
 *  SUBROUTINE GetSelectedOutputFileName(ID,FILENAME)
 *    INTEGER(KIND=4),   INTENT(IN)   :: ID
 *    CHARACTER(LEN=*),  INTENT(OUT)  :: FILENAME
 *  END SUBROUTINE GetSelectedOutputFileName
 *  </PRE>
 *  </CODE>
 *  @endhtmlonly
 */
	IPQ_DLL_EXPORT const char* GetSelectedOutputFileName(int id);


/**
 *  Retrieves the current selected-output file switch (see @ref SetCurrentSelectedOutputUserNumber).
 *  @param id                    The instance id returned from @ref CreateIPhreeqc.
 *  @return                      Non-zero if output is written to the selected-output (<B><I>selected_n.id.out</I></B> if unspecified) file, 0 (zero) otherwise.
 *  @see                         GetSelectedOutputColumnCount, GetSelectedOutputRowCount, GetSelectedOutputValue, SetCurrentSelectedOutputUserNumber, SetSelectedOutputFileOn
 *  @par Fortran90 Interface:
 *  @htmlonly
 *  <CODE>
 *  <PRE>
 *  FUNCTION GetSelectedOutputFileOn(ID)
 *    INTEGER(KIND=4),  INTENT(IN)  :: ID
 *    LOGICAL(KIND=4)               :: GetSelectedOutputFileOn
 *  END FUNCTION GetSelectedOutputFileOn
 *  </PRE>
 *  </CODE>
 *  @endhtmlonly
 */
	IPQ_DLL_EXPORT int         GetSelectedOutputFileOn(int id);


/**
 *  Retrieves the number of rows in the current selected-output buffer (see @ref SetCurrentSelectedOutputUserNumber).
 *  @param id            The instance id returned from @ref CreateIPhreeqc.
 *  @return              The number of rows.
 *  @see                 GetSelectedOutputFileOn, GetSelectedOutputColumnCount, GetSelectedOutputValue, SetCurrentSelectedOutputUserNumber, SetSelectedOutputFileOn
 *  @par Fortran90 Interface:
 *  @htmlonly
 *  <CODE>
 *  <PRE>
 *  FUNCTION GetSelectedOutputRowCount(ID)
 *    INTEGER(KIND=4),  INTENT(IN)  :: ID
 *    INTEGER(KIND=4)               :: GetSelectedOutputRowCount
 *  END FUNCTION GetSelectedOutputRowCount
 *  </PRE>
 *  </CODE>
 *  @endhtmlonly
 */
	IPQ_DLL_EXPORT int         GetSelectedOutputRowCount(int id);


/**
 *  Retrieves the string buffer containing the current <b>SELECTED_OUTPUT</b> (see @ref SetCurrentSelectedOutputUserNumber).
 *  @param id            The instance id returned from @ref CreateIPhreeqc.
 *  @return              A null terminated string containing <b>SELECTED_OUTPUT</b>.
 *  @pre                 @ref SetSelectedOutputStringOn must have been set to true (non-zero) in order to receive <b>SELECTED_OUTPUT</b>.
 *  @see                 GetSelectedOutputFileOn, GetSelectedOutputStringLine, GetSelectedOutputStringLineCount, GetSelectedOutputStringOn, SetSelectedOutputFileOn, SetCurrentSelectedOutputUserNumber SetSelectedOutputStringOn
 *  @par Fortran90 Interface:
 *  Not implemented. (see @ref GetSelectedOutputStringLineCount, @ref GetSelectedOutputStringLine)
 *
 *  @anchor GetSelectedOutputString_c
 *  @par  C Example:
 *  @include GetSelectedOutputString.c
 */
	IPQ_DLL_EXPORT const char* GetSelectedOutputString(int id);


/**
 *  Retrieves the given line of the current selected output string (see @ref SetCurrentSelectedOutputUserNumber).
 *  @param id            The instance id returned from @ref CreateIPhreeqc.
 *  @param n             The zero-based index of the line to retrieve.
 *  @return              A null terminated string containing the given line.
 *                       Returns an empty string if n is out of range.
 *  @pre                 @ref SetSelectedOutputStringOn must have been set to true (non-zero).
 *  @see                 GetSelectedOutputFileOn, GetSelectedOutputString, GetSelectedOutputStringLineCount, GetSelectedOutputStringOn, SetCurrentSelectedOutputUserNumber, SetSelectedOutputFileOn, SetSelectedOutputStringOn
 *  @par Fortran90 Interface:
 *  @htmlonly
 *  (Note: N is one-based for the Fortran interface.)
 *  <CODE>
 *  <PRE>
 *  SUBROUTINE GetSelectedOutputStringLine(ID,N,LINE)
 *    INTEGER(KIND=4),   INTENT(IN)   :: ID
 *    INTEGER(KIND=4),   INTENT(IN)   :: N
 *    CHARACTER(LEN=*),  INTENT(OUT)  :: LINE
 *  END SUBROUTINE GetSelectedOutputStringLine
 *  </PRE>
 *  </CODE>
 *  @endhtmlonly
 *
 *  @anchor GetSelectedOutputStringLine_f90
 *  @par  Fortran90 Example:
 *  @include F90GetSelectedOutputStringLine.f90
 */
	IPQ_DLL_EXPORT const char* GetSelectedOutputStringLine(int id, int n);


/**
 *  Retrieves the number of lines in the current selected output string buffer (see @ref SetCurrentSelectedOutputUserNumber).
 *  @param id            The instance id returned from @ref CreateIPhreeqc.
 *  @return              The number of lines.
 *  @pre                 @ref SetSelectedOutputStringOn must have been set to true (non-zero).
 *  @see                 GetSelectedOutputFileOn, GetSelectedOutputString, GetSelectedOutputStringLine, GetSelectedOutputStringOn, SetCurrentSelectedOutputUserNumber, SetSelectedOutputFileOn, SetSelectedOutputStringOn
 *  @par Fortran90 Interface:
 *  @htmlonly
 *  <CODE>
 *  <PRE>
 *  FUNCTION GetSelectedOutputStringLineCount(ID)
 *    INTEGER(KIND=4),  INTENT(IN)  :: ID
 *    INTEGER(KIND=4)               :: GetSelectedOutputStringLineCount
 *  END FUNCTION GetSelectedOutputStringLineCount
 *  </PRE>
 *  </CODE>
 *  @endhtmlonly
 *
 *  @par Fortran90 Example:
 *  see @ref GetSelectedOutputStringLine_f90 "GetSelectedOutputStringLine"
 */
	IPQ_DLL_EXPORT int         GetSelectedOutputStringLineCount(int id);


/**
 *  Retrieves the value of the current selected output string switch (see @ref SetCurrentSelectedOutputUserNumber).
 *  @param id            The instance id returned from @ref CreateIPhreeqc.
 *  @return              Non-zero if output defined by the <B>SELECTED_OUTPUT</B> keyword is stored, 0 (zero) otherwise.
 *  @see                 GetSelectedOutputFileOn, GetSelectedOutputString, GetSelectedOutputStringLine, GetSelectedOutputStringLineCount, SetCurrentSelectedOutputUserNumber, SetSelectedOutputFileOn, SetSelectedOutputStringOn
 *  @par Fortran90 Interface:
 *  @htmlonly
 *  <CODE>
 *  <PRE>
 *  FUNCTION GetSelectedOutputStringOn(ID)
 *    INTEGER(KIND=4),  INTENT(IN)  :: ID
 *    LOGICAL(KIND=4)               :: GetSelectedOutputStringOn
 *  END FUNCTION GetSelectedOutputStringOn
 *  </PRE>
 *  </CODE>
 *  @endhtmlonly
 */
	IPQ_DLL_EXPORT int         GetSelectedOutputStringOn(int id);


/**
 *  Returns the @c VAR associated with the specified row and column.  The current <b>SELECTED_OUTPUT</b> block is set using the @ref SetCurrentSelectedOutputUserNumber method.
 *  @param id                The instance id returned from @ref CreateIPhreeqc.
 *  @param row               The row index.
 *  @param col               The column index.
 *  @param pVAR              Pointer to the @c VAR to receive the requested data.
 *  @retval IPQ_OK           Success.
 *  @retval IPQ_INVALIDROW   The given row is out of range.
 *  @retval IPQ_INVALIDCOL   The given column is out of range.
 *  @retval IPQ_OUTOFMEMORY  Memory could not be allocated.
 *  @retval IPQ_BADINSTANCE  The given id is invalid.
 *  @see                     GetCurrentSelectedOutputUserNumber, GetSelectedOutputFileOn, GetSelectedOutputColumnCount, GetSelectedOutputRowCount, GetSelectedOutputValue2, SetCurrentSelectedOutputUserNumber, SetSelectedOutputFileOn
 *  @remarks
 *  Row 0 contains the column headings to the selected_ouput.
 *  @par Examples:
 *  The headings will include a suffix and/or prefix in order to differentiate the
 *  columns.
 *  @htmlonly
<p>
<table border=1>

<TR VALIGN="top">
<TH width=65%>
Input
</TH>
<TH width=35%>
Headings
</TH>
</TR>

<TR VALIGN="top">
<TD width=65%>
<CODE><PRE>
  SELECTED_OUTPUT
    -reset false
    -totals Ca Na
</PRE></CODE>
</TD>
<TD width=35%>
<CODE><PRE>
  Ca(mol/kgw)  Na(mol/kgw)
</PRE></CODE>
</TD>
</TR>

<TR VALIGN="top">
<TD width=65%>
<CODE><PRE>
  SELECTED_OUTPUT
    -reset false
    -molalities Fe+2 Hfo_sOZn+
</PRE></CODE>
</TD>
<TD width=35%>
<CODE><PRE>
  m_Fe+2(mol/kgw)  m_Hfo_sOZn+(mol/kgw)
</PRE></CODE>
</TD>
</TR>

<TR VALIGN="top">
<TD width=65%>
<CODE><PRE>
  SELECTED_OUTPUT
    -reset false
    -activities H+ Ca+2
</PRE></CODE>
</TD>
<TD width=35%>
<CODE><PRE>
  la_H+  la_Ca+2
</PRE></CODE>
</TD>
</TR>

<TR VALIGN="top">
<TD width=65%>
<CODE><PRE>
  SELECTED_OUTPUT
    -reset false
    -equilibrium_phases Calcite Dolomite
</PRE></CODE>
</TD>
<TD width=35%>
<CODE><PRE>
  Calcite  d_Calcite  Dolomite  d_Dolomite
</PRE></CODE>
</TD>
</TR>

<TR VALIGN="top">
<TD width=65%>
<CODE><PRE>
  SELECTED_OUTPUT
    -reset false
    -saturation_indices CO2(g) Siderite
</PRE></CODE>
</TD>
<TD width=35%>
<CODE><PRE>
  si_CO2(g)  si_Siderite
</PRE></CODE>
</TD>
</TR>

<TR VALIGN="top">
<TD width=65%>
<CODE><PRE>
  SELECTED_OUTPUT
    -reset false
    -gases CO2(g) N2(g)
</PRE></CODE>
</TD>
<TD width=35%>
<CODE><PRE>
  pressure "total mol" volume g_CO2(g) g_N2(g)
</PRE></CODE>
</TD>
</TR>

<TR VALIGN="top">
<TD width=65%>
<CODE><PRE>
  SELECTED_OUTPUT
    -reset false
    -kinetic_reactants CH2O Pyrite
</PRE></CODE>
</TD>
<TD width=35%>
<CODE><PRE>
  k_CH2O dk_CH2O k_Pyrite dk_Pyrite
</PRE></CODE>
</TD>
</TR>

<TR VALIGN="top">
<TD width=65%>
<CODE><PRE>
  SELECTED_OUTPUT
    -reset false
    -solid_solutions CaSO4 SrSO4
</PRE></CODE>
</TD>
<TD width=35%>
<CODE><PRE>
  s_CaSO4 s_SrSO4
</PRE></CODE>
</TD>
</TR>

</table>
 *  @endhtmlonly
 *  @par Fortran90 Interface:
 *  ROW is 1-based for the Fortran interface except that the column headings are stored in ROW=0.
 *  COL is 1-based for the Fortran interface.
 *  @htmlonly
 *  <CODE>
 *  <PRE>
 *  FUNCTION GetSelectedOutputValue(ID,ROW,COL,VTYPE,DVALUE,SVALUE,SLENGTH)
 *    INTEGER(KIND=4),   INTENT(IN)            :: ID
 *    INTEGER(KIND=4),   INTENT(IN)            :: ROW
 *    INTEGER(KIND=4),   INTENT(IN)            :: COL
 *    INTEGER(KIND=4),   INTENT(OUT)           :: VTYPE
 *    REAL(KIND=8),      INTENT(OUT)           :: DVALUE
 *    CHARACTER(LEN=*),  INTENT(OUT)           :: SVALUE
 *    INTEGER(KIND=4),   INTENT(OUT), OPTIONAL :: SLENGTH
 *    INTEGER(KIND=4)                          :: GetSelectedOutputValue
 *  END FUNCTION GetSelectedOutputValue
 *  </PRE>
 *  </CODE>
 *  @endhtmlonly
 *  @param ID                The instance id returned from @ref CreateIPhreeqc.
 *  @param ROW               The row index.
 *  @param COL               The column index.
 *  @param VTYPE             Returns the variable type.  See @ref VAR_TYPE.
 *  @param DVALUE            Returns the numeric value when (VTYPE=@ref TT_DOUBLE) or (VTYPE=@ref TT_LONG).
 *  @param SVALUE            Returns the string variable when (VTYPE=@ref TT_STRING).  When (VTYPE=@ref TT_DOUBLE) or (VTYPE=@ref TT_LONG) this variable is filled with a string equivalent of DVALUE.
 *  @param SLENGTH           Optional, if the length of SVALUE isn't sufficient to hold the entire string value, returns the required length, otherwise returns 0 (zero).
 *  @anchor GetSelectedOutputValue_c
 *  @par  C Example:
 *  @include GetSelectedOutputValue.c
 *
 *  @anchor F90GetSelectedOutputValue_f90
 *  @par  Fortran90 Example:
 *  @include F90GetSelectedOutputValue.f90
 */
	IPQ_DLL_EXPORT IPQ_RESULT  GetSelectedOutputValue(int id, int row, int col, VAR* pVAR);


/**
 *  Returns the associated data with the specified row and column.  The current <b>SELECTED_OUTPUT</b> block is set using the @ref SetCurrentSelectedOutputUserNumber method.
 *  @param id                The instance id returned from @ref CreateIPhreeqc.
 *  @param row               The row index.
 *  @param col               The column index.
 *  @param vtype             Receives the variable type.  See @ref VAR_TYPE.
 *  @param dvalue            Receives the numeric value when (VTYPE=@ref TT_DOUBLE) or (VTYPE=@ref TT_LONG).
 *  @param svalue            Receives the string variable when (VTYPE=@ref TT_STRING).  When (VTYPE=@ref TT_DOUBLE) or (VTYPE=@ref TT_LONG) this variable is filled with a string equivalent of DVALUE.
 *  @param svalue_length     The length of the svalue buffer.
 *  @retval IPQ_OK           Success.
 *  @retval IPQ_INVALIDROW   The given row is out of range.
 *  @retval IPQ_INVALIDCOL   The given column is out of range.
 *  @retval IPQ_OUTOFMEMORY  Memory could not be allocated.
 *  @retval IPQ_BADINSTANCE  The given id is invalid.
 *  @see                     GetCurrentSelectedOutputUserNumber, GetSelectedOutputFileOn, GetSelectedOutputColumnCount, GetSelectedOutputRowCount, GetSelectedOutputValue, SetCurrentSelectedOutputUserNumber, SetSelectedOutputFileOn
 *  @remarks
 *  Row 0 contains the column headings to the selected_ouput.
 *  @par Examples:
 *  The headings will include a suffix and/or prefix in order to differentiate the
 *  columns.
 *  @htmlonly
<p>
<table border=1>

<TR VALIGN="top">
<TH width=65%>
Input
</TH>
<TH width=35%>
Headings
</TH>
</TR>

<TR VALIGN="top">
<TD width=65%>
<CODE><PRE>
  SELECTED_OUTPUT
    -reset false
    -totals Ca Na
</PRE></CODE>
</TD>
<TD width=35%>
<CODE><PRE>
  Ca(mol/kgw)  Na(mol/kgw)
</PRE></CODE>
</TD>
</TR>

<TR VALIGN="top">
<TD width=65%>
<CODE><PRE>
  SELECTED_OUTPUT
    -reset false
    -molalities Fe+2 Hfo_sOZn+
</PRE></CODE>
</TD>
<TD width=35%>
<CODE><PRE>
  m_Fe+2(mol/kgw)  m_Hfo_sOZn+(mol/kgw)
</PRE></CODE>
</TD>
</TR>

<TR VALIGN="top">
<TD width=65%>
<CODE><PRE>
  SELECTED_OUTPUT
    -reset false
    -activities H+ Ca+2
</PRE></CODE>
</TD>
<TD width=35%>
<CODE><PRE>
  la_H+  la_Ca+2
</PRE></CODE>
</TD>
</TR>

<TR VALIGN="top">
<TD width=65%>
<CODE><PRE>
  SELECTED_OUTPUT
    -reset false
    -equilibrium_phases Calcite Dolomite
</PRE></CODE>
</TD>
<TD width=35%>
<CODE><PRE>
  Calcite  d_Calcite  Dolomite  d_Dolomite
</PRE></CODE>
</TD>
</TR>

<TR VALIGN="top">
<TD width=65%>
<CODE><PRE>
  SELECTED_OUTPUT
    -reset false
    -saturation_indices CO2(g) Siderite
</PRE></CODE>
</TD>
<TD width=35%>
<CODE><PRE>
  si_CO2(g)  si_Siderite
</PRE></CODE>
</TD>
</TR>

<TR VALIGN="top">
<TD width=65%>
<CODE><PRE>
  SELECTED_OUTPUT
    -reset false
    -gases CO2(g) N2(g)
</PRE></CODE>
</TD>
<TD width=35%>
<CODE><PRE>
  pressure "total mol" volume g_CO2(g) g_N2(g)
</PRE></CODE>
</TD>
</TR>

<TR VALIGN="top">
<TD width=65%>
<CODE><PRE>
  SELECTED_OUTPUT
    -reset false
    -kinetic_reactants CH2O Pyrite
</PRE></CODE>
</TD>
<TD width=35%>
<CODE><PRE>
  k_CH2O dk_CH2O k_Pyrite dk_Pyrite
</PRE></CODE>
</TD>
</TR>

<TR VALIGN="top">
<TD width=65%>
<CODE><PRE>
  SELECTED_OUTPUT
    -reset false
    -solid_solutions CaSO4 SrSO4
</PRE></CODE>
</TD>
<TD width=35%>
<CODE><PRE>
  s_CaSO4 s_SrSO4
</PRE></CODE>
</TD>
</TR>

</table>
 *  @endhtmlonly
 *  @anchor GetSelectedOutputValue2_c
 *  @par  C Example:
 *  @include GetSelectedOutputValue2.c
 */
	IPQ_DLL_EXPORT IPQ_RESULT  GetSelectedOutputValue2(int id, int row, int col, int *vtype, double* dvalue, char* svalue, unsigned int svalue_length);


/**
 *  Retrieves the string buffer containing the version in the form of X.X.X-XXXX.
 *  @return              A null terminated string containing the IPhreeqc version number.
 *  @par Fortran90 Interface:
 *  @htmlonly
 *  <CODE>
 *  <PRE>
 *  SUBROUTINE GetVersionString(VERSION)
 *    CHARACTER(LEN=*), INTENT(OUT) :: VERSION
 *  END SUBROUTINE GetVersionString
 *  </PRE>
 *  </CODE>
 *  @endhtmlonly
 *
 *  @par C Example:
 *  @include GetVersionString.c
 *
 *  @par Fortran90 Example:
 *  @include F90GetVersionString.f90
 */
	IPQ_DLL_EXPORT const char* GetVersionString(void);


/**
 *  Retrieves the warning messages from the last call to (@ref RunAccumulated, @ref RunFile, @ref RunString, @ref LoadDatabase, or @ref LoadDatabaseString).
 *  @param id            The instance id returned from @ref CreateIPhreeqc.
 *  @return              A null terminated string containing warning messages.
 *  @see                 GetWarningStringLine, GetWarningStringLineCount, OutputWarningString
 *  @par Fortran90 Interface:
 *  Not implemented. (see @ref GetWarningStringLineCount, @ref GetWarningStringLine, @ref OutputWarningString)
 */
	IPQ_DLL_EXPORT const char* GetWarningString(int id);


/**
 *  Retrieves the given warning line.
 *  @param id            The instance id returned from @ref CreateIPhreeqc.
 *  @param n             The zero-based index of the line to retrieve.
 *  @return              A null terminated string containing the given warning line message.
 *  @see                 GetWarningString, GetWarningStringLineCount, OutputWarningString
 *  @par Fortran90 Interface:
 *  (Note: N is one-based for the Fortran interface.)
 *  @htmlonly
 *  <CODE>
 *  <PRE>
 *  SUBROUTINE GetWarningStringLine(ID,N,LINE)
 *    INTEGER(KIND=4),  INTENT(IN)  :: ID
 *    LOGICAL(KIND=4),  INTENT(IN)  :: N
 *    CHARACTER(LEN=*), INTENT(OUT) :: LINE
 *  END SUBROUTINE GetWarningStringLine
 *  </PRE>
 *  </CODE>
 *  @endhtmlonly
 */
	IPQ_DLL_EXPORT const char* GetWarningStringLine(int id, int n);


/**
 *  Retrieves the number of lines in the current warning string buffer.
 *  @param id            The instance id returned from @ref CreateIPhreeqc.
 *  @return              The number of lines.
 *  @see                 GetWarningString, GetWarningStringLine, OutputWarningString
 *  @par Fortran90 Interface:
 *  @htmlonly
 *  <CODE>
 *  <PRE>
 *  FUNCTION GetWarningStringLineCount(ID)
 *    INTEGER(KIND=4),  INTENT(IN)  :: ID
 *    INTEGER(KIND=4)               :: GetWarningStringLineCount
 *  END FUNCTION GetWarningStringLineCount
 *  </PRE>
 *  </CODE>
 *  @endhtmlonly
 */
	IPQ_DLL_EXPORT int         GetWarningStringLineCount(int id);


/**
 *  Load the specified database file into phreeqc.
 *  @param id            The instance id returned from @ref CreateIPhreeqc.
 *  @param filename      The name of the phreeqc database to load.
 *                       The full path (or relative path with respect to the working directory)
 *                       must be given if the file is not in the current working directory.
 *  @return              The number of errors encountered.
 *  @see                 LoadDatabaseString
 *  @remarks
 *  All previous definitions are cleared.
 *  @par Fortran90 Interface:
 *  @htmlonly
 *  <CODE>
 *  <PRE>
 *  FUNCTION LoadDatabase(ID,FILENAME)
 *    INTEGER(KIND=4),   INTENT(IN)  :: ID
 *    CHARACTER(LEN=*),  INTENT(IN)  :: FILENAME
 *    INTEGER(KIND=4)                :: LoadDatabase
 *  END FUNCTION LoadDatabase
 *  </PRE>
 *  </CODE>
 *  @endhtmlonly
 *
 *  @par C Example:
 *  see @ref CreateIPhreeqc_c "CreateIPhreeqc"
 *
 *  @par Fortran90 Example:
 *  see @ref CreateIPhreeqc_f90 "CreateIPhreeqc"
 */
	IPQ_DLL_EXPORT int         LoadDatabase(int id, const char* filename);


/**
 *  Load the specified string as a database into phreeqc.
 *  @param id            The instance id returned from @ref CreateIPhreeqc.
 *  @param input         String containing data to be used as the phreeqc database.
 *  @return              The number of errors encountered.
 *  @see                 LoadDatabase
 *  @remarks
 *  All previous definitions are cleared.
 *  @par Fortran90 Interface:
 *  @htmlonly
 *  <CODE>
 *  <PRE>
 *  FUNCTION LoadDatabaseString(ID,INPUT)
 *    INTEGER(KIND=4),   INTENT(IN)  :: ID
 *    CHARACTER(LEN=*),  INTENT(IN)  :: INPUT
 *    INTEGER(KIND=4)                :: LoadDatabaseString
 *  END FUNCTION LoadDatabaseString
 *  </PRE>
 *  </CODE>
 *  @endhtmlonly
 */
	IPQ_DLL_EXPORT int         LoadDatabaseString(int id, const char* input);


/**
 *  Output the accumulated input buffer to stdout.  This input buffer can be run with a call to @ref RunAccumulated.
 *  @param id            The instance id returned from @ref CreateIPhreeqc.
 *  @see                 AccumulateLine, ClearAccumulatedLines, RunAccumulated
 *  @par Fortran90 Interface:
 *  @htmlonly
 *  <CODE>
 *  <PRE>
 *  SUBROUTINE OutputAccumulatedLines(ID)
 *    INTEGER(KIND=4),  INTENT(IN)  :: ID
 *  END SUBROUTINE OutputAccumulatedLines
 *  </PRE>
 *  </CODE>
 *  @endhtmlonly
 *
 *  @par Fortran90 Example:
 *  see @ref GetDumpStringLine_f90 "GetDumpStringLine"
 */
	IPQ_DLL_EXPORT void        OutputAccumulatedLines(int id);


/**
 *  Output the error messages normally stored in the <B><I>phreeqc.id.err</I></B> file to stdout.
 *  @param id            The instance id returned from @ref CreateIPhreeqc.
 *  @see                 GetErrorFileOn, GetErrorStringLine, GetErrorStringLineCount, SetErrorFileOn
 *  @par Fortran90 Interface:
 *  @htmlonly
 *  <CODE>
 *  <PRE>
 *  SUBROUTINE OutputErrorString(ID)
 *    INTEGER(KIND=4),  INTENT(IN)  :: ID
 *  END SUBROUTINE OutputErrorString
 *  </PRE>
 *  </CODE>
 *  @endhtmlonly
 *
 *  @par C Example:
 *  see @ref GetComponent_c "GetComponent"
 *
 *  @par Fortran90 Example:
 *  see @ref GetDumpStringLine_f90 "GetDumpStringLine"
 */
	IPQ_DLL_EXPORT void        OutputErrorString(int id);


/**
 *  Output the warning messages to stdout.
 *  @param id            The instance id returned from @ref CreateIPhreeqc.
 *  @see                 GetWarningString, GetWarningStringLine, GetWarningStringLineCount
 *  @par Fortran90 Interface:
 *  @htmlonly
 *  <CODE>
 *  <PRE>
 *  SUBROUTINE OutputWarningString(ID)
 *    INTEGER(KIND=4),  INTENT(IN)  :: ID
 *  END SUBROUTINE OutputWarningString
 *  </PRE>
 *  </CODE>
 *  @endhtmlonly
 */
	IPQ_DLL_EXPORT void        OutputWarningString(int id);


/**
 *  Runs the input buffer as defined by calls to @ref AccumulateLine.
 *  @param id            The instance id returned from @ref CreateIPhreeqc.
 *  @return              The number of errors encountered.
 *  @see                 AccumulateLine, ClearAccumulatedLines, OutputAccumulatedLines, RunFile, RunString
 *  @remarks
 *  The accumulated input is cleared at the next call to @ref AccumulateLine.
 *  @pre @ref LoadDatabase/@ref LoadDatabaseString must have been called and returned 0 (zero) errors.
 *  @par Fortran90 Interface:
 *  @htmlonly
 *  <CODE>
 *  <PRE>
 *  FUNCTION RunAccumulated(ID)
 *    INTEGER(KIND=4),  INTENT(IN)  :: ID
 *    INTEGER(KIND=4)               :: RunAccumulated
 *  END FUNCTION RunAccumulated
 *  </PRE>
 *  </CODE>
 *  @endhtmlonly
 *
 *  @par Fortran90 Example:
 *  see @ref GetDumpStringLine_f90 "GetDumpStringLine"
 */
	IPQ_DLL_EXPORT int         RunAccumulated(int id);


/**
 *  Runs the specified phreeqc input file.
 *  @param id            The instance id returned from @ref CreateIPhreeqc.
 *  @param filename      The name of the phreeqc input file to run.
 *  @return              The number of errors encountered during the run.
 *  @see                 RunAccumulated, RunString
 *  @pre                 (@ref LoadDatabase, @ref LoadDatabaseString) must have been called and returned 0 (zero) errors.
 *  @par Fortran90 Interface:
 *  @htmlonly
 *  <CODE>
 *  <PRE>
 *  FUNCTION RunFile(ID,FNAME)
 *    INTEGER(KIND=4),   INTENT(IN)  :: ID
 *    CHARACTER(LEN=*),  INTENT(IN)  :: FNAME
 *    INTEGER(KIND=4)                :: RunFile
 *  END FUNCTION RunFile
 *  </PRE>
 *  </CODE>
 *  @endhtmlonly
 *
 *  @par C Example:
 *  see @ref CreateIPhreeqc_c "CreateIPhreeqc"
 *
 *  @par Fortran90 Example:
 *  see @ref CreateIPhreeqc_f90 "CreateIPhreeqc"
 */
	IPQ_DLL_EXPORT int         RunFile(int id, const char* filename);


/**
 *  Runs the specified string as input to phreeqc.
 *  @param id            The instance id returned from @ref CreateIPhreeqc.
 *  @param input         String containing phreeqc input.
 *  @return              The number of errors encountered during the run.
 *  @see                 RunAccumulated, RunFile
 *  @pre                 (@ref LoadDatabase, @ref LoadDatabaseString) must have been called and returned 0 (zero) errors.
 *  @par Fortran90 Interface:
 *  @htmlonly
 *  <CODE>
 *  <PRE>
 *  FUNCTION RunString(ID,INPUT)
 *    INTEGER(KIND=4),  INTENT(IN)   :: ID
 *    CHARACTER(LEN=*),  INTENT(IN)  :: INPUT
 *    INTEGER(KIND=4)                :: RunString
 *  END FUNCTION RunString
 *  </PRE>
 *  </CODE>
 *  @endhtmlonly
 *
 *  @par C Example:
 *  see @ref GetDumpString_c "GetDumpString"
 *
 */
	IPQ_DLL_EXPORT int         RunString(int id, const char* input);

/**
 *  Sets a C callback function for Basic programs. The syntax for the Basic command is
 *  10 result = CALLBACK(x1, x2, string$)
 *  The syntax for the C function is
 *  double my_callback(double x1, double x2, const char * string)
 *  @param id               The instance id returned from @ref CreateIPhreeqc.
 *  @param fcn              The name of a user-defined function.
 *  @param cookie1          A user defined value to be passed to the callback function.
 *  @retval IPQ_OK          Success.
 *  @retval IPQ_BADINSTANCE The given id is invalid.
 *  @par Fortran90 Interface:
 *  see @ref SetBasicFortranCallback
 *  @par C Example:
 *  @include SetBasicCallback.c
 */
	IPQ_DLL_EXPORT IPQ_RESULT  SetBasicCallback(int id, double (*fcn)(double x1, double x2, const char *str, void *cookie), void *cookie1);

/**
 *  Sets a Fortran callback function for Basic programs. The syntax for the Basic command is
 *  10 result = CALLBACK(x1, x2, string$)
 *  The syntax for the Fortran function is
 *  REAL(KIND=C_DOUBLE) my_callback(x1, x2, string), where x1 and x2 are REAL(KIND=C_DOUBLE) and string is a CHARACTER(KIND=C_CHAR).
 *  @param id               The instance id returned from @ref CreateIPhreeqc.
 *  @param fcn              The name of a REAL(KIND=C_DOUBLE) Fortran function with three arguments (two REAL(KIND=C_DOUBLE), and one CHARACTER(KIND=C_CHAR)).
 *  @retval IPQ_OK          Success.
 *  @retval IPQ_BADINSTANCE The given id is invalid.
 *  @par Fortran90 Interface:
 *  @htmlonly
 *  <CODE>
 *  <PRE>
 *  !
 *  ! if using  include files (IPhreeqc.f.inc or IPhreeqc.f90.inc)
 *  !
 *  #ifdef IPHREEQC_NO_FORTRAN_MODULE
 *  FUNCTION SetBasicFortranCallback(ID,FCN)
 *    INTEGER(KIND=4),  INTENT(IN)  :: ID
 *    INTERFACE
 *      real(kind=8) FUNCTION FCN(x1, x2, str)
 *        real(kind=8), INTENT(in) :: x1
 *        real(kind=8), INTENT(in) :: x2
 *        CHARACTER(*), INTENT(in)     :: str
 *      END FUNCTION
 *    END INTERFACE
 *    INTEGER(KIND=4)               :: SetBasicFortranCallback
 *  END FUNCTION SetBasicFortranCallback
 *  #else
 *  !
 *  ! if using the fortran module (USE IPhreeqc)
 *  ! must also add IPhreeqc_interface.F90 to your project
 *  !
 *  FUNCTION SetBasicFortranCallback(ID,FCN)
 *    INTEGER(KIND=4),  INTENT(IN)  :: ID
 *    INTERFACE
 *      REAL(KIND=C_DOUBLE) FUNCTION fcn(x1, x2, str, l) BIND(C)
 *        USE ISO_C_BINDING
 *        IMPLICIT none
 *        REAL(KIND=C_DOUBLE),    INTENT(in)        :: x1, x2
 *        CHARACTER(KIND=C_CHAR), INTENT(in)        :: str(*)
 *        INTEGER(KIND=C_INT),    INTENT(in), VALUE :: l
 *      END FUNCTION fcn
 *    END INTERFACE
 *  END FUNCTION SetBasicFortranCallback
 *  #endif
 *  </PRE>
 *  </CODE>
 *  @endhtmlonly
 *  @par  Fortran90 Example:
 *  @include F90SetBasicFortranCallback.f90
 *  @par File ic :
 *  @include ic
 */
#ifdef IPHREEQC_NO_FORTRAN_MODULE
	IPQ_DLL_EXPORT IPQ_RESULT  SetBasicFortranCallback(int id, double (*fcn)(double *x1, double *x2, const char *str, size_t l));
#else
	IPQ_DLL_EXPORT IPQ_RESULT  SetBasicFortranCallback(int id, double (*fcn)(double *x1, double *x2, const char *str, int l));
#endif


/**
 *  Sets the current <B>SELECTED_OUTPUT</B> user number for use in subsequent calls to (@ref GetSelectedOutputColumnCount,
 *  @ref GetSelectedOutputFileName, @ref GetSelectedOutputRowCount, @ref GetSelectedOutputString, @ref GetSelectedOutputStringLine,
 *  @ref GetSelectedOutputStringLineCount, @ref GetSelectedOutputValue, @ref GetSelectedOutputValue2) routines.
 *  The initial setting after calling @ref CreateIPhreeqc is 1.
 *  @param id                   The instance id returned from @ref CreateIPhreeqc.
 *  @param n                    The user number specified in the <B>SELECTED_OUTPUT</B> block.
 *  @retval IPQ_OK              Success.
 *  @retval IPQ_BADINSTANCE     The given id is invalid.
 *  @retval IPQ_INVALIDARG      The given user number is invalid.
 *  @see                        GetSelectedOutputColumnCount, GetSelectedOutputFileName, GetSelectedOutputRowCount, GetSelectedOutputString, GetSelectedOutputStringLine, GetSelectedOutputStringLineCount, GetSelectedOutputValue
 *  @par Fortran90 Interface:
 *  @htmlonly
 *  <CODE>
 *  <PRE>
 *  FUNCTION SetCurrentSelectedOutputUserNumber(ID,N)
 *    INTEGER(KIND=4),  INTENT(IN)  :: ID
 *    INTEGER(KIND=4),  INTENT(IN)  :: N
 *    INTEGER(KIND=4)               :: SetCurrentSelectedOutputUserNumber
 *  END FUNCTION SetCurrentSelectedOutputUserNumber
 *  </PRE>
 *  </CODE>
 *  @endhtmlonly
 *
 *  @anchor SetCurrentSelectedOutputUserNumber_c
 *  @par C Example:
 *  @include SetCurrentSelectedOutputUserNumber.c
 *  @par File multi_punch :
 *  @include multi_punch
 */
	IPQ_DLL_EXPORT IPQ_RESULT  SetCurrentSelectedOutputUserNumber(int id, int n);

/**
 *  Sets the name of the dump file.  This file name is used if not specified within <B>DUMP</B> input.
 *  The default value is <B><I>dump.id.out</I></B>.
 *  @param id               The instance id returned from @ref CreateIPhreeqc.
 *  @param filename         The name of the file to write <B>DUMP</B> output to.
 *  @retval IPQ_OK          Success.
 *  @retval IPQ_BADINSTANCE The given id is invalid.
 *  @see                    GetDumpFileName, GetDumpFileOn, GetDumpString, GetDumpStringOn, GetDumpStringLine, GetDumpStringLineCount, SetDumpFileOn, SetDumpStringOn
 *  @par Fortran90 Interface:
 *  @htmlonly
 *  <CODE>
 *  <PRE>
 *  FUNCTION SetDumpFileName(ID,FILENAME)
 *    INTEGER(KIND=4),   INTENT(IN)   :: ID
 *    CHARACTER(LEN=*),  INTENT(OUT)  :: FILENAME
 *    INTEGER(KIND=4)                 :: SetDumpFileName
 *  END FUNCTION SetDumpFileName
 *  </PRE>
 *  </CODE>
 *  @endhtmlonly
 */
	IPQ_DLL_EXPORT IPQ_RESULT  SetDumpFileName(int id, const char* filename);


/**
 *  Sets the dump file switch on or off.  This switch controls whether or not phreeqc writes to the dump file.
 *  The initial setting after calling @ref CreateIPhreeqc is off.
 *  @param id                   The instance id returned from @ref CreateIPhreeqc.
 *  @param dump_on              If non-zero, turns on output to the <B>DUMP</B> (<B><I>dump.id.out</I></B> if unspecified) file;
 *                              if zero, turns off output to the <B>DUMP</B> file.
 *  @retval IPQ_OK              Success.
 *  @retval IPQ_BADINSTANCE     The given id is invalid.
 *  @see                        GetDumpFileOn, GetDumpString, GetDumpStringLine, GetDumpStringOn, GetDumpStringLineCount, SetDumpStringOn
 *  @par Fortran90 Interface:
 *  @htmlonly
 *  <CODE>
 *  <PRE>
 *  FUNCTION SetDumpFileOn(ID,DUMP_ON)
 *    INTEGER(KIND=4),  INTENT(IN)  :: ID
 *    LOGICAL(KIND=4),  INTENT(IN)  :: DUMP_ON
 *    INTEGER(KIND=4)               :: SetDumpFileOn
 *  END FUNCTION SetDumpFileOn
 *  </PRE>
 *  </CODE>
 *  @endhtmlonly
 */
	IPQ_DLL_EXPORT IPQ_RESULT  SetDumpFileOn(int id, int dump_on);


/**
 *  Sets the dump string switch on or off.  This switch controls whether or not the data normally sent
 *  to the dump file are stored in a buffer for retrieval.  The initial setting after calling
 *  @ref CreateIPhreeqc is off.
 *  @param id                   The instance id returned from @ref CreateIPhreeqc.
 *  @param dump_string_on       If non-zero, captures the output defined by the <B>DUMP</B> keyword into a string buffer;
 *                              if zero, output defined by the <B>DUMP</B> keyword is not captured to a string buffer.
 *  @retval IPQ_OK              Success.
 *  @retval IPQ_BADINSTANCE     The given id is invalid.
 *  @see                        GetDumpFileOn, GetDumpStringOn, GetDumpString, GetDumpStringLine, GetDumpStringLineCount, SetDumpFileOn
 *  @par Fortran90 Interface:
 *  @htmlonly
 *  <CODE>
 *  <PRE>
 *  FUNCTION SetDumpStringOn(ID,DUMP_STRING_ON)
 *    INTEGER(KIND=4),  INTENT(IN)  :: ID
 *    LOGICAL(KIND=4),  INTENT(IN)  :: DUMP_STRING_ON
 *    INTEGER(KIND=4)               :: SetDumpStringOn
 *  END FUNCTION SetDumpStringOn
 *  </PRE>
 *  </CODE>
 *  @endhtmlonly
 *
 *  @par C Example:
 *  see @ref GetDumpString_c "GetDumpString"
 *
 *  @par Fortran90 Example:
 *  see @ref GetDumpStringLine_f90 "GetDumpStringLine"
 */
	IPQ_DLL_EXPORT IPQ_RESULT  SetDumpStringOn(int id, int dump_string_on);

/**
 *  Sets the name of the error file.  The default value is <B><I>phreeqc.id.err</I></B>.
 *  @param id               The instance id returned from @ref CreateIPhreeqc.
 *  @param filename         The name of the error file.
 *  @retval IPQ_OK          Success.
 *  @retval IPQ_BADINSTANCE The given id is invalid.
 *  @see                    GetErrorFileName, GetErrorFileOn, GetErrorString, GetErrorStringOn, GetErrorStringLine, GetErrorStringLineCount, SetErrorFileOn, SetErrorStringOn
 *  @par Fortran90 Interface:
 *  @htmlonly
 *  <CODE>
 *  <PRE>
 *  FUNCTION SetErrorFileName(ID,FILENAME)
 *    INTEGER(KIND=4),   INTENT(IN)   :: ID
 *    CHARACTER(LEN=*),  INTENT(OUT)  :: FILENAME
 *    INTEGER(KIND=4)                 :: SetErrorFileName
 *  END FUNCTION SetErrorFileName
 *  </PRE>
 *  </CODE>
 *  @endhtmlonly
 */
	IPQ_DLL_EXPORT IPQ_RESULT  SetErrorFileName(int id, const char* filename);

/**
 *  Sets the error file switch on or off.  This switch controls whether or not
 *  error messages are written to the <B><I>phreeqc.id.err</I></B> file.  The initial setting after calling
 *  @ref CreateIPhreeqc is off.
 *  @param id                   The instance id returned from @ref CreateIPhreeqc.
 *  @param error_on             If non-zero, writes errors to the error file; if zero, no errors are written to the error file.
 *  @retval IPQ_OK              Success.
 *  @retval IPQ_BADINSTANCE     The given id is invalid.
 *  @see                        GetErrorFileOn, GetErrorStringLine, GetErrorStringLineCount, OutputErrorString
 *  @par Fortran90 Interface:
 *  @htmlonly
 *  <CODE>
 *  <PRE>
 *  FUNCTION SetErrorFileOn(ID,ERR_ON)
 *    INTEGER(KIND=4),  INTENT(IN)  :: ID
 *    LOGICAL(KIND=4),  INTENT(IN)  :: ERR_ON
 *    INTEGER(KIND=4)               :: SetErrorFileOn
 *  END FUNCTION SetErrorFileOn
 *  </PRE>
 *  </CODE>
 *  @endhtmlonly
 */
	IPQ_DLL_EXPORT IPQ_RESULT  SetErrorFileOn(int id, int error_on);

/**
 *  Sets the error switch on or off.  This switch controls whether or not
 *  error messages are generated and displayed.  The initial setting after calling
 *  @ref CreateIPhreeqc is on.
 *  @param id                   The instance id returned from @ref CreateIPhreeqc.
 *  @param error_on             If non-zero, writes errors to the error file and error string; if zero, no errors are written to the error file or stored in the error string.
 *  @retval IPQ_OK              Success.
 *  @retval IPQ_BADINSTANCE     The given id is invalid.
 *  @see                        GetErrorOn, GetErrorStringLine, GetErrorStringLineCount, OutputErrorString
 *  @par Fortran90 Interface:
 *  @htmlonly
 *  <CODE>
 *  <PRE>
 *  FUNCTION SetErrorOn(ID,ERR_ON)
 *    INTEGER(KIND=4),  INTENT(IN)  :: ID
 *    LOGICAL(KIND=4),  INTENT(IN)  :: ERR_ON
 *    INTEGER(KIND=4)               :: SetErrorOn
 *  END FUNCTION SetErrorOn
 *  </PRE>
 *  </CODE>
 *  @endhtmlonly
 */
	IPQ_DLL_EXPORT IPQ_RESULT  SetErrorOn(int id, int error_on);


/**
 *  Sets the error string switch on or off.  This switch controls whether or not the data normally sent
 *  to the error file are stored in a buffer for retrieval.  The initial setting after calling
 *  @ref CreateIPhreeqc is on.
 *  @param id                   The instance id returned from @ref CreateIPhreeqc.
 *  @param error_string_on      If non-zero, captures the error output into a string buffer;
 *                              if zero, error output is not captured to a string buffer.
 *  @retval IPQ_OK              Success.
 *  @retval IPQ_BADINSTANCE     The given id is invalid.
 *  @see                        GetErrorFileOn, GetErrorStringOn, GetErrorString, GetErrorStringLine, GetErrorStringLineCount, SetErrorFileOn
 *  @par Fortran90 Interface:
 *  @htmlonly
 *  <CODE>
 *  <PRE>
 *  FUNCTION SetErrorStringOn(ID,ERR_STRING_ON)
 *    INTEGER(KIND=4),  INTENT(IN)  :: ID
 *    LOGICAL(KIND=4),  INTENT(IN)  :: ERR_STRING_ON
 *    INTEGER(KIND=4)               :: SetErrorStringOn
 *  END FUNCTION SetErrorStringOn
 *  </PRE>
 *  </CODE>
 *  @endhtmlonly
 *
 *  @par C Example:
 *  see @ref GetErrorString_c "GetErrorString"
 *
 *  @par Fortran90 Example:
 *  see @ref GetErrorStringLine_f90 "GetErrorStringLine"
 */
	IPQ_DLL_EXPORT IPQ_RESULT  SetErrorStringOn(int id, int error_string_on);

/**
 *  Sets the name of the log file.  The default value is <B><I>phreeqc.id.log</I></B>.
 *  @param id               The instance id returned from @ref CreateIPhreeqc.
 *  @param filename         The name of the log file.
 *  @retval IPQ_OK          Success.
 *  @retval IPQ_BADINSTANCE The given id is invalid.
 *  @see                    GetLogFileName, GetLogFileOn, GetLogString, GetLogStringOn, GetLogStringLine, GetLogStringLineCount, SetLogFileOn, SetLogStringOn
 *  @par Fortran90 Interface:
 *  @htmlonly
 *  <CODE>
 *  <PRE>
 *  FUNCTION SetLogFileName(ID,FILENAME)
 *    INTEGER(KIND=4),   INTENT(IN)   :: ID
 *    CHARACTER(LEN=*),  INTENT(OUT)  :: FILENAME
 *    INTEGER(KIND=4)                 :: SetLogFileName
 *  END FUNCTION SetLogFileName
 *  </PRE>
 *  </CODE>
 *  @endhtmlonly
 */
	IPQ_DLL_EXPORT IPQ_RESULT  SetLogFileName(int id, const char* filename);

/**
 *  Sets the log file switch on or off.  This switch controls whether or not phreeqc
 *  writes log messages to the <B><I>phreeqc.id.log</I></B> file.  The initial setting after calling
 *  @ref CreateIPhreeqc is off.
 *  @param id            The instance id returned from @ref CreateIPhreeqc.
 *  @param log_on        If non-zero, log messages are written to the log file; if zero, no log messages are written to the log file.
 *  @retval IPQ_OK       Success.
 *  @retval              IPQ_BADINSTANCE The given id is invalid.
 *  @remarks
 *      Logging must be enabled through the use of the KNOBS -logfile option in order to receive any log messages.
 *  @see                 GetLogFileOn
 *  @par Fortran90 Interface:
 *  @htmlonly
 *  <CODE>
 *  <PRE>
 *  FUNCTION SetLogFileOn(ID,LOG_ON)
 *    INTEGER(KIND=4),  INTENT(IN)  :: ID
 *    LOGICAL(KIND=4),  INTENT(IN)  :: LOG_ON
 *    INTEGER(KIND=4)               :: SetLogFileOn
 *  END FUNCTION SetLogFileOn
 *  </PRE>
 *  </CODE>
 *  @endhtmlonly
 */
	IPQ_DLL_EXPORT IPQ_RESULT  SetLogFileOn(int id, int log_on);

/**
 *  Sets the log string switch on or off.  This switch controls whether or not the data normally sent
 *  to the log file are stored in a buffer for retrieval.  The initial setting after calling
 *  @ref CreateIPhreeqc is off.
 *  @param id                   The instance id returned from @ref CreateIPhreeqc.
 *  @param log_string_on        If non-zero, captures the log output into a string buffer;
 *                              if zero, log output is not captured to a string buffer.
 *  @retval IPQ_OK              Success.
 *  @retval IPQ_BADINSTANCE     The given id is invalid.
 *  @see                        GetLogFileOn, GetLogStringOn, GetLogString, GetLogStringLine, GetLogStringLineCount, SetLogFileOn
 *  @par Fortran90 Interface:
 *  @htmlonly
 *  <CODE>
 *  <PRE>
 *  FUNCTION SetLogStringOn(ID,LOG_STRING_ON)
 *    INTEGER(KIND=4),  INTENT(IN)  :: ID
 *    LOGICAL(KIND=4),  INTENT(IN)  :: LOG_STRING_ON
 *    INTEGER(KIND=4)               :: SetLogStringOn
 *  END FUNCTION SetLogStringOn
 *  </PRE>
 *  </CODE>
 *  @endhtmlonly
 *
 *  @par C Example:
 *  see @ref GetLogString_c "GetLogString"
 *
 *  @par Fortran90 Example:
 *  see @ref GetLogStringLine_f90 "GetLogStringLine"
 */
	IPQ_DLL_EXPORT IPQ_RESULT  SetLogStringOn(int id, int log_string_on);


/**
 *  Sets the name of the output file.  This file name is used if not specified within <B>DUMP</B> input.
 *  The default value is <B><I>phreeqc.id.out</I></B>.
 *  @param id               The instance id returned from @ref CreateIPhreeqc.
 *  @param filename         The name of the phreeqc output file.
 *  @retval IPQ_OK          Success.
 *  @retval IPQ_BADINSTANCE The given id is invalid.
 *  @see                    GetOutputFileName, GetOutputFileOn, GetOutputString, GetOutputStringOn, GetOutputStringLine, GetOutputStringLineCount, SetOutputFileOn, SetOutputStringOn
 *  @par Fortran90 Interface:
 *  @htmlonly
 *  <CODE>
 *  <PRE>
 *  FUNCTION SetOutputFileName(ID,FILENAME)
 *    INTEGER(KIND=4),   INTENT(IN)   :: ID
 *    CHARACTER(LEN=*),  INTENT(OUT)  :: FILENAME
 *    INTEGER(KIND=4)                 :: SetOutputFileName
 *  END FUNCTION SetOutputFileName
 *  </PRE>
 *  </CODE>
 *  @endhtmlonly
 */
	IPQ_DLL_EXPORT IPQ_RESULT  SetOutputFileName(int id, const char* filename);

/**
 *  Sets the output file switch on or off.  This switch controls whether or not phreeqc
 *  writes to the <B><I>phreeqc.id.out</I></B> file.  This is the output normally generated
 *  when phreeqc is run.  The initial setting after calling @ref CreateIPhreeqc is off.
 *  @param id               The instance id returned from @ref CreateIPhreeqc.
 *  @param output_on        If non-zero, writes output to the output file; if zero, no output is written to the output file.
 *  @retval IPQ_OK          Success.
 *  @retval IPQ_BADINSTANCE The given id is invalid.
 *  @see                    GetOutputFileOn
 *  @par Fortran90 Interface:
 *  @htmlonly
 *  <CODE>
 *  <PRE>
 *  FUNCTION SetOutputFileOn(ID,OUT_ON)
 *    INTEGER(KIND=4),  INTENT(IN)  :: ID
 *    LOGICAL(KIND=4),  INTENT(IN)  :: OUT_ON
 *    INTEGER(KIND=4)               :: SetOutputFileOn
 *  END FUNCTION SetOutputFileOn
 *  </PRE>
 *  </CODE>
 *  @endhtmlonly
 */
	IPQ_DLL_EXPORT IPQ_RESULT  SetOutputFileOn(int id, int output_on);

/**
 *  Sets the output string switch on or off.  This switch controls whether or not the data normally sent
 *  to the output file are stored in a buffer for retrieval.  The initial setting after calling
 *  @ref CreateIPhreeqc is off.
 *  @param id                   The instance id returned from @ref CreateIPhreeqc.
 *  @param output_string_on     If non-zero, captures the phreeqc output into a string buffer;
 *                              if zero, phreeqc output is not captured to a string buffer.
 *  @retval IPQ_OK              Success.
 *  @retval IPQ_BADINSTANCE     The given id is invalid.
 *  @see                        GetOutputFileOn, GetOutputStringOn, GetOutputString, GetOutputStringLine, GetOutputStringLineCount, SetOutputFileOn
 *  @par Fortran90 Interface:
 *  @htmlonly
 *  <CODE>
 *  <PRE>
 *  FUNCTION SetOutputStringOn(ID,OUTPUT_STRING_ON)
 *    INTEGER(KIND=4),  INTENT(IN)  :: ID
 *    LOGICAL(KIND=4),  INTENT(IN)  :: OUTPUT_STRING_ON
 *    INTEGER(KIND=4)               :: SetOutputStringOn
 *  END FUNCTION SetOutputStringOn
 *  </PRE>
 *  </CODE>
 *  @endhtmlonly
 *
 *  @par C Example:
 *  see @ref GetOutputString_c "GetOutputString"
 *
 *  @par Fortran90 Example:
 *  see @ref GetOutputStringLine_f90 "GetOutputStringLine"
 */
	IPQ_DLL_EXPORT IPQ_RESULT  SetOutputStringOn(int id, int output_string_on);


/**
 *  Sets the name of the current selected output file (see @ref SetCurrentSelectedOutputUserNumber).  This file name is used if not specified within <B>SELECTED_OUTPUT</B> input.
 *  The default value is <B><I>selected_n.id.out</I></B>.
 *  @param id               The instance id returned from @ref CreateIPhreeqc.
 *  @param filename         The name of the file to write <B>SELECTED_OUTPUT</B> output to.
 *  @retval IPQ_OK          Success.
 *  @retval IPQ_BADINSTANCE The given id is invalid.
 *  @see                    GetSelectedOutputFileName, GetSelectedOutputFileOn, GetSelectedOutputString, GetSelectedOutputStringOn, GetSelectedOutputStringLine, GetSelectedOutputStringLineCount, SetCurrentSelectedOutputUserNumber, SetSelectedOutputFileOn, SetSelectedOutputStringOn
 *  @par Fortran90 Interface:
 *  @htmlonly
 *  <CODE>
 *  <PRE>
 *  FUNCTION SetSelectedOutputFileName(ID,FILENAME)
 *    INTEGER(KIND=4),   INTENT(IN)   :: ID
 *    CHARACTER(LEN=*),  INTENT(OUT)  :: FILENAME
 *    INTEGER(KIND=4)                 :: SetSelectedOutputFileName
 *  END FUNCTION SetSelectedOutputFileName
 *  </PRE>
 *  </CODE>
 *  @endhtmlonly
 */
	IPQ_DLL_EXPORT IPQ_RESULT  SetSelectedOutputFileName(int id, const char* filename);


/**
 *  Sets the selected-output file switch on or off.  This switch controls whether or not phreeqc writes output to
 *  the current <B>SELECTED_OUTPUT</B> file (see @ref SetCurrentSelectedOutputUserNumber). The initial
 *  setting after calling @ref CreateIPhreeqc is off.
 *  @param id               The instance id returned from @ref CreateIPhreeqc.
 *  @param sel_on           If non-zero, writes output to the selected-output file; if zero, no output is written to the selected-output file.
 *  @retval IPQ_OK          Success.
 *  @retval IPQ_BADINSTANCE The given id is invalid.
 *  @see                    GetSelectedOutputFileOn, GetSelectedOutputColumnCount, GetSelectedOutputRowCount, GetSelectedOutputValue, SetCurrentSelectedOutputUserNumber
 *  @par Fortran90 Interface:
 *  @htmlonly
 *  <CODE>
 *  <PRE>
 *  FUNCTION SetSelectedOutputFileOn(ID,SEL_ON)
 *    INTEGER(KIND=4),  INTENT(IN) :: ID
 *    LOGICAL(KIND=4),  INTENT(IN) :: SEL_ON
 *    INTEGER(KIND=4)              :: SetSelectedOutputFileOn
 *  END FUNCTION SetSelectedOutputFileOn
 *  </PRE>
 *  </CODE>
 *  @endhtmlonly
 */
	IPQ_DLL_EXPORT IPQ_RESULT  SetSelectedOutputFileOn(int id, int sel_on);

/**
 *  Sets the current selected output string switch on or off.  This switch controls whether or not the data normally sent
 *  to the current selected output file (see @ref SetCurrentSelectedOutputUserNumber) are stored in a buffer for retrieval.  The initial setting after calling
 *  @ref CreateIPhreeqc is off.
 *  @param id                   The instance id returned from @ref CreateIPhreeqc.
 *  @param sel_string_on        If non-zero, captures the output defined by the <B>SELECTED_OUTPUT</B> keyword into a string buffer;
 *                              if zero, output defined by the <B>SELECTED_OUTPUT</B> keyword is not captured to a string buffer.
 *  @retval IPQ_OK              Success.
 *  @retval IPQ_BADINSTANCE     The given id is invalid.
 *  @see                        GetSelectedOutputFileOn, GetSelectedOutputStringOn, GetSelectedOutputString, GetSelectedOutputStringLine, GetSelectedOutputStringLineCount, SetCurrentSelectedOutputUserNumber, SetSelectedOutputFileOn
 *  @par Fortran90 Interface:
 *  @htmlonly
 *  <CODE>
 *  <PRE>
 *  FUNCTION SetSelectedOutputStringOn(ID,SELECTED_OUTPUT_STRING_ON)
 *    INTEGER(KIND=4),  INTENT(IN)  :: ID
 *    LOGICAL(KIND=4),  INTENT(IN)  :: SELECTED_OUTPUT_STRING_ON
 *    INTEGER(KIND=4)               :: SetSelectedOutputStringOn
 *  END FUNCTION SetSelectedOutputStringOn
 *  </PRE>
 *  </CODE>
 *  @endhtmlonly
 *
 *  @par C Example:
 *  see @ref GetSelectedOutputString_c "GetSelectedOutputString"
 *
 *  @par Fortran90 Example:
 *  see @ref GetSelectedOutputStringLine_f90 "GetSelectedOutputStringLine"
 */
	IPQ_DLL_EXPORT IPQ_RESULT  SetSelectedOutputStringOn(int id, int sel_string_on);

// TODO int RunWithCallback(PFN_PRERUN_CALLBACK pfn_pre, PFN_POSTRUN_CALLBACK pfn_post, void *cookie, int output_on, int error_on, int log_on, int selected_output_on);


// TODO int CatchErrors(PFN_CATCH_CALLBACK pfn, void *cookie);


#if defined(__cplusplus)
}
#endif

#endif // INC_IPHREEQC_H
